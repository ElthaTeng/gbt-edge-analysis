import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from reproject import reproject_interp
from astropy.wcs import WCS
import warnings
import astropy.units as u
from astropy.coordinates import Angle
from scipy.optimize import curve_fit

table_gbt = np.genfromtxt('galaxy_parameters_erik_v2_fix62.csv', delimiter=',', skip_header=1, dtype='str')
table_ecaliga = np.genfromtxt('ecalifa_global.csv', delimiter=',', skip_header=118, dtype='str')

def power_law(X, a, b):
    return a * X + b 

def radius_arcsec(shape, w, ra, dec, pa, incl,
                  incl_correction=False, cosINCL_limit=0.5):
    # All inputs assumed as Angle
    if incl_correction and (np.isnan(pa.rad + incl.rad)):
        pa = Angle(0 * u.rad)
        incl = Angle(0 * u.rad)
        # Not written to the header
        msg = '\n::z0mgs:: PA or INCL is NaN in ' + \
            'radius calculation \n' + \
            '::z0mgs:: Setting both to zero.'
        # Warning message ends
        warnings.warn(msg, UserWarning)
        # Warning ends
    cosPA, sinPA = np.cos(pa.rad), np.sin(pa.rad)
    cosINCL = np.cos(incl.rad)
    if incl_correction and (cosINCL < cosINCL_limit):
        cosINCL = cosINCL_limit
        # Not written to the header
        msg = '\n::z0mgs:: Large inclination encountered in ' + \
            'radius calculation \n' + \
            '::z0mgs:: Input inclination: ' + str(incl.deg) + \
            ' degrees. \n' + \
            '::z0mgs:: cos(incl) is set to ' + str(cosINCL_limit)
        # Warning message ends
        warnings.warn(msg, UserWarning)
        # Warning ends
    xcm, ycm = ra.rad, dec.rad

    dp_coords = np.zeros(list(shape) + [2])
    # Original coordinate is (y, x)
    # :1 --> x, RA --> the one needed to be divided by cos(incl)
    # :0 --> y, Dec
    dp_coords[:, :, 0], dp_coords[:, :, 1] = \
        np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    # Now, value inside dp_coords is (x, y)
    # :0 --> x, RA --> the one needed to be divided by cos(incl)
    # :1 --> y, Dec
    for i in range(shape[0]):
        dp_coords[i] = Angle(w.wcs_pix2world(dp_coords[i], 1) * u.deg).rad
    dp_coords[:, :, 0] = 0.5 * (dp_coords[:, :, 0] - xcm) * \
        (np.cos(dp_coords[:, :, 1]) + np.cos(ycm))
    dp_coords[:, :, 1] -= ycm
    # Now, dp_coords is (dx, dy) in the original coordinate
    # cosPA*dy-sinPA*dx is new y
    # cosPA*dx+sinPA*dy is new x
    radius = np.sqrt((cosPA * dp_coords[:, :, 1] +
                      sinPA * dp_coords[:, :, 0])**2 +
                     ((cosPA * dp_coords[:, :, 0] -
                       sinPA * dp_coords[:, :, 1]) / cosINCL)**2)
    radius = Angle(radius * u.rad).arcsec
    return radius


galaxies = list(table_gbt[:,0])
bestfit_Z2r = np.full((len(galaxies), 4), np.nan) #slope, offset, e_slope, e_offset

for i, gal in enumerate(galaxies):

    print(gal)
    dist = table_gbt[i,24].astype('float') * 1e6

    # Get resolved Mst and Sigma_st maps
    pipe3d = fits.open(gal+'.Pipe3D.cube.fits')
    Mst_dustcorr = 10**pipe3d[1].data[19]
    Sigma_st = Mst_dustcorr / (dist/206265)**2

    # Grab Curti O3N2 O/H maps for all gals except N2596 where PG16S is used
    if gal=='NGC2596':  
        flux_elines = pipe3d[3] 
        I_N2_6548 = flux_elines.data[47]
        I_N2_6584 = flux_elines.data[46]
        I_S2_6717 = flux_elines.data[49]
        I_S2_6731 = flux_elines.data[50]
        I_O3_4959 = flux_elines.data[27]
        I_O3_5007 = flux_elines.data[26]
        I_Hb = flux_elines.data[28]

        N2 = (I_N2_6548 + I_N2_6584) / I_Hb
        S2 = (I_S2_6717 + I_S2_6731) / I_Hb
        R3 = (I_O3_4959 + I_O3_5007) / I_Hb

        # Make metallicity map based on PG15S
        OH_u = 8.424 + 0.03*np.log10(R3/S2) + 0.751*np.log10(N2) + (-0.349+0.182*np.log10(R3/S2) + 0.508*np.log10(N2))*np.log10(S2)
        OH_l = 8.072 + 0.789*np.log10(R3/S2) + 0.726*np.log10(N2) + (1.069-0.17*np.log10(R3/S2) + 0.022*np.log10(N2))*np.log10(S2)

        Zmap = np.full(N2.shape, np.nan)
        Zmap[np.log10(N2) >= -0.6] = OH_u[np.log10(N2) >= -0.6]
        Zmap[np.log10(N2) < -0.6] = OH_l[np.log10(N2) < -0.6]
        Zprime = 10**(Zmap - 8.5) 
    else:
        OH_map = fits.open('data/EDGE-CALIFA/eDR_OH_cubes/'+gal+'.OH.cube.fits')[0].data[34]
        Zprime = 10**(OH_map - 8.69)  

    # Create WCS for pipe3d's 2D map
    mask_pipe3d = fits.open('masks/from_matlab/mask_'+gal+'_block_v3.fits')[0]
    pipe3d_header_2d = mask_pipe3d.header
    pipe3d_header_2d['NAXIS'] = 2
    keywords_to_remove = ['NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CDELT3', 'CUNIT3']
    for key in keywords_to_remove:
        del pipe3d_header_2d[key]
    #print(WCS(pipe3d_header_2d))

    # Get deprojected radii map (in arcsec)
    ra = table_ecaliga[table_ecaliga[:,1]==gal][0,2] + 'd'
    dec = table_ecaliga[table_ecaliga[:,1]==gal][0,3] + 'd'
    pa = table_gbt[0,2] + 'd'
    incl = table_gbt[0,3] + 'd'
    angles = Angle([ra, dec, pa, incl])

    wcs = WCS(pipe3d_header_2d)
    dshape = Zprime.shape
    radii_map = radius_arcsec(dshape, wcs, angles[0], angles[1], angles[2], angles[3])
    #np.save('maps/radii_maps/'+gal+'_radii_map.npy', radii_map)

    # fit Zprime-radius relation
    valid = ~np.isnan(Zprime)
    data_x = radii_map[valid]
    data_y = np.log10(Zprime[valid])

    ic = [-0.1, 1.]
    popt, pcov = curve_fit(power_law, data_x, data_y, p0 = ic)
    errors = np.sqrt(np.diagonal(pcov))
    bestfit_Z2r[i, :2] = popt
    bestfit_Z2r[i, 2:] = errors

    plt.scatter(data_x, data_y, facecolor='grey', s=25, marker='.')
    xrange = np.arange(np.nanmin(data_x), np.nanmax(data_x), 0.2)
    plt.title(gal+': y = '+str(round(popt[0],5))+'x + '+str(round(popt[1], 5)))
    plt.plot(xrange, xrange*popt[0] + popt[1], 'k--', lw=2)
    plt.ylim(-1, 1)

    plt.tick_params(axis="y", labelsize=14, labelleft=True)
    plt.tick_params(axis="x", labelsize=14, labelbottom=True)
    plt.xlabel('Radius (arcsec)', fontsize=16) 
    plt.ylabel(r'log (Z / Z$_\odot$)', fontsize=16)
    #plt.savefig('plots/fit_Zprime_gradient/'+gal+'_Zvsr_fitting_Curti.png', bbox_inches='tight', pad_inches=0.02)
    plt.clf()

np.save('bestfit_Zprime_gradient_Curti_fix62.npy', bestfit_Z2r)
np.savetxt('bestfit_Zprime_gradient_Curti_fix62.csv', bestfit_Z2r, delimiter=',')