import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MaxNLocator
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle
import re
 
table_gbt = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')
sampletable = np.genfromtxt('GBTEDGE.cat', skip_header=2, dtype='str')
galaxies = list(table_gbt[:,0])
sessions = table_gbt[:,5]
masks = table_gbt[:,7]

datadir = './maps/'  ### directory of the data
plotdir = './plots/moments/' ### folder in which the plot will be saved

fig = plt.figure(figsize=(15, 4))

# gal = 'NGC0001'
# i = 0

for i, gal in enumerate(galaxies):
    
    print(gal)

    ## Load moment maps
    if gal == 'Mrk1418':
        mom0 = fits.open(datadir + gal + '_12CO_mom0_' + masks[i] + '_se' + sessions[i] + '.fits')[0].data
        emom0 = fits.open(datadir + gal + '_12CO_emom0_' + masks[i] + '_se' + sessions[i] + '.fits')[0].data
    else:
        mom0 = fits.open(datadir + gal + '_12CO_mom0_' + masks[i] + '_se' + sessions[i] + '_noR25.fits')[0].data
        emom0 = fits.open(datadir + gal + '_12CO_emom0_' + masks[i] + '_se' + sessions[i] + '_noR25.fits')[0].data
    
    wcs = WCS(fits.open(datadir + gal + '_12CO_mom0_' + masks[i] + '_se' + sessions[i] + '.fits')[0].header)
    mom1 = fits.open(datadir + gal + '_12CO_mom1_' + masks[i] + '_se' + sessions[i] + '.fits')[0].data
    mom2 = fits.open(datadir + gal + '_12CO_vdisp_' + masks[i] + '_se' + sessions[i] + '.fits')[0].data 

    ra = sampletable[sampletable[:,0]==gal][0,2]
    dec = sampletable[sampletable[:,0]==gal][0,3]
    ra_num = np.array(re.findall('\d+\.\d+|\d+', ra), dtype='float64')
    dec_num = np.array(re.findall('\d+\.\d+|\d+', dec), dtype='float64')
    ra_deg = Angle(ra, unit=u.hour).degree
    dec_deg = Angle(dec, unit=u.degree).degree

    if gal == 'NGC0169':  
        vmax_mom0 = 100
    elif gal == 'CGCG536-030':
        vmax_mom0 = 200
    else:
        vmax_mom0 = 50  

    vmin_mom1 = np.nanpercentile(mom1, 5)  
    vmax_mom1 = np.nanpercentile(mom1, 95)  
    vmin_mom2 = 10 
    vmax_mom2 = 150
    
    if gal=='UGC04136': # ensure a minimum thickness for this galaxy due to its high inclination
        mask_r25 = np.sum(fits.open('masks/'+gal+'_mask_thickened_block.fits')[0].data, axis=0) > 0
    else:
        mask_r25 = fits.open('maps/'+gal+'_12CO_emom0_block_se'+sessions[i]+'.fits')[0].data > 0

    plt.suptitle(gal, x=0.05, y=0.95, weight='bold', fontsize=16)

    ## Moment 0
    ax = fig.add_subplot(1, 3, 1, projection=wcs)
    ra = ax.coords[0]
    ra.set_major_formatter('hh:mm:ss.s')

    #mom0[mom0==0] = np.nan
    cmapi = plt.cm.get_cmap('hot').copy()  
    cmapi.set_bad(color='dimgray')

    plt.imshow(mom0, cmap=cmapi, origin='lower', vmin=0, vmax=vmax_mom0)  
    plt.tick_params(axis="y", labelsize=12, labelleft=True)
    plt.tick_params(axis="x", labelsize=12, labelbottom=True)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    cb.ax.set_title(r'K km s$^{-1}$', fontsize=12, loc='right')
    
    sigma = np.mean(emom0[mask_r25])
    snr = mom0 / emom0
    levelK = [5] # S/N ratio
    ax.contour(snr, levels = levelK, colors = ['white'], alpha = 0.5, linewidths = 1.5)
    ax.contour(mask_r25, levels = [0.5], colors = ['cyan'], alpha = 0.5, linewidths = 1.5)
    ax.contourf(mask_r25, levels=[-0.5, 0.5], hatches=['xx', ''], alpha=0.15)

    # Set image FOV coverage
    fits_map = fits.open(datadir + gal + '_12CO_mom0_' + masks[i] + '_se' + sessions[i] + '.fits')[0]
    cornerleftx = 80*u.arcsec # including all GBT GOV
    cornerlefty = 80*u.arcsec
    offsetx = cornerleftx.to(u.deg)
    offsety = cornerlefty.to(u.deg)
    posx1, posy1 = wcs.all_world2pix((ra_deg*u.deg+offsetx).value, (dec_deg*u.deg + offsety).value, 1)
    posx2, posy2 = wcs.all_world2pix((ra_deg*u.deg-offsetx).value, (dec_deg*u.deg - offsety).value, 0)
    racheck, deccheck = wcs.all_world2pix(ra_deg, dec_deg, 0)
    cdelt1 = abs(fits_map.header['cdelt1']) * u.deg
    cdelt2 = abs(fits_map.header['cdelt2']) * u.deg
    crpix1 = abs(fits_map.header['crpix1'])
    crpix2 = abs(fits_map.header['crpix2'])
    deltapix_x = offsetx/cdelt1
    deltapix_y = offsety/cdelt2
    posx1_check = racheck - deltapix_x.value
    posx2_check = racheck + deltapix_x.value
    posy1_check = deccheck - deltapix_y.value
    posy2_check = deccheck + deltapix_y.value
    ax.set_xlim(posx1_check, posx2_check)
    ax.set_ylim(posy1_check, posy2_check)

    # Add beam size at the lower left corner
    b1 = abs(fits_map.header['bmaj']) * u.deg
    b2 = abs(fits_map.header['bmin']) * u.deg
    bpa = abs(fits_map.header['BPA']) * u.deg
    panel_width = posx2_check-posx1_check
    panel_heigth = posy2_check-posy1_check
    center_beamx  = posx1_check + 0.1*panel_width
    center_beamy  = posy1_check + 0.1*panel_heigth
    beam_sky = wcs.pixel_to_world(center_beamx, center_beamy)
    rect_width = 0.1 *panel_width
    rect_height = 0.1 * panel_heigth
    pos_rect_x = center_beamx - rect_width/2
    pos_rect_y = center_beamy - rect_width/2
    rect = Rectangle((pos_rect_x, pos_rect_y), rect_width, rect_height,
                edgecolor='black', facecolor='white')
    ax.add_patch(rect)
    beam = Ellipse((center_beamx, center_beamy), width = b1.value/cdelt1.value, height = b2.value/cdelt2.value, angle = bpa.value,
                        edgecolor = 'black', facecolor = 'black', alpha = 0.9)    
    ax.add_patch(beam)

    plt.title(r'Moment 0', fontsize=16) 
    plt.xlabel('R.A. (J2000) ', fontsize=14) # 
    plt.ylabel('Decl. (J2000)', fontsize=14) #


    ## Moment 1
    ax = fig.add_subplot(1, 3, 2, projection=wcs)
    ra = ax.coords[0]
    ra.set_major_formatter('hh:mm:ss.s')

    cmapi = plt.cm.get_cmap('coolwarm').copy()  
    cmapi.set_bad(color='dimgray')

    plt.imshow(mom1, cmap=cmapi, origin='lower', vmin=vmin_mom1, vmax=vmax_mom1)  
    plt.tick_params(axis="y", labelsize=12, labelleft=False)
    plt.tick_params(axis="x", labelsize=12, labelbottom=True)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    cb.ax.set_title(r'km s$^{-1}$', fontsize=12)

    ax.contour(mask_r25, levels = [0.5], colors = ['cyan'], alpha = 0.5, linewidths = 1.5)

    plt.title(r'Moment 1', fontsize=16)   
    plt.xlabel('R.A. (J2000) ', fontsize=14) # 

    # Inset plot for integrated spectrum within the mom1 panel 
    #ax_inset = inset_axes(ax, width="98%", height="24%", loc="lower center")
    ax_inset = ax.inset_axes([0.01, 0.01, 0.98, 0.24])
    cube = fits.open('data/' + gal + '_12CO_rebase7_smooth1.3_hanning2_se' + sessions[i] + '.fits')[0].data
    header = fits.open('data/' + gal + '_12CO_rebase7_smooth1.3_hanning2_se' + sessions[i] + '.fits')[0].header

    v_ref = header['CRVAL3']
    v_del = header['CDELT3']
    RefPix = header['CRPIX3']
    N_vel = header['NAXIS3']
    v_list = np.arange(v_ref-v_del*(RefPix-1), v_ref+v_del*(-RefPix+N_vel+0.5), v_del)

    mom0 = fits.open('maps/' + gal + '_12CO_mom0_' + masks[i] + '_se' + sessions[i] + '.fits')[0].data
    emom0 = fits.open('maps/' + gal + '_12CO_emom0_' + masks[i] + '_se' + sessions[i] + '.fits')[0].data
    snr = mom0 / emom0
    mask = np.copy(mom0)
    mask[np.isnan(mom0)] = 0
    if gal != 'Mrk1418':   
        mask[snr < 5] = 0  # Comment out this to include all pixels within R25 contour

    indices = np.stack((np.nonzero(mask)[0],np.nonzero(mask)[1]),axis=1)
    N_pix = indices.shape[0]

    spec_stacked_direct = np.zeros((cube.shape[0],))

    for i in range(N_pix):
        spec = cube[:, indices[i,0], indices[i,1]]
        spec_stacked_direct += spec

    spec_avg_direct = spec_stacked_direct/N_pix
    
    ax_inset.step(v_list,spec_avg_direct,c='#8E8E8E',lw=2,where='mid')
    ax_inset.tick_params(top=True, labeltop=True, labelleft=False, bottom=False, labelbottom=False, direction="in", pad=-15, labelsize=8)
    ax_inset.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.set_xticks(ax.get_xticks()[::5])  
    ax_inset.set_yticks([])
    ax_inset.xaxis.set_label_position('top')
    ax_inset.set_xlabel(r'velocity (km s$^{-1}$)', fontsize=8, loc='right')

    ax.contour(snr, levels = levelK, colors = ['white'], alpha = 0.5, linewidths = 1.5)


    ## Moment 2
    ax = fig.add_subplot(1, 3, 3, projection=wcs)
    ra = ax.coords[0]
    ra.set_major_formatter('hh:mm:ss.s')

    cmapi = plt.cm.get_cmap('coolwarm').copy()  
    cmapi.set_bad(color='dimgray')

    plt.imshow(mom2, cmap=cmapi, origin='lower', vmin=vmin_mom2, vmax=vmax_mom2)  
    plt.tick_params(axis="y", labelsize=12, labelleft=False)
    plt.tick_params(axis="x", labelsize=12, labelbottom=True)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    cb.ax.set_title(r'km s$^{-1}$', fontsize=12)

    ax.contour(mask_r25, levels = [0.5], colors = ['cyan'], alpha = 0.5, linewidths = 1.5)
    ax.contour(snr, levels = levelK, colors = ['white'], alpha = 0.5, linewidths = 1.5)

    plt.title(r'Moment 2', fontsize=16)   
    plt.xlabel('R.A. (J2000) ', fontsize=14) # 

    plt.savefig(plotdir + gal + '_moments.png', bbox_inches='tight', pad_inches=0.1)  
    plt.clf()  # Clear the figure for the next galaxy

    #plt.show()

