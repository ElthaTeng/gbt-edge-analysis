import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
### plot all moment maps
from astropy.table import Table
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from astropy.coordinates import Angle, SkyCoord
import matplotlib.colors as colors


plotdir = './plots/'

table = np.genfromtxt('galaxy_parameters_erik_v2_fix62.csv', delimiter=',', skip_header=1, dtype='str')
galaxies = list(table[:,0])

fig = plt.figure(figsize=(12, 14.5))

gs = fig.add_gridspec(nrows = 9, ncols = 7, height_ratios=(1,1,1,1,1,1,1,1,1), width_ratios = (1,1,1,1,1,1,1),
                      left=0.01, right=0.99, bottom=0.01, top=0.99,
                      wspace=0.03, hspace=0.03)

for i, gal in enumerate(galaxies):      

    img = mpimg.imread('cutouts/'+gal+'.jpg')

    ax = fig.add_subplot(gs[i//7, i%7])

    im = ax.imshow(img)

    ax.tick_params(axis="x", labelsize=14, bottom=False, labelbottom=False)
    ax.tick_params(axis="y", labelsize=14, left=False, labelleft=False)


    ### add contours if required
    # emom0 = fits.open('maps/'+row['Galaxy']+'_12CO_emom0_'+row['Mask']+'_se'+row['Session']+'.fits')[0].data
    # mask_r25 = fits.open('maps/'+row['Galaxy']+'_12CO_emom0_block_se'+row['Session']+'.fits')[0].data > 0
    # sigma = np.mean(emom0[mask_r25])
    # snr = mapdata / emom0
    # levelK = [2] # 2*sigma
    # ax.contour(snr, levels = levelK, colors = ['white'], alpha = 0.5, linewidths = 1.5)
    # ax.contour(mask_r25, levels = [1], colors = ['cyan'], alpha = 0.5, linewidths = 1.5)

    # ax.tick_params(axis = 'y', labelleft = False, labelright = False, left = False, right =False, direction = 'in', color = 'white')
    # ax.tick_params(axis = 'x',  labeltop = False, labelbottom = False, top = False, bottom = False, direction = 'in', color = 'white')# which = 'both', direction = 'in')
    # ax.set_xlabel(' ')
    # ax.set_ylabel(' ')

    ##################################
    #### add center position as plus:
    #ax.scatter(row['ra'], row['dec'], transform=ax.get_transform('fk5'), s=30, marker = 'x', alpha = 0.8, facecolor='red', zorder = 100)

    ##################################
    ########## CUT THE IMAGE #########

    # cornerleftx = 50*u.arcsec
    # cornerlefty = 50*u.arcsec

    # offsetx = cornerleftx.to(u.deg)
    # offsety = cornerlefty.to(u.deg)
    # posx1, posy1 = wcs_map.all_world2pix((row['ra']*u.deg+offsetx).value, (row['dec']*u.deg + offsety).value, 1)
    # posx2, posy2 = wcs_map.all_world2pix((row['ra']*u.deg-offsetx).value, (row['dec']*u.deg - offsety).value, 0)


    ### try different set of FOV measurements:
    # racheck, deccheck = wcs_map.all_world2pix((row['ra']), (row['dec']), 0)

    # cdelt1 = abs(mapi.header['cdelt1']) * u.deg
    # cdelt2 = abs(mapi.header['cdelt2']) * u.deg
    # crpix1 = abs(mapi.header['crpix1'])
    # crpix2 = abs(mapi.header['crpix2'])
    # deltapix_x = offsetx/cdelt1
    # deltapix_y = offsety/cdelt2

    # posx1_check = racheck - deltapix_x.value
    # posx2_check = racheck + deltapix_x.value
    # posy1_check = deccheck - deltapix_y.value
    # posy2_check = deccheck + deltapix_y.value

    # ax.set_xlim(posx1_check, posx2_check)
    # ax.set_ylim(posy1_check, posy2_check)


    labely = gal
    starneed = ['ngc0628', 'ngc4689', 'ngc3521', 'ngc4254', 'ngc5248', 'ngc4941', 'ngc4536', 'ngc4569']
    labely += r'$\star$' if gal in starneed else ''

    ### text position: do it in percent!
    #ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='w', fontweight = 'bold') #, bbox = dict(boxstyle='square', facecolor='white', alpha=0.7)

    if gal in ['NGC0932', 'NGC3406NED01', 'NGC5216', 'NGC5631', 'UGC02222', 'UGC03960', 'UGC08234', 'UGC09629']:
        ax.text(0.05*img.shape[1], 0.15*img.shape[0], s= labely, color='r', fontweight = 'bold')
    elif gal in ['NGC3106', 'NGC3619', 'NGC5157', 'NGC6154', 'NGC6338', 'UGC04136', 'UGC08322', 'UGC10097', 'UGC10905', 'IC0674', 'IC3598']:
        ax.text(0.05*img.shape[1], 0.15*img.shape[0], s= labely, color='lime', fontweight = 'bold')
    else:
        ax.text(0.05*img.shape[1], 0.15*img.shape[0], s= labely, color='w', fontweight = 'bold')

    ##################################
    #### add rectangle!
    ### FOV should be a green rectangle! PA is given and lengths in arcsec! OR a circle for some galaxies:

    # if not row['Galaxy'] in circlegals:
    #     width_rect = (row['fov1']*u.arcsec).to(u.deg)
    #     height_rect = (row['fov2']*u.arcsec).to(u.deg)
    #     posang_rect = (row['pa'])*u.deg

    #     #### do in pixel coord:
    #     width_rect_pix = (row['fov1']*u.arcsec).to(u.deg)/cdelt1
    #     height_rect_pix = (row['fov2']*u.arcsec).to(u.deg)/cdelt2
    #     reg = RectanglePixelRegion(PixCoord(x=racheck, y=deccheck), width=width_rect_pix.value,
    #                        height=height_rect_pix.value, angle=posang_rect)

    #     patch = reg.plot(ax=ax, facecolor='none', edgecolor='lime', lw=1,
    #                      label='Rectangle', zorder = 200)
    # else:
    #     r = SphericalCircle((row['ra'] * u.deg, row['dec']* u.deg), (25/2*u.arcsec).to(u.deg),
    #                          edgecolor='lime', facecolor='none',
    #                          transform=ax.get_transform('fk5'))

    #     ax.add_patch(r)


    ##################################################################
    ### add ellipse that indicates central 3kpc!! need distance here!

    # rad_3kpc  = 3 * u.kpc
    # rad_major = np.arctan(rad_3kpc/ (row['dist']*u.Mpc)).to(u.deg) ### should be the larger one!
    # rad_minor = abs(rad_major * np.cos(row['incl']*u.deg)) ### should be the smaller one!

    # if not do_other_geometry:

    #     ell = Ellipse((row['ra'], row['dec']), width = rad_minor.value, height = rad_major.value, angle = (360-row['pagal']),
    #                     edgecolor = 'red', facecolor = 'none', alpha = 0.9,
    #                     transform = ax.get_transform('fk5'))
    #     ax.add_patch(ell)

    # if do_other_geometry:
    #     #### better to use: EllipseSkyRegion
    #     center_sky = SkyCoord(row['ra'], row['dec'], unit='deg', frame='fk5')
    #     region_Ell = EllipseSkyRegion(center=center_sky,
    #                                   height=rad_major, width=rad_minor,
    #                                   angle=row['pagal']*u.deg)
    #     pixel_Ell = region_Ell.to_pixel(wcs_map)
    #     pixel_Ell.plot(edgecolor = 'red')


###############################
#### add legend
# axleg = fig.add_subplot(gs[1,8])
# axleg.spines['right'].set_visible(False)
# axleg.spines['left'].set_visible(False)
# axleg.spines['top'].set_visible(False)
# axleg.spines['bottom'].set_visible(False)
# axleg.tick_params(bottom = False, left = False, labelbottom = False, labelleft = False)
# axleg.patch.set_alpha(0.0)
# axleg.set_xlim(0,100)
# axleg.set_ylim(0,100)

# reg = RectanglePixelRegion(PixCoord(x=10, y=90), width=10,
#                    height=10, angle=0*u.deg)
# patch = reg.plot(ax=axleg, facecolor='none', edgecolor='lime', lw=1,
#                  label='Rectangle')
# axleg.text(25,87, s= 'FoV')

# ell = Ellipse((10,70), width = 10, height = 10, angle = (0),
#                 edgecolor = 'red', facecolor = 'none', alpha = 0.9)
# axleg.add_patch(ell)

# axleg.text(25,67, s= r'$3\,$kpc diameter')



#####################
### save the plot
savename = os.path.join(plotdir, 'sdss_gallery')

fig.savefig(savename + '_fix62_color.pdf', dpi = 300)
