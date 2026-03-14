import numpy as np
import matplotlib.pyplot as plt
import os
### plot all moment maps
from astropy.table import Table
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
### previously run the pvcreation script in CASA
from astropy.coordinates import Angle, SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as colors
import re
import matplotlib as mpl
#from regions import PixCoord, RectanglePixelRegion


### environment M51en

#########################################################
#### structure
#########################################################

############
### first load all info we need:

datadir = './maps/'  ### directory of the data
plotdir = './plots/' ### folder in which the plot will be saved
#os.system('mkdir ' + plotdir)

table = np.genfromtxt('galaxy_parameters_erik_v2_fix62.csv', delimiter=',', skip_header=1, dtype='str')

### all galaxies we want to plot.
galaxies = list(table[:,0])
posang = list(table[:,2])
sessions = list(table[:,5])
### SET THIS !!!
masks = list(table[:,7])  # [:,6] if blockpref, [:,7] if datapref

### units: arcsec
# FOV1     = [100,40,40,40,30,40,30,30,30,40,
#             40,50,30,50,40,40,40,40,25,40,
#             25,40,40,40,20,40,50,25,25,20]

# FOV2     = [100,40,40,40,30,40,30,30,30,40,
#             40,30,50,30,40,40,40,40,60,40,
#             60,40,40,40,40,40,50,60,50,40]

### for these galaxies, use circle with 25 arcsec diameter instead of rectangle:
circlegals = []

###########################
### additional info e.g. distance from sampletable
sampletable = np.genfromtxt('GBTEDGE.cat', skip_header=2, dtype='str')

distancelist = []
ras          = []
decs         = []
pa_gallist   = []
incllist     = []

for nameG in galaxies:
    ra             = sampletable[sampletable[:,0]==nameG][0,2]
    dec            = sampletable[sampletable[:,0]==nameG][0,3]
    ra_num = np.array(re.findall('\d+\.\d+|\d+', ra), dtype='float64')
    dec_num = np.array(re.findall('\d+\.\d+|\d+', dec), dtype='float64')
    # print(ra, dec, ra_num, dec_num)
    # ra_deg = Angle(tuple(ra_num), unit=u.hour).degree
    # dec_deg = Angle(tuple(dec_num), unit=u.degree).degree
    ra_deg = Angle(ra, unit=u.hour).degree
    dec_deg = Angle(dec, unit=u.degree).degree
    #print(ra_deg, dec_deg)

    # dist           = sampletable[sampletable[:,0]==nameG]['dist']
    # pa_disk        = sampletable[sampletable[:,0]==nameG]['orient_posang']
    # incl           = sampletable[sampletable[:,0]==nameG]['orient_incl']


    ras.append(ra_deg)
    decs.append(dec_deg)
    #distancelist.append(dist[0])
    #pa_gallist.append(pa_disk[0])
    #incllist.append(incl[0])


T = Table([galaxies, posang, ras, decs, masks, sessions], names = ('Galaxy', 'pa', 'ra', 'dec', 'Mask', 'Session'))
#T.write('galaxy_basics_all71.csv', format='csv')

##########################################################
#### a) Create gallerie of all tracers for mom0 and tpeak
##########################################################

### plot a sqrt colorscale:
do_norm = False

### to plot rectangle use matplotlib.patches vs regions
do_other_geometry = False

colors1 = plt.cm.binary(np.linspace(0., 1, 128))
colors2 = plt.cm.gist_heat(np.linspace(0, 1, 128))
combined_colors = np.vstack((colors1, colors2))
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', combined_colors)

cmapi = plt.cm.get_cmap('hot').copy()  # mymap  #
cmapi.set_bad(color='dimgray')

fig = plt.figure(figsize=(12, 14.5))
gs = fig.add_gridspec(nrows = 9, ncols = 7, height_ratios=(1,1,1,1,1,1,1,1,1), width_ratios = (1,1,1,1,1,1,1),
                      left=0.01, right=0.99, bottom=0.01, top=0.99,
                      wspace=0.03, hspace=0.03)

for i, row in enumerate(T):
    #filepath    = os.path.join(datadir, row['Galaxy'].lower() + '_12m+7m+tp_co21_strict_mom0.fits') if not 'ngc1808' in row['Galaxy'] else os.path.join(datadir, 'ngc1808/ngc1808_12m+7m+tp_co21lores_strict_mom0.fits')
    # if row['Galaxy'] in ['NGC2596', 'NGC5929', 'UGC08322', 'IC3598']:
    #     continue
    # else:
    filepath    = os.path.join(datadir, row['Galaxy'] + '_12CO_mom0_' + row['Mask'] + '_se' + row['Session'] + '_noR25.fits')  #_noR25
    emom0 = fits.open('maps/'+row['Galaxy']+'_12CO_emom0_'+row['Mask']+'_se'+row['Session']+'_noR25.fits')[0].data  #_noR25

    ## Uncomment this if running "noR25" version of mom0 map gallery
    if 'Mrk1418' in row['Galaxy']:
        filepath    = os.path.join(datadir, row['Galaxy'] + '_12CO_mom0_' + row['Mask'] + '_se' + row['Session'] + '.fits')
        emom0 = fits.open('maps/'+row['Galaxy']+'_12CO_emom0_'+row['Mask']+'_se'+row['Session']+'.fits')[0].data

    mapi        = fits.open(filepath)[0]
    wcs_map     = WCS(mapi.header)
    mapdata     = mapi.data        

    vmini = 0  #np.nanpercentile(mapdata, 1)
    mapdata[mapdata==0] = np.nan

    if do_norm:
        vmini = 0.01 if vmini <=0 else vmini
        vmininorm = np.nanpercentile(mapdata, 10)
        vmininorm = 0.01 if vmininorm <=0 else vmininorm

    if row['Galaxy'] in ['NGC0169']:  
        vmaxi = 100
    elif row['Galaxy'] in ['CGCG536-030']:
        vmaxi = 200
    else:
        vmaxi = 50  #40  #np.nanpercentile(mapdata, 99.9) 

    ax = fig.add_subplot(gs[i//7, i%7], projection = wcs_map)

    if do_norm: ### use a sqrt colorscale intsead of linear.
        def _forward(x):
            return np.sqrt(x)
        def _inverse(x):
            return x**2
        normsqrt = colors.FuncNorm((_forward,_inverse), vmin=vmininorm, vmax = vmaxi)
        im    = ax.imshow(mapdata, cmap = cmapi, norm=normsqrt)
    else: ### linear colorscale
        im    = ax.imshow(mapdata, cmap = cmapi, vmin = vmini, vmax = vmaxi)


    ### add contours if required
    mask_r25 = fits.open('maps/'+row['Galaxy']+'_12CO_emom0_block_se'+row['Session']+'.fits')[0].data > 0
    sigma = np.mean(emom0[mask_r25])
    snr = mapdata / emom0
    levelK = [5] # S/N ratio
    ax.contour(snr, levels = levelK, colors = ['white'], alpha = 0.5, linewidths = 1.5)
    ax.contour(mask_r25, levels = [0.5], colors = ['cyan'], alpha = 0.5, linewidths = 1.5)

    ax.contourf(mask_r25, levels=[-0.5, 0.5], hatches=['xxx', ''], alpha=0) 

    ax.tick_params(axis = 'y', labelleft = False, labelright = False, left = False, right =False, direction = 'in', color = 'white')
    ax.tick_params(axis = 'x',  labeltop = False, labelbottom = False, top = False, bottom = False, direction = 'in', color = 'white')# which = 'both', direction = 'in')
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')

    ##################################
    #### add center position as plus:
    #ax.scatter(row['ra'], row['dec'], transform=ax.get_transform('fk5'), s=30, marker = 'x', alpha = 0.8, facecolor='red', zorder = 100)

    ##################################
    ########## CUT THE IMAGE #########

    cornerleftx = 50*u.arcsec
    cornerlefty = 50*u.arcsec

    offsetx = cornerleftx.to(u.deg)
    offsety = cornerlefty.to(u.deg)
    posx1, posy1 = wcs_map.all_world2pix((row['ra']*u.deg+offsetx).value, (row['dec']*u.deg + offsety).value, 1)
    posx2, posy2 = wcs_map.all_world2pix((row['ra']*u.deg-offsetx).value, (row['dec']*u.deg - offsety).value, 0)


    ### try different set of FOV measurements:
    racheck, deccheck = wcs_map.all_world2pix((row['ra']), (row['dec']), 0)

    cdelt1 = abs(mapi.header['cdelt1']) * u.deg
    cdelt2 = abs(mapi.header['cdelt2']) * u.deg
    crpix1 = abs(mapi.header['crpix1'])
    crpix2 = abs(mapi.header['crpix2'])
    deltapix_x = offsetx/cdelt1
    deltapix_y = offsety/cdelt2

    posx1_check = racheck - deltapix_x.value
    posx2_check = racheck + deltapix_x.value
    posy1_check = deccheck - deltapix_y.value
    posy2_check = deccheck + deltapix_y.value

    ax.set_xlim(posx1_check, posx2_check)
    ax.set_ylim(posy1_check, posy2_check)


    ####################################
    ### add beam size!!!################
    b1 = abs(mapi.header['bmaj']) * u.deg
    b2 = abs(mapi.header['bmin']) * u.deg
    bpa = abs(mapi.header['BPA']) * u.deg

    ### define the position of the ellipse:
    panel_width = posx2_check-posx1_check
    panel_heigth = posy2_check-posy1_check

    center_beamx  = posx1_check + 0.1*panel_width
    center_beamy  = posy1_check + 0.1*panel_heigth

    beam_sky = wcs_map.pixel_to_world(center_beamx, center_beamy)

    ### width of the rectangle = 0.1 image width??
    rect_width = 0.15 *panel_width
    rect_height = 0.15 * panel_heigth
    pos_rect_x = center_beamx - rect_width/2
    pos_rect_y = center_beamy - rect_width/2


    rect = Rectangle((pos_rect_x, pos_rect_y), rect_width, rect_height,
                edgecolor='black', facecolor='white',
                )
    
    if i==0:  # Add beam size only on the first panel
        ax.add_patch(rect)

        # beam = Ellipse((beam_sky.ra.value, beam_sky.dec.value), width = b1.value, height = b2.value, angle = bpa.value,
        #                 edgecolor = 'black', facecolor = 'black', alpha = 0.9,
        #                 transform = ax.get_transform('fk5')) 
        
        beam = Ellipse((center_beamx, center_beamy), width = b1.value/cdelt1.value, height = b2.value/cdelt2.value, angle = bpa.value,
                        edgecolor = 'black', facecolor = 'black', alpha = 0.9
                        )
        
        ax.add_patch(beam)

    #print(beam_sky)

    labely = row['Galaxy'].upper()
    # starneed = ['NGC0169', 'NGC0495', 'UGC8909', 'CGCG536-030', 'NGC5954', 'UGC09777']
    # labely += r'$\star$' if row['Galaxy'] in starneed else ''

    ### text position: do it in percent!
    if row['Galaxy'] in ['NGC0932', 'NGC3406NED01', 'NGC5216', 'NGC5631', 'UGC02222', 'UGC03960', 'UGC08234', 'UGC09629']:
        ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='r', fontweight = 'bold')
    elif row['Galaxy'] in ['NGC3106', 'NGC3619', 'NGC5157', 'NGC6154', 'NGC6338', 'UGC04136', 'UGC08322', 'UGC10097', 'UGC10905', 'IC0674', 'IC3598']:
        ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='lime', fontweight = 'bold')
    else:
        ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='w', fontweight = 'bold') #, bbox = dict(boxstyle='square', facecolor='white', alpha=0.7)


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
axleg = fig.add_subplot(gs[8, 6])
axleg.spines['right'].set_visible(False)
axleg.spines['left'].set_visible(False)
axleg.spines['top'].set_visible(False)
axleg.spines['bottom'].set_visible(False)
axleg.tick_params(bottom = False, left = False, labelbottom = False, labelleft = False)
axleg.patch.set_alpha(0.0)
axleg.set_xlim(0,100)
axleg.set_ylim(0,100)

cmap = mpl.colormaps['hot']  # mymap  #
norm = mpl.colors.Normalize(vmin=0, vmax=50)
#norm = colors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=70)
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Dummy array is needed for proper initialization
cbar = fig.colorbar(sm, ax=axleg, orientation='horizontal')
cbar.ax.set_xscale('linear')

# reg = RectanglePixelRegion(PixCoord(x=10, y=90), width=10,
#                    height=10, angle=0*u.deg)
# patch = reg.plot(ax=axleg, facecolor='none', edgecolor='lime', lw=2.5,
#                  label='Rectangle')
ell = Ellipse((10,60), width = 10, height = 10, angle = (0),
                edgecolor = 'c', facecolor = 'none', alpha = 1, lw=2.5)
axleg.add_patch(ell)

axleg.text(25,57, s= r'R$_{25} \times$FoV', fontweight='bold', fontsize=12)

ell2 = Ellipse((10,30), width = 10, height = 10, angle = (0),
                edgecolor = 'silver', facecolor = 'none', alpha = 1, lw=2.5)
axleg.add_patch(ell2)

axleg.text(25,27, s= 'S/N > 5', fontweight='bold', fontsize=12)


#####################
### save the plot
savename = os.path.join(plotdir, 'masked_mom0_gallery')

if do_other_geometry:
    savename += '_geom'
if do_norm:
    savename += '_sqrt'

fig.savefig(savename + '_datapref_erik_v2_fix62_noR25_max50_color_hatch_new.pdf', dpi = 300)
