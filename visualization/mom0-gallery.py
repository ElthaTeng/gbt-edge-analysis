import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from astropy.coordinates import Angle
import matplotlib.colors as colors
import re
import matplotlib as mpl


#########################################################
#### structure
#########################################################

############
### first load all info we need:

datadir = './maps/'  ### directory of the data
plotdir = './plots/' ### folder in which the plot will be saved

table = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')

### all galaxies we want to plot.
galaxies = list(table[:,0])
posang = list(table[:,2])
sessions = list(table[:,5])
masks = list(table[:,7])  # [:,6] if blockpref, [:,7] if datapref


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
    ra_deg = Angle(ra, unit=u.hour).degree
    dec_deg = Angle(dec, unit=u.degree).degree

    ras.append(ra_deg)
    decs.append(dec_deg)

T = Table([galaxies, posang, ras, decs, masks, sessions], names = ('Galaxy', 'pa', 'ra', 'dec', 'Mask', 'Session'))
#T.write('galaxy_basics.csv', format='csv')

########################
#### a) Create gallery 
########################

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


    ### add contours and hatches
    if row['Galaxy']=='UGC04136': # ensure a minimum thickness for this galaxy due to its high inclination
         mask_r25 = np.sum(fits.open('masks/'+row['Galaxy']+'_mask_thickened_block.fits')[0].data, axis=0) > 0
    else:
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
        
        beam = Ellipse((center_beamx, center_beamy), width = b1.value/cdelt1.value, height = b2.value/cdelt2.value, angle = bpa.value,
                        edgecolor = 'black', facecolor = 'black', alpha = 0.9
                        )
        
        ax.add_patch(beam)

    labely = row['Galaxy'].upper()

    ### text position: do it in percent!
    if row['Galaxy'] in ['NGC0932', 'NGC3406NED01', 'NGC5216', 'NGC5631', 'UGC02222', 'UGC03960', 'UGC08234', 'UGC09629']:
        ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='r', fontweight = 'bold')
    elif row['Galaxy'] in ['NGC3106', 'NGC3619', 'NGC5157', 'NGC6154', 'NGC6338', 'UGC04136', 'UGC08322', 'UGC10097', 'UGC10905', 'IC0674', 'IC3598']:
        ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='lime', fontweight = 'bold')
    else:
        ax.text(posx1_check + 0.05*(posx2_check-posx1_check), posy1_check+0.88*(posy2_check-posy1_check), s= labely, color='w', fontweight = 'bold') #, bbox = dict(boxstyle='square', facecolor='white', alpha=0.7)


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
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Dummy array is needed for proper initialization
cbar = fig.colorbar(sm, ax=axleg, orientation='horizontal')
cbar.ax.set_xscale('linear')

ell = Ellipse((10,60), width = 10, height = 10, angle = (0),
                edgecolor = 'c', facecolor = 'none', alpha = 1, lw=2.5)
axleg.add_patch(ell)

axleg.text(25,57, s= r'R$_{25} \times$FoV', fontweight='bold', fontsize=12)

ell2 = Ellipse((10,30), width = 10, height = 10, angle = (0),
                edgecolor = 'silver', facecolor = 'none', alpha = 1, lw=2.5)
axleg.add_patch(ell2)

axleg.text(25, 27, s= 'S/N > 5', fontweight='bold', fontsize=12)


#####################
### save the plot
savename = os.path.join(plotdir, 'masked_mom0_gallery')

if do_other_geometry:
    savename += '_geom'
if do_norm:
    savename += '_sqrt'

fig.savefig(savename + '_datapref_fix62_noR25_max50_color_hatch.pdf', dpi = 300)
