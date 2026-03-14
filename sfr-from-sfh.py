import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS

table_gbt = np.genfromtxt('galaxy_parameters_erik_v2_fix62.csv', delimiter=',', skip_header=1, dtype='str')

idx_MS = table_gbt[:,14] == 'MS'
idx_GV = table_gbt[:,14] == 'GV'
idx_RG = table_gbt[:,14] == 'RG'

gals = table_gbt[:,0]#[idx_GV]
SFR_list = np.full((gals.shape[0],1), np.nan)

for i, gal in enumerate(list(gals)):

    pipe3d = fits.open(gal+'.Pipe3D.cube.fits')
    mask_pipe3d = fits.open('masks/from_matlab/mask_'+gal+'_block_v3.fits')[0]
    mask_2d = np.sum(mask_pipe3d.data, axis=0) > 0

    Vband_img = pipe3d[1].data[0]
    Mst_dustcorr = 10**pipe3d[1].data[19]
    Av_star = pipe3d[1].data[11]

    # Convert luminosity fractions to mass fractions for each age-metallicity bin
    if gal=='NGC2596':  # using old pipe3D
        ages = np.load('data/DR3_age-met-mass2light.npy')[0]
        ML_all = np.load('data/DR3_age-met-mass2light.npy')[2]
        lumfrac_3d = pipe3d[2].data[0:156]
    else:
        ages = np.load('data/eDR_age-met-mass2light.npy')[0]
        ML_all = np.load('data/eDR_age-met-mass2light.npy')[2]
        lumfrac_3d = pipe3d[2].data[0:273]
        #print(lumfrac_3d.shape)

    ML_temp = np.repeat(ML_all[:, np.newaxis], lumfrac_3d.shape[1], axis=1)
    ML_3d = np.repeat(ML_temp[:, :, np.newaxis], lumfrac_3d.shape[2], axis=2)
    mass_age_met_3d = lumfrac_3d * ML_3d * Vband_img
    mass_age_met_3d_dustcorr = mass_age_met_3d * 10**(Av_star/2.5)
    mass_integrated = np.sum(mass_age_met_3d_dustcorr * mask_2d) 

    # Get mass fraction with age <= 33 Myr (or <= 11.5 Myr)
    if gal=='NGC2596':  # using old pipe3D
        #f_Myoung = np.sum(mass_age_met_3d[:48] * mask_2d * 10**(Av_star/2.5)) / mass_integrated  # <= 31.5 Myr
        f_Myoung = np.sum(mass_age_met_3d[:40] * mask_2d * 10**(Av_star/2.5)) / mass_integrated  # <= 10 Myr
    else:
        #f_Myoung = np.sum(mass_age_met_3d[:70] * mask_2d * 10**(Av_star/2.5)) / mass_integrated  # <= 33 Myr 
        f_Myoung = np.sum(mass_age_met_3d[:42] * mask_2d * 10**(Av_star/2.5)) / mass_integrated  # <= 11.5 Myr 
    Mst_tot = np.nansum(Mst_dustcorr * mask_2d) 
    M_young_int = Mst_tot * f_Myoung

    SFR = M_young_int / 12e6  #  33e6  # 
    print(f_Myoung, Mst_tot, SFR)
    SFR_list[i, 0] = SFR

    # M_young map
    # ax = plt.subplot(111, projection=WCS(mask_pipe3d.header)[0])
    # ra = ax.coords[0]
    # ra.set_major_formatter('hh:mm:ss.s')

    # image = M_young * mask_valid
    # image[image==0] = np.nan
    # plt.imshow(image, origin='lower', vmin=0)
    # cb = plt.colorbar()

    # plt.tick_params(axis="y", labelsize=14, labelleft=True)
    # plt.tick_params(axis="x", labelsize=14, labelbottom=True)
    # cb.ax.tick_params(labelsize=14)
    # plt.xlabel('R.A. (J2000)', fontsize=16) 
    # plt.ylabel('Decl. (J2000)', fontsize=16) 
    #plt.savefig('plots/AHa_SFR_maps/'+gal+'_Myoung.pdf', bbox_inches='tight', pad_inches=0.02)

    #print(mask_2d.sum(), mask_valid.sum())
    print(f'{gal}: Total young mass = {M_young_int} Msun. SFR = {SFR} Msun/yr. \n')

np.savetxt('SFR_12Myr_block_dustcorr_Av.csv', SFR_list, delimiter=',')




