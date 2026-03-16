import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from spectral_cube import SpectralCube
from astropy import units as u

source_list = list(np.genfromtxt('galaxy_list.csv', delimiter=',', skip_header=1, dtype='str')[:, 0])
session_list = list(np.genfromtxt('galaxy_list.csv', delimiter=',', skip_header=1, dtype='str')[:, 6]) 

for n, (source, session) in enumerate(zip(source_list, session_list)):

    hdu_data = fits.open('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'.fits')[0]
    mask_Hav = fits.open('masks/'+source+'_mask_reprojected_Havfield_se'+session+'.fits')[0].data  
    mask_expand = fits.open('masks/'+source+'_mask_reprojected_datacube_expand_se'+session+'.fits')[0].data  
    mask_block = fits.open('masks/'+source+'_mask_reprojected_block_se'+session+'.fits')[0].data 
    block_zeros = np.zeros(len(mask_block), dtype=int)

    # identify the channels where 'block' mask is all zeros
    for i in range(len(mask_block)):
        block_zeros[i] = (mask_block[i]==0).all() 

    # chop off channels outside +/- vmaxg (use 'block' mask as reference)
    mask_expand[np.where(block_zeros)] = np.zeros((np.sum(block_zeros), hdu_data.shape[1], hdu_data.shape[2]))

    # apply a spatial R25 mask in addition to the datacube mask  
    mask_R25 = np.sum(mask_block, axis=0)
    mask_R25[mask_R25 > 0] = 1
    mask = mask_expand * mask_R25

    # Create HavORexpand mask based on GBT data grid
    mask_Havexpand = np.logical_or(mask_Hav, mask).astype(int)
    fits.writeto('masks/'+source+'_mask_reprojected_Havexpand_se'+session+'.fits', mask_Havexpand, hdu_data.header, overwrite=True)  #_noR25

    # Apply mask to data
    data_masked = hdu_data.data * mask_Havexpand
    fits.writeto('data/GBT-masked/'+source+'_12CO_masked_Havexpand_se'+session+'.fits', data_masked, hdu_data.header, overwrite=True)  #_noR25

    # Make maps
    cube = SpectralCube.read('data/GBT-masked/'+source+'_12CO_masked_Havexpand_se'+session+'.fits').with_spectral_unit(u.km/u.s)  #_noR25
    moment_0 = cube.moment(order=0)
    moment_1 = cube.moment(order=1)  
    sigma_map = cube.linewidth_sigma() 

    moment_0.write('maps/'+source+'_12CO_mom0_Havexpand_se'+session+'.fits', overwrite=True)   #_noR25
    moment_1.write('maps/'+source+'_12CO_mom1_Havexpand_se'+session+'.fits', overwrite=True)   #_noR25
    sigma_map.write('maps/'+source+'_12CO_vdisp_Havexpand_se'+session+'.fits', overwrite=True)   #_noR25

    # Show and save moment 0 plot
    fits_map = fits.open('maps/'+source+'_12CO_mom0_Havexpand_se'+session+'.fits')  #_noR25
    mom0 = fits_map[0].data
    wcs = WCS(fits_map[0].header)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=wcs)  
    ra = ax.coords[0]
    ra.set_major_formatter('hh:mm:ss.s')

    plt.imshow(mom0, origin='lower', cmap='hot', vmin=0) 
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=14)

    zeros = np.ma.masked_where(~(mom0 == 0), mom0)
    plt.imshow(zeros, origin='lower', cmap='Pastel2_r')
    plt.annotate('Peak: '+str(np.round(np.nanmax(mom0), 1))+' (K km/s)', (3, 0.9*mom0.shape[0]), xycoords='data', fontsize=14)  #(3, 130)  #(3, 88)
    plt.annotate('Total: '+str(np.round(np.nansum(mom0), 1))+' (K km/s)', (3, 0.8*mom0.shape[0]), xycoords='data', fontsize=14)  #(3, 120)  #(3, 80)

    plt.tick_params(axis="x", labelsize=14)
    plt.tick_params(axis="y", labelsize=14)

    title = plt.title('HavORexpand', fontsize=18)
    plt.xlabel('R.A. (J2000)', fontsize=16)
    plt.ylabel('Decl. (J2000)', fontsize=16)

    plt.savefig('plots/'+source+'_mom0_Havexpand_se'+session+'.pdf', bbox_extra_artists=(title,), bbox_inches='tight', pad_inches=0.02) 
    #plt.show()

    # Make emom0 map
    rms_map_new = fits.open('maps/'+source+'_12CO_rmsmap_se'+session+'.fits')[0].data
    new_header = fits.open('maps/'+source+'_12CO_rmsmap_se'+session+'.fits')[0].header
    N_masked = np.sum(mask_Havexpand, axis=0)
    emom0 = rms_map_new * abs(new_header['CDELT3']) * np.sqrt(N_masked)

    fits.writeto('maps/'+source+'_12CO_emom0_Havexpand_se'+session+'.fits', emom0, new_header, overwrite=True)  #_noR25

    print(f'{source} done.')
