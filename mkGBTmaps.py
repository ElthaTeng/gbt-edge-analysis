import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp
from astropy.stats import mad_std

import matplotlib.pyplot as plt
from astropy.wcs import WCS
from spectral_cube import SpectralCube
from astropy.convolution import Gaussian1DKernel
from astropy import units as u
import radio_beam

def get_mask(source, method, session, version, broad=False, write_fits=False):
    
    hdu_data = fits.open(get_pkg_data_filename('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'.fits'))[0]
    new_header = hdu_data.header.copy() 
    data = hdu_data.data

    rms_map = mad_std(data, axis = 0, ignore_nan = True)
    rms_cube = np.repeat(rms_map[np.newaxis, :, :], data.shape[0], axis=0)
    data_noise = np.copy(data)
    data_noise[data > 3 * rms_cube] = np.nan
        
    rms_map_new = mad_std(data_noise, axis = 0, ignore_nan = True)

    if method=='datacube':

        #start with a mask covering all > 2.5 sigma voxels
        rms_cube_new = np.repeat(rms_map_new[np.newaxis, :, :], data.shape[0], axis=0)
        
        mask = np.array(data > 2.5 * rms_cube_new, dtype = int) if broad else np.array(data > 3 * rms_cube_new, dtype = int)  

        #remove spikes
        mask = mask & (np.roll(mask, 1, 0)|np.roll(mask, -1, 0))
        
        new_mask = np.copy(mask)

        if write_fits:
            # rms map for a single channel (not for integrated map)
            fits.writeto('maps/'+source+'_12CO_rmsmap_se'+session+'.fits', rms_map_new, new_header, overwrite=True)
            # reprojected mask
            fits.writeto('masks/'+source+'_mask_reprojected_'+method+'_se'+session+'.fits', new_mask, new_header, overwrite=True) 
    
    else:
        hdu_mask = fits.open(get_pkg_data_filename('masks/from_matlab/mask_'+source+'_'+method+'_'+version+'.fits'))[0]
        new_mask, footprint = reproject_interp(hdu_mask, new_header, order='nearest-neighbor') 
        new_mask[np.isnan(new_mask)] = 0

        if write_fits:
            fits.writeto('masks/'+source+'_mask_reprojected_'+method+'.fits', new_mask, new_header, overwrite=True)  #'_se'+session+
    
    # Create mom0 error maps
    N_masked = np.sum(new_mask, axis=0)
    emom0 = rms_map_new * abs(new_header['CDELT3']) * np.sqrt(N_masked)

    if write_fits:
        # error map for the integrated map
        fits.writeto('maps/'+source+'_12CO_emom0_'+method+'_se'+session+'.fits', emom0, new_header, overwrite=True)

    return new_mask


def expand_mask(source, session, cutoff=0.05, write_fits=False): #only apply this on datacube masks

    #expand mask spatially by convolving to 2D Gaussian with 2x beam size (~16")
    mask = SpectralCube.read('masks/'+source+'_mask_reprojected_datacube_se'+session+'.fits')
    beam = radio_beam.Beam(major=mask.header['BMAJ']*2*u.deg, minor=mask.header['BMIN']*2*u.deg, pa=0*u.deg)
    mask_conv = mask.convolve_to(beam)

    #expand mask spectrally by convolving to 1D Gaussian with std = 2 pixels
    spectral_kernel = Gaussian1DKernel(stddev=2)
    mask_conv_spec = mask_conv.spectral_smooth(spectral_kernel)

    new_mask = mask_conv_spec.hdu.data
    new_mask[new_mask <= cutoff] = 0
    new_mask[new_mask > cutoff] = 1

    new_header = mask.header.copy()
    new_header['BMAJ'] = mask.header['BMAJ'] * 2
    new_header['BMIN'] = mask.header['BMIN'] * 2    

    rms_map = fits.open('maps/'+source+'_12CO_rmsmap_se'+session+'.fits')[0].data
    N_masked = np.sum(new_mask, axis=0)
    emom0 = rms_map * abs(new_header['CDELT3']) * np.sqrt(N_masked)

    if write_fits:
        fits.writeto('masks/'+source+'_mask_reprojected_datacube_expand_se'+session+'.fits', new_mask, new_header, overwrite=True)
        fits.writeto('maps/'+source+'_12CO_emom0_datacube_expand_se'+session+'.fits', emom0, new_header, overwrite=True)
    

def apply_mask(source, method, session, write_fits=False):  # method='datacube_expand' if expand_mask was applied on datacube masks
    
    data = fits.open('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'.fits')[0]
    
    if method=='datacube' or method=='datacube_expand': # 
        mask = fits.open('masks/'+source+'_mask_reprojected_'+method+'_se'+session+'.fits')[0].data
        mask_block = fits.open('masks/'+source+'_mask_reprojected_block.fits')[0].data  #_se'+session+'
        block_zeros = np.zeros(len(mask_block), dtype=int)

        # identify the channels where 'block' mask is all zeros
        for i in range(len(mask_block)):
            block_zeros[i] = (mask_block[i]==0).all() 

        # chop off channels outside +/- vmaxg (use 'block' mask as reference)
        mask[np.where(block_zeros)] = np.zeros((np.sum(block_zeros), data.shape[1], data.shape[2]))

        # apply a spatial R25 mask in addition to the datacube mask
        mask_R25 = np.sum(mask_block, axis=0)
        mask_R25[mask_R25 > 0] = 1
        mask = mask * mask_R25

    else:
        mask = fits.open('masks/'+source+'_mask_reprojected_'+method+'.fits')[0].data  #'_se'+session+

    data_masked = data.data * mask

    if write_fits:
        fits.writeto('data/GBT-masked/'+source+'_12CO_masked_'+method+'_se'+session+'.fits', data_masked, data.header, overwrite=True)  #_noR25
    
    return data_masked


def get_maps(source, method, session, write_fits=False):

    cube = SpectralCube.read('data/GBT-masked/'+source+'_12CO_masked_'+method+'_se'+session+'.fits').with_spectral_unit(u.km/u.s)  #_noR25
    moment_0 = cube.moment(order=0)
    moment_1 = cube.moment(order=1)  
    sigma_map = cube.linewidth_sigma() 

    if write_fits:
        moment_0.write('maps/'+source+'_12CO_mom0_'+method+'_se'+session+'.fits', overwrite=True)   #_noR25
        moment_1.write('maps/'+source+'_12CO_mom1_'+method+'_se'+session+'.fits', overwrite=True)   #_noR25
        sigma_map.write('maps/'+source+'_12CO_vdisp_'+method+'_se'+session+'.fits', overwrite=True)   #_noR25

    return moment_0


def compare_mom0(source, session, save_fig=False, interactive=True):
     
    methods = ['Havfield', 'rotnoHa+', 'rotnoHa-', 'flat', 'block', 'datacube_expand']  #'datacube'
    
    fig = plt.figure(figsize=(18,10))

    for i, method in enumerate(methods):

        fits_map = fits.open('maps/'+source+'_12CO_mom0_'+method+'_se'+session+'.fits')

        mom0 = fits_map[0].data
        wcs = WCS(fits_map[0].header)

        ax = fig.add_subplot(2,3,i+1, projection=wcs)  
        ra = ax.coords[0]
        ra.set_major_formatter('hh:mm:ss.s')
  
        plt.imshow(mom0, origin='lower', cmap='hot', vmin=0) 
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=14)

        zeros = np.ma.masked_where(~(mom0 == 0), mom0)
        plt.imshow(zeros, origin='lower', cmap='Pastel2_r')
        plt.annotate('Peak: '+str(np.round(np.nanmax(mom0), 1))+' (K km/s)', (3, 0.91*mom0.shape[0]), xycoords='data', fontsize=12)  
        plt.annotate('Total: '+str(np.round(np.nansum(mom0), 1))+' (K km/s)', (3, 0.85*mom0.shape[0]), xycoords='data', fontsize=12)  

        plt.tick_params(axis="x", labelsize=14)
        plt.tick_params(axis="y", labelsize=14)

        if i==0:
            title = plt.title(method, fontsize=18)
        else:
            plt.title(method, fontsize=16)

        if i//3 == 0:
            plt.xlabel(' ')
            plt.tick_params(axis="x", labelbottom=False)
        else:
            plt.xlabel('R.A. (J2000)', fontsize=16)
        if i%3 == 0:
            plt.ylabel('Decl. (J2000)', fontsize=16)
        else:
            plt.ylabel(' ')
            plt.tick_params(axis="y", labelleft=False)

    plt.subplots_adjust(hspace=0.1)
    plt.subplots_adjust(wspace=0.1)
    if save_fig:
        plt.savefig('plots/'+source+'_compare_mom0_se'+session+'.pdf', bbox_extra_artists=(title,), bbox_inches='tight', pad_inches=0.02)
    
    if interactive:
        plt.show()




