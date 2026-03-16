from mkGBTmaps import *
import time
import numpy as np
from astropy.io import fits

def make_2d_gaussian(size_xy, fwhm_xy):

    sigma_xy = fwhm_xy / 2.3548  
    ax = np.arange(-size_xy // 2 + 1., size_xy // 2 + 1.)
    xx, yy = np.meshgrid(ax, ax, indexing='ij')
    kernel = np.exp(-(xx**2 + yy**2) / (2. * sigma_xy**2))

    return kernel / kernel.max()


def make_1d_gaussian(size_v, fwhm_v):

    sigma_v = fwhm_v / 2.3548  
    vv = np.arange(-size_v // 2 + 1., size_v // 2 + 1.)
    kernel = np.exp(-(vv**2) / (2. * sigma_v**2))

    return kernel / kernel.max()


def Mmol(source, mask, session, aco, distance): 
    
    kmom0_fits = fits.open('maps/'+source+'_12CO_mom0_'+mask+'_se'+session+'.fits')
    mom0 = kmom0_fits[0]
    valid = mom0.data > 0  # Exclude negative pixels on mom0 map when computing total flux 
    Abeam = np.pi * (mom0.header['BMAJ']/2) * (mom0.header['BMIN']/2) / np.log(2)  #deg^2
    Npixperbeam = Abeam / abs(mom0.header['CDELT1'] * mom0.header['CDELT2'])
    print(f'# of pixels per beam: {Npixperbeam}')
    
    I_CO = np.nansum(mom0.data[valid]) / Npixperbeam  
    Omega_beam = np.pi * (np.pi/180)**2 * (mom0.header['BMAJ']/2) * (mom0.header['BMIN']/2) / np.log(2)  #rad^2
    L_CO = I_CO * Omega_beam * (1e6*distance)**2
    Mmol_kcube = L_CO * aco
    print(f'Omega_beam = {Omega_beam}')

    emom0 = fits.open('maps/'+source+'_12CO_emom0_'+mask+'_se'+session+'.fits')[0]
    I_CO_err = np.nansum(emom0.data[valid]) / Npixperbeam  #
    Mmol_err = I_CO_err * Omega_beam * (1e6*distance)**2 * aco
    
    print(f'I_tot = {I_CO} +/- {I_CO_err} K km/s') 
    print(f'L_tot = {L_CO} K km/s pc2')
    print(f'M_mol (using K cube) = {Mmol_kcube} +/- {Mmol_err} Msun')

    return I_CO, I_CO_err, L_CO


galaxy_list = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')[:,0]
gbt_session_list = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')[:,5]
create_file = True
broad_datamask = True
session = 'loop100'
num_run = 100
f_peak = 1
f_size = 1

N_gal = galaxy_list.shape[0]
Ico_ref_list = np.full((N_gal, num_run), np.nan)
Ico_list = np.full((N_gal, num_run), np.nan)
eIco_list = np.full((N_gal, num_run), np.nan)
Lco_list = np.full((N_gal, num_run), np.nan)

program_start = time.time() 

for n_run in range(num_run): 

    start_time = time.time() 

    for n_gal, galaxy in enumerate(galaxy_list):

        gbt_session = gbt_session_list[n_gal]
        hdu_data = fits.open('data/'+galaxy+'_12CO_rebase7_smooth1.3_hanning2_se'+gbt_session+'.fits')[0]
        print(f'Run {n_run+1}: {galaxy}')
        print(f'Total flux of original cube: {np.nansum(hdu_data.data)} K km/s')    

        ## Create datacube with added fake source
        # create a PSF: beam size (8.3"), 60 km/s in fwhm, Tpeak nomalized to 1 
        fwhm_xy = f_size * hdu_data.header['BMAJ']/hdu_data.header['CDELT2']  # ~5.2 pixels
        fwhm_v = 4  # 15.2 * 4 ~ 60 km/s 
        size_xy = 13  # ~2x fwhm_xy
        size_v = 9   # ~2x fwhm_
        psf = f_peak * make_2d_gaussian(size_xy, fwhm_xy)
        gauss1d = make_1d_gaussian(size_v, fwhm_v)
        fake_source = psf[None, :, :] * gauss1d[:, None, None]
        print(f'Fake source total flux: {fake_source.sum()} K km/s')

        # randomly select a ppv position within the block_3d mask that ensures full fit of the fake source
        mask_block = fits.open('masks/'+galaxy+'_mask_reprojected_block.fits')[0].data  #_se'+gbt_session+' 
        margin_xy = size_xy // 2
        margin_v = size_v // 2
        valid_mask = np.zeros_like(mask_block, dtype=bool)
        valid_mask[margin_v:-margin_v, margin_xy:-margin_xy, margin_xy:-margin_xy] = True
        eligible = np.argwhere((mask_block>0.5) & valid_mask)
        print(f'Total flux of original cube within block mask: {np.nansum(hdu_data.data * mask_block)} K km/s')

        if eligible.size == 0:
            raise ValueError("No eligible positions found.")
        else:
            i, j, k = eligible[np.random.choice(len(eligible))]
            print(f"Inserted fake source at (i, j, k) = ({i}, {j}, {k})")
        
        # place fake source at the selected ppv position
        v_slice = slice(i - margin_v, i + margin_v + 1)
        y_slice = slice(j - margin_xy, j + margin_xy + 1)
        x_slice = slice(k - margin_xy, k + margin_xy + 1)
        new_cube = hdu_data.data.copy()
        new_cube[v_slice, y_slice, x_slice] += fake_source

        fits.writeto('data/'+galaxy+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'.fits', new_cube, hdu_data.header, overwrite=True) #


        ## mkGBTmaps.py
        get_mask(galaxy, 'datacube', session, version=None, broad=broad_datamask, write_fits=create_file)
        cutoff = 0.05 #if broad_datamask else 0.05  #0.02 #0.01 #
        expand_mask(galaxy, session, cutoff=cutoff, write_fits=create_file)  
        mask_method = 'datacube_expand'

        apply_mask(galaxy, mask_method, session, write_fits=create_file)
        get_maps(galaxy, mask_method, session, write_fits=create_file)


        ## runORmasks.py   
        mask_expand = fits.open('masks/'+galaxy+'_mask_reprojected_datacube_expand_se'+session+'.fits')[0].data 
        mask_Hav = fits.open('masks/'+galaxy+'_mask_reprojected_Havfield.fits')[0].data  # _se'+gbt_session+' 
        block_zeros = np.zeros(len(mask_block), dtype=int)

        # identify the channels where 'block' mask is all zeros
        for idx in range(len(mask_block)):
            block_zeros[idx] = (mask_block[idx]==0).all() 

        # chop off channels outside +/- vmaxg (use 'block' mask as reference)
        mask_expand[np.where(block_zeros)] = np.zeros((np.sum(block_zeros), hdu_data.shape[1], hdu_data.shape[2]))

        # apply a spatial R25 mask in addition to the datacube mask  
        mask_R25 = np.sum(mask_block, axis=0)
        mask_R25[mask_R25 > 0] = 1
        mask = mask_expand * mask_R25

        # Create HavORexpand mask based on GBT data grid
        mask_Havexpand = np.logical_or(mask_Hav, mask).astype(int)
        fits.writeto('masks/'+galaxy+'_mask_reprojected_Havexpand_se'+session+'.fits', mask_Havexpand, hdu_data.header, overwrite=True)  #_noR25

        # Apply mask to data
        data_masked = hdu_data.data * mask_Havexpand
        fits.writeto('data/GBT-masked/'+galaxy+'_12CO_masked_Havexpand_se'+session+'.fits', data_masked, hdu_data.header, overwrite=True)  #_noR25

        # Make maps
        cube = SpectralCube.read('data/GBT-masked/'+galaxy+'_12CO_masked_Havexpand_se'+session+'.fits').with_spectral_unit(u.km/u.s)  #_noR25
        moment_0 = cube.moment(order=0)
        moment_1 = cube.moment(order=1)  
        sigma_map = cube.linewidth_sigma() 
        moment_0.write('maps/'+galaxy+'_12CO_mom0_Havexpand_se'+session+'.fits', overwrite=True)
        # emom0 map
        rms_map_new = fits.open('maps/'+galaxy+'_12CO_rmsmap_se'+session+'.fits')[0].data
        new_header = fits.open('maps/'+galaxy+'_12CO_rmsmap_se'+session+'.fits')[0].header
        N_masked = np.sum(mask_Havexpand, axis=0)
        emom0 = rms_map_new * abs(new_header['CDELT3']) * np.sqrt(N_masked)
        fits.writeto('maps/'+galaxy+'_12CO_emom0_Havexpand_se'+session+'.fits', emom0, new_header, overwrite=True)  #_noR25


        ## Get theoretical and measured integrated fluxes 
        # Compute theoretical line flux
        kmom0_fits = fits.open('maps/'+galaxy+'_12CO_mom0_Havexpand_se'+gbt_session+'.fits')
        mom0 = kmom0_fits[0]
        valid = mom0.data > 0  # Exclude negative pixels on mom0 map when computing total flux 
        Abeam = np.pi * (mom0.header['BMAJ']/2) * (mom0.header['BMIN']/2) / np.log(2)  #deg^2
        Npixperbeam = Abeam / abs(mom0.header['CDELT1'] * mom0.header['CDELT2'])
        Ico_ideal = (np.nansum(mom0.data[valid]) + fake_source.sum()) / Npixperbeam 
        Ico_ref_list[n_gal, n_run] = Ico_ideal

        # Compute measured line flux
        distance = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1)[n_gal, 24]
        Ico, eIco, Lco = Mmol(galaxy, 'Havexpand', session, 4.35, distance)
        Ico_list[n_gal, n_run] = Ico
        eIco_list[n_gal, n_run] = eIco
        Lco_list[n_gal, n_run] = Lco

    print(f'Run {n_run+1} done; took {np.round(time.time() - start_time, 1)} sec.\n')

# Save measured total flux (incl. fake source) for each run and each galaxy
np.save('Ico_list_allruns_'+session+'.npy', Ico_list)

mean_Ico_ref = np.mean(Ico_ref_list, axis=1)
mean_Ico = np.mean(Ico_list, axis=1)
mean_eIco = np.mean(eIco_list, axis=1)
mean_Lco = np.mean(Lco_list, axis=1)

# Save fake source fluxes, mean of the measured fluxes, and errors
np.savetxt('Ico_ref_list_'+session+'.csv', mean_Ico_ref.reshape(len(galaxy_list),1), delimiter=',')
np.savetxt('Ico_list_'+session+'.csv', mean_Ico.reshape(len(galaxy_list),1), delimiter=',')
np.savetxt('Ico_err_list_'+session+'.csv', mean_eIco.reshape(len(galaxy_list),1), delimiter=',')
np.savetxt('Lco_list_'+session+'.csv', mean_Lco.reshape(len(galaxy_list),1), delimiter=',')

print(f'Total elapsed time: {np.round(time.time() - program_start, 1)} sec.\n')

