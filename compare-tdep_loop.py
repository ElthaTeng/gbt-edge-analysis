import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import warnings


### SET THE ALPHA_CO PRESCRIPTION: 'const', 'B13_Zgrad', 'SL24_Zgrad', or 'T24' 
aco_choice = 'const'  


def Mmol(source, mask, session, aco): 
    
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
    
    tdep_kcube = Mmol_kcube / SFR_tot * 1e-9
    print(f't_dep (using K cube) = {tdep_kcube} Gyr \n')
    etdep_kcube = Mmol_err / SFR_tot * 1e-9

    return Mmol_kcube, Mmol_err, tdep_kcube, etdep_kcube, L_CO


def Mmol_B13(source, mask, session, metallicity=False, gradient=False):
    pipe3d = fits.open(source+'.Pipe3D.cube.fits')
    Mst_dustcorr = 10**pipe3d[1].data[19]
    Sigma_st = Mst_dustcorr / (1e6*distance/206265)**2

    if metallicity:
        if gradient:
            radii_map = np.load('maps/radii_maps/'+source+'_radii_map.npy')
            Zgrad_fit = np.load('bestfit_Zprime_gradient_Curti_fix62.npy')
            gal_id = np.where(table_gbt[:,0]==source)[0][0]
            slope = Zgrad_fit[gal_id, 0]
            offset = Zgrad_fit[gal_id, 1]
            Zprime = 10**(radii_map * slope + offset)

            plt.imshow(Zprime, origin='lower')
            plt.colorbar()
            plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_Zprime_gradfit.png', bbox_inches='tight', pad_inches=0.02)
            plt.clf()

            np.save('plots/fit_Zprime_gradient/npys/'+source+'_Zprime_B13_gradfit.npy', Zprime)
       
        else:
            flux_elines = pipe3d[3] if source=='NGC2596' else pipe3d[5]
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

            OH_u = 8.424 + 0.03*np.log10(R3/S2) + 0.751*np.log10(N2) + (-0.349+0.182*np.log10(R3/S2) + 0.508*np.log10(N2))*np.log10(S2)
            OH_l = 8.072 + 0.789*np.log10(R3/S2) + 0.726*np.log10(N2) + (1.069-0.17*np.log10(R3/S2) + 0.022*np.log10(N2))*np.log10(S2)

            Zmap = np.full(N2.shape, np.nan)
            Zmap[np.log10(N2) >= -0.6] = OH_u[np.log10(N2) >= -0.6]
            Zmap[np.log10(N2) < -0.6] = OH_l[np.log10(N2) < -0.6]
            Zprime = 10**(Zmap - 8.69)

            plt.imshow(Zprime, origin='lower', vmin=0.8, vmax=1.2)
            plt.colorbar()
            plt.savefig('plots/fit_Zprime_gradient/maps'+source+'_Zprime_pg16s.png', bbox_inches='tight', pad_inches=0.02)
            plt.clf()

            np.save('plots/fit_Zprime_gradient/npys/'+source+'_Zprime_B13_pg16s.npy', Zprime)

    else:
        if gradient:
            msg = 'gradient only works when metallicity=True. Proceeding with metallicity=False and gradient=False.'
            warnings.warn(msg, UserWarning)
        Zprime = 1

    aco_B13 = np.full(Sigma_st.shape, np.nan)
    Z_term = 2.9 * np.exp(0.4/Zprime)
    aco_B13 = Z_term * (Sigma_st/100)**-0.5 
    aco_B13[Sigma_st < 100] = Z_term[Sigma_st < 100] if metallicity else Z_term

    plt.imshow(aco_B13, origin='lower')  
    plt.colorbar()
    if metallicity:
        if gradient:
            plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_aco_B13_Zgradfit.png', bbox_inches='tight', pad_inches=0.02)
        else:
            plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_aco_B13_pg16s.png', bbox_inches='tight', pad_inches=0.02)
    else:
        plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_aco_B13.png', bbox_inches='tight', pad_inches=0.02)
    plt.clf()
    
    kmom0_fits = fits.open('maps/'+source+'_12CO_mom0_'+mask+'_se'+session+'.fits')
    mom0 = kmom0_fits[0].data
    gbt_header = kmom0_fits[0].header

    # Create pipe3d's 2D map header for regridding aco map onto GBT grid
    mask_pipe3d = fits.open('masks/from_matlab/mask_'+source+'_block_v3.fits')[0]
    pipe3d_header_2d = mask_pipe3d.header
    pipe3d_header_2d['NAXIS'] = 2
    keywords_to_remove = ['NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CDELT3', 'CUNIT3']
    for key in keywords_to_remove:
        del pipe3d_header_2d[key]

    aco_B13_regrid, footprint = reproject_interp((aco_B13, pipe3d_header_2d), gbt_header, order='nearest-neighbor')

    I_CO = np.copy(mom0)
    I_CO[mom0 <= 0] = np.nan  # Exclude negative pixels on mom0 map when computing total flux

    Abeam = np.pi * (kmom0_fits[0].header['BMAJ']/2) * (kmom0_fits[0].header['BMIN']/2) / np.log(2)  #deg^2
    Npixperbeam = Abeam / abs(kmom0_fits[0].header['CDELT1'] * kmom0_fits[0].header['CDELT2'])
    print(f'# of pixels per beam: {Npixperbeam}')

    Nvalid_mom0 = np.sum(np.isfinite(I_CO))
    Nvalid_aco = np.sum(np.isfinite(I_CO) * np.isfinite(aco_B13_regrid))
    completeness = Nvalid_aco / Nvalid_mom0
    print(Nvalid_mom0, Nvalid_aco)

    Omega_beam = (np.pi/180)**2 * Abeam  #rad^2
    L_CO = I_CO * Omega_beam * (1e6*distance)**2
    Mmol_kcube = L_CO * aco_B13_regrid
    print(f'Omega_beam = {Omega_beam}')

    emom0 = fits.open('maps/'+source+'_12CO_emom0_'+mask+'_se'+session+'.fits')[0]
    I_CO_err = np.copy(emom0.data)
    I_CO_err[mom0 <= 0] = np.nan
    Mmol_err = I_CO_err * Omega_beam * (1e6*distance)**2 * aco_B13_regrid
    
    if metallicity and not gradient:
        Mmol_tot = np.nansum(Mmol_kcube) / (Npixperbeam * completeness)
        Mmol_err_tot = np.nansum(Mmol_err) / (Npixperbeam * completeness)
    else:
        Mmol_tot = np.nansum(Mmol_kcube) / Npixperbeam 
        Mmol_err_tot = np.nansum(Mmol_err) / Npixperbeam       
    
    print(f'I_tot = {I_CO} +/- {I_CO_err} K km/s')  
    print(f'L_tot = {L_CO} K km/s pc2')
    print(f'M_mol (using K cube) = {Mmol_tot} +/- {Mmol_err_tot} Msun')
    
    tdep_tot = Mmol_tot / SFR_tot * 1e-9
    etdep_tot = Mmol_err_tot / SFR_tot * 1e-9
    print(f't_dep (using K cube) = {tdep_tot} +/- {etdep_tot} Gyr \n')

    return Mmol_tot, Mmol_err_tot, tdep_tot, etdep_tot, completeness


def Mmol_SL24(source, mask, session, gradient=False):
    pipe3d = fits.open(source+'.Pipe3D.cube.fits')
    Mst_dustcorr = 10**pipe3d[1].data[19]
    Sigma_st = Mst_dustcorr / (1e6*distance/206265)**2

    if gradient:
        radii_map = np.load('maps/radii_maps/'+source+'_radii_map.npy')
        Zgrad_fit = np.load('bestfit_Zprime_gradient_Curti_fix62.npy')
        gal_id = np.where(table_gbt[:,0]==source)[0][0]
        slope = Zgrad_fit[gal_id, 0]
        offset = Zgrad_fit[gal_id, 1]
        Zprime = 10**(radii_map * slope + offset)

        plt.imshow(Zprime, origin='lower')
        plt.colorbar()
        plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_Zprime_gradfit.png', bbox_inches='tight', pad_inches=0.02)
        plt.clf()
    else:
        flux_elines = pipe3d[3] if source=='NGC2596' else pipe3d[5]
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

        OH_u = 8.424 + 0.03*np.log10(R3/S2) + 0.751*np.log10(N2) + (-0.349+0.182*np.log10(R3/S2) + 0.508*np.log10(N2))*np.log10(S2)
        OH_l = 8.072 + 0.789*np.log10(R3/S2) + 0.726*np.log10(N2) + (1.069-0.17*np.log10(R3/S2) + 0.022*np.log10(N2))*np.log10(S2)

        Zmap = np.full(N2.shape, np.nan)
        Zmap[np.log10(N2) >= -0.6] = OH_u[np.log10(N2) >= -0.6]
        Zmap[np.log10(N2) < -0.6] = OH_l[np.log10(N2) < -0.6]
        Zprime = 10**(Zmap - 8.69)

        plt.imshow(Zprime, origin='lower', vmin=0.8, vmax=1.2)
        plt.colorbar()
        plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_Zprime_pg16s.png', bbox_inches='tight', pad_inches=0.02)
        plt.clf()

    Z_term = 4.35 * Zprime**-1.5
    aco_SL24 = Z_term * (Sigma_st/100)**-0.25 
    aco_SL24[Sigma_st < 100] = Z_term[Sigma_st < 100] 

    plt.imshow(aco_SL24, origin='lower')  
    plt.colorbar()
    if gradient:
        plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_aco_SL24_Zgradfit.png', bbox_inches='tight', pad_inches=0.02)
    else:
        plt.savefig('plots/fit_Zprime_gradient/maps/'+source+'_aco_SL24.png', bbox_inches='tight', pad_inches=0.02)
    plt.clf()
    
    kmom0_fits = fits.open('maps/'+source+'_12CO_mom0_'+mask+'_se'+session+'.fits')
    mom0 = kmom0_fits[0].data
    gbt_header = kmom0_fits[0].header

    # Create GBT 2D map header for regridding aco map onto GBT grid
    mask_pipe3d = fits.open('masks/from_matlab/mask_'+source+'_block_v3.fits')[0]
    pipe3d_header_2d = mask_pipe3d.header
    pipe3d_header_2d['NAXIS'] = 2
    keywords_to_remove = ['NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CDELT3', 'CUNIT3']
    for key in keywords_to_remove:
        del pipe3d_header_2d[key]

    aco_SL24_regrid, footprint = reproject_interp((aco_SL24, pipe3d_header_2d), gbt_header, order='nearest-neighbor')

    I_CO = np.copy(mom0)
    I_CO[mom0 <= 0] = np.nan  # Exclude negative pixels on mom0 map when computing total flux

    Abeam = np.pi * (kmom0_fits[0].header['BMAJ']/2) * (kmom0_fits[0].header['BMIN']/2) / np.log(2)  #deg^2
    Npixperbeam = Abeam / abs(kmom0_fits[0].header['CDELT1'] * kmom0_fits[0].header['CDELT2'])
    print(f'# of pixels per beam: {Npixperbeam}')

    Nvalid_mom0 = np.sum(np.isfinite(I_CO))
    Nvalid_aco = np.sum(np.isfinite(I_CO) * np.isfinite(aco_SL24_regrid))
    completeness = Nvalid_aco / Nvalid_mom0
    print(Nvalid_mom0, Nvalid_aco)

    Omega_beam = (np.pi/180)**2 * Abeam  #rad^2
    L_CO = I_CO * Omega_beam * (1e6*distance)**2
    Mmol_kcube = L_CO * aco_SL24_regrid
    print(f'Omega_beam = {Omega_beam}')

    emom0 = fits.open('maps/'+source+'_12CO_emom0_'+mask+'_se'+session+'.fits')[0]
    I_CO_err = np.copy(emom0.data)
    I_CO_err[mom0 <= 0] = np.nan
    Mmol_err = I_CO_err * Omega_beam * (1e6*distance)**2 * aco_SL24_regrid
    
    if gradient:
        Mmol_tot = np.nansum(Mmol_kcube) / Npixperbeam 
        Mmol_err_tot = np.nansum(Mmol_err) / Npixperbeam
    else:
        Mmol_tot = np.nansum(Mmol_kcube) / (Npixperbeam * completeness)
        Mmol_err_tot = np.nansum(Mmol_err) / (Npixperbeam * completeness)
    
    print(f'I_tot = {I_CO} +/- {I_CO_err} K km/s')  
    print(f'L_tot = {L_CO} K km/s pc2')
    print(f'M_mol (using K cube) = {Mmol_tot} +/- {Mmol_err_tot} Msun')
    
    tdep_tot = Mmol_tot / SFR_tot * 1e-9
    etdep_tot = Mmol_err_tot / SFR_tot * 1e-9
    print(f't_dep (using K cube) = {tdep_tot} +/- {etdep_tot} Gyr \n')

    return Mmol_tot, Mmol_err_tot, tdep_tot, etdep_tot, completeness


def Mmol_T24(source, mask, session):
    
    vdisp_raw = fits.open('maps/'+source+'_12CO_vdisp_'+mask+'_se'+session+'.fits')[0].data
    vdisp = vdisp_raw * np.sqrt(np.cos(inclination*np.pi/180))
    vdisp[vdisp<=0.1] = np.nan
    kmom0_fits = fits.open('maps/'+source+'_12CO_mom0_'+mask+'_se'+session+'.fits')
    mom0 = kmom0_fits[0].data

    aco_T24 = 10**(-0.96 * np.log10(vdisp) + 1.77)  
    I_CO = np.copy(mom0)
    I_CO[mom0 <= 0] = np.nan  # Mask out zero and negative pixels

    Nvalid_mom0 = np.sum(np.isfinite(I_CO))
    Nvalid_vdisp = np.sum(np.isfinite(I_CO) * np.isfinite(vdisp))
    completeness = Nvalid_vdisp / Nvalid_mom0
    print(Nvalid_mom0, Nvalid_vdisp)

    Abeam = np.pi * (kmom0_fits[0].header['BMAJ']/2) * (kmom0_fits[0].header['BMIN']/2) / np.log(2)  #deg^2
    Npixperbeam = Abeam / abs(kmom0_fits[0].header['CDELT1'] * kmom0_fits[0].header['CDELT2'])
    print(f'# of pixels per beam: {Npixperbeam}')

    Omega_beam = (np.pi/180)**2 * Abeam  #rad^2
    L_CO = I_CO * Omega_beam * (1e6*distance)**2
    Mmol_kcube = L_CO * aco_T24
    print(f'Omega_beam = {Omega_beam}')

    emom0 = fits.open('maps/'+source+'_12CO_emom0_'+mask+'_se'+session+'.fits')[0]
    I_CO_err = np.copy(emom0.data)
    I_CO_err[mom0 < 0] = 0
    Mmol_err = I_CO_err * Omega_beam * (1e6*distance)**2 * aco_T24
    
    Mmol_tot = np.nansum(Mmol_kcube) / (Npixperbeam * completeness)
    Mmol_err_tot = np.nansum(Mmol_err) / (Npixperbeam * completeness)
    
    print(f'I_tot = {I_CO} +/- {I_CO_err} K km/s') 
    print(f'L_tot = {L_CO} K km/s pc2')
    print(f'M_mol (using K cube) = {Mmol_tot} +/- {Mmol_err_tot} Msun')
    
    tdep_tot = Mmol_tot / SFR_tot * 1e-9
    etdep_tot = Mmol_err_tot / SFR_tot * 1e-9
    print(f't_dep (using K cube) = {tdep_tot} +/- {etdep_tot} Gyr \n')

    return Mmol_tot, Mmol_err_tot, tdep_tot, etdep_tot, completeness


table_gbt = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')
source_list = list(table_gbt[:, 0]) 
session_list = list(table_gbt[:, 5])
method_list = list(table_gbt[:, 7]) 

catalog = np.genfromtxt('edge_califa.csv', delimiter=',', skip_header=61, dtype='str')
inclinations = np.array(table_gbt[:, 3], dtype='float')

N_gal = np.array(source_list).shape[0]
Mmol_list = np.full((N_gal,), np.nan)
eMmol_list = np.full((N_gal,), np.nan)
tdep_list = np.full((N_gal,), np.nan)
etdep_list = np.full((N_gal,), np.nan)
complete_list = np.full((N_gal,), np.nan)

for n, (source, session, method) in enumerate(zip(source_list, session_list, method_list)):

    mask_pipe3d = fits.open('masks/from_matlab/mask_'+source+'_block_v3.fits')[0]       

    coord_ref = fits.open('masks/from_matlab/mask_'+source+'_Havfield_v1.fits')[0]

    inclination = inclinations[n]

    # Get Ha, Hb measurements from Pipe3D "FLUX_ELINES"
    elineflux = fits.open(source+'.Pipe3D.cube.fits')[3] if source== 'NGC2596' else fits.open(source+'.Pipe3D.cube.fits')[5]
    F_Ha = elineflux.data[45]
    F_Ha_err = elineflux.data[249] if source== 'NGC2596' else elineflux.data[261]
    EW_Ha = elineflux.data[198] if source== 'NGC2596' else elineflux.data[207]
    v_Ha = elineflux.data[96] if source== 'NGC2596' else elineflux.data[99]
    v_Ha_err = elineflux.data[300] if source== 'NGC2596' else elineflux.data[315]
    F_Hb = elineflux.data[28]
    F_Hb_err = elineflux.data[232] if source== 'NGC2596' else elineflux.data[244]

    nonsf = abs(EW_Ha) < 6
    snr = F_Ha / F_Ha_err
    snrcut = snr < 5
    F_Ha[nonsf] = np.nan
    F_Hb[nonsf] = np.nan
    v_Ha[nonsf] = np.nan

    pc2cm = 3.08567758128e18
    distance = catalog[catalog[:,1]==source][0,34].astype('float')  #c * z / H0  #Mpc
    A_sphere = 4 * np.pi * (distance * 1e6 * pc2cm)**2  #cm2

    A_Ha = 5.86 * np.log10(F_Ha / F_Hb / 2.86)
    A_Ha[A_Ha > 3] = 3
    SFR = 7.9e-42 * F_Ha * 10**(A_Ha/2.5) * 1e-16 * A_sphere

    # A_Ha map
    ax = plt.subplot(111, projection=WCS(coord_ref.header)[0])
    ra = ax.coords[0]
    ra.set_major_formatter('hh:mm:ss.s')

    plt.imshow(A_Ha, origin='lower', vmin=0, vmax=5)
    cb = plt.colorbar()

    plt.tick_params(axis="y", labelsize=14, labelleft=True)
    plt.tick_params(axis="x", labelsize=14, labelbottom=True)
    cb.ax.tick_params(labelsize=14)
    plt.xlabel('R.A. (J2000)', fontsize=16) 
    plt.ylabel('Decl. (J2000)', fontsize=16) 
    plt.savefig('plots/AHa_SFR_maps/'+source+'_AHa_block_se'+session+'.pdf', bbox_inches='tight', pad_inches=0.02)
    plt.clf()

    # SFR map
    ax = plt.subplot(111, projection=WCS(coord_ref.header)[0])
    ra = ax.coords[0]
    ra.set_major_formatter('hh:mm:ss.s')

    image = SFR * (np.sum(mask_pipe3d.data, axis=0) > 0)
    image[image==0] = np.nan
    plt.imshow(image, origin='lower', vmin=0)
    cb = plt.colorbar()

    plt.tick_params(axis="y", labelsize=14, labelleft=True)
    plt.tick_params(axis="x", labelsize=14, labelbottom=True)
    cb.ax.tick_params(labelsize=14)
    plt.xlabel('R.A. (J2000)', fontsize=16) 
    plt.ylabel('Decl. (J2000)', fontsize=16) 
    plt.savefig('plots/AHa_SFR_maps/'+source+'_SFR_block_se'+session+'.pdf', bbox_inches='tight', pad_inches=0.02)
    plt.clf()

    print(source + ':')
    print('Distance:', distance, 'Mpc')

    SFR_tot = np.nansum(image)
    print(f'SFR = {SFR_tot} Msun/yr')

    # Run M_mol, t_dep calculations
    if aco_choice == 'const':
        M_mol, M_mol_err, t_dep, t_dep_err, completeness = Mmol(source, method, session, 4.35) 
    elif aco_choice == 'B13_Zgrad':
        M_mol, M_mol_err, t_dep, t_dep_err, completeness = Mmol_B13(source, method, session, metallicity=True, gradient=True)
    elif aco_choice == 'SL24_Zgrad':
        M_mol, M_mol_err, t_dep, t_dep_err, completeness = Mmol_SL24(source, method, session, gradient=True)
    elif aco_choice == 'T24':
        M_mol, M_mol_err, t_dep, t_dep_err, completeness = Mmol_T24(source, method, session)
    else:
        print('Unrecognized choice of aco. Please ensure entering one of the following in aco_choice : const, B13, B13_Zgrad, SL24_Zgrad, T24.')

    Mmol_list[n] = M_mol
    eMmol_list[n] = M_mol_err
    tdep_list[n] = t_dep
    etdep_list[n] = t_dep_err
    complete_list[n] = completeness

# Save each output as an individual npy array and csv list
np.save('Mmol_galaxy_list_sfrblk_datapref_'+aco_choice+'.npy', Mmol_list)
np.save('Mmol_err_galaxy_list_sfrblk_datapref_'+aco_choice+'.npy', eMmol_list)
np.save('tdep_galaxy_list_sfrblk_datapref_'+aco_choice+'.npy', tdep_list)
np.save('tdep_err_galaxy_list_sfrblk_datapref_'+aco_choice+'.npy', tdep_list)

np.savetxt('Mmol_galaxy_list_sfrblk_datapref_'+aco_choice+'.csv', Mmol_list.reshape(len(source_list),1), delimiter=',')
np.savetxt('Mmol_err_galaxy_list_sfrblk_datapref_'+aco_choice+'.csv', eMmol_list.reshape(len(source_list),1), delimiter=',')
np.savetxt('tdep_galaxy_list_sfrblk_datapref_'+aco_choice+'.csv', tdep_list.reshape(len(source_list),1), delimiter=',')
np.savetxt('tdep_err_galaxy_list_sfrblk_datapref_'+aco_choice+'.csv', etdep_list.reshape(len(source_list),1), delimiter=',')

if aco_choice == 'T24':
    np.save('completeness_galaxy_list_sfrblk_'+aco_choice+'.npy', complete_list)
    np.savetxt('completeness_galaxy_list_sfrblk_'+aco_choice+'.csv', complete_list.reshape(len(source_list),1), delimiter=',')
elif aco_choice == 'const':
    np.save('Lco_galaxy_list_sfrblk_datapref.npy', complete_list)
    np.savetxt('Lco_galaxy_list_sfrblk_datapref.csv', complete_list.reshape(len(source_list),1), delimiter=',')
