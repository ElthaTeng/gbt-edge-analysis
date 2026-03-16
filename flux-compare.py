import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from Jy_K_convertor import *
from spectral_cube import SpectralCube
import radio_beam

def Mmol(source, mask, session, aco): 
    
    dir = 'acacompare_outputs/' if session=='7msmo12' else './'
    kmom0_fits = fits.open(dir+'maps/'+source+'_12CO_mom0_'+mask+'_se'+session+'.fits')
    mom0 = kmom0_fits[0]
    valid = mom0.data > 0  # Exclude negative pixels on mom0 map when computing total flux 
    Abeam = np.pi * (mom0.header['BMAJ']/2) * (mom0.header['BMIN']/2) / np.log(2)  #deg^2
    Npixperbeam = Abeam / abs(mom0.header['CDELT1'] * mom0.header['CDELT2'])
    #print(f'# of pixels per beam: {Npixperbeam}')
    
    I_CO = np.nansum(mom0.data[valid]) / Npixperbeam  
    Omega_beam = np.pi * (np.pi/180)**2 * (mom0.header['BMAJ']/2) * (mom0.header['BMIN']/2) / np.log(2)  #rad^2
    L_CO = I_CO * Omega_beam * (1e6*distance)**2
    Mmol_kcube = L_CO * aco
    #print(f'Omega_beam = {Omega_beam}')

    emom0 = fits.open(dir+'maps/'+source+'_12CO_emom0_'+mask+'_se'+session+'.fits')[0]
    I_CO_err = np.nansum(emom0.data[valid]) / Npixperbeam  #
    Mmol_err = I_CO_err * Omega_beam * (1e6*distance)**2 * aco
       
    print(f'I_tot = {I_CO} +/- {I_CO_err} K km/s')  
    print(f'L_tot = {L_CO} K km/s pc2')
    print(f'M_mol (using K cube) = {Mmol_kcube} +/- {Mmol_err} Msun')

    return I_CO, I_CO_err, Mmol_kcube, Mmol_err 


table_all = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')
source_list = ['NGC0001', 'NGC0169', 'NGC1056', 'NGC2540', 'NGC2596', 'UGC04245', 'UGC05396', 'UGC08322', 'UGC08781']
session_list = ['1+36+38+55', '2', '56', '7+48+59', 'all5', '8+32+35+37', '51+59', '41+54', '44+46+68']
method_list = ['block_noR25', 'block_noR25', 'block_noR25', 'block_noR25', 'block_noR25', 'block_noR25', 'Havexpand', 'block_noR25', 'block_noR25'] 

Ico21_list = np.full((len(source_list),2), np.nan)
Ico10_list = np.full((len(source_list),2), np.nan)
Ico_test_list = np.full((len(source_list),2), np.nan)

for n, (source, method, session) in enumerate(zip(source_list, method_list, session_list)):
    distance = float(table_all[table_all[:,0]==source][0, 24])
    z = 70 * distance / 299792.458
    print(f'{source}: D = {distance} Mpc, z = {z}')

    # Generate GBT cube with units of *Jy/beam* km/s
    GBT_natv = fits.open('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session_list[n]+'.fits')[0]
    freq = GBT_natv.header['RESTFRQ'] #Hz
    Bmaj = GBT_natv.header['BMAJ'] #deg
    Bmin = GBT_natv.header['BMIN'] #deg
    f_K2Jy = 1. / Jy2K(freq, Bmaj, Bmin)
    new_data = GBT_natv.data * f_K2Jy
    new_header = GBT_natv.header.copy()
    new_header['BUNIT'] = 'Jy beam-1'
    #fits.writeto('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session_list[n]+'_Jykms.fits', new_data, new_header, overwrite=True)  

    # Convolve GBT cubes that were run with pipeline masks to ACA 12" round beam 
    gbtcube = SpectralCube.read('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session_list[n]+'.fits')
    acacube = SpectralCube.read('data/ACA/'+source+'_7m_co21_smo12.pbcor_Kkms.fits')
    beam = radio_beam.Beam(major=acacube.header['BMAJ']*u.deg, minor=acacube.header['BMIN']*u.deg, pa=acacube.header['BPA']*u.deg)  
    cube_conv = gbtcube.convolve_to(beam)
    #cube_conv.write('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_segbtsmo12wR25.fits', overwrite=True) 

    # Generate convolved 12" GBT cube with units of *Jy/beam* km/s
    GBT_smo12wmask = fits.open('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_segbtsmo12wR25.fits')[0]
    freq = GBT_smo12wmask.header['RESTFRQ'] #Hz
    Bmaj = GBT_smo12wmask.header['BMAJ'] #deg
    Bmin = GBT_smo12wmask.header['BMIN'] #deg
    f_K2Jy = 1. / Jy2K(freq, Bmaj, Bmin)
    new_data = GBT_smo12wmask.data * f_K2Jy
    new_header = GBT_smo12wmask.header.copy()
    new_header['BUNIT'] = 'Jy beam-1'
    #fits.writeto('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_segbtsmo12wR25_Jykms.fits', new_data, new_header, overwrite=True) 

    # Generate ACA cube with units of Jy/beam *km/s*
    ACA_smo12 = SpectralCube.read('data/ACA/'+source+'_7m_co21_smo12.pbcor.fits').with_spectral_unit(u.km/u.s)
    #ACA_smo12.write('data/'+source+'_12CO_rebase5_smooth1.3_hanning2_se7msmo12_Jykms.fits', overwrite=True)

    print('- Convolved 12" GBT cube:')
    if source=='NGC0169':
        Ico_gbtsmo12, eIco_gbtsmo12, _, _ = Mmol(source, 'flat', 'gbtsmo12wR25', 4.35)
    else:
        Ico_gbtsmo12, eIco_gbtsmo12, _, _ = Mmol(source, 'Havexpand', 'gbtsmo12wR25', 4.35) #
    print('- Convolved 12" ACA cube:')
    Ico_7msmo12, eIco_7msmo12, _, _ = Mmol(source, method, '7msmo12', 6.7) #

    Ico10_list[n, 0] = Ico_gbtsmo12
    Ico21_list[n, 0] = Ico_7msmo12
    Ico10_list[n, 1] = eIco_gbtsmo12
    Ico21_list[n, 1] = eIco_7msmo12


co10_xerr_low = np.log10(Ico10_list[:,0]) - np.log10(Ico10_list[:,0] - Ico10_list[:,1])
co10_xerr_high = np.log10(Ico10_list[:,0] + Ico10_list[:,1]) - np.log10(Ico10_list[:,0])
co10_xerr = np.vstack((co10_xerr_low, co10_xerr_high))

co21_xerr_low = np.log10(Ico21_list[:,0]) - np.log10(Ico21_list[:,0] - Ico21_list[:,1])
co21_xerr_high = np.log10(Ico21_list[:,0] + Ico21_list[:,1]) - np.log10(Ico21_list[:,0])
co21_xerr = np.vstack((co21_xerr_low, co21_xerr_high))

plt.rc("axes", linewidth=1.5)
plt.rc("font", size=14)

plt.figure(figsize=(4.5,4.))

plt.errorbar(np.log10(Ico10_list[:,0]), np.log10(Ico21_list[:,0]), xerr=co10_xerr, yerr=co21_xerr, fmt='o', mfc='grey', mec='k', ecolor='grey', elinewidth=1.5)
plt.plot(np.arange(0.5,2.3,0.1), np.arange(0.5,2.3,0.1), 'k--', label=r'R$_{21}$ = 1')
plt.plot(np.arange(0.5,2.3,0.1), np.arange(0.5,2.3,0.1) + np.log10(0.5), 'k-.', label=r'R$_{21}$ = 0.5')
plt.plot(np.arange(0.5,2.3,0.1), np.arange(0.5,2.3,0.1) + np.log10(0.3), 'k:', label=r'R$_{21}$ = 0.3')

plt.text(1.1, 0.3, 'GBT CO(1-0): this work', fontsize=10)
plt.text(1.1, 0.1, 'ACA CO(2-1): Villanueva+24', fontsize=10)

plt.xlabel(r'$\log\ I_\mathrm{CO(1-0)}^\mathrm{GBT}$ (K km s$^{-1}$)', fontsize=16)  # w/mask
plt.ylabel(r'$\log\ I_\mathrm{CO(2-1)}^\mathrm{ACA}$ (K km s$^{-1}$)', fontsize=16)
plt.legend(fontsize=12)
plt.savefig('plots/flux_compare_errbar.pdf', bbox_inches='tight', pad_inches=0.02)

plt.show()