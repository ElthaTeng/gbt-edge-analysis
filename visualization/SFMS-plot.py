import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats

catalog = np.genfromtxt('edge_califa.csv', delimiter=',', skip_header=61, dtype='str')  #np.genfromtxt('GBTEDGE.cat', skip_header=2, dtype='str')

# GBT Sample
table_gbt = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')
galaxies_gbt = list(table_gbt[:,0])
logSFR_list_gbt = []
logMst_list_gbt = []

for nameG in galaxies_gbt:
    logSFR = catalog[catalog[:,1]==nameG][0,26] 
    logMst = catalog[catalog[:,1]==nameG][0,24] 

    logSFR_list_gbt.append(logSFR)
    logMst_list_gbt.append(logMst)

# CARMA Sample
table_carma = np.genfromtxt('carma-table1.csv', delimiter=',', skip_header=2, dtype='str')
galaxies_carma = list(table_carma[:,0])
logSFR_list_carma = []
logMst_list_carma = []

for nameG in galaxies_carma:
    galname = nameG.replace('"', '')
    logSFR = catalog[catalog[:,1]==galname][0,26] 
    logMst = catalog[catalog[:,1]==galname][0,24] 

    logSFR_list_carma.append(logSFR)
    logMst_list_carma.append(logMst)

# Write table to galaxy_parameters_SFMS.csv
dSFMS = np.array(logSFR_list_gbt, dtype='float') - (0.81 * np.array(logMst_list_gbt, dtype='float') - 8.34)
MS = dSFMS > -0.5
T = Table([galaxies_gbt, logSFR_list_gbt, logMst_list_gbt, list(dSFMS), list(MS)], names = ('Galaxy', 'logSFR', 'logMst', 'dSFMS', 'MS?'))

plt.rc("axes", linewidth=1.5)
plt.rc("font", size=14)
plt.figure(figsize=(6.4,5.8))

# SFMS plot
xmin=8.2; xmax=11.8; ymin=-4; ymax=1.5
X, Y = np.mgrid[xmin:xmax:30j, ymin:ymax:30j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([np.array(catalog[:,24], dtype='float'), np.array(catalog[:,26], dtype='float')])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)
plt.contour(X, Y, Z, colors='silver', linestyles='-') #, levels=6
print(Z.min(),Z.max())

plt.scatter(np.array(logMst_list_gbt, dtype='float'), np.array(logSFR_list_gbt, dtype='float'), c='k', marker='o', label='GBT-EDGE', zorder=5)
plt.scatter(np.array(logMst_list_carma, dtype='float'), np.array(logSFR_list_carma, dtype='float'), c='grey', marker='*', label='CARMA-EDGE')

xrange = np.arange(8, 12, 0.05)
f_MS = 0.81 * xrange - 8.34
plt.plot(xrange, f_MS, c='darkblue', ls='--')  #, label='SFMS (Cano-Diaz+16)'
plt.plot(xrange, 0.86 * xrange - 10.32, c='darkred', ls='--')  #, label='Red galaxies (Cano-Diaz+16)'

plt.fill_between(xrange, f_MS - 0.5, f_MS + 0.7, alpha=0.2, color='tab:blue', label='Main Sequence (43)')
plt.fill_between(xrange, f_MS - 1, f_MS - 0.5, alpha=0.2, color='tab:green', label='Green Valleys (11)')
plt.fill_between(xrange, f_MS - 2.5, f_MS - 1, alpha=0.2, color='tab:red', label='Red Galaxies (8)')

plt.legend(fontsize=12, loc='lower left')
plt.xlabel(r'$\log\ M_\mathrm{star}$ (M$_\odot$)', fontsize=16) 
plt.ylabel(r'$\log$ SFR (M$_\odot$ yr$^{-1}$)', fontsize=16)
plt.ylim(-4, 1.5)
plt.xlim(8., 11.8)
plt.savefig('plots/SFMS-plot_full.pdf', bbox_inches='tight', pad_inches=0.02)
plt.show()



## Histogram of depletion times
table_gbt_float = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1)
tdep_blockpref = table_gbt_float[:, 8] 
tdep_datapref = table_gbt_float[:, 9]
tdep_datapref_B13 = table_gbt_float[:, 27]
tdep_datapref_T24 = table_gbt_float[:, 31]
tdep_datapref_B13_Zgrad = table_gbt_float[:, 35]
tdep_datapref_SL24_Zgrad = table_gbt_float[:, 39]

error_add = table_gbt_float[:, 27]

etdep_blockpref = table_gbt_float[:, 15]
etdep_datapref = table_gbt_float[:, 16]
etdep_datapref_B13 = table_gbt_float[:, 28]
etdep_datapref_T24 = table_gbt_float[:, 32]
etdep_datapref_B13_Zgrad = table_gbt_float[:, 36]
etdep_datapref_SL24_Zgrad = table_gbt_float[:, 40]

Mmol_blockpref = table_gbt_float[:,17] 
eMmol_blockpref = table_gbt_float[:,18] 
Mmol_datapref = table_gbt_float[:,19] 
eMmol_datapref = table_gbt_float[:,20] 
Mmol_datapref_T24 = table_gbt_float[:,29] 
eMmol_datapref_T24 = table_gbt_float[:,30]  
Mmol_datapref_B13_Zgrad = table_gbt_float[:,33] 
eMmol_datapref_B13_Zgrad = table_gbt_float[:,34]
Mmol_datapref_SL24_Zgrad = table_gbt_float[:,37] 
eMmol_datapref_SL24_Zgrad = table_gbt_float[:,38]  

SFR_Ha = table_gbt_float[:,21] 
SFR_33M = table_gbt_float[:,22]

idx_GV = table_gbt[:,14] == 'GV'
idx_RG = table_gbt[:,14] == 'RG'
idx_MS = table_gbt[:,14] == 'MS'
idx_GVRG = table_gbt[:,13] == 'FALSE'

### SET THESE BASED ON ALPHA_CO CHOICE !! 
tdep = tdep_datapref#_B13_Zgrad
Mmol = Mmol_datapref#_B13_Zgrad
etdep_raw = etdep_datapref#_B13_Zgrad

etdep = np.sqrt((tdep * error_add * 0.01)**2 + etdep_raw**2) 

counts, bins = np.histogram(tdep, bins=10, range=(0.03, 2000))
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))

plt.figure(figsize=(6.4, 2.4))  
plt.hist(tdep[idx_MS], bins=logbins, histtype='step', color='C0', hatch='//', label='Main Sequence (43)')
plt.hist(tdep[idx_GV], bins=logbins, histtype='stepfilled', color='C2', alpha=0.5, label='Green Valleys (11)')
plt.hist(tdep[idx_RG], bins=logbins, histtype='stepfilled', color='C3', alpha=0.5, label='Red Galaxies (8)')
plt.hist(tdep, bins=logbins, histtype='step', color='k', label='Total (62)')
plt.xscale('log')

plt.axvline(x=np.median(tdep[idx_MS]), c='blue', ls='--')
plt.axvline(x=np.median(tdep[idx_GV]), c='green', ls='--')
plt.axvline(x=np.median(tdep[idx_RG]), c='red', ls='--')
print(f'MS: tdep = {np.median(tdep[idx_MS])} + {np.percentile(tdep[idx_MS], 84) - np.median(tdep[idx_MS])} - {np.median(tdep[idx_MS]) - np.percentile(tdep[idx_MS], 16)}')
print(f'GV: tdep = {np.median(tdep[idx_GV])} + {np.percentile(tdep[idx_GV], 84) - np.median(tdep[idx_GV])} - {np.median(tdep[idx_GV]) - np.percentile(tdep[idx_GV], 16)}')
print(f'RG: tdep = {np.median(tdep[idx_RG])} + {np.percentile(tdep[idx_RG], 84) - np.median(tdep[idx_RG])} - {np.median(tdep[idx_RG]) - np.percentile(tdep[idx_RG], 16)}')

plt.annotate(r'MW $\alpha_\mathrm{CO}$', weight='bold', fontsize=14, xy=(0.05, 0.8), xycoords='axes fraction', color='k')
plt.xlabel(r'$t_\mathrm{dep}$ (Gyr)', fontsize=16)  
plt.ylabel('N of galaxies', fontsize=16)
plt.legend(fontsize=12, loc='upper right')
plt.savefig('plots/hist_tdep_datapref_MW.pdf', bbox_inches='tight', pad_inches=0.02)
plt.show()



# dSFMS-dtdep plot 
logfmol = np.log10(Mmol) - np.array(logMst_list_gbt, dtype='float')
valid = tdep > etdep

# Uncomment below ONLY if running "nolab" for B13/SL24
#plt.figure(figsize=(5.4, 4.8))  

## xCOLDGASS data
logsfr_S17 = fits.open('xCOLDGASS_PubCat.fits')[1].data['LOGSFR_BEST']
logMst_S17 = fits.open('xCOLDGASS_PubCat.fits')[1].data['LOGMSTAR']
logMmol_S17 = fits.open('xCOLDGASS_PubCat.fits')[1].data['LOGMH2']
logMmol_lim_S17 = fits.open('xCOLDGASS_PubCat.fits')[1].data['LIM_LOGMH2']
aco_S17 = fits.open('xCOLDGASS_PubCat.fits')[1].data['XCO_A17']

# compute tdep using MW alpha_CO (reversing the Z-based alpha_CO used in S17)
tdep_S17 = 10**logMmol_S17 / aco_S17 * 4.35 / 10**logsfr_S17 / 1e9  
tdep_lim_S17 = 10**logMmol_lim_S17 / aco_S17 * 4.35 / 10**logsfr_S17 / 1e9  #
dSFMS_S17 = logsfr_S17 - (0.81 * logMst_S17 - 8.34)
valid_S17 = dSFMS_S17 > -3

plt.scatter(tdep_S17[valid_S17], dSFMS_S17[valid_S17], marker='.', c='lightgray', label='xCOLD GASS')

## iEDGE APEX data 
table_iedge = Table.read('iedge_v1.ecsv')
sfr_iedge = table_iedge['Glob_SFR']
Mst_iedge = table_iedge['Glob_Mstar']
Mmol_iedge = table_iedge['APEX_Glob_Mmol']
SNR_iedge = table_iedge['APEX_Glob_SNR']
tdep_iedge = Mmol_iedge / sfr_iedge / 1e9
dSFMS_iedge = np.log10(sfr_iedge) - (0.81 * np.log10(Mst_iedge) - 8.34)

# apply a S/N > 5 cut on the iEDGE-APEX sample
valid_iedge = SNR_iedge > 5  

plt.scatter(tdep_iedge[valid_iedge], dSFMS_iedge[valid_iedge], marker='+', c='lightgray', label='iEDGE-APEX')

# GBT-EDGE data
plt.errorbar(tdep[valid], dSFMS[valid], xerr=etdep[valid], fmt='o', mfc='grey', mec='k', ecolor='grey', elinewidth=1.5)
plt.scatter(tdep[valid], dSFMS[valid], marker='o', c=logfmol[valid], edgecolors='k', cmap='Spectral_r', zorder=10, vmin=-2.5, vmax=0)
plt.scatter(tdep[~valid], dSFMS[~valid], marker='<', c=logfmol[~valid], edgecolors='k', cmap='Spectral_r', zorder=8, vmin=-2.5, vmax=0)

# Showing colorbar for Mmol/Mstar (disable it ONLY for B13/SL24)
cb = plt.colorbar()
cb.set_label(r'$\log(M_\mathrm{mol}/M_\mathrm{star})$', rotation=270, labelpad=20)

plt.scatter(np.nan, np.nan, c='dimgrey', marker='<', label=r'CO upper limits')

plt.axhspan(- 0.5, 0.7, alpha=0.2, color='tab:blue', label='Main Sequence')
plt.axhspan(- 1, - 0.5, alpha=0.2, color='tab:green', label='Green Valleys')
plt.axhspan(- 2.5, - 1, alpha=0.2, color='tab:red', label='Red Galaxies')

plt.legend(fontsize=10, loc='lower left')
plt.xlabel(r'$t_\mathrm{dep}$ (Gyr)', fontsize=16) 
plt.ylabel(r'$\Delta$ SFMS (dex)', fontsize=16)
plt.xscale('log')
plt.xlim(0.005, 2000)  
plt.ylim(-2.7, 1.02)
plt.savefig('plots/dSFMS-tdep_datapref_errbar_fmol_S17MW_iedgesn5.pdf', bbox_inches='tight', pad_inches=0.02)  #_sfr33m
plt.show()

