import numpy as np
import matplotlib.pyplot as plt

version = 'datapref_B13_Zgrad' #

table_gbt = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1)
SFR_all = table_gbt[:,21]  # Msun/yr

if version == 'datapref':
    Mmol_all = table_gbt[:,19]
    eMmol_all = table_gbt[:,20]
    tdep_all = table_gbt[:,9]
elif version == 'blockpref':
    Mmol_all = table_gbt[:,17]
    eMmol_all = table_gbt[:,18]
    tdep_all = table_gbt[:,8]
elif version == 'datapref_T24':
    Mmol_all = table_gbt[:,29]
    eMmol_all = table_gbt[:,30]
    tdep_all = table_gbt[:,31]
elif version == 'datapref_B13_Zgrad':
    Mmol_all = table_gbt[:,33]
    eMmol_all = table_gbt[:,34]
    tdep_all = table_gbt[:,35]
elif version == 'datapref_SL24_Zgrad':
    Mmol_all = table_gbt[:,37]
    eMmol_all = table_gbt[:,38]
    tdep_all = table_gbt[:,39]

# Add uncertainty from masking and baseline variation
error_add = table_gbt[:, 27]
eMmol_all = np.sqrt((Mmol_all * error_add * 0.01)**2 + eMmol_all**2)  # eMmol_all + Mmol_all * error_add * 0.01  # 

### MAKE SURE TO USE THE CORRECT TABLE
idx_RG = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')[:,14] == 'RG'
idx_GV = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')[:,14] == 'GV'
idx_MS = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')[:,14] == 'MS'

logMmol = np.log10(Mmol_all)
logSFR = np.log10(SFR_all)     
xerr_low = logMmol - np.log10(Mmol_all - eMmol_all)
xerr_high = np.log10(Mmol_all + eMmol_all) - logMmol
valid = ~np.isnan(xerr_low)

galaxies = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')[:,0]
print(galaxies[~valid])

# plotting upper limits separately (preferred):
xerr = np.vstack((xerr_low[valid], xerr_high[valid]))

plt.rc("axes", linewidth=1.5)
plt.rc("font", size=14)
plt.figure(figsize=(6.4,5.8))

xerr_MS = np.vstack((xerr_low[valid * idx_MS], xerr_high[valid * idx_MS]))
plt.errorbar(logMmol[valid * idx_MS], logSFR[valid * idx_MS], xerr=xerr_MS, fmt='o', mfc='tab:blue', mec='k', ecolor='grey', elinewidth=1.5, label='Main Sequence')
plt.scatter(np.log10(Mmol_all)[~valid * idx_MS], logSFR[~valid * idx_MS], c='tab:blue', marker='<') # + eMmol_all

xerr_GV = np.vstack((xerr_low[valid * idx_GV], xerr_high[valid * idx_GV]))
plt.errorbar(logMmol[valid * idx_GV], logSFR[valid * idx_GV], xerr=xerr_GV, fmt='o', mfc='tab:green', mec='k', ecolor='grey', elinewidth=1.5, label='Green Valleys')
plt.scatter(np.log10(Mmol_all)[~valid * idx_GV], logSFR[~valid * idx_GV], c='darkgreen', marker='<') # + eMmol_all

xerr_RG = np.vstack((xerr_low[valid * idx_RG], xerr_high[valid * idx_RG]))
plt.errorbar(logMmol[valid * idx_RG], logSFR[valid * idx_RG], xerr=xerr_RG, fmt='o', mfc='tab:red', mec='k', ecolor='grey', elinewidth=1.5, label='Red Galaxies')
plt.scatter(np.log10(Mmol_all)[~valid * idx_RG], logSFR[~valid * idx_RG], c='darkred', marker='<') # + eMmol_all

plt.scatter(np.nan, np.nan, c='grey', marker='<', label=r'CO upper limits')

xrange = np.arange(5,12.5,0.05)
plt.plot(xrange, xrange-8, 'k:')
plt.plot(xrange, xrange-9, 'k:')
plt.plot(xrange, xrange-10, 'k:')

plt.text(6.1, -1.7, '0.1 Gyr', fontsize=11, color='k', rotation=40)
plt.text(6.2, -2.6, '1 Gyr', fontsize=11, color='k', rotation=40)
plt.text(6.3, -3.5, '10 Gyr', fontsize=11, color='k', rotation=40)

plt.xlabel(r'$\log\ M_\mathrm{mol}$ (M$_\odot$)', fontsize=16) 
plt.ylabel(r'$\log$ SFR (M$_\odot$ yr$^{-1}$)', fontsize=16)
plt.ylim(-4, 1.5)  
plt.xlim(6., 10.7) 
plt.legend(fontsize=13, loc='upper left') 
plt.savefig('plots/KS-plot_'+version+'_grouped_fix62_err1k1b_m1mean.pdf', bbox_inches='tight', pad_inches=0.02)
plt.show()
