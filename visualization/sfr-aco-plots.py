import numpy as np
import matplotlib.pyplot as plt

table_gbt = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')

idx_MS = table_gbt[:,14] == 'MS'
idx_GV = table_gbt[:,14] == 'GV'
idx_RG = table_gbt[:,14] == 'RG'

SFR_Ha = np.array(table_gbt[:,21], dtype='float')
SFR_33Myr = np.array(table_gbt[:,22], dtype='float')
SFR_12Myr = np.array(table_gbt[:,23], dtype='float')

plt.rc("axes", linewidth=1.5)
plt.rc("font", size=14)

# Compare SFR estimates using SSP vs. Ha 
SFR_SSP = SFR_33Myr

plt.figure(figsize=(4.2,4.))
plt.scatter(np.log10(SFR_Ha)[idx_MS], np.log10(SFR_SSP)[idx_MS], marker='o', c='C0', edgecolors='k', label='Main Sequence')
plt.scatter(np.log10(SFR_Ha)[idx_GV], np.log10(SFR_SSP)[idx_GV], marker='s', c='C2', edgecolors='k', label='Green Valleys')
plt.scatter(np.log10(SFR_Ha)[idx_RG], np.log10(SFR_SSP)[idx_RG], marker='D', c='C3', edgecolors='k', label='Red Galaxies')

plt.plot(np.arange(-2.5,1.5,0.1), np.arange(-2.5,1.5,0.1), 'k--')
plt.errorbar(1.3, -2.4, yerr=0.3, xerr=0.1, marker='.', c='k', elinewidth=1.5)

plt.grid(linestyle=':', linewidth=1)
plt.xlabel(r'$\log$ SFR (H$\alpha$)', fontsize=16)
plt.ylabel(r'$\log$ SFR (SSP$_\mathrm{< 33\,Myr}$)', fontsize=16)
plt.legend(fontsize=11)
plt.savefig('plots/SFR_compare_33Myr_fix62_arrow_Avcorr_final.pdf', bbox_inches='tight', pad_inches=0.02)
plt.show()


# Compare aco estimates using different prescriptions
Mmol_const = np.array(table_gbt[:,19], dtype='float')
Mmol_B13 = np.array(table_gbt[:,33], dtype='float')
Mmol_SL24 = np.array(table_gbt[:,37], dtype='float')
Mmol_T24 = np.array(table_gbt[:,29], dtype='float')
Lco = np.array(table_gbt[:,41], dtype='float')

aco_const = Mmol_const / Lco
aco_B13 = Mmol_B13 / Lco
aco_SL24 = Mmol_SL24 / Lco
aco_T24 = Mmol_T24 / Lco

plt.figure(figsize=(4.2,4.))
plt.scatter(aco_SL24[idx_MS], aco_B13[idx_MS], marker='o', c='C0', edgecolors='k', label='Main Sequence')
plt.scatter(aco_SL24[idx_GV], aco_B13[idx_GV], marker='s', c='C2', edgecolors='k', label='Green Valleys')
plt.scatter(aco_SL24[idx_RG], aco_B13[idx_RG], marker='D', c='C3', edgecolors='k', label='Red Galaxies')

plt.plot(np.arange(2.5,6.,0.1), np.arange(2.5,6.,0.1), 'k--')
plt.axhline(y=4.35, c='k', ls=':')
plt.axvline(x=4.35, c='k', ls=':')
plt.annotate(r'MW $\alpha_\mathrm{CO}$ = 4.35', weight=None, fontsize=11, xy=(0.5, 4.5), xycoords='data', color='k')

plt.grid(linestyle=':', linewidth=1)
plt.xlim(0,10)

plt.legend(fontsize=12, loc='lower right')
plt.xlabel(r'$\alpha_\mathrm{CO}$ (SL24)', fontsize=16)
plt.ylabel(r'$\alpha_\mathrm{CO}$ (B13)', fontsize=16)
#plt.savefig('plots/aco_SL24_vs_B13_Curti_fix62_final.pdf', bbox_inches='tight', pad_inches=0.02)
plt.show()


plt.figure(figsize=(4.2,4.))
plt.scatter(aco_T24[idx_MS], aco_B13[idx_MS], marker='o', c='C0', edgecolors='k', label='Main Sequence')
plt.scatter(aco_T24[idx_GV], aco_B13[idx_GV], marker='s', c='C2', edgecolors='k', label='Green Valleys')
plt.scatter(aco_T24[idx_RG], aco_B13[idx_RG], marker='D', c='C3', edgecolors='k', label='Red Galaxies')

plt.plot(np.arange(2.5,6.,0.1), np.arange(2.5,6.,0.1), 'k--')
plt.axhline(y=4.35, c='k', ls=':')
plt.axvline(x=4.35, c='k', ls=':')
plt.annotate(r'MW $\alpha_\mathrm{CO}$ = 4.35', weight=None, fontsize=11, xy=(6.2, 4.5), xycoords='data', color='k')

plt.grid(linestyle=':', linewidth=1)
plt.xlim(0, 10)
plt.xlabel(r'$\alpha_\mathrm{CO}$ (T24*)', fontsize=16)
#plt.ylabel(r'$\alpha_\mathrm{CO}$ (SL24)', fontsize=16)

#plt.savefig('plots/aco_T24_vs_B13_Curti_fix62_final.pdf', bbox_inches='tight', pad_inches=0.02)
plt.show()