import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.convolution import Gaussian1DKernel

## SET THE SMOOTHING FACTOR (originally 1, which corresponds to ~15 km/s channel width)
smooth_factor = 6


# Load source and session info 
table_gbt = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')
source_list = list(table_gbt[:, 0]) 
session_list = list(table_gbt[:, 5])

# Generate spectrally smoothed cubes for all galaxies
for n, (source, session) in enumerate(zip(source_list, session_list)):
    cube = SpectralCube.read('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'.fits')
    print(cube.spectral_axis[0], cube.spectral_axis[-1], cube.header['CDELT3'])
    
    new_axis = np.arange(cube.spectral_axis[0].value, cube.spectral_axis[-1].value, smooth_factor*cube.header['CDELT3']) * u.km/u.s
    fwhm_factor = np.sqrt(8*np.log(2))
    smcube = cube.spectral_smooth(Gaussian1DKernel(smooth_factor/fwhm_factor))
    new_cube = smcube.spectral_interpolate(new_axis, suppress_smooth_warning=True)
    
    new_cube.write('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'_specsm6.fits', overwrite=True)


# Generate Tpeak maps for all galaxies
for n, (source, session) in enumerate(zip(source_list, session_list)):  
    cube = fits.open('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'_specsm6.fits')[0].data 
    cube_header = fits.open('data/'+source+'_12CO_rebase7_smooth1.3_hanning2_se'+session+'_specsm6.fits')[0].header.copy() 
    map_header = fits.open('maps/'+source+'_12CO_mom0_block_se'+session+'.fits')[0].header
    map_header['BUNIT'] = 'K'
    
    tpeak = np.max(cube, axis=0)

    fits.writeto('maps/pplsquare_erik_v2/'+source+'_12CO_Tpeak_se'+session+'_specsm6.fits', tpeak, map_header, overwrite=True)  
