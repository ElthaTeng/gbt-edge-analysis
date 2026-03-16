import numpy as np
from astropy.io import fits

# Combine NGC169 Pipe3D cubes: use NGC0169_2 as base and combine it with NGC5929_1
pipe3d_1 = fits.open('data/EDGE-CALIFA/new-py-v2.3/NGC0169_1.Pipe3D.cube.fits')
pipe3d_2 = fits.open('data/EDGE-CALIFA/new-py-v2.3/NGC0169_2.Pipe3D.cube.fits')

base_data = pipe3d_2[5].data.copy()
base_novalue = np.logical_or(base_data == 0, np.isnan(base_data))
base_data[base_novalue] = pipe3d_1[5].data[base_novalue]

base_data_ssp = pipe3d_2[1].data.copy()
base_novalue_ssp = np.logical_or(base_data_ssp == 0, np.isnan(base_data_ssp))
base_data_ssp[base_novalue_ssp] = pipe3d_1[1].data[base_novalue_ssp]

pipe3d_2[5].data = base_data
pipe3d_2[1].data = base_data_ssp
pipe3d_2.writeto('data/EDGE-CALIFA/new-py-v2.3/NGC0169.Pipe3D.cube.fits')

# Combine NGC5929 Pipe3D cubes: use NGC5929_0 as base and combine it with NGC5929_1
pipe3d_1 = fits.open('data/EDGE-CALIFA/new-py-v2.3/NGC5929_1.Pipe3D.cube.fits')
pipe3d_2 = fits.open('data/EDGE-CALIFA/new-py-v2.3/NGC5929_0.Pipe3D.cube.fits')

base_data = pipe3d_2[5].data.copy()
base_novalue = np.logical_or(base_data == 0, np.isnan(base_data))
base_data[base_novalue] = pipe3d_1[5].data[base_novalue]

base_data_ssp = pipe3d_2[1].data.copy()
base_novalue_ssp = np.logical_or(base_data_ssp == 0, np.isnan(base_data_ssp))
base_data_ssp[base_novalue_ssp] = pipe3d_1[1].data[base_novalue_ssp]

pipe3d_2[5].data = base_data
pipe3d_2[1].data = base_data_ssp
pipe3d_2.writeto('data/EDGE-CALIFA/new-py-v2.3/NGC5929.Pipe3D.cube.fits')