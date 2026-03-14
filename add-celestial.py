import numpy as np
from astropy.io import fits

source_list = list(np.genfromtxt('galaxy_list.csv', delimiter=',', skip_header=1, dtype='str')[:, 0])

for i, source in enumerate(source_list):

    pipe3d = fits.open('data/'+source+'.Pipe3D.cube.fits') 

    pipe3d[5].header['CTYPE1'] = pipe3d[0].header['CTYPE1']
    pipe3d[5].header['CRVAL1'] = pipe3d[0].header['CRVAL1']
    pipe3d[5].header['CRPIX1'] = pipe3d[0].header['CRPIX1']
    pipe3d[5].header['CUNIT1'] = pipe3d[0].header['CUNIT1']

    pipe3d[5].header['CTYPE2'] = pipe3d[0].header['CTYPE2']
    pipe3d[5].header['CRVAL2'] = pipe3d[0].header['CRVAL2']
    pipe3d[5].header['CRPIX2'] = pipe3d[0].header['CRPIX2']
    pipe3d[5].header['CUNIT2'] = pipe3d[0].header['CUNIT2']

    pipe3d.writeto(source+'.Pipe3D.cube.fits') 
    print(f'{source} done.')