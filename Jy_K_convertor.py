import numpy as np

# input parameters:
# freq in Hz, bmaj and bmin in deg

def Jy2K(freq, bmaj, bmin):

    k = 1.380648e-23
    c = 299792458.
    lamb = c / freq

    beam_size = np.array((bmaj, bmin))  #deg 
    a = (beam_size[0]/2) * np.pi/180. 
    b = (beam_size[1]/2) * np.pi/180.
    omega = np.pi * a * b / np.log(2)

    f_Jy2K = lamb**2 / (2*k*omega * 10**26)

    print(f'Omega = {omega} (sr or radian^2)')
    print(f'Tb(Kelvin) to S(Jy/beam), factor = {1/f_Jy2K}')
    print(f'S(Jy/beam) to Tb(Kelvin), factor = {f_Jy2K}')

    return f_Jy2K
