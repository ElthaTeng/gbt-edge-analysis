from mkGBTmaps import *
import time
import numpy as np

galaxy_list = np.genfromtxt('galaxy_list_customfix.csv', delimiter=',', skip_header=1, dtype='str')[-1,:].reshape(1,8)  #[0:9]  #[-4:, :]
methods = ['Havfield', 'rotnoHa+', 'rotnoHa-', 'flat', 'block', 'datacube']  #['datacube']  #
create_file = True
broad_datamask = True

for i in range(len(galaxy_list)):
    
    start_time = time.time()

    galaxy = galaxy_list[i, 0]
    session = galaxy_list[i, 6]
    
    for method_n, mask_method in enumerate(methods):    
        if mask_method=='datacube':
            version = None 
            get_mask(galaxy, mask_method, session, version, broad=broad_datamask, write_fits=create_file)

            cutoff = 0.05 #if broad_datamask else 0.05  #0.02 #0.01 #
            expand_mask(galaxy, session, cutoff=cutoff, write_fits=create_file)  
            mask_method = mask_method+'_expand'

        else:
            version = galaxy_list[i, method_n+1]
            get_mask(galaxy, mask_method, session, version, write_fits=create_file)
        
        apply_mask(galaxy, mask_method, session, write_fits=create_file)
        get_maps(galaxy, mask_method, session, write_fits=create_file)

    
    #compare_mom0(galaxy, session, save_fig = create_file, interactive=True)

    print(f'{galaxy} done; took {np.round(time.time() - start_time, 1)} sec.')