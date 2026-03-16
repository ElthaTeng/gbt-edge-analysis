import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

plotdir = './plots/'

table = np.genfromtxt('galaxy_parameters.csv', delimiter=',', skip_header=1, dtype='str')
galaxies = list(table[:,0])

fig = plt.figure(figsize=(12, 14.5))

gs = fig.add_gridspec(nrows = 9, ncols = 7, height_ratios=(1,1,1,1,1,1,1,1,1), width_ratios = (1,1,1,1,1,1,1),
                      left=0.01, right=0.99, bottom=0.01, top=0.99,
                      wspace=0.03, hspace=0.03)

for i, gal in enumerate(galaxies):      

    img = mpimg.imread('cutouts/'+gal+'.jpg')

    ax = fig.add_subplot(gs[i//7, i%7])

    im = ax.imshow(img)

    ax.tick_params(axis="x", labelsize=14, bottom=False, labelbottom=False)
    ax.tick_params(axis="y", labelsize=14, left=False, labelleft=False)

    labely = gal
    starneed = ['ngc0628', 'ngc4689', 'ngc3521', 'ngc4254', 'ngc5248', 'ngc4941', 'ngc4536', 'ngc4569']
    labely += r'$\star$' if gal in starneed else ''


    if gal in ['NGC0932', 'NGC3406NED01', 'NGC5216', 'NGC5631', 'UGC02222', 'UGC03960', 'UGC08234', 'UGC09629']:
        ax.text(0.05*img.shape[1], 0.15*img.shape[0], s= labely, color='r', fontweight = 'bold')
    elif gal in ['NGC3106', 'NGC3619', 'NGC5157', 'NGC6154', 'NGC6338', 'UGC04136', 'UGC08322', 'UGC10097', 'UGC10905', 'IC0674', 'IC3598']:
        ax.text(0.05*img.shape[1], 0.15*img.shape[0], s= labely, color='lime', fontweight = 'bold')
    else:
        ax.text(0.05*img.shape[1], 0.15*img.shape[0], s= labely, color='w', fontweight = 'bold')

#####################
### save the plot
savename = os.path.join(plotdir, 'sdss_gallery')

fig.savefig(savename + '_fix62_color.pdf', dpi = 300)
