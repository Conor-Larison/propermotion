import glob
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
import pdb
import numpy as np
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
import sys

file = '/home/conor/.data/000717/c01/q2/zg/ztf_20190828158530_000717_zg_c01_o_q2_sciimg.fits'


a = fits.open(file)
image_file = get_pkg_data_filename(file)
fits.info(image_file)

image_data = fits.getdata(image_file, ext=0)

table = ascii.read('coords.txt')

x = [0]*len(table['X'])
y = [0]*len(table['Y'])

for i in range(len(table['X'])):
	x[i] = table['X'][i]
	y[i] = table['Y'][i]
plt.imshow(image_data)
plt.colorbar()
plt.scatter(x,y, s = .05)

plt.savefig('/home/conor/ZTF/fig1.png')

plt.show()






