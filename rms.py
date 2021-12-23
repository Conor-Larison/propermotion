import glob
import os
import csv
import pandas as pd
import matplotlib.pyplot as plt
from itertools import zip_longest
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
import pdb
import numpy as np
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
import sys
from scipy import stats
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u






file_list = os.listdir('/home/conor/.data/000717/c01/q2/zg')


files = []
coords = []

for file in file_list:
	if 'stamp.cat' in file:
		coords.append('/home/conor/.data/000717/c01/q2/zg/'+file)
	elif 'stamp.fits' in file:
		files.append('/home/conor/.data/000717/c01/q2/zg/' + file)
coords.sort()
files.sort()


i = 0
for image in files:
	f = fits.open(image)[0]
	if f.header['MAGLIM'] >= 20.5 and f.header['SEEING'] <= 2:
		base_image = image
		base_coord = coords[i]
		break
		i+=1

f = fits.open(base_coord)

ra_base = []
dec_base = []
mag_base = []

for i in range(len(f[2].data['XWIN_WORLD'])):
	ra_base.append(f[2].data['XWIN_WORLD'][i])
	dec_base.append(f[2].data['YWIN_WORLD'][i])
	mag_base.append(f[2].data['MAG_AUTO'])


files.remove(base_image)
coords.remove(base_coord)


ra = []
dec = []
mag = []


count = 0
for reg in coords:
	
	f = fits.open(reg)

	for i in range(len(f[2].data['XWIN_WORLD'])):
		ra.append(f[2].data['XWIN_WORLD'][i])
		dec.append(f[2].data['YWIN_WORLD'][i])
		mag.append(f[2].data['MAG_AUTO'][i])
	

x = []

stars = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=ra_base*u.degree,dec=dec_base*u.degree),0.000277778*u.degree)
	

ra_matched = []
dec_matched = []
ra_bmatched = []
mag_matched = []
for j in range(len(sep2d.arcsecond)):
	x.append(sep2d.arcsecond[j])
	ra_matched.append(ra[idxself[j]])
	dec_matched.append(dec[idxself[j]])
	ra_bmatched.append(ra_base[idxs[j]])
	mag_matched.append(mag[idxself[j]])


		
stds = []
index = 0

newer_mags = []

for i in range(len(x)):
	x_chunk = []
	ra_chunk = []
	dec_chunk = []
	new_mags = []
	if i!=0 and np.abs(ra_matched[i] - ra_matched[i-1]) < 0.000277778:
		continue
	elif i == 0:
		continue
	else:
		for l in range(index,i):
			x_chunk.append(x[l])
			ra_chunk.append(ra_matched[l])
			dec_chunk.append(dec_matched[l])
			new_mags.append(mag_matched[l])

		
		if len(x_chunk) != 1:
			stds.append(np.median(x_chunk))
			newer_mags.append(np.median(new_mags))


		index = i




fig, ax = plt.subplots()
ax.add_artist(plt.scatter(newer_mags,stds))
#ax.add_artist(plt.scatter(mag2,x2,color = 'r'))
plt.xlabel("Magnitude")
plt.ylabel("Combined RMS (arsec)")
plt.yscale('log')
plt.title('Scatter relative to ZTF vs relative to Gaia')
plt.show()
