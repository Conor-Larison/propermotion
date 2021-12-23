import zuds
import glob
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
import pdb
import numpy as np
import datetime as dt
import fakes



datadir = '/home/conor/.data/stackdata1/'



zuds.init_db()
refpaths = os.listdir(datadir)


refimpaths = []
refmskpaths = []


for file in refpaths:
	if 'sciimg.fits' in file:
		refimpaths.append('/home/conor/.data/stackdata1/' +file)
	elif 'mskimg.fits' in file:
		refmskpaths.append('/home/conor/.data/stackdata1/' +file)
    
refimpaths.sort()
refmskpaths.sort()



scis = []



for ipath,mpath in zip(refimpaths,refmskpaths):

	sci = zuds.ScienceImage.from_file(ipath)

	
	msk = zuds.MaskImage.from_file(mpath)
	msk.parent_image = sci



	scis.append(sci)

	zuds.DBSession().add_all([sci, msk])



index = 0
for i in range(len(scis)):

	chunk = []
	dates = []
	if i == len(scis) - 1:
		for l in range(index,i):
			chunk.append(scis[l])
			dates.append(scis[l].header['SHUTOPEN'][0:4])
	elif i!=0 and scis[i].header['SHUTOPEN'][0:4] == scis[i-1].header['SHUTOPEN'][0:4]:
		continue
	elif i == 0:
		continue
	else:
		for l in range(index,i):
			chunk.append(scis[l])
			dates.append(scis[l].header['SHUTOPEN'][0:4])

		index = i


	sci1 = zuds.ScienceCoadd.from_images(chunk, outname = '/home/conor/.data/stackdata1/' + str(dates[0]) + 'coadd.fits')

	psf = fakes.measure_psf('/home/conor/.data/stackdata1/' + str(dates[0]) + 'coadd.fits',persist = True)

	cat = zuds.PipelineFITSCatalog.from_image(sci1)

	zuds.calibrate_astrometry(sci1,inplace = True)

	catalog = zuds.PipelineFITSCatalog.from_image(sci1)