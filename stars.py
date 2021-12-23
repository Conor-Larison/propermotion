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


datadir = '/home/conor/.data/000717/c01/q2/zg'

'''coords = '/home/conor/.data/000717/c01/q2/zg/0190828158530_0coords.txt'



table = ascii.read(coords)
ra1 = [0]*len(table['RA'])
dec1 = [0]*len(table['DEC'])
rad = [0]*len(table['RA'])
decd = [0]*len(table['DEC'])





for i in range(1,len(table['RA'])):
    ra1[i] = table['RA'][i]
    dec1[i] = table['DEC'][i]
    rad[i] = table['RAD'][i]
    decd[i] = table['DECD'][i]'''


zuds.init_db()
refpaths = os.listdir(datadir)


refimpaths = []
refmskpaths = []


for file in refpaths:
	if 'stamp.fits' in file:
		refimpaths.append('/home/conor/.data/000717/c01/q2/zg/' +file)
	elif 'mask.fits' in file:
		refmskpaths.append('/home/conor/.data/000717/c01/q2/zg/' +file)
    
refimpaths.sort()
refmskpaths.sort()

fname = '/home/conor/ZTF/zuds-pipeline/coords.txt'

table = ascii.read(fname)


rad = []
decd = []
for i in range(len(table['RA'])):
    rad.append(table['RAD'][i])
    decd.append(table['DECD'][i])

print(rad)

scis = []
catalogs = []
ra = np.empty(1)
dec = np.empty(1)
a = np.empty(1)
b = np.empty(1)
sn = np.empty(1)
days = []

for ipath, mpath in zip(refimpaths, refmskpaths):

    

    # create a python object for each file

    psf = fakes.measure_psf(ipath,persist = True)


    


    f = fits.open('/home/conor/ZTF/zuds-pipeline/coadd.fits')


    w = WCS(f[0].header)

    hdu = fits.open(ipath)[0]

    hdu.header.update(w.to_header())

    hdu.writeto(ipath, overwrite=True)




    sci = zuds.ScienceImage.from_file(ipath)
    msk = zuds.MaskImage.from_file(mpath)

    msk.parent_image = sci

    msk.data = msk.data.astype('int32')
    msk.save() 
    msk.load()
    sci.header['CRVL1'] = rad[0]
    sci.header['CRVL2'] = decd[0]

    sci.save()
    sci.load()
   

    cat = zuds.run_sextractor(sci, sextractor_kws = {'PSF_NAME': ipath[0:-4] + 'psf'})

    


    # associate each mask image with a science image
    
    zuds.calibrate_astrometry(sci,inplace = True)
    # add the objects to the database session
    #zuds.DBSession().add_all([sci, msk])
    catalog = zuds.PipelineFITSCatalog.from_image(sci)
    
    scis.append(sci)
    w = WCS(ipath)
    hdulist = fits.open(ipath)
    reg = zuds.PipelineRegionFile.from_catalog(catalog)
    

    
    