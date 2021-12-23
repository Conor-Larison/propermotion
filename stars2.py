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
from astropy.utils.data import get_pkg_data_filename


datadir = '/home/conor/.data/000717/c01/q2/zg'


zuds.init_db()
refpaths = os.listdir(datadir)


refimpaths = []
refmskpaths = []

for file in refpaths:
	if 'sciimg.fits' in file:
		refimpaths.append('/home/conor/.data/000717/c01/q2/zg/' +
file)
	elif 'mskimg.fits' in file:
		refmskpaths.append('/home/conor/.data/000717/c01/q2/zg/' +
file)


refimpaths.sort()
refmskpaths.sort()


scis = []
catalogs = []
ra = np.empty(1)
dec = np.empty(1)
a = np.empty(1)
b = np.empty(1)
sn = np.empty(1)
days = []
i = 0
ipath = '/home/conor/.data/000717/c01/q2/zg/ztf_20190226359919_000717_zg_c01_o_q2_sciimg.fits'
mpath = '/home/conor/.data/000717/c01/q2/zg/ztf_20190226359919_000717_zg_c01_o_q2_mskimg.fits'

# create a python object for each file




sci = zuds.ScienceImage.from_file(ipath)
msk = zuds.MaskImage.from_file(mpath)

'''msk.data = msk.data.astype('int32')
msk.save() 
msk.load()'''

# associate each mask image with a science image
msk.parent_image = sci
cat = zuds.PipelineFITSCatalog.from_image(sci)
zuds.ScienceCoadd.from_images([sci],outname = 'coadd.fits', sci_swarp_kws={'CENTER_TYPE': 'ALL','PIXELSCALE_TYPE': 'MEDIAN','PIXEL_SCALE': 1.012,'IMAGE_SIZE': '768,770','RESAMPLE': 'N'},mask_swarp_kws={'CENTER_TYPE': 'ALL','PIXELSCALE_TYPE': 'MEDIAN','PIXEL_SCALE': 1.012,'IMAGE_SIZE': '768,770','RESAMPLE': 'N'})

'''zuds.calibrate_astrometry(sci,inplace = True)
# add the objects to the database session
#zuds.DBSession().add_all([sci, msk])
catalog = zuds.PipelineFITSCatalog.from_image(sci)
catalogs.append(catalog)
scis.append(sci)
w = WCS(ipath)
hdulist = fits.open(ipath)
reg = zuds.PipelineRegionFile.from_catalog(catalog)
    
ra = []
dec = []
mag = []
for i in range(len(catalog.data['X_WORLD'])):
    if catalog.data['MAG_AUTO'][i] + sci.header['MAGZP'] < 19:
        ra.append(catalog.data['X_WORLD'][i])
        dec.append(catalog.data['Y_WORLD'][i])
        mag.append(catalog.data['MAG_AUTO'][i]+sci.header['MAGZP'])

output = [ra,dec,mag]

ascii.write(output,'/home/conor/.data/000717/c01/q2/zg/' + 'hello' + 'coords.txt', delimiter = '\t', names=['RA','DEC','MAG'],overwrite = True)'''