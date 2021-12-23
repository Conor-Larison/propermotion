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
import tractor as tr


datadir = '/home/conor/.data/newdata/'




zuds.init_db()
refpaths = os.listdir(datadir)


refimpaths = []
refmskpaths = []


for file in refpaths:
    if 'sciimg.fits' in file:
        refimpaths.append('/home/conor/.data/newdata/' +file)
    elif 'mskimg.fits' in file:
        refmskpaths.append('/home/conor/.data/newdata/' +file)
    
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



for ipath, mpath in zip(refimpaths, refmskpaths):

    

    # create a python object for each file

    psf = fakes.measure_psf(ipath,persist = True)


    sci = zuds.ScienceImage.from_file(ipath)
    msk = zuds.MaskImage.from_file(mpath)

    msk.parent_image = sci

    cat = zuds.run_sextractor(sci, sextractor_kws = {'PSF_NAME': ipath[0:-4] + 'psf'})




    # associate each mask image with a science image
    
    zuds.calibrate_astrometry(sci,inplace = True)
    # add the objects to the database session
    #zuds.DBSession().add_all([sci, msk])
    catalog = zuds.PipelineFITSCatalog.from_image(sci)
    
    reg = zuds.PipelineRegionFile.from_catalog(catalog)
    
    '''ra = []
    dec = []
    mag = []
    for i in range(len(catalog.data['XWIN_WORLD'])):
        if catalog.data['MAG_AUTO'][i] + sci.header['MAGZP'] < 19:
            ra.append(catalog.data['XWIN_WORLD'][i])
            dec.append(catalog.data['YWIN_WORLD'][i])
            mag.append(catalog.data['MAG_AUTO'][i]+sci.header['MAGZP'])

        

    output = [ra,dec,mag]

    ascii.write(output,'/home/conor/.data/newdata/' + ipath[0:7] + 'truecoords.txt', delimiter = '\t', names=['RA','DEC','MAG'],overwrite = True)'''


    
    