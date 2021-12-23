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
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from scipy import stats
import pandas as pd
import fakes


datadir = '/home/conor/.data/newdata/'

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


sci = zuds.ScienceImage.from_file(refimpaths[1])

ra_g = []
dec_g = []

rag = sci.gaia_dr2_calibrators()['ra']
decg = sci.gaia_dr2_calibrators()['dec']
for i in range(len(rag)):
    ra_g.append(rag[i])
    dec_g.append(decg[i])



x_removed = []
mag_removed = []
x1_removed = []
mag1_removed = []

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
    '''catalogs.append(catalog)'''
    
    
    ra = []
    dec = []
    mag = []
    
    for i in range(len(catalog.data['XWIN_WORLD'])):
        ra.append(catalog.data['XWIN_WORLD'][i])
        dec.append(catalog.data['YWIN_WORLD'][i])
        mag.append(catalog.data['MAG_AUTO'][i]+sci.header['MAGZP'])

    x = []

    stars = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
    idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=ra_g*u.degree,dec=dec_g*u.degree),0.000277778*u.degree)
    

    ra_matched = []
    dec_matched = []
    ra_bmatched = []
    mag_matched = []
    for j in range(len(sep2d.arcsecond)):
        x.append(sep2d.arcsecond[j])
        ra_matched.append(ra[idxself[j]])
        dec_matched.append(dec[idxself[j]])
        mag_matched.append(mag[idxself[j]])

    gdindex = []
    for i in range(len(x)):
        if abs(x[i] - np.median(x)) < 5 * 1.4826 * stats.median_absolute_deviation(x-np.median(x)):
            gdindex.append(idxself[i])
        else:
            x_removed.append(x[i])
            mag_removed.append(mag[idxself[i]])
    gdindex.sort()
    out = catalog.data[gdindex]



    catalog.load()
    catalog.data = out

    catalog.save()
    catalog.load()

    # associate each mask image with a science imag



    zuds.calibrate_astrometry(sci,inplace = True)
    # add the objects to the database session
    #zuds.DBSession().add_all([sci, msk])

    catalog.load()
    catalog = zuds.PipelineFITSCatalog.from_image(sci)
    catalog.save()
    catalog.load()


    scis.append(sci)

    
    ra = []
    dec = []
    mag = []
    
    for i in range(len(catalog.data['XWIN_WORLD'])):
        ra.append(catalog.data['XWIN_WORLD'][i])
        dec.append(catalog.data['YWIN_WORLD'][i])
        mag.append(catalog.data['MAG_AUTO'][i]+sci.header['MAGZP'])

    x = []

    stars = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
    idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=ra_g*u.degree,dec=dec_g*u.degree),0.000277778*u.degree)
    

    ra_matched = []
    dec_matched = []
    ra_bmatched = []
    mag_matched = []
    for j in range(len(sep2d.arcsecond)):
        x.append(sep2d.arcsecond[j])
        ra_matched.append(ra[idxself[j]])
        dec_matched.append(dec[idxself[j]])
        mag_matched.append(mag[idxself[j]])


    gdindex = []
    for i in range(len(x)):
        if abs(x[i] - np.median(x)) < 5 * 1.4826 * stats.median_absolute_deviation(x-np.median(x)):
            gdindex.append(idxself[i])
        else:
            x1_removed.append(x[i])
            mag1_removed.append(mag[idxself[i]])
    gdindex.sort()
    out = catalog.data[gdindex]
    catalog.load()
    catalog.data = out

    catalog.save()
    catalog.load()

    zuds.calibrate_astrometry(sci,inplace = True)
    # add the objects to the database session
    #zuds.DBSession().add_all([sci, msk])
   
    catalog.load()
    catalog = zuds.PipelineFITSCatalog.from_image(sci)

    catalog.save()
    catalog.load()

    '''catalogs.append(catalog)'''
    
    
    ra = []
    dec = []
    mag = []
    
    for i in range(len(catalog.data['XWIN_WORLD'])):
        if catalog.data['MAG_AUTO'][i] + sci.header['MAGZP'] < 19:
            ra.append(catalog.data['XWIN_WORLD'][i])
            dec.append(catalog.data['YWIN_WORLD'][i])
            mag.append(catalog.data['MAG_AUTO'][i]+sci.header['MAGZP'])



    output = [ra,dec,mag]

    ascii.write(output,'/home/conor/.data/newdata/' + ipath[40:55] + 'truecoords.txt', delimiter = '\t', names=['RA','DEC','MAG'],overwrite = True)






file_list = os.listdir('/home/conor/.data/newdata/')


files = []
coords = []

for file in file_list:
    if 'truecoords.txt' in file:
        coords.append('/home/conor/.data/newdata/'+file)
    elif 'sciimg.fits' in file:
        files.append('/home/conor/.data/newdata/' + file)
coords.sort()
files.sort()

for coord,file in zip(coords,files):
    f = fits.open(file)
    hdu = f[0].header
    if hdu['MAGLIM'] >= 20.5 and hdu['SEEING'] <= 2.0:
        base_file = coord
        break


ra_base = []
dec_base = []
mag_base = []

data_base = ascii.read(base_file)
for i in range(len(data_base)):
    ra_base.append(data_base['RA'][i])
    dec_base.append(data_base['DEC'][i])
    mag_base.append(data_base['MAG'][i])


coords.remove(base_file)

ra = []
dec = []
mag = []


count = 0
for reg in coords:
    data = ascii.read(reg)


    for i in range(len(data['RA'])):
        ra.append(data['RA'][i])
        dec.append(data['DEC'][i])
        mag.append(data['MAG'][i])
    

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
    new_mags = []
    if i!=0 and np.abs(ra_matched[i] - ra_matched[i-1]) < 0.000277778:
        continue
    elif i == 0:
        continue
    else:
        for l in range(index,i):
            x_chunk.append(x[l])
            new_mags.append(mag_matched[l])
        
        if len(x_chunk) != 1:
            stds.append(np.median(x_chunk))
            newer_mags.append(np.median(new_mags))


        index = i


print(np.median(stds))



x2 = []
mag2 = []
stars = SkyCoord(ra=ra_base*u.degree,dec=dec_base*u.degree)
idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=ra_g*u.degree,dec=dec_g*u.degree),0.000277778*u.degree)
for i in range(len(sep2d.arcsecond)):
    x2.append(sep2d.arcsecond[i])
    mag2.append(mag_base[idxself[i]])



print(np.median(x2))
fig, ax = plt.subplots()
ax.add_artist(plt.scatter(newer_mags,stds))
#ax.add_artist(plt.scatter(mag_removed,x_removed,color = 'g'))
#ax.add_artist(plt.scatter(mag1_removed,x1_removed,color = 'y'))
ax.add_artist(plt.scatter(mag2,x2,color = 'r'))
plt.xlabel("Magnitude")
plt.ylabel("Combined RMS (arsec)")
plt.yscale('log')
plt.title('Scatter relative to ZTF vs relative to Gaia')
plt.show()
