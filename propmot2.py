import sys
import zuds
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.optimize import fmin_l_bfgs_b
from scipy.stats import linregress
from iminuit import Minuit
from iminuit.util import describe, make_func_code
from iminuit.cost import LeastSquares

PIX_SCALE = 1.012 * u.arcsec

# initialize the database
zuds.init_db()

# get the list of science images to read from the command line
image_list = sys.argv[1]

file_list = os.listdir(image_list)


files = []

for file in file_list:
    if 'sciimg.fits' in file:
        files.append('/home/conor/.data/newdata/' + file)
    

files.sort()


# load in image paths
image_paths = np.genfromtxt(files, encoding="ascii", dtype=None)

# load in images, filter out bad images, and sort by date (ascending)
images = list(
    filter(
        lambda img: img.header["MAGLIM"] > 20.5 and img.header["SEEING"] < 2.5,
        sorted(
            [zuds.ScienceImage.from_file(p) for p in image_paths],
            key=lambda img: img.mjd,
        ),
    )
)

# make catalogs
for i in images:
    i.mask_image = zuds.MaskImage.from_file(i.local_path.replace("sci", "msk"))
    i.catalog = zuds.PipelineFITSCatalog.from_image(i)

# define the base image as the chronologically first image, pop it off the stack
base_image = images[0]
images = images[1:]

# get proper motions from gaia
gaia = base_image.gaia_dr2_calibrators()

# keep only the objects in gaia that are in a reasonable magnitude range
gaia = gaia[(gaia["phot_rp_mean_mag"] > 15) & (gaia["phot_rp_mean_mag"] < 21)]

# match the base catalog to gaia
base_coords = SkyCoord(
    ra=base_image.catalog.data["XWIN_WORLD"],
    dec=base_image.catalog.data["YWIN_WORLD"],
    unit="deg",
)
gaia_coords = SkyCoord(ra=gaia["ra"], dec=gaia["dec"], unit="deg")
idx_gaia, idx_base, _, _ = base_coords.search_around_sky(gaia_coords, 1 * u.arcsec)

# prune the base and gaia catalogs, keeping only the matches in both catalogs
base_catalog = base_image.catalog.data[idx_base]
gaia_catalog = gaia[idx_gaia]
base_coords = base_coords[idx_base]

# match each catalog to the base catalog, storing matches in a dictionary of
# the following form:
#       base_star         : [(   matched_image,          matched_row         )]
#  base_catalog_row_index : [(image_list_index_other, other_catalog_row_index)]


# iterate over the remaining images
match_dictionary = {ind: [] for ind in range(len(base_coords))}
for img_ind, other in enumerate(images):
    match_coords = SkyCoord(
        ra=other.catalog.data["XWIN_WORLD"],
        dec=other.catalog.data["YWIN_WORLD"],
        unit="deg",
    )
    idx_other, idx_base, _, _ = base_coords.search_around_sky(
        match_coords, 1 * u.arcsec
    )

    for cat_idxother, cat_idxbase in zip(idx_other, idx_base):
        match_dictionary[cat_idxbase].append((img_ind, cat_idxother))

# fit a line to the motion of each star

pmra = []
pmraerr = []
ras = []

decs = []


def line(x, a, b):  # simple straight line model with explicit parameters
    return a + b * x

for base_index in match_dictionary:

    # fit RA and Dec separately, then combine
    for key in ["ra", "dec"]:
        sex_key = "XWIN_WORLD" if key == "ra" else "YWIN_WORLD"
        err_key = "ERRX2WIN_IMAGE" if key == "ra" else "ERRY2WIN_IMAGE"

        # get the ztf position and error along this axis for all matched stars
        y_data_deg = np.asarray(
            [
                images[img_idx].catalog.data[row_idx][sex_key]
                for img_idx, row_idx in match_dictionary[base_index]
            ]
        )

        var_pix2 = np.asarray(
            [
                images[img_idx].catalog.data[row_idx][err_key]
                for img_idx, row_idx in match_dictionary[base_index]
            ]
        )

        # get the mjd of the image of each matched star
        x = np.asarray(
            [images[img_idx].mjd for img_idx, _ in match_dictionary[base_index]]
        )

        # convert from variance to sigma
        sigma_deg = (np.sqrt(var_pix2) * PIX_SCALE).to("deg").value

        def objective_function(parameters):
            """The chi-squared of the model against the data. Minimize this."""
            slope, intercept = parameters
            y_mod_deg = slope * x + intercept
            chi = (y_mod_deg - y_data_deg) / sigma_deg
            return np.sum(chi * chi)

        # minimize the function
        guess = linregress(x, y_data_deg)
        minres = fmin_l_bfgs_b(objective_function, guess[:2], approx_grad=True)
        slope_fit, intercept_fit = minres[0]


        


        least_squares = LeastSquares(x, y_data_deg, sigma_deg,line)


        m = Minuit(least_squares, a=0,b=0)

        
        m.migrad()
        m.hesse()


        if key == 'ra':
            ras.append(np.median(y_data_deg))
            pmra.append(m.values['b']*365*1000*3600)
            pmraerr.append(m.errors['b']*365*1000*3600)
        else:
            decs.append(np.median(y_data_deg))


        if np.median(y_data_deg) == 257.90185999651567:
            plt.errorbar(x,y_data_deg,yerr = sigma_deg,ls = 'none')
            plt.plot(x,[slope_fit*x[i] + intercept_fit for i in range(len(x))])
            plt.title('Proper Motion of Example Star After Filter')
            plt.ylabel('Offset in RA (degrees)')
            plt.xlabel('MJD of Images')
            plt.show()







gmra_matched = []
pmra_matched = []
gmraerr_matched = []
pmraerr_matched = []
sn_matched = []
mag_matched = []



stars = SkyCoord(ra=ras*u.degree,dec=decs*u.degree)
idxs,idxself,sep2d,dist3d = stars.search_around_sky(SkyCoord(ra=gaia_catalog['ra'],dec=gaia_catalog['dec'], unit = 'deg'),0.000277778*u.degree)
for i in range(len(idxs)):
    pmra_matched.append(pmra[idxself[i]])
    pmraerr_matched.append(pmraerr[idxself[i]])
    gmra_matched.append(gaia_catalog['pmdec'][idxs[i]])
    gmraerr_matched.append(gaia_catalog['pmdec_error'][idxs[i]])

model = []

for i in range(len(pmra_matched)):
    model.append((pmra_matched[i]-gmra_matched[i])/pmraerr_matched[i])


count = 0
model1 = []
gmra_matched1 = []
for i in range(len(model)):
    if np.abs(model[i]) < 2:
        count +=1
    if np.abs(model[i]) < 1000:
        model1.append(model[i])
        gmra_matched1.append(gmra_matched[i])
        
print(count/len(idxs))



#plt.errorbar(gmra_matched,model,yerr = pmraerr_matched,xerr = gmraerr_matched,fmt = 'o',ecolor = 'b',capthick=1,capsize = 2)
plt.scatter(gmra_matched1,model1)
#ax.add_artist(plt.scatter(mag2,x2,color = 'r'))
#plt.plot(gmra_matched,gmra_matched)
plt.xlabel("Gaia Proper Motion in RA (mas/yr)")
plt.ylabel("Simga Offset from True Values")
plt.title('Sigma Difference Between Filtered Results and Gaia')
plt.show()
