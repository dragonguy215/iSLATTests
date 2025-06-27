import numpy as np

# Set-up default input parameters for model generation
min_lamb = 4.5
max_lamb = 28.
dist = 160.0
star_rv = 0.0
fwhm = 130.  # FWHM of the observed lines or instrument
pix_per_fwhm = 10  # number of pixels per fwhm element

intrinsic_line_width = 1.0
cc = 2.99792458e5  # speed of light in km/s
model_line_width = cc / fwhm
model_pixel_res = (np.mean ([min_lamb, max_lamb]) / cc * fwhm) / pix_per_fwhm

# Constants used in generating the rotation diagram
au = 1.496e11  # 1AU in m
pc = 3.08567758128e18  # From parsec to cm
ccum = 2.99792458e14  # speed of light in um/s
hh = 6.62606896e-27  # erg s

# Dictionary to store the initial values for each chemical
initial_values = {}

# define other defaults needed below
spanmol = "h2o"
specsep = .01  # default value for the separation to determine if line is single
fwhmtolerance = 5  # default value for the tolerance in FWHM for the de-blender (in km/s)
centrtolerance = 0.0001 # default value for the tolerance in centroid for the de-blender (in um)
line_threshold = 0.03  # percent value (where 0.01 = 1%) of the strongest line in the plot;
# lines below this this limit are ignored in the plot and in the single line selection