# import fits image and save as cropped image then delete old one

import matplotlib.pylab as plt
import numpy as np
import pyfits
from astropy.io import fits
from astropy import wcs
import os

# TO ADD #
# first estimate memory needed if want to rename them first then delete
# add change of binning
# make these data reduction things universal in class
# Load fits file                                                                                                          
night    = '061116'
objname  = 'polaris'
path2raw = '../raw/' + night + '/'
data = pyfits.open('%sfinalData_%s_%s.fits' %(path2raw,night,objname))[1].data

# Load Calibrations
if os.path.isfile(path2raw + "medFLAT_Blue.fits"):
    bias = medbiasfits = pyfits.open(path2raw + "medBIAS.fits")[0].data
else:

    bias = 0

if os.path.isfile(path2raw + "medFLAT_Blue.fits"):
    flats = {}
    flats['Luminance'] = pyfits.open(path2raw + "medFLAT_Luminance.fits")[0].data
    flats['Blue'] = pyfits.open(path2raw + "medFLAT_Blue.fits")[0].data
    flats['Filter 5'] = pyfits.open(path2raw + "medFLAT_Filter 5.fits")[0].data
else:
    flats = 1

xcent = data['Xcenter']
ycent = data['Ycenter']
rad = 150 # pixels - w, l for cropping box

datafiles = np.loadtxt(path2raw.replace('"','') + objname + 'list.txt',dtype='str')
i = 0
for f in datafiles:
    science = pyfits.open(path2raw + f)
    hdr = science[0].header
    xdim , ydim = np.shape(science[0].data)
    # Define new cutout being wary of the edges
    ylims = (0,rad) if ycent[i] < rad else (ycent[i] - rad, ycent[i] + rad)
    xlims = (0,rad) if xcent[i] < rad else (xcent[i] - rad, xcent[i] + rad)
    ylims = (ydim-2*rad,ydim) if ycent[i] + rad > ydim else  (int(ycent[i]) - rad, int(ycent[i]) + rad)
    xlims = (xdim-2*rad,xdim) if xcent[i] + rad > xdim else (int(xcent[i]) - rad, int(xcent[i]) + rad)
    # Trim & Note of it
    science[0].data = (science[0].data - bias)/flats[hdr['Filter']]
    science[0].data = science[0].data[ylims[0]:ylims[1],xlims[0]:xlims[1]]
    science[0].header['Notes'] = 'Trimmed: [%s:%s,%s:%s]' %(ylims[0],ylims[1],xlims[0],xlims[1])
    # Save (aka replace)
    if os.path.exists(path2raw + f) : os.remove(path2raw + f)
    science.writeto(path2raw + f)
    i+=1
    print f
