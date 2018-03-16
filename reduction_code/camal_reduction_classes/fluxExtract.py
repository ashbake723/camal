# CAMAL Reduction Class for CAMAL data
# 

import os
import matplotlib.pylab as plt
from photutils import datasets
from sigma_clipping import sigma_clip,sigma_clipped_stats
#from astropy.stats import sigma_clipped_stats #not sure why wont work
import numpy as np
from photutils import daofind, CircularAperture,CircularAnnulus,aperture_photometry
import astropy.io.fits as pyfits
import scipy
from scipy import signal
import camal
plt.ion()

__all__ = ['loadData','fluxExtract','sourceFinder'] 

    

def loadData(filename):
    ' Load data and return data array and header '
    dat = pyfits.open(filename)
    return dat[0].data, dat[0].header


def sourceFinder(science):
    ' Input: loaded data, Output: x,y center of aperture '
    # Convolve image with gaussian kernel
    kernel = np.outer(signal.gaussian(3,3), signal.gaussian(3,3))
    blurred = signal.fftconvolve(science, kernel, mode='same')
   
    # Take the normalized STD along x,y axes
    xstd = np.std(blurred,axis=0)
    ystd = np.std(blurred,axis=1)

    xsum = np.sum(blurred,axis=0)
    ysum = np.sum(blurred,axis=1)

    xstdn = (xstd - np.median(xstd))/max(xstd)
    ystdn = (ystd - np.median(ystd))/max(ystd)
    
    # Determine center by maximum. Eventually add check that there's only one source!
    try: x,y = np.where(xstdn == max(xstdn))[0][0], np.where(ystdn == max(ystdn))[0][0]
    except IndexError:
        x,y = 0,0
        
    return x,y



def fluxExtract(science,bias,dark,flats,hdr,plotap=False,rad=18,plotname=None):
    ' Returns final flux of source '
    # Calibrate Image. Flatten? add bias and flat to input
    flagg = 0
    
    # Defaul shape of sbig 2x2 binning
    xsize,ysize = 1266,1676
    if xsize != np.shape(bias)[0]:
        print 'WARNING: check that size of science file matches assumptions'

    if flats == 0:
        flat = np.ones((xsize,ysize))
    else:
        flat = flats[hdr['FILTER']]
    
    if dark == 0:
        dark = np.ones((xsize,ysize))
    else:
        dark = dark[str(hdr['EXPTIME'])]

    if np.mean(bias) == 0:
        bias = np.zeros((xsize,ysize))
    else:
        bias = bias

    # Size of science to trim calib frames (which are always full)
    # Science frames either full or sub: 400x400 around center ()
    centx, centy = int(1663/2), int(1252/2) # use config file in future, x and y are switched from thesky vs pyfits???
    subframe_size = np.shape(science)
    dx,dy         = subframe_size[0]/2,subframe_size[0]/2
    
    l,b,r,t = centx - dx, centy+dy, centx+dx, centy-dy # top and bottom switched

    data = (science - bias[t:b,l:r] - dark[t:b,l:r])/flat[t:b,l:r]
    
    if np.mean(data) < 0:
        flagg += 1
        
    # Get source x,y position
    x,y = sourceFinder(science)
    positions = (x,y)
    
    if x==0:
        flagg += 1
    
    if y==0:
        flagg += 1
        
    # Define Apertures
    apertures = CircularAperture(positions, r=rad)
    annulus_apertures = CircularAnnulus(positions, r_in=rad+5, r_out=rad+20)

    # Get fluxes
    rawflux_table = aperture_photometry(data, apertures)
    bkgflux_table = aperture_photometry(data, annulus_apertures)
    bkg_mean = bkgflux_table['aperture_sum'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    final_sum = rawflux_table['aperture_sum'] - bkg_sum

    # Plot
    if plotap == True:
        plt.ioff()
        plt.figure(-999)
        plt.clf()
        plt.imshow(np.log(data), origin='lower')
        apertures.plot(color='red', lw=1.5, alpha=0.5)
        annulus_apertures.plot(color='orange', lw=1.5, alpha=0.5)
        plt.savefig(plotname)

    return final_sum[0],flagg,x,y,bkg_mean[0]

