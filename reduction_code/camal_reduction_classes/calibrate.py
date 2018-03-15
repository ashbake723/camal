# Functions to Calibrate CAMAL Data

import numpy as np
import os
import astropy.io.fits as pyfits

import camal

class calibrate():
    """
    Calibration class to make calibration images for night
    """
    def __init__(self,path):
        self.path      = path
        self.medbias   = None
        self.flatfield = None
        self.medDark   = None

    def makeBias(self,biasfiles):
        """
        Create Median Bias files with biases found in night folder
    
        Note
        ----
        - To use older night's data, copy over the medBIAS.fits file
        - If no medBIAS.fits file found or if no biases in folder found
           then returns 0
        """
        if (np.size(biasfiles) > 0) & (not os.path.isfile(self.path + "\\medBIAS.fits")):
            nbias= 0
            for name in biasfiles:
                temp = pyfits.open(name)
                # Check CCD temperature
                if (abs(temp[0].header['CCD-TEMP'] - temp[0].header['CCD-TEMP'])< 0.3) & (nbias < 39):
                    try:
                        bias = np.dstack([bias,temp[0].data])
                    except NameError:
                        bias = temp[0].data
                    nbias += 1
                temp.close()
            medbias = np.median(bias ,axis=2)
            hdu=pyfits.PrimaryHDU(medbias)
            hdu.writeto(self.path + "\\medBIAS.fits")
            print 'Made Median Bias Frame and saved to ' + self.path + "\\medBIAS.fits"
            del bias #save room
            self.medbias = medbias
            return medbias

        elif os.path.isfile(self.path + "\\medBIAS.fits"):
            medbiasfits = pyfits.open(self.path + "\\medBIAS.fits")
            medbias = medbiasfits[0].data
            print 'Loaded Already Made Median Bias Frame'
            self.medbias = medbias
            return medbias
        else:
            print 'No Biases Found and no pre-made Median Bias Frame Found, \
                   returning 0'
            return 0


    def makeFlat(self,flatfiles):
        """
        Create the median, normalized Flat fields for each filter
    
        Notes:
        -----
        - If none exists, check if Blue one exists and then load and return all saved flats
          Must copy over these from another night...haven't implemented an automatic copying thing
        - If none exists, return 1
        """
        if (np.size(flatfiles) > 0) & (not os.path.isfile(self.path + "\\medFLAT_Blue.fits")):
            flats = {}
            flats['Blue'] = []
            flats['Luminance'] = []
            flats['Filter 5'] =[]
            for name in flatfiles:
                temp = pyfits.open(name)
                tempdat = temp[0].data
                temphdr = temp[0].header
                if temphdr['XBINNING'] == 1:
                    mleft,mright = camal.getBinMatrix(size)
                    tempdat = camal.rebin(temp[0].data,mleft = mleft,mright=mright)
                flats[temphdr['FILTER']].append((tempdat - self.medbias)/np.median(tempdat - self.medbias))
                print 'Added %s to flats for filter %s' %(name,temphdr['FILTER'])
            medflats = {}
            for keys in flats:
                flats[str(keys)] = np.median(flats[str(keys)],axis=0)
                hdu=pyfits.PrimaryHDU(flats[keys])
                hdu.writeto(self.path + "\\medFLAT_" + keys + ".fits")
                print 'Saved median flat for filter %s' % keys
            self.flatfield = flats
            return flats
        elif os.path.isfile(self.path + "\\medFLAT_Blue.fits"):
            flats = {}
            flats['Luminance'] = pyfits.open(self.path + "\\medFLAT_Luminance.fits")[0].data
            flats['Blue'] = pyfits.open(self.path + "\\medFLAT_Blue.fits")[0].data
            flats['Filter 5'] =pyfits.open(self.path + "\\medFLAT_Filter 5.fits")[0].data
            print 'Loaded Already Made Median Flat Fields'
            self.flatfield = flats
            return flats
        else:
            flats = 0
            return flats


    def makeDark(self,darkfiles):
        """
        Create median dark field for night
    
        Notes:
        -----
        - Right now, don't use darks. Have one for the polaris data.
    
        ToDo:
        ----
        - Implement check on the exposure time
        """
        if (np.size(darkfiles) > 0) & (not os.path.isfile(self.path + "medDark_2.3.fits")):
            dark = {}
            dark['2.3'] = []
            dark['5.0'] = []
            dark['2.1'] = []
            for name in darkfiles:
                temp = pyfits.open(name)
                exptime = temp[0].header['EXPTIME']
                dark[str(exptime)].append(temp[0].data - self.medbias)
            darks = {}
            for key in darks:
                darks[key] = np.median(dark[key],axis=0)[0]
                hdu=pyfits.PrimaryHDU(meddark[key])
                hdu.writeto(self.path + "medDARK_" + key + '.fits')
            return darks
        elif os.path.isfile(self.path + "medDark_2.3.fits"):
            darks = {}
            darks['2.3'] = pyfits.open(self.path + "medDARK_2.3.fits")[0].data
            darks['5.0'] = pyfits.open(self.path + "medDARK_5.0.fits")[0].data
            darks['2.1'] = pyfits.open(self.path + "medDARK_2.1.fits")[0].data
            print 'Loaded Already Made Median Darks'
            return darks
        else:
            print 'No Dark files found, nor existing file so returning 0'
            darks = 0
            return darks


def badPixels():
    """
    Define bad/hot pixels based on 2x2 binned pixel locations
    - Ashley picked these out pretty much by eye, finding the bad
       ones by looking at plots of convolved slices of the data.
    - Can try using a more reproducible way..
    """
    badpix_x = np.array([234,423,965, 1492, 1525, 1364])
    badpix_y = np.array([881, 1205, 808, 1043, 829, 746])
    return badpix_x, badpix_y



def badPixels_subframe():
    """
    need to do badpix of whole frame and then apply the
    subframe at the time. this is just a quick fix
    """
    badpix_x = np.array([274,281,334])
    badpix_y = np.array([217,262, 381])
    return badpix_x, badpix_y

