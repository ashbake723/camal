# Store information for reduction process, tools want to load once
from configobj import ConfigObj
import datetime, logging, ipdb
import ephem, time, math, os, sys, json, urllib2


class tonight:

    def __init__(self,site_name, night, configfile='camal_reduction_classes/reduce.ini'):
 
        # Reducation Things:
        self.night     = night
        self.dataPath = None
        self.finaldataPath = None
        self.schedulePath = None
        self.pwvPath   = None
        self.schedule   = None
        self.schedic   = None
        #self.mleft,self.mright = camal.getBinMatrix((2532, 3352))

        # Photometry things
        self.f780 = None
        self.f823 = None
        self.f860 = None
        self.c1   = None
        self.c2   = None
        self.alt  = None
        self.X    = None
        self.jd   = None
        self.hours= None
        self.iswitch = None
        self.objnames = None
        self.fnames = None
        self.Tgrid = None
        self.Xinds = None
        self.knots = None
        # ANalysis Things
        self.Temps  = None   #


        # Final things
        self.A1     = None
        self.A2     = None
        self.PWVs   = None
        self.PWVserr_lo = None
        self.PWVserr_hi = None
        self.pwvs_c1 = None
        self.pwvs_c2 = None
        self.p_best  = None
        self.samples = None
        self.m1best = None
        self.m2best = None

        
        # Site Things
        self.t_kitt    = None
        self.pwv_kitt  = None
        self.epwv_kitt = None
        
        # Load config file
        configObj  = ConfigObj(configfile)
        siteconfig = configObj[site_name]
        self.latitude   = siteconfig['Setup']['LATITUDE']
        self.longitude  = siteconfig['Setup']['LONGITUDE']        
        self.elevation  = float(siteconfig['Setup']['ELEVATION'])
        self.enclosures = siteconfig['Setup']['ENCLOSURES']

        self.currentStatusFile = 'current_' + site_name + '.log'
        
        logger_name = siteconfig['Setup']['LOGNAME']
        log_file    = 'logs/' + self.night + '/' + siteconfig['Setup']['LOGFILE']

