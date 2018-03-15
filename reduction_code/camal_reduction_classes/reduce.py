# Store information for reduction process, tools want to load once
from configobj import ConfigObj
import datetime, logging, ipdb
import ephem, time, math, os, sys, json, urllib2
import calibrate, camal

class tonight:

    def __init__(self,site_name, night, configfile='camal_reduction_classes/reduce.ini'):
 
        # Reducation Things:
        self.night     = night
        self.dataPath = None
        self.finaldataPath = None
        self.objects   = None
        self.schedic   = None
        self.badpix_x, self.badpix_y = calibrate.badPixels_subframe()
        self.mleft,self.mright = camal.getBinMatrix((2532, 3352))

        
        # Site Things
        # Load config file
        configObj  = ConfigObj(configfile)
        self.name  = site_name
        siteconfig = configObj[self.name]
        self.latitude   = siteconfig['Setup']['LATITUDE']
        self.longitude  = siteconfig['Setup']['LONGITUDE']        
        self.elevation  = float(siteconfig['Setup']['ELEVATION'])
        self.enclosures = siteconfig['Setup']['ENCLOSURES']

        self.currentStatusFile = 'current_' + site_name + '.log'
        
        logger_name = siteconfig['Setup']['LOGNAME']
        log_file    = 'logs/' + self.night + '/' + siteconfig['Setup']['LOGFILE']

        self.obs = ephem.Observer()
        self.obs.lat = ephem.degrees(str(self.latitude)) # N
        self.obs.lon = ephem.degrees(str(self.longitude)) # E
        self.obs.elevation = self.elevation # meters

    def sunrise(self, horizon=0):

        self.obs.horizon = str(horizon)
        sunrise = self.obs.next_rising(ephem.Sun(), start=self.startNightTime, use_center=True).datetime()
        return sunrise
    
    def sunset(self, horizon=0):

        self.obs.horizon = str(horizon)
        sunset = self.obs.next_setting(ephem.Sun(), start=self.startNightTime, use_center=True).datetime()
        return sunset
        
    def NautTwilBegin(self, horizon=-12):

        self.obs.horizon = str(horizon)
        NautTwilBegin = self.obs.next_rising(ephem.Sun(), start=self.startNightTime, use_center=True).datetime()
        return NautTwilBegin
    
    def NautTwilEnd(self, horizon=-12):

        self.obs.horizon = str(horizon)
        NautTwilEnd = self.obs.next_setting(ephem.Sun(), start=self.startNightTime, use_center=True).datetime()
        return NautTwilEnd

    def sunalt(self):

        self.obs.date = datetime.datetime.utcnow()
        sun = ephem.Sun()
        sun.compute(self.obs)
        return float(sun.alt)*180.0/math.pi

    def sunaz(self):

        self.obs.date = datetime.datetime.utcnow()
        sun = ephem.Sun()
        sun.compute(self.obs)
        return float(sun.az)*180.0/math.pi

    def moonpos(self):
        moon = ephem.Moon()
        moon.compute(datetime.datetime.utcnow())
        moonpos = (moon.ra,moon.dec)
        return moonpos
    
    def moonphase(self):
        moon = ephem.Moon()
        moon.compute(datetime.datetime.utcnow())
        moonphase = moon.phase/100.0
        return moonphase