# Store information for site, weather
from configobj import ConfigObj
import datetime, logging, ipdb
import ephem, time, math, os, sys, json, urllib2
class site:

    def __init__(self,site_name, night, configfile=''):

        self.name = site_name

        #set appropriate parameter based on aqawan_num
        #create configuration file object 
        configObj = ConfigObj(configfile)
        
        try:
            siteconfig = configObj[self.name]
        except:
            print('ERROR accessing ', self.name, ".", 
                self.name, " was not found in the configuration file", configfile)
            return 

        self.latitude = siteconfig['Setup']['LATITUDE']
        self.longitude = siteconfig['Setup']['LONGITUDE']        
        self.elevation = float(siteconfig['Setup']['ELEVATION'])
        self.enclosures = siteconfig['Setup']['ENCLOSURES']

        self.currentStatusFile = 'current_' + site_name + '.log'
        self.observing = True
        self.weather = -1
        self.startNightTime = -1
        self.night = night
        
        # touch a file in the current directory to enable cloud override
        self.cloudOverride = os.path.isfile('cloudOverride.txt') 
        self.sunOverride = os.path.isfile('sunOverride.txt')

        logger_name = siteconfig['Setup']['LOGNAME']
        log_file = 'logs/' + night + '/' + siteconfig['Setup']['LOGFILE']

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

 