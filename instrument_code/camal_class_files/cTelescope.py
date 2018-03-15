import time
import win32com.client
import logging
import mail

ERROR = True
NOERROR = False
 
 #http://www.rkblog.rk.edu.pl/w/p/ascom-end-user-application-developers/
 
##------------------------------------------------------------------------------
## Class: cTelescope
##------------------------------------------------------------------------------

__all__ =["cMyT"]

class cMyT:
    def __init__(self,night):
        self.name  = 'MyT'
        self.night = night
        self.connected = False
        self.initialized = False

        logger_name = self.name
        log_file = 'logs/' + night + '/' + logger_name
        self.currentStatusFile = 'current_' + self.name + '.log'

        # setting up telescope logger
        fmt = "%(asctime)s [%(filename)s:%(lineno)s - %(funcName)s()] %(levelname)s: %(message)s"
        datefmt = "%Y-%m-%dT%H:%M:%S"

        self.logger = logging.getLogger(logger_name)
        formatter = logging.Formatter(fmt,datefmt=datefmt)
        formatter.converter = time.gmtime
        
        fileHandler = logging.FileHandler(log_file, mode='a')
        fileHandler.setFormatter(formatter)

        console = logging.StreamHandler()
        console.setFormatter(formatter)
        console.setLevel(logging.INFO)
        
        self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(fileHandler)
        self.logger.addHandler(console)

        self.logger.info("Connecting to The Sky...")
        self.appl = win32com.client.Dispatch("TheSkyX.Application")
        self.__CHOOSER = win32com.client.Dispatch("ASCOM.Utilities.Chooser")
        self.__CHOOSER.DeviceType = 'Telescope'
        # Use this class to control the telescope
        self.tel = win32com.client.Dispatch("ASCOM.SoftwareBisque.Telescope") #make it __TELE once done testing       
 
        # Use this class to get the pointing model corrected sky coords
        self.skychart = win32com.client.Dispatch("TheSkyXAdaptor.RASCOMTheSky")


    def initialize(self):
        """
        Unpark and home telescope
        """
        # Dont need these anymore
        #self.CLS = win32com.client.Dispatch("TheSkyX.ClosedLoopSlew")
        #self.SC = win32com.client.Dispatch("TheSkyXAdaptor.StarChart")
        #self.OI = win32com.client.Dispatch("TheSkyXAdaptor.ObjectInformation")
 
        # Connect to telescope
        if self.tel.Connected:
            self.logger.info(" ->Telescope was already connected")
        else:
            self.tel.Connected = True
        if self.tel.Connected:
            self.logger.info(" Connected to telescope now")
        else:
            self.logger.error(" Unable to connect to telescope, expect exception")
 
        # Unpark
        if self.tel.AtPark:
            if self.tel.CanUnpark:
                self.tel.Unpark()
                self.logger.info('Telescope Unparked')
            else:
                self.logger.info('Telescope cannot unpark')

        # Home the scope
        self.home()

        self.initialized = True


    def gotoObject(self,objname):
        # Find object coordinates (PM corrected)
        try:
            self.skychart.GetObjectRaDec(objname) # load objname in sky chart
        except NameError:
            pass #eh just don't get the name wrong

        RA , DEC = self.skychart.dObjectRA, self.skychart.dObjectDEC

        # Slew to object coordinates
        if self.tel.CanSlew:
            self.tel.SlewToCoordinates(RA,DEC)     # !!!pick coords currently in sky!
            self.logger.info('Slewing to RA: %s, DEC: %s' %(RA,DEC))
        else:
            self.logger.warning('Warning: Cannot move axis. Did not slew')
        while(self.tel.Slewing):
            time.sleep(1)



    def goto(self,RA, DEC):
        """
        Go to RA DEC coordinates described in degrees
        """
        if self.tel.CanSlew:
            self.tel.SlewToCoordinates(RA,DEC)     # !!!pick coords currently in sky!
            self.logger.info('Slewing to RA: %s, DEC: %s' %(RA,DEC))
        else:
            self.logger.warning('Warning: Cannot move axis. Did not slew')
        if self.tel.CanSync:
            pass #eventually sync? make sure star is in center?
        while(self.tel.Slewing):
            time.sleep(1)

    def getRA(self):
        return self.tel.RightAscension

    def getDEC(self):
        return self.tel.Declination

    def park(self):
        self.tel.Park()
        while not self.tel.AtPark:
            time.sleep(1)
        if self.tel.AtPark:
            self.logger.info('Telescope Parked')

    def home(self):
        self.tel.FindHome()
        self.logger.info('Telescope Finding Home Position')
        while(self.tel.Slewing):
            time.sleep(1)  # adding this in case asynchronous is set


    def shutDown(self):
        # Park Telescope
        self.logger.info('Parking telescope')
        if not self.tel.AtPark:
            if self.tel.Connected == True:
                self.tel.Park()
                while not self.tel.AtPark:
                    time.sleep(1)
            if self.tel.Connected == False: # didn't connect for night
                self.tel.Connected = True
                self.tel.Unpark()
                self.home()
                self.tel.Park()
        if self.tel.AtPark:
            self.logger.info('Telescope Parked')
        else:
            self.logger.error('Telescope did not park')
            mail.send('CAMAL Telescope Did not Park!','Please log into CAMALPC and check on things')
        # Disconnect
        self.tel.Connected = False
        if not self.tel.Connected:
            self.logger.info('Disconnected from Telescope but not in the SkyX')



