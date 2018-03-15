import time
from win32com.client import Dispatch
import logging
import subprocess 
import sys
import mail

ERROR = True
NOERROR = False
 
##------------------------------------------------------------------------------
## Class: cCamera
##------------------------------------------------------------------------------
class cSBIGtheskyX:
    def __init__(self,night):
        # CCD control through the SkyX
        # Define attributes
        # hard code things since CAMAL is simple. can put into config file later
        self.dataPath  = "C:/camal/data/" + night + "/"
        self.finaldataPath =  "C:/camal/final/" + night + "/"
        self.night     = night
        self.name      = 'SBIG'
        self.program   = 'theskyX'
        self.filters   = {
                        '823'      : 'Blue',
                        '860'      : 'Clear',
                        '780'      : 'Lunar',
                        '0'        : 'Red'
                        }
        self.filterindices   = {
                        'Blue'      : 2,
                        'Clear'     : 3,
                        'Lunar'     : 4,
                        'Red'       : 0
                        }
        self.frametype = {
                        'LIGHT'     :  1,
                        'BIAS'      :  2,
                        'DARK'      :  3,
                        'FLAT'      :  4
                        }
        self.xbin = 2
        self.ybin = 2
        
        self.asynchronous = False

        # set up logger
        logger_name = self.name
        log_file = 'logs/' + night + '/' + self.name
            
        # setting up imager logger
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
        
  

    def connect(self, cooler=True):
        # Connect to an instance of Maxim's camera control.
        # (This launches the app if needed)
        self.logger.info('Connecting to theSkyX')
        thesky = Dispatch("theSkyX.Application") #should already be open bc of launching telescope..
        self.CAMERA = Dispatch("CCDSoft2XAdaptor.CCDSoft5Camera")
        self.IMAGE  = Dispatch("CCDSoft2XAdaptor.ccdsoft5Image")
        # Connect to the camera 
        self.logger.info('Connecting to camera')
        try:
            self.CAMERA.Connect()
            self.nfailed = 0
        except:
            self.recover()

        # Set binning
        self.CAMERA.BinX = self.xbin
        self.CAMERA.BinY = self.ybin
        self.logger.info("Camera binning set to %dx%d" % (self.xbin,self.ybin))

        # Asynchronous? aka return before command complete
        self.CAMERA.Asynchronous  = self.asynchronous
        # Set to full frame
        self.setFrame('full')

        # Turn AutoSave off
        self.CAMERA.AutoSaveOn = False

        # Cool CDD
        if cooler == True:
            self.coolCCD()  

    def Expose(self,exptime,exptype,filterlam):
        filt  = self.filterindices[self.filters[str(filterlam)]]
        if exptype == 1 or exptype == 4:
            # Science or flat image
            self.logger.info("Exposing light frame...")
            self.CAMERA.FilterIndexZeroBased = filt
            self.CAMERA.ExposureTime         = exptime
            self.CAMERA.Frame                = exptype
            # Take Image
            self.CAMERA.TakeImage()
            # Attach this image to IMAGE object for saving
            self.IMAGE.AttachToActive()
            while not self.CAMERA.IsExposureComplete:
                time.sleep(1)
            self.logger.info("Light frame exposure and download complete!")
        
        if exptype == 2 or exptype == 3:
            # Bias or dark image, shutter closed
            self.logger.info("Exposing Dark frame...")
            self.CAMERA.ExposureTime         = exptime
            self.CAMERA.Frame                = exptype
            self.CAMERA.TakeImage()
            # Attach this image to IMAGE object for saving
            self.IMAGE.AttachToActive()
            while not self.CAMERA.IsExposureComplete:
                time.sleep(1)
            self.logger.info("Dark frame exposure and download complete!")

    def saveImage(self,filename):
        try:
            self.IMAGE.Path = filename
            self.IMAGE.Save()
            self.logger.info("saved file to " + filename)
        except:
            self.logger.error("Cannot save file")
            raise EnvironmentError, 'Halting program'


    def setFrame(self,opt,l=0,b=0,r=0,t=0):
        if opt == 'full':
            self.CAMERA.Subframe = 0
            self.logger.info("Camera set to full-frame mode")
        elif opt == 'sub':
            "Left 0 , Bottom 0 , Top, Right"
            " Must give unbinned pixel coordinates "
            self.CAMERA.SubframeLeft = l
            self.CAMERA.SubframeRight = r
            self.CAMERA.SubframeBottom = b
            self.CAMERA.SubframeTop = t
            self.CAMERA.Subframe = 1           
            self.logger.info("Camera set to subframe. Note: Inputs must ignore binning.")
        else:
            self.logger.error("Set opt to either 'full' or 'sub'. If sub, give coordinates")

    def coolCCD(self):
        # Determine SetPoint Temp if cooler not on, make divisible by 5
        if self.CAMERA.RegulateTemperature == 0:
            ambTemp = self.CAMERA.Temperature
            if ambTemp < 15:
                settemp = -5
            elif (ambTemp >= 15) & (ambTemp < 20):
                settemp = 0
            elif (ambTemp >= 20) & (ambTemp < 25):
                settemp = 5
            else:
                settemp = round(ambTemp-20) - round(ambTemp)%5
        else: # make a set temp if cooler already on
            ambTemp = round(self.CAMERA.Temperature)
            if self.CAMERA.ThermalElectricCoolerPower > 30:
                # sign: +/- : Temperature larger/smaller than set point
                sign = (self.CAMERA.Temperature - self.CAMERA.TemperatureSetPoint)/abs(self.CAMERA.Temperature - self.CAMERA.TemperatureSetPoint)
                settemp = round(self.CAMERA.TemperatureSetPoint + sign * 5)
            else:
                settemp = round(self.CAMERA.TemperatureSetPoint)
        # Set Temp and turn on cooling
        self.CAMERA.TemperatureSetPoint = settemp
        self.CAMERA.RegulateTemperature = 1 # Turn on cooling
        time.sleep(5)
        self.logger.info("Cooling camera to " + str(settemp) + " C , Amb Temp= " + str(ambTemp))
        self.logger.info("Waiting for Camera Temperature to match Set Temp to within 0.4 Celcius" )
        time0 = time.time()
        it = 0
        while (abs(self.CAMERA.Temperature - self.CAMERA.TemperatureSetPoint) > 0.4) or (self.CAMERA.ThermalElectricCoolerPower > 50):
            time.sleep(10)   # sleep until temp and settemp match
            if time.time() - time0 > 600.0:
                # If Temperature still too high (low)..raise (lower) settemp
                sign = (self.CAMERA.Temperature - self.CAMERA.TemperatureSetPoint)/abs(self.CAMERA.Temperature - self.CAMERA.TemperatureSetPoint)
                settemp = settemp + sign * 5
                self.logger.info("Cooling Took Too Long, changing settemp to " +str(settemp))
                self.CAMERA.TemperatureSetPoint = settemp
                time0 = time.time()
                self.logger.warning('Camera failed to Cool! Time: ' + str(it))
                it += 1
                if it > 3:
                    self.logger.error('Camera failed to cool!')
                    mail.send("SBIG Camera failed to Cool","please do something")
                    sys.exit()



    def shutDown(self):
        if self.CAMERA.RegulateTemperature == 1:
            # Warm up cooler
            self.CAMERA.TemperatureSetPoint = self.CAMERA.Temperature + 10.0
            self.logger.info('Warming Up CCD to Amb. Temp.')
            time.sleep(25)
            # Turn Cooler Off
            self.CAMERA.RegulateTemperature = 0
        # Quit from camera, disconnect
        self.logger.info('Disconnecting and Quitting CAMERA')
        self.CAMERA.Disconnect()

    def restartthesky(self):
        self.logger.info('Killing maxim') 
        subprocess.call(['Taskkill','/IM','theSkyX.exe','/F'])

        time.sleep(15)

        self.logger.info('Reconnecting')
        self.connect()


    def recover(self):
        self.nfailed = self.nfailed + 1

        try:
            self.shutDown()
        except:
            pass

        if self.nfailed == 1:
            # attempt to reconnect
            self.logger.warning('Camera failed to connect; retrying') 
            self.connect()
        elif self.nfailed == 2:
            # then restart maxim
            self.logger.warning('Camera failed to connect; restarting maxim') 
            self.restartthesky()
        elif self.nfailed == 3:
            self.logger.error('Camera failed to connect!') 
            mail.send("SBIG Camera failed to connect","please do something",level="serious")
            sys.exit()




class cSBIGmaxim:
    def __init__(self,night):
        # Define attributes
        # hard code things since CAMAL is simple. can put into config file later
        self.dataPath  = None
        self.finaldataPath = None
        self.night     = night
        self.name      = 'SBIG'
        self.program   = 'maxim'
        self.filters   = {
                        '823'      : 'Luminance',
                        '860'      : 'Blue',
                        '780'      : 'Filter 5',
                        '0'        : 'Red'
                        }
        self.filterindices   = {
                        'Luminance'   : 2,
                        'Blue'        : 3,
                        'Filter 5'    : 4,
                        'Red'         : 0
                        }
        self.xbin = 2
        self.ybin = 2
        
        # set up logger
        logger_name = self.name
        log_file = 'logs/' + night + '/' + self.name
            
        # setting up imager logger
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

    def Expose(self,exptime,exptype,filterlam):
        """
        Take exposure time, exosure type (0 or 1), and 
        filter wavelength as string (780, 823, 860, or 0 for none)
        """

        filt = self.filterindices[self.filters[filterlam]]
        if exptype == 1:
            # Science or flat image
            self.logger.info("Exposing light frame...")
            self.CAMERA.Expose(exptime,1,filt)
            while not self.CAMERA.ImageReady:
                time.sleep(1)
            self.logger.info("Light frame exposure and download complete!")
        
        if exptype == 0:
            # Bias or dark image, shutter closed
            self.logger.info("Exposing Dark frame...")
            self.CAMERA.Expose(exptime,0,-1)
            while not self.CAMERA.ImageReady:
                time.sleep(1)
            self.logger.info("Dark frame exposure and download complete!")
        
        
    def setFrame(self,opt,l=0,b=0,r=0,t=0):
        """
        """
        x = l
        y = t
        wx = r-l
        wy = b-t
        if opt == 'full':
            self.CAMERA.SetFullFrame()
            self.logger.info("Camera set to full-frame mode")
        elif opt == 'sub':
            self.CAMERA.StartX = x
            self.CAMERA.StartY = y
            self.CAMERA.NumX = wx 
            self.CAMERA.NumY = wy 
            self.logger.info("Camera set to subframe. Note: Inputs must consider binning.")
        else:
            self.logger.error("Set opt to either 'full' or 'sub'")
        

    def setBinning(self,binmode):
        tup = (1,2,4)
        if binmode in tup:
            self.CAMERA.BinX = binmode
            self.CAMERA.BinY = binmode
            self.logger.info("Camera binning set to %dx%d" % (binmode,binmode))
        else:
            self.logger.error("ERROR: Invalid binning specified")


    def saveImage(self,filename):
        try:
            self.CAMERA.SaveImage(filename)
            self.logger.info("saved file to " + filename)
        except:
            self.logger.error("Cannot save file")
            raise EnvironmentError, 'Halting program'


 
    def coolCCD(self):
        if not self.CAMERA.CanSetTemperature:
            self.logger.error("ERROR: cannot change temperature")
        else:
            self.CAMERA.CoolerOn = True
            time.sleep(0.5)
            ambTemp = self.CAMERA.Temperature
            if ambTemp < 15:
                settemp = -5
            elif (ambTemp >= 15) & (ambTemp < 20):
                settemp = 0
            elif (ambTemp >= 20) & (ambTemp < 25):
                settemp = 5
            else:
                settemp = ambTemp - 20
                # logging.info('Temperature too hot..not cooling camera')

            self.CAMERA.TemperatureSetpoint = settemp
            time.sleep(5)
            self.logger.info("Cooling camera to " + str(settemp) + " C , Amb Temp= " + str(ambTemp))
            self.logger.info("Waiting for Cooler Power to Stabalize Below 35%" )
            tt = time.time()
            fails = 0
            while self.CAMERA.CoolerPower > 35:
                time.sleep(1)
                if time.time()  - tt > 300:
                    settemp += 5
                    self.CAMERA.TemperatureSetpoint = settemp
                    tt = time.time()
                    fails += 0
                    print 'failed to cool once'
                if fails > 2:
                    mail.send('CAMERA Cooler Power Never Settled','The temperature of the camera never'
                        'settled to a cooler power less than 35%. Check out why'
                        'Continuing anyways.')
                    print 'temperature didnt reach setpoint'
                    break
            tt = time.time()
            while abs(self.CAMERA.Temperature - self.CAMERA.TemperatureSetpoint) > 0.4:
                time.sleep(1)   # sleep more because it usually overshoots
                if time.time() - tt > 300:
                    mail.send('CAMERA NOT SETTLED','The temperature of the camera never'
                        'settled to within 0.4 degrees C within the setpoint. Check out why'
                        'Continuing anyways.')
                    break


    def shutDown(self):
        if self.CAMERA.CoolerOn:
            # Warm up cooler
            if self.CAMERA.TemperatureSetpoint < self.CAMERA.AmbientTemperature:
                self.CAMERA.TemperatureSetpoint = self.CAMERA.AmbientTemperature
                self.logger.info('Warming Up CCD to Amb. Temp.')
                time.sleep(25)
            # Turn Cooler Off
            self.CAMERA.CoolerOn = False
        # Quit from camera, disconnect
        self.logger.info('Disconnecting and Quitting CAMERA')
        self.CAMERA.Quit()

    def restartmaxim(self):
        self.logger.info('Killing maxim') 
        subprocess.call(['Taskkill','/IM','MaxIm_DL.exe','/F'])

        time.sleep(15)

        self.logger.info('Reconnecting')
        self.connect()


    def recover(self):

        self.nfailed = self.nfailed + 1

        try:
            self.shutDown()
        except:
            pass

        if self.nfailed == 1:
            # attempt to reconnect
            self.logger.warning('Camera failed to connect; retrying') 
            self.connect()
        elif self.nfailed == 2:
            # then restart maxim
            self.logger.warning('Camera failed to connect; restarting maxim') 
            self.restartmaxim()
        elif self.nfailed == 4:
            self.logger.error('Camera failed to connect!') 
            mail.send("Camera " + str(self.num) + " failed to connect","please do something",level="serious")
            sys.exit()

    def connect(self, cooler=True):
        settleTime = 1200
        oscillationTime = 120.0

        # Connect to an instance of Maxim's camera control.
        # (This launches the app if needed)
        self.logger.info('Connecting to Maxim') 
        self.CAMERA = Dispatch("MaxIm.CCDCamera")

        # Connect to the camera 
        self.logger.info('Connecting to camera')
        try:
            self.CAMERA.LinkEnabled = True
            self.nfailed = 0
        except:
            self.recover()

        # Prevent the camera from disconnecting when we exit
        self.logger.info('Preventing the camera from disconnecting when we exit') 
        self.CAMERA.DisableAutoShutdown = True

        # If we were responsible for launching Maxim, this prevents
        # Maxim from closing when our application exits
        self.logger.info('Preventing maxim from closing upon exit')
        maxim = Dispatch("MaxIm.Application")
        maxim.LockApp = True

        # Set binning
        self.setBinning(self.xbin)

        # Set to full frame
        self.setFrame('full')

        # Cool CDD
        if cooler == True:
            self.coolCCD()  
            
##
##    END OF 'cCamera' Class
##
