

import camal_class_files.mysite as camalsite
import camal_class_files.cCamera as camalimager
import camal_class_files.cTelescope as camaltelescope
import camal_class_files.cSchedule as camalscheduler
import camal_class_files.mail as mail
import numpy as np
from camal_class_files.get_centroid import *
import shutil, re

import datetime, logging, os, sys, time, subprocess, glob, math, json, copy
import ipdb
import socket, threading, ephem
import astropy.io.fits as pyfits
from scipy import stats, signal
import urllib

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive


def ten(string):
    array = string.split(':')
    if "-" in array[0]:
        return float(array[0]) - float(array[1])/60.0 - float(array[2])/3600.0
    return float(array[0]) + float(array[1])/60.0 + float(array[2])/3600.0


def endNight(site, telescope, imager):

    # park the scope
    try: logger.info("Parking/Disconnection from Telescope")
    except: pass
    telescope.shutDown()
    
    # Backup the data
    try: logger.info("Compressing data")
    except: pass
    #compressData(imager.dataPath) # does nothing

    # Turn off the camera cooler, disconnect
    try: logger.info("Disconnecting imager")
    except: pass
    try: imager.shutDown()
    except: pass

    #TODO: Back up the data
    try: logger.info("Backing Up data")
    except: pass
    #backupData(imager.dataPath) # does nothing


    # copy schedule to data directory
    try: logger.info("Copying schedule file from ./schedule/" + site.night + " to " + imager.dataPath)
    except: pass
    shutil.copyfile("C:\camal\schedule" + '//' + site.night, imager.dataPath + site.night)

    # copy logs to data directory
    logs = glob.glob("./logs/" + site.night + "/*.log")
    for log in logs:
        try: logger.info("Copying log file " + log + " to " + imager.dataPath)
        except: pass
        shutil.copyfile(log, imager.dataPath + os.path.basename(log))

    #### create an observing report ####

    # compose the observing report
    body = "Dear humans,\n\n" + \
           "While you were sleeping, I observed:\n\n" + \
           "Love,\n" + \
           "CAMAL"

    # email observing report
    mail.send('camal is done observing',body)


def prepNight(email=True):

    # reset the night at 10 am local
    today = datetime.datetime.utcnow()
    if datetime.datetime.now().hour >= 10 and datetime.datetime.now().hour <= 16:
        today = today + datetime.timedelta(days=1)
    night = 'n' + today.strftime('%Y%m%d')

    # make sure the log path exists
    logpath = 'logs/' + night + '/'
    if not os.path.exists(logpath):
        os.makedirs(logpath)
        
    hostname = socket.gethostname()

    site = camalsite.site('Mount_Hopkins',night,configfile='camal_class_files/site.ini')
    site.night = night

    site.startNightTime = datetime.datetime(today.year, today.month, today.day, 17) - datetime.timedelta(days=1)

    telescope = camaltelescope.cMyT(night)
    imager = camalimager.cSBIGtheskyX(night)
    #imager = camalimager.cSBIGmaxim(night)

    imager.dataPath = "C:/camal/data/" + night + "/"
    imager.finaldataPath = "C:/camal/final/" + night + "/"

    if not os.path.exists(imager.dataPath):
        os.makedirs(imager.dataPath)
    # email notice
    if email: mail.send('CAMAL Starting observing','Love,\nCAMAL')

    return site, telescope, imager


def backupData(dataPath,night):
    """
    backup raw images to ashleys google drive
    """
    # Collect data files to upload
    files = glob.glob(dataPath + "/*.fits")

    # Start pydrive instance
    gauth = GoogleAuth()
    auth_url = gauth.GetAuthUrl() # Create authentication url user needs to visit
    gauth.LocalWebserverAuth()  #will this require me to click a button every time?????
    #http://pythonhosted.org/PyDrive/oauth.html#customizing-authentication-with-settings-yaml
    drive = GoogleDrive(gauth)

    # Check if folder in CAMAL_data for night
    camalid   = '0BxAnuNxRgJ3vcWN5TVV4Q0J3NUE' # id of CAMAL data folder
    file_list = drive.ListFile({'q': "'%s' in parents and trashed=false" % camalid}).GetList()
    folders   = {}
    for f in file_list:
        if f['mimeType']=='application/vnd.google-apps.folder': # if folder
            folders[str(f['title'])] = str(f['id'])
    logging.info('Create Folder in Google Drive for Back Up')
    if night in folders.keys():
        nightid = folders[night]  # store night id
    else:
        # Create folder
        nightfolder = drive.CreateFile({'title': night, 
                                        "parents":  [{"id": camalid}], 
                                        "mimeType": "application/vnd.google-apps.folder"})
        nightfolder.Upload()

    # Upload Files to night's folder
    for filepath in files:
        fname = filepath.split('\\')[-1]
        f     = drive.CreateFile({'title':fname,"parents": [{"kind": "drive#fileLink", "id": nightid}]})
        f.SetContentFile( filepath )
        f.Upload()
 
    logging.info('Successfully Uploaded Files from %s' %night)
 


def doBias(site, telescope, imager, num=20):
    doDark(site, telescope, imager,exptime=0,num=num)

def doDark(site, telescope, imager, exptime, num):
    DARK = 0
    if exptime == 0:
        objectName = 'BIAS'
        for x in range(num):    
            logger.info('Taking ' + objectName + ' ' + str(x+1) + ' of ' + str(num) + ' (exptime = ' + str(exptime) + ')')
            takeImage(site, telescope, imager, exptime,objectName)
    else:
        objectName = 'DARK'
        # Take num Dark frames and loop over more than one exptime
        for time in exptime:
            for x in range(num):    
                logger.info('Taking ' + objectName + ' ' + str(x+1) + ' of ' + str(num) + ' (exptime = ' + str(time) + ')')
                takeImage(site, aqawan, telescope, imager, time,objectName)


def isSupersaturated(filename,satLevel=60000,size=1000):
    """
    Inputs: Filename   - name of file to check
            satLevel   - max pixel value over which returns True
                            default = 60000
            size       - max size of image to check (square)
                            default = 1000
    """
    image = pyfits.getdata(filename,0)
    # Take 1000x1000 field to avoid edges
    nx = len(image)
    ny = len(image[1])
    
    x1 = max(nx/2.0 - size/2.0,0)
    x2 = min(nx/2.0 + size/2.0,nx-1)
    y1 = max(ny/2.0 - size/2.0,0)
    y2 = min(ny/2.0 + size/2.0,ny-1)

    if np.max(image) > satLevel:
        return True
    
    return False


def takeImage(site, telescope, imager, exptime, objname,opt='full',frame=(0,0,0,0),filterlam=0):
    """
    Filterlam: options are 780, 823, 860, and 0 (no filter)
    opt: 'full' or 'sub'
    	if 'sub' then must provide frame=(l,b,r,t) tuple with frame boundaries
    """
    if imager.program == 'maxim':
    	exptypes = {
        	'DARK' : 0,
        	'BIAS' : 0,
        	'FLAT' : 1,
        	}

    if imager.program == 'theskyX':
    	exptypes = {
        	'DARK' : 3,
        	'BIAS' : 2,
        	'FLAT' : 4,
        	}
   	
    if objname in exptypes.keys():
        exptype = exptypes[objname]
    else: exptype = 1 # science exposure

    if str(filterlam) not in imager.filters.keys():
        logger.error("Requested filter (" + filterlam + ") not present")
        return
   
    # Take image
    logger.info("Taking " + str(exptime) + " second image")
    t0 = datetime.datetime.utcnow()
    l,b,r,t = frame
    imager.setFrame(opt,l=l,b=b,r=r,t=t)
    imager.Expose(exptime, exptype, str(filterlam))

    if exptype == 1:
        # Get status info for headers while exposing/reading out 
        telra  = telescope.getRA()
        teldec = telescope.getDEC()

        moonpos   = site.moonpos()
        moonra    = moonpos[0]
        moondec   = moonpos[1]
        moonsep   = ephem.separation((telra*math.pi/180.0,teldec*math.pi/180.0),moonpos)*180.0/math.pi
        moonphase = site.moonphase()

        # todo: use moon info or put in header at least

    # Save the image
    filename = imager.dataPath + "/" + getIndex(imager.dataPath) + objname + "." + str(filterlam) + ".fits"
    logger.info('Saving image: ' + filename)
    imager.saveImage(filename)


def getIndex(dirname):
    files = glob.glob(dirname + "/*.fits*")

    return str(len(files)+1).zfill(4)

    if len(files) == 0:
        return '0001'

    lastnum = files[-1][0:4]
    index = str(int(lastnum) + 1).zfill(4)
    return index

def loadSchedule(night):
    """
    Load Schedule assigned to that night
    returns schedule with target names
            and schedic with target info
            RA,DEC in degrees right now...do this better later?
            like read big file from online instead
            -Use imag to estimate exp times in each filter
    """
    f = open('C:/camal/schedule' + '\\' + night)
    lines = f.readlines()
    schedic  =  {}
    schedule =  []
    #logger.info('Reading Schedule File: ' + file)
    for line in lines:
        if not line.startswith('#'):
            stuff =  line.split()
            schedule.append(line.split()[0])
            # Fill in Dictionary
            schedic[stuff[0]]              = {}
            schedic[stuff[0]]['starttime'] = formatTime(stuff[1],night)
            schedic[stuff[0]]['endtime']   = formatTime(stuff[2],night)
            schedic[stuff[0]]['name']      = stuff[3]
            schedic[stuff[0]]['RA']        = stuff[4]
            schedic[stuff[0]]['DEC']       = stuff[5]
            schedic[stuff[0]]['imag']      = stuff[6]

    return schedule,schedic

def formatTime(string,night):
    strarr = string.split(':')
    hour   = int(strarr[0])
    minute = int(strarr[1])
    second = int(float(strarr[2]))
    day=int(night[7:])
    year=int(night[1:5])
    month=int(night[5:7])
    if hour < 14:
        return datetime.datetime(second=second, minute=minute, hour=hour,year=year,month=month,day=day)
    else:
        return datetime.datetime(second=second, minute=minute, hour=hour,year=year,month=month,day=day) - datetime.timedelta(1)

def CalcExpTime(imag,filt,area=0.07,snr=1000):
    """ 
    Calculate the exposure time of an astronomical object
    filter options: 780,823,860 (assuming i band flux)
    area meters squared for camal (rough, check later)
    Source: https://www.astro.umd.edu/~ssm/ASTR620/mags.html
    """
    if filt=='780':
        cent  = 0.78 # microns
        dll   = 0.00385 # dlambda/lambda
        flux0 = 2550 # Jy at top of atm
        ccdqe = 0.33
    if filt=='823':
        cent  = 0.823    
        dll   = 0.00389  
        flux0 = 2550
        ccdqe = 0.26
    if filt=='860':
        cent  = 0.86    
        dll   = 0.00221    
        flux0 = 2550
        ccdqe = 0.19
    # Calculate flux at top of atm
    f       = flux0 * 10**(-0.4*imag) * 1.51e7 * dll #photons/sec/m^2      
    # Calc again with eff
    eff     = 0.9 * ccdqe * 0.9 # atm, ccd, telescope throughput
    f_eff   = f * eff * area
    # Calculate exp time in seconds
    exptime = (snr**2/f_eff)
    return max(0.01,round(exptime,2))  # round to nearest milisecond

def doSkyFlat(site, telescope, imager, exptimes, nflat=30):
    # Calc exposure times
    exptime780, exptime823, exptime860 = exptimes

    logger.info('Looping through filters now for target')
    takeImage(site, telescope, imager, exptime823, 'FLAT',filterlam=823)
    takeImage(site, telescope, imager, exptime860, 'FLAT',filterlam=860)
    takeImage(site, telescope, imager, exptime780, 'FLAT',filterlam=780)


def doScience(site, telescope, imager, targetinfo):
    if datetime.datetime.utcnow() > targetinfo['endtime']:
        logger.info("Target " + targetinfo['name'] + " past its endtime (" + str(targetinfo['endtime']) + "); skipping")
        return

    # if before start time, wait
    if datetime.datetime.utcnow() < targetinfo['starttime']:
        waittime = (targetinfo['starttime']-datetime.datetime.utcnow()).total_seconds()
        logger.info("Target " + targetinfo['name'] + " is before its starttime (" + str(targetinfo['starttime']) + "); waiting " + str(waittime) + " seconds")
        time.sleep(waittime)

    # Slew to Target
    telescope.gotoObject(targetinfo['name'])
    telescope.tel.Tracking = True

    # Calc exposure times
    imag = float(targetinfo['imag'])
    exptime780 = CalcExpTime(imag,'780')
    exptime823 = CalcExpTime(imag,'823')
    exptime860 = CalcExpTime(imag,'860')

    # Find Star center and pick subframe
#    takeImage(site, telescope, imager, 0.1, targetinfo['name'],opt='sub',frame=(l,b,r,t),filterlam=)
    #^^ do this if centering still sucks

    # Define subframe
    centx, centy = 1663, 1252 # use config file
    boxwidth = 400
    l,b,r,t = centx - boxwidth, centy+boxwidth, centx+boxwidth, centy-boxwidth # top and bottom switched

    # Set time tracking
    timetracker = time.time()

    # Start Loop of Observations
    logger.info('Looping through filters now for target')
    while datetime.datetime.utcnow() < targetinfo['endtime']:
        if isDomeOpen():
            takeImage(site, telescope, imager, exptime780, targetinfo['name'],opt='sub',frame=(l,b,r,t),filterlam=780)
            takeImage(site, telescope, imager, exptime823, targetinfo['name'],opt='sub',frame=(l,b,r,t),filterlam=823)
            takeImage(site, telescope, imager, exptime860, targetinfo['name'],opt='sub',frame=(l,b,r,t),filterlam=860)
            time.sleep(1)
        else:
            print 'Dome Is not Open'
        if datetime.datetime.utcnow() > targetinfo['endtime']: return
        if time.time() - timetracker > 600:
            telescope.gotoObject(targetinfo['name']) # slew again to recenter
            downloadSkyJPG(imager.dataPath)
            timetracker = time.time()


def sourceFinder(science):
    ' Input: loaded data, Output: x,y center of aperture '
    # Convolve image with gaussian kernel
    kernel = np.outer(signal.gaussian(3,3), signal.gaussian(3,3))
    blurred = signal.fftconvolve(science, kernel, mode='same')
   
    # Take the normalized STD along x,y axes
    xstd = np.std(blurred,axis=0)
    ystd = np.std(blurred,axis=1)
    xstdn = (xstd - np.median(xstd))/max(xstd)
    ystdn = (ystd - np.median(ystd))/max(ystd)
    
    # Determine center by maximum. Eventually add check that there's only one source!
    try: x,y = np.where(xstdn == max(xstdn))[0][0], np.where(ystdn == max(ystdn))[0][0]
    except IndexError:
        x,y = 0,0
        
    return x,y


def downloadSkyJPG(path):
    strtime = time.strftime("%Y-%m-%dT%H.%M.%S")
    urllib.urlretrieve("http://linmax.sao.arizona.edu/weather/sky/last_sky.jpg", "%s//0000last_sky%s.jpg" %(path, strtime))

def isDomeOpen():
    os.system('scp -P 22222 minerva@minerva.sao.arizona.edu:/home/minerva/minerva-control/minerva_library/aqawan1.stat .')
    f = open('aqawan1.stat','r')
    line = f.readlines()
    f.close()
    if len(line) != 0:
        if line[0].split()[-1] == 'False':
            return False
        else:
            return True
    else:
        return True # if dome file is empty, keep observing

def updateReduction(night):
    """
    Write night's name to reduction file with 0 to indicate
    that it needs to be reduced
    Reduction code in analysis_code will edit the 0 to 1 when it's 
    successfully reduced or -1 if it failed
    """
    reductionfile = "C:/camal/reduction_code/" + 'reductionlist.txt'
    if os.path.isfile(reductionfile):
        with open(reductionfile,'a') as f:
            f.write("%s  0\n" %night)
        f.close()
    else:
        f = open(reductionfile,'w')
        f.write("%s  0\n" %night)
        f.close()


# --------------------------------------------------------#



if __name__ == '__main__':

    # this file is created at the beginning and deleted at the end
    # to enable automatic recovery
    with open('running.txt','w') as f:
        f.write(str(datetime.datetime.utcnow()))

    site, telescope, imager = prepNight()

    print 'night = ', site.night

    # setting up main logger
    fmt = "%(asctime)s [%(filename)s:%(lineno)s - %(funcName)s()] %(levelname)s: %(message)s"
    datefmt = "%Y-%m-%dT%H:%M:%S"

    logger = logging.getLogger('main')
    formatter = logging.Formatter(fmt,datefmt=datefmt)
    formatter.converter = time.gmtime
        
    fileHandler = logging.FileHandler('logs/' + site.night + '/main.log', mode='a')
    fileHandler.setFormatter(formatter)

    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)
        
    logger.setLevel(logging.DEBUG)
    logger.addHandler(fileHandler)
    logger.addHandler(console)

    # Wait until sunset to start biases (maybe in future start this earlier to have time for flats?)  
    timeUntilSunset = (site.sunset() - datetime.datetime.utcnow()).total_seconds()
    if timeUntilSunset > 0:
        logger.info('Waiting for sunset (' + str(timeUntilSunset) + 'seconds)')
        time.sleep(timeUntilSunset)

    # wait until it's darker to take biases/darks
    readtime = 1.0 #sec
    nbias = 30

    # Connect to imager
    imager.connect(cooler=True)

    # Take Biases
    logger.info('Taking %s Biases' %nbias)
    doBias(site, telescope, imager, num=nbias)

    # Wait until nautical twilight ends 
    timeUntilTwilEnd = (site.NautTwilEnd() - datetime.datetime.utcnow()).total_seconds()
    if timeUntilTwilEnd > 0:
        logger.info('Waiting for nautical twilight to end (' + str(timeUntilTwilEnd) + 'seconds)')
        time.sleep(timeUntilTwilEnd)

    # Load the night's schedule    
    try:
        logger.info('Making Schedule')
        setsched = camalscheduler.scheduler(site)
        setsched.run()
        setsched.save_schedule()
        # Load Schedule
        schedule, schedic = loadSchedule(site.night)
    except IOError:
        telescope.park()
        logger.info('Schedule did not work. Breaking.')
        mail.send('Forgot to Setup Schedule!!!', 'Hi there, please make a schedule'
                'for the night %s' %site.night)
 

    # Check if Dome is open. If not, wait 10min
    while not isDomeOpen():
        logger.info('Waiting For the Dome to Open at Beginning of Night')
        if (site.NautTwilBegin() - datetime.datetime.utcnow()).total_seconds() < 600:
            break
        time.sleep(60)


    if (site.NautTwilBegin() - datetime.datetime.utcnow()).total_seconds() < 600:
        logger.info('Past Nautical Twilight. Should close up')
    else:
        telescope.initialize()

    # read the target list
    for target in schedule:
        print target
        # Load Target Schedule Information
        targetinfo = schedic[target]

        # compute the rise/set times of the target
        site.obs.horizon = '20.0'
        body = ephem.FixedBody()
        body._ra  = str(ten(targetinfo['RA']))
        body._dec = str(ten(targetinfo['DEC']))
        body._epoch = '2000.0'
        body.compute()
        try:
            risetime = site.obs.next_rising(body,start=site.NautTwilEnd()).datetime()
        except ephem.AlwaysUpError:
            # if it's always up, don't modify the start time
            risetime = targetinfo['starttime']
        except ephem.NeverUpError:
            # if it's never up, skip the target
            risetime = targetinfo['endtime']
        try:
            settime = site.obs.next_setting(body,start=site.NautTwilEnd()).datetime()
        except ephem.AlwaysUpError:
            # if it's always up, don't modify the end time
            settime = targetinfo['endtime']
        except ephem.NeverUpError:
            # if it's never up, skip the target
            settime = targetinfo['starttime']

        # START OBSERVING
        if (site.NautTwilBegin() - datetime.datetime.utcnow()).total_seconds() > 600:
            logger.info(targetinfo['name']+ ' about to be observed')
            doScience(site, telescope, imager, targetinfo)
            logger.info(targetinfo['name']+ ' done; moving on to next target')
        else:
            pass
                
    # Park
    if telescope.initialized == True:
        logger.info('Parking Telescope (even if it never unparked)')
        telescope.park()

    # Take more biases
    logger.info('End Night')
    endNight(site, telescope, imager)

    # Backup data to ashley's google drive
    try:
        pass#backupData(imager.dataPath,site.night)
    except:
        mail.send('backup failed for night %s' %(site.night), 'from, \n \t CAMAL')

    # Put night's folder into 'To Reduce' list in analysis_code
    updateReduction(site.night)
    logger.info('Updating Reduction File list')

    if os.path.exists('running.txt'):
        os.remove('running.txt')


