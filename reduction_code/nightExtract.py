# To be copied to night of data and run there
# Extract fluxes from each observation
# images must be full (no subframe)
# if binning 1x1 by mistake, this is fixed with rebin1to2

import numpy as np
import os
import matplotlib.pylab as plt
import astropy.io.fits as pyfits
import glob, socket
from astropy.time import Time 
import datetime
import ephem
import urllib

from camal_reduction_classes import fluxExtract,camal,mail,calibrate
from camal_reduction_classes import reduce as camalreduce

def prepNight(email=True):
  night = editReductionList(mode='r')
  if night == None:
    print 'Will break if there are no nights to reduce'

  print 'night = ', night
  tonight = camalreduce.tonight('Mount_Hopkins',night,configfile='camal_reduction_classes/reduce.ini') 
  # Load target names from schedule but maybe should put this in reduction list
  # in case things don't get observed..

  schedule,schedic      = loadSchedule(night)
  tonight.objects       = schedule
  tonight.schedic       = schedic
  tonight.dataPath      = 'C:/camal/data/' + '\\' + night 
  tonight.finaldataPath = 'C:/camal/final/' + '\\' + night

  # make sure the log path exists
  logpath = 'logs/' + tonight.night + '/'
  if not os.path.exists(logpath):
    os.makedirs(logpath)
        
  hostname = socket.gethostname()
  if not os.path.isdir(tonight.finaldataPath):
    os.makedirs(tonight.finaldataPath)

  # email notice
  if email: mail.send('CAMAL Started reducing night %s' %tonight.night,'Love,\nCAMAL')

  return tonight

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
        return datetime.datetime(second=second, minute=minute, hour=hour,year=year,month=month,day=(day-1))

def get_alt(tonight,objname,jd):
    targetinfo  = tonight.schedic[objname]
    t=Time(jd, format='jd')                                                                  
    obs = ephem.Observer()
    obs.lat = ephem.degrees(str(tonight.latitude)) # N                                                        
    obs.lon = ephem.degrees(str(tonight.longitude)) # E                                                       
    obs.elevation = 2000 # meters                                                                
    obs.date  = t.iso.replace('-','/')
    obs.horizon = '20.0'
    body = ephem.FixedBody()
    body._ra  = str(ten(targetinfo['RA'])) 
    body._dec = str(ten(targetinfo['RA'])) 
    body._epoch = str(tonight.elevation)
    body.compute(obs)
    return body.alt

def ten(string):
    array = string.split(':')
    if "-" in array[0]:
        return float(array[0]) - float(array[1])/60.0 - float(array[2])/3600.0
    return float(array[0]) + float(array[1])/60.0 + float(array[2])/3600.0

def editReductionList(mode='r',night=None,error=False):
  """
  Read or write reductionlist.txt file
  if mode == 'r' then return the night to reduce
  if mode == 'w' then provide night and change 0 to 1
  if error==True in write mode, then will change 0 to -1
  """
  reductionfile = "C:/camal/reduction_code/" + 'reductionlist.txt'
  if mode == 'r':
    f = open(reductionfile,'r')
    lines = f.readlines()
    f.close()
    for line in lines:
      if line.split()[1] == '0':
        return line.split()[0]

  if mode == 'w':
    f = open(reductionfile,'r')
    lines = f.readlines()
    f.close()
    fnew = open(reductionfile,'w')
    for line in lines:
      if line.startswith(night):
        if error:
          fnew.write("%s  -1\n" %night)
        if not error:
          fnew.write("%s  1\n" %night)
      else:
        fnew.write(line)
    return


def updateAnalysis(night):
    """
    Write night's name to analysis file with 0 to indicate
    that it needs to be analyzed
    analysis code in analysis_code will edit the 0 to 1 when it's 
    successfully reduced or -1 if it failed
    """
    analysisfile = "C:/camal/analysis_code/" + 'analysislist.txt'
    if os.path.isfile(analysisfile):
        with open(analysisfile,'a') as f:
            f.write("%s  0\n" %night)
        f.close()
    else:
        f = open(analysisfile,'w')
        f.write("%s  0\n" %night)
        f.close()

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
            schedule.append(line.split()[3])
            # Fill in Dictionary
            schedic[stuff[3]]              = {}
            schedic[stuff[3]]['starttime'] = formatTime(stuff[1],night)
            schedic[stuff[3]]['endtime']   = formatTime(stuff[2],night)
            schedic[stuff[3]]['name']      = stuff[3]
            schedic[stuff[3]]['RA']        = stuff[4]
            schedic[stuff[3]]['DEC']       = stuff[5]
            schedic[stuff[3]]['imag']      = stuff[6]

    return schedule,schedic

def getWeather(night):
    """
    Downloads night and night before's weather report from Mt Hopkins Ridge data
    center to final data folder for the night
    """
    finalpath = 'C:\camal\\final' + '//' + night
    day   =int(night[7:])
    year  =int(night[1:5])
    month =int(night[5:7])
    date1 = datetime.datetime(year=year,month=month,day=day)
    date2 = datetime.datetime(year=year,month=month,day=day)  - datetime.timedelta(1)
    weather1 = 'http://www.sao.arizona.edu/weather/wview/img/Archive/ARC-%s.txt' %date1.strftime('%Y-%m-%d')
    weather2 = 'http://www.sao.arizona.edu/weather/wview/img/Archive/ARC-%s.txt' %date2.strftime('%Y-%m-%d')
    urllib.urlretrieve(weather1, finalpath + "//ridgeweather_%s.txt" %date1.strftime('%Y-%m-%d'))
    urllib.urlretrieve(weather2, finalpath + "//ridgeweather_%s.txt" %date2.strftime('%Y-%m-%d'))

def doReduction(fname,objname,bias=0,flat=0,dark=0,plotap=False):
  science, hdr = fluxExtract.loadData(fname)
  if hdr['XBINNING'] == 1:
    science = camal.rebin1to2(science,mleft = mleft,mright=mright)

  # Bad Pixels - change this so it medians neigboring pixels
  badpix_x, badpix_y = tonight.badpix_x, tonight.badpix_y #update this
  science[badpix_y,badpix_x] = science[badpix_y+2,badpix_x+2]
  flux, flags,x,y,bkg_mean = fluxExtract.fluxExtract(science,bias,dark,flat,hdr,plotap=plotap)
  altitude = hdr['CENTALT'] 
  #altitude = get_alt(tonight, objname, hdr['JD'])
  obsTime = Time(hdr['DATE-OBS'], format='isot', scale='utc')
  # Return things wanna save
  out = fname.split('\\')[-1], obsTime.jd, flux, altitude, hdr['EXPTIME'], x, y, bkg_mean, hdr['FILTER'],flags
  return out

# 'JD flux alt T_exp xcent ycent bkg_mean filter flags'
def backupDataGoogle(file,tonight):
    """
    backup raw images to ashleys google drive
    """
    # Collect data files to upload
    night = tonight.night
    files = glob.glob(tonight.dataPath + "/*.fits")

    # Start pydrive instance
    gauth = GoogleAuth()
    auth_url = gauth.GetAuthUrl() # Create authentication url user needs to visit
    gauth.LocalWebserverAuth()  #will this require me to click a button every time?????
    #http://pythonhosted.org/PyDrive/oauth.html#customizing-authentication-with-settings-yaml
    drive = GoogleDrive(gauth)

    # Check if folder in CAMAL_data for night
    camalid   = '' # id of CAMAL data folder
    file_list = drive.ListFile({'q': "'%s' in parents and trashed=false" % camalid}).GetList()
    folders   = {}
    for f in file_list:
        if f['mimeType']=='application/vnd.google-apps.folder': # if folder
            folders[str(f['title'])] = str(f['id'])
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

def uploadData(tonight):
  filename = tonight.finaldataPath + 'output.txt'
  if filename: # whi was there an indent here?
    gauth = GoogleAuth()
    auth_url = gauth.GetAuthUrl() # Create authentication url user needs to visit
    gauth.LocalWebserverAuth()  #will this require me to click a button every time?????
    #http://pythonhosted.org/PyDrive/oauth.html#customizing-authentication-with-settings-yaml
    drive = GoogleDrive(gauth)

    # Check if folder in CAMAL_data for night
    camalid   = '' # id of CAMAL data folder
    file_list = drive.ListFile({'q': "'%s' in parents and trashed=false" % camalid}).GetList()
    folders   = {}
    for f in file_list:
        if f['mimeType']=='application/vnd.google-apps.folder': # if folder
            folders[str(f['title'])] = str(f['id'])
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


def cropFits(file):
  """
  Given center of star, crop fits file around it
  """
  pass

def compressNight():
  """
  Create compression file of night's data (after cropping) saved in 
  final data folder...(then use backUptoMinerva function to transfer it to Minerva)
  """
  night = tonight.night
  files = glob.glob(tonight.dataPath + "/*.fits")


def backUptoMinerva(tonight):
  night = tonight.night
  files = glob.glob(tonight.dataPath + "/*")
  # Make Folder on minerva main for night
  out = os.system('ssh -p 22222 minerva@minerva.sao.arizona.edu ls /nas/camal/data/%s' %tonight.night )
  if out == 2:
    os.system('ssh -p 22222 minerva@minerva.sao.arizona.edu mkdir /nas/camal/data/%s' %tonight.night)
  # Move files over
  for i,f in enumerate(files):
    os.system('scp -P 22222 %s minerva@minerva.sao.arizona.edu:/nas/camal/data/%s/' %(f,tonight.night))

def saveFits(file):
  # Variable Names
  fname, time, flux, alt, texp, xcent, ycent, bkg, filt, flags = readOutput(file)

  # Save to Fits
  tbhdu = pyfits.BinTableHDU.from_columns(
       [pyfits.Column(name='filename', format='26A', array=np.array(fname)),
        pyfits.Column(name='JD',format='D',array=np.array(time)),
        pyfits.Column(name='flux', format='E', array=np.array(flux)),
        pyfits.Column(name='alt', format='16A', array=np.array(alt)),
        pyfits.Column(name='expTime',format='E',array=np.array(texp)),
        pyfits.Column(name='Xcenter',format='E',array=np.array(xcent)),
        pyfits.Column(name='Ycenter',format='E',array=np.array(ycent)),
        pyfits.Column(name='BkgMeanCount',format='E',array=np.array(bkg)),
        pyfits.Column(name='filter', format='20A',array=np.array(filt)),
        pyfits.Column(name='flags',format='E',array=np.array(flags)) ])

  if os.path.exists(file.replace('txt','fits')):
    os.system('rm %s' % file.replace('txt','fits'))
    tbhdu.writeto(file.replace('txt','fits'))
  else:
    tbhdu.writeto(file.replace('txt','fits'))

def readOutput(file):
  f = open(file,'r')
  lines = f.readlines()
  f.close()

  fname   = []
  flux    = []
  alt     = []
  filt    = []
  flags   = []
  texp    = []
  xcent   = []
  ycent   = []
  bkg     = []
  juldate = []

  for i in range(1,len(lines)):
    line = lines[i].split()
    fname.append(line[0])
    juldate.append(float(line[1]))
    flux.append(float(line[2]))
    alt.append(float(line[3]))
    texp.append(float(line[4]))
    xcent.append(float(line[5]))
    ycent.append(float(line[6]))
    bkg.append(float(line[7]))
    filt.append(line[8])
    flags.append(float(line[9]))

  return np.array(fname), np.array(juldate), np.array(flux),np.array(alt), np.array(texp), np.array(xcent),np.array(ycent),np.array(bkg),np.array(filt), np.array(flags)


def plotRawFluxes(tonight):
  file = tonight.finaldataPath + '\\' + 'raw_camal_%s.txt' %tonight.night
  params = readOutput(file)
  import matplotlib.pylab as plt
  plt.ion()
  plt.plot(params[1],params[2],'go')

if __name__ == '__main__':
  tonight = prepNight(email=False) # load schedule & reduction list, store things

  # Load Lists of filenames for Data & Calibration
  darkfiles = glob.glob(tonight.dataPath + '\\*DARK.0.fits')
  biasfiles = glob.glob(tonight.dataPath + '\\*BIAS.0.fits')
  flatfiles = glob.glob(tonight.dataPath + '\\*FLAT*.fits')

  # Make or Load Calibration Files, Save in raw folder
  cal  = calibrate.calibrate(tonight.dataPath)
  bias = cal.makeBias(biasfiles)
  flat = cal.makeFlat(flatfiles)
  dark = cal.makeDark(darkfiles)
  # Loop through science files, subtract bias, divide flat, extract flux
  fout         = open(tonight.finaldataPath + '\\' + \
                        'raw_camal_%s.txt' %tonight.night,'w')
  fout.write('# fname JD flux alt T_exp xcent ycent bkg_mean filter flags \n')
  for objname in tonight.objects:
    sciencefiles = glob.glob(tonight.dataPath + '\\*' + objname + '*.fits')
    for sciencefile in sciencefiles:
      out = doReduction(sciencefile,objname,bias=bias,flat=flat,dark=dark,plotap=False)
      str_out = '%s'
      for var in range(len(out)-1):
        str_out += '\t%s'
      fout.write(str_out %out + '\n')
      print sciencefile
  
  fout.close()

  saveFits(tonight.finaldataPath + '\\' + 'raw_camal_%s.txt' %tonight.night)

  # Update reduction list saying night is reduced
  editReductionList(mode='w',night=tonight.night,error=False)

  # Edit Analysis List saying this night needs PWV extraction
  updateAnalysis(tonight.night)

  # Download Weather from ridge, save to final folder
  getWeather(tonight.night)
  # crop all files & compress then backup onto minerva machine
  # copy log file (which need to add to) to final night data

  # Backup data
  backUptoMinerva(tonight)

  # Make Plots

