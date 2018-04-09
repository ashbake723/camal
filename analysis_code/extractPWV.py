#extract PWV for night of photometry

import numpy as np
import os
import matplotlib.pylab as plt
from astropy.stats import median_absolute_deviation
import astropy.io.fits as pyfits
import glob, socket
from astropy.time import Time 
from scipy import interpolate
import datetime
import ephem
import urllib
import cPickle as pickle
import scipy.optimize as op
import emcee, corner

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

from camal_analysis_classes import objectextract as obj_ext


def prepNight(email=False):
  night = editAnalysisList(mode='r') # import night in same way do reduction?

  if night == None:
    print 'There are no nights to analyze!'
    return None

  print 'night = ' + night

  oe = obj_ext.tonight('Mount_Hopkins',night,configfile='camal_analysis_classes/analyze.ini') 
  oe.night = night
  oe.dataPath = 'C:/camal/data/' + '\\' + night
  oe.finaldataPath = 'C:/camal/final/' + '\\' + night
  oe.schedulePath = 'C:/camal/schedule/' + '\\' + night
  oe.schedule  , oe.schedic = loadSchedule(oe.night)

  # Load Data
  f = pyfits.open(oe.finaldataPath + \
        '\\' + 'raw_camal_%s.fits' %oe.night)
  out = f[1].data
  f.close()

  oe.jdraw = out['JD']

  # make sure the log path exists
  logpath = 'logs/' + oe.night + '/'
  if not os.path.exists(logpath):
    os.makedirs(logpath)
        
  hostname = socket.gethostname()

  # email notice
  if email: mail.send('CAMAL Started extracting PWVs for night %s' %oe.night,'Love,\nCAMAL')

  return oe


def loadData(oe):
  fitsname = oe.finaldataPath + \
        '\\' + 'raw_camal_%s.fits' %oe.night
  oe.f780, oe.f823, oe.f860, \
       oe.alt, oe.jd, oe.objnames, oe.fnames = cleanData(oe.finaldataPath + \
       	'\\' + 'raw_camal_%s.fits' %oe.night)
  oe.X = 1/np.cos((np.pi/180.0) * (90.0 - oe.alt))
  oe.c1 = oe.f823/oe.f780
  oe.c2 = oe.f823/oe.f860


def editAnalysisList(mode='r',night=None,error=False):
  """
  Read or write reductionlist.txt file
  if mode == 'r' then return the night to reduce
  if mode == 'w' then provide night and change 0 to 1
  if error==True in write mode, then will change 0 to -1
  """
  analysisfile = "C:/camal/analysis_code/" + 'analysislist.txt'
  if mode == 'r':
    f = open(analysisfile,'r')
    lines = f.readlines()
    f.close()
    for line in lines:
      if line.split()[1] == '0':
        return line.split()[0]

  if mode == 'w':
    f = open(analysisfile,'r')
    lines = f.readlines()
    f.close()
    fnew = open(analysisfile,'w')
    for line in lines:
      if line.startswith(night):
        if error:
          fnew.write("%s  -1\n" %night)
        if not error:
          fnew.write("%s  1\n" %night)
      else:
        fnew.write(line)
    return

def formatTime(string,night):
  """
  Brute forcing the conversion of string time to datetime object
  """
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

def getScalar(oe,mrat):
  """
  - Take f860 to f780 ratio and correct f860 based on this
  - Use TRES data paper results to get f823/f780 ratio
  """
  # Define A2
  A2 = 0.979    # Scalar, comes from fitting TRES

  # Fit f860 to f780 (c1/c2)
  crat = oe.c1/oe.c2
  #  mrat = oe.m1[:,100]/oe.m2[:,100] # this is same for all PWV so take any 

  # Find where object switches
  # Create array of index of switch of targets
  iswitch1 = np.where((oe.objnames[1:] == oe.objnames[:-1]) == False)[0]
  if len(iswitch1) == 0:
    iswitch = np.array([0,len(oe.f780)])
  else:
    iswitch = np.concatenate(([0],iswitch1,[len(oe.f780)]))

  A1 = np.zeros((len(iswitch)-1,250,250))

  # Do for each star
  for ii in range(len(iswitch)-1):
    A1[ii,:,:] = np.median(crat[iswitch[ii]:iswitch[ii+1]])/np.median(mrat[iswitch[ii]:iswitch[ii+1]])

  return A1, A2

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
            schedic[stuff[3]]['temp']      = int(stuff[7])

    return schedule,schedic

def getKittPeak(oe):
  pass

def uploadData(oe):
  filename = oe.finaldataPath + '\\PWVout_%s.fits' %oe.night

  def ListFolder(parent):
    filelist=[]
    file_list = drive.ListFile({'q': "'%s' in parents and trashed=false" % parent}).GetList()
    for f in file_list:
      if f['mimeType']=='application/vnd.google-apps.folder': # if folder
        filelist.append({"id":f['id'],"title":f['title'],"list":ListFolder(f['id'])})
      else:
        filelist.append({"title":f['title'],"title1":f['alternateLink']})
    return filelist

  if filename: # whi was there an indent here?
    gauth = GoogleAuth()
    auth_url = gauth.GetAuthUrl() # Create authentication url user needs to visit
    gauth.LocalWebserverAuth()  #will this require me to click a button every time?????
    #http://pythonhosted.org/PyDrive/oauth.html#customizing-authentication-with-settings-yaml
    drive = GoogleDrive(gauth)

   # Check if folder in CAMAL_data for night
    out_root = ListFolder('root')
    for ifold in range(len(out_root)):
      if out_root[ifold]['title'] == u'Data':
        camalid  = out_root[ifold]['id']  
    #camalid   = '0B18tnyqgrwlpbFV0aXA4ckxXUlE' # id of CAMAL data folder
    #file_list = drive.ListFile({'q': "'%s' in parents and trashed=false" % camalid}).GetList()
    file_list = ListFolder(camalid)
    folders   = {}
    for f in file_list:
      folders[str(f['title'])] = str(f['id'])
    if oe.night in folders.keys():
        nightid = folders[oe.night]  # store night id
    else:
      # Create folder
      nightfolder = drive.CreateFile({'title': oe.night, 
                                      "parents":  [{"id": camalid}], 
                                      "mimeType": "application/vnd.google-apps.folder"})
      nightfolder.Upload()
      file_list = ListFolder(camalid)
      for f in file_list:
        folders[str(f['title'])] = str(f['id'])
      nightid = folders[oe.night]  # store night id

    files = glob.glob(oe.finaldataPath + '\\*')
    # Upload Files to night's folder
    for filepath in files:
        fname = filepath.split('\\')[-1]
        f     = drive.CreateFile({'title':fname,"parents": [{"kind": "drive#fileLink", "id": nightid}]})
        f.SetContentFile( filepath )
        f.Upload()

def readOutput(filename):
  f = open(filename,'r')
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
    alt.append(line[3])
    texp.append(float(line[4]))
    xcent.append(float(line[5]))
    ycent.append(float(line[6]))
    bkg.append(float(line[7]))
    filt.append(line[8])
    flags.append(float(line[9]))

  return np.array(fname), np.array(juldate), np.array(flux),np.array(alt), np.array(texp), np.array(xcent),np.array(ycent),np.array(bkg),np.array(filt), np.array(flags)

def cleanData(fitsname):
  """
  Load final data and output time and colors and Xs

  this is pretty bad - make a mask instead of 0s and 1s
  """
  #fitsname = oe.finaldataPath + \
  #      '\\' + 'raw_camal_%s.fits' %oe.night

  # Load Fits File (whole path)
  f = pyfits.open(fitsname)
  out = f[1].data
  f.close()

  # get alt, objnames
  alt = np.array(out['alt']).astype('float')
  objnames = []
  for fname in out['filename']:
    objnames.append(fname.split('.')[0][4:])
  objnames = np.array(objnames)

  # Pull out each color
  i780 = np.where(out['filter']=='Lunar')[0]
  i823 = np.where(out['filter']=='Blue')[0]
  i860 = np.where(out['filter']=='Clear')[0]

  # Add better check that there are no time delays
  if (len(i780) != len(i823)) or (len(i823) != len(i860)):
    raise TypeError('Number of 823 images not the' 
           'number of 860 of 780 images') # code solution when this happens

  # Proceed for perfect case
  f780 = (out['flux']/out['expTime'])[i780]
  f823 = (out['flux']/out['expTime'])[i823] 
  f860 = (out['flux']/out['expTime'])[i860]

  # FIRST THING ---------------------------
  # ------------------

  # Create array of index of switch of targets
  iswitch1 = np.where((objnames[i780][1:] == objnames[i780][:-1]) == False)[0]
  if len(iswitch1) == 0:
    iswitch = np.array([0,len(f780)])
  else:
    iswitch = np.concatenate(([0],iswitch1,[len(f780)]))

  #Take bkgmean over flux to be goodness indicator
  sf780 = (out['BkgMeanCount']/out['flux'])[i780]
  sf823 = (out['BkgMeanCount']/out['flux'])[i823]
  sf860 = (out['BkgMeanCount']/out['flux'])[i860]
  gdness = sf780 + sf823 + sf860

  # FInd median of goodness indicator for each star
  meds = np.zeros(len(f780))
  stds = np.zeros(len(f780))
  for i in range(len(iswitch)-1):
    meds[iswitch[i]:iswitch[i+1]] = np.median(gdness[iswitch[i]:iswitch[i+1]])
    stds[iswitch[i]:iswitch[i+1]] = 3*median_absolute_deviation(gdness[iswitch[i]:iswitch[i+1]])

  # Pick bad points
  ibad1 = np.where(np.abs(gdness - meds) > stds)[0]
 
  # SECOND THING ---------------------------
  # ------------------
  
  objnames2 = np.delete(objnames[i780],ibad1)

  # Create array of index of switch of targets
  iswitch1 = np.where((objnames2[1:] == objnames2[:-1]) == False)[0]
  if len(iswitch1) == 0:
    iswitch2 = np.array([0,len(objnames2)])
  else:
    iswitch2 = np.concatenate(([0],iswitch1,[len(objnames2)]))

  gdness2 = np.delete(f780+f823+f860,ibad1)

  # Loop through each star
  # FInd median of goodness indicator for each star
  meds = np.zeros(len(gdness2))
  stds = np.zeros(len(gdness2))
  for i in range(len(iswitch)-1):
    meds[iswitch2[i]:iswitch2[i+1]] = np.median(gdness2[iswitch2[i]:iswitch2[i+1]])
    stds[iswitch2[i]:iswitch2[i+1]] = 3*median_absolute_deviation(gdness2[iswitch2[i]:iswitch2[i+1]])

  ibad2 = np.where(np.abs(gdness2 - meds) > stds)[0]

  # THIRD THING ---------------------------
  # ------------------
  
  objnames3 = np.delete(objnames2,ibad2)

  # Create array of index of switch of targets
  iswitch1 = np.where((objnames3[1:] == objnames3[:-1]) == False)[0]
  if len(iswitch1) == 0:
    iswitch3 = np.array([0,len(objnames3)])
  else:
    iswitch3 = np.concatenate(([0],iswitch1,[len(objnames3)]))

  # REMOVE WHERE INDIVIDUAL FLUX DROPS BELOW 50%
  gdness3 = np.delete(np.delete(f860,ibad1),ibad2)
  max860 = np.zeros(len(gdness3))
  for i in range(len(iswitch3)-1):
    max860[iswitch3[i]:iswitch3[i+1]] = np.max(gdness3[iswitch3[i] + 2 :iswitch3[i+1] - 2])

  #ibad3   = np.where((gdness3 < (0.5*max860)) or (gdness < 50000))[0]
  ibad3 = np.arange(len(gdness3))[np.any((gdness3 < (0.5*max860),gdness3 < 10000.0),axis=0)]

  # Return wanted items minus bad indices
  # wallahi theres a better to do this ..
  return np.delete(np.delete(np.delete(f780,ibad1),ibad2),ibad3), \
  np.delete(np.delete(np.delete(f823,ibad1),ibad2),ibad3), \
  np.delete(np.delete(np.delete(f860,ibad1),ibad2),ibad3),\
  np.delete(np.delete(np.delete(alt[i780],ibad1),ibad2),ibad3), \
  np.delete(np.delete(np.delete(out['JD'][i780],ibad1),ibad2),ibad3), \
  np.delete(np.delete(np.delete(objnames[i780],ibad1),ibad2),ibad3),\
  np.delete(np.delete(np.delete(out['filename'][i780],ibad1),ibad2),ibad3)

def binData(oe):
  """
  Bin data inputs by 15 min? Reject outliers in bins & get error est.

  outputs: data, err, gtime binned
  """
  dt    = 5 # bin interval in minutes
  nbins = int(np.round( (oe.jd[-1] - oe.jd[0]) / (dt/60.0/24.0) ))
  
  # Binning part
  jdbin, c1bin, c1err = binfunc(oe.jd,oe.c1,nbins)
  jdbin, c2bin, c2err = binfunc(oe.jd,oe.c2,nbins)

  return c1bin,c2bin,c1err,c2err,jdbin

def binfunc(x,y,nbins):
    n, bins = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    bins = (bins[1:] + bins[:-1])/2
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    return bins,mean,std


def loadColorGrids(oe):
    # Load Spectral Grids
    specgrid780 = []
    specgrid823 = []
    specgrid860 = []

    # Extract Temperatures
    Temps = []
    for obj in oe.schedic:
      if obj in oe.objnames:
        Temps.append(np.round(oe.schedic[obj]['temp'],-3)) #round to 1000s


    for k,Temp in enumerate(Temps):
      specgrid780.append(pickle.load(open("C:/camal/analysis_code/specgrids/specgrid780_T%s.p" \
                                        %Temp,'rb')))
      specgrid823.append(pickle.load(open("C:/camal/analysis_code/specgrids/specgrid823_T%s.p" \
                                        %Temp,'rb')))
      specgrid860.append(pickle.load(open("C:/camal/analysis_code/specgrids/specgrid860_T%s.p" \
                                        %Temp,'rb')))

    specgrid780 = np.stack(specgrid780) # Stack so shape 
    specgrid823 = np.stack(specgrid823) # is nTemp, nPWVgrid, nXgrid
    specgrid860 = np.stack(specgrid860)

    # PWV, X grid arrays
    PWVgrid = np.arange(0,25,0.1)   # # PWV: 0-25 in steps of 0.1
    Xgrid   = np.arange(1,3.5,0.01) # airmass: 1:3.5 steps of 0.01

    # Load temperature grid for night
    Tgrid = np.zeros(len(oe.f823)).astype('int') # T: 0,1,2 etc matching data
    iT= 0
    for obj in oe.schedule:
      if obj in oe.objnames:
        Tgrid[np.where(oe.objnames==obj)[0]] = iT
        iT+=1
                       
    # Pick out airmass indices of observations
    Xinds   = ((np.round(oe.X,2) - 1)/0.01).astype(int)

    # define fluxes
    f1 = 0.965 * specgrid780[Tgrid,:,Xinds]  # 
    f2 = specgrid823[Tgrid,:,Xinds] # This gets A2
    f3 = specgrid860[Tgrid,:,Xinds] # this gets A1

    # Get calibration params A1, A2 
    oe.A1, oe.A2 = getScalar(oe,f3[:,100]/f1[:,100])

    oe.m1 = oe.A2 * specgrid823 / (0.965*specgrid780)
    oe.m2 = oe.A2 * specgrid823 / (oe.A1*specgrid860)

    oe.Tgrid = Tgrid
    oe.Xinds = Xinds

def setupSpline(gtime):
  gridsize = np.int(np.round((gtime[len(gtime)-1]-gtime[1])/(1.5/24.0))) + 1
  knots_in = np.linspace(gtime[1],gtime[len(gtime)-1],gridsize)
  N = len(knots_in) #doesnt include 4 repeats at endpoint and 4 at front    
  front = [0]*4 + gtime[0]- 0.02
  end = [0]*4 + gtime[-1] + 0.02
  knots = np.concatenate((front,knots_in,end))
  return knots


def doQuickExtract(oe):
  """
  Do quick PWV extraction
  """
  # Bin data for this
  c1,c2,c1err,c2err,jdbin = oe.c1,oe.c2, oe.c1*.01,oe.c2*0.01,oe.jd #binData(oe)
  # m1, m2
  m1 = oe.m1[oe.Tgrid,:,oe.Xinds]
  m2 = oe.m2[oe.Tgrid,:,oe.Xinds]
  # PWV, X grid arrays
  PWVgrid = np.arange(0,25,0.1)   # # PWV: 0-25 in steps of 0.1
  pwvsnow1 = np.zeros((len(c1)))
  pwvsnow2 = np.zeros((len(c2)))
  cmdiff   = np.zeros((len(c2)))
  for i in range(len(c1)):
    ipwv1 = np.where(abs(m1[i,:] - c1[i]) == min(abs((m1[i,:] - c1[i]))))[0]
    ipwv2 = np.where(abs(m2[i,:] - c2[i]) == min(abs((m2[i,:] - c2[i]))))[0]
    pwvsnow1[i] = np.mean(PWVgrid[ipwv1])
    pwvsnow2[i] = np.mean(PWVgrid[ipwv2])
    cmdiff[i] = np.mean((m1[i,:] - c1[i])[ipwv1])

  return jdbin,pwvsnow1,pwvsnow2

def lnlike(params, oe,data, err):
    """
    data: (oe.c1, oe.c2)
    err = error on data
    gtime: JD array, oe.jd
    knots: knots of spline fromg getknots
    """
    PWVcoeffs    = np.concatenate((params,np.zeros(4))) # pwv coefficients

    # Create PWV array for time sequence for star at each time step
    PWVs = interpolate.splev(oe.jd, (oe.knots,PWVcoeffs,3),ext=0)
    #PWVs[-1] = PWVs[-2]
    #PWVs[0]  = PWVs[1]

    # Get array indices of PWV .. 10 is bc of grid spacing by 0.1's
    PWVinds = np.round(PWVs*10,0).astype(int)

    # Model
    model = np.concatenate((oe.m1[oe.Tgrid,PWVinds,oe.Xinds],oe.m2[oe.Tgrid,PWVinds,oe.Xinds]))

    # Likelihood stuff
    inv_sig2 = 1.0/(err**2)
    return -0.5*(np.sum((data - model)**2 * inv_sig2))


def lnprior(params):
    n = len(params)
    if np.all((params > 0,params < 25)):
        return 0.0
    return -np.inf

def lnprob(params, oe,data, err):
    lp = lnprior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(params, oe, data, err)


def doExtract(oe):
  """
  Run emcee
  """
  def calc_err(a,b):
    "Compute error for a/b ratio assuming sigma_z is sqrt(z)"
    return (a/b) * np.sqrt( ( np.abs(a)/a**2 + np.abs(b)/b**2) )

  err1 = calc_err(oe.f823,oe.f780)
  err2 = calc_err(oe.f823,oe.f860)

  # Define Args
  data, err, oe.knots = np.concatenate([oe.c1,oe.c2]),\
    np.concatenate([err1,err2]),\
    setupSpline(oe.jd)

  # Starting knot coefficients
  p0 = np.ones(len(oe.knots)) + 5
  ndim, nwalkers = len(oe.knots), 80

  # Get first guess at coefficients using optimization thingy
  nll = lambda *args: -lnlike(*args)
  args = (oe,data,err)
  #  result = op.minimize(nll, p0, method='TNC',jac=False,args=args,bounds=[(0,25)]*ndim)
  #result = op.minimize(nll, p0, method='COBYLA',args=args)

  p0_new = p0 #result["x"]

  # Set up first guesses
  pos = [p0_new+ np.random.randn(ndim) for i in range(nwalkers)]

  # Sample with nwalkers, chain length..?
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)

  nsteps = 2000
  for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
      if (i+1) % 100 == 0:
          print("{0:5.1%}".format(float(i) / nsteps))

  # Without progress
  # sampler.run_mcmc(pos, nsteps)

  # Save samples
  samples = sampler.chain[:, 500:, :].reshape((-1, ndim))
  p_best = np.zeros(ndim)
  for i in range(ndim):
    p_best[i] = np.mean(samples[:,i])

  oe.samples = samples
  oe.p_best = p_best

  # Get PWVs from p_best
  PWVcoeffs    = np.concatenate((p_best,np.zeros(4))) # pwv coefficients

  # Create PWV array for time sequence for star at each time step
  oe.PWVs = interpolate.splev(oe.jd, (oe.knots,PWVcoeffs,3),ext=0)

  # Save Stuff
  PWVinds = np.round(oe.PWVs*10,0).astype(int)
  oe.m1_best = oe.m1[oe.Tgrid,PWVinds,oe.Xinds]
  oe.m2_best = oe.m2[oe.Tgrid,PWVinds,oe.Xinds]
  jdbin,oe.pwvs_c1,oe.pwvs_c2 = doQuickExtract(oe)

  # Make PWV errors
  subsamp = oe.samples[np.random.randint(len(oe.samples), size=1000)]
  PWVall = np.zeros((1000,len(oe.PWVs)))
  for i in range(1000):
    p_temp = np.concatenate((subsamp[i,:],np.zeros(4))) # pwv coefficients
    PWVall[i,:] = interpolate.splev(oe.jd, (oe.knots,p_temp,3),ext=0)


  pwvstd_lo = np.zeros((len(PWVall[0])))
  pwvstd_hi = np.zeros((len(PWVall[0])))
  medPWVs = np.zeros((len(PWVall[0])))
  for i in range(len(PWVall[0])):
    medPWVs[i]   = np.median(PWVall[:,i])
    los          = np.where(PWVall[:,i] < medPWVs[i])[0]
    his          = np.where(PWVall[:,i] > medPWVs[i])[0]
    pwvstd_lo[i] = np.std(PWVall[:,i][los])
    pwvstd_hi[i] = np.std(PWVall[:,i][his])

  oe.PWVserr_lo = pwvstd_lo
  oe.PWVserr_hi = pwvstd_hi

def getKP(oe):
  """
  Download kitt peak and azam PWV data

  Currently just saving kitt peak time, pwv,epwv
  """
  def suominet(filename):
    dat = np.genfromtxt(filename, dtype=None,
      names=('loc', 'date', 'PWV', 'PWVerr', 'null1', 'TZD',
      'null2','press', 'temp', 'hum', 'null3', 'null4', 'null5', 'null6'))

    # Convert times to julian date
    times = Time(dat['date'],scale='utc').jd

    # Return useful things
    return times, dat['PWV'], dat['PWVerr'], dat['temp'], dat['press'],dat['hum']

  year = oe.night[1:5]
  iworked = True
  # Try d/ling azam, kp from suominet site..don't want code to break if site is down
#  try:
#    url_azam = 'http://www.suominet.ucar.edu/data/staYrHr/AZAMnrt_%s.plot' %year
#    urllib.urlretrieve(url_azam,'%s/AZAMnrt_%s.txt' %(oe.finaldataPath,year))
#    azam = suominet('%s/AZAMnrt_%s.txt' %(oe.finaldataPath,year))
#  except IOError:
#    print 'Cannot D/L AZAM File from Suominet'

  url_kitt = 'http://www.suominet.ucar.edu/data/staYrHr/KITTnrt_%s.plot' %year
  urllib.urlretrieve(url_kitt,'%s/KITTnrt_%s.txt' %(oe.finaldataPath,year))

  try:
    kitt = suominet('%s/KITTnrt_%s.txt' %(oe.finaldataPath,year))
  except IOError:
    print 'Cannot D/L KITT File from Suominet'
    iworked = False

  # Get portion corresponding to night's time line
  if iworked == True:
    inight = np.where((kitt[0] > (oe.jd[0]- 0.02)) & (kitt[0] < oe.jd[-1]+0.02))[0]
  # Store
    oe.t_kitt = kitt[0][inight]
    oe.pwv_kitt = kitt[1][inight]
    oe.epwv_kitt = kitt[2][inight]
  # later can add more to save


def getWeather(oe):
  """
  Get the weather to save and compare..reformat it into PWVout fits file
  """
  pass


def saveFits(oe):
  # Variable Names
  tbhdu = pyfits.BinTableHDU.from_columns(
       [pyfits.Column(name='fnames', format='26A', array=np.array(oe.fnames)),
        pyfits.Column(name='JD',format='D',array=np.array(oe.jd)),
        pyfits.Column(name='alt', format='E', array=np.array(oe.alt)),
        pyfits.Column(name='objname', format='16A', array=np.array(oe.objnames)),
        pyfits.Column(name='PWVs',format='E',array=np.array(oe.PWVs)),
        pyfits.Column(name='PWVserr_lo',format='E',array=np.array(oe.PWVserr_lo)),
        pyfits.Column(name='PWVserr_hi',format='E',array=np.array(oe.PWVserr_hi)),
        pyfits.Column(name='PWV1',format='E',array=np.array(oe.pwvs_c1)),
        pyfits.Column(name='PWV2',format='E',array=np.array(oe.pwvs_c2)),
        pyfits.Column(name='m1_best', format='E',array=np.array(oe.m1_best)),
        pyfits.Column(name='m2_best', format='E',array=np.array(oe.m2_best)),
        pyfits.Column(name='c1', format='E',array=np.array(oe.c1)),
        pyfits.Column(name='c2', format='E',array=np.array(oe.c2)),
        pyfits.Column(name='A1', format='E',array=np.array(oe.A1[:,0,0])),
        pyfits.Column(name='A2', format='E',array=np.array([oe.A2])) ])

  filename = oe.finaldataPath + '\\PWVout_%s.fits' %oe.night
  if os.path.exists(filename):
    os.system('rm %s' % filename)
    tbhdu.writeto(filename)
  else:
    tbhdu.writeto(filename)



def PlotPWVs(oe):
  labels = []
  for i in range(len(oe.knots)):
    labels.append('p%s' %i)

  fig = corner.corner(oe.samples, labels=labels,
    truths=oe.p_best)
  fig.savefig(oe.finaldataPath + "\\triangle.png")

  plt.figure()
  plt.plot(oe.jd,oe.c1,'co',mec='None')
  plt.plot(oe.jd,oe.m1_best,'k',lw=3)
  plt.xlabel('JD')
  plt.ylabel('$c_1$')
  plt.savefig(oe.finaldataPath + '\\plot_c1.png')

  
  plt.figure()
  plt.plot(oe.jd,oe.c2,'co',mec='None')
  plt.xlabel('JD')
  plt.ylabel('$c_2$')
  plt.plot(oe.jd,oe.m2_best,'k',lw=3)

  plt.savefig(oe.finaldataPath + '\\plot_c2.png')

  # Get Kitt Peak points
  getKP(oe)

  # Plot PWVs
  plt.figure()
  plt.plot(oe.jd,oe.pwvs_c1,'o',c='grey',zorder=-1)
  plt.plot(oe.jd,oe.pwvs_c2,'o',c='grey',zorder=-1)
  #  plt.plot(oe.jd,oe.PWVs,'k-',lw=2)
  plt.fill_between(oe.jd,oe.PWVs-3*oe.PWVserr_lo,oe.PWVs+3*oe.PWVserr_hi,facecolor='gray',zorder=100,label='CAMAL')
  if len(oe.t_kitt) !=0:
    plt.errorbar(oe.t_kitt,oe.pwv_kitt,oe.epwv_kitt,fmt='o',c='r',zorder=200,label='Kitt Peak')
  plt.legend(loc='best')
  plt.xlabel('JD')
  plt.ylabel('PWV (mm)')
  plt.xlim(oe.jd[0],oe.jd[-1])
  plt.ylim(0,25)

  plt.savefig(oe.finaldataPath + '\\plot_PWV.png')



if __name__ == '__main__':
  # load schedule & reduction list, store things into oe object
  oe = prepNight(email=False) 

  try:
    oe.night
    noNightsLeft = False
  except AttributeError:
    noNightsLeft = True

  if not noNightsLeft:
    # Check if there were even any observations that night
    if len(oe.jdraw) != 0:
    # Load cleaned up reduction data
      loadData(oe)

    else:
      editAnalysisList(mode='w',night=oe.night,error=True)

  if not noNightsLeft:
    # Make sure after cleaning data, didnt nix all points
    if len(oe.jd) !=0:
      # Get Tind, Xind, load color grids
      loadColorGrids(oe)

      # Do Quick PWV fit
      doExtract(oe)

      # Write to file  
      saveFits(oe)

      # Make Plots
      PlotPWVs(oe)

      # Update reduction list saying night is reduced
      editAnalysisList(mode='w',night=oe.night,error=False)

      # Backup data
      uploadData(oe)

    else: #edit analysis file, mark error in that it was a bad night
      editAnalysisList(mode='w',night=oe.night,error=True)

  else:
    pass # if no nights to observe, just dont do anything