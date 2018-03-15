
import sys, os
import numpy as np
import time
import matplotlib.pylab as plt
import cPickle as pickle
from scipy import interpolate
# Load emcee - had to switch b/c mc3 is linux only
import emcee

import urllib

# ---------------------------
night    = 'n20170520'
filename = 'camal_' + night + '.txt'
dataPath = '/data3/planetgroup/ashbaker/camal/camal120716/final/'
A2 = 1.0
# ----------------------------

print 'night = ' + night

sys.path.append('../tools/')
import toolbox
import camal as cm

lf = toolbox.loadfile()
tl = toolbox.tools()

import loadCAMAL as lc
import loadCAMAL_old as lc # if date after june 2016

def getA1(objects):
  scalarr = np.ones(len(gtime))
  for jj, Ti in enumerate(np.unique(Tgrid)):
    scalarr[np.where(Tgrid == Ti)[0]] = params[jj]
  A1 = [scalarr]*250
  A1 = np.array(A1).T
  return A1


class fitPWV():
  """
  Find PWV best fit using emcee
  """
  def __init__(self,oe):
    self.oe = oe

    # Colors
    c1 = oe.f823/oe.f780
    c2 = oe.f823/oe.f860

    # Make Data Array
    ibreak = len(oe.f780)
    data = np.concatenate((c1,c2))
 
    # Make Uncertainties by binning data
    unc = np.ones(len(data)) # in future switch to binning

    # Load temperature grid for night
    Temps = []
    Tgrid = np.zeros(len(oe.f823)).astype('int') # T: 0,1,2 etc matching data
    iT= 0
    for obj in oe.objects:
      if obj in oe.objnames:
        Temps.append(oe.schedic[obj]['Temp'])
        Tgrid[np.where(oe.objnames==obj)[0]] = iT
        iT+=1

    # Load Spectral Grids
    specgrid780 = []
    specgrid823 = []
    specgrid860 = []
    for k,Temp in enumerate(Temps):
      specgrid780.append(pickle.load(open("../specgrids/specgrid780_T%s.p" \
                                        %Temp,'rb')))
      specgrid823.append(pickle.load(open("../specgrids/specgrid823_T%s.p" \
                                        %Temp,'rb')))
      specgrid860.append(pickle.load(open("../specgrids/specgrid860_T%s.p" \
                                        %Temp,'rb')))

    specgrid780 = np.stack(specgrid780) # Stack so shape 
    specgrid823 = np.stack(specgrid823) # is nTemp, nPWVgrid, nXgrid
    specgrid860 = np.stack(specgrid860)

    # PWV, X grid arrays
    PWVgrid = np.arange(0,25,0.1)   # # PWV: 0-25 in steps of 0.1
    Xgrid   = np.arange(1,3.5,0.01) # airmass: 1:3.5 steps of 0.01
                       
    # Pick out airmass indices of observations
    Xinds   = ((np.round(oe.X,2) - 1)/0.01).astype(int)

    # Define calibration parameters A1, A2
    A2 = 0.979    # Scalar
    A1 = getA1(oe.schedule)   # array length flux arrays

    # Make ratios for stars
    f1 = 0.965 * specgrid780[Tgrid,:,Xinds]  # 
    f2 = A2 * specgrid823[Tgrid,:,Xinds]
    f3 = A1 * specgrid860[Tgrid,:,Xinds]

    self.m1 = f2/f1
    self.m2 = f2/f3

# Make coefficients be close to one at first (will all be one?) then also give some
# PWV to multiply by the fxn of pwv of time
# keep knots as once every hour

def fxn(self,params,knots,gtime):
    """
    Airmass
    """
    NTs          = len(np.unique(Tgrid))
    PWVcoeffs    = np.concatenate((params,np.zeros(4))) # pwv coefficients
    # Create PWV array for time sequence for star at each time step
    PWVs = interpolate.splev(gtime, (knots,PWVcoeffs,3),ext=0)
    PWVs[-1] = PWVs[-2]
    PWVs[0]  = PWVs[1]
    # Get array indices of PWV
    PWVinds = np.round(PWVs*10,0).astype(int)
    # Define Scalar array
    scalarr = np.ones(len(gtime))
    return np.concatenate((self.m1[PWVinds], self.m2[PWVinds]))

# Set up Spline knots ( make past endpoints )
  def run(self,gtime):
    # DEFINE KNOTS
    gridsize = np.int(np.round((gtime[len(gtime)-1]-gtime[1])/(1.0/24.0))) + 1
    knots_in = np.linspace(gtime[1],gtime[len(gtime)-1],gridsize)
    N = len(knots_in) #doesnt include 4 repeats at endpoint and 4 at front    
    front = [0]*4 + gtime[0]- 0.02
    end = [0]*4 + gtime[-1] + 0.02
    knots = np.concatenate((front,knots_in,end))

    # DEFINE STARTING PARAMS
    params   = np.ones(N+4)*5.0
    pmin     = np.zeros(N+4)
    pmax     = np.ones(N+4)*24
    stepsize = np.zeros(N+4)+0.1
    prior    = np.ones(N+4)
    priorlow = np.ones(N+4)-0.8
    priorup  = np.ones(N+4)+16.0

    # Run the MCMC:
    bestp, CRlo, CRhi, stdp, posterior, Zchain = mc3.mcmc(data, unc,
             func=fxn, indparams=[Xinds,knots,gtime], \
             params=params, stepsize=stepsize, pmin=pmin, pmax=pmax, \
             prior=prior,priorlow=priorlow,priorup=priorup, \
             nsamples=10**5, burnin=10**4,plots=False)



    return bestp



class plotPWV():
    def __init__(self,results):
        pass

    def plotSpline(self):
        pass


# Overplot PWV from water vapor monitor at amado
year = night[1:5]
url_azam = 'http://www.suominet.ucar.edu/data/staYrHr/AZAMnrt_%s.plot' %year
url_kitt = 'http://www.suominet.ucar.edu/data/staYrHr/KITTnrt_%s.plot' %year
urllib.urlretrieve(url_azam,'suominet/AZAMnrt_%s.txt' %year) 
urllib.urlretrieve(url_kitt,'suominet/KITTnrt_%s.txt' %year)
azam = lf.suominet('suominet/AZAMnrt_%s.txt' %year)
kitt = lf.suominet('suominet/KITTnrt_%s.txt' %year)
ax[0].errorbar(24*(azam[0]-gtime[0]),*azam[1:3],label='Amado',fmt='o',color='blue',ecolor='blue',ms=0,elinewidth=4,capsize=6)
ax[0].errorbar(24*(kitt[0]-gtime[0]),*kitt[1:3],label='Kitt Peak',fmt='o',color='red',ecolor='red',ms=0,elinewidth=4,capsize=6)

if night == 'n20161204':
    leg = ax[0].legend(loc='best',prop={'size':20})
    #leg.get_frame().set_alpha(0)
    leg2 = ax[1].legend(loc=4,prop={'size':20})
    #leg2.get_frame().set_alpha(0)

# Formatting
plt.rcParams.update({'font.size': 22})
plt.tick_params('both', length=5, width=1, which='minor')
plt.tick_params('both', length=10, width=1, which='major')
ax[1].tick_params(axis='x', pad=8)
#yticks2 = [0.24,0.26,0.28]
yticks2 = [2.6,2.8,3.0]
ax[2].set_yticks(yticks2)
yticks1 = [0.60,0.65,0.70,0.75]
ax[1].set_yticks(yticks1)
yticks0 = [0,4,8,12,16]
ax[0].set_yticks(yticks0)
plt.subplots_adjust(bottom=0.14, left = 0.15,hspace=0)
#ax[0].set_title(night)
ax[0].text(0.5,14,night)
ax[0].set_ylim(0,18)
ax[0].grid(True)
ax[1].grid(True)
# Save Plot
plt.savefig('plots/PWV_withGPS_%s_noA2.eps' %(night),transparent=False)



# ------------ PWVnow ------------------- #
scalarr = np.ones(len(gtime))
for jj, Ti in enumerate(np.unique(Tgrid)):
    scalarr[np.where(Tgrid == Ti)[0]] = params[jj]

A1 = [scalarr]*250
A1 = np.array(A1).T
f1 = 0.965 * specgrid780[Tgrid,:,Xinds]
f2 = A2 * specgrid823[Tgrid,:,Xinds]
f3 = A1*specgrid860[Tgrid,:,Xinds]
m1 = f2/f1
m2 = f2/f3
c1 = np.roll(f823,-1)/f780
c2 = f823/f860
pwvsnow1 = np.zeros((len(c1)))
pwvsnow2 = np.zeros((len(c2)))
cmdiff   = np.zeros((len(c2)))
for i in range(len(c1)):
    ipwv1 = np.where(abs(m1[i,:] - c1[i]) == min(abs((m1[i,:] - c1[i]))))[0]
    ipwv2 = np.where(abs(m2[i,:] - c2[i]) == min(abs((m2[i,:] - c2[i]))))[0]
    pwvsnow1[i] = np.mean(PWVgrid[ipwv1])
    pwvsnow2[i] = np.mean(PWVgrid[ipwv2])
    cmdiff[i] = np.mean((m1[i,:] - c1[i])[ipwv1])

fig, ax = plt.subplots(1,figsize=(8,4),sharex=True)
ax.plot(hours,pwvsnow1,'.',c='gray',zorder=-10,label='CAMAL')
ax.plot(hours,pwvsnow2,'.',c='gray',zorder=-10)
ax.errorbar(24*(azam[0]-gtime[0]),azam[1],azam[2]/2.0,label='Amado',fmt='o',color='blue',ecolor='blue',ms=0,elinewidth=4,capsize=6)
ax.errorbar(24*(kitt[0]-gtime[0]),kitt[1],kitt[2]/2.0,label='Kitt Peak',fmt='o',color='red',ecolor='red',ms=0,elinewidth=4,capsize=6)
plt.xlim(min(hours),max(hours))
plt.ylim(min(pwvsnow1),24)
plt.xlabel('Time from start (hr)')
plt.ylabel('PWVs (mm)')
if night == 'n20161204':
    leg = ax.legend(loc='best',prop={'size':18})

ax.fill_between(hours,medPWVs + 3*pwvstd_hi,medPWVs-3*pwvstd_lo,facecolor='gray')

plt.subplots_adjust(bottom=0.19)
plt.text(1,20, night)
plt.savefig('paperplots/outputs/CAMALGPS_pwvnow_%s.eps' %night)

# Save PWVs
plt.savetxt('pwvs/PWVnow_%s_noA2.txt' %night, zip(gtime,pwvsnow1, pwvsnow2),header='time(JD)   PWVSnow1 (mm) PWVSnow2')


# Plot against e/o
#--------------------------------------
# CAMAL VS GPS
ilap = np.where(((kitt[0]-gtime[0]) < (gtime[-1]-gtime[0])) & ((kitt[0]-gtime[0]) > 0))[0]
pwv2 = interpolate.splev(kitt[0],(knots,PWVcoeffs,3))
ilap3 = np.where(((azam[0]-gtime[0]) < (gtime[-1]-gtime[0])) & ((azam[0]-gtime[0]) > 0))[0]
pwv3 = interpolate.splev(azam[0],(knots,PWVcoeffs,3))

plt.figure(100)
plt.plot(pwv2[ilap],kitt[1][ilap],'bo',label='Kitt Peak')
plt.plot(pwv3[ilap3],azam[1][ilap3],'g+',label='Amado')
plt.plot([0,20],[0,20],'k--')
#plt.legend()
plt.xlabel('PWV CAMAL (mm)')
plt.ylabel('GPS PWV (mm)')
plt.subplots_adjust(bottom=0.15)
#plt.plot(kitt[1][ilap],azam[1][ilap3],'ro')

# bin camal pwvsnow to match
#kittpeak
t_kitt = kitt[0][ilap]
pwvsnow_kitt = np.zeros(len(t_kitt))
for i in range(len(t_kitt)):
    selection = np.where((gtime > t_kitt[i]-0.005) & (gtime < t_kitt[i]+0.005))[0]
    pwvsnow_kitt[i] = np.median(pwvsnow1[selection])
    print np.std(pwvsnow1[selection])

t_azam = azam[0][ilap3]
pwvsnow_azam = np.zeros(len(t_azam))
for i in range(len(t_azam)):
    selection = np.where((gtime > t_azam[i]-0.005) & (gtime < t_azam[i]+0.005))[0]
    pwvsnow_azam[i] = np.median(pwvsnow1[selection])

np.savetxt('pwvs/both_pwvs_out_%s_noA2.txt' %night,zip(kitt[0][ilap],pwv2[ilap],kitt[1][ilap],pwvsnow_kitt, azam[0][ilap3],pwv3[ilap3],azam[1][ilap3], pwvsnow_azam),
           header='t_kitt pwv_camal pwv_kitt pwvnow_kitt t_azam  pwv_camal  pwv_azam pwvnow_azam')



# ------------- SAVE PWVs ------------------- #
# save PWVs pwvstd gtime to text file
f = open('pwvs/camal_out_%s_noA2.txt' %night,'wn')
f.write('# night: %s \n' %night)
f.write('# ncoeff: %s \n' %len(PWVcoeffs))
f.write('# coeffs: %s \n' %PWVcoeffs)
f.write('# Time (JD) \t PWV (mm) \t PWVunclo \t PWVunchi (mm) \n')
for i in range(0,len(gtime)):
    line = np.str(gtime[i])+'\t'+np.str(medPWVs[i])+'\t'+ \
                  np.str(pwvstd_lo[i])+'\t' + np.str(pwvstd_hi[i]) + '\n'
    f.write(line)

f.close()