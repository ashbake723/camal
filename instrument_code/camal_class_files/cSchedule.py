import numpy as np
import matplotlib.pylab as plt
import astropy
from astropy.io.votable import parse_single_table
import time
from astropy.time import Time
import os
import cTargets
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.utils import iers
import datetime
import jdcal
import astroplan
import astropy.units as u
import mail
#iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)


####### TO DO ############
# include way to recompute schedule in case of clouds (mask out target that didn't work?
# or give up and wait an hour?

class scheduler():
    def __init__(self,site):
        """
        Set up observing site using config file
        Have it take config info since whipple is probs summit
        Also sunset and sunrise should already be loaded for the site so
        ...just pass site to init and also config
        Switch to pyephem to match minerva
        """
        # Load Data
        self.dat = cTargets.loadTargets()
        self.night = site.night

        # Split time window into 30 increments
        # convert datetime into JD to do calc, then convert back using delta added to sunset to define times back into datetimes
        DAY  = datetime.timedelta(1)
        J2000_JD = datetime.timedelta(2451545)
        JULIAN_EPOCH = datetime.datetime(2000, 1, 1, 12)

        dt = site.sunset()
        self.JDsunset = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.0 + dt.minute/(24.*60) + dt.second/(24.*3600.0)

        dt = site.sunrise()
        self.JDsunrise = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) + dt.hour/24.0 + dt.minute/(24.*60) + dt.second/(24.*3600.0)

        self.time_window = self.JDsunset + (self.JDsunrise - self.JDsunset) * \
            np.linspace(0, 1, 30)


    def calcTargetLimits(self,X=1.5):
        """
        Calculate the limites of each target candidates
        
        Input:
        X -  airmass cutoff limit

        Computes:
        S - Start time target is above X (units: JD) [Ntargs]
        F - Set time of target below X   (units: JD) [Ntargs]
        alt - array of alts for target   [Ntimesteps, Ntargs]
        """
        Ntargs = len(self.dat.ids)
        # Set up arrays
        alt = np.zeros((len(self.time_window),Ntargs))   # altitudes of targets
        S = np.zeros(Ntargs)   # Start time above X
        F = np.zeros(Ntargs)   # Time sets below X
        for star in range(0,Ntargs):
            tempcoords = SkyCoord(self.dat.RA[star],self.dat.DEC[star],
                                      frame='icrs',unit=(u.hourangle,u.deg))
            temptarg = astroplan.FixedTarget(name=self.dat.ids[star], coord=tempcoords)
            self.whipple = astroplan.Observer.at_site('whipple')
            alt[:,star] = self.whipple.altaz(Time(self.time_window,format='jd'), temptarg).alt
            max_alt = np.arccos(1.0/X) * 180.0/np.pi
            star_up = np.where(alt[:,star] > 40.0)[0]
            if len(star_up) != 0:
                S[star] = self.time_window[star_up][0]   #star starts being up
                F[star] = self.time_window[star_up][-1]  #star sets below 30 deg
            else:
                S[star]=9.e9
                F[star]=-9.e9
        self.S = S
        self.F = F
        self.alt = alt

    def pickStar(self,t0, dt, badlist=None):
        """
        Pick best star above airmass 1.5 after t0 for at least dt hours
        
        Inputs:
        t0 :   Start time needs to be up by     [JD]
        dt :   Length of time needs to be up    [days]
        badlist: List of targets to not consider [list]  optional

        Uses:
        S  :   Rising times of targets          [JD]
        F  :   Setting times of targets         [JD]
        dat:   data class with target info      class structure
        
        Outputs: 
        index of the star to observe at that time
        if -1 then didn't find anything and you should reduce X limits or dt
        """
        # Remove targets because of moon or clouds somehow
        if badlist != None:
            pass #remove ids in badlist
        
        # Which targets are up at t0
        i_up     = np.where((self.S <= t0) & (self.F > t0))[0]
        index    = np.arange(len(self.S))[i_up] # keep track of index
        # mask : zero means not observable at start time for min dt
        mask     = np.zeros(len(self.S),dtype=int)

        # Which targets will be up for dt days longer (dt ~ 1-4 hours)
        longUp = np.where(self.F[index] > (t0 + dt))[0]
        index1  = index[longUp]
        mask[longUp] = 1

        if sum(mask) == 0:
            return -1 # try reducing dt or airmass constraint

        # Pick brightest star
        brightest = np.where(self.dat.imag[np.where(mask == 1)[0]] < 1.8)[0]
        if len(brightest) == 0:
            brightest = np.where(self.dat.imag[index1] == min(self.dat.imag[index1]))[0]
            index3 = index1[brightest]
        # Pick best spectral type
        elif len(brightest) > 1:
            index2 = index1[brightest]
            hottest = np.where(self.dat.specnum[index2] == \
                                   np.min(self.dat.specnum[index2]))
            index3 = index2[hottest]
        else:
            index3 = index1[brightest]
            
        return index3[0]


    def run(self):
        # Use pickstar
        self.calcTargetLimits()
        
        ischedule = []  # object index
        tstarts   = []  # times to start the object
        tends     = []
        dt = 0.15
        t0 = self.JDsunset + 0.025
        error = False
        while t0 < self.JDsunrise - 0.025: # if .025 days left just finish target
            ischedule_temp = self.pickStar(t0, dt)
            if ischedule_temp == -1:
                dt = dt - 0.025  # reduce time span and try again
                print dt
            else:
                tstarts.append(t0)
                ischedule.append(ischedule_temp)
                t0 = self.F[ischedule][-1] # redefine t0 as end for that target
                tends.append(t0)  # save end time for that target
                print ischedule
            # If statement for when close to end of night
            if t0 - (self.JDsunset - 0.025) < 0.05:
                tends[-1] = self.JDsunrise - 0.025
                break # Just finish on the previous target til sunrise
            # For when dt is too small, reduce X limit and try again
            if (dt < 0.05) & (error == False) & (ischedule_temp == -1):
                self.calcTargetLimits(X=1.8)
                dt = 0.15   # reduce airmass limits
                error=True
            elif (dt < 0.05) & (error == True):
                print 'eek shouldnt be here - make send email'
                body = 'Please check shcedule for night %s '\
                        'becuase the code couldnt finish '\
                        'running the scheduler to find targets'\
                        'for the entire night\n Love CAMAL'
                mail.send('Error in Scheduler',body)
                break


        self.ischedule = ischedule
        self.tstarts   = tstarts
        self.tends     = tends

        return ischedule, self.dat.ids[ischedule]
        
    def save_schedule(self):
        path  =  'C:\camal\schedule\\'
        f = open(path + self.night, 'w')

        f.write('# Object\t UTstart\t UTendtime\t name\t RA\t DEC\t Imag\t Temp\t \n')
        for i, idd in enumerate(self.ischedule):
            f.write(self.dat.ids[idd] + '\t')
            f.write(Time(self.tstarts[i],format='jd').isot.split('T')[1] + '\t')

            f.write(Time(self.tends[i],format='jd').isot.split('T')[1] + '\t')

            f.write(self.dat.names[idd] + '\t')
            f.write(self.dat.RA[idd]    + '\t')
            f.write(self.dat.DEC[idd]   + '\t')
            f.write(str(self.dat.imag[idd])  + '\t')
            f.write(str(self.dat.temps[idd]) + '\n')

        f.close()

    def plot_schedule(self,night):
        # Folder to save plot
        path = 'C:\camal\schedule\\'

        # Plot results
        plt.figure()
        for i in range(len(self.ischedule)):
            plt.plot(self.time_window,self.alt[:,ischedule[i]],\
                        label=str(self.dat.ids[ischedule[i]]))
            plt.plot([self.tstarts[i],self.tstarts[i]],[0,100],'k-',lw=2)
        plt.xlabel('Time (JD)')
        plt.ylabel('Altitude Angle (deg)')
        plt.legend(loc='best')
        plt.ylim(0,100)
        plt.savefig(path + '%s.png' %self.night)

#np.savetxt('schedule_test.txt',schedule.T,header = 'order , id, imag, ra, dec',fmt='%s13')








