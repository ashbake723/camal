##------------------------------------------------------------------------------
## Class: cSchedule
##------------------------------------------------------------------------------
import numpy as np
import time
import astroplan
from astropy.coordinates import SkyCoord
from astroplan import FixedTarget
from astropy.time import Time

__all__ = ["loadStars", "objectInfo", "loadSchedule"]

class loadStars():
    def __init__(self):
        # Load Full Target List...make a class                                                 
        filename     = 'Targets_short.txt'
        dat          = np.genfromtxt(filename, dtype=None,
                      names=('num', 'id','typ','coord1' ,'MagU' ,'MagB','MagV','MagR', 'MagI',
                     'spectype','bib','not'),delimiter='|')
        self.rawdat  = dat
        coords       = np.zeros((len(dat['coord1']),2))
        id           = []
        # Reformat Coordinates and ids
        for row in range(0,len(dat['coord1'])):
            coords[row] = dat['coord1'][row].strip().split(' ')
            id.append(dat['id'][row].strip()[2::].replace(' ', '_'))
        # FIll in rest of info
        self.coords   = coords
        self.RA       = coords[:,0]
        self.DEC      = coords[:,1]
        self.rmag     = dat['MagR']
        self.spectype = dat['spectype']
        self.id       = np.array(id)


class pickTargets():
    def __init__(self):
        self.info = ' eventually code up sweet algorithm from tanyas alg. class'
        # save the targets to file in night loaded by loadSchedule
        #np.savetxt('schedule_test.txt',schedule.T,header = 'order , id, rmag, ra, dec',fmt='%s13')


class objectInfo():
    def __init__(self,RA,DEC):
        ' Give coordinates of object, use astroplan to get alt and az at current time'
        whipple  = astroplan.Observer.at_site('whipple')
        t        = Time( time.strftime("%Y-%m-%d %H:%M:%S.000", time.gmtime()), format = 'iso')
        coords   = SkyCoord(RA, DEC, frame='icrs',unit="deg")
        temptarg = FixedTarget(name='', coord=coords)
        self.alt = whipple.altaz(t, temptarg).alt
        self.az  = whipple.altaz(t, temptarg).az
        
        
class loadSchedule():
    def __init__(self,filename='11_2017'):
    """
    Input: 
        filename:   name of file containing targets (e.g. 11_2017)
    Location:
        Schedule file exists as/in instrument_code\schedule\MM_YYYY
    Returns: (so Sched. file should contain):
        target:   target's ID (e.g. alf_Lyr)
        lststart: LST start time, when telescope should slew to it (e.g. 00:31:00)
        name:     common name (e.g. Vega)
    """
    to_file = '..\\schedule\\'
    f       = open(tofile + filename,'r')
    lines   = f.readlines()
    f.close()
    # Read lines and store
    target   = []
    lststart = []
    name     = [] 
     for l in lines:
        if l.startswith('#'):
            pass
        else:
            target.append(l.split()[0])
            lststart.append(l.split()[1])
            name.append(l.split()[2])
    self.target   =  target
    self.lststart =  lststart
    self.name     =  name