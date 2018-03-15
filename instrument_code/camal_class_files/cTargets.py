import numpy as np



class loadTargets():
    def __init__(self):
        # Load Full Target List...make a class                       
        filename= 'C:/camal/instrument_code/camal_class_files/TargetList_v1.1.txt'
        dat = np.genfromtxt(filename, dtype=None,
              names=('num','name', 'id','temps','typ','coord1' ,'MagU' ,'MagB','MagV','MagR', 'MagI',
                     'spectype','bib','not'),delimiter='|')
        self.rawdat = dat
        coords = np.zeros((len(dat['coord1']),2),dtype='|S14')
        ids = []
        names = []
        for row in range(0,len(dat['coord1'])):
            coords[row] = dat['coord1'][row].strip().split(' ')
            ids.append(dat['id'][row].strip()[2::].replace(' ', '_'))
            names.append(dat['name'][row].strip()[5:])
        self.coords = coords
        self.RA = coords[:,0]
        self.DEC= coords[:,1]
        self.imag = dat['MagI']
        self.ids = np.array(ids)
        self.spectype = dat['spectype']
        self.names = np.array(names)
        self.temps = dat['temps']
        
        # Convert spectypes to numbers
        specnum = np.zeros(len(dat['spectype']),dtype='int')
        specmap = {'O':0, 'B':1, 'A':2, 'F':3}
        for i,stype in enumerate(dat['spectype']):
            specnum[i] = specmap[stype.strip()[0]]
            print i
        self.specnum = specnum

