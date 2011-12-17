from pyNastran.op2.resultObjects.op2_Objects import scalarObject
from pyNastran.op2.op2Errors import *

class StrainEnergyObject(scalarObject):
    """
                               E L E M E N T   S T R A I N   E N E R G I E S
 
    ELEMENT-TYPE = QUAD4      * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   9.817708E+08
    SUBCASE               1   * TOTAL ENERGY OF ALL ELEMENTS IN SET       1 =   4.192036E+08
 
       ELEMENT-ID   STRAIN-ENERGY  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
               12   2.291087E+07        2.3336            2.291087E+02
               13   1.582968E+07        1.6124            1.055312E+02
               14   6.576075E+07        6.6982            3.288037E+02
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.energy  = {}
        self.percent = {}
        self.density = {}
        print self.dataCode
        print "numWide = %s %s"  %(self.dataCode['numWide'],type(self.dataCode['numWide']))
        if self.dataCode['numWide']==4:
            self.getLength    = self.getLength4
            self.add          = self.add4
            self.addTransient = self.addTransient4
        elif self.dataCode['numWide']==5:
            self.getLength    = self.getLength5
            self.add          = self.add5
            self.addTransient = self.addTransient5
        else:
            raise RuntimeError('invalid numWide=%s' %(self.numWide))
        ###
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
            self.isTransient = True
            self.add         = self.addTransient
        ###

    def getLength4(self):
        return(16,'ifff')

    def getLength5(self):
        return(16,'ssssfff')

    def updateDt(self,dataCode,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        print "updating dt...dt=%s" %(dt)
        print "dataCode = ",self.dataCode
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.energy:
            self.energy[self.dt] = {}
            self.percent[self.dt] = {}
            self.density[self.dt] = {}
        
    def add4(self,out):
        (grid,energy,percent,density) = out
        grid = (grid-self.deviceCode)/10
        assert grid not in self.energy
        if grid<=0:
            raise InvalidGridID_Error('grid=%s' %(grid))
        self.energy[grid]  = energy
        self.percent[grid] = percent
        self.density[grid] = density

    def add5(self,out):
        (a,b,c,d,energy,percent,density) = out
        #print "out = ",out
        word = a+b+c+d
        #assert word not in self.energy,'%s in energy...' %(word)
        #if grid<=0:
        #    raise InvalidGridID_Error('grid=%s' %(grid))
        self.energy[word]  = energy
        self.percent[word] = percent
        self.density[word] = density

    def addTransient4(self,out):
        dt = self.dt
        (grid,energy,percent,density) = out
        grid = (grid-self.deviceCode)/10
        #print str(self)
        #assert grid not in self.energy[dt],'grid=%s dt=%s energy=%s percent=%s density=%s' %(grid,dt,energy,percent,density)
        if grid<=0:
            raise InvalidGridID_Error('grid=%s' %(grid))

        self.energy[dt][grid]  = energy
        self.percent[dt][grid] = percent
        self.density[dt][grid] = density
    
    def addTransient5(self,out):
        dt = self.dt
        (a,b,c,d,energy,percent,density) = out
        word = a+b+c+d
        assert word not in self.energy[dt]
        #if grid<=0:
        #    raise InvalidGridID_Error('grid=%s' %(grid))
        self.energy[word]  = energy
        self.percent[word] = percent
        self.density[word] = density

    def __reprTransient__(self):
        msg  = '---Transient Strain Energy Object---\n'
        msg += "%-10s %-14s% -14s% -14s\n" %('EID','Energy','PercentTotal','density')
        for dt,Energy in sorted(self.energy.items()):
            msg += "%s = %s\n" %(self.dataCode['name'],dt)
            for eid,energy in sorted(Energy.items()):
                percent = self.percent[dt][eid]
                density = self.density[dt][eid]
                msg += "%-10i %-14g %-14g %-14g" %(eid,energy,percent,density)+'\n'
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg  = '---Strain Energy Object---\n'
        msg += "%-10s %-14s% -14s% -14s\n" %('EID','Energy','PercentTotal','density')
        for eid,energy in sorted(self.energy.items()):
            percent = self.percent[eid]
            density = self.density[eid]
            msg += "%-14i %-14g %-14g %-14g" %(eid,energy,percent,density)+'\n'
        return msg
