## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class eigenVectorObject(scalarObject): # approachCode=2, sortCode=0, thermal=0
    """
    EIGENVALUE =  6.158494E+07
        CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

    POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
           1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
        2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
        2003      G      0.0            0.0            0.0            0.0            0.0            0.0
    """
    def __init__(self,dataCode,iSubcase,mode):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.caseVal = mode
        self.updateDt = self.updateMode
        #print "mode = %s" %(mode)
        #print "dataCode = ",self.dataCode
        self.setDataMembers()
        
        #assert mode>=0.
        self.gridTypes = {}
        self.translations = {self.caseVal: {}}
        self.rotations    = {self.caseVal: {}}
    
    def readF06Data(self,dataCode,data):
        iMode = dataCode['mode']
        if iMode not in self.translations:
            self.updateMode(dataCode,iMode)
        
        for line in data:
            (nid,gridType,t1,t2,t3,r1,r2,r3) = line
            self.gridTypes[nid] = gridType
            self.translations[iMode][nid] = array([t1,t2,t3])
            self.rotations[iMode][nid]    = array([r1,r2,r3])
        ###

    def updateMode(self,dataCode,mode):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert mode>=0.
        self.dataCode = dataCode
        self.applyDataCode()
        self.caseVal = mode
        #print "mode = %s" %(str(mode))
        self.translations[self.caseVal] = {}
        self.rotations[self.caseVal] = {}
        self.setDataMembers()

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations

        #if gridType==0:
        #    Type = 'S'
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        elif gridType==7:
            Type = 'L'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.caseVall = %s' %(self.caseVal),type(self.caseVal)
        #print "d = ",self.translations
        self.translations[self.caseVal][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.caseVal][nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def eigenvalues(self):
        return self.eigrs

    def writeF06(self,header,pageStamp,pageNum=1):
        """
        EIGENVALUE =  6.158494E+07
            CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
            2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
            2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        """
        msg = []
        for i,(iMode,eigVals) in enumerate(sorted(self.translations.items())):
            msg += header
            freq = self.eigrs[i]
            msg.append('%16s = %13E\n' %('EIGENVALUE',freq))
            msg.append('%16s = %13E         R E A L   E I G E N V E C T O R   N O . %10i\n \n' %('CYCLES',self.modeCycle,iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx,dy,dz) = displacement
                (rx,ry,rz) = rotation
                
                vals = [dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum += 1
        ###
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'        
        msg += self.printDataMembers()
        name = self.dataCode['name']

        headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for i,(iMode,eigVals) in enumerate(sorted(self.translations.items())):
            freq = self.eigrs[i]
            msg += '%s = %g\n' %(name,iMode)
            msg += 'eigenvalueReal = %g\n' %(freq)
            #msg += 'eigenvalueReal = %f\n' %(freq)
            msg += headerLine
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx,dy,dz) = displacement
                (rx,ry,rz) = rotation

                msg += '%-8i %8s ' %(nodeID,gridType)
                vals = [dx,dy,dz,rx,ry,rz]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'
            ###
            msg += '\n'
            #print msg
            #return msg
        return msg

class realEigenVectorObject(scalarObject): # approachCode=2, sortCode=0, thermal=0
    """
                                         R E A L   E I G E N V E C T O R   N O .          1
      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
             1      G      0.0            0.0            0.0            0.0            1.260264E-01   0.0
             7      G      9.999849E-01   0.0            6.728968E-03   0.0            8.021386E-03   0.0
    """
    def __init__(self,dataCode,iSubcase,mode):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.caseVal = mode
        #print "mode = %s" %(mode)
        self.caseVal = self.getUnsteadyValue()
        self.setDataMembers()
        
        #assert mode>=0.
        self.gridTypes = {}
        self.translations = {self.caseVal: {}}
        self.rotations    = {self.caseVal: {}}

    def updateDt(self,dataCode,dt):
        #print " self.dataCode = ",self.dataCode
        self.dataCode = dataCode
        self.applyDataCode()
        self.setDataMembers()
        self.caseVal = dt
        
        #print "*self.dataCode = ",self.dataCode
        self.translations[self.caseVal] = {}
        self.rotations[self.caseVal]    = {}
        #print "dt = ",dt
        #raise Exception(self.dataCode)
        
    def deleteTransient(self,dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        msg += "           v4=%s v5=%s v6=%s" %(v4,v5,v6)
        #print msg
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations

        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        elif gridType==7:
            Type = 'L'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.caseVal = %s' %(self.caseVal),type(self.caseVal)
        #print "d = ",self.translations
        self.translations[self.caseVal][nodeID] = [v1,v2,v3]
        self.rotations[self.caseVal][nodeID]    = [v4,v5,v6]
    ###

    def modes(self):
        return sorted(self.translations.keys())

    def eigenvalues(self):
        return self.eigrs

    def writeF06(self,header,pageStamp,pageNum=1):
        """
        EIGENVALUE =  6.158494E+07
                                           R E A L   E I G E N V E C T O R   N O .          1

        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
            2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
            2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        """
        msg = []
        #print self.dataCode
        for i,(iMode,eigVals) in enumerate(sorted(self.translations.items())):
            msg += header
            freq = self.eigrs[i]
            msg.append('%16s = %12E\n' %('EIGENVALUE',freq))
            msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' %(iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID,translation in sorted(eigVals.items()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation
                
                vals = [dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum += 1
        ###
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        msg = '---REAL EIGENVECTORS---\n'
        msg += self.printDataMembers()
        name = self.dataCode['name']

        headers = ['T']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for iMode,eigVals in sorted(self.translations.items()):
            msg += '%s = %s\n' %(name,iMode)
            msg += headerLine
            for nodeID,translation in sorted(eigVals.items()):
                Type = self.gridTypes[nodeID]
                
                rotation = self.rotations[iMode][nodeID]
                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation
                
                vals = [dx,dy,dz,rx,ry,rz]
                msg += '%-8i %8s ' %(nodeID,Type)
                for v in vals:
                    if abs(v)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3f ' %(v)
                    ###
                ###
                msg += '\n'
            ###
            msg += '\n'
        return msg

class complexEigenVectorObject(scalarObject): # approachCode=2, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,mode):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.caseVal = mode
        self.updateDt = self.updateMode
        self.setDataMembers()
        
        #print "mode = %s" %(mode)
        
        #assert mode>=0.
        self.gridTypes = {}
        self.translations = {self.caseVal: {}}
        self.rotations    = {self.caseVal: {}}

    def updateMode(self,dataCode,mode):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert mode>=0.
        self.caseVal = mode
        self.dataCode = dataCode
        self.applyDataCode()
        self.setDataMembers()
        #print "mode = %s" %(str(mode))
        self.translations[self.caseVal] = {}
        self.rotations[self.caseVal] = {}

    def add(self,nodeID,gridType,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i):
        #msg = "nodeID=%s v1=%s v2=%s v3=%s v4=%s v5=%s v6=%s" %(nodeID,v1,v2,v3,v4,v5,v6)
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations

        #if gridType==0:
        #    Type = 'S??'
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.caseVal = %s' %(self.caseVal),type(self.caseVal)
        #print "d = ",self.translations
        self.translations[self.caseVal][nodeID] = [v1r,v1i,v2r,v2i,v3r,v3i] # dx,dy,dz
        self.rotations[self.caseVal][nodeID]    = [v4r,v4i,v5r,v5i,v6r,v6i] # rx,ry,rz
    ###

    def eigenvalues(self):
        return sorted(self.translations.keys())

    def writeF06(self,header,pageStamp,pageNum=1):
        msg = []
        #print self.dataCode
        for i,(iMode,eigVals) in enumerate(sorted(self.translations.items())):
            msg += header
            freq = self.eigrs[i]
            #freq = 0.0
            msg.append('%16s = %12E\n' %('EIGENVALUE',freq))
            msg.append('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' %('CYCLES',self.modeCycle,iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (v1r,v2r,v3r,v1i,v2i,v3i) = displacement
                (v4r,v5r,v6r,v4i,v5i,v6i) = rotation
                
                vals = [v1r,v2r,v3r,v1i,v2i,v3i,v4r,v5r,v6r,v4i,v5i,v6i]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [v1r,v2r,v3r,v1i,v2i,v3i,v4r,v5r,v6r,v4i,v5i,v6i] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,v1r,v2r,v3r,v4r,v5r,v6r.rstrip()))
                msg.append('%14s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %('','',          v1i,v2i,v3i,v4i,v5i,v6i.rstrip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum += 1
        ###
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += self.printDataMembers()

        headers = ['T1','T2','T3','R1','R2','R3']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'
        name = self.dataCode['name']

        for i,(iMode,eigVals) in enumerate(sorted(self.translations.items())):
            msg += '%s = %g\n' %(name,iMode)
            msg += headerLine
            for nodeID,translation in sorted(eigVals.items()):
                rotation = self.rotations[iMode][nodeID]
                Type = self.gridTypes[nodeID]
                #(dx,dy,dz) = displacement
                #(rx,ry,rz) = rotation

                msg += '%-8i %8s ' %(nodeID,Type)
                vals = translation+rotation
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'
            ###
            msg += '\n'
        return msg

