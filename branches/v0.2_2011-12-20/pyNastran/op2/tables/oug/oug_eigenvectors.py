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
    def __init__(self,dataCode,iSubcase,eigReal):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        self.updateDt = self.updateEigReal
        #print "eigReal = %s" %(eigReal)
        
        #assert eigReal>=0.
        self.gridTypes = {}
        self.displacements = {self.eigReal: {}}
        self.rotations     = {self.eigReal: {}}

    def updateEigReal(self,dataCode,eigReal):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert eigReal>=0.
        #self.eigReal = int(eigReal)
        self.dataCode = dataCode
        self.applyDataCode()
        self.eigReal = eigReal
        #print "eigReal = %s" %(str(eigReal))
        self.displacements[self.eigReal] = {}
        self.rotations[self.eigReal] = {}

        self.appendDataMember('modes','mode')
        self.appendDataMember('eigrs','eigr')
        self.appendDataMember('eigis','eigi')

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements

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
        #print 'self.eigReal = %s' %(self.eigReal),type(self.eigReal)
        #print "d = ",self.displacements
        self.displacements[self.eigReal][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.eigReal][nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def eigenvalues(self):
        return sorted(self.displacements.keys())

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += '-eigenvalues-\n'
        for i,eigenvalue in enumerate(self.eigenvalues()):
            msg += '%-2s %f\n' %(i,eigenvalue)
        msg += '\n'

        #if self.eigReal is not None:
        #    msg += 'eigReal = %g\n' %(self.eigReal)
        headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for freq,eigVals in sorted(self.displacements.items()):
            msg += 'eigenvalueReal = %f\n' %(freq)
            msg += headerLine
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[freq][nodeID]
                Type = self.gridTypes[nodeID]
                (dx,dy,dz) = displacement
                (rx,ry,rz) = rotation

                msg += '%-8i %8s ' %(nodeID,Type)
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
        return msg

class realEigenVectorObject(scalarObject): # approachCode=2, sortCode=0, thermal=0
    """
                                         R E A L   E I G E N V E C T O R   N O .          1
      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
             1      G      0.0            0.0            0.0            0.0            1.260264E-01   0.0
             7      G      9.999849E-01   0.0            6.728968E-03   0.0            8.021386E-03   0.0
    """
    def __init__(self,dataCode,iSubcase,eigReal):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        #print "eigReal = %s" %(eigReal)
        
        #assert eigReal>=0.
        self.gridTypes = {}
        self.T = {self.eigReal: {}}

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #print "self.dataCode = ",self.dataCode
        #print "dt = ",dt
        #raise Exception(self.dataCode)
        
    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements

        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        elif gridType==7:
            Type = 'L'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.eigReal = %s' %(self.eigReal),type(self.eigReal)
        #print "d = ",self.displacements
        self.T[self.eigReal][nodeID] = v1
    ###

    def eigenvalues(self):
        return sorted(self.T.keys())

    def __repr__(self):
        msg = '---REAL EIGENVECTORS---\n'
        msg += '-eigenvalues-\n'
        for i,eigenvalue in enumerate(self.eigenvalues()):
            msg += '%-2s %f\n' %(i,eigenvalue)
        msg += '\n'

        #if self.eigReal is not None:
        #    msg += 'eigReal = %g\n' %(self.eigReal)
        headers = ['T']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for freq,eigVals in sorted(self.T.items()):
            msg += 'eigenvalueReal = %f\n' %(freq)
            msg += headerLine
            for nodeID,T in sorted(eigVals.items()):
                Type = self.gridTypes[nodeID]
                
                msg += '%-8i %8s ' %(nodeID,Type)

                if abs(T)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3f ' %(T)
                ###
                msg += '\n'
            ###
            msg += '\n'
        return msg

class complexEigenVectorObject(scalarObject): # approachCode=2, sortCode=0, thermal=0
    """
    @todo add table data
    """
    def __init__(self,dataCode,iSubcase,eigReal):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        self.updateDt = self.updateEigReal
        #print "eigReal = %s" %(eigReal)
        
        #assert eigReal>=0.
        self.gridTypes = {}
        self.displacements = {self.eigReal: {}}
        self.rotations     = {self.eigReal: {}}

    def updateEigReal(self,dataCode,eigReal):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert eigReal>=0.
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        self.dataCode = dataCode
        self.applyDataCode()
        #print "eigReal = %s" %(str(eigReal))
        self.displacements[self.eigReal] = {}
        self.rotations[self.eigReal] = {}

        #self.appendDataMember('modes','mode')
        #self.appendDataMember('eigrs','eigr')
        #self.appendDataMember('eigis','eigi')

    def add(self,nodeID,gridType,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i):
        #msg = "nodeID=%s v1=%s v2=%s v3=%s v4=%s v5=%s v6=%s" %(nodeID,v1,v2,v3,v4,v5,v6)
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements

        #if gridType==0:
        #    Type = 'S??'
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.eigReal = %s' %(self.eigReal),type(self.eigReal)
        #print "d = ",self.displacements
        self.displacements[self.eigReal][nodeID] = [[v1r,v1i],[v2r,v2i],[v3r,v3i]] # dx,dy,dz
        self.rotations[self.eigReal][nodeID]     = [[v4r,v4i],[v5r,v5i],[v6r,v6i]] # rx,ry,rz
    ###

    def eigenvalues(self):
        return sorted(self.displacements.keys())

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += '-eigenvalues-\n'
        for i,eigenvalue in enumerate(self.eigenvalues()):
            msg += '%-2s %f\n' %(i,eigenvalue)
        msg += '\n'

        #if self.eigReal is not None:
        #    msg += 'eigReal = %g\n' %(self.eigReal)
        headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for freq,eigVals in sorted(self.displacements.items()):
            msg += 'eigenvalueReal = %f\n' %(freq)
            msg += headerLine
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[freq][nodeID]
                Type = self.gridTypes[nodeID]
                (dx,dy,dz) = displacement
                (rx,ry,rz) = rotation

                msg += '%-8i %8s ' %(nodeID,Type)
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
        return msg

