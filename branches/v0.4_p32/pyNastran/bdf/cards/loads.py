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
#import sys
from numpy import array,cross
from numpy.linalg import norm

from .baseCard import BaseCard

class Load(BaseCard):
    """defines the DefaultLoad class"""
    type = 'DefLoad'
    def __init__(self,card,data):
        pass
        #self.type = card[0]

    #def normalize(self,v):
    #    #print "v = ",v
    #    return v/norm(v)

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        else:
            return self.cid.cid
        ###

    def nodeIDs(self,nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
           nodes = self.nodes
        if isinstance(nodes[0],int):
            #print 'if'
            return [node     for node in nodes]
        else:
            #print 'else'
            return [node.nid for node in nodes]
        ###

    def rawFields(self):
        fields = [self.type,self.lid]
        return fields

    def reprFields(self):
        return self.rawFields()

class GRAV(BaseCard):
    """
    Defines acceleration vectors for gravity or other acceleration loading
    GRAV SID CID A     N1  N2 N3    MB
    GRAV 1   3   32.2 0.0 0.0 -1.0
    """
    type = 'GRAV'
    def __init__(self,card=None,data=None):
        if card:
            ## Set identification number
            self.sid = card.field(1)
            ## Coordinate system identification number.
            self.cid = card.field(2,0)
            ## scale factor
            self.a   = card.field(3)
            ## Acceleration vector components measured in coordinate system CID
            self.N   = array(card.fields(4,7,[0.,0.,0.]))
            ## Indicates whether the CID coordinate system is defined in the main Bulk
            ## Data Section (MB = -1) or the partitioned superelement Bulk Data
            ## Section (MB = 0). Coordinate systems referenced in the main Bulk Data
            ## Section are considered stationary with respect to the assembly basic
            ## coordinate system. See Remark 10. (Integer; Default = 0)
            self.mb  = card.field(7,0)
        else:
            self.sid = data[0]
            self.cid = data[1]
            self.a   = data[2]
            self.N   = data[3:6]
            self.mb  = data[6]
            assert len(data)==7
        ###

    def crossReference(self,model):
        self.cid = model.Cid(self.cid)
    
    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def GravityVector(self):
        """returns the gravity vector in absolute coordinates"""
        p,matrix = self.cid.transformToGlobal(self.N)
        return p

    def rawFields(self):
        N = []
        for n in self.N:
            N.append(self.setBlankIfDefault(n,0.0))
        fields = ['GRAV',self.sid,self.Cid(),self.a]+N+[self.mb]
        return fields

class LSEQ(BaseCard): # Requires LOADSET in case control deck
    """
    Defines a sequence of static load sets
    @todo how does this work...
    """
    type = 'LSEQ'
    def __init__(self,card=None,data=None):
        self.sid  = card.field(1)
        self.exciteID = card.field(2)
        self.lid = card.field(3)
        self.tid = card.field(4)

    def nodeIDs(self,nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
           nodes = self.nodes
        if isinstance(nodes[0],int):
            #print 'if'
            return [node     for node in nodes]
        else:
            #print 'else'
            return [node.nid for node in nodes]
        ###

    def crossReference(self,model):
        self.lid = model.Load(self.lid)
        if self.tid:
            self.tid = model.Load(self.tid)
        ###
    
    def Lid(self):
        if isinstance(self.lid,int):
            return self.lid
        return self.lid.lid
        
    def Tid(self):
        if self.tid is None:
            return None
        if isinstance(self.tid,int):
            return self.tid
        return self.tid.tid

    def rawFields(self):
        fields = ['LSEQ',self.sid,self.exciteID,self.Lid(),self.Tid()]
        return fields

    def reprFields(self):
        return self.rawFields()

class LOAD(Load):
    type = 'LOAD'
    def __init__(self,card=None,data=None):
        """
        @todo parse the loads data to have scale factor and load
        """
        if card:
            #fields   = card.fields()

            ## load ID
            self.lid = card.field(1)
            #self.id  = self.lid

            ## overall scale factor
            self.s  = card.field(2)

            #print "nFields = ",nFields
            #print "nLoads  = ",nLoads

            loads = card.fields(3) # temp list
            nLoadFields = len(loads)
            nLoads  = nLoadFields/2
            assert nLoadFields%2==0

            ## individual scale factors (corresponds to loadIDs)
            self.scaleFactors = []

            ## individual loadIDs (corresponds to scaleFactors)
            self.loadIDs = []

            for i in range(0,nLoadFields,2):   # alternating of scale factor & load set ID
                self.scaleFactors.append(loads[i  ])
                self.loadIDs.append(     loads[i+1])
        else:
            self.lid = data[0]
            self.s   = data[1]
            self.scaleFactors = data[2]
            self.loadIDs      = data[3]
        ###

    def crossReference(self,model):
        loadIDs2 = []
        for loadID in self.loadIDs:
            loadID2 = model.Load(loadID)
            loadIDs2.append(loadID2)
        self.loadIDs = loadIDs2
    
    def LoadID(self,ID):
        if isinstance(ID,int):
            return ID
        return ID.lid
    
    def getLoads(self):
        loads = []
        for load in self.loadIDs:
            loads += loadID.getLoads()
        ###
        return loads

    def writeCodeAster(self):
        msg = '# Loads\n'
        (forceLoads,momentLoads,forceConstraints,momentConstraints) = self.organizeLoadsForCodeAster()
        # stuff
        msg += "load_bc=AFFE_CHAR_MECA(MODELE=modmod,\n"
        msg += "                      DDL_IMPO=(_F(GROUP_MA='Lleft',\n"
        msg += "                                   DX=0.0,\n"
        msg += "                                   DY=0.0,\n"
        msg += "                                   DZ=0.0,),\n"
        msg += "                                _F(GROUP_MA='Lright',\n"
        msg += "                                   DZ=0.0,),),\n"

        mags = {}
        for node,load in sorted(forceLoads.items()):
            pass
        return msg

    def organizeLoadsForCodeAster(self):
        forceLoads  = {} # spc enforced displacement (e.g. FORCE=0)
        momentLoads = {}
        forceConstraints  = {}
        momentConstraints = {}
        allLoads = array([0.,0.,0.])
        for loadID in self.loadIDs():
            loads = loadID.getLoads()
            for load in loads:
                isLoad,node,vector = load.transformLoad()
                if isinstance(load,Force):
                    if isLoad:  #load
                        if node not in forceLoads:
                            forceLoads[node]  = vector
                        else:
                            forceLoads[node] += vector
                        ###
                    else: # constraint
                        if node not in forceLoads:
                            forceConstraints[node]  = vector
                        else:
                            forceConstraints[node] += vector
                        ###
                    ###
                else:  # Moment
                    #if 
                    if isLoad: # load
                        if node not in momentLoads:
                            momentLoads[node]  = vector
                        else:
                            momentLoads[node] += vector
                        ###
                    else: # constraint
                        if node not in momentLoads:
                            momentConstraints[node]  = vector
                        else:
                            momentConstraints[node] += vector
                        ###
                    ###
                ###
            ###
        ###
        return (forceLoads,momentLoads,forceConstraints,momentConstraints)

    def rawFields(self):
        fields = ['LOAD',self.lid,self.s]
        for scaleFactor,loadID in zip(self.scaleFactors,self.loadIDs):
            fields += [scaleFactor,self.LoadID(loadID)]
        return fields

    def reprFields(self):
        return self.rawFields()

class DLOAD(Load):
    type = 'DLOAD'
    def __init__(self,card=None,data=None):
        ## load ID
        self.lid   = card.field(1)
        self.scale = card.field(2)

        fields = card.fields(3)
        n = len(fields)//2
        if len(fields)%2==1:
            n+=1
            raise Exception('missing last magnitude on DLOAD card=%s' %(card.fields()) )

        self.sids = []
        self.mags = []
        for i in range(n):
            j = 2*i
            self.mags.append(fields[j  ])
            self.sids.append(fields[j+1])  # RLOADx,TLOADx,ACSRC
        ###

    def crossReference(self,model):
        for (i,sid) in enumerate(self.sids):
            self.sids[i] = model.Load(sid)
        ###

    def Sid(self,sid):
        if isinstance(sid,int):
            return sid
        return sid.lid

    def rawFields(self):
        fields = ['DLOAD',self.lid,self.scale]
        for (mag,sid) in zip(self.mags,self.sids):
            fields += [mag,self.Sid(sid)]
        return fields

    def reprFields(self):
        return self.rawFields()

class SLOAD(Load):
    type = 'SLOAD'
    def __init__(self,card=None,data=None):
        ## load ID
        self.lid = card.field(1)
        
        fields = card.fields(2)
        n = len(fields)/2
        if len(fields)%2==1:
            n+=1
            raise Exception('missing last magnitude on SLOAD card=%s' %(card.fields()) )

        self.sids = []
        self.mags = []
        for i in range(n):
            j = 2*i
            self.sids.append(fields[j  ])
            self.mags.append(fields[j+1])
        ###

    def crossReference(self,model):
        for (i,sid) in enumerate(self.sids):
            self.sids[i] = model.Load(sid)
        ###

    def Sid(self,sid):
        if isinstance(sid,int):
            return sid
        return sid.lid

    def rawFields(self):
        fields = ['SLOAD',self.lid]
        for sid,mag in zip(self.sids,self.mags):
            fields += [self.Sid(sid),mag]
        return fields

    def reprFields(self):
        return self.rawFields()

class OneDeeLoad(Load): # FORCE/MOMENT
    type = '1D_Load'
    def __init__(self,card,data):
        Load.__init__(self,card,data)

    def getLoads(self):
        return [self]

    def transformLoad(self):
        xyxz = self.cid.transformToGlobal(self.xyz)
        if self.mag>0.:
            return (True,self.node,self.mag*xyz) # load
        return (False,self.node,xyz) # enforced displacement

    def normalize(self):
        """
        adjust the vector to a unit length
        scale up the magnitude of the vector
        """
        if self.mag != 0.0:  # enforced displacement
            normXYZ = norm(self.xyz)
            mag = self.mag*normXYZ
            self.mag *= normXYZ
            self.xyz = self.xyz/normXYZ

class Force(OneDeeLoad):
    type = '1D_Load'
    def __init__(self,card,data):
        OneDeeLoad.__init__(self,card,data)

    def F(self):
        return self.xyz*self.mag

class Moment(OneDeeLoad):
    type = 'Moment'
    def __init__(self,card,data):
        OneDeeLoad.__init__(self,card,data)

    def M(self):
        return self.xyz*self.mag

class FORCE(Force):
    type = 'FORCE'
    def __init__(self,card=None,data=None):
        """
        FORCE          3       1            100.      0.      0.      1.
        """
        Force.__init__(self,card,data)
        if card:
            self.lid  = card.field(1)
            self.node = card.field(2)
            self.cid  = card.field(3,0)
            self.mag  = card.field(4)
            xyz = card.fields(5,8,[0.,0.,0.])
        else:
            self.lid  = data[0]
            self.node = data[1]
            self.cid  = data[2]
            self.mag  = data[3]
            xyz  = data[4:7]

        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def F(self):
        return {self.node:  self.mag*self.xyz}

    #def nodeID(self):
    #    return self.node

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def rawFields(self):
        fields = ['FORCE',self.lid,self.node,self.Cid(),self.mag] + list(self.xyz)
        return fields

    def reprFields(self):
        cid = self.setBlankIfDefault(self.Cid(),0)
        fields = ['FORCE',self.lid,self.node,cid,self.mag] + list(self.xyz)
        return fields

class FORCE1(Force):
    """
    Defines a static concentrated force at a grid point by specification of a magnitude and
    two grid points that determine the direction.
    """
    type = 'FORCE1'
    def __init__(self,card=None,data=None):
        Force.__init__(self,card,data)
        if card:
            self.lid  = card.field(1)
            self.node = card.field(2)
            self.mag  = card.field(3)
            self.g1   = card.field(4)
            self.g2   = card.field(5)
        else:
            self.lid  = data[0]
            self.node = data[1]
            self.mag  = data[2]
            self.g1   = data[3]
            self.g2   = data[4]
        ###

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        self.node = mdodel.Node(self.node)
        self.g1 = model.Node(self.g1)
        self.g2 = model.Node(self.g2)
        self.xyz = self.g2.Position()-self.g1.Position()
        #self.xyz = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        self.Normalize()
    
    def G1(self):
        if isinstance(self.g1,int) or isinstance(self.g1,float):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2,int) or isinstance(self.g1,float):
            return self.g2
        return self.g2.nid

    def NodeID(self):
        if isinstance(self.node,int):
            return self.node
        return self.node.nid

    def rawFields(self):
        (node,g1,g2) = self.nodeIDs([self.node,self.G1(),self.G2()])
        fields = ['FORCE1',self.lid,node,self.mag,g1,g2]
        return fields

    def reprFields(self):
        return self.rawFields()

class FORCE2(Force):
    """
    Defines a static concentrated force at a grid point by specification of a magnitude and
    four grid points that determine the direction.
    """
    type = 'FORCE2'
    def __init__(self,card=None,data=None):
        """
        FORCE2 SID G F G1 G2 G3 G4
        """
        Force.__init__(self,card,data)
        if card:
            self.lid  = card.field(1)
            self.node = card.field(2)
            self.mag  = card.field(3)
            self.g1   = card.field(4)
            self.g2   = card.field(5)
            self.g3   = card.field(5)
            self.g4   = card.field(5)
        else:
            self.lid  = data[0]
            self.node = data[1]
            self.mag  = data[2]
            self.g1   = data[3]
            self.g2   = data[4]
            self.g3   = data[5]
            self.g4   = data[6]
        ###

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        self.node = mdodel.Node(self.node)

        v12 = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        v34 = model.Node(self.g4).Position() - model.Node(self.g3).Position()
        v12 = v12/norm(v12)
        v34 = v34/norm(v34)
        self.xyz = cross(v12,v34)
        self.Normalize()

    def NodeID(self):
        if isinstance(self.node,int):
            return self.node
        return self.node.nid

    def rawFields(self):
        (node,g1,g2,g3,g4) = self.nodeIDs([self.node,self.g1,self.g2,self.g3,self.g4])
        fields = ['FORCE2',self.lid,node,self.mag,g1,g2,g3,g4]
        return fields

    def reprFields(self):
        return self.rawFields()

class MOMENT(Moment):  # can i copy the force init without making the MOMENT a FORCE ???
    type = 'MOMENT'
    def __init__(self,card=None,data=None):
        """
        Defines a static concentrated moment at a grid point by specifying a scale factor and
        a vector that determines the direction.

        MOMENT SID G CID M    N1  N2  N3
        MOMENT 2   5   6 2.9 0.0 1.0 0.0
        """
        Moment.__init__(self,card,data)
        self.lid  = card.field(1)
        self.node = card.field(2)
        self.cid  = card.field(3,0)
        self.mag  = card.field(4)

        xyz = card.fields(5,8,[0.,0.,0.])
        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def rawFields(self):
        fields = ['MOMENT',self.lid,self.node,self.Cid(),self.mag] + list(self.xyz)
        return fields

    def reprFields(self):
        cid = self.setBlankIfDefault(self.Cid(),0)
        fields = ['MOMENT',self.lid,self.node,cid,self.mag] + list(self.xyz)
        return fields

class MOMENT1(Moment):
    type = 'MOMENT1'
    def __init__(self,card=None,data=None):
        """
        Defines a static concentrated moment at a grid point by specifying a magnitude and
        two grid points that determine the direction

        MOMENT1 SID G M G1 G2
        """
        Moment.__init__(self,card,data)
        if card:
            self.lid  = card.field(1)
            self.node = card.field(2)
            self.mag  = card.field(3)
            self.g1   = card.field(4)
            self.g2   = card.field(5)
            self.g3   = card.field(6)
            self.g4   = card.field(7)
            xyz = card.fields(5,8,[0.,0.,0.])
        else:
            self.lid  = data[0]
            self.node = data[1]
            self.mag  = data[2]
            self.g1   = data[3]
            self.g2   = data[4]
            self.g3   = data[5]
            self.g4   = data[6]
            xyz       = data[7:10]
        ###

        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        self.node = mdodel.Node(self.node)
        self.xyz = model.Node(self.g2).Position() - model.Node(self.g1).Position()
        self.Normalize()

    def rawFields(self):
        (node,g1,g2) = self.nodeIDs([self.node,self.g1,self.g2])
        fields = ['MOMENT1',self.lid,node,self.mag,g1,g2]
        return fields

    def reprFields(self):
        return self.rawFields()


class MOMENT2(Moment):
    type = 'MOMENT2'
    def __init__(self,card=None,data=None):
        """
        Defines a static concentrated moment at a grid point by specification of a magnitude
        and four grid points that determine the direction.

        MOMENT2 SID G M G1 G2 G3 G4
        """
        Moment.__init__(self,card,data)
        if card:
            self.lid  = card.field(1)
            self.node = card.field(2)
            self.mag  = card.field(3)
            self.g1   = card.field(4)
            self.g2   = card.field(5)
            self.g3   = card.field(6)
            self.g4   = card.field(7)
            xyz = card.fields(5,8,[0.,0.,0.])
        else:
            self.lid  = data[0]
            self.node = data[1]
            self.mag  = data[2]
            self.g1   = data[3]
            self.g2   = data[4]
            self.g3   = data[5]
            self.g4   = data[6]
            xyz       = data[7:10]
        ###
        assert len(xyz)==3,'xyz=%s' %(xyz)
        self.xyz = array(xyz)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        (self.g1,self.g2,self.g3,self.g4) = model.Nodes(self.g1,self.g2,self.g3,self.g4)
        v12 = g2.Position()-g1.Position()
        v34 = g4.Position()-g3.Position()
        v12 = v12/norm(v12)
        v34 = v34/norm(v34)
        self.xyz = cross(v12,v34)

    def rawFields(self):
        (node,g1,g2,g3,g4) = self.nodeIDs([self.node,self.g1,self.g2,self.g3,self.g4])
        fields = ['MOMENT2',self.lid,node,self.mag,g1,g2,g3,g4]
        return fields

    def reprFields(self):
        return self.rawFields()

class PLOAD(Load):
    type = 'PLOAD'
    def __init__(self,card=None,data=None):
        if card:
            self.lid   = card.field(1)
            self.p     = card.field(2)
            nodes      = card.fields(3,7)
            self.nodes = self.wipeEmptyFields(nodes)
        else:
            self.lid   = data[0]
            self.p     = data[1]
            self.nodes = data[2:]
            print("PLOAD = ",data)
            raise NotImplementedError('PLOAD')
        assert len(self.nodes) in [3,4],'nodes=%s' %(self.nodes)
    
    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def rawFields(self):
        fields = ['PLOAD',self.lid,self.p]+self.nodeIDs()
        return fields

    def reprFields(self):
        return self.rawFields()

class PLOAD1(Load):
    type = 'PLOAD1'
    validTypes = ['FX','FY','FZ','FXE','FYE','FZE',
                  'MX','MY','MZ','MXE','MYE','MZE']
    validScales = ['LE','FR','LEPR','FRPR']
    def __init__(self,card=None,data=None):
        if card:
            self.lid   = card.field(1)
            self.eid   = card.field(2)
            self.Type  = card.field(3)
            self.scale = card.field(4)
            self.x1    = card.field(5)
            self.p1    = card.field(6)
            self.x2    = card.field(7)
            self.p2    = card.field(8)
        else:
            self.lid   = data[0]
            self.eid   = data[1]
            self.Type  = data[2]
            self.scale = data[3]
            self.x1    = data[4]
            self.p1    = data[5]
            self.x2    = data[6]
            self.p2    = data[7]
        ###
        assert self.Type  in self.validTypes, '%s is an invalid type on the PLOAD1 card' %(self.Type)
        assert self.scale in self.validScales,'%s is an invalid scale on the PLOAD1 card' %(self.scale)

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def rawFields(self):
        fields = ['PLOAD1',self.lid,self.eid,self.Type,self.scale,self.x1,self.p1,self.x2,self.p2]
        return fields

    def reprFields(self):
        return self.rawFields()

class PLOAD2(Load):
    type = 'PLOAD2'
    def __init__(self,card=None,data=None):
        if card:
            self.lid = card.field(1)
            self.p   = card.field(2)
            eids = card.fields(3,9)

            if card.field(4)=='THRU':
                #print "PLOAD2 %s %s" %(eids[0],eids[-1])
                eids = [i for i in range(eids[0],eids[2]+1)]
                #print "found a THRU on PLOAD2"
                #raise NotImplementedError('PLOAD2')
            ###
            self.eids = eids
        else:
            self.lid   = data[0]
            self.p     = data[1]
            self.eids = list(data[2:])
            #print "PLOAD2 = ",data
        ###

    def crossReference(self,model):
        """@todo cross reference and fix repr function"""
        pass

    def rawFields(self):
        fields = ['PLOAD2',self.lid,self.p]
        if len(self.eids)>6:
            fields += [self.eids[0],'THRU',self.eids[-1]]
        else:
            fields += self.eids
        return fields

    def reprFields(self):
        return self.rawFields()

class PLOAD4(Load):
    """
    @todo needs work on g1
    """
    type = 'PLOAD4'
    def __init__(self,card=None,data=None):
        if card:
            #print "card.fields() = ",card.fields()
            self.lid = card.field(1)
            self.eid = card.field(2)
            p1 = card.field(3)
            p  = card.fields(4,7,[p1,p1,p1])
            self.p = [p1]+p

            if card.field(7)=='THRU' and card.field(8):
                #print "found a THRU on PLOAD4"
                pass
                eid2 = card.field(8)
                self.eids= self.expandThru([self.eid,'THRU',eid2])
                self.g3 = None
                self.g4 = None
            else:   # used for CPENTA, CHEXA
                self.eids = None
                self.g3 = card.field(7)
                self.g4 = card.field(8)
            ###

            ## Coordinate system identification number. See Remark 2. (Integer >= 0;Default=0)
            self.cid     = card.field(9,0)
            #print "PLOAD4 cid = ",self.cid
            self.NVector = card.fields(10,13,[0.,0.,0.])
            self.sorl    = card.field(13,'SURF')
        else:
            #print "PLOAD4 = ",data
            self.lid     = data[0]
            self.eid     = data[1]
            self.p       = data[2]

            self.g1      = data[3]
            self.g34     = data[4]
            self.cid     = data[5]
            self.NVector = data[6]

            self.sorl    = data[7]
            #assert len(data)==8
            
            self.g3 = self.g1
            self.g4 = self.g34
            self.eids = []
        ###

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def crossReference(self,model):
        self.cid = model.Coord(self.cid)
        if self.g1: self.g1 = model.Node(self.g1)
        if self.g3: self.g3 = model.Node(self.g3)
        if self.g4: self.g4 = model.Node(self.g4)
        if self.eids:
            self.eids = model.Elements(self.eids)

    def Eid(self,element):
        if isinstance(element,int):
            return element
        return element.eid

    def getElementIDs(self,eid=None):
        if eid:
            return self.Eid(eid)
        eids = []
        for element in self.eids:
            eids.append(self.Eid(element))
        return eids
            
    def rawFields(self):
        cid  = self.setBlankIfDefault(self.Cid(),0)
        sorl = self.setBlankIfDefault(self.sorl,'SURF')
        p1   = self.p[0]
        p2   = self.setBlankIfDefault(self.p[1],p1)
        p3   = self.setBlankIfDefault(self.p[2],p1)
        p4   = self.setBlankIfDefault(self.p[3],p1)
        fields = ['PLOAD4',self.lid,self.eid,self.p[0],p2,p3,p4]

        #print "g3=|%s| g4=%s eids=|%s|" %(self.g3,self.g4,self.eids)
        if self.g3 is not None:
            (g3,g4) = self.nodeIDs([self.g3,self.g4])
            fields.append(g3)
            fields.append(g4)
        else:
            #print "eids = %s" %(self.eids)
            if not self.eids==None:
                try:
                    fields.append('THRU')
                    eid = self.eids[-1]
                except:
                    print("g3 = ",self.g3)
                    print("g4 = ",self.g4)
                    print("self.eids = ",self.eids)
                    raise
                ###
                fields.append(self.getElementIDs(eid) )
            else:
                fields += [None,None]
            ###
        fields.append(cid)
        
        n1 = self.setBlankIfDefault(self.NVector[0],0.0)
        n2 = self.setBlankIfDefault(self.NVector[1],0.0)
        n3 = self.setBlankIfDefault(self.NVector[2],0.0)
        fields += [n1,n2,n3]
        fields.append(sorl)
        #print "fields = ",fields
        return fields

    def reprFields(self):
        return self.rawFields()
