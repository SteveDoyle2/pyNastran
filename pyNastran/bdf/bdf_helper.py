import sys
import copy
from pyNastran.bdf.fieldWriter import printCard

from bdf_writeMesh import writeMesh
from bdf_cardMethods import cardMethods
from crossReference import XrefMesh

class getMethods(object):
    def __init__(self):
        pass

    def getNodeIDs(self):
        raise Exception('use self.nodeIDs()')
        return sorted(self.nodes.keys())

    def nodeIDs(self):
        return self.nodes.keys()
        
    def elementIDs(self):
        return self.elements.keys()

    def propertyIDs(self):
        return self.properties.keys()

    def materialIDs(self):
        return self.materials.keys()

    def getElementIDs(self):
        raise Exception('use self.elementIDs() and sort it...')
        return sorted(self.elements.keys())

    def getPropertyIDs(self):
        raise Exception('use self.propertyIDs() and sort it...')
        return sorted(self.properties.keys())

    def getMaterialIDs(self):
        raise Exception('use self.materialIDs() and sort it...')
        return sorted(self.materials.keys())

    def getNodes(self):
        nodes = []
        for nid,node in sorted(self.nodes.items()):
            nodes.append(node)
        return nodes

    def getNodeIDsWithElement(self,eid):
        return self.getNodeIDsWithElements([eid])

    def getNodeIDsWithElements(self,eids):
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid)
            self.log().info("element.pid = %s" %(element.pid))
            nids = set(element.nids)
            nids2 = nids2.union(nids)
        ###
        return nids2

    def getElementIDsWithPID(self,pid):
        return self.getElementIDsWithPIDs([pid])

    def getElementIDsWithPIDs(self,pids):
        #self.log().info("pids = %s" %(pids))

        eids = self.elementIDs()
        eids2 = []
        #print "eids = ",eids
        for eid in eids:
            element = self.Element(eid)
            if element.Pid() in pids:
                eids2.append(eid)
            ###
            #print ""
        ###
        return (eids2)

    def Node(self,nid):
        return self.nodes[nid]

    def Nodes(self,nids):
        """
        Returns a series of node objects given a list of node IDs
        """
        #print "nids",nids
        nodes = []
        for nid in nids:
            nodes.append(self.nodes[nid])
        return nodes

    def Element(self,eid):
        return self.elements[eid]

    def Elements(self,eids):
        elements = []
        for eid in eids:
            elements.append(self.elements[eid])
        return elements

    def Property(self,pid):
        return self.properties[pid]

    def Properties(self,pids):
        properties = []
        for pid in pids:
            properties.append(self.properties[pid])
        return properties

    def Phbdy(self,pid):
        return self.phbdys[pid]

    def Material(self,mid):
        return self.materials[mid]

    def Materials(self,mids):
        materials = []
        for mid in mids:
            materials.append(self.materials[mid])
        return materials

    def Load(self,lid):
        return self.loads[lid]

    def Coord(self,cid):
        return self.coords[cid]


    def Aero(self,acsid):
        return self.aeros[acsid]
    def Flfact(self,sid):
        return self.flfacts[sid]
    def Flutter(self,fid):
        return self.flutters[fid]
    def Gust(self,sid):
        return self.gusts[sid]
    def CAero(self,eid):
        return self.caeros[eid]
    def Aero(self,acsid):
        return self.aeros[acsid]
    def Spline(self,eid):
        return self.splines[eid]


    def DConstr(self,oid):
        return self.dconstrs[oid]
    def Desvar(self,oid):
        return self.desvars[oid]
    def DDVal(self,oid):
        return self.ddvals[oid]

    def NLParml(self,nid):
        return self.nlparms[nid]

    def sumForces(self):
        for key,loadCase in self.loads.items():
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,FORCE):
                    f = load.mag*load.xyz
                    print "f = ",f
                    F += f
                ###
            self.log().info("case=%s F=%s\n\n" %(key,F))
        ###

    def sumMoments(self):
        p = array([0.,0.5,0.])
        for key,loadCase in self.loads.items():
            M = array([0.,0.,0.])
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,FORCE):
                    f = load.mag*load.xyz
                    node = self.Node(load.node)
                    #print "node = ",node
                    r = node.Position() - p
                    m = cross(r,f)
                    #print "m    = ",m
                    M += m
                    F += f
                #elif isinstance(load,MOMENT):
                #    m = load.mag*load.xyz
                #    M += m
                ###
            print "case=%s F=%s M=%s\n\n" %(key,F,M)
        ###

class addMethods(object):
    def __init__(self):
        pass

    def addParam(self,param):
        assert param.key not in self.params
        self.params[param.key] = param

    def addNode(self,node):
        #print node
        assert node.nid not in self.nodes
        assert node.nid>0
        self.nodes[node.nid] = node

    def addElement(self,elem):
        assert elem.eid not in self.elements,'eid=%s oldElement=\n%s' %(elem.eid,self.elements[elem.eid])
        assert elem.eid>0
        self.elements[elem.eid] = elem

    def addThermalElement(self,elem):  # same function at the moment...
        self.addElement(elem)
        #assert elem.eid not in self.elements
        #assert elem.eid>0
        #self.elements[elem.eid] = elem

    def addProperty(self,prop):
        assert prop.pid not in self.properties,'eid=%s oldElement=\n%s' %(prop.pid,self.properties[prop.pid])
        assert prop.pid>0
        self.properties[prop.pid] = prop

    def addMaterial(self,material):
        assert material.mid not in self.materials
        assert material.mid>0
        self.materials[material.mid] = material

    def addCoord(self,coord):
        assert coord.cid not in self.coords
        assert coord.cid>-1
        self.coords[coord.cid] = coord

    def addLoad(self,load):
        key = load.lid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addPHBDY(self,prop):
        assert prop.pid>0
        assert prop.pid not in self.phbdys
        self.phbdys[prop.pid] = prop

    def addConvectionProperty(self,prop):
        assert prop.pconid>0
        assert prop.pconid not in self.convectionProperties
        self.convectionProperties[prop.pconid] = prop

    #def addThermalProperty(self,prop):
    #    assert prop.pconid not in self.thermalProperties
    #    self.thermalProperties[prop.pconid] = prop

    def addThermalBC(self,bc,key):
        assert key>0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]

    def addThermalLoad(self,load):  # same function at the moment...
        key = load.sid
        assert key>0
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addConstraint_MPCADD(self,constraint):
        self.mpcObject.add(constraint)

    def addConstraint_MPC(self,constraint):
        self.mpcObject.append(constraint)

    def addConstraint_SPCADD(self,constraint):
        self.spcObject.add(constraint)

    def addConstraint_SPC(self,constraint):
        self.spcObject.append(constraint)
        #key = constraint.cid
        #if self.constraints.has_key(key):
        #    self.constraints[key].append(constraint)
        #else:
        #    self.constraints[key] = [constraint]

    def addConstraint(self,constraint):
        #self.spcObject.append(constraint)
        key = constraint.cid
        assert key>0
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def addSUPORT(self,suport):
        self.suports.append(suport)


    def addDArea(self,darea):
        assert darea.sid not in self.dareas,'eid=%s oldElement=\n%s' %(darea.sid,self.dareas[darea.sid])
        assert darea.sid>0
        self.dareas[darea.sid] = darea

    def addAero(self,aero):
        assert aero.acsid not in self.aeros
        assert aero.acsid>0
        self.aeros[aero.acsid] = aero

    def addCAero(self,caero):
        assert caero.eid not in self.caeros,'self.caeros=|%s| caero.eid=|%s|' %(self.caeros,caero.eid)
        assert caero.eid>0
        self.caeros[caero.eid] = caero

    def addGrav(self,grav):
        assert grav.sid not in self.gravs
        assert grav.sid>0
        self.gravs[grav.sid] = grav

    def addSpline(self,spline):
        assert spline.eid not in self.splines
        assert spline.eid>0
        self.splines[spline.eid] = spline

    def addGust(self,gust):
        assert gust.sid not in self.gusts
        assert gust.sid>0
        self.gusts[gust.sid] = gust

    def addFlutter(self,flutter):
        assert flutter.sid not in self.flutters
        assert flutter.sid>0
        self.flutters[flutter.sid] = flutter

    def addFLFACT(self,flfact):
        #assert flfact.sid not in self.flfacts
        self.flfacts[flfact.sid] = flfact # set id...
        #print "added flfact...flflact =\n"+flfact


    def addDConstr(self,dconstr):
        assert dconstr.oid not in self.dconstrs
        assert dconstrs.oid>0
        self.dconstrss[dconstrs.oid] = dconstrs
    def addDesvar(self,desvar):
        assert desvar.oid not in self.desvars
        assert desvar.oid>0
        self.desvars[desvar.oid] = desvar
    def addDDVal(self,ddval):
        assert ddval.oid not in self.ddvals
        assert ddval.oid>0
        self.ddvals[ddval.oid] = ddval


    def addNLParm(self,nlparm):
        assert nlparm.nid not in self.nlparms
        assert nlparm.nid>0
        self.nlparms[nlparm.nid] = nlparm
