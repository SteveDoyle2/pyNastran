import sys
import copy
from pyNastran.bdf.fieldWriter import printCard
#from BDF_Card import BDF_Card
from bdf_writeMesh import writeMesh
from bdf_cardMethods import cardMethods

class getMethods(object):
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
        for eid in eids:
            element = self.Element(eid)
            #print "element = ",element
            #print "element.pid = ",element.pid
            if element.pid in pids:
                eids2.append(eid)
            ###
            #print ""
        ###
        return (eids2)

    def Node(self,nid):
        return self.nodes[nid]

    def Element(self,eid):
        return self.elements[eid]

    def Property(self,pid):
        return self.properties[pid]

    def Material(self,mid):
        return self.materials[mid]

    def Load(self,lid):
        return self.loads[lid]

    def Coord(self,cid):
        return self.coords[cid]

    def FLFACT(self,fid):
        self.flfacts[param.key] = flfact

class addMethods(object):
    def addParam(self,param):
        assert param.key not in self.params
        self.params[param.key] = param

    def addNode(self,node):
        self.nodes[node.nid] = node

    def addElement(self,elem):
        #assert elem.eid not in self.elements
        self.elements[elem.eid] = elem

    def addProperty(self,prop):
        assert prop.pid not in self.properties
        self.properties[prop.pid] = prop

    def addMaterial(self,material):
        assert material.mid not in self.materials
        self.materials[material.mid] = material

    def addCoord(self,coord):
        assert coord.cid not in self.coords
        self.coords[coord.cid] = coord

    def addLoad(self,load):
        key = load.lid
        if self.loads.has_key(key):
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

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
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

    def addAero(self,aero):
        assert aero.acsid not in self.aeros
        self.aeros[aero.acsid] = aero

    def addGust(self,gust):
        assert gust.sid not in self.gusts
        self.gusts[gust.sid] = gust

    def addFlutter(self,flutter):
        assert flutter.sid not in self.flutters
        self.flutters[flutter.sid] = flutter

    def addFLFACT(self,flfact):
        assert flfact.sid not in self.flfacts
        self.flfacts[flfact.sid] = flfact # set id...
        #print "added flfact...flflact =\n"+flfact

