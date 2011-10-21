import sys
import copy
from pyNastran.bdf.fieldWriter import printCard
#from BDF_Card import BDF_Card
from bdf_writeMesh import writeMesh
from bdf_cardMethods import cardMethods

class getMethods(object):
    def getNodeIDs(self):
        return sorted(self.nodes.keys())

    def propertyIDs(self):
        return self.properties.keys()

    def nodeIDs(self):
        return self.nodes.keys()
        
    def getElementIDs(self):
        return sorted(self.elements.keys())

    def elementIDs(self):
        return self.elements.keys()

    def getPropertyIDs(self):
        return sorted(self.properties.keys())

    def getMaterialIDs(self):
        return sorted(self.materials.keys())

    def getNodes(self):
        nodes = []
        for nid,node in sorted(self.nodes.items()):
            nodes.append(node)
        return nodes

    def getNodeIDsWithElements(self,eids):
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid)
            self.log().info("element.pid = %s" %(element.pid))
            nids = set(element.nids)
            nids2 = nids2.union(nids)
        ###
        return nids2

    def getElementIDsWithPIDs(self,pids):
        self.log().info("pids = %s" %(pids))

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
    def Load(self,lid):
        return self.loads[lid]

    def Property(self,pid):
        return self.properties[pid]

    def Material(self,mid):
        return self.materials[mid]

    def Coord(self,cid):
        return self.coords[cid]

    def FLFACT(self,fid):
        self.flfacts[param.key] = flfact

class addMethods(object):
    def addFLFACT(self,flfact):
        self.flfacts[flfact.sid] = flfact # set id...
        asdf
        print "added flfact...flflact = ",flfact

    def addParam(self,param):
        self.params[param.key] = param

    def addNode(self,node):
        self.nodes[node.nid] = node

    def addElement(self,elem):
        self.elements[elem.eid] = elem

    def addProperty(self,prop):
        self.properties[prop.pid] = prop

    def addMaterial(self,material):
        self.materials[material.mid] = material

    def addCoord(self,coord):
        self.coords[coord.cid] = coord

    def addLoad(self,load):
        key = load.lid
        #print "type(self.loads) = ",type(self.loads)
        if self.loads.has_key(key):
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addConstraint(self,constraint):
        key = constraint.cid
        #print "type(self.constraints) = ",type(self.constraints)
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

