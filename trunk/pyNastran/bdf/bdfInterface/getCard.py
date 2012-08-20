# pylint: disable=E1101,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
#import copy

from pyNastran.bdf.cards.nodes import SPOINT


class GetMethods(object):
    def __init__(self):
        pass

    #--------------------
    # NODE CARDS

    def nNodes(self):
        return len(self.nodes)

    def nodeIDs(self):
        return self.nodes.keys()

    def getNodes(self):
        nodes = []
        for (nid, node) in sorted(self.nodes.iteritems()):
            nodes.append(node)
        return nodes

    def getNodeIDsWithElement(self, eid):
        return self.getNodeIDsWithElements([eid])

    def getNodeIDsWithElements(self, eids):
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid)
            self.log.debug("element.pid = %s" % (element.pid))
            nids = set(element.nids)
            nids2 = nids2.union(nids)
        ###
        return nids2

    def Node(self, nid, allowEmptyNodes=False):
        if nid == 0 and allowEmptyNodes:
            return None
        elif nid in self.nodes:
            return self.nodes[nid]
        elif self.spoints and nid in self.spoints.spoints:
            return SPOINT(nid)
        else:
            raise RuntimeError('nid=%s is not a GRID or SPOINT' % (nid))
        ###

    def Nodes(self, nids, allowEmptyNodes=False):
        """
        Returns a series of node objects given a list of node IDs
        """
        #print "nids",nids
        nodes = []
        for nid in nids:
            nodes.append(self.Node(nid, allowEmptyNodes))
        return nodes

    #--------------------
    # ELEMENT CARDS

    def nElements(self):
        return len(self.elements)

    def elementIDs(self):
        return self.elements.keys()

    def getElementIDsWithPID(self, pid):
        return self.getElementIDsWithPIDs([pid])

    def getElementIDsWithPIDs(self, pids):
        eids = self.elementIDs()
        eids2 = []
        #print "eids = ",eids
        for eid in eids:
            element = self.Element(eid)
            if element.Pid() in pids:
                eids2.append(eid)
            ###
        ###
        return (eids2)

    def getNodeIDToElementIDsMap(self):
        """
        Returns a dictionary that maps a node ID to a list of elemnents
        @todo
          support 0d or 1d elements
        @todo
          support elements with missing nodes (e.g. CQUAD8 with missing nodes)
        """
        nidToElementsMap = {}
        for nid in self.nodes:  # initalize the mapper
            nidToElementsMap[nid] = []

        if self.spoints:  # SPOINTs
            for nid in sorted(self.spoints.spoints):  # SPOINTs
                nidToElementsMap[nid] = []
            ###
        ###
        for (eid, element) in self.elements.iteritems():  # load the mapper
            try:
                # not supported for 0-D and 1-D elements
                nids = element.nodeIDs()
                for nid in nids:  # (e.g. CQUAD8 with missing node)
                    nidToElementsMap[nid].append(eid)
            except:
                pass
            ###
        ###
        return nidToElementsMap

    def getPropertyIDToElementIDsMap(self):
        """
        Returns a dictionary that maps a property ID to a list of elemnents
        """
        pidToEidsMap = {}
        pids = self.propertyIDs()
        for pid in pids:
            pidToEidsMap[pid] = []

        for eid in self.elementIDs():
            element = self.Element(eid)
            #print dir(element)
            if hasattr(element, 'pid'):
                pid = element.Pid()
                if pid == 0:  # CONM2
                    continue
                pidToEidsMap[pid].append(eid)
            ###
        ###
        return (pidToEidsMap)

    def getMaterialIDToPropertyIDsMap(self):
        """
        Returns a dictionary that maps a material ID to a list of properties
        @note
          all properties require an mid to be counted (except for PCOMP,
          which has multiple mids)
        """
        midToPidsMap = {}
        for mid in self.materialIDs():
            midToPidsMap[mid] = []

        for pid in self.propertyIDs():
            prop = self.Property(pid)
            if prop.type == 'PCOMP':
                mids = prop.Mids()

                for mid in mids:
                    if pid not in midToPidsMap[mid]:
                        midToPidsMap[mid].append(pid)
            else:  # PCOMP
                if hasattr(prop, 'mid') and prop.Mid() in mids:
                    if pid not in midToPidsMap[mid]:
                        midToPidsMap[mid].append(pid)
        return (midToPidsMap)

    def Element(self, eid):
        return self.elements[eid]

    def Elements(self, eids):
        elements = []
        for eid in eids:
            elements.append(self.elements[eid])
        return elements

    def RigidElement(self, eid):
        return self.rigidElements[eid]

    #--------------------
    # PROPERTY CARDS

    def propertyIDs(self):
        return self.properties.keys()

    def Property(self, pid):
        try:
            return self.properties[pid]
        except KeyError:
            raise KeyError('pid=%s not found.  Allowed Pids=%s'
                           % (pid, self.propertyIDs()))

    def Properties(self, pids):
        properties = []
        #print "pids = ",pids
        for pid in pids:
            properties.append(self.properties[pid])
        return properties

    def Phbdy(self, pid):
        return self.phbdys[pid]

    #--------------------
    # MATERIAL CARDS

    def structuralMaterialIDs(self):
        return self.materials.keys()

    def materialIDs(self):
        return self.materials.keys() + self.thermalMaterials.keys()

    def thermalMaterialIDs(self):
        return self.thermalMaterials.keys()

    def Material(self, mid):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermalMaterials:
            return self.thermalMaterials[mid]
        else:
            raise KeyError('Invalid Material ID:  mid=%s' % (mid))

    def StructuralMaterial(self, mid):
        return self.materials[mid]

    def ThermalMaterial(self, mid):
        return self.thermalMaterials[mid]

    def Materials(self, mids):
        materials = []
        for mid in mids:
            materials.append(self.Material(mid))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid):
        #print 'sid=%s self.loads=%s' %(sid,(self.loads.keys()))
        assert isinstance(sid, int), 'sid=%s is not an integer\n' % (sid)
        if sid in self.loads:
            load = self.loads[sid]
        else:
            raise KeyError('cannot find LoadID=|%s|.' % (sid))
        return load

    def Grav(self, sid):
        raise DeprecationWarning('use Load(sid) instead of Grav(sid)')
        return self.Load(sid)

    #--------------------
    # SPCs

    def SPC(self, conid):
        #print 'conid=%s self.spcs=%s' %(conid,(self.spcs.keys()))
        assert isinstance(conid, int), 'conid=%s is not an integer\n' % (conid)
        if conid in self.spcs:
            constraint = self.spcs[conid]
        else:
            raise KeyError('cannot find ConstraintID=|%s|.' % (conid))
        return constraint

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid):
        return self.coords[cid]

    #--------------------
    # AERO CARDS

    def nCAeros(self):
        return len(self.caeros)

    def Aero(self, acsid):
        return self.aero[acsid]

    def Aeros(self, acsid):
        return self.aeros[acsid]

    def Spline(self, eid):
        return self.splines[eid]

    def CAero(self, eid):
        return self.caeros[eid]

    def PAero(self, pid):
        return self.paeros[pid]

    def Gust(self, sid):
        return self.gusts[sid]

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid):
        return self.aestats[aid]

    def AELink(self, linkID):
        return self.aelinks[linkID]

    def AEParam(self, aid):
        return self.aeparams[aid]

    #--------------------
    # FLUTTER CARDS

    def Flfact(self, sid):
        return self.flfacts[sid]

    def Flutter(self, fid):
        return self.flutters[fid]

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid):
        return self.dconstrs[oid]

    def Desvar(self, oid):
        return self.desvars[oid]

    def DDVal(self, oid):
        return self.ddvals[oid]

    #--------------------
    # SET CARDS

    def Set(self, sid):
        return self.sets[sid]

    def SetSuper(self, seid):
        return self.setsSuper[seid]

    #--------------------
    # METHOD CARDS
    def Method(self, sid):
        return self.methods[sid]

    def CMethod(self, sid):
        return self.cMethods[sid]

    #--------------------
    # TABLE CARDS
    def Table(self, tid):
        return self.tables[tid]

    def RandomTable(self, tid):
        return self.randomTables[tid]

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid):
        return self.nlparms[nid]

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname):
        return self.dmig[dname]
    #--------------------
