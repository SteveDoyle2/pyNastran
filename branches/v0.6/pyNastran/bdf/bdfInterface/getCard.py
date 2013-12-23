# pylint: disable=E1101,C0103,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.cards.nodes import SPOINT

class GetMethodsDeprecated(object):

    def getElementIDsWithPID(self, pid):
        """
        Gets all the element IDs with a specific property ID

        :param pid: property ID
        :returns elementIDs: as a list

        .. deprecated:: will be removed in version 0.8

        The same functionality may be used by calling
          >>> self.getElementIDsWithPIDs([pid], mode='list')
        """
        warnings.warn('getElementIDsWithPID has been deprecated; use '
                      'getElementIDsWithPIDs', DeprecationWarning, stacklevel=2)
        return self.getElementIDsWithPIDs([pid], mode='list')

class GetMethods(GetMethodsDeprecated):
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
        return nids2

    def Node(self, nid, allowEmptyNodes=False, msg=''):
        if (nid == 0 or nid is None) and allowEmptyNodes:
            return None
        elif nid in self.nodes:
            return self.nodes[nid]
        elif self.spoints and nid in self.spoints.spoints:
            return SPOINT(nid)
        else:
            assert isinstance(nid, int), 'nid should be an integer; not %s' % type(nid)
            raise RuntimeError('nid=%s is not a GRID or SPOINT%s' % (nid, msg))

    def Nodes(self, nids, allowEmptyNodes=False, msg=''):
        """
        Returns a series of node objects given a list of node IDs
        """
        #print("nids",nids)
        nodes = []
        for nid in nids:
            #print("nid = %s" %(nid))
            nodes.append(self.Node(nid, allowEmptyNodes, msg))
        return nodes

    #--------------------
    # ELEMENT CARDS

    def nElements(self):
        return len(self.elements)

    def elementIDs(self):
        return self.elements.keys()

    def getElementIDsWithPIDs(self, pids, mode='list'):
        """
        Gets all the element IDs with a specific property ID
        .
        :param self: the BDF object
        :param pids: list of property ID
        :param mode:  'list' - returns the data as one list (default)
                      'dict' - returns the data as a dictionary of lists based on the property ID

        :returns elementIDs: as a list or dictionary of lists by property based on the mode
        """
        if mode not in ['list', 'dict']:
            msg = "mode=%r is not supported.  Use 'list' or 'dict'\n" % mode
            raise ValueError(msg)
        if mode == 'list':
            eids2 = []
            for eid, element in self.elements.iteritems():
                if element.Pid() in pids:
                    eids2.append(eid)
        else:
            eids2 = {}
            for pid in pids:
                eids2[pid] = []
            for eid, element in self.elements.iteritems():
                pid = element.Pid()
                if pid in pids:
                    eids2[pid].append(eid)
        return eids2

    def getNodeIDToElementIDsMap(self):
        """
        Returns a dictionary that maps a node ID to a list of elemnents

        .. todo:: support 0d or 1d elements
        .  todo:: support elements with missing nodes
                  (e.g. CQUAD8 with missing nodes)
        """
        nidToElementsMap = {}
        for nid in self.nodes:  # initalize the mapper
            nidToElementsMap[nid] = []

        if self.spoints:  # SPOINTs
            for nid in sorted(self.spoints.spoints):  # SPOINTs
                nidToElementsMap[nid] = []

        for (eid, element) in self.elements.iteritems():  # load the mapper
            try:
                # not supported for 0-D and 1-D elements
                nids = element.nodeIDs()
                for nid in nids:  # (e.g. CQUAD8 with missing node)
                    nidToElementsMap[nid].append(eid)
            except:
                pass
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
        return (pidToEidsMap)

    def getMaterialIDToPropertyIDsMap(self):
        """
        Returns a dictionary that maps a material ID to a list of properties

        .. note:: all properties require an mid to be counted (except for
                  PCOMP, which has multiple mids)
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

    def elementIDs(self):
        return self.elements.keys()

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

    def Property(self, pid, msg):
        try:
            return self.properties[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Pids=%s'
                           % (pid, msg, self.propertyIDs()))

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

    def Material(self, mid, msg=''):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermalMaterials:
            return self.thermalMaterials[mid]
        else:
            msg = '\n' + msg
            raise KeyError('Invalid Material ID:  mid=%s%s' % (mid, msg))

    def StructuralMaterial(self, mid, msg=''):
        return self.materials[mid]

    def ThermalMaterial(self, mid, msg=''):
        return self.thermalMaterials[mid]

    def Materials(self, mids, msg=''):
        materials = []
        for mid in mids:
            materials.append(self.Material(mid))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid, msg=''):
        #print 'sid=%s self.loads=%s' %(sid,(self.loads.keys()))
        assert isinstance(sid, int), 'sid=%s is not an integer\n' % sid
        if sid in self.loads:
            load = self.loads[sid]
        else:
            raise KeyError('cannot find LoadID=%r%s.' % (sid, msg))
        return load

    def Grav(self, sid, msg=''):
        raise DeprecationWarning('use Load(sid) instead of Grav(sid)')
        return self.Load(sid)

    #--------------------
    # SPCs

    def SPC(self, conid, msg=''):
        #print 'conid=%s self.spcs=%s' %(conid,(self.spcs.keys()))
        assert isinstance(conid, int), 'conid=%s is not an integer\n' % (conid)
        if conid in self.spcs:
            constraint = self.spcs[conid]
        else:
            raise KeyError('cannot find ConstraintID=%r.\n%s' % (conid, msg))
        return constraint

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid, msg=''):
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coordIDs()))

    def coordIDs(self):
        return sorted(self.coords.keys())
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
    def Table(self, tid, msg=''):
        return self.tables[tid]

    def RandomTable(self, tid, msg=''):
        return self.randomTables[tid]

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid):
        return self.nlparms[nid]

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname):
        return self.dmig[dname]