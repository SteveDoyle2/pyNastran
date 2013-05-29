# pylint: disable=E1101,C0103,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.cards.nodes import SPOINT

class DeprecatedGetMethods(object):

    def Grav(self, sid, msg):
        raise DeprecationWarning('use Load(sid) instead of Grav(sid)')
        return self.Load(sid, msg)

    def Flfact(self, sid, msg):
        return self.FLFACT(sid, msg)


class GetMethods(DeprecatedGetMethods):
    def __init__(self):
        pass

    #--------------------
    # NODE CARDS

    def nNodes(self):
        return len(self.nodes.nids)

    def nodeIDs(self):
        return self.nodes.nids.keys()

    def getNodeIDsWithElement(self, eid):
        return self.getNodeIDsWithElements([eid])

    def getNodeIDsWithElements(self, eids, msg=''):
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid, msg=msg)
            self.log.debug("element.pid = %s" % (element.pid))
            nids = set(element.nids)
            nids2 = nids2.union(nids)
        return nids2

    def Node(self, nid, allowEmptyNodes=False, msg=''):
        try:
            if (nid == 0 or nid is None) and allowEmptyNodes:
                return None
            elif nid in self.nodes.nid:
                return nid
            elif self.spoints and nid in self.spoints.spoints:
                return SPOINT(nid)
            else:
                raise KeyError('nid=%s is not a GRID or SPOINT%s' % (nid, msg))
        except:
            if self._xref == 2: # 'partial'
                return nid
            raise

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

    def getElementIDsWithPID(self, pid):
        return self.getElementIDsWithPIDs([pid])

    def getElementIDsWithPIDs(self, pids, msg=''):
        eids = self.elementIDs()
        eids2 = []
        #print "eids = ",eids
        for eid in eids:
            element = self.Element(eid, msg)
            if element.Pid() in pids:
                eids2.append(eid)
        return (eids2)

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

    def Element(self, eid, msg):
        try:
            return self.elements[eid]
        except KeyError:
            if self._xref == 2:
                return eid
            raise KeyError('eid=%s not found%s.  Allowed elements=%s'
                           % (eid, msg, self.elements.keys()))

    def Elements(self, eids, msg):
        elements = []
        for eid in eids:
            elements.append(self.Element(eid, msg))
        return elements

    def RigidElement(self, eid, msg):
        try:
            return self.rigidElements[eid]
        except KeyError:
            if self._xref == 2:
                return eid
            raise KeyError('eid=%s not found%s.  Allowed rigidElements=%s'
                           % (pid, msg, self.rigidElements.keys()))

    #--------------------
    # PROPERTY CARDS

    def propertyIDs(self):
        return self.properties.keys()

    def Property(self, pid, msg):
        try:
            return self.properties[pid]
        except KeyError:
            if self._xref == 2:
                return pid
            raise KeyError('pid=%s not found%s.  Allowed Pids=%s'
                           % (pid, msg, self.propertyIDs()))

    def Properties(self, pids, msg):
        properties = []
        #print "pids = ",pids
        for pid in pids:
            properties.append(self.Property(pid, msg))
        return properties

    def Phbdy(self, pid, msg):
        try:
            return self.phbdys[pid]
        except KeyError:
            if self._xref == 2:
                return pid
            raise KeyError('pid=%s not found%s.  Allowed PHBDY Pids=%s'
                           % (pid, msg, self.phbdys.keys()))

    #--------------------
    # MATERIAL CARDS

    def structuralMaterialIDs(self):
        return self.materials.keys()

    def materialIDs(self):
        return self.materials.keys() + self.thermalMaterials.keys()

    def thermalMaterialIDs(self):
        return self.thermalMaterials.keys()

    def Material(self, mid, msg):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermalMaterials:
            return self.thermalMaterials[mid]
        else:
            if self._xref == 2:
                return mid
            msg = '\n' + msg
            raise KeyError('Invalid Material ID:  mid=%s%s' % (mid, msg))

    def StructuralMaterial(self, mid, msg):
        try:
            mat = self.materials[mid]
        except KeyError:            
            if self._xref == 2:
                return mid
            msg = '\n' + msg
            raise KeyError('Invalid Material ID:  mid=%s%s' % (mid, msg))

    def ThermalMaterial(self, mid, msg):
        try:
            mat = self.thermalMaterials[mid]
        except KeyError:            
            if self._xref == 2:
                return mid
            msg = '\n' + msg
            raise KeyError('Invalid Material ID:  mid=%s%s' % (mid, msg))

    def Materials(self, mids, msg):
        materials = []
        for mid in mids:
            materials.append(self.Material(mid, msg))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid, msg):
        #print 'sid=%s self.loads=%s' %(sid,(self.loads.keys()))
        assert isinstance(sid, int), 'sid=%s is not an integer\n' % sid
        if sid in self.loads:
            load = self.loads[sid]
        else:
            if self._xref == 2:
                return sid
            raise KeyError('cannot find LoadID=|%s|%s.' % (sid, msg))
        return load

    #--------------------
    # SPCs

    def SPC(self, conid, msg):
        #print 'conid=%s self.spcs=%s' %(conid,(self.spcs.keys()))
        assert isinstance(conid, int), 'conid=%s is not an integer\n' % conid
        if conid in self.spcs:
            constraint = self.spcs[conid]
        else:
            if self._xref == 2:
                return conid
            raise KeyError('cannot find ConstraintID=|%s|.' % conid)
        return constraint

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid, msg):
        try:
            return self.coords[cid]
        except KeyError:
            if self._xref == 2:
                return cid
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coordIDs()))
    
    def coordIDs(self):
        return sorted(self.coords.keys())
    #--------------------
    # AERO CARDS

    def nCAeros(self):
        return len(self.caeros)

    def Aero(self, acsid, msg):
        try:
            return self.aero[acsid]
        except KeyError:
            if self._xref == 2:
                return acsid
            raise KeyError('acsid=%s not found%s.  Allowed AERO=%s'
                           % (acsid, msg, self.aero.keys()))

    def Aeros(self, acsid, msg):
        try:
            return self.aeros[acsid]
        except KeyError:
            if self._xref == 2:
                return acsid
            raise KeyError('acsid=%s not found%s.  Allowed AEROS=%s'
                           % (acsid, msg, self.aeros.keys()))

    def Spline(self, eid, msg):
        try:
            return self.splines[eid]
        except KeyError:
            if self._xref == 2:
                return eid
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, self.splines.keys()))

    def CAero(self, eid, msg):
        try:
            return self.caeros[eid]
        except KeyError:
            if self._xref == 2:
                return eid
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, self.caeros.keys()))

    def PAero(self, pid, msg):
        try:
            return self.paeros[pid]
        except KeyError:
            if self._xref == 2:
                return pid
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, self.paeros.keys()))

    def Gust(self, sid, msg):
        try:
            return self.gusts[sid]
        except KeyError:
            if self._xref == 2:
                return sid
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, self.gusts.keys()))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid, msg):
        try:
            return self.aestats[aid]
        except KeyError:
            if self._xref == 2:
                return aid
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, self.aestats.keys()))

    def AELink(self, linkID, msg):
        try:
            return self.aelinks[linkID]
        except KeyError:
            if self._xref == 2:
                return linkID
            raise KeyError('linkID=%s not found%s.  Allowed AELINKs=%s'
                           % (linkID, msg, self.aelinks.keys()))

    def AEParam(self, aid, msg):
        try:
            return self.aeparams[aid]
        except KeyError:
            if self._xref == 2:
                return aid
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, self.aeparams.keys()))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid, msg):
        try:
            return self.flfacts[sid]
        except KeyError:
            if self._xref == 2:
                return sid
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, self.flfacts.keys()))

    def Flutter(self, fid, msg):
        try:
            return self.flutters[fid]
        except KeyError:
            if self._xref == 2:
                return fid
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, self.flutters.keys()))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid, msg):
        try:
            return self.dconstrs[oid]
        except KeyError:
            if self._xref == 2:
                return oid
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, self.dconstrs.keys()))

    def Desvar(self, oid, msg):
        try:
            return self.desvars[oid]
        except KeyError:
            if self._xref == 2:
                return oid
            raise KeyError('oid=%s not found%s.  Allowed DESVARs=%s'
                           % (oid, msg, self.desvars.keys()))

    def DDVal(self, oid, msg):
        try:
            return self.ddvals[oid]
        except KeyError:
            if self._xref == 2:
                return oid
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, self.ddvals.keys()))

    #--------------------
    # SET CARDS

    def Set(self, sid, msg):
        try:
            return self.sets[sid]
        except KeyError:
            if self._xref == 2:
                return sid
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, self.sets.keys()))

    def SetSuper(self, seid, msg):
        try:
            return self.setsSuper[seid]
        except KeyError:
            if self._xref == 2:
                return seid
            raise KeyError('seid=%s not found%s.  Allowed SETx=%s'
                           % (seid, msg, self.setsSuper.keys()))

    #--------------------
    # METHOD CARDS
    def Method(self, sid, msg):
        try:
            return self.methods[sid]
        except KeyError:
            if self._xref == 2:
                return sid
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, self.methods.keys()))

    def CMethod(self, sid, msg):
        try:
            return self.cmethods[sid]
        except KeyError:
            if self._xref == 2:
                return sid
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, self.cmethods.keys()))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg):
        try:
            return self.tables[tid]
        except KeyError:
            if self._xref == 2:
                return tid
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, self.tables.keys()))

    def RandomTable(self, tid, msg):
        try:
            return self.randomTables[tid]
        except KeyError:
            if self._xref == 2:
                return tid
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, self.randomTables.keys()))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid, msg):
        try:
            return self.nlparms[nid]
        except KeyError:
            if self._xref == 2:
                return nid
            raise KeyError('nid=%s not found%s.  Allowed NLPARMs=%s'
                           % (nid, msg, self.nlparms.keys()))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname, msg):
        try:
            return self.dmigs[dname]
        except KeyError:
            if self._xref == 2:
                return dname
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, self.dmigs.keys()))