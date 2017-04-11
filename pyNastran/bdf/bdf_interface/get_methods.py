from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import numpy as np

from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.cards.nodes import SPOINT, EPOINT
from pyNastran.utils import integer_types


class GetMethods(BDFAttributes):
    def __init__(self):
        BDFAttributes.__init__(self)

    def Node(self, nid, allow_empty_nodes=False, msg=''):
        if (nid == 0 or nid is None) and allow_empty_nodes:
            return None
        elif nid in self.nodes:
            return self.nodes[nid]
        elif self.spoints and nid in self.spoints.points:
            return SPOINT(nid)
        elif self.epoints and nid in self.epoints.points:
            return EPOINT(nid)
        else:
            assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
            nid_list = np.unique(list(self.nodes.keys()))
            msg = 'nid=%s is not a GRID, SPOINT, or EPOINT%s\n' % (nid, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % list(self.spoints.points)
            if self.epoints:
                msg += 'epoints=%s\n' % list(self.epoints.points)
            raise KeyError(msg)

    def Nodes(self, nids, allow_empty_nodes=False, msg=''):
        """
        Returns a series of node objects given a list of IDs
        """
        nodes = []
        for nid in nids:
            nodes.append(self.Node(nid, allow_empty_nodes=allow_empty_nodes, msg=msg))
        return nodes

    def Point(self, nid, msg=''):
        if nid in self.points:
            return self.points[nid]
        else:
            assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
            nid_list = np.unique(list(self.points.keys()))
            raise KeyError('nid=%s is not a POINT%s\n%s' % (nid, msg, nid_list))

    def Points(self, nids, msg=''):
        """
        Returns a series of POINT objects given a list of IDs
        """
        points = []
        for nid in nids:
            points.append(self.Point(nid, msg=msg))
        return points

    def Element(self, eid, msg=''):
        """
        Gets an element

        Doesn't get rigid (RROD, RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE)
        or mass (CMASS1, CONM2))
        """
        try:
            return self.elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed elements=%s'
                           % (eid, msg, np.unique(list(self.elements.keys()))))

    def Elements(self, eids, msg=''):
        elements = []
        for eid in eids:
            elements.append(self.Element(eid, msg))
        return elements

    def Mass(self, eid, msg=''):
        """gets a mass element (CMASS1, CONM2)"""
        try:
            return self.masses[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed masses=%s'
                           % (eid, msg, np.unique(list(self.masses.keys()))))

    def RigidElement(self, eid, msg=''):
        """gets a rigid element (RBAR, RBE2, RBE3, RBAR, RBAR1, RROD, RSPLINE)"""
        try:
            return self.rigid_elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed rigid_elements=%s'
                           % (eid, msg, np.unique(list(self.rigid_elements.keys()))))

    #--------------------
    # PROPERTY CARDS

    def Property(self, pid, msg=''):
        """
        gets an elemental property (e.g. PSOLID, PLSOLID, PCOMP, PSHELL, PSHEAR);
        not mass property (PMASS)
        """
        try:
            return self.properties[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Pids=%s'
                           % (pid, msg, self.property_ids))

    def Properties(self, pids, msg=''):
        properties = []
        for pid in pids:
            properties.append(self.Property(pid, msg))
        return properties

    def PropertyMass(self, pid, msg=''):
        """
        gets a mass property (PMASS)
        """
        try:
            return self.properties_mass[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Mass Pids=%s'
                           % (pid, msg, np.unique(list(self.mass_property.keys()))))

    def Phbdy(self, pid, msg=''):
        """gets a PHBDY"""
        try:
            return self.phbdys[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PHBDY Pids=%s'
                           % (pid, msg, np.unique(list(self.phbdys.keys()))))

    #--------------------
    # MATERIAL CARDS

    def get_structural_material_ids(self):
        return self.materials.keys()

    def get_material_ids(self):
        return (self.materials.keys() + self.thermal_materials.keys() +
                self.hyperelastic_materials.keys())

    def get_thermal_material_ids(self):
        return self.thermal_materials.keys()

    def Material(self, mid, msg=''):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermal_materials:
            return self.thermal_materials[mid]
        else:
            msg = '\n' + msg
            raise KeyError('Invalid Material ID:  mid=%s%s' % (mid, msg))

    def StructuralMaterial(self, mid, msg=''):
        try:
            mat = self.materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Structural Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def ThermalMaterial(self, mid, msg=''):
        try:
            mat = self.thermal_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Thermal Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def HyperelasticMaterial(self, mid, msg=''):
        try:
            mat = self.hyperelastic_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Hyperelastic Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def Materials(self, mids, msg=''):
        materials = []
        for mid in mids:
            materials.append(self.Material(mid, msg))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid, msg=''):
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid
        if sid in self.loads:
            load = self.loads[sid]
        else:
            raise KeyError('cannot find LoadID=%r%s.\nLoadIDs=%s\n' % (
                sid, msg, np.unique(list(self.loads.keys()))))
        return load

    def DLoad(self, sid, msg=''):
        """gets a DLOAD"""
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid
        if sid in self.dloads:
            load = self.dloads[sid]
        else:
            raise KeyError('cannot find DLoadID=%r%s.\nDLoadIDs=%s\n' % (
                sid, msg, np.unique(list(self.dloads.keys()))))
        return load

    def get_dload_entries(self, sid, msg=''):
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid
        if sid in self.dload_entries:
            load = self.dload_entries[sid]
        else:
            raise KeyError('cannot find DLoad Entry ID=%r%s.\nDLoadEntryIDs=%s\n' % (
                sid, msg, np.unique(list(self.dload_entries.keys()))))
        return load

    def DELAY(self, delay_id, msg=''):
        """gets a DELAY"""
        assert isinstance(delay_id, integer_types), delay_id
        try:
            return self.delays[delay_id]
        except KeyError:
            raise KeyError('delay_id=%s not found%s.  Allowed DELAY=%s'
                           % (delay_id, msg, list(self.delays.keys())))

    def DPHASE(self, dphase_id, msg=''):
        """gets a DPHASE"""
        assert isinstance(dphase_id, integer_types), dphase_id
        try:
            return self.dphases[dphase_id]
        except KeyError:
            raise KeyError('dphase_id=%s not found%s.  Allowed DPHASE=%s'
                           % (dphase_id, msg, list(self.dphases.keys())))

    #--------------------
    def MPC(self, mpc_id, msg=''):
        """gets an MPC"""
        assert isinstance(mpc_id, integer_types), 'mpc_id=%s is not an integer\n' % mpc_id
        if mpc_id in self.mpcs:
            constraint = self.mpcs[mpc_id]
        else:
            raise KeyError('cannot find MPC ID=%r%s.\nAllowed MPCs=%s' % (
                mpc_id, msg, np.unique(list(self.mpcs.keys()))))
        return constraint

    def SPC(self, spc_id, msg=''):
        """gets an SPC"""
        assert isinstance(spc_id, integer_types), 'spc_id=%s is not an integer\n' % spc_id
        if spc_id in self.spcs:
            constraint = self.spcs[spc_id]
        else:
            raise KeyError('cannot find SPC ID=%r%s.\nAllowed SPCs=%s' % (
                spc_id, msg, np.unique(list(self.spcs.keys()))))
        return constraint

    #--------------------
    # Sets
    def SET1(self, set_id, msg=''):
        """gets a SET1"""
        assert isinstance(set_id, integer_types), 'set_id=%s is not an integer\n' % set_id
        if set_id in self.sets:
            set1 = self.sets[set_id]
        else:
            raise KeyError('cannot find SET1 ID=%r.\n%s' % (set_id, msg))
        return set1

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid, msg=''):
        """gets an COORDx"""
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coord_ids))

    #--------------------
    # AERO CARDS

    def AEList(self, aelist, msg=''):
        """gets an AELIST"""
        try:
            return self.aelists[aelist]
        except KeyError:
            raise KeyError('aelist=%s not found%s.  Allowed AELIST=%s'
                           % (aelist, msg, np.unique(list(self.aelists.keys()))))

    def AEFact(self, aefact, msg=''):
        """gets an AEFACT"""
        try:
            return self.aefacts[aefact]
        except KeyError:
            raise KeyError('aefact=%s not found%s.  Allowed AEFACT=%s'
                           % (aefact, msg, np.unique(list(self.aefacts.keys()))))

    def AESurf(self, aesurf_id, msg=''):
        """gets an AESURF"""
        try:
            return self.aesurf[aesurf_id]
        except KeyError:
            raise KeyError('aesurf=%s not found%s.  Allowed AESURF=%s'
                           % (aesurf_id, msg, np.unique(list(self.aesurf.keys()))))

    def Acsid(self, msg=''):
        """gets the aerodynamic system coordinate"""
        if self.aero is not None:
            acsid_aero = self.aero.Acsid()
        if self.aeros is not None:
            acsid_aeros = self.aeros.Acsid()

        if self.aero is not None:
            if self.aeros is not None:
                assert acsid_aero == acsid_aeros, 'AERO acsid=%s, AEROS acsid=%s' % (acsid_aero,
                                                                                     acsid_aeros)
            cid = self.Coord(acsid_aero, msg=msg)
        elif self.aeros is not None:
            cid = self.Coord(acsid_aeros, msg=msg)
        else:
            msg = 'neither AERO nor AEROS cards exist.'
            raise RuntimeError(msg)
        return cid

    def Aero(self, msg=''):
        """gets the AERO"""
        if self.aero is not None:
            return self.aero
        else:
            raise RuntimeError('no AERO card found%s.' % (msg))

    def Aeros(self, msg=''):
        """gets the AEROS"""
        if self.aeros is not None:
            return self.aeros
        else:
            raise RuntimeError('no AEROS card found%s.' % (msg))

    def Spline(self, eid, msg=''):
        """gets a SPLINEx"""
        try:
            return self.splines[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, np.unique(list(self.splines.keys()))))

    def CAero(self, eid, msg=''):
        """gets an CAEROx"""
        try:
            return self.caeros[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, np.unique(list(self.caero_ids))))

    def PAero(self, pid, msg=''):
        """gets a PAEROx"""
        try:
            return self.paeros[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, np.unique(list(self.paeros.keys()))))

    def Gust(self, sid, msg=''):
        """gets a GUST"""
        try:
            return self.gusts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, np.unique(list(self.gusts.keys()))))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid, msg=''):
        """gets an AESTAT"""
        try:
            return self.aestats[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, np.unique(list(self.aestats.keys()))))

    def AELIST(self, aid, msg=''):
        """gets an AELIST"""
        try:
            return self.aelists[aid]
        except KeyError:
            raise KeyError('id=%s not found%s.  Allowed AELISTs=%s'
                           % (aid, msg, np.unique(list(self.aelists.keys()))))

    def AELink(self, link_id, msg=''):
        """gets an AELINK"""
        try:
            return self.aelinks[link_id]
        except KeyError:
            raise KeyError('link_id=%s not found%s.  Allowed AELINKs=%s'
                           % (link_id, msg, np.unique(list(self.aelinks.keys()))))

    def AEParam(self, aid, msg=''):
        """gets an AEPARM"""
        try:
            return self.aeparams[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, np.unique(list(self.aeparams.keys()))))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid, msg=''):
        """gets an FLFACT"""
        try:
            return self.flfacts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, np.unique(list(self.flfacts.keys()))))

    def Flutter(self, fid, msg=''):
        """gets a FLUTTER"""
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, np.unique(list(self.flutters.keys()))))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid, msg=''):
        """gets a DCONSTR"""
        try:
            return self.dconstrs[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, np.unique(list(self.dconstrs.keys()))))

    def DResp(self, dresp_id, msg=''):
        """gets a DRESPx"""
        try:
            return self.dresps[dresp_id]
        except KeyError:
            raise KeyError('dresp_id=%s not found%s.  Allowed DRESPx=%s'
                           % (dresp_id, msg, np.unique(list(self.dresps.keys()))))

    def Desvar(self, desvar_id, msg=''):
        """gets a DESVAR"""
        try:
            return self.desvars[desvar_id]
        except KeyError:
            raise KeyError('desvar_id=%s not found%s.  Allowed DESVARs=%s'
                           % (desvar_id, msg, np.unique(list(self.desvars.keys()))))

    def DDVal(self, oid, msg=''):
        """gets a DDVAL"""
        try:
            return self.ddvals[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, np.unique(list(self.ddvals.keys()))))

    def DVcrel(self, dv_id, msg=''):
        """gets a DVCREL1/DVCREL2"""
        try:
            return self.dvcrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVCRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvcrels.keys()))))

    def DVmrel(self, dv_id, msg=''):
        """gets a DVMREL1/DVMREL2"""
        try:
            return self.dvmrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVMRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvmrels.keys()))))

    def DVprel(self, dv_id, msg=''):
        """gets a DVPREL1/DVPREL2"""
        try:
            return self.dvprels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVPRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvprels.keys()))))

    #--------------------
    # SET CARDS

    def Set(self, sid, msg=''):
        try:
            return self.sets[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, np.unique(list(self.sets.keys()))))

    #--------------------
    # METHOD CARDS
    def Method(self, sid, msg=''):
        """gets a METHOD (EIGR, EIGRL)"""
        try:
            return self.methods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, np.unique(list(self.methods.keys()))))

    def CMethod(self, sid, msg=''):
        """gets a METHOD (EIGC)"""
        try:
            return self.cmethods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, np.unique(list(self.cmethods.keys()))))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg=''):
        """gets a TABLES1, TABLEST, ???"""
        try:
            return self.tables[tid]
        except KeyError:
            table_keys = np.unique(list(self.tables.keys()))
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, table_keys))

    def TableD(self, tid, msg=''):
        """gets a TABLEDx (TABLED1, TABLED2, TABLED3, TABLED4)"""
        try:
            return self.tables_d[tid]
        except KeyError:
            table_keys = np.unique(list(self.tables.keys()))
            tabled_keys = np.unique(list(self.tables_d.keys()))
            tablem_keys = np.unique(list(self.tables_m.keys()))
            raise KeyError('tid=%s not found%s.  Allowed TABLEDs=%s; TABLEs=%s; TABLEMs=%s'
                           % (tid, msg, tabled_keys, table_keys, tablem_keys))

    def TableM(self, tid, msg=''):
        """gets a TABLEx (TABLEM1, TABLEM2, TABLEM3, TABLEM4)"""
        try:
            return self.tables_m[tid]
        except KeyError:
            table_keys = np.unique(list(self.tables.keys()))
            tabled_keys = np.unique(list(self.tables_d.keys()))
            tablem_keys = np.unique(list(self.tables_m.keys()))
            raise KeyError('tid=%s not found%s.  Allowed TABLEMs=%s; TABLEs=%s; TABLEDs=%s'
                           % (tid, msg, tablem_keys, table_keys, tabled_keys))

    def RandomTable(self, tid, msg=''):
        try:
            return self.random_tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, np.unique(list(self.random_tables.keys()))))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid, msg=''):
        """gets an NLPARM"""
        try:
            return self.nlparms[nid]
        except KeyError:
            raise KeyError('nid=%s not found%s.  Allowed NLPARMs=%s'
                           % (nid, msg, np.unique(list(self.nlparms.keys()))))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname, msg=''):
        """gets a DMIG"""
        try:
            return self.dmig[dname]
        except KeyError:
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, np.unique(list(self.dmig.keys()))))

    def DEQATN(self, equation_id, msg=''):
        """gets a DEQATN"""
        try:
            return self.dequations[equation_id]
        except KeyError:
            raise KeyError('equation_id=%s not found%s.  Allowed DEQATNs=%s'
                           % (equation_id, msg, np.unique(list(self.dequations.keys()))))
