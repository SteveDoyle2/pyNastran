"""defines various methods to access low level BDF data"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import chain
import numpy as np

from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.utils.numpy_utils import integer_types

class GetMethods(BDFAttributes):
    """defines various methods to access low level BDF data"""
    def __init__(self):
        # type: () -> None
        BDFAttributes.__init__(self)

    def EmptyNode(self, nid, msg=''):
        # type: (int, str) -> Optional[Any]
        """
        Gets a GRID/SPOINT/EPOINT object, but allows for empty nodes
        (i.e., the CTETRA10 supports empty nodes, but the CTETRA4 does not).

        Parameters
        ----------
        nid : int / None
            the node id
            0, None : indicate blank
        msg : str; default=''
            a debugging message

        """
        if nid == 0 or nid is None:
            return None
        elif nid in self.nodes:
            return self.nodes[nid]
        elif nid in self.spoints:
            return self.spoints[nid]
        elif nid in self.epoints:
            return self.epoints[nid]
        else:
            assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
            nid_list = np.unique(list(self.nodes.keys()))
            msg = 'nid=%s is not a GRID, SPOINT, or EPOINT%s\n' % (nid, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % np.unique(list(self.spoints.keys()))
            if self.epoints:
                msg += 'epoints=%s\n' % np.unique(list(self.epoints.keys()))
            raise KeyError(msg)

    def Node(self, nid, msg=''):
        # type: (int, str) -> Any
        """
        Gets a GRID/SPOINT/EPOINT object.  This method does not allow for empty nodes
        (i.e., the CTETRA10 supports empty nodes, but the CTETRA4 does not).

        Parameters
        ----------
        nid : int
            the node id
        msg : str; default=''
            a debugging message

        """
        #assert isinstance(nid, integer_types), 'nid=%s' % str(nid)
        if nid in self.nodes:
            return self.nodes[nid]
        elif nid in self.spoints:
            return self.spoints[nid]
        elif nid in self.epoints:
            return self.epoints[nid]
        else:
            assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
            nid_list = np.unique(list(self.nodes.keys()))
            msg = 'nid=%s is not a GRID, SPOINT, or EPOINT%s\n' % (nid, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % np.unique(list(self.spoints.keys()))
            if self.epoints:
                msg += 'epoints=%s\n' % np.unique(list(self.epoints.keys()))
            raise KeyError(msg)

    def EmptyNodes(self, nids, msg=''):
        # type: (List[int], str) -> List[Optional[Any]]
        """
        Returns a series of node objects given a list of IDs

        """
        nodes = []
        bad_nids = []
        for nid in nids:
            try:
                nodes.append(self.EmptyNode(nid, msg=msg))
            except KeyError:
                bad_nids.append(nid)

        if bad_nids:
            nid_list = np.unique(list(self.nodes.keys()))
            msg = 'nids=%s are not a GRID, SPOINT, or EPOINT%s\n' % (bad_nids, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % np.unique(list(self.spoints.keys()))
            if self.epoints:
                msg += 'epoints=%s\n' % np.unique(list(self.epoints.keys()))
            raise KeyError(msg)

        return nodes

    def Nodes(self, nids, msg=''):
        # type: (List[int], str) -> List[Any]
        """
        Returns a series of node objects given a list of IDs

        """
        nodes = []
        #self.axic
        if self._is_axis_symmetric and self.axif is not None:
            # GRIDB is only active if AXIF exists
            for nid in nids:
                try:
                    gridb = self.gridb[nid]
                except KeyError:
                    assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
                    nid_list = np.unique(list(self.gridb.keys()))
                    raise KeyError('nid=%s is not a GRIDB%s\n%s' % (nid, msg, nid_list))
                nodes.append(gridb)
        else:
            for nid in nids:
                nodes.append(self.Node(nid, msg=msg))
        return nodes

    def Point(self, nid, msg=''):
        # type: (int, str) -> Any
        """Returns a POINT card"""
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
        # type: (int, str) -> Any
        """
        Gets an element

        Doesn't get rigid (RROD, RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE, RSSCON)
        or mass (CMASS1, CONM2)
        """
        try:
            return self.elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed elements=%s'
                           % (eid, msg, np.unique(list(self.elements.keys()))))

    def Elements(self, eids, msg=''):
        """
        Gets an series of elements

        Doesn't get rigid (RROD, RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE, RSSCON)
        or mass (CMASS1, CONM2)

        """
        elements = []
        bad_eids = []
        for eid in eids:
            try:
                elements.append(self.Element(eid, msg))
            except KeyError:
                bad_eids.append(eid)
        if bad_eids:
            msg = 'eids=%s not found%s.  Allowed elements=%s' % (
                bad_eids, msg, np.unique(list(self.elements.keys())))
            raise KeyError(msg)
        return elements

    def Mass(self, eid, msg=''):
        # type: (int, str) -> Any
        """gets a mass element (CMASS1, CONM2)"""
        try:
            return self.masses[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed masses=%s'
                           % (eid, msg, np.unique(list(self.masses.keys()))))

    def RigidElement(self, eid, msg=''):
        # type: (int, str) -> Any
        """gets a rigid element (RBAR, RBE2, RBE3, RBAR, RBAR1, RROD, RSPLINE, RSSCON)"""
        try:
            return self.rigid_elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed rigid_elements=%s'
                           % (eid, msg, np.unique(list(self.rigid_elements.keys()))))

    #--------------------
    # PROPERTY CARDS

    def Property(self, pid, msg=''):
        # type: (int, str) -> Any
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
        """
        gets one or more elemental property (e.g. PSOLID, PLSOLID,
        PCOMP, PSHELL, PSHEAR); not mass property (PMASS)

        """
        properties = []
        for pid in pids:
            properties.append(self.Property(pid, msg))
        return properties

    def PropertyMass(self, pid, msg=''):
        # type: (int, str) -> Any
        """
        gets a mass property (PMASS)
        """
        try:
            return self.properties_mass[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Mass Pids=%s'
                           % (pid, msg, np.unique(list(self.properties_mass.keys()))))

    def Phbdy(self, pid, msg=''):
        # type: (int, str) -> Any
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
        """gets the material ids"""
        keys = chain(
            self.materials.keys(),
            self.thermal_materials.keys(),
            self.hyperelastic_materials.keys(),
        )
        return keys

    def get_thermal_material_ids(self):
        """gets the thermal material ids"""
        return self.thermal_materials.keys()

    def Material(self, mid, msg=''):
        # type: (int, str) -> Any
        """gets a structural or thermal material"""
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermal_materials:
            return self.thermal_materials[mid]
        else:
            msg = '\n' + msg
            msg2 = (
                'Invalid Material ID:  mid=%s%s\nAllowed=%s' % (
                    mid, msg, np.unique(list(self.materials.keys()))
                )
            )
            raise KeyError(msg2)

    def StructuralMaterial(self, mid, msg=''):
        # type: (int, str) -> Any
        """gets a structural material"""
        try:
            mat = self.materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Structural Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def ThermalMaterial(self, mid, msg=''):
        # type: (int, str) -> Any
        """gets a thermal material"""
        try:
            mat = self.thermal_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Thermal Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def HyperelasticMaterial(self, mid, msg=''):
        # type: (int, str) -> Any
        """gets a hyperelastic material"""
        try:
            mat = self.hyperelastic_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Hyperelastic Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def Materials(self, mids, msg=''):
        """gets one or more Materials"""
        if isinstance(mids, integer_types):
            mids = [mids]
        materials = []
        for mid in mids:
            materials.append(self.Material(mid, msg))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid, consider_load_combinations=True, msg=''):
        """
        Gets an LOAD or FORCE/PLOAD4/etc.

        Parameters
        ----------
        sid : int
            the LOAD id
        consider_load_combinations : bool; default=True
            LOADs should not be considered when referenced from an LOAD card
            from a case control, True should be used.
        msg : str
            additional message to print when failing

        """
        assert isinstance(sid, integer_types), 'sid=%s is not an integer; type=%s\n' % (sid, type(sid))
        if consider_load_combinations and sid in self.load_combinations:
            load = self.load_combinations[sid]
        elif sid in self.loads:
            load = self.loads[sid]
        else:
            loads_ids = list(self.loads.keys())
            load_combination_ids = list(self.load_combinations.keys())
            raise KeyError('cannot find LOAD ID=%r%s.\nAllowed loads (e.g., FORCE)=%s; LOAD=%s' % (
                sid, msg, np.unique(loads_ids), np.unique(load_combination_ids)))
        return load

    def DLoad(self, sid, consider_dload_combinations=True, msg=''):
        """
        Gets a DLOAD, TLOAD1, TLOAD2, etc. associcated with the
        Case Control DLOAD entry

        """
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid

        if consider_dload_combinations and sid in self.dloads:
            dload = self.dloads[sid]
        elif sid in self.dload_entries:
            dload = self.dload_entries[sid]
        else:
            dloads_ids = list(self.dload_entries.keys())
            dload_combination_ids = list(self.dloads.keys())
            raise KeyError('cannot find DLOAD ID=%r%s.\nAllowed dloads '
                           '(e.g., TLOAD1)=%s; DLOAD=%s' % (
                               sid, msg, np.unique(dloads_ids),
                               np.unique(dload_combination_ids)))
        return dload

    def get_dload_entries(self, sid, msg=''):
        """gets the dload entries (e.g., TLOAD1, TLOAD2)"""
        self.deprecated(
            "get_dload_entries(sid, msg='')",
            "DLoad(sid, consider_dload_combinations=False, msg='')",
            '1.1')
        return self.DLoad(sid, consider_dload_combinations=False, msg=msg)

    def DAREA(self, darea_id, msg=''):
        """gets a DAREA"""
        assert isinstance(darea_id, integer_types), darea_id
        try:
            return self.dareas[darea_id]
        except KeyError:
            raise KeyError('darea_id=%s not found%s.  Allowed DAREA=%s'
                           % (darea_id, msg, list(self.dareas.keys())))

    def DELAY(self, delay_id, msg=''):
        # type: (int, str) -> Any
        """gets a DELAY"""
        assert isinstance(delay_id, integer_types), delay_id
        try:
            return self.delays[delay_id]
        except KeyError:
            raise KeyError('delay_id=%s not found%s.  Allowed DELAY=%s'
                           % (delay_id, msg, list(self.delays.keys())))

    def DPHASE(self, dphase_id, msg=''):
        # type: (int, str) -> Any
        """gets a DPHASE"""
        assert isinstance(dphase_id, integer_types), dphase_id
        try:
            return self.dphases[dphase_id]
        except KeyError:
            raise KeyError('dphase_id=%s not found%s.  Allowed DPHASE=%s'
                           % (dphase_id, msg, list(self.dphases.keys())))

    #--------------------
    def MPC(self, mpc_id, consider_mpcadd=True, msg=''):
        """
        Gets an MPCADD or MPC

        Parameters
        ----------
        mpc_id : int
            the MPC id
        consider_mpcadd : bool; default=True
            MPCADDs should not be considered when referenced from an MPCADD
            from a case control, True should be used.
        msg : str
            additional message to print when failing

        """
        assert isinstance(mpc_id, integer_types), 'mpc_id=%s is not an integer\n' % mpc_id
        if consider_mpcadd and mpc_id in self.mpcadds:
            constraint = self.mpcadds[mpc_id]
        elif mpc_id in self.mpcs:
            constraint = self.mpcs[mpc_id]
        else:
            mpc_ids = list(self.mpcs.keys())
            mpcadd_ids = list(self.mpcadds.keys())
            raise KeyError('cannot find MPC ID=%r%s.\nAllowed MPCs=%s; MPCADDs=%s' % (
                mpc_id, msg, np.unique(mpc_ids), np.unique(mpcadd_ids)))
        return constraint

    def SPC(self, spc_id, consider_spcadd=True, msg=''):
        """
        Gets an SPCADD or SPC

        Parameters
        ----------
        spc_id : int
            the SPC id
        consider_spcadd : bool; default=True
            SPCADDs should not be considered when referenced from an SPCADD
            from a case control, True should be used.
        msg : str
            additional message to print when failing

        """
        assert isinstance(spc_id, integer_types), 'spc_id=%s is not an integer\n' % spc_id
        if consider_spcadd and spc_id in self.spcadds:
            constraint = self.spcadds[spc_id]
        elif spc_id in self.spcs:
            constraint = self.spcs[spc_id]
        else:
            spc_ids = list(self.spcs.keys())
            spcadd_ids = list(self.spcadds.keys())
            raise KeyError('cannot find SPC ID=%r%s.\nAllowed SPCs=%s; SPCADDs=%s' % (
                spc_id, msg, np.unique(spc_ids), np.unique(spcadd_ids)))
        return constraint

    def NSM(self, nsm_id, consider_nsmadd=True, msg=''):
        """
        Gets an LOAD or FORCE/PLOAD4/etc.

        Parameters
        ----------
        nsm_id : int
            the LOAD id
        consider_nsmadd : bool; default=True
            NSMADDs should not be considered when referenced from an NSM card
            from a case control, True should be used.
        msg : str
            additional message to print when failing

        """
        assert isinstance(nsm_id, integer_types), 'nsm_id=%s is not an integer\n' % nsm_id
        if consider_nsmadd and nsm_id in self.nsmadds:
            nsm = self.nsmadds[nsm_id]
        elif nsm_id in self.nsms:
            nsm = self.nsms[nsm_id]
        else:
            nsm_ids = list(self.nsms.keys())
            nsmadd_ids = list(self.nsmadds.keys())
            raise KeyError('cannot find NSM ID=%r%s.\nAllowed NSMs (e.g., NSM1)=%s; '
                           'NSMADDs=%s; consider_nsmadd=%s' % (
                               nsm_id, msg, np.unique(nsm_ids), np.unique(nsmadd_ids),
                               consider_nsmadd))
        return nsm

    #--------------------
    # Sets
    def SET1(self, set_id, msg=''):
        # type: (int, str) -> Any
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
        # type: (int, str) -> Any
        """gets an COORDx"""
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coord_ids))

    #--------------------
    # AERO CARDS

    def AEList(self, aelist, msg=''):
        # type: (int, str) -> Any
        """gets an AELIST"""
        try:
            return self.aelists[aelist]
        except KeyError:
            raise KeyError('aelist=%s not found%s.  Allowed AELIST=%s'
                           % (aelist, msg, np.unique(list(self.aelists.keys()))))

    def AEFact(self, aefact, msg=''):
        # type: (int, str) -> Any
        """gets an AEFACT"""
        try:
            return self.aefacts[aefact]
        except KeyError:
            raise KeyError('aefact=%s not found%s.  Allowed AEFACT=%s'
                           % (aefact, msg, np.unique(list(self.aefacts.keys()))))

    def AESurf(self, aesurf_id, msg=''):
        # type: (int, str) -> Any
        """gets an AESURF"""
        try:
            return self.aesurf[aesurf_id]
        except KeyError:
            raise KeyError('aesurf=%s not found%s.  Allowed AESURF=%s'
                           % (aesurf_id, msg, np.unique(list(self.aesurf.keys()))))

    def Acsid(self, msg=''):
        # type: (int, str) -> Any
        """gets the aerodynamic coordinate system"""
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

    def safe_acsid(self, msg=''):
        # type: (str) -> Any
        """gets the aerodynamic coordinate system"""
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
            ## TODO: consider changing this...
            self.log.error('neither AERO nor AEROS cards exist; assuming global (cid=0).')
            return self.Coord(0, msg=msg)
        return cid

    def Aero(self, msg=''):
        # type: (str) -> Any
        """gets the AERO"""
        if self.aero is not None:
            return self.aero
        raise RuntimeError('no AERO card found%s.' % (msg))

    def Aeros(self, msg=''):
        # type: (str) -> Any
        """gets the AEROS"""
        if self.aeros is not None:
            return self.aeros
        raise RuntimeError('no AEROS card found%s.' % (msg))

    def Spline(self, eid, msg=''):
        # type: (int, str) -> Any
        """gets a SPLINEx"""
        try:
            return self.splines[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, np.unique(list(self.splines.keys()))))

    def CAero(self, eid, msg=''):
        # type: (int, str) -> Any
        """gets an CAEROx"""
        try:
            return self.caeros[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, np.unique(list(self.caero_ids))))

    def PAero(self, pid, msg=''):
        # type: (int, str) -> Any
        """gets a PAEROx"""
        try:
            return self.paeros[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, np.unique(list(self.paeros.keys()))))

    def Gust(self, sid, msg=''):
        # type: (int, str) -> Any
        """gets a GUST"""
        try:
            return self.gusts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, np.unique(list(self.gusts.keys()))))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid, msg=''):
        # type: (int, str) -> Any
        """gets an AESTAT"""
        try:
            return self.aestats[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, np.unique(list(self.aestats.keys()))))

    def AELIST(self, aid, msg=''):
        # type: (int, str) -> Any
        """gets an AELIST"""
        try:
            return self.aelists[aid]
        except KeyError:
            raise KeyError('id=%s not found%s.  Allowed AELISTs=%s'
                           % (aid, msg, np.unique(list(self.aelists.keys()))))

    def AELink(self, link_id, msg=''):
        # type: (int, str) -> Any
        """gets an AELINK"""
        try:
            return self.aelinks[link_id]
        except KeyError:
            raise KeyError('link_id=%s not found%s.  Allowed AELINKs=%s'
                           % (link_id, msg, np.unique(list(self.aelinks.keys()))))

    def AEParam(self, aid, msg=''):
        # type: (int, str) -> Any
        """gets an AEPARM"""
        try:
            return self.aeparams[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, np.unique(list(self.aeparams.keys()))))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid, msg=''):
        # type: (int, str) -> Any
        """gets an FLFACT"""
        try:
            return self.flfacts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, np.unique(list(self.flfacts.keys()))))

    def Flutter(self, fid, msg=''):
        # type: (int, str) -> Any
        """gets a FLUTTER"""
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, np.unique(list(self.flutters.keys()))))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid, msg=''):
        # type: (int, str) -> List[Any]
        """gets a DCONSTR"""
        try:
            return self.dconstrs[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, np.unique(list(self.dconstrs.keys()))))

    def DResp(self, dresp_id, msg=''):
        # type: (int, str) -> Any
        """gets a DRESPx"""
        try:
            return self.dresps[dresp_id]
        except KeyError:
            raise KeyError('dresp_id=%s not found%s.  Allowed DRESPx=%s'
                           % (dresp_id, msg, np.unique(list(self.dresps.keys()))))

    def Desvar(self, desvar_id, msg=''):
        # type: (int, str) -> Any
        """gets a DESVAR"""
        try:
            return self.desvars[desvar_id]
        except KeyError:
            raise KeyError('desvar_id=%s not found%s.  Allowed DESVARs=%s'
                           % (desvar_id, msg, np.unique(list(self.desvars.keys()))))

    def DDVal(self, oid, msg=''):
        # type: (int, str) -> Any
        """gets a DDVAL"""
        try:
            return self.ddvals[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, np.unique(list(self.ddvals.keys()))))

    def DVcrel(self, dv_id, msg=''):
        # type: (int, str) -> Any
        """gets a DVCREL1/DVCREL2"""
        try:
            return self.dvcrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVCRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvcrels.keys()))))

    def DVmrel(self, dv_id, msg=''):
        # type: (int, str) -> Any
        """gets a DVMREL1/DVMREL2"""
        try:
            return self.dvmrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVMRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvmrels.keys()))))

    def DVprel(self, dv_id, msg=''):
        # type: (int, str) -> Any
        """gets a DVPREL1/DVPREL2"""
        try:
            return self.dvprels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVPRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvprels.keys()))))

    #--------------------
    # SET CARDS

    def Set(self, sid, msg=''):
        # type: (int, str) -> Any
        """gets a SET, SET1, SET2, or SET3 card"""
        try:
            return self.sets[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, np.unique(list(self.sets.keys()))))

    #--------------------
    # METHOD CARDS
    def Method(self, sid, msg=''):
        # type: (int, str) -> Any
        """gets a METHOD (EIGR, EIGRL)"""
        try:
            return self.methods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, np.unique(list(self.methods.keys()))))

    def CMethod(self, sid, msg=''):
        # type: (int, str) -> Any
        """gets a METHOD (EIGC)"""
        try:
            return self.cMethods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, np.unique(list(self.cMethods.keys()))))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg=''):
        # type: (int, str) -> Any
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
        # type: (int, str) -> Any
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
        # type: (int, str) -> Any
        """gets a TABRND1 / TABRNDG"""
        try:
            return self.random_tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, np.unique(list(self.random_tables.keys()))))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid, msg=''):
        # type: (int, str) -> Any
        """gets an NLPARM"""
        try:
            return self.nlparms[nid]
        except KeyError:
            raise KeyError('nid=%s not found%s.  Allowed NLPARMs=%s'
                           % (nid, msg, np.unique(list(self.nlparms.keys()))))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname, msg=''):
        # type: (str, str) -> Any
        """gets a DMIG"""
        try:
            return self.dmigs[dname]
        except KeyError:
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, np.unique(list(self.dmigs.keys()))))

    def DEQATN(self, equation_id, msg=''):
        # type: (int, str) -> Any
        """gets a DEQATN"""
        try:
            return self.dequations[equation_id]
        except KeyError:
            raise KeyError('equation_id=%s not found%s.  Allowed DEQATNs=%s'
                           % (equation_id, msg, np.unique(list(self.dequations.keys()))))
