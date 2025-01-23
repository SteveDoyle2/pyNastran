"""defines various methods to access low level BDF data"""
from __future__ import annotations
from itertools import chain
from collections import defaultdict
from typing import Optional, Iterable, Any, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.utils.numpy_utils import integer_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import (
        BDF, Element, Property, Material, ThermalMaterial, RigidElement,
        CAEROs, PAEROs, SPLINEs)
    from pyNastran.bdf.cards.coordinate_systems import Coord
    from pyNastran.bdf.cards.nodes import POINT, GRID, SPOINT, EPOINT # , SPOINTs, EPOINTs, SEQGP, GRIDB
    from pyNastran.bdf.cards.aero.aero import (
        #AECOMP, AECOMPL
        AEFACT, AELINK,
        AELIST,
        AEPARM, AESURF, # AESURFS,
        CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
        PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
        #MONPNT1, MONPNT2, MONPNT3,
        SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5,
    )
    from pyNastran.bdf.cards.aero.static_loads import AESTAT, AEROS, TRIM # CSSCHD, TRIM2, DIVERG
    from pyNastran.bdf.cards.aero.dynamic_loads import AERO, FLFACT, FLUTTER, GUST
    # MKAERO1, MKAERO2
    # ----------------------------------------------------

    from pyNastran.bdf.cards.elements.mass import CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4
    from pyNastran.bdf.cards.properties.mass import PMASS, NSM1, NSML, NSML1, NSM, NSMADD
    from pyNastran.bdf.cards.constraints import (SPCADD, SPC1, # SPC,
                                                 #MPC,
                                                 MPCADD) # SUPORT1, SUPORT
    #from pyNastran.bdf.cards.deqatn import DEQATN
    from pyNastran.bdf.cards.dynamic import (
        #DELAY, DPHASE,
        NLPARM)
    from pyNastran.bdf.cards.elements.rigid import (
        RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE, RSSCON)
    from pyNastran.bdf.cards.loads.loads import SLOAD, DAREA
    from pyNastran.bdf.cards.loads.dloads import DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
    from pyNastran.bdf.cards.loads.static_loads import (LOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                        FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                        PLOAD, PLOAD1, PLOAD2, PLOAD4)

    from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                               MAT8, MAT9, MAT10, MAT11, MAT3D,
                                               MATG, MATHE, MATHP, EQUIV)

    from pyNastran.bdf.cards.methods import EIGC, EIGR, EIGRL
    #from pyNastran.bdf.cards.nodes import GRID, SPOINTs, EPOINTs, POINT

    from pyNastran.bdf.cards.aero.zona import (
        #ACOORD, AEROZ, AESURFZ, BODY7, CAERO7, MKAEROZ, PAFOIL7, PANLST1, PANLST3,
        #SEGMESH, SPLINE1_ZONA, SPLINE2_ZONA, SPLINE3_ZONA, TRIMLNK, TRIMVAR, TRIM_ZONA,
        ZONA)

    from pyNastran.bdf.cards.optimization import (
        DCONADD, DCONSTR, DESVAR, DDVAL, # DOPTPRM, DLINK,
        DRESP1, DRESP2, DRESP3,
        DVCREL1, DVCREL2,
        DVMREL1, DVMREL2,
        DVPREL1, DVPREL2)
    from pyNastran.bdf.cards.bdf_sets import SET1, SET3
    from pyNastran.bdf.cards.thermal.thermal import PHBDY
    #from pyNastran.bdf.cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, TEMPB3,
                                                   #TEMPRB, QVOL, QVECT)
    from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                                TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                                TABLES1, TABLEST, TABLEHT, TABLEH1,
                                                TABRND1, TABRNDG)

class GetMethods(BDFAttributes):
    """defines various methods to access low level BDF data"""
    def __init__(self):
        BDFAttributes.__init__(self)

    def EmptyNode(self, nid: int, msg: str='') -> Optional[GRID | SPOINT | EPOINT]:
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
            nid_list = _unique_keys(self.nodes)
            msg = 'nid=%s is not a GRID, SPOINT, or EPOINT%s\n' % (nid, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % _unique_keys(self.spoints)
            if self.epoints:
                msg += 'epoints=%s\n' % _unique_keys(self.epoints)
            raise KeyError(msg)

    def Node(self, nid: int, msg: str='') -> GRID | SPOINT | EPOINT:
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
            nid_list = _unique_keys(self.nodes)
            msg = 'nid=%s is not a GRID, SPOINT, or EPOINT%s\n' % (nid, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % _unique_keys(self.spoints)
            if self.epoints:
                msg += 'epoints=%s\n' % _unique_keys(self.epoints)
            raise KeyError(msg)

    def EmptyNodes(self, nids: list[int],
                   msg: str='') -> list[Optional[GRID | SPOINT | EPOINT]]:
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
            nid_list = _unique_keys(self.nodes)
            msg = 'nids=%s are not a GRID, SPOINT, or EPOINT%s\n' % (bad_nids, msg)
            msg += 'nids=%s\n' % nid_list
            if self.spoints:
                msg += 'spoints=%s\n' % _unique_keys(self.spoints)
            if self.epoints:
                msg += 'epoints=%s\n' % _unique_keys(self.epoints)
            raise KeyError(msg)
        return nodes

    def Nodes(self, nids: list[int], msg: str='') -> list[GRID | SPOINT | EPOINT]:
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
                    nid_list = _unique_keys(self.gridb)
                    raise KeyError('nid=%s is not a GRIDB%s\n%s' % (nid, msg, nid_list))
                nodes.append(gridb)
        else:
            try:
                for nid in nids:
                    nodes.append(self.Node(nid, msg=msg))
            except AssertionError:
                print(msg)
                print(nids)
                raise
        return nodes

    def Point(self, nid: int, msg: str='') -> POINT:
        """Returns a POINT card"""
        if nid in self.points:
            return self.points[nid]
        else:
            assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
            nid_list = _unique_keys(self.points)
            raise KeyError('nid=%s is not a POINT%s\n%s' % (nid, msg, nid_list))

    def Points(self, nids: list[int], msg: str='') -> list[POINT]:
        """
        Returns a series of POINT objects given a list of IDs
        """
        points = []
        for nid in nids:
            points.append(self.Point(nid, msg=msg))
        return points

    def Element(self, eid: int, msg: str='') -> Any:
        """
        Gets an element

        Doesn't get rigid (RROD, RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE, RSSCON)
        or mass (CMASS1, CONM2)
        """
        try:
            return self.elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed elements=%s'
                           % (eid, msg, _unique_keys(self.elements)))

    def Elements(self, eids: list[int], msg: str='') -> list[Any]:
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
                bad_eids, msg, _unique_keys(self.elements))
            raise KeyError(msg)
        return elements

    def Mass(self, eid: int, msg: str='') -> CMASS1 | CMASS2 | CMASS3 | CMASS4 | CONM1 | CONM2:
        """gets a mass element (CMASS1, CONM2)"""
        try:
            return self.masses[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed masses=%s'
                           % (eid, msg, _unique_keys(self.masses)))

    def RigidElement(self, eid: int, msg: str='') -> RigidElement:
        """gets a rigid element (RBAR, RBE1, RBE2, RBE3, RBAR, RBAR1, RROD, RSPLINE, RSSCON)"""
        try:
            return self.rigid_elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed rigid_elements=%s'
                           % (eid, msg, _unique_keys(self.rigid_elements)))

    #--------------------
    # PROPERTY CARDS
    def Property(self, pid: int, msg: str='') -> Property:
        """
        gets an elemental property (e.g. PSOLID, PLSOLID, PCOMP, PSHELL, PSHEAR);
        not mass property (PMASS)

        """
        try:
            return self.properties[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Pids=%s'
                           % (pid, msg, self.property_ids))

    def Properties(self, pids: list[int], msg: str='') -> list[Any]:
        """
        gets one or more elemental property (e.g. PSOLID, PLSOLID,
        PCOMP, PSHELL, PSHEAR); not mass property (PMASS)

        """
        properties = []
        for pid in pids:
            properties.append(self.Property(pid, msg))
        return properties

    def PropertyMass(self, pid: int, msg: str='') -> PMASS:
        """
        gets a mass property (PMASS)
        """
        try:
            return self.properties_mass[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Mass Pids=%s'
                           % (pid, msg, _unique_keys(self.properties_mass)))

    def Phbdy(self, pid: int, msg: str='') -> PHBDY:
        """gets a PHBDY"""
        try:
            return self.phbdys[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PHBDY Pids=%s'
                           % (pid, msg, _unique_keys(self.phbdys)))

    #--------------------
    # MATERIAL CARDS

    def get_structural_material_ids(self) -> Iterable[int]:
        return self.materials.keys()

    def get_material_ids(self) -> Iterable[int]:
        """gets the material ids"""
        keys = chain(
            self.materials.keys(),
            self.thermal_materials.keys(),
            self.hyperelastic_materials.keys(),
        )
        return keys

    def get_thermal_material_ids(self) -> Iterable[int]:
        """gets the thermal material ids"""
        return self.thermal_materials.keys()

    def Material(self, mid: int, msg: str='') -> Material | ThermalMaterial:
        """gets a structural or thermal material"""
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermal_materials:
            return self.thermal_materials[mid]
        else:
            msg = '\n' + msg
            msg2 = (
                'Invalid Material ID:  mid=%s%s\nAllowed=%s' % (
                    mid, msg, _unique_keys(self.materials)
                )
            )
            raise KeyError(msg2)

    def StructuralMaterial(self, mid, msg='') -> Material:
        """gets a structural material"""
        try:
            mat = self.materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Structural Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def ThermalMaterial(self, mid: int, msg: str='') -> MAT4 | MAT5:
        """gets a thermal material"""
        try:
            mat = self.thermal_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Thermal Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def HyperelasticMaterial(self, mid: int, msg: str='') -> MATHE | MATHP:
        """gets a hyperelastic material"""
        try:
            mat = self.hyperelastic_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Hyperelastic Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def Materials(self, mids, msg='') -> list[Material]:
        """gets one or more Materials"""
        if isinstance(mids, integer_types):
            mids = [mids]
        materials = []
        for mid in mids:
            materials.append(self.Material(mid, msg))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid: int, consider_load_combinations: bool=True,
             msg: str='') -> list[LOAD | GRAV | ACCEL | ACCEL1 | SLOAD |
                                  FORCE | FORCE1 | FORCE2 |
                                  MOMENT | MOMENT1 | MOMENT2 |
                                  PLOAD | PLOAD1 | PLOAD2 | PLOAD4]:
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
        load = []
        if consider_load_combinations and sid in self.load_combinations:
            load = self.load_combinations[sid]
        elif sid in self.loads:
            load = self.loads[sid]

        if sid in self.tempds:
            load.append(self.tempds[sid])

        if len(load) == 0:
            loads_ids = list(self.loads.keys())
            load_combination_ids = list(self.load_combinations.keys())
            raise KeyError('cannot find LOAD ID=%r%s.\nAllowed loads (e.g., FORCE)=%s; LOAD=%s' % (
                sid, msg, np.unique(loads_ids), np.unique(load_combination_ids)))
        return load

    def DLoad(self, sid: int, consider_dload_combinations: bool=True, msg: str='') -> DLOAD:
        """
        Gets a DLOAD, TLOAD1, TLOAD2, etc. associated with the
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

    def get_dload_entries(self, sid: int, msg: str='') -> TLOAD1 | TLOAD2 | RLOAD1 | RLOAD2:
        """gets the dload entries (e.g., TLOAD1, TLOAD2)"""
        self.deprecated(
            "get_dload_entries(sid, msg='')",
            "DLoad(sid, consider_dload_combinations=False, msg='')",
            '1.1')
        return self.DLoad(sid, consider_dload_combinations=False, msg=msg)

    def DAREA(self, darea_id: int, msg: str='') -> DAREA:
        """gets a DAREA"""
        assert isinstance(darea_id, integer_types), darea_id
        try:
            return self.dareas[darea_id]
        except KeyError:
            raise KeyError('darea_id=%s not found%s.  Allowed DAREA=%s'
                           % (darea_id, msg, list(self.dareas.keys())))

    def DELAY(self, delay_id: int, msg: str='') -> DELAY:
        """gets a DELAY"""
        assert isinstance(delay_id, integer_types), delay_id
        try:
            return self.delays[delay_id]
        except KeyError:
            raise KeyError('delay_id=%s not found%s.  Allowed DELAY=%s'
                           % (delay_id, msg, list(self.delays.keys())))

    def DPHASE(self, dphase_id: int, msg: str='') -> DPHASE:
        """gets a DPHASE"""
        assert isinstance(dphase_id, integer_types), dphase_id
        try:
            return self.dphases[dphase_id]
        except KeyError:
            raise KeyError('dphase_id=%s not found%s.  Allowed DPHASE=%s'
                           % (dphase_id, msg, list(self.dphases.keys())))

    #--------------------
    def MPC(self, mpc_id: int, consider_mpcadd: bool=True, msg: str='') -> MPC | MPCADD:
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

    def SPC(self, spc_id: int, consider_spcadd: bool=True, msg: str='') -> SPC | SPC1 | SPCADD:
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

    def NSM(self, nsm_id: int, consider_nsmadd: bool=True,
            msg: str='') -> NSM | NSM1 | NSML | NSML1 | NSMADD:
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
    def SET1(self, set_id: int, msg: str='') -> SET1:
        """gets a SET1"""
        assert isinstance(set_id, integer_types), 'set_id=%s is not an integer\n' % set_id
        if set_id in self.sets:
            set1 = self.sets[set_id]
        else:
            raise KeyError(f'cannot find SET1 ID={set_id}{msg}.\n{msg}')
        return set1

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid: int, msg: str='') -> Coord:
        """gets an COORDx"""
        try:
            return self.coords[cid]
        except KeyError:
            cids = np.array(list(self.coord_ids))
            raise KeyError(f'cid={cid} not found{msg}.  Allowed Cids={cids}')

    #--------------------
    # AERO CARDS
    def AEList(self, aelist: int, msg: str='') -> AELIST:
        """gets an AELIST"""
        try:
            return self.aelists[aelist]
        except KeyError:
            aelists = _unique_keys(self.aelists)
            raise KeyError(f'aelist={aelist} not found{msg}.  Allowed AELIST={aelists}')

    def AEFact(self, aefact: int, msg: str='') -> AEFACT:
        """gets an AEFACT"""
        try:
            return self.aefacts[aefact]
        except KeyError:
            aefacts = _unique_keys(self.aefacts)
            raise KeyError(f'aefact={aefact} not found{msg}.  Allowed AEFACT={aefacts}')

    def Trim(self, trim_id: int, msg: str='') -> TRIM:
        """gets an TRIM"""
        try:
            return self.trims[trim_id]
        except KeyError:
            raise KeyError('TRIM=%s not found%s.  Allowed TRIM=%s'
                           % (trim_id, msg, _unique_keys(self.trims)))

    def AESurf(self, aesurf_name: str, msg: str='') -> AESURF:
        """gets an AESURF"""
        #if isinstance(aesurf_name, integer_types):
            #aesurf_id = aesurf_name
            #try:
                #return self.aesurf[aesurf_id]
            #except KeyError:
                #raise KeyError('aesurf=%s not found%s.  Allowed AESURF=%s'
                               #% (aesurf_id, msg, _unique_keys(self.aesurf)))
        #else:
        assert isinstance(aesurf_name, str), f'aesurf_name={aesurf_name!r}'

        for aesurf_int, aesurf in self.aesurf.items():
            if aesurf.label == aesurf_name:
                return aesurf

        aesurf_names = list(aesurf.label for aesurf in self.aesurf.values())
        aeparam_names = list(aeparm.label for aeparm in self.aeparams.values())
        raise KeyError('aesurf=%r not found%s.  Allowed AESURF=%s.  Allowed AEPARM=%s.' % (
            aesurf_name, msg, aesurf_names, aeparam_names))

    def AESurf_int(self, aesurf_id: int, msg: str='') -> AESURF:
        """gets an AESURF"""
        try:
            return self.aesurf[aesurf_id]
        except KeyError:
            raise KeyError('aesurf=%s not found%s.  Allowed AESURF=%s'
                           % (aesurf_id, msg, _unique_keys(self.aesurf)))

    def Acsid(self, msg: str='') -> Coord:
        """gets the aerodynamic coordinate system"""
        if self.aero is not None:
            acsid_aero = self.aero.Acsid()
        if self.aeros is not None:
            acsid_aeros = self.aeros.Acsid()

        if self.aero is not None:
            if self.aeros is not None:
                assert acsid_aero == acsid_aeros, f'AERO acsid={acsid_aero}, AEROS acsid={acsid_aeros}'
            coord = self.Coord(acsid_aero, msg=msg)
        elif self.aeros is not None:
            coord = self.Coord(acsid_aeros, msg=msg)
        else:
            msg = 'neither AERO nor AEROS cards exist%s.' % msg
            raise RuntimeError(msg)
        return coord

    def safe_acsid(self, msg: str='') -> Optional[Coord]:
        """gets the aerodynamic coordinate system"""
        if self.aero is not None:
            acsid_aero = self.aero.Acsid()
        if self.aeros is not None:
            acsid_aeros = self.aeros.Acsid()

        if self.aero is not None:
            if self.aeros is not None:
                assert acsid_aero == acsid_aeros, f'AERO acsid={acsid_aero}, AEROS acsid={acsid_aeros}'
            coord = self.Coord(acsid_aero, msg=msg)
        elif self.aeros is not None:
            coord = self.Coord(acsid_aeros, msg=msg)
        else:
            ## TODO: consider changing this...
            self.log.error(f'neither AERO nor AEROS cards exist; assuming global (cid=0){msg}.')
            return self.Coord(0, msg=msg)
        return coord

    def Aero(self, msg: str='') -> AERO:
        """gets the AERO"""
        if self.aero is not None:
            return self.aero
        raise RuntimeError(f'no AERO card found{msg}.')

    def Aeros(self, msg: str='') -> AEROS:
        """gets the AEROS"""
        if self.aeros is not None:
            return self.aeros
        raise RuntimeError(f'no AEROS card found{msg}.')

    def Spline(self, eid: int, msg: str='') -> SPLINEs:
        """gets a SPLINEx"""
        try:
            return self.splines[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, _unique_keys(self.splines)))

    def CAero(self, eid: int, msg: str='') -> CAEROs:
        """gets an CAEROx"""
        try:
            return self.caeros[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, _unique_keys(self.caeros)))

    def PAero(self, pid: int, msg: str='') -> PAEROs:
        """gets a PAEROx"""
        try:
            return self.paeros[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, _unique_keys(self.paeros)))

    def Gust(self, sid: int, msg: str='') -> GUST:
        """gets a GUST"""
        try:
            return self.gusts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, _unique_keys(self.gusts)))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aestat_name: str, msg: str='') -> AESTAT:
        """gets an AESTAT"""
        #if isinstance(aesurf_name, integer_types):
            #aesurf_id = aesurf_name
            #try:
                #return self.aesurf[aesurf_id]
            #except KeyError:
                #raise KeyError('aesurf=%s not found%s.  Allowed AESURF=%s'
                               #% (aesurf_id, msg, _unique_keys(self.aesurf)))
        #else:
        assert isinstance(aestat_name, str), f'aestat_name={aestat_name!r}'

        for aestat_int, aestat in self.aestats.items():
            if aestat.label == aestat_name:
                return aestat

        aesurf_names = list(aesurf.label for aesurf in self.aesurf.values())
        aeparam_names = list(aeparm.label for aeparm in self.aeparams.values())
        aestat_names = list(aeparm.label for aeparm in self.aeparams.values())
        raise KeyError('aestat=%r not found%s.  Allowed AESURF=%s.  Allowed AEPARM=%s.  Allowed AESTAT=%s' % (
            aestat_name, msg, aesurf_names, aeparam_names, aestat_names))

    def AEStat_int(self, aid: int, msg: str='') -> AESTAT:
        """gets an AESTAT"""
        try:
            return self.aestats[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, _unique_keys(self.aestats)))

    def AELIST(self, aid: int, msg: str='') -> AELIST:
        """gets an AELIST"""
        try:
            return self.aelists[aid]
        except KeyError:
            raise KeyError('id=%s not found%s.  Allowed AELISTs=%s'
                           % (aid, msg, _unique_keys(self.aelists)))

    def AELink(self, link_id: int, msg: str='') -> AELINK:
        """gets an AELINK"""
        try:
            return self.aelinks[link_id]
        except KeyError:
            raise KeyError('link_id=%s not found%s.  Allowed AELINKs=%s'
                           % (link_id, msg, _unique_keys(self.aelinks)))

    def AEParam(self, aeparm_name: str, msg: str='') -> AEPARM:
        """gets an AEPARM"""
        #if isinstance(aeparam_name, integer_types):
            #aeparam_id = aeparam_name
            #try:
                #return self.aeparams[aesurf_id]
            #except KeyError:
                #raise KeyError('aesurf=%s not found%s.  Allowed AEPARM=%s'
                               #% (aesurf_id, msg, _unique_keys(self.aesurf)))
        #else:
        assert isinstance(aeparm_name, str), f'aeparm_name={aeparm_name!r}'

        for aeparam_int, aeparam in self.aeparams.items():
            if aeparam.label == aeparm_name:
                return aeparam

        aesurf_names = list(aesurf.label for aesurf in self.aesurf.values())
        aeparam_names = list(aeparm.label for aeparm in self.aeparams.values())
        raise KeyError('aeparam=%r not found%s.  Allowed AESURF=%s.  Allowed AEPARM=%s.' % (
            aeparm_name, msg, aesurf_names, aeparam_names))

    def AEParam_int(self, aid: int, msg: str='') -> AEPARM:
        """gets an AEPARM"""
        try:
            return self.aeparams[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, _unique_keys(self.aeparams)))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid: int, msg: str='') -> FLFACT:
        """gets an FLFACT"""
        try:
            return self.flfacts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, _unique_keys(self.flfacts)))

    def Flutter(self, fid: int, msg: str='') -> FLUTTER:
        """gets a FLUTTER"""
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, _unique_keys(self.flutters)))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid: int, msg: str='') -> list[DCONSTR | DCONADD]:
        """gets a DCONSTR"""
        try:
            return self.dconstrs[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, _unique_keys(self.dconstrs)))

    def DResp(self, dresp_id: int, msg: str='') -> DRESP1 | DRESP2 | DRESP3:
        """gets a DRESPx"""
        try:
            return self.dresps[dresp_id]
        except KeyError:
            raise KeyError('dresp_id=%s not found%s.  Allowed DRESPx=%s'
                           % (dresp_id, msg, _unique_keys(self.dresps)))

    def Desvar(self, desvar_id: int, msg: int='') -> DESVAR:
        """gets a DESVAR"""
        try:
            return self.desvars[desvar_id]
        except KeyError:
            raise KeyError('desvar_id=%s not found%s.  Allowed DESVARs=%s'
                           % (desvar_id, msg, _unique_keys(self.desvars)))

    def DDVal(self, oid: int, msg: str='') -> DDVAL:
        """gets a DDVAL"""
        try:
            return self.ddvals[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, _unique_keys(self.ddvals)))

    def DVcrel(self, dv_id: int, msg: str='') -> DVCREL1 | DVCREL2:
        """gets a DVCREL1/DVCREL2"""
        try:
            return self.dvcrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVCRELx=%s'
                           % (dv_id, msg, _unique_keys(self.dvcrels)))

    def DVmrel(self, dv_id: int, msg: str='') -> DVMREL1 | DVMREL2:
        """gets a DVMREL1/DVMREL2"""
        try:
            return self.dvmrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVMRELx=%s'
                           % (dv_id, msg, _unique_keys(self.dvmrels)))

    def DVprel(self, dv_id: int, msg: str='') -> DVPREL1 | DVPREL2:
        """gets a DVPREL1/DVPREL2"""
        try:
            return self.dvprels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVPRELx=%s'
                           % (dv_id, msg, _unique_keys(self.dvprels)))

    #--------------------
    # SET CARDS

    def Set(self, sid: int, msg: str='') -> SET1 | SET3:
        """gets a SET, SET1, SET2, or SET3 card"""
        try:
            return self.sets[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, _unique_keys(self.sets)))

    #--------------------
    # METHOD CARDS
    def Method(self, sid: int, msg: str='') -> EIGR | EIGRL:
        """gets a METHOD (EIGR, EIGRL)"""
        try:
            return self.methods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, _unique_keys(self.methods)))

    def CMethod(self, sid: int, msg: str='') -> EIGC:
        """gets a METHOD (EIGC)"""
        try:
            return self.cMethods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, _unique_keys(self.cMethods)))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg='') -> TABLES1 | TABLEST | TABLEH1 | TABLEHT:
        """gets a TABLES1, TABLEST, TABLEH1, TABLEHT"""
        try:
            return self.tables[tid]
        except KeyError:
            table_keys = _unique_keys(self.tables)
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, table_keys))

    def TableD(self, tid: int, msg: str='') -> TABLED1 | TABLED2 | TABLED3 | TABLED4:
        """gets a TABLEDx (TABLED1, TABLED2, TABLED3, TABLED4)"""
        try:
            return self.tables_d[tid]
        except KeyError:
            table_keys = _unique_keys(self.tables)
            tabled_keys = _unique_keys(self.tables_d)
            tablem_keys = _unique_keys(self.tables_m)
            raise KeyError('tid=%s not found%s.  Allowed TABLEDs=%s; TABLEs=%s; TABLEMs=%s'
                           % (tid, msg, tabled_keys, table_keys, tablem_keys))

    def TableM(self, table_id: int, msg: str='') -> TABLEM1 | TABLEM2 | TABLEM3 | TABLEM4:
        """gets a TABLEx (TABLEM1, TABLEM2, TABLEM3, TABLEM4)"""
        try:
            return self.tables_m[table_id]
        except KeyError:
            table_keys = _unique_keys(self.tables)
            tabled_keys = _unique_keys(self.tables_d)
            tablem_keys = _unique_keys(self.tables_m)
            raise KeyError('table_id=%s not found%s.  Allowed TABLEMs=%s; TABLEs=%s; TABLEDs=%s'
                           % (table_id, msg, tablem_keys, table_keys, tabled_keys))

    def RandomTable(self, table_id: int, msg: str='') -> TABRND1 | TABRNDG:
        """gets a TABRND1 / TABRNDG"""
        try:
            return self.random_tables[table_id]
        except KeyError:
            raise KeyError('table_id=%s not found%s.  Allowed TABLEs=%s'
                           % (table_id, msg, _unique_keys(self.random_tables)))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nlparm_id: int, msg: str='') -> NLPARM:
        """gets an NLPARM"""
        try:
            return self.nlparms[nlparm_id]
        except KeyError:
            raise KeyError('nlparm_id=%s not found%s.  Allowed NLPARMs=%s'
                           % (nlparm_id, msg, _unique_keys(self.nlparms)))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname: str, msg: str='') -> DMIG:
        """gets a DMIG"""
        try:
            return self.dmig[dname]
        except KeyError:
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, _unique_keys(self.dmig)))

    def DEQATN(self, equation_id: int, msg: str='') -> DEQATN:
        """gets a DEQATN"""
        try:
            return self.dequations[equation_id]
        except KeyError:
            raise KeyError('equation_id=%s not found%s.  Allowed DEQATNs=%s'
                           % (equation_id, msg, _unique_keys(self.dequations)))

def get_pid_to_eid_map(model: BDF) -> dict[int, list[int]]:
    pid_to_eid_map = defaultdict(set)
    for eid, elem in model.elements.items():
        pid = elem.pid
        pid_to_eid_map[pid].add(eid)
    return pid_to_eid_map

def get_pid_to_nid_map(model: BDF) -> dict[int, list[int]]:
    """TODO: doesn't support CONROD"""
    from collections import defaultdict
    get_pid_to_eid_map(model)
    property_to_nodes_map = defaultdict(set)
    for eid, elem in model.elements.items():
        pid = elem.pid
        property_to_nodes_map[pid].update(elem.node_ids)

    property_to_nodes_map2 = {}
    for pid, nodes in property_to_nodes_map.items():
        nodes_list = list(nodes)
        nodes_list.sort()
        property_to_nodes_map2[pid] = nodes_list
    return property_to_nodes_map2

def _unique_keys(mydict: dict[int, Any]) -> str:
    """helper method"""
    return np.unique(list(mydict.keys()))
