"""
Creates safe cross referencing

Safe cross-referencing skips failed xref's

"""
from __future__ import annotations
from typing import Any, TYPE_CHECKING

import numpy as np
from numpy import zeros, argsort, arange, array_equal
from pyNastran.bdf.bdf_interface.cross_reference import XrefMesh
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import (
        BDF,
        Element, Property, Material, PAEROs, CAEROs, Coord,
        AEFACT, AELIST, TABLEDs, TABLEMs)
    from pyNastran.bdf.cards.nodes import GRID
    from pyNastran.bdf.bdf_interface.cross_reference_obj import CrossReference


class SafeXrefMesh(XrefMesh):
    """
    Links up the various cards in the BDF.
    """
    def __init__(self) -> None:
        """The main BDF class defines all the parameters that are used."""
        XrefMesh.__init__(self)

    # def geom_check(self):
        # """
        # Performs various geometry checks
          # 1.  nodal uniqueness on elements
        # """
        # for elem in model.elements:
            # elem.check_unique_nodes()

    def safe_cross_reference(self, xref: bool=True,
                             xref_nodes: bool=True,
                             xref_elements: bool=True,
                             xref_nodes_with_elements: bool=False,
                             xref_properties: bool=True,
                             xref_masses: bool=True,
                             xref_materials: bool=True,
                             xref_loads: bool=True,
                             xref_constraints: bool=True,
                             xref_aero: bool=True,
                             xref_sets: bool=True,
                             xref_optimization: bool=True,
                             create_superelement_geometry: bool=False,
                             debug=True,
                             word: str='') -> None:
        """
        Performs cross referencing in a way that skips data gracefully.

        .. warning:: not fully implemented
        """
        if not xref:
            return
        self.log.debug('Safe Cross Referencing{word}...')

        xref_obj: CrossReference = self.xref_obj
        if xref_nodes:
            xref_obj.safe_cross_reference_nodes()
            xref_obj.safe_cross_reference_coordinates()
        if xref_elements:
            xref_obj.safe_cross_reference_elements()
        if xref_properties:
            xref_obj.safe_cross_reference_properties()
        if xref_masses:
            xref_obj.safe_cross_reference_masses()
        if xref_materials:
            xref_obj.safe_cross_reference_materials()

        if xref_sets:
            xref_obj.safe_cross_reference_sets()
        if xref_aero:
            xref_obj.safe_cross_reference_aero()
        if xref_constraints:
            xref_obj.safe_cross_reference_constraints()
        if xref_loads:
            xref_obj.safe_cross_reference_loads()
        if xref_optimization:
            xref_obj.safe_cross_reference_optimization()
        if xref_nodes_with_elements:
            xref_obj.cross_reference_nodes_with_elements()

        xref_obj.safe_cross_reference_contact()
        xref_obj.safe_cross_reference_superelements(create_superelement_geometry)

        self.pop_xref_errors()
        for superelement_tuple, superelement in sorted(self.superelement_models.items()):
            if isinstance(superelement_tuple, int):
                word = f' (Superelement {superelement_tuple:d})'
            else:
                wordi, value, label = superelement_tuple
                if label:
                    word = f'BEGIN {wordi}={value:d} LABEL={label}\n'
                else:
                    word = f'BEGIN {wordi}={value:d}\n'

            superelement.safe_cross_reference(
                xref=xref, xref_nodes=xref_nodes, xref_elements=xref_elements,
                xref_nodes_with_elements=xref_nodes_with_elements,
                xref_properties=xref_properties, xref_masses=xref_masses,
                xref_materials=xref_materials, xref_loads=xref_loads,
                xref_constraints=xref_constraints, xref_aero=xref_aero,
                xref_sets=xref_sets, xref_optimization=xref_optimization,
                word=word)

    def safe_empty_nodes(self, nids: list[int],
                         msg: str='') -> tuple[list[GRID], list[int]]:
        """safe xref version of self.Nodes(nid, msg='')"""
        nodes = []
        missing_nodes = []
        for nid in nids:
            try:
                node = self.EmptyNode(nid)
            except KeyError:
                node = nid
                missing_nodes.append(nid)
            nodes.append(node)
        if missing_nodes:
            missing_nodes.sort()
            self.log.warning('Nodes %s are missing%s' % (str(missing_nodes), msg))
        return nodes, missing_nodes

    def safe_get_nodes(self, nids: list[int], msg: str='') -> tuple[list[GRID], str]:
        """safe xref version of self.Nodes(nid, msg='')"""
        nodes = []
        error_nodes = []
        msgi = ''
        for nid in nids:
            try:
                node = self.Node(nid)
            except KeyError:
                error_nodes.append(str(nid))
                node = nid
            nodes.append(nid)
        if error_nodes:
            msgi += 'Could not find nodes %s%s\n' % (', '.join(error_nodes), msg)
        return nodes, msgi

    def safe_get_points(self, point_ids: list[int], msg: str=''):
        """safe xref version of self.Points(point_ids, msg='')"""
        points = []
        error_points = []
        msgi = ''
        for point_id in point_ids:
            try:
                point = self.Point(point_id)
            except KeyError:
                error_points.append(str(point_id))
                point = point_id
            points.append(point)
        if error_points:
            msgi += 'Could not find POINTs %s%s\n' % (', '.join(error_points), msg)
        return points, msgi

    def safe_get_elements(self, eids: list[int], msg: str=''):
        """safe xref version of self.Elements(eid, msg='')"""
        elements = []
        msgi = ''
        for eid in eids:
            try:
                element = self.Element(eid)
            except KeyError:
                element = eid
                msgi += msg % eid
            elements.append(element)
        return elements, msgi

    def safe_element(self, eid: int, ref_id: int, xref_errors: dict[str, Any],
                     msg: str='') -> Element:
        """
        Gets an element card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., a load references an element)

        ref_id = 10 # PLOAD4
        pid = 42  # CQUAD4
        xref_errors = {'eid' : []}
        self.safe_element(eid, ref_id, xref_errors)

        """
        try:
            eid_ref = self.Element(eid, msg=msg)
        except KeyError:
            eid_ref = None
            #self.log.error('cant find Element=%s%s' % (mid, msg))
            xref_errors['eid'].append((ref_id, eid))
        return eid_ref

    def safe_elements(self, eids: list[int], ref_id: int,
                      xref_errors: dict[str, tuple[int, int]],
                      msg: str='') -> list[Element]:
        """
        Gets an series of elements

        Doesn't get rigid (RROD, RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE, RSSCON)
        or mass (CMASS1, CONM2)

        Parameters
        ----------
        ref_id: int
            typically a load_id

        """
        elements = []
        bad_eids = []
        for eid in eids:
            try:
                # elements.append(self.safe_element(eid, ref_id, xref_errors, msg))
                elements.append(self.Element(eid, msg))
            except KeyError:
                bad_eids.append(eid)
                elements.append(None)
                xref_errors['eid'].append((ref_id, eid))
        #if bad_eids:
            #msg = 'eids=%s not found%s.  Allowed elements=%s' % (
                #bad_eids, msg, _unique_keys(self.elements.keys())))
            #self.log.error(msg)
            #raise KeyError(msg)
        return elements

    def safe_property(self, pid: int, ref_id: int, xref_errors,
                      msg: str='') -> Property:
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)

        ref_id = 10 # CQUAD4
        pid = 42  # PSHELL
        xref_errors = {'pid' : []}
        self.safe_property(pid, ref_id, xref_errors)
        """
        try:
            pid_ref = self.Property(pid, msg=msg)
        except KeyError:
            pid_ref = None
            #self.log.error('cant find Property=%s%s' % (mid, msg))
            xref_errors['pid'].append((ref_id, pid))
        return pid_ref

    def safe_property_mass(self, pid, ref_id, xref_errors, msg=''):
        """
        Gets a mass_property card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)
        """
        try:
            pid_ref = self.PropertyMass(pid, msg=msg)
        except KeyError:
            pid_ref = None
            #self.log.error('cant find Property=%s%s' % (mid, msg))
            xref_errors['pid'].append((ref_id, pid))
        return pid_ref

    def safe_material(self, mid: int, ref_id: int,
                      xref_errors, msg='') -> Material:
        """
        Gets a material card

        Parameters
        ----------
        mid : int
            the material_id
        ref_id : int
            the referencing value (e.g., an property references a material, so use self.pid)
        """
        try:
            mid_ref = self.Material(mid, msg=msg)
        except KeyError:
            mid_ref = None
            #self.log.error('cant find Material=%s%s' % (mid, msg))
            xref_errors['mid'].append((ref_id, mid))
        return mid_ref

    def safe_hyperelastic_material(self, mid: int, ref_id: int,
                                   xref_errors, msg='') -> Material:
        """
        Gets a material card

        Parameters
        ----------
        mid : int
            the material_id
        ref_id : int
            the referencing value (e.g., an property references a material, so use self.pid)
        """
        try:
            mid_ref = self.HyperelasticMaterial(mid, msg=msg)
        except KeyError:
            mid_ref = None
            #self.log.error('cant find Material=%s%s' % (mid, msg))
            xref_errors['mid'].append((ref_id, mid))
        return mid_ref

    def safe_coord(self, cid: int, ref_id: int,
                   xref_errors: dict[str, tuple[int, int]], msg: str='') -> Coord:
        """
        Gets a Coord card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an node and element references a coord)

        """
        try:
            cid_ref = self.Coord(cid, msg=msg)
        except KeyError:
            cid_ref = None
            #msgi = 'cant find cid=%s%s' % (cid, msg)
            #self.log.error(msgi)
            xref_errors['cid'].append((ref_id, cid))
        return cid_ref

    def safe_paero(self, paero_id: int, ref_id: int,
                   xref_errors: dict[str, Any], msg: str='') -> PAEROs:
        """
        Gets a PAEROx card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., a load references an element)

        ref_id = 10 # CAERO1
        pid = 42  # PAERO1
        xref_errors = {'paero' : []}
        self.safe_element(pid, ref_id, xref_errors)

        """
        try:
            paero_ref = self.PAero(paero_id, msg=msg)
        except KeyError:
            paero_ref = None
            xref_errors['paero'].append((ref_id, paero_id))
        return paero_ref

    def safe_aefact(self, aefact_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> AEFACT:
        """
        Gets an AEFACT card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an CAERO eid references a AEFACT)

        """
        try:
            aefact_ref = self.AEFact(aefact_id, msg=msg)
        except KeyError:
            aefact_ref = None
            #self.log.error('cant find AFEACT=%s%s' % (aefact_id, msg))
            xref_errors['aefact'].append((ref_id, aefact_id))
        return aefact_ref

    def safe_aelist(self, aelist_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> AELIST:
        """
        Gets an AELIST card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an AESURF eid references a AELIST)

        """
        try:
            aefact_ref = self.AELIST(aelist_id, msg=msg)
        except KeyError:
            aefact_ref = None
            xref_errors['aelist'].append((ref_id, aelist_id))
        return aefact_ref

    def safe_caero(self, caero_id: int, ref_id: int,
                   xref_errors: dict[str, Any], msg: str='') -> CAEROs:
        try:
            caero_ref = self.CAero(caero_id, msg=msg)
        except KeyError:
            caero_ref = None
            xref_errors['caero'].append((ref_id, caero_id))
        return caero_ref

    def safe_tabled(self, tabled_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> TABLEDs:
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an TLOAD1 eid references a TABLED1)
        """
        try:
            tabled_ref = self.TableD(tabled_id, msg=msg)
        except KeyError:
            tabled_ref = None
            xref_errors['tabled'].append((ref_id, tabled_id))
        return tabled_ref

    def safe_tablem(self, table_id: int, ref_id: int,
                    xref_errors: dict[str, tuple[int, int]], msg: str='') -> TABLEMs:
        """
        Gets a Table card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an node and element references a coord)

        """
        try:
            table_ref = self.TableM(table_id, msg=msg)
        except KeyError:
            table_ref = None
            # msgi = 'cant find cid=%s%s' % (cid, msg)
            # self.log.error(msgi)
            xref_errors['tablem'].append((ref_id, table_id))
        return table_ref

    def safe_tableh(self, tableh_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an MATT1 eid references a TABLEH1)
        """
        try:
            tableh_ref = self.TableH(tableh_id, msg=msg)
        except KeyError:
            tableh_ref = None
            xref_errors['tableh'].append((ref_id, tableh_id))
        return tableh_ref

    def safe_desvar(self, desvar_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> DESVAR:
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an DVPREL1 eid references a DESVAR)
        """
        try:
            desvar_ref = self.Desvar(desvar_id, msg=msg)
        except KeyError:
            desvar_ref = None
            xref_errors['desvar'].append((ref_id, desvar_id))
        return desvar_ref
