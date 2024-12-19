"""
Creates safe cross referencing

Safe cross-referencing skips failed xref's

"""
from __future__ import annotations
from typing import Callable, Any, TYPE_CHECKING

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
        self.log.debug(f'Safe Cross Referencing{word}...')

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
        return _safe_attr(self, eid, ref_id, xref_errors,
                          func=self.Element, word='eid', msg=msg)

    def safe_nodes(self, nids: list[int], ref_id: int,
                   xref_errors: dict[str, tuple[int, int]],
                   msg: str='') -> list[Element]:
        """
        Gets an series of GRIDs/SPOINTs/EPOINTs

        Parameters
        ----------
        ref_id: int
            typically an element_id

        """
        return _safe_attrs(self, nids, ref_id, xref_errors,
                          func=self.Node, word='nid', msg=msg)

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
        return _safe_attrs(self, eids, ref_id, xref_errors,
                          func=self.Element, word='eid', msg=msg)
    def safe_node(self, nid: int, ref_id: int, xref_errors,
                  msg: str='') -> Property:
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)

        ref_id = 10 # CQUAD4
        nid = 42
        xref_errors = {'nid' : []}
        self.safe_node(pid, ref_id, xref_errors)
        """
        return _safe_attr(self, nid, ref_id, xref_errors,
                          func=self.Node, word='nid', msg=msg)

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
        return _safe_attr(self, pid, ref_id, xref_errors,
                          func=self.Property, word='pid', msg=msg)

    def safe_property_mass(self, pid, ref_id, xref_errors, msg=''):
        """
        Gets a mass_property card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an element references a property)
        """
        return _safe_attr(self, pid, ref_id, xref_errors,
                          func=self.PropertyMass, word='pid', msg=msg)

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
        return _safe_attr(self, mid, ref_id, xref_errors,
                          func=self.Material, word='mid', msg=msg)

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
        return _safe_attr(self, mid, ref_id, xref_errors,
                          func=self.HyperelasticMaterial, word='mid', msg=msg)

    def safe_coord(self, cid: int, ref_id: int,
                   xref_errors: dict[str, tuple[int, int]], msg: str='') -> Coord:
        """
        Gets a Coord card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an node and element references a coord)

        """
        return _safe_attr(self, cid, ref_id, xref_errors,
                          func=self.Coord, word='cid', msg=msg)

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
        return _safe_attr(self, paero_id, ref_id, xref_errors,
                          func=self.PAero, word='paero', msg=msg)

    def safe_aefact(self, aefact_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> AEFACT:
        """
        Gets an AEFACT card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an CAERO eid references a AEFACT)

        """
        return _safe_attr(self, aefact_id, ref_id, xref_errors,
                          func=self.AEFact, word='aefact', msg=msg)

    def safe_aelist(self, aelist_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> AELIST:
        """
        Gets an AELIST card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an AESURF eid references a AELIST)

        """
        return _safe_attr(self, aelist_id, ref_id, xref_errors,
                          func=self.AELIST, word='aelist', msg=msg)
        return aefact_ref

    def safe_caero(self, caero_id: int, ref_id: int,
                   xref_errors: dict[str, Any], msg: str='') -> CAEROs:
        return _safe_attr(self, caero_id, ref_id, xref_errors,
                          func=self.CAero, word='caero', msg=msg)

    def safe_tabled(self, tabled_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> TABLEDs:
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an TLOAD1 eid references a TABLED1)
        """
        return _safe_attr(self, tabled_id, ref_id, xref_errors,
                          func=self.TableD, word='tabled', msg=msg)

    def safe_tablem(self, tablem_id: int, ref_id: int,
                    xref_errors: dict[str, tuple[int, int]], msg: str='') -> TABLEMs:
        """
        Gets a Table card

        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an node and element references a coord)

        """
        return _safe_attr(self, tablem_id, ref_id, xref_errors,
                          func=self.TableM, word='tablem', msg=msg)

    def safe_tableh(self, tableh_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str=''):
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an MATT1 eid references a TABLEH1)
        """
        return _safe_attr(self, tableh_id, ref_id, xref_errors,
                          func=self.TableH, word='tableh', msg=msg)

    def safe_desvar(self, desvar_id: int, ref_id: int,
                    xref_errors: dict[str, Any], msg: str='') -> DESVAR:
        """
        Parameters
        ----------
        ref_id : int
            the referencing value (e.g., an DVPREL1 eid references a DESVAR)
        """
        return _safe_attr(self, desvar_id, ref_id, xref_errors,
                          func=self.Desvar, word='desvar', msg=msg)

def _safe_attr(model: BDF, idi: int, ref_id: int,
               xref_errors: dict[str, Any],
               func: Callable,
               word: str, msg: str):
    try:
        id_ref = func(idi, msg=msg)
    except KeyError:
        id_ref = None
        # self.log.error('cant find Property=%s%s' % (mid, msg))
        xref_errors[word].append((ref_id, idi))
    return id_ref

def _safe_attrs(model: BDF, ids: list[int], ref_id: int,
                xref_errors: dict[str, Any],
                func: Callable,
                word: str, msg: str):
    elements = []
    bad_eids = []
    for eid in ids:
        try:
            # elements.append(self.safe_element(eid, ref_id, xref_errors, msg))
            elements.append(func(eid, msg))
        except KeyError:
            bad_eids.append(eid)
            elements.append(None)
            xref_errors[word].append((ref_id, eid))
    #if bad_eids:
        #msg = 'eids=%s not found%s.  Allowed elements=%s' % (
            #bad_eids, msg, _unique_keys(self.elements.keys())))
        #self.log.error(msg)
        #raise KeyError(msg)
    return elements
