"""
Defines:
 - SURF tests
"""
from __future__ import print_function
import os
import unittest

import pyNastran
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.converters.nastran.nastran_to_surf import nastran_to_surf, read_bdf
from pyNastran.converters.aflr.surf.surf_io import SurfIO
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.utils.log import get_logger

PKG_PATH = pyNastran.__path__[0]
model_path = os.path.join(PKG_PATH, 'converters', 'tecplot', 'models')
nastran_path = os.path.join(PKG_PATH, '..', 'models')


def remap_cards(model,
                remap_nodes=True, remap_elements=True, remap_properties=True,
                remap_materials=True):
    model.uncross_reference()

    if remap_nodes:
        nodes = {}
        for node in model.nodes.values():
            nid = node.nid
            nodes[nid] = node
        model.nodes = nodes

    if remap_elements:
        elements = {}
        for element in model.elements.values():
            eid = element.eid
            elements[eid] = element
        model.elements = elements

    if remap_properties:
        properties = {}
        for prop in model.properties.values():
            pid = prop.pid
            properties[pid] = prop
        model.properties = properties

    if remap_materials:
        materials = {}
        for material in model.materials.values():
            mid = material.mid
            materials[mid] = material
        model.materials = materials
    model.cross_reference()

def delete_properties(bdf_model, property_types_to_save=None):
    pids_to_delete = set([])
    if property_types_to_save:
        for pid, prop in bdf_model.properties.items():
            ptype = prop.type
            if ptype not in property_types_to_save:
                pids_to_delete.add(pid)
                bdf_model._type_to_id_map[ptype].remove(pid)
    for pid in pids_to_delete:
        del bdf_model.properties[pid]

def delete_elements(bdf_model, element_types_to_save=None):
    eids_to_delete = set([])
    if element_types_to_save:
        for eid, element in bdf_model.elements.items():
            etype = element.type
            if etype not in element_types_to_save:
                eids_to_delete.add(eid)
                bdf_model._type_to_id_map[etype].remove(eid)
    for eid in eids_to_delete:
        del bdf_model.elements[eid]

#def delete_forces(bdf_model, eids_to_delete=None):
    #eids_to_delete = set([])
    #loads = {}
    #for load_id, load_set in bdf_model.loads.items():
        #load_set = []
        #for load in load_set:
            #print(load)
    #if element_types_to_save:
        #for eid, element in bdf_model.elements.items():
            #if element.type not in element_types_to_save:
                #eids_to_delete.add(eid)
    #for eid in eids_to_delete:
        #del bdf_model.elements[eid]

class SurfGui(SurfIO, FakeGUIMethods):
    """defines the UGRID 2D/3D interface"""
    def __init__(self):
        FakeGUIMethods.__init__(self)
        SurfIO.__init__(self, self)
        self.build_fmts(['surf'], stop_on_failure=True)

class TestSurfGui(unittest.TestCase):
    """defines *.surf tests"""
    def test_surf_gui_01(self):
        """tests two_blade_wake_sym_extended.surf"""
        ugrid_filename = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models',
                                      'two_blade_wake_sym_extended.surf')
        log = get_logger(level='warning')
        test = SurfGui()
        test.log = log
        test.on_load_geometry(ugrid_filename, geometry_format='surf', raise_error=True)
        #test.load_surf_geometry(ugrid_filename)

    def test_surf_01(self):
        """tests two_blade_wake_sym_extended.surf"""
        MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
        bdf_filename = os.path.join(MODEL_PATH, 'iSat', 'ISat_Launch_Sm_Rgd.dat')
        surf_filename = os.path.join(MODEL_PATH, 'iSat', 'ISat_Launch_Sm_Rgd.surf')
        bdf_model = read_bdf(bdf_filename)

        #ugrid_filename = os.path.join(PKG_PATH, 'converters', 'aflr', 'ugrid', 'models',
                                      #'two_blade_wake_sym_extended.surf')
        #log = get_logger(level='warning')

        pid_to_element_flags = {}
        for pid, prop in bdf_model.properties.items():
            if prop.type in ['PSHELL', 'PCOMP']:
                # name, initial_normal_spacing, bl_thickness, grid_bc
                pid_to_element_flags[pid] = ['na;me', 0.01, 0.1, 1]

        with self.assertRaises(RuntimeError):
            nastran_to_surf(bdf_model, pid_to_element_flags, surf_filename,
                            renumber_pids=None,
                            line_map=None, scale=1.0,
                            tol=1e-10, xref=True)

        delete_elements(
            bdf_model,
            element_types_to_save=['CTRIA3', 'CQUAD4'])
        delete_properties(
            bdf_model,
            property_types_to_save=['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'])
        #print(bdf_model.properties)

        bdf_model.uncross_reference()
        remove_unused(bdf_model, remove_nids=True, remove_cids=True,
                      remove_pids=True, remove_mids=True)
        #delete_forces(bdf_model)
        bdf_model.case_control_deck = None

        #bdf_filename_re = os.path.join(MODEL_PATH, 'iSat', 'ISat_Launch_Sm_Rgd_re.dat')
        bdf_filename_re = None
        bdf_model.cross_reference()
        bdf_model_re = bdf_renumber(
            bdf_model, bdf_filename_re,
            #size=8, is_double=False,
            #starting_id_dict=None,
            round_ids=False, cards_to_skip=None,
            log=bdf_model.log, debug=False)[0]

        remap_cards(bdf_model_re)
        #print(bdf_model_re.properties)
        #print(bdf_model_re.elements)
        #aaa

        #bdf_model_re = read_bdf(bdf_filename_re)
        #print(bdf_model_re.get_bdf_stats())

        pid_to_element_flags = {}
        for pid, prop in bdf_model_re.properties.items():
            if prop.type in ['PSHELL', 'PCOMP']:
                # name, initial_normal_spacing, bl_thickness, grid_bc
                pid_to_element_flags[pid] = ['na;me', 0.01, 0.1, 1]

        nastran_to_surf(bdf_model_re,
                        pid_to_element_flags, surf_filename,
                        renumber_pids=None,
                        line_map=None, scale=1.0,
                        tol=1e-10, xref=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
