"""tests the Nastran converters"""
import os
import unittest
import numpy as np
from pathlib import Path
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import read_op2

from pyNastran.converters.format_converter import cmd_line_format_converter
from pyNastran.converters.nastran.nastran_to_cart3d import nastran_to_cart3d, nastran_to_cart3d_filename
from pyNastran.converters.nastran.nastran_to_stl import nastran_to_stl, nastran_to_stl_filename
from pyNastran.converters.nastran.nastran_to_surf import nastran_to_surf, clear_out_solids
from pyNastran.converters.nastran.nastran_to_tecplot import nastran_to_tecplot, nastran_to_tecplot_filename
from pyNastran.converters.nastran.nastran_to_ugrid import nastran_to_ugrid

from pyNastran.converters.tecplot.tecplot_to_nastran import tecplot_to_nastran
from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.converters.cart3d.cart3d import read_cart3d
from pyNastran.bdf.mesh_utils.skin_solid_elements import write_skin_solid_faces

import pyNastran.converters.nastran.nastran_to_ugrid3d
from pyNastran.converters.tecplot.tecplot_to_nastran import nastran_tables_to_tecplot_filenames

from pyNastran.converters.nastran.gui.result_objects.simple_table_results import SimpleTableResults
from pyNastran.converters.nastran.gui.result_objects.layered_table_results import LayeredTableResults
from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2, _composite_method_map

from pyNastran.gui.gui_objects.force_results import ForceResults2
from pyNastran.gui.gui_objects.displacement_results import DisplacementResults2
from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2, DERIVATION_METHODS as shell_derivation_methods
from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'
DIRNAME = Path(os.path.dirname(__file__))
RED = (1., 0., 0.)


class TestNastranGUIObjects(unittest.TestCase):
    """tests:
    - real
     - CompositeStrainStressResults2
     - LayeredTableResults (old)
    - complex
     - DisplacementResults2
     - ForceResults2
    """
    def test_displacement_results(self):
        """tests:
        - complex
         - DisplacementResults2
         - ForceResults2
        """
        bdf_filename = MODEL_PATH / 'aero' / '2_mode_flutter' / '0012_flutter.bdf'
        op2_filename = MODEL_PATH / 'aero' / '2_mode_flutter' / '0012_flutter.op2'
        model = read_bdf(bdf_filename, debug=False)

        out = model.get_displacement_index_xyz_cp_cd()
        icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        node_ids = nid_cp_cd[:, 0]
        xyz_cid0 = xyz_cp

        results = read_op2(op2_filename, debug=None)
        real_mode = (1, 2, 1, 0, 0, '', '')
        complex_mode = (1, 9, 1, 0, 0, '', '')
        subcase_id = 1
        #print(list(results.eigenvectors.keys()))
        disp = results.eigenvectors[complex_mode]
        dxyz = disp #.data

        obj = DisplacementResults2(
            subcase_id,
            node_ids,
            xyz_cid0,
            dxyz,
            title='title',
            t123_offset=0,
            dim_max=1.0,
            data_format='%g',
            is_variable_data_format=False, # ???
            nlabels=None, labelsize=None, ncolors=None,
            colormap='',
            set_max_min=False,
            uname='DisplacementResults2',
        )
        str(obj)
        scale = 1.0
        phase = 0.0
        itime = 0
        res_name = 'asdf'
        xyz, deflected_xyz = obj.get_vector_result_by_scale_phase(
            itime, res_name, scale, phase)
        xyz, deflected_xyz = obj.get_vector_result(itime, res_name, return_dense=True)
        xyz, deflected_xyz = obj.get_vector_result(itime, res_name, return_dense=False)
        assert obj.deflects(itime, res_name) == True, obj.deflects(itime, res_name)
        obj.has_output_checks(itime, res_name)
        assert obj.is_complex == True, obj.is_complex
        xyz, deflected_xyz = obj.get_force_vector_result(itime, res_name)
        mini, maxi = obj.get_default_min_max(itime, res_name)
        imini, imaxi = obj.get_default_min_max(itime, res_name)
        mini, maxi = obj.get_min_max(itime, res_name)
        methods = obj.get_methods(itime, res_name)
        arrow_scale = obj.get_default_arrow_scale(itime, res_name)
        assert np.allclose(arrow_scale, 1.3849555), arrow_scale
        #-------------------------------------------------------------

        methods_txyz_rxyz = [
            'a', 'b', 'c',
            'd', 'e', 'f',
        ]
        #dict[int, tuple[str, str]],
        index_to_base_title_annotation = {
           # 0 : ('abc', 'def'),
            0: {
                'title':' mytitle',
                'corner': 'mycorner',
            },
        }
        force_obj = ForceResults2(
            subcase_id,
            node_ids,
            xyz_cid0,
            dxyz,
            methods_txyz_rxyz=methods_txyz_rxyz,
            index_to_base_title_annotation=index_to_base_title_annotation,
            title='title',
            t123_offset=0,
            dim_max=1.0,
            data_format='%g',
            is_variable_data_format=False, # ???
            nlabels=None, labelsize=None, ncolors=None,
            colormap='',
            set_max_min=False,
            uname='ForceResults2',
        )
        str(force_obj)
        #fxyz, deflected_xyz = force_obj.get_vector_result_by_scale_phase(
        #    itime, res_name, scale, phase)
        fxyz, deflected_xyz = force_obj.get_vector_result(itime, res_name)
        #fxyz, deflected_xyz = force_obj.get_vector_result(itime, res_name, return_dense=False)
        assert force_obj.deflects(itime, res_name) is False, force_obj.deflects(itime, res_name)
        force_obj.has_output_checks(itime, res_name)
        assert force_obj.is_complex is True, force_obj.is_complex
        fxyz, deflected_xyz = force_obj.get_force_vector_result(itime, res_name)
        mini, maxi = force_obj.get_default_min_max(itime, res_name)
        imini, imaxi = force_obj.get_default_min_max(itime, res_name)
        mini, maxi = force_obj.get_min_max(itime, res_name)
        methods = force_obj.get_methods(itime, res_name)
        arrow_scale = force_obj.get_default_arrow_scale(itime, res_name)
        assert np.allclose(arrow_scale, 1.3849555), arrow_scale

        force_obj.get_imin_imax(itime, res_name)

        force_obj.get_default_phase(itime, res_name)
        force_obj.set_phase(itime, res_name, 90.)

        force_obj.set_data_format(itime, res_name, '%g')
        force_obj.get_default_data_format(itime, res_name)

        force_obj.get_nlabels_labelsize_ncolors_colormap(itime, res_name)
        force_obj.set_nlabels_labelsize_ncolors_colormap(
            itime, res_name, nlabels=8, labelsize=8, ncolors=8, colormap='jet')

        force_obj.get_default_legend_title(itime, res_name)
        force_obj.get_annotation(itime, res_name)

        force_obj.get_fringe_result_dense(itime, res_name)
        force_obj.get_fringe_vector_result(itime, res_name)

        force_obj.get_scale(itime, res_name)
        force_obj.get_arrow_scale(itime, res_name)
        force_obj.set_arrow_scale(itime, res_name, 4.0)

    def test_plate_wingbox(self):
        dirname = MODEL_PATH / 'wingbox'
        bdf_filename = dirname / 'wingbox_stitched_together-000.bdf'
        op2_filename = dirname / 'wingbox_stitched_together-000.op2'
        model = read_bdf(bdf_filename, debug=False)
        element_id = np.array(list(model.elements), dtype='int32')

        model_results = read_op2(op2_filename, debug=None)
        subcase_id = 1

        eid_to_nid_map = {}
        for eid, elem in model.elements.items():
            eid_to_nid_map[eid] = elem.nodes

        is_stress = False
        obj2 = PlateStrainStressResults2.load_from_code(
            subcase_id, model, model_results, element_id,
            is_stress, eid_to_nid_map,
            #is_variable_data_format=False,
            require_results=True,
        )
        str(obj2)

        is_stress = True
        obj = PlateStrainStressResults2.load_from_code(
            subcase_id, model, model_results, element_id,
            is_stress, eid_to_nid_map,
            #is_variable_data_format=False,
            require_results=True,
        )
        str(obj)
        assert obj.is_complex == False
        itime = 0
        res_name = (itime, 2, 'header')
        obj.get_fringe_vector_result(itime, res_name)
        obj.get_default_min_max(itime, res_name)
        obj.get_min_max(itime, res_name)
        obj.get_methods(itime, res_name)
        obj.get_phase(itime, res_name)
        obj.get_scale(itime, res_name)
        obj.get_annotation(itime, res_name)
        obj.get_case_flag(itime, res_name)
        obj.get_default_legend_title(itime, res_name)
        obj.get_default_phase(itime, res_name)
        obj.get_default_scale(itime, res_name)
        obj.get_default_arrow_scale(itime, res_name)
        obj.get_case_flag(itime, res_name)
        all_ids, ids = obj.get_location_arrays()

        for method in shell_derivation_methods:
            obj.set_centroid(top_bottom_both='Both',
                             min_max_method=method)
            assert obj.nodal_combine == 'Centroid', obj.nodal_combine
            obj.get_fringe_result(itime, res_name)

        obj.set_centroid(top_bottom_both='Top',
                         min_max_method='Absolute Max')  # doesn't matter
        assert obj.nodal_combine == 'Centroid', obj.nodal_combine
        obj.get_fringe_result(itime, res_name)

        obj.set_centroid(top_bottom_both='Bottom',
                         min_max_method='Absolute Max')  # doesn't matter
        assert obj.nodal_combine == 'Centroid', obj.nodal_combine
        obj.get_fringe_result(itime, res_name)

        obj.set_corner(top_bottom_both='both',
                       min_max_method='absolute max',
                       nodal_combine_method='mean')
        assert obj.nodal_combine == 'Mean', obj.nodal_combine
        obj.get_fringe_result(itime, res_name)

        for method in shell_derivation_methods:
            obj.set_corner(top_bottom_both='Both',
                           min_max_method=method,
                           nodal_combine_method=method)
            obj.get_fringe_result(itime, res_name)

        obj.set_sidebar_args(
             itime, res_name,
             min_max_method='', # Absolute Max
             transform='', # Material
             methods_keys=None,
             # unused
             nodal_combine='', # Centroid
        )

    def test_solid_bending(self):
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.op2')
        model = read_bdf(bdf_filename, debug=False)
        node_id = np.array(list(model.nodes), dtype='int32')
        element_id = np.array(list(model.elements), dtype='int32')

        model_results = read_op2(op2_filename, debug=None)
        subcase_id = 1

        cases = [
            model_results.op2_results.stress.ctetra_stress[subcase_id],
        ]
        result = {2 : 'asdf'}
        title = '???'
        obj = SolidStrainStressResults2(
            subcase_id,
            model,
            node_id,
            element_id,
            cases,
            result,
            title,
            #data_format: str = '%g',
            #is_variable_data_format=False,
            #nlabels=None, labelsize=None, ncolors=None,
            #colormap='',
            #set_max_min=False,
            #uname='SolidStressStrainResults2',
        )
        str(obj)
        is_stress = True
        obj = SolidStrainStressResults2.load_from_code(
            subcase_id, model, model_results, element_id,
            is_stress,
            #is_variable_data_format=False,
            require_results=True,
        )

        assert obj.is_complex == False
        itime = 0
        res_name = (itime, 2, 'header')
        obj.get_fringe_vector_result(itime, res_name)
        obj.get_default_min_max(itime, res_name)
        obj.get_min_max(itime, res_name)
        obj.get_methods(itime, res_name)
        obj.get_phase(itime, res_name)
        obj.get_scale(itime, res_name)
        obj.get_annotation(itime, res_name)
        obj.get_case_flag(itime, res_name)
        obj.get_default_legend_title(itime, res_name)
        obj.get_default_phase(itime, res_name)
        obj.get_default_scale(itime, res_name)
        obj.get_default_arrow_scale(itime, res_name)
        obj.get_case_flag(itime, res_name)
        all_ids, ids = obj.get_location_arrays()
        obj.get_fringe_result(itime, res_name)

    def test_layered_table2(self):
        """tests CompositeStrainStressResults2"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        model = read_bdf(bdf_filename, debug=False)

        #out = model.get_displacement_index_xyz_cp_cd()
        #icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        #node_ids = nid_cp_cd[:, 0]
        #xyz_cid0 = xyz_cp

        model_results = read_op2(op2_filename, debug=None)
        title = 'title'
        subcase_id = 1
        element_id = np.array(list(model.elements))
        result = {
            0: 'aaa',
            1: 'bbb',
        }
        case = model_results.op2_results.strain.cquad4_composite_strain[subcase_id]
        obj = CompositeStrainStressResults2(
            subcase_id,
            model_results,
            element_id, # : np.ndarray
            case, # RealCompositePlateArray
            result, # str
            title, # : str
            data_format='%g',
            is_variable_data_format=False,
            nlabels=None, labelsize=None, ncolors=None,
            colormap='',
            set_max_min=False,
            uname='CompositeStressResults2')
        str(obj)

        itime = 0
        ilayer = 1
        #imethod = 1
        header = 'asdf'
        i = 0
        res_name = (itime, ilayer, header)
        obj.get_fringe_result(i, res_name)

        obj.get_fringe_vector_result(itime, res_name)
        obj.get_default_min_max(itime, res_name)
        obj.get_min_max(itime, res_name)
        obj.get_methods(itime, res_name)
        obj.get_phase(itime, res_name)
        obj.get_scale(itime, res_name)
        obj.get_annotation(itime, res_name)
        obj.get_default_legend_title(itime, res_name)
        obj.get_default_phase(itime, res_name)
        obj.get_default_scale(itime, res_name)
        obj.get_default_arrow_scale(itime, res_name)


    def test_layered_table(self):
        """tests LayeredTableResults"""

        word, method_map = _composite_method_map(is_stress=True)
        word, method_map = _composite_method_map(is_stress=False)

        subcase_id = 1
        eid_max = 5
        nlayers = 2
        headers = ['asdf']
        eids = np.arange(1, eid_max + 1)
        ntimes = 1
        neids = eid_max
        nmethods = 2
        scalars = np.random.random((ntimes, neids, nlayers, nmethods))
        methods = ['a', 'b']
        obj = LayeredTableResults(subcase_id, headers, eids, eid_max, scalars,
                 methods,
                 data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 set_max_min=False, uname='LayeredTableResults')
        str(obj)

        itime = 0
        ilayer = 1
        imethod = 1
        header = 'asdf'
        i = 0
        res_name = (itime, ilayer, imethod, header)
        obj.get_fringe_result(i, res_name)

        obj.get_fringe_vector_result(itime, res_name)
        obj.get_default_min_max(itime, res_name)
        obj.get_min_max(itime, res_name)
        obj.get_methods(itime, res_name)
        obj.get_phase(itime, res_name)
        obj.get_scale(itime, res_name)
        obj.get_annotation(itime, res_name)
        obj.get_default_legend_title(itime, res_name)
        obj.get_default_phase(itime, res_name)
        obj.get_default_scale(itime, res_name)
        obj.get_default_arrow_scale(itime, res_name)
        #ll_ids, ids = obj.get_location_arrays()

class FakeCase:
    def __init__(self, times: np.ndarray):
        self._times = times
        self.headers = ['a', 'b']
        self.title = 'title'
        self.subtitle = 'subtitle'
        ntimes = 4
        nresults = len(self.headers)
        self.data = np.zeros((ntimes, 36, 2))


class TestNastran(unittest.TestCase):
    def test_nastran_to_cart3d_se2(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 2, [1, 2, 3, 4])
        model.add_ctria3(2, 2, [2, 3, 4])

        model.add_cquadr(3, 2, [1, 2, 3, 4])
        model.add_ctriar(4, 2, [2, 3, 4])
        model.add_pshell(2, mid1=100, t=0.1)
        model.add_mat1(100, 3.0e7, None, 0.3)
        model.add_cbar(101, 200, [3, 2], [1., 1., 1.], None)
        model.add_pbarl(200, 100, 'ROD', [1.])
        model.card_count = {
            'GRID': 4,
            'CQUAD4': 1,
            'CQUADR': 1,
            'CTRIA3': 1,
            'CTRIAR': 1,
            'CBAR': 1,
        }
        cart3d = nastran_to_cart3d(model)
        cart3d.flip_model()
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 6
        cart3d.remove_elements([0], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 5
        cart3d.keep_elements([2], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 3
        assert len(cart3d.elements) == 1

        bdf_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.bdf')
        model.write_bdf(bdf_filename)
        cart3d_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.tri')
        nastran_to_cart3d_filename(bdf_filename, cart3d_filename)

    def test_nastran_to_cart3d(self):
        model = BDF(debug=False)
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        model.add_grid(4, [1., 1., 0.])
        model.add_grid(5, [0., 1., 0.])
        model.add_cquad4(10, 2, [2, 3, 4, 5])
        model.add_ctria3(12, 2, [2, 3, 4])

        model.add_cquadr(100, 2, [2, 3, 4, 5])
        model.add_ctriar(121, 2, [2, 3, 4])
        model.add_pshell(2, mid1=100, t=0.1)
        model.add_mat1(100, 3.0e7, None, 0.3)
        model.add_cbar(101, 200, [3, 2], [1., 1., 1.], None)
        model.add_pbarl(200, 100, 'ROD', [1.])
        model.card_count = {
            'GRID': 4,
            'CQUAD4': 1,
            'CQUADR': 1,
            'CTRIA3': 1,
            'CTRIAR': 1,
            'CBAR': 1,
        }
        cart3d = nastran_to_cart3d(model)
        cart3d.flip_model()
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 6
        cart3d.remove_elements([0], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 4
        assert len(cart3d.elements) == 5
        cart3d.keep_elements([2], remove_associated_nodes=True)
        assert len(cart3d.nodes) == 3
        assert len(cart3d.elements) == 1

        bdf_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.bdf')
        bdf_filename2 = os.path.join(DIRNAME, 'nastran_to_cart3d_2.bdf')
        model.write_bdf(bdf_filename)
        cart3d_filename = os.path.join(DIRNAME, 'nastran_to_cart3d.tri')
        nastran_to_cart3d_filename(bdf_filename, cart3d_filename)

        args = ['format_converter', 'nastran', bdf_filename, 'cart3d', cart3d_filename, '--scale', '2.0']
        cmd_line_format_converter(args, quiet=True)

        args = ['format_converter', 'nastran', bdf_filename, 'nastran', bdf_filename2]
        cmd_line_format_converter(args, quiet=True)

    def test_nastran_to_tecplot(self):
        """tests a large number of elements and results in SOL 101"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        tecplot_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.plt')
        tecplot_filename2 = os.path.join(MODEL_PATH, 'elements', 'static_elements2.plt')
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, log=log)
        with self.assertRaises(RuntimeError):
            nastran_to_tecplot(model)
        nastran_to_tecplot_filename(bdf_filename, tecplot_filename, log=log)

        argv = ['format_converter', 'nastran', bdf_filename, 'tecplot', tecplot_filename2]
        with self.assertRaises(RuntimeError):
            cmd_line_format_converter(argv=argv, quiet=True)

    # def test_nastran_to_tecplot_crod(self):
    #     model = BDF(debug=False)
    #     model.add_grid(1, [0., 0., 0.])
    #     model.add_grid(2, [0., 1., 0.])
    #     model.add_crod(10, 100, [1, 2])
    #     tecplot = nastran_to_tecplot(model)
    #     zone = tecplot.zones[0]
    #     assert len(zone.tri_elements) == 0, zone
    #     assert len(zone.quad_elements) == 0, zone
    #     assert len(zone.tet_elements) == 0, zone
    #     assert len(zone.hexa_elements) == 0, zone

    def test_nastran_to_tecplot_tri(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 1., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_ctria3(10, 100, [1, 2, 3])
        tecplot = nastran_to_tecplot(model)
        zone = tecplot.zones[0]
        assert len(zone.tri_elements) == 1, zone
        assert len(zone.quad_elements) == 0, zone
        assert len(zone.tet_elements) == 0, zone
        assert len(zone.hexa_elements) == 0, zone

        bdf_filename = DIRNAME / 'tri.bdf'
        tecplot_to_nastran(tecplot, bdf_filename)
        model = read_bdf(bdf_filename)
        assert len(model.elements) == 1, model.elements
        elem = model.elements[1]
        assert elem.type == 'CTRIA3', elem
        os.remove(bdf_filename)

    def test_nastran_to_tecplot_cquad4(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 1., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_cquad4(10, 100, [1, 2, 3, 4])
        tecplot = nastran_to_tecplot(model)
        zone = tecplot.zones[0]
        assert len(zone.tri_elements) == 0, zone
        assert len(zone.quad_elements) == 1, zone
        assert len(zone.tet_elements) == 0, zone
        assert len(zone.hexa_elements) == 0, zone

        bdf_filename = DIRNAME / 'quad.bdf'
        tecplot_to_nastran(tecplot, bdf_filename)
        model = read_bdf(bdf_filename)
        assert len(model.elements) == 1, model.elements
        elem = model.elements[1]
        assert elem.type == 'CQUAD4', elem
        os.remove(bdf_filename)

    def test_nastran_to_tecplot_cpenta(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 1., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 0., 1.])
        model.add_grid(5, [0., 1., 1.])
        model.add_grid(6, [1., 1., 1.])
        model.add_cpenta(10, 100, [1, 2, 3, 4, 5, 6])
        tecplot = nastran_to_tecplot(model)
        zone = tecplot.zones[0]
        assert len(zone.tri_elements) == 0, zone
        assert len(zone.quad_elements) == 0, zone
        assert len(zone.tet_elements) == 0, zone
        assert len(zone.hexa_elements) == 1, zone

        bdf_filename = DIRNAME / 'penta.bdf'
        tecplot_to_nastran(tecplot, bdf_filename)
        model = read_bdf(bdf_filename)
        assert len(model.elements) == 1, model.elements
        elem = model.elements[1]
        assert elem.type == 'CPENTA', elem
        os.remove(bdf_filename)

    def test_nastran_to_tecplot_chexa(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 1., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [0., 1., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [1., 0., 1.])
        model.add_chexa(10, 100, [1, 2, 3, 4, 5, 6, 7, 8])

        bdf_filename = DIRNAME / 'hexa.bdf'
        tecplot = nastran_to_tecplot(model)
        zone = tecplot.zones[0]
        assert len(zone.tri_elements) == 0, zone
        assert len(zone.quad_elements) == 0, zone
        assert len(zone.tet_elements) == 0, zone
        assert len(zone.hexa_elements) == 1, zone
        tecplot_to_nastran(tecplot, bdf_filename)
        model = read_bdf(bdf_filename)
        assert len(model.elements) == 1, model.elements
        elem = model.elements[1]
        assert elem.type == 'CHEXA', elem
        os.remove(bdf_filename)

    def test_nastran_to_tecplot_chexa_cpenta(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 1., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [0., 1., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [1., 0., 1.])
        model.add_cpenta(11, 100, [1, 2, 3, 4, 5, 6])
        model.add_chexa(10, 100, [1, 2, 3, 4, 5, 6, 7, 8])

        tecplot = nastran_to_tecplot(model)
        zone = tecplot.zones[0]
        assert len(zone.tri_elements) == 0, zone
        assert len(zone.quad_elements) == 0, zone
        assert len(zone.tet_elements) == 0, zone
        assert len(zone.hexa_elements) == 2, zone

        bdf_filename = DIRNAME / 'penta_hexa.bdf'
        tecplot_to_nastran(tecplot, bdf_filename)
        model = read_bdf(bdf_filename)
        assert len(model.elements) == 2, model.elements
        elem1 = model.elements[1]
        elem2 = model.elements[2]
        assert elem1.type == 'CPENTA', elem1
        assert elem2.type == 'CHEXA', elem2
        os.remove(bdf_filename)

    def test_nastran_to_ugrid_01(self):
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')

        size = 8
        debug = False
        log = SimpleLogger(level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, log=log, debug=debug)
        #log = model.log
        #model.get_element_faces()
        skin_bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin.bdf')
        write_skin_solid_faces(model, skin_bdf_filename, write_solids=True,
                               write_shells=True,
                               size=size, is_double=False, encoding=None)

        bdf_model = read_bdf(skin_bdf_filename, log=log, debug=debug)
        ugrid_filename_out = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin.b8.ugrid')
        ugrid_filename_out2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin2.b8.ugrid')
        nastran_to_ugrid(bdf_model, ugrid_filename_out, properties=None,
                         check_shells=True, check_solids=True)
        ugrid = read_ugrid(ugrid_filename_out, encoding=None, log=log,
                           debug=debug)

        skin_bdf_filename2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin2.bdf')
        skin_cart3d_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin2.tri')
        skin_cart3d_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin3.tri')
        skin_stl_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_skin3.stl')

        #msg += "  format_converter nastran   <INPUT> <format2> <OUTPUT> [-o <OP2>] --no_xref\n"
        #msg += "  format_converter <format1> <INPUT> tecplot   <OUTPUT> [-r RESTYPE...] [-b] [--block] [-x <X>] [-y <Y>] [-z <Z>] [--scale SCALE]\n"
        #msg += "  format_converter <format1> <INPUT> stl       <OUTPUT> [-b]  [--scale SCALE]\n"
        #msg += "  format_converter cart3d    <INPUT> <format2> <OUTPUT> [-b]  [--scale SCALE]\n"
        #msg += "  format_converter <format1> <INPUT> <format2> <OUTPUT> [--scale SCALE]\n"
        argv = ['format_converter', 'nastran', bdf_filename, 'ugrid', ugrid_filename_out2]
        with self.assertRaises(RuntimeError):
            cmd_line_format_converter(argv=argv, quiet=True)

        #argv = ['format_converter', 'nastran', bdf_filename, 'cart3d', skin_cart3d_filename3]
        #cmd_line_format_converter(argv=argv)

        #argv = ['format_converter', 'nastran', bdf_filename, 'stl', skin_stl_filename3]
        #cmd_line_format_converter(argv=argv)

        ugrid.write_bdf(skin_bdf_filename2, include_shells=True, include_solids=True,
                        convert_pyram_to_penta=True, encoding=None,
                        size=size, is_double=False)
        read_bdf(skin_bdf_filename2, log=log, debug=debug)

        with self.assertRaises(NotImplementedError):
            nastran_to_cart3d_filename(skin_bdf_filename2, skin_cart3d_filename)

        ugrid.write_bdf(skin_bdf_filename2, include_shells=True, include_solids=False,
                        convert_pyram_to_penta=True, encoding=None,
                        size=size, is_double=False)

        nastran_to_cart3d_filename(skin_bdf_filename2, skin_cart3d_filename)
        read_cart3d(skin_cart3d_filename, log=log)

        os.remove(ugrid_filename_out)
        os.remove(skin_bdf_filename)
        os.remove(skin_bdf_filename2)
        os.remove(skin_cart3d_filename)

    def test_nastran_to_stl(self):
        """tests nastran_to_stl"""
        bdf_filename = MODEL_PATH / 'plate' / 'plate.bdf'
        stl_filename = MODEL_PATH / 'plate' / 'plate.stl'
        log = SimpleLogger(level='warning', encoding='utf-8')
        nastran_to_stl_filename(bdf_filename, stl_filename, is_binary=False, log=log)

        model = BDF(debug=None)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_mat1(1000, 3.0e7, None, 0.3)
        model.add_pshell(100, 1000, 0.1)
        model.add_conrod(10, 1000, [1, 2])
        model.add_ctria3(1, 100, [1, 2, 3])
        model.add_ctriar(2, 100, [1, 2, 3])
        model.add_cquad4(3, 100, [1, 2, 3, 4])
        model.add_cquadr(4, 100, [1, 2, 3, 4])
        stl_filename2 = MODEL_PATH / 'plate' / 'spike.stl'
        nastran_to_stl(model, stl_filename2)

    def test_nastran_to_tecplot_case(self):
        bdf_filename = MODEL_PATH / 'plate' / 'plate.bdf'
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_model = read_bdf(bdf_filename, log=log, debug=False)

        times = np.arange(2)
        case = FakeCase(times)
        tecplot_filename_base = 'cat%d'
        nastran_tables_to_tecplot_filenames(tecplot_filename_base, bdf_model, case,
                                            variables=None, ivars=None)

    def test_format_converter(self):
        """tests nastran_to_stl"""
        bdf_filename = os.path.join(MODEL_PATH, 'plate', 'plate.bdf')
        bdf_filename2 = os.path.join(MODEL_PATH, 'plate', 'plate2.bdf')

        stl_filename = os.path.join(MODEL_PATH, 'plate', 'plate.stl')
        ugrid_filename = os.path.join(MODEL_PATH, 'plate', 'plate.b8.ugrid')
        cart3d_filename = os.path.join(MODEL_PATH, 'plate', 'plate.tri')
        tecplot_filename = os.path.join(MODEL_PATH, 'plate', 'plate.plt')

        argv = ['format_converter', 'nastran', bdf_filename, 'stl', stl_filename]
        cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'nastran', bdf_filename, 'tecplot', tecplot_filename]
        cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'nastran', bdf_filename, 'ugrid', ugrid_filename]
        with self.assertRaises(RuntimeError):
            cmd_line_format_converter(argv=argv, quiet=True)


        #argv = ['format_converter', 'nastran', bdf_filename, 'cart3d', cart3d_filename]
        #cmd_line_format_converter(argv=argv, quiet=True)
        #os.remove(stl_filename)
        #os.remove(cart3d_filename)
        #os.remove(tecplot_filename)
        # -------------------------
        tecplot_filename2 = os.path.join(MODEL_PATH, 'plate', 'plate2.plt')
        argv = ['format_converter', 'stl', stl_filename, 'nastran', bdf_filename2]
        cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'stl', stl_filename, 'tecplot', tecplot_filename2]
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'stl', tecplot_filename, 'ugrid', ugrid_filename]
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(bdf_filename2)
        #os.remove(stl_filename)
        #os.remove(cart3d_filename)
        #os.remove(tecplot_filename)
        # -------------------------
        argv = ['format_converter', 'tecplot', tecplot_filename, 'nastran', bdf_filename2]
        cmd_line_format_converter(argv=argv, quiet=True)

        #argv = ['format_converter', 'tecplot', tecplot_filename, 'stl', stl_filename]
        #cmd_line_format_converter(argv=argv, quiet=True)

        argv = ['format_converter', 'tecplot', tecplot_filename, 'ugrid', ugrid_filename]
        with self.assertRaises(AssertionError):
            cmd_line_format_converter(argv=argv, quiet=True)

        os.remove(bdf_filename2)
        os.remove(stl_filename)
        #os.remove(cart3d_filename)
        os.remove(tecplot_filename)

    def test_clear_out_solids(self):
        """tests clear_out_solids"""
        deck = (
            "$ pyNastran: punch=True\n"
            "GRID,1\n"
            "GRID,2\n"
            "GRID,3\n"
            "GRID,4\n"
            "GRID,5\n"
            "GRID,6\n"
            "GRID,7\n"
            "GRID,8\n"

            "GRID,9\n"
            "GRID,10\n"
            "GRID,11\n"
            "GRID,12\n"
            "CHEXA,1,1, 5,6,7,8,9,10,\n"
            ",7,8\n"
            "CQUAD4,2,200, 1,2,3,4\n"
            # doesn't work
            #"CHEXA,1,1, 1,2,3,4,5,6,\n"
            #",7,8\n"
            #"CQUAD4,2,200, 8,9,10,11\n"
            "PSHELL,200,1000,0.1\n"
            "PSOLID,100,1000\n"
            "MAT1,1000,3.0e7,,0.3\n"
        )

        bdf_filename = 'deck.bdf'
        bdf_clean_filename = 'clean.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(deck)

        log = SimpleLogger(level='warning', encoding='utf-8')
        model = read_bdf(bdf_filename, xref=False, log=log)
        clear_out_solids(model, bdf_clean_filename, renumber=True,
                         equivalence=False, equivalence_tol=0.01)

        model = read_bdf(bdf_clean_filename, log=log)
        assert len(model.nodes) == 4, len(model.nodes)
        assert len(model.elements) == 1, len(model.elements)
        os.remove(bdf_filename)
        os.remove(bdf_clean_filename)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
