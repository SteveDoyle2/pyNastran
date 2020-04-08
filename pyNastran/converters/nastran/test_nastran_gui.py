"""tests the NastranIO class"""
import os
from copy import deepcopy
import unittest

import numpy as np
try:
    import matplotlib
    matplotlib.use('Agg')
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pyparsing is missing
    IS_MATPLOTLIB = False
#except ImportError:
    #pass
import vtk

from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.test_aero import get_zona_model
from pyNastran.bdf.errors import DuplicateIDsError
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.gui.nastran_io import NastranIO
RED = (1., 0., 0.)


class NastranGUI(NastranIO, FakeGUIMethods):
    def __init__(self, inputs=None):
        FakeGUIMethods.__init__(self, inputs=inputs)
        NastranIO.__init__(self)
        self.build_fmts(['nastran'], stop_on_failure=True)

PKG_PATH = pyNastran.__path__[0]
STL_PATH = os.path.join(PKG_PATH, 'converters', 'stl')
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestNastranGUI(unittest.TestCase):

    def test_settings(self):
        from qtpy import QtCore
        settings = QtCore.QSettings()
        test = NastranGUI()
        is_loaded = test.settings.load(settings)
        assert is_loaded is True
        test.settings.save(settings, is_testing=True)
        is_loaded = test.settings.load(settings)
        assert is_loaded is True

        test.settings.set_annotation_size_color(size=10, color=None)
        #test.settings.set_annotation_size_color(size=10, color=RED)

        test.settings.set_coord_scale(2.0, render=True)
        test.settings.set_coord_text_scale(10, render=True)

        test.settings.update_coord_scale(coord_scale=None, render=True)
        test.settings.set_background_color_to_white(render=True)

        color = RED
        opacity = 0.4
        test.settings.set_background_color(color, render=True)
        test.settings.set_background_color2(color, render=True)
        test.settings.set_highlight_color(color)
        test.settings.set_highlight_opacity(opacity)
        test.settings.set_text_color(color, render=True)
        test.settings.set_text_size(10)
        test.settings.set_magnify(magnify=4)
        #self.settings.s

    def test_solid_shell_bar_obj(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        obj_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.obj')
        model = BDF()
        model.read_bdf(bdf_filename)
        model.save(obj_filename, unxref=True)

        test = NastranGUI()
        test.load_nastran_geometry(obj_filename)

    @unittest.skipIf(IS_MATPLOTLIB is False, 'No matplotlib')
    def test_solid_shell_bar_01(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)
        test.cycle_results()
        test.on_rcycle_results()

        #print('test.result_cases', test.result_cases)
        #gpforce = test.model.grid_point_forces[1]

        icase_gpforce = None
        for icase, (case, dummy) in test.result_cases.items():
            if hasattr(case, 'gpforce_array'):
                icase_gpforce = icase
                break
        else:
            raise RuntimeError('missing gpforce')

        case, (unused_i, unused_name) = test.result_cases[icase_gpforce]
        str(case)
        gpforce = case.gpforce_array
        model_name = 'main'


        p1 = [0., 0., 0.]
        p3 = [1., 0., 0.]

        p2 = [0., 1., 0.]
        zaxis = [0., 0., 1.]
        with self.assertRaises(ValueError):
            test.shear_moment_torque_obj.plot_shear_moment_torque(
                model_name, gpforce,
                p1, p2, p3, zaxis,
                method='Z-Axis Projection',
                cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
                nplanes=5, plane_color=None, plane_opacity=0.5,
                csv_filename=None, show=False, stop_on_failure=True)

        p1 = np.array([0., 0., 0.]) # origin
        p2 = np.array([1., 0., 0.]) # xaxis
        p3 = np.array([1., 0., 0.]) # end
        zaxis = np.array([0., 0., 1.])
        #idir = 0
        test.shear_moment_torque_obj.plot_shear_moment_torque(
            model_name, gpforce,
            p1, p2, p3, zaxis,
            method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
            nplanes=5, plane_color=None, plane_opacity=0.5,
            csv_filename=None, show=False, stop_on_failure=True)

        with self.assertRaises(TypeError):
            # we need to set the case to a grid point force result
            test.cutting_plane_obj.make_cutting_plane(
                model_name,
                p1, p2, zaxis,
                method='Z-Axis Projection',
                cid_p1=0, cid_p2=0, cid_zaxis=0,
                ytol=1., plane_atol=1e-5,
                plane_color=None, plane_opacity=0.5,
                csv_filename=None, show=False, stop_on_failure=True)

        # setting the case to a grid point force result
        test.icase_fringe = icase_gpforce
        test._cycle_results(icase_gpforce)
        test.cutting_plane_obj.make_cutting_plane(
            model_name,
            p1, p2, zaxis,
            method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_zaxis=0,
            ytol=1., plane_atol=1e-5,
            plane_color=None, plane_opacity=0.5,
            csv_filename=None, show=False, stop_on_failure=True)

        test.icase_fringe = 0
        #with self.assertRaises(RuntimeError):
        test.cutting_plane_obj.make_cutting_plane(
            model_name,
            p1, p2, zaxis,
            method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_zaxis=0,
            ytol=1., plane_atol=1e-5,
            plane_color=None, plane_opacity=0.5,
            csv_filename=None, show=False, stop_on_failure=True)

    def test_solid_shell_bar_02(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'mode_solid_shell_bar.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'mode_solid_shell_bar.op2')

        test = NastranGUI()
        test.legend_obj.set_legend_menu()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)
        assert len(test.models['main'].elements) > 0

        test.on_rcycle_results()
        test.on_update_legend(
            title='Title', min_value=0., max_value=1.,
            scale=0.0, phase=0.0,
            arrow_scale=1.,
            data_format='%.0f',
            is_low_to_high=True, is_discrete=True, is_horizontal=True,
            nlabels=None, labelsize=None, ncolors=None, colormap=None,
            is_shown=True, render=True)
        test.on_update_legend(
            title='Title', min_value=0., max_value=1.,
            scale=0.0, phase=0.0,
            arrow_scale=1.,
            data_format='%.0f',
            is_low_to_high=True, is_discrete=True, is_horizontal=False,
            nlabels=None, labelsize=None, ncolors=None, colormap='viridis',
            is_shown=True, render=True)
        test.legend_obj.set_legend_menu()
        test.on_set_camera_data(
            {'distance': 15.23729238729831,
             'prallel_proj': None,
             'view_angle': 30.0,
             'parallel_scale': 3.9437014656284517,
             'position': (-8.279127062822164, 4.306812025814127, 11.191236382055052),
             'view_up': (0.14388395111701072, 0.9587296714789404, -0.245224031523912),
             'clip_range': (7.44295814719721, 25.085506595796954),
             'focal_point': (0.49249999999999994, 0.0, -0.5)}
        )
        test.settings.reset_settings()
        test.on_set_font_size(8)
        test.on_increase_font_size()
        test.on_decrease_font_size()

        labels_list = []
        text = 'text'
        x, y, z = 0., 0., 0.
        labels_list.append(test.create_annotation(text, x, y, z))

        cell_id = 1
        world_position = [0., 0., 1.]
        res_name, result_values, xyz = test.get_result_by_cell_id(
            cell_id, world_position,
            icase=0)
        assert res_name == 'NodeID', 'res_name=%r' % res_name
        assert result_values == 2, 'result_values=%r' % result_values
        assert isinstance(xyz, list), xyz

        #node_xyz = None
        cell_id = 5
        #out = test.mark_actions.get_result_by_xyz_cell_id(node_xyz, cell_id)
        #result_name, result_values, node_id, xyz = out

        eids = [1, 2]
        icase_result = 2
        icase_to_apply = 3
        test.label_actors[2] = []
        test.label_actors[3] = []
        test.mark_actions.mark_elements_by_different_case(
            eids, icase_result, icase_to_apply, stop_on_failure=True, )

        #eids = [1, 2]
        with self.assertRaises(NotImplementedError):
            test.mark_actions.highlight_elements(eids, model_name='main')

        nids = [1, 2]
        icase = 1
        test.label_actors[1] = []
        text = 'cat'
        test.mark_actions.mark_nodes(nids, icase, text)

        with self.assertRaises(RuntimeError): # icase_to_apply=166 doesn't exist
            test.mark_elements(eids, stop_on_failure=True, show_command=True)
        #test.mark_elements_by_case(eids, stop_on_failure=True, show_command=True)

        test.icase = 2 # PropertyID
        test.mark_elements(eids, stop_on_failure=True, show_command=True)
        test.mark_elements_by_case(eids, stop_on_failure=True, show_command=True)

        #for key, obj in test.result_cases.items():
            #print(key)
            #print(obj)

        # fail mapping strain energy because we're on NodeID
        test.icase = 0 # NodeID
        test.icase_fringe = 0 # NodeID
        is_passed = test.map_element_centroid_to_node_fringe_result(update_limits=True, show_msg=True)

        obj, (itime, name) = test.result_cases[test.icase]
        str(obj)
        assert is_passed == False, f'map_element_centroid_to_node_fringe_result should fail for NodeID\n{obj}'

        # map strain energy
        keys = list(test.result_cases.keys())
        assert len(keys) == 689, len(keys)
        icase = keys[-1]
        obj, (itime, name) = test.result_cases[icase]
        test.icase_fringe = icase
        str(obj)

        title = obj.get_title(itime, name)
        assert title == 'Strain Energy Density', str(obj)
        is_passed = test.map_element_centroid_to_node_fringe_result(update_limits=True, show_msg=False)
        assert is_passed == True, 'map_element_centroid_to_node_fringe_result failed'

    def test_solid_shell_bar_02b(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'mode_solid_shell_bar.bdf')

        test = NastranGUI()
        test.on_load_geometry(infile_name=bdf_filename, geometry_format='nastran', name='main',
                              plot=True, raise_error=True)

    def test_solid_shell_bar_03(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'buckling_solid_shell_bar.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'buckling_solid_shell_bar.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_solid_bending(self):
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        #op2_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_ogs.op2')
        deflection_filename1 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_deflection_node.txt')
        deflection_filename2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_deflection_node_short.txt')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        #test.load_nastran_results(op2_filename)
        nresult_cases = len(test.result_cases)
        icase = max(test.result_cases)

        # these are the two cases we're checking were added
        test.on_load_custom_results(out_filename=deflection_filename1, restype='Deflection')
        test.on_load_custom_results(out_filename=deflection_filename1, restype='Force')
        dresult_cases = len(test.result_cases) - nresult_cases
        icase_final = max(test.result_cases)
        dcase = icase_final - icase
        assert dresult_cases == 2, dresult_cases
        assert dcase == 2, dcase
        assert (icase_final - 1) in test.label_actors
        assert icase_final in test.label_actors
        assert len(test.label_actors[icase_final]) == 0

        nids = [1, 2, 3, 5]
        icase = icase_final
        text = 'word'
        test.mark_nodes(nids, icase, text)
        assert len(test.label_actors[icase_final]) == 4, len(test.label_actors[icase_final])

        # test nodal results
        #'node', 'element', 'deflection', 'force', 'patran_nod',
        csv_filename1 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node.csv')
        csv_filename2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node_extra.txt')
        csv_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node_short.txt')

        # missing/extra nodes
        test.on_load_custom_results(out_filename=csv_filename1, restype='node', stop_on_failure=True)
        test.on_load_custom_results(out_filename=csv_filename2, restype='node', stop_on_failure=True)
        test.on_load_custom_results(out_filename=csv_filename3, restype='node', stop_on_failure=True)

        # missing nodes
        test.on_load_custom_results(out_filename=deflection_filename2, restype='Deflection')


    def test_beam_modes_01(self):
        """CBAR/CBEAM - PARAM,POST,-1"""
        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        op2_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes_m1.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_beam_modes_02(self):
        """CBAR/CBEAM - PARAM,POST,-2"""
        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        op2_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes_m2.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_beam_modes_03(self):
        dirname = os.path.join(MODEL_PATH, 'beam_modes')
        bdf_filename = os.path.join(dirname, 'beam_modes.dat')
        op2_filename = os.path.join(dirname, 'beam_modes_m1.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        #test.load_nastran_results(op2_filename)

        test.load_nastran_geometry(bdf_filename)
        #test.load_nastran_results(op2_filename)

        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_beam_modes_04(self):
        dirname = os.path.join(MODEL_PATH, 'beam_modes')
        bdf_filename = os.path.join(dirname, 'beam_modes.dat')
        op2_filename = os.path.join(dirname, 'beam_modes_m2.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

        test.load_nastran_geometry(bdf_filename)


    #@unittest.expectedFailure
    #def test_contact(self):
        #"""this test fails because of a misparsed card"""
        #bdf_filename = os.path.join(MODEL_PATH, 'contact', 'contact.bdf')
        #op2_filename = os.path.join(MODEL_PATH, 'contact', 'contact.op2')

        #test = NastranGUI()
        #test.load_nastran_geometry(bdf_filename)
        #test.load_nastran_results(op2_filename)

    def test_fsi(self):
        """tests -1 coordinate systems (flag for a fluid contact face)"""
        bdf_filename = os.path.join(MODEL_PATH, 'fsi', 'fsi.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'fsi', 'fsi.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_thermal_01(self):
        dirname = os.path.join(MODEL_PATH, 'thermal')
        bdf_filename = os.path.join(dirname, 'thermal_test_153.bdf')
        op2_filename = os.path.join(dirname, 'thermal_test_153.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_bwb_gui(self):
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        test = NastranGUI()
        #test.log = get_logger2()
        test.load_nastran_geometry(bdf_filename)
        test.group_actions.create_groups_by_property_id()
        test.group_actions.create_groups_by_visible_result(nlimit=50)
        test.toggle_conms()

    def test_femap_rougv1_01(self):
        """tests the exhaust manifold and it's funny eigenvectors"""
        dirname = os.path.join(MODEL_PATH, 'femap_exhaust')
        #bdf_filename = os.path.join(dirname, 'modal_example.bdf')
        op2_filename = os.path.join(dirname, 'modal_example.op2')

        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_aero(self):
        """tests the bah_plane"""
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)
        out_datai = deepcopy(test.geometry_properties)
        test.on_update_geometry_properties_override_dialog(out_datai)

        out_data = {
            'clicked_ok' : True,
            'Global XYZ' : out_datai['Global XYZ'],
            'conm2' : out_datai['conm2'],
            'bar_z' : out_datai['bar_z'],
            'caero' : out_datai['caero'],
        }

        #print(test.geometry_properties)
        coord = out_data['Global XYZ']
        coord.is_visible = False
        str(coord)
        #print('coord = %r' % coord)

        conm2 = out_data['conm2']
        conm2.point_size = 10

        barz = out_data['bar_z']
        barz.bar_scale = 0.5
        barz.is_visible = True
        #print(barz)

        caero = test.geometry_properties['caero']
        str(caero)
        caero.color = (255, 0, 0)
        caero.line_width = 10
        caero.opacity = 0.8
        caero.is_visible = False
        #print(caero)
        #print(out_data)
        test.on_update_geometry_properties(out_data, name='caero',
                                           write_log=True)

    def test_gui_elements_01(self):
        """tests forces/pressure in SOL 101"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

        idisp = None
        iforce_xyz = None
        for key, case_data in test.result_cases.items():
            case, data = case_data
            #print(key, case)
            if idisp is None and case.uname == 'Displacement':
                idisp = key
            elif idisp is not None and iforce_xyz is None and case.uname == 'LoadVectors':
                iforce_xyz = key
                break
            elif key > 70:
                break

        ifringe = len(test.result_cases) - 1 # Strain Energy Density
        test.on_fringe(icase=ifringe, stop_on_failure=True)
        with self.assertRaises(ValueError):
            test.on_vector(icase=ifringe, stop_on_failure=True)
        with self.assertRaises(ValueError):
            test.on_disp(icase=ifringe, stop_on_failure=True) # disp

        test.on_fringe(icase=iforce_xyz, stop_on_failure=True)
        test.on_vector(icase=iforce_xyz, stop_on_failure=True)
        test.on_disp(icase=idisp, stop_on_failure=True) # disp
        test.on_clear_results()

        test.on_fringe(icase=iforce_xyz, stop_on_failure=True)
        test.on_vector(icase=iforce_xyz, stop_on_failure=True) # force_xyz
        test.on_disp(icase=idisp, stop_on_failure=True) # disp
        test.on_fringe(icase=37, update_legend_window=True, show_msg=True)  # normal

    def test_gui_elements_01b(self):
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        #model = read_bdf(bdf_filename)
        model = BDF()
        model.disable_cards(['CHEXA', 'CTETRA', 'CPENTA',
                             'CROD', 'PLOTEL', 'CBAR', 'CBEAM', 'CTRIA3', 'CQUAD4', 'CQUADR', 'CTRIAR',
                             'CQUAD8', 'CTRIA6', 'CSHEAR', 'CTUBE',
                             'CONM2', 'CVISC', #'CONROD',
                             'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',

                             'PLOAD1', 'PLOAD2', 'PLOAD4',
                             ])
        model.read_bdf(bdf_filename)
        #print(model.elements)
        test2 = NastranGUI()
        test2.load_nastran_geometry(model)
        test2.load_nastran_results(op2_filename)

    def test_gui_elements_02(self):
        """tests a large number of elements and results in SOL 101"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

        #test = NastranGUI()
        test.settings.nastran_create_coords = False
        test.settings.nastran_is_bar_axes = False
        test.settings.nastran_is_3d_bars = False
        test.settings.nastran_is_3d_bars_update = False
        test.settings.nastran_is_element_quality = False
        test.settings.nastran_is_properties = False
        test.load_nastran_geometry(op2_filename)

    def test_gui_elements_03(self):
        """tests a large number of elements and results in SOL 103-modes"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)
        #test.create_groups_by_property_id()
        test.create_groups_by_visible_result()

    def test_gui_elements_04(self):
        """tests a large number of elements and results in SOL 108-freq"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)
        icase = 33
        name = 'Normal'
        subcase_id = -1
        test.set_normal_result(icase, name, subcase_id)

        test.setup_fake_text_actors()
        icase = 0
        icase2 = icase + 1
        while icase2 < len(test.result_cases):
            #test.on_cycle_results(case=icase2, show_msg=True)
            unused_result_name = 'dummy'
            test._set_case(unused_result_name, icase2, explicit=False, cycle=False,
                           skip_click_check=False, min_value=None, max_value=None,
                           is_legend_shown=None, show_msg=True)
            icase2 += 1

    def test_gui_elements_05(self):
        """tests a large number of elements and results in SOL 108-freq"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_elements_06(self):
        """tests a large number of elements and results in SOL 106-loadstep"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_elements_07(self):
        """tests a large number of elements and results in SOL 107-complex modes"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_elements_08(self):
        """tests a large number of elements and results in SOL 109-linear time"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_pload_01(self):
        """tests a PLOAD4/CTETRA"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'ctetra.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'ctetra.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_pload_02(self):
        """tests a PLOAD4/CHEXA"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'chexa.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'chexa.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_pload_03(self):
        """tests a PLOAD4/CPENTA"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'cpenta.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'cpenta.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_pload_04(self):
        """tests a PLOAD4/CQUAD4"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'cquad4.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'cquad4.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_pload_05(self):
        """tests a PLOAD4/CTRIA3"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'ctria3.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'ctria3.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    #def test_gui_pload_06(self):
        #"""tests a PLOAD1/CBAR"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'pload1.bdf')
        #op2_filename = os.path.join(MODEL_PATH, 'unit', 'pload4', 'pload1.op2')
        #test = NastranGUI()
        #test.load_nastran_geometry(op2_filename)
        #test.load_nastran_results(op2_filename)

    #def test_gui_bar_rod(self):
        #"""tests a PBARL/ROD"""
        #bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_rod.bdf')
        #test = NastranGUI()
        #test.load_nastran_geometry(bdf_filename)

    #def test_gui_bar_tube2(self):
    def test_gui_bar_tube(self):
        """tests a PBARL/TUBE"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_tube.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_chan(self):
        """tests a PBARL/CHAN"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_chan.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.on_pan_left(None)
        test.on_pan_right(None)
        test.on_pan_up(None)
        test.on_pan_down(None)
        test.on_increase_magnification()
        test.on_decrease_magnification()
        test.zoom(1.2)
        test.on_rotate_clockwise()
        test.on_rotate_cclockwise()
        test.rotate(15.0)
        test.set_focal_point([0., 1., 2.])
        test.export_case_data(icases=[0, 1])

        test.update_camera('+x')
        test.update_camera('-x')
        test.update_camera('+y')
        test.update_camera('-y')
        test.update_camera('+z')
        test.update_camera('-z')
        test._update_camera()
        camera_data = test.get_camera_data()
        test.on_set_camera_data(camera_data, show_log=True)

        csv_filename = os.path.join(MODEL_PATH, 'custom_geom.csv')
        test.on_load_user_geom(csv_filename=csv_filename, name=None, color=None)

        stl_filename = os.path.join(STL_PATH, 'sphere.stl')
        test.on_load_user_geom(csv_filename=stl_filename, name=None, color=None)
        test.clear_labels()
        test.reset_labels()

        with open('xyz1.csv', 'w') as xyz_file:
            xyz_file.write('1., 2., 3.\n')
            xyz_file.write('4., 5., 6.\n')
        csv_filename = 'xyz1.csv' # os.path.join(MODEL_PATH, 'xyz1.csv')
        test.on_load_csv_points(csv_filename=csv_filename, name=None, color=None)
        os.remove(csv_filename)

        with open('xyz2.csv', 'w') as xyz_file:
            xyz_file.write('10., 20., 30.')
        csv_filename = 'xyz2.csv' # os.path.join(MODEL_PATH, 'xyz2.csv')
        test.on_load_csv_points(csv_filename=csv_filename, name=None, color=None)
        os.remove(csv_filename)

        #test.on_wireframe()
        #test.on_surface()

        os.remove('0_NodeID.csv')
        os.remove('1_ElementID.csv')

        with open('rotate.py', 'w') as pyfile:
            pyfile.write('self.rotate(20.)\n')
        test.on_run_script('rotate.py')
        os.remove('rotate.py')

    def test_gui_screenshot(self):
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_chan.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

        magnify = None

        render_large = vtk.vtkRenderLargeImage()
        test.run_vtk = True
        #test.create_corner_axis()

        # faking coordinate system
        axes_actor = vtk.vtkAxesActor()
        test.corner_axis = vtk.vtkOrientationMarkerWidget()
        test.corner_axis.SetOrientationMarker(axes_actor)

        #test.on_take_screenshot(fname='chan.png', magnify=None, show_msg=True)
        out = test.tool_actions._screenshot_setup(magnify, render_large)
        line_widths0, point_sizes0, coord_scale0, coord_text_scale0, linewidth0, fake_axes_actor, magnify = out
        test.tool_actions._screenshot_teardown(
            line_widths0, point_sizes0, coord_scale0, coord_text_scale0, linewidth0, axes_actor)

    def test_gui_bar_chan1(self):
        """tests a PBARL/CHAN1"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_chan1.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
    #def test_gui_bar_chan2(self):

    def test_gui_bar_bar(self):
        """tests a PBARL/BAR"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_bar.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_box(self):
        """tests a PBARL/BOX"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_box.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_z(self):
        """tests a PBARL/Z"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_z.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_t(self):
        """tests a PBARL/T"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_t.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_t1(self):
        """tests a PBARL/T1"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_t1.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        str(test.geometry_properties)

        T1z = deepcopy(test.geometry_properties['T1_z'])
        T1z.line_width = 4
        T1z.color = (255, 0, 0)
        T1z.opacity = 0.6
        T1z.bar_scale = 0.20
        test.edit_geometry_properties_obj.on_update_geometry_properties({'T1_z' : T1z}, name=None, write_log=True)
        test.edit_geometry_properties_obj.on_update_geometry_properties({'T1_z' : T1z}, name='T1_z', write_log=True)

    def test_gui_bar_t2(self):
        """tests a PBARL/T2"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_t2.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_hexa(self):
        """tests a PBARL/HEXA"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_hexa.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_hat(self):
        """tests a PBARL/HAT"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_hat.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_i(self):
        """tests a PBARL/I"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_i.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_i1(self):
        """tests a PBARL/I1"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_i1.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_bar_h(self):
        """tests a PBARL/H"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_h.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_beam_l(self):
        """tests a PBEAML/L"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbeaml_l.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        total_length = test.model.get_length_breakdown()[100]
        assert np.allclose(total_length, 100.)

    def test_gui_thermal_01(self):
        """tests thermal"""
        #bdf_filename = os.path.join(MODEL_PATH, 'thermal', 'thermal_test_153.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'thermal', 'thermal_test_153.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_thermal_02(self):
        """tests thermal"""
        bdf_filename = os.path.join(MODEL_PATH, 'thermal', 'hd15901.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'thermal', 'hd15901.op2')
        test = NastranGUI()
        with self.assertRaises(DuplicateIDsError):
            test.load_nastran_geometry(op2_filename)
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_thermal_03(self):
        """tests thermal"""
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_superelement_1(self):
        """tests flyswatter"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'flyswatter', 'flyswatter_renumber.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        #test.load_nastran_results(op2_filename)

    def test_gui_superelement_2(self):
        """tests superelement mirror/shift/renumber"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'see103q4.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        os.remove('spike.bdf')
        os.remove('super_12.bdf')
        os.remove('super_13.bdf')
        os.remove('super_15.bdf')

    def test_gui_dvprel(self):
        """tests dvprel"""
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'dofm12.bdf')
        #op2_filename = os.path.join(MODEL_PATH, 'other', 'dofm12.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        #test.load_nastran_results(op2_filename)

    def test_gui_optimization_mcpads4(self):
        """tests mcpads4.bdf, which tests *.des and convergence"""
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'mcpads4.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'mcpads4.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_patran(self):
        """tests patran format"""
        bdf_filename = os.path.join(MODEL_PATH, 'patran_fmt', '0012_20.bdf')
        nod_filename = os.path.join(MODEL_PATH, 'patran_fmt', 'normals.nod')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(nod_filename)

    def test_gui_patran2(self):
        """tests patran format"""
        bdf_filename = os.path.join(MODEL_PATH, 'patran_fmt', '0012_20.bdf')
        nod_filename = os.path.join(MODEL_PATH, 'patran_fmt', 'normals.nod')
        test = NastranGUI()
        test.on_load_geometry(bdf_filename, geometry_format='nastran', raise_error=True)
        test.on_load_custom_results(out_filename=nod_filename, restype='Patran_nod')

    def test_gui_axi(self):
        """tests axisymmetric elements"""
        bdf_filename = os.path.join(MODEL_PATH, 'axisymmetric', 'model.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)

    def test_gui_ogs(self):
        """
        tests ogs.op2:
         - GridPointSurfaceStressesArray
        """
        bdf_filename = os.path.join(MODEL_PATH, 'ogs', 'ogs.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'ogs', 'ogs.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_bdf_op2_other_23(self):
        """
        tests ehbus69.op2:
         - RealBush1DStressArray
         - GridPointStressesVolumeDirectArray
         - GridPointStressesVolumePrincipalArray
         - GridPointStressesSurfaceDiscontinutiesArray
        """
        # TODO: support these results...
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'ehbus69.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ehbus69.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)


    def test_gui_zona_model_1(self):
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'f16_ma41.bdf')
        test = NastranGUI()
        test.log = SimpleLogger(level='error', encoding='utf-8', log_func=None)
        test.load_nastran_geometry(bdf_filename)

    def test_gui_zona_model_2(self):
        bdf_file = get_zona_model()
        test = NastranGUI()
        test.log = SimpleLogger(level='error', encoding='utf-8', log_func=None)
        test.load_nastran_geometry(bdf_file)

#def test_bottle():  # pragma: no cover
    #"""
    #Tests Nastran GUI loading
    #"""
    #test = NastranGUI()
    #test.load_nastran_geometry('bottle_shell_w_holes_pmc.bdf', '')
    #test.load_nastran_results('bottle_shell_w_holes_pmc.op2', '')

    #keys = test.result_cases.keys()
    #assert (1, 'Stress1', 1, 'centroid', '%.3f') in keys, keys

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
