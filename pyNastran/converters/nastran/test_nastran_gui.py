"""tests the NastranIO class"""
# encoding: utf8
import os
import sys
from pathlib import Path
from copy import deepcopy
import getpass
from collections import defaultdict
import unittest

import numpy as np
from pyNastran.converters.neu.neu import read_neu

try:
    import matplotlib
    matplotlib.use('Agg')
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pyparsing is missing
    IS_MATPLOTLIB = False

import vtkmodules
from vtk import vtkRenderLargeImage, vtkAxesActor, vtkOrientationMarkerWidget

from cpylog import SimpleLogger

import pyNastran
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.bdf.cards.test.test_aero import get_zona_model
from pyNastran.bdf.errors import DuplicateIDsError

from pyNastran.gui import (
    USE_NEW_SIDEBAR_OBJS_ as USE_NEW_SIDEBAR_OBJS,
    USE_NEW_TERMS_ as USE_NEW_TERMS)

USE_OLD_TERMS = not USE_NEW_TERMS
from pyNastran.gui.testing_methods import FakeGUIMethods

from pyNastran.converters.nastran.gui.nastran_io import NastranIO
from pyNastran.converters.nastran.nastran_to_vtk import nastran_to_vtk, save_nastran_results
from pyNastran.converters.nastran.gui.stress import get_composite_sort

from pyNastran.gui.qt_files.gui_attributes import IS_CUTTING_PLANE

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult, GridPointForceResult
from pyNastran.gui.gui_objects.displacements import (
    #DisplacementResults,
    ForceTableResults, ElementalTableResults)

from pyNastran.gui.gui_objects.force_results import ForceResults2
from pyNastran.gui.gui_objects.displacement_results import DisplacementResults2

from pyNastran.converters.nastran.gui.result_objects.simple_table_results import SimpleTableResults
from pyNastran.converters.nastran.gui.result_objects.layered_table_results import LayeredTableResults
from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2

RED_FLOAT = (1., 0., 0.)


class NastranGUI(NastranIO, FakeGUIMethods):
    def __init__(self, inputs=None):
        FakeGUIMethods.__init__(self, inputs=inputs)
        NastranIO.__init__(self)
        self.build_fmts(['nastran'], stop_on_failure=True)
        self.stop_on_failure = True

    def load_nastran_geometry(self, bdf_filename: PathLike | BDF,
                         name: str='main',
                         plot: bool=True,
                         stop_on_failure: bool=False):
        super().load_nastran_geometry(
            bdf_filename, name=name,
            plot=plot, stop_on_failure=stop_on_failure)
        self.validate_result_object_methods()

    def load_nastran_results(self, op2_filename: PathLike | OP2):
        self.stop_on_failure = True
        super().load_nastran_results(op2_filename)
        self.validate_result_object_methods()
        vtk_ugrid = self.grid
        save_nastran_results(self.gui, vtk_ugrid)

    def write_result_cases(self):  # pramga: no cover
        case_id0 = 0
        for case_id, (result_case, flag) in self.result_cases.items():

            if result_case.class_name == 'SimpleTableResults':
                aflag, bflag = flag
                bflag2 = list(bflag)
                imethod = bflag2[1]
                #flag2 = [aflag, bflag2]
                bflag2.append(result_case.methods[imethod])
                print(f"{case_id}, {str(bflag2)}, '{result_case.uname}'", result_case.class_name)
            elif hasattr(result_case, 'headers'):
                if result_case.headers:
                    print(case_id, flag, result_case.headers, result_case.class_name)
                else:
                    print(case_id, flag, result_case.__class__.__name__, result_case.methods)
            elif hasattr(result_case, 'header'):
                print(case_id, flag, result_case.header, result_case.class_name)
            else:
                raise RuntimeError(result_case)
            assert case_id == case_id0, (case_id, case_id0)
            case_id0 += 1

    def validate_result_object_methods(self):
        scale = 10.
        checks = defaultdict(bool)
            #'GuiResult': True,
            #'NormalResult': True,
            #'GridPointForceResult': True,
            #'CompositeStrainStressResults2': False,
            #'PlateStrainStressResults2': False,
            #'CompositeStrainStressResults2': False,
            #'SolidStrainStressResults2': False,
            #'PlateStrainStressResults2': False,
            #'PlateStrainStressResults2': False,
        #}
        for icase, (res, (itime, res_name)) in self.result_cases.items():
            # make one check per result type.
            # We check each term (e.g., for SolidStress, check oxx, oyy, von Mises, ...)
            is_complex = res.is_complex
            legend_title = res.get_legend_title(itime, res_name)
            key = (res.class_name, is_complex, legend_title)

            is_analyzed = checks[key]
            if is_analyzed:  # don't analyze the same object type twice
                continue

            res.get_data_format(itime, res_name)
            is_min = 'Min' in legend_title
            is_max = 'Max' in legend_title
            if isinstance(res, (GuiResult, SimpleTableResults, LayeredTableResults)):
                res.get_fringe_result(itime, res_name)
                res.get_fringe_vector_result(itime, res_name)
                res.get_annotation(itime, res_name)
                res.get_default_min_max(itime, res_name)
                res.get_imin_imax(itime, res_name)
            elif isinstance(res, NormalResult):
                res.get_annotation(itime, res_name)
            elif isinstance(res, GridPointForceResult):
                pass
            elif isinstance(res, ForceTableResults):
                # DisplacementResults, , ElementalTableResults
                res.get_annotation(itime, res_name)
                res.get_arrow_scale(itime, res_name)
                #res.get_case_flag(itime, res_name)

                #res.get_data_type(itime, res_name)
                res.get_fringe_result(itime, res_name)
                res.get_fringe_vector_result(itime, res_name)
                res.get_vector_result(itime, res_name)
                deflects = res.deflects(itime, res_name)
                if is_complex and deflects:
                    res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                res.get_default_arrow_scale(itime, res_name)
                res.get_default_min_max(itime, res_name)
                res.get_imin_imax(itime, res_name)

            elif isinstance(res, (ForceResults2, DisplacementResults2)):
                force_disp_flags = [
                    # transform, method_keys, nodal_combine
                    ('Global', None, ''),
                    ('Global', [0], ''),
                    ('Global', [0, 1], ''),
                    ('Global', [1, 2], ''),
                    ('Global', [0, 1, 2], ''),
                ]
                res.get_annotation(itime, res_name)
                res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)

                for transform, method_keys, nodal_combine in force_disp_flags:
                    res.set_sidebar_args(itime, res_name,
                                         #min_max_method=min_max_method,
                                         transform=transform,
                                         #methods_keys=method_keys,
                                         nodal_combine='')
                    #res.get_data_type(itime, res_name)
                    res.get_fringe_result(itime, res_name)
                    res.get_fringe_vector_result(itime, res_name)
                    res.get_vector_result(itime, res_name)
                    deflects = res.deflects(itime, res_name)
                    if is_complex and deflects:
                        res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    res.get_default_arrow_scale(itime, res_name)
                    res.get_default_min_max(itime, res_name)
                    res.get_imin_imax(itime, res_name)

            elif isinstance(res, CompositeStrainStressResults2):
                composite_flags = [
                    # min_max_method, transform, method_keys, nodal_combine
                    ('Absolute Max', 'Material', None, ''),
                    ('Mean', 'Material', None, ''),
                    ('Std. Dev.', 'Material', None, ''),
                    ('Difference', 'Material', None, ''),
                    ('Max', 'Material', [1, 2], ''),
                    ('Min', 'Material', [1, 2], ''),
                ]
                res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)
                #res.get_data_type(itime, res_name)
                for min_max_method, transform, method_keys, nodal_combine in composite_flags:
                    if (is_min, min_max_method) == (True, 'Max') or (is_max, min_max_method) == (True, 'Min'):
                        continue
                    res.set_sidebar_args(itime, res_name,
                                         min_max_method=min_max_method,
                                         transform=transform,
                                         methods_keys=method_keys,
                                         nodal_combine='')
                    res.get_annotation(itime, res_name)
                    res.get_fringe_result(itime, res_name)
                    res.get_fringe_vector_result(itime, res_name)
                    #res.get_vector_result(itime, res_name)
                    if is_complex:
                        res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    res.get_default_arrow_scale(itime, res_name)
                    res.get_default_min_max(itime, res_name)
                    res.get_imin_imax(itime, res_name)

            elif isinstance(res, PlateStrainStressResults2):
                plate_flags = [
                    # min_max_method, transform, method_keys, nodal_combine
                    ('Absolute Max', 'Material', None, 'Absolute Max'),
                    ('Mean', 'Material', None, 'Mean'),
                    ('Min', 'Material', None, 'Min'),
                    ('Max', 'Material', None, 'Max'),
                    ('Std. Dev.', 'Material', None, 'Std. Dev.'),
                    ('Difference', 'Material', None, 'Difference'),
                    ('Max', 'Material', [0], 'Max'),
                    ('Max', 'Material', [1], 'Max'),
                    ('Max', 'Material', [2], 'Max'),
                ]
                res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)
                #res.get_data_type(itime, res_name)
                for min_max_method, transform, method_keys, nodal_combine in plate_flags:
                    if (is_min, min_max_method) == (True, 'Max') or (is_max, min_max_method) == (True, 'Min'):
                        continue
                    res.set_sidebar_args(itime, res_name,
                                         min_max_method=min_max_method,
                                         transform=transform,
                                         methods_keys=method_keys,
                                         nodal_combine='')
                    res.get_annotation(itime, res_name)
                    res.get_fringe_result(itime, res_name)
                    res.get_fringe_vector_result(itime, res_name)
                    #res.get_vector_result(itime, res_name)
                    if is_complex:
                        res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    res.get_default_arrow_scale(itime, res_name)
                    res.get_default_min_max(itime, res_name)
                    res.get_imin_imax(itime, res_name)

            elif isinstance(res, SolidStrainStressResults2):
                solid_flags = [
                    # min_max_method, transform, method_keys, nodal_combine
                    ('Absolute Max', 'Material', None, 'Absolute Max'),
                    ('Mean', 'Material', None, 'Mean'),
                    ('Min', 'Material', None, 'Min'),
                    ('Max', 'Material', None, 'Max'),
                    ('Std. Dev.', 'Material', None, 'Std. Dev.'),
                    ('Difference', 'Material', None, 'Difference'),
                    ('Max', 'Material', [0], 'Max'),
                    ('Max', 'Material', [1], 'Max'),
                ]
                res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)
                #res.get_data_type(itime, res_name)
                for min_max_method, transform, method_keys, nodal_combine in solid_flags:
                    if (is_min, min_max_method) == (True, 'Max') or (is_max, min_max_method) == (True, 'Min'):
                        continue
                    res.set_sidebar_args(itime, res_name,
                                         min_max_method=min_max_method,
                                         transform=transform,
                                         methods_keys=method_keys,
                                         nodal_combine='')
                    res.get_annotation(itime, res_name)
                    res.get_fringe_result(itime, res_name)
                    res.get_fringe_vector_result(itime, res_name)
                    #res.get_vector_result(itime, res_name)
                    if is_complex:
                        res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    res.get_default_arrow_scale(itime, res_name)
                    res.get_default_min_max(itime, res_name)
                    res.get_imin_imax(itime, res_name)
            else:
                raise NotImplementedError(res)
            checks[key] = True
        return

    def get_results_by_type(self) -> dict[str, tuple]:
        results_by_type = defaultdict(list)
        for icase, (res, (itime, res_name)) in self.result_cases.items():
            #print(res.class_name)
            #class_name = res.__class__.__name__

            # make one check per result type.
            # We check each term (e.g., for SolidStress, check oxx, oyy, von Mises, ...)
            is_complex = res.is_complex
            legend_title = ''
            #key = (res.class_name, is_complex, legend_title)

            # legend_title = res.get_legend_title(itime, res_name)
            # res.get_data_format(itime, res_name)
            # is_min = 'Min' in legend_title
            # is_max = 'Max' in legend_title
            results_by_type[res.class_name].append((icase, (itime, res_name), res))

            if isinstance(res, (GuiResult, SimpleTableResults, LayeredTableResults)):
                res.get_fringe_result(itime, res_name)
                res.get_fringe_vector_result(itime, res_name)
                res.get_annotation(itime, res_name)
                res.get_default_min_max(itime, res_name)
                res.get_imin_imax(itime, res_name)
            elif isinstance(res, NormalResult):
                res.get_annotation(itime, res_name)
            elif isinstance(res, GridPointForceResult):
                pass
            elif isinstance(res, ForceTableResults):
                # DisplacementResults, , ElementalTableResults
                res.get_annotation(itime, res_name)
                res.get_arrow_scale(itime, res_name)
                #res.get_case_flag(itime, res_name)

                #res.get_data_type(itime, res_name)
                res.get_fringe_result(itime, res_name)
                res.get_fringe_vector_result(itime, res_name)
                res.get_vector_result(itime, res_name)
                deflects = res.deflects(itime, res_name)
                if is_complex and deflects:
                    scale = 1.0
                    res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                res.get_default_arrow_scale(itime, res_name)
                res.get_default_min_max(itime, res_name)
                res.get_imin_imax(itime, res_name)

            elif isinstance(res, (ForceResults2, DisplacementResults2)):
                continue
            elif isinstance(res, CompositeStrainStressResults2):
                composite_flags = [
                    # min_max_method, transform, method_keys, nodal_combine
                    ('Absolute Max', 'Material', None, ''),
                    ('Mean', 'Material', None, ''),
                    ('Std. Dev.', 'Material', None, ''),
                    ('Difference', 'Material', None, ''),
                    ('Max', 'Material', [1, 2], ''),
                    ('Min', 'Material', [1, 2], ''),
                ]
                res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)
                # res.get_data_type(itime, res_name)
                for min_max_method, transform, method_keys, nodal_combine in composite_flags:
                    if (is_min, min_max_method) == (True, 'Max') or (is_max, min_max_method) == (True, 'Min'):
                        continue
                    res.set_sidebar_args(itime, res_name,
                                         min_max_method=min_max_method,
                                         transform=transform,
                                         methods_keys=method_keys,
                                         nodal_combine='')
                    res.get_annotation(itime, res_name)
                    # res.get_fringe_result(itime, res_name)
                    # res.get_fringe_vector_result(itime, res_name)
                    # res.get_vector_result(itime, res_name)
                    # if is_complex:
                    #     res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    # res.get_default_arrow_scale(itime, res_name)
                    # res.get_default_min_max(itime, res_name)
                    # res.get_imin_imax(itime, res_name)

            elif isinstance(res, PlateStrainStressResults2):
                plate_flags = [
                    # min_max_method, transform, method_keys, nodal_combine
                    ('Absolute Max', 'Material', None, 'Absolute Max'),
                    ('Mean', 'Material', None, 'Mean'),
                    ('Min', 'Material', None, 'Min'),
                    ('Max', 'Material', None, 'Max'),
                    ('Std. Dev.', 'Material', None, 'Std. Dev.'),
                    ('Difference', 'Material', None, 'Difference'),
                    ('Max', 'Material', [0], 'Max'),
                    ('Max', 'Material', [1], 'Max'),
                    ('Max', 'Material', [2], 'Max'),
                ]
                # res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)
                # res.get_data_type(itime, res_name)
                for min_max_method, transform, method_keys, nodal_combine in plate_flags:
                    if (is_min, min_max_method) == (True, 'Max') or (is_max, min_max_method) == (True, 'Min'):
                        continue
                    res.set_sidebar_args(itime, res_name,
                                         min_max_method=min_max_method,
                                         transform=transform,
                                         methods_keys=method_keys,
                                         nodal_combine='')
                    an = res.get_annotation(itime, res_name)
                    print(an)
                    asdf
                    # res.get_fringe_result(itime, res_name)
                    # res.get_fringe_vector_result(itime, res_name)
                    # res.get_vector_result(itime, res_name)
                    # if is_complex:
                    #     res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    # res.get_default_arrow_scale(itime, res_name)
                    # res.get_default_min_max(itime, res_name)
                    # res.get_imin_imax(itime, res_name)

            elif isinstance(res, SolidStrainStressResults2):
                solid_flags = [
                    # min_max_method, transform, method_keys, nodal_combine
                    ('Absolute Max', 'Material', None, 'Absolute Max'),
                    ('Mean', 'Material', None, 'Mean'),
                    ('Min', 'Material', None, 'Min'),
                    ('Max', 'Material', None, 'Max'),
                    ('Std. Dev.', 'Material', None, 'Std. Dev.'),
                    ('Difference', 'Material', None, 'Difference'),
                    ('Max', 'Material', [0], 'Max'),
                    ('Max', 'Material', [1], 'Max'),
                ]
                # res.get_arrow_scale(itime, res_name)
                res.get_case_flag(itime, res_name)
                # res.get_data_type(itime, res_name)
                for min_max_method, transform, method_keys, nodal_combine in solid_flags:
                    # if (is_min, min_max_method) == (True, 'Max') or (is_max, min_max_method) == (True, 'Min'):
                    #     continue
                    res.set_sidebar_args(itime, res_name,
                                         min_max_method=min_max_method,
                                         transform=transform,
                                         methods_keys=method_keys,
                                         nodal_combine='')
                    an = res.get_annotation(itime, res_name)

                    # res.get_fringe_result(itime, res_name)
                    # res.get_fringe_vector_result(itime, res_name)
                    # res.get_vector_result(itime, res_name)
                    # if is_complex:
                    #     res.get_vector_result_by_scale_phase(itime, res_name, scale, phase=45.)
                    # res.get_default_arrow_scale(itime, res_name)
                    # res.get_default_min_max(itime, res_name)
                    # res.get_imin_imax(itime, res_name)
            else:
                raise NotImplementedError(res)
        return dict(results_by_type)


PKG_PATH = Path(pyNastran.__path__[0])
STL_PATH = PKG_PATH / 'converters' / 'stl'
MODEL_PATH = PKG_PATH / '..' / 'models'
FLUTTER_PATH = PKG_PATH / 'bdf' / 'cards' / 'aero' / 'examples' / 'flutter'


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

        color = RED_FLOAT
        opacity = 0.4
        test.settings.set_background_color(color, render=True)
        test.settings.set_background_color2(color, render=True)
        test.settings.set_highlight_color(color)
        test.settings.set_highlight_opacity(opacity)
        test.settings.set_corner_text_color(color, render=True)
        test.settings.set_corner_text_size(10)
        test.settings.set_magnify(magnify=4)
        #self.settings.s

    def test_nastran_f16_aero(self):
        dirname = FLUTTER_PATH / 'case5'
        neu_filename = dirname / 'f16-aero.neu'
        test = NastranGUI()
        test.load_nastran_geometry(neu_filename)
        test.load_nastran_results(neu_filename)

    def test_nastran_flut_anti(self):
        dirname = FLUTTER_PATH / 'case5'
        neu_filename = dirname / 'flut_anti.neu'
        test = NastranGUI()
        test.load_nastran_geometry(neu_filename)
        test.load_nastran_results(neu_filename)

    def test_nastran_cp_anti(self):
        dirname = FLUTTER_PATH / 'case5'
        neu_filename = dirname / 'cp2anti.neu'
        test = NastranGUI()
        test.load_nastran_geometry(neu_filename)
        test.load_nastran_results(neu_filename)

    def test_solid_shell_bar_obj(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        obj_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.obj')
        model = BDF()
        model.read_bdf(bdf_filename)
        model.save(obj_filename, unxref=True)

        test = NastranGUI()
        test.load_nastran_geometry(obj_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 56, len(test.result_cases)
        else:
            assert len(test.result_cases) == 56, len(test.result_cases)

    @unittest.skipIf(IS_MATPLOTLIB is False, 'No matplotlib')
    def test_solid_shell_bar_01(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 56, len(test.result_cases)
        else:
            assert len(test.result_cases) == 56, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 190, len(test.result_cases)
            #assert len(test.result_cases) == 206, len(test.result_cases)  # new? faked
        else:
            assert USE_OLD_TERMS
            # assert len(test.result_cases) == 278, len(test.result_cases)
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

        p1 = [0, 0, 0]
        p3 = [1, 0, 0]

        p2 = [0, 1, 0]
        zaxis = [0, 0, 1]

        test.shear_moment_torque_obj.setup_model_data(model_name)
        force_sum, moment_sum = test.shear_moment_torque_obj.plot_shear_moment_torque(
            icase_gpforce,
            p1, p2, p3, zaxis,
            method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
            nplanes=5,
            #plane_color=None, plane_opacity=0.5,
            csv_filename=None, show=False, stop_on_failure=True)
        assert np.allclose(np.abs(force_sum).max(), 0.000732421875), np.abs(force_sum).max()
        assert np.allclose(np.abs(moment_sum).max(), 0.000244140625), np.abs(moment_sum).max()

        p1 = np.array([0, 0, 0])  # origin
        p2 = np.array([1, 0, 0])  # xaxis
        p3 = np.array([1, 0, 0])  # end
        zaxis = np.array([0, 0, 1])
        test.shear_moment_torque_obj.plot_shear_moment_torque(
            icase_gpforce,
            p1, p2, p3, zaxis,
            method='Z-Axis Projection',
            cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
            nplanes=5,
            #plane_color=None, plane_opacity=0.5,
            csv_filename=None, show=False, stop_on_failure=True)

        if IS_CUTTING_PLANE:
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
        if IS_CUTTING_PLANE:
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
        if IS_CUTTING_PLANE:
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
        assert len(test.result_cases) == 55, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        assert len(test.models['main'].elements) > 0
        #test.write_result_cases()

        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 430, len(test.result_cases)
        else:
            assert USE_OLD_TERMS
            # assert len(test.result_cases) == 694, len(test.result_cases)

        #assert_result_cases(test, ncases=694)

        test.on_rcycle_results()
        test.on_update_legend(
            title='Title', min_value=0., max_value=1.,
            filter_value=1.,
            scale=0.0, phase=0.0,
            arrow_scale=1.,
            data_format='%.0f',
            is_low_to_high=True, is_discrete=True, is_horizontal=True,
            nlabels=None, labelsize=None, ncolors=None, colormap=None,
            is_shown=True, render=True)
        test.on_update_legend(
            title='Title', min_value=0., max_value=1.,
            filter_value=None,
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

        with self.assertRaises(RuntimeError):  # icase_to_apply=166 doesn't exist
            test.mark_elements(eids, stop_on_failure=True, show_command=True)
        #test.mark_elements_by_case(eids, stop_on_failure=True, show_command=True)

        test.icase = 2  # PropertyID
        test.mark_elements(eids, stop_on_failure=True, show_command=True)
        test.mark_elements_by_case(eids, stop_on_failure=True, show_command=True)

        #for key, obj in test.result_cases.items():
            #print(key)
            #print(obj)

        # fail mapping strain energy because we're on NodeID
        test.icase = 0  # NodeID
        test.icase_fringe = 0  # NodeID
        is_passed = test.map_element_centroid_to_node_fringe_result(update_limits=True, show_msg=True)

        obj, (itime, name) = test.result_cases[test.icase]
        str(obj)
        assert is_passed is False, f'map_element_centroid_to_node_fringe_result should fail for NodeID\n{obj}'

        # map strain energy
        keys = list(test.result_cases.keys())
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 430, len(test.result_cases)
            #assert len(test.result_cases) == 478, len(test.result_cases)  # new?; faked
        else:
            assert USE_OLD_TERMS
            # assert len(test.result_cases) == 694, len(test.result_cases)

        #assert_result_cases(test, ncases=694)
        icase = keys[-1]
        obj, (itime, name) = test.result_cases[icase]
        test.icase_fringe = icase
        str(obj)

        title = obj.get_legend_title(itime, name)
        assert title == 'Strain Energy Density', str(obj)
        is_passed = test.map_element_centroid_to_node_fringe_result(update_limits=True, show_msg=False)
        assert is_passed == True, 'map_element_centroid_to_node_fringe_result failed'

    def test_solid_shell_bar_02b(self):
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'mode_solid_shell_bar.bdf')

        test = NastranGUI()
        test.on_load_geometry(infile_name=bdf_filename, geometry_format='nastran', name='main',
                              plot=True, stop_on_failure=True)

    def test_aero_02(self):
        """checks 0012_flutter.op2"""
        bdf_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.op2')
        #log = SimpleLogger(level='error', encoding='utf-8')
        test = NastranGUI()
        test.on_load_geometry(infile_name=bdf_filename, geometry_format='nastran', name='main',
                              plot=True, stop_on_failure=True)
        test.on_load_results(op2_filename)
        test.validate_result_object_methods()

    def test_stack_composites(self):
        e1 = np.array([
            [18,  1],
            [18,  2],
            [18,  3],
            [18,  4],
            [19,  1],
            [19,  2],
            [19,  3],
            [19,  4],
            [20,  1],
            [20,  2],
            [20,  3],
            [20,  4],
            [20,  5],
            [21,  1],
            [21,  2],
            [21,  3],
            [21,  4],
            [21,  5],
        ])

        e2 = np.array([
            [16,  1],
            [16,  2],
            [16,  3],
            [16,  4],
            [17,  1],
            [17,  2],
            [17,  3],
            [17,  4],
            [17,  5],
        ])
        e1_e2 = np.vstack([e1, e2])
        out, isort = get_composite_sort(e1_e2)
        assert out.shape == e1_e2.shape

    def test_solid_shell_bar_03(self):
        bdf_filename = MODEL_PATH / 'sol_101_elements' / 'buckling_solid_shell_bar.bdf'
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'buckling_solid_shell_bar.op2'

        test = NastranGUI()
        test.stop_on_failure = True
        test.load_nastran_geometry(bdf_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            # we lost 3 cases for SPCD
            assert len(test.result_cases) == 57, len(test.result_cases)
        else:
            assert len(test.result_cases) == 57, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            # we lost 3 cases for SPCD
            assert len(test.result_cases) == 640, len(test.result_cases)
            #assert len(test.result_cases) == 684, len(test.result_cases)  # new terms?
        else:
            assert len(test.result_cases) == 640, len(test.result_cases)
            assert USE_OLD_TERMS

    def test_solid_bending(self):
        dirname = MODEL_PATH / 'solid_bending'
        bdf_filename = dirname / 'solid_bending.bdf'
        op2_filename = dirname / 'solid_bending.op2'
        #op2_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_ogs.op2')
        deflection_filename1 = dirname / 'solid_bending_multi_deflection_node.txt'
        deflection_filename2 = dirname / 'solid_bending_multi_deflection_node_short.txt'

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 10, len(test.result_cases)
        else:
            assert len(test.result_cases) == 10, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        #results_by_type = test.get_results_by_type()
        #for res_type, results in results_by_type.items():
            #for icase, (itime, res_name), obj in results:
                #an = obj.get_annotation(itime, res_name)
                #print(f'{icase}: {res_type} {repr(res_name)} {an!r}')

        nresults = get_nreal_nresults(
            test,
            ndisplacement=1,
            nspc_force=1, nmpc_force=0,
            nload_vectors=0,
            #neigenvectors=0, nspring_stress=0, nspring_strain=0, nspring_force=0,
            nsolid_stress=1, nsolid_strain=0,
            nabs_stress=1, nabs_strain=0,
            nstrain_energy=0, ngrid_point_forces=0)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert nresults == 24, nresults
            assert len(test.result_cases) == 34, len(test.result_cases)
        else:
            assert USE_OLD_TERMS
            # assert nresults == 24, nresults # 49-10; ???
            # assert len(test.result_cases) == 39, len(test.result_cases)

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
        csv_filename1 = dirname / 'solid_bending_multi_node.csv'
        csv_filename2 = dirname / 'solid_bending_multi_node_extra.txt'
        csv_filename3 = dirname / 'solid_bending_multi_node_short.txt'

        # missing/extra nodes
        test.on_load_custom_results(out_filename=csv_filename1, restype='node', stop_on_failure=True)
        test.on_load_custom_results(out_filename=csv_filename2, restype='node', stop_on_failure=True)
        test.on_load_custom_results(out_filename=csv_filename3, restype='node', stop_on_failure=True)

        # missing nodes
        test.on_load_custom_results(out_filename=deflection_filename2, restype='Deflection')

    def test_solid_bending_missing_eids(self):
        """
        same as the nominal version, but:
         - remove a solid element
        """
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.op2')
        model = read_bdf(bdf_filename)

        # make the problem a little harder
        del model.elements[1]

        test = NastranGUI()
        test.load_nastran_geometry(model)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 10, len(test.result_cases)
        else:
            assert len(test.result_cases) == 10, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 34, len(test.result_cases)
        else:
            assert USE_NEW_TERMS
            assert len(test.result_cases) == 39, len(test.result_cases)

    @unittest.skipIf(getpass.getuser() != 'sdoyle', 'local test')
    def test_solid_bending_missing_nodes(self):
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        op2_filename1 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.op2')
        op2_filename2 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_extra_nodes.op2')

        model = read_op2(op2_filename=op2_filename1, debug=False, combine=False)
        print(list(model.displacements.keys()))
        disp = model.displacements[(1, 1, 1, 0, 0, '', '')]
        # print(disp.object_attributes())

        # ['acoustic_flag', 'analysis_code', 'analysis_fmt', 'approach_code', 'class_name', 'data', 'data_code',
        #  'data_frame', 'data_names', 'dataframe', 'device_code', 'dt', 'format_code', 'gridtype_str', 'h5_file',
        #  'headers', 'is_built', 'is_cid', 'is_complex', 'is_msc', 'is_real', 'is_sort1', 'is_sort2',
        #  'isubcase', 'itime', 'itotal', 'label', 'load_as_h5', 'lsdvmn', 'lsdvmns', 'name', 'node_gridtype',
        #  'nonlinear_factor', 'ntimes', 'ntotal', 'num_wide', 'ogs', 'pval_step', 'random_code', 'result_name', 'size',
        #  'sort_bits', 'sort_code', 'sort_method', 'subtitle', 'subtitle_original', 'superelement_adaptivity_index',
        #  'tCode', 'table_code', 'table_name', 'table_name_str', 'thermal', 'thermal_bits', 'title', 'words']

        nids = disp.node_gridtype[:, 0]
        gridtype = disp.node_gridtype[:, 1]
        nnids = len(nids)
        nid = nids[-1] + 1
        nids2 = np.arange(nid, nid+nnids, dtype=nids.dtype)
        node_gridtype2 = disp.node_gridtype.copy()
        node_gridtype2[:, 0] = nids2

        datai = disp.data[0, :, :]
        data = np.vstack([datai, datai])
        disp.node_gridtype = np.vstack([disp.node_gridtype, node_gridtype2])
        print(disp.node_gridtype)
        disp.data = data.reshape(1, 2*nnids, 6)

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 10, len(test.result_cases)
        else:
            assert len(test.result_cases) == 10, len(test.result_cases)

        test.load_nastran_results(model)

    def test_beam_modes_01(self):
        """CBAR/CBEAM - PARAM,POST,-1"""
        bdf_filename = MODEL_PATH / 'beam_modes' / 'beam_modes.dat'
        op2_filename = MODEL_PATH / 'beam_modes' / 'beam_modes_m1.op2'

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        assert len(test.result_cases) == 7, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        nmodes = 10
        nresults = get_nreal_nresults(
            test,
            neigenvectors=nmodes,
            nbar_stress=nmodes,
            nbeam_stress=nmodes,  # beam stress is dropped
            nbar_force=nmodes,
            nbeam_force=nmodes)  # beam force is dropped
        #assert nresults == 231, nresults  # 238-7
        #test.write_result_cases()
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 238, len(test.result_cases)
        else:
            assert len(test.result_cases) == 238, len(test.result_cases)

    @unittest.skipIf(getpass.getuser() != 'sdoyle', 'local test')
    def test_beam_modes_01_missing_eids(self):
        """
        same as test_beam_modes_01 except:
         - missing CBAR eid=1
         - missing GRID nid=1
        """
        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        op2_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes_m1.op2')
        model = read_bdf(bdf_filename)
        del model.elements[1]
        del model.nodes[1]
        model._type_to_id_map['CBAR'].remove(1)

        test = NastranGUI()
        test.stop_on_failure = False
        test.load_nastran_geometry(model)
        assert len(test.result_cases) == 7, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        #test.write_result_cases()
        assert len(test.result_cases) == 238, len(test.result_cases)

    def test_beam_modes_02(self):
        """CBAR/CBEAM - PARAM,POST,-2"""
        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        op2_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes_m2.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

        nmodes = 10
        nresults = get_nreal_nresults(
            test,
            neigenvectors=nmodes,
            nbar_stress=nmodes,
            nbeam_stress=nmodes,  # beam stress is dropped
            nbar_force=nmodes,
            nbeam_force=nmodes)  # beam force is dropped
        #assert nresults == 231, nresults  # 238-7
        assert len(test.result_cases) == 238, len(test.result_cases)

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
        nmodes = 10
        nresults = get_nreal_nresults(
            test,
            neigenvectors=nmodes,
            nbar_stress=nmodes,
            nbeam_stress=nmodes,  # beam stress is dropped
            nbar_force=nmodes,
            nbeam_force=nmodes)  # beam force is dropped
        #assert nresults == 231, nresults  # 238-7
        assert len(test.result_cases) == 238, len(test.result_cases)

    def test_beam_modes_04(self):
        dirname = os.path.join(MODEL_PATH, 'beam_modes')
        bdf_filename = os.path.join(dirname, 'beam_modes.dat')
        op2_filename = os.path.join(dirname, 'beam_modes_m2.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        assert len(test.result_cases) == 7, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        nmodes = 10
        nresults = get_nreal_nresults(
            test,
            neigenvectors=nmodes,
            nbar_stress=nmodes,
            nbeam_stress=nmodes,  # beam stress is dropped
            nbar_force=nmodes,
            nbeam_force=nmodes)  # beam force is dropped
        #assert nresults == 231, nresults  # 238-7
        assert len(test.result_cases) == 238, len(test.result_cases)

        test.load_nastran_geometry(bdf_filename)
        assert len(test.result_cases) == 7, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        assert len(test.result_cases) == 238, len(test.result_cases)

        test.load_nastran_geometry(bdf_filename)
        assert len(test.result_cases) == 7, len(test.result_cases)

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
        assert len(test.result_cases) == 33, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        assert len(test.result_cases) == 53, len(test.result_cases)

    def test_thermal_01(self):
        """runs models/thermal/thermal_test_153"""
        dirname = os.path.join(MODEL_PATH, 'thermal')
        bdf_filename = os.path.join(dirname, 'thermal_test_153.bdf')
        op2_filename = os.path.join(dirname, 'thermal_test_153.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        assert len(test.result_cases) == 9, len(test.result_cases)
        test.load_nastran_results(op2_filename)
        assert len(test.result_cases) == 10, len(test.result_cases)

    def test_bwb_gui(self):
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.op2')
        test = NastranGUI()

        model = read_bdf(bdf_filename, debug=None)
        #CTRIA3      8043  901512    7571    7569    7572
        #CQUAD4      8044  901512    7569    7570    7573    7572
        if 8043 in model.elements:
            del model.elements[8043]
            del model.elements[8044]
        #model._type_to_id_map['CTRIA3'].remove(8043)
        #model._type_to_id_map['CQUAD4'].remove(8044)

        test.load_nastran_geometry(model)
        assert len(test.result_cases) == 95, len(test.result_cases)
        if os.path.exists(op2_filename) and 0:  # pragma: no cover
            nresults = get_nreal_nresults(
                test,
                nspc_force=1,
                nmpc_force=1,
                ndisplacement=1,
                #nspring_stress=0, nspring_strain=0,
                #nrod_stress=0, ctube_stress=0, nconrod_stress=0,
                #nrod_strain=0, ctube_strain=0, nconrod_strain=0,
                nbar_stress=1, nbar_strain=1,
                #nbeam_stress=0,
                #nplate_stress=0,
                #nshear_stress=0,
                ncomposite_layers=10, ncomposite_plate_stress=1, ncomposite_plate_strain=1,
                #nsolid_stress=0,
                #nbar_force=0,
                #nplate_force=0
            )
            assert nresults == 85, (len(test.result_cases), nresults) #  -95
            raise RuntimeError(nresults)

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
        assert len(test.result_cases) == 36, len(test.result_cases)
        test.load_nastran_results(op2_filename)

        assert len(test.result_cases) == 56, len(test.result_cases)

    def test_aero_op2(self):
        """tests the freedlm model (OP2 with aero)"""
        #bdf_filename = os.path.join(MODEL_PATH, 'aero', 'freedlm', 'freedlm.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'aero', 'freedlm', 'freedlm.op2')
        test = NastranGUI()
        #test.log.level = 'debug'
        #test.load_nastran_geometry(bdf_filename)
        test.load_nastran_geometry(op2_filename)
        #assert_result_cases(test, ncases=150)
        assert len(test.result_cases) == 150, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        #test.write_result_cases()

        #print(test.result_cases[154])
        #print(test.result_cases[160])
        #assert_result_cases(test, ncases=236)
        if USE_OLD_TERMS:
            assert len(test.result_cases) == 220, len(test.result_cases)  # old terms
        else:
            assert len(test.result_cases) == 224, len(test.result_cases)  # new terms
        #print(test.result_cases)

    def test_vba1(self):
        """vibroacoustics"""
        test = NastranGUI()

        bdf_filename = MODEL_PATH / 'nx' / 'test_vba' / 'test_vba.bdf'
        test.load_nastran_geometry(bdf_filename)

        bdf_filename = MODEL_PATH / 'nx' / 'test_vba' / 'ac108vatv5tc.bdf'
        test.load_nastran_geometry(bdf_filename)

        bdf_filename = MODEL_PATH / 'nx' / 'test_vba' / 'acssn108presvar.bdf'
        test.load_nastran_geometry(bdf_filename)

    def test_aero(self):
        """tests the bah_plane"""
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        assert_result_cases(test, ncases=7)
        test.load_nastran_results(op2_filename)
        assert_result_cases(test, ncases=7)

        out_datai = deepcopy(test.geometry_properties)
        test.on_update_geometry_properties_override_dialog(out_datai)

        out_data = {
            'clicked_ok': True,
            'Global XYZ': out_datai['Global XYZ'],
            'conm2': out_datai['conm2'],
            'bar_z': out_datai['bar_z'],
            'caero': out_datai['caero'],
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
        ngeometry = 62
        assert len(test.result_cases) == 62, len(test.result_cases)

        test.load_nastran_results(op2_filename)

        #nastran_settings = test.settings.nastran_settings
        #nastran_settings.displacement = False
        #nastran_settings.spc_force = False
        #nastran_settings.mpc_force = False
        #nastran_settings.grid_point_force = False
        #nastran_settings.strain_energy = False
        #nastran_settings.applied_load = False
        # 12-base good
        # ----------------------------------------------
        # 18-force good
        #nastran_settings.force = False
        # ----------------------------------------------
        # 10-abs
        # 4-rod
        # 13-bar
        # 8-plate** shold be 10, but for now, it's 8 (we added two parameters)
        # 9-comp
        # 10-solid
        # 1-spring
        # =55
        #nastran_settings.stress = True
        # ----------------------------------------------
        # 10-abs
        # 9-comp
        # 4-rod
        # 13-bar
        # 10-solid
        # 1-spring
        # =47 good
        #nastran_settings.strain = False
        # ----------------------------------------------
        nresults = get_nreal_nresults(
            test,
            nspc_force=1, nmpc_force=1, ndisplacement=1,
            nload_vectors=1,
            #neigenvectors=0,
            nspring_stress=1, nspring_strain=1, nspring_force=1,
            ncrod_stress=1,  # ctube_stress=0, nconrod_stress=0,
            ncrod_strain=1,  # ctube_strain=0, nconrod_strain=0,
            nbar_stress=1, nbar_strain=1, nbar_force=1,
            #nbeam_stress=0, nbeam_strain=0, nbeam_force=0,
            nplate_stress=1, nplate_strain=1, nplate_force=1,
            #nshear_stress=0, nshear_strain=0, nshear_force=0,
            ncomposite_layers=5, ncomposite_plate_stress=1, ncomposite_plate_strain=1,
            nsolid_stress=1, nsolid_strain=1,
            nabs_stress=1, nabs_strain=1,
            #----------------------
            nstrain_energy=1, ngrid_point_forces=1,
        )
        test.write_result_cases()
        #assert nresults == 139, nresults  # 202-139; alt is 196-63=133
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 203, len(test.result_cases)
        else:
            assert USE_OLD_TERMS
            assert len(test.result_cases) == -1, len(test.result_cases) # new terms; ???

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
            elif key > 69:
                break

        ifringe = len(test.result_cases) - 1  # Strain Energy Density
        test.on_fringe(icase=ifringe, stop_on_failure=True)
        with self.assertRaises(ValueError):
            test.on_vector(icase=ifringe, stop_on_failure=True)
        with self.assertRaises(ValueError):
            test.on_disp(icase=ifringe, stop_on_failure=True)  # disp

        test.on_fringe(icase=iforce_xyz, stop_on_failure=True)
        test.on_vector(icase=iforce_xyz, stop_on_failure=True)
        test.on_disp(icase=idisp, stop_on_failure=True)  # disp
        test.on_clear_results()

        test.on_fringe(icase=iforce_xyz, stop_on_failure=True)
        test.on_vector(icase=iforce_xyz, stop_on_failure=True)  # force_xyz
        test.on_disp(icase=idisp, stop_on_failure=True)  # disp
        test.on_fringe(icase=37, update_legend_window=True, show_msg=True)  # normal

        #op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        vtu_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.vtu')
        nastran_to_vtk(op2_filename, op2_filename, vtu_filename, log_level='error')

        assert os.path.exists(vtu_filename), vtu_filename

    def _test_gui_vtk(self):  # pragma: no cover
        dirname = r'pyNastran\pyNastran\converters\nastran\models'
        bdf_filename = os.path.join(dirname, 'demo.bdf')
        op2_filename = os.path.join(dirname, 'demo.op2')
        vtu_filename = os.path.join(dirname, 'demo.vtu')
        nastran_to_vtk(bdf_filename, op2_filename, vtu_filename, log_level='error')

    def test_gui_vtk1_elements_01(self):
        """tests forces/pressure in SOL 101 using an op2 model object"""
        op2_filename = MODEL_PATH / 'elements' / 'static_elements.op2'
        vtu_filename = MODEL_PATH / 'elements' / 'static_elements0.vtu'
        vtk_filename = MODEL_PATH / 'elements' / 'static_elements.vtk'
        model = read_op2(op2_filename, load_geometry=True, combine=False, debug=False)
        nastran_to_vtk(model, model, vtu_filename, log_level='error')
        nastran_to_vtk(model, model, vtk_filename, log_level='error')

    def test_gui_vtk2_elements_01(self):
        """tests forces/pressure in SOL 101 using a Path object"""
        op2_filename = MODEL_PATH / 'elements' / 'static_elements.op2'
        vtu_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements9.vtu')
        nastran_to_vtk(op2_filename, op2_filename, vtu_filename, compression_level=9)

        vtu_filename = MODEL_PATH / 'elements' / 'static_elements0.vtu'
        nastran_to_vtk(op2_filename, op2_filename, vtu_filename, compression_level=0)

    def test_bdf_op2_64_bit(self):
        """
        checks d173.bdf, which tests MSC Nastran 64-bit without the
        op2.is_interlaced flag
        """
        dirname = MODEL_PATH / 'msc' / '64_bit'
        bdf_filename = os.path.join(dirname, 'd173.bdf')
        op2_filename = os.path.join(dirname, 'd173.op2')
        vtk_filename = os.path.join(dirname, 'd173.vtu')
        nastran_to_vtk(bdf_filename, op2_filename, vtk_filename)

    def test_gui_elements_01_missing_eids(self):
        """
        same as test_gui_elements_01 except missing a single:
         -
        """
        bdf_filename = MODEL_PATH / 'elements' / 'static_elements.bdf'
        op2_filename = MODEL_PATH / 'elements' / 'static_elements.op2'
        model = read_bdf(bdf_filename)
        test = NastranGUI()
        test.load_nastran_geometry(model)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 62, len(test.result_cases)
        else:
            assert len(test.result_cases) == 62, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 203, len(test.result_cases)
        else:
            assert USE_OLD_TERMS
            assert len(test.result_cases) == -1, len(test.result_cases)  # new terms; ???

    def _test_gui_elements_01b(self):  # pragma: no cover
        bdf_filename = MODEL_PATH / 'elements' / 'static_elements.bdf'
        op2_filename = MODEL_PATH / 'elements' / 'static_elements.op2'
        #model = read_bdf(bdf_filename)
        model = BDF()
        model.disable_cards(['CHEXA', 'CTETRA', 'CPENTA',
                             'CROD', 'PLOTEL', 'CBAR', 'CBEAM', 'CTRIA3', 'CQUAD4', 'CQUADR', 'CTRIAR',
                             'CQUAD8', 'CTRIA6', 'CSHEAR', 'CTUBE',
                             'CONM2', 'CVISC',  # 'CONROD',
                             'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',

                             'PLOAD1', 'PLOAD2', 'PLOAD4',
                             ])
        model.read_bdf(bdf_filename)
        #print(model.elements)
        test2 = NastranGUI()
        test2.stop_on_failure = False
        test2.load_nastran_geometry(model)
        assert len(test2.result_cases) == 5, len(test2.result_cases)
        test2.load_nastran_results(op2_filename)
        #print(test2.result_cases[27])
        #print(test2.write_result_cases())

        nresults = get_nreal_nresults(
            test2,
            nspc_force=1, nmpc_force=1, nload_vectors=1,
            ndisplacement=1, neigenvectors=0,
            #nspring_stress=0, nspring_strain=0, nspring_force=0,
            #ncrod_stress=0, ctube_stress=0, nconrod_stress=0,
            ncrod_strain=0, ctube_strain=1, nconrod_strain=1,
            nbar_stress=0, nbar_strain=1, nbar_force=0,
            #nbeam_stress=1, nbeam_strain=1, nbeam_force=0,
            nplate_stress=1, nplate_strain=1, nplate_force=0,
            nshear_stress=1, nshear_strain=0, nshear_force=0,
            #ncomposite_layers=0, ncomposite_plate_stress=1, ncomposite_plate_strain=1,
            nsolid_stress=0, nsolid_strain=1,
            nabs_stress=1, nabs_strain=1,
            nstrain_energy=1, ngrid_point_forces=1)

        if USE_OLD_TERMS:
            assert len(test2.result_cases) == 105, len(test2.result_cases)
        else:
            assert len(test2.result_cases) == 125, len(test2.result_cases)  # wrong

    def test_gui_elements_02(self):
        """tests a large number of elements and results in SOL 101"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        test = NastranGUI()
        test.load_nastran_geometry(op2_filename)

        assert USE_NEW_SIDEBAR_OBJS
        assert len(test.result_cases) == 57, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        assert USE_NEW_SIDEBAR_OBJS
        assert len(test.result_cases) == 198, len(test.result_cases)

        #test = NastranGUI()
        test.settings.nastran_create_coords = False
        test.settings.nastran_is_bar_axes = False
        test.settings.nastran_is_3d_bars = False
        test.settings.nastran_is_3d_bars_update = False
        test.settings.nastran_is_element_quality = False
        test.settings.nastran_is_properties = False
        test.load_nastran_geometry(op2_filename)

    def test_gui_elements_missing_pcomp_subcase(self):
        op2_filename = MODEL_PATH / 'bugs' / 'flat_plate_pcomp_test' / 'flat_plate_tip_loads_mixed_2cases.op2'
        test = NastranGUI()
        #test.stop_on_failure = False
        test.load_nastran_geometry(op2_filename)
        assert USE_NEW_SIDEBAR_OBJS
        assert len(test.result_cases) == 58, len(test.result_cases)

        test.load_nastran_results(op2_filename)
        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 116, len(test.result_cases)  # old terms
        else:
            assert len(test.result_cases) == 120, len(test.result_cases)  # new terms

    def test_gui_elements_03(self):
        """tests a large number of elements and results in SOL 103-modes"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.bdf')
        op2_filename = MODEL_PATH / 'elements' / 'modes_elements.op2'
        test = NastranGUI()
        #test.stop_on_failure = False
        test.load_nastran_geometry(op2_filename)
        assert USE_NEW_SIDEBAR_OBJS
        assert len(test.result_cases) == 56, len(test.result_cases)

        test.load_nastran_results(op2_filename)

        #junk_filename = r'C:\NASA\m4\formats\git\pyNastran\models\bwb\old.txt'
        #junk_filename = r'C:\NASA\m4\formats\git\pyNastran\models\bwb\new.txt'
        #with open(junk_filename, 'w') as fobj:
            #for key, value in test.result_cases.items():
                #fobj.write(f'{key}: {value}\n')

        if USE_NEW_SIDEBAR_OBJS and USE_OLD_TERMS:
            assert len(test.result_cases) == 503, len(test.result_cases)  # old terms
        else:
            assert USE_OLD_TERMS

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
        test.stop_on_failure = False
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

        render_large = vtkRenderLargeImage()
        test.run_vtk = True
        #test.create_corner_axis()

        # faking coordinate system
        axes_actor = vtkAxesActor()
        test.corner_axis = vtkOrientationMarkerWidget()
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
        test.edit_geometry_properties_obj.on_update_geometry_properties({'T1_z': T1z}, name=None, write_log=True)
        test.edit_geometry_properties_obj.on_update_geometry_properties({'T1_z': T1z}, name='T1_z', write_log=True)

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

    def test_gui_pcomp_01(self):
        """tests composite von mises cquad4"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pcomp', 'pcomp_cquad4.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pcomp', 'pcomp_cquad4.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
        test.load_nastran_results(op2_filename)

    def test_gui_pcomp_02(self):
        """tests composite von mises ctria3"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'pcomp', 'pcomp_ctria3.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'pcomp', 'pcomp_ctria3.op2')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename)
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
        test.on_load_geometry(bdf_filename, geometry_format='nastran', stop_on_failure=True)
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


class TestZonaGui(unittest.TestCase):
    def test_gui_zona_model_1(self):
        bdf_filename = MODEL_PATH / 'aero' / 'zona' / 'f16_ma41.bdf'
        test = NastranGUI()
        test.log = SimpleLogger(level='error', encoding='utf-8')
        test.load_nastran_geometry(bdf_filename)

    def test_gui_zona_model_2(self):
        bdf_file = get_zona_model()
        test = NastranGUI()
        test.log = SimpleLogger(level='error', encoding='utf-8')
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


def assert_result_cases(test: NastranGUI, ncases: int,
                        debug: bool=False) -> None:  # pragma: no cover
    msg = ''
    for i, (obj, name) in test.result_cases.items():
        #if i > 328:
            #break
        msg += f'{i}: {str(obj)}'
    assert isinstance(msg, str)
    ncases_actual = len(test.result_cases)
    if ncases_actual != ncases:
        #print(msg)
        print(f'ncases_actual={ncases_actual}; ncases={ncases}')
        raise RuntimeError(msg)
    if debug:
        assert isinstance(msg, str)
        if hasattr(sys.stdin, 'reconfigure'):
            sys.stdin.reconfigure(encoding='utf-8')
        if hasattr(sys.stdout, 'reconfigure'):
            sys.stdout.reconfigure(encoding='utf-8')
        print(msg)


def get_nreal_nresults_from_model(
        model: BDF,
        ntimes: int=1,
        ndisplacement=0, neigenvectors=0, nspc=0, nmpc=0,
        nstress=0, nstrain=0, nforce=0):  # pragma: no cover
    #cc = model.case_control_deck
    nspring = model.card_count['CELAS1'] + model.card_count['CELAS2'] + model.card_count['CELAS3'] + model.card_count['CELAS4']
    nshear = model.card_count['CSHEAR']
    ncrod = model.card_count['CROD']
    ntube = model.card_count['CTUBE']
    nconrod = model.card_count['CONROD']
    #ntri = model.card_count['CTRIA3'] + model.card_count['CTRIA6'] + model.card_count['CTRIAR']
    #nquad = model.card_count['CQUAD4'] + model.card_count['CQUAD8'] + model.card_count['CQUADR']
    nbar = model.card_count['CBAR']
    nbeam = model.card_count['CBEAM']
    ncomposite_layers = 0
    ncomposite_plates = 0
    nplate = 0
    nsolid = 0
    for prop in model.properties.items():
        if prop.type in {'PCOMP', 'PCOMPG'}:
            ncomposite_layers = max(ncomposite_layers, prop.nlayers)
            ncomposite_plates = 1
        elif prop.type == 'PSHELL':
            nplate = 1
        elif prop.type == 'PSOLID':
            nsolid = 1

    test = False
    nresults = get_nreal_nresults(
        test,
        nspc_force=nspc, nmpc_force=nmpc, ndisplacement=ndisplacement,
        neigenvectors=neigenvectors,
        nspring_stress=nspring*nstress, nspring_strain=nspring*nstrain,
        ncrod_stress=ncrod*nstress, ctube_stress=ntube*nstress, nconrod_stress=nconrod*nstress,
        ncrod_strain=ncrod*nstrain, ctube_strain=ntube*nstrain, nconrod_strain=nconrod*nstrain,
        nbar_stress=nbar*nstress, nbar_strain=nbar*nstrain,
        nbeam_stress=nbeam*nstress, nbeam_strain=nbeam*nstrain,
        nplate_stress=nplate*nstress, nplate_strain=nplate*nstrain,
        nshear_stress=nstrain*nstress, nshear_strain=nstrain*nstrain,
        ncomposite_layers=ncomposite_layers,
        ncomposite_plate_stress=ncomposite_plates*nstress,
        ncomposite_plate_strain=ncomposite_plates*nstrain,
        nsolid_stress=nsolid*nstress, nsolid_strain=nsolid*nstrain,
        nbar_force=nbar*nforce, nbeam_force=nbeam*nforce,
        nplate_force=nplate*nforce, nshear_force=nshear*nforce)
    return nresults


def get_nreal_nresults(
        test,
        nspc_force=0,
        nmpc_force=0,
        nload_vectors=0,
        ndisplacement=0,
        neigenvectors=0,
        nspring_stress=0, nspring_strain=0, nspring_force=0,
        ncrod_stress=0, ctube_stress=0, nconrod_stress=0,
        ncrod_strain=0, ctube_strain=0, nconrod_strain=0,
        nbar_stress=0, nbar_strain=0, nbar_force=0,
        nbeam_stress=0, nbeam_strain=0, nbeam_force=0,
        nplate_stress=0, nplate_strain=0, nplate_force=0,
        nshear_stress=0, nshear_strain=0, nshear_force=0,
        ncomposite_layers=0, ncomposite_plate_stress=0, ncomposite_plate_strain=0,
        nsolid_stress=0, nsolid_strain=0,
        nabs_stress=0, nabs_strain=0,
        #-------------------------
        nstrain_energy=0, ngrid_point_forces=0) -> int:
    """the beginnings of trying to calculate the number of results vs. guessing"""
    # txyz, rxyz
    nastran_settings = test.settings.nastran_settings
    if not nastran_settings.displacement:
        ndisplacement = 0
    if not nastran_settings.eigenvector:
        neigenvectors = 0
    if not nastran_settings.applied_load:
        nload_vectors = 0
    if not nastran_settings.spc_force:
        nspc_force = 0
    if not nastran_settings.mpc_force:
        nmpc_force = 0
    if not nastran_settings.grid_point_force:
        ngrid_point_forces = 0
    if not nastran_settings.strain_energy:
        nstrain_energy = 0
    ntables = (
        # *2 is for translation/rotation
        nspc_force + nmpc_force + nload_vectors +
        ndisplacement + neigenvectors) * 2

    if not nastran_settings.force:
        nspring_force = 0
        nbar_force = 0
        nbeam_force = 0
        nplate_force = 0
        nshear_force = 0
    if not nastran_settings.strain:
        nspring_strain = 0
        ncrod_strain = 0
        ctube_strain = 0
        nconrod_strain = 0
        nbar_strain = 0
        nbeam_strain = 0
        nplate_strain = 0
        nshear_strain = 0
        ncomposite_plate_strain = 0
        nsolid_strain = 0
        nabs_strain = 0
    if not nastran_settings.stress:
        nspring_stress = 0
        nconrod_stress = 0
        ncrod_stress = 0
        ctube_stress = 0
        nbar_stress = 0
        nbeam_stress = 0
        nplate_stress = 0
        nshear_stress = 0
        ncomposite_plate_stress = 0
        nsolid_stress = 0
        nabs_stress = 0

    nplate_stress_strain = nplate_stress + nplate_strain
    ncomposite_stress_strain = ncomposite_plate_stress + ncomposite_plate_strain
    nstress_strain = (
        # oxx
        (nspring_stress + nspring_strain) +

        # [oxx, MS_axial, txy, MS_torsion]
        (ncrod_stress + ctube_stress + nconrod_stress +
         ncrod_strain + ctube_strain + nconrod_strain) * 4 +

        #[s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
        # s1b, s2b, s3b, s4b,        smaxb, sminb, MS_compression]
        # no MS_tension/compression
        (nbar_stress + nbar_strain) * 13 +
        (nbeam_stress + nbeam_strain) * 0 +  # sxc, sxd, sxe, sxf, smax, smin

        # [max_shear, avg_shear, margin]; no margin
        (nshear_stress + nshear_strain) * 2 +

        #[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]; # 8; old
        nplate_stress_strain * 8 * USE_NEW_SIDEBAR_OBJS * 1 * USE_OLD_TERMS +  # 1 layer

        #[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises, oabs, max_shear]  # 10; new
        nplate_stress_strain * 10 * USE_NEW_SIDEBAR_OBJS * 1 * USE_NEW_TERMS +  # 1 layer
        #added abs_princpal and von_mises/max_shear, so +2 results; disabld for now for testing

        # [o11, o22, t12, t1z, t2z, angle, major, minor, max_shea]
        ncomposite_stress_strain * 9 * USE_NEW_SIDEBAR_OBJS +

        # [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, von_mises]
        (nsolid_stress + nsolid_strain) * 10  # OLD
    )

    nforce = (
        #[Fspring]
        nspring_force +
        # order is different, but close enough
        # [bending_moment_a1, bending_moment_a2, shear1,
        #  bending_moment_b1, bending_moment_b2, shear2,
        #  axial, torque]
        nbar_force * 8 +
        nbeam_force * 0 +

        # [Fx, Fy, Fxy, Mx, My, Mxy, Vx, Vy]
        nplate_force * 8
    )
    nresults = (
        ntables +
        nstress_strain + nforce +
        #[Strain Energy, Percent, Strain Energy Density]
        nstrain_energy*3 +
        ngrid_point_forces
    )
    if nbar_force:
        # isBarOn flag
        nresults += 1

    if nbeam_stress:
        # is_stress_off
        nresults += 1
    if nbeam_strain:
        # is_strain_off
        nresults += 1
    if nbeam_force:
        # is_force_off
        nresults += 1

    nplate_stress_strain = nplate_stress + nplate_strain
    nabs = nabs_stress + nabs_strain
    if nabs > 0:
        # absolute max corner stress
        nresults += nabs * 10
    assert nresults > 0, nresults
    return nresults


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
