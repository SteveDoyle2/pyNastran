"""
defines GuiQtCommon

This file defines functions related to the result updating that are VTK specific

"""
# coding: utf-8
# pylint: disable=C0111
import sys
from copy import deepcopy
from collections import namedtuple

import numpy as np
from numpy import issubdtype
from numpy.linalg import norm  # type: ignore
import vtk

#import pyNastran
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.aero.utils import points_elements_from_quad_points

from pyNastran.gui.gui_objects.names_storage import NamesStorage
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.qt_files.gui_attributes import GuiAttributes
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk, VTK_VERSION
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points
#from pyNastran.gui import IS_DEV
IS_TESTING = 'test' in sys.argv[0]

WHITE = (1., 1., 1.)
BLUE = (0., 0., 1.)
RED = (1., 0., 0.)

FringeData = namedtuple(
    'FringeData',
    'icase, result_type, location, min_value, max_value, norm_value,'
    'data_format, scale, methods,'
    'subcase_id, subtitle, label,'
    'nlabels, labelsize, ncolors, colormap,'
    'imin, imax')

DispData = namedtuple(
    'DispData',
    'icase, result_type, location, min_value, max_value, norm_value,'
    'data_format, scale, phase, methods,'
    'subcase_id, subtitle, label,'
    'nlabels, labelsize, ncolors, colormap,'
    'xyz_nominal, vector_data,'
    'is_checked,'
    'imin, imax')


class GuiQtCommon(GuiAttributes):
    def __init__(self, **kwds):
        inputs = kwds['inputs']
        kwds['res_widget'] = None
        super(GuiQtCommon, self).__init__(**kwds)
        self.is_groups = inputs['is_groups']

        #self.groups = set()
        self._group_elements = {}
        self._group_coords = {}
        self._group_shown = {}
        self._names_storage = NamesStorage()

        self.vtk_version = VTK_VERSION
        #if not IS_TESTING or not pyNastran.is_pynastrangui_exe:  # pragma: no cover
            #print('vtk_version = %s' % (self.vtk_version))
        if self.vtk_version[0] < 7:  # TODO: should check for 7.1
            raise RuntimeError('VTK %s is no longer supported' % vtk.VTK_VERSION)

    def _cycle_results(self, icase=0):
        """used for testing"""
        self.icase = icase
        self.on_cycle_results(show_msg=True)
        i = 0
        while self.icase != 0:
            self.on_cycle_results(show_msg=False)
            if i > 10000:   # pragma: no cover
                raise RuntimeError('max cycle count; i=%s' % i)
            i += 1

    def on_rcycle_results(self, show_msg=True):
        """the reverse of on_cycle_results"""
        if len(self.case_keys) <= 1:
            return

        #ncases = len(self.case_keys)
        #icasei = self.case_keys.index(self.icase)
        #icasei2 = ncases + icasei - 1
        #icase = (self.case_keys + self.case_keys)[icasei2]
        #self.cycle_results(icase)
        #self.icase = icase

        icase = self.icase - 1
        while 1:  # TODO: speed this up
            if icase == -1:
                icase = self.ncases - 1
            try:
                self.cycle_results(icase, show_msg=show_msg)
                break
            except IndexError:
                icase -= 1
        self.update_icase()

    def on_cycle_results(self, show_msg=True):
        """the gui method for calling cycle_results"""
        if len(self.case_keys) <= 1:
            return

        #ncases = len(self.case_keys)
        #icasei = self.case_keys.index(self.icase)
        #icasei2 = icasei + 1
        #icase = (self.case_keys + self.case_keys)[icasei2]
        #self.cycle_results(icase)
        #self.icase = icase

        icase = self.icase + 1
        while 1:  # TODO: speed this up
            if icase == self.ncases:
                icase = 0
            try:
                self.cycle_results(icase, show_msg=show_msg)
                break
            except IndexError:
                icase += 1
        self.update_icase()

    def update_icase(self):
        """updates the Qt sidebar"""
        self.res_widget.update_icase(self.icase)

    def cycle_results(self, case=None, show_msg=True):
        """
        Selects the next case

        Parameters
        ----------
        case : int; default=None
            selects the icase
            None : defaults to self.icase+1

        """
        #print('-----------------')
        #print('real-cycle_results(case=%r)' % case)
        if case is None:
            #print('adding 1 to case (%s)' % (self.icase))
            case = self.icase + 1

        assert case is not False, case
        ncases = len(self.case_keys)
        if ncases <= 1:
            self.log_warning('cycle_results(result_name=%r); ncases=%i' % (
                case, ncases))
            if self.ncases == 0:
                self.scalar_bar_actor.SetVisibility(False)
            return
        case = self.cycle_results_explicit(case, explicit=False, show_msg=show_msg)
        assert case is not False, case
        if show_msg:
            self.log_command('cycle_results(case=%r)' % self.icase)

    def get_new_icase(self):
        if len(self.result_cases):
            return max(self.result_cases) + 1
        return 0

    def update_result_cases(self, cases):
        """acts like result_cases.update(cases)"""
        for key, case in cases.items():
            assert key not in self.result_cases, 'key=%r is already used' % key
            self.result_cases[key] = case

    def get_subtitle_label(self, subcase_id):
        try:
            subtitle, label = self.isubcase_name_map[subcase_id]
        except TypeError:
            subtitle = 'case=NA'
            label = 'label=NA'
        except KeyError:
            subtitle = 'case=NA'
            label = 'label=NA'
        return subtitle, label

    def on_clear_results(self, show_msg=True):
        """clears the model of all results"""
        (obj, (i, name)) = self.result_cases[self.icase]
        unused_location = obj.get_location(i, name)

        unused_name_str = self._names_storage.get_name_string(name)
        grid = self.grid

        if self._is_displaced:
            if self._xyz_nominal is None:
                self.log_error('the model is displaced, but self._xyz_nominal is None...')
                return
            self._is_displaced = False
            self._update_grid(self._xyz_nominal)
        if show_msg:
            self.log_command(f'on_clear_results(show_msg={show_msg})')

    #def clear_grid_fringe(grid):
        point_data = grid.GetPointData()
        npoint_arrays = point_data.GetNumberOfArrays()
        for i in range(npoint_arrays):
            name = point_data.GetArrayName(i)
            if name is not None:
                #print("pointArray ", i, ": ", name)
                point_data.RemoveArray(name)

        cell_data = grid.GetCellData()
        ncell_arrays = cell_data.GetNumberOfArrays()
        for i in range(ncell_arrays):
            name = cell_data.GetArrayName(i)
            if name is not None:
                #print("cellArray ", i, ": ", name)
                cell_data.RemoveArray(name)

        cell_data.SetActiveScalars(None)
        point_data.SetActiveScalars(None)
        cell_data.SetActiveVectors(None)
        point_data.SetActiveVectors(None)
        self.arrow_actor.SetVisibility(False)
        self.arrow_actor_centroid.SetVisibility(False)

        prop = self.geom_actor.GetProperty()
        prop.SetColor(WHITE)

        # the backface property could be null
        back_prop = vtk.vtkProperty()
        back_prop.SetColor(WHITE)
        self.geom_actor.SetBackfaceProperty(back_prop)
        self.geom_actor.Modified()

        self.hide_legend()
        self.scalar_bar.is_shown = False

        self.clear_legend()
        self.vtk_interactor.Render()

        self.res_widget.result_case_window.treeView.fringe.setChecked(False)
        self.res_widget.result_case_window.treeView.disp.setChecked(False)
        self.res_widget.result_case_window.treeView.vector.setChecked(False)
        self.icase = -1
        self.icase_fringe = None
        self.icase_disp = None
        self.icase_vector = None

    def _get_fringe_data(self, icase, scale=None):
        """helper for ``on_fringe``"""
        is_valid = False
        # (grid_result, name_tuple, name_str, data)
        failed_data = (None, None, None, None)

        if not isinstance(icase, integer_types):
            self.log_error('icase=%r is not an integer; type=%s' % (icase, type(icase)))
            return is_valid, failed_data
            #assert isinstance(icase, integer_types), icase

        try:
            case = self.result_cases[icase]
        except KeyError:
            self.log_error('icase=%r is not a valid result_case' % icase)
            return is_valid, failed_data

        label2 = ''
        (obj, (i, name)) = self.result_cases[icase]
        subcase_id = obj.subcase_id

        case = obj.get_result(i, name)
        if scale is None:
            scale = 1.0
        else:
            # we can have ints...
            #case *= scale
            case = np.multiply(case, scale, casting="unsafe")
        #else:
            # phase is not None
            #xyz_nominal, vector_data = obj.get_vector_result(i, name, phase)

        if case is None:
            # normal result
            self.res_widget.result_case_window.treeView.fringe.setChecked(False)
            #self.res_widget.result_case_window.treeView.disp.setChecked(False)
            #self.res_widget.result_case_window.treeView.vector.setChecked(False)

            #is_legend_shown = False
            self.set_normal_result(icase, name, subcase_id)
            # this is valid, but we want to skip out
            return is_valid, failed_data

        result_type = obj.get_title(i, name)
        data_format = obj.get_data_format(i, name)
        vector_size = obj.get_vector_size(i, name)
        location = obj.get_location(i, name)
        methods = obj.get_methods(i)
        #scale = obj.get_scale(i, name)
        phase = obj.get_phase(i, name)
        label2 = obj.get_header(i, name)
        out = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
        nlabels, labelsize, ncolors, colormap = out

        #normi = _get_normalized_data(self.result_cases[icase])
        normi = _get_normalized_data(case)
        imin = np.nanargmin(normi)
        imax = np.nanargmax(normi)

        #if min_value is None and max_value is None:
            #max_value = normi.max()
            #min_value = normi.min()

        #if min_value is None and max_value is None:
        min_value, max_value = obj.get_min_max(i, name)
        subtitle, label = self.get_subtitle_label(subcase_id)
        if label2:
            label += '; ' + label2

        #================================================
        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)

        vector_size = 1
        name_tuple = (vector_size, subcase_id, result_type, label, min_value, max_value, scale)
        name_str = self._names_storage.get_name_string(name)
        #return name, normi, vector_size, min_value, max_value, norm_value

        grid_result = self.set_grid_values(name_tuple, normi, vector_size, phase)

        data = FringeData(
            icase, result_type, location, min_value, max_value, norm_value,
            data_format, scale, methods,
            subcase_id, subtitle, label,
            nlabels, labelsize, ncolors, colormap,
            imin, imax,
        )
        is_valid = True
        return is_valid, (grid_result, name_tuple, name_str, data)

    def get_mapping_for_location(self, location):
        """helper method for ``export_case_data``"""
        if location == 'centroid':
            word = 'ElementID'
            eids_nids = self.element_ids
        elif location == 'node':
            word = 'NodeID'
            eids_nids = self.node_ids
        else:
            raise NotImplementedError(location)
        return word, eids_nids

    def _get_disp_data(self, icase, is_disp, stop_on_failure=False):
        """helper for ``on_disp``"""
        is_valid = False
        failed_data = (None, None, None, None)
        #is_checked = self.res_widget.get_checked()
        is_checked = False

        if not isinstance(icase, integer_types):
            self.log_error('icase=%r is not an integer; type=%s' % (icase, type(icase)))
            #self.res_widget.set_checked('normals', is_checked)

            return is_valid, failed_data
            #assert isinstance(icase, integer_types), icase

        try:
            case = self.result_cases[icase]
        except KeyError:
            self.log_error('icase=%r is not a valid result_case' % icase)
            return is_valid, failed_data

        label2 = ''
        (obj, (i, name)) = case
        subcase_id = obj.subcase_id
        #case = obj.get_result(i, name)  # TODO: remove
        #if case is None:
            ## normal result
            #self.log_error('icase=%r is not a displacement/force' % icase)
            #return is_valid, failed_data

        result_type = obj.get_title(i, name)
        vector_size = obj.get_vector_size(i, name)
        if vector_size == 1:
            msg = 'icase=%r is not a displacement/force' % icase
            self.log_error(msg)
            if stop_on_failure:
                raise ValueError(msg)
            return is_valid, failed_data
        location = obj.get_location(i, name)
        if is_disp and location != 'node':
            self.log_error('icase=%r is not a displacement-like (nodal vector) result' % icase)
            return is_valid, failed_data

        xyz_nominal, vector_data = obj.get_vector_result(i, name)
        if is_disp and xyz_nominal is None:
            self.log_error('icase=%r is not a displacement-like (nodal vector) result' % icase)
            return is_valid, failed_data

        methods = obj.get_methods(i)
        data_format = obj.get_data_format(i, name)
        scale = obj.get_scale(i, name)
        phase = obj.get_phase(i, name)
        label2 = obj.get_header(i, name)
        out = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
        nlabels, labelsize, ncolors, colormap = out

        # setup
        #-----
        # make results
        #if len(case.shape) == 1:
            #normi = case
        #else:
            #normi = norm(case, axis=1)

        #if min_value is None and max_value is None:
            #max_value = normi.max()
            #min_value = normi.min()

        #if min_value is None and max_value is None:
        min_value, max_value = obj.get_min_max(i, name)
        unused_subtitle, label = self.get_subtitle_label(subcase_id)
        if label2:
            label += '; ' + label2

        #================================================
        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)

        vector_size = 3
        name_tuple = (vector_size, subcase_id, result_type, label, scale)
        name_str = self._names_storage.get_name_string(name)
        #return name, normi, vector_size, min_value, max_value, norm_value

        #grid_result = self.set_grid_values(name_tuple, normi, vector_size)
        grid_result = None
        min_value = None
        max_value = None
        norm_value = None


        imin = None
        imax = None
        subtitle = None
        data = DispData(
            icase, result_type, location, min_value, max_value, norm_value,
            data_format, scale, phase, methods,
            subcase_id, subtitle, label,
            nlabels, labelsize, ncolors, colormap,
            xyz_nominal, vector_data,
            is_checked,
            imin, imax,
        )
        is_valid = True
        return is_valid, (grid_result, name_tuple, name_str, data)

    def on_fringe(self, icase, update_legend_window=True, show_msg=True,
                  stop_on_failure=False):
        """
        Sets the icase data to the active fringe

        Parameters
        ----------
        icase : int; default=None
            selects the icase
            None : defaults to self.icase+1

        """
        self.icase = icase
        is_valid, data = self._update_vtk_fringe(icase)
        if not is_valid:
            return is_valid

        icase = data.icase
        result_type = data.result_type
        location = data.location
        min_value = data.min_value
        max_value = data.max_value
        #norm_value = data.norm_value
        data_format = data.data_format
        scale = data.scale
        methods = data.methods
        nlabels = data.nlabels
        labelsize = data.labelsize
        ncolors = data.ncolors
        colormap = data.colormap
        imin = data.imin
        imax = data.imax
        #subcase_id = data.subcase_id
        #subtitle = data.subtitle
        #label = data.label

        #is_legend_shown = True
        #if is_legend_shown is None:
        self.show_legend()
        self.scalar_bar.is_shown = True
        is_legend_shown = self.scalar_bar.is_shown

        # TODO: normal -> fringe screws up legend
        #print('is_legend_shown = ', is_legend_shown)
        if not is_legend_shown:
            #print('showing')
            self.show_legend()

        self.update_contour_filter(nlabels, location, min_value, max_value)

        self.update_scalar_bar(result_type, min_value, max_value,
                               data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_shown=is_legend_shown)

        icase_fringe = icase
        icase_disp = self.icase_disp
        icase_vector = self.icase_vector

        phase = 0.0
        arrow_scale = 0.0
        if update_legend_window:
            self.legend_obj.update_legend(
                icase_fringe, icase_disp, icase_vector,
                result_type, min_value, max_value, data_format, scale, phase,
                arrow_scale,
                nlabels, labelsize, ncolors, colormap,
                use_disp_internal=True, use_vector_internal=True)
            self._set_legend_fringe(True)
        self.res_widget.update_method(methods)

        self.icase_fringe = icase
        self.grid.Modified()
        self.grid_selected.Modified()
        self._update_min_max_actors(location, icase_fringe,
                                    imin, min_value,
                                    imax, max_value)
        self.vtk_interactor.Render()
        self.res_widget.result_case_window.treeView.fringe.setChecked(True)
        is_valid = True
        return is_valid

    def show_hide_min_actor(self, render=True):
        """flips the status of the min label actor"""
        actor = self.min_max_actors[0]
        show_hide_actor(actor)
        if render:
            self.Render()

    def show_hide_max_actor(self, render=True):
        """flips the status of the max label actor"""
        actor = self.min_max_actors[1]
        show_hide_actor(actor)
        if render:
            self.Render()

    def _update_min_max_actors(self, location, icase_fringe,
                               imin, min_value, imax, max_value):
        """updates the values for the min and max actors"""
        settings = self.settings
        self.is_min_actor = True
        self.is_max_actor = True

        if location == 'node':
            if hasattr(self, 'xyz_cid0') and self.xyz_cid0 is None:
                return
            points = self.grid.GetPoints()
            xyz_min = np.array(points.GetPoint(imin), dtype='float64')
            xyz_max = np.array(points.GetPoint(imax), dtype='float64')
            xyzs = [xyz_min, xyz_max]
            #xyzs = [
            #    self.xyz_cid0[imin, :], # min
            #    self.xyz_cid0[imax, :], # max
            #]
        elif location == 'centroid':
            #print(self.models['main'].elements)
            #print('neids', len(self.models['main'].elements))
            xyzs = [
                self.cell_centroid_grid(self.grid, imin),
                self.cell_centroid_grid(self.grid, imax),
            ]
        else:  # pragma: no cover
            raise NotImplementedError(location)

        min_maxs = [
            (imin, min_value, xyzs[0]),
            (imax, max_value, xyzs[1]),
        ]

        for (imin_max, unused_value, xyz), text_actor in zip(min_maxs, self.min_max_actors):
            self.imin = imin
            self.imax = imax
            text_actor.SetPosition(*xyz)
            text_prop = text_actor.GetTextProperty()
            text_prop.SetFontSize(settings.annotation_size)
            text_prop.SetColor(settings.annotation_color)
            text_actor.Modified()

    def _update_vtk_fringe(self, icase, scale=None):
        """helper method for ``on_fringe``"""
        is_valid, (grid_result, name, name_str, data) = self._get_fringe_data(icase, scale)
        #print("is_valid=%s scale=%s" % (is_valid, scale))
        if not is_valid:
            return is_valid, data

        icase = data.icase
        location = data.location
        min_value = data.min_value
        max_value = data.max_value
        imin = data.imin
        imax = data.imax
        subcase_id = data.subcase_id
        subtitle = data.subtitle
        label = data.label

        #-----------------------------------
        grid = self.grid

        grid_result.SetName(name_str)
        self._names_storage.add(name)

        cell_data = grid.GetCellData()
        point_data = grid.GetPointData()
        if location == 'centroid':
            #cell_data.RemoveArray(name_str)
            self._names_storage.remove(name)
            cell_data.AddArray(grid_result)

            #if location != obj_location:
            point_data.SetActiveScalars(None)
            cell_data.SetActiveScalars(name_str)

        elif location == 'node':
            #point_data.RemoveArray(name_str)
            self._names_storage.remove(name)
            point_data.AddArray(grid_result)

            #if location != obj_location:
            cell_data.SetActiveScalars(None)
            point_data.SetActiveScalars(name_str)
        else:
            raise RuntimeError(location)

        self.tool_actions.update_text_actors(subcase_id, subtitle,
                                             imin, min_value,
                                             imax, max_value, label, location)
        is_valid = True
        return is_valid, data

    def on_disp(self, icase, apply_fringe=False, update_legend_window=True, show_msg=True,
                stop_on_failure=False):
        """Sets the icase data to the active displacement"""
        is_disp = True
        self._on_disp_vector(icase, is_disp, apply_fringe, update_legend_window, show_msg=show_msg,
                             stop_on_failure=stop_on_failure)
        self.res_widget.result_case_window.treeView.disp.setChecked(True)

    def on_vector(self, icase, apply_fringe=False, update_legend_window=True, show_msg=True,
                  stop_on_failure=False):
        """Sets the icase data to the active vector"""
        is_disp = False
        self._on_disp_vector(icase, is_disp, apply_fringe, update_legend_window, show_msg=show_msg,
                             stop_on_failure=stop_on_failure)
        self.res_widget.result_case_window.treeView.vector.setChecked(True)

    def _on_disp_vector(self, icase, is_disp, apply_fringe=False,
                        update_legend_window=True, show_msg=True, stop_on_failure=False):
        """
        Sets the icase data to the active displacement/vector

        Parameters
        ----------
        case : int; default=None
            selects the icase
            None : defaults to self.icase+1

        """
        self.icase = icase
        is_valid, (unused_grid_result, unused_name, unused_name_str, data) = self._get_disp_data(
            icase, is_disp, stop_on_failure=stop_on_failure)

        if not is_valid:
            return

        icase = data.icase
        scale = data.scale
        location = data.location
        min_value = data.min_value
        max_value = data.max_value
        xyz_nominal = data.xyz_nominal
        vector_data = data.vector_data
        #imin = data.imin
        #imax = data.imax
        #unused_subcase_id = data.subcase_id
        #subtitle = data.subtitle
        #unused_label = data.label

        #deflects = obj.deflects(i, res_name)

        #(obj, (obj_i, obj_name)) = self.result_cases[self.icase]
        #obj_location = obj.get_location(obj_i, obj_name)
        #obj_location = ''
        #-----------------------------------
        grid = self.grid

        #print('disp=%s location=%r' % (is_disp, location))
        if is_disp: # or obj.deflects(i, res_name):
            #grid_result1 = self.set_grid_values(name, case, 1)
            #point_data.AddArray(grid_result1)

            self._is_displaced = True
            self._is_forces = False
            self._xyz_nominal = xyz_nominal
            self._update_grid(vector_data)

            grid.Modified()
            self.grid_selected.Modified()
            self.icase_disp = icase
        else:
            self.icase_vector = icase
            if location == 'node':
                #self._is_displaced = False
                self._is_forces = True
                #xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
                self._update_forces(vector_data, set_scalars=apply_fringe, scale=scale)
            else:
                #self._is_displaced = False
                self._is_forces = True
                #xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
                self._update_elemental_vectors(vector_data, set_scalars=apply_fringe, scale=scale)
            self.icase_vector = icase

        icase_fringe = self.icase_fringe
        icase_disp = self.icase_disp
        icase_vector = self.icase_vector

        result_type = None
        max_value = None
        min_value = None
        data_format = None
        nlabels = None
        ncolors = None
        colormap = None
        labelsize = None

        scale = None
        phase = None
        arrow_scale = None
        if update_legend_window:
            self.legend_obj.update_legend(
                icase_fringe, icase_disp, icase_vector,
                result_type, min_value, max_value, data_format, scale, phase,
                arrow_scale,
                nlabels, labelsize, ncolors, colormap, use_fringe_internal=True,
                use_disp_internal=True, use_vector_internal=True,
                external_call=False)

        self.vtk_interactor.Render()


    def cycle_results_explicit(self, case=None, explicit=True, min_value=None, max_value=None,
                               show_msg=True):
        """
        Forces the result to cycle regardless of whether or not the
        icase value is the same.  You'd do this when you've just
        loaded a model.

        """
        assert case is not False, case
        #if explicit:
            #self.log_command('cycle_results(case=%r)' % case)
        found_cases = self.increment_cycle(case)
        if found_cases:
            icase = self._set_case(case, self.icase, explicit=explicit, cycle=True,
                                   min_value=min_value, max_value=max_value,
                                   show_msg=show_msg)
            assert icase is not False, case
        else:
            icase = None
        return icase

    def get_name_result_data(self, icase):
        key = self.case_keys[icase]
        assert isinstance(key, integer_types), key
        (obj, (i, name)) = self.result_cases[key]
        #subcase_id = obj.subcase_id
        case = obj.get_result(i, name)
        return name, case

    def delete_cases(self, icases_to_delete, ask=True):
        """
        Used by the Sidebar to delete results

        Parameters
        ----------
        icases_to_delete : List[int]
            the result cases to delete
        ask : bool; default=True
            TODO: does nothing...

        """
        for icase in icases_to_delete:
            if icase not in self.case_keys:
                continue
            self.case_keys.remove(icase)
            del self.result_cases[icase]

            # we leave this, so we can still cycle the results without renumbering the cases
            #self.ncases -= 1
        self.res_widget.set_case_keys(self.case_keys)
        self.log_command('delete_cases(icases_to_delete=%s, ask=%s)' % (icases_to_delete, ask))

    def _get_sidebar_data(self, unused_name: str):
        """
        gets the form for the selected name

        Parameters
        ----------
        name : str
            the name that was selected

        Returns
        -------
        form : List[tuple]
            the form data

        """
        return []

    def _set_case(self, unused_result_name, icase, explicit=False, cycle=False,
                  skip_click_check=False, min_value=None, max_value=None,
                  is_legend_shown=None, show_msg=True):
        """
        Internal method for doing results updating

        Parameters
        ----------
        unused_result_name : str
            the name of the case for debugging purposes
        icase : int
            the case id
        cycle : bool; default=False
            ???
        skip_click_check : bool; default=False
            There is no reason to update if the case didn't change on the
            Results Sidebar
            True  : Legend Menu
            False : Results Sidebar
        min_value : float; default=None
            the min value
            None : use the default
        max_value  : float; default=None
            the max value
            None : use the default
        is_legend_shown : bool; default=None
            is the scalar bar shown
            None : stick with the current value
            True : show the legend
            False : hide the legend
        explicit : bool; default=False
            show the command when we're doing in the log
        show_msg : bool; default=True
            ???

        """
        _update_icase = (
            self.icase != self.icase_fringe and self.icase_fringe is not None or
            self.icase != self.icase_disp  and self.icase_disp is not None or
            self.icase != self.icase_vector and self.icase_vector is not None
        )
        if _update_icase:
            skip_click_check = True

        if not skip_click_check:
            if not cycle and icase == self.icase:
                # don't click the button twice
                # cycle=True means we're cycling
                # cycle=False skips that check
                return None

        try:
            key = self.case_keys[icase]
        except (KeyError, TypeError):
            print('icase=%r case_keys=%s' % (icase, str(self.case_keys)))
            raise
        icase = key
        self.icase = icase
        #assert self.icase >= 0, 'icase=%i key=%s case_keys=%s' % (icase, key, self.case_keys)


        # these will be fixed later in this function
        self.icase_fringe = None
        self.icase_disp = None
        self.icase_vector = None

        case = self.result_cases[key]
        label2 = ''
        assert isinstance(key, integer_types), key
        (obj, (i, name)) = self.result_cases[key]
        subcase_id = obj.subcase_id
        case = obj.get_result(i, name)
        result_type = obj.get_title(i, name)
        vector_size0 = obj.get_vector_size(i, name)
        location = obj.get_location(i, name)
        methods = obj.get_methods(i)
        data_format = obj.get_data_format(i, name)
        scale = obj.get_scale(i, name)
        phase = obj.get_phase(i, name)
        label2 = obj.get_header(i, name)
        out = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
        nlabels, labelsize, ncolors, colormap = out
        #default_max, default_min = obj.get_default_min_max(i, name)
        if min_value is None or max_value is None:
            min_valuei, max_valuei = obj.get_min_max(i, name)
            if min_value is None:
                min_value = min_valuei
            if max_value is None:
                max_value = max_valuei

        complex_types = ['complex64']
        if hasattr(max_value, 'dtype') and max_value.dtype.name in complex_types:
            raise TypeError(max_value)
        if hasattr(min_value, 'dtype') and min_value.dtype.name in complex_types:
            raise TypeError(min_value)

        subtitle, label = self.get_subtitle_label(subcase_id)
        if label2:
            label += '; ' + label2
        #print("subcase_id=%s result_type=%r subtitle=%r label=%r"
              #% (subcase_id, result_type, subtitle, label))

        #================================================
        if case is None:
            self.icase_fringe = icase
            self.set_normal_result(icase, name, subcase_id)
            return icase

        elif not self._is_fringe:
            self.icase_fringe = icase
            # we maybe hacked the scalar bar to turn off for Normals/Clear Results
            # so we turn it back on
            self.show_legend()
            self._set_legend_fringe(True)

        if len(case.shape) == 1:
            normi = case
        else:
            normi = norm(case, axis=1)

        #print(case, len(case))
        imin = np.nanargmin(normi)
        imax = np.nanargmax(normi)
        #print(imin, imax, normi[imin], normi[imax])

        #if min_value is None and max_value is None:
            #max_value = normi.max()
            #min_value = normi.min()

        #================================================
        # flips sign to make colors go from blue -> red
        #norm_value = float(max_value - min_value)

        vector_size = 1
        name = (vector_size, subcase_id, result_type, label, scale)
        if self._names_storage.has_exact_name(name):
            grid_result = None
        else:
            grid_result = self.set_grid_values(name, normi, vector_size, phase)

        if vector_size0 == 1:
            name_vector = None
            grid_result_vector = None
        else:
            vector_size = 3
            name_vector = (vector_size, subcase_id, result_type, label, scale)
            if self._names_storage.has_exact_name(name_vector):
                grid_result_vector = None
            else:
                grid_result_vector = self.set_grid_values(name_vector, case, vector_size, phase)

        self.final_grid_update(icase, name, grid_result,
                               name_vector, grid_result_vector,
                               key, subtitle, label,
                               min_value, max_value, show_msg)

        if is_legend_shown is None:
            is_legend_shown = self.scalar_bar.is_shown

        self.update_scalar_bar(result_type, min_value, max_value,
                               data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_shown=is_legend_shown)

        icase_fringe = icase
        icase_disp = self.icase_disp
        icase_vector = self.icase_vector

        #print('norm =', normi, location)
        self._update_min_max_actors(location, icase_fringe,
                                    imin, min_value,
                                    imax, max_value)

        self.update_text_actors(subcase_id, subtitle,
                                imin, min_value,
                                imax, max_value, label, location)


        arrow_scale = 0.0
        self.legend_obj.update_legend(
            icase_fringe, icase_disp, icase_vector,
            result_type, min_value, max_value, data_format, scale, phase,
            arrow_scale,
            nlabels, labelsize, ncolors, colormap, use_fringe_internal=True,
            external_call=False)

        # updates the type of the result that is displayed
        # method:
        #     for a nodeID, the method is [node]
        #     for an elementID, the method is [centroid]
        #     for a displacement, the methods are [magnitude, tx, ty, tz, rx, ry, rz]
        # location:
        #     for a nodeID, the location is [node]
        #     for an elementID, the location is [centroid]
        #     for a displacement, the location is [node]

        #location = self.get_case_location(key)
        self.res_widget.update_method(methods)
        if explicit:
            self.log_command('cycle_results(case=%r)' % self.icase)
        assert self.icase is not False, self.icase
        return self.icase

    def set_normal_result(self, icase, name, unused_subcase_id):
        """plots a NormalResult"""
        unused_name_str = self._names_storage.get_name_string(name)
        prop = self.geom_actor.GetProperty()

        prop.SetColor(RED)

        # the backface property is null
        #back_prop = self.geom_actor.GetBackfaceProperty()
        back_prop = vtk.vtkProperty()
        back_prop.SetColor(BLUE)
        self.geom_actor.SetBackfaceProperty(back_prop)

        grid = self.grid
        if self._is_displaced:
            self._is_displaced = False
            self._update_grid(self._xyz_nominal)
            self.icase_disp = None

        if self._is_forces:
            self.arrow_actor.SetVisibility(False)
            self.icase_vector = None

        cell_data = grid.GetCellData()
        cell_data.SetActiveScalars(None)

        point_data = grid.GetPointData()
        point_data.SetActiveScalars(None)
        self.icase_fringe = icase

        #if is_legend_shown is None:
            #is_legend_shown = self.scalar_bar.is_shown
        #self.update_scalar_bar(result_type, min_value, max_value,
                               #data_format,
                               #nlabels=nlabels, labelsize=labelsize,
                               #ncolors=ncolors, colormap=colormap,
                               #is_shown=is_legend_shown)
        #scale = 0.0
        #phase = None

        #min_value = -1.
        #max_value = 1.
        #icase_fringe = icase
        #self.legend_obj.update_legend(
            #icase_fringe, icase_disp, icase_vector,
            #result_type, min_value, max_value, data_format, scale, phase,
            #nlabels, labelsize, ncolors, colormap, external_call=False)
        self.hide_legend()
        self.scalar_bar.is_shown = False
        self._set_legend_fringe(False)
        #min_value = 'Front Face'
        #max_value = 'Back Face'
        #self.update_text_actors(subcase_id, subtitle,
                                #min_value, max_value, label)
        self.vtk_interactor.Render()

    def set_grid_values(self, name, case, vector_size: int, phase: float):
        """
        https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
        """
        if self._names_storage.has_exact_name(name):
            return None
        #print('name=%r case=%r' % (name, case))

        if not hasattr(case, 'dtype'):
            raise RuntimeError('name=%s case=%s' % (name, case))

        if issubdtype(case.dtype, np.integer):
            data_type = vtk.VTK_INT
            self.grid_mapper.InterpolateScalarsBeforeMappingOn()
        elif issubdtype(case.dtype, np.floating):
            data_type = vtk.VTK_FLOAT
            self.grid_mapper.InterpolateScalarsBeforeMappingOff()
        elif case.dtype.name in ['complex64']:
            if phase:
                phaser = np.radians(phase)
                case = (np.cos(phaser) * case.real + np.sin(phaser) * case.imag).real
            else:
                case = case.real
            data_type = vtk.VTK_FLOAT
            self.grid_mapper.InterpolateScalarsBeforeMappingOff()
        else:
            raise NotImplementedError(case.dtype.type)


        #if 0: # nan testing
            #if case.dtype.name == 'float32':
                #case[50] = np.float32(1) / np.float32(0)
            #else:
                #case[50] = np.int32(1) / np.int32(0)

        if vector_size == 1:
            if case.flags.contiguous:
                case2 = case
            else:
                case2 = deepcopy(case)
            grid_result = numpy_to_vtk(
                num_array=case2,
                deep=True,
                array_type=data_type
            )
            #print('grid_result = %s' % grid_result)
            #print('max2 =', grid_result.GetRange())
        else:
            # vector_size=3
            if case.flags.contiguous:
                case2 = case
            else:
                case2 = deepcopy(case)
            grid_result = numpy_to_vtk(
                num_array=case2,
                deep=True,
                array_type=data_type
            )
        return grid_result

    def update_grid_by_icase_scale_phase(self, icase, scale, phase=0.0):
        """
        Updates to the deflection state defined by the cases

        Parameters
        ----------
        icase : int
            result number in self.result_cases
        scale : float
            deflection scale factor; true scale
        phase : float; default=0.0
            phase angle (degrees); unused for real results

        Returns
        -------
        xyz : (nnodes, 3) float ndarray
            the nominal state
        deflected_xyz : (nnodes, 3) float ndarray
            the deflected state

        """
        #print('update_grid_by_icase_scale_phase')
        (obj, (i, res_name)) = self.result_cases[icase]
        xyz_nominal, vector_data = obj.get_vector_result_by_scale_phase(
            i, res_name, scale, phase)

        #grid_result1 = self.set_grid_values(name, case, 1,
            #min_value, max_value, norm_value)
        #point_data.AddArray(grid_result1)

        self._is_displaced = True
        self._xyz_nominal = xyz_nominal
        self._update_grid(vector_data)

    def update_forces_by_icase_scale_phase(self, icase, arrow_scale, phase=0.0):
        """
        Updates to the force state defined by the cases

        Parameters
        ----------
        icase : int
            result number in self.result_cases
        arrow_scale : float
            force scale factor; ??? scale
        phase : float; default=0.0
            phase angle (degrees); unused for real results

        """
        #print('update_grid_by_icase_scale_phase')
        (obj, (i, res_name)) = self.result_cases[icase]
        unused_xyz_nominal, vector_data = obj.get_vector_result_by_scale_phase(
            i, res_name, arrow_scale, phase)

        #grid_result1 = self.set_grid_values(name, case, 1)
        #point_data.AddArray(grid_result1)

        self._is_forces = True
        ## TODO: support elemental forces
        self._update_forces(vector_data, set_scalars=False, scale=arrow_scale)
        #self._update_elemental_vectors(forces_array, set_scalars=True, scale=None)

    def final_grid_update(self, icase, name, grid_result,
                          name_vector, grid_result_vector,
                          key, subtitle, label,
                          min_value, max_value, show_msg):
        assert isinstance(key, integer_types), key
        (obj, (i, res_name)) = self.result_cases[key]
        subcase_id = obj.subcase_id
        result_type = obj.get_title(i, res_name)
        vector_size = obj.get_vector_size(i, res_name)
        location = obj.get_location(i, res_name)
        out = obj.get_nlabels_labelsize_ncolors_colormap(i, res_name)
        nlabels, unused_labelsize, unused_ncolors, unused_colormap = out

        #if vector_size == 3:
            #print('name, grid_result, vector_size=3', name, grid_result)
        self._final_grid_update(icase, name, grid_result, None, None, None,
                                1, subcase_id, result_type, location, subtitle, label,
                                revert_displaced=True, show_msg=show_msg)
        self.update_contour_filter(nlabels, location, min_value, max_value)

        if vector_size == 3:
            self._final_grid_update(icase, name_vector, grid_result_vector, obj, i, res_name,
                                    vector_size, subcase_id, result_type, location, subtitle, label,
                                    revert_displaced=False, show_msg=show_msg)
            #xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
            #self._update_grid(vector_data)

    def _final_grid_update(self, icase, name, grid_result, obj, i, res_name,
                           vector_size, subcase_id, result_type, location, subtitle, label,
                           revert_displaced=True, show_msg=True):
        if name is None:
            return
        assert icase >= 0

        # the result type being currently shown
        # for a Nastran NodeID/displacement, this is 'node'
        # for a Nastran ElementID/PropertyID, this is 'element'
        self.result_location = location

        grid = self.grid
        name_str = self._names_storage.get_name_string(name)
        if not self._names_storage.has_exact_name(name):
            grid_result.SetName(name_str)
            self._names_storage.add(name)

            if self._is_displaced and revert_displaced:
                self._is_displaced = False
                self._update_grid(self._xyz_nominal)

            if self._is_forces:
                self.arrow_actor.SetVisibility(False)

            if location == 'centroid':
                cell_data = grid.GetCellData()
                if self._names_storage.has_close_name(name):
                    cell_data.RemoveArray(name_str)
                    self._names_storage.remove(name)

                cell_data.AddArray(grid_result)
                if show_msg:
                    self.log_info('centroidal plotting vector=%s - subcase_id=%s '
                                  'result_type=%s subtitle=%s label=%s'
                                  % (vector_size, subcase_id, result_type, subtitle, label))
            elif location == 'node':
                point_data = grid.GetPointData()
                if self._names_storage.has_close_name(name):
                    point_data.RemoveArray(name_str)
                    self._names_storage.remove(name)

                if vector_size == 1:
                    if show_msg:
                        self.log_info('node plotting vector=%s - subcase_id=%s '
                                      'result_type=%s subtitle=%s label=%s"'
                                      % (vector_size, subcase_id, result_type, subtitle, label))
                    point_data.AddArray(grid_result)
                elif vector_size == 3:
                    #print('vector_size3; get, update')

                    xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
                    if obj.deflects(i, res_name):
                        #grid_result1 = self.set_grid_values(name, case, 1,
                                                            #min_value, max_value, norm_value)
                        #point_data.AddArray(grid_result1)

                        self._is_displaced = True
                        self._is_forces = False
                        self._xyz_nominal = xyz_nominal
                        self._update_grid(vector_data)
                        self.icase_disp = icase
                    else:
                        self._is_displaced = False
                        self._is_forces = True
                        scale = obj.get_scale(i, res_name)
                        xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
                        self._update_forces(vector_data, scale)
                        self.icase_vector = icase

                    if show_msg:
                        self.log_info('node plotting vector=%s - subcase_id=%s '
                                      'result_type=%s subtitle=%s label=%s'
                                      % (vector_size, subcase_id, result_type, subtitle, label))
                    #point_data.AddVector(grid_result) # old
                    #point_data.AddArray(grid_result)
                else:
                    raise RuntimeError(vector_size)
            else:
                raise RuntimeError(location)

        if location == 'centroid':
            self.icase_fringe = icase
            cell_data = grid.GetCellData()
            cell_data.SetActiveScalars(name_str)

            point_data = grid.GetPointData()
            point_data.SetActiveScalars(None)
            if vector_size == 1:
                #point_data.SetActiveVectors(None)   # I don't think I need this
                pass
            else:
                raise RuntimeError(vector_size)
        elif location == 'node':
            cell_data = grid.GetCellData()
            cell_data.SetActiveScalars(None)

            point_data = grid.GetPointData()
            if vector_size == 1:
                self.icase_fringe = icase
                point_data.SetActiveScalars(name_str)  # TODO: None???
            elif vector_size == 3:
                pass
                #point_data.SetActiveVectors(name_str)
            else:
                raise RuntimeError(vector_size)
            #print('name_str=%r' % name_str)
        else:
            raise RuntimeError(location)

        grid.Modified()
        self.grid_selected.Modified()
        #self.contour_filter.Modified()
        #self.update_all()
        #self.update_all()
        if len(self.groups):
            self.post_group_by_name(self.group_active)

        self.vtk_interactor.Render()

        self.hide_labels(show_msg=False)
        self.show_labels(case_keys=[self.icase], show_msg=False)

    def _update_forces(self, forces_array, set_scalars=True, scale=None):
        """changes the glyphs"""
        grid = self.grid
        if scale is not None:
            self.glyphs.SetScaleFactor(self.glyph_scale_factor * scale)
        new_forces, mag = normalize_forces(forces_array)

        vtk_vectors = numpy_to_vtk(new_forces, deep=1)
        grid.GetPointData().SetVectors(vtk_vectors)
        if set_scalars:
            vtk_mag = numpy_to_vtk(mag, deep=1)
            grid.GetPointData().SetScalars(vtk_mag)
            grid.GetCellData().SetScalars(None)
        self.arrow_actor.SetVisibility(True)
        grid.Modified()
        self.grid_selected.Modified()

    def _update_elemental_vectors(self, forces_array, set_scalars=True, scale=None):
        """changes the glyphs"""
        grid = self.grid
        if scale is not None:
            # TODO: glyhs_centroid?
            self.glyphs_centroid.SetScaleFactor(self.glyph_scale_factor * scale)
        new_forces, mag = normalize_forces(forces_array)

        vtk_vectors = numpy_to_vtk(new_forces, deep=1)
        grid.GetPointData().SetVectors(None)
        #print('_update_elemental_vectors; shape=%s' % (str(new_forces.shape)))
        grid.GetCellData().SetVectors(vtk_vectors)
        if set_scalars:
            vtk_mag = numpy_to_vtk(mag, deep=1)
            grid.GetPointData().SetScalars(None)
            grid.GetCellData().SetScalars(vtk_mag)
        self.arrow_actor_centroid.SetVisibility(True)
        grid.Modified()
        self.grid_selected.Modified()

    def _update_grid(self, nodes):
        """deflects the geometry"""
        grid = self.grid
        points = grid.GetPoints()
        #inan = np.where(nodes.ravel() == np.nan)[0]
        #if len(inan) > 0:
            #raise RuntimeError('nan in nodes...')
        numpy_to_vtk_points(nodes, points=points, dtype='<f', deep=1)
        grid.Modified()
        self.grid_selected.Modified()
        self._update_follower_grids(nodes)
        self._update_follower_grids_complex(nodes)

    def _update_follower_grids(self, nodes):
        """updates grids that use the same ids as the parent model"""
        for name, nids in self.follower_nodes.items():
            grid = self.alt_grids[name]
            points = grid.GetPoints()
            for j, nid in enumerate(nids):
                i = self.nid_map[nid]
                points.SetPoint(j, *nodes[i, :])
            grid.Modified()

    def _update_follower_grids_complex(self, nodes):
        """updates grids that use a complicated update method"""
        for name, follower_function in self.follower_functions.items():
            grid = self.alt_grids[name]
            points = grid.GetPoints()
            follower_function(self.nid_map, grid, points, nodes)
            grid.Modified()

    def _get_icase(self, result_name):
        if not self.result_cases:
            raise IndexError('result_name=%r not in self.result_cases' % result_name)

        #print('result_cases.keys() =', self.result_cases.keys())
        i = 0
        icase = None
        for icase in sorted(self.result_cases.keys()):
            #cases = self.result_cases[icase]
            if result_name == icase[1]:
                unused_found_case = True
                icase = i
                break
            i += 1
        return icase

    def increment_cycle(self, icase=None):
        #print('1-icase=%r self.icase=%s ncases=%r' % (icase, self.icase, self.ncases))
        #print(type(icase))
        if isinstance(icase, integer_types):
            self.icase = icase
            if self.icase >= self.ncases:
                self.icase = 0
        else:
            self.icase += 1
            if self.icase == self.ncases:
                self.icase = 0
        #print('2-icase=%r self.icase=%s ncases=%r' % (icase, self.icase, self.ncases))

        if self.ncases > 0:
            #key = self.case_keys[self.icase]
            #try:
                #key = self.case_keys[self.icase]
            #except IndexError:
                #found_cases = False
                #return found_cases
            #except TypeError:
                #msg = 'type(case_keys)=%s\n' % type(self.case_keys)
                #msg += 'icase=%r\n' % str(self.icase)
                #msg += 'case_keys=%r' % str(self.case_keys)
                #print(msg)
                #raise TypeError(msg)
            msg = 'icase=%r\n' % str(self.icase)
            msg += 'case_keys=%r' % str(self.case_keys)
            #print(msg)
            found_cases = True
        else:
            # key = self.case_keys[self.icase]
            # location = self.get_case_location(key)
            self.log_error('No results found.')
            self.scalar_bar_actor.SetVisibility(False)
            found_cases = False
        #print("next icase=%s key=%s" % (self.icase, key))
        return found_cases

    def get_result_name(self, icase):
        assert isinstance(icase, integer_types), icase
        (unused_obj, (unused_i, res_name)) = self.result_cases[icase]
        return res_name

    def get_case_location(self, icase):
        assert isinstance(icase, integer_types), icase
        (obj, (i, res_name)) = self.result_cases[icase]
        return obj.get_location(i, res_name)

    #---------------------------------------------------------------------------
    def hide_labels(self, case_keys=None, show_msg=True):
        if case_keys is None:
            names = 'None)  # None -> all'
            case_keys = sorted(self.label_actors.keys())
        else:
            mid = '%s,' * len(case_keys)
            names = '[' + mid[:-1] + '])'

        count = 0
        for icase in case_keys:
            actors = self.label_actors[icase]
            for actor in actors:
                actor.VisibilityOff()
                #prop = actor.GetProperty()
                count += 1
        if count and show_msg:
            self.log_command('hide_labels(%s)' % names)

    def show_labels(self, case_keys=None, show_msg=True):
        if case_keys is None:
            names = 'None)  # None -> all'
            case_keys = sorted(self.label_actors.keys())
        else:
            mid = '%s,' * len(case_keys)
            names = mid[:-1] % case_keys + ')'

        count = 0
        for icase in case_keys:
            assert icase >= 0, case_keys
            try:
                actors = self.label_actors[icase]
            except KeyError:
                keys = list(self.label_actors.keys())
                msg = 'Cant find label_actors for icase=%r; keys=%s' % (
                    icase, keys)
                self.log.error(msg)
                continue
            for actor in actors:
                actor.VisibilityOn()
                count += 1
        if count and show_msg:
            # yes the ) is intentionally left off because it's already been added
            self.log_command('show_labels(%s)' % names)

    def remove_alt_grid(self, name, remove_geometry_property=False):
        if name in self.alt_grids:
            del self.alt_grids[name]
        if remove_geometry_property and name in self.geometry_properties:
            slot = self.geometry_properties[name].label_actors
            for sloti in slot:
                self.rend.RemoveActor(sloti)
            del self.geometry_properties[name]

    def reset_label_actors(self, name):
        slot = self.geometry_properties[name].label_actors
        for sloti in slot:
            self.rend.RemoveActor(sloti)
        return slot

    def create_alternate_vtk_grid(self, name, color=None, line_width=5, opacity=1.0, point_size=1,
                                  bar_scale=0.0, representation=None, display=None, is_visible=True,
                                  follower_nodes=None, follower_function=None,
                                  is_pickable=False, ugrid=None):
        """
        Creates an AltGeometry object

        Parameters
        ----------
        line_width : int
            the width of the line for 'surface' and 'main'
        color : [int, int, int]
            the RGB colors
        opacity : float
            0.0 -> solid
            1.0 -> transparent
        point_size : int
            the point size for 'point'
        bar_scale : float
            the scale for the CBAR / CBEAM elements
        representation : str
            main - change with main mesh
            wire - always wireframe
            point - always points
            surface - always surface
            bar - can use bar scale
        is_visible : bool; default=True
            is this actor currently visable
        is_pickable : bool; default=False
            can you pick a node/cell on this actor
        follower_nodes : List[int]
            the nodes that are brought along with a deflection
        follower_function : function
            a custom follower_node update function
        ugrid : vtk.vtkUnstructuredGrid(); default=None
            the grid object; one will be created that you can fill
            if None is passed in

        """
        if ugrid is None:
            ugrid = vtk.vtkUnstructuredGrid()
        self.alt_grids[name] = ugrid

        if name not in self.geometry_properties:
            self.geometry_properties[name] = AltGeometry(
                self, name, color=color,
                line_width=line_width, opacity=opacity,
                point_size=point_size, bar_scale=bar_scale,
                representation=representation, display=display,
                is_visible=is_visible, is_pickable=is_pickable)
        if follower_nodes is not None:
            self.follower_nodes[name] = follower_nodes
        if follower_function is not None:
            self.follower_functions[name] = follower_function

    def duplicate_alternate_vtk_grid(self, name, name_duplicate_from, color=None, line_width=5,
                                     opacity=1.0, point_size=1, bar_scale=0.0, is_visible=True,
                                     follower_nodes=None, is_pickable=False):
        """
        Copies the VTK actor

        Parameters
        ----------
        line_width : int
            the width of the line for 'surface' and 'main'
        color : [int, int, int]
            the RGB colors
        opacity : float
            0.0 -> solid
            1.0 -> transparent
        point_size : int
            the point size for 'point'
        bar_scale : float
            the scale for the CBAR / CBEAM elements
        is_visible : bool; default=True
            is this actor currently visable
        is_pickable : bool; default=False
            can you pick a node/cell on this actor
        follower_nodes : List[int]
            the nodes that are brought along with a deflection

        """
        self.alt_grids[name] = vtk.vtkUnstructuredGrid()
        if name_duplicate_from == 'main':
            grid_copy_from = self.grid
            representation = 'toggle'
        else:
            grid_copy_from = self.alt_grids[name_duplicate_from]
            props = self.geometry_properties[name_duplicate_from]
            representation = props.representation
        self.alt_grids[name].DeepCopy(grid_copy_from)

        #representation : str
            #main - change with main mesh
            #wire - always wireframe
            #point - always points
            #surface - always surface
            #bar - can use bar scale
        self.geometry_properties[name] = AltGeometry(
            self, name, color=color, line_width=line_width,
            opacity=opacity, point_size=point_size,
            bar_scale=bar_scale, representation=representation,
            is_visible=is_visible, is_pickable=is_pickable)

        if follower_nodes is not None:
            self.follower_nodes[name] = follower_nodes

    def _create_plane_actor_from_points(self, p1, p2, i, k, dim_max,
                                        actor_name='plane'):
        """
        This is used by the cutting plane tool and the shear/moment/torque tool.

           4+------+3
            |      |
            p1     p2
            |      |
           1+------+2

        """
        shift = 1.1
        dshift = (shift - 1) / 2.
        half_shift = 0.5 + dshift
        delta = half_shift * dim_max
        #dim_xy = shift * dim_max

        #n1 = 1 - dim_max * (dshift * i + half_shift * k)
        #n2 = n1 + shift * dim_max * i
        #n3 = n2 + shift * dim_max * k
        #n4 = n1 + shift * dim_max * k
        pcenter = (p1 + p2) / 2
        n1 = pcenter - delta * i - delta * k
        n2 = pcenter + delta * i - delta * k
        n3 = pcenter + delta * i + delta * k
        n4 = pcenter - delta * i + delta * k

        x = np.linspace(0., 1., num=10)
        y = x
        if actor_name in self.alt_grids:
            plane_actor = self.plane_actor
            add = False
            #alt_grid =
            #plane_source = vtk.vtkPlaneSource()
            #self.rend.AddActor(plane_actor)
            #self.plane_actor = plane_actor
        else:
            add = True
            alt_grid = vtk.vtkUnstructuredGrid()
            self.alt_grids[actor_name] = alt_grid

            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputData(alt_grid)
            plane_actor = vtk.vtkActor()
            plane_actor.SetMapper(mapper)

            #plane_source = self.plane_source
            #plane_actor = self.plane_actor
            self.plane_actor = plane_actor
            self.rend.AddActor(plane_actor)

        nodes, elements = points_elements_from_quad_points(n1, n2, n3, n4, x, y)
        self.set_quad_grid(actor_name, nodes, elements, color=RED,
                           line_width=1, opacity=1., representation='surface',
                           add=add)
        #plane_actor.Modified()
        return plane_actor

    def _create_point_actor_from_points(self, points, point_size=8,
                                        actor_name='plane_points'):  # pragma: no cover
        """
        This is used by the shear/moment/torque tool.

            p1------p2

        """
        points = np.asarray(points)

        if actor_name in self.alt_grids:
            point_actor = self.point_actor
            #alt_grid =
            #plane_source = vtk.vtkPlaneSource()
            #self.rend.AddActor(point_actor)
            #self.point_actor = point_actor
        else:
            alt_grid = vtk.vtkUnstructuredGrid()
            self.alt_grids[actor_name] = alt_grid

            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputData(alt_grid)
            point_actor = vtk.vtkActor()
            point_actor.SetMapper(mapper)

            #plane_source = self.plane_source
            #point_actor = self.point_actor
            self.point_actor = point_actor
            self.rend.AddActor(point_actor)

        ## TODO: not done...
        #nodes, elements = points_elements_from_quad_points(n1, n2, n3, n4, x, y)
        #color = RED
        #self.set_quad_grid(actor_name, nodes, elements, color,
                           #line_width=1, opacity=1., add=False)
        return point_actor

    def _make_contour_filter(self):  # pragma: no cover
        """trying to make model lines...doesn't work"""
        if not self.make_contour_filter:
            return

        self.contour_filter = vtk.vtkContourFilter()

        #if 0:
            # doesn't work...in progress
            #geometry_filter = vtk.vtkGeometryFilter()
            #geometry_filter.SetInputData(self.grid_selected)
            #geometry_filter.Update()
            #poly_data = geometry_filter.GetOutput()

            #self.contour_filter.SetInputData(poly_data)
        #elif 0:  # pragma: no cover
            # doesn't work
            #self.contour_filter.SetInputData(self.grid_selected)
        if 1:
            # https://blog.kitware.com/cell-set-as-unstructured-grid/
            self.contour_filter.SetInputData(self.grid)
        else:
            raise RuntimeError('invalid contour_filter option')
        #self.contour_filter.GenerateValues(1, 10, 10)
        #self.contour_filter.SetComputeScalars(1)
        #contour_filter.SetInputConnection(self.grid_selected.GetOutputPort())
        #self.contour_filter.SetInputData(None)
        self.contour_filter.ComputeScalarsOff()
        self.contour_filter.ComputeNormalsOff()


        # Connect the segments of the conours into polylines
        contour_stripper = vtk.vtkStripper()
        contour_stripper.SetInputConnection(self.contour_filter.GetOutputPort())
        contour_stripper.Update()

        number_of_contour_lines = contour_stripper.GetOutput().GetNumberOfLines()
        print('There are %s contours lines.' % number_of_contour_lines)

        include_labels = False
        if include_labels:
            label_points = contour_stripper.GetOutput().GetPoints()
            unused_cells = contour_stripper.GetOutput().GetLines()
            label_scalars = contour_stripper.GetOutput().GetPointData().GetScalars()

            label_poly_data.SetPoints(label_points)
            label_poly_data.GetPointData().SetScalars(label_scalars)

            # The labeled data mapper will place labels at the points
            label_mapper = vtk.vtkLabeledDataMapper()
            label_mapper.SetFieldDataName("Isovalues")
            label_mapper.SetInputData(label_poly_data)

            label_mapper.SetLabelModeToLabelScalars()
            label_mapper.SetLabelFormat("%6.2f")

            label_mapper.SetLabelModeToLabelScalars()
            label_mapper.SetLabelFormat("%6.2f")

            isolabels_actor = vtk.vtkActor2D()
            isolabels_actor.SetMapper(label_mapper)

        contour_mapper = vtk.vtkPolyDataMapper()
        contour_mapper.SetInputConnection(contour_stripper.GetOutputPort())
        contour_mapper.ScalarVisibilityOff()

        isolines_actor = vtk.vtkActor()
        isolines_actor.SetMapper(contour_mapper)
        isolines_actor.GetProperty().SetColor(0., 0., 0.)

        # Add the actors to the scene
        self.rend.AddActor(isolines_actor)
        if include_labels:
            self.rend.AddActor(isolabels_actor)

        self.contour_mapper = contour_mapper
        self.contour_stripper = contour_stripper
        self.contour_lines_actor = isolines_actor
        self.contour_lines_actor.VisibilityOff()

    def update_contour_filter(self, nlabels, location, min_value=None, max_value=None):
        """update the contour lines"""
        if not self.make_contour_filter:
            return
        if nlabels is None:
            nlabels = 11

        if location == 'centroid': # node/centroid
            self.contour_lines_actor.VisibilityOff()
            #self.contour_mapper.SetScalarModeToUseCellData()
            #res_data = self.grid.GetCellData().GetScalars()
            return
        elif location == 'node':
            #self.contour_mapper.SetScalarModeToUsePointData()
            res_data = self.grid.GetPointData().GetScalars()
        else:
            raise RuntimeError('location=%r' % location)

        self.contour_lines_actor.VisibilityOn()
        number_of_cuts = nlabels - 1
        if min_value is None or max_value is None:
            data_range = res_data.GetRange()
            if min_value is None:
                min_value = data_range[0]
            if max_value is None:
                max_value = data_range[1]

        self.contour_filter.GenerateValues(
            number_of_cuts,
            0.99 * min_value,
            0.99 * max_value)

        self.contour_filter.Modified()
        #self.contour_stripper.Modified()
        self.contour_mapper.Modified()
        #self.contour_filter.Update()
        #self.contour_stripper.Update()
        #self.contour_mapper.Update()

        number_of_contour_lines = self.contour_stripper.GetOutput().GetNumberOfLines()
        self.log.info('There are %s contours lines.' % number_of_contour_lines)


def _get_normalized_data(case):
    """helper method for ``_get_fringe_data``"""
#def _get_normalized_data(result_case):
    #(obj, (i, name)) = result_case
    #case = obj.get_result(i, name)
    if case is None:
        return None

    if len(case.shape) == 1:
        normi = case
    else:
        assert isinstance(case, np.ndarray), case
        normi = norm(case, axis=1)
    return normi

def normalize_forces(forces_array):
    """normalizes the forces"""
    mag = np.linalg.norm(forces_array, axis=1)
    #assert len(forces_array) == len(mag)

    mag_max = mag.max()
    if mag_max > 0.:
        new_forces = np.copy(forces_array / mag_max)
    else:
        new_forces = np.copy(forces_array)
    #mag /= mag_max

    #inonzero = np.where(mag > 0)[0]
    #print('new_forces_max =', new_forces.max())
    #print('new_forces =', new_forces[inonzero])
    #print('mag =', mag[inonzero])
    return new_forces, mag

def show_hide_actor(actor):
    """flips the visibility for an actor"""
    is_visible = actor.GetVisibility()
    actor.SetVisibility(not is_visible)
    actor.Modified()
