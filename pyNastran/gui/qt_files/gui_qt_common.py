# -*- coding: utf-8 -*-
# pylint: disable=C0111
from __future__ import print_function
from six import iteritems
from copy import deepcopy

import numpy as np
from numpy import ndarray, full, issubdtype
from numpy.linalg import norm
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from pyNastran.utils import integer_types
from pyNastran.gui.gui_objects.names_storage import NamesStorage
from pyNastran.gui.testing_methods import GuiAttributes


class GuiCommon(GuiAttributes):
    def __init__(self, inputs):
        GuiAttributes.__init__(self, inputs, res_widget=None)
        self.is_groups = inputs['is_groups']

        #self.groups = set([])
        self._group_elements = {}
        self._group_coords = {}
        self._group_shown = {}
        self._names_storage = NamesStorage()

        self.vtk_version = [int(i) for i in vtk.VTK_VERSION.split('.')[:1]]
        print('vtk_version = %s' % (self.vtk_version))

    #def nCells(self):
        #try:
            #cell_data = self.grid.GetCellData()
            #return cell_data.GetNumberOfCells()
        #except AttributeError:
            #return 0

    #def nPoints(self):
        #try:
            #point_data = self.grid.GetPointData()
            #return point_data.GetNumberOfPoints()
        #except AttributeError:
            #return 0

    def update_axes_length(self, dim_max):
        """
        sets the driving dimension for:
          - picking?
          - coordinate systems
          - label size
        """
        self.dim_max = dim_max
        dim = self.dim * 0.10
        self.on_set_axes_length(dim)

    def on_set_axes_length(self, dim=None):
        """
        scale coordinate system based on model length
        """
        if dim is None:
            dim = self.dim_max * 0.10
        if hasattr(self, 'axes'):
            for cid, axes in iteritems(self.axes):
                axes.SetTotalLength(dim, dim, dim)

    def update_text_actors(self, subcase_id, subtitle, min_value, max_value, label):
        self.text_actors[0].SetInput('Max:  %g' % max_value)  # max
        self.text_actors[1].SetInput('Min:  %g' % min_value)  # min
        self.text_actors[2].SetInput('Subcase: %s Subtitle: %s' % (subcase_id, subtitle))  # info

        if label:
            self.text_actors[3].SetInput('Label: %s' % label)  # info
            self.text_actors[3].VisibilityOn()
        else:
            self.text_actors[3].VisibilityOff()

    def cycle_results(self, result_name=None):
        if self.ncases <= 1:
            self.log.warning('cycle_results(result_name=%r); ncases=%i' % (result_name, self.ncases))
            if self.ncases == 0:
                self.scalarBar.SetVisibility(False)
            return
        result_type = self.cycle_results_explicit(result_name, explicit=False)
        self.log_command('cycle_results(result_name=%r)' % result_type)

    def get_subtitle_label(self, subcase_id):
        try:
            subtitle, label = self.iSubcaseNameMap[subcase_id]
        except KeyError:
            subtitle = 'case=NA'
            label = 'label=NA'
        return subtitle, label

    def cycle_results_explicit(self, result_name=None, explicit=True):
        #if explicit:
            #self.log_command('cycle_results(result_name=%r)' % result_name)
        found_cases = self.increment_cycle(result_name)
        if found_cases:
            result_type = self._set_case(result_name, self.icase, explicit=explicit, cycle=True)
        else:
            result_type = None
        #else:
            #self.log_command(""didn't find case...")
        return result_type

    def get_name_result_data(self, icase):
        key = self.case_keys[icase]
        if isinstance(key, integer_types):
            (obj, (i, name)) = self.result_cases[key]
            #subcase_id = obj.subcase_id
            case = obj.get_result(i, name)
        else:
            assert len(key) == 7, key
            #(subcase_id, j, result_type, vector_size, location, data_format, label2) = key
        return name, case

    def _set_case(self, result_name, icase, explicit=False, cycle=False, skip_click_check=False,
                  min_value=None, max_value=None, is_legend_shown=None):
        if not skip_click_check:
            if not cycle and icase == self.icase:
                # don't click the button twice
                # cycle=True means we're cycling
                # cycle=False skips that check
                return

        try:
            key = self.case_keys[icase]
        except:
            print('icase=%s case_keys=%s' % (icase, str(self.case_keys)))
            raise
        self.icase = icase
        case = self.result_cases[key]
        label2 = ''
        if isinstance(key, integer_types):
            (obj, (i, name)) = self.result_cases[key]
            subcase_id = obj.subcase_id
            case = obj.get_result(i, name)
            result_type = obj.get_title(i, name)
            vector_size = obj.get_vector_size(i, name)
            location = obj.get_location(i, name)
            data_format = obj.get_data_format(i, name)
            scale = obj.get_scale(i, name)
            label2 = obj.get_header(i, name)
            nlabels, labelsize, ncolors, colormap = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
            #default_max, default_min = obj.get_default_min_max(i, name)
            if min_value is None and max_value is None:
                min_value, max_value = obj.get_min_max(i, name)
        else:
            assert len(key) == 7, key
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
            nlabels = None
            labelsize = None
            ncolors = None
            colormap = 'jet'
            #normi = case
            scale = 0.0
            if min_value is None and max_value is None:
                max_value = case.max()
                min_value = case.min()
            if np.isnan(max_value):
                inotnan = not np.isnan(case)
                max_value = case[inotnan].max()
                min_value = case[inotnan].min()
                print('max_value = ', max_value)

        subtitle, label = self.get_subtitle_label(subcase_id)
        if label2:
            label += '; ' + label2
        print("subcase_id=%s result_type=%r subtitle=%r label=%r"
              % (subcase_id, result_type, subtitle, label))

        #================================================

        if isinstance(case, ndarray):
            if len(case.shape) == 1:
                normi = case
            else:
                normi = norm(case, axis=1)
        else:
            msg = 'list-based results have been removed; use numpy.array; name=%s case=%s' % (
                result_name, case)
            raise RuntimeError(msg)

        #if min_value is None and max_value is None:
            #max_value = normi.max()
            #min_value = normi.min()

        #================================================
        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)

        vector_size = 1
        name = (vector_size, subcase_id, result_type, label, min_value, max_value, scale)
        if self._names_storage.has_exact_name(name):
            grid_result = None
        else:
            grid_result = self.set_grid_values(name, normi, vector_size,
                                               min_value, max_value, norm_value)

        vector_size = 3
        name_vector = (vector_size, subcase_id, result_type, label, min_value, max_value, scale)
        if self._names_storage.has_exact_name(name_vector):
            grid_result_vector = None
        else:
            grid_result_vector = self.set_grid_values(name_vector, case, vector_size,
                                                      min_value, max_value, norm_value)

        self.update_text_actors(subcase_id, subtitle,
                                min_value, max_value, label)

        self.final_grid_update(name, grid_result,
                               name_vector, grid_result_vector,
                               key, subtitle, label)

        is_low_to_high = True
        if is_legend_shown is None:
            is_legend_shown = self.scalar_bar.is_shown
        self.update_scalar_bar(result_type, min_value, max_value, norm_value,
                               data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_low_to_high=is_low_to_high,
                               is_horizontal=self.is_horizontal_scalar_bar,
                               is_shown=is_legend_shown)
        self.update_legend(icase,
                           result_type, min_value, max_value, data_format, scale,
                           nlabels, labelsize, ncolors, colormap,
                           is_low_to_high, self.is_horizontal_scalar_bar)
        location = self.get_case_location(key)
        self.res_widget.update_method(location)
        if explicit:
            self.log_command('cycle_results(result_name=%r)' % result_type)
        return result_type

    def set_grid_values(self, name, case, vector_size, min_value, max_value, norm_value,
                        is_low_to_high=True):
        """
        https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
        """
        if self._names_storage.has_exact_name(name):
            return
        #print('name, case =', name, case)

        if not hasattr(case, 'dtype'):
            raise RuntimeError('name=%s case=%s' % (name, case))

        if issubdtype(case.dtype, np.integer):
            data_type = vtk.VTK_INT
            self.grid_mapper.InterpolateScalarsBeforeMappingOn()
        elif issubdtype(case.dtype, np.float):
            data_type = vtk.VTK_FLOAT
            self.grid_mapper.InterpolateScalarsBeforeMappingOff()
        else:
            raise NotImplementedError(case.dtype.type)


        if 0: # nan testing
            if case.dtype.name == 'float32':
                case[50] = np.float32(1) / np.float32(0)
            else:
                case[50] = np.int32(1) / np.int32(0)

        if vector_size == 1:
            if is_low_to_high:
                if norm_value == 0:
                    nvalues = len(case)
                    case2 = full((nvalues), 1.0 - min_value, dtype='float32')
                else:
                    case2 = 1.0 - (case - min_value) / norm_value
            else:
                if norm_value == 0:
                    case2 = full((nvalues), min_value, dtype='float32')
                else:
                    case2 = (case - min_value) / norm_value

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

    def get_result_data_from_icase(self, icase):
        """
        Gets the data stored in the object in a more generic way.

        Parameters
        ----------
        icase : int
            the case ID

        Returns
        -------
           see ``self.get_result_data_from_key(key)``
        """
        key = self.case_keys[icase]
        return self.get_result_data_from_key(key)

    def get_result_data_from_key(self, key):
        """
        You should probably use ``self.get_result_data_from_icase(icase)``
        instead of this.

        Parameters
        ----------
        key : tuple(varies)
            the result key

        Returns
        -------
        obj : varies
            the object that stores the data, if it exists (or None)
        i : int/None
            the object array index
            None : obj is None
        j : int/None
            does ???
            None : obj is not None
        res_name : str
            the scalar bar default title (???)
        result_type : str
            the scalar bar title (???)
        subcase_id : int
            the subcase ID
        vector_size : int
            1 - scalar quantity
            3 - vector quantity
        data_format : str
            Python string formatter (e.g. '%.3f')
        label2 : str
            something to append to the text at the bottom of the viewport
        """
        obj = None
        i = None
        j = None
        if isinstance(key, int):
            (obj, (i, res_name)) = self.result_cases[key]
            subcase_id = obj.subcase_id
            #case = obj.get_result(i, name)
            result_type = obj.get_title(i, res_name)
            vector_size = obj.get_vector_size(i, res_name)
            location = obj.get_location(i, res_name)
            data_format = obj.get_data_format(i, res_name)
            label2 = ''
        else:
            assert len(key) == 7, key
            # j is icase? and is used to...
            # label2 defaults to ''
            res_name = result_type # ???
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key

        return obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2

    def final_grid_update(self, name, grid_result,
                          name_vector, grid_result_vector,
                          key, subtitle, label):
        obj = None
        if isinstance(key, int):
            (obj, (i, res_name)) = self.result_cases[key]
            subcase_id = obj.subcase_id
            #case = obj.get_result(i, name)
            result_type = obj.get_title(i, res_name)
            vector_size = obj.get_vector_size(i, res_name)
            #print('res_name=%s vector_size=%s' % (res_name, vector_size))
            location = obj.get_location(i, res_name)
            #data_format = obj.get_data_format(i, res_name)
        else:
            assert len(key) == 7, key
            # j is icase? and is used to...
            # label2 defaults to ''
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key

        #if vector_size == 3:
            #print('name, grid_result, vector_size=3', name, grid_result)
        self._final_grid_update(name, grid_result, None, None, None,
                                1, subcase_id, result_type, location, subtitle, label,
                                revert_displaced=True)
        if obj is None:
            return
        if vector_size == 3:
            self._final_grid_update(name_vector, grid_result_vector, obj, i, res_name,
                                    vector_size, subcase_id, result_type, location, subtitle, label,
                                    revert_displaced=False)
            #xyz_nominal, vector_data = obj.get_vector_result(i, res_name)
            #self._update_grid(vector_data)

    def _final_grid_update(self, name, grid_result, obj, i, res_name,
                           vector_size, subcase_id, result_type, location, subtitle, label,
                           revert_displaced=True):
        if name is None:
            return
        name_str = self._names_storage.get_name_string(name)
        if not self._names_storage.has_exact_name(name):
            grid_result.SetName(name_str)
            self._names_storage.add(name)

            if self._is_displaced and revert_displaced:
                self._is_displaced = False
                self._update_grid(self._xyz_nominal)

            if location == 'centroid':
                cell_data = self.grid.GetCellData()
                if self._names_storage.has_close_name(name):
                    cell_data.RemoveArray(name_str)
                    self._names_storage.remove(name)

                cell_data.AddArray(grid_result)
                self.log_info("centroidal plotting vector=%s - subcase_id=%s result_type=%s subtitle=%s label=%s"
                              % (vector_size, subcase_id, result_type, subtitle, label))
            elif location == 'node':
                point_data = self.grid.GetPointData()
                if self._names_storage.has_close_name(name):
                    point_data.RemoveArray(name_str)
                    self._names_storage.remove(name)

                if vector_size == 1:
                    self.log_info("node plotting vector=%s - subcase_id=%s result_type=%s subtitle=%s label=%s"
                                  % (vector_size, subcase_id, result_type, subtitle, label))
                    point_data.AddArray(grid_result)
                elif vector_size == 3:
                    #print('vector_size3; get, update')
                    xyz_nominal, vector_data = obj.get_vector_result(i, res_name)

                    #grid_result1 = self.set_grid_values(name, case, 1,
                                                        #min_value, max_value, norm_value)
                    #point_data.AddArray(grid_result1)

                    self._is_displaced = True
                    self._xyz_nominal = xyz_nominal
                    self._update_grid(vector_data)
                    self.log_info("node plotting vector=%s - subcase_id=%s result_type=%s subtitle=%s label=%s"
                                  % (vector_size, subcase_id, result_type, subtitle, label))
                    #point_data.AddVector(grid_result) # old
                    #point_data.AddArray(grid_result)
                else:
                    raise RuntimeError(vector_size)
            else:
                raise RuntimeError(location)

        if location == 'centroid':
            cell_data = self.grid.GetCellData()
            cell_data.SetActiveScalars(name_str)

            point_data = self.grid.GetPointData()
            point_data.SetActiveScalars(None)
        elif location == 'node':
            cell_data = self.grid.GetCellData()
            cell_data.SetActiveScalars(None)

            point_data = self.grid.GetPointData()
            if vector_size == 1:                point_data.SetActiveScalars(name_str)
            elif vector_size == 3:                point_data.SetActiveVectors(name_str)
            else:
                raise RuntimeError(vector_size)
            #print('name_str=%r' % name_str)
        else:
            raise RuntimeError(location)

        self.grid.Modified()
        self.grid_selected.Modified()
        #self.update_all()
        #self.update_all()
        if len(self.groups):
            self.post_group_by_name(self.group_active)
        self.vtk_interactor.Render()

        self.hide_labels(show_msg=False)
        self.show_labels(result_names=[result_type], show_msg=False)

    def _update_grid(self, vector_data):
        nnodes = vector_data.shape[0]
        points = self.grid.GetPoints()
        for j in range(nnodes):
            points.SetPoint(j, *vector_data[j, :])
        self.grid.Modified()
        self.grid_selected.Modified()

    def _get_icase(self, result_name):
        found_case = False
        print('result_cases.keys() =', self.result_cases.keys())
        i = 0
        for icase, cases in sorted(iteritems(self.result_cases)):
            if result_name == icase[1]:
                found_case = True
                icase = i
                break
            i += 1
        assert found_case == True, 'result_name=%r' % result_name
        return icase

    def increment_cycle(self, result_name=False):
        found_case = False
        if result_name is not False and result_name is not None:
            for icase, cases in sorted(iteritems(self.result_cases)):
                if result_name == cases[1]:
                    found_case = True
                    self.icase = icase  # no idea why this works...if it does...

        if not found_case:
            if self.icase is not self.ncases:
                self.icase += 1
            else:
                self.icase = 0
        if self.icase == len(self.case_keys):
            self.icase = 0

        if len(self.case_keys) > 0:
            try:
                key = self.case_keys[self.icase]
            except IndexError:
                found_cases = False
                return found_cases
            except TypeError:
                msg = 'type(case_keys)=%s\n' % type(self.case_keys)
                msg += 'icase=%r\n' % str(self.icase)
                msg += 'case_keys=%r' % str(self.case_keys)
                print(msg)
                raise TypeError(msg)
            msg = 'icase=%r\n' % str(self.icase)
            msg += 'case_keys=%r' % str(self.case_keys)
            #print(msg)

            location = self.get_case_location(key)
            #print("key_increment_cycle = %s" % str(key))
            #if key[2] == 3:  # vector size=3 -> vector, skipping ???
                #self.increment_cycle()
            found_cases = True
        else:
            # key = self.case_keys[self.icase]
            # location = self.get_case_location(key)
            location = 'N/A'
            #result_type = 'centroidal' if location == 'centroid' else 'nodal'
            result_type = '???'
            self.log_error("No Results found.  Many results are not supported in the GUI.\nTry using %s results."
                           % result_type)
            self.scalarBar.SetVisibility(False)
            found_cases = False
        #print("next icase=%s key=%s" % (self.icase, key))
        return found_cases

    def get_result_name(self, key):
        if isinstance(key, int):
            (obj, (i, name)) = self.result_cases[key]
            return name
        else:
            assert len(key) == 7, '%s = (subcase_id, j, result_type, vector_size, location, data_format, label2)' % str(key)
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
        return result_type

    def get_case_location(self, key):
        if isinstance(key, int):
            (obj, (i, name)) = self.result_cases[key]
            return obj.get_location(i, name)
        else:
            assert len(key) == 7, key
            try:
                (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
            except:
                self.log.error(key)
                return
        return location

