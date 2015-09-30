# -*- coding: utf-8 -*-
# pylint: disable=C0111
from __future__ import print_function
from six import iteritems
from copy import deepcopy

import numpy
from numpy import ndarray, asarray, hstack, searchsorted, ones, full, int32, float32, issubdtype
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from pyNastran.gui.names_storage import NamesStorage


class GuiCommon(object):
    def __init__(self):
        self._is_displaced = False
        self._xyz_nominal = None

        self.nvalues = 9
        self.groups = set([])
        self._group_elements = {}
        self._group_coords = {}
        self._group_shown = {}
        self._names_storage = NamesStorage()

        self.dim_max = 1.0
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
        scale coordinate system based on model length
        """
        self.dim_max = dim_max
        dim_max *= 0.10
        if hasattr(self, 'axes'):
            for cid, axes in iteritems(self.axes):
                axes.SetTotalLength(dim_max, dim_max, dim_max)

    def update_text_actors(self, subcase_id, subtitle, min_value, max_value, label):
        self.text_actors[0].SetInput('Max:  %g' % max_value)  # max
        self.text_actors[1].SetInput('Min:  %g' % min_value)  # min
        self.text_actors[2].SetInput('Subcase: %s Subtitle: %s' % (subcase_id, subtitle))  # info

        if label:
            self.text_actors[3].SetInput('Label: %s' % label)  # info
            self.text_actors[3].VisibilityOn()
        else:
            self.text_actors[3].VisibilityOff()

    def cycleResults(self, result_name=None):
        if self.nCases <= 1:
            self.log.warning('cycleResults(result_name=%r); nCases=%i' % (result_name, self.nCases))
            if self.nCases == 0:
                self.scalarBar.SetVisibility(False)
            return
        result_type = self.cycleResults_explicit(result_name, explicit=False)
        self.log_command('cycleResults(result_name=%r)' % result_type)

    def get_subtitle_label(self, subcase_id):
        try:
            subtitle, label = self.iSubcaseNameMap[subcase_id]
        except KeyError:
            subtitle = 'case=NA'
            label = 'label=NA'
        return subtitle, label

    def cycleResults_explicit(self, result_name=None, explicit=True):
        #if explicit:
            #self.log_command('cycleResults(result_name=%r)' % result_name)
        found_cases = self.incrementCycle(result_name)
        if found_cases:
            result_type = self._set_case(result_name, self.iCase, explicit=explicit, cycle=True)
        else:
            result_type = None
        #else:
            #self.log_command(""didn't find case...")
        return result_type

    def _set_case(self, result_name, icase, explicit=False, cycle=False):
        if not cycle and icase == self.iCase:
            # don't click the button twice
            return

        try:
            key = self.caseKeys[icase]
        except:
            print('icase=%s caseKeys=%s' % (icase, str(self.caseKeys)))
            raise
        self.iCase = icase
        case = self.resultCases[key]
        label2 = ''
        if isinstance(key, (int, int32)):
            (obj, (i, name)) = self.resultCases[key]
            subcase_id = obj.subcase_id
            case = obj.get_result(i, name)
            result_type = obj.get_title(i, name)
            vector_size = obj.get_vector_size(i, name)
            location = obj.get_location(i, name)
            data_format = obj.get_data_format(i, name)
            scale = obj.get_scale(i, name)
        elif len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
            scale = 0.0
        elif len(key) == 6:
            (subcase_id, j, result_type, vector_size, location, data_format) = key
            scale = 0.0
        else:
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
            scale = 0.0

        subtitle, label = self.get_subtitle_label(subcase_id)
        label += label2
        print("subcase_id=%s result_type=%r subtitle=%r label=%r"
              % (subcase_id, result_type, subtitle, label))

        #================================================
        if isinstance(case, ndarray):
            max_value = case.max()
            min_value = case.min()
        else:
            raise RuntimeError('list-based results have been removed; use numpy.array')

        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)

        vector_size = 1
        name = (vector_size, subcase_id, result_type, label, min_value, max_value, scale)
        if self._names_storage.has_exact_name(name):
            grid_result = None
        else:
            grid_result = self.set_grid_values(name, case, vector_size,
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

        # TODO: results can only go from centroid->node and not back to centroid
        self.final_grid_update(name, grid_result,
                               name_vector, grid_result_vector,
                               key, subtitle, label)

        self.update_scalar_bar(result_type, min_value, max_value, norm_value,
                               data_format, is_blue_to_red=True, is_horizontal=self.is_horizontal_scalar_bar)

        location = self.get_case_location(key)
        self.res_widget.update_method(location)
        if explicit:
            self.log_command('cycleResults(result_name=%r)' % result_type)
        return result_type

    def set_grid_values(self, name, case, vector_size, min_value, max_value, norm_value,
                        is_blue_to_red=True):
        """
        https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
        """
        if self._names_storage.has_exact_name(name):
            return

        if issubdtype(case.dtype, numpy.integer):
            data_type = vtk.VTK_INT
            self.aQuadMapper.InterpolateScalarsBeforeMappingOn()
        elif issubdtype(case.dtype, numpy.float):
            data_type = vtk.VTK_FLOAT
            self.aQuadMapper.InterpolateScalarsBeforeMappingOff()
        else:
            raise NotImplementedError(case.dtype.type)


        if 0: # nan testing
            from numpy import float32, int32
            if case.dtype.name == 'float32':
                case[50] = float32(1) / float32(0)
            else:
                case[50] = int32(1) / int32(0)

        if vector_size == 1:
            if is_blue_to_red:
                if norm_value == 0:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue(1 - min_value)
                    nvalues = len(case)
                    case2 = full((nvalues), 1.0 - min_value, dtype='float32')
                    #case2 = 1 - ones(nvalues) * min_value
                else:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue(1.0 - (value - min_value) / norm_value)
                    case2 = 1.0 - (case - min_value) / norm_value
            else:
                if norm_value == 0:
                    # how do you even get a constant nodal result on a surface?
                    # nodal normals on a constant surface, but that's a bit of an edge case
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue(min_value)
                    case2 = full((nvalues), min_value, dtype='float32')
                    #case2 = case
                else:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue((value - min_value) / norm_value)
                    case2 = (case - min_value) / norm_value

            #scalar_range = self.grid.GetScalarRange()
            #print(scalar_range)
            #self.aQuadMapper.SetScalarRange(scalar_range)

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
            #for value in case:
                #grid_result.InsertNextTuple3(*value)  # x, y, z
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

    def final_grid_update(self, name, grid_result,
                          name_vector, grid_result_vector,
                          key, subtitle, label):
        if isinstance(key, int):
            (obj, (i, res_name)) = self.resultCases[key]
            subcase_id = obj.subcase_id
            #case = obj.get_result(i, name)
            result_type = obj.get_title(i, res_name)
            vector_size = obj.get_vector_size(i, res_name)
            #print('res_name=%s vector_size=%s' % (res_name, vector_size))
            location = obj.get_location(i, res_name)
            data_format = obj.get_data_format(i, res_name)
        elif len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
        elif len(key) == 6:
            (subcase_id, j, result_type, vector_size, location, data_format) = key
        else:
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key

        self._final_grid_update(name, grid_result, None, None, None,
                                1, subcase_id, result_type, location, subtitle, label)
        if vector_size == 3:
            self._final_grid_update(name_vector, grid_result_vector, obj, i, res_name,
                                    vector_size, subcase_id, result_type, location, subtitle, label)

    def _final_grid_update(self, name, grid_result, obj, i, res_name,
                           vector_size, subcase_id, result_type, location, subtitle, label):
        name_str = self._names_storage.get_name_string(name)
        if not self._names_storage.has_exact_name(name):
            grid_result.SetName(name_str)
            self._names_storage.add(name)

            if self._is_displaced:
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
                    xyz_nominal, vector_data, norm = obj.get_vector_result(i, res_name)

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
            if vector_size == 1:
                point_data.SetActiveScalars(name_str)
            elif vector_size == 3:
                point_data.SetActiveVectors(name_str)
            else:
                raise RuntimeError(vector_size)
        else:
            raise RuntimeError(location)

        self.grid.Modified()
        self.vtk_interactor.Render()

        self.hide_labels(show_msg=False)
        self.show_labels(result_names=[result_type], show_msg=False)

    def _update_grid(self, vector_data):
        nnodes = vector_data.shape[0]
        points = self.grid.GetPoints()
        for j in range(nnodes):
            points.SetPoint(j, *vector_data[j,:])
        self.grid.Modified()

    def _get_icase(self, result_name):
        found_case = False
        print(self.resultCases.keys())
        i = 0
        for icase, cases in sorted(iteritems(self.resultCases)):
            if result_name == icase[1]:
                found_case = True
                iCase = i
                break
            i += 1
        assert found_case == True, 'result_name=%r' % result_name
        return iCase

    def incrementCycle(self, result_name=False):
        found_case = False
        if result_name is not False and result_name is not None:
            for icase, cases in sorted(iteritems(self.resultCases)):
                if result_name == cases[1]:
                    found_case = True
                    self.iCase = icase  # no idea why this works...if it does...

        if not found_case:
            if self.iCase is not self.nCases:
                self.iCase += 1
            else:
                self.iCase = 0
        if self.iCase == len(self.caseKeys):
            self.iCase = 0

        if len(self.caseKeys) > 0:
            try:
                key = self.caseKeys[self.iCase]
            except IndexError:
                found_cases = False
                return found_cases

            location = self.get_case_location(key)
            print("key = %s" % str(key))
            #if key[2] == 3:  # vector size=3 -> vector, skipping ???
                #self.incrementCycle()
            found_cases = True
        else:
            # key = self.caseKeys[self.iCase]
            # location = self.get_case_location(key)
            location = 'N/A'
            #result_type = 'centroidal' if location == 'centroid' else 'nodal'
            result_type = '???'
            self.log_error("No Results found.  Many results are not supported in the GUI.\nTry using %s results."
                           % result_type)
            self.scalarBar.SetVisibility(False)
            found_cases = False
        #print("next iCase=%s key=%s" % (self.iCase, key))
        return found_cases

    def get_result_name(self, key):
        if isinstance(key, int):
            (obj, (i, name)) = self.resultCases[key]
            return name
        elif len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
        elif len(key) == 6:
            (subcase_id, i, result_type, vector_size, location, data_format) = key
        else:
            (subcase_id, i, result_type, vector_size, location, data_format, label2) = key
        return result_type

    def get_case_location(self, key):
        if isinstance(key, int):
            (obj, (i, name)) = self.resultCases[key]
            return obj.get_location(i, name)
        elif len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
        elif len(key) == 6:
            (subcase_id, i, result_type, vector_size, location, data_format) = key
        else:
            (subcase_id, i, result_type, vector_size, location, data_format, label2) = key
        return location

