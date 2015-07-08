# -*- coding: utf-8 -*-
# pylint: disable=C0111
from __future__ import print_function
from six import iteritems
from copy import deepcopy

from numpy import ndarray, asarray, hstack, searchsorted, ones
import vtk
from vtk.util.numpy_support import numpy_to_vtk


class GuiCommon(object):
    def __init__(self):
        self.nvalues = 9
        self.groups = set([])
        self._group_elements = {}
        self._group_coords = {}
        self._group_shown = {}
        self.dim_max = 1.0
        self.vtk_version = [int(i) for i in vtk.VTK_VERSION.split('.')[:1]]
        print('vtk_version = %s' % (self.vtk_version))

    def nCells(self):
        try:
            cell_data = self.grid.GetCellData()
            return cell_data.GetNumberOfCells()
        except AttributeError:
            return 0

    def nPoints(self):
        try:
            point_data = self.grid.GetPointData()
            return point_data.GetNumberOfPoints()
        except AttributeError:
            return 0

    def _is_int_result(self, data_format):
        if 'i' in data_format:
            return True
        return False

    def update_axes_length(self, dim_max):
        """
        scale coordinate system based on model length
        """
        self.dim_max = dim_max
        dim_max *= 0.10
        if hasattr(self, 'axes'):
            for cid, axes in iteritems(self.axes):
                axes.SetTotalLength(dim_max, dim_max, dim_max)

    def update_text_actors(self, case, subcase_id, subtitle, min_value, max_value, label):
        self.textActors[0].SetInput('Max:  %g' % max_value)  # max
        self.textActors[1].SetInput('Min:  %g' % min_value)  # min
        self.textActors[2].SetInput('Subcase: %s Subtitle: %s' % (subcase_id, subtitle))  # info

        if label:
            self.textActors[3].SetInput('Label: %s' % label)  # info
            self.textActors[3].VisibilityOn()
        else:
            self.textActors[3].VisibilityOff()

    def cycleResults(self, result_name=None):
        if self.nCases <= 1:
            self.log.warning('cycleResults(result_name=%r); nCases=%i' % (result_name, self.nCases))
            if self.nCases == 0:
                self.scalarBar.SetVisibility(False)
            return
        result_type = self.cycleResults_explicit(result_name, explicit=False)
        self.log_command('cycleResults(result_name=%r)' % result_type)

    def cycleResults_explicit(self, result_name=None, explicit=True):
        #if explicit:
            #self.log_command('cycleResults(result_name=%r)' % result_name)
        print("is_nodal=%s is_centroidal=%s" % (self.is_nodal, self.is_centroidal))

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
        if len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
        elif len(key) == 6:
            (subcase_id, j, result_type, vector_size, location, data_format) = key
        else:
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key

        try:
            case_name = self.iSubcaseNameMap[subcase_id]
        except KeyError:
            case_name = ('case=NA', 'label=NA')
        (subtitle, label) = case_name
        label += label2
        print("subcase_id=%s result_type=%r subtitle=%r label=%r"
              % (subcase_id, result_type, subtitle, label))

        #================================================
        grid_result = self.build_grid_result(vector_size, location)
        #================================================
        if isinstance(case, ndarray):
            max_value = case.max()
            min_value = case.min()
        else:
            raise RuntimeError('list-based results have been disabled; use numpy.array')
            print('resultType=%r should really use numpy arrays...' % result_type)
            max_value = case[0]
            min_value = case[0]
            for value in case:
                max_value = max(value, max_value)
                min_value = min(value, min_value)

        norm_value, nvalues_set, grid_result = self.set_grid_values(grid_result, case, vector_size, min_value, max_value)
        self.update_text_actors(case, subcase_id, subtitle, min_value, max_value, label)
        self.UpdateScalarBar(result_type, min_value, max_value, norm_value, data_format, is_blue_to_red=True)

        # TODO: results can only go from centroid->node and not back to centroid
        self.final_grid_update(grid_result, key, subtitle, label)
        if explicit:
            self.log_command('cycleResults(result_name=%r)' % result_type)
        return result_type

    def set_grid_values(self, grid_result, case, vector_size,
                        min_value, max_value,
                        is_blue_to_red=True):
        """
        https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
        """
        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)
        #print('max_value=%s min_value=%r norm_value=%r' % (max_value, min_value, norm_value))
        #print("case = ", case)
        #if norm_value == 0.: # avoids division by 0.
        #    norm_value = 1.

        #warp_vector = vtk.vtkWarpVector()
        #warp_vector.setInput(grid_result.GetOuput())

        value_set = set()
        if vector_size == 1:
            if is_blue_to_red:
                if norm_value == 0:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue(1 - min_value)
                    nvalues = len(case)
                    case2 = 1 - ones(nvalues) * min_value
                else:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue(1.0 - (value - min_value) / norm_value)
                    case2 = 1.0 - (case - min_value) / norm_value
            else:
                if norm_value == 0:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue(min_value)
                    case2 = ones(nvalues) * min_value
                else:
                    #for i, value in enumerate(case):
                        #grid_result.InsertNextValue((value - min_value) / norm_value)
                    case2 = (case - min_value) / norm_value
            grid_result = numpy_to_vtk(
                num_array=case2,
                deep=False,
                array_type=vtk.VTK_FLOAT
            )
        else:
            # vector_size=3
            #for value in case:
                #grid_result.InsertNextTuple3(*value)  # x, y, z
            if case.flags.contiguous:
                case2 = case
                deep = False
            else:
                case2 = deepcopy(case)
                deep = True
            grid_result = numpy_to_vtk(
                num_array=case2,
                deep=deep,
                array_type=vtk.VTK_FLOAT
            )

        nvalues_set = len(value_set)
        return norm_value, nvalues_set, grid_result

    def final_grid_update(self, grid_result, key, subtitle, label):
        if len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
        elif len(key) == 6:
            (subcase_id, j, result_type, vector_size, location, data_format) = key
        else:
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key
        npoints = self.nPoints()
        ncells = self.nCells()

        if location == 'centroid':
        #if location == 'centroid' and self.is_centroidal:
            if npoints:
                point_data = self.grid.GetPointData()
                point_data.Reset()
            self.grid.GetCellData().SetScalars(grid_result)
            self.log_info("centroidal plotting vector=%s - subcase_id=%s result_type=%s subtitle=%s label=%s"
                          % (vector_size, subcase_id, result_type, subtitle, label))
        elif location == 'node':
        #elif location == 'node' and self.is_nodal:
            if ncells:
                cell_data = self.grid.GetCellData()
                #print(dir(cell_data))
                cell_data.Reset()
            if vector_size == 1:
                self.log_info("node plotting vector=%s - subcase_id=%s result_type=%s subtitle=%s label=%s"
                              % (vector_size, subcase_id, result_type, subtitle, label))
                self.grid.GetPointData().SetScalars(grid_result)
            elif vector_size == 3:
                self.log_info("node plotting vector=%s - subcase_id=%s result_type=%s subtitle=%s label=%s"
                              % (vector_size, subcase_id, result_type, subtitle, label))
                self.grid.GetPointData().SetScalars(grid_result)
            else:
                #print("***node skipping - subcase_id=%s result_type=%s subtitle=%s label=%s"
                      #% (subcase_id, result_type, subtitle, label))
                raise RuntimeError(vector_size)
        else:
            raise RuntimeError(location)
            #self.log_info("***D%s skipping - subcase_id=%s result_type=%s subtitle=%s label=%s"
                          #% (location, subcase_id, result_type, subtitle, label))
            #self.scalarBar.SetVisibility(False)
        self.grid.Modified()
        self.vtk_interactor.Render()

        self.hide_labels(show_msg=False)
        self.show_labels(result_names=[result_type], show_msg=False)

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

            print("key = %s" % str(key))
            #if key[2] == 3:  # vector size=3 -> vector, skipping ???
                #self.incrementCycle()
            found_cases = True
        else:
            result_type = 'centroidal' if self.is_centroidal else 'nodal'
            self.log_error("No Results found.  Many results are not supported in the GUI.\nTry using %s results."
                           % result_type)
            self.scalarBar.SetVisibility(False)
            found_cases = False
        #print("next iCase=%s key=%s" % (self.iCase, key))
        return found_cases

    def build_grid_result(self, vector_size, location):
        grid_result = vtk.vtkFloatArray()
        self.rend.Modified()

        grid_result.SetNumberOfComponents(vector_size)
        if location == 'centroid':
        #if location == 'centroid' and self.is_centroidal:
            #allocationSize = vector_size*location (where location='centroid'-> self.nElements)
            grid_result.Allocate(self.nElements, 1000)
        #elif location == 'node' and self.is_nodal:
        elif location == 'node':
            #allocationSize = vector_size*self.nNodes # (where location='node'-> self.nNodes)
            grid_result.Allocate(self.nNodes * vector_size, 1000)
            #grid_result.SetNumberOfComponents(vector_size)
        else:
            raise RuntimeError(location)
            #print("***%s skipping" % location)
        return grid_result

    def UpdateScalarBar(self, title, min_value, max_value, norm_value, data_format,
                        is_blue_to_red=True):
        """
        :param title:       the scalar bar title
        :param min_value:   the blue value
        :param max_value:   the red value
        :param data_format: '%g','%f','%i', etc.
        """
        print("UpdateScalarBar min=%s max=%s norm=%s" % (min_value, max_value, norm_value))
        self.colorFunction.RemoveAllPoints()

        if is_blue_to_red:
            self.colorFunction.AddRGBPoint(min_value, 0.0, 0.0, 1.0)  # blue
            self.colorFunction.AddRGBPoint(max_value, 1.0, 0.0, 0.0)  # red
        else:
            self.colorFunction.AddRGBPoint(min_value, 1.0, 0.0, 0.0)  # red
            self.colorFunction.AddRGBPoint(max_value, 0.0, 0.0, 1.0)  # blue

        #self.scalarBar.SetLookupTable(self.colorFunction)
        self.scalarBar.SetTitle(title)

        nvalues = 11
        data_format_display = data_format
        if data_format == '%i':
            data_format_display = '%.0f'
            nvalues = int(max_value - min_value) + 1
            if nvalues < 7:
                nvalues = 7
            elif nvalues > 30:
                nvalues = 11
        self.scalarBar.SetLabelFormat(data_format_display)

        #if title in ['ElementID', 'Eids', 'Region'] and norm_value < 11:
            #nvalues = int(max_value - min_value) + 1
            #print("need to adjust axes...max_value=%s" % max_value)

        if self.nvalues is not None:
            if not self._is_int_result(data_format):
                # don't change nvalues for int results
                nvalues = self.nvalues

        self.scalarBar.SetNumberOfLabels(nvalues)
        self.scalarBar.SetMaximumNumberOfColors(nvalues)
        self.scalarBar.Modified()


