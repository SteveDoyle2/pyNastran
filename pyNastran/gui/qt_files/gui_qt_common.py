# -*- coding: utf-8 -*-
from __future__ import print_function
from six import iteritems
import vtk
from numpy import ndarray, asarray, hstack, searchsorted


class GuiCommon(object):
    def __init__(self):
        self.nvalues = 9
        self.groups = set([])
        self._group_elements = {}
        self._group_coords = {}
        self._group_shown = {}
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
        # scale based on model length
        dim_max *= 0.10
        if hasattr(self, 'axes'):
            for cid, axes in iteritems(self.axes):
                axes.SetTotalLength(dim_max, dim_max, dim_max)

    def update_text_actors(self, case, subcaseID, subtitle, min_value, max_value, label):
        self.textActors[0].SetInput('Max:  %g' % max_value)  # max
        self.textActors[1].SetInput('Min:  %g' % min_value)  # min
        self.textActors[2].SetInput('Subcase: %s Subtitle: %s' % (subcaseID, subtitle))  # info

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
        self.cycleResults_explicit(result_name, explicit=False)

    def cycleResults_explicit(self, result_name=None, explicit=True):
        #if explicit:
            #self.log_command('cycleResults(result_name=%r)' % result_name)
        print("is_nodal=%s is_centroidal=%s" % (self.is_nodal, self.is_centroidal))

        found_cases = self.incrementCycle(result_name)
        if found_cases:
            self._set_case(result_name, self.iCase, explicit=explicit, cycle=True)
        #else:
            #self.log_command(""didn't find case...")

    def _set_case(self, result_name, iCase, explicit=False, cycle=False):
        if not cycle and iCase == self.iCase:
            #print('double...')
            return  # don't click the button twice

        try:
            key = self.caseKeys[iCase]
        except:
            print('icase=%s caseKeys=%s' % (iCase, str(self.caseKeys)))
            raise
        self.iCase = iCase
        case = self.resultCases[key]
        label2 = ''
        if len(key) == 5:
            (subcase_id, result_type, vector_size, location, data_format) = key
        elif len(key) == 6:
            (subcase_id, j, result_type, vector_size, location, data_format) = key
        else:
            (subcase_id, j, result_type, vector_size, location, data_format, label2) = key

        try:
            caseName = self.iSubcaseNameMap[subcase_id]
        except KeyError:
            caseName = ('case=NA', 'label=NA')
        (subtitle, label) = caseName
        label += label2
        print("subcaseID=%s resultType=%r subtitle=%r label=%r" % (subcase_id, result_type, subtitle, label))

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

        norm_value, nValueSet = self.set_grid_values(grid_result, case, vector_size, min_value, max_value)
        self.update_text_actors(case, subcase_id, subtitle, min_value, max_value, label)
        self.UpdateScalarBar(result_type, min_value, max_value, norm_value, data_format, is_blue_to_red=True)
        #self.scalarBar.SetNumberOfLabels(nValueSet)
        #self.scalarBar.SetMaximumNumberOfColors(nValueSet)
        #prop = self.scalarBar.GetLabelTextProperty()
        #fontSize = prop.GetFontSize()
        #print("fontSize = %s" % fontSize)
        #prop.SetFontSize(40)

        # TODO results can only go from centroid->node and not back to
        # centroid
        #print(dir(self.grid))
        #self.grid.Reset()
        self.final_grid_update(grid_result, key, subtitle, label)
        if explicit:
            self.log_command('cycleResults(result_name=%r)' % result_type)

    def set_grid_values(self, gridResult, case, vectorSize, min_value, max_value, is_blue_to_red=True):
        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)
        #print('max_value=%s min_value=%r norm_value=%r' % (max_value, min_value, norm_value))
        #print("case = ", case)
        #if norm_value == 0.: # avoids division by 0.
        #    norm_value = 1.

        valueSet = set()
        if vectorSize == 1:
            if is_blue_to_red:
                if norm_value == 0:
                    for i, value in enumerate(case):
                        gridResult.InsertNextValue(1 - min_value)
                else:
                    for i, value in enumerate(case):
                        gridResult.InsertNextValue(1.0 - (value - min_value) / norm_value)
            else:
                if norm_value == 0:
                    for i, value in enumerate(case):
                        gridResult.InsertNextValue(min_value)
                else:
                    for i, value in enumerate(case):
                        gridResult.InsertNextValue((value - min_value) / norm_value)
        else:  # vectorSize=3
            for value in case:
                gridResult.InsertNextTuple3(value)  # x,y,z

        nValueSet = len(valueSet)
        return norm_value, nValueSet

    def final_grid_update(self, gridResult, key, subtitle, label):
        if len(key) == 5:
            (subcaseID, resultType, vectorSize, location, data_format) = key
        elif len(key) == 6:
            (subcaseID, j, resultType, vectorSize, location, data_format) = key
        else:
            (subcaseID, j, resultType, vectorSize, location, data_format, label2) = key
        npoints = self.nPoints()
        ncells = self.nCells()

        if location == 'centroid':
        #if location == 'centroid' and self.is_centroidal:
            if npoints:
                point_data = self.grid.GetPointData()
                point_data.Reset()
            self.grid.GetCellData().SetScalars(gridResult)
            self.log_info("centroidal plotting vector=%s - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
        elif location == 'node':
        #elif location == 'node' and self.is_nodal:
            if ncells:
                cell_data = self.grid.GetCellData()
                #print(dir(cell_data))
                cell_data.Reset()
            if vectorSize == 1:
                self.log_info("node plotting vector=%s - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                self.grid.GetPointData().SetScalars(gridResult)
            elif vectorSize == 3:
                self.log_info("node plotting vector=%s - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                self.grid.GetPointData().SetScalars(gridResult)
            else:
                #print("***node skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID, resultType, subtitle, label))
                raise RuntimeError(vectorSize)
        else:
            raise RuntimeError(location)
            self.log_info("***D%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (location, subcaseID, resultType, subtitle, label))
            self.scalarBar.SetVisibility(False)
        self.grid.Modified()
        self.vtk_interactor.Render()

    def _get_icase(self, result_name):
        found_case = False
        print(self.resultCases.keys())
        i = 0
        for icase, cases in sorted(iteritems(self.resultCases)):
            #print(cases[1])
            if result_name == icase[1]:
                found_case = True
                iCase = i
                break
            i += 1
        assert found_case == True, 'result_name=%r' % result_name
        #print('***icase = %s' % iCase)
        return iCase

    def incrementCycle(self, result_name=False):
        found_case = False
        if result_name is not False and result_name is not None:
            for icase, cases in sorted(iteritems(self.resultCases)):
                if result_name == cases[1]:
                    found_case = True
                    self.iCase = icase  # no idea why this works...if it does...

        if not found_case:
            #print('iCase=%s nCases=%s' % (self.iCase, self.nCases))
            if self.iCase is not self.nCases:
                self.iCase += 1
            else:
                self.iCase = 0
        if self.iCase == len(self.caseKeys):
            self.iCase = 0

        if len(self.caseKeys) > 0:
            #print('caseKeys =', self.caseKeys)
            try:
                key = self.caseKeys[self.iCase]
            except IndexError:
                foundCases = False
                return foundCases

            print("key = %s" % str(key))
            #if key[2] == 3:  # vector size=3 -> vector, skipping ???
                #self.incrementCycle()
            foundCases = True
        else:
            result_type = 'centroidal' if self.is_centroidal else 'nodal'
            self.log_error("No Results found.  Many results are not supported in the GUI.\nTry using %s results." % result_type)
            self.scalarBar.SetVisibility(False)
            foundCases = False
        #print("next iCase=%s key=%s" % (self.iCase, key))
        return foundCases

    def build_grid_result(self, vectorSize, location):
        #gridResult.Reset()
        gridResult = vtk.vtkFloatArray()
        #gridResult.Reset()
        #gridResult.Modified()
        self.rend.Modified()
        #emptyResult = vtk.vtkFloatArray()

        gridResult.SetNumberOfComponents(vectorSize)
        if location == 'centroid':
        #if location == 'centroid' and self.is_centroidal:
            #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
            gridResult.Allocate(self.nElements, 1000)
        #elif location == 'node' and self.is_nodal:
        elif location == 'node':
            #allocationSize = vectorSize*self.nNodes # (where location='node'-> self.nNodes)
            gridResult.Allocate(self.nNodes * vectorSize, 1000)
            #gridResult.SetNumberOfComponents(vectorSize)
        else:
            raise RuntimeError(location)
            #print("***%s skipping" % location)
        return gridResult

    def UpdateScalarBar(self, Title, min_value, max_value, norm_value, data_format, is_blue_to_red=True):
        """
        :param Title the:   scalar bar title
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
        self.scalarBar.SetTitle(Title)

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

        #if (Title in ['ElementID', 'Eids', 'Region'] and norm_value < 11):
            #nvalues = int(max_value - min_value) + 1
            #print("need to adjust axes...max_value=%s" % max_value)

        if self.nvalues is not None:
            if not self._is_int_result(data_format):  # don't change nvalues for int results
                nvalues = self.nvalues

        self.scalarBar.SetNumberOfLabels(nvalues)
        self.scalarBar.SetMaximumNumberOfColors(nvalues)
        self.scalarBar.Modified()


#========================================================================================
# Groups - not done

    def clear_groups(self):
        all_groups = self.groups # set
        self.remove_groups(all_groups)

    def remove_groups(self, groups):
        assert isinstance(groups, list), type(groups)
        assert len(groups) >= 0, 'groups is empty'
        assert len(all_groups) >= 0, 'all_groups is empty'

        all_groups = self.groups # set
        for group in all_groups:
            if group in groups:
                self.remove_group(group)

    def remove_group(self, group):
        if group not in all_groups:
            raise RuntimeError('group=%r not found' % group)

    def show_group(self, name):
        self._group_shown[name] = True

    def hide_group(self, name):
        self._group_shown[name] = False

    def post_groups(self, groups):
        assert isinstance(groups, list), type(groups)
        assert len(groups) >= 0, 'groups is empty'

        all_groups = self.groups # set
        assert len(all_groups) >= 0, 'all_groups is empty'

        for group in all_groups:
            if group in groups:
                self.show_group(group)
            else:
                self.show_group(group)

    def _check_add(self, Format, name, element_id=None, property_id=None, coord_id=None):
        if element_id is None and property_id is None:
            raise RuntimeError('either element_id or property_id must be set')
        if isinstance(element_id, int):
            element_id = [element_id]
        if isinstance(property_id, int):
            property_id = [property_id]

        if Format == 'nastran':
            if property_id:
                element_id = self.model.get_element_id_by_property_id(property_id)

        elif Format == 'cart3d':
            if property_id:
                element_id = self.model.get_gelement_id_by_region_id(property_id)
        elif Format == 'panair':
            if element_id is None:
                raise RuntimeError('element_id must be set for panair')
        else:
            msg = "Format=%r is not supported; use 'nastran', 'cart3d', 'panair'" % Format
            raise NotImplementedError(msg)

        if coord_id is not None and Format != 'nastran':
            raise RuntimeError('coord_id must be None for format=%r' % Format)

        element_id = asarray(element_id)
        return element_id

    def _add_coord_id(self, name, coord_id):
        if coord_id is None:
            coord_id = set([0])
        elif isinstance(coord_id, int):
            coord_id = set([coord_id])
        else:
            for cid in coord_id:
                assert isinstance(cid, int), type(cid)
        if name in self._group_coords:
            self._group_coords[name].union(set(coord_id))
        else:
            self._group_coords[name] = set(coord_id)

    def _create_grid_mapper(self, name):
        self.grid = vtk.vtkUnstructuredGrid()
        self.aQuadMapper = vtk.vtkDataSetMapper()
        self.aQuadMapper.SetInput(self.grid)

        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(self.aQuadMapper)
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(geometryActor)

class Groups(object):
    def __init__(self):
        self.nNodes = None
        self.model = None
        self.grid = None

    def add_to_group(self, Format, name, element_id=None, property_id=None, coord_id=None):
        assert name in self._group_elements
        element_id = self._check_add(Format, name,
                                     element_id=element_id, property_id=property_id,
                                     coord_id=coord_id)
        self._group_elements[name] = hstack([self._group_elements[name],
                                             element_id])
        self._add_coord_id(name, coord_id)

    def create_group(self, Format, name,
                     element_id=None, property_id=None, coord_id=None, show=True):
        element_id = self._check_add(Format, name,
                                     element_id=element_id, property_id=property_id,
                                     coord_id=coord_id)

        self.groups.add(name)
        self._group_elements[name] = element_id
        self._add_coord_id(name, coord_id)
        self._group_shown[name] = show
        if Format == 'nastran':
            self._create_nastran_group(name, self.model, element_id)
        elif Format == 'cart3d':
            self._create_cart3d_group(name, self.model, element_id)
        #elif Format == 'panair':
            #self._create_panair_group(name, self.model, element_id)

    def _create_nastran_group(self, name, model, element_id):
        pass

    def _create_cart3d_group(self, name, model, element_id):
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        nelements = len(element_id)
        self.grid.Allocate(nelements, 1000)

        nodes = self.model.get_nodes_associated_with_elements(element_id)
        nodes.sort()

        nid = 0
        all_nodes = self.nodes
        for i in all_nodes:
            #if nid in nodes:
                points.InsertPoint(nid, all_nodes[i, :])
                nid += 1

        from vtk import vtkTriangle
        for eid in element_id:
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem_nodes = searchsorted(nodes, node_ids)
            elem.GetPointIds().SetId(0, elem_nodes[0])
            elem.GetPointIds().SetId(1, elem_nodes[1])
            elem.GetPointIds().SetId(2, elem_nodes[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())

        self.grid[name].SetPoints(points)
        self.grid[name].Modified()
        self.grid[name].Update()

