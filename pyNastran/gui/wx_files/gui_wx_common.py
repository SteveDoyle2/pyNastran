# -*- coding: utf-8 -*-
import vtk
from numpy import ndarray

class GuiCommon(object):
    def __init__(self):
        self.nvalues = 9

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
            self.axes.SetTotalLength(dim_max, dim_max, dim_max)

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
        self.cycleResults_explicit(result_name)

    def cycleResults_explicit(self, result_name=None):
        self.log_command('cycleResults(result_name=%r)' % result_name)
        #print('cycling...')
        print("is_nodal=%s is_centroidal=%s" % (self.is_nodal, self.is_centroidal))

        foundCases = self.incrementCycle(result_name)
        if foundCases:
            try:
                key = self.caseKeys[self.iCase]
            except:
                print('icase=%s caseKeys=%s' % (self.iCase, str(self.caseKeys)))
                raise
            case = self.resultCases[key]
            (subcaseID, resultType, vectorSize, location, data_format) = key

            try:
                caseName = self.iSubcaseNameMap[subcaseID]
            except KeyError:
                caseName = ('case=NA', 'label=NA')
            (subtitle, label) = caseName

            print("subcaseID=%s resultType=%r subtitle=%r label=%r" % (subcaseID, resultType, subtitle, label))

            #================================================
            gridResult = self.build_grid_result(vectorSize, location)
            #================================================
            if isinstance(case, ndarray):
                max_value = case.max()
                min_value = case.min()
            else:
                print('resultType=%r should really use numpy arrays...' % resultType)
                max_value = case[0]
                min_value = case[0]
                for value in case:
                    max_value = max(value, max_value)
                    min_value = min(value, min_value)

            norm_value, nValueSet = self.set_grid_values(gridResult, case, vectorSize, min_value, max_value)
            self.update_text_actors(case, subcaseID, subtitle, min_value, max_value, label)
            self.UpdateScalarBar(resultType, min_value, max_value, norm_value, data_format, is_blue_to_red=True)
            #self.scalarBar.SetNumberOfLabels(nValueSet)
            #self.scalarBar.SetMaximumNumberOfColors(nValueSet)
            #prop = self.scalarBar.GetLabelTextProperty()
            #fontSize = prop.GetFontSize()
            #print("fontSize = %s" % fontSize)
            #prop.SetFontSize(40)

            # TODO results can only go from centroid->node and not back to
            ## centroid
            #print(dir(self.grid))
            #self.grid.Reset()
            self.final_grid_update(gridResult, key, subtitle, label)

    def set_grid_values(self, gridResult, case, vectorSize, min_value, max_value, is_blue_to_red=True):
        # flips sign to make colors go from blue -> red
        norm_value = float(max_value - min_value)
        print('max_value=%s min_value=%r norm_value=%r' % (max_value, min_value, norm_value))
        #print("case = ", case)
        #if normValue == 0.: # avoids division by 0.
        #    normValue = 1.

        valueSet = set()
        if vectorSize == 1:
            for value in case:
                gridResult.InsertNextValue(value)
                #if len(valueSet) < 20:
                    #valueSet.add(value)
        else:  # vectorSize=3
            for value in case:
                print(value)
                gridResult.InsertNextTuple3(*value)  # x,y,z

        nValueSet = len(valueSet)
        return norm_value, nValueSet

    def final_grid_update(self, gridResult, key, subtitle, label):
        (subcaseID, resultType, vectorSize, location, data_format) = key
        npoints = self.nPoints()
        ncells = self.nCells()

        if location == 'centroid' and self.is_centroidal:
            #self.grid.GetPointData().Reset()
            self.grid.GetCellData().SetScalars(gridResult)
            print("***centroidal plotting vector=%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
            self.grid.Modified()
        elif location == 'node' and self.is_nodal:
            self.grid.GetCellData().Reset()
            if vectorSize == 1:
                print("***node plotting vector=%s skipping (centroid/node) - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                self.grid.GetPointData().SetScalars(gridResult)
                self.grid.Modified()
            elif vectorSize == 3:
                #print("***node vector=%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                self.grid.GetPointData().SetVectors(gridResult)
                self.grid.Modified()
                #print("***node skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID,resultType,subtitle,label))
            else:
                #print("***node skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID, resultType, subtitle, label))
                raise RuntimeError(vectorSize)
        else:
            raise RuntimeError(location)
            self.log_info("***D%s skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" % (location, subcaseID, resultType, subtitle, label))
            self.scalarBar.SetVisibility(False)

    def incrementCycle(self, result_name=False):
        found_case = False
        if result_name is not False and result_name is not None:
            for icase, cases in sorted(self.resultCases.iteritems()):
                if result_name == cases[1]:
                    found_case = True
                    self.iCase = icase

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
            self.log_error("No Results found.  Many results are not supported in the GUI.\n")
            self.scalarBar.SetVisibility(False)
            foundCases = False
        #print("next key = %s" % key)
        return foundCases

    def build_grid_result(self, vectorSize, location):
        gridResult = vtk.vtkFloatArray()
        #emptyResult = vtk.vtkFloatArray()

        gridResult.SetNumberOfComponents(vectorSize)
        if location == 'centroid' and self.is_centroidal:
            #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
            gridResult.Allocate(self.nElements, 1000)
        elif location == 'node' and self.is_nodal:
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
        self.scalarBar.SetLabelFormat(data_format)

        nvalues = 11
        if (Title in ['Element_ID', 'Eids', 'Region'] and (max_value - min_value + 1) < 11):
            nvalues = int(max_value - min_value) + 1
            #ncolors = nvalues
            #if nvalues < 5:
                #ncolors = 5
            #print("need to adjust axes...max_value=%s" % max_value)
        #if data_format == '%.0f' and maxValue>

        self.scalarBar.SetNumberOfLabels(nvalues)
        self.scalarBar.SetMaximumNumberOfColors(nvalues)
        self.scalarBar.Modified()
