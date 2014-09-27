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
        plotNodal = self.is_nodal
        plotCentroidal = self.is_centroidal
        #print("plotNodal=%s plotCentroidal=%s" %(plotNodal,plotCentroidal))
        #print("nCases = %i" %(self.nCases+1))
        if self.nCases == 0:
            self.scalarBar.SetVisibility(False)
            return

        foundCases = self.incrementCycle()
        if foundCases:
        #if 1:
            print("incremented case")
            #gridResult.Reset()
            gridResult = vtk.vtkFloatArray()
            emptyResult = vtk.vtkFloatArray()

            try:
                key = self.caseKeys[self.iCase]
            except:
                print('icase=%s caseKeys=%s' % (self.iCase, str(self.caseKeys)))
                raise
            case = self.resultCases[key]
            print("len(case) = %i" % len(case))
            (subcaseID, resultType, vectorSize, location, dataFormat) = key

            if location == 'centroid' and plotCentroidal:
                #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
                gridResult.Allocate(self.nElements, 1000)
            elif location == 'node' and plotNodal:
                #allocationSize = vectorSize*self.nNodes # (where location='node'-> self.nNodes)
                gridResult.Allocate(self.nNodes * vectorSize, 1000)
                gridResult.SetNumberOfComponents(vectorSize)
            else:
                print("***%s skipping" % location)

            try:
                #self.iSubcaseNameMap[self.isubcase] = [Subtitle, Label]
                caseName = self.iSubcaseNameMap[subcaseID]
            except KeyError:
                #print "cant find subcaseID=%s" % subcaseID
                caseName = ('case=NA', 'label=NA')
            (subtitle, label) = caseName

            print("subcaseID=%s resultType=%s subtitle=%r label=%r" % (subcaseID, resultType, subtitle, label))

            if isinstance(case, ndarray):
                max_value = case.max()
                min_value = case.min()
            else:
                max_value = case[0]
                min_value = case[0]
                for value in case:
                    max_value = max(value, max_value)
                    min_value = min(value, min_value)

            # flips sign to make colors go from blue -> red
            try:
                norm_value = max_value - min_value
            except:
                raise RuntimeError(resultType)
            #print "case = ",case
            #if normValue==0.: # avoids division by 0.
            #    normValue = 1.

            #valueSet = set()
            if vectorSize == 1:
                #print "minValue = ",min(case)
                for value in case:
                    gridResult.InsertNextValue(value)
                    #if len(valueSet) < 20:
                        #valueSet.add(value)
            else:  # vectorSize=3
                pass
                #for value in case:
                #    gridResult.InsertNextTuple3(value)  # x,y,z

            print("max=%g min=%g norm=%g\n" % (max_value, min_value, norm_value))

            #nValueSet = len(valueSet)
            self.update_text_actors(case, subcaseID, subtitle, min_value, max_value, label)

            self.UpdateScalarBar(resultType, min_value, max_value, dataFormat)
            #self.scalarBar.SetNumberOfLabels(nValueSet)
            #self.scalarBar.SetMaximumNumberOfColors(nValueSet)
            #prop = self.scalarBar.GetLabelTextProperty()
            #fontSize = prop.GetFontSize()
            #print "fontSize = ",fontSize
            #prop.SetFontSize(40)

            # TODO results can only go from centroid->node and not back to
            ## centroid
            #print dir(self.grid)
            #self.grid.Reset()
            if location == 'centroid' and plotCentroidal:
                #self.grid.GetPointData().Reset()
                self.grid.GetCellData().SetScalars(gridResult)
                print("***plotting vector=%s skipping (centroid/node) - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                self.grid.Modified()
            elif location == 'node' and plotNodal:
                self.grid.GetCellData().Reset()
                if vectorSize == 1:
                    print("***plotting vector=%s skipping (centroid/node) - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                    self.grid.GetPointData().SetScalars(gridResult)
                    self.grid.Modified()
                else:
                    print("***node vector=%s skipping (centroid/node) - subcaseID=%s resultType=%s subtitle=%s label=%s" % (vectorSize, subcaseID, resultType, subtitle, label))
                    self.grid.GetPointData().SetVectors(gridResult)
                    #print("***node skipping - subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID,resultType,subtitle,label))
            else:
                print("***%s skipping (centroid/node) - subcaseID=%s resultType=%s subtitle=%s label=%s" % (location, subcaseID, resultType, subtitle, label))
                self.scalarBar.SetVisibility(False)

    def UpdateScalarBar(self, Title, min_value, max_value, dataFormat):
        """
        @param Title the scalar bar title
        @param min_value the blue value
        @param max_value the red value
        @param dataFormat '%g','%f','%i',etc.
        """
        #drange = [10.,20.]
        self.colorFunction.RemoveAllPoints()
        self.colorFunction.AddRGBPoint(min_value, 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(max_value, 1.0, 0.0, 0.0)
        #self.scalarBar.SetLookupTable(self.colorFunction)

        self.scalarBar.SetTitle(Title)
        self.scalarBar.SetLabelFormat(dataFormat)

        nvalues = 11
        if (Title in ['Element_ID', 'Eids', 'Region'] and (max_value - min_value + 1) < 11):
            nvalues = int(max_value - min_value) + 1
            #ncolors = nvalues
            #if nvalues < 5:
                #ncolors = 5
            #print "need to adjust axes...max_value=%s" %(max_value)
        #if dataFormat=='%.0f' and maxValue>

        #print("Title =", Title)
        #print("nvalues =", nvalues)
        self.scalarBar.SetNumberOfLabels(nvalues)
        self.scalarBar.SetMaximumNumberOfColors(nvalues)
        self.scalarBar.Modified()
