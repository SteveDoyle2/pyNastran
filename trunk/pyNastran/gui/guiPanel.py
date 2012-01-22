import os
import wx
import vtk
from vtk import (vtkTriangle,vtkQuad,vtkTetra,vtkWedge,vtkHexahedron,
                 vtkQuadraticTriangle,vtkQuadraticQuad,vtkQuadraticTetra,vtkQuadraticWedge,vtkQuadraticHexahedron)
from numpy import zeros,ones

import pyNastran
version = pyNastran.__version__

from pyNastran.bdf.bdf import *
from pyNastran.op2.op2 import OP2
from mouseStyle import MouseStyle
from actionsControl import pyWidget

def getScreenCorner(x,y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen,yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen-x)//2
    yCorner = (yScreen-y)//2
    return(xCorner,yCorner)

class Pan(wx.Panel):
    def __init__(self, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        isEdges = True
        self.isEdges = isEdges # surface wireframe
        self.widget = pyWidget(self, -1)

        window = self.widget.GetRenderWindow()
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(window)
        self.iText = 0
        self.textActors = {}

    def DisplayEdges(self,event):
        self.isEdges = not(self.isEdges)
        if 0:
            self.getEdges()

        if 1:
            prop = self.edgeActor.GetProperty()
            print "dir(prop) = ",dir(prop)
            print "visible = ",prop.GetEdgeVisibility()
            if self.isEdges:
                prop.EdgeVisibilityOn()
                print "edges are now on\n"
            else:
                prop.EdgeVisibilityOff()
                print "edges are now off\n"
            prop.Modified()
            #prop.Update()
            #self.edgeActor.Modified()
            #self.edgeActor.Update()
        if 0:
            print "dir(widget) = ",dir(self.widget)

            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                print "iactor = ",i
                actor = actors.GetNextItem()
                
                try:
                    if self.isEdges:
                        actor.getProperty().EdgeVisibilityOn()
                    else:
                        actor.getProperty().EdgeVisibilityOff()
                    print "set the edges..."
                except:
                    pass
                #actor.GetProperty().SetLineWidth(0.0)
        self.widget.Render()
        self.widget.Update()


    def SetToWireframe(self,event):
        if self.bdfFileName is not None:
            self.widget.Wireframe()
            self.widget.Update()

    def SetToSurface(self,event):
        if self.bdfFileName is not None:
            self.widget.Surface()
            self.widget.Update()

    def SetToFlatShading(self,event): # Flat
        if self.bdfFileName is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToFlat()
            self.widget.Render()
        ###

    def SetToGouraudShading(self,event): # Gouraud
        if self.bdfFileName is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToGouraud()
            self.widget.Render()
        ###

    def SetToPhongShading(self,event): # Phong
        if self.bdfFileName is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToPhong()
            self.widget.Render()
        ###

    def WireframeTemp(self):
        """Sets the current actor representation as wireframe.
        """
        actors = self.widget._CurrentRenderer.GetActors()
        numActors = actors.GetNumberOfItems()
        actors.InitTraversal()
        for i in range(0,numActors):
            actor = actors.GetNextItem()
            actor.GetProperty().SetRepresentationToWireframe()
        self.Render()

    def getColors(self):
        pass
        #Merger = vtk.vtkMergeFilter()
        #Filter = vtk.vtkGeometryFilter()
        #Filter.SetInput(self.aQuadGrid)
        
        #Filter = vtk.vtkRotationFilter()
        #Filter.SetCenter(0.,0.,0.)
        #Filter.SetNumberOfCopies(1)
        #Filter.SetInput(aQuadGrid)
        #Filter.Update()
        
    def getEdges(self):
        edges = vtk.vtkExtractEdges()
        edges.SetInput(self.grid)
        self.edgeMapper = vtk.vtkPolyDataMapper()
        self.edgeMapper.SetInput(edges.GetOutput())
        #self.edgeMapper.EdgeVisibilityOff()

        self.edgeActor = vtk.vtkActor()
        self.edgeActor.SetMapper(self.edgeMapper)
        self.edgeActor.GetProperty().SetColor(0,0,0)
        prop = self.edgeActor.GetProperty()
        #prop.SetLineWidth(0.0)
        if self.isEdges:
            prop.EdgeVisibilityOn()
        else:
            prop.EdgeVisibilityOff()

        self.rend.AddActor(self.edgeActor)

        print "visible = ",prop.GetEdgeVisibility()
        #self.edgeActor.Update()


    def addGeometry(self):
        print "addGeometry"
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid)
        
        #lut = vtk.vtkLookupTable()
        #aQuadMapper.SetLookupTable(lut)

        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0) # red
        #geometryActor.GetProperty().SetBackfaceProperty(1, 0, 0) # red
        #geometryActor.GetProperty().BackfaceCullingOn()  # hidges elements that have normals not facing camera
        #geometryActor.GetProperty().SetLineWidth(0.5)
        self.rend.AddActor(geometryActor)

    def addAltGeometry(self):
        print "addAltGeometry"
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid2)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 1, 0) # green
        self.rend.AddActor(geometryActor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

    def buildLookupTable2(self):

        scalarBar = vtkScalarBarActor()
        scalarBar.SetLookupTable(mapper.GetLookupTable())
        scalarBar.SetTitle("Title")
        scalarBar.SetNumberOfLabels(4)
 
        # Create a lookup table to share between the mapper and the scalarbar
        hueLut =vtkLookupTable()

        hueLut.SetTableRange(0, 1)
        hueLut.SetHueRange(0, 1)
        hueLut.SetSaturationRange(1, 1)
        hueLut.SetValueRange(1, 1)
        hueLut.Build()
        mapper.SetLookupTable(hueLut)
        scalarBar.SetLookupTable(hueLut)
        mapper.ScalarVisibilityOn()
        mapper.SetScalarModeToUsePointData()
        mapper.SetColorModeToMapScalars()

    def buildLookupTable(self):
        self.colorFunction = vtk.vtkColorTransferFunction()
        self.colorFunction.SetColorSpaceToHSV()
        self.colorFunction.HSVWrapOff()
        
        drange = [10.,20.]
        self.colorFunction.AddRGBPoint(drange[0], 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(drange[1], 1.0, 0.0, 0.0)

        self.scalarBar = vtk.vtkScalarBarActor()
        self.scalarBar.SetTitle("Title1")
        self.scalarBar.SetLookupTable(self.colorFunction)
        self.scalarBar.SetOrientationToVertical()

        self.scalarBar.SetHeight(0.9)
        self.scalarBar.SetWidth(0.20) # the width is set first
        self.scalarBar.SetPosition(0.77, 0.1) # after the width is set, this is adjusted
        #self.scalarBar.SetPosition2(0.1, 0.3)
        print self.scalarBar.GetPosition()

        propTitle = vtk.vtkTextProperty()
        propTitle.SetFontFamilyToArial()
        #propTitle.ItalicOff()
        propTitle.BoldOn()
        propTitle.ShadowOn()

        propLabel = vtk.vtkTextProperty()
        propLabel.BoldOff()
        propLabel.ShadowOn()

        #self.scalarBar.SetTitleTextProperty(propTitle);
        #self.scalarBar.SetLabelTextProperty(propLabel);
        self.scalarBar.SetLabelFormat("%i")
        
        # allows 0-1 to be nice number when ranging values (gotta pick something)
        self.scalarBar.SetNumberOfLabels(11)
        self.scalarBar.SetMaximumNumberOfColors(11)
        
        visibility = True
        if visibility:
            self.scalarBar.VisibilityOn()
        else:
            self.scalarBar.VisibilityOff()
        #self.scalarBar.ShadowOn()
        #self.scalarBar.RepositionableOn()
        self.rend.AddActor(self.scalarBar)
        #return scalarBar

    def UpdateScalarBar(self,Title,minValue,maxValue,dataFormat):
        """
        @param Title the scalar bar title
        @param minValue the blue value
        @param maxValue the red value
        @param dataFormat '%g','%f','%i',etc.
        """
        #drange = [10.,20.]
        self.colorFunction.RemoveAllPoints()
        self.colorFunction.AddRGBPoint(minValue, 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(maxValue, 1.0, 0.0, 0.0)
        #self.scalarBar.SetLookupTable(self.colorFunction)

        self.scalarBar.SetTitle(Title)
        self.scalarBar.SetLabelFormat(dataFormat)
        
        nValues = 11
        if Title=='Element_ID' and (maxValue-minValue+1)<11:
            nValues = int(maxValue-minValue)+1
            #print "need to adjust axes...maxValue=%s" %(maxValue)
        #if dataFormat=='%.0f' and maxValue>

        self.scalarBar.SetNumberOfLabels(nValues)
        self.scalarBar.SetMaximumNumberOfColors(nValues)
        self.scalarBar.Modified()

    def Update(self):
        #print '\n'.join(dir(self.widget))
        window = self.widget.GetRenderWindow()
        self.rend.ResetCamera()
        window.Render()
        self.widget.Update()
        #self.rend.Update()
        #self.renWin.Render()
        
    def main(self):
        print "main builder"
        window = self.widget.GetRenderWindow()
        window.AddRenderer(self.rend)

        self.addGeometry()
        self.addAltGeometry()
        
        textSize = 15
        self.createText([5,35],'Max: ', textSize) # text actor 0
        self.createText([5,20],'Min: ', textSize) # text actor 1
        self.createText([5,5 ],'Word: ',textSize) # text actor 2
        #self.createText([5,35],'Yet Again',  textSize)

        # Create the usual rendering stuff.
        if self.isEdges:
            self.getEdges()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(.1,.2,.4)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        self.rend.ResetCamera()

        #iren.Start()

    def createText(self,position,label,textSize=18,movable=False):
        # create a text actor
        txt = vtk.vtkTextActor()
        txt.SetInput(label)
        txtprop=txt.GetTextProperty()
        #txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(textSize)
        txtprop.SetColor(1,1,1)
        txt.SetDisplayPosition(*position)
        #txt.SetDisplayPosition(5,5) # bottom left
        #txt.SetDisplayPosition(5,95)

        #print dir(txt)
        #txt.SetPosition(0.1,0.5)

        # assign actor to the renderer
        self.rend.AddActor(txt)
        self.textActors[self.iText] = txt
        self.iText+=1
        
    def getWindow(self):
        return self.widget.GetRenderWindow()

    def getWindowName(self,winName=''):
        return "pyNastran v%s - %s" %(version,self.bdfFileName)
        return "pyNastran v%s - %s" %(version,'solid.bdf')

    def setWindowName(self,winName=''):
        window = self.getWindow()
        window.SetWindowName("pyNastran v%s - %s" %(version,self.bdfFileName))
    
    def isStress(self,op2,ID):
        if ID in op2.solidStress or ID in op2.plateStress or ID in op2.compositePlateStress or ID in op2.barStress  or ID in op2.beamStress or ID in op2.rodStress:
            return True
        return False
        
    def loadResults(self,op2FileName):
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        op2 = OP2(op2FileName)
        op2.readOP2()
        #print op2.printResults()
        
        #case = op2.displacements[1]
        #print "case = ",case
        #for nodeID,translation in sorted(case.translations.items()):
        #    print "nodeID=%s t=%s" %(nodeID,translation)
        #self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = op2.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = op2.iSubcaseNameMap
        
        nidsSet = False # set to True to disable nodeIDs
        eidsSet = False
        for ID in subcaseIDs:
            if not nidsSet:
                nids = zeros(self.nNodes,'d')
                for nid,nid2 in self.nidMap.items():
                    nids[nid2] = nid
                cases[(ID,'Node_ID',1,'node','%.0f')] = nids
                nidsSet = True

            if not eidsSet:
                eids = zeros(self.nElements,'d')
                for eid,eid2 in self.eidMap.items():
                    eids[eid2] = eid
               
                eKey = (ID,'isElementOn',1,'centroid','%.0g')
                cases[(ID,'Element_ID',1,'centroid','%.0f')] = eids
                cases[eKey] = zeros(self.nElements) # is the element supported
                eidsSet = True

            if ID in op2.displacements: # not correct?
                case = op2.displacements[ID]
                key = (ID,'DisplacementX',3,'node','%g')
                cases[key] = case.translations

            if ID in op2.temperatures:
                case = op2.temperatures[ID]
                print case
                temps = zeros(self.nNodes)
                key = (ID,'Temperature',1,'node','%g')
                for nid,T in case.temperatures.items():
                    print T
                    nid2 = self.nidMap[nid]
                    temps[nid2] = T
                cases[key] = temps

            if self.isStress(op2,ID):
                cases = self.fillStressCase(cases,op2,ID,eKey)
            ###
        ###
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases)-1 # number of keys in dictionary
        self.cycleResults() # start at nCase=0

    def fillStressCase(self,cases,op2,ID,eKey):
        oxx = zeros(self.nElements)
        oyy = zeros(self.nElements)
        ozz = zeros(self.nElements)

        o1  = zeros(self.nElements)
        o2  = zeros(self.nElements)
        o3  = zeros(self.nElements)
        ovm = zeros(self.nElements)

        vmWord = 'N/A'
        if ID in op2.rodStress:
            case = op2.rodStress[ID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.
                #print "bar eid=%s" %(eid)

                axial   = case.axial[eid]
                torsion = case.torsion[eid]

                oxx[eid2] = axial
                oyy[eid2] = torsion

                o1[eid2]  = max(axial,torsion)  # not really
                #oyy[eid2] = torsion
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = min(axial,torsion)

        if ID in op2.barStress:
            #self.s1    = {}
            #self.s2    = {}
            #self.s3    = {}
            #self.s4    = {}
            case = op2.barStress[ID]
            for eid in case.axial:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "bar eid=%s" %(eid)
                oxxi = case.axial[eid]
                o1i  = max(case.smax[eid])
                o3i  = min(case.smin[eid])
                
                oxx[eid2] = oxxi
                #oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                #ovm[eid2] = ovmi

        if ID in op2.plateStress:
            #self.txy    = {}
            case = op2.plateStress[ID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "plate eid=%s" %(eid)
                oxxi = case.oxx[eid]['C']
                #self.oyy[eid][nid][iLayer]

                oxxi = max(case.oxx[eid]['C'])
                oyyi = max(case.oyy[eid]['C'])
                ozzi = min(case.oxx[eid]['C'],min(case.oyy[eid]['C']))

                o1i = max(case.majorP[eid]['C'])
                o2i = max(case.minorP[eid]['C'])
                o3i = min(case.majorP[eid]['C'],min(case.minorP[eid]['C']))
                ovmi = max(case.ovmShear[eid]['C'])


                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        if ID in op2.solidStress:
            case = op2.solidStress[ID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 1.

                #print "solid eid=%s" %(eid)
                oxxi = case.oxx[eid]['C']
                oyyi = case.oyy[eid]['C']
                ozzi = case.ozz[eid]['C']
                o1i = case.o1[eid]['C']
                o2i = case.o2[eid]['C']
                o3i = case.o3[eid]['C']
                ovmi = case.ovmShear[eid]['C']
                #if ID==1:
                #    print "ovm[%s] = %s" %(eid,ovmi)
                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                ozz[eid2] = ozzi

                o1[eid2] = o1i
                o2[eid2] = o2i
                o3[eid2] = o3i
                ovm[eid2] = ovmi
            ###

        cases[(ID,'StressXX',1,'centroid','%.3f')] = oxx
        cases[(ID,'StressYY',1,'centroid','%.3f')] = oyy
        cases[(ID,'StressZZ',1,'centroid','%.3f')] = ozz

        cases[(ID,'Stress1',1,'centroid','%.3f')] = o1
        cases[(ID,'Stress2',1,'centroid','%.3f')] = o2
        cases[(ID,'Stress3',1,'centroid','%.3f')] = o3
        cases[(ID,vmWord   ,1,'centroid','%.3f')] = ovm
        return cases

    def incrementCycle(self):
        if self.iCase is not self.nCases:
            self.iCase +=1
        else:
            self.iCase = 0
        key = self.caseKeys[self.iCase]
        if key[2] ==3: # vector size=3 -> vector
            self.incrementCycle()
        
        #print "next key = ",key
        return

    def cycleResults(self):
        self.incrementCycle()
        self.gridResult = vtk.vtkFloatArray()

        key = self.caseKeys[self.iCase]
        case = self.resultCases[key]
        (subcaseID,resultType,vectorSize,location,dataFormat) = key

        if location=='centroid':
            #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
            self.gridResult.Allocate(self.nElements,1000)
        else: # node
            #allocationSize = vectorSize*location (where location='node'-> self.nNodes)
            self.gridResult.Allocate(self.nNodes,1000)

        #self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]
        caseName = self.iSubcaseNameMap[subcaseID]
        subtitle,label = caseName

        print "subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID,resultType,subtitle,label)

        for value in case:
            maxValue = value
            minValue = value
            break

        for value in case:
            maxValue = max(value,maxValue)
            minValue = min(value,minValue)
        
        # flips sign to make colors go from blue -> red
        normValue = maxValue-minValue
        print "case = ",case
        if normValue==0.: # avoids division by 0.
            normValue = 1.
        for value in case:
            self.gridResult.InsertNextValue((value-minValue)/normValue)
        print "max=%g min=%g norm=%g\n" %(maxValue,minValue,normValue)

        self.textActors[0].SetInput('Max:  %g' %(maxValue)) # max
        self.textActors[1].SetInput('Min:  %g' %(minValue)) # min
        self.textActors[2].SetInput('subcase=%s subtitle=%s' %(subcaseID,subtitle)) # info
        self.UpdateScalarBar(resultType,minValue,maxValue,dataFormat)

        # @todo results can only go from centroid->node and not back to centroid
        if location=='centroid':
            self.grid.GetCellData().SetScalars(self.gridResult)
            #self.grid.GetPointData().SetScalars(self.emptyResult) # causes a crash
        else:
            #self.grid.GetCellData().SetScalars(self.emptyResult)
            self.grid.GetPointData().SetScalars(self.gridResult)
        self.grid.Modified()

    def loadGeometry(self,bdfFileName,dirname):

        if bdfFileName is None:
            self.grid       = vtk.vtkUnstructuredGrid()
            self.gridResult = vtk.vtkFloatArray()
            self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.grid2      = vtk.vtkUnstructuredGrid()
            self.scalarBar.VisibilityOff()
            return
        else:
            self.grid.Reset()
            self.gridResult.Reset()
            self.grid2.Reset()
            try:
                del self.resultCases
                del self.caseKeys
                del self.iCase
                del self.nCases
                del self.iSubcaseNameMap
            except:
                print "cant cleanup..."
            ###
        self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()

        model = BDF()
        model.readBDF(bdfFileName,includeDir=dirname)

        nNodes    = model.nNodes()
        nElements = model.nElements()
        nCAeros   = model.nCAeros()
        self.nNodes = nNodes
        self.nElements = nElements

        #print "nNodes = ",nNodes
        print "nElements = ",nElements

        #self.aQuadGrid.Allocate(nElements+nNodes, 1000)

        if 'CONM2' in model.cardCount:
            nCONM2 = model.cardCount['CONM2']
        else:
            nCONM2 = 0
        self.grid.Allocate(nElements, 1000)
        self.grid2.Allocate(nCAeros +nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nNodes)
        self.gridResult.Allocate(nNodes,1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        i=0
        #elem.SetNumberOfPoints(nNodes)
        fraction = 1./nNodes
        for nid,node in sorted(model.nodes.items()):
            #print "i = ",i
            point = node.Position()
            #print "point = ",point
            #sys.stdout.flush()
            #aVoxel = vtk.vtkPixel()
            #print "made voxel"; sys.stdout.flush()
            points.InsertPoint(i, *point)
            self.gridResult.InsertNextValue(i*fraction)
            #print str(element)

            #elem = vtk.vtkVertex()
            #elem.GetPointIds().SetId(0, i)
            #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

            self.nidMap[nid] = i
            i+=1

        j = 0
        points2 = vtk.vtkPoints()
        points2.SetNumberOfPoints(nCAeros*4+nCONM2)
        for eid,element in sorted(model.caeros.items()):
            if isinstance(element,CAERO1):
                cpoints = element.Points()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j+1)
                elem.GetPointIds().SetId(2, j+2)
                elem.GetPointIds().SetId(3, j+3)
                points2.InsertPoint(j,   *cpoints[0])
                points2.InsertPoint(j+1, *cpoints[1])
                points2.InsertPoint(j+2, *cpoints[2])
                points2.InsertPoint(j+3, *cpoints[3])
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j+=4
        self.mapElements(points,points2,self.nidMap,model,j)

    def mapElements(self,points,points2,nidMap,model,j):
        self.eidMap = {}
        i = 0
        for eid,element in sorted(model.elements.items()):
            self.eidMap[eid] = i
            if isinstance(element,CTRIA3):
                #print "ctria3"
                elem = vtkTriangle()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTRIA6):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CQUAD4):
                nodeIDs = element.nodeIDs()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CQUAD8):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticQuad()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                else:
                    elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA4):
                elem = vtkTetra()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA10):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTetra()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                else:
                    elem = vtkTetra()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CPENTA6):
                elem = vtkWedge()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element,CPENTA15):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6,  nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7,  nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8,  nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9,  nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                else:
                    elem = vtkWedge()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CHEXA8):
                nodeIDs = element.nodeIDs()
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CHEXA20):
                nodeIDs = element.nodeIDs()
                #print "nodeIDs = ",nodeIDs
                if None not in nodeIDs:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8,  nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9,  nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                    elem.GetPointIds().SetId(15, nidMap[nodeIDs[15]])
                    elem.GetPointIds().SetId(16, nidMap[nodeIDs[16]])
                    elem.GetPointIds().SetId(17, nidMap[nodeIDs[17]])
                    elem.GetPointIds().SetId(18, nidMap[nodeIDs[18]])
                    elem.GetPointIds().SetId(19, nidMap[nodeIDs[19]])
                else:
                    elem = vtkHexahedron()

                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,LineElement) or isinstance(element,SpringElement):
                elem = vtk.vtkLine()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            ###
            elif isinstance(element,CONM2): # not perfectly located
                nid  = element.Nid()
                c    = element.Centroid()
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()
                #elem.SetRadius(1.0)
                #print str(element)

                points2.InsertPoint(j,     *c)
                elem.GetPointIds().SetId(0, j)
                #elem.SetCenter(points.GetPoint(nidMap[nid]))
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j+=1
            else:
                print "skipping %s" %(element.type)
            i+=1
        ###
        self.grid.SetPoints(points)
        self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        print "updated grid"

    def OnKeyPress(self,obj,event):
        rwi = obj
        key = rwi.GetKeySym()
        print "*Pressed %s" %(key)

    def buildVTK(self,bdfFileName=None,dirname=None):
        self.bdfFileName = bdfFileName

        self.rend = vtk.vtkRenderer()
        self.buildLookupTable()
        self.loadGeometry(self.bdfFileName,dirname)
        self.main()

        #self.renWin = vtk.vtkRenderWindow()
        #renWin = window
        
        cam = self.rend.GetActiveCamera()
        mouseArgs = {'pipeline':self,'camera':cam}

        xSize = 800
        ySize = 600
        (x,y) = getScreenCorner(xSize,ySize)
        #window.SetSize(xSize,ySize)
        #window.SetPosition(x,y)

        #iren = wxVTKRenderWindowInteractor(self,-1)
        iren = vtk.vtkRenderWindowInteractor()
        mouseStyle = MouseStyle(mouseArgs,iren)
        iren.SetInteractorStyle(mouseStyle)
        window = self.widget.GetRenderWindow()
        iren.SetRenderWindow(window)

        iren.AddObserver("KeyPressEvent", self.OnKeyPress)

        #print "type(ren) = ",type(self.rend)
        #self.rend.GetActiveCamera().Zoom(2.0)

        #self.setWindowName()
        # Render the scene and start interaction.
        iren.Initialize()
        #iren.Start()
        #window.Render()

