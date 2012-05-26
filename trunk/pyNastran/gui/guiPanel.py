import os
import platform

import wx
import vtk
from numpy import zeros,ones

import pyNastran
version = pyNastran.__version__

from mouseStyle import MouseStyle
from actionsControl import pyWidget
from nastranIO import NastranIO

def getScreenCorner(x,y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen,yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen-x)//2
    yCorner = (yScreen-y)//2
    return(xCorner,yCorner)

class Pan(wx.Panel,NastranIO):
    def __init__(self, *args, **kwargs):
        isEdges=kwargs['isEdges']
        del kwargs['isEdges']
        wx.Panel.__init__(self, *args, **kwargs)
        NastranIO.__init__(self)
        #isEdges = False
        print "isEdges = ",isEdges
        self.isEdges = isEdges # surface wireframe
        self.widget = pyWidget(self, -1)

        window = self.widget.GetRenderWindow()
        self.iren = vtk.vtkRenderWindowInteractor()

        if platform.system=='Windows':
            self.iren.SetRenderWindow(window)
        
        self.iText = 0
        self.textActors = {}

    def createTriAxes(self):
        pass

    def DisplayEdges(self,event):
        self.isEdges = not(self.isEdges)
        if 0:
            self.getEdges()
        
        try:
            hasEdgeActor = hasattr(self,edgeActor)
        except:
            return
        
        if hasEdgeActor:
            prop = self.edgeActor.GetProperty()
            #print "dir(prop) = ",dir(prop)
            print "visible = ",prop.GetEdgeVisibility()
            if self.isEdges:
                prop.EdgeVisibilityOn()
                print "edges are now on\n"
            else:
                prop.EdgeVisibilityOff()
                print "edges are now off\n"
            prop.Modified()
            #prop.Update()
            self.edgeActor.Modified()
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
                    raise
                #actor.GetProperty().SetLineWidth(0.0)
        window = self.widget.GetRenderWindow()
        window.Render()
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
        self.aQuadMapper = vtk.vtkDataSetMapper()
        self.aQuadMapper.SetInput(self.grid)
        
        #lut = vtk.vtkLookupTable()
        #aQuadMapper.SetLookupTable(lut)

        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        prop = geometryActor.GetProperty()
        #prop = geometryActor.GetProperty()
        #prop.SetColor(0.9,0.9,0.9)
        #prop.SetColor(1.,1.,1.)
        #prop.SetAmbient(0.9)
        #prop.SetDiffuse(0.1)
        #prop.SetSpecular(0.1)

        prop.SetDiffuseColor(1, 0, 0) # red
        #prop = geometryActor.SetBackfaceProperty(prop)
        
        #geometryActor.GetProperty().SetBackfaceProperty(1, 0, 0) # red
        #geometryActor.GetProperty().BackfaceCullingOn()  # hidges elements that have normals not facing camera
        #geometryActor.GetProperty().SetLineWidth(0.5)
        self.rend.AddActor(geometryActor)
        #self.rend.TwoSidedLightingOn()
        #self.rend.LightFollowCameraOn()

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

        self.scalarBar.SetTitle("Title1")
        self.scalarBar.SetLookupTable(self.colorFunction)
        self.scalarBar.SetOrientationToVertical()

        self.scalarBar.SetHeight(0.9)
        self.scalarBar.SetWidth(0.20) # the width is set first
        self.scalarBar.SetPosition(0.77, 0.1) # after the width is set, this is adjusted
        #self.scalarBar.SetPosition2(0.1, 0.3)
        #print self.scalarBar.GetPosition()

        propTitle = vtk.vtkTextProperty()
        propTitle.SetFontFamilyToArial()
        #propTitle.ItalicOff()
        propTitle.BoldOn()
        propTitle.ShadowOn()

        propLabel = vtk.vtkTextProperty()
        propLabel.BoldOff()
        propLabel.ShadowOn()

        scalar_range = self.grid.GetScalarRange()
        self.aQuadMapper.SetScalarRange(scalar_range)
        self.aQuadMapper.SetLookupTable(self.colorFunction)

        #self.scalarBar.SetTitleTextProperty(propTitle);
        #self.scalarBar.SetLabelTextProperty(propLabel);
        self.scalarBar.SetLabelFormat("%i")
        
        # allows 0-1 to be nice number when ranging values (gotta pick something)
        self.scalarBar.SetNumberOfLabels(11)
        self.scalarBar.SetMaximumNumberOfColors(11)
        
        #self.scalarBar.VisibilityOff()  # first load -> scalar bar off
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

    def startWireframeLight(self):
        #light = vtk.vtkLight()
        #light.SetFocalPoint(self.grid.GetOutput().GetCenter())
        #self.rend.AddLight(light)
        sphere1.GetProperty().SetColor(1,1,1)
        sphere1.GetProperty().SetAmbient(1.0)
        sphere1.GetProperty().SetDiffuse(0)
        sphere1.GetProperty().SetSpecular(0)

    def main(self):
        print "main builder"
        window = self.widget.GetRenderWindow()
        window.AddRenderer(self.rend)

        self.addGeometry()
        self.addAltGeometry()
        self.buildLookupTable()
        #self.startWireframeLight()
        self.createTriAxes()
        
        textSize = 15
        self.createText([5,50],'Max  ', textSize) # text actor 0
        self.createText([5,35],'Min  ', textSize) # text actor 1
        self.createText([5,20],'Word1',textSize) # text actor 2
        self.createText([5,5 ],'Word2',textSize) # text actor 3
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
        
        #print "dir(text) = ",dir(txt)
        txt.VisibilityOff()

        #txt.SetDisplayPosition(5,5) # bottom left
        #txt.SetDisplayPosition(5,95)

        #print dir(txt)
        #txt.SetPosition(0.1,0.5)

        # assign actor to the renderer
        self.rend.AddActor(txt)
        self.textActors[self.iText] = txt
        self.iText+=1
        
    def TurnTextOff(self):
        for i,text in self.textActors.iteritems():
            text.VisibilityOff()
        
    def TurnTextOn(self):
        for i,text in self.textActors.iteritems():
            text.VisibilityOn()

    def getWindow(self):
        return self.widget.GetRenderWindow()

    def getWindowName(self,winName=''):
        return "pyNastran v%s - %s" %(version,self.bdfFileName)

    def setWindowName(self,winName=''):
        window = self.getWindow()
        window.SetWindowName("pyNastran v%s - %s" %(version,self.bdfFileName))
    
    def isStress(self,op2,ID):
        if ID in op2.solidStress or ID in op2.plateStress or ID in op2.compositePlateStress or ID in op2.barStress  or ID in op2.beamStress or ID in op2.rodStress:
            return True
        return False

    def incrementCycle(self):
        if self.iCase is not self.nCases:
            self.iCase +=1
        else:
            self.iCase = 0

        if len(self.caseKeys)>0:
            key = self.caseKeys[self.iCase]
        
            if key[2] ==3: # vector size=3 -> vector
                self.incrementCycle()
            foundCases = True
        else:
            print "No Results found.  Many results are not supported in the GUI.\n"
            foundCases = False
        #print "next key = ",key
        return foundCases

    def cycleResults(self):
        if self.nCases==0:
            return

        foundCases = self.incrementCycle()
        if foundCases:
            self.gridResult = vtk.vtkFloatArray()

            key = self.caseKeys[self.iCase]
            case = self.resultCases[key]
            (subcaseID,resultType,vectorSize,location,dataFormat) = key

            if location=='centroid':
                #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
                self.gridResult.Allocate(self.nElements,1000)
            #else: # node
                #allocationSize = vectorSize*location (where location='node'-> self.nNodes)
                #self.gridResult.Allocate(self.nNodes*vectorSize,1000)
                #self.gridResult.SetNumberOfComponents(vectorSize)

            #self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]
            caseName = self.iSubcaseNameMap[subcaseID]
            (subtitle,label) = caseName

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
            #print "case = ",case
            #if normValue==0.: # avoids division by 0.
            #    normValue = 1.

            if vectorSize==1:
                for value in case:
                    self.gridResult.InsertNextValue(value)
                ###
            else: # vectorSize=3
                pass
                #for value in case:
                #    self.gridResult.InsertNextTuple3(value)  # x,y,z
                ###
            ###
            print "max=%g min=%g norm=%g\n" %(maxValue,minValue,normValue)

            self.textActors[0].SetInput('Max:  %g' %(maxValue)) # max
            self.textActors[1].SetInput('Min:  %g' %(minValue)) # min
            self.textActors[2].SetInput('Subcase=%s Subtitle: %s' %(subcaseID,subtitle)) # info
            self.textActors[3].SetInput('Label: %s'               %(label)) # info
            self.UpdateScalarBar(resultType,minValue,maxValue,dataFormat)

            # @todo results can only go from centroid->node and not back to centroid
            if location=='centroid':
                self.grid.GetCellData().SetScalars(self.gridResult)
                #self.grid.GetPointData().SetScalars(self.emptyResult) # causes a crash
                self.grid.Modified()
            else:
                #self.grid.GetCellData().SetScalars(self.emptyResult)
                #if vectorSize==1:
                #    self.grid.GetPointData().SetScalars(self.gridResult)
                #else:
                #    self.grid.GetPointData().SetScalars(self.gridResult)
                #self.grid.Modified()
                pass
            ###
        ### if results       

    def OnKeyPress(self,obj,event):
        rwi = obj
        key = rwi.GetKeySym()
        #print "*Pressed %s" %(key)

    def buildVTK(self,bdfFileName=None,dirname=None):
        self.bdfFileName = bdfFileName

        self.rend = vtk.vtkRenderer()
        self.scalarBar = vtk.vtkScalarBarActor()
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
        #iren = vtk.vtkRenderWindowInteractor()
        #mouseStyle = MouseStyle(mouseArgs,iren)
        #iren.SetInteractorStyle(mouseStyle)
        #window = self.widget.GetRenderWindow()
        #iren.SetRenderWindow(window)

        #iren.AddObserver("KeyPressEvent", self.OnKeyPress)

        #print "type(ren) = ",type(self.rend)
        #self.rend.GetActiveCamera().Zoom(2.0)

        #self.setWindowName()
        # Render the scene and start interaction.
        #iren.Initialize()
        #iren.Start()
        #window.Render()

