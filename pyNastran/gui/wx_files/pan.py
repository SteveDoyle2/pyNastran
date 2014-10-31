from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import range
import platform

import wx
import vtk
#from numpy import zeros, ones
from numpy import ndarray, amax, amin

import pyNastran
version = pyNastran.__version__

#from pyNastran.gui.mouseStyle import MouseStyle
from pyNastran.gui.wx_files.actionsControl import pyWidget
from pyNastran.gui.wx_files.gui_wx_common import GuiCommon

from pyNastran.gui.formats import (NastranIO, Cart3dIO, PanairIO, LaWGS_IO, STL_IO, TetgenIO, Usm3dIO, Plot3d_io,
    is_nastran, is_cart3d, is_panair, is_lawgs, is_stl, is_tetgen, is_usm3d, is_plot3d)

def getScreenCorner(x, y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen, yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen - x) // 2
    yCorner = (yScreen - y) // 2
    return(xCorner, yCorner)


class Pan(wx.Panel, GuiCommon, NastranIO, Cart3dIO, LaWGS_IO, PanairIO, STL_IO, TetgenIO, Usm3dIO, Plot3d_io):
    def __init__(self, *args, **kwargs):

        self.grid = vtk.vtkUnstructuredGrid()
        #gridResult = vtk.vtkFloatArray()
        #self.emptyResult = vtk.vtkFloatArray()
        #self.vectorResult = vtk.vtkFloatArray()
        self.grid2 = vtk.vtkUnstructuredGrid()

        # edges
        self.edgeActor = vtk.vtkActor()
        self.edgeMapper = vtk.vtkPolyDataMapper()

        self.is_edges = kwargs['is_edges']
        self.is_nodal = kwargs['is_nodal']
        self.is_centroidal = kwargs['is_centroidal']
        self.magnify = kwargs['magnify']
        self.debug = True

        gui_parent = kwargs['gui_parent']
        if 'log' in kwargs:
            self.log = kwargs['log']
            del kwargs['log']
        del kwargs['gui_parent']
        del kwargs['is_edges']
        del kwargs['is_nodal']
        del kwargs['is_centroidal']
        del kwargs['magnify']
        #del kwargs['rotation']
        wx.Panel.__init__(self, *args, **kwargs)
        GuiCommon.__init__(self)
        NastranIO.__init__(self)
        Cart3dIO.__init__(self)
        LaWGS_IO.__init__(self)
        PanairIO.__init__(self)
        STL_IO.__init__(self)
        TetgenIO.__init__(self)
        Usm3dIO.__init__(self)


        self.nCases = 0

        #print("is_edges = %r" % self.is_edges)
        self.widget = pyWidget(self, -1)

        window = self.widget.GetRenderWindow()
        self.iren = vtk.vtkRenderWindowInteractor()

        if platform.system == 'Windows':
            self.iren.SetRenderWindow(window)

        self.iText = 0
        self.textActors = {}

    def log_command(self, msg):
        pass

    def get_renderer(self):
        return self.rend

    def createTriAxes(self):
        pass

    def onSetToWireframe(self, event):
        if self.infile_name is not None:
            self.widget.Wireframe()
            self.widget.Update()

    def onSetToSurface(self, event):
        if self.infile_name is not None:
            self.widget.Surface()
            self.widget.Update()

    def onSetToFlatShading(self, event):  # Flat
        if self.infile_name is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0, numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToFlat()
            self.widget.Render()

    def onSetToGouraudShading(self, event):  # Gouraud
        if self.infile_name is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0, numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToGouraud()
            self.widget.Render()

    def onSetToPhongShading(self, event):  # Phong
        if self.infile_name is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0, numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToPhong()
            self.widget.Render()

    def getColors(self):
        pass
        #Merger = vtk.vtkMergeFilter()
        #Filter = vtk.vtkGeometryFilter()
        #Filter.SetInput(self.aQuadGrid)

        #Filter = vtk.vtkRotationFilter()
        #Filter.SetCenter(0., 0., 0.)
        #Filter.SetNumberOfCopies(1)
        #Filter.SetInput(aQuadGrid)
        #Filter.Update()

    def onFlipEdges(self, event):
        self.is_edges = not(self.is_edges)
        self.edgeActor.SetVisibility(self.is_edges)
        #self.edgeActor.GetProperty().SetColor(0, 0, 0)  # cart3d edge color isn't black...
        self.edgeActor.Modified()
        self.widget.Update()
        self.Refresh()

    def get_edges(self):
        edges = vtk.vtkExtractEdges()
        edges.SetInput(self.grid)
        self.edgeMapper.SetInput(edges.GetOutput())

        self.edgeActor.SetMapper(self.edgeMapper)
        self.edgeActor.GetProperty().SetColor(0, 0, 0)

        prop = self.edgeActor.GetProperty()
        self.edgeActor.SetVisibility(self.is_edges)
        self.rend.AddActor(self.edgeActor)

    def addGeometry(self):
        print("addGeometry")
        self.aQuadMapper = vtk.vtkDataSetMapper()
        self.aQuadMapper.SetInput(self.grid)

        #lut = vtk.vtkLookupTable()
        #aQuadMapper.SetLookupTable(lut)

        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        #geometryActor.GetOutput().ReleaseDataFlagOn()
        geometryActor.SetMapper(self.aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        prop = geometryActor.GetProperty()
        #prop = geometryActor.GetProperty()
        #prop.SetColor(0.9, 0.9, 0.9)
        #prop.SetColor(1., 1., 1.)
        #prop.SetAmbient(0.9)
        #prop.SetDiffuse(0.1)
        #prop.SetSpecular(0.1)

        prop.SetDiffuseColor(1, 0, 0)  # red
        #prop = geometryActor.SetBackfaceProperty(prop)

        #geometryActor.GetProperty().SetBackfaceProperty(1, 0, 0) # red
        #geometryActor.GetProperty().BackfaceCullingOn()  # hidges elements that have normals not facing camera
        #geometryActor.GetProperty().SetLineWidth(0.5)
        self.rend.AddActor(geometryActor)
        #self.rend.TwoSidedLightingOn()
        #self.rend.LightFollowCameraOn()

    def addAltGeometry(self):
        print("addAltGeometry")
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid2)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 1, 0)  # green
        self.rend.AddActor(geometryActor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

    def buildLookupTable2(self):
        scalarBar = vtkScalarBarActor()
        scalarBar.SetLookupTable(mapper.GetLookupTable())
        scalarBar.SetTitle("Title")
        scalarBar.SetNumberOfLabels(4)
        #scalarBar.GetLabelTextProperty().SetFontSize(8)

        # Create a lookup table to share between the mapper and the scalarbar
        hueLut = vtkLookupTable()

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

        drange = [10., 20.]
        self.colorFunction.AddRGBPoint(drange[0], 0.0, 0.0, 1.0)
        self.colorFunction.AddRGBPoint(drange[1], 1.0, 0.0, 0.0)

        self.scalarBar.SetTitle("Title1")
        self.scalarBar.SetLookupTable(self.colorFunction)
        self.scalarBar.SetOrientationToVertical()
        #self.scalarBar.GetLabelTextProperty().SetFontSize(8)

        self.scalarBar.SetHeight(0.9)
        self.scalarBar.SetWidth(0.20)  # the width is set first
        # after the width is set, this is adjusted
        self.scalarBar.SetPosition(0.77, 0.1)
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
        #propLabel.SetFontSize(8)

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
        sphere1.GetProperty().SetColor(1, 1, 1)
        sphere1.GetProperty().SetAmbient(1.0)
        sphere1.GetProperty().SetDiffuse(0)
        sphere1.GetProperty().SetSpecular(0)

    def main(self):
        print("main builder")
        window = self.widget.GetRenderWindow()
        window.AddRenderer(self.rend)

        self.addGeometry()
        self.addAltGeometry()
        self.buildLookupTable()
        #self.startWireframeLight()
        self.createTriAxes()

        textSize = 14 * self.magnify
        self.createText([5, 50], 'Max  ', textSize)  # text actor 0
        self.createText([5, 35], 'Min  ', textSize)  # text actor 1
        self.createText([5, 20], 'Word1', textSize)  # text actor 2
        self.createText([5, 5], 'Word2', textSize)  # text actor 3
        #self.createText([5,35],'Yet Again',  textSize)

        # Create the usual rendering stuff.
        self.get_edges()
        if self.is_edges:
            prop = self.edgeActor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edgeActor.GetProperty()
            prop.EdgeVisibilityOff()

        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(.1, .2, .4)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        self.rend.ResetCamera()

        #iren.Start()

    def createText(self, position, label, textSize=18, movable=False):
        # create a text actor
        txt = vtk.vtkTextActor()
        txt.SetInput(label)
        txtprop = txt.GetTextProperty()
        #txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(textSize)
        txtprop.SetColor(1, 1, 1)
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
        self.iText += 1

    def TurnTextOff(self):
        for (i, text) in iteritems(self.textActors):
            text.VisibilityOff()

    def TurnTextOn(self):
        for (i, text) in iteritems(self.textActors):
            text.VisibilityOn()

    def getWindow(self):
        return self.widget.GetRenderWindow()

    def getWindowName(self, winName=''):
        if self.infile_name:
            return "pyNastran v%s - %s" % (version, self.infile_name)
        else:
            return "pyNastran v%s" % (version)

    def setWindowName(self, winName=''):
        window = self.getWindow()
        window.SetWindowName(
            "pyNastran v%s - %s" % (version, self.infile_name))

    def log_info(self, msg):
        print(msg)

    def log_debug(self, msg):
        print(msg)

    def onKeyPress(self, obj, event):
        rwi = obj
        key = rwi.GetKeySym()
        #print "*Pressed %s" %(key)

    def buildVTK(self, infile_name, dirname=None):
        self.infile_name = infile_name

        self.rend = vtk.vtkRenderer()
        self.scalarBar = vtk.vtkScalarBarActor()
        self.load_nastran_geometry(self.infile_name, dirname)
        self.main()

        #self.renWin = vtk.vtkRenderWindow()
        #renWin = window

        cam = self.rend.GetActiveCamera()
        mouseArgs = {'pipeline': self, 'camera': cam}

        xSize = 800
        ySize = 600
        (x, y) = getScreenCorner(xSize, ySize)
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
