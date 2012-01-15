#!/usr/bin/python

import wx
import vtk
from vtk.wx.wxVTKRenderWindow import wxVTKRenderWindow
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

import pyNastran
from pyNastran.bdf.bdf import *
from mouseStyle import MouseStyle
from vtk import (vtkTriangle,vtkQuad,vtkTetra,vtkWedge,vtkHexahedron,
                 vtkQuadraticTriangle,vtkQuadraticQuad,vtkQuadraticTetra,vtkQuadraticWedge,vtkQuadraticHexahedron)

#ID_OPEN = 801
ID_SAVEAS = 803
ID_ABOUT = 3

def getScreenCorner(x,y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen,yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen-x)//2
    yCorner = (yScreen-y)//2
    return(xCorner,yCorner)

version = pyNastran.__version__

class pyWidget(wxVTKRenderWindow):
    def __init__(self,*args,**kwargs):
        wxVTKRenderWindow.__init__(self, *args, **kwargs)
        self.OnChar = self.onChar2

    def ResetCamera(self):
        self.Reset()

    def GetCamera(self):
        return self._CurrentCamera

    def onChar2(self,event):
        #print "onChar2 = ",event.GetKeyCode()
        camera = self.GetCamera()
        code = event.GetKeyCode()
        if code == ord('m'): # zooming in
            camera.Zoom(1.1)
        elif code == ord('M'):  # zooming out
            camera.Zoom(0.9)

        elif code == ord('o'): # counter-clockwise
            camera.Roll(5.)
        elif code == ord('O'): # clockwise
            camera.Roll(-5.)

        elif code == ord('x'): # set x-axis
            camera.SetFocalPoint(0.,0., 0.)
            camera.SetViewUp(    0.,0., 1.)
            camera.SetPosition(  1.,0., 0.)
            self.ResetCamera()
        elif code == ord('X'): # set x-axis
            camera.SetFocalPoint(0.,0., 0.)
            camera.SetViewUp(    0.,0.,-1.)
            camera.SetPosition( -1.,0., 0.)
            self.ResetCamera()


        elif code == ord('y'): # set y-axis
            camera.SetFocalPoint(0.,0.,0.)
            camera.SetViewUp(    0.,0.,1.)
            camera.SetPosition(  0.,1.,0.)
            self.ResetCamera()
        elif code == ord('Y'): # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(    0., 0.,-1.)
            camera.SetPosition(  0.,-1., 0.)
            self.ResetCamera()

        elif code == ord('z'): # set z-axis
            camera.SetFocalPoint(0.,0.,0.)
            camera.SetViewUp(    0.,1.,0.)
            camera.SetPosition(  0.,0.,1.)
            self.ResetCamera()
        elif code == ord('Z'): # set z-axis
            camera.SetFocalPoint(0.,0., 0.)
            camera.SetViewUp(   0., -1.,0.)
            camera.SetPosition( 0., 0.,-1.)
            self.ResetCamera()

        elif code == ord('i'): # picture taking doesnt work
            self.takePicture()

        self.Update()
        self.Render()
        ###

    def takePicture(self):
        renderLarge = vtk.vtkRenderLargeImage()
        renderLarge.SetInput(self.getRenderer())
        renderLarge.SetMagnification(4)

        # We write out the image which causes the rendering to occur. If you
        # watch your screen you might see the pieces being rendered right
        # after one another.
        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(renderLarge.GetOutputPort())
        writer.SetFileName("largeImage.png")
        writer.Write()

    def getRenderer(self):
        return self.GetCurrentRenderer()
        window = self.GetRenderWindow()
        renderers = window.GetRenderers()
        return renderers[0]

class Pan(wx.Panel):
    def __init__(self, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        isEdges = False
        self.isEdges = isEdges # surface wireframe
        self.widget = pyWidget(self, -1)
        

        window = self.widget.GetRenderWindow()
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(window)

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
        self.rend = vtk.vtkRenderer()
        window.AddRenderer(self.rend)

        #self.rend = vtk.vtkRenderer()
        self.addGeometry()
        self.addAltGeometry()

        # Create the usual rendering stuff.
        if self.isEdges:
            self.getEdges()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(.1,.2,.4)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        self.rend.ResetCamera()

        #iren.Start()

    def getWindow(self):
        return self.widget.GetRenderWindow()

    def getWindowName(self,winName=''):
        return "pyNastran v%s - %s" %(version,self.bdfFileName)

    def setWindowName(self,winName=''):
        window = self.getWindow()
        window.SetWindowName("pyNastran v%s - %s" %(version,self.bdfFileName))
    
    def loadGeometry(self,bdfFileName):

        if bdfFileName is None:
            self.grid  = vtk.vtkUnstructuredGrid()
            self.grid2 = vtk.vtkUnstructuredGrid()
            return
        else:
            self.grid.Reset()
            self.grid2.Reset()

        model = BDF()
        model.readBDF(bdfFileName)

        nNodes    = model.nNodes()
        nElements = model.nElements()
        nCAeros   = model.nCAeros()

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
        nidMap = {}
        i=0
        #elem.SetNumberOfPoints(nNodes)
        for nid,node in sorted(model.nodes.items()):
            #print "i = ",i
            point = node.Position()
            #print "point = ",point
            #sys.stdout.flush()
            #aVoxel = vtk.vtkPixel()
            #print "made voxel"; sys.stdout.flush()
            points.InsertPoint(i, *point)

            #print str(element)

            #elem = vtk.vtkVertex()
            #elem.GetPointIds().SetId(0, i)
            #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            nidMap[nid] = i
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
        self.mapElements(points,points2,nidMap,model,j)

    def mapElements(self,points,points2,nidMap,model,j):
        for eid,element in sorted(model.elements.items()):
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

        ###
        self.grid.SetPoints(points)
        self.grid2.SetPoints(points2)
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        print "updated grid"

    def OnKeyPress(self,obj,event):
        rwi = obj
        key = rwi.GetKeySym()
        print "*Pressed %s" %(key)

    def buildVTK(self,bdfFileName=None):
        self.bdfFileName = bdfFileName
        self.loadGeometry(self.bdfFileName)
        self.main()

        #self.renWin = vtk.vtkRenderWindow()
        #renWin = window
        
        cam = self.rend.GetActiveCamera()
        mouseArgs = {'pipeline':self,'camera':cam}

        xSize = 500
        ySize = 400
        (x,y) = getScreenCorner(xSize,ySize)
        #renWin.SetSize(xSize,ySize)
        #renWin.SetPosition(x,y)


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


    

#------------------------------------------------------------------------------

class AppFrame( wx.Frame ) :

    def __init__( self) :
        
        wx.Frame.__init__( self, None, -1, title='pyNastran' )
        self.bdfFileName = None
        self.dirname = ''

        #bdfFileName = 'box_cylindrical_coord_sys.bdf'

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        # Must call before any event handler is referenced.
        self.eventsHandler = EventsHandler( self )

        self.buildMenuBar()
        #self.buildToolBar()
        #self.buildToolBar2()
        self.buildStatusBar()
        self.frmPanel = Pan(self)
        
        self.SetMenuBar(self.menubar)

        self.frmPanel.bdfFileName = self.bdfFileName
        self.frmPanel.buildVTK(self.bdfFileName)

        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)


        #self.vbox.Add(self.frmPanel.widget, 0, wx.EXPAND)
        #self.frmPanel.BackgroundColour = (200, 240, 250)    # light blue
        #-----

        # Add them to sizer.
        hbox = wx.BoxSizer( wx.HORIZONTAL )
        hbox.Add( self.frmPanel.widget, 1, wx.EXPAND|wx.ALL, 1 )

        # Add buttons in their own sizer
        if 1:
            self.redBtn   = wx.Button( self.frmPanel, label='Red' )
            self.greenBtn = wx.Button( self.frmPanel, label='Green' )
            self.exitBtn  = wx.Button( self.frmPanel, label='Exit' )

            buttonSizer = wx.BoxSizer( wx.VERTICAL )
            buttonSizer.AddStretchSpacer()
            buttonSizer.Add( self.redBtn,   proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            buttonSizer.Add( self.greenBtn, proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            buttonSizer.Add( self.exitBtn,  proportion=0, flag=wx.EXPAND|wx.ALL, border=5 )
            buttonSizer.AddStretchSpacer()

            hbox.Add( buttonSizer, 0, wx.EXPAND| wx.ALL, 5 )

        # SetSizer both sizers in the most senior control that has sizers in it.
        self.vbox.AddStretchSpacer()
        #self.vbox.Add(hbox, 0, wx.EXPAND|wx.ALL, 5 )
        self.vbox.Add(hbox)
        #self.frmPanel.SetSizer( self.vbox )
        self.frmPanel.SetSizer( hbox )
        self.frmPanel.Layout()
        #self.toolbar1.Realize()


        #-----

        if 0:
            # Bind event handlers to all controls that have one.
            self.redBtn.  Bind( wx.EVT_BUTTON, self.eventsHandler.OnRedBtn )
            self.greenBtn.Bind( wx.EVT_BUTTON, self.eventsHandler.OnGreenBtn )
            self.exitBtn. Bind( wx.EVT_BUTTON, self.eventsHandler.OnExit )

            # Create more convenient ways to close this app.
            # Adding these makes a total of 5 separate ways to exit.
            self.frmPanel .Bind( wx.EVT_LEFT_DCLICK, self.eventsHandler.OnExit )
            self.colourPnl.Bind( wx.EVT_LEFT_DCLICK, self.eventsHandler.OnExit )
        ###
        
        events = self.eventsHandler
        # Bind Controls
        #self.Bind(wx.EVT_RIGHT_DOWN, events.OnRightDown)

        # Bind File Menu
        self.Bind(wx.EVT_MENU, events.OnExit, self.exitButton)
        
        # Bind View Menu
        self.Bind(wx.EVT_MENU, events.OnBackgroundColor, self.bkgColorView)
        self.Bind(wx.EVT_MENU, events.ToggleStatusBar, self.showStatusBar)
        self.Bind(wx.EVT_MENU, events.ToggleToolBar, self.showToolBar)
        
        # Bind Help Menu
        self.Bind(wx.EVT_MENU, events.OnAbout, id=ID_ABOUT)



    #end __init__

    def UpdateWindowName(self,bdfFileName):
        self.bdfFileName = bdfFileName
        self.frmPanel.bdfFileName = bdfFileName
        windowName = self.frmPanel.getWindowName()
        self.SetTitle(windowName)
    
    def buildStatusBar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText('Ready')

    def buildToolBar2(self):
        self.toolbar1 = wx.ToolBar(self)
        topen = self.toolbar1.AddLabelTool(wx.ID_OPEN, '', wx.Bitmap('icons/topen.png'))
        qtool = self.toolbar1.AddLabelTool(wx.ID_EXIT, '', wx.Bitmap('icons/texit.png'))
        self.vbox.Add(self.toolbar1)

    def buildToolBar(self):
        events = self.eventsHandler

        #toolbar1.AddSeparator()
        #toolbar1.AddSeparator()
        #tnew  = toolbar1.AddLabelTool(wx.ID_ANY,  '', wx.Bitmap('icons/new.png'))
        #tsave = toolbar1.AddLabelTool(ID_SAVEAS,  '', wx.Bitmap('icons/tsave.png'))
        #tundo = toolbar1.AddLabelTool(wx.ID_UNDO, '', wx.Bitmap('icons/tundo.png'))
        #tredo = toolbar1.AddLabelTool(wx.ID_REDO, '', wx.Bitmap('icons/tredo.png'))

        # toolbar at top - toggles
        toolbar1 = wx.ToolBar(self)
        topen = toolbar1.AddLabelTool(wx.ID_OPEN, '', wx.Bitmap('icons/topen.png'))
        qtool = toolbar1.AddLabelTool(wx.ID_EXIT, '', wx.Bitmap('icons/texit.png'))
        toolbar1.EnableTool(wx.ID_REDO, False)

        self.frmPanel.toolbar1 = toolbar1

        self.vbox.Add(self.toolbar1, 0, wx.EXPAND)

        self.Bind(wx.EVT_TOOL, events.OnLoadBDF,  id=wx.ID_OPEN)
        self.Bind(wx.EVT_TOOL, events.OnExit,     qtool)
        #self.Bind(wx.EVT_TOOL, events.OnSaveAsFile, id=ID_SAVEAS)
        #self.Bind(wx.EVT_TOOL, events.OnUndo, tundo)
        #self.Bind(wx.EVT_TOOL, events.OnRedo, tredo)

    def buildMenuBar(self):
        events = self.eventsHandler

        menubar = wx.MenuBar()
        # file menu
        fileMenu = wx.Menu()
        #fileMenu.Append(wx.ID_NEW,  '&New','does nothing')
        fileMenu.Append(wx.ID_OPEN, '&Load BDF','Loads a BDF')
        #fileMenu.Append(wx.ID_RES, 'Load OP2 &Results','Loads a OP2 - does nothing')
        #fileMenu.Append(wx.ID_SAVE, '&Save','does nothing')
        fileMenu.AppendSeparator()

        
        # dummy import submenu
        imp = wx.Menu()
        imp.Append(wx.ID_ANY, 'Import newsfeed list...')
        imp.Append(wx.ID_ANY, 'Import bookmarks...')
        imp.Append(wx.ID_ANY, 'Import mail...')

        #fileMenu.AppendMenu(wx.ID_ANY, 'I&mport', imp)

        self.exitButton = wx.MenuItem(fileMenu, wx.ID_EXIT, 'Exit','Exits pyNastran')
        self.exitButton.SetBitmap(wx.Image('icons/texit.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        fileMenu.AppendItem(self.exitButton)

        # view menu
        # status bar at bottom - toggles
        viewMenu = wx.Menu()
        self.bkgColorView  = viewMenu.Append(wx.ID_ANY, 'Change Background Color','Change Background Color')
        self.showStatusBar = viewMenu.Append(wx.ID_ANY, 'Show statusbar', 'Show Statusbar', kind=wx.ITEM_CHECK)
        self.showToolBar   = viewMenu.Append(wx.ID_ANY, 'Show toolbar',   'Show Toolbar',   kind=wx.ITEM_CHECK)
        viewMenu.Check(self.showStatusBar.GetId(), True)
        viewMenu.Check(self.showToolBar.GetId(), True)

        self.Bind(wx.EVT_TOOL, events.OnLoadBDF,  id=wx.ID_OPEN)
        # help/about menu
        helpMenu = wx.Menu()
        helpMenu.Append(ID_ABOUT, '&About', 'About pyNastran')

        # menu bar
        menubar.Append(fileMenu, '&File')
        menubar.Append(viewMenu, '&View')
        menubar.Append(helpMenu, '&Help')
        self.menubar = menubar



#end AppFrame class

#------------------------------------------------------------------------------

class EventsHandler(object) :

    def __init__(self,parent):
        self.parent = parent

    # File Menu
    def OnLoadBDF(self, event):
        """ Open a file"""
        #print "OnOpen..."

        if 0:
            if self.parent.bdfFileName=='test_tet10.bdf':
                bdfFileName = 'box_cylindrical_coord_sys.bdf'
            else:
                bdfFileName = 'test_tet10.bdf'
            self.parent.bdfFileName = bdfFileName
            #self.parent.frmPanel.buildVTK(bdfFileName)
            self.parent.frmPanel.loadGeometry(bdfFileName)
            #self.parent.frmPanel.getWindow.Render()
            self.parent.frmPanel.Update()
            return
        
        wildcard = "Nastran BDF (*.bdf; *.dat; *.nas)|*.bdf;*.dat;*.nas|" \
         "All files (*.*)|*.*"

        dlg = wx.FileDialog(None, "Choose a file", self.parent.dirname, "", wildcard, wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            bdfFileName         = dlg.GetFilename()
            self.parent.dirname = dlg.GetDirectory()
            fname = os.path.join(self.parent.dirname, bdfFileName)
            print "fname = ",fname
            self.parent.UpdateWindowName(bdfFileName)
            self.parent.frmPanel.loadGeometry(bdfFileName)
            self.parent.frmPanel.Update()
        dlg.Destroy()

    def OnExit(self,event):
        self.parent.Destroy()

    # View Menu
    def OnBackgroundColor(self,event):
        dlg = wx.ColourDialog(self.parent)
        dlg.GetColourData().SetChooseFull(True)
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            #print 'You selected: %s\n' % str(data.GetColour().Get())
            self.ShowColourAndDialog(data.GetColour().Get())
        dlg.Destroy()

    def ShowColourAndDialog(self,pnlColor):
        rend = self.parent.frmPanel.widget.GetCurrentRenderer()
        if not rend:
            rend = self.parent.frmPanel.widget.GetCurrentRenderer()
        
        color = [pnlColor[0]/255.,pnlColor[1]/255.,pnlColor[2]/255.]
        rend.SetBackground(color)
        self.parent.frmPanel.widget.Render()

    def ToggleStatusBar(self, e):
        if self.parent.showStatusBar.IsChecked():
            self.parent.statusbar.Show()
        else:
            self.parent.statusbar.Hide()

    def OnKeyPress(self,obj,event):
        rwi = obj
        key = rwi.GetKeySym()
        print "*Pressed %s" %(key)

    def ToggleToolBar(self, e):
        if self.parent.showToolBar.IsChecked():
            self.parent.toolbar1.Show()
        else:
            self.parent.toolbar1.Hide()

    # Help Menu
    def OnAbout(self, event):
        about = [
            'pyNastran v0.3.0',
            'Copyright Steven P. Doyle 2011-2012\n',
            'code.google.com/p/pynastran/',
            '',
            'Controls',
              'X/x - snap to x axis',
              'Y/y - snap to y axis',
              'Z/z - snap to z axis',
              '',
              'm/M scale up/scale down',
              'o/O rotate counter-clockwise/clockwise 5 degrees',
              'w   wireframe',
              's   surface',
              'i   take a screenshot (image)',]
        
        # not done
              #'',
              #'left arrow  - pan left  (not done)',
              #'right arrow - pan right (not done)',
              #'up arrow    - pan up    (not done)',
              #'down arrow  - pan down  (not done)',
              #'',
              #'p   project point (not done)',
              #'f   fly to rotation point (not done)',

        dlg = wx.MessageDialog(None, '\n'.join(about), 'About',
                 wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

#end Events class

#------------------------------------------------------------------------------

def Main():
    app = wx.App( redirect=False )
    appFrm = AppFrame()
    appFrm.Show()
    app.MainLoop()

#end class

#==============================================================================

if __name__ == '__main__' :

    Main()