#!/usr/bin/env python

# This example shows how to manually construct unstructured grids
# using Python.  Unstructured grids require explicit point and cell
# representations, so every point and cell must be created, and then
# added to the vtkUnstructuredGrid instance.

import wx
import vtk
import sys
import pyNastran
from pyNastran.bdf.bdf import BDF,CTRIA3,CQUAD4,CTETRA4,CPENTA6,CHEXA8,LineElement,CONM2,SpringElement
from mouseStyle import MouseStyle


def getScreenCorner(x,y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen,yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen-x)//2
    yCorner = (yScreen-y)//2
    return(xCorner,yCorner)

version = pyNastran.__version__
bdfFileName = sys.argv[1]


class Pipeline(object):
    def __init__(self):
        model = BDF()
        model.readBDF(bdfFileName)

        nNodes = model.nNodes()
        nElements = model.nElements()
        print "nNodes = ",nNodes
        print "nElements = ",nElements

        aQuadGrid = vtk.vtkUnstructuredGrid()
        #aQuadGrid.Allocate(nElements+nNodes, 1000)
        aQuadGrid.Allocate(nElements, 1000)


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
            #aVoxel.GetPointIds().SetId(i, i)
            points.InsertPoint(i, *point)

            #print str(element)

            #elem = vtk.vtkVertex()
            #elem.GetPointIds().SetId(0, i)
            #aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())


            nidMap[nid] = i
            i+=1

        if 0:
            for eid,element in sorted(model.caeros.items()):
                if isinstance(element,CAERO1):
                    nodeIDs = element.nodeIDs()
                    elem = vtk.vtkQuad()
                    elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                    elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                    elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                    aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())


        for eid,element in sorted(model.elements.items()):
            if isinstance(element,CTRIA3):
                #print "ctria3"
                elem = vtk.vtkTriangle()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CQUAD4):
                #print "cquad4"
                nodeIDs = element.nodeIDs()
                elem = vtk.vtkQuad()
                #print nodeIDs
                #print "n1=%s n2=%s n3=%s n4=%s" %(nidMap[nodeIDs[0]], nidMap[nodeIDs[1]], nidMap[nodeIDs[2]], nidMap[nodeIDs[3]])
                sys.stdout.flush()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA4):
                #print "ctetra"
                elem = vtk.vtkTetra()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
            elif isinstance(element,CPENTA6):
                #print "cpenta"
                elem = vtk.vtkWedge()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
            elif isinstance(element,CHEXA8):
                #print "chexa"
                elem = vtk.vtkHexahedron()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,LineElement) or isinstance(element,SpringElement):
                elem = vtk.vtkLine()
                #print str(element)
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            ###
            elif isinstance(element,CONM2):
                nid = element.Nid()
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()
                #elem.SetRadius(1.0)
                #print str(element)
                elem.GetPointIds().SetId(0, nidMap[nid])
                #elem.SetCenter(points.GetPoint(nidMap[nid]))
                aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                print "skipping %s" %(element.type)

        ###

        aQuadGrid.SetPoints(points)

        #Filter = vtk.vtkRotationFilter()
        #Filter.SetCenter(0.,0.,0.)
        #Filter.SetNumberOfCopies(1)
        #Filter.SetInput(aQuadGrid)
        #Filter.Update()


        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(aQuadGrid)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0) # red



        # Create the usual rendering stuff.
        self.rend    = vtk.vtkRenderer()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(.1, .2, .4)
        self.rend.AddActor(geometryActor)
        self.rend.ResetCamera()

        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer(self.rend)
        
        cam = self.rend.GetActiveCamera()
        mouseArgs = {'pipeline':self,'camera':cam}

        xSize = 500
        ySize = 400
        k = wx.App(False) # required to get screen resolution; will be True in final gui
        (x,y) = getScreenCorner(xSize,ySize)
        self.renWin.SetSize(xSize,ySize)
        self.renWin.SetPosition(x,y)


        iren = vtk.vtkRenderWindowInteractor()
        mouseStyle = MouseStyle(mouseArgs,iren)
        iren.SetInteractorStyle(mouseStyle)
        iren.SetRenderWindow(self.renWin)

        iren.AddObserver("KeyPressEvent", mouseStyle.OnKeyPress)

        print "type(ren) = ",type(self.rend)
        #self.rend.GetActiveCamera().Zoom(2.0)

        # Render the scene and start interaction.
        iren.Initialize()
        self.renWin.Render()
        self.renWin.SetWindowName("pyNastran v%s - %s" %(version,bdfFileName))
        iren.Start()

    def takePicture(self):
        """
        doesnt work...
        """
        #imageName = 'image.tiff'
        imageName = 'image.png'
        self.captureImage(imageName)
        print "took picture %s" %(imageName)

    def captureImage(self,imageName):
        """
        doesnt work...
        """
        w2i = vtk.vtkWindowToImageFilter()
        #writer = vtk.vtkTIFFWriter()
        writer = vtk.vtkPNGWriter()
        w2i.SetInput(self.renWin)
        w2i.Update()
        writer.SetInputConnection(w2i.GetOutputPort())
        self.renWin.Render()
        writer.SetFileName(imageName)
        writer.Write()

if __name__=='__main__':
    p = Pipeline()
