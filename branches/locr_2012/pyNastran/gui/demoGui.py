#!/usr/bin/env python

# This example shows how to manually construct unstructured grids
# using Python.  Unstructured grids require explicit point and cell
# representations, so every point and cell must be created, and then
# added to the vtkUnstructuredGrid instance.

import wx
import vtk
from vtk import (vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
                 vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra, vtkQuadraticWedge, vtkQuadraticHexahedron)


import sys
import pyNastran
from pyNastran.bdf.bdf import (BDF, CAERO1, CQUAD4, CQUAD8, CTRIA3, CTRIA6,
                               CTETRA4, CTETRA10, CPENTA6, CPENTA15,
                               CHEXA8, CHEXA20, CONM2, LineElement,
                               SpringElement)
from mouseStyle import MouseStyle


def getScreenCorner(x, y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen, yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen - x) // 2
    yCorner = (yScreen - y) // 2
    return(xCorner, yCorner)

version = pyNastran.__version__
bdfFileName = sys.argv[1]


class FrameVTK(object):
    def __init__(self, isEdges):
        self.isEdges = isEdges  # surface wireframe
        self.loadGeometry(bdfFileName)
        self.main()

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
        self.edgeActor.GetProperty().SetColor(0, 0, 0)
        prop = self.edgeActor.GetProperty()
        #prop.SetLineWidth(0.0)
        if self.isEdges:
            prop.EdgeVisibilityOn()
        else:
            prop.EdgeVisibilityOff()

        self.rend.AddActor(self.edgeActor)

        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        print "visible = ", prop.GetEdgeVisibility()
        #self.edgeActor.Update()

    def addGeometry(self):
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(geometryActor)

    def addAltGeometry(self):
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

    def main(self):
        self.rend = vtk.vtkRenderer()
        self.addGeometry()
        self.addAltGeometry()

        # Create the usual rendering stuff.
        if self.isEdges:
            self.getEdges()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(.1, .2, .4)
        self.rend.ResetCamera()

        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer(self.rend)

        cam = self.rend.GetActiveCamera()
        mouseArgs = {'pipeline': self, 'camera': cam}

        xSize = 500
        ySize = 400
        k = wx.App(False)  # required to get screen resolution; will be True in final gui
        (x, y) = getScreenCorner(xSize, ySize)
        self.renWin.SetSize(xSize, ySize)
        self.renWin.SetPosition(x, y)

        iren = vtk.vtkRenderWindowInteractor()
        mouseStyle = MouseStyle(mouseArgs, iren)
        iren.SetInteractorStyle(mouseStyle)
        iren.SetRenderWindow(self.renWin)

        iren.AddObserver("KeyPressEvent", mouseStyle.OnKeyPress)

        print "type(ren) = ", type(self.rend)
        #self.rend.GetActiveCamera().Zoom(2.0)

        # Render the scene and start interaction.
        iren.Initialize()
        self.renWin.Render()
        self.renWin.SetWindowName(
            "pyNastran v%s - %s" % (version, bdfFileName))
        iren.Start()

    def loadGeometry(self, bdf_filename):
        model = BDF()
        model.readBDF(bdf_filename)

        nNodes = model.nNodes()
        nElements = model.nElements()
        nCAeros = model.nCAeros()

        #print "nNodes = ",nNodes
        print "nElements = ", nElements

        self.grid = vtk.vtkUnstructuredGrid()
        self.grid2 = vtk.vtkUnstructuredGrid()
        #self.aQuadGrid.Allocate(nElements+nNodes, 1000)

        if 'CONM2' in model.cardCount:
            nCONM2 = model.cardCount['CONM2']
        else:
            nCONM2 = 0
        self.grid.Allocate(nElements, 1000)
        self.grid2.Allocate(nCAeros + nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nNodes)
        nidMap = {}
        i = 0
        #elem.SetNumberOfPoints(nNodes)
        for nid, node in sorted(model.nodes.iteritems()):
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
            i += 1

        j = 0
        points2 = vtk.vtkPoints()
        points2.SetNumberOfPoints(nCAeros * 4 + nCONM2)
        for eid, element in sorted(model.caeros.iteritems()):
            if isinstance(element, CAERO1):
                cpoints = element.Points()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j + 1)
                elem.GetPointIds().SetId(2, j + 2)
                elem.GetPointIds().SetId(3, j + 3)
                points2.InsertPoint(j, *cpoints[0])
                points2.InsertPoint(j + 1, *cpoints[1])
                points2.InsertPoint(j + 2, *cpoints[2])
                points2.InsertPoint(j + 3, *cpoints[3])
                self.grid2.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
                j += 4
        self.mapElements(points, points2, nidMap, model, j)

    def mapElements(self, points, points2, nidMap, model, j):
        for eid, element in sorted(model.elements.iteritems()):
            if isinstance(element, CTRIA3):
                #print "ctria3"
                elem = vtkTriangle()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTRIA6):
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
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CQUAD4):
                nodeIDs = element.nodeIDs()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CQUAD8):
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
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTETRA4):
                elem = vtkTetra()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CTETRA10):
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
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CPENTA6):
                elem = vtkWedge()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element, CPENTA15):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
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
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CHEXA8):
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
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, CHEXA20):
                nodeIDs = element.nodeIDs()
                #print "nodeIDs = ",nodeIDs
                if None not in nodeIDs:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
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
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element, LineElement) or isinstance(element, SpringElement):
                elem = vtk.vtkLine()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                self.grid.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
            ###
            elif isinstance(element, CONM2):  # not perfectly located
                nid = element.Nid()
                c = element.Centroid()
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()
                #elem.SetRadius(1.0)
                #print str(element)

                points2.InsertPoint(j, *c)
                elem.GetPointIds().SetId(0, j)
                #elem.SetCenter(points.GetPoint(nidMap[nid]))
                self.grid2.InsertNextCell(
                    elem.GetCellType(), elem.GetPointIds())
                j += 1
            else:
                print "skipping %s" % (element.type)

        ###
        self.grid.SetPoints(points)
        self.grid2.SetPoints(points2)
        self.grid.Update()
        self.grid2.Update()

    def takePicture(self):
        """
        doesnt work...
        """
        #imageName = 'image.tiff'
        imageName = 'image.png'
        self.captureImage(imageName)
        print "took picture %s" % (imageName)

    def captureImage(self, imageName):
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


def main():
    isEdges = False
    p = FrameVTK(isEdges)

if __name__ == '__main__':
    main()
