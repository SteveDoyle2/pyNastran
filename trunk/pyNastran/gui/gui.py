#!/usr/bin/env python

# This example shows how to manually construct unstructured grids
# using Python.  Unstructured grids require explicit point and cell
# representations, so every point and cell must be created, and then
# added to the vtkUnstructuredGrid instance.

import vtk
import sys
import pyNastran
from pyNastran.bdf.bdf import BDF,CTRIA3,CQUAD4,CTETRA4,CPENTA6,CHEXA8,LineElement,CONM2,SpringElement
from mouseStyle import MouseStyle

version = pyNastran.__version__
bdfFileName = sys.argv[1]


model = BDF()
model.readBDF(bdfFileName)

nNodes = model.nNodes()
nElements = model.nElements()
print "nNodes = ",nNodes
print "nElements = ",nElements

aQuadGrid = vtk.vtkUnstructuredGrid()
#aQuadGrid.Allocate(nElements+nNodes, 1000)
aQuadGrid.Allocate(nElements, 1000)

#voxelPoints.InsertPoint(0, 0, 0, 0)
#voxelPoints.InsertPoint(1, 1, 0, 0)
#voxelPoints.InsertPoint(2, 0, 1, 0)
#voxelPoints.InsertPoint(3, 1, 1, 0)
#voxelPoints.InsertPoint(4, 0, 0, 1)
#voxelPoints.InsertPoint(5, 1, 0, 1)
#voxelPoints.InsertPoint(6, 0, 1, 1)
#voxelPoints.InsertPoint(7, 1, 1, 1)
#aVoxel.GetPointIds().SetId(1, 1)
#aVoxel.GetPointIds().SetId(2, 2)
#aVoxel.GetPointIds().SetId(3, 3)
#aVoxel.GetPointIds().SetId(4, 4)
#aVoxel.GetPointIds().SetId(5, 5)
#aVoxel.GetPointIds().SetId(6, 6)
#aVoxel.GetPointIds().SetId(7, 7)
#aVoxelGrid = vtk.vtkUnstructuredGrid()
#aVoxelGrid.Allocate(1, 1)
#aVoxelGrid.InsertNextCell(aVoxel.GetCellType(), aVoxel.GetPointIds())
#aVoxelGrid.SetPoints(voxelPoints)



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
ren    = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

xSize = 500
ySize = 400
renWin.SetSize(xSize,ySize)
x = 350
y = 150
renWin.SetPosition(x,y)




iren = vtk.vtkRenderWindowInteractor()
style = MouseStyle(iren)
iren.SetInteractorStyle(style)
iren.SetRenderWindow(renWin)

ren.SetBackground(.1, .2, .4)
ren.AddActor(geometryActor)
ren.ResetCamera()
ren.GetActiveCamera().ParallelProjectionOn()

# Render the scene and start interaction.
iren.Initialize()
renWin.Render()
renWin.SetWindowName("pyNastran v%s - %s" %(version,bdfFileName))
iren.Start()