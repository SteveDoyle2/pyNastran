import numpy as np
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

def box_faces(n1, n2):
    """
    0----3
    |    |
    |    |
    1----2
    """
    w = 0.5
    h = 0.5

    points = np.array([
        [0., -w/2, h/2],   # 0
        [0., -w/2, -h/2],  # 1
        [0., w/2, -h/2],   # 2
        [0., w/2, h/2],    # 3

        [0., w/2, h/2],    # 3
        [0., w/2, -h/2],   # 2
        [0., -w/2, -h/2],  # 1
        [0., -w/2, h/2],   # 0
    ])  # 16 x 3

    xform = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ]) # 3x3 unit matrix

    points2 = points @ xform  # 16x3
    points2[:12, :] += np.array(n1)
    points2[12:, :] += np.array(n2)

    element = vtk.vtkHexahedron()
    return element, points2

def i_faces(n1, n2):
    """
         ^z
         |
    0----------11 - z3
    |          |
    1---2  9---10 - z2
        |  |
        |  |
        |  |
    4---3  8---7 - z1
    |          |
    5-----+----6-----> y
    """
    wflange_top = 0.5
    wflange_btm = 0.5
    hweb = 1.
    #ysc = 0.
    tflange_top = 0.2
    tflange_btm = 0.2
    tweb = 0.1
    z0 = 0.
    z1 = tflange_btm
    z2 = tflange_btm + hweb
    z3 = tflange_btm + hweb + tflange_top

    faces = [
        # around the I (counter-clockwise)
        [0, 12, 13, 1],
        [1, 13, 14, 2],
        [2, 14, 15, 3],
        [3, 15, 16, 4],
        [4, 16, 17, 5],
        [5, 17, 18, 6],
        [6, 18, 19, 7],
        [7, 19, 20, 8],
        [8, 20, 21, 9],
        [9, 21, 22, 10,],
        [10, 22, 23, 11,],
        [11, 23, 12, 0],

        # front face
        [0, 1, 2, 9, 10, 11],
        [2, 3, 8, 9],
        [4, 5, 6, 7, 8, 3],

        # back face
        [12, 13, 14, 21, 22, 23],
        [14, 15, 20, 21],
        [16, 17, 18, 19, 20, 15],
    ]

    points = np.array([
        [0., -wflange_top/2, z3],   # 0
        [0., -wflange_top/2, z2],   # 1
        [0., -tweb/2,        z2],   # 2
        [0., -tweb/2,        z1],   # 3
        [0., -wflange_btm/2, z1],   # 4
        [0., -wflange_btm/2, z0],   # 5

        [0., wflange_btm/2, z0],   # 6
        [0., wflange_btm/2, z1],   # 7
        [0., tweb/2,        z1],   # 8
        [0., tweb/2,        z2],   # 9
        [0., wflange_top/2, z2],   # 10
        [0., wflange_top/2, z3],   # 11

        [0., -wflange_top/2, z3],   # 0
        [0., -wflange_top/2, z2],   # 1
        [0., -tweb/2,        z2],   # 2
        [0., -tweb/2,        z1],   # 3
        [0., -wflange_btm/2, z1],   # 4
        [0., -wflange_btm/2, z0],   # 5
        [0., wflange_btm/2, z0],   # 6
        [0., wflange_btm/2, z1],   # 7
        [0., tweb/2,        z1],   # 8
        [0., tweb/2,        z2],   # 9
        [0., wflange_top/2, z2],   # 10
        [0., wflange_top/2, z3],   # 11
    ])  # 16 x 3

    xform = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ]) # 3x3 unit matrix

    points2 = points @ xform  # 16x3
    points2[:12, :] += np.array(n1)
    points2[12:, :] += np.array(n2)
    return faces, points2

def chan_faces(n1, n2):
    """
 ^z x0
 |  0-------7
 |  |       |
 |  |  5----6
 }  |  |
 +--|--|------> y
    |  |
    |  4----3
    |       |
    1-------2
    """
    faces = [
        # around the T (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],

        # front face
        [0, 1, 4, 5],
        [1, 2, 3, 4],
        [0, 5, 6, 7],

        # back face
        [8, 9, 12, 13],
        [9, 10, 11, 12],
        [8, 13, 14, 15],

        # alternate method for front/back face
        #[0, 1, 2,  3,  4,  5,  6,  7],
        #[8, 9, 10, 11, 12, 13, 14, 15],
    ]

    wflange = 0.5
    hweb = 1.
    #ysc = 0.
    tflange = 0.2
    tweb = 0.1
    points = np.array([
        [0.,  hweb/2, 0.],                  # 0
        [0., -hweb/2, 0.],                  # 1
        [0., -hweb/2, wflange],             # 2
        [0., -hweb/2 + tflange, wflange],   # 3
        [0., -hweb/2 + tflange, tweb],      # 4
        [0.,  hweb/2 - tflange, tweb],      # 5
        [0.,  hweb/2 - tflange, wflange],   # 6
        [0.,  hweb/2,           wflange],   # 7

        [0.,  hweb/2, 0.],                  # 0
        [0., -hweb/2, 0.],                  # 1
        [0., -hweb/2, wflange],             # 2
        [0., -hweb/2 + tflange, wflange],   # 3
        [0., -hweb/2 + tflange, tweb],      # 4
        [0.,  hweb/2 - tflange, tweb],      # 5
        [0.,  hweb/2 - tflange, wflange],   # 6
        [0.,  hweb/2,           wflange],   # 7
    ])  # 16 x 3

    xform = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ]) # 3x3 unit matrix

    points2 = points @ xform  # 16x3
    points2[:8, :] += np.array(n1)
    points2[8:, :] += np.array(n2)
    return faces, points2

def t_faces(n1, n2):
    """
         ^ z
         |
    0----------7
    |    +     |----> y
    1---2  5---6
        |  |
        }  |
        3--4

    8----------15
    |    +     |
    9---10 13--14
        |  |
        }  |
        11-12
    """
    faces = [
        # around the T (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],

        # front face
        [0, 1, 2, 5, 6, 7],
        [2, 3, 4, 5],

    # back face
    [8, 9, 10, 13, 14, 15],
    [10, 11, 12, 13],

        # front/back face
        #[0, 1, 2,  3,  4,  5,  6,  7],
        #[8, 9, 10, 11, 12, 13, 14, 15],
    ]

    wflange = 0.5
    wweb = 0.1
    hweb = 1.
    t = 0.1
    points = np.array([
        [0.,-wflange/2,  t/2], # 0
        [0.,-wflange/2, -t/2], # 1
        [0.,-wweb/2, -t/2], # 2
        [0.,-wweb/2, -hweb + t/2], # 3
        [0., wweb/2, -hweb + t/2], # 4
        [0., wweb/2, -t/2], # 5
        [0., wflange/2, -t/2], # 6
        [0., wflange/2,  t/2], # 7

        [0.,-wflange/2,  t/2], # 0
        [0.,-wflange/2, -t/2], # 1
        [0.,-wweb/2, -t/2], # 2
        [0.,-wweb/2, -hweb + t/2], # 3
        [0., wweb/2, -hweb + t/2], # 4
        [0., wweb/2, -t/2], # 5
        [0., wflange/2, -t/2], # 6
        [0., wflange/2,  t/2], # 7
    ])  # 16 x 3
    xform = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ]) # 3x3 unit matrix

    points2 = points @ xform  # 16x3
    points2[:8, :] += np.array(n1)
    points2[8:, :] += np.array(n2)
    return faces, points2

def t2_faces(n1, n2):
    """
         ^ z
         |
    0----------7
    |    +     |----> y
    1---2  5---6
        |  |
        }  |
        3--4

    8----------15
    |    +     |
    9---10 13--14
        |  |
        }  |
        11-12
    """
    wflange = 0.5
    wweb = 0.1
    hweb = 1.
    t = 0.1
    #y = [0., 1., 0]
    #z = [0., 1., 0]

    points = [
        [n1[0], -wflange/2,  t/2], # 0
        [n1[0], -wflange/2, -t/2], # 1
        [n1[0], -wweb/2, -t/2], # 2
        [n1[0], -wweb/2, -hweb + t/2], # 3
        [n1[0],  wweb/2, -hweb + t/2], # 4
        [n1[0], wweb/2, -t/2], # 5
        [n1[0], wflange/2, -t/2], # 6
        [n1[0], wflange/2,  t/2], # 7

        [n2[0], -wflange/2,  t/2], # 0
        [n2[0], -wflange/2, -t/2], # 1
        [n2[0], -wweb/2, -t/2], # 2
        [n2[0], -wweb/2, -hweb + t/2], # 3
        [n2[0],  wweb/2, -hweb + t/2], # 4
        [n2[0],  wweb/2, -t/2], # 5
        [n2[0],  wflange/2, -t/2], # 6
        [n2[0],  wflange/2,  t/2], # 7
    ]
    faces = [
        # around the T (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],

        # front/back face
        [0, 1, 2,  3,  4,  5,  6,  7,  0],
        [8, 9, 10, 11, 12, 13, 14, 15, 8],
    ]
    return faces, points

def hollow_box_faces(n1, n2):
    """
         ^ y
         |
    0---------3
    |         |
    |  4---7  |
    |  | + |  |---> z
    |  5---6  |
    |         |
    1---------2

    8---------11
    |         |
    |  12--15 |
    |  | + |  |
    |  13--14 |
    |         |
    9---------10
    """
    hbox = 1.0
    wbox = 2.0
    tw = 0.1
    th = 0.1

    hbox_in = hbox - 2*th
    wbox_in = wbox - 2*tw
    faces = [
        # front face
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],

        # back face
        [8, 9, 13, 12],
        [9, 10, 14, 13],
        [10, 11, 15, 14],
        [11, 8, 12, 15],

        # front/back face
       #[0, 1,  2,  3, 0,  4,  7,  6,  5,  4, 0],
       #[8, 9, 10, 11, 8, 12, 15, 14, 13, 12, 8],

       #[0, 1,  2,  3,  4,  5,  6,  7,  4],
       #[8, 9, 10, 11, 12, 13, 14, 15, 12]

        # around the box (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 8, 0],

        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 12, 4],
    ]
    points = [
        [n1[0], -wbox/2., hbox/2], # 0
        [n1[0], -wbox/2, -hbox/2], # 1
        [n1[0],  wbox/2, -hbox/2], # 2
        [n1[0],  wbox/2,  hbox/2], # 3

        [n1[0], -wbox_in/2., hbox_in/2], # 4
        [n1[0], -wbox_in/2, -hbox_in/2], # 5
        [n1[0],  wbox_in/2, -hbox_in/2], # 6
        [n1[0],  wbox_in/2,  hbox_in/2], # 7

        # faces2
        [n2[0], -wbox/2., hbox/2], # 8
        [n2[0], -wbox/2, -hbox/2], # 9
        [n2[0],  wbox/2, -hbox/2], # 10
        [n2[0],  wbox/2,  hbox/2], # 11

        [n2[0], -wbox_in/2., hbox_in/2], # 12
        [n2[0], -wbox_in/2, -hbox_in/2], # 13
        [n2[0],  wbox_in/2, -hbox_in/2], # 14
        [n2[0],  wbox_in/2,  hbox_in/2], # 15
    ]
    return faces, points

import vtk
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk, numpy_to_vtkIdTypeArray
from pyNastran.gui.utils.vtk.base_utils import VTK_VERSION


# Create polyhedron (cube)



# List of pointIds that make up the face

n1 = [0., 0., 0.]
n2 = [5., 0., 0.]
faceList, nodes = t_faces(n1, n2)
faceList, nodes = hollow_box_faces(n1, n2)
#faceList, nodes = chan_faces(n1, n2)
faceList, nodes = i_faces(n1, n2)

nnodes = len(nodes)
pointIds = vtk.vtkIdList()
pointIds.SetNumberOfIds(nnodes)

for i in range(nnodes):
    pointIds.InsertNextId(i)
    points = vtk.vtkPoints()

for node in nodes:
    #print(node)
    points.InsertNextPoint(*node)

#faceList = [
#    [0, 3, 2, 1],
#    [0, 4, 7, 3],
#    [4, 5, 6, 7],
#    [5, 1, 2, 6],
#    [0, 1, 5, 4],
#    [2, 3, 7, 6],
#]

def faces_to_element_facelist(faces):
    face_idlist = vtk.vtkIdList()

    nfaces = len(faces)
    face_idlist.InsertNextId(nfaces) # Number faces that make up the cell.
    for face in faces: # Loop over all the faces
        face_idlist.InsertNextId(len(face)) # Number of points in face
        for i in face:
            face_idlist.InsertNextId(i)
        #[faceId.InsertNextId(i) for i in face] # Insert the pointIds for the face

    return face_idlist

face_idlist = faces_to_element_facelist(faceList)
ugrid = vtkUnstructuredGrid()
ugrid.SetPoints(points)
ugrid.InsertNextCell(vtk.VTK_POLYHEDRON, faceId)

#-------------------------------------
if 1:
    ug = ugrid
    grid_mapper = vtk.vtkDataSetMapper()
    vtk_version = int(VTK_VERSION[0])
    grid_mapper.SetInputData(ug)

#-------------------------------------

#Create a mapper and actor
#mapper = vtk.vtkPolyDataMapper()
#mapper.SetInputConnection(text_source.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(grid_mapper)

#Create a renderer, render window, and interactor
renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

#Add the actor to the scene
renderer.AddActor(actor)
renderer.SetBackground(0, 0, 0) # Background color white

#Render and interact
render_window.Render()
render_window_interactor.Start()
