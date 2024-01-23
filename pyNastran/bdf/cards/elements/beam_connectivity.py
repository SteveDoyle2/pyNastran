# pylint: disable=C0326
"""
Defines the faces for representing the PBARL/PBEAML cards as 3D objects
instead of just 1D elements.

These functions intentionally do not handle offsets.
Node 1/2 should be the offset/actual coordinates.

The following shapes are fully supported:
 - ROD
 - TUBE
 - BAR
 - BOX
 - L
 - H
 - I1  (I1 is handled correctly because it's doubly symmetric)
 - T
 - T1
 - T2
 - Z
 - HEXA

The following shapes don't do the shear center correctly:
 - CHAN
 - CHAN1
 - I
 - HAT

The following shapes aren't supported:
 - CROSS
 - TUBE2
 - CHAN2 (will be a shear center issue)
 - BOX1  (will be a shear center issue)
 - HAT1  (will be a shear center issue)
 - DBOX  (will be a shear center issue)

"""
from typing import Union
import numpy as np
from pyNastran.nptyping_interface import NDArray3float, NDArray33float, NDArrayN3float
from pyNastran.bdf.cards.aero.utils import elements_from_quad, tri_cap

Faces = list[list[int]]
Dim1 = Union[np.ndarray, tuple[float]]
Dim2 = Union[np.ndarray, tuple[float, float]]
Dim3 = Union[np.ndarray, tuple[float, float, float]]
Dim4 = Union[np.ndarray, tuple[float, float, float, float]]
Dim6 = Union[np.ndarray, tuple[float, float, float, float, float, float]]

def _transform_points(n1: NDArray3float,
                      n2: NDArray3float,
                      points1: NDArrayN3float,
                      points2: NDArrayN3float,
                      xform: NDArray33float) -> NDArrayN3float:
    points_list = []
    for nid, points in [(n1, points1), (n2, points2)]:
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return np.vstack(points_list)

def rod_setup(dim1: Dim1, dim2: Dim1,
              ) -> tuple[Faces, np.ndarray, np.ndarray]: # validated
    """defines points in a circle with triangle based end caps"""
    # 4,8,12,16,... becomes 5,9,13,17,...
    ntheta = 17
    thetas = np.radians(np.linspace(0., 360., ntheta))

    nfaces = 0
    all_faces = []
    points_list = []
    x = np.zeros(ntheta)
    for dim in (dim1, dim2):
        radius, = dim
        y = radius * np.cos(thetas)
        z = radius * np.sin(thetas)
        xyz = np.vstack([x, y, z]).T
        assert xyz.shape == (ntheta, 3), xyz.shape

        # the tri_cap buils triangles that fan out from the first node
        tris = tri_cap(ntheta)

        # we need to use the tolist because we're going to
        # combine quads and tris (the elements have different
        # lengths)
        all_faces += (nfaces + tris).tolist()
        nfaces += tris.shape[0]
        points_list.append(xyz)

    # the main cylinder uses the points defined independent
    # of the points n1/n2
    faces = elements_from_quad(2, ntheta)
    all_faces += faces.tolist()
    return all_faces, points_list[0], points_list[1]

def rod_faces(n1: NDArray3float,
              n2: NDArray3float,
              xform: NDArray33float,
              dim1: Dim1,
              dim2: Dim1,
              ) -> tuple[Faces, np.ndarray, int]: # validated
    """defines points in a circle with triangle based end caps"""
    faces, points1, points2 = rod_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array, points_array.shape[0]

def tube_setup(dim1: Dim2, dim2: Dim2,
               ) -> tuple[np.ndarray, NDArrayN3float, NDArrayN3float]:  # validated
    """defines a rod with a hole"""
    # 4,8,12,16,... becomes 5,9,13,17,...
    ntheta = 17
    thetas = np.radians(np.linspace(0., 360., ntheta))
    npoints = ntheta * 2

    x = np.zeros(ntheta)

    points = []
    #points_list_in = []
    #points_list_out = []
    for dim in (dim1, dim2):
        radius_out, radius_in = dim

        # inner rod
        y = radius_in * np.cos(thetas)
        z = radius_in * np.sin(thetas)
        xyz2 = np.vstack([x, y, z]).T
        #points_list_in.append(xyz2)

        # outer rod
        y = radius_out * np.cos(thetas)
        z = radius_out * np.sin(thetas)
        xyz1 = np.vstack([x, y, z]).T
        #points_list_out.append(xyz1)
        xyz = np.vstack([xyz1, xyz2])
        points.append(xyz)


    # the main cylinder uses the points defined independent
    # of the inner/outer faces
    faces_n1 = elements_from_quad(2, ntheta)
    faces_n2 = faces_n1 + npoints

    #faces_n1
    #[[ 0 17 18  1]
     #[ 1 18 19  2]
     #[ 2 19 20  3]
     #[ 3 20 21  4]
     #[ 4 21 22  5]
     #[ 5 22 23  6]

    # we use slicing to link the outer surface to the
    # inner surface using quads
    #
    # we'd like to hstack the column arrays, but hstack and a
    # transpose is easier.  Presumably, a "newaxis" index would work...
    faces_out = np.vstack([
        faces_n1[:, 0],
        faces_n1[:, 3],
        faces_n2[:, 3],
        faces_n2[:, 0]]).T
    faces_in = np.vstack([
        faces_n1[:, 1],
        faces_n1[:, 2],
        faces_n2[:, 2],
        faces_n2[:, 1]]).T

    # combine everything together
    all_faces = np.vstack([faces_n1, faces_n2, faces_out, faces_in])
    return all_faces, points[0], points[1]

def tube_faces(n1: NDArray3float,
               n2: NDArray3float,
               xform: NDArray33float,
               dim1: Dim2,
               dim2: Dim2) -> tuple[np.ndarray, np.ndarray, int]:  # validated
    """defines a rod with a hole"""
    faces, points1, points2 = tube_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array, points_array.shape[0]

def bar_setup(dim1: Dim2, dim2: Dim2,
              ) -> tuple[Faces, NDArrayN3float, NDArrayN3float]:  # validated
    """
       ^y
       |
    0-----3
    |     |
    |     |----> z
    |     |
    1-----2

    """
    faces = [
        # front/back
        [0, 1, 2, 3],
        [4, 5, 6, 7],

        # 4-sides
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    points = []
    for dim in (dim1, dim2):
        w, h = dim
        pointsi = np.array([
            [0.,  h/2, -w/2],   # 0
            [0., -h/2, -w/2],   # 1
            [0., -h/2,  w/2],   # 2
            [0.,  h/2,  w/2],   # 3
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def bar_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
              dim1: Dim2, dim2: Dim2) -> np.ndarray:  # validated
    """builds the BAR faces"""
    unused_faces, points1, points2 = bar_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return points_array

def box_setup(dim1: Dim4, dim2: Dim4,
              ) -> tuple[Faces, NDArrayN3float, NDArrayN3float]:
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

    """
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

        # outer
        # around the box (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 8, 0],

        # inner
        # backwards - I don't think this matters
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 12, 4],
    ]
    points = []
    for dim in (dim1, dim2):
        wbox, hbox, th, tw = dim
        hbox_in = hbox - 2*th
        wbox_in = wbox - 2*tw
        pointsi = np.array([
            [0., -hbox/2., wbox/2], # 0
            [0., -hbox/2, -wbox/2], # 1
            [0.,  hbox/2, -wbox/2], # 2
            [0.,  hbox/2,  wbox/2], # 3

            [0., -hbox_in/2., wbox_in/2], # 4
            [0., -hbox_in/2, -wbox_in/2], # 5
            [0.,  hbox_in/2, -wbox_in/2], # 6
            [0.,  hbox_in/2,  wbox_in/2], # 7
        ])
        points.append(pointsi)
    return faces, points[0], points[1]

def box_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
              dim1: Dim4,
              dim2: Dim4) -> tuple[Faces, np.ndarray]:  # validated
    """builds the BOX faces"""
    faces, points1, points2 = box_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def i_setup(dim1: Dim6, dim2: Dim6):   # validated
    """
         ^y
         |
    0----------11 - y3
    |          |
    1---2  9---10 - y2
        |  |
        |  |
        |  |
    4---3  8---7 - y1
    |          |
    5-----+----6-----> z
    """
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
    points = []
    for dim in (dim1, dim2):
        hall, bflange_btm, bflange_top, tweb, tflange_btm, tflange_top = dim
        hweb = hall - tflange_top - tflange_btm
        ysc = -hall / 2.  # TODO: fix the shear center formula

        y0 = ysc
        y1 = y0 + tflange_btm
        y2 = y1 + hweb
        y3 = y2 + tflange_top
        pointsi = np.array([
            [0., y3, -bflange_top/2],   # 0
            [0., y2, -bflange_top/2],   # 1
            [0., y2, -tweb/2,      ],   # 2
            [0., y1, -tweb/2,      ],   # 3
            [0., y1, -bflange_btm/2],   # 4
            [0., y0, -bflange_btm/2],   # 5

            [0., y0, bflange_btm/2],   # 6
            [0., y1, bflange_btm/2],   # 7
            [0., y1, tweb/2,      ],   # 8
            [0., y2, tweb/2,      ],   # 9
            [0., y2, bflange_top/2],   # 10
            [0., y3, bflange_top/2],   # 11
        ])  # 24 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def i_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
            dim1: Dim6, dim2: Dim6):  # validated
    """builds the I faces"""
    faces, points1, points2 = i_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def i1_setup(dim1: Dim4,
             dim2: Dim4) -> tuple[Faces,
                                  NDArrayN3float, NDArrayN3float]:
    """builds the I1 faces"""
    """
         ^y
         |
    0----------11 - y3
    |          |
    1---2  9---10 - y2
        |  |
        |  |
        |  |
    4---3  8---7 - y1
    |          |
    5-----+----6-----> z
    """
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

    points = []
    for dim in (dim1, dim2):
        bfoot_2x, tweb, hin, hall = dim
        bflange = bfoot_2x + tweb
        bflange_top = bflange
        bflange_btm = bflange
        tflange_top = (hall - hin) / 2.
        tflange_btm = tflange_top

        hweb = hall - tflange_top - tflange_btm
        ysc = -hall / 2.  # TODO: fix the shear center formula

        y0 = ysc
        y1 = y0 + tflange_btm
        y2 = y1 + hweb
        y3 = y2 + tflange_top
        pointsi = np.array([
            [0., y3, -bflange_top/2],   # 0
            [0., y2, -bflange_top/2],   # 1
            [0., y2, -tweb/2,      ],   # 2
            [0., y1, -tweb/2,      ],   # 3
            [0., y1, -bflange_btm/2],   # 4
            [0., y0, -bflange_btm/2],   # 5

            [0., y0, bflange_btm/2],   # 6
            [0., y1, bflange_btm/2],   # 7
            [0., y1, tweb/2,      ],   # 8
            [0., y2, tweb/2,      ],   # 9
            [0., y2, bflange_top/2],   # 10
            [0., y3, bflange_top/2],   # 11
        ])  # 24 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def i1_faces(n1: NDArray3float,
             n2: NDArray3float,
             xform: NDArray33float,
             dim1: Dim4,
             dim2: Dim4):
    """builds the I1 faces"""
    """
         ^y
         |
    0----------11 - y3
    |          |
    1---2  9---10 - y2
        |  |
        |  |
        |  |
    4---3  8---7 - y1
    |          |
    5-----+----6-----> z
    """
    faces, points1, points2 = i1_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def h_setup(dim1: Dim4,
            dim2: Dim4):
    """
            ^y
            |
    0----11 |  8----7
    |    |  |  |    |
    |    10----9    |
    |    |  |  |    |
    |    |  +--|----|---> z
    |    |     |    |
    |    3-----4    |
    |    |     |    |
    |    |     |    |
    1----2     5----6

          <---> dim1
    <----> 0.5*dim2
    <---------------> w_all
    """
    faces = [
        # front face
        [0, 1, 2, 3, 10, 11],
        [3, 4, 9, 10],
        [4, 5, 6, 7, 8, 9],

        # back face
        [12, 13, 14, 15, 22, 23],
        [15, 16, 21, 22],
        [16, 17, 18, 19, 20, 21],

        # around the H (counter-clockwise)
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
    ]
    points = []
    for dim in (dim1, dim2):
        winner, wouter, hall, hinner = dim
        w_all = winner + wouter

        pointsi = np.array([
            [0.,  hall/2, -w_all/2],        # 0 - upper far left
            [0., -hall/2, -w_all/2],        # 1 - lower far left
            [0., -hall/2,   -winner/2],     # 2 - lower near left
            [0., -hinner/2, -winner/2],     # 3 - lower H left
            [0., -hinner/2,  winner/2],     # 4 - lower H right
            [0., -hall/2,    winner/2],     # 5 - lower near right
            [0., -hall/2,  w_all/2],        # 6 - lower far right

            [0., hall/2,  w_all/2],        # 7  - upper far right
            [0., hall/2,    winner/2],     # 8  - upper near right
            [0., hinner/2,  winner/2],     # 9  - upper H right
            [0., hinner/2, -winner/2],     # 10 - upper H left
            [0., hall/2,   -winner/2],     # 11 - upper near left
        ])  # 24 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def h_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
               dim1: Dim4,
               dim2: Dim4):  # validated
    """builds the H faces"""
    faces, points1, points2 = h_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def chan_setup(dim1: Dim4,
               dim2: Dim4):
    """
    ^y
    |  0--------7
    |  |        |
    |  |  5-----6
    |  |  |
    +--|--|-------> z
       |  |
       |  4-----3
       |        |
       1--------2
    <----> e/zsc

       <--------> bflange
    """
    faces = [
        # front face
        [0, 5, 6, 7],
        [0, 1, 4, 5],
        [1, 2, 3, 4],

        # back face
        [8, 13, 14, 15],
        [8, 9, 12, 13],
        [9, 10, 11, 12],

        # around the C (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],
    ]
    points = []
    for dim in (dim1, dim2):
        bflange, hall, tweb, tflange = dim
        # distance from shear center to neutral axis
        #zsc_na = 3 * bflange ** 2 / (6 * bflange + h) # per msc 2018 refman

        zsc = 0. ## TODO: consider the shear center
        #if 0:  # pragma: no cover
            # msc 2018 refman; p.670
            #h = hall - tflange
            #tw = tweb
            #tf = tflange
            #b = bflange - tw / 2.
            #bf = bflange - tw
            #hw = hall - 2. * tf
            #A = 2 * tf * bf + (h + tf) * tw
            #zc = bf * tf * (bf + tw) / A
            #zsc = b**2 * tf / (2*b*tw + 1/3. * h * tf)
            #E = zs - tw / 2.
            #zna = zc + zsc

        pointsi = np.array([
            [0., hall/2,  zsc], # 0
            [0., -hall/2, zsc], # 1

            [0., -hall/2,           zsc + bflange], # 2
            [0., -hall/2 + tflange, zsc + bflange], # 3
            [0., -hall/2 + tflange, zsc + tweb], # 4

            [0.,  hall/2 - tflange, zsc + tweb], # 5
            [0.,  hall/2 - tflange, zsc + bflange], # 6
            [0.,  hall/2,           zsc + bflange], # 7
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def chan_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
               dim1: Dim4,
               dim2: Dim4):  # validated
    """builds the CHAN faces"""
    faces, points1, points2 = chan_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def chan1_setup(dim1: Dim4,
                dim2: Dim4) -> tuple[Faces,
                                     NDArrayN3float, NDArrayN3float]:
    """
    ^y
    |  0--------7      ^ hall
    |  |        |      |
    |  |  5-----6      |  ^
    |  |  |            |  | hweb
    +--|--|-------> z  |  |
       |  |            |  v
       |  4-----3      |  ^
       |        |      |  | tflange
       1--------2      v  v

       <--> tweb
          <-----> bflange_out
       <--------> bflange

    """
    faces = [
        # front face
        [0, 5, 6, 7],
        [0, 1, 4, 5],
        [1, 2, 3, 4],

        # back face
        [8, 13, 14, 15],
        [8, 9, 12, 13],
        [9, 10, 11, 12],

        # around the C (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],
    ]

    zsc = 0.  # TODO: consider the shear center
    points = []
    for dim in (dim1, dim2):
        bflange_out, tweb, hin, hall = dim
        bflange = bflange_out + tweb
        tflange = (hall - hin) / 2.
        pointsi = np.array([
            [0., hall/2,  zsc], # 0
            [0., -hall/2, zsc], # 1

            [0., -hall/2,           zsc + bflange], # 2
            [0., -hall/2 + tflange, zsc + bflange], # 3
            [0., -hall/2 + tflange, zsc + tweb], # 4

            [0.,  hall/2 - tflange, zsc + tweb], # 5
            [0.,  hall/2 - tflange, zsc + bflange], # 6
            [0.,  hall/2,           zsc + bflange], # 7
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def chan1_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
                dim1: Dim4,
                dim2: Dim4):
    """builds the CHAN1 faces"""
    faces, points1, points2 = chan1_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def z_setup(dim1: Dim4,
            dim2: Dim4) -> tuple[Faces,
                                 NDArrayN3float, NDArrayN3float]:
    """
           ^ y
           |
    0--------7        ^ hall
    |      | |        |
    1----2 | |        |  ^
         | | |   z    |  | hin
         | +-|--->    |  |
         |   |        |  v
         |   6-----5  |
         |         |  |
         3---------4  v
         <---> tweb
             <-----> bfoot
       <-----------> bflange

    """
    faces = [
        # front face
        [0, 1, 2, 7],
        [2, 3, 6, 7],
        [3, 4, 5, 6],

        # back face
        [8, 9, 10, 15],
        [10, 11, 14, 15],
        [11, 12, 13, 14],

        # around the C (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],
    ]
    points = []
    for dim in (dim1, dim2):
        bfoot, tweb, hin, hall = dim
        wall = bfoot * 2. + tweb
        tflange = (hall - hin) / 2.
        pointsi = np.array([
            [0., hall/2,           -wall/2],  # 0
            [0., hall/2 - tflange, -wall/2],  # 1
            [0., hall/2 - tflange, -tweb/2],  # 2
            [0., -hall/2,          -tweb/2],  # 3
            [0., -hall/2,           wall/2], # 4
            [0., -hall/2 + tflange, wall/2], # 5
            [0., -hall/2 + tflange, tweb/2],  # 6
            [0., hall/2, tweb/2],             # 7
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def z_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
            dim1: Dim4,
            dim2: Dim4):
    """builds the Z faces"""
    faces, points1, points2 = z_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def hexa_setup(dim1: Dim3,
               dim2: Dim3) -> tuple[Faces,
                                    NDArrayN3float, NDArrayN3float]:
    r"""
          ^ y
          |
       0-----5        ^ hall
      /   |   \       |
     /    |    \      |
    1     +-----4---z |
    \           /     |
     \         /      |
      2_______3       v
              <-> wtri
    <-----------> wall

    """
    faces = [
        # front face
        [0, 1, 2, 5],
        [5, 2, 3, 4],

        # back face
        [6, 7, 8, 11],
        [11, 8, 9, 10],

        # around the C (counter-clockwise)
        [0, 6, 7, 1],
        [1, 7, 8, 2],
        [2, 8, 9, 3],
        [3, 9, 10, 4],
        [4, 10, 11, 5],
        [5, 11, 6, 0],
    ]
    points = []
    for dim in (dim1, dim2):
        wtri, wall, hall = dim
        pointsi = np.array([
            [0., hall/2,  -wall/2 + wtri],  # 0
            [0., 0.,      -wall/2],  # 1
            [0., -hall/2, -wall/2 + wtri],  # 2
            [0., -hall/2,  wall/2 - wtri],  # 3
            [0., 0.,       wall/2], # 4
            [0., hall/2,   wall/2 - wtri], # 5
        ])  # 12 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def hexa_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
               dim1: Dim3,
               dim2: Dim3) -> tuple[Faces, np.ndarray]:
    """builds the HEXA faces"""
    faces, points1, points2 = hexa_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def l_setup(dim1: Dim4,
            dim2: Dim4,
            ) -> tuple[Faces, NDArrayN3float, NDArrayN3float]:
    """
       ^y
       |
    0-----5               ^ hall
    |  |  |               |
    |  |  |               |
    |  |  |               |
    |  |  4---3           |  ^ tflange
    |  +------|-------> z |  |
    1---------2           v  v

    <----> tweb
    <---------> bflange
    """
    faces = [
        # front face
        [0, 1, 4, 5],
        [1, 2, 3, 4],

        # back face
        [6, 7, 10, 11],
        [7, 8, 9, 10],

        # around the C (counter-clockwise)
        [0, 6, 7, 1],
        [1, 7, 8, 2],
        [2, 8, 9, 3],
        [3, 9, 10, 4],
        [4, 10, 11, 5],
        [5, 11, 6, 0],
    ]
    points = []
    for dim in (dim1, dim2):
        bflange, hall, tflange, tweb = dim
        pointsi = np.array([
            [0., hall - tflange/2,  -tweb/2], # 0
            [0., -tflange/2,        -tweb/2], # 1
            [0., -tflange/2, bflange - tweb/2], # 2
            [0., tflange/2,  bflange - tweb/2], # 3
            [0., tflange/2,        tweb/2], # 4
            [0., hall - tflange/2, tweb/2], # 5
        ])  # 12 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def l_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
            dim1: Dim4,
            dim2: Dim4):
    faces, points1, points2 = l_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def t_setup(dim1: Dim4,
            dim2: Dim4) -> tuple[Faces,
                                 NDArrayN3float, NDArrayN3float]:
    """
         ^ y
         |
    0----------7
    |    +     |----> z
    1---2  5---6
        |  |
        |  |
        3--4

      bflange
    <---------->

    8----------15
    |    +     |
    9---10 13--14
        |  |
        |  |
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
    points = []
    for dim in (dim1, dim2):
        bflange, htotal, tflange, tweb = dim
        #htotal = tflange + hweb
        pointsi = np.array([
            [0., tflange/2,  -bflange/2], # 0
            [0., -tflange/2, -bflange/2], # 1
            [0., -tflange/2, -tweb/2], # 2
            [0., -htotal+tflange/2,-tweb/2], # 3
            [0., -htotal+tflange/2, tweb/2], # 4
            [0., -tflange/2, tweb/2], # 5
            [0., -tflange/2, bflange/2], # 6
            [0.,  tflange/2, bflange/2], # 7
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def t_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
            dim1: Dim4,
            dim2: Dim4):
    """builds the T faces"""
    faces, points1, points2 = t_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def t1_setup(dim1: Dim4,
             dim2: Dim4) -> tuple[Faces,
                                  NDArrayN3float, NDArrayN3float]:
    """
                 ^ y
                 |
               6---5       ^ hall
               | | |       |
               | | |       |
    0----------7 | |       | ^
    |            +-|---> z | | tflange
    1----------2   |       | v
               |   |       |
               |   |       |
               3---4       v
    <---------><--->
     bfoot      tweb

      bflange
    <---------->

    """
    faces = [
        # front face
        [0, 1, 2, 7],
        [2, 3, 4, 5, 6, 7],

        # back face
        [8, 9, 10, 15],
        [10, 11, 12, 13, 14, 15],

        # around the T (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],
    ]
    points = []
    for dim in (dim1, dim2):
        hall, bfoot, tweb, tflange = dim
        #htotal = tflange + hweb
        pointsi = np.array([
            [0.,  tflange/2,  -bfoot - tweb/2], # 0
            [0., -tflange/2, -bfoot - tweb/2], # 1
            [0., -tflange/2, -tweb/2], # 2
            [0., -hall/2,    -tweb/2], # 3
            [0., -hall/2,     tweb/2], # 4
            [0.,  hall/2,     tweb/2], # 5
            [0.,  hall/2,    -tweb/2], # 6
            [0.,  tflange/2, -tweb/2], # 7
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def t1_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
             dim1: Dim4,
             dim2: Dim4):  # validated
    """builds the T1 faces"""
    faces, points1, points2 = t1_setup(dim1, dim2)
    points_list = []
    for nid, points in [(n1, points1), (n2, points2)]:
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def t2_setup(dim1: Dim4,
             dim2: Dim4) -> tuple[Faces,
                                  NDArrayN3float, NDArrayN3float]:
    """
       <-->  tweb
         ^y
         |
       0--7    ^
       |  |    | hweb
       |  |    V
    2--1  6--5 ^ hflange
    |        |-|-->z
    3--------4 v

      bflange
    <-------->

    """
    faces = [
        # front face
        [0, 1, 6, 7],
        [1, 2, 3, 4, 5, 6],

        # back face
        [8, 9, 14, 15],
        [9, 10, 11, 12, 13, 14],

        # around the T (counter-clockwise)
        [0, 8, 9, 1],
        [1, 9, 10, 2],
        [2, 10, 11, 3],
        [3, 11, 12, 4],
        [4, 12, 13, 5],
        [5, 13, 14, 6],
        [6, 14, 15, 7],
        [7, 15, 8, 0],
    ]
    points = []
    for dim in (dim1, dim2):
        bflange, htotal, tflange, tweb = dim
        hweb = htotal - tflange
        pointsi = np.array([
            # x, y, z
            [0., hweb + tflange/2, -tweb/2],  # 0
            [0., tflange/2,        -tweb/2],  # 1
            [0., tflange/2,     -bflange/2],  # 2
            [0., -tflange/2,    -bflange/2],  # 3

            [0., -tflange/2,    bflange/2],  # 4
            [0., tflange/2,     bflange/2],  # 5
            [0., tflange/2,        tweb/2],  # 6
            [0., hweb + tflange/2, tweb/2],  # 7
        ])  # 16 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def t2_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
             dim1: Dim4,
             dim2: Dim4):  # validated
    """builds the T2 faces"""
    faces, points1, points2 = t2_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

def hat_setup(dim1: Dim4,
              dim2: Dim4) -> tuple[Faces,
                                   NDArrayN3float, NDArrayN3float]:
    """
            <--------d3------->

            0----------------11              ^   --y3
       d4   |        A        |   d4         |
     <----> +-d2-5-------6-d2-+ <---->       |   --y2
            | B  |       |  B |              | d1
     2------1----+       +----10-----9       |   --y1
     |     C     |       |     C     | t=d2  |
     3-----------4       7-----------8       v   --y0
     |
     -z0    |-z1 |-z2    |z2  |z1    |z0
    """
    faces = [
        # front face
        [1, 2, 3, 4],
        [0, 1, 4, 5],
        [0, 5, 6, 11],
        [6, 7, 10, 11],
        [7, 8, 9, 10],

        # back face
        [13, 14, 15, 16],
        [12, 13, 16, 17],
        [12, 17, 18, 23],
        [18, 19, 22, 23],
        [19, 20, 21, 22],

        # around the T (counter-clockwise)
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
    ]
    points = []
    for dim in (dim1, dim2):
        d1, d2, d3, d4 = dim
        #hall, that, wmid, bfoot
        z0 = d3/2 + d4
        z1 = d3/2.
        z2 = d3/2 - d2
        y0 = 0.
        y1 = y0 + d2
        y2 = d1 - d2
        y3 = d1
        pointsi = np.array([
            # x, y, z
            [0., y3, -z1],  # 0
            [0., y1, -z1],  # 1
            [0., y1, -z0],  # 2
            [0., y0, -z0],  # 3

            [0., y0, -z2],  # 4
            [0., y2, -z2],  # 5
            [0., y2,  z2],  # 6
            [0., y0,  z2],  # 7

            [0., y0, z0],  # 3
            [0., y1, z0],  # 2
            [0., y1, z1],  # 1
            [0., y3, z1],  # 0
        ])  # 24 x 3
        points.append(pointsi)
    return faces, points[0], points[1]

def hat_faces(n1: NDArray3float, n2: NDArray3float, xform: NDArray33float,
              dim1: Dim4,
              dim2: Dim4):
    """builds the HAT faces"""
    faces, points1, points2 = hat_setup(dim1, dim2)
    points_array = _transform_points(n1, n2, points1, points2, xform)
    return faces, points_array

