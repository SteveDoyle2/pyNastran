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
from typing import Tuple
import numpy as np
from pyNastran.bdf.cards.aero.utils import elements_from_quad, tri_cap


def rod_faces(n1: int, n2: int, xform,
              dim1: Tuple[float, float],
              dim2: Tuple[float, float]): # validated
    """defines points in a circle with triangle based end caps"""
    # 4,8,12,16,... becomes 5,9,13,17,...
    ntheta = 17
    thetas = np.radians(np.linspace(0., 360., ntheta))

    nfaces = 0
    all_faces = []
    points_list = []
    x = np.zeros(ntheta)
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        radius, = dim
        y = radius * np.cos(thetas)
        z = radius * np.sin(thetas)
        xyz = np.vstack([x, y, z]).T
        assert xyz.shape == (ntheta, 3), xyz.shape

        pointsi = xyz @ xform + nid
        points_list.append(pointsi)

        # the tri_cap buils triangles that fan out from the first node
        tris = tri_cap(ntheta)

        # we need to use the tolist because we're going to
        # combine quads and tris (the elements have different
        # lengths)
        all_faces += (nfaces + tris).tolist()
        nfaces += tris.shape[0]

    # the main cylinder uses the points defined independent
    # of the points n1/n2
    faces = elements_from_quad(2, ntheta)
    all_faces += faces.tolist()

    # used by the tri_caps
    points_list.append(n1)
    points_list.append(n2)
    points = np.vstack(points_list)
    return all_faces, points, points.shape[0]

def tube_faces(n1: int, n2: int, xform,
               dim1: Tuple[float, float],
               dim2: Tuple[float, float]):  # validated
    """defines a rod with a hole"""
    # 4,8,12,16,... becomes 5,9,13,17,...
    ntheta = 17
    thetas = np.radians(np.linspace(0., 360., ntheta))
    npoints = ntheta * 2

    points_list1 = []
    points_list2 = []
    x = np.zeros(ntheta)

    for nid, dim in [(n1, dim1), (n2, dim2)]:
        radius_out, radius_in = dim

        # outer rod
        y = radius_out * np.cos(thetas)
        z = radius_out * np.sin(thetas)
        xyz1 = np.vstack([x, y, z]).T
        points1i = xyz1 @ xform + nid
        points_list2.append(points1i)

        # inner rod
        y = radius_in * np.cos(thetas)
        z = radius_in * np.sin(thetas)
        xyz2 = np.vstack([x, y, z]).T
        points2i = xyz2 @ xform + nid
        points_list1.append(points2i)

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
    points = np.vstack(points_list1 + points_list2)
    return all_faces, points, points.shape[0]

def bar_faces(n1: int, n2: int, xform,
              dim1: Tuple[float, float],
              dim2: Tuple[float, float]):  # validated
    """
       ^y
       |
    0-----3
    |     |
    |     |----> z
    |     |
    1-----2

    """
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        w, h = dim
        points = np.array([
            [0.,  h/2, -w/2],   # 0
            [0., -h/2, -w/2],   # 1
            [0., -h/2,  w/2],   # 2
            [0.,  h/2,  w/2],   # 3
        ])  # 16 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return np.vstack(points_list)

def box_faces(n1: int, n2: int, xform,
              dim1: Tuple[float, float, float, float],
              dim2: Tuple[float, float, float, float]):  # validated
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        wbox, hbox, th, tw = dim
        hbox_in = hbox - 2*th
        wbox_in = wbox - 2*tw
        points = np.array([
            [0., -hbox/2., wbox/2], # 0
            [0., -hbox/2, -wbox/2], # 1
            [0.,  hbox/2, -wbox/2], # 2
            [0.,  hbox/2,  wbox/2], # 3

            [0., -hbox_in/2., wbox_in/2], # 4
            [0., -hbox_in/2, -wbox_in/2], # 5
            [0.,  hbox_in/2, -wbox_in/2], # 6
            [0.,  hbox_in/2,  wbox_in/2], # 7
        ])
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def i_faces(n1: int, n2: int, xform,
            dim1: Tuple[float, float, float, float, float, float],
            dim2: Tuple[float, float, float, float, float, float]):   # validated
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        hall, bflange_btm, bflange_top, tweb, tflange_btm, tflange_top = dim
        hweb = hall - tflange_top - tflange_btm
        ysc = -hall / 2.  # TODO: fix the shear center formula

        y0 = ysc
        y1 = y0 + tflange_btm
        y2 = y1 + hweb
        y3 = y2 + tflange_top
        points = np.array([
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
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def i1_faces(n1: int, n2: int, xform,
             dim1: Tuple[float, float, float, float],
             dim2: Tuple[float, float, float, float]):
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

    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
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
        points = np.array([
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
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def h_faces(n1: int, n2: int, xform,
            dim1: Tuple[float, float, float, float],
            dim2: Tuple[float, float, float, float]):
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        winner, wouter, hall, hinner = dim
        w_all = winner + wouter

        points = np.array([
            [0.,  hall/2, -w_all/2],        # 0 - upper far left
            [0., -hall/2, -w_all/2],        # 1 - lower far left
            [0., -hall/2,   -winner/2],     # 2 - lower near left
            [0., -hinner/2, -winner/2],     # 3 - lower H left
            [0., -hinner/2,  winner/2],     # 4 - lower H right
            [0., -hall/2,    winner/2],     # 5 - lower near right
            [0., -hall/2,  w_all/2],         # 6 - lower far right

            [0., hall/2,  w_all/2],         # 7  - upper far right
            [0., hall/2,    winner/2],     # 8  - upper near right
            [0., hinner/2,  winner/2],     # 9  - upper H right
            [0., hinner/2, -winner/2],     # 10 - upper H left
            [0., hall/2,   -winner/2],     # 11 - upper near left
        ])  # 24 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def chan_faces(n1: int, n2: int, xform,
               dim1: Tuple[float, float, float, float],
               dim2: Tuple[float, float, float, float]):
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
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

        points = np.array([
            [0., hall/2,  zsc], # 0
            [0., -hall/2, zsc], # 1

            [0., -hall/2,           zsc + bflange], # 2
            [0., -hall/2 + tflange, zsc + bflange], # 3
            [0., -hall/2 + tflange, zsc + tweb], # 4

            [0.,  hall/2 - tflange, zsc + tweb], # 5
            [0.,  hall/2 - tflange, zsc + bflange], # 6
            [0.,  hall/2,           zsc + bflange], # 7
        ])  # 16 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def chan1_faces(n1: int, n2: int, xform,
                dim1: Tuple[float, float, float, float],
                dim2: Tuple[float, float, float, float]):
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
    zsc = 0.  # TODO: consider the shear center
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        bflange_out, tweb, hin, hall = dim
        bflange = bflange_out + tweb
        tflange = (hall - hin) / 2.
        points = np.array([
            [0., hall/2,  zsc], # 0
            [0., -hall/2, zsc], # 1

            [0., -hall/2,           zsc + bflange], # 2
            [0., -hall/2 + tflange, zsc + bflange], # 3
            [0., -hall/2 + tflange, zsc + tweb], # 4

            [0.,  hall/2 - tflange, zsc + tweb], # 5
            [0.,  hall/2 - tflange, zsc + bflange], # 6
            [0.,  hall/2,           zsc + bflange], # 7
        ])  # 16 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def z_faces(n1: int, n2: int, xform,
            dim1: Tuple[float, float, float, float],
            dim2: Tuple[float, float, float, float]):
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        bfoot, tweb, hin, hall = dim
        wall = bfoot * 2. + tweb
        tflange = (hall - hin) / 2.
        points = np.array([
            [0., hall/2,           -wall/2],  # 0
            [0., hall/2 - tflange, -wall/2],  # 1
            [0., hall/2 - tflange, -tweb/2],  # 2
            [0., -hall/2,          -tweb/2],  # 3
            [0., -hall/2,           wall/2], # 4
            [0., -hall/2 + tflange, wall/2], # 5
            [0., -hall/2 + tflange, tweb/2],  # 6
            [0., hall/2, tweb/2],             # 7
        ])  # 16 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def hexa_faces(n1: int, n2: int, xform,
               dim1: Tuple[float, float, float],
               dim2: Tuple[float, float, float]):
    """
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        wtri, wall, hall = dim
        points = np.array([
            [0., hall/2,  -wall/2 + wtri],  # 0
            [0., 0.,      -wall/2],  # 1
            [0., -hall/2, -wall/2 + wtri],  # 2
            [0., -hall/2,  wall/2 - wtri],  # 3
            [0., 0.,       wall/2], # 4
            [0., hall/2,   wall/2 - wtri], # 5
        ])  # 12 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def l_faces(n1: int, n2: int, xform,
            dim1: Tuple[float, float, float, float],
            dim2: Tuple[float, float, float, float]):
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        bflange, hall, tflange, tweb = dim
        points = np.array([
            [0., hall - tflange/2,  -tweb/2], # 0
            [0., -tflange/2,        -tweb/2], # 1
            [0., -tflange/2, bflange - tweb/2], # 2
            [0., tflange/2,  bflange - tweb/2], # 3
            [0., tflange/2,        tweb/2], # 4
            [0., hall - tflange/2, tweb/2], # 5
        ])  # 12 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)

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
    return faces, np.vstack(points_list)

def t_faces(n1: int, n2: int, xform,
            dim1: Tuple[float, float, float, float],
            dim2: Tuple[float, float, float, float]):  # validated
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

    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        bflange, htotal, tflange, tweb = dim
        #htotal = tflange + hweb
        points = np.array([
            [0., tflange/2,  -bflange/2], # 0
            [0., -tflange/2, -bflange/2], # 1
            [0., -tflange/2, -tweb/2], # 2
            [0., -htotal+tflange/2,-tweb/2], # 3
            [0., -htotal+tflange/2, tweb/2], # 4
            [0., -tflange/2, tweb/2], # 5
            [0., -tflange/2, bflange/2], # 6
            [0.,  tflange/2, bflange/2], # 7
        ])  # 16 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def t1_faces(n1: int, n2: int, xform,
             dim1: Tuple[float, float, float, float],
             dim2: Tuple[float, float, float, float]):  # validated
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        hall, bfoot, tweb, tflange = dim
        #htotal = tflange + hweb
        points = np.array([
            [0.,  tflange/2,  -bfoot - tweb/2], # 0
            [0., -tflange/2, -bfoot - tweb/2], # 1
            [0., -tflange/2, -tweb/2], # 2
            [0., -hall/2,    -tweb/2], # 3
            [0., -hall/2,     tweb/2], # 4
            [0.,  hall/2,     tweb/2], # 5
            [0.,  hall/2,    -tweb/2], # 6
            [0.,  tflange/2, -tweb/2], # 7
        ])  # 16 x 3
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def t2_faces(n1: int, n2: int, xform,
             dim1: Tuple[float, float, float, float],
             dim2: Tuple[float, float, float, float]):  # validated
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        bflange, htotal, tflange, tweb = dim
        hweb = htotal - tflange
        points = np.array([
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
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)


def hat_faces(n1: int, n2: int, xform,
              dim1: Tuple[float, float],
              dim2: Tuple[float, float]):
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
    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        d1, d2, d3, d4 = dim
        #hall, that, wmid, bfoot
        z0 = d3/2 + d4
        z1 = d3/2.
        z2 = d3/2 - d2
        y0 = 0.
        y1 = y0 + d2
        y2 = d1 - d2
        y3 = d1
        points = np.array([
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
        pointsi = points @ xform + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)
