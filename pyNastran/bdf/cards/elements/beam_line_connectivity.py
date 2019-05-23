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
import numpy as np
from pyNastran.bdf.cards.aero.utils import elements_from_quad, tri_cap


def rod_faces(n1, n2, xform, dim1, dim2): # validated
    """
    defines points in a circle with triangle based end caps
    """
    # 4,8,12,16,... becomes 5,9,13,17,...
    thetas = np.radians(np.linspace(0., 360., 17))
    ntheta = len(thetas)

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

        # the tri_cap is made from points that aren't defined yet
        # (the n1/n2 end points)
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

def tube_faces(n1, n2, xform, dim1, dim2):  # validated
    """
    defines a rod with a hole
    """
    # 4,8,12,16,... becomes 5,9,13,17,...
    thetas = np.radians(np.linspace(0., 360., 17))
    ntheta = len(thetas)
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

def bar_lines(dim):  # validaeted
    """
       ^y
       |
    0----3
    |    |
    |    |----> z
    |    |
    1----2
    """
    lines = [[0, 1, 2, 3, 0]]
    w, h = dim
    points = np.array([
        [h/2, -w/2],   # 0
        [-h/2, -w/2],   # 1
        [-h/2,  w/2],   # 2
        [h/2,  w/2],   # 3
        ])  # 16 x 3
    return lines, points

def box_lines(dim):  # validated
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
    lines = [
        [0, 1, 2, 3, 0],
        [4, 5, 6, 7, 4],
    ]

    points_list = []
    wbox, hbox, th, tw = dim
    hbox_in = hbox - 2*th
    wbox_in = wbox - 2*tw
    points = np.array([
        [-hbox/2., wbox/2], # 0
        [-hbox/2, -wbox/2], # 1
        [hbox/2, -wbox/2], # 2
        [hbox/2,  wbox/2], # 3

        [-hbox_in/2., wbox_in/2], # 4
        [-hbox_in/2, -wbox_in/2], # 5
        [hbox_in/2, -wbox_in/2], # 6
        [hbox_in/2,  wbox_in/2], # 7
    ])
    return lines, points

def i_lines(dim):   # validated
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0]]

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
    return lines, points

def i1_lines(dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0]]
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
    return lines, points

def h_lines(dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0]]
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
    return lines, points

def chan_faces(n1, n2, xform, dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]

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
    return lines, points

def chan1_lines(n1, n2, xform, dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]
    zsc = 0.  # TODO: consider the shear center
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
    return lines, points

def z_lines(dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]
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
    return lines, points

def hexa_lines(dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 0]]
    wtri, wall, hall = dim
    points = np.array([
        [0., hall/2,  -wall/2 + wtri],  # 0
        [0., 0.,      -wall/2],  # 1
        [0., -hall/2, -wall/2 + wtri],  # 2
        [0., -hall/2,  wall/2 - wtri],  # 3
        [0., 0.,       wall/2], # 4
        [0., hall/2,   wall/2 - wtri], # 5
        ])  # 12 x 3
    return lines, points

def l_lines(dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 0]]
    bflange, hall, tflange, tweb = dim
    points = np.array([
        [0., hall - tflange/2,  -tweb/2], # 0
        [0., -tflange/2,        -tweb/2], # 1
        [0., -tflange/2, bflange - tweb/2], # 2
        [0., tflange/2,  bflange - tweb/2], # 3
        [0., tflange/2,        tweb/2], # 4
        [0., hall - tflange/2, tweb/2], # 5
        ])  # 12 x 3
    return lines, points

def t_lines(dim):  # validated
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
    """
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]
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
    return lines, points

def t1_lines(dim):  # validated
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]
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
    return lines, points

def t2_lines(n1, n2, xform, dim):  # validated
    """
         ^y
         |
       0--7
       |  |
       |  |
    2--1  6--5
    |        |---->z
    3--------4

      bflange
    <---------->
    """
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]
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
    return lines, points


def hat_lines(dim):
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
    lines = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0]]
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
    return lines, points
