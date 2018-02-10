import numpy as np

def bar_faces(n1, n2, ihat, yhat, zhat, dim1, dim2):  # validaeted
    """
       ^y
       |
    0----3
    |    |
    |    |----> z
    |    |
    1----2
    """
    xform = np.vstack([ihat, yhat, zhat]) # 3x3 unit matrix

    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        w, h = dim
        points = np.array([
            [0.,  h/2, -w/2],   # 0
            [0., -h/2, -w/2],   # 1
            [0., -h/2,  w/2],   # 2
            [0.,  h/2,  w/2],   # 3
        ])  # 16 x 3
        pointsi = np.dot(points, xform) + nid
        points_list.append(pointsi)
    return np.vstack(points_list)

def box_faces(n1, n2, ihat, yhat, zhat, dim1, dim2):  # validaeted
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
    xform = np.vstack([ihat, yhat, zhat]) # 3x3 unit matrix
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
            [0., -wbox/2., hbox/2], # 0
            [0., -wbox/2, -hbox/2], # 1
            [0.,  wbox/2, -hbox/2], # 2
            [0.,  wbox/2,  hbox/2], # 3

            [0., -wbox_in/2., hbox_in/2], # 4
            [0., -wbox_in/2, -hbox_in/2], # 5
            [0.,  wbox_in/2, -hbox_in/2], # 6
            [0.,  wbox_in/2,  hbox_in/2], # 7
        ])
        pointsi = np.dot(points, xform) + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def i_faces(n1, n2, ihat, yhat, zhat, dim1, dim2):
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

    xform = np.vstack([ihat, yhat, zhat]) # 3x3 unit matrix

    ## TODO: most of the faces are missing; I don't know why
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
        ysc =  -hall / 2.  # TODO: fix the shear center formula

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
        ])  # 16 x 3
        pointsi = np.dot(points, xform) + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def t_faces(n1, n2, ihat, yhat, zhat, dim1, dim2):  # validated
    """
         ^ y
         |
    0----------7
    |    +     |----> z
    1---2  5---6
        |  |
        |  |
        3--4

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
    xform = np.vstack([ihat, yhat, zhat]) # 3x3 unit matrix

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
        pointsi = np.dot(points, xform) + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)

def t2_faces(n1, n2, ihat, yhat, zhat, dim1, dim2):  # validated
    """
         ^y
         |
       0--7
       |  |
       |  |
    2--1  6--5
    |        |---->z
    3--------4
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

    xform = np.vstack([ihat, yhat, zhat]) # 3x3 unit matrix

    points_list = []
    for nid, dim in [(n1, dim1), (n2, dim2)]:
        bflange, htotal, tflange, tweb = dim
        hweb = htotal - tflange
        points = np.array([
            # x, y, z
            [0., hweb + tflange/2, -tflange/2],  # 0
            [0., tflange/2, -tflange/2],  # 1
            [0., tflange/2, -bflange/2],   # 2
            [0., -tflange/2,-bflange/2], # 3

            [0., -tflange/2,  bflange/2], # 4
            [0., tflange/2, bflange/2],   # 5
            [0., tflange/2, tflange/2],  # 6
            [0., hweb + tflange/2, tflange/2],  # 7
        ])  # 16 x 3
        pointsi = np.dot(points, xform) + nid
        points_list.append(pointsi)
    return faces, np.vstack(points_list)
