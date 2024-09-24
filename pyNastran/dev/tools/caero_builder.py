import numpy as np
import matplotlib.pyplot as plt
from pyNastran.bdf.cards.aero.aero import CAERO1

def build_caero_fore_aft():
    a_point_location12 = 1
    x12 = 1.0
    x43 = 1.0
    pa = [0., 0., 0.]
    pb = [0., 10., 0.]
    b_point_location43 = 3
    assert a_point_location12 in [1, 2], a_point_location12
    assert b_point_location43 in [3, 4], b_point_location43

    dxyz12 = np.array([x12, 0., 0.])
    dxyz43 = np.array([x43, 0., 0.])
    point_pa = np.array(pa)
    point_pb = np.array(pb)
    assert len(point_pa) == 3, point_pa
    assert len(point_pb) == 3, point_pb
    if a_point_location12 == 1:
        p1 = point_pa
        p2 = p1 + dxyz12
    else:
        assert a_point_location12 == 2, a_point_location12
        p2 = point_pa
        p1 = p2 - dxyz12

    if b_point_location43 == 3:
        p3 = point_pb
        p4 = p3 + dxyz43
    else:
        assert b_point_location43 == 3, b_point_location43
        p4 = point_pb
        p3 = p4 - dxyz43

    eid = 1
    pid = 1
    #eid: int, pid: int, igroup: int,
                 #p1: NDArray3float, x12: float,
                 #p4: NDArray3float, x43: float,
                 #cp: int=0,
                 #nspan: int=0, lspan: int=0,
                 #nchord: int=0, lchord: int=0, 
    igroup = 1
    elem = CAERO1(
        eid, pid, igroup,
        p1, x12,
        p4, x43,
        cp=0, nspan=10, nchord=10)
    elem.validate()
    return elem

def build_caero_dimensions():
    pa = [0., 0., 0.]
    x12 = 5.0
    x43 = 1.0
    b = 11.

    percent_chord1_sweep = 0.25
    percent_chord2_sweep = 0.25

    sweep_deg = 10.
    dihedral_deg = 20.

    dihedral = np.radians(dihedral_deg)
    sweep = np.radians(sweep_deg)

    #-----------------------------
    chord1 = x12
    chord2 = x43

    point_pa = np.array(pa)
    p1 = point_pa
    dxyz12 = np.array([x12, 0., 0.])
    dxyz43 = np.array([x43, 0., 0.])
    d12_34 = np.array([b*np.tan(sweep), b, 0.])
    d1_p12 = np.array([chord1*percent_chord1_sweep, 0., 0.])
    d4_p43 = np.array([chord2*percent_chord2_sweep, 0., 0.])
    p12 = p1 + d1_p12
    p43 = p12 + d12_34
    p4 = p43 - d4_p43
    p3 = p4 + dxyz43
    p2 = p1 + dxyz12

    dsweep = p4 - p1
    sweep_LE = np.arctan2(dsweep[0], dsweep[1])
    sweep_LE_deg = np.degrees(sweep_LE)
    print(f'sweepLE (deg) = {sweep_LE_deg:g}')
    print('dsweep = ', dsweep)

    # stack points in ring
    points = np.vstack([p1, p2, p3, p4, p1])
    

    # transform points from xy plane to xyz
    cosd = np.cos(dihedral)
    sind = np.sin(dihedral)
    T = np.array([
        [1., 0., 0.],
        [0., cosd, sind],
        [0., -sind, cosd],
    ])
    points2 = points @ T
    
    fig = plt.figure()

    x = points2[:, 0]
    y = points2[:, 1]
    z = points2[:, 2]
    #print('z/b = ', z)
    dz = z.max() - z.min()
    #print('dz/b = ', dz/b)

    p1a = points2[0, :]
    p2a = points2[1, :]
    p3a = points2[2, :]
    p4a = points2[3, :]
    dn = 0.2
    centroid2 = (p1a + p2a + p3a + p4a) / 4
    p1b = p1a + (centroid2 - p1a)*dn
    p2b = p2a + (centroid2 - p2a)*dn
    p3b = p3a + (centroid2 - p3a)*dn
    p4b = p4a + (centroid2 - p4a)*dn

    span_split = 0.20
    bi = b*span_split
    p4_prime = p1 + np.array([bi*span_split, bi, 0.])
    p3_prime = p2 + np.array([bi*span_split, bi, 0.])
    points_split = np.vstack([p3_prime, p4_prime])
    points2_split = points_split @ T

    if dz/b > 0.1 or 1:
        ax = fig.add_subplot(projection='3d')
        ax.plot(x, y, z, label='caero panel')
        ax.set_zlabel('z')

        # ax.text(x, y, z, label, zdir)
        #ax.text(centroid2[0], centroid2[1], centroid2[2], 'centroid')
        ax.text(*centroid2, 'centroid')
        ax.text(*p1b, 'p1')
        ax.text(*p2b, 'p2')
        ax.text(*p3b, 'p3')
        ax.text(*p4b, 'p4')
    else:
        ax = fig.gca()
        ax.grid(True)
        ax.plot(x, y, label='caero panel')
        #ax.set_text(p1a[:1],"p1")
        # n = (centroid - p1)
        # p1' = p1 + n*scale
        ax.annotate('centroid', centroid2[:2], ha='center')
        ax.annotate('p1', p1b[:2], ha='center')
        ax.annotate('p2', p2b[:2], ha='center')
        ax.annotate('p3', p3b[:2], ha='center')
        ax.annotate('p4', p4b[:2], ha='center')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')

    plt.show()


if __name__ == '__main__':
    build_caero_dimensions()
    #build_caero_fore_aft()

