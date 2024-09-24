import numpy as np
import matplotlib.pyplot as plt
from pyNastran.bdf.cards.aero.aero import CAERO1

def build_caero_fore_aft(
        x12: float, x43: float,
        a_point_location12: int=1,
        b_point_location43: int=3):
    pa = [0., 0., 0.]
    pb = [0., 10., 0.]
    assert a_point_location12 in [1, 2], a_point_location12
    assert b_point_location43 in [3, 4], b_point_location43

    dxyz12 = np.array([x12, 0., 0.])
    dxyz43 = np.array([x43, 0., 0.])
    point_pa = np.asarray(pa)
    point_pb = np.asarray(pb)
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

def build_caero_dimensions(x12: float, x43: float, b: float,
                           sweep_deg: float=0.0,
                           dihedral_deg: float=0.0,
                           percent_chord1_sweep: float=0.25,
                           percent_chord2_sweep: float=0.25):
    pa = [0., 0., 0.]
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

    pa_1a = points2[0, :]
    pa_2a = points2[1, :]
    pa_3a = points2[2, :]
    pa_4a = points2[3, :]
    dn = 0.2
    centroid2 = (pa_1a + pa_2a + pa_3a + pa_4a) / 4
    pa_1b = pa_1a + (centroid2 - pa_1a)*dn
    pa_2b = pa_2a + (centroid2 - pa_2a)*dn
    pa_3b = pa_3a + (centroid2 - pa_3a)*dn
    pa_4b = pa_4a + (centroid2 - pa_4a)*dn

    # create points that are at the 20% spanwise location
    span_split = 0.20
    p4_prime = p1 + (p4 - p1) * span_split
    p3_prime = p2 + (p3 - p2) * span_split
    points_split = np.vstack([p3_prime, p4_prime])
    points2_split = points_split @ T
    pb_1a = points2_split[0, :]
    pb_2a = points2_split[1, :]

    if dz/b > 0.1 or 1:
        ax = fig.add_subplot(projection='3d')
        ax.plot(x, y, z, label='caero panel')
        ax.set_zlabel('z')

        # ax.text(x, y, z, label, zdir)
        #ax.text(centroid2[0], centroid2[1], centroid2[2], 'centroid')
        ax.text(*centroid2, 'centroid')
        ax.text(*pa_1b, 'p1')
        ax.text(*pa_2b, 'p2')
        ax.text(*pa_3b, 'p3')
        ax.text(*pa_4b, 'p4')


        ax.text(*pb_1a, 'pa')
        ax.text(*pb_2a, 'pb')
    else:
        ax = fig.gca()
        ax.grid(True)
        ax.plot(x, y, label='caero panel')
        #ax.set_text(p1a[:1],"p1")
        # n = (centroid - p1)
        # p1' = p1 + n*scale
        ax.annotate('centroid', centroid2[:2], ha='center')
        ax.annotate('p1', pa_1b[:2], ha='center')
        ax.annotate('p2', pa_2b[:2], ha='center')
        ax.annotate('p3', pa_3b[:2], ha='center')
        ax.annotate('p4', pa_4b[:2], ha='center')

        ax.annotate('pa', pb_1a[:2], ha='center')
        ax.annotate('pb', pb_2a[:2], ha='center')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')

    plt.show()


if __name__ == '__main__':
    x12 = 5.0
    x43 = 1.0
    b = 11.
    sweep_deg = 10.
    dihedral_deg = 20.
    build_caero_dimensions(x12, x43, b,
                           sweep_deg=sweep_deg,
                           dihedral_deg=dihedral_deg)

    x12 = 1.0
    x43 = 1.0
    #build_caero_fore_aft(x12, x43)
