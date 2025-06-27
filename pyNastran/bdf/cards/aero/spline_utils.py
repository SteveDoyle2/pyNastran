"""

{w_aero} = [xK] [C]^-1 {w_structure}
{w_aero} = [xK] [C]^-1 {w_structure}

Eq 1-56
-------
w = [C] {P}

Eq 1-41
-------
{w_structure} = [C]{P}
so:
{P} = [C]^-1 {w_structure}

Eq 1-42
-------
           [1   x1A y1A K1A,1, K2A,2 ... K1nA,n]        {0 }
{w_aero} = [1   x2A y2A K2A,1, K2A,2 ... K2nA,n] [C]^-1 {0 }
           [... ... ... ...    ...   ... ...   ]        {0 } = [xK] [C]^-1 {w_structure}
           [1   xnA ynA KnA,1, K2A,2 ... KnnA,n]        {w1}
                                                        {w2}
                                                        {wn}
so:
{w_aero} = [xK] [C]^-1 {w_structure} = [xK] {Cws} = [xK] {P}

"""
#from typing import TYPE_CHECKING
import numpy as np
#if TYPE_CHECKING:
from pyNastran.bdf.bdf import BDF
from cpylog import SimpleLogger

natural_log = np.log


def map_deflections(structure_model: BDF, node_list: list[int],
                    deflections: dict[int, np.ndarray],
                    apoints: dict[int, np.ndarray],
                    log: SimpleLogger):
    #nodeList = remove_duplicate_nodes(node_list, structure_model)
    C = get_c_matrix(node_list, structure_model, log)
    w_structure = get_w_structure(node_list, deflections, log)
    del deflections

    #aPoints = read_half_cart3d_points(cart3d)
    w_aero = get_w_aero(node_list, C, w_structure, structure_model, apoints, log)
    del C
    del structure_model

    #write_new_cart3d_mesh(cart3d, cart3d2, wA)
    return w_aero, w_structure


def get_w_aero(node_list: list[int],
               C: np.ndarray,
               w_structure: np.ndarray,
               structure_model: BDF,
               apoints: dict[int, np.ndarray],
               log: SimpleLogger) -> dict[int, np.ndarray]:
    log.info('---starting get_w_aero---')
    #print printMatrix(C)

    Cws = np.linalg.inv(C) @ w_structure  # Cws matrix, P matrix
    #P = solve(C, wS)
    #C*P=wS
    #P = C^-1*wS

    w_aero = get_XK_matrix(Cws, node_list, structure_model, apoints, log)
    #{w_aero} = [xK] [C] {w_structure}
    log.info('---finished get_w_aero---')
    return w_aero


def get_XK_matrix(Cws: np.ndarray,
                  node_list: list[int],
                  structure_model: BDF,
                  apoints: dict[int, np.ndarray],
                  log: SimpleLogger) -> dict[int, np.ndarray]:
    """


    Parameters
    ----------
    Cws
    node_list
    structure_model
    apoints
    log

    Returns
    -------

    """
    log.info('---starting get_XK_matrix---')
    D = 1.
    piD16 = np.pi * D * 16.

    nnodes = len(node_list)
    #npoints = len(apoints.keys())
    w_aero = {}
    i = 0
    for (iaero, anode) in sorted(apoints.items()):
        xK = np.zeros(nnodes+3, dtype='float64')
        #nodeI = mesh.Node(iNode)

        xa, ya, za = anode
        xK[0] = 1.
        xK[1] = xa
        xK[2] = ya

        j = 3
        for jNode in node_list:
            sNode = structure_model.Node(jNode)
            (xs, ys, zs) = sNode.get_position()

            Rij2 = (xa-xs)**2. + (ya-ys)**2  # Rij^2
            if Rij2 == 0.:
                xK[j] = 0.
            else:
                Kij = Rij2 * natural_log(Rij2) / piD16
                xK[j] = Kij
            j += 1

        wai = xK @ Cws
        w_aero[iaero] = wai  # [0, 0]
        #print "w[%s]=%s" % (iAero, wi[0,0])
        i += 1
    log.info('---finished get_XK_matrix---')
    return w_aero


def get_w_structure(node_list: list[int],
                    deflections: dict[int, np.ndarray],
                    log: SimpleLogger) -> np.ndarray:
    log.info('---staring get_w_structure---')
    nnodes = len(node_list)
    w_structure = np.zeros((3+nnodes, 1), dtype='float64')
    i = 3
    for inode in node_list:
        dx, dy, dz = deflections[inode]
        w_structure[i, 0] = dz
        log.info("w_structure[node=%s; i=%s]=%s" % (inode, i, dz))
        i += 1
    #print(max(w_column))
    log.info('---finished get_w_structure---')

    print(f'w_structure = {w_structure}')
    w_structure_max = max(w_structure)
    print(f'w_structure_max = {w_structure_max}')
    return w_structure


def get_c_matrix(node_list: list[int], structural_model: BDF, log: SimpleLogger):
    """
    Parameters
    ----------
    node_list
    structural_model
    log

    Returns
    -------

    """
    log.info('---starting get_c_matrix---')
    D = 1.
    piD16 = np.pi * D * 16.

    nnodes = len(node_list)
    i = 3
    log.info('nnodes=%s' % nnodes)
    C = np.zeros((3+nnodes, 3+nnodes), dtype='float64')
    for iNode in node_list:
        nodeI = structural_model.Node(iNode)
        #i = iNode+3
        (xi, yi, zi) = nodeI.get_position()
        #x,y,z = p

        C[0, i] = 1.
        C[1, i] = xi
        C[2, i] = yi

        C[i, 0] = 1.
        C[i, 1] = xi
        C[i, 2] = yi

        j = 3
        for jNode in node_list:
            #j = 3+jNode
            nodeJ = structural_model.Node(jNode)
            xj, yj, zj = nodeJ.get_position()
            if i == j:
                C[i, j] = 0.
            else:
                Rij2 = (xi-xj)**2. + (yi-yj)**2  # Rij^2
                if Rij2 == 0.:
                    C[i, j] = 0.
                else:
                    Kij = Rij2 * natural_log(Rij2) / piD16
                    C[i, j] = Kij
                    #msg = "i=%s j=%s xi=%s xj=%s yi=%s yj=%s Rij2=%s Kij=%s" %(i,j,xi,xj,yi,yj,Rij2,Kij)
                    #assert isinstance(Kij,float64), msg
            j += 1
        i += 1
    log.info('---finished get_c_matrix---')
    return C


def main():
    structure_model = BDF()
    caero_id = 1001
    paero_id = 1
    p1 = [0., 0., 0.]
    p4 = [0., 10., 0.]
    x12 = 1.
    x43 = 1.
    node_list = [1, 2, 3, 4]
    structure_model.add_grid(1, [0., 0., 0.])
    structure_model.add_grid(2, [1., 0., 0.])
    structure_model.add_grid(3, [1., 10., 0.])
    structure_model.add_grid(4, [0., 10., 0.])
    igroup = 1
    structure_model.add_caero1(caero_id, paero_id, igroup, p1, x12, p4, x43)
    log = structure_model.log
    deflections = {
        1: np.array([0., 0, 1.]),
        2: np.array([0., 0., 0.]),
        3: np.array([0., 0., 10.]),
        4: np.array([0., 0., 9.]),
    }
    apoints = {
        1: np.array([0., 0., 0.]),
        2: np.array([1., 0., 0.]),
        3: np.array([1., 10., 0.]),
        4: np.array([0., 10., 0.]),
        5: np.array([0.5, 0., 0.]),
        6: np.array([0.5, 5., 0.]),
        7: np.array([0.5, 10., 0.]),
    }
    map_deflections(structure_model, node_list,
                    deflections, apoints, log)


if __name__ == '__main__':
    main()
