import os
import sys
import numpy as np
from pyNastran.bdf import read_bdf, BDF
from pyNastran.op2 import read_op2
from pyNastran.nptyping_interface import NDArray3float


def remove_duplicate_nodes(node_list, model, log=None):
    """
    Removes nodes that have the same (x,y) coordinate.
    Note that if 2 nodes with different z values are found, only 1 is returned.
    This is intentional.
    """
    node_list.sort()
    log.info("node_list_a = %s" % node_list)
    node_dict = {}
    for inode in node_list:
        x, y, unused_z = model.Node(inode).get_position()
        node_dict[(x, y)] = inode
    node_list = node_dict.values()
    node_list.sort()
    log.info("node_list_b = %s" % node_list)
    sys.stdout.flush()
    return node_list


def run_map_deflections(node_list, bdf_filename, out_filename,
                        cart3d, cart3d2, log=None):
    """
    Runs the spline deflection mapping method to morph the Cart3d model
    to a deformed mesh.

    Parameters
    ----------
    node_list : list[int]
        the list of nodes from the BDF to spline???
        is this just the SPLINE1 card???
    bdf_filename : str
        the name of the undeformed BDF model
    out_filename : str
        the name of the deformed result (OP2) model
    cart3d : str
        the name of the undeformed Cart3d model
    cart3d2 : str
        the name of the deformed Cart3d model
    log : Log(); default=None
        a python logging object

    """
    ext = os.path.splitext(out_filename)[1]
    if ext == ".op2":
        deflections = read_op2(out_filename, log=log)
    # elif ext == '.f06':
    # deflections = read_f06(out_filename)
    else:
        raise NotImplementedError("out_filename = %r" % out_filename)

    model = read_bdf(bdf_filename, xref=True, punch=False, log=log, debug=True)

    node_list = remove_duplicate_nodes(node_list, model, log=log)
    C = get_c_matrix(node_list, model, log=log)
    wS = get_ws(node_list, deflections, log=log)
    del deflections

    aero_points = [1, 2, 3]
    wA = get_wa(node_list, C, wS, model, aero_points, log=log)
    del C
    return wA, wS

def get_wa(node_list, C: NDArray3float, wS: NDArray3float,
           model: BDF,
           aero_points: dict[int, NDArray3float],
           log=None):
    """
    Cannot use solve

    [C]*[P] = [wS]
    [P] = [C]^-1*[wS]
    [wA] = [xK] [Cws]

    Parameters
    ----------
    node_list : list[int]
        list of node ids
    C : (3+nnodes, 3+nnodes) float ndarray
        ???
    wS : (???, ???) float ndarray
        [wS] = [C]*[P]
    aero_points : ???
        ???
    aero_points : : dict[int, NDArray3float]
        nid -> xyz mapping

    Returns
    -------
    wA : dict[cart3d_nid] = zi
        cart3d_nid : int
            the cart3d node id
        zi : float
            the delta z location

    """
    log.info("---starting get_wa---")

    Cws = np.linalg.inv(C) @ wS  # Cws matrix, P matrix

    wA = get_xk_matrix(Cws, node_list, model, aero_points, log=log)
    # wA = xK*C*wS
    log.info("---finished get_wa---")
    sys.stdout.flush()
    return wA


def get_xk_matrix(Cws, node_list, model: BDF,
                  aero_points: dict[int, NDArray3float], log=None):
    """
    Calculates the XK matrix to

    xK = Rij^2 = (xa-xs)^2. + (ya-ys)^2
    xK = Rij^2 * ln(Rij^2) / piD16
    """
    log.info("---starting get_xk_matrix---")
    D = 1.0
    pi_d_16 = np.pi * D * 16.0

    nnodes = len(node_list)
    # npoints = len(aero_points.keys())
    wa = {}
    for iaero, aero_node in sorted(aero_points.items()):
        xK = np.zeros(nnodes + 3, dtype='float64')
        # nodeI = mesh.Node(iNode)

        xa, ya, unused_za = aero_node

        xK[0] = 1.0
        xK[1] = xa
        xK[2] = ya

        j = 3
        for jnode in node_list:
            structural_node = model.Node(jnode)
            (xs, ys, unused_zs) = structural_node.get_position()

            Rij2 = (xa - xs) ** 2.0 + (ya - ys) ** 2  # Rij^2
            if Rij2 == 0.0:
                xK[j] = 0.0
            else:
                Kij = Rij2 * np.log(Rij2) / pi_d_16
                xK[j] = Kij
            j += 1

        wai = xK @ Cws
        wa[iaero] = wai[0, 0]
        #print("w[%s]=%s" % (iaero, wi[0, 0]))
    #print('---wa---')
    #print('wa = ', wa)
    log.info("---finished getXK_matrix---")
    sys.stdout.flush()
    return wa


def get_ws(node_list: list[int], deflections, log=None):
    """
    Parameters
    ----------
    node_list : list[int]
        list of node ids

    """
    log.info("---staring get_ws---")
    nnodes = len(node_list)
    w_column = np.zeros((3 + nnodes, 1), dtype="float64")
    i = 3
    nodes = deflections.node_grid[:, 0]
    for inode in node_list:
        inodei = np.searchsorted(nodes, inode)
        unused_dx, unused_dy, dz = deflections.data[0, inodei, :2]  # deflections[inode]
        w_column[i] = dz
        log.info(f"wS[{inode}={i}]={dz}")
        i += 1
    print('w_column =', w_column.max())
    log.info("---finished get_ws---")
    sys.stdout.flush()

    ws_max = np.max(w_column)
    log.debug(f"ws_max = {ws_max}")
    return w_column


def get_c_matrix(node_list: list[int], model: BDF, log=None):
    """
    Parameters
    ----------
    node_list : list[int]
        list of node ids

    Returns
    -------
    C : (3+nnodes, 3+nnodes) float ndarray
        ???
    """
    log.info("---starting get_c_matrix---")
    D = 1.0
    pi_d_16 = np.pi * D * 16.0

    nnodes = len(node_list)
    i = 3
    log.info("nnodes=%s" % nnodes)
    sys.stdout.flush()

    C = np.zeros((3 + nnodes, 3 + nnodes), dtype='float64')
    for inode in node_list:
        node_i = model.Node(inode)
        # i = inode+3
        (xi, yi, unused_zi) = node_i.get_position()
        # x,y,z = p

        C[0, i] = 1.0
        C[1, i] = xi
        C[2, i] = yi

        C[i, 0] = 1.0
        C[i, 1] = xi
        C[i, 2] = yi

        j = 3
        for jnode in node_list:
            # j = 3+jnode
            node_j = model.Node(jnode)
            xj, yj, unused_zj = node_j.get_position()
            if i == j:
                C[i, j] = 0.0
            else:
                Rij2 = (xi - xj) ** 2.0 + (yi - yj) ** 2  # Rij^2
                if Rij2 == 0.0:
                    C[i, j] = 0.0
                else:
                    Kij = Rij2 * np.log(Rij2) / pi_d_16
                    C[i, j] = Kij
                    # msg = "i=%s j=%s xi=%s xj=%s yi=%s yj=%s Rij2=%s Kij=%s" %(
                    # i, j, xi, xj, yi, yj, Rij2, Kij)
                    # assert isinstance(Kij,float64), msg
            j += 1
        i += 1
    log.info("---finished getCmatrix---")
    sys.stdout.flush()
    return C
