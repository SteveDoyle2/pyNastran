"""
defines:
    * nids_close = find_closest_nodes(nodes_xyz, nids, xyz_compare, neq_max, tol)
    * ieq = find_closest_nodes_index(nodes_xyz, xyz_compare, neq_max, tol)

"""
from itertools import count
from typing import List, Optional
import numpy as np

from pyNastran.bdf.mesh_utils.bdf_equivalence import (
    _get_tree)

from pyNastran.nptyping import NDArray3float, NDArrayNint

def find_closest_nodes(nodes_xyz: NDArray3float, nids: NDArrayNint,
                       xyz_compare: NDArray3float, neq_max: int=1, tol: Optional[float]=None,
                       msg: str='') -> NDArrayNint:
    """
    Finds the closest nodes to an arbitrary set of xyz points

    Parameters
    ----------
    nodes_xyz : (Nnodes, 3) float ndarray
        the source points (e.g., xyz_cid0)
    nids : (Nnodes, ) int ndarray
        the source node ids (e.g.; nid_cp_cid[:, 0])
    xyz_compare : (Ncompare, 3) float ndarray
        the xyz points to compare to; xyz_to_find
    tol : float; default=None
        the max spherical tolerance
        None : the whole model
    neq_max : int; default=1.0
        the number of "close" points
    msg : str; default=''
        custom message used for errors

    Returns
    -------
    nids_close: (Ncompare, ) int ndarray
        the close node ids

    """
    if not isinstance(neq_max, int):
        msgi = 'neq_max=%r must be an int; type=%s\n%s' % (
            neq_max, type(neq_max), msg)
        raise TypeError(msgi)
    #ieq = find_closest_nodes_index(nodes_xyz, xyz_compare, neq_max, tol)
    if tol is None:
        xyz_max = nodes_xyz.max(axis=0)
        xyz_min = nodes_xyz.min(axis=0)
        assert len(xyz_max) == 3, xyz_max
        dxyz = np.linalg.norm(xyz_max - xyz_min)
        tol = 2. * dxyz

    ieq = _not_equal_nodes_build_tree(nodes_xyz, xyz_compare, tol,
                                      neq_max=neq_max, msg=msg)[1]
    ncompare = xyz_compare.shape[0]
    assert len(ieq) == ncompare, 'increase the tolerance so you can find nodes; tol=%r' % tol
    try:
        nids_out = nids[ieq]
    except IndexError:
        # if you get a crash while trying to create the error message
        # check to see if your nodes are really far from each other
        #
        nnids = len(nids)
        msgi = 'Cannot find:\n'
        for i, ieqi, nid in zip(count(), ieq, nids):
            if ieqi == nnids:
                xyz = xyz_compare[i, :]
                msgi += '  nid=%s xyz=%s\n' % (nid, xyz)
        msgi += msg
        raise IndexError(msgi)
    return nids_out


def find_closest_nodes_index(nodes_xyz, xyz_compare, neq_max, tol, msg=''):
    """
    Finds the closest nodes to an arbitrary set of xyz points

    Parameters
    ----------
    nodes_xyz : (Nnodes, 3) float ndarray
        the source points
    xyz_compare : (Ncompare, 3) float ndarray
        the xyz points to compare to
    neq_max : int
        the number of "close" points (default=4)
    tol : float
        the max spherical tolerance
    msg : str; default=''
        error message

    Returns
    -------
    slots : (Ncompare, ) int ndarray
        the indices of the close nodes corresponding to nodes_xyz

    """
    #nodes_xyz, model, nids, inew = _eq_nodes_setup(
        #bdf_filename, tol, renumber_nodes=renumber_nodes,
        #xref=xref, node_set=node_set, debug=debug)
    ieq, slots = _not_equal_nodes_build_tree(nodes_xyz, xyz_compare, tol,
                                             neq_max=neq_max, msg=msg)[1:3]
    return ieq


def _not_equal_nodes_build_tree(nodes_xyz, xyz_compare, tol, neq_max=4, msg=''):
    # type: (np.ndarray, np.ndarray, float, int, str) -> (Any, np.ndarray, np.ndarray)
    """
    helper function for `bdf_equivalence_nodes`

    Parameters
    ----------
    nodes_xyz : (Nnodes, 3) float ndarray
         the source points
    xyz_compare : (Ncompare, 3) float ndarray
         the xyz points to compare to
    tol : float
        the max spherical tolerance
    neq_max : int; default=4
        the number of close nodes
    msg : str; default=''
        error message

    Returns
    -------
    kdt : cKDTree()
        the kdtree object
    ieq : int ndarray
        The indices of nodes_xyz where the nodes in xyz_compare are close???
        neq_max = 1:
            (N, ) int ndarray
        neq_max > 1:
            (N, N) int ndarray
    slots : int ndarray
        The indices of nodes_xyz where the nodes in xyz_compare are close???
        neq_max = 1:
            (N, ) int ndarray
        neq_max > 1:
            (N, N) int ndarray
    msg : str; default=''
        error message

    """
    assert isinstance(nodes_xyz, np.ndarray), type(nodes_xyz)
    assert isinstance(xyz_compare, np.ndarray), type(xyz_compare)
    if len(nodes_xyz.shape) != len(xyz_compare.shape) or nodes_xyz.shape[1] != xyz_compare.shape[1]:
        msgi = 'nodes_xyz.shape=%s xyz_compare.shape=%s%s' % (
            str(nodes_xyz.shape), str(xyz_compare.shape), msg)
        raise RuntimeError(msgi)
    kdt = _get_tree(nodes_xyz, msg=msg)
    # check the closest 10 nodes for equality
    deq, ieq = kdt.query(xyz_compare, k=neq_max, distance_upper_bound=tol)
    #print(deq)
    #print('ieq =', ieq)
    #print('neq_max = %s' % neq_max)

    # get the ids of the duplicate nodes
    nnodes = nodes_xyz.shape[0]
    if neq_max == 1:
        assert len(deq.shape) == 1, deq.shape
        slots = np.where(ieq < nnodes)
    else:
        assert len(deq.shape) == 2, deq.shape
        slots = np.where(ieq[:, :] < nnodes)
    #print('slots =', slots)
    return kdt, ieq, slots
