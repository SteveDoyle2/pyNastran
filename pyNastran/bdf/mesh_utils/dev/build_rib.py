from pyNastran.bdf.bdf import BDF
import numpy as np

def build_rib(model: BDF,
              nids_top: list[int],
              nids_mid: list[int],
              nids_bottom: list[int],
              pid: int=1,
              nid0: int=0,
              eid0: int=0):

    """
    T1---T2----T3----3
                |
    *     *    M1----2
                |
    *     *    M2----1
                |
    B1---B2----B3----0

    3-0 = 3 / 3 elemnts = 1
    Parameters
    ----------
    nids_top
    nids_mid
    nids_bottom

    Returns
    -------

    """
    if nid0 == 0:
        nid0 = max(model.nodes) + 1
    assert len(nids_top) == len(nids_bottom), (nids_top, nids_bottom)

    ncols = len(nids_top)
    nrows = len(nids_mid) + 2
    nids_bottom2 = nids_bottom[::-1]
    A = np.zeros((nrows, ncols), dtype='int32')
    A[0, :] = nids_top
    print(A)
    print(nids_mid)
    A[1:1+len(nids_mid), -1] = nids_mid
    A[-1, :] = nids_bottom2
    print(A)

    # assume linear
    nnodes_mid = len(nids_mid)
    nelements_tall = nnodes_mid + 1
    for jcol in range(ncols-1):
        nid_top = nids_top[jcol]
        nid_btm = nids_bottom2[jcol]
        print('*', nid_top, nid_btm)
        xyz_top = model.nodes[nid_top].get_position()
        xyz_btm = model.nodes[nid_btm].get_position()
        dz = (xyz_top - xyz_btm) / nelements_tall

        for inode in range(nnodes_mid):
            xyz_mid = xyz_top + dz * (inode + 1)
            model.add_grid(nid0, xyz=xyz_mid)
            nid0 += 1
            A[inode+1, jcol] = nid0
    eid0 = nid_matrix_to_quads(A, eid0, pid)


def nid_matrix_to_quads(nids_array: np.ndarray, eid0: int, pid: int) -> int:
    print(nids_array)
    n1 = nids_array[:-1, :-1].ravel()
    n2 = nids_array[:-1, 1:].ravel()
    n3 = nids_array[1:, 1:].ravel()
    n4 = nids_array[1:, :-1].ravel()

    for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
        print(n1i, n2i, n3i, n4i)
        nodes = [n1i, n2i, n3i, n4i]
        model.add_cquad4(eid0, pid, nodes)
        eid0 += 1
    return eid0
    #print(n1)
    #print(n2)
    #print(n3)
    #print(n4)

def main():
    # A = np.array([
    #     [1, 2, 3],
    #     [4, 5, 6],
    #     [7, 8, 9],
    # ])
    # nid_matrix_to_quads(A, 0)
    # asdf
    model = BDF()
    nids_top = [1, 2, 3]
    nids_mid = [4, 5]
    nids_bottom = [6, 7, 8]
    nid_xyzs = [
        (1, 0., 0., 0.),
        (2, 1., 0., 0.),
        (3, 2., 0., 0.),
        # mid
        (4, 2., 1., 0.),
        (5, 2., 2., 0.),
        # btm
        (6, 2., 3., 0.),
        (7, 1., 0., 0.),
        (8, 0., 0., 0.),
    ]
    for nid, x, y, z in nid_xyzs:
        xyz = [x, y, z]
        model.add_grid(nid, xyz=xyz)
    build_rib(model, nids_top, nids_mid, nids_bottom)
    nids = [(6001145, 6001146, 6011896),
            (6011897, 6011903, 6011908, 6011960),
            (6011963, 6011968, 6011971),
            ]


if __name__ == '__main__':  # pragma: no cover
    main()
