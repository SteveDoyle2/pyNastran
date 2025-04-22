import numpy as np

from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.cards.elements.solid import CTETRA

def fill_volume(model: BDF,
                mass=100.,  # lb
                density=0.0294,  # lb/in^3; 6.8 lb/gal
                ) -> tuple[float, float]:

    # add a fuel centroid. we'll change the position later
    grid = model.grid
    nid_centroid = max(grid.node_id) + 1
    xyz_centroid = np.zeros(3)
    model.add_grid(nid_centroid, xyz_centroid)
    model.setup()

    nids = grid.node_id
    xyz_cid0 = grid.xyz_cid0()
    z_fill, mass_fill = get_mass_z_level(
        model, nids, xyz_cid0,
        mass, density)

    ctria3 = model.ctria3
    nfaces = len(ctria3)
    normal = ctria3.normal()
    nz = normal[:, 2]
    inz = np.where(nz > 0)[0]
    print(f'inz = {inz}')
    n_inz = len(inz)

    area = ctria3.area()[inz]
    centroid = ctria3.centroid()[inz, :]
    xyz_fill = centroid.copy()
    xyz_fill[:, 2] = z_fill
    assert len(area) == n_inz, (len(area), n_inz)
    assert len(centroid) == n_inz, (centroid.shape, n_inz)

    print(f'xyz_fill.shape = {xyz_fill.shape}')
    dxyz = xyz_fill - centroid
    print(f'dxyz.shape = {dxyz.shape}')
    height = np.linalg.norm(dxyz, axis=1)

    print(f'n_inz = {len(inz)}')
    print(f'area = {area}; n={len(area)}')
    print(f'height = {height}; n={len(height)}')
    assert len(height) == len(dxyz), (height.shape, dxyz.shape)
    assert len(dxyz) == len(xyz_fill)
    volume = area * height
    mass_fill_array = volume * density

    eid = 1
    nid = 1
    for ieid, massi in enumerate(mass_fill_array):
        if massi > 0:
            continue
        model.add_conm2(eid, nid, massi)
        eid += 1
        nid += 1

    z = xyz_cid0[:, 2]
    z1 = z.min()
    z2 = z.max()
    return z_fill, mass_fill

def get_mass_z_level(model: BDF,
                     nids: np.ndarray,
                     xyz_cid0: np.ndarray,
                     mass: float,
                     density: float) -> tuple[float, float]:
    # assume there already are some tris

    z = xyz_cid0[:, 2]
    z1 = z.min()
    z2 = z.max()

    #-------------------------------------
    # convert quads to tris
    cquad4 = model.cquad4
    eid0 = max(model.shell_element_ids) + 1
    if len(cquad4):
        pid = cquad4.property_id[0]
        quads = model.cquad4.nodes
        tris = cquad4.split_to_ctria3()
        for ieid, tri in enumerate(tris):
            eid = eid0 + ieid
            model.add_ctria3(eid, pid, tri)
    cquad4.clear()

    #-------------------------------------
    # convert tris to tets
    model.setup()
    tris = model.ctria3.nodes
    assert len(tris) > 0, tris.shape

    itris = np.searchsorted(nids, tris)

    # take the average triangle location as the "centroid"
    # z value is what is important; ignore it for now
    xyz_centroids = model.ctria3.centroid()
    xyz_centroid = xyz_centroids.mean(axis=0)
    assert len(xyz_centroid) == 3

    #
    mass1 = 0.  # get_mass(tris, xyz_cid0, xyz_centroid,
                           #z1, density)
    mass2 = get_mass(model, tris, nids, xyz_cid0, xyz_centroid,
                     z2, density)
    mass_total = mass2

    # placeholder values; fully fueled
    z_interp = z2
    mass_interp = mass_total

    # kind of unused
    tol = 0.1  # %
    nmax = 20
    n = 0
    print(f'{n}: z1={z1:g} mass1={mass1}; z2={z2:g} mass2={mass2}')
    while abs(mass2 - mass1) / mass_total > tol and n < nmax:
        slope = (z2 - z1) / (mass2 - mass1)
        z_interp = slope * (mass - mass1) + z1

        mass_interp = get_mass(model,
                               tris, nids, xyz_cid0, xyz_centroid,
                               z_interp, density)
        dz1 = abs(z_interp - z1)
        dz2 = abs(z_interp - z2)

        if dz1 <= dz2:
            z1 = z_interp
            mass1 = mass_interp
        else:
            z2 = z_interp
            mass2 = mass_interp
        n += 1
        print(f'{n}: z1={z1:g} mass1={mass1:g}; z2={z2:g} mass2={mass2:g}; dz1={dz1:g} dz2={dz2:g}')

    # if abs(z_interp - z_old) > tol:
    #     raise RuntimeError(f'Did not converge; Check your units; n=nmax={nmax}\n'
    #                        f'target alt={alt_final} alt_current={alt1}')
    #-------------------------------------
    return z_interp, mass_interp

def get_mass(model: BDF,
             tris: np.ndarray,
             nids: np.ndarray,
             xyz_cid0: np.ndarray,
             xyz_centroid: np.ndarray,
             z: float, density: float) -> float:
    xyz_centroid[-1] = z
    nelements = len(tris)
    assert len(tris) > 0, tris.shape
    inid_max = len(nids)

    grid = model.grid
    assert len(nids) == len(grid)
    grid.xyz[:, :] = xyz_centroid
    #grid.xyz[:, -1] = z

    tetra_elem = CTETRA(model)
    tetra_elem.element_id = np.arange(1, nelements+1, dtype='int32')
    tetra_elem.ifile = np.zeros(nelements, dtype='int32')
    tetra_elem.property_id = np.ones(nelements, dtype='int32')
    tetra_elem.nodes = np.column_stack([tris, tetra_elem.property_id * inid_max])
    tetra_elem.n = nelements
    #print(tetra_elem.get_stats())
    #model.setup()
    # tetra_elem.set_xyz_cid0(xyz_cid0)

    volume = tetra_elem.volume()
    print(f'density = {density}')
    print(f'volume = {volume}')
    mass = volume * density
    massi = mass.sum()
    assert isinstance(massi, float), massi
    return massi


def main():
    model = BDF()
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])

    model.add_grid(5, [0., 0., 1.])
    model.add_grid(6, [1., 0., 1.])
    model.add_grid(7, [1., 1., 1.])
    model.add_grid(8, [0., 1., 1.])

    model.add_chexa(1, 1, [1, 2, 3, 4, 5, 6, 7, 8])
    model.add_psolid(1, 1)
    model.add_mat1(1, 3.0e7, None, 0.3)
    model.add_pshell(2, 1, 0.1)
    model.setup()

    chexa = model.chexa
    faces = chexa.faces()
    nfaces = len(faces)

    eid = 2

    cquad4 = model.cquad4
    cquad4.element_id = np.arange(1, nfaces+1, dtype='int32')
    cquad4.property_id = cquad4.element_id * 2
    cquad4.nodes = faces

    cquad4.ifile = np.zeros(nfaces, dtype='int32')
    cquad4.mcid = np.ones(nfaces, dtype='int32') * -1
    cquad4.theta = np.zeros(nfaces, dtype='float64')
    cquad4.zoffset = np.zeros(nfaces, dtype='float64')
    cquad4.tflag = np.zeros(nfaces, dtype='int32')
    cquad4.T = np.zeros((nfaces, 4), dtype='float64')

    cquad4.n = nfaces
    #print(cquad4.get_stats())
    #print(cquad4)
    # cquad4.validate()

    model.setup()

    # volume = 1
    # density = 1
    # mass = volume * density
    z_fill, mass_fill = fill_volume(model, mass=0.5, density=1.0)
    print(f'z_fill={z_fill:g} mass_fill={mass_fill:g}')

if __name__ == '__main__':
    main()
