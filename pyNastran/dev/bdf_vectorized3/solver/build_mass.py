import numpy as np
from cpylog import SimpleLogger
from pyNastran.nptyping_interface import (
    NDArrayNbool,
    NDArrayNint, NDArrayN2int,
    NDArrayNfloat, NDArrayNNfloat)

from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from pyNastran.dev.bdf_vectorized3.bdf_interface.breakdowns import NO_MASS

from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    consistent_mass, lumped_mass, beam_transform,)
from .elements.solids import build_mbb_solids

from .build_stiffness import _COOAccumulator
from .utils import DOF_MAP


def build_Mbb(model: BDF, subcase: Subcase,
              dof_map: DOF_MAP, ndof: int,
              fdtype: str="float64") -> NDArrayNNfloat:
    """builds the mass matrix in the basic frame, [Mbb]"""
    log = model.log
    log.info("starting build_Mbb")
    wtmass = 1.0
    coo_m = _COOAccumulator(ndof)
    str(model)
    str(subcase)
    # no_mass = {
    #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
    # }

    mass_rod_2x2 = 1/3 * np.array([
        [2, 0, 1, 0],
        [0, 2, 0, 1],
        [1, 0, 2, 0],
        [0, 1, 0, 2],
    ], dtype="float64")
    mass_tri = 1/6 * np.array([
        [4, 0, 2, 0, 1, 0],
        [0, 4, 0, 2, 0, 1],
        [2, 0, 4, 0, 2, 0],
        [0, 2, 0, 4, 0, 2],
        [1, 0, 2, 0, 4, 0],
        [0, 1, 0, 2, 0, 4],
    ], dtype="float64")
    # mass_quad_1x1 = np.array([
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # ], dtype='float64') / 36.

    mass_quad_2x2 = 1/36 * np.array([
        [4, 0, 2, 0, 1, 0, 2, 0],
        [0, 4, 0, 2, 0, 1, 0, 2],
        [2, 0, 4, 0, 2, 0, 1, 0],
        [0, 2, 0, 4, 0, 2, 0, 1],
        [1, 0, 2, 0, 4, 0, 2, 0],
        [0, 1, 0, 2, 0, 4, 0, 2],
        [2, 0, 1, 0, 2, 0, 4, 0],
        [0, 2, 0, 1, 0, 2, 0, 4],
    ], dtype="float64")

    mass_total = 0.0
    if model.conm1:
        nid = elem.nid
        nid_ref = elem.nid_ref
        if nid_ref.type == "GRID":
            i1 = dof_map[(nid, 1)]

            # TODO: support CID
            if nid_ref.cd != elem.cid:
                log.warning(
                    f"  CONM1 eid={eid} nid={nid} CD={nid_ref.cd} to cid={elem.cid} is not supported")
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)
        coo_m.add_matrix(list(range(i1, i1 + 6)), elem.mass_matrix)

    # PARAM,COUPMASS: 1 => consistent mass, -1 or absent => lumped
    use_consistent_shells = False
    if hasattr(model, "params") and "COUPMASS" in model.params:
        coupmass_val = model.params["COUPMASS"].values[0]
        if coupmass_val >= 1:
            use_consistent_shells = True

    # Consistent mass matrices for shells (m/36 * [4,2,1,2;...])
    mass_quad_consistent = 1/36 * np.array([
        [4, 2, 1, 2],
        [2, 4, 2, 1],
        [1, 2, 4, 2],
        [2, 1, 2, 4],
    ], dtype="float64")
    mass_tri_consistent = 1/12 * np.array([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ], dtype="float64")
    
    # has possibility of mass
    has_mass = False

    for elem in model.element_cards:
        etype = elem.type
        if etype in NO_MASS:
            continue
        if elem.n == 0:
            continue

        has_mass = True
        if etype == "CONM2":
            mass_total = conm2_fill_Mbb(model, mass_total, coo_m, dof_map)

        elif etype in ["CROD", "CONROD", "CTUBE"]:
            # verified
            mass = elem.mass()
            if mass == 0.0:
                log.warning(f"  no mass for {etype} eid={eid}")
                continue

            nids1 = elem.nodes[:, 0]
            nids2 = elem.nodes[:, 1]
            for nid1, nid2 in zip(nids1, nids2):
                i1 = dof_map[(nid1, 1)]
                j1 = dof_map[(nid2, 1)]
                ii = [i1, i1 + 1, j1, j1 + 1]
                coo_m.add_matrix(ii, mass_rod_2x2 * mass)
        elif etype in {"CBAR", "CBEAM"}:
            # PARAM,COUPMASS,1 => consistent; default (0 or absent) => lumped
            use_consistent = False
            if hasattr(model, 'params') and 'COUPMASS' in model.params:
                coupmass_val = model.params['COUPMASS'].values[0]
                if coupmass_val >= 1:
                    use_consistent = True

            area = elem.area()
            inertia = elem.inertia()
            xyz1, xyz2 = elem.get_xyz()
            lengths = np.linalg.norm(xyz2 - xyz1, axis=1)
            v, ihat_arr, yhat_arr, zhat_arr, wa_arr, wb_arr = elem.get_axes(xyz1, xyz2)
            k_arr = elem.k()
            e_g_nus = elem.e_g_nu()

            elem_masses = elem.mass()
            mass_total += elem_masses.sum()
            mass_per_length_total = elem_masses / lengths

            for (nid1, nid2), areai, inertiai, lengthi, ki, e_g_nu, ihati, jhati, khati, mpl in zip(
                elem.nodes, area, inertia, lengths, k_arr, e_g_nus,
                ihat_arr, yhat_arr, zhat_arr, mass_per_length_total):

                i1_inertia, i2_inertia, i12, j = inertiai
                e, g, nu = e_g_nu
                k1, k2 = ki

                rho_eff = mpl / areai if areai > 0 else 0.0
                if use_consistent:
                    Me = consistent_mass(
                        areai, lengthi, rho_eff,
                        i1_inertia, i2_inertia, j,
                        k1, k2, nsm=0.0)
                else:
                    Me = lumped_mass(
                        areai, lengthi, rho_eff,
                        i1_inertia, i2_inertia, j,
                        nsm=0.0)

                Teb = beam_transform(ihati, jhati, khati)
                M_basic = Teb.T @ Me @ Teb

                gi1 = dof_map[(nid1, 1)]
                gi2 = dof_map[(nid2, 1)]
                n_ijv = [
                    gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
                    gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
                ]
                coo_m.add_matrix(n_ijv, M_basic)
        elif etype == "CTRIA3":
            masses = elem.mass()
            if masses.sum() == 0.0:
                log.warning(f"  no mass for CTRIA3 eid={elem.element_id}")
                continue
            for (nid1, nid2, nid3), massi in zip(elem.nodes, masses):
                if massi == 0.0:
                    continue
                if use_consistent_shells:
                    # Consistent: m * [2,1,1;1,2,1;1,1,2]/12
                    nids_e = [nid1, nid2, nid3]
                    M_consist = mass_tri_consistent * massi
                    for dof_offset in range(3):
                        ii = [dof_map[(n, 1)] + dof_offset for n in nids_e]
                        coo_m.add_matrix(ii, M_consist)
                else:
                    # Lumped: m/3 per node on Tx, Ty, Tz
                    m_node = massi / 3.0
                    for nid in [nid1, nid2, nid3]:
                        i1 = dof_map[(nid, 1)]
                        coo_m.add_scalar(i1, i1, m_node)
                        coo_m.add_scalar(i1 + 1, i1 + 1, m_node)
                        coo_m.add_scalar(i1 + 2, i1 + 2, m_node)
            # Mbb[i1, i1] = Mbb[i1+1, i1+1] = Mbb[i1+2, i1+2] = \
            # Mbb[i2, i2] = Mbb[i2+1, i2+1] = Mbb[i2+2, i2+2] = \
            # Mbb[i3, i3] = Mbb[i3+1, i3+1] = Mbb[i3+2, i3+2] = mass / 3
        elif etype == "CQUAD4":
            masses = elem.mass()
            if masses.sum() == 0.0:
                log.warning(f"  no mass for CQUAD4 eid={elem.element_id}")
                continue
            for (nid1, nid2, nid3, nid4), massi in zip(elem.nodes, masses):
                if massi == 0.0:
                    continue
                if use_consistent_shells:
                    # Consistent: m * [4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]/36
                    # Applied to each translational DOF (Tx, Ty, Tz) independently
                    nids_e = [nid1, nid2, nid3, nid4]
                    M_consist = mass_quad_consistent * massi
                    for dof_offset in range(3):
                        ii = [dof_map[(n, 1)] + dof_offset for n in nids_e]
                        coo_m.add_matrix(ii, M_consist)
                else:
                    # Lumped: m/4 per node on Tx, Ty, Tz
                    m_node = massi / 4.0
                    for nid in [nid1, nid2, nid3, nid4]:
                        i1 = dof_map[(nid, 1)]
                        coo_m.add_scalar(i1, i1, m_node)
                        coo_m.add_scalar(i1 + 1, i1 + 1, m_node)
                        coo_m.add_scalar(i1 + 2, i1 + 2, m_node)
            # if 0:  # pragma: no cover
            # mass4 = mass / 9. # 4/36
            # mass2 = mass / 18. # 2/36
            # mass1 = mass / 36.
            # print(mass1, mass2, mass4)
            # Mbb[i1, i1] += mass4
            # Mbb[i1+1, i1+1] += mass4
            # Mbb[i2, i2] += mass4
            # Mbb[i2+1, i2+1] += mass4
            # Mbb[i3, i3] += mass4
            # Mbb[i3+1, i3+1] += mass4
            # Mbb[i4, i4] += mass4
            # Mbb[i4+1, i4+1] += mass4

            # Mbb[i1, i3] += mass1
            # Mbb[i1, i2] += mass2
            # Mbb[i1, i4] += mass2
            # Mbb[i3, i1] += mass1
            # Mbb[i2, i1] += mass2
            # Mbb[i4, i1] += mass2

            # Mbb[i1+1, i3+1] += mass1
            # Mbb[i1+1, i2+1] += mass1
            # Mbb[i1+1, i4+1] += mass2
            # Mbb[i3+1, i1+1] += mass1
            # Mbb[i2+1, i1+1] += mass2
            # Mbb[i4+1, i1+1] += mass2

            # Mbb[i2, i4] += mass1
            # Mbb[i2+1, i4+1] += mass1
            # Mbb[i2, i3] += mass2
            # Mbb[i2+1, i3+1] += mass2

            # Mbb[i4, i2] += mass1
            # Mbb[i4+1, i2+1] += mass1
            # Mbb[i3, i2] += mass2
            # Mbb[i3+1, i2+1] += mass2

            # Mbb[i3, i4] += mass2
            # Mbb[i3+1, i4+1] = mass2
            # Mbb[i4, i3] += mass2
            # Mbb[i4+1, i3+1] += mass2

            # Mbb[i2+1, i3+1] = 1
            # Mbb[i2+1, i2+1] = Mbb[i1+1, i4+1] = 2

            # Mbb[i2, i1+1] = Mbb[i3, i1+1] = Mbb[i4, i1+1] = 1
            # Mbb[i1, i2+1] = Mbb[i3, i2+1] = Mbb[i4, i2+1] = 1
            # Mbb[i1, i3+1] = Mbb[i2, i3+1] = Mbb[i4, i3+1] = 1
            # Mbb[i1, i4+1] = Mbb[i2, i4+1] = Mbb[i3, i4+1] = 1
            # print(Mbb)
            # print(Mbb[ii, :][:, ii])
        elif etype == "CSHEAR":
            masses = elem.mass()
            if masses.sum() == 0.0:
                continue
            for (nid1, nid2, nid3, nid4), massi in zip(elem.nodes, masses):
                if massi == 0.0:
                    continue
                i1 = dof_map[(nid1, 1)]
                i2 = dof_map[(nid2, 1)]
                i3 = dof_map[(nid3, 1)]
                i4 = dof_map[(nid4, 1)]
                ii = [
                    i1,
                    i1 + 1,
                    i2,
                    i2 + 1,
                    i3,
                    i3 + 1,
                    i4,
                    i4 + 1,
                ]
                coo_m.add_matrix(ii, mass_quad_2x2 * massi)
        elif etype in {"CHEXA", "CTETRA", "CPENTA"}:
            pass  # handled below
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)

    # Solid elements (CHEXA, CTETRA, CPENTA) — consistent mass
    mass_total += build_mbb_solids(model, coo_m, dof_map)

    # Convert COO accumulator to sparse CSC
    Mbb = coo_m.to_csc()

    if wtmass != 1.0:
        Mbb *= wtmass

    has_special_points = "SPOINT" in model.card_count or "EPOINT" in model.card_count
    is_all_grids = not has_special_points
    unused_can_dof_slice = is_all_grids and not has_mass
    if Mbb.nnz > 0:
        i = np.arange(0, ndof).reshape(ndof // 6, 6)[:, :3].ravel()
        massi = sum(Mbb[ii, ii] for ii in i)
        log.info(f"finished build_Mbb; M={massi:.6g}; mass_total={mass_total:.6g}")
    else:
        return None
    return Mbb



def conm2_fill_Mbb(model: BDF, mass_total: float,
                   coo_m, dof_map: dict[tuple[int, int], int]) -> float:
    eye3 = np.eye(3, dtype="float64")
    conm2 = model.conm2
    log = model.log
    # mass = elem.Mass()
    # nid = elem.nid
    # nid_ref = elem.nid_ref
    inid = model.grid.index(conm2.node_id)
    cds = model.grid.cd[inid]
    for eid, nid, cd, cid, mass, elem_x, elem_i in zip(
        conm2.element_id,
        conm2.node_id,
        cds,
        conm2.coord_id,
        conm2.mass(),
        conm2.xyz_offset,
        conm2.inertia,
    ):
        i1 = dof_map[(nid, 1)]
        if cd != cid:
            log.warning(f"  CONM2 eid={eid} nid={nid} CD={cd} to CONM2 cid={cid} is not supported")
        # Mbb[i1, i1] = mass
        # Mbb[i1+1, i1+1] = mass
        # Mbb[i1+2, i1+2] = mass
        # TODO: support CID
        I11, I21, I22, I31, I32, I33 = elem_i
        x1, x2, x3 = elem_x
        mxx = (
            np.array(
                [
                    [x1 * x1, -x1 * x2, -x1 * x3],
                    [-x2 * x1, x2 * x2, -x2 * x3],
                    [-x3 * x1, x3 * x2, x3 * x3],
                ]
            )
            * mass
        )
        Tr = np.array(
            [
                [0, x3, -x2],
                [-x3, 0, x1],
                [x2, -x1, 0],
            ],
            dtype="float64",
        )
        mx = Tr * mass
        I = (
            np.array(
                [
                    [I11, -I21, I31],
                    [-I21, I22, -I32],
                    [-I31, -I32, I33],
                ]
            )
            + mxx
        )

        # [mass, 01, 02, 03, mass * X3, -mass * X2]
        # [10, mass, 12, -mass * X3, 14, mass * X1]
        # [20, 21, mass, mass * X2, -mass * X1, 25]
        # [30, -mass * X3, mass * X2,        I11 + mass * X2 * X2 + mass * X3 * X3, -I21 - mass * X2 * X1,                  -I31 - mass * X3 * X1]
        # [mass * X3, 41, -mass * X1,       -I21 - mass * X2 * X1,                   I22 + mass * X1 * X1 + mass * X3 * X3, -I32 - mass * X3 * X2]
        # [-mass * X2, mass * X1, 52,       -I31 - mass * X3 * X1,                  -I32 - mass * X3 * X2,                   I33 + mass * X2 * X2 + mass * X1 * X1]

        M6 = np.zeros((6, 6), dtype='float64')
        M6[:3, :3] = eye3 * mass
        M6[:3, 3:] = mx
        M6[3:, :3] = mx.T
        M6[3:, 3:] = I
        coo_m.add_matrix(list(range(i1, i1 + 6)), M6)
        mass_total += mass
    return mass_total


def build_Mbb_og(model: BDF,
              subcase: Subcase,
              dof_map: DOF_MAP,
              ndof: int, fdtype='float64') -> NDArrayNNfloat:
    """builds the mass matrix in the basic frame, [Mbb]"""
    log = model.log
    log.info('starting build_Mbb')
    wtmass = model.get_param('WTMASS', 1.0)
    #Mbb = np.eye(ndof, dtype=fdtype)
    Mbb = np.zeros((ndof, ndof), dtype=fdtype)
    #Mbb = np.eye(ndof, dtype='int32')
    str(model)
    str(subcase)
    no_mass = {
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
    }

    eye3 = np.eye(3, dtype='float64')
    mass_rod_2x2 = np.array([
        [2, 0, 1, 0],
        [0, 2, 0, 1],
        [1, 0, 2, 0],
        [0, 1, 0, 2],
    ], dtype='float64') / 3.
    mass_tri = np.array([
        [4, 0, 2, 0, 1, 0],
        [0, 4, 0, 2, 0, 1],
        [2, 0, 4, 0, 2, 0],
        [0, 2, 0, 4, 0, 2],
        [1, 0, 2, 0, 4, 0],
        [0, 1, 0, 2, 0, 4],
    ], dtype='float64') / 6.
    #mass_quad_1x1 = np.array([
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
    #], dtype='float64') / 36.

    mass_quad_2x2 = np.array([
        [4, 0, 2, 0, 1, 0, 2, 0],
        [0, 4, 0, 2, 0, 1, 0, 2],
        [2, 0, 4, 0, 2, 0, 1, 0],
        [0, 2, 0, 4, 0, 2, 0, 1],
        [1, 0, 2, 0, 4, 0, 2, 0],
        [0, 1, 0, 2, 0, 4, 0, 2],
        [2, 0, 1, 0, 2, 0, 4, 0],
        [0, 2, 0, 1, 0, 2, 0, 4],
    ], dtype='float64') / 36.

    mass_total = 0.
    for eid, elem in model.masses.items():
        etype = elem.type
        if etype in no_mass:
            continue

        if etype == 'CONM1':
            nid = elem.nid
            nid_ref = elem.nid_ref
            if nid_ref.type == 'GRID':
                i1 = dof_map[(nid, 1)]

                # TODO: support CID
                if nid_ref.cd != elem.cid:
                    log.warning(f'  CONM1 eid={eid} nid={nid} CD={nid_ref.cd} to cid={elem.cid} is not supported')
            else:  # pragma: no cover
                print(elem.get_stats())
                raise NotImplementedError(elem)
            Mbb[i1:i1+6, i1:i1+6] += elem.mass_matrix

        if etype == 'CONM2':
            mass = elem.Mass()
            nid = elem.nid
            nid_ref = elem.nid_ref
            if nid_ref.type == 'GRID':
                i1 = dof_map[(nid, 1)]
                if nid_ref.cd != elem.cid:
                    log.warning(f'  CONM2 eid={eid} nid={nid} CD={nid_ref.cd} to cid={elem.cid} is not supported')
                #Mbb[i1, i1] = mass
                #Mbb[i1+1, i1+1] = mass
                #Mbb[i1+2, i1+2] = mass
                # TODO: support CID
                I11, I21, I22, I31, I32, I33 = elem.I
                x1, x2, x3 = elem.X
                mxx = np.array([
                    [x1 * x1, -x1 * x2, -x1 * x3],
                    [-x2 * x1, x2 * x2, -x2 * x3],
                    [-x3 * x1, x3 * x2, x3 * x3],
                ]) * mass
                Tr = np.array([
                    [0, x3, -x2],
                    [-x3, 0, x1],
                    [x2, -x1, 0],
                ], dtype='float64')
                mx = Tr * mass
                I = np.array([
                    [I11, -I21, I31],
                    [-I21, I22, -I32],
                    [-I31, -I32, I33],
                ]) + mxx

                #[mass, 01, 02, 03, mass * X3, -mass * X2]
                #[10, mass, 12, -mass * X3, 14, mass * X1]
                #[20, 21, mass, mass * X2, -mass * X1, 25]
                #[30, -mass * X3, mass * X2,        I11 + mass * X2 * X2 + mass * X3 * X3, -I21 - mass * X2 * X1,                  -I31 - mass * X3 * X1]
                #[mass * X3, 41, -mass * X1,       -I21 - mass * X2 * X1,                   I22 + mass * X1 * X1 + mass * X3 * X3, -I32 - mass * X3 * X2]
                #[-mass * X2, mass * X1, 52,       -I31 - mass * X3 * X1,                  -I32 - mass * X3 * X2,                   I33 + mass * X2 * X2 + mass * X1 * X1]

                Mbb[i1:i1+3, i1:i1+3] += eye3 * mass
                Mbb[i1:i1+3, i1+3:i1+6] += mx
                Mbb[i1+3:i1+6, i1:i1+3] += mx.T
                mass_total += mass
                Mbb[i1+3:i1+6, i1+3:i1+6] += I
            else:  # pragma: no cover
                print(elem.get_stats())
                raise NotImplementedError(elem)
            del i1
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)

    # has possibility of mass
    has_mass = False
    for eid, elem in model.elements.items():
        etype = elem.type
        if etype in no_mass:
            continue

        has_mass = True
        if etype in ['CROD', 'CONROD', 'CTUBE']:
            # verified
            mass = elem.Mass()
            if mass == 0.0:
                log.warning(f'  no mass for {etype} eid={eid}')
                continue

            nid1, nid2 = elem.nodes
            i1 = dof_map[(nid1, 1)]
            j1 = dof_map[(nid2, 1)]
            ii = [i1, i1 + 1,
                  j1, j1 + 1]
            #Mbb[i1, i1] = Mbb[i1+1, i1+1] = \
            #Mbb[j1, j1] = Mbb[j1+1, j1+1] = mass / 3

            #Mbb[i1, j1] = Mbb[j1, i1] = \
            #Mbb[i1+1, j1+1] = Mbb[j1+1, i1+1] = mass / 6
            iii, jjj = np.meshgrid(ii, ii)
            Mbb[iii, jjj] += mass_rod_2x2 * mass
            #ii = [i1, i1 + 1, j1, j1 + 1]
            #print(Mbb[ii, :][:, ii])
        elif etype in ['CBAR', 'CBEAM']:
            # TODO: verify
            # TODO: add rotary inertia
            mass = elem.Mass()
            if mass == 0.0:
                log.warning(f'  no mass for {etype} eid={eid}')
                continue
            nid1, nid2 = elem.nodes
            i1 = dof_map[(nid1, 1)]
            j1 = dof_map[(nid2, 1)]

            #Mbb[i1, i1] = Mbb[i1+1, i1+1] = \
            #Mbb[j1, j1] = Mbb[j1+1, j1+1] = mass / 3

            #Mbb[i1, j1] = Mbb[j1, i1] = \
            #Mbb[i1+1, j1+1] = Mbb[j1+1, i1+1] = mass / 6
            ii = [i1, i1 + 1,
                  j1, j1 + 1]
            iii, jjj = np.meshgrid(ii, ii)
            Mbb[iii, jjj] += mass_rod_2x2 * mass
        elif etype == 'CTRIA3':
            # TODO: verify
            # TODO: add rotary inertia
            mass = elem.Mass()
            nid1, nid2, nid3 = elem.nodes
            i1 = dof_map[(nid1, 1)]
            i2 = dof_map[(nid2, 1)]
            i3 = dof_map[(nid3, 1)]
            ii = [
                i1, i1 + 1,
                i2, i2 + 1,
                i3, i3 + 1,
            ]
            iii, jjj = np.meshgrid(ii, ii)
            Mbb[iii, jjj] += mass_tri * mass
            #Mbb[i1, i1] = Mbb[i1+1, i1+1] = Mbb[i1+2, i1+2] = \
            #Mbb[i2, i2] = Mbb[i2+1, i2+1] = Mbb[i2+2, i2+2] = \
            #Mbb[i3, i3] = Mbb[i3+1, i3+1] = Mbb[i3+2, i3+2] = mass / 3
        elif etype == 'CQUAD4':
            # TODO: verify
            # TODO: add rotary inertia
            try:
                mass = elem.Mass()
            except AttributeError:
                pid_ref = elem.pid_ref
                if pid_ref.mid1 is None and pid_ref.mid2 is None:
                    log.warning(f'  no mass for CQUAD4 eid={eid}')
                    continue
                #raise
                #mid_ref = elem.mid_ref
                #rho = mid_ref.Rho()
            nid1, nid2, nid3, nid4 = elem.nodes
            i1 = dof_map[(nid1, 1)]
            i2 = dof_map[(nid2, 1)]
            i3 = dof_map[(nid3, 1)]
            i4 = dof_map[(nid4, 1)]
            if mass == 0.:
                pid_ref = elem.pid_ref
                ptype = pid_ref.type
                if ptype == 'PSHELL':
                    mid_ref = pid_ref.mid_ref
                    rho = mid_ref.rho
                    log.warning(f'  no mass for CQUAD4 eid={eid} ptype={ptype} rho={rho}')
                else:
                    log.warning(f'  no mass for CQUAD4 eid={eid} ptype={ptype}')
            else:
                pid_ref = elem.pid_ref
                ptype = pid_ref.type
                log.info(f'  mass={mass} for CQUAD4 eid={eid} ptype={ptype}')

            ii = [
                i1, i1 + 1,
                i2, i2 + 1,
                i3, i3 + 1,
                i4, i4 + 1,
            ]
            iii, jjj = np.meshgrid(ii, ii)
            Mbb[iii, jjj] += mass_quad_2x2 * mass
            #if 0:  # pragma: no cover
                #mass4 = mass / 9. # 4/36
                #mass2 = mass / 18. # 2/36
                #mass1 = mass / 36.
                #print(mass1, mass2, mass4)
                #Mbb[i1, i1] += mass4
                #Mbb[i1+1, i1+1] += mass4
                #Mbb[i2, i2] += mass4
                #Mbb[i2+1, i2+1] += mass4
                #Mbb[i3, i3] += mass4
                #Mbb[i3+1, i3+1] += mass4
                #Mbb[i4, i4] += mass4
                #Mbb[i4+1, i4+1] += mass4

                #Mbb[i1, i3] += mass1
                #Mbb[i1, i2] += mass2
                #Mbb[i1, i4] += mass2
                #Mbb[i3, i1] += mass1
                #Mbb[i2, i1] += mass2
                #Mbb[i4, i1] += mass2

                #Mbb[i1+1, i3+1] += mass1
                #Mbb[i1+1, i2+1] += mass1
                #Mbb[i1+1, i4+1] += mass2
                #Mbb[i3+1, i1+1] += mass1
                #Mbb[i2+1, i1+1] += mass2
                #Mbb[i4+1, i1+1] += mass2

                #Mbb[i2, i4] += mass1
                #Mbb[i2+1, i4+1] += mass1
                #Mbb[i2, i3] += mass2
                #Mbb[i2+1, i3+1] += mass2

                #Mbb[i4, i2] += mass1
                #Mbb[i4+1, i2+1] += mass1
                #Mbb[i3, i2] += mass2
                #Mbb[i3+1, i2+1] += mass2

                #Mbb[i3, i4] += mass2
                #Mbb[i3+1, i4+1] = mass2
                #Mbb[i4, i3] += mass2
                #Mbb[i4+1, i3+1] += mass2

                #Mbb[i2+1, i3+1] = 1
                #Mbb[i2+1, i2+1] = Mbb[i1+1, i4+1] = 2

                #Mbb[i2, i1+1] = Mbb[i3, i1+1] = Mbb[i4, i1+1] = 1
                #Mbb[i1, i2+1] = Mbb[i3, i2+1] = Mbb[i4, i2+1] = 1
                #Mbb[i1, i3+1] = Mbb[i2, i3+1] = Mbb[i4, i3+1] = 1
                #Mbb[i1, i4+1] = Mbb[i2, i4+1] = Mbb[i3, i4+1] = 1
                #print(Mbb)
            #print(Mbb[ii, :][:, ii])
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)

    if wtmass != 1.0:
        Mbb *= wtmass

    # TODO: working on a hack to NOTt fake the mass when there are no SPOINTs
    #       as mass can be easily added to a spring model
    #
    # TODO: non-zero mass case doesn't handle SPOINTs
    #
    has_special_points = 'SPOINT' in model.card_count or 'EPOINT' in model.card_count
    is_all_grids = not has_special_points
    unused_can_dof_slice = is_all_grids and not has_mass
    if Mbb.sum() != 0.0:
    #if Mbb.sum() != 0.0 or can_dof_slice:
        #print(f'is_all_grids={is_all_grids} has_mass={has_mass}; can_dof_slice={can_dof_slice} Mbb.shape={Mbb.shape}')
        i = np.arange(0, ndof).reshape(ndof//6, 6)[:, :3].ravel()
        #print(Mbb[i, i])
        massi = Mbb[i, i].sum()
        log.info(f'finished build_Mbb; M={massi:.6g}; mass_total={mass_total:.6g}')
    else:
        Mbb = np.eye(ndof, dtype=fdtype)
        log.error(f'finished build_Mbb; faking mass; M={Mbb.sum()} ndof={ndof}')
    return Mbb
