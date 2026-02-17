from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d
from pyNastran.dev.bdf_vectorized3.solver.utils import get_ieids_eids, get_element
from pyNastran.op2.op2_interface.op2_classes import (
    RealSpringForceArray, RealRodForceArray, RealCBarForceArray,
)
# from pyNastran.dev.solver.build_stiffness import ke_cbar
from .static_spring import _recover_force_celas
from .utils import get_plot_request

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, Subcase, CBAR, PBAR, PBARL
    DOF_MAP = dict[tuple[int, int], int]


def recover_force_101(f06_file: TextIO, op2: OP2,
                      model: BDF, dof_map: DOF_MAP,
                      subcase: Subcase, xb: np.ndarray, fdtype: str='float32',
                      title: str='', subtitle: str='', label: str='',
                      page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - FORCE = ALL

    """
    if 'FORCE' not in subcase:
        asdf
        return

    eid_str = 'ALL'
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, 'FORCE')
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS3', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS4', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)

    nelements += _recover_force_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CONROD', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_rod(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CTUBE', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_force_cbar(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBAR', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no force output...{model.card_count}; {model.bdf_filename}')
        # raise RuntimeError(f'nelements={nelements} cards={model.card_count}')


def _recover_force_rod(f06_file: TextIO, op2: OP2,
                       model: BDF, dof_map: DOF_MAP, isubcase: int,
                       xb: np.ndarray, eids_str,
                       element_name: str, fdtype='float32',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp='PAGE %s') -> int:
    """recovers static rod force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    forces = np.full((neids, 2), np.nan, dtype=fdtype)

    if element_name == 'CONROD':
        prop = elem
    elif element_name == 'CROD':
        pid = elem.property_id
        prop = model.prod.slice_card_by_property_id(pid)
    elif element_name == 'CTUBE':
        pid = elem.property_id
        prop = model.ptube.slice_card_by_property_id(pid)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    mat1 = model.mat1.slice_card_by_material_id(prop.material_id)
    E = mat1.E
    G = mat1.G
    if prop.type == 'PTUBE':
        A = prop.area()
        J = prop.J()
    else:
        A = prop.A
        J = prop.J

    assert isinstance(J, np.ndarray), (prop.type, J)
    assert isinstance(A, np.ndarray), (prop.type, A)
    assert isinstance(E, np.ndarray), (prop.type, E)
    for ieid, eid, nodes, xyz1i, xyz2i, Gi, Ji, Ai, Ei in zip(
            ieids, eids, elem.nodes, xyz1, xyz2, G, J, A, E):
        forces[ieid, :] = _recover_forcei_rod(
            xb, dof_map, nodes, xyz1i, xyz2i, Gi, Ji, Ai, Ei)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealRodForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    if element_name == 'CONROD':
        force.conrod_force[isubcase] = force_obj
    elif element_name == 'CROD':
        force.crod_force[isubcase] = force_obj
    elif element_name == 'CTUBE':
        force.ctube_force[isubcase] = force_obj
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return neids

def _recover_forcei_rod(xb: np.ndarray,
                        dof_map: DOF_MAP,
                        nodes: np.ndarray,
                        xyz1: np.ndarray, xyz2: np.ndarray,
                        G, J, A, E):
    """get the static rod force"""
    nid1, nid2 = nodes
    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_axial = np.array([
        xb[i1], xb[i1+1], xb[i1+2],
        xb[i2], xb[i2+1], xb[i2+2]
    ])
    q_torsion = np.array([
        xb[i1+3], xb[i1+4], xb[i1+5],
        xb[i2+3], xb[i2+4], xb[i2+5]
    ])
    dxyz12 = xyz1 - xyz2
    Lambda = lambda1d(dxyz12, debug=False)

    u_axial = Lambda @ q_axial
    u_torsion = Lambda @ q_torsion
    du_axial = u_axial[0] - u_axial[1]
    du_torsion = u_torsion[0] - u_torsion[1]
    #headers = ['axial', 'SMa', 'torsion', 'SMt']
    #C = prop.c
    L = np.linalg.norm(dxyz12)
    axial_strain = du_axial / L
    #torsional_strain = du_torsion * C / L

    axial_stress = E * axial_strain
    #torsional_stress = G * torsional_strain
    axial_force = axial_stress * A
    torsional_moment = du_torsion * G * J / L
    return axial_force, torsional_moment

def _recover_force_cbar(f06_file: TextIO, op2,
                        model: BDF, dof_map, isubcase, xb, eids_str,
                        element_name, fdtype='float32',
                        title: str='', subtitle: str='', label: str='',
                        page_num: int=1, page_stamp='PAGE %s') -> int:
    """
    Recovers static CBAR force.

    .. todo:: doesn't support CBAR-100

    """
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    xyz1 = model.grid.get_position_by_node_id(elem.nodes[:, 0])
    xyz2 = model.grid.get_position_by_node_id(elem.nodes[:, 1])

    # mat1 = model.mat1.slice_card_by_material_id(prop.material_id)
    # E = mat1.E
    # G = mat1.G
    # if prop.type == 'PTUBE':
    #     A = prop.area()
    #     J = prop.J()
    # else:
    #     A = prop.A
    #     J = prop.J

    AIJEG = elem.stiffness_info()
    A = AIJEG[:, 0]
    I = AIJEG[:, [1, 2, 3]]
    J = AIJEG[:, 4]
    E = AIJEG[:, 5]
    G = AIJEG[:, 6]
    assert isinstance(J, np.ndarray), (elem.type, J)
    assert isinstance(A, np.ndarray), (elem.type, A)
    assert isinstance(E, np.ndarray), (elem.type, E)

    forces = np.full((neids, 8), np.nan, dtype=fdtype)

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
    for (ieid, eid, nodes, xyz1i, xyz2i,
         Ai, I1, Ji, Ei, Gi,
         vi, ihati, yhati, zhati, wai, wbi) in zip(
            ieids, eids, elem.nodes, xyz1, xyz2, A, I, J, E, G,
            v, ihat, yhat, zhat, wa, wb):
        forces[ieid, :] = _recover_forcei_cbar(model, xb, dof_map, nodes,
                                               xyz1i, xyz2i, Ai, I1, Ji, Ei, Gi,
                                               vi, ihati, yhati, zhati, wai, wbi)

    data = forces.reshape(1, *forces.shape)
    table_name = 'OEF1'
    force_obj = RealCBarForceArray.add_static_case(
        table_name, 'CBAR', eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    force.cbar_force[isubcase] = force_obj

    force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
                        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return len(elem)

def ke_cbar(xyz1, xyz2,
            A, I, J, E, G,
            ihat, jhat, khat, wa, wb,
            fdtype: str='float64'):
    """get the elemental stiffness matrix in the basic frame"""
    I1, I2, I12 = I
    pa = pb = 0
    z = np.zeros((3, 3), dtype='float64')
    T = z
    dxyz = xyz2 - xyz1
    L = np.linalg.norm(dxyz)
    T = np.vstack([ihat, jhat, khat])
    z = np.zeros((3, 3), dtype=fdtype)
    Teb = np.block([
        [T, z, z, z],
        [z, T, z, z],
        [z, z, T, z],
        [z, z, z, T],
    ])
    k1 = k2 = 0.
    # k1 = pid_ref.k1
    # k2 = pid_ref.k2
    Ke = _beami_stiffness(L, A,
                          I1, I2, I12, J,
                          E, G,
                          k1=k1, k2=k2, pa=pa, pb=pb)
    K = Teb.T @ Ke @ Teb
    is_passed = True
    return is_passed, K

def _beami_stiffness(L: float, A: float,
                     Iy: float, Iz: float, Iyz: float, J: float,
                     E: float, G: float,
                     pa: int=0, pb: int=0,
                     k1: Optional[float]=None,
                     k2: Optional[float]=None):
    """gets the ith Euler-Bernoulli beam stiffness"""
    kaxial = E * A / L
    ktorsion = G * J / L

    L2 = L * L
    if k1 is not None:
        phiy = 12. * E * Iy / (k1 * G * A * L2)
    if k2 is not None:
        phiz = 12. * E * Iy / (k2 * G * A * L2)

    phiy = 1.0
    phiz = 1.0
    ky = E * Iy / (L * phiy)
    kz = E * Iz / (L * phiz)

    K = np.zeros((12, 12), dtype='float64')
    # axial
    K[0, 0] = K[6, 6] = kaxial
    K[6, 0] = K[0, 6] = -kaxial

    # torsion
    K[3, 3] = K[9, 9] = ktorsion
    K[9, 3] = K[3, 9] = -ktorsion

    #Fx - 0, 6
    #Fy - 1, 7**
    #Fz - 2, 8
    #Mx - 3, 9
    #My - 4, 10
    #Mz - 5, 11**
    # bending z (Fy/1/7, Mz/5/11)
    #      1     5       7   11
    # 1  [12  & 6L   & -12 & 6L
    # 5  [6L  & 4L^2 & -6L & 2L^2
    # 7  [-12 &-6L   &  12 & -6L
    # 11 [6L  & 2L^2 & -6L & 4L^2
    K[1, 1] = K[7, 7] = 12. * kz
    K[1, 7] = K[1, 7] = -12. * kz
    K[1, 5] = K[5, 1] = K[11, 1] = K[1, 11] = 6. * L * kz

    K[5, 7] = K[7, 5] = K[7, 11] = K[11, 7] = -6. * L * kz
    K[5, 11] = K[11, 5] = 2. * L2 * kz * (2 - phiz)
    K[5, 5] = K[11, 11] = 4. * L2 * kz * (4 + phiz)

    #Fx - 0, 6
    #Fy - 1, 7
    #Fz - 2, 8**
    #Mx - 3, 9
    #My - 4, 10**
    #Mz - 5, 11
    # bending y (Fz/2/8, My/4/10)
    #      2     4       8   10
    # 2  [12  & 6L   & -12 & 6L
    # 4  [6L  & 4L^2 & -6L & 2L^2
    # 8  [-12 &-6L   &  12 & -6L
    # 10 [6L  & 2L^2 & -6L & 4L^2
    K[2, 2] = K[8, 8] = 12. * ky
    K[2, 8] = K[2, 8] = -12. * ky
    K[2, 4] = K[4, 2] = K[10, 2] = K[2, 10] = 6. * L * ky

    K[4, 8] = K[8, 4] = K[8, 10] = K[10, 8] = -6. * L * ky
    K[4, 10] = K[10, 4] = 2. * L * L * ky * (2. - phiy)
    K[4, 4] = K[10, 10] = 4. * L * L * ky * (4. + phiy)

    if pa != 0:
        assert pa > 0
        for pas in str(pa): # 123456
            pai = int(pas) - 1 # 012345 (0-5)
            K[pai, :] = 0
            K[:, pai] = 0
    if pb != 0:
        assert pb > 0
        for pbs in str(pb): # 123456
            pbi = int(pbs) + 5 #  (6 - 11)
            K[pbi, :] = 0
            K[:, pbi] = 0
    return K

def _recover_forcei_cbar(model: BDF,
                         xb: np.ndarray,
                         dof_map: DOF_MAP,
                         nodes: np.ndarray,
                         xyz1: np.ndarray,
                         xyz2: np.ndarray,
                         A, I, J, E, G,
                         v, ihat, jhat, khat, wa, wb,
                         fdtype: str='float64'):
    """get the static CBAR force"""
    #words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
             #'0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
             #'       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
    nid1, nid2 = nodes
    # prop = elem.pid_ref
    # mat = prop.mid_ref
    I1, I2, I12 = I

    i1 = dof_map[(nid1, 1)]
    i2 = dof_map[(nid2, 1)]

    q_all = np.hstack([
        xb[i1:i1+6],
        xb[i2:i2+6],
    ])
    #q_axial = np.array([
        #xb[i1], xb[i1+1], xb[i1+2],
        #xb[i2], xb[i2+1], xb[i2+2]
    #])
    #q_torsion = np.array([
        #xb[i1+3], xb[i1+4], xb[i1+5],
        #xb[i2+3], xb[i2+4], xb[i2+5]
    #])

    #u_axial = Lambda @ q_axial
    #u_torsion = Lambda @ q_torsion

    is_passed, Ke = ke_cbar(
        xyz1, xyz2,
        A, I, J, E, G,
        ihat, jhat, khat, wa, wb,
        fdtype=fdtype)
    assert is_passed

    #pid_ref = elem.pid_ref
    #mat = pid_ref.mid_ref

    # ------------------
    #dxyz = xyz2 - xyz1
    #L = np.linalg.norm(dxyz)
    #pid_ref = elem.pid_ref
    #mat = pid_ref.mid_ref
    T = np.vstack([ihat, jhat, khat])
    z = np.zeros((3, 3), dtype=fdtype)
    Teb = np.block([
        [T, z, z, z],
        [z, T, z, z],
        [z, z, T, z],
        [z, z, z, T],
    ]) # 12x12
    q_element = Teb @ q_all
    u_e = q_element.reshape(12, 1)
    Fe = Ke @ q_element

    # ---------------------------------
    f_e = Fe
    #c, d, e, f = prop.get_cdef()
    #C1, C2 = c
    #D1, D2 = d
    #E1, E2 = e
    #F1, F2 = f
    # C1, C2, D1, D2, E1, E2, F1, F2 = cdef
    C1, C2, D1, D2, E1, E2, F1, F2 = [0.]*8
    # C1, C2, D1, D2, E1, E2, F1, F2 = prop.get_cdef().ravel()

    #f_e = obj.k_e * u_e
    force2stress = np.array([
        [1/A, 0, 0, 0, C2/I2, -C1/I1],
        [1/A, 0, 0, 0, D2/I2, -D1/I1],
        [1/A, 0, 0, 0, E2/I2, -E1/I1],
        [1/A, 0, 0, 0, F2/I2, -F1/I1],
    ])

    #print(force2stress.shape)
    #print(f_e.shape)
    stress = np.hstack([
        -force2stress @ f_e[:6],
        force2stress @ f_e[6:],
    ])
    #print(E, stress)

    # [End A Long. Stress or Strain at Point C;
    #  End A Long. Stress or Strain at Point D;
    #  End A Long. Stress or Strain at Point E;
    #  End A Long. Stress or Strain at Point F;
    #  End B Long. Stress or Strain at Point C;
    #  End B Long. Stress or Strain at Point D;
    #  End B Long. Stress or Strain at Point E;
    #  End B Long. Stress or Strain at Point F]
    strain = 1 / E * stress

    strain_energy = 0.5 * np.diag(u_e.T @ f_e).T
    #print('strain_energy =', strain_energy)
    # ---------------------------------
    #k1 = pid_ref.k1
    #k2 = pid_ref.k2
    #Ke = _beami_stiffness(pid_ref, mat, L, I1, I2, k1=k1, k2=k2, pa=pa, pb=pb)
    #K = Teb.T @ Ke @ Teb

    #is_passed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
    #T = np.vstack([ihat, jhat, khat])
    #z = np.zeros((3, 3), dtype='float64')
    #I1 = prop.I11()
    #I2 = prop.I22()
    #A = prop.Area()
    #J = prop.J()
    #unused_I12 = prop.I12()

    (Fx1, Fy1, Fz1, Mx1, My1, Mz1,
    Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fe

    axial = Fx1
    torque = Mx1
    shear1 = Fy1
    shear2 = Fz1
    bending_moment_a1 = My1
    bending_moment_a2 = Mz1

    bending_moment_b1 = My2
    bending_moment_b2 = Mz2

    out = (
        bending_moment_a1, bending_moment_a2,
        bending_moment_b1, bending_moment_b2,
        shear1, shear2,
        axial, torque)
    return out
