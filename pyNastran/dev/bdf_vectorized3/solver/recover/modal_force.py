from typing import TextIO
import numpy as np
from pyNastran.bdf.bdf import BDF, Subcase
from .utils import (
    get_f06_op2_pch_set, #get_mag_phase_from_options,
    get_ieids_eids)
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealSpringForceArray


def recover_force_103(f06_file: TextIO, op2,
                      model: BDF, dof_map,
                      subcase: Subcase,
                      phig: np.ndarray,
                      phib: np.ndarray,
                      modes, eigns, #freqs,
                      cfdtype: str='complex64',
                      title: str='', subtitle: str='', label: str='',
                      page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - FORCE = ALL

    """
    omegas = np.sqrt(np.abs(eigns))
    freqs = omegas / (2 * np.pi)

    log = model.log
    if 'FORCE' not in subcase:
        return
    write_f06, write_op2, write_pch, options, eids_write = get_f06_op2_pch_set(
        subcase, 'FORCE')
    write_data = np.any((write_f06, write_op2))
    if not write_data:
        log.warning(f'skipping modal force')
        return page_num

    isubcase = subcase.id

    ngrid = len(model.grid)
    nspoint = len(model.spoint)
    assert nspoint == 0, nspoint
    nmode = phig.shape[0]
    xg = phig.reshape(nmode, ngrid, 6)
    nelements = 0
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xg, eids_write, options,
        modes, eigns, freqs,
        'CELAS1', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xg, eids_write, options,
        modes, eigns, freqs,
        'CELAS2', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xg, eids_write, options,
        modes, eigns, freqs,
        'CELAS3', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xg, eids_write, options,
        modes, eigns, freqs,
        'CELAS4', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)

    # nelements += _recover_force_rod(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CROD', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    # nelements += _recover_force_rod(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CONROD', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    # nelements += _recover_force_rod(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CTUBE', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    # nelements += _recover_force_cbar(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CBAR', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        log.warning(f'no force output...{model.card_count}; {model.bdf_filename}')
        # asdf


def _recover_force_celas(f06_file: TextIO, op2,
                         model: BDF,
                         isubcase: int,
                         xg: np.ndarray,
                         eids_write: np.ndarray,  options: list[str],
                         modes, eigns, freqs,
                         element_name: str, cfdtype='complex64',
                         title: str='', subtitle: str='', label: str='',
                         page_num: int=1, page_stamp='PAGE %s',
                         nelements: int=0,
                         write_f06: bool=False, write_op2: bool=False,
                         ) -> tuple[int, int]:
    """recovers static spring force"""
    neid, elem, ieid, element_id = get_ieids_eids(model, element_name, eids_write)
    nelements += neid
    if not neid:
        return nelements, page_num

    nmode = xg.shape[0]

    forces = op2.op2_results.force
    slot = getattr(forces, element_name.lower() + '_force')

    if element_name in {'CELAS1', 'CELAS3'}:
        pid = elem.property_id[ieid]
        pelas = model.pelas.slice_card_by_property_id(pid)
        k = pelas.k.ravel()
    else:
        k = elem.k[ieid].ravel()
    assert len(k) == len(elem)

    if element_name in {'CELAS1', 'CELAS2'}:
        force = _celas2_force(model, elem, k, xg)
    # elif element_name in {'CELAS3', 'CELAS4'}:
    #     force = _celas4_force(model, elem, k, xg)
    #     force = _recover_forcei_celas34(xg, dof_map, elem, k)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    table_name = 'OEF1'
    data = force.reshape(nmode, neid, 1)
    assert data.shape == (nmode, neid, 1), data.shape
    obj = RealSpringForceArray.add_modal_case(
        table_name, element_name, element_id, data, isubcase,
        modes, eigns, freqs,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    header = ['', '', '']
    if write_f06:
        page_num = obj.write_f06(
            f06_file, header=header,
            page_stamp=page_stamp, page_num=page_num,
            is_mag_phase=False, is_sort1=True)
        f06_file.write('\n')
    else:
        raise RuntimeError('write the f06')

    slot[isubcase] = obj
    # if write_op2:
    #     asdf
    # write_f06_force = True
    # _save_spring_force(
    #     op2, f06_file, page_num, page_stamp,
    #     element_name,
    #     force, eids, write_f06_force,
    #     isubcase, title, subtitle, label)
    return nelements, page_num


def _celas2_force(model: BDF,
                  elem, k: np.ndarray,
                  xg: np.ndarray) -> np.ndarray:
    if elem.nodes.min() == 0 or elem.components.min() == 0:
        inid = np.full(elem.nodes.shape)
        i0 = np.where(elem.nodes > 0)
        raise RuntimeError(f'SPOINT {elem.type}')
    else:
        inid = np.searchsorted(model.grid.node_id, elem.nodes)
        inid1 = inid[:, 0]
        inid2 = inid[:, 1]
        comp1 = elem.components[:, 0] - 1
        comp2 = elem.components[:, 0] - 1
        # xg1 = xg[:, inid1, :][:, comp1]
        # xg2 = xg[:, inid2, :][:, comp2]
        xg1 = xg[:, inid1, comp1]
        xg2 = xg[:, inid2, comp2]
        dx = xg2 - xg1

    force = k * dx  # (nmode, neid)
    return force

# def _recover_force_rod(f06_file, op2,
#                        model: BDF, dof_map, isubcase, phib, eids_str,
#                        element_name, fdtype='float32',
#                        title: str='', subtitle: str='', label: str='',
#                        page_num: int=1, page_stamp='PAGE %s') -> None:
#     """recovers static rod force"""
#     neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
#     if not neids:
#         return neids
#     forces = np.full((neids, 2), np.nan, dtype=fdtype)
#     if element_name == 'CONROD':
#         for ieid, eid in zip(irod, eids):
#             elem = model.elements[eid]
#             forces[ieid, :] = _recover_forcei_rod(xb, dof_map, elem, elem)
#     elif element_name == 'CROD':
#         for ieid, eid in zip(irod, eids):
#             elem = model.elements[eid]
#             forces[ieid, :] = _recover_forcei_rod(xb, dof_map, elem, elem.pid_ref)
#     elif element_name == 'CTUBE':
#         for ieid, eid in zip(irod, eids):
#             elem = model.elements[eid]
#             forces[ieid, :] = _recover_forcei_rod(xb, dof_map, elem, elem.pid_ref)
#     else:  # pragma: no cover
#         raise NotImplementedError(element_name)
#
#     data = forces.reshape(1, *forces.shape)
#     table_name = 'OEF1'
#     force_obj = RealRodForceArray.add_static_case(
#         table_name, element_name, eids, data, isubcase,
#         is_sort1=True, is_random=False, is_msc=True,
#         random_code=0, title=title, subtitle=subtitle, label=label)
#
#     force = op2.op2_results.force
#     if element_name == 'CONROD':
#         force.conrod_force[isubcase] = force_obj
#     elif element_name == 'CROD':
#         force.crod_force[isubcase] = force_obj
#     elif element_name == 'CTUBE':
#         force.ctube_force[isubcase] = force_obj
#     else:  # pragma: no cover
#         raise NotImplementedError(element_name)
#
#     force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
#                         page_num=page_num, is_mag_phase=False, is_sort1=True)
#     return neids
#
# def _recover_forcei_rod(xb, dof_map, elem, prop):
#     """get the static rod force"""
#     nid1, nid2 = elem.nodes
#
#     i1 = dof_map[(nid1, 1)]
#     i2 = dof_map[(nid2, 1)]
#
#     q_axial = np.array([
#         xb[i1], xb[i1+1], xb[i1+2],
#         xb[i2], xb[i2+1], xb[i2+2]
#     ])
#     q_torsion = np.array([
#         xb[i1+3], xb[i1+4], xb[i1+5],
#         xb[i2+3], xb[i2+4], xb[i2+5]
#     ])
#     xyz1 = elem.nodes_ref[0].get_position()
#     xyz2 = elem.nodes_ref[1].get_position()
#     dxyz12 = xyz1 - xyz2
#     Lambda = lambda1d(dxyz12, debug=False)
#
#     u_axial = Lambda @ q_axial
#     u_torsion = Lambda @ q_torsion
#     du_axial = u_axial[0] - u_axial[1]
#     du_torsion = u_torsion[0] - u_torsion[1]
#     #headers = ['axial', 'SMa', 'torsion', 'SMt']
#
#     #C = prop.c
#     mat = prop.mid_ref
#
#     L = np.linalg.norm(dxyz12)
#     G = mat.G()
#     J = elem.J()
#     A = elem.Area()
#     E = elem.E()
#
#     axial_strain = du_axial / L
#     #torsional_strain = du_torsion * C / L
#
#     axial_stress = E * axial_strain
#     #torsional_stress = G * torsional_strain
#     axial_force = axial_stress * A
#     torsional_moment = du_torsion * G * J / L
#
#     return axial_force, torsional_moment
#
# def _recover_force_cbar(f06_file, op2,
#                         model: BDF, dof_map, isubcase, xb, eids_str,
#                         element_name, fdtype='float32',
#                         title: str='', subtitle: str='', label: str='',
#                         page_num: int=1, page_stamp='PAGE %s') -> None:
#     """
#     Recovers static CBAR force.
#
#     .. todo:: doesn't support CBAR-100
#
#     """
#     neids, irod, eids = get_ieids_eids(model, element_name, eids_str)
#     if not neids:
#         return neids
#     forces = np.full((neids, 8), np.nan, dtype=fdtype)
#
#     for ieid, eid in zip(irod, eids):
#         elem = model.elements[eid]
#         forces[ieid, :] = _recover_forcei_cbar(model, xb, dof_map, elem, elem.pid_ref)
#
#     data = forces.reshape(1, *forces.shape)
#     table_name = 'OEF1'
#     force_obj = RealCBarForceArray.add_static_case(
#         table_name, 'CBAR', eids, data, isubcase,
#         is_sort1=True, is_random=False, is_msc=True,
#         random_code=0, title=title, subtitle=subtitle, label=label)
#
#     force = op2.op2_results.force
#     force.cbar_force[isubcase] = force_obj
#
#     force_obj.write_f06(f06_file, header=None, page_stamp=page_stamp,
#                         page_num=page_num, is_mag_phase=False, is_sort1=True)
#     return neids
#
# def _recover_forcei_cbar(model: BDF,
#                          xb, dof_map, elem: CBAR,
#                          prop: PBAR | PBARL, fdtype: str='float64'):
#     """get the static CBAR force"""
#     #words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
#              #'0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
#              #'       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
#     nid1, nid2 = elem.nodes
#     prop = elem.pid_ref
#     mat = prop.mid_ref
#
#     i1 = dof_map[(nid1, 1)]
#     i2 = dof_map[(nid2, 1)]
#
#     q_all = np.hstack([
#         xb[i1:i1+6],
#         xb[i2:i2+6],
#     ])
#     #q_axial = np.array([
#         #xb[i1], xb[i1+1], xb[i1+2],
#         #xb[i2], xb[i2+1], xb[i2+2]
#     #])
#     #q_torsion = np.array([
#         #xb[i1+3], xb[i1+4], xb[i1+5],
#         #xb[i2+3], xb[i2+4], xb[i2+5]
#     #])
#
#     #u_axial = Lambda @ q_axial
#     #u_torsion = Lambda @ q_torsion
#
#     nid1, nid2 = elem.nodes
#     is_passed, Ke = ke_cbar(model, elem, fdtype=fdtype)
#     assert is_passed
#
#     #pid_ref = elem.pid_ref
#     #mat = pid_ref.mid_ref
#
#
#     # ------------------
#     is_failed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
#     assert is_failed is False
#     #print(wa, wb)
#     #xyz1 = elem.nodes_ref[0].get_position() + wa
#     #xyz2 = elem.nodes_ref[1].get_position() + wb
#     #dxyz = xyz2 - xyz1
#     #L = np.linalg.norm(dxyz)
#     #pid_ref = elem.pid_ref
#     #mat = pid_ref.mid_ref
#     T = np.vstack([ihat, jhat, khat])
#     z = np.zeros((3, 3), dtype=fdtype)
#     Teb = np.block([
#         [T, z, z, z],
#         [z, T, z, z],
#         [z, z, T, z],
#         [z, z, z, T],
#     ]) # 12x12
#     q_element = Teb @ q_all
#     u_e = q_element.reshape(12, 1)
#     Fe = Ke @ q_element
#
#     # ---------------------------------
#     f_e = Fe
#     #c, d, e, f = prop.get_cdef()
#     #C1, C2 = c
#     #D1, D2 = d
#     #E1, E2 = e
#     #F1, F2 = f
#     C1, C2, D1, D2, E1, E2, F1, F2 = prop.get_cdef().ravel()
#     A = prop.A
#     E = mat.E()
#     I1 = prop.I11()
#     I2 = prop.I22()
#
#     #f_e = obj.k_e * u_e
#     force2stress = np.array([
#         [1/A, 0, 0, 0, C2/I2, -C1/I1],
#         [1/A, 0, 0, 0, D2/I2, -D1/I1],
#         [1/A, 0, 0, 0, E2/I2, -E1/I1],
#         [1/A, 0, 0, 0, F2/I2, -F1/I1],
#     ])
#
#     #print(force2stress.shape)
#     #print(f_e.shape)
#     stress = np.hstack([
#         -force2stress @ f_e[:6],
#         force2stress @ f_e[6:],
#     ])
#     #print(E, stress)
#
#     # [End A Long. Stress or Strain at Point C;
#     #  End A Long. Stress or Strain at Point D;
#     #  End A Long. Stress or Strain at Point E;
#     #  End A Long. Stress or Strain at Point F;
#     #  End B Long. Stress or Strain at Point C;
#     #  End B Long. Stress or Strain at Point D;
#     #  End B Long. Stress or Strain at Point E;
#     #  End B Long. Stress or Strain at Point F]
#     strain = 1 / E * stress
#
#     strain_energy = 0.5 * np.diag(u_e.T @ f_e).T
#     #print('strain_energy =', strain_energy)
#     # ---------------------------------
#     #k1 = pid_ref.k1
#     #k2 = pid_ref.k2
#     #Ke = _beami_stiffness(pid_ref, mat, L, I1, I2, k1=k1, k2=k2, pa=pa, pb=pb)
#     #K = Teb.T @ Ke @ Teb
#
#     #is_passed, (v, ihat, jhat, khat, wa, wb) = elem.get_axes(model)
#     #T = np.vstack([ihat, jhat, khat])
#     #z = np.zeros((3, 3), dtype='float64')
#     #I1 = prop.I11()
#     #I2 = prop.I22()
#     #A = prop.Area()
#     #J = prop.J()
#     #unused_I12 = prop.I12()
#
#     (Fx1, Fy1, Fz1, Mx1, My1, Mz1,
#     Fx2, Fy2, Fz2, Mx2, My2, Mz2) = Fe
#
#     axial = Fx1
#     torque = Mx1
#     shear1 = Fy1
#     shear2 = Fz1
#     bending_moment_a1 = My1
#     bending_moment_a2 = Mz1
#
#     bending_moment_b1 = My2
#     bending_moment_b2 = Mz2
#
#     out = (
#         bending_moment_a1, bending_moment_a2,
#         bending_moment_b1, bending_moment_b2,
#         shear1, shear2,
#         axial, torque)
#     return out
