"""
defines:
 - export_caero_mesh(model, caero_bdf_filename='caero.bdf', is_subpanel_model=True)

"""
from __future__ import annotations
import math
from typing import TextIO, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.utils import PathLike
    from pyNastran.bdf.bdf import BDF, CAERO2, Coord, AELIST
from pyNastran.bdf.field_writer_8 import print_card_8


def export_caero_mesh(model: BDF,
                      caero_bdf_filename: PathLike='caero.bdf',
                      is_subpanel_model: bool=True,
                      pid_method: str='aesurf',
                      rotate_panel_angle_deg: float=0.0,
                      write_panel_xyz: bool=True) -> None:
    """
    Write the CAERO cards as CQUAD4s that can be visualized

    model: BDF
        a valid geometry
    caero_bdf_filename : str
        the file to write
    is_subpanel_model : bool; default=True
        True : write the subpanels as CQUAD4s
        False : write the macro elements as CQUAD4s
    pid_method : str; default='aesurf'
        'aesurf' : write the referenced AESURF as the property ID
                   main structure will be pid=1
        'caero' : write the CAERO1 as the property id
        'paero' : write the PAERO1 as the property id
    rotate_panel_angle_deg : float; default=0.0
        panel angle to rotate (e.g., rotate all surfaces by 30 degrees)
    write_panel_xyz : bool; default=True
        write the following table...
        $$  CAEROID      EID       XLE       YLE       ZLE     CHORD      SPAN   XLE+C/4
        $$        1        1    0.0000    0.2500    0.0000    0.0988    0.5000    0.0247
        $$        1        2    0.0988    0.2500    0.0000    0.0988    0.5000    0.1234
    """
    rotate_panel_angle = np.radians(rotate_panel_angle_deg)
    log = model.log
    if pid_method not in {'aesurf', 'caero', 'paero'}:
        raise RuntimeError(f'pid_method={pid_method!r} is not [aesurf, caero, paero]')

    inid = 1
    mid = 1
    log.info(f'export_caero_mesh -> {caero_bdf_filename}')

    #all_points = []
    aero_eid_map = {}
    #if is_subpanel_model:
    isubpanel_ieid = 0
    model.xref_obj.safe_cross_reference_aero()
    for caero_eid, caero in sorted(model.caeros.items()):
        if caero.type == 'CAERO2':
            log.warning('CAERO2 will probably cause issues...put it at the max id')
            continue
        points, elements = caero.panel_points_elements()
        for isubpanel_eid in range(len(elements)):
            aero_eid_map[isubpanel_ieid] = caero_eid + isubpanel_eid
            isubpanel_ieid += 1
    log.debug(f'  nsubpanels = {len(aero_eid_map)}')
    subcases, loads = _write_subcases_loads(model, aero_eid_map, is_subpanel_model)
    coords_to_write_dict = _get_coords_to_write_dict(model)

    eids_to_rotate_dict = {}
    aesurf_subpanel_eid_list = []
    for aesurf_id, aesurf in model.aesurf.items():
        cid1_ref: Coord = aesurf.cid1_ref
        aelist1_ref: AELIST = aesurf.aelist_id1_ref
        panel1_eids = aelist1_ref.elements
        eids_to_rotate_dict[(aesurf.label, 1)] = (cid1_ref, panel1_eids)
        aesurf_subpanel_eid_list.extend(panel1_eids)
        if aesurf.aelist_id2_ref is not None:
            cid2_ref: Coord = aesurf.cid2_ref
            aelist2_ref: AELIST = aesurf.aelist_id2_ref
            panel2_eids = aelist2_ref.elements
            eids_to_rotate_dict[(aesurf.label, 2)] = (cid2_ref, panel2_eids)
            aesurf_subpanel_eid_list.extend(panel2_eids)
    aesurf_subpanel_eids = np.array(aesurf_subpanel_eid_list, dtype='int32')

    with open(caero_bdf_filename, 'w') as bdf_file:
        #bdf_file.write('$ pyNastran: punch=True\n')
        bdf_file.write('SOL 101\n')
        bdf_file.write('CEND\n')
        bdf_file.write(subcases)
        bdf_file.write('BEGIN BULK\n')

        bdf_file.write(loads)
        _write_properties(model, bdf_file, pid_method=pid_method)

        for caero_eid, caero in sorted(model.caeros.items()):
            #assert caero_eid != 1, 'CAERO eid=1 is reserved for non-flaps'
            scaero = str(caero).rstrip().split('\n')
            if is_subpanel_model:
                if caero.type == 'CAERO2':
                    _write_caero2_subpanel(bdf_file, caero)
                    continue

                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
                if hasattr(caero, 'lspan'):
                    assert caero.type in {'CAERO1', 'CAERO4', 'CAERO5'}, caero
                    if caero.lspan_ref:
                        aefact_chord = str(caero.lspan_ref).rstrip().split('\n')
                        bdf_file.write('$ ' + '\n$ '.join(aefact_chord) + '\n')
                if hasattr(caero, 'lchord'):
                    assert caero.type in {'CAERO1', }, caero
                    if caero.lchord_ref:
                        aefact_span = str(caero.lchord_ref).rstrip().split('\n')
                        bdf_file.write('$ ' + '\n$ '.join(aefact_span) + '\n')

                #bdf_file.write("$   CAEROID       ID       XLE      YLE      ZLE     CHORD      SPAN\n")
                points, elements = caero.panel_points_elements()
                if write_panel_xyz:
                    _write_subpanel_strips(bdf_file, model, caero_eid, points, elements)

                box_ids = caero.box_ids.flatten()
                eids_aesurf = np.union1d(box_ids, aesurf_subpanel_eids)

                # verify eids_aesurf is sorted
                assert np.allclose(eids_aesurf, np.unique(eids_aesurf))
                if len(eids_aesurf) and rotate_panel_angle != 0.0:
                    nid_all = np.unique(elements.ravel())
                    # get the aesurf and fixed ids
                    # eids_aesurf = np.intersect1d(box_ids, aesurf_subpanel_eids)
                    eids_fixed = np.setdiff1d(box_ids, aesurf_subpanel_eids)

                    # get the index for each element
                    assert isinstance(box_ids, np.ndarray), box_ids
                    assert box_ids.ndim == 1, box_ids.shape
                    assert eids_fixed.ndim == 1, eids_fixed.shape
                    ieid_fixed = np.searchsorted(box_ids, eids_fixed)
                    # ieid_aesurf = np.searchsorted(box_ids, eids_aesurf)

                    # get the node ids associated with each element
                    nid_fixed = np.unique(elements[ieid_fixed, :])
                    # nid_aesurf = np.unique(elements[ieid_aesurf, :])

                    # get their index (so we can write the nodes in a rotated frame)
                    inid_fixed = np.searchsorted(nid_all, nid_fixed)

                    # get the base points and the points we're going to rotate
                    # note that they might overlap or be used by multiple flaps
                    #
                    # we'll write points_fixed later
                    points_fixed = points[inid_fixed, :]
                    elements_fixed = elements[ieid_fixed, :]

                    # now get do the same thing for each flap
                    for (label, idi), (cid_ref, eid_surface) in eids_to_rotate_dict.items():
                        eids_aesurf = np.intersect1d(box_ids, eid_surface)
                        if len(eids_aesurf) == 0:
                            continue

                        # print(f'box_ids     = {box_ids}; n={len(box_ids)}')
                        # print(f'eids_aesurf = {eids_aesurf}; n={len(eids_aesurf)}')
                        ieid_aesurf = np.searchsorted(box_ids, eids_aesurf)
                        # print(f'ieid_aesurf = {ieid_aesurf}')
                        # print('elements:\n', elements)
                        nid_aesurf = np.unique(elements[ieid_aesurf, :])

                        inid_aesurf = np.searchsorted(nid_all, nid_aesurf)
                        points_to_rotate = points[inid_aesurf, :]
                        elements_to_rotate = elements[ieid_aesurf, :]
                        xyz_rotated = rodriguez_rotate(
                            points_to_rotate, rotate_panel_angle,
                            cid_ref, iaxis=1)  # y-axis

                        # TODO: still need to renumber the points
                        #       and figure out what the offset is
                    raise NotImplementedError(f'rotate_panel_angle_deg={rotate_panel_angle_deg:g}; '
                                              f'rotate_panel_angle={rotate_panel_angle:g}')
                else:
                    npoints = points.shape[0]
                    #nelements = elements.shape[0]
                    for ipoint, point in enumerate(points):
                        x, y, z = point
                        bdf_file.write(print_card_8(['GRID', inid+ipoint, None, x, y, z]))

                #pid = caero_eid
                #mid = caero_eid
                jeid = 0
                for elem in elements + inid:
                    p1, p2, p3, p4 = elem
                    eid2 = jeid + caero_eid
                    pidi = _get_subpanel_property(
                        model, caero_eid, eid2, pid_method=pid_method)
                    fields = ['CQUAD4', eid2, pidi, p1, p2, p3, p4]
                    bdf_file.write(print_card_8(fields))
                    jeid += 1
            else:
                # macro model
                if caero.type == 'CAERO2':
                    continue
                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
                points = caero.get_points()
                npoints = 4
                for ipoint, point in enumerate(points):
                    x, y, z = point
                    bdf_file.write(print_card_8(['GRID', inid+ipoint, None, x, y, z]))

                pid = _get_subpanel_property(
                    model, caero_eid, caero_eid, pid_method=pid_method)
                p1 = inid
                p2 = inid + 1
                p3 = inid + 2
                p4 = inid + 3
                bdf_file.write(print_card_8(['CQUAD4', caero_eid, pid, p1, p2, p3, p4]))
            inid += npoints

        for cid, coord in sorted(coords_to_write_dict.items()):
            bdf_file.write(str(coord))

        # aluminum
        E = 350e9  # 350 GPa
        #G = None
        nu = 0.3
        rho = 2700.  # 2700 kg/m^3
        bdf_file.write(f'MAT1,{mid},{E},,{nu},{rho}\n')
        bdf_file.write('ENDDATA\n')
    log.debug(f'  ---finished export_caero_mesh of {caero_bdf_filename}---')
    return


def rodriguez_rotate(xyz: np.ndarray,
                     theta: float,
                     coord: Coord,
                     iaxis: int=0) -> np.ndarray:
    """
    v_rotated = v * cos(θ) + (k x v) * sin(θ) + k * (k · v) * (1 - cos(θ))
    Where:
      v is the vector to be rotated.
      k is the unit vector representing the axis of rotation.
      θ is the angle of rotation.
      x denotes the cross product.
      · denotes the dot product.

    :return:
    """
    if iaxis == 0:
        k = coord.i
    elif iaxis == 1:
        k = coord.j
    elif iaxis == 2:
        k = coord.k
    else:
        raise RuntimeError(f'iaxis={iaxis!r} and must be [0, 1, 2] for [x, y, z]')
    assert xyz.ndim == 2, f'xyz.shape={str(xyz.shape)} and must be 2d'

    nxyz = len(xyz)
    # thetas = np.full(nxyz, theta, dtype='float64')
    kmat = np.vstack([k] * nxyz, dtype='float64')
    # print(f'k; shape={str(k.shape)} = {k}')
    # kxv = np.cross(k[np.newaxis, :], xyz, axis=0)
    # kov = np.dot(k[np.newaxis, :], xyz, axis=0)
    try:
        kxv = np.cross(kmat, xyz, axis=1)
    except ValueError:  # pragma: no cover
        print(f'xyz; shape={str(xyz.shape)}:\n{xyz}')
        print(f'kmat; shape={str(kmat.shape)}:\n{kmat}')
        raise
    assert kxv.shape == (nxyz, 3), (kxv.shape, nxyz)
    # kov = np.dot(kmat, xyz, axis=0)
    kov = np.einsum('ij,ij->i', kmat, xyz)
    kov2 = kov.reshape(nxyz, 1)
    assert len(kov) == nxyz, (len(kov), nxyz)
    xyz_rotated = (
        xyz * np.cos(theta) +
        kxv * np.sin(theta) +
        kmat * kov2 * (1 - np.cos(theta))
    )
    return xyz_rotated

def _get_coords_to_write_dict(model: BDF) -> dict[int, Coord]:
    coords_to_write_dict = {}
    if model.aeros:
        rcsid_ref = model.aeros.rcsid_ref
        coords_to_write_dict[rcsid_ref.cid] = rcsid_ref
    for aesurf_id, aesurf in model.aesurf.items():
        #print(aesurf)
        if aesurf.cid1_ref is not None:
            #print(aesurf.cid1_ref)
            coords_to_write_dict[aesurf.cid1_ref.cid] = aesurf.cid1_ref
        if aesurf.cid2_ref is not None:
            coords_to_write_dict[aesurf.cid2_ref.cid] = aesurf.cid2_ref
    _add_traced_coords(model, coords_to_write_dict)
    return coords_to_write_dict


def _add_traced_coords(model: BDF,
                       coords_to_write_dict: dict[int, Coord]) -> None:
    nrids = 1
    while nrids > 0:
        cids = set(list(coords_to_write_dict))
        coords_to_add = set([])
        for cid in cids:
            coord = model.coords[cid]
            coords_to_add.add(coord.rid)
            coords_to_add.update(coord.rid_trace)  # seems to be busted
        coords_to_add_filtered = coords_to_add - cids
        coords_to_add_next = set([])
        for cid in coords_to_add_filtered:
            if cid == 0:
                continue
            coords_to_add_next.add(cid)
            coords_to_write_dict[cid] = model.coords[cid]
        nrids = len(coords_to_add_next)


def _write_caero2_subpanel(bdf_file: TextIO, caero: CAERO2):
    scaero = str(caero).rstrip().split('\n')
    bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
    if caero.lsb_ref:
        aefact_lsb = str(caero.lsb_ref).rstrip().split('\n')
        bdf_file.write('$ ' + '\n$ '.join(aefact_lsb) + '\n')
    if caero.lint_ref:
        aefact_lint = str(caero.lint_ref).rstrip().split('\n')
        bdf_file.write('$ ' + '\n$ '.join(aefact_lint) + '\n')
    #points, elements = caero.panel_points_elements()
    #raise NotImplementedError('caero2')
    return


def _write_subcases_loads(model: BDF,
                          aero_eid_map: dict[int, int],
                          is_subpanel_model: bool) -> tuple[str, str]:
    """writes the DMI, DMIJ, DMIK cards to a series of load cases"""
    nsubpanels = len(aero_eid_map)
    if len(model.dmi) == 0 and len(model.dmij) == 0 and len(model.dmik) == 0 and len(model.dmiji) == 0:
        loads = ''
        subcases = ''
        return subcases, loads

    log = model.log
    isubcase, loads, subcases = _write_dmi(model, aero_eid_map)

    for name, dmij in model.dmij.items():
        data, rows, cols = dmij.get_matrix(is_sparse=False, apply_symmetry=True)
        log.info(f'  {name}: shape={data.shape}')
        msg = f'{name}:\n'
        msg += str(data)

        if name in {'W2GJ', 'FA2J'}:
            assert len(rows) == len(data)
            assert data.shape[1] == 1, f'name={name}; shape={data.shape}'
            if name == 'W2GJ':
                data *= 180 / np.pi
                subtitle = f'DMI {name} (degrees)'
            else:
                subtitle = f'DMI {name}'
            subcases += (
                f'SUBCASE {isubcase}\n'
                f'  SUBTITLE = {subtitle}\n'
                f'  LOAD = {isubcase}\n')
            loads += f'$ {subtitle}\n'
            loads += '$ PLOAD2 SID P EID1\n'
            # aero_eid_map[isubpanel_ieid] = caero_eid + isubpanel_eid
            # raise NotImplementedError(msg)
            for irow, value in zip(rows, data):
                row = rows[irow]   # row = (1000,3)
                idi = row[0]
                eid = aero_eid_map[idi]
                loads += f'PLOAD2,{isubcase},{value[0]},{eid}\n'
            isubcase += 1
        else:  # pragma: no cover
            raise NotImplementedError(msg)

    for name, dmiji in model.dmiji.items():
        data, rows, cols = dmiji.get_matrix(is_sparse=False, apply_symmetry=True)
        log.info(f'  {name}: shape={data.shape}')
        msg = f'{name}:\n'
        msg += str(data)
        raise NotImplementedError(msg)

    for name, dmik in model.dmik.items():
        data, rows, cols = dmik.get_matrix(is_sparse=False, apply_symmetry=True)
        log.info(f'  {name}: shape={data.shape}')
        msg = f'{name}:\n'
        msg += str(data)
        if name == 'WKK':
            # column matrix of (neids*2,1)
            #assert data.shape[1] == 1, f'name={name}; shape={data.shape}'  # (112,1)
            if data.shape[1] != 1 and 0:
                # nsubpanels = 40
                # WKK = (80,80)
                log.warning(f'WKK is the wrong shape; shape={data.shape}')
                continue

            isubcase5 = isubcase + 1
            is_dof3 = False
            is_dof5 = False
            loads += f'$ {name} - FORCE/MOMENT\n'
            loads += '$ PLOAD2 SID P EID1\n'
            for irow, eid_dof1 in rows.items():
                eid1, dof1 = eid_dof1
                for icol, eid_dof2 in cols.items():
                    eid2, dof2 = eid_dof2
                    if eid1 == eid2 and dof1 == dof2:
                        value = data[irow, icol]
                        if dof1 == 3:
                            is_dof3 = True
                            loads += f'PLOAD2,{isubcase},{value},{eid1}\n'
                        else:
                            is_dof5 = True
                            assert dof1 == 5, dof1
                            loads += f'PLOAD2,{isubcase5},{value},{eid1}\n'
            #nrows = data.shape[0] // 2
            #data = data.reshape(nrows, 2)
            #subcases += (
                #f'SUBCASE {isubcase}\n'
                #f'  SUBTITLE = DMI {name} - FORCE\n'
                #f'  LOAD = {isubcase}\n')
            if is_dof3:
                subcases += (
                    f'SUBCASE {isubcase}\n'
                    f'  SUBTITLE = DMI {name} - FORCE\n'
                    f'  LOAD = {isubcase}\n'
                )
            if is_dof5:
                subcases += (
                    f'SUBCASE {isubcase5}\n'
                    f'  SUBTITLE = DMI {name} - MOMENT\n'
                    f'  LOAD = {isubcase5}\n'
                )
            isubcase += 2

            ## TODO: assume first column is forces & second column is moments...verify
            #loads += f'$ {name} - FORCE\n'
            #loads += '$ PLOAD2 SID P EID1\n'
            #for ieid, value in enumerate(data[:, 0].ravel()):
                #eid = aero_eid_map[ieid]
                #loads += f'PLOAD2,{isubcase},{value},{eid}\n'
            #isubcase += 1

            #subcases += (
                #f'SUBCASE {isubcase}\n'
                #f'  SUBTITLE = DMI {name} - MOMENT\n'
                #f'  LOAD = {isubcase}\n')
            #loads += f'$ {name} - MOMENT\n'
            #loads += '$ PLOAD2 SID P EID1\n'
            #for irow, value in enumerate(data[:, 1].ravel()):
                #row = rows[irow]
                #eid = row[0]
                #loads += f'PLOAD2,{isubcase},{value},{eid}\n'
        else:
            raise NotImplementedError(msg)

    if not is_subpanel_model:
        # we put this here to test
        #model.log.warning('cannot export "loads" because not a subpanel model')
        subcases = ''
        loads = ''
    return subcases, loads


def _write_dmi(model: BDF,
               aero_eid_map: dict[int, int]) -> tuple[int, str, str]:
    """writes the DMI cards to a series of load cases"""
    nsubpanels = len(aero_eid_map)
    isubcase = 1
    loads = ''
    subcases = ''
    log = model.log
    for name, dmi in model.dmi.items():
        form_str = dmi.matrix_form_str
        data, rows, cols = dmi.get_matrix(is_sparse=False, apply_symmetry=True)
        log.debug(f'  {name}: shape={data.shape}; form_str={form_str}')

        if name == 'WTFACT':
            # square matrix of (neids,neids) that has values on only? the diagonal
            subcases += (
                f'SUBCASE {isubcase}\n'
                f'  SUBTITLE = DMI {name} - FORCE\n'
                f'  LOAD = {isubcase}\n')

            # diagonal matrix
            diag = data.diagonal()
            data = data.copy()
            np.fill_diagonal(data, 0.)
            assert np.abs(data).max() == 0.0, 'WTFACT has cross terms'

            loads += f'$ {name} - FORCE\n'
            loads += '$ PLOAD2 SID P EID1\n'
            for ieid, value in enumerate(diag[::2]):
                eid = aero_eid_map[ieid]
                loads += f'PLOAD2,{isubcase},{value},{eid}\n'
            isubcase += 1

            subcases += (
                f'SUBCASE {isubcase}\n'
                f'  SUBTITLE = DMI {name} - MOMENT\n'
                f'  LOAD = {isubcase}\n')
            loads += f'$ {name} - MOMENT\n'
            loads += '$ PLOAD2 SID P EID1\n'
            for ieid, value in enumerate(diag[1::2]):
                eid = aero_eid_map[ieid]
                loads += f'PLOAD2,{isubcase},{value},{eid}\n'

        elif name in {'W2GJ', 'FA2J'}:
            # column matrix of (neids,1)
            # boxid and not dof in the k-set
            assert data.shape[1] == 1, f'name={name}; shape={data.shape}'  # (56,1)
            subcases += (
                f'SUBCASE {isubcase}\n'
                f'  SUBTITLE = DMI {name}\n'
                f'  LOAD = {isubcase}\n')

            # (56,1)
            if name == 'W2GJ':
                data *= 180 / np.pi
                loads += f'$ {name} (degrees)\n'
            else:
                loads += f'$ {name}\n'
            loads += '$ PLOAD2 SID P EID1\n'
            for ieid, value in enumerate(data.ravel()):
                eid = aero_eid_map[ieid]
                loads += f'PLOAD2,{isubcase},{value},{eid}\n'
        elif name == 'WKK':
            # column matrix of (neids*2,1)

            # 1: 'square',
            # 2: 'rectangular',  # 9 ???
            # 3: 'diagonal',
            # 6: 'symmetric',
            # 9: 'identity',
            matrix_form_str = dmi.matrix_form_str
            if matrix_form_str == 'diagonal':
                pass
            elif matrix_form_str in {'square', 'rectangular'}:
                assert data.shape[0] == data.shape[1], f'name={name}; shape={data.shape} matrix_form={dmi.matrix_form}={matrix_form_str}'  # (100,100)
                assert data.shape == (2*nsubpanels, 2*nsubpanels), f'name={name}; shape={data.shape} expected=({2*nsubpanels},{2*nsubpanels}); matrix_form={dmi.matrix_form}={matrix_form_str}'
                matrix_form_str = 'square'
            else:
                matrix_form_str = 'column'
                assert data.shape[1] == 1, f'name={name}; shape={data.shape} matrix_form={dmi.matrix_form}={matrix_form_str}'  # (112,1)

            if matrix_form_str in {'diagonal', 'column'}:
                assert data.shape[1] == 1, f'name={name}; shape={data.shape}'  # (112,1)
                assert data.shape == (2*nsubpanels, 1), f'name={name}; shape={data.shape}  expected=({2*nsubpanels},{2*nsubpanels}); matrix_form={dmi.matrix_form}={matrix_form_str}'
                nrows = data.shape[0] // 2
                data = data.reshape(nrows, 2)
                force_correction = data[:, 0]
                moment_correction = data[:, 1]
                corrections = [
                    ('FORCE', force_correction),
                    ('MOMENT', moment_correction),
                ]
            elif matrix_form_str == 'square':
                assert data.shape[0] == data.shape[1]
                force = data[::2, ::2]
                moment = data[1::2, 1::2]
                force_correction_row = force.sum(axis=0)
                force_correction_col = force.sum(axis=1)
                moment_correction_row = moment.sum(axis=0)
                moment_correction_col = moment.sum(axis=1)

                corrections = [
                    ('FORCE_ROW', force_correction_row),
                    ('FORCE_COL', force_correction_col),
                    ('MOMENT_ROW', moment_correction_row),
                    ('MOMENT_COL', moment_correction_col),
                ]
            else:  # pragma: no cover
                raise NotImplementedError(matrix_form_str)

            for result_name, correction_column in corrections:
                #print(result_name, isubcase)
                subcases += (
                    f'SUBCASE {isubcase}\n'
                    f'  SUBTITLE = DMI {name} - {result_name}\n'
                    f'  LOAD = {isubcase}\n')

                loads += f'$ {name} - {result_name}\n'
                loads += '$ PLOAD2 SID P EID1\n'
                for ieid, value in enumerate(correction_column.ravel()):
                    eid = aero_eid_map[ieid]
                    loads += f'PLOAD2,{isubcase},{value},{eid}\n'
                isubcase += 1
        else:  # pragma: no cover
            msg = f'{name} matrix_form_str={matrix_form_str}:\n'
            msg += str(data)
            #continue
            raise NotImplementedError(msg)
        isubcase += 1
    return isubcase, loads, subcases


def _write_subpanel_strips(bdf_file: TextIO, model: BDF,
                           caero_eid: int,
                           points: np.ndarray, elements: np.ndarray) -> None:
    """writes the strips for the subpanels"""
    #bdf_file.write("$   CAEROID       ID       XLE      YLE      ZLE     CHORD      SPAN\n")
    bdf_file.write('$$\n$$ XYZ_LE is taken at the center of the leading edge; (p1+p4)/2\n$$\n')
    bdf_file.write('$$ %8s %8s %9s %9s %9s %9s %9s %9s %9s\n' % (
        'CAEROID', 'EID', 'XLE', 'YLE', 'ZLE', 'CHORD', 'SPAN', 'XLE+C/4', 'XLE+C/2'))

    for i in range(elements.shape[0]):
        # The point numbers here are consistent with the CAERO1
        p1 = points[elements[i, 0], :]
        p4 = points[elements[i, 1], :]
        p2 = points[elements[i, 2], :]
        p3 = points[elements[i, 3], :]
        le: list[float] = (p1 + p4)*0.5
        te: list[float] = (p2 + p3)*0.5
        dy = (p4 - p1)[1]
        dz = (p4 - p1)[2]
        span = math.sqrt(dy**2 + dz**2)
        chord: float = te[0] - le[0]
        xqc: float = le[0] + chord / 4.
        xmid: float = le[0] + chord / 2.
        bdf_file.write("$$ %8d %8d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n" % (
            caero_eid, caero_eid+i, le[0], le[1], le[2], chord, span, xqc, xmid))


def _get_subpanel_property(model: BDF, caero_id: int, eid: int,
                           pid_method: str='aesurf') -> int:
    """gets the property id for the subpanel"""
    pid = None
    if pid_method == 'aesurf':
        for aesurf_id, aesurf in model.aesurf.items():
            aelist_id = aesurf.Aelist_id1()
            aelist = model.aelists[aelist_id]
            if eid in aelist.elements:
                pid = aesurf_id
                break
    elif pid_method == 'caero':
        pid = caero_id
    elif pid_method == 'paero':
        caero = model.caeros[caero_id]
        pid = caero.pid
    else:  # pragma: no cover
        raise RuntimeError(f'pid_method={pid_method!r} is not [aesurf, caero, paero]')

    if pid is None:
        pid = 1
    return pid


def _write_properties(model: BDF, bdf_file: TextIO,
                      pid_method: str='aesurf') -> None:
    """writes the PSHELL with a different comment depending on the pid_method flag"""
    if pid_method == 'aesurf':
        _write_aesurf_properties(model, bdf_file)
    elif pid_method == 'paero':
        for paero_id, paero in sorted(model.paeros.items()):
            spaero = str(paero).rstrip().split('\n')
            bdf_file.write('$ ' + '\n$ '.join(spaero) + '\n')
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (paero_id, 1))
    elif pid_method == 'caero':
        for caero_eid, caero in sorted(model.caeros.items()):
            scaero = str(caero).rstrip().split('\n')
            bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (caero_eid, 1))
    else:  # pragma: no cover
        raise RuntimeError(f'pid_method={repr(pid_method)} is not [aesurf, caero, paero]')


def _write_aesurf_properties(model: BDF, bdf_file: TextIO) -> None:
    """the AESURF ID will be the PSHELL ID"""
    aesurf_mid = 1
    for aesurf_id, aesurf in model.aesurf.items():
        #cid = aesurf.cid1

        #aesurf_mid = aesurf_id
        saesurf = str(aesurf).rstrip().split('\n')
        bdf_file.write('$ ' + '\n$ '.join(saesurf) + '\n')
        bdf_file.write('PSHELL,%s,%s,0.1\n' % (aesurf_id, aesurf_mid))
        #print(cid)
        #ax, ay, az = cid.i
        #bx, by, bz = cid.j
        #cx, cy, cz = cid.k
        #bdf_file.write('CORD2R,%s,,%s,%s,%s,%s,%s,%s\n' % (
            #cid, ax, ay, az, bx, by, bz))
        #bdf_file.write(',%s,%s,%s\n' % (cx, cy, cz))
        #print(cid)
        #aesurf.elements
    # dummy property
    bdf_file.write('PSHELL,%s,%s,0.1\n' % (1, 1))
    return
