from collections import defaultdict

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF, CaseControlDeck
from pyNastran.converters.abaqus.abaqus_cards import Part, Material, Step, ShellSection, SolidSection

from pyNastran.converters.abaqus.abaqus import (
    Abaqus, read_abaqus, get_nodes_nnodes_nelements)

def _add_part_to_nastran(nastran_model: BDF, elements, pid: int, nid_offset: int) -> None:
    for etype, eids_nids in elements.element_types.items():
        eids, part_nids = eids_nids
        if eids is None and part_nids is None:
            continue

        if nid_offset > 0:
            # don't use += or it's an inplace operation
            part_nids = part_nids + nid_offset

        if etype == 'b31h':
            for eid, nids in zip(eids, part_nids):
                x = None
                g0 = nids[2]
                nids = [nids[0], nids[1]]
                nastran_model.add_cbeam(
                    eid, pid, nids, x, g0, offt='GGG', bit=None,
                    pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')
                pass
        elif etype == 'cpe4':
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, T4=None, comment='')
        elif etype == 's8r':
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cquad8(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, T4=None, comment='')
        elif etype == 'c3d4':
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_ctetra(eid, pid, nids, comment='')
        else:
            raise NotImplementedError(etype)

    #add_lines(grid, nidsi, part.r2d2, nid_offset)

    #add_tris(grid, nidsi, part.cps3, nid_offset)
    #add_tris(grid, nidsi, part.cpe3, nid_offset)

    #add_quads(grid, nidsi, part.cpe4, nid_offset)
    #add_quads(grid, nidsi, part.cpe4r, nid_offset)

    #add_quads(grid, nidsi, part.cps4, nid_offset)
    #add_quads(grid, nidsi, part.cps4r, nid_offset)

    #add_quads(grid, nidsi, part.coh2d4, nid_offset)
    #add_quads(grid, nidsi, part.cohax4, nid_offset)

    #add_tris(grid, nidsi, part.cax3, nid_offset)
    #add_quads(grid, nidsi, part.cax4, nid_offset)
    #add_quads(grid, nidsi, part.cax4r, nid_offset)

    ## solids
    #add_tetras(grid, nidsi, part.c3d10h, nid_offset)
    #add_hexas(grid, nidsi, part.c3d8r, nid_offset)


def _create_nastran_nodes_elements(model: Abaqus, nastran_model: BDF) -> None:
    log = model.log
    nid_offset = 0

    pid = 1
    if model.nids is not None and len(model.nids):
        nnodesi = model.nodes.shape[0]
        elements = model.elements
        _add_part_to_nastran(nastran_model, elements, pid, nid_offset)
        nid_offset += nnodesi

    for unused_part_name, part in model.parts.items():
        #log.info('part_name = %r' % unused_part_name)
        nnodesi = part.nodes.shape[0]
        #nidsi = part.nids

        elements = part.elements
        _add_part_to_nastran(nastran_model, elements, pid, nid_offset)
        #nids.append(nidsi)

        nid_offset += nnodesi
        for shell_section in part.shell_sections:
            log.info('shell')
        for shell_section in part.shell_sections:
            log.info('solid')
    #nids = np.hstack(nids)

    pid = 1
    mid = 1
    for solid_section in model.solid_sections:
        #print(solid_section)
        mat_name = solid_section.material_name
        mat = model.materials[mat_name]
        element_set = solid_section.elset
        log.warning(f'element_set = {element_set}')
        #print(mat)
        nastran_model.add_psolid(pid, mid, cordm=0,
                                 integ=None, stress=None, isop=None, fctn='SMECH', comment='')
        G = None
        elastic = mat.sections['elastic']
        E = elastic[0]
        nu = elastic[1]
        nastran_model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0,
                               ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        log.info('shell2')
        mid += 1
    for shell_section in model.shell_sections:
        #print(shell_section)
        mat_name = shell_section.material_name
        mat = model.materials[mat_name]
        #print(mat)
        t = shell_section.thickness
        nastran_model.add_pshell(pid, mid1=mid, t=t, mid2=mid, twelveIt3=1.0, mid3=None, tst=0.833333,
                                 nsm=0.0, z1=None, z2=None, mid4=None, comment='')
        G = None
        elastic = mat.sections['elastic']
        E = elastic[0]
        nu = elastic[1]
        nastran_model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0,
                               ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        log.info('shell2')
        mid += 1

def abaqus_to_nastran_filename(abaqus_inp_filename: str,
                               nastran_filename_out: str,
                               log=None) -> BDF:
    if isinstance(abaqus_inp_filename, Abaqus):
        model = abaqus_inp_filename
    else:
        model = read_abaqus(abaqus_inp_filename, log=log, debug=True)
    nnodes, nids, nodes, nelements = get_nodes_nnodes_nelements(
        model, stop_for_no_elements=True)
    assert nnodes > 0, nnodes
    assert nelements > 0, nelements

    nastran_model = BDF(debug=True, log=log, mode='msc')
    for nid, xyz in zip(nids, nodes):
        nastran_model.add_grid(nid, xyz)
    _create_nastran_nodes_elements(model, nastran_model)
    #for step in model.steps:
        #print(step)

    _create_nastran_loads(model, nastran_model)
    try:
        str(nastran_model.case_control_deck)
    except AssertionError:
        log.error('No case control deck found...skipping')
        #print(f'{nastran_model.case_control_deck}')
        nastran_model.case_control_deck = None
    nastran_model.write_bdf(nastran_filename_out)
    x = 1
    return nastran_model

def _write_boundary_as_nastran(model: Abaqus,
                               nastran_model: BDF,
                               case_control_deck: CaseControlDeck) -> None:
    if not model.boundaries:
        return
    for iboundary, boundary in enumerate(model.boundaries):
        fixed_spcs = defaultdict(list)
        for (nid, dof), value in boundary.nid_dof_to_value.items():
            assert isinstance(nid, int), nid
            if value == 0.0:
                fixed_spcs[dof].append(nid)
            else:
                raise RuntimeError((nid, dof, value))

        subcase_id = iboundary + 1
        spc_id = iboundary + 1
        subcase = case_control_deck.create_new_subcase(subcase_id)
        subcase.add('SPC', spc_id, [], 'STRESS-type')
        for dof, nids in fixed_spcs.items():
            #  spc_id: int, components: str, nodes
            nastran_model.add_spc1(spc_id, str(dof), nids)

def _create_nastran_loads(model: Abaqus, nastran_model: BDF):
    def xyz():
        return np.zeros(3, dtype='float32')

    if nastran_model.case_control_deck is None:
        nastran_model.case_control_deck = CaseControlDeck([], log=nastran_model.log)

    case_control_deck = nastran_model.case_control_deck
    _write_boundary_as_nastran(model, nastran_model, case_control_deck)

    output_map = {
        'U': 'DISP',
        'RF': 'SPCFORCE',
        'S': 'STRESS',
        'E': 'STRAIN',
        'ENER': 'ESE',
    }
    for istep, step in enumerate(model.steps):
        subcase_id = istep + 1
        load_id = subcase_id
        if subcase_id in nastran_model.subcases:
            subcase = nastran_model.subcases[subcase_id]
        else:
            subcase = case_control_deck.create_new_subcase(subcase_id)

        for output in step.node_output + step.element_output:

            try:
                base_output = output_map[output.upper()]
            except KeyError:
                raise KeyError(f'output={output!r} is not in [u, rf, s, e, ener]')
            subcase.add(base_output, 'ALL', ['PRINT', 'PLOT'], 'STRESS-type')

        for cload in step.cloads:
            assert nastran_model.sol is None or nastran_model.sol == 101, nastran_model.sol
            nastran_model.sol = 101
            subcase.add('LOAD', load_id, [], 'STRESS-type')
            #subcase['LOAD'] = load_id
            forces = defaultdict(xyz)
            moments = defaultdict(xyz)
            for cloadi in cload:
                nid, dof, mag = cloadi
                assert dof in [1, 2, 3], cload
                if dof in {1, 2, 3}:
                    forces[nid][dof - 1] = mag
                elif dof in {4, 5, 6}:
                    moments[nid][dof - 4] = mag
                else:
                    raise NotImplementedError(cloadi)

            mag = 1.0
            if len(forces):
                for nid, xyz in forces.items():
                    nastran_model.add_force(load_id, nid, mag, xyz, cid=0, comment='')
            if len(moments):
                for nid, xyz in moments.items():
                    nastran_model.add_moment(load_id, nid, mag, xyz, cid=0, comment='')


            #print(step.cloads)
        #step.cloads


if __name__ == '__main__':   # pragma: no cover
    nastran_filename = r'C:\NASA\m4\formats\git\pyNastran\models\plate\plate.bdf'
    abaqus_inp_filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\abaqus\plate.inp'
    nastran_to_abaqus_filename(nastran_filename, abaqus_inp_filename)

    #model = read_abaqus(abaqus_filename, debug=True)
    nastran_filename_out = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\abaqus\plate2.bdf'
    abaqus_to_nastran_filename(abaqus_inp_filename, nastran_filename_out)
    x = 1
    #nastran_filename = 'junk.bdf'

    #abaqus_to_nastran_filename(abaqus_filename, nastran_filename)
    #nastran_filename = 'junk.bdf'
