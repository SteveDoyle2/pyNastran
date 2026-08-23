from __future__ import annotations
from itertools import count
from collections import defaultdict
import warnings
from typing import TYPE_CHECKING

import numpy as np
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

try:
    from tables import File, Int64Col, Float64Col, StringCol
    IS_PYTABLES = True
except ImportError:
    IS_PYTABLES = False


def get_h5_elemental_nodal(model: OP2):
    nodal_dicts = []
    elemental_dicts = []

    stress = model.op2_results.stress
    stress.get_h5_elemental_tables(elemental_dicts)

    strain = model.op2_results.strain
    strain.get_h5_elemental_tables(elemental_dicts)

    force = model.op2_results.force
    force.get_h5_elemental_tables(elemental_dicts)

    modal_contribution = model.op2_results.modal_contribution
    modal_contribution.get_h5_nodal_tables(nodal_dicts)
    modal_contribution.get_h5_elemental_tables(elemental_dicts)

    # only tested for real
    split_table_by_type(nodal_dicts, model.displacements,
                        'DISPLACEMENT', 'DISPLACEMENT_CPLX', 'DISPLACEMENT_RANDOM')
    split_table_by_type(nodal_dicts, model.velocities,
                        'VELOCITY', 'VELOCITY_CPLX', 'VELOCITY_RANDOM')
    split_table_by_type(nodal_dicts, model.accelerations,
                        'ACCELERATION', 'ACCELERATION_CPLX', 'ACCELERATION_RANDOM')
    split_table_by_type(nodal_dicts, model.eigenvectors,
                        'EIGENVECTOR', 'EIGENVECTOR_CPLX', 'EIGENVECTOR_RANDOM')
    split_table_by_type(nodal_dicts, model.load_vectors,
                        'APPLIED_LOAD', 'APPLIED_LOAD_CPLX', 'APPLIED_LOAD_RANDOM')
    split_table_by_type(nodal_dicts, model.mpc_forces,
                        'MPC_FORCE', 'MPC_FORCE_CPLX', 'MPC_FORCE_RANDOM')
    split_table_by_type(nodal_dicts, model.spc_forces,
                        'SPC_FORCE', 'SPC_FORCE_CPLX', 'SPC_FORCE_RANDOM')
    split_table_by_type(nodal_dicts, model.grid_point_forces,
                        'GRID_POINT_FORCE', '', '')

    # assert len(model.displacements) + len(model.eigenvectors) > 0, len(nodal_dicts)
    assert len(nodal_dicts) > 0, nodal_dicts

    elemental_dicts = [(name, dicti, table_dicti)
                       for name, dicti, table_dicti in elemental_dicts if len(dicti)]

    key_to_id_map = []
    for name, key_obj_tuple, table_dicti in elemental_dicts:
        for obj_key, obj in key_obj_tuple:
            keys = obj_to_domain_key(obj)
            for key in keys:
                if key not in key_to_id_map:
                    key_to_id_map.append(key)
    for name, key_obj_tuple, table_dicti in nodal_dicts:
        for obj_key, obj in key_obj_tuple:
            # print('nodal', obj)
            keys = obj_to_domain_key(obj)
            for key in keys:
                if key not in key_to_id_map:
                    key_to_id_map.append(key)
    return elemental_dicts, nodal_dicts, key_to_id_map

def split_table_by_type(nodal_dicts: list[tuple],
                        tables_dict: dict,
                        name_real: str | tuple[str, str]='',
                        name_imag: str | tuple[str, str]='',
                        name_random: str | tuple[str, str]='') -> list[tuple]:
    """breaks the tables into separate blocks based on result type (e.g., real vs. imag)"""
    reals = []
    imags = []
    randoms = []
    for key, table in tables_dict.items():
        # print(key, table.analysis_code)
        if table.analysis_code in {1, 2, 6}:
            # 1: statics
            # 2: modes
            # 6: time
            reals.append((key, table))
        elif table.analysis_code == 5:
            # 5: freq
            # print('5, freq, imag')
            imags.append((key, table))
        elif table.analysis_code == 8:
            # 5: post-buckling
            reals.append((key, table))
        else:  # pragma: no cover
            raise NotImplementedError(table.analysis_code)
            # 7: pre-buckling
            # 8: post-buckling
            # 9: complex eigenvalues

    if reals:
        assert len(name_real), name_real
        key0, table0 = reals[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_real, reals, h5_table_dict))
    elif reals:
        warnings.warn(f'missing {name_real}')

    if imags and len(name_imag):
        assert len(name_imag), name_imag
        key0, table0 = imags[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_imag, imags, h5_table_dict))
    elif imags:
        warnings.warn(f'missing {name_imag}')

    if randoms and len(name_random):
        assert len(name_random), name_random
        key0, table0 = randoms[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_random, randoms, h5_table_dict))
    return nodal_dicts

def split_quad_table_by_type(elemental_dicts,
                             my_dict: dict, result_group: str, word: str):
    # QUAD_CN vs ???
    reals_centroid = []
    reals_corner = []
    for key, table in my_dict.items():
        if table.nnodes_per_element == 1:
            reals = reals_centroid
        else:
            reals = reals_corner

        if table.analysis_code in {1, 2, 6, 8}:
            # 1: statics
            # 2: modes
            # 6: time
            # 8: post-buckling
            reals.append((key, table))
        else:  # pragma: no cover
            raise NotImplementedError(table.analysis_code)

    if len(reals_centroid):
        name_real = (result_group, 'QUAD_CEN')
        key0, table0 = reals_corner[0]
        h5_table_dict = table0.h5_table_dict()
        elemental_dicts.append((name_real, reals_corner, h5_table_dict))
    if len(reals_corner):
        name_real = (result_group, 'QUAD_CN')  # corner
        key0, table0 = reals_corner[0]
        h5_table_dict = table0.h5_table_dict()
        elemental_dicts.append((name_real, reals_corner, h5_table_dict))
    return elemental_dicts

def obj_to_domain_key(obj) -> list[tuple]:
    # print(obj.get_stats())
    subcase_id = obj.isubcase
    analysis_code = obj.analysis_code
    ndomains = obj.data.shape[0]
    step = 0.
    time = 0.
    eigi = 0.
    mode = 0
    design_cycle = 0
    random = 0
    se = 0
    afpm = 0
    trmc = 0
    instance = 0
    module = 0
    substep = 0
    impfid = 0
    keys = []
    if analysis_code == 1:  # statics
        key = (subcase_id, step, analysis_code, time, eigi, mode,
               design_cycle, random, se,
               afpm, trmc, instance, module, substep, impfid, ndomains)
        keys.append(key)
    elif analysis_code == 2:  # modes
        for mode, eign in zip(obj.modes, obj.eigns):
            key = (subcase_id, step, analysis_code, float(eign), eigi, int(mode),
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
    elif analysis_code == 5:  # freq
        for freq in obj.freqs:
            key = (subcase_id, step, analysis_code, float(freq), eigi, mode,
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
    elif analysis_code == 6:  # time
        for time in obj._times:
            key = (subcase_id, step, analysis_code, float(time), eigi, mode,
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
    elif analysis_code == 8:  # post-buckling
        for mode, eigr in zip(count(), obj.eigrs):
            key = (subcase_id, step, analysis_code, float(eigr), eigi, int(mode),
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
        # raise NotImplementedError(obj.get_stats())
    elif analysis_code == 9:  # complex modes
        for mode, eigr, eigi in zip(obj.modes, obj.eigrs, obj.eigis):
            key = (subcase_id, step, analysis_code, float(eigr), float(eigi), int(mode),
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
    else:
        raise NotImplementedError(obj.get_stats())
    # ID, SUBCASE, STEP, ANALYSIS, TIME_FREQ_EIGR, EIGI, MODE, DESIGN_CYCLE, RANDOM, SE,
    #     AFPM, TRMC, INSTANCE, MODULE, SUBSTEP, IMPFID,
    return keys


def write_h5_domain(h5file: File, result_group, key_to_id_map):
    domain_dict = {
        'ID': Int64Col(pos=0),
        'SUBCASE': Int64Col(pos=1),
        'STEP': Int64Col(pos=2),
        'ANALYSIS': Int64Col(pos=3),
        'TIME_FREQ_EIGR': Float64Col(pos=4),
        'EIGI': Float64Col(pos=5),
        'MODE': Int64Col(pos=6),
        'DESIGN_CYCLE': Int64Col(pos=7),
        'RANDOM': Int64Col(pos=8),
        'SE': Int64Col(pos=9),
        # AFPM   Indicates the beginning of an acoustic field point mesh Bulk Data Section.
        # AFPMID Acoustic field point mesh identification number
        'AFPM': Int64Col(pos=10),
        # TRMC stands for Trim Component
        'TRMC': Int64Col(pos=11),
        'INSTANCE': Int64Col(pos=12),
        'MODULE': Int64Col(pos=13),
        'SUBSTEP': Int64Col(pos=14),
        # IMPFID stands for Imperfection Case ID
        'IMPFID': Int64Col(pos=15),
    }
    ndomain = len(key_to_id_map)

    domain_table = h5file.create_table(result_group, 'DOMAINS', domain_dict)
    arr = np.empty(ndomain, dtype=domain_table.dtype)
    arr["ID"] = np.arange(ndomain) + 1
    arr["SUBCASE"] = np.array([val[0] for val in key_to_id_map])
    arr["STEP"] = np.array([val[1] for val in key_to_id_map])
    arr["ANALYSIS"] = np.array([val[2] for val in key_to_id_map])
    arr["TIME_FREQ_EIGR"] = np.array([val[3] for val in key_to_id_map])
    arr["EIGI"] = np.array([val[4] for val in key_to_id_map])
    arr["MODE"] = np.array([val[5] for val in key_to_id_map])
    arr["DESIGN_CYCLE"] = np.array([val[6] for val in key_to_id_map])
    arr["RANDOM"] = np.array([val[7] for val in key_to_id_map])
    arr["SE"] = np.array([val[8] for val in key_to_id_map])
    arr["AFPM"] = np.array([val[9] for val in key_to_id_map])
    arr["TRMC"] = np.array([val[10] for val in key_to_id_map])
    arr["INSTANCE"] = np.array([val[11] for val in key_to_id_map])
    arr["MODULE"] = np.array([val[12] for val in key_to_id_map])
    arr["SUBSTEP"] = np.array([val[13] for val in key_to_id_map])
    arr["IMPFID"] = np.array([val[14] for val in key_to_id_map])
    # ID, SUBCASE, STEP, ANALYSIS, TIME_FREQ_EIGR, EIGI, MODE, DESIGN_CYCLE, RANDOM, SE,
    #     AFPM, TRMC, INSTANCE, MODULE, SUBSTEP, IMPFID,
    domain_table.append(arr)


def write_h5_results(model: OP2, h5file: File,
                     nastran_group, key_to_id_map,
                     elemental_dicts, nodal_dicts,
                     root: str='/'):
    """
    supports:
     - domains support
     - nodal/elemental results
     - nodal/elemental index support

    doesn't handle:
     - modal/transient/buckling/freq for grid_point_forces/strain_energy
     - imaginary/random elemental results (stress/strain/force/strain_energy)
     - optimization
     - matrices
     - trim
     - flutter

    not sure if supported:
     - multiple subcases

    real result types supported:
     - nodal: displacement, velocity, acceleration, load_vector, spc/mpc forces, grid point forces
     - elemental stress/strain/force:
       - crod
     - elemental stress/strain
       - ctria3, cquad4 (corner), composite ctria3/cquad4
       - ctetra, cpenta, chexa
     - strain_energy: N/A
    """
    # nastran_group = h5file.create_group('/', 'NASTRAN')
    result_group = h5file.create_group(nastran_group, 'RESULT')

    index_group = h5file.create_group(root, 'INDEX')
    nastran_index_group = h5file.create_group(index_group, 'NASTRAN')
    nastran_index_result_group = h5file.create_group(nastran_index_group, 'RESULT')
    write_h5_domain(h5file, result_group, key_to_id_map)

    write_elemental_dicts(elemental_dicts, key_to_id_map, h5file, result_group, nastran_index_result_group)
    write_nodal_dicts(nodal_dicts, key_to_id_map, h5file, result_group, nastran_index_result_group)
    write_summary(model, h5file, result_group, nastran_index_result_group)

def write_summary(model: OP2, h5file: File, result_group, index_group):
    is_summary = False
    domain_table_dicti = {
        "DOMAIN_ID": Int64Col(pos=0),
        "POSITION": Int64Col(pos=1),
        "LENGTH": Int64Col(pos=2),
    }

    if len(model.eigenvalues):
        assert len(model.eigenvalues) == 1, model.eigenvalues
        for key, obj in model.eigenvalues.items():
            if not is_summary:
                summary = h5file.create_group(result_group, 'SUMMARY')
                summary_index = h5file.create_group(index_group, 'SUMMARY')
                is_summary = True

            name = 'EIGENVALUE'
            h5_table_dict = obj.h5_table_dict()
            table = h5file.create_table(summary, name, h5_table_dict)
            table_index = h5file.create_table(summary_index, name, domain_table_dicti)

            nmode = len(obj.mode)
            arr = np.empty(nmode, dtype=table.dtype)
            arr_index = np.empty(nmode, dtype=table_index.dtype)

            idomain0 = 1
            domains = idomain0 + obj.mode

            obj.add_to_h5_array(arr)
            arr["DOMAIN_ID"] = domains

            # domain
            arr_index["DOMAIN_ID"] = 0
            arr_index["POSITION"] = 0
            arr_index["LENGTH"] = nmode


def write_elemental_dicts(elemental_dicts: list[tuple],
                          key_to_id_map,
                          h5file: File,
                          result_group, index_group):
    if len(elemental_dicts) == 0:
        return

    domain_table_dicti = {
        "DOMAIN_ID": Int64Col(pos=0),
        "POSITION": Int64Col(pos=1),
        "LENGTH": Int64Col(pos=2),
    }
    elemental_group_ = h5file.create_group(result_group, 'ELEMENTAL')
    elemental_index_group_ = h5file.create_group(index_group, 'ELEMENTAL')

    elemental_groups = defaultdict(list)
    for (group_name, name), key_obj_tuple, table_dicti in elemental_dicts:
        #print(f'adding group={group_name} name={name}')
        elemental_groups[group_name].append((name, key_obj_tuple, table_dicti))

    for group_name, element_groups_ in elemental_groups.items():
        elemental_group = h5file.create_group(elemental_group_, group_name)
        elemental_index_group = h5file.create_group(elemental_index_group_, group_name)
        for (name, key_obj_tuple, table_dicti) in element_groups_:
            print(f'adding {group_name} / {name}')
            flag = (group_name, name)

            # the table has to be out here in order to handle multi-subcase
            table = h5file.create_table(elemental_group, name, table_dicti)
            table_index = h5file.create_table(elemental_index_group, name, domain_table_dicti)

            ntime, ntime_neid = get_ntime_neid(name, key_obj_tuple)
            arr = np.empty(ntime_neid, dtype=table.dtype)
            arr_index = np.empty(ntime, dtype=table_index.dtype)

            ntime_neid0 = 0
            for keyi, obj in key_obj_tuple:
                domain_keys = obj_to_domain_key(obj)
                domain_key = domain_keys[0]
                idomain0 = key_to_id_map.index(domain_key) + 1

                data = obj.data
                ntime = data.shape[0]
                neid = get_neid(obj)
                for itime in range(ntime):
                    idomain = idomain0 + itime
                    ntime_neid1 = ntime_neid0 + neid
                    # print(f'idomain={idomain} position={ntime_neid0} length={ntime_neid1-ntime_neid0} neid={neid}')
                    # print(f'ntime_neid0={ntime_neid0} ntime_neid1={ntime_neid1} nelements={len(obj.element)}')
                    assert ntime_neid1 > ntime_neid0
                    obj.add_to_h5_array(arr, ntime_neid0, ntime_neid1, itime)
                    arr["DOMAIN_ID"][ntime_neid0:ntime_neid1] = np.full(neid, idomain0, dtype='int64')

                    # domain
                    arr_index["DOMAIN_ID"][itime] = idomain
                    arr_index["POSITION"][itime] = ntime_neid0
                    arr_index["LENGTH"][itime] = ntime_neid1 - ntime_neid0

                table.append(arr)
                table.flush()
                table_index.append(arr_index)
                table_index.flush()
    return

def write_nodal_dicts(nodal_dicts: list[tuple],
                      key_to_id_map,
                      h5file: File,
                      result_group, index_group):
    if len(nodal_dicts) == 0:
        return
    domain_table_dicti = {
        "DOMAIN_ID": Int64Col(pos=0),
        "POSITION": Int64Col(pos=1),
        "LENGTH": Int64Col(pos=2),
    }
    nodal_group = h5file.create_group(result_group, 'NODAL')
    nodal_index_group = h5file.create_group(index_group, 'NODAL')
    for name, key_obj_tuple, table_dicti in nodal_dicts:
        # print(f'adding {name}: {len(key_obj_tuple)}')
        # for key, obj in key_obj_tuple:
        #     print(f'key = {key}')
        #     print(f'obj = {obj}')
        #     print('---------------------------------------')

        # the table has to be out here in order to handle multi-subcase
        table = h5file.create_table(nodal_group, name, table_dicti)
        table_index = h5file.create_table(nodal_index_group, name, domain_table_dicti)

        ntime, ntime_nnode = get_ntime_nnode(name, key_obj_tuple)

        ntime_nnode0 = 0
        arr = np.empty(ntime_nnode, dtype=table.dtype)
        arr_index = np.empty(ntime, dtype=table_index.dtype)
        for key, obj in key_obj_tuple:
            domain_keys = obj_to_domain_key(obj)
            domain_key = domain_keys[0]
            idomain0 = key_to_id_map.index(domain_key) + 1
            data = obj.data
            ntime, nnode = data.shape[:2]

            if name == 'GRID_POINT_FORCE':
                for itime in range(ntime):
                    idomain = idomain0 + itime
                    # print(f'idomain: {idomain}')
                    # assert idomain < 20, idomain

                    ntime_nnode1 = ntime_nnode0 + nnode
                    obj.add_to_h5_array(arr, ntime_nnode0, ntime_nnode1, itime)
                    arr["DOMAIN_ID"][ntime_nnode0:ntime_nnode1] = np.full(nnode, idomain, dtype='int64')

                    # domain
                    arr_index["DOMAIN_ID"][itime] = idomain
                    arr_index["POSITION"][itime] = ntime_nnode0
                    arr_index["LENGTH"][itime] = ntime_nnode1 - ntime_nnode0
                table.append(arr)
                table.flush()
                table_index.append(arr_index)
                table_index.flush()
            else:
                for itime in range(ntime):
                    idomain = idomain0 + itime
                    # assert idomain < 20, idomain
                    ntime_nnode1 = ntime_nnode0 + nnode
                    obj.add_to_h5_array(arr, ntime_nnode0, ntime_nnode1, itime)
                    arr["DOMAIN_ID"][ntime_nnode0:ntime_nnode1] = np.full(nnode, idomain, dtype='int64')

                    # domain
                    # print(f'idomain={idomain} position={ntime_nnode0} length={ntime_nnode1-ntime_nnode0+1} nnode={nnode}')
                    arr_index["DOMAIN_ID"][itime] = idomain
                    arr_index["POSITION"][itime] = ntime_nnode0
                    arr_index["LENGTH"][itime] = ntime_nnode1 - ntime_nnode0
                    ntime_nnode0: ntime_nnode1
                    ntime_nnode0 += nnode
                # print('DOMAIN_ID', arr_index["DOMAIN_ID"])
                # print('POSITION', arr_index["POSITION"])
                table.append(arr)
                table.flush()
                table_index.append(arr_index)
                table_index.flush()
    return

def get_ntime_nnode(name: str, key_obj_tuple) -> tuple[int, int]:
    ntime = 0
    ntime_nnode = 0
    if name == 'GRID_POINT_FORCE':
        for key, obj in key_obj_tuple:
            data = obj.data
            ntimei, nnodei = data.shape[:2]
            # assert ntimei == 1, data.shape  # TODO: limited to statics
            ntime += ntimei
            ntime_nnode += nnodei * ntimei
    else:
        # this block handles multiple subcases
        for key, obj in key_obj_tuple:
            data = obj.data
            ntimei, nnodei = data.shape[:2]
            ntime += ntimei
            ntime_nnode += nnodei * ntimei
    return ntime, ntime_nnode


def get_ntime_neid(name: str, key_obj_tuple) -> tuple[int, int]:
    ntime = 0
    ntime_neid = 0
    for key, obj in key_obj_tuple:
        neidi = get_neid(obj)
        data = obj.data
        ntimei = data.shape[0]
        ntime += ntimei
        ntime_neid += neidi * ntimei
    return ntime, ntime_neid


def get_neid(obj) -> int:
    if hasattr(obj, "get_neid"):
        neid = obj.get_neid()
    elif hasattr(obj, "element"):
        neid = obj.element.shape[0]
    elif hasattr(obj, "element_cid"):
        neid = obj.element_cid.shape[0]
    elif hasattr(obj, "element_layer"):
        neid = obj.element_layer.shape[0]
    else:
        return obj.get_neid()
    return neid
