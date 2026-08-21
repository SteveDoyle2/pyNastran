from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
if TYPE_CHECKING:  # pramga: no cover
    from pyNastran.op2.op2 import OP2

try:
    from tables import File, Int64Col, Float64Col, StringCol
    IS_PYTABLES = True
except ImportError:
    IS_PYTABLES = False


def get_h5_elemental_nodal(model: OP2):
    real_h5_stress_rod_dict = {
        'EID': Int64Col(pos=0),
        'A': Float64Col(pos=1),
        'MSA': Float64Col(pos=2),
        'T': Float64Col(pos=3),
        'MST': Float64Col(pos=4),
        'DOMAIN_ID': Int64Col(pos=5),
    }

    stress = model.op2_results.stress
    elemental_dicts = [
        # TODO: doesn't support complex/random
        (('STRESS', 'ROD'), stress.crod_stress, real_h5_stress_rod_dict),
    ]

    nodal_dicts = []
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
                        'GRID_POINT_FORCE', 'GRID_POINT_FORCE_CPLX', '')

    assert len(model.displacements) + len(model.eigenvectors) > 0, len(nodal_dicts)
    assert len(nodal_dicts) > 0, nodal_dicts

    elemental_dicts = [(name, dicti, table_dicti)
                       for name, dicti, table_dicti in elemental_dicts if len(dicti)]

    key_to_id_map = []
    for name, dicts, table_dicti in elemental_dicts:
        for obj_key, obj in dicts.items():
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
                        name_real: str='',
                        name_imag: str='',
                        name_random: str='') -> list[tuple]:
    reals = []
    imags = []
    randoms = []
    for key, table in tables_dict.items():
        if table.analysis_code in {1, 2, 6}:
            # 1: statics
            # 2: modes
            # 6: time
            reals.append((key, table))
        elif table.analysis_code == 5:
            # 5: freq
            imags.append((key, table))
        else:
            raise NotImplementedError(table)
            # 7: pre-buckling
            # 8: post-buckling
            # 9: complex eigenvalues

    if reals:
        assert len(name_real), name_real
        key0, table0 = reals[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_real, reals, h5_table_dict))
    if imags:
        assert len(name_imag), name_imag
        key0, table0 = imags[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_imag, imags, h5_table_dict))
    if randoms:
        assert len(name_random), name_random
        key0, table0 = randoms[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_random, randoms, h5_table_dict))
    return nodal_dicts

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
        for freq in obj.freq:
            key = (subcase_id, step, analysis_code, float(freq), eigi, mode,
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
    elif analysis_code == 6:  # time
        for time in obj.times:
            key = (subcase_id, step, analysis_code, float(time), eigi, mode,
                   design_cycle, random, se,
                   afpm, trmc, instance, module, substep, impfid, ndomains)
            keys.append(key)
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
     - single subcase nodal results (not grid_point_forces)
     - index for nodal results

    doesn't handle:
     - multiple subcases
     - index for elemental results
     - modal/transient/buckling/freq for stress/strain/forces/grid_point_forces/strain_energy
     - optimization
    """
    # nastran_group = h5file.create_group('/', 'NASTRAN')
    result_group = h5file.create_group(nastran_group, 'RESULT')

    index_group = h5file.create_group(root, 'INDEX')
    nastran_index_group = h5file.create_group(index_group, 'NASTRAN')
    nastran_index_result_group = h5file.create_group(nastran_index_group, 'RESULT')
    write_h5_domain(h5file, result_group, key_to_id_map)

    write_elemental_dicts(elemental_dicts, key_to_id_map, h5file, result_group, nastran_index_result_group)
    write_nodal_dicts(nodal_dicts, key_to_id_map, h5file, result_group, nastran_index_result_group)

def write_elemental_dicts(elemental_dicts: list[tuple],
                          key_to_id_map,
                          h5file: File,
                          result_group, index_group):
    if len(elemental_dicts) == 0:
        return
    elemental_group_ = h5file.create_group(result_group, 'ELEMENTAL')
    elemental_groups = defaultdict(list)
    for (group_name, name), objs, table_dicti in elemental_dicts:
        elemental_groups[group_name].append((name, objs, table_dicti))
    for group_name, element_groups_ in elemental_groups.items():
        for (name, objs, table_dicti) in element_groups_:
            # print(f'adding {group_name} / {name}')
            elemental_group = h5file.create_group(elemental_group_, group_name)
            table = h5file.create_table(elemental_group, name, table_dicti)
            flag = (group_name, name)

            for keyi, obj in objs.items():
                domain_keys = obj_to_domain_key(obj)
                domain_key = domain_keys[0]
                idomain0 = key_to_id_map.index(domain_key) + 1

                data = obj.data
                ntime, neid = data.shape[:2]
                assert ntime == 1, ntime
                arr = np.empty(neid, dtype=table.dtype)
                if flag == ('STRESS', 'ROD'):
                    arr["EID"] = obj.element
                    arr["A"] = data[0, :, 0]
                    arr["MSA"] = data[0, :, 1]
                    arr["T"] = data[0, :, 2]
                    arr["MST"] = data[0, :, 3]
                else:  # pragma: no cover0
                    raise NotImplementedError(flag)
                arr["DOMAIN_ID"] = np.full(neid, idomain0, dtype='int64')
                table.append(arr)
                table.flush()
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

        table = h5file.create_table(nodal_group, name, table_dicti)
        table_index = h5file.create_table(nodal_index_group, name, domain_table_dicti)

        ntime = 0
        ntime_nnode = 0
        if name == 'GRID_POINT_FORCE':
            for key, obj in key_obj_tuple:
                data = obj.data
                ntimei, nnodei = data.shape[:2]
                assert ntimei == 1, data.shape  # TODO: limited to statics
                ntime += ntimei
                ntime_nnode += nnodei * ntimei
        else:
            for key, obj in key_obj_tuple:
                data = obj.data
                ntimei, nnodei = data.shape[:2]
                ntime += ntimei
                ntime_nnode += nnodei * ntimei

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
                itime = 0
                idomain = idomain0 + itime
                # print(f'idomain: {idomain}')
                # assert idomain < 20, idomain

                ntime_nnode1 = ntime_nnode0 + nnode
                assert ntime == 1, data.shape  # TODO: limited to statics
                arr["ID"] = obj.node_element[0, :, 0]
                arr["EID"] = obj.node_element[0, :, 1]
                arr["ELNAME"] = obj.element_names[0, :]
                arr["F1"] = data[0, :, 0]
                arr["F2"] = data[0, :, 1]
                arr["F3"] = data[0, :, 2]
                arr["M1"] = data[0, :, 3]
                arr["M2"] = data[0, :, 4]
                arr["M3"] = data[0, :, 5]
                arr["DOMAIN_ID"] = np.full(nnode, idomain, dtype='int64')

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
                    arr["ID"][ntime_nnode0:ntime_nnode1] = obj.node_gridtype[:, 0]
                    arr["X"][ntime_nnode0:ntime_nnode1] = data[itime, :, 0]
                    arr["Y"][ntime_nnode0:ntime_nnode1] = data[itime, :, 1]
                    arr["Z"][ntime_nnode0:ntime_nnode1] = data[itime, :, 2]
                    arr["RX"][ntime_nnode0:ntime_nnode1] = data[itime, :, 3]
                    arr["RY"][ntime_nnode0:ntime_nnode1] = data[itime, :, 4]
                    arr["RZ"][ntime_nnode0:ntime_nnode1] = data[itime, :, 5]
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
