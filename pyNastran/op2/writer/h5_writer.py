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

    split_table_by_type(nodal_dicts, model.displacements,
                        'DISPLACEMENT', 'DISPLACEMENT_CPLX', '')
    split_table_by_type(nodal_dicts, model.velocities,
                        'VELOCITY', 'VELOCITY_CPLX', '')
    split_table_by_type(nodal_dicts, model.accelerations,
                        'ACCELERATION', 'ACCELERATION_CPLX', '')
    split_table_by_type(nodal_dicts, model.eigenvectors,
                        'EIGENVECTOR', 'EIGENVECTOR_CPLX', '')
    split_table_by_type(nodal_dicts, model.load_vectors,
                        'APPLIED_LOAD', 'APPLIED_LOAD_CPLX', '')
    split_table_by_type(nodal_dicts, model.mpc_forces,
                        'MPC_FORCE', 'MPC_FORCE_CPLX', '')
    split_table_by_type(nodal_dicts, model.spc_forces,
                        'SPC_FORCE', 'SPC_FORCE_CPLX', '')
    split_table_by_type(nodal_dicts, model.grid_point_forces,
                        'GRID_POINT_FORCE', 'GRID_POINT_FORCE_CPLX', '')
    split_table_by_type(nodal_dicts, model.temperatures,
                        'TEMPERATURE', '', '')

    # assert len(nodal_dicts) > 0, nodal_dicts

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

    strain_energy = model.op2_results.strain_energy
    se_reals, se_imags = strain_energy.get_h5_strain_energy_tables()

    strain_energy_dicts = []
    if len(se_reals):
        strain_energy_dicts.append(('REAL', se_reals, {}))
    if len(se_imags):
        assert len(se_imags) == 0, se_imags
        strain_energy_dicts.append(('IMAG', se_imags, {}))

    for name, key_obj_tuple, table_dicti in strain_energy_dicts:
        for obj_key, obj in key_obj_tuple:
            # print('se', obj)
            keys = obj_to_domain_key(obj)
            for key in keys:
                if key not in key_to_id_map:
                    key_to_id_map.append(key)
    
    return elemental_dicts, nodal_dicts, strain_energy_dicts, key_to_id_map

def split_table_by_type(nodal_dicts: list[tuple],
                        tables_dict: dict,
                        name_real: str | tuple[str, str]='',
                        name_imag: str | tuple[str, str]='',
                        name_random: str | tuple[str, str]='',
                        ) -> list[tuple]:
    """
    breaks the tables into separate blocks based on result
    type (e.g., real vs. imag)
    """
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
            # 8: post-buckling
            reals.append((key, table))
        elif table.analysis_code == 9:
            # 9: complex eigenvalue
            imags.append((key, table))
        else:  # pragma: no cover
            raise NotImplementedError(table.analysis_code)
        # 7: pre-buckling

    if reals and len(name_real):
        assert len(name_real), name_real
        key0, table0 = reals[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_real, reals, h5_table_dict))
    elif reals:
        warnings.warn(f'missing real {name_real}')

    if imags and len(name_imag):
        assert len(name_imag), name_imag
        key0, table0 = imags[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_imag, imags, h5_table_dict))
    elif imags:
        warnings.warn(f'missing {name_imag} (real={name_real})')

    if randoms and len(name_random):
        assert len(name_random), name_random
        key0, table0 = randoms[0]
        h5_table_dict = table0.h5_table_dict()
        nodal_dicts.append((name_random, randoms, h5_table_dict))
    elif randoms:
        warnings.warn(f'missing random {name_random} (real={name_real})')
    return nodal_dicts

def split_quad_table_by_type(elemental_dicts,
                             my_dict: dict, result_group: str, word: str):
    # QUAD_CN vs ???
    reals_centroid = []
    reals_corner = []

    imags_centroid = []
    imags_corner = []

    for key, table in my_dict.items():
        if table.nnodes_per_element == 1:
            reals = reals_centroid
            imags = imags_centroid
        else:
            reals = reals_corner
            imags = imags_corner

        if table.analysis_code in {1, 2, 6, 8}:
            # 1: statics
            # 2: modes
            # 6: time
            # 8: post-buckling
            reals.append((key, table))
        elif table.analysis_code in {9}:
            # complex eigenvalues
            imags.append((key, table))
        else:  # pragma: no cover
            raise NotImplementedError(table.analysis_code)

    if len(reals_centroid):
        name_real = (result_group, 'QUAD_CEN')
        key0, table0 = reals_centroid[0]
        h5_table_dict = table0.h5_table_dict()
        elemental_dicts.append((name_real, reals_corner, h5_table_dict))
    if len(reals_corner):
        name_real = (result_group, 'QUAD_CN')  # corner
        key0, table0 = reals_corner[0]
        h5_table_dict = table0.h5_table_dict()
        elemental_dicts.append((name_real, reals_corner, h5_table_dict))

    if len(imags_centroid):
        name_imag = (result_group, 'QUAD_CEN_CPLX')
        key0, table0 = imags_centroid[0]
        h5_table_dict = table0.h5_table_dict()
        elemental_dicts.append((name_imag, reals_corner, h5_table_dict))
    if len(imags_corner):
        name_imag = (result_group, 'QUAD_CN_CPLX')  # corner
        key0, table0 = imags_corner[0]
        h5_table_dict = table0.h5_table_dict()
        elemental_dicts.append((name_imag, reals_corner, h5_table_dict))
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
    else:  # pragma: no cover
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
                     elemental_dicts, nodal_dicts, se_dicts,
                     root: str='/'):
    """
    supports:
     - static, modal, transient, buckling, freq, thermal
     - domains support
     - nodal/elemental results
     - eigenvalues
       - modal, buckling, complex
     - strain energy

    doesn't handle:
     - thermal forces
     - elastic/plastic/thermal stress/strain
     - random elemental results (stress/strain/force/strain_energy)
     - optimization
     - matrices
     - trim
     - flutter
     - grid point weight

    not sure if supported:
     - multiple subcases

    real/complex result types supported:
     - nodal:
       - displacement, velocity, acceleration
       - load_vector, spc/mpc forces, eigenvector
       - grid point forces
         - might have transient issues
     - elemental stress/strain/force:
       - 1d
         - crod, ctube, conrod, cvisc
         - celas1-4, cdamp1-4
         - no cbush
         - no cbar, cbeam, cbend
         - no cshear
       - shells
         - isotropic/composite
         - ctria3, cquad4 (corner), ctria6, ctriar, cquad8, cquadr
         - no cquad4 (centroid)
       - solids:
         - ctetra, cpenta, chexa
    """
    assert isinstance(root, str), root
    # nastran_group = h5file.create_group('/', 'NASTRAN')
    result_group = h5file.create_group(nastran_group, 'RESULT')

    index_group = h5file.create_group(root, 'INDEX')
    nastran_index_group = h5file.create_group(index_group, 'NASTRAN')
    nastran_index_result_group = h5file.create_group(nastran_index_group, 'RESULT')
    write_h5_domain(h5file, result_group, key_to_id_map)

    write_elemental_dicts(elemental_dicts, key_to_id_map, h5file, result_group, nastran_index_result_group)
    write_nodal_dicts(nodal_dicts, key_to_id_map, h5file, result_group, nastran_index_result_group)
    write_summary(model, h5file, result_group, nastran_index_result_group)
    write_strain_energy(
        model, h5file, result_group, nastran_index_result_group,
        se_dicts, key_to_id_map)


def write_strain_energy(model: OP2, h5file: File,
                        result_group, index_group,
                        se_dicts, key_to_id_map):
    """
    done
     - sort like MSC (loop over modes and then element types)
    TODO:
     - verify ident
     - verify domain
    """
    if len(se_dicts) == 0:
        return

    energy = h5file.create_group(result_group, 'ENERGY')
    energy_index = h5file.create_group(index_group, 'ENERGY')

    domain_table_dicti = {
        "DOMAIN_ID": Int64Col(pos=0),
        "POSITION": Int64Col(pos=1),
        "LENGTH": Int64Col(pos=2),
    }
    ident_dict = {
        "IDENT": Int64Col(pos=0),
        "ELNAME": StringCol(8, pos=1),
        "ETOTAL": Float64Col(pos=2),
        "CVALRES": Float64Col(pos=3),
        "ESUBT": Float64Col(pos=4),
        "ETOTPOS": Float64Col(pos=5),
        "ETOTNEG": Float64Col(pos=6),
    }
    energy_dict = {
        "ID": Int64Col(pos=0),
        "ENERGY": Float64Col(pos=1),
        "PCT": Float64Col(pos=2),
        "DEN": Float64Col(pos=3),
        "IDENT": Int64Col(pos=4),
        "DOMAIN_ID": Int64Col(pos=5),
        #"ELNAME": StringCol(8, pos=6),
    }

    assert len(se_dicts) == 1, se_dicts

    for group_name, key_obj_tuple, fake_dict in se_dicts:
        table_ident  = h5file.create_table(energy, 'IDENT', ident_dict)
        table_energy = h5file.create_table(energy, 'STRAIN_ELEM', energy_dict)
        table_index  = h5file.create_table(energy_index, 'STRAIN_ELEM', domain_table_dicti)

        neids = 0
        nelem_types_nmodes = 0
        
        ntimes = 0
        for key, obj in key_obj_tuple:
            #print(f'nelem_types_nmodes {obj.element_name} = {nelem_types_nmodes}')
            #print(obj.get_stats())
            ntime, neid, _nresult = obj.data.shape
            if ntimes == 0:
                ntimes = ntime
            else:
                assert ntime == ntimes, 'verifying same number of times/modes'
            
            assert isinstance(obj.element_name, str), obj.element_name
            for itime in range(ntime):
                element = obj.element[itime, :]
                #ilast = np.where(element == 100000000)[0][0]
                #neids += ilast
                neids += neid
            nelem_types_nmodes += ntime

        #print(f'{group_name} strain_energy neids={neids}')
        arr_ident  = np.empty(nelem_types_nmodes, dtype=table_ident.dtype)
        arr_energy = np.empty(neids, dtype=table_energy.dtype)
        arr_index  = np.empty(nelem_types_nmodes, dtype=table_index.dtype)

        i = 0
        eid0 = 0
        iposition = 1

        key0, obj0 = key_obj_tuple[0]
        domain_keys = obj_to_domain_key(obj0)
        domain_key = domain_keys[0]
        idomain0 = key_to_id_map.index(domain_key) + 1

        nelement_types = len(key_obj_tuple)
        for itime in range(ntime):
            arr_index["DOMAIN_ID"][itime] = idomain0 + itime
            arr_index["POSITION"][itime] = iposition
            arr_index["LENGTH"][itime] = nelement_types

            for key, obj in key_obj_tuple:
                #print(f'nelem_types_nmodes {obj.element_name} = {i}')
                #print(obj.get_stats())
                ntime, neid, nresult = obj.data.shape

                #----------------
                element = obj.element[itime, :]
                #print(f'element = {element}')
                ilast = np.where(element == 100000000)[0][0]
                eid1 = eid0 + neid
                #print(f'ilast = {ilast}')
                #print(obj.data[itime, :, :])
                assert isinstance(i, int), i
                arr_ident["IDENT"][i] = i + 1
                arr_ident["ELNAME"][i] = obj.element_name
                arr_ident["ETOTAL"][i] = obj.data[itime, ilast, 0]
                arr_ident["CVALRES"][i] = obj.cvalres
                arr_ident["ESUBT"][i] = 0.0   # TODO: what is this?
                arr_ident["ETOTPOS"][i] = obj.etotpos
                arr_ident["ETOTNEG"][i] = obj.etotneg

                #print('eids', element)
                arr_energy["ID"][eid0:eid1] = element
                arr_energy["ENERGY"][eid0:eid1] = obj.data[itime, :, 0]
                arr_energy["PCT"][eid0:eid1] = obj.data[itime, :, 1]
                arr_energy["DEN"][eid0:eid1] = obj.data[itime, :, 2]
                arr_energy["IDENT"][eid0:eid1] = i + 1
                arr_energy["DOMAIN_ID"][eid0:eid1] = idomain0 + itime
                #arr_energy["ELNAME"][eid0:eid1] = obj.element_name
                eid0 = eid1
                i += 1
            iposition += nelement_types
        
        assert i == nelem_types_nmodes, f'nelem_types_nmodes={nelem_types_nmodes} i={i}'
        assert eid0 == neids, f'neids={neids} eid0={eid0}'

        table_ident.append(arr_ident)
        table_ident.flush()
        table_energy.append(arr_energy)
        table_energy.flush()
        table_index.append(arr_index)
        table_index.flush()

        #print(list(ident_dict))
        #print(arr_ident)
        #print(list(ident_dict))
        #print('------------------------')
        #print(list(energy_dict))
        #print(arr_energy)
        #print(list(energy_dict))
    return

def write_summary(model: OP2, h5file: File, result_group, index_group):
    """
    - modal eigenvalues
    - buckling eigenvalues
    - complex eigenvalues
    TODO: modal/buckling use the same flag, but different headers
      - need an optimization case to test
    """
    if len(model.eigenvalues) == 0:
        return

    domain_table_dicti = {
        "DOMAIN_ID": Int64Col(pos=0),
        "POSITION": Int64Col(pos=1),
        "LENGTH": Int64Col(pos=2),
    }
    modes = []
    buckling = []
    imags = []
    assert len(model.eigenvalues) == 1, model.eigenvalues
    for key, obj in model.eigenvalues.items():
        if obj.is_real:
            modes.append((key, obj))
        if obj.is_buckling:
            buckling.append((key, obj))
        if obj.is_complex:
            imags.append((key, obj))

    eigenvalue_objs = []
    if len(modes):
        #print('eigenvalue modes')
        assert len(buckling) == 0, buckling  # shares common name with modes
        eigenvalue_objs.append(('EIGENVALUE', modes))
    if len(buckling):
        #print('eigenvalue buckling')
        assert len(modes) == 0, modes  # shares common name with buckling
        eigenvalue_objs.append(('EIGENVALUE', buckling))
    if len(imags):
        #print('eigenvalue imags')
        eigenvalue_objs.append(('EIGENVALUE_CPLX', imags))

    is_summary = False
    for name, key_obj_tuple in eigenvalue_objs:
        #print(name)
        assert len(key_obj_tuple) == 1, key_obj_tuple
        key, obj = key_obj_tuple[0]
        if not is_summary:
            summary = h5file.create_group(result_group, 'SUMMARY')
            summary_index = h5file.create_group(index_group, 'SUMMARY')
            is_summary = True

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

        table.append(arr)
        table.flush()
        table_index.append(arr_index)
        table_index.flush()
    return

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
                    #print(f'idomain={idomain} position={ntime_neid0} length={ntime_neid1-ntime_neid0} neid={neid}')
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
