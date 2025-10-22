from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
from itertools import count
#from collections import defaultdict

import numpy as np
from tables import Group, Node

from pyNastran.op2.op2_interface.op2_classes import (
    RealDisplacementArray, RealVelocityArray, RealAccelerationArray,
    RealLoadVectorArray, RealSPCForcesArray, RealMPCForcesArray, RealEigenvectorArray,

    ComplexDisplacementArray, ComplexVelocityArray, ComplexAccelerationArray, ComplexEigenvectorArray,
    ComplexLoadVectorArray, ComplexSPCForcesArray, # ComplexMPCForcesArray
)
from .utils import break_domain_by_case, get_name, get_analysis_code_times
from .stress_strain import get_idomain

from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.utils import get_group_name, get_attributes
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.op2_vectorized3.op2_geom import OP2 # , OP2Geom
    import pandas as pd


def read_nodal_result(model: OP2, domains: np.ndarray,
                      nodal: Group, nodal_index: Group):
    iresult = 0
    results = {}
    #geom_model = None
    ids = None

    for h5_node_ in nodal._f_iter_nodes():
        if isinstance(h5_node_, Group):
            print(h5_node_)
            name = get_group_name(h5_node_)
            if name == 'ELEMENTAL':
                group1
            else:
                raise NotImplementedError(name)

        elif isinstance(h5_node_, Node):
            attrs = get_attributes(h5_node_)
            assert len(attrs) == 1, attrs
            version = attrs['version'][0]

            #data = h5_node_.read()
            #assert isinstance(data, np.ndarray), data
            name = h5_node_.name

            group = nodal[name].read()
            index = nodal_index[name].read()
            if name in {'GRID_FORCE', 'GRID_WEIGHT', 'TEMPERATURE', 'KINETIC_ENERGY'}:
                pass
            #--------------------------------------------------------------------------
            elif name == 'DISPLACEMENT':
                basename = 'displacements'
                load_displacement(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'VELOCITY':
                basename = 'velocities'
                load_displacement(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'APPLIED_LOAD':
                basename = 'load_vectors'
                load_displacement(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'SPC_FORCE':
                basename = 'spc_forces'
                load_displacement(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'MPC_FORCE':
                basename = 'mpc_forces'
                #load_displacement(
                    #basename, iresult, results, domains, group, index, version,
                    #ids, model, subcases=None)
                #--------------------------------------------------------------------------
            elif name == 'EIGENVECTOR':
                basename = 'eigenvectors'
                load_displacement(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            #--------------------------------------------------------------------------
            elif name == 'EIGENVECTOR_CPLX':
                basename = 'eigenvectors'
                load_eigenvector_complex(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'DISPLACEMENT_CPLX':
                basename = 'displacements'
                load_eigenvector_complex(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'APPLIED_LOAD_CPLX':
                basename = 'load_vectors'
                load_eigenvector_complex(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'SPC_FORCE_CPLX':
                basename = 'spc_forces'
                load_eigenvector_complex(
                    basename, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            #elif name == 'MPC_FORCE_CPLX':
                #basename = 'mpc_forces'
                #load_eigenvector_complex(
                    #basename, iresult, results, domains, group, index, ids,
                    #model, subcases=None)
            else:
                raise NotImplementedError(name)
        else:
            raise NotImplementedError(h5_node_)

def load_displacement(result_name: str,
                      iresult: int,
                      results: dict[int, Function],
                      domains_df: pd.DataFrame,
                      group: h5py._hl.dataset.Dataset,
                      index: h5py._hl.dataset.Dataset,
                      version: int,
                      ids: np.ndarray,
                      model: OP2,
                      subcases: Optional[list[int]]=None) -> int:
    assert version == 1, version
    if result_name == 'displacements':
        class_obj = RealDisplacementArray
        table_name = 'OUGV1'
    elif result_name == 'velocities':
        class_obj = RealVelocityArray
        table_name = 'OVG1'
    elif result_name == 'eigenvectors':
        class_obj = RealEigenvectorArray
        table_name = 'OUGV1'
    elif result_name == 'load_vectors':
        class_obj = RealLoadVectorArray
        table_name = 'OPG1'
    elif result_name == 'spc_forces':
        class_obj = RealSPCForcesArray
        table_name = 'OQG1'
    elif result_name == 'mpc_forces':
        class_obj = RealMPCForcesArray
        table_name = 'OQG1' # OQG1/OQMG1?
    else:
        #table_name = 'OUGV1'
        raise RuntimeError(result_name)

    basename = result_name + ' (real)'

    # dtype([('ID', '<i8'), ('X', '<f8'), ('Y', '<f8'), ('Z', '<f8'),
    #       ('RX', '<f8'), ('RY', '<f8'), ('RZ', '<f8'), ('DOMAIN_ID', '<i8')])
    NID = group['ID']
    nnodes = len(NID)
    #names = group.dtype.names
    TXYZ_RXYZ = np.stack([group['X'], group['Y'], group['Z'],
                          group['RX'], group['RY'], group['RZ'],
    ], axis=1)
    #DOMAIN = group['DOMAIN_ID']
    # what are the subcases...
    #udomains = np.unique(DOMAIN)
    #all_data = np.stack([X, Y, Z, RX, RY, RZ], axis=1, out=None)

    mode, eigr, eigi, data_by_group = _get_data_by_group_node(
        basename, domains_df, index,
        NID, TXYZ_RXYZ)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            node_gridtype = []
            datas = []
            idomain = []
            nnode_gridtype = []
            ntimes = len(group)

            #is_modes, is_freq,
            for itime, (idomaini, node_gridtypei, datai) in group.items():
                idomain.append(idomaini)
                node_gridtype = node_gridtypei
                nnode_gridtypei = node_gridtypei.shape[0]
                nnode_gridtype.append(nnode_gridtypei)
                datas.append(datai)

            assert max(nnode_gridtype) == min(nnode_gridtype)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nnodes = node_gridtype.shape[0]
            assert data.shape == (ntimes, nnodes, 6), data.shape
            #freqs = np.zeros(ntimes, dtype=TX.dtype)
            if is_static:
                disp = class_obj.add_static_case(
                    table_name, node_gridtype, data, isubcase, # freqs,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
                #model.displacements[isubcase] = disp
                #model.subcase_key[isubcase].append(isubcase)
            elif is_modes:
                modes, eigenvalues, mode_cycles = modes_data
                disp = class_obj.add_modal_case(
                    table_name, node_gridtype, data, isubcase,
                    modes, eigenvalues, mode_cycles,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                if not hasattr(class_obj, 'add_buckling_case'):
                    model.log.warning(f'skipping {result_name} - add_buckling_case')
                    continue
                modes, eigenvalues, mode_cycles = buckling_data
                disp = class_obj.add_buckling_case(
                    table_name, node_gridtype, data, isubcase,
                    modes, eigenvalues, mode_cycles,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                disp = class_obj.add_transient_case(
                    table_name, node_gridtype, data, isubcase, times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            else:
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            getattr(model, result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def load_eigenvector_complex(result_name: str,
                             iresult: int,
                             results: dict[int, Function],
                             domains_df: pd.DataFrame,
                             group: h5py._hl.dataset.Dataset,
                             index: h5py._hl.dataset.Dataset,
                             version: int,
                             ids: np.ndarray,
                             model: OP2,
                             subcases: Optional[list[int]]=None) -> int:
    assert version == 1, version
    if result_name == 'displacements':
        class_obj = ComplexDisplacementArray
        table_name = 'OUGV1'
    elif result_name == 'eigenvectors':
        class_obj = ComplexEigenvectorArray
        table_name = 'OUGV1'
    elif result_name == 'load_vectors':
        class_obj = ComplexLoadVectorArray
        table_name = 'OPGV1'
    elif result_name == 'spc_forces':
        class_obj = ComplexSPCForcesArray
        table_name = 'OQG1'
    else:
        #table_name = 'OUGV1'
        raise RuntimeError(result_name)

    basename = result_name + ' (complex)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'

    #dtype=[('ID', '<i8'),
    # ('XR', '<f8'), ('YR', '<f8'), ('ZR', '<f8'),
    #('RXR', '<f8'), ('RYR', '<f8'), ('RZR', '<f8'),
    # ('XI', '<f8'), ('YI', '<f8'), ('ZI', '<f8'),
    #('RXI', '<f8'), ('RYI', '<f8'), ('RZI', '<f8'), ('DOMAIN_ID', '<i8')])

    NID = group['ID']
    nnodes = len(NID)
    #names = group.dtype.names
    TX = group['XR'] + group['XI'] * 1j
    TY = group['YR'] + group['YI'] * 1j
    TZ = group['ZR'] + group['ZI'] * 1j
    RX = group['RXR'] + group['RXI'] * 1j
    RY = group['RYR'] + group['RYI'] * 1j
    RZ = group['RZR'] + group['RZI'] * 1j
    TXYZ_RXYZ = np.stack([TX, TY, TZ,
                          RX, RY, RZ], axis=1)
    mode, eigr, eigi, data_by_group = _get_data_by_group_node(
        basename, domains_df, index,
        NID, TXYZ_RXYZ)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for ianalysis_code, group in data_analysis.items():
            node_gridtype = []
            datas = []
            #modes = []
            #eigrs = []
            #eigis = []
            idomain = []
            nnode_gridtype = []
            ntimes = len(group)
            for itime, (idomaini, node_gridtypei, datai) in group.items():
                #(modei, eigri, eigii) = key
                #modes.append(modei)
                #eigrs.append(eigri)
                #eigis.append(eigii)
                idomain.append(idomaini)
                node_gridtype = node_gridtypei
                nnode_gridtypei = node_gridtypei.shape[0]
                nnode_gridtype.append(nnode_gridtypei)
                datas.append(datai)

            #is_static = True
            #is_modes = True
            #is_static = False
            #is_modes = False
            #is_freq = False
            #if ianalysis_code == 5: # frequency
                #times = eigr.loc[idomain].values
                #is_freq = True
            #else:
                #raise NotImplementedError(ianalysis_code)
            assert max(nnode_gridtype) == min(nnode_gridtype)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes_data) = get_analysis_code_times(
                ianalysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nnodes = node_gridtype.shape[0]
            assert data.shape == (ntimes, nnodes, 6), data.shape
            #freqs = np.zeros(ntimes, dtype=TX.dtype)
            if is_freq:
                freq = freq_data
                disp = class_obj.add_freq_case(
                    table_name, node_gridtype, data, isubcase, freq,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_complex_modes:
                #model.log.warning(f'skipping {result_name} - add_complex_modes_case')
                freqs, eigrs, eigis = complex_modes_data
                disp = class_obj.add_complex_modes_case(
                    table_name, node_gridtype, data, isubcase, freqs, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            else:
                raise NotImplementedError(ianalysis_code)
            model.subcase_key[isubcase].append(isubcase)
            getattr(model, result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _get_data_by_group_node(basename: str,
                            domains_df: pd.Dataframe, index: Group,
                            NID, DATA: np.ndarray) -> dict[Any, Any]:
    INDEX_DOMAIN = index['DOMAIN_ID']
    INDEX_POSITION = index['POSITION']
    INDEX_LENGTH = index['LENGTH']

    #DOMAIN = group['DOMAIN_ID']
    # what are the subcases...
    #udomains = np.unique(DOMAIN)
    #all_data = np.stack([X, Y, Z, RX, RY, RZ], axis=1, out=None)

    #idi, subcase, step, analysis, time_freq_eigr, eigi, mode, design_cycle,
    #random, se, afpm, trmc, instance, module, substep, impfid
    #DOMAIN_ID = domains_df[:, 0] # ['ID']
    #DOMAIN_ID = domains_df['ID']
    grouped_df = break_domain_by_case(domains_df, INDEX_DOMAIN)
    for grouped_dfi in grouped_df:
        no_domain, idomain, dfi = get_idomain(grouped_dfi, INDEX_DOMAIN)
        if no_domain:
            continue

    #for domain, position, length in zip(INDEX_DOMAIN, INDEX_POSITION, INDEX_LENGTH):
        #idomain = (DOMAIN_ID == domain)
        #mycase = domains_df.loc[idomain]
        position = INDEX_POSITION[idomain]
        length = INDEX_LENGTH[idomain]
        # ID=1; Analysis=0 -> eigenvalues
        #    ID  SUBCASE  STEP  ANALYSIS  TIME_FREQ_EIGR  EIGI  MODE  DESIGN_CYCLE  RANDOM  SE  AFPM  TRMC  INSTANCE  MODULE
        # 0   1        0     0         0    0.000000e+00   0.0     0             0       0   0     0     0         0       0
        # 1   2        1     0         2   -3.087735e-10   0.0     1             0       0   0     0     0         0       0
        # 2   3        1     0         2   -2.082743e-10   0.0     2             0       0   0     0     0         0       0
        # 3   4        1     0         2   -1.514309e-10   0.0     3             0       0   0     0     0         0       0
        subcase = dfi['SUBCASE']
        analysis_code = dfi['ANALYSIS']
        eigr = dfi['TIME_FREQ_EIGR']
        eigi = dfi['EIGI']
        mode = dfi['MODE']

        step = dfi['STEP'].values[0]
        design_cycle = dfi['DESIGN_CYCLE'].values[0]
        random = dfi['RANDOM'].values[0]
        se = dfi['SE'].values[0]
        afpm = dfi['AFPM'].values[0]
        trmc = dfi['TRMC'].values[0]
        inst = dfi['INSTANCE'].values[0]
        module = dfi['MODULE'].values[0]
        assert step in [0, 1], step
        assert design_cycle == 0, design_cycle
        assert random == 0, random
        assert se == 0, se
        assert afpm == 0, afpm
        assert trmc == 0, trmc
        assert inst == 0, inst
        assert module == 0, module

        #is_static = False
        #is_modes = np.abs(mode).max() != 0
        #is_freq = np.abs(eigr).max() != 0 or np.abs(eigi).max() != 0

        #analysis_codei = np.uniqu
        #if analysis_code == 1:
            #is_static = True
        #else:
            #raise NotImplementedError(analysis_code)

        #iresults = np.full(len(idomain), np.nan, dtype='int32')
        #assert is_modes or is_freq or is_static
        #assert xor(is_static, is_modes, is_freq), f'is_static={is_static} is_modes={is_modes} is_freq={is_freq} are both True/False...'  # xor

        data_by_group = {}
        #data_by_group[subcase][analysis_codei]
        for subcasei in subcase:
            data_by_group[subcasei] = {}
        for subcasei, analysis_codei in zip(subcase, analysis_code):
            data_by_group[subcasei][analysis_codei] = {}

        for itime, subcasei, analysis_codei, modei, eigri, eigii, idomaini, positioni, lengthi in zip(
                count(), subcase, analysis_code, mode, eigr, eigi, idomain, position, length):
            nnodesi = lengthi
            i0 = positioni
            i1 = positioni + lengthi
            #name = get_name(basename, is_freq, subcasei, analysis_codei, modei, eigri, eigii)
            #print('  ', name)
            #print(itime, domain, positioni, lengthi)

            # make it so we can determine the other "times"
            #iresults[itime] = iresult
            #_domain = DOMAIN[i0]

            nodei = NID[i0:i1]
            #gridtypei = np.zeros(nnodesi.shape, dtype=nodei.dtype)
            node_gridtype = np.ones((nnodesi, 2), dtype=nodei.dtype)
            node_gridtype[:, 0] = nodei
            datai = DATA[i0:i1]
            #op2.displacements[subcasei]

            #key = (modei, eigri, eigii)

            data_by_group[subcasei][analysis_codei][itime] = (idomaini, node_gridtype, datai) # is_modes, is_freq,

            #results[iresult] = RealVectorTable(
                #name, itime, iresult, iresults,
                #_domain, position, length,
                #NID[i0:i1],
                #TX[i0:i1], TY[i0:i1], TZ[i0:i1],
                #RX[i0:i1], RY[i0:i1], RZ[i0:i1],
                #DOMAIN[i0:i1], location='node')
            #iresult += 1
    return mode, eigr, eigi, data_by_group

def xor(*values):
    assert len(values) > 0, values
    return sum(values) == 1
