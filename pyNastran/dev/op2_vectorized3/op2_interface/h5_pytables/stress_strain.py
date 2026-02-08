from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
from itertools import count

import numpy as np
#from numpy.lib import recfunctions as rfn  # random numpy package to merge arrays...
from tables import Group, Node

from pyNastran.op2.op2_interface.op2_classes import (
    RealRodForceArray, RealRodStressArray, RealRodStrainArray,
    ComplexRodForceArray, ComplexRodStressArray, ComplexRodStrainArray,

    RealCBarForceArray, RealBarStressArray, RealBarStrainArray,
    RealCBeamForceArray, RealBeamStressArray, RealBeamStrainArray,

    RealSpringForceArray, RealSpringStressArray, RealSpringStrainArray,
    ComplexSpringForceArray, ComplexSpringStressArray, ComplexSpringStrainArray,
    RealDamperForceArray, ComplexDamperForceArray,

    RealViscForceArray, ComplexViscForceArray,

    RealCShearForceArray, RealShearStressArray, RealShearStrainArray,
    ComplexShearStressArray, ComplexShearStrainArray, ComplexCShearForceArray,

    RealPlateStressArray, RealPlateStrainArray, RealPlateForceArray, RealPlateBilinearForceArray,
    ComplexPlateStressArray, ComplexPlateStrainArray,
    #ComplexPlateVMStressArray, ComplexPlateVMStrainArray,

    RealCompositePlateStressArray, RealCompositePlateStrainArray,


    RealSolidStressArray, RealSolidStrainArray,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.utils import get_group_name, get_attributes
#from pyNastran2.op2.op2_interface.h5_pytables.utils import get_analysis_code_times
from .utils import break_domain_by_case, get_analysis_code_times # get_name
if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd
    import h5py
    from pyNastran.dev.op2_vectorized3.op2_geom import OP2
    GroupData = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]

def read_elemental_force(model: OP2, domains: np.ndarray,
                         elemental: Node, elemental_index: Node):
    return read_elemental_stress_strain(model, domains,
                                        elemental, elemental_index, 'force')

def read_elemental_stress(model: OP2, domains: np.ndarray,
                          elemental: Node, elemental_index: Node):
    return read_elemental_stress_strain(model, domains,
                                        elemental, elemental_index, 'stress')

def read_elemental_strain(model: OP2, domains: np.ndarray,
                          elemental: Node, elemental_index: Node):
    return read_elemental_stress_strain(model, domains,
                                        elemental, elemental_index, 'strain')


def read_elemental_stress_strain(model: OP2, domains: np.ndarray,
                                 elemental: Node, elemental_index: Node,
                                 stress_strain: str):
    """
    Paramters
    ---------
    stress_strain : str
       'stress', 'strain', 'force'
    """
    iresult = 0
    results = {}
    #geom_model = None
    ids = None

    log = model.log
    for h5_node_ in elemental._f_iter_nodes():
        cresult_name = ''
        if isinstance(h5_node_, Group):
            print(h5_node_)
            name = get_group_name(h5_node_)
            if name == 'TEMP_JUNK':
                temp
            else:
                raise NotImplementedError(name)

        elif isinstance(h5_node_, Node):
            attrs = get_attributes(h5_node_)
            assert len(attrs) == 1, attrs
            version = attrs['version'][0]

            #data = h5_node_.read()
            name = h5_node_.name
            group = elemental[name].read()
            index = elemental_index[name].read()
            #print(h5_node_)
            if name == 'BAR':
                cresult_name = f'{stress_strain}.cbar_{stress_strain}'
                iresult = _bar_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'BEAM':
                #log.warning(f'skipping {name} {stress_strain}')
                cresult_name = f'{stress_strain}.cbeam_{stress_strain}'
                iresult = _beam_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'ROD':
                cresult_name = f'{stress_strain}.crod_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TUBE':
                cresult_name = f'{stress_strain}.ctube_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'CONROD':
                cresult_name = f'{stress_strain}.conrod_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name in {'DAMP1', 'DAMP2', 'DAMP3', 'DAMP4'}:
                n = name[4]
                cresult_name = f'{stress_strain}.cdamp{n}_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name in {'ELAS1', 'ELAS2', 'ELAS3', 'ELAS4'}:
                n = name[4]
                cresult_name = f'{stress_strain}.celas{n}_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'VISC':
                cresult_name = f'{stress_strain}.cvisc_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'SHEAR':
                cresult_name = f'{stress_strain}.cshear_{stress_strain}'
                iresult = _real_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'QUAD4':
                cresult_name = f'{stress_strain}.cquad4_{stress_strain}'
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUAD8':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquad8_{stress_strain}'
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUADR':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquadr_{stress_strain}'
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIA3':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctria3_{stress_strain}'
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIA6':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctria6_{stress_strain}'
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIAR':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctriar_{stress_strain}'
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'QUAD4_CN':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquad4_{stress_strain}'
                iresult = _real_shell_element(
                    cresult_name, name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'QUAD4_COMP':
                cresult_name = f'{stress_strain}.cquad4_composite_{stress_strain}'
                iresult = _real_composite_shell(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIA3_COMP':
                cresult_name = f'{stress_strain}.ctria3_composite_{stress_strain}'
                iresult = _real_composite_shell(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUADR_COMP':
                cresult_name = f'{stress_strain}.cquadr_composite_{stress_strain}'
                iresult = _real_composite_shell(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIAR_COMP':
                cresult_name = f'{stress_strain}.ctriar_composite_{stress_strain}'
                iresult = _real_composite_shell(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'TETRA':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctetra_{stress_strain}'
                iresult = _real_solid_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'PENTA':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cpenta_{stress_strain}'
                iresult = _real_solid_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'HEXA':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.chexa_{stress_strain}'
                iresult = _real_solid_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            # -----------------------------------------
            #complex       01234
            elif name in {'ELAS1_CPLX', 'ELAS2_CPLX', 'ELAS3_CPLX', 'ELAS4_CPLX'}:
                n = name[4]
                cresult_name = f'{stress_strain}.celas{n}_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name in {'DAMP1_CPLX', 'DAMP2_CPLX', 'DAMP3_CPLX', 'DAMP4_CPLX'}:
                n = name[4]
                cresult_name = f'{stress_strain}.cdamp{n}_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'BAR_CPLX':
                pass
            elif name == 'BEAM_CPLX':
                pass
            elif name == 'CONROD_CPLX':
                cresult_name = f'{stress_strain}.conrod_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'ROD_CPLX':
                cresult_name = f'{stress_strain}.crod_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TUBE_CPLX':
                cresult_name = f'{stress_strain}.ctube_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'SHEAR_CPLX':
                cresult_name = f'{stress_strain}.cshear_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'BUSH_CPLX':
                cresult_name = f'{stress_strain}.cbush_{stress_strain}'
                iresult = load_complex_element(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'TRIA3_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctria3_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIA6_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctria6_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'TRIAR_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.ctriar_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUAD4_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquad4_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUAD8_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquad8_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUADR_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquadr_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'QUAD_CN':
                #if stress_strain == 'force':
                    #log.warning(f'skipping {name} {stress_strain}')
                    #return
                log.warning(f'skipping {name} {stress_strain}')
                #cresult_name = f'cquad_{stress_strain}'
                #iresult = _real_shell_element(
                    #cresult_name, iresult, results, domains, group, index, version,
                    #ids, model, subcases=None)
            elif name == 'QUAD4_CN_CPLX':
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return
                cresult_name = f'{stress_strain}.cquad4_{stress_strain}'
                iresult = _complex_shell_stress(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)
            elif name == 'QUAD_CN_CPLX':
                log.warning(f'skipping {name} {stress_strain}')

            elif name in {'TRIA3_COMP_CPLX', 'TRIAR_COMP_CPLX',
                          'QUAD4_COMP_CPLX', 'QUADR_COMP_CPLX'}:
                if stress_strain == 'force':
                    log.warning(f'skipping {name} {stress_strain}')
                    return

                if name == 'TRIA3_COMP_CPLX':
                    cresult_name = f'{stress_strain}.ctria3_composite_{stress_strain}'
                elif name == 'TRIAR_COMP_CPLX':
                    cresult_name = f'{stress_strain}.ctriar_composite_{stress_strain}'
                elif name == 'QUAD4_COMP_CPLX':
                    cresult_name = f'{stress_strain}.cquad4_composite_{stress_strain}'
                elif name == 'QUADR_COMP_CPLX':
                    cresult_name = f'{stress_strain}.cquadr_composite_{stress_strain}'

                iresult = _complex_composite_shell(
                    cresult_name, iresult, results, domains, group, index, version,
                    ids, model, subcases=None)

            elif name == 'HEXA_CPLX':
                log.warning(f'skipping {name} {stress_strain}')
            elif name == 'PENTA_CPLX':
                log.warning(f'skipping {name} {stress_strain}')
            elif name == 'TETRA_CPLX':
                log.warning(f'skipping {name} {stress_strain}')

            elif name == 'GRAD_FLUX':
                log.warning(f'skipping {name} {stress_strain}')
            elif name == 'HBDYE':
                log.warning(f'skipping {name} {stress_strain}')
            else:
                print(h5_node_)
                #raise NotImplementedError(name)
        else:
            print(h5_node_)
            raise NotImplementedError(h5_node_)

def load_complex_element(result_name: str,
                         iresult: int,
                         results: dict[int, Function],
                         domains_df: pd.DataFrame,
                         group: h5py._hl.dataset.Dataset,
                         index: h5py._hl.dataset.Dataset,
                         version: int,
                         ids: np.ndarray,
                         model: OP2,
                         subcases: Optional[list[int]]=None) -> int:

    if result_name in {'stress.crod_stress', 'stress.conrod_stress', 'stress.ctube_stress'}:
        class_obj = ComplexRodStressArray
        table_name = 'OES1'
    elif result_name  in {'force.crod_force', 'force.conrod_force', 'force.ctube_force'}:
        class_obj = ComplexRodForceArray
        table_name = 'OEF1'
    elif result_name in {'strain.crod_strain', 'strain.conrod_strain', 'strain.ctube_strain'}:
        class_obj = ComplexRodStrainArray
        table_name = 'OSTR1'

    elif result_name in {'stress.celas1_stress', 'stress.celas2_stress', 'stress.celas3_stress'}:
        class_obj = ComplexSpringStressArray
        table_name = 'OES1'
    elif result_name in {'strain.celas1_strain', 'strain.celas2_strain', 'strain.celas3_strain'}:
        class_obj = ComplexSpringStrainArray
        table_name = 'OSTR1'
    elif result_name in {'force.celas1_force', 'force.celas2_force', 'force.celas3_force', 'force.celas4_force'}:
        class_obj = ComplexSpringForceArray
        table_name = 'OEF1'
    elif result_name in {'force.cdamp1_force', 'force.cdamp2_force', 'force.cdamp3_force', 'force.cdamp4_force'}:
        class_obj = ComplexDamperForceArray
        table_name = 'OEF1'

    elif result_name == 'stress.cshear_stress':
        class_obj = ComplexShearStressArray
        table_name = 'OES1'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)
    #is_strain = ('STRAIN' == stress_strain)
    is_force = ('FORCE' == stress_strain)

    basename = result_name + ' (complex)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'

    # dtype([('EID', '<i8'), ('AR', '<f8'), ('AI', '<f8'), ('TR', '<f8'), ('TI', '<f8'), ('DOMAIN_ID', '<i8')])  # CONROD/CROD
    #dtype([('EID', '<i8'), ('ASR', '<f8'), ('ASI', '<f8'), ('TSR', '<f8'), ('TSI', '<f8'), ('DOMAIN_ID', '<i8')])  # CTUBE
    EID = group['EID']
    nelements = len(EID)
    #names = group.dtype.names

    if is_force:
        if result_name == 'force.conrod_force':
            #('EID', 'AFR', 'AFI', 'TRQR', 'TRQI', 'DOMAIN_ID')
            assert version == 0, version
            nresults = 2
            A = group['AFR'] + group['AFI'] * 1j
            T = group['TRQR'] + group['TRQI'] * 1j
            DATA = hstack_shape(nelements, A, T)
        elif 'celas' in result_name or 'cdamp' in result_name:
            #('EID', 'FR', 'FI', 'DOMAIN_ID')
            assert version == 0, version
            nresults = 1
            F = group['FR'] + group['FI'] * 1j
            DATA = hstack_shape(nelements, F)
        else:  # pragma: no cover
            raise NotImplementedError(result_name)
    elif 'crod_' in result_name or 'conrod_' in result_name:
        assert version == 0, version
        nresults = 2
        A = group['AR'] + group['AI'] * 1j
        T = group['TR'] + group['TI'] * 1j
        DATA = hstack_shape(nelements, A, T)
    elif 'ctube_' in result_name:
        assert version == 0, version
        nresults = 2
        A = group['ASR'] + group['ASI'] * 1j
        T = group['TSR'] + group['TSI'] * 1j
        DATA = hstack_shape(nelements, A, T)
    elif 'cshear_' in result_name:
        assert version == 0, version
        nresults = 2
        tmax = group['TMAXR'] + group['TMAXI'] * 1j
        tavg = group['TAVGR'] + group['TAVGI'] * 1j
        DATA = hstack_shape(nelements, tmax, tavg)
    elif 'celas' in result_name:
        assert version == 0, version
        nresults = 1
        DATA = (group['SR'] + group['SI'] * 1j).reshape(nelements, 1)
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    mode, eigr, eigi, data_by_group = _get_data_by_group_element(
        basename, domains_df, index,
        EID, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            ntimes = len(group)
            for itime, (idomaini, elementi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(elementi)
                nelement.append(nelementi)
                datas.append(datai)
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nelements = len(elementi)
            if data.shape != (ntimes, nelements, nresults):
                raise RuntimeError(f'result_name={result_name}; nresults={nresults} data.shape={data.shape}')

            if is_freq:
                freq = freq_data
                assert isinstance(element_name, str), element_name
                if not hasattr(class_obj, 'add_freq_case'):
                    model.log.warning(f'skipping {result_name} - add_freq_case')
                    continue

                res = class_obj.add_freq_case(
                    table_name,
                    elementi, data, isubcase,
                    freq,
                    element_name,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')

            elif is_complex_modes:
                if not hasattr(class_obj, 'add_complex_modes_case'):
                    model.log.warning(f'skipping {result_name} - add_complex_modes_case')
                    continue
                modes, eigr, eigi = complex_modes
                res = class_obj.add_complex_modes_case(
                    table_name, elementi, data, isubcase,
                    modes, eigr, eigi,
                    element_name,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = res
            #model.displacements[isubcase] = disp

    return iresult

def hstack_shape(nelements: int, *arrays):
    DATA = np.hstack([array.reshape(nelements, 1) for array in arrays])
    return DATA

#result_name_to_class_table_name_map = {}
def _real_element(result_name: str,
                  iresult: int,
                  results: dict[int, Function],
                  domains_df: pd.DataFrame,
                  group: h5py._hl.dataset.Dataset,
                  index: h5py._hl.dataset.Dataset,
                  version: int,
                  ids: np.ndarray,
                  model: OP2,
                  subcases: Optional[list[int]]=None) -> int:

    if result_name in {'stress.crod_stress', 'stress.conrod_stress', 'stress.ctube_stress'}:
        class_obj = RealRodStressArray
        table_name = 'OES1'
    elif result_name  in {'force.crod_force', 'force.conrod_force', 'force.ctube_force'}:
        class_obj = RealRodForceArray
        table_name = 'OEF1'
    elif result_name in {'strain.crod_strain', 'strain.conrod_strain', 'strain.ctube_strain'}:
        class_obj = RealRodStrainArray
        table_name = 'OSTR1'

    elif result_name in {'stress.celas1_stress', 'stress.celas2_stress', 'stress.celas3_stress'}:
        class_obj = RealSpringStressArray
        table_name = 'OES1'
    elif result_name in {'strain.celas1_strain', 'strain.celas2_strain', 'strain.celas3_strain'}:
        class_obj = RealSpringStrainArray
        table_name = 'OSTR1'
    elif result_name in {'force.celas1_force', 'force.celas2_force',
                         'force.celas3_force', 'force.celas4_force'}:
        class_obj = RealSpringForceArray
        table_name = 'OEF1'
    elif result_name in {'force.cdamp1_force', 'force.cdamp2_force',
                         'force.cdamp3_force', 'force.cdamp4_force'}:
        class_obj = RealDamperForceArray
        table_name = 'OEF1'
    elif result_name == 'force.cvisc_force':
        class_obj = RealViscForceArray
        table_name = 'OEF1'

    elif result_name == 'stress.cshear_stress':
        class_obj = RealShearStressArray
        table_name = 'OES1'
    #elif result_name == 'cshear_strain':
        #class_obj = RealShearStrainArray
        #table_name = 'OSTR1'
    elif result_name == 'force.cshear_force':
        class_obj = RealCShearForceArray
        table_name = 'OEF1'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)
    #is_strain = ('STRAIN' == stress_strain)
    #is_force = ('FORCE' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    EID = group['EID']
    nelements = len(EID)
    nlength = index['LENGTH'].sum()

    if result_name in {'force.conrod_force', 'force.crod_force', 'force.ctube_force', 'force.cvisc_force'}:
        #('EID', 'AF', 'TRQ', 'DOMAIN_ID')
        assert version == 0, version
        axial = group['AF']
        torque = group['TRQ']
        DATA = hstack_shape(nelements, axial, torque)
    elif result_name in {'force.cdamp1_force', 'force.cdamp2_force', 'force.cdamp3_force', 'force.cdamp4_force',
                         'force.celas1_force', 'force.celas2_force', 'force.celas3_force', 'force.celas4_force',}:
        #('EID', 'F', 'DOMAIN_ID')
        assert version == 0, version
        force = group['F']
        DATA = hstack_shape(nelements, force)
    elif result_name in {'strain.conrod_strain', 'strain.crod_strain',
                         'stress.conrod_stress', 'stress.crod_stress'}:
        #('EID', 'A', 'T', 'MSA', 'MST', 'DOMAIN_ID')
        assert version == 0, version
        axial = group['A']
        torque = group['T']
        axial_margin = group['MSA']
        torsion_margin = group['MST']
        DATA = hstack_shape(nelements, axial, torque, axial_margin, torsion_margin)
    elif result_name in {'stress.ctube_stress', 'strain.ctube_strain'}:
        #('EID', 'AS', 'MSA', 'TS', 'MST', 'DOMAIN_ID')
        assert version == 0, version
        axial = group['AS']
        torque = group['TS']
        axial_margin = group['MSA']
        torsion_margin = group['MST']
        DATA = hstack_shape(nelements, axial, torque, axial_margin, torsion_margin)

    elif result_name in {'stress.celas1_stress', 'stress.celas2_stress', 'stress.celas3_stress', 'stress.celas4_stress',
                         'strain.celas1_strain', 'strain.celas2_strain', 'strain.celas3_strain', 'strain.celas4_strain'}:
        #('EID', 'S', 'DOMAIN_ID')
        assert version == 0, version
        stress = group['S']
        DATA = hstack_shape(nelements, stress)
    elif result_name in {'stress.cshear_stress', 'strain.cshear_strain'}:
        # dtype=[('EID', '<i8'), ('TMAX', '<f8'), ('TAVG', '<f8'), ('MS', '<f8'), ('DOMAIN_ID', '<i8')])
        assert version == 0, version
        tmax = group['TMAX']
        tavg = group['TAVG']
        margin = group['MS']
        DATA = hstack_shape(nelements, tmax, tavg, margin)

    elif result_name == 'force.cshear_force':
        assert version == 0, version
        # dtype=[('EID', '<i8'), ('F41', '<f8'), ('F21', '<f8'), ('F12', '<f8'), ('F32', '<f8'),
        #        ('F23', '<f8'), ('F43', '<f8'), ('F34', '<f8'), ('F14', '<f8'), ('KF1', '<f8'),
        #        ('S12', '<f8'), ('KF2', '<f8'), ('S23', '<f8'), ('KF3', '<f8'), ('S34', '<f8'),
        #        ('KF4', '<f8'), ('S41', '<f8'), ('DOMAIN_ID', '<i8')])
        #'force41', 'force21', 'force12', 'force32', 'force23', 'force43',
        #'force34', 'force14',
        #'kick_force1', 'shear12', 'kick_force2', 'shear23',
        #'kick_force3', 'shear34', 'kick_force4', 'shear41',
        force41 = group['F41']
        force21 = group['F21']
        force12 = group['F12']
        force32 = group['F32']
        force23 = group['F23']
        force43 = group['F43']
        force34 = group['F34']
        force14 = group['F14']
        kick_force1 = group['KF1']
        kick_force2 = group['KF2']
        kick_force3 = group['KF3']
        kick_force4 = group['KF4']
        shear12 = group['S12']
        shear23 = group['S23']
        shear34 = group['S34']
        shear41 = group['S41']
        DATA = hstack_shape(
            nelements,
            force41, force21, force12, force32,
            force23, force43, force34, force14,
            kick_force1, shear12, kick_force2, shear23,
            kick_force3, shear34, kick_force4, shear41)
    else:  # pragma: no cover
        raise NotImplementedError(result_name)

    assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'
    mode, eigr, eigi, data_by_group = _get_data_by_group_element(
        basename, domains_df, index,
        EID, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            #ntimes = len(group)
            for itime, (idomaini, elementi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(elementi)
                nelement.append(nelementi)
                datas.append(datai)
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nelements = len(elementi)
            #assert data.shape == (ntimes, nelements, 9), data.shape

            if is_static:
                disp = class_obj.add_static_case(
                    table_name, element_name, elementi, data, isubcase,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                #if not hasattr(class_obj, 'add_transient_case'):
                    #model.log.warning(f'skipping {result_name} - add_transient_case')
                    #continue
                disp = class_obj.add_transient_case(
                    table_name, element_name, elementi, data, isubcase,
                    times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_modes:
                modes, eigenvalues, freqs = modes_data
                #if not hasattr(class_obj, 'add_modal_case'):
                    #model.log.warning(f'skipping {result_name} - add_modal_case')
                    #continue
                disp = class_obj.add_modal_case(
                    table_name, element_name, elementi, data, isubcase,
                    modes, eigenvalues, freqs,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                if not hasattr(class_obj, 'add_post_buckling_case'):
                    model.log.warning(f'skipping {result_name} - add_post_buckling_case')
                    continue
                modes, eigrs, eigis = buckling_data
                disp = class_obj.add_post_buckling_case(
                    table_name, element_name, elementi, data, isubcase,
                    modes, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _bar_element(result_name: str,
                  iresult: int,
                  results: dict[int, Function],
                  domains_df: pd.DataFrame,
                  group: h5py._hl.dataset.Dataset,
                  index: h5py._hl.dataset.Dataset,
                  version: int,
                  ids: np.ndarray,
                  model: OP2,
                  subcases: Optional[list[int]]=None) -> int:

    if result_name == 'stress.cbar_stress':
        class_obj = RealBarStressArray
        table_name = 'OES1X'
    elif result_name == 'strain.cbar_strain':
        class_obj = RealBarStrainArray
        table_name = 'OSTR1X'
    elif result_name == 'force.cbar_force':
        class_obj = RealCBarForceArray
        table_name = 'OEF1X'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)

    #is_strain = ('STRAIN' == stress_strain)
    #is_force = ('FORCE' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    EID = group['EID']
    nelements = len(EID)
    nlength = index['LENGTH'].sum()

    if result_name == 'force.cbar_force':
        # ('EID', 'BM1A', 'BM2A', 'BM1B', 'BM2B', 'TS1', 'TS2', 'AF', 'TRQ', 'DOMAIN_ID')
        assert version == 0, version
        bm1a = group['BM1A']
        bm2a = group['BM2A']
        bm1b = group['BM1B']
        bm2b = group['BM2B']
        ts1  = group['TS1']
        ts2  = group['TS2']
        axial = group['AF']
        torque = group['TRQ']
        DATA = hstack_shape(nelements, bm1a, bm2a, bm1b, bm2b, ts1, ts2, axial, torque)
    elif result_name in {'stress.cbar_stress', 'strain.cbar_strain'}:
        # ('EID', 'X1A', 'X2A', 'X3A', 'X4A', 'MAXA', 'MINA', 'MST', 'AX',
        #        'X1B', 'X2B', 'X3B', 'X4B', 'MAXB', 'MINB', 'MSC', 'DOMAIN_ID')
        #headers = ['s1a', 's2a', 's3a', 's4a', 'axial', 'smaxa', 'smina', 'MS_tension',
                   #'s1b', 's2b', 's3b', 's4b',          'smaxb', 'sminb', 'MS_compression']
        s1a = group['X1A']
        s2a = group['X2A']
        s3a = group['X3A']
        s4a = group['X4A']

        s1b = group['X1B']
        s2b = group['X2B']
        s3b = group['X3B']
        s4b = group['X4B']
        smaxa = group['MAXA']
        smina = group['MINA']

        smaxb = group['MAXB']
        sminb = group['MINB']
        mst = group['MST']
        msc = group['MSC']
        axial = group['AX']
        DATA = hstack_shape(nelements,
                            s1a, s2a, s3a, s4a, axial, smaxa, smina, mst,
                            s1b, s2b, s3b, s4b,        smaxb, sminb, msc)
    else:  # pragma: no cover
        raise NotImplementedError(result_name)

    assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'
    mode, eigr, eigi, data_by_group = _get_data_by_group_element(
        basename, domains_df, index,
        EID, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            #ntimes = len(group)
            for itime, (idomaini, elementi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(elementi)
                nelement.append(nelementi)
                datas.append(datai)
                #x = 1
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nelements = len(elementi)
            #assert data.shape == (ntimes, nelements, 9), data.shape

            if is_static:
                disp = class_obj.add_static_case(
                    table_name, element_name, elementi, data, isubcase,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                #if not hasattr(class_obj, 'add_transient_case'):
                    #model.log.warning(f'skipping {result_name} - add_transient_case')
                    #continue
                disp = class_obj.add_transient_case(
                    table_name, element_name, elementi, data, isubcase,
                    times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_modes:
                modes, eigenvalues, freqs = modes_data
                #if not hasattr(class_obj, 'add_modal_case'):
                    #model.log.warning(f'skipping {result_name} - add_modal_case')
                    #continue
                disp = class_obj.add_modal_case(
                    table_name, element_name, elementi, data, isubcase,
                    modes, eigenvalues, freqs,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                #if not hasattr(class_obj, 'add_post_buckling_case'):
                    #model.log.warning(f'skipping {result_name} - add_post_buckling_case')
                    #continue
                modes, eigrs, eigis = buckling_data
                disp = class_obj.add_post_buckling_case(
                    table_name, element_name, elementi, data, isubcase,
                    modes, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _beam_element(result_name: str,
                  iresult: int,
                  results: dict[int, Function],
                  domains_df: pd.DataFrame,
                  group: h5py._hl.dataset.Dataset,
                  index: h5py._hl.dataset.Dataset,
                  version: int,
                  ids: np.ndarray,
                  model: OP2,
                  subcases: Optional[list[int]]=None) -> int:

    if result_name == 'stress.cbeam_stress':
        class_obj = RealBeamStressArray
        table_name = 'OES1X'
    elif result_name == 'strain.cbeam_strain':
        class_obj = RealBeamStrainArray
        table_name = 'OSTR1X'
    elif result_name == 'force.cbeam_force':
        class_obj = RealCBeamForceArray
        table_name = 'OEF1X'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)

    #is_strain = ('STRAIN' == stress_strain)
    #is_force = ('FORCE' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    EID = group['EID']
    nelements = len(EID)
    nlength = index['LENGTH'].sum()

    nstations_per_element = 11
    nstations = nelements * nstations_per_element
    if result_name == 'force.cbeam_force':
        #('EID', 'GRID', 'SD', 'BM1', 'BM2', 'TS1', 'TS2', 'AF', 'TTRQ', 'WTRQ', 'DOMAIN_ID')
        assert version == 0, version
        #'sd', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2',
        #'axial_force', 'total_torque', 'warping_torque']
        element_nodes = np.zeros((nstations, 2), dtype=EID.dtype)
        element_nodes[:, 0] = np.repeat(EID, nstations_per_element)
        element_nodes[:, 1] = group['GRID'].ravel()
        xxb = group['SD']
        bm1 = group['BM1']
        bm2 = group['BM2']
        shear1  = group['TS1']
        shear2  = group['TS2']
        axial = group['AF']
        torque = group['TTRQ']
        warping_torque = group['WTRQ']
        nresults = 8
        DATA = np.stack([xxb, bm1, bm2, shear1, shear2,
                         axial, torque, warping_torque], axis=2)
    elif result_name in {'stress.cbeam_stress', 'strain.cbeam_strain'}:
        assert version == 0, version
        element_nodes = np.zeros((nstations, 2), dtype=EID.dtype)
        element_nodes[:, 0] = np.repeat(EID, nstations_per_element)
        element_nodes[:, 1] = group['GRID'].ravel()

        xxb = group['SD']
        #('EID', 'GRID', 'SD', 'XC', 'XD', 'XE', 'XF', 'MAX', 'MIN', 'MST', 'MSC', 'DOMAIN_ID')
        ##'grid', 'xxb',
        #'sxc', 'sxd', 'sxe', 'sxf',
        #'smax', 'smin', 'MS_tension', 'MS_compression'

        sxc = group['XC']
        sxd = group['XD']
        sxe = group['XE']
        sxf = group['XF']

        smax = group['MAX']
        smin = group['MIN']

        mst = group['MST']
        msc = group['MSC']
        nresults = 8
        DATA = np.stack([sxc, sxd, sxe, sxf,
                         smax, smin, mst, msc], axis=2)
    else:  # pragma: no cover
        raise NotImplementedError(result_name)
    element_nodes = element_nodes.reshape(nelements, nstations_per_element, 2)

    assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'

    #DATA = DATA.reshape(nelements, nstations_per_element, nresults)
    mode, eigr, eigi, data_by_group = _get_data_by_group_fiber_element(
        basename, domains_df, index,
        EID, element_nodes, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            #ntimes = len(group)
            for itime, (idomaini, elementi, element_nodei, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(elementi)
                nelement.append(nelementi)
                nstationsi = nstations_per_element * nelementi
                datas.append(datai.reshape(nstationsi, nresults))
            element_node = element_nodei[0, :, :]

            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nelements = len(elementi)
            #assert data.shape == (ntimes, nelements, 9), data.shape

            if is_static:
                disp = class_obj.add_static_case(
                    table_name, element_name, element_node, xxb, data, isubcase,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                #if not hasattr(class_obj, 'add_transient_case'):
                    #model.log.warning(f'skipping {result_name} - add_transient_case')
                    #continue
                disp = class_obj.add_transient_case(
                    table_name, element_name, element_node, xxb, data, isubcase,
                    times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_modes:
                modes, eigenvalues, freqs = modes_data
                #if not hasattr(class_obj, 'add_modal_case'):
                    #model.log.warning(f'skipping {result_name} - add_modal_case')
                    #continue
                disp = class_obj.add_modal_case(
                    table_name, element_name, element_node, xxb, data, isubcase,
                    modes, eigenvalues, freqs,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                #if not hasattr(class_obj, 'add_post_buckling_case'):
                    #model.log.warning(f'skipping {result_name} - add_post_buckling_case')
                    #continue
                #model.log.warning(f'skipping {result_name} - add_buckling_case')
                #continue
                modes, eigrs, eigis = buckling_data
                disp = class_obj.add_modal_case(
                    table_name, element_name, element_node, xxb, data, isubcase,
                    modes, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _real_shell_element(result_name: str,
                        group_name: str,
                        iresult: int,
                        results: dict[int, Function],
                        domains_df: pd.DataFrame,
                        group: h5py._hl.dataset.Dataset,
                        index: h5py._hl.dataset.Dataset,
                        version: int,
                        ids: np.ndarray,
                        model: OP2,
                        subcases: Optional[list[int]]=None) -> int:

    if result_name in {'stress.ctria3_stress', 'stress.ctria6_stress', 'stress.ctriar_stress',
                       'stress.cquad4_stress', 'stress.cquad8_stress', 'stress.cquadr_stress', 'stress.cquad_stress'}:
        class_obj = RealPlateStressArray
        table_name = 'OES1'
    elif result_name in {'strain.ctria3_strain', 'strain.ctria6_strain', 'strain.ctriar_strain',
                         'strain.cquad4_strain', 'strain.cquad8_strain', 'strain.cquadr_strain', 'strain.cquad_strain'}:
        class_obj = RealPlateStrainArray
        table_name = 'OSTR1'
    else:
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)
    #is_strain = ('STRAIN' == stress_strain)
    #is_force = ('FORCE' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    EID = group['EID']
    nelements = len(EID)
    nlength = index['LENGTH'].sum()

    if 'ctria3' in result_name:
        nnodes = 1
        # ('EID', 'FD1', 'X1', 'Y1', 'TXY1', 'FD2', 'X2', 'Y2', 'TXY2', 'DOMAIN_ID')
        assert version == 0, version
        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]

        nlayers_half = nelements
        nlayers = nlayers_half * 2

        element_node = np.zeros((nlayers_half, 4), dtype=EID.dtype)
        element_node[:, 0] = EID
        element_node[:, 2] = EID

        FIBER = np.stack([group['FD1'], group['FD2']], axis=1)
        assert FIBER.shape[0] == nlength, f'fiber.shape={FIBER.shape}; nlength={nlength}'
        X1 = group['X1']
        Y1 = group['Y1']
        TXY1 = group['TXY1']
        X2 = group['X2']
        Y2 = group['Y2']
        TXY2 = group['TXY2']
        #fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm
        angle = np.nan * X1
        major = np.nan * X1
        minor = np.nan * X1
        ovm = np.nan * X1

        DATA = np.stack([
            group['FD1'], X1, Y1, TXY1, angle, major, minor, ovm,
            group['FD2'], X2, Y2, TXY2, angle, major, minor, ovm,
        ], axis=1)
        assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'
    elif group_name == 'QUAD4':
        # why is there no T????
        #('EID', 'FD1', 'X1', 'Y1', 'XY1', 'FD2', 'X2', 'Y2', 'XY2', 'DOMAIN_ID')
        assert version == 0, version
        nnodes = 1
        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]

        nlayers_half = nelements
        #nlayers = nlayers_half * 2
        element_node = np.zeros((nlayers_half, 4), dtype=EID.dtype)
        element_node[:, 0] = EID
        element_node[:, 2] = EID

        FIBER = np.stack([group['FD1'], group['FD2']], axis=1)
        assert FIBER.shape[0] == nlength, f'fiber.shape={FIBER.shape}; nlength={nlength}'
        X1 = group['X1']
        Y1 = group['Y1']
        TXY1 = group['XY1']
        X2 = group['X2']
        Y2 = group['Y2']
        TXY2 = group['XY2']
        #fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm
        angle = np.nan * X1
        major = np.nan * X1
        minor = np.nan * X1
        ovm = np.nan * X1

        DATA = np.stack([
            group['FD1'], X1, Y1, TXY1, angle, major, minor, ovm,
            group['FD2'], X2, Y2, TXY2, angle, major, minor, ovm,
        ], axis=1)
        assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'
    else:# if 'cquad8' in result_name:
        # ('EID', 'FD1', 'X1', 'Y1', 'TXY1', 'FD2', 'X2', 'Y2', 'TXY2', 'DOMAIN_ID')
        assert version == 0, version
        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]

        X1 = group['X1']
        Y1 = group['Y1']
        TXY1 = group['TXY1']
        X2 = group['X2']
        Y2 = group['Y2']
        TXY2 = group['TXY2']

        nnodes = X1.shape[1]

        nlayers_half = X1.size
        nlayers = nlayers_half * 2
        element_node = np.zeros((nlayers_half, 2), dtype=EID.dtype)

        grid = group['GRID']
        grid[:, 0] = 0  #  why was the centroid not 0???
        element_node[:, 0] = np.repeat(EID, nnodes)
        element_node[:, 1] = grid.ravel()
        element_node = element_node.reshape((nelements, 2*nnodes))
        #element_node = np.zeros((nelements, 2*nnodes), dtype=EID.dtype)

        #irange = np.arange(0, 2*nnodes, 2)
        #element_node[:, irange] = EID[:, np.newaxis]
        #for ii in irange:

        FIBER = np.stack([
            group['FD1'], # .reshape(nlayers_half, 1),
            group['FD2'], # .reshape(nlayers_half, 1),
        ], axis=2) # .reshape(nlayers)
        assert FIBER.shape[0] == nlength, f'fiber.shape={FIBER.shape}; nlength={nlength}'
        assert len(FIBER) == nelements
        #fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm
        angle = np.nan * X1
        major = np.nan * X1
        minor = np.nan * X1
        ovm = np.nan * X1
        DATA = np.stack([
            group['FD1'], X1, Y1, TXY1, angle, major, minor, ovm,
            group['FD2'], X2, Y2, TXY2, angle, major, minor, ovm,
        ], axis=2)
        assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'
    #else:  # pragma: no cover
        #raise NotImplementedError(result_name)

    if nnodes == 5 and element_name == 'CQUAD4':
        element_name = 'CQUAD4-144'
    elif nnodes == 1 and element_name in {'CQUAD8', 'CQUADR', 'CTRIA6', 'CTRIAR'}:
        element_name += '_LINEAR'
        raise RuntimeError(element_name)
    elif element_name in {'CQUAD8', 'CTRIA6'}:
        model.log.warning(f'skipping shell {element_name}')
        return iresult

    #assert element_node.shape[1] == 4, element_node.shape
    mode, eigr, eigi, data_by_group = _get_data_by_group_fiber_element(
        basename, domains_df, index,
        element_node, FIBER, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            #ntimes = len(group)
            for itime, (idomaini, element_nodei, fiberi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(fiberi)
                nelement.append(nelementi)
                datas.append(datai)
            fiber_distance = fiberi.flatten()
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)

            # reshape it from:
            #  (N, 4) to (2N, 2); nlayers=2N
            #  (N, 4*5)=(N, 10*2) to (10N, 2); nlayers=10N
            #
            nlayersi = element_nodei.shape[0] * element_nodei.shape[1] // 2
            element_nodei = element_nodei.reshape((nlayersi, 2))

            if data.ndim == 3:
                # ctria3
                ntimes, _nelements, nresults2 = data.shape
                data = data.reshape(ntimes, 2*nelementi, 8)
            else:
                ntimes, nelements, nnodes, nresults2 = data.shape
                data = data.reshape(ntimes, 2*nelementi*nnodes, 8)
            assert nresults2 == 16

            #nelements = len(elementi)
            #assert data.shape == (ntimes, nelements, 9), data.shape

            if is_static:
                disp = class_obj.add_static_case(
                    table_name, element_name, nnodes, element_nodei, fiber_distance, data, isubcase,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                if not hasattr(class_obj, 'add_transient_case'):
                    model.log.warning(f'skipping {result_name} - add_transient_case')
                    continue
                disp = class_obj.add_transient_case(
                    table_name, element_name, nnodes, element_nodei, fiber_distance, data, isubcase,
                    times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_modes:
                modes, eigenvalues, freqs = modes_data
                if not hasattr(class_obj, 'add_modal_case'):
                    model.log.warning(f'skipping {result_name} - add_modal_case')
                    continue
                disp = class_obj.add_modal_case(
                    table_name, element_name, nnodes, element_nodei, fiber_distance, data, isubcase,
                    modes, eigenvalues, freqs,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                model.log.warning(f'skipping {result_name} - add_buckling_case')
                continue
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _real_solid_element(result_name: str,
                        iresult: int,
                        results: dict[int, Function],
                        domains_df: pd.DataFrame,
                        group: h5py._hl.dataset.Dataset,
                        index: h5py._hl.dataset.Dataset,
                        version: int,
                        ids: np.ndarray,
                        model: OP2,
                        subcases: Optional[list[int]]=None) -> int:

    if result_name in {'stress.chexa_stress', 'stress.cpenta_stress', 'stress.ctetra_stress'}:
        class_obj = RealSolidStressArray
        table_name = 'OES1'
    elif result_name in {'strain.chexa_strain', 'strain.cpenta_strain', 'strain.ctetra_strain'}:
        class_obj = RealSolidStrainArray
        table_name = 'OSTR1'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)
    #is_strain = ('STRAIN' == stress_strain)
    #is_force = ('FORCE' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    EID = group['EID']
    nelements = len(EID)

    if 'chexa_' in result_name or 1:
    #elif 'ctetra_' in result_name:
    #elif 'cpenta_' in result_name:
        # dtype({'names': ['EID', 'CID', 'CTYPE', 'NODEF', 'GRID', 'X', 'Y', 'Z',
        #                  'TXY', 'TYZ', 'TZX', 'DOMAIN_ID'],
        #        'formats': ['<i8', '<i8', 'S4', '<i8', ('<i8', (9,)), ('<f8', (9,)), ('<f8', (9,)),
        #                    ('<f8', (9,)), ('<f8', (9,)), ('<f8', (9,)), ('<f8', (9,)), '<i8'],
        #        'offsets': [0, 8, 16, 24, 32, 104, 176, 248, 320, 392, 464, 536], 'itemsize': 544})
        assert version == 0, version
        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        cid = group['CID']
        ctype = group['CTYPE']
        nodef = group['NODEF']
        grid = group['GRID'].ravel()
        ngrid = len(grid)

        oxx = group['X']
        oyy = group['Y']
        ozz = group['Z']
        txy = group['TXY']
        tyz = group['TYZ']
        txz = group['TZX']
        o1 = np.nan * oxx
        o2 = np.nan * oxx
        o3 = np.nan * oxx
        ovm_shear = np.nan * oxx

        nnodes_per_element = ngrid // nelements

        element_node = np.zeros((ngrid, 2), dtype=EID.dtype)
        element_node[:, 0] = np.repeat(EID, nnodes_per_element)
        element_node[:, 1] = grid

        element_cid = np.zeros((nelements, 2), dtype=EID.dtype)
        element_cid[:, 0] = EID
        element_cid[:, 1] = cid
        DATA = np.stack([oxx, oyy, ozz, txy, tyz, txz,
                         o1, o2, o3, ovm_shear], axis=2)
    else:  # pragma: no cover
        raise NotImplementedError(result_name)
    element_node = element_node.reshape(nelements, nnodes_per_element, 2)
    #element_cid = element_cid.reshape(nelements, 2)

    mode, eigr, eigi, data_by_group = _get_data_by_group_fiber_element(
        basename, domains_df, index,
        element_node, element_cid, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            ntimes = len(group)
            for itime, (idomaini, element_nodei, element_cidi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(element_cidi)
                nelement.append(nelementi)
                datas.append(datai)
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            #nelements = len(elementi)
            #assert data.shape == (ntimes, nelements, 9), data.shape
            element_nodei = element_nodei.reshape(nelementi*nnodes_per_element, 2)
            data = data.reshape(ntimes, nelementi*nnodes_per_element, 10)

            if is_static:
                disp = class_obj.add_static_case(
                    table_name, element_name, element_nodei, element_cidi, data, isubcase,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                if not hasattr(class_obj, 'add_transient_case'):
                    model.log.warning(f'skipping {result_name} - add_transient_case')
                    continue
                disp = class_obj.add_transient_case(
                    table_name, element_name, element_nodei, element_cidi, data, isubcase,
                    times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_modes:
                modes, eigenvalues, freqs = modes_data
                if not hasattr(class_obj, 'add_modal_case'):
                    model.log.warning(f'skipping {result_name} - add_modal_case')
                    continue
                disp = class_obj.add_modal_case(
                    table_name, element_name, element_nodei, element_cidi, data, isubcase,
                    modes, eigenvalues, freqs,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                modes, eigrs, eigis = buckling_data
                disp = class_obj.add_post_buckling_case(
                    table_name, element_name, element_nodei, element_cidi, data, isubcase,
                    modes, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
                #model.log.warning(f'skipping {result_name} - add_post_buckling_case')
                continue
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _real_composite_shell(result_name: str,
                          iresult: int,
                          results: dict[int, Function],
                          domains_df: pd.DataFrame,
                          group: h5py._hl.dataset.Dataset,
                          index: h5py._hl.dataset.Dataset,
                          version: int,
                          ids: np.ndarray,
                          model: OP2,
                          subcases: Optional[list[int]]=None) -> int:
    if result_name in {'stress.ctria3_composite_stress', 'stress.ctria6_composite_stress', 'stress.ctriar_composite_stress',
                       'stress.cquad4_composite_stress', 'stress.cquad8_composite_stress', 'stress.cquadr_composite_stress'}:
        class_obj = RealCompositePlateStressArray
        table_name = 'OES1C'
    elif result_name in {'strain.ctria3_composite_strain', 'strain.ctria6_composite_strain', 'strain.ctriar_composite_strain',
                         'strain.cquad4_composite_strain', 'strain.cquad8_composite_strain', 'strain.cquadr_composite_strain'}:
        class_obj = RealCompositePlateStrainArray
        table_name = 'OSTR1C'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, composite, stress_strain = result_name.upper().split('.')[1].split('_')
    #is_strain = ('STRAIN' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    assert version == 0, version
    EID = group['EID']
    nlayers = len(EID)

    element_layer = np.zeros((nlayers, 2), dtype=EID.dtype)
    element_layer[:, 0] = EID
    element_layer[:, 1] = group['PLY']
    X1 = group['X1']
    Y1 = group['Y1']
    T1 = group['T1']
    L1 = group['L1']
    L2 = group['L2']
    angle = X1 * np.nan
    major = X1 * np.nan
    minor = X1 * np.nan
    ovm = X1 * np.nan

    #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
    #if X1.ndim == 1:
        #nnodes_per_element = 1
    #else:
        #nnodes_per_element = X1.shape[1]

    DATA = hstack_shape(nlayers,
                        X1, Y1, T1, L1, L2,
                        angle, major, minor, ovm)
    mode, eigr, eigi, data_by_group = _get_data_by_group_element(
        basename, domains_df, index,
        element_layer, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            ntimes = len(group)
            for itime, (idomaini, elementi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(elementi)
                nelement.append(nelementi)
                datas.append(datai)
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nelements = len(elementi)
            assert data.shape == (ntimes, nelements, 9), data.shape

            if is_static:
                disp = class_obj.add_static_case(
                    table_name, elementi, data, isubcase, # freq,
                    element_name,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_transient:
                times = transient_data
                #if not hasattr(class_obj, 'add_transient_case'):
                    #model.log.warning(f'skipping {result_name} - add_transient_case')
                    #continue
                disp = class_obj.add_transient_case(
                    table_name, elementi, data, isubcase,
                    element_name, times,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')

            elif is_modes:
                #model.log.warning(f'skipping {result_name} - add_modal_case')
                modes, eigns, cycles = modes_data
                disp = class_obj.add_modal_case(
                    table_name, elementi, data, isubcase,
                    element_name, modes, eigns, cycles,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_post_buckling:
                modes, eigrs, eigis = buckling_data
                #if not hasattr(class_obj, 'add_post_buckling_case'):
                    #model.log.warning(f'skipping {result_name} - add_post_buckling_case')
                    #continue
                disp = class_obj.add_post_buckling_case(
                    table_name, elementi, data, isubcase,
                    element_name, modes, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
                continue
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _complex_composite_shell(result_name: str,
                             iresult: int,
                             results: dict[int, Function],
                             domains_df: pd.DataFrame,
                             group: h5py._hl.dataset.Dataset,
                             index: h5py._hl.dataset.Dataset,
                             version: int,
                             ids: np.ndarray,
                             model: OP2,
                             subcases: Optional[list[int]]=None) -> int:
    if result_name in {'stress.ctria3_composite_stress', 'stress.ctria6_composite_stress', 'stress.ctriar_composite_stress',
                       'stress.cquad4_composite_stress', 'stress.cquad8_composite_stress', 'stress.cquadr_composite_stress'}:
        #class_obj = RealCompositePlateStressArray
        table_name = 'OES1C'
    elif result_name in {'strain.ctria3_composite_strain', 'strain.ctria6_composite_strain', 'strain.ctriar_composite_strain',
                         'strain.cquad4_composite_strain', 'strain.cquad8_composite_strain', 'strain.cquadr_composite_strain'}:
        #class_obj = RealCompositePlateStrainArray
        table_name = 'OSTR1C'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, composite, stress_strain = result_name.upper().split('_')
    #is_strain = ('STRAIN' == stress_strain)

    basename = result_name + ' (real)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'


    #dtype=[('EID', '<i8'), ('PLY', '<i8'), ('X1', '<f8'), ('Y1', '<f8'),
    #       ('T1', '<f8'), ('L1', '<f8'), ('L2', '<f8'), ('DOMAIN_ID', '<i8')])
    assert version == 0, version
    EID = group['EID']
    nlayers = len(EID)

    element_layer = np.zeros((nlayers, 2), dtype=EID.dtype)
    element_layer[:, 0] = EID
    element_layer[:, 1] = group['PLY']

    ('EID', 'PLY', 'X1R', 'Y1R', 'T1R', 'L1R', 'L2R', 'X1I', 'Y1I', 'T1I', 'L1I', 'L2I', 'DOMAIN_ID')
    X1 = group['X1R'] + 1j * group['X1I']
    Y1 = group['Y1R'] + 1j * group['Y1I']
    T1 = group['T1R'] + 1j * group['T1I']  # T12

    # interlaminar shear
    L1 = group['L1R'] + 1j * group['L1I']  # TXZ-mat
    L2 = group['L2R'] + 1j * group['L2I']  # TYZ-mat
    #angle = X1 * np.nan
    #major = X1 * np.nan
    #minor = X1 * np.nan
    #ovm = X1 * np.nan

    #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
    #if X1.ndim == 1:
        #nnodes_per_element = 1
    #else:
        #nnodes_per_element = X1.shape[1]

    DATA = hstack_shape(nlayers,
                        X1, Y1, T1, L1, L2)
    mode, eigr, eigi, data_by_group = _get_data_by_group_element(
        basename, domains_df, index,
        element_layer, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            idomain = []
            nelement = []
            ntimes = len(group)
            for itime, (idomaini, element_layeri, datai) in group.items():
                idomain.append(idomaini)
                nelementi = len(element_layeri)
                nelement.append(nelementi)
                datas.append(datai)
            assert max(nelement) == min(nelement)
            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            data = np.stack(datas, axis=0)
            nelements = len(element_layeri)
            assert data.shape == (ntimes, nelements, 5), data.shape
            assert element_layeri.shape == (nelements, 2), element_layeri.shape

            if is_freq:
                freq = freq_data
                #model.log.warning(f'skipping {result_name} - add_freq_case')
                #nelements_, nlayers_, two_ = element_layeri.shape
                #element_layeri = element_layeri.reshape(nelements_*nlayers_, 2)
                model.log.warning(f'skipping {result_name} - add_freq_case')
                return
                disp = class_obj.add_freq_case(
                    table_name, element_name, element_layeri, data, isubcase, freq,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_complex_modes:
                model.log.warning(f'skipping {result_name} - add_complex_modes_case')
                continue
            #if is_static:
                #disp = class_obj.add_static_case(
                    #table_name, elementi, data, isubcase, # freq,
                    #element_name,
                    #is_sort1=True, is_random=False, is_msc=True,
                    #random_code=0, title='', subtitle='', label='')
            #elif is_transient:
                #times = transient_data
                ##if not hasattr(class_obj, 'add_transient_case'):
                    ##model.log.warning(f'skipping {result_name} - add_transient_case')
                    ##continue
                #disp = class_obj.add_transient_case(
                    #table_name, elementi, data, isubcase,
                    #element_name, times,
                    #is_sort1=True, is_random=False, is_msc=True,
                    #random_code=0, title='', subtitle='', label='')

            #elif is_modes:
                ##model.log.warning(f'skipping {result_name} - add_modal_case')
                #modes, eigns, cycles = modes_data
                #disp = class_obj.add_modal_case(
                    #table_name, elementi, data, isubcase,
                    #element_name, modes, eigns, cycles,
                    #is_sort1=True, is_random=False, is_msc=True,
                    #random_code=0, title='', subtitle='', label='')
            #elif is_buckling:
                #model.log.warning(f'skipping {result_name} - add_buckling_case')
                #continue
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _complex_shell_stress(result_name: str,
                          iresult: int,
                          results: dict[int, Function],
                          domains_df: pd.DataFrame,
                          group: h5py._hl.dataset.Dataset,
                          index: h5py._hl.dataset.Dataset,
                          version: int,
                          ids: np.ndarray,
                          model: OP2,
                          subcases: Optional[list[int]]=None) -> int:
    if result_name in {'stress.ctria3_stress', 'stress.ctria6_stress', 'stress.ctriar_stress',
                       'stress.cquad8_stress', 'stress.cquadr_stress'}:
        class_obj = ComplexPlateStressArray
        table_name = 'OES1'
    elif result_name in {'strain.ctria3_strain', 'strain.ctria6_strain', 'strain.ctriar_strain',
                         'strain.cquad8_strain', 'strain.cquadr_strain'}:
        class_obj = ComplexPlateStrainArray
        table_name = 'OSTR1'
    else:  # pragma: no cover
        raise RuntimeError(result_name)

    element_name, stress_strain = get_element_name_from_result_name(result_name)
    #is_strain = ('STRAIN' == stress_strain)

    basename = result_name + ' (complex)'
    # TODO: check real/imaginary or magnitude/phase
    #'ID', 'X', 'Y', 'Z', 'RX', 'RY', 'RZ', 'DOMAIN_ID',
    #'DOMAIN_ID', 'POSITION', 'LENGTH'

    # CTRIA3_CPLX / CQUAD4_CPLX???
    # dtype=[('EID', '<i8'), ('FD1', '<f8'), ('X1R', '<f8'), ('X1I', '<f8'), ('Y1R', '<f8'), ('Y1I', '<f8'), ('TXY1R', '<f8'), ('TXY1I', '<f8'),
    #                        ('FD2', '<f8'), ('X2R', '<f8'), ('X2I', '<f8'), ('Y2R', '<f8'), ('Y2I', '<f8'), ('TXY2R', '<f8'), ('TXY2I', '<f8'), ('DOMAIN_ID', '<i8')])
    #
    assert version == 0, version
    EID = group['EID']
    nelements = len(EID)
    nlength = index['LENGTH'].sum()
    #nueid = len(np.unique(EID))

    #element_node = np.zeros((nelements, 2), dtype=EID.dtype)
    #element_node[:, 0] = EID

    #names = group.dtype.names
    X1 = group['X1R'] + group['X1I'] * 1j
    Y1 = group['Y1R'] + group['Y1I'] * 1j
    TXY1 = group['TXY1R'] + group['TXY1I'] * 1j

    X2 = group['X2R'] + group['X2I'] * 1j
    Y2 = group['Y2R'] + group['Y2I'] * 1j
    TXY2 = group['TXY1R'] + group['TXY2I'] * 1j

    if X1.ndim == 1:
        nnodes_per_element = 1
    else:
        nnodes_per_element = X1.shape[1]
    #element_node = element_node.reshape(nelements, nnodes_per_element, 2)

    nlayers_half = nelements * nnodes_per_element
    nlayers = nlayers_half * 2


    #FIBER = np.hstack([
        #group['FD1'].reshape(nlayers_half, 1),
        #group['FD2'].reshape(nlayers_half, 1),
    #]).reshape(nlayers, 1)
    if element_name == 'CTRIA3':
        nnodes = 1
        element_node = np.zeros((nlayers, 2), dtype=EID.dtype)
        element_node[:, 0] = np.repeat(EID, 2)

        #nfiber_group = len(group['FD1'])
        FIBER = np.hstack([
            group['FD1'].reshape(nlayers_half, 1),
            group['FD2'].reshape(nlayers_half, 1),
        ])
        DATA = np.hstack([
            X1.reshape(nlayers_half, 1),
            Y1.reshape(nlayers_half, 1),
            TXY1.reshape(nlayers_half, 1),
            X2.reshape(nlayers_half, 1),
            Y2.reshape(nlayers_half, 1),
            TXY2.reshape(nlayers_half, 1),
        ])# .reshape(nlayers, 3)
    else:
        nnodes = group['GRID'].shape[1]
        # ('EID', 'TERM', 'GRID',
        #  'FD1', 'X1R', 'X1I', 'Y1R', 'Y1I', 'TXY1R', 'TXY1I',
        #  'FD2', 'X2R', 'X2I', 'Y2R', 'Y2I', 'TXY2R', 'TXY2I', 'DOMAIN_ID')
        element_node = np.full((nlayers, 2), -1, dtype=EID.dtype)
        element_node[:, 0] = np.repeat(EID, 2*nnodes)

        ## TODO: element_node seems very wrong...why are there no zeros in it...
        ## for a CQUAD8 with 5 times...where is the 0?  why is the 7 not 0?
        ## array([
        ## [ 7,  4, 15, 61, 60],
        ## [ 7,  4, 15, 61, 60],
        ## [ 7,  4, 15, 61, 60],
        ## [ 7,  4, 15, 61, 60],
        ## [ 7,  4, 15, 61, 60]], dtype=int64)

        grid = group['GRID']
        grid[:, 0] = 0  # why is the first column not 0???
        #element_node[:, 1] = np.repeat([grid.ravel(), grid.ravel()])  #old
        element_node[:, 1] = np.repeat(grid.ravel(), 2)  #old

        #grid = np.hstack([
            #np.zeros(nelements, 1),
            #group['GRID'],
        #])
        #element_node[:, 1] = np.vstack([grid, grid])

        # 2 layers
        # 2 columns (eid, nid)
        element_node = element_node.reshape((nelements, 4*nnodes))
        FIBER = np.stack([
            group['FD1'],
            group['FD2'],
        ], axis=2)
        DATA = np.stack([
            X1, Y1, TXY1,
            X2, Y2, TXY2,
        ], axis=2)
    assert FIBER.shape[0] == nlength, f'fiber.shape={FIBER.shape}; nlength={nlength}'
    element_node = element_node.reshape(nelements, nnodes*2, 2)
    nresults = 3

    assert DATA.shape[0] == nlength, f'DATA.shape={DATA.shape}; nlength={nlength}'
    mode, eigr, eigi, data_by_group = _get_data_by_group_fiber_element(
        basename, domains_df, index,
        element_node, FIBER, DATA)

    for isubcase, data_analysis in sorted(data_by_group.items()):
        assert len(data_analysis) == 1, data_analysis
        for analysis_code, group in data_analysis.items():
            #element = []
            datas = []
            fibers = []
            idomain = []
            nelement = []
            ntimes = len(group)
            for itime, (idomaini, element_nodei, fiberi, datai) in group.items():
                idomain.append(idomaini)
                nelementi = element_nodei.shape[0]
                nnodei = element_nodei.shape[1]
                nelement.append(nelementi)
                fibers.append(fiberi.flatten())
                datas.append(datai)
            assert max(nelement) == min(nelement)
            element_nodei = element_nodei.reshape(nelementi*nnodei, 2)
            assert element_nodei[0, 0] == element_nodei[1, 0]
            assert element_nodei[0, 1] == element_nodei[1, 1]

            #data = data.reshape(ntimes, 2*nelementi, 8) # real
            data = np.stack(datas, axis=0)
            if element_name == 'CTRIA3':
                ntimes_, nlayers_per_element_over_2, nresults_times_2 = data.shape
                nlayers_ = nlayers_per_element_over_2 * 2
            else:
                ntimes_, nelements_, nlayers_per_element_over_2, nresults_times_2 = data.shape
                nlayers_ = nelements_ * nlayers_per_element_over_2 * 2
            data = data.reshape((ntimes, nlayers_, nresults))
            assert data.ndim == 3, f'data.shape={data.shape}'

            (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
             static_data, modes_data, freq_data, buckling_data,
             transient_data, complex_modes) = get_analysis_code_times(
                analysis_code, mode, eigr, eigi, idomain)

            fiber = fibers[0]
            #nelements = len(elementi)
            nfiber = len(fiber)
            assert data.shape == (ntimes, nfiber, 3), data.shape

            #nelements_, nlayers_, two_ = element_nodei.shape
            #element_nodei = element_nodei.reshape(nelements_*nlayers_, 2)
            if is_freq:
                freq = freq_data
                #model.log.warning(f'skipping {result_name} - add_freq_case')
                disp = class_obj.add_freq_case(
                    table_name, element_name, element_nodei, fiber, data, isubcase, freq,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
            elif is_complex_modes:
                #model.log.warning(f'skipping {result_name} - add_complex_modes_case')
                modes, eigrs, eigis = complex_modes
                #model.log.warning(f'skipping {result_name} - add_freq_case')
                disp = class_obj.add_complex_modes_case(
                    table_name, element_name, element_nodei, fiber, data, isubcase,
                    modes, eigrs, eigis,
                    is_sort1=True, is_random=False, is_msc=True,
                    random_code=0, title='', subtitle='', label='')
                continue
            else:  # pragma: no cover
                raise NotImplementedError(analysis_code)
            model.subcase_key[isubcase].append(isubcase)
            model.get_result(result_name)[isubcase] = disp
            #model.displacements[isubcase] = disp
    return iresult

def _get_data_by_group_element(basename: str,
                               domains_df: pd.Dataframe,
                               index: Group,
                               EID: np.ndarray,
                               DATA: np.ndarray) -> tuple[
                                   pd.Series, pd.Series, pd.Series,
                                   dict[int,
                                        dict[int,
                                             dict[int,
                                                  tuple[int, np.ndarray, np.ndarray]]]]]:
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

        #idomain = []
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

        #iresults = np.full(len(idomain), np.nan, dtype='int32')
        #is_modes = np.abs(mode).max() != 0
        #is_freq = np.abs(eigr).max() != 0 or np.abs(eigi).max() != 0
        #assert is_modes or is_freq
        #assert is_modes ^ is_freq, f'is_modes={is_modes} is_freq={is_freq} are both True/False...'  # xor

        data_by_group = _preallocate_data_by_group(subcase, analysis_code)
        for itime, subcasei, analysis_codei, modei, eigri, eigii, idomaini, positioni, lengthi in zip(
                count(), subcase, analysis_code, mode, eigr, eigi, idomain, position, length):
            #nelementsi = lengthi
            i0 = positioni
            i1 = positioni + lengthi
            #name = get_name(basename, is_freq, subcasei, analysis_codei, modei, eigri, eigii)

            #print('  ', name)
            #print(itime, domain, positioni, lengthi)

            # make it so we can determine the other "times"
            #iresults[itime] = iresult
            #_domain = DOMAIN[i0]

            elementi = EID[i0:i1]
            datai = DATA[i0:i1]
            #print(elementi[:, 0])
            #op2.displacements[subcasei]
            #DataByGroup =
            data_by_group[subcasei][analysis_codei][itime] = (idomaini, elementi, datai)

            #results[iresult] = RealVectorTable(
                #name, itime, iresult, iresults,
                #_domain, position, length,
                #NID[i0:i1],
                #TX[i0:i1], TY[i0:i1], TZ[i0:i1],
                #RX[i0:i1], RY[i0:i1], RZ[i0:i1],
                #DOMAIN[i0:i1], location='node')
            #iresult += 1
    return mode, eigr, eigi, data_by_group

def get_idomain(grouped_dfi: pd.DataFrame, INDEX_DOMAIN: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    no_domain = True
    indexi, dfi = grouped_dfi
    DOMAINs = dfi['ID']
    #print(indexi)
    #print(dfi)
    #print('---------------------')
    idomain = np.searchsorted(DOMAINs, INDEX_DOMAIN)

    # doesn't handle lower out of range
    exists = (idomain < len(INDEX_DOMAIN))

    # condition #2 out of range
    #exists2 = ((idomain < len(INDEX_DOMAIN) & (DOMAINs.values[idomain] == INDEX_DOMAIN))

    #  why doesn't this work???
    #exists2 = np.array([idomaini if existi for idomaini, existi in zip(idomain, exists)])
    exists2 = []
    for idomaini, existi in zip(idomain, exists):
        if existi:
            exists2.append(idomaini)
    exists2 = np.array(exists2, dtype='bool')

    if not np.all(exists):
        #if np.any(exists2):
            #raise RuntimeError(idomain)
        return no_domain, idomain, dfi

    no_domain = False
    return no_domain, idomain, dfi

def _get_data_by_group_fiber_element(basename: str,
                                     domains_df: pd.Dataframe, index: Group,
                                     EID: np.ndarray,
                                     FIBER: np.ndarray,
                                     DATA: np.ndarray) -> dict[int, dict[int, dict[int, GroupData]]]:
    INDEX_DOMAIN = index['DOMAIN_ID']
    INDEX_POSITION = index['POSITION']
    INDEX_LENGTH = index['LENGTH']
    total_length = INDEX_LENGTH.sum()
    assert len(EID) == total_length, f'EID.shape={EID.shape} total_length={total_length}'
    assert len(FIBER) == total_length, f'FIBER.shape={FIBER.shape} total_length={total_length}'
    assert len(DATA) == total_length, f'DATA.shape={DATA.shape} total_length={total_length}'

    ntotal = 0

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
        domains: np.ndarray = grouped_dfi[1]['ID'].values
        no_domain, idomain, dfi = get_idomain(grouped_dfi, INDEX_DOMAIN)
        if no_domain:
            continue

        #idomain = []
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

        #iresults = np.full(len(idomain), np.nan, dtype='int32')
        #is_modes = np.abs(mode).max() != 0
        #is_freq = np.abs(eigr).max() != 0 or np.abs(eigi).max() != 0
        #assert is_modes or is_freq
        #assert is_modes ^ is_freq, f'is_modes={is_modes} is_freq={is_freq} are both True/False...'  # xor
        data_by_group = _preallocate_data_by_group(subcase, analysis_code)

        for itime, subcasei, analysis_codei, modei, eigri, eigii, idomaini, positioni, lengthi in zip(
                count(), subcase, analysis_code, mode, eigr, eigi, idomain, position, length):
            #nelementsi = lengthi
            i0 = positioni
            i1 = positioni + lengthi
            #name = get_name(basename, is_freq, subcasei, analysis_codei, modei, eigri, eigii)

            #print('  ', name)
            #print(itime, domain, positioni, lengthi)

            # make it so we can determine the other "times"
            #iresults[itime] = iresult
            #_domain = DOMAIN[i0]

            elementi = EID[i0:i1]
            #print(f'i0={i0} i1={i1}')
            #print(elementi[:, 0])
            #print(f'element[0]={elementi[0,0]} element[-1]={elementi[-1,0]}; '
                  #f'min={elementi[:,0].min()}; max={elementi[:,0].max()}')
            #print()
            fiberi = FIBER[i0:i1]
            datai = DATA[i0:i1]
            #op2.displacements[subcasei]

            data_by_group[subcasei][analysis_codei][itime] = (idomaini, elementi, fiberi, datai)
            ntotal += len(elementi)

            #results[iresult] = RealVectorTable(
                #name, itime, iresult, iresults,
                #_domain, position, length,
                #NID[i0:i1],
                #TX[i0:i1], TY[i0:i1], TZ[i0:i1],
                #RX[i0:i1], RY[i0:i1], RZ[i0:i1],
                #DOMAIN[i0:i1], location='node')
            #iresult += 1
    assert ntotal == total_length, f'total_length={total_length} ntotal={ntotal}'
    return mode, eigr, eigi, data_by_group

def _preallocate_data_by_group(subcase: pd.Series,
                               analysis_code: int) -> dict[int, dict[int, dict[Any, Any]]]:
    data_by_group = {}
    #data_by_group[subcase][analysis_codei]
    for subcasei in subcase:
        data_by_group[subcasei] = {}
    for subcasei, analysis_codei in zip(subcase, analysis_code):
        data_by_group[subcasei][analysis_codei] = {}
    return data_by_group

def get_element_name_from_result_name(result_name: str) -> tuple[str, str]:
    element_name, stress_strain = result_name.upper().split('_')
    if '.' in element_name:
        element_name = element_name.split('.')[1]
    return element_name, stress_strain
