import os
import numpy as np
import h5py
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray, ComplexDisplacementArray
from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray, ComplexVelocityArray
from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray, ComplexAccelerationArray
from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray, ComplexEigenvectorArray

from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray, ComplexSPCForcesArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import RealMPCForcesArray, ComplexMPCForcesArray
#from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import RealTemperatureGradientAndFluxArray
from pyNastran.utils import print_bad_path

def cast(h5_result_attr):
    if len(h5_result_attr.shape) == 0:
        return np.array(h5_result_attr).tolist()
        #raise NotImplementedError(h5_result_attr.dtype)
    else:
        return np.array(h5_result_attr)

TABLE_OBJ_MAP = {
    'displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'velocities' : (RealVelocityArray, ComplexVelocityArray),
    'accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'eigenvectors' : (RealEigenvectorArray, ComplexEigenvectorArray),
    'spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
}
TABLE_OBJ_KEYS = list(TABLE_OBJ_MAP.keys())

def load_table(model, h5_result, real_obj, complex_obj, debug=False):
    """loads a RealEigenvectorArray/ComplexEigenvectorArray"""
    is_real = cast(h5_result.get('is_real'))
    is_complex = cast(h5_result.get('is_complex'))
    nonlinear_factor = cast(h5_result.get('nonlinear_factor'))

    data_names = cast(h5_result.get('data_names')).tolist()
    sdata_names = [data_name + 's' for data_name in data_names]
    data_code = {
        'nonlinear_factor' : nonlinear_factor,
        'sort_bits' : cast(h5_result.get('sort_bits')).tolist(),
        'sort_method' : cast(h5_result.get('sort_method')),
        'is_msc' : cast(h5_result.get('is_msc')),
        'format_code' : cast(h5_result.get('format_code')),
        'device_code' : cast(h5_result.get('device_code')),
        'approach_code' : cast(h5_result.get('approach_code')),
        'analysis_code' : cast(h5_result.get('analysis_code')),
        'table_code' : cast(h5_result.get('table_code')),
        'tCode' : cast(h5_result.get('tCode')),
        'sort_code' : cast(h5_result.get('sort_code')),
        'thermal' : cast(h5_result.get('thermal')),
        'subtitle' : cast(h5_result.get('subtitle')),
        'acoustic_flag' : cast(h5_result.get('acoustic_flag')),
        'data_names' : data_names,
        'name' : data_names[0],
    }
    is_sort1 = cast(h5_result.get('is_sort1'))
    isubcase = cast(h5_result.get('isubcase'))
    dt = nonlinear_factor

    if is_real:
        obj = real_obj(
            data_code, is_sort1, isubcase, dt,
            #f06_flag=False,
        )
    else:
        return None

    keys_to_skip = [
        'class_name', 'headers', 'is_real', 'is_complex',
        'is_sort1', 'is_sort2', 'table_name_str']
    #time_keys = ['eigns', 'mode_cycles', 'modes']
    for key in h5_result.keys():
        if key in keys_to_skip:
            continue
        elif key in sdata_names:
            if debug:  # pragma: no cover
                print('  *****key=%r' % key)
            datai = cast(h5_result.get(key))
            setattr(obj, key, datai)
            setattr(obj, '_times', datai)
        elif key not in data_code:
            if debug:  # pragma: no cover
                print('  **key=%r' % key)
            datai = cast(h5_result.get(key))
            setattr(obj, key, datai)
    return obj

def load_op2_from_h5(h5_filename, log=None):
    assert os.path.exists(h5_filename), print_bad_path(h5_filename)
    model = OP2(log=None)
    model.op2_filename = h5_filename

    h5_file = h5py.File(h5_filename, 'r')
    debug = False
    for key in h5_file.keys():
        if key.startswith('Subcase'):
            h5_subcase = h5_file.get(key)
            log.info('subcase:')
            for result_name in h5_subcase.keys():
                if result_name in ['eigenvalues']:
                    log.warning('    skipping %r...' % result_name)
                elif result_name in TABLE_OBJ_KEYS:
                    if debug:
                        log.info('  %s:' % result_name)
                    real_obj, complex_obj = TABLE_OBJ_MAP[result_name]
                    h5_result = h5_subcase.get(result_name)
                    obj = load_table(model, h5_result, real_obj, complex_obj, debug=debug)
                    if obj is None:
                        log.warning('    skipping %r...' % result_name)
                        continue

                    # isubcase, analysis_code, sort_method,
                    #  count, ogs, superelement_adaptivity_index
                    opt_count = 0
                    ogs = 0
                    superelement_adaptivity_index = ''
                    # (1, 2, 1, 0, 0, u'')
                    key = (obj.isubcase, obj.analysis_code, obj.sort_method,
                           opt_count, ogs,
                           superelement_adaptivity_index)
                    model.eigenvectors[key] = obj
                    log.info('  loaded %r' % result_name)
                else:
                    log.warning('skipping %r...' % result_name)

            #print(h5_subcase.keys())
        elif key == 'info':
            pass
        elif key == 'matrices':
            log.info('matrices:')
            h5_matrix_group = h5_file.get(key)
            for matrix_name in h5_matrix_group.keys():
                h5_matrix = h5_matrix_group.get(matrix_name)
                nkeys = len(h5_matrix.keys())
                if not nkeys:
                    log.warning('  %s is empty...skipping' % h5_matrix)
                else:
                    log.warning('  skipping %r...' % matrix_name)
                    #for attr in h5_matrix.keys():
                        #print('    attr=%s' % attr)

        else:
            print('key = %r' % key)

    return model