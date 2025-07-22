from __future__ import annotations
import numpy as np
import h5py
from typing import TYPE_CHECKING

from pyNastran.utils.dict_to_h5py import _cast # , _cast_array
from pyNastran.bdf.bdf_interface.subcase.cards import OBJ_MAP
from pyNastran.bdf.bdf_interface.case_control_cards import CLASS_MAP
from pyNastran.utils.numpy_utils import integer_types, integer_float_types
from pyNastran.utils import object_attributes

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.subcase import Subcase


def load_hdf5_file(subcase: Subcase, hdf5_file, encoding: str):
    log: SimpleLogger = subcase.log

    keys = list(hdf5_file.keys())
    for key in keys:
        #print(key)
        group = hdf5_file[key]
        if key in ['id']: # scalars
            value = _cast(group)
            setattr(subcase, key, value)
        elif key == 'params':
            #print(group)
            #print(group.keys())
            for group_key in group.keys():
                #log.debug('%s %s' % (group_key, key))
                value, options, param_type = _load_hdf5_param(group, group_key, encoding, log)
                if param_type == 'OBJ-type':
                    subcase.params[group_key] = value
                    str(subcase)
                    continue

                #log.debug('%s (%s, %s, %s)' % (key, value, options, param_type))
                if isinstance(options, list):
                    options = [
                        option.decode(encoding) if isinstance(option, bytes) else option
                        for option in options]

                subcase.params[group_key] = (value, options, param_type)
                str(subcase)
        else:  # pragma: no cover
            raise RuntimeError(f'failed exporting Subcase/{key}')

def _load_hdf5_param(group, key: str, encoding: str, log: SimpleLogger) -> tuple[str, list[str], str]:
    """('ALL', ['SORT2'], 'STRESS-type')"""
    #print('-----------------------------------------')
    #print(type(key), key)
    sub_group = group[key]
    keys = list(sub_group.keys())
    #print(group)
    #print('subgroup.keys() =', sub_group.keys())

    if key == 'blank':
        key = ''

    if 'options' in sub_group:
        keys.remove('options')
        options = _cast(sub_group['options'])
        if isinstance(options, (integer_types, str)):
            pass
        elif isinstance(options, bytes):
            options = options.decode(encoding)
        else:
            options = [
                option.decode(encoding) if isinstance(option, bytes) else option
                for option in options]
    else:
        options = None

    param_type = None
    if 'param_type' in sub_group:
        param_type = _cast(sub_group['param_type'])
        if isinstance(param_type, bytes):
            param_type = param_type.decode(encoding)
        keys.remove('param_type')

    value = None
    if param_type == 'OBJ-type':
        #sub_groupi.attrs['type'] = key
        Type = sub_group.attrs['type']
        try:
            cls = OBJ_MAP[Type]
        except KeyError:
            log.error(f'*param_type={param_type}: Type={Type} key={key}')
            raise

        obj = cls.load_hdf5(sub_group, encoding)
        str(obj)
        return obj, None, param_type

    if 'value' in sub_group:
        keys.remove('value')
        value = _cast(sub_group['value'])
        if isinstance(value, bytes):
            value = value.decode(encoding)
        elif isinstance(value, (integer_types, str)):
            pass
        elif isinstance(value, np.ndarray):
            value = value.tolist()
        #else:
            #value = value.tolist()

    elif 'object' in sub_group:
        value, options = _load_hdf5_object(key, keys, sub_group, encoding)


    if len(keys) > 0:
        #keyi = _cast(sub_group['key'])
        #print('keyi = %r' % keyi)
        raise RuntimeError(f'keys = {keys}')

    #print(value, options, param_type)
    return value, options, param_type

def _load_hdf5_object(key, keys, sub_group, encoding):
    keys.remove('object')
    sub_groupi = sub_group['object']

    Type = sub_groupi.attrs['type']
    cls = CLASS_MAP[Type]

    if hasattr(cls, 'load_hdf5'):
        class_obj, options = cls.load_hdf5(sub_groupi, 'utf8')
        value = class_obj
        return value, options

    use_data = True
    if 'options' in sub_groupi:
        options2 = _cast(sub_groupi['options'])
        value = _cast(sub_groupi['value'])
        #print('sub_keys =', sub_groupi, sub_groupi.keys())

        options_str = [
            option.decode(encoding) if isinstance(option, bytes) else option
            for option in options2]
        use_data = False
    if isinstance(value, bytes):
        value = value.decode(encoding)

    data_group = sub_groupi['data']
    keys2 = _cast(data_group['keys'])

    h5_values = data_group['values']
    if isinstance(h5_values, h5py._hl.group.Group):
        values2 = [None] * len(keys2)
        for ih5 in h5_values.keys():
            ih5_int = int(ih5)
            h5_value = _cast(h5_values[ih5])
            values2[ih5_int] = h5_value
    else:
        values2 = _cast(h5_values)
    #print('data_keys =', data_group, data_group.keys())

    unused_keys_str = [
        keyi.decode(encoding) if isinstance(keyi, bytes) else keyi
        for keyi in keys2]
    unused_values_str = [
        valuei.decode(encoding) if isinstance(valuei, bytes) else valuei
        for valuei in values2]

    #print('keys2 =', keys2)
    #print('values2 =', values2)
    #print('options2 =', options2)

    if use_data:
        #print('keys2 =', keys2)
        #print('values2 =', values2)
        data = []
        for keyi, valuei in zip(keys2, values2):
            data.append((keyi, valuei))
        class_obj = cls(data)
        #assert options is None, options
    else:
        class_obj = cls(key, value, options_str)
        options = options_str
    value = class_obj
    #print(class_obj)
    #class_obj.load_hdf5_file(hdf5_file, encoding)
    return value, options

def export_to_hdf5(subcase: Subcase, hdf5_file, encoding: str) -> None:
    keys_to_skip = ['log', 'solCodeMap', 'allowed_param_types']
    h5attrs = object_attributes(subcase, mode='both', keys_to_skip=keys_to_skip)
    #print(f'Subcase {subcase.id:d}')
    for h5attr in h5attrs:
        value = getattr(subcase, h5attr)
        if h5attr in ['id']: # scalars
            # simple export
            hdf5_file.create_dataset(h5attr, data=value)
        elif h5attr in ['sol']: # scalars/None
            if value is None:
                continue
            hdf5_file.create_dataset(h5attr, data=value)
        elif h5attr in ['params']:
            if len(value) == 0:
                continue
            keys = list(subcase.params.keys())
            params_group = hdf5_file.create_group('params')
            #print('keys =', keys)
            unused_keys_bytes = [key.encode(encoding) for key in keys]
            #params_group.create_dataset('keys', data=keys_bytes)
            for key, (value, options, param_type) in subcase.params.items():
                #print('  %-14s: %-8r %r %r' % (key, value, options, param_type))
                _export_hdf5_param(params_group, key, value, options, param_type, encoding)

        #if h5attr in ['_begin_count', 'debug', 'write_begin_bulk']: # scalars
            ## do nothing on purpose
            #hdf5_file.create_dataset(h5attr, data=value)
        #elif h5attr in ['reject_lines', 'begin_bulk', 'lines', 'output_lines']:
            # lists of strings
            #if len(value) == 0:
                #continue
            #value_bytes = [line.encode(encoding) for line in value]
            ##print(value_bytes)
            #hdf5_file.create_dataset(h5attr, data=value_bytes)
        #elif h5attr == 'subcases':
            #keys = list(self.subcases.keys())
            #subcase_group = hdf5_file.create_group('subcases')
            #subcase_group.create_dataset('keys', data=keys)
            #for key, subcase in self.subcases.items():
                #sub_group = subcase_group.create_group(str(key))
                #subcase.export_to_hdf5(subcase_group, encoding)
        else:  # pragma: no cover
            print(key, value)
            raise RuntimeError(f'cant export to hdf5 Subcase/{h5attr}')

def _export_hdf5_param(params_group,
                       key: str, value, options, param_type: str,
                       encoding: str):
    """
    Export a Subcase param like:
    K2GG = 1.0*K2GG1
    SET 10 = 1 THRU 10
    """
    if param_type == 'OBJ-type':
        #print('**********', key, param_type)
        sub_groupi = params_group.create_group(key)
        sub_groupi.attrs['type'] = key
        sub_groupi['param_type'] = param_type
        value.export_to_hdf5(sub_groupi, encoding)
        return

    if key == '':
        sub_group = params_group.create_group('blank')
        sub_group.create_dataset('value', data=value)

    else:
        #print('key = %r' % key)
        sub_group = params_group.create_group(key)
        if value is not None:
            if isinstance(value, list):
                value_bytes = [
                    valuei.encode(encoding) if isinstance(valuei, str) else valuei
                    for valuei in value]
                sub_group.create_dataset('value', data=value_bytes)
            elif isinstance(value, (integer_float_types, str)):
                sub_group.create_dataset('value', data=value)
            elif hasattr(value, 'export_to_hdf5'):
                sub_groupi = sub_group.create_group('object')
                sub_groupi.attrs['type'] = key
                value.export_to_hdf5(sub_groupi, encoding)
            else:
                print(f'value = {value!r}')
                raise NotImplementedError(value)

    if param_type is not None:
        sub_group.create_dataset('param_type', data=param_type)

    if options is not None:
        if isinstance(options, list):
            options_bytes = [
                option.encode(encoding) if isinstance(option, str) else option
                for option in options]
            #print('options_bytes', options_bytes)
            sub_group.create_dataset('options', data=options_bytes)
        else:
            sub_group.create_dataset('options', data=options)
