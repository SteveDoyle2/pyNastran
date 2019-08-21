# coding: utf-8
"""
defines:
 - mydict = load_obj_from_hdf5(hdf5_filename, log=None, debug=False)
 - mydict = load_obj_from_hdf5_file(mydict, h5_file, log=None, debug=False)
 - export_obj_to_hdf5(hdf5_filename, mydict)
 - export_obj_to_hdf5_file(hdf5_file, mydict)

Supports:
- integers, floats, None, strings, unicode, lists, tuple,
- numpy arrays (including NaN),
- objects (including custom objects),
- scikit-learn StandardScalar
- minimal dependencies (e.g., no scikit-learn)

Limitations:
- Dictionary keys must be strings/unicode
- May run into problems if you have two classes with the same name,
  but point to different locations.  There is some support for this,
  but hopefully you aren’t using it.

"""
#from types import MethodType, FunctionType

import h5py
import numpy as np
from cpylog import get_logger2

from pyNastran.utils import object_attributes, check_path
from pyNastran.utils.numpy_utils import integer_types, float_types

string_types = (str, bytes)
#integer_types = (int, np.int32, np.int64)
#float_types = (float, np.float32, np.float64)


#---------------------------------------------------------------------------------------------

def export_obj_to_hdf5(hdf5_filename, obj, user_custom_types=None, log=None, debug=False):
    """exports an object to an HDF5 file"""
    log = get_logger2(log=log, debug=debug, encoding='utf-8')
    with h5py.File(hdf5_filename, 'w') as hdf5_file:
        log.info('starting export_op2_to_hdf5_file of %r' % hdf5_filename)
        export_obj_to_hdf5_file(
            hdf5_file, obj, user_custom_types=user_custom_types, log=log)

def export_obj_to_hdf5_file(hdf5_file, obj, user_custom_types=None, log=None, debug=False):
    """exports an object to an HDF5 file object"""
    exporter = HDF5Exporter(hdf5_file, user_custom_types=user_custom_types, log=log, debug=debug)
    exporter._create_dict_group(hdf5_file, obj, exporter.user_custom_types, nlevels=0)

class HDF5Exporter:
    def __init__(self, hdf5_file, user_custom_types=None, log=None, debug=False):
        log = get_logger2(log=log, debug=debug, encoding='utf-8')
        if user_custom_types is None:
            user_custom_types = []
        self.user_custom_types = user_custom_types
        self.hdf5_file = hdf5_file
        self.log = log

    def _create_dict_group(self, hdf5_file, mydict, user_custom_types, nlevels):
        """creates the info HDF5 group"""

        key_types = [key.__class__.__name__ for key in mydict.keys()]
        ukey_types = set(list(key_types))
        #if ukey_types:
            #print('ukey_types =', ukey_types)

        is_int = 'int' in ukey_types
        #if not len(ukey_types) in [0, 1]:
            #print('%skey_types = %s' % (nlevels*' ',  list(key.__class__.__name__ for key in mydict.keys())))
        #if is_int:
            #print('%sis_int = %s' % (nlevels*' ', is_int))

        #is_int = any(isinstance(key, (integer_types, float_types)) for key in mydict.keys())
        if is_int:
            sub_group = hdf5_file.create_group('dict_keys')
            if len(ukey_types):
                sub_group.attrs['key_type'] = list(ukey_types)[0]
            else:
                raise NotImplementedError('variable type key')

            #key_types2 = (str(key) for key in mydict.keys())
            for key, value in sorted(mydict.items()):
                skey = str(key)
                self._add_dataset(hdf5_file, skey, value, user_custom_types, nlevels+1)
        else:
            for key, value in sorted(mydict.items()):
                self._add_dataset(hdf5_file, key, value, user_custom_types, nlevels+1)


    def _add_attrs(self, sub_group, obj, attrs, user_custom_types, nlevels):
        #print(nlevels*' ', sub_group, obj)
        for attr in attrs:
            #print('%sattr %s' % ((nlevels+1)*' ', attr))
            if not hasattr(obj, attr):
                # an attribute is not required, but may exist
                continue
            value = getattr(obj, attr)
            #if not isinstance(value, (dict, list)):
            print('%sattr %s %s' % ((nlevels+1)*' ', attr, type(value)))
            self._add_dataset(sub_group, attr, value, user_custom_types, nlevels+1)

    def _add_dataset(self, hdf5_file, key, value, user_custom_types, nlevels):
        #print(key, type(key))
        #print(key, type(key), value, type(value))
        if value is None:
            # Nones can't be stored, so we create a custom type
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = 'None'
            return

        if isinstance(key, (integer_types, float_types)):
            raise TypeError('key=%r; key must be a string, not %s\nvalue:\n%r' % (key, type(key), value))

        custom_types_list = user_custom_types + ['BDF', 'OP2', 'OP2Geom', 'StandardScaler', 'lil_matrix']
        class_name = value.__class__.__name__

        if isinstance(value, dict):
            try:
                sub_group = hdf5_file.create_group(key)
            except:
                print('key = %s; type=%s' % (key, type(key)))
                raise
            sub_group.attrs['type'] = 'dict'
            self._create_dict_group(sub_group, value, user_custom_types, nlevels+1)

        elif isinstance(value, (integer_types, float_types, string_types, np.ndarray)):
            try:
                hdf5_file.create_dataset(key, data=value)
            except (TypeError, RuntimeError):
                print('key=%r value=%s type=%s' % (key, str(value), type(value)))
                raise

        elif isinstance(value, tuple):
            self._add_list_tuple(hdf5_file, key, value, 'tuple')
        elif isinstance(value, list):
            self._add_list_tuple(hdf5_file, key, value, 'list')
        elif isinstance(value, set):
            self._add_list_tuple(hdf5_file, key, value, 'set')
        elif hasattr(value, 'export_hdf5_file'):
            #print('export_hdf5_file', key)
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = value.__class__.__name__
            value.export_hdf5_file(sub_group, exporter=self)
        elif hasattr(value, 'get_h5attrs'):
            h5attrs = value.get_h5attrs()
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = value.__class__.__name__
            self._add_attrs(sub_group, value, h5attrs, user_custom_types, nlevels+1)
        elif hasattr(value, 'object_attributes'):
            h5attrs = value.object_attributes(mode='both')
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = value.__class__.__name__
            self._add_attrs(sub_group, value, h5attrs, user_custom_types, nlevels+1)
        elif class_name == 'lil_matrix':
            h5attrs = ['dtype', 'shape', 'ndim', 'nnz'] # 'data', 'rows'
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = value.__class__.__name__
            self._add_attrs(sub_group, value, h5attrs, user_custom_types, nlevels+1)
        elif class_name == 'dtype':
            h5attrs = ['dtype']
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = value.__class__.__name__
            self._add_attrs(sub_group, value, h5attrs, user_custom_types, nlevels+1)
        elif class_name in custom_types_list:
            attrs = object_attributes(value, mode='both', keys_to_skip=None)
            sub_group = hdf5_file.create_group(key)
            sub_group.attrs['type'] = class_name
            self._add_attrs(sub_group, value, attrs, user_custom_types, nlevels+1)
        else:
            print('string_types =', string_types)
            print('value =', value)
            msg = (
                'key=%r Type=%r is not in custom_types=%s and does not have:\n'
                ' - export_hdf5_file(h5_file)\n'
                ' - object_attributes()\n'
                ' - get_h5attrs(self)' % (
                key, class_name, custom_types_list))
            raise TypeError(msg)
            #raise TypeError(type(value))

    def _add_list_tuple(self, hdf5_file, key, value, Type):
        """
        tuples are indistinguishable from lists as a dataset,
        so we'll store it as a numpy array, list it, and then tuple it back

        lists/tuples with numpy unicode are special
        """
        _add_list_tuple(hdf5_file, key, value, Type, self.log)

def _add_list_tuple(hdf5_file, key, value, Type, log):
    """
    tuples are indistinguishable from lists as a dataset,
    so we'll store it as a numpy array, list it, and then tuple it back

    lists/tuples with numpy unicode are special
    """
    try:
        sub_group = hdf5_file.create_group(key)
    except ValueError:  # pragma: no cover
        print('key=%s is duplicated' % key)
        print('value = ', value)
        raise
    sub_group.attrs['type'] = Type
    try:
        sub_group.create_dataset('value', data=value)
    except TypeError:
        # contains unicode
        sub_group.attrs['type'] = Type
        for i, valuei in enumerate(value):
            try:
                sub_group.create_dataset(str(i), data=valuei)
            except TypeError:
                log.error('key=%r Type=%r' % (key, Type))
                print('value = %s' % str(value))
                print('value[%i] = %s' % (i, str(valuei)))
                raise

def load_obj_from_hdf5(hdf5_filename, custom_types_dict=None, log=None, debug=False):
    """
    loads an hdf5 file into an object

    Parameters
    ----------
    hdf5_filename : str
       the h5 filename to load
    custom_types_dict : dict[key] : function()
        the custom mapper
    """
    check_path(hdf5_filename, 'hdf5_filename')
    log = get_logger2(log=log, debug=debug, encoding='utf-8')
    log.info('hdf5_filename = %r' % hdf5_filename)

    model = {}
    with h5py.File(hdf5_filename, 'r') as h5_file:
        load_obj_from_hdf5_file(model, h5_file, custom_types_dict=custom_types_dict, log=log, debug=debug)
    return model

def load_obj_from_hdf5_file(model, h5_file, log=None, custom_types_dict=None, debug=False):
    """loads an h5 file object into an dict object"""
    importer = HDF5Importer(h5_file, custom_types_dict=custom_types_dict, log=log, debug=debug)
    importer.load(model, h5_file)



class HDF5Importer:
    def __init__(self, h5_file, custom_types_dict=None, log=None, debug=False):
        self.log = get_logger2(log=log, debug=debug, encoding='utf-8')
        if custom_types_dict is None:
            custom_types_dict = {}
        self.custom_types_dict = custom_types_dict
        self.h5_file = h5_file

    def load(self, model, h5_file, self_obj=None):
        for key in h5_file.keys():
            self._load_value(model, h5_file, key, self.custom_types_dict, self_obj, nlevels=1)

    def _load_value(self, model, h5_file, key, custom_types_dict, self_obj, nlevels):
        value = h5_file.get(key)
        keys = None
        if not hasattr(value, 'attrs'):
            #print('%s%s %s' % ((nlevels)*'  ', key, value))
            #print("%sno_attrs_cast" % ((nlevels)*'  '))
            value2 = self.cast(h5_file, key, value, nlevels)
            model[key] = value2
            return

        # group
        attrs = value.attrs
        keys = list(attrs.keys())
        if not keys:
            #print('%s%s %s' % ((nlevels)*'  ', key, value))
            #print("%sno_keys_cast" % ((nlevels)*'  '))
            value2 = self.cast(h5_file, key, value, nlevels)
            model[key] = value2
            return

        # keys exist
        #print('%s%s %s %s' % ((nlevels)*'  ', key, value, keys))
        Type = None
        if 'type' in keys:
            Type = value.attrs['type']
            #print('%sType=%s' % ((nlevels)*'  ', Type))
            keys.remove('type')

        key_type = None
        if 'key_type' in keys:
            self.log.warning('not handling key_type for %s' % value)
            key_type = value.attrs['key_type']
            #print('%sType=%s' % ((nlevels)*'  ', key_type))
            keys.remove('key_type')
            return

        function = None
        if 'function' in keys:
            function = value.attrs['function']
            keys.remove('function')

        if keys:
            raise NotImplementedError(keys)


        if function is not None:
            _function_data = getattr(self_obj, function)(self, model, h5_file, key, value)
            model[key] = _function_data
            return

        if Type == 'dict':
            self._load_dict(model, h5_file, key, value, custom_types_dict, self_obj,
                            nlevels+1, print_dict=False)
        elif Type == 'None':
            model[key] = None

        elif Type == 'list':
            _list = self._load_mixed_tuple_list(h5_file, key, value, custom_types_dict, self_obj, nlevels+1)
            model[key] = _list
        elif Type == 'tuple':
            _list = self._load_mixed_tuple_list(h5_file, key, value, custom_types_dict, self_obj, nlevels+1)
            model[key] = tuple(_list)
        elif Type == 'set':
            _list = self._load_mixed_tuple_list(h5_file, key, value, custom_types_dict, self_obj, nlevels+1)
            model[key] = set(_list)

        elif Type in custom_types_dict:
            try:
                obj = self._load_custom_type(h5_file, Type, key, value, custom_types_dict,
                                             self_obj, nlevels)
            except:
                msg = ('Cannot load custom type: %s.  Try setting:\n'
                       ' - load_hdf5_file\n'
                       ' - function\n' % (Type))
                self.log.error(msg)
                raise
            model[key] = obj
        else:
            print('%s%s %s %s' % ((nlevels)*'  ', key, value, keys))
            custom_type_keys = list(custom_types_dict.keys())
            custom_type_keys.sort()
            raise TypeError('Type=%r is not in custom_types_dict=%s' % (Type, custom_type_keys))
            #print("%stype_cast Type=%s" % ((nlevels)*'  ', Type))
            #value2 = self.cast(h5_file, key, value, nlevels)
            #model[key] = value2

    @classmethod
    def cast(cls, h5_file, key, value, nlevels):
        """casts a value"""
        # value
        #print('%s****castingA' % (nlevels*'  '))
        #print(key, value)
        try:
            value2 = _cast(h5_file.get(key))
        except AttributeError:
            print(key)
            raise
        #print('%s****%s' % (nlevels*'  ', value2))
        #print('%s  %r : %s %s' % (nlevels*'  ', key, value2, type(value2)))
        return value2

    def _load_mixed_tuple_list(self, h5_file, key, value, custom_types_dict, self_obj, nlevels):
        """
        Lists/tuples are stored as lists if the data doesn't contain unicode.
        Otherwise, they're stored like dictionaries, with string indices that are
        integer values.
        """
        keys = value.keys()
        is_unicode_list = '0' in keys
        if is_unicode_list:
            mylist = self._load_unicode_list(h5_file, key, value, custom_types_dict, self_obj, nlevels)
        else:
            temp_dict = {}
            sub_h5 = value
            self._load_value(temp_dict, sub_h5, 'value', custom_types_dict, self_obj, nlevels+2)
            mylist = temp_dict['value']
        return mylist

    def _load_unicode_list(self, h5_file, key, value, custom_types_dict, self_obj, nlevels):
        """
        We have a dictionary like:
        data = {
           '1' : value1,
           '2' : value2,
           '3' : value3,
        }
        We do this because we need to worry about unicode

        """
        temp_dict = {}
        sub_h5 = value
        self._load_dict(temp_dict, h5_file, key, value, custom_types_dict, self_obj, nlevels, print_dict=False)

        mydict = temp_dict[key]
        nvalues = len(mydict)
        mylist = [None] * nvalues
        for int_key, valuei in mydict.items():
            i = int(int_key)
            mylist[i] = valuei
        return mylist

    def _load_custom_type(self, h5_file, Type, key, value, custom_types, self_obj, nlevels):
        """
        The following custom methods can/should be defined in a class:
         - init_from_empty()
         - _init_from_self(parent)
         - get_custom_types()

        """
        class_instance = custom_types[Type]
        #print('******Type=%r' % Type)
        if hasattr(class_instance, '_init_from_empty'):
            obj = class_instance._init_from_empty()
        elif hasattr(class_instance, '_init_from_self'):
            #print('self_obj', self_obj)
            obj = class_instance._init_from_self(self_obj)
        else:
            try:
                obj = class_instance()
            except:
                self.log.error('%s cannot load with 0 arguments' % Type)
                raise

        temp_dict = {}
        local_custom_types = custom_types
        if hasattr(obj, 'load_hdf5_file'):
            keys = list(value.keys())
            obj.load_hdf5_file(value)
            return obj

        if hasattr(obj, 'get_custom_types'):
            local_custom_types = obj.get_custom_types()
            #print('local_custom_types =', local_custom_types)
        print(key, value)
        self._load_dict(temp_dict, h5_file, key, value, local_custom_types,
                        obj, nlevels+1, print_dict=False)
        #print('!!!!!!', temp_dict)
        for keyi, valuei in sorted(temp_dict[key].items()):
            #print('&& %s %s' % (keyi, valuei))
            try:
                setattr(obj, keyi, valuei)
            except AttributeError:
                self.log.error('cant set %r as %s; is this a property?' % (keyi, valuei))
                continue
        return obj

    def _load_dict(self, model, h5_file, key, value, custom_types, self_obj, nlevels, print_dict=True):
        #print('%svalue = %s' % (nlevels*'  ', value))
        keys = value.keys()
        new_dict = {}
        sub_h5 = value
        for keyi in sorted(keys):
            self._load_value(new_dict, sub_h5, keyi, custom_types, self_obj, nlevels+2)
        if print_dict:
            print('%s%s' % (nlevels*'  ', str(new_dict)))
        model[key] = new_dict


def _cast(h5_result_attr):
    """converts the h5py type back into the actual type"""
    if h5_result_attr is None:
        return None

    if len(h5_result_attr.shape) == 0:
        return np.array(h5_result_attr).tolist()
        #raise NotImplementedError(h5_result_attr.dtype)
    return np.array(h5_result_attr)
