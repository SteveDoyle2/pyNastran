# coding: utf-8
"""tests for dict_to_h5py"""
import unittest
from six import text_type, binary_type
import numpy as np

try:
    import h5py  # pylint: disable=unused-import
    IS_H5PY = True
except ImportError:  # pragma: no cover
    IS_H5PY = False

if IS_H5PY:
    from pyNastran.utils.dict_to_h5py import load_obj_from_hdf5, export_obj_to_hdf5
from pyNastran.bdf.bdf import BDF


class TestDictToH5(unittest.TestCase):

    @unittest.skipIf(not IS_H5PY, "No h5py")
    def test_dict_to_h5py(self):
        model = BDF()
        obj = {
            'bdf' : model,
            'key' : 'value',
            #1 : 2,
            #3.0 : 4.0,
            'int_list' : [1, 2, 3],
            'float_list' : [1.1, 2.2, 3.3],
            'mydict' : {'one' : 1},
            'five' : np.zeros(5),
            'None' : None,
            'nan' : np.nan,
            'variable_type_list' : [1, 2.2, b'3.14s', u'4.4u'],
            'variable_type_tuple' : (1, 2.2, b'3.14s', u'4.4u'),
            'str_key_unicode_value' : u'helló wörld from two',
            'helló wörld from two' : b'cat',
        }

        custom_types = {
            'BDF' : BDF,
        }
        export_obj_to_hdf5('test.h5', obj, log=None)
        #export_obj_to_hdf5_file('test.h5', ap, log)
        new_dict = load_obj_from_hdf5('test.h5', custom_types, log=None)
        #print('new_dict[ap]', new_dict['ap'])
        #print('variable_type_list', new_dict['variable_type_list'])
        #print('variable_type_tuple', new_dict['variable_type_tuple'])
        export_obj_to_hdf5('test_new.h5', new_dict, log=None)


        #os.remove('test.h5')
        #os.remove('test_new.h5')

        #obj = {
            #'key' : 'value',
            ##1 : 2,
            ##3.0 : 4.0,
            #'mydict' : {'one' : 1},
            #'five' : np.zeros(5),
            #'ap' : analytic_predictor,
            #'None' : None,
            #'nan' : np.nan,
            #'variable_type_list' : [1, 2.2, '3.14s'],
            #'variable_type_tuple' : (1, 2.2, '3.14s'),
            #'breaking my code' : u'helló wörld from two',
            #'helló wörld from two' : 'cat',
        #}

        #print('new_dict[ap]', new_dict['ap'])
        assert isinstance(new_dict['variable_type_list'], list)
        assert isinstance(new_dict['variable_type_tuple'], tuple)
        assert isinstance(new_dict['five'], np.ndarray)
        assert len(new_dict['five']) == 5
        assert isinstance(new_dict['str_key_unicode_value'], text_type)
        assert isinstance(new_dict[u'helló wörld from two'], binary_type), type(new_dict[u'helló wörld from two'])
        assert new_dict['None'] is None, new_dict['None']
        assert np.isnan(new_dict['nan']), new_dict['nan']
        #str_key_unicode_value

if __name__ == '__main__':
    unittest.main()
