# coding: utf-8
"""tests for dict_to_h5py"""
import os
import unittest
import numpy as np
from cpylog import get_logger

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
        log = get_logger(log=None, level='warning', encoding='utf-8')
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
        export_obj_to_hdf5('test.h5', obj, log=log)
        #export_obj_to_hdf5_file('test.h5', ap, log)
        new_dict = load_obj_from_hdf5('test.h5', custom_types, log=log)
        #print('new_dict[ap]', new_dict['ap'])
        #print('variable_type_list', new_dict['variable_type_list'])
        #print('variable_type_tuple', new_dict['variable_type_tuple'])
        export_obj_to_hdf5('test_new.h5', new_dict, log=log)


        os.remove('test.h5')
        os.remove('test_new.h5')

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
        assert isinstance(new_dict['str_key_unicode_value'], str)
        assert isinstance(new_dict[u'helló wörld from two'], bytes), type(new_dict[u'helló wörld from two'])
        assert new_dict['None'] is None, new_dict['None']
        assert np.isnan(new_dict['nan']), new_dict['nan']
        #str_key_unicode_value

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
