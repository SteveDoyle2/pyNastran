"""
tests rainflow.py
"""
from __future__ import print_function
import os
import unittest

import numpy as np
from pyNastran.utils import remove_files
from pyNastran.applications.rainflow import (
    rainflow_from_arrays, rainflow_from_csv, write_min_max_stress_dict)

class TestRainflow(unittest.TestCase):
    def test_rainflow_csv(self):
        """tests rainflow functions on a contrived problem"""
        input_csv = 'test.csv'
        n = 700
        n1 = n // 3
        casenames = {
            ('normal', 0, n1),
            ('impulse', n1, n - 1),
        }

        x = np.linspace(0., 3.14*5, num=n)
        y = np.sin(x) * np.cos(201 * x)
        z = np.sin(x) * np.cos(201 * x) * np.tan(x)
        A = np.vstack([y, z])
        np.savetxt(input_csv, A.T, delimiter=',')
        features = {
            0 : 'fillet',
            1 : 'groove',
        }
        rainflow_from_csv(input_csv, casenames, features,
                          write_csvs=True, delimiter=',',
                          xmax=None, legend_alpha=1.0, debug=False)

        stress_filename = 'stress.csv'
        eids = list(features.keys())
        min_max_stresses_dict = rainflow_from_arrays(A.T, features=None, casenames=None,
                                                     nfeatures_max=None,
                                                     filter_stresses=True)
        write_min_max_stress_dict(stress_filename, min_max_stresses_dict, eids)

        min_max_stresses_dict = rainflow_from_arrays(A.T, features=features, casenames=casenames,
                                                     nfeatures_max=None,
                                                     filter_stresses=True)
        write_min_max_stress_dict(stress_filename, min_max_stresses_dict, eids)


        failed = remove_files(
            'test.csv',
            'feature0_impulse_fillet.csv', 'feature0_normal_fillet.csv',
            'feature1_impulse_groove.csv', 'feature1_normal_groove.csv',
            'fillet.png', 'groove.png',
            'stress.csv',
        )
        if failed:  # pragma: no cover
            print(failed)


if __name__ == '__main__':
    unittest.main()


