"""Interface to various numpy utilities"""
from __future__ import absolute_import

from six import PY2
import numpy as np

from pyNastran.femutils.io import loadtxt_nice, savetxt_nice

from pyNastran.femutils.nan import (
    isfinite, isfinite_and_greater_than, isfinite_and_nonzero)

from pyNastran.femutils.utils import (
    cross2d, unique2d, unique_rows, augmented_identity)

if PY2:
    integer_types = (int, long, np.int32, np.int64)
    integer_string_types = (int, long, np.int32, np.int64, str, unicode)
    integer_float_types = (int, long, np.int32, np.int64, float, np.float32)
else:
    integer_types = (int, np.int32, np.int64)
    integer_string_types = (int, np.int32, np.int64, bytes, str)
    integer_float_types = (int, np.int32, np.int64, float, np.float32)
float_types = (float, np.float32, np.float64)
