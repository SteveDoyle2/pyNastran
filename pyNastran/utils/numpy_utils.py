"""Interface to various numpy utilities"""
from __future__ import absolute_import
from pyNastran.femutils.io import loadtxt_nice, savetxt_nice

from pyNastran.femutils.matrix3d import (
    norm2d, axes_stack, normalize_vector2d,
    dot3d, transpose3d, triple, triple_transpose,)

from pyNastran.femutils.nan import (
    isfinite, isfinite_and_greater_than, isfinite_and_nonzero)

from pyNastran.femutils.utils import (
    cross2d, unique2d, unique_rows, augmented_identity,
    perpendicular_vector, perpendicular_vector2d,)
