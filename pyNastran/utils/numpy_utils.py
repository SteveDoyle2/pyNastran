"""Interface to various numpy utilities"""
from __future__ import absolute_import
from .numpy_functions.io import loadtxt_nice, savetxt_nice

from .numpy_functions.matrix3d import (
    norm2d, axes_stack, normalize_vector2d,
    dot3d, transpose3d, triple, triple_transpose,)

from .numpy_functions.nan import (
    isfinite, isfinite_and_greater_than, isfinite_and_nonzero)

from .numpy_functions.numpy_utils import (
    cross2d, unique2d, unique_rows, augmented_identity,
    perpendicular_vector, perpendicular_vector2d,)
