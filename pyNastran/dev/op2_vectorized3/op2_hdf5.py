from pathlib import PurePath
import warnings
from typing import Union

try:
    import pandas
    IS_PANDAS = True
except ImportError:  # pragma: no cover
    IS_PANDAS = False


from .op2_geom import OP2, OP2Geom
if IS_PANDAS:
    from pyNastran.dev.op2_vectorized3.op2_interface.h5_pytables.h5_results import read_h5_result
    from pyNastran.dev.op2_vectorized3.op2_interface.h5_pytables.h5_results import read_h5_geometry_result
else:  # pragma: no cover
    warnings.warn('cannot find pandas.  Disabling h5 support')


class Results(OP2):
    def read_h5(self, h5_filename: str | PurePath, combine=None):
        read_h5_result(self, h5_filename, root_path='/')

class ResultsGeom(Results, OP2Geom):
    def read_h5(self, h5_filename: str | PurePath, combine=None):
        """TODO: should support geometry"""
        read_h5_geometry_result(self, h5_filename, root_path='/')
