from pathlib import PurePath
import warnings

try:
    import pandas
    IS_PANDAS = True
except ImportError:  # pragma: no cover
    IS_PANDAS = False

try:
    from tables import open_file, Group, Node, File
    IS_TABLES = True
except ImportError:
    print('pytables was not found; no h5 support.  Run ">>> pip install tables"\n'
          'Do you have h5py installed?  That can cause conflicts.')
    IS_TABLES = False


from .op2_geom import OP2, OP2Geom
if IS_PANDAS and IS_TABLES:
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
