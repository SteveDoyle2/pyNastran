from pathlib import PurePath
from .op2_geom import OP2, OP2Geom
from pyNastran.dev.op2_vectorized3.op2_interface.h5_pytables.h5_results import read_h5_result
from pyNastran.dev.op2_vectorized3.op2_interface.h5_pytables.h5_results import read_h5_geometry_result

class Results(Results):
    def read_h5(self, h5_filename: Union[str, PurePath], combine=None):
        read_h5_result(self, h5_filename, root_path='/')

class ResultsGeom(Results, OP2Geom):
    def read_h5(self, h5_filename: Union[str, PurePath], combine=None):
        """TODO: should support geometry"""
        read_h5_geometry_result(self, h5_filename, root_path='/')
