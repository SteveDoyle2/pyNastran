#pylint:  disable=C0111
import copy
from struct import unpack, Struct, pack
from collections import defaultdict
from typing import Optional

import numpy as np
import scipy
from scipy.spatial import KDTree

from pyNastran.utils import is_binary_file
from cpylog import __version__ as CPYLOG_VERSION, SimpleLogger
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger


def read_fld(fld_filename: str,
             log: Optional[SimpleLogger]=None, debug: bool=False):
    """
    Reads an STL file

    Parameters
    ----------
    fld_filename : str
        the filename to read

    Returns
    -------
    model : FLD()
       the fld model

    """
    model = FLD(log=log, debug=debug)
    model.read_fld(fld_filename)
    return model


class FLD:
    #model_type = 'fld'
    #is_structured = False
    #is_outward_normals = True

    def __init__(self, log: Optional[SimpleLogger]=None, debug: bool=False):
        """
        Initializes the FLD object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        """
        self.log = get_logger(log, debug)

        #self.nodes = None
        #self.elements = None
        #self.header = ''
        self.infilename = None

    def write_fld(self, fld_filename: str,
                  stop_on_failure: bool=True) -> None:
        """
        Writes an FLD file

        Parameters
        ----------
        fld_filename : str
            the filename to write
        """
        self.log.info(f'---writing FLD...{stop_on_failure!r}---')
        raise RuntimeError('write-fld')

    def read_fld(self, fld_filename: str) -> None:
        """
        Reads an FLD file

        Parameters
        ----------
        fld_filename : str
            the filename to read

        FIELD: [loadcase_name] : [TABLE]
        FIELD LOCK STATE: [NO]
        DUPLICATE_VALUE_OPTION: [0]
        PARAMETERIZE INDEPENDENT DOMAIN: [NO]
        PERSIST INTERPOL: [NO]
        CREATE INTERPOLATION: [NO]
        FALLBACK DEFAULT INTERPOLATOR: [YES]
        INTERPOL [10]
        VALUES OUTSIDE: [0]
        REMOVE DELAUNAY SLIVERS: [NO]
        MAP: [1]
        INDEP VAR: [x] : [Length] ; [in] : [0]
        BOUNDS: [-100.000[ : [YES] : [100.000] : [YES] : [1000] : [207.000]
        INDEP VAR: [y] : [Length] : [in] : [0]
        BOUNDS: [-100.000[ : [YES] : [100.000] : [YES] : [1000] : [-10.000]
        INDEP VAR: [z] : [Length] : [in] : [0]
        BOUNDS: [-100.000[ : [YES] : [100.000] : [YES] : [1000] : [37.000]
        DEP VAR: [pressure] : [Pressure] : [lbf/in^2(psi)] : [0]
        START DATA
        1.11177, 1.50000, 0.05111, 1.0
        1.50000, 1.88823, 0.05111, 2.0
        1.50000, 1.50000, 0.00000, 3.0
        END DATA
        """

        self.infilename = fld_filename
        self.log.info(f'---reading FLD...{self.infilename}---')

        with open(fld_filename, 'r') as infile:
            lines = infile.readlines()

        header_lines = []
        data_lines = []

        for i, line in enumerate(lines):
            line = line.strip()
            if line.startswith('START'):
                break
            header_lines.append(line)

        i += 1
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('END'):
                break
            sline = line.split(',')
            data_lines.append(sline)
            i += 1
        self.xyzp = np.array(data_lines, dtype='float64')

        ind_vars = []
        dep_vars = []
        for line in header_lines:
            if line.startswith('INDEP VAR'):
                sline = line.split(':')
                print(sline)
                indep_, var__, var_, unit_, zero_ = sline
                var = var_.strip('[] ')
                unit = unit_.strip('[] ')
                ind_vars.append((var, unit))
            elif line.startswith('DEP VAR'):
                sline = line.split(':')
                indep_, var__, var_, unit_, zero_ = sline
                var = var_.strip('[] ')
                unit = unit_.strip('[] ')
                ind_vars.append((var, unit))
        for var in ind_vars:
            print(f'I: {var}')
        for var in dep_vars:
            print(f'D: {var}')
