from pyNastran.op2.tables.geom.geom1 import GEOM1
from pyNastran.op2.tables.geom.geom2 import GEOM2
from pyNastran.op2.tables.geom.geom3 import GEOM3
from pyNastran.op2.tables.geom.geom4 import GEOM4

from pyNastran.op2.tables.geom.ept import EPT
from pyNastran.op2.tables.geom.mpt import MPT

from pyNastran.op2.tables.geom.dit import DIT
from pyNastran.op2.tables.geom.dynamics import DYNAMICS

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

def read_op2_geom(op2_filename=None, combine=True,
             log=None, debug=True, debug_file=None, build_dataframe=False,
             skip_undefined_matrices=True, mode='msc', encoding=None):
    """
    Creates the OP2 object without calling the OP2 class.

    Parameters
    ----------
    op2_filename : str (default=None -> popup)
        the op2_filename
    combine : bool; default=True
        True : objects are isubcase based
        False : objects are (isubcase, subtitle) based;
                will be used for superelements regardless of the option
    build_dataframe : bool; default=False
        builds a pandas DataFrame for op2 objects
    skip_undefined_matrices : bool; default=False
         True : prevents matrix reading crashes
    debug : bool; default=False
        enables the debug log and sets the debug in the logger
    log : Log()
        a logging object to write debug messages to
     (.. seealso:: import logging)
    debug_file : str; default=None (No debug)
        sets the filename that will be written to
    encoding : str
        the unicode encoding (default=None; system default)

    Returns
    -------
    model : OP2()
        an OP2 object

    .. todo:: creates the OP2 object without all the read methods

    .. note :: this method will change in order to return an object that
               does not have so many methods
    """
    model = OP2Geom(log=log, debug=debug, debug_file=debug_file, mode=mode)
    model.read_op2(op2_filename=op2_filename, build_dataframe=build_dataframe,
                   skip_undefined_matrices=skip_undefined_matrices, combine=combine,
                   encoding=encoding)
    return model

class OP2Geom(OP2, BDF,
              GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT, DIT, DYNAMICS):
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None, mode='msc'):
        """
        Initializes the OP2 object

        Parameters
        ----------
        make_geom : bool; default=False
            reads the BDF tables
        debug : bool; default=False
            enables the debug log and sets the debug in the logger
        log: log()
            a logging object to write debug messages to
            (.. seealso:: import logging)
        debug_file : default=None -> no debug
            sets the filename that will be written to
        mode : str; default='msc'
            {msc, nx}
        """
        # make_geom=False, debug=True, log=None, debug_file=None

        BDF.__init__(self, debug=debug, log=log)
        GEOM1.__init__(self)
        GEOM2.__init__(self)
        GEOM3.__init__(self)
        GEOM4.__init__(self)

        EPT.__init__(self)
        MPT.__init__(self)
        DIT.__init__(self)
        DYNAMICS.__init__(self)

        OP2.__init__(self, debug, log=log, debug_file=debug_file, mode=mode)
        self.make_geom = True

    def _get_table_mapper(self):
        table_mapper = OP2._get_table_mapper(self)

        # geometry
        table_mapper[b'GEOM1'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4'] = [self._read_geom4_4, self._read_geom4_4]

        # superelements
        table_mapper[b'GEOM1S'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2S'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3S'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4S'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper[b'GEOM1N'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2N'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3N'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4N'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper[b'GEOM1OLD'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2OLD'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3OLD'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4OLD'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper[b'EPT'] = [self._read_ept_4, self._read_ept_4]
        table_mapper[b'EPTS'] = [self._read_ept_4, self._read_ept_4]
        table_mapper[b'EPTOLD'] = [self._read_ept_4, self._read_ept_4]

        table_mapper[b'MPT'] = [self._read_mpt_4, self._read_mpt_4]
        table_mapper[b'MPTS'] = [self._read_mpt_4, self._read_mpt_4]

        table_mapper[b'DYNAMIC'] = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper[b'DYNAMICS'] = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper[b'DIT'] = [self._read_dit_4, self._read_dit_4]   # table objects (e.g. TABLED1)
        return table_mapper
