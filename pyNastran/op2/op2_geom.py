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
from pyNastran.op2.op2_vectorized import OP2_Vectorized

class OP2Geom(BDF,
              GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT, DIT, DYNAMICS,
              OP2):
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param debug: enables the debug log and sets the debug in the logger (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        :param debug_file: sets the filename that will be written to (default=None -> no debug)
        """
        assert make_geom == True, make_geom

        BDF.__init__(self, debug=debug, log=log)
        GEOM1.__init__(self)
        GEOM2.__init__(self)
        GEOM3.__init__(self)
        GEOM4.__init__(self)

        EPT.__init__(self)
        MPT.__init__(self)
        DIT.__init__(self)
        DYNAMICS.__init__(self)

        OP2.__init__(self, debug, log=log, debug_file=debug_file)
        self.make_geom = make_geom

    def _get_table_mapper(self):
        table_mapper = OP2._get_table_mapper(self)

        # geometry
        table_mapper['GEOM1'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4'] = [self._read_geom4_4, self._read_geom4_4]

        # superelements
        table_mapper['GEOM1S'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2S'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3S'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4S'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper['GEOM1N'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2N'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3N'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4N'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper['GEOM1OLD'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2OLD'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3OLD'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4OLD'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper['EPT']    = [self._read_ept_4, self._read_ept_4]
        table_mapper['EPTS']   = [self._read_ept_4, self._read_ept_4]
        table_mapper['EPTOLD'] = [self._read_ept_4, self._read_ept_4]

        table_mapper['MPT']  = [self._read_mpt_4, self._read_mpt_4]
        table_mapper['MPTS'] = [self._read_mpt_4, self._read_mpt_4]

        table_mapper['DYNAMIC']  = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper['DYNAMICS'] = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper['DIT']      = [self._read_dit_4, self._read_dit_4]   # table objects (e.g. TABLED1)


class OP2Geom_Vectorized(OP2Geom, OP2_Vectorized):
    def __init__(self, debug=False, log=None, debug_file=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param debug: enables the debug log and sets the debug in the logger (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        :param debug_file: sets the filename that will be written to (default=None -> no debug)
        """
        # make_geom=False, debug=True, log=None, debug_file=None
        OP2_Vectorized.__init__(self, debug=debug, log=log, debug_file=debug_file)
        OP2Geom.__init__(self, debug=debug, log=log, debug_file=debug_file)
        self.make_geom = True
