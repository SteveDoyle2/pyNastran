#pylint: disable=C0103,W0201,W0223,R0901,R0902,R0904
"""
Main OP2 class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
#import warnings
from numpy import array

#from pyNastran.utils import is_binary
#from pyNastran.utils.gui_io import load_file_dialog
from pyNastran.op2.dev_explicit import OP2

class OP2_Vectorized(OP2):

    def __init__(self, make_geom=False, save_skipped_cards=False,
                 debug=True, log=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param save_skipped_cards: creates the skippedCards.out file (default=False)
        :param debug: prints data about how the OP2 was parsed (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        """
        OP2.__init__(self, make_geom=False, save_skipped_cards=False,
                 debug=True, log=None, debug=debug, log=log)
        self.ask = None

    def set_as_vectorized(ask=False):
        """
        Enables vectorization (currently subject to change and break frequently)

        The code will degenerate to dictionary based results when
        a result does not support vectorization.

        Vectorization is always True here.
        :param ask:  Do you want to see a GUI of result types.

        +--------+---------------+---------+------------+
        | Case # | Vectorization |  Ask    | Read Modes |
        +========+===============+=========+============+
        |    1   | True          |  True   |  1, 2      |
        |    2   | True          |  False  |  1, 2      |
        |    3   | False         |  True   |  1, 2      |
        |    4   | False         |  False  |  0         | <------ < v0.8
        +--------+---------------+---------+------------|

        Definitions
        ===========
          Vectorization - A storage structure that allows for faster read/access
                          speeds and better memory usage, but comes with a more
                          difficult to use data structure.

                          It limits the node IDs to all be integers (e.g. element
                          centroid).  Composite plate elements (even for 1 type)
                          with an inconsistent number of layers will have a more
                          difficult data structure.  Models with solid elements of
                          mixed type will also be more complicated (or potentially
                          split up).
          Scanning   - a quick check used to figure out how many results to process
                      that takes almost no time
          Reading    - process the op2 data
          Build      - call the __init__ on a results object (e.g. DisplacementObject)
          Start Over - Go to the start of the op2 file
          Ask        - launch a GUI dialog to let the user click which results to load

        Read Mode Definitions
        =====================
          0.   The default OP2 dictionary based-approach with no asking GUI
          1.   The first read of a result to get the shape of the data
               (and result types if vectorization=True or result types
               if vectorization=False)
          2.   The second read of a result to get the results

        Cases
        ======
          1.   Scan the block to get the size, build the object (read_mode=1),
               ask the user, start over, fill the data (read_mode=2).
               Degenerate to read_mode=0 when read_mode=2 cannot be used based
               upon the value of ask.
          2.   Scan the block to get the size, build the object (read_mode=1),
               start over, read the OP2 to fill the objects (read_mode=2).
          3.   Scan the block to get the object types (read_mode=1), ask the user,
               build the object & fill it (read_mode=2)
          4.   Read the block to get the size, build the object & fill it (read_mode=0)

        """
        self.ask = ask

    def read_op2(self, op2_filename=None):
        """
        Starts the OP2 file reading
        """
        assert self.ask in [True, False], self.ask
        if self.ask:
            self.read_mode = 0
            self._close_op2 = False

            # get GUI object names, build objects, but don't read data
            self.read_mode = 0
            self.OP2.read_op2(op2_filename)

            # TODO: stuff to figure out objects
            # TODO: stuff to show gui of table names
            # TODO: clear out objects the user doesn't want
            self.read_mode = 1
            self._close_op2 = True
            self.OP2.read_op2(self.op2_filename)
        else:
            self.read_mode = 2
            self.OP2.read_op2(self.op2_filename)
        self.op2.close()
        self.skippedCardsFile.close()

if __name__ == '__main__':
    import pyNastran
    pkg_path = pyNastran.__file__[0]

    # we don't want the variable name to get picked up by the class
    _op2_filename = os.path.join(pyNastran, '..', 'models', 'solid_bending', 'solid_bending.op2')

    model = OP2_Vectorized()
    model.read_op2(_op2_filename)