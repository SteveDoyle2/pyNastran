import numpy as np

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

class DELAY(VectorizedCard):
    """
    +-------+-----+-----------+-----+--------+------+-----+--------+-----+
    |   1   |  2  |     3     |  4  |   5    |  6   |  7  |   8    |  9  |
    +=======+=====+===========+=====+========+======+=====+========+=====+
    | DELAY | SID | POINT ID1 | C1  |   T1   | P2   | C2  |   T2   |     |
    +-------+-----+-----------+-----+--------+------+-----+--------+-----+

    referenced by:
      - TLOAD1, TLOAD2 (time domain loads)
      - RLOAD1, TLOAD2 (frequency domain loads)
      - ACSRCE
    """
    type = 'DELAY'
    def __init__(self, model):
        """
        Defines the DELAY object.

        Parameters
        ----------
        model : BDF
           the BDF object

        .. todo:: collapse loads
        """
        VectorizedCard.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt

            #: Identification number of DELAY entry. (Integer > 0)
            self.sid = np.zeros(ncards, 'int32')
            #: Grid, extra, or scalar point identification number. (Integer > 0)
            self.nodes = np.zeros(ncards, 'int32')
            #: Component number. (Integers 1 through 6 for grid points; zero or blank for extra
            #: or scalar points)
            self.component = np.zeros(ncards, 'int32')
            #: Time delay (tau) for designated point Pi and component Ci. (Real)
            self.delay = np.zeros(ncards, float_fmt)

    #def __getitem__(self, i):
        #unique_lid = unique(self.load_id)
        #if len(i):
            #f = FORCE1(self.model)
            #f.load_id = self.load_id[i]
            #f.node_id = self.node_id[i]
            #f.coord_id = self.coord_id[i]
            #f.mag = self.mag[i]
            #f.xyz = self.xyz[i]
            #f.n = len(i)
            #return f
        #raise RuntimeError('len(i) = 0')

    def add_card(self, card, comment=''):
        i = self.i

        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        component = integer(card, 3, 'component')
        delay = double_or_blank(card, 4, 'delay')

        assert component in [0, 1, 2, 3, 4, 5, 6], component

        #assert len(card) <= 13, 'len(TLOAD2 card) = %i\ncard=%s' % (len(card), card)

        self.sid[i] = sid
        self.node_id[i] = node
        self.component[i] = component
        self.delay[i] = delay
        self.i += 1

    def build(self):
        if self.n:
            i = self.sid.argsort()
            self.sid = self.sid[i]
            self.node_id = self.node_id[i]
            self.component = self.component[i]
            self.delay = self.delay[i]

    #def get_load_ids(self):
        #return unique(self.load_id)

    #@property
    #def load_id(self):
        #return self.sid

    def get_delay_index_by_delay_id(self, sid):
        if sid is None:
            return np.arange(self.n)
        #msg = ''
        assert isinstance(sid, int), sid
        return np.where(self.sid == sid)[0]

    def write_card(self, bdf_file, size=8, is_double=False, sid=None):
        if self.n:
            i = self.get_delay_index_by_delay_id(sid)
            for (sid, nid, component, delay) in zip(
                    self.sid[i], self.node_id[i], self.component[i], self.delay[i],):

                list_fields = ['DELAY', sid, nid, component, delay]

                if size == 8:
                    bdf_file.write(print_card_8(list_fields))
                else:
                    bdf_file.write(print_card_16(list_fields))
