import numpy as np

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, components_or_blank,
    # integer_or_blank, double_or_blank, string_or_blank
)

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

#class DAREA:
    #type = 'DAREA'
    #def __init__(self, model):
        #self.n = 0
    #def allocate(self, card_count):
        #pass
    #def build(self):
        #pass
    #def parse(self, card_obj, icard=0, comment=''):
        #card = (1, 2., icard)
        #return card, comment
    #def add_card(self, card, comment=''):
        #return

class DAREA(VectorizedCard):
    """
    Defines scale (area) factors for static and dynamic loads. In dynamic
    analysis, DAREA is used in conjunction with ACSRCE, RLOADi and TLOADi
    entries.

    RLOAD1 -> DAREA by SID

    +-------+-----+----+----+-----+----+----+------+
    |   1   |  2  | 3  |  4 |  5  | 6  |  7 |  8   |
    +=======+=====+====+====+=====+====+====+======+
    | DAREA | SID | P1 | C1 | A1  | P2 | C2 | A2   |
    +-------+-----+----+----+-----+----+----+------+
    | DAREA | 3   | 6  | 2  | 8.2 | 15 | 1  | 10.1 |
    +-------+-----+----+----+-----+----+----+------+

    referenced by:
      - TLOAD1, TLOAD2, RLOAD1, RLOAD2 and ??? entries
    """
    type = 'DAREA'
    def __init__(self, model):
        """
        Defines the DAREA object.

        Parameters
        ----------
        model : BDF
           the BDF object

        .. todo:: collapse loads
        """
        VectorizedCard.__init__(self, model)

    def parse(self, card, icard=0, comment=''):
        noffset = 3 * icard
        sid = integer(card, 1, 'sid')
        node = integer(card, 2 + noffset, 'node_p')
        component = int(components_or_blank(card, 3 + noffset, 'c', '0'))
        scale = double(card, 4 + noffset, 'scale')
        data = (sid, node, component, scale)
        return data, comment

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt

            #: Identification number of DELAY entry. (Integer > 0)
            self.sid = np.zeros(ncards, 'int32')
            #: Grid, extra, or scalar point identification number. (Integer > 0)
            self.node_id = np.zeros(ncards, 'int32')
            #: Component number. (Integers 1 through 6 for grid points; zero or blank for extra
            #: or scalar points)
            self.component = np.zeros(ncards, 'int32')
            #: Scale (area) factor. (Real)
            self.scale = np.zeros(ncards, float_fmt)

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

    def add_card(self, data, comment=''):
        i = self.i
        print(data)
        (sid, node, component, scale) = data

        #sid = integer(card, 1, 'sid')
        #node = integer(card, 2, 'node')
        #component = integer(card, 3, 'component')
        #delay = double_or_blank(card, 4, 'delay')

        assert component in [0, 1, 2, 3, 4, 5, 6], 'component=%r' % component

        #assert len(card) <= 13, 'len(TLOAD2 card) = %i\ncard=%s' % (len(card), card)

        self.sid[i] = sid
        self.node_id[i] = node
        self.component[i] = component
        self.scale[i] = scale
        self.i += 1

    def build(self):
        if self.n > 1:
            i = self.sid.argsort()
            self.sid = self.sid[i]
            self.node_id = self.node_id[i]
            self.component = self.component[i]
            self.scale = self.scale[i]

    def update(self, maps):
        """
        maps = {
            'node_id' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            ## TODO: support DAREA id
            nid_map = maps['node']
            for i, nid in enumerate(self.node_id):
                print(self.print_card(i))
                self.node_ids[i] = nid_map[nid]

    #def get_load_ids(self):
        #return unique(self.load_id)

    #@property
    #def load_id(self):
        #return self.sid

    def get_darea_index_by_darea_id(self, sid):
        if sid is None:
            return np.arange(self.n)
        ##msg = ''
        assert isinstance(sid, int), sid
        return np.where(self.sid == sid)[0]

    def write_card(self, bdf_file, size=8, is_double=False, sid=None):
        if self.n:
            i = self.get_darea_index_by_darea_id(sid)
            for (sid, nid, component, scale) in zip(
                    self.sid[i], self.node_id[i], self.component[i], self.scale[i],):

                list_fields = ['DAREA', sid, nid, component, scale]
                if size == 8:
                    bdf_file.write(print_card_8(list_fields))
                else:
                    bdf_file.write(print_card_16(list_fields))
