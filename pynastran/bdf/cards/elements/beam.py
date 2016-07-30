# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from six import string_types
import numpy as np
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.elements.bars import CBAR, LineElement
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_string_or_blank)
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class CBEAM(CBAR):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    | CBEAM | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    | CBEAM | EID | PID | GA  | GB  | G0  |     |     | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    """
    type = 'CBEAM'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', #5:'x_g0', 6:'g1', 7:'g2',
        #8:'offt',
        9:'pa', 10:'pb',
        17:'sa', 18:'sb',
    }

    def _update_field_helper(self, n, value):
        if n == 11:
            self.wa[0] = value
        elif n == 12:
            self.wa[1] = value
        elif n == 13:
            self.wa[2] = value

        elif n == 14:
            self.wb[0] = value
        elif n == 15:
            self.wb[1] = value
        elif n == 16:
            self.wb[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:  # offt
                    msg = 'Field %r=%r is an invalid %s entry or is unsupported.' % (
                        n, value, self.type)
                    raise KeyError(msg)
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    msg = 'Field %r=%r is an invalid %s entry or is unsupported.' % (
                        n, value, self.type)
                    raise KeyError(msg)

    def __init__(self, eid, pid, ga, gb, x, g0, is_offt, offt, bit,
                 pa, pb, wa, wb, sa, sb, comment=''):
        LineElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = ga
        self.gb = gb
        self.x = x
        self.g0 = g0
        self.is_offt = is_offt
        self.offt = offt
        self.bit = bit
        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb
        self.sa = sa
        self.sb = sb
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')

        x, g0 = cls._init_x_g0(card, eid)
        is_offt, offt, bit = cls._init_offt_bit(card, eid)# offt doesn't exist in NX nastran
        pa = integer_or_blank(card, 9, 'pa', 0)
        pb = integer_or_blank(card, 10, 'pb', 0)

        wa = np.array([double_or_blank(card, 11, 'w1a', 0.0),
                       double_or_blank(card, 12, 'w2a', 0.0),
                       double_or_blank(card, 13, 'w3a', 0.0)], 'float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], 'float64')

        sa = integer_or_blank(card, 17, 'sa', 0)
        sb = integer_or_blank(card, 18, 'sb', 0)
        assert len(card) <= 19, 'len(CBEAM card) = %i\ncard=%s' % (len(card), card)
        return CBEAM(eid, pid, ga, gb, x, g0, is_offt, offt, bit,
                     pa, pb, wa, wb, sa, sb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, f, comment=''):
        #: .. todo:: verify
        assert len(data) == 2, 'data=%s len(data)=%s' % (data, len(data))
        #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
        #        [f,g0]]
        #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
        #        [f,x1,x2,x3]]

        main, aft = data
        flag = aft[0]
        assert f == flag, 'f=%s flag=%s' % (f, flag)
        if flag == 0:
            # basic cid
            #data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, x1, x2, x3]]
            assert len(aft) == 4, 'f=%s aft=%s len(aft)=%s' % (f, aft, len(aft))
            x1, x2, x3 = aft[1:]
            g0 = None
            x = np.array([x1, x2, x3], dtype='float64')
        elif flag == 1:
            # global cid
            #data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, x1, x2, x3]]
            assert len(aft) == 4, 'f=%s aft=%s len(aft)=%s' % (f, aft, len(aft))
            g0 = None
            x1, x2, x3 = aft[1:]
            x = np.array([x1, x2, x3], dtype='float64')
        elif flag == 2:
            # grid option
            #data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, g0]]
            assert len(aft) == 2, 'f=%s aft=%s len(aft)=%s' % (f, aft, len(aft))
            g0 = data[1][1]
            x = None
        else:
            raise NotImplementedError()

        eid = main[0]
        pid = main[1]
        ga = main[2]
        gb = main[3]
        sa = main[4]
        sb = main[5]

        is_offt = True  #: .. todo:: is this correct???
        #offt = str(data[6]) # GGG
        bit = None # ????
        offt = 'GGG'  #: .. todo:: is this correct???

        pa = main[6]
        pb = main[7]

        wa = np.array([main[8], main[9], main[10]], 'float64')
        wb = np.array([main[11], main[12], main[13]], 'float64')
        return CBEAM(eid, pid, ga, gb, x, g0, is_offt, offt, bit,
                     pa, pb, wa, wb, sa, sb, comment=comment)

    def _validate_input(self):
        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

    def Nodes(self):
        return [self.ga, self.gb]

    @classmethod
    def _init_offt_bit(cls, card, eid):
        """
        offt doesn't exist in NX nastran
        """
        field8 = integer_double_string_or_blank(card, 8, 'field8')
        if isinstance(field8, float):
            is_offt = False
            offt = None
            bit = field8
        elif field8 is None:
            is_offt = True
            offt = 'GGG'  # default
            bit = None
        elif isinstance(field8, string_types):
            is_offt = True
            bit = None
            offt = field8
            msg = 'invalid offt parameter of CBEAM...offt=%s' % offt
            assert offt[0] in ['G', 'B', 'O', 'E'], msg
            assert offt[1] in ['G', 'B', 'O', 'E'], msg
            assert offt[2] in ['G', 'B', 'O', 'E'], msg
        else:
            msg = ('field8 on %s card is not a string(offt) or bit '
                   '(float)...field8=%s\n' % (cls.type, field8))
            raise RuntimeError("Card Instantiation: %s" % msg)
        return is_offt, offt, bit

    def Mid(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Mid()

    def Area(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Area()

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Nsm()

    def get_offt_bit_defaults(self):
        """
        offt doesn't exist in NX nastran
        """
        if self.is_offt:
            field8 = field8 = set_blank_if_default(self.offt, 'GGG')
        else:
            field8 = set_blank_if_default(self.bit, 0.0)
        return field8

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        self.gb = model.Node(self.gb, msg=msg)
        self.gb_ref = self.gb
        self.nodes = model.Nodes([self.ga.nid, self.gb.nid], msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid
        if self.g0:
            g0 = model.nodes[self.g0]
            self.g0_vector = g0.get_position() - self.ga.get_position()
        else:
            self.g0_vector = self.x
        if model.is_nx:
            assert self.offt == 'GGG', 'NX only support offt=GGG; offt=%r' % self.offt

    def safe_cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)

        self.ga_ref = self.ga
        self.gb_ref = self.gb
        try:
            self.pid = model.Property(self.pid, msg=msg)
            self.pid_ref = self.pid
        except KeyError:
            model.log.warning('pid=%s%s' % (self.pid, msg))

        if self.g0:
            try:
                g0 = model.nodes[self.g0]
                self.g0_vector = g0.get_position() - self.ga.get_position()
            except KeyError:
                model.log.warning('Node=%s%s' % (self.g0, msg))
        else:
            self.g0_vector = self.x

    def uncross_reference(self):
        self.pid = self.Pid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        del self.ga_ref, self.gb_ref, self.pid_ref

    def raw_fields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga, gb = self.node_ids
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                       self.pa, self.pb] + list(self.wa) + list(self.wb) + [self.sa, self.sb]
        return list_fields

    def repr_fields(self):
        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)
        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)

        sa = set_blank_if_default(self.sa, 0)
        sb = set_blank_if_default(self.sb, 0)
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.get_offt_bit_defaults()
        ga, gb = self.node_ids
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                       self.pa, self.pb, w1a, w2a, w3a,
                       w1b, w2b, w3b, sa, sb]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)
