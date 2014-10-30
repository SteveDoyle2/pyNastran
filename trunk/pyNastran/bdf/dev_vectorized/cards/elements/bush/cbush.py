from six.moves import zip
from numpy import array, dot, arange, zeros, unique, searchsorted, full, nan, asarray, isnan
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, string_or_blank, blank)
from pyNastran.bdf.cards.elements.bars import CBAROR

from pyNastran.bdf.fieldWriter import set_string8_blank_if_default
from pyNastran.bdf.fieldWriter16 import set_string16_blank_if_default

from pyNastran.bdf.dev_vectorized.cards.elements.element import Element


class CBUSH(Element):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR  | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR  | EID | PID | GA  | GB  | G0  |     |     | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |  CBAR | 2     |  39 | 7     | 6     |  105   |       |       |  GGG  |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |       |       | 513 | 0.0+0 | 0.0+0 |    -9. | 0.0+0 | 0.0+0 |   -9. |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    """
    type = 'CBUSH'
    def __init__(self, model):
        """
        Defines the CBUSH object.

        :param self: the CBUSH object
        :param model: the BDF object
        """
        Element.__init__(self, model)

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.element_id = zeros(ncards, 'int32')
        self.property_id = zeros(ncards, 'int32')
        self.node_ids = zeros((ncards, 2), 'int32')
        self.is_g0 = zeros(ncards, 'bool')
        self.g0 = full(ncards, nan, 'int32')
        self.x = full((ncards, 3), nan, float_fmt)
        self.cid = full(ncards, nan, 'int32')
        self.s = full(ncards, nan, float_fmt)
        self.ocid = full(ncards, nan, 'int32')
        self.si = full((ncards, 3), nan, float_fmt)

    def build(self):
        """
        :param self: the CBUSH object
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards

        if ncards:
            float_fmt = self.model.float

            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.node_ids = zeros((ncards, 2), 'int32')

            self.is_g0 = zeros(ncards, 'bool')
            self.g0 = full(ncards, nan, 'int32')
            self.x = full((ncards, 3), nan, float_fmt)

            self.cid = full(ncards, nan, 'int32')
            self.s = full(ncards, nan, float_fmt)
            self.ocid = full(ncards, nan, 'int32')
            self.si = full((ncards, 3), nan, float_fmt)

            for i, card in enumerate(cards):
                eid = integer(card, 1, 'element_id')
                self.element_id[i] = eid
                self.property_id[i] = integer_or_blank(card, 2, 'property_id', eid)
                self.node_ids[i] = [integer(card, 3, 'GA'),
                                    integer(card, 4, 'GB')]

                #---------------------------------------------------------
                # x / g0
                field5 = integer_double_or_blank(card, 5, 'x1_g0')
                if isinstance(field5, int):
                    self.is_g0[i] = True
                    self.g0[i] = field5
                elif isinstance(field5, float):
                    self.is_g0[i] = False
                    x = array([field5,
                               double_or_blank(card, 6, 'x2', x2_default),
                               double_or_blank(card, 7, 'x3', x3_default)], dtype='float64')
                    self.x[i, :] = x
                    if norm(x) == 0.0:
                        msg = 'G0 vector defining plane 1 is not defined on %s %s.\n' % (self.type, eid)
                        msg += 'G0 = %s\n' % field5
                        msg += 'X  = %s\n' % x
                        msg += '%s' % card
                        raise RuntimeError(msg)
                #else:
                    #msg = ('field5 on %s (G0/X1) is the wrong type...eid=%s field5=%s '
                           #'type=%s' % (self.type, eid, field5, type(field5)))
                    #raise RuntimeError(msg)

                #---------------------------------------------------------
                #: Element coordinate system identification. A 0 means the basic
                #: coordinate system. If CID is blank, then the element coordinate
                #: system is determined from GO or Xi.
                #: (default=blank=element-based)
                cid = integer_or_blank(card, 8, 'cid')
                if cid is not None:
                    self.cid[i] = cid
                #: Location of spring damper (0 <= s <= 1.0)
                self.s[i] = double_or_blank(card, 9, 's', 0.5)
                #: Coordinate system identification of spring-damper offset. See
                #: Remark 9. (Integer > -1; Default = -1, which means the offset
                #: point lies on the line between GA and GB
                self.ocid[i] = integer_or_blank(card, 10, 'ocid', -1)
                #: Components of spring-damper offset in the OCID coordinate system
                #: if OCID > 0.
                self.si[i, :] = [double_or_blank(card, 11, 's1'),
                           double_or_blank(card, 12, 's2'),
                           double_or_blank(card, 13, 's3')]
                assert len(card) <= 14, 'len(CBUSH card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]

            self.is_g0 = self.is_g0[i]
            self.g0 = self.g0[i]
            self.x = self.x[i, :]

            self.cid = self.cid[i]
            self.s = self.s[i]
            self.ocid = self.ocid[i]
            self.si = self.si[i]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CBUSH IDs...')
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')


    #=========================================================================
    def write_bdf(self, f, size=8, element_ids=None):
        if self.n:
            if element_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.element_id)

            for (eid, pid, n, is_g0, g0, x, cid, s, ocid, si) in zip(
                self.element_id[i], self.property_id[i], self.node_ids[i],
                self.is_g0[i], self.g0[i], self.x[i], self.cid[i], self.s[i], self.ocid[i], self.si[i]):

                x1 = g0 if is_g0 else x[0]
                x2 = 0 if is_g0 else x[1]
                x3 = 0 if is_g0 else x[2]
                if isnan(x1):
                    x1 = ''
                if isnan(x2):
                    x2 = ''
                if isnan(x3):
                    x3 = ''

                si1, si2, si3 = si
                if isnan(si1):
                    si1 = ''
                if isnan(si2):
                    si2 = ''
                if isnan(si3):
                    si3 = ''

                ocid = set_blank_if_default(ocid, -1)
                s = set_blank_if_default(s, 0.5)
                card = ['CBUSH', eid, pid, n[0], n[1], x1, x2, x3,
                        cid, s, ocid, si1, si2, si3 ]
                f.write(print_card(card))

    def slice_by_index(self, i):
        i = asarray(i)
        obj = CBUSH(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        obj.is_g0 = self.is_g0[i]
        obj.g0 = self.g0[i]
        obj.x = self.x[i, :]
        obj.s = self.s[i, :]
        obj.si = self.si[i]
        obj.cid = self.cid[i]
        obj.ocid = self.ocid[i]
        return obj
