from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.dev.bdf_vectorized2.cards.elements.bars import init_x_g0
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class BeamElement:
    """base class for CBEAM"""
    card_name = ''
    def __init__(self, model):
        """initializes the BeamElement"""
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='float64')
        self.offt = np.array([], dtype='|U8')
        self.x = np.array([], dtype='float64')
        self.g0 = np.array([], dtype='int32')
        self.pin_flags = np.array([], dtype='int32')
        self.wa_offset = np.array([], dtype='float64')
        self.wb_offset = np.array([], dtype='float64')
        self.sab_warping = np.array([], dtype='int32')
        self._eid = []
        self._pid = []
        self._nids = []
        self._offt = []
        self._x = []
        self._g0 = []
        self._pin_flags = []
        self._wa_offset = []
        self._wb_offset = []
        self._sab_warping = []
        self.comment = defaultdict(str)

    def check_if_current(self, nid, nids):
        """we split this up to reason about it easier"""
        if self.is_current:
            if nid in nids:
                # card exists, so we use that slot
                add_card = False
            else:
                add_card = True
        else:
            add_card = True
        return add_card

    #def get_element_by_eid(self, eid):
        #self.make_current()
        #ieid = np.searchsorted(eid, self.eid)
        #return self[ieid]

    def make_current(self):
        """creates an array of the GRID points"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                self.offt = np.hstack([self.offt, self._offt])
                self.x = np.hstack([self.x, self._x])
                self.g0 = np.hstack([self.g0, self._g0])
                self.pin_flags = np.hstack([self.pin_flags, self._pin_flags])
                self.wa_offset = np.hstack([self.wa_offset, self._wa_offset])
                self.wb_offset = np.hstack([self.wb_offset, self._wb_offset])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.offt = np.array(self._offt, dtype='|U8')
                self.x = np.array(self._x, dtype='float64')
                self.g0 = np.array(self._g0, dtype='int32')
                self.pin_flags = np.array(self._pin_flags, dtype='int32')
                self.wa_offset = np.array(self._wa_offset, dtype='float64')
                self.wb_offset = np.array(self._wb_offset, dtype='float64')
            assert len(self.eid) == len(np.unique(self.eid))

            isort = np.argsort(self.eid)
            self.eid = self.eid[isort]
            self.pid = self.pid[isort]
            self.nids = self.nids[isort, :]
            self.offt = self.offt[isort]
            self.x = self.x[isort, :]
            self.g0 = self.g0[isort]
            self.pin_flags = self.pin_flags[isort, :]
            self.wa_offset = self.wa_offset[isort, :]
            self.wb_offset = self.wb_offset[isort, :]

            self._eid = []
            self._pid = []
            self._nids = []
            self._offt = []
            self._x = []
            self._g0 = []
            self._pin_flags = []
            self._wa_offset = []
            self._wb_offset = []
            self._sab_warping = []
            self.is_current = True

    def cross_reference(self, model: BDF) -> None:
        """does this do anything?"""
        self.make_current()

    def __len__(self):
        """returns the number of elements"""
        return len(self.eid) + len(self._eid)

    def repr_indent(self, indent=''):
        self.make_current()
        neids = len(self.eid)
        if neids == 0:
            return '%s%sv; nelements=%s' % (indent, self.card_name, neids)
        msg = '%s%sv; nelements=%s:\n' % (indent, self.card_name, neids)
        msg += '%s  eid = %s\n' % (indent, self.eid)

        upid = np.unique(self.pid)
        if len(upid) == 1 and upid[0] == 0:
            msg += '%s  upid = %s\n' % (indent, upid)
        else:
            msg += '%s  pid = %s\n' % (indent, self.pid)

        #msg += '  nid =\n%s' % self.nid
        return msg

    def __repr__(self):
        return self.repr_indent('')


class CBEAMv(BeamElement):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEAM | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEAM | EID | PID | GA  | GB  | G0  |     |     | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    offt/bit are MSC specific fields

    """
    card_name = 'CBEAM'

    def add(self, eid, pid, nids, x, g0, offt='GGG', bit=None,
            pin_flags=None, wa=None, wb=None, sa=0, sb=0, comment=''):
        """
        Adds a CBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        x : list[float, float, float]
            Components of orientation vector, from GA, in the displacement
            coordinate system at GA (default), or in the basic coordinate system
        g0 : int
            Alternate method to supply the orientation vector using grid
            point G0. Direction of is from GA to G0. is then transferred
            to End A
        offt : str; default='GGG'
            Offset vector interpretation flag
            None : bit is active
        bit : float; default=None
            Built-in twist of the cross-sectional axes about the beam axis
            at end B relative to end A.
            For beam p-elements ONLY!
            None : offt is active
        pin_flags : list[int, int]; default=None
            None : [0, 0]; don't release the DOFs
            Pin Flag at End A/B.  Releases the specified DOFs
        wa / wb : list[float, float, float]
            Components of offset vectors from the grid points to the end
            points of the axis of the shear center
        sa / sb : int; default=0
            Scalar or grid point identification numbers for the ends A and B,
            respectively. The degrees-of-freedom at these points are the
            warping variables . SA and SB cannot be specified for
            beam p-elements
        comment : str; default=''
            a comment for the card

        offt/bit are MSC specific fields

        """
        if g0 is None:
            g0 = -1
        else:
            x = [np.nan, np.nan, np.nan]
        if pin_flags is None:
            pin_flags = [0, 0]

        self.model.bars.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        self._x.append(x)
        self._g0.append(g0)
        self._offt.append(offt)
        self._wa_offset.append(wa)
        self._wb_offset.append(wb)
        self._sab_warping.append([sa, sb])
        self._pin_flags.append(pin_flags)
        #self._offset.append(wa_offset)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CBEAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')

        x, g0 = init_x_g0(card, eid)
        offt, bit = init_offt_bit(card, eid)# offt doesn't exist in NX nastran
        pin_flag_a = integer_or_blank(card, 9, 'pa', 0)
        pin_flag_b = integer_or_blank(card, 10, 'pb', 0)

        wa = np.array([double_or_blank(card, 11, 'w1a', 0.0),
                       double_or_blank(card, 12, 'w2a', 0.0),
                       double_or_blank(card, 13, 'w3a', 0.0)], 'float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], 'float64')

        sa = integer_or_blank(card, 17, 'sa', 0)
        sb = integer_or_blank(card, 18, 'sb', 0)
        assert len(card) <= 19, 'len(CBEAM card) = %i\ncard=%s' % (len(card), card)
        return self.add(eid, pid, [ga, gb], x, g0, offt, bit,
                        [pin_flag_a, pin_flag_b], wa, wb, sa, sb, comment=comment)

    #def update(self, grid):
        #"""functions like a dictionary"""
        #nid = grid.nid
        #add_card = self.check_if_current(eid, self.eid)
        #if add_card:
            #self.add(nid, grid.xyz, cp=grid.cp, cd=grid.cd,  # add_cquad4
                     #ps=grid.ps, seid=grid.seid, comment=grid.comment)
            #self.is_current = False
        #else:
            #inid = np.where(nid == self.nid)[0]
            #self.nid[inid] = grid.nid
            #self.xyz[inid] = grid.xyz
            #self.cp[inid] = grid.cp
            #self.cd[inid] = grid.cd
            #self.ps[inid] = grid.ps
            #self.seid[inid] = grid.seid
            #self.comment[nid] = comment
            #self.is_current = True  # implicit

    #def __iter__(self):
        #pass
    #def __next__(self):
        #pass
    #def __items__(self):
        #pass
    #def __keys__(self):
        #pass
    #def __values__(self):
        #pass
    #def __getitem__(self, i):
        #"""this works on index"""
        #self.make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass

    @classmethod
    def get_x_g0_defaults(cls, x, g0):
        """
        X and G0 compete for the same fields, so the method exists to
        make it easier to write the card

        Returns
        -------
        x_g0 : varies
            g0 : list[int, None, None]
            x : list[float, float, float]

        """
        if g0 is not None:
            return (g0, None, None)
        else:
            #print('x =', self.x)
            #print('g0 =', self.g0)
            #x1 = set_blank_if_default(self.x[0], 0.0)
            #x2 = set_blank_if_default(self.x[1], 0.0)
            #x3 = set_blank_if_default(self.x[2], 0.0)
            return list(x)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes, x, g0, offt, pin_flags, wa_offset, wb_offset in zip(
            self.eid, self.pid, self.nids, self.x, self.g0, self.offt, self.pin_flags, self.wa_offset, self.wb_offset):
            x1, x2, x3 = self.get_x_g0_defaults(x, g0)
            #pa = set_blank_if_default(self.pa, 0)
            #pb = set_blank_if_default(self.pb, 0)

            #w1a = set_blank_if_default(self.wa[0], 0.0)
            #w2a = set_blank_if_default(self.wa[1], 0.0)
            #w3a = set_blank_if_default(self.wa[2], 0.0)

            #w1b = set_blank_if_default(self.wb[0], 0.0)
            #w2b = set_blank_if_default(self.wb[1], 0.0)
            #w3b = set_blank_if_default(self.wb[2], 0.0)
            ga, gb = nodes
            pin_flag_a, pin_flag_b = pin_flags
            w1a, w2a, w3a = wa_offset
            w1b, w2b, w3b = wb_offset

            #sa = set_blank_if_default(self.sa, 0)
            #sb = set_blank_if_default(self.sb, 0)
            sa = sb = 0
            #(x1, x2, x3) = self.get_x_g0_defaults()

            # offt doesn't exist in NX nastran
            offt = set_blank_if_default(offt, 'GGG')

            list_fields = ['CBEAM', eid, pid, ga, gb, x1, x2, x3, offt,
                           pin_flag_a, pin_flag_b, w1a, w2a, w3a, w1b, w2b, w3b, sa, sb]
            msg += self.comment[eid] + print_card_8(list_fields)
        bdf_file.write(msg)
        return msg


class Beams:
    """Stores CBEAM elements that exist in 3D space"""
    def __init__(self, model):
        self.model = model
        self.cbeam = model.cbeam
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.cbeam):
            self.cbeam.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.cbeam.make_current()

    def __len__(self):
        return len(self.cbeam)

    def repr_indent(self, indent=''):
        msg = '%s<Beams> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CBEAM: %s\n' % (indent, len(self.cbeam))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')

def init_offt_bit(card, eid):
    """offt doesn't exist in NX nastran"""
    field8 = integer_double_string_or_blank(card, 8, 'field8')
    if isinstance(field8, float):
        offt = None
        bit = field8
    elif field8 is None:
        offt = 'GGG'  # default
        bit = None
    elif isinstance(field8, str):
        bit = None
        offt = field8
        msg = 'invalid offt parameter of CBEAM...offt=%s' % offt
        assert offt[0] in ['G', 'B', 'O', 'E'], msg
        assert offt[1] in ['G', 'B', 'O', 'E'], msg
        assert offt[2] in ['G', 'B', 'O', 'E'], msg
    else:
        msg = ('field8 on %s card is not a string(offt) or bit '
               '(float)...field8=%s\n' % (card.field(0), field8))
        raise RuntimeError("Card Instantiation: %s" % msg)
    return offt, bit
