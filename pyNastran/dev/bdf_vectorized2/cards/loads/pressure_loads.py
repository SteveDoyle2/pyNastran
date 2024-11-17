"""
defines:
 - PLOADv
 - PLOAD1v
 - PLOAD2v
 - PLOAD4v

"""
from collections import defaultdict
from itertools import count
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string_or_blank,
    integer_string_or_blank, string, fields, components_or_blank)
from pyNastran.bdf.field_writer_8 import (
    print_float_8, print_card_8, set_blank_if_default, set_string8_blank_if_default)
from pyNastran.bdf.field_writer_16 import (
    print_float_16, set_string16_blank_if_default)
from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.base_card import _format_comment
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.dev.bdf_vectorized2.cards.loads.loads import BaseLoad


class PLOADv(BaseLoad):
    """
    Static Pressure Load

    Defines a uniform static pressure load on a triangular or quadrilateral surface
    comprised of surface elements and/or the faces of solid elements.

    +-------+-----+------+----+----+----+----+
    |   1   |  2  |  3   | 4  | 5  | 6  | 7  |
    +=======+=====+======+====+====+====+====+
    | PLOAD | SID |  P   | G1 | G2 | G3 | G4 |
    +-------+-----+------+----+----+----+----+
    | PLOAD |  1  | -4.0 | 16 | 32 | 11 |    |
    +-------+-----+------+----+----+----+----+
    """
    card_name = 'PLOAD'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.pressure = np.array([], dtype='float64')
        self.nids = np.array([], dtype='int32')

        #self._sid = []
        self._pressure = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, sid, pressure, nids, comment=''):
        """
        Creates a PLOAD card, which defines a uniform pressure load on a
        shell/solid face or arbitrarily defined quad/tri face.

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply
        nodes : list[int]
            The nodes that are used to define the normal are defined
            using the same method as the CTRIA3/CQUAD4 normal.
            n = 3 or 4
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._pressure.append(pressure)
        self._nids.append(nids)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a PLOAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'pressure')
        nodes = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', 0),
        ]
        assert len(card) <= 7, 'len(PLOAD card) = %i\ncard=%s' % (len(card), card)
        self.add(sid, pressure, nodes, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, pressure, node_ids in zip(count(), self.sid, self.pressure, self.nids):
            list_fields = ['PLOAD', sid, pressure] + node_ids.tolist()
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.nids = np.vstack([self.nids, self._nids])

                self.pressure = np.hstack([self.pressure, self._pressure])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.pressure = np.array(self._pressure, dtype='float64')

            self._sid = []
            self._nids = []
            self._pressure = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sPLOADv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class PLOAD1v(BaseLoad):
    """
    Applied Load on CBAR, CBEAM or CBEND Elements

    Defines concentrated, uniformly distributed, or linearly distributed
    applied loads to the CBAR or CBEAM elements at user-chosen points
    along the axis. For the CBEND element, only distributed loads over
    an entire length may be defined.

    +--------+-----+------+------+-------+-----+-------+-----+-------+
    |   1    |  2  |  3   |  4   |   5   |  6  |   7   |  8  |   9   |
    +========+=====+======+======+=======+=====+=======+=====+=======+
    | PLOAD1 | SID | EID  | TYPE | SCALE | X1  |  P1   |  X2 |  P2   |
    +--------+-----+------+------+-------+-----+-------+-----+-------+
    | PLOAD1 | 25  | 1065 |  MY  | FRPR  | 0.2 | 2.5E3 | 0.8 | 3.5E3 |
    +--------+-----+------+------+-------+-----+-------+-----+-------+
    """
    card_name = 'PLOAD1'
    valid_types = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                   'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']

    # LE: length-based; FR: fractional; PR:projected
    valid_scales = ['LE', 'FR', 'LEPR', 'FRPR']

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.eid = np.array([], dtype='int32')
        self.load_type = np.array([], dtype='|U4')
        self.scale = np.array([], dtype='|U4')
        self.p12 = np.array([], dtype='float64')
        self.x12 = np.array([], dtype='float64')

        self._sid = []
        self._eid = []
        self._p12 = []
        self._x12 = []
        self._scale = []
        self._load_type = []
        self.comment = defaultdict(str)

    def add(self, sid, eid, load_type, scale, x1, p1, x2=None, p2=None, comment=''):
        """
        Creates a PLOAD1 card, which may be applied to a CBAR/CBEAM

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element to apply the load to
        load_type : str
            type of load that's applied
            valid_types = {FX, FY, FZ, FXE, FYE, FZE,
                           MX, MY, MZ, MXE, MYE, MZE}
        scale : float
            local pressure scaling factor
        x1 / x2 : float / float
            the starting/end position for the load application
            the default for x2 is x1
        p1 / p2 : float / float
            the magnitude of the load at x1 and x2
            the default for p2 is p1
        comment : str; default=''
            a comment for the card

        Point Load       : x1 == x2
        Distributed Load : x1 != x2
        """

        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._eid.append(eid)
        self._scale.append(scale)
        self._load_type.append(load_type)
        self._p12.append((p1, p2))
        self._x12.append((x1, x2))

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a PLOAD1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        load_type = string(card, 3, 'Type ("%s")' % '",  "'.join(self.valid_types))
        scale = string(card, 4, 'scale ("%s")' % '", "'.join(self.valid_scales))
        x1 = double(card, 5, 'x1')
        p1 = double(card, 6, 'p1')
        x2 = double_or_blank(card, 7, 'x2', x1)
        p2 = double_or_blank(card, 8, 'p2', p1)
        assert len(card) <= 9, 'len(PLOAD1 card) = %i\ncard=%s' % (len(card), card)
        self.add(sid, eid, load_type, scale, x1, p1, x2, p2, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, eid, load_typei, scalei, x12, p12 in zip(
            count(), self.sid, self.eid, self.load_type, self.scale, self.x12, self.p12):
            list_fields = ['PLOAD1', sid, eid, load_typei, scalei,
                           x12[0], p12[0], x12[1], p12[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.eid = np.vstack([self.eid, self._eid])

                self.p12 = np.hstack([self.p12, self._p12])
                self.x12 = np.hstack([self.x12, self._x12])
                self.scale = np.hstack([self.scale, self._scale])
                self.load_type = np.hstack([self.load_type, self._load_type])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.eid = np.array(self._eid, dtype='int32')
                self.load_type = np.array(self._load_type, dtype='|U4')
                self.scale = np.array(self._scale, dtype='|U4')
                self.p12 = np.array(self._p12, dtype='float64')
                self.x12 = np.array(self._x12, dtype='float64')

            self._sid = []
            self._eid = []
            self._p12 = []
            self._x12 = []
            self._scale = []
            self._load_type = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sPLOAD1v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class PLOAD2v(BaseLoad):
    """
    +--------+-----+------+------+------+------+------+------+------+
    |    1   |   2 |  3   |  4   |   5  |   6  |   7  |   8  |   9  |
    +========+=====+======+======+======+=============+======+======+
    | PLOAD2 | SID |  P   | EID1 | EID2 | EID3 | EID4 | EID5 | EID6 |
    +--------+-----+------+------+------+------+------+------+------+
    | PLOAD2 | 21  | -3.6 |  4   |  16  |  2   |      |      |      |
    +--------+-----+------+------+------+------+------+------+------+
    | PLOAD2 | SID |  P   | EID1 | THRU | EID2 |      |      |      |
    +--------+-----+------+------+------+------+------+------+------+
    """
    card_name = 'PLOAD2'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.pressure = np.array([], dtype='float64')
        self.eid = np.array([], dtype='int32')

        self._pressure = []
        self._eid = []
        self.comment = defaultdict(str)

    def add(self, sid, pressure, eids, comment=''):
        """
        Creates a PLOAD2 card, which defines an applied load normal to the quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply to the elements
        eids : list[int]
            the elements to apply pressure to
            n < 6 or a continouus monotonic list of elements (e.g., [1, 2, ..., 1000])
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._pressure.append(pressure)
        self._eid.append(eids)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a PLOAD2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'p')

        if integer_string_or_blank(card, 4, 'THRU') == 'THRU':
            e1 = integer(card, 3, 'Element1')
            e2 = integer(card, 5, 'Element1')
            eids = [i for i in range(e1, e2 + 1)]
            assert len(card) == 6, 'len(PLOAD2 card) = %i\ncard=%s' % (len(card), card)
        else:
            eids = fields(integer, card, 'eid', i=3, j=len(card))
        for eid in eids:
            self.add(sid, pressure, eid, comment=comment)
            comment = ''

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, pressure, eid in zip(count(), self.sid, self.pressure, self.eid):
            list_fields = ['PLOAD2', sid, pressure, eid]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.eid = np.vstack([self.eid, self._eid])

                self.pressure = np.hstack([self.pressure, self._pressure])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.eid = np.array(self._eid, dtype='int32')
                self.pressure = np.array(self._pressure, dtype='float64')

            self._sid = []
            self._eid = []
            self._pressure = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sPLOAD2v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class PLOAD4v(BaseLoad):
    """
    ``Solid Format``

    Defines a pressure load on a face of a CHEXA, CPENTA, or CTETRA element.

    +--------+-----+-----+----+----+------+------+------+-------+
    |   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    +========+=====+=====+====+====+======+======+======+=======+
    | PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  |  G1  | G3/G4 |
    +--------+-----+-----+----+----+------+------+------+-------+
    |        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    +--------+-----+-----+----+----+------+------+------+-------+

    ``Shell Format``

    Defines a pressure load on a face of a CTRIA3, CTRIA6, CTRIAR,
    CQUAD4, CQUAD8, or CQUADR element.

    +--------+-----+-----+----+----+------+------+------+-------+
    |   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    +========+=====+=====+====+====+======+======+======+=======+
    | PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  | THRU | EID2  |
    +--------+-----+-----+----+----+------+------+------+-------+
    |        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    +--------+-----+-----+----+----+------+------+------+-------+

    .. warning:: NX does not support SORL and LDIR, MSC does
    """
    card_name = 'PLOAD4'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.eid = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.pressures = np.array([], dtype='float64')
        self.nvector = np.array([], dtype='float64')
        self.g1 = np.array([], dtype='int32')
        self.g34 = np.array([], dtype='int32')
        self.surf_or_line = np.array([], dtype='|U8')
        self.line_load_dir = np.array([], dtype='|U8')

        self._eid = []
        self._cid = []
        self._pressures = []
        self._g1 = []
        self._g34 = []
        self._nvector = []
        self._surf_or_line = []
        self._line_load_dir = []
        self.comment = defaultdict(str)

    def add(self, sid, eid, pressures, g1, g34,
            cid=0, nvector=None, surf_or_line='SURF',
            line_load_dir='NORM', comment=''):
        """
        Creates a PLOAD4 card

        Parameters
        ----------
        sid : int
            the load id
        eids : list[int, ...]
            shells : the range of element ids; must be sequential
            solids : must be length 1
        pressures : list[float, float, float, float]
            tri : must be length 4 (the last value should be the same as the 0th value)
            quad : must be length 4
        g1 : int/None
            only used for solid elements
        g34 : int / None
            only used for solid elements
        cid : int; default=0
            the coordinate system for ???
        nvector : (3, ) float ndarray
           blank : load acts normal to the face
           the local pressure vector
        surf_or_line : str; default='SURF'
           SURF : surface load
           LINE : line load    (only defined for QUADR, TRIAR)
           not supported
        line_load_dir : str; default='NORM'
           direction of the line load (see surf_or_line); {X, Y, Z, TANG, NORM}
           not supported
        comment : str; default=''
            a comment for the card

        TODO: fix the way "pressures" works
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._eid.append(eid)
        self._cid.append(cid)
        self._pressures.append(pressures)
        self._g1.append(g1)
        self._g34.append(g34)
        self._nvector.append(nvector)
        self._surf_or_line.append(surf_or_line)
        self._line_load_dir.append(line_load_dir)

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a PLOAD4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        p1 = double_or_blank(card, 3, 'p1', 0.0)
        pressures = [
            p1,
            double_or_blank(card, 4, 'p2', p1),
            double_or_blank(card, 5, 'p3', p1),
            double_or_blank(card, 6, 'p4', p1)]

        eids = [eid]
        g1_thru = integer_string_or_blank(card, 7, 'g1/THRU')
        if g1_thru == 'THRU' and integer_or_blank(card, 8, 'eid2'):
            # alternate form
            eid2 = integer(card, 8, 'eid2')
            if eid2:
                eids = list(np.unique(
                    expand_thru([eid, 'THRU', eid2], set_fields=False, sort_fields=False)
                ))
            g1 = 0
            g34 = 0
        else:
            # standard form
            eids = [eid]
            g1 = integer_or_blank(card, 7, 'g1', -1)
            g34 = integer_or_blank(card, 8, 'g34', -1)

        cid = integer_or_blank(card, 9, 'cid', 0)
        nvector = np.array([double_or_blank(card, 10, 'N1', 0.0),
                            double_or_blank(card, 11, 'N2', 0.0),
                            double_or_blank(card, 12, 'N3', 0.0)])
        surf_or_line = string_or_blank(card, 13, 'sorl', 'SURF')
        line_load_dir = string_or_blank(card, 14, 'ldir', 'NORM')
        assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)
        for eid in eids:
            self.add(sid, eid, pressures, g1, g34, cid, nvector,
                     surf_or_line, line_load_dir, comment=comment)
            comment = ''

    @property
    def is_solid(self):
        return np.where(self.g1 == 0)[0]

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for sid, eid, cid, g1, g34, pressures, nvector, surf_or_line, line_load_dir in zip(
            self.sid, self.eid, self.cid, self.g1, self.g34, self.pressures,
            self.nvector, self.surf_or_line, self.line_load_dir):
            p1 = pressures[0]
            p2 = set_blank_if_default(pressures[1], p1)
            p3 = set_blank_if_default(pressures[2], p1)
            p4 = set_blank_if_default(pressures[3], p1)
            list_fields = ['PLOAD4', sid, eid, pressures[0], p2, p3, p4]

            if g1 != -1:
                # is it a SOLID element
                list_fields += [g1, g34]
            else:
                list_fields.extend([eid, None])
                #if len(eids) > 1:
                    #try:
                        #list_fields.append('THRU')
                        #eidi = eids[-1]
                    #except Exception:
                        #print("g1  = %s" % g1)
                        #print("g34 = %s" % g34)
                        #print("eid = %s" % eid)
                        #raise
                    #list_fields.append(eidi)
                #else:
                    #list_fields += [None, None]

            if cid or np.linalg.norm(nvector) > 0.0:
                n1 = nvector[0]
                n2 = nvector[1]
                n3 = nvector[2]
                list_fields.append(cid)
                list_fields += [n1, n2, n3]
                surf_or_line = surf_or_line
                line_load_dir = line_load_dir
            else:
                list_fields += [None, None, None, None]
                surf_or_line = set_blank_if_default(surf_or_line, 'SURF')
                line_load_dir = set_blank_if_default(line_load_dir, 'NORM')
            list_fields.append(surf_or_line)
            if surf_or_line == 'LINE':
                list_fields.append(line_load_dir)

            msgi = print_card_8(list_fields)
            #msg += self.comment[eid] + msgi.rstrip() + '\n'
            msg += msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.eid = np.vstack([self.eid, self._eid])
                self.cid = np.hstack([self.cid, self._cid])

                self.pressures = np.hstack([self.pressures, self._pressures])
                self.g1 = np.hstack([self.g1, self._g1])
                self.g34 = np.hstack([self.g34, self._g34])
                self.nvector = np.hstack([self.nvector, self._nvector])
                self.surf_or_line = np.hstack([self.surf_or_line, self._surf_or_line])
                self.line_load_dir = np.hstack([self.line_load_dir, self._line_load_dir])

                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.eid = np.array(self._eid, dtype='int32')
                self.cid = np.array(self._cid, dtype='int32')
                self.pressures = np.array(self._pressures, dtype='float64')
                self.g1 = np.array(self._g1, dtype='int32')
                self.g34 = np.array(self._g34, dtype='int32')
                self.nvector = np.array(self._nvector, dtype='float64')
                self.surf_or_line = np.array(self._surf_or_line, dtype='|U8')
                self.line_load_dir = np.array(self._line_load_dir, dtype='|U8')

            self._sid = []
            self._eid = []
            self._cid = []
            self._pressures = []
            self._g1 = []
            self._g34 = []
            self._nvector = []
            self._surf_or_line = []
            self._line_load_dir = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sPLOAD4v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg
