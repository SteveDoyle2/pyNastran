from __future__ import print_function
from collections import defaultdict
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.cards.base_card import _format_comment


class Loads(object):
    """intializes the Loads"""
    def __init__(self, model):
        self.model = model
        self.pload2 = model.pload2
        self.pload4 = model.pload4
        self.force = model.force
        self.force1 = model.force1
        self.force2 = model.force2

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.pload4):
            self.pload4.write_card(size, is_double, bdf_file)
        if len(self.force):
            self.force.write_card(size, is_double, bdf_file)
        if len(self.force1):
            self.force1.write_card(size, is_double, bdf_file)
        if len(self.force2):
            self.force2.write_card(size, is_double, bdf_file)

    def __len__(self):
        return len(self.pload4) + len(self.force)

    def repr_indent(self, indent='  '):
        msg = '%s<Rods> : nelements=%s\n' % (indent, len(self))
        msg += '%s  PLOAD4:  %s\n' % (indent, len(self.pload4))
        msg += '%s  FORCE :  %s\n' % (indent, len(self.force))
        msg += '%s  FORCE1:  %s\n' % (indent, len(self.force1))
        msg += '%s  FORCE2:  %s\n' % (indent, len(self.force2))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class BaseLoad(object):
    """base class for FORCE, PLOAD4"""
    card_name = ''
    def __init__(self, model):
        self.model = model
        self.is_current = False
        self.sid = np.array([], dtype='int32')
        self._sid = []

    def check_if_current(self, eid, eids):
        """we split this up to reason about it easier"""
        if self.is_current:
            if eid in eids:
                # card exists, so we use that slot
                add_card = False
            else:
                add_card = True
        else:
            add_card = True
        return add_card

    def cross_reference(self, model):
        """does this do anything?"""
        self._make_current()

    def __len__(self):
        """returns the number of elements"""
        return len(self.sid) + len(self._sid)

    def __repr__(self):
        return self.repr_indent(indent='')

    def _make_current(self):
        raise NotImplementedError(self.card_name)
    def repr_indent(self, indent=''):
        raise NotImplementedError(self.card_name)


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
        #self.model = model
        self.is_current = False
        #self.sid = np.array([], dtype='int32')
        self.pressure = np.array([], dtype='float64')
        self.eid = np.array([], dtype='int32')

        #self._sid = []
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
        eids : List[int]
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

    def add_card(self, card, comment=''):
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
        self._make_current()
        msg = ''
        for i, sid, pressure, eid in zip(count(), self.sid, self.pressure, self.eid):
            list_fields = ['PLOAD2', sid, pressure, eid]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi
            msg += msgi
        bdf_file.write(msg)
        return msg

    def _make_current(self):
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
        #else:
            #print('no GRIDs')

    def repr_indent(self, indent=''):
        msg = '%sPLOAD2v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class PLOAD4v(BaseLoad):
    """
    Solid Format
    ============
    Defines a pressure load on a face of a CHEXA, CPENTA, or CTETRA element.

    +--------+-----+-----+----+----+------+------+------+-------+
    |   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
    +========+=====+=====+====+====+======+======+======+=======+
    | PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  |  G1  | G3/G4 |
    +--------+-----+-----+----+----+------+------+------+-------+
    |        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
    +--------+-----+-----+----+----+------+------+------+-------+

    Shell Format
    ============
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
        #self.model = model
        self.is_current = False
        #self.sid = np.array([], dtype='int32')
        self.eid = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.pressures = np.array([], dtype='float64')
        self.nvector = np.array([], dtype='float64')
        self.g1 = np.array([], dtype='int32')
        self.g34 = np.array([], dtype='int32')
        self.surf_or_line = np.array([], dtype='|U8')
        self.line_load_dir = np.array([], dtype='|U8')

        #self._sid = []
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
        eids : List[int, ...]
            shells : the range of element ids; must be sequential
            solids : must be length 1
        pressures : List[float, float, float, float]
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

    def add_card(self, card, comment=''):
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
                eids = list(unique(
                    expand_thru([eid, 'THRU', eid2], set_fields=False, sort_fields=False)
                ))
            g1 = None
            g34 = None
        else:
            # standard form
            eids = [eid]
            g1 = integer_or_blank(card, 7, 'g1')
            g34 = integer_or_blank(card, 8, 'g34')

        cid = integer_or_blank(card, 9, 'cid', 0)
        nvector = array([double_or_blank(card, 10, 'N1', 0.0),
                         double_or_blank(card, 11, 'N2', 0.0),
                         double_or_blank(card, 12, 'N3', 0.0)])
        surf_or_line = string_or_blank(card, 13, 'sorl', 'SURF')
        line_load_dir = string_or_blank(card, 14, 'ldir', 'NORM')
        assert len(card) <= 15, 'len(PLOAD4 card) = %i\ncard=%s' % (len(card), card)
        for eid in eids:
            self.add(sid, eid, pressures, g1, g34, cid, nvector,
                     surf_or_line, line_load_dir, comment=comment)
            comment = ''

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for sid, eid, cid, g1, g34, pressures, nvector in zip(self.sid, self.eid, self.cid,
                                                              self.g1, self.g34, self.pressures,
                                                              self.nvector):
            p1 = self.pressures[0]
            p2 = set_blank_if_default(self.pressures[1], p1)
            p3 = set_blank_if_default(self.pressures[2], p1)
            p4 = set_blank_if_default(self.pressures[3], p1)
            list_fields = ['PLOAD4', self.sid, eid, self.pressures[0], p2, p3, p4]

            if g1 is not None:
                # is it a SOLID element
                list_fields += [g1, g34]
            else:
                if len(eids) > 1:
                    try:
                        list_fields.append('THRU')
                        eidi = eids[-1]
                    except:
                        print("g1  = %s" % g1)
                        print("g34 = %s" % g34)
                        print("eid = %s" % eid)
                        raise
                    list_fields.append(eidi)
                else:
                    list_fields += [None, None]

            if cid or norm(nvector) > 0.0:
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

    def _make_current(self):
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
        #else:
            #print('no GRIDs')

    def repr_indent(self, indent=''):
        msg = '%sPLOAD4v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class FORCEv(BaseLoad):
    """
    +-------+-----+------+-------+------+------+------+------+
    |   1   |  2  |  3   |   4   |  5   |  6   |   7  |   8  |
    +=======+=====+======+=======+======+======+======+======+
    | FORCE | SID | NODE |  CID  | MAG  |  FX  |  FY  |  FZ  |
    +-------+-----+------+-------+------+------+------+------+

    +-------+-----+------+-------+------+------+------+------+
    | FORCE |  3  |  1   |       | 100. |  0.  |  0.  |  1.  |
    +-------+-----+------+-------+------+------+------+------+
    """
    card_name = 'FORCE'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        #self.model = model
        #self.is_current = False
        #self.sid = np.array([], dtype='int32')
        self.nid = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.xyz = np.array([], dtype='float64')

        #self._sid = []
        self._nid = []
        self._cid = []
        self._mag = []
        self._xyz = []
        self.comment = defaultdict(str)

    def add(self, sid, node, mag, xyz, cid=0, comment=''):
        """
        Creates a FORCE card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        xyz : (3, ) float ndarray
            the load direction in the cid frame
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._nid.append(nid)
        self._cid.append(cid)
        self._mag.append(mag)
        self._xyz.append(xyz)

    def add_card(self, card, comment=''):
        """
        Adds a FORCE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mag = double(card, 4, 'mag')
        xyz = array([double_or_blank(card, 5, 'X1', 0.0),
                     double_or_blank(card, 6, 'X2', 0.0),
                     double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(FORCE card) = %i\ncard=%s' % (len(card), card)
        #return FORCE(sid, node, mag, xyz, cid=cid, comment=comment)
        self.add(sid, node, mag, xyz, cid=cid, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for sid, node_id, cid, mag, xyz in zip(self.sid, self.nid, self.cid, self.mag, self.xyz):
            if size == 8:
                cids = set_string8_blank_if_default(cid, 0)
                msgi = 'FORCE   %8i%8i%8s%8s%8s%8s%8s\n' % (
                    sid, node_id,
                    cids, print_float_8(mag), print_float_8(xyz[0]),
                    print_float_8(xyz[1]), print_float_8(xyz[2]))
            else:
                cids = set_string16_blank_if_default(cid, 0)
                if is_double:
                    msgi = ('FORCE*  %16i%16i%16s%s\n'
                            '*       %16s%16s%16s\n') % (
                                sid, node_id,
                                cids, print_scientific_double(mag),
                                print_scientific_double(xyz[0]),
                                print_scientific_double(xyz[1]),
                                print_scientific_double(xyz[2]))
                else:
                    msgi = ('FORCE*  %16i%16i%16s%s\n'
                            '*       %16s%16s%16s\n') % (
                                sid, node_id,
                                cids, print_float_16(mag), print_float_16(xyz[0]),
                                print_float_16(xyz[1]), print_float_16(xyz[2]))
            #msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def _make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.nid = np.vstack([self.nid, self._nid])
                self.cid = np.hstack([self.cid, self._cid])
                self.mag = np.hstack([self.mag, self._mag])
                self.xyz = np.hstack([self.xyz, self._xyz])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.cid = np.array(self._cid, dtype='int32')
                self.mag = np.array(self._mag, dtype='float64')
                self.xyz = np.array(self._xyz, dtype='float64')

            self._sid = []
            self._nid = []
            self._cid = []
            self._mag = []
            self._xyz = []
            self.is_current = True
        #else:
            #print('no GRIDs')

    def repr_indent(self, indent=''):
        msg = '%sFORCEv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class FORCE1v(BaseLoad):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.

    +--------+-----+----+-------+----+----+
    |   1    |  2  | 3  |   4   | 5  | 6  |
    +========+=====+====+=======+====+====+
    | FORCE1 | SID | G  |   F   | G1 | G2 |
    +--------+-----+----+-------+----+----+
    | FORCE1 |  6  | 13 | -2.93 | 16 | 13 |
    +--------+-----+----+-------+----+----+
    """
    card_name = 'FORCE1'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        #self.model = model
        #self.is_current = False
        #self.sid = np.array([], dtype='int32')
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.g12 = np.array([], dtype='int32')

        #self._sid = []
        self._nid = []
        self._mag = []
        self._g12 = []
        self.comment = defaultdict(str)

    def add(self, sid, node, mag, g1, g2, comment=''):
        """
        Creates a FORCE1 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        n1 / n2 : int / int
            defines the load direction
            n = n2 - n1
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._nid.append(nid)
        self._mag.append(mag)
        self._g12.append([g1, g2])

    def add_card(self, card, comment=''):
        """
        Adds a FORCE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        assert len(card) == 6, 'len(FORCE1 card) = %i\ncard=%s' % (len(card), card)
        #return FORCE1(sid, node, mag, g1, g2, comment=comment)
        self.add(sid, node, mag, g1, g2, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for i, sid, node_id, mag, g12 in zip(count(), self.sid, self.nid, self.mag, self.g12):
            nid1, nid2 = g12
            if size == 8:
                msgi = print_card_8(list_fields)
                msgi = 'FORCE1  %8i%8i%8s%8i%8i\n' % (
                    sid, node_id,
                    print_float_8(mag), nid1, nid2)
            else:
                if is_double:
                    msgi = ('FORCE1* %16i%16i%16s%i\n'
                            '*       %16i\n') % (
                                sid, node_id,
                                print_scientific_double(mag),
                                nid1, nid2)
                else:
                    msgi = ('FORCE1* %16i%16i%16s%i\n'
                            '*       %16i\n') % (
                                sid, node_id,
                                print_float_16(mag),
                                nid1, nid2)
            msg += self.comment[i] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def _make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.nid = np.vstack([self.nid, self._nid])
                self.mag = np.hstack([self.mag, self._mag])
                self.g12 = np.hstack([self.g12, self._g12])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.mag = np.array(self._mag, dtype='float64')
                self.g12 = np.array(self._g12, dtype='int32')

            self._sid = []
            self._nid = []
            self._mag = []
            self._g12 = []
            self.is_current = True
        #else:
            #print('no GRIDs')

    def repr_indent(self, indent=''):
        msg = '%sFORCE1v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg

class FORCE2v(BaseLoad):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and four grid points that determine the direction.

    +--------+-----+---+---+----+----+----+----+
    |   1    |  2  | 3 | 4 |  5 |  6 |  7 |  8 |
    +========+=====+===+===+====+====+====+====+
    | FORCE2 | SID | G | F | G1 | G2 | G3 | G4 |
    +--------+-----+---+---+----+----+----+----+
    """
    card_name = 'FORCE2'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        #self.model = model
        #self.is_current = False
        #self.sid = np.array([], dtype='int32')
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.g1234 = np.array([], dtype='int32')

        #self._sid = []
        self._nid = []
        self._mag = []
        self._g1234 = []
        self.comment = defaultdict(str)

    def add(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        """
        Creates a FORCE2 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        g1 / g2 / g3 / g4 : int / int / int / int
            defines the load direction
            n = (g2 - g1) x (g4 - g3)
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._nid.append(nid)
        self._mag.append(mag)
        self._g1234.append([g1, g2, g3, g4])

    def add_card(self, card, comment=''):
        """
        Adds a FORCE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        g3 = integer(card, 6, 'g3')
        g4 = integer(card, 7, 'g4')
        assert len(card) == 8, 'len(FORCE2 card) = %i\ncard=%s' % (len(card), card)
        self.add(sid, node, mag, g1, g2, g3, g4, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for i, sid, node_id, mag, g1234 in zip(count(), self.sid, self.nid, self.mag, self.g1234):
            nid1, nid2, nid3, nid4 = g1234
            list_fields = ['FORCE2', sid, node_id, mag, nid1, nid2, nid3, nid4]
            #if size == 8:
                #msgi = print_card_8(list_fields)
                #msgi = 'FORCE2  %8i%8i%8s%8i%8i%8i%8i\n' % (
                    #sid, node_id,
                    #print_float_8(mag), nid1, nid2, nid3, nid4)
            #else:
                #if is_double:
                    #msgi = ('FORCE2* %16i%16i%16s%i\n'
                            #'*       %16i%16i%16i\n') % (
                                #sid, node_id,
                                #print_scientific_double(mag),
                                #nid1, nid2, nid3, nid4)
                #else:
                    #msgi = ('FORCE2* %16i%16i%16s%i\n'
                            #'*       %16i%16i%16i\n') % (
                                #sid, node_id,
                                #print_float_16(mag),
                                #nid1, nid2, nid3, nid4)
            msg += self.comment[i] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def _make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.nid = np.vstack([self.nid, self._nid])
                self.mag = np.hstack([self.mag, self._mag])
                self.g1234 = np.hstack([self.g1234, self._g1234])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.mag = np.array(self._mag, dtype='float64')
                self.g1234 = np.array(self._g1234, dtype='int32')

            self._sid = []
            self._nid = []
            self._mag = []
            self._g1234 = []
            self.is_current = True
        #else:
            #print('no GRIDs')

    def repr_indent(self, indent=''):
        msg = '%sFORCE2v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg
