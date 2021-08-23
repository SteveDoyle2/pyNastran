from __future__ import annotations
from typing import TYPE_CHECKING
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


class Loads:
    """initializes the Loads"""
    def __init__(self, model):
        """
        This is mostly me thinking about the problem, not what the code does...

        In the BDF, we create multiple Loads objects that is stored in a
        dictionary called self.loads[sid].  Thus, we don't need to store
        the sid locally (we do though).

        We also store a self.load_combinations, which stores only LOAD cards.
        Then, in the case control deck, there are a few cases:

          - The LOAD=1 is set:
                1) no LOAD card exists
                2) a FORCE is applied
                Currently handled in the non-vectorized BDF

          - The LOAD=1 is set:
                1) a LOAD card exists
                2) a FORCE is applied
                3) the LOAD doesn't reference the FORCE
            Currently handled in the non-vectorized BDF

          - The LOAD=1 is set:
                1) a LOAD card exists
                2) a FORCE is applied
                3) the LOAD doesn't reference the FORCE
            Currently **not** handled in the non-vectorized BDF


        By using a self.load_combinations and self.loads (or some other similar
        names), we need to:

          - Prevent recursion.
            (Currently handled in the non-vectorized BDF)

          - Prevent secondary LOAD cards from being called.
            (Currently **not** handled in the non-vectorized BDF)

          - Allow for LOAD=1 to reference and scale FORCE=1, which is bizarre,
            but allowed.
            (Currently **not** handled in the non-vectorized BDF)
        """
        self.model = model
        self.sload = model.sload
        self.pload = model.pload
        self.pload1 = model.pload1
        self.pload2 = model.pload2
        self.pload4 = model.pload4
        self.force = model.force
        self.force1 = model.force1
        self.force2 = model.force2
        self.moment = model.moment
        self.moment1 = model.moment1
        self.moment2 = model.moment2
        self.spcd = model.spcd
        self.temp = model.temp
        self.tempd = model.tempd
        self.lseq = model.lseq
        #self.rforce = model.rforce
        #self.rforce1 = model.rforce1
        self.grav = model.grav
        #self.presax = model.presax
        #self.gmload = model.gmload
        self.unhandled = []

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.sload):
            self.sload.write_card(size, is_double, bdf_file)
        if len(self.pload):
            self.pload.write_card(size, is_double, bdf_file)
        if len(self.pload1):
            self.pload1.write_card(size, is_double, bdf_file)
        if len(self.pload2):
            self.pload2.write_card(size, is_double, bdf_file)
        if len(self.pload4):
            self.pload4.write_card(size, is_double, bdf_file)

        if len(self.force):
            self.force.write_card(size, is_double, bdf_file)
        if len(self.force1):
            self.force1.write_card(size, is_double, bdf_file)
        if len(self.force2):
            self.force2.write_card(size, is_double, bdf_file)

        if len(self.moment):
            self.moment.write_card(size, is_double, bdf_file)
        if len(self.moment1):
            self.moment1.write_card(size, is_double, bdf_file)
        if len(self.moment2):
            self.moment2.write_card(size, is_double, bdf_file)
        if len(self.spcd):
            self.spcd.write_card(size, is_double, bdf_file)

        if len(self.temp):
            self.temp.write_card(size, is_double, bdf_file)
        if len(self.tempd):
            self.tempd.write_card(size, is_double, bdf_file)
        if len(self.lseq):
            self.lseq.write_card(size, is_double, bdf_file)

        #if len(self.rforce):
            #self.rforce.write_card(size, is_double, bdf_file)
        #if len(self.rforce1):
            #self.rforce1.write_card(size, is_double, bdf_file)
        if len(self.grav):
            self.grav.write_card(size, is_double, bdf_file)
        #if len(self.presax):
            #self.presax.write_card(size, is_double, bdf_file)
        #if len(self.gmload):
            #self.gmload.write_card(size, is_double, bdf_file)

    #def make_current(self):
        #"""calls make_current() for each group"""
        #self.eids = []
        #for group in self.groups:
            #group.make_current()

    @property
    def groups(self):
        """gets the sub-load groups"""
        groups = [
            self.sload,
            self.pload, self.pload1, self.pload2, self.pload4,
            self.force, self.force1, self.force2,
            self.moment, self.moment1, self.moment2,
            self.spcd, self.temp, self.tempd,
            self.grav, self.lseq,#self.rforce, self.rforce1,
            #self.presax,
        ]
        return groups

    #@property
    #def loads(self):
        #"""gets all the loads in the expected sorted order"""
        #elems = []
        #for group in self.groups:
            #elems += group.elements
        #return elems

    def __len__(self):
        """gets the number of loads"""
        return sum([len(group) for group in self.groups])

    def repr_indent(self, indent='  '):
        msg = '%s<Loads> : nloads=%s\n' % (indent, len(self))
        msg += '%s  SLOAD  :  %s\n' % (indent, len(self.sload))
        msg += '%s  PLOAD  :  %s\n' % (indent, len(self.pload4))
        msg += '%s  PLOAD1 :  %s\n' % (indent, len(self.pload1))
        msg += '%s  PLOAD2 :  %s\n' % (indent, len(self.pload2))
        msg += '%s  PLOAD4 :  %s\n' % (indent, len(self.pload4))
        msg += '%s  FORCE  :  %s\n' % (indent, len(self.force))
        msg += '%s  FORCE1 :  %s\n' % (indent, len(self.force1))
        msg += '%s  FORCE2 :  %s\n' % (indent, len(self.force2))
        msg += '%s  MOMENT :  %s\n' % (indent, len(self.moment))
        msg += '%s  MOMENT1:  %s\n' % (indent, len(self.moment1))
        msg += '%s  MOMENT2:  %s\n' % (indent, len(self.moment2))
        msg += '%s  SPCD   :  %s\n' % (indent, len(self.spcd))
        msg += '%s  TEMP   :  %s\n' % (indent, len(self.temp))
        msg += '%s  TEMPD  :  %s\n' % (indent, len(self.tempd))
        msg += '%s  GRAV   :  %s\n' % (indent, len(self.grav))
        #msg += '%s  RFORCE:  %s\n' % (indent, len(self.rforce))
        #msg += '%s  RFORCE1:  %s\n' % (indent, len(self.rforce1))
        #msg += '%s  GMLOAD:  %s\n' % (indent, len(self.gmload))
        #msg += '%s  PRESAX:  %s\n' % (indent, len(self.presax))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class BaseLoad:
    """base class for FORCE, PLOAD4"""
    card_name = ''
    def __init__(self, model):
        self.model = model
        self.is_current = True
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

    def cross_reference(self, model: BDF) -> None:
        """does this do anything?"""
        self.make_current()

    def __len__(self):
        """returns the number of elements"""
        return len(self.sid) + len(self._sid)

    def __repr__(self):
        return self.repr_indent(indent='')

    def make_current(self):
        raise NotImplementedError(self.card_name)
    def repr_indent(self, indent=''):
        raise NotImplementedError(self.card_name)


class FORCEv(BaseLoad):
    """
    Defines a static concentrated force at a grid point by specifying a
    scale factor and a vector that determines the direction.

    +-------+-----+------+-------+------+------+------+------+
    |   1   |  2  |  3   |   4   |  5   |  6   |   7  |   8  |
    +=======+=====+======+=======+======+======+======+======+
    | FORCE | SID | NODE |  CID  | MAG  |  FX  |  FY  |  FZ  |
    +-------+-----+------+-------+------+------+------+------+
    | FORCE |  3  |  1   |       | 100. |  0.  |  0.  |  1.  |
    +-------+-----+------+-------+------+------+------+------+
    """
    card_name = 'MOMENT'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.xyz = np.array([], dtype='float64')

        self._nid = []
        self._cid = []
        self._mag = []
        self._xyz = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, xyz, cid=0, comment=''):
        """
        Creates a FORCE card

        Parameters
        ----------
        sid : int
            load id
        nid : int
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
        xyz = np.array([double_or_blank(card, 5, 'X1', 0.0),
                        double_or_blank(card, 6, 'X2', 0.0),
                        double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(FORCE card) = %i\ncard=%s' % (len(card), card)
        #return FORCE(sid, node, mag, xyz, cid=cid, comment=comment)
        self.add(sid, node, mag, xyz, cid=cid, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, node_id, cid, mag, xyz in zip(count(), self.sid, self.nid, self.cid, self.mag, self.xyz):
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
            msg += self.comment[i] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
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
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.g12 = np.array([], dtype='int32')

        self._nid = []
        self._mag = []
        self._g12 = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, g1, g2, comment=''):
        """
        Creates a FORCE1 card

        Parameters
        ----------
        sid : int
            load id
        nid : int
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
        self.make_current()
        msg = ''
        for i, sid, node_id, mag, g12 in zip(count(), self.sid, self.nid, self.mag, self.g12):
            nid1, nid2 = g12
            if size == 8:
                #msgi = print_card_8(list_fields)
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

    def make_current(self):
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
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.g1234 = np.array([], dtype='int32')

        self._nid = []
        self._mag = []
        self._g1234 = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, g1, g2, g3, g4, comment=''):
        """
        Creates a FORCE2 card

        Parameters
        ----------
        sid : int
            load id
        nid : int
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
        self.make_current()
        msg = ''
        for i, sid, node_id, mag, g1234 in zip(count(), self.sid, self.nid, self.mag, self.g1234):
            nid1, nid2, nid3, nid4 = g1234
            list_fields = ['FORCE2', sid, node_id, mag, nid1, nid2, nid3, nid4]
            msgi = print_card_8(list_fields)
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

    def make_current(self):
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

    def repr_indent(self, indent=''):
        msg = '%sFORCE2v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class MOMENTv(BaseLoad):
    """
    Defines a static concentrated moment at a grid point by specifying a
    scale factor and a vector that determines the direction.

    +--------+-----+---+-----+-----+-----+-----+-----+
    |   1    |  2  | 3 |  4  |  5  |  6  |  7  |  8  |
    +========+=====+===+=====+=====+=====+=====+=====+
    | MOMENT | SID | G | CID |  M  |  N1 |  N2 |  N3 |
    +--------+-----+---+-----+-----+-----+-----+-----+
    | MOMENT |  2  | 5 |  6  | 2.9 | 0.0 | 1.0 | 0.0 |
    +--------+-----+---+-----+-----+-----+-----+-----+
    """
    card_name = 'FORCE'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.xyz = np.array([], dtype='float64')

        self._nid = []
        self._cid = []
        self._mag = []
        self._xyz = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, xyz, cid=0, comment=''):
        """
        Creates a MOMENT card

        Parameters
        ----------
        sid : int
            load id
        nid : int
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
        Adds a MOMENT card from ``BDF.add_card(...)``

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
        xyz = np.array([double_or_blank(card, 5, 'X1', 0.0),
                        double_or_blank(card, 6, 'X2', 0.0),
                        double_or_blank(card, 7, 'X3', 0.0)])
        assert len(card) <= 8, 'len(MOMENT card) = %i\ncard=%s' % (len(card), card)
        self.add(sid, node, mag, xyz, cid=cid, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, node_id, cid, mag, xyz in zip(count(), self.sid, self.nid, self.cid, self.mag, self.xyz):
            if size == 8:
                cids = set_string8_blank_if_default(cid, 0)
                msgi = 'MOMENT  %8i%8i%8s%8s%8s%8s%8s\n' % (
                    sid, node_id,
                    cids, print_float_8(mag), print_float_8(xyz[0]),
                    print_float_8(xyz[1]), print_float_8(xyz[2]))
            else:
                cids = set_string16_blank_if_default(cid, 0)
                if is_double:
                    msgi = ('MOMENT* %16i%16i%16s%s\n'
                            '*       %16s%16s%16s\n') % (
                                sid, node_id,
                                cids, print_scientific_double(mag),
                                print_scientific_double(xyz[0]),
                                print_scientific_double(xyz[1]),
                                print_scientific_double(xyz[2]))
                else:
                    msgi = ('MOMENT* %16i%16i%16s%s\n'
                            '*       %16s%16s%16s\n') % (
                                sid, node_id,
                                cids, print_float_16(mag), print_float_16(xyz[0]),
                                print_float_16(xyz[1]), print_float_16(xyz[2]))
            msg += self.comment[i] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
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

    def repr_indent(self, indent=''):
        msg = '%sMOMENTv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class MOMENT1v(BaseLoad):
    """
    Defines a static concentrated moment at a grid point by specifying a
    magnitude and two grid points that determine the direction.

    +---------+-----+----+-------+----+----+
    |    1    |  2  | 3  |   4   | 5  | 6  |
    +=========+=====+====+=======+====+====+
    | MOMENT1 | SID | G  |   M   | G1 | G2 |
    +---------+-----+----+-------+----+----+
    | MOMENT1 |  6  | 13 | -2.93 | 16 | 13 |
    +---------+-----+----+-------+----+----+
    """
    card_name = 'MOMENT1'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.g12 = np.array([], dtype='int32')

        self._nid = []
        self._mag = []
        self._g12 = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, g1, g2, comment=''):
        """
        Creates a MOMENT1 card

        Parameters
        ----------
        sid : int
            load id
        nid : int
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
        Adds a MOMENT1 card from ``BDF.add_card(...)``

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
        assert len(card) == 6, 'len(MOMENT1 card) = %i\ncard=%s' % (len(card), card)
        self.add(sid, node, mag, g1, g2, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, node_id, mag, g12 in zip(count(), self.sid, self.nid, self.mag, self.g12):
            nid1, nid2 = g12
            if size == 8:
                msgi = 'MOMENT1 %8i%8i%8s%8i%8i\n' % (
                    sid, node_id,
                    print_float_8(mag), nid1, nid2)
            else:
                if is_double:
                    msgi = ('MOMENT1*%16i%16i%16s%i\n'
                            '*       %16i\n') % (
                                sid, node_id,
                                print_scientific_double(mag),
                                nid1, nid2)
                else:
                    msgi = ('MOMENT1*%16i%16i%16s%i\n'
                            '*       %16i\n') % (
                                sid, node_id,
                                print_float_16(mag),
                                nid1, nid2)
            msg += self.comment[i] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
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

    def repr_indent(self, indent=''):
        msg = '%sMOMENTv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class MOMENT2v(BaseLoad):
    """
    Defines a static concentrated moment at a grid point by specification
    of a magnitude and four grid points that determine the direction.

    +---------+-----+---+---+----+----+----+----+
    |    1    |  2  | 3 | 4 |  5 |  6 |  7 |  8 |
    +=========+=====+===+===+====+====+====+====+
    | MOMENT2 | SID | G | M | G1 | G2 | G3 | G4 |
    +---------+-----+---+---+----+----+----+----+
    """
    card_name = 'MOMENT2'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='int32')
        self.g1234 = np.array([], dtype='int32')

        self._nid = []
        self._mag = []
        self._g1234 = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, g1, g2, g3, g4, comment=''):
        """
        Creates a MOMENT2 card

        Parameters
        ----------
        sid : int
            load id
        nid : int
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
        Adds a MOMENT2 card from ``BDF.add_card(...)``

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
        assert len(card) == 8, 'len(MOMENT2 card) = %i\ncard=%s' % (len(card), card)
        self.add(sid, node, mag, g1, g2, g3, g4, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, node_id, mag, g1234 in zip(count(), self.sid, self.nid, self.mag, self.g1234):
            nid1, nid2, nid3, nid4 = g1234
            list_fields = ['MOMENT2', sid, node_id, mag, nid1, nid2, nid3, nid4]
            msgi = print_card_8(list_fields)
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

    def make_current(self):
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

    def repr_indent(self, indent=''):
        msg = '%sMOMENT2v:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg

class SLOADv(BaseLoad):
    """
    Static Scalar Load
    Defines concentrated static loads on scalar or grid points.

    +-------+-----+----+-----+----+------+----+-------+
    |   1   |  2  | 3  |  4  |  5 |  6   |  7 |   8   |
    +=======+=====+====+=====+====+======+====+=======+
    | SLOAD | SID | S1 | F1  | S2 |  F2  | S3 |   F3  |
    +-------+-----+----+-----+----+------+----+-------+
    | SLOAD | 16  | 2  | 5.9 | 17 | -6.3 | 14 | -2.93 |
    +-------+-----+----+-----+----+------+----+-------+

    .. note:: Can be used in statics OR dynamics.

    If Si refers to a grid point, the load is applied to component T1 of the
    displacement coordinate system (see the CD field on the GRID entry).
    """
    card_name = 'SLOAD'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.mag = np.array([], dtype='float64')

        self._nid = []
        self._mag = []
        self.comment = defaultdict(str)

    def add(self, sid, nid, mag, comment=''):
        """
        Creates a SLOAD card, which applies a load to DOF=1 of the GRID/SPOINT

        Parameters
        ----------
        sid : int
            load id
        nid : int
            the GRID/SPOINT id
        mag : float
            load magnitude
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._nid.append(nid)
        self._mag.append(mag)

    def add_card(self, card, comment=''):
        """
        Adds a SLOAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')

        nfields = len(card) - 2
        ngroups = nfields // 2
        if nfields % 2 == 1:
            ngroups += 1
            msg = 'Missing last magnitude on SLOAD card=%s' % card.fields()
            raise RuntimeError(msg)

        nodes = []
        mags = []
        for i in range(ngroups):
            j = 2 * i + 2
            nodes.append(integer(card, j, 'nid' + str(i)))
            mags.append(double(card, j + 1, 'mag' + str(i)))

        for nid, mag in zip(nodes, mags):
            self.add(sid, nid, mag, comment=comment)
            comment = ''

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, nid, mag in zip(count(), self.sid, self.nid, self.mag):
            list_fields = ['SLOAD', sid, nid, mag]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.nid = np.vstack([self.nid, self._nid])
                self.mag = np.vstack([self.mag, self._mag])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.mag = np.array(self._mag, dtype='float64')

            self._sid = []
            self._nid = []
            self._mag = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sSLOADv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg

class GRAVv(BaseLoad):
    """
    Defines acceleration vectors for gravity or other acceleration loading.

    +------+-----+-----+------+-----+-----+------+-----+
    |  1   |  2  |  3  |   4  |  5  |  6  |   7  |  8  |
    +======+=====+=====+======+=====+=====+======+=====+
    | GRAV | SID | CID |  A   | N1  | N2  |  N3  |  MB |
    +------+-----+-----+------+-----+-----+------+-----+
    | GRAV | 1   | 3   | 32.2 | 0.0 | 0.0 | -1.0 |     |
    +------+-----+-----+------+-----+-----+------+-----+
    """
    def __init__(self, model):
        BaseLoad.__init__(self, model)
        #self.nid = np.array([], dtype='int32')
        self.mb = np.array([], dtype='int32')
        self.cid = np.array([], dtype='int32')
        self.scale = np.array([], dtype='float64')
        self.N = np.array([], dtype='float64')

        #self._nid = []
        self._cid = []
        self._mb = []
        self._scale = []
        self._N = []
        self.comment = defaultdict(str)

    def add(self, sid, scale, N, cid=0, mb=0, comment=''):
        """
        Creates an GRAV card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            scale factor for load
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        cid : int; default=0
            the coordinate system for the load
        mb : int; default=0
            ???
        comment : str; default=''
            a comment for the card
        """
        #if comment:
            #self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._cid.append(cid)
        self._scale.append(scale)
        self._N.append(N)
        self._mb.append(mb)

    def add_card(self, card, comment=''):
        """
        Adds a GRAV card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', 0)
        scale = double(card, 3, 'scale')
        N = np.array([double_or_blank(card, 4, 'N1', 0.0),
                      double_or_blank(card, 5, 'N2', 0.0),
                      double_or_blank(card, 6, 'N3', 0.0)])
        mb = integer_or_blank(card, 7, 'mb', 0)
        assert len(card) <= 8, 'len(GRAV card) = %i\ncard=%s' % (len(card), card)
        return self.add(sid, scale, N, cid=cid, mb=mb, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, cid, scale, N, mb in zip(count(), self.sid, self.cid, self.scale, self.N, self.mb):
            list_fields = ['GRAV', sid, cid, scale] + list(N) + [mb]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.cid = np.vstack([self.cid, self._cid])
                self.mb = np.vstack([self.mb, self._mb])
                self.scale = np.vstack([self.scale, self._scale])
                self.N = np.vstack([self.N, self._N])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.cid = np.array(self._cid, dtype='int32')
                self.mb = np.array(self._mb, dtype='int32')
                self.N = np.array(self._N, dtype='float64')
                self.scale = np.array(self._scale, dtype='float64')

            self._sid = []
            self._cid = []
            self._scale = []
            self._mb = []
            self._N = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sGRAVv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class SPCDv(BaseLoad):
    """
    Defines an enforced displacement value for static analysis and an
    enforced motion value (displacement, velocity or acceleration) in
    dynamic analysis.

     +------+-----+-----+-----+------+----+----+----+
     |   1  |  2  |  3  |  4  |   5  |  6 | 7  |  8 |
     +======+=====+=====+=====+======+====+====+====+
     | SPCD | SID |  G1 |  C1 |  D1  | G2 | C2 | D2 |
     +------+-----+-----+-----+------+----+----+----+
     | SPCD | 100 | 32  | 436 | -2.6 | 5  | 2  | .9 |
     +------+-----+-----+-----+------+----+----+----+
    """
    card_name = 'SPCD'

    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.comp = np.array([], dtype='int32')
        self.mag = np.array([], dtype='float64')

        self._nid = []
        self._comp = []
        self._mag = []
        self.comment = defaultdict(str)

    def add(self, sid, nodes, constraints, mag, comment=''):
        """
        Creates an SPCD card, which defines the degree of freedoms to be
        set during enforced motion

        Parameters
        ----------
        sid : List[int]
            constraint ids
        nodes : List[int]
            GRID/SPOINT ids
        constraints : List[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        enforced : List[float]
            the constrained value for the given node (typically 0.0)

        .. note:: len(nodes) == len(constraints) == len(enforced)
        .. warning:: Non-zero enforced deflection requires an SPC/SPC1 as well.
                     Yes, you really want to constrain the deflection to 0.0
                     with an SPC1 card and then reset the deflection using an
                     SPCD card.
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid += sid
        self._nid += nodes
        self._comp += constraints
        self._mag += mag

    def add_card(self, card, comment=''):
        """
        Adds a SPCD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = [integer(card, 1, 'sid')]
        if card.field(5) in [None, '']:
            nodes = [integer(card, 2, 'G1'),]
            constraints = [components_or_blank(card, 3, 'C1', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0)]
        else:
            sid.append(sid[0])
            nodes = [
                integer(card, 2, 'G1'),
                integer(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            constraints = [components_or_blank(card, 3, 'C1', 0),
                           components_or_blank(card, 6, 'C2', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0),
                        double_or_blank(card, 7, 'D2', 0.0)]
        self.add(sid, nodes, constraints, enforced, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, nid, comp, mag in zip(count(), self.sid, self.nid, self.comp, self.mag):
            list_fields = ['SPCD', sid]
            list_fields += [nid, comp, mag]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.nid = np.vstack([self.nid, self._nid])
                self.comp = np.vstack([self.comp, self._comp])
                self.mag = np.vstack([self.mag, self._mag])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.comp = np.array(self._comp, dtype='int32')
                self.mag = np.array(self._mag, dtype='float64')

            self._sid = []
            self._nid = []
            self._mag = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sSPCDv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class LSEQv(BaseLoad):
    """
    Defines a sequence of static load sets

    .. todo:: how does this work...
    +------+-----+----------+-----+-----+
    |   1  |  2  |     3    |  4  |  5  |
    +======+=====+==========+=====+=====+
    | LSEQ | SID | EXCITEID | LID | TID |
    +------+-----+----------+-----+-----+

    ACSRCE : If there is no LOADSET Case Control command, then EXCITEID
             may reference DAREA and SLOAD entries. If there is a LOADSET
             Case Control command, then EXCITEID may reference DAREA
             entries as well as SLOAD entries specified by the LID field
             in the selected LSEQ entry corresponding to EXCITEID.

    DAREA :  Refer to RLOAD1, RLOAD2, TLOAD1, TLOAD2, or ACSRCE entries
             for the formulas that define the scale factor Ai in dynamic
             analysis.

    DPHASE :

    SLOAD :  In the static solution sequences, the load set ID (SID) is
             selected by the Case Control command LOAD. In the dynamic
             solution sequences, SID must be referenced in the LID field
             of an LSEQ entry, which in turn must be selected by the Case
             Control command LOADSET.

    LSEQ LID : Load set identification number of a set of static load
               entries such as those referenced by the LOAD Case Control
               command.


    LSEQ,  SID, EXCITEID, LID, TID

    #--------------------------------------------------------------
    # F:\\Program Files\\Siemens\\NXNastran\\nxn10p1\\nxn10p1\\nast\\tpl\\cube_iter.dat

    DLOAD       1001     1.0     1.0   55212
    sid = 1001
    load_id = [55212] -> RLOAD2.SID

    RLOAD2,     SID, EXCITEID, DELAYID, DPHASEID,   TB,     TP,  TYPE
    RLOAD2     55212   55120              55122   55123   55124
    EXCITEID = 55120 -> DAREA.SID
    DPHASEID = 55122 -> DPHASE.SID

    DARA        SID      NID    COMP  SCALE
    DAREA      55120     913    3     9.9E+9
    SID = 55120 -> RLOAD2.SID

    DPHASE      SID     POINTID   C1    TH1
    DPHASE     55122     913       3   -90.0
    SID = 55122
    POINTID = 913 -> GRID.NID

    GRID       NID       X     Y     Z
    GRID       913      50.  0.19  -39.9
    """
    type = 'LSEQ'
    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.sid = np.array([], dtype='int32')
        self.excite_id = np.array([], dtype='int32')
        self.lid = np.array([], dtype='int32')
        self.tid = np.array([], dtype='int32')

        self._excite_id = []
        self._lid = []
        self._tid = []
        self.comment = defaultdict(str)

    def add(self, sid, excite_id, lid, tid=0, comment=''):
        """
        Creates a LSEQ card

        Parameters
        ----------
        sid : int
            loadset id; LOADSET points to this
        excite_id : int
            set id assigned to this static load vector
        lid : int
            load set id of a set of static load entries;
            LOAD in the Case Control
        tid : int; default=None
            temperature set id of a set of thermal load entries;
            TEMP(LOAD) in the Case Control
        comment : str; default=''
            a comment for the card
        """
        #if comment:
            #self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid.append(sid)
        self._excite_id.append(excite_id)
        self._lid.append(lid)
        self._tid.append(tid)

    def add_card(self, card, comment=''):
        """
        Adds an LSEQ card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        load_id = integer_or_blank(card, 3, 'lid', 0)
        temp_id = integer_or_blank(card, 4, 'tid', 0)
        if load_id != 0 and temp_id != 0:
            msg = 'LSEQ load_id/temp_id are both 0; load_id=%s temp_id=%s' % (load_id, temp_id)
            raise RuntimeError(msg)
        assert len(card) <= 5, 'len(LSEQ card) = %i\ncard=%s' % (len(card), card)
        return self.add(sid, excite_id, load_id, tid=temp_id, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, excite_id, lid, tid in zip(count(), self.sid, self.excite_id, self.lid, self.tid):
            list_fields = ['LSEQ', sid, excite_id, lid, tid]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.excite_id = np.vstack([self.excite_id, self._excite_id])
                self.lid = np.vstack([self.lid, self._lid])
                self.tid = np.vstack([self.tid, self._tid])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.excite_id = np.array(self._excite_id, dtype='int32')
                self.lid = np.array(self._lid, dtype='int32')
                self.tid = np.array(self._tid, dtype='int32')

            self._sid = []
            self._excite_id = []
            self._lid = []
            self._tid = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sLSEQv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg
