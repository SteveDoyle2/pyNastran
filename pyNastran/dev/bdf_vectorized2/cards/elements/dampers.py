from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Dampers:
    """intializes the Dampers"""
    def __init__(self, model):
        self.model = model
        self.cdamp1 = model.cdamp1
        self.cdamp2 = model.cdamp2
        self.cdamp3 = model.cdamp3
        self.cdamp4 = model.cdamp4
        self.cvisc = model.cvisc    # TODO: temp
        self.plotel = model.plotel  # TODO: temp
        #self.cdamp5 = model.cdamp5
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.cdamp1):
            self.cdamp1.write_card(size, is_double, bdf_file)
        if len(self.cdamp2):
            self.cdamp2.write_card(size, is_double, bdf_file)
        if len(self.cdamp3):
            self.cdamp3.write_card(size, is_double, bdf_file)
        if len(self.cdamp4):
            self.cdamp4.write_card(size, is_double, bdf_file)
        #if len(self.cdamp5):
            #self.cdamp5.write_card(size, is_double, bdf_file)
        if len(self.cvisc):
            self.cvisc.write_card(size, is_double, bdf_file)
        if len(self.plotel):
            self.plotel.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.cdamp1.make_current()
        self.cdamp2.make_current()
        self.cdamp3.make_current()
        self.cdamp4.make_current()
        #self.cdamp5.make_current()
        self.cvisc.make_current()
        self.plotel.make_current()

    def __len__(self):
        return(len(self.cdamp1) + len(self.cdamp2) +
               len(self.cdamp3) + len(self.cdamp4) +
               len(self.cvisc) + len(self.plotel)) #+ len(self.cdamp5)

    def repr_indent(self, indent='  '):
        msg = '%s<Dampers> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CDAMP1:  %s\n' % (indent, len(self.cdamp1))
        msg += '%s  CDAMP2:  %s\n' % (indent, len(self.cdamp2))
        msg += '%s  CDAMP3:  %s\n' % (indent, len(self.cdamp3))
        msg += '%s  CDAMP4:  %s\n' % (indent, len(self.cdamp4))
        #msg += '%s  CDAMP5:  %s\n' % (indent, len(self.cdamp5))
        msg += '%s  CVISC:  %s\n' % (indent, len(self.cvisc))
        msg += '%s  PLOTEL:  %s\n' % (indent, len(self.plotel))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class DamperElement:
    """base class for CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5"""
    card_name = ''

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
        return len(self.eid) + len(self._eid)

    def repr_indent(self, indent=''):
        self.make_current()
        neids = len(self.eid)
        if neids == 0:
            return '%s%sv; nelements=%s' % (indent, self.card_name, neids)
        msg = '%s%sv; nelements=%s\n' % (indent, self.card_name, neids)
        msg += '%s  eid = %s\n' % (indent, self.eid)

        if hasattr(self, 'pid'):
            upid = np.unique(self.pid)
            if len(upid) == 1:
                msg += '%s  upid = %s' % (indent, upid)
            else:
                msg += '%s  pid = %s' % (indent, self.pid)
        else:
            msg += '%s  b = %s' % (indent, self.b)
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class CDAMP1(DamperElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CDAMP1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    card_name = 'CDAMP1'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')
        self.dofs = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self._dofs = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, dofs, comment=''):
        """
        Creates a CDAMP1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int]
            GRID/SPOINT ids; n=2
        dofs : List[int]
            node ids; n=2
        comment : str; default=''
            a comment for the card
        """
        self.model.solids.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        self._dofs.append(dofs)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CDAMP1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'g1'), integer_or_blank(card, 5, 'g2', 0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 7, 'len(CDAMP1 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, [c1, c2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes, dofs in zip(self.eid, self.pid, self.nids, self.dofs):
            list_fields = ['CDAMP1', eid, pid, nodes[0], dofs[0], nodes[1], dofs[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                self.dofs = np.hstack([self.dofs, self._dofs])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.dofs = np.array(self._dofs, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            #print(self.nid)
            self._eid = []
            self._pid = []
            self._nids = []
            self._dofs = []
            self.is_current = True
        #else:
            #print('no GRIDs')


class CDAMP2(DamperElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CDAMP2 | EID |  B  | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    card_name = 'CDAMP2'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.b = np.array([], dtype='float64')
        self.nids = np.array([], dtype='int32')
        self.dofs = np.array([], dtype='int32')

        self._eid = []
        self._b = []
        self._nids = []
        self._dofs = []
        self.comment = defaultdict(str)

    def add(self, eid, b, nids, dofs, comment=''):
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            spring damping
        nids : List[int]
            GRID/SPOINT ids; n=2
        dofs : List[int]
            node ids; n=2
        ge : float
            ???
        s : float
            ???
        comment : str; default=''
            a comment for the card
        """
        self.model.solids.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._b.append(b)
        self._nids.append(nids)
        self._dofs.append(dofs)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CELAS2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'b')
        nids = [integer_or_blank(card, 3, 'g1', 0),
                integer_or_blank(card, 5, 'g2', 0)]
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 9, 'len(CDAMP2 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, b, nids, [c1, c2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, b, nodes, dofs in zip(self.eid, self.b, self.nids, self.dofs):
            list_fields = ['CDAMP2', eid, b, nodes[0], dofs[0], nodes[1], dofs[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.b = np.vstack([self.b, self._b])
                self.nids = np.hstack([self.nids, self._nids])
                self.dofs = np.hstack([self.dofs, self._dofs])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.b = np.array(self._b, dtype='float64')
                self.nids = np.array(self._nids, dtype='int32')
                self.dofs = np.array(self._dofs, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            #print(self.nid)
            self._eid = []
            self._b = []
            self._nids = []
            self._dofs = []
            self.is_current = True
        #else:
            #print('no GRIDs')


class CDAMP3(DamperElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    card_name = 'CDAMP3'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CELAS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int]
            SPOINT ids; n=2
        comment : str; default=''
            a comment for the card
        """
        self.model.solids.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CELAS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)

        s1 = integer_or_blank(card, 3, 's1', 0)
        s2 = integer_or_blank(card, 4, 's2', 0)
        assert len(card) <= 5, 'len(CDAMP3 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, [s1, s2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            list_fields = ['CDAMP3', eid, pid, nodes[0], nodes[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            #print(self.nid)
            self._eid = []
            self._pid = []
            self._nids = []
            self.is_current = True
        #else:
            #print('no GRIDs')

class CDAMP5(DamperElement):  # TODO: temp
    def __init__(self, model):
        self.model = model
        self.is_current = True
    def __len__(self):
        return 0

class CDAMP4(DamperElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP4 | EID |  B  | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    card_name = 'CDAMP4'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.b = np.array([], dtype='float64')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._b = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, b, nids, comment=''):
        """
        Creates a CDAMP4 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            spring damping
        nids : List[int]
            SPOINT ids; n=2
        comment : str; default=''
            a comment for the card
        """
        self.model.springs.add(eid)
        self.is_current = False
        self._eid.append(eid)
        assert isinstance(b, float), b
        self._b.append(b)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CDAMP4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'k')
        s1 = integer_or_blank(card, 3, 's1', 0)
        s2 = integer_or_blank(card, 4, 's2', 0)
        assert len(card) <= 5, 'len(CDAMP4 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, b, [s1, s2], comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, b, nodes in zip(self.eid, self.b, self.nids):
            list_fields = ['CDAMP4', eid, b, nodes[0], nodes[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.b = np.vstack([self.b, self._b])
                self.nids = np.hstack([self.nids, self._nids])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.b = np.array(self._b, dtype='float64')
                self.nids = np.array(self._nids, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            #print(self.nid)
            self._eid = []
            self._b = []
            self._nids = []
            self.is_current = True
        #else:
            #print('no GRIDs')


class CVISCv(DamperElement):
    """
    Viscous Damper Connection
    Defines a viscous damper element.

    +-------+-----+-----+----+----+
    |   1   |  2  |  3  | 4  | 5  |
    +=======+=====+=====+====+====+
    | CVISC | EID | PID | G1 | G2 |
    +-------+-----+-----+----+----+
    """
    card_name = 'CVISC'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CVISC card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PVISC)
        nids : List[int, int]
            GRID ids
        comment : str; default=''
            a comment for the card
        """
        self.model.bushes.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CVISC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 4, 'n2', 0)]
        assert len(card) <= 5, 'len(CVISC card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            list_fields = ['CVISC', eid, pid, nodes[0], nodes[1]]
            msgi = print_card_8(list_fields)
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            self._eid = []
            self._pid = []
            self._nids = []
            self.is_current = True

class PLOTELv(DamperElement):
    """
    Defines a 1D dummy element used for plotting.

    This element is not used in the model during any of the solution
    phases of a problem. It is used to simplify plotting of
    structures with large numbers of colinear grid points, where the
    plotting of each grid point along with the elements connecting
    them would result in a confusing plot.

    +--------+-----+-----+-----+
    |   1    |  2  |  3  |  4  |
    +========+=====+=====+=====+
    | PLOTEL | EID | G1  | G2  |
    +--------+-----+-----+-----+
    """
    card_name = 'PLOTEL'

    def __init__(self, model):
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')

        self._eid = []
        self._nids = []
        self.comment = defaultdict(str)

    def add(self, eid, nids, comment=''):
        """
        Adds a PLOTEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        self.model.bushes.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._nids.append(nids)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a PLOTEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
        ]
        assert len(card) <= 4, 'len(PLOTEL card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, nodes, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, nodes in zip(self.eid, self.nids):
            msgi = 'PLOTEL  %8i%8i%8i\n' % (eid, nodes[0], nodes[1])
            msg += self.comment[eid] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.nids = np.hstack([self.nids, self._nids])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))
            self._eid = []
            self._nids = []
            self.is_current = True
