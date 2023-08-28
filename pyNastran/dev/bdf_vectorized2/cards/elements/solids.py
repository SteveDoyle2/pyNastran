from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank)
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Solids:
    def __init__(self, model):
        self.model = model
        self.ctetra4 = model.ctetra4
        self.ctetra10 = model.ctetra10
        self.chexa8 = model.chexa8
        self.chexa20 = model.chexa20
        self.cpenta6 = model.cpenta6
        self.cpenta15 = model.cpenta15
        self.cpyram5 = model.cpyram5
        self.cpyram13 = model.cpyram13
        self._eids = set()

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.ctetra4):
            self.ctetra4.write_card(size, is_double, bdf_file)
        if len(self.cpenta6):
            self.cpenta6.write_card(size, is_double, bdf_file)
        if len(self.chexa8):
            self.chexa8.write_card(size, is_double, bdf_file)
        if len(self.cpyram5):
            self.cpyram5.write_card(size, is_double, bdf_file)

        if len(self.ctetra10):
            self.ctetra10.write_card(size, is_double, bdf_file)
        if len(self.cpenta15):
            self.cpenta15.write_card(size, is_double, bdf_file)
        if len(self.chexa20):
            self.chexa20.write_card(size, is_double, bdf_file)
        if len(self.cpyram13):
            self.cpyram13.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.ctetra4.make_current()
        self.ctetra10.make_current()
        self.chexa8.make_current()
        self.chexa20.make_current()
        self.cpenta6.make_current()
        self.cpenta15.make_current()
        self.cpyram5.make_current()
        self.cpyram13.make_current()

    def __len__(self):
        return(len(self.ctetra4) + len(self.ctetra10) +
               len(self.cpenta6) + len(self.cpenta15) +
               len(self.chexa8) + len(self.chexa20) +
               len(self.cpyram5) + len(self.cpyram13))

    def repr_indent(self, indent='  '):
        msg = '%s<Solids> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CTETRA4:  %s\n' % (indent, len(self.ctetra4))
        msg += '%s  CPENTA6:  %s\n' % (indent, len(self.cpenta6))
        msg += '%s  CHEXA8:   %s\n' % (indent, len(self.chexa8))
        msg += '%s  CPYRAM5:  %s\n' % (indent, len(self.cpyram5))
        msg += '%s  CTETRA10: %s\n' % (indent, len(self.ctetra10))
        msg += '%s  CPENTA15: %s\n' % (indent, len(self.cpenta15))
        msg += '%s  CHEXA20:  %s\n' % (indent, len(self.chexa20))
        msg += '%s  CPYRAM13: %s\n' % (indent, len(self.cpyram13))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class SolidElement:
    """base class for CTETRA4, CPENTA6, CHEXA8, CPYRAM5, CTETRA10, CPENTA15, CHEXA20, CPYRAM13"""
    card_name = ''
    nnodes = 0
    def __init__(self, model):
        self.model = model
        self.is_current = False
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')
        #self.theta = np.array([], dtype='int32')  # np.nan if undefined
        #self.mcid = np.array([], dtype='int32') # -1 if undefined
        #self.theta_mcid_flag = np.array([], dtype='bool')
        #self.thickness
        #self.thickness_flag

        self._eid = []
        self._pid = []
        self._nids = []
        #self._theta = []
        #self._mcid = []
        #self._theta_mcid_flag = []
        #self.thickness
        #self.thickness_flag
        self.comment = defaultdict(str)
        self.is_current = True

    def add(self, eid, pid, nids, comment=''):
        """
        Creates a CTETRA4, CPENTA6, CHEXA8, CPYRAM5, CTETRA10, CPENTA15, CHEXA20, or CPYRAM13

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : list[int]
            node ids; n=4
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

            isort = np.argsort(self.eid)
            self.eid = self.eid[isort]
            self.pid = self.pid[isort]
            self.nids = self.nids[isort, :]

            self._eid = []
            self._pid = []
            self._nids = []
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
            return '%s%s%sv; nelements=%s' % (indent, self.card_name, self.nnodes, neids)
        msg = '%s%s%sv; nelements=%s\n' % (indent, self.card_name, self.nnodes, neids)
        msg += '%s  eid = %s\n' % (indent, self.eid)

        upid = np.unique(self.pid)
        if len(upid) == 1:
            msg += '%s  upid = %s\n' % (indent, upid)
        else:
            msg += '%s  pid = %s' % (indent, self.pid)

        #umcid = np.unique(self.mcid)
        #if len(umcid) == 1 and umcid[0] == 0:
            #msg += '  umcid = %s\n' % umcid
        #else:
            #msg += '  umcid = %s\n' % umcid
            #msg += '  mcid = %s\n' % self.mcid

        #utheta = np.unique(self.theta)
        #if len(utheta) == 1 and umcid[0] == 0:
            #msg += '  utheta = %s\n' % utheta
        #else:
            #msg += '  theta = %s\n' % self.theta
        #msg += '  is_theta = %s\n' % self.is_theta
        #msg += '  nid =\n%s' % self.nid
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')


class CTETRA4v(SolidElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CTETRA | EID | PID | G1 | G2 | G3 | G4 |
    +--------+-----+-----+----+----+----+----+
    """
    card_name = 'CTETRA'
    nnodes = 4

    def add_card(self, card, comment=''):
        """
        Adds a CTETRA4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'nid1'),
                integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'),
                integer(card, 6, 'nid4'), ]
        assert len(card) == 7, 'len(CTETRA4 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nids in zip(self.eid, self.pid, self.nids):
            data = [eid, pid] + nids.tolist()
            msgi = 'CTETRA  %8i%8i%8i%8i%8i%8i\n' % tuple(data)
            msg += self.comment[eid] + msgi
        bdf_file.write(msg)
        return msg

class CTETRA10v(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7 |  8  |  9  |
    +========+=====+=====+=====+=====+=====+====+=====+=====+
    | CTETRA | EID | PID | G1  | G2  | G3  | G4 | G5  | G6  |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |        | G7  |  G8 | G9  | G10 |     |    |     |     |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    | CTETRA | 1   | 1   | 239 | 229 | 516 | 99 | 335 | 103 |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |        | 265 | 334 | 101 | 102 |     |    |     |     |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    """
    card_name = 'CTETRA'
    nnodes = 10
    nrequired = 4

    def add_card(self, card, comment=''):
        """
        Adds a CTETRA10 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'nid1'),
                integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'),
                integer(card, 6, 'nid4'),
                integer_or_blank(card, 7, 'nid5'),
                integer_or_blank(card, 8, 'nid6'),
                integer_or_blank(card, 9, 'nid7'),
                integer_or_blank(card, 10, 'nid8'),
                integer_or_blank(card, 11, 'nid9'),
                integer_or_blank(card, 12, 'nid10'), ]
        assert len(card) <= 13, 'len(CTETRA10 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            #data = [eid, pid] + nids.tolist()
            nodes2 = ['' if node is None else '%8i' % node for node in nodes[4:].tolist()]
            data = [eid, pid] + nodes[:4].tolist() + nodes2
            msgi = ('CTETRA  %8i%8i%8i%8i%8i%8i%8s%8s\n'
                    '        %8s%8s%8s%8s' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg


class CPENTA6v(SolidElement):
    r"""
    +--------+-----+-----+----+----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +========+=====+=====+====+====+====+====+====+====+
    | CPENTA | EID | PID | G1 | G2 | G3 | G4 | G5 | G6 |
    +--------+-----+-----+----+----+----+----+----+----+

    ::
         3         6
        *----------*
       / \        / \
      / A \      / c \
      *---*-----*-----*
      1    2    4      5
      V = (A1+A2)/2  * norm(c1-c2)
      C = (c1-c2)/2
    """
    card_name = 'CPENTA'
    nnodes = 6

    def add_card(self, card, comment=''):
        """
        Adds a CPENTA6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
        ]
        assert len(card) == 9, 'len(CPENTA6 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nids in zip(self.eid, self.pid, self.nids):
            data = [eid, pid] + nids.tolist()
            msgi = 'CPENTA  %8i%8i%8i%8i%8i%8i%8i%8i\n' % tuple(data)
            msg += self.comment[eid] + msgi
        bdf_file.write(msg)
        return msg


class CPENTA15v(SolidElement):
    """
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    |    1    |  2  |  3  |  4 |  5  |  6  |  7  |  8  |  9  |
    +=========+=====+=====+====+=====+=====+=====+=====+=====+
    |  CPENTA | EID | PID | G1 | G2  | G3  | G4  | G5  | G6  |
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    |         | G7  | G8  | G9 | G10 | G11 | G12 | G13 | G14 |
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    |         | G15 |     |    |     |     |     |     |     |
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    """
    card_name = 'CPENTA'
    nnodes = 15
    nrequired = 6

    def add_card(self, card, comment=''):
        """
        Adds a CPENTA15 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
            integer_or_blank(card, 9, 'nid7'),
            integer_or_blank(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            integer_or_blank(card, 15, 'nid13'),
            integer_or_blank(card, 16, 'nid14'),
            integer_or_blank(card, 17, 'nid15'),
        ]
        assert len(card) <= 18, 'len(CPENTA15 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            nodes2 = ['' if node is None else '%8i' % node for node in nodes[6:].tolist()]
            data = [eid, pid] + nodes[:6].tolist() + nodes2
            msgi = ('CPENTA  %8i%8i%8i%8i%8i%8i%8i%8i\n'
                    '        %8s%8s%8s%8s%8s%8s%8s%8s\n'
                    '        %8s' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg


class CHEXA8v(SolidElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CTETRA | EID | PID | G1 | G2 | G3 | G4 |
    +--------+-----+-----+----+----+----+----+
    """
    card_name = 'CHEXA'
    nnodes = 8

    def add_card(self, card, comment=''):
        """
        Adds a CHEXA8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'),
            integer(card, 10, 'nid8')
        ]
        assert len(card) == 11, 'len(CHEXA8 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nids in zip(self.eid, self.pid, self.nids):
            data = [eid, pid] + nids.tolist()
            msgi = ('CHEXA   %8i%8i%8i%8i%8i%8i%8i%8i\n'
                    '        %8i%8i\n' % tuple(data))
            msg += self.comment[eid] + msgi
        bdf_file.write(msg)
        return msg


class CHEXA20v(SolidElement):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+=====+
    | CHEXA | EID | PID | G1  | G2  | G3  | G4  | G5  | G6  |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |       | G7  | G8  | G9  | G10 | G11 | G12 | G13 | G14 |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |       | G15 | G16 | G17 | G18 | G19 | G20 |     |     |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    """
    card_name = 'CHEXA'
    nnodes = 20
    nrequired = 8

    def add_card(self, card, comment=''):
        """
        Adds a CHEXA20 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            integer_or_blank(card, 15, 'nid13'),
            integer_or_blank(card, 16, 'nid14'),
            integer_or_blank(card, 17, 'nid15'),
            integer_or_blank(card, 18, 'nid16'),
            integer_or_blank(card, 19, 'nid17'),
            integer_or_blank(card, 20, 'nid18'),
            integer_or_blank(card, 21, 'nid19'),
            integer_or_blank(card, 22, 'nid20'),
        ]
        assert len(card) <= 23, 'len(CHEXA20 card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            nodes2 = ['' if node is None else '%8i' % node for node in nodes[8:].tolist()]
            data = [eid, pid] + nodes[:8].tolist() + nodes2
            msgi = ('CHEXA   %8i%8i%8i%8i%8i%8i%8i%8i\n'
                    '        %8i%8i%8s%8s%8s%8s%8s%8s\n'
                    '        %8s%8s%8s%8s%8s%8s' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg


class CPYRAM5v(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
    +========+=====+=====+=====+=====+=====+=====+=====+
    | CPYRAM | EID | PID | G1  | G2  | G3  | G4  | G5  |
    +--------+-----+-----+-----+-----+-----+-----+-----+
    """
    card_name = 'CPYRAM'
    nnodes = 5

    def add_card(self, card, comment=''):
        """
        Adds a CPYRAM5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
                integer(card, 7, 'nid5')]
        assert len(card) == 8, 'len(CPYRAM5 1card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            data = [eid, pid] + nodes.tolist()
            msgi = ('CPYRAM  %8i%8i%8i%8i%8i%8i%8i\n' % tuple(data))
            msg += self.comment[eid] + msgi
        bdf_file.write(msg)
        return msg

class CPYRAM13v(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+=====+=====+=====+=====+=====+=====+=====+
    | CPYRAM | EID | PID | G1  | G2  | G3  | G4  | G5  | G6  |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |        | G7  | G8  | G9  | G10 | G11 | G12 |     |     |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    """
    card_name = 'CPYRAM'
    nnodes = 13
    nrequired = 5

    def add_card(self, card, comment=''):
        """
        Adds a CPYRAM5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'nid1'),
                integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'),
                integer(card, 6, 'nid4'),
                integer(card, 7, 'nid5')]
        assert len(card) == 8, 'len(CPYRAM5 1card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nodes in zip(self.eid, self.pid, self.nids):
            nodes2 = ['' if node is None else '%8i' % node for node in nodes[5:].tolist()]
            data = [eid, pid] + nodes[:5].tolist() + nodes2
            msgi = ('CPYRAM  %8i%8i%8i%8i%8i%8i%8i%8s\n'
                    '        %8s%8s%8s%8s%8s%8s%s' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg
