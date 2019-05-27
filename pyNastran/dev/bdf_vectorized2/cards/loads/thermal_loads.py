"""
defines:
 - TEMPv
 - TEMPDv

"""
from collections import defaultdict
from itertools import count
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import integer, double
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import _format_comment
from pyNastran.dev.bdf_vectorized2.cards.loads.loads import BaseLoad


class TEMPv(BaseLoad):
    """
    Defines temperature at grid points for determination of thermal loading,
    temperature-dependent material properties, or stress recovery.

    +------+-----+----+-------+----+-------+----+----+
    |   1  |  2  |  3 |   4   |  5 |   6   |  7 |  8 |
    +======+=====+====+=======+====+=======+====+====+
    | TEMP | SID | G1 |  T1   | G2 |  T2   | G3 | T3 |
    +------+-----+----+-------+----+-------+----+----+
    | TEMP |  3  | 94 | 316.2 | 49 | 219.8 |    |    |
    +------+-----+----+-------+----+-------+----+----+
    """
    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.nid = np.array([], dtype='int32')
        self.temperatures = np.array([], dtype='int32')

        self._nid = []
        self._temperatures = []
        self.comment = defaultdict(str)

    def add(self, sid, nodes, temperatures, comment=''):
        """
        Creates an TEMP card

        Parameters
        ----------
        sid : List[int]
            constraint ids
        nodes : List[int]
            GRID/SPOINT ids
        temperatures : List[float]
            the nodal temperature
        enforced : List[float]
            the constrained value for the given node (typically 0.0)
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid += sid
        self._nid += nodes
        self._temperatures += temperatures

    def add_card(self, card, comment=''):
        """
        Adds a TEMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')

        nfields = len(card)
        assert nfields <= 8, 'len(card)=%i card=%s' % (len(card), card)

        ntemps = (nfields -2) // 2
        assert nfields % 2 == 0, card
        assert nfields // 2 > 1, card

        nodes = []
        temperatures = []
        for i in range(ntemps):
            n = i * 2 + 2
            node = integer(card, n, 'g' + str(i))
            temperaturei = double(card, n + 1, 'T' + str(i))
            nodes.append(node)
            temperatures.append(temperaturei)

        sids = [sid] * len(temperatures)
        self.add(sids, nodes, temperatures, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, nid, temperature in zip(count(), self.sid, self.nid, self.temperatures):
            list_fields = ['TEMP', sid, nid, temperature]
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
                self.temperatures = np.vstack([self.temperatures, self._temperatures])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.nid = np.array(self._nid, dtype='int32')
                self.temperatures = np.array(self._temperatures, dtype='float64')

            self._sid = []
            self._nid = []
            self._temperatures = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sTEMPv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg


class TEMPDv(BaseLoad):
    """
    Defines a temperature value for all grid points of the structural model
    that have not been given a temperature on a TEMP entry

    +-------+------+----+------+----+------+----+------+----+
    |   1   |  2   | 3  |  4   |  5 |  6   | 7  |  8   | 9  |
    +=======+======+====+======+====+======+====+======+====+
    | TEMPD | SID1 | T1 | SID2 | T2 | SID3 | T3 | SID4 | T4 |
    +-------+------+----+------+----+------+----+------+----+
    """
    def __init__(self, model):
        BaseLoad.__init__(self, model)
        self.temperatures = np.array([], dtype='int32')

        self._temperatures = []
        self.comment = defaultdict(str)

    def add(self, sid, temperatures, comment=''):
        """
        Creates an TEMPD card

        Parameters
        ----------
        sid : List[int]
            constraint ids
        temperatures : List[float]
            the nodal temperature
        """
        if comment:
            self.comment[len(self)] = _format_comment(comment)
        self.is_current = False
        self._sid += sid
        self._temperatures += temperatures

    def add_card(self, card, comment=''):
        """
        Adds a TEMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        #sids = []
        #temperatures = []
        sids = [integer(card, 1, 'sid')]
        temperatures = [double(card, 2, 'temp')]

        i = 3
        while i < len(card):
            sids.append(integer(card, i, 'sid'))
            temperatures.append(double(card, i + 1, 'temp'))
            i += 2
        assert len(sids) <= 4, sids

        self.add(sids, temperatures, comment=comment)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for i, sid, temperature in zip(count(), self.sid, self.temperatures):
            list_fields = ['TEMPD', sid, temperature]
            msgi = print_card_8(list_fields)
            msg += self.comment[i] + msgi
        bdf_file.write(msg)
        return msg

    def make_current(self):
        """creates an array of the elements"""
        if not self.is_current:
            if len(self.sid) > 0: # there are already elements in self.eid
                self.sid = np.hstack([self.sid, self._sid])
                self.temperatures = np.vstack([self.temperatures, self._temperatures])
                # TODO: need to handle comments
            else:
                self.sid = np.array(self._sid, dtype='int32')
                self.temperatures = np.array(self._temperatures, dtype='float64')

            self._sid = []
            self._temperatures = []
            self.is_current = True

    def repr_indent(self, indent=''):
        msg = '%sTEMPDv:\n' % indent
        msg += '%s  sid = %s\n' % self.sid
        return msg
