from __future__ import print_function
from collections import defaultdict
import numpy as np
from pyNastran.utils import integer_types

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, blank,
    integer_double_or_blank)
from pyNastran.bdf.field_writer_8 import print_field_8, print_card_8, set_blank_if_default
from pyNastran.bdf.cards.base_card import _format_comment


class ShellElement(object):
    """base class for CTRIA3, CQUAD4"""
    card_name = ''
    def __init__(self, model):
        """intializes the ShellElement"""
        self.model = model
        self.is_current = False
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='float64')
        self.theta = np.array([], dtype='int32')  # np.nan if undefined
        self.mcid = np.array([], dtype='int32') # -1 if undefined
        #self.theta_mcid_flag = np.array([], dtype='bool')
        self.zoffset = np.array([], dtype='float64')
        self.thickness = np.array([], dtype='float64')  # np.nan is undefined

        # 0 : Ti are actual user specified thicknesses
        # 1 : Ti are fractions relative to the T value of the PSHELL
        self.thickness_flag = np.array([], dtype='int32')

        self._eid = []
        self._pid = []
        self._nids = []
        self._theta = []
        self._mcid = []
        #self._theta_mcid_flag = []
        self._zoffset = []
        self._thickness = []
        self._thickness_flag = []
        self.comment = defaultdict(str)

    @property
    def is_theta(self):
        is_theta = np.full(len(self.eid), False, dtype='bool', order='C')
        i = np.where(self.mcid == -1)[0]
        is_theta[i] = True
        return is_theta

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
        #self._make_current()
        #ieid = np.searchsorted(eid, self.eid)
        #return self[ieid]

    def _make_current(self):
        """creates an array of the GRID points"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                self.theta = np.hstack([self.theta, self._theta])
                self.mcid = np.hstack([self.mcid, self._mcid])
                self.thickness = np.vstack([self.thickness, self._thickness])
                self.thickness_flag = np.hstack([self.thickness_flag, self._thickness_flag])
                self.zoffset = np.hstack([self.zoffset, self._zoffset])
                #self.seid = np.hstack(self.seid, self._seid)
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = np.array(self._nids, dtype='int32')
                self.theta = np.array(self._theta, dtype='float64')
                self.mcid = np.array(self._mcid, dtype='int32')
                self.thickness = np.array(self._thickness, dtype='float64')
                self.thickness_flag = np.array(self._thickness_flag, dtype='int32')
                self.zoffset = np.array(self._zoffset, dtype='float64')
                #self.ps = np.array(self._ps)
                #self.seid = np.array(self._seid)
            assert len(self.eid) == len(np.unique(self.eid))

            isort = np.argsort(self.eid)
            self.eid = self.eid[isort]
            self.pid = self.pid[isort]
            self.nids = self.nids[isort, :]
            self.theta = self.theta[isort]
            self.mcid = self.mcid[isort]
            self.thickness = self.thickness[isort, :]
            self.thickness_flag = self.thickness_flag[isort]
            self.zoffset = self.zoffset[isort]

            #print(self.nid)
            self._eid = []
            self._pid = []
            self._nids = []
            self._theta = []
            self._mcid = []
            self._zoffset = []
            self._thickness = []
            self._thickness_flag = []
            #self._cd = []
            self.is_current = True
        #else:
            #print('no GRIDs')

    def cross_reference(self, model):
        """does this do anything?"""
        self._make_current()

    def __len__(self):
        """returns the number of elements"""
        return len(self.eid) + len(self._eid)

    def repr_indent(self, indent=''):
        self._make_current()
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

        umcid = np.unique(self.mcid)
        if len(umcid) == 1 and umcid[0] == -1:
            show_theta = True
        elif len(umcid) == 1 and umcid[0] == 0:
            msg += '%s  umcid = %s\n' % (indent, umcid)
            show_theta = False
        else:
            msg += '%s  umcid = %s\n' % (indent, umcid)
            msg += '%s  mcid = %s\n' % (indent, self.mcid)
            show_theta = True

        if show_theta:
            utheta = np.unique(self.theta)
            if len(utheta) == 1 and umcid[0] == 0:
                msg += '%s  utheta = %s\n' % (indent, utheta)
            else:
                msg += '%s  theta = %s\n' % (indent, self.theta)
            if umcid[0] != -1:
                msg += '%s  is_theta = %s' % (indent, self.is_theta)
        #msg += '  nid =\n%s' % self.nid
        return msg

    def __repr__(self):
        return self.repr_indent('')


class CTRIA3v(ShellElement):
    """
    +--------+-------+-------+----+----+----+------------+---------+-----+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |  9  |
    +========+=======+=======+=====+===+====+============+=========+=====+
    | CTRIA3 |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |     |
    +--------+-------+-------+----+----+----+------------+---------+-----+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |     |
    +--------+-------+-------+----+----+----+------------+---------+-----+
    """
    card_name = 'CTRIA3'
    nnodes = 3
    nthickness = 3

    def add(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
            thickness_flag=0, thickness=None, comment=''):
        """
        Creates a CTRIA3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        thickness_flag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        thickness : List[float, float, float]; default=None
            If a thickness is not supplied, then the thickness will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card
        """
        self.model.shells.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if isinstance(theta_mcid, integer_types):
            self._mcid.append(theta_mcid)
            self._theta.append(np.nan)
        else:
            self._theta.append(theta_mcid)
            self._mcid.append(-1)
        self._zoffset.append(zoffset)
        self._thickness_flag.append(thickness_flag)
        self._thickness.append(thickness)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        # type: (Any, str) -> CTRIA3
        """
        Adds a CTRIA3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        if len(card) > 5:
            theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 7, 'zoffset', 0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            thickness_flag = integer_or_blank(card, 10, 'tflag', 0)
            thickness = [
                double_or_blank(card, 11, 'T1'),
                double_or_blank(card, 12, 'T2'),
                double_or_blank(card, 13, 'T3'),
            ]
            assert len(card) <= 14, 'len(CTRIA3 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            thickness_flag = 0
            thickness = [1.0, 1.0, 1.0]
        self.add(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                 thickness_flag=thickness_flag, thickness=thickness)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for eid, pid, nids, theta, mcid, thickness_flag, thickness in zip(
            self.eid, self.pid, self.nids, self.theta, self.mcid, self.thickness_flag, self.thickness):
            #zoffset = set_blank_if_default(self.zoffset, 0.0)
            thickness_flag = set_blank_if_default(thickness_flag, 0)
            #theta_mcid = set_blank_if_default(self.Theta_mcid(), 0.0)
            zoffset = 0.
            #T1 = set_blank_if_default(self.T1, 1.0)
            #T2 = set_blank_if_default(self.T2, 1.0)
            #T3 = set_blank_if_default(self.T3, 1.0)

            if mcid != -1:
                theta = mcid
            row2_data = [theta, zoffset, # theta is theta_mcid
                         thickness_flag] + thickness.tolist()
            row2 = [print_field_8(field) for field in row2_data]
            data = [eid, pid] + nids.tolist() + row2
            msgi = ('CTRIA3  %8i%8i%8i%8i%8i%8s%8s\n'
                    '                %8s%8s%8s%8s\n' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

class CTRIA6v(ShellElement):
    """
    +--------+------------+---------+----+----+----+----+----+-----+
    |   1    |      2     |    3    |  4 |  5 |  6 | 7  | 8  |  9  |
    +========+============+=========+=====+===+====+====+====+=====+
    | CTRIA3 |    EID     |   PID   | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+------------+---------+----+----+----+----+----+-----+
    |        | THETA/MCID | ZOFFSET | T1 | T2 | T3 |    |    |     |
    +--------+------------+---------+----+----+----+----+----+-----+
    """
    card_name = 'CTRIA6'
    nnodes = 6
    nthickness = 3

    def add(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
            thickness_flag=0, thickness=None, comment=''):
        """
        Creates a CTRIA3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        thickness_flag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        thickness : List[float, float, float]; default=None
            If a thickness is not supplied, then the thickness will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card
        """
        self.model.shells.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if isinstance(theta_mcid, integer_types):
            self._mcid.append(theta_mcid)
            self._theta.append(np.nan)
        else:
            self._theta.append(theta_mcid)
            self._mcid.append(-1)
        self._zoffset.append(zoffset)
        self._thickness_flag.append(thickness_flag)
        self._thickness.append(thickness)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        # type: (Any, str) -> CTRIA3
        """
        Adds a CTRIA6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', 0),
            integer_or_blank(card, 7, 'n5', 0),
            integer_or_blank(card, 8, 'n6', 0),
        ]
        if len(card) > 9:
            theta_mcid = integer_double_or_blank(card, 9, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 10, 'zoffset', 0.0)

            thickness = [
                double_or_blank(card, 11, 'T1'),
                double_or_blank(card, 12, 'T2'),
                double_or_blank(card, 13, 'T3'),
            ]
            thickness_flag = integer_or_blank(card, 14, 'tflag', 0)
            assert len(card) <= 15, 'len(CTRIA6 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            thickness = [1.0, 1.0, 1.0]
            thickness_flag = 0
        self.add(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                 thickness_flag=thickness_flag, thickness=thickness)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for eid, pid, nids, theta, mcid, thickness_flag, thickness in zip(
            self.eid, self.pid, self.nids, self.theta, self.mcid, self.thickness_flag, self.thickness):
            #zoffset = set_blank_if_default(self.zoffset, 0.0)
            thickness_flag = set_blank_if_default(thickness_flag, 0)
            #theta_mcid = set_blank_if_default(self.Theta_mcid(), 0.0)
            if mcid == -1:
                mcid = theta # theta_mcid
            zoffset = 0.
            #T1 = set_blank_if_default(self.T1, 1.0)
            #T2 = set_blank_if_default(self.T2, 1.0)
            #T3 = set_blank_if_default(self.T3, 1.0)

            if mcid != -1:
                theta = mcid
            list_fields = (['CTRIA6', eid, pid] + nids.tolist() +
                           [mcid, zoffset] + thickness.tolist() + [thickness_flag])
            msg += self.comment[eid] + print_card_8(list_fields)
        bdf_file.write(msg)
        return msg


class CQUAD4v(ShellElement):
    """
    +--------+-------+-------+----+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |     8      |    9    |
    +========+=======+=======+=====+===+====+====+============+=========+
    | CQUAD4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 |            |         |
    +--------+-------+-------+----+----+----+----+------------+---------+
    """
    card_name = 'CQUAD4'
    nnodes = 4
    nthickness = 4

    def add(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
            thickness_flag=0, thickness=None, comment=''):
        """
        Creates a CQUAD4 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        thickness_flag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        thickness : List[float, float, float, float]; default=None
            If a thickness is not supplied, then the thickness will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card
        """
        self.model.shells.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if isinstance(theta_mcid, integer_types):
            self._mcid.append(theta_mcid)
            self._theta.append(np.nan)
        else:
            self._theta.append(theta_mcid)
            self._mcid.append(-1)
        self._zoffset.append(zoffset)
        self._thickness_flag.append(thickness_flag)
        self._thickness.append(thickness)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CQUAD4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4')]
        if len(card) > 6:
            theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 8, 'zoffset', 0.0)
            blank(card, 9, 'blank')
            thickness_flag = integer_or_blank(card, 10, 'tflag', 0)
            thickness = [
                double_or_blank(card, 11, 'T1'),
                double_or_blank(card, 12, 'T2'),
                double_or_blank(card, 13, 'T3'),
                double_or_blank(card, 14, 'T4'),
            ]
            assert len(card) <= 15, 'len(CQUAD4 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            thickness_flag = 0
            thickness = [1.0, 1.0, 1.0, 1.0]
        self.add(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                 thickness_flag=thickness_flag, thickness=thickness)

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
        #self._make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for eid, pid, nids, theta, mcid, thickness_flag, thickness in zip(
            self.eid, self.pid, self.nids, self.theta, self.mcid, self.thickness_flag, self.thickness):
            #zoffset = set_blank_if_default(self.zoffset, 0.0)
            thickness_flag = set_blank_if_default(thickness_flag, 0)
            #theta_mcid = set_blank_if_default(self.Theta_mcid(), 0.0)
            zoffset = 0.
            #T1 = set_blank_if_default(self.T1, 1.0)
            #T2 = set_blank_if_default(self.T2, 1.0)
            #T3 = set_blank_if_default(self.T3, 1.0)
            #T4 = set_blank_if_default(self.T4, 1.0)

            if mcid == -1:
                mcid = theta # theta_mcid
            row2_data = [mcid, zoffset,
                         thickness_flag] + thickness.tolist()
            row2 = [print_field_8(field) for field in row2_data]
            data = [eid, pid] + nids.tolist() + row2
            msgi = ('CQUAD4  %8i%8i%8i%8i%8i%8i%8s%8s\n'
                    '                %8s%8s%8s%8s%8s\n' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg

class CQUAD8v(ShellElement):
    """
    +--------+-------+-----+----+----+----+----+------------+-------+
    |    1   |   2   |  3  |  4 |  5 |  6 |  7 |      8     |   9   |
    +========+=======+=====+====+====+====+====+============+=======+
    | CQUAD8 |  EID  | PID | G1 | G2 | G3 | G4 |     G5     |  G6   |
    +--------+-------+-----+----+----+----+----+------------+-------+
    |        |   G7  | G8  | T1 | T2 | T3 | T4 | THETA/MCID | ZOFFS |
    +--------+-------+-----+----+----+----+----+------------+-------+
    |        | TFLAG |     |    |    |    |    |            |       |
    +--------+-------+-----+----+----+----+----+------------+-------+
    """
    card_name = 'CQUAD8'
    nnodes = 8
    nthickness = 4

    def add(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
            thickness_flag=0, thickness=None, comment=''):
        """
        Creates a CQUAD8 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int, int/None, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        thickness_flag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        thickness : List[float, float, float, float]; default=None
            If a thickness is not supplied, then the thickness will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card
        """
        self.model.shells.add(eid)
        self.is_current = False
        self._eid.append(eid)
        self._pid.append(pid)
        self._nids.append(nids)
        if isinstance(theta_mcid, integer_types):
            self._mcid.append(theta_mcid)
            self._theta.append(np.nan)
        else:
            self._theta.append(theta_mcid)
            self._mcid.append(-1)
        self._zoffset.append(zoffset)
        self._thickness_flag.append(thickness_flag)
        self._thickness.append(thickness)
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CQUAD8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5', 0),
                integer_or_blank(card, 8, 'n6', 0),
                integer_or_blank(card, 9, 'n7', 0),
                integer_or_blank(card, 10, 'n8', 0)]
        if len(card) > 11:
            thickness = [
                double_or_blank(card, 11, 'T1'),
                double_or_blank(card, 12, 'T2'),
                double_or_blank(card, 13, 'T3'),
                double_or_blank(card, 14, 'T4'),
            ]
            theta_mcid = integer_double_or_blank(card, 15, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 16, 'zoffset', 0.0)
            thickness_flag = integer_or_blank(card, 17, 'tflag', 0)
            assert len(card) <= 18, 'len(CQUAD8 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            thickness_flag = 0
            thickness = [1.0, 1.0, 1.0, 1.0]
        self.add(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                 thickness_flag=thickness_flag, thickness=thickness)

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
        #self._make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self._make_current()
        msg = ''
        for eid, pid, nids, theta, mcid, thickness_flag, thickness in zip(
            self.eid, self.pid, self.nids, self.theta, self.mcid, self.thickness_flag, self.thickness):
            #zoffset = set_blank_if_default(self.zoffset, 0.0)
            thickness_flag = set_blank_if_default(thickness_flag, 0)
            #theta_mcid = set_blank_if_default(self.Theta_mcid(), 0.0)

            if mcid == -1:
                mcid = theta # theta_mcid
            zoffset = 0.
            #T1 = set_blank_if_default(self.T1, 1.0)
            #T2 = set_blank_if_default(self.T2, 1.0)
            #T3 = set_blank_if_default(self.T3, 1.0)
            #T4 = set_blank_if_default(self.T4, 1.0)

            if mcid != -1:
                theta = mcid

            list_fields = (['CQUAD8', eid, pid] + nids.tolist() + thickness.tolist() + [
                mcid, zoffset, thickness_flag])
            msg += self.comment[eid] + print_card_8(list_fields)
        bdf_file.write(msg)
        return msg


class Shells(object):
    """
    Stores shell elements that exist in 3D space
    (e.g., not axisysmmetric elements).
    """
    def __init__(self, model):
        self.model = model
        self.ctria3 = model.ctria3
        self.cquad4 = model.cquad4
        self.ctria6 = model.ctria6
        self.cquad8 = model.cquad8
        self.cquad = model.cquad
        self._eids = set([])

    def add(self, eid):
        if eid not in self._eids:
            self._eids.add(eid)
        else:
            raise RuntimeError('eid=%s is duplicated' % eid)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        if len(self.ctria3):
            msg = self.ctria3.write_card(size, is_double, bdf_file)
        if len(self.cquad4):
            self.cquad4.write_card(size, is_double, bdf_file)
        if len(self.ctria6):
            self.ctria6.write_card(size, is_double, bdf_file)
        if len(self.cquad8):
            self.cquad8.write_card(size, is_double, bdf_file)
        if len(self.cquad):
            self.cquad.write_card(size, is_double, bdf_file)

    def __len__(self):
        return(len(self.cquad4) + len(self.cquad8) +
               len(self.ctria3) + len(self.ctria6) + len(self.cquad))

    def repr_indent(self, indent=''):
        msg = '%s<Shells> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CTRIA3: %s\n' % (indent, len(self.ctria3))
        msg += '%s  CQUAD4: %s\n' % (indent, len(self.cquad4))
        msg += '%s  CTRIA6: %s\n' % (indent, len(self.ctria6))
        msg += '%s  CQUAD8: %s\n' % (indent, len(self.cquad8))
        msg += '%s  CQUAD : %s\n' % (indent, len(self.cquad))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')
