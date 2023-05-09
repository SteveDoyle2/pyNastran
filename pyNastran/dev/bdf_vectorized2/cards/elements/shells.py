from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, blank,
    integer_double_or_blank)
from pyNastran.bdf.field_writer_8 import print_field_8, print_card_8, set_blank_if_default
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, BDFCard


class ShellElement:
    """base class for CTRIA3, CQUAD4"""
    card_name = ''
    def __init__(self, model):
        """initializes the ShellElement"""
        self.model = model
        self.is_current = True
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
    nrequired = 3

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
        nids : list[int, int, int]
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
        thickness : list[float, float, float]; default=None
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

    def add_card(self, card: BDFCard, comment: str=''):
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
        self.make_current()
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

    def quality(self, eids=None):
        """gets the quality metrics for a tri"""
        assert eids is None
        self.make_current()
        nelements = len(self)
        xyz_cid0 = self.model.nodes.xyz_cid0
        i123 = self.model.get_node_index(self.nids)
        quality = ctria3_quality(nelements, xyz_cid0, i123)
        areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai, min_edge_length = quality
        return areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai, min_edge_length

def ctria3_quality(nelements, xyz_cid0, i123):
    piover3 = np.pi / 3.
    #print('i123 = ', i123)
    p1 = xyz_cid0[i123[:, 0], :]
    p2 = xyz_cid0[i123[:, 1], :]
    p3 = xyz_cid0[i123[:, 2], :]
    #print('p1 = ', p1)
    e1 = (p1 + p2) / 2.
    e2 = (p2 + p3) / 2.
    e3 = (p3 + p1) / 2.

    #    3
    #    / \
    # e3/   \ e2
    #  /    /\
    # /    /  \
    # 1---/----2
    #    e1
    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2

    e3_p2 = e3 - p2
    e2_p1 = e2 - p1
    e1_p3 = e1 - p3

    v21 = p2 - p1
    v32 = p3 - p2
    v13 = p1 - p3
    #print(v21)
    length21 = np.linalg.norm(v21, axis=1)
    length32 = np.linalg.norm(v32, axis=1)
    length13 = np.linalg.norm(v13, axis=1)
    assert len(length13) == nelements, 'len(length13)=%s nelements=%s' % (len(length13), nelements)

    length = np.vstack([length21, length32, length13])
    #print(length)

    #min_edge_length = min(length21, length32, length13)
    min_edge_length = length.min(axis=0)
    assert len(min_edge_length) == nelements, 'len(min_edge_length)=%s nelements=%s' % (len(min_edge_length), nelements)

    areai = 0.5 * np.linalg.norm(np.cross(v21, v13), axis=1)
    assert len(areai) == nelements, 'len(areai)=%s nelements=%s' % (len(areai), nelements)
    assert len(np.linalg.norm(e2_p1, axis=1)) == nelements
    assert len(np.linalg.norm(e31, axis=1)) == nelements


    # TODO: the stuff below needs work...
    #cos_skew1 = (e2_p1 @  e31) / (np.linalg.norm(e2_p1, axis=1) * np.linalg.norm(e31, axis=1))
    #cos_skew2 = (e2_p1 @ -e31) / (np.linalg.norm(e2_p1, axis=1) * np.linalg.norm(e31, axis=1))
    #cos_skew3 = (e3_p2 @  e21) / (np.linalg.norm(e3_p2, axis=1) * np.linalg.norm(e21, axis=1))
    #cos_skew4 = (e3_p2 @ -e21) / (np.linalg.norm(e3_p2, axis=1) * np.linalg.norm(e21, axis=1))
    #cos_skew5 = (e1_p3 @  e32) / (np.linalg.norm(e1_p3, axis=1) * np.linalg.norm(e32, axis=1))
    #cos_skew6 = (e1_p3 @ -e32) / (np.linalg.norm(e1_p3, axis=1) * np.linalg.norm(e32, axis=1))
    #max_skew = np.pi / 2. - np.abs(np.arccos(np.clip([
        #cos_skew1, cos_skew2, cos_skew3,
        #cos_skew4, cos_skew5, cos_skew6], -1., 1.))).min()
    max_skew = None
    lengths = np.linalg.norm([length21, length32, length13], axis=0)
    #assert len(lengths) == 3, lengths
    #print('lengths = ', lengths)
    aspect_ratio = lengths.max() / lengths.min()
    #assert len(aspect_ratio) == nelements, 'len(aspect_ratio)=%s nelements=%s' % (len(aspect_ratio), nelements)
    #assert len(cos_skew1) == nelements, 'len(cos_skew1)=%s nelements=%s' % (len(cos_skew1), nelements)
    cos_skew1 = None

    #cos_theta1 = (v21 @ -v13) / (length21 * length13)
    #cos_theta2 = (v32 @ -v21) / (length32 * length21)
    #cos_theta3 = (v13 @ -v32) / (length13 * length32)
    #thetas = np.arccos(np.clip([cos_theta1, cos_theta2, cos_theta3], -1., 1.))
    thetas = np.array([60.])
    min_thetai = thetas.min()
    max_thetai = thetas.max()
    dideal_thetai = max(max_thetai - piover3, piover3 - min_thetai)
    #assert len(cos_theta1) == nelements, 'len(cos_theta1)=%s nelements=%s' % (len(cos_theta1), nelements)

    #theta_deg = np.degrees(np.arccos(max_cos_theta))
    #if theta_deg < 60.:
        #print('p1=%s' % xyz_cid0[p1, :])
        #print('p2=%s' % xyz_cid0[p2, :])
        #print('p3=%s' % xyz_cid0[p3, :])
        #print('theta1=%s' % np.degrees(np.arccos(cos_theta1)))
        #print('theta2=%s' % np.degrees(np.arccos(cos_theta2)))
        #print('theta3=%s' % np.degrees(np.arccos(cos_theta3)))
        #print('max_theta=%s' % theta_deg)
        #asdf
    return areai, max_skew, aspect_ratio, min_thetai, max_thetai, dideal_thetai, min_edge_length

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
    nrequired = 3

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
        nids : list[int, int, int, int/None, int/None, int/None]
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
        thickness : list[float, float, float]; default=None
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

    def add_card(self, card: Any, comment: str='') -> None:
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
        self.make_current()
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


class CTRIARv(ShellElement):
    """
    +--------+-------+-------+----+----+----+------------+---------+-----+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |  9  |
    +========+=======+=======+=====+===+====+============+=========+=====+
    | CTRIAR |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |     |
    +--------+-------+-------+----+----+----+------------+---------+-----+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |     |
    +--------+-------+-------+----+----+----+------------+---------+-----+
    """
    card_name = 'CTRIAR'
    nnodes = 6
    nrequired = 3

    def add(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.0,
            thickness_flag=0, thickness=None, comment=''):
        """
        Creates a CTRIAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int]
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
        thickness : list[float, float, float]; default=None
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
        # type: (Any, str) -> CTRIAR
        """
        Adds a CTRIAR card from ``BDF.add_card(...)``

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
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]

        theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', 0.0)
        zoffset = double_or_blank(card, 7, 'zoffset', 0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')

        thickness_flag = integer_or_blank(card, 10, 'tflag', 0)
        thickness = [
            double_or_blank(card, 11, 'T1'),
            double_or_blank(card, 12, 'T2'),
            double_or_blank(card, 13, 'T3'),]
        assert len(card) <= 14, 'len(CTRIAR card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                 thickness_flag=thickness_flag, thickness=thickness)

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
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
            list_fields = (['CTRIAR', eid, pid] + nids.tolist() +
                           [mcid, zoffset, None, None, thickness_flag] + thickness.tolist())
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
    nrequired = 4

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
        nids : list[int, int, int, int]
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
        thickness : list[float, float, float, float]; default=None
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
        #self.make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nids, theta, mcid, zoffset, thickness_flag, thickness in zip(
            self.eid, self.pid, self.nids, self.theta, self.mcid, self.zoffset,
            self.thickness_flag, self.thickness):

            theta_mcid = mcid
            if mcid == -1:
                theta_mcid = theta

            row2_data = [theta_mcid, zoffset,
                         thickness_flag] + thickness.tolist()

            #if row2_data == [0.0, 0.0, 0, 1.0, 1.0, 1.0, 1.0]:
            if row2_data == [0.0, 0.0, 0, np.nan, np.nan, np.nan, np.nan]:
                data = [eid, pid] + nids.tolist()
                msgi = ('CQUAD4  %8i%8i%8i%8i%8i%8i\n' % tuple(data))
            else:
                theta_mcid = set_blank_if_default(theta_mcid, 0.0)
                zoffset = set_blank_if_default(zoffset, 0.0)
                thickness_flag = set_blank_if_default(thickness_flag, 0)
                row2_data = [theta_mcid, zoffset,
                             thickness_flag] + thickness.tolist()
                #T1 = set_blank_if_default(self.T1, 1.0)
                #T2 = set_blank_if_default(self.T2, 1.0)
                #T3 = set_blank_if_default(self.T3, 1.0)
                #T4 = set_blank_if_default(self.T4, 1.0)
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
    nrequired = 4

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
        nids : list[int, int, int, int, int/None, int/None, int/None, int/None]
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
        thickness : list[float, float, float, float]; default=None
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
        #self.make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
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


class CQUADv(ShellElement):
    """
    +-------+-------+-----+----+------------+----+----+----+----+
    |    1  |   2   |  3  |  4 |     5      |  6 |  7 | 8  |  9 |
    +=======+=======+=====+====+============+====+====+====+====+
    | CQUAD |  EID  | PID | G1 |     G2     | G3 | G4 | G5 | G6 |
    +-------+-------+-----+----+------------+----+----+----+----+
    |       |   G7  | G8  | G9 | THETA/MCID |    |    |    |    |
    +-------+-------+-----+----+------------+----+----+----+----+

    theta_mcid is an MSC specific variable
    """
    card_name = 'CQUAD'
    nnodes = 9
    nrequired = 4

    def __init__(self, model):
        """initializes the ShellElement"""
        self.model = model
        self.is_current = True
        self.eid = np.array([], dtype='int32')
        self.pid = np.array([], dtype='int32')
        self.nids = np.array([], dtype='float64')
        self.theta = np.array([], dtype='int32')  # np.nan if undefined
        self.mcid = np.array([], dtype='int32') # -1 if undefined

        self._eid = []
        self._pid = []
        self._nids = []
        self._theta = []
        self._mcid = []
        self.comment = defaultdict(str)

    def add(self, eid, pid, nids, theta_mcid=0., comment=''):
        """
        Creates a CQUAD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int, int/None, int/None,
                    int/None, int/None, int/None]
            node ids
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
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
        if comment:
            self.comment[eid] = _format_comment(comment)

    def add_card(self, card, comment=''):
        """
        Adds a CQUAD card from ``BDF.add_card(...)``

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
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5'),
                integer_or_blank(card, 8, 'n6'),
                integer_or_blank(card, 9, 'n7'),
                integer_or_blank(card, 10, 'n8'),
                integer_or_blank(card, 11, 'n9'),]
        theta_mcid = integer_double_or_blank(card, 12, 'theta_mcid', 0.)
        assert len(card) <= 13, 'len(CQUAD card) = %i\ncard=%s' % (len(card), card)
        self.add(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)

    def make_current(self):
        """creates an array of the GRID points"""
        if not self.is_current:
            if len(self.eid) > 0: # there are already elements in self.eid
                self.eid = np.hstack([self.eid, self._eid])
                self.pid = np.vstack([self.pid, self._pid])
                self.nids = np.hstack([self.nids, self._nids])
                self.theta = np.hstack([self.theta, self._theta])
                self.mcid = np.hstack([self.mcid, self._mcid])
                # don't need to handle comments
            else:
                self.eid = np.array(self._eid, dtype='int32')
                self.pid = np.array(self._pid, dtype='int32')
                self.nids = positive_int_array(self._nids, dtype='int32')
                self.theta = np.array(self._theta, dtype='float64')
                self.mcid = np.array(self._mcid, dtype='int32')
            assert len(self.eid) == len(np.unique(self.eid))

            isort = np.argsort(self.eid)
            self.eid = self.eid[isort]
            self.pid = self.pid[isort]
            self.nids = self.nids[isort, :]
            self.theta = self.theta[isort]
            self.mcid = self.mcid[isort]

            #print(self.nid)
            self._eid = []
            self._pid = []
            self._nids = []
            self._theta = []
            self._mcid = []
            self.is_current = True

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
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
        msg = ''
        for eid, pid, nids, theta, mcid in zip(
            self.eid, self.pid, self.nids, self.theta, self.mcid):
            nodes2 = ['' if node is None else '%8i' % node for node in nids[4:].tolist()]
            if mcid == -1:
                mcid = theta

            data = [eid, pid] + nids[:4].tolist() + nodes2 + [mcid]
            msgi = ('CQUAD   %8i%8i%8i%8i%8i%8i%8s%8s\n'  # 6 nodes
                    '        %8s%8s%8s%8s\n' % tuple(data))
            msg += self.comment[eid] + msgi.rstrip() + '\n'
        bdf_file.write(msg)
        return msg


class CQUADRv(ShellElement):
    """
    +--------+-------+-------+----+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |     8      |    9    |
    +========+=======+=======+=====+===+====+====+============+=========+
    | CQUADR |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 |            |         |
    +--------+-------+-------+----+----+----+----+------------+---------+
    """
    card_name = 'CQUADR'
    nnodes = 8
    nrequired = 4

    def add(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
            thickness_flag=0, thickness=None, comment=''):
        """
        Creates a CQUADR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int]
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
        thickness : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
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
        Adds a CQUADR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),]

        theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', 0.0)
        zoffset = double_or_blank(card, 8, 'zoffset', 0.0)

        thickness_flag = integer_or_blank(card, 10, 'tflag', 0)
        thickness = [
            double_or_blank(card, 11, 'T1'),
            double_or_blank(card, 12, 'T2'),
            double_or_blank(card, 13, 'T3'),
            double_or_blank(card, 14, 'T4'),
        ]
        assert len(card) <= 15, 'len(CQUADR card) = %i\ncard=%s' % (len(card), card)
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
        #self.make_current()
        #eid = self.eid[i]
        #return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    #ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    #def __setitem__(self, i, value):
        #pass
    #def __delitem__(self, i):
        #pass
    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.make_current()
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


class Shells:
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
        self.ctriar = model.ctriar
        self.cquadr = model.cquadr
        self.cquad = model.cquad
        self._eids = set()

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
        if len(self.cquadr):
            self.cquadr.write_card(size, is_double, bdf_file)
        if len(self.ctriar):
            self.ctriar.write_card(size, is_double, bdf_file)

    def make_current(self):
        self.ctria3.make_current()
        self.cquad4.make_current()
        self.ctria6.make_current()
        self.cquad8.make_current()
        self.cquad.make_current()
        self.cquadr.make_current()
        self.ctriar.make_current()

    @property
    def elements(self):
        """gets the sub-shell element groups"""
        return self.groups

    @property
    def groups(self):
        """gets the sub-shell element groups"""
        return [self.ctria3, self.cquad4,
                self.ctria6, self.cquad8,
                self.cquad, self.cquadr, self.ctriar,]

    def __len__(self):
        """gets the number of shell elements"""
        return sum([len(group) for group in self.groups])

    def repr_indent(self, indent=''):
        msg = '%s<Shells> : nelements=%s\n' % (indent, len(self))
        msg += '%s  CTRIA3: %s\n' % (indent, len(self.ctria3))
        msg += '%s  CQUAD4: %s\n' % (indent, len(self.cquad4))
        msg += '%s  CTRIA6: %s\n' % (indent, len(self.ctria6))
        msg += '%s  CQUAD8: %s\n' % (indent, len(self.cquad8))
        msg += '%s  CQUAD : %s\n' % (indent, len(self.cquad))
        msg += '%s  CQUADR: %s\n' % (indent, len(self.cquadr))
        msg += '%s  CTRIAR: %s\n' % (indent, len(self.ctriar))
        return msg

    def __repr__(self):
        return self.repr_indent(indent='')

def positive_int_array(array_, dtype='int32'):
    """used for a CQUAD8 to support None -> 0"""
    try:
        new_array = np.array(array_, dtype=dtype)
    except TypeError:
        # TODO: potentially risky...check...
        float_array = np.array(array_, dtype='float64')
        inan = np.isnan(float_array)
        float_array[inan] = 0.
        new_array = float_array.astype(dtype)
        #i0 = np.where(new_array < 0)
        #new_array[i0] = 0
    return new_array
