# coding: utf-8
# pylint: disable=C0103
"""
All beam properties are defined in this file.  This includes:
 *   PBEAM
 *   PBEAML
 *   PBCOMP
 *   PBRSECT

All beams are LineProperty objects.
Multi-segment beams are IntegratedLineProperty objects.

"""
from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING
import numpy as np
from numpy import array, unique, argsort, allclose, ndarray

from pyNastran.bdf.bdf_interface.utils import to_fields
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.cards.properties.bars import (
    IntegratedLineProperty, LineProperty, _bar_areaL,
    get_beam_sections, split_arbitrary_thickness_section, write_arbitrary_beam_section,
    plot_arbitrary_section)
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    string, string_or_blank, double_string_or_blank)
from pyNastran.utils.mathematics import integrate_unit_line, integrate_positive_unit_line
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import (
    break_word_by_trailing_parentheses_integer_ab, break_word_by_trailing_integer)
from pyNastran.utils.numpy_utils import integer_types, float_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class PBEAM(IntegratedLineProperty):
    """
    Defines the properties of a beam element (CBEAM entry). This element may be
    used to model tapered beams.


    +-------+-------+-------+-------+-------+-------+--------+-------+--------+
    | PBEAM |  PID  |  MID  | A(A)  | I1(A) | I2(A) | I12(A) | J(A)  | NSM(A) |
    +-------+-------+-------+-------+-------+-------+--------+-------+--------+
    |       | C1(A) | C2(A) | D1(A) | D2(A) | E1(A) | E2(A)  | F1(A) | F2(A)  |
    +-------+-------+-------+-------+-------+-------+--------+-------+--------+

    The next two continuations are repeated for each intermediate station as
    described in Remark 5. and SO and X/XB must be specified.

    +----+------+----+----+----+-----+----+-----+
    | SO | X/XB | A  | I1 | I2 | I12 | J  | NSM |
    +----+------+----+----+----+-----+----+-----+
    | C1 |  C2  | D1 | D2 | E1 | E2  | F1 | F2  |
    +----+------+----+----+----+-----+----+-----+

    The last two continuations are:
    +-------+-------+-------+-------+--------+--------+-------+-------+
    |   K1  |   K2  |   S1  |   S2  | NSI(A) | NSI(B) | CW(A) | CW(B) |
    +-------+-------+-------+-------+--------+--------+-------+-------+
    | M1(A) | M2(A) | M1(B) | M2(B) | N1(A)  | N2(A)  | N1(B) | N2(B) |
    +-------+-------+-------+-------+--------+--------+-------+-------+
    """
    type = 'PBEAM'

    #_opt_map = {
        #'I1(B)' : 'i1',
        #'I1(B)',
    #}
    def update_by_pname_fid(self, pname_fid, value):
        if isinstance(pname_fid, integer_types):
            if pname_fid < 0:
                # convert fid to pname
                pname_fid = update_pbeam_negative_integer(pname_fid)

                # call this function again to update pname
                self.update_by_pname_fid(pname_fid, value)
            else:  # pragma: no cover
                msg = "property_type='PBEAM' has not implemented %r in pname_map" % pname_fid
                raise NotImplementedError(msg)

        #elif pname_fid.startswith('DIM'):
            #word, num = break_word_by_trailing_integer(pname_fid)
            #ndim = len(self.dim[0])

            #num = int(num) - 1
            #idim = num % ndim
            #istation = num // ndim
            #try:
                #dim_station = self.dim[istation]
                #dim_station[idim] = value
            #except:
                #print('pname_fid=%r num=%r ndim=%r' % (pname_fid, num, ndim))
                #print('istation=%r idim=%r' % (istation, idim))
                #print(self)
                #raise
        elif isinstance(pname_fid, str):
            if '(A)' in pname_fid:
                i = 0
                end = '(A)'
            elif '(' not in pname_fid:
                i = 0
                end = ''
            elif '(B)' in pname_fid:
                i = -1
                end = '(B)'
            elif '(' in pname_fid:
                word, num = break_word_by_trailing_parentheses_integer_ab(pname_fid)
                end = '(%i)' % num
                pname_fid = word + end
                i = num - 1
            else:
                raise NotImplementedError(pname_fid)

            if pname_fid.startswith('I12'):
                self.i12[i] = value
                word = 'I12'
            elif pname_fid.startswith('I1'):
                self.i1[i] = value
                word = 'I1'
            elif pname_fid.startswith('I2'):
                self.i2[i] = value
                word = 'I2'
            elif pname_fid.startswith('J'):
                self.j[i] = value
                word = 'J'
            elif pname_fid.startswith('A'):
                self.A[i] = value
                word = 'A'
            elif pname_fid.startswith('C1'):
                self.c1[i] = value
                word = 'C1'
            elif pname_fid.startswith('C2'):
                self.c2[i] = value
                word = 'C2'
            elif pname_fid.startswith('D1'):
                self.d1[i] = value
                word = 'D1'
            elif pname_fid.startswith('D2'):
                self.d2[i] = value
                word = 'D2'
            elif pname_fid.startswith('E1'):
                self.e1[i] = value
                word = 'E1'
            elif pname_fid.startswith('E2'):
                self.e2[i] = value
                word = 'E2'
            elif pname_fid.startswith('F1'):
                self.f1[i] = value
                word = 'F1'
            elif pname_fid.startswith('F2'):
                self.f2[i] = value
                word = 'F2'
            else:  # pragma: no cover
                msg = "property_type='PBEAM' has not implemented %r in pname_map" % (
                    pname_fid)
                raise NotImplementedError(msg)

            expected_word = word + end
            if pname_fid != expected_word:
                raise RuntimeError('pname_fid=%r expected_word=%r is invalid' % (
                    pname_fid, expected_word))
        else:  # pragma: no cover
            msg = "property_type='PBEAM' has not implemented %r in pname_map; type=%s" % (
                pname_fid, type(pname_fid))
            raise NotImplementedError(msg)

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        xxb = [0.]
        so = ['YES']
        area = [0.]
        i1 = [0.]
        i2 = [0.]
        i12 = [0.]
        j = [0.]
        return PBEAM(pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                     c1=None, c2=None, d1=None, d2=None,
                     e1=None, e2=None, f1=None, f2=None,
                     k1=1., k2=1., s1=0., s2=0.,
                     nsia=0., nsib=None, cwa=0., cwb=None,
                     m1a=0., m2a=0., m1b=None, m2b=None,
                     n1a=0., n2a=0., n1b=None, n2b=None,
                     comment='')

    def __init__(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                 c1=None, c2=None, d1=None, d2=None,
                 e1=None, e2=None, f1=None, f2=None,
                 k1=1., k2=1., s1=0., s2=0.,
                 nsia=0., nsib=None, cwa=0., cwb=None,
                 m1a=0., m2a=0., m1b=None, m2b=None,
                 n1a=0., n2a=0., n1b=None, n2b=None,
                 comment=''):
        """
        .. todo:: fix 0th entry of self.so, self.xxb

        Creates a PBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        xxb : List[float]
            The percentage locations along the beam [0., ..., 1.]
        so : List[str]
            YES, YESA, NO
        area : List[float]
            area
        i1, i2, i12, j : List[float]
            moments of inertia
        nsm : List[float]
            nonstructural mass per unit length
        c1/c2, d1/d2, e1/e2, f1/f2 : List[float]; default=None -> [0.]*nxxb
           the y/z locations of the stress recovery points
           c1 - point C.y
           c2 - point C.z
        k1 / k2 : float; default=1.
            Shear stiffness factor K in K*A*G for plane 1/2.
        s1 / s2 : float; default=0.
            Shear relief coefficient due to taper for plane 1/2.
        nsia / nsia : float; default=0. / nsia
            non structural mass moment of inertia per unit length
            about nsm center of gravity at Point A/B.
        cwa / cwb : float; default=0. / cwa
            warping coefficient for end A/B.
        m1a / m2a : float; default=0. / 0.
            y/z coordinate of center of gravity of
            nonstructural mass for end A.
        m1b / m2b : float; default=m1a / m2a
            y/z coordinate of center of gravity of
            nonstructural mass for end B.
        n1a / n2a : float; default=0. / 0.
            y/z coordinate of neutral axis for end A.
        n1b / n2b : float; default=n1a / n2a
            y/z coordinate of neutral axis for end B.
        comment : str; default=''
            a comment for the card

        """
        IntegratedLineProperty.__init__(self)
        if comment:
            self.comment = comment
        if nsib is None:
            nsib = nsia
        if cwb is None:
            cwb = cwa

        # A / B - the grid points
        # 1 / 2 - the y/z coords

        assert m1a is not None, m1a
        assert m2a is not None, m2a

        if m1b is None:
            m1b = m1a
        if m2b is None:
            m2b = m2a
        assert m1b is not None, m1b
        assert m2b is not None, m2b

        # A / B - the grid points
        # 1 / 2 - the y/z coords
        assert n1a is not None, n1a
        assert n2a is not None, n2a
        if n1b is None:
            n1b = n1a
        if n2b is None:
            n2b = n2a
        assert n1b is not None, n1b
        assert n2b is not None, n2b

        nxxb = len(xxb)
        if nsm is None:
            nsm = [0.] * nxxb
        if c1 is None:
            c1 = [None] * nxxb
        if c2 is None:
            c2 = [None] * nxxb
        if d1 is None:
            d1 = [None] * nxxb
        if d2 is None:
            d2 = [None] * nxxb
        if e1 is None:
            e1 = [None] * nxxb
        if e2 is None:
            e2 = [None] * nxxb
        if f1 is None:
            f1 = [None] * nxxb
        if f2 is None:
            f2 = [None] * nxxb

        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid

        assert isinstance(xxb, list), xxb
        assert isinstance(so, list), so
        #assert min(so) in [0., 1.], so  # YES, NO
        #assert max(so) == 1.0, so
        #print('xxb', xxb)
        assert np.allclose(min(xxb), 0.), 'pid=%s min(xxb)=%s xxb=%s' % (pid, min(xxb), xxb)  # x/L
        assert 0. <= max(xxb) <= 1.0, 'pid=%s max(xxb)=%s xxb=%s' % (pid, max(xxb), xxb)  # x/L
        assert isinstance(area, list), area
        assert isinstance(i1, list), i1
        assert isinstance(i2, list), i2
        assert isinstance(i12, list), i12
        assert isinstance(j, list), j
        assert isinstance(nsm, list), nsm
        # at least one cross section is required; even if it's empty
        # xxb[0] isn't explicitly used
        #: Section position
        #self.xxb = xxb
        #: Output flag
        #self.so = so
        #: Area
        #self.A = area
        #: Moment of Inertia about the 1 axis :math:`I_1`
        #self.i1 = i1
        #: Moment of Inertia about the 2 axis  :math:`I_2`
        #self.i2 = i2
        #: Moment of Inertia about the 12 axis  :math:`I_{12}`
        #self.i12 = i12
        #: Polar Moment of Inertia :math:`J`
        #self.j = j
        #: Non-structural mass :math:`nsm`
        #self.nsm = nsm

        assert isinstance(c1, list), c1
        assert isinstance(c2, list), c2
        assert isinstance(d1, list), d1
        assert isinstance(d2, list), d2
        assert isinstance(e1, list), e1
        assert isinstance(e2, list), e2
        assert isinstance(f1, list), f1
        assert isinstance(f2, list), f2
        #self.c1 = c1
        #self.c2 = c2
        #self.d1 = d1
        #self.d2 = d2
        #self.e1 = e1
        #self.e2 = e2
        #self.f1 = f1
        #self.f2 = f2


        #: Shear stiffness factor K in K*A*G for plane 1.
        self.k1 = k1
        #: Shear stiffness factor K in K*A*G for plane 2.
        self.k2 = k2

        #: Shear relief coefficient due to taper for plane 1.
        self.s1 = s1
        #assert self.s1 == 0., self.s1
        #: Shear relief coefficient due to taper for plane 2.
        self.s2 = s2
        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point A.
        self.nsia = nsia
        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point B.
        self.nsib = nsib

        #: warping coefficient for end A.
        self.cwa = cwa
        #: warping coefficient for end B.
        self.cwb = cwb

        #: y coordinate of center of gravity of
        #: nonstructural mass for end A.
        self.m1a = m1a
        #: z coordinate of center of gravity of
        #: nonstructural mass for end A.
        self.m2a = m2a

        #: y coordinate of center of gravity of
        #: nonstructural mass for end B.
        self.m1b = m1b
        #: z coordinate of center of gravity of
        #: nonstructural mass for end B.
        self.m2b = m2b

        #: y coordinate of neutral axis for end A.
        self.n1a = n1a
        #: z coordinate of neutral axis for end A.
        self.n2a = n2a
        #: y coordinate of neutral axis for end B.
        self.n1b = n1b
        #: z coordinate of neutral axis for end B.
        self.n2b = n2b

        nxxb = len(xxb)
        # sort xxb
        ixxb = argsort(xxb)
        self.so = array(so, dtype='|U8')[ixxb]
        self.xxb = array(xxb, dtype='float64')[ixxb]
        #print('ixxb = %s' % ixxb)
        #print('i12 = %s' % i12)

        assert len(area) == nxxb, 'pid=%s len(xxb)=%s len(A =)=%s' % (pid, nxxb, len(area))
        assert len(i1) == nxxb, 'pid=%s len(xxb)=%s len(i1 )=%s' % (pid, nxxb, len(i1))
        assert len(i2) == nxxb, 'pid=%s len(xxb)=%s len(i2 )=%s' % (pid, nxxb, len(i2))
        assert len(i12) == nxxb, 'pid=%s len(xxb)=%s len(i12)=%s' % (pid, nxxb, len(i12))
        assert len(j) == nxxb, 'pid=%s len(xxb)=%s len(j =)=%s' % (pid, nxxb, len(j))
        assert len(nsm) == nxxb, 'pid=%s len(xxb)=%s len(nsm)=%s' % (pid, nxxb, len(nsm))

        self.A = array(area, dtype='float64')[ixxb]
        self.i1 = array(i1, dtype='float64')[ixxb]
        self.i2 = array(i2, dtype='float64')[ixxb]
        self.i12 = array(i12, dtype='float64')[ixxb]
        self.j = array(j, dtype='float64')[ixxb]
        self.nsm = array(nsm, dtype='float64')[ixxb]

        #assert len(area) == nxxb, 'pid=%s len(xxb)=%s len(area)=%s' % (nxxb, len(area))
        assert len(c1) == nxxb, 'pid=%s len(xxb)=%s len(c1)=%s' % (pid, nxxb, len(c1))
        assert len(c2) == nxxb, 'pid=%s len(xxb)=%s len(c2)=%s' % (pid, nxxb, len(c2))
        assert len(d1) == nxxb, 'pid=%s len(xxb)=%s len(d1)=%s' % (pid, nxxb, len(d1))
        assert len(d2) == nxxb, 'pid=%s len(xxb)=%s len(d2)=%s' % (pid, nxxb, len(d2))
        assert len(e1) == nxxb, 'pid=%s len(xxb)=%s len(e1)=%s' % (pid, nxxb, len(e1))
        assert len(e2) == nxxb, 'pid=%s len(xxb)=%s len(e2)=%s' % (pid, nxxb, len(e2))
        assert len(f1) == nxxb, 'pid=%s len(xxb)=%s len(f1)=%s' % (pid, nxxb, len(f1))
        assert len(f2) == nxxb, 'pid=%s len(xxb)=%s len(f2)=%s' % (pid, nxxb, len(f2))

        self.c1 = array(c1, dtype='float64')[ixxb]
        self.c2 = array(c2, dtype='float64')[ixxb]
        self.d1 = array(d1, dtype='float64')[ixxb]
        self.d2 = array(d2, dtype='float64')[ixxb]
        self.e1 = array(e1, dtype='float64')[ixxb]
        self.e2 = array(e2, dtype='float64')[ixxb]
        self.f1 = array(f1, dtype='float64')[ixxb]
        self.f2 = array(f2, dtype='float64')[ixxb]
        self._interpolate_sections()
        self.mid_ref = None

    def _interpolate_sections(self):
        """
        now we interpolate to fix up missing data
        from arrays that were potentially out of order
        (they're sorted now)

        we've also already checked xxb=0.0 and xxb=1.0 for I1, I2, I12, J
        """
        loop_vars = (
            count(), self.xxb, self.A, self.i1, self.i2, self.i12, self.j, self.nsm,
            self.c1, self.c2, self.d1, self.d2, self.e1, self.e2, self.f1, self.f2
        )
        for interp_data in zip(*loop_vars):
            i, xxb, a, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2 = interp_data
            if xxb not in [0., 1.]:
                if a == 0.0:
                    self.A[i] = self.A[-1] + self.A[0] * (1 - xxb)
                if i1 == 0.0:
                    self.i1[i] = self.i1[-1] + self.i1[0] * (1 - xxb)
                if i2 == 0.0:
                    self.i2[i] = self.i2[-1] + self.i2[0] * (1 - xxb)
                if i12 == 0.0:
                    self.i12[i] = self.i12[-1] + self.i12[0] * (1 - xxb)
                if j == 0.0:
                    self.j[i] = self.j[-1] + self.j[0] * (1 - xxb)

                assert self.A[i] >= 0., self.A
                assert self.i1[i] >= 0., self.i1
                assert self.i2[i] >= 0., self.i2
                assert self.j[i] >= 0., self.j  # we check warping later

                if nsm == 0.0:
                    #print('iterpolating nsm; i=%s xxb=%s' % (i, xxb))
                    self.nsm[i] = self.nsm[-1] + self.nsm[0] * (1 - xxb)

                if c1 == 0.0:
                    self.c1[i] = self.c1[-1] + self.c1[0] * (1 - xxb)
                if c2 == 0.0:
                    self.c2[i] = self.c2[-1] + self.c2[0] * (1 - xxb)

                if d1 == 0.0:
                    self.d1[i] = self.d1[-1] + self.d1[0] * (1 - xxb)
                if d2 == 0.0:
                    self.d2[i] = self.d2[-1] + self.d2[0] * (1 - xxb)

                if e1 == 0.0:
                    self.e1[i] = self.e1[-1] + self.e1[0] * (1 - xxb)
                if e2 == 0.0:
                    self.e2[i] = self.e2[-1] + self.e2[0] * (1 - xxb)

                if f1 == 0.0:
                    self.f1[i] = self.f1[-1] + self.f1[0] * (1 - xxb)
                if f2 == 0.0:
                    self.f2[i] = self.f2[-1] + self.f2[0] * (1 - xxb)


        if self.cwa or self.cwb:  # if either is non-zero
            for i, xxbi, ji in zip(count(), self.xxb, self.j):
                if ji < 0.:
                    msg = 'Warping Check Error; j[%i] must be greater than 0.0' % i
                    msg += '  cwa=%s cwb=%s\n' % (self.cwa, self.cwb)
                    msg += '  i=%s xxb=%s j=%s; j[%i]=%s\n' % (i, xxbi, self.j, i, ji)
                    raise ValueError(msg)

    def validate(self):
        nstations = len(self.A)
        for ilayer in range(nstations):
            di12 = self.i1[ilayer] * self.i2[ilayer] - self.i12[ilayer] ** 2
            if not di12 > 0.:
                msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End B\n'
                msg += 'pid=%s xxb=%s i1=%s i2=%s i12=%s i1*i2-i12^2=%s'  % (
                    self.pid, self.xxb[ilayer], self.i1[ilayer], self.i2[ilayer],
                    self.i12[ilayer], di12)
                raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBEAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'property_id')
        mid = integer(card, 2, 'material_id')

        area0 = double(card, 3, 'Area')
        i1a = double_or_blank(card, 4, 'I1', 0.0)
        i2a = double_or_blank(card, 5, 'I2', 0.0)
        i12a = double_or_blank(card, 6, 'I12', 0.0)
        ja = double_or_blank(card, 7, 'J', 0.0)
        nsma = double_or_blank(card, 8, 'nsm', 0.0)
        area = [area0]
        i1 = [i1a]
        i2 = [i2a]
        i12 = [i12a]
        j = [ja]
        nsm = [nsma]

        assert area[0] >= 0., 'PBEAM pid=%s area=%s' % (pid, area)
        assert i1[0] >= 0., 'PBEAM pid=%s i1=%s' % (pid, i1)
        assert i2[0] >= 0., 'PBEAM pid=%s i2=%s' % (pid, i2)

        # we'll do a check for warping later; cwa/cwb -> j > 0.0
        assert j[0] >= 0., 'PBEAM pid=%s j=%s' % (pid, j)

        if i1a * i2a - i12a ** 2 <= 0.:
            msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End A\n'
            msg += 'i1=%s i2=%s i12=%s i1*i2-i12^2=%s'  % (i1a, i2a, i12a, i1a*i2a-i12a**2)
            raise ValueError(msg)

        # TODO: can you have a single lined PBEAM...I think so...
        # the second line is blank, so all values would be None (End A)
        # the NO xxb would be implicitly defined as 1.0
        # the End B values would try to use End A values, but because they're not set,
        # the defaults would get applied
        # the final 2 lines will default
        # finally, there would be no output at End A, but there would be output at End A.
        ifield = 9
        field9 = double_string_or_blank(card, 9, 'field9', 0.0)
        if isinstance(field9, float):
            # C/D/E/F
            c1a = double_or_blank(card, 9, 'c1', 0.0)
            c2a = double_or_blank(card, 10, 'c2', 0.0)
            d1a = double_or_blank(card, 11, 'd1', 0.0)
            d2a = double_or_blank(card, 12, 'd2', 0.0)
            e1a = double_or_blank(card, 13, 'e1', 0.0)
            e2a = double_or_blank(card, 14, 'e2', 0.0)
            f1a = double_or_blank(card, 15, 'f1', 0.0)
            f2a = double_or_blank(card, 16, 'f2', 0.0)
            c1 = [c1a]
            c2 = [c2a]
            d1 = [d1a]
            d2 = [d2a]
            e1 = [e1a]
            e2 = [e2a]
            f1 = [f1a]
            f2 = [f2a]
            so = ['YES']
            ifield += 8 # 9 + 8 = 17
        else:
            c1a = c2a = d1a = d2a = e1a = e2a = f1a = f2a = 0.0
            c1 = [None]
            c2 = [None]
            d1 = [None]
            d2 = [None]
            e1 = [None]
            e2 = [None]
            f1 = [None]
            f2 = [None]
            so = ['NO']
            if field9 not in ['YES', 'YESA', 'NO']:
                msg = ('field9=%r on the PBEAM pid=%s must be [YES, YESA, NO] '
                       'because C/D/E/F at A is not specified' % field9)
                raise ValueError(field9)
        xxb = [0.]

        irow = 0
        nrows_max = 10
        for irow in range(nrows_max):
            nrepeated = irow + 1
            SOi_k1 = double_string_or_blank(card, ifield, 'SO_%i/K1' % nrepeated)
            if isinstance(SOi_k1, float) or SOi_k1 is None:
                # we found K1
                break
            else:
                soi = string(card, ifield, 'SO%i' % nrepeated)
                xxbi = double(card, ifield + 1, 'x/xb%i' % nrepeated)
                if xxbi == 1.0:
                    # these have already been checked such that they're greater than 0
                    # so when we interpolate, our values will be correct
                    areai = double_or_blank(card, ifield + 2, 'Area%i' % nrepeated, area0)
                    i1i = double_or_blank(card, ifield + 3, 'I1 %i' % nrepeated, i1a)
                    i2i = double_or_blank(card, ifield + 4, 'I2 %i' % nrepeated, i2a)
                    i12i = double_or_blank(card, ifield + 5, 'I12 %i' % nrepeated, i12a)

                    assert area[-1] >= 0., area
                    assert i1[-1] >= 0., i1
                    assert i2[-1] >= 0., i2

                    # we'll do a check for warping later; cwa/cwb -> j > 0.0
                    assert j[-1] >= 0., j
                    if i1i * i2i - i12i ** 2 <= 0.:
                        msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End B\n'
                        msg += 'xxb=1.0 i1=%s i2=%s i12=%s'  % (i1i, i2i, i12i)
                        raise ValueError(msg)

                    ji = double_or_blank(card, ifield + 6, 'J%i' % nrepeated, ja)
                    nsmi = double_or_blank(card, ifield + 7, 'nsm%i' % nrepeated, nsma)
                else:
                    # we'll go through and do linear interpolation afterwards
                    areai = double_or_blank(card, ifield + 2, 'Area%i' % nrepeated, 0.0)
                    i1i = double_or_blank(card, ifield + 3, 'I1 %i' % nrepeated, 0.0)
                    i2i = double_or_blank(card, ifield + 4, 'I2 %i' % nrepeated, 0.0)
                    i12i = double_or_blank(card, ifield + 5, 'I12 %i' % nrepeated, 0.0)
                    ji = double_or_blank(card, ifield + 6, 'J%i' % nrepeated, 0.0)
                    nsmi = double_or_blank(card, ifield + 7, 'nsm%i' % nrepeated, 0.0)

                so.append(soi)
                xxb.append(xxbi)
                area.append(areai)
                i1.append(i1i)
                i2.append(i2i)
                i12.append(i12i)
                j.append(ji)
                nsm.append(nsmi)

                if soi == 'YES':
                    c1i = double_or_blank(card, ifield + 8, 'c1 %i' % nrepeated, 0.0)
                    c2i = double_or_blank(card, ifield + 9, 'c2 %i' % nrepeated, 0.0)
                    d1i = double_or_blank(card, ifield + 10, 'd1 %i' % nrepeated, 0.0)
                    d2i = double_or_blank(card, ifield + 11, 'd2 %i' % nrepeated, 0.0)
                    e1i = double_or_blank(card, ifield + 12, 'e1 %i' % nrepeated, 0.0)
                    e2i = double_or_blank(card, ifield + 13, 'e2 %i' % nrepeated, 0.0)
                    f1i = double_or_blank(card, ifield + 14, 'f1 %i' % nrepeated, 0.0)
                    f2i = double_or_blank(card, ifield + 15, 'f2 %i' % nrepeated, 0.0)
                    ifield += 16
                elif soi == 'YESA':
                    c1i = c1a
                    c2i = c2a
                    d1i = d1a
                    d2i = d2a
                    e1i = e1a
                    e2i = e2a
                    f1i = f1a
                    f2i = f2a
                    ifield += 8
                elif soi == 'NO':
                    c1i = None
                    c2i = None
                    d1i = None
                    d2i = None
                    e1i = None
                    e2i = None
                    f1i = None
                    f2i = None
                    ifield += 8
                else:
                    raise RuntimeError('so=%r and not [YES, YESA, NO]' % soi)
                c1.append(c1i)
                c2.append(c2i)
                d1.append(d1i)
                d2.append(d2i)
                e1.append(e1i)
                e2.append(e2i)
                f1.append(f1i)
                f2.append(f2i)
        if irow != 0:
            assert min(xxb) == 0.0, 'pid=%s x/xb=%s' % (pid, xxb)
            assert max(xxb) == 1.0, 'pid=%s x/xb=%s' % (pid, xxb)
            assert len(xxb) == len(unique(xxb)), xxb

        # calculate:
        #    k1, k2, s1, s2
        #    m1a, m2a, n1a, n2a, etc.

        # footer fields
        #: Shear stiffness factor K in K*A*G for plane 1.
        k1 = double_or_blank(card, ifield, 'k1', 1.0)
        #: Shear stiffness factor K in K*A*G for plane 2.
        k2 = double_or_blank(card, ifield + 1, 'k2', 1.0)

        #: Shear relief coefficient due to taper for plane 1.
        s1 = double_or_blank(card, ifield + 2, 's1', 0.0)
        #: Shear relief coefficient due to taper for plane 2.
        s2 = double_or_blank(card, ifield + 3, 's2', 0.0)

        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point A.
        nsia = double_or_blank(card, ifield + 4, 'nsia', 0.0)
        #: non structural mass moment of inertia per unit length
        #: about nsm center of gravity at Point B.
        nsib = double_or_blank(card, ifield + 5, 'nsib', nsia)

        #: warping coefficient for end A.
        cwa = double_or_blank(card, ifield + 6, 'cwa', 0.0)
        #: warping coefficient for end B.
        cwb = double_or_blank(card, ifield + 7, 'cwb', cwa)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end A.
        m1a = double_or_blank(card, ifield + 8, 'm1a', 0.0)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end A.
        m2a = double_or_blank(card, ifield + 9, 'm2a', 0.0)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end B.
        m1b = double_or_blank(card, ifield + 10, 'm1b', m1a)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end B.
        m2b = double_or_blank(card, ifield + 11, 'm2b', m2a)

        #: y coordinate of neutral axis for end A.
        n1a = double_or_blank(card, ifield + 12, 'n1a', 0.0)
        #: z coordinate of neutral axis for end A.
        n2a = double_or_blank(card, ifield + 13, 'n2a', 0.0)

        #: y coordinate of neutral axis for end B.
        n1b = double_or_blank(card, ifield + 14, 'n1a', n1a)
        #: z coordinate of neutral axis for end B.
        n2b = double_or_blank(card, ifield + 15, 'n2b', n2a)


        ifield += 16
        if len(card) > ifield:
            msg = 'len(card)=%s is too long; max=%s\n' % (len(card), ifield)
            msg += 'You probably have a empty line after the YESA/NO line.\n'
            msg += 'The next line must have K1.\n'
            msg += 'pid = %s\n' % pid
            msg += 'mid = %s\n' % mid
            msg += 's0 = %s\n' % so
            msg += 'xxb = %s\n' % xxb

            msg += 'A = %s\n' % area
            msg += 'i1 = %s\n' % i1
            msg += 'i2 = %s\n' % i2
            msg += 'i12 = %s\n' % i12
            msg += 'j = %s\n' % j
            msg += 'nsm = %s\n\n' % nsm

            msg += 'c1 = %s\n' % c1
            msg += 'c2 = %s\n' % c2
            msg += 'd1 = %s\n' % d1
            msg += 'd2 = %s\n' % d2
            msg += 'e1 = %s\n' % e1
            msg += 'e2 = %s\n' % e2
            msg += 'f1 = %s\n' % f1
            msg += 'f2 = %s\n\n' % f2

            msg += 'k1 = %s\n' % k1
            msg += 'k2 = %s\n' % k2
            msg += 's1 = %s\n' % s1
            msg += 's2 = %s\n' % s2
            msg += 'nsia = %s\n' % nsia
            msg += 'nsib = %s\n\n' % nsib

            msg += 'cwa = %s\n' % cwa
            msg += 'cwb = %s\n' % cwb
            msg += 'm1a = %s\n' % m1a
            msg += 'm2a = %s\n' % m2a
            msg += 'mb1 = %s\n' % m1b
            msg += 'm2b = %s\n' % m2b
            msg += 'n1a = %s\n' % n1a
            msg += 'n2a = %s\n' % n2a
            msg += 'n1b = %s\n' % n1b
            msg += 'n2b = %s\n' % n2b
            raise RuntimeError(msg)

        return PBEAM(
            pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
            c1, c2, d1, d2, e1, e2, f1, f2,
            k1, k2, s1, s2,
            nsia, nsib, cwa, cwb, m1a,
            m2a, m1b, m2b, n1a, n2a, n1b, n2b,
            comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PBEAM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        (pid, mid, unused_nsegs, unused_ccf, unused_x) = data[:5]

        rows = data[5:]
        if len(rows) != 12:
            msg = 'PBEAM: len(rows)=%s expected=13\n' % len(rows)
            for datai in data:
                msg += '    %s\n' % str(datai)
            raise SyntaxError(msg)
        area = []
        so = []
        xxb = []
        i1 = []
        i2 = []
        i12 = []
        j = []
        nsm = []
        c1 = []
        c2 = []
        d1 = []
        d2 = []
        e1 = []
        e2 = []
        f1 = []
        f2 = []
        ## TODO: PBEAM:op2 handle repeated x/xb = 0.0, 1.0

        iis = []
        for i, pack in enumerate(rows[:-1]):
            (soi, xxbi, areai, i1i, i2i, i12i, ji, nsmi, c1i, c2i,
             d1i, d2i, e1i, e2i, f1i, f2i) = pack
            if i > 0 and allclose(xxbi, 0.0):
                #print('PBEAM - skipping i=%s x/xb=%s' % (i, xxbi))
                continue
            if i > 0 and i != 10 and allclose(xxbi, 1.0):
                #print('PBEAM - skipping i=%s x/xb=%s' % (i, xxbi))
                continue
            area.append(areai)
            xxb.append(xxbi)
            iis.append(i)
            so.append(soi)
            i1.append(i1i)
            i2.append(i2i)
            i12.append(i12i)
            j.append(ji)
            nsm.append(nsmi)
            c1.append(c1i)
            c2.append(c2i)
            d1.append(d1i)
            d2.append(d2i)
            e1.append(e1i)
            e2.append(e2i)
            f1.append(f1i)
            f2.append(f2i)
        #print('*i = %s' % iis)
        #print('*xxb = %s' % xxb)
        #print('*i1 = %s' % i1)

        (k1, k2, s1, s2, nsia, nsib, cwa, unused_cwb,
         m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b) = data[-1]
        return PBEAM(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                     c1, c2, d1, d2, e1, e2, f1, f2,
                     k1=k1, k2=k2, s1=s1, s2=s2,
                     nsia=nsia, nsib=nsib, cwa=cwa, cwb=None,
                     m1a=m1a, m2a=m2a, m1b=m1b, m2b=m2b,
                     n1a=n1a, n2a=n2a, n1b=n1b, n2b=n2b,
                     comment=comment)

    def set_optimization_value(self, name_str, value):
        if name_str == 'I1(A)':
            self.i1[0] = value
        elif name_str == 'I1(B)':
            self.i1[-1] = value

        elif name_str == 'I2(A)':
            self.i1[-1] = value
        elif name_str == 'I2(B)':
            self.i2[-1] = value
        else:
            raise NotImplementedError(name_str)

    def get_optimization_value(self, name_str):
        if name_str == 'I1(A)':
            out = self.i1[0]
        elif name_str == 'I1(B)':
            out = self.i1[-1]

        elif name_str == 'I2(A)':
            out = self.i2[-1]
        elif name_str == 'I2(B)':
            out = self.i2[-1]
        else:
            raise NotImplementedError(name_str)
        return out

    #def Area(self):
       #""".. warning:: area field not supported fully on PBEAM card"""
       ##raise RuntimeError(self.A[0])
       #return self.A[0]

    #def Nsm(self):
       #""".. warning:: nsm field not supported fully on PBEAM card"""
       ##raise RuntimeError(self.nsm[0])
       #return self.nsm[0]

    def I1_I2_I12(self):
        #assert self.i1  is not None, 'I1=%r' % self.i1
        #assert self.i2  is not None, 'I2=%r' % self.i2
        #assert self.i12 is not None, 'I12=%r' % self.i12
        return self.i1[0], self.i2[0], self.i12[0]

    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho = self.Rho()
        mass_per_lengths = []
        for (area, nsm) in zip(self.A, self.nsm):
            #print('area=%s rho=%s nsm=%s' % (area, rho, nsm))
            mass_per_lengths.append(area * rho + nsm)
        mass_per_length = integrate_positive_unit_line(self.xxb, mass_per_lengths)
        return mass_per_length

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBEAM mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)
        #if model.sol != 600:
            #assert max(self.j) == 0.0, self.j
            #assert min(self.j) == 0.0, self.j

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def _verify(self, xref):
        pid = self.Pid()
        mid = self.Mid()
        A = self.Area()
        try:
            J = self.J()
        except NotImplementedError:
            msg = "J is not implemented for pid.type=%s pid.Type=%s" % (self.type, self.beam_type)
            print(msg)
            J = 0.0
        nsm = self.Nsm()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'pid=%r' % A
        assert isinstance(J, float), 'cid=%r' % J
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        if xref:
            assert self.mid_ref.type in ['MAT1', 'MAT4', 'MAT5'], 'pid.type=%s; mid_ref.type=%s' % (
                self.type, self.mid_ref.type)
            #self.MassPerLength()

    def raw_fields(self):
        list_fields = ['PBEAM', self.pid, self.Mid()]

        i = 0
        loop_vars = (
            self.so, self.xxb, self.A, self.i1, self.i2, self.i12, self.j, self.nsm,
            self.c1, self.c2, self.d1, self.d2, self.e1, self.e2, self.f1, self.f2
        )
        for out in zip(*loop_vars):
            (so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2) = out

            if i == 0:  # the first 2 fields aren't written
                list_fields += [A, i1, i2, i12, j, nsm]
                if any([isinstance(cdefi, float) for cdefi in [c1, c2, d1, d2, e1, e2, f1, f2]]):
                    list_fields += [c1, c2, d1, d2, e1, e2, f1, f2]
            else:
                if so in ['YES']:
                    list_fields += [
                        'YES', xxb, A, i1, i2, i12, j, nsm,
                        c1, c2, d1, d2, e1, e2, f1, f2
                    ]
                elif so in ['NO']:
                    list_fields += ['NO', xxb, A, i1, i2, i12, j, nsm]
                elif so in ['YESA']:
                    list_fields += ['YESA', xxb, A, i1, i2, i12, j, nsm]
                else:
                    raise RuntimeError('so=%r type(so)=%s' % (so, type(so)))
            i += 1

        footer = [
            self.k1, self.k2, self.s1, self.s2, self.nsia, self.nsib, self.cwa, self.cwb,
            self.m1a, self.m2a, self.m1b, self.m2b, self.n1a, self.n2a, self.n1b, self.n2b
        ]
        list_fields += footer
        return list_fields

    def repr_fields(self):
        list_fields = ['PBEAM', self.pid, self.Mid()]
        #print("c1=%r c2=%r" % (self.c1, self.c2))
        #print("d1=%r d2=%r" % (self.d1, self.d2))
        #print("e1=%r e2=%r" % (self.e1, self.e2))
        #print("f1=%r f2=%r" % (self.f1, self.f2))
        #print('i1  = %r' % self.i1)
        #print('i2  = %r' % self.i2)
        #print('i12 = %r' % self.i12)
        i = 0
        for (so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1,
             f2) in zip(self.so, self.xxb, self.A, self.i1, self.i2, self.i12,
                        self.j, self.nsm, self.c1, self.c2, self.d1, self.d2,
                        self.e1, self.e2, self.f1, self.f2):
            i1 = set_blank_if_default(i1, 0.0)
            i2 = set_blank_if_default(i2, 0.0)
            i12 = set_blank_if_default(i12, 0.0)
            j = set_blank_if_default(j, 0.0)

            nsm = set_blank_if_default(nsm, 0.0)
            c1 = set_blank_if_default(c1, 0.0)
            d1 = set_blank_if_default(d1, 0.0)
            e1 = set_blank_if_default(e1, 0.0)
            f1 = set_blank_if_default(f1, 0.0)

            c2 = set_blank_if_default(c2, 0.0)
            d2 = set_blank_if_default(d2, 0.0)
            e2 = set_blank_if_default(e2, 0.0)
            f2 = set_blank_if_default(f2, 0.0)

            if i == 0:  # the first 2 fields aren't written
                list_fields += [A, i1, i2, i12, j, nsm,
                                c1, c2, d1, d2, e1, e2, f1, f2]
            else:
                if so in ['YES']:
                    list_fields += ['YES', xxb, A, i1, i2, i12, j, nsm,
                                    c1, c2, d1, d2, e1, e2, f1, f2]
                elif so in ['NO']:
                    list_fields += ['NO', xxb, A, i1, i2, i12, j, nsm]
                elif so in ['YESA']:
                    list_fields += ['YESA', xxb, A, i1, i2, i12, j, nsm]
                else:
                    raise RuntimeError('so=%r type(so)=%s' % (so, type(so)))

            i += 1
        k1 = set_blank_if_default(self.k1, 1.0)
        k2 = set_blank_if_default(self.k2, 1.0)

        s1 = set_blank_if_default(self.s1, 0.0)
        s2 = set_blank_if_default(self.s2, 0.0)
        #k1 = self.k1
        #k2 = self.k2
        #s1 = self.s1
        #s2 = self.s2

        nsia = set_blank_if_default(self.nsia, 0.0)
        nsib = set_blank_if_default(self.nsib, self.nsia)

        cwa = set_blank_if_default(self.cwa, 0.0)
        cwb = set_blank_if_default(self.cwb, self.cwa)

        #m1a = self.m1a
        #m2a = self.m2a
        #m1b = self.m1b
        #m2b = self.m2b
        # Point A/B
        # Directions 1/2
        m1a = set_blank_if_default(self.m1a, 0.0)
        m2a = set_blank_if_default(self.m2a, 0.0)
        m1b = set_blank_if_default(self.m1b, self.m1a)
        m2b = set_blank_if_default(self.m2b, self.m2a)

        n1a = set_blank_if_default(self.n1a, 0.0)
        n2a = set_blank_if_default(self.n2a, 0.0)
        n1b = set_blank_if_default(self.n1b, self.n1a)
        n2b = set_blank_if_default(self.n2b, self.n2b)

        footer = [k1, k2, s1, s2, nsia, nsib, cwa, cwb,
                  m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b]

        list_fields += footer
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)

def update_pbeam_negative_integer(pname_fid):
    """
    Converts the negative PBEAM value to a positive one

    Parameters
    ----------
    pname_fid : int
        for a PBEAM this should be between [-5, -167]
        the negative values correspond to the numbers in the MPT OP2 table

    Returns
    -------
    pname_fid : str
        a pname is of 'J(3)' is far more clear than -37

    TODO: only handles istation=0 for now (e.g., 'J(1)')
    """
    # shift to divisible by 16
    if not (-167 <= pname_fid <= -6):  # pragma: no cover
        msg = "A-property_type='PBEAM' has not implemented %r in pname_map" % (
            pname_fid)
        raise NotImplementedError(msg)
    ioffset = -pname_fid - 6
    istation = ioffset // 16
    iterm = ioffset % 16

    # 0    1   2   3   4   5   6   7    8  9
    #(soi, xxb, a, i1, i2, i12, j, nsm, c1, c2,
    #
     #10  11  12  13  14  15
     #d1, d2, e1, e2, f1, f2) = pack
    assert istation == 0, istation
    if iterm == 2:
        word = 'A'
    elif iterm == 3:
        word = 'I1'
    elif iterm == 4:
        word = 'I2'
    elif iterm == 5:
        word = 'I12'
    elif iterm == 6:
        word = 'J'
    #elif iterm == 7:
        #self.nsm[istation] = value # 13-6 = 6
    elif iterm == 8:
        word = 'C1'
    elif iterm == 9:
        word = 'C2'
    elif iterm == 10:
        word = 'D1'
    elif iterm == 11:
        word = 'D2'
    elif iterm == 12:
        word = 'E1'
    elif iterm == 13:
        word = 'E2'
    elif iterm == 14:
        word = 'F1'
    elif iterm == 15:
        word = 'F2'
    else:  # pragma: no cover
        #print('istation=%s iterm=%s' % (istation, iterm))
        msg = ("property_type='PBEAM' has not implemented %r (istation=%r, iterm=%r)"
               " in pname_map" % (pname_fid, istation, iterm))
        raise NotImplementedError(msg)
    word = '%s(%i)' % (word, istation + 1)
    return word

class PBEAML(IntegratedLineProperty):
    """
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |   1    |    2    |    3    |    4    |    5    |    6    |    7    |    8    |    9    |
    +========+=========+=========+=========+=========+=========+=========+=========+=========+
    | PBEAML |   PID   |   MID   |  GROUP  |  TYPE   |         |         |         |         |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | DIM1(A) | DIM2(A) |  etc.   | DIMn(A) | NSM(A)  | SO(1)   | X(1)/XB | DIM1(1) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | DIM2(1) |  etc.   | DIMn(1) | NSM(1)  | SO(2)   | X(2)/XB | DIM1(2) | DIM2(2) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        |  etc.   | DIMn(2) | NSM(m)  |  etc.   | SO(m)   | X(m)/XB | DIM1(m) |  etc.   |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | DIMn(m) |  NSM(m) | SO(B)   |   1.0   | DIM1(B) | DIM2(B) |  etc.   | DIMn(B) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | NSM(B)  |         |         |         |         |         |         |         |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    """
    type = 'PBEAML'
    _properties = ['valid_types', 'Type']
    valid_types = {
        "ROD": 1,
        "TUBE": 2,
        "TUBE2": 2,
        "L": 4,
        "I": 6,
        "CHAN": 4,
        "T": 4,
        "BOX": 4,
        "BAR": 2,
        "CROSS": 4,
        "H": 4,
        "T1": 4,
        "I1": 4,
        "CHAN1": 4,
        "Z": 4,
        "CHAN2": 4,
        "T2": 4,
        "BOX1": 6,
        "HEXA": 3,
        "HAT": 4,
        "HAT1": 5,
        "DBOX": 10,  # TODO: was 12???
    }  # for GROUP="MSCBML0"
    def update_by_pname_fid(self, pname_fid, value):
        if isinstance(pname_fid, int):
            raise NotImplementedError('property_type=%r has not implemented %r in pname_map' % (
                self.type, pname_fid))
        elif pname_fid.startswith('DIM'):
            # DIM1, DIM1(A), DIM1(10), DIM2(B)
            #
            # (DIM, num, station)
            #
            station_str = 'A' #  the default, A
            if '(' in pname_fid:
                assert pname_fid.endswith(')'), pname_fid
                pname_fid, station_str = pname_fid[:-1].split('(', 1)

            ndim = len(self.dim[0])
            if station_str == 'A':
                istation = 0
            elif station_str == 'B':
                istation = ndim - 1
            else:
                istation = int(station_str) - 1

            unused_dim_word, num = break_word_by_trailing_integer(pname_fid)

            idim = int(num) - 1

            # in Nastran syntax
            #print('DIM%i(%i)' % (idim+1, istation+1))
            try:
                dim_station = self.dim[istation]
                dim_station[idim] = value
            except:
                print('pname_fid=%r num=%r ndim=%r' % (pname_fid, num, ndim))
                print('istation=%r idim=%r' % (istation, idim))
                print(self)
                raise
        elif isinstance(pname_fid, str):
            if not pname_fid[-1].isdigit():
                param_name = pname_fid
                idim = -1
            else:
                param_name, num = break_word_by_trailing_integer(pname_fid)
                idim = int(num) - 1

            if param_name == 'NSM':
                self.nsm[idim] = value
            else:
                raise NotImplementedError('property_type=%r param_name=%r idim=%s '
                                          'has not implemented %r in pname_map' % (
                                              self.type, param_name, idim, pname_fid))
        else:
            raise NotImplementedError('property_type=%r has not implemented %r in pname_map' % (
                self.type, pname_fid))

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        beam_type = 'ROD'
        xxb = [1.]
        dims = [[1.]]
        return PBEAML(pid, mid, beam_type, xxb, dims,
                      so=None, nsm=None, group='MSCBML0', comment='')

    def __init__(self, pid, mid, beam_type, xxb, dims, so=None, nsm=None,
                 group='MSCBML0', comment=''):
        """
        Creates a PBEAML card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        beam_type : str
            the section profile
        xxb : List[float]
            The percentage locations along the beam [0., ..., 1.]
        dims : List[dim]
            dim : List[float]
                The dimensions for each section
        so : List[str]; default=None
            YES, YESA, NO
            None : [0.] * len(xxb)
        nsm : List[float]; default=None
            nonstructural mass per unit length
            None : [0.] * len(xxb)
        group : str; default='MSCBML0'
            this parameter can lead to a very broken deck with a very
            bad error message; don't touch it!
        comment : str; default=''
            a comment for the card
        """
        IntegratedLineProperty.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.group = group
        #: Section Type (e.g. 'ROD', 'TUBE', 'I', 'H')
        self.beam_type = beam_type
        ndim = self.valid_types[self.beam_type]

        nxxb = len(xxb)
        if nxxb == 0:
            raise IndexError('pid=%s; len(xxb)=0; at least 1 station must be defined' % pid)
        if nsm is None:
            nsm = [0.] * nxxb
        elif not isinstance(nsm, (list, tuple, ndarray)):
            msg = 'pid=%s; nsm=%s and must be a list/tuple/ndarray; type=%s' % (
                pid, nsm, type(nsm))
            raise TypeError(msg)

        if so is None:
            so = ['YES'] * nxxb
        elif not isinstance(so, (list, tuple, ndarray)):
            msg = 'pid=%s; so=%s and must be a list/tuple/ndarray; type=%s' % (
                pid, so, type(so))
            raise TypeError(msg)

        for istation, xxbi, nsmi, dim in zip(count(), xxb, nsm, dims):
            if not isinstance(dim, (list, ndarray)):
                msg = 'dims = List[dim]; dim=List[floats]; type(dim)=%s' % (type(dim))
                raise TypeError(msg)
            assert len(dim) == ndim, 'beam_type=%s ndim=%s len(dim)=%s xxb=%s dim=%s' % (
                beam_type, ndim, len(dim), xxbi, dim)

            if not isinstance(xxbi, float_types):
                raise TypeError('istation=%i xxb=%s and must be a float' % (
                    istation, xxbi))
            if not isinstance(nsmi, float_types):
                raise TypeError('istation=%i nsm=%s and must be a float' % (
                    istation, nsmi))

            for idim, dimi in enumerate(dim):
                if not isinstance(dimi, float_types):
                    raise TypeError('istation=%i dim%i=%s and must be a float' % (
                        istation, idim+1, dimi))

        self.dim = np.asarray(dims)
        self.xxb = np.asarray(xxb)
        self.so = so
        self.nsm = np.asarray(nsm)
        self.mid_ref = None
        unused_A = self.Area()

    def validate(self):
        uxxb = np.unique(self.xxb)
        if len(self.xxb) != len(uxxb):
            raise ValueError('xxb=%s unique(xxb)=%s' % (self.xxb, uxxb))

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        if isinstance(self.dim, list):
            self.dim = np.asarray(self.dim)
            self.xxb = np.asarray(self.xxb)
            self.nsm = np.asarray(self.nsm)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBEAML card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        group = string_or_blank(card, 3, 'group', 'MSCBML0')
        beam_type = string(card, 4, 'Type')

        # determine the number of required dimensions on the PBEAM
        ndim = cls.valid_types[beam_type]

        #: dimension list
        dims = []
        dim = []

        #: Section position
        xxb = [0.]

        #: Output flag
        so = ['YES']  # station 0

        #: non-structural mass :math:`nsm`
        nsm = []

        i = 9
        n = 0

        #n_so = (len(card) - 9) // (ndim + 2) #- 1
        #n_extra = (len(card) - 9) % (ndim + 2)
        xxbi = 0.0
        while i < len(card):
            if n > 0:
                soi = string_or_blank(card, i, 'so_n=%i' % n, 'YES')
                xxbi = double_or_blank(card, i + 1, 'xxb_n=%i' % n, 1.0)
                so.append(soi)
                xxb.append(xxbi)
                i += 2

            # PBARL
            # 9. For DBOX section, the default value for DIM5 to DIM10 are
            #    based on the following rules:
            #     a. DIM5, DIM6, DIM7 and DIM8 have a default value of
            #        DIM4if not provided.
            #     b. DIM9 and DIM10 have a default value of DIM6 if not
            #        provided.

            #If any of the fields NSM(B), DIMi(B) are blank on the
            #continuation entry for End B, the values are set to the
            #values given for end A. For the continuation entries that
            #have values of X(j)/XB between 0.0 and 1.0 and use the
            #default option (blank field), a linear interpolation between
            #the values at ends A and B is performed to obtain the
            #missing field.
            dim = []
            if beam_type == 'DBOX':
                for ii in range(ndim):
                    field_name = 'istation=%s; ndim=%s; dim%i' % (n, ndim, ii+1)
                    if ii in [4, 5, 6, 7]:
                        dim4 = dim[3]
                        dimi = double_or_blank(card, i, field_name, default=dim4)
                    elif ii in [8, 9]:
                        dim6 = dim[5]
                        dimi = double_or_blank(card, i, field_name, default=dim6)
                    else:
                        dimi = double(card, i, field_name)
                    dim.append(dimi)
                    i += 1
            else:
                for ii in range(ndim):
                    field_name = 'istation=%s; ndim=%s; dim%i' % (n, ndim, ii+1)
                    if xxbi == 0.0:
                        dimi = double(card, i, field_name)
                    elif xxbi == 1.0:
                        dims0 = dims[0]
                        dimi = double_or_blank(card, i, field_name, dims0[ii])
                    else:
                        ## TODO: use linear interpolation
                        dimi = double(card, i, field_name)

                    dim.append(dimi)
                    i += 1
            dims.append(dim)

            nsmi = double_or_blank(card, i, 'nsm_n=%i' % n, 0.0)
            nsm.append(nsmi)
            n += 1
            i += 1
        assert len(card) > 5, card
        return PBEAML(pid, mid, beam_type, xxb, dims, group=group,
                      so=so, nsm=nsm, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        nsm = self.Nsm()
        area = self.Area()
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))
        assert isinstance(nsm, float), 'nsm=%r\n%s' % (nsm, str(self))
        assert isinstance(area, float), 'area=%r\n%s' % (area, str(self))
        if xref:
            rho = self.Rho()
            mass_per_length = self.MassPerLength()
            assert isinstance(rho, float), 'rho=%r\n%s' % (rho, str(self))
            assert isinstance(mass_per_length, float), 'mass/L=%r\n%s' % (mass_per_length, str(self))

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PBEAML card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        TODO: this doesn't work right for the calculation of area
              the card is all messed up
        """
        (pid, mid, group, beam_type, fvalues) = data
        group = group.strip()
        beam_type = beam_type.strip()
        ndim = cls.valid_types[beam_type]
        nfvalues = len(fvalues)
        nsections = nfvalues // (3 + ndim)
        sections = fvalues.reshape(nsections, ndim+3)
        #print('sections = \n%s' % sections)

        xxb = []
        so = []
        dims = []
        nsm = []
        # XXB, SO, NSM, and dimensions
        isections = []
        for i, section in enumerate(sections):
            xxbi = section[1]
            #print('PBEAML - i=%s x/xb=%s' % (i, xxbi))
            if i > 0 and allclose(xxbi, 0.0):
                #print('  PBEAML - skipping i=%s x/xb=%s' % (i, xxbi))
                continue
            if xxbi in xxb:
                #print("  PBEAML - skipping i=%s x/xb=%s because it's a duplicate" % (i, xxbi))
                continue
            isections.append(i)

            so_float = section[0]
            if so_float == 0.:
                so_string = 'YES'
            elif so_float == 1.:
                so_string = 'NO'
            else:
                msg = 'so_float=%r; expected 0.0 or 1.0; data=%s' % (so_float, data)
                raise NotImplementedError(msg)

            dim = list(section[2:-1])
            nsmi = section[-1]
            #print(dim, section)
            xxb.append(xxbi)
            so.append(so_string)
            dims.append(dim)
            nsm.append(nsmi)
        #print('sections2 = \n%s' % sections[isections])
        return PBEAML(pid, mid, beam_type, xxb, dims, group=group,
                      so=so, nsm=nsm, comment=comment)

    @property
    def Type(self):
        """gets Type"""
        return self.beam_type
    @Type.setter
    def Type(self, beam_type):
        """sets Type"""
        self.beam_type = beam_type

    def get_mass_per_lengths(self):
        """
        helper method for MassPerLength
        a*rho + nsm
        """
        rho = self.Rho()
        mass_per_lengths = []
        for (dim, nsm) in zip(self.dim, self.nsm):
            a = _bar_areaL('PBEAML', self.beam_type, dim, self)
            try:
                mass_per_lengths.append(a * rho + nsm)
            except:
                msg = "PBEAML a*rho+nsm a=%s rho=%s nsm=%s" % (a, rho, nsm)
                raise RuntimeError(msg)
        return mass_per_lengths

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the PBEAML.

        .. math:: \frac{m}{L} = A(x) \rho + nsm

        .. math:: \frac{m}{L} = nsm L + \rho \int \, A(x) dx
        """
        mass_per_lengths = self.get_mass_per_lengths()
        mass_per_length = integrate_positive_unit_line(self.xxb, mass_per_lengths)
        return mass_per_length

    def Area(self):
        r"""
        Gets the Area :math:`A` of the PBEAML.

        .. math:: A = \int \, A(x) dx

        Notes
        -----
        a spline is fit to :math:`A(x)` and then integrated.
        """
        areas = []
        for dim in self.dim:
            areas.append(_bar_areaL('PBEAML', self.beam_type, dim, self))
        try:
            A = integrate_unit_line(self.xxb, areas)
        except ValueError:
            print('PBEAML integration error; pid=%s x/xb=%s areas=%s' % (
                self.pid, self.xxb, areas))
            assert len(self.xxb) == len(areas)
            raise
            #A = mean(areas)
        return A

    #def Mid(self):
    #    return self.mid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. warning:: For structural problems, PBEAML entries must
                     reference a MAT1 material entry
        .. warning:: For heat-transfer problems, the MID must
                     reference a MAT4 or MAT5 material entry.
        .. todo:: What happens when there are 2 subcases?
        """
        msg = ', which is required by PBEAML mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def verify(self, model, isubcase):
        if model.is_thermal_solution(isubcase):
            assert self.mid_ref.type in ['MAT4', 'MAT5']
        else:
            assert self.mid_ref.type in ['MAT1']

    def J(self):
        #j = []
        #for unused_dims in self.dim:
            # calculate J for the station
        #Js = self._J()
        #j = integrate_positive_unit_line(self.xxb, Js)
        j = None
        return j

    def I11(self):
        #i1 = integrate_positive_unit_line(self.xxb,self.i1)
        i1 = None
        return i1

    def I22(self):
        #i2 = integrate_positive_unit_line(self.xxb,self.i2)
        i2 = None
        return i2

    def I12(self):
        #i12 = integrate_unit_line(self.xxb,self.i12)
        i12 = None
        return i12

    def raw_fields(self):
        list_fields = ['PBEAML', self.pid, self.Mid(), self.group, self.beam_type,
                       None, None, None, None]
        #print("xxb=%s so=%s dim=%s nsm=%s" % (
            #self.xxb,self.so, self.dim,self.nsm))
        for (i, xxb, so, dim, nsm) in zip(count(), self.xxb, self.so,
                                          self.dim, self.nsm):
            if i == 0:
                list_fields += dim.tolist() + [nsm]
            else:
                list_fields += [so, xxb] + dim.tolist() + [nsm]
        #raise NotImplementedError('verify PBEAML...')
        return list_fields

    def repr_fields(self):
        group = set_blank_if_default(self.group, 'MSCBML0')
        list_fields = self.raw_fields()
        list_fields[3] = group
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """.. todo:: having bug with PBEAML"""
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)  #is this allowed???


class PBMSECT(LineProperty):
    """
    not done
    """
    type = 'PBMSECT'
    _properties = ['outp_id']

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 2
        form = 'FORM'
        options = [('OUTP', 10)]
        return PBMSECT(pid, mid, form, options, comment='')

    def _finalize_hdf5(self, encoding):
        self.brps = {key : value for key, value in zip(*self.brps)}
        self.ts = {key : value for key, value in zip(*self.ts)}
        self.inps = {key : value for key, value in zip(*self.inps)}
        self.core = {key : value for key, value in zip(*self.core)}

    def __init__(self, pid, mid, form, options, comment=''):
        LineProperty.__init__(self)
        if comment:
            self.comment = comment

        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.form = form

        # integer
        self.outp = None

        # float
        self.nsm = 0.

        # int : int
        self.brps = {}
        self.inps = {}
        self.core = {}

        # int : floats
        self.ts = {}
        assert isinstance(options, list), options
        for key_value in options:
            try:
                key, value = key_value
            except ValueError:
                print(key_value)
                raise
            key = key.upper()
            if key == 'NSM':
                self.nsm = float(value)

            elif 'INP' in key:
                if key.startswith('INP('):
                    assert key.endswith(')'), 'key=%r' % key
                    key_id = int(key[4:-1])
                    self.inps[key_id] = int(value)
                else:
                    if 1 not in self.inps:
                        self.inps[1] = []
                    self.inps[1].append(int(value))

            elif key == 'OUTP':
                self.outp = int(value)

            elif key.startswith('BRP'):
                if key.startswith('BRP('):
                    assert key.endswith(')'), 'key=%r' % key
                    key_id = int(key[4:-1])
                    self.brps[key_id] = int(value)
                else:
                    self.brps[0] = int(value)

            elif key.startswith('T('):
                index, out = split_arbitrary_thickness_section(key, value)
                self.ts[index] = out
            elif key == 'T':
                self.ts[1] = float(value)
            elif key == 'CORE':
                key = 'CORE(1)'
                index, out = split_arbitrary_thickness_section(key, value)
                self.core[index] = out
            elif key.startswith(('CORE(', 'C(')):
                index, out = split_arbitrary_thickness_section(key, value)
                self.core[index] = out
            else:
                raise NotImplementedError('PBMSECT.pid=%s key=%r value=%r' % (pid, key, value))

        assert self.outp is not None, 'options=%s' % str(options)
        self.mid_ref = None
        self.outp_ref = None
        self.brps_ref = {}

    def validate(self):
        assert self.form in ['GS', 'OP', 'CP'], 'pid=%s form=%r' % (self.pid, self.form)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBMSECT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : List[str]
            this card is special and is not a ``BDFCard`` like other cards
        comment : str; default=''
            a comment for the card
        """
        line0 = card[0]
        if '\t' in line0:
            line0 = line0.expandtabs()
        bdf_card = BDFCard(to_fields([line0], 'PBMSECT'))
        #line0_eq = line0[16:]
        lines_joined = ','.join(card[1:]).replace(' ', '').replace(',,', ',')

        if lines_joined:
            fields = get_beam_sections(lines_joined)
            options = [field.split('=', 1) for field in fields]
        else:
            options = []

        pid = integer(bdf_card, 1, 'pid')
        mid = integer(bdf_card, 2, 'mid')
        form = string_or_blank(bdf_card, 3, 'form')

        return PBMSECT(pid, mid, form, options, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):  # pragma: no cover
        #pid = data[0]
        #mid = data[1]
        #group = data[2].strip()
        #Type = data[3].strip()
        #dim = list(data[4:-1])
        #nsm = data[-1]
        #print("group = %r" % group)
        #print("Type  = %r" % Type)
        #print("dim = ",dim)
        #print(str(self))
        #print("*PBMSECT = ", data)
        #raise NotImplementedError('PBMSECT not finished...data=%s' % str(data))
        #return PBMSECT(pid, mid, group, Type, dim, nsm, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBMSECT mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        self.outp_ref = model.Set(self.outp)
        self.outp_ref.cross_reference_set(model, 'Point', msg=msg)

        self.brps_ref = {}
        if len(self.brps):
            for key, brpi in self.brps.items():
                brpi_ref = model.Set(brpi, msg=msg)
                brpi_ref.cross_reference_set(model, 'Point', msg=msg)
                self.brps_ref[key] = brpi_ref

    def plot(self, model, figure_id=1, show=False):
        """
        Plots the beam section

        Parameters
        ----------
        model : BDF()
            the BDF object
        figure_id : int; default=1
            the figure id
        show : bool; default=False
            show the figure when done
        """
        class_name = self.__class__.__name__
        form_map = {
            'GS' : 'General Section',
            'OP' : 'Open Profile',
            'CP' : 'Closed Profile',
        }
        formi = ' form=%s' % form_map[self.form]
        plot_arbitrary_section(
            model, self,
            self.inps, self.ts, self.brps_ref, self.nsm, self.outp_ref,
            figure_id=figure_id,
            title=class_name + ' pid=%s' % self.pid + formi,
            show=show)

    @property
    def outp_id(self):
        if self.outp_ref is not None:
            return self.outp_ref.sid
        return self.outp

    #@property
    #def brp1_id(self):
        #if self.brp1_ref is not None:
            #return self.brp1_ref.sid
        #return self.brp1

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None
        self.outp_ref = None
        self.brps_ref = {}

    def _verify(self, xref):
        pid = self.Pid()
        mid = self.Mid()
        #A = self.Area()
        #J = self.J()
        #nsm = self.Nsm()
        #mpl = self.MassPerLength()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        #assert isinstance(A, float), 'pid=%r' % A
        #assert isinstance(J, float), 'cid=%r' % J
        #assert isinstance(nsm, float), 'nsm=%r' % nsm
        #assert isinstance(mpl, float), 'mass_per_length=%r' % mpl

    def Area(self):
        """
        Gets the area :math:`A` of the CBEAM.
        """
        return 0.
        #raise NotImplementedError('Area is not implemented for PBMSECT')

    def Nsm(self):
        """
        Gets the non-structural mass :math:`nsm` of the CBEAM.
        """
        return 0.
        #raise NotImplementedError('Nsm is not implemented for PBMSECT')

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CBEAM.

        .. math:: \frac{m}{L} = A \rho + nsm
        """
        rho = self.Rho()
        area = self.Area()
        nsm = self.Nsm()
        return area * rho + nsm

    def I11(self):
        raise NotImplementedError('I11 is not implemented for PBMSECT')

    #def I12(self):
        #return self.I12()

    def J(self):
        raise NotImplementedError('J is not implemented for PBMSECT')

    def I22(self):
        raise NotImplementedError('I22 is not implemented for PBMSECT')

    def raw_fields(self):
        """not done..."""
        list_fields = ['PBMSECT', self.pid, self.Mid(), self.form]
        return list_fields

    def repr_fields(self):
        """not done..."""
        list_fields = ['PBMSECT', self.pid, self.Mid(), self.form]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = ['PBMSECT', self.pid, self.Mid(), self.form]
        end = write_arbitrary_beam_section(
            self.inps, self.ts, self.brps, self.nsm, self.outp_id, self.core)
        out = self.comment + print_card_8(card) + end
        return out

    def __repr__(self):
        return self.write_card()


class PBCOMP(LineProperty):
    """
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |   1    |   2  |  3  |  4  |   5  |  6 |  7  |   8    |  9  |
    +========+======+=====+=====+======+====+=====+========+=====+
    | PBCOMP | PID  | MID | A   |  I1  | I2 | I12 |   J    | NSM |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        |  K1  | K2  | M1  |  M2  | N1 | N2  | SYMOPT |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        |  Y1  | Z1  | C1  | MID1 |    |     |        |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        |  Y2  | Z2  | C2  | MID2 |    |     |        |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    |        | ...  | ... | ... |      |    |     |        |     |
    +--------+------+-----+-----+------+----+-----+--------+-----+
    """
    type = 'PBCOMP'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        y = [1.]
        z = [1.]
        c = [1.]
        mids = [1]
        return PBCOMP(pid, mid, y, z, c, mids,
                      area=0.0, i1=0.0, i2=0.0, i12=0.0, j=0.0, nsm=0.0,
                      k1=1.0, k2=1.0, m1=0.0, m2=0.0, n1=0.0, n2=0.0,
                      symopt=0, comment='')

    def __init__(self, pid, mid, y, z, c, mids,
                 area=0.0, i1=0.0, i2=0.0, i12=0.0, j=0.0, nsm=0.0,
                 k1=1.0, k2=1.0, m1=0.0, m2=0.0, n1=0.0, n2=0.0,
                 symopt=0, comment=''):
        """
        Creates a PBCOMP card

        Parameters
        ----------
        pid : int
            Property ID
        mid : int
            Material ID
        mids : List[int]
            Material ID for the i-th integration point
        y / z : List[float]
            The (y,z) coordinates of the lumped areas in the element
            coordinate system
        c : List[float]; default=0.0
            Fraction of the total area for the i-th lumped area
            default not supported...
        area : float
            Area of beam cross section
        i1 / i2 : float; default=0.0
            Area moment of inertia about plane 1/2 about the neutral axis
        i12 : float; default=0.0
           area product of inertia
        j : float; default=0.0
            Torsional moment of interia
        nsm : float; default=0.0
            Nonstructural mass per unit length
        k1 / k2 : float; default=1.0
            Shear stiffness factor K in K*A*G for plane 1/2
        m1 / m2 : float; default=0.0
            The (y,z) coordinates of center of gravity of nonstructural mass
        n1 / n2 : float; default=0.0
            The (y,z) coordinates of neutral axis
        symopt : int; default=0
            Symmetry option to input lumped areas for the beam cross section
            0 < Integer < 5
        comment : str; default=''
            a comment for the card

        """
        LineProperty.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        # Area
        self.A = area
        self.i1 = i1
        self.i2 = i2
        self.i12 = i12
        #: Polar Moment of Inertia :math:`J`
        self.j = j
        #: Non-structural mass per unit length :math:`\frac{m}{L}`
        self.nsm = nsm
        self.k1 = k1
        self.k2 = k2
        self.m1 = m1
        self.m2 = m2
        self.n1 = n1
        self.n2 = n2
        self.symopt = symopt
        self.y = y
        self.z = z
        self.c = c
        self.mids = mids
        assert 0 <= self.symopt <= 5, 'symopt=%i is invalid; ' % self.symopt
        self.mid_ref = None
        self.mids_ref = None

    def validate(self):
        assert isinstance(self.mids, list), 'mids=%r type=%s' % (self.mids, type(self.mids))
        assert isinstance(self.y, list), 'y=%r type=%s' % (self.y, type(self.y))
        assert isinstance(self.z, list), 'z=%r type=%s' % (self.z, type(self.z))
        assert isinstance(self.c, list), 'c=%r type=%s' % (self.c, type(self.c))

        nmids = len(self.mids)
        assert nmids == len(self.y), 'len(mids)=%s len(y)=%s' % (nmids, len(self.y))
        assert nmids == len(self.z), 'len(mids)=%s len(z)=%s' % (nmids, len(self.z))
        assert nmids == len(self.c), 'len(mids)=%s len(c)=%s' % (nmids, len(self.c))
        assert self.symopt in [0, 1, 2, 3, 4, 5], 'symopt=%r' % self.symopt

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBCOMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        area = double_or_blank(card, 3, 'Area', 0.0)
        i1 = double_or_blank(card, 4, 'I1', 0.0)
        i2 = double_or_blank(card, 5, 'I2', 0.0)
        i12 = double_or_blank(card, 6, 'I12', 0.0)
        j = double_or_blank(card, 7, 'J', 0.0)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)
        k1 = double_or_blank(card, 9, 'k1', 1.0)
        k2 = double_or_blank(card, 10, 'k2', 1.0)
        m1 = double_or_blank(card, 11, 'm1', 0.0)
        m2 = double_or_blank(card, 12, 'm2', 0.0)
        n1 = double_or_blank(card, 13, 'n1', 0.0)
        n2 = double_or_blank(card, 14, 'n2', 0.0)
        symopt = integer_or_blank(card, 15, 'symopt', 0)
        y = []
        z = []
        c = []
        mids = []

        nfields = len(card) - 17
        nrows = nfields // 8
        if nfields % 8 > 0:
            nrows += 1

        for row in range(nrows):
            i = 8 * row + 17
            yi = double(card, i, 'y' + str(row))
            zi = double(card, i + 1, 'z' + str(row))
            ci = double_or_blank(card, i + 2, 'c' + str(row), 0.0)
            mid = integer_or_blank(card, i + 3, 'mid' + str(row), mid)
            y.append(yi)
            z.append(zi)
            c.append(ci)
            mids.append(mid)
        return PBCOMP(pid, mid, y, z, c, mids,
                      area, i1, i2, i12, j, nsm,
                      k1, k2, m1, m2, n1, n2,
                      symopt, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        data1, data2 = data
        (pid, mid, area, i1, i2, i12, j, nsm, k1, k2, m1, m2, n1, n2, unused_nsections) = data1
        y = []
        z = []
        c = []
        mids = []

        symopt = 5
        for yi, zi, ci, midi in data2:
            y.append(yi)
            z.append(zi)
            c.append(ci)
            mids.append(midi)
        return PBCOMP(pid, mid, y, z, c, mids,
                      area, i1, i2, i12, j, nsm,
                      k1, k2, m1, m2, n1, n2,
                      symopt, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, int)

    def MassPerLength(self):
        return self.nsm + self.mid_ref.Rho() * self.A

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PBCOMP mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)
        self.mids_ref = model.Materials(self.mids, msg=msg)

    def Mids(self):
        if self.mids_ref is None:
            return self.mids
        return [mid.mid for mid in self.mids_ref]

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mids = self.Mids()
        self.mid_ref = None
        self.mids_ref = None

    def raw_fields(self):
        list_fields = ['PBCOMP', self.pid, self.Mid(), self.A, self.i1,
                       self.i2, self.i12, self.j, self.nsm, self.k1, self.k2,
                       self.m1, self.m2, self.n1, self.n2, self.symopt, None]
        for (yi, zi, ci, mid) in zip(self.y, self.z, self.c, self.Mids()):
            list_fields += [yi, zi, ci, mid, None, None, None, None]
        return list_fields

    def repr_fields(self):
        area = set_blank_if_default(self.A, 0.0)
        j = set_blank_if_default(self.j, 0.0)
        i1 = set_blank_if_default(self.i1, 0.0)
        i2 = set_blank_if_default(self.i2, 0.0)
        i12 = set_blank_if_default(self.i12, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)

        k1 = set_blank_if_default(self.k1, 1.0)
        k2 = set_blank_if_default(self.k2, 1.0)

        m1 = set_blank_if_default(self.m1, 0.0)
        m2 = set_blank_if_default(self.m2, 0.0)

        n1 = set_blank_if_default(self.n1, 0.0)
        n2 = set_blank_if_default(self.n2, 0.0)

        symopt = set_blank_if_default(self.symopt, 0)

        list_fields = ['PBCOMP', self.pid, self.Mid(), area, i1, i2, i12, j,
                       nsm, k1, k2, m1, m2, n1, n2, symopt, None]

        for (yi, zi, ci, mid) in zip(self.y, self.z, self.c, self.mids):
            ci = set_blank_if_default(ci, 0.0)
            list_fields += [yi, zi, ci, mid, None, None, None, None]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
