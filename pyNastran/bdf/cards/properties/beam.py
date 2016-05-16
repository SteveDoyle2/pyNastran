# coding: utf-8
# pylint: disable=C0103
"""
All beam properties are defined in this file.  This includes:
 *   PBEAM
 *   PBEAML
 *   PBCOMP

All beams are LineProperty objects.
Multi-segment beams are IntegratedLineProperty objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count
from six.moves import zip, range
from numpy import array, unique, argsort, mean, allclose

from pyNastran.bdf.cards.properties.bars import IntegratedLineProperty, LineProperty, _bar_areaL
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    string, string_or_blank, double_string_or_blank)
from pyNastran.utils.mathematics import integrate_line, integrate_positive_line
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class PBEAM(IntegratedLineProperty):
    type = 'PBEAM'

    #_opt_map = {
        #'I1(B)' : 'i1',
        #'I1(B)',
    #}

    def __init__(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                 c1, c2, d1, d2, e1, e2, f1, f2,
                 k1, k2, s1, s2, nsia, nsib, cwa, cwb,
                 m1a, m2a, m1b, m2b,
                 n1a, n2a, n1b, n2b,
                 comment=''):
        """
        .. todo:: fix 0th entry of self.so, self.xxb
        """
        IntegratedLineProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid

        assert isinstance(xxb, list), xxb
        assert isinstance(so, list), so
        #assert min(so) in [0., 1.], so  # YES, NO
        #assert max(so) == 1.0, so
        #print('xxb', xxb)
        assert 0. <= min(xxb) <= 0.0, xxb  # x/L
        assert 0. <= max(xxb) <= 1.0, xxb
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

        # sort xxb
        ixxb = argsort(xxb)
        self.so = array(so, dtype='|U8')[ixxb]
        self.xxb = array(xxb, dtype='float64')[ixxb]
        #print('ixxb = %s' % ixxb)
        #print('i12 = %s' % i12)

        self.A = array(area, dtype='float64')[ixxb]
        self.i1 = array(i1, dtype='float64')[ixxb]
        self.i2 = array(i2, dtype='float64')[ixxb]
        self.i12 = array(i12, dtype='float64')[ixxb]
        self.j = array(j, dtype='float64')[ixxb]
        self.nsm = array(nsm, dtype='float64')[ixxb]

        self.c1 = array(c1, dtype='float64')[ixxb]
        self.c2 = array(c2, dtype='float64')[ixxb]
        self.d1 = array(d1, dtype='float64')[ixxb]
        self.d2 = array(d2, dtype='float64')[ixxb]
        self.e1 = array(e1, dtype='float64')[ixxb]
        self.e2 = array(e2, dtype='float64')[ixxb]
        self.f1 = array(f1, dtype='float64')[ixxb]
        self.f2 = array(f2, dtype='float64')[ixxb]

        # now we interpolate to fix up missing data
        # from arrays that were potentially out of order
        # (they're sorted now)
        #
        # we've also already checked xxb=0.0 and xxb=1.0 for I1, I2, I12, J
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
                    self.i12[i] = self.i2[-1] + self.i2[0] * (1 - xxb)
                if j == 0.0:
                    self.j[i] = self.j[-1] + self.j[0] * (1 - xxb)

                assert self.A[i] >= 0., self.A
                assert self.i1[i] >= 0., self.i1
                assert self.i2[i] >= 0., self.i2
                assert self.j[i] >= 0., self.j  # we check warping later
                di12 = self.i1[i] * self.i2[i] - self.i12[i] ** 2
                if not di12 > 0.:
                    msg = 'I1 * I2 - I12^2=0 and must be greater than 0.0 at End B\n'
                    msg += 'xxb=%s i1=%s i2=%s i12=%s i1*i2-i12^2=%s'  % (
                        self.xxb[i], self.i1[i], self.i2[i], self.i12[i], di12)
                    raise ValueError(msg)

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
            for i, xxb, j in zip(count(), self.xxb, self.j):
                if self.j[i] < 0.:
                    ji = self.j[i]
                    msg = 'Warping Check Error; j[%i] must be greater than 0.0' % i
                    msg += '  cwa=%s cwb=%s\n' % (self.cwa, self.cwb)
                    msg += '  i=%s xxb=%s j=%s; j[%i]=%s\n' % (i, xxb, self.j, i, j)
                    raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
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

        assert area[0] >= 0., area
        assert i1[0] >= 0., i1
        assert i2[0] >= 0., i2

        # we'll do a check for warping later; cwa/cwb -> j > 0.0
        assert j[0] >= 0., j

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
        m2a = double_or_blank(card, ifield + 9, 'm2a', m1a)

        #: y coordinate of center of gravity of
        #: nonstructural mass for end B.
        m1b = double_or_blank(card, ifield + 10, 'm1b', 0.0)
        #: z coordinate of center of gravity of
        #: nonstructural mass for end B.
        m2b = double_or_blank(card, ifield + 11, 'm2b', m1b)

        #: y coordinate of neutral axis for end A.
        n1a = double_or_blank(card, ifield + 12, 'n1a', 0.0)
        #: z coordinate of neutral axis for end A.
        n2a = double_or_blank(card, ifield + 13, 'n2a', n1a)

        #: y coordinate of neutral axis for end B.
        n1b = double_or_blank(card, ifield + 14, 'n1a', 0.0)
        #: z coordinate of neutral axis for end B.
        n2b = double_or_blank(card, ifield + 15, 'n2b', n1b)


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
        (pid, mid, nsegs, ccf, x) = data[:5]

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

        (k1, k2, s1, s2, nsia, nsib, cwa, cwb,
         m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b) = data[-1]
        return PBEAM(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                     c1, c2, d1, d2, e1, e2,
                     f1, f2, k1, k2, s1, s2,
                     nsia, nsib, cwa, cwb, m1a, m2a, m1b,
                     m2b, n1a, n2a, n1b, n2b, comment=comment)

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
            return self.i1[0]
        elif name_str == 'I1(B)':
            return self.i1[-1]

        elif name_str == 'I2(A)':
            return self.i2[-1]
        elif name_str == 'I2(B)':
            return self.i2[-1]
        else:
            raise NotImplementedError(name_str)

    #def Area(self):
    #    """.. warning:: area field not supported fully on PBEAM card"""
    #    #raise RuntimeError(self.A[0])
    #    return self.A[0]

    #def Nsm(self):
    #    """.. warning:: nsm field not supported fully on PBEAM card"""
    #    #raise RuntimeError(self.nsm[0])
    #    return self.nsm[0]

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
            mass_per_lengths.append(area * rho + nsm)
        mass_per_length = integrate_positive_line(self.xxb, mass_per_lengths)
        return mass_per_length

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PBEAM mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid
        #if model.sol != 600:
            #assert max(self.j) == 0.0, self.j
            #assert min(self.j) == 0.0, self.j

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        A = self.Area()
        try:
            J = self.J()
        except NotImplementedError:
            msg = "J is not implemented for pid.type=%s pid.Type=%s" % (self.type, self.Type)
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

    def _write_code_aster(self):  # PBEAM
        a = self.Area()
        iy = self.I11()
        iz = self.I22()
        j = self.J()
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='P%s', # PBEAM\n" % (self.pid)
        msg += "              SECTION='GENERALE',\n"
        msg += "              CARA=('A','IY','IZ','JX'), # area, moments of inertia\n"
        msg += "              VALE=(%g,  %g,  %g,  %g),\n" % (a, iy, iz, j)

        msg += "              ORIENTATION=_F( \n"
        ## .. todo:: is this correct
        msg += "                  CARA=('VECT_Y'), # direction of beam ???\n"
        msg += "                  VALE=(1.0,0.0,0.0,)"

        if [self.n1a, self.n1b] != [0., 0.]:
            msg += "              \n),\n"
            msg += "              CARA=('AX','AY'), # shear centers\n"
            msg += "              VALE=(%g, %g),\n" % (self.n1a, self.n1b)
            msg += "             ),\n"
        else:
            msg += " )),\n"
        return msg

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
        #k1 = set_blank_if_default(self.k1, 1.0)
        #k2 = set_blank_if_default(self.k2, 1.0)
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
        m1a = set_blank_if_default(self.m1a, 0.0)
        m2a = set_blank_if_default(self.m2a, self.m1a)
        m1b = set_blank_if_default(self.m1b, 0.0)
        m2b = set_blank_if_default(self.m2b, self.m1b)

        n1a = set_blank_if_default(self.n1a, 0.0)
        n2a = set_blank_if_default(self.n2a, self.n1a)
        n1b = set_blank_if_default(self.n1b, 0.0)
        n2b = set_blank_if_default(self.n2b, self.n1b)

        #footer = [k1, k2, s1, s2, nsia, nsib, cwa, cwb,
                  #m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b]
        footer = [self.k1, self.k2, s1, s2, nsia, nsib, cwa, cwb,
                  m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b]

        list_fields += footer
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBEAML(IntegratedLineProperty):
    """
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    | PBEAML | PID     | MID     | GROUP   | TYPE    |         |         |         |         |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | DIM1(A) | DIM2(A) | -etc.-  | DIMn(A) | NSM(A)  | SO(1)   | X(1)/XB | DIM1(1) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | DIM2(1) | -etc.-  | DIMn(1) | NSM(1)  | SO(2)   | X(2)/XB | DIM1(2) | DIM2(2) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | -etc.-  | DIMn(2) | NSM(m)  | -etc.-  | SO(m)   | X(m)/XB | DIM1(m) | -etc.-  |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | DIMn(m) |  NSM(m) | SO(B)   | 1.0     | DIM1(B) | DIM2(B) | -etc.-  | DIMn(B) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | NSM(B)  |
    +--------+---------+
    """
    type = 'PBEAML'
    valid_types = {
        "ROD": 1,
        "TUBE": 2,
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

    def __init__(self, pid, mid, group, Type, xxb, so, dims, nsm, comment=''):
        IntegratedLineProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.group = group
        #: Section Type (e.g. 'ROD', 'TUBE', 'I', 'H')
        self.Type = Type
        ndim = self.valid_types[self.Type]

        self.dim = dims
        for xxbi, dim in zip(xxb, dims):
            assert len(dim) == ndim, 'Type=%s ndim=%s len(dim)=%s xxb=%s dim=%s' % (
                Type, ndim, len(dim), xxbi, dim)
        self.xxb = xxb
        self.so = so
        self.nsm = nsm

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        group = string_or_blank(card, 3, 'group', 'MSCBMLO')
        Type = string(card, 4, 'Type')

        # determine the number of required dimensions on the PBEAM
        ndim = cls.valid_types[Type]

        #: dimension list
        dims = []
        dim = []

        #: Section position
        xxb = [0.]

        #: Output flag
        so = ['YES']

        #: non-structural mass :math:`nsm`
        nsm = []

        i = 9
        n = 0
        while i < len(card):
            if n > 0:
                soi = string_or_blank(card, i, 'so_n=%i' % n, 'YES')
                xxbi = double_or_blank(card, i + 1, 'xxb_n=%i' % n, 1.0)
                so.append(soi)
                xxb.append(xxbi)
                i += 2

            dim = []
            for ii in range(ndim):
                dimi = double(card, i, 'dim_n=%i_ii=%i' % (n, ii))
                dim.append(dimi)
                i += 1
            dims.append(dim)

            nsmi = double_or_blank(card, i, 'nsm_n=%i' % n, 0.0)
            nsm.append(nsmi)
            n += 1
            i += 1
        return PBEAML(pid, mid, group, Type, xxb, so, dims, nsm, comment=comment)

    #def add_op2_data(self, data, comment=''):
        #if comment:
            #self._comment = comment
        #raise NotImplementedError(data)

    def _verify(self, xref=False):
        pid = self.Pid()
        rho = self.Rho()
        nsm = self.Nsm()
        area = self.Area()
        mass_per_length = self.MassPerLength()
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))
        assert isinstance(rho, float), 'rho=%r\n%s' % (rho, str(self))
        assert isinstance(nsm, float), 'nsm=%r\n%s' % (nsm, str(self))
        assert isinstance(area, float), 'area=%r\n%s' % (area, str(self))
        assert isinstance(mass_per_length, float), 'mass/L=%r\n%s' % (mass_per_length, str(self))

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the PBEAML.

        .. math:: \frac{m}{L} = A(x) \rho + nsm

        .. math:: \frac{m}{L} = nsm L + \rho \int \, A(x) dx
        """
        rho = self.Rho()
        massPerLs = []
        for (dim, nsm) in zip(self.dim, self.nsm):
            a = _bar_areaL('PBEAML', self.Type, dim)
            try:
                massPerLs.append(a * rho + nsm)
            except:
                msg = "PBEAML a*rho+nsm a=%s rho=%s nsm=%s" % (a, rho, nsm)
                raise RuntimeError(msg)
        massPerL = integrate_positive_line(self.xxb, massPerLs)
        return massPerL

    def Area(self):
        r"""
        Gets the Area :math:`A` of the PBEAML.

        .. math:: A = \int \, A(x) dx

        .. note:: a spline is fit to :math:`A(x)` and then integrated.
        """
        areas = []
        for dim in self.dim:
            areas.append(_bar_areaL('PBEAML', self.Type, dim))
        try:
            A = integrate_line(self.xxb, areas)
        except ValueError:
            print('PBEAML integration error; pid=%s x/xb=%s areas=%s' % (self.pid, self.xxb, areas))
            A = mean(areas)
        return A

    #def Mid(self):
    #    return self.mid

    def cross_reference(self, model):
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
        msg = ' which is required by PBEAML mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

    def verify(self, model, isubcase):
        if model.is_thermal_solution(isubcase):
            assert self.mid_ref.type in ['MAT4', 'MAT5']
        else:
            assert self.mid_ref.type in ['MAT1']

    def _J(self):
        j = []
        for dims in self.dim:
            pass
            #print("dims = ",dims)
            #IAreaL()
        return j

    def J(self):
        #raise NotImplementedError()
        #Js = self._J()
        #j = integrate_positive_line(self.xxb, Js)
        j = None
        return j

    def I11(self):
        #i1 = integrate_positive_line(self.xxb,self.i1)
        i1 = None
        return i1

    def I22(self):
        #i2 = integrate_positive_line(self.xxb,self.i2)
        i2 = None
        return i2

    def I12(self):
        #i12 = integrate_line(self.xxb,self.i12)
        i12 = None
        return i12

    def _write_code_aster(self, icut=0, iface=0, istart=0):  # PBEAML
        msg = ''
        msg2 = 'Cut_%s = geompy.MakeCut(' % (icut + 1)
        for xxb, dim, nsm in zip(self.xxb, self.dim, self.nsm):
            msg += self.CA_Section(iface, istart, self.dim)
            msg2 += 'Face_%i, ' % (iface + 1)
            iface += 1
            istart += len(self.dim)
        msg2 = msg2[-2:]
        msg2 += ')\n'

        msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (
            icut + 1, icut + 1)
        icut += 1
        return (msg + msg2, icut, iface, istart)

    def raw_fields(self):
        list_fields = ['PBEAML', self.pid, self.Mid(), self.group, self.Type,
                       None, None, None, None]
        #print("self.nsm = ",self.nsm)
        #print("xxb=%s so=%s dim=%s nsm=%s" %(self.xxb,self.so,
        #                                     self.dim,self.nsm))
        for (i, xxb, so, dim, nsm) in zip(count(), self.xxb, self.so,
                                          self.dim, self.nsm):
            if i == 0:
                list_fields += dim + [nsm]
            else:
                list_fields += [so, xxb] + dim + [nsm]
        #raise NotImplementedError('verify PBEAML...')
        return list_fields

    def repr_fields(self):
        #group = set_blank_if_default(self.group, 'MSCBMLO')
        list_fields = self.raw_fields()
        #list_fields[3] = group
        return list_fields

    def write_card(self, size=8, is_double=False):
        """.. todo:: having bug with PBEAML"""
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)  #is this allowed???


class PBCOMP(LineProperty):
    type = 'PBCOMP'

    def __init__(self, pid, mid, area, i1, i2, i12, j, nsm,
                 k1, k2, m1, m2, n1, n2, symopt, y, z, c,
                 mids, comment=''):
        LineProperty.__init__(self)
        if comment:
            self._comment = comment
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

    @classmethod
    def add_card(cls, card, comment=''):
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
        return PBCOMP(pid, mid, area, i1, i2, i12, j, nsm,
                      k1, k2, m1, m2, n1, n2,
                      symopt, y, z, c, mids,
                      comment=comment)

    def add_op2_data(self, data, comment=''):
        raise NotImplementedError()

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int)

    def MassPerLength(self):
        return self.nsm + self.mid_ref.Rho() * self.A

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PBCOMP mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

    def raw_fields(self):
        list_fields = ['PBCOMP', self.pid, self.Mid(), self.A, self.i1,
                       self.i2, self.i12, self.j, self.nsm, self.k1, self.k2,
                       self.m1, self.m2, self.n1, self.n2, self.symopt, None]
        for (yi, zi, ci, mid) in zip(self.y, self.z, self.c, self.mids):
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
