# pylint: disable=C0103,R0902,R0904,R0914,C0111
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
#import sys
from six.moves import zip, range
from itertools import count
from numpy import pi, array

from pyNastran.bdf.cards.properties.bars import IntegratedLineProperty, LineProperty, _bar_areaL
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank,
    double_string_or_blank, integer_double_string_or_blank)
from pyNastran.utils.mathematics import integrate_line, integrate_positive_line
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class PBEAM(IntegratedLineProperty):
    type = 'PBEAM'

    def __init__(self, card=None, data=None, comment=''):
        """
        .. todo:: fix 0th entry of self.so, self.xxb
        """
        IntegratedLineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            #: Property ID
            self.pid = integer(card, 1, 'property_id')
            #: Material ID
            self.mid = integer(card, 2, 'material_id')

            # at least one cross section are required
            # so[0] and xxb[0] aren't used
            #: Output flag
            self.so = ['YES']
            #: Section position
            self.xxb = [0.]
            A = double(card, 3, 'Area')
            #: Area
            self.A = [A]
            #: Moment of Inertia about the 1 axis :math:`I_1`
            self.i1 = [double_or_blank(card, 4, 'I1', 0.0)]
            #: Moment of Inertia about the 2 axis  :math:`I_2`
            self.i2 = [double_or_blank(card, 5, 'I2', 0.0)]
            #: Moment of Inertia about the 12 axis  :math:`I_{12}`
            self.i12 = [double_or_blank(card, 6, 'I12', 0.0)]
            #: Polar Moment of Inertia :math:`J`
            self.j = [double_or_blank(card, 7, 'J', 0.0)]
            #: Non-structural mass :math:`nsm`
            self.nsm = [double_or_blank(card, 8, 'nsm', 0.0)]

            assert self.A[0] >= 0., self.A
            assert self.i1[0] >= 0., self.i1
            assert self.i2[0] >= 0., self.i2


            is_cdef = False
            field9 = double_string_or_blank(card, 9, 'field9')
            field17 = double_string_or_blank(card, 17, 'field17')
            try:
                isinstance(field17, float)
                is_cdef = True
                is_footer = True
            except SyntaxError:
                pass
            #print("f9=%s f17=%s" % (field9, field17))

            #nlines = nfields // 8

            if field9 in ['YES', 'YESA', 'NO']:
                is_cdef = False
                is_continue = True
            elif field17 in ['YES', 'YESA', 'NO']:
                is_cdef = True
                is_continue = True
            else:
                is_continue = False
                is_cdef = True
                #if nlines == 2:
                #    isCDEF = True
                #elif nlines == 3:
                #    isCDEF = Tru
                #else:


            #print("isCDEF=%s isContinue=%s" % (isCDEF, isContinue))
            #if isCDEF:
            self.c1 = [double_or_blank(card, 9, 'c1', 0.0)]
            self.c2 = [double_or_blank(card, 10, 'c2', 0.0)]
            self.d1 = [double_or_blank(card, 11, 'd1', 0.0)]
            self.d2 = [double_or_blank(card, 12, 'd2', 0.0)]
            self.e1 = [double_or_blank(card, 13, 'e1', 0.0)]
            self.e2 = [double_or_blank(card, 14, 'e2', 0.0)]
            self.f1 = [double_or_blank(card, 15, 'f1', 0.0)]
            self.f2 = [double_or_blank(card, 16, 'f2', 0.0)]
            #else:
                #msg = ('On PBEAM %s, only YES format is supported.\n'
                #       'All C, D, E, and F fields must be specified.'% self.pid)
                #raise RuntimeError(msg)
                #self.c1 = [None]
                #self.c2 = [None]
                #self.d1 = [None]
                #self.d2 = [None]
                #self.e1 = [None]
                #self.e2 = [None]
                #self.f1 = [None]
                #self.f2 = [None]

            # ----------------------------------------------------------------
            # figure out how many YES/YESA/NO fields there are
            # and if there is a footer

            # counting continuation cards
            nfields = len(card) - 1  # -1 for PBEAM field
            if is_cdef:
                nmajor = nfields // 16
                nleftover = nfields % 16
                if nleftover == 0:
                    nmajor -= 1

                if nmajor == 0:
                    nmajor = 1

                # jump to the last block of 16
                x = nmajor * 16 + 1

                # If it's an SO field, we don't read the footer
                # remark 6:
                # The fourth and fifth continuation entries, which
                # contain fields K1 through N2(B), are optional
                # and may be omitted if the default values are appropriate.
                val = integer_double_string_or_blank(card, x, 'YES/YESA/NO')
                if val in ['YES', 'YESA', 'NO']:  # there is no footer
                    nmajor += 1
                    x += 16
                else:
                    # read the footer
                    pass
            else:
                nmajor = nfields // 8
                nleftover = nfields % 8
                if nleftover == 0:
                    nmajor -= 1
                if nmajor == 0:
                    nmajor = 1
                x = nmajor * 8 + 1

                val = integer_double_string_or_blank(card, x, 'YES/YESA/NO')
                if val in ['YES', 'YESA', 'NO']:  # there is no footer
                    nmajor += 1
                    x += 8
                else:
                    # read the footer
                    pass

            # ----------------------------------------------------------------
            for nrepeated in range(1, nmajor): # start at 1 to drop the header
                # field 17 is the first possible so
                if is_cdef:
                    nstart = nrepeated * 16 + 1
                else:
                    nstart = nrepeated * 8 + 1

                so = string(card, nstart, 'SO%i' % nrepeated)
                xxb = double(card, nstart + 1, 'x/xb%i' % nrepeated)
                A = double_or_blank(card, nstart + 2, 'Area%i' % nrepeated, 0.0)
                i1 = double_or_blank(card, nstart + 3, 'I1 %i' % nrepeated, 0.0)
                i2 = double_or_blank(card, nstart + 4, 'I2 %i' % nrepeated, 0.0)
                i12 = double_or_blank(card, nstart + 5, 'I12 %i' % nrepeated, 0.0)
                j = double_or_blank(card, nstart + 6, 'J%i' % nrepeated, 0.0)
                nsm = double_or_blank(card, nstart + 7, 'nsm%i' % nrepeated, 0.0)

                self.so.append(so)
                self.xxb.append(xxb)
                self.A.append(A)
                self.i1.append(i1)
                self.i2.append(i2)
                self.i12.append(i12)
                self.j.append(j)
                self.nsm.append(nsm)

                if is_cdef:
                    c1 = double_or_blank(card, nstart + 8, 'c1 %i' % nrepeated, 0.0)
                    c2 = double_or_blank(card, nstart + 9, 'c2 %i' % nrepeated, 0.0)
                    d1 = double_or_blank(card, nstart + 10, 'd1 %i' % nrepeated, 0.0)
                    d2 = double_or_blank(card, nstart + 11, 'd2 %i' % nrepeated, 0.0)
                    e1 = double_or_blank(card, nstart + 12, 'e1 %i' % nrepeated, 0.0)
                    e2 = double_or_blank(card, nstart + 13, 'e2 %i' % nrepeated, 0.0)
                    f1 = double_or_blank(card, nstart + 14, 'f1 %i' % nrepeated, 0.0)
                    f2 = double_or_blank(card, nstart + 15, 'f2 %i' % nrepeated, 0.0)
                    self.c1.append(c1)
                    self.c2.append(c2)
                    self.d1.append(d1)
                    self.d2.append(d2)
                    self.e1.append(e1)
                    self.e2.append(e2)
                    self.f1.append(f1)
                    self.f2.append(f2)
                else:
                    # YESA or NO, values MUST be omitted; remark 5
                    self.c1.append(None)
                    self.c2.append(None)
                    self.d1.append(None)
                    self.d2.append(None)
                    self.e1.append(None)
                    self.e2.append(None)
                    self.f1.append(None)
                    self.f2.append(None)
            if len(self.xxb) > 1:
                assert min(self.xxb) == 0.0, 'min=%s, but should be 0.0\nxxb=%s' % (min(self.xxb), self.xxb)
                assert max(self.xxb) == 1.0, 'max=%s, but should be 1.0\nxxb=%s' % (max(self.xxb), self.xxb)


            # footer fields
            #: Shear stiffness factor K in K*A*G for plane 1.
            self.k1 = double_or_blank(card, x, 'k1', 1.0)
            #: Shear stiffness factor K in K*A*G for plane 2.
            self.k2 = double_or_blank(card, x + 1, 'k2', 1.0)

            #: Shear relief coefficient due to taper for plane 1.
            self.s1 = double_or_blank(card, x + 2, 's1', 0.0)
            #: Shear relief coefficient due to taper for plane 2.
            self.s2 = double_or_blank(card, x + 3, 's2', 0.0)

            #: non structural mass moment of inertia per unit length
            #: about nsm center of gravity at Point A.
            self.nsia = double_or_blank(card, x + 4, 'nsia', 0.0)
            #: non structural mass moment of inertia per unit length
            #: about nsm center of gravity at Point B.
            self.nsib = double_or_blank(card, x + 5, 'nsib', self.nsia)

            #: warping coefficient for end A.
            self.cwa = double_or_blank(card, x + 6, 'cwa', 0.0)
            #: warping coefficient for end B.
            self.cwb = double_or_blank(card, x + 7, 'cwb', self.cwa)

            #: y coordinate of center of gravity of
            #: nonstructural mass for end A.
            self.m1a = double_or_blank(card, x + 8, 'm1a', 0.0)
            #: z coordinate of center of gravity of
            #: nonstructural mass for end A.
            self.m2a = double_or_blank(card, x + 9, 'm2a', self.m1a)

            #: y coordinate of center of gravity of
            #: nonstructural mass for end B.
            self.m1b = double_or_blank(card, x + 10, 'm1b', 0.0)
            #: z coordinate of center of gravity of
            #: nonstructural mass for end B.
            self.m2b = double_or_blank(card, x + 11, 'm2b', self.m1b)

            #: y coordinate of neutral axis for end A.
            self.n1a = double_or_blank(card, x + 12, 'n1a', 0.0)
            #: z coordinate of neutral axis for end A.
            self.n2a = double_or_blank(card, x + 13, 'n2a', self.n1a)

            #: y coordinate of neutral axis for end B.
            self.n1b = double_or_blank(card, x + 14, 'n1a', 0.0)
            #: z coordinate of neutral axis for end B.
            self.n2b = double_or_blank(card, x + 15, 'n2b', self.n1b)
        else:
            raise NotImplementedError(data)

    #def Area(self):
    #    """.. warning:: area field not supported fully on PBEAM card"""
    #    #raise RuntimeError(self.A[0])
    #    return self.A[0]

    #def Nsm(self):
    #    """.. warning:: nsm field not supported fully on PBEAM card"""
    #    #raise RuntimeError(self.nsm[0])
    #    return self.nsm[0]

    def I1_I2_I12(self):
        assert self.i1  is not None, 'I1=%r' % self.i1
        assert self.i2  is not None, 'I2=%r' % self.i2
        assert self.i12 is not None, 'I12=%r' % self.i12
        return self.i1[0], self.i2[0], self.i12[0]

    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho = self.Rho()
        massPerLs = []
        for (area, nsm) in zip(self.A, self.nsm):
            massPerLs.append(area * rho + nsm)
        massPerL = integrate_positive_line(self.xxb, massPerLs)
        return massPerL

    def cross_reference(self, model):
        msg = ' which is required by PBEAM mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        #if model.sol != 600:
            #assert max(self.j) == 0.0, self.j
            #assert min(self.j) == 0.0, self.j

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
            assert self.mid.type in ['MAT1', 'MAT4', 'MAT5'], 'pid.type=%s; mid.type=%s' % (self.type, self.mid.type)
            #self.MassPerLength()

    def writeCodeAster(self):  # PBEAM
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
        msg += "                  CARA=('VECT_Y'), # direction of beam ???\n"  ## .. todo:: is this correct
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
        for (so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2,
             f1, f2) in zip(self.so, self.xxb, self.A, self.i1, self.i2,
                            self.i12, self.j, self.nsm, self.c1, self.c2,
                            self.d1, self.d2, self.e1, self.e2, self.f1,
                            self.f2):

            if i == 0:  # the first 2 fields aren't written
                list_fields += [A, i1, i2, i12, j, nsm,
                                c1, c2, d1, d2, e1, e2, f1, f2]
            else:
                list_fields += [so, xxb, A, i1, i2, i12, j, nsm,
                                c1, c2, d1, d2, e1, e2, f1, f2]
            i += 1

        footer = [self.k1, self.k2, self.s1, self.s2, self.nsia, self.nsib,
                  self.cwa, self.cwb, self.m1a, self.m2a, self.m1b, self.m2b,
                  self.n1a, self.n2a, self.n1b, self.n2b]
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
                list_fields += [so, xxb, A, i1, i2, i12, j, nsm,
                                c1, c2, d1, d2, e1, e2, f1, f2]
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
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


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
    validTypes = {
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
        "DBOX": 10,  # was 12
    }  # for GROUP="MSCBML0"

    def __init__(self, card=None, data=None, comment=''):
        IntegratedLineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            self.group = string_or_blank(card, 3, 'group', 'MSCBMLO')
            #: Section Type (e.g. 'ROD', 'TUBE', 'I', 'H')
            self.Type = string(card, 4, 'Type')

            # determine the number of required dimensions on the PBEAM
            ndim = self.validTypes[self.Type]
            #nAll = ndim + 1

            #: dimension list
            self.dim = []
            Dim = []

            #: Section position
            self.xxb = [0.]

            #: Output flag
            self.so = ['YES']

            #: non-structural mass :math:`nsm`
            self.nsm = []

            i = 9
            n = 0
            while i < len(card):
                if n > 0:
                    so = string_or_blank(card, i, 'so_n=%i' % n, 'YES')
                    xxb = double_or_blank(card, i + 1, 'xxb_n=%i' % n, 1.0)
                    self.so.append(so)
                    self.xxb.append(xxb)
                    i += 2

                Dim = []
                for ii in range(ndim):
                    dim = double(card, i, 'dim_n=%i_ii=%i' % (n, ii))
                    Dim.append(dim)
                    i += 1
                self.dim.append(Dim)

                nsm = double_or_blank(card, i, 'nsm_n=%i' % n, 0.0)
                self.nsm.append(nsm)
                n += 1
                i += 1

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
        Areas = []
        for dim in self.dim:
            Areas.append(_bar_areaL('PBEAML', self.Type, dim))
        A = integrate_line(self.xxb, Areas)
        return A

    #def Mid(self):
    #    return self.mid

    def cross_reference(self, model):
        """
        .. warning:: For structural problems, PBEAML entries must
                     reference a MAT1 material entry
        .. warning:: For heat-transfer problems, the MID must
                     reference a MAT4 or MAT5 material entry.
        .. todo:: What happens when there are 2 subcases?
        """
        msg = ' which is required by PBEAML mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

    def verify(self, model, isubcase):
        if model.is_thermal_solution(isubcase):
            assert self.mid.type in ['MAT4', 'MAT5']
        else:
            assert self.mid.type in ['MAT1']

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

    def writeCodeAster(self, iCut=0, iFace=0, iStart=0):  # PBEAML
        msg = ''
        msg2 = 'Cut_%s = geompy.MakeCut(' % (iCut + 1)
        for (xxb, so, dim, nsm) in zip(self.xxb, self.so, self.dim, self.nsm):
            msg += self.CA_Section(iFace, iStart, self.dim)
            msg2 += 'Face_%i, ' % (iFace + 1)
            iFace += 1
            iStart += len(self.dim)
        msg2 = msg2[-2:]
        msg2 += ')\n'

        msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (
            iCut + 1, iCut + 1)
        iCut += 1
        return (msg + msg2, iCut, iFace, iStart)

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

    def write_card(self, size, is_double):
        """.. todo:: having bug with PBEAML"""
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)  #is this allowed???


class PBCOMP(LineProperty):
    type = 'PBCOMP'

    def __init__(self, card=None, data=None, comment=''):
        LineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            # Area
            self.A = double_or_blank(card, 3, 'A', 0.0)
            self.i1 = double_or_blank(card, 4, 'I1', 0.0)
            self.i2 = double_or_blank(card, 5, 'I2', 0.0)
            self.i12 = double_or_blank(card, 6, 'I12', 0.0)
            #: Polar Moment of Inertia :math:`J`
            self.j = double_or_blank(card, 7, 'J', 0.0)
            #: Non-structural mass per unit length :math:`\frac{m}{L}`
            self.nsm = double_or_blank(card, 8, 'nsm', 0.0)
            self.k1 = double_or_blank(card, 9, 'k1', 1.0)
            self.k2 = double_or_blank(card, 10, 'k2', 1.0)
            self.m1 = double_or_blank(card, 11, 'm1', 0.0)
            self.m2 = double_or_blank(card, 12, 'm2', 0.0)
            self.n1 = double_or_blank(card, 13, 'n1', 0.0)
            self.n2 = double_or_blank(card, 14, 'n2', 0.0)
            self.symopt = integer_or_blank(card, 15, 'symopt', 0)
            assert 0 <= self.symopt <= 5, 'symopt=%i is invalid; ' % self.symopt
            self.y = []
            self.z = []
            self.c = []
            self.mids = []

            nfields = len(card) - 17
            nrows = nfields // 8
            if nfields % 8 > 0:
                nrows += 1

            for row in range(nrows):
                i = 8 * row + 17
                yi = double(card, i, 'y' + str(row))
                zi = double(card, i + 1, 'z' + str(row))
                ci = double_or_blank(card, i + 2, 'c' + str(row), 0.0)
                mid = integer_or_blank(card, i + 3, 'mid' + str(row), self.mid)
                self.y.append(yi)
                self.z.append(zi)
                self.c.append(ci)
                self.mids.append(mid)

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int)

    def MassPerLength(self):
        return self.nsm + self.mid.Rho() * self.A

    def cross_reference(self, model):
        msg = ' which is required by PBCOMP mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

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
            return self.comment() + print_card_8(card)
        #if is_double:
            #return self.comment() + print_card_double(card)
        return self.comment() + print_card_16(card)
