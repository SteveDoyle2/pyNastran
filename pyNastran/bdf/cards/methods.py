# pylint: disable=C0103,R0902,R0904,R0914
"""
All method cards are defined in this file.  This includes:

 * EIGB
 * EIGC
 * EIGR
 * EIGP
 * EIGRL

All cards are Method objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from six.moves import zip, range

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string, string_or_blank,
    components, components_or_blank, integer_double_string_or_blank, blank,
    interpret_value)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class Method(BaseCard):
    """
    Generic class for all methods.
    Part of self.methods
    """
    def __init__(self):
        pass


class EIGB(Method):
    """
    Defines data needed to perform buckling analysis
    """
    type = 'EIGB'

    def __init__(self, sid, method, L1, L2, nep, ndp, ndn, norm, G, C, comment=''):
        Method.__init__(self)
        if comment:
            self._comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid

        #: Method of eigenvalue extraction. (Character: 'INV' for inverse
        #: power method or 'SINV' for enhanced inverse power method.)
        #: apparently it can also be blank...
        self.method = method

        #: Eigenvalue range of interest. (Real, L1 < L2)
        self.L1 = L1
        self.L2 = L2

        #: Estimate of number of roots in positive range not used for
        #: METHOD = 'SINV'. (Integer > 0)
        self.nep = nep

        #: Desired number of positive and negative roots.
        #: (Integer>0; Default = 3*NEP)
        self.ndp = ndp
        self.ndn = ndn

        #: Method for normalizing eigenvectors.
        #: ('MAX' or 'POINT';Default='MAX')
        self.norm = norm
        self.G = G
        self.C = C
        if not self.L1 < self.L2:
            msg = 'L1=%s L2=%s; L1<L2 is requried' % (self.L1, self.L2)
            raise RuntimeError(msg)
        if self.method not in ['INV', 'SINV', None]:
            msg = 'method must be INV or SINV.  method=%r' % self.method
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        method = string_or_blank(card, 2, 'method')

        L1 = double(card, 3, 'L1')
        L2 = double(card, 4, 'L2')

        nep = integer_or_blank(card, 5, 'nep', 0)
        ndp = integer_or_blank(card, 6, 'ndp', 3 * nep)
        ndn = integer_or_blank(card, 7, 'ndn', 3 * nep)

        norm = string_or_blank(card, 9, 'norm', 'MAX')
        if norm == 'POINT':
            G = integer(card, 10, 'G')
            C = components(card, 11, 'C')
        else:
            G = integer_or_blank(card, 10, 'G')
            C = components_or_blank(card, 11, 'C')
        assert len(card) <= 12, 'len(EIGB card) = %i\ncard=%s' % (len(card), card)
        return EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                    comment=comment)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        list_fields = ['EIGB', self.sid, self.method, self.L1, self.L2, self.nep,
                       self.ndp, self.ndn, None, self.norm, self.G, self.C]
        return list_fields

    def repr_fields(self):
        #method = set_blank_if_default(self.method,'INV')
        nep = set_blank_if_default(self.nep, 0)
        ndp = set_blank_if_default(self.ndp, 3 * self.nep)
        ndn = set_blank_if_default(self.ndn, 3 * self.nep)

        norm = set_blank_if_default(self.norm, 'MAX')
        list_fields = ['EIGB', self.sid, self.method, self.L1, self.L2, nep, ndp,
                       ndn, None, norm, self.G, self.C]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGC(Method):
    """
    Defines data needed to perform complex eigenvalue analysis
    .. todo: not done
    """
    type = 'EIGC'

    def __init__(self, sid, method, norm, G, C, E, ndo, # common
                 mblkszs=None, iblkszs=None, ksteps=None, NJIs=None, # CLAN
                 alphaAjs=None, omegaAjs=None, alphaBjs=None, omegaBjs=None, # HESS/INV
                 LJs=None, NEJs=None, NDJs=None, # HESS/INV
                 shift_r1=None, shift_i1=None, isrr_flag=None, nd1=None, # ISRR
                 comment=''):
        Method.__init__(self)
        if comment:
            self._comment = comment

        #: Set identification number. (Unique Integer > 0)
        self.sid = sid

        #: Method of complex eigenvalue extraction
        #:   MSC 2014 = [INV, HESS, CLAN, IRAM]
        #:   NX 8.5 = [INV, HESS, CLAN, ISRR]
        #:   Autodesk 2015 = [ARNO, HESS, CLAN]
        self.method = method

        #: Method for normalizing eigenvectors
        self.norm = norm
        #: Grid or scalar point identification number. Required only if
        #: NORM='POINT'. (Integer>0)
        self.G = G


        #: Component number. Required only if NORM='POINT' and G is a
        #: geometric grid point. (1<Integer<6)
        self.C = C

        #: Convergence criterion. (Real > 0.0. Default values are:
        #: 10^-4 for METHOD = "INV",
        #: 10^-8 for METHOD = "CLAN",
        #: 10^-8 for METHOD = "ISRR",
        #: 10^-15 for METHOD = "HESS",
        #: E is machine dependent for METHOD = "CLAN".)
        self.E = E
        self.ndo = ndo

        # CLAN
        self.mblkszs = mblkszs
        self.iblkszs = iblkszs
        self.ksteps = ksteps
        self.NJIs = NJIs

        # HESS
        self.alphaBjs = alphaBjs
        self.omegaBjs = omegaBjs
        self.LJs = LJs
        self.NEJs = NEJs
        self.NDJs = NDJs

        self.alphaAjs = alphaAjs
        self.omegaAjs = omegaAjs
        #self.alphaBjs = []
        self.omegaBjs = []
        #self.LJs = []
        #self.NEJs = []
        #self.NDJs = []

        #----------
        # ISRR
        self.shift_r1 = shift_r1
        self.shift_i1 = shift_i1
        self.isrr_flag = isrr_flag
        self.nd1 = nd1

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        method = string(card, 2, 'method')
        assert method in ['ARNO', 'INV', 'HESS', 'CLAN', 'ISRR', 'IRAM'], (
            'method=%s is not ARNO, INV, HESS, CLAN, ISRR, IRAM' % method)

        norm = string_or_blank(card, 3, 'norm')
        if norm == 'POINT':
            G = integer(card, 4, 'G')
            C = components(card, 5, 'C')
        else:
            G = blank(card, 4, 'G')
            C = blank(card, 5, 'C')

        E = double_or_blank(card, 6, 'E')
        ndo = integer_double_string_or_blank(card, 7, 'ND0')

        # ALPHAAJ OMEGAAJ ALPHABJ OMEGABJ LJ NEJ NDJ
        fields = [interpret_value(field) for field in card[9:]]

        #-------CLAN--------------
        mblkszs = []
        iblkszs = []
        ksteps = []
        NJIs = []
        #-------CLAN--------------

        #-------HESS--------------
        alphaAjs = []
        alphaBjs = []
        omegaAjs = []
        omegaBjs = []
        mblkszs = []
        iblkszs = []
        ksteps = []
        LJs = []
        NEJs = []
        NDJs = []
        #-------HESS--------------

        #-------ISRR--------------
        shift_r1 = 0.0
        shift_i1 = 0.0
        isrr_flag = 0
        nd1 = None
        #-------ISRR--------------
        nfields = len(fields)
        nrows = nfields // 8
        if nfields % 8 > 0:
            nrows += 1
        #if nrows == 0:
            #msg = 'invalid row count=0; nfields=%s \ncard=%s\nfields=%s' % (
                #nfields, card, fields)
            #raise RuntimeError(msg)

        if method == 'CLAN':
            alphaAjs, omegaAjs, mblkszs, iblkszs, ksteps, NJIs = cls._load_clan(nrows, card)
        elif method in ['HESS', 'INV']:  # HESS, INV
            alphaAjs, omegaAjs, alphaBjs, omegaBjs, LJs, NEJs, NDJs = cls._load_hess_inv(
                nrows, method, card)
        elif method == 'ISRR':
            shift_r1, shift_i1, isrr_flag, nd1 = cls._load_isrr(nrows, card)
        else:
            msg = 'invalid EIGC method...method=%r' % method
            raise RuntimeError(msg)
        #assert card.nFields() < 8, 'card = %s' % card
        return EIGC(sid, method, norm, G, C, E, ndo,
                    mblkszs, iblkszs, ksteps, NJIs, # CLAN
                    alphaAjs, omegaAjs, alphaBjs, omegaBjs, LJs, NEJs, NDJs, # HESS/INV
                    shift_r1, shift_i1, isrr_flag, nd1, # ISRR
                    comment=comment)

    @staticmethod
    def _load_isrr(nrows, card):
        shift_r1 = []
        shift_i1 = []
        isrr_flag = []
        nd1 = []
        for irow in range(nrows):
            i = 9 + 8 * irow
            shift_r1i = double_or_blank(card, i, 'SHIFT_R1', 0.0)
            shift_i1i = double_or_blank(card, i + 1, 'SHIFT_I1', 0.0)
            #2
            #3
            #4
            isrr_flagi = integer_or_blank(card, i + 5, 'ISRR_FLAG', 0)
            nd1i = integer(card, i + 6, 'ND1')
            shift_r1.append(shift_r1i)
            shift_i1.append(shift_i1i)
            isrr_flag.append(isrr_flagi)
            nd1.append(nd1i)
        return shift_r1, shift_i1, isrr_flag, nd1

    @staticmethod
    def _load_clan(nrows, card):
        alphaAjs = []
        omegaAjs = []
        mblkszs = []
        iblkszs = []
        ksteps = []
        NJIs = []
        for irow in range(nrows):
            #NDJ_default = None
            i = 9 + 8 * irow
            alphaAjs.append(
                double_or_blank(card, i, 'alpha' + str(irow), 0.0))
            omegaAjs.append(
                double_or_blank(card, i + 1, 'omega' + str(irow), 0.0))
            mblkszs.append(
                double_or_blank(card, i + 2, 'mblock' + str(irow), 7))

            iblkszs.append(
                integer_or_blank(card, i + 3, 'iblksz' + str(irow), 2))
            ksteps.append(
                integer_or_blank(card, i + 4, 'kstep' + str(irow), 5))
            NJIs.append(
                integer(card, i + 6, 'NJI' + str(irow)))
        return alphaAjs, omegaAjs, mblkszs, iblkszs, ksteps, NJIs

    @staticmethod
    def _load_hess_inv(nrows, method, card):
        alphaOmega_default = None
        LJ_default = None
        if method == 'INV':
            alphaOmega_default = 0.0
            LJ_default = 1.0

        alphaAjs = []
        alphaBjs = []
        omegaAjs = []
        omegaBjs = []
        #mblkszs = []
        #iblkszs = []
        #ksteps = []
        LJs = []
        NEJs = []
        NDJs = []
        for irow in range(nrows):
            NEj = integer_or_blank(card, 9 + 7 * irow + 5, 'NE%s' % str(irow), 0)
            NDJ_default = None
            if method == 'INV':
                NDJ_default = 3 * NEj

            i = 9 + 8 * irow
            alphaAjs.append(
                double_or_blank(card, i, 'alphaA' + str(irow), alphaOmega_default))
            omegaAjs.append(
                double_or_blank(card, i + 1, 'omegaA' + str(irow), alphaOmega_default))
            alphaBjs.append(
                double_or_blank(card, i + 2, 'alphaB' + str(irow), alphaOmega_default))
            omegaBjs.append(
                double_or_blank(card, i + 3, 'omegaB' + str(irow), alphaOmega_default))
            LJs.append(
                double_or_blank(card, i + 4, 'LJ' + str(irow), LJ_default))
            NEJs.append(NEj)
            NDJs.append(
                integer_or_blank(card, i + 6, 'NDJ' + str(irow), NDJ_default))
        return alphaAjs, omegaAjs, alphaBjs, omegaBjs, LJs, NEJs, NDJs

    def cross_reference(self, model):
        pass

    def raw_method(self):
        list_fields = []
        if self.method in ['HESS', 'INV']:
            for (alphaA, omegaA, alphaB, omegaB, Lj, NEj, NDj) in zip(
                    self.alphaAjs, self.omegaAjs, self.alphaBjs, self.omegaBjs,
                    self.LJs, self.NEJs, self.NDJs):
                alphaA = set_blank_if_default(alphaA, 0.0)
                omegaA = set_blank_if_default(omegaA, 0.0)
                alphaB = set_blank_if_default(alphaB, 0.0)
                omegaB = set_blank_if_default(omegaB, 0.0)
                list_fields += [alphaA, omegaA, alphaB, omegaB, Lj, NEj, NDj, None]

        elif self.method == 'CLAN':
            assert len(self.alphaAjs) == len(self.omegaAjs)
            assert len(self.alphaAjs) == len(self.mblkszs)
            assert len(self.alphaAjs) == len(self.iblkszs)
            assert len(self.alphaAjs) == len(self.ksteps)
            assert len(self.alphaAjs) == len(self.NJIs)
            for (alphaA, omegaA, mblksz, iblksz, kstep, Nj) in zip(
                    self.alphaAjs, self.omegaAjs, self.mblkszs, self.iblkszs,
                    self.ksteps, self.NJIs):
                alphaA = set_blank_if_default(alphaA, 0.0)
                omegaA = set_blank_if_default(omegaA, 0.0)
                mblksz = set_blank_if_default(mblksz, 7)
                iblksz = set_blank_if_default(iblksz, 2)
                kstep = set_blank_if_default(kstep, 5)

                list_fields += [alphaA, omegaA, mblksz, iblksz,
                                kstep, None, Nj, None]
        elif self.method == 'ISRR':
            for shift_r1i, shift_i1i, isrr_flagi, nd1i in zip(
                self.shift_r1, self.shift_i1, self.isrr_flag, self.nd1):
                list_fields += [shift_r1i, shift_i1i, None, None, None, isrr_flagi, nd1i, None]

        else:
            msg = 'invalid EIGC method...method=%r' % self.method
            raise RuntimeError(msg)
        return list_fields

    def repr_method(self):
        return self.raw_method()

    def raw_fields(self):
        list_fields = ['EIGC', self.sid, self.method, self.norm, self.G, self.C,
                       self.E, self.ndo, None]
        list_fields += self.raw_method()
        return list_fields

    def repr_fields(self):
        if self.E is None:
            E = None
        else:
            E = str(self.E)
        list_fields = ['EIGC', self.sid, self.method, self.norm, self.G, self.C,
                       E, self.ndo, None]
        list_fields += self.repr_method()
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGP(Method):
    """
    Defines poles that are used in complex eigenvalue extraction by the
    Determinant method.
    """
    type = 'EIGP'

    def __init__(self, sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=''):
        Method.__init__(self)
        if comment:
            self._comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid
        #: Coordinates of point in complex plane. (Real)
        self.alpha1 = alpha1
        #: Coordinates of point in complex plane. (Real)
        self.omega1 = omega1
        #: Multiplicity of complex root at pole defined by point at ALPHAi
        #: and OMEGAi
        self.m1 = m1

        #: Coordinates of point in complex plane. (Real)
        self.alpha2 = alpha2
        #: Coordinates of point in complex plane. (Real)
        self.omega2 = omega2
        #: Multiplicity of complex root at pole defined by point at ALPHAi
        #: and OMEGAi
        self.m2 = m2

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')

        alpha1 = double(card, 2, 'alpha1')
        omega1 = double(card, 3, 'omega1')
        m1 = integer(card, 4, 'm1')

        alpha2 = double(card, 5, 'alpha2')
        omega2 = double(card, 6, 'omega2')
        m2 = integer(card, 7, 'm2')
        assert len(card) == 8, 'len(EIGP card) = %i\ncard=%s' % (len(card), card)
        return EIGP(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=comment)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        list_fields = ['EIGP', self.alpha1, self.omega1, self.m1,
                       self.alpha2, self.omega2, self.m2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGR(Method):
    """
    Defines data needed to perform real eigenvalue analysis
    """
    type = 'EIGR'
    allowed_methods = ['LAN', 'AHOU', 'INV', 'SINV', 'GIV', 'MGIV',
                       'HOU', 'MHOU', 'AGIV']

    def __init__(self, sid, method, f1, f2, ne, nd, norm, G, C, comment=''):
        Method.__init__(self)
        if comment:
            self._comment = comment

        #: Set identification number. (Unique Integer > 0)
        self.sid = sid

        #: Method of eigenvalue extraction. (Character: 'INV' for inverse
        #: power method or 'SINV' for enhanced inverse power method.)
        self.method = method

        #: Frequency range of interest
        self.f1 = f1
        self.f2 = f2

        #: Estimate of number of roots in range (Required for
        #: METHOD = 'INV'). Not used by 'SINV' method.
        self.ne = ne

        #: Desired number of roots (default=600 for SINV 3*ne for INV)
        self.nd = nd

        #: Method for normalizing eigenvectors. ('MAX' or 'POINT';
        #: Default='MAX')
        self.norm = norm

        #: Grid or scalar point identification number. Required only if
        #: NORM='POINT'. (Integer>0)
        self.G = G

        #: Component number. Required only if NORM='POINT' and G is a
        #: geometric grid point. (1<Integer<6)
        self.C = C

        if self.method not in self.allowed_methods:
            msg = 'method=%s; allowed_methods=[%s]' % (
                self.method, ', '.join(self.allowed_methods))
            raise ValueError(msg)
        assert norm in ['POINT', 'MASS', 'MAX']

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        method = string_or_blank(card, 2, 'method', 'LAN')

        f1 = double_or_blank(card, 3, 'f1')
        f2 = double_or_blank(card, 4, 'f2')
        ne = integer_or_blank(card, 5, 'ne')

        if method not in cls.allowed_methods:
            msg = 'method=%s; allowed_methods=[%s]' % (
                method, ', '.join(cls.allowed_methods))
            raise ValueError(msg)

        if method == 'SINV':
            nd = integer_or_blank(card, 6, 'nd', 600)
        elif method == 'INV':
            nd = integer_or_blank(card, 6, 'nd', 3 * ne)
        elif method in ['GIV', 'MGIV', 'HOU', 'MHOU']:
            nd = integer_or_blank(card, 6, 'nd', 0)
        else:
            nd = integer(card, 6, 'nd')
        norm = string_or_blank(card, 9, 'norm', 'MASS')

        if method == 'POINT':
            G = integer(card, 10, 'G')
            C = components(card, 11, 'C')
        else:
            G = blank(card, 10, 'G')
            C = blank(card, 11, 'C')
        assert len(card) <= 12, 'len(EIGR card) = %i\ncard=%s' % (len(card), card)
        return EIGR(sid, method, f1, f2, ne, nd, norm, G, C, comment=comment)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        list_fields = ['EIGR', self.sid, self.method, self.f1, self.f2, self.ne,
                       self.nd, None, None, self.norm, self.G, self.C]
        return list_fields

    def repr_fields(self):
        method = set_blank_if_default(self.method, 'LAN')
        norm = set_blank_if_default(self.norm, 'MASS')
        list_fields = ['EIGR', self.sid, method, self.f1, self.f2, self.ne,
                       self.nd, None, None, norm, self.G, self.C]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGRL(Method):
    """
    Defines data needed to perform real eigenvalue (vibration or buckling)
    analysis with the Lanczos method
    """
    type = 'EIGRL'

    def __init__(self, sid, v1, v2, nd, msglvl, maxset, shfscl, norm,
                 options, values, comment=''):
        Method.__init__(self)
        if comment:
            self._comment = comment

        #: Set identification number. (Unique Integer > 0)
        self.sid = sid

        #: For vibration analysis: frequency range of interest. For
        #: buckling analysis: eigenvalue range of interest. See Remark 4.
        #: (Real or blank, -5 10e16 <= V1 < V2 <= 5.10e16)
        self.v1 = v1
        self.v2 = v2

        #: Number of roots desired
        self.nd = nd

        #: Diagnostic level. (0 < Integer < 4; Default = 0)
        self.msglvl = msglvl

        #: Number of vectors in block or set. Default is machine dependent
        self.maxset = maxset

        #: Estimate of the first flexible mode natural frequency
        #: (Real or blank)
        self.shfscl = shfscl

        #: Method for normalizing eigenvectors (Character: 'MASS' or 'MAX')
        self.norm = norm
        self.options = options
        self.values = values

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        v1 = double_or_blank(card, 2, 'v1')
        v2 = double_or_blank(card, 3, 'v2')
        nd = integer_or_blank(card, 4, 'nd')
        msglvl = integer_or_blank(card, 5, 'msglvl', 0)
        maxset = integer_or_blank(card, 6, 'maxset')
        shfscl = double_or_blank(card, 7, 'shfscl')
        norm = string_or_blank(card, 8, 'norm')

        option_values = [interpret_value(field) for field in card[9:]]
        options = []
        values = []
        for option_value in option_values:
            try:
                (option, value) = option_value.split('=')
            except AttributeError:
                msg = 'parsing EIGRL card incorrectly; option_values=%s\ncard=%s' % (
                    option_values, card)
                raise RuntimeError(msg)
            options.append(option)
            values.append(value)

        #: Method for normalizing eigenvectors
        #if sol in [103, 115, 146]:
            ## normal modes,cyclic normal modes, flutter
            #self.norm = string_or_blank(card, 8, 'norm', 'MASS')
        #elif sol in [105, 110, 111, 116]:
            ## buckling, modal complex eigenvalues,
            ## modal frequency response,cyclic buckling
            #self.norm = string_or_blank(card, 8, 'norm', 'MAX')
        #else:
        norm = string_or_blank(card, 8, 'norm')

        #assert len(card) <= 9, 'len(EIGRL card) = %i\ncard=%s' % (len(card), card)
        assert len(card) <= 10, 'len(EIGRL card) = %i\ncard=%s' % (len(card), card)

        #msg = 'norm=%s sol=%s' % (self.norm, sol)
        #assert self.norm in ['MASS', 'MAX'],msg
        #assert card.nFields()<9,'card = %s' %(card.fields(0))
        return EIGRL(sid, v1, v2, nd, msglvl, maxset, shfscl, norm,
                     options, values, comment=comment)


    def cross_reference(self, model):
        pass
        #if self.norm is None:
            #if model.is_modal_solution():
                #self.norm = 'MASS'
            #elif mdoel.is_buckling_solution():
                #self.norm = 'MAX'

    def raw_fields(self):
        list_fields = ['EIGRL', self.sid, self.v1, self.v2, self.nd,
                       self.msglvl, self.maxset, self.shfscl, self.norm]
        for (option, value) in zip(self.options, self.values):
            list_fields += [option + '=' + str(value)]
        return list_fields

    def repr_fields(self):
        msglvl = set_blank_if_default(self.msglvl, 0)
        list_fields = ['EIGRL', self.sid, self.v1, self.v2, self.nd, msglvl,
                       self.maxset, self.shfscl, self.norm]
        for (option, value) in zip(self.options, self.values):
            list_fields += [option + '=' + str(value)]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
