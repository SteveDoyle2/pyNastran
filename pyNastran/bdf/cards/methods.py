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
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank, components,
    components_or_blank, integer_double_string_or_blank, blank, interpret_value)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class Method(BaseCard):
    """
    Generic class for all methods.
    Part of self.methods
    """
    def __init__(self, card, data):
        pass


class EIGB(Method):
    """
    Defines data needed to perform buckling analysis
    """
    type = 'EIGB'

    def __init__(self, card=None, data=None, comment=''):
        Method.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Set identification number. (Unique Integer > 0)
            self.sid = integer(card, 1, 'sid')

            #: Method of eigenvalue extraction. (Character: 'INV' for inverse
            #: power method or 'SINV' for enhanced inverse power method.)
            #: apparently it can also be blank...
            self.method = string_or_blank(card, 2, 'method')

            if self.method not in ['INV', 'SINV', None]:
                msg = 'method must be INV or SINV.  method=|%s|' % self.method
                raise RuntimeError(msg)

            #: Eigenvalue range of interest. (Real, L1 < L2)
            self.L1 = double(card, 3, 'L1')
            self.L2 = double(card, 4, 'L2')
            if not self.L1 < self.L2:
                msg = 'L1=%s L2=%s; L1<L2 is requried' % (self.L1, self.L2)
                raise RuntimeError(msg)

            #: Estimate of number of roots in positive range not used for
            #: METHOD = 'SINV'. (Integer > 0)
            self.nep = integer_or_blank(card, 5, 'nep', 0)

            #: Desired number of positive and negative roots.
            #: (Integer>0; Default = 3*NEP)
            self.ndp = integer_or_blank(card, 6, 'ndp', 3 * self.nep)
            self.ndn = integer_or_blank(card, 7, 'ndn', 3 * self.nep)

            #: Method for normalizing eigenvectors.
            #: ('MAX' or 'POINT';Default='MAX')
            self.norm = string_or_blank(card, 9, 'norm', 'MAX')
            if self.norm == 'POINT':
                self.G = integer(card, 10, 'G')
                self.C = components(card, 11, 'C')
            else:
                self.G = integer_or_blank(card, 10, 'G')
                self.C = components_or_blank(card, 11, 'C')
            assert len(card) <= 12, 'len(EIGB card) = %i' % len(card)
        else:
            raise NotImplementedError('EIGB')

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

    def __init__(self, card=None, data=None, comment=''):
        Method.__init__(self, card, data)
        if comment:
            self._comment = comment
        # CLAN
        self.mblkszs = []
        self.iblkszs = []
        self.ksteps = []
        self.NJIs = []

        # HESS
        self.alphaBjs = []
        self.omegaBjs = []
        self.LJs = []
        self.NEJs = []
        self.NDJs = []

        if card:
            #: Set identification number. (Unique Integer > 0)
            self.sid = integer(card, 1, 'sid')
            #: Method of complex eigenvalue extraction
            self.method = string(card, 2, 'method')
            assert self.method in ['INV', 'HESS', 'CLAN', 'ISRR'], (
                'method=%s is not INV, HESS, CLAN, ISRR' % self.method)
            #: Method for normalizing eigenvectors
            self.norm = string_or_blank(card, 3, 'norm')
            if self.norm == 'POINT':
                #: Grid or scalar point identification number. Required only if
                #: NORM='POINT'. (Integer>0)
                self.G = integer(card, 4, 'G')

                #: Component number. Required only if NORM='POINT' and G is a
                #: geometric grid point. (1<Integer<6)
                self.C = components(card, 5, 'C')
            else:
                self.G = blank(card, 4, 'G')
                self.C = blank(card, 5, 'C')

            #: Convergence criterion. (Real > 0.0. Default values are:
            #: 10^-4 for METHOD = "INV",
            #: 10^-8 for METHOD = "CLAN",
            #: 10^-8 for METHOD = "ISRR",
            #: 10^-15 for METHOD = "HESS",
            #: E is machine dependent for METHOD = "CLAN".)
            self.E = double_or_blank(card, 6, 'E')
            self.ndo = integer_double_string_or_blank(card, 7, 'ND0')

            # ALPHAAJ OMEGAAJ ALPHABJ OMEGABJ LJ NEJ NDJ
            fields = [interpret_value(field) for field in card[9:]]
            self.alphaAjs = []
            self.omegaAjs = []
            nfields = len(fields)
            nrows = nfields // 8
            if nfields % 8 > 0:
                nrows += 1
            #if nrows == 0:
                #raise RuntimeError('invalid row count=0; nfields=%s \ncard=%s\nfields=%s' % (nfields, card, fields))

            if self.method == 'CLAN':
                self.loadCLAN(nrows, card)
            elif self.method in ['HESS', 'INV']:  # HESS, INV
                self.loadHESS_INV(nrows, card)
            elif self.method == 'ISRR':
                self._load_isrr(nrows, card)
            else:
                msg = 'invalid EIGC method...method=%r' % self.method
                raise RuntimeError(msg)
            #assert card.nFields() < 8, 'card = %s' % card
        else:
            raise NotImplementedError('EIGC')

    def _load_isrr(self, nrows, card):
        assert nrows == 1, card
        for irow in range(nrows):
            i = 9 + 8 * irow
            self.shift_r1 = double_or_blank(card, i, 'SHIFT_R1', 0.0)
            self.shift_i1 = double_or_blank(card, i + 1, 'SHIFT_I1', 0.0)
            #2
            #3
            #4
            self.isrr_flag = integer_or_blank(card, i + 5, 'ISRR_FLAG', 0)
            self.nd1 = integer(card, i + 6, 'ND1')

    def loadCLAN(self, nrows, card):
        for irow in range(nrows):
            #NDJ_default = None
            i = 9 + 8 * irow
            self.alphaAjs.append(
                double_or_blank(card, i, 'alpha' + str(irow), 0.0))
            self.omegaAjs.append(
                double_or_blank(card, i + 1, 'omega' + str(irow), 0.0))
            self.mblkszs.append(
                double_or_blank(card, i + 2, 'mblock' + str(irow), 7))

            self.iblkszs.append(
                integer_or_blank(card, i + 3, 'iblksz' + str(irow), 2))
            self.ksteps.append(
                integer_or_blank(card, i + 4, 'kstep' + str(irow), 5))
            self.NJIs.append(
                integer(card, i + 6, 'NJI' + str(irow)))

    def loadHESS_INV(self, nrows, card):
        alphaOmega_default = None
        LJ_default = None
        if self.method == 'INV':
            alphaOmega_default = 0.0
            LJ_default = 1.0

        for irow in range(nrows):
            NEj = integer(card, 9 + 7 * irow + 5, 'NE%s' % str(irow))
            NDJ_default = None
            if self.method == 'INV':
                NDJ_default = 3 * NEj

            i = 9 + 8 * irow
            self.alphaAjs.append(
                double_or_blank(card, i, 'alphaA' + str(irow), alphaOmega_default))
            self.omegaAjs.append(
                double_or_blank(card, i + 1, 'omegaA' + str(irow), alphaOmega_default))
            self.alphaBjs.append(
                double_or_blank(card, i + 2, 'alphaB' + str(irow), alphaOmega_default))
            self.omegaBjs.append(
                double_or_blank(card, i + 3, 'omegaB' + str(irow), alphaOmega_default))
            self.LJs.append(
                double_or_blank(card, i + 4, 'LJ' + str(irow), LJ_default))
            self.NEJs.append(
                integer(card, i + 5, 'NEJ' + str(irow)))
            self.NDJs.append(
                integer_or_blank(card, i + 6, 'NDJ' + str(irow), NDJ_default))

    def cross_reference(self, model):
        pass

    def rawMethod(self):
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
        else:
            msg = 'invalid EIGC method...method=%r' % self.method
            raise RuntimeError(msg)
        return list_fields

    def reprMethod(self):
        return self.rawMethod()

    def raw_fields(self):
        list_fields = ['EIGC', self.sid, self.method, self.norm, self.G, self.C,
                       self.E, self.ndo, None]
        list_fields += self.rawMethod()
        return list_fields

    def repr_fields(self):
        if self.E is None:
            E = None
        else:
            E = str(self.E)
        list_fields = ['EIGC', self.sid, self.method, self.norm, self.G, self.C,
                       E, self.ndo, None]
        list_fields += self.reprMethod()
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

    def __init__(self, card=None, data=None, comment=''):
        Method.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Set identification number. (Unique Integer > 0)
            self.sid = integer(card, 1, 'sid')

            #: Coordinates of point in complex plane. (Real)
            self.alpha1 = double(card, 2, 'alpha1')
            #: Coordinates of point in complex plane. (Real)
            self.omega1 = double(card, 3, 'omega1')
            #: Multiplicity of complex root at pole defined by point at ALPHAi
            #: and OMEGAi
            self.m1 = integer(card, 4, 'm1')

            #: Coordinates of point in complex plane. (Real)
            self.alpha2 = double(card, 5, 'alpha2')
            #: Coordinates of point in complex plane. (Real)
            self.omega2 = double(card, 6, 'omega2')
            #: Multiplicity of complex root at pole defined by point at ALPHAi
            #: and OMEGAi
            self.m2 = integer(card, 7, 'm2')
            assert len(card) == 8, 'len(EIGP card) = %i' % len(card)
        else:
            raise NotImplementedError('EIGP')

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

    def __init__(self, card=None, data=None, comment=''):
        Method.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Set identification number. (Unique Integer > 0)
            self.sid = integer(card, 1, 'sid')

            #: Method of eigenvalue extraction. (Character: 'INV' for inverse
            #: power method or 'SINV' for enhanced inverse power method.)
            self.method = string_or_blank(card, 2, 'method', 'LAN')
            assert self.method in ['LAN', 'AHOU', 'INV', 'SINV', 'GIV', 'MGIV', 'HOU', 'MHOU', 'AGIV'], 'method=%s' % self.method

            #: Frequency range of interest
            self.f1 = double_or_blank(card, 3, 'f1')
            self.f2 = double_or_blank(card, 4, 'f2')

            #: Estimate of number of roots in range (Required for
            #: METHOD = 'INV'). Not used by 'SINV' method.
            self.ne = integer_or_blank(card, 5, 'ne')

            #: Desired number of roots (default=600 for SINV 3*ne for INV)
            if self.method in ['SINV']:
                self.nd = integer_or_blank(card, 6, 'nd', 600)
            if self.method in ['INV']:
                self.nd = integer_or_blank(card, 6, 'nd', 3 * self.ne)
            elif self.method in ['GIV', 'MGIV', 'HOU', 'MHOU']:
                self.nd = integer_or_blank(card, 6, 'nd', 0)
            else:
                self.nd = integer(card, 6, 'nd')
            #: Method for normalizing eigenvectors. ('MAX' or 'POINT';
            #: Default='MAX')
            self.norm = string_or_blank(card, 9, 'norm', 'MASS')
            assert self.norm in ['POINT', 'MASS', 'MAX']

            if self.method == 'POINT':
                #: Grid or scalar point identification number. Required only if
                #: NORM='POINT'. (Integer>0)
                self.G = integer(card, 10, 'G')

                #: Component number. Required only if NORM='POINT' and G is a
                #: geometric grid point. (1<Integer<6)
                self.C = components(card, 11, 'C')
            else:
                self.G = blank(card, 10, 'G')
                self.C = blank(card, 11, 'C')
            assert len(card) <= 12, 'len(EIGR card) = %i' % len(card)
        else:
            raise NotImplementedError('EIGR')

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

    def __init__(self, card=None, data=None, sol=None, comment=''):
        Method.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Set identification number. (Unique Integer > 0)
            self.sid = integer(card, 1, 'sid')
            #: For vibration analysis: frequency range of interest. For
            #: buckling analysis: eigenvalue range of interest. See Remark 4.
            #: (Real or blank, -5 10e16 <= V1 < V2 <= 5.10e16)
            self.v1 = double_or_blank(card, 2, 'v1')
            self.v2 = double_or_blank(card, 3, 'v2')
            #: Number of roots desired
            self.nd = integer_or_blank(card, 4, 'nd')
            #: Diagnostic level. (0 < Integer < 4; Default = 0)
            self.msglvl = integer_or_blank(card, 5, 'msglvl', 0)
            #: Number of vectors in block or set. Default is machine dependent
            self.maxset = integer_or_blank(card, 6, 'maxset')
            #: Estimate of the first flexible mode natural frequency
            #: (Real or blank)
            self.shfscl = double_or_blank(card, 7, 'shfscl')
            #: Method for normalizing eigenvectors (Character: 'MASS' or 'MAX')
            self.norm = string_or_blank(card, 8, 'norm')

            optionValues = [interpret_value(field) for field in card[9:]]
            self.options = []
            self.values = []
            #print "optionValues = ",optionValues
            for optionValue in optionValues:
                #print "optionValue = ",optionValue
                (option, value) = optionValue.split('=')
                self.options.append(option)
                self.values.append(value)

            #: Method for normalizing eigenvectors
            if sol in [103, 115, 146]:
                # normal modes,cyclic normal modes, flutter
                self.norm = string_or_blank(card, 8, 'norm', 'MASS')
            elif sol in [105, 110, 111, 116]:
                # buckling, modal complex eigenvalues,
                # modal frequency response,cyclic buckling
                self.norm = string_or_blank(card, 8, 'norm', 'MAX')
            else:
                self.norm = string_or_blank(card, 8, 'norm')

            #assert len(card) <= 9, 'len(EIGRL card) = %i' % len(card)
            assert len(card) <= 10, 'len(EIGRL card) = %i' % len(card)

            #msg = 'norm=%s sol=%s' % (self.norm, sol)
            #assert self.norm in ['MASS', 'MAX'],msg
            #assert card.nFields()<9,'card = %s' %(card.fields(0))
        else:
            raise NotImplementedError('EIGRL')

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
