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
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string, string_or_blank,
    parse_components, components_or_blank, integer_double_string_or_blank, blank,
    interpret_value)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Method(BaseCard):
    """
    Generic class for all methods.
    Part of self.methods

    """
    def __init__(self):
        pass


class EIGB(Method):
    """Defines data needed to perform buckling analysis"""
    type = 'EIGB'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        method = 'INV'
        G = 1
        C = 1
        norm = 'MAX'
        L1 = 1.0
        L2 = 2.0
        nep = 10
        ndp = 20
        ndn = 30
        return EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, G, C, comment='')

    def __init__(self, sid, method, L1, L2, nep, ndp, ndn, norm, G, C, comment=''):
        Method.__init__(self)
        if comment:
            self.comment = comment
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
        """
        Adds an EIGB card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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
            C = parse_components(card, 11, 'C')
        else:
            G = integer_or_blank(card, 10, 'G')
            C = components_or_blank(card, 11, 'C')
        assert len(card) <= 12, 'len(EIGB card) = %i\ncard=%s' % (len(card), card)
        return EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                    comment=comment)

    def cross_reference(self, model: BDF) -> None:
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGC(Method):
    """
    Defines data needed to perform complex eigenvalue analysis
    .. todo: not done

    ``inverse power``

    +------+---------+---------+---------+---------+---------+---------+-----+
    |   1  |    2    |    3    |    4    |    5    |    6    |   7     |  8  |
    +======+=========+=========+=========+=========+=========+=========+=====+
    | EIGC |   SID   | METHOD  |         |         |         |   EPS   | ND0 |
    +------+---------+---------+---------+---------+---------+---------+-----+
    |      | ALPHAAj | OMEGAAj | ALPHABj | OMEGABj |   Lj    |   NEj   | NDj |
    +------+---------+---------+---------+---------+---------+---------+-----+

    ``complex Lanczos``

    +------+---------+---------+---------+---------+---------+---------+-----+
    |   1  |    2    |    3    |    4    |    5    |    6    |   7     |  8  |
    +======+=========+=========+=========+=========+=========+=========+=====+
    |      | SHIFTRj | SHIFTIj | MBLKSZj | IBLKSZj | KSTEPSj |   NDj   |     |
    +------+---------+---------+---------+---------+---------+---------+-----+

    ``iterative Schur-Rayleigh-Ritz``

    +------+---------+---------+---------+---------+---------+---------+-----+
    |   1  |    2    |    3    |    4    |    5    |    6    |   7     |  8  |
    +======+=========+=========+=========+=========+=========+=========+=====+
    |      | SHIFTR1 | SHIFTI1 |         |         |         | ISRRFLG | ND1 |
    +------+---------+---------+---------+---------+---------+---------+-----+

    """
    type = 'EIGC'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        method = 'CLAN'
        grid = 1
        component = 1
        epsilon = 0.1
        neigenvalues = 10
        return EIGC(sid, method, grid, component, epsilon, neigenvalues,
                    norm='MAX', mblkszs=None, iblkszs=None, ksteps=None,
                    NJIs=None, alphaAjs=None, omegaAjs=None, alphaBjs=None,
                    omegaBjs=None, LJs=None, NEJs=None, NDJs=None,
                    shift_r1=None, shift_i1=None, isrr_flag=None, nd1=None, comment='')

    def __init__(self, sid, method, grid, component, epsilon, neigenvalues,
                 norm='MAX', # common
                 mblkszs=None, iblkszs=None, ksteps=None, NJIs=None, # CLAN
                 alphaAjs=None, omegaAjs=None, alphaBjs=None, omegaBjs=None, # HESS/INV
                 LJs=None, NEJs=None, NDJs=None, # HESS/INV
                 shift_r1=None, shift_i1=None, isrr_flag=None, nd1=None, # ISRR
                 comment=''):
        """
        Creates a EIGC card, which is required for a SOL 107 analysis

        Parameters
        ----------
        sid : int
            CMETHOD id in the case control deck
        method : str
           Method of complex eigenvalue extraction
             MSC 2014 = [INV, HESS, CLAN, IRAM]
             NX 8.5 = [INV, HESS, CLAN, ISRR]
             Autodesk 2015 = [ARNO, HESS, CLAN]
             INV  : Inverse Power
             IRAM : Implicitly Restarted Arnoldi method
             ISRR : Iterative Schur-Rayleigh-Ritz method
             CLAN : Complex Lanczos.  For linear perturbation of ANALYSIS=DCEIG
                    with large displacement, CLAN is recommended.
             HESS : Upper Hessenberg. For linear perturbation of ANALYSIS=DCEIG
                    with large displacement, please don't use HESS.
             ARNO: ???
        norm : str; default='MAX'
            Method for normalizing eigenvectors
            valid_norm = {MAX, POINT}
        grid : int
            GRID/SPOINT id
            Required if norm='POINT'
        component : int
            Required if norm='POINT'
        epsilon : float
        neigenvalues : int
            Number of Eigenvalues
        mblkszs : List[float]; default=None
            used by CLAN
        iblkszs : List[int]; default=None
            used by CLAN
        ksteps : List[int]; default=None
            used by CLAN
        NJIs : List[int]; default=None
            used by CLAN
        alphaAjs : List[float]; default=None
            used by HESS/INV
        omegaAjs : List[float]; default=None
            used by HESS/INV
        alphaBjs : List[float]; default=None
            used by HESS/INV
        omegaBjs : List[float]; default=None
            used by HESS/INV
        LJs : List[float]; default=None
            used by HESS/INV
        NEJs : List[int]; default=None
            used by HESS/INV
        NDJs : List[int]; default=None
            used by HESS/INV
        shift_r1 : List[float]; default=None
            used by ISSR
        shift_i1 : List[float]; default=None
            used by ISSR
        isrr_flag : List[int]; default=None
            used by ISSR
        nd1 : List[int]; default=None
            used by ISSR
        comment : str; default=''
            a comment for the card

        """
        Method.__init__(self)
        if comment:
            self.comment = comment

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
        self.G = grid


        #: Component number. Required only if NORM='POINT' and G is a
        #: geometric grid point. (1<Integer<6)
        self.C = component

        #: Convergence criterion. (Real > 0.0. Default values are:
        #: 10^-4 for METHOD = "INV",
        #: 10^-8 for METHOD = "CLAN",
        #: 10^-8 for METHOD = "ISRR",
        #: 10^-15 for METHOD = "HESS",
        #: E is machine dependent for METHOD = "CLAN".)
        self.epsilon = epsilon

        #Number of eigenvalues and/or eigenvectors desired. See Remark
        #3. (Integer > 0 or blank; No default)
        self.neigenvalues = neigenvalues

        # CLAN
        if mblkszs is None:
            mblkszs = []
        if iblkszs is None:
            iblkszs = []
        if ksteps is None:
            ksteps = []
        if NJIs is None:
            NJIs = []

        self.mblkszs = mblkszs
        self.iblkszs = iblkszs
        self.ksteps = ksteps
        self.NJIs = NJIs

        # HESS
        if alphaBjs is None:
            alphaBjs = []
        if omegaBjs is None:
            omegaBjs = []
        self.alphaBjs = alphaBjs
        self.omegaBjs = omegaBjs

        if LJs is None:
            LJs = []
        self.LJs = LJs
        if NEJs is None:
            NEJs = []
        self.NEJs = NEJs
        if NDJs is None:
            NDJs = []
        self.NDJs = NDJs

        if alphaAjs is None:
            alphaAjs = []
        if omegaAjs is None:
            omegaAjs = []
        self.alphaAjs = alphaAjs
        self.omegaAjs = omegaAjs

        #----------
        # ISRR
        self.shift_r1 = shift_r1
        self.shift_i1 = shift_i1
        self.isrr_flag = isrr_flag
        self.nd1 = nd1

    def validate(self):
        assert self.norm in ['MAX', 'POINT'], 'norm=%r' % self.norm
        nalpha_a = len(self.alphaAjs)
        assert nalpha_a == len(self.omegaAjs), 'alphaAjs=%s omegaAj=%s' % (self.alphaAjs, self.omegaAjs)
        if self.method in ['HESS', 'INV']:
            assert nalpha_a == len(self.alphaBjs), 'alphaAjs=%s alphaBj=%s' % (self.alphaAjs, self.alphaBjs)
            #assert nalpha_a == len(self.omegaBjs), 'alphaAjs=%s omegaBjs=%s' % (self.alphaAjs, self.omegaBjs)
            assert nalpha_a == len(self.LJs), 'alphaAjs=%s LJs=%s' % (self.alphaAjs, self.LJs)
            assert nalpha_a == len(self.NEJs), 'alphaAjs=%s NEJs=%s' % (self.alphaAjs, self.NEJs)
            assert nalpha_a == len(self.NDJs), 'alphaAjs=%s NDJs=%s' % (self.alphaAjs, self.NDJs)
        elif self.method == 'CLAN':
            if nalpha_a == len(self.alphaBjs):
                assert nalpha_a == len(self.alphaBjs), f'nalpha_a={nalpha_a} nalpha_b={nalpha_b}'
                assert nalpha_a == len(self.omegaBjs), f'nalpha_a={nalpha_a} nomega_b={len(self.omegaBjs)}'
                assert nalpha_a == len(self.LJs)
                assert nalpha_a == len(self.NEJs)
                assert nalpha_a == len(self.NDJs)
            else:
                assert nalpha_a == len(self.omegaAjs)
                assert nalpha_a == len(self.mblkszs), 'alphaAjs=%s mblkszs=%s' % (self.alphaAjs, self.mblkszs)
                assert nalpha_a == len(self.iblkszs)
                assert nalpha_a == len(self.ksteps)
                assert nalpha_a == len(self.NJIs)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an EIGC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        method = string(card, 2, 'method')
        assert method in ['ARNO', 'INV', 'HESS', 'CLAN', 'ISRR', 'IRAM', 'DET'], (
            'method=%s is not ARNO, INV, HESS, CLAN, ISRR, IRAM, DET' % method)

        norm = string_or_blank(card, 3, 'norm', 'MAX')
        if norm == 'POINT':
            grid = integer(card, 4, 'G')
            component = parse_components(card, 5, 'C')
        else:
            grid = blank(card, 4, 'G')
            component = blank(card, 5, 'C')

        epsilon = double_or_blank(card, 6, 'epsilon')
        neigenvalues = integer_double_string_or_blank(card, 7, 'ND0/neigenvalues')

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
            out = _load_clan(nrows, card)
            (alphaAjs, omegaAjs, mblkszs, iblkszs, ksteps, NJIs,
             alphaBjs, omegaBjs, LJs, NEJs, NDJs) = out
        elif method in ['HESS', 'INV', 'DET']:  # HESS, INV
            alphaAjs, omegaAjs, alphaBjs, omegaBjs, LJs, NEJs, NDJs = _load_hess_inv(
                nrows, method, card)
        elif method == 'ISRR':
            shift_r1, shift_i1, isrr_flag, nd1 = _load_isrr(nrows, card)
        else:
            raise RuntimeError(f'invalid EIGC method...method={method!r}')
        #assert card.nfields() < 8, 'card = %s' % card
        return EIGC(sid, method, grid, component, epsilon, neigenvalues,
                    norm, # common
                    mblkszs, iblkszs, ksteps, NJIs, # CLAN
                    alphaAjs, omegaAjs, alphaBjs, omegaBjs, LJs, NEJs, NDJs, # HESS/INV
                    shift_r1, shift_i1, isrr_flag, nd1, # ISRR
                    comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def raw_method(self):
        list_fields = []
        if self.method in ['HESS', 'INV', 'DET']:
            for (alphaA, omegaA, alphaB, omegaB, Lj, NEj, NDj) in zip(
                    self.alphaAjs, self.omegaAjs, self.alphaBjs, self.omegaBjs,
                    self.LJs, self.NEJs, self.NDJs):
                alphaA = set_blank_if_default(alphaA, 0.0)
                omegaA = set_blank_if_default(omegaA, 0.0)
                alphaB = set_blank_if_default(alphaB, 0.0)
                omegaB = set_blank_if_default(omegaB, 0.0)
                list_fields += [alphaA, omegaA, alphaB, omegaB, Lj, NEj, NDj, None]

        elif self.method == 'CLAN':
            nalpha_a = len(self.alphaAjs)
            assert nalpha_a == len(self.omegaAjs)
            if nalpha_a == len(self.alphaBjs):  # pragma:no cover
                assert nalpha_a == len(self.alphaBjs), f'nalpha_a={nalpha_a} nalpha_b={nalpha_b}'
                assert nalpha_a == len(self.omegaBjs), f'nalpha_a={nalpha_a} nomega_b={len(self.omegaBjs)}'
                assert nalpha_a == len(self.LJs)
                assert nalpha_a == len(self.NEJs)
                assert nalpha_a == len(self.NDJs)
                for (alphaA, omegaA, alphaB, omegaB, Lj, Nej, Ndj) in zip(
                        self.alphaAjs, self.omegaAjs,
                        self.alphaBjs, self.omegaBjs,
                        self.LJs, self.NEJs, self.NDJs):
                    #alphaA = set_blank_if_default(alphaA, 0.0)
                    #omegaA = set_blank_if_default(omegaA, 0.0)
                    #mblksz = set_blank_if_default(mblksz, 7)
                    #iblksz = set_blank_if_default(iblksz, 2)
                    #kstep = set_blank_if_default(kstep, 5)
                    list_fields += [alphaA, omegaA, alphaB, omegaB, Lj,
                                    Nej, Ndj, None]
            else:
                assert nalpha_a == len(self.mblkszs)
                assert nalpha_a == len(self.iblkszs)
                assert nalpha_a == len(self.ksteps)
                assert nalpha_a == len(self.NJIs)
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
            assert self.shift_r1 is not None, self.get_stats()
            assert len(self.shift_r1) > 0, self.get_stats()
            for shift_r1i, shift_i1i, isrr_flagi, nd1i in zip(
                self.shift_r1, self.shift_i1, self.isrr_flag, self.nd1):
                list_fields += [shift_r1i, shift_i1i, None, None, None, isrr_flagi, nd1i, None]

        else:
            raise RuntimeError(f'invalid EIGC method.  method={self.method!r} '
                               'expected=[HESS, INV, DET, CLAN, ISRR]')
        return list_fields

    def repr_method(self):
        return self.raw_method()

    def raw_fields(self):
        list_fields = ['EIGC', self.sid, self.method, self.norm, self.G, self.C,
                       self.epsilon, self.neigenvalues, None]
        list_fields += self.raw_method()
        return list_fields

    def repr_fields(self):
        if self.epsilon is None:
            epsilon = None
        else:
            epsilon = self.epsilon
        list_fields = ['EIGC', self.sid, self.method, self.norm, self.G, self.C,
                       epsilon, self.neigenvalues, None]
        list_fields += self.repr_method()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGP(Method):
    """
    Defines poles that are used in complex eigenvalue extraction by the
    Determinant method.

    +------+-------+--------+--------+-------+--------+--------+-----+
    |   1  |   2   |   3    |   4    |   5   |   6    |   7    |  8  |
    +======+=======+========+========+=======+========+========+=====+
    | EIGP |  SID  | ALPHA1 | OMEGA1 |   M1  | ALPHA2 | OMEGA2 |  M2 |
    +------+-------+--------+--------+-------+--------+--------+-----+
    | EIGP |  15   |  -5.2  |  0.0   |   2   |  6.3   |  5.5   |  3  |
    +------+-------+--------+--------+-------+--------+--------+-----+

    """
    type = 'EIGP'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        alpha1 = 1.
        omega1 = 1.
        m1 = 1.
        alpha2 = 1.
        omega2 = 1.
        m2 = 1.
        return EIGP(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment='')

    def __init__(self, sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=''):
        Method.__init__(self)
        if comment:
            self.comment = comment
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
        """
        Adds an EIGPX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')

        alpha1 = double(card, 2, 'alpha1')
        omega1 = double(card, 3, 'omega1')
        m1 = integer(card, 4, 'm1')

        alpha2 = double(card, 5, 'alpha2')
        omega2 = double(card, 6, 'omega2')
        m2 = integer(card, 7, 'm2')
        assert len(card) == 8, 'len(EIGP card) = %i\ncard=%s' % (len(card), card)
        return EIGP(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def raw_fields(self):
        list_fields = ['EIGP', self.sid, self.alpha1, self.omega1, self.m1,
                       self.alpha2, self.omega2, self.m2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGR(Method):
    """
    Defines data needed to perform real eigenvalue analysis
    """
    type = 'EIGR'
    allowed_methods = [
        'LAN', 'AHOU', # recommended
        'INV', 'SINV', 'GIV', 'MGIV', 'HOU', 'MHOU', 'AGIV' # obsolete
    ]

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        return EIGR(sid, method='LAN', f1=None, f2=None, ne=None, nd=None,
                    norm='MASS', G=None, C=None, comment='')

    def __init__(self, sid, method='LAN', f1=None, f2=None, ne=None, nd=None,
                 norm='MASS', G=None, C=None, comment=''):
        """
        Adds a EIGR card

        Parameters
        ----------
        sid : int
            method id
        method : str; default='LAN'
            eigenvalue method
            recommended: {LAN, AHOU}
            obsolete : {INV, SINV, GIV, MGIV, HOU, MHOU, AGIV}
        f1 / f2 : float; default=None
            lower/upper bound eigenvalue
        f2 : float; default=None
            upper bound eigenvalue
        ne : int; default=None
            estimate of number of roots (used for INV)
        nd : int; default=None
            desired number of roots
        msglvl : int; default=0
            debug level; 0-4
        maxset : int; default=None
            Number of vectors in block or set
        shfscl : float; default=None
            estimate of first flexible mode natural frequency
        norm : str; default=None
            {MAX, MASS, AF, POINT}
            default=MASS (NX)
        G : int; default=None
            node id for normalization; only for POINT
        C : int; default=None
            component for normalization (1-6); only for POINT
        comment : str; default=''
            a comment for the card

        """
        Method.__init__(self)
        if comment:
            self.comment = comment
        if G == 0:
            G = None
        if C == 0:
            C = None

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
        """
        Adds an EIGR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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
            ne = integer(card, 5, 'ne')
            nd = integer_or_blank(card, 6, 'nd', 3 * ne)
        elif method in ['GIV', 'MGIV', 'HOU', 'MHOU']:
            nd = integer_or_blank(card, 6, 'nd', 0)
        else:
            nd = integer(card, 6, 'nd')
        norm = string_or_blank(card, 9, 'norm', 'MASS')

        if norm == 'POINT':
            G = integer(card, 10, 'G')
            C = parse_components(card, 11, 'C')
        else:
            G = blank(card, 10, 'G')
            C = blank(card, 11, 'C')
        assert len(card) <= 12, 'len(EIGR card) = %i\ncard=%s' % (len(card), card)
        return EIGR(sid, method, f1, f2, ne, nd, norm, G, C, comment=comment)

    def cross_reference(self, model: BDF) -> None:
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EIGRL(Method):
    """
    Defines data needed to perform real eigenvalue (vibration or buckling)
    analysis with the Lanczos method

    +-------+-----+----+----+----+--------+--------+--------+------+
    |   1   |  2  |  3 |  4 |  5 |    6   |    7   |    8   |   9  |
    +=======+=====+====+====+====+========+========+========+======+
    | EIGRL | SID | V1 | V2 | ND | MSGLVL | MAXSET | SHFSCL | NORM |
    +-------+-----+----+----+----+--------+--------+--------+------+
    |        option_1 = value_1 option_2 = value_2, etc.           |
    +--------------------------------------------------------------+

    """
    type = 'EIGRL'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        return EIGRL(sid, v1=None, v2=None, nd=None, msglvl=0, maxset=None,
                     shfscl=None, norm=None, options=None, values=None, comment='')

    def __init__(self, sid, v1=None, v2=None, nd=None, msglvl=0, maxset=None, shfscl=None,
                 norm=None, options=None, values=None, comment=''):
        """
        Adds an EIGRL card

        Parameters
        ----------
        sid : int
            method id
        v1 : float; default=None
            lower bound eigenvalue
        v2 : float; default=None
            upper bound eigenvalue
        nd : int
            number of roots
        msglvl : int; default=0
            debug level; 0-4
        maxset : int; default=None
            Number of vectors in block or set
        shfscl : float; default=None
            estimate of first flexible mode natural frequency
        norm : str; default=None
            {MAX, MASS, AF}
        options : ???; default=None -> []
            ???
        values : ???; default=None -> []
            ???
        comment : str; default=''
            a comment for the card

        """
        Method.__init__(self)
        if comment:
            self.comment = comment
        if options is None:
            options = []
        if values is None:
            values = []

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

    def validate(self):
        assert self.norm in [None, 'MAX', 'MASS', 'AF'], 'norm=%r' % self.norm
        assert self.msglvl in [0, 1, 2, 3, 4], 'msglvl=%r' % self.msglvl
        if len(self.options) != len(self.values):
            raise RuntimeError('len(options) != len(values); noptions=%s nvalues=%s\n'
                               'options=%s values=%s' % (len(self.options), len(self.values),
                                                         self.options, self.values))
        for option, value in zip(self.options, self.values):
            if option == 'NORM':
                assert value in ['MAX', ], 'option=%r value=%r' % (option, value)
            elif option == 'ALPH':
                # float
                pass
            elif option == 'NUMS':
                # integer
                pass
            else:
                raise NotImplementedError('option=%r value=%r' % (option, value))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an EIGRL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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


    def cross_reference(self, model: BDF) -> None:
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


def _load_isrr(nrows, card):
    """loads the iterative Schur-Rayleigh-Ritz"""
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

def _load_clan(nrows, card):
    """loads complex Lanczos"""
    alphaAjs = []
    omegaAjs = []
    mblkszs = []
    iblkszs = []
    ksteps = []
    NJIs = []

    alphaBjs = []
    omegaBjs = []
    ljs = []
    nejs = []
    ndjs = []
    is_nej = None
    for irow in range(nrows):
        #NDJ_default = None
        i = 9 + 8 * irow
        alphaAjs.append(
            double_or_blank(card, i, 'alpha' + str(irow), 0.0))
        omegaAjs.append(
            double_or_blank(card, i + 1, 'omega' + str(irow), 0.0))

        nej_blank = integer_or_blank(card, i + 6, 'NEJ_blank')
        if nej_blank is not None and 0:  # pragma: no cover
            assert is_nej in [True, None], is_nej
            is_nej = True
            # ALPHAAJ OMEGAAJ ALPHABJ OMEGABJ LJ NEJ NDJ
            assert isinstance(nej_blank, int), nej_blank
            alpha_bj = double(card, i + 2, 'alpha_bj' + str(irow))
            omega_bj = double(card, i + 3, 'omega_bj' + str(irow))
            lj = double_or_blank(card, i + 4, 'LJ' + str(irow), 1.0)
            nej = integer_or_blank(card, i + 5, 'NEJ' + str(irow))
            ndj = integer(card, i + 6, 'NDJ' + str(irow))
            alphaBjs.append(alpha_bj)
            omegaBjs.append(omega_bj)
            ljs.append(lj)
            nejs.append(nej)
            ndjs.append(ndj)
        else:
            assert is_nej in [False, None], is_nej
            is_nej = False
            # ALPHAAJ OMEGAAJ MBLKSZ IBLKSZ KSTEPS blank NJi
            mblock_size = double_or_blank(card, i + 2, 'mblock' + str(irow), 7)

            # iblkszs is an integer, but entered as a float...
            iblock_size = double_or_blank(card, i + 3, 'iblksz' + str(irow), 2.0)
            kstep = integer_or_blank(card, i + 4, 'kstep' + str(irow), 5)
            nji = integer(card, i + 6, 'NJI' + str(irow))

            mblkszs.append(mblock_size)
            iblkszs.append(iblock_size)
            ksteps.append(kstep)
            NJIs.append(nji)

    out = (
        alphaAjs, omegaAjs, mblkszs, iblkszs, ksteps, NJIs,
        alphaBjs, omegaBjs, ljs, nejs, ndjs,
    )
    return out

def _load_hess_inv(nrows, method, card):
    """loads inverse power"""
    alpha_omega_default = None
    LJ_default = None
    if method == 'INV':
        alpha_omega_default = 0.0
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
            double_or_blank(card, i, 'alphaA' + str(irow), alpha_omega_default))
        omegaAjs.append(
            double_or_blank(card, i + 1, 'omegaA' + str(irow), alpha_omega_default))
        alphaBjs.append(
            double_or_blank(card, i + 2, 'alphaB' + str(irow), alpha_omega_default))
        omegaBjs.append(
            double_or_blank(card, i + 3, 'omegaB' + str(irow), alpha_omega_default))
        LJs.append(
            double_or_blank(card, i + 4, 'LJ' + str(irow), LJ_default))
        NEJs.append(NEj)
        NDJs.append(integer_or_blank(card, i + 6, 'NDJ' + str(irow), NDJ_default))
    return alphaAjs, omegaAjs, alphaBjs, omegaBjs, LJs, NEJs, NDJs


class MODTRAK(BaseCard):
    """
    MODTRAK SID LOWRNG HIGHRNG MTFILTER
    MODTRAK 100   1      26      0.80
    """
    def __init__(self, sid, low_range, high_range, mt_filter, comment=''):
        BaseCard.__init__(self)
        self.sid = sid
        self.low_range = low_range
        self.high_range = high_range
        self.mt_filter = mt_filter

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        low_range = integer_or_blank(card, 2, 'low_range', 0)
        high_range = integer(card, 3, 'high_range')
        mt_filter = double_or_blank(card, 4, 'mt_filter', 0.9)
        return MODTRAK(sid, low_range, high_range, mt_filter, comment=comment)

    def raw_fields(self) -> List[Any]:
        list_fields = ['MODTRAK', self.sid, self.low_range, self.high_range, self.mt_filter]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        fields = self.raw_fields()
        #if size == 8:
        return self.comment + print_card_8(fields)
        #return self.comment + print_card_16(fields)
