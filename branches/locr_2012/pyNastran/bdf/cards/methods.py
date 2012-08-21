## pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from itertools import izip

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard


class Method(BaseCard):
    """
    Generic class for all methods.
    Part of self.methods
    """
    def __init__(self, card=None, data=None):
        pass


class EIGB(Method):
    """
    Defines data needed to perform buckling analysis
    """
    type = 'EIGB'

    def __init__(self, card=None, data=None):
        Method.__init__(self, card, data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## Method of eigenvalue extraction. (Character: 'INV' for inverse
            ## power method or 'SINV' for enhanced inverse power method.)
            self.method = card.field(2)
            assert self.method in ['INV', 'SINV'], 'method must be INV or SINV.  method=|%s|' % (self.method)
            ## Eigenvalue range of interest. (Real, L1 < L2)
            self.L1 = card.field(3)
            self.L2 = card.field(4)
            assert self.L1 < self.L2, 'L1=%s L2=%s; L1<L2 is requried' % (
                self.L1, self.L2)
            ## Estimate of number of roots in positive range not used for
            ## METHOD = 'SINV'. (Integer > 0)
            self.nep = card.field(5, 0)
            ## Desired number of positive and negative roots.
            ## (Integer>0; Default = 3*NEP)
            self.ndp = card.field(6, 3 * self.nep)
            self.ndn = card.field(7, 3 * self.nep)
            ## Method for normalizing eigenvectors.
            ## ('MAX' or 'POINT';Default='MAX')
            self.norm = card.field(9, 'MAX')
            self.G = card.field(10)
            self.C = card.field(11)
        else:
            raise NotImplementedError('EIGB')
        ###
        #print self.rawFields()
        #print self.reprFields()
        #print self

    def cross_reference(self, model):
        pass

    def rawFields(self):
        fields = ['EIGB', self.sid, self.method, self.L1, self.L2, self.nep,
                  self.ndp, self.ndn, None, self.norm, self.G, self.C]
        return fields

    def reprFields(self):
        #method = set_blank_if_default(self.method,'INV')
        nep = set_blank_if_default(self.nep, 0)
        ndp = set_blank_if_default(self.ndp, 3 * self.nep)
        ndn = set_blank_if_default(self.ndn, 3 * self.nep)
        #print "nep = ",self.nep,ndn
        norm = set_blank_if_default(self.norm, 'MAX')
        fields = ['EIGB', self.sid, self.method, self.L1, self.L2, nep, ndp, ndn, None,
                  norm, self.G, self.C]
        return fields


class EIGC(Method):  ## not done
    """
    Defines data needed to perform complex eigenvalue analysis
    """
    type = 'EIGC'

    def __init__(self, card=None, data=None):
        Method.__init__(self, card, data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## Method of complex eigenvalue extraction
            self.method = card.field(2)
            assert self.method in ['INV', 'HESS', 'CLAN'],(
                    'method=%s is not INV, HESS, CLAN' % (self.method))
            ## Method for normalizing eigenvectors
            self.norm = card.field(3)

            ## Grid or scalar point identification number. Required only if
            ## NORM='POINT'. (Integer>0)
            self.G = card.field(4)
            ## Component number. Required only if NORM='POINT' and G is a
            ## geometric grid point. (1<Integer<6)
            self.C = card.field(5)
            ## Convergence criterion. (Real > 0.0. Default values are:
            # 10^-4 for METHOD = "INV",
            # 10^-15 for METHOD = "HESS",
            ## E is machine dependent for METHOD = "CLAN".)
            self.E = card.field(6)
            self.ndo = card.field(7)

            # ALPHAAJ OMEGAAJ ALPHABJ OMEGABJ LJ NEJ NDJ
            fields = card.fields(9)
            self.alphaAjs = []
            self.omegaAjs = []
            nFields = len(fields)
            nRows = nFields // 8
            if nFields % 7 > 0:
                nRows += 1

            if self.method == 'CLAN':
                self.loadCLAN(nRows, card)
            elif self.method in ['HESS', 'INV']:  # HESS, INV
                self.loadHESS_INV(nRows, card)
            else:
                msg = 'invalid EIGC method...method=|%r|' % (self.method)
                raise RuntimeError(msg)
            ###
            #assert card.nFields()<8,'card = %s' %(card.fields(0))
        else:
            raise NotImplementedError('EIGC')
        ###

    def loadCLAN(self, nRows, card):
        self.mblkszs = []
        self.iblkszs = []
        self.ksteps = []
        self.NJIs = []
        for iRow in xrange(nRows):
            #NDJ_default = None
            self.alphaAjs.append(card.field(9 + 8 * iRow, 0.0))
            self.omegaAjs.append(card.field(9 + 8 * iRow + 1, 0.0))
            self.mblkszs.append(card.field(9 + 8 * iRow + 2, 7))

            #self.alphaAjs.append(card.field(9+8*iRow  ,'ALPHA%s'%(iRow)))
            #self.omegaAjs.append(card.field(9+8*iRow+1,'OMEGA%s'%(iRow)))
            #self.mblkszs.append( card.field(9+8*iRow+2,'MBLOCK%s'%(iRow)))

            self.iblkszs.append(card.field(9 + 8 * iRow + 3, 2))
            self.ksteps.append(card.field(9 + 8 * iRow + 4, 5))
            self.NJIs.append(card.field(9 + 8 * iRow + 6))

    def loadHESS_INV(self, nRows, card):
        self.alphaBjs = []
        self.omegaBjs = []
        self.LJs = []
        self.NEJs = []
        self.NDJs = []

        alphaOmega_default = None
        LJ_default = None
        if self.method == 'INV':
            alphaOmega_default = 0.0
            LJ_default = 1.0

        for iRow in xrange(nRows):
            NEj = card.field(9 + 7 * iRow + 5)
            NDJ_default = None
            if self.method == 'INV':
                NDJ_default = 3 * NEj

            self.alphaAjs.append(
                card.field(9 + 8 * iRow, alphaOmega_default))
            self.omegaAjs.append(
                card.field(9 + 8 * iRow + 1, alphaOmega_default))
            self.alphaBjs.append(
                card.field(9 + 8 * iRow + 2, alphaOmega_default))
            self.omegaBjs.append(
                card.field(9 + 8 * iRow + 3, alphaOmega_default))
            self.LJs.append(card.field(9 + 8 * iRow + 4, LJ_default))
            self.NEJs.append(card.field(9 + 8 * iRow + 5))
            self.NDJs.append(card.field(9 + 8 * iRow + 6, NDJ_default))
        ###

    def cross_reference(self, model):
        pass

    def rawMethod(self):
        fields = []
        if self.method in ['HESS', 'INV']:
            for (alphaA, omegaA, alphaB, omegaB, Lj, NEj, NDj) in izip(
                self.alphaAjs, self.omegaAjs, self.alphaBjs, self.omegaBjs,
                self.LJs, self.NEJs, self.NDJs):
                alphaA = set_blank_if_default(alphaA, 0.0)
                omegaA = set_blank_if_default(omegaA, 0.0)
                alphaB = set_blank_if_default(alphaB, 0.0)
                omegaB = set_blank_if_default(omegaB, 0.0)
                fields += [alphaA, omegaA, alphaB, omegaB, Lj, NEj, NDj, None]

        elif self.method == 'CLAN':
            for (alphaA, omegaA, mblksz, iblksz, kstep, Nj) in izip(
                 self.alphaAjs, self.omegaAjs, self.mblkszs, self.iblkszs,
                 self.ksteps, self.NJIs):
                alphaA = set_blank_if_default(alphaA, 0.0)
                omegaA = set_blank_if_default(omegaA, 0.0)
                mblksz = set_blank_if_default(mblksz, 7)
                iblksz = set_blank_if_default(iblksz, 2)
                kstep = set_blank_if_default(kstep, 5)

                fields += [alphaA, omegaA, mblksz, iblksz,
                    kstep, None, Nj, None]
        else:
            msg = 'invalid EIGC method...method=|%r|' % (self.method)
            raise RuntimeError(msg)
        return fields

    def reprMethod(self):
        return self.rawMethod()

    def rawFields(self):
        fields = ['EIGC', self.sid, self.method, self.norm, self.G,
            self.C, self.E, self.ndo, None]
        fields += self.rawMethod()
        return fields

    def reprFields(self):
        if self.E is None:
            E = None
        else:
            E = str(self.E)
        fields = ['EIGC', self.sid, self.method, self.norm, self.G,
            self.C, E, self.ndo, None]
        fields += self.reprMethod()
        return fields


class EIGR(Method):
    """
    Defines data needed to perform real eigenvalue analysis
    """
    type = 'EIGR'

    def __init__(self, card=None, data=None):
        Method.__init__(self, card, data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## Method of eigenvalue extraction. (Character: 'INV' for inverse power
            ## method or 'SINV' for enhanced inverse power method.)
            self.method = card.field(2, 'LAN')
            ## Frequency range of interest
            self.f1 = card.field(3)
            self.f2 = card.field(4)
            ## Estimate of number of roots in range (Required for
            ## METHOD = 'INV'). Not used by 'SINV' method.
            self.ne = card.field(5)
            ## Desired number of roots (default=600 for SINV 3*ne for INV)
            self.nd = card.field(6)
            ## Method for normalizing eigenvectors. ('MAX' or 'POINT';Default='MAX')
            self.norm = card.field(9, 'MASS')
            assert self.norm in ['POINT', 'MASS', 'MAX']
            ## Grid or scalar point identification number. Required only if NORM='POINT'. (Integer>0)
            self.G = card.field(10)
            ## Component number. Required only if NORM='POINT' and G is a geometric grid point. (1<Integer<6)
            self.C = card.field(11)
        else:
            raise NotImplementedError('EIGR')
        ###

    def cross_reference(self, model):
        pass

    def rawFields(self):
        fields = ['EIGR', self.sid, self.method, self.f1, self.f2, self.ne,
                  self.nd, None, None, self.norm, self.G, self.C]
        return fields

    def reprFields(self):
        method = set_blank_if_default(self.method, 'LAN')
        norm = set_blank_if_default(self.norm, 'MASS')
        fields = ['EIGR', self.sid, method, self.f1, self.f2, self.ne,
                  self.nd, None, None, norm, self.G, self.C]
        return fields


class EIGP(Method):
    """
    Defines poles that are used in complex eigenvalue extraction by the Determinant method.
    """
    type = 'EIGP'

    def __init__(self, card=None, data=None):
        Method.__init__(self, card, data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)

            ## Coordinates of point in complex plane. (Real)
            self.alpha1 = card.field(2)
            ## Coordinates of point in complex plane. (Real)
            self.omega1 = card.field(3)
            ## Multiplicity of complex root at pole defined by point at ALPHAi
            ## and OMEGAi
            self.m1 = card.field(4)

            ## Coordinates of point in complex plane. (Real)
            self.alpha2 = card.field(5)
            ## Coordinates of point in complex plane. (Real)
            self.omega2 = card.field(6)
            ## Multiplicity of complex root at pole defined by point at ALPHAi
            ## and OMEGAi
            self.m2 = card.field(7)
        else:
            raise NotImplementedError('EIGP')
        ###

    def cross_reference(self, model):
        pass

    def rawFields(self):
        fields = ['EIGP', self.alpha1, self.omega1, self.m1,
                  self.alpha2, self.omega2, self.m2]
        return fields

    def reprFields(self):
        return self.rawFields()


class EIGRL(Method):
    """
    Defines data needed to perform real eigenvalue (vibration or buckling)
    analysis with the Lanczos method
    """
    type = 'EIGRL'

    def __init__(self, card=None, data=None, sol=None):
        Method.__init__(self, card, data)
        if card:
            ## Set identification number. (Unique Integer > 0)
            self.sid = card.field(1)
            ## For vibration analysis: frequency range of interest. For buckling
            ## analysis: eigenvalue range of interest. See Remark 4.
            ## (Real or blank, -5 10e16 <= V1 < V2 <= 5.10e16)
            self.v1 = card.field(2)
            self.v2 = card.field(3)
            ## Number of roots desired
            self.nd = card.field(4)
            ## Diagnostic level. (0 < Integer < 4; Default = 0)
            self.msglvl = card.field(5, 0)
            ## Number of vectors in block or set. Default is machine dependent
            self.maxset = card.field(6)
            ## Estimate of the first flexible mode natural frequency
            ## (Real or blank)
            self.shfscl = card.field(7)
            ## Method for normalizing eigenvectors (Character: 'MASS' or 'MAX')
            self.norm = card.field(8)

            optionValues = card.fields(9)
            self.options = []
            self.values = []
            #print "optionValues = ",optionValues
            for optionValue in optionValues:
                #print "optionValue = ",optionValue
                (option, value) = optionValue.split('=')
                self.options.append(option)
                self.values.append(value)

            ## Method for normalizing eigenvectors
            if sol in [103, 115, 146]:  # normal modes,cyclic normal modes, flutter
                self.norm = card.field(8, 'MASS')
            elif sol in [105, 110, 111, 116]:  # buckling, modal complex eigenvalues,modal frequency response,cyclic buckling
                self.norm = card.field(8, 'MAX')
            else:
                self.norm = card.field(8)
            #assert self.norm in ['MASS', 'MAX'],'norm=%s sol=%s' %(self.norm,sol)
            #assert card.nFields()<9,'card = %s' %(card.fields(0))
        else:
            raise NotImplementedError('EIGRL')
        ###

    def cross_reference(self, model):
        pass

    def rawFields(self):
        fields = ['EIGRL', self.sid, self.v1, self.v2, self.nd,
                  self.msglvl, self.maxset, self.shfscl, self.norm]
        for (option, value) in izip(self.options, self.values):
            fields += [option + '=' + str(value)]
        return fields

    def reprFields(self):
        msglvl = set_blank_if_default(self.msglvl, 0)
        fields = ['EIGRL', self.sid, self.v1, self.v2, self.nd,
                  msglvl, self.maxset, self.shfscl, self.norm]
        for (option, value) in izip(self.options, self.values):
            fields += [option + '=' + str(value)]
        return fields
