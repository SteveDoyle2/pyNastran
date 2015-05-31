# pylint: disable=C0103,R0902,R0904,R0914
"""
All dynamic control cards are defined in this file.  This includes:

 * FREQ
 * FREQ1
 * FREQ2 (not implemented)
 * FREQ3
 * FREQ4
 * FREQ5 (not implemented)
 * NLPCI
 * NLPARM
 * TSTEP
 * TSTEPNL

All cards are BaseCard objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import log, exp, ceil
from six.moves import zip, range

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string_or_blank, blank, fields)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class FREQ(BaseCard):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems.

    +-----+-----+-----+-----+------+-----+-----+-----+-----+
    |  1  |  2  |  3  |  4  |  5   |  6  |  7  |  8  |  9  |
    +=====+=====+=====+=====+======+=====+=====+=====+=====+
    |FREQ | SID | F1  | F2  | etc. |     |     |     |     |
    +-----+-----+-----+-----+------+-----+-----+-----+-----+
    """
    type = 'FREQ'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        self.freqs = fields(double, card, 'freq', i=2, j=len(card))
        self.cleanFreqs()

    def cleanFreqs(self):
        self.freqs = list(set(self.freqs))
        self.freqs.sort()

    def getFreqs(self):
        return self.freqs

    def add_frequencies(self, freqs):
        """
        Combines the frequencies from 1 FREQx object with another.
        All FREQi entries with the same frequency set identification numbers
        will be used. Duplicate frequencies will be ignored.

        :param self:  the object pointer
        :param freqs: the frequencies for a FREQx object
        """
        #print("self.freqs = ",self.freqs)
        #print("freqs = ",freqs)
        self.freqs += freqs
        self.cleanFreqs()

    def add_frequency_object(self, freq):
        """
        :param self: the object pointer
        :param freq: a FREQx object

        .. seealso:: :func:`addFrequencies`
        """
        self.add_frequencies(freq.freqs)

    def raw_fields(self):
        list_fields = ['FREQ', self.sid] + self.freqs
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class FREQ1(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems by specification of a starting frequency, frequency
    increment, and the number of increments desired.

    +------+-----+-----+-----+-----+-----+-----+-----+-----+
    |  1   |  2  | 3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +======+=====+=====+=====+=====+=====+=====+=====+=====+
    |FREQ1 | SID |  F1 | DF  | NDF |     |     |     |     |
    +------+-----+-----+-----+-----+-----+-----+-----+-----+

    .. note:: this card rewrites as a FREQ card
    """
    type = 'FREQ1'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        f1 = double_or_blank(card, 2, 'f1', 0.0)
        df = double(card, 3, 'df')
        ndf = integer_or_blank(card, 4, 'ndf', 1)
        assert len(card) <= 5, 'len(FREQ card) = %i' % len(card)

        self.freqs = []
        for i in range(ndf):
            self.freqs.append(f1 + i * df)
        self.cleanFreqs()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class FREQ2(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems by specification of a starting frequency, final
    frequency, and the number of logarithmic increments desired.

    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |   1   |  2  | 3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+=====+
    | FREQ2 | SID |  F1 | F2  | NDF |     |     |     |     |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+

    .. note:: this card rewrites as a FREQ card
    """
    type = 'FREQ2'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        f1 = double(card, 2, 'f1')  # default=0.0 ?
        f2 = double(card, 3, 'f2')
        nf = integer_or_blank(card, 4, 'nf', 1)
        assert len(card) <= 5, 'len(FREQ2 card) = %i' % len(card)

        d = 1. / nf * log(f2 / f1)
        self.freqs = []
        for i in range(nf):
            self.freqs.append(f1 * exp(i * d))  # 0 based index
        self.cleanFreqs()


class FREQ3(FREQ):
    type = 'FREQ3'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        raise NotImplementedError()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class FREQ4(FREQ):
    """
    Defines a set of frequencies used in the solution of modal frequency
    response problems by specifying the amount of 'spread' around each natural
    frequency and the number of equally spaced excitation frequencies within
    the spread.

    +------+-----+-----+-----+------+-----+-----+-----+-----+
    |  1   |  2  | 3   |  4  |  5   |  6  |  7  |  8  |  9  |
    +======+=====+=====+=====+======+=====+=====+=====+=====+
    |FREQ4 | SID |  F1 | F2  | FSPD | NFM |     |     |     |
    +------+-----+-----+-----+------+-----+-----+-----+-----+

    .. note:: this card rewrites as a FREQ card
    .. todo:: not done...
    """
    type = 'FREQ4'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        self.f1 = double_or_blank(card, 2, 'f1', 0.0)
        self.f2 = double_or_blank(card, 3, 'f2', 1.e20)
        self.fspd = double_or_blank(card, 4, 'fspd', 0.1)
        self.nfm = integer_or_blank(card, 5, 'nfm', 3)
        assert len(card) <= 6, 'len(FREQ card) = %i' % len(card)

    def raw_fields(self):
        list_fields = ['FREQ4', self.sid, self.f1, self.f2, self.fspd,
                       self.nfm]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class FREQ5(FREQ):
    type = 'FREQ5'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        raise NotImplementedError()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class NLPARM(BaseCard):
    """
    Defines a set of parameters for nonlinear static analysis iteration
    strategy.

    +--------+--------+------+------+---------+-------+---------+---------+--------+
    |    1   |  2     |  3   |  4   |    5    |   6   |   7     |  8      |   9    |
    +========+========+======+======+=========+=======+=========+=========+========+
    | NLPARM |   ID   | NINC |  DT  | KMETHOD | KSTEP | MAXITER |  CONV   | INTOUT |
    +--------+--------+------+------+---------+-------+---------+---------+--------+
    |        |  ESPU  | EPSP | EPSW | MAXDIV  | MAXQN | MAXLS   | FSTRESS | LSTOL  |
    +--------+--------+------+------+---------+-------+---------+---------+--------+
    |        | MAXBIS |      |      |         | MAXR  |         | RTOLB   | CONV   |
    +--------+--------+------+------+---------+-------+---------+---------+--------+
    """
    type = 'NLPARM'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.nlparm_id = integer(card, 1, 'nlparm_id')
            self.ninc = integer_or_blank(card, 2, 'ninc', 10)
            self.dt = double_or_blank(card, 3, 'dt', 0.0)
            self.kMethod = string_or_blank(card, 4, 'kMethod', 'AUTO')
            self.kStep = integer_or_blank(card, 5, 'kStep', 5)
            self.maxIter = integer_or_blank(card, 6, 'maxIter', 25)
            self.conv = string_or_blank(card, 7, 'conv', 'PW')
            self.intOut = string_or_blank(card, 8, 'intOut', 'NO')

            # line 2
            self.epsU = double_or_blank(card, 9, 'epsU', 0.01)
            self.epsP = double_or_blank(card, 10, 'epsP', 0.01)
            self.epsW = double_or_blank(card, 11, 'epsW', 0.01)
            self.maxDiv = integer_or_blank(card, 12, 'maxDiv', 3)

            if self.kMethod == 'PFNT':
                self.maxQn = integer_or_blank(card, 13, 'maxQn', 0)
            else:
                self.maxQn = integer_or_blank(card, 13, 'maxQn', self.maxIter)

            self.maxLs = integer_or_blank(card, 14, 'maxLs', 4)
            self.fStress = double_or_blank(card, 15, 'fStress', 0.2)
            self.lsTol = double_or_blank(card, 16, 'lsTol', 0.5)

            # line 3
            self.maxBisect = integer_or_blank(card, 17, '', 5)
            self.maxR = double_or_blank(card, 21, 'maxR', 20.)
            self.rTolB = double_or_blank(card, 23, 'rTolB', 20.)
            assert len(card) <= 24, 'len(NLPARM card) = %i' % len(card)
        else:
            (nlparm_id, ninc, dt, kMethod, kStep, maxIter, conv, intOut, epsU, epsP,
             epsW, maxDiv, maxQn, maxLs, fStress, lsTol, maxBisect, maxR,
             rTolB) = data
            self.nlparm_id = nlparm_id
            self.ninc = ninc
            self.dt = dt
            self.kMethod = kMethod
            self.kStep = kStep
            self.maxIter = maxIter
            self.conv = conv
            self.intOut = intOut

            # line 2
            self.epsU = epsU
            self.epsP = epsP
            self.epsW = epsW
            self.maxDiv = maxDiv
            self.maxQn = maxQn
            self.maxLs = maxLs
            self.fStress = fStress
            self.lsTol = lsTol

            # line 3
            self.maxBisect = maxBisect
            self.maxR = maxR
            self.rTolB = rTolB

    def raw_fields(self):
        list_fields = ['NLPARM', self.nlparm_id, self.ninc, self.dt, self.kMethod,
                       self.kStep, self.maxIter, self.conv, self.intOut, self.epsU,
                       self.epsP, self.epsW, self.maxDiv, self.maxQn, self.maxLs,
                       self.fStress, self.lsTol, self.maxBisect, None, None, None,
                       self.maxR, None, self.rTolB]
        return list_fields

    def repr_fields(self):
        ninc = set_blank_if_default(self.ninc, 10)
        dt = set_blank_if_default(self.dt, 0.0)
        kMethod = set_blank_if_default(self.kMethod, 'AUTO')
        kStep = set_blank_if_default(self.kStep, 5)
        maxIter = set_blank_if_default(self.maxIter, 25)
        conv = set_blank_if_default(self.conv, 'PW')
        intOut = set_blank_if_default(self.intOut, 'NO')
        epsU = set_blank_if_default(self.epsU, 0.01)
        epsP = set_blank_if_default(self.epsP, 0.01)
        epsW = set_blank_if_default(self.epsW, 0.01)
        maxDiv = set_blank_if_default(self.maxDiv, 3)
        maxQn = set_blank_if_default(self.maxQn, self.maxIter)
        maxLs = set_blank_if_default(self.maxLs, 4)
        fStress = set_blank_if_default(self.fStress, 0.2)
        lsTol = set_blank_if_default(self.lsTol, 0.5)
        maxBisect = set_blank_if_default(self.maxBisect, 5)
        maxR = set_blank_if_default(self.maxR, 20.)
        rTolB = set_blank_if_default(self.rTolB, 20.)

        list_fields = ['NLPARM', self.nlparm_id, ninc, dt, kMethod, kStep, maxIter,
                       conv, intOut, epsU, epsP, epsW, maxDiv, maxQn, maxLs,
                       fStress, lsTol, maxBisect, None, None, None, maxR, None,
                       rTolB]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card) # having trouble with double precision...
        return self.comment() + print_card_16(card)


class NLPCI(BaseCard):
    type = 'NLPCI'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.nlpci_id = integer(card, 1, 'nlpci_id')
        self.Type = string_or_blank(card, 2, 'Type', 'CRIS')
        self.minalr = double_or_blank(card, 3, 'minalr', 0.25)
        self.maxalr = double_or_blank(card, 4, 'maxalr', 4.0)
        self.scale = double_or_blank(card, 5, 'scale', 0.0)
        blank(card, 6, 'blank')
        self.desiter = integer_or_blank(card, 7, 'desiter', 12)
        self.mxinc = integer_or_blank(card, 8, 'mxinc', 20)

    def raw_fields(self):
        list_fields = ['NLPCI', self.nlpci_id, self.Type, self.minalr,
                       self.maxalr, self.scale, None, self.desiter, self.mxinc]
        return list_fields

    def repr_fields(self):
        #minalr = set_blank_if_default(self.minalr, 0.25)
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class TSTEP(BaseCard):
    """
    Transient Time Step
    Defines time step intervals at which a solution will be generated and
    output in transient analysis.

    +-------+------+-----+-----+-----+-----+-----+-----+-----+
    |   1   |   2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +-------+------+-----+-----+-----+-----+-----+-----+-----+
    | TSTEP |  N1  | DT1 | NO1 |     |     |     |     |     |
    +-------+------+-----+-----+-----+-----+-----+-----+-----+
    |       |  N2  | DT2 | NO2 |     |     |     |     |     |
    +-------+------+-----+-----+-----+-----+-----+-----+-----+
    |       | etc. |     |     |     |     |     |     |     |
    +-------+------+-----+-----+-----+-----+-----+-----+-----+
    """
    type = 'TSTEP'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.sid = integer(card, 1, 'sid')
        self.N = []
        self.DT = []
        self.NO = []

        nrows = int(ceil((len(card) - 1.) / 8.))
        for i in range(nrows):
            n = 8 * i + 1
            N = integer_or_blank(card, n + 1, 'N' + str(i), 1)
            dt = double_or_blank(card, n + 2, 'dt' + str(i), 0.)
            no = integer_or_blank(card, n + 3, 'NO' + str(i), 1)
            self.N.append(N)
            self.DT.append(dt)
            self.NO.append(no)

    def raw_fields(self):
        list_fields = ['TSTEP', self.sid]
        for (N, dt, no) in zip(self.N, self.DT, self.NO):
            list_fields += [N, dt, no, None, None, None, None, None]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class TSTEPNL(BaseCard):
    """
    Defines parametric controls and data for nonlinear transient structural or
    heat transfer analysis. TSTEPNL is intended for SOLs 129, 159, and 600.
    Parameters for Nonlinear Transient Analysis.

    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    |    1    |   2    |   3    |   4   |   5    |   6    |   7   |    8    |  9   |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    | TSTEPNL |   ID   |  NDT   |  DT   |   NO   | METHOD | KSTEP | MAXITER | CONV |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    |         |  ESPU  |  EPSP  |  EPSW | MAXDIV | MAXQN  | MAXLS | FSTRESS |      |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    |         | MAXBIS | ADJUST | MSTEP |   RB   | MAXR   | UTOL  | RTOLB   |      |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    """
    type = 'TSTEPNL'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.ndt = integer(card, 2, 'ndt')
            assert self.ndt >= 3
            self.dt = double(card, 3, 'dt')
            assert self.dt > 0.
            self.no = integer_or_blank(card, 4, 'no', 1)

            #: .. note:: not listed in all QRGs
            self.method = string_or_blank(card, 5, 'method', 'ADAPT')
            if self.method == 'ADAPT':
                self.kStep = integer_or_blank(card, 6, 'kStep', 2)
            elif self.method == 'ITER':
                self.kStep = integer_or_blank(card, 6, 'kStep', 10)
            elif self.method in ['AUTO', 'TSTEP']:
                self.kStep = None
                #self.kStep = blank(card, 6, 'kStep') #: .. todo:: not blank
            else:
                msg = 'invalid TSTEPNL Method.  method=|%s|' % (self.method)
                raise RuntimeError(msg)
            self.maxIter = integer_or_blank(card, 7, 'maxIter', 10)
            self.conv = string_or_blank(card, 8, 'conv', 'PW')

            # line 2
            self.epsU = double_or_blank(card, 9, 'epsU', 1.E-2)
            self.epsP = double_or_blank(card, 10, 'epsP', 1.E-3)
            self.epsW = double_or_blank(card, 11, 'epsW', 1.E-6)
            self.maxDiv = integer_or_blank(card, 12, 'maxDiv', 2)
            self.maxQn = integer_or_blank(card, 13, 'maxQn', 10)
            self.MaxLs = integer_or_blank(card, 14, 'MaxLs', 2)
            self.fStress = double_or_blank(card, 15, 'fStress', 0.2)

            # line 3
            self.maxBisect = integer_or_blank(card, 17, 'maxBisect', 5)
            self.adjust = integer_or_blank(card, 18, 'adjust', 5)
            self.mStep = integer_or_blank(card, 19, 'mStep')
            self.rb = double_or_blank(card, 20, 'rb', 0.6)
            self.maxR = double_or_blank(card, 21, 'maxR', 32.)
            self.uTol = double_or_blank(card, 22, 'uTol', 0.1)
            self.rTolB = double_or_blank(card, 23, 'rTolB', 20.)

            # not listed in all QRGs
            self.minIter = integer_or_blank(card, 24, 'minIter')
            assert len(card) <= 25, 'len(TSTEPNL card) = %i' % len(card)
        else:
            (sid, ndt, dt, no, method, kStep, maxIter, conv, epsU, epsP, epsW,
             maxDiv, maxQn, maxLs, fStress, maxBisect,
             adjust, mStep, rb, maxR, uTol, rTolB) = data
            self.sid = sid
            self.ndt = ndt
            self.dt = dt
            self.no = no
            self.method = method
            self.kStep = kStep
            self.maxIter = maxIter
            self.conv = conv

            # line 2
            self.epsU = epsU
            self.epsP = epsP
            self.epsW = epsW
            self.maxDiv = maxDiv
            self.maxQn = maxQn
            self.MaxLs = maxLs
            self.fStress = fStress

            # line 3
            self.maxBisect = maxBisect
            self.adjust = adjust
            self.mStep = mStep
            self.rb = rb
            self.maxR = maxR
            self.uTol = uTol
            self.rTolB = rTolB
            self.minIter = None  # not listed in DMAP 2005

    def raw_fields(self):
        list_fields = ['TSTEPNL', self.sid, self.ndt, self.dt, self.no,
                       self.method, self.kStep, self.maxIter, self.conv, self.epsU,
                       self.epsP, self.epsW, self.maxDiv, self.maxQn, self.MaxLs,
                       self.fStress, None, self.maxBisect, self.adjust, self.mStep,
                       self.rb, self.maxR, self.uTol, self.rTolB, self.minIter]
        return list_fields

    def repr_fields(self):
        #no = set_blank_if_default(self.no,1)
        no = self.no
        method = set_blank_if_default(self.method, 'ADAPT')

        kStep = self.kStep
        #if self.method=='ADAPT':
            #kStep = set_blank_if_default(self.kStep, 2)
        #elif self.method=='ITER':
            #kStep = set_blank_if_default(self.kStep, 10)
        #else:
            #msg = 'invalid TSTEPNL Method.  method=|%s|' %(self.method)
            #raise RuntimeError(msg)

        #maxIter = set_blank_if_default(self.maxIter, 10)
        conv = set_blank_if_default(self.conv, 'PW')

        epsU = set_blank_if_default(self.epsU, 1e-2)
        epsP = set_blank_if_default(self.epsP, 1e-3)
        epsW = set_blank_if_default(self.epsW, 1e-6)
        maxDiv = set_blank_if_default(self.maxDiv, 2)
        maxQn = set_blank_if_default(self.maxQn, 10)
        MaxLs = set_blank_if_default(self.MaxLs, 2)
        fStress = set_blank_if_default(self.fStress, 0.2)

        maxBisect = set_blank_if_default(self.maxBisect, 5)
        adjust = set_blank_if_default(self.adjust, 5)
        rb = set_blank_if_default(self.rb, 0.6)
        maxR = set_blank_if_default(self.maxR, 32.)
        uTol = set_blank_if_default(self.uTol, 0.1)
        rTolB = set_blank_if_default(self.rTolB, 20.)

        list_fields = ['TSTEPNL', self.sid, self.ndt, self.dt, no, method,
                       kStep, self.maxIter, conv, epsU, epsP, epsW, maxDiv, maxQn,
                       MaxLs, fStress, None, maxBisect, adjust, self.mStep, rb,
                       maxR, uTol, rTolB, self.minIter]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
