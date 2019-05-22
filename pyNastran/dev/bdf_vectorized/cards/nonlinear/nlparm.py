from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank, string_or_blank)


class NLPARM:
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

    def __init__(self):
        pass

    def add_card(self, card=None, comment=''):
        if comment:
            self.comment = comment

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
        assert len(card) <= 24, 'len(NLPARM card) = %i\ncard=%s' % (len(card), card)

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
        max_bisect = set_blank_if_default(self.maxBisect, 5)
        maxR = set_blank_if_default(self.maxR, 20.)
        rTolB = set_blank_if_default(self.rTolB, 20.)

        list_fields = ['NLPARM', self.nlparm_id, ninc, dt, kMethod, kStep, maxIter,
                       conv, intOut, epsU, epsP, epsW, maxDiv, maxQn, maxLs,
                       fStress, lsTol, max_bisect, None, None, None, maxR, None,
                       rTolB]
        return list_fields

    def write_card(self, bdf_file, size=8):
        card = self.raw_fields()
        if size == 8:
            bdf_file.write(print_card_8(card))
        else:
            bdf_file.write(print_card_16(card))
