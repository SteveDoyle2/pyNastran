from numpy import zeros, where, arange, searchsorted, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    string, string_or_blank, blank)


class MAT1(object):
    """
    Defines the material properties for linear isotropic materials.

    +-----+-----+-----+-----+-------+-----+------+------+-----+
    |  1  |  2  | 3   | 4   |   5   |  6  |  7   |  8   |  9  |
    +=====+=====+=====+=====+=======+=====+======+======+=====+
    |MAT1 | MID |  E  |  G  |  NU   | RHO |  A   | TREF | GE  |
    +-----+-----+-----+-----+-------+-----+------+------+-----+
    |     | ST  | SC  | SS  | MCSID |     |      |      |     |
    +-----+-----+-----+-----+-------+-----+------+------+-----+
    """
    type = 'MAT1'
    def __init__(self, model):
        self.model = model
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)
        
    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            self.material_id = zeros(ncards, 'int32')
            self.rho = zeros(ncards, float_fmt)
            self.E = zeros(ncards, float_fmt)
            self.G = zeros(ncards, float_fmt)
            self.nu = zeros(ncards, float_fmt)
            self.a = zeros(ncards, float_fmt)
            self.TRef = zeros(ncards, float_fmt)
            self.ge = zeros(ncards, float_fmt)
            self.St = zeros(ncards, float_fmt)
            self.Sc = zeros(ncards, float_fmt)
            self.Ss = zeros(ncards, float_fmt)
            self.mcsid = zeros(ncards, 'int32')

            for (i, card) in enumerate(cards):
                #if comment:
                #    self._comment = comment
                self.material_id[i] = integer(card, 1, 'mid')
                self.set_E_G_nu(i, card)
                self.rho[i] = double_or_blank(card, 5, 'rho', 0.)
                self.a[i] = double_or_blank(card, 6, 'a', 0.0)
                self.TRef[i] = double_or_blank(card, 7, 'TRef', 0.0)
                self.ge[i] = double_or_blank(card, 8, 'ge', 0.0)
                self.St[i] = double_or_blank(card, 9, 'St', 0.0)
                self.Sc[i] = double_or_blank(card, 10, 'Sc', 0.0)
                self.Ss[i] = double_or_blank(card, 11, 'Ss', 0.0)
                self.mcsid[i] = integer_or_blank(card, 12, 'Mcsid', 0)
                assert len(card) <= 13, 'len(MAT1 card) = %i' % len(card)

            i = self.material_id.argsort()
            self.material_id = self.material_id[i]
            self.E = self.E[i]
            self.G = self.G[i]
            self.rho = self.rho[i]
            self.a = self.a[i]
            self.TRef = self.TRef[i]
            self.ge = self.ge[i]
            self.St = self.St[i]
            self.Sc = self.Sc[i]
            self.Ss = self.Ss[i]
            self.mcsid = self.mcsid[i]
            self._cards = []
            self._comments = []
    
    def _G_default(self, E, G, nu):
        if G == 0.0 or nu == 0.0:
            pass
        else:
            G = E / 2. / (1 + nu)
        return G

    def write_bdf(self, f, size=8, material_ids=None):
        if self.n:
            if material_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.material_id, material_ids)
            
            assert material_ids is None
            #print "n = ", self.n
            #print "mids MAT1", self.material_id
            Rho  = ['' if rhoi  == 0.0 else rhoi  for rhoi  in self.rho[i]]
            A    = ['' if ai    == 0.0 else ai    for ai    in self.a[i]]
            TRef = ['' if trefi == 0.0 else trefi for trefi in self.TRef[i]]
            ge   = ['' if gei   == 0.0 else gei   for gei   in self.ge[i]]
            
            St   = ['' if st    == 0.0 else st    for st    in self.St[i]]
            Sc   = ['' if sc    == 0.0 else sc    for sc    in self.Sc[i]]
            Ss   = ['' if ss    == 0.0 else ss    for ss    in self.Ss[i]]

            card = ['$MAT1', 'mid', 'E', 'G', 'nu', 'rho', 'a', 'TRef', 'ge']
            f.write(print_card(card, size=size))
            card = ['$', 'st', 'sc', 'ss', 'mcsid']
            f.write(print_card(card, size=size))
            for (mid, E, G, nu, rho, a, TRef, ge, st, sc, ss, mcsid) in zip(
                 self.material_id[i], self.E[i], self.G[i], self.nu[i], Rho, A,
                 TRef, ge, St, Sc, Ss, self.mcsid[i]):

                #Gdefault = self.getG_default()
                Gdefault = self._G_default(E, G, nu)
                G = set_blank_if_default(G, Gdefault)
                #rho = set_blank_if_default(rho, 0.)
                #a = set_blank_if_default(a, 0.)
                #TRef = set_blank_if_default(TRef, 0.)
                #ge = set_blank_if_default(ge, 0.)
                #st = set_blank_if_default(st, 0.)
                #sc = set_blank_if_default(sc, 0.)
                #ss = set_blank_if_default(ss, 0.)
                mcsid = set_blank_if_default(mcsid, 0)
                card = ['MAT1', mid, E, G, nu, rho, a, TRef, ge, st, sc, ss, mcsid]
                f.write(print_card(card, size=size))
        
    def reprFields(self, material_id):
        i = where(self.material_id == material_id)[0]
        i = i[0]
        card = ['MAT1', self.material_id[i], self.E[i], self.G[i], self.nu[i],
                        self.rho[i], self.a[i], self.TRef[i], self.ge[i],
                        self.St[i], self.Sc[i], self.Ss[i], self.mcsid[i]]
        return card
        
    def _verify(self, xref=True):
        pass

    def set_E_G_nu(self, i, card):
        r"""
        \f[ G = \frac{E}{2 (1+\nu)} \f]
        """
        E = double_or_blank(card, 2, 'E')
        G = double_or_blank(card, 3, 'G')
        nu = double_or_blank(card, 4, 'nu')

        if G is None and E is None:  # no E,G
            raise RuntimeError('G=%s E=%s cannot both be None' % (G, E))
        elif E is not None and G is not None and nu is not None:
            pass
        elif E is not None and nu is not None:
            G = E / 2. / (1 + nu)
        elif G is not None and nu is not None:
            E = 2 * (1 + nu) * G
        elif G is not None and E is not None:
            nu = E / (2 * G) - 1.
        elif G is None and nu is None:
            G = 0.0
            nu = 0.0
        elif E is None and nu is None:
            E = 0.0
            nu = 0.0
        else:
            msg = 'G=%s E=%s nu=%s' % (G, E, nu)
            raise RuntimeError(msg)
        #print 'G =', G
        self.E[i] = E
        self.G[i] = G
        self.nu[i] = nu