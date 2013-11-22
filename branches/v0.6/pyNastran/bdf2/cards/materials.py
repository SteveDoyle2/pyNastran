from numpy import zeros, where, arange, searchsorted, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    string, string_or_blank, blank)


class Materials(object):
    def __init__(self, model):
        self.model = model

        self._mat1 = []
        self._mat1_comments = []
        self.mat1 = MAT1(model)
    
    def add_mat1(self, card, comment):
        self._mat1.append(card)
        self._mat1_comments.append(comment)

    def build(self):
        self.mat1.build(self._mat1)
        #self.mat1s.build(self._mat1s)
        #self.mat2.build(self._mat2)
        #self.mat8.build(self._mat8)

    def get_stats(self):
        msg = []
        types = self._get_types()
        for mat in types:
            nmat = len(mat.mid)
            if nmat:
                msg.append('  %-8s: %i' % (mat.type, nmat))
        return msg

    def __getitem__(self, mid):
        assert isinstance(mid, int), 'mid=%r' % mid
        types = self._get_types()
        for mat in types:
            if mid in mat.mid:
                return mat
        raise RuntimeError('mid=%r does not exist' % mid)

    def __len__(self):
        types = self._get_types()
        n = 0
        for mat in types:
            n += len(mat.mid)
        return n

    def __iter__(self):
        types = self._get_types()
        for mat in types:
            for mid in mat.mid:
                yield mid

    def _get_types(self):
        return [self.mat1]

    def _verify(self, xref=True):
        self.mat1._verify(xref)

    def write_bdf(self, f, size=8):
        f.write('$MATERIALS\n')
        self.mat1.write_bdf(f, size=size)
        #self.mat1s.write_bdf(f, size=size)
        #self.mat2.write_bdf(f, size=size)
        #self.mat8.write_bdf(f, size=size)

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

    def build(self, cards):
        ncards = len(cards)
        assert ncards > 0
        
        self.mid = zeros(ncards, 'int32')
        self.rho = zeros(ncards, 'float32')
        self.E = zeros(ncards, 'float32')
        self.G = zeros(ncards, 'float32')
        self.nu = zeros(ncards, 'float32')
        self.a = zeros(ncards, 'float32')
        self.TRef = zeros(ncards, 'float32')
        self.ge = zeros(ncards, 'float32')
        self.St = zeros(ncards, 'float32')
        self.Sc = zeros(ncards, 'float32')
        self.Ss = zeros(ncards, 'float32')
        self.mcsid = zeros(ncards, 'int32')
        
        for (i, card) in enumerate(cards):
            #if comment:
            #    self._comment = comment
            self.mid[i] = integer(card, 1, 'mid')
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

        i = self.mid.argsort()
        self.mid = self.mid[i]
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
    
    def _G_default(self, E, G, nu):
        if G == 0.0 or self.nu == 0.0:
            pass
        else:
            G = E / 2. / (1 + nu)
        return G

    def write_bdf(self, f, size=8, mids=None):
        assert mids is None
        Rho  = ['' if rhoi  == 0.0 else rhoi for rhoi  in self.rho]
        A    = ['' if ai    == 0.0 else ai   for ai    in self.a]
        TRef = ['' if trefi == 0.0 else tefi for trefi in self.TRef]
        ge   = ['' if gei   == 0.0 else gei  for gei   in self.ge]
        
        card = ['$MAT1', 'mid', 'E', 'G', 'nu', 'rho', 'a', 'TRef', 'ge']
        f.write(print_card(card, size=size))
        card = ['$', 'st', 'sc', 'ss', 'mcsid']
        f.write(print_card(card, size=size))
        for (mid, E, G, nu, rho, a, TRef, ge, st, sc, ss, mcsid) in zip(
             self.mid, self.E, self.G, self.nu, Rho, A,
             TRef, ge, self.St, self.Sc, self.Ss, self.mcsid):

            #Gdefault = self.getG_default()
            #G = set_blank_if_default(self.g, Gdefault)
            G = self._G_default(E, G, nu)
            #rho = set_blank_if_default(rho, 0.)
            #a = set_blank_if_default(a, 0.)
            #TRef = set_blank_if_default(TRef, 0.)
            #ge = set_blank_if_default(ge, 0.)
            st = set_blank_if_default(st, 0.)
            sc = set_blank_if_default(sc, 0.)
            ss = set_blank_if_default(ss, 0.)
            mcsid = set_blank_if_default(mcsid, 0)
            card = ['MAT1', mid, E, G, nu, rho, a, TRef, ge, st, sc, ss, mcsid]
            f.write(print_card(card, size=size))
        
    def reprFields(self, mid):
        i = where(self.mid == mid)[0]
        assert len(i) == 1, i
        i = i[0]
        card = ['MAT1', self.mid[i], self.E[i], self.G[i], self.nu[i],
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
        self.E[i] = E
        self.G[i] = G
        self.nu[i] = nu
    