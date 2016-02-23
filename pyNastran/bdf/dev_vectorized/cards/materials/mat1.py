from __future__ import print_function
from six.moves import zip
from numpy import zeros, where, arange, searchsorted

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)

#from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard
from pyNastran.bdf.dev_vectorized.cards.materials.material import Material


class MAT1(Material):
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
        Material.__init__(self, model)

    def allocate(self, ncards):
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
        self.n = ncards

    def add(self, card, comment=''):
        assert self.n > 0, 'self.n=%s self.i=%s' % (self.n, self.i)
        i = self.i
        print('i=%s' % i)
        mid = integer(card, 1, 'mid')
        if comment:
            self._comments[mid] = comment
        self.material_id[i] = mid
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
        assert self.material_id[i] > 0, self.material_id
        self.i += 1


    def build(self):
        if self.n:
            print('MAT1.materialsA = %s' % self.material_id)
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
            print('MAT1.materialsB = %s' % self.material_id)
            assert self.material_id.min() > 0, 'MAT1.materials = %s' % self.material_id

    def get_D_matrix(self):
        """
        // The isotropic Elasticity matrix D is given by
        //               -                                             -     -                        -
        //               | 1-nu  nu   nu      0         0        0     |     | D0  D1  D1  0   0   0  |
        //    young      |  nu  1-nu  nu      0         0        0     |     | D1  D0  D1  0   0   0  |
        // ------------- |  nu   nu  1-nu     0         0        0     |     | D1  D1  D0  0   0   0  |
        // (1+nu)(1-2nu) |  0    0   0    (1-2nu)/2     0        0     |  =  | 0   0   0   D2  0   0  |
        //               |  0    0   0        0     (1-2nu)/2    0     |     | 0   0   0   0   D2  0  |
        //               |  0    0   0        0         0    (1-2nu)/2 |     | 0   0   0   0   0   D2 |
        //               -                                             -     -                        -
        http://image.diku.dk/svn/OpenTissue/archieve/silcowitz/OpenTissue/dynamics/fem/fem_compute_isotropic_elasticity.h
        """
        poisson2 = 2.0*poisson
        scale = young / ((1.0 + poisson) * (1.0 - poisson2))
        D[0] = (1.0 - poisson) * scale
        D[1] = poisson * scale
        D[2] = young / (2.0  + poisson2)

    def _G_default(self, E, G, nu):
        if G == 0.0 or nu == 0.0:
            pass
        else:
            G = E / 2. / (1 + nu)
        return G

    def get_density_by_material_index(self, i=None):
        if i is None:
            i = slice(None)
        return self.rho[i]

    def get_E_by_material_index(self, i=None):
        if i is None:
            i = slice(None)
        return self.E[i]

    def get_E_by_material_id(self, material_id=None):
        i = self.get_material_index_by_material_id(material_id)
        return self.get_E_by_material_index(i)

    def get_density_by_material_id(self, material_id=None):
        i = self.get_material_index_by_material_id(material_id)
        return self.get_density_by_material_index(i)

    def get_G_by_material_id(self, material_id=None):
        i = self.get_material_index_by_material_id(material_id)
        return self.G[i]

    def write_card(self, f, size=8, material_id=None):
        if self.n:
            if material_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.material_id, material_id)

            assert material_id is None
            #self.model.log.debug"n = %s" % self.n)
            #self.model.log.debug"mids MAT1 %s" % self.material_id)

            Rho = ['' if rhoi == 0.0 else rhoi for rhoi in self.rho[i]]
            A = ['' if ai == 0.0 else ai for ai in self.a[i]]
            TRef = ['' if trefi == 0.0 else trefi for trefi in self.TRef[i]]
            ge = ['' if gei == 0.0 else gei for gei in self.ge[i]]
            St = ['' if st == 0.0 else st for st in self.St[i]]
            Sc = ['' if sc == 0.0 else sc for sc in self.Sc[i]]
            Ss = ['' if ss == 0.0 else ss for ss in self.Ss[i]]

            card_a = ['$MAT1', 'mid', 'E', 'G', 'nu', 'rho', 'a', 'TRef', 'ge']
            card_b = ['$', 'st', 'sc', 'ss', 'mcsid']
            if size == 8:
                fmt_card = print_card_8
            else:
                fmt_card = print_card_16

            f.write(fmt_card(card_a))
            f.write(fmt_card(card_b))

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
                f.write(fmt_card(card))

    def repr_fields(self, material_id):
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
        #self.model.log.debug('G = %s' % G)
        self.E[i] = E
        self.G[i] = G
        self.nu[i] = nu

    def slice_by_index(self, i):
        #assert not isinstance(i, int), type(i)
        #self.model.log.info('i=%s; mids=%s type(i)=%s' % (i, self.material_id, type(i)))
        i = self._validate_slice(i)
        #print('MAT1.slice.i = %s' % i, type(i))

        obj = MAT1(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]

        obj.material_id = self.material_id[i]
        obj.rho = self.rho[i]
        obj.E = self.E[i]
        obj.G = self.G[i]
        obj.nu = self.nu[i]
        obj.a = self.a[i]
        obj.TRef = self.TRef[i]
        obj.ge = self.ge[i]
        obj.St = self.St[i]
        obj.Sc = self.Sc[i]
        obj.Ss = self.Ss[i]
        obj.mcsid = self.mcsid[i]
        return obj
