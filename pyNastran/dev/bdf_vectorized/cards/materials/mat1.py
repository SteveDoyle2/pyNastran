import numpy as np

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)

#from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard
from pyNastran.dev.bdf_vectorized.cards.materials.material import Material


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
        self.material_id = []
        self.rho = []
        self.E = []
        self.G = []
        self.nu = []
        self.a = []
        self.tref = []
        self.ge = []
        self.St = []
        self.Sc = []
        self.Ss = []
        self.mcsid = []

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            float_fmt = self.model.float_fmt
            self.material_id = np.zeros(ncards, 'int32')
            self.rho = np.zeros(ncards, float_fmt)
            self.E = np.zeros(ncards, float_fmt)
            self.G = np.zeros(ncards, float_fmt)
            self.nu = np.zeros(ncards, float_fmt)
            self.a = np.zeros(ncards, float_fmt)
            self.tref = np.zeros(ncards, float_fmt)
            self.ge = np.zeros(ncards, float_fmt)
            self.St = np.zeros(ncards, float_fmt)
            self.Sc = np.zeros(ncards, float_fmt)
            self.Ss = np.zeros(ncards, float_fmt)
            self.mcsid = np.zeros(ncards, 'int32')
            self.n = ncards

    def add_card(self, card, comment=''):
        assert self.n > 0, 'self.n=%s self.i=%s' % (self.n, self.i)
        i = self.i
        #self.model.log.debug('i=%s' % i)
        mid = integer(card, 1, 'mid')
        if comment:
            self.set_comment(mid, comment)
        self.material_id[i] = mid
        self.set_E_G_nu(i, card)
        self.rho[i] = double_or_blank(card, 5, 'rho', 0.)
        self.a[i] = double_or_blank(card, 6, 'a', 0.0)
        self.tref[i] = double_or_blank(card, 7, 'tref', 0.0)
        self.ge[i] = double_or_blank(card, 8, 'ge', 0.0)
        self.St[i] = double_or_blank(card, 9, 'St', 0.0)
        self.Sc[i] = double_or_blank(card, 10, 'Sc', 0.0)
        self.Ss[i] = double_or_blank(card, 11, 'Ss', 0.0)
        self.mcsid[i] = integer_or_blank(card, 12, 'Mcsid', 0)
        #if mid == 5:
        #print(card)
        #print('i=%-2s; E=%s; G=%s; nu=%s\n' % (self.i, self.E[self.i], self.G[self.i], self.nu[self.i]))
        #print(self.print_card(self.i) + '\n')

        assert len(card) <= 13, 'len(MAT1 card) = %i\ncard=%s' % (len(card), card)
        assert self.material_id[i] > 0, self.material_id
        self.i += 1


    def build(self):
        if self.n:
            self.model.log.debug('MAT1.materialsA = %s' % self.material_id)
            #print('G =', self.G)
            #print('nu =', self.nu)
            i = self.material_id.argsort()
            #if not np.array_equal(np.arange(i.size), i):
            self.material_id = self.material_id[i]
            self.E = self.E[i]
            self.G = self.G[i]
            self.nu = self.nu[i]
            self.rho = self.rho[i]
            self.a = self.a[i]
            self.tref = self.tref[i]
            self.ge = self.ge[i]
            self.St = self.St[i]
            self.Sc = self.Sc[i]
            self.Ss = self.Ss[i]
            self.mcsid = self.mcsid[i]
            self.model.log.debug('MAT1.materialsB = %s' % self.material_id)
            assert self.material_id.min() > 0, 'MAT1.materials = %s' % self.material_id
            #print('G =', self.G)
            #print('nu =', self.nu)

    def update(self, maps):
        """
        maps = {
            'material' : mid_map,
        }
        """
        if self.n:
            mid_map = maps['material']
            for i, mid in enumerate(self.material_id):
                self.material_id[i] = mid_map[mid]

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
        D = np.array((6, 6), dtype='float64')
        D[0] = (1.0 - poisson) * scale
        D[1] = poisson * scale
        D[2] = young / (2.0  + poisson2)
        return D

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

    def write_card(self, bdf_file, material_id=None, size=8, is_double=False):
        if self.n:
            if material_id is None:
                i = np.arange(self.n)
            else:
                if isinstance(material_id, int):
                    material_id = [material_id]
                i = np.searchsorted(self.material_id, material_id)

            #print('imat1 = ', i)
            #assert material_id is None, 'i=%i' % i
            #self.model.log.debug"n = %s" % self.n)
            #self.model.log.debug"mids MAT1 %s" % self.material_id)

            Rho = ['' if rhoi == 0.0 else rhoi for rhoi in self.rho[i]]
            A = ['' if ai == 0.0 else ai for ai in self.a[i]]
            tref = ['' if trefi == 0.0 else trefi for trefi in self.tref[i]]
            ge = ['' if gei == 0.0 else gei for gei in self.ge[i]]
            St = ['' if st == 0.0 else st for st in self.St[i]]
            Sc = ['' if sc == 0.0 else sc for sc in self.Sc[i]]
            Ss = ['' if ss == 0.0 else ss for ss in self.Ss[i]]

            if size == 8:
                fmt_card = print_card_8
                bdf_file.write('$MAT1        mid       E       G      nu     rho       a    tref      ge\n')
                bdf_file.write('$             st      sc      ss   mcsid\n')
            else:
                bdf_file.write('$MAT1*               mid               E               G              nu\n')
                bdf_file.write('$*                   rho               a            tref              ge\n')
                bdf_file.write('$*                   st              sc              ss           mcsid\n')
                if is_double:
                    fmt_card = print_card_double
                else:
                    fmt_card = print_card_16


            for (mid, E, G, nu, rho, a, tref, ge, st, sc, ss, mcsid) in zip(
                    self.material_id[i], self.E[i], self.G[i], self.nu[i], Rho, A,
                    tref, ge, St, Sc, Ss, self.mcsid[i]):
                if mid in self._comments:
                    bdf_file.write(self._comments[mid])

                #Gdefault = self.getG_default()
                Gdefault = self._G_default(E, G, nu)
                G = set_blank_if_default(G, Gdefault)
                #rho = set_blank_if_default(rho, 0.)
                #a = set_blank_if_default(a, 0.)
                #tref = set_blank_if_default(tref, 0.)
                #ge = set_blank_if_default(ge, 0.)
                #st = set_blank_if_default(st, 0.)
                #sc = set_blank_if_default(sc, 0.)
                #ss = set_blank_if_default(ss, 0.)
                mcsid = set_blank_if_default(mcsid, 0)
                card = ['MAT1', mid, E, G, nu, rho, a, tref, ge, st, sc, ss, mcsid]
                bdf_file.write(fmt_card(card))

    def repr_fields(self, material_id):
        i = np.where(self.material_id == material_id)[0]
        i = i[0]
        card = ['MAT1', self.material_id[i], self.E[i], self.G[i], self.nu[i],
                self.rho[i], self.a[i], self.tref[i], self.ge[i],
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
        obj.tref = self.tref[i]
        obj.ge = self.ge[i]
        obj.St = self.St[i]
        obj.Sc = self.Sc[i]
        obj.Ss = self.Ss[i]
        obj.mcsid = self.mcsid[i]
        return obj
