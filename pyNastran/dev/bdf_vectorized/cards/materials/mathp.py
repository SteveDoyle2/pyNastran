from numpy import zeros, where, arange, searchsorted, argsort

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard
from pyNastran.dev.bdf_vectorized.cards.materials.material import Material


class MATHP(Material):
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
    type = 'MATHP'
    def __init__(self, model):
        Material.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            float_fmt = self.model.float_fmt
            self.material_id = zeros(ncards, 'int32')

            self.a10 = zeros(ncards, dtype=float_fmt)
            self.a01 = zeros(ncards, dtype=float_fmt)
            self.d1 = zeros(ncards, dtype=float_fmt)
            self.rho = zeros(ncards, dtype=float_fmt)
            self.av = zeros(ncards, dtype=float_fmt)
            self.tref = zeros(ncards, dtype=float_fmt)
            self.ge = zeros(ncards, dtype=float_fmt)

            self.na = zeros(ncards, dtype='int32')
            self.nd = zeros(ncards, dtype='int32')

            self.a20 = zeros(ncards, dtype=float_fmt)
            self.a11 = zeros(ncards, dtype=float_fmt)
            self.a02 = zeros(ncards, dtype=float_fmt)
            self.d2 = zeros(ncards, dtype=float_fmt)

            self.a30 = zeros(ncards, dtype=float_fmt)
            self.a21 = zeros(ncards, dtype=float_fmt)
            self.a12 = zeros(ncards, dtype=float_fmt)
            self.a03 = zeros(ncards, dtype=float_fmt)
            self.d3 = zeros(ncards, dtype=float_fmt)

            self.a40 = zeros(ncards, dtype=float_fmt)
            self.a31 = zeros(ncards, dtype=float_fmt)
            self.a22 = zeros(ncards, dtype=float_fmt)
            self.a13 = zeros(ncards, dtype=float_fmt)
            self.a04 = zeros(ncards, dtype=float_fmt)
            self.d4 = zeros(ncards, dtype=float_fmt)

            self.a50 = zeros(ncards, dtype=float_fmt)
            self.a41 = zeros(ncards, dtype=float_fmt)
            self.a32 = zeros(ncards, dtype=float_fmt)
            self.a23 = zeros(ncards, dtype=float_fmt)
            self.a14 = zeros(ncards, dtype=float_fmt)
            self.a05 = zeros(ncards, dtype=float_fmt)
            self.d5 = zeros(ncards, dtype=float_fmt)

            self.tab1 = zeros(ncards, dtype='int32')
            self.tab2 = zeros(ncards, dtype='int32')
            self.tab3 = zeros(ncards, dtype='int32')
            self.tab4 = zeros(ncards, dtype='int32')
            self.tabd = zeros(ncards, dtype='int32')
            self.n = ncards

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        if comment:
            self.set_comment(mid, comment)
        self.material_id[i] = integer(card, 1, 'material_id')
        self.model.log.debug('add MATHP.material_id = %s' % self.material_id)
        a10 = double_or_blank(card, 2, 'a10', 0.)
        a01 = double_or_blank(card, 3, 'a01', 0.)
        self.a10[i] = a10
        self.a01[i] = a01
        self.d1[i] = double_or_blank(card, 4, 'd1', (a10 + a01) * 1000)
        self.rho[i] = double_or_blank(card, 5, 'rho', 0.)
        self.av[i] = double_or_blank(card, 6, 'av', 0.)
        self.tref[i] = double_or_blank(card, 7, 'tref', 0.)
        self.ge[i] = double_or_blank(card, 8, 'ge', 0.)

        self.na[i] = integer_or_blank(card, 10, 'na', 1)
        self.nd[i] = integer_or_blank(card, 11, 'nd', 1)

        self.a20[i] = double_or_blank(card, 17, 'a20', 0.)
        self.a11[i] = double_or_blank(card, 18, 'a11', 0.)
        self.a02[i] = double_or_blank(card, 19, 'a02', 0.)
        self.d2[i] = double_or_blank(card, 20, 'd2', 0.)

        self.a30[i] = double_or_blank(card, 25, 'a30', 0.)
        self.a21[i] = double_or_blank(card, 26, 'a21', 0.)
        self.a12[i] = double_or_blank(card, 27, 'a12', 0.)
        self.a03[i] = double_or_blank(card, 28, 'a03', 0.)
        self.d3[i] = double_or_blank(card, 29, 'd3', 0.)

        self.a40[i] = double_or_blank(card, 33, 'a40', 0.)
        self.a31[i] = double_or_blank(card, 34, 'a31', 0.)
        self.a22[i] = double_or_blank(card, 35, 'a22', 0.)
        self.a13[i] = double_or_blank(card, 36, 'a13', 0.)
        self.a04[i] = double_or_blank(card, 37, 'a04', 0.)
        self.d4[i] = double_or_blank(card, 38, 'd4', 0.)

        self.a50[i] = double_or_blank(card, 41, 'a50', 0.)
        self.a41[i] = double_or_blank(card, 42, 'a41', 0.)
        self.a32[i] = double_or_blank(card, 43, 'a32', 0.)
        self.a23[i] = double_or_blank(card, 44, 'a23', 0.)
        self.a14[i] = double_or_blank(card, 45, 'a14', 0.)
        self.a05[i] = double_or_blank(card, 46, 'a05', 0.)
        self.d5[i] = double_or_blank(card, 47, 'd5', 0.)

        self.tab1[i] = integer_or_blank(card, 49, 'tab1', 0)
        self.tab2[i] = integer_or_blank(card, 50, 'tab2', 0)
        self.tab3[i] = integer_or_blank(card, 51, 'tab3', 0)
        self.tab4[i] = integer_or_blank(card, 52, 'tab4', 0)
        self.tabd[i] = integer_or_blank(card, 56, 'tabd', 0)
        assert len(card) <= 57, 'len(MATHP card) = %i\ncard=%s' % (len(card), card)
        self.i += 1


    def build(self):
        if self.n:
            i = argsort(self.material_id)
            self.model.log.debug('build1 MATHP.material_id = %s' % self.material_id)
            self.material_id = self.material_id[i]
            self.model.log.debug('build2 MATHP.material_id = %s' % self.material_id)
            self.a10 = self.a10[i]
            self.a01 = self.a01[i]
            self.d1 = self.d1[i]
            self.rho = self.rho[i]
            self.av = self.av[i]
            self.tref = self.tref[i]
            self.ge = self.ge[i]

            self.na = self.na[i]
            self.nd = self.nd[i]

            self.a20 = self.a20[i]
            self.a11 = self.a11[i]
            self.a02 = self.a02[i]
            self.d2 = self.d2[i]

            self.a30 = self.a30[i]
            self.a21 = self.a21[i]
            self.a12 = self.a12[i]
            self.a03 = self.a03[i]
            self.d3 = self.d3[i]

            self.a40 = self.a40[i]
            self.a31 = self.a31[i]
            self.a22 = self.a22[i]
            self.a13 = self.a13[i]
            self.a04 = self.a04[i]
            self.d4 = self.d4[i]

            self.a50 = self.a50[i]
            self.a41 = self.a41[i]
            self.a32 = self.a32[i]
            self.a23 = self.a23[i]
            self.a14 = self.a14[i]
            self.a05 = self.a05[i]
            self.d5 = self.d5[i]

            self.tab1 = self.tab1[i]
            self.tab2 = self.tab2[i]
            self.tab3 = self.tab3[i]
            self.tab4 = self.tab4[i]
            self.tabd = self.tabd[i]
            self.model.log.debug('MATHP.material_id = %s' % self.material_id)

    #def get_D_matrix(self):
        #"""
        #// The isotropic Elasticity matrix D is given by
        #//               -                                             -     -                        -
        #//               | 1-nu  nu   nu      0         0        0     |     | D0  D1  D1  0   0   0  |
        #//    young      |  nu  1-nu  nu      0         0        0     |     | D1  D0  D1  0   0   0  |
        #// ------------- |  nu   nu  1-nu     0         0        0     |     | D1  D1  D0  0   0   0  |
        #// (1+nu)(1-2nu) |  0    0   0    (1-2nu)/2     0        0     |  =  | 0   0   0   D2  0   0  |
        #//               |  0    0   0        0     (1-2nu)/2    0     |     | 0   0   0   0   D2  0  |
        #//               |  0    0   0        0         0    (1-2nu)/2 |     | 0   0   0   0   0   D2 |
        #//               -                                             -     -                        -
        #http://image.diku.dk/svn/OpenTissue/archieve/silcowitz/OpenTissue/dynamics/fem/fem_compute_isotropic_elasticity.h
        #"""
        #poisson2 = 2.0 * poisson
        #scale = young / ((1.0 + poisson) * (1.0 - poisson2))
        #D[0] = (1.0 - poisson) * scale
        #D[1] = poisson * scale
        #D[2] = young / (2.0  + poisson2)

    def _G_default(self, E, G, nu):
        if G == 0.0 or nu == 0.0:
            pass
        else:
            G = E / 2. / (1 + nu)
        return G

    def get_density_by_index(self, i):
        return self.rho[i]

    def get_density_by_material_id(self, material_id):
        i = self.get_material_index_by_material_id(material_id)
        return self.get_density_by_index(i)

    def write_card(self, bdf_file, size=8, material_id=None):
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
            tref = ['' if trefi == 0.0 else trefi for trefi in self.tref[i]]
            ge = ['' if gei == 0.0 else gei for gei in self.ge[i]]
            St = ['' if st == 0.0 else st for st in self.St[i]]
            Sc = ['' if sc == 0.0 else sc for sc in self.Sc[i]]
            Ss = ['' if ss == 0.0 else ss for ss in self.Ss[i]]

            card = ['$MAT1', 'mid', 'E', 'G', 'nu', 'rho', 'a', 'tref', 'ge']
            bdf_file.write(print_card_8(card))
            card = ['$', 'st', 'sc', 'ss', 'mcsid']
            bdf_file.write(print_card_8(card))
            for(mid, E, G, nu, rho, a, tref, ge, st, sc, ss, mcsid) in zip(
                self.material_id[i], self.E[i], self.G[i], self.nu[i], Rho, A,
                tref, ge, St, Sc, Ss, self.mcsid[i]):

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
                bdf_file.write(print_card_8(card))

    def repr_fields(self, material_id):
        i = where(self.material_id == material_id)[0]
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
        i = self._validate_slice(i)
        #self.model.log.debug('i = %s' % i)
        obj = MATHP(self.model)
        n = len(i)
        obj.n = n
        obj.i = n
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.material_id = self.material_id[i]

        obj.a10 = self.a10[i]
        obj.a01 = self.a01[i]
        obj.d1 = self.d1[i]
        obj.rho = self.rho[i]
        obj.av = self.av[i]
        obj.tref = self.tref[i]
        obj.ge = self.ge[i]

        obj.na = self.na[i]
        obj.nd = self.nd[i]

        obj.a20 = self.a20[i]
        obj.a11 = self.a11[i]
        obj.a02 = self.a02[i]
        obj.d2 = self.d2[i]

        obj.a30 = self.a30[i]
        obj.a21 = self.a21[i]
        obj.a12 = self.a12[i]
        obj.a03 = self.a03[i]
        obj.d3 = self.d3[i]

        obj.a40 = self.a40[i]
        obj.a31 = self.a31[i]
        obj.a22 = self.a22[i]
        obj.a13 = self.a13[i]
        obj.a04 = self.a04[i]
        obj.d4 = self.d4[i]

        obj.a50 = self.a50[i]
        obj.a41 = self.a41[i]
        obj.a32 = self.a32[i]
        obj.a23 = self.a23[i]
        obj.a14 = self.a14[i]
        obj.a05 = self.a05[i]
        obj.d5 = self.d5[i]

        obj.tab1 = self.tab1[i]
        obj.tab2 = self.tab2[i]
        obj.tab3 = self.tab3[i]
        obj.tab4 = self.tab4[i]
        obj.tabd = self.tabd[i]
        self.model.log.debug("obj = %s" % obj.material_id)
        return obj
