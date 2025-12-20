from numpy import zeros, where, arange, searchsorted

from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank, string_or_blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.materials.material import Material


class MAT8(Material):
    """
    Defines the material property for an orthotropic material for isoparametric
    shell elements.

    +-----+-----+-----+------+------+-----+-----+-----+-----+
    |  1  |  2  |  3  |  4   |  5   |  6  |  7  |  8  |  9  |
    +=====+=====+=====+======+======+=====+=====+=====+=====+
    |MAT8 | MID | E1  |  E2  | NU12 | G12 | G1Z | G2Z | RHO |
    +-----+-----+-----+------+------+-----+-----+-----+-----+
    |     | A1  |  A2 | TREF |  Xt  |  Xc |  Yt |  Yc |  S  |
    +-----+-----+-----+------+------+-----+-----+-----+-----+
    |     | GE1 | F12 | STRN |      |     |     |     |     |
    +-----+-----+-----+------+------+-----+-----+-----+-----+
    """
    type = 'MAT8'
    def __init__(self, model):
        Material.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.material_id = zeros(ncards, dtype='int32')
            self.e11 = zeros(ncards, dtype=float_fmt)
            self.e22 = zeros(ncards, dtype=float_fmt)
            self.nu12 = zeros(ncards, dtype=float_fmt)

            self.g12 = zeros(ncards, dtype=float_fmt)
            self.g1z = zeros(ncards, dtype=float_fmt)
            self.g2z = zeros(ncards, dtype=float_fmt)
            self.rho = zeros(ncards, dtype=float_fmt)
            self.a1 = zeros(ncards, dtype=float_fmt)
            self.a2 = zeros(ncards, dtype=float_fmt)
            self.tref = zeros(ncards, dtype=float_fmt)
            self.Xt = zeros(ncards, dtype=float_fmt)
            self.Xc = zeros(ncards, dtype=float_fmt)
            self.Yt = zeros(ncards, dtype=float_fmt)
            self.Yc = zeros(ncards, dtype=float_fmt)
            self.S = zeros(ncards, dtype=float_fmt)
            self.ge = zeros(ncards, dtype=float_fmt)
            self.F12 = zeros(ncards, dtype=float_fmt)
            self.strn = zeros(ncards, dtype=float_fmt)

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i

        mid = integer(card, 1, 'mid')
        if comment:
            self.set_comment(mid, comment)

        self.material_id[i] = mid
        self.e11[i] = double(card, 2, 'E11')    #: .. todo:: is this the correct default
        self.e22[i] = double(card, 3, 'E22')    #: .. todo:: is this the correct default
        self.nu12[i] = double_or_blank(card, 4, 'nu12', 0.0)

        self.g12[i] = double_or_blank(card, 5, 'g12', 0.0)
        self.g1z[i] = double_or_blank(card, 6, 'g1z', 1e8)
        self.g2z[i] = double_or_blank(card, 7, 'g2z', 1e8)
        self.rho[i] = double_or_blank(card, 8, 'rho', 0.0)
        self.a1[i] = double_or_blank(card, 9, 'a1', 0.0)
        self.a2[i] = double_or_blank(card, 10, 'a2', 0.0)
        self.tref[i] = double_or_blank(card, 11, 'tref', 0.0)
        self.Xt[i] = double_or_blank(card, 12, 'Xt', 0.0)
        self.Xc[i] = double_or_blank(card, 13, 'Xc', self.Xt[i])
        self.Yt[i] = double_or_blank(card, 14, 'Yt', 0.0)
        self.Yc[i] = double_or_blank(card, 15, 'Yc', self.Yt[i])
        self.S[i] = double_or_blank(card, 16, 'S', 0.0)
        self.ge[i] = double_or_blank(card, 17, 'ge', 0.0)
        self.F12[i] = double_or_blank(card, 18, 'F12', 0.0)
        self.strn[i] = double_or_blank(card, 19, 'strn', 0.0)
        assert len(card) <= 20, 'len(MAT8 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1


    def build(self):
        if self.n:
            i = self.material_id.argsort()
            self.material_id = self.material_id[i]
            self.e11 = self.e11[i]
            self.e22 = self.e22[i]
            self.nu12 = self.nu12[i]

            self.g12 = self.g12[i]
            self.g1z = self.g1z[i]
            self.g2z = self.g2z[i]
            self.rho = self.rho[i]
            self.a1 = self.a1[i]
            self.a2 = self.a2[i]
            self.tref = self.tref[i]
            self.Xt = self.Xt[i]
            self.Xc = self.Xc[i]
            self.Yt = self.Yt[i]
            self.Yc = self.Yc[i]
            self.S = self.S[i]
            self.ge = self.ge[i]
            self.F12 = self.F12[i]
            self.strn = self.strn[i]

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

            #card = ['$MAT1', 'mid', 'E', 'G', 'nu', 'rho', 'a', 'tref', 'ge']
            #bdf_file.write(print_card_8(card))
            #card = ['$', 'st', 'sc', 'ss', 'mcsid']
            #bdf_file.write(print_card_8(card))
            for (mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, tref,
                 Xt, Xc, Yt, Yc, S, ge, F12, strn) in zip(
                               self.material_id[i], self.e11[i], self.e22[i], self.nu12[i], self.g12[i],
                               self.g1z[i], self.g2z[i], self.rho[i], self.a1[i], self.a2[i], self.tref[i],
                               self.Xt[i], self.Xc[i], self.Yt[i], self.Yc[i], self.S[i], self.ge[i],
                               self.F12[i], self.strn[i]):
                if mid in self._comments:
                    bdf_file.write(self._comments[mid])

                g12 = set_blank_if_default(g12, 0.)
                g1z = set_blank_if_default(g1z, 1e8)
                g2z = set_blank_if_default(g2z, 1e8)

                rho = set_blank_if_default(rho, 0.0)
                a1 = set_blank_if_default(a1, 0.0)
                a2 = set_blank_if_default(a2, 0.0)
                tref = set_blank_if_default(tref, 0.0)

                Xt = set_blank_if_default(Xt, 0.)
                Yt = set_blank_if_default(Yt, 0.)

                Xc = set_blank_if_default(Xc, Xt)
                Yc = set_blank_if_default(Yc, Yt)

                S = set_blank_if_default(S, 0.0)
                ge = set_blank_if_default(ge, 0.0)
                F12 = set_blank_if_default(F12, 0.0)
                strn = set_blank_if_default(strn, 0.0)

                card = ['MAT8', mid, e11, e22, nu12, g12, g1z, g2z,
                        rho, a1, a2, tref,
                        Xt, Xc, Yt, Yc, S, ge, F12, strn]
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
        obj = MAT8(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]

        obj.material_id = self.material_id[i]
        obj.e11 = self.e11[i]
        obj.e22 = self.e22[i]
        obj.nu12 = self.nu12[i]

        obj.g12 = self.g12[i]
        obj.g1z = self.g1z[i]
        obj.g2z = self.g2z[i]
        obj.rho = self.rho[i]
        obj.a1 = self.a1[i]
        obj.a2 = self.a2[i]
        obj.tref = self.tref[i]
        obj.Xt = self.Xt[i]
        obj.Xc = self.Xc[i]
        obj.Yt = self.Yt[i]
        obj.Yc = self.Yc[i]
        obj.S = self.S[i]
        obj.ge = self.ge[i]
        obj.F12 = self.F12[i]
        obj.strn = self.strn[i]

        return obj
