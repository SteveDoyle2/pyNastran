from __future__ import annotations
from itertools import zip_longest
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank, string_or_blank)
from pyNastran.bdf.cards.materials import mat1_E_G_nu, get_G_default, set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import Material
from pyNastran.dev.bdf_vectorized3.cards.write_utils import get_print_card, array_str, array_default_int

if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class MAT1(Material):
    """
    Defines the material properties for linear isotropic materials.

    +------+-----+-----+-----+-------+-----+------+------+-----+
    |   1  |  2  | 3   | 4   |   5   |  6  |  7   |  8   |  9  |
    +======+=====+=====+=====+=======+=====+======+======+=====+
    | MAT1 | MID |  E  |  G  |  NU   | RHO |  A   | TREF | GE  |
    +------+-----+-----+-----+-------+-----+------+------+-----+
    |      | ST  | SC  | SS  | MCSID |     |      |      |     |
    +------+-----+-----+-----+-------+-----+------+------+-----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n = 0
        self.material_id = np.array([], dtype='int32')
        self.E = np.array([], dtype='float64')
        self.G = np.array([], dtype='float64')
        self.nu = np.array([], dtype='float64')
        self.rho = np.array([], dtype='float64')
        self.alpha = np.array([], dtype='float64')
        self.tref = np.array([], dtype='float64')
        self.ge = np.array([], dtype='float64')
        self.St = np.array([], dtype='float64')
        self.Sc = np.array([], dtype='float64')
        self.Ss = np.array([], dtype='float64')
        self.mcsid = np.array([], dtype='int32')

    def add(self, mid: int, E: float, G: float, nu: float,
            rho: float=0.0, alpha: float=0.0, tref: float=0.0, ge: float=0.0,
            St: float=0.0, Sc: float=0.0, Ss: float=0.0,
            mcsid: int=0, comment: str=''):
        """
        Creates a MAT1 card

        Parameters
        ----------
        mid : int
            material id
        E : float / None
            Young's modulus
        G : float / None
            Shear modulus
        nu : float / None
            Poisson's ratio
        rho : float; default=0.
            density
        alpha : float; default=0.
            coefficient of thermal expansion
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            damping coefficient
        St / Sc / Ss : float; default=0.
            tensile / compression / shear allowable
        mcsid : int; default=0
            material coordinate system id
            used by PARAM,CURV
        comment : str; default=''
            a comment for the card

        If E, G, or nu is None (only 1), it will be calculated

        """
        E, G, nu = mat1_E_G_nu(E, G, nu)

        self.cards.append((mid, E, G, nu, rho, alpha, tref, ge, St, Sc, Ss,
                           mcsid, comment))
        self.n += 1
        #self.material_id = np.hstack([self.material_id, [mid]])
        #self.E = np.hstack([self.E, [E]])
        #self.G = np.hstack([self.G, [G]])
        #self.nu = np.hstack([self.nu, [nu]])
        #self.rho = np.hstack([self.rho, [rho]])
        #self.alpha = np.hstack([self.alpha, [alpha]])
        #self.tref = np.hstack([self.tref, [tref]])
        #self.St = np.hstack([self.St, [St]])
        #self.Sc = np.hstack([self.Sc, [Sc]])
        #self.Ss = np.hstack([self.Ss, [Ss]])
        #self.mcsid = np.hstack([self.mcsid, [mcsid]])
        #self.n += 1

    def add_card(self, card: BDFCard, comment: str=''):
        mid = integer(card, 1, 'mid')
        E = double_or_blank(card, 2, 'E')
        G = double_or_blank(card, 3, 'G')
        nu = double_or_blank(card, 4, 'nu')
        if E is None and G is None and nu is None:
            E = G = nu = 0.
            self.model.log.warning(f'MAT1; mid={mid} E=G=nu=0.')
        else:
            E, G, nu = mat1_E_G_nu(E, G, nu)

        rho = double_or_blank(card, 5, 'rho', default=0.)
        alpha = double_or_blank(card, 6, 'a', default=0.0)
        tref = double_or_blank(card, 7, 'tref', default=0.0)
        ge = double_or_blank(card, 8, 'ge', default=0.0)
        St = double_or_blank(card, 9, 'St', default=0.0)
        Sc = double_or_blank(card, 10, 'Sc', default=0.0)
        Ss = double_or_blank(card, 11, 'Ss', default=0.0)
        mcsid = integer_or_blank(card, 12, 'mcsid', default=0)
        assert len(card) <= 13, f'len(MAT1 card) = {len(card):d}\ncard={card}'
        self.cards.append((mid, E, G, nu, rho, alpha, tref, ge, St, Sc, Ss,
                           mcsid, comment))
        self.n += 1

    def add_op2_data(self, data, comment: str=''):
        """
        Adds a MAT1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        mid = data[0]
        E = data[1]
        G = data[2]
        nu = data[3]
        rho = data[4]
        alpha = data[5]
        tref = data[6]
        ge = data[7]
        St = data[8]
        Sc = data[9]
        Ss = data[10]
        mcsid = data[11]
        #return MAT1(mid, e, g, nu, rho, a, tref, ge,
                    #St, Sc, Ss, mcsid, comment=comment)
        assert rho is not None, rho
        assert mcsid is not None, mcsid
        self.cards.append((mid, E, G, nu, rho, alpha, tref, ge, St, Sc, Ss,
                           mcsid, comment))
        self.n += 1

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        material_id = np.zeros(ncards, dtype='int32')
        E = np.zeros(ncards, dtype='float64')
        G = np.zeros(ncards, dtype='float64')
        nu = np.zeros(ncards, dtype='float64')
        rho = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')
        tref = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')
        Ss = np.zeros(ncards, dtype='float64')
        St = np.zeros(ncards, dtype='float64')
        Sc = np.zeros(ncards, dtype='float64')
        mcsid = np.zeros(ncards, dtype='int32')

        for i, card in enumerate(self.cards):
            (mid, Ei, Gi, nui, rhoi, alphai, trefi, gei,
             Sti, Sci, Ssi,
             mcsidi, comment) = card
            material_id[i] = mid
            E[i] = Ei
            G[i] = Gi
            nu[i] = nui
            rho[i] = rhoi
            alpha[i] = alphai
            tref[i] = trefi
            ge[i] = gei
            Ss[i] = Ssi
            St[i] = Sti
            Sc[i] = Sci
            mcsid[i] = mcsidi
        self._save(material_id, E, G, nu, rho, alpha, tref, ge,
                   Ss, St, Sc, mcsid)
        self.sort()
        self.cards = []

    def _save(self, material_id, E, G, nu,
              rho, alpha, tref, ge,
              Ss, St, Sc, mcsid):
        if len(self.material_id) > 0:
            material_id = np.hstack([self.material_id, material_id])
            E = np.hstack([self.E, E])
            G = np.hstack([self.G, G])
            nu = np.hstack([self.nu, nu])
            rho = np.hstack([self.rho, rho])
            alpha = np.hstack([self.alpha, alpha])
            tref = np.hstack([self.tref, tref])
            ge = np.hstack([self.ge, ge])
            material_id = np.hstack([self.material_id, material_id])
            Ss = np.hstack([self.Ss, Ss])
            St = np.hstack([self.St, St])
            Sc = np.hstack([self.Sc, Sc])
            mcsid = np.hstack([self.mcsid, mcsid])
        self.material_id = material_id
        self.E = E
        self.G = G
        self.nu = nu
        self.rho = rho
        self.alpha = alpha
        self.tref = tref
        self.ge = ge
        self.Ss = Ss
        self.St = St
        self.Sc = Sc
        self.mcsid = mcsid

    def __apply_slice__(self, mat: MAT1, i: np.ndarray) -> None:
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.E = self.E[i]
        mat.G = self.G[i]
        mat.nu = self.nu[i]
        mat.rho = self.rho[i]
        mat.alpha = self.alpha[i]
        mat.tref = self.tref[i]
        mat.ge = self.ge[i]
        mat.Ss = self.Ss[i]
        mat.St = self.St[i]
        mat.Sc = self.Sc[i]
        mat.mcsid = self.mcsid[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def write(self, size: int=8) -> str:
        if len(self.material_id) == 0:
            return ''

        max_int = max(self.material_id.max(), self.mcsid.max())
        print_card = get_print_card(size, max_int)
        lines = []
        for mid, e, g, nu, rho, alpha, tref, ge, Ss, St, Sc, \
            mcsid in zip_longest(self.material_id, self.E, self.G, self.nu, self.rho,
                                 self.alpha, self.tref, self.ge,
                                 self.Ss, self.St, self.Sc, self.mcsid):
            g_default = get_G_default(e, g, nu)
            G = set_blank_if_default(g, g_default)

            rho = set_blank_if_default(rho, 0.)
            a = set_blank_if_default(alpha, 0.)
            tref = set_blank_if_default(tref, 0.)
            ge = set_blank_if_default(ge, 0.)

            if [St, Sc, Ss, mcsid] == [0., 0., 0., 0]:
                list_fields = ['MAT1', mid, e, G, nu, rho, a, tref, ge]
            else:
                St = set_blank_if_default(St, 0.)
                Sc = set_blank_if_default(Sc, 0.)
                Ss = set_blank_if_default(Ss, 0.)
                mcsid = set_blank_if_default(mcsid, 0)
                list_fields = ['MAT1', mid, e, G, nu, rho, a, tref, ge,
                               St, Sc, Ss, mcsid]
            lines.append(print_card(list_fields))
        return ''.join(lines)

    def s33(self):
        """
        ei2 = [
            [  1 / e, -nu / e,    0.],
            [-nu / e,   1 / e,    0.],
            [     0.,      0., 1 / g],
        ]
        """
        nmaterial = len(self.material_id)
        s33 = np.zeros((nmaterial, 3, 3), dtype='float64')
        s33[:, 0, 0] = s33[:, 1, 1] = 1 / self.E
        s33[:, 1, 0] = s33[:, 0, 1] = -self.nu / self.E
        s33[:, 2, 2] = 1 / self.G
        return s33


class MAT2(Material):
    """
    Defines the material properties for linear anisotropic materials for
    two-dimensional elements.

    +------+-------+-----+-----+------+-----+------+-----+-----+
    |   1  |   2   |  3  |  4  |  5   |  6  |  7   | 8   |  9  |
    +======+=======+=====+=====+======+=====+======+=====+=====+
    | MAT2 |  MID  | G11 | G12 | G13  | G22 | G23  | G33 | RHO |
    +------+-------+-----+-----+------+-----+------+-----+-----+
    |      |  A1   | A2  | A3  | TREF | GE  |  ST  | SC  | SS  |
    +------+-------+-----+-----+------+-----+------+-----+-----+
    |      | MCSID |     |     |      |     |      |     |     |
    +------+-------+-----+-----+------+-----+------+-----+-----+

    """
    def add(self, mid: float,
            G11: float, G12: float, G13: float,
            G22: float, G23: float, G33: float, rho: float=0.,
            a1: Optional[float]=None, a2: Optional[float]=None, a3: Optional[float]=None,
            tref: float=0., ge: float=0.,
            St: Optional[float]=None, Sc: Optional[float]=None, Ss: Optional[float]=None,
            mcsid: Optional[int]=None, comment: str='') -> MAT2:
        """Creates an MAT2 card"""
        if a1 is None:
            a1 = np.nan
        if a2 is None:
            a2 = np.nan
        if a3 is None:
            a3 = np.nan

        if St is None:
            St = np.nan
        if Sc is None:
            Sc = np.nan
        if Ss is None:
            Ss = np.nan

        if mcsid is None:
            mcsid = -1
        ge_matrix = [np.nan] * 6
        self.cards.append((mid, G11, G12, G13, G22, G23, G33, rho,
                           [a1, a2, a3], tref, ge, St, Sc, Ss,
                           mcsid, ge_matrix, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str=''):
        mid = integer(card, 1, 'mid')
        G11 = double_or_blank(card, 2, 'G11', default=0.0)
        G12 = double_or_blank(card, 3, 'G12', default=0.0)
        G13 = double_or_blank(card, 4, 'G13', default=0.0)
        G22 = double_or_blank(card, 5, 'G22', default=0.0)
        G23 = double_or_blank(card, 6, 'G23', default=0.0)
        G33 = double_or_blank(card, 7, 'G33', default=0.0)

        rho = double_or_blank(card, 8, 'rho', default=0.0)
        a1 = double_or_blank(card, 9, 'a1') # blank?
        a2 = double_or_blank(card, 10, 'a2') # blank?
        a3 = double_or_blank(card, 11, 'a3') # blank?
        tref = double_or_blank(card, 12, 'tref', default=0.0)
        ge = double_or_blank(card, 13, 'ge', default=0.0)
        St = double_or_blank(card, 14, 'St') # or blank?
        Sc = double_or_blank(card, 15, 'Sc') # or blank?
        Ss = double_or_blank(card, 16, 'Ss') # or blank?
        mcsid = integer_or_blank(card, 17, 'mcsid', default=-1)

        if len(card) > 18:
            ge_matrix = [
                double_or_blank(card, 18, 'ge11', default=0.0),
                double_or_blank(card, 19, 'ge12', default=0.0),
                double_or_blank(card, 20, 'ge13', default=0.0),
                double_or_blank(card, 21, 'ge22', default=0.0),
                double_or_blank(card, 22, 'ge23', default=0.0),
                double_or_blank(card, 23, 'ge33', default=0.0),
            ]
            assert len(card) <= 24, f'len(MAT2 card) = {len(card):d}\ncard={card}'
            #self.ge_matrix[i, :] = ge_matrix
        else:
            ge_matrix = [np.nan] * 6
            assert len(card) <= 18, f'len(MAT2 card) = {len(card):d}\ncard={card}'
        self.cards.append((mid, G11, G12, G13, G22, G23, G33, rho,
                           [a1, a2, a3], tref, ge, St, Sc, Ss, mcsid, ge_matrix, comment))
        self.n += 1

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        material_id = np.zeros(ncards, dtype='int32')
        G11 = np.zeros(ncards, dtype='float64')
        G12 = np.zeros(ncards, dtype='float64')
        G13 = np.zeros(ncards, dtype='float64')
        G22 = np.zeros(ncards, dtype='float64')
        G23 = np.zeros(ncards, dtype='float64')
        G33 = np.zeros(ncards, dtype='float64')

        rho = np.zeros(ncards, dtype='float64')
        alpha = np.zeros((ncards, 3), dtype='float64')
        tref = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')
        Ss = np.zeros(ncards, dtype='float64')
        St = np.zeros(ncards, dtype='float64')
        Sc = np.zeros(ncards, dtype='float64')
        mcsid = np.zeros(ncards, dtype='int32')
        ge_matrix = np.full((ncards, 6), np.nan, dtype='float64')

        for i, card in enumerate(self.cards):
            (mid, G11i, G12i, G13i, G22i, G23i, G33i, rhoi,
             alphas, trefi, gei, Sti, Sci, Ssi, mcsidi, ge_matrixi, comment) = card

            material_id[i] = mid
            G11[i] = G11i
            G12[i] = G12i
            G13[i] = G13i
            G22[i] = G22i
            G23[i] = G23i
            G33[i] = G33i

            rho[i] = rhoi
            alpha[i] = alphas
            tref[i] = trefi
            ge[i] = gei
            Ss[i] = Ssi
            St[i] = Sti
            Sc[i] = Sci
            mcsid[i] = mcsidi
            ge_matrix[i] = ge_matrixi
        self._save(material_id, G11, G12, G13, G22, G23, G33,
                   rho, alpha, tref, ge, Ss, St, Sc,
                   mcsid, ge_matrix)
        self.sort()
        self.cards = []

    def _save(self, material_id, G11, G12, G13, G22, G23, G33,
              rho, alpha, tref, ge, Ss, St, Sc,
              mcsid, ge_matrix):
        #print('calling MAT2 save')
        nmaterial = len(material_id)
        assert nmaterial > 0, nmaterial
        self.material_id = material_id
        self.G11 = G11
        self.G12 = G12
        #print(G12.shape, G11.shape)
        self.G13 = G13
        self.G22 = G22
        self.G23 = G23
        self.G33 = G33

        self.rho = rho
        self.alpha = alpha
        self.tref = tref
        self.ge = ge
        self.Ss = Ss
        self.St = St
        self.Sc = Sc
        self.mcsid = mcsid
        if ge_matrix is None:
            ge_matrix = np.full((nmaterial, 6), np.nan, dtype='float64')
        self.ge_matrix = ge_matrix
        self.n = len(material_id)

    def validate(self):
        assert isinstance(self.material_id, np.ndarray), self.material_id
        assert isinstance(self.G11, np.ndarray), self.G11
        assert isinstance(self.G12, np.ndarray), self.G12
        assert isinstance(self.G13, np.ndarray), self.G13
        assert isinstance(self.G22, np.ndarray), self.G22
        assert isinstance(self.G23, np.ndarray), self.G23
        assert isinstance(self.G33, np.ndarray), self.G33

        assert isinstance(self.rho, np.ndarray), self.rho
        assert isinstance(self.alpha, np.ndarray), self.alpha
        assert isinstance(self.tref, np.ndarray), self.tref
        assert isinstance(self.ge, np.ndarray), self.ge
        assert isinstance(self.Ss, np.ndarray), self.Ss
        assert isinstance(self.St, np.ndarray), self.St
        assert isinstance(self.Sc, np.ndarray), self.Sc
        assert isinstance(self.mcsid, np.ndarray), self.mcsid
        #if self.ge_matrix is None:
        assert isinstance(self.ge_matrix, np.ndarray), self.ge_matrix

        #print(self.material_id.dtype.name)
        #print(self.G11.dtype.name)
        #print(self.G12.dtype.name)
        #print(self.G13.dtype.name)
        #print(self.G22.dtype.name)
        #print(self.G23.dtype.name)
        #print(self.G33.dtype.name)

        #print(self.rho.dtype.name)
        #print(self.alpha.dtype.name)
        #print(self.tref.dtype.name)
        #print(self.ge.dtype.name)
        #print(self.Ss.dtype.name)
        #print(self.St.dtype.name)
        #print(self.Sc.dtype.name)
        #print(self.mcsid.dtype.name)

    def __apply_slice__(self, mat: MAT2, i: np.ndarray) -> None:
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.G11 = self.G11[i]
        mat.G12 = self.G12[i]
        mat.G22 = self.G22[i]
        mat.G33 = self.G33[i]
        mat.G13 = self.G13[i]
        mat.G23 = self.G23[i]
        mat.ge_matrix = self.ge_matrix[i, :]

        mat.rho = self.rho[i]
        mat.alpha = self.alpha[i, :]
        mat.tref = self.tref[i]
        mat.ge = self.ge[i]
        mat.Ss = self.Ss[i]
        mat.St = self.St[i]
        mat.Sc = self.Sc[i]
        mat.mcsid = self.mcsid[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def write(self, size: int=8) -> str:
        if len(self.material_id) == 0:
            return ''

        max_int = max(self.material_id.max(), self.mcsid.max())
        print_card = get_print_card(size, max_int)
        lines = []
        self.validate()
        #print(self.material_id)
        #print(self.G11)
        #print('G12', self.G12)
        #print(self.G13)
        #print(self.G22)
        #print(self.G23)
        #print(self.G33)
        #print(self.rho)
        #print(self.alpha)
        #print(self.tref)
        #print(self.ge)
        #print(self.Ss)
        #print(self.St)
        #print(self.Sc)
        #print(self.mcsid)
        for mid, G11, G12, G13, G22, G23, G33, \
            rho, (a1, a2, a3), tref, \
            ge, Ss, St, Sc, mcsid in zip_longest(self.material_id,
                                                 self.G11, self.G12, self.G13,
                                                 self.G22, self.G23, self.G33,
                                                 self.rho,
                                                 self.alpha, self.tref, self.ge,
                                                 self.Ss, self.St, self.Sc, self.mcsid):
            G11 = set_blank_if_default(G11, 0.0)
            G12 = set_blank_if_default(G12, 0.0)
            G13 = set_blank_if_default(G13, 0.0)
            G22 = set_blank_if_default(G22, 0.0)
            G23 = set_blank_if_default(G23, 0.0)
            G33 = set_blank_if_default(G33, 0.0)
            rho = set_blank_if_default(rho, 0.0)
            tref = set_blank_if_default(tref, 0.0)
            ge = set_blank_if_default(ge, 0.0)
            list_fields = [
                'MAT2', mid, G11, G12, G13, G22, G23, G33, rho,
                a1, a2, a3, tref, ge,
                St, Sc, Ss, mcsid]
            #if ge_matrix != [0., 0., 0., 0., 0., 0.]:
                #ge11 = set_blank_if_default(ge_matrix[0], 0.0)
                #ge12 = set_blank_if_default(ge_matrix[1], 0.0)
                #ge13 = set_blank_if_default(ge_matrix[2], 0.0)
                #ge22 = set_blank_if_default(ge_matrix[3], 0.0)
                #ge23 = set_blank_if_default(ge_matrix[4], 0.0)
                #ge33 = set_blank_if_default(ge_matrix[5], 0.0)
                #list_fields += [ge11, ge12, ge13, ge22, ge23, ge33]
            lines.append(print_card(list_fields))
        return ''.join(lines)


class MAT8(Material):
    """
    Defines the material property for an orthotropic material for
    isoparametric shell elements.

    +------+-----+-----+------+------+-----+-----+-----+-----+
    |  1   |  2  |  3  |  4   |  5   |  6  |  7  |  8  |  9  |
    +======+=====+=====+======+======+=====+=====+=====+=====+
    | MAT8 | MID | E1  |  E2  | NU12 | G12 | G1Z | G2Z | RHO |
    +------+-----+-----+------+------+-----+-----+-----+-----+
    |      | A1  |  A2 | TREF |  Xt  |  Xc |  Yt |  Yc |  S  |
    +------+-----+-----+------+------+-----+-----+-----+-----+
    |      | GE1 | F12 | STRN |      |     |     |     |     |
    +------+-----+-----+------+------+-----+-----+-----+-----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n = 0
        self.material_id = np.array([], dtype='int32')
        self.E11 = np.array([], dtype='float64')
        self.E22 = np.array([], dtype='float64')
        self.G12 = np.array([], dtype='float64')
        self.G13 = np.array([], dtype='float64')
        self.G23 = np.array([], dtype='float64')
        self.nu12 = np.array([], dtype='float64')

        self.rho = np.array([], dtype='float64')
        self.alpha = np.zeros((0, 2), dtype='float64')
        self.tref = np.array([], dtype='float64')
        self.ge = np.array([], dtype='float64')

        self.Xt = np.array([], dtype='float64')
        self.Xc = np.array([], dtype='float64')
        self.Yt = np.array([], dtype='float64')
        self.Yc = np.array([], dtype='float64')
        self.S = np.array([], dtype='float64')
        self.f12 = np.array([], dtype='float64')
        self.strn = np.array([], dtype='float64')

    def add(self, mid: int, e11: float, e22: float, nu12: float,
            g12: float=0.0, g1z: float=1e8, g2z: float=1e8,
            rho: float=0., a1: float=0., a2: float=0., tref: float=0.,
            xt: float=0., xc: float=None, yt: float=0., yc: float=None,
            s: float=0., ge: float=0., f12: float=0., strn: float=0.,
            comment: str='') -> MAT8:
        """Creates a MAT8 card"""
        xc = xc if xc is not None else xt
        yc = yc if yc is not None else yt
        #self.material_id = np.hstack([self.material_id, [mid]])
        #self.E11 = np.hstack([self.E11, [e11]])
        #self.E22 = np.hstack([self.E22, [e22]])
        #self.G12 = np.hstack([self.G12, [g12]])
        #self.G13 = np.hstack([self.G13, [g1z]])
        #self.G23 = np.hstack([self.G23, [g2z]])
        #self.nu12 = np.hstack([self.nu12, [nu12]])
        #self.rho = np.hstack([self.rho, [rho]])
        #self.alpha = np.vstack([self.alpha, [a1, a2]])
        #self.tref = np.hstack([self.tref, [tref]])
        #self.ge = np.hstack([self.ge, [ge]])

        #self.Xt = np.hstack([self.Xt, [Xt]])
        #self.Xc = np.hstack([self.Xc, [Xc]])
        #self.Yt = np.hstack([self.Yt, [Yt]])
        #self.Yc = np.hstack([self.Yc, [Yc]])
        #self.S = np.hstack([self.S, [S]])
        #self.f12 = np.hstack([self.f12, [F12]])
        #self.strn = np.hstack([self.strn, [strn]])
        self.cards.append((mid, e11, e22, nu12, g12, g1z, g2z, rho, [a1, a2], tref,
                           xt, xc, yt, yc, s, ge, f12, strn, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str=''):
        mid = integer(card, 1, 'mid')
        e11 = double(card, 2, 'E11')    #: .. todo:: is this the correct default
        e22 = double(card, 3, 'E22')    #: .. todo:: is this the correct default

        nu12 = double_or_blank(card, 4, 'nu12', default=0.0)

        g12 = double_or_blank(card, 5, 'g12', default=0.0)
        g1z = double_or_blank(card, 6, 'g1z', default=1e8)
        g2z = double_or_blank(card, 7, 'g2z', default=1e8)
        rho = double_or_blank(card, 8, 'rho', default=0.0)
        a1 = double_or_blank(card, 9, 'a1', default=0.0)
        a2 = double_or_blank(card, 10, 'a2', default=0.0)
        tref = double_or_blank(card, 11, 'tref', default=0.0)
        xt = double_or_blank(card, 12, 'Xt', default=0.0)
        xc = double_or_blank(card, 13, 'Xc', default=xt)
        yt = double_or_blank(card, 14, 'Yt', default=0.0)
        yc = double_or_blank(card, 15, 'Yc', default=yt)
        s = double_or_blank(card, 16, 'S', default=0.0)
        ge = double_or_blank(card, 17, 'ge', default=0.0)
        f12 = double_or_blank(card, 18, 'F12', default=0.0)
        strn = double_or_blank(card, 19, 'strn', default=0.0)
        assert len(card) <= 20, f'len(MAT8 card) = {len(card):d}\ncard={card}'
        self.cards.append((mid, e11, e22, nu12, g12, g1z, g2z, rho, [a1, a2], tref,
                           xt, xc, yt, yc, s, ge, f12, strn, comment))
        self.n += 1

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        material_id = np.zeros(ncards, dtype='int32')

        E11 = np.zeros(ncards, dtype='float64')
        E22 = np.zeros(ncards, dtype='float64')
        G12 = np.zeros(ncards, dtype='float64')
        G13 = np.zeros(ncards, dtype='float64')
        G23 = np.zeros(ncards, dtype='float64')
        nu12 = np.zeros(ncards, dtype='float64')
        rho = np.zeros(ncards, dtype='float64')
        alpha = np.zeros((ncards, 2), dtype='float64')
        tref = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')

        Xt = np.zeros(ncards, dtype='float64')
        Xc = np.zeros(ncards, dtype='float64')
        Yt = np.zeros(ncards, dtype='float64')
        Yc = np.zeros(ncards, dtype='float64')
        S = np.zeros(ncards, dtype='float64')
        f12 = np.zeros(ncards, dtype='float64')
        strn = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (mid, e11, e22, nu12i, g12, g1z, g2z, rhoi, alphai, trefi,
             xt, xc, yt, yc, s, gei, f12i, strni, comment) = card
            material_id[icard] = mid
            E11[icard] = e11
            E22[icard] = e22
            G12[icard] = g12
            G13[icard] = g1z
            G23[icard] = g2z
            nu12[icard] = nu12i
            rho[icard] = rhoi
            alpha[icard, :] = alphai  ## thermal expansion in 1 and 2 directions
            tref[icard] = trefi       ## reference temperature
            ge[icard] = gei           ## damping
            Xt[icard] = xt            ## x tension allowable
            Xc[icard] = xc            ## x compression allowable
            Yt[icard] = yt            ## y tension allowable
            Yc[icard] = yc            ## y compression allowable
            S[icard] = s              ## shear allowable
            f12[icard] = f12i
            strn[icard] = strni
        self._save(material_id, E11, E22, G12, G13, G23, nu12,
                   rho, alpha, tref, ge, Xt, Xc, Yt, Yc, S, f12, strn)
        self.sort()
        self.cards = []

    def _save(self, material_id, E11, E22, G12, G13, G23, nu12,
              rho, alpha, tref, ge,
              Xt, Xc, Yt, Yc, S, f12, strn):
        nmaterials = len(material_id)
        self.material_id = material_id
        self.E11 = E11
        self.E22 = E22
        self.G12 = G12
        self.G13 = G13
        self.G23 = G23
        self.nu12 = nu12
        self.rho = rho
        self.alpha = alpha
        self.tref = tref
        self.ge = ge

        self.Xt = Xt
        self.Xc = Xc
        self.Yt = Yt
        self.Yc = Yc
        self.S = S
        self.f12 = f12
        self.strn = strn
        self.n = nmaterials

    def __apply_slice__(self, mat: MAT8, i: np.ndarray) -> None:
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.E11 = self.E11[i]
        mat.E22 = self.E22[i]
        mat.G12 = self.G12[i]
        mat.G13 = self.G13[i]
        mat.G23 = self.G23[i]

        mat.nu12 = self.nu12[i]
        mat.rho = self.rho[i]
        mat.alpha = self.alpha[i, :]
        mat.tref = self.tref[i]
        mat.ge = self.ge[i]

        mat.Xt = self.Xt[i]
        mat.Xc = self.Xc[i]
        mat.Yt = self.Yt[i]
        mat.Yc = self.Yc[i]
        mat.S = self.S[i]
        mat.f12 = self.f12[i]
        mat.strn = self.strn[i]

    @property
    def a1(self) -> np.ndarray:
        return self.alpha[:, 0]
    @property
    def a2(self) -> np.ndarray:
        return self.alpha[:, 1]

    @a1.setter
    def a1(self, a1: np.ndarray) -> None:
        self.alpha[:, 0] = a1
    @a2.setter
    def a2(self, a2: np.ndarray) -> None:
        self.alpha[:, 1] = a2

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def write(self, size: int=8) -> str:
        if len(self.material_id) == 0:
            return ''
        lines = []
        for mid, e11, e22, g12, g13, g23, nu12, rho, \
            alpha, tref, ge, xt, xc, yt, yc, S, f12, strn, in zip_longest(self.material_id, self.E11, self.E22,
                                                                          self.G12, self.G13, self.G13, self.nu12,
                                                                          self.rho,
                                                                          self.alpha, self.tref, self.ge,
                                                                          self.Xt, self.Xc, self.Yt, self.Yc, self.S,
                                                                          self.f12, self.strn):
            a1, a2 = alpha
            G12 = set_blank_if_default(g12, 0.)
            G1z = set_blank_if_default(g13, 1e8)
            G2z = set_blank_if_default(g23, 1e8)

            rho = set_blank_if_default(rho, 0.0)
            a1 = set_blank_if_default(a1, 0.0)
            a2 = set_blank_if_default(a2, 0.0)
            tref = set_blank_if_default(tref, 0.0)

            xc = set_blank_if_default(xc, xt)
            yc = set_blank_if_default(yc, yt)

            xt = set_blank_if_default(xt, 0.)
            yt = set_blank_if_default(yt, 0.)

            S = set_blank_if_default(S, 0.0)
            ge = set_blank_if_default(ge, 0.0)
            f12 = set_blank_if_default(f12, 0.0)
            strn = set_blank_if_default(strn, 0.0)

            list_fields = ['MAT8', mid, e11, e22, nu12, G12, G1z,
                           G2z, rho, a1, a2, tref, xt, xc, yt, yc, S, ge, f12, strn]
            lines.append(print_card_8(list_fields))
        return ''.join(lines)


class MAT9(Material):
    """
    Defines the material properties for linear, temperature-independent,
    anisotropic materials for solid isoparametric elements

    .. seealso::  PSOLID entry description

    +------+-----+-----+-----+-----+-----+------+-----+-----+
    |   1  |  2  | 3   | 4   |  5  |  6  |  7   | 8   |  9  |
    +======+=====+=====+=====+=====+=====+======+=====+=====+
    | MAT9 | MID | G11 | G12 | G13 | G14 | G15  | G16 | G22 |
    +------+-----+-----+-----+-----+-----+------+-----+-----+
    |      | G23 | G24 | G25 | G26 | G33 | G34  | G35 | G36 |
    +------+-----+-----+-----+-----+-----+------+-----+-----+
    |      | G44 | G45 | G46 | G55 | G56 | G66  | RHO | A1  |
    +------+-----+-----+-----+-----+-----+------+-----+-----+
    |      | A2  | A3  | A4  | A5  | A6  | TREF | GE  |     |
    +------+-----+-----+-----+-----+-----+------+-----+-----+

    .. warning:: MSC 2020: gelist is not supported.
    """
    def add(self, mid: int,
            G11=0., G12=0., G13=0., G14=0., G15=0., G16=0.,
            G22=0., G23=0., G24=0., G25=0., G26=0.,
            G33=0., G34=0., G35=0., G36=0.,
            G44=0., G45=0., G46=0.,
            G55=0., G56=0., G66=0.,
            rho=0., A=None, tref=0., ge=0., comment='') -> None:
        if A is None:
            A = [0.] * 6
        self.cards.append((mid, G11, G12, G13, G14, G15, G16,
                           G22, G23, G24, G25, G26,
                           G33, G34, G35, G36,
                           G44, G45, G46,
                           G55, G56, G66, rho, A, tref, ge, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        mid = integer(card, 1, 'mid')
        G11 = double_or_blank(card, 2, 'G11', default=0.0)
        G12 = double_or_blank(card, 3, 'G12', default=0.0)
        G13 = double_or_blank(card, 4, 'G13', default=0.0)
        G14 = double_or_blank(card, 5, 'G14', default=0.0)
        G15 = double_or_blank(card, 6, 'G15', default=0.0)
        G16 = double_or_blank(card, 7, 'G16', default=0.0)
        G22 = double_or_blank(card, 8, 'G22', default=0.0)
        G23 = double_or_blank(card, 9, 'G23', default=0.0)
        G24 = double_or_blank(card, 10, 'G24', default=0.0)
        G25 = double_or_blank(card, 11, 'G25', default=0.0)
        G26 = double_or_blank(card, 12, 'G26', default=0.0)
        G33 = double_or_blank(card, 13, 'G33', default=0.0)
        G34 = double_or_blank(card, 14, 'G34', default=0.0)
        G35 = double_or_blank(card, 15, 'G35', default=0.0)
        G36 = double_or_blank(card, 16, 'G36', default=0.0)
        G44 = double_or_blank(card, 17, 'G44', default=0.0)
        G45 = double_or_blank(card, 18, 'G45', default=0.0)
        G46 = double_or_blank(card, 19, 'G46', default=0.0)
        G55 = double_or_blank(card, 20, 'G55', default=0.0)
        G56 = double_or_blank(card, 21, 'G56', default=0.0)
        G66 = double_or_blank(card, 22, 'G66', default=0.0)
        rho = double_or_blank(card, 23, 'rho', default=0.0)
        alpha = [double_or_blank(card, 24, 'A1', default=0.0),
                 double_or_blank(card, 25, 'A2', default=0.0),
                 double_or_blank(card, 26, 'A3', default=0.0),
                 double_or_blank(card, 27, 'A4', default=0.0),
                 double_or_blank(card, 28, 'A5', default=0.0),
                 double_or_blank(card, 29, 'A6', default=0.0)]
        tref = double_or_blank(card, 30, 'tref', default=0.0)
        ge = double_or_blank(card, 31, 'ge', default=0.0)
        assert len(card) <= 32, f'len(MAT9 card) = {len(card):d}\ncard={card}'
        self.cards.append((mid, G11, G12, G13, G14, G15, G16,
                           G22, G23, G24, G25, G26,
                           G33, G34, G35, G36,
                           G44, G45, G46,
                           G55, G56, G66, rho, alpha, tref, ge, comment))
        self.n += 1

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        material_id = np.zeros(ncards, dtype='int32')
        G11 = np.zeros(ncards, dtype='float64')
        G12 = np.zeros(ncards, dtype='float64')
        G13 = np.zeros(ncards, dtype='float64')
        G14 = np.zeros(ncards, dtype='float64')
        G15 = np.zeros(ncards, dtype='float64')
        G16 = np.zeros(ncards, dtype='float64')
        G22 = np.zeros(ncards, dtype='float64')
        G23 = np.zeros(ncards, dtype='float64')
        G24 = np.zeros(ncards, dtype='float64')
        G25 = np.zeros(ncards, dtype='float64')
        G26 = np.zeros(ncards, dtype='float64')

        G33 = np.zeros(ncards, dtype='float64')
        G34 = np.zeros(ncards, dtype='float64')
        G35 = np.zeros(ncards, dtype='float64')
        G36 = np.zeros(ncards, dtype='float64')

        G44 = np.zeros(ncards, dtype='float64')
        G45 = np.zeros(ncards, dtype='float64')
        G46 = np.zeros(ncards, dtype='float64')

        G55 = np.zeros(ncards, dtype='float64')
        G56 = np.zeros(ncards, dtype='float64')
        G66 = np.zeros(ncards, dtype='float64')

        rho = np.zeros(ncards, dtype='float64')
        alpha = np.zeros((ncards, 6), dtype='float64')
        tref = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (mid, G11i, G12i, G13i, G14i, G15i, G16i,
             G22i, G23i, G24i, G25i, G26i,
             G33i, G34i, G35i, G36i,
             G44i, G45i, G46i,
             G55i, G56i, G66i, rhoi, alphai, trefi, gei, comment) = card
            material_id[icard] = mid
            G11[icard] = G11i
            G12[icard] = G12i
            G13[icard] = G13i
            G14[icard] = G14i
            G15[icard] = G15i
            G16[icard] = G16i
            G22[icard] = G22i
            G23[icard] = G23i
            G24[icard] = G24i
            G25[icard] = G25i
            G26[icard] = G26i

            G33[icard] = G33i
            G34[icard] = G34i
            G35[icard] = G35i
            G36[icard] = G36i

            G44[icard] = G44i
            G45[icard] = G45i
            G46[icard] = G46i

            G55[icard] = G55i
            G56[icard] = G56i
            G66[icard] = G66i

            rho[icard] = rhoi
            alpha[icard] = alphai
            tref[icard] = trefi
            ge[icard] = gei
        self._save(material_id, G11, G12, G13, G14, G15, G16,
                   G22, G23, G24, G25, G26,
                   G33, G34, G35, G36,
                   G44, G45, G46,
                   G55, G56,
                   G66, rho, alpha, tref, ge, ge_list=None)
        self.sort()
        self.cards = []

    def _save(self, material_id, G11, G12, G13, G14, G15, G16,
              G22, G23, G24, G25, G26,
              G33, G34, G35, G36,
              G44, G45, G46,
              G55, G56,
              G66, rho, alpha, tref, ge, ge_list):
        nmaterials = len(material_id)
        self.material_id = material_id
        self.G11 = G11
        self.G12 = G12
        self.G13 = G13
        self.G14 = G14
        self.G15 = G15
        self.G16 = G16
        self.G22 = G22
        self.G23 = G23
        self.G24 = G24
        self.G25 = G25
        self.G26 = G26

        self.G33 = G33
        self.G34 = G34
        self.G35 = G35
        self.G36 = G36

        self.G44 = G44
        self.G45 = G45
        self.G46 = G46

        self.G55 = G55
        self.G56 = G56
        self.G66 = G66

        self.rho = rho
        self.alpha = alpha
        self.tref = tref
        self.ge = ge
        if ge_list is None:
            ge_list = np.full((nmaterials, 21), np.nan, dtype=ge.dtype)
        self.ge_list = ge_list
        self.n = nmaterials

    def __apply_slice__(self, mat: MAT9, i: np.ndarray) -> None:
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.G11 = self.G11[i]
        mat.G12 = self.G12[i]
        mat.G13 = self.G13[i]
        mat.G14 = self.G14[i]
        mat.G15 = self.G15[i]
        mat.G16 = self.G16[i]
        mat.G22 = self.G22[i]
        mat.G23 = self.G23[i]
        mat.G24 = self.G24[i]
        mat.G25 = self.G25[i]
        mat.G26 = self.G26[i]

        mat.G33 = self.G33[i]
        mat.G34 = self.G34[i]
        mat.G35 = self.G35[i]
        mat.G36 = self.G36[i]
        mat.G44 = self.G44[i]
        mat.G45 = self.G45[i]
        mat.G46 = self.G46[i]
        mat.G55 = self.G55[i]
        mat.G56 = self.G56[i]
        mat.G66 = self.G66[i]

        mat.rho = self.rho[i]
        mat.alpha = self.alpha[i, :]
        mat.tref = self.tref[i]
        mat.ge = self.ge[i]
        mat.ge_list = self.ge_list[i, :]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def write(self, size: int=8) -> str:
        if len(self.material_id) == 0:
            return ''
        lines = []
        # TODO: ge_list
        for mid, G11, G12, G13, G14, G15, G16, \
        G22, G23, G24, G25, G26, \
        G33, G34, G35, G36, \
        G44, G45, G46, \
        G55, G56, \
        G66, rho, alpha, tref, \
            ge in zip_longest(self.material_id,
                              self.G11, self.G12, self.G13, self.G14, self.G15, self.G16,
                              self.G22, self.G23, self.G24, self.G25, self.G26,
                              self.G33, self.G34, self.G35, self.G36,
                              self.G44, self.G45, self.G46,
                              self.G55, self.G56,
                              self.G66,
                              self.rho, self.alpha, self.tref, self.ge):
            A = []
            for a in alpha:
                a = set_blank_if_default(a, 0.0)
                A.append(a)

            rho = set_blank_if_default(rho, 0.0)
            tref = set_blank_if_default(tref, 0.0)
            ge = set_blank_if_default(ge, 0.0)
            list_fields = (['MAT9', mid, G11, G12, G13, G14,
                            G15, G16, G22, G23, G24, G25,
                            G26, G33, G34, G35, G36, G44,
                            G45, G46, G55, G56, G66, rho]
                           + A + [tref, ge])
            lines.append(print_card_8(list_fields))
        return ''.join(lines)
