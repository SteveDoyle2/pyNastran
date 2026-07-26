from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, string)

from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points,)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, parse_check, save_ifile_comment)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, get_print_card_size)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class CAERO7(VectorizedBaseCard):
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')

    def add(self, eid: int, label: str,
            p1: np.ndarray, x12: float,
            p4: np.ndarray, x43: float,
            lchord_root: int=0, attach_root: int=0, achord_root: int=0,
            lchord_tip: int=0, attach_tip: int=0, achord_tip: int=0,
            cp: int=0,
            nspan: int=0, lspan: int=0,
            nchord: int=0,
            p_airfoil: int=0, ztaic: int=0,
            ifile: int=0, comment: str='') -> int:
        p_airfoil = p_airfoil if p_airfoil is not None else 0
        ztaic = ztaic if ztaic is not None else 0
        card = (eid, label,
                p1, x12, lchord_root, attach_root, achord_root,
                p4, x43, lchord_tip, attach_tip, achord_tip,
                cp, nspan, lspan, nchord,
                p_airfoil, ztaic, ifile, comment)
        assert len(card) == 20, card
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a CAERO7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        name = string(card, 2, 'name')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nspan = integer_or_blank(card, 4, 'nspan', default=0)
        nchord = integer_or_blank(card, 5, 'nchord', default=0)
        lspan = integer_or_blank(card, 6, 'aefact_lchord', default=0)
        ztaic = integer_or_blank(card, 7, 'ztaic', default=0)
        p_airfoil = integer_or_blank(card, 8, 'aefact', default=0)
        #assert cp == 0
        #igroup = integer(card, 8, 'igid')

        x1 = double_or_blank(card, 9, 'x1', default=0.0)
        y1 = double_or_blank(card, 10, 'y1', default=0.0)
        z1 = double_or_blank(card, 11, 'z1', default=0.0)
        p1 = np.array([x1, y1, z1])
        x12 = double_or_blank(card, 12, 'x12', default=0.)
        lchord_root = integer_or_blank(card, 13, 'lchord_root', default=0)
        attach_root = integer_or_blank(card, 14, 'attach_root', default=0)
        achord_root = integer_or_blank(card, 15, 'achord_root', default=0)

        x4 = double_or_blank(card, 17, 'x4', default=0.0)
        y4 = double_or_blank(card, 18, 'y4', default=0.0)
        z4 = double_or_blank(card, 19, 'z4', default=0.0)
        p4 = np.array([x4, y4, z4])
        x43 = double_or_blank(card, 20, 'x43', default=0.)

        lchord_tip = integer_or_blank(card, 21, 'lchord_tip', default=0)
        attach_tip = integer_or_blank(card, 22, 'attach_tip', default=0)
        achord_tip = integer_or_blank(card, 23, 'achord_tip', default=0)

        assert len(card) <= 23, f'len(CAERO7 card) = {len(card):d}\ncard={card}'
        #return CAERO7(eid, name, p1, x12, p4, x43,
                      #cp=cp, nspan=nspan, nchord=nchord, lspan=lspan,
                      #p_airfoil=p_airfoil, ztaic=ztaic,
                      #comment=comment)
        card = (eid, name,
                p1, x12, lchord_root, attach_root, achord_root,
                p4, x43, lchord_tip, attach_tip, achord_tip,
                cp, nspan, lspan, nchord,
                p_airfoil, ztaic, ifile, comment)
        assert len(card) == 20, card
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        p1 = np.zeros((ncards, 3), dtype='float64')
        p4 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        x43 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        nspan = np.zeros(ncards, dtype='int32')
        nchord = np.zeros(ncards, dtype='int32')
        lspan = np.zeros(ncards, dtype='int32')
        lchord = np.zeros(ncards, dtype='int32')
        ztaic = np.zeros(ncards, dtype='int32')
        p_airfoil = np.zeros(ncards, dtype='int32')
        comment = {}
        for icard, card in enumerate(self.cards):
            (eid, labeli,
             p1i, x12i, lchord_rooti, attach_rooti, achord_rooti,
             p4i, x43i, lchord_tipi, attach_tipi, achord_tipi,
             cpi, nspani, lspani, nchordi,
             p_airfoili, ztaici, ifilei, commenti) = card
            assert len(labeli) <= 8, f'label={labeli!r}'
            ifile[icard] = ifilei
            if commenti:
                comment[eid] = commenti
            element_id[icard] = eid
            label[icard] = labeli
            p1[icard, :] = p1i
            x12[icard] = x12i
            p4[icard, :] = p4i
            x43[icard] = x43i
            cp[icard] = cpi
            nspan[icard] = nspani
            lspan[icard] = lspani
            nchord[icard] = nchordi
            #lchord[icard] = lchordi
            ztaic[icard] = ztaici
            p_airfoil[icard] = p_airfoili
        self._save(element_id, label, p1, p4, x12, x43, cp,
                   nspan, lspan, nchord, lchord, ztaic, p_airfoil,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, element_id, label, p1, p4, x12, x43, cp,
              nspan, lspan, nchord, lchord, ztaic, p_airfoil,
              ifile=None, comment=None):
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.element_id):
            ifile = np.hstack([self.ifile, ifile])
            element_id = np.hstack([self.element_id, element_id])
            label = np.hstack([self.label, label])
            p1 = np.vstack([self.p1, p1])
            p4 = np.vstack([self.p4, p4])
            x12 = np.vstack([self.x12, x12])
            x43 = np.vstack([self.x43, x43])
            cp = np.hstack([self.cp, cp])
            nspan = np.hstack([self.nspan, nspan])
            lspan = np.hstack([self.lspan, lspan])
            nchord = np.hstack([self.nchord, nchord])
            lchord = np.hstack([self.lchord, lchord])
            ztaic = np.hstack([self.ztaic, ztaic])
            p_airfoil = np.hstack([self.p_airfoil, p_airfoil])

        self.ifile = ifile
        self.element_id = element_id
        self.label = label
        self.p1 = p1
        self.p4 = p4
        self.x12 = x12
        self.x43 = x43
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.nchord = nchord
        self.lchord = lchord
        self.ztaic = ztaic
        self.p_airfoil = p_airfoil
        self.n = len(element_id)

    @property
    def max_id(self) -> int:
        arrays = [
            self.element_id, # self.property_id,
            self.lchord, self.lspan,
            self.ztaic, self.p_airfoil,
        ]
        # assert isinstance(self.property_id[0], integer_types), self.property_id
        return max([array.max() for array in arrays])

    def __apply_slice__(self, elem: CAERO7, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.label = self.label[i]
        elem.p1 = self.p1[i, :]
        elem.p4 = self.p4[i, :]
        elem.x12 = self.x12[i]
        elem.x43 = self.x43[i]
        elem.cp = self.cp[i]
        elem.nspan = self.nspan[i]
        elem.lspan = self.lspan[i]
        elem.nchord = self.nchord[i]
        elem.lchord = self.lchord[i]
        elem.ztaic = self.ztaic[i]
        elem.p_airfoil = self.p_airfoil[i]

    @property
    def shape(self) -> tuple[np.ndarray, np.ndarray]:
        nchord = self.nchord.copy()
        nspan = self.nspan.copy()

        ichord = (nchord == 0)
        ispan = (nspan == 0)

        aefact = self.model.aefact
        if np.any(ichord):
            jchord = np.where(ichord)[0]
            lchord = self.lchord[jchord]
            aefacti = aefact.slice_card_by_aefact_id(lchord)
            nchord[ichord] = aefacti.nfractions - 1  # points -> boxes

        if np.any(ispan):
            jspan = np.where(ispan)[0]
            lspan = self.lspan[jspan]
            aefacti = aefact.slice_card_by_aefact_id(lspan)
            nspan[ispan] = aefacti.nfractions - 1  # points -> boxes
        return nchord, nspan

    def get_panel_npoints_nelements(self) -> tuple[int, int]:
        nchord, nspan = self.shape

        # nchord, nspan are number of boxes in the different directions
        npoints = (nchord + 1) * (nspan + 1)
        nelements = nchord * nspan
        return npoints, nelements

    def panel_points_elements(self) -> tuple[np.ndarray, np.ndarray]:
        points = []
        elements = []

        icp = (self.cp != 0)
        if np.any(icp):
            coord = self.model.coord
            raise NotImplementedError('CAERO7: cp')
        else:
            p1 = self.p1
            p4 = self.p4
            p2 = p1.copy()
            p3 = p4.copy()
        p2[:, 0] = p2[:, 0] + self.x12
        p4[:, 0] = p4[:, 0] + self.x43

        x, y = self.xy()
        ipoint = 0
        for p1i, p2i, p3i, p4i, xi, yi in zip(p1, p2, p3, p4, x, y):
            pointsi, elementsi = points_elements_from_quad_points(p1i, p4i, p3i, p2i, yi, xi, dtype='int32')
            points.append(pointsi)
            elements.append(elementsi + ipoint)
            ipoint += len(pointsi)
        return points, elements

    def validate(self):
        msg = ''
        is_failed = False
        #if not isinstance(self.p1, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True
        #if not isinstance(self.p4, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True

        element_id = self.element_id

        ibad = (self.x12 <= 0.)
        if np.any(ibad):
            msg += 'X12 and must be greater than or equal to 0\n'
            msg += f'  element_id={element_id[ibad]}\n   x12={self.x12[ibad]}'
            is_failed = True
        #if self.x43 <= 0.:
            #msg += 'X43=%s and must be greater than or equal to 0\n' % (self.x43)
            #is_failed = True

        ibad = (self.nspan == 0) & (self.lspan == 0)
        if np.any(ibad):
            msg += 'NSPAN or LSPAN must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nspan != 0) & (self.lspan != 0)
        if np.any(ibad):
            msg += 'Either NSPAN or LSPAN must 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nchord == 0) & (self.lchord == 0)
        if np.any(ibad):
            msg += 'NCHORD or LCHORD must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nchord={self.nchord[ibad]}\n  lchord={self.lchord}\n'
            is_failed = True

        ibad = (self.nchord != 0) & (self.lchord != 0)
        if np.any(ibad):
            msg += 'Either NCHORD or LCHORD must 0\n'
            msg += f'  element_id={element_id}\n  nchord={self.nchord[ibad]}\n  lchord={self.lchord}\n'
            is_failed = True
        if is_failed:
            msg += str(self)
            #msg += CAERO7_MSG
            raise ValueError(msg)

        #assert self.p1.shape[1] == 3, 'p1=%s' % self.p1.shape
        #assert self.p4.shape[1] == 3, 'p4=%s' % self.p4.shape

        ## calculating area; assuming coordinate transformations don't matter
        #if 1:
            #p1 = self.p1
            #p4 = self.p4
            #d12 = np.zeros(p1.shape)
            #d12[:, 0] = self.x12
            #d43 = np.zeros(p1.shape)
            #d43[:, 0] = self.x43
            #p2 = p1 + d12
            #p3 = p4 + d43

            #a = p3 - p1
            #b = p4 - p2
            #area = np.linalg.norm(np.cross(a, b), axis=1)
            #assert len(area) == p1.shape[0]
            #ibad = (area < 0.0)
            #if np.any(ibad):
                #msg += 'Either NCHORD or LCHORD must 0\n'
                #msg += f'eid={self.element_id[ibad,:]} p1={p1[ibad,:]} p2={p2[ibad,:]} p3={p3[ibad,:]} p4={p4[ibad,:]} area={area[ibad]}\n'
                #raise RuntimeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def flip_normal_by_element_id(self, element_id=None):
        """flips the CAERO1 normal vector"""
        i = self.index(element_id)
        self.flip_normal_by_index(i)

    def flip_normal_by_index(self, i: np.ndarray) -> None:
        """flips the CAERO1 normal vector"""
        self.p1[i, :], self.p4[i, :] = self.p4[i, :], self.p1[i, :]
        self.x12[i, :], self.x43[i, :] = self.x43[i, :], self.x12[i, :]

    def xy(self) -> tuple[list[np.ndarray], list[np.ndarray]]:
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        #xy = []
        x = []
        y = []
        for eid, nchord, lchord, nspan, lspan in zip(self.element_id,
                                                     self.nchord, self.lchord,
                                                     self.nspan, self.lspan):
            if nchord == 0:
                lchord_ref = self.model.aefact.slice_card_by_id(lspan)
                xi = lchord_ref.fractions
                nchord = len(xi) - 1
            else:
                xi = np.linspace(0., 1., nchord + 1)

            if nspan == 0:
                lspan_ref = self.model.aefact.slice_card_by_id(lspan)
                yi = lspan_ref.fractions
                nspan = len(yi) - 1
            else:
                yi = np.linspace(0., 1., nspan + 1)

            if nchord < 1 or nspan < 1:
                msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                    eid, nchord, nspan, lchord, lspan)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    # def xy(self):
    #     """
    #     Returns
    #     -------
    #     x : (nchord,) ndarray
    #         The percentage x location in the chord-wise direction of each panel
    #     y : (nspan,) ndarray
    #         The percentage y location in the span-wise direction of each panel
    #
    #     """
    #     #xy = []
    #     x = []
    #     y = []
    #     for eid, nchord, lchord, nspan, lspan in zip(self.element_id, self.nchord, self.lchord, self.nspan, self.lspan):
    #         if nchord == 0:
    #             xi = self.lchord_ref.fractions
    #             nchord = len(x) - 1
    #         else:
    #             xi = np.linspace(0., 1., nchord + 1)
    #
    #         if nspan == 0:
    #             yi = self.lspan_ref.fractions
    #             nspan = len(y) - 1
    #         else:
    #             yi = np.linspace(0., 1., nspan + 1)
    #
    #         if nchord < 1 or nspan < 1:
    #             msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
    #                 eid, nchord, nspan, lchord, lspan)
    #             raise RuntimeError(msg)
    #         x.append(xi)
    #         y.append(yi)
    #     return x, y

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        nspan_ = array_default_int(self.nspan, default=0, size=size)
        lspan_ = array_default_int(self.lspan, default=0, size=size)
        nchord_ = array_default_int(self.nchord, default=0, size=size)
        ztaic_ = array_default_int(self.ztaic, default=0, size=size)
        p_airfoil_ = array_default_int(self.p_airfoil, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        assert self.p4.shape[1] == 3, self.p4.shape
        p1_ = self.p1.tolist()
        p4_ = self.p4.tolist()
        for eid, label, p1, x12, p4, x43, cp, nspan, lspan, nchord, ztaic, p_airfoil, in zip(
                element_id, self.label, p1_, self.x12, p4_, self.x43, cp_,
                nspan_, lspan_, nchord_, ztaic_, p_airfoil_):
            list_fields = [
                'CAERO7', eid, label, cp, nspan, nchord, lspan, ztaic, p_airfoil,] + \
                p1 + [x12, None, None, None, None] + \
                p4 + [x43, None, None, None, None]
            bdf_file.write(print_card(list_fields))
        return
