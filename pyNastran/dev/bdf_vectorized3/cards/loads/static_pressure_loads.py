from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8 # , print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_string_or_blank, fields)
from pyNastran.utils.numpy_utils import (
    integer_types,
)

from pyNastran.dev.bdf_vectorized3.cards.elements.solid import (
    CTETRA, CTETRA_FACE_MAPPER,
    CPENTA, CPENTA_FACE_MAPPER,
    CHEXA, CHEXA_FACE_MAPPER,
    CPYRAM
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    hslice_by_idim, make_idim, searchsorted_filter, get_print_card_8_16, parse_load_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_float, array_str, array_default_int, array_default_float
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from .static_loads import Load

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf.types import TextIOLike
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from .dynamic_loads import LOADSET


class PLOAD(Load):
    """
    Static Pressure Load

    Defines a uniform static pressure load on a triangular or quadrilateral surface
    comprised of surface elements and/or the faces of solid elements.

    +-------+-----+------+----+----+----+----+
    |   1   |  2  |  3   | 4  | 5  | 6  | 7  |
    +=======+=====+======+====+====+====+====+
    | PLOAD | SID |  P   | G1 | G2 | G3 | G4 |
    +-------+-----+------+----+----+----+----+
    | PLOAD |  1  | -4.0 | 16 | 32 | 11 |    |
    +-------+-----+------+----+----+----+----+

    """
    _id_name = 'load_id'

    def add(self, sid: int, pressure: float, nodes: list[int],
            comment: str='') -> int:
        """
        Creates a PLOAD card, which defines a uniform pressure load on a
        shell/solid face or arbitrarily defined quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply
        nodes : list[int]
            The nodes that are used to define the normal are defined
            using the same method as the CTRIA3/CQUAD4 normal.
            n = 3 or 4
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, pressure, nodes, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'pressure')
        nodes = [integer(card, 3, 'n1'),
                 integer(card, 4, 'n2'),
                 integer(card, 5, 'n3'),
                 integer_or_blank(card, 6, 'n4', default=0)]
        assert len(card) <= 7, f'len(PLOAD card) = {len(card):d}\ncard={card}'
        self.cards.append((sid, pressure, nodes, comment))
        self.n += 1
        return self.n

    @Load.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        load_id = np.zeros(ncards, dtype='int32')
        pressure = np.zeros(ncards, dtype='float64')
        node_id = np.zeros((ncards, 4), dtype='int32')

        for icard, card in enumerate(self.cards):
            (sid, pressurei, nodesi, comment) = card
            load_id[icard] = sid
            pressure[icard] = pressurei
            if len(nodesi) == 3:
                node_id[icard, :-1] = nodesi
            else:
                node_id[icard, :] = nodesi
        self._save(load_id, pressure, node_id)
        self.cards = []

    def _save(self, load_id, pressure, node_id):
        if len(self.load_id) != 0:
            raise NotImplementedError()
        nloads = len(load_id)
        self.load_id = load_id
        self.pressure = pressure
        self.node_id = node_id
        self.n = nloads

    def __apply_slice__(self, load: PLOAD, i: np.ndarray) -> None:  # ignore[override]
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.pressure = self.pressure[i]
        load.node_id = self.node_id[i, :]

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        node_id = self.node_id.flatten()
        node_id = node_id[node_id != 0]
        geom_check(self,
                   missing,
                   node=(nid, node_id),)

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        #print_card = get_print_card_8_16(size)

        load_ids = array_default_int(self.load_id, size=size)
        node_ids = array_default_int(self.node_id, default=0, size=size)
        pressures = array_float(self.pressure, size=8, is_double=False)
        for load_id, pressure, nodes in zip(load_ids, pressures, node_ids):
            #list_fields = ['PLOAD', load_id, pressure] + nodes
            msg = 'PLOAD   %8s%8s%8s%8s%8s%8s\n' % (load_id, pressure,
                                                    nodes[0], nodes[1], nodes[2], nodes[3])
            bdf_file.write(msg)
        return

    def sum_forces_moments(self):
        nloads = len(self.load_id)
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()
        inode = np.searchsorted(grid.node_id, self.node_id)

        is_tri = (self.node_id[:, 3] == 0)
        is_quad = ~is_tri

        in1 = inode[:, 0]
        in2 = inode[:, 1]
        in3 = inode[:, 2]
        xyz1 = xyz_cid0[in1, :]
        xyz2 = xyz_cid0[in2, :]
        xyz3 = xyz_cid0[in3, :]

        area_normal = np.full((nloads, 3), np.nan, dtype='float64')
        centroid = np.full((nloads, 3), np.nan, dtype='float64')
        if np.any(is_tri):
            area_normali = np.cross(xyz2[is_tri, :]-xyz1[is_tri, :], xyz3[is_tri, :]-xyz1[is_tri, :])
            centroidi = (xyz1[is_tri, :] + xyz2[is_tri, :] + xyz3[is_tri, :]) / 3.
            area_normal[is_tri] = area_normali
            centroid[is_tri] = centroidi

        if np.any(is_quad):
            in4 = inode[:, 3]
            xyz4 = xyz_cid0[in4, :]
            area_normali = np.cross(xyz2[is_quad, :]-xyz1[is_quad, :], xyz3[is_quad, :]-xyz1[is_quad, :])
            centroidi = (xyz1[is_quad, :] + xyz2[is_quad, :] + xyz3[is_quad, :] + xyz4[is_quad, :]) / 4.
            area_normal[is_quad] = area_normali
            centroid[is_quad] = centroidi

        force = self.pressure[:, np.newaxis] * area_normal
        moment = np.cross(centroid, force)
        force_moment = np.hstack([force, moment])
        assert force_moment.shape == (nloads, 6)
        return force_moment


class PLOAD1(Load):
    """
    Applied Load on CBAR, CBEAM or CBEND Elements

    Defines concentrated, uniformly distributed, or linearly distributed
    applied loads to the CBAR or CBEAM elements at user-chosen points
    along the axis. For the CBEND element, only distributed loads over
    an entire length may be defined.

    +--------+-----+------+------+-------+-----+-------+-----+-------+
    |   1    |  2  |  3   |  4   |   5   |  6  |   7   |  8  |   9   |
    +========+=====+======+======+=======+=====+=======+=====+=======+
    | PLOAD1 | SID | EID  | TYPE | SCALE | X1  |  P1   |  X2 |  P2   |
    +--------+-----+------+------+-------+-----+-------+-----+-------+
    | PLOAD1 | 25  | 1065 |  MY  | FRPR  | 0.2 | 2.5E3 | 0.8 | 3.5E3 |
    +--------+-----+------+------+-------+-----+-------+-----+-------+

    """
    valid_types = ['FX', 'FY', 'FZ', 'FXE', 'FYE', 'FZE',
                   'MX', 'MY', 'MZ', 'MXE', 'MYE', 'MZE']

    # LE: length-based; FR: fractional; PR:projected
    valid_scales = ['LE', 'FR', 'LEPR', 'FRPR']

    _id_name = 'load_id'
    def add(self, sid: int, eid: int, load_type: str, scale: float,
            x1: float, p1: float,
            x2: Optional[float]=None, p2: Optional[float]=None,
            comment: str='') -> PLOAD1:
        """
        Creates a PLOAD1 card, which may be applied to a CBAR/CBEAM

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element to apply the load to
        load_type : str
            type of load that's applied
            valid_types = {FX, FY, FZ, FXE, FYE, FZE,
                           MX, MY, MZ, MXE, MYE, MZE}
        scale : str
            Determines scale factor for X1, X2.
            {LE, FR, LEPR, FRPR}
        x1 / x2 : float / float
            the starting/end position for the load application
            the default for x2 is x1
        p1 / p2 : float / float
            the magnitude of the load at x1 and x2
            the default for p2 is p1
        comment : str; default=''
            a comment for the card

        Point Load       : x1 == x2
        Distributed Load : x1 != x2

        """
        x2 = x2 if x2 is not None else x1
        p2 = p2 if p2 is not None else p1
        self.cards.append((sid, eid, load_type, scale, [x1, x2], [p1, p2], comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        load_type = string(card, 3, 'Type ("%s")' % '",  "'.join(self.valid_types))
        scale = string(card, 4, 'scale ("%s")' % '", "'.join(self.valid_scales))
        x1 = double(card, 5, 'x1')
        p1 = double(card, 6, 'p1')
        x2 = double_or_blank(card, 7, 'x2', default=x1)
        p2 = double_or_blank(card, 8, 'p2', default=p1)
        assert len(card) <= 9, f'len(PLOAD1 card) = {len(card):d}\ncard={card}'
        self.cards.append((sid, eid, load_type, scale, [x1, x2], [p1, p2], comment))
        self.n += 1
        return self.n

    @Load.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        load_id = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype='int32')
        load_type = np.zeros(ncards, dtype='|U8')
        scale = np.zeros(ncards, dtype='|U8')
        x = np.zeros((ncards, 2), dtype='float64')
        pressure = np.zeros((ncards, 2), dtype='float64')

        for icard, card in enumerate(self.cards):
            (sid, eid, load_typei, scalei, x12, pressurei, comment) = card
            load_id[icard] = sid
            element_id[icard] = eid
            load_type[icard] = load_typei
            scale[icard] = scalei
            x[icard] = x12
            pressure[icard, :] = pressurei
        self._save(load_id, element_id, load_type, scale, x, pressure)
        self.cards = []

    def _save(self, load_id, element_id, load_type, scale, x, pressure):
        if len(self.load_id) != 0:
            raise NotImplementedError()
        nloads = len(load_id)
        self.load_id = load_id
        self.element_id = element_id
        self.load_type = load_type
        self.scale = scale
        self.x = x
        self.pressure = pressure
        self.n = nloads


    def __apply_slice__(self, load: PLOAD1, i: np.ndarray) -> None:  # ignore[override]
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.element_id = self.element_id[i]
        load.load_type = self.load_type[i]
        load.scale = self.scale[i]
        load.x = self.x[i, :]
        load.pressure = self.pressure[i, :]

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> str:
        print_card = get_print_card_8_16(size)
        load_ids = array_str(self.load_id, size=size)
        element_ids = array_str(self.element_id, size=size)
        xs = array_float(self.x, size=size, is_double=is_double)
        pressures = array_float(self.pressure, size=size, is_double=is_double)
        for sid, eid, load_type, scale, x, pressure in zip(load_ids, element_ids, self.load_type,
                                                           self.scale, xs, pressures):
            x1, x2 = x
            p1, p2 = pressure
            list_fields = ['PLOAD1', sid, eid, load_type, scale,
                           x1, p1, x2, p2]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def line_elements(self) -> list[Any]:
        model = self.model
        return [card for card in [model.cbar, model.cbeam] if card.n]

    def geom_check(self, missing: dict[str, np.ndarray]):
        eids = hstack_msg([elem.element_id for elem in self.line_elements],
                          msg=f'no bar/beam elements for {self.type}')
        eids.sort()
        geom_check(self,
                   missing,
                   element_id=(eids, self.element_id))

    def sum_forces_moments(self):
        log = self.model.log
        load_ids = np.unique(self.load_id)
        log.warning(f'skipping PLOAD1 for load_id={load_ids}')
        nloads = len(self.load_id)
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        #load.load_id = self.load_id[i]
        #load.element_id = self.element_id[i]
        #load.scale = self.scale[i]
        #load.x = self.x[i, :]
        #load.pressure = self.pressure[i, :]
        scale = self.scale

        # could do proper integration
        #
        # ^
        # |
        # |     p1
        # |     *
        # |     | \
        # |     |  \
        # |     |   \
        # |     +----+ p2
        # |     |    |      phigh-plow = p1 - p2
        # |     |    |
        # +-----+----+------->
        #       x1   x2
        #
        x = self.x
        pressure = self.pressure.copy()
        load_type = self.load_type

        x1 = self.x[:, 0]
        x2 = self.x[:, 1]

        # basic frame
        ifx = np.where(load_type == 'FX')[0]
        ify = np.where(load_type == 'FY')[0]
        ifz = np.where(load_type == 'FZ')[0]
        imx = np.where(load_type == 'MX')[0]
        imy = np.where(load_type == 'MY')[0]
        imz = np.where(load_type == 'MZ')[0]

        is_basic = np.zeros(nloads, dtype='bool')
        is_basic_force = np.zeros(nloads, dtype='bool')
        is_basic_force[ifx] = True
        is_basic_force[ify] = True
        is_basic_force[ifz] = True
        is_basic_moment = np.zeros(nloads, dtype='bool')
        is_basic_moment[imx] = True
        is_basic_moment[imy] = True
        is_basic_moment[imz] = True
        is_basic[is_basic_force] = True
        is_basic[is_basic_moment] = True

        #is_basic_force = np.hstack([ifx, ify, ifz])
        #is_basic_moment = np.hstack([imx, imy, imz])
        #is_basic = np.hstack([is_basic_force, is_basic_moment])

        vector = np.zeros((nloads, 3), dtype='float64')
        if np.any(ifx):
            vector[ifx, :] = [1., 0., 0.]
        if np.any(ify):
            vector[ify, :] = [0., 1., 0.]
        if np.any(ifz):
            vector[ifz, :] = [0., 0., 1.]

        if np.any(imx):
            vector[imx, :] = [1., 0., 0.]
        if np.any(imy):
            vector[imy, :] = [0., 1., 0.]
        if np.any(imy):
            vector[imy, :] = [0., 0., 1.]

        # element coordinate system
        ifxe = np.where(load_type == 'FXE')[0]
        ifye = np.where(load_type == 'FYE')[0]
        ifze = np.where(load_type == 'FZE')[0]
        imxe = np.where(load_type == 'MXE')[0]
        imye = np.where(load_type == 'MYE')[0]
        imze = np.where(load_type == 'MZE')[0]
        is_element = np.zeros(nloads, dtype='bool')
        is_element_force = np.zeros(nloads, dtype='bool')
        is_element_force[ifxe] = True
        is_element_force[ifye] = True
        is_element_force[ifze] = True
        is_element_moment = np.zeros(nloads, dtype='bool')
        is_element_moment[imxe] = True
        is_element_moment[imye] = True
        is_element_moment[imze] = True
        is_element[is_element_force] = True
        is_element[is_element_moment] = True
        # -----------------------------------------------------------------------------

        #force_moment = np.full((nloads, 6), np.nan, dtype='float64')
        #force = force_moment[:, 3:]
        #moment = force_moment[:, :3]

        length = np.full(nloads, np.nan, dtype='float64')
        centroid = np.full((nloads, 3), np.nan, dtype='float64')
        #load_direction = np.full((nloads, 3), np.nan, dtype='float64')
        #x_fractional_local = np.full((nloads, 2), np.nan, dtype='float64')
        ivector = np.full((nloads, 3), np.nan, dtype='float64')
        is_length = (self.scale == 'LE')
        is_fractional = (self.scale == 'FR')
        is_length_projected = (self.scale == 'LEPR')
        is_fractional_projected = (self.scale == 'FRPR')

        is_concentrated = (self.x[:, 0] == self.x[:, 1])
        is_distributed = ~is_concentrated
        for card in self.line_elements:
            i_lookup, i_all = searchsorted_filter(card.element_id, self.element_id, msg=f'{card.type} eids', debug=True)
            if len(i_lookup) == 0:
                continue

            # we're at least using some loads
            #pressure[i_lookup] = load_pressure[i_lookup, :nnids]
            ivectori, lengthi = card.line_vector_length()
            centroidi = card.centroid()
            #area[i_lookup] = areai[i_all]
            length[i_lookup] = lengthi[i_all]
            ivector[i_lookup, :] = ivectori[i_all, :]
            centroid[i_lookup, :] = centroidi[i_all, :]


        # LE: length-based; FR: fractional; PR:projected
        #valid_scales = ['LE', 'FR', 'LEPR', 'FRPR']

        # concentrated loads
        #if np.any(is_concentrated):
            #concentrated_load
        #if np.any(is_distributed):
            #distributed_load

        # 6. If SCALE = 'LE' (length), the xi values are actual distances along the element axis,
        # and, if X2 != X1, then Pi are load intensities per unit length of the element.
        # 7. If SCALE = 'FR' (fractional), the xi values are ratios of the distance along the
        # axis to the total length, and (if X2 != X1) Pi are load intensities per unit length of
        # the element.
        # 8. If SCALE = 'LEPR' (length projected), the xi values are actual distances along
        # the element axis, and (if X2 != X1) the distributed load is input in terms of the
        # projected length of the element.

        #normal[i_lookup, :] = normali[i_all, :]
        #x = 1


        # convert to per length; get rid of fractional
        # just multiply by length
        if np.any(is_fractional):
            pressure[is_fractional] *= length[is_fractional, np.newaxis]
            is_length = is_length | is_fractional
        if np.any(is_fractional_projected):
            pressure[is_fractional_projected] *= length[is_fractional_projected, np.newaxis]
            is_length_projected = is_length_projected | is_fractional_projected

        if 0:
            p1 = self.pressure[:, 0]
            p2 = self.pressure[:, 1]
            ptri = p2 - p1
            psquare = p2

            #unit_load = psquare * (x2 - x1) + ptri * (x2 - x1) / 2.
            #x_avg = self.x.mean(axis=1)
            #x_avg = (x2 * (x2 - x1) + x1 * (x2 - x1) / 2. ) # / x + avg
            p_avg = (2 * p1 + p2) / 3.
            x_avg = (x1 + 2 * x2) / 3.
            # p_avg = self.pressure.mean(axis=1)
            assert len(x_avg) == nloads

        p1 = pressure[:, 0]
        p2 = pressure[:, 1]
        #ptri = p2 - p1
        #psquare = p2

        #unit_load = psquare * (x2 - x1) + ptri * (x2 - x1) / 2.
        #x_avg = self.x.mean(axis=1)
        #x_avg = (x2 * (x2 - x1) + x1 * (x2 - x1) / 2. ) # / x + avg
        # simple centroid formula with:
        # (x1,p1),(x2,p1),(x2,p2)
        p_avg = (2 * p1 + p2) / 3.
        x_avg = (x1 + 2 * x2) / 3.

        mag = p_avg * x_avg
        global_force = np.full((nloads, 3), np.nan, dtype='float64')
        global_moment = np.full((nloads, 3), np.nan, dtype='float64')
        if np.any(is_length):
            # just take the raw value
            iforce = is_basic_force & is_length
            imoment = is_basic_moment & is_length
            global_force[iforce, :] = vector[iforce, :]
            global_moment[imoment, :] = vector[imoment, :]

        if np.any(is_length_projected):
            # doing the example in the QRG...
            # a 30-60-90 triangle has sides of length:
            #  - opposite to 30: 1
            #  - opposite to 60: sqrt(3)
            #  - opposite to 90: 2
            # the projection is calculated as:
            #   Fy = vector of <0., -1., 0.>
            #   ivector = <1/sqrt(3), -1, 0.>
            #   cos(30) = sqrt(3) / 2
            # proj of v onto u
            #   proj = (u @ v) / v^2 * v
            # proj of Fy onto i
            #   proj = (i @ Fy) / Fy^2 * Fy
            #
            # not 100%...
            iforce = is_basic_force & is_length_projected
            imoment = is_basic_moment & is_length_projected
            if np.any(iforce):
                forcei = vector[iforce, :] @ ivector[iforce, :].T
                assert len(forcei) == iforce.sum()
                global_force[iforce, :] = forcei
            if np.any(imoment):
                global_moment[imoment, :] = vector[imoment, :] @ ivector[imoment, :].T

        global_force[is_basic_force, :] *= mag[is_basic_force, np.newaxis]
        global_force[is_basic_moment, :] = 0.

        # should be the local coordinatate, aka xavg
        global_moment[is_basic_force, :] = np.cross(centroid[is_basic_force, :], global_force[is_basic_force, :])
        global_moment[is_basic_moment, :] *= mag[is_basic_moment, np.newaxis]
        is_basic
        is_element
        #inode = np.searchsorted(grid.node_id, self.node_id)
        force_moment = np.hstack([global_force, global_moment])
        return force_moment


class PLOAD2(Load):
    """
    +--------+-----+------+------+------+------+------+------+------+
    |    1   |   2 |  3   |  4   |   5  |   6  |   7  |   8  |   9  |
    +========+=====+======+======+======+=============+======+======+
    | PLOAD2 | SID |  P   | EID1 | EID2 | EID3 | EID4 | EID5 | EID6 |
    +--------+-----+------+------+------+------+------+------+------+
    | PLOAD2 | 21  | -3.6 |  4   |  16  |  2   |      |      |      |
    +--------+-----+------+------+------+------+------+------+------+
    | PLOAD2 | SID |  P   | EID1 | THRU | EID2 |      |      |      |
    +--------+-----+------+------+------+------+------+------+------+

    """
    _id_name = 'load_id'
    #def __init__(self, model: BDF):
        #super().__init__(model)
        #self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> PLOAD2:
        #if len(i) == self.n:
            #return self
        #load = PLOAD2(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def add(self, sid: int, pressure: float, eids: list[int],
            comment: str='') -> int:
        """
        Creates a PLOAD2 card, which defines an applied load normal
        to the quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply to the elements
        eids : list[int]
            the elements to apply pressure to
            n < 6 or a continouus monotonic list of elements (e.g., [1, 2, ..., 1000])
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, pressure, eids, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'p')

        if integer_string_or_blank(card, 4, 'THRU') == 'THRU':
            e1 = integer(card, 3, 'Element1')
            e2 = integer(card, 5, 'Element1')
            eids = [i for i in range(e1, e2 + 1)]
            assert len(card) == 6, f'len(PLOAD2 card) = {len(card):d}\ncard={card}'
        else:
            eids = fields(integer, card, 'eid', i=3, j=len(card))
            assert len(eids) <= 6, f'A maximum of 6 eids may be on the PLOAD2; n={len(eids)}\ncard={card}'
        self.cards.append((sid, pressure, eids, comment))
        self.n += 1
        return self.n

    @Load.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        load_id = np.zeros(ncards, dtype='int32')
        pressure = np.zeros(ncards, dtype='float64')
        nelement = np.zeros(ncards, dtype='int32')

        element_ids_ = []
        for icard, card in enumerate(self.cards):
            (sid, pressurei, eids, comment) = card
            load_id[icard] = sid
            pressure[icard] = pressurei
            nelement[icard] = len(eids)
            element_ids_.extend(eids)
        element_ids = np.array(element_ids_, dtype='int32')
        self. _save(load_id, pressure, element_ids, nelement)
        self.cards = []

    def _save(self, load_id, pressure, element_ids, nelement):
        if len(self.load_id):
            load_id = np.hstack([self.load_id, load_id])
            pressure = np.hstack([self.pressure, pressure])
            element_ids = np.hstack([self.element_ids, element_ids])
            nelement = np.hstack([self.nelement, nelement])
        nloads = len(load_id)
        self.load_id = load_id
        self.pressure = pressure
        self.element_ids = element_ids
        self.nelement = nelement
        self.n = nloads

    def __apply_slice__(self, load: PLOAD2, i: np.ndarray) -> None:  # ignore[override]
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.pressure = self.pressure[i]
        load.nelement = self.nelement[i]

        idim = self.idim
        load.element_ids = hslice_by_idim(i, idim, self.element_ids)

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def shell_elements(self) -> list[Any]:
        model = self.model
        elements = [card for card in [model.ctria3, model.ctria6, model.ctriar,
                                      model.cquad4, model.cquad8, model.cquadr]
                if card.n]
        return elements

    def geom_check(self, missing: dict[str, np.ndarray]):
        eids = hstack_msg([elem.element_id for elem in self.shell_elements],
                          msg=f'no shell elements for {self.type}')
        eids.sort()
        geom_check(self,
                   missing,
                   element_id=(eids, self.element_ids))

    @property
    def is_small_field(self):
        return max(self.load_id.max(), self.element_ids.max()) < 99_999_999

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        #print_card = get_print_card_8_16(size)
        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        for load_id, pressure, idim in zip(self.load_id, self.pressure, self.idim):
            idim0, idim1 = idim
            eids = self.element_ids[idim0:idim1]
            list_fields = ['PLOAD2', load_id, pressure]
            eids = self.element_ids[idim0:idim1].tolist()
            if len(eids) <= 6:
                list_fields += eids
            else:
                eids.sort()
                delta_eid = eids[-1] - eids[0] + 1
                if delta_eid != len(eids):
                    msg = 'eids=%s len(eids)=%s delta_eid=%s must be continuous' % (
                        eids, len(eids), delta_eid)
                    raise RuntimeError(msg)
                #list_fields += eids
                list_fields += [eids[0], 'THRU', eids[-1]]
            bdf_file.write(print_card(list_fields))
        return

    def sum_forces_moments(self) -> np.ndarray:
        #log = self.model.log
        load_ids = np.unique(self.load_id)
        nloads = len(self.load_id)
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        ithru_load, ithru_result = _get_ithru(self.nelement)
        pressure = self.pressure[ithru_result]

        nelement = self.nelement.sum()
        area = np.full(nelement, np.nan, dtype='float64')
        centroid = np.full((nelement, 3), np.nan, dtype='float64')
        normal = np.full((nelement, 3), np.nan, dtype='float64')
        #pressure = np.full(nelement, np.nan, dtype='float64')
        for card in self.shell_elements:
            i_lookup, i_all = searchsorted_filter(card.element_id, self.element_ids, msg=f'{card.type} eids', debug=True)
            if len(i_lookup) == 0:
                continue

            # we're at least using some loads
            nnids = card.base_nodes.shape[1]
            #pressure[i_lookup] = load_pressure[i_lookup, :nnids]
            areai, centroidi, normali = card.area_centroid_normal()
            area[i_lookup] = areai[i_all]
            centroid[i_lookup, :] = centroidi[i_all, :]
            normal[i_lookup, :] = normali[i_all, :]

        #inode = np.searchsorted(grid.node_id, self.node_id)

        #self.load_id = np.zeros(ncards, dtype='int32')
        #self.pressure = np.zeros(ncards, dtype='float64')
        #self.nelement = np.zeros(ncards, dtype='int32')

        pa =pressure * area
        force = pa[:, np.newaxis] * normal
        moment = np.cross(centroid, force)
        force_moment = np.hstack([force, moment])
        assert force_moment.shape == (nelement, 6)
        return force_moment


class PLOAD4(Load):
    _id_name = 'load_id'
    #def slice_card_by_index(self, i: np.ndarray) -> PLOAD4:
        #load = PLOAD4(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def add(self,
            sid: int,
            eids: list[int],
            pressures: list[float],
            g1=-1, g34=-1, cid: int=0,
            nvector=None,
            surf_or_line: str='SURF', line_load_dir: str='NORM',
            comment: str='') -> int:
        """
        Creates a PLOAD4 card

        Parameters
        ----------
        sid : int
            the load id
        eids : list[int, ...]
            shells : the range of element ids; must be sequential
            solids : must be length 1
        pressures : list[float, float, float, float]
            tri : must be length 4 (the last value should be the same as the 0th value)
            quad : must be length 4
        g1 : int/None
            only used for solid elements
        g34 : int / None
            only used for solid elements
        cid : int; default=0
            the coordinate system for ???
        nvector : (3, ) float ndarray
           blank : load acts normal to the face
           the local pressure vector
        surf_or_line : str; default='SURF'
           SURF : surface load
           LINE : line load (only defined for QUADR, TRIAR)
           not supported
        line_load_dir : str; default='NORM'
           direction of the line load (see surf_or_line); {X, Y, Z, TANG, NORM}
           not supported
        comment : str; default=''
            a comment for the card

        TODO: fix the way "pressures" works

        """
        if nvector is None:
            nvector = [0., 0., 0.]
        if isinstance(eids, integer_types):
            eids = [eids]
        eid = eids[0]

        g1 = g1 if g1 is not None else -1
        g34 = g34 if g34 is not None else -1
        self.cards.append((sid, eid, pressures, eids, g1, g34,
                           cid, nvector, surf_or_line, line_load_dir, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        p1 = double_or_blank(card, 3, 'p1', default=0.0)
        pressures = [
            p1,
            double_or_blank(card, 4, 'p2', default=np.nan),
            double_or_blank(card, 5, 'p3', default=np.nan),
            double_or_blank(card, 6, 'p4', default=np.nan),
        ]

        eids = [eid]
        g1_thru = integer_string_or_blank(card, 7, 'g1/THRU')
        if g1_thru == 'THRU' and integer_or_blank(card, 8, 'eid2'):
            # alternate form
            eid2 = integer(card, 8, 'eid2')
            if eid2:
                eids = list(np.unique(
                    expand_thru([eid, 'THRU', eid2], set_fields=False, sort_fields=False)
                ))
            g1 = -1
            g34 = -1
        else:
            # standard form
            eids = [eid]
            g1 = integer_or_blank(card, 7, 'g1', default=-1)
            g34 = integer_or_blank(card, 8, 'g34', default=-1)

        # If both (CID, N1, n2, N3) and LDIR are blank, then the default is LDIR=NORM.
        cid = integer_or_blank(card, 9, 'cid', default=-1)
        n1 = double_or_blank(card, 10, 'N1', default=0.0)
        n2 = double_or_blank(card, 11, 'N2', default=0.0)
        n3 = double_or_blank(card, 12, 'N3', default=0.0)

        surf_or_line = string_or_blank(card, 13, 'sorl', default='SURF')
        line_load_dir = string_or_blank(card, 14, 'ldir', default='NORM')
        assert len(card) <= 15, f'len(PLOAD4 card) = {len(card):d}\ncard={card}'
        assert isinstance(g1, integer_types), g1
        assert isinstance(g34, integer_types), g34
        self.cards.append((sid, eid, pressures, eids, g1, g34,
                           cid, [n1, n2, n3], surf_or_line, line_load_dir, comment))
        self.n += 1
        return self.n

    @Load.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        load_id = np.zeros(ncards, dtype='int32')
        coord_id = np.full(ncards, -1, dtype='int32')
        pressure = np.zeros((ncards, 4), dtype='float64')
        nodes_g1_g34 = np.full((ncards, 2), -1, dtype='int32')

        surf_or_line = np.full(ncards, '', dtype='|U8')
        line_load_dir = np.full(ncards, '', dtype='|U8')
        nvector = np.zeros((ncards, 3), dtype='float64')
        nelement = np.zeros(ncards, dtype='int32')
        element_ids = []
        for icard, card in enumerate(self.cards):
            (sid, eid, pressures, eids, g1, g34,
             cid, n123, surf_or_linei, line_load_diri, comment) = card
            load_id[icard] = sid
            pressure[icard, :] = pressures
            if cid is not None:
                coord_id[icard] = cid
            #element_id[icard] = eid
            nvector[icard, :] = n123
            surf_or_line[icard] = surf_or_linei
            line_load_dir[icard] = line_load_diri
            nodes_g1_g34[icard, :] = [g1, g34]
            nelement[icard] = len(eids)
            element_ids.extend(eids)

        element_ids_array = np.array(element_ids, dtype='int32')
        self._save(load_id, element_ids_array, coord_id, pressure, nodes_g1_g34,
                   surf_or_line, line_load_dir, nvector, nelement)
        self.cards = []
        self.nvector

    def _save(self, load_id, element_ids, coord_id, pressure, nodes_g1_g34,
              surf_or_line, line_load_dir, nvector, nelement):

        nloads = len(load_id)
        if surf_or_line is None:
            surf_or_line = np.full(nloads, 'SURF', dtype='|U8')
        if line_load_dir is None:
            line_load_dir = np.full(nloads, 'NOMR', dtype='|U8')

        if len(self.load_id) != 0:
            load_id = np.hstack([self.load_id, load_id])
            element_ids = np.hstack([self.element_ids, element_ids])
            coord_id = np.hstack([self.coord_id, coord_id])
            pressure = np.vstack([self.pressure, pressure])
            nodes_g1_g34 = np.vstack([self.nodes_g1_g34, nodes_g1_g34])
            surf_or_line = np.hstack([self.surf_or_line, surf_or_line])
            line_load_dir = np.hstack([self.line_load_dir, line_load_dir])
            nvector = np.vstack([self.nvector, nvector])
            nelement = np.hstack([self.nelement, nelement])

        nloads = len(load_id)
        self.load_id = load_id
        self.element_ids = element_ids
        self.coord_id = coord_id
        self.pressure = pressure
        assert pressure.shape == (nloads, 4), pressure
        self.nodes_g1_g34 = nodes_g1_g34

        self.surf_or_line = surf_or_line
        self.line_load_dir = line_load_dir
        assert nvector is not None
        self.nvector = nvector
        self.nelement = nelement
        self.n = nloads

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id

        shell_elements = self.shell_elements

        is_solid = self.is_solid
        is_shell = self.is_shell
        if np.any(is_shell):
            eids_shell = hstack_msg([elem.element_id for elem in shell_elements],
                                    msg=f'no shell elements for {self.type}')
            eids_shell.sort()
            shell_element_ids = np.hstack([self.element_ids[idim0:idim1]
                                           for is_solidi, (idim0, idim1) in zip(is_shell, self.ielement)])
            geom_check(self,
                       missing,
                       element_id=(eids_shell, shell_element_ids))

        if np.any(is_solid):
            solid_elements = self.solid_elements
            eids_solid = hstack_msg([elem.element_id for elem in solid_elements],
                                    msg=f'no solid elements for {self.type}')
            eids_solid.sort()
            solid_element_ids = np.hstack([self.element_ids[idim0:idim1]
                                           for is_solidi, (idim0, idim1) in zip(is_solid, self.ielement)])
            solid_nodes = self.nodes_g1_g34[is_solid, :]
            solid_element_ids = np.array([]),
            geom_check(self,
                       missing,
                       node=(nid, solid_nodes),
                       element_id=(eids_solid, solid_element_ids),)
        #pids.sort()
        # shells
        self.nodes_g1_g34
        shell_elements

    def __apply_slice__(self, load: PLOAD4, i: np.ndarray) -> None:  # ignore[override]
        nloads = len(self.load_id)
        assert self.nvector.shape == (nloads, 3), self.nvector.shape
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.coord_id = self.coord_id[i]
        load.pressure = self.pressure[i]
        load.nodes_g1_g34 = self.nodes_g1_g34[i, :]
        if len(i) == 0:
            load.element_ids = load.element_ids[i]
        else:
            load.element_ids = hslice_by_idim(i, self.ielement, self.element_ids)
        load.surf_or_line = self.surf_or_line[i]
        load.line_load_dir = self.line_load_dir[i]
        load.nvector = self.nvector[i, :]
        load.nelement = self.nelement[i]

        nloads = len(load.load_id)
        assert load.nvector.shape == (nloads, 3), load.nvector.shape

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> str:
        print_card = get_print_card_8_16(size)
        assert self.nvector is not None
        #nvectors = array_default_float(self.nvector)
        for load_id, cid, pressures, g1_g34, nvector, \
            surf_or_line, line_load_dir, ielement in zip(self.load_id, self.coord_id, self.pressure,
                                                         self.nodes_g1_g34, self.nvector,
                                                         self.surf_or_line, self.line_load_dir,
                                                         self.ielement):
            idim0, idim1 = ielement
            eids = self.element_ids[idim0:idim1]
            eid = eids[0]
            p1 = pressures[0]
            p2 = set_blank_if_default(pressures[1], p1)
            p3 = set_blank_if_default(pressures[2], p1)
            p4 = set_blank_if_default(pressures[3], p1)
            list_fields = ['PLOAD4', load_id, eid, pressures[0], p2, p3, p4]

            g1, g34 = g1_g34
            if g1 != -1:
                # is it a SOLID element
                list_fields += [g1, g34]
            else:
                if len(eids) > 1:
                    try:
                        list_fields.append('THRU')
                        eidi = eids[-1]
                    except Exception:
                        print("g1  = %s" % g1)
                        print("g34 = %s" % g34)
                        print("self.eids = %s" % eids)
                        raise
                    list_fields.append(eidi)
                else:
                    list_fields += [None, None]

            #+--------+-----+-----+----+----+------+------+------+-------+
            #|   1    |  2  |  3  |  4 |  5 |  6   |   7  |   8  |   9   |
            #+========+=====+=====+====+====+======+======+======+=======+
            #| PLOAD4 | SID | EID | P1 | P2 |  P3  |  P4  | THRU | EID2  |
            #+--------+-----+-----+----+----+------+------+------+-------+
            #|        | CID | N1  | N2 | N3 | SORL | LDIR |      |       |
            #+--------+-----+-----+----+----+------+------+------+-------+

            if cid != -1 or np.abs(nvector).max() > 0.:
                n1, n2, n3 = nvector
                list_fields.append(cid)
                list_fields += [n1, n2, n3]
                surf_or_line = set_blank_if_default(surf_or_line, 'SURF')
                line_load_dir = set_blank_if_default(line_load_dir, 'NORM')
            else:
                list_fields += [None, None, None, None]
                surf_or_line = set_blank_if_default(surf_or_line, 'SURF')
                line_load_dir = set_blank_if_default(line_load_dir, 'NORM')
            list_fields.append(surf_or_line)
            if surf_or_line == 'LINE':
                list_fields.append(line_load_dir)
            bdf_file.write(print_card(list_fields))
        return

    @property
    def shell_elements(self) -> list[Any]:
        model = self.model
        elements = [card for card in [model.ctria3, model.ctria6, model.ctriar,
                                      model.cquad4, model.cquad8, model.cquadr]
                    if card.n]
        return elements
    @property
    def solid_elements(self) -> list[Any]:
        model = self.model
        elements = [card for card in [model.ctetra, model.cpenta, model.cpyram, model.chexa]
                    if card.n]
        return elements

    @property
    def applied_pressure(self) -> np.ndarray:
        pressure = self.pressure
        applied_pressure = pressure.copy()
        for j in range(1, 4):
            pressurej = pressure[:, j]
            inan = np.isnan(pressurej)
            applied_pressure[inan, j] = pressure[inan, 0]
        return applied_pressure

    @property
    def is_shell(self) -> np.ndarray:
        nloads = len(self.load_id)
        min_nid = self.nodes_g1_g34.min(axis=1)
        assert len(min_nid) == nloads, min_nid.shape
        is_shell = (min_nid == -1)
        #assert len(is_shell) == self.element_ids
        return is_shell

    @property
    def is_solid(self) -> np.ndarray:
        return ~self.is_shell

    def geom_check(self, missing: dict[str, np.ndarray]):
        available_eids = hstack_msg([elem.element_id for elem in self.shell_elements + self.solid_elements],
                                    msg=f'no shell/solid elements for {self.type}')
        available_eids.sort()
        geom_check(self,
                   missing,
                   element_id=(available_eids, self.element_ids))

        #available_shell_eids = hstack_msg([elem.element_id for elem in self.shell_elements],
                                          #msg=f'no shell elements for {self.type}')
        #available_solid_eids = hstack_msg([elem.element_id for elem in self.solid_elements],
                                          #msg=f'no solid elements for {self.type}')
        #available_shell_eids.sort()
        #available_solid_eids.sort()

        #shell_eids = self.element_ids[self.is_shell]
        #solid_eids = self.element_ids[self.is_solid]
        #geom_check(self,
                   #missing,
                   #element_id=(available_shell_eids, shell_eids))
        #geom_check(self,
                   #missing,
                   #element_id=(available_solid_eids, solid_eids))

    def area_centroid_normal_pressure(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """gets the area, centroid, normal, and pressure"""
        log = self.model.log

        xyz_cid0 = self.model.grid.xyz_cid0()
        nids = self.model.grid.node_id
        #self.load_id[icard] = sid
        #self.pressure[icard, :] = pressures
        #self.coord_id[icard] = cid
        #self.nvector[icard, :] = [n1, n2, n3]
        #self.surf_or_line[icard] = surf_or_line
        #self.line_load_dir[icard] = line_load_dir
        #self.nodes_g1_g34[icard, :] = [g1, g34]
        nloads = len(self.load_id)
        assert self.nvector.shape == (nloads, 3), self.nvector.shape

        #grid = self.model.grid
        #xyz_cid0 = grid.xyz_cid0()
        #nid = grid.node_id

        is_shell_ = self.is_shell
        is_solid_ = ~is_shell_

        # nvector defined in coord_id
        is_normal_ = self.nvector.max(axis=1) == 0.0
        assert len(is_normal_) == nloads

        nvector_norm_ = np.linalg.norm(self.nvector, axis=1)
        assert len(nvector_norm_) == nloads

        with np.errstate(divide='ignore', invalid='ignore'):
            nvector_ = self.nvector / nvector_norm_[:, np.newaxis]
        assert nvector_.shape == (nloads, 3)

        nelement = self.nelement.sum()
        area = np.full(nelement, np.nan, dtype='float64')
        centroid = np.full((nelement, 3), np.nan, dtype='float64')
        normal = np.full((nelement, 3), np.nan, dtype='float64')
        #pressure = np.full((nelement, 4), np.nan, dtype='float64')
        mean_pressure = np.full(nelement, np.nan, dtype='float64')

        ithru_load, ithru_result = _get_ithru(self.nelement)

        # transform these to the thru numbering system
        #load_id = self.load_id[ithru_result]
        #coord_id = self.coord_id[ithru_result]
        nvector = nvector_[ithru_result]
        #nvector_norm = nvector_norm_[ithru_result]
        is_normal = is_normal_[ithru_result]
        is_shell = is_shell_[ithru_result]
        is_solid = is_solid_[ithru_result]

        applied_pressure = self.applied_pressure
        #applied_pressure2 = applied_pressure.mean(axis=1)
        load_pressure = applied_pressure[ithru_result, :]

        if np.all(is_shell_) and 0:
            for card in self.shell_elements:
                # TODO: use is_shell
                i_lookup, i_all = searchsorted_filter(card.element_id, self.element_ids, msg=f'{card.type} eids', debug=True)
                if len(i_lookup) == 0:
                    continue

                # we're at least using some loads
                #print(card.type)
                nnids = card.base_nodes.shape[1]
                pressurei = load_pressure[i_lookup, :nnids].copy()
                remove_nan_pressure(pressurei)
                #pressure[i_lookup, :nnids] = pressurei
                mean_pressure[i_lookup] = pressurei.mean(axis=1)
                areai, centroidi, normali = card.area_centroid_normal()
                area[i_lookup] = areai[i_all]
                centroid[i_lookup, :] = centroidi[i_all, :]
                normal[i_lookup, :] = normali[i_all, :]

        elif np.all(is_solid) and 0:
            for card in self.solid_elements:
                # TODO: use is_shell
                i_lookup, i_all = searchsorted_filter(card.element_id, self.element_ids, msg=f'{card.type} eids', debug=True)
                if len(i_lookup) == 0:
                    continue

                load_pressurei = load_pressure[i_lookup]
                element_nodes = card.base_nodes[i_all, :]
                #print(nodes)
                #assert nodes.shape[0] == 1, nodes.shape

                g1 = self.nodes_g1_g34[i_lookup, 0]
                g34 = self.nodes_g1_g34[i_lookup, 1]
                if load_pressurei.shape[0] != len(g1):
                    msg = f'ng1 = {len(g1)}\nload_pressurei.shape={load_pressurei.shape}'
                    raise ValueError(msg)

                #print(self)
                areai, centroidi, normali, mean_pressurei = get_solid_face_area(
                    card.type, card,
                    nids, xyz_cid0, element_nodes,
                    g1, g34, load_pressurei)
                area[i_lookup] = areai
                mean_pressure[i_lookup] = mean_pressurei
                assert len(areai) == len(mean_pressurei)
                assert len(normali) == len(mean_pressurei)
                normal[i_lookup, :] = normali
                centroid[i_lookup, :] = centroidi

            if np.any(np.isnan(mean_pressure)):
                raise RuntimeError('there are nan pressures', mean_pressure)

        else:
            # mixed or maybe mixed
            #print(self)
            #mixed
            #print('is_shell =', is_shell)
            for card in self.shell_elements:
                # TODO: use is_shell
                self.model.log.debug(card.type)
                i_lookup, i_all = searchsorted_filter(card.element_id, self.element_ids, msg=f'{card.type} eids', debug=True)
                if len(i_lookup) == 0:
                    continue
                log.debug(f'ilookup = {i_lookup}')

                # we're at least using some loads
                nnids = card.base_nodes.shape[1]
                pressurei = load_pressure[i_lookup, :nnids].copy()
                remove_nan_pressure(pressurei)
                #pressure[i_lookup, :nnids] = pressurei
                mean_pressure[i_lookup] = pressurei.mean(axis=1)
                areai, centroidi, normali = card.area_centroid_normal()
                area[i_lookup] = areai[i_all]
                centroid[i_lookup, :] = centroidi[i_all, :]
                normal[i_lookup, :] = normali[i_all, :]

            for card in self.solid_elements:
                # TODO: use is_shell
                i_lookup, i_all = searchsorted_filter(card.element_id, self.element_ids, msg=f'{card.type} eids', debug=True)
                if len(i_lookup) == 0:
                    continue

                load_pressurei = load_pressure[i_lookup]
                element_nodes = card.base_nodes[i_all, :]
                #print(nodes)
                #assert nodes.shape[0] == 1, nodes.shape

                g1 = self.nodes_g1_g34[i_lookup, 0]
                g34 = self.nodes_g1_g34[i_lookup, 1]
                if load_pressurei.shape[0] != len(g1):
                    msg = f'ng1 = {len(g1)}\nload_pressurei.shape={load_pressurei.shape}'
                    raise ValueError(msg)

                #print(self)
                areai, centroidi, normali, mean_pressurei = get_solid_face_area(
                    card.type, card,
                    nids, xyz_cid0, element_nodes,
                    g1, g34, load_pressurei)
                area[i_lookup] = areai
                mean_pressure[i_lookup] = mean_pressurei
                assert len(areai) == len(mean_pressurei)
                assert len(normali) == len(mean_pressurei)
                normal[i_lookup, :] = normali
                centroid[i_lookup, :] = centroidi

        # apply nvector to the normal for places that the normal vector isn't used
        normal[~is_normal, :] = nvector[~is_normal, :]

        inan = np.isnan(mean_pressure)
        if np.any(inan):
            assert len(mean_pressure) == len(self.element_ids)
            eids_nan = self.element_ids[inan]
            raise RuntimeError(f'there are nan PLOAD4 pressures for eids={eids_nan}')
        assert len(area) == nelement, f'narea={len(area)}; nelement={nelement}'
        assert len(mean_pressure) == nelement, f'npressure={len(mean_pressure)}; nelement={nelement}'

        normal_norm = np.linalg.norm(normal)
        if np.any(np.isnan(normal_norm)):
            raise RuntimeError('there are nan normals', normal)


        centroid_norm = np.linalg.norm(centroid)
        if np.any(np.isnan(centroid_norm)):
            raise RuntimeError('there are nan centroid', centroid)
        assert centroid.shape == normal.shape

        return area, centroid, normal, mean_pressure

    def sum_forces_moments(self):
        area, centroid, normal, mean_pressure = self.area_centroid_normal_pressure()
        nelement = len(area)

        # global or default
        if self.coord_id.max() in {-1, 0} and self.coord_id.min() in {-1, 0}:
            # F = p*A
            pa = mean_pressure * area
            force = pa[:, np.newaxis] * normal
            moment = np.cross(centroid, force)
        else:
            ucoords = np.unique(self.coord_id)
            msg = f'coords={ucoords}'
            raise NotImplementedError(msg)

        #print('mean_pressure =', mean_pressure)
        #print('area =', area)
        #print('pa =', pa)
        #print('normal=', normal)
        force_moment = np.hstack([force, moment])
        assert force_moment.shape == (nelement, 6)
        #force_moment = np.zeros((nelement, 6), dtype='float64')
        #force_moment[:, :3] = force
        #force_moment[:, 3:] = moment
        return force_moment

def _get_ithru(nelement: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nelementi = nelement.sum()
    is_thru = (nelement.max() > 1)

    ithru_load = np.arange(nelementi, dtype='int32')
    if is_thru:
        ithru_result_ = []
        for i, nelement_ in enumerate(nelement):
            ithru_result_.extend([i]*nelement_)
        ithru_result = np.array(ithru_result_, dtype='int32')
        #del ithru_result_
    else:
        ithru_result = ithru_load
    return ithru_load, ithru_result

def remove_nan_pressure(pressure: np.ndarray):
    nnodes = pressure.shape[1]
    for i in range(1, nnodes):
        inan = np.isnan(pressure[:, i])
        pressure[inan, i] = pressure[inan, 0]

def _solid_mean_pressure(face, pressure):
    nnodes_on_face = len(face)
    pi = pressure[:, :nnodes_on_face]
    #remove_nan_pressure(pi)
    mean_pressures = pi.mean(axis=1)
    return mean_pressures

CTETRA_FACE_MAPPER = {
    #(3, 1, 0),
    #(1, 2):
    #(2, 3): (1, 2, 4),
    #(1, 3): (1, 2, 4),
    #(1, 4): ,
    #(3, 2, 1),
    #(0, 2, 3),
    (1, 0): (0, 2, 3),
}


def get_solid_face_area(element_type: str, element: Union[CTETRA, CHEXA, CPENTA, CPYRAM],
                        nids, xyz_cid0,
                        element_nodes,
                        g1: np.ndarray, g34: np.ndarray,
                        pressure: np.ndarray) -> list[Any]:
    # we're at least using some loads
    nnodes = element.base_nodes.shape[1]
    nloads = len(g1)

    area = np.full(nloads, np.nan, dtype=xyz_cid0.dtype)
    normal = np.full((nloads, 3), np.nan, dtype=xyz_cid0.dtype)
    centroid = np.full((nloads, 3), np.nan, dtype=xyz_cid0.dtype)

    iblank34 = (g34 == -1)
    i34 = ~iblank34
    all_ids = np.arange(nloads, dtype=element_nodes.dtype)
    ids_with_34 = all_ids[i34]
    ids_no_34 = all_ids[iblank34]
    log = element.model.log
    log.debug(f'g1 = {g1}')
    log.debug(f'g3/4 = {g34}')

    #ig1_list = []
    #ig34_list = []
    # both nodes are valid?
    #for enodes, g1i, g34i in zip(nodes[i34, :].tolist(), g1[i34], g34[i34]):
        #ig1i = enodes.index(g1i)
        #ig34i = enodes.index(g34i)
        #print(f'{element_type}: nnodes={nnodes}; ig1={ig1i}; ig3/4={ig34i}')
        ##card.map_node_index_to_face_number(ig1, ig34)
        #ig1_list.append(ig1i)
        #ig34_list.append(ig34i)
    #ig1 = np.array(ig1_list, dtype=nodes.dtype)
    #ig34 = np.array(ig34_list, dtype=nodes.dtype)


    face_list = []
    #print(ig1)
    #print(ig34)
    if len(g1) == 0:
        area_list = []
        return area_list
    #print('nsum =', iblank34.sum())
    #from itertools import count

    # i1 and i34; CHEXA/CPYRAM/CPENTA/CTETRA only
    #face = np.full((nloads, 4), -1, dtype='int32')
    CPYRAM_FACE_MAPPER = {
        # quad bottom
        (0, 2): np.array([0, 1, 2, 3]),
        (1, 3): np.array([0, 1, 2, 3]),
        (2, 0): np.array([0, 1, 2, 3]),  # flipped
        (3, 1): np.array([0, 1, 2, 3]),

        # tri faces???
        #(0, 1) : np.array([0, 1, 4]),
        #(1, 2) : np.array([1, 2, 4]),
        #(2, 3) : np.array([2, 3, 4]),
        #(3, 0) : np.array([3, 0, 4]),
    }
    for i, enodes, enodes_list, g1i, g34i in zip(ids_with_34,
                                                 element_nodes[i34, :],
                                                 element_nodes[i34, :].tolist(),
                                                 g1[i34], g34[i34]):
        ig1i = enodes_list.index(g1i)
        ig34i = enodes_list.index(g34i)

        assert element_type in ['CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM'], element_type
        #print(enodes_list, g1i)
        face_list = []
        if element_type == 'CHEXA':
            #For CHEXA, CPYRAM, or CPENTA quadrilateral faces, G3 is the
            #identification number of a grid point connected to a corner diagonally
            #opposite to G1. Required for quadrilateral faces of CHEXA,
            facei = CHEXA_FACE_MAPPER[(ig1i, ig34i)]
            #print(f'  face={facei}')
            #print('enodes =', enodes)
            enodes2 = enodes[facei]
            #print(f'  enodes2={enodes2}')

        elif element_type == 'CPYRAM':
            facei = CPYRAM_FACE_MAPPER[(ig1i, ig34i)]
            enodes2 = enodes[facei]
        elif element_type == 'CPENTA':
            facei = CPENTA_FACE_MAPPER[(ig1i, ig34i)]
            enodes2 = enodes[facei]
        elif element_type == 'CTETRA':
            #print('ig1=', g1i, ig1i)
            #
            #
            #
            #
            # g1 is irrelevant? - on the face
            # g4 not on the face
            # g1 loaded, g2 not loaded
            _ctetra_faces = {
                2: np.array([3, 1, 0], dtype='int32'),
                3: np.array([0, 1, 2], dtype='int32'),
                0: np.array([3, 2, 1], dtype='int32'),
                1: np.array([0, 2, 3], dtype='int32'),
            }
            for face_ in _ctetra_faces.values():
                if ig1i in face_ and ig34i not in face_:
                    found_face = face_
            #print((ig1i, ig34i), found_face)


            #face = CTETRA_FACE_MAPPER[(ig1i, ig34i)]
            face = _ctetra_faces[ig34i]
            facei = found_face
            del found_face
            #print(f'  face={facei}; type={type(face)}')
            #print(f'   enodes={enodes}; type={type(enodes)}')
            enodes2 = enodes[facei]
            #print(f'  enodes2={enodes2}')
        else:
            raise NotImplementedError(element_type)
        face_list.append(enodes2)
        assert len(face_list) > 0, face_list
        face = np.array(face_list, dtype=element_nodes.dtype)
        #print('face =', face)
        #element_nodes = element.nodes[:, face]
        #print('element_nodes =', element_nodes)
        if element_type == 'CTETRA':
            areai, centroidi, normali = _tri_area(nids, xyz_cid0, face)
        else:
            areai, centroidi, normali = _quad_area(nids, xyz_cid0, face)
        #n1 = enodes[ig1i]
        #n1 = nodes[]
        area[i] = areai
        centroid[i, :] = centroidi
        normal[i, :] = normali

    for i, enodes, enodes_list, g1i, g34i in zip(ids_no_34,
                                                 element_nodes[iblank34, :],
                                                 element_nodes[iblank34, :].tolist(),
                                                 g1[iblank34],
                                                 g34[iblank34]):
        ig1 = enodes_list.index(g1i)
        log.debug(f'{element_type}: nnodes={nnodes}; ig1={ig1}; g3/4={g34}')
        ig34 = -1
        #card.map_node_index_to_face_number(ig1, ig34)

        nids = [ig1, ig34]
        nids.sort()
        if element_type == 'CPENTA':
            #if load.g34 is None:
                #face_acn = elem.get_face_area_centroid_normal(g1)
                #nface = 3
            #else:
                #face_acn = elem.get_face_area_centroid_normal(g1, load.g34_ref.nid)
                #nface = 4
            facei = CPENTA_FACE_MAPPER[(ig1, ig34)]
            log.debug(f'  face={facei}')
        elif element_type == 'CHEXA':
            facei = CHEXA_FACE_MAPPER[(ig1, ig34)]
            log.debug(f'  face={facei}')
        elif element_type == 'CTETRA':
            log.debug('(ig1, ig34)=', (ig1, ig34))
            #
            #
            #
            #
            # g1 is irrelevant? - on the face
            # g4 not on the face
            # g1 loaded, g2 not loaded
            _ctetra_faces = (
                (3, 1, 0),
                (0, 1, 2),
                (3, 2, 1),
                (0, 2, 3),
            )
            for face in _ctetra_faces:
                if ig1 in face and ig34 not in face:
                    found_face = face
            log.debug(f'(ig1={ig1} ig34={ig34} found_face={found_face}')

            face = CTETRA_FACE_MAPPER[(ig1, ig34)]
            facei = found_face
            del found_face
            log.debug(f'  face={facei}')
            #adsf
        else:
            raise NotImplementedError(element)
        #x = 1
        #print(enodes, type(enodes))
        #print(facei, type(facei))
        enodes2 = enodes[facei]
        #print(enodes2)
        face_list.append(enodes2)

        assert len(facei) == 3, facei
        areai = 1.
        centroidi = [0., 0., 0.]
        normali = [0., 0., 1.]
        #areai, centroidi, normali = _tri_area(nids, xyz_cid0, face)
        area[i] = areai
        centroid[i, :] = centroidi
        normal[i, :] = normali

        #print('area =', area)
        del facei
    assert len(face_list) > 0, face_list

    face = np.array(face_list)
    log.debug('face = {face}')

    #if element_type == 'CHEXA':
        #face = np.array(face_list)
        #print(face)
        #n1 = face[:, 0]
        #n2 = face[:, 1]
        #n3 = face[:, 2]
        #n4 = face[:, 3]
        #crossi = np.cross(n3-n1, n4-n2)
    #elif element_type == 'CTETRA':
        #face = np.array(face_list)
        #print(face)
        #n1 = face[:, 0]
        #n2 = face[:, 1]
        #n3 = face[:, 2]
        #crossi = np.cross(n2-n1, n3-n1)
    #elif element_type == 'CTETRA':
        #face = np.array(face_list)
        #print(face)
        #n1 = face[:, 0]
        #n2 = face[:, 1]
        #n3 = face[:, 2]
        #crossi = np.cross(n2-n1, n3-n1)
    #else:
        #raise NotImplementedError(element)
    #print(crossi.shape)
    #area = np.linalg.norm(crossi, axis=0)
    #assert len(area) == len(n1)

    mean_pressure = _solid_mean_pressure(face, pressure)
    #mean_pressure_list.extend(mean_pressures)
    #n1 = face
    if np.any(np.isnan(area)):
        raise RuntimeError(f'there is nan area; area={area}')
    if np.any(np.isnan(centroid)):
        raise RuntimeError(f'there is nan centroid; centroid={centroid}')
    if np.any(np.isnan(normal)):
        raise RuntimeError(f'there is nan normal; normal={normal}')
    if np.any(np.isnan(mean_pressure)):
        raise RuntimeError(f'there is nan pressure; mean_pressure={mean_pressure}')
    return area, centroid, normal, mean_pressure



def _tri_area(nids: np.ndarray, xyz_cid0: np.ndarray, face):
    iface = np.searchsorted(nids, face)
    in1 = iface[:, 0]
    in2 = iface[:, 1]
    in3 = iface[:, 2]
    nelements = len(in1)

    n1 = xyz_cid0[in1, :]
    n2 = xyz_cid0[in2, :]
    n3 = xyz_cid0[in3, :]
    centroid = (n1 + n2 + n3) / 3.
    crossi = np.cross(n2-n1, n3-n1)
    area = np.linalg.norm(crossi, axis=1)
    assert len(area) == nelements, f'len(area)={len(area)}; nelmements={nelements}'
    normal = crossi / area[:, np.newaxis]
    assert normal.shape == (nelements, 3), normal.shape
    return area, centroid, normal

def _quad_area(nids: np.ndarray, xyz_cid0: np.ndarray, face):
    iface = np.searchsorted(nids, face)
    in1 = iface[:, 0]
    in2 = iface[:, 1]
    in3 = iface[:, 2]
    in4 = iface[:, 3]
    n1 = xyz_cid0[in1, :]
    n2 = xyz_cid0[in2, :]
    n3 = xyz_cid0[in3, :]
    n4 = xyz_cid0[in4, :]
    centroid = (n1 + n2 + n3 + n4) / 4.
    crossi = np.cross(n4-n2, n3-n1)
    area = np.linalg.norm(crossi, axis=1)
    normal = crossi / area[:, np.newaxis]
    assert len(area) == len(n1)
    return area, centroid, normal
