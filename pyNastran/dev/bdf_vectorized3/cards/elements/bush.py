from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, blank,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_double_or_blank,
    fields)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, get_print_card_8_16,
    parse_element_check, parse_property_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    get_print_card_size, array_str,
    array_default_int, array_default_float, array_float_nan)
from .rod import line_length_nan, line_centroid, line_centroid_with_spoints
from .bar import get_bar_vector, safe_normalize
from .utils import get_mass_from_property
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class CBUSH(Element):
    """
    Generalized Spring-and-Damper Connection

    Defines a generalized spring-and-damper structural element that
    may be nonlinear or frequency dependent.

    +-------+-----+------+----+----+-------+----+----+-----+
    |   1   |  2  |  3   |  4 |  5 |   6   |  7 |  8 |  9  |
    +=======+=====+======+====+====+=======+====+====+=====+
    | CBUSH | EID | PID  | GA | GB | GO/X1 | X2 | X3 | CID |
    +-------+-----+------+----+----+-------+----+----+-----+
    |       |  S  | OCID | S1 | S2 |   S3  |    |    |     |
    +-------+-----+------+----+----+-------+----+----+-----+
    """
    @Element.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            x: Optional[list[float]], g0: Optional[int], cid=None,
            s: float=0.5, ocid: int=-1, si: Optional[list[float]]=None, comment='') -> int:
        """
        Creates a CBUSH card

        Parameters
        ----------
        eid : int
            Element id
        pid : int
            Property id (PBUSH)
        nids : list[int, int]
            node ids; connected grid points at ends A and B
            The nodes may be coincident, but then cid is required.
        x : list[float, float, float]; None
            list : the directional vector used to define the stiffnesses
                   or damping from the PBUSH card
            None : use g0
        g0 : int/None
            int : the directional vector used to define the stiffnesses
                  or damping from the PBUSH card
            None : use x
        cid : int; default=None
            Element coordinate system identification. A 0 means the basic
            coordinate system. If CID is blank, then the element coordinate
            system is determined from GO or Xi.
        s: float; default=0.5
            Location of spring damper (0 <= s <= 1.0)
        ocid : int; default=-1
            Coordinate system identification of spring-damper offset.
            (Integer > -1; Default = -1, which means the offset
            point lies on the line between GA and GB)
        si : list[float, float, float]; default=None
            Components of spring-damper offset in the OCID coordinate system
            if OCID > 0.
            None : [None, None, None]
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, cid, g0, x, s, ocid, si, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb', default=0)

        #: Element coordinate system identification. A 0 means the basic
        #: coordinate system. If CID is blank, then the element coordinate
        #: system is determined from GO or Xi.
        #: (default=blank=element-based)
        cid = integer_or_blank(card, 8, 'cid', default=-1)

        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0')
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x1 = x1_g0
            x2 = double_or_blank(card, 6, 'x2', default=0.0)
            x3 = double_or_blank(card, 7, 'x3', default=0.0)
            x = [x1, x2, x3]
            if not isinstance(cid, integer_types):
                assert max(x) != min(x), 'x=%s' % x
        else:
            g0 = -1
            x = [None, None, None]

        #: Location of spring damper (0 <= s <= 1.0)
        s = double_or_blank(card, 9, 's', default=0.5)

        #: Coordinate system identification of spring-damper offset. See
        #: Remark 9. (Integer > -1; Default = -1, which means the offset
        #: point lies on the line between GA and GB
        ocid = integer_or_blank(card, 10, 'ocid', default=-1)

        #: Components of spring-damper offset in the OCID coordinate system
        #: if OCID > 0.
        si = [double_or_blank(card, 11, 's1'),
              double_or_blank(card, 12, 's2'),
              double_or_blank(card, 13, 's3'), ]
        assert len(card) <= 14, f'len(CBUSH card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, [ga, gb], cid, g0, x, s, ocid, si, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        g0 = np.full(ncards, -1, dtype=idtype)
        x = np.full((ncards, 3), np.nan, dtype='float64')

        coord_id = np.zeros(ncards, dtype='int32')
        s = np.zeros(ncards, dtype='float64')
        ocid = np.zeros(ncards, dtype='int32')
        ocid_offset = np.zeros((ncards, 3), dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nodesi, cidi, g0i, xi, si, ocidi, ocid_offseti, comment) = card
            if cidi is None:
                cidi = -1

            s[icard] = si
            coord_id[icard] = cidi
            ocid[icard] = ocidi
            ocid_offset[icard, :] = ocid_offseti

            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nodesi
            if g0i is None:
                x[icard, :] = xi
            else:
                g0[icard] = g0i
        self._save(element_id, property_id, nodes, x, g0,
                   s, coord_id, ocid, ocid_offset)
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              x: np.ndarray, g0: np.ndarray,
              s: np.ndarray,
              coord_id: np.ndarray,
              ocid: np.ndarray, ocid_offset: np.ndarray) -> None:
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.x = x
        self.g0 = g0
        self.s = s
        self.coord_id = coord_id
        self.ocid = ocid
        self.ocid_offset = ocid_offset

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)
        coords = np.unique(np.hstack([self.coord_id, self.ocid]))
        coords = coords[coords >= 0]
        used_dict['coord_id'].append(coords)

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        """
        s: float; default=0.5
            Location of spring damper (0 <= s <= 1.0)
        si : list[float, float, float]; default=None
            Components of spring-damper offset in the OCID coordinate system
            if OCID > 0.
            None : [None, None, None]
        """
        icoord = (self.coord_id != -1)
        ncoord = icoord.sum()

        is_x = self.is_x

        #xmax = self.x.max(axis=1)
        #xnan = np.isnan(xmax)
        nx = is_x.sum()
        if nx:
            ga = self.nodes[is_x, 0]
            grid_ga = self.model.grid.slice_card_by_id(ga, assume_sorted=True, sort_ids=False)
            cp = grid_ga.cp
            coords = self.model.coord.slice_card_by_id(cp, assume_sorted=True, sort_ids=False)
            xi = self.x[is_x, :].copy()
            ixyz = (coords.coord_type == 'R')
            irtz = (coords.coord_type == 'C')
            irtp = (coords.coord_type == 'S')
            xi[ixyz, :] *= xyz_scale
            xi[irtz, 0] *= xyz_scale
            xi[irtz, 2] *= xyz_scale
            xi[irtp, 0] *= xyz_scale
            self.x[is_x] = xi

        #if ncoord:
            #asfd
            #ib = self.offt
        self.s *= xyz_scale
        self.si *= xyz_scale

    def __apply_slice__(self, elem: CBUSH, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.g0 = self.g0[i]
        elem.x = self.x[i, :]
        elem.s = self.s[i]
        elem.coord_id = self.coord_id[i]
        elem.ocid = self.ocid[i]
        elem.ocid_offset = self.ocid_offset[i, :]
        elem.n = len(i)

    @property
    def is_x(self) -> np.ndarray:
        return (self.g0 == -1)

    @property
    def is_g0(self) -> np.ndarray:
        return ~self.is_x

    @property
    def si(self) -> np.ndarray:
        return self.ocid_offset
    @si.setter
    def si(self, ocid_offset: np.ndarray) -> None:
        self.ocid_offset = ocid_offset

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodess = array_str(self.nodes, size=size)
        ocids = array_default_int(self.ocid, default=-1, size=size)
        cids = array_default_int(self.coord_id, default=-1, size=size)
        ss = array_default_float(self.s, default=0.5, size=size, is_double=False)
        for eid, pid, nodes, g0, x, cid, s, ocid, si in zip(element_ids, property_ids, nodess,
                                                            self.g0, self.x, cids, ss, ocids, self.si):
            n1, n2 = nodes

            if g0 == -1:
                x1 = ''
                x2 = ''
                x3 = ''
            elif g0 > 0:
                x1 = g0
                x2 = ''
                x3 = ''
            else:
                x1, x2, x3 = x # self.get_x_g0_defaults()
                assert not np.any(np.isnan(x)), x

            si1, si2, si3 = '', '', ''
            if not np.all(np.isnan(si)):
                si1, si2, si3 = si
            list_fields = ['CBUSH', eid, pid, n1, n2,
                           x1, x2, x3,
                           cid, s, ocid, si1, si2, si3]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no bush properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pbush]
                if prop.n > 0]

    def mass(self) -> np.ndarray:
        mass = get_mass_from_property(self.property_id, self.allowed_properties)
        return mass

    def get_xyz(self) -> tuple[np.ndarray, np.ndarray]:
        #neids = len(self.element_id)
        grid = self.model.grid
        xyz = grid.xyz_cid0()
        nid = grid.node_id
        inode = np.searchsorted(nid, self.nodes)
        assert np.array_equal(nid[inode], self.nodes)
        in1 = inode[:, 0]
        in2 = inode[:, 1]
        xyz1 = xyz[in1, :]
        xyz2 = xyz[in2, :]
        return xyz1, xyz2

    def get_bar_vector(self, xyz1: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        idefault = (self.coord_id == -1)
        icoord = ~idefault

        i = np.full(self.x.shape, np.nan, self.x.dtype)
        j = np.full(self.x.shape, np.nan, self.x.dtype)
        k = np.full(self.x.shape, np.nan, self.x.dtype)
        xyz1, xyz2 = self.get_xyz()
        i_vector = xyz2 - xyz1
        nelement = self.n
        maxs = np.abs(i_vector).max(axis=1)
        assert len(maxs) == nelement, (len(maxs), nelement)

        if np.any(idefault):
            sub_i = i_vector[idefault, :]
            ihat = safe_normalize(sub_i)

            index = np.where(idefault)[0]
            sub_cbush = self.slice_card_by_index(index)
            v, cd = get_bar_vector(sub_cbush, xyz1[idefault, :])

            ki = np.cross(ihat, v, axis=1)
            khat = safe_normalize(ki)
            jhat = np.cross(khat, ihat, axis=1)
            i[idefault, :] = ihat
            j[idefault, :] = jhat
            k[idefault, :] = khat

        if np.any(icoord):
            coord_ids = self.coord_id[icoord]
            assert len(coord_ids) == icoord.sum()
            coords = self.model.coord.slice_card_by_id(coord_ids)
            ihat = coords.i
            jhat = coords.j
            khat = coords.k
            i[icoord, :] = ihat
            j[icoord, :] = jhat
            k[icoord, :] = khat
        #if 0:
            #v, cd = get_bar_vector(self, xyz1)
        return i, j, k

    def get_axes(self, xyz1: np.ndarray, xyz2: np.ndarray,
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                            np.ndarray, np.ndarray, np.ndarray]:
        log = self.model.log
        coords = self.model.coord
        #xyz1, xyz2 = self.get_xyz()

        neids = xyz1.shape[0]
        #i = xyz2 - xyz1
        #ihat_norm = np.linalg.norm(i, axis=1)
        #assert len(ihat_norm) == neids
        #if min(ihat_norm) == 0.:
            #msg = 'xyz1=%s xyz2=%s\n%s' % (xyz1, xyz2, self)
            #log.error(msg)
            #raise ValueError(msg)
        #i_offset = i / ihat_norm[:, np.newaxis]

        #log.info(f'x =\n{self.x}')
        #log.info(f'g0   = {self.g0}')
        ihat, yhat, zhat = self.get_bar_vector(xyz1)

        v = yhat
        wa = np.zeros(ihat.shape, ihat.dtype)
        wb = wa
        #ihat = xform[0, :]
        #yhat = xform[1, :]
        #zhat = xform[2, :]
        #wa, wb, _ihat, jhat, khat = out

        # we finally have the nodal coordaintes!!!! :)
        return v, ihat, yhat, zhat, wa, wb

    def length(self) -> np.ndarray:
        length = line_length_nan(self.model, self.nodes, default_node=-1)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid_with_spoints(self.model, self.nodes)
        return centroid


class PBUSH(Property):
    """
    Generalized Spring-and-Damper Property
    Defines the nominal property values for a generalized spring-and-damper
    structural element.

    +-------+-----+-------+------+-------+-----+-----+-----+----+
    |   1   |  2  |   3   |   4  |   5   |  6  |  7  |  8  | 9  |
    +=======+=====+=======+======+=======+=====+=====+=====+====+
    | PBUSH | PID |   K   |  K1  |   K2  |  K3 |  K4 |  K5 | K6 |
    +-------+-----+-------+------+-------+-----+-----+-----+----+
    |       |  B  |  B1   |  B2  |   B3  |  B4 |  B5 |  B6 |    |
    +-------+-----+-------+------+-------+-----+-----+-----+----+
    |       | GE  |  GE1  |  GE2 |  GE3  | GE4 | GE5 | GE6 |    |
    +-------+-----+-------+------+-------+-----+-----+-----+----+
    |       | RCV |  SA   |  ST  |   EA  |  ET |     |     |    |
    +-------+-----+-------+------+-------+-----+-----+-----+----+
    |       |  M  |  MASS |      |       |     |     |     |    |
    +-------+-----+-------+------+-------+-----+-----+-----+----+
    |       |  T  | ALPHA | TREF | COINL |     |     |     |    |
    +-------+-----+-------+------+-------+-----+-----+-----+----+

    RCV was added <= MSC 2016
    MASS was added <= MSC 2016
    T/ALPHA/TREF/COINL was added in MSC 2021
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.k_fields = np.zeros((0, 6), dtype='float64')
        self.b_fields = np.zeros((0, 6), dtype='float64')
        self.ge_fields = np.zeros((0, 6), dtype='float64')
        self.rcv_fields = np.zeros((0, 4), dtype='float64')
        self._mass = np.array([], dtype='float64')
        self.alpha = np.array([], dtype='float64')
        self.tref = np.array([], dtype='float64')
        self.coincident_length = np.array([], dtype='float64')

    def add(self, pid: int, k: list[float], b: list[float], ge: list[float],
            rcv: Optional[list[float]]=None, mass: Optional[float]=None,
            alpha: float=0., tref: float=0., coincident_length=None,
            comment: str='') -> int:
        """
        Creates a PBUSH card, which defines a property for a PBUSH

        Parameters
        ----------
        pid : int
            property id
        k : list[float]
            Nominal stiffness values in directions 1 through 6.
            len(k) = 6
        b : list[float]
            Nominal damping coefficients in direction 1 through 6 in units of
            force per unit velocity
            len(b) = 6
        ge : list[float]
            Nominal structural damping constant in directions 1 through 6.
            len(ge) = 6
        rcv : list[float]; default=None -> (None, None, None, None)
            [sa, st, ea, et] = rcv
            length(rcv) = 4
        mass : float; default=None
            lumped mass of the CBUSH
            This is an MSC only parameter.
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, k, b, ge, rcv,
                           mass, alpha, tref, coincident_length, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        k_fields = []
        b_fields = []
        ge_fields = []
        t_fields = [np.nan, np.nan, np.nan]
        rcv_fields = [1., 1., 1., 1.]
        mass = None

        pid = integer(card, 1, 'pid')

        nfields = card.nfields
        istart = 2
        while istart < nfields:
            pname = string(card, istart, 'pname')

            if pname == 'K':
                # Flag indicating that the next 1 to 6 fields are stiffness values in
                # the element coordinate system.
                #self.k = string(card, istart, 'k')

                #: Nominal stiffness values in directions 1 through 6.
                #: See Remarks 2 and 3. (Real; Default = 0.0)
                k_fields = self._read_var(card, 'Ki', istart + 1, istart + 7)

            elif pname == 'B':
                # Flag indicating that the next 1 to 6 fields are force-per-velocity
                # damping.
                #self.b = string(card, istart, 'b')

                #: Force per unit velocity (Real)
                #: Nominal damping coefficients in direction 1 through 6 in units of
                #: force per unit velocity. See Remarks 2, 3, and 9. (Real; Default=0.)
                b_fields = self._read_var(card, 'Bi', istart + 1, istart + 7)

            elif pname == 'GE':
                # Flag indicating that the next fields, 1 through 6 are structural
                # damping constants. See Remark 7. (Character)
                #self.ge = string(card, istart, 'ge')

                #: Nominal structural damping constant in directions 1 through 6. See
                #: Remarks 2. and 3. (Real; Default = 0.0)
                ge_fields = self._read_var(card, 'GEi', istart + 1, istart + 7)
            elif pname == 'RCV':
                rcv_fields = self._read_rcv(card, istart)
            elif pname == 'M':
                # Lumped mass of the cbush; default=0.0
                mass = double_or_blank(card, istart + 1, 'mass', default=0.)
            elif pname == 'T':
                t_fields = self._read_var(card, 'Ti', istart + 1, istart + 4)
                assert len(t_fields) == 3, t_fields
            else:
                raise RuntimeError(f'unsupported PBUSH type; pname={pname!r}')
                #break #  old version...
            istart += 8

        alpha, tref, coincident_length = t_fields
        self.cards.append((pid, k_fields, b_fields, ge_fields, rcv_fields,
                           mass, alpha, tref, coincident_length, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        property_id = np.zeros(ncards, dtype=idtype)

        k_fields = np.zeros((ncards, 6), dtype='float64')
        b_fields = np.zeros((ncards, 6), dtype='float64')
        ge_fields = np.zeros((ncards, 6), dtype='float64')
        rcv_fields = np.zeros((ncards, 4), dtype='float64')
        _mass = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')
        tref = np.zeros(ncards, dtype='float64')
        coincident_length = np.zeros(ncards, dtype='float64')
        for icard, card in enumerate(self.cards):
            (pid, k_fieldsi, b_fieldsi, ge_fieldsi, rcv_fieldsi,
             massi, alphai, trefi, coincident_lengthi, comment) = card

            property_id[icard] = pid
            if k_fieldsi:
                k_fields[icard, :] = k_fieldsi
            if b_fieldsi:
                b_fields[icard, :] = b_fieldsi
            if ge_fieldsi:
                ge_fields[icard, :] = ge_fieldsi
            #t_fields[icard, :] = t_fields
            rcv_fields[icard, :] = rcv_fieldsi
            _mass[icard] = massi
            alpha[icard] = alphai
            tref[icard] = trefi
            coincident_length[icard] = coincident_lengthi
        self._save(property_id, k_fields, b_fields, ge_fields, rcv_fields,
                   _mass, alpha, tref, coincident_length)
        self.cards = []

    def _save(self, property_id, k_fields, b_fields, ge_fields, rcv_fields,
              _mass, alpha, tref, coincident_length) -> None:
        self.property_id = property_id
        self.k_fields = k_fields
        self.b_fields = b_fields
        self.ge_fields = ge_fields
        self.rcv_fields = rcv_fields
        self._mass = _mass
        self.alpha = alpha
        self.tref = tref
        self.coincident_length = coincident_length

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['pbusht_id'].append(self.property_id)

    def convert(self,
                xyz_scale: float=1.0,
                mass_scale: float=1.0,
                temperature_scale: float=1.0,
                linear_stiffness_scale: float= 1.0,
                rotational_stiffness_scale: float= 1.0,
                linear_damping_scale: float= 1.0,
                rotational_damping_scale: float= 1.0,
                alpha_scale: float=1.0, **kwargs) -> None:
        self._mass *= mass_scale
        self.k_fields[:, [0, 1, 2]] *= linear_stiffness_scale
        self.k_fields[:, [3, 4, 5]] *= rotational_stiffness_scale

        self.b_fields[:, [0, 1, 2]] *= linear_damping_scale
        self.b_fields[:, [3, 4, 5]] *= rotational_damping_scale
        self.coincident_length *= temperature_scale
        self.tref *= xyz_scale
        self.alpha *= alpha_scale

    @classmethod
    def _read_var(cls, card, var_prefix, istart, iend):
        ki = fields(double_or_blank, card, var_prefix, istart, iend)
        return ki

    @classmethod
    def _read_rcv(cls, card, istart):
        """
        Flag indicating that the next 1 to 4 fields are stress or strain
        coefficients. (Character)
        """
        #self.rcv = string(card, istart, 'rcv')
        sa = double_or_blank(card, istart + 1, 'sa', default=1.)
        st = double_or_blank(card, istart + 2, 'st', default=1.)
        ea = double_or_blank(card, istart + 3, 'ea', default=1.)
        et = double_or_blank(card, istart + 4, 'et', default=1.)
        return [sa, st, ea, et]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass
        #nid = self.model.grid.node_id
        #pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          #msg=f'no spring properties for {self.type}')
        #pids.sort()
        #geom_check(self,
                   #missing,
                   #node=(nid, self.nodes),
                   #property_id=(pids, self.property_id))

    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #RCV was added <= MSC 2016
        #MASS was added <= MSC 2016
        #T/ALPHA/TREF/COINL was added in MSC 2021
        no_rcv = self.rcv_fields.max() == 1. and self.rcv_fields.min() == 1.
        no_t = np.any(np.isfinite(self.alpha)) | np.any(np.isfinite(self.tref)) | np.any(np.isfinite(self.coincident_length))
        no_mass = np.any(np.isfinite(self._mass))
        if no_rcv and no_t and no_mass:
            # MSC 2005r2
            for pid, k_fields, b_fields, ge_fields in zip(self.property_id, self.k_fields, self.b_fields,
                                                          self.ge_fields):
                list_fields = ['PBUSH', pid]
                _set_fields_pbush(list_fields, 'K', k_fields)
                _set_fields_pbush(list_fields, 'B', b_fields)
                _set_fields_pbush(list_fields, 'GE', ge_fields)
                bdf_file.write(print_card(list_fields))
        #elif no_t:
        else:
            for pid, k_fields, b_fields, ge_fields, rcv_fields, alpha, tref, coincident_length, mass in zip(
                    self.property_id, self.k_fields, self.b_fields,
                    self.ge_fields, self.rcv_fields,
                    self.alpha, self.tref, self.coincident_length, self._mass):
                #self.k_fields[icard, :] = k_fields
                #if b_fields:
                    #self.b_fields[icard, :] = b_fields
                #if ge_fields:
                    #self.ge_fields[icard, :] = ge_fields
                #self.rcv_fields[icard, :] = rcv_fields
                #self.mass[icard] = mass
                #self.alpha[icard] = alpha
                #self.tref[icard] = tref
                #self.coinl[icard] = coinl

                mass = None if mass == 0. else mass
                rcv_fields = [None if rcv == 1. else rcv
                              for rcv in rcv_fields]
                list_fields = ['PBUSH', pid]
                _set_fields_pbush(list_fields, 'K', k_fields)
                _set_fields_pbush(list_fields, 'B', b_fields)
                _set_fields_pbush(list_fields, 'GE', ge_fields)
                _set_fields_pbush(list_fields, 'RCV', rcv_fields)
                _set_fields_pbush(list_fields, 'M', [mass])
                _set_fields_pbush(list_fields, 'T', [alpha, tref, coincident_length])

                bdf_file.write(print_card(list_fields))
        return

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass(self) -> np.ndarray:
        mass = self._mass
        return mass


class PBUSHT(Property):
    #def __init__(self, pid, k_tables, b_tables,
                 #ge_tables, kn_tables, comment=''):
        #BushingProperty.__init__(self)
        #if comment:
            #self.comment = comment
        #self.pid = pid

        #_append_nones(k_tables, 6)
        #_append_nones(b_tables, 6)
        #_append_nones(ge_tables, 6)
        #_append_nones(kn_tables, 6)

        #self.k_tables = k_tables
        #self.b_tables = b_tables
        #self.ge_tables = ge_tables
        #self.kn_tables = kn_tables

    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.k_tables = np.zeros((0, 6), dtype='int32')
        self.b_tables = np.zeros((0, 6), dtype='int32')
        self.ge_tables = np.zeros((0, 6), dtype='int32')
        self.kn_tables = np.zeros((0, 6), dtype='int32')

    def add(self, pid: int, k_tables: list[int], b_tables: list[int],
            ge_tables: list[int], kn_tables: list[int], comment: str='') -> int:
        self.cards.append((pid, k_tables, b_tables, ge_tables, kn_tables, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PBUSHT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        k_tables = [0] * 6
        b_tables = [0] * 6
        ge_tables = [0] * 6
        kn_tables = [0] * 6

        pid = integer(card, 1, 'pid')
        nfields = len(card) - 1
        nrows = nfields // 8
        if nfields % 8 != 0:
            nrows += 1

        for irow in range(nrows):
            ifield = 1 + irow * 8
            param = string(card, ifield + 1, 'param_type')
            table = []
            for j in range(6):
                table_value = integer_or_blank(card, ifield + j + 2, param + '%i' % (j+1),
                                               default=0)
                table.append(table_value)
            if param == 'K':
                k_tables = table
            elif param == 'B':
                b_tables = table
            elif param == 'GE':
                ge_tables = table
            elif param == 'KN':
                kn_tables = table
            else:
                raise ValueError(param)
        self.cards.append((pid, k_tables, b_tables, ge_tables, kn_tables, comment))
        self.n += 1
        return self.n - 1
        #return PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables,
                      #comment=comment)

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        k_tables = np.zeros((ncards, 6), dtype='int32')
        b_tables = np.zeros((ncards, 6), dtype='int32')
        ge_tables = np.zeros((ncards, 6), dtype='int32')
        kn_tables = np.zeros((ncards, 6), dtype='int32')
        for icard, card in enumerate(self.cards):
            pid, k_tablesi, b_tablesi, ge_tablesi, kn_tablesi, comment = card
            property_id[icard] = pid
            k_tables[icard, :] = k_tablesi
            b_tables[icard, :] = b_tablesi
            ge_tables[icard, :] = ge_tablesi
            kn_tables[icard, :] = kn_tablesi
        self._save(property_id, k_tables, b_tables, ge_tables, kn_tables)

    def _save(self, property_id, k_tables, b_tables, ge_tables, kn_tables):
        ncards_existing = len(self.property_id)
        if ncards_existing != 0:
            property_id = np.hstack([self.property_id, property_id])
            k_tables = np.vstack([self.k_tables, k_tables])
            b_tables = np.vstack([self.b_tables, b_tables])
            ge_tables = np.vstack([self.ge_tables, ge_tables])
            kn_tables = np.vstack([self.kn_tables, kn_tables])

        self.property_id = property_id
        self.k_tables = k_tables
        self.b_tables = b_tables
        self.ge_tables = ge_tables
        self.kn_tables = kn_tables

    def __apply_slice__(self, prop: PBUSHT, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop.k_tables = self.k_tables[i, :]
        prop.b_tables = self.b_tables[i, :]
        prop.ge_tables = self.ge_tables[i, :]
        prop.kn_tables = self.kn_tables[i, :]
        prop.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['tabled_id'].append(self.k_tables.ravel())
        used_dict['tabled_id'].append(self.b_tables.ravel())
        used_dict['tabled_id'].append(self.ge_tables.ravel())
        used_dict['tabled_id'].append(self.kn_tables.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        k_tables = array_default_int(self.k_tables, default=0, size=size).tolist()
        b_tables = array_default_int(self.b_tables, default=0, size=size).tolist()
        ge_tables = array_default_int(self.ge_tables, default=0, size=size).tolist()
        kn_tables = array_default_int(self.kn_tables, default=0, size=size).tolist()
        is_k = self.k_tables.max(axis=1) > 0
        is_b = self.b_tables.max(axis=1) > 0
        is_ge = self.ge_tables.max(axis=1) > 0
        is_kn = self.kn_tables.max(axis=1) > 0
        assert len(is_k) == len(property_id)
        for pid, ik, ib, ige, ikn, k_tablesi, b_tablesi, ge_tablesi, kn_tablesi in \
                zip(property_id,
                    is_k, is_b, is_ge, is_kn,
                    k_tables, b_tables, ge_tables, kn_tables):
            list_fields = ['PBUSHT', pid]
            if ik:
                list_fields += ['K'] + k_tablesi + [None]
            if ib:
                list_fields += ['B'] + b_tablesi + [None]
            if ige:
                list_fields += ['GE'] + ge_tablesi + [None]
            if ikn:
                list_fields += ['KN'] + kn_tablesi + [None]
            bdf_file.write(print_card(list_fields))
        return

def _set_fields_pbush(list_fields: list[Any], var: str, fields: list[Any]):
    """helper method for PBUSH"""
    fields2 = []
    write_fields = False
    for field in fields:
        if field is None or np.isnan(field):
            fields2.append(None)
        else:
            fields2.append(field)
            write_fields = True
    if not write_fields:
        return
    list_fields += [var] + fields2

    nspaces = 8 - (len(list_fields) - 1) % 8
    if nspaces == 8:
        list_fields += [None]
    elif nspaces < 8:
        list_fields += [None] * (nspaces + 1)


class CBUSH1D(Element):
    @Element.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int], cid: int=-1,
            comment: str='') -> int:
        if cid is None:
            cid = -1
        self.cards.append((eid, pid, nids, cid, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb', default=-1)
        cid = integer_or_blank(card, 5, 'cid', default=-1)
        assert len(card) <= 6, f'len(CBUSH1D card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, [ga, gb], cid, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        coord_id = np.full(ncards, -1, dtype='int32')
        for icard, card in enumerate(self.cards):
            (eid, pid, nodesi, cid, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nodesi
            coord_id[icard] = cid
        self._save(element_id, property_id, nodes, coord_id)
        self.cards = []

    def _save(self, element_id, property_id, nodes, coord_id) -> None:
        if len(self.element_id) != 0:
            asdf
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.coord_id = coord_id

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)
        coords = np.unique(self.coord_id)
        coords = coords[coords >= 0]
        used_dict['coord_id'].append(coords)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no pbush1d properties for {self.type}')
        pids.sort()
        cids = self.model.coord.coord_id
        coord_id = self.coord_id[self.coord_id != -1]
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id),
                   coord=(cids, coord_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodess = array_str(self.nodes, size=size)
        cids = array_default_int(self.coord_id, default=-1, size=size)
        for eid, pid, nodes, cid in zip(element_ids, property_ids, nodess, cids):
            n1, n2 = nodes
            list_fields = ['CBUSH1D', eid, pid, n1, n2, cid]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pbush1d]
                if prop.n > 0]

    #def mass(self) -> np.ndarray:
        #mass = get_mass_from_property(self.property_id, self.allowed_properties)
        #return mass

    def length(self) -> np.ndarray:
        length = line_length_nan(self.model, self.nodes, default_node=-1)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid


class PBUSH1D(Property):
    """
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |   1     |    2   |   3   |    4   |    5   |   6   |   7   |   8   |
    +=========+========+=======+========+========+=======+=======+=======+
    | PBUSH1D |   PID  |   K   |    C   |    M   |       |   SA  |   SE  |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | SHOCKA | TYPE  |   CVT  |   CVC  | EXPVT | EXPVC |  IDTS |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | IDETS  | IDECS | IDETSD | IDECSD |       |       |       |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | SPRING |  TYPE |   IDT  |   IDC  | IDTDU | IDCDU |       |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | DAMPER |  TYPE |   IDT  |   IDC  | IDTDV | IDCDV |       |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    |         | GENER  |  IDT  |   IDC  |  IDTDU | IDCDU | IDTDV | IDCDV |
    +---------+--------+-------+--------+--------+-------+-------+-------+
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.k = np.array([], dtype='float64')
        self.c = np.array([], dtype='float64')
        self.sa = np.array([], dtype='float64')
        self.se = np.array([], dtype='float64')
        self.mass = np.array([], dtype='float64')

        self.spring_type = np.array([], dtype='|U8')
        self.spring_table = np.zeros((0, 4), dtype='int32')
        self.spring_equation = np.zeros((0, 4), dtype='int32')

        self.damper_type = np.array([], dtype='|U8')
        self.damper_table = np.zeros((0, 4), dtype='int32')
        self.damper_equation = np.zeros((0, 4), dtype='int32')

        self.shock_type = np.array([], dtype='|U8')
        self.shock_table = np.zeros((0, 5), dtype='int32')
        self.shock_equation = np.zeros((0, 4), dtype='int32')

        # TYPE = EQUAT
        self.gener_equation = np.zeros((0, 6), dtype='int32')

    def add(self, pid: int,
            k: float=0., c: float=0., m: float=0.,
            sa: float=0., se: float=0., optional_vars=None,
            comment: str='') -> int:
        """Creates a PBUSH1D card"""

        if optional_vars is None:
            optional_vars = {}
        self.cards.append((pid, k, c, m, sa, se, optional_vars, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        k = double_or_blank(card, 2, 'k', default=0.0)
        c = double_or_blank(card, 3, 'c', default=0.0)
        m = double_or_blank(card, 4, 'm', default=0.0)

        sa = double_or_blank(card, 6, 'sa', default=0.0)
        se = double_or_blank(card, 7, 'se', default=0.0)

        nfields = card.nfields
        optional_vars = {}
        istart = 9
        while istart < nfields:
            pname = string(card, istart, 'pname')
            if pname == 'SHOCKA':
                istart, out = self._read_shock(card, istart)
                optional_vars['SHOCKA'] = out

            elif pname == 'SPRING':
                out = self._read_spring(card, istart)
                optional_vars['SPRING'] = out

            elif pname == 'DAMPER':
                out = self._read_damper(card, istart)
                optional_vars['DAMPER'] = out

            elif pname == 'GENER':
                out = self._read_gener(card, istart)
                optional_vars['GENER'] = out
            else:
                break
            istart += 8

        #return PBUSH1D(pid, k=k, c=c, m=m, sa=sa, se=se,
                       #optional_vars=optional_vars, comment=comment)
        self.cards.append((pid, k, c, m, sa, se, optional_vars, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        property_id = np.zeros(ncards, dtype=idtype)

        k = np.zeros(ncards, dtype='float64')
        c = np.zeros(ncards, dtype='float64')
        sa = np.zeros(ncards, dtype='float64')
        se = np.zeros(ncards, dtype='float64')
        mass = np.zeros(ncards, dtype='float64')

        spring_type = np.zeros(ncards, dtype='|U8')
        spring_table = np.zeros((ncards, 4), dtype='int32')
        spring_equation = np.zeros((ncards, 4), dtype='int32')

        damper_type = np.zeros(ncards, dtype='|U8')
        damper_table = np.zeros((ncards, 4), dtype='int32')
        damper_equation = np.zeros((ncards, 4), dtype='int32')

        shock_type = np.zeros(ncards, dtype='|U8')
        shock_table = np.zeros((ncards, 5), dtype='int32')
        shock_equation = np.zeros((ncards, 4), dtype='int32')

        # TYPE = EQUAT
        gener_equation = np.zeros((ncards, 6), dtype='int32')
        #self.gener_idt = gener_idt
        #self.gener_idc = gener_idc
        #self.gener_idtdu = gener_idtdu
        #self.gener_idcdu = gener_idcdu
        #self.gener_idtdv = gener_idtdv
        #self.gener_idcdv = gener_idcdv

        for icard, card in enumerate(self.cards):
            (pid, ki, ci, massi, sai, sei, optional_vars, comment) = card
            property_id[icard] = pid
            k[icard] = ki
            c[icard] = ci
            mass[icard] = massi
            sa[icard] = sai
            se[icard] = sei

            self.model.log.info(f'PBUSH pid={pid}')
            for key, values in optional_vars.items():
                self.model.log.info(f'  {key} values={values}')
                if key == 'SPRING':
                    spring_typei, *spring_data = values
                    spring_type[icard] = spring_typei
                    if spring_typei == 'TABLE':
                        spring_table[icard, :] = spring_data
                    elif spring_typei == 'EQUAT':
                        spring_equation[icard, :] = spring_data
                    else:
                        raise RuntimeError(values)

                elif key == 'DAMPER':
                    damper_typei, *damper_data = values
                    damper_type[icard] = damper_typei
                    if damper_typei == 'TABLE':
                        damper_table[icard, :] = damper_data
                    elif damper_typei == 'EQUAT':
                        damper_equation[icard, :] = damper_data
                    else:
                        raise RuntimeError(values)
                elif key == 'SHOCKA':
                    #shock_cvt = double(card, istart + 2, 'shockCVT')
                    #shock_cvc = double_or_blank(card, istart + 3, 'shockCVC', default=shock_cvt)
                    #shock_exp_vt = double_or_blank(card, istart + 4, 'shockExpVT', default=1.0)
                    #shock_exp_vc = double_or_blank(card, istart + 5,

                    shock_typei, shock_cvt, scock_cvc, shock_exp_vt, shock_exp_vc, *shock_data = values
                    shock_type[icard] = shock_typei
                    if shock_typei == 'TABLE':
                        self.model.log.debug(f'shock_data={shock_data}')
                        shock_table[icard, :] = shock_data
                    elif shock_typei == 'EQUAT':
                        shock_equation[icard, :] = shock_data
                    else:
                        raise RuntimeError(values)
                elif key == 'GENER':
                    #  F(u, v) = Ft(u, v)
                    #IDT IDC IDTDU IDCDU IDTDV IDCDV
                    gener_idt, gener_idc, gener_idtdu, gener_idcdu, gener_idtdv, gener_idcdv = values
                    print('values =', values)

                    #gener_idt   = integer(card, istart + 2, 'generIDT')
                    #gener_idc   = integer_or_blank(card, istart + 3, 'generIDC', gener_idt)
                    #gener_idtdu = integer(card, istart + 4, 'generIDTDU')
                    #gener_idcdu = integer_or_blank(card, istart + 5, 'generIDCDU', gener_idtdu)
                    #gener_idtdv = integer(card, istart + 6, 'generIDTDV')
                    #gener_idcdv = integer_or_blank(card, istart + 7, 'generIDCDV', gener_idtdv)
                    #shock_cvt = double(card, istart + 2, 'shockCVT')
                    #shock_cvc = double_or_blank(card, istart + 3, 'shockCVC', default=shock_cvt)
                    #shock_exp_vt = double_or_blank(card, istart + 4, 'shockExpVT', default=1.0)
                    #shock_exp_vc = double_or_blank(card, istart + 5,

                    self.model.log.debug(f'gener_data={values}')
                    gener_equation[icard, :] = values
                    #self.gener_idt[icard, :] = gener_idt
                    #self.gener_idc[icard, :] = gener_idc
                    #self.gener_idtdu[icard, :] = gener_idtdu
                    #self.gener_idcdu[icard, :] = gener_idcdu
                    #self.gener_idtdv[icard, :] = gener_idtdv
                    #self.gener_idcdv[icard, :] = gener_idcdv
                else:
                    raise NotImplementedError(key)
            #assert len(optional_vars) == 0, optional_vars
        self._save(property_id, k, c, sa, se, mass,
                   spring_type, spring_table, spring_equation,
                   damper_type, damper_table, damper_equation,
                   shock_type, shock_table, shock_equation,
                   gener_equation,
                   #gener_idt, gener_idc,
                   #gener_idtdu, gener_idcdu, gener_idtdv, gener_idcdv,
                   )
        self.sort()
        assert len(self.property_id) > 0, self.property_id
        self.cards = []

    def _save(self, property_id, k, c, sa, se, mass,
              spring_type, spring_table, spring_equation,
              damper_type, damper_table, damper_equation,
              shock_type, shock_table, shock_equation,
              gener_equation,
              ):
        if len(self.property_id) != 0:
            asdf
        self.property_id = property_id
        self.k = k
        self.c = c
        self.sa = sa
        self.se = se
        self.mass = mass

        self.spring_type = spring_type
        self.spring_table = spring_table
        self.spring_equation = spring_equation

        self.damper_type = damper_type
        self.damper_table = damper_table
        self.damper_equation = damper_equation

        self.shock_type = shock_type
        self.shock_table = shock_table
        self.shock_equation = shock_equation

        self.gener_equation = gener_equation
        #self.gener_idt = gener_idt
        #self.gener_idc = gener_idc
        #self.gener_idtdu = gener_idtdu
        #self.gener_idcdu = gener_idcdu
        #self.gener_idtdv = gener_idtdv
        #self.gener_idcdv = gener_idcdv
        self.n = len(property_id)

    def __apply_slice__(self, prop: PBUSH1D, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop.k = self.k[i]
        prop.c = self.c[i]
        prop.sa = self.sa[i]
        prop.se = self.se[i]
        prop.mass = self.mass[i]

        prop.spring_type = self.spring_type[i]
        prop.spring_table = self.spring_table[i, :]
        prop.spring_equation = self.spring_equation[i, :]

        prop.damper_type = self.damper_type[i]
        prop.damper_table = self.damper_table[i, :]
        prop.damper_equation = self.damper_equation[i, :]

        prop.shock_type = self.shock_type[i]
        prop.shock_table = self.shock_table[i, :]
        prop.shock_equation = self.shock_equation[i, :]

        prop.gener_equation = self.gener_equation[i, :]
        prop.n = len(i)

    @staticmethod
    def _read_spring(card, istart: int) -> tuple[str, int, int, int, int]:
        """
        F(u) = Ft(u)
        """
        spring_type = string_or_blank(card, istart + 1, 'springType', default='TABLE')
        spring_idt = integer(card, istart + 2, 'springIDT')

        if spring_type == 'TABLE':
            spring_idc = integer_or_blank(card, istart + 3, 'springIDC', default=spring_idt)
            spring_idtdu = integer_or_blank(card, istart + 4, 'springIDTDU', default=0)
            spring_idcdu = integer_or_blank(card, istart + 5, 'springIDCDU', default=spring_idtdu)
        elif spring_type == 'EQUAT':
            spring_idc = integer_or_blank(card, istart + 3, 'springIDC', default=spring_idt)
            spring_idtdu = integer(card, istart + 4, 'springIDTDU')
            spring_idcdu = integer_or_blank(card, istart + 5, 'springIDCDU', default=spring_idtdu)
        else:
            msg = 'Invalid spring_type=%r on card\n%s' % (spring_type, card)
            raise RuntimeError(msg)

        return spring_type, spring_idt, spring_idc, spring_idtdu, spring_idcdu

    @staticmethod
    def _read_damper(card: BDFCard, istart: int) -> tuple[str, int, int, int, int]:
        """
        F(v) = Ft(u)
        """
        damper_type = string_or_blank(card, istart + 1, 'damperType', default='TABLE')
        damper_idt = integer(card, istart + 2, 'damperIDT')
        if damper_type == 'TABLE':
            damper_idc = integer_or_blank(card, istart + 3, 'damperIDC', default=damper_idt)
            damper_idtdv = integer_or_blank(card, istart + 4, 'damperIDTDV', default=0)
            damper_idcdv = integer_or_blank(card, istart + 5, 'damperIDCDV', default=damper_idtdv)
        elif damper_type == 'EQUAT':
            damper_idc = integer_or_blank(card, istart + 3, 'damperIDC', default=damper_idt)
            damper_idtdv = integer(card, istart + 4, 'damperIDTDV')
            damper_idcdv = integer_or_blank(card, istart + 5, 'damperIDCDV', default=damper_idtdv)
        else:
            msg = 'Invalid damper_type=%r on card\n%s' % (damper_type, card)
            raise RuntimeError(msg)

        return damper_type, damper_idt, damper_idc, damper_idtdv, damper_idcdv

    @staticmethod
    def _read_shock(card: BDFCard, istart: int) -> tuple[int,
                                                         tuple[str, float, float, float, float]]:
        """
        F(u, v) = Cv * S(u) * sign(v) * |v|^ev
        """
        #CVT   Viscous damping coefficient CV for tension v > 0,
        #      force per unit velocity. (Real).
        #CVC   Viscous damping coefficient CV for compression v > 0, force
        #      per unit velocity. (Real); default = CVT
        #EXPVT Exponent of velocity EXPV for tension v > 0. (Real); default = 1.
        #EXPVC Exponent of velocity EXPV for compression v < 0. (Real). EXPVT

        # TYPE=TABLE
        # IDTS   TABLEDi id for tension and compression. Defines the
        #        scale factor S, versus displacement u.
        #
        # Type=EQUAT
        # IDETS DEQATN id for tension. Defines the scale
        #       factor S, versus displacement u, for tension u > 0.
        # IDECS DEQATN id for compression. Defines the scale
        #       factor S, versus displacement u, for compression u < 0; default = IDETS
        # IDETSD DEQATN id for tension. Defines the
        #        derivative of the scale factor S, with respect
        shock_type = string_or_blank(card, istart + 1, 'shockType')
        shock_cvt = double(card, istart + 2, 'shockCVT')
        shock_cvc = double_or_blank(card, istart + 3, 'shockCVC', default=shock_cvt)
        shock_exp_vt = double_or_blank(card, istart + 4, 'shockExpVT', default=1.0)
        shock_exp_vc = double_or_blank(card, istart + 5,
                                       'shockExpVC', default=shock_exp_vt)

        if shock_type == 'TABLE':
            #shock_idts = None
            shock_idets = 0
            shock_idecs = 0
            shock_idetsd = 0
            shock_idecsd = 0
            shock_idts = integer(card, istart + 6, 'shockIDTS')
            #shock_idets = blank(card, istart + 9, 'shockIDETS')
            #shock_idecs = blank(card, istart + 10, 'shockIDECS')
            #shock_idetsd = blank(card, istart + 11, 'shockIDETSD')
            #shock_idecsd = blank(card, istart + 12, 'shockIDECSD')
        elif shock_type == 'EQUAT':
            shock_idts = blank(card, istart + 6, 'shockIDTS')
            shock_idets = integer(card, istart + 9, 'shockIDETS')
            shock_idecs = integer_or_blank(card, istart + 10,
                                           'shockIDECS', default=shock_idets)
            shock_idetsd = integer(card, istart + 11, 'shockIDETSD')
            shock_idecsd = integer_or_blank(card, istart + 11,
                                            'shockIDECSD', default=shock_idetsd)
            #def DEquation(self):
                #if isinstance(self.dequation, int):
                    #return self.dequation
                #return self.dequation.equation_id
        else:
            msg = 'Invalid shockType=%r on card\n%s' % (shock_type, card)
            raise RuntimeError(msg)

        out = (
            shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
            shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
        )
        istart += 8
        return istart, out

    @staticmethod
    def _read_gener(card, istart: int) -> tuple[int,
                                                tuple[int, int, int, int, int, int]]:
        """
        F(u, v) = Ft(u, v)
        """
        gener_idt = integer(card, istart + 2, 'generIDT')
        gener_idc = integer_or_blank(card, istart + 3,
                                     'generIDC', gener_idt)
        gener_idtdu = integer(card, istart + 4, 'generIDTDU')
        gener_idcdu = integer_or_blank(card, istart + 5,
                                       'generIDCDU', gener_idtdu)
        gener_idtdv = integer(card, istart + 6, 'generIDTDV')
        gener_idcdv = integer_or_blank(card, istart + 7,
                                       'generIDCDV', gener_idtdv)
        out = (
            gener_idt, gener_idc,
            gener_idtdu, gener_idcdu,
            gener_idtdv, gener_idcdv,
        )
        return out

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass
        #nid = self.model.grid.node_id
        #pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          #msg=f'no spring properties for {self.type}')
        #pids.sort()
        #geom_check(self,
                   #missing,
                   #node=(nid, self.nodes),
                   #property_id=(pids, self.property_id)
        #)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass
    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        sas = array_float_nan(self.sa, size=size, is_double=False)
        ses = array_float_nan(self.se, size=size, is_double=False)
        for pid, k, c, mass, sa, se in zip(property_id, self.k, self.c,
                                           self.mass, sas, ses):
            list_fields = ['PBUSH1D', pid, k, c, mass, None,
                           sa, se, None]
            self_vars = {}
            for var in self_vars:
                if var == 'SHOCKA':
                    list_fields += self._shock_fields()
                elif var == 'SPRING':
                    list_fields += self._spring_fields()
                elif var == 'DAMPER':
                    list_fields += self._damper_fields()
                elif var == 'GENER':
                    list_fields += self._gener_fields()
                else:
                    msg = f'var={var!r} not supported PBUSH1D field...'
                    raise RuntimeError(msg)
                nspaces = 8 - (len(list_fields) - 1) % 8

                if nspaces < 8:
                    list_fields += [None] * (nspaces)

            bdf_file.write(print_card(list_fields))
        return

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass(self) -> np.ndarray:
        mass = self.mass
        return mass

CBUSH2D = CBUSH
PBUSH2D = PBUSH
