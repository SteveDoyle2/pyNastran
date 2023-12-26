from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

#from pyNastran.bdf.field_writer_8 import set_blank_if_default, set_string8_blank_if_default
from pyNastran.bdf.cards.base_card import expand_thru, expand_thru_by # BaseCard, expand_thru_by #  _node_ids,
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_or_string, fields,
    #integer_types,
    float_types)
from pyNastran.bdf.cards.utils import wipe_empty_fields

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, hslice_by_idim, make_idim,
    parse_load_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int, get_print_card_size
from .static_loads import Load

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class QHBDY(Load):
    """Defines a uniform heat flux into a set of grid points."""
    flag_to_nnodes = {
        'POINT' : (1, 1),
        'LINE' : (2, 2),
        'REV' : (2, 2),
        'AREA3' : (3, 3),
        'AREA4' : (4, 4),
        'AREA6' : (4, 6), # 4-6
        'AREA8' : (5, 8), # 5-8
    }
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> PLOAD1:
        #load = PLOAD1(self.model)
        #load.n = len(i)
        #load.load_id = self.load_id[i]
        #load.element_id = self.element_id[i]
        #load.scale = self.scale[i]
        #load.x = self.x[i, :]
        #load.pressure = self.pressure[i, :]
        #return load

    def add(self, sid: int, flag: int, q0: float, grids: list[int],
            area_factor: Optional[float]=None, comment: str='') -> int:
        """
        Creates a QHBDY card

        Parameters
        ----------
        sid : int
            load id
        flag : str
            valid_flags = {POINT, LINE, REV, AREA3, AREA4, AREA6, AREA8}
        q0 : float
            Magnitude of thermal flux into face. Q0 is positive for heat
            into the surface
        area_factor : float; default=None
            Area factor depends on type
        grids : list[int]
            Grid point identification of connected grid points
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, flag, q0, area_factor, grids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a QHBDY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'eid')
        flag = string(card, 2, 'flag')

        q0 = double(card, 3, 'q0')
        af = double_or_blank(card, 4, 'af')
        nnodes_required, nnodes_max = self.flag_to_nnodes[flag]

        grids = []
        if nnodes_required == nnodes_max:
            for i in range(nnodes_required):
                grid = integer(card, 5 + i, 'grid%i' % (i + 1))
                grids.append(grid)
        else:
            int_node_count = 0
            for i in range(nnodes_max):
                grid = integer_or_blank(card, 5 + i, 'grid%i' % (i + 1))
                if grid is not None:
                    int_node_count += 1
                grids.append(grid)
            if int_node_count < nnodes_required:
                msg = 'int_node_count=%s nnodes_required=%s' % (int_node_count, nnodes_required)
                raise RuntimeError(msg)

        self.cards.append((sid, flag, q0, af, grids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #: Load set identification number. (Integer > 0)
        load_id = np.zeros(ncards, dtype='int32')

        #: Magnitude of thermal flux into face. Q0 is positive for heat
        #: into the surface. (Real)
        q0 = np.zeros(ncards, dtype='float64')

        #: Area factor depends on type. (Real > 0.0 or blank)
        area_factor = np.zeros(ncards, dtype='float64')

        #: Grid point identification of connected grid points.
        #: (Integer > 0 or blank)
        ngrid = np.zeros(ncards, dtype='int32')
        all_grids = []

        #assert flag in ['POINT', 'LINE', 'REV', 'AREA3', 'AREA4', 'AREA6', 'AREA8'], str(self)
        flag = np.zeros(ncards, dtype='|U5')

        assert ncards > 0, ncards

        for icard, card in enumerate(self.cards):
            (sid, flagi, q0i, af, gridsi, comment) = card
            load_id[icard] = sid
            q0[icard] = q0i
            flag[icard] = flagi
            area_factor[icard] = af
            ngrid[icard] = len(gridsi)
            all_grids.extend(gridsi)

        grids = np.array(all_grids, dtype='int32')
        self._save(load_id, q0, flag, area_factor, grids, ngrid)
        self.cards = []

    def _save(self, load_id, q0, flag, area_factor, grids, ngrid):
        assert len(self.load_id) == 0, self.load_id
        nloads = len(load_id)
        self.load_id = load_id
        self.q0 = q0
        self.flag = flag
        self.area_factor = area_factor
        self.grids = grids
        self.ngrid = ngrid
        self.n = nloads

    #def slice_card_by_index(self, i: np.ndarray) -> QHBDY:
        #load = QHBDY(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: QHBDY, i: np.ndarray) -> None:
        load.load_id = self.load_id[i]
        load.q0 = self.q0[i]
        load.area_factor = self.area_factor[i]
        load.grids = hslice_by_idim(i, self.inode, self.grids)
        load.ngrid = self.ngrid[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.grids), filter_node0=False,
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        grids = self.grids[self.grids > 0]
        used_dict['node_id'].append(grids)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.grids
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def inode(self) -> np.ndarray:
        return make_idim(self.n, self.ngrid)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.grids.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        load_ids = array_str(self.load_id, size=size)
        for sid, flag, q0, af, inid in zip(load_ids, self.flag, self.q0, self.area_factor, self.inode):
            inid0, inid1 = inid
            grids = self.grids[inid0:inid1].tolist()
            #assert len(grids) > 0, grids
            #eids.sort()
            list_fields = ['QHBDY', sid, flag, q0, af] + grids
            bdf_file.write(print_card(list_fields))
        return


class QBDY1(VectorizedBaseCard):
    """
    Defines a uniform heat flux into CHBDYj elements.

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    def add(self, sid: int, qflux: float, eids: list[int],
            comment: str='') -> int:
        """Creates a QBDY1 card"""
        self.cards.append((sid, qflux, eids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a QBDY1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        qflux = double(card, 2, 'qflux')
        eids = []
        j = 1
        for i in range(3, len(card)):
            eid = integer_or_string(card, i, 'eid%d' % j)
            eids.append(eid)
            j += 1
        self.cards.append((sid, qflux, eids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #: Load set identification number. (Integer > 0)
        load_id = np.zeros(ncards, dtype='int32')

        qflux = np.zeros(ncards, dtype='float64')
        nelement = np.zeros(ncards, dtype='int32')

        assert ncards > 0, ncards
        elements = []
        all_elements = []
        for icard, card in enumerate(self.cards):
            (sid, qfluxi, elementsii, comment) = card

            load_id[icard] = sid
            qflux[icard] = qfluxi
            elementsi = expand_thru(elementsii)
            nelement[icard] = len(elementsi)
            all_elements.extend(elementsi)

        #: CHBDYj element identification numbers
        elements = np.array(all_elements, dtype='int32')
        self._save(load_id, qflux, elements, nelement)
        self.cards = []

    def _save(self, load_id, qflux, elements, nelement):
        assert len(self.load_id) == 0, self.load_id
        nloads = len(load_id)
        self.load_id = load_id
        self.qflux = qflux
        self.elements = elements
        self.nelement = nelement
        self.n = nloads

    #def slice_card_by_index(self, i: np.ndarray) -> QBDY1:
        #load = QBDY1(self.model)
        #self.__apply_slice__(load, i)
        #self.n = len(o)
        #return load

    def __apply_slice__(self, load: QBDY1, i: np.ndarray) -> None:
        load.load_id = self.load_id[i]
        load.qflux = self.qflux[i]
        load.elements = hslice_by_idim(i, self.ielement, self.elements)
        load.nelement = self.nelement[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        element_id = self.model.shell_element_ids
        geom_check(self,
                   missing,
                   element_id=(element_id, self.elements),
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['element_id'].append(self.elements)

    @property
    def ielement(self) -> np.ndarray:
        #print('ielement =', self.nelement, self.elements)
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.elements.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        load_ids = array_str(self.load_id, size=size)
        for sid, qflux, ieid in zip(load_ids, self.qflux, self.ielement):
            ieid0, ieid1 = ieid
            eids = self.elements[ieid0:ieid1].tolist()
            assert len(eids) > 0, eids
            #eids.sort()
            list_fields = ['QBDY1', sid, qflux] + eids
            bdf_file.write(print_card(list_fields))
        return


class QBDY2(VectorizedBaseCard):
    """
    Defines a uniform heat flux load for a boundary surface.

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    def add(self, sid: int, eid: int, qfluxs: list[float],
            comment: str='') -> int:
        """Creates a QBDY1 card"""
        if isinstance(qfluxs, float_types):
            qfluxs = [qfluxs]
        self.cards.append((sid, eid, qfluxs, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a QBDY2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')

        qfluxs = []
        j = 1
        for i in range(3, len(card)):
            q = double_or_blank(card, i, 'qFlux%d' % j)
            qfluxs.append(q)
            j += 1

        assert len(qfluxs) > 0, qfluxs
        qfluxs = wipe_empty_fields(qfluxs)
        self.cards.append((sid, eid, qfluxs, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #: Load set identification number. (Integer > 0)
        load_id = np.zeros(ncards, dtype='int32')

        all_qfluxs = []
        nflux = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype='int32')
        assert ncards > 0, ncards

        for icard, card in enumerate(self.cards):
            (sid, eid, qfluxs, comment) = card
            load_id[icard] = sid
            element_id[icard] = eid
            all_qfluxs.extend(qfluxs)
            nflux[icard] = len(qfluxs)
        qflux = np.array(all_qfluxs, dtype='float64')
        self._save(load_id, element_id, qflux, nflux)
        self.cards = []

    def _save(self, load_id, element_id, qflux, nflux):
        assert len(self.load_id) == 0, self.load_id
        nloads = len(load_id)
        self.load_id = load_id
        self.element_id = element_id
        self.qflux = qflux
        self.nflux = nflux
        self.n = nloads

    def slice_card_by_index(self, i: np.ndarray) -> QBDY2:
        load = QBDY2(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: QBDY2, i: np.ndarray) -> None:
        load.load_id = self.load_id[i]
        load.element_id = self.element_id[i]
        load.qflux = hslice_by_idim(i, self.iflux, self.qflux)
        load.nflux = self.nflux[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        element_id = self.model.shell_element_ids
        geom_check(self,
                   missing,
                   element_id=(element_id, self.element_id), filter_node0=False,
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['element_id'].append(self.element_id)

    @property
    def iflux(self) -> np.ndarray:
        return make_idim(self.n, self.nflux)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.element_id.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        load_ids = array_str(self.load_id, size=size)
        element_ids = array_str(self.element_id, size=size)
        for sid, eid, (iflux0, iflux1) in zip(load_ids, element_ids, self.iflux):
            qflux = self.qflux[iflux0:iflux1].tolist()
            assert len(qflux) > 0, qflux
            #eids.sort()
            list_fields = ['QBDY2', sid, eid, ] + qflux
            bdf_file.write(print_card(list_fields))
        return


class QBDY3(Load):
    """
    Defines a uniform heat flux load for a boundary surface.

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')
        self.q0 = np.array([], dtype='float64')
        self.control_node = np.array([], dtype='int32')
        self.elements = np.array([], dtype='int32')
        self.nelement = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> PLOAD1:
        #load = PLOAD1(self.model)
        #load.n = len(i)
        #load.load_id = self.load_id[i]
        #load.element_id = self.element_id[i]
        #load.scale = self.scale[i]
        #load.x = self.x[i, :]
        #load.pressure = self.pressure[i, :]
        #return load

    def add(self, sid: int, q0, cntrlnd: int, eids: list[int],
            comment: str='') -> int:
        """
        Creates a QBDY3 card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        q0 : float; default=None
            Magnitude of thermal flux vector into face
        control_id : int; default=0
            Control point
        eids : list[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, q0, cntrlnd, eids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        q0 = double(card, 2, 'q0')
        cntrlnd = integer_or_blank(card, 3, 'cntrlnd', default=0)

        nfields = card.nfields
        eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
        eids = expand_thru_by(eids)

        self.cards.append((sid, q0, cntrlnd, eids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #: Load set identification number. (Integer > 0)
        load_id = np.zeros(ncards, dtype='int32')

        #: Heat flux into element
        q0 = np.zeros(ncards, dtype='float64')

        #: Control point for thermal flux load. (Integer > 0; Default = 0)
        control_node = np.zeros(ncards, dtype='int32')

        nelement = np.zeros(ncards, dtype='int32')
        assert ncards > 0, ncards
        elements_list = []
        for icard, card in enumerate(self.cards):
            (sid, q0i, cntrlnd, eids, comment) = card
            nelementi = len(eids)

            load_id[icard] = sid
            q0[icard] = q0i
            control_node[icard] = cntrlnd
            nelement[icard] = nelementi
            elements_list.extend(eids)

        #: CHBDYj element identification numbers
        elements = np.array(elements_list, dtype='int32')
        self._save(load_id, q0, control_node, elements, nelement)
        self.cards = []

    def _save(self, load_id, q0, control_node, elements, nelement):
        assert len(self.load_id) == 0, self.load_id
        nloads = len(load_id)
        self.load_id = load_id
        self.q0 = q0
        self.control_node = control_node
        self.elements = elements
        self.nelement = nelement
        self.n = nloads

    #def slice_card_by_index(self, i: np.ndarray) -> QBDY3:
        #load = QBDY3(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: QBDY3, i: np.ndarray) -> None:
        load.load_id = self.load_id[i]
        load.q0 = self.q0[i]
        load.control_node = self.control_node[i]
        load.elements = hslice_by_idim(i, self.ielement, self.elements)
        load.nelement = self.nelement[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.control_node), filter_node0=False,
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.control_node)
        used_dict['element_id'].append(self.elements)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.control_node
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.control_node.max(), self.elements.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        load_ids = array_str(self.load_id, size=size)
        control_nodes = array_default_int(self.control_node, default=0, size=size)
        for sid, q0, control_node, ieid in zip(load_ids, self.q0, control_nodes, self.ielement):
            ieid0, ieid1 = ieid
            eids = self.elements[ieid0:ieid1]
            assert len(eids) > 0, eids
            eids.sort()
            list_fields = ['QBDY3', sid, q0, control_node] + collapse_thru_by(eids)
            bdf_file.write(print_card(list_fields))
        return


class QVOL(Load):
    """
    Defines a rate of volumetric heat addition in a conduction element.

    +------+------+------+---------+------+------+------+------+------+
    |  1   |   2  |   3  |    4    |   5  |   6  |   7  |   8  |   9  |
    +======+======+======+=========+======+======+======+======+======+
    | QVOL | SID  | QVOL | CNTRLND | EID1 | EID2 | EID3 | EID4 | EID5 |
    +------+------+------+---------+------+------+------+------+------+
    |      | EID6 | etc. |         |      |      |      |      |      |
    +------+------+------+---------+------+------+------+------+------+

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')
        self.qvol = np.array([], dtype='float64')
        self.control_node = np.array([], dtype='int32')
        self.elements = np.array([], dtype='int32')
        self.nelement = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> PLOAD1:
        #load = PLOAD1(self.model)
        #load.n = len(i)
        #load.load_id = self.load_id[i]
        #load.element_id = self.element_id[i]
        #load.scale = self.scale[i]
        #load.x = self.x[i, :]
        #load.pressure = self.pressure[i, :]
        #return load

    def add(self, sid: int, qvol: float,
            control_point: int, elements: list[int],
            comment: str='') -> int:
        """Creates a QVOL card"""
        self.cards.append((sid, qvol, control_point, elements, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        qvol = double(card, 2, 'qvol')
        control_point = integer_or_blank(card, 3, 'control_id', default=0)

        i = 1
        eids = []
        for ifield in range(4, len(card)):
            eid = integer_or_string(card, ifield, 'eid_%d' % i)
            eids.append(eid)
            i += 1
        elements = expand_thru_by(eids)
        self.cards.append((sid, qvol, control_point, elements, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #: Load set identification number. (Integer > 0)
        load_id = np.zeros(ncards, dtype='int32')

        #: volumetric heat load
        qvol = np.zeros(ncards, dtype='float64')

        #: Control point for thermal flux load. (Integer > 0; Default = 0)
        control_node = np.zeros(ncards, dtype='int32')

        nelement = np.zeros(ncards, dtype='int32')
        #self.ielement = np.zeros((ncards, 2), dtype='int32')

        elements = []

        assert ncards > 0, ncards
        all_elements = []
        for icard, card in enumerate(self.cards):
            (sid, qvoli, control_nodei, elements, comment) = card
            nelementi = len(elements)

            load_id[icard] = sid
            qvol[icard] = qvoli
            control_node[icard] = control_nodei
            nelement[icard] = nelementi
            all_elements.extend(elements)

        #: CHBDYj element identification numbers
        elements = np.array(all_elements, dtype='int32')
        self._save(load_id, qvol, control_node, elements, nelement)
        self.cards = []

    def _save(self, load_id, qvol, control_node, elements, nelement):
        nloads = len(load_id)
        self.load_id = load_id
        self.qvol = qvol
        self.control_node = control_node
        self.elements = elements
        assert elements is not None
        self.nelement = nelement
        assert nelement.dtype.name in {'int32', 'int64'}, nelement.dtype.name
        self.n = nloads

    #def slice_card_by_index(self, i: np.ndarray) -> QVOL:
        #load = QVOL(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: QVOL, i: np.ndarray) -> None:
        load.load_id = self.load_id[i]
        load.qvol = self.qvol[i]
        load.control_node = self.control_node[i]
        load.elements = hslice_by_idim(i, self.ielement, self.elements)
        load.nelement = self.nelement[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        element_id = self.model.shell_element_ids
        geom_check(self,
                   missing,
                   node=(nid, self.control_node), filter_node0=False,
                   element_id=(element_id, self.elements),
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.control_node)
        used_dict['element_id'].append(self.elements)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.control_node
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.control_node.max(), self.elements.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        load_ids = array_str(self.load_id, size=size)
        control_nodes = array_default_int(self.control_node, default=0, size=size)
        #elements = array_str(self.elements, size=size)
        for sid, qvol, control_node, ieid in zip(load_ids, self.qvol, control_nodes, self.ielement):
            ieid0, ieid1 = ieid
            eids = self.elements[ieid0:ieid1].tolist()
            assert len(eids) > 0, eids
            eids.sort()
            list_fields = ['QVOL', sid, qvol, control_node] + collapse_thru_by(eids)
            bdf_file.write(print_card(list_fields))
        return


class TEMPBC(VectorizedBaseCard):
    _id_name = 'spc_id'
    def clear(self) -> None:
        self.n = 0
        self.spc_id = np.array([], dtype='int32')
        #self.control_node = np.array([], dtype='int32')
        self.temperature = np.array([], dtype='float64')
        self.nodes = np.array([], dtype='int32')
        self.nnode = np.zeros((0, 2), dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> TEMPBC:
        #load = TEMPBC(self.model)
        #load.n = len(i)
        #load.load_id = self.load_id[i]
        #load.element_id = self.element_id[i]
        #load.scale = self.scale[i]
        #load.x = self.x[i, :]
        #load.pressure = self.pressure[i, :]
        #return load

    def add_card(self, card: BDFCard, comment: str='') -> None:
        sid = integer(card, 1, 'sid')
        bc_type = string_or_blank(card, 2, 'Type', default='STAT')
        nfields_left = len(card) - 3
        assert nfields_left > 0, card
        assert nfields_left % 2 == 0, card

        temps = []
        nodes = []
        for i in range(nfields_left // 2):
            ifield = 3 + i*2
            temp = double(card, ifield, 'temp_%d'%  ((i+1)))
            nid = integer(card, ifield+1, 'temp_%d'%  ((i+1)))
            temps.append(temp)
            nodes.append(nid)
        self.cards.append((sid, bc_type, temps, nodes, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        spc_id = np.zeros(ncards, dtype='int32')
        bc_type = np.zeros(ncards, dtype='|U4')
        #control_node = np.zeros(ncards, dtype='int32')
        nnode = np.zeros(ncards, dtype='int32')

        assert ncards > 0, ncards

        all_temps = []
        all_nodes = []
        for icard, card in enumerate(self.cards):
            (sid, bc_type, temps, nodes, comment) = card
            spc_id[icard] = sid
            bc_type[icard] = bc_type
            #control_node[icard] = control_point
            nnode[icard] = len(nodes)
            all_nodes.extend(nodes)
            all_temps.extend(temps)

        nodes = np.array(all_nodes, dtype='int32')
        temperature = np.array(all_temps, dtype='float64')
        self._save(spc_id, bc_type, nnode, nodes, temperature)
        self.cards = []

    #def slice_card_by_index(self, i: np.ndarray) -> TEMPBC:
        #load = TEMPBC(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: TEMPBC, i: np.ndarray) -> None:
        load.spc_id = self.spc_id[i]
        load.bc_type = self.bc_type[i]
        #load.control_node = self.control_node[i]
        load.temperature = hslice_by_idim(i, self.inode, self.temperature)
        load.nodes = hslice_by_idim(i, self.inode, self.nodes)
        load.nnode = self.nnode[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.nodes), filter_node0=False,
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def inode(self) -> np.ndarray:
        return make_idim(self.n, self.nnode)

    @property
    def max_id(self) -> int:
        return max(self.spc_id.max(), self.nodes.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spc_id) == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        spc_ids = array_str(self.spc_id, size=size)
        nodes_ = array_str(self.nodes, size=size)
        for sid, bc_type, inode in zip(spc_ids, self.bc_type, self.inode):
            inode0, inode1 = inode
            nodes = nodes_[inode0:inode1]
            temps = self.temperature[inode0:inode1]

            list_fields = ['TEMPBC', sid, bc_type]
            for temp, node in zip(temps, nodes):
                list_fields.extend([temp, node])
            bdf_file.write(print_card(list_fields))
        return


class RADM(VectorizedBaseCard):
    """
    Defines the radiation properties of a boundary element for heat transfer
    analysis

    """
    _id_name = 'rad_mid'
    def clear(self) -> None:
        self.n = 0
        self.rad_mid = np.array([], dtype='int32')
        self.absorptivity = np.array([], dtype='float64')
        self.nemissivity = np.array([], dtype='int32')
        self.emissivity = np.array([], dtype='float64')

    def add(self, radmid: int, absorb: float, emissivity: list[float], comment: str='') -> int:
        if isinstance(emissivity, float_types):
            emissivity = [emissivity]
        self.cards.append((radmid, absorb, emissivity, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card, comment='') -> int:
        """
        Adds a RADM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nfields = card.nfields
        radmid = integer(card, 1, 'radmid')
        absorb = double_or_blank(card, 2, 'absorb')
        emissivity = fields(double, card, 'emissivity', i=3, j=nfields)
        assert isinstance(emissivity, list), emissivity
        self.cards.append((radmid, absorb, emissivity, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        rad_mid = np.zeros(ncards, dtype='int32')
        absorptivity = np.zeros(ncards, dtype='float64')
        nemissivity = np.zeros(ncards, dtype='int32')

        assert ncards > 0, ncards

        all_emissivity = []
        for icard, card in enumerate(self.cards):
            (radmid, absorb, emissivity, comment) = card
            rad_mid[icard] = radmid
            absorptivity[icard] = absorb
            nemissivity[icard] = len(emissivity)
            all_emissivity.extend(emissivity)
        emissivity = np.array(all_emissivity, dtype='float64')
        self._save(rad_mid, absorptivity, nemissivity, emissivity)
        self.cards = []

    def _save(self, rad_mid, absorptivity, nemissivity, emissivity):
        self.rad_mid = rad_mid
        self.absorptivity = absorptivity
        self.nemissivity = nemissivity
        self.emissivity = emissivity

    def __apply_slice__(self, load: QHBDY, i: np.ndarray) -> None:
        load.rad_mid = self.rad_mid[i]
        load.absorptivity = self.absorptivity[i]
        load.emissivity = hslice_by_idim(i, self.iemissivity, self.emissivity)
        load.nemissivity = self.nemissivity[i]
        load.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    @property
    def iemissivity(self) -> np.ndarray:
        return make_idim(self.n, self.nemissivity)

    @property
    def max_id(self) -> int:
        return self.rad_mid.max()

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.rad_mid) == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)
        rad_mids = array_str(self.rad_mid, size=size)
        for rad_mid, absorb, (iemissivity0, iemissivity1) in zip(
            rad_mids, self.absorptivity, self.iemissivity):
            emissivity = self.emissivity[iemissivity0:iemissivity1].tolist()
            list_fields = ['RADM', rad_mid, absorb] + emissivity
            bdf_file.write(print_card(list_fields))
        return


class RADBC(VectorizedBaseCard):
    """
    Specifies an CHBDYi element face for application of radiation boundary
    conditions

    """
    _id_name = 'node_id'
    def clear(self) -> None:
        self.n = 0
        self.node_id = np.array([], dtype='int32')
        self.factor_ambient = np.array([], dtype='float64')
        self.control_node = np.array([], dtype='int32')
        self.elements = np.array([], dtype='int32')
        self.nelement = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> RADBC:
        #load = RADBC(self.model)
        #load.n = len(i)
        #load.load_id = self.load_id[i]
        #load.element_id = self.element_id[i]
        #load.scale = self.scale[i]
        #load.x = self.x[i, :]
        #load.pressure = self.pressure[i, :]
        #return load

    def add(self, node_amb, famb, control_node, eids, comment='') -> int:
        assert len(eids) > 0, eids
        self.cards.append((node_amb, famb, control_node, eids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        node_amb = integer(card, 1, 'nodamb')
        famb = double(card, 2, 'famb')
        control_node = integer_or_blank(card, 3, 'cntrlnd', default=0)
        nfields = card.nfields
        eids = fields(integer_or_string, card, 'eid', i=4, j=nfields)
        eids_expand = expand_thru_by(eids)
        self.cards.append((node_amb, famb, control_node, eids_expand, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')
        node_id = np.zeros(ncards, dtype='int32')
        factor_ambient = np.zeros(ncards, dtype='float64')
        control_node = np.zeros(ncards, dtype='int32')
        nelement = np.zeros(ncards, dtype='int32')

        assert ncards > 0, ncards

        all_elements = []
        for icard, card in enumerate(self.cards):
            (node_amb, famb, control_nodei, eids, comment) = card
            node_id[icard] = node_amb
            factor_ambient[icard] = famb
            control_node[icard] = control_nodei
            nelement[icard] = len(eids)
            all_elements.extend(eids)

        elements = np.array(all_elements, dtype='int32')
        self._save(node_id, factor_ambient, control_node, nelement, elements)
        self.cards = []

    def _save(self, node_id, factor_ambient, control_node, nelement, elements):
        self.node_id = node_id
        self.factor_ambient = factor_ambient
        self.control_node = control_node
        self.nelement = nelement
        self.elements = elements

    def slice_card_by_index(self, i: np.ndarray) -> RADBC:
        load = RADBC(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: RADBC, i: np.ndarray) -> None:
        load.node_id = self.node_id[i]
        load.factor_ambient = self.factor_ambient[i]
        load.control_node = self.control_node[i]
        load.elements = hslice_by_idim(i, self.ielement, self.elements)
        load.nelement = self.nelement[i]
        load.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        element_id = self.model.shell_element_ids
        card_nodes = np.hstack([self.node_id, self.control_node])
        geom_check(self,
                   missing,
                   node=(nid, card_nodes), filter_node0=False,
                   element_id=(element_id, self.elements),
                   )

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)
        used_dict['node_id'].append(self.control_node)
        used_dict['element_id'].append(self.elements)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.control_node
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return max(self.node_id.max(), self.control_node.max(), self.elements.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.node_id) == 0:
            return ''
        print_card, size = get_print_card_size(size, self.max_id)
        node_ids = array_str(self.node_id, size=size)
        control_nodes = array_default_int(self.control_node, default=0, size=size)
        for nid, famb, (ielement0, ielement1), control_nodei in zip(node_ids, self.factor_ambient, self.ielement, control_nodes):
            eids = self.elements[ielement0:ielement1].tolist()
            eids2 = collapse_thru_by(eids)
            list_fields = ['RADBC', nid, famb, control_nodei] + eids2

            bdf_file.write(print_card(list_fields))
        return
