from __future__ import annotations
from itertools import zip_longest
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

# from pyNastran.bdf.field_writer_8 import print_card_8
# from pyNastran.utils.numpy_utils import integer_types # , cast_ints
from pyNastran.bdf.cards.expand_card import expand_thru_by

# from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
# from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer,
    double,
    string,
    # integer_or_blank, double_or_blank,
    integer_or_string,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.breakdowns import NO_MASS
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    hslice_by_idim,
    make_idim,
    VectorizedBaseCard,
    parse_check,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str,
    array_float,
    get_print_card_size,
)  # , array_default_int
from pyNastran.dev.bdf_vectorized3.cards.constraints import ADD
from pyNastran.op2.result_objects.scalar6_table_object import float_types

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    # from pyNastran.dev.bdf_vectorized3.bdf import BDF

AREA_ETYPES = {"CTRIA3", "CQUAD4", "CTRIA6", "CQUAD8", "CQUAD"}
LENGTH_ETYPES = {"CONROD", "CROD", "PROD", "CBAR", "CBEAM"}


class NSMi(VectorizedBaseCard):
    """
    Defines a set of non structural mass.

    +-----+-----+------+----+-------+----+-------+----+-------+
    |  1  |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +=====+=====+======+====+=======+====+=======+====+=======+
    | NSM | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +-----+-----+------+----+-------+----+-------+----+-------+
    """

    _id_name = "nsm_id"

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.nsm_id = np.array([], dtype="int32")
        self.nsm_type = np.array([], dtype="|U4")
        self.pid_eid = np.array([], dtype="int32")
        self.value = np.array([], dtype="int32")

    def add_card(self, card: BDFCard, ifile: int, comment: str = "") -> int:
        """
        Adds an NSM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, "sid")
        nsm_type = string(card, 2, "Type")
        pid_eid = integer(card, 3, "pid/eid")
        value = double(card, 4, "value")
        pid_eids = [pid_eid]
        values = [value]
        # self.cards.append((sid, nsm_type, pid_eid, value, comment))
        self.n += 1
        ifield = 5
        while card.field(ifield):
            pid_eid = integer(card, ifield, "pid/eid")
            value = double(card, ifield + 1, "value")
            pid_eids.append(pid_eid)
            values.append(value)
            ifield += 2

        self.cards.append((sid, nsm_type, pid_eids, values, comment))
        # return cls(sid, nsm_type, pid_eid, value, comment=comment)
        return self.n - 1

    def add(self, sid: int, nsm_type: str, pid_eid: int, value: float, comment: str = "") -> int:
        """
        Creates an NSM card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        pid_eid : list[int]; int
            property id or element id depending on nsm_type
        value : list[float]; float
            the non-structural pass per unit length/area
            same length as pid_eid
        comment : str; default=''
            a comment for the card

        """
        NSM_TYPES = {
            "PSHELL",
            "PCOMP",
            "PSHEAR",
            "PBAR",
            "PBARL",
            "PBEAM",
            "PBEAML",
            "PBCOMP",
            "PBEND",
            "PROD",
            "CONROD",
            "PTUBE",
            "PCONEAX",
            "PRAC2D",
            "ELEMENT",
        }
        assert nsm_type in NSM_TYPES, "nsm_type={nsm_type!r}"
        # if isinstance(pid_eid, integer_types):
        self.cards.append((sid, nsm_type, pid_eid, value, comment))
        self.n += 1
        # else:
        ##for pidi, valuei in zip(pid_eid, value):
        # self.cards.append((sid, nsm_type, pidi, valuei, comment))
        # comment = ''
        # self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        nsm_id = np.zeros(ncards, dtype="int32")
        nsm_type = np.zeros(ncards, dtype="|U7")  # ELEMENT
        nvalues = np.zeros(ncards, dtype="int32")
        pid_eids = []
        values = []
        # pid_eid = np.zeros(ncards, dtype='int32')
        # value = np.zeros(ncards, dtype='float64')
        # I11, I21, I22, I31, I32, I33 = I
        # inertia = np.zeros((ncards, 6), dtype='float64')
        for icard, card in enumerate(self.cards):
            (sid, nsm_typei, pid_eidi, valuei, comment) = card
            nsm_id[icard] = sid
            nsm_type[icard] = nsm_typei
            assert len(nsm_typei) <= 7, f"nsm_type={nsm_typei!r} len={len(nsm_typei)}"
            # model.log.debug(pid_eidi)
            if isinstance(pid_eidi, int):
                pid_eidi = [pid_eidi]
                valuei = [valuei]
            pid_eids.extend(pid_eidi)
            values.extend(valuei)
            nvalues[icard] = len(valuei)
            # pid_eid[icard] = pid_eidi
            # value[icard] = valuei
        pid_eid = np.array(pid_eids, dtype="int32")
        value = np.array(values, dtype="float64")
        self._save(nsm_id, nsm_type, pid_eid, value, nvalues)
        self.sort()
        self.cards = []

    def _save(self, nsm_id, nsm_type, pid_eid, value, nvalue):
        assert len(self.nsm_id) == 0, self.nsm_id
        self.nsm_id = nsm_id
        self.nsm_type = nsm_type
        self.pid_eid = pid_eid
        self.value = value
        self.nvalue = nvalue

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        properties = {
            "PROD",
            "PTUBE",
            "PBAR",
            "PBARL",
            "PBEAM",
            "PBEAML",
            "PSHELL",
            "PSHEAR",
        }
        elements = {
            #'CQUAD4', 'CQUAD8', 'CQUADR',
            #'CTRIA3', 'CTRIA6', 'CTRIAR',
            "ELEMENT",
            "CONROD",
        }
        pids_used = []
        eids_used = []
        for nsm_type, pid_eid in zip_longest(self.nsm_type, self.pid_eid):
            if nsm_type in properties:
                pids_used.append(pid_eid)
            elif nsm_type in elements:
                eids_used.append(pid_eid)
            else:  # pragma: no cover
                raise NotImplementedError(nsm_type)
        if len(pids_used):
            used_dict["property_id"] = np.unique(pids_used)
        if len(eids_used):
            used_dict["element_id"] = np.unique(eids_used)

    @property
    def ivalue(self) -> np.ndarray:
        return make_idim(self.n, self.nvalue)

    def slice_card_by_id(self, ids: np.ndarray | int,
                         assume_sorted: bool=True,
                         sort_ids: bool=False) -> NSMi:
        card = slice_duplicate_card_by_id(
            self, ids, assume_sorted=assume_sorted, sort_ids=sort_ids)
        return card

    def __apply_slice__(self, nsm: NSMi, i: np.ndarray) -> None:
        nsm.nsm_id = self.nsm_id[i]
        nsm.nsm_type = self.nsm_type[i]
        # nsm.pid_eid = self.pid_eid[i]
        # nsm.value = self.value[i]
        nsm.nvalue = self.nvalue[i]
        idim = self.ivalue
        nsm.pid_eid = hslice_by_idim(i, idim, self.pid_eid)
        nsm.value = hslice_by_idim(i, idim, self.value)
        nsm.n = len(i)

    @property
    def ivalue(self) -> np.ndarray:
        return make_idim(self.n, self.nvalue)

    @property
    def max_id(self) -> int:
        return max(self.nsm_id.max(), self.pid_eid.max())

    # def geom_check(self, missing: dict[str, np.ndarray]):
    # nid = self.model.grid.node_id
    # cid = self.model.coord.coord_id
    # ucid = np.unique(self.coord_id)
    # if ucid[0] == -1:
    # ucid = ucid[1:]
    # geom_check(self,
    # missing,
    # node=(nid, self.node_id),
    # coord=(cid, ucid))

    # def mass(self) -> np.ndarray:
    # return self._mass

    # def centroid(self) -> np.ndarray:
    # nid = self.model.grid.node_id
    # xyz = self.model.grid.xyz_cid0()
    # inode = np.searchsorted(nid, self.node_id)
    # assert np.array_equal(nid[inode], self.node_id)
    # centroid = xyz[inode, :] + self.xyz_offset

    ## handle cid=-1
    ##ucid = np.unique(self.coord_id)
    # icoord = np.where(self.coord_id == -1)
    # centroid[icoord, :] = self.xyz_offset[icoord, :]
    # return centroid


class NSM1i(VectorizedBaseCard):
    """
    Defines a set of non structural mass.

    +------+-----+------+-------+-----+----+----+----+----+
    |  1   |  2  |  3   |   4   |  5  | 6  | 7  | 8  | 9  |
    +======+=====+======+=======+=====+====+====+====+====+
    | NSM1 | SID | TYPE | VALUE | ID  | ID | ID | ID | ID |
    +------+-----+------+-------+-----+----+----+----+----+
    |      |  ID |  ID  |  ID   | etc |    |    |    |    |
    +------+-----+------+-------+-----+----+----+----+----+
    """

    _id_name = "nsm_id"

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.nsm_id = np.array([], dtype="int32")
        self.nsm_type = np.array([], dtype="|U7")
        self.pid_eid = np.array([], dtype="int32")
        self.npid_eid = np.array([], dtype="int32")
        self.value = np.array([], dtype="int32")

    def add(self, sid: int, nsm_type: str, value: float, ids: list[int], comment: str = "") -> int:
        """
        Creates an NSM1/NSML1 card, which defines non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        value : float
            the non-structural pass per unit length/area
        ids : list[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card

        """
        # assert isinstance(ids, (list, np.ndarray)), ids
        assert isinstance(value, float), value
        self.cards.append((sid, nsm_type, ids, value, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str = "") -> int:
        """
        Adds a NSM1/NSML1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, "sid")
        nsm_type = string(card, 2, "Type")
        value = double(card, 3, "value")

        # TODO: doesn't support 1 THRU 11 BY 2
        ids = []
        # _id = 1
        nfields = len(card)
        if nfields == 5:
            id1 = integer_or_string(card, 4, "ID_1")
            if id1 != "ALL" and not isinstance(id1, int):
                msg = "*ID_1 = %r (field #4) on card must be an integer or ALL.\ncard=%s" % (
                    id1,
                    card,
                )
                raise SyntaxError(msg)
            ids = id1
        else:
            # we'll handle expansion in the init
            ids = card[4:]
        # return cls(sid, nsm_type, value, ids, comment=comment)
        self.cards.append((sid, nsm_type, ids, value, comment))
        self.n += 1

        # return cls(sid, nsm_type, pid_eid, value, comment=comment)
        assert len(card) >= 5, f"len(NSM1 card) = {len(card):d}\ncard={card}"
        return self.n - 1

    def slice_card_by_id(self, ids: np.ndarray | int,
                         assume_sorted: bool=True,
                         sort_ids: bool=False) -> NSM1i:
        card = slice_duplicate_card_by_id(
            self, ids, assume_sorted=assume_sorted, sort_ids=sort_ids)
        return card

    def __apply_slice__(self, elem: NSM1i, i: np.ndarray) -> None:
        elem.nsm_id = self.nsm_id[i]
        elem.nsm_type = self.nsm_type[i]
        elem.value = self.value[i]

        ielement = self.ielement  # [i, :]
        elem.pid_eid = hslice_by_idim(i, ielement, self.pid_eid)
        elem.npid_eid = self.npid_eid[i]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        properties = {
            "PROD",
            "PTUBE",
            "PBAR",
            "PBARL",
            "PBEAM",
            "PBEAML",
            "PSHELL",
            "PSHEAR",
        }
        elements = {
            #'CQUAD4', 'CQUAD8', 'CQUADR',
            #'CTRIA3', 'CTRIA6', 'CTRIAR',
            "ELEMENT",
            "CONROD",
        }
        pids_used = []
        eids_used = []
        insm = self.ielement
        for nsm_type, (insm0, insm1) in zip_longest(self.nsm_type, insm):
            ids = self.pid_eid[insm0:insm1].tolist()
            if nsm_type in properties:
                pids_used.extend(ids)
            elif nsm_type in elements:
                eids_used.extend(ids)
            else:
                raise NotImplementedError(nsm_type)
        if len(pids_used):
            used_dict["property_id"] = np.unique(pids_used)
        if len(eids_used):
            used_dict["element_id"] = np.unique(eids_used)

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.npid_eid)

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        nsm_id = np.zeros(ncards, dtype="int32")

        # PSHELL, ELEMENT
        nsm_type = np.zeros(ncards, dtype="|U7")
        pid_eid_list = []
        npid_eid = np.zeros(ncards, dtype="int32")
        value = np.zeros(ncards, dtype="float64")

        for icard, card in enumerate(self.cards):
            (sid, nsm_typei, pid_eidi, valuei, comment) = card
            nsm_id[icard] = sid
            nsm_type[icard] = nsm_typei
            assert len(nsm_typei) <= 7, f"nsm_type={nsm_typei!r} len={len(nsm_typei)}"
            if isinstance(pid_eidi, int):
                pid_eidi = [pid_eidi]
            elif isinstance(pid_eidi, str):  #  pragma: no cover
                if pid_eidi == "ALL":
                    pid_eidi = [-1]
                else:
                    raise RuntimeError(pid_eidi)
            elif isinstance(pid_eidi, list):
                pid_eidi = expand_thru_by(
                    pid_eidi,
                    set_fields=True,
                    sort_fields=True,
                    require_int=True,
                    allow_blanks=False,
                )
                assert "THRU" not in pid_eidi
            elif isinstance(pid_eidi, np.ndarray):
                pass
            else:
                raise NotImplementedError(pid_eidi)
            pid_eid_list.extend(pid_eidi)
            npid_eid[icard] = len(pid_eidi)
            value[icard] = valuei
        pid_eid = np.array(pid_eid_list, dtype="int32")
        self._save(nsm_id, nsm_type, pid_eid, npid_eid, value)
        self.sort()
        self.cards = []

    def _save(self, nsm_id, nsm_type, pid_eid, npid_eid, value) -> None:
        assert len(self.nsm_id) == 0, self.nsm_id
        self.nsm_id = nsm_id
        self.nsm_type = nsm_type
        self.pid_eid = pid_eid
        self.npid_eid = npid_eid
        self.value = value

    @property
    def max_id(self) -> int:
        return max(self.nsm_id.max(), self.pid_eid.max())

    # def geom_check(self, missing: dict[str, np.ndarray]):
    # nid = self.model.grid.node_id
    # cid = self.model.coord.coord_id
    # ucid = np.unique(self.coord_id)
    # if ucid[0] == -1:
    # ucid = ucid[1:]
    # geom_check(self,
    # missing,
    # node=(nid, self.node_id),
    # coord=(cid, ucid))

    # def mass(self) -> np.ndarray:
    # return self._mass

    # def centroid(self) -> np.ndarray:
    # nid = self.model.grid.node_id
    # xyz = self.model.grid.xyz_cid0()
    # inode = np.searchsorted(nid, self.node_id)
    # assert np.array_equal(nid[inode], self.node_id)
    # centroid = xyz[inode, :] + self.xyz_offset

    ## handle cid=-1
    ##ucid = np.unique(self.coord_id)
    # icoord = np.where(self.coord_id == -1)
    # centroid[icoord, :] = self.xyz_offset[icoord, :]
    # return centroid


class NSM1(NSM1i):
    def inertia(
        self,
        element_id: np.ndarray | list[int] | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        self.model.log.info("NSM1 - inertia")
        element_id, mass, centroid, inertia = inertia1_func(self, element_id=element_id)
        return element_id, mass, centroid, inertia

    @parse_check
    def write_file(
        self,
        bdf_file: TextIOLike,
        size: int = 8,
        is_double: bool = False,
        write_card_headers: bool = False,
    ) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        nsm_str = array_str(self.nsm_id, size=size)
        # pid_eid_str = array_str(self.pid_eid, size=size)
        # insm = self.insm
        ielement = self.ielement
        for nsm_id, nsm_type, value, (insm0, insm1) in zip_longest(
            nsm_str, self.nsm_type, self.value, ielement
        ):
            ids = self.pid_eid[insm0:insm1].tolist()
            list_fields = [
                "NSM1",
                nsm_id,
                nsm_type,
                value,
            ] + ids
            bdf_file.write(print_card(list_fields))
        return


class NSML1(NSM1i):
    def inertia(
        self,
        element_id: np.ndarray | list[int] | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        element_id, mass, centroid, inertia = inertia1_func(self, element_id=element_id)
        return element_id, mass, centroid, inertia

    @parse_check
    def write_file(
        self,
        bdf_file: TextIOLike,
        size: int = 8,
        is_double: bool = False,
        write_card_headers: bool = False,
    ) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        nsm_str = array_str(self.nsm_id, size=size)
        # pid_eid_str = array_str(self.pid_eid, size=size)
        insm = self.ielement
        for nsm_id, nsm_type, value, (insm0, insm1) in zip_longest(
            nsm_str, self.nsm_type, self.value, insm
        ):
            ids = self.pid_eid[insm0:insm1].tolist()
            list_fields = [
                "NSML1",
                nsm_id,
                nsm_type,
                value,
            ] + ids
            bdf_file.write(print_card(list_fields))
        return


def get_pids_list(cards: list, ids) -> np.ndarray:
    pids_list = []
    for card in cards:
        pids_list.append(card.property_id)
    pids = np.hstack(pids_list)
    if len(ids) == 1 and ids[0] != -1:
        pids = np.intersect1d(pids, ids)
    assert len(pids) > 0, pids
    return pids


def _compute_inertia(mass: np.ndarray, centroid: np.ndarray) -> np.ndarray:
    """Compute [Ixx, Iyy, Izz, Ixy, Ixz, Iyz] from point masses at centroids."""
    dx = centroid[:, 0]
    dy = centroid[:, 1]
    dz = centroid[:, 2]
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz
    inertia = np.empty((len(mass), 6), dtype=centroid.dtype)
    inertia[:, 0] = mass * (dy2 + dz2)
    inertia[:, 1] = mass * (dx2 + dz2)
    inertia[:, 2] = mass * (dx2 + dy2)
    inertia[:, 3] = mass * (dx * dy)
    inertia[:, 4] = mass * (dx * dz)
    inertia[:, 5] = mass * (dy * dz)
    return inertia


def _nsml_property_vectorized(
    ids: np.ndarray,
    value: np.ndarray,
    valid_pids: np.ndarray,
    element_cards: list,
    is_area: bool,
    is_length: bool,
    elementi_id_list: list,
    centroidi_list: list,
    mass_area_length_list: list,
    total_area_length_list: list,
) -> tuple[bool, bool]:
    """Batch all PIDs: get elements once, apply per-PID values via searchsorted."""
    pids_common = np.intersect1d(ids, valid_pids)
    if len(pids_common) == 0:
        return is_area, is_length

    all_eids = []
    all_pids_elem = []
    all_area = []
    all_centroid = []
    for card in element_cards:
        if len(card) == 0:
            continue
        pids_in_card = np.intersect1d(pids_common, card.property_id)
        if len(pids_in_card) == 0:
            continue
        card2 = card.slice_card_by_property_id(pids_in_card)
        all_eids.append(card2.element_id)
        all_pids_elem.append(card2.property_id)
        all_centroid.append(card2.centroid())
        if card2.type in {"CROD", "CTUBE", "CBAR", "CBEAM", "CBEND"}:
            is_length = True
            all_area.append(card2.length())
        else:
            is_area = True
            all_area.append(card2.area())

    if len(all_eids) == 0:
        return is_area, is_length

    elem_ids = np.hstack(all_eids)
    elem_pids = np.hstack(all_pids_elem)
    area = np.hstack(all_area)
    centroid = np.vstack(all_centroid)

    # map per-PID values to each element via its property_id
    isort = np.argsort(ids)
    ids_sorted = ids[isort]
    value_sorted = value[isort]
    idx = np.searchsorted(ids_sorted, elem_pids)
    elem_value = value_sorted[idx]

    mass = area * elem_value
    elementi_id_list.append(elem_ids)
    centroidi_list.append(centroid)
    mass_area_length_list.append(mass)
    total_area_length_list.append(area)
    return is_area, is_length


def _nsml_property(
    nsm_id: int,
    value: float | np.floating,
    pids: np.ndarray,
    element_cards: list,
    is_area: bool,
    is_length: bool,
    elementi_id_list,
    centroidi_list,
    mass_area_length_list,
    total_area_length_list,
    divide_by_sum: bool,
) -> tuple[bool, bool]:
    assert isinstance(value, (float, np.floating)), value

    element_cards2 = []
    for card in element_cards:
        if len(card) == 0:
            continue
        pids_common = np.intersect1d(pids, card.property_id)
        if len(pids_common) == 0:
            continue
        card2 = card.slice_card_by_property_id(pids_common)
        element_cards2.append(card2)
    assert len(element_cards2) > 0, element_cards2

    card0 = element_cards2[0]
    if card0.type in {"CROD", "CTUBE", "CBAR", "CBEAM", "CBEND"}:
        is_length = True
    elif card0.type in {"CTRIA3", "CTRIA6", "CTRIAR", "CQUAD4", "CQUAD8", "CQUADR"}:
        is_area = True
    else:  # pragma: no cover
        raise NotImplementedError(card0.type)

    element_id_list = []
    area_list = []
    centroid_list = []
    for card in element_cards2:
        element_id_list.append(card.element_id)
        centroid_list.append(card.centroid())
        if is_length:
            area_list.append(card.length())
        else:
            area_list.append(card.area())

    element_id = np.hstack(element_id_list)
    area = np.hstack(area_list)
    centroid = np.vstack(centroid_list)
    mass = area * value

    elementi_id_list.append(element_id)
    centroidi_list.append(centroid)
    mass_area_length_list.append(mass)
    total_area_length_list.append(area)
    return is_area, is_length


def _nsml_element(
    value,
    ids: np.ndarray,
    card,
    is_area: bool,
    is_length: bool,
    elementi_id_list,
    centroidi_list,
    mass_area_length_list,
    total_area_length_list,
    divide_by_sum: bool,
) -> tuple[bool, bool]:
    if len(ids) == 1 and ids[0] == -1:  # ALL
        eid_common = card.element_id
    else:
        assert ids.min() > 0, ids
        eid_common = np.intersect1d(ids, card.element_id)
    if len(eid_common) == 0:
        return is_area, is_length

    if card.type in AREA_ETYPES:
        card2 = card.slice_card_by_id(eid_common)
        area = card2.area()
        is_area = True
    elif card.type in LENGTH_ETYPES:
        card2 = card.slice_card_by_id(eid_common)
        area = card2.length()
        is_length = True
    else:
        raise NotImplementedError(card.type)
    centroid = card2.centroid()

    if isinstance(value, np.ndarray) and len(value) > 1:
        icommon = np.searchsorted(ids, eid_common)
        value_filtered = value[icommon]
        mass_area_length = area * value_filtered
    else:
        v = value if isinstance(value, (float, np.floating)) else value[0]
        mass_area_length = area * v

    elementi_id_list.append(eid_common)
    centroidi_list.append(centroid)
    mass_area_length_list.append(mass_area_length)
    total_area_length_list.append(area)
    return is_area, is_length


class NSM(NSMi):
    def write_file(
        self,
        bdf_file: TextIOLike,
        size: int = 8,
        is_double: bool = False,
        write_card_headers: bool = False,
    ) -> None:
        if len(self.nsm_id) == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        values = array_float(self.value, size=size, is_double=False)
        for nsm_id, nsm_type, (ivalue0, ivalue1) in zip_longest(
            nsm_str, self.nsm_type, self.ivalue
        ):
            pid_eid_value = alternate_values_in_list(
                pid_eid_str[ivalue0:ivalue1], values[ivalue0:ivalue1]
            )
            # pid_eid = pid_eid_str[ivalue0:ivalue1]
            # value = self.value[ivalue0:ivalue1]
            list_fields = ["NSM", nsm_id, nsm_type] + pid_eid_value
            bdf_file.write(print_card(list_fields))
        return

    def inertia(
        self,
        element_id: np.ndarray | list[int] | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        element_id, mass, centroid, inertia = inertia_func(self, element_id=element_id)
        return element_id, mass, centroid, inertia


class NSML(NSMi):
    def inertia(
        self,
        element_id: np.ndarray | list[int] | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        element_id, mass, centroid, inertia = inertia_func(self, element_id=element_id)
        return element_id, mass, centroid, inertia

    @parse_check
    def write_file(
        self,
        bdf_file: TextIOLike,
        size: int = 8,
        is_double: bool = False,
        write_card_headers: bool = False,
    ) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        values = array_float(self.value, size=size, is_double=False)
        for nsm_id, nsm_type, (ivalue0, ivalue1) in zip_longest(
            nsm_str, self.nsm_type, self.ivalue
        ):
            pid_eid_value = alternate_values_in_list(
                pid_eid_str[ivalue0:ivalue1], values[ivalue0:ivalue1]
            )
            list_fields = ["NSML", nsm_id, nsm_type] + pid_eid_value
            bdf_file.write(print_card(list_fields))
        return


NSMs = NSM | NSM1 | NSML | NSML1


def alternate_values_in_list(pid_eid_str, values):
    pid_eid_value = []
    for pid_eid, value in zip(pid_eid_str, values):
        pid_eid_value.append(pid_eid)
        pid_eid_value.append(value)
    return pid_eid_value


class NSMADD(ADD):
    """
    Defines an NSM combination set as a union of NSM cards.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | NSMADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """

    _id_name = "nsm_id"

    @property
    def nsm_id(self):
        return self.sid

    @nsm_id.setter
    def nsm_id(self, nsm_id: np.ndarray):
        self.sid = nsm_id

    @property
    def nsm_ids(self):
        return self.sids

    @nsm_ids.setter
    def nsm_ids(self, nsm_ids: np.ndarray):
        self.sids = nsm_ids

    @property
    def nnsm(self):
        return self.nsids

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict["nsm_id"] = self.nsm_ids

    def write_file(
        self,
        bdf_file: TextIOLike,
        size: int = 8,
        is_double: bool = False,
        write_card_headers: bool = False,
    ) -> None:
        if len(self.nsm_id) == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)
        nsm_id = array_str(self.nsm_id, size=size)
        nsm_ids = array_str(self.nsm_ids, size=size)
        for nsm_id, idim in zip(nsm_id, self.idim):
            idim0, idim1 = idim
            assert idim1 > idim0, self.idim
            nsm_idsi = nsm_ids[idim0:idim1].tolist()
            assert len(nsm_idsi) > 0, self.idim
            list_fields = ["NSMADD", nsm_id] + nsm_idsi
            bdf_file.write(print_card(list_fields))
        return

    def get_nsms_by_nsm_id(self) -> dict[int, NSMs]:
        model = self.model
        """"""
        # unsm_ids = np.unique(self.nsm_id)
        nsm_by_nsm_id = defaultdict(list)
        # for nsm_id in unsm_ids:
        # nsm_by_nsm_id[nsm_id] = []

        for nsm in model.nsm_cards:
            if nsm.type in {"NSMADD"}:
                continue

            unsm_idsi = np.unique(nsm.nsm_id)
            for unsm_id in unsm_idsi:
                if unsm_id not in nsm.nsm_id:
                    continue

                i = np.where(unsm_id == nsm.nsm_id)[0]
                if len(i) == 0:
                    continue
                # print(nsm.type, unsm_id, nsm.nsm_id, i)
                self._ids
                nsmi = nsm.slice_card_by_index(i)
                nsm_by_nsm_id[unsm_id].append(nsmi)

        if self.n > 0:
            # build the nsmadds
            unsm_ids_list = []
            for nsm_id, (idim0, idim1) in zip(self.nsm_id, self.idim):
                nsm_ids = self.nsm_ids[idim0:idim1]
                unsm_ids_list.append(nsm_ids)
            unsm_ids = np.hstack(unsm_ids_list)
            for nsm_idi in unsm_ids:
                nsm_by_nsm_id[nsm_id].extend(nsm_by_nsm_id[nsm_idi])
        return dict(nsm_by_nsm_id)

    # def slice_card_by_nsm_id(self, nsm_id: np.ndarray) -> NSMADD:
    # assert self.n > 0, self.n
    # assert len(self.element_id) > 0, self.element_id
    # i = self.index(element_id)
    ##cls_obj = cls(self.model)
    ##cls_obj.__apply_slice__(self, i)
    # cls_obj = self.slice_card_by_index(i)
    # assert cls_obj.n > 0, cls_obj
    # return cls_obj

    def get_reduced_nsms(
        self,
        # resolve_load_card: bool=False,
        stop_on_failure: bool = True,
    ) -> dict[int, NSMs]:
        """
        Parameters
        ----------
        resolve_load_card : bool; default=False
            ???
        """
        stop_on_failure = True
        nsm_by_nsm_id = self.get_nsms_by_nsm_id()
        log = self.model.log

        reduced_nsms = {}
        for sid, idim in zip(self.nsm_id, self.idim):
            reduced_nsmsi = []
            idim0, idim1 = idim

            nsm_ids = self.nsm_ids[idim0:idim1]
            for nsm_idi in nsm_ids:
                nsms_found = nsm_by_nsm_id[nsm_idi]
                if len(nsms_found) == 0:
                    msg = f"No referenced NSMs found for nsm_id={nsm_idi} on NSMADD nsm_id={sid}"
                    log.error(msg)
                    if stop_on_failure:
                        raise RuntimeError(msg)
                reduced_nsmsi.append(nsms_found)
            if sid in reduced_nsms:
                reduced_nsms[sid].extend(reduced_nsmsi)
            else:
                reduced_nsms[sid] = reduced_nsmsi
        return reduced_nsms

    @property
    def _ids(self) -> np.ndarray:
        return getattr(self, self._id_name)

    @_ids.setter
    def _ids(self, ids: np.ndarray) -> None:
        return setattr(self, self._id_name, ids)

    def __apply_slice__(self, nsm: NSMADD, i: np.ndarray) -> None:
        nsm.n = len(i)
        nsm.sid = self.sid[i]
        nsm.nsm_id = self.nsm_id[i]
        nsm.nsids = self.nsids[i]
        # nsm.components = self.components[i]

        # nsm.nnodes = self.nnodes[i]
        idim = self.idim
        nsm.sids = hslice_by_idim(i, idim, self.sids)
        return nsm

    def slice_card_by_id(self, ids: np.ndarray | int,
                         assume_sorted: bool=True,
                         sort_ids: bool=False) -> NSMADD:
        card = slice_duplicate_card_by_id(
            self, ids, assume_sorted=assume_sorted, sort_ids=sort_ids)
        return card

    # def __apply_slice__(self, spc: NSMADD, i: np.ndarray) -> None:
    # spc.n = len(i)
    # spc.spc_id = self.spc_id[i]
    # spc.components = self.components[i]

    # spc.nnodes = self.nnodes[i]
    # idim = self.inode
    # spc.node_id = hslice_by_idim(i, idim, self.node_id)
    # return spc


def slice_duplicate_card_by_id(self,
                               ids: np.ndarray | int,
                               assume_sorted: bool = True,
                               sort_ids: bool = False):
    ids = np.atleast_1d(np.asarray(ids, dtype=self.nsm_id.dtype))
    i = np.where(np.isin(self.nsm_id, ids))[0]
    return self.slice_card_by_index(i)


def inertia_func(
    nsm_obj: NSM | NSML,
    element_id: np.ndarray | list[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    divide_by_sum = False
    if nsm_obj.type in ["NSML1", "NSML"]:
        divide_by_sum = True

    element_id_list = []
    mass_list = []
    centroid_list = []
    inertia_list = []
    model = nsm_obj.model

    for nsm_id, nsm_type, (ivalue0, ivalue1) in zip_longest(
        nsm_obj.nsm_id, nsm_obj.nsm_type, nsm_obj.ivalue
    ):
        # log.debug(f'nsm_id={nsm_id} nsm_type={nsm_type}')
        ids = nsm_obj.pid_eid[ivalue0:ivalue1]
        value = nsm_obj.value[ivalue0:ivalue1]
        assert len(ids) == len(value), value
        ids.sort()
        is_area = False
        is_length = False

        elementi_id_list = []
        centroidi_list = []
        # inertiai_list = []
        mass_area_length_list = []
        total_area_length_list = []
        assert ids.min() > 0, ids
        if nsm_type == "CONROD":
            card = model.conrod
            is_area, is_length = _nsml_element(
                value,
                ids,
                card,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )
            assert is_length, (is_area, is_length)
        elif nsm_type == "PSHELL":
            cards = [
                model.ctria3,
                model.ctriar,
                model.ctria6,
                model.cquad4,
                model.cquadr,
                model.cquad8,
                model.cquad,
            ]
            is_area, is_length = _nsml_property_vectorized(
                ids,
                value,
                model.pshell.property_id,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
            )

        elif nsm_type == "PCOMP":
            cards = [
                model.ctria3,
                model.ctriar,
                model.ctria6,
                model.cquad4,
                model.cquadr,
                model.cquad8,
                model.cquad,
            ]
            pcomp_pids = np.hstack([model.pcomp.property_id, model.pcompg.property_id])
            is_area, is_length = _nsml_property_vectorized(
                ids,
                value,
                pcomp_pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
            )

        elif nsm_type in {"PROD", "PTUBE"}:
            cards = [model.crod, model.ctube]
            prod_pids = np.hstack([model.prod.property_id, model.ptube.property_id])
            is_area, is_length = _nsml_property_vectorized(
                ids,
                value,
                prod_pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
            )

        elif nsm_type in {"PBAR", "PBARL"}:
            cards = [model.cbar]
            pbar_pids = np.hstack([model.pbar.property_id, model.pbarl.property_id])
            is_area, is_length = _nsml_property_vectorized(
                ids,
                value,
                pbar_pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
            )

        elif nsm_type in {"PBEAM", "PBEAML"}:
            cards = [model.cbeam]
            pbeam_pids = np.hstack([model.pbeam.property_id, model.pbeaml.property_id])
            is_area, is_length = _nsml_property_vectorized(
                ids,
                value,
                pbeam_pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
            )

        elif nsm_type == "ELEMENT":
            for card in model.element_cards:
                if len(card) == 0 or card.type in NO_MASS:
                    continue
                is_area, is_length = _nsml_element(
                    value,
                    ids,
                    card,
                    is_area,
                    is_length,
                    elementi_id_list,
                    centroidi_list,
                    mass_area_length_list,
                    total_area_length_list,
                    divide_by_sum,
                )
        else:
            raise NotImplementedError(nsm_type)

        if len(elementi_id_list) == 0:
            continue

        if is_area and is_length:
            raise RuntimeError("Area and length are both used")
        assert is_area or is_length, (is_area, is_length)

        elementi_id = np.hstack(elementi_id_list)
        centroidi = np.vstack(centroidi_list)
        mass_area_length = np.hstack(mass_area_length_list)
        if divide_by_sum:
            if nsm_type == "ELEMENT":
                total_area_length = np.hstack(total_area_length_list)
                massi = mass_area_length / total_area_length
            else:
                total_area_length = np.hstack(total_area_length_list).sum()
                massi = mass_area_length / total_area_length
            del mass_area_length, total_area_length
        else:
            massi = mass_area_length
            del mass_area_length

        neidi = len(elementi_id)
        if neidi == 0:
            continue
        assert len(elementi_id) == neidi, (massi, neidi)
        assert centroidi.shape == (neidi, 3), f"neidi={neidi}, centroidi.shape={centroidi.shape}"

        if element_id is not None:
            ikeep = np.isin(elementi_id, element_id)
            elementi_id = elementi_id[ikeep]
            centroidi = centroidi[ikeep]
            massi = massi[ikeep]
            neidi = len(elementi_id)
            assert centroidi.shape == (neidi, 3), (
                f"neidi={neidi}, centroidi.shape={centroidi.shape}"
            )

        if len(elementi_id) == 0:
            continue
        element_id_list.append(elementi_id)
        mass_list.append(massi)
        centroid_list.append(centroidi)
        inertia_list.append(_compute_inertia(massi, centroidi))

    if len(element_id_list) == 0:
        empty_int = np.array([], dtype="int32")
        empty_float = np.array([], dtype="float64")
        return (
            empty_int,
            empty_float,
            np.zeros((0, 3), dtype="float64"),
            np.zeros((0, 6), dtype="float64"),
        )

    element_id_out = np.hstack(element_id_list)
    mass_out = np.hstack(mass_list)
    centroid_out = np.vstack(centroid_list)
    inertia_out = np.vstack(inertia_list)

    neidi = len(element_id_out)
    assert len(mass_out) == neidi, (mass_out, neidi)
    assert centroid_out.shape == (neidi, 3), f"neidi={neidi}, centroid.shape={centroid_out.shape}"
    assert inertia_out.shape == (neidi, 6), f"neidi={neidi}, inertia.shape={inertia_out.shape}"
    return element_id_out, mass_out, centroid_out, inertia_out


def inertia1_func(
    nsm_obj: NSM1 | NSML1,
    element_id: np.ndarray | list[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    divide_by_sum = False
    if nsm_obj.type in ["NSML1", "NSML"]:
        divide_by_sum = True

    element_id_list = []
    mass_list = []
    centroid_list = []
    inertia_list = []
    model = nsm_obj.model

    log = model.log
    log.debug(f"  {nsm_obj.type}: divide_by_mass={divide_by_sum}")
    for nsm_id, nsm_type, value, (insm0, insm1) in zip_longest(
        nsm_obj.nsm_id, nsm_obj.nsm_type, nsm_obj.value, nsm_obj.ielement
    ):
        log.debug(f"  nsm_id={nsm_id} nsm_type={nsm_type}")
        ids = nsm_obj.pid_eid[insm0:insm1]
        assert isinstance(value, float_types), value
        ids.sort()
        is_area = False
        is_length = False

        elementi_id_list = []
        centroidi_list = []
        # inertiai_list = []
        mass_area_length_list = []
        total_area_length_list = []
        # assert ids.min() > 0, (ids, nsm_type)
        if nsm_type == "CONROD":
            card = model.conrod
            log.debug(f"ids = {ids}")
            if len(ids) == 1 and ids[0] == -1:
                ids = card.element_id
            is_area, is_length = _nsml_element(
                value,
                ids,
                card,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )
            assert is_length, (is_area, is_length)

        elif nsm_type in {"PROD", "PTUBE"}:
            pids = get_pids_list([model.prod, model.ptube], ids)
            cards = [model.crod, model.ctube]

            is_area, is_length = _nsml_property(
                nsm_id,
                value,
                pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )
        elif nsm_type in {"PBAR", "PBARL"}:
            pids = get_pids_list([model.pbar, model.pbarl], ids)
            cards = [model.cbar]
            is_area, is_length = _nsml_property(
                nsm_id,
                value,
                pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )
        elif nsm_type in {"PBEAM", "PBEAML"}:
            pids = get_pids_list([model.pbeam, model.pbeaml], ids)
            cards = [model.cbeam]
            is_area, is_length = _nsml_property(
                nsm_id,
                value,
                pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )

        elif nsm_type == "PSHELL":
            pids = get_pids_list([model.pshell], ids)
            cards = [model.ctria3, model.cquad4, model.ctria6, model.cquad8, model.cquad]
            is_area, is_length = _nsml_property(
                nsm_id,
                value,
                pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )

        elif nsm_type == "PCOMP":
            pids = get_pids_list([model.pcomp, model.pcompg], ids)
            cards = [model.ctria3, model.cquad4, model.ctria6, model.cquad8, model.cquad]
            is_area, is_length = _nsml_property(
                nsm_id,
                value,
                pids,
                cards,
                is_area,
                is_length,
                elementi_id_list,
                centroidi_list,
                mass_area_length_list,
                total_area_length_list,
                divide_by_sum,
            )

        elif nsm_type == "ELEMENT":
            for card in model.element_cards:
                if len(card) == 0 or card.type in NO_MASS:
                    continue
                is_area, is_length = _nsml_element(
                    value,
                    ids,
                    card,
                    is_area,
                    is_length,
                    elementi_id_list,
                    centroidi_list,
                    mass_area_length_list,
                    total_area_length_list,
                    divide_by_sum,
                )
        else:  # pragma: no cover
            raise NotImplementedError(nsm_type)

        if len(elementi_id_list) == 0:
            continue

        if is_area and is_length:
            raise RuntimeError("Area and length are are both used")
        assert is_area or is_length, (is_area, is_length)

        assert len(elementi_id_list) == len(centroidi_list), (
            nsm_type,
            len(elementi_id_list),
            len(centroidi_list),
        )
        elementi_id = np.hstack(elementi_id_list)
        centroidi = np.vstack(centroidi_list)
        mass_area_length = np.hstack(mass_area_length_list)
        neidi = len(elementi_id)
        assert centroidi.shape == (neidi, 3), f"neidi={neidi}, centroidi.shape={centroidi.shape}"

        if divide_by_sum:
            total_area_length = np.hstack(total_area_length_list).sum()
            massi = mass_area_length / total_area_length
            del mass_area_length, total_area_length
        else:
            # assert total_area_length_list[0] is None, total_area_length_list
            massi = mass_area_length
            del mass_area_length

        if element_id is not None:
            ikeep = np.isin(elementi_id, element_id)
            elementi_id = elementi_id[ikeep]
            centroidi = centroidi[ikeep]
            massi = massi[ikeep]
            neidi = len(elementi_id)
            assert centroidi.shape == (neidi, 3), (
                f"neidi={neidi}, centroidi.shape={centroidi.shape}"
            )

        if len(elementi_id) == 0:
            continue
        element_id_list.append(elementi_id)
        mass_list.append(massi)
        centroid_list.append(centroidi)
        inertia_list.append(_compute_inertia(massi, centroidi))

    if len(element_id_list) == 0:
        empty_int = np.array([], dtype="int32")
        empty_float = np.array([], dtype="float64")
        return (
            empty_int,
            empty_float,
            np.zeros((0, 3), dtype="float64"),
            np.zeros((0, 6), dtype="float64"),
        )

    element_id_out = np.hstack(element_id_list)
    mass_out = np.hstack(mass_list)
    centroid_out = np.vstack(centroid_list)
    inertia_out = np.vstack(inertia_list)

    neidi = len(element_id_out)
    assert len(mass_out) == neidi, (mass_out, neidi)
    assert centroid_out.shape == (neidi, 3), f"neidi={neidi}, centroid.shape={centroid_out.shape}"
    assert inertia_out.shape == (neidi, 6), f"neidi={neidi}, inertia.shape={inertia_out.shape}"
    return element_id_out, mass_out, centroid_out, inertia_out
