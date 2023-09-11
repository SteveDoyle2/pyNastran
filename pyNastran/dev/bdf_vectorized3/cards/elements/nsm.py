from __future__ import annotations
from itertools import zip_longest
from collections import defaultdict
from typing import Union, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.numpy_utils import integer_types, cast_ints
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    integer_or_blank, double_or_blank,
    integer_or_string,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import Element, get_print_card_8_16
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.cards.constraints import ADD
if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class NSMi(Element):
    """
    Defines a set of non structural mass.

    +-----+-----+------+----+-------+----+-------+----+-------+
    |  1  |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +=====+=====+======+====+=======+====+=======+====+=======+
    | NSM | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +-----+-----+------+----+-------+----+-------+----+-------+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nsm_id = np.array([], dtype='int32')
        self.nsm_type = np.array([], dtype='|U4')
        self.pid_eid = np.array([], dtype='int32')
        self.value = np.array([], dtype='int32')

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
        sid = integer(card, 1, 'sid')
        nsm_type = string(card, 2, 'Type')
        pid_eid = integer(card, 3, 'pid/eid')
        value = double(card, 4, 'value')
        self.cards.append((sid, nsm_type, pid_eid, value, comment))
        self.n += 1
        if card.field(5):
            pid_eid = integer(card, 5, 'pid/eid')
            value = double(card, 6, 'value')
            self.cards.append((sid, nsm_type, pid_eid, value, comment))
            self.n += 1
        if card.field(7):
            pid_eid = integer(card, 7, 'pid/eid')
            value = double(card, 8, 'value')
            self.cards.append((sid, nsm_type, pid_eid, value, comment))
            self.n += 1

        #return cls(sid, nsm_type, pid_eid, value, comment=comment)
        assert len(card) <= 9, f'len(NSM card) = {len(card):d}\ncard={card}'
        return self.n

    def add(self, sid: int, nsm_type: str, pid_eid: int, value: float,
            comment: str='') -> int:
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
            'PSHELL', 'PCOMP', 'PSHEAR',
            'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP', 'PBEND',
            'PROD', 'CONROD', 'PTUBE', 'PCONEAX', 'PRAC2D',
            'ELEMENT'}
        assert nsm_type in NSM_TYPES, 'nsm_type={nsm_type!r}'
        if isinstance(pid_eid, integer_types):
            self.cards.append((sid, nsm_type, pid_eid, value, comment))
            self.n += 1
        else:
            for pidi, valuei in zip(pid_eid, value):
                self.cards.append((sid, nsm_type, pidi, valuei, comment))
                comment = ''
                self.n += 1
        return self.n

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        nsm_id = np.zeros(ncards, dtype='int32')
        nsm_type = np.zeros(ncards, dtype='|U7') # ELEMENT
        pid_eid = np.zeros(ncards, dtype='int32')
        value = np.zeros(ncards, dtype='float64')
        #I11, I21, I22, I31, I32, I33 = I
        #inertia = np.zeros((ncards, 6), dtype='float64')
        for icard, card in enumerate(self.cards):
            (sid, nsm_typei, pid_eidi, valuei, comment) = card
            nsm_id[icard] = sid
            nsm_type[icard] = nsm_typei
            assert len(nsm_typei) <= 7, f'nsm_type={nsm_typei!r} len={len(nsm_typei)}'
            print(pid_eidi)
            pid_eid[icard] = pid_eidi
            value[icard] = valuei
        self._save(nsm_id, nsm_type, pid_eid, value)
        self.sort()
        self.cards = []

    def _save(self, nsm_id, nsm_type, pid_eid, value):
        assert len(self.nsm_id) == 0, self.nsm_id
        self.nsm_id = nsm_id
        self.nsm_type = nsm_type
        self.pid_eid = pid_eid
        self.value = value

    def __apply_slice__(self, elem: NSMi, i: np.ndarray) -> None:
        elem.nsm_id = self.nsm_id[i]
        elem.nsm_type = self.nsm_type[i]
        elem.pid_eid = self.pid_eid[i]
        elem.value = self.value[i]
        elem.n = len(i)

    #def geom_check(self, missing: dict[str, np.ndarray]):
        #nid = self.model.grid.node_id
        #cid = self.model.coord.coord_id
        #ucid = np.unique(self.coord_id)
        #if ucid[0] == -1:
            #ucid = ucid[1:]
        #geom_check(self,
                   #missing,
                   #node=(nid, self.node_id),
                   #coord=(cid, ucid))

    #def mass(self) -> np.ndarray:
        #return self._mass

    #def centroid(self) -> np.ndarray:
        #nid = self.model.grid.node_id
        #xyz = self.model.grid.xyz_cid0()
        #inode = np.searchsorted(nid, self.node_id)
        #assert np.array_equal(nid[inode], self.node_id)
        #centroid = xyz[inode, :] + self.xyz_offset

        ## handle cid=-1
        ##ucid = np.unique(self.coord_id)
        #icoord = np.where(self.coord_id == -1)
        #centroid[icoord, :] = self.xyz_offset[icoord, :]
        #return centroid


class NSM1i(Element):
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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nsm_id = np.array([], dtype='int32')
        self.nsm_type = np.array([], dtype='|U4')
        self.pid_eid = np.array([], dtype='int32')
        self.value = np.array([], dtype='int32')

    def add(self, sid: int, nsm_type: str, value: float, ids: list[int],
            comment: str='') -> int:
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
        self.cards.append((sid, nsm_type, ids, value, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a NSM1/NSML1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nsm_type = string(card, 2, 'Type')
        value = double(card, 3, 'value')

        # TODO: doesn't support 1 THRU 11 BY 2
        ids = []
        #_id = 1
        nfields = len(card)
        if nfields == 5:
            id1 = integer_or_string(card, 4, 'ID_1')
            if id1 != 'ALL' and not isinstance(id1, int):
                msg = ('*ID_1 = %r (field #4) on card must be an integer or ALL.\n'
                       'card=%s' % (id1, card))
                raise SyntaxError(msg)
            ids = id1
        else:
            # we'll handle expansion in the init
            ids = card[4:]
        #return cls(sid, nsm_type, value, ids, comment=comment)
        self.cards.append((sid, nsm_type, ids, value, comment))
        self.n += 1

        #return cls(sid, nsm_type, pid_eid, value, comment=comment)
        assert len(card) <= 9, f'len(NSM1 card) = {len(card):d}\ncard={card}'
        return self.n

    def __apply_slice__(self, elem: NSMi, i: np.ndarray) -> None:
        elem.nsm_id = self.nsm_id[i]
        elem.nsm_type = self.nsm_type[i]
        elem.pid_eid = self.pid_eid[i]
        elem.value = self.value[i]
        elem.n = len(i)

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        nsm_id = np.zeros(ncards, dtype='int32')
        nsm_type = np.zeros(ncards, dtype='|U4')
        pid_eid_list = []
        npid_eid = np.zeros(ncards, dtype='int32')
        value = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (sid, nsm_typei, pid_eidi, valuei, comment) = card
            nsm_id[icard] = sid
            nsm_type[icard] = nsm_typei
            assert len(nsm_typei) <= 4, nsm_typei
            pid_eid_list.extend(pid_eidi)
            npid_eid[icard] = len(pid_eidi)
            value[icard] = valuei
        pid_eid = np.array(pid_eid_list, dtype='int32')
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

    #def geom_check(self, missing: dict[str, np.ndarray]):
        #nid = self.model.grid.node_id
        #cid = self.model.coord.coord_id
        #ucid = np.unique(self.coord_id)
        #if ucid[0] == -1:
            #ucid = ucid[1:]
        #geom_check(self,
                   #missing,
                   #node=(nid, self.node_id),
                   #coord=(cid, ucid))

    #def mass(self) -> np.ndarray:
        #return self._mass

    #def centroid(self) -> np.ndarray:
        #nid = self.model.grid.node_id
        #xyz = self.model.grid.xyz_cid0()
        #inode = np.searchsorted(nid, self.node_id)
        #assert np.array_equal(nid[inode], self.node_id)
        #centroid = xyz[inode, :] + self.xyz_offset

        ## handle cid=-1
        ##ucid = np.unique(self.coord_id)
        #icoord = np.where(self.coord_id == -1)
        #centroid[icoord, :] = self.xyz_offset[icoord, :]
        #return centroid

class NSM1(NSM1i):
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        for nsm_id, nsm_type, pid_eid, value, (insm0, insm1) in zip_longest(
                nsm_str, self.nsm_type, pid_eid_str, self.value, self.insm):
            ids = self.pid_eid[insm0:insm1]
            list_fields = ['NSM1', nsm_id, nsm_type, value, ] + ids
            bdf_file.write(print_card(list_fields))
        return

class NSML1(NSM1i):
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        for nsm_id, nsm_type, pid_eid, value, (insm0, insm1) in zip_longest(
                nsm_str, self.nsm_type, pid_eid_str, self.value, self.insm):
            ids = self.pid_eid[insm0:insm1]
            list_fields = ['NSML1', nsm_id, nsm_type, value, ] + ids
            bdf_file.write(print_card(list_fields))
        return

class NSM(NSMi):
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        for nsm_id, nsm_type, pid_eid, value in zip_longest(nsm_str, self.nsm_type, pid_eid_str, self.value):
            list_fields = (['NSM', nsm_id, nsm_type, pid_eid, value])
            bdf_file.write(print_card(list_fields))
        return

class NSML(NSMi):
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        for nsm_id, nsm_type, pid_eid, value in zip_longest(nsm_str, self.nsm_type, pid_eid_str, self.value):
            list_fields = (['NSML', nsm_id, nsm_type, pid_eid, value])
            bdf_file.write(print_card(list_fields))
        return


NSMs = Union[NSM, NSM1,
             NSML, NSML1]

class NSMADD(ADD):
    """
    Defines an NSM combination set as a union of NSM cards.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | NSMADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    #def __init__(self, model: BDF):
        #super().__init__(model)
        #self.sid = np.array([], dtype='int32')
        #self.spc_ids = np.array([], dtype='int32')
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
    def spc_ids(self, nsm_ids: np.ndarray):
        self.sids = nsm_ids

    @property
    def nnsm(self):
        return self.nsids

    def write(self, size: int=8) -> str:
        if len(self.nsm_id) == 0:
            return ''
        lines = []
        nsm_id = array_str(self.nsm_id, size=size)
        nsm_ids = array_str(self.nsm_ids, size=size)
        for nsm_id, idim in zip(nsm_id, self.idim):
            idim0, idim1 = idim
            assert idim1 > idim0, self.idim
            nsm_idsi = nsm_ids[idim0:idim1].tolist()
            assert len(nsm_idsi) > 0, self.idim
            list_fields = ['NSMADD', nsm_id] + nsm_idsi
            lines.append(print_card_8(list_fields))
        return ''.join(lines)

    def get_nsms_by_nsm_id(self) -> dict[int, NSMs]:
        model = self.model
        """"""
        #unsm_ids = np.unique(self.nsm_id)
        nsm_by_nsm_id = defaultdict(list)
        #for nsm_id in unsm_ids:
            #nsm_by_nsm_id[nsm_id] = []

        for nsm in model.nsms:
            if nsm.type in {'NSMADD'}:
                continue

            unsm_idsi = np.unique(nsm.nsm_id)
            for unsm_id in unsm_idsi:
                i = np.where(unsm_id == nsm.nsm_id)[0]
                if len(i) == 0:
                    continue
                nsmi = nsm.slice_card_by_index(i)
                nsm_by_nsm_id[unsm_id].append(nsmi)
        return dict(nsm_by_nsm_id)

    def get_reduced_nsms(self,
                         #resolve_load_card: bool=False,
                         stop_on_failure: bool=True) -> dict[int, SPCs]:
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
                    msg = f'No referenced NSMs found for nsm_id={nsm_idi} on NSMADD nsm_id={sid}'
                    log.error(msg)
                    if stop_on_failure:
                        raise RuntimeError(msg)
                reduced_nsmsi.append(nsms_found)
            reduced_nsms[sid] = reduced_nsmsi
        return reduced_nsms


