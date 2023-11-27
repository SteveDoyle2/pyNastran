from __future__ import annotations
from itertools import zip_longest
from collections import defaultdict
from typing import Union, TYPE_CHECKING
import numpy as np
#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.numpy_utils import integer_types # , cast_ints
from pyNastran.bdf.cards.expand_card import expand_thru_by
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    #integer_or_blank, double_or_blank,
    integer_or_string,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import get_print_card_8_16, hslice_by_idim, make_idim, VectorizedBaseCard
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str # , array_default_int
from pyNastran.dev.bdf_vectorized3.cards.constraints import ADD
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class NSMi(VectorizedBaseCard):
    """
    Defines a set of non structural mass.

    +-----+-----+------+----+-------+----+-------+----+-------+
    |  1  |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +=====+=====+======+====+=======+====+=======+====+=======+
    | NSM | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +-----+-----+------+----+-------+----+-------+----+-------+
    """
    _id_name = 'nsm_id'
    def clear(self) -> None:
        self.n = 0
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
        pid_eids = [pid_eid]
        values = [value]
        #self.cards.append((sid, nsm_type, pid_eid, value, comment))
        self.n += 1
        ifield = 5
        while card.field(ifield):
            pid_eid = integer(card, ifield, 'pid/eid')
            value = double(card, ifield+1, 'value')
            pid_eids.append(pid_eid)
            values.append(value)
            ifield += 2

        self.cards.append((sid, nsm_type, pid_eids, values, comment))
        #return cls(sid, nsm_type, pid_eid, value, comment=comment)
        return self.n - 1

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
        #if isinstance(pid_eid, integer_types):
        self.cards.append((sid, nsm_type, pid_eid, value, comment))
        self.n += 1
        #else:
        ##for pidi, valuei in zip(pid_eid, value):
        #self.cards.append((sid, nsm_type, pidi, valuei, comment))
        #comment = ''
        #self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        nsm_id = np.zeros(ncards, dtype='int32')
        nsm_type = np.zeros(ncards, dtype='|U7') # ELEMENT
        nvalues = np.zeros(ncards, dtype='int32')
        pid_eids = []
        values = []
        #pid_eid = np.zeros(ncards, dtype='int32')
        #value = np.zeros(ncards, dtype='float64')
        #I11, I21, I22, I31, I32, I33 = I
        #inertia = np.zeros((ncards, 6), dtype='float64')
        for icard, card in enumerate(self.cards):
            (sid, nsm_typei, pid_eidi, valuei, comment) = card
            nsm_id[icard] = sid
            nsm_type[icard] = nsm_typei
            assert len(nsm_typei) <= 7, f'nsm_type={nsm_typei!r} len={len(nsm_typei)}'
            #model.log.debug(pid_eidi)
            if isinstance(pid_eidi, int):
                pid_eidi = [pid_eidi]
                valuei = [valuei]
            pid_eids.extend(pid_eidi)
            values.extend(valuei)
            nvalues[icard] = len(valuei)
            #pid_eid[icard] = pid_eidi
            #value[icard] = valuei
        pid_eid = np.array(pid_eids, dtype='int32')
        value = np.array(values, dtype='float64')
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
            'PROD', 'PTUBE',
            'PBAR', 'PBARL',
            'PBEAM', 'PBEAML',
            'PSHELL', 'PSHEAR',
        }
        elements = {
            #'CQUAD4', 'CQUAD8', 'CQUADR',
            #'CTRIA3', 'CTRIA6', 'CTRIAR',
            'ELEMENT', 'CONROD',
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
            used_dict['property_id'] = np.unique(pids_used)
        if len(eids_used):
            used_dict['element_id'] = np.unique(eids_used)

    @property
    def ivalue(self) -> np.ndarray:
        return make_idim(self.n, self.nvalue)

    def __apply_slice__(self, nsm: NSMi, i: np.ndarray) -> None:
        nsm.nsm_id = self.nsm_id[i]
        nsm.nsm_type = self.nsm_type[i]
        #nsm.pid_eid = self.pid_eid[i]
        #nsm.value = self.value[i]
        nsm.nvalue = self.nvalue[i]
        idim = self.ivalue
        nsm.pid_eid = hslice_by_idim(i, idim, self.pid_eid)
        nsm.value = hslice_by_idim(i, idim, self.value)
        nsm.n = len(i)

    @property
    def ivalue(self) -> np.ndarray:
        return make_idim(self.n, self.nvalue)

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
    _id_name = 'nsm_id'
    def clear(self) -> None:
        self.n = 0
        self.nsm_id = np.array([], dtype='int32')
        self.nsm_type = np.array([], dtype='|U7')
        self.pid_eid = np.array([], dtype='int32')
        self.npid_eid = np.array([], dtype='int32')
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
        return self.n - 1

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
        return self.n - 1

    def __apply_slice__(self, elem: NSMi, i: np.ndarray) -> None:
        elem.nsm_id = self.nsm_id[i]
        elem.nsm_type = self.nsm_type[i]
        elem.value = self.value[i]

        ielement = self.ielement # [i, :]
        elem.pid_eid = hslice_by_idim(i, ielement, self.pid_eid)
        elem.npid_eid = self.npid_eid[i]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        properties = {
            'PROD', 'PTUBE',
            'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
            'PSHELL', 'PSHEAR',
        }
        elements = {
            #'CQUAD4', 'CQUAD8', 'CQUADR',
            #'CTRIA3', 'CTRIA6', 'CTRIAR',
            'ELEMENT', 'CONROD',
        }
        pids_used = []
        eids_used = []
        insm = self.ielement
        for nsm_type, (insm0, insm1) in zip_longest(
                self.nsm_type, insm):
            ids = self.pid_eid[insm0:insm1].tolist()
            if nsm_type in properties:
                pids_used.extend(ids)
            elif nsm_type in elements:
                eids_used.extend(ids)
            else:
                raise NotImplementedError(nsm_type)
        if len(pids_used):
            used_dict['property_id'] = np.unique(pids_used)
        if len(eids_used):
            used_dict['element_id'] = np.unique(eids_used)

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.npid_eid)

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        nsm_id = np.zeros(ncards, dtype='int32')

        # PSHELL, ELEMENT
        nsm_type = np.zeros(ncards, dtype='|U7')
        pid_eid_list = []
        npid_eid = np.zeros(ncards, dtype='int32')
        value = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (sid, nsm_typei, pid_eidi, valuei, comment) = card
            nsm_id[icard] = sid
            nsm_type[icard] = nsm_typei
            assert len(nsm_typei) <= 7, f'nsm_type={nsm_typei!r} len={len(nsm_typei)}'
            if isinstance(pid_eidi, int):
                pid_eidi = [pid_eidi]
            elif isinstance(pid_eidi, str):  #  pragma: no cover
                if pid_eidi == 'ALL':
                    pid_eidi = [-1]
                else:
                    raise RuntimeError(pid_eidi)
            elif isinstance(pid_eidi, list):
                pid_eidi = expand_thru_by(pid_eidi, set_fields=True, sort_fields=True, require_int=True, allow_blanks=False)
                assert 'THRU' not in pid_eidi
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
        if len(self.nsm_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        nsm_str = array_str(self.nsm_id, size=size)
        pid_eid_str = array_str(self.pid_eid, size=size)
        #insm = self.insm
        ielement = self.ielement
        for nsm_id, nsm_type, value, (insm0, insm1) in zip_longest(
                nsm_str, self.nsm_type, self.value, ielement):
            ids = self.pid_eid[insm0:insm1].tolist()
            list_fields = ['NSM1', nsm_id, nsm_type, value, ] + ids
            bdf_file.write(print_card(list_fields))
        return

class NSML1(NSM1i):
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.nsm_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        nsm_str = array_str(self.nsm_id, size=size)
        #pid_eid_str = array_str(self.pid_eid, size=size)
        insm = self.ielement
        for nsm_id, nsm_type, value, (insm0, insm1) in zip_longest(
                nsm_str, self.nsm_type, self.value, insm):
            ids = self.pid_eid[insm0:insm1].tolist()
            list_fields = ['NSML1', nsm_id, nsm_type, value, ] + ids
            bdf_file.write(print_card(list_fields))
        return

class NSM(NSMi):
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.nsm_id) == 0:
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
        if len(self.nsm_id) == 0:
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
    _id_name = 'nsm_id'
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
        used_dict['nsm_id'] = self.nsm_ids

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_headers: bool=False) -> None:
        if len(self.nsm_id) == 0:
            return
        print_card = get_print_card_8_16(size)
        nsm_id = array_str(self.nsm_id, size=size)
        nsm_ids = array_str(self.nsm_ids, size=size)
        for nsm_id, idim in zip(nsm_id, self.idim):
            idim0, idim1 = idim
            assert idim1 > idim0, self.idim
            nsm_idsi = nsm_ids[idim0:idim1].tolist()
            assert len(nsm_idsi) > 0, self.idim
            list_fields = ['NSMADD', nsm_id] + nsm_idsi
            bdf_file.write(print_card(list_fields))
        return

    def get_nsms_by_nsm_id(self) -> dict[int, NSMs]:
        model = self.model
        """"""
        #unsm_ids = np.unique(self.nsm_id)
        nsm_by_nsm_id = defaultdict(list)
        #for nsm_id in unsm_ids:
            #nsm_by_nsm_id[nsm_id] = []

        for nsm in model.nsm_cards:
            if nsm.type in {'NSMADD'}:
                continue

            unsm_idsi = np.unique(nsm.nsm_id)
            for unsm_id in unsm_idsi:
                if unsm_id not in nsm.nsm_id:
                    continue

                i = np.where(unsm_id == nsm.nsm_id)[0]
                if len(i) == 0:
                    continue
                #print(nsm.type, unsm_id, nsm.nsm_id, i)
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

    #def slice_card_by_nsm_id(self, nsm_id: np.ndarray) -> NSMADD:
        #assert self.n > 0, self.n
        #assert len(self.element_id) > 0, self.element_id
        #i = self.index(element_id)
        ##cls_obj = cls(self.model)
        ##cls_obj.__apply_slice__(self, i)
        #cls_obj = self.slice_card_by_index(i)
        #assert cls_obj.n > 0, cls_obj
        #return cls_obj

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
        #nsm.components = self.components[i]

        #nsm.nnodes = self.nnodes[i]
        idim = self.idim
        nsm.sids = hslice_by_idim(i, idim, self.sids)
        return nsm

    #def slice_card_by_id(self, spc_id: int) -> SPC1:
        #i = np.where(spc_id == self.spc_id)[0]
        #spc = SPC1(self.model)
        #self.__apply_slice__(spc, i)
        #return spc

    #def slice_card_by_index(self, i: np.ndarray) -> SPC1:
        #spc = SPC1(self.model)
        #self.__apply_slice__(spc, i)
        #return spc

    #def __apply_slice__(self, spc: NSMADD, i: np.ndarray) -> None:
        #spc.n = len(i)
        #spc.spc_id = self.spc_id[i]
        #spc.components = self.components[i]

        #spc.nnodes = self.nnodes[i]
        #idim = self.inode
        #spc.node_id = hslice_by_idim(i, idim, self.node_id)
        #return spc


