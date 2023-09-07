from __future__ import annotations
from typing import Union, TYPE_CHECKING
from collections import defaultdict
from itertools import count
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard, hslice_by_idim, make_idim
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    components_or_blank)
#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16 # print_float_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.dev.bdf_vectorized3.cards.base_card import get_print_card_8_16
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_default_str, array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.utils import cast_int_array, print_card_8

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class SPC(VectorizedBaseCard):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).

     +-----+-----+----+----+------+----+----+----+
     |  1  |  2  |  3 |  4 |   5  |  6 |  7 |  8 |
     +=====+=====+====+====+======+====+====+====+
     | SPC | SID | G1 | C1 |  D1  | G2 | C2 | D2 |
     +-----+-----+----+----+------+----+----+----+
     | SPC |  2  | 32 | 3  | -2.6 |  5 |    |    |
     +-----+-----+----+----+------+----+----+----+
    """
    _id_name = 'spc_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.spc_id = np.array([])

    def slice_card_by_index(self, i: np.ndarray) -> SPC:
        spc = SPC(self.model)
        self.__apply_slice__(spc, i)
        return spc

    def __apply_slice__(self, spc: SPC, i: np.ndarray) -> None:
        spc.n = len(i)
        spc.spc_id = self.spc_id[i]
        spc.node_id = self.node_id[i]
        spc.components = self.components[i]
        spc.enforced = self.enforced[i]
        return spc

    def add(self, spc_id: int, nodes: list[int], components: list[str],
            enforced: list[float], comment: str=''):
        """
        Creates an SPC card, which defines the degree of freedoms to be
        constrained

        Parameters
        ----------
        conid : int
            constraint id
        nodes : list[int]
            GRID/SPOINT ids
        components : list[int]
            the degree of freedoms to constrain (e.g., 1, 123)
        enforced : list[float]
            the constrained value for the given node (typically 0.0)
        comment : str; default=''
            a comment for the card

        Notes
        -----
        len(nodes) == len(components) == len(enforced)

        .. warning:: non-zero enforced deflection requires an SPCD as well

        """
        if isinstance(spc_id, integer_types):
            spc_id = [spc_id]
        self.cards.append((spc_id, nodes, components, enforced, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        spc_idi = integer(card, 1, 'sid')
        if card.field(5) in [None, '']:
            spc_id = [spc_idi]
            nodes = [integer(card, 2, 'G1'),]
            components = [components_or_blank(card, 3, 'C1', default=0)]
            enforced = [double_or_blank(card, 4, 'D1', default=0.0)]
        else:
            spc_id = [spc_idi, spc_idi]
            nodes = [
                integer(card, 2, 'G1'),
                integer_or_blank(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            components = [
                components_or_blank(card, 3, 'C1', default=0),
                components_or_blank(card, 6, 'C2', default=0),
            ]
            enforced = [
                double_or_blank(card, 4, 'D1', default=0.0),
                double_or_blank(card, 7, 'D2', default=0.0),
            ]

        #nnodes = len(nodesi)
        #inid1 = inid0 + nnodes
        self.cards.append((spc_id, nodes, components, enforced, comment))
        self.n += 1

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        spc_id_ = []
        components_ = []
        node_id_ = []
        enforced_ = []
        for icard, card in enumerate(self.cards):
            (spc_idi, nodesi, componentsi, enforcedi, comment) = card
            spc_id_.extend(spc_idi)
            node_id_.extend(nodesi)
            enforced_.extend(enforcedi)
            components_.extend(componentsi)
        spc_id = np.array(spc_id_, dtype='int32')
        node_id = np.array(node_id_, dtype='int32')
        components = np.array(components_, dtype='int32')
        enforced = np.array(enforced_, dtype='float64')
        self._save(spc_id, node_id, components, enforced)
        assert len(self.spc_id) == len(self.node_id)
        assert len(self.spc_id) == len(self.components)
        assert len(self.spc_id) == len(self.enforced)
        self.cards = []

    def _save(self, spc_id, node_id, components, enforced):
        nspcs = len(spc_id)
        self.spc_id = spc_id
        self.node_id = node_id
        self.components = components
        self.enforced = enforced
        self.n = nspcs

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spc_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        spc_str = array_str(self.spc_id, size=size)
        node_str = array_str(self.node_id, size=size)
        components_str = array_default_int(self.components, default=0, size=size)

        no_enforced = self.enforced.max() == 0. and self.enforced.min() == 0.
        if no_enforced:
            for spc_id, node_id, components in zip(spc_str, node_str, components_str):
                msg = 'SPC     %8s%8s%8s\n' % (spc_id, node_id, components)
                bdf_file.write(msg)
        else:
            for spc_id, node_id, components, enforced in zip(spc_str, node_str, components_str, self.enforced):
                list_fields = ['SPC', spc_id, node_id, components, enforced]
                bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id))


class SPC1(VectorizedBaseCard):
    """
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    |  1   |  2  |   3   |    4   |   5    |    6   |    7   |    8   |  9 |
    +======+=====+=======+========+========+========+========+========+====+
    | SPC1 | SID |   C   |   G1   |   G2   |   G3   |   G4   |   G5   | G6 |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    |      |  G7 |   G8  |   G9   |  etc.  |        |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 |  3  |  246  | 209075 | 209096 | 209512 | 209513 | 209516 |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 |  3  |   2   |   1    |   3    |   10   |   9    |   6    | 5  |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    |      |  2  |   8   |        |        |        |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 | SID |   C   |   G1   |  THRU  |   G2   |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 | 313 | 12456 |    6   |  THRU  |   32   |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+

    """
    _id_name = 'spc_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.spc_id = np.array([])

    def slice_card_by_id(self, spc_id: int) -> SPC1:
        i = np.where(spc_id == self.spc_id)[0]
        spc = SPC1(self.model)
        self.__apply_slice__(spc, i)
        return spc


    def slice_card_by_index(self, i: np.ndarray) -> SPC1:
        spc = SPC1(self.model)
        self.__apply_slice__(spc, i)
        return spc

    def __apply_slice__(self, spc: SPC1, i: np.ndarray) -> None:
        spc.n = len(i)
        spc.spc_id = self.spc_id[i]
        spc.components = self.components[i]

        spc.nnodes = self.nnodes[i]
        idim = self.inode
        spc.node_id = hslice_by_idim(i, idim, self.node_id)
        return spc

    def add(self, spc_id: int, components: int, nodes: list[int], comment: str='') -> None:
        """
        Creates an SPC1 card, which defines the degree of freedoms to be
        constrained to a value of 0.0

        Parameters
        ----------
        conid : int
            constraint id
        components : int
            the degree of freedoms to constrain (e.g., 1, 123)
        nodes : list[int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        self.cards.append((spc_id, components, nodes, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str=''):
        spc_id = integer(card, 1, 'conid')
        components = components_or_blank(card, 2, 'components', default=0)  # 246 = y; dx, dz dir
        #nodes = [node for node in card.fields(3) if node is not None]
        node_fields = card.fields(3)
        nnodes = len(node_fields)
        assert nnodes > 0, node_fields
        nodes = expand_thru(node_fields, set_fields=True, sort_fields=True)
        assert 'THRU' not in nodes, nodes
        self.cards.append((spc_id, components, nodes, comment))
        self.n += 1

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        spc_id = np.zeros(ncards, dtype='int32')
        components = np.zeros(ncards, dtype='int32')
        nnodes = np.zeros(ncards, dtype='int32')
        node_id = []
        for icard, card in enumerate(self.cards):
            (spc_idi, componentsi, nodesi, comment) = card
            nnodesi = len(nodesi)
            spc_id[icard] = spc_idi
            nnodes[icard] = nnodesi
            node_id.extend(nodesi)
            components[icard] = componentsi
        node_id_ = cast_int_array(node_id)
        self._save(spc_id, node_id_, components, nnodes)
        self.cards = []

    def _save(self, spc_id, node_id, components, nnodes):
        nspcs = len(spc_id)
        self.spc_id = spc_id
        self.node_id = node_id
        self.components = components
        self.nnodes = nnodes
        self.n = nspcs

    @property
    def inode(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spc_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        spc_str = array_str(self.spc_id, size=size)
        #node_str = array_str(self.node_id, size=size)
        components_str = array_default_int(self.components, default=0, size=size)

        for spc_id, components, inode in zip(spc_str, components_str, self.inode):
            inode0, inode1 = inode
            nodes = self.node_id[inode0 : inode1].tolist()
            assert len(nodes) > 0, nodes
            list_fields = ['SPC1', spc_id, components] + nodes
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),)

class MPC(VectorizedBaseCard):
    """
    Multipoint Constraint
    Defines a multipoint constraint equation of the form:
      sum(A_j * u_j) = 0

    where:
      uj represents degree-of-freedom Cj at grid or scalar point Gj.
      Aj represents the scale factor

    +-----+-----+----+----+-----+----+----+----+-----+
    |  1  |  2  |  3 |  4 |  5  |  6 |  7 |  8 |  9  |
    +=====+=====+====+====+=====+====+====+====+=====+
    | MPC | SID | G1 | C1 |  A1 | G2 | C2 | A2 |     |
    +-----+-----+----+----+-----+----+----+----+-----+
    |     |  G3 | C3 | A3 | ... |    |    |    |     |
    +-----+-----+----+----+-----+----+----+----+-----+
    """
    _id_name = 'mpc_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.mpc_id = np.array([])

    def add(self, mpc_id: int,
            nodes: list[int], components: list[str], coefficients: list[float],
            comment: str=''):
        self.cards.append((mpc_id, nodes, components, coefficients, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str=''):
        mpc_id = integer(card, 1, 'conid')
        nodes = []
        components = []
        coefficients = []

        nfields = len(card)

        i = 1
        for ifield in range(2, nfields, 8):
            nid = integer(card, ifield, 'G%d' % i)
            component = components_or_blank(card, ifield + 1, 'constraint%d' % i, '0')  # scalar point
            if i == 1:
                coefficient = double(card, ifield + 2, 'coefficient%d' % i)
                if coefficient == 0.0:
                    raise RuntimeError('coefficient1 must be nonzero; coefficient=%r' % coefficient)
            else:
                coefficient = double_or_blank(card, ifield + 2, 'coefficient%d' % i, 0.0)
            nodes.append(nid)
            components.append(component)
            coefficients.append(coefficient)
            i += 1

            if ifield + 4 > nfields and i != 2:
                # if G2 is empty (it's ifield+4 because nfields is length based
                # and not loop friendly)
                break
            nid = integer(card, ifield + 3, 'G%d' % i)
            component = int(components_or_blank(card, ifield + 4, 'constraint%d' % i, default=0))  # scalar point
            coefficient = double_or_blank(card, ifield + 5, 'coefficient%d' % i, default=0.0)
            nodes.append(nid)
            components.append(component)
            coefficients.append(coefficient)
            i += 1
        self.cards.append((mpc_id, nodes, components, coefficients, comment))
        self.n += 1

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.mpc_id = np.zeros(ncards, dtype='int32')
        self.idim = np.zeros((ncards, 2), dtype='int32')
        all_nodes = []
        all_components = []
        all_coefficients = []
        idim0 = 0
        for icard, card in enumerate(self.cards):
            (mpc_id, nodes, components, coefficients, comment) = card

            nnodes = len(nodes)
            idim1 = idim0 + nnodes
            self.mpc_id[icard] = mpc_id
            self.idim[icard, :] = [idim0, idim1]
            all_nodes.extend(nodes)
            all_components.extend(components)
            all_coefficients.extend(coefficients)
            idim0 = idim1
        self.node_id = np.array(all_nodes, dtype='int32')
        self.components = np.array(all_components, dtype='int32')
        self.coefficients = np.array(all_coefficients, dtype='float64')
        #assert len(self.mpc_id) == len(self.node_id)
        #assert len(self.mpc_id) == len(self.components)
        #assert len(self.mpc_id) == len(self.coefficients)
        self.cards = []

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.mpc_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        mpc_str = array_str(self.mpc_id, size=size)
        node_str = array_str(self.node_id, size=size)
        components_str = array_default_int(self.components, default=0, size=size)

        for mpc_id, idim in zip(mpc_str, self.idim):
            idim0, idim1 = idim
            nodes = node_str[idim0:idim1]
            components = components_str[idim0:idim1]
            coeffs = self.coefficients[idim0:idim1]

            list_fields = ['MPC', mpc_id]
            for i, gid, component, coefficient in zip(count(), nodes, components, coeffs):
                list_fields += [gid, component, coefficient]
                if i % 2 == 1 and i > 0:
                    list_fields.append(None)
                    list_fields.append(None)
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        all_spoint_ids = self.model.spoint.spoint_id

        ispoint = np.where(self.components == 0)
        spoint_id = self.node_id[ispoint]

        inode = where_not(self.components, ispoint)
        geom_check(self,
                   missing,
                   node=(nid, self.node_id[inode]),
                   spoint=(all_spoint_ids, spoint_id),
                   )

def where_not(base_vector: np.ndarray, iwhere: np.ndarray) -> np.ndarray:
    mask = np.full(base_vector.shape, True, dtype='bool')
    mask[iwhere] = False
    return mask


class ADD(VectorizedBaseCard):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.sid = np.array([], dtype='int32')
        self.sids = np.array([], dtype='int32')
        self.nsids = np.array([], dtype='int32')

    def add(self, sid, sets, comment='') -> ADD:
        """Creates an MPCADD card"""
        if isinstance(sets, integer_types):
            sets = [sets]
        self.cards.append((sid, sets, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        sid = integer(card, 1, 'conid')
        sets = np.unique(card.fields(2)).tolist()
        nset_cards = len(sets)
        assert nset_cards > 0, f'nset_cards={nset_cards} card={card}'

        self.cards.append((sid, sets, comment))
        self.n += 1

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        sid = np.zeros(ncards, dtype='int32')
        nsids = np.zeros(ncards, dtype='int32')

        sids_list = []
        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sidi, setsi, comment) = card
            sid[icard] = sidi

            nset_cards = len(setsi)
            assert nset_cards > 0, f'nset_cards={nset_cards} card={card}'
            sids_list.extend(setsi)
            sid[icard] = sidi
            nsids[icard] = nset_cards
        sids = np.array(sids_list, dtype='int32')
        self._save(sid, sids, nsids)
        #self.sort()
        self.cards = []

    def _save(self, sid, sids, nsids):
        if len(self.sid):
            sid = np.hstack([self.sid, sids])
            sids = np.hstack([self.sids, sids])
            nsids = np.hstack([self.nsids, nsids])
        nsid = len(sid)
        self.sid = sid
        self.sids = sids
        self.nsids = nsids
        self.n = nsid

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nsids)


class SPCADD(ADD):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | SPCADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    _id_name = 'spc_id'
    #def __init__(self, model: BDF):
        #super().__init__(model)
        #self.sid = np.array([], dtype='int32')
        #self.spc_ids = np.array([], dtype='int32')
    @property
    def spc_id(self):
        return self.sid
    @spc_id.setter
    def spc_id(self, spc_id: np.ndarray):
        self.sid = spc_id

    @property
    def spc_ids(self):
        return self.sids
    @spc_ids.setter
    def spc_ids(self, spc_ids: np.ndarray):
        self.sids = spc_ids

    @property
    def nspcs(self):
        return self.nsids
    @nspcs.setter
    def nspcs(self):
        return self.nsids

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool = False,
                   write_card_header: bool=False) -> None:
        if len(self.spc_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        spc_ids = array_str(self.spc_ids, size=size)
        for spc_id, idim in zip(self.spc_id, self.idim):
            idim0, idim1 = idim
            assert idim1 > idim0, self.idim
            spc_idsi = spc_ids[idim0:idim1].tolist()
            assert len(spc_idsi) > 0, self.idim
            list_fields = ['SPCADD', spc_id] + spc_idsi
            bdf_file.write(print_card(list_fields))
        return

    def get_spcs_by_spc_id(self) -> dict[int, SPCs]:
        model = self.model
        """"""
        #uspc_ids = np.unique(self.spc_id)
        spc_by_spc_id = defaultdict(list)
        #for spc_id in uspc_ids:
            #spc_by_spc_id[spc_id] = []

        for spc in model.spcs:
            if spc.type in {'SPCADD'}:
                continue

            uspc_idsi = np.unique(spc.spc_id)
            for uspc_id in uspc_idsi:
                i = np.where(uspc_id == spc.spc_id)[0]
                if len(i) == 0:
                    continue
                spci = spc.slice_card_by_index(i)
                spc_by_spc_id[uspc_id].append(spci)
        return dict(spc_by_spc_id)

    def get_reduced_spcs(self,
                         #resolve_load_card: bool=False,
                         stop_on_failure: bool=True) -> dict[int, SPCs]:
        """
        Parameters
        ----------
        resolve_load_card : bool; default=False
            ???
        """
        stop_on_failure = True
        spc_by_spc_id = self.get_spcs_by_spc_id()
        log = self.model.log

        reduced_spcs = {}
        for sid, idim in zip(self.spc_id, self.idim):
            reduced_spcsi = []
            idim0, idim1 = idim

            spc_ids = self.spc_ids[idim0:idim1]
            for spc_idi in spc_ids:
                spcs_found = spc_by_spc_id[spc_idi]
                if len(spcs_found) == 0:
                    msg = f'No referenced SPCs found for spc_id={spc_idi} on SPCADD spc_id={sid}'
                    log.error(msg)
                    if stop_on_failure:
                        raise RuntimeError(msg)
                reduced_spcsi.append(spcs_found)
            reduced_spcs[sid] = reduced_spcsi
        return reduced_spcs


class MPCADD(ADD):
    _id_name = 'mpc_id'
    @property
    def mpc_id(self):
        return self.sid
    @mpc_id.setter
    def mpc_id(self, spc_id: np.ndarray):
        self.sid = spc_id

    @property
    def mpc_ids(self):
        return self.sids
    @mpc_ids.setter
    def mpc_ids(self, spc_ids: np.ndarray):
        self.sids = spc_ids

    @property
    def nmpcs(self):
        return self.nsids
    @nmpcs.setter
    def nmpcs(self):
        return self.nsids

    @property
    def is_small_field(self):
        return max(self.mpc_id.max(), self.mpc_ids.max()) < 99_999_999

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> str:
        if len(self.mpc_id) == 0:
            return
        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16
        #self.get_reduced_spcs()
        mpc_ids = array_str(self.mpc_ids, size=size)
        for mpc_id, idim in zip(self.mpc_id, self.idim):
            idim0, idim1 = idim
            mpc_idsi = mpc_ids[idim0:idim1].tolist()
            list_fields = ['MPCADD', mpc_id] + mpc_idsi
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        mpc_id = np.unique(self.model.mpc.mpc_id)
        geom_check(self,
                   missing,
                   mpc=(mpc_id, self.mpc_ids))

SPCs = Union[SPC, SPC1]

