from __future__ import annotations
from typing import TYPE_CHECKING
from collections import defaultdict
from itertools import count
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.cards.collpase_card import collapse_thru
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard, hslice_by_idim, make_idim
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    components_or_blank, parse_components,)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    remove_unused_primary, remove_unused_duplicate, parse_check)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    #array_default_str,
    array_str, array_default_int,
    get_print_card_size)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF


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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spc_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self.components = np.array([], dtype='int32')
        self.enforced = np.array([], dtype='float64')

    #def slice_card_by_index(self, i: np.ndarray, **kwargs) -> SPC:
        #spc = SPC(self.model)
        #self.__apply_slice__(spc, i)
        #return spc

    def __apply_slice__(self, spc: SPC, i: np.ndarray) -> None:
        spc.n = len(i)
        spc.spc_id = self.spc_id[i]
        spc.node_id = self.node_id[i]
        spc.components = self.components[i]
        spc.enforced = self.enforced[i]
        return spc

    def add(self, spc_id: int, nodes: list[int], components: list[str],
            enforced: list[float],
            ifile: int=0, comment: str='') -> int:
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
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        nnodes = len(nodes)

        if isinstance(components, integer_types):
            components = [components]
        if isinstance(enforced, float_types):
            enforced = [enforced] * nnodes

        if isinstance(spc_id, integer_types):
            spc_id = [spc_id] * nnodes
        assert nnodes == len(components)
        assert nnodes == len(enforced)
        self.cards.append((spc_id, nodes, components, enforced, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        self.cards.append((spc_id, nodes, components, enforced, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        idtype = self.model.idtype
        spc_id_ = []
        components_ = []
        node_id_ = []
        enforced_ = []
        for icard, card in enumerate(self.cards):
            (spc_idi, nodesi, componentsi, enforcedi, ifilei, comment) = card
            spc_id_.extend(spc_idi)
            node_id_.extend(nodesi)
            enforced_.extend(enforcedi)
            components_.extend(componentsi)

        spc_id = np.array(spc_id_, dtype=idtype)
        node_id = np.array(node_id_, dtype=idtype)
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

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.node_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        #node_ids = used_dict['node_id']
        spc_ids = used_dict['spc_id']

        #nids_to_remove = np.intersect1d
        ncards_removed = remove_unused_duplicate(self, spc_ids, self.spc_id, 'spc_id')
        return ncards_removed

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        if np.abs(self.enforced).max() == 0.:
            return
        nids_list = []
        components_list = []
        enforced_list = []
        rotation_dof_set = set('456')
        translation_dof_set = set('123')
        for nidi, compi, enforcedi in zip(self.node_id, self.components, self.enforced):
            if enforcedi == 0.0:
                nids_list.append(nidi)
                components_list.append(compi)
                enforced_list.append(enforcedi)
                continue

            comp_str = str(compi)
            comp_str = ''.join(list(sorted(comp_str)))  # sort the string
            if comp_str in {'1', '2', '3', '12', '13', '23', '123'}:
                # pure translation - displacement
                nids_list.append(nidi)
                components_list.append(compi)
                enforced_list.append(enforcedi * xyz_scale)
            elif comp_str in {'4', '5', '6', '45', '46', '56', '456'}:
                # pure rotation - radians?
                nids_list.append(nidi)
                components_list.append(compi)
                enforced_list.append(enforcedi)
            else:
                dofs = set(comp_str)
                translation_dofs = dofs.difference(rotation_dof_set)
                rotation_dofs = dofs.difference(translation_dof_set)
                if len(translation_dofs):
                    components_list.append(nidi)
                    components_list.append(translation_dofs)
                    enforced_list.append(enforcedi)
                if len(rotation_dofs):
                    components_list.append(nidi)
                    components_list.append(rotation_dofs)
                    enforced_list.append(enforcedi)
        spc_id = np.array(self.spc_id, dtype=self.spc_id.dtype)
        node_id = np.array(self.node_id, dtype=self.node_id.dtype)
        components = np.array(self.components, dtype=self.components.dtype)
        enforced = np.array(self.enforced, dtype=self.enforced.dtype)
        self._save(spc_id, node_id, components, enforced)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),
                   )
    @property
    def max_id(self) -> int:
        return max(self.spc_id.max(), self.node_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        spc_str = array_str(self.spc_id, size=size)
        node_str = array_str(self.node_id, size=size)
        components_str = array_default_int(self.components, default=0, size=size)

        no_enforced = self.enforced.max() == 0. and self.enforced.min() == 0.
        if no_enforced:
            for spc_id, node_id, components in zip(spc_str, node_str, components_str):
                list_fields = ['SPC', spc_id, node_id, components]
                #msg = 'SPC     %8s%8s%8s\n' % (spc_id, node_id, components)
                #bdf_file.write(msg)
                bdf_file.write(print_card(list_fields))
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spc_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self.components = np.array([], dtype='int32')
        self.nnodes = np.array([], dtype='int32')
        self.enforced = np.array([], dtype='float64')

    def slice_card_by_id(self, spc_id: int, **kwargs) -> SPC1:
        i = np.where(spc_id == self.spc_id)[0]
        spc = SPC1(self.model)
        self.__apply_slice__(spc, i)
        return spc

    #def slice_card_by_index(self, i: np.ndarray) -> SPC1:
        #spc = SPC1(self.model)
        #self.__apply_slice__(spc, i)
        #return spc

    def __apply_slice__(self, spc: SPC1, i: np.ndarray) -> None:
        spc.n = len(i)
        spc.spc_id = self.spc_id[i]
        spc.components = self.components[i]

        spc.nnodes = self.nnodes[i]
        idim = self.inode
        spc.node_id = hslice_by_idim(i, idim, self.node_id)
        return spc

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.node_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        #node_ids = used_dict['node_id']
        spc_ids = used_dict['spc_id']
        ncards_removed = remove_unused_duplicate(self, spc_ids, self.spc_id, 'spc_id')
        return ncards_removed

    def add(self, spc_id: int, components: int, nodes: list[int],
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((spc_id, components, nodes, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        spc_id = integer(card, 1, 'conid')
        components = components_or_blank(card, 2, 'components', default=0)  # 246 = y; dx, dz dir
        #nodes = [node for node in card.fields(3) if node is not None]
        node_fields = card.fields(3)
        nnodes = len(node_fields)
        assert nnodes > 0, node_fields
        nodes = expand_thru(node_fields, set_fields=True, sort_fields=True)
        assert 'THRU' not in nodes, nodes
        self.cards.append((spc_id, components, nodes, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        spc_id = np.zeros(ncards, dtype='int32')
        components = np.zeros(ncards, dtype='int32')
        nnodes = np.zeros(ncards, dtype='int32')
        node_id_list = []
        for icard, card in enumerate(self.cards):
            (spc_idi, componentsi, nodesi, ifilei, comment) = card
            nnodesi = len(nodesi)
            spc_id[icard] = spc_idi
            nnodes[icard] = nnodesi
            node_id_list.extend(nodesi)
            components[icard] = componentsi
        node_id = np.array(node_id_list, dtype=idtype)
        self._save(spc_id, node_id, components, nnodes)
        self.cards = []

    def _save(self, spc_id, node_id, components, nnodes):
        nspcs = len(spc_id)
        self.spc_id = spc_id
        self.node_id = node_id
        self.components = components
        self.nnodes = nnodes
        self.n = nspcs

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def inode(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    @property
    def max_id(self) -> int:
        return max(self.spc_id.max(), self.node_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        spc_str = array_str(self.spc_id, size=size)
        node_str = array_str(self.node_id, size=size)
        components_str = array_default_int(self.components, default=0, size=size)

        for spc_id, components, inode in zip(spc_str, components_str, self.inode):
            inode0, inode1 = inode
            nodes = node_str[inode0 : inode1].tolist()
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.mpc_id = np.array([], dtype='int32')

    def add(self, mpc_id: int,
            nodes: list[int], components: list[str], coefficients: list[float],
            ifile: int=0, comment: str=''):
        self.cards.append((mpc_id, nodes, components, coefficients, ifile, comment))
        self.n += 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        self.cards.append((mpc_id, nodes, components, coefficients, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        all_nodes = []
        all_components = []
        all_coefficients = []
        idim0 = 0

        mpc_id = np.zeros(ncards, dtype='int32')
        nnode = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (mpc_idi, nodes, components, coefficients, ifilei, comment) = card

            nnodes = len(nodes)
            idim1 = idim0 + nnodes
            mpc_id[icard] = mpc_idi
            #idim[icard, :] = [idim0, idim1]
            nnode[icard] = len(nodes)
            all_nodes.extend(nodes)
            all_components.extend(components)
            all_coefficients.extend(coefficients)
            idim0 = idim1

        node_id = np.array(all_nodes, dtype=idtype)
        components = np.array(all_components, dtype='int32')
        coefficients = np.array(all_coefficients, dtype='float64')
        self._save(mpc_id, nnode, node_id, components, coefficients)
        #assert len(self.mpc_id) == len(self.node_id)
        #assert len(self.mpc_id) == len(self.components)
        #assert len(self.mpc_id) == len(self.coefficients)
        self.cards = []

    def _save(self, mpc_id, nnode, node_id, components, coefficients) -> None:
        if len(self.mpc_id) != 0:
            asdf
        self.mpc_id = mpc_id
        self.nnode = nnode
        self.node_id = node_id
        self.components = components
        self.coefficients = coefficients

    def __apply_slice__(self, mpc: MPC, i: np.ndarray):
        """this kind of doesn't make sense..."""
        mpc.mpc_id = self.mpc_id[i]
        impc = self.idim
        mpc.node_id = hslice_by_idim(i, impc, self.node_id)
        mpc.components = hslice_by_idim(i, impc, self.components)
        mpc.coefficients = hslice_by_idim(i, impc, self.coefficients)
        mpc.nnode = self.nnode[i]
        mpc.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnode)

    @property
    def max_id(self) -> int:
        return max(self.mpc_id.max(), self.node_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

def where_not(base_vector: np.ndarray, iwhere: np.ndarray) -> np.ndarray:
    mask = np.full(base_vector.shape, True, dtype='bool')
    mask[iwhere] = False
    return mask


class ADD(VectorizedBaseCard):
    def clear(self) -> None:
        self.sid = np.array([], dtype='int32')
        self.sids = np.array([], dtype='int32')
        self.nsids = np.array([], dtype='int32')
        self.n = 0

    @property
    def max_id(self) -> int:
        return max(self.sid.max(), self.sids.max())

    def add(self, sid: int, sets,
            ifile: int=0, comment: str='') -> int:
        """Creates an MPCADD card"""
        if isinstance(sets, integer_types):
            sets = [sets]
        self.cards.append((sid, sets, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        sid = integer(card, 1, 'conid')
        sets = np.unique(card.fields(2)).tolist()
        nset_cards = len(sets)
        assert nset_cards > 0, f'nset_cards={nset_cards} card={card}'

        self.cards.append((sid, sets, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        sid = np.zeros(ncards, dtype='int32')
        nsids = np.zeros(ncards, dtype='int32')

        sids_list = []
        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sidi, setsi, ifilei, comment) = card
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
        assert nsids.min() > 0, nsids
        self.n = nsid

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nsids)

    def __apply_slice__(self, obj: ADD, i: np.ndarray) -> None:
        obj.n = len(i)
        obj.sid = self.sid[i]

        obj.nsids = self.nsids[i]
        idim = self.idim
        obj.sids = hslice_by_idim(i, idim, self.sids)


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

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['spc_id'].append(self.spc_id)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool = False,
                   write_card_header: bool=False) -> None:
        if len(self.spc_id) == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

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

        for spc in model.spc_cards:
            if spc.type in {'SPCADD', 'SPCOFF', 'BNDFIX', 'BNDFREE'}:
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
                reduced_spcsi.extend(spcs_found)
            reduced_spcs[sid] = reduced_spcsi
        return reduced_spcs

    def get_reduced_node_component(self,
                                   reduced_spc_dict: dict[int, list[SPCs]],
                                   ) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        """Gets the node/component dict for an SPCADD"""
        compressed_spc_dict = {}
        for spcadd_id, cards2 in sorted(reduced_spc_dict.items()):
            all_spc_ids_list = []
            for card in cards2:
                #all_spc_ids_list.append(card.components)
                all_spc_ids_list.append(card.spc_id)
            uspc_ids = np.unique(np.hstack(all_spc_ids_list))

            nids_all_list = []
            comp_all_list = []
            for spc_idi in uspc_ids:
                is_failed, nidsi, compsi = spc_cards_to_nid_dof(spc_idi, cards2)
                if is_failed:
                    continue
                nids_all_list.append(nidsi)
                comp_all_list.append(compsi)
            nids = np.hstack(nids_all_list)
            comp = np.hstack(comp_all_list)
            compressed_spc_dict[spcadd_id] = (nids, comp)
        return compressed_spc_dict

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

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['mpc_id'].append(self.mpc_id)

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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


class CommonSet(VectorizedBaseCard):
    _id_name = 'node_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.node_id = np.array([], dtype='int32')
        self.component = np.array([], dtype='int32')

    def add_set(self, nids: list[int], components: list[int],
                comment: str='') -> int:
        assert isinstance(nids, (list, np.ndarray, tuple))
        assert isinstance(components, (list, np.ndarray, tuple))
        nnodes = len(nids)
        ncomp = len(components)
        assert nnodes == ncomp, (nnodes, ncomp)
        suport_id = 0
        self.cards.append((suport_id, nids, components, ifile, comment))
        #if comment:
            #self.comment[nid] = _format_comment(comment)
        self.n += nnodes
        return self.n - 1

    def add_set1(self, suport_id: int, nids: list[int], component: list[int],
                  comment: str='') -> int:
        assert isinstance(component, (str, integer_types)), component
        nids = expand_thru(nids, set_fields=True, sort_fields=False)
        nnodes = len(nids)
        components = [component] * nnodes
        self.cards.append((suport_id, nids, components, ifile, comment))
        #if comment:
            #self.comment[nid] = _format_comment(comment)
        self.n += nnodes
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str=''):
        card_name = card[0].upper()
        msg = f'add_card(...) has been removed for {card_name}.  Use add_set_card or add_set1_card'
        raise AttributeError(msg)

    def add_set_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a SPCOFF/BNDFIX/BNDFREE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        nodes = []
        components = []
        nfields = len(card) - 1
        nconstraints = nfields // 2
        if nfields % 2 == 1:
            nconstraints += 1
        for counter in range(nconstraints):
            igrid = counter + 1
            ifield = counter * 2 + 1
            node = integer(card, ifield, 'G%i' % igrid)
            component = components_or_blank(card, ifield+1, 'C%i' % igrid, default='0')
            nodes.append(node)
            components.append(component)
        assert len(card) > 1, f'len({self.type} card) = {len(card):d}\ncard={card}'
        self.cards.append((nodes, components, ifile, comment))
        self.n += len(nodes)
        return self.n - 1

    def add_set1_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a SPCOFF1/BNDFIX1/BNDFREE1/BNDGRID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        component = parse_components(card, 1, 'components')  # 246 = y; dx, dz dir
        nodes = card.fields(2)
        assert len(card) > 2, f'len({self.type}1 card) = {len(card):d}\ncard={card}'
        #return cls(components, nodes, comment=comment)
        nodes = expand_thru(nodes)
        nnodes = len(nodes)
        components = [component] * nnodes
        self.cards.append((nodes, components, ifile, comment))
        #if comment:
            #self.comment[nid] = comment
        self.n += nnodes
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        idtype = self.model.idtype
        node_id_list = []
        component_list = []
        #comment = {}
        for icard, card in enumerate(self.cards):
            (nidi, componenti, ifilei, commenti) = card
            assert isinstance(nidi, list), nidi
            assert isinstance(componenti, list), componenti
            #nnodes = len(nidi)
            node_id_list.extend(nidi)
            component_list.extend(componenti)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti
        node_id = np.array(node_id_list, dtype=idtype)
        component = np.array(component_list, dtype=idtype)
        self._save(node_id, component, comment=None)
        #self.sort()
        self.cards = []

    def _save(self,
              node_id: np.ndarray,
              component: np.ndarray,
              comment: dict[int, str]=None) -> None:
        #ncards = len(node_id)
        ncards_existing = len(self.node_id)

        if ncards_existing != 0:
            node_id = np.hstack([self.node_id, node_id])
            component = np.hstack([self.component, component])
        #if comment:
            #self.comment.update(comment)
        self.node_id = node_id
        self.component = component
        #print(node_id, component)
        self.n = len(node_id)
        #self.sort()
        #self.cards = []

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id # .ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.node_id)

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    #def slice_card_by_node_id(self, node_id: np.ndarray) -> GRID:
        #"""uses a node_ids to extract GRIDs"""
        #inid = self.index(node_id)
        ##assert len(self.node_id) > 0, self.node_id
        ##i = np.searchsorted(self.node_id, node_id)
        #grid = self.slice_card_by_index(inid)
        #return grid

    #def slice_card_by_index(self, i: np.ndarray) -> GRID:
        #"""uses a node_index to extract GRIDs"""
        #assert self.xyz.shape == self._xyz_cid0.shape
        #assert len(self.node_id) > 0, self.node_id
        #i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        #i.sort()
        #grid = GRID(self.model)
        #self.__apply_slice__(grid, i)
        #return grid

    def __apply_slice__(self, spcoff: SPCOFF, i: np.ndarray) -> None:
        self._slice_comment(spcoff, i)
        spcoff.n = len(i)
        spcoff.node_id = self.node_id[i]
        spcoff.component = self.component[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),)

    @property
    def max_id(self) -> int:
        return self.node_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        component_to_nodes = defaultdict(list)
        for nid, comp in zip(self.node_id, self.component):
            component_to_nodes[comp].append(nid)

        # BNDGRID1 doesn't exist, but BNDGRID has that format
        class_name = self.type if self.type == 'BNDGRID' else f'{self.type}1'
        for component, nodes in component_to_nodes.items():
            nodes_collapsed = collapse_thru(nodes, nthru=None)
            list_fields = [class_name, component] + nodes_collapsed
            bdf_file.write(print_card(list_fields))
        return

    #def index(self, node_id: np.ndarray, safe: bool=False) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #if safe:
            #ibad = inid >= len(self.node_id)
            #if sum(ibad):
                ##self.model.log.error(f'bad nids; node_id={node_id[ibad]}')
                #raise RuntimeError(f'bad nids; node_id={node_id[ibad]}')
            #inids_leftover = inid[~ibad]
            #if len(inids_leftover):
                #actual_nids = self.node_id[inids_leftover]
                #assert np.array_equal(actual_nids, node_id)
        #return inid

class SPCOFF(CommonSet):
    pass
class BNDFIX(CommonSet):
    pass
class BNDFREE(CommonSet):
    pass
class BNDGRID(CommonSet):
    pass
    #def add_card(self, card: BDFCard, ifile: int, comment: str=''):
        #thru = card.field(3).strip()
        #if thru == 'THRU':
            #self.add_set_card(card, comment='')


def spc_cards_to_nid_dof(spc_id: int,
                         cards: list[SPC | SPC1],
                         ) -> tuple[bool, np.ndarray, np.ndarray]:
    """helper for making SPCs/SPCADD node/component"""
    comp_list = []
    nids_list = []
    is_failed = True
    for card in cards:
        if spc_id not in card.spc_id:
            continue
        spc = card.slice_card_by_id(spc_id)

        if len(spc.components) == len(spc.node_id):
            component = spc.components
            comp_list.append(component)
            nids_list.append(spc.node_id)
        else:
            assert spc.type == 'SPC1', spc
            for i, nnode, (inode0, inode1) in zip(count(), spc.nnodes, spc.inode):
                componenti = spc.components[i]
                component = np.ones(nnode, dtype='int32') * componenti
                node_id = spc.node_id[inode0:inode1]
                comp_list.append(component)
                nids_list.append(node_id)

        comp = np.hstack(comp_list)
        nids = np.hstack(nids_list)
        assert len(comp) == len(nids)
        del comp, nids

        #x = 1
    if len(comp_list) == 0 and len(nids_list) == 0:
        return is_failed, nids_list, comp_list

    comp = np.hstack(comp_list)
    nids = np.hstack(nids_list)
    assert len(comp) == len(nids)
    is_failed = False
    return is_failed, nids, comp

SPCs = SPC | SPC1
