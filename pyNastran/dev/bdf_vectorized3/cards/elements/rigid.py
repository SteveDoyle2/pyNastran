from __future__ import annotations
from itertools import count, zip_longest
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import cast_ints # integer_types
#from pyNastran.bdf import MAX_INT
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, hslice_by_idim, make_idim,
    parse_element_check, get_print_card_8_16,
    remove_unused_primary)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, string, blank,
    integer_or_blank, double_or_blank,
    integer_or_double, integer_double_or_blank, integer_double_or_string,
    parse_components, components_or_blank,
)
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16 # print_float_16,
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, array_default_float,
    get_print_card_size)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class RigidElement(VectorizedBaseCard):
    _id_name = 'element_id'
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        pass

    def slice_card_by_element_id(self, element_id: np.ndarray) -> RigidElement:
        assert self.n > 0, self.n
        assert len(self.element_id) > 0, self.element_id
        i = self.index(element_id)
        #cls_obj = cls(self.model)
        #cls_obj.__apply_slice__(self, i)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def sort(self) -> None:
        i = np.argsort(self.element_id)
        self.__apply_slice__(self, i)

    def convert(self, alpha_scale: float=1., **kwargs):
        self.alpha *= alpha_scale


def _str_to_int(value: str) -> int:
    if isinstance(value, str):
        if value == '':
            value = -1
        else:
            value = int(value)
        return value
    assert isinstance(value, int), value
    return value

class RBAR(RigidElement):
    """
    Defines a rigid bar with six degrees-of-freedom at each end.

    +------+-----+----+----+--------+-----+-----+-----+-------+
    |  1   |  2  |  3 |  4 |    5   |  6  |  7  |  8  |   9   |
    +======+=====+====+====+========+=====+=====+=====+=======+
    | RBAR | EID | GA | GB |  CNA   | CNB | CMA | CMB | ALPHA |
    +------+-----+----+----+--------+-----+-----+-----+-------+
    | RBAR |  5  | 1  |  2 | 123456 |     |     |     | 6.5-6 |
    +------+-----+----+----+--------+-----+-----+-----+-------+

    TREF was added in MSC 2021
    """
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def add(self, eid: int, nids: list[int],
            cna: int, cnb: int,
            cma: int, cmb: int,
            alpha: float=0.,
            tref: float=0.,
            comment: str='') -> int:
        """
        Creates a RBAR element

        Parameters
        ----------
        eid : int
            element id
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        cna / cnb : str
            independent DOFs in '123456'
        cma / cmb : str
            dependent DOFs in '123456'
        alpha : float; default=0.0
            coefficient of thermal expansion
        tref : float; default=0.0
            reference temperature
        comment : str; default=''
            a comment for the card

        """
        cna = _str_to_int(cna)
        cnb = _str_to_int(cnb)
        cma = _str_to_int(cma)
        cmb = _str_to_int(cmb)
        self.cards.append((eid, nids, [cna, cnb], [cma, cmb], alpha, tref, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        ga = integer(card, 2, 'ga')
        gb = integer(card, 3, 'gb')
        cna = components_or_blank(card, 4, 'cna', default='0')
        cnb = components_or_blank(card, 5, 'cnb', default='0')
        cma = components_or_blank(card, 6, 'cma', default='0')
        cmb = components_or_blank(card, 7, 'cmb', default='0')
        alpha = double_or_blank(card, 8, 'alpha', default=0.0)
        tref = double_or_blank(card, 9, 'tref', default=0.0)
        assert len(card) <= 10, f'len(RBAR card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, [ga, gb], [cna, cnb], [cma, cmb], alpha, tref, comment))
        self.n += 1
        return self.n - 1

    @RigidElement.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        dependent_dof = np.zeros((ncards, 2), dtype='int32')
        independent_dof = np.zeros((ncards, 2), dtype='int32')
        tref = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, nodesi, cn, cm, alphai, trefi, comment) = card
            element_id[icard] = eid
            nodes[icard, :] = nodesi
            independent_dof[icard, :] = cn
            dependent_dof[icard, :] = cm
            alpha[icard] = alphai
            tref[icard] = trefi
        self._save(element_id, nodes, independent_dof, dependent_dof, alpha, tref)
        self.cards = []

    def _save(self, element_id, nodes, independent_dof, dependent_dof, alpha, tref):
        nelements = len(element_id)
        if alpha is None:
            alpha = np.zeros(nelements, dtype='float64')
        if tref is None:
            tref = np.zeros(nelements, dtype='float64')

        if len(self.element_id):
            element_id = np.hstack([self.element_id, element_id])
            nodes = np.vstack([self.nodes, nodes])
            independent_dof = np.vstack([self.independent_dof, independent_dof])
            dependent_dof = np.vstack([self.dependent_dof, dependent_dof])
            alpha = np.hstack([self.alpha, alpha])
            tref = np.hstack([self.tref, tref])
        self.element_id = element_id
        self.nodes = nodes
        self.independent_dof = independent_dof
        self.dependent_dof = dependent_dof
        self.alpha = alpha
        self.tref = tref
        self.n = len(element_id)

    def __apply_slice__(self, element: RBAR, i: np.ndarray):
        element.element_id = self.element_id[i]
        element.nodes = self.nodes[i, :]
        element.independent_dof = self.independent_dof[i]
        element.dependent_dof = self.dependent_dof[i]
        element.alpha = self.alpha[i]
        element.tref = self.tref[i]
        element.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        eids = array_str(self.element_id)
        nodes = array_str(self.nodes, size=size)
        ind_nodes = array_default_int(self.independent_dof, default=-1, size=size)
        dep_nodes = array_default_int(self.dependent_dof, default=-1, size=size)

        alphas = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
        trefs = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
        for eid, nodes, independent_dof, dependent_dof, alpha, tref in zip(eids, nodes, ind_nodes, dep_nodes, alphas, trefs):
            ga, gb = nodes
            cna, cnb = independent_dof
            cma, cmb = dependent_dof

            list_fields = ['RBAR', eid, ga, gb, cna, cnb,
                           cma, cmb, alpha, tref]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        #all_nodes = np.intersect1d(self.independent_node, self.dependent_nodes, assume_unique=False)
        #geom_check(self,
                   #missing,
                   #node=(nid, all_nodes))


class RROD(RigidElement):
    """
    Rigid Pin-Ended Element Connection
    Defines a pin-ended element that is rigid in translation

    +------+-----+----+----+-----+-----+-------+
    |   1  |  2  | 3  | 4  |  5  |  6  |   7   |
    +======+=====+====+====+=====+=====+=======+
    | RROD | EID | GA | GB | CMA | CMB | ALPHA |
    +------+-----+----+----+-----+-----+-------+
    | RROD | 5   | 1  |  2 |     |     | 6.5-6 |
    +------+-----+----+----+-----+-----+-------+
    """
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def add(self, eid, nids, cma='', cmb='', alpha=0.0, comment='') -> int:
        if cma is None or cma == '':
            cma = 0
        else:
            cma = int(cma)

        if cmb is None or cmb == '':
            cmb = 0
        else:
            cmb = int(cmb)
        self.cards.append((eid, nids, [cma, cmb], alpha, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a RROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        ga = integer(card, 2, 'ga')
        gb = integer(card, 3, 'gb')
        cma = components_or_blank(card, 4, 'cma', default=0)
        cmb = components_or_blank(card, 5, 'cmb', default=0)
        alpha = double_or_blank(card, 6, 'alpha', default=0.0)
        assert len(card) <= 7, f'len(RROD card) = {len(card):d}\ncard={card}'
        #return RROD(eid, [ga, gb], cma, cmb, alpha, comment=comment)
        self.cards.append((eid, [ga, gb], [cma, cmb], alpha, comment))
        self.n += 1
        return self.n - 1

    @RigidElement.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        dependent_dof = np.zeros((ncards, 2), dtype='int32')
        #independent_dof = np.zeros((ncards, 2), dtype='int32')
        #tref = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, nodesi, cm, alphai, comment) = card
            element_id[icard] = eid
            nodes[icard, :] = nodesi
            dependent_dof[icard] = cm
            alpha[icard] = alphai
            #self.tref[icard] = tref
        self._save(element_id, nodes, dependent_dof, alpha)
        self.cards = []

    def _save(self, element_id, nodes, cm, alpha) -> None:
        if len(self.element_id):
            asdf
        nelements = len(element_id)
        assert element_id.min() > 0, element_id
        assert nodes.min() > 0, nodes
        if alpha is None:
            alpha = np.zeros(nelements, dtype='float64')

        self.element_id = element_id
        self.nodes = nodes
        self.dependent_dof = cm
        self.alpha = alpha
        self.n = nelements

    def __apply_slice__(self, element: RROD, i: np.ndarray):
        element.element_id = self.element_id[i]
        element.nodes = self.nodes[i, :]
        element.dependent_dof = self.dependent_dof[i]
        element.alpha = self.alpha[i]
        element.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        eids = array_str(self.element_id)
        nodes = array_str(self.nodes, size=size)
        #ind_nodes = array_default_int(self.independent_dof, default=0, size=size)
        dep_nodes = array_default_int(self.dependent_dof, default=0, size=size)
        alphas = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
        for eid, nodes, dependent_dof, alpha in zip(eids, nodes, dep_nodes, alphas):
            ga, gb = nodes
            #cna, cnb = independent_dof
            cma, cmb = dependent_dof

            list_fields = ['RROD', eid, ga, gb,
                           cma, cmb, alpha]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        #all_nodes = np.intersect1d(self.independent_node, self.dependent_nodes, assume_unique=False)
        #geom_check(self,
                   #missing,
                   #node=(nid, all_nodes))


class RBAR1(RigidElement):
    """
    +-------+-----+----+----+-----+-------+
    |   1   |  2  |  3 |  4 |  5  |   6   |
    +=======+=====+====+====+=====+=======+
    | RBAR1 | EID | GA | GB | CB  | ALPHA |
    +-------+-----+----+----+-----+-------+
    | RBAR1 | 5   |  1 |  2 | 123 | 6.5-6 |
    +-------+-----+----+----+-----+-------+
    """
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def add(self, eid: int, nids: list[int],
            cb: str='123456', alpha: float=0.0, comment: str='') -> int:
        #if cma is None or cma == '':
            #cma = 0
        #else:
            #cma = int(cma)

        #if cmb is None or cmb == '':
            #cmb = 0
        #else:
            #cmb = int(cmb)
        self.cards.append((eid, nids, cb, alpha, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a RBAR1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        ga = integer(card, 2, 'ga')
        gb = integer(card, 3, 'gb')
        cb = components_or_blank(card, 4, 'cb', default='123456')
        alpha = double_or_blank(card, 5, 'alpha', default=0.0)
        assert len(card) <= 6, f'len(RBAR1 card) = {len(card):d}\ncard={card}'
        #return RBAR1(eid, [ga, gb], cb, alpha=alpha, comment=comment)
        self.cards.append((eid, [ga, gb], cb, alpha, comment))
        self.n += 1
        return self.n - 1

    @RigidElement.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        dependent_dof = np.zeros(ncards, dtype='int32')
        #independent_dof = np.zeros((ncards, 2), dtype='int32')
        #tref = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, nodesi, cb, alphai, comment) = card
            element_id[icard] = eid
            nodes[icard, :] = nodesi
            dependent_dof[icard] = cb
            alpha[icard] = alphai
            #self.tref[icard] = tref
        self._save(element_id, nodes, dependent_dof, alpha)
        self.cards = []

    def _save(self, element_id, nodes, cb, alpha) -> None:
        if len(self.element_id):
            asdf
        nelements = len(element_id)
        assert element_id.min() > 0, element_id
        assert nodes.min() > 0, nodes
        if alpha is None:
            alpha = np.zeros(nelements, dtype='float64')

        self.element_id = element_id
        self.nodes = nodes
        self.dependent_dof = cb
        self.alpha = alpha
        self.n = nelements

    def __apply_slice__(self, element: RBAR1, i: np.ndarray):
        element.element_id = self.element_id[i]
        element.nodes = self.nodes[i, :]
        element.dependent_dof = self.dependent_dof[i]
        element.alpha = self.alpha[i]
        element.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        eids = array_str(self.element_id)
        nodes = array_str(self.nodes, size=size)
        dep_nodes = array_default_int(self.dependent_dof, default=0, size=size)
        alphas = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
        for eid, nodes, cb, alpha in zip(eids, nodes, dep_nodes, alphas):
            ga, gb = nodes
            list_fields = ['RBAR1', eid, ga, gb, cb, alpha]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        #all_nodes = np.intersect1d(self.independent_node, self.dependent_nodes, assume_unique=False)
        #geom_check(self,
                   #missing,
                   #node=(nid, all_nodes))


class RBE1(RigidElement):
    """
    +------+-----+-----+-----+-------+-----+-----+-----+
    |   1  |  2  |  3  |  4  |   5   |  6  |  7  |  8  |
    +======+=====+=====+=====+=======+=====+=====+=====+
    | RBE1 | EID | GN1 | CN1 |  GN2  | CN2 | GN3 | CN3 |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      |     | GN4 | CN4 |  GN5  | CN5 | GN6 | CN6 |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      | UM  | GM1 | CM1 |  GM2  | CM2 | GM3 | CM3 |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      | GM4 | CM4 | etc | ALPHA |     |     |     |
    +------+-----+-----+-----+-------+-----+-----+-----+

    +------+-----+-----+-----+-------+-----+-----+-----+
    | RBE1 | 59  | 59  | 123 |  60   | 456 |     |     |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      | UM  | 61  | 246 |       |     |     |     |
    +------+-----+-----+-----+-------+-----+-----+-----+
    """
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def add(self, eid: int,
            Gni: list[int], Cni: list[int],
            Gmi: list[int], Cmi: list[int],
            alpha: float=0., comment: str='') -> int:
        self.cards.append((eid, alpha, Gni, Cni, Gmi, Cmi, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        ium = card.index('UM')
        if ium > 0:
            assert string(card, ium, 'UM') == 'UM', 'RBE1=%s must contain UM' % str(card)

        #assert isinstance(card[-1], str), 'card[-1]=%r type=%s' %(card[-1], type(card[-1]))
        alpha_last = integer_double_or_string(card, -1, 'alpha_last')
        if isinstance(alpha_last, float):
            alpha = alpha_last
            card.pop()  # remove the last field so len(card) will not include alpha
        else:
            alpha = 0.

        # loop till UM, no field9,field10
        n = 1
        i = 0
        offset = 2
        Gni = []
        Cni = []
        Gmi = []
        Cmi = []
        while offset + i < ium - 1:
            #print('field(%s) = %s' % (offset + i, card.field(offset + i)))
            gni = integer_or_blank(card, offset + i, 'gn%i' % n)
            cni = components_or_blank(card, offset + i + 1, 'cn%i' % n)

            if gni:
                #print("gni=%s cni=%s" % (gni ,cni))
                Gni.append(gni)
                Cni.append(cni)
                n += 1
            else:
                assert cni is None
            i += 2

        # loop till alpha, no field9,field10
        n = 1
        offset = ium + 1
        i = 0

        # dont grab alpha
        while offset + i < len(card):
            gmi = integer_or_blank(card, offset + i, 'gm%i' % n)
            cmi = components_or_blank(card, offset + i + 1, 'cm%i' % n)
            if gmi:
                Gmi.append(gmi)
                Cmi.append(cmi)
                n += 1
            else:
                assert cmi is None
            i += 2
        self.cards.append((eid, alpha, Gni, Cni, Gmi, Cmi, comment))
        self.n += 1
        return self.n - 1

    @RigidElement.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        ndependent = np.zeros(ncards, dtype='int32')
        nindependent = np.zeros(ncards, dtype='int32')
        tref = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')

        all_dependent_nodes = []
        all_dependent_dofs = []
        all_independent_nodes = []
        all_independent_dofs = []
        for icard, card in enumerate(self.cards):
            (eid, alphai, Gni, Cni, Gmi, Cmi, comment) = card

            #def independent_nodes(self):
                #return self.Gni_node_ids
            #def dependent_nodes(self):
                #return self.Gmi_node_ids

            element_id[icard] = eid
            ndependent[icard] = len(Gni)
            nindependent[icard] = len(Gmi)

            all_dependent_nodes.extend(Gni)
            all_dependent_dofs.extend(Cni)
            all_independent_nodes.extend(Gmi)
            all_independent_dofs.extend(Cmi)
            alpha[icard] = alphai

        dependent_node = np.array(all_dependent_nodes, dtype=idtype)
        dependent_dof = np.array(all_dependent_dofs, dtype='int32')
        independent_node = np.array(all_independent_nodes, dtype=idtype)
        independent_dof = np.array(all_independent_dofs, dtype='int32')
        self._save(element_id,
                   ndependent, dependent_node, dependent_dof,
                   nindependent, independent_node, independent_dof,
                   tref, alpha)
        self.cards = []

    def _save(self, element_id,
              ndependent, dependent_node, dependent_dof,
              nindependent, independent_node, independent_dof,
              tref, alpha):
        self.element_id = element_id

        self.ndependent = ndependent
        self.dependent_node = dependent_node
        self.dependent_dof = dependent_dof

        self.nindependent = nindependent
        self.independent_node = independent_node
        self.independent_dof = independent_dof
        self.tref = tref
        self.alpha = alpha

    def __apply_slice__(self, element: RBAR, i: np.ndarray):
        element.element_id = self.element_id[i]
        idependent = self.idependent
        element.dependent_node = hslice_by_idim(i, idependent, self.dependent_node)
        element.dependent_dof = hslice_by_idim(i, idependent, self.dependent_dof)
        element.ndependent = self.ndependent[i]

        iindependent = self.iindependent
        element.independent_node = hslice_by_idim(i, iindependent, self.independent_node)
        element.independent_dof = hslice_by_idim(i, iindependent, self.independent_dof)
        element.nindependent = self.nindependent[i]

        element.tref = self.tref[i]
        element.alpha = self.alpha[i]
        element.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.dependent_node)
        used_dict['node_id'].append(self.independent_node)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        for nodes in (self.dependent_node, self.independent_node):
            for i, nid1 in enumerate(nodes):
                nid2 = nid_old_to_new.get(nid1, nid1)
                nodes[i] = nid2

    @property
    def idependent(self) -> np.ndarray:
        return make_idim(self.n, self.ndependent)

    @property
    def iindependent(self) -> np.ndarray:
        return make_idim(self.n, self.nindependent)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(),
                   self.independent_node.max(),
                   self.dependent_node.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        eids = array_str(self.element_id)
        #ind_nodes = array_default_int(self.independent_dof, default=0, size=size)
        #dep_nodes = array_default_int(self.dependent_dof, default=0, size=size)

        for eid, idependent, iindependent, alpha in zip(eids, self.idependent, self.iindependent, self.alpha):
            idependent0, idependent1 = idependent
            iindependent0, iindependent1 = iindependent
            Gmi = self.independent_node[iindependent0:iindependent1]
            Cmi = self.independent_dof[iindependent0:iindependent1]

            Gni = self.dependent_node[idependent0:idependent1]
            Cni = self.dependent_dof[idependent0:idependent1]

            list_fields = ['RBE1', eid]
            for (i, gn, cn) in zip(count(), Gni, Cni):
                list_fields += [gn, cn]
                if i > 0 and i % 3 == 0:
                    list_fields += [None]

            nspaces = 8 - (len(list_fields) - 1) % 8  # puts UM/ALPHA onto next line
            if nspaces < 8:
                list_fields += [None] * nspaces

            # overly complicated loop to print the UM section
            list_fields += ['UM']
            j = 1
            for (i, gm, cm) in zip(count(), Gmi, Cmi):
                list_fields += [gm, cm]
                if i > 0 and j % 3 == 0:
                    list_fields += [None, None]
                    j -= 3
                j += 1

            if alpha > 0.:  # handles default alpha value
                nspaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
                if nspaces == 1:
                    list_fields += [None, None]
                list_fields += [alpha]

            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        #all_nodes = np.intersect1d(self.independent_node, self.dependent_nodes, assume_unique=False)
        #geom_check(self,
                   #missing,
                   #node=(nid, all_nodes))


class RBE2(RigidElement):
    """
    +-------+-----+-----+-----+------+-------+------+-----+-----+
    |   1   |  2  |  3  |  4  |  5   |   6   |  7   |  8  |  9  |
    +=======+=====+=====+=====+======+=======+======+=====+=====+
    |  RBE2 | EID | GN  | CM  | GM1  |  GM2  | GM3  | GM4 | GM5 |
    +-------+-----+-----+-----+------+-------+------+-----+-----+
    |       | GM6 | GM7 | GM8 | etc. | ALPHA | TREF |     |     |
    +-------+-----+-----+-----+------+-------+------+-----+-----+

    TREF was added in MSC 2021
    """
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def __apply_slice__(self, element: RBE2, i: np.ndarray):
        element.element_id = self.element_id[i]
        element.independent_node = self.independent_node[i]
        element.independent_dof = self.independent_dof[i]
        element.tref = self.tref[i]
        element.dependent_nodes = hslice_by_idim(i, self.idim, self.dependent_nodes)
        element.nnode = self.nnode[i]
        element.n = len(i)

    def add(self, eid: int,
            gn: int, # independent
            cm: str, Gmi: list[int],  # dependent
            alpha: float=0.0, tref: float=0.0, comment: str='') -> int:
        self.cards.append((eid, gn, cm, Gmi, alpha, tref, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        gn = integer(card, 2, 'gn')
        cm = components_or_blank(card, 3, 'cm', default=0)

        tref = 0.0
        gm_alpha_tref = integer_or_double(card, len(card) - 1, 'gm/alpha/tref')
        if isinstance(gm_alpha_tref, float):
            gm_alpha = integer_double_or_blank(card, len(card) - 2, 'gm/alpha', default=0.0)
            if isinstance(gm_alpha, float):
                # alpha/tref is correct
                # the last field is tref
                n = 2
                alpha = gm_alpha
                tref = gm_alpha_tref
            else:
                # alpha is correct
                # the last field is alpha
                n = 1
                alpha = gm_alpha_tref
        else:
            # the last field is Gm
            n = 0
            alpha = 0.0

        j = 4
        Gmi = []
        for i in range(len(card) - 4 - n):
            gmi = integer(card, j + i, 'Gm%i' % (i + 1))
            Gmi.append(gmi)
        self.cards.append((eid, gn, cm, Gmi, alpha, tref, comment))
        self.n += 1
        return self.n - 1

    @RigidElement.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        independent_node = np.zeros(ncards, dtype=idtype)
        independent_dof = np.zeros(ncards, dtype='int32')
        tref = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')
        nnode = np.zeros(ncards, dtype='int32')

        dependent_nodes = []
        for icard, card in enumerate(self.cards):
            eid, gn, cm, Gmi, alphai, trefi, comment = card
            dependent_nodes.extend(Gmi)
            nnode[icard] = len(Gmi)

            element_id[icard] = eid
            independent_node[icard] = gn
            independent_dof[icard] = cm
            alpha[icard] = alphai
            tref[icard] = trefi
        dependent_nodes = np.array(dependent_nodes, dtype=idtype)

        self._save(element_id, independent_node, independent_dof, dependent_nodes, nnode, alpha, tref)
        self.sort()
        self.cards = []

    def _save(self, element_id, independent_node, independent_dof, dependent_nodes, nnode, alpha, tref):
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            independent_node = np.hstack([self.independent_node, independent_node])
            independent_dof = np.hstack([self.independent_dof, independent_dof])
            dependent_nodes = np.hstack([self.dependent_nodes, dependent_nodes])
            nnode = np.hstack([self.nnode, nnode])
            alpha = np.hstack([self.alpha, alpha])
            tref = np.hstack([self.tref, tref])

        self.element_id = element_id
        self.independent_dof = independent_dof
        self.independent_node = independent_node
        self.dependent_nodes = dependent_nodes
        self.nnode = nnode
        self.alpha = alpha
        self.tref = tref

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.dependent_nodes.ravel())
        used_dict['node_id'].append(self.independent_node)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        for nodes in (self.independent_node, self.dependent_nodes):
            for i, nid1 in enumerate(nodes):
                nid2 = nid_old_to_new.get(nid1, nid1)
                nodes[i] = nid2

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnode)

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        eid_str = array_default_int(self.element_id, default=0, size=size)
        ind_node_str = array_default_int(self.independent_node, default=0, size=size)
        ind_components_str = array_default_int(self.independent_dof, default=0, size=size)

        dep_node_str = array_str(self.dependent_nodes, size=size)

        no_alpha = self.alpha.max() == 0. and self.alpha.min() == 0.
        no_tref = self.tref.max() == 0. and self.tref.min() == 0.
        if no_tref and no_alpha:
            for eid, ind_node, ind_component, idim in zip_longest(eid_str, ind_node_str, ind_components_str, self.idim):
                idim0, idim1 = idim
                dep_nodes = dep_node_str[idim0:idim1].tolist()  # Gmi
                list_fields = ['RBE2', eid, ind_node, ind_component] + dep_nodes
                bdf_file.write(print_card(list_fields))
        elif no_tref:
            alphas = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
            for eid, ind_node, ind_component, idim, alpha in zip_longest(eid_str, ind_node_str, ind_components_str,
                                                                         self.idim, alphas):
                idim0, idim1 = idim
                dep_nodes = dep_node_str[idim0:idim1].tolist()  # Gmi

                list_fields = ['RBE2', eid, ind_node, ind_component] + dep_nodes + [alpha]
                bdf_file.write(print_card(list_fields))
        else:
            alphas = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
            trefs = array_default_float(self.alpha, default=0., size=size, is_double=is_double)
            for eid, ind_node, ind_component, idim, alpha, tref in zip_longest(eid_str, ind_node_str, ind_components_str,
                                                                               self.idim, alphas, trefs):
                idim0, idim1 = idim
                dep_nodes = dep_node_str[idim0:idim1].tolist()  # Gmi

                list_fields = ['RBE2', eid, ind_node, ind_component] + dep_nodes + [alpha, tref]
                bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        all_nodes = np.intersect1d(self.independent_node, self.dependent_nodes, assume_unique=False)
        geom_check(self,
                   missing,
                   node=(nid, all_nodes))


class RBE3(RigidElement):
    """
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |   1  |    2    |    3    |    4    |  5   |    6   |    7   |   8  |    9   |
    +======+=========+=========+=========+======+========+========+======+========+
    | RBE3 |   EID   |         | REFGRID | REFC |  WT1   |   C1   | G1,1 |  G1,2  |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   G1,3  |   WT2   |   C2    | G2,1 |  G2,2  |  etc.  | WT3  |   C3   |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   G3,1  |   G3,2  |  etc.   | WT4  |  C4    |  G4,1  | G4,2 |  etc.  |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   'UM'  |   GM1   |   CM1   | GM2  |  CM2   |  GM3   | CM3  |        |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   GM4   |   CM4   |   GM5   | CM5  |  etc.  |        |      |        |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      | 'ALPHA' |   ALPHA |         |      |        |        |      |        |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    """
    def clear(self) -> None:
        self.n = 0
        self.element_id = np.array([])

    def __apply_slice__(self, element: RBE2, i: np.ndarray):
        element.element_id = self.element_id[i]
        element.independent_node = self.independent_node[i]
        element.independent_dof = self.independent_dof[i]
        element.tref = self.tref[i]
        element.nodes = hslice_by_idim(i, self.inode, self.nodes)
        element.nnode = self.nnode[i]

    def add(self, eid: int, refgrid: int, refc: str,
            weights: list[float], comps: list[str], Gijs: list[int],
            Gmi=None, Cmi=None,
            alpha: float=0.0, tref: float=0.0,
            comment: str='') -> int:
        if Gmi is None or Cmi is None:
            assert Gmi is None and Cmi is None, 'Gmi and Cmi should both be None if either is None'
            Gmi = []
            Cmi = []
        self.cards.append((eid, refgrid, refc, alpha, tref,
                           weights, Gijs, comps,
                           Gmi, Cmi, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        blank(card, 2, 'blank')
        refgrid = integer(card, 3, 'refgrid')
        refc = components_or_blank(card, 4, 'refc')

        fields = [field.upper() if isinstance(field, str) else field for field in card[5:]]
        ioffset = 5
        iwt_max = len(fields) + ioffset
        try:
            ialpha = fields.index('ALPHA') + ioffset
            iwt_max = ialpha  # the index to start parsing UM
            ium_stop = ialpha  # the index to stop  parsing UM
        except ValueError:
            ialpha = None
            ium_stop = iwt_max

        try:
            ium = fields.index('UM') + ioffset
            iwt_max = ium
        except ValueError:
            ium = None

        i = ioffset
        n = 1
        weights = []
        comps = []
        Gijs = []

        while i < iwt_max:
            Gij = []
            wtname = 'wt' + str(n)
            wt = double_or_blank(card, i, wtname)
            if wt is not None:
                cname = 'c' + str(n)
                compi = components_or_blank(card, i + 1, cname)

                #print("%s=%s %s=%s" % (wtname, wt, cname, compi))
                i += 2
                gij = 0

                j = 0
                while isinstance(gij, int) and i < iwt_max:
                    j += 1
                    gij_name = 'g%s,%s' % (n, j)
                    gij = integer_double_or_blank(card, i, gij_name)
                    if isinstance(gij, float):
                        break
                    #print("%s = %s" % (gij_name, gij))
                    if gij is not None:
                        Gij.append(gij)
                    i += 1
                assert compi is not None
                assert len(Gij) > 0, Gij
                assert Gij[0] is not None, Gij
                weights.append(wt)

                # original
                comps.append(compi)
                Gijs.append(Gij)
                # new
                #print('----finished a group=%r----' % weight_cg_group)
            else:
                i += 1

        Gmi = []
        Cmi = []
        if ium:
            #print('UM = %s' % card.field(ium))  # UM
            i = ium + 1
            n = 1
            #print("i=%s iUmStop=%s" % (i, iUmStop))
            for j in range(i, ium_stop, 2):

                gm_name = 'gm' + str(n)
                cm_name = 'cm' + str(n)
                gmi = integer_or_blank(card, j, gm_name)
                if gmi is not None:
                    cmi = parse_components(card, j + 1, cm_name)
                    #print("gmi=%s cmi=%s" % (gmi, cmi))
                    Gmi.append(gmi)
                    Cmi.append(cmi)

        if ialpha:
            alpha = double_or_blank(card, ialpha + 1, 'alpha', default=0.0)
            tref = double_or_blank(card, ialpha + 2, 'tref', default=0.0)
        else:
            #: thermal expansion coefficient
            alpha = 0.0
            tref = 0.0

        #def independent_nodes(self) -> list[int]:
            #for gij in self.Gijs:

        #def dependent_nodes(self) -> list[int]:
            #nodes = [self.ref_grid_id]
            #nodes += self.Gmi_node_ids

        #self.nnode[icard] = len(Gmi)
        self.cards.append((eid, refgrid, refc, alpha, tref,
                           weights, Gijs, comps,
                           Gmi, Cmi, comment))
        self.n += 1
        return self.n - 1

    @RigidElement.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        # basic
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        ref_grid = np.zeros(ncards, dtype=idtype)
        ref_component = np.zeros(ncards, dtype='int32')
        tref = np.zeros(ncards, dtype='float64')
        alpha = np.zeros(ncards, dtype='float64')
        #self.nnode = np.zeros(ncards, dtype='int32')

        # weight/components/Gij
        nweight = np.zeros(ncards, dtype='int32')
        ndependent = np.zeros(ncards, dtype='int32')
        #self.independent_node = np.zeros(ncards, dtype='int32')
        #self.independent_dof = np.zeros(ncards, dtype='int32')

        # dependent
        #dependent_nodes = []
        #dependent_dofs = []

        all_weights = []
        independent_nodes = [] # Gij
        independent_dofs = []  # comps

        dependent_nodes = [] # Gmi
        dependent_dofs = []  # Cmi

        ngrid_per_weight = []
        for icard, card in enumerate(self.cards):
            (eid, refgrid, refc, alphai, trefi,
             weights, Gijs, comps,
             Gmi, Cmi, comment) = card
            # basic
            element_id[icard] = eid
            ref_grid[icard] = refgrid
            ref_component[icard] = refc
            alpha[icard] = alphai
            tref[icard] = trefi

            #assert len(Gij) > 0
            #assert len(compi) > 0
            for Gij in Gijs:
                independent_nodes.extend(Gij)
                ngrid_per_weight.append(len(Gij))
            for compi in comps:
                independent_dofs.append(compi)

            # Gij
            nweight[icard] = len(weights)
            all_weights.extend(weights)

            ndependent[icard] = len(Gmi)
            dependent_nodes.extend(Gmi)
            dependent_dofs.extend(Cmi)

        #element_id = cast_ints(element_id, dtype='int32')
        #ref_grid = cast_ints(ref_grid, dtype='int32')

        weight = np.array(all_weights, dtype='float64')
        independent_nodes = cast_ints(independent_nodes, dtype='int32')
        independent_dofs = np.array(independent_dofs, dtype='int32')
        ngrid_per_weight = np.array(ngrid_per_weight, dtype='int32')

        dependent_nodes = cast_ints(dependent_nodes, dtype='int32')
        dependent_dofs = np.array(dependent_dofs, dtype='int32')

        #self.element_id = np.zeros(ncards, dtype='int32')
        #self.ref_grid = np.zeros(ncards, dtype='int32')
        #self.ref_component = np.zeros(ncards, dtype='int32')
        #self.tref = np.zeros(ncards, dtype='float64')
        #self.alpha = np.zeros(ncards, dtype='float64')
        #self.nweight = np.zeros(ncards, dtype='int32')
        #self.ndependent = np.zeros(ncards, dtype='int32')
        self._save(element_id, ref_grid, ref_component, tref, alpha,
                   nweight, ndependent,
                   weight, independent_nodes, independent_dofs,
                   ngrid_per_weight,
                   dependent_nodes, dependent_dofs)
        self.cards = []

    def _save(self, element_id,
              ref_grid, ref_component,
              tref, alpha,
              nweight, ndependent,
              weight, independent_nodes, independent_dofs,
              ngrid_per_weight, dependent_nodes, dependent_dofs) -> None:
        if len(self.element_id):
            asdf

        self.element_id = element_id
        self.ref_grid = ref_grid
        self.ref_component = ref_component
        self.tref = tref
        self.alpha = alpha

        self.nweight = nweight
        self.weight = weight

        self.independent_nodes = independent_nodes
        self.independent_dofs = independent_dofs
        self.ngrid_per_weight = ngrid_per_weight

        self.ndependent = ndependent
        self.dependent_nodes = dependent_nodes
        self.dependent_dofs = dependent_dofs

    def __apply_slice__(self, element: RBE3, i: np.ndarray):
        element.element_id = self.element_id[i]
        element.ref_grid = self.ref_grid[i]
        element.elemeref_componentt_id = self.ref_component[i]
        element.tref = self.tref[i]
        element.alpha = self.alpha[i]

        iweight = self.iweight
        element.weight = hslice_by_idim(i, iweight, self.weight)
        element.independent_dofs = hslice_by_idim(i, iweight, self.independent_dofs)
        element.independent_nodes = hslice_by_idim(i, iweight, self.independent_nodes)
        element.ngrid_per_weight = self.ngrid_per_weight[i]
        element.nweight = self.nweight[i]

        idependent = self.idependent
        element.dependent_dofs = hslice_by_idim(i, idependent, self.dependent_dofs)
        element.dependent_nodes = hslice_by_idim(i, idependent, self.dependent_nodes)

        element.ndependent = self.ndependent[i]
        element.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.dependent_nodes)
        used_dict['node_id'].append(self.independent_nodes)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        for nodes in (self.ref_grid, self.independent_nodes, self.dependent_nodes):
            for i, nid1 in enumerate(nodes):
                nid2 = nid_old_to_new.get(nid1, nid1)
                nodes[i] = nid2

    @property
    def iweight(self) -> np.ndarray:
        return make_idim(self.n, self.nweight)

    @property
    def idependent(self) -> np.ndarray:
        return make_idim(self.n, self.ndependent)

    @property
    def is_small_field(self):
        max_dependent_node = 0
        if self.ndependent.max() > 0:
            max_dependent_node = self.dependent_nodes.max()

        return max(self.element_id.max(),
                   self.independent_nodes.max(),
                   max_dependent_node) < 99_999_999

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        #print_card = get_print_card_8_16(size)
        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        eid_str = array_default_int(self.element_id, default=0, size=size)
        #ind_node_str = array_default_int(self.independent_node, default=0, size=size)
        #ind_components_str = array_default_int(self.independent_dof, default=0, size=size)

        #dep_node_str = array_str(self.dependent_nodes, size=size)

        #no_alpha = self.alpha.max() == 0. and self.alpha.min() == 0.
        #no_tref = self.tref.max() == 0. and self.tref.min() == 0.
        iweight0_for_grids = 0
        #iGij0 = 0
        igrid0 = 0
        for eid, ref_grid, refc, \
                nweight, iweight, \
                idependent, \
                alpha, tref in zip_longest(eid_str, self.ref_grid, self.ref_component,
                                           self.nweight, self.iweight,
                                           self.idependent,
                                           self.alpha, self.tref):
            #idim0, idim1 = idim
            #dep_nodes = dep_node_str[idim0:idim1].tolist()  # Gmi

            # independent
            iweight0, iweight1 = iweight
            weights = self.weight[iweight0:iweight1]
            comps = self.independent_dofs[iweight0:iweight1]
            Gijs = []
            for i in range(nweight):
                ngrid = self.ngrid_per_weight[iweight0_for_grids+i]
                igrid1 = igrid0 + ngrid
                Gij = self.independent_nodes[igrid0:igrid1].tolist()
                Gijs.append(Gij)
                igrid0 = igrid1
            assert len(weights) == len(Gijs)
            assert len(weights) == len(comps)
            iweight0_for_grids += nweight
            #Gij = self.independent_nodes[iweight0:iweight1]
            #igrid = make_idim(len(Gij), ngrid_per_weight)

            # dependent
            idependent0, idependent1 = idependent
            Gmi = self.dependent_nodes[idependent0:idependent1]
            Cmi = self.dependent_dofs[idependent0:idependent1]

            list_fields = ['RBE3', eid, None, ref_grid, refc]
            for (wt, ci, Gij) in zip(weights, comps, Gijs):
                list_fields += [wt, ci] + Gij
            nspaces = 8 - (len(list_fields) - 1) % 8  # puts UM onto next line

            if nspaces < 8:
                list_fields += [None] * nspaces

            if len(Gmi):
                list_fields += ['UM']
                for (gmi, cmi) in zip(Gmi, Cmi):
                    list_fields += [gmi, cmi]

            nspaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
            if nspaces < 8:
                list_fields += [None] * nspaces

            is_alpha = (alpha != 0.0)
            is_tref = (tref != 0.0)
            if is_alpha or is_tref:  # handles the default value
                list_fields += ['ALPHA', alpha]
                if is_tref:
                    list_fields.append(tref)
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        all_nodes = np.intersect1d(self.independent_nodes, self.dependent_nodes, assume_unique=False)
        geom_check(self,
                   missing,
                   node=(nid, all_nodes))

RSSCON = RBAR
