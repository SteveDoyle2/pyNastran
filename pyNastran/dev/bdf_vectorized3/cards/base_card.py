from __future__ import annotations
from abc import abstractmethod
from functools import wraps
from io import StringIO
from typing import Callable, Optional, Any, TYPE_CHECKING
import numpy as np
from pyNastran.dev.bdf_vectorized3.utils import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
#from pyNastran.bdf.bdf_interface.assign_type import (
    #integer, integer_or_blank, double_or_blank, components_or_blank)
#from pyNastran.dev.bdf_vectorized3.cards.line_elements import CBAR, CBEAM, PBAR, PBARL, PBEAM, PBEAML
#from pyNastran.dev.bdf_vectorized3.cards.shell_elements import CQUAD4, CTRIA3, PSHELL, PCOMP, PCOMPG
#from pyNastran.dev.bdf_vectorized3.cards.solid_elements import CTETRA, CHEXA, CPENTA, CPYRAM
#from pyNastran.dev.bdf_vectorized3.cards.static_loads import PLOAD4, FORCE
from pyNastran.utils import object_stats, object_attributes # , object_methods,
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


def make_idim(n: int, ndim: np.ndarray) -> np.ndarray:
    #assert n > 0, n
    idim = np.zeros((n, 2), dtype=ndim.dtype)
    assert ndim.dtype.name in {'int32', 'int64'}, ndim.dtype.name
    if len(ndim) == 0:
        return idim

    assert isinstance(ndim, np.ndarray), ndim

    # if this fails, self.n probably wasn't set right; should be the length of self.sid
    assert n == len(ndim), f'n={n:d} ndim={ndim}'

    csum = np.cumsum(ndim)
    idim[1:, 0] = csum[:-1]
    idim[:, 1] = csum
    assert idim[0, 0] == 0, idim
    assert idim.ndim == 2, idim
    return idim

def hslice_by_idim(i: np.ndarray, idim: np.ndarray, elements: np.ndarray) -> np.ndarray:
    elements2 = []
    assert len(i) > 0, i
    for ii in i:
        idim0, idim1 = idim[ii]
        ndim = idim1 - idim0
        eids = elements[idim0:idim1]
        assert len(eids) == ndim
        elements2.append(eids)
    return np.hstack(elements2)

def vslice_by_idim(i: np.ndarray, idim: np.ndarray, elements: np.ndarray) -> np.ndarray:
    elements2 = []
    for ii in i:
        idim0, idim1 = idim[ii]
        ndim = idim1 - idim0
        eids = elements[idim0:idim1, :]
        assert len(eids) == ndim
        elements2.append(eids)
    return np.vstack(elements2)

def searchsorted_filter_(all_ids, lookup_ids, msg: str='') -> np.ndarray:
    all_ids_array = np.array(all_ids)
    lookup_ids_array = np.array(lookup_ids)
    return searchsorted_filter(all_ids_array, lookup_ids_array, msg=msg)

def searchsorted_filter(all_ids: np.ndarray,
                        lookup_ids: np.ndarray,
                        msg: str='',
                        debug: bool=False) -> tuple[np.ndarray, np.ndarray]:
    """finds the index (like searchsorted) and filters out the incorrect ids"""
    if not np.array_equal(np.unique(all_ids), all_ids):
        raise RuntimeError(f'{msg}={all_ids} is unsorted')
    if lookup_ids.dtype.name == 'object':
        raise TypeError(lookup_ids)
    i_lookup_full = np.arange(len(lookup_ids))
    i_sorted = np.searchsorted(all_ids, lookup_ids)
    i_lookup_bool = (i_sorted != len(all_ids))

    # filter out of range values
    i_lookup = i_lookup_full[i_lookup_bool]
    i_all = i_sorted[i_lookup_bool]

    # filter values in range but not found
    expected = lookup_ids[i_lookup]
    actual = all_ids[i_all]
    i_lookup2_bool = (expected == actual)
    if not np.any(i_lookup2_bool):
        null = np.array([], dtype='int32')
        return null, null
    i_lookup2 = i_lookup[i_lookup2_bool]
    i_all2 = i_all[i_lookup2_bool]

    # check
    expected2 = lookup_ids[i_lookup2]
    actual2 = all_ids[i_all2]
    assert np.array_equal(expected2, actual2)
    return i_lookup2, i_all2


#DOUBLE_CARDS = {'GRID', 'DMIG', 'DMIK', 'DMIJ', 'DMIJI', 'DMIAX'}
class VectorizedBaseCard:
    _skip_equality_check = True
    def __init__(self, model: BDF):
        self.model = model
        self.cards: list[tuple] = []
        self.comment: dict[int, str] = {}
        self.debug = False
        self.write_default_fields = True
        self.id = np.array([], dtype='int32')
        self.ifile = np.array([], dtype='int32')
        self.comment: dict[int, str] = {}
        #self.n = 0
        if hasattr(self, 'clear'):
            self.clear()
        #else:
            #self.n = 0
        assert self.n >= 0

    def __len__(self) -> int:
        return self.n

    def __contains__(self, key) -> bool:
        return key in self._ids

    def clear_check(func):
        @wraps(func)
        def wrapper(self):
            self.n = 0
            self.cards = []
            self.comment = {}
            return func(self)
        return wrapper

    def parse_cards_check(func):
        @wraps(func)
        def wrapper(self):
            if self.n == 0:
                return
            ncards = len(self.cards)
            if ncards == 0:
                return
            return func(self)
        return wrapper

    def _slice_comment(self, new_obj: VectorizedBaseCard, i: np.ndarray) -> None:
        """
        only slices comments when ids are unique
         - yes: GRID, Coord, CBAR, PSHELL, MAT1, RBE2
         - no:  MPC, SPC, SPC1, ASET, FREQ
        """
        assert isinstance(self.comment, dict), self.comment
        if not self.comment:
            return

        # new method
        #new_ids = self.node_id[i]
        try:
            # intersect the ids we're slicing with the existing ones
            if not i.max() <= self._ids.max():
                raise RuntimeError(f'_slice_comment: i={i}; {self._id_name}={self._ids}')
            new_ids = self._ids[i]
            common_ids = np.intersect1d(list(self.comment.keys()), new_ids)
            new_obj.comment = {nid: self.comment[nid] for nid in common_ids}
        except:
            print(self.comment)
            raise
        #----
        #old method
        #keys = list(self.comment.keys())
        #ikeys = [ii for ii in i
                 #if ii in keys]
        #new_obj.comment = {self.comment[i] for i in ikeys}


    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        self.model.log.warning(f'adding card {card}')
        self.cards.append((card, comment))
        i = len(self.cards) - 1
        self.n += 1
        return i

    @property
    def type(self) -> str:
        return self.__class__.__name__

    def rstrip(self) -> str:
        msg = self.write(size=8, write_card_header=False)
        return msg.rstrip()

    def __repr__(self) -> str:
        msg = f'{self.type}: n={self.n}'
        return msg

    def validate(self):
        self.model.log.warning(f'{self.type} has not implemented validate')

    @property
    def _ids(self) -> np.ndarray:
        return getattr(self, self._id_name)
    @_ids.setter
    def _ids(self, ids: np.ndarray) -> None:
        return setattr(self, self._id_name, ids)

    def remove_duplicates(self, inplace: bool=True) -> Any:
        """
        Removes duplicates

        Parameters
        ----------
        inplace: bool; default=True
            True: edit the current object
            False: return a new object

        Returns
        -------
        obj : parent class
            a class of the same type as the original

        """
        ids = self._ids
        uids = np.unique(ids)
        use_new_method = True
        if use_new_method:
            # get the sort index
            iarg = np.argsort(ids)

            # sort the sort index
            uarg = np.unique(iarg)

            # check if the sort index is sorted
            # True: no sorting is required
            # False: we need to sort
            if not np.array_equal(uarg, iarg):
                if len(iarg) == len(uarg):
                    # if the lengths are the same, we can use dumb sorting
                    jarg = iarg
                #else:
                    # we need to filter terms
                    #x = 1

                    # inplace
                    #self.__apply_slice__(card, i)
                # inplace
                self.__apply_slice__(self, jarg)
        else:
            # TODO: do we really need to sort?
            #       -> yes b/c we're going to call index
            self.sort()

        i = self.index(uids, inverse=False)
        if inplace:
            card = self
            self.__apply_slice__(card, i)
        else:
            card = self.slice_card_by_index(i, sort_ids=True)
        return card

    def slice_card_by_id(self, ids: np.ndarray,
                         assume_sorted: bool=True,
                         sort_ids: bool=False) -> Any:
        """uses a node_id to extract Elements, Properties, etc.

        Parameters
        ----------
        ids : (n,) int array
            the ids to extract (should be a subset of elem.element_id)
        assume_sorted: bool; default=True
            assume the parent array (e.g., elem.element_id is sorted)
        sort_ids: bool; default=False
            True: sort the input ids
            False: output is unsorted, which could cause issues

        """
        i = self.index(ids, assume_sorted=assume_sorted)
        cls_obj = self.slice_card_by_index(i, sort_index=sort_ids) # , assume_sorted=assume_sorted)
        #if 1:
            #assert np.array_equal(self._ids, np.unique(self._ids))
            #ids = np.atleast_1d(np.asarray(ids, dtype=self._ids.dtype))
        return cls_obj

    def remove_card_by_id(self, ids: np.ndarray) -> Any:
        """inplace operation"""
        index = self.index(ids, inverse=True)
        #if len(index) == 0:
        card = self
        self.__apply_slice__(card, index)
        return card
        #card = self.slice_card_by_index(index)
        #return card

    def slice_card_by_index(self, i: np.ndarray,
                            sort_index: bool=False,
                            #assume_sorted: bool=True,
                            ) -> Any:
        """
        Uses a node_index to extract Elements, Properties, etc.

        Parameters
        ----------
        i : (n,) int array
            the indices to extract (should be a subset of elem.element_id)
        sort_index: bool; default=False
            True: sort the input ids (i)
            False: output is unsorted, which could cause issues
        #assume_sorted: bool; default=True
            #assume the parent array (e.g., elem.element_id is sorted)

        """
        assert self.n > 0, self
        self_ids = self._ids
        assert len(self_ids), self_ids

        # i needs to be a numpy array of integers (vs. a list of integers)
        # in general, we'll use the dtype from:
        #  - self.node_id
        #  - self.element_id
        # but some cards are special snowflakes and don't have ids
        non_integer_based_ids = {
            'AECOMP', 'AECOMPL', 'MONPNT1', 'MONPNT2', 'MONPNT3',
            'MONDSP1', 'USET'}
        idtype = self_ids.dtype
        if idtype.name not in {'int32', 'int64'}:
            if self.type in non_integer_based_ids:
                idtype = self.model.idtype
            else:  # pragma: no cover
                raise RuntimeError(self.get_stats())
        i = np.atleast_1d(np.asarray(i, dtype=idtype))

        assert len(i) > 0, (i, self.type)
        assert i.min() >= 0, (i, self.type)
        if sort_index:
            i.sort()
            imax = i[-1]
        else:
            imax = i.max()
        imax_allowable = len(self_ids) - 1
        #if imax_allowable == -1:
            #imax_allowable = 0
        if imax > imax_allowable:
            raise RuntimeError(f'{self.type}: imax_allowable={imax_allowable}; ids={self_ids}; len(i)={imax}')
        cls = self.__class__
        card = cls(self.model)
        #card = CQUAD4(self.model)
        assert self.n > 0, self
        self.__apply_slice__(card, i)
        try:
            card.max_id
        except ValueError:
            raise ValueError(str(self) + self.write())
        assert card.n > 0, card
        return card

    def index(self, ids: np.ndarray,
              assume_sorted: bool=True,
              check_index: bool=True,
              inverse: bool=False) -> np.ndarray:
        """
        Parameters
        ----------
        ids: (n,) int array
            the node/element/property/material/etc. ids
        assume_sorted: bool; default=True
            assume the parent array (e.g., elem.element_id is sorted)
        check_index: bool; default=True
            validate the lookup
        inverse: bool; default=False
            False: get the indices for the ids
            True: get the inverse indices for the ids

        Returns
        -------
        index: (n,) int array
            the indicies in the node_ids array (or other array)

        Example
        -------
        >>> all_ids   = [1, 2, 3, 4, 5]
        >>> all_index = [0, 1, 2, 3, 4]
        >>> ids = [3, 4]
        >>> index(all_ids, ids, inverse=False)
        [2, 3]
        >>> index(all_ids, ids, inverse=True)
        [0, 1, 4]

        """
        if not assume_sorted:
            self.sort()
        self_ids = self._ids
        assert len(self_ids) > 0, f'{self.type}: {self._id_name}={self_ids}'
        if ids is None:
            return None # np.arange(len(ids), dtype='int32')
        ids = np.atleast_1d(np.asarray(ids, dtype=self_ids.dtype))
        ielem = np.searchsorted(self_ids, ids)
        if check_index:
            try:
                actual_ids = self_ids[ielem]
            except IndexError:
                raise
                missing_ids = np.setdiff1d(ids, self_ids)
                raise IndexError(f'{self.type}: missing {self._id_name}={missing_ids}')

            if not np.array_equal(actual_ids, ids):
                print(self.type)
                print(f'self.{self._id_name}= {self_ids}')
                print(f'{self._id_name}     = {ids}')
                missing = np.setdiff1d(self_ids, ids)
                assert len(missing) > 0, missing

                raise KeyError(f'{self.type}: expected_{self._id_name}={ids}\n'
                               f'actual_{self._id_name}={actual_ids}\n'
                               f'missing_{self._id_name}={missing}')
        if inverse:
            i = np.arange(len(self_ids), dtype=self_ids.dtype)
            index = np.setdiff1d(i, ielem)
            return index
        return ielem

    def sort(self) -> None:
        """sorts the card by node_id"""
        if hasattr(self, '_show_attributes'):
            sort_duplicates(self)
            self._is_sorted = True
            #if self.type == 'PSHELL':
                #assert len(self.property_id) == 1, self.property_id
            return

        self_ids = self._ids
        uid = np.unique(self_ids)
        if np.array_equal(uid, self_ids):
            return
        i = np.argsort(self_ids)
        assert len(i) > 0, self_ids
        self.__apply_slice__(self, i)

    @abstractmethod
    def __apply_slice__(self, prop: VectorizedBaseCard, i: np.ndarray) -> None:
        #...
        raise NotImplementedError(f'{self.type}: add __apply_slice__')

    #def write(self, size: int=8, is_double: bool=False,
              #write_card_header: bool=False) -> str:
        ##...
        #raise NotImplementedError(f'{self.type}: add write')

    def write(self, size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> str:
        """write typically exists and is probably overwritten; otherwise write_file is used"""
        stringio = StringIO()
        self.write_file(stringio, size=size, is_double=is_double, write_card_header=write_card_header)
        #if size == 8:
            #return self.write_8()
        #else:
            #return self.write_16(is_double=is_double)
        msg = stringio.getvalue()
        return msg

    def write_8(self, is_double: bool=False,
                write_card_header: bool=False) -> str:
        """write is the base function"""
        #stringio = StringIO()
        msg = self.write(size=8, is_double=is_double, write_card_header=write_card_header)
        #if hasattr(self, 'write_file_8'):
           #self.write_file_8(stringio, is_double=is_double,
                             #write_card_header=write_card_header)
        #else:
           #self.write_file(stringio, size=8, is_double=is_double,
                           #write_card_header=write_card_header)
        #msg = stringio.getvalue()
        return msg

    def write_16(self, is_double: bool=False,
                 write_card_header: bool=False) -> str:
        """write is the base function"""
        #stringio = StringIO()
        msg = self.write(size=16, is_double=is_double, write_card_header=write_card_header)
        #if hasattr(self, 'write_file_16'):
        #    self.write_file_16(stringio, is_double=is_double, write_card_header=write_card_header)
        #else:
        #    self.write_file(stringio, size=8, is_double=is_double, write_card_header=write_card_header)
        #msg = stringio.getvalue()
        return msg

    #def write_file_8(self, file_obj: TextIOLike,
    #                  is_double: bool=False,
    #                  write_card_header: bool=False) -> None:
    #    if not hasattr(self, 'write_file'):
    #        raise NotImplementedError(f'{self.type}: add write_file')
    #    self.write_file(file_obj, size=8, write_card_header=write_card_header)

    #def write_file_16(self, bdf_file: TextIOLike,
    #                  is_double: bool=False,
    #                  write_card_header: bool=False) -> None:
    #    if not hasattr(self, 'write_file'):
    #        raise NotImplementedError(f'{self.type}: add write_file')
        #self.write_file(bdf_file, size=16, is_double=is_double, write_card_header=write_card_header)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if not hasattr(self, 'write_file_8'):
            raise NotImplementedError(f'{self.type}: add write_file')

        if size == 8:
            self.write_file_8(bdf_file, write_card_header=write_card_header)
        else:
            self.write_file_16(bdf_file, is_double=is_double,
                               write_card_header=write_card_header)
        return

    def get_attributes(self, mode: str='public',
                       keys_to_skip: Optional[list[str]]=None,
                       filter_properties: bool=False) -> list[str]:
        if hasattr(self, '_show_attributes'):
            return self._show_attributes

        keys_to_skip = keys_to_skip if keys_to_skip is not None else []
        keys_to_skip.append('model')
        #if hasattr(self, '_remove_attributes'):
            #keys_to_skip += self._remove_attributes
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip,
                                 filter_properties=filter_properties)
    def get_stats(self,
                  mode: str='public',
                  keys_to_skip: Optional[list[str]]=None,
                  filter_properties: bool=False) -> str:
        keys_to_skip = keys_to_skip if keys_to_skip is not None else []
        keys_to_skip.append('model')
        return object_stats(self, mode=mode,
                            keys_to_skip=keys_to_skip,
                            filter_properties=filter_properties)


class Element(VectorizedBaseCard):
    _id_name = 'element_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n: int = 0
        self.element_id: np.ndarray = np.array([], dtype='int32')
        if hasattr(self, 'clear'):
            self.clear()

    def slice_card_by_element_id(self, element_id: np.ndarray,
                                 sort_ids: bool=False) -> Element:
        assert self.n > 0, self.n
        assert len(self.element_id) > 0, self.element_id
        i = self.index(element_id)
        #cls_obj = cls(self.model)
        #cls_obj.__apply_slice__(self, i)
        cls_obj = self.slice_card_by_index(i, sort_index=sort_ids)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        return 0
        #element_id_all = used_dict['element_id']
        #element_id = np.intersect1d(element_id_all, self.element_id)

        #ncards_removed = len(self.element_id) - len(element_id)
        #if ncards_removed:
            #if len(element_id) == 0:
                #self.clear()
            #else:
                #try:
                    #self.slice_card_by_id(element_id, assume_sorted=True, sort_ids=False)
                #except IndexError:
                    #raise RuntimeError(self.get_stats())
                #except ValueError:
                    #raise RuntimeError(f'{self.type} element_id is empty...n={self.n}')
        #return ncards_removed

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2
        if hasattr(self, 'g0'):
            nodes = self.g0
            for i, nid1 in enumerate(nodes):
                if nid1 == 0:
                    continue
                nid2 = nid_old_to_new.get(nid1, nid1)
                nodes[i] = nid2

    def filter_by_ifile_filter(self,
                               ifile_filter: np.ndarray) -> tuple[bool, Any]:
        log = self.model.log
        is_removed = True
        element = self
        nelement = element.n
        if len(ifile_filter):
            if hasattr(element, 'ifile'):
                ifile = element.ifile
                if len(element.ifile) == 0:
                    log.error(f'{element.type} has no ifile; skipping')
                    is_removed = False
                    return is_removed, element
            else:
                log.error(f'{element.type} has no ifile; skipping')
                is_removed = False
                return is_removed, element
                #log.error(f'{element.type} has no ifile; assuming ifile=0')
                #ifile = np.zeros(nelement, dtype='int32')

            is_in = np.isin(ifile, ifile_filter)
            if not np.any(is_in):
                log.error(f'skipping {element.type} because it is empty after filtering; '
                          f'ifile={np.unique(element.ifile)}')
                return is_removed, element
            if not np.all(is_in):
                #print(f'is_in = {is_in.tolist()}')
                assert len(is_in) == nelement
                # assert False not in is_in, is_in
                index = np.zeros(len(is_in), dtype='int32')[is_in]
                element2 = element.slice_card_by_index(index)
                assert element2.n <= element.n, element2.n
                element = element2
                # assert element.n == element.n
                # nelement = element2.n
                # assert element == element2
                del element2, index

        elif not hasattr(element, 'ifile'):
            log.error(f'{element.type} has no ifile')
        is_removed = False
        return is_removed, element

    def filter(self,
               ifile_filter: np.ndarray,
               property_filter: np.ndarray,
               material_filter: np.ndarray) -> tuple[bool, Any]:
        """filers """
        log = self.model.log
        element = self
        if len(ifile_filter) == 0 or len(property_filter) == 0 or len(material_filter) == 0:
            is_removed = False
            return is_removed, element

        is_removed = True
        nelement = element.n

        is_in_file = np.isin(element.ifile, ifile_filter)
        is_ins = [is_in_file]
        if hasattr(element, 'property_id'):
            log.debug(f'property_id')
            property_id = element.property_id
            #upids = np.unique(property_id)
            is_in_prop = np.isin(element.property_id, property_filter)
            if not np.all(is_in_prop):
                log.debug(f'  all')
            elif not np.all(is_in_prop):
                log.debug(f'  adding')
                is_ins.append(is_in_prop)
            else:
                log.debug(f'  else...n={is_in_prop.sum()}/{self.n}')


        if hasattr(element, 'material_id'):
            log.debug(f'material_id')
            is_in_mat = np.isin(element.material_id, material_filter)
            if not np.all(is_in_mat):
                log.debug(f'  all')
            elif not np.all(is_in_mat):
                log.debug(f'  adding')
                is_ins.append(is_in_mat)
            else:
                log.debug(f'  else...n={is_in_mat.sum()}/{self.n}')
        else:
            if -1 not in material_filter:  # PSHELLs
                material_filter_ = np.hstack([-1], material_filter)
            #is_in_mat = get_material_ids_from_property(
                #self.property_id, self.allowed_materials,
                #material_filter)
            is_in_mat = np.any(is_in_mats, axis=1)
            assert len(is_in_mat) == nelement, len(is_in_mat)
            if not np.all(is_in_mat):
                log.debug(f'  all')
            elif not np.all(is_in_mat):
                log.debug(f'  adding')
                is_ins.append(is_in_mat)
            else:
                log.debug(f'  else...n={is_in_mat.sum()}/{self.n}')

        log.debug(f'n is_ins={len(is_ins)}')
        if len(is_ins) == 0:
            is_removed = True
            return is_removed, element
        elif len(is_ins) == 1:
            is_in = is_ins[0]
        elif len(is_ins) == 2:
            is_in = is_ins[0] & is_ins[1]
        elif len(is_ins) == 3:
            is_in = is_ins[0] & is_ins[1] & is_ins[2]
        else:  # pragma: no cover
            #is_in = np.logical_and(*is_ins)
            raise RuntimeError(len(is_ins))

        assert len(is_in) == nelement, len(is_in)
        #print(f'prop is_in = {is_in.tolist()}')
        if not np.any(is_in):
            log.error(f'skipping {element.type} because it is empty after filtering; '
                      f'property_ids={upids}')
            return is_removed, element

        is_removed = False
        if np.all(is_in):
            # keep all elements
            return is_removed, element

        nexpected = is_in.sum()
        #print(f'is_in = {is_in.tolist()}')
        assert len(is_in) == nelement
        # assert False not in is_in, is_in
        for pid, flag in zip(element.property_id, is_in):
            is_in2 = pid in property_filter
            assert flag == is_in2, (pid, flag, is_in2)
        index = np.arange(len(is_in), dtype='int32')[is_in]
        #print('pid', element.property_id[is_in])
        #print('upid', np.unique(element.property_id[is_in]))
        #print('slicing...')
        element2 = element.slice_card_by_index(index)
        #print('back from slice')
        assert element2.n <= element.n, element2.n
        #print(f'check n; {element.type}; nexpected={nexpected} n={element2.n}')
        #print('setting...')
        element = element2
        # assert element.n == element.n
        # nelement = element2.n
        # assert element == element2
        # del element2, index
        # print(element.write())
        # assert element.element_id[0] != element.element_id[1]
        return is_removed, element

    def filter_by_property_filter(self,
                                  property_filter: np.ndarray) -> tuple[bool, Any]:
        """property filter is inclusive"""
        # property_filter = np.arange(20000)
        log = self.model.log
        element = self
        if not hasattr(element, 'property_id') or len(property_filter) == 0:
            is_removed = False
            return is_removed, element

        is_removed = True
        nelement = element.n

        property_id = element.property_id
        upids = np.unique(property_id)

        is_in = np.isin(property_id, property_filter)
        #print(f'prop is_in = {is_in.tolist()}')
        if not np.any(is_in):
            log.error(f'skipping {element.type} because it is empty after filtering; '
                      f'property_ids={upids}')
            return is_removed, element

        is_removed = False
        if np.all(is_in):
            # keep all elements
            return is_removed, element

        nexpected = is_in.sum()
        #print(f'is_in = {is_in.tolist()}')
        assert len(is_in) == nelement
        # assert False not in is_in, is_in
        for pid, flag in zip(element.property_id, is_in):
            is_in2 = pid in property_filter
            assert flag == is_in2, (pid, flag, is_in2)
        index = np.arange(len(is_in), dtype='int32')[is_in]
        #print('pid', element.property_id[is_in])
        #print('upid', np.unique(element.property_id[is_in]))
        #print('slicing...')
        element2 = element.slice_card_by_index(index)
        #print('back from slice')
        assert element2.n <= element.n, element2.n
        #print(f'check n; {element.type}; nexpected={nexpected} n={element2.n}')
        #print('setting...')
        element = element2
        # assert element.n == element.n
        # nelement = element2.n
        # assert element == element2
        # del element2, index
        # print(element.write())
        # assert element.element_id[0] != element.element_id[1]
        return is_removed, element


class Property(VectorizedBaseCard):
    _id_name = 'property_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n = 0
        self.property_id: np.ndarray = np.array([], dtype='int32')

    def slice_card_by_property_id(self, property_id: np.ndarray,
                                 sort_ids: bool=False) -> Property:
        assert self.n > 0, self.n
        assert len(self.property_id) > 0, self.property_id
        i = self.index(property_id)
        cls_obj = self.slice_card_by_index(i, sort_index=sort_ids)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        property_id_all = used_dict['property_id']
        property_id = np.intersect1d(property_id_all, self.property_id)
        removed_id = np.setdiff1d(self.property_id, property_id_all)
        ncards_removed = len(self.property_id) - len(property_id)
        if ncards_removed:
            if len(property_id) == 0:
                self.model.log.info(f'removing {self.type} removed property_id={removed_id}; ncards={ncards_removed}')
                self.clear()
            else:
                try:
                    self.slice_card_by_id(property_id, assume_sorted=True, sort_ids=False)
                except IndexError:
                    print(self.get_stats())
                    raise
                except ValueError:
                    raise RuntimeError(f'{self.type} property_id is empty...n={self.n}')
        return ncards_removed


class Material(VectorizedBaseCard):
    _id_name = 'material_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n = 0
        self.material_id: np.ndarray = np.array([], dtype='int32')
        if hasattr(self, 'clear'):
            self.clear()

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        pass

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        material_id_all = used_dict['material_id']
        material_id = np.intersect1d(material_id_all, self.material_id)
        removed_id = np.setdiff1d(self.material_id, material_id_all)
        ncards_removed = len(self.material_id) - len(material_id)
        if ncards_removed:
            if len(material_id) == 0:
                self.model.log.info(f'removing {self.type} removed material_id={removed_id}; ncards={ncards_removed}')
                self.clear()
            else:
                try:
                    self.slice_card_by_id(material_id, assume_sorted=True, sort_ids=False)
                except IndexError:
                    raise RuntimeError(self.get_stats())
                except ValueError:
                    raise RuntimeError(f'{self.type} material_id is empty...n={self.n}')
        return ncards_removed

    def slice_card_by_material_id(self, material_id: np.ndarray,
                                  sort_ids: bool=False) -> Material:
        assert self.n > 0, self.n
        i = self.index(material_id)
        cls_obj = self.slice_card_by_index(i, sort_index=sort_ids)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def get_density(self) -> np.ndarray:
        return self.rho

    @abstractmethod
    def __apply_slice__(self, mat, i: np.ndarray) -> None:  # ignore[override]
        #...
        raise NotImplementedError(f'{self.type}: add __apply_slice__')


def parse_check(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        ids = getattr(self, self._id_name)
        if len(ids) == 0:
            if self.n == 0:
                return
            self.parse_cards()
        return func(self, *args, **kwargs)
    return wrapper

def get_print_card_8_16(size: int) -> Callable[[int], list]:
    if size == 8:
        print_card = print_card_8
    else:
        print_card = print_card_16
    return print_card


def remove_unused_primary(card, used_ids: np.ndarray, ids: np.ndarray,
                          ids_str: str) -> int:
    """
    spc_ids = used_dict['spc_id']
    ncards_removed = remove_unused_primary(self, spc_ids, self.spc_id, 'spc_id')
    """
    assert len(ids) == len(np.unique(ids)), f'{card.type}: did you mean to use remove_unused_duplicate?'
    spc_id = np.intersect1d(used_ids, ids)
    removed_id = np.setdiff1d(ids, used_ids)
    ncards_removed = len(ids) - len(spc_id)
    if ncards_removed:
        if len(spc_id) == 0:
            card.model.log.info(f'removing {card.type} removed spc_id={removed_id}; ncards={ncards_removed}')
            card.clear()
        else:
            try:
                card.slice_card_by_id(spc_id, assume_sorted=True, sort_ids=False)
            except IndexError:
                raise RuntimeError(card.get_stats())
            except ValueError:
                raise RuntimeError(f'{card.type} {ids_str} is empty...n={card.n}')
    return ncards_removed

def remove_unused_duplicate(card, used_ids: np.ndarray, ids: np.ndarray,
                            ids_str: str) -> int:
    """FORCE, FORCE1, LOAD, ..."""
    index = np.array([i for i, load_id in enumerate(ids)])
    if len(index):
        card.__apply_slice__(card, index)
        ncards_removed = len(index)
    else:
        card.clear()
        ncards_removed = 0
    return ncards_removed

# -----------------------------------------------------------

def is_equal_by_index(card: VectorizedBaseCard,
                      inode1: np.ndarray,
                      inode2: np.ndarray) -> bool:
    if card._skip_equality_check:
        return False
    is_equal = True
    for key in card._show_attributes:
        values = getattr(card, key)
        val1 = values[inode1]
        val2 = values[inode2]
        if not np.array_equal(val1, val2):
            print(f'key={key!r} are not equal')
            is_equal = False
            break
    return is_equal


def compare_by_indexs(card: VectorizedBaseCard,
                      iarg1: np.ndarray,
                      iarg2: np.ndarray) -> tuple[str]:
    attributes = card.get_attributes(mode='public', keys_to_skip=None,
                                     filter_properties=False)
    name = card.type.lower()
    card1 = card.slice_card_by_index(iarg1)
    card2 = card.slice_card_by_index(iarg2)
    msg = ''
    for value_name in attributes:
        value1 = getattr(card1, value_name)
        value2 = getattr(card2, value_name)
        if value1.dtype.name in {'int32', 'int64', 'float32', 'float64'}:
            dvalue = value1 - value2
            if np.abs(dvalue).max() > 0:
                msg += f'd{value_name}:\n{dvalue}\n'
        else:
            # strings
            valuei = value1[0]
            assert isinstance(valuei, str), valuei
            is_diff = (value1 != value2)
            if np.any(is_diff):
                msg += (
                    f'{value_name}:\n'
                    f'  {value_name}1: {value1[is_diff]}\n'
                    f'  {value_name}2: {value2[is_diff]}\n')

    if msg:
        msg = (f'arg1={iarg1}\n'
               f'arg2={iarg2}\n'
               f'{name}1:\n{card1.write()}'
               f'{name}2:\n{card2.write()}') + msg
    return msg
    #dxyz = node1.xyz - node2.xyz
    #dcp = node1.cp - node2.cp
    #dcd = node1.cd - node2.cd
    #if np.abs(dxyz).max() > 0.0:
        #msg += f'dxyz:\n{dxyz}\n'
    #if np.abs(dcp).max() > 0:
        #msg += f'dcp={dcp}\n'
    #if np.abs(dcd).max() > 0:
        #msg += f'dcd={dcd}\n'
    #raise RuntimeError(msg)


def sort_duplicates(card: VectorizedBaseCard) -> None:
    """
    sort_duplicates is sort, but you can have duplicate values
    the duplicate values will be removed only if the values are unique
    (so all GRID NID=1 are identical).  Otherwise, there will be a crash
    """
    if not hasattr(card, '_id_name'):
        raise NotImplementedError('_id_name is required for sort_duplicates')
    if not hasattr(card, '_show_attributes'):
        raise NotImplementedError('_show_attributes is required for sort_duplicates')
    # get the sort index
    #ids = self.node_id
    ids = getattr(card, card._id_name)
    iarg = np.argsort(ids)
    rarg = np.flip(iarg)

    # sort the sort index
    uarg, idx, rdx, ucounts = np.unique(
        ids, return_index=True, return_inverse=True,
        return_counts=True)

    ids_sorted = ids[iarg]
    uids_sorted = np.unique(ids_sorted)
    if not np.array_equal(uarg, iarg) or len(ids) != len(uids_sorted):
        #print(f'card.type = {card.type}')
        model = card.model
        is_duplicate_id = ucounts.max() > 1
        if is_duplicate_id:
            if card.type in model.allow_overwrites_set:
                isort = np.searchsorted(ids_sorted, uarg, side='right') - 1
                #model.log.warning('picking lower indicies (first added to model)')
                card.__apply_slice__(card, isort)
                card._is_sorted = True
                return
            else:
                # we have duplicate/unsorted nodes
                dnode = ids_sorted[1:] - ids_sorted[:-1]  # high - low
                # dnode = np.diff(ids_sorted)
                # inode = (dnode == 0)
                inode1 = np.hstack([False, dnode == 0])
                inode2 = np.hstack([dnode == 0, False])
                # ids_unique1 = ids_sorted[inode1] # high node
                # ids_unique2 = ids_sorted[inode2] # low  node

                iarg1 = iarg[inode1]
                iarg2 = iarg[inode2]
                if not is_equal_by_index(card, iarg1, iarg2):
                    msg = compare_by_indexs(card, iarg1, iarg2)
                    if msg != '' or not card._skip_equality_check:
                        raise RuntimeError(msg)

                uids_sorted, idx = np.unique(ids, return_index=True)
                nunique_nodes = len(uids_sorted)
                assert nunique_nodes == len(idx)

                card.__apply_slice__(card, idx)
                card._is_sorted = True
                return
                # raise RuntimeError(f'ids={ids}')
        else:
            card.__apply_slice__(card, idx)
            card._is_sorted = True
            return

        # we have duplicate/unsorted nodes
        # print(f'{card._id_name} = {ids} -> ids[iarg]={ids[iarg]}')
        # print(f'idx = {idx}; iarg={iarg}; uarg={uarg}')
        # print(f'rdx = {rdx}; rarg={rarg}')
        # print(f'ucounts = {ucounts}')
        # mapper = {i: nid for i, nid in enumerate(ids.tolist())}
        # print(f'mapper = {mapper}')
        #
        #print(iarg.tolist())
        #print('->  ', ids_sorted.tolist())
        #print('--> ', uids_sorted.tolist())
        assert (iarg - uarg).sum() == 0, (iarg - uarg).sum()
        #if len(iarg) == len(uarg):
            ## if the lengths are the same, we can use dumb sorting

            ## check on dumb sort; no missing entries and no duplicates
            #assert (iarg - uarg).sum() == 0, (iarg - uarg).sum()

        is_duplicate_ids = (len(ids_sorted) != len(uids_sorted))
        # print(f'ids = {ids}')
        # print(f'ids_sorted = {ids_sorted}')
        # print(f'is_duplicate_ids = {is_duplicate_ids}')
        model = card.model
        if is_duplicate_ids and card.type in model.allow_overwrites_set:
            #print(f'ids_sorted = {ids_sorted}')
            #print(f'uids_sorted = {uids_sorted}')
            idx = np.searchsorted(ids_sorted, uids_sorted)
            #idx = np.searchsorted(ids_sorted, uids_sorted, side='right') - 1
            # print(f'idx = {idx}')
            assert np.array_equal(ids_sorted[idx], uids_sorted)
            card.__apply_slice__(card, idx)
            # card._is_sorted = True
            # return

        if is_duplicate_ids:
            # Now that we've sorted the ids, we'll take the delta id with the next id.
            # We check the two neighoring values when there is a duplicate,
            # so if they're they same in xyz, cp, cd, etc., the results are the same
            # we don't need to check unique ids
            dnode = ids_sorted[1:] - ids_sorted[:-1] # high - low
            #dnode = np.diff(ids_sorted)
            #inode = (dnode == 0)
            inode1 = np.hstack([False, dnode == 0])
            inode2 = np.hstack([dnode == 0, False])
            #ids_unique1 = ids_sorted[inode1] # high node
            #ids_unique2 = ids_sorted[inode2] # low  node

            iarg1 = iarg[inode1]
            iarg2 = iarg[inode2]
            if not is_equal_by_index(card, iarg1, iarg2):
                msg = compare_by_indexs(card, iarg1, iarg2)
                if msg != '' or not card._skip_equality_check:
                    raise RuntimeError(msg)

            # now that we know all the nodes are unique
            # (i.e., 2x node_id=1 is OK, but they have the same xyz),
            # we can slice ids to get the unique ids. Then take
            # the mapping (idx) and slice the data to make the nodes
            # unique
            uids_sorted, idx = np.unique(ids, return_index=True)
            nunique_nodes = len(uids_sorted)
            assert nunique_nodes == len(idx)

            card.__apply_slice__(card, idx)

            # check we found the unique ids
            node_id = getattr(card, card._id_name)
            assert np.array_equal(node_id, uids_sorted)
            #assert np.array_equal(self.node_id, np.unique(self.node_id))
        else:
            # no duplicate nodes
            #
            # we should still sort using iargsort
            card.__apply_slice__(card, iarg)

    ids2 = getattr(card, card._id_name)
    unid = np.unique(ids2)

    if len(ids2) != len(unid):
        # we need to filter nodes
        msg = f'len(self.{ids2})={len(ids2)} != len(unid)={len(unid)}'
        print(msg)

    #self.remove_duplicates(inplace=False)
    card._is_sorted = True

def save_ifile_comment(self, ifile, comment):
    assert len(ifile) > 0, ifile
    if comment is not None and len(comment):
        #print(f'comment = {comment!r}')
        #print(f'self.comment = {self.comment!r}')
        self.comment.update(comment)
    self.ifile = ifile
