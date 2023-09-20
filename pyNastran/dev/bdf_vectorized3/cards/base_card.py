from __future__ import annotations
from abc import abstractmethod
from functools import wraps
from io import StringIO
from typing import Callable, Optional, Any, TYPE_CHECKING
import numpy as np
from pyNastran.dev.bdf_vectorized3.utils import print_card_8 # , print_float_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
#from pyNastran.bdf.bdf_interface.assign_type import (
    #integer, integer_or_blank, double_or_blank, components_or_blank)
#from pyNastran.dev.bdf_vectorized3.cards.line_elements import CBAR, CBEAM, PBAR, PBARL, PBEAM, PBEAML
#from pyNastran.dev.bdf_vectorized3.cards.shell_elements import CQUAD4, CTRIA3, PSHELL, PCOMP, PCOMPG
#from pyNastran.dev.bdf_vectorized3.cards.solid_elements import CTETRA, CHEXA, CPENTA, CPYRAM
#from pyNastran.dev.bdf_vectorized3.cards.static_loads import PLOAD4, FORCE
from pyNastran.utils import object_stats # object_attributes, object_methods,
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
    #if debug:
        #x = 2
    return i_lookup2, i_all2


#DOUBLE_CARDS = {'GRID', 'DMIG', 'DMIK', 'DMIJ', 'DMIJI', 'DMIAX'}
class VectorizedBaseCard:
    def __init__(self, model: BDF):
        self.model = model
        self.cards: list[tuple] = []
        self.comment: dict[int, str] = {}
        self.n = 0
        self.debug = False
        self.write_default_fields = True
        self.id = np.array([], dtype='int32')

    def __len__(self) -> int:
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
            card = self.slice_card_by_index(i)
        return card

    def slice_card_by_id(self, ids: np.ndarray,
                         assume_sorted: bool=True) -> Any:
        """uses a node_id to extract Elements, Properties, etc.

        Parameters
        ----------
        ids : (n,) int array
            the ids to extract (should be a subset of elem.element_id)
        assume_sorted: bool; default=True
            assume the parent array (e.g., elem.element_id is sorted)

        """
        i = self.index(ids, assume_sorted=assume_sorted)
        cls_obj = self.slice_card_by_index(i) # , assume_sorted=assume_sorted)
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
                            #assume_sorted: bool=True,
                            ) -> Any:
        """
        Uses a node_index to extract Elements, Properties, etc.

        Parameters
        ----------
        i : (n,) int array
            the indices to extract (should be a subset of elem.element_id)
        assume_sorted: bool; default=True
            assume the parent array (e.g., elem.element_id is sorted)

        """
        assert self.n > 0, self
        self_ids = self._ids
        assert len(self_ids), self_ids
        i = np.atleast_1d(np.asarray(i, dtype=self_ids.dtype))
        i.sort()
        imax = i[-1]
        imax_allowable = len(self_ids) - 1
        #if imax_allowable == -1:
            #imax_allowable = 0
        if imax > imax_allowable:
            raise RuntimeError(f'imax_allowable={imax_allowable}; ids={self_ids}; len(i)={imax}')
        cls = self.__class__
        card = cls(self.model)
        #card = CQUAD4(self.model)
        self.__apply_slice__(card, i)
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
            actual_ids = self_ids[ielem]
            if not np.array_equal(actual_ids, ids):
                raise KeyError(f'{self.type}: expected_{self._id_name}={ids}; actual_{self._id_name}={actual_ids}')
        if inverse:
            i = np.arange(len(self_ids), dtype=self_ids.dtype)
            index = np.setdiff1d(i, ielem)
            return index
        return ielem

    def sort(self) -> None:
        self_ids = self._ids
        uid = np.unique(self_ids)
        if np.array_equal(uid, self_ids):
            return
        i = np.argsort(self_ids)
        self.__apply_slice__(self, i)

    @abstractmethod
    def __apply_slice__(self, prop: VectorizedBaseCard, i: np.ndarray) -> None:
        #...
        raise NotImplementedError(f'{self.type}: add __apply_slice__')

    #def write(self, size: int=8, is_double: bool=False, write_card_header: bool=False) -> str:
        ##...
        #raise NotImplementedError(f'{self.type}: add write')

    def write(self, size: int=8, is_double: bool=False, write_card_header: bool=False) -> str:
        """write typically exists and is probably overwritten; otherwise write_file is used"""
        stringio = StringIO()
        self.write_file(stringio, size=size, is_double=is_double, write_card_header=write_card_header)
        #if size == 8:
            #return self.write_8()
        #else:
            #return self.write_16(is_double=is_double)
        msg = stringio.getvalue()
        return msg

    def write_8(self, is_double: bool=False, write_card_header: bool=False) -> str:
        """write is the base function"""
        #stringio = StringIO()
        msg = self.write(size=8, is_double=is_double, write_card_header=write_card_header)
        #if hasattr(self, 'write_file_8'):
        #    self.write_file_8(stringio, is_double=is_double, write_card_header=write_card_header)
        #else:
        #    self.write_file(stringio, size=8, is_double=is_double, write_card_header=write_card_header)
        #msg = stringio.getvalue()
        return msg

    def write_16(self, is_double: bool=False, write_card_header: bool=False) -> str:
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
            self.write_file_16(bdf_file, is_double=is_double, write_card_header=write_card_header)
        return

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

    def slice_card_by_element_id(self, element_id: np.ndarray) -> Element:
        assert self.n > 0, self.n
        assert len(self.element_id) > 0, self.element_id
        i = self.index(element_id)
        #cls_obj = cls(self.model)
        #cls_obj.__apply_slice__(self, i)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj


class Property(VectorizedBaseCard):
    _id_name = 'property_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n = 0
        self.property_id: np.ndarray = np.array([], dtype='int32')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> Property:
        assert self.n > 0, self.n
        assert len(self.property_id) > 0, self.property_id
        i = self.index(property_id)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj


class Material(VectorizedBaseCard):
    _id_name = 'material_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards = []
        self.n = 0
        self.material_id: np.ndarray = np.array([], dtype='int32')

    def slice_card_by_material_id(self, material_id: np.ndarray) -> Material:
        assert self.n > 0, self.n
        i = self.index(material_id)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def get_density(self) -> np.ndarray:
        return self.rho

    @abstractmethod
    def __apply_slice__(self, mat, i: np.ndarray) -> None:  # ignore[override]
        #...
        raise NotImplementedError(f'{self.type}: add __apply_slice__')


def parse_element_check(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if len(self.element_id) == 0:
            if self.n == 0:
                return
            self.parse_cards()
        return func(self, *args, **kwargs)
    return wrapper

def parse_property_check(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if len(self.property_id) == 0:
            if self.n == 0:
                return
            self.parse_cards()
        return func(self, *args, **kwargs)
    return wrapper

def parse_material_check(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if len(self.material_id) == 0:
            if self.n == 0:
                return
            self.parse_cards()
        return func(self, *args, **kwargs)
    return wrapper

def get_print_card_8_16(size: int) -> Callable[list]:
    if size == 8:
        print_card = print_card_8
    else:
        print_card = print_card_16
    return print_card
