from io import StringIO
from numpy import (array, searchsorted, array_equal, setdiff1d, int64, argsort,
                   arange, ndarray, asarray)
from pyNastran.utils import object_attributes, object_methods
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import _format_comment


class BaseMethods:
    def object_attributes(self, mode='public', keys_to_skip=None,
                          filter_properties=False):
        """.. seealso:: `pyNastran.utils.object_attributes(...)`"""
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode='public', keys_to_skip=None):
        """.. seealso:: `pyNastran.utils.object_methods(...)`"""
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)


class VectorizedCard:
    type = 'VectorizedCard'
    def __init__(self, model):
        self.model = model
        self.n = 0
        self.i = 0
        self._comments = {}
        if 0:
            if self.type in model._element_name_to_element_type_mapper:
                self.op2_id = model._element_name_to_element_type_mapper[self.type]
            else:
                if self.type.startswith('C'):
                    self.model.log.warning('there is no op2_id to apply for element=%r' % self.type)

    def set_comment(self, idi, new_comment):
        """sets a comment"""
        self._comments[idi] = _format_comment(new_comment)

    def __len__(self):
        return self.n

    def shrink(self, refcheck=True):
        raise NotImplementedError()

    def resize(self, n, refcheck=True):
        names = object_attributes(self, mode="public")
        for name in names:
            attr = getattr(self, name)
            if isinstance(attr, ndarray):
                #self.model.log.info('resizing %r; shape=%s; size=%s' % (
                    #name, attr.shape, attr.size))
                # resize the array
                shape2 = list(attr.shape)
                shape2[0] = n
                attr.resize(tuple(shape2), refcheck=refcheck)

                if n > self.n:
                    # TODO: fill the data with nan values ideally, but it's not working
                    if attr.ndim == 1:
                        attr[self.n:] = 0
                    elif attr.ndim == 2:
                        attr[self.n:, :] = 0
                    elif attr.ndim == 3:
                        attr[self.n:, :, :] = 0
                    else:
                        raise NotImplementedError(attr.shape)
                    #print(attr)
            else:
                # metadata
                pass
        if self.i >= n:
            self.i = n
        self.n = n

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % (self.type, self.n))
        return msg

    def print_card(self, i=None, size=8):
        #i = self._validate_slice(i)
        string_io = StringIO()
        self.write_card(string_io, i, size=size)
        return string_io.getvalue().rstrip()

    def _verify(self, xref=True):
        raise NotImplementedError(self.type)

    def _validate_slice(self, i):
        if self.n == 0:
            raise RuntimeError('%s has not been allocated or built' % self.type)
        if isinstance(i, (int, int64)):
            i2 = array([i], dtype='int32')
        elif isinstance(i, list):
            i2 = asarray(i)
        elif i is None:
            i2 = None  # slice(None)
        elif len(i.shape) == 1:
            i2 = i
        else:
            #print(i, type(i), i.shape)
            i2 = i
        #print('shape=%s' % str(i2.shape))
        return i2

    def _set_as_array(self, i):
        if isinstance(i, (int, int64)):
            i = array([i], dtype='int32')
        elif isinstance(i, list):
            i = array(i, dtype='int32')
        elif i is None:
            raise RuntimeError(i)
        elif len(i.shape) == 1:
            pass
        else:
            self.model.log.warning('???', i, type(i))
        return i

    def _get_sorted_index(self, sorted_array, unsorted_array, field_name, msg, check=True):
        if field_name == 'coord_id':
            if sorted_array.min() < 0:
                msg = '%s.min()=%s; did you allocate & build it?\n' % (
                    field_name, sorted_array.min())
                raise RuntimeError(msg)
        else:
            if sorted_array.min() < 1:
                msg = '%s.min()=%s; did you allocate & build it?\n' % (
                    field_name, sorted_array.min())
                raise RuntimeError(msg)

        if not array_equal(argsort(sorted_array), arange(len(sorted_array))):
            msg2 = '%s is not sorted\nsorted_array=%s' % (msg, sorted_array)
            raise RuntimeError(msg2)
        if not isinstance(self.n, (int64, int)):
            raise TypeError('field_name=%s n=%s type=%s' % (field_name, self.n, type(self.n)))
        if not isinstance(check, bool):
            raise TypeError('check should be a bool; check=%r' % check)

        if unsorted_array is None:
            i = slice(None)
            if self.n == 1:
                i = array([0], dtype='int32')
            return i
        else:
            i = searchsorted(sorted_array, unsorted_array)
            i.astype('int32')
            if check:
                try:
                    sorted_array_i = sorted_array[i]
                except IndexError:
                    #msg = 'sorted_array_i = sorted_array[i] - %s\n' % field_name
                    msg = 'sorted_%s=%s\n' % (field_name, sorted_array)
                    msg += 'i=%s' % (i)
                    raise IndexError(msg)

                if not array_equal(sorted_array_i, unsorted_array):
                    # undefined nodes/elements
                    #print('unsorted %s' % unsorted_array)
                    #print('sorted %s' % sorted_array)
                    #pass
                    msg2 = 'Undefined %s\n' % msg
                    msg2 += 'diff=%s\n' % setdiff1d(unsorted_array, sorted_array)
                    msg2 += 'sorted %s= %s; n=%s\n' % (field_name, str(sorted_array), len(sorted_array))
                    msg2 += 'unsorted %s = %s\n' % (field_name, unsorted_array)
                    msg2 += 'sorted %s[i]= %s\n' % (field_name, sorted_array[i])
                    msg2 += 'i=%s\n' % i
                    raise RuntimeError(msg2)
        if isinstance(i, (int64, integer_types)):
            i = array([i], dtype='int32')
        return i

    #def __repr__(self):
        #f = StringIO()
        #self.write_card(f)
        #return f.getvalue().rstrip()

    def __repr__(self):
        return '<%s object; n=%s>' % (self.type, self.n)

def by_converter(value, default):
    """
    For use in:
        - get_index_by_?
        - get_?_by_index

    INPUT:
      - list
      - int
      - 1d-array
      - None
    OUTPUT
      - 1d-array

    Assumes dtype='int32'
    """
    if value is None:
        return default
    if isinstance(value, integer_types):
        return array([value], dtype='int32')
    else:
        return asarray(value)
