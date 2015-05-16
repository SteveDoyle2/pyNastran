from __future__ import print_function
from six.moves import StringIO
from numpy import (array, searchsorted, array_equal, setdiff1d, int64, argsort,
                   arange, ndarray, nan)
from pyNastran.utils import object_attributes

class VectorizedCard(object):
    def __init__(self, model):
        self.model = model
        self.n = 0
        self.i = 0
        self._comments = {}
        if self.type in model._element_name_to_element_type_mapper:
            self.op2_id = model._element_name_to_element_type_mapper[self.type]
        else:
            if self.type.startswith('C'):
                print('there is no op2_id to apply for element=%r' % self.type)

    def __len__(self):
        return self.n

    def shrink(self, refcheck=True):
        raise NotImplementedError()

    def resize(self, n, refcheck=True):
        names = object_attributes(self, mode="public")
        for name in names:
            attr = getattr(self, name)
            if isinstance(attr, ndarray):
                #self.model.log.info('resizing %r; shape=%s; size=%s' % (name, attr.shape, attr.size))
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

    def __repr__(self):
        f = StringIO()
        self.write_bdf(f)
        return f.getvalue().rstrip()

    def _get_sorted_index(self, sorted_array, unsorted_array, n, field_name, msg, check=True):
        if not array_equal(argsort(sorted_array), arange(len(sorted_array))):
            msg2 = '%s is not sorted\nsorted_array=%s' % (msg, sorted_array)
            raise RuntimeError(msg2)
        assert isinstance(n, int64) or isinstance(n, int), 'field_name=%s n=%s type=%s' % (field_name, n, type(n))
        assert isinstance(check, bool)
        if unsorted_array is None:
            i = slice(None)
            if n == 1:
                i = array([0], dtype='int32')
            return i
        else:
            i = searchsorted(sorted_array, unsorted_array)
            i.astype('int32')
            if check:
                if not array_equal(sorted_array[i], unsorted_array):
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
        if isinstance(i, int64) or isinstance(i, int):
            i = array([i], dtype='int32')
        return i

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
    if isinstance(value, int):
        return array([value], dtype='int32')
    else:
        return asarray(value)
