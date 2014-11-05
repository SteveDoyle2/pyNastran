from six.moves import StringIO
from numpy import array, searchsorted, array_equal, setdiff1d, int64

class VectorizedCard(object):
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []
        if self.type in model._element_name_to_element_type_mapper:
            self.op2_id = model._element_name_to_element_type_mapper[self.type]
        else:
            if self.type.startswith('C'):
                print('there is no op2_id to apply for element=%r' % self.type)

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % (self.type, self.n))
        return msg

    def __repr__(self):
        f = StringIO()
        self.write_bdf(f)
        return f.getvalue().rstrip()

    def _get_sorted_index(self, sorted_array, unsorted_array, n, msg, check=True):
        assert isinstance(n, int)
        assert isinstance(check, bool)
        if unsorted_array is None:
            i = None
            if n == 1:
                i = array([0], dtype='int32')
        else:
            i = searchsorted(sorted_array, unsorted_array)
            i.astype('int32')
            if check:
                if not array_equal(sorted_array[i], unsorted_array):
                    # undefined nodes/elements
                    #print('unsorted %s' % unsorted_array)
                    #print('sorted %s' % sorted_array)
                    #pass
                    raise RuntimeError('Undefined %s\n%s' % (msg, setdiff1d(unsorted_array, sorted_array)))
        if isinstance(i, int64) or isinstance(i, int):
            i = array([i], dtype='int32')
        return i
