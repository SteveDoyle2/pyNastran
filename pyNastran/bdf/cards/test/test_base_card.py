from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import unittest
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.cards.expand_card import expand_thru, expand_thru_by, setup_data #, expand_thru_exclude

def get_except(data, i, ndata, end_value):
    removed = set()
    while i < len(data):
        value = data[i]
        #print('  exclude?', i, value)
        if isinstance(value, int):
            ivalue = value
        else:
            ivalue = int(value)
        #print('  ivalue=%s > end_value=%s' %  (ivalue, end_value))
        if ivalue > end_value:
            #print('    break')
            break
        removed.add(ivalue)
        i += 1
    #print('removed...', removed)
    #print('*datai =', data[i])
    return i, removed

def expand(data_in):
    """new expand method"""
    data, stype = setup_data(data_in)
    #print('***************************')
    #print('data =', data)
    assert stype == 'int', data
    ndata = len(data) - 1
    out = []
    removed_set = set()

    i = 0
    is_thru = False
    while i < len(data):
        value = data[i]
        #print('value =', value)
        if isinstance(value, int):
            out.append(value)
            i += 1
            continue

        if is_thru:
            is_thru = False
            continue_after_range = False

            out2 = []
            #print('*', i, value)
            value0 = out.pop()
            end_value = int(value)
            i += 1

            by_value = 1
            if i < ndata:
                next_svalue = data[i]
                #print('found by?', next_svalue)
                if next_svalue == 'BY':
                    i += 1
                    next_svalue2 = data[i]
                    if isinstance(next_svalue2, int):
                        by_value = next_svalue2
                    else:
                        next_svalue2 = next_svalue2.upper()
                        #print('next_svalue2 =', next_svalue2)
                        by_value = int(next_svalue2)
                    #print('by_value =', by_value)
                elif next_svalue == 'EXCEPT':
                    i, removed_seti = get_except(data, i+1, ndata, end_value)
                    removed_set.update(removed_seti)
                    continue_after_range = True
                else:
                    raise RuntimeError('next_svalue=%s data=%s' % (next_svalue,  data))
            rangei = range(value0, end_value + 1, by_value)
            #print('range(%r, %r, %r) = %s' % (value0, end_value + 1, by_value, rangei))
            out2.extend(rangei)
            if end_value != out2[-1]:
                out2.append(end_value)

            if continue_after_range:
                #print('continue_after_range')
                out.extend(out2)
                continue

            i += 1
            if i < ndata:
                exclude_svalue = data[i]
                #print('found exclude?', exclude_svalue)
                if exclude_svalue == 'EXCEPT':
                    i, removed = get_except(data, i, ndata, end_value)
                    removed_set.update(removed_seti)
                    asdf
                    continue
                else:
                    not_except
            else:
                i -= 1

            out.extend(out2)
            #print(i, value0, 'THRU', end_value)
            #print('------------')
            del by_value
            i += 1
        else:
            #print('add', i, value)
            assert is_thru is False, data
            if value == 'THRU':
                is_thru = True
            else:
                ivalue = int(value)
                out.append(ivalue)
            i += 1
    if len(removed_set):
        #print('data =', data)
        #print('out =', out)
        #print('removed_set =', removed_set)
        out_set = set(out) - removed_set
        out = list(out_set)
    return out

#expand_thru = expand
class TestBaseCard(unittest.TestCase):
    """Tests methods used by ``BaseCard``"""
    def test_expand_thru(self):
        """tests expand_thru"""
        values1 = expand_thru(['1', 'THRU', '10'])
        values2 = expand(['1', 'THRU', '10'])
        values3 = expand(['1 THRU 10'])
        assert values1 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values1
        assert values2 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values2
        assert values3 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values3


        values1 = expand_thru(['1', 'thru', '10'])
        values2 = expand(['1', 'thru', '10'])
        values3 = expand(['1 thru 10'])
        values4 = expand(['1, thru, 10'])
        assert values1 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values1
        assert values2 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values2
        assert values3 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values3
        assert values4 == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], values4

        values1 = expand_thru([1, 2, 3, 4, 5])
        values2 = expand([1, 2, 3, 4, 5])
        values3 = expand(['1, 2, 3, 4, 5'])
        values4 = expand(['1 2 3 4 5'])
        assert values1 == [1, 2, 3, 4, 5], values1
        assert values2 == [1, 2, 3, 4, 5], values2
        assert values3 == [1, 2, 3, 4, 5], values3
        assert values4 == [1, 2, 3, 4, 5], values4

        values1 = expand_thru_by(['1', 'THRU', '10', 'BY', 2])
        values2 = expand(['1', 'THRU', '10', 'BY', 2])
        values3 = expand(['1 THRU 10 BY 2'])
        assert values1 == [1, 3, 5, 7, 9, 10], values1
        assert values2 == [1, 3, 5, 7, 9, 10], values2
        assert values3 == [1, 3, 5, 7, 9, 10], values3

        values1 = expand_thru_by(['1', 'thru', '10', 'by', 2])
        values2 = expand(['1', 'thru', '10', 'by', 2])
        values3 = expand(['1 thru 10 by 2'])
        values4 = expand(['1 thru 10, by, 2'])
        assert values1 == [1, 3, 5, 7, 9, 10], values1
        assert values2 == [1, 3, 5, 7, 9, 10], values2
        assert values3 == [1, 3, 5, 7, 9, 10], values3
        assert values4 == [1, 3, 5, 7, 9, 10], values4

        values1 = expand_thru_by([1, 2, 3, 4, 5])
        values2 = expand([1, 2, 3, 4, 5])
        assert values1 == [1, 2, 3, 4, 5], values1
        assert values2 == [1, 2, 3, 4, 5], values2

        #values = expand_thru_exclude(['1', 'thru', '5', 'exclude', 2])
        #assert values == [1, 3, 4, 5], values
        values = expand(['1,5,THRU,11,EXCEPT,7,8,13'])
        assert values == [1, 5, 6, 9, 10, 11, 13], values

    def test_collapse_thru(self):
        """
        tests collapse_thru method used by SETx cards
        """
        data = [1, 2, 3, 4, 5, 10]
        expected = [1, u'THRU', 5, 10]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 17]
        expected = [1, 3, 4, 5, 6, 17]
        msg = 'expected=%s actual=%s' % (expected, collapse_thru_by(data))
        self.assertEqual(collapse_thru_by(data), expected, msg)

        data = [1, 3, 4, 5, 6, 7, 17]
        expected = [1, 3, 4, 'THRU', 7, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 6, 8, 10, 12, 14, 17]
        expected = [1, 3, 4, 'THRU', 14, 'BY', 2, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 101]
        expected = [1, 3, 4, 5, 6, 8, 'THRU', 22, 'BY', 2, 101]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5]
        expected = [1, 'THRU', 5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [5]
        expected = [5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 7, 9, 11, 12, 14, 16]
        expected = [1, 'THRU', 5,
                    7, 9, 11,
                    12, 14, 16]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2]
        expected = [1, 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 5, 7, 9, 11]
        expected = [1, 'THRU', 11, 'BY', 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4]
        expected = [1, 'THRU', 4]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3]
        expected = [1, 2, 3]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = [1, 'THRU', 8]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))


if __name__ == '__main__':
    unittest.main()
