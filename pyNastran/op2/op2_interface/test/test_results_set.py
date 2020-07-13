import unittest
from pyNastran.op2.op2_interface.result_set import ResultSet


class TestResultsSet(unittest.TestCase):

    def test_result_set_print(self):
        allowed_results = ['a', 'b', 'c']
        results_map = {}
        unused_log = None
        res = ResultSet(allowed_results, results_map, unused_log)
        msg_expected = (
            'ResultSet:\n'
            ' results:\n'
            '  a\n'
            '  b\n'
            '  c\n'
        )
        for a, b in zip(msg_expected.split('\n'), str(res).split('\n')):
            assert a == b, f'expected={a!r} actual={b}'
        assert str(res) == msg_expected
        #print(res)
        res.remove('a')

        msg_expected = (
            'ResultSet:\n'
            ' results:\n'
            '  a (disabled)\n'
            '  b\n'
            '  c\n'
        )
        assert str(res) == msg_expected, str(res)

    def test_result_set_1(self):
        allowed_results = ['chexa_stress', 'strain.chexa_strain', 'displacement', 'cquad4_strain',
                           'strain.ctetra_strain']
        results_map = {}
        unused_log = None
        res = ResultSet(allowed_results, results_map, unused_log)
        with self.assertRaises(RuntimeError):
            matched = res._get_matched_results('chexa_')
        #assert len(matched) == [], matched
        #print(matched)

        matched = res._get_matched_results('chexa_*')
        assert matched == ['chexa_stress'], matched

        matched = res._get_matched_results('chexa_stress')
        assert matched == ['chexa_stress'], matched

        matched = res._get_matched_results('*_strain')
        assert set(matched) == set(['strain.chexa_strain', 'strain.ctetra_strain', 'cquad4_strain']), matched

        matched = res._get_matched_results('strain*')
        assert set(matched) == set(['strain.chexa_strain', 'strain.ctetra_strain']), matched

        matched = res._get_matched_results('.*chexa_*')
        assert set(matched) == set(['chexa_stress', 'strain.chexa_strain']), matched

        matched = res._get_matched_results('.*chexa_.*')
        assert set(matched) == set(['chexa_stress', 'strain.chexa_strain']), matched

        #matched = res._get_matched_results('*chexa_.*')
        #assert set(matched) == set(['chexa_stress', 'strain.chexa_strain']), matched

        matched = res._get_matched_results('strain.ch*')
        assert set(matched) == set(['strain.chexa_strain']), matched

        matched = res._get_matched_results('strain.ch.*')
        assert set(matched) == set(['strain.chexa_strain']), matched

        #print(matched)
        results = ['chexa*']
        res.remove(results)
        #print(res)


if __name__ == '__main__':
    unittest.main()
