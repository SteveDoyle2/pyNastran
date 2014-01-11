from pyNastran.utils.mathematics import get_numpy_version

assert get_numpy_version('1.8.0') == (1, 8, 0)
assert get_numpy_version('1.7.0') == (1, 7, 0)
assert get_numpy_version('1.8.1') == (1, 8, 1)
assert get_numpy_version('1.8.2') == (1, 8, 2)
assert get_numpy_version('1.10.11') == (1, 10, 11)
assert get_numpy_version('1.9.0.dev-d1dbf8e') == (1, 9, 0)

assert (2, 0, 0) > (1, 8, 4)
assert (2, 1) > (1, 8, 4)

assert (1, 8, 0) > (1, 8)

print "passed"