#!/usr/bin/env python

import numpy as np
import op4

fh = op4.File('../test/mat_t_dn.op4', 'r')
a = op4.load(10)

print 'The array created is %s' % a
print 'It carries a reference to our deallocator: %s ' % a.base
# np.testing.assert_allclose(a, np.arange(10))
