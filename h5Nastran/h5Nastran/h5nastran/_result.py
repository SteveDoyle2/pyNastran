from __future__ import print_function, absolute_import

from collections import OrderedDict

import numpy as np
import tables
from six import iteritems

from ._punch import H5NastranResultPunch
from ._op2 import H5NastranResultOP2


class H5NastranResult(H5NastranResultPunch, H5NastranResultOP2):
    def __init__(self, *args, **kwargs):
        super(H5NastranResult, self).__init__(*args, **kwargs)
