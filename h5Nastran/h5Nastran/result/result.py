from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range


import tables
import numpy as np

from .elemental import Elemental
from .nodal import Nodal


class Result(object):
    def __init__(self, h5n):
        self._h5n = h5n  # type: H5Nastran

        self.acoustic = None
        self.contact = None
        self.domains = None
        self.elemental = Elemental(self._h5n, self)
        self.fatigue = None
        self.fatigue_vibration = None
        self.matrix = None
        self.modal = None
        self.nlout = None
        self.nodal = Nodal(self._h5n, self)
        self.optimization = None
        self.summary = None
        self.monitor = None

    def path(self):
        return self._h5n.path() + ['RESULT']
