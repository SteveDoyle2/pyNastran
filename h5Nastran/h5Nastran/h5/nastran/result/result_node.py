from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range


import tables
import numpy as np

from .elemental import Elemental
from .nodal import Nodal


class ResultNode(object):
    def __init__(self, h5n, nastran):
        self._h5n = h5n  # type: H5Nastran
        self._nastran = nastran

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
        return self._nastran.path() + ['RESULT']

    @property
    def subcases(self):
        # TODO: temporary for now, should be somewhere else, say under h5nastran.result
        return self._h5n.h5f.get_node(self._h5n.table_paths.subcase).read()

