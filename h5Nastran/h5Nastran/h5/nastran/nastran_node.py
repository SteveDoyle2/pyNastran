from __future__ import print_function, absolute_import

from six import iteritems

from h5Nastran.h5nastrannode import H5NastranNode
from .input import InputNode
from .result import ResultNode


class NastranNode(H5NastranNode):
    def __init__(self, h5n):
        self._h5n = h5n  # type: H5Nastran
        
        self.input = InputNode(self._h5n, self)
        self.result = ResultNode(self._h5n, self)

    def path(self):
        return self._h5n.path() + ['NASTRAN']

    def update(self):
        self.input.update()

