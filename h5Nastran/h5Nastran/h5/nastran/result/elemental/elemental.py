from __future__ import print_function, absolute_import

from h5Nastran.h5nastrannode import H5NastranNode
from .element_force import ElementForce
from .energy import Energy
from .strain import Strain
from .stress import Stress


class Elemental(H5NastranNode):
    def __init__(self, h5n, result):
        self._h5n = h5n
        self._result = result

        self.element_force = ElementForce(self._h5n, self)
        self.energy = Energy(self._h5n, self)
        self.strain = Strain(self._h5n, self)
        self.stress = Stress(self._h5n, self)

    def path(self):
        return self._h5n.path() + ['ELEMENTAL']
