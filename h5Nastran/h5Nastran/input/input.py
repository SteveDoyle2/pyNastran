from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range


import tables
import numpy as np


from .constraint import Constraint
from .coordinate_system import CoordinateSystem
from .element import Element
from .load import Load
from .property import Property
from .node import Node
from .material import Material
from .parameter import Parameter


class Input(object):
    def __init__(self, h5n):
        self._h5n = h5n  # type: H5Nastran

        self.constraint = Constraint(self._h5n, self)
        # self.contact = Contact(self.h5n, self)
        self.coordinate_system = CoordinateSystem(self._h5n, self)
        # self.design = None
        # self.domains = None
        # self.dynamic = None
        self.element = Element(self._h5n, self)
        # self.fatigue = None
        self.load = Load(self._h5n, self)
        self.material = Material(self._h5n, self)
        # self.matrix = None
        # self.modules = None
        self.node = Node(self._h5n, self)
        self.parameter = Parameter(self._h5n, self)
        # self.partition = None
        self.property = Property(self._h5n, self)
        # self.table = None
        # self.uds = None

    def path(self):
        return self._h5n.path() + ['INPUT']

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass
            
    def update(self):
        self.coordinate_system.update()
