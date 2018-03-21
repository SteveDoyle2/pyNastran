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
from .partition import Partition
from .table import Table
from .dynamic import Dynamic
from .design import Design
from ..h5nastrannode import H5NastranNode


class Input(H5NastranNode):
    def __init__(self, h5n):
        self._h5n = h5n  # type: H5Nastran

        self.constraint = Constraint(self._h5n, self)
        # self.contact = Contact(self.h5n, self)
        self.coordinate_system = CoordinateSystem(self._h5n, self)
        self.design = Design(self._h5n, self)
        # self.domains = None
        self.dynamic = Dynamic(self._h5n, self)
        self.element = Element(self._h5n, self)
        # self.fatigue = None
        self.load = Load(self._h5n, self)
        self.material = Material(self._h5n, self)
        # self.matrix = None
        # self.modules = None
        self.node = Node(self._h5n, self)
        self.parameter = Parameter(self._h5n, self)
        self.partition = Partition(self._h5n, self)
        self.property = Property(self._h5n, self)
        self.table = Table(self._h5n, self)
        # self.uds = None

    def path(self):
        return self._h5n.path() + ['INPUT']
            
    def update(self):
        self.coordinate_system.update()

    def to_bdf(self, bdf):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            if hasattr(item, 'to_bdf'):
                item.to_bdf(bdf)

