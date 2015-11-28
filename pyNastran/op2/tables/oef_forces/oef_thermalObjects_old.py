#pylint disable=C0103,C0301
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject


class HeatFlux_2D_3D(ScalarObject):  # 33-QUAD4, 39-TETRA, 53-TRIAX6,64-QUAD8, 67-HEXA, 68-PENTA, 74-TRIA3, 75-TRIA6
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.grad = {}
        self.flux = {}

        # TODO if dt=None, handle SORT1 case
        if is_sort1:
            self.add = self.add_sort1
        else:
            self.add = self.add_sort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self.eType)
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.grad)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, grad, flux\n')
        return msg

    def add_new_transient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add_sort1(self, dt, data):
        [eid, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xgrad, ygrad, zgrad]
        self.flux[dt][eid] = [xflux, yflux, zflux]

    def add_sort2(self, eid, data):
        [dt, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xgrad, ygrad, zgrad]
        self.flux[dt][eid] = [xflux, yflux, zflux]

