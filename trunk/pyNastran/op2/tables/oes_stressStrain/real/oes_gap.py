from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject


class NonlinearGapStress(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.compX = {}
        self.shearY = {}
        self.shearZ = {}
        self.axialU = {}
        self.shearV = {}
        self.shearW = {}
        self.slipV = {}
        self.slipW = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)
        eTypes = list(set(self.eType.values()))
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.compX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW\n')
        msg.append('  eTypes = %s\n' %(', '.join(eTypes)))
        return msg

    def getLength(self):
        return (8, 'f')

    def delete_transient(self, dt):
        del self.compX[dt]

    def get_transients(self):
        k = self.compX.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.element_name = self.data_code['element_name']
        self.dt = dt
        self.compX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.axialU[dt] = {}
        self.shearV[dt] = {}
        self.shearW[dt] = {}
        self.slipV[dt] = {}
        self.slipW[dt] = {}


    def add_new_eid(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        (stress,) = out
        self.eType[eid] = self.element_name
        self.compX[eid] = cpx
        self.shearY[eid] = shy
        self.shearZ[eid] = shz
        self.axialU[eid] = au
        self.shearV[eid] = shv
        self.shearW[eid] = shw
        self.slipV[eid] = slv
        self.slipW[eid] = slp

    def add_new_eid_sort1(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.add_new_transient(dt)
        self.eType[eid] = self.element_name
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp

    def add_new_eid_sort2(self, eid, dt, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.add_new_transient(dt)
        (stress,) = out
        self.eType[eid] = self.element_name
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp