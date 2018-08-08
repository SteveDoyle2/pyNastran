from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import numpy as np
try:
    import pandas as pd  # type: ignore
except ImportError:
    pass

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e

ints = (int, np.int32)


class ComplexBendArray(OES_Object):
    """
    Common class for:
     - ComplexBendStressArray
     - ComplexBendStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    #def get_nnodes(self):
        #return get_nnodes(self)

    def build(self):
        """sizes the vectorized attributes of the ComplexBendArray"""
        #print('data_code = %s' % self.data_code)
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = np.zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        self.element_node = np.zeros((self.ntotal, 2), 'int32')

        # the number is messed up because of the offset for the element's properties
        if not self.nelements * nnodes * 2 == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes,
                self.ntotal)
            raise RuntimeError(msg)

        # [angle, sc, sd, se, sf]
        self.data = np.zeros((self.ntimes, self.ntotal, 5), 'complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        headers = self.headers
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values,
                                   major_axis=self.element, minor_axis=self.headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, (eid, nid) in enumerate(self.element_node):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (angle1, sc1, sd1, se1, sf1) = t1
                        (angle2, sc2, sd2, se2, sf2) = t2
                        d = t1 - t2
                        if not np.allclose([angle1.real, sc1.real, sc1.imag, sd1.real, sd1.imag, se1.real, se1.imag, sf1.real, sf1.imag, ],
                                           [angle2.real, sc2.real, sc2.imag, sd2.real, sd2.imag, se2.real, se2.imag, sf2.real, sf2.imag, ], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                sc1.real, sc1.imag, sd1.real, sd1.imag,
                                sc2.real, sc2.imag, sd2.real, sd2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, grid, angle, sc, sd, se, sf):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, int) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [angle, sc, sd, se, sf]
        self.element_node[self.itotal] = [eid, grid]
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  data: [ntimes, nnodes, 5] where 5=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    @property
    def headers(self):
        return self._get_headers()

    def get_headers(self):
        return self.headers

class ComplexBendStressArray(ComplexBendArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['angle', 'sc', 'sd', 'se', 'sf']

class ComplexBendStrainArray(ComplexBendArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBendArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        return ['angle', 'sc', 'sd', 'se', 'sf']
