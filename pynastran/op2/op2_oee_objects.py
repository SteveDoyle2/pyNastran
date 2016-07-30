from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import itervalues
from math import isnan
from collections import OrderedDict
from numpy import zeros, empty, array_equal
import numpy as np

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import get_key0, _eigenvalue_header, write_float_13e
try:
    import pandas as pd
except ImportError:
    pass


class RealStrainEnergy(ScalarObject):
    """
    ::

                                 E L E M E N T   S T R A I N   E N E R G I E S

      ELEMENT-TYPE = QUAD4      * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   9.817708E+08
      SUBCASE               1   * TOTAL ENERGY OF ALL ELEMENTS IN SET       1 =   4.192036E+08

         ELEMENT-ID   STRAIN-ENERGY  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
                 12   2.291087E+07        2.3336            2.291087E+02
                 13   1.582968E+07        1.6124            1.055312E+02
                 14   6.576075E+07        6.6982            3.288037E+02
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.energy = {}
        self.percent = {}
        self.density = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.energy)
            time0 = get_key0(self.energy)
            nelements = len(self.energy[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.energy)
            msg.append('  type=%s nelements=%s\n'
                       % (self.__class__.__name__, nelements))
        msg.append('  energy, percent, density\n  ')
        msg += self.get_data_code()
        return msg

    def update_dt(self, data_code, dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.data_code = data_code
        self.apply_data_code()
        #assert dt >= 0.
        self.log.debug("updating %s...%s=%s  isubcase=%s" % (
            self.name, self.name, dt, self.isubcase))
        if dt is not None:
            self.dt = dt
            self.add_new_transient(dt)
        self.updateNumWide()

    def delete_transient(self, dt):
        del self.energy[dt]
        del self.percent[dt]
        del self.density[dt]

    def get_transients(self):
        k = self.energy.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.energy[dt] = {}
        self.percent[dt] = {}
        self.density[dt] = {}

    def add(self, dt, eid, energy, percent, density):
        #print "energyGridIDs = %s" % (self.energy.keys())
        #assert grid not in self.energy,'grid=%s out=%s' % (grid, out)
        if isinstance(eid, int) and eid <= 0:
            raise ValueError("Invalid Grid ID: eid=%s" % eid)
        self.energy[eid] = energy
        self.percent[eid] = percent
        self.density[eid] = density

    def add_sort1(self, dt, eid, energy, percent, density):
        if dt not in self.energy:
            self.add_new_transient(dt)

        #print str(self)
        #assert grid not in self.energy[dt],'grid=%s dt=%s energy=%s percent=%s density=%s' % (grid, dt, energy, percent, density)
        if eid <= 0:
            raise ValueError("Invalid Grid ID: eid=%s" % eid)

        self.energy[dt][eid] = energy
        self.percent[dt][eid] = percent
        self.density[dt][eid] = density


class RealStrainEnergyArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itime = None
        self.itotal2 = 0
        #self.element_name_count = OrderedDict()
        self.dt_temp = None

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'strain_energy', 'percent', 'strain_energy_density'
        ]
        return headers

    def build(self):
        if self.is_built:
            return
        del self.dt_temp

        #print(self._ntotals)
        self.ntotal = max(self._ntotals)  # TODO: is this correct???

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements = self.ntotal // self.ntimes
        self.nelements = self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.itotal2 = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')
        #self.element_data_type = empty(self.nelements, dtype='|U8')

        #[energy, percent, density]
        assert isinstance(self.ntimes, int), self.ntimes
        assert isinstance(self.ntotal, int), self.ntotal
        self.data = zeros((self.ntimes, self.nelements, 3), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    #def finalize(self):
        #self._times = self._times[::self.nelement_types]
        #ntimes = len(self._times)
        #print('times =', self._times)
        #ntypes_ntimes, _nelements, nheaders = self.data.shape
        #nelements = ntypes_nelements // self.nelement_types
        #nelements = self.nelements
        #nheaders = len(self.get_headers())

        #print('self.nelements =', self.nelements)
        #print('data1.shape =', self.data.shape)
        #print(self.data[:10, :, 0])
        #assert nelements == 29
        #assert ntimes == 3
        #print('ntimes=%s nelements=%s nheaders=%s' % (
            #ntimes, nelements, nheaders))
        #self.data = self.data.reshape(ntimes, nelements, nheaders)
        #print('data2.shape =', self.data.shape)

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'class_name=%s times=%s table.times=%s' % (
                self.class_name, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for (eid, eid2) in zip(self.element, table.element):
                msg += '(%s), (%s)\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (energyi1, percenti1, densityi1) = t1
                    (energyi2, percenti2, densityi2) = t2

                    if np.isnan(densityi1) or not np.isfinite(densityi1):
                        if not np.array_equal(t1[:2], t2[:2]):
                            msg += (
                                '%s (%s, %s)\n'
                                '%s (%s, %s)\n' % (
                                    eid, energyi1, percenti1,
                                    ' ' * len(str(eid)),
                                    energyi2, percenti2,
                                ))
                            i += 1
                            if i > 10:
                                print(msg)
                                raise ValueError(msg)
                    elif not np.array_equal(t1, t2):
                        msg += (
                            '%s (%s, %s, %s)\n'
                            '%s (%s, %s, %s)\n' % (
                                eid, energyi1, percenti1, densityi1,
                                ' ' * len(str(eid)),
                                energyi2, percenti2, densityi2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, energyi, percenti, densityi):
        #itime = self.itime // self.nelement_types
        itime = self.itime
        self._times[itime] = dt
        try:
            self.element[self.ielement] = eid
            #self.element_data_type[self.ielement] = etype
            self.data[itime, self.ielement, :] = [energyi, percenti, densityi]
        except IndexError:
            print('RealStrainEnergyArray', dt, eid, energyi, percenti, densityi)
            raise
        self.ielement += 1
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        """
        '      EIGENVALUE =  2.005177E+05'
        '      CYCLES =  7.126832E+01'
        '                                           E L E M E N T   S T R A I N   E N E R G I E S'
        ' '
        '                ELEMENT-TYPE = TETRA               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   1.002589E+05'
        '                   MODE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   1.002589E+05'
        '0'
        '                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY'
        '                                             4          3.247409E+00                 0.0032              1.948445E+01'
        '                                             5          3.977916E+00                 0.0040              2.386749E+01'
        ''
        '                        TYPE = TETRA    SUBTOTAL        7.225325E+00                 0.0072'
        """
        msg_temp = (
            '                                           E L E M E N T   S T R A I N   E N E R G I E S\n'
            ' \n'
            '                ELEMENT-TYPE = TETRA               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   %s\n'
            '                   MODE        %8i            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =   %s\n'
            '0\n'
            '                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY\n'
        )
        ntimes = self.data.shape[0]

        eids = self.element
        #etype = self.element_data_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            total_energy = 0.
            total_set_energy = 0.
            msg_temp2 = [msg_temp % (total_energy, itime, total_set_energy)]
            f.write(''.join(header + msg_temp2))

            # energy, percent, density
            energy = self.data[itime, :, 0]
            percent = self.data[itime, :, 1]
            density = self.data[itime, :, 2]

            for (eid, energyi, percenti, densityi) in zip(eids, energy, percent, density):
                senergyi = write_float_13e(energyi)
                sdensityi = write_float_13e(densityi)

                f.write(' %8i  %-13s %.4fs %s\n' % (
                    eid, senergyi, percenti, sdensityi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
