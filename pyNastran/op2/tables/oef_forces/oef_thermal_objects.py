#pylint disable=C0103,C0301
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import zeros, empty, array_equal
from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import get_key0, write_float_13e, write_floats_13e, _eigenvalue_header
from pyNastran.op2.result_objects.element_table_object import RealElementTableArray
import numpy as np
try:
    import pandas as pd
except ImportError:
    pass


class Real1DHeatFluxArray(ScalarObject):  # 1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itotal = 0
        self.ielement = 0

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux'
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')
        self.element_data_type = empty(self.nelements, dtype='|U8')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        assert 0 not in self.element
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, EType\n'
            for (eid, etype, eid2, etype2) in zip(self.element, self.element_data_type, table.element, table.element_data_type):
                msg += '(%s, %s), (%s, %s)\n' % (eid, etype, eid2, etype2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (xgrad1, ygrad1, zgrad1, xflux1, yflux1, zflux1) = t1
                    (xgrad2, ygrad2, zgrad2, xflux2, yflux2, zflux2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                xgrad1, ygrad1, zgrad1, xflux1, yflux1, zflux1,
                                xgrad2, ygrad2, zgrad2, xflux2, yflux2, zflux2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux):
        self.add_sort1(dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux)

    def add_sort1(self, dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.element_data_type[self.ielement] = etype
        self.data[self.itime, self.ielement, :] = [xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.ielement += 1

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
        msg_temp = [
            '                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S  \n'
            ' \n'
            '    ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n'
            #'            10    ROD         -1.889713E+02                                       3.779427E+04'
        ]

        ntimes = self.data.shape[0]

        eids = self.element
        etype = self.element_data_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            xgrad = self.data[itime, :, 0]
            #ygrad = self.data[itime, :, 1]
            #zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 1]
            #yflux = self.data[itime, :, 4]
            #zflux = self.data[itime, :, 5]

            for (eid, etypei, xgradi, xfluxi) in zip(eids, etype, xgrad, xflux):
                (sxgradi, sxfluxi) = write_floats_13e([xgradi, xfluxi])

                # TODO: hopa is probably the wrong type
                f.write(' %8i  %8s %-13s %-13s %-13s %s\n' % (
                        eid, etypei, sxgradi, '', '', sxfluxi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealHeatFluxVU3DArray(ScalarObject):  # 189-VUQUAD 190-VUTRIA,191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.ielement = 0
        self.itotal = 0
        self.itime = 0

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux',
        ]
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_parent = zeros((self.nelements, 2), dtype='int32')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.vugrid = zeros((self.ntimes, self.ntotal, 1), dtype='int32')
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        # TODO: fix me
        headers = self.get_headers()
        #assert 0 not in self.element
        element_node = [
            self.element_node[:, 0],
            self.element_node[:, 1],
        ]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Node', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.element_parent, table.element_parent):
            assert self.element_parent.shape == table.element_parent.shape, 'element_parent shape=%s table.shape=%s' % (
                self.element_parent.shape, table.element_parent.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Parent, Coord, iCoord\n'
            for (eid1, parent1, coord1, icord1), (eid2, parent2, coord2, icord2) in zip(self.element_parent, table.element_parent_coord_icord):
                msg += '(%s, %s, %s, %s) (%s, %s, %s, %s)\n' % (eid1, parent1, coord1, icord1, eid2, parent2, coord2, icord2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            #eids = self.element_node[:, 0]
            #ntotal = self.data.shape[2]
            for itime in range(self.ntimes):
                vugrids = self.int_data[itime, :, :]
                for j, vugrid in enumerate(vugrids):
                    t1 = self.data[itime, j, :]
                    t2 = table.data[itime, j, :]
                    if not np.array_equal(t1, t2):
                        (xgrad1, ygrad1, zgrad1, xflux1, yflux1, zflux1) = t1
                        (xgrad2, ygrad2, zgrad2, xflux2, yflux2, zflux2) = t2
                        msg += (
                            '(%s, %s)   (%s, %s, %s, %s, %s, %s) (%s, %s, %s, %s, %s, %s)\n' % (
                                j, vugrid,
                                xgrad1, ygrad1, zgrad1, xflux1, yflux1, zflux1,
                                xgrad2, ygrad2, zgrad2, xflux2, yflux2, zflux2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, parent, grad_fluxes):
        self._times[self.itime] = dt
        #icord,
        #print([eid, parent])
        self.element_parent[self.ielement, :] = [eid, parent]
        for grad_flux in grad_fluxes:
            #print(self.itime, self.itotal, grad_flux)
            self.vugrid[self.itime, self.itotal, :] = grad_flux[0]
            self.data[self.itime, self.itotal, :] = grad_flux[1:]
            self.itotal += 1
        self.ielement += 1

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
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        ## TODO: add the f06 header
        msg_temp = [
            '          T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   T R I A N G U L A R   P - E L E M E N T S\n'
            '                 VU-ELEMENT ID=  100005001, P-ELEMENT ID =       5, OUTPUT COORD. ID= (LOCAL), P OF EDGES =  2  2  2\n'  # TODO: wrong
            '                       LOCAL X DIR. = PROJECTED +X DIR.,  LOCAL NORMAL = COUNTER-CLOCKWISE,  ANGLE =    0.0000\n'  # TODO: wrong
            '\n'
            '             VUGRID      X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX    \n'
            #'          111005001     2.000000E+01    -4.799646E-14     0.000000E+00    -4.080000E+03     9.791279E-12     0.000000E+00\n'
        ]
        ntimes = self.data.shape[0]

        #eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # [xgrad, ygrad, zgrad, xflux, yflux, zflux]
            #nids = self.int_data[itime, :, 0]
            vugrids = self.int_data[itime, :, 0]
            print(vugrids)
            xgrad = self.data[itime, :, 0]
            ygrad = self.data[itime, :, 1]
            zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 3]
            yflux = self.data[itime, :, 4]
            zflux = self.data[itime, :, 5]

            for (vugrid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi) in zip(
                 vugrids, xgrad, ygrad, zgrad, xflux, yflux, zflux):
                f.write('         %10i    %-13E    %-13E    %-13E    %-13E    %-13E    %-13E\n' % (
                    vugrid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class HeatFlux_VU_3D(ScalarObject):  # 146-VUPENTA, 147-VUTETRA, 148-VUPENTA
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        #self.eType = {}
        self.parent = {}

        self.grad = {}
        self.flux = {}

        # TODO if dt=None, handle SORT1 case
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self. parent)
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.grad)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, grad, flux\n')
        return msg

    def add_new_transient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, dt, eid, parent, grad_fluxes):
        self.parent[eid] = parent
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[eid][nid] = [xflux, yflux, zflux]

    def add_sort1(self, dt, eid, parent, grad_fluxes):
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]

    def add_sort2(self, eid, dt, parent, grad_fluxes):
        if dt not in self.flux:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        #self.coord[eid] = coord
        #self.icord[eid] = icord
        #self.theta[eid] = theta
        #self.eType[eid] = etype

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]


class RealHeatFluxVUArray(ScalarObject):  # 189-VUQUAD 190-VUTRIA,191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itotal = 0
        self.ielement = 0

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux',
        ]
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_parent_coord_icord = zeros((self.nelements, 4), dtype='int32')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.int_data = zeros((self.ntimes, self.ntotal, 1), dtype='int32')
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def _build_dataframe(self):
        # TODO: fix me
        headers = self.get_headers()
        #assert 0 not in self.element
        element_node = [
            self.element_node[:, 0],
            self.element_node[:, 1],
        ]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Node', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.element_parent_coord_icord, table.element_parent_coord_icord):
            assert self.element_parent_coord_icord.shape == table.element_parent_coord_icord.shape, 'element_parent_coord_icord shape=%s table.shape=%s' % (
                self.element_parent_coord_icord.shape, table.element_parent_coord_icord.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Parent, Coord, iCoord\n'
            for (eid1, parent1, coord1, icord1), (eid2, parent2, coord2, icord2) in zip(self.element_parent_coord_icord, table.element_parent_coord_icord):
                msg += '(%s, %s, %s, %s) (%s, %s, %s, %s)\n' % (eid1, parent1, coord1, icord1, eid2, parent2, coord2, icord2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            #eids = self.element_node[:, 0]
            #ntotal = self.data.shape[2]
            for itime in range(self.ntimes):
                vugrids = self.int_data[itime, :, :]
                for j, vugrid in enumerate(vugrids):
                    t1 = self.data[itime, j, :]
                    t2 = table.data[itime, j, :]
                    if not np.array_equal(t1, t2):
                        (xgrad1, ygrad1, zgrad1, xflux1, yflux1, zflux1) = t1
                        (xgrad2, ygrad2, zgrad2, xflux2, yflux2, zflux2) = t2
                        msg += (
                            '(%s, %s)   (%s, %s, %s, %s, %s, %s) (%s, %s, %s, %s, %s, %s)\n' % (
                                j, vugrid,
                                xgrad1, ygrad1, zgrad1, xflux1, yflux1, zflux1,
                                xgrad2, ygrad2, zgrad2, xflux2, yflux2, zflux2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, parent, coord, icord, theta, grad_fluxes):
        self._times[self.itime] = dt
        #icord,
        #print([eid, parent, coord, theta])
        self.element_parent_coord_icord[self.ielement, :] = [eid, parent, coord, theta]
        for grad_flux in grad_fluxes:
            self.int_data[self.itime, self.itotal, :] = grad_flux[0]
            self.data[self.itime, self.itotal, :] = grad_flux[1:]
            self.itotal += 1
        self.ielement += 1

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
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        ## TODO: add the f06 header
        msg_temp = [
            '          T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   T R I A N G U L A R   P - E L E M E N T S\n'
            '                 VU-ELEMENT ID=  100005001, P-ELEMENT ID =       5, OUTPUT COORD. ID= (LOCAL), P OF EDGES =  2  2  2\n'  # TODO: wrong
            '                       LOCAL X DIR. = PROJECTED +X DIR.,  LOCAL NORMAL = COUNTER-CLOCKWISE,  ANGLE =    0.0000\n'  # TODO: wrong
            '\n'
            '             VUGRID      X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX    \n'
            #'          111005001     2.000000E+01    -4.799646E-14     0.000000E+00    -4.080000E+03     9.791279E-12     0.000000E+00\n'
        ]
        ntimes = self.data.shape[0]

        #eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # [xgrad, ygrad, zgrad, xflux, yflux, zflux]
            #nids = self.int_data[itime, :, 0]
            vugrids = self.int_data[itime, :, 0]
            print(vugrids)
            xgrad = self.data[itime, :, 0]
            ygrad = self.data[itime, :, 1]
            zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 3]
            yflux = self.data[itime, :, 4]
            zflux = self.data[itime, :, 5]

            for (vugrid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi) in zip(
                 vugrids, xgrad, ygrad, zgrad, xflux, yflux, zflux):
                f.write('         %10i    %-13E    %-13E    %-13E    %-13E    %-13E    %-13E\n' % (
                    vugrid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class HeatFlux_VU(ScalarObject):  # 189-VUQUAD 190-VUTRIA,191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}

        self.grad = {}
        self.flux = {}

        # TODO if dt=None, handle SORT1 case
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self. parent)
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.grad)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, coord, icord, theta, grad, flux\n')
        return msg

    def add_new_transient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, nnodes, dt, eid, parent, coord, icord, theta, grad_fluxes):
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

        self.grad[eid] = {}
        self.flux[eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[eid][nid] = [xflux, yflux, zflux]

    def add_sort1(self, nnodes, dt, eid, parent, coord, icord, theta, grad_fluxes):
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]

    def add_sort2(self, nnodes, eid, dt, parent, coord, icord, theta, grad_fluxes):
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]


class RealHeatFluxVUBeamArray(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itotal = 0
        self.ielement = 0
        self.itime = 0

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux',
        ]
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_parent_coord = zeros((self.nelements, 3), dtype='int32')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.vugrid = zeros((self.ntimes, self.ntotal, 1), dtype='int32')
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        #assert 0 not in self.element
        element_node = [
            self.element_parent_coord[:, 0],
            #self.element_parent_coord[:, 1],
            #self.element_parent_coord[:, 2],
        ]
        #print(pd.DataFrame(self.element_parent_coord))
        element_parent_coord = [
            np.vstack([self.element_parent_coord[:, 0], self.element_parent_coord[:, 0]]).T.ravel(),
            np.vstack([self.element_parent_coord[:, 1], self.element_parent_coord[:, 1]]).T.ravel(),
            np.vstack([self.element_parent_coord[:, 2], self.element_parent_coord[:, 2]]).T.ravel(),
        ]
        #print(element_parent_coord)
        if self.nonlinear_factor is not None:
            # TODO: rework
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Node', 'Item']
        else:
            df1 = pd.DataFrame(element_parent_coord).T
            df1.columns = ['Element', 'Parent', 'Coord']
            df2 = pd.DataFrame(self.vugrid[0])
            df2.columns = ['VU_Grid']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            self.data_frame = df1.join([df2, df3])
        #print(self.data_frame)

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.element_parent_coord, table.element_parent_coord):
            assert self.element_parent_coord.shape == table.element_parent_coord.shape, 'element_parent_coord shape=%s table.shape=%s' % (
                self.element_parent_coord.shape, table.element_parent_coord.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Parent, Coord\n'
            for (eid1, parent1, coord1), (eid2, parent2, coord2) in zip(self.element_parent_coord, table.element_parent_coord):
                msg += '(%s, %s, %s)  (%s, %s, %s)\n' % (eid1, parent1, coord1, eid2, parent2, coord2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            eids = self.element_node[:, 0]
            for itime in range(self.ntimes):
                for ie, e in enumerate(eids):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (free_conv1, free_conv_k1) = t1
                    (free_conv2, free_conv_k2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s) (%s, %s)\n' % (
                                eid,
                                free_conv1, free_conv_k1,
                                free_conv2, free_conv_k2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, parent, coord, icord, grad_fluxes):
        self._times[self.itime] = dt
        #icord
        self.element_parent_coord[self.ielement, :] = [eid, parent, coord]
        for grad_flux in grad_fluxes:
            self.vugrid[self.itime, self.itotal, :] = grad_flux[0]
            self.data[self.itime, self.itotal, :] = grad_flux[1:]
            self.itotal += 1
        self.ielement += 1

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
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        asdf
        msg_temp = [
            '                T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   B E A M   P - E L E M E N T S\n'
            '                    VU-ELEMENT ID=  100005001, P-ELEMENT ID =       5, OUTPUT COORD. ID= (LOCAL), P OF EDGES =  2\n'
            '\n'
            '             VUGRID      X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX    \n'
            #'          111005001    -2.000000E+01     0.000000E+00     0.000000E+00     4.080000E+03     0.000000E+00     0.000000E+00\n'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        #eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            vugrids = self.vugrid[itime, :, 0]
            # [xgrad, ygrad, zgrad, xflux, yflux, zflux]
            xgrad = self.data[itime, :, 0]
            ygrad = self.data[itime, :, 1]
            zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 3]
            yflux = self.data[itime, :, 4]
            zflux = self.data[itime, :, 5]

            for (nid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi) in zip(
                 vugrids, xgrad, ygrad, zgrad, xflux, yflux, zflux):
                vals2 = write_floats_13e(
                    [fappliedi, free_convi, force_convi, fradi, ftotali])
                [sfapplied, sfree_conv, sforce_conv, sfrad, sftotal] = vals2

                f.write('         %10i    %13E    %13E    %13E    %13E    %13E    %13E\n' % (
                    eid, sfapplied, sfree_conv, sforce_conv, sfrad, sftotal))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class HeatFlux_VUBEAM(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}

        self.grad = {}
        self.flux = {}

        # TODO if dt=None, handle SORT1 case
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self. parent)
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.grad)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, coord, icord, theta, grad, flux\n')
        return msg

    def add_new_transient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, nnodes, dt, data):
        [eid, parent, coord, icord, grad_fluxes] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[eid][nid] = [xflux, yflux, zflux]

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, coord, icord, grad_fluxes] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]

    def add_sort2(self, nnodes, eid, data):
        [dt, parent, coord, icord, grad_fluxes] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]


class HeatFlux_2D_3DArray(RealElementTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealElementTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def build_dataframe(self):
        headers = self.get_headers()

        #nelements = self.element.shape[0]# // 2
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            df1 = pd.DataFrame(self.element)
            df1.columns = ['ElementID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join(df2)
        #print(self.data_frame)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = [
            '                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S  \n \n',
            '    ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)

    def get_headers(self):
        return ['grad1', 'grad2', 'grad3', 'flux1', 'flux2', 'flux3']


class RealConvHeatFluxArray(ScalarObject):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.ielement = 0
        self.itotal = 0

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'free_conv', 'free_conv_k',
        ]
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.nelements, 2), dtype='int32')

        #[free_conv, free_conv_k]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='float32')

    def build_dataframe(self):
        # TODO: fix me
        headers = self.get_headers()
        #assert 0 not in self.element
        element_node = [
            self.element_node[:, 0],
            self.element_node[:, 1],
        ]
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Node', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            eids = self.element_node[:, 0]
            for itime in range(self.ntimes):
                for ie, e in enumerate(eids):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (free_conv1, free_conv_k1) = t1
                    (free_conv2, free_conv_k2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s) (%s, %s)\n' % (
                                eid,
                                free_conv1, free_conv_k1,
                                free_conv2, free_conv_k2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, cntl_node, free_conv, free_conv_k):
        self.add_sort1(dt, eid, cntl_node, free_conv, free_conv_k)

    def add_sort1(self, dt, eid, cntl_node, free_conv, free_conv_k):
        self._times[self.itime] = dt
        self.element_node[self.ielement, :] = [eid, cntl_node]
        self.data[self.itime, self.ielement, :] = [free_conv, free_conv_k]
        self.ielement += 1

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
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = [
            #'                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S  '
            #' '
            #'    ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX'
            #'             1    QUAD4       -8.372393E-01     1.776357E-15                      8.372393E-01    -1.776357E-15'
            '                                RealConvHeatFluxArray\n'
            '               ELEMENT-ID      FREE-CONVECTION   CONTROL-NODE   FREE-CONVECTION-K\n'
        ]
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # [free_conv, free_conv_k]
            free_conv = self.data[itime, :, 0]
            free_conv_k = self.data[itime, :, 1]

            for (eid, nid, free_convi, free_conv_ki) in zip(eids, nids, free_conv, free_conv_k):
                f.write(' %8i  %-13s %-13s %s\n' % (
                    eid,
                    write_float_13e(free_convi),
                    nid,
                    write_float_13e(free_conv_ki)
                ))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealChbdyHeatFluxArray(ScalarObject):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.ielement = 0
        self.itotal = 0

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'fapplied', 'free_conv', 'force_conv', 'frad', 'ftotal',
        ]
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')
        self.element_type = empty(self.nelements, dtype='|U8')

        #[fapplied, free_conv, force_conv, frad, ftotal]
        self.data = zeros((self.ntimes, self.ntotal, 5), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        assert 0 not in self.element
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    eid = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fapplied1, free_conv1, force_conv1, frad1, ftotal1) = t1
                    (fapplied2, free_conv2, force_conv2, frad2, ftotal2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fapplied1, free_conv1, force_conv1, frad1, ftotal1,
                                fapplied2, free_conv2, force_conv2, frad2, ftotal2,
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal):
        self.add_sort1(dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal)

    def add_sort1(self, dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.element_type[self.ielement] = etype
        self.data[self.itime, self.ielement, :] = [fapplied, free_conv, force_conv, frad, ftotal]
        self.ielement += 1

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
        msg.append('  data.shape = %s\n  ' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        assert self.is_sort1() == True, self.is_sort1()


        msg_temp = [
            '                                H E A T   F L O W   I N T O   H B D Y   E L E M E N T S   (CHBDY)\n'
            ' \n'
            '               ELEMENT-ID      APPLIED-LOAD   FREE-CONVECTION   FORCED-CONVECTION     RADIATION           TOTAL\n'
            #'                       60      0.000000E+00      1.641941E+02      0.000000E+00      0.000000E+00      1.641941E+02'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # [fapplied, free_conv, force_conv, frad, ftotal]
            fapplied = self.data[itime, :, 0]
            free_conv = self.data[itime, :, 1]
            force_conv = self.data[itime, :, 2]
            frad = self.data[itime, :, 3]
            ftotal = self.data[itime, :, 4]

            for (eid, fappliedi, free_convi, force_convi, fradi, ftotali) in zip(
                eids, fapplied, free_conv, force_conv, frad, ftotal):
                #vals2 = write_floats_13e(
                    #[fappliedi, free_convi, force_convi, fradi, ftotali])
                #[sfapplied, sfree_conv, sforce_conv, sfrad, sftotal] = vals2

                f.write('                 %8i     %13E     %13E     %13E     %13E     %13E\n' % (
                    eid, fappliedi, free_convi, force_convi, fradi, ftotali))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
