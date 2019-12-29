#pylint disable=C0103,C0301
from typing import List

import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import BaseElement
from pyNastran.f06.f06_formatting import (
    write_float_13e, write_floats_13e, _eigenvalue_header)
from pyNastran.op2.result_objects.element_table_object import RealElementTableArray


class Real1DHeatFluxArray(BaseElement):
    """1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itotal = 0
        self.ielement = 0

        if not is_sort1:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        """is the result real?"""
        return True

    @property
    def is_complex(self) -> bool:
        """is the result complex?"""
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux'
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the Real1DHeatFluxArray"""
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element = np.zeros(self.nelements, dtype='int32')
        self.element_data_type = np.empty(self.nelements, dtype='|U8')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        assert 0 not in self.element
        if self.nonlinear_factor not in (None, np.nan):
            #LoadStep                  1.0
            #ElementID Item
            #14        xgrad  0.000000e+00
            #          ygrad  1.401298e-45
            #          zgrad  1.401298e-45
            #          xflux -0.000000e+00
            #          yflux  1.401298e-45
            #          zflux  1.401298e-45
            #15        xgrad -2.842171e-14
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            data_frame = pd.Panel(self.data,
                                  major_axis=self.element,
                                  minor_axis=headers).to_frame()
            data_frame.columns.names = ['Static']
            data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, EType\n'
            for (eid, etype, eid2, etype2) in zip(self.element, self.element_data_type,
                                                  table.element, table.element_data_type):
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

    def add_sort1(self, dt, eid, etype, xgrad, ygrad, zgrad, xflux, yflux, zflux):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.element_data_type[self.ielement] = etype
        self.data[self.itime, self.ielement, :] = [xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
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
            f06_file.write(''.join(header + msg_temp))

            xgrad = self.data[itime, :, 0]
            #ygrad = self.data[itime, :, 1]
            #zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 1]
            #yflux = self.data[itime, :, 4]
            #zflux = self.data[itime, :, 5]

            for (eid, etypei, xgradi, xfluxi) in zip(eids, etype, xgrad, xflux):
                (sxgradi, sxfluxi) = write_floats_13e([xgradi, xfluxi])

                # TODO: hopa is probably the wrong type
                f06_file.write(' %8i  %8s %-13s %-13s %-13s %s\n' % (
                    eid, etypei, sxgradi, '', '', sxfluxi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealHeatFluxVU3DArray(BaseElement):
    """189-VUQUAD 190-VUTRIA,191-VUBEAM"""
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
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

    def get_headers(self) -> List[str]:
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux',
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealHeatFluxVU3DArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element_parent = np.zeros((self.nelements, 2), dtype='int32')

        self.vugrid = np.zeros((self.ntimes, self.ntotal), dtype='int32')
        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def _build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        # TODO: fix me
        headers = self.get_headers()
        #assert 0 not in self.element
        element_parent = [
            self.element_parent[:, 0],
            self.element_parent[:, 1],
        ]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element_parent,
                                       minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
        else:
            self.data_frame = pd.Panel(self.data,
                                       major_axis=element_parent,
                                       minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
        self.data_frame.index.names = ['ElementID', 'Parent', 'Item']

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.element_parent, table.element_parent):
            assert self.element_parent.shape == table.element_parent.shape, 'element_parent shape=%s table.shape=%s' % (
                self.element_parent.shape, table.element_parent.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Parent, Coord, iCoord\n'
            for (eid1, parent1, coord1, icord1), (eid2, parent2, coord2, icord2) in zip(
                    self.element_parent, table.element_parent_coord_icord):
                msg += '(%s, %s, %s, %s) (%s, %s, %s, %s)\n' % (
                    eid1, parent1, coord1, icord1,
                    eid2, parent2, coord2, icord2)
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
        """unvectorized method for adding SORT1 transient data"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        #icord,
        self.element_parent[self.ielement, :] = [eid, parent]
        #try:
            #self.element_parent[self.ielement, :] = [eid, parent]
            #print([self.ielement, eid, parent])
        #except:
            #print(['*', self.ielement, eid, parent])

        for grad_flux in grad_fluxes:
            #print(self.itime, self.itotal, grad_flux)
            self.vugrid[self.itime, self.itotal] = grad_flux[0]
            self.data[self.itime, self.itotal, :] = grad_flux[1:]
            self.itotal += 1
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
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
        #vu3d
        ntimes = self.data.shape[0]

        #eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            # [xgrad, ygrad, zgrad, xflux, yflux, zflux]
            #nids = self.int_data[itime, :, 0]
            #self.element_parent = np.zeros((self.nelements, 2), dtype='int32')
            #self.vugrid = np.zeros((self.ntimes, self.ntotal), dtype='int32')
            vugrids = self.vugrid[itime, :]
            #print(vugrids)
            xgrad = self.data[itime, :, 0]
            ygrad = self.data[itime, :, 1]
            zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 3]
            yflux = self.data[itime, :, 4]
            zflux = self.data[itime, :, 5]

            for (vugrid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi) in zip(
                    vugrids, xgrad, ygrad, zgrad, xflux, yflux, zflux):
                f06_file.write(
                    '         %10i    %-13E    %-13E    %-13E    %-13E    %-13E    %-13E\n' % (
                        vugrid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealHeatFluxVUBeamArray(BaseElement):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.itotal = 0
        self.ielement = 0
        self.itime = 0

        if not is_sort1:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux',
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealHeatFluxVUBeamArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element_parent_coord = np.zeros((self.nelements, 3), dtype='int32')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.vugrid = np.zeros((self.ntimes, self.ntotal, 1), dtype='int32')
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
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
        if self.nonlinear_factor not in (None, np.nan):
            # TODO: rework
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element_node,
                                       minor_axis=headers).to_frame()
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

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.element_parent_coord, table.element_parent_coord):
            assert self.element_parent_coord.shape == table.element_parent_coord.shape, 'element_parent_coord shape=%s table.shape=%s' % (
                self.element_parent_coord.shape, table.element_parent_coord.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Parent, Coord\n'
            for (eid1, parent1, coord1), (eid2, parent2, coord2) in zip(self.element_parent_coord, table.element_parent_coord):
                msg += '(%s, %s, %s)  (%s, %s, %s)\n' % (
                    eid1, parent1, coord1, eid2, parent2, coord2)
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

    def add_sort1(self, dt, eid, parent, coord, unused_icord, grad_fluxes):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element_parent_coord[self.ielement, :] = [eid, parent, coord]
        for grad_flux in grad_fluxes:
            self.vugrid[self.itime, self.itotal, :] = grad_flux[0]
            self.data[self.itime, self.itotal, :] = grad_flux[1:]
            self.itotal += 1
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #vubeam
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
            f06_file.write(''.join(header + msg_temp))

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
                vals2 = write_floats_13e([xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi])
                [sxgradi, sygradi, szgradi, sxfluxi, syfluxi, szfluxi] = vals2
                f06_file.write('         %10i    %13s    %13s    %13s    %13s    %13s    %s\n' % (
                    nid, sxgradi, sygradi, szgradi, sxfluxi, syfluxi, szfluxi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealHeatFlux_2D_3DArray(RealElementTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealElementTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()

        #nelements = self.element.shape[0]# // 2
        if self.nonlinear_factor not in (None, np.nan):
            #Time            0.0           10.0
            #ElementID Item
            #1         grad1   0.0 -1.734723e-18
            #          grad2   0.0 -1.301043e-18
            #          grad3   0.0  1.951564e-18
            #          flux1  -0.0  3.538836e-16
            #          flux2  -0.0  2.654127e-16
            #          flux3  -0.0 -3.981190e-16
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            df1 = pd.DataFrame(self.element)
            df1.columns = ['ElementID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            data_frame = df1.join(df2)
        #print(self.data_frame)
        self.data_frame = data_frame

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = [
            '                   F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S  \n \n',
            '    ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX           Z-FLUX\n']
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)

    def get_headers(self) -> List[str]:
        return ['grad1', 'grad2', 'grad3', 'flux1', 'flux2', 'flux3']


class RealConvHeatFluxArray(BaseElement):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.ielement = 0
        self.itotal = 0

        if not is_sort1:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        """is the result real?"""
        return True

    @property
    def is_complex(self) -> bool:
        """is the result complex?"""
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'free_conv', 'free_conv_k',
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealConvHeatFluxArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element_node = np.zeros((self.nelements, 2), dtype='int32')

        #[free_conv, free_conv_k]
        self.data = np.zeros((self.ntimes, self.ntotal, 2), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        # TODO: fix me
        headers = self.get_headers()
        #assert 0 not in self.element
        element_node = [
            self.element_node[:, 0],
            self.element_node[:, 1],
        ]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            #data_frame = self._build_pandas_transient_elements(
                #column_values, column_names,
                #headers, self.element, self.data)
            #print(data_frame)
            #asdf
            data_frame = pd.Panel(self.data, items=column_values,
                                  major_axis=element_node,
                                  minor_axis=headers).to_frame()
            data_frame.columns.names = column_names
        else:
            # >=25.0
            #Static            free_conv  free_conv_k
            #ElementID NodeID
            #1         0       -0.166667         10.0
            #2         0       -0.166667         10.0
            #3         0       -0.166667         10.0
            #4         0       -0.166667         10.0
            #5         0       -0.166667         10.0
            #6         0       -0.166667         10.0
            # <v24.2
            #Static                              0
            #ElementID Node Item
            #1         0    free_conv    -0.166667
            #               free_conv_k  10.000000
            #2         0    free_conv    -0.166667
            #               free_conv_k  10.000000
            index = pd.MultiIndex.from_arrays(self.element_node.T, names=['ElementID', 'NodeID'])
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=index)
            data_frame.columns.names = ['Static']
            #data_frame = pd.Panel(self.data,
                                  #major_axis=element_node,
                                  #minor_axis=headers).to_frame()
            #data_frame.columns.names = ['Static']
            #data_frame.index.names = ['ElementID', 'Node', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
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

    def add_sort1(self, dt, eid, cntl_node, free_conv, free_conv_k):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element_node[self.ielement, :] = [eid, cntl_node]
        self.data[self.itime, self.ielement, :] = [free_conv, free_conv_k]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
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
            f06_file.write(''.join(header + msg_temp))

            # [free_conv, free_conv_k]
            free_conv = self.data[itime, :, 0]
            free_conv_k = self.data[itime, :, 1]

            for (eid, nid, free_convi, free_conv_ki) in zip(eids, nids, free_conv, free_conv_k):
                f06_file.write(' %8i  %-13s %-13s %s\n' % (
                    eid,
                    write_float_13e(free_convi),
                    nid,
                    write_float_13e(free_conv_ki)
                ))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealChbdyHeatFluxArray(BaseElement):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.ielement = 0
        self.itotal = 0

        if not is_sort1:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        """is the result real?"""
        return True

    @property
    def is_complex(self) -> bool:
        """is the result complex?"""
        return False

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'fapplied', 'free_conv', 'force_conv', 'frad', 'ftotal',
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the RealChbdyHeatFluxArray"""
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element = np.zeros(self.nelements, dtype='int32')
        self.element_type = np.empty(self.nelements, dtype='|U8')

        #[fapplied, free_conv, force_conv, frad, ftotal]
        self.data = np.zeros((self.ntimes, self.ntotal, 5), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        assert 0 not in self.element
        if self.nonlinear_factor not in (None, np.nan):
            #Time                 0.0         10.0
            #ElementID Item
            #10        fapplied     0.0    0.000000
            #          free_conv    0.0  499.376068
            #          force_conv   0.0    0.000000
            #          frad         0.0    0.000000
            #          ftotal       0.0  499.376068
            #20        fapplied     0.0    0.000000
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            # >=25.0
            #Static     fapplied  free_conv  force_conv  frad  ftotal
            #ElementID
            #1          0.166667  -0.166667         0.0   0.0     0.0
            #2          0.166667  -0.166667         0.0   0.0     0.0
            #3          0.166667  -0.166667         0.0   0.0     0.0
            #4          0.166667  -0.166667         0.0   0.0     0.0
            #5          0.166667  -0.166667         0.0   0.0     0.0
            #6          0.166667  -0.166667         0.0   0.0     0.0
            #
            # <=24.2
            #Static                       0
            #ElementID Item
            #1         fapplied    0.166667
            #          free_conv  -0.166667
            #          force_conv  0.000000
            #          frad        0.000000
            #          ftotal      0.000000
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
            #data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            #data_frame.columns.names = ['Static']
            #data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
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
                            '%s   (%s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s)\n' % (
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

    def add_sort1(self, dt, eid, etype, fapplied, free_conv, force_conv, frad, ftotal):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.element_type[self.ielement] = etype
        self.data[self.itime, self.ielement, :] = [fapplied, free_conv, force_conv, frad, ftotal]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        assert self.is_sort1 == True, self.is_sort1


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
            f06_file.write(''.join(header + msg_temp))

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

                f06_file.write('                 %8i     %13E     %13E     %13E     %13E     %13E\n' % (
                    eid, fappliedi, free_convi, force_convi, fradi, ftotali))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

class RealHeatFluxVUShellArray(BaseElement):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = np.nan
        self.table_name = None
        self.approach_code = None
        self.analysis_code = None
        BaseElement.__init__(self, data_code, isubcase, apply_data_code=True)  # no double inheritance
        unused_sort1 = self.is_sort1
        #self.dt = dt
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self.nelements = 0  # result specific

    @property
    def is_real(self) -> bool:
        """is the result real?"""
        return True

    @property
    def is_complex(self) -> bool:
        """is the result complex?"""
        return False

    def data_type(self):
        return 'float32'

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
        #ngrids = len(self.gridTypes)
        msg = []

        unused_ntimesi, ntotal = self.data.shape[:2]
        ntimes = len(self._times)
        nelements = self.element.shape[0]

        nmajor = self.ntimes
        nminor = self.ntotal
        if self.is_sort1:
            assert nmajor == ntimes, 'ntimes=%s expected=%s' % (nmajor, ntimes)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, nelements)
        else:
            assert nmajor == nelements, 'nelements=%s expected=%s' % (nmajor, nelements)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, ntimes)

        msg.append('  isubcase = %s\n' % self.isubcase)
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n'
                       % (self.__class__.__name__, nelements))
        headers = ', '.join(self._get_headers())
        #msg.append('  data: [%s] shape=%s dtype=%s\n'
                   #% (headers, [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  data: [%s] shape=%s dtype=%s\n'
                   % (headers,
                      [int(i) for i in self.data.shape], self.data.dtype))
        msg += self.get_data_code()
        return msg

    @property
    def headers(self):
        return ['xgrad', 'ygrad', 'zgrad', 'xflux', 'yflux', 'zflux']

    def _get_headers(self):
        return self.headers

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        """sizes the vectorized attributes of the ElementTableArray"""
        #print('nelements=%s ntimes=%s sort1?=%s ntotal=%s -> _nelements=%s' % (
            #self.nelements, self.ntimes, self.is_sort1,
            #self.ntotal, self.nelements))

        self.nelements //= self.ntimes
        self.itime = 0
        self.itotal = 0
        self.is_built = True

        if self.is_sort1:
            ntimes = self.ntimes
            nelements = self.ntotal
            nx = ntimes
            ny = self.ntotal
            #print("ntimes=%s nelements=%s" % (ntimes, nelements))
        if self.is_sort2:
            #unused_ntotal = self.ntotal
            nelements = self.ntimes
            ntimes = self.ntotal
            nx = nelements
            ny = ntimes
            #print("ntotal=%s nelements=%s ntimes=%s" % (ntotal, nelements, ntimes))

        self._times = np.zeros(ntimes, dtype=self._times_dtype)
        #self.types = array(self.nelements, dtype='|S1')

        self.element = np.zeros(nelements, dtype='int32')
        self.element_parent_coord_icord = np.zeros((nelements, 4), dtype='int32')
        #self.element_data_type = empty(nelements, dtype='|U8')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.data = np.zeros((nx, ny, 6), self.data_type())

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (self.nonlinear_factor is not None and
                  np.isnan(self.nonlinear_factor) and
                  np.isnan(table.nonlinear_factor))
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not is_nan:
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)

        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'eid:'
            for (eid, eid2) in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)

        if not np.array_equal(self.element_parent_coord_icord, table.element_parent_coord_icord):
            assert self.element_parent_coord_icord.shape == table.element_parent_coord_icord.shape, 'shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'element_parent_coord_icord:'
            for (epci1, epci2) in zip(self.element_parent_coord_icord, table.element_parent_coord_icord):
                msg += '%s, %s\n' % (epci1, epci2)
            print(msg)
            raise ValueError(msg)

        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                tx1, ty1, tz1, rx1, ry1, rz1,
                                tx2, ty2, tz2, rx2, ry2, rz2)
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

    def add_sort1(self, dt, eid, parent, coord, unused_icord, unused_theta,
                  xgrad, ygrad, zgrad, xflux, yflux, zflux):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        #print(eid, parent, coord, icord)
        # icord is a string?
        self.element_parent_coord_icord[self.itotal] = [eid, parent, coord, 0]
        #self.element_data_type[self.itotal] = etype
        self.data[self.itime, self.itotal, :] = [xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.itotal += 1

    #def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  #page_num=1, is_mag_phase=False, is_sort1=True):
        #pass
