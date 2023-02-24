from itertools import count

import numpy as np
from numpy import zeros, searchsorted, ravel

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class NonlinearGapStressArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def is_stress(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        msgs = [
            '                      S T R E S S E S   ( F O R C E S )   I N   G A P   E L E M E N T S      ( C G A P )'
            ' '
            '    ELEMENT   - F O R C E S   I N   E L E M   S Y S T -       - D I S P L A C E M E N T S   I N   E L E M   S Y S T -'
            '       ID       COMP-X       SHEAR-Y       SHEAR-Z       AXIAL-U       TOTAL-V       TOTAL-W       SLIP-V        SLIP-W    STATUS'
            #'       3801   3.71080E+05   0.0           0.0           2.37879E-01   9.51516E-01  -5.55112E-17   9.51516E-01  -5.55112E-17 SLIDE   '
        ]
        return msgs

    def get_headers(self) -> list[str]:
        headers = ['compX', 'shearY', 'shearZ', 'axialU', 'shearV', 'shearW', 'slipV', 'slipW']
        return headers

    def build(self):
        """sizes the vectorized attributes of the NonlinearGapStressArray"""
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        #if self.element_type == None?:
            #nnodes_per_element = 1
        #else:
            #raise NotImplementedError(self.element_type)

        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        if self.is_sort1:
            ntimes = self.ntimes
            ntotal = self.ntotal
        else:
            #print("NonlinearGapStressArray: name=%s type=%s ntimes=%s nelements=%s ntotal=%s" % (
                #self.element_name, self.element_type,
                #self.ntimes, self.nelements, self.ntotal))
            ntotal = self.ntimes
            ntimes = self.ntotal
            nelements = self.ntimes
            self.ntimes = ntimes
            self.nelements = nelements
            #print("-> ntimes=%s nelements=%s ntotal=%s" % (
                #ntimes, nelements, ntotal))

        _times = zeros(ntimes, dtype=self.analysis_fmt)
        element = zeros(ntotal, dtype=idtype)
        form = zeros(ntotal, dtype='|U8')

        # [comp_x, shear_y, shear_z, axial_u, shear_v, shear_w, slip_v, slip_w]
        data = zeros((ntimes, ntotal, 8), dtype=fdtype)

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element = group.create_dataset('element', data=element)
            self.data = group.create_dataset('data', data=data)
            self.form = group.create_dataset('form', data=form)
        else:
            self._times = _times
            self.element = element
            self.data = data
            self.form = form

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Time                 0.025     0.050
            #ElementID Item
            #899       compX   0.537281  0.851454
            #          shearY  0.000000  0.000000
            #          shearZ  0.000000  0.000000
            #          axialU  0.000005  0.000009
            #          shearV  0.000000  0.000000
            #          shearW  0.000000  0.000000
            #          slipV   0.000000  0.000000
            #          slipW   0.000000  0.000000
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            data_frame = pd.Panel(self.data,
                                  major_axis=self.element, minor_axis=headers).to_frame()
            data_frame.columns.names = ['Static']
            data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.form, table.form):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += f'cur_form: {self.form}\n'
            msg += f'old_form: {table.form}\n'
            msg += '%s\n' % str(self.code_information())

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

                        (axial_stress1, equiv_stress1, total_strain1,
                         effective_plastic_creep_strain1, effective_creep_strain1,
                         linear_torsional_stress1) = t1

                        (axial_stress2, equiv_stress2, total_strain2,
                         effective_plastic_creep_strain2, effective_creep_strain2,
                         linear_torsional_stress2) = t2

                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += ('%s\n  (%s, %s, %s, %s, %s, %s)\n  '
                                    '(%s, %s, %s, %s, %s, %s)\n' % (
                                        eid,
                                        axial_stress1, equiv_stress1, total_strain1,
                                        effective_plastic_creep_strain1, effective_creep_strain1,
                                        linear_torsional_stress1,

                                        axial_stress2, equiv_stress2, total_strain2,
                                        effective_plastic_creep_strain2, effective_creep_strain2,
                                        linear_torsional_stress2))
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

    def add_sort1(self, dt, eid, comp_xi, shear_yi, shear_zi, axial_ui,
                  shear_vi, shear_wi, slip_vi, slip_wi, form):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.form[self.itotal] = form
        self.data[self.itime, self.itotal, :] = [comp_xi, shear_yi, shear_zi, axial_ui,
                                                 shear_vi, shear_wi, slip_vi, slip_wi]
        self.itotal += 1
        self.ielement += 1

    def add_sort2(self, dt, eid, comp_xi, shear_yi, shear_zi, axial_ui,
                  shear_vi, shear_wi, slip_vi, slip_wi, form):
        """unvectorized method for adding SORT2 transient data"""
        assert self.sort_method == 2, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        itime = self.itotal
        itotal = self.itime
        ntimes = len(self._times)
        ntotal = len(self.element)
        #print(f'NonlinearGapStressArray: itime={itime}/{ntimes} itime={itotal}/{ntotal} -> eid={eid} dt={dt}')
        self._times[itime] = dt
        self.element[itotal] = eid
        self.form[itotal] = form
        self.data[itime, itotal, :] = [comp_xi, shear_yi, shear_zi, axial_ui,
                                       shear_vi, shear_wi, slip_vi, slip_wi]
        self.itotal += 1
        #self.ielement += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.ntotal
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append(f'  type={self.__class__.__name__} ntimes={ntimes:d} nelements={nelements:d}; table_name={self.table_name}\n')
            ntimes_word = 'ntimes'
        else:
            msg.append(f'  type={self.__class__.__name__} nelements={nelements:d}; table_name={self.table_name}\n')
            ntimes_word = '1'
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element.shape = {self.element.shape}\n')
        uform = np.unique(self.form)
        msg.append(f'  form.shape = {self.form.shape}; unique={uform}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1,
                  is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_msgs()
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #comp_x, shear_y, shear_z, axial_u, shear_v, shear_w, slip_v, slip_w
            comp_x = self.data[itime, :, 0]
            shear_y = self.data[itime, :, 1]
            shear_z = self.data[itime, :, 2]
            axial_u = self.data[itime, :, 3]
            shear_v = self.data[itime, :, 4]
            shear_w = self.data[itime, :, 5]
            slip_v = self.data[itime, :, 6]
            slip_w = self.data[itime, :, 7]
            for (i, eid, comp_xi, shear_yi, shear_zi, axial_ui, shear_vi, shear_wi, slip_vi, slip_wi) in zip(
                count(), eids, comp_x, shear_y, shear_z, axial_u, shear_v, shear_w, slip_v, slip_w):

                vals = [comp_xi, shear_yi, shear_zi, axial_ui, shear_vi, shear_wi, slip_vi, slip_wi]
                vals2 = write_floats_13e(vals)
                [comp_xi, shear_yi, shear_zi, axial_ui,
                 shear_vi, shear_wi, slip_vi, slip_wi] = vals2
                f06_file.write(
                    '0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                    % (eid, comp_xi, shear_yi, shear_zi, axial_ui,
                       shear_vi, shear_wi, slip_vi, slip_wi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor in (None, np.nan):
            page_num -= 1
        return page_num
