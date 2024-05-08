import numpy as np
from numpy import zeros

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_float_13e, _eigenvalue_header

#BASIC_TABLES = {
    #'OESATO1', 'OESCRM1', 'OESPSD1', 'OESRMS1', 'OESNO1',
    #'OESATO2', 'OESCRM2', 'OESPSD2', 'OESRMS2', 'OESNO2',
    #'OSTRATO1', 'OSTRCRM1', 'OSTRPSD1', 'OSTRRMS1', 'OSTRNO1',
    #'OSTRATO2', 'OSTRCRM2', 'OSTRPSD2', 'OSTRRMS2', 'OSTRNO2',
#}
#VM_TABLES = {'OESXRMS1', 'OESXNO1'}


class RandomPlateVMArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        self.fiber_distance = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return get_nnodes(self)

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def build(self):
        r"""sizes the vectorized attributes of the RandomPlateVMArray

        SORT1:
         - etype     ndata numwide size  -> nelements     nnodes nlayers
         - CQUAD-144 376   47      4        376/(47*4)=2  5      2*2*5=20  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\plate_111.op2

        """
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']

        nnodes = self.nnodes_per_element

        #self.names = []
        self.nelements //= nnodes
        #print('element_type=%r ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, self.ntotal, self.subtitle))
        #self.nelements //= self.ntimes
        #self.ntotal = self.nelements # * 2
        #if self.element_name == 'CTRIA3':
        #print('element_type=%r ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, self.ntotal, self.subtitle))
        #print(self)
        if self.is_sort1:
            # old
            #ntimes = self.ntimes
            #nelements = self.nelements
            ntimes = len(self._ntotals)
            ntotal = self._ntotals[0]
            nelements = ntotal // (2 * nnodes)

            #ntotal = self.ntotal
            #ny = nelements * 2
            nelements_nnodes = nelements * 2
            nlayers = nelements * 2 * nnodes
            if nlayers != ntotal:
                msg = f'nlayers={nlayers} ntotal={ntotal}\n'
                msg += f'SORT1 {self.element_name}-{self.element_type} ntimes={ntimes} nelements={nelements} ntotal={ntotal}\n'
                msg = f'_ntotals={self._ntotals}\n'
                raise RuntimeError(msg)
            #nx = ntimes
            #ny = nelements_nnodes
            #ntotal = nelements * 2
            #if self.element_name in ['CTRIA3', 'CQUAD8']:
            #print(f"SORT1 ntimes={ntimes} nelements={nelements} ntotal={ntotal}")
        elif self.is_sort2:
            # ntotal=164
            # len(_ntotals) = 4580 -> nelements=4580
            # nfreqs=82
            # flip this to sort1?
            #ntimes = self.ntotal
            #nnodes = self.ntimes
            #ntotal = nnodes
            #nelements = self.ntimes
            #ntimes = self.nelements # // nelements
            #ntotal = self.ntotal
            nelements = len(self._ntotals)
            ntotal = self._ntotals[0]
            ntimes = ntotal // 2 // nnodes

            #print(self._ntotals)
            ntotal = self._ntotals[0]
            #nelements = len(self._ntotals)
            #ntimes = ntotal // 2
            #ntimes, nelements = nelements_real, ntimes_real
            #ntotal = self.ntotal
            nlayers = nelements * 2 * nnodes
            #ny = nlayers
            #nx = ntimes
            #if self.element_name in ['CTRIA3', 'CQUAD8']:
            #if self.element_name in ['CQUAD4']:
                #print(f"SORT2 ntimes={ntimes} nelements={nelements} ntotal={ntotal} nnodes={nnodes} nlayers={nlayers}")
        else:  # pragma: no cover
            raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())

        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self.build_data(ntimes, nelements, nlayers, nnodes, self._times_dtype)
        #print(''.join(self.get_stats()))
        #try:
            #name = self.data_code['name']
        #except KeyError:
            #print(''.join(self.get_stats()))
            #raise

    def build_data(self, ntimes, nelements, nlayers, nnodes, dtype):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements
        #ntotal = nelements * 2
        self.ntotal = nlayers
        #_times = zeros(ntimes, dtype=self.analysis_fmt)
        #element = zeros(nelements, dtype=idtype)

        self._times = zeros(ntimes, dtype=self.analysis_fmt)
        #self.ntotal = self.nelements * nnodes

        #print(f'***nelements={nelements} nlayers={nlayers} ntimes={ntimes}')
        self.element_node = zeros((nlayers, 2), 'int32')

        # the number is messed up because of the offset for the element's properties
        #if not self.nelements * 2 == self.ntotal:
            #msg = 'ntimes=%s nelements=%s nlayers=%s ntotal=%s' % (
                #self.ntimes, self.nelements, self.nelements * 2, self.ntotal)
            #raise RuntimeError(msg)

        self.fiber_distance = zeros(nlayers, 'float32')

        # [oxx, oyy, txy]
        #print(f'ntimes={self.ntimes} nelements={self.nelements} ntotal={self.ntotal}')
        self.data = zeros((ntimes, nlayers, 3), 'float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        #print(f'column_names = {column_names} column_values={column_values}')

        #print(self.element_node)
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\psdo7026.op2
        #self.data_frame = pd.Panel(self.data, items=column_values,
                                   #major_axis=self.element_node, minor_axis=headers).to_frame()
        #self.data_frame.columns.names = column_names
        #self.data_frame.index.names = ['ElementID', 'Item']
        #return

        names = ['ElementID', 'NodeID']
        #ipos = np.where(self.element_node[:, 0] > 0)
        ipos = None

        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_element_node(
            column_values, column_names,
            headers, self.element_node, self.data)
        #print(data_frame)
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        assert self.element_node[:, 0].min() > 0, self.element_node
        assert table.element_node[:, 0].min() > 0, table.element_node
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (
                self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid, Nid\n' % str(self.code_information())
            for eid1, eid2 in zip(self.element, table.element):
                msg += '(%s, %s), (%s, %s)\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
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
                        (oxx1, oyy1, txy1) = t1
                        (oxx2, oyy2, txy2) = t2
                        #d = t1 - t2
                        if not np.allclose(
                                [oxx1, oyy1, txy1], # atol=0.0001
                                [oxx2, oyy2, txy2], atol=0.075):
                            ni = len(str(eid)) + len(str(eid))
                        #if not np.array_equal(t1, t2):
                            msg += ('%s  (%s, %s, %s)\n'
                                    '%s  (%s, %s, %s)\n' % (
                                        eid,
                                        oxx1, oyy1, txy1,
                                        ' ' * ni,
                                        oxx2, oyy2, txy2,
                                    ))
                            msg += ('%s  (%s, %s, %s)\n'
                                    % (
                                        ' ' * ni,
                                        oxx1 - oxx2,
                                        oyy1 - oyy2,
                                        txy1 - txy2,
                                    ))

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

    def finalize(self):
        """Calls any OP2 objects that need to do any post matrix calcs"""
        if self.is_sort1:
            return
        #print('finalize random plate')
        self.set_as_sort1()
        #print(self.get_stats())
        #print(self._times, self._times.dtype)
        #print(self.element_node)

    #def add_sort1(self, dt, eid, nid,
                  #fd1, oxx1, oyy1, txy1,
                  #fd2, oxx2, oyy2, txy2):
        #assert self.is_sort1, self.sort_method
        ##assert self.element_node.max() == 0, self.element_node
        ##if self.element_name in ['CTRIA3', 'CQUAD8']:
        ##print(f'SORT1 {self.element_name}: itime={self.itime} ielement={self.ielement} itotal={self.itotal} dt={dt} eid={eid} nid={nid} fd={fd1} oxx={oxx1}')
        ##print(f'SORT1 {self.element_name}: itime={self.itime} ielement={self.ielement} itotal={self.itotal+1} dt={dt} eid={eid} nid={nid} fd={fd2} oxx={oxx2}')
            ##print('%s: itime=%s itotal=%s dt=%s eid=%s fd=%s oxx=%s' % (self.element_name, self.itime, self.itotal, dt, eid, fd, oxx))
        #self._times[self.itime] = dt
        ##print(self.element_types2, element_type, self.element_types2.dtype)

        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #self.data[self.itime, self.itotal, :] = [oxx1, oyy1, txy1]
        #self.element_node[self.itotal, :] = [eid, nid]  # 0 is center
        #self.fiber_distance[self.itotal] = fd1
        ##self.ielement += 1
        #self.itotal += 1

        #self.data[self.itime, self.itotal, :] = [oxx2, oyy2, txy2]
        #self.element_node[self.itotal, :] = [eid, nid]  # 0 is center
        #self.fiber_distance[self.itotal] = fd2
        #self.itotal += 1

    #def _get_sort2_itime_ilower_iupper_from_itotal(self):
        #nnodes = self.nnodes_per_element
        #itime = self.ielement // nnodes
        #inid = self.ielement % nnodes
        #itotal = self.itotal
        ##if itime >= self.data.shape[0]:# or itotal >= self.element_node.shape[0]:
            ##print(f'*SORT2 {self.element_name}: itime={itime} ielement={self.itime} inid={inid} itotal={itotal} dt={dt} eid={eid} nid={nid} fd={fd1:.2f} oxx={oxx1:.2f}')
            ##print(f'*SORT2 {self.element_name}: itime={itime} ielement={self.itime} inid={inid} itotal={itotal+1} dt={dt} eid={eid} nid={nid} fd={fd2:.2f} oxx={oxx2:.2f}')
            ##print(self.data.shape)
            ##print(self.element_node.shape)
        ##else:
        #ielement = self.itime

        ##ibase = 2 * ielement # ctria3/cquad4-33
        #ibase = 2 * (ielement * nnodes + inid)
        #ie_upper = ibase
        #ie_lower = ibase + 1

        ##if self.element_name == 'CTRIAR': # and self.table_name == 'OESATO2':
        ##debug = False
        ##if self.element_name == 'CTRIA3' and self.table_name in ['OSTRRMS1', 'OSTRRMS2']:
            ##debug = True
            ##print(f'SORT2 {self.table_name} {self.element_name}: itime={itime} ie_upper={ie_upper} ielement={self.itime} inid={inid} nid={nid} itotal={itotal} dt={dt} eid={eid} nid={nid} fd={fd1:.2f} oxx={oxx1:.2f}')
            ##print(f'SORT2 {self.table_name} {self.element_name}: itime={itime} ie_lower={ie_lower} ielement={self.itime} inid={inid} nid={nid} itotal={itotal+1} dt={dt} eid={eid} nid={nid} fd={fd2:.2f} oxx={oxx2:.2f}')
        #return itime, ie_upper, ie_lower

    #def add_sort2(self, dt, eid, nid,
                  #fd1, oxx1, oyy1, txy1,
                  #fd2, oxx2, oyy2, txy2):
        #assert self.is_sort2, self
        ##if self.element_name == 'CTRIA3':
        ##assert self.element_node.max() == 0, self.element_node
        ##print(self.element_node, nid)
        #itime, ie_upper, ie_lower = self._get_sort2_itime_ilower_iupper_from_itotal()
        #self._times[itime] = dt
        ##print(self.element_types2, element_type, self.element_types2.dtype)

        ##itime = self.ielement
        ##itime = self.itime
        ##ielement = self.itime
        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)

        #if itime == 0:
            #self.element_node[ie_upper, :] = [eid, nid]  # 0 is center
            #self.element_node[ie_lower, :] = [eid, nid]  # 0 is center
            #self.fiber_distance[ie_upper] = fd1
            #self.fiber_distance[ie_lower] = fd2
            ##if self.element_name == 'CQUAD4':
                ##print(self.element_node)

        #self.data[itime, ie_upper, :] = [oxx1, oyy1, txy1]
        #self.data[itime, ie_lower, :] = [oxx2, oyy2, txy2]

        #self.itotal += 2
        #self.ielement += 1
        #debug = False
        #if debug:
            #print(self.element_node)
    #---------------------------------------------------------------------------

    def add_ovm_sort1(self, dt, eid, nid,
                      fd1, oxx1, oyy1, txy1, ovm1,
                      fd2, oxx2, oyy2, txy2, ovm2):
        """unvectorized method for adding SORT1 transient data"""
        assert self.is_sort1, self.sort_method
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #if self.element_name == 'CTRIA3':
        #print('%s itotal=%s dt=%s eid=%s nid=%-5s oxx=%s' % (
            #self.element_name, self.itotal, dt, eid, nid, oxx1))

        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.data[self.itime, self.itotal, :] = [oxx1, oyy1, txy1, ovm1]
        self.element_node[self.itotal, :] = [eid, nid]  # 0 is center
        self.fiber_distance[self.itotal] = fd1
        #self.ielement += 1
        self.itotal += 1

        #print(self.data.shape)
        self.data[self.itime, self.itotal, :] = [oxx2, oyy2, txy2, ovm2]
        self.element_node[self.itotal, :] = [eid, nid]  # 0 is center
        self.fiber_distance[self.itotal] = fd2
        self.itotal += 1
        #print(self.element_node)

    def add_ovm_sort2(self, dt, eid, nid,
                      fd1, oxx1, oyy1, txy1, ovm1,
                      fd2, oxx2, oyy2, txy2, ovm2):
        """unvectorized method for adding SORT2 transient data"""
        #self.add_sort2(dt, eid, nid, fd1, oxx1, oyy1, txy1, fd2, oxx2, oyy2, txy2)
        self.add_ovm_sort1(dt, eid, nid,
                           fd1, oxx1, oyy1, txy1, ovm1,
                           fd2, oxx2, oyy2, txy2, ovm2)

    #---------------------------------------------------------------------------

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append(f'  type={self.class_name} ntimes={ntimes} nelements={nelements:d} nnodes={nnodes:d} table_name={self.table_name}\n')
        else:
            msg.append(f'  type={self.class_name} nelements={nelements:d} nnodes={nnodes:d} {self.table_name}\n')
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 3] where 3=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        #print(self.table_name, type(self.table_name))
        has_von_mises = False
        msg_temp, unused_nnodes, unused_is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1, has_von_mises)

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))
            #dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            #header[1] = dt_line
            #msg = header + msg_temp
            #f06_file.write('\n'.join(msg))

            if self.element_type == 144: # CQUAD4 bilinear
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 4, has_von_mises)
            elif self.element_type == 33:  # CQUAD4 linear
                #assert has_von_mises is False, has_von_mises
                self._write_f06_tri3_transient(f06_file, itime)
            elif self.element_type == 74: # CTRIA3
                #assert has_von_mises is False, has_von_mises
                self._write_f06_tri3_transient(f06_file, itime)
            elif self.element_type == 64:  #CQUAD8
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 8, has_von_mises)
            elif self.element_type == 82:  # CQUADR
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 4, has_von_mises)
            elif self.element_type == 70:  # CTRIAR
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 6, has_von_mises)
            elif self.element_type == 75:  # CTRIA6
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 6, has_von_mises)
            else:
                raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_f06_tri3_transient(self, f06_file, itime):
        """
        CQUAD4 linear
        CTRIA3
        """
        fds = self.fiber_distance
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        #nids = self.element_node[:, 1]

        ilayer0 = True
        for eid, fd, oxx, oyy, txy in zip(eids, fds, oxx, oyy, txy):
            sfd = write_float_13e(fd)
            soxx = write_float_13e(oxx)
            soyy = write_float_13e(oyy)
            stxy = write_float_13e(txy)

            if ilayer0:    # TODO: assuming 2 layers?
                f06_file.write('0  %-13s  %6s   %-13s  %-13s  %s\n' % (
                    eid, sfd, soxx, soyy, stxy))
            else:
                f06_file.write('   %-13s  %6s   %-13s  %-13s  %s\n' % (
                    '', sfd, soxx, soyy, stxy))
            ilayer0 = not ilayer0

    def _write_f06_quad4_bilinear_transient(self, f06_file, itime: int, nnodes: int, has_von_mises: bool) -> None:
        """
        CQUAD4 bilinear
        CQUAD8
        CTRIAR
        CTRIA6
        """
        fds = self.fiber_distance
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        ilayer0 = True
        for eid, nid, fd, doxx, doyy, dtxy in zip(eids, nids, fds, oxx, oyy, txy):
            sfd = write_float_13e(fd)
            soxx = write_float_13e(doxx)
            soyy = write_float_13e(doyy)
            stxy = write_float_13e(dtxy)

            if ilayer0:    # TODO: assuming 2 layers?
                f06_file.write('0  %-13s  %6s   %-13s  %-13s  %s\n' % (
                    eid, sfd, soxx, soyy, stxy))
            else:
                f06_file.write('   %-13s  %6s   %-13s  %-13s  %s\n' % (
                    '', sfd, soxx, soyy, stxy))
            ilayer0 = not ilayer0

class RandomPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        self.fiber_distance = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return get_nnodes(self)

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def build(self):
        r"""sizes the vectorized attributes of the RandomPlateArray

        SORT1:
         - etype     ndata numwide size  -> nelements     nnodes nlayers
         - CQUAD-144 376   47      4        376/(47*4)=2  5      2*2*5=20  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\plate_111.op2

        """
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']

        nnodes = self.nnodes_per_element

        #self.names = []
        self.nelements //= nnodes
        #print('element_type=%r ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, self.ntotal, self.subtitle))
        #self.nelements //= self.ntimes
        #self.ntotal = self.nelements # * 2
        #if self.element_name == 'CTRIA3':
        #print('element_type=%r ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.element_type, self.ntimes, self.nelements, self.ntotal, self.subtitle))
        #print(self)
        if self.is_sort1:
            # old
            #ntimes = self.ntimes
            #nelements = self.nelements
            ntimes = len(self._ntotals)
            ntotal = self._ntotals[0]
            nelements = ntotal // (2 * nnodes)

            #ntotal = self.ntotal
            #ny = nelements * 2
            nelements_nnodes = nelements * 2
            nlayers = nelements * 2 * nnodes
            if nlayers != ntotal:
                msg = f'nlayers={nlayers} ntotal={ntotal}\n'
                msg += f'SORT1 {self.element_name}-{self.element_type} ntimes={ntimes} nelements={nelements} ntotal={ntotal}\n'
                msg = f'_ntotals={self._ntotals}\n'
                raise RuntimeError(msg)
            #nx = ntimes
            #ny = nelements_nnodes
            #ntotal = nelements * 2
            #if self.element_name in ['CTRIA3', 'CQUAD8']:
            #print(f"SORT1 ntimes={ntimes} nelements={nelements} ntotal={ntotal}")
        elif self.is_sort2:
            # ntotal=164
            # len(_ntotals) = 4580 -> nelements=4580
            # nfreqs=82
            # flip this to sort1?
            #ntimes = self.ntotal
            #nnodes = self.ntimes
            #ntotal = nnodes
            #nelements = self.ntimes
            #ntimes = self.nelements # // nelements
            #ntotal = self.ntotal
            nelements = len(self._ntotals)
            ntotal = self._ntotals[0]
            ntimes = ntotal // 2 // nnodes

            #print(self._ntotals)
            ntotal = self._ntotals[0]
            #nelements = len(self._ntotals)
            #ntimes = ntotal // 2
            #ntimes, nelements = nelements_real, ntimes_real
            #ntotal = self.ntotal
            nlayers = nelements * 2 * nnodes
            #ny = nlayers
            #nx = ntimes
            #if self.element_name in ['CTRIA3', 'CQUAD8']:
            #if self.element_name in ['CQUAD4']:
                #print(f"SORT2 ntimes={ntimes} nelements={nelements} ntotal={ntotal} nnodes={nnodes} nlayers={nlayers}")
        else:  # pragma: no cover
            raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())

        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self.build_data(ntimes, nelements, nlayers, nnodes, self._times_dtype)
        #print(''.join(self.get_stats()))
        #try:
            #name = self.data_code['name']
        #except KeyError:
            #print(''.join(self.get_stats()))
            #raise

    def build_data(self, ntimes, nelements, nlayers, nnodes, dtype):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements
        #ntotal = nelements * 2
        self.ntotal = nlayers
        #_times = zeros(ntimes, dtype=self.analysis_fmt)
        #element = zeros(nelements, dtype='int32')

        self._times = zeros(ntimes, dtype=self.analysis_fmt)
        #self.ntotal = self.nelements * nnodes

        #print(f'***nelements={nelements} nlayers={nlayers} ntimes={ntimes}')
        self.element_node = zeros((nlayers, 2), 'int32')

        # the number is messed up because of the offset for the element's properties
        #if not self.nelements * 2 == self.ntotal:
            #msg = 'ntimes=%s nelements=%s nlayers=%s ntotal=%s' % (
                #self.ntimes, self.nelements, self.nelements * 2, self.ntotal)
            #raise RuntimeError(msg)

        self.fiber_distance = zeros(nlayers, 'float32')

        # [oxx, oyy, txy]
        #print(f'ntimes={self.ntimes} nelements={self.nelements} ntotal={self.ntotal}')
        self.data = zeros((ntimes, nlayers, 3), 'float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        #print(f'column_names = {column_names} column_values={column_values}')

        #print(self.element_node)
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\psdo7026.op2
        #self.data_frame = pd.Panel(self.data, items=column_values,
                                   #major_axis=self.element_node, minor_axis=headers).to_frame()
        #self.data_frame.columns.names = column_names
        #self.data_frame.index.names = ['ElementID', 'Item']
        #return

        names = ['ElementID', 'NodeID']
        #ipos = np.where(self.element_node[:, 0] > 0)
        ipos = None

        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_element_node(
            column_values, column_names,
            headers, self.element_node, self.data)
        #print(data_frame)
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        assert self.element_node[:, 0].min() > 0, self.element_node
        assert table.element_node[:, 0].min() > 0, table.element_node
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (
                self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid, Nid\n' % str(self.code_information())
            for eid1, eid2 in zip(self.element, table.element):
                msg += '(%s, %s), (%s, %s)\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)
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
                        (oxx1, oyy1, txy1) = t1
                        (oxx2, oyy2, txy2) = t2
                        #d = t1 - t2
                        if not np.allclose(
                                [oxx1, oyy1, txy1], # atol=0.0001
                                [oxx2, oyy2, txy2], atol=0.075):
                            ni = len(str(eid)) + len(str(eid))
                        #if not np.array_equal(t1, t2):
                            msg += ('%s  (%s, %s, %s)\n'
                                    '%s  (%s, %s, %s)\n' % (
                                        eid,
                                        oxx1, oyy1, txy1,
                                        ' ' * ni,
                                        oxx2, oyy2, txy2,
                                    ))
                            msg += ('%s  (%s, %s, %s)\n'
                                    % (
                                        ' ' * ni,
                                        oxx1 - oxx2,
                                        oyy1 - oyy2,
                                        txy1 - txy2,
                                    ))

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

    def finalize(self):
        """Calls any OP2 objects that need to do any post matrix calcs"""
        if self.is_sort1:
            return
        #print('finalize random plate')
        self.set_as_sort1()
        #print(self.get_stats())
        #print(self._times, self._times.dtype)
        #print(self.element_node)

    def add_sort1(self, dt, eid, nid,
                  fd1, oxx1, oyy1, txy1,
                  fd2, oxx2, oyy2, txy2):
        assert self.is_sort1, self.sort_method
        #assert self.element_node.max() == 0, self.element_node
        #if self.element_name in ['CTRIA3', 'CQUAD8']:
        #print(f'SORT1 {self.element_name}: itime={self.itime} ielement={self.ielement} itotal={self.itotal} dt={dt} eid={eid} nid={nid} fd={fd1} oxx={oxx1}')
        #print(f'SORT1 {self.element_name}: itime={self.itime} ielement={self.ielement} itotal={self.itotal+1} dt={dt} eid={eid} nid={nid} fd={fd2} oxx={oxx2}')
            #print('%s: itime=%s itotal=%s dt=%s eid=%s fd=%s oxx=%s' % (self.element_name, self.itime, self.itotal, dt, eid, fd, oxx))
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)

        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.data[self.itime, self.itotal, :] = [oxx1, oyy1, txy1]
        self.element_node[self.itotal, :] = [eid, nid]  # 0 is center
        self.fiber_distance[self.itotal] = fd1
        #self.ielement += 1
        self.itotal += 1

        self.data[self.itime, self.itotal, :] = [oxx2, oyy2, txy2]
        self.element_node[self.itotal, :] = [eid, nid]  # 0 is center
        self.fiber_distance[self.itotal] = fd2
        self.itotal += 1

    def _get_sort2_itime_ilower_iupper_from_itotal(self):
        nnodes = self.nnodes_per_element
        itime = self.ielement // nnodes
        inid = self.ielement % nnodes
        itotal = self.itotal
        #if itime >= self.data.shape[0]:# or itotal >= self.element_node.shape[0]:
            #print(f'*SORT2 {self.element_name}: itime={itime} ielement={self.itime} inid={inid} itotal={itotal} dt={dt} eid={eid} nid={nid} fd={fd1:.2f} oxx={oxx1:.2f}')
            #print(f'*SORT2 {self.element_name}: itime={itime} ielement={self.itime} inid={inid} itotal={itotal+1} dt={dt} eid={eid} nid={nid} fd={fd2:.2f} oxx={oxx2:.2f}')
            #print(self.data.shape)
            #print(self.element_node.shape)
        #else:
        ielement = self.itime

        #ibase = 2 * ielement # ctria3/cquad4-33
        ibase = 2 * (ielement * nnodes + inid)
        ie_upper = ibase
        ie_lower = ibase + 1

        #if self.element_name == 'CTRIAR': # and self.table_name == 'OESATO2':
        #debug = False
        #if self.element_name == 'CTRIA3' and self.table_name in ['OSTRRMS1', 'OSTRRMS2']:
            #debug = True
            #print(f'SORT2 {self.table_name} {self.element_name}: itime={itime} ie_upper={ie_upper} ielement={self.itime} inid={inid} nid={nid} itotal={itotal} dt={dt} eid={eid} nid={nid} fd={fd1:.2f} oxx={oxx1:.2f}')
            #print(f'SORT2 {self.table_name} {self.element_name}: itime={itime} ie_lower={ie_lower} ielement={self.itime} inid={inid} nid={nid} itotal={itotal+1} dt={dt} eid={eid} nid={nid} fd={fd2:.2f} oxx={oxx2:.2f}')
        return itime, ie_upper, ie_lower

    def add_sort2(self, dt, eid, nid,
                  fd1, oxx1, oyy1, txy1,
                  fd2, oxx2, oyy2, txy2):
        assert self.is_sort2, self
        #if self.element_name == 'CTRIA3':
        #assert self.element_node.max() == 0, self.element_node
        #print(self.element_node, nid)
        itime, ie_upper, ie_lower = self._get_sort2_itime_ilower_iupper_from_itotal()
        self._times[itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)

        #itime = self.ielement
        #itime = self.itime
        #ielement = self.itime
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)

        if itime == 0:
            self.element_node[ie_upper, :] = [eid, nid]  # 0 is center
            self.element_node[ie_lower, :] = [eid, nid]  # 0 is center
            self.fiber_distance[ie_upper] = fd1
            self.fiber_distance[ie_lower] = fd2
            #if self.element_name == 'CQUAD4':
                #print(self.element_node)

        self.data[itime, ie_upper, :] = [oxx1, oyy1, txy1]
        self.data[itime, ie_lower, :] = [oxx2, oyy2, txy2]

        self.itotal += 2
        self.ielement += 1
        #debug = False
        #if debug:
            #print(self.element_node)

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append(f'  type={self.class_name} ntimes={ntimes} nelements={nelements:d} nnodes={nnodes:d} table_name={self.table_name}\n')
        else:
            msg.append(f'  type={self.class_name} nelements={nelements:d} nnodes={nnodes:d} {self.table_name}\n')
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 3] where 3=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append(f'  element_node.shape = {self.element_node.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        #print(self.table_name, type(self.table_name))
        has_von_mises = False
        msg_temp, unused_nnodes, unused_is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1, has_von_mises)

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))
            #dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            #header[1] = dt_line
            #msg = header + msg_temp
            #f06_file.write('\n'.join(msg))

            if self.element_type == 144: # CQUAD4 bilinear
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 4, has_von_mises)
            elif self.element_type == 33:  # CQUAD4 linear
                #assert has_von_mises is False, has_von_mises
                self._write_f06_tri3_transient(f06_file, itime)
            elif self.element_type == 74: # CTRIA3
                #assert has_von_mises is False, has_von_mises
                self._write_f06_tri3_transient(f06_file, itime)
            elif self.element_type == 64:  #CQUAD8
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 8, has_von_mises)
            elif self.element_type == 82:  # CQUADR
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 4, has_von_mises)
            elif self.element_type == 70:  # CTRIAR
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 6, has_von_mises)
            elif self.element_type == 75:  # CTRIA6
                self._write_f06_quad4_bilinear_transient(f06_file, itime, 6, has_von_mises)
            else:
                raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_f06_tri3_transient(self, f06_file, itime):
        """
        CQUAD4 linear
        CTRIA3
        """
        fds = self.fiber_distance
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        #nids = self.element_node[:, 1]

        ilayer0 = True
        for eid, fd, oxx, oyy, txy in zip(eids, fds, oxx, oyy, txy):
            sfd = write_float_13e(fd)
            soxx = write_float_13e(oxx)
            soyy = write_float_13e(oyy)
            stxy = write_float_13e(txy)

            if ilayer0:    # TODO: assuming 2 layers?
                f06_file.write('0  %-13s  %6s   %-13s  %-13s  %s\n' % (
                    eid, sfd, soxx, soyy, stxy))
            else:
                f06_file.write('   %-13s  %6s   %-13s  %-13s  %s\n' % (
                    '', sfd, soxx, soyy, stxy))
            ilayer0 = not ilayer0

    def _write_f06_quad4_bilinear_transient(self, f06_file, itime: int, nnodes: int, has_von_mises: bool) -> None:
        """
        CQUAD4 bilinear
        CQUAD8
        CTRIAR
        CTRIA6
        """
        fds = self.fiber_distance
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        ilayer0 = True
        for eid, nid, fd, doxx, doyy, dtxy in zip(eids, nids, fds, oxx, oyy, txy):
            sfd = write_float_13e(fd)
            soxx = write_float_13e(doxx)
            soyy = write_float_13e(doyy)
            stxy = write_float_13e(dtxy)

            if ilayer0:    # TODO: assuming 2 layers?
                f06_file.write('0  %-13s  %6s   %-13s  %-13s  %s\n' % (
                    eid, sfd, soxx, soyy, stxy))
            else:
                f06_file.write('   %-13s  %6s   %-13s  %-13s  %s\n' % (
                    '', sfd, soxx, soyy, stxy))
            ilayer0 = not ilayer0

def _get_plate_msg(self, is_mag_phase=True, is_sort1=True, has_von_mises: bool=False):
    #if self.is_von_mises:
        #von_mises = 'VON MISES'
    #else:
        #von_mises = 'MAX SHEAR'

    if self.is_stress:
        if self.is_fiber_distance:
            vm_msg_temp = ['    ELEMENT              FIBER                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            msg_temp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                        '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            vm_msg_temp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            msg_temp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                        '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
    else:
        if self.is_fiber_distance:
            vm_msg_temp = ['    ELEMENT              FIBER                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            msg_temp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                        '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            vm_msg_temp = ['    ELEMENT              FIBRE                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            msg_temp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                        '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']


    #if is_mag_phase:
        #mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
    #else:
        #mag_real = ['                                                          (REAL/IMAGINARY)\n', ' \n']
    ## TODO: validation on header formatting...

    if self.is_stress:
        cquad4_bilinear = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n']
        cquad4_linear = ['                S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad8 = ['                S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria3 = ['                   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        ctria6 = ['                   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']
    else:
        cquad4_bilinear = ['                S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n']
        cquad4_linear = ['                S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']
        cquad8 = ['                S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria3 = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']
        ctria6 = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

    #msg = []
    is_bilinear = False
    if self.element_type == 144: # CQUAD4
        is_bilinear = True
        #msg += cquad4_linear + mag_real + grid_msg_temp
        msg0 = cquad4_linear
    elif self.element_type == 33: # CQUAD4
        is_bilinear = False
        #msg += cquad4_bilinear + mag_real + fiber_msg_temp
        msg0 = cquad4_bilinear
    elif self.element_type == 64:  #CQUAD8
        #msg += cquad8 + mag_real + grid_msg_temp
        msg0 = cquad8
        is_bilinear = True
    elif self.element_type == 82:  # CQUADR
        #msg += cquadr + mag_real + grid_msg_temp
        msg0 = cquadr
        is_bilinear = True

    elif self.element_type == 74: # CTRIA3
        #msg += ctria3 + mag_real + fiber_msg_temp
        msg0 = ctria3
    elif self.element_type == 75:  # CTRIA6
        #msg += ctria6 + mag_real + grid_msg_temp
        msg0 = ctria6
        is_bilinear = True
    elif self.element_type == 70:  # CTRIAR
        #msg += ctriar + mag_real + grid_msg_temp
        msg0 = ctriar
        is_bilinear = True
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
    msgi = msg0

    nnodes = get_nnodes(self)
    # TODO: STRESSES IN ELEMENT  COORDINATE SYSTEM???
    if self.table_name in ['OESPSD1', 'OSTRPSD1']:
        #assert self.is_fiber_distance, self.table_name
        msg = msgi + [
            #'                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
            '                                                ( POWER SPECTRAL DENSITY FUNCTION )\n'
            ' \n'
            #'                    FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n'
            #'    FREQUENCY      DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n'
            #'0  2.000000E+01 -5.000000E-02              1.925767E-05                      1.404795E-04                     1.097896E-03'
            #'                 5.000000E-02              1.925766E-05                      1.404794E-04                     1.097896E-03'
        ] + vm_msg_temp
    elif self.table_name in ['OESNO1', 'OSTRNO1']:
        #assert self.is_fiber_distance is False, self.table_name
        msg = msgi + [
            #'                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S     ( Q U A D 4 )        OPTION = BILIN  '
            '                                                ( NUMBER OF ZERO CROSSINGS )\n'
            ' \n'
            #'                      STRAIN                                  -  STRAINS IN ELEMENT  COORDINATE SYSTEM -\n'
            #'   ELEMENT-ID  GRID ID   CURVATURE               NORMAL-X                       NORMAL-Y                          SHEAR-XY\n'
            #'0           2  CEN/4   0.0                     5.492230E+02                   5.786333E+02                      5.789550E+02'
            #'                      -1.000000E+00            6.206278E+01                   5.312986E+01                      5.217435E+01'
        ] + vm_msg_temp
    elif self.table_name in ['OESRMS1', 'OSTRRMS1']:  # OESRMS1
        #assert self.is_fiber_distance is False, self.table_name
        msg = msgi + [
            #'                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S     ( Q U A D 4 )        OPTION = BILIN  '
            '                                    ( ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  X.XXE+XX )\n'
            ' \n'
            #'                      STRAIN                                  -  STRAINS IN ELEMENT  COORDINATE SYSTEM -\n'
            #'   ELEMENT-ID  GRID ID   CURVATURE               NORMAL-X                       NORMAL-Y                          SHEAR-XY\n'
            #'0           2  CEN/4   0.0                     1.728666E-16                   3.661332E-16                      3.067126E-15'
            #'                      -1.000000E+00            6.043902E-01                   6.198736E-02                      1.028200E-01'
        ] + vm_msg_temp
    elif self.table_name in ['OSTRCRM1', 'OESCRM1']:
        msg = msgi + [
            '                               ( CUMULATIVE ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  X.XXE+XX )\n'
            ' \n'
            ] + msg_temp
    elif self.table_name in ['OESATO1', 'OSTRATO1']:
        msg = msgi + [
            '                                    ( ATO ? )\n'
            ' \n'
            ] + msg_temp
    else:
        raise NotImplementedError(self.table_name)

    return msg, nnodes, is_bilinear

def get_nnodes(self) -> int:
    if self.element_type in [64, 82, 144]:  # ???, CQUADR, CQUAD4 bilinear
        nnodes = 4 + 1 # centroid
    elif self.element_type in [70, 75]:   #???, CTRIA6
        nnodes = 3 + 1 # centroid
    elif self.element_type in [144, 74, 33]:  # CTRIA3, CQUAD4 linear
        nnodes = 1
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
    return nnodes

class RandomPlateStressArray(RandomPlateArray, StressObject):
    """
    NX
    2 FD1  RS Z1 = Fibre Distance
    3 SX1  RS Normal in x at Z1
    4 SY1  RS Normal in y at Z1
    5 TXY1 RS Shear in xy at Z1
    6 FD2  RS Z2 = Fibre Distance
    7 SX2  RS Normal in x at Z2
    8 SY2  RS Normal in y at Z2
    9 TXY2 RS Shear in xy at Z2
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        headers = ['oxx', 'oyy', 'txy']
        return headers

    def get_headers(self) -> list[str]:
        return self._get_headers()


class RandomPlateStrainArray(RandomPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        headers = ['exx', 'eyy', 'exy']
        return headers

    def get_headers(self) -> list[str]:
        return self._get_headers()

class RandomPlateStressVMArray(RandomPlateVMArray, StressObject):
    """
    NX
    2 FD1  RS Z1 = Fibre Distance
    3 SX1  RS Normal in x at Z1
    4 SY1  RS Normal in y at Z1
    5 TXY1 RS Shear in xy at Z1
    6 FD2  RS Z2 = Fibre Distance
    7 SX2  RS Normal in x at Z2
    8 SY2  RS Normal in y at Z2
    9 TXY2 RS Shear in xy at Z2
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomPlateVMArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        headers = ['oxx', 'oyy', 'txy']
        return headers

    def get_headers(self) -> list[str]:
        return self._get_headers()


class RandomPlateStrainVMArray(RandomPlateVMArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomPlateVMArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        headers = ['exx', 'eyy', 'exy']
        return headers

    def get_headers(self) -> list[str]:
        return self._get_headers()
