import warnings
from typing import List
import numpy as np

from pyNastran.op2.result_objects.op2_objects import ScalarObject, get_times_dtype
from pyNastran.f06.f06_formatting import (
    write_floats_10e, _eigenvalue_header)


class GridPointSurfaceArray(ScalarObject):
    """
    '                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
    '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
    '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
    '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
    '0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
    '      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
    '      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'

    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        self.node_element = None
        self.location = None
        self._times = None

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self) -> bool:
        return True
    @property
    def is_complex(self) -> bool:
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        if self.is_built:
            return
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self.node_element = np.zeros((self.ntotal, 2), dtype=idtype)
        #oxx, oyy, txy, angle, major, minor, ovm
        self.data = np.zeros((self.ntimes, self.nelements, 8), dtype=fdtype)
        self.location = np.empty(self.ntotal, dtype='U8')

        self._times = np.zeros(self.ntimes, dtype=dtype)

    def _write_table_3(self, op2_file, op2_ascii, new_result, itable, itime): #, itable=-3, itime=0):
        import inspect
        from struct import pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        #if itable == -3:
        #print('*writing itable=%s' % itable)
        if new_result and itable != -3:
            header = [
                4, 146, 4,
            ]
        else:
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
                4, 146, 4,
            ]
        op2_file.write(pack(b'%ii' % len(header), *header))
        op2_ascii.write('table_3_header = %s\n' % header)
        #op2_file.write(pack('12i', *header))
        #else:
            #print('***writing itable=%s' % itable)
            #op2_file.write(pack('3i', *[
                ##4, itable, 4,
                ##4, 1, 4,
                ##4, 0, 4,
                #4, 146, 4,
            #]))
        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        #[
            #'aCode', 'tCode', 'element_type', 'isubcase',
            #'???', '???', '???', 'load_set'
            #'format_code', 'num_wide', 's_code', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', 'Title', 'subtitle', 'label']
        #random_code = self.random_code
        ogs = self.ogs
        if ogs is None:
            #print(''.join(self.get_stats()))
            warnings.warn('ogs=0...')
            ogs = 0

        format_code = self.format_code
        s_code = self.sCode
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'50i 128s 128s 128s'
        unused_oCode = 0

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        field6 = 0
        field7 = 0
        if self.analysis_code == 1:
            field5 = self.lsdvmns[itime]
            if np.isnan(field5):  # poor sort2 -> sort1
                raise RuntimeError('field5 in a static case is nan...; do you have SORT2?')
                #field5 = 1

        elif self.analysis_code == 2:
            field5 = self.modes[itime]
            field6 = self.eigns[itime]
            field7 = self.cycles[itime]
            assert isinstance(field6, float), type(field6)
            assert isinstance(field7, float), type(field7)
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7

        #elif self.analysis_code == 3:
            #field5 = self.freqs[itime]
        elif self.analysis_code == 5:
            field5 = self.freqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            field5 = self.dts[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 7:  # pre-buckling
            field5 = self.lsdvmns[itime] # load set number
        elif self.analysis_code == 8:  # post-buckling
            field5 = self.lsdvmns[itime] # load set number
            #if hasattr(self, 'eigns'):
            if hasattr(self, 'eigens'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
        elif self.analysis_code == 9:  # complex eigenvalues
            field5 = self.modes[itime]
            if hasattr(self, 'eigns'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                print(self.get_stats())
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=9')

            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            field7 = self.eigis[itime]
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.lftsfqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            field5 = self.lsdvmns[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)


        #self.ogs = self.add_data_parameter(data, 'ogs_id', b'i', 3, False)
        #self.refid = self.add_data_parameter(data, 'refid', b'i', 8, False)
        #self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)
        #self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)
        #self.sCode = self.add_data_parameter(data, 'sCode', b'i', 11, False)
        #self.oCoord = self.add_data_parameter(data, 'oCoord', b'i', 12, False)
        #self.axis = self.add_data_parameter(data, 'axis', b'i', 13, False)
        #self.normal = self.add_data_parameter(data, 'normal', b'i', 14, False)

        table3 = [
            approach_code, table_code, ogs, isubcase, field5,
            field6, field7, self.refid, format_code, num_wide,
            s_code, self.oCoord, self.axis, self.normal, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert table3[22] == thermal

        table3 = fix_table3_types(table3, size=4)
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2_file.write(pack(fmt, *data))

    #def build_dataframe(self):
        #"""creates a pandas dataframe"""
        #import pandas as pd
        #headers = self.get_headers()
        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        #if self.nonlinear_factor not in (None, np.nan):
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
        #else:
            #self.data_frame = pd.Panel(self.data, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
        #self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']

    def add_sort1(self, dt, nid, eid, fiber, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        """unvectorized method for adding SORT1 transient data"""
        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.node_element[self.itotal, :] = [nid, eid]
        self.location[self.itotal] = fiber
        self.data[self.itime, self.itotal, :] = [nx, ny, txy, angle, majorP, minorP, tmax, ovm]
        self.itotal += 1

    def get_stats(self, short: bool=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

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
        msg.append(f'  node_element.shape = {self.node_element.shape}\n')
        msg.append(f'  location.shape = {self.location.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg += self.get_data_code()
        return msg


    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []

        cid = self.refid
        axis_int = self.oCoord
        axis_map = {0 : 'X', 1 : 'Y', 2 : 'Z'}
        axis = axis_map[axis_int]
        msg = self._get_f06_message(self.ogs_id, cid, axis)

        ntimes = self.data.shape[0]

        nids = self.node_element[:, 0]
        eids = self.node_element[:, 1]
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            nx = self.data[itime, :, 0]
            ny = self.data[itime, :, 1]
            txy = self.data[itime, :, 2]
            angle = self.data[itime, :, 3]
            majorp = self.data[itime, :, 4]
            minorp = self.data[itime, :, 5]
            tmax = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]
            fibers = self.location
            nid_old = -1
            for (nid, eid, fiber, nxi, nyi, txyi, anglei, majorpi, minorpi, tmaxi, ovmi) in zip(
                nids, eids, fibers, nx, ny, txy, angle, majorp, minorp, tmax, ovm):
                [nxi, nyi, txyi, majorpi, minorpi, tmaxi, ovmi] = write_floats_10e([
                    nxi, nyi, txyi, majorpi, minorpi, tmaxi, ovmi])
                if nid > nid_old:
                    f06_file.write(
                        '0%8s  %8s   %4s    %-10s %-10s %-10s  %8.4f %10s %10s %10s  %s\n' % (
                            nid, eid, fiber, nxi, nyi, txyi, anglei, majorpi, minorpi,
                            tmaxi, ovmi))
                else:
                    f06_file.write(
                        ' %8s  %8s   %4s    %-10s %-10s %-10s  %8.4f %10s %10s %10s  %s\n' % (
                            '', '', fiber, nxi, nyi, txyi, anglei, majorpi, minorpi,
                            tmaxi, ovmi))
                nid_old = nid
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            #print('***************', itable)
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        #eids2 = self.element_node[:, 0]
        #nodes = self.element_node[:, 1]
        #nelements_nodes = len(nodes)

        #eids3 = self.element_cid[:, 0]
        #cids3 = self.element_cid[:, 1]

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        #nelements = len(np.unique(eids2))

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)
        #nnodes_centroid = self.nnodes_per_element
        #nnodes_no_centroid = self.nnodes_per_element_no_centroid
        nnodes = self.data.shape[1]
        #ntotali = 11
        ntotali = self.num_wide
        assert ntotali == 11, ntotali
        ntotal = ntotali * nnodes


        #print('shape = %s' % str(self.data.shape))
        #assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')
        ntimes = self.ntimes

        #print('ntotal=%s' % (ntotal))
        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #op2_format = endian + b'2i6f'

        #idtype = self.element_cid.dtype
        fdtype = self.data.dtype
        #print(self.size)
        if self.size == 4:
            grid_bytes = b'GRID'
        else:
            warnings.warn(f'downcasting {self.class_name}...')
            idtype = np.int32(1)
            fdtype = np.float32(1.0)
            grid_bytes = b'GRID'

        #[nids, eids, fibers, nx, ny, txy, angle, majorp, minorp, tmax, ovm]
        nids = self.node_element[:, 0]
        eids = self.node_element[:, 1]
        nids_device = nids * 10 + self.device_code

        nids_device

        # speed up transient cases, but slightly slows down static cases
        data_out = np.empty((nnodes, 11), dtype=fdtype)
        # setting:
        #  - [nid_device, eids, location_bytes]
        data_out[:, 0] = nids_device
        data_out[:, 1] = eids
        location_bytes = np.array([loc.encode('ascii') for loc in self.location])
        data_out[:, 2] = location_bytes.view(fdtype)


        #nx = self.data[itime, :, 0]
        #ny = self.data[itime, :, 1]
        #txy = self.data[itime, :, 2]
        #angle = self.data[itime, :, 3]
        #majorp = self.data[itime, :, 4]
        #minorp = self.data[itime, :, 5]
        #tmax = self.data[itime, :, 6]
        #ovm = self.data[itime, :, 7]
        #fibers = self.location

        #cen_array = np.full(nelements, grid_bytes, dtype='|S4')
        #nnodes_no_centroid_array = np.full(nelements, nnodes_no_centroid, dtype=idtype)

        #element_wise_data = to_column_bytes([
            #element_device, # ints
            #cids3, # ints
            #cen_array, # bytes
            #nnodes_no_centroid_array, # ints
        #], fdtype, debug=False)

        # we could tack the nodes on, so we don't have to keep stacking it
        # but we run into issues with datai
        #
        # total=nelements_nodes
        #nodes_view = nodes.view(fdtype).reshape(nelements, nnodes_centroid)
        #inode = np.arange(nnodes_centroid)
        #data_out[:, 4+inode*21] = nodes_view[:, inode]

        op2_ascii.write(f'nnodes={nnodes:d}\n')
        struct_i = Struct('i')
        struct_13i = Struct('13i')
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(struct_13i.pack(*header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')


            # stack each output by columns and fix any dtypes
            #datai2 = datai.reshape(nelements, 21*nnodes_centroid)
            #data_out = np.hstack([element_wise_data, datai2])
            #data_out[:, 4:] = datai2

            # switch datai to element format and put it in the output buffer
            data_out[:, 3:] = self.data[itime, :, :]
            assert data_out.size == ntotal
            op2_file.write(data_out)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(struct_i.pack(*header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for inid, (nid, eid) in enumerate(self.node_element):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (nx1, ny1, txy1, majorp1, minorp1, tmax1, ovm1) = t1
                        (nx2, ny2, txy2, majorp2, minorp2, tmax2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s %s\n  (%s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                nid, eid,
                                nx1, ny1, txy1, majorp1, minorp1, tmax1, ovm1,
                                nx2, ny2, txy2, majorp2, minorp2, tmax2, ovm2)
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

    def _get_f06_message(self, ogs_id: int, cid: int, axis: str) -> List[str]:
        raise NotImplementedError()

class GridPointSurfaceStressesArray(GridPointSurfaceArray):

    def get_headers(self) -> List[str]:
        headers = ['nx', 'ny', 'txy', 'angle', 'majorP', 'minorP', 'tmax', 'ovm']
        return headers

    def _get_f06_message(self, ogs_id: int, cid: int, axis: str) -> List[str]:
        msg = [
            f'                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E    {ogs_id:d}\n',
            f'0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  {axis}         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        {cid}\n',
            '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
            '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
           #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
           #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
           #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        return msg



class GridPointSurfaceStrainsArray(GridPointSurfaceArray):

    def get_headers(self) -> List[str]:
        headers = ['nx', 'ny', 'exy', 'angle', 'majorP', 'minorP', 'emax', 'evm']
        return headers

    def _get_f06_message(self, ogs_id: int, cid: int, axis: str) -> List[str]:
        msg = [
            f'                                    S T R A I N S   A T   G R I D   P O I N T S   - -     S U R F A C E    {ogs_id:d}\n',
           #f'                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E    {ogs_id:d}\n',
            f'0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  {axis}         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        {cid}\n',
           #'     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
            '     GRID      ELEMENT            STRAINS  IN SURFACE SYSTEM           PRINCIPAL STRAINS             MAX             \n',
            '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
           #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
           #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
           #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        return msg




class GridPointStressesVolumePrincipalArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        self._times = None

    def get_headers(self) -> List[str]:
        headers = [
            'lxa', 'lxb', 'lxc',
            'lya', 'lyb', 'lyc',
            'lza', 'lzb', 'lzc',
            'sa', 'sb', 'sc',
            'epr', 'ovm']
        return headers

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for inid, nid in enumerate(self.node):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (lxa1, lxb1, lxc1, lya1, lyb1, lyc1, lza1, lzb1, lzc1, sa1, sb1, sc1, epr1, ovm1) = t1
                        (lxa2, lxb2, lxc2, lya2, lyb2, lyc2, lza2, lzb2, lzc2, sa2, sb2, sc2, epr2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                                nid,
                                lxa1, lxb1, lxc1, lya1, lyb1, lyc1, lza1,
                                lxa2, lxb2, lxc2, lya2, lyb2, lyc2, lza2)
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

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self) -> bool:
        return True
    @property
    def is_complex(self) -> bool:
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #print('self.IDs', self.data)
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self.node = np.zeros(self.ntotal, dtype=idtype)
        #lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm
        self.data = np.zeros((self.ntimes, self.ntotal, 14), dtype=fdtype)
        self.location = np.empty(self.ntotal, dtype='U8')

        self._times = np.zeros(self.ntimes, dtype=dtype)

    def get_stats(self, short: bool=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

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
        msg.append(f'  node.shape = {self.node.shape}\n')
        msg.append(f'  location.shape = {self.location.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg += self.get_data_code()
        return msg

    def add_sort1(self, dt, nid, lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm):
        assert isinstance(nid, int) and nid > 0, 'dt=%s nid=%s' % (dt, nid)
        self._times[self.itime] = dt
        self.node[self.itotal] = nid
        self.data[self.itime, self.itotal, :] = [lxa, lxb, lxc, lya, lyb, lyc, lza, lzb, lzc, sa, sb, sc, epr, ovm]
        self.itotal += 1

    #def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  #page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        #pass


class GridPointStressesVolumeDirectArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        self._times = None

    def get_headers(self) -> List[str]:
        headers = ['ox', 'oy', 'oz', 'txy', 'tyz', 'txz', 'pressure', 'ovm']
        return headers

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self) -> bool:
        return True
    @property
    def is_complex(self) -> bool:
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #print('self.IDs', self.data)
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        self.nelements //= self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self.node = np.zeros(self.ntotal, dtype=idtype)
        #oxx, oyy, txy, angle, major, minor, ovm
        self.data = np.zeros((self.ntimes, self.ntotal, 8), dtype=fdtype)
        self.location = np.empty(self.ntotal, dtype='U8')
        self._times = np.zeros(self.ntimes, dtype=dtype)

    def get_stats(self, short: bool=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

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
        msg.append(f'  node.shape = {self.node.shape}\n')
        msg.append(f'  location.shape = {self.location.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg += self.get_data_code()
        return msg

    def add_sort1(self, dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm):
        assert isinstance(nid, int) and nid > 0, 'dt=%s nid=%s' % (dt, nid)
        self._times[self.itime] = dt
        self.node[self.itotal] = nid
        self.data[self.itime, self.itotal, :] = [nx, ny, nz, txy, tyz, txz, pressure, ovm]
        self.itotal += 1

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        """
        '    D I R E C T   S T R E S S E S   A T   G R I D   P O I N T S   - -       V O L U M E      101'
        '        OUTPUT COORDINATE SYSTEM =       0  BASIC   '
        '    GRID            NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX        MEAN      VON MISES'
        '    ID                                                                                           PRESSURE'
        '        1           1.455E+03  -1.548E+02  -2.927E+02    -1.573E+01   3.326E+01  -3.438E+03     -3.357E+02   6.188E+03'
        '        2           1.093E+03  -1.996E+02  -1.682E+02     1.542E+02   5.962E+01  -4.104E+03     -2.417E+02   7.227E+03'
        """
        if header is None:
            header = []


        cid = self.refid
        #axis_int = self.oCoord
        #axis_map = {0 : 'X', 1 : 'Y', 2 : 'Z'}
        #axis = axis_map[axis_int]
        msg = [
            '                    D I R E C T   S T R E S S E S   A T   G R I D   P O I N T S   - -       V O L U M E      %3i\n'
            '                              OUTPUT COORDINATE SYSTEM = %7i  ELEMENT \n'
            '     GRID            NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX        MEAN      VON MISES\n'
            '     ID                                                                                           PRESSURE\n' % (
            #'     8086           6.136E-02   2.131E-01   8.353E-02    -2.268E+00  -2.274E-13   1.525E-13     -1.193E-01   3.930E+00'
            self.ogs_id, cid)
        ]

        ntimes = self.data.shape[0]

        nids = self.node
        zero = ' '
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            nx = self.data[itime, :, 0]
            ny = self.data[itime, :, 1]
            nz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]
            pressure = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]
            for (nid, nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi) in zip(
                nids, nx, ny, nz, txy, tyz, txz, pressure, ovm):
                [nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi] = write_floats_10e([
                    nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi])

                f06_file.write('%s%8s          %-10s  %-10s  %-10s    %-10s  %-10s  %-10s     %-10s  %-s\n' % (
                    zero, nid, nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi.rstrip()))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for inid, nid in enumerate(self.node):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (nx1, ny1, nz1, txy1, tyz1, txz1, pressure1, ovm1) = t1
                        (nx2, ny2, nz2, txy2, tyz2, txz2, pressure2, ovm2) = t2
                        if not np.allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                nid,
                                nx1, ny1, nz1, txy1, tyz1, txz1, pressure1, ovm1,
                                nx2, ny2, nz2, txy2, tyz2, txz2, pressure2, ovm2)
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

#msg = [
    #' P R I N C I P A L   G R I D   P O I N T   S T R E S S   D I S C O N T I N U I T I E S  - -       V O L U M E       %s\n'
    #'                              OUTPUT COORDINATE SYSTEM = %7i  ELEMENT \n'
    #'                              GRID       PRINCIPAL STRESS DISCONTINUITY     MEAN      VON MISES      ERROR\n'
    #'                              ID           A          B          C          PRESSURE                 EST.\n' % (
        #ivolume, cid)
    #'                              8086         5.448E-09  9.886E-08  2.026E-15  2.484E-09  1.086E-07   5.716E-08'
#]
# not sure what result this is for
#zero = '                              '
#f06_file.write('%s%8s  %-10s %-10s %-10s   %-10s %-10s %-10s %-10s  %-s\n' % (
    #zero, nid, nxi, nyi, nzi, txyi, tyzi, txzi, pressurei, ovmi.rstrip()))

GridPointStressesVolumeDiscontinutiesArray = None # tCode=34

class GridPointStressesSurfaceDiscontinutiesArray(ScalarObject): # tCode=35
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.ntotal = 0
        self.ntimes = 0
        self.nelements = 0
        self.itotal = 0
        self.ielement = 0
        self.data = None
        self.itime = None
        #self.node_element = None
        self._times = None

    def get_headers(self) -> List[str]:
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'pressure']
        return headers

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    @property
    def is_real(self) -> bool:
        return True
    @property
    def is_complex(self) -> bool:
        return False

    def build(self):
        """sizes the vectorized attributes of the GridPointStressesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #print('self.IDs', self.data)
        self.itime = 0
        self.ielement = 0
        self.itotal = 0

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes

        self.node = np.zeros(self.ntotal, dtype='int32')
        #oxx, oyy, ozz, txy, pressure
        self.data = np.zeros((self.ntimes, self.ntotal, 5), dtype='float32')
        self.location = np.empty(self.ntotal, dtype='U8')
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)

        self._times = np.zeros(self.ntimes, dtype=dtype)

    def get_stats(self, short: bool=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

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
        msg.append(f'  node.shape = {self.node.shape}\n')
        msg.append(f'  location.shape = {self.location.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg += self.get_data_code()
        return msg

    def add_sort1(self, dt, nid, oxx, oyy, ozz, txy, pressure):
        assert isinstance(nid, int) and nid > 0, 'dt=%s nid=%s' % (dt, nid)
        self._times[self.itime] = dt
        self.node[self.itotal] = nid
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, pressure]
        self.itotal += 1

class GridPointStrainsVolumePrincipalArray(GridPointStressesVolumePrincipalArray):
    pass

class GridPointStrainsVolumeDirectArray(GridPointStressesVolumeDirectArray):
    pass

GridPointStrainsVolumeDiscontinutiesArray = None

class GridPointStrainsSurfaceDiscontinutiesArray(GridPointStressesSurfaceDiscontinutiesArray):
    pass
