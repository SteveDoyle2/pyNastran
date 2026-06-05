# pylint: disable=E1101
from typing import TextIO, BinaryIO
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_13e_long)
import numpy as np
from pyNastran.op2.result_objects.scalar6_table_object import RealScalarTableArray

class RealTemperatureArray(RealScalarTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_csv(self, csv_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        Displacement Table
        ------------------
        Flag, SubcaseID,  iTime, NID,        T,  dy,  dz,  rx,  ry,  rz,  cd,  PointType
        1,            1,  0,     101, 0.014159, nan, nan, nan, nan, nan,   0,  1
        uses cd=-1 for unknown cd

        """
        name = str(self.__class__.__name__)
        if write_header:
            csv_file.write('%s\n' % name)
            headers = ['Flag', 'Subcase', 'iTime', 'Node', 'Temperature', 'dy', 'dz', 'rx', 'ry', 'rz', 'cd', 'PointType']
            csv_file.write('# ' + ','.join(headers) + '\n')
        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]

        #unused_times = self._times
        isubcase = self.isubcase

        flag_map = {
            # 'RealDisplacementArray': 1,
            # 'RealVelocityArray': 2,
            # 'RealAccelerationArray': 3,
            # 'RealEigenvectorArray': 4,

            # 'RealLoadVectorArray': 5,
            #'RealAppliedLoadsArray': 5,
            # 'RealSPCForcesArray': 6,
            # 'RealMPCForcesArray': 7,
            'RealTemperatureArray' : 8,

            #'HeatFlux' : 9,
            # 'RealTemperatureGradientAndFluxArray': 9,
        }
        flag = flag_map[name]
        # sort1 as sort1
        assert is_sort1 is True, is_sort1
        nid_len = '%d' % len(str(node.max()))
        cd = -1
        for itime in range(self.ntimes):
            #dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            for node_id, gridtypei, t1i in zip(node, gridtype, t1):
                if is_exponent_format:
                    vals2 = write_floats_13e_long([t1i])
                    (t1i, ) = vals2

                csv_file.write(f'{flag}, {isubcase}, {itime}, {node_id:{nid_len}d}, '
                               f'{t1i}, nan, nan, nan, nan, nan, -1, {gridtypei}\n')
        return

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient(header, page_stamp, page_num, f06_file,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words=False)

    def _write_f06_transient(self, header, page_stamp, page_num=1, f06_file=None,
                             is_mag_phase=False, is_sort1=True):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words=False)
