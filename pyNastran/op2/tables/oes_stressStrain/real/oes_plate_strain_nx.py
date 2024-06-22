# coding: utf-8
#pylint disable=C0103
from itertools import count
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_13e_long, _eigenvalue_header)
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object


class RealCPLSTRNPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if not is_sort1:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def is_bilinear(self):
        #if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            #return False
        #elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            #return True
        return False

    def build(self):
        """sizes the vectorized attributes of the RealCPLSTRNPlateArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 271:  # CPLSTN3
            nnodes_per_element = 1
        elif self.element_type == 272:  # CPLSTN4
            nnodes_per_element = 5
        elif self.element_type == 273:  # CPLSTN6
            nnodes_per_element = 4 # 7
        elif self.element_type == 274:  # CPLSTN8
            nnodes_per_element = 5 # 9
        else:  # pragma: no cover
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

        self.nnodes = nnodes_per_element
        self.nelements //= nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element,
            #self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)
        self.element = np.zeros(self.nelements, dtype=idtype)
        self.element_node = np.zeros((self.ntotal, 2), dtype=idtype)

        #[oxx, oyy, ozz, txy, ovm]
        self.data = np.zeros((self.ntimes, self.ntotal, 5), dtype=fdtype)

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, (eid, nid) in enumerate(self.element_node):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (oxx1, oyy1, ozz1, txy1, ovm1) = t1
                    (oxx2, oyy2, ozz2, txy2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            oxx1, ozz1, oyy1, txy1, ovm1,
                            oxx2, ozz2, oyy2, txy2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.nnodes
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // self.nnodes // 2

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, nelements, nnodes, nlayers, ntotal)
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = np.ravel([np.searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ##ind = np.searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def add_new_eid_sort1(self, dt, eid, node_id,
                          oxx, oyy, ozz, txy, ovm):
        assert isinstance(eid, integer_types), eid
        assert isinstance(node_id, integer_types), node_id
        assert eid > 0, f'eid={eid}, nid={node_id}'
        self._times[self.itime] = dt
        # assert self.itotal == 0, oxx
        self.element[self.ielement] = eid
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, ovm]
        self.itotal += 1
        self.ielement += 1

    def add_sort1(self, dt, eid, node_id,
                  oxx, oyy, ozz, txy, ovm):
        assert self.sort_method == 1, self
        assert eid is not None, eid
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        assert isinstance(node_id, integer_types), node_id
        assert eid > 0, f'eid={eid}, nid={node_id}'
        #print(self.element_node.shape, self.itotal, eid, node_id)
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, ovm]
        self.itotal += 1
        #self.ielement += 1

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg, nnodes, cen = _get_clpstsn_msg(self)

        ## write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        #'    ELEMENT'
        #'      ID      GRID-ID       NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES'
        #'0        72    CEN 0     -2.420339E-02      0.0               1.499626E-02     -1.436477E-02      2.429825E-02'
        #'                  53     -3.316031E-02      0.0               2.198507E-02     -1.747013E-02      3.360480E-02'
        #'                  39     -2.398538E-02      0.0               1.464934E-02     -8.596349E-03      2.306218E-02'
        #'                  38     -1.546447E-02      0.0               8.354364E-03     -1.702782E-02      1.706980E-02'
        #cen_word = 'CEN/%i' % nnodes
        #cen_word = cen
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[oxx, oyy, ozz, shear, ovm]
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            ovm = self.data[itime, :, 4]

            for (i, eid, nid, oxxi, oyyi, ozzi, txyi, ovmi) in zip(
                 count(), eids, nids, oxx, oyy, ozz, txy, ovm):
                [oxxi, oyyi, ozzi, txyi, ovmi] = write_floats_13e(
                [oxxi, oyyi, ozzi, txyi, ovmi])

                if self.element_type in [271, 272, 273, 274]:  # CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8
                    if nid == 0:
                        #'    ELEMENT'
                        #'      ID      GRID-ID       NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES'
                        #'0        72    CEN 0     -2.420339E-02      0.0               1.499626E-02     -1.436477E-02      2.429825E-02')
                        #'                  53     -3.316031E-02      0.0               2.198507E-02     -1.747013E-02      3.360480E-02'
                        #'                  39     -2.398538E-02      0.0               1.464934E-02     -8.596349E-03      2.306218E-02'
                        #'                  38     -1.546447E-02      0.0               8.354364E-03     -1.702782E-02      1.706980E-02'
                        f06_file.write(f'0  {eid}    CEN 0     {oxxi}     {oyyi}     {ozzi}     {txyi}     {ovmi}\n')
                    else:
                        f06_file.write(f'            {nid}     {oxxi}     {oyyi}     {ozzi}     {txyi}     {ovmi}\n')
                else:
                    raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

class RealCPLSTRNPlateStressNXArray(RealCPLSTRNPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCPLSTRNPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'von_mises']
        return headers


class RealCPLSTRNPlateStrainNXArray(RealCPLSTRNPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCPLSTRNPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['exx', 'eyy', 'ezz', 'exy', 'von_mises']
        return headers


def _get_clpstsn_msg(self):
    if self.is_von_mises:
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    assert von_mises == 'VON MISES', von_mises
    quad_msg_temp = [
        '    ELEMENT\n',
        '      ID      GRID-ID       NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES\n']
    tri_msg_temp = [
        '    ELEMENT\n',
        '      ID           NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES\n']
    if self.is_stress:
        cplstn3_msg = ['              S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T N 3 )\n'] + tri_msg_temp
        cplstn6_msg = ['              S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T N 6 )\n'] + tri_msg_temp
        cplstn8_msg = ['            S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( C P L S T N 8 )\n'] + quad_msg_temp
        cplstn4_msg = ['            S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( C P L S T N 4 )\n'] + quad_msg_temp
    else:
        # '              S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( C P L S T N 8 )'
        # '    ELEMENT'
        # '      ID      GRID-ID       NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES'
        # '0        41    CEN 0     -5.869190E-03      0.0               3.099033E-03      2.562248E-03      5.463578E-03'
        # '                  42     -5.577235E-03      0.0               2.865523E-03      8.868802E-04      4.983902E-03'
        # '                  40     -7.421685E-03      0.0               5.138225E-03      2.677066E-03      7.453323E-03'
        # '                  41     -4.968971E-03      0.0               2.495173E-03      5.169902E-03      5.306638E-03'
        # '                  43     -5.148691E-03      0.0               2.261392E-03      1.509432E-03      4.470889E-03'
        # '                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T N 3 )'
        # '    ELEMENT'
        # '      ID           NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES'
        # '0        36     -1.018986E-02      0.0               1.542302E-03     -7.523935E-04      7.374198E-03'

        cplstn3_msg = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T N 3 )\n'] + tri_msg_temp
        cplstn6_msg = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T N 6 )\n'] + tri_msg_temp
        cplstn8_msg = ['              S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( C P L S T N 8 )\n'] + quad_msg_temp
        cplstn4_msg = ['              S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( C P L S T N 4 )\n'] + quad_msg_temp


    #is_bilinear = False
    if self.element_type == 271:
        msg = cplstn3_msg
        nnodes = 3
        cen = 'CEN/3'
    elif self.element_type == 272:
        msg = cplstn4_msg
        nnodes = 4
        cen = 'CEN/4'
    elif self.element_type == 273:
        msg = cplstn6_msg
        nnodes = 6
        cen = 'CEN/6'
    elif self.element_type == 274:
         msg = cplstn8_msg
         nnodes = 8
         cen = 'CEN/8'
    else:  # pragma: no cover
        raise NotImplementedError(f'name={self.element_name} type={self.element_type}')
    return msg, nnodes, cen
