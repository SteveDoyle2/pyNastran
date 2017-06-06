from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import isnan
from itertools import cycle
from six import iteritems, integer_types
from six.moves import range
import numpy as np
from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header, write_float_13e
try:
    import pandas as pd
except ImportError:
    pass


class RealNonlinearPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        #asdfasdfasdf
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def is_stress(self):
        return True

    def get_headers(self):
        headers = [
            #[fiber_dist, oxx, oyy, ozz, txy, es, eps, ecs, exx, eyy, ezz, etxy]
            'fiber_distance', 'oxx', 'oyy', 'ozz', 'txy',
            'eff_plastic_strain', 'eff_plastic_strain', 'eff_creep_strain',
            'exx', 'eyy', 'ezz', 'exy',
        ]
        return headers

    #def is_bilinear(self):
        #if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            #return False
        #elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            #return True
        #else:
            #raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))

    def build(self):
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 88:  # Tria3-nonlinear
            nnodes_per_element = 1
        elif self.element_type == 90:  # Quad4-nonlinear
            nnodes_per_element = 1
        #elif self.element_type == 144:
            #nnodes_per_element = 5
        #elif self.element_type == 64:  # CQUAD8
            #nnodes_per_element = 5
        #elif self.element_type == 82:  # CQUADR
            #nnodes_per_element = 5
        #elif self.element_type == 70:  # CTRIAR
            #nnodes_per_element = 4
        #elif self.element_type == 75:  # CTRIA6
            #nnodes_per_element = 4
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        nnodes_per_element = 1
        self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element

        if self.nelements % self.ntimes != 0:
            msg = 'nelements=%s ntimes=%s nelements/ntimes=%s'  % (
                self.nelements, self.ntimes, self.nelements / float(self.ntimes))
            #return
            raise RuntimeError(msg)

        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[fiber_dist, oxx, oyy, ozz, txy, es, eps, ecs, exx, eyy, ezz, etxy]
        self.data = zeros((self.ntimes, self.ntotal, 12), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()[1:]

        nelements = self.element.shape[0]
        if self.is_fiber_distance():
            fiber_distance = ['Top', 'Bottom'] * nelements
        else:
            fiber_distance = ['Mean', 'Curvature'] * nelements
        fd = np.array(fiber_distance, dtype='unicode')
        element = np.vstack([self.element, self.element]).T.ravel()
        element_fd = [element, fd]

        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data[:, :, 1:], items=column_values, major_axis=element_fd, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Location', 'Item']
        else:
            # option B - nice!
            df1 = pd.DataFrame(element_fd).T
            df1.columns = ['ElementID', 'Location']
            df2 = pd.DataFrame(self.data[0, :, 1:])
            df2.columns = headers
            self.data_frame = df1.join(df2)
        self.data_frame = self.data_frame.reset_index().set_index(['ElementID', 'Location'])
        #print(self.data_frame)

    #def add_new_eid(self, dt, eid, etype, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy):
        #self.add_sort1(dt, eid, etype, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy)

    def add_new_eid_sort1(self, dt, eid, etype, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy):
        self.element[self.ielement] = eid
        self.ielement += 1
        self.add_sort1(dt, eid, etype, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy)

    def add_sort1(self, dt, eid, etype, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy):
        """unvectorized method for adding SORT1 transient data"""
        if isnan(fd):
            fd = 0.
        if isnan(sz):
            sz = 0.
        if isnan(ez):
            ez = 0.
        self._times[self.itime] = dt
        #if self.ielement == 10:
            #print(self.element_node[:10, :])
            #raise RuntimeError()
        #[fiber_dist, oxx, oyy, ozz, txy, es, eps, ecs, exx, eyy, ezz, etxy]
        assert eid == self.element[self.ielement - 1], 'eid=%s self.element[i-1] = ' % (eid, self.element[self.ielement - 1])
        self.data[self.itime, self.itotal, :] = [fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy]
        self.itotal += 1

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]

                    # TODO: this name order is wrong
                    #[fiber_dist, oxx, oyy, ozz, txy, es, eps, ecs, exx, eyy, ezz, etxy]
                    (fiber_distance1, oxx1, oyy1, ozz1, txy1, exx1, eyy1, ezz1, exy1, es1, eps1, ecs1) = t1
                    (fiber_distance2, oxx2, oyy2, ozz2, txy2, exx2, eyy2, ezz2, exy2, es2, eps2, ecs2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1, t2):
                        msg += (
                            '%s    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '%s    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                fiber_distance1, oxx1, oyy1, ozz1, txy1, exx1, eyy1, ezz1, exy1, es1, eps1, ecs1,
                                ' ' * (len(str(eid)),
                                fiber_distance2, oxx2, oyy2, ozz2, txy2, exx2, eyy2, ezz2, exy2, es2, eps2, ecs2)))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.nnodes
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // self.nnodes // 2

        msg = []
        if self.nonlinear_factor is not None:  # transient
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
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg, nnodes, cen = _get_plate_msg(self)
        if self.element_type == 88:
            msg = [
                '                   N O N L I N E A R   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S      ( T R I A 3 )\n'
                ' \n'
                '    ELEMENT      FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP\n'
                '       ID      DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN\n'
            ]
        elif self.element_type == 90:
            msg = [
                '               N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D 4 )\n'
                ' \n'
                '    ELEMENT      FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP\n'
                '       ID      DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN\n'
                #'0         1  -2.500000E-02  -4.829193E+00  -1.640651E-05                 -1.907010E-04   4.829185E+00   0.0            0.0\n'
                #'                            -4.829188E-05   1.448741E-05                 -4.958226E-09\n'
                #'              2.500000E-02   4.770547E+00   1.493975E-04                  1.907012E-04   4.770473E+00   0.0            0.0\n'
                #'                             4.770502E-05  -1.431015E-05                  4.958231E-09\n'
            ]
        else:
            raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

        #msg = [
        #'               N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D 4 )\n'
        #' \n'
        #'    ELEMENT      FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP\n'
        #'       ID      DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN\n'
        ##'0         1  -2.500000E-02  -4.829193E+00  -1.640651E-05                 -1.907010E-04   4.829185E+00   0.0            0.0\n'
        ##'                            -4.829188E-05   1.448741E-05                 -4.958226E-09\n'
        ##'              2.500000E-02   4.770547E+00   1.493975E-04                  1.907012E-04   4.770473E+00   0.0            0.0\n'
        ##'                             4.770502E-05  -1.431015E-05                  4.958231E-09\n'
        #]


        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        #cen_word = 'CEN/%i' % nnodes
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[fiber_dist, oxx, oyy, ozz, txy, es, eps, ecs, exx, eyy, ezz, etxy]
            fiber_dist = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            ozz = self.data[itime, :, 3]
            txy = self.data[itime, :, 4]
            es = self.data[itime, :, 5]
            eps = self.data[itime, :, 6]
            ecs = self.data[itime, :, 7]
            exx = self.data[itime, :, 8]
            eyy = self.data[itime, :, 9]
            ezz = self.data[itime, :, 10]
            exy = self.data[itime, :, 11]

            for (i, eid, nid, fdi, oxxi, oyyi, ozzi, txyi, exxi, eyyi, ezzi, exyi, esi, epsi, ecsi) in zip(
                 cycle([0, 1]), eids, nids, fiber_dist, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs):
                #[fdi, oxxi, oyyi, txyi, major, minor, ovmi] = write_floats_13e(
                    #[fdi, oxxi, oyyi, txyi, major, minor, ovmi])

                #'    ELEMENT      FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP\n'
                #'       ID      DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN\n'
                #'0         1  -2.500000E-02  -4.829193E+00  -1.640651E-05                 -1.907010E-04   4.829185E+00   0.0            0.0\n'
                #'                            -4.829188E-05   1.448741E-05                 -4.958226E-09\n'
                #'              2.500000E-02   4.770547E+00   1.493975E-04                  1.907012E-04   4.770473E+00   0.0            0.0\n'
                #'                             4.770502E-05  -1.431015E-05                  4.958231E-09\n'
                if i == 0:
                    f.write(
                        '0  %8i  %-13s  %-13s  %-13s                 %-13s  %-13s  %-13s  %s\n'
                        '                            %-13s  %-13s                 %s\n' % (
                            # A
                            #eid, write_float_13e(fdi),
                            #write_float_13e(oxxi), write_float_13e(oyyi),
                            ##write_float_13e(ozzi),
                            #write_float_13e(txyi),
                            #write_float_13e(esi), write_float_13e(epsi),
                            #write_float_13e(ecsi),

                            #write_float_13e(exxi), write_float_13e(eyyi),
                            ##write_float_13e(ezzi),
                            #write_float_13e(exyi),

                            #    ELEMENT  FIBER  XYZ STRESS     EQUIVALENT  EFF.STRAIN  EFF.CREEP\n'
                            eid, write_float_13e(fdi),
                            write_float_13e(oxxi), write_float_13e(oyyi),
                            #write_float_13e(ozzi),
                            write_float_13e(txyi),

                            write_float_13e(esi), write_float_13e(epsi),
                            write_float_13e(ecsi),

                            write_float_13e(exxi), write_float_13e(eyyi),
                            #write_float_13e(ezzi),
                            write_float_13e(exyi),
                        ))
                else:
                    f.write(
                        '             %-13s  %-13s  %-13s                 %-13s  %-13s  %-13s  %s\n'
                        '                            %-13s  %-13s                 %s\n' % (
                            write_float_13e(fdi),
                            write_float_13e(oxxi), write_float_13e(oyyi),
                            #write_float_13e(ozzi),
                            write_float_13e(txyi),

                            write_float_13e(esi), write_float_13e(epsi),
                            write_float_13e(ecsi),

                            write_float_13e(exxi), write_float_13e(eyyi),
                            #write_float_13e(ezzi),
                            write_float_13e(exyi),
                    ))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

class NonlinearQuad(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        #self.eType = 'QUAD4FD' # or CTRIA3

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.eType = {}
        self.fiberDistance = {}
        self.oxx = {}
        self.oyy = {}
        self.ozz = {}
        self.txy = {}

        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}

        self.es = {}
        self.eps = {}
        self.ecs = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self, short=False):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberDistance, oxx, oyy, ozz, txy, '
                   'exx, eyy, ezz, exy, es, eps, ecs\n')
        return msg

    def delete_transient(self, dt):
        del self.fiberDistance[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.ozz[dt]
        del self.txy[dt]

        del self.exx[dt]
        del self.eyy[dt]
        del self.ezz[dt]
        del self.exy[dt]

        del self.es[dt]
        del self.eps[dt]
        del self.ecs[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        self.fiberDistance[dt] = {}
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.ozz[dt] = {}
        self.txy[dt] = {}

        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.ezz[dt] = {}
        self.exy[dt] = {}

        self.es[dt] = {}
        self.eps[dt] = {}
        self.ecs[dt] = {}

    def add_new_eid_sort1(self, dt, eid, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy):
        if dt not in self.oxx:
            self.add_new_transient(dt)
        self.fiberDistance[dt][eid] = [fd]
        if isnan(sz):
            sz = 0.
        if isnan(ez):
            ez = 0.
        self.oxx[dt][eid] = [sx]
        self.oyy[dt][eid] = [sy]
        self.ozz[dt][eid] = [sz]
        self.txy[dt][eid] = [txy]

        self.exx[dt][eid] = [ex]
        self.eyy[dt][eid] = [ey]
        self.ezz[dt][eid] = [ez]
        self.exy[dt][eid] = [exy]

        self.es[dt][eid] = [es]
        self.eps[dt][eid] = [eps]
        self.ecs[dt][eid] = [ecs]

    def add_sort1(self, dt, eid, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy):
        """unvectorized method for adding SORT1 transient data"""
        self.fiberDistance[dt][eid].append(fd)
        if isnan(sz):
            sz = 0.
        if isnan(ez):
            ez = 0.

        self.oxx[dt][eid].append(sx)
        self.oyy[dt][eid].append(sy)
        self.ozz[dt][eid].append(sz)
        self.txy[dt][eid].append(txy)

        self.exx[dt][eid].append(ex)
        self.eyy[dt][eid].append(ey)
        self.ezz[dt][eid].append(ez)
        self.exy[dt][eid].append(exy)

        self.es[dt][eid].append(es)
        self.eps[dt][eid].append(eps)
        self.ecs[dt][eid].append(ecs)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_start = [
            '      ELEMENT-ID =     129\n'
            '               N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D 4 )\n'
            ' \n',
            '    TIME         FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP\n'
            '               DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN\n'
        ]
        #0 5.000E-05  -5.000000E-01  -4.484895E+01  -1.561594E+02                 -2.008336E-02   1.392609E+02   0.0            0.0
        msg_element = {}
        msg_time = {}
        for (dt, Oxxs) in sorted(iteritems(self.oxx)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)

            for (eid, oxxs) in sorted(iteritems(Oxxs)):
                msg_element[eid] = header + ['      ELEMENT-ID = %8i\n' % (eid)]
                if eid not in msg_time:
                    msg_time[eid] = []
                for i, oxx in enumerate(oxxs):
                    fd = self.fiberDistance[dt][eid][i]
                    oxx = self.oxx[dt][eid][i]
                    oyy = self.oyy[dt][eid][i]
                    ozz = self.ozz[dt][eid][i]
                    txy = self.txy[dt][eid][i]

                    exx = self.exx[dt][eid][i]
                    eyy = self.eyy[dt][eid][i]
                    ezz = self.ezz[dt][eid][i]
                    exy = self.exy[dt][eid][i]

                    es = self.es[dt][eid][i]
                    eps = self.eps[dt][eid][i]
                    ecs = self.ecs[dt][eid][i]
                    [oxx, oyy, ozz, txy, exx, eyy, es, eps, ecs, exx, eyy, ezz, exy] = write_floats_13e([oxx, oyy, ozz, txy, exx, eyy, es, eps, ecs, exx, eyy, ezz, exy])
                    if i == 0:
                        msg_time[eid].append('0 %9.3E %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (dt, fd, oxx, oyy, ozz, txy, es, eps, ecs))
                    else:
                        msg_time[eid].append('     %9s %-13s  %-13s  %-13s  %-13s  %-13s\n' % ('', '', exx, eyy, ezz, exy))

        msg = []
        for eid, e in sorted(iteritems(msg_element)):
            msg += header + e + msg_start + msg_time[eid]
            msg.append(page_stamp % page_num)
            page_num += 1
        f.write(''.join(msg))
        return page_num - 1
