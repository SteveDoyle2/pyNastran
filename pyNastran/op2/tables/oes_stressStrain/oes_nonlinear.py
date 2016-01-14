from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import range
from math import isnan
from numpy import zeros, array_equal

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
try:
    import pandas as pd
except ImportError:
    pass


class RealNonlinearRodArray(OES_Object):
    """
    ::

      ELEMENT-ID =     102
                               N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                             STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
      2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
      3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        headers = ['axial_stress', 'equiv_stress', 'total_strain',
        'effective_plastic_creep_strain', 'effective_creep_strain',
        'linear_torsional_stress']
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

        #[axial_stress, equiv_stress, total_strain, effective_plastic_creep_strain,
        # effective_creep_strain, linear_torsional_stress]
        self.data = zeros((self.ntimes, self.nelements, 6), dtype='float32')

    def build_dataframe(self):
        headers = self.get_headers()
        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            df1 = pd.DataFrame(self.element).T
            df1.columns = ['ElementID']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join([df2])
        #print(self.data_frame)

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        (axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1) = t1
                        (axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2) = t2
                        if not allclose(t1, t2):
                        #if not array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                axial_stress1, equiv_stress1, total_strain1, effective_plastic_creep_strain1, effective_creep_strain1, linear_torsional_stress1,
                                axial_stress2, equiv_stress2, total_strain2, effective_plastic_creep_strain2, effective_creep_strain2, linear_torsional_stress2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial_stress, equiv_stress, total_strain,
                  effective_plastic_creep_strain, effective_creep_strain, linear_torsional_stress):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            axial_stress, equiv_stress, total_strain, effective_plastic_creep_strain,
            effective_creep_strain, linear_torsional_stress
        ]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if is_sort1:
            msg = [
                '                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                ' \n',
                '    ELEMENT-ID    AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL\n',
                '                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS\n'
            ]
        else:
            msg = [
                '                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                ' \n',
                '    TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL\n',
                '                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS\n'
            ]

        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f, msg)
        else:
            raise NotImplementedError('RealNonlinearRodArray')
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f, msg_temp):
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            eqs = self.data[itime, :, 1]
            total = self.data[itime, :, 2]
            epcs = self.data[itime, :, 3]
            ecs = self.data[itime, :, 4]
            lts = self.data[itime, :, 5]

            #print "dt=%s axials=%s eqs=%s ts=%s epcs=%s ecs=%s lts=%s" %(dt,axial,eqs,ts,epcs,ecs,lts)
            #msgE[eid] = '      ELEMENT-ID = %8i\n' % (eid)
            #if eid not in msgT:
                #msgT[eid] = []
            #msgT[eid].append('  %9.3E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E\n' % (dt, axial, eqs, ts, epcs, ecs, lts))

            for eid, axiali, eqsi, totali, epcsi, ecsi, ltsi in zip(eids, axial, eqs, total, epcs, ecs, lts):
                ([saxial, seqs, stotal, sepcs, secs, slts]) = write_floats_13e(
                    [axiali, eqsi, totali, epcsi, ecsi, ltsi])

            f.write('  %8i       %-13s       %-13s       %-13s       %-13s       %-13s       %s\n' % (
                eid, saxial, seqs, stotal, sepcs, secs, slts))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealNonlinearPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'fiberDistance', 'oxx', 'oyy', 'ozz', 'txy',
            'exx', 'eyy', 'ezz', 'exy', 'es', 'eps', 'ecs'
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
        #if self.element_type in [33, 74]:
            #nnodes_per_element = 1
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
        #else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

        self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        #[fiberDistance, oxx, oyy, ozz, txy, exx, eyy, ezz, exy, es, eps, ecs]
        self.data = zeros((self.ntimes, self.ntotal, 12), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'element_node shape=%s table.shape=%s' % (self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Eid, Nid)\n'
            for (eid, nid), (eid2, nid2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s)    (%s, %s)\n' % (eid, nid, eid2, nid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiberDistance1, oxx1, oyy1, ozz1, txy1, exx1, eyy1, ezz1, exy1, es1, eps1, ecs1) = t1
                    (fiberDistance2, oxx2, oyy2, ozz2, txy2, exx2, eyy2, ezz2, exy2, es2, eps2, ecs2) = t2

                    # vm stress can be NaN for some reason...
                    if not array_equal(t1[:-1], t2[:-1]):
                        msg += ('(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                                '            (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiberDistance1, oxx1, oyy1, ozz1, txy1, exx1, eyy1, ezz1, exy1, es1, eps1, ecs1,
                            fiberDistance2, oxx2, oyy2, ozz2, txy2, exx2, eyy2, ezz2, exy2, es2, eps2, ecs2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    #def add_new_eid(self, etype, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #self.add_new_eid_sort1(etype, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    #def add_new_eid_sort1(self, etype, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert isinstance(eid, int), eid
        #assert isinstance(node_id, int), node_id
        #self._times[self.itime] = dt
        ##assert self.itotal == 0, oxx
        #self.element_node[self.itotal, :] = [eid, node_id]
        #self.data[self.itime, self.itotal, :] = [fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        #self.itotal += 1
        #self.ielement += 1

    #def addNewNode(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert isinstance(node_id, int), node_id
        #self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    #def add(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert isinstance(node_id, int), node_id
        #self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    #def addNewNodeSort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    #def add_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None, eid
        #assert isinstance(node_id, int), node_id
        #self.element_node[self.itotal, :] = [eid, node_id]
        #self.data[self.itime, self.itotal, :] = [fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        #self.itotal += 1

    def get_stats(self):
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
            ntimes_word = 1
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        real_nonlinear_plate_array
        msg, nnodes, cen = _get_plate_msg(self)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        #cen_word = 'CEN/%i' % nnodes
        cen_word = cen
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            fiber_dist = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            angle = self.data[itime, :, 4]
            majorP = self.data[itime, :, 5]
            minorP = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]

            # loop over all the elements
            for (i, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi) in zip(
                 count(), eids, nids, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
                [fdi, oxxi, oyyi, txyi, major, minor, ovmi] = write_floats_13e(
                    [fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                ilayer = i % 2
                # tria3
                if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
                    if ilayer == 0:
                        f.write('0  %6i   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    else:
                        f.write('   %6s   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % ('', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))

                elif self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                    # bilinear
                    if nid == 0 and ilayer == 0:  # CEN
                        f.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (eid, cen_word, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 0:
                        f.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % ('', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 1:
                        f.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % ('', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                else:
                    raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

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

    def get_stats(self):
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

    def add_new_eid_sort1(self, eType, dt, data):
        if dt not in self.oxx:
            self.add_new_transient(dt)
        (eid, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy) = data
        self.fiberDistance[dt][eid] = [fd]
        if isnan(sz):
            sz = 0.
        if isnan(ez):
            ez = 0.
        self.eType[eid] = eType
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

    def add_sort1(self, dt, data):
        (eid, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy) = data
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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
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


class HyperelasticQuad(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'QUAD4FD'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.Type = {}
        self.IDs = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  Type, oxx, oyy, txy, angle, majorP, minorP\n')
        return msg

    def delete_transient(self, dt):
        del self.fiberDistance[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]

        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.txy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}

    def add_new_eid_sort1(self, dt, data):
        if dt not in self.oxx:
            self.add_new_transient(dt)
        (eid, Type, oxx, oyy, txy, angle, majorP, minorP) = data
        self.Type[eid] = Type
        self.oxx[dt] = {eid: [oxx]}
        self.oyy[dt] = {eid: [oyy]}
        self.txy[dt] = {eid: [txy]}
        self.angle[dt] = {eid: [angle]}
        self.majorP[dt] = {eid: [majorP]}
        self.minorP[dt] = {eid: [minorP]}

    def add_sort1(self, dt, eid, data):
        (ID, oxx, oyy, txy, angle, majorP, minorP) = data
        self.oxx[dt][eid].append(oxx)
        self.oyy[dt][eid].append(oyy)
        self.txy[dt][eid].append(txy)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):  # .. todo:: doesnt support CTRIA3NL (calls them CQUAD4s)
        msg = ['           S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )\n',
               '  ELEMENT     GRID/    POINT       ---------CAUCHY STRESSES--------             PRINCIPAL STRESSES (ZERO SHEAR)\n',
               '     ID       GAUSS      ID      NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR\n', ]
               #0       1     GAUS         1   7.318995E+00   6.367099E-01  -6.551054E+00   -31.4888    1.133173E+01   -3.376026E+00
               #                           2   1.097933E+01   4.149028E+00   6.278160E+00    30.7275    1.471111E+01    4.172537E-01

        for dt, Oxxs in sorted(iteritems(self.oxx)):
            #header[-1] = '     LOAD STEP = %12.5E' %(dt)
            msg += header
            for eid, oxxs in sorted(iteritems(Oxxs)):
                gauss = self.Type[eid]
                oxx = self.oxx[dt][eid]
                oyy = self.oyy[dt][eid]
                txy = self.txy[dt][eid]
                angle = self.angle[dt][eid]
                majorP = self.majorP[dt][eid]
                minorP = self.minorP[dt][eid]

                for i in range(4):  # 1,2,3,4
                    if i == 0:
                        msg.append('0%8i %8s  %8i  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6\n' % (eid, gauss, i + 1, oxx[i], oyy[i], txy[i], angle[i], majorP[i], minorP[i]))
                    else:
                        msg.append(' %8s %8s  %8i  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6\n' % ('', '', i + 1, oxx[i], oyy[i], txy[i], angle[i], majorP[i], minorP[i]))
        f.write(''.join(msg))
        return page_num


#class NonlinearRod(StressObject):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StressObject.__init__(self, data_code, isubcase)
        ##self.eType = 'CROD'
        #self.eTypeMap = {89: 'CRODNL', 92: 'CONRODNL'}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.eType = {}
        #self.axialStress = {}
        #self.equivStress = {}
        #self.totalStrain = {}
        #self.effectivePlasticCreepStrain = {}
        #self.effectiveCreepStrain = {}
        #self.linearTorsionalStress = {}

        #self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                ##self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            ##self.add = self.add_sort2
            ##self.add_new_eid = self.add_new_eid_sort2

    #def get_stats(self):
        #nelements = len(self.eType)
        #msg = self.get_data_code()
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.axialStress)
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, axialStress, equivStress, totalStrain, '
                   #'effectivePlasticCreepStrain, effectiveCreepStrain, '
                   #'linearTorsionalStress\n')
        #return msg

    #def delete_transient(self, dt):
        #del self.axialStress[dt]
        #del self.equivStress[dt]
        #del self.totalStrain[dt]
        #del self.effectivePlasticCreepStrain[dt]

        #del self.effectiveCreepStrain[dt]
        #del self.linearTorsionalStress[dt]

    #def get_transients(self):
        #k = self.axialStress.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #self.axialStress[dt] = {}
        #self.equivStress[dt] = {}
        #self.totalStrain[dt] = {}
        #self.effectivePlasticCreepStrain[dt] = {}
        #self.effectiveCreepStrain[dt] = {}
        #self.linearTorsionalStress[dt] = {}

    #def add_sort1(self, eType, dt, data):
        #if dt not in self.axialStress:
            #self.add_new_transient(dt)
        #eid = data[0]
        #self.eType[eid] = eType
        #self.axialStress[dt][eid] = data[1]
        #self.equivStress[dt][eid] = data[2]
        #self.totalStrain[dt][eid] = data[3]
        #self.effectivePlasticCreepStrain[dt][eid] = data[4]
        #self.effectiveCreepStrain[dt][eid] = data[5]
        #self.linearTorsionalStress[dt][eid] = data[6]
        ##print data

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):  # .. todo:: doesnt support CONROD/CTUBE (calls them CRODs)
        #"""
        #::

          #ELEMENT-ID =     102
                                   #N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
            #TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                                 #STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
          #2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
          #3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
        #"""
        #msg = []
        #msg_start = [
            #'                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
            #' \n',
            #'    TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL\n',
            #'                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS\n'
        #]
        #msg_element = {}
        #msg_time = {}
        #for dt, axials in sorted(iteritems(self.axialStress)):
            #for eid, axial in sorted(iteritems(axials)):
                #eqs = self.equivStress[dt][eid]
                #ts = self.totalStrain[dt][eid]
                #epcs = self.effectivePlasticCreepStrain[dt][eid]
                #ecs = self.effectiveCreepStrain[dt][eid]
                #lts = self.linearTorsionalStress[dt][eid]
                ##print "dt=%s axials=%s eqs=%s ts=%s epcs=%s ecs=%s lts=%s" %(dt,axial,eqs,ts,epcs,ecs,lts)
                #msg_element[eid] = '      ELEMENT-ID = %8i\n' % (eid)
                #if eid not in msg_time:
                    #msg_time[eid] = []
                #msg_time[eid].append('  %9.3E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E\n' % (dt, axial, eqs, ts, epcs, ecs, lts))

        #for eid, e in sorted(iteritems(msg_element)):
            #msg += header + [e] + msg_start + msg_time[eid]
            #msg.append(page_stamp % page_num)
        #f.write(''.join(msg))
        #return page_num
