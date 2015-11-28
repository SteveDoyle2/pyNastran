#pylint disable=C0103,C0301
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import zeros, empty, array_equal
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import get_key0, writeFloats13E, _eigenvalue_header
from pyNastran.op2.resultObjects.element_table_object import RealElementTableArray



class Real1DHeatFluxArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

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
        self.element_type = empty(self.nelements, dtype='|U8')

        #[xgrad, ygrad, zgrad, xflux, yflux, zflux]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for (eid, nid), (eid2, nid2) in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
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

                    if not array_equal(t1, t2):
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
        self.element_type[self.ielement] = etype
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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        onedeeheatflux
        msg_temp = [
            '             F O R C E S   I N   A X I S - S Y M M E T R I C   C O N I C A L   S H E L L   E L E M E N T S   (CCONEAX)\n'
            ' \n'
            '  ELEMENT     HARMONIC    POINT           BEND-MOMENT       BEND-MOMENT      TWIST-MOMENT           SHEAR            SHEAR\n'
            '   ID.         NUMBER     ANGLE               V                 U                                     V                U\n'
            #'      101        0                       5.864739E-09      1.759422E-09      0.0               0.0               0.0'
            #'      101                0.0000          5.864739E-09      1.759422E-09      0.0               0.0               0.0'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            xgrad = self.data[itime, :, 0]
            ygrad = self.data[itime, :, 1]
            zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 3]
            yflux = self.data[itime, :, 4]
            zflux = self.data[itime, :, 5]

            for (eid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi) in zip(
                eids, xgrad, ygrad, zgrad, xflux, yflux, zflux):
                (vals2, is_all_zeros) = writeFloats13E(
                    [xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi])
                [xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi] = vals2

                # TODO: hopa is probably the wrong type
                f.write(' %8i  %-13s %-13s %-13s %-13s %-13s %s\n' % (
                    eid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi))
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

    def add(self, nnodes, dt, data):
        [eid, parent, grad_fluxes] = data
        self.parent[eid] = parent
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[eid][nid] = [xflux, yflux, zflux]

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, grad_fluxes] = data
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

    def add_sort2(self, nnodes, eid, data):
        (dt, parent, grad_fluxes) = data
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


class HeatFlux_VU(ScalarObject):  # 189-VUQUAD 190-VUTRIA,191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        #self.eType = {}
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

    def add(self, nnodes, dt, data):
        [eid, parent, coord, icord, theta, grad_fluxes] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[eid][nid] = [xflux, yflux, zflux]

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, coord, icord, theta, grad_fluxes] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]

    def add_sort2(self, nnodes, eid, data):
        [dt, parent, coord, icord, theta, grad_fluxes] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for grad_flux in grad_fluxes:
            [nid, xgrad, ygrad, zgrad, xflux, yflux, zflux] = grad_flux
            self.grad[dt][eid][nid] = [xgrad, ygrad, zgrad]
            self.flux[dt][eid][nid] = [xflux, yflux, zflux]


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


class HeatFlux_1D(ScalarObject):  # 1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.eType = {}
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
        nelements = len(self.eType)
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.grad)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, grad, flux\n')
        return msg

    def add_new_transient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, dt, data):
        [eid, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux] = data
        self.eType[eid] = eType
        self.grad[eid] = [xgrad, ygrad, zgrad]
        self.flux[eid] = [xflux, yflux, zflux]

    def add_sort1(self, dt, data):
        [eid, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux] = data
        if dt not in self.grad:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xgrad, ygrad, zgrad]
        self.flux[dt][eid] = [xflux, yflux, zflux]

    def add_sort2(self, eid, data):
        [dt, eType, xgrad, ygrad, zgrad, xflux, yflux, zflux] = data
        if dt not in self.fApplied:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xgrad, ygrad, zgrad]
        self.flux[dt][eid] = [xflux, yflux, zflux]


class HeatFlux_2D_3DArray(RealElementTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealElementTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #words = ['                                             D I S P L A C E M E N T   V E C T O R\n', ]
                 ##' \n',
                 ##'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        ##words += self.get_table_marker()
        #write_words = True
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words,
                                                   #is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                                         #is_mag_phase=is_mag_phase, is_sort1=is_sort1)

    def _get_headers(self):
        return ['grad1', 'grad2', 'grad3', 'flux1', 'flux2', 'flux3']


class HeatFlux_CONV(ScalarObject):  # 110-CONV
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.cntlNode = {}
        self.freeConv = {}
        self.freeConvK = {}

        # TODO if dt=None, handle SORT1 case
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.cntlNode)
            times0 = get_key0(self.cntlNode)
            nelements = len(self.cntlNode[times0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.cntlNode)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  cntlNode, freeConv, freeConvK\n')
        return msg

    def add_new_transient(self, dt):
        self.cntlNode[dt] = {}
        self.freeConv[dt] = {}
        self.freeConvK[dt] = {}

    def add(self, dt, data):
        [eid, cntl_node, freeConv, freeConvK] = data
        #self.eType[eid]     = eType
        self.cntlNode[eid] = cntl_node
        self.freeConv[eid] = freeConv
        self.freeConvK[eid] = freeConvK

    def add_sort1(self, dt, data):
        [eid, cntl_node, freeConv, freeConvK] = data
        if dt not in self.freeConv:
            self.add_new_transient(dt)
        #self.eType[eid]     = eType
        self.cntlNode[dt][eid] = cntl_node
        self.freeConv[dt][eid] = freeConv
        self.freeConvK[dt][eid] = freeConvK

    def add_sort2(self, eid, data):
        [dt, eType, fApplied, freeConv, forceConv, fRad, fTotal] = data
        if dt not in self.freeConv:
            self.add_new_transient(dt)
        #self.eType[eid]     = eType
        self.fApplied[dt][eid] = fApplied
        self.freeConv[dt][eid] = freeConv
        self.forceConv[dt][eid] = forceConv
        self.fRad[dt][eid] = fRad
        self.fTotal[dt][eid] = fTotal


class RealChbdyHeatFluxArray(ScalarObject):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

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
        self.element_type = empty(self.nelements, dtype='|U8')

        #[fapplied, free_conv, force_conv, frad, ftotal]
        self.data = zeros((self.ntimes, self.ntotal, 5), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for (eid, nid), (eid2, nid2) in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
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

                    if not array_equal(t1, t2):
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
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        #msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        onedeeheatflux
        msg_temp = [
            '             F O R C E S   I N   A X I S - S Y M M E T R I C   C O N I C A L   S H E L L   E L E M E N T S   (CCONEAX)\n'
            ' \n'
            '  ELEMENT     HARMONIC    POINT           BEND-MOMENT       BEND-MOMENT      TWIST-MOMENT           SHEAR            SHEAR\n'
            '   ID.         NUMBER     ANGLE               V                 U                                     V                U\n'
            #'      101        0                       5.864739E-09      1.759422E-09      0.0               0.0               0.0'
            #'      101                0.0000          5.864739E-09      1.759422E-09      0.0               0.0               0.0'
        ]

        #(elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            xgrad = self.data[itime, :, 0]
            ygrad = self.data[itime, :, 1]
            zgrad = self.data[itime, :, 2]
            xflux = self.data[itime, :, 3]
            yflux = self.data[itime, :, 4]
            zflux = self.data[itime, :, 5]

            for (eid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi) in zip(
                eids, xgrad, ygrad, zgrad, xflux, yflux, zflux):
                (vals2, is_all_zeros) = writeFloats13E(
                    [xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi])
                [xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi] = vals2

                # TODO: hopa is probably the wrong type
                f.write(' %8i  %-13s %-13s %-13s %-13s %-13s %s\n' % (
                    eid, xgradi, ygradi, zgradi, xfluxi, yfluxi, zfluxi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class HeatFlux_CHBDYx(ScalarObject):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.fApplied = {}
        self.freeConv = {}
        self.forceConv = {}
        self.fRad = {}
        self.fTotal = {}

        # TODO if dt=None, handle SORT1 case
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self.eType)
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.fApplied)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fApplied, freeConv, forceConv, fRad, fTotal\n')
        return msg

    def add_new_transient(self, dt):
        self.fApplied[dt] = {}
        self.freeConv[dt] = {}
        self.forceConv[dt] = {}
        self.fRad[dt] = {}
        self.fTotal[dt] = {}

    def add(self, dt, data):
        [eid, etype, fapplied, free_conv, force_conv, frad, ftotal] = data

        self.eType[eid] = etype
        self.fApplied[eid] = fapplied
        self.freeConv[eid] = free_conv
        self.forceConv[eid] = force_conv
        self.fRad[eid] = frad
        self.fTotal[eid] = ftotal

    def add_sort1(self, dt, data):
        [eid, etype, fapplied, free_conv, force_conv, frad, ftotal] = data
        if dt not in self.fApplied:
            self.add_new_transient(dt)
        self.eType[eid] = etype
        self.fApplied[dt][eid] = fapplied
        self.freeConv[dt][eid] = free_conv
        self.forceConv[dt][eid] = force_conv
        self.fRad[dt][eid] = frad
        self.fTotal[dt][eid] = ftotal

    def add_sort2(self, eid, data):
        [dt, etype, fapplied, free_conv, force_conv, frad, ftotal] = data
        if dt not in self.fApplied:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.fApplied[dt][eid] = fapplied
        self.freeConv[dt][eid] = free_conv
        self.forceConv[dt][eid] = force_conv
        self.fRad[dt][eid] = frad
        self.fTotal[dt][eid] = ftotal
