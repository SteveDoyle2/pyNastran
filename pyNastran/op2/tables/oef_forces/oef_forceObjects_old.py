from numpy import array
from six import iteritems
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e, get_key0, writeFloats12E
from pyNastran.f06.f06_formatting import _eigenvalue_header
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import _write_f06_springs


class RealRodForce(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)

        self.axial_force = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def add_new_transient(self, dt):
        self.dt = dt
        self.axial_force[dt] = {}
        self.torque[dt] = {}

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            #if ntimes == 0:
            #    time0 = None
            #    nelements = None
            #else:
            time0 = get_key0(self.torque)
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append(' axial_force, torque\n')
        return msg

    def add(self, dt, eid, axial_force, torque):
        self.axial_force[eid] = axial_force
        self.torque[eid] = torque

    def add_sort1(self, dt, eid, axial_force, torque):
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, data):
        [dt, axial_force, torque] = data
        if dt not in self.axial_force:
            self.add_new_transient(dt)

        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torque) = line
                self.axial_force[eid] = axial
                self.torque[eid] = torque
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        #if dt not in self.axial_force:
            #self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial_force[dt][eid] = axial
            self.torque[dt][eid] = torsion

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        msg = []
        itime = 0
        ntimes = len(self.axial_force)
        for dt, axials in sorted(iteritems(self.axial_force)):
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            #dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            #header[2] = dtLine
            msg += header
            out = []
            for eid in sorted(axials):
                axial = self.axial_force[dt][eid]
                torsion = self.torque[dt][eid]
                out.append([eid, axial, torsion])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in range(0, nWrite, 2):


                #            14       -3.826978E+05   2.251969E+03                       15       -4.054566E+05   8.402607E+02


                out_line = '      %8i       %13.6E  %13.6E                 %8i       %13.6E  %13.6E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(out_line)

            if nOut % 2 == 1:
                out_line = '      %8i       %13.6E  %13.6E\n' % (
                    tuple(out[-1]))
                msg.append(out_line)
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            page_num += 1
            itime += 1
        return page_num - 1

    def _write_f06(self, msg, page_stamp, page_num, f):
        out = []
        for eid in sorted(self.axial_force):
            axial = self.axial_force[eid]
            torsion = self.torque[eid]
            vals2 = write_floats_13e([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, torsion])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in range(0, nWrite, 2):
            out_line = '      %8i   %-13s  %-13s  %8i   %-13s  %-13s\n' % tuple(out[i] + out[i + 1])
            msg.append(out_line)

        if nOut % 2 == 1:
            out_line = '      %8i   %-13s  %-13s\n' % tuple(out[-1])
            msg.append(out_line)
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num


class RealCTubeForce(RealRodForce):  # 3-TUBE
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealRodForce.__init__(self, data_code, is_sort1, isubcase, dt)
        self.elementType = 'CTUBE'

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = header + ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C T U B E )\n',
                          '       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n',
                          '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(words, page_stamp, page_num, f, is_mag_phase, is_sort1)
        return self._write_f06(words, page_stamp, page_num, f)


class RealConrodForce(RealRodForce):  # 10-CONROD
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealRodForce.__init__(self, data_code, is_sort1, isubcase, dt)
        self.elementType = 'CONROD'

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = header + ['                                           F O R C E S   I N   R O D   E L E M E N T S     ( C O N R O D )\n',
                          '       ELEMENT           AXIAL                                     ELEMENT           AXIAL\n',
                          '         ID.             FORCE          TORQUE                       ID.             FORCE          TORQUE\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(words, page_stamp, page_num, f, is_mag_phase, is_sort1)
        return self._write_f06(words, page_stamp, page_num, f)


class RealCRodForce(RealRodForce):  # 1-CROD
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        RealRodForce.__init__(self, data_code, is_sort1, isubcase, dt)
        self.elementType = 'CROD'

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = header + ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                          '       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n',
                          '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(words, page_stamp, page_num, f, is_mag_phase, is_sort1)
        return self._write_f06(words, page_stamp, page_num, f)


class RealPlateBilinearForce(ScalarObject):  # 64-CQUAD8, 75-CTRIA6, 82-CQUADR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.term = {}
        self.ngrids = {}
        self.mx = {}
        self.my = {}
        self.mxy = {}
        self.bmx = {}
        self.bmy = {}
        self.bmxy = {}
        self.tx = {}
        self.ty = {}
        self.eid_old = None

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add_new_element = self.addNewElementSort1
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add_new_element = self.addNewElementSort2
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.mx)
            if ntimes == 0:
                nelements = 0
            else:
                time0 = get_key0(self.mx)
                nelements = len(self.mx[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.mx)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  term, ngrids, mx, my, mxy, bmx, bmy, bmxy, tx, ty\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.ngrids = {}
        self.mx[dt] = {}
        self.my[dt] = {}
        self.mxy[dt] = {}
        self.bmx[dt] = {}
        self.bmy[dt] = {}
        self.bmxy[dt] = {}
        self.tx[dt] = {}
        self.ty[dt] = {}

    def add_f06_data(self, dt, data, element_name, ngrids):
        if dt is None:
            term = None
            eid_old = self.eid_old
            for (eid, grid, fx, fy, fxy, mx, my, mxy, qx, qy) in data:
                if isinstance(eid, int):
                    self.term[eid] = term
                    self.ngrids[eid] = [grid]
                    self.mx[eid] = [fx]
                    self.my[eid] = [fy]
                    self.mxy[eid] = [fxy]
                    self.bmx[eid] = [mx]
                    self.bmy[eid] = [my]
                    self.bmxy[eid] = [mxy]
                    self.tx[eid] = [qx]
                    self.ty[eid] = [qy]
                    eid_old = eid
                else:
                    eid = eid_old
                    self.ngrids[eid].append(grid)
                    self.mx[eid].append(fx)
                    self.my[eid].append(fy)
                    self.mxy[eid].append(fxy)
                    self.bmx[eid].append(mx)
                    self.bmy[eid].append(my)
                    self.bmxy[eid].append(mxy)
                    self.tx[eid].append(qx)
                    self.ty[eid].append(qy)
            self.eid_old = eid_old
        else:
            if dt not in self.mx:
                pass
                #raise NotImplementedError()

    def add_new_element(self, eid, dt, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.term[eid] = term
        self.ngrids[eid] = [nid]
        self.mx[eid] = [mx]
        self.my[eid] = [my]
        self.mxy[eid] = [mxy]
        self.bmx[eid] = [bmx]
        self.bmy[eid] = [bmy]
        self.bmxy[eid] = [bmxy]
        self.tx[eid] = [tx]
        self.ty[eid] = [ty]

    def add(self, eid, dt, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.ngrids[eid].append(nid)
        self.mx[eid].append(mx)
        self.my[eid].append(my)
        self.mxy[eid].append(mxy)
        self.bmx[eid].append(bmx)
        self.bmy[eid].append(bmy)
        self.bmxy[eid].append(bmxy)
        self.tx[eid].append(tx)
        self.ty[eid].append(ty)

    def addNewElementSort1(self, eid, dt, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)

        self.term[eid] = term
        self.ngrids[eid] = [nid]
        self.mx[dt][eid] = [mx]
        self.my[dt][eid] = [my]
        self.mxy[dt][eid] = [mxy]
        self.bmx[dt][eid] = [bmx]
        self.bmy[dt][eid] = [bmy]
        self.bmxy[dt][eid] = [bmxy]
        self.tx[dt][eid] = [tx]
        self.ty[dt][eid] = [ty]

    def add_sort1(self, eid, dt, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)

        self.ngrids[eid].append(nid)
        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def addNewElementSort2(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)

        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mx[dt][eid] = [mx]
        self.my[dt][eid] = [my]
        self.mxy[dt][eid] = [mxy]
        self.bmx[dt][eid] = [bmx]
        self.bmy[dt][eid] = [bmy]
        self.bmxy[dt][eid] = [bmxy]
        self.tx[dt][eid] = [tx]
        self.ty[dt][eid] = [ty]

    def add_sort2(self, dt, eid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)

        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=False, is_sort1=is_sort1)

        words = header + [
            '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n'
            ' \n'
            '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
            '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
        ]
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        f.write(''.join(words))
        for eid in sorted(self.mx):
            mxii = self.mx[eid]
            term = self.term[eid]
            #self.term[eid] = term
            for i in range(len(mxii)):
                node_id = self.ngrids[eid][i]
                mxi = self.mx[eid][i]
                myi = self.my[eid][i]
                mxyi = self.mxy[eid][i]
                bmxi = self.bmx[eid][i]
                bmyi = self.bmy[eid][i]
                bmxyi = self.bmxy[eid][i]
                txi = self.tx[eid][i]
                tyi = self.ty[eid][i]
                vals2 = write_floats_13e([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                [mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = vals2
                if i == 0:
                    f.write('0  %8i    CEN/4 %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                else:
                    f.write('            %8i %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (node_id, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
#            1     2.504029E+06  9.728743E+06   5.088001E+05  1.976808E+06   1.995229E+06  7.751935E+06  -3.684978E-07  -1.180941E-07

        f.write(page_stamp % page_num)
        return page_num


class RealPlateForce(ScalarObject):  # 33-CQUAD4, 74-CTRIA3
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.mx = {}
        self.my = {}
        self.mxy = {}
        self.bmx = {}
        self.bmy = {}
        self.bmxy = {}
        self.tx = {}
        self.ty = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.mx)
            time0 = get_key0(self.mx)
            nelements = len(self.mx[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.mx)
            msg.append('  type=%s nelements=%s\n'
                       % (self.__class__.__name__, nelements))
        msg.append('  mx, my, mxy, bmx, bmy, bmxy, tx, ty\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.mx[dt] = {}
        self.my[dt] = {}
        self.mxy[dt] = {}
        self.bmx[dt] = {}
        self.bmy[dt] = {}
        self.bmxy[dt] = {}
        self.tx[dt] = {}
        self.ty[dt] = {}

    def add_f06_data(self, dt, data):
        if dt is None:
            for (eid, grid, fx, fy, fxy, mx, my, mxy, qx, qy) in data:
                self.mx[eid] = mx

        else:
            if dt not in self.mx:
                pass
                #raise NotImplementedError()

    def add(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.mx[eid] = mx
        self.my[eid] = my
        self.mxy[eid] = mxy
        self.bmx[eid] = bmx
        self.bmy[eid] = bmy
        self.bmxy[eid] = bmxy
        self.tx[eid] = tx
        self.ty[eid] = ty

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def add_sort2(self, eid, dt, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            f.write('%s._write_f06_transient is not implemented\n' % self.__class__.__name__)
            return page_num
            #return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = [
            '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
            ' \n'
            '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
            '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
        ]
        #f.write('      ID      FX                            FY                            FXY                           MX                            MY                            MXY                           QX                            QY\n')
        f.write(''.join(words))
        for eid in self.mx:
            mx = self.mx[eid]
            my = self.my[eid]
            mxy = self.mxy[eid]

            bmx = self.bmx[eid]
            bmy = self.bmy[eid]
            bmxy = self.bmxy[eid]

            tx = self.tx[eid]
            ty = self.ty[eid]
            #Fmt = '% 8i   ' + '%27.20E   ' * 8 + '\n'

            [mx, my, mxy, bmx, bmy, bmxy, tx, ty] = write_floats_13e(
                [mx, my, mxy, bmx, bmy, bmxy, tx, ty])
            Fmt = '% 8i   ' + '%s ' * 8 + '\n'
            f.write(Fmt % (eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty))
        return page_num


class RealCBushForce(ScalarObject):  # 102-CBUSH
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.force = {}
        self.moment = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = get_key0(self.force)
            nelements = len(self.force[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.force)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force, moment\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.force[dt] = {}
        self.moment[dt] = {}

    def add(self, dt, eid, fx, fy, fz, mx, my, mz):
        self.force[eid] = array([fx, fy, fz], dtype='float32')
        self.moment[eid] = array([mx, my, mz], dtype='float32')

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = array([fx, fy, fz], dtype='float32')
        self.moment[dt][eid] = array([mx, my, mz], dtype='float32')

    def add_sort2(self, eid, dt, fx, fy, fz, mx, my, mz):
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = array([fx, fy, fz], dtype='float32')
        self.moment[dt][eid] = array([mx, my, mz], dtype='float32')


class RealSpringForce(ScalarObject):  # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.force = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = get_key0(self.force)
            nelements = len(self.force[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.force)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force\n')
        return msg

    def add_f06_data(self, data, dt):
        if dt is not None:
            for datai in data:
                (eid, forcei) = datai
                self.force[dt][eid] = forcei
            return

        for datai in data:
            (eid, forcei) = datai
            self.force[eid] = forcei

    def add_new_transient(self, dt):
        self.dt = dt
        self.force[dt] = {}

    def add(self, dt, data):
        [eid, force] = data
        self.force[eid] = force

    def add_sort1(self, dt, data):
        [eid, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def add_sort2(self, eid, data):
        [dt, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                 '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n',
                 '        ID.                              ID.                              ID.                              ID.\n',
                 ]
        msg = []
        for dt, Force in sorted(iteritems(self.force)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words

            forces = []
            #elements = []
            line = '   '
            for eid, force in sorted(iteritems(Force)):
                #elements.append(eid)
                forces.append(force)
                line += '%13s  %13s     ' % (eid, f)
                if len(forces) == 3:
                    forces = []
                    msg.append(line.rstrip() + '\n')
                    line = '   '

            if forces:
                msg.append(line.rstrip() + '\n')

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1

        return page_num - 1

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)
        msg = header + ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n',
                        '        ID.                              ID.                              ID.                              ID.\n',
                        ]
        f.write(''.join(msg))
        _write_f06_springs(f, self.force)
        f.write(page_stamp % page_num)

        return page_num


class RealDamperForce(ScalarObject):  # 20-CDAMP1,21-CDAMP2,22-CDAMP3,23-CDAMP4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.force = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = get_key0(self.force)
            nelements = len(self.force[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.force)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.force[dt] = {}

    def add(self, dt, data):
        [eid, force] = data
        self.force[eid] = force

    def add_sort1(self, dt, data):
        [eid, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def add_sort2(self, eid, data):
        [dt, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                 ' \n',
                 '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        msg = []
        for dt, Force in sorted(self.force.items()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words

            #packs = []
            forces = []
            elements = []
            line = '   '
            for eid, force in sorted(Force.items()):
                elements.append(eid)
                forces.append(force)
                #pack.append(eid)
                #pack.append(f)
                line += '%13s  %13s     ' % (eid, f)
                if len(forces) == 3:
                    msg.append(line.rstrip() + '\n')

            if forces:
                msg.append(line.rstrip() + '\n')
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1

        return page_num - 1

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)
        msg = header + ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        ' \n',
                        '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        #packs = []
        forces = []
        elements = []
        line = '   '
        for eid, force in sorted(self.force.items()):
            elements.append(eid)
            forces.append(force)
            #pack.append(eid)
            #pack.append(f)
            line += '%13s  %13s     ' % (eid, force)
            if len(forces) == 3:
                msg.append(line.rstrip() + '\n')

        if forces:
            msg.append(line.rstrip() + '\n')

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num


class RealViscForce(ScalarObject):  # 24-CVISC
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.axial_force = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = get_key0(self.torque)
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append(' axial_force, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.axial_force[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, axial_force, torque] = data
        self.axial_force[eid] = axial_force
        self.torque[eid] = torque

    def add_sort1(self, dt, data):
        [eid, axial_force, torque] = data
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, data):
        [dt, axial_force, torque] = data
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque


class RealCBarForce(ScalarObject):  # 34-CBAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.bendingMomentA = {}
        self.bendingMomentB = {}
        self.shear = {}
        self.axial = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = get_key0(self.torque)
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  bendingMomentA, bendingMomentB, shear, axial, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMomentA[dt] = {}
        self.bendingMomentB[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        self.bendingMomentA[eid] = [bm1a, bm2a]
        self.bendingMomentB[eid] = [bm1b, bm2b]
        self.shear[eid] = [ts1, ts2]
        self.axial[eid] = af
        self.torque[eid] = trq

    def add_sort1(self, dt, data):
        [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        if dt not in self.axial:
            self.add_new_transient(dt)

        self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def add_sort2(self, eid, data):
        [dt, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        if dt not in self.axial:
            self.add_new_transient(dt)

        self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq


    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=False, is_sort1=is_sort1)

        words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        f.write(''.join(words))
        for eid in sorted(self.bendingMomentA):
            bm1a, bm2a = self.bendingMomentA[eid]
            bm1b, bm2b = self.bendingMomentB[eid]
            ts1, ts2 = self.shear[eid]
            af = self.axial[eid]
            trq = self.torque[eid]
            vals2 = write_floats_13e([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq])
            [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = vals2
            f.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
            #1     2.504029E+06  9.728743E+06   5.088001E+05  1.976808E+06   1.995229E+06  7.751935E+06  -3.684978E-07  -1.180941E-07

            f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num, f, is_mag_phase=False, is_sort1=True):
        assert f is not None
        words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']

        for dt, bm in sorted(iteritems(self.bendingMomentA)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for eid in sorted(bm):
                bm1a, bm2a = self.bendingMomentA[dt][eid]
                bm1b, bm2b = self.bendingMomentB[dt][eid]
                ts1, ts2 = self.shear[dt][eid]
                af = self.axial[dt][eid]
                trq = self.torque[dt][eid]
                vals2 = write_floats_13e([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq])
                [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = vals2
                f.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
                #1     2.504029E+06  9.728743E+06   5.088001E+05  1.976808E+06   1.995229E+06  7.751935E+06  -3.684978E-07  -1.180941E-07

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCShearForce(ScalarObject):  # 4-CSHEAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.force41 = {}
        self.force14 = {}
        self.force21 = {}
        self.force12 = {}
        self.force32 = {}
        self.force23 = {}
        self.force43 = {}
        self.force34 = {}
        self.kickForce1 = {}
        self.kickForce2 = {}
        self.kickForce3 = {}
        self.kickForce4 = {}
        self.shear12 = {}
        self.shear23 = {}
        self.shear34 = {}
        self.shear41 = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.shear12)
            time0 = get_key0(self.shear12)
            nelements = len(self.shear12[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.shear12)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force41, force14, force21, force12, force32, force23, '
                   '  force43, force34, kickForce1, kickForce2, kickForce3, '
                   '  kickForce4, shear12, shear23, shear34, shear41\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.force41[dt] = {}
        self.force14[dt] = {}
        self.force21[dt] = {}
        self.force12[dt] = {}
        self.force32[dt] = {}
        self.force23[dt] = {}
        self.force43[dt] = {}
        self.force34[dt] = {}
        self.kickForce1[dt] = {}
        self.kickForce2[dt] = {}
        self.kickForce3[dt] = {}
        self.kickForce4[dt] = {}
        self.shear12[dt] = {}
        self.shear23[dt] = {}
        self.shear34[dt] = {}
        self.shear41[dt] = {}

    def add(self, dt, data):
        [eid, f41, f21, f12, f32, f23, f43, f34, f14,
         kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data
        self.force41[eid] = f41
        self.force14[eid] = f14
        self.force21[eid] = f21
        self.force12[eid] = f12
        self.force32[eid] = f32
        self.force23[eid] = f23
        self.force43[eid] = f43
        self.force34[eid] = f34
        self.kickForce1[eid] = kf1
        self.kickForce2[eid] = kf2
        self.kickForce3[eid] = kf3
        self.kickForce4[eid] = kf4
        self.shear12[eid] = s12
        self.shear23[eid] = s23
        self.shear34[eid] = s34
        self.shear41[eid] = s41

    def add_sort1(self, dt, data):
        [eid, f41, f21, f12, f32, f23, f43, f34, f14,
         kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data
        self._fill_object(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                          kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def add_sort2(self, eid, data):
        [dt, f41, f21, f12, f32, f23, f43, f34, f14,
            kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data

        self._fill_object(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                          kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def add_f06_data(self, data, dt=None):
        if dt:
            raise NotImplementedError(dt)

        for d in data:
            [
                eid,
                f41, f21, tau12, kick1,
                f12, f32, tau23, kick2,
                f23, f43, tau34, kick3,
                f34, f14, tau41, kick4,
            ] = d
            #print('eid, nid, sd', eid, nid, sd)
            self.force41[eid] = f41
            self.force14[eid] = f14

            self.force21[eid] = f21
            self.force12[eid] = f12

            self.force32[eid] = f32
            self.force23[eid] = f23

            self.force43[eid] = f43
            self.force34[eid] = f34

            self.kickForce1[eid] = kick1
            self.kickForce2[eid] = kick2
            self.kickForce3[eid] = kick3
            self.kickForce4[eid] = kick4

            self.shear12[eid] = tau12
            self.shear23[eid] = tau23
            self.shear34[eid] = tau34
            self.shear41[eid] = tau41

    def _fill_object(self, dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                    kf1, s12, kf2, s23, kf3, s34, kf4, s41):
        if dt not in self.force41:
            self.add_new_transient(dt)
        self.force41[dt][eid] = f41
        self.force14[dt][eid] = f14
        self.force21[dt][eid] = f21
        self.force12[dt][eid] = f12
        self.force32[dt][eid] = f32
        self.force23[dt][eid] = f23
        self.force43[dt][eid] = f43
        self.force34[dt][eid] = f34
        self.kickForce1[dt][eid] = kf1
        self.kickForce2[dt][eid] = kf2
        self.kickForce3[dt][eid] = kf3
        self.kickForce4[dt][eid] = kf4
        self.shear12[dt][eid] = s12
        self.shear23[dt][eid] = s23
        self.shear34[dt][eid] = s34
        self.shear41[dt][eid] = s41

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = [
            '                           F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n'
            ' \n'
            '                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n'
            '   ELEMENT        F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n'
            '         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n'
        ]
        msg = header + words
        for eid, forcei in sorted(iteritems(self.force14)):
            f41 = self.force41[eid]
            f14 = self.force14[eid]
            f21 = self.force21[eid]
            f12 = self.force12[eid]
            f32 = self.force32[eid]
            f23 = self.force23[eid]
            f43 = self.force43[eid]
            f34 = self.force34[eid]
            kick1 = self.kickForce1[eid]
            kick2 = self.kickForce2[eid]
            kick3 = self.kickForce3[eid]
            kick4 = self.kickForce4[eid]
            tau12 = self.shear12[eid]
            tau23 = self.shear23[eid]
            tau34 = self.shear34[eid]
            tau41 = self.shear41[eid]

            vals = [
                f14, f12,
                f21, f23,
                f32, f34,
                f43, f41,
                kick1, tau12,
                kick2, tau23,
                kick3, tau34,
                kick4, tau41,
            ]
            (vals2, is_all_zeros) = writeFloats12E(vals)
            [f14, f12,
             f21, f23,
             f32, f34,
             f43, f41,
             kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41,
            ] = vals2
            msg.append('0%13i%-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, f14, f12, f21, f23, f32, f34, f43, f41))
            msg.append('                     %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num




class RealCGapForce(ScalarObject):  # 38-CGAP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.fx = {}
        self.sfy = {}
        self.sfz = {}
        self.u = {}
        self.v = {}
        self.w = {}
        self.sv = {}
        self.sw = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.fx)
            time0 = get_key0(self.fx)
            nelements = len(self.fx[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.fx)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  fx, sfy, sfz, u, v, w, sv, sw\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.fx[dt] = {}
        self.sfy[dt] = {}
        self.sfz[dt] = {}
        self.u[dt] = {}
        self.v[dt] = {}
        self.w[dt] = {}
        self.sv[dt] = {}
        self.sw[dt] = {}

    def add(self, dt, eid, fx, sfy, sfz, u, v, w, sv, sw):
        self.fx[eid] = fx
        self.sfy[eid] = sfy
        self.sfz[eid] = sfz
        self.u[eid] = u
        self.v[eid] = v
        self.w[eid] = w
        self.sv[eid] = sv
        self.sw[eid] = sw

    def add_sort1(self, dt, eid, fx, sfy, sfz, u, v, w, sv, sw):
        if dt not in self.fx:
            self.add_new_transient(dt)
        self.fx[dt][eid] = fx
        self.sfy[dt][eid] = sfy
        self.sfz[dt][eid] = sfz
        self.u[dt][eid] = u
        self.v[dt][eid] = v
        self.w[dt][eid] = w
        self.sv[dt][eid] = sv
        self.sw[dt][eid] = sw

    def add_sort2(self, eid, dt, fx, sfy, sfz, u, v, w, sv, sw):
        if dt not in self.fx:
            self.add_new_transient(dt)
        self.fx[dt][eid] = fx
        self.sfy[dt][eid] = sfy
        self.sfz[dt][eid] = sfz
        self.u[dt][eid] = u
        self.v[dt][eid] = v
        self.w[dt][eid] = w
        self.sv[dt][eid] = sv
        self.sw[dt][eid] = sw


class RealConeAxForce(ScalarObject):  # 35-CCONEAX
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.hopa = {}
        self.bmu = {}
        self.bmv = {}
        self.tm = {}
        self.su = {}
        self.sv = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.hopa)
            time0 = get_key0(self.hopa)
            nelements = len(self.hopa[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.hopa)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  hopa, bmu, bmv, tm, su, sv\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.hopa[dt] = {}
        self.bmu[dt] = {}
        self.bmv[dt] = {}
        self.tm[dt] = {}
        self.su[dt] = {}
        self.sv[dt] = {}

    def add(self, dt, data):
        [eid, hopa, bmu, bmv, tm, su, sv] = data
        self.hopa[eid] = hopa
        self.bmu[eid] = bmu
        self.bmv[eid] = bmv
        self.tm[eid] = tm
        self.su[eid] = su
        self.sv[eid] = sv

    def add_sort1(self, dt, data):
        [eid, hopa, bmu, bmv, tm, su, sv] = data
        if dt not in self.hopa:
            self.add_new_transient(dt)
        self.hopa[dt][eid] = hopa
        self.bmu[dt][eid] = bmu
        self.bmv[dt][eid] = bmv
        self.tm[dt][eid] = tm
        self.su[dt][eid] = su
        self.sv[dt][eid] = sv

    def add_sort2(self, eid, data):
        [dt, hopa, bmu, bmv, tm, su, sv] = data
        if dt not in self.hopa:
            self.add_new_transient(dt)
        self.hopa[dt][eid] = hopa
        self.bmu[dt][eid] = bmu
        self.bmv[dt][eid] = bmv
        self.tm[dt][eid] = tm
        self.su[dt][eid] = su
        self.sv[dt][eid] = sv


class RealPentaPressureForce(ScalarObject):  # 77-PENTA_PR,78-TETRA_PR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.acceleration = {}
        self.velocity = {}
        self.pressure = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.acceleration)
            time0 = get_key0(self.acceleration)
            nelements = len(self.acceleration[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.acceleration)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  acceleration, velocity, pressure\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.acceleration[dt] = {}
        self.velocity[dt] = {}
        self.pressure[dt] = {}

    def add(self, dt, data):
        [eid, ename, ax, ay, az, vx, vy, vz, pressure] = data
        self.acceleration[eid] = [ax, ay, az]
        self.velocity[eid] = [vx, vy, vz]
        self.pressure[eid] = pressure

    def add_sort1(self, dt, data):
        [eid, ename, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure

    def add_sort2(self, eid, data):
        [dt, ename, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #words = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
        #         ' \n',
        #         '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)
        f.write('%s write_f06 not implemented...\n' % self.__class__.__name__)
        #raise NotImplementedError()
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                 ' \n',
                 '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        msg = []
        aaa
        for dt, acc in sorted(iteritems(self.acceleration)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid in sorted(acc):
                ax, ay, az = self.acceleration[dt][eid]
                vx, vy, vz = self.velocity[dt][eid]
                pressure = self.pressure[dt][eid]
                vals = [ax, ay, az, pressure]
                vals2 = write_floats_13e(vals)
                [ax, ay, az, pressure] = vals2
                etype = 'PENPR'
                msg.append('0%13s    %5s               %-13s             %-13s             %-13s             %s\n' % (eid, etype, ax, ay, az, pressure))
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
