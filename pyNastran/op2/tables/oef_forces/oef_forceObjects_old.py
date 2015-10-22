class RealRodForce(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)

        self.axialForce = {}
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
        self.axialForce[dt] = {}
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
        msg.append('  axialForce, torque\n')
        return msg

    def add(self, dt, eid, axialForce, torque):
        self.axialForce[eid] = axialForce
        self.torque[eid] = torque

    def add_sort1(self, dt, eid, axialForce, torque):
        if dt not in self.axialForce:
            self.add_new_transient(dt)
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, data):
        [dt, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)

        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torque) = line
                self.axialForce[eid] = axial
                self.torque[eid] = torque
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        #if dt not in self.axialForce:
            #self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axialForce[dt][eid] = axial
            self.torque[dt][eid] = torsion

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        msg = []
        itime = 0
        ntimes = len(self.axialForce)
        for dt, axials in sorted(iteritems(self.axialForce)):
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            #dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            #header[2] = dtLine
            msg += header
            out = []
            for eid in sorted(axials):
                axial = self.axialForce[dt][eid]
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
        for eid in sorted(self.axialForce):
            axial = self.axialForce[eid]
            torsion = self.torque[eid]
            (vals2, is_all_zeros) = writeFloats13E([axial, torsion])
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
                (vals2, is_all_zeros) = writeFloats13E([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
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

            ([mx, my, mxy, bmx, bmy, bmxy, tx, ty],
             is_all_zeros) = writeFloats13E([mx, my, mxy, bmx, bmy, bmxy, tx, ty])
            Fmt = '% 8i   ' + '%s ' * 8 + '\n'
            f.write(Fmt % (eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty))
        return page_num

