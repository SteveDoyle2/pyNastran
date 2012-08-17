from pyNastran.op2.resultObjects.op2_Objects import scalarObject


class GridPointStressesObject(scalarObject):
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.nx = {}
        self.ny = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.tmax = {}
        self.ovm = {}

        self.elemName = {}
        self.eids = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def addNewTransient(self, dt):  # eKey
        """initializes the transient variables"""
        self.nx[dt] = {}
        self.ny[dt] = {}
        self.txy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.tmax[dt] = {}
        self.ovm[dt] = {}

        self.elemName = {}
        self.eids = {}

    def add(self, dt, eKey, eid, elemName, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        if eKey not in self.nx:
            self.eids[eKey] = []
            self.elemName[eKey] = []
            self.nx[eKey] = []
            self.ny[eKey] = []
            self.txy[eKey] = []
            self.angle[eKey] = []
            self.majorP[eKey] = []
            self.minorP[eKey] = []
            self.tmax[eKey] = []
            self.ovm[eKey] = []
        self.nx[eKey].append(nx)
        self.ny[eKey].append(ny)
        self.txy[eKey].append(txy)
        self.angle[eKey].append(angle)
        self.majorP[eKey].append(majorP)
        self.minorP[eKey].append(minorP)
        self.tmax[eKey].append(tmax)
        self.ovm[eKey].append(ovm)

        self.elemName[eKey].append(elemName)
        self.eids[eKey].append(eid)

    def addSort1(self, dt, eKey, eid, elemName, nx, ny, txy, angle, majorP, minorP, tmax, ovm):
        if dt not in self.nx:
            self.addNewTransient(dt)

        #print "%s=%s eKey=%s eid=%s elemName=%s f1=%s" %(self.dataCode['name'],dt,eKey,eid,elemName,f1)
        if eKey not in self.nx[dt]:
            self.eids[eKey] = []
            self.elemName[eKey] = []
            self.nx[dt][eKey] = []
            self.ny[dt][eKey] = []
            self.txy[dt][eKey] = []
            self.angle[dt][eKey] = []
            self.majorP[dt][eKey] = []
            self.minorP[dt][eKey] = []
            self.tmax[dt][eKey] = []
            self.ovm[dt][eKey] = []
        self.eids[eKey].append(eid)
        self.elemName[eKey].append(elemName)

        self.nx[dt][eKey].append(nx)
        self.ny[dt][eKey].append(ny)
        self.txy[dt][eKey].append(txy)
        self.angle[dt][eKey].append(angle)
        self.majorP[dt][eKey].append(majorP)
        self.minorP[dt][eKey].append(minorP)
        self.tmax[dt][eKey].append(tmax)
        self.ovm[dt][eKey].append(ovm)

    def deleteTransient(self, dt):
        del self.nx[dt]
        del self.ny[dt]
        del self.txy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.tmax[dt]
        del self.ovm[dt]

    def getTransients(self):
        k = self.nx.keys()
        k.sort()
        return k

    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for eKey, nxs in sorted(self.nx.iteritems()):
            eKey2 = eKey
            zero = '0'
            for iLoad, nx in enumerate(nxs):
                ny = self.ny[eKey][iLoad]
                txy = self.txy[eKey][iLoad]
                angle = self.angle[eKey][iLoad]
                majorP = self.majorP[eKey][iLoad]
                minorP = self.minorP[eKey][iLoad]
                tmax = self.tmax[eKey][iLoad]
                ovm = self.ovm[eKey][iLoad]

                (elemName) = self.elemName[eKey][iLoad]
                eid = self.eids[eKey][iLoad]
                vals = [nx, ny, txy, majorP, minorP, tmax, ovm]
                (vals2, isAllZeros) = self.writeFloats10E(vals)
                [nx, ny, txy, majorP, minorP, tmax, ovm] = vals2
                if eid == 0:
                    eid = zero
                msg.append('%s%8s  %8s   %4s    %s %s %s   %7.4f %s %s %s  %-s\n' % (zero, eKey2, eid, elemName, nx, ny, txy, angle, majorP, minorP, tmax, ovm.rstrip()))
                zero = ' '
                eKey2 = ' '
            ###
        ###
        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'GridPointStressesObject writeF06 is not implemented...', pageNum
        #raise NotImplementedError()
        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for dt, Forces in sorted(self.forces.iteritems()):
            for eKey, force in sorted(Forces.iteritems()):
                zero = '0'
                for iLoad, f in enumerate(force):
                    (f1, f2, f3) = f
                    (m1, m2, m3) = self.moments[dt][eKey][iLoad]
                    (elemName) = self.elemName[eKey][iLoad]
                    eid = self.eids[eKey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    (vals2, isAllZeros) = self.writeFloats13E(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''

                    msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' % (zero, eKey, eid, elemName, f1, f2, f3, m1, m2, m3))
                    zero = ' '

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        return self.writeF06([], 'PAGE ', 1)[0]
        #return '---gridPointStressesObject---'


class GridPointStressesVolumeObject(scalarObject):
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.nx = {}
        self.ny = {}
        self.nz = {}
        self.txy = {}
        self.tyz = {}
        self.txz = {}
        self.pressure = {}
        self.ovm = {}

        self.elemName = {}
        self.eids = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def addNewTransient(self, dt):  # eKey
        """initializes the transient variables"""
        self.nx[dt] = {}
        self.ny[dt] = {}
        self.nz[dt] = {}
        self.txy[dt] = {}
        self.tyz[dt] = {}
        self.txz[dt] = {}
        self.pressure[dt] = {}
        self.ovm[dt] = {}

        self.elemName = {}
        self.eids = {}

    def add(self, dt, eKey, nx, ny, nz, txy, tyz, txz, pressure, ovm):
        if eKey not in self.nx:
            #self.eids[eKey] = []
            #self.elemName[eKey] = []
            self.nx[eKey] = []
            self.ny[eKey] = []
            self.nz[eKey] = []
            self.txy[eKey] = []
            self.tyz[eKey] = []
            self.txz[eKey] = []
            self.pressure[eKey] = []
            self.ovm[eKey] = []
        self.nx[eKey].append(nx)
        self.ny[eKey].append(ny)
        self.nz[eKey].append(nz)
        self.txy[eKey].append(txy)
        self.tyz[eKey].append(tyz)
        self.txz[eKey].append(txz)
        self.pressure[eKey].append(pressure)
        self.ovm[eKey].append(ovm)

        #self.elemName[eKey].append(elemName)
        #self.eids[eKey].append(eid)

    def addSort1(self, dt, eKey, nx, ny, nz, txy, tyz, txz, pressure, ovm):
        if dt not in self.nx:
            self.addNewTransient(dt)

        #print "%s=%s eKey=%s eid=%s elemName=%s f1=%s" %(self.dataCode['name'],dt,eKey,eid,elemName,f1)
        if eKey not in self.nx[dt]:
            #self.eids[eKey] = []
            #self.elemName[eKey] = []
            self.nx[dt][eKey] = []
            self.ny[dt][eKey] = []
            self.nz[dt][eKey] = []
            self.txy[dt][eKey] = []
            self.tyz[dt][eKey] = []
            self.txz[dt][eKey] = []
            self.pressure[eKey] = []
            self.ovm[dt][eKey] = []
        self.eids[eKey].append(eid)
        #self.elemName[eKey].append(elemName)

        self.nx[dt][eKey].append(nx)
        self.ny[dt][eKey].append(ny)
        self.nz[dt][eKey].append(nz)
        self.txy[dt][eKey].append(txy)
        self.tyz[dt][eKey].append(tyz)
        self.txz[dt][eKey].append(txz)
        self.pressure[dt][eKey].append(pressure)
        self.ovm[dt][eKey].append(ovm)

    def deleteTransient(self, dt):
        del self.nx[dt]
        del self.ny[dt]
        del self.nz[dt]
        del self.txy[dt]
        del self.tyz[dt]
        del self.txz[dt]
        del self.pressure[dt]
        del self.ovm[dt]

    def getTransients(self):
        k = self.nx.keys()
        k.sort()
        return k

    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'GridPointStressesVolumeObject writeF06 is not implemented...', pageNum
        #raise NotImplementedError()
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for eKey, nxs in sorted(self.nx.iteritems()):
            eKey2 = eKey
            zero = '0'
            for iLoad, nx in enumerate(nxs):
                ny = self.ny[eKey][iLoad]
                nz = self.nz[eKey][iLoad]
                txy = self.txy[eKey][iLoad]
                tyz = self.tyz[eKey][iLoad]
                txz = self.txz[eKey][iLoad]
                pressure = self.pressure[eKey][iLoad]
                ovm = self.ovm[eKey][iLoad]

                #(elemName) = self.elemName[eKey][iLoad]
                #eid = self.eids[eKey][iLoad]
                vals = [nx, ny, nz, txy, tyz, txz, pressure, ovm]
                (vals2, isAllZeros) = self.writeFloats10E(vals)
                [nx, ny, nz, txy, tyz, txz, pressure, ovm] = vals2
                msg.append('%s%8s  %s %s %s   %s %s %s %s  %-s\n' % (zero, eKey, nx, ny, nz, txy, tyz, txz, pressure, ovm.rstrip()))
                zero = ' '
                eKey2 = ' '
            ###
        ###
        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'GridPointStressesVolumeObject writeF06Transient is not implemented...', pageNum
        #raise NotImplementedError()
        msg = header + ['                                  S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5\n',
                        '0                       SURFACE X-AXIS X  NORMAL(Z-AXIS)  Z         REFERENCE COORDINATE SYSTEM FOR SURFACE DEFINITION CID        0\n',
                        '     GRID      ELEMENT            STRESSES IN SURFACE SYSTEM           PRINCIPAL STRESSES            MAX             \n',
                        '     ID          ID    FIBRE   NORMAL-X   NORMAL-Y   SHEAR-XY     ANGLE      MAJOR      MINOR      SHEAR     VON MISES\n']
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for dt, Forces in sorted(self.forces.iteritems()):
            for eKey, force in sorted(Forces.iteritems()):
                zero = '0'
                for iLoad, f in enumerate(force):
                    (f1, f2, f3) = f
                    (m1, m2, m3) = self.moments[dt][eKey][iLoad]
                    (elemName) = self.elemName[eKey][iLoad]
                    eid = self.eids[eKey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    (vals2, isAllZeros) = self.writeFloats13E(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''

                    msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' % (zero, eKey, eid, elemName, f1, f2, f3, m1, m2, m3))
                    zero = ' '

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        return self.writeF06([], 'PAGE ', 1)[0]
        #return '---gridPointStressesVolumeObject---'
