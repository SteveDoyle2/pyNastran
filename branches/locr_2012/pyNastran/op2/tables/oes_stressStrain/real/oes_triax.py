from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys

from .oes_objects import stressObject, strainObject


class TriaxStressObject(stressObject):
    """
    @code
    # formatCode=1 sortCode=0 stressCode=0
                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S
    ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES
       ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
       5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CTRIAX6'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.radial = {}
        self.azimuthal = {}
        self.axial = {}
        self.shear = {}
        self.omax = {}
        self.oms = {}
        self.ovm = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def addF06Data(self, data, transient):
        raise Exception('Not Implemented')
        if transient is None:
            for line in data:
                (eid, axial, MSa, torsion, MSt) = line
                if MSa is None:
                    MSa = 0.
                if MSt is None:
                    MSt = 0.
                self.axial[eid] = axial
                self.MS_axial[eid] = MSa
                self.torsion[eid] = torsion
                self.MS_torsion[eid] = MSt
            return

        (dtName, dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.s1:
            self.updateDt(self.dataCode, dt)
            self.isTransient = True

        for line in data:
            (eid, axial, MSa, torsion, MSt) = line
            if MSa is None:
                MSa = 0.
            if MSt is None:
                MSt = 0.
            self.axial[dt][eid] = axial
            self.MS_axial[dt][eid] = MSa
            self.torsion[dt][eid] = torsion
            self.MS_torsion[dt][eid] = MSt

    def deleteTransient(self, dt):
        del self.radial[dt]
        del self.azimuthal[dt]
        del self.axial[dt]
        del self.shear[dt]
        del self.omax[dt]
        del self.oms[dt]
        del self.ovm[dt]

    def getTransients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.radial[dt] = {}
        self.azimuthal[dt] = {}
        self.axial[dt] = {}
        self.shear[dt] = {}
        self.omax[dt] = {}
        self.oms[dt] = {}
        self.ovm[dt] = {}

    def addNewEid(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid] = {nid: rs}
        self.azimuthal[eid] = {nid: azs}
        self.axial[eid] = {nid: As}
        self.shear[eid] = {nid: ss}
        self.omax[eid] = {nid: maxp}
        self.oms[eid] = {nid: tmax}
        self.ovm[eid] = {nid: octs}

    def add(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid][nid] = rs
        self.azimuthal[eid][nid] = azs
        self.axial[eid][nid] = As
        self.shear[eid][nid] = ss
        self.omax[eid][nid] = maxp
        self.oms[eid][nid] = tmax
        self.ovm[eid][nid] = octs

    def addNewEidSort1(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #assert isinstance(eid,int)
        #assert eid >= 0
        #print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        if dt not in self.radial:
            self.addNewTransient(dt)
        self.radial[dt][eid] = {nid: rs}
        self.azimuthal[dt][eid] = {nid: azs}
        self.axial[dt][eid] = {nid: As}
        self.shear[dt][eid] = {nid: ss}
        self.omax[dt][eid] = {nid: maxp}
        self.oms[dt][eid] = {nid: tmax}
        self.ovm[dt][eid] = {nid: octs}

    def addSort1(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid][nid] = rs
        self.azimuthal[dt][eid][nid] = azs
        self.axial[dt][eid][nid] = As
        self.shear[dt][eid][nid] = ss
        self.omax[dt][eid][nid] = maxp
        self.oms[dt][eid][nid] = tmax
        self.ovm[dt][eid][nid] = octs

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        msg = header + ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
                        '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                        '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #out = []
        for eid, radial in sorted(self.radial.iteritems()):
            for nid in sorted(radial):
                rad = self.radial[eid][nid]
                azimuth = self.azimuthal[eid][nid]
                axial = self.axial[eid][nid]
                shear = self.shear[eid][nid]
                omax = self.omax[eid][nid]
                oms = self.oms[eid][nid]
                ovm = self.ovm[eid][nid]
                if nid == 0:
                    Eid = eid
                else:
                    Eid = ''
                ([rad, azimuth, axial, shear, omax, oms, ovm], isAllZeros) = self.writeFloats13E([rad, azimuth, axial, shear, omax, oms, ovm])
                msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, radial, azimuth, axial, shear, omax, oms, ovm.rstrip()))
            msg.append('\n')

        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return(''.join(msg), pageNum)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
                 '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                 '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        msg = []
        for dt, Radial in sorted(self.radial.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.dataCode['name'], dt)
            msg += header + words
            for eid, radial in sorted(Radial.iteritems()):
                for nid in sorted(radial):
                    rad = self.radial[dt][eid][nid]
                    azimuth = self.azimuthal[dt][eid][nid]
                    axial = self.axial[dt][eid][nid]
                    shear = self.shear[dt][eid][nid]
                    omax = self.omax[dt][eid][nid]
                    oms = self.oms[dt][eid][nid]
                    ovm = self.ovm[dt][eid][nid]
                    if nid == 0:
                        Eid = eid
                    else:
                        Eid = ''
                    ([rad, azimuth, axial, shear, omax, oms, ovm], isAllZeros) = self.writeFloats13E([rad, azimuth, axial, shear, omax, oms, ovm])
                    msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, rad, azimuth, axial, shear, omax, oms, ovm.rstrip()))
                msg.append('\n')

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return(''.join(msg), pageNum - 1)

    def __repr__(self):
        return self.writeF06(['', ''], 'PAGE ', 1)[0]
        if self.nonlinearFactor is not None:
            pass


class TriaxStrainObject(strainObject):
    """
    @code
    # formatCode=1 sortCode=0 stressCode=0
                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S
    ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES
       ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
       5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        strainObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CTRIAX6'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.radial = {}
        self.azimuthal = {}
        self.axial = {}
        self.shear = {}
        self.emax = {}
        self.ems = {}
        self.evm = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def addF06Data(self, data, transient):
        raise Exception('Not Implemented')

    def deleteTransient(self, dt):
        del self.radial[dt]
        del self.azimuthal[dt]
        del self.axial[dt]
        del self.shear[dt]
        del self.emax[dt]
        del self.ems[dt]
        del self.evm[dt]

    def getTransients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.radial[dt] = {}
        self.azimuthal[dt] = {}
        self.axial[dt] = {}
        self.shear[dt] = {}
        self.emax[dt] = {}
        self.ems[dt] = {}
        self.evm[dt] = {}

    def addNewEid(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid] = {nid: rs}
        self.azimuthal[eid] = {nid: azs}
        self.axial[eid] = {nid: As}
        self.shear[eid] = {nid: ss}
        self.emax[eid] = {nid: maxp}
        self.ems[eid] = {nid: emax}
        self.evm[eid] = {nid: ects}

    def add(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid][nid] = rs
        self.azimuthal[eid][nid] = azs
        self.axial[eid][nid] = As
        self.shear[eid][nid] = ss
        self.emax[eid][nid] = maxp
        self.ems[eid][nid] = emax
        self.evm[eid][nid] = ects

    def addNewEidSort1(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        #assert isinstance(eid,int)
        #assert eid >= 0
        #print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid] = {nid: rs}
        self.azimuthal[dt][eid] = {nid: azs}
        self.axial[dt][eid] = {nid: As}
        self.shear[dt][eid] = {nid: ss}
        self.emax[dt][eid] = {nid: maxp}
        self.ems[dt][eid] = {nid: emax}
        self.evm[dt][eid] = {nid: ects}

    def addSort1(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid][nid] = rs
        self.azimuthal[dt][eid][nid] = azs
        self.axial[dt][eid][nid] = As
        self.shear[dt][eid][nid] = ss
        self.emax[dt][eid][nid] = maxp
        self.ems[dt][eid][nid] = emax
        self.evm[dt][eid][nid] = ects

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        msg = header + ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
                        '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                        '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #out = []
        for eid, radial in sorted(self.radial.iteritems()):
            for nid in sorted(radial):
                rad = self.radial[eid][nid]
                azimuth = self.azimuthal[eid][nid]
                axial = self.axial[eid][nid]
                shear = self.shear[eid][nid]
                emax = self.emax[eid][nid]
                ems = self.ems[eid][nid]
                evm = self.evm[eid][nid]
                if nid == 0:
                    Eid = eid
                else:
                    Eid = ''
                ([rad, azimuth, axial, shear, emax, ems, evm], isAllZeros) = self.writeFloats13E([rad, azimuth, axial, shear, emax, ems, evm])
                msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, radial, azimuth, axial, shear, emax, ems, evm.rstrip()))
            msg.append('\n')

        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return(''.join(msg), pageNum)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
                 '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                 '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        msg = []
        for dt, Radial in sorted(self.radial.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.dataCode['name'], dt)
            msg += header + words
            for eid, radial in sorted(Radial.iteritems()):
                for nid in sorted(radial):
                    rad = self.radial[dt][eid][nid]
                    azimuth = self.azimuthal[dt][eid][nid]
                    axial = self.axial[dt][eid][nid]
                    shear = self.shear[dt][eid][nid]
                    emax = self.emax[dt][eid][nid]
                    ems = self.ems[dt][eid][nid]
                    evm = self.evm[dt][eid][nid]
                    if nid == 0:
                        Eid = eid
                    else:
                        Eid = ''
                    ([rad, azimuth, axial, shear, emax, ems, evm], isAllZeros) = self.writeFloats13E([rad, azimuth, axial, shear, emax, ems, evm])
                    msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, rad, azimuth, axial, shear, emax, ems, evm.rstrip()))
                msg.append('\n')

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return(''.join(msg), pageNum - 1)

    def __repr__(self):
        return self.writeF06(['', ''], 'PAGE ', 1)[0]
        if self.nonlinearFactor is not None:
            pass
