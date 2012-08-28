from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys

from ..real.oes_objects import stressObject, strainObject


class ComplexPlateStressObject(stressObject):
    """
    @code
                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )
                                                          (REAL/IMAGINARY)

    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -
      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY
0       100    CEN/8  -2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0
                       2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        #print "making complex plate stress obj"
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = {}

        #self.appendDataMember('sCodes','sCode')
        #print "self.sCode = ",self.sCode
        self.code = [self.formatCode, self.sortCode, self.sCode]

        self.fiberCurvature = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
                self.addNewNode = self.addNewNodeSort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2
            self.addNewNode = self.addNewNodeSort2

    def addF06Data(self, data, transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ox1, oy1, txy1,
                     f2, ox2, oy2, txy2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'C': [f1, f2]}
                    self.oxx[eid] = {'C': [ox1, ox2]}
                    self.oyy[eid] = {'C': [oy1, oy2]}
                    self.txy[eid] = {'C': [txy1, txy2]}
                elif eType == 'CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, nid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        if nid == 'CEN/4':
                            nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                    elif len(line) == 17:  # Bilinear
                        #print line
                        (nid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.oxx[eid][nid] = [ox1, ox2]
                        self.oyy[eid][nid] = [oy1, oy2]
                        self.txy[eid][nid] = [txy1, txy2]
                    else:
                        assert len(line) == 19, 'len(line)=%s' % (len(line))
                        raise NotImplementedError()
                else:
                    raise NotImplementedError('line=%s not supported...' %
                                              (line))
            return
        #for line in data:
            #print line
        raise NotImplementedError('transient results not supported')

    def deleteTransient(self, dt):
        #del self.fiberCurvature[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]

    def getTransients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        #self.fiberCurvature[dt] = {}
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.txy[dt] = {}

    def addNewEid(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy):
        #print "***addNewEid..."
        #if eid in self.oxx:
            #return self.add(dt,eid,nodeID,fd,oxx,oyy,txy)

        #assert eid not in self.oxx
        #print self.oxx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [fd]}
        assert isinstance(eid, int)
        self.oxx[eid] = {nodeID: [oxx]}
        self.oyy[eid] = {nodeID: [oyy]}
        self.txy[eid] = {nodeID: [txy]}
        msg = "eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            eid, nodeID, fd, oxx, oyy, txy)
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def addNewEidSort1(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy):
        #print "***addNewEidSort1..."
        #msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" %(dt,eid,nodeID,fd,oxx,oyy,txy)
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%s" % (
            dt, eid, nodeID, fd, oxx)
        #print msg
        #if eid in self.oxx[dt]:
        #    return self.add(eid,nodeID,fd,oxx,oyy,txy)

        if dt in self.oxx and eid in self.oxx[dt]:  # SOL200, erase the old result
            nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid])))
            self.deleteTransient(dt)
            self.addNewTransient(dt)

        assert isinstance(eid, int)
        self.eType[eid] = eType
        if dt not in self.oxx:
            self.addNewTransient(dt)
        self.fiberCurvature[eid] = {nodeID: [fd]}
        self.oxx[dt][eid] = {nodeID: [oxx]}
        self.oyy[dt][eid] = {nodeID: [oyy]}
        self.txy[dt][eid] = {nodeID: [txy]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add(self, dt, eid, nodeID, fd, oxx, oyy, txy):
        #print "***add"
        msg = "eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            eid, nodeID, fd, oxx, oyy, txy)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        assert isinstance(eid, int)
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[eid][nodeID].append(oxx)
        self.oyy[eid][nodeID].append(oyy)
        self.txy[eid][nodeID].append(txy)
        if nodeID == 0:
            raise Exception(msg)

    def addSort1(self, dt, eid, nodeID, fd, oxx, oyy, txy):
        #print "***addSort1"
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            dt, eid, nodeID, fd, oxx, oyy, txy)
        #print msg
        #print self.oxx
        #print self.fiberCurvatrure
        assert eid is not None
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[dt][eid][nodeID].append(oxx)
        self.oyy[dt][eid][nodeID].append(oyy)
        self.txy[dt][eid][nodeID].append(txy)
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNode(self, dt, eid, nodeID, fd, oxx, oyy, txy):
        #print "***addNewNode"
        #print self.oxx
        assert eid is not None
        #assert nodeID not in self.oxx[eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[eid][nodeID] = [oxx]
        self.oyy[eid][nodeID] = [oyy]
        self.txy[eid][nodeID] = [txy]
        msg = "eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            eid, nodeID, fd, oxx, oyy, txy)
        #print msg
        if nodeID == 0:
            raise Exception(msg)

    def addNewNodeSort1(self, dt, eid, nodeID, fd, oxx, oyy, txy):
        #print "***addNewNodeTransient_Sort1"
        #print self.oxx
        assert eid is not None
        msg = "eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            eid, nodeID, fd, oxx, oyy, txy)
        #print msg
        #assert nodeID not in self.oxx[dt][eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[dt][eid][nodeID] = [oxx]
        self.oyy[dt][eid][nodeID] = [oyy]
        self.txy[dt][eid][nodeID] = [txy]
        if nodeID == 0:
            raise ValueError(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
        headers += ['oxx', 'oyy', 'txy']
        if self.isVonMises():
            headers.append('oVonMises')
        else:
            headers.append('maxShear')
        return headers

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        #print self.oxx.keys()
        for dt, oxxs in sorted(self.oxx.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid, oxxNodes in sorted(oxxs.iteritems()):
                eType = self.eType[eid]
                k = oxxNodes.keys()
                k.remove('C')
                k.sort()
                for nid in ['C'] + k:
                    for iLayer in xrange(len(self.oxx[dt][eid][nid])):
                        fd = self.fiberCurvature[eid][nid][iLayer]
                        oxx = self.oxx[dt][eid][nid][iLayer]
                        oyy = self.oyy[dt][eid][nid][iLayer]
                        txy = self.txy[dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s %10g ' % (eid, eType,
                                                           nid, iLayer + 1, fd)
                        vals = [oxx, oyy, txy]
                        for val in vals:
                            if abs(val) < 1e-6:
                                msg += '%10s ' % ('0')
                            else:
                                try:
                                    msg += '%10i %10i' % (int(val.real), int(val.imag))
                                except:
                                    print("bad val = %s" % (val))
                                    raise
                            ###
                        msg += '\n'
                    ###
                ###
            ###
        ###
        return msg

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        assert f is not None
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f, isMagPhase)

        #if self.isVonMises():
        #    vonMises = 'VON MISES'
        #else:
        #    vonMises = 'MAX SHEAR'

        if isMagPhase:
            pass
        else:
            magReal = ['                                                          (REAL/IMAGINARY)\n', ' \n']

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -',
                           '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY']
        else:
            pass

#'0       100    CEN/8  -2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'
#'                       2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'
        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None

        quadMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'] + formWord + quadMsgTemp
        quadrMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )'] + formWord + quadMsgTemp
        quad8Msg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )'] + formWord + quadMsgTemp

        eTypes = self.eType.values()
        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                eids.sort()
                #print "eType = ",eType
                #print "eids = ",eids
                #print "eType = ",eType
                msgPack = msgPacks[eType]

                msg += header + msgPack
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for eid in eids:
                            out = self.writeF06_Quad4_Bilinear(eid, 4)
                            msg.append(out)
                    else:
                        for eid in eids:
                            out = self.writeF06_Tri3(eid)
                            msg.append(out)
                    ###
                elif eType in ['CTRIA3']:
                    for eid in eids:
                        out = self.writeF06_Tri3(eid)
                        msg.append(out)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 8)
                        msg.append(out)
                elif eType in ['CTRIAR', 'CTRIA6']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 3)
                        msg.append(out)
                else:
                    raise NotImplementedError('eType = |%r|' %
                                              (eType))  # CQUAD8, CTRIA6
                msg.append(pageStamp + str(pageNum) + '\n')
                if f is not None:
                    f.write(''.join(msg))
                    msg = ['']
                pageNum += 1
            ###
        ###
        return (''.join(msg), pageNum - 1)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if isMagPhase:
            magReal = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            magReal = ['                                                          (REAL/IMAGINARY)\n \n']

        if self.isFiberDistance():
            gridMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                            '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']

        else:
            pass

        #if self.isFiberDistance():
            #quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
            #               '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            #triMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
            #              '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
        #else:
            #quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
            #               '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            #triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
            #              '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        #quadMsg  = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )']+magReal+fiberMsgTemp
        #quadrMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )']+magReal+fiberMsgTemp
        #quad8Msg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )']+magReal+fiberMsgTemp
        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None

        eTypes = self.eType.values()
        dts = self.oxx.keys()
        #print self.oxx
        #print "dts = ",dts
        dt = dts[0]
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            #print "qkey=%s kkey=%s" %(qkey,kkey)
            ekey = self.oxx[dt][kkey].keys()
            isBilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + magReal + gridMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + magReal + fiberMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + magReal + gridMsgTemp

        if 'CQUADR' in eTypes:
            quadrMsg = header + ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + magReal + gridMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + magReal + fiberMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + magReal + gridMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + magReal + gridMsgTemp

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        dts = self.oxx.keys()
        dts.sort()
        for eType in typesOut:
            #print "eType = ",eType
            eids = orderedETypes[eType]
            #print "eids = ",eids
            if eids:
                msgPack = msgPacks[eType]
                eids.sort()
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.dataCode[
                                'name'], dt)
                            msg += header + msgPack
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt, eid, 4, isMagPhase)
                                msg.append(out)
                            msg.append(pageStamp + str(pageNum) + '\n')
                            pageNum += 1
                            if f is not None:
                                f.write(''.join(msg))
                                msg = ['']
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.dataCode[
                                'name'], dt)
                            msg += header + msgPack
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid, isMagPhase)
                                msg.append(out)
                            msg.append(pageStamp + str(pageNum) + '\n')
                            pageNum += 1
                            if f is not None:
                                f.write(''.join(msg))
                                msg = ['']
                    ###
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.dataCode['name'], dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self.writeF06_Tri3Transient(dt,
                                                              eid, isMagPhase)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        pageNum += 1
                        if f is not None:
                            f.write(''.join(msg))
                            msg = ['']
                elif eType in ['CTRIA6', 'CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.dataCode['name'], dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 3, isMagPhase)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        pageNum += 1
                        if f is not None:
                            f.write(''.join(msg))
                            msg = ['']
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.dataCode['name'], dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 8, isMagPhase)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        pageNum += 1
                        if f is not None:
                            f.write(''.join(msg))
                            msg = ['']
                else:
                    raise NotImplementedError('eType = |%r|' %
                                              (eType))  # CQUAD8, CTRIA6
                if f is not None:
                    f.write(''.join(msg))
                    msg = ['']
            ### eids
        ### all eTypes
        return (''.join(msg), pageNum - 1)

    def writeF06_Quad4_Bilinear(self, eid, n, isMagPhase):
        msg = ''
        k = self.oxx[eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[eid][nid][iLayer]
                oyy = self.oyy[eid][nid][iLayer]
                txy = self.txy[eid][nid][iLayer]
                ([fd, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], isAllZeros) = self.writeImagFloats13E([fd, oxx, oyy, txy], isMagPhase)

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % (eid, 'CEN/' + str(n), fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % ('', nid, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n\n' % ('', '', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi.rstrip())
                else:
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Quad4_BilinearTransient(self, dt, eid, n, isMagPhase):
        msg = ''
        k = self.oxx[dt][eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.oxx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[dt][eid][nid][iLayer]
                oyy = self.oyy[dt][eid][nid][iLayer]
                txy = self.txy[dt][eid][nid][iLayer]
                ([fd, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], isAllZeros) = self.writeImagFloats13E([fd, oxx, oyy, txy], isMagPhase)

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % (eid, 'CEN/' + str(n), fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % ('', nid, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n\n' % ('', '', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi.rstrip())
                else:
                    #msg += '   %8s %8s  %13E  %13E %13E %13E   %8.4F  %13E %13E %13E\n' %('','',  fd,oxx,oyy,txy)
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Tri3(self, eid, isMagPhase):
        msg = ''
        k = self.oxx[eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[eid][nid][iLayer]
                oyy = self.oyy[eid][nid][iLayer]
                txy = self.txy[eid][nid][iLayer]
                ([fd, oxx, oyy, txy], isAllZeros) = self.writeFloats13E(
                    [fd, oxx, oyy, txy])

                if iLayer == 0:
                    msg += '0G  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, oxx, oyy, txy)
                else:
                    msg += ' H  %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, oxx, oyy, txy)
        return msg

    def writeF06_Tri3Transient(self, dt, eid, isMagPhase):
        msg = ''
        k = self.oxx[dt][eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.oxx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[dt][eid][nid][iLayer]
                oyy = self.oyy[dt][eid][nid][iLayer]
                txy = self.txy[dt][eid][nid][iLayer]
                ([fd, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], isAllZeros) = self.writeImagFloats13E([fd, oxx, oyy, txy], isMagPhase)

                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s / %13s     %13s / %13s     %13s / %-s\n' % (eid, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                else:
                    msg += '   %6s   %13s     %13s / %13s     %13s / %13s     %13s / %-s\n' % ('', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
        return msg

    def __repr__(self):
        #print "sCodes = ",self.sCodes
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        #print self.oxx.keys()
        for eid, oxxNodes in sorted(self.oxx.iteritems()):
            eType = self.eType[eid]
            k = oxxNodes.keys()
            k.remove('C')
            k.sort()
            for nid in ['C'] + k:
                for iLayer in xrange(len(self.oxx[eid][nid])):
                    fd = self.fiberCurvature[eid][nid][iLayer]
                    oxx = self.oxx[eid][nid][iLayer]
                    oyy = self.oyy[eid][nid][iLayer]
                    txy = self.txy[eid][nid][iLayer]

                    msg += '%-6i %6s %8s %7s %10g ' % (
                        eid, eType, nid, iLayer + 1, fd)
                    vals = [oxx, oyy, txy]
                    for val in vals:
                        if abs(val) < 1e-6:
                            msg += '%10s ' % ('0')
                        else:
                            try:
                                msg += '%10i %10i' % (val.real, val.imag)
                            except:
                                print("bad val = %s" % (val))
                                raise
                    msg += '\n'
        return msg


class ComplexPlateStrainObject(strainObject):
    """
    @code
    # ??? - is this just 11
    ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
      ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

    # sCode=11
                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)
      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES

    # sCode=15
                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
    ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

    # sCode=10
                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)          MAX
      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        strainObject.__init__(self, dataCode, iSubcase)
        self.eType = {}

        self.code = [self.formatCode, self.sortCode, self.sCode]

        #print self.dataCode
        self.fiberCurvature = {}
        self.exx = {}
        self.eyy = {}
        self.exy = {}

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
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ex1, ey1, exy1,
                     f2, ex2, ey2, exy2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'C': [f1, f2]}
                    self.exx[eid] = {'C': [ex1, ex2]}
                    self.eyy[eid] = {'C': [ey1, ey2]}
                    self.exy[eid] = {'C': [exy1, exy2]}
                elif eType == 'CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    #print line
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, nid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        if nid == 'CEN/4':
                            nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                    elif len(line) == 17:  # Bilinear node
                        #print line
                        (nid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.exx[eid][nid] = [ex1, ex2]
                        self.eyy[eid][nid] = [ey1, ey2]
                        self.txy[eid][nid] = [exy1, exy2]
                    else:
                        assert len(line) == 19, 'len(line)=%s' % (len(line))
                        raise NotImplementedError()
                else:
                    raise NotImplementedError('line=%s not supported...' %
                                              (line))
            return
        raise NotImplementedError('transient results not supported')

    def deleteTransient(self, dt):
        #del self.fiberCurvature[dt]
        del self.exx[dt]
        del self.eyy[dt]
        del self.exy[dt]

    def getTransients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberCurvature[dt] = {}
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.exy[dt] = {}

    def addNewEid(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)

        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        #if eid in self.exx:
            #return self.add(eid,nodeID,curvature,exx,eyy,exy)
        #assert eid not in self.exx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[eid] = {nodeID: [exx]}
        self.eyy[eid] = {nodeID: [eyy]}
        self.exy[eid] = {nodeID: [exy]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def addNewEidSort1(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)
        #print msg

        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        if dt not in self.exx:
            self.addNewTransient(dt)
        #if eid in self.exx[dt]:  # SOL200, erase the old result
            #nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid]))
            #self.deleteTransient(dt)
            #self.addNewTransient()

        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[dt][eid] = {nodeID: [exx]}
        self.eyy[dt][eid] = {nodeID: [eyy]}
        self.exy[dt][eid] = {nodeID: [exy]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add(self, dt, eid, nodeID, curvature, exx, eyy, exy):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[eid][nodeID].append(exx)
        self.eyy[eid][nodeID].append(eyy)
        self.exy[eid][nodeID].append(exy)
        if nodeID == 0:
            raise ValueError(msg)

    def addSort1(self, dt, eid, nodeID, curvature, exx, eyy, exy):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)

        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[dt][eid][nodeID].append(exx)
        self.eyy[dt][eid][nodeID].append(eyy)
        self.exy[dt][eid][nodeID].append(exy)
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNode(self, dt, eid, nodeID, curvature, exx, eyy, exy):
        #print "***addNewNode"
        #print self.oxx
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)
        assert nodeID not in self.exx[eid], msg
        self.fiberCurvature[eid][nodeID] = [curvature]
        self.exx[eid][nodeID] = [exx]
        self.eyy[eid][nodeID] = [eyy]
        self.exy[eid][nodeID] = [exy]
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']

        headers += ['exx', 'eyy', 'exy']
        if self.isVonMises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        assert f is not None
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        triMsg = None
        tri6Msg = None
        trirMsg = None

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear = True
            quadMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            qkey = eTypes.index('CQUAD8')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear = True
            quad8Msg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quad8Msg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadrMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                #print "eids = ",eids
                #print "eType = ",eType
                msgPack = msgPacks[eType]
                eids.sort()
                msg += header + msgPack
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for eid in eids:
                            out = self.writeF06_Quad4_Bilinear(eid, 4)
                    else:
                        for eid in eids:
                            out = self.writeF06_Tri3(eid)
                    ###
                elif eType in ['CTRIA3']:
                    for eid in eids:
                        out = self.writeF06_Tri3(eid)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 5)
                elif eType in ['CTRIA6', 'CTRIAR']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 3)
                else:
                    raise NotImplementedError('eType = |%r|' %
                                              (eType))  # CQUAD8, CTRIA6

                msg.append(out)
                msg.append(pageStamp + str(pageNum) + '\n')
                if f is not None:
                    f.write(''.join(msg))
                    msg = ['']
                pageNum += 1
            ###
        ###
        return (''.join(msg), pageNum - 1)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        triMsg = None
        tri6Msg = None
        trirMsg = None

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            ElemKey = eTypes.index('CQUAD4')
            #print qkey
            eid = self.eType.keys()[ElemKey]
            #print "self.oxx = ",self.oxx
            #print "eid=%s" %(eid)
            dt = self.exx.keys()[0]
            #print "dt=%s" %(dt)
            nLayers = len(self.exx[dt][eid])
            #print "elementKeys = ",elementKeys
            isBilinear = True
            quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if nLayers == 1:
                isBilinear = False
                quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msg = []
        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        dts = self.exx.keys()
        dts.sort()
        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                eids.sort()
                eType = self.eType[eid]
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.dataCode['name'], dt)
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt, eid, 4, isMagPhase)
                                msg.append(out)
                            msg.append(pageStamp + str(pageNum) + '\n')
                            if f is not None:
                                f.write(''.join(msg))
                                msg = ['']
                            pageNum += 1
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.dataCode[
                                'name'], dt)
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid, isMagPhase)
                                msg.append(out)
                            msg.append(pageStamp + str(pageNum) + '\n')
                            if f is not None:
                                f.write(''.join(msg))
                                msg = ['']
                            pageNum += 1
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.dataCode['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Tri3Transient(dt,
                                                              eid, isMagPhase)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        if f is not None:
                            f.write(''.join(msg))
                            msg = ['']
                        pageNum += 1
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.dataCode['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 5, isMagPhase)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        if f is not None:
                            f.write(''.join(msg))
                            msg = ['']
                        pageNum += 1
                elif eType in ['CTRIA6', 'CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.dataCode['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 3, isMagPhase)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        if f is not None:
                            f.write(''.join(msg))
                            msg = ['']
                        pageNum += 1
                else:
                    raise NotImplementedError('eType = |%r|' %
                                              (eType))  # CQUAD8, CTRIA6
                ###
            ###
        ###
        return (''.join(msg), pageNum - 1)

    def writeF06_Quad4_Bilinear(self, eid, n, isMagPhase):
        msg = ''
        k = self.exx[eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.exx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[eid][nid][iLayer]
                eyy = self.eyy[eid][nid][iLayer]
                exy = self.exy[eid][nid][iLayer]
                ([fd, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], isAllZeros) = self.writeImagFloats13E([fd, exx, eyy, exy], isMagPhase)

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % (eid, 'CEN/' + str(n), fd, exxr, exxi, eyyr, eyyi, exyr, exyi.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % ('', nid, fd, exxr, exxi, eyyr, eyyi, exyr, exyi.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n\n' % ('', '', fd, exxr, exxi, eyyr, eyyi, exyr, exyi.rstrip())
                else:
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Quad4_BilinearTransient(self, dt, eid, n, isMagPhase):
        msg = ''
        k = self.exx[dt][eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[dt][eid][nid][iLayer]
                eyy = self.eyy[dt][eid][nid][iLayer]
                exy = self.exy[dt][eid][nid][iLayer]
                ([fd, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], isAllZeros) = self.writeImagFloats13E([fd, exx, eyy, exy], isMagPhase)

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % (eid, 'CEN/' + str(n), fd, exxr, exxi, eyyr, eyyi, exyr, exyi.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n' % ('', nid, fd, exxr, exxi, eyyr, eyyi, exyr, exyi.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s   %13s / %13s   %13s / %13s   %13s /   %-s\n\n' % ('', '', fd, exxr, exxi, eyyr, eyyi, exyr, exyi.rstrip())
                else:
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Tri3(self, eid, isMagPhase):
        msg = ''
        k = self.exx[eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.exx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[eid][nid][iLayer]
                eyy = self.eyy[eid][nid][iLayer]
                exy = self.exy[eid][nid][iLayer]

                ([fd, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], isAllZeros) = self.writeImagFloats13E([fd, exx, eyy, exy], isMagPhase)
                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s / %13s     %13s / %13s     %13s / %-s\n' % (eid, fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
                else:
                    msg += '   %6s   %13s     %13s / %13s     %13s / %13s     %13s / %-s\n' % ('', fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
        return msg

    def writeF06_Tri3Transient(self, dt, eid, isMagPhase):
        msg = ''
        k = self.exx[dt][eid].keys()
        k.remove('C')
        k.sort()
        for nid in ['C'] + k:
            for iLayer in xrange(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[dt][eid][nid][iLayer]
                eyy = self.eyy[dt][eid][nid][iLayer]
                exy = self.exy[dt][eid][nid][iLayer]

                ([fd, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], isAllZeros) = self.writeImagFloats13E([fd, exx, eyy, exy], isMagPhase)
                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s / %13s     %13s / %13s     %13s / %-s\n' % (eid, fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
                else:
                    msg += '   %6s   %13s     %13s / %13s     %13s / %13s     %13s / %-s\n' % ('', fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
        return msg

    def __repr__(self):
        #print "dt = ",self.dt
        #print "nf = ",self.nonlinearFactor
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for eid, exxNodes in sorted(self.exx.iteritems()):
            eType = self.eType[eid]
            k = exxNodes.keys()
            k.remove('C')
            k.sort()
            for nid in ['C'] + k:
                for iLayer in xrange(len(self.exx[eid][nid])):
                    fd = self.fiberCurvature[eid][nid][iLayer]
                    exx = self.exx[eid][nid][iLayer]
                    eyy = self.eyy[eid][nid][iLayer]
                    exy = self.exy[eid][nid][iLayer]

                    msg += '%-6i %6s %8s %7s %10g ' % (
                        eid, eType, nid, iLayer + 1, fd)
                    vals = [exx, eyy, exy]
                    for val in vals:
                        if abs(val) < 1e-6:
                            msg += '%10s ' % ('0.')
                        else:
                            msg += '%10.3g ' % (val)
                    msg += '\n'

                    #msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy)
        return msg

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, exx in sorted(self.exx.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid, exxNodes in sorted(exx.iteritems()):
                eType = self.eType[eid]
                k = exxNodes.keys()
                k.remove('C')
                k.sort()
                for nid in ['C'] + k:
                    for iLayer in xrange(len(self.exx[dt][eid][nid])):
                        fd = self.fiberCurvature[eid][nid][iLayer]
                        exx = self.exx[dt][eid][nid][iLayer]
                        eyy = self.eyy[dt][eid][nid][iLayer]
                        exy = self.exy[dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s %10g ' % (eid, eType,
                                                           nid, iLayer + 1, fd)
                        vals = [exx, eyy, exy]
                        for val in vals:
                            if abs(val) < 1e-6:
                                msg += '%10s ' % ('0.')
                            else:
                                msg += '%10.3g ' % (val)
                            ###
                        msg += '\n'

                        #msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy)
        return msg
