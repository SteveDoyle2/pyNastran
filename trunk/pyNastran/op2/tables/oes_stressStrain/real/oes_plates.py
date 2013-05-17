from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import zeros, array

from pyNastran import as_array
from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeFloats8p4F


class PlateStressObject(StressObject):
    """
    ::
    
      ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)
        ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
            6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                         1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03
  
  
                           S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)          MAX
        ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
            6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                         1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        #self.append_data_member('sCodes','s_code')
        #print "self.s_code = ",self.s_code
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.fiberCurvature = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.ovmShear = {}
        
        self.ovmShear2 = None

        # eid_layer = [eid1, ilayer1]
        #             [eid1, ilayer2]
        #             [eid2, ilayer1]
        #             [eid2, ilayer2]
        
        # ovm[1] = [1234.5]
        #    [2] = [1234.5]
        #    [3]   [5678.8]
        #    [4]   [5678.8]
        #
        # >>> elements2, col = numpy.where(eid_layer[:,0] == eid2)
        # >>> print elements
        # array([2, 3])
        # >>> ovm[elements2]
        # [5678.8, 5678.8]

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.addNewNodeSort1
                
                self.add_array = self.add_array_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2
            self.addNewNode = self.addNewNodeSort2
            self.add_array = self.add_array_sort2

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
        msg.append('  eType, fiberCurvature, oxx, oyy, txy, angle, '
                   'majorP, minorP, ovmShear\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            #print(data)
            
            eType = data[0][0]
            print('eType = %s' % eType)
            
            n = 0
            line2 = data[0]
            if eType == 'CTRIA3':
                n = 1
            elif eType == 'CQUAD4':
                if len(line2) == 19:  # Centroid - bilinear
                    n += 5
                elif len(line2) == 18:  # Centroid
                    n += 1
            
            n *= 2 # 2 layers - top & bottom
            #n = 10
            eTypes = zeros(n, dtype='string')
            ovmShear = zeros(n, dtype='float32')
            index_to_elementNodeLayer = zeros((n, 3), dtype='int32')  # layer=0 -> center
            
            #print(index_to_elementNodeLayer)
            #import sys
            #print(data)
            #sys.exit()
            i = 0
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                     f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line

                    if as_array:
                        index_to_elementNodeLayer[i]   = [eid, 0, 0]
                        index_to_elementNodeLayer[i+1] = [eid, 0, 1]
                        #eTypes[i:i+1]   = [eType, eType]
                        #ovmShear[i:i+1] = [ovm1, ovm2]

                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'C': [f1, f2]}
                    self.oxx[eid] = {'C': [ox1, ox2]}
                    self.oyy[eid] = {'C': [oy1, oy2]}
                    self.txy[eid] = {'C': [txy1, txy2]}
                    self.angle[eid] = {'C': [angle1, angle2]}
                    self.majorP[eid] = {'C': [o11, o12]}
                    self.minorP[eid] = {'C': [o21, o22]}
                    self.ovmShear[eid] = {'C': [ovm1, ovm2]}
                    i += 2
                elif eType == 'CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, nid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                          f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        if as_array:
                            index_to_elementNodeLayer[i]   = [eid, 0, 0]
                            index_to_elementNodeLayer[i+1] = [eid, 0, 1]
                            #print(eTypes)
                            #a = array([eType, eType])
                            #print(a.shape)
                            #print(eTypes.shape)
                            #print(i)
                            #eTypes[i:i+1] = a
                            #print(eTypes)
                            #ovmShear[i:i+1] = [ovm1, ovm2]

                        if nid == 'CEN/4':
                            nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [o11, o12]}
                        self.minorP[eid] = {nid: [o21, o22]}
                        self.ovmShear[eid] = {nid: [ovm1, ovm2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                     f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [o11, o12]}
                        self.minorP[eid] = {nid: [o21, o22]}
                        self.ovmShear[eid] = {nid: [ovm1, ovm2]}
                    elif len(line) == 17:  # Bilinear
                        #print line
                        (nid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                              f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.oxx[eid][nid] = [ox1, ox2]
                        self.oyy[eid][nid] = [oy1, oy2]
                        self.txy[eid][nid] = [txy1, txy2]
                        self.angle[eid][nid] = [angle1, angle2]
                        self.majorP[eid][nid] = [o11, o12]
                        self.minorP[eid][nid] = [o21, o22]
                        self.ovmShear[eid][nid] = [ovm1, ovm2]
                    else:
                        assert len(line) == 19, 'line=%s len(line)=%s' % (line, len(line))
                        raise NotImplementedError()
                else:
                    msg = 'line=%s not supported...' % (line)
                    raise NotImplementedError(msg)

            self.eTypes2 = eTypes
            self.ovmShear2 = ovmShear
            return
        #for line in data:
            #print line
        raise NotImplementedError('transient results not supported')

    def delete_transient(self, dt):
        #del self.fiberCurvature[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.ovmShear[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        #self.fiberCurvature[dt] = {}
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.txy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.ovmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print "***add_new_eid..."
        #if eid in self.oxx:
            #return self.add(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)

        #assert eid not in self.oxx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [fd]}
        assert isinstance(eid, int)
        self.oxx[eid] = {nodeID: [oxx]}
        self.oyy[eid] = {nodeID: [oyy]}
        self.txy[eid] = {nodeID: [txy]}
        self.angle[eid] = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.ovmShear[eid] = {nodeID: [ovm]}
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        if nodeID == 0:
            raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print "***addNewEidTransient..."
        #msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g major=%g vm=%g" % (
            dt, eid, nodeID, fd, oxx, majorP, ovm)
        #print msg
        #if eid in self.ovmShear[dt]:
        #    return self.add(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)

        if dt in self.oxx and eid in self.oxx[dt]:  # SOL200, erase the old result
            #nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s major=%s vm=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid]),str(self.majorP[dt][eid][nid]),str(self.ovmShear[dt][eid][nid]))
            self.delete_transient(dt)
            self.add_new_transient(dt)

        assert isinstance(eid, int)
        self.eType[eid] = eType
        if dt not in self.oxx:
            self.add_new_transient(dt)
        self.fiberCurvature[eid] = {nodeID: [fd]}
        self.oxx[dt][eid] = {nodeID: [oxx]}
        self.oyy[dt][eid] = {nodeID: [oyy]}
        self.txy[dt][eid] = {nodeID: [txy]}
        self.angle[dt][eid] = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.ovmShear[dt][eid] = {nodeID: [ovm]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print "***add"
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        assert isinstance(eid, int)
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[eid][nodeID].append(oxx)
        self.oyy[eid][nodeID].append(oyy)
        self.txy[eid][nodeID].append(txy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.ovmShear[eid][nodeID].append(ovm)
        if nodeID == 0:
            raise ValueError(msg)

    def add_array_sort1(self, dt, sadf):
        raise NotImplementedError()

    def add_array_sort2(self, dt, sadf):
        raise NotImplementedError()

    def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print "***addTransient"
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" % (dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        #print self.oxx
        #print self.fiberCurvatrure
        assert eid is not None
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[dt][eid][nodeID].append(oxx)
        self.oyy[dt][eid][nodeID].append(oyy)
        self.txy[dt][eid][nodeID].append(txy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.ovmShear[dt][eid][nodeID].append(ovm)
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNode(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print "***addNewNode"
        assert eid is not None
        #assert nodeID not in self.oxx[eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[eid][nodeID] = [oxx]
        self.oyy[eid][nodeID] = [oyy]
        self.txy[eid][nodeID] = [txy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.ovmShear[eid][nodeID] = [ovm]
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNodeSort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print "***addNewNodeTransient"
        #print self.oxx
        assert eid is not None
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        #assert nodeID not in self.oxx[dt][eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[dt][eid][nodeID] = [oxx]
        self.oyy[dt][eid][nodeID] = [oyy]
        self.txy[dt][eid][nodeID] = [txy]
        self.angle[dt][eid][nodeID] = [angle]
        self.majorP[dt][eid][nodeID] = [majorP]
        self.minorP[dt][eid][nodeID] = [minorP]
        self.ovmShear[dt][eid][nodeID] = [ovm]
        if nodeID == 0:
            raise ValueError(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
        headers += ['oxx', 'oyy', 'txy', 'majorP', 'minorP']
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
            msg += '%10s ' % header
        msg += '\n'

        #print self.oxx.keys()
        for dt, oxxs in sorted(self.oxx.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid, oxxNodes in sorted(oxxs.iteritems()):
                eType = self.eType[eid]
                for nid in sorted(oxxNodes):
                    for iLayer in xrange(len(self.oxx[dt][eid][nid])):
                        fd = self.fiberCurvature[eid][nid][iLayer]
                        oxx = self.oxx[dt][eid][nid][iLayer]
                        oyy = self.oyy[dt][eid][nid][iLayer]
                        txy = self.txy[dt][eid][nid][iLayer]
                       #angle = self.angle[ dt][eid][nid][iLayer]
                        major = self.majorP[dt][eid][nid][iLayer]
                        minor = self.minorP[dt][eid][nid][iLayer]
                        ovm = self.ovmShear[dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s %10g ' % (eid, eType,
                                                           nid, iLayer + 1, fd)
                        vals = [oxx, oyy, txy, major, minor, ovm]
                        for val in vals:
                            if abs(val) < 1e-6:
                                msg += '%10s ' % '0'
                            else:
                                try:
                                    msg += '%10i ' % int(val)
                                except:
                                    print("bad val = %s" % val)
                                    raise
                        msg += '\n'
        return msg

    def write_matlab(self, name, isubcase, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

        #if self.isVonMises():
        #    vonMises = 'vonMises'
        #else:
        #    vonMises = 'maxShear'

        #if self.isFiberDistance():
        #    fiberCurvature = 'fiberDistance'
        #else:
        #    fiberCurvature = 'fiberCurvature'

        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

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
                msgPack = msgPacks[eType]

                msg += header + msgPack
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for eid in eids:
                            out = self._write_matlab_quad4_bilinear(eid, 4)
                            msg.append(out)
                    else:
                        for eid in eids:
                            out = self._write_matlab_tri3(eid)
                            msg.append(out)

                elif eType in ['CTRIA3']:
                    a = 'fem.plateStress(%i).tri3.elementIDs = %s\n' % (
                        isubcase, eids)
                    b = 'fem.plateStress(%i).tri3.oxx = [' % isubcase
                    for eid in eids:
                        out = self._write_matlab_tri3(eid)
                        msg.append(out)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self._write_matlab_quad4_bilinear(eid, 5)
                        msg.append(out)
                elif eType in ['CTRIAR', 'CTRIA6']:
                    for eid in eids:
                        out = self._write_matlab_quad4_bilinear(eid, 3)
                        msg.append(out)
                else:
                    raise NotImplementedError('eType = |%r|' % eType)  # CQUAD8, CTRIA6
                f.write(''.join(msg))
                msg = []

    def _write_matlab_tri3(self, eid):
        msg = ''
        oxxNodes = self.oxx[eid].keys()
        for nid in sorted(oxxNodes):
            for iLayer in xrange(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[eid][nid][iLayer]
                oyy = self.oyy[eid][nid][iLayer]
                txy = self.txy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
        return msg

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

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
                elif eType in ['CTRIA3']:
                    for eid in eids:
                        out = self.writeF06_Tri3(eid)
                        msg.append(out)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 5)
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
        return (''.join(msg), pageNum - 1)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None

        eTypes = self.eType.values()
        dts = self.oxx.keys()
        dt = dts[0]
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            #print "qkey=%s kkey=%s" %(qkey,kkey)
            ekey = self.oxx[dt][kkey].keys()
            isBilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

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
                            header[1] = ' %s = %10.4E\n' % (self.data_code[
                                'name'], dt)
                            msg += header + msgPack
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt,
                                                                            eid, 4)
                                msg.append(out)
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (
                                self.data_code['name'], dt)
                            msg += header + msgPack
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid)
                                msg.append(out)
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self.writeF06_Tri3Transient(dt, eid)
                            msg.append(out)
                elif eType in ['CTRIA6', 'CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(
                                dt, eid, 3)
                            msg.append(out)
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'],
                                                        dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(
                                dt, eid, 5)
                            msg.append(out)
                else:
                    raise NotImplementedError('eType = |%r|' % eType)  # CQUAD8, CTRIA6

                msg.append(pageStamp + str(pageNum) + '\n')
                if f is not None:
                    f.write(''.join(msg))
                    msg = ['']
                pageNum += 1

        return (''.join(msg), pageNum - 1)

    def writeF06_Quad4_Bilinear(self, eid, n):
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
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % (eid, 'CEN/' + str(n), fd, oxx, oyy, txy, angle, major, minor, ovm)
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % ('', nid, fd, oxx, oyy, txy, angle, major, minor, ovm)
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' % ('', '', fd, oxx, oyy, txy, angle, major, minor, ovm)
                else:
                    raise Exception('Invalid option for cquad4')

        return msg

    def writeF06_Quad4_BilinearTransient(self, dt, eid, n):
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
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                ovm = self.ovmShear[dt][eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % (eid, 'CEN/' + str(n), fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % ('', nid, fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' % ('', '', fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                else:
                    #msg += '   %8s %8s  %13E  %13E %13E %13E   %8.4F  %13E %13E %13E\n' %('','',  fd,oxx,oyy,txy,angle,major,minor,ovm)
                    raise RuntimeError('Invalid option for cquad4')
        return msg

    def writeF06_Tri3(self, eid):
        msg = ''
        oxxNodes = self.oxx[eid].keys()
        for nid in sorted(oxxNodes):
            for iLayer in xrange(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[eid][nid][iLayer]
                oyy = self.oyy[eid][nid][iLayer]
                txy = self.txy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
        return msg

    def writeF06_Tri3Transient(self, dt, eid):
        msg = ''
        oxxNodes = self.oxx[dt][eid].keys()
        for nid in sorted(oxxNodes):
            for iLayer in xrange(len(self.oxx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[dt][eid][nid][iLayer]
                oyy = self.oyy[dt][eid][nid][iLayer]
                txy = self.txy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                ovm = self.ovmShear[dt][eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, oxx, oyy, txy, angle, major, minor, ovm)
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, oxx, oyy, txy, angle, major, minor, ovm)
        return msg

    def __repr__(self):
        #print "sCodes = ",self.sCodes
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        #print self.oxx.keys()
        for eid, oxxNodes in sorted(self.oxx.iteritems()):
            eType = self.eType[eid]
            for nid in sorted(oxxNodes):
                for iLayer in xrange(len(self.oxx[eid][nid])):
                    fd = self.fiberCurvature[eid][nid][iLayer]
                    oxx = self.oxx[eid][nid][iLayer]
                    oyy = self.oyy[eid][nid][iLayer]
                    txy = self.txy[eid][nid][iLayer]
                   #angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    ovm = self.ovmShear[eid][nid][iLayer]

                    msg += '%-6i %6s %8s %7s %10g ' % (
                        eid, eType, nid, iLayer + 1, fd)
                    vals = [oxx, oyy, txy, major, minor, ovm]
                    for val in vals:
                        if abs(val) < 1e-6:
                            msg += '%10s ' % '0'
                        else:
                            try:
                                msg += '%10i ' % val
                            except:
                                print("bad val = %s" % val)
                                raise
                    msg += '\n'
        return msg


class PlateStrainObject(StrainObject):
    """
    ::

      # ??? - is this just 11
      ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
        ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
  
      # s_code=11
                             S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)
        ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES
  
      # s_code=15
                             S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
      ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
        ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
  
      # s_code=10
                             S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)          MAX
        ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]

        #print self.data_code
        self.fiberCurvature = {}
        self.exx = {}
        self.eyy = {}
        self.exy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.evmShear = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.exx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberCurvature, exx, eyy, exy, angle, '
                   'majorP, minorP, evmShear\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                     f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'C': [f1, f2]}
                    self.exx[eid] = {'C': [ex1, ex2]}
                    self.eyy[eid] = {'C': [ey1, ey2]}
                    self.exy[eid] = {'C': [exy1, exy2]}
                    self.angle[eid] = {'C': [angle1, angle2]}
                    self.majorP[eid] = {'C': [e11, e12]}
                    self.minorP[eid] = {'C': [e21, e22]}
                    self.evmShear[eid] = {'C': [evm1, evm2]}
                elif eType == 'CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    #print line
                    if len(line) == 19:  # Centroid - bilinear
                        (
                            eType, eid, nid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                            f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        if nid == 'CEN/4':
                            nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [e11, e12]}
                        self.minorP[eid] = {nid: [e21, e22]}
                        self.evmShear[eid] = {nid: [evm1, evm2]}
                    elif len(line) == 18:  # Centroid
                        (
                            eType, eid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                            f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        nid = 'C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [e11, e12]}
                        self.minorP[eid] = {nid: [e21, e22]}
                        self.evmShear[eid] = {nid: [evm1, evm2]}
                    elif len(line) == 17:  # Bilinear node
                        #print line
                        (nid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                         f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.exx[eid][nid] = [ex1, ex2]
                        self.eyy[eid][nid] = [ey1, ey2]
                        self.txy[eid][nid] = [exy1, exy2]
                        self.angle[eid][nid] = [angle1, angle2]
                        self.majorP[eid][nid] = [e11, e12]
                        self.minorP[eid][nid] = [e21, e22]
                        self.evmShear[eid][nid] = [evm1, evm2]
                    else:
                        assert len(line) == 19, 'len(line)=%s' % len(line)
                        raise NotImplementedError()
                else:
                    msg = 'line=%s not supported...' % (line)
                    raise NotImplementedError(msg)
            return
        raise NotImplementedError('transient results not supported')

    def delete_transient(self, dt):
        #del self.fiberCurvature[dt]
        del self.exx[dt]
        del self.eyy[dt]
        del self.exy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.evmShear[dt]

    def get_transients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberCurvature[dt] = {}
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.exy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.evmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)

        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        #if eid in self.exx:
            #return self.add(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #assert eid not in self.exx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[eid] = {nodeID: [exx]}
        self.eyy[eid] = {nodeID: [eyy]}
        self.exy[eid] = {nodeID: [exy]}
        self.angle[eid] = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.evmShear[eid] = {nodeID: [evm]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #print msg

        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        if dt not in self.exx:
            self.add_new_transient(dt)
        #if eid in self.evmShear[dt]:  # SOL200, erase the old result
            #nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s major=%s vm=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid]),str(self.majorP[dt][eid][nid]),str(self.ovmShear[dt][eid][nid]))
            #self.delete_transient(dt)
            #self.add_new_transient()

        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[dt][eid] = {nodeID: [exx]}
        self.eyy[dt][eid] = {nodeID: [eyy]}
        self.exy[dt][eid] = {nodeID: [exy]}
        self.angle[dt][eid] = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.evmShear[dt][eid] = {nodeID: [evm]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[eid][nodeID].append(exx)
        self.eyy[eid][nodeID].append(eyy)
        self.exy[eid][nodeID].append(exy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.evmShear[eid][nodeID].append(evm)
        if nodeID == 0:
            raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        if nodeID != 'C':  # centroid
            assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)

        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[dt][eid][nodeID].append(exx)
        self.eyy[dt][eid][nodeID].append(eyy)
        self.exy[dt][eid][nodeID].append(exy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.evmShear[dt][eid][nodeID].append(evm)
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNode(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "***addNewNode"
        #print self.oxx
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        assert nodeID not in self.exx[eid], msg
        self.fiberCurvature[eid][nodeID] = [curvature]
        self.exx[eid][nodeID] = [exx]
        self.eyy[eid][nodeID] = [eyy]
        self.exy[eid][nodeID] = [exy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.evmShear[eid][nodeID] = [evm]
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']

        headers += ['exx', 'eyy', 'exy', 'eMajor', 'eMinor']
        if self.isVonMises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)

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
                                              eType)  # CQUAD8, CTRIA6
                msg.append(out)
                msg.append(pageStamp + str(pageNum) + '\n')
                if f is not None:
                    f.write(''.join(msg))
                    msg = ['']
                pageNum += 1
        return (''.join(msg), pageNum - 1)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % vonMises]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % vonMises]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % vonMises]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % vonMises]

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
                            header[1] = ' %s = %10.4E\n' % (
                                self.data_code['name'], dt)
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt,
                                                                            eid, 4)
                                msg.append(out)
                            msg.append(pageStamp + str(pageNum) + '\n')
                            pageNum += 1
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (
                                self.data_code['name'], dt)
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid)
                                msg.append(out)
                            msg.append(pageStamp + str(pageNum) + '\n')
                            pageNum += 1
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'],
                                                        dt)
                        for eid in eids:
                            out = self.writeF06_Tri3Transient(dt, eid)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        pageNum += 1
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'],
                                                        dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(
                                dt, eid, 5)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        pageNum += 1
                elif eType in ['CTRIA6', 'CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'],
                                                        dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(
                                dt, eid, 3)
                            msg.append(out)
                        msg.append(pageStamp + str(pageNum) + '\n')
                        pageNum += 1
                else:
                    raise NotImplementedError('eType = |%r|' %
                                              eType)  # CQUAD8, CTRIA6
                if f is not None: # stream write the file to save memory
                    f.write(''.join(msg))
                    msg = ['']
        return (''.join(msg), pageNum - 1)

    def writeF06_Quad4_Bilinear(self, eid, n):
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
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                evm = self.evmShear[eid][nid][iLayer]
                ([fd, exx, eyy, exy, major, minor, evm], isAllZeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % (eid, 'CEN/' + str(n), fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % ('', nid, fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' % ('', '', fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                else:
                    raise RuntimeError('Invalid option for cquad4')
        return msg

    def writeF06_Quad4_BilinearTransient(self, dt, eid, n):
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
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                evm = self.evmShear[dt][eid][nid][iLayer]

                ([fd, exx, eyy, exy, major, minor, evm], isAllZeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])

                if nid == 'C' and iLayer == 0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % (eid, 'CEN/' + str(n), fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                elif iLayer == 0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % ('', nid, fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                elif iLayer == 1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' % ('', '', fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                else:
                    raise RuntimeError('Invalid option for cquad4')
        return msg

    def writeF06_Tri3(self, eid):
        msg = ''
        k = self.exx[eid].keys()
        for nid in sorted(k):
            for iLayer in xrange(len(self.exx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[eid][nid][iLayer]
                eyy = self.eyy[eid][nid][iLayer]
                exy = self.exy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                evm = self.evmShear[eid][nid][iLayer]

                ([fd, exx, eyy, exy, major, minor, evm], isAllZeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])
                if iLayer == 0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, exx, eyy, exy, angle, major, minor, evm.rstrip())
        return msg

    def writeF06_Tri3Transient(self, dt, eid):
        msg = ''
        exxNodes = self.exx[dt][eid]
        #k = exxNodes.keys()
        for nid in sorted(exxNodes):
            for iLayer in xrange(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[dt][eid][nid][iLayer]
                eyy = self.eyy[dt][eid][nid][iLayer]
                exy = self.exy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                evm = self.evmShear[dt][eid][nid][iLayer]

                ([fd, exx, eyy, exy, major, minor, evm], isAllZeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], isAllZeros) = writeFloats8p4F([angle])
                if iLayer == 0:
                    msg += ('0  %6i   %13s     %13s  %13s  %13s   %8s   '
                            '%13s   %13s  %-s\n' % (eid, fd, exx, eyy, exy,
                                                    angle, major, minor, evm))
                else:
                    msg += ('   %6s   %13s     %13s  %13s  %13s   %8s   '
                            '%13s   %13s  %-s\n' % ('', fd, exx, eyy, exy,
                                                    angle, major, minor, evm))
        return msg

    def __repr__(self):
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for eid, exxNodes in sorted(self.exx.iteritems()):
            eType = self.eType[eid]
            for nid in sorted(exxNodes):
                for iLayer in xrange(len(self.exx[eid][nid])):
                    fd = self.fiberCurvature[eid][nid][iLayer]
                    exx = self.exx[eid][nid][iLayer]
                    eyy = self.eyy[eid][nid][iLayer]
                    exy = self.exy[eid][nid][iLayer]
                   #angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    evm = self.evmShear[eid][nid][iLayer]

                    msg += '%-6i %6s %8s %7s %10g ' % (
                        eid, eType, nid, iLayer + 1, fd)
                    vals = [exx, eyy, exy, major, minor, evm]
                    for val in vals:
                        if abs(val) < 1e-6:
                            msg += '%10s ' % '0.'
                        else:
                            msg += '%10.3g ' % val
                    msg += '\n'

                    #msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g evm=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy,evm)
        return msg

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' % ('EID', 'eType', 'nodeID', 'iLayer')
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, exx in sorted(self.exx.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid, exxNodes in sorted(exx.iteritems()):
                eType = self.eType[eid]
                for nid in sorted(exxNodes):
                    for iLayer in xrange(len(self.exx[dt][eid][nid])):
                        fd = self.fiberCurvature[eid][nid][iLayer]
                        exx = self.exx[dt][eid][nid][iLayer]
                        eyy = self.eyy[dt][eid][nid][iLayer]
                        exy = self.exy[dt][eid][nid][iLayer]
                       #angle = self.angle[ dt][eid][nid][iLayer]
                        major = self.majorP[dt][eid][nid][iLayer]
                        minor = self.minorP[dt][eid][nid][iLayer]
                        evm = self.evmShear[dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s %10g ' % (eid, eType,
                                                           nid, iLayer + 1, fd)
                        vals = [exx, eyy, exy, major, minor, evm]
                        for val in vals:
                            if abs(val) < 1e-6:
                                msg += '%10s ' % '0.'
                            else:
                                msg += '%10.3g ' % val
                        msg += '\n'

                        #msg += ('eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g '
                        #        'eyy=%-9.3g exy=%-9.3g evm=%-9.3g\n'
                        #      % (eid, eType, nid, iLayer, exx, eyy, exy, evm))
        return msg