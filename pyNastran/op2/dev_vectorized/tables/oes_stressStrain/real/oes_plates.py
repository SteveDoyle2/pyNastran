from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import zeros, array
import pandas as pd

from pyNastran import update_index

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeFloats8p4F

class RealPlateResults(object):
    def __init__(self):
        self.shape = {}
        self._inode_start = None
        self._inode_end = None
        self._ielement_start = None
        self._ielement_end = None
        self._ncount = 0
        self.data = None
        self.element_data = None

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

    def _is_full(self, nnodes, nelements):
        self._inode_start += nnodes
        self._ielement_start += nelements

        #print('check....inode_start=%s size_node_end=%s' %(self._inode_start, self._size_node_end))
        #print('check...._ielement_start=%s size_element_end=%s' %(self._ielement_start, self._size_element_end))
        if self._inode_start == self._size_node_end:
            return True
        elif self._ielement_start == self._size_element_end:
            self._ielement_start = 0
        return False

    def _increase_size(self, dt, nnodes, nelements):
        if dt in self.shape:  # default dictionary
            self.shape[dt][0] += nnodes
            self.shape[dt][1] += nelements
        else:
            self.shape[dt] = [nnodes, nelements]
        #print("shape =", self.shape)

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nnodes = self.shape[shape0][0]
        nelements = self.shape[shape0][1]
        #print("ndt=%s nnodes=%s dts=%s" % (str(ndt), str(nnodes), str(dts)))
        return ndt, nnodes, nelements, dts

    #def _increment(self, nnodes, nelements):
        #self._inode_start += nnodes
        #self._inode_end += nnodes
        #self._ielement_start += nelements
        #self._ielement_end += nelements
        #return self._inode_start, self._inode_end, self._ielement_start, self._ielement_end

    def _preallocate(self, dt, nnodes, nelements):
        """
        nodes will create self.data
        elements will create self.element_data

        data is primary, so nodes are input first
        """
        #print('---preallocate nnodes=%s nelements=%s' % (nnodes, nelements))
        if self.shape is None:
            self._inode_start += nnodes
            self._inode_end += nnodes
            self._ielement_start += nelements
            self._ielement_end += nelements
        #elif self.data is not None:
            #ndt, nnodes_size, nelements_size, dts = self._get_shape()
            #inode_end = min(self._inode_end, nnodes_size)
            #ielement_end = min(self._ielement_end, nelements_size)
            #return (self._inode_start, inode_end, self._ielement_start, ielement_end)
            #print "data =", self.data
            #pass
        else:
            ndt, nnodes_size, nelements_size, dts = self._get_shape()
            #print("ndt=%s nelements_size=%s nnodes_size=%s dts=%s" % (ndt, nelements_size, nnodes_size, str(dts)))

            if self._inode_start is not None:
                return (self._inode_start, self._inode_start + nnodes,
                        self._ielement_start, self._ielement_start + nelements)

            n = ndt * nnodes_size
            if self._ncount != 0:
                asfd
            self._ncount += 1
            self._inode_start = 0
            self._inode_end = nnodes

            self._ielement_start = 0
            self._ielement_end = nelements

            data = {}
            element_data = {}
            columns = []
            if dts[0] is not None:
                name = self.data_code['name']
                #print('***name=%r***' % name)
                if isinstance(dt, int):
                    data[name] = pd.Series(zeros((n), dtype='int32'))
                else:
                    data[name] = pd.Series(zeros((n), dtype='float32'))
                columns.append(name)

            element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))
            element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))
            element_data['nlayers'] = pd.Series(zeros(nelements_size, dtype='int32'))

            data['element_id'] = pd.Series(zeros((n), dtype='int32'))
            data['node_id'] = pd.Series(zeros((n), dtype='int32'))
            data['layer'] = pd.Series(zeros((n), dtype='int32'))

            #columns.append('element_type')

            #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
            #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
            #print('n =', n)


            headers = self._get_headers()
            (fd, oxx, oyy, txy, omax, omin, ovm) = headers

            # fiber_distance / fiber_curvature
            data[fd] = pd.Series(zeros((n), dtype='float32'))

            data[oxx] = pd.Series(zeros((n), dtype='float32'))
            data[oyy] = pd.Series(zeros((n), dtype='float32'))
            data[txy] = pd.Series(zeros((n), dtype='float32'))

            data['angle'] = pd.Series(zeros((n), dtype='float32'))
            data[omax] = pd.Series(zeros((n), dtype='float32'))
            data[omin] = pd.Series(zeros((n), dtype='float32'))
            data[ovm] = pd.Series(zeros((n), dtype='float32'))

            columns += ['element_id', 'node_id', 'layer', fd, oxx, oyy, txy, 'angle', omax, omin, ovm]

            self.data = pd.DataFrame(data, columns=columns)
            self.element_data = pd.DataFrame(element_data, columns=['element_id', 'element_type', 'nnodes'])
            size_end = n

            # max sizes
            self._size_node_start = 0
            self._size_node_end = n
            self._size_element_start = 0
            self._size_element_end = nelements_size

            self._inode_start = 0
            self._inode_end = nnodes
            self._ielement_start = 0
            self._ielement_end = nelements
            #print('_inode=%s _ielement=%s' %(nnodes, nelements))
        return (self._inode_start, self._inode_end, self._ielement_start, self._ielement_end)

    def _finalize(self, dt):
        ndt, nnodes, nelements, dts = self._get_shape()

        if dt is not None and dt != dts[-1]:
            return
        #print("----finalize----")

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        update_index = True
        #print("final A\n", self.data.to_string())
        if update_index:
            if dts[0] is not None:
                name = self.data_code['name']
                self.data = self.data.set_index([name, 'element_id', 'node_id', 'layer'])
            else:
                self.data = self.data.set_index(['element_id', 'node_id', 'layer'])
            self.element_data = self.element_data.set_index(['element_id'])

        #print(self.element_data.to_string())
        #print("final\n", self.data.to_string())
        #print(self.element_data.to_string())
        del self._inode_start
        del self._inode_end
        #print('---PlateStressObject---')
        #print(self.data.to_string())

    def get_element_types(self):
        etypes = self.element_data['element_type']
        return list(set(etypes))

    def get_stats(self):
        ndt, nnodes, nelements, dts = self._get_shape()
        msg = self._get_data_code()
        if self.nonlinear_factor is not None:  # transient
            name = self.data_code['name']
            dt_string = '%s, ' % name
            msg.append('  real type=%s n%ss=%s nelements=%s\n'
                       % (self.__class__.__name__, name, ndt, nelements))
        else:
            dt_string = ''
            msg.append('  real type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        headers = self._get_headers()
        (fd, oxx, oyy, txy, omax, omin, ovm) = headers

        etypes = self.get_element_types()
        msg.append('  element data: index :  element_id\n')
        msg.append('              : result:  element_type, nnodes\n')
        msg.append('  data        : index :  %selement_id, node_id, layer\n' % dt_string)
        msg.append('              : result:  %s, %s, %s, %s, '
                                            '%s, %s, %s, angle\n' % (fd, oxx, oyy, txy,
                                                              omax, omin, ovm) )
        msg.append('                element_types: %s\n' %(', '.join(etypes)))
        return msg

    def _get_first_type_index(self):
        first_type_index = {}
        #print(self.element_data.to_string())
        for i in xrange(len(self.element_data)):
            eid = self.element_data.index[i]
            etype = self.element_data['element_type'][eid]
            if etype not in first_type_index:
                first_type_index[etype] = eid
            #print("i=%s eid=%s" % (i, eid))
        return first_type_index

    def _write_f06_helper(self, msgPacks, pageStamp, f, pageNum):
        msg = []
        isBilinear = True

        (kfd, koxx, koyy, ktxy, komax, komin, kovm) = self._get_headers()
        i = 0
        j = 0
        (eid, nid, ilayer) = self.data.index[i]
        etype_old = self.element_data.ix[eid]['element_type']

        ndata = len(self.data)
        try:
            while i < ndata:
                (eid, nid, ilayer) = self.data.index[i]
                #print('eid=%s nid=%s ilayer=%s' % (eid, nid, ilayer))
                etype = self.element_data.ix[eid]['element_type']
                msg = msgPacks[etype]
                if etype in ['CQUAD4']:
                    if isBilinear:
                        n = 4
                        cen4 = 'CEN/' + str(n)
                        eid2 = eid
                        while eid2 == eid and i < ndata:
                            (eid, nid, ilayer) = self.data.index[i]
                            index = self.data.index[i]
                            data = self.data.ix[index]

                            fd = data[kfd]
                            angle = data['angle']
                            oxx = data[koxx]
                            oyy = data[koyy]
                            txy = data[ktxy]
                            major = data[komax]
                            minor = data[komin]
                            ovm = data[kovm]

                            ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy,
                                                                                                   major, minor, ovm])
                            ([angle], isAllZeros) = writeFloats8p4F([angle])

                            # ..note:: nid==0 is the center of the CQUAD4
                            if nid == 0 and ilayer == 1:
                                msg.append('0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % (eid, cen4, fd, oxx, oyy, txy, angle, major, minor, ovm))
                            elif ilayer == 1:
                                msg.append('   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % ('', nid, fd, oxx, oyy, txy, angle, major, minor, ovm))
                            elif ilayer == 2:
                                msg.append('   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' % ('', '', fd, oxx, oyy, txy, angle, major, minor, ovm))
                            else:
                                raise Exception('Invalid option for cquad4; nid=%r ilayer=%r' % (nid, ilayer))
                            i += 1
                            eid2 = self.data.index[i][0]
                        #print("eid=%s eid2=%s" % (eid, eid2))
                        #assert eid != eid2, 'i=%s eid=%s eid2=%s msg=\n%s' % (i, eid, eid2,''.join(msg))
                    else:
                        raise NotImplementedError('CQUAD4 not-billinear')
                elif etype == 'CTRIA3':
                    eid2 = eid
                    while eid2 == eid and i < ndata:
                        (eid, nid, ilayer) = self.data.index[i]
                        index = self.data.index[i]
                        data = self.data.ix[index]

                        fd = data[kfd]
                        angle = data['angle']
                        oxx = data[koxx]
                        oyy = data[koyy]
                        txy = data[ktxy]
                        major = data[komax]
                        minor = data[komin]
                        ovm = data[kovm]

                        ([fd, oxx, oyy, txy, major, minor, ovm], isAllZeros) = writeFloats13E([fd, oxx, oyy, txy,
                                                                                               major, minor, ovm])
                        ([angle], isAllZeros) = writeFloats8p4F([angle])

                        # ..note:: nid==0 is the center of the CTRIA3
                        if ilayer == 1:
                            msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                        else:
                            msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip())
                        i += 1
                        eid2 = self.data.index[i][0]
                else:
                    raise NotImplementedError(etype)

            j += 1
            if etype != etype_old:
                etype_old = etype
                msg.append(pageStamp + str(pageNum) + '\n\n')
                f.write(''.join(msg))
                msg = ['']
                pageNum += 1
        except IndexError:
            pass

        if etype == etype_old:
            msg.append(pageStamp + str(pageNum) + '\n\n')
            f.write(''.join(msg))
            pageNum += 1
        return pageNum - 1


    def _write_f06_transient_helper(self, msgPacks, header, pageStamp, f, pageNum=1, is_mag_phase=False):

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        #(typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        return pageNum - 1

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
            raise NotImplementedError('eType = %r' % eType)  # CQUAD8, CTRIA6
        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        msg = ['']
        pageNum += 1
        return pageNum - 1


class PlateStressObject(RealPlateResults, StressObject):
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
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealPlateResults.__init__(self)
        StressObject.__init__(self, data_code, isubcase, read_mode)

    def getOrderedETypes(self, valid_types):
        """
        Groups element IDs by element type

        :param validTypes: list of valid element types
                           e.g. ['CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8']
        :returns TypesOut:      the ordered list of types
        :returns orderedETypes: dictionary of key=Type; value=list of IDs; to write
        """
        raise RuntimeError('removed...')
        ordered_etypes = {}
        types_out = []
        ordered_etypes = {}

        #print("valid_types =", valid_types)
        #validTypes = ['CTRIA3','CTRIA6','CQUAD4','CQUAD8']
        for etype in valid_types:
            ordered_etypes[etype] = []

        #eids = self.element_data['element_id']
        #print(self.element_data)
        eids = list(self.element_data.index)
        etypes = list(self.element_data['element_type'])
        #print("eids =", eids)
        for i in xrange(len(self.element_data)):
            #print("i =", i)
            etype = etypes[i]
            eid = eids[i]
            #print("eid =", eid, etype)
            assert etype in valid_types, 'unsupported eType=%s' % etype
            ordered_etypes[etype].append(eid)
            if etype not in types_out:
                types_out.append(etype)
        #print('ordered_etypes', ordered_etypes)
        return types_out, ordered_etypes

    def add_f06_data(self, data, transient):
        if transient is None:
            #print(data)

            eType = data[0][0]
            #print('eType = %s' % eType)

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

    def _get_headers(self):
        if self.isFiberDistance():
            fd = 'fd'  # 'fiber_distance'
        else:
            fd = 'fiber_curvature'
        return(fd, 'oxx', 'oyy', 'txy', 'omax', 'omin', 'ovm')


    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum, is_mag_phase)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % vonMises]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % vonMises]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % vonMises]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % vonMises]

        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None

        #eTypes = []
        #eTypes = self.get_element_types()
        #print(self.element_data)

        first_type_index = self._get_first_type_index()
        #print("first_type_index", first_type_index)
        #element_types = self.element_data.ix['element_type']

        #self.element_data.set_index(['element_id'])
        for element_type, eid in first_type_index.iteritems():
            if 'CQUAD4' == element_type:
                #print(self.data)
                #print('self.element_data =\n', self.element_data)
                #i_first_quad = element_types.index('CQUAD4')
                #i_first_quad = eid
                #eid = self.data['element_id'][i_first_quad]
                nnodes = self.element_data['nnodes'][eid]
                #print('nnodes =', nnodes)
                isBilinear = True
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
                if nnodes == 1:
                    isBilinear = False
                    quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp
            elif 'CTRIA3' == element_type:
                triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp
            else:
                raise NotImplementedError(element_type)

        if 'CQUAD8' in first_type_index:
            quad8Msg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in first_type_index:
            qkey = element_types .index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.data['oxx'][kkey]
            isBilinear = True
            quadrMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA6' in first_type_index:
            tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in first_type_index:
            trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        return self._write_f06_helper(msgPacks, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
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

        first_type_index = self._get_first_type_index()
        #print("first_type_index", first_type_index)
        #element_types = self.element_data.ix['element_type']

        #self.element_data.set_index(['element_id'])
        for element_type, eid in first_type_index.iteritems():
            if 'CQUAD4' == element_type:
                #print('self.element_data =\n', self.element_data)
                i_first_quad = eid
                nnodes = self.element_data['nnodes'][eid]
            elif 'CTRIA3' == element_type:
                triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp
            elif 'CTRIA6' == element_type:
                tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp
            elif 'CTRIAR' == element_type:
                trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp
            else:
                raise NotImplementedError(msg)

        eTypes = self.get_element_types()
        dts = self.oxx.keys()
        dt = dts[0]

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

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        return self._write_f06_transient_helper(msgPacks, header, pageStamp, f, pageNum, is_mag_phase)


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


class PlateStrainObject(RealPlateResults, StrainObject):
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
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealPlateResults.__init__(self)
        StrainObject.__init__(self, data_code, isubcase, read_mode)

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

    def _get_headers(self):
        if self.isFiberDistance():
            fd = 'fiber_distance'
        else:
            fd = 'fiber_curvature'
        return(fd, 'exx', 'eyy', 'exy', 'emax' , 'emin', 'evm')

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

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

        eTypes = self.get_element_types()

        first_type_index = self._get_first_type_index()
        #print("first_type_index", first_type_index)
        #element_types = self.element_data.ix['element_type']

        #self.element_data.set_index(['element_id'])
        isBilinear = None
        for element_type, eid in first_type_index.iteritems():
            if 'CQUAD4' == element_type:
                #print('self.element_data =\n', self.element_data)
                i_first_quad = eid
                nnodes = self.element_data['nnodes'][eid]
                quadMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
                if nnodes == 1:
                    isBilinear = False
                    quadMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

            elif 'CTRIA3' == element_type:
                triMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp
            elif 'CTRIA6' in eTypes:
                tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp
            elif 'CTRIAR' in eTypes:
                trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

            else:
                raise NotImplementedError('element_type = %r' % element_type)

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

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        return self._write_f06_helper(msgPacks, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
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


        first_type_index = self._get_first_type_index()
        #print('first_type_index =', first_type_index)
        for element_type, eid in first_type_index.iteritems():
            if element_type == 'CQUAD4':
                #print('self.element_data =\n', self.element_data)
                i_first_quad = eid
                nnodes = self.element_data['nnodes'][eid]
                quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
                if nnodes == 1:
                    isBilinear = False
                    quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp
            elif element_type == 'CTRIA3':
                triMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp
            elif element_type == 'CTRIA6':
                tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp
            elif element_type == 'CTRIAR':
                trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp
            else:
                raise NotImplementedError('element_type = %r' % element_type)

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        return self._write_f06_transient_helper(msgPacks, header, pageStamp, f, pageNum, is_mag_phase)


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