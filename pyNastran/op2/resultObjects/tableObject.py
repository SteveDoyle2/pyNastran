#from struct import pack
import numpy as np
import numpy

import pandas as pd
from numpy import array, sqrt, abs, angle, zeros, ones  # dot,

from pyNastran import as_array
from pyNastran.op2.resultObjects.op2_Objects import scalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E

try:
    from pylab import xlabel, ylabel, show, grid, legend, plot, title
    import matplotlib.pyplot as plt
except ImportError:
    pass


class TableObject(scalarObject):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        self.shape = {}
        self._inode_start = None
        self._inode_end = None
        scalarObject.__init__(self, data_code, isubcase, read_mode)

        # new method, not finished
        self.nodeIDs_to_index = None
        self.gridTypes2 = {}

        #self.gridTypes = {}
        #self.translations = {}
        #self.rotations = {}
        
        self.dt = dt
        if is_sort1:
            if dt is not None:
                pass
        else:
            assert dt is not None
            self.add = self.add_sort2
            self.add_array = self.add_sort2

    def get_stats(self):
        ndt, nnodes, dts = self._get_shape()

        msg = self.get_data_code()
        if dts[0] is not None:
            msg.append('  type=%s ntimes=%s nnodes=%s\n'
                       % (self.__class__.__name__, ndt, nnodes))
        else:
            msg.append('  type=%s nnodes=%s\n' % (self.__class__.__name__, nnodes))
        msg.append('  T1, T2, T3, R1, R2, R3, gridTypes\n')
        return msg

    def isImaginary(self):
        return False

    def add_array_f06_data(self, data, transient):
        nnodes = len(data)

        nodeIDs_to_index = zeros(nnodes, dtype='int32')
        gridTypes = zeros(nnodes, dtype='string')
        translations = zeros((nnodes, 6), dtype='float32')

        for i, line in enumerate(data):
            (nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
            nodeIDs_to_index[i] = nodeID
            gridTypes[i] = gridType
            translations[i, :] = array([t1, t2, t3, r1, r2, r3])

        if transient is None:
            self.add_array(None, nodeIDs_to_index, gridTypes, translations)
            #print('%s gridTypes=%s'  % (self.__class__.__name__, gridTypes))
            return
        else:
            (dtName, dt) = transient
            self.data_code['name'] = dtName
            if dt not in self.translations:
                self.update_dt(self.data_code, dt, read_mode)
            self.add_array_sort1(dt, nodeIDs_to_index, gridTypes, translations)

    def add_f06_data(self, data, transient):
        if as_array:
            return self.add_array_f06_data(data, transient)

        #raise RuntimeError('this should be commented out')
        if transient is None:
            for line in data:
                (nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
                self.gridTypes[nodeID] = gridType
                self.translations[nodeID] = array([t1, t2, t3])
                self.rotations[nodeID] = array([r1, r2, r3])
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt, read_mode)

        for line in data:
            (nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[nodeID] = gridType
            self.translations[dt][nodeID] = array([t1, t2, t3])
            self.rotations[dt][nodeID] = array([r1, r2, r3])

    def update_dt(self, data_code, dt, read_mode):
        if read_mode == 1:
            return
        self.data_code = data_code
        self.apply_data_code()
        if dt is not None:
            print("updating %s...%s=%s  isubcase=%s" % (self.data_code['name'],
                self.data_code['name'], dt, self.isubcase))
            self.dt = dt

    def get_transients(self):
        if as_array:
            k = self.translations.keys()
        else:
            k = self.translations[:, 0, 0]  # itime, inode, icomponent
        k.sort()
        return k

    def _preallocate(self, dt, nnodes):
        if self.shape is None:
            self._inode_start += nnodes
            self._inode_end += nnodes
        else:
            ndt, nnodes, dts = self._get_shape()
            #print "ndt=%s nnodes=%s dts=%s" % (ndt, nnodes, str(dts))
            
            data = {}
            columns = []
            indexs = []
            if dts[0] is not None:
                if isinstance(dt, int):
                    data['dt'] = pd.Series(zeros((ndt * nnodes), dtype='int32'))
                else:
                    data['dt'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
                columns.append('dt')
                indexs = ['dt']
            data['node_id'] = pd.Series(zeros((ndt * nnodes), dtype='int32'))
            indexs.append('node_id')

            #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
            #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
            data['T1'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
            data['T2'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
            data['T3'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
            data['R1'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
            data['R2'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
            data['R3'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
            columns += ['node_id', 'T1', 'T2', 'T3', 'R1', 'R2', 'R3']

            self.data = pd.DataFrame(data, columns=columns)
            self._inode_start = 0
            self._inode_end = nnodes
        return self._inode_start, self._inode_end

    def _finalize(self):
        ndt, nnodes, dts = self._get_shape()
        
        mapper = {
            1 : 'G',
        }
        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append(mapper[grid_type])
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if dts[0] is not None:
            self.data = self.data.set_index(['dt', 'node_id'])
        else:
            self.data = self.data.set_index('node_id')
        #print "final\n", self.data
        del self._inode_start
        del self._inode_end

    def _increase_size(self, dt, nnodes):
        #self.shape += 1
        if dt in self.shape:  # default dictionary
            self.shape[dt] += nnodes
        else:
            self.shape[dt] = nnodes

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nnodes = self.shape[shape0]
        print "ndt=%s nnodes=%s dts=%s" % (ndt, nnodes, str(dts))
        return ndt, nnodes, dts

    def add_array(self, dt, nodeIDs_to_index, gridTypes):
        dts = self.shape.keys()
        if dts[0] == dt:
            for nid, gridType in zip(nodeIDs_to_index, gridTypes):
                self.gridTypes2[nid] = gridType
            #if self.gridTypes2 is None:
                #self.gridTypes2 = gridTypes
            #else:
                #self.gridTypes2 = numpy.concatenate((self.gridTypes2, gridTypes))
            print 'gridTypes2', self.gridTypes2

    def get_as_sort1(self):
        return (self.translations, self.rotations)

    def get_as_sort2(self):
        """returns translations and rotations in sort2 format"""
        translations2 = {}
        rotations2 = {}
        if self.dt is not None:
            #return self.__reprTransient__()

            for dt, translations in sorted(self.translations.iteritems()):
                nodeIDs = translations.keys()
                for nodeID in nodeIDs:
                    translations2[nodeID] = {}
                    rotations2[nodeID] = {}

            for dt, translations in sorted(self.translations.iteritems()):
                for nodeID, translation in sorted(translations.iteritems()):
                    rotation = self.rotations[dt][nodeID]
                    translations2[nodeID][dt] = translation
                    rotations2[nodeID][dt] = rotation
        else:
            return (self.translations, self.rotations)
            #for nodeID,translation in sorted(self.translations.iteritems()):
            #    rotation = self.rotations[nodeID]
            #    translations2[nodeID] = translation
            #    rotations2[nodeID]    = rotation
        return (translations2, rotations2)

    def _write_matlab(self, name, isubcase, f=None, is_mag_phase=False):
        """
        name = displacements
        """
        if self.nonlinear_factor is not None:
            self._write_matlab_transient(name, isubcase, f, is_mag_phase)
        #print "static!!!!"
        #msg = []
        #magPhase = 0
        #if is_mag_phase:
        #    magPhase = 1
        #msg.append('fem.%s.is_mag_phase = %s' %(name,magPhase))

        nodes = self.translations.keys()

        nodes.sort()
        #nNodes = len(nodes)
        self._write_matlab_args(name, isubcase, f)
        f.write('fem.%s(%i).nodes = %s;\n' % (name, isubcase, nodes))

        msgG = "fem.%s(%i).gridTypes    = ['" % (name, isubcase)
        msgT = 'fem.%s(%i).translations = [' % (name, isubcase)
        msgR = 'fem.%s(%i).rotations    = [' % (name, isubcase)
        spaceT = ' ' * len(msgT)
        spaceR = ' ' * len(msgR)
        i = 0
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]
            msgG += '%s' % (gridType)

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, isAllZeros) = writeFloats13E(vals)
            msgT += '[%s,%s,%s];' % (dx, dy, dz)
            msgR += '[%s,%s,%s];' % (rx, ry, rz)
            i += 1
            if i == 100:
                msgT += '\n%s' % (spaceT)
                msgR += '\n%s' % (spaceR)
                i = 0

        msgG += "'];\n"
        msgT += '];\n'
        msgR += '];\n'
        f.write(msgG)
        f.write(msgT)
        f.write(msgR)

    def _write_matlab_transient(self, name, isubcase, f=None, is_mag_phase=False):
        """
        name = displacements
        """
        times = self.translations.keys()
        nodes = self.translations[times[0]].keys()
        times.sort()
        nodes.sort()
        #nNodes = len(nodes)
        self._write_matlab_args(name, isubcase, f)
        dtName = '%s' % (self.data_code['name'])
        f.write('fem.%s(%i).nodes = %s;\n' % (name, isubcase, nodes))
        f.write('fem.%s(%i).%s = %s;\n' % (name, isubcase, dtName, times))

        msgG = "fem.%s(%i).gridTypes = ['" % (name, isubcase)

        for nodeID, gridType in sorted(self.gridTypes.iteritems()):
            msgG += '%s' % gridType
        msgG += "'];\n"
        f.write(msgG)
        del msgG

        nDt = len(self.translations)
        msgT = 'fem.%s(%i).translations.%s = cell(1,%i);\n' % (
            name, isubcase, dtName, nDt)
        msgR = 'fem.%s(%i).rotations.%s    = cell(1,%i);\n' % (
            name, isubcase, dtName, nDt)
        #msgT = 'fem.%s(%i).translations.%s = zeros(%i,3);\n' %(name,isubcase,dtName,n+1,nNodes)
        #msgR = 'fem.%s(%i).rotations.%s    = zeros(%i,3);\n' %(name,isubcase,dtName,n+1,nNodes)
        
        node_ids = self.data['node_id']
        for n, (dt, translations) in enumerate(sorted(self.translations.iteritems())):
            #msgT = 'fem.%s(%i).translations.%s(%i) = zeros(%i,3);\n' %(name,isubcase,dtName,n+1,nNodes)
            #msgR = 'fem.%s(%i).rotations.%s(%i)    = zeros(%i,3);\n' %(name,isubcase,dtName,n+1,nNodes)

            msgT += 'fem.%s(%i).translations.%s(%i) = [' % (
                name, isubcase, dtName, n + 1)
            msgR += 'fem.%s(%i).rotations.%s(%i)    = [' % (
                name, isubcase, dtName, n + 1)

            i = 0
            for nodeID in node_ids:
                translation = self.translations[dt][nodeID]
                rotation = self.rotations[dt][nodeID]
                #gridType = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation
                #vals = [dx,dy,dz,rx,ry,rz]
                msgT += '[%s,%s,%s];' % (dx, dy, dz)
                msgR += '[%s,%s,%s];' % (rx, ry, rz)
                i += 1
                if i == 100:
                    msgT += '\n'
                    msgR += '\n'
                    i = 0

            msgT += '];\n'
            msgR += '];\n'
            f.write(msgT)
            f.write(msgR)
            msgT = ''
            msgR = ''

    def _write_f06_block(self, words, header, pageStamp, pageNum=1, f=None):
        msg = words
        #assert f is not None # remove
        
        print self.data.to_string()
        node_ids = self.data['node_id']
        for nodeID in node_ids:
            translation = self.translations[nodeID]
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, isAllZeros) = writeFloats13E(vals)
            if not isAllZeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n'
                        % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def _write_f06_transient_block(self, words, header, pageStamp, pageNum=1, f=None):
        msg = []
        #assert f is not None # remove
        for dt, translations in sorted(self.translations.iteritems()):
            if isinstance(dt, float):  # fix
                header[1] = ' %s = %10.4E float %s\n' % (self.data_code[
                    'name'], dt, self.analysis_code)
            else:
                header[1] = ' %s = %10i integer %s\n' % (self.data_code[
                    'name'], dt, self.analysis_code)
            msg += header + words
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation
                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                if not isAllZeros:
                    [dx, dy, dz, rx, ry, rz] = vals2
                    msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def get_table_marker(self):
        if self.isATO():
            words = self.ATO_words()
        elif self.isCRM():
            words = self.CRM_words()
        elif self.isPSD():
            words = self.PSD_words()
        elif self.isRMS():
            words = self.RMS_words()
        elif self.isZERO():
            return self.ZERO_words()
        else:
            words = ['']
        return words

    def isATO(self):
        """Auto-Correlation Function"""
        if 'ATO' in self.table_name:
            return True
        return False

    def isCRM(self):
        """Correlated Root-Mean Square"""
        if 'CRM' in self.table_name:
            return True
        return False

    def isPSD(self):
        """Power Spectral Density"""
        if 'PSD' in self.table_name:
            return True
        return False

    def isRMS(self):
        """Root-Mean Square"""
        if 'RMS' in self.table_name:
            return True
        return False

    def isZERO(self):
        """Zero Crossings"""
        if 'NO' in self.table_name:
            return True
        return False

    def ATO_words(self):
        words = ['                                                 ( AUTO-CORRELATION FUNCTION )\n', ' \n']
        return words

    def CRM_words(self):
        words = ['                                               ( CUMULATIVE ROOT MEAN SQUARE )\n', ' \n']
        return words

    def PSD_words(self):
        words = ['                                             ( POWER SPECTRAL DENSITY FUNCTION )\n', ' \n']
        return words

    def RMS_words(self):
        words = ['                                                     ( ROOT MEAN SQUARE )\n', ' \n']
        return words

    def ZERO_words(self):
        words = ['                                                 ( NUMBER OF ZERO CROSSINGS )\n', ' \n']
        return words

    def plot(self, nodeList=None, resultType='Translation', coord=3, markers=None,
             Title=None, hasLegend=True, Legend=None, xLabel=None, yLabel=None, alphaLegend=0.5):
        """
        :param nodeList: a list of the node IDs to plot vs the
               independent variable (default=None; all nodes)
        :param resultType: the variable to plot ('Translation','Rotation')
        :param coord: the coordinate to plot (for <x,y,z>, x=0,y=1,z=2,Mag=3);
               default=Magnitude
        :param markers:  a list of colors/marker shapes for each line
        :param Title: title of the plot (default=the object name)
        :param hasLegend: should a legend be shown (default=False)
        :param Legend: the list of the legend titles (default=No Legend)
        :param xLabel: the name of the xAxis (default=the name of
               the independent variable; string)
        :param yLabel: the name of the xAxis (default=the name of
               the dependent variable; string)
        :param alphaLegend: the transparency of the legend;
               (0.0=solid; 1.0=transparent; default=0.5)

        .. todo:: fix alphaLegend; test options more...
        """
        (results, nodeList, markers, Title, xLabel, yLabel) = getPlotData(
            self, nodeList, resultType, coord, markers, Title, hasLegend,
            Legend, xLabel, yLabel)

        i = 0
        Labels = []
        #print "nodeList = ",nodeList
        node0 = nodeList[0]
        Xs = sorted(results[node0].keys())

        (fig, ax) = plt.subplots(1)
        #leg = plt.legend(loc='best', fancybox=True,alpha=0.5)
        for nodeID in nodeList:  # translation/rotation
            result = results[nodeID]
            Ys = []
            for dt, res in sorted(result.iteritems()):
                if coord == 3:
                    val = sqrt(res.dot(res))
                else:
                    val = res[coord]
                Ys.append(val)
            Label = 'Node %s' % (nodeID)
            Labels.append(Label)
            #ax.plot(Xs,Ys,markers[i],label=Label)
            #plot(Xs,Ys,Label)
            plot(Xs, Ys)
            i += 1
        plt.legend(Labels, loc='best', fancybox=True)
        #print dir(leg)
        #leg.get_frame().set_alpha(0.5)
        ax.set_title(Title)
        #ax.set_legend(Labels)
        #if hasLegend and Legend is None:
        #    plt.legend(Labels)
        #elif hasLegend:
        #    plt.legend(Legend)
        #plt.axes
        #print dir(plt)
        ax.grid(True)
        ax.set_ylabel(yLabel)
        ax.set_xlabel(xLabel)
        #alpha(alphaLegend)
        show()


class ComplexTableObject(scalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        scalarObject.__init__(self, data_code, isubcase)
        self.gridTypes = {}
        self.translations = {}
        self.rotations = {}

        self.dt = dt
        if is_sort1:
            pass
            #if dt is not None:
                #self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        ngrids = len(self.gridTypes)
        msg = self.get_data_code()

        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.translations)
            msg.append('  imaginary type=%s ntimes=%s ngrids=%s\n'
                       % (self.__class__.__name__, ntimes, ngrids))
        else:
            msg.append('  imaginary type=%s ngrids=%s\n'
                       % (self.__class__.__name__, ngrids))
        msg.append('  translations, rotations, gridTypes\n')
        return msg
            
    def isImaginary(self):
        return True

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                try:
                    (nodeID, gridType, v1, v2, v3, v4, v5, v6) = line
                except:
                    print('line = %r' % line)
                    raise
                self.gridTypes[nodeID] = gridType
                self.translations[self.dt][nodeID] = [v1, v2, v3]  # dx,dy,dz
                self.rotations[self.dt][nodeID] = [v4, v5, v6]  # rx,ry,rz
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt, read_mode)
        
        for line in data:
            try:
                (nodeID, gridType, v1, v2, v3, v4, v5, v6) = line
            except:
                print('line = %r' % line)
                raise
            self.gridTypes[nodeID] = gridType
            self.translations[self.dt][nodeID] = [v1, v2, v3]  # dx,dy,dz
            self.rotations[self.dt][nodeID] = [v4, v5, v6]  # rx,ry,rz

    def add_complex_f06_data(self, data, transient):
        raise NotImplementedError()

    def update_dt(self, data_code, dt, read_mode):
        self.data_code = data_code
        self.apply_data_code()

    def delete_transient(self, dt):
        del self.translations[dt]
        del self.rotations[dt]

    def get_transients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add(self, dt, out):
        (nodeID, gridType, v1, v2, v3, v4, v5, v6) = out
        #msg = "dt=%s nodeID=%s v1=%s v2=%s v3=%s" %(dt,nodeID,v1,v2,v3)
        #assert isinstance(nodeID,int),nodeID
        msg = "nodeID=%s v1=%s v2=%s v3=%s\n" % (nodeID, v1, v2, v3)
        msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print msg
        assert 0 < nodeID < 1000000000, msg  # -1
        assert -0.5 < dt, msg              # remove
        assert isinstance(nodeID, int), msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[nodeID] = self.recastGridType(gridType)
        self.translations[nodeID] = [v1, v2, v3]  # dx,dy,dz
        self.rotations[nodeID] = [v4, v5, v6]  # rx,ry,rz

    def _write_matlab_transient(self, name, isubcase, f=None, is_mag_phase=False):
        """
        name = displacements
        """
        #msg = []
        times = self.translations.keys()
        nodes = self.translations[times[0]].keys()
        times.sort()
        nodes.sort()
        #nNodes = len(nodes)
        self._write_matlab_args(name, isubcase, f)
        dtName = '%s' % (self.data_code['name'])
        f.write('fem.%s(%i).nodes = %s;\n' % (name, isubcase, nodes))
        f.write('fem.%s(%i).%s = %s;\n' % (name, isubcase, dtName, times))

        msgG = "fem.%s(%i).gridTypes = ['" % (name, isubcase)
        for nodeID, gridType in sorted(self.gridTypes.iteritems()):
            msgG += '%s' % gridType
        msgG += "'];\n"
        f.write(msgG)
        del msgG

        for n, (dt, translations) in enumerate(sorted(self.translations.iteritems())):
            #msgT = 'fem.%s(%i).translations.%s(%i) = zeros(%i,3);\n' %(name,isubcase,dtName,n+1,nNodes)
            #msgR = 'fem.%s(%i).rotations.%s(%i)    = zeros(%i,3);\n' %(name,isubcase,dtName,n+1,nNodes)

            msgT = 'fem.%s(%i).translations.%s(%i) = [' % (
                name, isubcase, dtName, n + 1)
            msgR = 'fem.%s(%i).rotations.%s(%i)    = [' % (
                name, isubcase, dtName, n + 1)

            i = 0
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                msgT += '[%s+%sj,%s+%sj,%s+%sj];' % (dx.real, dx.imag,
                                                     dy.real, dy.imag, dz.real, dz.imag)
                msgR += '[%s+%sj,%s+%sj,%s+%sj];' % (rx.real, rx.imag,
                                                     ry.real, ry.imag, rz.real, rz.imag)
                i += 1
                if i == 100:
                    msgT += '\n'
                    msgR += '\n'
                    i = 0

            msgT += '];\n'
            msgR += '];\n'
            f.write(msgT)
            f.write(msgR)

    def _write_f06_block(self, words, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        #words += self.getTableMarker()
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        msg = words
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation

            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, isAllZeros) = writeImagFloats13E(vals)
            [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi, dzi, rxi,
                ryi, rzi] = vals2
            msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
            msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %
                       ('', '', dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def _write_f06_transient_block(self, words, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        msg = []
        #assert f is not None
        for dt, translations in sorted(self.translations.iteritems()):
            #print "dt = ",dt
            #sys.stdout.flush()
            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for nodeID, translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeImagFloats13E(vals, is_mag_phase)
                [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi,
                    dzi, rxi, ryi, rzi] = vals2
                if not isAllZeros:
                    msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
                    msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % ('', '', dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def plot(self, nodeList=None, resultType='Translation', displayType='Real Imag', coord=3, markers=None,
             Title=None, hasLegend=True, Legend=None, xLabel=None, yLabel=None, alphaLegend=0.5):
        """
        :param nodeList: a list of the node IDs to plot vs the
               independent variable (default=None; all nodes)
        :param resultType: the variable to plot ('Translation','Rotation')
        :param displayType: 'Real Imag' or 'Mag Phase'
        :param coord: the coordinate to plot (for <x,y,z>, x=0,y=1,z=2,Mag=3);
               default=Magnitude
        :param markers: a list of colors/marker shapes for each line
        :param Title: title of the plot (default=the object name)
        :param hasLegend: should a legend be shown (default=False)
        :param Legend: the list of the legend titles (default=No Legend)
        :param xLabel: the name of the xAxis (default=the name of
               the independent variable; string)
        :param yLabel: the name of the xAxis (default=the name of
               the dependent variable; string)
        :param alphaLegend: the transparency of the legend;
               (0.0=solid; 1.0=transparent; default=0.5)

        .. todo:: fix alphaLegend; test options more...
        """
        displayType = displayType.title()
        (results, nodeList, markers, Title, xLabel, yLabel) = getPlotData(
            self, nodeList, resultType, coord, markers, Title, hasLegend, Legend, xLabel, yLabel)

        i = 0
        Labels = []
        node0 = nodeList[0]
        Xs = sorted(results[node0].keys())

        (fig, ax) = plt.subplots(1)
        leg = plt.legend(loc='best', fancybox=True, alpha=0.5)
        for nodeID in nodeList:  # translation/rotation
            result = results[nodeID]
            Ys = []
            for dt, res in sorted(result.iteritems()):
                if coord == 3:
                    val = sqrt(res.dot(res))
                else:
                    val = res[coord]

                mag = abs(val)
                phase = angle(val, deg=True)
                Ys.append(val)
            Label = 'Node %s' % (nodeID)
            Labels.append(Label)
            ax.plot(Xs, Ys, markers[i])
            #plot(Xs,Ys,Label)
            i += 1
        #print dir(leg)
        #leg.get_frame().set_alpha(0.5)
        ax.set_title(Title)
        #ax.set_legend(Labels)
        if hasLegend and Legend is None:
            plt.legend(Labels)
        elif hasLegend:
            plt.legend(Legend)
        #plt.axes
        #print dir(plt)
        ax.grid(True)
        ax.set_ylabel(yLabel)
        ax.set_xlabel(xLabel)
        #alpha(alphaLegend)
        show()


def getPlotData(obj, nodeList, resultType, coord, markers, Title, hasLegend, Legend, xLabel, yLabel):
    if Title is None:
        label = ''
        subtitle = ''
        if obj.label:
            label = ' - %s' % (obj.label)
        if obj.subtitle:
            subtitle = ' - %s' % (obj.subtitle)
        Title = '%s%s%s - Subcase %s' % (
            obj.__class__.__name__, label, subtitle, obj.isubcase)

    resultType = resultType.title()
    assert coord in [0, 1, 2,
                     3], 'invalid coord...options=[0,1,2,3].  coord=%s' % (coord)
    assert resultType in ['Translation', 'Rotation'], "invalid resultType...options=['Translation','Rotation']."

    if xLabel is None:
        xLabel = obj.data_code['name'].title()
    if yLabel is None:
        if   coord == 0:
            yLabel = 'X %s' % (resultType)
        elif coord == 1:
            yLabel = 'Y %s' % (resultType)
        elif coord == 2:
            yLabel = 'Z %s' % (resultType)
        elif coord == 3:
            yLabel = '%s (Magnitude)' % (resultType)
        else:
            raise RuntimeError('invalid coord...options=[0,1,2,3].  choice=%s' % (coord))

    (translations, rotations) = obj.getAsSort2()

    if resultType == 'Translation':
        results = translations
    else:
        results = rotations

    nodeListAll = results.keys()
    if nodeList is None:
        nodeList = nodeListAll

    if hasLegend and Legend is not None:
        assert len(nodeList) == len(Legend), 'len(nodeList)=%s len(legend)=%s' % (len(nodeList), len(legend))
    if markers is None:
        markers = ['-'] * len(nodeList)
    #print "subcase = ",obj.isubcase
    return (results, nodeList, markers, Title, xLabel, yLabel)