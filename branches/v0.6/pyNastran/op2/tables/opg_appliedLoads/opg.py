from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack

from .opg_Objects import AppliedLoadsObject  # ComplexAppliedLoadsObject
from .opg_loadVector import LoadVectorObject, ComplexLoadVectorObject, ThermalLoadVectorObject
from .opnl_forceVector import ForceVectorObject, ComplexForceVectorObject

# OGS table # TODO move this...
from ..ogf_gridPointForces.ogs_surfaceStresses import GridPointStressesObject, GridPointStressesVolumeObject


class OPG(object):
    """Table of element forces"""
    def readTable_OPG(self):
        table3 = self.readTable_OPG_3
        table4Data = self.readOPG_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OPG()

    def _delete_attributes_OPG(self):
        params = ['lsdvm', 'mode', 'eigr', 'eign', 'eigi', 'mode_cycle', 'freq',
                  'time', 'lftsfq', 'dLoadID', 'format_code', 'num_wide', 'oCode']
        self._delete_attributes(params)

    def readTable_OPG_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack(b'i', data)
        data = self.get_data(4 * 50)

        #self.print_block(data)

        aCode = self.get_block_int_entry(data, 1)
        #print "aCode = ",aCode
        self.parse_approach_code(data)
        #isubcase = self.get_values(data,'i',4)

        ## dynamic load set ID/random code
        self.add_data_parameter(data, 'dLoadID', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## undefined in DMAP...
        self.add_data_parameter(data, 'oCode', 'i', 11, False)
        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        #print "dLoadID(8)=%s format_code(9)=%s num_wide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.format_code,self.num_wide,self.oCode,self.thermal)
        if not self.is_sort1():
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eign', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        #elif self.analysis_code == 4: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'time', 'f', 5)
            self.apply_data_code_value('dataNames', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               (self.analysis_code))

        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set, Mode number

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)
        #print self.code_information()

        #self.print_block(data)
        self.read_title()

    def readOPG_Data(self):
        #print "self.analysis_code=%s table_code(1)=%s thermal(23)=%g" %(self.analysis_code,self.table_code,self.thermal)

        if self.table_code == 19:
            assert self.table_name in [None], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self.readOPG_Data_table19()  # grid point force balance
        elif self.table_code == 2:  # load vector
            assert self.table_name in ['OPG1', 'OPGV1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self.readOPG_Data_table2()
        elif self.table_code == 12:  # nonlinear force vector
            assert self.table_name in ['OPNL1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self.readOPG_Data_table12()
        elif self.table_code == 26:  # OGS1 - grid point stresses - surface
            assert self.table_name in ['OGS1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self.readOGS1_Data_table26()
        elif self.table_code == 27:  # OGS1 - grid point stresses - volume direct
            assert self.table_name in ['OGS1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self.readOGS1_Data_table27()

        #elif self.table_code == 28:  # OGS1- grid point stresses - principal
            #assert self.table_name in ['OGS1'],'table_name=%s table_code=%s' %(self.table_name,self.table_code)
            #self.readOGS1_Data_table28()
            #self.not_implemented_or_skip()

        #elif self.table_code == 35:  # OGS - Grid point stress discontinuities (plane strain)
            #self.not_implemented_or_skip()

        # OFMPF2M - does this belong here?
        #elif tfsCode==[51,3,3]:
        #    self.readOPG_Data_format3_sort3()

        # OSMPF2M - does this belong here?
        #elif tfsCode==[52,3,3]:
        #    self.readOPG_Data_format3_sort3()

        # OPMPF2M - does this belong here?
        #elif tfsCode==[53,3,3]:
        #    self.readOPG_Data_format3_sort3()

        # OLMPF2M - does this belong here?
        #elif tfsCode==[54,3,3]:
        #    self.readOPG_Data_format3_sort3()

        # OGMPF2M - does this belong here?
        #elif tfsCode==[55,3,3]:
        #    self.readOPG_Data_format3_sort3()
        else:
            self.not_implemented_or_skip('bad OPG table')

        #print self.obj

    def readOPG_Data_table2(self):  # Load Vector
        #is_sort1 = self.is_sort1()
        #print "********\n",self.code_information()
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'loadVectors'
                self.create_transient_object(self.loadVectors, LoadVectorObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            elif self.thermal == 1:
                resultName = 'thermalLoadVectors'
                self.create_transient_object(self.thermalLoadVectors,
                                           ThermalLoadVectorObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'loadVectors'
                self.create_transient_object(self.loadVectors, ComplexLoadVectorObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            raise RuntimeError('only num_wide=8 or 14 is allowed  num_wide=%s' %
                               (self.num_wide))

    def readOPG_Data_table12(self):  # Nonlinear Force Vector (in progress)
        #is_sort1 = self.is_sort1()
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'forceVectors'
                self.create_transient_object(self.forceVectors, ForceVectorObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'forceVectors'
                self.create_transient_object(self.forceVectors, ComplexForceVectorObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            msg = 'only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide)
            raise RuntimeError(msg)

    def readOPG_Data_table19(self):  # Applied Loads
        #is_sort1 = self.is_sort1()
        if self.num_wide == 8:  # real/random
            resultName = 'appliedLoads'
            if self.thermal == 0:
                self.create_transient_object(self.appliedLoads,
                                           AppliedLoadsObject)
                self.handle_results_buffer(self.readOPGForces, resultName)
            else:
                self.not_implemented_or_skip()
        #elif self.num_wide == 14:  # real/imaginary or mag/phase
        #    if self.thermal == 0:
        #        self.create_transient_object(self.appliedLoads,ComplexAppliedLoadsObject) # complex
        #        raise NotImplementedError('can this use a OUG_Complex table???')
        #    else:
        #        raise NotImplementedError(self.code_information())
        #    #self.handle_results_buffer(self.OUG_ComplexTable,resultName)
        #    raise NotImplementedError(self.code_information())
        else:
            msg = 'only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide)
            raise RuntimeError(msg)

    def readOGS1_Data_table26(self):  # OGS1 - grid point stresses - surface
        #is_sort1 = self.is_sort1()
        resultName = 'gridPointStresses'
        if self.num_wide == 11:  # real/random
            self.create_transient_object(self.gridPointStresses,
                                       GridPointStressesObject)
            self.handle_results_buffer(
                self.readOGS1_table26_numWide11, resultName)
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % (self.num_wide)
            raise RuntimeError(msg)

    def readOGS1_table26_numWide11(self):  # surface stresses
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i4s8f'
        format1 = bytes(format1)

        while len(self.data) >= 44:
            eData = self.data[0:44]
            self.data = self.data[44:]  # 11*4
            out = unpack(format1, eData)
            (eKey, eid, fiber, nx, ny, txy, angle, major,
                minor, tmax, ovm) = out
            eKey = extract(eKey, dt)
            fiber = fiber.decode('utf-8').strip()
            self.obj.add(dt, eKey, eid, fiber, nx, ny, txy,
                         angle, major, minor, tmax, ovm)
        #print len(self.data)

    def readOGS1_Data_table27(self):  # OGS1 - grid point stresses - volume direct
        #is_sort1 = self.is_sort1()
        #print(self.code_information())
        if self.num_wide == 9:  # real/random
            resultName = 'gridPointVolumeStresses'
            self.create_transient_object(self.gridPointVolumeStresses,
                                       GridPointStressesVolumeObject)
            self.handle_results_buffer(
                self.readOGS1_table27_numWide9, resultName)
        else:
            msg = 'only num_wide=9 is allowed  num_wide=%s' % (self.num_wide)
            raise RuntimeError(msg)

    def readOGS1_table27_numWide9(self):  # surface stresses
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i7f'
        format1 = bytes(format1)

        while len(self.data) >= 36:
            eData = self.data[0:36]
            self.data = self.data[36:]  # 9*4
            out = unpack(format1, eData)
            (eKey, nx, ny, nz, txy, tyz, txz, pressure, ovm) = out
            eKey = extract(eKey, dt)

            self.obj.add(dt, eKey, nx, ny, nz, txy, tyz, txz, pressure, ovm)
        #print len(self.data)

    def readOPGForces(self):
        """
        .. todo:: needs some work...
        """
        raise RuntimeError('this should never been called...')
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i'
        format1 = bytes(format1)

        #nTotal = self.num_wide*4
        nTotal = 40  # same as dn
        while len(data) > dn:
            #print "len(data) = ",len(data)
            eData = self.data[0:dn]
            #self.print_block(data[:dn])
            (grid_device, eid) = unpack(format1, data[0:8])
            nodeID = extract(grid_device, dt)

            source = unpack(b'8s', data[8:16])
            (dx, dy, dz, rx, ry, rz) = unpack(b'6f', data[16:40])
            #print "source = |%s|" %(source)

            #print "nodeID=%s eid=%s source=|%s| dx=%-4i dy=%-4i dz=%-4i rx=%-4i ry=%-4i rz=%-4i" %(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            source2 = source.replace('*', '').replace('-', '').strip()
            assert source2.isalnum(), 'source=|%s| contains invalid characters...' % (source)

            self.obj.add(nodeID, eid, source, dx, dy, dz, rx, ry, rz)
            #print "grid_device = ",grid_device
            #print "device_code = ",device_code
            #print "nodeID=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(nodeID,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            self.data = self.data[dn:]

        #print "***********"
        #print self.obj
