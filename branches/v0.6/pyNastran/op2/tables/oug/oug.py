## GNU Lesser General Public License
##
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
##
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
##
## This file is part of pyNastran.
##
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
##
# http://www.cadfamily.com/online-help/I-DEAS/SDRCHelp/LANG/English/slv_ug/NAS_results_imported.htm
import sys
from struct import unpack

from pyNastran.op2.tables.oug.oug_displacements import (
    DisplacementObject,              # table_code=1     format_code=1 sort_code=0
    ComplexDisplacementObject)       # analysis_code=5  format_code=3 sort_code=1

# table_code=10 format_code=1 sort_code=0
from pyNastran.op2.tables.oug.oug_velocities import (
    VelocityObject, ComplexVelocityObject)

# table_code=11 format_code=1 sort_code=0
from pyNastran.op2.tables.oug.oug_accelerations import (
    AccelerationObject, ComplexAccelerationObject)

# table_code=1 format_code=1 sort_code=0
from pyNastran.op2.tables.oug.oug_temperatures import (
    TemperatureObject)

from pyNastran.op2.tables.oug.oug_eigenvectors import (
     EigenVectorObject,                     # analysis_code=2, sort_code=0 format_code   table_code=7
     ComplexEigenVectorObject,              # analysis_code=5, sort_code=1 format_code=1 table_code=7
    #RealEigenVectorObject,                 # analysis_code=9, sort_code=1 format_code=1 table_code=7
)
from pyNastran.op2.tables.opg_appliedLoads.opg_loadVector import ThermalVelocityVectorObject
from pyNastran.op2.op2_helper import polar_to_real_imag


class OUG(object):
    """Table of displacements/velocities/acceleration/heat flux/temperature"""

    def readTable_OUG2(self):  # OUGPSD
        #self.table_name = 'OUG'
        table3 = self.readTable_OUG2_3
        table4Data = self.readOUG2_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OUG()

        sys.exit('stopping...')

    def readOUG2_Data(self):
        self.readOUG_Data_table1(debug=True)

    def readTable_OUG2_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        if self.make_op2_debug:
            self.op2Debug.write('buffer_words=%s\n' % (str(buffer_words)))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        data = self.get_data(4 * 50)
        #print self.print_block(data)

        (three) = self.parse_approach_code2(data)

        ## random code
        self.add_data_parameter(data, 'randomCode', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)
        ## number of words per entry in record; .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## acoustic pressure flag
        self.add_data_parameter(data, 'acousticFlag', 'f', 13, False)
        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)
        self.isFlipped = False

        eidDevice = self.get_values(data, 'i', 5)
        floatVal = self.get_values(data, 'f', 5)
        eid = (eidDevice-self.device_code)//10
        #print("EID = %s" %(eidDevice))
        #print("floatVal = %s" %(floatVal))

        if self.table_name == 'OUGRMS2' and self.analysis_code == 1:
            self.add_data_parameter(data, 'nodeID', 'i',
                                  5, fixDeviceCode=True)  # frequency
            self.apply_data_code_value('dataNames', ['nodeID'])
            #print("nodeID = %s" %(self.nodeID))

        #self.isRegular = False
        elif self.analysis_code in [1, 5]:  # 5 # freq
            self.add_data_parameter(data, 'nodeID', 'i',
                                  5, fixDeviceCode=True)  # frequency
            self.apply_data_code_value('dataNames', ['nodeID'])
            #print("nodeID = %s" %(self.nodeID))
            #sys.exit(self.nodeID)
        elif self.analysis_code == 6:  # transient dt
            self.add_data_parameter(data, 'nodeID', 'i',
                                  5, fixDeviceCode=True)  # time step
            self.apply_data_code_value('dataNames', ['nodeID'])
        elif self.analysis_code == 10:  # freq/time step fqts
            self.add_data_parameter(data, 'nodeID', 'i', 5, fixDeviceCode=True)  # frequency / time step
            self.apply_data_code_value('dataNames', ['nodeID'])
        else:
            self.isRegular = True
            self.add_data_parameter(data, 'nodeID', 'i',
                                  5, fixDeviceCode=True)  # node ID
            self.apply_data_code_value('dataNames', ['nodeID'])
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)
        print(self.code_information())

        #self.print_block(data)
        self.read_title()


    def readTable_OUG(self):
        #self.table_name = 'OUG'
        table3 = self.readTable_OUG_3
        table4Data = self.readOUG_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OUG()

    def _delete_attributes_OUG(self):
        params = ['lsdvm', 'mode', 'eigr', 'mode_cycle', 'freq', 'dt', 'lftsfq',
                  'thermal', 'randomCode', 'fCode', 'num_wide', 'acousticFlag']
        self._delete_attributes(params)

    def readTable_OUG_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        if self.make_op2_debug:
            self.op2Debug.write('buffer_words=%s\n' % (str(buffer_words)))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        data = self.get_data(4 * 50)
        #print self.print_block(data)

        (three) = self.parse_approach_code(data)

        ## random code
        self.add_data_parameter(data, 'randomCode', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)
        ## number of words per entry in record; .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## acoustic pressure flag
        self.add_data_parameter(data, 'acousticFlag', 'f', 13, False)
        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)
        self.isFlipped = False
        if self.is_sort1():
            ## assuming tCode=1
            if self.analysis_code == 1:   # statics / displacement / heat flux
                self.add_data_parameter(data, 'lsdvmn', 'i',
                                      5, False)  # load set number
                self.apply_data_code_value('dataNames', ['lsdvmn'])
                self.setNullNonlinearFactor()
            elif self.analysis_code == 2:  # real eigenvalues
                self.add_data_parameter(data, 'mode',
                                      'i', 5)  # mode number
                self.add_data_parameter(data, 'eigr',
                                      'f', 6, False)  # real eigenvalue
                self.add_data_parameter(data, 'mode_cycle', 'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
                self.apply_data_code_value('dataNames', [
                    'mode', 'eigr', 'mode_cycle'])
            #elif self.analysis_code==3: # differential stiffness
                #self.lsdvmn = self.get_values(data,'i',5) ## load set number
                #self.data_code['lsdvmn'] = self.lsdvmn
            #elif self.analysis_code==4: # differential stiffness
                #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            elif self.analysis_code == 5:   # frequency
                self.add_data_parameter(data, 'freq', 'f', 5)  # frequency
                self.apply_data_code_value('dataNames', ['freq'])
            elif self.analysis_code == 6:  # transient
                self.add_data_parameter(data, 'dt', 'f', 5)  # time step
                self.apply_data_code_value('dataNames', ['dt'])
            elif self.analysis_code == 7:  # pre-buckling
                self.add_data_parameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.apply_data_code_value('dataNames', ['lsdvmn'])
            elif self.analysis_code == 8:  # post-buckling
                self.add_data_parameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.add_data_parameter(data, 'eigr', 'f',
                                      6, False)  # real eigenvalue
                self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
            elif self.analysis_code == 9:  # complex eigenvalues
                self.add_data_parameter(data, 'mode', 'i', 5)  # mode number
                self.add_data_parameter(data, 'eigr', 'f', 6,
                                      False)  # real eigenvalue
                self.add_data_parameter(data, 'eigi', 'f', 7,
                                      False)  # imaginary eigenvalue
                self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
            elif self.analysis_code == 10:  # nonlinear statics
                self.add_data_parameter(data, 'lftsfq', 'f', 5)  # load step
                self.apply_data_code_value('dataNames', ['lftsfq'])
            elif self.analysis_code == 11:  # old geometric nonlinear statics
                self.add_data_parameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.apply_data_code_value('dataNames', ['lsdvmn'])
            elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
                self.add_data_parameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.apply_data_code_value('dataNames', ['lsdvmn'])
            else:
                raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                                   (self.analysis_code))
        else:  # sort2
            eidDevice = self.get_values(data, 'i', 5)
            floatVal = self.get_values(data, 'f', 5)
            #eid = (eidDevice-self.device_code)//10
            #print("EID = %s" %(eidDevice))
            #print("floatVal = %s" %(floatVal))

            if self.table_name == 'OUGRMS2' and self.analysis_code == 1:
                self.add_data_parameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # frequency
                self.apply_data_code_value('dataNames', ['nodeID'])

            #self.isRegular = False
            elif self.analysis_code in [1, 5]:  # 5 # freq
                self.add_data_parameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # frequency
                self.apply_data_code_value('dataNames', ['nodeID'])
                #print("nodeID = %s" %(self.nodeID))
                #sys.exit(self.nodeID)
            elif self.analysis_code == 6:  # transient dt
                self.add_data_parameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # time step
                self.apply_data_code_value('dataNames', ['nodeID'])
            elif self.analysis_code == 10:  # freq/time step fqts
                self.add_data_parameter(data, 'nodeID', 'i', 5, fixDeviceCode=True)  # frequency / time step
                self.apply_data_code_value('dataNames', ['nodeID'])
            else:
                self.isRegular = True
                self.add_data_parameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # node ID
                self.apply_data_code_value('dataNames', ['nodeID'])
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)
        #print(self.code_information())

        if not self.is_sort1():
            raise NotImplementedError('sort2...')

        #self.print_block(data)
        self.read_title()

    def getOUG_FormatStart(self):
        """
        Returns an i or an f depending on if it's SORT2 or not.
        Also returns an extraction function that is called on the first argument
        """
        is_sort1 = self.is_sort1()
        #print('is_sort1 = %s' %is_sort1)
        if self.table_name == 'OUGRMS2':
            if self.analysis_code == 1:
                format1 = 'f'  # SORT2
                extract = self.extractSort2
            else:
                raise KeyError('invalid analysis_code...analysis_code=%s' %
                               (self.analysis_code))
        else:
            if is_sort1:
                #print "SORT1 - %s" %(self.get_element_type(self.element_type))
                #print "SORT1"
                format1 = 'i'  # SORT1
                extract = self.extractSort1
                #if self.analysis_code in [5]:
                    #extract==self.extractSort2
            else:  # values from IDENT   #: .. todo:: test this...
                #print "SORT2"
                #print "SORT2 - %s" %(self.get_element_type(self.element_type))
                if self.analysis_code in [2, 3, 4, 6, 7, 8, 11]:
                    format1 = 'f'  # SORT2
                    extract = self.extractSort2
                elif self.analysis_code in [5]:
                    format1 = 'f'
                    extract = self.extractSort2
                elif self.analysis_code in [1, 9, 10, 12]:
                    format1 = 'f'  # SORT1
                    extract = self.extractSort2
                else:
                    raise KeyError('invalid analysis_code...analysis_code=%s' %
                                   (self.analysis_code))
                #eid = self.nonlinear_factor
        return (format1, extract)

    def readOUG_Data(self):
        #print "self.analysis_code=%s table_code(1)=%s thermal(23)=%g" %(self.analysis_code,self.table_code,self.thermal)
        #tfsCode = [self.table_code,self.format_code,self.sort_code]

        #print self.data_code
        #print "tfsCode=%s" %(tfsCode)

        if self.table_code == 1 and self.table_name in ['OUGV1', 'OUG1', 'OUPV1']:    # displacement
            if self.table_name in ['OUGV1', 'OUG1']:
                assert self.table_name in ['OUGV1', 'OUG1'], 'table_name=%s table_code=%s\n%s' % (self.table_name, self.table_code, self.code_information())
                self.readOUG_Data_table1()
            else:  # 'OUPV1'
                msg = ('bad approach_code=%s, table_code=%s, format_code-%s '
                       'sort_code=%s on %s-OUG table' % (self.analysis_code,
                        self.table_code, self.format_code, self.sort_code,
                        self.table_name))
                self.not_implemented_or_skip(msg)
        elif self.table_code == 1 and self.table_name in ['OUGATO2', 'OUGCRM2', 'OUGPSD2', 'OUGRMS2', 'OUGNO2', ]:    # displacement
            #assert self.table_name in ['OUGATO2','OUGCRM2','OUGPSD2','OUGRMS2','OUGNO2',],'table_name=%s table_code=%s\n%s' %(self.table_name,self.table_code,self.code_information())
            self.readOUG_Data_table1()
        elif self.table_code == 7:  # modes
            if self.table_name in ['OUGV1', 'OUG1']:
                assert self.table_name in ['OUGV1', 'OUG1'], 'table_name=%s table_code=%s\n%s' % (self.table_name, self.table_code, self.code_information())
                self.readOUG_Data_table7()
            else:
                self.not_implemented_or_skip('bad OUG table')
        elif self.table_code == 10:  # velocity
            if self.table_name in ['OUGV1', 'OUG1']:
                assert self.table_name in ['OUGV1', 'OUG1'], 'table_name=%s table_code=%s\n%s' % (self.table_name, self.table_code, self.code_information())
                self.readOUG_Data_table10()
            else:
                self.not_implemented_or_skip('bad OUG table')
        elif self.table_code == 11:  # Acceleration vector
            if self.table_name in ['OUGV1', 'OUG1']:
                assert self.table_name in ['OUGV1', 'OUG1'], 'table_name=%s table_code=%s\n%s' % (self.table_name, self.table_code, self.code_information())
                self.readOUG_Data_table11()
            else:
                self.not_implemented_or_skip('bad OUG table')
        else:
            self.not_implemented_or_skip('bad OUG table')
        #print self.obj

    def readThermal4(self):  # used on self.thermal in [2,4,8]:
        #print self.code_information()
        #print self.print_block(self.data)
        n = 0
        nEntries = len(self.data) // 32
        for i in xrange(nEntries):
            eData = self.data[n:n + 32]
            out = unpack('2i6f', eData)
            #nid = (out[0]-self.device_code)//10    # TODO fix the device_code

            #print out
            n += 32
            #print "nid = ",nid
        #sys.exit('thermal4...')

    def readOUG_Data_table1(self,debug=False):  # displacement / temperature OUGV1, OUPV1
        """
        OUGV1   - global coordinate system in sort 1
        OUPV1   - scaled response spectra in sort 1
        OUGPSD2 - PSD in sort 2
        OUGATO2 - auto-correlated in sort 2
        """
        self.debug = debug
        isSkip = False
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                #print self.data_code
                if self.table_name in ['OUGV1', 'OUG1']:
                    resultName = 'displacements'
                    self.create_transient_object(self.displacements,
                                               DisplacementObject)
                elif self.table_name in ['OUGATO2']:
                    resultName = 'displacementsATO'
                    self.create_transient_object(self.displacementsATO,
                                               DisplacementObject)
                elif self.table_name in ['OUGCRM2']:
                    resultName = 'displacementsCRM'
                    self.create_transient_object(self.displacementsCRM,
                                               DisplacementObject)
                elif self.table_name in ['OUGPSD2']:
                    resultName = 'displacementsPSD'
                    self.create_transient_object(self.displacementsPSD,
                                               DisplacementObject)
                elif self.table_name in ['OUGRMS2']:
                    resultName = 'displacementsRMS'
                    self.create_transient_object(self.displacementsRMS,
                                               DisplacementObject)
                elif self.table_name in ['OUGNO2']:
                    resultName = 'displacementsNO'
                    self.create_transient_object(self.displacementsNO,
                                               DisplacementObject)
                else:
                    isSkip = True
                    self.not_implemented_or_skip('***table=%s***\n%s' % (self.table_name, self.code_information()))
                if not isSkip:
                    self.handle_results_buffer(self.OUG_RealTable, resultName)
            elif self.thermal == 1:
                resultName = 'temperatures'
                self.create_transient_object(
                    self.temperatures, TemperatureObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            #elif self.thermal == 8:
                #resultName = 'scaledDisplacements'
                #self.create_transient_object(self.scaledDisplacements,displacementObject)
                #self.handle_results_buffer(self.OUG_RealTable,resultName)
            else:
                self.not_implemented_or_skip('***thermal=%s***\n%s' %
                                          (self.thermal, self.code_information()))
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'displacements'
                self.create_transient_object(self.displacements,
                                           ComplexDisplacementObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip('only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide))

    def readOUG_Data_table7(self):  # eigenvector
        #is_sort1 = self.is_sort1()
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'eigenvectors'
                self.create_transient_object(self.eigenvectors,
                                           EigenVectorObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'eigenvectors'
                self.create_transient_object(self.eigenvectors,
                                           ComplexEigenVectorObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip('only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide))

    def readOUG_Data_table10(self):  # velocity
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'velocities'
                self.create_transient_object(
                    self.velocities, VelocityObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            elif self.thermal == 1:
                resultName = 'velocities'
                self.create_transient_object(self.velocities,
                                           ThermalVelocityVectorObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'velocities'
                self.create_transient_object(self.velocities,
                                           ComplexVelocityObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip('only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide))

    def readOUG_Data_table11(self):  # acceleration
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'accelerations'
                self.create_transient_object(self.accelerations,
                                           AccelerationObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'accelerations'
                self.create_transient_object(self.accelerations,
                                           ComplexAccelerationObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip('only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide))

    def OUG_RealTable(self):
        print("---starting OUG_RealTable---")
        dt = self.nonlinear_factor
        (format1, extract) = self.getOUG_FormatStart()
        #format1 = 'f'
        format1 += 'i6f'
        #print("format1 =",format1)

        nnodes = len(self.data) // 32
        istart = 0
        iend = 32
        for inode in xrange(nnodes):
            #print('istart=%i iend=%i' % (istart, iend))
            print('inode=%i nnodes=%i' % (inode, nnodes))
            eData = self.data[istart:iend]

            out = unpack(format1, eData)
            (eid, gridType, tx, ty, tz, rx, ry, rz) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, gridType, tx, ty, tz, rx, ry, rz]
            print('eid=%g gridType=%g tx=%g ty=%g tz=%g rx=%g ry=%g rz=%g' %(eid2, gridType, tx, ty, tz, rx, ry, rz))
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #print "%s" %(self.table_name),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            istart = iend
            iend += 32
        iend2 = 32 * (nnodes + 1)
        print("iend=%i iend2=%i" % (iend, iend2))
        self.data = self.data[iend:]
        print("len(self.data)=%i" % len(self.data))

    def OUG_ComplexTable(self):
        dt = self.nonlinear_factor

        (format1, extract) = self.getOUG_FormatStart()

        format1 += 'i12f'
        #print "format1 = ",format1
        is_magnitude_phase = self.is_magnitude_phase()

        nnodes = len(self.data) // 56
        #print('nnodes =', nnodes)
        istart = 0
        iend = 56  # 14 * 4
        for i in range(nnodes):
            eData = self.data[istart:iend]

            out = unpack(format1, eData)
            (eid, gridType, txr, tyr, tzr, rxr, ryr, rzr,
             txi, tyi, tzi, rxi, ryi, rzi) = out

            if is_magnitude_phase:
                tx = polar_to_real_imag(txr, txi)
                rx = polar_to_real_imag(rxr, rxi)
                ty = polar_to_real_imag(tyr, tyi)
                ry = polar_to_real_imag(ryr, ryi)
                tz = polar_to_real_imag(tzr, tzi)
                rz = polar_to_real_imag(rzr, rzi)
            else:
                tx = complex(txr, txi)
                rx = complex(rxr, rxi)
                ty = complex(tyr, tyi)
                ry = complex(ryr, ryi)
                tz = complex(tzr, tzi)
                rz = complex(rzr, rzi)

            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, gridType, tx, ty, tz, rx, ry, rz]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            istart = iend
            iend += 56
        self.data = self.data[iend:]

