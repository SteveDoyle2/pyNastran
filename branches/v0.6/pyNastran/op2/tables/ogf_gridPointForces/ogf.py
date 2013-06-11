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
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack

# pyNastran
from pyNastran.op2.op2_helper import polar_to_real_imag
from .ogf_Objects import gridPointForcesObject, complexGridPointForcesObject


class OGF(object):
    """Table of Grid Point Forces"""

    def readTable_OGF(self):
        table3 = self.readTable_OGF_3
        table4Data = self.readOGF_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OGF()

    def _delete_attributes_OGF(self):
        params = ['format_code', 'appCode', 'num_wide', 'value1', 'value2', ]
        self._delete_attributes(params)

    #def add_data_parameter(self,data,Name,Type,FieldNum,applyNonlinearFactor=True):
    #    #self.mode = self.get_values(data, 'i', 5) ## mode number
    #    value = self.get_values(data, Type, FieldNum)
    #    setattr(self, Name, value)
    #    self.data_code[Name] = value
    #
    #    if applyNonlinearFactor:
    #        self.nonlinear_factor = value
    #        self.data_code['nonlinear_factor'] = value
    #        self.data_code['name'] = Name

    def apply_data_code_value(self, Name, value):
        self.data_code[Name] = value

    def readTable_OGF_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        if self.make_op2_debug:
            self.op2Debug.write('buffer_words=%s\n' % (str(buffer_words)))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack(b'i', data)
        data = self.get_data(4 * 50)
        #print self.print_block(data)

        (three) = self.parse_approach_code(data)
        ## format code
        self.add_data_parameter(data,'format_code', 'i', 9, False)
        ## approach code ???
        self.add_data_parameter(data,'appCode', 'i', 9, False)
        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.add_data_parameter(data,'num_wide', 'i', 10, False)
        ## Data Value 1
        self.add_data_parameter(data,'value1', 'i', 11, False)
        ## Data Value 2
        self.add_data_parameter(data,'value2', 'f', 12, False)

        #self.print_block(data) # on
        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            #self.extractDt = self.extractInt
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            self.apply_data_code_value('dataNames', ['mode'])
        #elif self.analysis_code==3: # differential stiffness
            ## load set number
            #self.lsdvmn = self.get_values(data,'i',5)
            #self.extractDt = self.extractInt
        #elif self.analysis_code==4: # differential stiffness
            #self.extractDt = self.extractInt
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'time', 'f', 5)
            self.apply_data_code_value('dataNames', ['time'])
        #elif self.analysis_code==7: # pre-buckling
            #self.extractDt = self.extractInt
            #self.apply_data_code_value('dataNames',['lsdvmn'])
        #elif self.analysis_code==8: # post-buckling
            #self.extractDt = self.extractInt
            #self.apply_data_code_value('dataNames',['lsdvmn','eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            #self.add_data_parameter(data,'eigr','f',6,False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load factor
            self.add_data_parameter(data, 'loadFactor', 'f', 5)
            self.apply_data_code_value('dataNames', ['loadFactor'])
        #elif self.analysis_code==11: # old geometric nonlinear statics
            #self.extractDt = self.extractInt
            #self.apply_data_code_value('dataNames',['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## time step
            self.add_data_parameter(data, 'time', 'f', 5)
            self.apply_data_code_value('dataNames', ['time'])
            #self.extractDt = self.extractInt
            #self.apply_data_code_value('dataNames',['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % (self.analysis_code)
            raise RuntimeError(msg)

        if not self.is_sort1():
            raise NotImplementedError('sort2...')

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)
        #self.print_block(data)
        self.read_title()

    def readOGF_Data(self):
        #print "self.analysis_code=%s table_code(1)=%s thermal(23)=%g" %(self.analysis_code,self.table_code,self.thermal)
        #tfsCode = [self.table_code,self.format_code,self.sort_code]
        #print self.data_code
        #print "tfsCode=%s" %(tfsCode)

        if self.table_code == 19:  # grid point forces
            assert self.table_name in ['OGPFB1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self.readOGF_Data_table19()
        else:
            self.not_implemented_or_skip('bad OGF table')

    def readOGF_Data_table19(self):  # grid point forces
        #is_sort1 = self.is_sort1()
        if self.num_wide == 10:  # real/random
            #if self.thermal==0:
            self.create_transient_object(self.gridPointForces,
                                       gridPointForcesObject)  # real
            self.handle_results_buffer(self.readOGF_numWide10,
                                      resultName='gridPointForces')
            #else:
                #self.not_implemented_or_skip()
            #self.handle_results_buffer(self.OUG_RealTable)
        elif self.num_wide == 16:  # real/imaginary or mag/phase
            #if self.thermal==0:
            self.create_transient_object(self.gridPointForces,
                                       complexGridPointForcesObject)  # complex
            self.handle_results_buffer(self.readOGF_numWide16,
                                      resultName='gridPointForces')
            #else:
                #self.not_implemented_or_skip()
            #self.handle_results_buffer(self.OUG_ComplexTable)
        else:
            raise NotImplementedError('only num_wide=10 or 16 is allowed  num_wide=%s' % (self.num_wide))

    def readOGF_numWide10(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i8s6f'
        format1 = bytes(format1)

        while len(self.data) >= 40:
            eData = self.data[0:4 * 10]
            self.data = self.data[4 * 10:]
            out = unpack(format1, eData)
            (eKey, eid, elemName, f1, f2, f3, m1, m2, m3) = out
            eKey = extract(eKey, dt)
            elemName = elemName.strip()
            #data = (eid,elemName,f1,f2,f3,m1,m2,m3)
            self.obj.add(dt, eKey, eid, elemName, f1, f2, f3, m1, m2, m3)
            #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%g f2=%g f3=%g m1=%g m2=%g m3=%g" %(ekey,eid,elemName,f1,f2,f3,m1,m2,m3)
        #print len(self.data)

    def readOGF_numWide16(self):
        dt = self.nonlinear_factor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i8s12f'
        format1 = bytes(format1)
        is_magnitude_phase = self.is_magnitude_phase()

        while len(self.data) >= 64:
            eData = self.data[0:4 * 16]
            self.data = self.data[4 * 16:]
            out = unpack(format1, eData)
            (eKey, eid, elemName, f1r, f2r, f3r, m1r, m2r, m3r,
                f1i, f2i, f3i, m1i, m2i, m3i) = out
            eKey = extract(eKey, dt)

            if is_magnitude_phase:
                f1 = polar_to_real_imag(f1r, f1i)
                m1 = polar_to_real_imag(m1r, m1i)
                f2 = polar_to_real_imag(f2r, f2i)
                m2 = polar_to_real_imag(m2r, m2i)
                f3 = polar_to_real_imag(f3r, f3i)
                m3 = polar_to_real_imag(m3r, m3i)
            else:
                f1 = complex(f1r, f1i)
                m1 = complex(m1r, m1i)
                f2 = complex(f2r, f2i)
                m2 = complex(m2r, m2i)
                f3 = complex(f3r, f3i)
                m3 = complex(m3r, m3i)

            elemName = elemName.strip()
            #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%s f2=%s f3=%s m1=%s m2=%s m3=%s" %(ekey,eid,elemName,f1r+f1i,f2r+f2i,f3r+f3i,m1r+m1i,m2r+m2i,m3r+m3i)
            self.obj.add(dt, eKey, eid, elemName, f1, f2, f3, m1, m2, m3)

    def readThermal4(self):
        #print self.code_information()
        #print self.print_block(self.data)
        n = 0
        nEntries = len(self.data) // 32
        for i in xrange(nEntries):
            eData = self.data[n:n + 32]
            out = unpack(b'2i6f', eData)
            #nid = (out[0]-self.device_code)//10  # TODO update...
            #print out
            n += 32
            #print "nid = ",nid
