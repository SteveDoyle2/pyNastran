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

from struct import unpack

from pyNastran.op2.tables.oes_stressStrain.real.elementsStressStrain import RealElementsStressStrain

class OESNLXR(RealElementsStressStrain):  ##todo real or complex?? see r767
    """Table of stresses/strains"""

    def readTable_OESNLXR(self):
        table3 = self.readTable_OESNLXR_3
        table4Data = self.readTable_OES_4_Data
        self.read_results_table(
            table3, table4Data, flag=1)  # flag=1 defines old style
        self._delete_attributes_OES()

    def _delete_attributes_OESNLXR(self):
        params = ['s_code', 'element_type', 'obj', 'markerStart', 'load_set',
                  'format_code', 's_code', 'thermal', 'lsdvmn', 'mode', 'eign',
                  'mode_cycle', 'freq', 'mode', 'eigr', 'eigi', 'dt']
        self._delete_attributes(params)

    def readTable_OESNLXR_3(self, iTable):
        #print "*iTable3 = ",iTable
        #if 0:
            #markers = self.read_markers([0,2])
            #print "markers=%s" %(markers)
            #block = self.read_block()
            #print "block = ",block
            #markers = self.read_markers([-1,7])
            #print "markers=%s" %(markers)
            #print self.print_section(200)

        buffer_words = self.get_buffer_words()

        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        if self.make_op2_debug:
            self.op2Debug.write('buffer_size=|%s|\n' % str(buffer_size))

        data = self.get_data(4 * 50)
        #self.print_block(data)
        if self.make_op2_debug:
            self.op2Debug.write('block3header\n')

        self.parse_approach_code(data)  # 3
        ## element type
        self.add_data_parameter(data, 'element_type', 'i', 3, False)
        ## load set ID
        self.add_data_parameter(data, 'load_set', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)
        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## stress/strain codes
        self.add_data_parameter(data, 's_code', 'i', 11, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        #print "loadset=%s format_code=%s numWordsEntry=%s s_code=%s" %(self.load_set,self.format_code,self.num_wide,self.s_code)
        #print "thermal(23)=%s element_type(3)=%s" %(self.thermal,self.element_type)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle
            # TODO confused on the type - F1???
            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
        #elif self.analysis_code==3: # differential stiffness
        #    ## load set number
        #    self.lsdvmn = self.get_values(data,'i',5)
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number

        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'dt', 'f', 5)
            #print "DT(5)=%s" %(self.dt)
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            #print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            #print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            #print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.add_data_parameter(data, 'dt', 'f', 5)
        else:
            raise RuntimeError('Invalid Analysis Code: analysis_code=%s' % (self.analysis_code))

        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        self.read_title()
        #print "n4 = ",self.n

    def readTable_OESNLXR_4_Data(self, iTable):
        isTable4Done = False
        isBlockDone = False
        #print self.print_section(100)

        data = self.get_data(16)
        #print self.print_block(data) # on
        #print "16 block..."
        #self.print_block(data)
        #self.print_block(data)
        buffer_words, = unpack('i', data[4:8])
        #print "buffer_words = ",buffer_words
        if self.make_op2_debug:
            self.op2Debug.write('buffer_words=|%s|\n' % (str(buffer_words)))

        #print "*********************"
        #buffer_words = self.get_marker() # 87 - buffer
        #print "OES4 buffer_words = ",buffer_words,buffer_words*4
        #self.verifyBufferSize(buffer_words)

        isBlockDone = not(buffer_words)
        #print "self.firstPass = ",self.firstPass

        ## table -4 is done, restarting table -3
        if self.isBufferDone:  # table is done when the buffer is done
            isTable4Done = True
            #print "exitA"
            return isTable4Done, isBlockDone
        if buffer_words == 0:
            #print "buffer_words 0 - done with Table4"
            isTable4Done = True
            #isBlockDone  = True
            self.print_section(40)
            #print "exitB"
            return isTable4Done, isBlockDone

        self.readOESNLXR_ElementTable()
        return isTable4Done, isBlockDone

    def readOESNLXR_ElementTable(self):
        #print "**self.readElementTable"
        #print "*element_type = ",self.element_type
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)

        self.rewind(4)
        self.data = self.read_block()  # 348
        #print "len(self.data) = ",len(self.data)

        if self.make_op2_debug:
            self.op2Debug.write('reading big data block\n')
        #print self.print_block(self.data)

        #msg = 'element_type=%s -> %s' %(self.element_type,self.get_element_type(self.element_type))
        tfsCode = [self.table_code, self.format_code, self.sort_code]

        self.parse_stress_code()

        if not self.is_valid_subcase():  # lets the user skip a certain subcase
            self.log.debug("***skipping table=%s isubcase=%s" %
                           (self.table_name, self.isubcase))
            self.skipOES_Element()
        elif self.thermal == 0:
            # Stress / Strain
            self.data_code['element_name'] = self.get_element_type(self.element_type)
            if   tfsCode == [5, 1, 0]:
                self.readOESNLXR_Data_format1_sort0()
            else:
                raise RuntimeError('invalid atfsCode=%s' % (self.atfsCode))
                self.skipOES_Element()
                pass
        else:
            raise RuntimeError('invalid thermal option...')

        #print self.obj

    def readOESNLXR_Data_format1_sort0(self):
        raise NotImplementedError('code this...')
        #msg = 'OES element_type=%-3s -> %-6s\n' %(self.element_type,self.get_element_type(self.element_type))
        msg = ''
        
        if self.element_type in [1, 3, 10]:  # crod/ctube/conrod
            self.makeOES_Object(self.rodStress, rodStressObject,
                                self.rodStrain, rodStrainObject)
            self.basicElement()
        elif self.element_type == 2:   # cbeam
            #print "    found cbeam_2"
            self.data_code['element_name'] = 'CBEAM'
            self.makeOES_Object(self.beamStress, beamStressObject,
                                self.beamStrain, beamStrainObject)
            self.CBEAM_2()

        elif self.element_type in [4]:  # cshear
            #print "    found crod_1"
            self.data_code['element_name'] = 'CSHEAR'
            self.makeOES_Object(self.shearStress, shearStressObject,
                                self.shearStrain, shearStrainObject)
            self.basicElement()
        elif self.element_type in [11, 12, 13]:   # celas1/celas2/celas3
            self.makeOES_Object(self.celasStress, celasStressObject,
                                self.celasStrain, celasStrainObject)
            self.basicElement()
        elif self.element_type == 34:   # cbar
            #print "    found cbar_34"
            self.data_code['element_name'] = 'CBAR'
            self.makeOES_Object(self.barStress, barStressObject,
                                self.barStrain, barStrainObject)
            self.CBAR_34()

        elif self.element_type == 33:  # cquad4_33
            #print "    found cquad_33"
            self.data_code['element_name'] = 'CQUAD4'
            self.makeOES_Object(self.plateStress, plateStressObject,
                                self.plateStrain, plateStrainObject)
            self.CQUAD4_33()
        elif self.element_type == 74:  # ctria
            #print "    found ctria_74"
            self.data_code['element_name'] = 'CTRIA3'
            self.makeOES_Object(self.plateStress, plateStressObject,
                                self.plateStrain, plateStrainObject)
            self.CTRIA3_74()  # ctria3
        elif self.element_type in [64, 144, 70, 75]:  # cquad8/cquad4/ctriar/ctria6
            #print "    found cquad_144"
            if     self.element_type == 64:
                self.data_code['element_name'] = 'CQUAD8'
            elif   self.element_type == 144:
                self.data_code['element_name'] = 'CQUAD4'
            elif   self.element_type == 70:
                self.data_code['element_name'] = 'CTRIAR'
            elif   self.element_type == 75:
                self.data_code['element_name'] = 'CTRIA6'
            else:
                raise NotImplementedError(self.element_type)
            self.makeOES_Object(self.plateStress, plateStressObject,
                                self.plateStrain, plateStrainObject)
            self.CQUAD4_144()

        elif self.element_type in [39, 67, 68]:   # ctetra/chexa/cpenta
            #print "    found ctetra_39 / hexa_67 / cpenta_68"
            self.makeOES_Object(self.solidStress, solidStressObject,
                                self.solidStrain, solidStrainObject)
            self.CSOLID_67()

        elif self.element_type in [85]:   # ctetra/chexa/cpenta (91,93)
            #print "    found ctetra_85 / hexa_93 / cpenta_91"
            self.makeOES_Object(self.solidStress, solidStressObject,
                                self.solidStrain, solidStrainObject)
            self.CSOLID_85()

        elif self.element_type in [89, 92]:   # RODNL
            #print "    found RODNL_89"
            self.makeOES_Object(self.rodStress, nonlinearRodObject,
                                self.rodStrain, nonlinearRodObject)
            self.RODNL_89_92()

        #elif self.element_type in [91]: # CPENTANL
        #    #print "hexa_93"
        #    self.CPENTANL_91()
        #elif self.element_type in [93]: # CHEXANL
        #    #print "hexa_93"
        #    self.CHEXANL_93()

        elif self.element_type in [95, 96, 97, 98]:  # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            self.eid2 = None  # stores the previous elementID
            self.makeOES_Object(
                self.compositePlateStress, compositePlateStressObject,
                self.compositePlateStrain, compositePlateStrainObject)
            self.CQUAD4_95()
            del self.eid2

        elif self.element_type in [94]:  # CBEAM (nonlinear)
            print "    found a 94!"
            #self.eid2 = None # stores the previous elementID
            self.makeOES_Object(self.beamStress, beamStressObject,
                                self.beamStrain, beamStrainObject)
            self.CBEAM_94()
            #sys.exit('stoping at end of CBEAM_94')
            #del self.eid2

        elif self.element_type in [139]:   # QUAD4FD
            #print "    found QUAD4FD_139"
            self.makeOES_Object(self.plateStress, nonlinearQuadObject,
                                self.plateStrain, nonlinearQuadObject)
            self.QUAD4FD_139()

        else:
            #self.print_block(self.data[0:100])
            self.skipOES_Element()
            msg = 'OES format1_sort0 element_type=%-3s -> %s is not supported - fname=%s\n' % (self.element_type, self.get_element_type(self.element_type), self.op2FileName)
            self.log.debug(msg)
            #msg = 'OES format1_sort0 element_type=%-3s -> %s is not supported' %(self.element_type,self.get_element_type(self.element_type))
            #raise RuntimeError(msg)
            self.skippedCardsFile.write(msg)
