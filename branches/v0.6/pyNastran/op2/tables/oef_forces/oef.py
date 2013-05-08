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

from pyNastran import isRelease
from .realForces import RealForces
from .complexForces import ComplexForces


from .oef_forceObjects import (RealRodForce, RealCBeamForce, RealCShearForce,
                               RealSpringForce, RealDamperForce, RealViscForce,
                               RealPlateForce, RealConeAxForce, RealPlate2Force,
                               RealCBar100Force, RealCGapForce, RealBendForce,
                               RealPentaPressureForce, RealCBushForce,
                               RealForce_VU_2D, RealCBarForce, RealForce_VU)
from .oef_complexForceObjects import (ComplexRodForce, ComplexCBeamForce,
                                      ComplexCShearForce, ComplexSpringForce,
                                      ComplexDamperForce, ComplexViscForce,
                                      ComplexPlateForce, ComplexPlate2Force,
                                      ComplexBendForce,
                                      ComplexPentaPressureForce,
                                      ComplexCBushForce, ComplexForce_VU_2D,
                                      ComplexCBarForce, ComplexForce_VU)
from .thermal_elements import ThermalElements


class OEF(ThermalElements, RealForces, ComplexForces):
    """Table of element forces"""
    def readTable_OEF(self):
        table3 = self.readTable_OEF_3
        table4Data = self.readOEF_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OEF()

    def _delete_attributes_OEF(self):
        params = ['element_type', 'dLoadID', 'loadID', 'obj', 'markerStart', 'oCode',
                  'eigr', 'eigi', 'eign', 'mode', 'freq', 'time', 'thermal', ]
        self._delete_attributes(params)
        #print self.obj

    def readTable_OEF_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        data = self.get_data(4 * 50)
        #self.print_block(data)

        aCode = self.get_block_int_entry(data, 1)

        self.parse_approach_code(data)
        #: element type
        self.add_data_parameter( data, 'element_type', 'i', 3, False)  

        # dynamic load set ID/random code
        #self.add_data_parameter(data, 'dLoadID', 'i', 8, False)

        #: format code
        self.add_data_parameter( data, 'format_code', 'i', 9, False)
        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        #: undefined in DMAP...
        self.add_data_parameter(data, 'oCode', 'i', 11, False)
        #: thermal flag; 1 for heat ransfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        #print "dLoadID(8)=%s format_code(9)=%s numwde(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.format_code,self.num_wide,self.oCode,self.thermal)
        #print "thermal(23)=%s element_type(3)=%s" %(self.thermal,self.element_type)

        #if not self.is_sort1():
            #raise NotImplementedError('sort2...')

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            self.add_data_parameter(data, 'loadID', 'i', 5,
                                  False)  # load set ID number
            #self.apply_data_code_value('dataNames',['lsdvmn'])
            self.apply_data_code_value('dataNames', ['loadID'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            #: eigenvalue
            self.add_data_parameter(data, 'eign', 'f', 6, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
            #print "mode(5)=%s eigr(6)=%s mode_cycle(7)=%s" %(self.mode,self.eigr,self.mode_cycle)
        elif self.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            self.add_data_parameter(data, 'loadID', 'i', 5)
            self.apply_data_code_value('dataNames', ['loadID'])
        elif self.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            self.add_data_parameter(data, 'loadID', 'i', 5)
            self.apply_data_code_value('dataNames', ['loadID'])
        elif self.analysis_code == 5:   # frequency
            self.add_data_parameter(data, 'freq', 'f', 5)  # frequency
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            self.add_data_parameter(data, 'time', 'f', 5)  # time step
            self.apply_data_code_value('dataNames', ['time'])
        elif self.analysis_code == 7:  # pre-buckling
            #: load set ID number
            self.add_data_parameter(data, 'loadID', 'i', 5)
            #self.apply_data_code_value('dataNames',['lsdvmn'])
            self.apply_data_code_value('dataNames', ['loadID'])
        elif self.analysis_code == 8:  # post-buckling
            #: load set ID number
            self.add_data_parameter(data, 'loadID', 'i', 5)
            #: real eigenvalue
            self.add_data_parameter( data, 'eigr', 'f', 6, False)
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
            #print "loadID(5)=%s  eigr(6)=%s" %(self.loadID,self.eigr)
        elif self.analysis_code == 9:  # complex eigenvalues
            #: mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            #: real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            #: imaginary eigenvalue
            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            #: load step
            self.add_data_parameter(data, 'load_step', 'f', 5)
            self.apply_data_code_value('dataNames', ['load_step'])
        elif self.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            self.add_data_parameter(data, 'loadID', 'i', 5)
            self.apply_data_code_value('dataNames', ['loadID'])
            #print "loadID(5)=%s" %(self.loadID)
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % (str(self.analysis_code) + '\n' + self.code_information()))

        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.loadID = self.get_values(data,'i',5) ## load set ID number

        if not self.is_sort1():
            if isRelease:
                self.not_implemented_or_skip('skipping OES SORT2')
            else:
                raise NotImplementedError('SORT2')

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)
        #print self.code_information()

        #self.print_block(data)
        #print '-'*80
        self.read_title()

    def OEF_ForceCode(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        realMapper = {
            1: 3,    # CROD
            2: 1 + (10 - 1) * 11,  # CBEAM
            3: 3,    # CTUBE
            4: 17,   # CSHEAR
            10: 3,    # CONROD
            11: 2,    # CELAS1
            12: 2,    # CELAS2
            13: 2,    # CELAS3
            14: 2,    # CELAS4

            20: 2,    # CDAMP1
            21: 2,    # CDAMP2
            22: 2,    # CDAMP3
            23: 2,    # CDAMP4
            24: 3,    # CVISC
            33: 9,    # CQUAD4
            34: 9,    # CBAR
            35: 7,    # CCONEAX
            38: 9,    # CGAP
            40: 8,    # CBUSH1D ???
            64: 2 + (11 - 2) * 5,  # CQUAD8
            69: 1 + (8 - 1) * 2,  # CBEND
            70: 2 + (11 - 2) * 4,  # CTRIAR
            74: 9,    # CTRIA3
            75: 2 + (11 - 2) * 4,  # CTRIA6


            #76:  16,   # Acoustic Velocity/Pressure CHEXA ???
            76: None,  # dummy so it doesnt go into the real results
            77: 10,   # Acoustic Velocity/Pressure CPENTA
            78: 10,   # Acoustic Velocity/Pressure CTETRA

            82: 2 + (11 - 2) * 5,  # CQUADR
            95: 9,    # composite CQUAD4 ???
            96: 9,    # composite CQUAD8 ???
            97: 9,    # composite CTRIA3 ???
            98: 9,    # composite CTRIA6 ???
            100: 8,    # BARS
            102: 7,    # CBUSH
            144: 2 + (11 - 2) * 5,  # bilinear CQUAD4
            189: 6 + (19 - 6) * 4,  # VUQUAD
            190: 6 + (19 - 6) * 3,  # VUTRIA
            191: 4 + (12 - 4) * 2,  # VUBEAM
            200: 9,    # CWELD
            232: 9,    # composite CQUADR ???
            233: 9,    # composite TRIAR ???
            235: 9,    # punch CQUADR...num_wide in DMAP is wrong...left out first entry...
            236: 8,    # punch CTRIAR
        }
        imagMapper = {
            1: 5,    # CROD
            2: 1 + (17 - 1) * 11,  # CBEAM
            3: 5,    # CTUBE
            4: 33,   # CSHEAR
            10: 5,    # CONROD

            11: 3,    # CELAS1
            12: 3,    # CELAS2
            13: 3,    # CELAS3
            14: 3,    # CELAS4

            20: 3,    # CDAMP1
            21: 3,    # CDAMP2
            22: 3,    # CDAMP3
            23: 3,    # CDAMP4
            24: 5,    # CVISC
            33: 17,   # CQUAD4
            34: 17,   # CBAR
            35: 7,    # CCONEAX # needed to not crash the code...
            38: 9,    # CGAP
            40: 8,    # CBUSH1D ???
            64: 2 + (19 - 2) * 5,  # CQUAD8
            69: 1 + (14 - 1) * 2,  # CBEND
            70: 2 + (19 - 2) * 4,  # CTRIAR
            74: 17,   # CTRIA3
            75: 2 + (19 - 2) * 4,  # CTRIA6

            76: 16,   # Acoustic Velocity/Pressure CHEXA_PR
            77: 16,   # Acoustic Velocity/Pressure CPENTA_PR
            78: 16,   # Acoustic Velocity/Pressure CTETRA_PR

            82: 2 + (19 - 2) * 5,  # CQUADR
            95: 9,    # composite CQUAD4 ???
            96: 9,    # composite CQUAD8 ???
            97: 9,    # composite CTRIA3 ???
            98: 9,    # composite CTRIA6 ???
            100: 14,   # BARS
            102: 13,   # CBUSH
            144: 2 + (19 - 2) * 5,  # bilinear CQUAD4
            189: 6 + (31 - 6) * 4,  # VUQUAD
            190: 6 + (31 - 6) * 3,  # VUTRIA
            191: 4 + (18 - 4) * 2,  # VUBEAM
            200: 17,   # CWELD
            232: 9,    # composite CQUADR ???
            233: 9,    # composite TRIAR ???
            235: 17,   # punch CQUADR...num_wide in DMAP is wrong...left out first entry...
            236: 16,   # punch CTRIAR
        }

        Real = realMapper[self.element_type]
        Imag = imagMapper[self.element_type]
        return (Real, Imag)

    def readOEF_Data(self):
        """
        OEF1X -
        DOEF1 -
        """
        if self.table_code == 4 and self.table_name in ['OEF1X', 'DOEF1']:  # Forces/Heat Flux
            assert self.table_name in ['OEF1X', 'DOEF1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self.readOEF_Data_table4()
        #elif self.table_name in ['OEFATO2','OEFCRM2','OEFPSD2','OEFRMS2','OEFNO2',]:
            #self.skipOES_Element() # skipping entire table
        else:
            self.not_implemented_or_skip()

    def readOEF_Data_table4(self):  # Forces/Heat Flux
        if self.thermal in [0, 8]:
            self.readOEF_Forces()
        elif self.thermal == 1:
            self.readOEF_Thermal()
        else:
            self.not_implemented_or_skip('thermal=%s' % (self.thermal))

        #if self.thermal==8:
        #    self.not_implemented_or_skip('thermal=%s' %(self.thermal))

    def readOEF_Forces(self):
        try:
            (numWideReal, numWideImag) = self.OEF_ForceCode()
        except KeyError:
            self.not_implemented_or_skip()

        #print(self.code_information())

        if self.element_type in [1, 3, 10]:  # CROD,CTUBE,CONROD
            resultName = 'rodForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.rodForces, RealRodForce)
                self.handle_results_buffer(self.OEF_Rod, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.rodForces, ComplexRodForce)
                self.handle_results_buffer(self.OEF_Rod_alt, resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [2]:  # CBEAM
            resultName = 'beamForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.beamForces, RealCBeamForce)
                self.handle_results_buffer(self.OEF_Beam, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.beamForces, ComplexCBeamForce)
                self.handle_results_buffer(self.OEF_Beam_alt, resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [4]:  # CSHEAR
            resultName = 'shearForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.shearForces, RealCShearForce)
                self.handle_results_buffer(self.OEF_Shear, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(
                    self.shearForces, ComplexCShearForce)
                self.handle_results_buffer(self.OEF_Shear_alt, resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [11, 12, 13, 14]:  # CELAS1,CELAS2,CELAS3,CELAS4
            resultName = 'springForces'
            if self.num_wide == numWideReal:  # .. todo:: is this correct or is DMAP wrong (CELAS1)???
                self.create_transient_object(self.springForces, RealSpringForce)
                self.handle_results_buffer(self.OEF_Spring, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(
                    self.springForces, ComplexSpringForce)
                self.handle_results_buffer(self.OEF_Spring_alt, resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [20, 21, 22, 23]:  # CDAMP1,CDAMP2,CDAMP3,CDAMP4
            resultName = 'damperForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.damperForces, RealDamperForce)
                self.handle_results_buffer(self.OEF_Spring,
                                          resultName)  # same reader as for springs
            elif self.num_wide == numWideImag:
                self.create_transient_object(
                    self.damperForces, ComplexDamperForce)
                self.handle_results_buffer(self.OEF_Spring_alt, resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [24]:  # CVISC
            resultName = 'viscForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.viscForces, RealViscForce)
                self.handle_results_buffer(self.OEF_CVisc, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.viscForces, ComplexViscForce)
                self.handle_results_buffer(self.OEF_CVisc_alt, resultName)
            else:
                self.not_implemented_or_skip()

        elif self.element_type in [33, 74, 235]:  # CQUAD4,CTRIA3,CQUADR
            resultName = 'plateForces'
            if self.num_wide == numWideReal:
                #print(self.code_information())
                self.create_transient_object(self.plateForces, RealPlateForce)
                self.handle_results_buffer(self.OEF_Plate, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.plateForces, ComplexPlateForce)
                self.handle_results_buffer(self.OEF_Plate_alt, resultName)
            else:
                self.not_implemented_or_skip()

        #elif self.element_type in [235]: # CQUADR
            #if self.num_wide == numWideReal:
                #self.OEF_Plate()
            #elif self.num_wide == numWideImag:
                #self.OEF_Plate_alt()
            #else:
                #raise NotImplementedError(self.code_information())
        elif self.element_type in [35]:  # CCONEAX
            resultName = 'coneAxForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.coneAxForces, RealConeAxForce)
                self.handle_results_buffer(self.OEF_ConeAx, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8,CTRIAR,CTRIA6,CQUADR,CQUAD4-bilinear
            resultName = 'plateForces2'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.plateForces2, RealPlate2Force)
                self.handle_results_buffer(self.OEF_Plate2, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(
                    self.plateForces2, ComplexPlate2Force)
                self.handle_results_buffer(self.OEF_Plate2_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [34]:  # CBAR
            resultName = 'barForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.barForces, RealCBarForce)
                self.handle_results_buffer(self.OEF_CBar, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.barForces, ComplexCBarForce)
                self.handle_results_buffer(self.OEF_CBar_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [100]:  # CBAR
            resultName = 'barForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.bar100Forces, RealCBar100Force)
                self.handle_results_buffer(self.OEF_CBar100, resultName)
            elif self.num_wide == numWideImag:
                self.handle_results_buffer(self.OEF_CBar100_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [38]:  # CGAP
            resultName = 'gapForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.gapForces, RealCGapForce)
                self.handle_results_buffer(self.OEF_CGap, resultName)
            elif self.num_wide == numWideImag:
                self.handle_results_buffer(self.OEF_CGap_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [69]:  # CBEND
            resultName = 'bendForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.bendForces, RealBendForce)
                self.handle_results_buffer(self.OEF_Bend, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.bendForces, ComplexBendForce)
                self.handle_results_buffer(self.OEF_Bend_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [76, 77, 78]:  # CHEXA_PR,PENTA_PR,CTETRA_PR
            resultName = 'solidPressureForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.solidPressureForces,
                                           RealPentaPressureForce)
                self.handle_results_buffer(self.OEF_PentaPressure, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.solidPressureForces,
                                           ComplexPentaPressureForce)
                self.handle_results_buffer(
                    self.OEF_PentaPressure_alt, resultName)
            else:
                self.not_implemented_or_skip()
        #elif self.element_type in [95,96,97,98]: # composite CQUAD4,CQUAD8,CTRIA3,CTRIA6
        #    resultName = '???' # TODO why is there no object...
        #    if self.num_wide == numWideReal:
        #        self.handle_results_buffer(self.OEF_CompositePlate,resultName)
        #    elif self.num_wide == numWideImag:
        #        self.handle_results_buffer(self.OEF_CompositePlate_alt)
        #    else:
        #        self.not_implemented_or_skip()
        elif self.element_type in [102]:  # CBUSH
            resultName = 'bushForces'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.bushForces, RealCBushForce)
                self.handle_results_buffer(self.OEF_CBush, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.bushForces, ComplexCBushForce)
                self.handle_results_buffer(self.OEF_CBush_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [189, 190]:  # VUQUAD,VUTRIA
            resultName = 'force_VU_2D'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.force_VU_2D, RealForce_VU_2D)
                self.handle_results_buffer(self.OEF_Force_VUTRIA, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(
                    self.force_VU_2D, ComplexForce_VU_2D)
                self.handle_results_buffer(
                    self.OEF_Force_VUTRIA_alt, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.element_type in [191]:  # VUBEAM
            resultName = 'force_VU'
            if self.num_wide == numWideReal:
                self.create_transient_object(self.force_VU, RealForce_VU)
                self.handle_results_buffer(self.OEF_Force_VU, resultName)
            elif self.num_wide == numWideImag:
                self.create_transient_object(self.force_VU, ComplexForce_VU)
                self.handle_results_buffer(self.OEF_Force_VU_alt, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip()
