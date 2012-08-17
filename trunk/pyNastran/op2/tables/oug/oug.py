# http://www.cadfamily.com/online-help/I-DEAS/SDRCHelp/LANG/English/slv_ug/NAS_results_imported.htm
#import sys
from struct import unpack

from pyNastran.op2.tables.oug.oug_displacements import (
    DisplacementObject,              # tableCode=1     formatCode=1 sortCode=0
    ComplexDisplacementObject)       # analysisCode=5  formatCode=3 sortCode=1

# tableCode=10 formatCode=1 sortCode=0
from pyNastran.op2.tables.oug.oug_velocities import (
    VelocityObject, ComplexVelocityObject)

# tableCode=11 formatCode=1 sortCode=0
from pyNastran.op2.tables.oug.oug_accelerations import (
    AccelerationObject, ComplexAccelerationObject)

# tableCode=1 formatCode=1 sortCode=0
from pyNastran.op2.tables.oug.oug_temperatures import (
    TemperatureObject)

from pyNastran.op2.tables.oug.oug_eigenvectors import (
     EigenVectorObject,                     # analysisCode=2, sortCode=0 formatCode   tableCode=7
     ComplexEigenVectorObject,              # analysisCode=5, sortCode=1 formatCode=1 tableCode=7
    #RealEigenVectorObject,                 # analysisCode=9, sortCode=1 formatCode=1 tableCode=7
)
from pyNastran.op2.tables.opg_appliedLoads.opg_loadVector import ThermalVelocityVectorObject
from pyNastran.op2.op2_helper import polarToRealImag


class OUG(object):
    """Table of displacements/velocities/acceleration/heat flux/temperature"""

    def readTable_OUG(self):
        #self.tableName = 'OUG'
        table3 = self.readTable_OUG_3
        table4Data = self.readOUG_Data
        self.readResultsTable(table3, table4Data)
        self.deleteAttributes_OUG()

    def deleteAttributes_OUG(self):
        params = ['lsdvm', 'mode', 'eigr', 'modeCycle', 'freq', 'dt', 'lftsfq',
                  'thermal', 'randomCode', 'fCode', 'numWide', 'acousticFlag']
        self.deleteAttributes(params)

    def readTable_OUG_3(self, iTable):  # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' % (str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i', data)
        data = self.getData(4 * 50)
        #print self.printBlock(data)

        (three) = self.parseApproachCode(data)

        ## random code
        self.addDataParameter(data, 'randomCode', 'i', 8, False)
        ## format code
        self.addDataParameter(data, 'formatCode', 'i', 9, False)
        ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data, 'numWide', 'i', 10, False)
        ## acoustic pressure flag
        self.addDataParameter(data, 'acousticFlag', 'f', 13, False)
        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.addDataParameter(data, 'thermal', 'i', 23, False)
        self.isFlipped = False
        if self.isSort1():
            ## assuming tCode=1
            if self.analysisCode == 1:   # statics / displacement / heat flux
                self.addDataParameter(data, 'lsdvmn', 'i',
                                      5, False)  # load set number
                self.applyDataCodeValue('dataNames', ['lsdvmn'])
                self.setNullNonlinearFactor()
            elif self.analysisCode == 2:  # real eigenvalues
                self.addDataParameter(data, 'mode',
                                      'i', 5)  # mode number
                self.addDataParameter(data, 'eigr',
                                      'f', 6, False)  # real eigenvalue
                self.addDataParameter(data, 'modeCycle', 'i', 7, False)  # mode or cycle @todo confused on the type - F1???
                self.applyDataCodeValue('dataNames', [
                    'mode', 'eigr', 'modeCycle'])
            #elif self.analysisCode==3: # differential stiffness
                #self.lsdvmn = self.getValues(data,'i',5) ## load set number
                #self.dataCode['lsdvmn'] = self.lsdvmn
            #elif self.analysisCode==4: # differential stiffness
                #self.lsdvmn = self.getValues(data,'i',5) ## load set number
            elif self.analysisCode == 5:   # frequency
                self.addDataParameter(data, 'freq', 'f', 5)  # frequency
                self.applyDataCodeValue('dataNames', ['freq'])
            elif self.analysisCode == 6:  # transient
                self.addDataParameter(data, 'dt', 'f', 5)  # time step
                self.applyDataCodeValue('dataNames', ['dt'])
            elif self.analysisCode == 7:  # pre-buckling
                self.addDataParameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.applyDataCodeValue('dataNames', ['lsdvmn'])
            elif self.analysisCode == 8:  # post-buckling
                self.addDataParameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.addDataParameter(data, 'eigr', 'f',
                                      6, False)  # real eigenvalue
                self.applyDataCodeValue('dataNames', ['lsdvmn', 'eigr'])
            elif self.analysisCode == 9:  # complex eigenvalues
                self.addDataParameter(data, 'mode', 'i', 5)  # mode number
                self.addDataParameter(data, 'eigr', 'f', 6,
                                      False)  # real eigenvalue
                self.addDataParameter(data, 'eigi', 'f', 7,
                                      False)  # imaginary eigenvalue
                self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'eigi'])
            elif self.analysisCode == 10:  # nonlinear statics
                self.addDataParameter(data, 'lftsfq', 'f', 5)  # load step
                self.applyDataCodeValue('dataNames', ['lftsfq'])
            elif self.analysisCode == 11:  # old geometric nonlinear statics
                self.addDataParameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.applyDataCodeValue('dataNames', ['lsdvmn'])
            elif self.analysisCode == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
                self.addDataParameter(data, 'lsdvmn', 'i',
                                      5)  # load set number
                self.applyDataCodeValue('dataNames', ['lsdvmn'])
            else:
                raise RuntimeError('invalid analysisCode...analysisCode=%s' %
                                   (self.analysisCode))
        else:  # sort2

            eidDevice = self.getValues(data, 'i', 5)
            floatVal = self.getValues(data, 'f', 5)
            #eid = (eidDevice-self.deviceCode)//10
            #print("EID = %s" %(eidDevice))
            #print("floatVal = %s" %(floatVal))

            if self.tableName == 'OUGRMS2' and self.analysisCode == 1:
                self.addDataParameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # frequency
                self.applyDataCodeValue('dataNames', ['nodeID'])

            #self.isRegular = False
            elif self.analysisCode in [1, 5]:  # 5 # freq
                self.addDataParameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # frequency
                self.applyDataCodeValue('dataNames', ['nodeID'])
                #print("nodeID = %s" %(self.nodeID))
                #sys.exit(self.nodeID)
            elif self.analysisCode == 6:  # transient dt
                self.addDataParameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # time step
                self.applyDataCodeValue('dataNames', ['nodeID'])
            elif self.analysisCode == 10:  # freq/time step fqts
                self.addDataParameter(data, 'nodeID', 'i', 5, fixDeviceCode=True)  # frequency / time step
                self.applyDataCodeValue('dataNames', ['nodeID'])
            else:
                self.isRegular = True
                self.addDataParameter(data, 'nodeID', 'i',
                                      5, fixDeviceCode=True)  # node ID
                self.applyDataCodeValue('dataNames', ['nodeID'])
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)

        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        if not self.isSort1():
            raise NotImplementedError('sort2...')

        #self.printBlock(data)
        self.readTitle()

    def getOUG_FormatStart(self):
        """
        Returns an i or an f depending on if it's SORT2 or not.
        Also returns an extraction function that is called on the first argument
        """
        isSort1 = self.isSort1()
        if self.tableName == 'OUGRMS2' and self.analysisCode == 1:
            format1 = 'i'  # SORT2
            extract = self.extractSort2

        elif isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            #print "SORT1"
            format1 = 'i'  # SORT1
            extract = self.extractSort1
            #if self.analysisCode in [5]:
                #extract==self.extractSort2
        else:  # values from IDENT   #@todo test this...
            #print "SORT2"
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            if self.analysisCode in [2, 3, 4, 6, 7, 8, 11]:
                format1 = 'f'  # SORT2
                extract = self.extractSort2
            elif self.analysisCode in [5]:
                format1 = 'f'
                extract = self.extractSort2
            elif self.analysisCode in [1, 9, 10, 12]:
                format1 = 'f'  # SORT1
                extract = self.extractSort2
            else:
                raise KeyError('invalid analysisCode...analysisCode=%s' %
                               (self.analysisCode))
            ###
            #eid = self.nonlinearFactor
        return (format1, extract)

    def readOUG_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        #tfsCode = [self.tableCode,self.formatCode,self.sortCode]

        #print self.dataCode
        #print "tfsCode=%s" %(tfsCode)

        if self.tableCode == 1 and self.tableName in ['OUGV1', 'OUPV1']:    # displacement
            if self.tableName == 'OUGV1':
                assert self.tableName in ['OUGV1'], 'tableName=%s tableCode=%s\n%s' % (self.tableName, self.tableCode, self.codeInformation())
                self.readOUG_Data_table1()
            else:  # 'OUPV1'
                self.NotImplementedOrSkip('bad approach/table/format/sortCode=%s on %s-OUG table' % (self.atfsCode, self.tableName))
        elif self.tableCode == 1 and self.tableName in ['OUGATO2', 'OUGCRM2', 'OUGPSD2', 'OUGRMS2', 'OUGNO2', ]:    # displacement
            #assert self.tableName in ['OUGATO2','OUGCRM2','OUGPSD2','OUGRMS2','OUGNO2',],'tableName=%s tableCode=%s\n%s' %(self.tableName,self.tableCode,self.codeInformation())
            self.readOUG_Data_table1()
        elif self.tableCode == 7:  # modes
            if self.tableName == 'OUGV1':
                assert self.tableName in ['OUGV1'], 'tableName=%s tableCode=%s\n%s' % (self.tableName, self.tableCode, self.codeInformation())
                self.readOUG_Data_table7()
            else:
                self.NotImplementedOrSkip('bad approach/table/format/sortCode=%s on %s-OUG table' % (self.atfsCode, self.tableName))
        elif self.tableCode == 10:  # velocity
            if self.tableName == 'OUGV1':
                assert self.tableName in ['OUGV1'], 'tableName=%s tableCode=%s\n%s' % (self.tableName, self.tableCode, self.codeInformation())
                self.readOUG_Data_table10()
            else:
                self.NotImplementedOrSkip('bad approach/table/format/sortCode=%s on %s-OUG table' % (self.atfsCode, self.tableName))
        elif self.tableCode == 11:  # Acceleration vector
            if self.tableName == 'OUGV1':
                assert self.tableName in ['OUGV1'], 'tableName=%s tableCode=%s\n%s' % (self.tableName, self.tableCode, self.codeInformation())
                self.readOUG_Data_table11()
            else:
                self.NotImplementedOrSkip('bad approach/table/format/sortCode=%s on %s-OUG table' % (self.atfsCode, self.tableName))
        else:
            #self.log.debug('skipping approach/table/format/sortCode=%s on %s-OUG table' %(self.atfsCode,self.tableName))
            self.NotImplementedOrSkip('bad approach/table/format/sortCode=%s on %s-OUG table' % (self.atfsCode, self.tableName))
        ###
        #print self.obj

    def readThermal4(self):  # used on self.thermal in [2,4,8]:
        #print self.codeInformation()
        #print self.printBlock(self.data)
        n = 0
        nEntries = len(self.data) // 32
        for i in xrange(nEntries):
            eData = self.data[n:n + 32]
            out = unpack('2i6f', eData)
            #nid = (out[0]-self.deviceCode)//10    ## @todo fix the deviceCode

            #print out
            n += 32
            #print "nid = ",nid
        #sys.exit('thermal4...')

    def readOUG_Data_table1(self):  # displacement / temperature OUGV1, OUPV1
        """
        OUGV1   - global coordinate system in sort 1
        OUPV1   - scaled response spectra in sort 1
        OUGPSD2 - PSD in sort 2
        OUGATO2 - auto-correlated in sort 2
        """
        isSkip = False
        if self.numWide == 8:  # real/random
            if self.thermal == 0:
                #print self.dataCode
                if self.tableName in ['OUGV1']:
                    resultName = 'displacements'
                    self.createTransientObject(self.displacements,
                                               DisplacementObject)
                elif self.tableName in ['OUGATO2']:
                    resultName = 'displacementsATO'
                    self.createTransientObject(self.displacementsATO,
                                               DisplacementObject)
                elif self.tableName in ['OUGCRM2']:
                    resultName = 'displacementsCRM'
                    self.createTransientObject(self.displacementsCRM,
                                               DisplacementObject)
                elif self.tableName in ['OUGPSD2']:
                    resultName = 'displacementsPSD'
                    self.createTransientObject(self.displacementsPSD,
                                               DisplacementObject)
                elif self.tableName in ['OUGRMS2']:
                    resultName = 'displacementsRMS'
                    self.createTransientObject(self.displacementsRMS,
                                               DisplacementObject)
                elif self.tableName in ['OUGNO2']:
                    resultName = 'displacementsNO'
                    self.createTransientObject(self.displacementsNO,
                                               DisplacementObject)
                else:
                    isSkip = True
                    self.NotImplementedOrSkip('***table=%s***\n%s' % (self.tableName, self.codeInformation()))
                if not isSkip:
                    self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            elif self.thermal == 1:
                resultName = 'temperatures'
                self.createTransientObject(
                    self.temperatures, TemperatureObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            #elif self.thermal == 8:
                #resultName = 'scaledDisplacements'
                #self.createTransientObject(self.scaledDisplacements,displacementObject)
                #self.handleResultsBuffer3(self.OUG_RealTable,resultName)
            else:
                self.NotImplementedOrSkip('***thermal=%s***\n%s' %
                                          (self.thermal, self.codeInformation()))
        elif self.numWide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'displacements'
                self.createTransientObject(self.displacements,
                                           ComplexDisplacementObject)
                self.handleResultsBuffer3(self.OUG_ComplexTable, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip('only numWide=8 or 14 is allowed  numWide=%s' % (self.numWide))
        ###

    def readOUG_Data_table7(self):  # eigenvector
        #isSort1 = self.isSort1()
        if self.numWide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'eigenvectors'
                self.createTransientObject(self.eigenvectors,
                                           EigenVectorObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.numWide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'eigenvectors'
                self.createTransientObject(self.eigenvectors,
                                           ComplexEigenVectorObject)
                self.handleResultsBuffer3(self.OUG_ComplexTable, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip('only numWide=8 or 14 is allowed  numWide=%s' % (self.numWide))

    def readOUG_Data_table10(self):  # velocity
        if self.numWide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'velocities'
                self.createTransientObject(
                    self.velocities, VelocityObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            elif self.thermal == 1:
                resultName = 'velocities'
                self.createTransientObject(self.velocities,
                                           ThermalVelocityVectorObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.numWide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'velocities'
                self.createTransientObject(self.velocities,
                                           ComplexVelocityObject)
                self.handleResultsBuffer3(self.OUG_ComplexTable, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip('only numWide=8 or 14 is allowed  numWide=%s' % (self.numWide))

    def readOUG_Data_table11(self):  # acceleration
        if self.numWide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'accelerations'
                self.createTransientObject(self.accelerations,
                                           AccelerationObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.numWide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'accelerations'
                self.createTransientObject(self.accelerations,
                                           ComplexAccelerationObject)
                self.handleResultsBuffer3(self.OUG_ComplexTable, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip('only numWide=8 or 14 is allowed  numWide=%s' % (self.numWide))

    def OUG_RealTable(self):
        dt = self.nonlinearFactor
        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i6f'

        #print "len(data) = ",len(self.data)
        while len(self.data) >= 32:  # 8*4
            eData = self.data[0:32]
            self.data = self.data[32:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, gridType, tx, ty, tz, rx, ry, rz) = out
            eid2 = extract(eid, dt)
            #print "eType=%s" %(eType)

            dataIn = [eid2, gridType, tx, ty, tz, rx, ry, rz]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #print "%s" %(self.tableName),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)

    def OUG_ComplexTable(self):
        dt = self.nonlinearFactor

        (format1, extract) = self.getOUG_FormatStart()
        format1 += 'i12f'
        #print "format1 = ",format1
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data) >= 56:  # 14*4
            eData = self.data[0:56]
            self.data = self.data[56:]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid, gridType, txr, tyr, tzr, rxr, ryr, rzr,
             txi, tyi, tzi, rxi, ryi, rzi) = out

            if isMagnitudePhase:
                tx = polarToRealImag(txr, txi)
                rx = polarToRealImag(rxr, rxi)
                ty = polarToRealImag(tyr, tyi)
                ry = polarToRealImag(ryr, ryi)
                tz = polarToRealImag(tzr, tzi)
                rz = polarToRealImag(rzr, rzi)
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
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt, dataIn)
            #print "len(data) = ",len(self.data)

