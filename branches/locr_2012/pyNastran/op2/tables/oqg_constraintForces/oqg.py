import sys
from struct import unpack

from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import(
    SPCForcesObject, ComplexSPCForcesObject)
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import(
    MPCForcesObject, ComplexMPCForcesObject)


class OQG(object):
    """Table of spc/mpc forces/momenets"""

    def readTable_OQG(self):
        table3 = self.readTable_OQG_3
        table4Data = self.readOQG_Data
        self.readResultsTable(table3, table4Data)
        self.deleteAttributes_OQG()

    def deleteAttributes_OQG(self):
        #print self.obj
        params = ['lsdvm', 'mode', 'eigr', 'modeCycle', 'freq', 'dt', 'lftsfq', 'thermal', 'rCode', 'fCode', 'numWide', 'acousticFlag', 'thermal']
        self.deleteAttributes(params)

    def readTable_OQG_3(self, iTable):  # iTable=-3
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
        self.addDataParameter( data, 'randomCode', 'i', 8, False)
        ## format code        self.addDataParameter( data, 'formatCode', 'i', 9, False)
        ## number of words per entry in record
        ## @note is this needed for this table ???        self.addDataParameter(data, 'numWide', 'i', 10, False)
        ## acoustic pressure flag
        self.addDataParameter(data, 'acousticFlag', 'f', 13, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.addDataParameter(data, 'thermal', 'i', 23, False)

        if not self.isSort1():
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode == 1:   # statics / displacement / heat flux
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5, False)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysisCode == 2:  # real eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle @todo confused on the type - F1???
            self.addDataParameter(data, 'modeCycle', 'f', 7, False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'modeCycle'])
        #elif self.analysisCode==3: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
            #self.dataCode['lsdvmn'] = self.lsdvmn
        #elif self.analysisCode==4: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.analysisCode == 5:   # frequency
            ## frequency
            self.addDataParameter(data, 'freq', 'f', 5)
            self.applyDataCodeValue('dataNames', ['freq'])
        elif self.analysisCode == 6:  # transient
            ## time step
            self.addDataParameter(data, 'dt', 'f', 5)
            self.applyDataCodeValue('dataNames', ['dt'])
        elif self.analysisCode == 7:  # pre-buckling
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
        elif self.analysisCode == 8:  # post-buckling
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            self.applyDataCodeValue('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysisCode == 9:  # complex eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.addDataParameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.addDataParameter(data, 'eigi', 'f', 7, False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysisCode == 10:  # nonlinear statics
            ## load step
            self.addDataParameter(data, 'lftsfq', 'f', 5)
            self.applyDataCodeValue('dataNames', ['lftsfq'])
        elif self.analysisCode == 11:  # old geometric nonlinear statics
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
        elif self.analysisCode == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.addDataParameter(data, 'lsdvmn', 'i', 5)
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
        else:
            msg = 'invalid analysisCode...analysisCode=%s' % (self.analysisCode)
            raise RuntimeError(msg)
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)

        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOQG_Data(self):
        #tfsCode = [self.tableCode,self.formatCode,self.sortCode]

        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        assert self.thermal in [0, 1, 8], self.codeInformation()

        if   self.tableCode == 3:   # SPC Forces
            assert self.tableName in ['OQG1', 'OQGV1', 'OQP1'], 'tableName=%s tableCode=%s' % (self.tableName, self.tableCode)
            self.readOQG_Data_table3()
        elif self.tableCode == 39:  # MPC Forces
            assert self.tableName in ['OQMG1'], 'tableName=%s tableCode=%s' % (
                self.tableName, self.tableCode)
            self.readOQG_Data_table3()
        else:
            self.NotImplementedOrSkip('bad analysis/table/format/sortCode=%s' %
                                      (self.atfsCode))
        #print self.obj

    def readOQG_Data_table3(self):  # SPC Forces
        #isSort1 = self.isSort1()
        #print(self.codeInformation())
        magPhase = self.isMagnitudePhase()
        if magPhase or self.numWide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'spcForces'
                self.createTransientObject(self.spcForces,
                                           ComplexSPCForcesObject)
                self.handleResultsBuffer3(self.OUG_ComplexTable, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.numWide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'spcForces'
                self.createTransientObject(
                    self.spcForces, SPCForcesObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip('only numWide=8 or 14 is allowed  numWide=%s' % (self.numWide))

        #if self.thermal not in [0,1]:
            #print self.obj
            #raise RuntimeError('check the printout for thermal...')

    def readOQG_Data_table39(self):  # MPC Forces
        #isSort1 = self.isSort1()
        if self.numWide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'mpcForces'
                self.createTransientObject(
                    self.mpcForces, MPCForcesObject)
                self.handleResultsBuffer3(self.OUG_RealTable, resultName)
            else:
                self.NotImplementedOrSkip()
        elif self.numWide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'mpcForces'
                self.createTransientObject(self.mpcForces,
                                           ComplexMPCForcesObject)
                self.handleResultsBuffer3(self.OUG_ComplexTable, resultName)
            else:
                self.NotImplementedOrSkip()
        else:
            self.NotImplementedOrSkip('only numWide=8 or 14 is allowed  numWide=%s' % (self.numWide))
