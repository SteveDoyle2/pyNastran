from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack

# pyNastran
from pyNastran.op2.op2_helper import polarToRealImag
from .ogf_Objects import gridPointForcesObject, complexGridPointForcesObject


class OGF(object):
    """Table of Grid Point Forces"""

    def readTable_OGF(self):
        table3 = self.readTable_OGF_3
        table4Data = self.readOGF_Data
        self.readResultsTable(table3, table4Data)
        self.deleteAttributes_OGF()

    def deleteAttributes_OGF(self):
        params = ['formatCode', 'appCode', 'numWide', 'value1', 'value2', ]
        self.deleteAttributes(params)

    #def addDataParameter(self,data,Name,Type,FieldNum,applyNonlinearFactor=True):
    #    #self.mode = self.getValues(data,'i',5) ## mode number
    #    value = self.getValues(data,Type,FieldNum)
    #    setattr(self,Name,value)
    #    self.dataCode[Name] = value
    #
    #    if applyNonlinearFactor:
    #        self.nonlinearFactor = value
    #        self.dataCode['nonlinearFactor'] = value
    #        self.dataCode['name'] = Name

    def applyDataCodeValue(self, Name, value):
        self.dataCode[Name] = value

    def readTable_OGF_3(self, iTable):  # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' % (str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack(b'i', data)
        data = self.getData(4 * 50)
        #print self.printBlock(data)

        (three) = self.parseApproachCode(data)
        ## format code
        self.addDataParameter(data,'formatCode', 'i', 9, False)
        ## approach code ???
        self.addDataParameter(data,'appCode', 'i', 9, False)
        ## number of words per entry in record
        ## @note is this needed for this table ???
        self.addDataParameter(data,'numWide', 'i', 10, False)
        ## Data Value 1
        self.addDataParameter(data,'value1', 'i', 11, False)
        ## Data Value 2
        self.addDataParameter(data,'value2', 'f', 12, False)
        
        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode == 1:   # statics / displacement / heat flux
            #self.extractDt = self.extractInt
            self.applyDataCodeValue('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysisCode == 2:  # real eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            self.applyDataCodeValue('dataNames', ['mode'])
        #elif self.analysisCode==3: # differential stiffness
            ## load set number
            #self.lsdvmn = self.getValues(data,'i',5)
            #self.extractDt = self.extractInt
        #elif self.analysisCode==4: # differential stiffness
            #self.extractDt = self.extractInt
        elif self.analysisCode == 5:   # frequency
            ## frequency
            self.addDataParameter(data, 'freq', 'f', 5)
            self.applyDataCodeValue('dataNames', ['freq'])
        elif self.analysisCode == 6:  # transient
            ## time step
            self.addDataParameter(data, 'time', 'f', 5)
            self.applyDataCodeValue('dataNames', ['time'])
        #elif self.analysisCode==7: # pre-buckling
            #self.extractDt = self.extractInt
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
        #elif self.analysisCode==8: # post-buckling
            #self.extractDt = self.extractInt
            #self.applyDataCodeValue('dataNames',['lsdvmn','eigr'])
        elif self.analysisCode == 9:  # complex eigenvalues
            ## mode number
            self.addDataParameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            #self.addDataParameter(data,'eigr','f',6,False)
            self.applyDataCodeValue('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysisCode == 10:  # nonlinear statics
            ## load factor
            self.addDataParameter(data, 'loadFactor', 'f', 5)
            self.applyDataCodeValue('dataNames', ['loadFactor'])
        #elif self.analysisCode==11: # old geometric nonlinear statics
            #self.extractDt = self.extractInt
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## time step
            self.addDataParameter(data, 'time', 'f', 5)
            self.applyDataCodeValue('dataNames', ['time'])
            #self.extractDt = self.extractInt
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
        else:
            msg = 'invalid analysisCode...analysisCode=%s' % (self.analysisCode)
            raise RuntimeError(msg)

        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #self.printBlock(data)
        self.readTitle()

    def readOGF_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        #tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        #print self.dataCode
        #print "tfsCode=%s" %(tfsCode)

        if self.tableCode == 19:  # grid point forces
            assert self.tableName in ['OGPFB1'], 'tableName=%s tableCode=%s' % (self.tableName, self.tableCode)
            self.readOGF_Data_table19()
        else:
            #self.log.debug('skipping approach/table/format/sortCode=%s on %s-OGF table' %(self.tableName,self.atfsCode))
            self.NotImplementedOrSkip('bad approach/table/format/sortCode=%s on %s-OGF table' % (self.atfsCode, self.tableName))
        ###
        #print self.obj

    def readOGF_Data_table19(self):  # grid point forces
        #isSort1 = self.isSort1()
        if self.numWide == 10:  # real/random
            #if self.thermal==0:
            self.createTransientObject(self.gridPointForces,
                                       gridPointForcesObject)  # real
            self.handleResultsBuffer3(self.readOGF_numWide10,
                                      resultName='gridPointForces')
            #else:
                #self.NotImplementedOrSkip()
            #self.handleResultsBuffer3(self.OUG_RealTable)
        elif self.numWide == 16:  # real/imaginary or mag/phase
            #if self.thermal==0:
            self.createTransientObject(self.gridPointForces,
                                       complexGridPointForcesObject)  # complex
            self.handleResultsBuffer3(self.readOGF_numWide16,
                                      resultName='gridPointForces')
            #else:
                #self.NotImplementedOrSkip()
            #self.handleResultsBuffer3(self.OUG_ComplexTable)
        else:
            raise NotImplementedError('only numWide=10 or 16 is allowed  numWide=%s' % (self.numWide))
        ###

    def readOGF_numWide10(self):
        dt = self.nonlinearFactor
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
        dt = self.nonlinearFactor
        (format1, extract) = self.getOEF_FormatStart()
        format1 += 'i8s12f'
        format1 = bytes(format1)
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data) >= 64:
            eData = self.data[0:4 * 16]
            self.data = self.data[4 * 16:]
            out = unpack(format1, eData)
            (eKey, eid, elemName, f1r, f2r, f3r, m1r, m2r, m3r,
                f1i, f2i, f3i, m1i, m2i, m3i) = out
            eKey = extract(eKey, dt)

            if isMagnitudePhase:
                f1 = polarToRealImag(f1r, f1i)
                m1 = polarToRealImag(m1r, m1i)
                f2 = polarToRealImag(f2r, f2i)
                m2 = polarToRealImag(m2r, m2i)
                f3 = polarToRealImag(f3r, f3i)
                m3 = polarToRealImag(m3r, m3i)
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
        #print self.codeInformation()
        #print self.printBlock(self.data)
        n = 0
        nEntries = len(self.data) // 32
        for i in xrange(nEntries):
            eData = self.data[n:n + 32]
            out = unpack(b'2i6f', eData)
            #nid = (out[0]-self.deviceCode)//10  ## @todo update...
            #print out
            n += 32
            #print "nid = ",nid
