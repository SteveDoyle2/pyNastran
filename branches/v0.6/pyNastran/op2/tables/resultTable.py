from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import copy
from numpy import array

from pyNastran.op2.tables.oug.oug import OUG
from pyNastran.op2.tables.oes_stressStrain.oes import OES
#from pyNastran.op2.tables.oes_stressStrain.oesnlxr import OESNLXR
#from pyNastran.op2.tables.oes_stressStrain.oesnlxd import OESNLXD

from pyNastran.op2.tables.oqg_constraintForces.oqg import OQG
from pyNastran.op2.tables.oef_forces.oef import OEF
from pyNastran.op2.tables.opg_appliedLoads.opg import OPG
from pyNastran.op2.tables.oee_energy.oee import OEE
from pyNastran.op2.tables.ogf_gridPointForces.ogf import OGF
#from pyNastran.op2.tables.hisadd import HISADD - combined with R1TAB for now
from pyNastran.op2.tables.r1tab import R1TAB
from pyNastran.op2.tables.destab import DESTAB
from pyNastran.op2.tables.lama_eigenvalues.lama import LAMA


class ResultTable(OQG, OUG, OEF, OPG, OES, OEE, OGF, R1TAB, DESTAB, LAMA):  # OESNLXR,OESNLXD,

    def __init__(self):
        pass

    def readTableA_DUMMY(self):
        """reads a dummy geometry table"""
        self.iTableMap = {}
        self.readRecordTable('DUMMY')

    def readTableB_DUMMY(self):
        """reads a dummy results table"""
        self.table_name = 'DUMMY'
        table3 = self.readTable_DUMMY_3
        table4Data = self.readDUMMY_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OPG()

    def readTable_DUMMY_3(self, iTable):
        """sets dummy parameters"""
        self.analysis_code = None
        self.table_code = None
        self.format_code = None
        self.sort_code = None

    def readDUMMY_Data(self):
        """creates a dummy object and skips the results"""
        self.obj = None
        self.readOES_Element()

    def updateSort1(self):
        extract = self.extractSort1
        return 'i', extract

    def updateSort2(self):
        extract = self.extractSort2
        return 'f'

    def extractSort1(self, eidDevice, dt):
        #eidDevice, = unpack(b'i',data)
        #print "eidDevice=%s dt=%s eid-dev=%s out=%s" %(eidDevice,dt,eidDevice-self.device_code,(eidDevice-self.device_code)/10.)
        return (eidDevice - self.device_code) // 10

    def extractSort2(self, timeFreq, eid):
        #print "timeFreq=%s eid=%s" %(timeFreq,eid)
        #grid_device, = unpack(b'i',data)
        return timeFreq

    def add_data_parameter(self, data, name, Type, field_num, applyNonlinearFactor=True,
                           fixDeviceCode=False):
        """
        Extracts the binary value and adds it to the self.data_code dictionary.

        :param data: the binary stream
        :param name: the name of the parameter
        :param Type: the type of the value 'i'=integer, 'f'=float, 'ssss'=string
        :param field_num: the field where the parameter is located
        :param applyNonlinearFactor: this parameter is the major nonlinear
               parameter (e.g. time, frequency, mode) and will be applied
               to self.nonlinear_factor; default=True (only 1 per nonlinear subcase)
        :param  fixDeviceCode: removes the device_code from the parameter
                        (e.g. eid = 201 -> 20); default=False

        >>> self.mode = self.get_values(data,'i',5) ## mode number
        >>>
        """
        value = self.get_values(data, Type, field_num)
        if self.make_op2_debug:
            self.op2_debug.write('add_data_parameter - %s=%r\n' % (name, value))
        if fixDeviceCode:
            value = (value - self.device_code) // 10
        setattr(self, name, value)
        self.data_code[name] = value
        #if name == 'mode':
            #print("table=%s name=%s Type=%s field_num=%s aCode=%s value=%s" %(self.table_name, name, Type, field_num, self.analysis_code, value))

        if applyNonlinearFactor:
            self.nonlinear_factor = value
            self.data_code['nonlinear_factor'] = value
            self.data_code['name'] = name

    def setNullNonlinearFactor(self):
        self.nonlinear_factor = None
        self.data_code['nonlinear_factor'] = None

    def apply_data_code_value(self, Name, value):
        self.data_code[Name] = value

    def create_transient_object(self, storageObj, classObj, debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).

        :param storageObj:  the dictionary to store the object in (e.g. self.bars)
        :param classObj:    the class object to instantiate
        :param debug:       developer debug

        .. note:: dt can also be load_step depending on the class
        """
        if debug:
            print("create Transient Object")
            print("***NF = %s" % self.nonlinear_factor)
            #print "DC = ",self.data_code

        if hasattr(self,'isubcase'):
            if self.isubcase in storageObj:
                #print "updating dt..."
                self.obj = storageObj[self.isubcase]
                #print "obj = ",self.obj.__class__.__name__
                #print self.obj.write_f06(['',''],'PAGE ',1)[0]

                try:
                    self.obj.update_data_code(self.data_code)
                    #self.obj.update_dt(self.data_code,self.nonlinear_factor)
                except:
                    #try:
                        #print "objName = ",self.obj.name()
                    #except:
                        #print "objName = ",self.obj
                    raise
            else:
                #if self.isRegular:
                    #self.obj = classObj(self.data_code,not(self.isRegular),self.isubcase,self.nonlinear_factor)
                #else:
                self.obj = classObj(self.data_code, self.is_sort1(), self.isubcase, self.nonlinear_factor)
                #print "obj2 = ",self.obj.__class__.__name__
            storageObj[self.isubcase] = self.obj
        else:
            if self.ID in storageObj:
                #print "updating dt..."
                self.obj = storageObj[self.ID]
            else:
                storageObj[self.ID] = self.obj

    def createThermalTransientObject(self, resultName, objClass, is_sort1):
        #print resultName
        if self.isubcase in resultName:
            self.obj = resultName[self.isubcase]
            #print "returning isubcase result=%s" %(self.isubcase)
        else:
            self.obj = objClass(self.data_code, is_sort1,
                                self.isubcase, self.nonlinear_factor)
            resultName[self.isubcase] = self.obj
            #print "creating isubcase result=%s" %(self.isubcase)
        #return self.obj

    def read_results_table(self, table3, table4Data, flag=0):
        self.dtMap = {}
        table_name = self.read_table_name(rewind=False)  # OEF
        self.table_init(table_name)
        if self.make_op2_debug:
            self.op2_debug.write("***start of %s table***\n" % table_name)
        #print("table_name = |%r|" % table_name)

        self.read_markers([-1, 7], table_name)
        ints = self.read_int_block()
        #print "*ints = ",ints

        self.read_markers([-2, 1, 0], table_name)  # 7
        buffer_words = self.get_marker()
        #print "1-buffer_words = ",buffer_words,buffer_words*4
        ints = self.read_int_block()
        #print "*ints = ",ints

        markerA = -4
        markerB = 0

        iTable = -3
        self.read_markers([iTable, 1, 0], table_name)

        exitFast = False
        while [markerA, markerB] != [0, 2]:
            self.isBufferDone = False
            #print self.print_section(140)
            #print "reading iTable3=%s" %(iTable)
            #self.obj = None

            ## the results object
            self.obj = None
            ## dt/loadFactor/frequency/load_step value (or None for static)
            self.nonlinear_factor = None
            self.data_code = {}

            n = self.op2.tell()
            marker = self.get_marker()
            self.goto(n)
            if marker != 146:
                self.log.debug("marker = %s" % marker)
                exitFast = True
                break

            table3(iTable)
            self.data_code['table_name'] = self.table_name
            ## developer parameter - Analysis/Table/Format/Sort Codes
            #self.atfsCode = [self.analysis_code, self.table_code,
            #                 self.format_code, self.sort_code]
            #print "self.tellA = ",self.op2.tell()

            self.isMarker = False

            isBlockDone = self.readTable4(table4Data, flag, iTable - 1)
            #self.firstPass = False

            #print "self.tellB = ",self.op2.tell()
            iTable -= 2
            #print "isBlockDone = ",isBlockDone
            #sys.exit('stopping')
            if isBlockDone:
                #print "iTable = ",iTable
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            n = self.n
            #print self.print_section(100)
            self.read_markers([iTable, 1, 0], table_name)
            #self.log.debug("")

        nOld = self.op2.tell()
        #try:
        if not exitFast:
            #print self.print_section(100000)
            self.read_markers([iTable, 1, 0], table_name)
            #self.get_marker()
            #self.get_marker()
            #self.get_marker()
        #except SyntaxError: # Invalid Markers
        #    self.goto(nOld)
            #print self.print_block(self.data)
            #print self.print_section(100)
            #markerZero = self.get_marker()
            #assert markerZero==0
            #self.goto(nOld)
            #print "finished markerZero"
            #return

        #print str(self.obj)
        if self.make_op2_debug:
            self.op2_debug.write("***end of %s table***\n" % table_name)
        del self.dtMap

    def readTable4(self, table4Data, flag, iTable):
        """loops over repeated table -4s"""
        #self.read_markers([iTable,1,0])
        markerA = 4

        while markerA is not None:
            self.markerStart = copy.deepcopy(self.n)
            #self.print_section(180)
            self.read_markers([iTable, 1, 0])
            #print "starting OEF table 4..."
            if flag:
                isTable4Done, isBlockDone = table4Data(iTable)
            else:
                isTable4Done, isBlockDone = self.readTable4DataSetup(
                    table4Data, iTable)
            if isTable4Done:
                #print "done with OEF4"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            #print "finished reading oef table..."
            markerA = self.get_marker('A')
            self.n -= 12
            self.op2.seek(self.n)

            self.n = self.op2.tell()
            #print "***markerA = ",markerA

            iTable -= 1
            #print "isBlockDone = ",isBlockDone
        #print "isBlockDone = ",isBlockDone
        return isBlockDone

    def readTable4DataSetup(self, table4Data, iTable):  # iTable=-4
        """
        Checks to see if table 4 is done, loads the data, and handles skipping
        """
        isTable4Done = False
        isBlockDone = False

        buffer_words = self.get_marker(self.table_name)
        #print "buffer_words = ",buffer_words
        #print len(buffer_words)
        self.data = self.read_block()
        #self.print_block(data)

        if buffer_words == 146:  # table -4 is done, restarting table -3
            isTable4Done = True
            return isTable4Done, isBlockDone
        elif buffer_words == 0:
            #print "buffer_words 0 - done with Table4"
            isTable4Done = True
            isBlockDone = True
            return isTable4Done, isBlockDone

        isBlockDone = not(buffer_words)

        if self.is_valid_subcase():  # lets the user skip a certain subcase
            table4Data()
        else:
            self.log.debug("***skipping table=%s isubcase=%s" % (self.table_name, self.isubcase))
            self.skipOES_Element()
        return (isTable4Done, isBlockDone)

    def update_dt_map(self):
        """
        Interfaces with setTransientTimes(times) to limit the amount
        of output that is in the transient results.  This allows
        for models with 1000s of time steps that would otherwise
        crash with memory errors to run.  Every result may be extracted
        if the OP2 is read multiple times.

        While not ideal, this function prevents having to drastically
        change the code to support large models, which would
        make the OP2 reader not as useful for reasonably sized models.

        The code works by taking the user-provided array of values
        for a given subcase and b/c it's an array can be subtracted
        from the current value of dt.  Then the minimum absolute value
        is found and that is mapped back to the original value.  If a
        closer value to the target user-specified value is found, the
        previous result will be deleted, which keeps memory usage down.
        The new result will then be read.  If a dt is not closer than
        an existing value to a given case, that case will not be read.
        """
        #numArray = array([1.,2.])
        isubcase = self.isubcase
        numArray = self.expected_times[isubcase]
        #nums = [0.9,1.11,  1.89,2.1]
        # they're sorted so the last value is the current dt
        num = self.obj.get_transients()[-1]

        readCase = True
        delta = numArray - num
        absDelta = list(abs(delta))
        closest = min(absDelta)
        iclose = absDelta.index(closest)
        actualValue = numArray[iclose]

        if isubcase not in self.dtMap:
            self.dtMap[isubcase] = {}
        if iclose in self.dtMap[isubcase]:
            v1 = self.dtMap[isubcase][iclose]
            vact = get_close_num(v1, num, actualValue)

            if vact != self.dtMap[isubcase][iclose]:
                del self.dtMap[isubcase][iclose]
                self.obj.delete_transient(v1)
                #print "num=%s closest=%s iclose=%s" %(num,actualValue,iclose)
                #print "***deleted v1=%s num=%s vact=%s actual=%s" %(v1,num,vact,actualValue)
                self.dtMap[isubcase][iclose] = vact
                readCase = True
                #print self.dtMap
                #print "A"
            else:  # cleanup previous creation of empty dt case (happened in update_dt outside this function)
                readCase = False
                #print self.dtMap
                #print "num=%s closest=%s iclose=%s" %(num,actualValue,iclose)
                #print "B"
                self.obj.delete_transient(num)
        else:  # read case
            self.dtMap[isubcase][iclose] = num
            readCase = True
            #print "num=%s closest=%s iclose=%s" %(num,actualValue,iclose)
            #print self.dtMap
            #print "C"
        #print "delta = ",delta,'\n'
        #print "readCase = ",readCase
        #if num>=0.14:
            #print self.obj.get_transients()
            #sys.exit('OUG !!!')

        return readCase

    def not_implemented_or_skip(self, msg=''):
        """stops if code is in development and continues otherwise"""
        if False:
            raise NotImplementedError(msg)
        else:
            #self.log.info("skipping...")
            #self.log.info("\n" + self.code_information())
            self.skipOES_Element()

    def handle_results_buffer(self, func, resultName, name=None, debug=False):
        """
        Method for getting results without recursion, which for large OP2s
        is necessary to prevent going past the recursion limit.

        Works by knowing that:

        * the end of an unbuffered table has a:
          - [4]
        * the end of an table with a buffer has a:
          - [4,4,x,4] where x is the next buffer size, which may have another
            buffer
        * the end of the final buffer block has
          - nothing!

        :param self:  the object pointer
        :param func:  the function to recursively call
                      (the function that called this)
        :param result_name: the name of the result object the output will go to
        :param name:  the unique name of the case, which will indicate to the
                      user what data is being saved using a GUI (currently unused)
        :param debug: developer debug

        .. note:: The code knows that the large buffer is the default size
          and the only way there will be a smaller buffer is if there are no
          more buffers.  So, the op2 is shifted by 1 word (4 bytes) to account
          for this end shift.  An extra marker value is read, but no big deal.
          Beyond that it's just appending some data to the binary string and
          calling the function that's passed in.

          We append the data since we read the first 10.5 blocks of data into
          memory, but could only extract the first 10.0 blocks.  We take the
          0.5 block and combine it with the next set of blocks until we have
          no more blocks to read.

        .. warning:: Dont modify this without LOTS of testing.
                     It's a VERY important function
        """
        #if resultName not in self.allowedResultNames:
        #    return self.self.skipOES_Element()

        #stopBuffer = False
        i = 0
        if self.make_op2_debug:
            self.op2_debug.write("***start of %s results***\n" % resultName)
            self.op2_debug.write(self.code_information())

        while not self.isBufferDone:
            #print "n=%s len(data)=%s" %(self.n,len(self.data))
            #sys.stdout.flush()
            func()
            nOld = self.n
            markers = self.read_header()
            #print "nOld=%s markers=%s" %(nOld,markers)

            if markers < 0:  # not a buffer, the table may be done
                self.goto(nOld)
                #print "markers%%2 = %s" %(markers%2)
                if markers is not None and markers % 2 == 1:
                    self.isBufferDone = True
            else:
                data = self.read_block()
                if type(data) != type(self.data):
                    msg = 'The function f=%s has a unicode error\n' % f.__name__
                    msg += ("type(self.data)=%s type(data)=%s"
                            % (type(self.data), type(data)))
                sys.stdout.flush()
                self.data += data
            i += 1
            if i == 2000:
                raise RuntimeError('Infinite Loop or a really big model...')
            #print "isBufferDone=%s" %(self.isBufferDone)
        if self.make_op2_debug:
            self.op2_debug.write("***end of %s results***\n" % resultName)




def get_close_num(v1, v2, closePoint):
    numList = [v1, v2]
    delta = array([v1, v2]) - closePoint
    #print "**delta=%s" %(delta)
    absDelta = list(abs(delta))
    closest = min(absDelta)
    iclose = absDelta.index(closest)
    actualValue = numList[iclose]
    return actualValue
