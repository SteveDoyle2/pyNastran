import sys

class Subcase(object):
    solCodeMap = {
                  1:  101,
                 21:  101,
                 24:  101,
                 26:  101,
                 61:  101,
                 64:  106,  # correct
                 66:  106,  # correct
                 68:  106,  # correct
                 76:  101,
                 99:  129,  # correct
                144:  101,  # correct
                187:  101,
             }

    def __init__(self, id=0):
        self.id = id
        self.params = {}
        self.sol = None
        #print "\n***adding subcase %s***" %(self.id)

    def getStressCode(self, key, options, value):
        """@note the individual element must take
        the stressCode and reduce it to what the element can
        return.  For example, for an isotropic CQUAD4
        the fiber field doesnt mean anything.
        
        BAR       - no von mises/fiber
        ISOTROPIC - no fiber
        
        @todo how does the MATERIAL bit get turned on?  I'm assuming it's element dependent...
        """
        stressCode = 0
        if 'VONMISES'  in options:  stressCode += 1
        if key=='STRAIN':           stressCode += 10 # 2+8=10 - fields 2 and 4
        if 'FIBER'     in options:  stressCode += 4
        #if 'MATERIAL' in options:  stressCode += 16  material coord (1) vs element (0)
        return stressCode

    def getFormatCode(self, options, value):
        """@todo not done..."""
        formatCode = 0
        if 'REAL'  in options:  formatCode += 1
        if 'IMAG'  in options:  formatCode += 2
        if 'PHASE' in options:  formatCode += 4
        formatCode = max(formatCode, 1)
        return formatCode

    def getSortCode(self, options, value):
        sortCode = 0
        if 'COMPLEX' in options:  sortCode += 1
        if 'SORT2'   in options:  sortCode += 2
        if 'RANDOM'  in options:  sortCode += 4
        return sortCode

    def getDeviceCode(self, options, value):
        deviceCode = 0
        if 'PRINT' in options:  deviceCode += 1
        if 'PLOT'  in options:  deviceCode += 2
        if 'PUNCH' in options:  deviceCode += 4
        deviceCode = max(deviceCode,1)
        #if deviceCode==0:
        #    deviceCode=1  # PRINT
        return deviceCode
        
    def getAnalysisCode(self, sol):
        """
        8 - post-buckling (maybe 7 depending on NLPARM???)

        # not important
        3/4 - differential stiffness (obsolete)
        11  - old geometric nonlinear statics
        12  - contran (???)
        @todo verify
        """
        codes = {
                 101: 1,  # staics
                 103: 2,  # modes
                 105: 7,  # pre-buckling
                 106: 10, # nonlinear statics
                 107: 9,  # complex eigenvalues
                 108: 5,  # frequency
                 111: 5,
                 112: 6,
                 114: 1,
                 115: 2,
                 116: 7,
                 118: 5,
                 129: 6,
                 144: 1,
                 145: 1,
                 146: 1,
                 153: 10,
                 159: 6,  # transient
                 }
        #print "sol=%s" %(sol)
        approachCode = codes[sol]
        #print 'approachCode = ',approachCode
        return approachCode
        
    def getTableCode(self, sol, tableName, options):
        if tableName in ['VECTOR', 'PRESSURE']:
            tableName = 'DISPLACEMENT' # equivalent tables...

        key = (sol,tableName)
        tables = { #SOL, tableName      tableCode

                  (101,'ACCELERATION'): 11,
                  (103,'ACCELERATION'): 11,
                  (106,'ACCELERATION'): 11,
                  (107,'ACCELERATION'): 11,
                  (108,'ACCELERATION'): 11,
                  (129,'ACCELERATION'): 11,
                 #(144,'ACCELERATION'): 11,
                  (145,'ACCELERATION'): 11,
                  (146,'ACCELERATION'): 11,


                  (101,'DISPLACEMENT'): 1,
                  (103,'DISPLACEMENT'): 7, # VECTOR
                  (105,'DISPLACEMENT'): 7,
                  (106,'DISPLACEMENT'): 1,
                  (107,'DISPLACEMENT'): 7,
                  (108,'DISPLACEMENT'): 1,
                  (109,'DISPLACEMENT'): 1,
                  (111,'DISPLACEMENT'): 7,
                  (112,'DISPLACEMENT'): 1,
                  (129,'DISPLACEMENT'): 7,
                 #(144,'DISPLACEMENT'): 1,
                  (145,'DISPLACEMENT'): 1,
                  (146,'DISPLACEMENT'): 1,

                  (101,'ESE'):          18, # energy
                  (103,'ESE'):          18, # energy
                  (105,'ESE'):          18, # energy
                  (106,'ESE'):          18, # energy
                  (107,'ESE'):          18, # energy
                  (108,'ESE'):          18, # energy
                  (109,'ESE'):          18, # energy
                  (110,'ESE'):          18, # energy
                  (111,'ESE'):          18, # energy
                  (112,'ESE'):          18, # energy
                  (145,'ESE'):          18, # energy
                  (146,'ESE'):          18, # energy

                  (101,'FORCE'):        3, # ???
                  (103,'FORCE'):        3, # ???
                  (105,'FORCE'):        3, # ???
                  (106,'FORCE'):        3, # ???
                  (107,'FORCE'):        4, # ???
                  (108,'FORCE'):        3, # ???
                  (111,'FORCE'):        3, # ???
                  (112,'FORCE'):        3, # ???
                  (129,'FORCE'):        3, # ???
                  (145,'FORCE'):        3, # ???
                  (146,'FORCE'):        3, # ???

                  (101,'GPFORCE'):      19,
                  (105,'GPFORCE'):      19,
                  (106,'GPFORCE'):      19,
                  (107,'GPFORCE'):      19,
                  (108,'GPFORCE'):      19,
                  (111,'GPFORCE'):      19,
                  (112,'GPFORCE'):      19,
                  (129,'GPFORCE'):      19,
                  (145,'GPFORCE'):      19,
                  (146,'GPFORCE'):      19,

                  (101,'GPSTRESS'):     20,
                  (105,'GPSTRESS'):     20,
                  (106,'GPSTRESS'):     20,
                  (107,'GPSTRESS'):     20,
                  (108,'GPSTRESS'):     20,
                  (111,'GPSTRESS'):     20,
                  (112,'GPSTRESS'):     20,
                  (129,'GPSTRESS'):     20,
                  (145,'GPSTRESS'):     20,
                  (146,'GPSTRESS'):     20,

                  (101,'GPSTRAIN'):     21,
                  (105,'GPSTRAIN'):     21,
                  (106,'GPSTRAIN'):     21,
                  (107,'GPSTRAIN'):     21,
                  (108,'GPSTRAIN'):     21,
                  (111,'GPSTRAIN'):     21,
                  (112,'GPSTRAIN'):     21,
                  (129,'GPSTRAIN'):     21,
                  (145,'GPSTRAIN'):     21,
                  (146,'GPSTRAIN'):     21,

                  (101,'MPCFORCES'):    3,
                  (103,'MPCFORCES'):    3,
                  (106,'MPCFORCES'):    3,
                  (108,'MPCFORCES'):    3,
                  (112,'MPCFORCES'):    3,
                  (129,'MPCFORCES'):    3,
                 #(144,'MPCFORCES'):    3,
                  (145,'MPCFORCES'):    3,
                  (146,'MPCFORCES'):    3,

                  (101,'OLOAD'):        2,
                  (103,'OLOAD'):        2,
                  (105,'OLOAD'):        2,
                  (106,'OLOAD'):        2,
                  (107,'OLOAD'):        2,
                  (108,'OLOAD'):        2,
                  (111,'OLOAD'):        2,
                  (112,'OLOAD'):        2,
                  (129,'OLOAD'):        2,
                 #(144,'OLOAD'):        2,
                  (145,'OLOAD'):        2,
                  (146,'OLOAD'):        2,

                  (101,'SPCFORCES'):    3,
                  (103,'SPCFORCES'):    3,
                  (105,'SPCFORCES'):    3,
                  (106,'SPCFORCES'):    3,
                  (107,'SPCFORCES'):    3,
                  (108,'SPCFORCES'):    3,
                  (110,'SPCFORCES'):    3,
                  (111,'SPCFORCES'):    3,
                  (112,'SPCFORCES'):    3,
                  (129,'SPCFORCES'):    3,
                 #(144,'SPCFORCES'):    3,
                  (145,'SPCFORCES'):    3,
                  (146,'SPCFORCES'):    3,

                  (101,'STRAIN'):       5,# 5/20/21 ???
                  (105,'STRAIN'):       5,
                  (106,'STRAIN'):       5,
                  (107,'STRAIN'):       5,
                  (108,'STRAIN'):       5,
                  (110,'STRAIN'):       5,
                  (111,'STRAIN'):       5,
                  (112,'STRAIN'):       5,
                  (129,'STRAIN'):       5,
                  (145,'STRAIN'):       5,
                  (146,'STRAIN'):       5,

                  (101,'STRESS'):       5,# 5/20/21 ???
                  (103,'STRESS'):       5,
                  (105,'STRESS'):       5,
                  (106,'STRESS'):       5,
                  (107,'STRESS'):       5,
                  (108,'STRESS'):       5,
                  (111,'STRESS'):       5,
                  (112,'STRESS'):       5,
                  (129,'STRESS'):       5,
                  (145,'STRESS'):       5,
                  (146,'STRESS'):       5,

                  (145,'SVECTOR'):      14,

                  (101,'FLUX'):         4,
                  (103,'FLUX'):         4,
                  (106,'FLUX'):         4,
                  (112,'FLUX'):         4,
                  (108,'FLUX'):         4,
                  (153,'FLUX'):         4,
                  (159,'FLUX'):         4,

                  (101,'THERMAL'):      3, # 3/4 ???
                  (159,'THERMAL'):      3, # 3/4 ???

                  (101,'VELOCITY'):    10,
                  (103,'VELOCITY'):    10,
                  (106,'VELOCITY'):    10,
                  (107,'VELOCITY'):    10,
                  (108,'VELOCITY'):    10,
                  (111,'VELOCITY'):    10,
                  (112,'VELOCITY'):    10,
                  (129,'VELOCITY'):    10,
                 #(144,'VELOCITY'):    10,
                  (145,'VELOCITY'):    10,
                  (146,'VELOCITY'):    10,

                  (101,'VUGRID'):    10,

                 }
        print("key=%s" %(str(key)))
        if key not in tables:
            raise KeyError(key)
        tableCode = tables[key]
        return tableCode
        
    def hasParameter(self, paramName):
        if paramName in self.params:
            return True
        return False

    def getParameter(self, paramName):
        paramName = self.updateParamName(paramName)
        if paramName not in self.params:
            raise KeyError('%s doesnt exist in subcase=%s in the case control deck.' %(paramName,self.id))
        return self.params[paramName][0:2]

    def updateParamName(self, paramName):
        """
        takes an abbreviated name and expands it so the user can type DISP or 
        DISPLACEMENT and get the same answer
        @todo not a complete list
        @warning fully tested yet...
        """
        #print 'paramName  = ',paramName
        if   paramName.startswith('ACCE'):  paramName = 'ACCELERATION'
        elif paramName.startswith('DESO'):  paramName = 'DESOBJ'
        elif paramName.startswith('DESS'):  paramName = 'DESSUB'
        elif paramName.startswith('DISP'):  paramName = 'DISPLACEMENT'
        elif paramName.startswith('EXPO'):  paramName = 'EXPORTLID'
        elif paramName.startswith('ELFO'):  paramName = 'FORCE'
        elif paramName.startswith('FORC'):  paramName = 'FORCE'
        elif paramName.startswith('FREQ'):  paramName = 'FREQUENCY'
        elif paramName.startswith('GPFO'):  paramName = 'GPFORCE'
        elif paramName=='GPST':             raise SyntaxError('invalid GPST stress/strain')
        elif paramName.startswith('METH'):  paramName = 'METHOD'
        elif paramName.startswith('MPCF'):  paramName = 'MPCFORCES'
        elif paramName.startswith('OLOA'):  paramName = 'OLOAD'
        elif paramName.startswith('PRES'):  paramName = 'PRESSURE'

        elif paramName.startswith('SPCF'):  paramName = 'SPCFORCES'

        elif paramName.startswith('STRA'):  paramName = 'STRAIN'
        elif paramName.startswith('STRE'):  paramName = 'STRESS'
        elif paramName.startswith('SUPO'):  paramName = 'SUPORT1'
        elif paramName.startswith('SVEC'):  paramName = 'SVECTOR'
        elif paramName.startswith('THER'):  paramName = 'THERMAL'
        elif paramName.startswith('VECT'):  paramName = 'VECTOR'
        elif paramName.startswith('VELO'):  paramName = 'VELOCITY'



        #elif paramName.startswith('DFRE'):  paramName = 'D'
#DFRE
        #elif paramName.startswith('TEMP'):  paramName = 'TEMPERATURE'  # handled in caseControlDeck.py
        #print '*paramName = ',paramName
        return  paramName

    def _addData(self, key, value, options, paramType):
        key = self.updateParamName(key)
        #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id, key, value, options, paramType)
        if isinstance(value,str) and value.isdigit():
            value = int(value)

        (key, value, options) = self._simplifyData(key, value, options, paramType)
        self.params[key] = [value, options, paramType]

    def _simplifyData(self,key, value, options, paramType):
        if paramType=='SET-type':
            #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id, key, value, options, paramType)
            values2 =[]
            for i,ivalue in enumerate(value):
                ivalue = ivalue.strip()
                if ivalue.isdigit():
                    values2.append(int(ivalue))
                else:
                    if value is 'EXCLUDE':
                        raise RuntimeError('EXCLUDE is not supported on CaseControlDeck SET card\n')
                    values2.append(ivalue)
                ###
            ###
            
            ## @todo expand values with THRU and EXCLUDE
            ## @todo sort values
            ## @todo collapse values when printing
            
            #print "values2 = ",values2
            options = int(options)
            return (key, values2, options)

        elif paramType=='CSV-type':
            #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id, key, value, options, paramType)
            if value.isdigit():  # PARAM,DBFIXED,-1
                value = value
            ###
        else:
            #a = 'key=|%s|'       %(key)
            #b = 'value=|%s|'     %(value)
            #c = 'options=|%s|'   %(options)
            #d = 'paramType=|%s|' %(paramType)
            #print("_adding iSubcase=%s %-18s %-12s %-12s %-12s" %(self.id, a, b, c, d))
            if isinstance(value, int) or value is None:
                pass
            elif value.isdigit():  # STRESS = ALL
                value = value
            #else: pass
        return (key,value,options)

    def getOp2Data(self,sol,solmap_toValue):
        self.sol = sol
        label = 'SUBCASE %s' %(self.id)
        op2Params = {'iSubcase':None, 'tables':[], 'analysisCodes':[], 
                     'deviceCodes':[], 'sortCodes':[], 'tableCodes':[],
                     'label':label,'subtitle':None, 'title':None,
                     'formatCodes':[], 'stressCodes':[], 'thermal':None}
        
        results = ['DISPLACEMENT', 'EKE', 'EDE', 'ELSDCON', 'ENTHALPY',
                   'EQUILIBRIUM', 'ESE', 'FLUX', 'FORCE', 'GPFORCE', 'GPKE',
                   'GPSDCON', 'GPSTRAIN', 'GPSTRESS', 'HOUTPUT', 'MODALKE',
                   'MODALSE', 'MPCFORCES', 'NLSTRESS', 'NOUTPUT', 'OLOAD',
                   'PFGRID', 'PFMODE', 'PFPANEL', 'RCROSS', 'RESVEC',
                   'SACCELERATION', 'SDISPACEMENT', 'SPCFORCES', 'STRAIN',
                   'STRESS', 'SVECTOR', 'SVELOCITY', 'THERMAL', 'VECTOR',
                   'VELOCITY', 'VUGRID', 'WEIGHTCHECK']

        # converts from solution 200 to solution 144
        if self.sol==200 or 'ANALYSIS' in self.params:
            param = self.params['ANALYSIS']
            (value, options, paramType) = param
            
            sol = solmap_toValue[value.upper()]
            print("***value=%s sol=%s" %(value, sol))
        else:  # leaves SOL the same
            sol  = self.sol
        ###
        if sol in self.solCodeMap:  # reduces SOL 144 to SOL 101
            sol = self.solCodeMap[sol]

        for key,param in self.params.iteritems():
            key = key.upper()
            (value, options, paramType) = param
            msg = "  -key=|%s| value=|%s| options=%s paramType=|%s|" %(key, value, options, paramType)

        thermal = 0
        for key,param in self.params.iteritems():
            key = key.upper()
            (value, options, paramType) = param
            msg = "  *key=|%s| value=|%s| options=%s paramType=|%s|" %(key, value, options, paramType)
            print(msg)
            #msg += self.printParam(key,param,printBeginBulk=False)
            if paramType=='SUBCASE-type':
                op2Params['iSubcase'].append(value)
            elif key in ['BEGIN','ECHO','ANALYSIS'] or 'SET' in key:
                pass
            elif key=='TEMPERATURE':
                thermal = 1
            elif key in results:
                sortCode   = self.getSortCode(options, value)
                deviceCode = self.getDeviceCode(options, value)
                
                if key in ['STRESS','STRAIN']:
                    stressCode = self.getStressCode(key, options, value)
                    op2Params['stressCodes'].append(stressCode)
                else:
                    op2Params['stressCodes'].append(0)
                ###
                
                formatCode   = self.getFormatCode(options,value)
                tableCode    = self.getTableCode(sol,key,options)
                analysisCode = self.getAnalysisCode(sol)
                
                approachCode = analysisCode*10 + deviceCode
                tCode        = tableCode*1000 + sortCode
                op2Params['tables'].append(key)
                                
                op2Params['analysisCodes'].append(analysisCode)
                #op2Params['approachCodes'].append(approachCode)
                op2Params['deviceCodes'].append(deviceCode)
                op2Params['formatCodes'].append(formatCode)
                op2Params['sortCodes'].append(sortCode)
                op2Params['tableCodes'].append(tableCode)
                #analysisMethod = value

            #elif key in ['ADACT','ADAPT','AERCONFIG',
            #            'TITLE','SUBTITLE','LABEL',
            #            'LOAD','SUPORT','SUPORT1','MPC','SPC',
            #            'TSTEPNL','NLPARM','TRIM','GUST','METHOD','DESOBJ',
            #            'DESSUB','FMETHOD','SEALL']:
            else:
                op2Params[key.lower()] = value
            ###
            #else:
            #    raise Exception('unsupported entry...\n%s' %(msg))
            ###
        ###
        op2Params['thermal'] = thermal
        
        print("\nThe estimated results...")
        for (key, value) in sorted(op2Params.iteritems()):
            if value is not None:
                print("   key=|%s| value=|%s|" %(key, value))
        #sys.exit('exit subcase.py in bdf/subcase.py')
        

    def printParam(self, key, param, printBeginBulk=True):
        """
        Prints a single entry of the a subcase from the global or local
        subcase list.
        """
        msg = ''
        #msg += 'id=%s   ' %(self.id)
        (value, options, paramType) = param
        
        spaces = ''
        if self.id>0:
            spaces = '    '

        if paramType=='SUBCASE-type':
            if self.id>0:
                msg += 'SUBCASE %s\n' %(self.id)
            ###
            #else:  global subcase ID=0 and is not printed
            #    pass
        elif paramType=='KEY-type':
            #print "KEY-TYPE:  |%s|" %(value)
            assert value is not None, param
            msg += spaces + '%s\n' %(value)
        elif paramType=='STRING-type':
            msg += spaces+'%s = %s\n' %(key, value)
        elif paramType=='CSV-type':
            msg += spaces+'%s,%s,%s\n' %(key, value, options)
        elif paramType=='STRESS-type':
            sOptions = ','.join(options)
            #print "sOptions = |%s|" %(sOptions)
            #print "STRESSTYPE key=%s value=%s options=%s" %(key,value,options)
            if len(sOptions)>0:
                msg += '%s(%s) = %s\n' %(key, sOptions, value)
            else:
                msg += '%s = %s\n' %(key, value)
            msg = spaces + msg

        elif paramType=='BEGIN_BULK-type':
            msg += '%s %s\n' %(key, value)
            if 'BEGIN BULK' not in msg:
                msg = spaces + msg
            elif printBeginBulk:
                pass
            else:
                msg = ''
            ###
        elif paramType=='SET-type':
            ## @todo collapse data...not written yet
            starter = 'SET %s = ' %(options)
            msg2 = spaces + starter
            nChars = len(msg2)
            
            i = 0
            while i<len(value):
                #print "type(value[i]) = ",type(value[i])
                newString = '%s, ' %(value[i])
                #print "newString[%i] = |%s|" %(i,newString)
                if len(msg2+newString)>70:
                    msg += msg2+'\n'
                    msg2 = ' '*nChars+newString
                else:
                    msg2 += newString
                ###
                i+=1
            ###
            msg += msg2.rstrip(' \n,')+'\n'
        else:
            raise NotImplementedError((key,param))  # SET-type is not supported yet...
        ###
        #print "msg = |%r|" %(msg)
        return msg

    def crossReference(self,mesh):
        """
        @note this is not integrated and probably never will be as it's not really that necessary
        """
        print("keys = %s" %(sorted(self.params.keys())))
        if 'LOAD' in self.params:
            loadID = self.params['LOAD'][0]
            loadObj = mesh.loads[loadID]
            loadObj.crossReference(mesh)
        if 'SUPORT' in self.params:
            pass
        if 'MPC' in self.params:
            #mpcID = self.params['MPC'][0]
            #mpcObj = mesh.mpcs[mpcID]
            #mpcObj.crossReference(mesh)
            pass
        if 'SPC' in self.params:
            #spcID = self.params['SPC'][0]
            #print "SPC ID = ",spcID
            #spcObj = mesh.spcObject
            #spcObj.crossReference(spcID,mesh)
            pass
        if 'TSTEPNL' in self.params:
            tstepnlID = self.params['TSTEPNL'][0]
            tstepnlObj = mesh.tstepnl[tstepnlID]
            tstepnlObj.crossReference(mesh)
        if 'NLPARM' in self.params:
            nlparmID = self.params['NLPARM'][0]
            nlparmObj = mesh.nlparms[nlparmID]
            nlparmObj.crossReference(mesh)
        if 'TRIM' in self.params:
            trimID = self.params['TRIM'][0]
            trimObj = mesh.trims[trimID]
            trimObj.crossReference(mesh)
        if 'GUST' in self.params:
            gustID = self.params['GUST'][0]
            gustObj = mesh.gusts[gustID]
            gustObj.crossReference(mesh)
        if 'DLOAD' in self.params: # ???
            pass
    
    def finishSubcase(self):
        """
        removes the subcase parameter from the subcase to avoid printing it in a funny spot
        """
        if 'SUBCASE' in self.params:
            del self.params['SUBCASE']
        #print "self.params %s = %s" %(self.id,self.params)

    def writeSubcase(self,subcase0):
        """
        internal method to print a subcase
        @param self the subcase object
        @param subcase0 the global subcase
        """
        if self.id==0:
            msg = str(self)
        else:
            msg = 'SUBCASE %s\n' %(self.id)
            for (key,param) in self.subcaseSorted(self.params.items()):
                if key in subcase0.params and subcase0.params[key]==param: # dont write global subcase parameters
                    pass
                else:
                    if 'key'=='BEGIN':
                        continue
                    else:
                        #print "key=%s param=%s" %(key,param)
                        (value,options,paramType) = param
                        #print "  *key=|%s| value=|%s| options=%s paramType=|%s|" %(key,value,options,paramType)
                        msg += self.printParam(key, param, printBeginBulk=False)
                        #print ""
                    ###
                ###
        if self.id>0 and 'BEGIN' in self.params: # self.id>0 and 'BEGIN' in self.params... why was that there???
            msg += self.printParam('BEGIN', self.params['BEGIN'])
        return msg

    def subcaseSorted(self,listA):
        """
        does a "smart" sort on the keys such that SET cards increment in
        numerical order.
        @param self the subcase object
        @param listA the list of subcase list objects
        @retval listB the sorted version of listA
        """
        # presort the list to put all the SET cards next to each other
        listA.sort()
        
        i=0 # needed in case the loop doesnt execute
        iSet=None # index of first SET card in the deck
        setDict    = {} 
        listBefore = []
        listAfter  = []
        setKeys    = []
        for (i, entry) in enumerate(listA):
            key = entry[0]
            if 'SET' in key[0:3]:
                if key=='SET': # handles "SET = ALL"
                    key = 0
                else: # handles "SET 100 = 1,2,3"
                    sline = key.split(' ')
                    key = int(sline[1])

                # store the integer ID and the SET-type list
                setDict[key] = entry
                setKeys.append(key)
            else:
                # only store the entries before the SET cards
                listBefore.append(entry)
                if iSet:
                    break
        
        # grab the other entries
        listAfter = listA[i:]
        
        # write the SET cards in a sorted order
        setList = []
        for key in sorted(setKeys):
            setList.append(setDict[key])
        
        # combine all the cards
        listB = listBefore+setList+listAfter
        return listB

    def __repr__(self):
        """
        prints out every entry in the subcase
        """
        #msg = "-------SUBCASE %s-------\n" %(self.id)
        msg = ''
        if self.id>0:
            msg += 'SUBCASE %s\n' %(self.id)
        ###

        for (key,param) in self.subcaseSorted(self.params.items()):
            if 'key'=='BEGIN':
                continue
            else:
                #print "key=%s param=%s" %(key,param)
                (value,options,paramType) = param
                #print "  *key=|%s| value=|%s| options=%s paramType=|%s|" %(key,value,options,paramType)
                msg += self.printParam(key, param, printBeginBulk=False)
                #print ""
            ###
        ###
        if self.id>0 and 'BEGIN' in self.params: # prevents 2 BEGIN BULKs
            msg += self.printParam('BEGIN', self.params['BEGIN'])
        return msg

