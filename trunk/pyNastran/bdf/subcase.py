import sys
class Subcase(object):
    solCodeMap = {
                144:  101,
             }

    def __init__(self,id=0):
        self.id = id
        self.params = {}
        self.sol = None
        #print "\n***adding subcase %s***" %(self.id)

    def getStressCode(self,key,options,value):
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

    def getFormatCode(self,options,value):
        """@todo not done..."""
        formatCode = 0
        if 'REAL'  in options:  formatCode += 1
        if 'IMAG'  in options:  formatCode += 2
        if 'PHASE' in options:  formatCode += 4
        formatCode = max(formatCode,1)
        return formatCode

    def getSortCode(self,options,value):
        sortCode = 0
        if 'COMPLEX' in options:  sortCode += 1
        if 'SORT2'   in options:  sortCode += 2
        if 'RANDOM'  in options:  sortCode += 4
        return sortCode

    def getDeviceCode(self,options,value):
        deviceCode = 0
        if 'PRINT' in options:  deviceCode += 1
        if 'PLOT' in options:   deviceCode += 2
        if 'PUNCH' in options:  deviceCode += 4
        deviceCode = max(deviceCode,1)
        #if deviceCode==0:
        #    deviceCode=1  # PRINT
        return deviceCode
        
    def getAnalysisCode(self,sol):
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
        print "sol=%s" %(sol)
        approachCode = codes[sol]
        print 'approachCode = ',approachCode
        return approachCode
        
    def getTableCode(self,sol,tableName,options):
        if tableName in ['VECTOR','PRESSURE']:
            tableName = 'DISPLACEMENT' # equivalent tables...

        key = (sol,tableName)
        tables = { #SOL, tableName      tableCode
                  (101,'DISPLACEMENT'): 1,
                  (103,'DISPLACEMENT'): 7, # VECTOR
                 #(144,'DISPLACEMENT'): 1,
                  (145,'DISPLACEMENT'): 1,
                  (146,'DISPLACEMENT'): 1,

                  (101,'FORCE'):        3, # ???
                  (145,'FORCE'):        3, # ???
                  (146,'FORCE'):        3, # ???

                  (101,'STRESS'):       5,# 5/20/21 ???
                  (145,'STRESS'):       5,# 5/20/21 ???
                  (146,'STRESS'):       5,# 5/20/21 ???

                  (101,'STRAIN'):       5,# 5/20/21 ???
                  (145,'STRAIN'):       5,# 5/20/21 ??? flutter
                  (146,'STRAIN'):       5,# 5/20/21 ??? saero

                  (101,'SPCFORCES'):    3,
                  (103,'SPCFORCES'):    3,
                 #(144,'SPCFORCES'):    3,
                  (145,'SPCFORCES'):    3,
                  (146,'SPCFORCES'):    3,

                  (101,'MPCFORCES'):    3,
                  (103,'MPCFORCES'):    3,
                 #(144,'MPCFORCES'):    3,
                  (145,'MPCFORCES'):    3,
                  (146,'MPCFORCES'):    3,

                  (145,'SVECTOR'): 14,

                  (101,'FLUX'):         4,
                  (103,'FLUX'):         4,
                  (159,'FLUX'):         4,
                  (159,'THERMAL'):      3, # 3/4 ???


#STRESS(PLOT) = ALL
#FORCE(PLOT) = ALL
#DISPLACEMENT(PLOT) = ALL
#DESOBJ(MIN) = 200

                 }
        print "key=%s" %(str(key))
        tableCode = tables[key]
        return tableCode
        
    def hasParameter(self,paramName):
        if paramName in self.params:
            return True
        return False

    def getParameter(self,paramName):
        paramName = self.updateParamName(paramName)
        return self.params[paramName][0:2]

    def updateParamName(self,paramName):
        """
        takes an abbreviated name and expands it so the user can type DISP or 
        DISPLACEMENT and get the same answer
        @todo not a complete list
        @warning not tested yet...
        """
        #print 'paramName  = ',paramName
        if   paramName.startswith('DESO'):  paramName = 'DESOBJ'
        elif paramName.startswith('DESS'):  paramName = 'DESSUB'
        elif paramName.startswith('DISP'):  paramName = 'DISPLACEMENT'
        elif paramName.startswith('EXPO'):  paramName = 'EXPORTLID'
        elif paramName.startswith('FREQ'):  paramName = 'FREQUENCY'
        elif paramName.startswith('PRES'):  paramName = 'PRESSURE'
        elif paramName.startswith('SUPO'):  paramName = 'SUPORT1'
        elif paramName.startswith('SVEC'):  paramName = 'SVECTOR'
        #elif paramName.startswith('TEMP'):  paramName = 'TEMPERATURE'  # handled in caseControlDeck.py
        #print '*paramName = ',paramName
        return  paramName

    def _addData(self,key,value,options,paramType):
        key = self.updateParamName(key)
        #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id,key,value,options,paramType)
        if isinstance(value,str) and value.isdigit():
            value = int(value)

        (key,value,options) = self._simplifyData(key,value,options,paramType)
        self.params[key] = [value,options,paramType]

    def _simplifyData(self,key,value,options,paramType):
        if paramType=='SET-type':
            #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id,key,value,options,paramType)
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
            return (key,values2,options)

        elif paramType=='CSV-type':
            #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id,key,value,options,paramType)
            if value.isdigit():  # PARAM,DBFIXED,-1
                value = value
            ###
        else:
            a = 'key=|%s|'       %(key)
            b = 'value=|%s|'     %(value)
            c = 'options=|%s|'   %(options)
            d = 'paramType=|%s|' %(paramType)
            #print "_adding iSubcase=%s %-18s %-12s %-12s %-12s" %(self.id,a,b,c,d)
            if isinstance(value,int) or value is None:
                pass
            elif value.isdigit():  # STRESS = ALL
                value = value
            #else: pass
        return (key,value,options)

    def getOp2Data(self,sol,solmap_toValue):
        self.sol = sol
        label = 'SUBCASE %s' %(self.id)
        op2Params = {'iSubcase':None,'tables':[],'analysisCodes':[],'deviceCodes':[],
                     'sortCodes':[],'tableCodes':[],'label':label,'subtitle':None,
                     'title':None,'formatCodes':[],
                     
                     'stressCodes':[],'thermal':None}
        
        results = ['DISPLACEMENT','EKE','EDE','ELSDCON','ENTHALPY','EQUILIBRIUM','ESE',
                  'FLUX','FORCE','GPFORCE','GPKE','GPSDCON','GPSTRAIN','GPSTRESS',
                   'HOUTPUT','MODALKE','MODALSE','MPCFORCES','NLSTRESS','NOUTPUT',
                   'OLOAD','PFGRID','PFMODE','PFPANEL','RCROSS','RESVEC','SACCELERATION',
                   'SDISPACEMENT','SPCFORCES','STRAIN','STRESS','SVECTOR','SVELOCITY',
                   'THERMAL','VECTOR','VELOCITY','VUGRID','WEIGHTCHECK']
                   
#   SPC = 2
#   LOAD = 123458
#   DISPLACEMENT(SORT1,REAL)=ALL
#   SPCFORCES(SORT1,REAL)=ALL
#   STRESS(SORT1,REAL,VONMISES,BILIN)=ALL
        
        if self.sol==200: # converts from solution 200 to solution 144
            param = self.params['ANALYSIS']
            (value,options,paramType) = param
            
            sol = solmap_toValue[value.upper()]
            print "***value=%s sol=%s" %(value,sol)
        else:  # leaves SOL the same
            sol  = self.sol
        ###
        if sol in self.solCodeMap:  # reduces SOL 144 to SOL 101
            sol = self.solCodeMap[sol]

        thermal = 0
        for key,param in self.params.items():
            key = key.upper()
            (value,options,paramType) = param
            msg = "  *key=|%s| value=|%s| options=%s paramType=|%s|" %(key,value,options,paramType)
            print msg
            #msg += self.printParam(key,param,printBeginBulk=False)
            if paramType=='SUBCASE-type':
                op2Params['iSubcase'].append(value)
            elif key in ['BEGIN','ECHO','ANALYSIS'] or 'SET' in key:
                pass
            elif key=='TEMPERATURE':
                thermal = 1
            elif key in results:
                sortCode   = self.getSortCode(options,value)
                deviceCode = self.getDeviceCode(options,value)
                
                if key in ['STRESS','STRAIN']:
                    stressCode = self.getStressCode(key,options,value)
                    op2Params['stressCodes'].append(stressCode)
                else:
                    op2Params['stressCodes'].append(0)
                ###
                
                formatCode   = self.getFormatCode(options,value)
                tableCode    = self.getTableCode(sol,key,options)
                analysisCode = self.getAnalysisCode(sol)
                
                approachCode = analysisCode*10    + deviceCode
                tCode        = tableCode*1000 + sortCode
                op2Params['tables'].append(key)
                                
                op2Params['analysisCodes'].append(analysisCode)
                #op2Params['approachCodes'].append(approachCode)
                op2Params['deviceCodes'].append(deviceCode)
                op2Params['formatCodes'].append(formatCode)
                op2Params['sortCodes'].append(sortCode)
                op2Params['tableCodes'].append(tableCode)
                #analysisMethod = value

            elif key in ['TITLE','SUBTITLE','LABEL',
                        'LOAD','SUPORT','SUPORT1','MPC','SPC',
                        'TSTEPNL','NLPARM','TRIM','GUST','METHOD','DESOBJ',
                        'DESSUB','FMETHOD',]:
                op2Params[key.lower()] = value
            ###
            else:
                raise Exception('unsupported entry...\n%s' %(msg))
        ###
        op2Params['thermal'] = thermal
        
        print "\nThe estimated results..."
        for key,value in sorted(op2Params.items()):
            if value is not None:
                print "   key=|%s| value=|%s|" %(key,value)
        #sys.exit('exit subcase.py in bdf/subcase.py')
        

    def printParam(self,key,param,printBeginBulk=True):
        """
        Prints a single entry of the a subcase from the global or local
        subcase list.
        """
        msg = ''
        #msg += 'id=%s   ' %(self.id)
        (value,options,paramType) = param
        
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
            msg += '%s\n' %(value)
        elif paramType=='STRING-type':
            msg += spaces+'%s = %s\n' %(key,value)
        elif paramType=='CSV-type':
            msg += spaces+'%s,%s,%s\n' %(key,value,options)
        elif paramType=='STRESS-type':
            sOptions = ','.join(options)
            #print "sOptions = |%s|" %(sOptions)
            #print "STRESSTYPE key=%s value=%s options=%s" %(key,value,options)
            if len(sOptions)>0:
                msg += '%s(%s) = %s\n' %(key,sOptions,value)
            else:
                msg += '%s = %s\n' %(key,value)
            msg = spaces + msg

        elif paramType=='BEGIN_BULK-type':
            msg += '%s %s\n' %(key,value)
            if 'BEGIN BULK' not in msg:
                msg = spaces + msg
            elif printBeginBulk:
                pass
            else:
                msg = ''
            ###
        elif paramType=='SET-type':
            ## @todo collapse data...not written yet
            msg2 = 'SET %s = ' %(options)
            nChars = len(msg2)
            
            i = 0
            while i<len(value):
                #print "type(value[i]) = ",type(value[i])
                newString = '%s, ' %(value[i])
                #print "newString[%i] = |%s|" %(i,newString)
                if len(msg2+newString)>70:
                    msg += msg2+'\n'
                    msg2 = '        '+newString
                else:
                    msg2 += newString
                ###
                i+=1
            ###
            msg += msg2.rstrip(' \n,')+'\n'
        else:
            raise Exception((key,param))  # SET-type is not supported yet...
        ###
        #print "msg = |%r|" %(msg)
        return msg

    def crossReference(self,mesh):
        print "keys = ",sorted(self.params.keys())
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
        if self.id==0:
            msg = str(self)
        else:
            msg = 'SUBCASE %s\n' %(self.id)
            for (key,param) in sorted(self.params.items()):
                if key in subcase0.params and subcase0.params[key]==param:
                    pass
                else:
                    if 'key'=='BEGIN':
                        continue
                    else:
                        #print "key=%s param=%s" %(key,param)
                        (value,options,paramType) = param
                        #print "  *key=|%s| value=|%s| options=%s paramType=|%s|" %(key,value,options,paramType)
                        msg += self.printParam(key,param,printBeginBulk=False)
                        #print ""
                    ###
                ###
        if self.id>0 and 'BEGIN' in self.params:
            msg += self.printParam('BEGIN',self.params['BEGIN'])
        return msg

    def __repr__(self):
        """
        prints out every entry in the subcase
        """
        #msg = "-------SUBCASE %s-------\n" %(self.id)
        msg = ''
        if self.id>0:
            msg += 'SUBCASE %s\n' %(self.id)
        ###

        for (key,param) in sorted(self.params.items()):
            if 'key'=='BEGIN':
                continue
            else:
                #print "key=%s param=%s" %(key,param)
                (value,options,paramType) = param
                #print "  *key=|%s| value=|%s| options=%s paramType=|%s|" %(key,value,options,paramType)
                msg += self.printParam(key,param,printBeginBulk=False)
                #print ""
            ###
        ###
        if self.id>0 and 'BEGIN' in self.params:
            msg += self.printParam('BEGIN',self.params['BEGIN'])
        return msg

