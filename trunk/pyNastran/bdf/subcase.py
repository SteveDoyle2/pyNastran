import sys
class Subcase(object):
    def __init__(self,id=0):
        self.id = id
        self.params = {}
        self.sol = None
        #print "\n***adding subcase %s***" %(self.id)

    def hasParameter(self,paramName):
        if paramName in self.params:
            return True
        return False

    def getParameter(self,paramName):
        return self.params[paramName][0:2]

    def updateParamName(self,paramName):
        """
        takes an abbreviated name and expands it so the user can type DISP or 
        DISPLACEMNT and get the same answer
        @todo not a complete list
        @warning not implemented yet...
        """
        if   paramName.startswith('SUPO'):  paramName = 'SUPORT1'
        elif paramName.startswith('DISP'):  paramName = 'DISPLACEMENT'
        elif paramName.startswith('TEMP'):  paramName = 'TEMPERATURE'
        elif paramName.startswith('FREQ'):  paramName = 'FREQUENCY'
        return  paramName

    def _addData(self,key,value,options,paramType):
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
            #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id,key,value,options,paramType)
            if isinstance(value,int):
                pass
            elif value.isdigit():  # STRESS = ALL
                value = value
            #else: pass
        return (key,value,options)

    def printParam(self,key,param,printBeginBulk=True):
        """
        Prints a single entry of the a subcase from the global or local
        subcase list.
        """
        msg = ''
        (value,options,paramType) = param
        
        spaces = '    '*self.id
        if paramType=='SUBCASE-type':
            if self.id>0:
                msg += 'SUBCASE %s\n' %(self.id)
            ###
            #else:  global subcase ID=0 and is not printed
            #    pass
        elif paramType=='CSV-type':
            msg += spaces+'%s,%s,%s\n' %(key,value,options)
        elif paramType=='STRESS-type':
            sOptions = ','.join(options)
            #print "sOptions = |%s|" %(sOptions)
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

    def finishSubcase(self):
        """
        removes the subcase parameter from the subcase to avoid printing it in a funny spot
        """
        if 'SUBCASE' in self.params:
            del self.params['SUBCASE']
        #print "self.params %s = %s" %(self.id,self.params)

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

