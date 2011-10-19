class Subcase(object):
    def __init__(self,id=0):
        self.id = id
        self.params = {}
        self.sol = None
        #print "\n***adding subcase %s***" %(self.id)

    def addData(self,key,value,options,paramType):
        #print "adding iSubcase=%s key=|%s| value=|%s| options=|%s| paramType=%s" %(self.id,key,value,options,paramType)
        self.params[key] = [value,options,paramType]

    def printParam(self,key,param,printBeginBulk=True):
        """
        Prints a single entry of the a subcase from the global or local
        subcase list.
        @todo SET-type is not supported yet...
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
        elif paramType=='PARAM-type':
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
        else:
            raise Exception((key,param))  # SET-type is not supported yet...
        ###
        #print "msg = |%r|" %(msg)
        return msg

    def finishSubcase(self):
        if 'SUBCASE' in self.params:
            del self.params['SUBCASE']
        #print "self.params %s = %s" %(self.id,self.params)

    def __repr__(self):
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

