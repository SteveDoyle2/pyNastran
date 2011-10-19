from subcase import Subcase

class CaseControlDeck(object):
    def __init__(self,lines,log=None):
        self.log = log
        self.lines = lines
        self.subcases = {0:Subcase(id=0)}
        self.iSubcase = 0
        self.read()

    def read(self):
        i = 0
        lines = self.lines
        while i < len(lines):
            line = lines[i].strip()
            #print "rawLine = |%s|" %(line)
            options = []
            value = None
            key = None
            paramType = None

            if 'SUBCASE' in line:
                (key,iSubcase) = line.split(' ')
                #print "key=|%s| iSubcase=|%s|" %(key,iSubcase)
                self.iSubcase = int(iSubcase)
                paramType = 'SUBCASE-type'
            elif '=' in line: # TITLE, STRESS
                (key,value) = line.strip().split('=')
                key   = key.strip()
                value = value.strip()
                #print "key=|%s| value=|%s|" %(key,value)
                paramType = 'STRESS-type'

                if '(' in key:  # comma may be in line - STRESS
                    sline = key.strip(')').split('(')
                    key = sline[0]
                    options = sline[1].split(',')
                    #print "key=|%s| options=%s" %(key,options)
                    paramType = 'STRESS-type'
                elif ' ' in key and ',' in value: # set
                    (key,ID) = key.split()
                    fivalues = value.split(',') # float/int values
                    
                    # read more lines....
                    if line[-1]==',':
                        i+=1
                        while lines[i]==',':
                            fivalues += lines.split(',')
                            i+=1
                        ###
                    ###
                    value = fivalues
                    paramType = 'SET-type'
                elif ',' in value: # special TITLE = stuffA,stuffB
                    print 'A ??? line = ',line
                    raise Exception(line)
                else:  # TITLE = stuff
                    #print 'B ??? line = ',line
                    pass
                ###
            ### = in line
            elif ',' in line: # param
                (key,value,options) = line.strip().split(',')
                paramType = 'PARAM-type'
            elif ' ' in line: # begin bulk
                (key,value) = line.strip().split(' ')
                paramType = 'BEGIN_BULK-type'
            else:
                print 'C ??? line = ',line
                raise Exception(line)
            ###
            i+=1
            #print "key=|%s| value=|%s| options=|%s| paramType=%s" %(key,value,options,paramType)
            self.addParameterToSubcase(key,value,options,paramType)
            #print "--------------"
        ###
        print "done with while loop...\n"
        
        #print str(self)
        #sys.exit('stopping...')
        self.finishSubcases()
    ###

    def finishSubcases(self):
        """
        removes any unwanted data in the subcase
        """
        for (iSubcase,subcase) in sorted(self.subcases.items()):
            subcase.finishSubcase()
        ###
    ###

    def addParameterToSubcase(self,key,value,options,paramType):
        #print "adding key=|%s| value=|%s| options=|%s| paramType=%s" %(key,value,options,paramType)
        if self.iSubcase not in self.subcases: # initialize new subcase
            #self.iSubcase += 1
            self.subcases[self.iSubcase] = Subcase(id=self.iSubcase)

        subcase = self.subcases[self.iSubcase]
        subcase.addData(key,value,options,paramType)

    def __repr__(self):
        msg = ''
        for (iSubcase,subcase) in sorted(self.subcases.items()):
            msg += str(subcase)
            #print "\n"
            #break
        return msg
    ###
###
