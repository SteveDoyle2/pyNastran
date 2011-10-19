import sys
import copy
from nastranParser.bdf.fieldWriter import printCard
from BDF_Card import BDF_Card


class getMethods(object):
    def getNodeIDs(self):
        return sorted(self.nodes.keys())

    def propertyIDs(self):
        return self.properties.keys()

    def nodeIDs(self):
        return self.nodes.keys()
        
    def getElementIDs(self):
        return sorted(self.elements.keys())

    def elementIDs(self):
        return self.elements.keys()

    def getPropertyIDs(self):
        return sorted(self.properties.keys())

    def getMaterialIDs(self):
        return sorted(self.materials.keys())

    def getNodes(self):
        nodes = []
        for nid,node in sorted(self.nodes.items()):
            nodes.append(node)
        return nodes

    def getNodeIDsWithElements(self,eids):
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid)
            self.log().info("element.pid = %s" %(element.pid))
            nids = set(element.nids)
            nids2 = nids2.union(nids)
        ###
        return nids2

    def getElementIDsWithPIDs(self,pids):
        self.log().info("pids = %s" %(pids))

        eids = self.elementIDs()
        eids2 = []
        for eid in eids:
            element = self.Element(eid)
            #print "element = ",element
            #print "element.pid = ",element.pid
            if element.pid in pids:
                eids2.append(eid)
            ###
            #print ""
        ###
        return (eids2)

    def Node(self,nid):
        return self.nodes[nid]

    def Element(self,eid):
        return self.elements[eid]
    def Load(self,lid):
        return self.loads[lid]

    def Property(self,pid):
        return self.properties[pid]

    def Material(self,mid):
        return self.materials[mid]

    def Coord(self,cid):
        return self.coords[cid]

class addMethods(object):
    def addParam(self,param):
        self.params[param.key] = param

    def addNode(self,node):
        self.nodes[node.id] = node

    def addElement(self,elem):
        self.elements[elem.eid] = elem

    def addProperty(self,prop):
        self.properties[prop.pid] = prop

    def addMaterial(self,material):
        self.materials[material.mid] = material

    def addCoord(self,coord):
        self.coords[coord.cid] = coord

    def addLoad(self,load):
        key = load.lid
        #print "type(self.loads) = ",type(self.loads)
        if self.loads.has_key(key):
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]

    def addConstraint(self,constraint):
        key = constraint.id
        #print "type(self.constraints) = ",type(self.constraints)
        if self.constraints.has_key(key):
            self.constraints[key].append(constraint)
        else:
            self.constraints[key] = [constraint]

class writeMesh(object):
    def writeHeader(self):
        msg = '$EXECUTIVE CONTROL DECK\n'
        for line in self.executiveControlLines:
            msg += line

        msg += '$CASE CONTROL DECK\n'
        msg += str(self.caseControlDeck)
        #for line in self.caseControlLines:
        #    msg += line
        return msg

    def writeParams(self):
        msg = '$PARAMS\n'
        #print "self.nodes = ",self.nodes
        for key,param in sorted(self.params.items()):
            #print "param = ",param
            msg += str(param)
        return msg

    def writeNodes(self):
        msg = '$NODES\n'
        #print "self.nodes = ",self.nodes
        print "nNodes = ",len(self.nodes)
        for key,node in sorted(self.nodes.items()):
            #print "node = ",node
            msg += str(node)
        return msg

    def writeElements(self):
        msg = '$ELEMENTS\n'
        for key,element in sorted(self.elements.items()):
            msg += str(element)
        return msg

    def writeProperties(self):
        msg = '$PROPERTIES\n'
        for key,prop in sorted(self.properties.items()):
            msg += str(prop)
        return msg

    def writeMaterials(self):
        msg = '$MATERIALS\n'
        for key,material in sorted(self.materials.items()):
            msg += str(material)
        return msg

    def writeConstraints(self):
        msg = '$CONSTRAINTS\n'
        for key,loadcase in sorted(self.constraints.items()):
            for constraint in loadcase:
                msg += str(constraint)
        return msg

    def writeLoads(self):
        msg = '$LOADS\n'
        for key,loadcase in sorted(self.loads.items()):
            for load in loadcase:
                msg += str(load)
        return msg

    def writeCoords(self):
        #print "output coords..."
        msg = '$COORDS\n'
        coordKeys = self.coords.keys()
        self.log().info("coordKeys = %s" %(coordKeys))
        for ID,coord in sorted(self.coords.items()):
            if ID!=0:
                msg += str(coord)
        return msg

    def writeRejects(self):
        msg = '$REJECTS\n'
        for rejectCard in self.rejectCards:
            #print "rejectCard = ",rejectCard
            #print ""
            msg += printCard(rejectCard)
        for rejectLines in self.rejects:
            if rejectLines[0][0]==' ':
                continue
            else:
                #print "rejectLines = ",rejectLines
                for reject in rejectLines:
                    reject2 = reject.rstrip()
                    if reject2:
                        msg += str(reject2)+'\n'
                    ###
                ###
            ###
        ###
        return msg

class cardMethods(object):

    def getCard(self,debug=False):
        """gets a single unparsed card"""
        #debug = False
        
        emptyLines=0
        while len(self.lines)<20:
            line=self.infile.readline()
            line = line.split('$')[0]
            if('$' not in line and len(line)>0):
                if debug:
                    print "line = |%s|" %(line)
                self.lines.append(line)
            else:
                emptyLines += 1
            ###
            if emptyLines==50:
                break

        tempcard = [self.lines[0]]
        #if 'CQUAD4' in tempcard[0]:
        #    debug = True
        i=1
        #if emptyLines<50:
        try:
            tempcard = self.getMultiLineCard(i,tempcard,debug=False)
        except IndexError:
            #try:
            #    tempcard = self.getMultiLineCard(i,tempcard,debug=True)
            #except IndexError:
            #    print "workscard = ",tempcard
            #    print ""
            #    #raise
            #    pass
            pass
            #raise
            #print "error..."
            #print "line = ",self.lines[i]
            self.doneReading=True
            #raise
        ###
        
        
        #try:
        self.lines = self.lines[i:]
        #except IndexError:
        #    self.lines = []
        
        #print "tempcard = ",''.join(tempcard)
        
        upperCard = [line.upper() for line in tempcard]
        cardName = self.getCardName(upperCard)
        #print "|%s|" %(cardName)
        #if cardName=='CTRIA3':
        if debug:
            self.log().debug("cardName  = |%s|" %(cardName))
            self.log().debug("upperCard = |%s|" %(upperCard))
            self.log().debug("tempcard  = |%s|" %(tempcard))
            self.log().debug("-------\n")
        return (upperCard,cardName)
    
    def getMultiLineCard(self,i,tempcard,isCSV=False,debug=False):
        if debug:
            print "tempcard1 = ",tempcard
        iline = self.lines[i].rstrip()
        #while iline=='':
        #    i+=1
        #    iline = self.lines[i].rstrip()

        sCardName = iline[0:8].strip()  # trying to find if it's blank...
        if debug:
            self.log().debug("  iline = |%s|" %(iline))
            self.log().debug("  sCardName = |%s|" %(sCardName))
            self.log().debug("  len(iline) = |%s|" %(len(iline)))
            #print "  iline[0] = ",iline[0]
            #print ""

        isNotDone = len(iline)>0 and ((iline[0]=='*') or (iline[0]=='+') or (iline[0]==',') or sCardName=='')
        while(isNotDone):
            tempcard.append(iline)
            i+=1
            #if debug:
            #    print card
            if 'ENDDATA' in iline:
                print "********"
                self.doneReading=True
                break
            if i+1==len(self.lines):
                break
            iline = self.lines[i]
            #try:
            #    iline = self.lines[i]
            #except IndexError:
            #    iline = ''

            sCardName = iline[0:8].strip()  # trying to find if it's blank...
            if debug:
                print "sCardName = ",sCardName
                print "len(iline) = ",len(iline)
            isNotDone = len(iline)>0 and ((iline[0]=='*') or (iline[0]=='+') or sCardName=='')
            #print "isNotDone = ",isNotDone
        ###
        if debug:
            self.log().debug("tempcard2 = |%s|" %(tempcard))
            #print ""
        return tempcard

    def nastranSplit(self,line,isLargeField):
        fields = []
        nChars = len(line)
        iStart = 0
        iEnd = 8
        lenField = 8
        nFields = int(ceil(len(line)/8.))
        #print "nFields = ",nFields
        for iField in range(nFields):
            field = line[iStart:iEnd]
            fields.append(field)
            iStart += lenField
            iEnd   += lenField
        ###
        #print "fields2 = ",fields
        return fields
    
    def nastranSplit2(self,line,isLargeField,debug=False):
        if isLargeField:
            #print "large"
            #print "line = |%s|" %(line)
            fields = [line[0:8],line[8:24],line[24:40],line[40:56],line[56:72]]
        else:
            #print "small"
            #fields = [line[0 :8 ],line[8 :16],line[16:24],line[24:32],line[32:40],
            #          line[40:48],line[48:56],line[56:64],line[64:72],line[72:80]]
            fields = [line[0 :8 ],line[8 :16],line[16:24],line[24:32],line[32:40],
                      line[40:48],line[48:56],line[56:64],line[64:72],line[72:80]]
        if debug:
            print "  fields = ",collapse(fields)
        
        fields2 = []
        for field in fields:
            field = field.strip()
            if '*'==field:
                pass
            else:
                fields2.append(field)
        return fields2

    def isLargeField(self,card):
        """returns True if the card is in 16-character width fields"""
        starField = ['*' in field for field in card]
        #print "starField = ",starField
        
        if any(starField):
            return True
        return False

    def processCard(self,tempcard,debug=False):
        """
        takes a list of strings and returns a list with the 
        proper value in the fields of the list
        """
        #debug = True
        card = []
        isLargeField = self.isLargeField(tempcard)
        #print "*** isLargeField = ",isLargeField

        #print "tempcard = ",tempcard
        for line in tempcard:
            if debug:
                self.log().debug("  line = %s" %(line))
            sline = line.strip('\r\n')
            if not(sline):
                break
            if debug:
                self.log().debug("  sline = %s" %(sline))
            if ',' in sline:
                sline = sline.split(',')
            else:
                
                sline = self.nastranSplit2(sline,isLargeField,debug=debug)
            #name = sline[0]
            #nFields = len(sline)
            #print "sline = ",sline
            
            cardName = sline[0].strip(' *')
            card.append(cardName)
            for valueIn in sline[1:]:
                value = self.getValue(valueIn,debug=debug)
                if debug:
                    print ""
                card.append(value)
            ###
        ###
        if debug:
            self.log().debug("  sline2 = %s" %(card))
            self.log().debug("  sline2 = %s" %(collapse(card)))
        return self.makeSingleStreamedCard(card)
        
    def makeSingleStreamedCard(self,card):
        """
        takes a card that has been split b/c it's a multiline card
        and gets rid of the required blanks in it.
        """
        cardOut = []
        n=0
        #print "card = ",card
        for i,field in enumerate(card):
            if n-9==0:
                pass
            elif n-10==0:
                n-=10
            else:
                cardOut.append(field)
            n+=1
        ###
        #print "cardOut = ",cardOut
        return cardOut
        return collapse(cardOut)
        
    def getValue(self,valueRaw,debug=False):
        """converts a value from nastran format into python format."""
        if debug:
            print "v1 = |%s|" %(valueRaw)
        valueIn = valueRaw.strip(' \n\r').rstrip('*').upper()
        
        if debug:
            pass
            #print "v2 = |%s|" %(valueIn)
        if len(valueIn)==0:
            if debug:
                print "BLANK!"
            return ''

        if '=' in valueIn or '(' in valueIn or '*' in valueRaw:
            if debug:
                print "=(! - special formatting"
            return valueRaw.strip()
        #valueIn = valueIn.upper()
        # int, float, string, exponent
        valuePositive = valueIn.strip('+-')
        if debug:
            print "isDigit = ",valuePositive.isdigit()
        if valuePositive.isdigit():
            if debug:
                print "INT!"
            return int(valueIn)
        try:
            value = float(valueIn)
            if debug:
                print "FLOAT!"
            return value
        except:
             pass

        #if('=' in valueIn or '(' in valueIn or ')' in valueIn):
        #    print "=()!"
        #    return valueIn

        noED = list(set(valueIn)-set('ED 1234567890+-'))
        word = ''.join(noED)
        #print "word=|%s|" %word
        if word.isalpha():
            if debug:
                print "WORD!"
            return valueIn

        v0 = valueIn[0]
        vm = valueIn.find('-',1)
        vp = valueIn.find('+',1)
        #print "valueIn = ",valueIn
        if '-'==v0 or '+'==v0:
            valueLeft = valueIn[1:]
        else:
            valueLeft = valueIn


        #print "valueIn = |%s|" %(valueIn)
        #print "v0 = |%s|" %v0
        if v0=='-':
            vFactor=-1.
        elif v0=='+' or v0.isdigit():
            vFactor=1.
        else:
            raise Exception('huh?_1')
        
        if vm>0:
            sline = valueLeft.split('-')
            expFactor = -1.
        elif vp>0:
            sline = valueLeft.split('+')
            expFactor = 1.
        else:
            self.log().info("*valueIn = %s" %(valueIn))
            raise Exception('huh?_2')

        s0 = vFactor*float(sline[0])
        s1 = expFactor*int(sline[1])

        value = s0*10**(s1)
        #print "valueOut = |%s|" %value
        
        if debug:
            print "SCIENTIFIC!"
        return value
    

def stringParser(stringIn):
    typeCheck = ''
    n=0
    for i,s in enumerate(stringIn):
        if s in "+-":   state='+'
        elif s==" ":      state=' '
        elif s==".":      state='.'
        elif s in "eEdD": state='e'
        elif s.isdigit(): state='1'
        elif s.isalpha() or s in "()*/=]['\"": return 'string' # string character
        else:
            msg = "s=|%r|" %(s)
            raise Exception(msg)
        ###
        
        #print "s=%s stringIn[i-1]=%s" %(state,typeCheck[i-1])
        #print "i=%s s=%s typeCheck=%s" %(i,s,typeCheck)
        if i==0:
            typeCheck += state
            n+=1
        elif typeCheck[n-1]!=state:
            typeCheck += state
            n+=1
        elif state in 'e .+': # double e, space, dot, plus
            return 'string'
        ###
    ###
    if typeCheck==' ':
        return None

    typeCheck = typeCheck.strip()
    if typeCheck in ['1','+1']: # integer
        return int(stringIn)

    elif typeCheck in [ '1.', '1.1', '.1',  # float
                       '+1.','+1.1','+.1']:
        return float(stringIn)

    elif typeCheck in [ '1.1e1', '1.1e+1', '1.e1', '1.e+1', # python scientific
                       '+1.1e1','+1.1e+1','+1.e1','+1.e+1',
                       '.1e1','.1e+1','+.1e1','+.1e+1',]:
        return float(stringIn)

    elif typeCheck in ['1+1','+1+1','.1+1','+.1+1']: # nastran scientific
        stringReversed = stringIn[::-1]
        i = stringReversed.index('+')
        lString = list(stringIn)
        lString.insert(-i-1,'e')
        #print "lString = ",lString
        out = ''.join(lString)
        print "out = ",out
        return float(out)
    else:
        #print "string = ",stringIn
        #print "typeCheck = ",typeCheck
        return 'string'
    
    print "typeCheck = |%s|" %(typeCheck)
    raise Exception('this should never happen...')

if __name__=='__main__':
    print stringParser('123')
    print stringParser('+123')
    print stringParser('.234')
    print stringParser('+.234')
    print stringParser('-.234')
    print stringParser('1+5')
    print "abc = |%s|" %(stringParser('abc'))
    print "eeg = |%s|" %(stringParser('eeg'))
    #print "e1 = |%s|" %(stringParser('\T'))
    print stringParser('.e1')
 