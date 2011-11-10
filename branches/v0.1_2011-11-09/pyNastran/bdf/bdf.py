import os
import sys
import copy
from math import ceil

# 3rd party
import numpy
from numpy import any,cross

# my code
from pyNastran.general.general import ListPrint
#from mathFunctions import *

from cards import * # reads all the card types - GRID, CQUAD4, FORCE, PSHELL, etc.
from BDF_Card        import BDF_Card
from bdf_helper      import getMethods,addMethods,writeMesh,cardMethods,XrefMesh
from caseControlDeck import CaseControlDeck
from fieldWriter     import printCard

class BDF(getMethods,addMethods,writeMesh,cardMethods,XrefMesh):
    modelType = 'nastran'
    isStructured = False
    
    #def setCardsToInclude():
    #    pass

    def __init__(self,debug=True,log=None):
        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            log = loggerObj.startLog('debug') # or info
        self.log = log

        ## allows the BDF variables to be scoped properly (i think...)
        getMethods.__init__(self)
        addMethods.__init__(self)
        writeMesh.__init__(self)
        cardMethods.__init__(self)
        XrefMesh.__init__(self)

        ## useful in debugging errors in input
        self.debug = debug

        self._initSolution()

        ## lines that were rejected b/c they were for a card
        ## that isnt supported
        self.rejects = []
        ## cards that were created, but not processed
        self.rejectCards = []
        ## list of execive control deck lines
        self.executiveControlLines = []
        ## list of case control deck lines
        self.caseControlLines = []
        ## the analysis type

        ## the list of possible cards that will be parsed
        self.cardsToRead = set([
        'PARAM','INCLUDE',  # '='
        'GRID','GRDSET', #'RINGAX',

        'CONM2','CMASS1','CMASS2','CMASS3','CMASS4',
        'CELAS1','CELAS2','CELAS3','CELAS4',
        'CBAR','CROD','CTUBE','CBEAM','CONROD',
        'CTRIA3','CTRIA6','CTRIAX6',
        'CQUAD4','CQUAD8','CQUADR','CQUADX',
        'CHEXA','CPENTA','CTETRA',
        'CSHEAR','CVISC','CRAC2D','CRAC3D',
        'RBAR','RBAR1','RBE1','RBE2','RBE3',

        'PMASS',
        'PELAS',
        'PROD','PBAR','PBEAM',#'PBEAM3','PBEAML'
        'PSHELL','PCOMP', # 'PCOMPG',
        'PSOLID','PLSOLID',
        'MAT1','MAT2','MAT3','MAT4','MAT5','MAT8','MAT9','MAT10',
         #'MATT1','MATT2','MATT3','MATT4','MATT5','MATT8','MATT9',
         #'MATS1',

        'SPC','SPC1','SPCD','SPCADD','SPCAX',
        'MPC','MPCADD',
        'SUPORT','SUPORT1',

        'LOAD',
        'FORCE','FORCE1','FORCE2',
        'PLOAD','PLOAD1','PLOAD2','PLOAD4',
        'MOMENT','MOMENT1','MOMENT2',

        # dynamic
        'DAREA','NLPARM',

        # aero
        'FLFACT','AERO','AEROS','GUST','FLUTTER','GRAV',
        'CAERO1',#'CAERO2','CAERO3','CAERO4','CAERO5',
        'SPLINE1',#'SPLINE2','SPLINE3','SPLINE4','SPLINE5','SPLINE6','SPLINE7',

        #'CORD1R','CORD1C','CORD1S',
        'CORD2R',#'CORD2C','CORD2S',
        'ENDDATA',
        
        'TEMP',#'TEMPD',
        'QBDY1','QBDY2','QBDY3','QHBDY',
        'CHBDYE','CHBDYG','CHBDYP',
        'PCONV','PCONVM','PHBDY',
        'RADBC','CONV',  #'RADM',
        
        # optimization
        'DCONSTR','DESVAR','DDVAL',
        ])
        ## was playing around with an idea...does nothing for now...
        self.cardsToWrite = self.cardsToRead

    def _initSolution(self):
        self.infilename = None
        self.autoReject = False
        self.solmap_toValue = {
                        'SESTSTATIC': 101,
                        'SEMODES'   : 103,
                        'BUCKLING'  : 105,
                        'SEBUCKL'   : 105,
                        'NLSTATIC'  : 106,
                        'SEDCEIG'   : 107,
                        'SEDFREQ'   : 108,
                        'SEDTRAN'   : 109,
                        'SEMCEIG'   : 110,
                        'SEMFREQ'   : 111,
                        'SEMTRAN'   : 112,
                        'CYCSTATX'  : 114,
                        'CYCMODE'   : 115,
                        'CYCBUCKL'  : 116,
                        'CYCFREQ'   : 118,
                        'NLTRAN'    : 129,
                        'AESTAT'    : 144,
                        'FLUTTR'    : 145,
                        'SEAERO'    : 146,
                        'NLSCSH'    : 153,
                        'NLTCSH'    : 159,
                        'DBTRANS'   : 190,
                        'DESOPT'    : 200,
                       }


        self.rsolmap_toStr = {
                         101 : 'SESTSTATIC',
                         103 : 'SEMODES'   ,
                         105 : 'BUCKLING'  ,  # SEBUCKL
                         106 : 'NLSTATIC'  ,
                         107 : 'SEDCEIG'   ,
                         108 : 'SEDFREQ'   ,
                         109 : 'SEDTRAN'   ,
                         110 : 'SEMCEIG'   ,
                         111 : 'SEMFREQ'   ,
                         112 : 'SEMTRAN'   ,
                         114 : 'CYCSTATX'  ,
                         115 : 'CYCMODE'   ,
                         116 : 'CYCBUCKL'  ,
                         118 : 'CYCFREQ'   ,
                         129 : 'NLTRAN'    ,
                         144 : 'AESTAT'    ,
                         145 : 'FLUTTR'    ,
                         146 : 'SEAERO'    ,
                         153 : 'NLSCSH'    ,
                         159 : 'NLTCSH'    ,
                         190 : 'DBTRANS'   ,
                         200 : 'DESOPT'    ,
                       }
        self._initStructuralDefaults()
        self._initThermalDefaults()

    def _initStructuralDefaults(self):
        self.sol = None
        ## used in solution 600
        self.solMethod = None
        self.iSolLine = None
        self.caseControlDeck = CaseControlDeck([],self.log)
        #self.executiveControlLines = [self.sol]

        # main structural block
        ## store the PARAM cards
        self.params = {}
        ## stores SPOINT, GRID cards
        self.nodes = {}
        ## stores GRIDSET card
        self.gridSet = None
        ## stores LOTS of elements
        self.elements = {}
        ## stores LOTS of elements
        self.properties = {}
        ## stores MAT1,MAT2,...,MAT10 materials
        self.materials = {}
        ## stores LOAD,FORCE,MOMENT
        self.loads = {}
        ## stores coordinate systems
        self.coords = {0: CORD2R() }

        # constraints
        ## stores SUPORT1s
        self.constraints = {} # suport1, anything else???
        self.suports = [] # suport

        ## stores SPCADD,SPC,SPC1,SPCD,SPCAX
        self.spcObject = constraintObject()
        ## stores MPCADD,MPC
        self.mpcObject = constraintObject()

        # dynamic cards
        ## stores DAREA
        self.dareas  = {}
        self.nlparms = {}

        # aero cards
        ## stores CAERO1
        self.caeros   = {}  # can this be combined with self.elements???
        ## stores AERO & AEROS
        ## @warning possible bug
        self.aeros    = {}
        ## stores GUST cards
        self.gusts    = {}  # can this be simplified ???
        ## stores FLFACT
        self.flfacts  = {}  # can this be simplified ???
        ## stores FLUTTER
        self.flutters = {}
        ## store SPLINE1
        self.splines  = {} # maybe put into self.elements???
        ## stores GRAV
        self.gravs = {}
        
        # optimization
        self.dconstrs = {}
        self.desvars  = {}
        self.ddvals   = {}

    def _initThermalDefaults(self):
        # BCs
        ## stores thermal boundary conditions - CONV,RADBC
        self.bcs   = {}  # e.g. RADBC
        
        # elements
        # see self.elements

        # properties
        ## stores other thermal properties - unused ???
        #self.thermalProperties    = {}
        ## stores PHBDY
        self.phbdys               = {}
        ## stores convection properties - PCONV, PCONVM ???
        self.convectionProperties = {}

    def _setInfile(self,infilename,includeDir=None):
        """
        sets up the basic file/lines/cardCounting operations
        """
        ## automatically rejects every parsable card (default=False)
        self.autoReject   = False
        ## is the active file done reading
        self.doneReading  = False
        ## was an ENDDATA card found
        self.foundEndData = False

        if includeDir is None:
            includeDir = os.path.dirname(infilename)
        ## the active filename (string)
        self.infilename = infilename
        ## the directory of the 1st BDF (include BDFs are relative to this one)
        self.includeDir = includeDir
        ## list of infile objects (needed for INCLUDE files)
        self.infilesPack     = []
        ## list of lines from self.activeFilename that are stored
        self.linesPack       = []
        ## list of all open filenames
        self.activeFileNames = []
        ## stores the line number of self.activefilename that the parser is on
        ## very helpful when debugging
        self.lineNumbers     = []
        ## dictionary that says whether self.infilename is open/close (boolean0
        self.isOpened = {self.infilename: False}
        ## list of all read in cards - useful in determining if
        ## entire BDF was read & really useful in debugging
        self.cardCount = {}

    def openFile(self,infileName):
        """
        Takes a filename and opens the file.  
        This method is used in order to support INCLUDE files.
        """
        #print self.isOpened
        if self.isOpened[infileName]==False:
            self.activeFileNames.append(infileName)
            #self.log().info("*openFile bdf=|%s|  pwd=|%s|" %(infileName,os.getcwd()))
            assert os.path.exists(infileName),"infileName=|%s| does not exist..." %(infileName)
            infile = open(infileName,'r')
            self.infilesPack.append(infile)
            self.lineNumbers.append(0)
            self.isOpened[infileName]=True
            self.linesPack.append([])
        ###
        else:
            pass
            #print "is already open...skipping"
        ###

    def getFileStats(self):
        filename   = self.activeFileNames[-1]
        return (filename,self.getLineNumber())

    def getLineNumber(self):
        lineNumber = self.lineNumbers[-1]
        return lineNumber

    def getNextLine(self,debug=False):
        self.lineNumbers[-1]+=1
        linesPack = self.makeLinesPack(debug=False)
        #print "len(linesPack) = ",len(linesPack)
        #for line in linesPack:
        #    print "$  |%r|" %(line)

        if len(linesPack)==0:
            self.closeFile()
            return None
            #linesPack = self.makeLinesPack(debug=debug)
            #return lastLine
        #print linesPack[0]
        return linesPack.pop(0)

    #def getNextLine2(self):
    #    self.lineNumbers[-1]+=1
    #    return self.infilesPack[-1].readline()

    #    infile = self.infilesPack[-1]
    #    #print "infile = |%s|" %(infile),type(infile)
    #    line = infile.readline()
    #    print "line = |%r|" %(line)
    #    return line

    def closeFile(self,debug=False):
        """
        Closes the active file object.
        If no files are open, the function is skipped.
        This method is used in order to support INCLUDE files.
        """
        if len(self.infilesPack)==0:
            return
        if debug:
            print "*closing"
        infile = self.infilesPack.pop()
        infile.close()

        if debug:
            print [os.path.relpath(fname) for fname in self.activeFileNames]
        lineNumbers = self.lineNumbers.pop()
        activeFileName = self.activeFileNames.pop()
        linesPack = self.linesPack.pop()
        self.isOpened[activeFileName] = False
        nlines = len(self.linesPack[-1])
        ## determines if self.activefilename should be closed at the next opportunity
        self.doneReading = False
        if debug:
            print "activeFileName=|%s| infilename=%s len(pack)=%s\n" %(os.path.relpath(activeFileName),os.path.relpath(self.infilename),nlines)
        #print "\n\n"

    def read(self,infilename,includeDir=None,xref=True):
        """
        main read method for the bdf
        @param infilename the input bdf
        @param includeDir the relative path to any include files (default=None if no include files)
        @param xref should the bdf be cross referenced (default=True)
        @param log a logger object (default=None -> a simplified one will be created)
        """
        self._setInfile(infilename,includeDir)

        self.log().info('---starting BDF.read of %s---' %(os.path.relpath(self.infilename)))
        sys.stdout.flush()

        #self.debug = True
        if self.debug:
            self.log().info("*BDF.read")
        self.readExecutiveControlDeck()
        self.readCaseControlDeck(self.infilename)
        self.readBulkDataDeck()
        #self.closeFile()
        self.crossReference(xref=xref)
        if self.debug:
            self.log().debug("***BDF.read")
        self.log().info('---finished BDF.read of %s---' %(os.path.relpath(self.infilename)))
        sys.stdout.flush()

        isDone = self.foundEndData
        return ('BulkDataDeck',isDone)

    def isExecutiveControlDeck(self,line):
        """@todo code this..."""
        return True

    def readExecutiveControlDeck(self):
        """
        reads the executive control deck
        """
        self.openFile(self.infilename)
        line = ''
        #self.executiveControlLines = []
        while len(self.activeFileNames)>0: # keep going until finished
            lineIn = self.getNextLine()
            if lineIn==None: # file was closed and a 2nd readCaseControl was called
                return

            line = lineIn.strip()
            if self.debug:
                (n) = self.getLineNumber()
                self.log().debug("line[%s]*= |%r|" %(n,line.upper()))
            self.executiveControlLines.append(lineIn)
            if 'CEND' in line:
                break
            ###
        ###
        self.parseExecutiveControlDeck()

    def parseExecutiveControlDeck(self):
        """
        extracts the solution from the executive control deck
        """
        for i,eline in enumerate(self.executiveControlLines):
            #print 'eLine = |%r|' %(eline)
            uline = eline.strip().upper()
            if 'SOL ' in uline:
                #print "line = ",uline
                if ',' in uline:
                    sline = uline.split(',') # SOL 600,method
                    solValue = sline[0]
                    method = sline[1]

                    #print "sline    = |%s|" %(sline)
                    #print "sline2   = |%s|" %(sline2)
                else:
                    solValue = uline
                    method = None

                #print "solValue = |%s|" %(solValue)
                sol = solValue[3:].strip()
                    
                assert self.sol==None,'cannot overwrite solution existing=|SOL %s| new =|%s|' %(self.sol,sline2)
                self.iSolLine = i

                try:
                    self.updateSolution(sol,method)
                except:
                    #msg = 'updateSolution failed...sline2=%s sline=%s' %(sline2,sline)
                    #raise RuntimeError(msg)
                    raise

            ###
        ###

    def updateSolution(self,sol,method=None):
        """
        updates the overall solution type (e.g. 101,200,600)
        @param self   the object pointer
        @param sol    the solution type (101,103, etc)
        @param method the solution method (only for SOL=600), default=None
        """

        ## the integer of the solution type (e.g. SOL 101)
        try:
            self.sol = int(sol)
            #print "sol = |%s|" %(sol)
        except:
            #print "sol = |%r|" %(sol)
            self.sol = self.solmap_toValue[sol]
            #print "sol = ",self.sol

        if self.sol==600:
            ## solution 600 method modifier
            self.solMethod = method.strip()
            print "sol=%s method=%s" %(self.sol,self.solMethod)
        else: # very common
            self.solMethod = None
        ###
        #print "sol=%s method=%s" %(self.sol,self.solMethod)
        

    def isCaseControlDeck(self,line):
        """@todo not done..."""
        return True

    def readCaseControlDeck(self,infilename):
        #print "opening |%s|" %(infilename)
        self.openFile(infilename)
        #self.log().info("reading Case Control Deck...")
        line = ''
        #self.caseControlControlLines = []

        i = 0
        while len(self.activeFileNames)>0: # keep going until finished
        #while 'BEGIN BULK' not in line:
            lineIn = self.getNextLine()
            #print "lineIn = ",lineIn
            if lineIn==None: # file was closed and a 2nd readCaseControl was called
                return
            line = lineIn.strip().split('$')[0].strip()

            if 'INCLUDE' in line.upper():
                nextLine = self.getNextLine().strip().split('$')[0].strip()
                includeLines = [line]
                #print "^&*1",nextLine
                while '\\' in nextLine or '/' in nextLine: # more includes
                    includeLines.append(nextLine)
                    nextLine = self.getNextLine().strip().split('$')[0].strip()
                    #print "^&*2",nextLine

                    
                #print "include lines = |%s|" %(includeLines)
                filename = self.getIncludeFileName(includeLines,'INCLUDE')

                self.addIncludeFile(filename)
                #self.openFile(filename)
                self.readCaseControlDeck(filename)
                line = nextLine
                #print "appending |%r|" %(nextLine)
                self.caseControlLines.append(nextLine)
            else:
                #print "appending |%r|" %(lineIn)
                self.caseControlLines.append(lineIn)
            ###

            #print "*line = |%s|" %(line)
            if 'BEGIN BULK' in line:
                #print "breaking"
                break
            if i>200:
                raise RuntimeError('there are too many lines in the Case Control Deck < 200')
            i+=1
        #self.log().info("finished with Case Control Deck..")
        #print "self.caseControlLines = ",self.caseControlLines
        
        #for line in self.caseControlLines:
        #    print "** line=|%r|" %(line)

        self.caseControlDeck = CaseControlDeck(self.caseControlLines,self.log)
        #print "done w/ case control..."
        #print '***********************'
        return self.caseControlLines

    def Is(self,card,cardCheck):
        """
        a method for determining the cardName
        """
        #print "card=%s" %(card)
        #return cardCheck in card[0][0:8]
        return any([cardCheck in field[0:8].lstrip().rstrip(' *') for field in card])

    def isPrintable(self,cardName):
        """can the card be printed"""
        #cardName = self.getCardName(card)
        
        if cardName in self.cardsToWrite:
            #print "*card = ",card
            #print "WcardName = |%s|" %(cardName)
            return False
        return True

    def getCardName(self,cardLines):
        """
        Given a list of card lines, determines the cardName.
        @param self      the object pointer
        @param cardLines the list of lines that define the card
        @retval cardName the name of the card
        @note
            Parses the first 8 characters of a card, splits bassed on a comma,
            pulls off any spaces or asterisks and returns what is left.
        """
        #self.log().debug("getting cardName...")
        cardName = cardLines[0][0:8].strip()
        if ',' in cardName:
            cardName = cardName.split(',')[0].strip()

        cardName = cardName.lstrip().rstrip(' *')
        #self.log().debug("getCardName cardName=|%s|" %(cardName))
        return cardName
    
    def isReject(self,cardName):
        """can the card be read"""
        #cardName = self.getCardName(card)
        if cardName.startswith('='):
            return False
        elif cardName in self.cardsToRead:
            #print "*card = ",card
            #print "RcardName = |%s|" %(cardName)
            return False
        if cardName.strip():
            self.log().info("RejectCardName = |%s|" %(cardName))
        return True

    def getIncludeFileName(self,cardLines,cardName):
        """parses an INCLUDE file split into multiple lines (as a list)
        @param self the object poitner
        @param cardLines the list of lines in the include card (all the lines!)
        @param cardName
            INCLUDE or include (needed to strip it off without converting the case
        """
        cardLines2 = []
        for line in cardLines:
            line = line.strip('\t\r\n ')
            cardLines2.append(line)

        #print "cardLinesRaw = ",cardLines2
        cardLines2[0] = cardLines2[0][7:].strip() # truncate the cardname
        filename = ''.join(cardLines2)
        filename = filename.strip('"').strip("'")
        #print 'filename = |%s|' %(filename)
        filename = os.path.join(self.includeDir,filename)
        return filename

    def addIncludeFile(self,infileName):
        """
        This method must be called before opening an INCLUDE file.
        Identifies the new file as being opened.
        """
        self.isOpened[infileName] = False

    def readBulkDataDeck(self):
        """parses the Bulk Data Deck"""
        debug = self.debug
        #debug = False
        
        if self.debug:
            self.log().debug("*readBulkDataDeck")
        self.openFile(self.infilename)
        #self.nodes = {}
        #self.elements = {}
        #self.rejects = []
        
        #oldCardObj = BDF_Card()
        while len(self.activeFileNames)>0: # keep going until finished
            (rawCard,card,cardName) = self.getCard(debug=False) # gets the cardLines
            #print "outcard = ",card
            #if cardName=='CQUAD4':
            #    print "card = ",card
            
            if cardName=='INCLUDE':
                #print "rawCard = ",rawCard
                #print "card    = ",card
                filename = self.getIncludeFileName(rawCard,cardName)
                #print 'filename = ',os.path.relpath(filename)
                self.addIncludeFile(filename)
                self.openFile(filename)
                reject = '$ INCLUDE processed:  %s\n' %(filename)
                self.rejects.append([reject])
                continue

            if not self.isReject(cardName):
                #print ""
                #print "not a reject"
                card = self.processCard(card) # parse the card into fields
                #print "processedCard = ",card
            elif card[0].strip()=='':
                #print "funny strip thing..."
                pass
            else:
                #print "reject!"
                self.rejects.append(card)
                continue
                #print " rejecting card = ",card
                #card = self.processCard(card)
                #sys.exit()
            

            #print "card2 = ",ListPrint(card)
            #print "card = ",card
            cardName = self.getCardName(card)
            
            if 'ENDDATA' in cardName:
                #print cardName
                break # exits while loop
            #self.log().debug('cardName = |%s|' %(cardName))
            
            #cardObj = BDF_Card(card,oldCardObj)
            cardObj = BDF_Card(card)

            nCards = 1
            #special = False
            if '=' in cardName:
                nCards = cardName.strip('=()')
                if nCards:
                    nCards = int(nCards)
                else:
                    nCards = 1
                    #special = True
                #print "nCards = ",nCards
                #cardName = oldCardObj.field(0)
            ###

            for iCard in range(nCards):
                #print "----------------------------"
                #if special:
                #    print "iCard = ",iCard
                self.addCard(card,cardName,iCard=0,oldCardObj=None)
                #if self.foundEndData:
                #    break
            ### iCard
            if self.doneReading or len(self.linesPack[-1])==0:
                #print "doneReading=%s len(pack)=%s" %(self.doneReading,len(self.linesPack[-1]))
                self.closeFile()
            ###
            #oldCardObj = copy.deepcopy(cardObj) # used for =(*1) stuff
            #print ""
            
            #print "self.linesPack[-1] = ",len(self.linesPack[-1])
            #print "self.activeFileNames = ",self.activeFileNames
        ### end of while loop

        #self.debug = True
        if self.debug:
            #for nid,node in self.nodes.items():
            #    print node
            #for eid,element in self.elements.items():
            #    print element
            
            self.log().debug("\n$REJECTS")
            #for reject in self.rejects:
                #print printCard(reject)
                #print ''.join(reject)
            self.log().debug("***readBulkDataDeck")
        ###

    def addCard(self,card,cardName,iCard=0,oldCardObj=None):
        """
        adds a card object to the BDF object. 
        @param self the object pointer
        @param card the list of the card fields -> ['GRID',1,2,]
        @param cardName the cardName -> 'GRID'
        @param iCard used when reading Nastran Free-Format (disabled)
        @param oldCardObj the last card object that was returned (type=BDF_Card or None; default=None)
        @retval cardObject the card object representation of card
        @note
            this is a very useful method for interfacing with the code
        @note
            the cardObject is not a card-type object...so not a GRID card
            or CQUAD4 object.  It's a BDF_Card Object.  However, you know the type (assuming a GRID),
            so just call the mesh.Node(nid) to get the Node object that was just created.
        @warning cardObject is not returned
        """
        #if cardName != 'CQUAD4':
        #    print cardName
        #print "card = ",card
        cardObj = BDF_Card(card,oldCardObj=None)
        if self.debug:
            self.log().debug("*oldCardObj = \n%s" %(oldCardObj))
            self.log().debug("*cardObj = \n%s" %(cardObj))
        #cardObj.applyOldFields(iCard)

        try:
            if self.autoReject==True:
                print 'rejecting processed %s' %(card)
                self.rejectCards.append(card)
            elif card==[] or cardName=='':
                pass
            elif cardName=='PARAM':
                param = PARAM(cardObj)
                self.addParam(param)
            elif cardName=='GRDSET':
                self.gridSet = GRDSET(cardObj)
            elif cardName=='GRID':
                node = GRID(cardObj)
                self.addNode(node)
            #elif cardName=='RINGAX':
            #    node = RINGAX(cardObj)
            #    self.addNode(node)
            #elif cardName=='SPOINT':  # not done
            #    node = SPOINT(cardObj)
            #    self.addNode(node)

            elif cardName=='CQUAD4':
                elem = CQUAD4(cardObj)
                self.addElement(elem)
            elif cardName=='CQUAD8':
                elem = CQUAD8(cardObj)
                self.addElement(elem)

            elif cardName=='CTRIA3':
                elem = CTRIA3(cardObj)
                self.addElement(elem)
            elif cardName=='CTRIA6':
                elem = CTRIA6(cardObj)
                self.addElement(elem)

            elif cardName=='CTRIAX6':
                elem = CTRIAX6(cardObj)
                self.addElement(elem)
            elif cardName=='CQUAD':
                elem = CQUAD(cardObj)
                self.addElement(elem)
            elif cardName=='CQUADR':
                elem = CQUADR(cardObj)
                self.addElement(elem)
            elif cardName=='CQUADX':
                elem = CQUADX(cardObj)
                self.addElement(elem)
            elif cardName=='CVISC':
                elem = CVISC(cardObj)
                self.addElement(elem)
            elif cardName=='CSHEAR':
                elem = CSHEAR(cardObj)
                self.addElement(elem)

            elif cardName=='CRAC2D':
                elem = CRAC2D(cardObj)
                self.addElement(elem)
            elif cardName=='CRAC3D':
                elem = CRAC3D(cardObj)
                self.addElement(elem)



            elif cardName=='CTETRA':
                nFields = cardObj.nFields()
                if   nFields==7:    elem = CTETRA4(cardObj) # 4+3
                else:               elem = CTETRA10(cardObj)# 10+3
                #elif nFields==13:   elem = CTETRA10(cardObj)# 10+3
                #else: raise Exception('invalid number of CTETRA nodes=%s card=%s' %(nFields-3,str(cardObj)))
                self.addElement(elem)
            elif cardName=='CHEXA':
                nFields = cardObj.nFields()
                if   nFields==11: elem = CHEXA8(cardObj)  # 8+3
                else:             elem = CHEXA20(cardObj) # 20+3
                #elif nFields==23: elem = CHEXA20(cardObj) # 20+3
                #else: raise Exception('invalid number of CPENTA nodes=%s card=%s' %(nFields-3,str(cardObj)))
                self.addElement(elem)
            elif cardName=='CPENTA': # 6/15
                nFields = cardObj.nFields()
                if   nFields==9:  elem = CPENTA6(cardObj)  # 6+3
                else:             elem = CPENTA15(cardObj) # 15+3
                #elif nFields==18: elem = CPENTA15(cardObj) # 15+3
                #else: raise Exception('invalid number of CPENTA nodes=%s card=%s' %(nFields-3,str(cardObj)))
                self.addElement(elem)

            elif cardName=='CBAR':
                elem = CBAR(cardObj)
                self.addElement(elem)
            elif cardName=='CBEAM':
                elem = CBEAM(cardObj)
                self.addElement(elem)
            elif cardName=='CROD':
                elem = CROD(cardObj)
                self.addElement(elem)
            elif cardName=='CONROD':
                elem = CONROD(cardObj)
                self.addElement(elem)
            elif cardName=='CTUBE':
                elem = CTUBE(cardObj)
                self.addElement(elem)

            elif cardName=='CELAS1':
                elem = CELAS1(cardObj)
                self.addElement(elem)
            elif cardName=='CELAS2':
                (elem) = CELAS2(cardObj)
                self.addElement(elem)
            elif cardName=='CELAS3':
                (elem) = CELAS3(cardObj)
                self.addElement(elem)
            elif cardName=='CELAS4':
                (elem) = CELAS4(cardObj)
                self.addElement(elem)
            elif cardName=='CONM2': # not done...
                elem = CONM2(cardObj)
                self.addElement(elem)

            elif cardName=='CMASS1':
                elem = CMASS1(cardObj)
                self.addElement(elem)
            elif cardName=='CMASS2':
                elem = CMASS2(cardObj)
                self.addElement(elem)
            elif cardName=='CMASS3':
                elem = CMASS3(cardObj)
                self.addElement(elem)
            elif cardName=='CMASS4':
                elem = CMASS4(cardObj)
                self.addElement(elem)


            elif cardName=='RBAR':
                (elem) = RBAR(cardObj)
                self.addElement(elem)
            elif cardName=='RBAR1':
                (elem) = RBAR1(cardObj)
                self.addElement(elem)

            elif cardName=='RBE1':
                (elem) = RBE1(cardObj)
                self.addElement(elem)
            elif cardName=='RBE2':
                (elem) = RBE2(cardObj)
                self.addElement(elem)
            elif cardName=='RBE3':
                (elem) = RBE3(cardObj)
                self.addElement(elem)

            elif cardName=='PELAS':
                prop = PELAS(cardObj)
                if cardObj.field(5):
                    prop = PELAS(cardObj,1) # makes 2nd PELAS card
                self.addProperty(prop)

            elif cardName=='PBAR':
                prop = PBAR(cardObj)
                self.addProperty(prop)
            #elif cardName=='PBARL':
            #    prop = PBARL(cardObj)
            #    self.addProperty(prop)
            elif cardName=='PBEAM':
                prop = PBEAM(cardObj)
                self.addProperty(prop)
            #elif cardName=='PBEAM3':
            #    prop = PBEAM3(cardObj)
            #    self.addProperty(prop)
            #elif cardName=='PBEAML':
            #    prop = PBEAML(cardObj)
            #    self.addProperty(prop)
            elif cardName=='PROD':
                prop = PROD(cardObj)
                self.addProperty(prop)
            elif cardName=='PTUBE':
                prop = PTUBE(cardObj)
                self.addProperty(prop)
            elif cardName=='PMASS':
                prop = PMASS(cardObj,nOffset=0)
                self.addProperty(prop)

                if card.field(3):
                    prop = PMASS(cardObj,nOffset=1)
                    self.addProperty(prop)
                if card.field(5):
                    prop = PMASS(cardObj,nOffset=2)
                    self.addProperty(prop)
                if card.field(7):
                    prop = PMASS(cardObj,nOffset=3)
                    self.addProperty(prop)
                ###
            elif cardName=='PSHELL':
                prop = PSHELL(cardObj)
                self.addProperty(prop)
            elif cardName=='PCOMP':
                prop = PCOMP(cardObj)
                self.addProperty(prop)
            #elif cardName=='PCOMPG':
            #    prop = PCOMPG(cardObj)
            #    self.addProperty(prop)

            elif cardName=='PSOLID':
                prop = PSOLID(cardObj)
                self.addProperty(prop)
            elif cardName=='PLSOLID':
                prop = PLSOLID(cardObj)
                self.addProperty(prop)

            elif cardName=='MAT1':
                material = MAT1(cardObj)
                self.addMaterial(material)
            elif cardName=='MAT2':
                material = MAT2(cardObj)
                self.addMaterial(material)
            elif cardName=='MAT3':
                material = MAT3(cardObj)
                self.addMaterial(material)
            elif cardName=='MAT4':
                material = MAT4(cardObj)
                self.addMaterial(material) # maybe addThermalMaterial
            elif cardName=='MAT5':
                material = MAT5(cardObj)
                self.addMaterial(material) # maybe addThermalMaterial
            elif cardName=='MAT8':  # note there is no MAT6 or MAT7
                material = MAT8(cardObj)
                self.addMaterial(material)
            elif cardName=='MAT9':
                material = MAT9(cardObj)
                self.addMaterial(material)
            elif cardName=='MAT10':
                material = MAT10(cardObj)
                self.addMaterial(material)

            #elif cardName=='MATS1':
            #    material = MATS1(cardObj)
            #    self.addStressMaterial(material)
            #elif cardName=='MATT1':
            #    material = MATT1(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName=='MATT2':
            #    material = MATT2(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName=='MATT3':
            #    material = MATT3(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName=='MATT4':
            #    material = MATT4(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName=='MATT5':
            #    material = MATT5(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName=='MATT8':
            #    material = MATT8(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName=='MATT9':
            #    material = MATT9(cardObj)
            #    self.addTempMaterial(material)

            elif cardName=='FORCE':
                force = FORCE(cardObj)
                self.addLoad(force)
            elif cardName=='FORCE1':
                force = FORCE1(cardObj)
                self.addLoad(force)
            elif cardName=='FORCE2':
                force = FORCE2(cardObj)
                self.addLoad(force)
            elif cardName=='MOMENT':
                moment = MOMENT(cardObj)
                self.addLoad(moment)
            elif cardName=='MOMENT1':
                moment = MOMENT1(cardObj)
                self.addLoad(moment)
            elif cardName=='MOMENT2':
                moment = MOMENT2(cardObj)
                self.addLoad(moment)
            elif cardName=='LOAD':
                load = LOAD(cardObj)
                self.addLoad(load)
            elif cardName=='PLOAD':
                load = PLOAD(cardObj)
                self.addLoad(load)
            elif cardName=='PLOAD1':
                load = PLOAD1(cardObj)
                self.addLoad(load)
            elif cardName=='PLOAD2':
                load = PLOAD2(cardObj)
                self.addLoad(load)
            elif cardName=='PLOAD4':
                load = PLOAD4(cardObj)
                self.addLoad(load)

            # thermal loads
            elif cardName=='TEMP':
                load = TEMP(cardObj)
                self.addThermalLoad(load)
            #elif cardName=='TEMPD':
            #    load = TEMPD(cardObj)
            #    self.addThermalLoad(load)
            elif cardName=='QBDY1':
                load = QBDY1(cardObj)
                self.addThermalLoad(load)
            elif cardName=='QBDY2':
                load = QBDY2(cardObj)
                self.addThermalLoad(load)
            elif cardName=='QBDY3':
                load = QBDY3(cardObj)
                self.addThermalLoad(load)
            elif cardName=='QHBDY':
                load = QHBDY(cardObj)
                self.addThermalLoad(load)

            # thermal elements
            elif cardName=='CHBDYE':
                element = CHBDYE(cardObj)
                self.addThermalElement(element)
            elif cardName=='CHBDYG':
                element = CHBDYG(cardObj)
                self.addThermalElement(element)
            elif cardName=='CHBDYP':
                element = CHBDYP(cardObj)
                self.addThermalElement(element)

            # thermal properties
            elif cardName=='PCONV':
                prop = PCONV(cardObj)
                self.addConvectionProperty(prop)
            elif cardName=='PCONVM':
                prop = PCONVM(cardObj)
                self.addConvectionProperty(prop)
            elif cardName=='PHBDY':
                prop = PHBDY(cardObj)
                self.addPHBDY(prop)

            # thermal BCs
            elif cardName=='CONV':
                bc = CONV(cardObj)
                self.addThermalBC(bc,bc.eid)
            #elif cardName=='RADM':
            #    bc = RADM(cardObj)
            #    self.addThermalBC(bc,bc.nodamb)
            elif cardName=='RADBC':
                bc = RADBC(cardObj)
                self.addThermalBC(bc,bc.nodamb)

            #elif cardName=='TABLEH1':
            #    load = TABLEH1(cardObj)
            #    self.addTable(load)

            # constraints
            elif cardName=='MPC':
                constraint = MPC(cardObj)
                self.addConstraint_MPC(constraint)
            elif cardName=='MPCADD':
                constraint = MPCADD(cardObj)
                assert not isinstance(constraint,SPCADD)
                self.addConstraint_MPCADD(constraint)

            # constraints
            elif cardName=='SPC':
                constraint = SPC(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPC1':
                constraint = SPC1(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPCAX':
                constraint = SPCAX(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPCD':
                constraint = SPCD(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPCADD':
                constraint = SPCADD(cardObj)
                assert not isinstance(constraint,MPCADD)
                self.addConstraint_SPCADD(constraint)
            elif cardName=='SUPORT':
                suport = SUPORT(cardObj)
                self.addSUPORT(suport)
            elif cardName=='SUPORT1':
                suport1 = SUPORT1(cardObj)
                self.addConstraint(suport1)

            # dynamic
            elif cardName=='DAREA':
                darea = DAREA(cardObj)
                self.addDArea(darea)
            elif cardName=='NLPARM':
                nlparmObj = NLPARM(cardObj)
                self.addNLParm(nlparmObj)

            # aero
            elif cardName=='SPLINE1':
                aero = SPLINE1(cardObj)
                self.addSpline(aero)
            elif cardName=='CAERO1':
                aero = CAERO1(cardObj)
                self.addCAero(aero)
            #elif cardName=='CAERO2':
            #    aero = CAERO2(cardObj)
            #    self.addCAero(aero)
            elif cardName=='AERO':
                aero = AERO(cardObj)
                self.addAero(aero)
            elif cardName=='AEROS':
                aeros = AEROS(cardObj)
                self.addAero(aeros)
            elif cardName=='FLFACT':
                flfact = FLFACT(cardObj)
                self.addFLFACT(flfact)
            elif cardName=='GUST':
                gust = GUST(cardObj)
                self.addGust(gust)
            elif cardName=='GRAV':
                grav = GRAV(cardObj)
                self.addGrav(grav)
            elif cardName=='FLUTTER':
                flutter = FLUTTER(cardObj)
                self.addFlutter(flutter)

            elif cardName=='DCONSTR':
                flutter = DCONSTR(cardObj)
                self.addDConstr(flutter)
            elif cardName=='DESVAR':
                desvar = DESVAR(cardObj)
                self.addDesvar(desvar)
            elif cardName=='DDVAL':
                ddval = DDVAL(cardObj)
                self.addDDVal(ddval)


            # coordinate systems
            elif cardName=='CORD2R':
                coord = CORD2R(cardObj)
                self.addCoord(coord)
            #elif cardName=='CORD2C':
            #    coord = CORD2C(cardObj)
            #    self.addCoord(coord)
            #elif cardName=='CORD2S':
            #    coord = CORD2S(cardObj)
            #    self.addCoord(coord)

            #elif cardName=='CORD1R':
            #    coord = CORD1R(cardObj)
            #    self.addCoord(coord)
            #    if cardObj.field(5):
            #        coord = CORD1R(cardObj,nCoord=1)
            #        self.addCoord(coord)
            #    ###
            #elif cardName=='CORD1C':
            #    coord = CORD1C(cardObj)
            #    self.addCoord(coord)
            #    if cardObj.field(5):
            #        coord = CORD1C(cardObj,nCoord=1)
            #        self.addCoord(coord)
            #    ###
            #elif cardName=='CORD1S':
            #    coord = CORD1S(cardObj)
            #    self.addCoord(coord)
            #    if cardObj.field(5):
            #        coord = CORD1S(cardObj,nCoord=1)
            #        self.addCoord(coord)
            #    ###
            #elif cardName=='CORD3G':
            #    coord = CORD3G(cardObj)
            #    self.addCoord(coord)


            elif 'ENDDATA' in cardName:
                self.foundEndData = True
                #break
            else:
                if '=' not in card[0]:  ## warning cards with = signs in them are not announced when they are rejected
                    print 'rejecting processed %s' %(card)
                self.rejectCards.append(card)
            ###
        except:
            print "cardName = |%r|" %(cardName)
            print "failed! Unreduced Card=%s\n" %(ListPrint(card))
            print "filename = %s\n" %(self.infilename)
            raise
        ### try-except block

        return cardObj
    
