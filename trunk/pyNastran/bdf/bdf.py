import os
import sys
import copy
from math import ceil
from pyNastran.general.general import ListPrint

# 3rd party
import numpy
from numpy import any,cross

# my code
from fieldWriter import printCard
from cards import * # reads all the card types - GRID, CQUAD4, FORCE, PSHELL, etc.
#from mathFunctions import *

from BDF_Card import BDF_Card
from bdf_helper import getMethods,addMethods,writeMesh,cardMethods
from caseControlDeck import CaseControlDeck

class BDF(getMethods,addMethods,writeMesh,cardMethods):
    modelType = 'nastran'
    isStructured = False
    
    def setCardsToInclude():
        pass

    def __init__(self,infilename,log=None):
        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            log = loggerObj.startLog('debug') # or info

        self.debug = False
        self.log = log
        self.infilename = infilename
        self.isOpened=False
        #self.n = 0
        #self.nCards = 0
        self.doneReading = False
        self.foundEndData = False

        self.params = {}

        self.nodes = {}
        self.gridSet = None

        self.elements = {}
        self.properties = {}
        self.materials = {}
        self.loads = {}
        self.coords = {0: CORD2R() }

        self.constraints = {} # suport1, anything else???
        self.spcObject = constraintObject()
        self.mpcObject = constraintObject()

        # aero cards
        self.aeros   = {}
        self.gusts   = {}  # can this be simplified ???
        self.flfacts = {}  # can this be simplified ???
        self.flutters = {}

        self.rejects = []
        self.rejectCards = []
        self.executiveControlLines = []
        self.caseControlLines = []

        self.cardsToRead = set([
        'PARAM','=',
        'GRID','GRDSET',
        
        'CONM2',
        'CELAS1','CELAS2',
        'CBAR','CROD','CTUBE','CBEAM',
        'CTRIA3','CQUAD4',
        'CHEXA','CPENTA','CTETRA',
        'RBAR','RBAR1','RBE1','RBE2','RBE3',
        
        'PELAS',
        'PROD','PBAR','PBEAM',#'PBEAM3','PBEAML'
        'PSHELL','PCOMP', # 'PCOMPG',
        'PSOLID','PLSOLID',
        'MAT1','MAT2','MAT3','MAT4','MAT5','MAT8','MAT9','MAT10',
         #'MATT1','MATT2','MATT3','MATT4','MATT5','MATT8','MATT9',
         #'MATS1',

        'SPC','SPC1','SPCD','SPCADD','SUPORT1',
        'MPC','MPCADD',

        'LOAD',
        'FORCE',#'FORCE1','FORCE2',
        'PLOAD','PLOAD2','PLOAD4',#'PLOAD1',
        #'MOMENT','MOMENT1','MOMENT2',

        'FLFACT','AERO','AEROS','GUST','FLUTTER',
        #'CAERO1','CAERO2','CAERO3','CAERO4','CAERO5',
        #'SPLINE1','SPLINE2','SPLINE3','SPLINE4','SPLINE5','SPLINE6','SPLINE7',
        #'NLPARM',

        'CORD1R','CORD1C','CORD1S',
        'CORD2R','CORD2C','CORD2S',
        'ENDDATA',
        ])
        self.cardsToWrite = self.cardsToRead

    def openFile(self):
        if self.isOpened==False:
            self.log().info("*FEM_Mesh bdf=|%s|  pwd=|%s|" %(self.infilename,os.getcwd()))
            self.infile = open(self.infilename,'r')
            self.isOpened=True
            self.lines = []

    def closeFile(self):
        self.infile.close()

    def read(self,debug=False):
        self.log().info('---starting FEM_Mesh.read of %s---' %(self.infilename))
        sys.stdout.flush()
        self.debug = debug
        if self.debug:
            self.log().info("*FEM_Mesh.read")
        self.readExecutiveControlDeck()
        self.readCaseControlDeck()
        self.readBulkDataDeck()
        self.crossReference()
        self.closeFile()
        if self.debug:
            self.log().debug("***FEM_Mesh.read")
        self.log().info('---finished FEM_Mesh.read of %s---' %(self.infilename))
        sys.stdout.flush()

        isDone = self.foundEndData
        return ('BulkDataDeck',isDone)

    def crossReference(self):
        #print "cross Reference is a temp function"
        #for key,e in self.elements.items():
        #    print(e)
        
        #self.spcObject.crossReference(self)
        #self.caseControlDeck.crossReference(self)
        self.crossReference_Nodes()
        #self.crossReference_Elements()
        #self.crossReference_Properties()

    def crossReference_Nodes(self):
        gridSet = self.gridSet
        for nid,n in self.nodes.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            n.crossReference(self,gridSet)
        ###

    def crossReference_Elements(self):
        for eid,e in self.elements.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            e.crossReference(self)
        ###

    def crossReference_Properties(self):
        for pid,p in self.properties.items():
            #print "n.cid = ",n.cid
            #coord = self.Coord(n.cid)
            #print "*",str(coord)
            p.crossReference(self)
        ###

    def readExecutiveControlDeck(self):
        self.openFile()
        line = ''
        #self.executiveControlLines = []
        while 'CEND' not in line:
            lineIn = self.infile.readline()
            line = lineIn.strip()
            self.executiveControlLines.append(lineIn)
        return self.executiveControlLines

    def readCaseControlDeck(self):
        self.openFile()
        self.log().info("reading Case Control Deck...")
        line = ''
        #self.caseControlControlLines = []
        while 'BEGIN BULK' not in line:
            lineIn = self.infile.readline()
            line = lineIn.strip().split('$')[0].strip()
            #print "*line = |%s|" %(line)
            self.caseControlLines.append(lineIn)
        self.log().info("finished with Case Control Deck..")
        #print "self.caseControlLines = ",self.caseControlLines
        
        self.caseControlDeck = CaseControlDeck(self.caseControlLines,self.log)
        #print "done w/ case control..."
        return self.caseControlLines

    def Is(self,card,cardCheck):
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

    def getCardName(self,card):
        #self.log().debug("getting cardName...")
        cardName = card[0][0:8].strip()
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
            print "RcardName = |%s|" %(cardName)
        return True

    def readBulkDataDeck(self):
        debug = self.debug
        #debug = False
        
        if self.debug:
            self.log().debug("*readBulkDataDeck")
        self.openFile()
        #self.nodes = {}
        #self.elements = {}
        #self.rejects = []
        
        #oldCardObj = BDF_Card()
        while 1: # keep going until finished
            (card,cardName) = self.getCard(debug=debug) # gets the cardLines
            #print "outcard = ",card
            #if cardName=='CQUAD4':
            #    print "card = ",card

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
                print cardName
                break
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
                cardName = oldCardObj.field(0)
            ###

            for iCard in range(nCards):
                #print "----------------------------"
                #if special:
                #    print "iCard = ",iCard
                self.addCard(card,cardName,iCard=0,oldCardObj=None)
                #if self.foundEndData:
                #    break
            ### iCard
            if self.doneReading or len(self.lines)==0:
                break
            ###
            #oldCardObj = copy.deepcopy(cardObj) # used for =(*1) stuff
            #print ""
        
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
    
    def addCard(self,card,cardName,iCard=0,oldCardObj=None):
        #if cardName != 'CQUAD4':
        #    print cardName
        if self.debug:
            print "*oldCardObj = \n",oldCardObj
            print "*cardObj = \n",cardObj
        cardObj = BDF_Card(card,oldCardObj=None)
        #cardObj.applyOldFields(iCard)
        try:
            if card==[] or cardName=='':
                pass
            elif cardName=='PARAM':
                param = PARAM(cardObj)
                self.addParam(param)
            elif cardName=='GRDSET':
                self.gridSet = GRDSET(cardObj)
            elif cardName=='GRID':
                node = GRID(cardObj)
                self.addNode(node)
            #elif cardName=='SPOINT':
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
                #print str(elem).strip()
            elif cardName=='CTUBE':
                elem = CBAR(cardObj)
                self.addElement(elem)

            elif cardName=='CELAS1':
                elem = CELAS1(cardObj)
                self.addElement(elem)
            elif cardName=='CELAS2':
                (elem) = CELAS2(cardObj)  # removed prop from outputs...
                self.addElement(elem)
                #self.addProperty(prop)
            elif cardName=='CONM2': # not done...
                elem = CONM2(cardObj)
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
                material = MAT9(cardObj)
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
                #print "fcard = ",card 
                force = FORCE(cardObj)
                self.addLoad(force)
            elif cardName=='LOAD':
                #print "fcard = ",card 
                load = LOAD(cardObj)
                self.addLoad(load)

            elif cardName=='SPC':
                constraint = SPC(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPC1':
                constraint = SPC1(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPCD':
                constraint = SPC1(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SPCADD':
                constraint = SPCADD(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName=='SUPORT1':
                #print "card = ",card
                constraint = SUPORT1(cardObj)
                self.addConstraint(constraint)
                #print "constraint = ",constraint

            elif cardName=='AERO':
                #print "card = ",card
                aero = AERO(cardObj)
                self.addAero(aero)
            elif cardName=='AEROS':
                #print "card = ",card
                aeros = AEROS(cardObj)
                self.addAero(aeros)
            elif cardName=='FLFACT':
                #print "card = ",card
                flfact = FLFACT(cardObj)
                self.addFLFACT(flfact)
            elif cardName=='GUST':
                #print "card = ",card
                gust = GUST(cardObj)
                self.addGust(gust)
            elif cardName=='FLUTTER':
                #print "card = ",card
                flutter = FLUTTER(cardObj)
                self.addFlutter(flutter)

            elif cardName=='CORD2R':
                coord = CORD2R(cardObj)
                self.addCoord(coord)
            #elif cardName=='CORD2C':
            #    coord = CORD2C(cardObj)
            #    self.addCoord(coord)
            #elif cardName=='CORD2S':
            #    coord = CORD2S(cardObj)
            #    self.addCoord(coord)
            #elif 'CORD' in cardName:
            #    raise Exception('unhandled coordinate system...cardName=%s' %(cardName))
            elif 'ENDDATA' in cardName:
                self.foundEndData = True
                #break
            else:
                print 'rejecting processed %s' %(card)
                self.rejectCards.append(card)
            ###
        except:
            print "failed! Unreduced Card=%s" %(ListPrint(card))
            raise
        ### try-except block

    def sumForces(self):
        for key,loadCase in self.loads.items():
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,FORCE):
                    f = load.mag*load.xyz
                    print "f = ",f
                    F += f
                ###
            self.log().info("case=%s F=%s\n\n" %(key,F))
        ###

    def sumMoments(self):
        p = array([0.,0.5,0.])
        for key,loadCase in self.loads.items():
            M = array([0.,0.,0.])
            F = array([0.,0.,0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load,FORCE):
                    f = load.mag*load.xyz
                    node = self.Node(load.node)
                    #print "node = ",node
                    r = node.Position() - p
                    m = cross(r,f)
                    #print "m    = ",m
                    M += m
                    F += f
                ###
            print "case=%s F=%s M=%s\n\n" %(key,F,M)
        ###

### FEM_Mesh

if __name__=='__main__':
    import sys
    #basepath = os.getcwd()
    #configpath = os.path.join(basepath,'inputs')
    #workpath   = os.path.join(basepath,'outputs')
    #os.chdir(workpath)
    

    #bdfModel   = os.path.join(configpath,'fem.bdf.txt')
    #bdfModel   = os.path.join(configpath,'aeroModel.bdf')
    #bdfModel   = os.path.join('aeroModel_mod.bdf')
    bdfModel   = os.path.join('aeroModel_2.bdf')
    #bdfModel   = os.path.join('hard.bdf')
    #bdfModel   = os.path.join(configpath,'aeroModel_Loads.bdf')
    #bdfModel   = os.path.join(configpath,'test_mesh.bdf')
    #bdfModel   = os.path.join(configpath,'test_tet10.bdf')
    assert os.path.exists(bdfModel),'|%s| doesnt exist' %(bdfModel)
    fem = BDF(bdfModel,log=None)
    fem.read()
    #fem.sumForces()
    #fem.sumMoments()
    
    #print "----------"
    #deck = fem.caseControlDeck #.subcases[1]
    #print "deck = \n",deck
    #print str()
    fem.write('fem.out.bdf')
    #fem.writeAsCTRIA3('fem.out.bdf')


