import os
import sys
import traceback

from pyNastran.bdf.bdf import BDF

import pyNastran.bdf.test
testPath = pyNastran.bdf.test.__path__[0]
print "testPath = ",testPath

def runBDF(folder,bdfFilename,debug=False):
    bdfModel = os.path.join(testPath,folder,bdfFilename)
    assert os.path.exists(bdfModel),'|%s| doesnt exist' %(bdfModel)

    fem1 = BDF(bdfModel,log=None,debug=debug)
    try:
        fem1.read()
        #fem1.sumForces()
        #fem1.sumMoments()
        outModel = bdfModel+'_out'
        fem1.writeAsPatran(outModel)
        #fem1.writeAsCTRIA3(outModel)

        fem2 = BDF(outModel,log=None,debug=debug)
        fem2.read(debug=debug)
        #fem2.sumForces()
        #fem2.sumMoments()
        outModel2 = bdfModel+'_out2'
        fem2.writeAsPatran(outModel2)
        #fem2.writeAsCTRIA3(outModel2)
        compare(fem1,fem2)

    #except KeyboardInterrupt:
    #    sys.exit()
    except:
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #print "\n"
        traceback.print_exc(file=sys.stdout)
        #print msg
        print "-"*80
        raise
    ###
    print "-"*80

def divide(value1,value2):
    if value1==value2:  # good for 0/0
        return 1.0
    else:
        try:
            v = value1/float(value2)
        except:
            v = 0.
        ###
    ###
    return v

def compareCardCount(fem1,fem2):
    cards1 = fem1.cardCount
    cards2 = fem2.cardCount
    computeInts(cards1,cards2)

def computeInts(cards1,cards2):
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys  = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))
    
    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    msg = 'diffKeys1=%s diffKeys2=%s' %(diffKeys2,diffKeys2)
    print msg
    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1: value1 = cards1[key]
        else:                value1 = 0

        if key in listKeys2: value2 = cards2[key]
        else:                value2 = 0
        
        diff = abs(value1-value2)
        star = ' '
        if diff:
            star = '*'

        factor1 = divide(value1,value2)
        factor2 = divide(value2,value1)
        factorMsg = ''
        if factor1 != factor2:
            factorMsg = 'diff=%s factor1=%g factor2=%g' %(diff,factor1,factor2)

        msg += '  %skey=%-7s value1=%-4s value2=%-4s' %(star,key,value1,value2)+factorMsg #+'\n'
        print msg
    ###

def compute(cards1,cards2):
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys  = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))
    
    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    msg = 'diffKeys1=%s diffKeys2=%s' %(diffKeys2,diffKeys2)
    #print msg
    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1: value1 = cards1[key]
        else:                value1 = 0

        if key in listKeys2: value2 = cards2[key]
        else:                value2 = 0
        msg += '   *key=%-7s value1=%-7s value2=%-7s' %(key,value1,value2) #+'\n'
        print msg
    ###

from pyNastran.bdf.bdf import ShellElement, SolidElement,LineElement

def getElementStats(fem1,fem2):
    for key,e in sorted(fem1.elements.items()):
        if isinstance(e,ShellElement):
            a = e.area()
            m = e.mass()
        elif isinstance(e,SolidElement):
            v = e.volume()
            m = e.mass()
        elif isinstance(e,LineElement):
            L = e.length()
            m = e.mass()
        else:
            print "e.type = ",e.type
            #try:
            #    print "e.type = ",e.type
            #except:
            #    print str(e)
            ###
        ###
    ###

def compare(fem1,fem2):
    compareCardCount(fem1,fem2)
    getElementStats(fem1,fem2)
    #compareParams(fem1,fem2)
    #printPoints(fem1,fem2)


def compareParams(fem1,fem2):
    compute(fem1.params,fem2.params)
    
def printPoints(fem1,fem2):
    for nid,n1 in sorted(fem1.nodes.items()):
        print "%s   xyz=%s  n1=%s  n2=%s" %(nid,n1.xyz,n1.Position(True),  fem2.Node(nid).Position())
        break
    ###
    coord = fem1.Coord(5)
    print coord
    #print coord.Stats()

