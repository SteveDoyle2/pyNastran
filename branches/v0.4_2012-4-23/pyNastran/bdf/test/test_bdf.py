## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import os
import sys
import numpy
numpy.seterr(all='raise')
import traceback

import pyNastran
from pyNastran.bdf.errors import *
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf import ShellElement,SolidElement,LineElement,RigidElement,SpringElement,PointElement,DamperElement
from pyNastran.bdf.test.compareCardContent import compareCardContent

import pyNastran.bdf.test
testPath = pyNastran.bdf.test.__path__[0]
#print "testPath = ",testPath

def runAllFilesInFolder(folder,debug=False,xref=True,check=True,cid=None):
    #debug = True
    print "folder = ",folder
    filenames  = os.listdir(folder)
    filenames2 = []
    diffCards = []
    for filename in filenames:
        if filename.endswith('.bdf') or filename.endswith('.dat') or filename.endswith('.nas') or filename.endswith('.nas'):
            filenames2.append(filename)
        ###
    ###
    for filename in filenames2:
        print "filename = ",os.path.abspath(os.path.join(folder,filename))
        try:
            (fem1,fem2,diffCards2) = runBDF(folder,filename,debug=debug,xref=xref,check=check,cid=cid,isFolder=True)
            diffCards += diffCards
        except KeyboardInterrupt:
            sys.exit('KeyboardInterrupt...sys.exit()')
        #except MissingFileError:
        #    pass
        #except TabCharacterError:
        #    pass
        #except ScientificParseError:
        #    pass
        #except ClosedBDFError:
        #    pass
        #except BDF_SyntaxError:
        #    pass
        except SystemExit:
            sys.exit('sys.exit...')
        except:
            traceback.print_exc(file=sys.stdout)
            #raise
        ###
        print '-'*80
    ###
    print '*'*80
    try:
        print "diffCards1 = ",list(set(diffCards))
    except TypeError:
        #print "type(diffCards) =",type(diffCards)
        print "diffCards2 = ",diffCards
###

def runBDF(folder,bdfFilename,debug=False,xref=True,check=True,cid=None,meshForm='combined',isFolder=False):
    bdfModel = str(bdfFilename)
    print "bdfModel = ",bdfModel
    if isFolder:
        bdfModel = os.path.join(testPath,folder,bdfFilename)
    
    assert os.path.exists(bdfModel),'|%s| doesnt exist' %(bdfModel)

    fem1 = BDF(debug=debug,log=None)
    fem1.log.info('starting fem1')
    sys.stdout.flush()
    fem2 = None
    diffCards = []
    try:
        #print "xref = ",xref
        try:
            fem1.readBDF(bdfModel,xref=xref)
        except:
            print "failed reading |%s|" %(bdfModel)
            raise
        #fem1.sumForces()
        #fem1.sumMoments()
        outModel = bdfModel+'_out'
        if cid is not None and xref:
            fem1.resolveGrids(cid=cid)
        if meshForm=='combined':
            fem1.writeBDFAsPatran(outModel)
        elif meshForm=='separate':
            fem1.writeBDF(outModel)
        else:
            raise NotImplementedError("meshForm=|%r| allowedForms=['combined','separate']" %(meshForm))
        #fem1.writeAsCTRIA3(outModel)

        fem2 = BDF(debug=debug,log=None)
        fem2.log.info('starting fem2')
        sys.stdout.flush()
        try:
            fem2.readBDF(outModel,xref=xref)
        except:
            print "failed reading |%s|" %(outModel)
            raise
        
        #fem2.sumForces()
        #fem2.sumMoments()
        outModel2 = bdfModel+'_out2'
        fem2.writeBDFAsPatran(outModel2)
        #fem2.writeAsCTRIA3(outModel2)
        diffCards = compare(fem1,fem2,xref=xref,check=check)
        os.remove(outModel2)

    except KeyboardInterrupt:
        sys.exit('KeyboardInterrupt...sys.exit()')
    #except MissingFileError:
    #    pass
    #except TabCharacterError:
    #    pass
    #except ScientificParseError:
    #    pass
    #except BDF_SyntaxError:
    #    pass
    #except ClosedBDFError:
    #    pass
    except SystemExit:
        sys.exit('sys.exit...')
    except:
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #print "\n"
        traceback.print_exc(file=sys.stdout)
        #print msg
        print "-"*80
        raise
    ###
    print "-"*80
    return (fem1,fem2,diffCards)

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
    return computeInts(cards1,cards2,fem1)

def computeInts(cards1,cards2,fem1):
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys  = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))
    
    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    msg = ' diffKeys1=%s diffKeys2=%s' %(diffKeys2,diffKeys2)
    print msg
    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1: value1 = cards1[key]
        else:                value1 = 0

        if key in listKeys2: value2 = cards2[key]
        else:                value2 = 0
        
        diff = abs(value1-value2)
        star = ' '
        if diff and key not in ['INCLUDE']:
            star = '*'
        if key not in fem1.cardsToRead:
            star = '-'

        factor1 = divide(value1,value2)
        factor2 = divide(value2,value1)
        factorMsg = ''
        if factor1 != factor2:
            factorMsg = 'diff=%s factor1=%g factor2=%g' %(diff,factor1,factor2)

        msg += '  %skey=%-7s value1=%-7s value2=%-7s' %(star,key,value1,value2)+factorMsg #+'\n'
        msg = msg.rstrip()
        print msg
    ###
    return listKeys1+listKeys2

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
        
        if key=='INCLUDE':
            msg += '    key=%-7s value1=%-7s value2=%-7s' %(key,value1,value2)
        else:
            msg += '   *key=%-7s value1=%-7s value2=%-7s' %(key,value1,value2)
        msg = msg.rstrip()
        print msg
    ###

def getElementStats(fem1,fem2):
    for key,e in sorted(fem1.elements.items()):
        try:
            if isinstance(e,ShellElement):
                a   = e.Area()
                t   = e.Thickness()
                nsm = e.Nsm()
                mA  = e.MassPerArea()
                m   = e.Mass()
                c   = e.Centroid()
                #mid = e.Mid()
                pid = e.Pid()
                n    = e.Normal()
                a,c,n = e.AreaCentroidNormal()
                
            elif isinstance(e,SolidElement):
                v   = e.Volume()
                m   = e.Mass()
                c   = e.Centroid()
                mid = e.Mid()
                pid = e.Pid()
            elif isinstance(e,LineElement): # ROD/BAR/BEAM
                L   = e.Length()
                nsm = e.Nsm()
                A   = e.Area()
                mL  = e.MassPerLength()
                m   = e.Mass()
                I22 = e.I22()
                I11 = e.I11()
                I12 = e.I12()
                J   = e.J()
                c   = e.Centroid()
                mid = e.Mid()
                pid = e.Pid()
                if J is None:
                    print "Moment of Inertia not available - e.type=%s e.eid=%i" %(e.type,e.eid)
            elif isinstance(e,RigidElement):
                pass
            elif isinstance(e,DamperElement):
                b = e.B()
            elif isinstance(e,SpringElement):
                L = e.Length()
                K = e.K()
                pid = e.Pid()
            elif isinstance(e,PointElement):
                m = e.Mass()
                c = e.Centroid()
            else:
                print "statistics not available - e.type=%s e.eid=%s" %(e.type,e.eid)
                #try:
                #    print "e.type = ",e.type
                #except:
                #    print str(e)
                ###
            ###
        except:
            print "*stats - e.type=%s eid=%s  element=\n%s" %(e.type,e.eid,str(e))
            raise
        ###
    ###

def compare(fem1,fem2,xref=True,check=True):
    diffCards = compareCardCount(fem1,fem2)
    if xref and check:
        getElementStats(fem1,fem2)
    compareCardContent(fem1,fem2)
    #compareParams(fem1,fem2)
    #printPoints(fem1,fem2)
    return diffCards

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

def main():
    import argparse

    ver = str(pyNastran.__version__)
    parser = argparse.ArgumentParser(description='Tests to see if a BDF will work with pyNastran.',add_help=True)
    parser.add_argument('bdfFileName', metavar='bdfFileName', type=str, nargs=1,
                       help='path to BDF/DAT file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-q','--quiet',  dest='quiet',action='store_true',  help='Prints   debug messages (default=False)')
    parser.add_argument('-x','--xref',   dest='xref', action='store_false', help='Disables cross-referencing and checks of the BDF')
    parser.add_argument('-c','--checks', dest='check',action='store_false', help='Disables BDF checks.  Checks run the methods on every element/property to test them.  May fails if a card is not supported.')
    parser.add_argument('-v','--version',action='version',version=ver,help="Shows pyNastran's version number and exits")
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    print "bdfFile     = ",args.bdfFileName[0]
    print "xref        = ",args.xref
    print "check       = ",args.check
    print "debug       = ",not(args.quiet)

    xref        = args.xref
    check       = args.check
    debug       = not(args.quiet)
    bdfFileName = args.bdfFileName[0]

    runBDF('.',bdfFileName,debug=debug,xref=xref,check=check)

if __name__=='__main__':
    main()

