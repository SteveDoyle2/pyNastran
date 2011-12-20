import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.op2.op2    import OP2,EndOfFileError
from pyNastran.bdf.errors import *
from pyNastran.op2.op2Errors import *

def parseTableNamesFromF06(f06Name):
   """gets the op2 names from the f06"""
   infile = open(f06Name,'r')
   marker = 'NAME OF DATA BLOCK WRITTEN ON FORTRAN UNIT IS'
   names = []
   for line in infile:

       if marker in line:
           word = line.replace(marker,'').strip().strip('.')
           names.append(word)
       ###
   ###
   infile.close()
   return names

def getOp2Files(dirname,files):
    files2 = []
    for fname in files:
        if '.op2' in fname:
            files2.append(os.path.join(dirname,fname))
        ###
    ###
    return files2

def getFailedFiles(filename):
    infile = open(filename,'r')
    lines = infile.readlines()
    infile.close()
    
    files = []
    for line in lines:
        files.append(line.strip())
    return files

def runOP2(op2file,makeGeom=False,writeBDF=False,debug=False,stopOnFailure=True):
    isPassed = False
    try:
        nTotal += 1
        op2 = OP2(op2file,makeGeom=makeGeom,debug=debug)
        op2.setSubcases(iSubcases)
        op2.readBDF(op2.bdfFileName,includeDir=None,xref=False)
        #op2.writeBDFAsPatran()
        op2.readOP2()
        if writeBDF:
            op2.writeBDFAsPatran()
        #tableNamesF06 = parseTableNamesFromF06(op2.f06FileName)
        tableNamesOP2 = op2.getTableNamesFromOP2()
        print op2.printResults()
        #print "subcases = ",op2.subcases

        #assert tableNamesF06==tableNamesOP2,'tableNamesF06=%s tableNamesOP2=%s' %(tableNamesF06,tableNamesOP2)
        pass
        #op2.caseControlDeck.sol = op2.sol
        #print op2.caseControlDeck.getOp2Data()
        #print op2.printResults()
        #print op2.caseControlDeck.getOp2Data()
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' %(op2file))
        sys.exit('keyboard stop...')
    except AddNewElementError:
        raise
    except TapeCodeError: # the op2 is bad, not my fault
        isPassed = True

    #except AssertionError:
    #    isPassed = True

    #except InvalidFormatCodeError:
    #    isPassed = True
    #except InvalidAnalysisCodeError:
    #    isPassed = True
    #except InvalidMarkersError:
    #    isPassed = True

    #except TabCharacterError:
    #    isPassed = True
    #except EndOfFileError:
    #    isPassed = True
    except SystemExit:
        #print_exc(file=sys.stdout)
        #sys.exit('stopping on sys.exit')
        raise
    #except NameError:  # variable isnt defined
    #    if stopOnFailure:
    #        raise
    #    else:
    #        isPassed = True
    #except AttributeError:  # missing function
    #    if stopOnFailure:
    #        raise
    #    else:
    #        isPassed = True
    #except KeyError:
    #    raise
    #except TypeError:  # numpy error
    #    isPassed = True
    #except IndexError: # bad bdf
    #    isPassed = True
    #except MissingFileError: # missing bdf file
    #    isPassed = True
    #except InvalidSubcaseParseError:
    #    isPassed = True
    #except ScientificParseError:  # bad value parsing
    #    isPassed = True
    #except ParamParseError:
    #    isPassed = True
    #except NotImplementedMethodError:
    #    isPassed = True
    #except InvalidFieldError: # bad bdf field
    #    isPassed = True
    except:
        #print e
        print_exc(file=sys.stdout)
        if stopOnFailure:
            raise
        else:
            isPassed = True
        ###
    return isPassed
    ###

if __name__=='__main__':  # op2
def main(op2FileName):
    import argparse

    ver = str(pyNastran.__version__)
    parser = argparse.ArgumentParser(description='Tests to see if an OP2 will work with pyNastran.',add_help=True) #,version=ver)
    parser.add_argument('op2FileName', metavar='op2FileName', type=str, nargs=1,
                       help='path to OP2 file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-q','--quiet',    dest='quiet',    action='store_true',help='Prints   debug messages (default=True)')

    parser.add_argument('-g','--geometry', dest='geometry', action='store_true',help='Reads the OP2 for geometry, which can be written out')
    parser.add_argument('-w','--writeBDF', dest='writeBDF', action='store_true',help='Writes the bdf to fem.bdf.out (requires --geometry)')
    parser.add_argument('-v','--version',action='version',version=ver)
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    print "op2FileName = ",args.op2FileName[0]
    print "debug       = ",not(args.quiet)

    debug       = not(args.quiet)
    makeGeom    = args.geometry
    writeBDF    = args.writeBDF
    op2FileName = args.op2FileName[0]

    runOP2(op2file,makeGeom=makeGeom,writeBDF=writeBDF,debug=debug)

