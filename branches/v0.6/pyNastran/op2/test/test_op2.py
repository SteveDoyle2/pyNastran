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
import time
from traceback import print_exc

import pyNastran
from pyNastran.op2.op2 import OP2

def parse_table_names_from_F06(f06Name):
    """gets the op2 names from the f06"""
    infile = open(f06Name,'r')
    marker = 'NAME OF DATA BLOCK WRITTEN ON FORTRAN UNIT IS'
    names = []
    for line in infile:
        if marker in line:
            word = line.replace(marker,'').strip().strip('.')
            names.append(word)

    infile.close()
    return names

def get_failed_files(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    
    files = []
    for line in lines:
        files.append(line.strip())
    return files

def run_lots_of_files(files ,make_geom=True, write_bdf=False, write_f06=True,
                   write_matlab=True, delete_f06=True, print_results=True,
                   debug=True, saveCases=True, skipFiles=[],
                   stopOnFailure=False, nStart=0, nStop=1000000000):
    n = ''
    iSubcases = []
    failedCases = []
    nFailed = 0
    nTotal  = 0
    nPassed = 0
    t0 = time.time()
    for (i, op2file) in enumerate(files[nStart:nStop], nStart):  # 149
        baseName = os.path.basename(op2file)
        #if baseName not in skipFiles and not baseName.startswith('acms') and i not in nSkip:
        if baseName not in skipFiles and '#' not in op2file:
            print("%"*80)
            print('file=%s\n' % op2file)
            n = '%s ' %(i)
            sys.stderr.write('%sfile=%s\n' %(n, op2file))
            nTotal += 1
            isPassed = run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf,
                               write_f06=write_f06, write_matlab=write_matlab,
                               delete_f06=delete_f06, print_results=print_results,
                               iSubcases=iSubcases, debug=debug,
                               stopOnFailure=stopOnFailure) # True/False
            if not isPassed:
                sys.stderr.write('**file=%s\n' % op2file)
                failedCases.append(op2file)
                nFailed +=1
            else:
                nPassed +=1
            #sys.exit('end of test...test_op2.py')

    if saveCases:
        f = open('failedCases.in','wb')
        for op2file in failedCases:
            f.write('%s\n' % op2file)
        f.close()
    
    seconds = time.time()-t0
    minutes = seconds/60.
    print("dt = %s seconds = %s minutes" % (seconds, minutes))
    
    #op2 = OP2('test_tet10_subcase_1.op2')
    #op2.read_op2()
    
    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' %(nPassed,nTotal,100.*nPassed/float(nTotal),nTotal-nPassed)
    print(msg)
    sys.exit(msg)

def run_op2(op2FileName, make_geom=False, write_bdf=False, write_f06=True,
            write_matlab=True, is_mag_phase=False, delete_f06=False,
            print_results=True, iSubcases=[], debug=False, stopOnFailure=True):
    assert '.op2' in op2FileName.lower(), 'op2FileName=%s is not an OP2' %(op2FileName)
    isPassed = False
    stopOnFailure = False
    #debug = True
    write_matlab = False
    try:
        op2 = OP2(op2FileName, make_geom=make_geom, debug=debug)
        op2.set_subcases(iSubcases)

        #op2.read_bdf(op2.bdfFileName,includeDir=None,xref=False)
        #op2.write_bdf_as_patran()
        op2.read_op2()
        print("---stats for %s---" % op2FileName)
        #op2.get_op2_stats()
        print(op2.get_op2_stats())
        if write_bdf:
            op2.write_bdf(interspersed=True)
        #tableNamesF06 = parse_table_names_from_F06(op2.f06FileName)
        #tableNamesOP2 = op2.getTableNamesFromOP2()
        if write_f06:
            (model, ext) = os.path.splitext(op2FileName)
            op2.write_f06(model+'.f06.out', is_mag_phase=is_mag_phase)
            if delete_f06:
                os.remove(model+'.f06.out')

        if write_matlab:
            (model, ext) = os.path.splitext(op2FileName)
            op2.write_matlab(model+'.m', is_mag_phase=is_mag_phase)

        #print op2.print_results()
        if print_results:
            op2.print_results()
        #print "subcases = ",op2.subcases

        #assert tableNamesF06==tableNamesOP2,'tableNamesF06=%s tableNamesOP2=%s' %(tableNamesF06,tableNamesOP2)
        #op2.caseControlDeck.sol = op2.sol
        #print op2.caseControlDeck.get_op2_data()
        #print op2.print_results()
        #print op2.caseControlDeck.get_op2_data()
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % op2FileName)
        sys.exit('keyboard stop...')
    #except AddNewElementError:
    #    raise
    #except RuntimeError: # the op2 is bad, not my fault
    #    isPassed = True
    #    if stopOnFailure:
    #        raise
    #    else:
    #        isPassed = True

    #except IOError: # missing file
        #pass
    #except AssertionError:
    #    isPassed = True

    #except InvalidFormatCodeError:
    #    isPassed = True
    #except RuntimeError: #invalid analysis code
    #    isPassed = True
    #except SyntaxError: # Invalid Marker
        #isPassed = True

    #except EOFError:
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
    #except IOError: # missing file
        #isPassed = False
        #raise
    #except SyntaxError: #Param Parse
        #isPassed = True
    except:
        #print e
        print_exc(file=sys.stdout)
        if stopOnFailure:
            raise
        else:
            isPassed = False

    return isPassed

def run_arg_parse():
    import argparse

    ver = str(pyNastran.__version__)
    parser = argparse.ArgumentParser(description='Tests to see if an OP2 will work with pyNastran.',add_help=True) #,version=ver)
    parser.add_argument('op2FileName', metavar='op2FileName', type=str, nargs=1,
                       help='path to OP2 file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-q','--quiet',    dest='quiet',    action='store_true',help='Prints   debug messages (default=True)')

    #group2 = parser.add_mutually_exclusive_group()  # should this be exclusive???
    parser.add_argument('-g','--geometry',     dest='geometry',       action='store_true', help='Reads the OP2 for geometry, which can be written out')
    parser.add_argument('-w','--write_bdf',    dest='write_bdf',      action='store_true', help='Writes the bdf to fem.bdf.out')
    parser.add_argument('-f','--write_f06',    dest='write_f06',      action='store_true', help='Writes the f06 to fem.f06.out')
    parser.add_argument('-z','--is_mag_phase', dest='is_mag_phase',   action='store_true', help='F06/Matlab Writer writes Magnitude/Phase instead of Real/Imaginary (still stores Real/Imag)')
    parser.add_argument('-p','--print_results',dest='print_results',  action='store_true', help='Prints objects to screen which can require lots of memory')
    parser.add_argument('-v','--version',action='version',version=ver)
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    print("op2FileName = %s" %(args.op2FileName[0]))
    print("debug       = %s" %(not(args.quiet)))

    debug         = not(args.quiet)
    make_geom     = args.geometry
    write_bdf     = args.write_bdf
    write_f06     = args.write_f06
    is_mag_phase  = args.is_mag_phase
    write_matlab  = False
    print_results = args.print_results
    op2FileName   = args.op2FileName[0]

    return (op2FileName,make_geom,write_bdf,write_f06,write_matlab,is_mag_phase,print_results,debug)

def main():
    (op2FileName,make_geom,write_bdf,write_f06,write_matlab,is_mag_phase,print_results,debug) = run_arg_parse()

    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')
    run_op2(op2FileName,make_geom=make_geom,write_bdf=write_bdf,write_f06=write_f06,write_matlab=write_matlab,is_mag_phase=is_mag_phase,print_results=print_results,debug=debug)

if __name__=='__main__':  # op2
    main()
