import os
import sys
import time
from traceback import print_exc

import pyNastran
from pyNastran.pch.pch import PCH

def parse_table_names_from_F06(f06Name):
    """gets the pch names from the f06"""
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

def run_lots_of_files(files ,makeGeom=True, writeBDF=False, write_f06=True,
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
    for (i, pchfile) in enumerate(files[nStart:nStop], nStart):  # 149
        baseName = os.path.basename(pchfile)
        #if baseName not in skipFiles and not baseName.startswith('acms') and i not in nSkip:
        if baseName not in skipFiles and '#' not in pchfile:
            print("%"*80)
            print('file=%s\n' % pchfile)
            n = '%s ' %(i)
            sys.stderr.write('%sfile=%s\n' %(n, pchfile))
            nTotal += 1
            isPassed = run_pch(pchfile, makeGeom=makeGeom, writeBDF=writeBDF,
                               write_f06=write_f06, write_matlab=write_matlab,
                               delete_f06=delete_f06, print_results=print_results,
                               iSubcases=iSubcases, debug=debug,
                               stopOnFailure=stopOnFailure) # True/False
            if not isPassed:
                sys.stderr.write('**file=%s\n' % pchfile)
                failedCases.append(pchfile)
                nFailed +=1
            else:
                nPassed +=1
            #sys.exit('end of test...test_pch.py')

    if saveCases:
        f = open('failedCases.in','wb')
        for pchfile in failedCases:
            f.write('%s\n' % pchfile)
        f.close()
    
    seconds = time.time()-t0
    minutes = seconds/60.
    print("dt = %s seconds = %s minutes" %(seconds,minutes))
    
    #pch = PCH('test_tet10_subcase_1.pch')
    #pch.read_pch()
    
    msg = '-----done with all models %s/%s=%.2f%%  nFailed=%s-----' %(nPassed,nTotal,100.*nPassed/float(nTotal),nTotal-nPassed)
    print(msg)
    sys.exit(msg)

def run_pch(pchFileName, makeGeom=False, writeBDF=False, write_f06=True,
            write_matlab=True, isMagPhase=False, delete_f06=False,
            print_results=True, iSubcases=[], debug=False, stopOnFailure=True):
    assert '.pch' in pchFileName.lower(), 'pchFileName=%s is not an PCH' % pchFileName
    isPassed = False
    stopOnFailure = False
    #debug = True
    try:
        pch = PCH(pchFileName, makeGeom=makeGeom, debug=debug)
        pch.set_subcases(iSubcases)

        #pch.read_bdf(pch.bdfFileName,includeDir=None,xref=False)
        #pch.write_bdf_as_patran()
        pch.read_pch()
        print("---stats for %s---" % pchFileName)
        pch.get_pch_stats()
        #print(pch.get_pch_stats())
        if writeBDF:
            pch.write_bdf_as_patran()
        #tableNamesF06 = parse_table_names_from_F06(pch.f06FileName)
        #tableNamesPCH = pch.getTableNamesFromPCH()
        if write_f06:
            (model, ext) = os.path.splitext(pchFileName)
            pch.write_f06(model+'.f06.out', isMagPhase=isMagPhase)
            if delete_f06:
                os.remove(model+'.f06.out')

        if write_matlab:
            (model, ext) = os.path.splitext(pchFileName)
            pch.write_matlab(model+'.m', isMagPhase=isMagPhase)

        #print pch.print_results()
        if print_results:
            pch.print_results()
        #print "subcases = ",pch.subcases

        #assert tableNamesF06==tableNamesPCH,'tableNamesF06=%s tableNamesPCH=%s' %(tableNamesF06,tableNamesPCH)
        #pch.caseControlDeck.sol = pch.sol
        #print pch.caseControlDeck.get_pch_data()
        #print pch.print_results()
        #print pch.caseControlDeck.get_pch_data()
        isPassed = True
    except KeyboardInterrupt:
        sys.stdout.flush()
        print_exc(file=sys.stdout)
        sys.stderr.write('**file=%s\n' % pchFileName)
        sys.exit('keyboard stop...')
    #except AddNewElementError:
    #    raise
    #except RuntimeError: # the pch is bad, not my fault
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
    except IOError: # missing file
        isPassed = False
        raise
    except SyntaxError: #Param Parse
        isPassed = True
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
    parser = argparse.ArgumentParser(description='Tests to see if an PCH will work with pyNastran.',add_help=True) #,version=ver)
    parser.add_argument('pchFileName', metavar='pchFileName', type=str, nargs=1,
                       help='path to PCH file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-q','--quiet',    dest='quiet',    action='store_true',help='Prints   debug messages (default=True)')

    #group2 = parser.add_mutually_exclusive_group()  # should this be exclusive???
    parser.add_argument('-g','--geometry', dest='geometry',    action='store_true', help='Reads the PCH for geometry, which can be written out')
    parser.add_argument('-w','--writeBDF', dest='writeBDF',    action='store_true', help='Writes the bdf to fem.bdf.out')
    parser.add_argument('-f','--write_f06', dest='write_f06',    action='store_true', help='Writes the f06 to fem.f06.out')
    parser.add_argument('-z','--magPhase', dest='isMagPhase',  action='store_true', help='F06/Matlab Writer writes Magnitude/Phase instead of Real/Imaginary (still stores Real/Imag)')
    parser.add_argument('-m','--matlab',   dest='write_matlab', action='store_true', help='Matlab Writer is enabled (fails for transient; limited support)')
    parser.add_argument('-p','--print_results',dest='print_results',  action='store_true', help='Prints objects to screen which can require lots of memory')
    parser.add_argument('-v','--version',action='version',version=ver)
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    print("pchFileName = %s" %(args.pchFileName[0]))
    print("debug       = %s" %(not(args.quiet)))

    debug         = not(args.quiet)
    makeGeom      = args.geometry
    writeBDF      = args.writeBDF
    write_f06     = args.write_f06
    isMagPhase    = args.isMagPhase
    write_matlab  = args.write_matlab
    print_results = args.print_results
    pchFileName   = args.pchFileName[0]

    return (pchFileName,makeGeom,writeBDF,write_f06,write_matlab,isMagPhase,print_results,debug)

def main():
    (pchFileName,makeGeom,writeBDF,write_f06,write_matlab,isMagPhase,print_results,debug) = run_arg_parse()

    if os.path.exists('skippedCards.out'):
        os.remove('skippedCards.out')
    run_pch(pchFileName,makeGeom=makeGeom,writeBDF=writeBDF,write_f06=write_f06,write_matlab=write_matlab,isMagPhase=isMagPhase,print_results=print_results,debug=debug)

if __name__=='__main__':  # pch
    main()
