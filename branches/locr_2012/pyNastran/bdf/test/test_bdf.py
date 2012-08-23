# pylint: disable=W0612,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
import numpy
import warnings
warnings.simplefilter('always')

numpy.seterr(all='raise')
import traceback

from pyNastran.general.general import printBadPath
from pyNastran.bdf.bdf import BDF, CTRIAX, CTRIAX6
from pyNastran.bdf.bdf import (ShellElement, SolidElement, LineElement,
                               RigidElement, SpringElement, PointElement,
                               DamperElement, RodElement, NastranMatrix)
from pyNastran.bdf.test.compareCardContent import compare_card_content

import pyNastran.bdf.test
test_path = pyNastran.bdf.test.__path__[0]
#print "test_path = ",test_path


def run_all_files_in_folder(folder, debug=False, xref=True, check=True,
                            cid=None):
    print("folder = %s" % (folder))
    filenames = os.listdir(folder)
    run_lots_of_files(filenames, debug=debug, xref=xref, check=check, cid=cid)


def run_lots_of_files(filenames, folder='', debug=False, xref=True, check=True,
                      cid=None):
    filenames = list(set(filenames))
    filenames.sort()

    #debug = True
    filenames2 = []
    diffCards = []
    for filename in filenames:
        if (filename.endswith('.bdf') or filename.endswith('.dat') or
            filename.endswith('.nas') or filename.endswith('.nas')):
            filenames2.append(filename)

    failedFiles = []
    for filename in filenames2:
        absFilename = os.path.abspath(os.path.join(folder, filename))
        if folder != '':
            print("filename = %s" % (absFilename))
        isPassed = False
        try:
            (fem1, fem2, diffCards2) = runBDF(folder, filename, debug=debug,
                                              xref=xref, check=check, cid=cid,
                                              isFolder=True)
            diffCards += diffCards
            isPassed = True
        except KeyboardInterrupt:
            sys.exit('KeyboardInterrupt...sys.exit()')
        #except SyntaxError:
        #    pass
        except SystemExit:
            sys.exit('sys.exit...')
        except:
            traceback.print_exc(file=sys.stdout)
            #raise
        ###
        print('-' * 80)
        if not isPassed:
            failedFiles.append(absFilename)
    ###
    print('*' * 80)
    try:
        print("diffCards1 = %s" % (list(set(diffCards))))
    except TypeError:
        #print "type(diffCards) =",type(diffCards)
        print("diffCards2 = %s" % (diffCards))
    return failedFiles
###


def runBDF(folder, bdfFilename, debug=False, xref=True, check=True, cid=None,
           meshForm='combined', isFolder=False):
    bdfModel = str(bdfFilename)
    print("bdfModel = %s" % (bdfModel))
    if isFolder:
        bdfModel = os.path.join(test_path, folder, bdfFilename)

    assert os.path.exists(bdfModel), '|%s| doesnt exist' % (bdfModel)

    fem1 = BDF(debug=debug, log=None)
    fem1.log.info('starting fem1')
    sys.stdout.flush()
    fem2 = None
    diffCards = []
    try:
        #print("xref = ", xref)
        (outModel) = run_fem1(fem1, bdfModel, meshForm, xref, cid)

        (fem2) = run_fem2(bdfModel, outModel, xref, debug=debug, log=None)
        (diffCards) = compare(fem1, fem2, xref=xref, check=check)

    except KeyboardInterrupt:
        sys.exit('KeyboardInterrupt...sys.exit()')
    #except IOError:
        #pass
    except SystemExit:
        sys.exit('sys.exit...')
    except:
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #print "\n"
        traceback.print_exc(file=sys.stdout)
        #print msg
        print("-" * 80)
        raise
    ###
    print("-" * 80)
    return (fem1, fem2, diffCards)


def run_fem1(fem1, bdfModel, meshForm, xref, cid):
    assert os.path.exists(bdfModel), printBadPath(bdfModel)
    try:
        if '.pch' in bdfModel:
            fem1.readBDF_Punch(bdfModel, xref=False)
        else:
            fem1.readBDF(bdfModel, xref=xref)
    except:
        print("failed reading |%s|" % (bdfModel))
        raise
    #fem1.sumForces()
    #fem1.sumMoments()
    outModel = bdfModel + '_out'
    if cid is not None and xref:
        fem1.resolveGrids(cid=cid)
    if meshForm == 'combined':
        fem1.writeBDFAsPatran(outModel)
    elif meshForm == 'separate':
        fem1.writeBDF(outModel)
    else:
        msg = "meshForm=|%r| allowedForms=['combined','separate']" % (meshForm)
        raise NotImplementedError(msg)
    #fem1.writeAsCTRIA3(outModel)
    return (outModel)


def run_fem2(bdfModel, outModel, xref, debug=False, log=None):
    assert os.path.exists(bdfModel), bdfModel
    assert os.path.exists(outModel), outModel
    fem2 = BDF(debug=debug, log=log)
    fem2.log.info('starting fem2')
    sys.stdout.flush()
    try:
        fem2.readBDF(outModel, xref=xref)
    except:
        print("failed reading |%s|" % (outModel))
        raise

    #fem2.sumForces()
    #fem2.sumMoments()
    outModel2 = bdfModel + '_out2'
    fem2.writeBDFAsPatran(outModel2)
    #fem2.writeAsCTRIA3(outModel2)
    os.remove(outModel2)
    return (fem2)


def divide(value1, value2):
    if value1 == value2:  # good for 0/0
        return 1.0
    else:
        try:
            v = value1 / float(value2)
        except ZeroDivisionError:
            v = 0.
    return v


def compare_card_count(fem1, fem2):
    cards1 = fem1.cardCount
    cards2 = fem2.cardCount
    return compute_ints(cards1, cards2, fem1)


def compute_ints(cards1, cards2, fem1):
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))

    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    if diffKeys1 or diffKeys2:
        print(' diffKeys1=%s diffKeys2=%s' % (diffKeys1, diffKeys2))

    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1:
            value1 = cards1[key]
        else:
            value1 = 0

        if key in listKeys2:
            value2 = cards2[key]
        else:
            value2 = 0

        diff = abs(value1 - value2)
        star = ' '
        if diff and key not in ['INCLUDE']:
            star = '*'
        if key not in fem1.cardsToRead:
            star = '-'

        factor1 = divide(value1, value2)
        factor2 = divide(value2, value1)
        factorMsg = ''
        if factor1 != factor2:
            factorMsg = 'diff=%s factor1=%g factor2=%g' % (diff, factor1,
                                                          factor2)
        msg += '  %skey=%-7s value1=%-7s value2=%-7s' % (star, key, value1,
                                                       value2) + factorMsg  # +'\n'
        msg = msg.rstrip()
        print(msg)
    ###
    return listKeys1 + listKeys2


def compute(cards1, cards2):
    cardKeys1 = set(cards1.keys())
    cardKeys2 = set(cards2.keys())
    allKeys = cardKeys1.union(cardKeys2)
    diffKeys1 = list(allKeys.difference(cardKeys1))
    diffKeys2 = list(allKeys.difference(cardKeys2))

    listKeys1 = list(cardKeys1)
    listKeys2 = list(cardKeys2)
    msg = ''
    if diffKeys1 or diffKeys2:
        msg = 'diffKeys1=%s diffKeys2=%s' % (diffKeys1, diffKeys2)

    for key in sorted(allKeys):
        msg = ''
        if key in listKeys1:
            value1 = cards1[key]
        else:
            value2 = 0

        if key in listKeys2:
            value2 = cards2[key]
        else:
            value2 = 0

        if key == 'INCLUDE':
            msg += '    key=%-7s value1=%-7s value2=%-7s' % (key,
                                                             value1, value2)
        else:
            msg += '   *key=%-7s value1=%-7s value2=%-7s' % (key,
                                                             value1, value2)
        msg = msg.rstrip()
        print(msg)
    ###


def get_element_stats(fem1, fem2):
    """verifies that the various element methods work"""
    for (key, loads) in sorted(fem1.loads.iteritems()):
        for load in loads:
            try:
                allLoads = load.getLoads()
                if not isinstance(allLoads, list):
                    raise TypeError('allLoads should return a list...%s'
                                    % (type(allLoads)))
            except:
                print("load statistics not available - load.type=%s "
                      "load.sid=%s" % (load.type, load.sid))
                raise

    for (key, e) in sorted(fem1.elements.iteritems()):
        try:
            if isinstance(e, ShellElement):
                a = e.Area()
                if (isinstance(e, CTRIAX) or isinstance(e, CTRIAX6)):
                    pass
                else:
                    t = e.Thickness()
                    nsm = e.Nsm()
                    mA = e.MassPerArea()
                    m = e.Mass()
                    c = e.Centroid()
                    #mid = e.Mid()
                    pid = e.Pid()
                    n = e.Normal()
                    a, c, n = e.AreaCentroidNormal()
            elif isinstance(e, SolidElement):
                v = e.Volume()
                m = e.Mass()
                c = e.Centroid()
                mid = e.Mid()
                pid = e.Pid()
            elif isinstance(e, LineElement):  # ROD/BAR/BEAM
                L = e.Length()
                nsm = e.Nsm()
                A = e.Area()
                mL = e.MassPerLength()
                m = e.Mass()
                I22 = e.I22()
                I11 = e.I11()
                I12 = e.I12()
                J = e.J()
                c = e.Centroid()
                mid = e.Mid()
                pid = e.Pid()
                if J is None:
                    print("Moment of Inertia not available - e.type=%s "
                          "e.eid=%i" % (e.type, e.eid))
            elif isinstance(e, RodElement):  # CROD, CONROD, CTUBE
                L = e.Length()
                nsm = e.Nsm()
                A = e.Area()
                mL = e.MassPerLength()
                m = e.Mass()
                c = e.Centroid()
                mid = e.Mid()
                pid = e.Pid()
            elif isinstance(e, RigidElement):
                pass
            elif isinstance(e, DamperElement):
                b = e.B()
            elif isinstance(e, SpringElement):
                #L = e.Length()
                K = e.K()
                pid = e.Pid()
            elif isinstance(e, PointElement):
                m = e.Mass()
                c = e.Centroid()
            else:
                print("statistics not available - e.type=%s e.eid=%s"
                    % (e.type, e.eid))
                #try:
                #    print("e.type = ",e.type)
                #except:
                #    print(str(e))
                ###
            ###
        except:
            print("*stats - e.type=%s eid=%s  element=\n%s"
                % (e.type, e.eid, str(e)))
            raise
        ###
    ###


def get_matrix_stats(fem1, fem2):
    for (key, dmig) in sorted(fem1.dmigs.iteritems()):
        try:
            if isinstance(dmig, NastranMatrix):
                dmig.getMatrix()
            else:
                print("statistics not available - "
                      "matrix.type=%s matrix.name=%s" % (dmig.type, dmig.name))
        except:
            print("*stats - matrix.type=%s name=%s  matrix=\n%s"
                % (dmig.type, dmig.name, str(dmig)))
            raise
        ###


def compare(fem1, fem2, xref=True, check=True):
    diffCards = compare_card_count(fem1, fem2)
    if xref and check:
        get_element_stats(fem1, fem2)
        get_matrix_stats(fem1, fem2)
    compare_card_content(fem1, fem2)
    #compare_params(fem1,fem2)
    #print_points(fem1,fem2)
    return diffCards


def compare_params(fem1, fem2):
    compute(fem1.params, fem2.params)


def print_points(fem1, fem2):
    for (nid, n1) in sorted(fem1.nodes.iteritems()):
        print("%s   xyz=%s  n1=%s  n2=%s" % (nid, n1.xyz, n1.Position(True),
                                            fem2.Node(nid).Position()))
        break
    coord = fem1.Coord(5)
    print(coord)
    #print coord.Stats()


def main():
    import argparse

    ver = str(pyNastran.__version__)
    parser = argparse.ArgumentParser(description='Tests to see if a BDF will '
                                     'work with pyNastran.', add_help=True)
    parser.add_argument('bdfFileName', metavar='bdfFileName', type=str,
                        nargs=1, help='path to BDF/DAT file')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                       help='Prints   debug messages (default=False)')
    parser.add_argument('-x', '--xref', dest='xref', action='store_false',
                       help='Disables cross-referencing and checks of the BDF')
    parser.add_argument('-c', '--checks', dest='check', action='store_false',
                       help='Disables BDF checks.  Checks run the methods on '
                       'every element/property to test them.  May fails if a '
                       'card is not supported.')
    parser.add_argument('-v', '--version', action='version', version=ver,
                       help="Shows pyNastran's version number and exits")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    print("bdfFile     = %s" % (args.bdfFileName[0]))
    print("xref        = %s" % (args.xref))
    print("check       = %s" % (args.check))
    print("debug       = %s" % (not(args.quiet)))

    xref = args.xref
    check = args.check
    debug = not(args.quiet)
    bdf_filename = args.bdfFileName[0]

    runBDF('.', bdf_filename, debug=debug, xref=xref, check=check)

if __name__ == '__main__':
    main()
