import os
import sys
from numpy import pi, zeros, array, matrix #, float64, memmap
from numpy import log as naturalLog
from numpy.linalg import inv

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader

#from mapping.f06 import F06Reader
#from grom.py2matlab import MatPrint

from pyNastran.utils.log import get_logger
debug = True
log = get_logger(None, 'debug' if debug else 'info')

def readOP2(infilename):
    log.info('---starting deflectionReader.init of %s---' % infilename)
    op2 = OP2(infilename)
    terms = ['force','stress','stress_comp','strain','strain_comp','displacement','grid_point_forces']
    op2.read_op2()

    subcase0 = op2.displacements.keys()[0]  # get the 0th subcase
    displacment_obj = op2.displacements[subcase0]

    #op2.nastranModel.printDisplacement()
    #displacements = convertDisplacements(displacements)
    log.info('---finished deflectionReader.init of %s---' % infilename)
    return displacment_obj.translations


def readCart3dPoints(cfdGridFile):
    """return half model points to shrink xK matrix"""
    cart = Cart3DReader()
    (points, elements, regions, loads) = cart.read_cart3d(cfdGridFile)
    (points, elements, regions, loads) = cart.make_half_model(points, elements, regions, loads)
    return points

def writeNewCart3dMesh(cfdGridFile, cfdGridFile2, wA):
    """takes in half model wA, and baseline cart3d model, updates full model grids"""
    log.info("---starting writeNewCart3dMesh---")

    # make half model
    cart3d = Cart3DReader()
    result_names = ['Cp']
    (points, elements, regions, loads) = cart.read_cart3d(cfdGridFile, result_names=result_names) # reading full model
    (points, elements, regions, loads) = cart.make_half_model(points, elements, regions, loads)

    # adjusting points
    points2 = {}
    for (iPoint, point) in sorted(points.iteritems()):
        wai = wA[iPoint]
        (x, y, z) = point
        points2[iPoint] = [x, y, z + wai]

    (points, elements, regions, loads) = cart.make_mirror_model(points2, elements, regions, loads)  # mirroring model
    cart.writeOutfile(cfdGridFile2, points, elements, regions) # writing half model  (cleans up leftover parameters)

    log.info("---finished writeNewCart3dMesh---")
    sys.stdout.flush()


def removeDuplicateNodes(nodeList,mesh):
    """
    Removes nodes that have the same (x,y) coordinate.
    Note that if 2 nodes with different z values are found, only 1 is returned.
    This is intentional.   splineSurface = f(x,y)
    """
    nodeList.sort()
    log.info("nodeListA = %s" %(nodeList))
    nodeDict = {}
    for iNode in nodeList:
        (x,y,z) = mesh.Node(iNode).Position()
        nodeDict[(x,y)] = iNode
    nodeList = nodeDict.values()
    nodeList.sort()
    log.info("nodeListB = %s" %(nodeList))
    sys.stdout.flush()
    return nodeList

def runMapDeflections(nodeList, bdf_filename, out_filename, cart3d, cart3d2, log=None):
    fbase, ext = os.path.splitext(out_filename)
    if ext == '.op2':
        deflections = readOP2(out_filename)
    elif ext == '.f06':
        deflections = readF06(out_filename)
    else:
        raise NotImplementedError('out_filename = %r' % out_filename)

    mesh = BDF(debug=True, log=log)
    mesh.read_bdf(bdf_filename, include_dir=None, xref=True, punch=False)

    nodeList = removeDuplicateNodes(nodeList, mesh)
    C = getCmatrix(nodeList, mesh)
    wS = getWS(nodeList, deflections)
    del deflections

    aPoints = readCart3dPoints(cart3d)
    wA = getWA(nodeList, C, wS, mesh, aPoints)
    del C
    #del wS
    del mesh

    writeNewCart3dMesh(cart3d, cart3d2, wA)
    return (wA,wS)

def getWA(nodeList, C, wS, mesh, aPoints):
    log.info('---starting getWA---')
    MatPrint(sys.stdout,C)

    C  = inv(C) * wS  # Cws matrix, P matrix
    #P = solve(C, wS)
    #C*P=wS
    #P = C^-1*wS

    wA = getXK_Matrix(C,nodeList,mesh,aPoints)
    #wA = xK*C*wS
    log.info('---finished getWA---')
    sys.stdout.flush()
    return wA

def getXK_Matrix(Cws, nodeList, mesh, aPoints):
    log.info('---starting getXK_matrix---')
    D = 1.
    piD16 = pi*D*16.

    nNodes = len(nodeList)
    nPoints = len(aPoints.keys())
    wa = {}
    i = 0
    for (iAero, aNode) in sorted(aPoints.iteritems()):
        xK = zeros(nNodes+3, 'd')
        #nodeI = mesh.Node(iNode)

        xa,ya,za = aNode

        xK[0] = 1.
        xK[1] = xa
        xK[2] = ya

        j = 3
        for jNode in nodeList:
            sNode = mesh.Node(jNode)
            (xs, ys, zs) = sNode.Position()

            Rij2 = (xa-xs)**2. + (ya-ys)**2  # Rij^2
            if Rij2==0.:
                xK[j] = 0.
            else:
                Kij = Rij2 * naturalLog(Rij2) / piD16
                xK[j] = Kij
            j += 1

        wai = xK*Cws
        wa[iAero] = wai[0,0]
        #print "w[%s]=%s" % (iAero, wi[0,0])
        i += 1
    #print '---wa---'
    #print 'wa = ',wa
    log.info('---finished getXK_matrix---')
    sys.stdout.flush()
    return wa

def getWS(nodeList,deflections):
    log.info('---staring WS---')
    nNodes = len(nodeList)
    Wcolumn = matrix(zeros((3+nNodes,1), 'd'))
    i = 3
    for iNode in nodeList:
        (dx,dy,dz) = deflections[iNode]
        Wcolumn[i]=dz
        log.info("wS[%s=%s]=%s" %(iNode,i,dz))
        i+=1
    print max(Wcolumn)
    log.info('---finished getWS---')
    sys.stdout.flush()

    wSmax = max(Wcolumn)
    print "wSmax = ",wSmax[0,0]
    return Wcolumn

def getCmatrix(nodeList, mesh):
    log.info('---starting getCmatrix---')
    D = 1.
    piD16 = pi*D*16.

    nNodes = len(nodeList)
    i = 3
    #C = memmap('Cmatrix.map',dtype='float64',mode='write',shape=(3+nNodes,3+nNodes) )
    log.info('nNodes=%s' % nNodes)
    sys.stdout.flush()
    C = matrix(zeros((3+nNodes,3+nNodes), 'float64'))
    for iNode in nodeList:
        nodeI = mesh.Node(iNode)
        #i = iNode+3
        (xi, yi, zi) = nodeI.Position()
        #x,y,z = p

        C[0,i] = 1.
        C[1,i] = xi
        C[2,i] = yi

        C[i,0] = 1.
        C[i,1] = xi
        C[i,2] = yi

        j = 3
        for (jNode) in nodeList:
            #j = 3+jNode
            nodeJ = mesh.Node(jNode)
            (xj,yj,zj) = nodeJ.Position()
            if i==j:
                C[i,j] = 0.
            else:
                Rij2 = (xi-xj)**2. + (yi-yj)**2  # Rij^2
                if Rij2==0.:
                    C[i,j] = 0.
                else:
                    Kij = Rij2 * naturalLog(Rij2) / piD16
                    C[i,j] = Kij
                    #msg = "i=%s j=%s xi=%s xj=%s yi=%s yj=%s Rij2=%s Kij=%s" %(i,j,xi,xj,yi,yj,Rij2,Kij)
                    #assert isinstance(Kij,float64), msg
            j += 1
        i += 1
    log.info('---finished getCmatrix---')
    sys.stdout.flush()
    return C


if __name__=='__main__':
    basepath = os.getcwd()
    configpath = os.path.join(basepath, 'inputs')
    workpath   = os.path.join(basepath, 'outputsFinal')

    bdf_name = os.path.join(configpath, 'fem3.bdf')
    f06_name = os.path.join(configpath, 'fem3.f06')
    op2_name = os.path.join(workpath, 'fem3.f06')
    cart3d = os.path.join(configpath, 'Cart3d_bwb.i.tri')
    nodeList = [20037, 21140, 21787, 21028, 1151, 1886, 2018, 1477, 1023, 1116, 1201, 1116, 1201, 1828, 2589, 1373, 1315, 1571, 1507, 1532, 1317, 1327, 2011, 1445, 2352, 1564, 1878, 1402, 1196, 1234, 1252, 1679, 1926, 1274, 2060, 2365, 21486, 20018, 20890, 20035, 1393, 2350, 1487, 1530, 1698, 1782]
    #nodeList = [1001,1002,1003,1004,1005,1006]  # these are the hard points
    #nodeList = mesh.getNodeIDs() # [0:200]
    cart3d2 = cart3d + '_deflected'

    (wA,wS) = runMapDeflections(nodeList, bdf_name, op2_name, cart3d, cart3d2, log=log)
    print "wAero=", wA
    wSmax = max(wS)
    print "wSmax = ", wSmax[0,0]
