from __future__ import print_function
from six import iteritems
import os
import sys
from numpy import pi, zeros, matrix, searchsorted #, float64, memmap
from numpy import log as naturalLog
from numpy.linalg import inv

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.f06.f06 import F06
from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader
#from pyNastran.applications.cart3d_nastran_fsi.mathFunctions import printMatrix

from pyNastran.utils.log import get_logger
debug = True
log = get_logger(None, 'debug' if debug else 'info')


def read_op2(op2_filename, isubcase=1):
    log.info('---starting deflectionReader.init of %s---' % op2_filename)
    op2 = OP2()
    op2.set_results('displacements')
    op2.read_op2(op2_filename)
    displacment_obj = op2.displacements[isubcase]

    log.info('---finished deflectionReader.init of %s---' % op2_filename)
    return displacment_obj


def read_f06(f06_filename, isubcase=1):
    log.info('---starting deflectionReader.init of %s---' % f06_filename)
    f06 = F06()
    #terms = ['force','stress','stress_comp','strain','strain_comp','displacement','grid_point_forces']
    f06.set_results('displacements')
    f06.read_f06(f06_filename)
    displacment_obj = f06.displacements[isubcase]

    #op2.nastranModel.printDisplacement()
    #displacements = convertDisplacements(displacements)
    log.info('---finished deflectionReader.init of %s---' % f06_filename)
    return displacment_obj.translations


def read_half_cart3d_points(cfdGridFile):
    """return half model points to shrink xK matrix"""
    cart = Cart3DReader()
    points, elements, regions, loads = cart.read_cart3d(cfdGridFile)
    points, elements, regions, loads = cart.make_half_model(points, elements, regions, loads)
    return points


def write_new_cart3d_mesh(cfdGridFile, cfdGridFile2, wA):
    """takes in half model wA, and baseline cart3d model, updates full model grids"""
    log.info("---starting write_new_cart3d_mesh---")

    # make half model
    cart = Cart3DReader()
    result_names = ['Cp']
    points, elements, regions, loads = cart.read_cart3d(cfdGridFile, result_names=result_names) # reading full model
    points, elements, regions, loads = cart.make_half_model(points, elements, regions, loads)

    # adjusting points
    points2 = {}
    for (iPoint, point) in sorted(iteritems(points)):
        wai = wA[iPoint]
        (x, y, z) = point
        points2[iPoint] = [x, y, z + wai]

    points, elements, regions, loads = cart.make_mirror_model(points2, elements, regions, loads)  # mirroring model
    cart.write_cart3d(cfdGridFile2, points, elements, regions) # writing half model; no loads (cleans up leftover parameters)

    log.info("---finished write_new_cart3d_mesh---")
    sys.stdout.flush()


def remove_duplicate_nodes(node_list, mesh):
    """
    Removes nodes that have the same (x,y) coordinate.
    Note that if 2 nodes with different z values are found, only 1 is returned.
    This is intentional.
    """
    node_list.sort()
    log.info("nodeListA = %s" % node_list)
    node_dict = {}
    for iNode in node_list:
        x, y, z = mesh.Node(iNode).get_position()
        node_dict[(x, y)] = iNode
    node_list = node_dict.values()
    node_list.sort()
    log.info("nodeListB = %s" % node_list)
    sys.stdout.flush()
    return node_list

def run_map_deflections(node_list, bdf_filename, out_filename, cart3d, cart3d2, log=None):
    fbase, ext = os.path.splitext(out_filename)
    if ext == '.op2':
        deflections = read_op2(out_filename)
    #elif ext == '.f06':
        #deflections = read_f06(out_filename)
    else:
        raise NotImplementedError('out_filename = %r' % out_filename)

    mesh = BDF(debug=True, log=log)
    mesh.read_bdf(bdf_filename, xref=True, punch=False)

    node_list = remove_duplicate_nodes(node_list, mesh)
    C = getCmatrix(node_list, mesh)
    wS = get_WS(node_list, deflections)
    del deflections

    aero_points = read_half_cart3d_points(cart3d)
    wA = get_WA(node_list, C, wS, mesh, aero_points)
    del C
    del mesh

    write_new_cart3d_mesh(cart3d, cart3d2, wA)
    return (wA, wS)

def get_WA(nodeList, C, wS, mesh, aero_points):
    log.info('---starting get_WA---')
    #print printMatrix(C)

    C = inv(C) * wS  # Cws matrix, P matrix
    #P = solve(C, wS)
    #C*P=wS
    #P = C^-1*wS

    wA = getXK_Matrix(C, nodeList, mesh, aero_points)
    #wA = xK*C*wS
    log.info('---finished get_WA---')
    sys.stdout.flush()
    return wA

def getXK_Matrix(Cws, nodeList, mesh, aero_points):
    log.info('---starting getXK_matrix---')
    D = 1.
    piD16 = pi * D * 16.

    nNodes = len(nodeList)
    #nPoints = len(aero_points.keys())
    wa = {}
    #i = 0
    for iAero, aNode in sorted(iteritems(aero_points)):
        xK = zeros(nNodes+3, 'd')
        #nodeI = mesh.Node(iNode)

        xa, ya, za = aNode

        xK[0] = 1.
        xK[1] = xa
        xK[2] = ya

        j = 3
        for jNode in nodeList:
            sNode = mesh.Node(jNode)
            (xs, ys, zs) = sNode.get_position()

            Rij2 = (xa-xs)**2. + (ya-ys)**2  # Rij^2
            if Rij2 == 0.:
                xK[j] = 0.
            else:
                Kij = Rij2 * naturalLog(Rij2) / piD16
                xK[j] = Kij
            j += 1

        wai = xK * Cws
        wa[iAero] = wai[0, 0]
        #print("w[%s]=%s" % (iAero, wi[0, 0]))
        #i += 1
    #print('---wa---')
    #print('wa = ', wa)
    log.info('---finished getXK_matrix---')
    sys.stdout.flush()
    return wa

def get_WS(node_list, deflections):
    log.info('---staring get_WS---')
    nNodes = len(node_list)
    Wcolumn = matrix(zeros((3 + nNodes, 1), 'd'))
    i = 3
    nodes = deflections.node_grid[:, 0]
    for inode in node_list:
        inodei = searchsorted(nodes, inode)
        dx, dy, dz = deflections.data[0, inodei, :2] #deflections[inode]
        Wcolumn[i] = dz
        log.info("wS[%s=%s]=%s" % (inode, i, dz))
        i += 1
    print(max(Wcolumn))
    log.info('---finished get_WS---')
    sys.stdout.flush()

    wSmax = max(Wcolumn)
    print("wSmax = %s" % wSmax[0, 0])
    return Wcolumn

def getCmatrix(node_list, mesh):
    log.info('---starting getCmatrix---')
    D = 1.
    piD16 = pi * D * 16.

    nNodes = len(node_list)
    i = 3
    #C = memmap('Cmatrix.map', dtype='float64', mode='write', shape=(3+nNodes, 3+nNodes) )
    log.info('nNodes=%s' % nNodes)
    sys.stdout.flush()
    C = matrix(zeros((3 + nNodes, 3 + nNodes), 'float64'))
    for inode in node_list:
        nodeI = mesh.Node(inode)
        #i = inode+3
        (xi, yi, zi) = nodeI.get_position()
        #x,y,z = p

        C[0, i] = 1.
        C[1, i] = xi
        C[2, i] = yi

        C[i, 0] = 1.
        C[i, 1] = xi
        C[i, 2] = yi

        j = 3
        for jNode in node_list:
            #j = 3+jNode
            nodeJ = mesh.Node(jNode)
            xj, yj, zj = nodeJ.get_position()
            if i == j:
                C[i, j] = 0.
            else:
                Rij2 = (xi-xj)**2. + (yi-yj)**2  # Rij^2
                if Rij2 == 0.:
                    C[i, j] = 0.
                else:
                    Kij = Rij2 * naturalLog(Rij2) / piD16
                    C[i, j] = Kij
                    #msg = "i=%s j=%s xi=%s xj=%s yi=%s yj=%s Rij2=%s Kij=%s" %(i,j,xi,xj,yi,yj,Rij2,Kij)
                    #assert isinstance(Kij,float64), msg
            j += 1
        i += 1
    log.info('---finished getCmatrix---')
    sys.stdout.flush()
    return C


def main():
    basepath = os.getcwd()
    configpath = os.path.join(basepath, 'inputs')
    workpath = os.path.join(basepath, 'outputsFinal')

    bdf_filename = os.path.join(configpath, 'fem3.bdf')
    f06_filename = os.path.join(configpath, 'fem3.f06')
    op2_filename = os.path.join(workpath, 'fem3.op2')
    cart3d = os.path.join(configpath, 'Cart3d_bwb.i.tri')
    node_list = [
        20037, 21140, 21787, 21028, 1151, 1886, 2018, 1477, 1023, 1116, 1201,
        1116, 1201, 1828, 2589, 1373, 1315, 1571, 1507, 1532, 1317, 1327, 2011,
        1445, 2352, 1564, 1878, 1402, 1196, 1234, 1252, 1679, 1926, 1274, 2060,
        2365, 21486, 20018, 20890, 20035, 1393, 2350, 1487, 1530, 1698, 1782,
    ]
    #node_list = [1001, 1002, 1003, 1004, 1005, 1006]  # these are the hard points
    #node_list = mesh.getNodeIDs() # [0:200]
    cart3d2 = cart3d + '_deflected'

    wA, wS = run_map_deflections(node_list, bdf_filename, op2_filename, cart3d, cart3d2, log=log)
    print("wAero = %s" % wA)
    wSmax = max(wS)
    print("wSmax = %s" % wSmax[0, 0])


if __name__ == '__main__':
    main()
