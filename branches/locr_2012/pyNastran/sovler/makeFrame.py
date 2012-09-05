import sys
from itertools import izip
from numpy import array, zeros, matrix
from numpy.linalg import solve

from pyNastran.bdf.bdf import BDF
from pyNastran.general.general import ListPrint
from pyNastran.general.mathematics import reduce_matrix


def makeTruss2():
    """
    @code
                  ^
                  |1000  lb
                  |
     12           |
    *
    *1-----A------2  34
     | D        B |
     |   \    /   C
     F    /  \    |   20"
     |  /       \ |
    *3-----E------4  78
    *56
          40"
    @endcode
    
    where:
      * \* indicates a constraint
      * 1234 are nodes
      * ABCDEF are rods
      * the truss is 20" by 40"
      * E = 1e7 psi
      * B and D do not touch, but cross at the diagonal
      * 12,34,56,78 are the degrees of freedom
    """
    model = BDF('')
    #model.executiveControlLines = ['SOL 101']
    #model.caseControlDeck = None
    #model.iSolLine = 0

    p1 = array([0., 20., 0.])
    p2 = array([40., 20., 0.])
    p3 = array([0., 0., 0.])
    p4 = array([40., 0., 0.])

    E = 1e7  # psi
    nu = 0.3
    G = None  # default is used
    J = 0.  # not needed

    nodes = [p1, p2, p3, p4]
    for i, node in enumerate(nodes):
        spc = None
        if node[0] == 0.:  # if x==0
            spc = str(123456)
        grid = ['GRID', i + 1, None, node[0], node[1], node[2], None, spc]
        model.add_card(grid, 'GRID')

          #id
    cA = [1, 2]
    cB = [2, 3]
    cC = [2, 4]
    cD = [1, 4]
    cE = [3, 4]
    cF = [1, 3]
    cRods = [cA, cB, cC, cD, cE, cF]

    mid1 = 199
    mat1 = ['MAT1', mid1, E, G, nu]
    A = [1.0, 0.5, 1.0, 0.5, 1.0, 1.0]

    i = 0
    i = addRods(model, cRods, mat1, i, A, J)

    scaleFactor = 1.
    loadID = 123
    load = ['LOAD', loadID, scaleFactor, 1., loadID]
    model.add_card(load, 'LOAD')

    force = ['FORCE', loadID, 2, None, 1000., 0., 1., 0.]
    model.add_card(force, 'FORCE')
    model.writeBDF('conrod.bdf')
    print "done"
    return model


def makeTruss():
    model = BDF('')

    p1 = array([0., 0., 0.]) * 12.
    p2 = array([10., 0., 0.]) * 12.
    p3 = array([10., 10., 0.]) * 12.
    p4 = array([0., 10., 0.]) * 12.

    p5 = array([0., 0., 10.]) * 12.
    p6 = array([10., 0., 10.]) * 12.
    p7 = array([10., 10., 10.]) * 12.
    p8 = array([0., 10., 10.]) * 12.

    nodes = [p1, p2, p3, p4, p5, p6, p7, p8]
    for i, node in enumerate(nodes):
        spc = None
        if node[2] == 0.:
            spc = str(123456)
        grid = ['GRID', i + 1, None, node[0], node[1], node[2], None, spc]
        model.add_card(grid, 'GRID')

    # c conrods
    Ac = 4.   # in^2
    Jc = 60.  # in4
    Ixc = 650.  # in4
    Iyc = 65.  # in4
    Ec = 29000.
    Gc = None
    nub = 0.3

          #id
    c1 = [5, 1]
    c2 = [6, 2]
    c3 = [7, 3]
    c4 = [8, 4]
    cRods = [c1, c2, c3, c4]

    mid1 = 99
    mat1 = ['MAT1', mid1, Ec, Gc, nub]

    i = 1
    i = addRods(model, cRods, mat1, i, Ac, Jc)

    # b conrods
    Ab = 3.2
    Jb = 43.
    Ixb = 450.
    Iyb = 32.
    Eb = 29000.
    Gb = None
    nub = 0.3

    b1 = [1, 2]
    b2 = [2, 3]
    b3 = [3, 4]
    b4 = [4, 1]

    b5 = [5, 6]
    b6 = [6, 7]
    b7 = [7, 8]
    b8 = [8, 5]
    bRods = [b1, b2, b3, b4, b6, b7, b8]

    mid2 = 22
    mat1 = ['MAT1', mid2, Eb, Gb, nub]
    i = addRods(model, bRods, mat1, i, Ab, Jb)

    scaleFactor = 1.
    loadID = 123
    load = ['LOAD', loadID, scaleFactor, 1., loadID]
    model.add_card(load, 'LOAD')

    force = ['FORCE', loadID, 3, None, 100., 1., 0., 0.]
    model.add_card(force, 'FORCE')
    model.writeBDF('frame.bdf')
    print "done"
    return model


def addRods(model, rods, mat1, i, A, J):
    if isinstance(A, float):
        A = [A] * len(rods)

    model.add_card(mat1, 'MAT1')
    for (area, rod) in izip(A, rods):
        (n1, n2) = rod
        #print "rod = ",rod           #mid
        conrod = ['CONROD', i + 1, n1, n2, mat1[1], area, J]
        model.add_card(conrod, 'CONROD')
        #print conrod
        #Ke = model.Element(i+1).Stiffness(model)
        i += 1

    return i


def buildGlobalStiffness(model):
    nodeIDs = model.nodeIDs()
    nDOFperNode = 2

    Dofs = {}
    nDOF = 0
    for nodeID in nodeIDs:
        Dofs[nodeID] = [nDOF, nDOF + 1, nDOF + 2]
        nDOF += 3

    nElements = model.nElements()
    Kg = matrix(zeros((nDOF, nDOF), 'd'))  # K_global

    for id, element in sorted(model.elements.iteritems()):
        #nodes = element.nodes
        print element
        Ke = element.Stiffness(model)
        #print "K_element[%s] = \n%s\n" %(id,Ke)

        nodes = element.nodeIDs()
        print "nodes = ", nodes

        dofs = []
        for i, iNode in enumerate(nodes):
            dofs += Dofs[iNode]
        print "allDOFs = ", dofs

        # put in global stiffness matrix
        for j, dof in enumerate(dofs):
            print "nid=%s dofs=%s" % (iNode, dof)
            for k, dof2 in enumerate(dofs):
                print "Kg[%s][%s] = Ke[%s][%s] = %s" % (
                    dof, dof2, j, k, Ke[j, k])
                #print "Kg[%s][%s] = Ke[%s][%s] = %s" %(dof2,dof, k,j,Ke[k,j])
                #print "Kg[%s][%s] = Ke[%s][%s] = %s" %(dof2,dof2,k,k,Ke[k,k])
                #print "Kg[%s][%s] = Ke[%s][%s] = %s" %(dof, dof, j,j,Ke[j,j])
                Kg[dof, dof2] += Ke[j, k]
                #Kg[dof2,dof ] += Ke[k,j]
                #Kg[dof2,dof2] += Ke[k,k]
                #Kg[dof ,dof ] += Ke[j,j]

        print "K_global =\n", ListPrint(Kg)
    print "K_global =\n", ListPrint(Kg)
    #print "K_global = \n",Kg
    return Kg, Dofs

# B - strain displacement matrix
# Uo = 1/2*d.T*B.T*E*B
# Ke = integral(B.T*E*B,dV)


def applyBoundaryConditions(model, Kg):
    nids = []
    allNodes = []
    for (eid, element) in model.elements.iteritems():
        nodes = element.nodeIDs()
        allNodes += nodes
        for nid in nodes:
            node = model.Node(nid)
            BCs = str(node.ps)
            if '1' in BCs:
                Kg = setRow(Kg, nid - 1, 0.)
                Kg = setCol(Kg, nid - 1, 0.)
            else:
                nids.append(nid - 1)

    notNids = [nid - 1 for nid in allNodes if nid not in nids]
    #print "notNids =\n",notNids
    print "nids = ", nids
    print "Kg2 =\n", ListPrint(Kg)
    Kr = reduce_matrix(Kg, nids)
    #Kr = reduce_matrix(Kg,notNids)
    print "Kreduced =\n", ListPrint(Kr)
    return Kr


def setRow(Kg, nid, value):
    #Kg[:][nid] = 10.
    nrows, ncols = Kg.shape
    for irow in xrange(nrows):
        Kg[irow, nid] = value
    return Kg


def setCol(Kg, nid, value):
    #Kg[nid][:] = 9.
    nrows, ncols = Kg.shape
    for icol in xrange(ncols):
        Kg[nid, icol] = value
    return Kg


def getForces(model, Dofs):
    Fvector = zeros(model.nNodes() * 3, 'd')
    print model.loads
    for loadSet, loads in model.loads.iteritems():
        ## @todo if loadset in required loadsets...
        #print loads
        for load in loads:
            if load.type == 'FORCE':
                loadDir = load.F()
                #nid = load.nodeID()
                for nodeID, F in loadDir.iteritems():
                    dof = Dofs[nodeID]
                    print "dof[%s] = %s" % (nodeID, dof)
                    Fvector[dof[0]] = F[0]
                    Fvector[dof[1]] = F[1]
                    Fvector[dof[2]] = F[2]
                    print "FVector = ", Fvector

    print "Fvector = ", Fvector
    return Fvector


def runTruss():
    model = makeTruss2()
    model.cross_reference()
    #for id,e in model.elements.iteritems():
    #    print "K = \n",e.Stiffness(model),'\n'

    Kglobal, Dofs = buildGlobalStiffness(model)
    #Kreduced = applyBoundaryConditions(model,Kglobal)

    F = getForces(model, Dofs)

    q = solveKF(model, Kglobal, F, Dofs)

    print q
    for eid, elem in model.elements.iteritems():
        elem.displacementStress(model, q)

    #x = solve(Kreduced,F) # deflections
    #print "x = ",x


def solveKF(model, Kg, F, Dofs):
    freeDOFs = {}
    constrainedDOFs = {}

    dofmap = {}
    i = 0
    for nid, node in model.nodes.iteritems():
        #constrained=[]; unconstrained=[]
        BCs = str(node.ps)
        dofs = Dofs[nid]

        constrainedDOFs[nid] = []
        freeDOFs[nid] = []
        print "BCs = ", BCs

        #dofmap[(nid,1)] = dofs[0]
        #dofmap[(nid,2)] = dofs[1]

        if '1' in BCs:
            constrainedDOFs[nid].append([dofs[0], i])
        else:
            freeDOFs[nid].append([dofs[0], i])
        i += 1

        if '2' in BCs:
            constrainedDOFs[nid].append([dofs[1], i])
        else:
            freeDOFs[nid].append([dofs[1], i])

        i += 1

        if '3' in BCs:
            constrainedDOFs[nid].append([dofs[1], i])
        else:
            freeDOFs[nid].append([dofs[1], i])

        i += 1

    nC = len(constrainedDOFs)
    nUC = len(freeDOFs)

    Kff = matrix(zeros((nUC, nUC), 'd'))  # K-free-free
    #Kfs = matrix(zeros((nUC,nC ),'d')))  # K-free-supported
    Qf = zeros(nUC, 'd')

    nDOF = 0
    allPacks = []
    for nodeID in model.nodes:
        aDOF = freeDOFs[nodeID]

        for dofPack in aDOF:
            (dofs, i) = dofPack
            print "dofs=%s i=%s" % (dofs, i)
            allPacks.append(i)
            #allDOFs.append(dofs)

    print allPacks

    for j, dof in enumerate(allPacks):
        for k, dof2 in enumerate(allPacks):
            print "Kff[%s][%s] = Kg[%s][%s] = %s" % (
                j + 1, k + 1, dof, dof2, Kg[dof, dof2])
            Kff[j, k] = Kg[dof, dof2]
        Qf[j] = F[dof]

    qf = solve(Kff, Qf)
    q = zeros(model.nNodes() * 3, 'd')

    print "Kff = \n%s" % (Kff)
    print "Qf  = %s" % (Qf)
    print "qf  = %s" % (qf)

    for dof, qi in izip(allPacks, qf):
        q[dof] = qi

    print "q   = %s" % (q)
    return q


def fKx(K, x):
    r"""
    \f[ {F} = [K]{x} \f]
    \f[ {x} = [K]^-1 {F} \f]
    """
    f = solve(K, x)
    return f


if __name__ == '__main__':
    runTruss()
