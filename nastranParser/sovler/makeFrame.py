from numpy import array,zeros,matrix

from nastranParser.bdf.FEM_Mesh import FEM_Mesh
from nastranParser.general.general import ListPrint
from nastranParser.general.generalMath import reduceMatrix


def makeTruss():
    p1 = array([ 0., 0., 0.]) *12.
    p2 = array([10., 0., 0.]) *12.
    p3 = array([10.,10., 0.]) *12.
    p4 = array([ 0.,10., 0.]) *12.
    
    p5 = array([ 0., 0.,10.]) *12.
    p6 = array([10., 0.,10.]) *12.
    p7 = array([10.,10.,10.]) *12.
    p8 = array([ 0.,10.,10.]) *12.
    
    mesh = FEM_Mesh('')
    nodes = [p1,p2,p3,p4,p5,p6,p7,p8]
    for i,node in enumerate(nodes):
        spc = None
        if node[2]==0.:
            spc = str(123456)
        grid = ['GRID',i+1,None,node[0],node[1],node[2],None,spc]
        mesh.addCard(grid,'GRID')


    Ac  = 4. # in^2
    Jc  = 60. # in4
    Ixc = 650. # in4
    Iyc = 65. # in4
    Ec  = 29000.
    
    Ab  = 3.2
    Jb  = 43.
    Ixb = 450.
    Iyb = 32.
    Eb  = 29000.
    
          #id
    c1 = [5,1]
    c2 = [6,2]
    c3 = [7,3]
    c4 = [8,4]
    c = [c1,c2,c3,c4]

    b1 = [1,2]
    b2 = [2,3]
    b3 = [3,4]
    b4 = [4,1]
    
    b5 = [5,6]
    b6 = [6,7]
    b7 = [7,8]
    b8 = [8,5]
    b = [b1,b2,b3,b4,b6,b7,b8]

    mid1 = 99
    mat1 = ['MAT1',mid1,Ec,None,0.3]
    mesh.addCard(mat1,'MAT1')
    for i,rod in enumerate(c):
        #print "e[%s] = %s" %(i,rod)
        (n1,n2) = rod
        conrod = ['CONROD',i+1,n1,n2,mid1,Ac,Jc]
        mesh.addCard(conrod,'CONROD')

    mid2 = 22
    mat1 = ['MAT1',mid2,Eb,None,0.3]
    mesh.addCard(mat1,'MAT1')
    for i,rod in enumerate(c):
        (n1,n2) = rod
        #print "rod = ",rod
        conrod = ['CONROD',i+1,n1,n2,mid2,Ab,Jb]
        mesh.addCard(conrod,'CONROD')
        #print conrod
        Ke = mesh.Element(i+1).Stiffness(mesh)
    
    scaleFactor = 1.
    loadID = 123
    load = ['LOAD',loadID,scaleFactor,1.,loadID]
    mesh.addCard(load,'LOAD')

    force = ['FORCE',loadID,3,None,100.,  1.,0.,0.]
    mesh.addCard(force,'FORCE')
    mesh.write('frame.bdf')
    print "done"
    return mesh

def doProblem(elements):
    Kg = buildGlobalStiffness()     # global  stiffness (Kg)
    Kr = applyBoundaryConditions()  # reduced stiffness (Kr)

    #F = K*x
    F = getForces()
    x = solve(Kr,F) # deflections
    
    (stress,strain) = recoverStressStrain(elements,x)

def buildGlobalStiffness(mesh):
    nElements = len(mesh.elements)+4
    Kg = matrix(zeros((nElements,nElements),'d') ) # K_global
    
    for id,element in sorted(mesh.elements.items()):
        #nodes = element.nodes
        print element
        Ke = element.Stiffness(mesh)
        #print "K_element[%s] = \n%s\n" %(id,Ke)

        nodes = element.getNodeIDs()

        # put in global stiffness matrix
        for (i,iNode) in enumerate(nodes): # row
            for (j,jNode) in enumerate(nodes): # column
                print "Kg[%s][%s] = Ke[%s][%s] = %s" %(iNode,jNode,i,j,Ke[i,j])

                Kg[iNode-1,jNode-1] = Ke[i,j]
            ###
        ###
    ###
    print "K_global =\n",ListPrint(Kg)
    #print "K_global = \n",Kg
    return Kg

# B - strain displacement matrix
# Uo = 1/2*d.T*B.T*E*B
# Ke = integral(B.T*E*B,dV)

def applyBoundaryConditions(mesh,Kg):
    #for (nid,node) in sorted(mesh.nodes.items()):
    nids = []
    allNodes = []
    for (eid,element) in mesh.elements.items():
        nodes = element.getNodeIDs()
        allNodes+=nodes
        for nid in nodes:
            node = mesh.Node(nid)
            BCs = str(node.ps)
            if '1' in BCs:
                Kg = setRow(Kg,nid-1,0.)
                Kg = setCol(Kg,nid-1,0.)
            else:
                nids.append(nid-1)
        ###
    ###
    notNids = [nid-1 for nid in allNodes if nid not in nids]
    #print "notNids =\n",notNids
    print "nids = ",nids
    print "Kg2 =\n",ListPrint(Kg)
    Kr = reduceMatrix(Kg,nids)
    #Kr = reduceMatrix(Kg,notNids)
    print "Kreduced =\n",ListPrint(Kr)
    return Kr

def setRow(Kg,nid,value):
    #Kg[:][nid] = 10.
    nrows,ncols = Kg.shape
    for irow in range(nrows):
        Kg[irow,nid] = value
    return Kg

def setCol(Kg,nid,value):
    #Kg[nid][:] = 9.
    nrows,ncols = Kg.shape
    for icol in range(ncols):
        Kg[nid,icol] = value
    return Kg

def runTruss():
    mesh = makeTruss()

    #for id,e in mesh.elements.items():
    #    print "K = \n",e.Stiffness(mesh),'\n'
    
    Kglobal  = buildGlobalStiffness(mesh)
    Kreduced = applyBoundaryConditions(mesh,Kglobal)
    

def fKx(K,x):
    """
    {F} = [K]{x}
    {x} = [K]^-1 {F}
    """
    f = solve(K,x)
    return f


runTruss()
