# pylint: disable=E0602,C0103
from __future__ import print_function
import sys
from itertools import izip
from numpy import array, zeros, ones
from numpy.linalg import solve

from pyNastran.general.mathematics import print_matrix, print_annotated_matrix
from pyNastran.bdf.bdf import BDF, SPC, SPC1
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2

from pyNastran.op2.tables.oug.oug_displacements import DisplacementObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject


def partition_sparse(Is, Js, Vs):
    I2 = []
    J2 = []
    V2 = []
    for (i, j, v) in (Is, Js, Vs):
        if abs(v) >= 1e-8:
            I2.append(i)
            J2.append(j)
            V2.append(v)
    return(I2, J2, V2)


def getDOF_Set(nAll, dofs):
    dofsAll = set([i for i in xrange(nAll)])
    dofs = list(dofsAll.difference(set(dofs)))
    return dofs


def partition_dense_symmetric(A, dofs):
    nAll = A.shape[0]
    dofs = getDOF_Set(nAll, dofs)
    dofs.sort()
    n = len(dofs)
    A2 = zeros((n, n), 'float64')
    for (i, dofI) in enumerate(dofs):
        for (j, dofJ) in enumerate(dofs):
            v = A[dofI, dofJ]
            if abs(v) >= 1e-8:
                A2[i, j] = v
    return(A2)


def partition_dense_vector(F, dofs):
    nAll = F.shape[0]
    #print("dofs = ", dofs)
    dofs = getDOF_Set(nAll, dofs)
    dofs.sort()
    #print("dofs = ", dofs)
    n = len(dofs)
    F2 = zeros(n, 'float64')
    for (i, dofI) in enumerate(dofs):
        v = F[dofI]
        if abs(v) >= 1e-8:
            F2[i] = v
    return(F2)


def partition_sparse_vector(F, dofs):
    dofs.sort()
    #n = len(dofs)
    F2i = []
    F2v = []
    for (i, dofI) in enumerate(dofs):
        v = F2v[dofI]
        if abs(v) >= 1e-8:
            F2i.append(i)
            F2v.append(v)
    return(F2i, F2v)


def departition_dense_vector(n, IsVs):
    V = zeros(n)
    for IV in IVs:
        (Is, Vs) = IV
        for (i, v) in izip(Is, Vs):
            V[i] = v
    return(V)


def reverseDict(A):
    B = {}
    for (key, value) in A.iteritems():
        B[value] = key
    return B


class Solver(F06, OP2):
    """
    Solves SOL 101

    Progress:
      - Solves a 2D problem using matrix partitioning.of the SPC set
      - 3D works, but is disabled for testing

    TODO:
      - Need to combine solved displacements with original known
        displacements to create displacement set
      - Calculate Stress/Strain
      - Write the OP2

    Case Control
      LOAD,SPC,TITLE,LABEL

    Bulk Data:
      GRID,CORDx
      CONROD, CROD, PROD
      MAT1
      LOAD, FORCE
      SPC, SPC1

    Results:
      @todo DISPLACEMENT solver results
      @todo STRESS solver results
      @todo STRAIN solver results
    """
    def __init__(self):
        F06.__initAlt__(self)
        OP2.__objectsInit__(self)
        self.nU = 0
        self.nUs = 0
        self.nUm = 0

        # displacments
        self.U = []
        self.Us = []
        self.Um = []

        # indicies in U that correspond to Us/Um
        self.iUs = []
        self.iUm = []

    def solve(self, K, F):  # can be overwritten
        r"""solves \f$ [K]{x} = {F}\f$ for \f${x}\f$"""
        print("--------------")
        print("Kaa_norm = \n" + str(K / 250000.))
        print("--------------")
        print("Fa = ", F)
        assert max(F) != min(F), 'no load is applied...'
        print("--------------")

        #asdf
        return solve(K, F)

    def run(self, bdfName):
        model = BDF()
        model.cardsToRead = get_cards()
        model.readBDF(bdfName)
        cc = model.caseControlDeck
        #print cc.subcases
        analysisCases = []
        for (isub, subcase) in sorted(cc.subcases.iteritems()):
            if subcase.hasParameter('LOAD'):
                analysisCases.append(subcase)

        #print analysisCases
        for case in analysisCases:
            print(case)
            (value, options) = case.get_parameter('STRESS')
            print("STRESS value   = %s" % (value))
            print("STRESS options = %s" % (options))

            if case.hasParameter('TEMPERATURE(INITIAL)'):
                (value, options) = case.get_parameter('TEMPERATURE(INITIAL)')
                print('value   = %s' % (value))
                print('options = %s' % (options))
                raise NotImplementedError('TEMPERATURE(INITIAL) not supported')
                #integrate(B.T*E*alpha*dt*Ads)
            #sys.exit('starting case')
            self.runCase(model, case)

    def runCase(self, model, case):
        sols = {101: self.runSOL101}

        iSubcase = case.id
        if model.sol in sols:
            if case.hasParameter('TITLE'):
                (self.Title, options) = case.get_parameter('TITLE')
            else:
                self.Title = 'pyNastran Default Title'
            if case.hasParameter('SUBTITLE'):
                (self.Subtitle, options) = case.get_parameter('SUBTITLE')
            else:
                self.Subtitle = 'DEFAULT'
            self.iSubcaseNameMap[iSubcase] = [self.Title, self.Subtitle]

            sols[model.sol](model, case)
        else:
            raise NotImplementedError('model.sol=%s not in %s' %
                                      (model.sol, sols.keys()))

    def buildNidComponentToID(self, model):
        i = 0
        nidComponentToID = {}
        for (nid, node) in sorted(model.nodes.iteritems()):  # GRIDs
            for ps in node.ps:
                if self.is3D or ps in ['1', '2']:
                    self.iUs.append(i + int(ps) - 1)
                    self.Us.append(0.0)
                elif ps == '5' and not(self.is3D):  # 4 or 5
                    self.iUs.append(i + 2)
                    self.Us.append(0.0)

            if self.is3D:
                nidComponentToID[(nid, 1)] = i
                nidComponentToID[(nid, 2)] = i + 1
                nidComponentToID[(nid, 3)] = i + 2
                nidComponentToID[(nid, 4)] = i + 3
                nidComponentToID[(nid, 5)] = i + 4
                nidComponentToID[(nid, 6)] = i + 5
                i += 6
            else:
                nidComponentToID[(nid, 1)] = i
                nidComponentToID[(nid, 2)] = i + 1
               #nidComponentToID[(nid, 4)] = i+2 # torsion
                nidComponentToID[(nid, 5)] = i + \
                    2  # bending - can pick one moment
                i += 3
            #print('iUs[%i] = %s' %(nid,self.iUs))

        if model.spoints:
            for nid in sorted(model.spoints.spoints):  # SPOINTS
                nidComponentToID[(nid, 1)] = i
                i += 1
        return(nidComponentToID, i)

    def build_Kgg_Fg(self, model, case, nidComponentToID, i):
        (isSPC) = self.applySPCs(model, case, nidComponentToID)
        (isMPC) = self.applyMPCs(model, case, nidComponentToID)

        #spcDOFs = self.iUs
        #mpcDOFs = self.iUm

        Ug = ones(i)
        Fg = zeros(i, 'float64')

        Kgg = zeros((i, i), 'float64')
        #Mgg = zeros((i, i))

        Fg = self.assembleForces(model, case, Fg, nidComponentToID)
        Kgg = self.assembleGlobalStiffness(model, Kgg, nidComponentToID, Fg)
        return(Kgg, Fg, isSPC, isMPC)

    def runSOL101(self, model, case):
        #print("case = ", case)
        assert model.sol == 101, 'model.sol=%s is not 101' % (model.sol)

        ## define IDs of grid point components in matrices

        self.is3D = True
        #self.is3D = False

        (Ua, n) = self.setupSOL101(model, case)
        print("------------------------\n")
        print("Ua = ", Ua)
        print("Us = ", self.Us)

        dofsA = getDOF_Set(n, self.iUs)
        dofsA.sort()
        U = zeros(n, 'float64')
        print("U = ", U)
        print("iUs   = ", self.iUs)

        ## @todo handle MPCs
        for (i, iu) in enumerate(self.iUs):
            U[iu] = self.Us[i]
        for (i, iu) in enumerate(dofsA):
            U[iu] = Ua[i]
        print("*U = ", U)
        print("dofsA = ", dofsA)

        self.storeDisplacements(model, U, case)
        self.writeDisplacements(case)
        q = U
        for eid, elem in model.elements.iteritems():
            elem.displacementStress(model, q, self.nidComponentToID)

    def storeDisplacements(self, model, U, case):
        """
        fills the displacement object
        """
        self.iSubcases = []
        #self.log = None
        analysisCode = 1
        transient = False
        iSubcase = case.id
        isSort1 = False
        dt = None

        dataCode = {'log': self.log, 'analysisCode': analysisCode,
                    'deviceCode': 1, 'tableCode': 1, 'sortCode': 0,
                    'sortBits': [0, 0, 0], 'numWide': 8, 'tableName': 'OUG',
                    'nonlinearFactor': None}
        #print "headers = %s" %(headers)

        disp = DisplacementObject(dataCode, isSort1, iSubcase, dt=False)

        data = []

        i = 0
        #(nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
        for (nid, node) in sorted(model.nodes.iteritems()):
            line = [nid]
            if node.type == 'GRID':
                line.append('G')
            else:
                raise NotImplementedError('node.type=%s' % (node.type))
            if self.is3D:
                line += U[i:i + 6]  # 1,2,3,4,5,6
                i += 6
            else:
                line += [U[i], U[i + 1], 0., 0., U[i + 2], 0.]  # 1,2,5
                i += 3
            #print("line = ",line)
            data.append(line)
        disp.addF06Data(data, dt)
        self.displacements[iSubcase] = disp
        self.iSubcases.append(iSubcase)

    def writeDisplacements(self, case):
        self.writeF06('2dTruss.f06')

    def setupSOL101(self, model, case):
        # the (GridID,componentID) -> internalID
        (self.nidComponentToID, i) = self.buildNidComponentToID(model)
        (Kgg, Fg, isSPC, isMPC) = self.build_Kgg_Fg(model, case,
                                                    self.nidComponentToID, i)

        (self.IDtoNidComponents) = reverseDict(self.nidComponentToID)
        print("IDtoNidComponents = ", self.IDtoNidComponents)
        print("Kgg =\n" + print_annotated_matrix(Kgg, self.IDtoNidComponents))
        #print("Kgg = \n", Kgg)
        #print("iSize = ", i)

        #(Kaa,Fa) = self.Partition(Kgg)
        #sys.exit('verify Kgg')

        Kaa = partition_dense_symmetric(Kgg, self.iUs)
        print("Kaa = \n%s" % (print_matrix(Kaa)))
        #print("Kaa.shape = ",Kaa.shape)

        #sys.exit('verify Kaa')
        Fa = partition_dense_vector(Fg, self.iUs)
        #print("Kaa = \n%s" %(print_matrix(Kaa)))

        print("Fg = ", Fg)
        print("Fa = ", Fa)
        print("Us = ", self.Us)

        self.Us = array(self.Us, 'float64')
        self.Um = array(self.Um, 'float64')

        zero = array([])
        MPCgg = zero
        Ksa = Kas = Kss = Cam = Cma = Kma = Kam = Kmm = Kms = Ksm = zero
        Kaa1 = Kaa2 = zero
        Fm = Fs = zero

        #Kaa = partition_dense_matrix(Kgg,iFree)
        Kaa0 = Kaa
        Fa0 = Fa

        if isSPC:
           #Fs  = partition_dense_vector(Fg,self.iUs)
            Ksa = partition_dense_matrix(Kgg, self.iUs, iFree)
            Kas = Ksa.transpose()
            Kss = partition_dense_matrix(Kgg, self.iUs)

        if isMPC:
            Fm = partition_dense_vector(Fg, self.iUm)
            Cam = partition_dense_matrix(MPCgg, iFree)
            Cma = partition_dense_matrix(MPCgg, self.iUm)

            Kaa1 = Cam * Kmm * Cma
            Kaa2 = Cam * Kma + Kam * Cma
            assert Cam.transpose() == Cma

            Kma = partition_dense_matrix(Kgg, self.iUm, iFree)
            Kam = Kma.transpose()
            Kmm = partition_dense_matrix(Kgg, self.iUm)
            if isSPC:
                Kms = partition_dense_matrix(Kgg, self.iUm, self.iUs)
                Ksm = Kms.transpose()

        Fa = Fa0  # + Cam*Fm
        Kaa = Kaa0  # +Kaa1+Kaa2

        Ua = self.solve(Kaa, Fa)
        #self.Um = Kma*Ua

        return(Ua, i)

    def assembleGlobalStiffness(self, model, Kgg, Dofs, F):
        for (eid, elem) in sorted(model.elements.iteritems()):  # CROD, CONROD

            # nIJV is the position of the values of K in the dof
            (K, nIJV, Fg, nGrav) = elem.Stiffness(model, self.gravLoad,
                                                  self.is3D)
            print("K[%s] = \n%s" % (eid, K))
            (Ki, Kj) = K.shape
            ij = 0
            nij2 = []
            for ij in nIJV:
                nij2.append(Dofs[ij])
            print('nij2', nij2)

            for (fg, dof) in izip(Fg, nGrav):
                #print("dof = ",dof)
                if dof in Dofs:  # is3D
                    F[Dofs[dof]] += fg

            for i in xrange(Ki):
                for j in xrange(Kj):
                    kij = K[i, j]
                    #if abs(kij)>1e-8:
                    #print('niJV = ',nIJV[ij])
                    ii = nij2[i]
                    jj = nij2[j]
                    #dof = Dofs[n]
                    Kgg[ii, jj] += kij

                    #ij += 1

        if 0:
            # n is (nid,componentID), IJV is the (ith,jth,value) in K
            #(n,IJV) = elem.nIJV()
            for (ni, ijv) in izip(n, IJV):
                i = nidComponentToID(ni)
                j = nidComponentToID(ji)
                #(i,j,v) = ijv
                Kgg[i, j] = v
                #KggI.append(i)
                #KggJ.append(j)
                #KggV.append(v)

        return Kgg

    def applySPCs2(self, model, case, nidComponentToID):
        if case.hasParameter('SPC'):
            # get the value, 1 is the options (SPC has no options)
            spcID = case.get_parameter('SPC')[0]
            SpcSet = model.SPC(spcID)

            for spcSet in SpcSet:
                (typesFound, positionSPCs) = spcSet.organizeConstraints(model)
        return

    def applySPCs(self, model, case, nidComponentToID):
        return self.applySPCs2(model, case, nidComponentToID)
        isSPC = False
        print('*Us', self.Us)
        print('*iUs', self.iUs)
        if case.hasParameter('SPC'):
            isSPC = True
            spcs = model.spcObject2.constraints
            # get the value, 1 is the options (SPC has no options)
            spcID = case.get_parameter('SPC')[0]
            print("SPC = ", spcID)
            #print model.spcObject2.constraints
            spcset = spcs[spcID]

            for spc in spcset:
                print(spc)
                if isinstance(spc, SPC):
                    ps = spc.constraints
                    print("spc.constraints = ", spc.constraints)
                    print("spc.enforced    = ", spc.enforced)
                    raise NotImplementedError('no support for SPC...')
                    self.Us.append(self.enforced)
                    for ips in ps:
                        key = (nid, int(ips))
                        i = nidComponentToID(key)
                        self.iUs.append(i)
                        self.Us.append(0.0)

                elif isinstance(spc, SPC1):
                    ps = spc.constraints
                    #print "ps = |%s|" %(ps)
                    nodes = spc.nodes
                    for nid in nodes:
                        for ips in ps:
                            ips = int(ips)
                            if not(self.is3D) and ips not in ['1', '2', '5']:
                                continue

                            key = (nid, ips)
                            i = nidComponentToID[key]
                            #print "i=%s Us=%s" %(i,0.0)
                            if i not in self.iUs:
                                self.iUs.append(i)
                                self.Us.append(0.0)

                else:
                    raise NotImplementedError('Invalid SPC...spc=\n%s' % (spc))

            print("iUs = ", self.iUs)
            print("Us  = ", self.Us)
            sys.exit('stopping in applySPCs')

        print("iUs = ", self.iUs)
        print("Us  = ", self.Us)
        return (isSPC)

    def applyMPCs(self, model, case, nidComponentToID):
        isMPC = False
        if case.hasParameter('MPC'):
            isMPC = True
            mpcs = model.mpcObject2.constraints
            # get the value, 1 is the options (MPC has no options)
            mpcID = case.get_parameter('MPC')[0]
            print("******")
            print(model.mpcObject2.constraints)
            print("mpcID = ", mpcID)
            mpcset = mpcs[mpcID]

            for mpc in mpcset:
                print(mpc, type(mpc))

            msg = 'MPCs are not supported...stopping in applyMPCs'
            raise NotImplementedError(msg)
        return (isMPC)

    def assembleForces(self, model, case, Fg, Dofs):
        """very similar to writeCodeAster loads"""
        #print(model.loads)
        (loadID, junk) = model.caseControlDeck.get_subcase_parameter(case.id,
                                                                   'LOAD')
        print("loadID = ", loadID)
        LoadSet = model.Load(loadID)

        self.gravLoad = array([0., 0., 0.])
        for loadSet in LoadSet:
            print(loadSet)
            (typesFound, forceLoads, momentLoads,
                         forceConstraints, momentConstraints,
                         gravityLoad) = loadSet.organizeLoads(model)

            nids = []
            for nid in forceLoads:
                nids.append(nid)
            for nid in momentLoads:
                nids.append(nid)

            if gravityLoad != []:
                print("gravityLoad = ", gravityLoad)
                self.gravLoad += gravityLoad
            for nid in sorted(nids):
                #print("nid = ", nid)
                #print("load = ", load)

                if nid in forceLoads:
                    force = forceLoads[nid]
                    if abs(force[0]) > 0.:
                        Fg[Dofs[(nid, 1)]] += force[0]
                    if abs(force[1]) > 0.:
                        Fg[Dofs[(nid, 2)]] += force[1]
                    if abs(force[2]) > 0. and self.is3D:
                        Fg[Dofs[(nid, 3)]] += force[2]

                if nid in momentLoads:
                    moment = momentLoads[nid]
                    if abs(moment[0]) > 0. and self.is3D:
                        Fg[Dofs[(nid, 4)]] += moment[2]
                    if abs(moment[1]) > 0.:
                        Fg[Dofs[(nid, 5)]] += moment[2]
                    if abs(moment[2]) > 0. and self.is3D:
                        Fg[Dofs[(nid, 6)]] += moment[2]

        if not self.is3D:
            self.gravLoad = array([self.gravLoad[0], self.gravLoad[1]])
        print("Fg  = ", Fg)
        return Fg

    def writeResults(self, case):
        Us = self.Us
        Um = self.Um
        Ua = self.Ua

        iUs = self.iUs
        iUm = self.iUm
        iUa = self.iUa
        pageNum = 1

        if case.hasParameter('DISPLACEMENT'):
            (value, options) = case.get_parameter('DISPLACEMENT')
            if options is not []:
                UgSeparate = [[Ua, iUa], [Us, iUs], [Um, iUm]]
                Ug = departition_dense_vector(UgSeparate)

                result = DisplacementObject(dataCode, transient)
                result.addF06Data()

                if 'PRINT' in options:
                    f06.write(result.writeF06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(self.Title, self.Subtitle))

        if case.hasParameter('SPCFORCES'):
            (value, options) = case.get_parameter('SPCFORCES')
            if options is not []:
                SPCForces = Ksa * Ua + Kss * Us
                if isMPC:
                    SPCForces += Ksm * Um

                result = SPCForcesObject(dataCode, transient)
                result.addF06Data()
                if 'PRINT' in options:
                    f06.write(result.writeF06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(Title, Subtitle))

        if case.hasParameter('MPCFORCES'):
            if options is not []:
                (value, options) = case.get_parameter('MPCFORCES')
                MPCForces = Kma * Ua + Kmm * Um
                if isSPC:
                    MPCForces += Kms * Us

                result = MPCForcesObject(dataCode, transient)
                result.addF06Data()
                if 'PRINT' in options:
                    f06.write(result.writeF06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    f06.write(result.writeOP2(Title, Subtitle))

        if case.hasParameter('GPFORCE'):
            if options is not []:
                (value, options) = case.get_parameter('GPFORCE')
                AppliedLoads = Kaa * Ua
                if isSPC:
                    AppliedLoads += Kas * Us
                if isMPC:
                    AppliedLoads += Kam * Um

                result = AppliedLoadsObject(dataCode, transient)
                result.addF06Data()
                if 'PRINT' in options:
                    f06.write(result.writeF06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(Title, Subtitle))

        if case.hasParameter('STRAIN'):
            if options is not []:
                (value, options) = case.get_parameter('STRAIN')

                for (eid, elem) in sorted(model.elements()):
                    pass

                result = xxxObject(dataCode, transient)
                result.addF06Data()
                if 'PRINT' in options:
                    f06.write(result.writeF06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(Title, Subtitle))


def get_cards():
    cardsToRead = set([
                      'PARAM',
                      'GRID', 'GRDSET',

                      # elements
                      'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
                      'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',

                      'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CONROD', 'CBEND',
                      'CTRIA3', 'CTRIA6',
                      'CQUAD4', 'CQUAD8',
                      'CTETRA', 'CPENTA', 'CHEXA',
                      'CSHEAR',

                      # rigid elements - represent as MPCs???
                      #'RBAR','RBAR1','RBE1','RBE2','RBE3',

                      # properties
                      'PELAS',
                      'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PTUBE', 'PBEND',
                      'PSHELL', 'PCOMP', 'PSHEAR',  # 'PCOMPG',
                      'PSOLID',

                      # materials
                      'MAT1', 'MAT2', 'MAT8',

                      # spc/mpc constraints
                      'SPC', 'SPC1', 'SPCADD',
                      #'MPC','MPCADD',

                      # loads
                      'LOAD',
                      'FORCE', 'FORCE1', 'FORCE2',
                      'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
                      'MOMENT', 'MOMENT1', 'MOMENT2',
                      'GRAV',

                      # coords
                      'CORD1R', 'CORD1C', 'CORD1S',
                      'CORD2R', 'CORD2C', 'CORD2S',

                      # other
                      'INCLUDE',
                      'ENDDATA',
                      ])
    return cardsToRead

if __name__ == '__main__':
    bdfName = sys.argv[1]
    s = Solver()
    s.run(bdfName)
