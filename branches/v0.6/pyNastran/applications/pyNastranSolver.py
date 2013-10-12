# pylint: disable=E0602,C0103
from __future__ import print_function
import os
import sys
from itertools import izip

# 3rd party
from numpy import array, zeros, ones
from numpy import searchsorted
from numpy.linalg import solve

# pyNastran
from pyNastran.utils.mathematics import print_matrix, print_annotated_matrix
from pyNastran.bdf.bdf import BDF, SPC, SPC1
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2

from pyNastran.f06.f06Writer import make_grid_point_singularity_table

# Tables
from pyNastran.op2.tables.oug.oug_displacements import DisplacementObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject

# rods
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RodStressObject, RodStrainObject
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealRodForce
RealRodForce

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
    print("dofs = ", dofs)
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
    for IV in IsVs:
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
        OP2.__init__(self, '')

        self.is3D = True

        # normalization of stiffness matrix
        self.fnorm = 1.0

        self.iSubcases = []
        self.nU = 0
        self.nUs = 0
        self.nUm = 0

        # displacements
        self.U = []
        self.Us = []
        self.Um = []

        # indicies in U that correspond to Us/Um
        self.iUs = []
        self.iUm = []

        self.case_result_flags = {}

    def solve(self, K, F):  # can be overwritten
        r"""solves \f$ [K]{x} = {F}\f$ for \f${x}\f$"""
        print("--------------")
        print("Kaa_norm = \n" + str(K / self.fnorm))
        print("--------------")
        print("Fa = ", F)
        assert max(F) != min(F), 'no load is applied...'
        print("--------------")

        #asdf
        return solve(K, F)

    def run_solver(self, bdfName):
        bdfName = os.path.abspath(bdfName)
        bdf_base, ext = os.path.splitext(bdfName)
        self.f06_name = bdf_base + '.f06'
        self.op2_name = bdf_base + '.op2'
        print("**is3D =", self.is3D)

        model = BDF()
        model.cardsToRead = get_cards()
        model.readBDF(bdfName)
        cc = model.caseControlDeck
        #print cc.subcases
        analysisCases = []
        for (isub, subcase) in sorted(cc.subcases.iteritems()):
            if subcase.has_parameter('LOAD'):
                analysisCases.append(subcase)

        #print analysisCases
        for case in analysisCases:
            print(case)
            (value, options) = case.get_parameter('STRESS')
            print("STRESS value   = %s" % (value))
            print("STRESS options = %s" % (options))

            if case.has_parameter('TEMPERATURE(INITIAL)'):
                (value, options) = case.get_parameter('TEMPERATURE(INITIAL)')
                print('value   = %s' % (value))
                print('options = %s' % (options))
                raise NotImplementedError('TEMPERATURE(INITIAL) not supported')
                #integrate(B.T*E*alpha*dt*Ads)
            #sys.exit('starting case')
            self.run_case(model, case)

    def run_case(self, model, case):
        sols = {101: self.run_sol_101}

        isubcase = case.id
        if model.sol in sols:
            if case.has_parameter('TITLE'):
                (self.Title, options) = case.get_parameter('TITLE')
            else:
                self.Title = 'pyNastran Default Title'
            if case.has_parameter('SUBTITLE'):
                (self.Subtitle, options) = case.get_parameter('SUBTITLE')
            else:
                self.Subtitle = 'DEFAULT'
            self.iSubcaseNameMap[isubcase] = [self.Title, self.Subtitle]

            # really should be is_f06_stress, is_op2_stress, etc.
            # also should have SET support
            #print(case.params)
            #if case.has_parameter('DISPLACEMENT'):
                #value, options = case.get_parameter('DISPLACEMENT')
                #if value == 'ALL':
                    #self.is_displacement = True
                #elif value == 'NONE':
                    #self.is_displacement = False
                #else:
                    #raise NotImplementedError('DISPLACEMENT = %r is not supported' % value)

            self.get_case_parameter(case, 'DISPLACEMENT')
            self.get_case_parameter(case, 'STRESS')
            self.get_case_parameter(case, 'STRAIN')
            self.get_case_parameter(case, 'FORCE')

            #if case.has_parameter('STRESS'):
                #value, options = case.get_parameter('STRESS')
                #if value == 'ALL':
                    #self.is_stress = True
                #elif value == 'NONE':
                    #self.is_stress = False
                #else:
                    #raise NotImplementedError('STRESS = %r is not supported' % value)

            #if case.has_parameter('STRAIN'):
                #value, options = case.get_parameter('STRAIN')
                #if value == 'ALL':
                    #self.is_strain = True
                #elif value == 'NONE':
                    #self.is_strain = False
                #else:
                    #raise NotImplementedError('STRAIN = %r is not supported' % value)

            #if case.has_parameter('FORCE'):
                #value, options = case.get_parameter('FORCE')
                #if value == 'ALL':
                    #self.is_force = True
                #elif value == 'NONE':
                    #self.is_force = False
                #else:
                    #raise NotImplementedError('FORCE = %r is not supported' % value)

            self.is_displacement = True
            self.is_stress = True
            self.is_strain = True
            self.is_force = True
            if not self.case_result_flags['DISPLACEMENT'][0]:
                self.is_displacement = False
            if not self.case_result_flags['STRESS'][0]:
                self.is_stress = False
            if not self.case_result_flags['STRAIN'][0]:
                self.is_strain = False
            if not self.case_result_flags['FORCE'][0]:
                self.is_force = False

            if not(self.is_displacement or self.is_stress or self.is_strain or self.is_force):
                msg = 'No results selected...'
                raise RuntimeError(msg)
            sols[model.sol](model, case)
        else:
            raise NotImplementedError('model.sol=%s not in %s' % (model.sol, sols.keys()))

    def get_case_parameter(self, case, param_name):
        if case.has_parameter(param_name):
            value, options = case.get_parameter(param_name)
            if value == 'ALL':
                save_results = True
                is_bool = True
            elif value == 'NONE':
                save_results = False
                is_bool = True
            else:
                save_results = True
                is_bool = False
                raise NotImplementedError('%s = %r is not supported' % (param_name, value))
        self.case_result_flags[param_name] = [save_results, is_bool]

    def build_nid_component_to_id(self, model):
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
            #sys.exit('3D %r' % self.is3D)
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
               #nidComponentToID[(nid, 4)] = i + 2 # torsion
                nidComponentToID[(nid, 5)] = i + 2 # bending - can pick one moment
                i += 3
            #print('iUs[%i] = %s' % (nid, self.iUs))

        if model.spoints:
            for nid in sorted(model.spoints.spoints):  # SPOINTS
                nidComponentToID[(nid, 1)] = i
                i += 1
        return(nidComponentToID, i)

    def build_Kgg_Fg(self, model, case, nidComponentToID, i):
        (isSPC) = self.apply_SPCs(model, case, nidComponentToID)
        (isMPC) = self.apply_MPCs(model, case, nidComponentToID)

        #spcDOFs = self.iUs
        #mpcDOFs = self.iUm

        Ug = ones(i)
        Fg = zeros(i, 'float64')

        Kgg = zeros((i, i), 'float64')
        #Mgg = zeros((i, i))

        Fg = self.assemble_forces(model, case, Fg, nidComponentToID)
        Kgg = self.assemble_global_stiffness(model, i, Kgg, nidComponentToID, Fg)
        return(Kgg, Fg, isSPC, isMPC)

    def run_sol_101(self, model, case):
        #print("case = ", case)
        assert model.sol == 101, 'model.sol=%s is not 101' % (model.sol)

        ## define IDs of grid point components in matrices
        if 1:
            # analysis
            (Ua, n) = self.setup_sol_101(model, case)
            print("------------------------\n")
            print("Ua = ", Ua)
            print("Us = ", self.Us)

            dofsA = getDOF_Set(n, self.iUs)
            dofsA.sort()
            U = zeros(n, 'float64')
            print("U = ", U)
            print("iUs   = ", self.iUs)

            # TODO handle MPCs
            for (i, iu) in enumerate(self.iUs):
                U[iu] = self.Us[i]
            for (i, iu) in enumerate(dofsA):
                U[iu] = Ua[i]
            print("*U = ", U)
            print("dofsA = ", dofsA)

            if self.is_displacement:
                self.store_displacements(model, U, case)
            q = U
        else:
            n = len(model.nodes)
            q = ones(n, 'float64')

        # results
        rods = []
        beams = []
        bars = []
        quads = []
        tris = []
        type_map = {
            'CONROD' : rods,
            'CROD'   : rods,
            'CBEAM'  : beams,
            'CBAR'  : bars,
        }
        elements = model.elements
        for eid, element in elements.iteritems():
            try:
                a = type_map[element.type]
            except KeyError:
                raise NotImplementedError(element.type)
            a.append(eid)

        #=========================
        nnodes = len(model.nodes)
        nrods = len(rods)
        nbars = len(bars)
        nbeams = len(beams) # half implemented

        # not implemented
        nbars = len(bars)
        nquads = len(quads)
        ntris = len(tris)

        #=========================
        rods = array(rods)
        #=========================
        # RODS
        if self.is_stress or self.is_strain or self.is_force:
            if nrods:
                ox = zeros(nrods, 'float64')
                ex = zeros(nrods, 'float64')
                fx = zeros(nrods, 'float64')
                for i, eid in enumerate(rods):
                    element = elements[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID, is3D=self.is3D)
                    ox[i] = oxi
                    ex[i] = exi
                    fx[i] = fxi

                if self.is_strain:
                    self.store_rod_oes(model, rods, ex, case, Type='strain')
                if self.is_stress:
                    self.store_rod_oes(model, rods, ox, case, Type='stress')
                if self.is_force:
                    self.store_rod_oef(model, rods, fx, case)

            if nbeams:
                ox = zeros(nbeams, 'float64')
                ex = zeros(nbeams, 'float64')
                fx = zeros(nbeams, 'float64')
                for i, eid in enumerate(beams):
                    element = elements[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID, is3D=self.is3D)
                    ox[i] = oxi
                    ex[i] = exi
                    fx[i] = fxi

                if self.is_strain:
                    self.store_beam_oes(model, beams, ex, case, Type='strain')
                if self.is_stress:
                    self.store_beam_oes(model, beams, ox, case, Type='stress')
                if self.is_force:
                    self.store_beam_oef(model, beams, fx, case)

        #=========================
        self.card_count = model.card_count
        self.write_f06(self.f06_name, end_flag=True)


    def store_beam_oes(self, model, beams, ex, case, Type='strain'):
        return
        if Type == 'stress':
            stress = BeamStressObject(data_code, is_sort1, isubcase, dt=False)
        elif Type == 'strain':
            stress = BeamStrainObject(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(Type)

    def store_beam_oef(self, model, beams, fx, case):
        pass

    def store_rod_oef(self, model, eids, axial, case):
        """
        fills the displacement object
        """
        #print('eids =', eids)
        if len(eids) == 0:
            return
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = False
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OEF',
                    'element_name': 'CONROD', 'format_code':format_code,
                    #'s_code': s_code,
                    'nonlinear_factor': None}
        forces = RealRodForce(data_code, is_sort1, isubcase, dt=False)

        data = []
        i = 0
        #(elementID, axial, torsion) = line
        for (eid, axiali) in zip(eids, axial):
            line = [eid, axiali, 0.]
            data.append(line)
        forces.add_f06_data(data, dt)

        self.rodForces[isubcase] = forces
        #stress.dt = None
        self.iSubcases.append(isubcase)

    def store_rod_oes(self, model, eids, axial, case, Type='stress'):
        """
        fills the displacement object
        """
        #print('eids =', eids)
        if len(eids) == 0:
            return
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = False
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
                    'element_name': 'CONROD', 'format_code':format_code,
                    's_code': s_code,
                    'nonlinear_factor': None}
        if Type == 'stress':
            stress = RodStressObject(data_code, is_sort1, isubcase, dt=False)
        elif Type == 'strain':
            stress = RodStrainObject(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(Type)

        data = []
        i = 0
        #(elementID, axial, torsion, margin_axial, margin_torsion) = line
        for (eid, axiali) in zip(eids, axial):
            line = [eid, axiali, 0., 0., 0.]
            data.append(line)
        stress.add_f06_data(data, dt)

        if Type == 'stress':
            self.rodStress[isubcase] = stress
        elif Type == 'strain':
            self.rodStrain[isubcase] = stress
            stress.dt = None
        else:
            raise NotImplementedError(Type)
        self.iSubcases.append(isubcase)

    def store_displacements(self, model, U, case):
        """
        fills the displacement object
        """
        self.iSubcases = []
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = False
        dt = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OUG',
                    'nonlinear_factor': None}
        disp = DisplacementObject(data_code, is_sort1, isubcase, dt=False)

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
            data.append(line)
        disp.add_f06_data(data, dt)
        self.displacements[isubcase] = disp
        self.iSubcases.append(isubcase)

    def setup_sol_101(self, model, case):
        # the (GridID,componentID) -> internalID
        (self.nidComponentToID, i) = self.build_nid_component_to_id(model)
        (Kgg, Fg, isSPC, isMPC) = self.build_Kgg_Fg(model, case, self.nidComponentToID, i)

        (self.IDtoNidComponents) = reverseDict(self.nidComponentToID)
        print("IDtoNidComponents = ", self.IDtoNidComponents)
        #print("Kgg =\n" + print_annotated_matrix(Kgg, self.IDtoNidComponents))
        #print("Kgg = \n", Kgg)
        #print("iSize = ", i)

        #(Kaa,Fa) = self.Partition(Kgg)
        #sys.exit('verify Kgg')

        print("Kgg = \n%s" % (print_matrix(Kgg)))
        Kaa = partition_dense_symmetric(Kgg, self.iUs)
        print("Kaa = \n%s" % (print_matrix(Kaa)))
        #print("Kaa.shape = ",Kaa.shape)

        #sys.exit('verify Kaa')
        Fa = partition_dense_vector(Fg, self.iUs)
        #print("Kaa = \n%s" % (print_matrix(Kaa)))

        print("Fg = ", Fg)
        print("Fa = ", Fa)
        #print("Us = ", self.Us)

        #self.Us = array(self.Us, 'float64')  # SPC
        #self.Um = array(self.Um, 'float64')  # MPC
        Ua = self.solve_sol_101(Kaa, Fa)
        return Ua, i

    def solve_sol_101(self, Kaa, Fa):
        if 0:
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

        return Ua

    def element_dof_start(self, elem, nids, is3D=True):
        node_ids = elem.nodeIDs()
        index0s = searchsorted(nids, node_ids)
        if is3D:
            index0s *= 6
        else:
            index0s *= 3
        return node_ids, index0s

    def assemble_global_stiffness(self, model, i, Kgg, Dofs, F):
        dof_mapper = []
        print("Kgg.shape", Kgg.shape)

        nnodes = len(model.nodes)
        try:
            nspoints = len(model.spoints)
        except TypeError:
            nspoints = 0
        #ndofs = 6 * nnodes + nspoints

        i = 0
        nids = model.nodes.keys()
        nids.sort()
        print('nids =', nids, i)
        #dof_ids = arange(0, 6*nnodes, 6)

        #dofs_0 = [nid=2, 1] -> searchsorted(nids, nid)[0]

        #for nid in sorted(self.nodes.iterkeys()):
            #nid_dof_mapper[]

        for (eid, elem) in sorted(model.elements.iteritems()):  # CROD, CONROD
            print("----------------------------")
            node_ids, index0s = self.element_dof_start(elem, nids)
            print("node_ids=%s index0s=%s" % (node_ids, index0s))

            # nIJV is the position of the values of K in the dof
            try:
                (K, Ke, dofs, nIJV, Fg, nGrav) = elem.Stiffness(model, node_ids, index0s, self.gravLoad, self.is3D, self.fnorm)
            except TypeError as e:
                msg = 'elem %s must take:\n>>>Stiffness(model, node_ids, index0s, self.gravLoad, self.is3D)' % elem.type
                print(msg)
                raise
            print("K[%s] = \n%s" % (eid, K/self.fnorm))
            print("Ke[%s] = \n%s" % (eid, Ke/self.fnorm))
            #print("dofs = %s" % dofs)

            if 0:
                ndofsi = len(dofs)
                (Ki, Kj) = K.shape
                assert ndofsi==Ki, 'ndofsi=%s Ki=%s' % (ndofsi, Ki)
                assert ndofsi==Kj, 'ndofsi=%s Kj=%s' % (ndofsi, Kj)
                #Kgg[dofs, dofs] = K
                for i, dof1 in enumerate(dofs):
                    for j, dof2 in enumerate(dofs):
                        if abs(K[i, j]) > 0.0:
                            print('i=%s j=%s dof1=%s dof2=%s Ke[i,j]=%s' % (i,j,dof1,dof2, K[i,j]/self.fnorm))
                            Kgg[dof1, dof2] = K[i, j]
                continue

            (Ki, Kj) = K.shape
            ij = 0
            nij2 = []
            print('**Dofs', Dofs.keys())
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

    def apply_SPCs2(self, model, case, nidComponentToID):
        if case.has_parameter('SPC'):
            # get the value, 1 is the options (SPC has no options)
            spcID = case.get_parameter('SPC')[0]
            SpcSet = model.SPC(spcID)

            print(SpcSet)
            for spcSet in SpcSet:
                (typesFound, positionSPCs) = spcSet.organizeConstraints(model)
        return

    def apply_SPCs(self, model, case, nidComponentToID):
        return self.apply_SPCs2(model, case, nidComponentToID)
        isSPC = False
        print('*Us', self.Us)
        print('*iUs', self.iUs)
        if case.has_parameter('SPC'):
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
                    #print("ps = |%s|" % ps)
                    nodes = spc.nodes
                    for nid in nodes:
                        for ips in ps:
                            ips = int(ips)
                            if not(self.is3D) and ips not in ['1', '2', '5']:
                                continue

                            key = (nid, ips)
                            i = nidComponentToID[key]
                            #print("i=%s Us=%s" % (i, 0.0))
                            if i not in self.iUs:
                                self.iUs.append(i)
                                self.Us.append(0.0)

                else:
                    raise NotImplementedError('Invalid SPC...spc=\n%s' % spc)

            print("iUs = ", self.iUs)
            print("Us  = ", self.Us)
            sys.exit('stopping in apply_SPCs')

        print("iUs = ", self.iUs)
        print("Us  = ", self.Us)
        return (isSPC)

    def apply_MPCs(self, model, case, nidComponentToID):
        isMPC = False
        if case.has_parameter('MPC'):
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

            msg = 'MPCs are not supported...stopping in apply_MPCs'
            raise NotImplementedError(msg)
        return (isMPC)

    def assemble_forces(self, model, case, Fg, Dofs):
        """very similar to writeCodeAster loads"""
        #print(model.loads)
        (loadID, junk) = model.caseControlDeck.get_subcase_parameter(case.id, 'LOAD')
        print("loadID = ", loadID)
        LoadSet = model.Load(loadID, 'loadID=%s' % loadID)

        self.gravLoad = array([0., 0., 0.])
        for load_set in LoadSet:
            print("load_set = %r" % str(load_set))
            #print("type", type(load_set))
            (typesFound, forceLoads, momentLoads,
             forceConstraints, momentConstraints,
             gravityLoad) = load_set.organizeLoads(model)

            print('typesFound', typesFound)
            assert isinstance(typesFound, list), type(typesFound)
            assert isinstance(forceLoads, dict), type(forceLoads)
            assert isinstance(momentLoads, dict), type(momentLoads)
            assert isinstance(forceConstraints, dict), type(forceConstraints)
            assert isinstance(momentConstraints, dict), type(momentConstraints)
            assert isinstance(gravityLoad, list), type(gravityLoad)

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

    def write_results(self, case):
        Us = self.Us
        Um = self.Um
        Ua = self.Ua

        iUs = self.iUs
        iUm = self.iUm
        iUa = self.iUa
        pageNum = 1

        if case.has_parameter('DISPLACEMENT'):
            (value, options) = case.get_parameter('DISPLACEMENT')
            if options is not []:
                UgSeparate = [[Ua, iUa], [Us, iUs], [Um, iUm]]
                Ug = departition_dense_vector(UgSeparate)

                result = DisplacementObject(data_code, transient)
                result.add_f06_data()

                if 'PRINT' in options:
                    f06.write(result.write_f06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(self.Title, self.Subtitle))

        if case.has_parameter('SPCFORCES'):
            (value, options) = case.get_parameter('SPCFORCES')
            if options is not []:
                SPCForces = Ksa * Ua + Kss * Us
                if isMPC:
                    SPCForces += Ksm * Um

                result = SPCForcesObject(data_code, transient)
                result.add_f06_data()
                if 'PRINT' in options:
                    f06.write(result.write_f06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(Title, Subtitle))

        if case.has_parameter('MPCFORCES'):
            if options is not []:
                (value, options) = case.get_parameter('MPCFORCES')
                MPCForces = Kma * Ua + Kmm * Um
                if isSPC:
                    MPCForces += Kms * Us

                result = MPCForcesObject(data_code, transient)
                result.add_f06_data()
                if 'PRINT' in options:
                    f06.write(result.write_f06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    f06.write(result.writeOP2(Title, Subtitle))

        if case.has_parameter('GPFORCE'):
            if options is not []:
                (value, options) = case.get_parameter('GPFORCE')
                AppliedLoads = Kaa * Ua
                if isSPC:
                    AppliedLoads += Kas * Us
                if isMPC:
                    AppliedLoads += Kam * Um

                result = AppliedLoadsObject(data_code, transient)
                result.add_f06_data()
                if 'PRINT' in options:
                    f06.write(result.write_f06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.writeOP2(Title, Subtitle))

        if case.has_parameter('STRAIN'):
            if options is not []:
                (value, options) = case.get_parameter('STRAIN')

                for (eid, elem) in sorted(model.elements()):
                    pass

                result = xxxObject(data_code, transient)
                result.add_f06_data()
                if 'PRINT' in options:
                    f06.write(result.write_f06(header, pageStamp, pageNum))
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
                      'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PTUBE',
                      'PBEND',
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
    inputs = sys.argv[1:]
    fnorm = 1.0
    is3D = False
    if len(inputs) == 1:
        bdfName = inputs[0]
    elif len(inputs) == 2:
        bdfName, is3D= inputs
        is3D = bool(is3D)
    elif len(inputs) == 3:
        bdfName, is3D, fnorm = inputs
        is3D = bool(is3D)
        fnorm = float(fnorm)
    s = Solver()
    s.is3D = is3D
    s.fnorm = fnorm
    s.run_solver(bdfName)