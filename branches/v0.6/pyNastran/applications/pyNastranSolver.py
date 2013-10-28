# pylint: disable=E0602,C0103
from __future__ import print_function
import os
import sys
from itertools import izip
from struct import pack

# 3rd party
from numpy import array, zeros, ones
from numpy import searchsorted
from numpy.linalg import solve

from scipy.sparse import coo_matrix

# pyNastran
from pyNastran.utils import list_print
from pyNastran.utils.mathematics import print_matrix, print_annotated_matrix
from pyNastran.bdf.bdf import BDF, SPC, SPC1
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2

# Tables
from pyNastran.op2.tables.oug.oug_displacements import DisplacementObject
#from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject
#from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject

# springs
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import CelasStressObject, CelasStrainObject
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealSpringForce

# rods
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RodStressObject, RodStrainObject, ConrodStressObject, ConrodStrainObject, CtubeStressObject, CtubeStrainObject
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealRodForce, RealConrodForce, RealCtubeForce

# shear
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import ShearStressObject, ShearStrainObject
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealCShearForce

# beams
from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import BeamStressObject, BeamStrainObject
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealCBeamForce


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


def partition_dense_symmetric(A, dofs_in):
    nAll = A.shape[0]
    dofs = getDOF_Set(nAll, dofs_in)
    dofs.sort()
    n = len(dofs)
    A2 = zeros((n, n), 'float64')
    for (i, dofI) in enumerate(dofs):
        for (j, dofJ) in enumerate(dofs):
            v = A[dofI, dofJ]
            if abs(v) >= 1e-8:
                A2[i, j] = v
    return (A2, dofs)


def partition_dense_vector(F, dofs_in):
    nAll = F.shape[0]
    #print("dofs = ", dofs)
    dofs = getDOF_Set(nAll, dofs_in)
    dofs.sort()
    print("dofs = ", dofs)
    n = len(dofs)
    F2 = zeros(n, 'float64')
    for (i, dofI) in enumerate(dofs):
        v = F[dofI]
        if abs(v) >= 1e-8:
            F2[i] = v
    return (F2, dofs)


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


def reverse_dict(A):
    B = {}
    for (key, value) in A.iteritems():
        B[value] = key
    return B


class Solver(F06, OP2):
    """

    Goals:
      - Solves SOL 101
        - calculate Kgg, Fg matrix
        - calculate SPC, MPC, Coord constraints
        - {F} = [K]{x}
        - {x} = [K]^-1 {F}
      - Solve SOL 103 eigenvalue
        - no discussion in limitations and progress
        - calculate Kgg, Mgg, Fg matrix
        - calculate SPC, MPC, Coord constraints
        - solve for eigenvectors
        - {F} = [K]{x} + [M]{xdd}
        -       [K]{x} + [M]{λ^2}{x}
        -       ([K] + [M]{λ^2}){x}
        - let F = 0 and assume {x} != 0
          - [K] + [M]{λ^2} = {0}
          - {λ^2} = -[K][M]^-1

    Progress:
      - Solves a structural problem using matrix partitioning of the SPC set

    TODO:
      - Need to combine solved displacements with original known
        displacements to create displacement set
      - Calculate Stress/Strain for all elements
      - Write the OP2

    Case Control
       - LOAD, SPC, TITLE, LABEL
       - DISP, STRESS, STRAIN
      @todo SETx not supported for output requests; ALL/NONE

    Bulk Data:
      GRID,CORDx
       - Position
       - ps constraint
       - @todo analysis & output coordinate system != 0

      CONROD, CROD/PROD, CTUBE/PTUBE
       - @todo CROD/CTUBE not tested

      CELAS1, CELAS2, CELAS3, CELAS4, PELAS
       - ge not supported (it's damping, so not an issue)
       - @todo test CELAS3/CELAS4

      CSHEAR/PSHEAR
       - @todo what loads can be applied?
       - @todo non-constant a, b not tested
       - @todo f1, f2 not supported and not checked
         @todo support PLOAD2
      MAT1

      LOAD, FORCE, MOMENT
       - coord 0
       - @todo coordx not tested
      PLOAD1
       - @todo validate...
       - distributed load; forces/moments
       - LE: x1=0.0*L; x2=1.0*L
       - FR: x1=0.0; x2=1.0
       - LE_PR:  @todo add this...
       - FR_PR:  @todo add this...
       - @todo support alternate coordinate system
       - @todo static load at x1=0.5

      SPC, SPC1
       - @todo not supported
       - @todo constraints in coordinate system

      MPC/RBE2/RBE3/RBAR
       - @todo not supported

     CQUAD4/CTRIA3
       - @todo not supported
       - @todo PLOAD, PLOAD2, PLOAD4

     CHEXA, CPENTA, CTETRA
       - @todo not supported
       - @todo PLOAD, PLOAD4

    Results:
      CONROD, CROD, PROD
       - FORCE
       - STRESS/STRAIN
         - margins not supported
         - ge not supported (it's damping, so not an issue)

      CELAS1, CELAS2, CELAS3, CELAS4, PELAS
       - STRESS/STRAIN/FORCE
       - @todo test CELAS3/CELAS4

      CBEAM
          - PLOAD1
        - @todo not done... but close
      CSHEAR
        - FORCE/STRESS/STRAIn
          - not calculated
          - no tables created
    """
    def __init__(self):
        F06.__initAlt__(self)
        OP2.__init__(self, '')

        self.pageNum = 1

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

    def _solve(self, K, F, dofs):  # can be overwritten
        r"""solves \f$ [K]{x} = {F}\f$ for \f${x}\f$"""
        print("--------------")
        print("Kaa_norm / %s = \n" % self.fnorm + list_print(K / self.fnorm))
        print("--------------")
        print("Fa = ", F)
        if F[0] == 0.0:
            assert max(F) != min(F), 'no load is applied...'
        print("--------------")

        try:
            U = solve(K, F)
        except:
            failed = []
            faileds = []
            for i, iu in enumerate(dofs):
                absF = abs(F)
                nid, dof = self.IDtoNidComponents[iu]
                #if absF[iu] == 0.0 and ??:
                if K[i, i] == 0.0:
                    failed.append([nid, dof])
                    faileds.append(i)
            msg = self.make_grid_point_singularity_table(failed)
            self.f06_file.write(msg)
            self.f06_file.flush()

            #if 'AUTOSPC' in self.model.params:
            if 1:
                #autospc = self.model.params['AUTOSPC']
                #value = autospc.values[0]
                value = 1

                if value in [1, 'YES']:
                    # figure out what are the DOFs that are removed
                    ilist = list(set(xrange(len(dofs))).difference(set(faileds)))
                    ilist.sort()

                    # remove the DOFs and solve
                    K2 = K[ilist, :][:, ilist]
                    F2 = F[ilist]
                    U2 = solve(K2, F2)

                    # put the removed DOFs back in and set their displacement to 0.0
                    U = zeros(len(F), 'float64')
                    U[ilist] = U2
            else:
                self.f06_file.close()
                raise

        return U

    def run_solver(self, bdfName):
        bdfName = os.path.abspath(bdfName)
        bdf_base, ext = os.path.splitext(bdfName)
        self.f06_name = bdf_base + '.f06'
        self.op2_name = bdf_base + '.op2'
        self.op2_pack_name = bdf_base + '_pack.op2'

        self.model = BDF()
        self.model.cardsToRead = get_cards()
        self.model.readBDF(bdfName)
        cc = self.model.caseControlDeck
        #print cc.subcases
        analysisCases = []
        for (isub, subcase) in sorted(cc.subcases.iteritems()):
            if subcase.has_parameter('LOAD'):
                analysisCases.append(subcase)

        self.f06_file = open(self.f06_name, 'wb')
        self.op2_file = open(self.op2_name, 'wb')
        self.op2_pack_file = open(self.op2_pack_name, 'w')

        self.f06_file.write(self.make_f06_header())
        pageStamp = self.make_stamp(self.Title)
        self.write_summary(self.f06_file, card_count=self.model.card_count)

        #print analysisCases
        for case in analysisCases:
            print(case)
            (value, options) = case.get_parameter('STRESS')
            print("STRESS value   = %s" % value)
            print("STRESS options = %s" % options)

            if case.has_parameter('TEMPERATURE(INITIAL)'):
                (value, options) = case.get_parameter('TEMPERATURE(INITIAL)')
                print('value   = %s' % value)
                print('options = %s' % options)
                raise NotImplementedError('TEMPERATURE(INITIAL) not supported')
                #integrate(B.T*E*alpha*dt*Ads)
            #sys.exit('starting case')
            self.run_case(self.model, case)

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
        else:
            save_results = False
            is_bool = False
        self.case_result_flags[param_name] = [save_results, is_bool]

    def build_nid_component_to_id(self, model):
        i = 0
        nidComponentToID = {}
        for (nid, node) in sorted(model.nodes.iteritems()):  # GRIDs
            for ps in node.ps:
                self.iUs.append(i + int(ps) - 1)
                self.Us.append(0.0)

            nidComponentToID[(nid, 1)] = i
            nidComponentToID[(nid, 2)] = i + 1
            nidComponentToID[(nid, 3)] = i + 2
            nidComponentToID[(nid, 4)] = i + 3
            nidComponentToID[(nid, 5)] = i + 4
            nidComponentToID[(nid, 6)] = i + 5
            i += 6
            #print('iUs[%i] = %s' % (nid, self.iUs))

        if model.spoints:
            for nid in sorted(model.spoints.spoints):  # SPOINTS
                nidComponentToID[(nid, 1)] = i
                i += 1
        return(nidComponentToID, i)

    def run_sol_101(self, model, case):
        #print("case = ", case)
        assert model.sol == 101, 'model.sol=%s is not 101' % (model.sol)

        if "GRDPNT" in model.params and model.params["GRDPNT"] >= 0:
            g0 = model.parameter["GRDPNT"]
            reference_point = None
            if g0 in model.nodes:
                reference_point = model.nodes[g0].Position()
            (mass, cg, I) = model.mass_properties(reference_point, sym_axis=None, num_cpus=1)
            # calculate mass - really should use the Mass matrix...

        ## define IDs of grid point components in matrices
        if 1:
            # analysis
            (Kgg, Fg, n) = self.setup_sol_101(model, case)
            print("------------------------\n")
            print("solving...")
            Ua = self.solve_sol_101(Kgg, Fg)
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
                self._store_displacements(model, U, case)
            q = U
        else:
            n = len(model.nodes)
            q = ones(n, 'float64')

        # =====================================================================
        # results
        #spc_forces = Ksa*Ua + Kss*Us + Ksm*Um
        #mpc_forces = Kma*Ua + Kms*Us + Kmm*Um

        if case.has_parameter('OLOAD'):
            val, options = case.get_parameter('OLOAD')
        del Fg, Kgg

        # =====================================================================
        # springs
        celas1s = []
        celas2s = []
        celas3s = []
        celas4s = []

        # bars
        crods = []
        conrods = []
        ctubes = []

        # bars
        cbeams = []
        cbars = []

        # shears
        cshears = []

        # shells
        cquad4s = []
        ctria3s = []

        # solids
        ctetra4s = []
        cpenta5s = []
        chexa8s = []

        type_map = {
            # springs
            'CELAS1' : celas1s,
            'CELAS2' : celas2s,
            'CELAS3' : celas3s,
            'CELAS4' : celas4s,

            # rods
            'CONROD' : conrods,
            'CROD'   : crods,
            'CTUBE'  : ctubes,

            # bars
            'CBEAM'  : cbeams,
            'CBAR'   : cbars,

            # shear
            'CSHEAR'  : cshears,

            # shells
            'CQUAD4'  : cquad4s,
            'CTRIA3'  : ctria3s,

            'CTETRA'  : ctetra4s,
            'CPENTA'  : cpenta5s,
            'CHEXA'  : chexa8s,
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

        # springs
        ncelas1s = len(celas1s)  # doesn't support cid != 0
        ncelas2s = len(celas2s)  # doesn't support cid != 0
        ncelas3s = len(celas3s)  # not tested
        ncelas4s = len(celas4s)  # not tested

        # rods
        ncrods = len(crods)    # not tested
        nconrods = len(conrods)
        nctubes = len(ctubes)  # not tested

        # bars
        ncbars = len(cbars)   # not tested
        ncbeams = len(cbeams) # half implemented

        # shear
        ncshears = len(cshears) # half implemented

        # not implemented - shells
        ncquad4s = len(cquad4s)
        nctria3s = len(ctria3s)

        # not implemented - shells
        nctetra4s = len(ctetra4s)
        ncpenta5s = len(cpenta5s)
        nchexa8s = len(chexa8s)

        #=========================
        # rods
        crods = array(crods)
        conrods = array(conrods)
        ctubes = array(ctubes)

        # shears
        cshears = array(cshears)

        # shells
        ctria3s = array(ctria3s)
        cquad4s = array(cquad4s)

        #=========================
        if self.is_stress or self.is_strain or self.is_force:
            # SPRINGS
            elementTypes = ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']
            for elementType in elementTypes:
                if elementType == 'CELAS1':
                    n = ncelas1s
                    eids = celas1s
                elif elementType == 'CELAS2':
                    n = ncelas2s
                    eids = celas2s
                elif elementType == 'CELAS3':
                    n = ncelas3s
                    eids = celas3s
                elif elementType == 'CELAS4':
                    n = ncelas4s
                    eids = celas4s
                else:
                    raise NotImplementedError(elementType)
                if n:
                    o1 = zeros(n, 'float64')
                    e1 = zeros(n, 'float64')
                    f1 = zeros(n, 'float64')

                    for i, eid in enumerate(eids):
                        element = elements[eid]
                        (e1i, o1i, f1i) = element.displacement_stress(model, q, self.nidComponentToID)
                        o1[i] = o1i
                        e1[i] = e1i
                        f1[i] = f1i
                    #if self.is_strain:
                    self._store_spring_oes(model, eids, e1, case, elementType, Type='strain')
                    #if self.is_stress:
                    self._store_spring_oes(model, eids, o1, case, elementType, Type='stress')
                    #if self.is_force:
                    self._store_spring_oef(model, eids, f1, case, elementType)
                    del e1
                    del o1
                    del f1

            # RODS
            elementTypes = ['CROD', 'CONROD', 'CTUBE']
            for elementType in elementTypes:
                if elementType == 'CROD':
                    n = ncrods
                    eids = crods
                elif elementType == 'CONROD':
                    n = nconrods
                    eids = conrods
                elif elementType == 'CTUBE':
                    n = nctubes
                    eids = ctubes
                else:
                    raise NotImplementedError(elementType)

                if n:
                    o1 = zeros(n, 'float64')
                    e1 = zeros(n, 'float64')
                    f1 = zeros(n, 'float64')

                    o4 = zeros(n, 'float64')
                    e4 = zeros(n, 'float64')
                    f4 = zeros(n, 'float64')

                    #margin_1 =
                    #margin_2 =
                    #margin_12 =
                    for i, eid in enumerate(eids):
                        element = elements[eid]
                        (e1i, e4i,
                         o1i, o4i,
                         f1i, f4i) = element.displacement_stress(model, q, self.nidComponentToID)
                        o1[i] = o1i
                        e1[i] = e1i
                        f1[i] = f1i

                        o4[i] = o4i
                        e4[i] = e4i
                        f4[i] = f4i
                    if self.is_strain:
                        self._store_rod_oes(model, eids, e1, e4, case, elementType, Type='strain')
                    del e1, e4
                    if self.is_stress:
                        self._store_rod_oes(model, eids, o1, o4, case, elementType, Type='stress')
                    del o1, o4
                    if self.is_force:
                        self._store_rod_oef(model, eids, f1, f4, case, elementType)
                    del f1, f4

            #=========================
            # CHSEAR
            if ncshears:
                print("ncshears", ncshears)
                print("cshears", cshears)
                stress = zeros((ncshears, 3), 'float64')
                strain = zeros((ncshears, 3), 'float64')
                force  = zeros((ncshears, 1), 'float64')
                for i, eid in enumerate(cshears):
                    element = elements[eid]
                    (stressi, straini, forcei) = element.displacement_stress(model, q, self.nidComponentToID)
                    stress[i, :] = stressi
                    strain[i, :] = straini
                    #force[i, :] = forcei
                #if self.is_strain:
                self._store_cshear_oes(model, cshears, strain, case, 'CSHEAR', Type='strain')
                #if self.is_stress:
                self._store_cshear_oes(model, cshears, stress, case, 'CSHEAR', Type='stress')
                #if self.is_force:
                #self._store_cshear_oef(model, cshears, force, case)
                del stress
                del strain
                del force

            #=========================
            # BARS / BEAMS
            if ncbeams:
                print("ncbeams", ncbeams)
                print("cbeams", cbeams)
                o1 = zeros(ncbeams, 'float64')
                e1 = zeros(ncbeams, 'float64')
                f1 = zeros(ncbeams, 'float64')
                for i, eid in enumerate(cbeams):
                    element = elements[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID)
                    o1[i] = oxi
                    e1[i] = exi
                    f1[i] = fxi
                #if self.is_strain:
                self._store_beam_oes(model, cbeams, e1, case, Type='strain')
                #if self.is_stress:
                self._store_beam_oes(model, cbeams, o1, case, Type='stress')
                #if self.is_force:
                self._store_beam_oef(model, cbeams, f1, case)
                del e1
                del o1
                del f1

            #=========================
            # SHELLS
            if nctria3s or ncquad4s:
                print("nctria3", nctria3s)
                stress = zeros((nctria3s+ncquad4s, 3), 'float64')
                strain = zeros((nctria3s+ncquad4s, 3), 'float64')
                force  = zeros((nctria3s+ncquad4s, 3), 'float64')

            i0 = i
            if nctria3s:
                for i, eid in enumerate(ctria3s):
                    element = elements[eid]
                    (stressi, straini, forcei) = element.displacement_stress(model, q, self.nidComponentToID)
                    stress[i, :] = stressi
                    strain[i, :] = straini
                    force[i, :] = forcei
                i0 = i

            if ncquad4s:
                for i, eid in enumerate(ncquad4s):
                    element = elements[eid]
                    (stressi, straini, forcei) = element.displacement_stress(model, q, self.nidComponentToID)
                    stress[i0+i, :] = stressi
                    strain[i0+i, :] = straini
                    force[i0+i, :] = forcei

            if nctria3s or ncquad4s:
                #if self.is_strain:
                self._store_plate_oes(model, cbeams, stress, case, Type='strain')
                #if self.is_stress:
                self._store_plate_oes(model, cbeams, strain, case, Type='stress')
                #if self.is_force:
                self._store_plate_oef(model, cbeams, force, case)
                del stress, strain, force

            # SOLIDS
        #=========================
        self.write_f06(self.f06_file, end_flag=True)
        self.write_op2(self.op2_file, packing=True)
        self.write_op2(self.op2_pack_file, packing=False)

    def _op2_header(self, f, packing=True):
        data = [4,3,4,
         1,28,12,
         4,7,4,   # 7
         'NASTRAN FORT TAPE ID CODE - ',   # 28 = 7*4
         4,2,4,   # 7
         4,-1,4,
         4,0,4,

         4,2,4,
         4,0,4,
         4,2,4]
        if packing:
            f.write(pack('9i28s18i', *data))
        if not packing:
            f.write(str(data)+'\n')

    def write_op2(self, f, packing=False):
        results = [self.displacements]
        header = None
        pageStamp = None

        self._op2_header(f, packing=packing)

        for result in results:
            for subcase, case in sorted(result.iteritems()):
                case.write_op2(header, pageStamp, f, is_mag_phase=False, packing=packing)
                if not packing:
                    f.write('\n')
        marker1 = [4,  0, 4]
        marker2 = [4,  0, 4]
        marker3 = [4,  0, 4]
        marker = marker1 + marker2 + marker3
        if packing:
            nmarker = len(marker)
            p = pack('%ii' % nmarker, *marker)
            f.write(p)
        else:
            f.write(str(marker)+'\n')
        f.close()

    def _store_beam_oes(self, model, eids, axial, case, elementType='CBEAM', Type='strain'):
        #print('eids =', eids)
        if len(eids) == 0:
            return
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = False
        dt = None
        format_code = 1  # ???
        s_code = 1

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
                    'element_name': elementType, 'format_code':format_code,
                    's_code': s_code,
                    'nonlinear_factor': None}
        if Type == 'stress':
            if elementType == 'CBEAM':
                stress = BeamStressObject(data_code, is_sort1, isubcase, dt=False)
        elif Type == 'strain':
            if elementType == 'CBEAM':
                stress = BeamStrainObject(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(Type)

        data = []
        i = 0

        for (eid, axiali) in zip(eids, axial):
            element = model.Element(eid)
            n1, n2 = element.nodeIDs()
            print(n1, n2)
            #      (eid, grid, sd,  sxc,   sxd, sxe, sxf,  smax, smin, mst, msc) = out
            line = [eid, n1,   0.0, axiali, 0., 0.0,  0.0, 0.0, 0.0,  0.0,  0.0]
            data.append(line)

            line = [eid, n2,   1.0, axiali, 0., 0.0,  0.0, 0.0, 0.0,  0.0,  0.0]
            data.append(line)
        stress.add_f06_data(data, dt)

        if elementType == 'CBEAM' and Type == 'stress':
            self.beamStress[isubcase] = stress
        elif elementType == 'CBEAM' and Type == 'strain':
            self.beamStrain[isubcase] = stress
        else:
            raise NotImplementedError('elementType=%r Type=%r' % (elementType, Type))
        stress.dt = None


    def _store_beam_oef(self, model, eids, fx, case, elementType='CBEAM'):
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
                    'element_name': elementType, 'format_code':format_code,
                    #'s_code': s_code,
                    'nonlinear_factor': None}

        if elementType == 'CBEAM':
            forces = RealCBeamForce(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(elementType)

        data = []
        i = 0
        for (eid, fxi) in zip(eids, fx):
            element = model.Element(eid)
            n1, n2 = element.nodeIDs()
            print('***(*', n1, n2)
            #      [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
            line = [eid, n1, 0.0, 0.,   0.,   0., 0., 0.,   0.,  0.0]
            data.append(line)
            line = [eid, n1, 1.0, 0.,   0.,   0., 0., 0.,   0.,  0.0]
            data.append(line)
            line = [eid, n2, 0.0, 0.,   0.,   0., 0., 0.,   0.,  0.0]
            data.append(line)
            line = [eid, n2, 1.0, 0.,   0.,   0., 0., 0.,   0.,  0.0]
            #data.append(line)
        print(data)
        forces.add_f06_data(data, dt)

        if elementType == 'CBEAM':
            self.beamForces[isubcase] = forces
        else:
            raise NotImplementedError(elementType)
        #stress.dt = None

    def _OEF_f06_header(self, case, elementType):
        analysis_code = 1
        transient = False
        is_sort1 = False
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OEF',
                    'element_name': elementType, 'format_code':format_code,
                    #'s_code': s_code,
                    'nonlinear_factor': None}
        return data_code

    def _store_cshear_oes(self, model, eids, results, case, elementType, Type='strain'):
        """
        self.shearForces = {}
        self.shearStress = {}
        self.shearStrain = {}

                                     S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )
      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
          328        1.721350E+03   1.570314E+03   7.2E+01
        """
        if len(eids) == 0:
            return
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = True
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
                    'element_name': elementType, 'format_code':format_code,
                    's_code': s_code,
                    'nonlinear_factor': None}

        if Type == 'stress':
            #if elementType == 'CELAS2':
            stress = ShearStressObject(data_code, is_sort1, isubcase, dt=None)
            #else:
                #raise NotImplementedError(elementType)
        elif Type == 'strain':
            #if elementType == 'CELAS2':
            stress = ShearStrainObject(data_code, is_sort1, isubcase, dt=None)
            #else:
                #raise NotImplementedError(elementType)
        else:
            raise NotImplementedError(Type)

        data = []
        i = 0
        #(elementID, max_shear, avg_shear, margin) = line
        for i, eid in enumerate(eids):
            resultsi = results[i, :]
            line = [eid] + list(resultsi)
            data.append(line)
        stress.add_f06_data(data, dt)

        if Type == 'stress':
            self.shearStress[isubcase] = stress
        elif Type == 'strain':
            self.shearStrain[isubcase] = stress
        else:
            raise NotImplementedError('elementType=%r Type=%r' % (elementType, Type))
        #stress.dt = None

    def _store_spring_oes(self, model, eids, axial, case, elementType, Type='stress'):
        if len(eids) == 0:
            return
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = True
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
                    'element_name': elementType, 'format_code':format_code,
                    's_code': s_code,
                    'nonlinear_factor': None}

        if Type == 'stress':
            #if elementType == 'CELAS2':
            stress = CelasStressObject(data_code, is_sort1, isubcase, dt=None)
            #else:
                #raise NotImplementedError(elementType)
        elif Type == 'strain':
            #if elementType == 'CELAS2':
            stress = CelasStrainObject(data_code, is_sort1, isubcase, dt=None)
            #else:
                #raise NotImplementedError(elementType)
        else:
            raise NotImplementedError(Type)

        data = []
        i = 0
        #(elementID, stress) = line
        for (eid, axiali) in zip(eids, axial):
            line = [eid, axiali]
            data.append(line)
        stress.add_f06_data(data, dt)

        if Type == 'stress':
            self.celasStress[isubcase] = stress
        elif Type == 'strain':
            self.celasStrain[isubcase] = stress
        else:
            raise NotImplementedError('elementType=%r Type=%r' % (elementType, Type))
        #stress.dt = None

    def _store_spring_oef(self, model, eids, axial, case, elementType):
        if len(eids) == 0:
            return
        data_code = self._OEF_f06_header(case, elementType)

        is_sort1 = True
        isubcase = case.id
        dt = None
        #if elementType == 'CELAS2':
        forces = RealSpringForce(data_code, is_sort1, isubcase, dt=None)
        #else:
            #raise NotImplementedError(elementType)

        data = []
        i = 0
        #(elementID, axial) = line
        for (eid, axiali) in zip(eids, axial):
            line = [eid, axiali]
            data.append(line)
        forces.add_f06_data(data, dt)

        #if elementType == 'CELAS2':
        self.springForces[isubcase] = forces
        #else:
            #raise NotImplementedError(elementType)
        #stress.dt = None

    def _store_rod_oef(self, model, eids, axial, torsion, case, elementType):
        if len(eids) == 0:
            return
        data_code = self._OEF_f06_header(case, elementType)
        is_sort1 = True
        isubcase = case.id

        if elementType == 'CROD':
            forces = RealRodForce(data_code, is_sort1, isubcase, dt=False)
        elif elementType == 'CONROD':
            forces = RealConrodForce(data_code, is_sort1, isubcase, dt=False)
        elif elementType == 'CTUBE':
            forces = RealCtubeForce(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(elementType)

        data = []
        i = 0
        #(elementID, axial, torsion) = line
        for (eid, axiali, torsioni) in zip(eids, axial, torsion):
            line = [eid, axiali, torsioni]
            data.append(line)

        dt = None
        forces.add_f06_data(data, dt)

        if elementType == 'CROD':
            self.rodForces[isubcase] = forces
        elif elementType == 'CONROD':
            self.conrodForces[isubcase] = forces
        elif elementType == 'CTUBE':
            self.ctubeForces[isubcase] = forces
        else:
            raise NotImplementedError(elementType)
        #stress.dt = None

    def _store_rod_oes(self, model, eids, axial, torsion, case, elementType, Type='stress'):
        if len(eids) == 0:
            return
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = True
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
                    'element_name': elementType, 'format_code':format_code,
                    's_code': s_code,
                    'nonlinear_factor': None}
        if Type == 'stress':
            if elementType == 'CROD':
                stress = RodStressObject(data_code, is_sort1, isubcase, dt=False)
            elif elementType == 'CONROD':
                stress = ConrodStressObject(data_code, is_sort1, isubcase, dt=False)
            elif elementType == 'CTUBE':
                stress = CtubeStressObject(data_code, is_sort1, isubcase, dt=False)
        elif Type == 'strain':
            if elementType == 'CROD':
                stress = RodStrainObject(data_code, is_sort1, isubcase, dt=False)
            elif elementType == 'CONROD':
                stress = ConrodStrainObject(data_code, is_sort1, isubcase, dt=False)
            elif elementType == 'CTUBE':
                stress = CtubeStrainObject(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(Type)

        data = []
        i = 0
        #(elementID, axial, torsion, margin_axial, margin_torsion) = line
        for (eid, axiali, torsioni) in zip(eids, axial, torsion):
            line = [eid, axiali, 0., torsioni, 0.]
            data.append(line)
        stress.add_f06_data(data, dt)

        if elementType == 'CROD' and Type == 'stress':
            self.rodStress[isubcase] = stress
        elif elementType == 'CROD' and Type == 'strain':
            self.rodStrain[isubcase] = stress
        elif elementType == 'CONROD' and Type == 'stress':
            self.conrodStress[isubcase] = stress
        elif elementType == 'CONROD' and Type == 'strain':
            self.conrodStrain[isubcase] = stress
        elif elementType == 'CTUBE' and Type == 'stress':
            self.ctubeStress[isubcase] = stress
        elif elementType == 'CTUBE' and Type == 'strain':
            self.ctubeStrain[isubcase] = stress
        else:
            raise NotImplementedError('elementType=%r Type=%r' % (elementType, Type))
        stress.dt = None

    def _store_displacements(self, model, U, case):
        """
        fills the displacement object
        """
        self.iSubcases = []
        analysis_code = 1
        transient = False
        isubcase = case.id
        is_sort1 = True
        dt = None
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OUG',
                    'nonlinear_factor': None}
        disp = DisplacementObject(data_code, is_sort1, isubcase, dt=None)

        data = []

        i = 0
        #(nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
        for (nid, node) in sorted(model.nodes.iteritems()):
            if node.type == 'GRID':
                line = [nid, 'G']
                xyz = U[i:i + 6]  # 1,2,3,4,5,6
                line += xyz
                i += 6
            #elif node.type == 'SPOINT':
                #line = [nid, 'S', U[i], 0., 0., 0., 0., 0.]
                #i += 1
            else:
                raise NotImplementedError('node.type=%s' % node.type)
            data.append(line)
        disp.add_f06_data(data, dt)
        self.displacements[isubcase] = disp
        self.iSubcases.append(isubcase)

    def setup_sol_101(self, model, case):
        # the (GridID,componentID) -> internalID
        (self.nidComponentToID, i) = self.build_nid_component_to_id(model)

        #=====================
        isSPC = self.apply_SPCs(model, case, self.nidComponentToID)
        isMPC = self.apply_MPCs(model, case, self.nidComponentToID)

        #spcDOFs = self.iUs
        #mpcDOFs = self.iUm

        #Mgg = zeros((i, i), 'float64')

        Kgg, Kgg_sparse = self.assemble_global_stiffness(model, i, self.nidComponentToID)
        Fg = self.assemble_forces(model, i, case, self.nidComponentToID)
        return Kgg, Fg, i

    def solve_sol_101(self, Kgg, Fg):
        #=====================

        (self.IDtoNidComponents) = reverse_dict(self.nidComponentToID)
        print("IDtoNidComponents = ", self.IDtoNidComponents)
        #print("Kgg =\n" + print_annotated_matrix(Kgg, self.IDtoNidComponents))
        #print("Kgg = \n", Kgg)
        #print("iSize = ", i)

        #(Kaa,Fa) = self.Partition(Kgg)
        #sys.exit('verify Kgg')

        #print("Kgg = \n%s" % (print_matrix(Kgg / self.fnorm)))
        Kaa, dofs2 = partition_dense_symmetric(Kgg, self.iUs)
        print("Kaa = \n%s" % (print_matrix(Kaa / self.fnorm)))
        #print("Kaa.shape = ",Kaa.shape)

        #sys.exit('verify Kaa')
        Fa, _dofs2 = partition_dense_vector(Fg, self.iUs)
        #print("Kaa = \n%s" % (print_matrix(Kaa)))

        print("Fg = ", Fg)
        print("Fa = ", Fa)
        #print("Us = ", self.Us)

        #self.Us = array(self.Us, 'float64')  # SPC
        #self.Um = array(self.Um, 'float64')  # MPC

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

        Ua = self._solve(Kaa, Fa, dofs2)
        #self.Um = Kma*Ua

        return Ua

    def element_dof_start(self, elem, nids):
        node_ids = elem.nodeIDs()
        index0s = searchsorted(nids, node_ids)
        index0s *= 6
        return node_ids, index0s

    def assemble_global_stiffness(self, model, i, Dofs):
        Kgg = zeros((i, i), 'float64')

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
                (K, dofs, nIJV) = elem.Stiffness(model, node_ids, index0s, self.fnorm)
            except TypeError as e:
                msg = 'elem %s must take:\n>>>Stiffness(model, node_ids, index0s, self.fnorm)' % elem.type
                print(msg)
                raise
            #(Fg, nGrav) = elem.Fg(model, self.gravLoad, fnorm)
            #for (fg, dof) in izip(Fg, nGrav):
                ##print("dof = ",dof)
                #if dof in Dofs:
                    #F[Dofs[dof]] += fg

            print("K[%s] = \n%s" % (eid, list_print(K / self.fnorm)))
            #print("Ke[%s] = \n%s" % (eid, Ke / self.fnorm))
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
                            Kgg[dof1, dof2] += K[i, j]
                continue

            (Ki, Kj) = K.shape
            ij = 0
            nij2 = []
            print('**Dofs', Dofs.keys())
            for ij in nIJV:
                nij2.append(Dofs[ij])
            print('nij2', nij2)

            for i in xrange(Ki):
                for j in xrange(Kj):
                    #kij = K[i, j]
                    #if abs(kij)>1e-8:
                    #print('niJV = ',nIJV[ij])
                    ii = nij2[i]
                    jj = nij2[j]
                    #dof = Dofs[n]
                    if K[i, j] / fnorm > 1000.0:
                        asdf
                    Kgg[ii, jj] += K[i, j]

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

        #Kgg_sparse = coo_matrix((entries, (rows, cols)), shape=(i, i))
        Kgg_sparse = None
        return Kgg, Kgg_sparse

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

    def assemble_forces(self, model, i, case, Dofs):
        """very similar to writeCodeAster loads"""
        Fg = zeros(i, 'float64')
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
            if not (isinstance(typesFound, list) or  isinstance(typesFound, set)):
                raise RuntimeError(type(typesFound))
            assert isinstance(forceLoads, dict), type(forceLoads)
            assert isinstance(momentLoads, dict), type(momentLoads)
            assert isinstance(forceConstraints, dict), type(forceConstraints)
            assert isinstance(momentConstraints, dict), type(momentConstraints)
            assert isinstance(gravityLoad, list), type(gravityLoad)

            nids = set([])
            for nid in forceLoads:
                nids.add(nid)
            for nid in momentLoads:
                nids.add(nid)

            if gravityLoad != []:
                print("gravityLoad = ", gravityLoad)
                self.gravLoad += gravityLoad

            for nid in nids:
                print("nid = ", nid)

                if nid in forceLoads:
                    force = forceLoads[nid]
                    #print("force", force)
                    if abs(force[0]) > 0.:
                        Fg[Dofs[(nid, 1)]] += force[0]
                    if abs(force[1]) > 0.:
                        Fg[Dofs[(nid, 2)]] += force[1]
                    if abs(force[2]) > 0.:
                        Fg[Dofs[(nid, 3)]] += force[2]

                if nid in momentLoads:
                    moment = momentLoads[nid]
                    #print("moment", moment)
                    if abs(moment[0]) > 0.:
                        Fg[Dofs[(nid, 4)]] += moment[0]
                    if abs(moment[1]) > 0.:
                        Fg[Dofs[(nid, 5)]] += moment[1]
                    if abs(moment[2]) > 0.:
                        Fg[Dofs[(nid, 6)]] += moment[2]

        if sum(abs(self.gravLoad)) > 0.0:
            for (eid, elem) in sorted(model.elements.iteritems()):  # CROD, CONROD
                print("----------------------------")
                node_ids, index0s = self.element_dof_start(elem, nids)
                print("node_ids=%s index0s=%s" % (node_ids, index0s))

                # nIJV is the position of the values of K in the dof
                (Fgi, nGrav) = elem.Fg(model, self.gravLoad, fnorm)
                for (fg, dof) in izip(Fgi, nGrav):
                    #print("dof = ",dof)
                    if dof in Dofs:
                        Fg[Dofs[dof]] += fg
        print("Fg  = ", Fg)
        return Fg

    def write_results(self, case):
        Us = self.Us  # SPC set
        Um = self.Um  # MPC set
        Ua = self.Ua  # constrained set
        #Ug - global set

        iUs = self.iUs
        iUm = self.iUm
        iUa = self.iUa
        pageNum = 1

        if case.has_parameter('DISPLACEMENT'):
            (value, options) = case.get_parameter('DISPLACEMENT')
            if options is not []:
                Ug_separate = [[Ua, iUa], [Us, iUs], [Um, iUm]]
                Ug = departition_dense_vector(Ug_separate)

                result = DisplacementObject(data_code, transient)
                result.add_f06_data()

                if 'PRINT' in options:
                    f06.write(result.write_f06(header, pageStamp, pageNum))
                if 'PLOT' in options:
                    op2.write(result.write_op2(self.Title, self.Subtitle))

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
                    op2.write(result.write_op2(Title, Subtitle))

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
                    f06.write(result.write_op2(Title, Subtitle))

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
                    op2.write(result.write_op2(Title, Subtitle))

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
                    op2.write(result.write_op2(Title, Subtitle))


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
    if 1:
        inputs = sys.argv[1:]
        fnorm = 1.0
        if len(inputs) == 1:
            bdfName = inputs[0]
        elif len(inputs) == 2:
            bdfName, fnorm = inputs
            fnorm = float(fnorm)
        s = Solver()
        s.fnorm = fnorm
        s.run_solver(bdfName)

        #if os.path.exists(s.op2_name):
            #op2_name = s.op2_name
            #op2 = OP2(op2FileName=op2_name, make_geom=False, debug=True, log=None)
    else:
        op2_name = 'solid_bending.op2'
        op2 = OP2(op2_name)
        op2.make_op2_debug = True
        op2.read_op2()
