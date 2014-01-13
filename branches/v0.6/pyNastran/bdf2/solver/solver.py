# -*- coding: utf-8 -*-
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
from pyNastran.f06.f06Writer import sorted_bulk_data_header
from pyNastran.utils import list_print
from pyNastran.utils.mathematics import print_matrix, print_annotated_matrix
#from pyNastran.bdf.bdf import BDF
from pyNastran.bdf2.bdf import BDF, SPC, SPC1
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

def remove_dofs(dofsAll, dofs_remove):
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
        -       [K]{x} + [M]{s^2}{x}
        -       ([K] + [M]{s^2}){x}
        - let F = 0 and assume {x} != 0
          - [K] + [M]{s^2} = {0}
          - {s^2} = -[K][M]^-1

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
       - @todo CTUBE not tested

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
       - @todo constraints in alternate coordinate system (specified by GRID cards)

      MPC/RBE2/RBE3/RBAR
       - @todo not supported

     CQUAD4/CTRIA3
       - @todo not supported
       - @todo PLOAD, PLOAD2, PLOAD4

     CHEXA, CPENTA, CTETRA
       - @todo not supported
       - @todo PLOAD, PLOAD4

    Results:
      CONROD, CROD/PROD
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
        - FORCE/STRESS/STRAIN
          - not calculated
          - no tables created
    """
    def __init__(self, fargs):
        F06.__init_data__(self)
        OP2.__init__(self, '')

        self.pageNum = 1
        self.fargs = fargs

        # normalization of stiffness matrix
        self.knorm = fargs['--k']
        # normalization of load vector
        self.fnorm = fargs['--f']
        # normalization of mass matrix
        self.mnorm = fargs['--m']

        self.iSubcases = []
        self.nU = 0
        self.nUs = 0
        self.nUm = 0

        #==============================
        #: displacements
        self.U = None


        #==============================
        #: Degrees-of-freedom eliminated by single-point constraints that are
        #: included in boundary condition changes and by the AUTOSPC feature.
        self.iUsb = []
        self.Usb = []
        #==============================
        #: Degrees-of-freedom eliminated by single-point constraints that are
        #: specified on the PS field on GRID Bulk Data entries.
        self.iUsg = []
        self.Usg = []
        #==============================
        #: s = sb + sg
        #: all degrees-of-freedom eliminated by single point constraints
        self.Us = None
        self.iUs = None
        #==============================

        #==============================


        #==============================
        #: Degrees-of-freedom eliminated by multipoint constraints.
        self.Ump = []
        self.iUmp = []
        self.jUmp = []

        #: Degrees-of-freedom eliminated by multipoint constraints created by the
        #: rigid elements using the LGELIM method on the Case Control command
        #: RIGID.
        self.Umr = []
        self.iUmr = []
        self.jUmr = []

        # m = mp + mr
        # all degrees-of-freedom eliminated by multiple constraints
        self.iUm = None
        self.jUm = None
        self.Um = None

        #==============================

        self.Ub = []
        self.iUb = []
        self.Uc = []
        self.iUc = []
        self.Ue = []
        self.iUe = []
        self.Uj = []
        self.iUj = []
        self.Uk = []
        self.iUk = []
        self.Ulm = []
        self.iUlm = []

        self.Uq = []
        self.iUq = []
        self.Ur = []
        self.iUr = []
        self.Uo = []
        self.iUo = []
        self.Usa = []
        self.iUsa = []
        #==============================
        self.case_result_flags = {}

    def _solve(self, K, F, dofs):  # can be overwritten
        r"""solves \f$ [K]{x} = {F}\f$ for \f${x}\f$"""
        print("--------------")
        print("Kaa_norm / %s = \n" % self.knorm + list_print(K / self.knorm))
        print("--------------")
        print("Fa/%g = %s" % (self.fnorm, F / self.fnorm))
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

    def run_solver(self):
        fargs = self.fargs
        bdf_filename = os.path.abspath(fargs['BDFNAME'])
        bdf_base = os.path.abspath(fargs['BDFBASE'])

        #bdf_base, ext = os.path.splitext(bdfName)
        self.f06_name = bdf_base + '.f06'
        self.op2_name = bdf_base + '.op2'
        self.op2_pack_name = bdf_base + '_pack.op2'

        self.f06_file = open(self.f06_name, 'wb')
        #self.op2_file = open(self.op2_name, 'wb')
        #self.op2_pack_file = open(self.op2_pack_name, 'w')
        self.op2_file = None
        self.op2_pack_file = None


        self.f06_file.write(self.make_f06_header())
        #self.f06_file.write(sorted_bulk_data_header())
        pageStamp = self.make_stamp(self.Title)

        #------------------------------------------
        # start of analysis

        self.model = BDF()
        self.model.cards_to_read = get_cards()
        self.model.f06 = self.f06_file

        if 1:
            data = {
                'bar1_a': 1.0,
                'bar2_a': 1.0,
                'bar3_a': 1.0,
                'youngs': 5e6,
                'loadmag': 1000.0,
                'loadx': 0.5,
                'loady': 1.0,
                'rho' : 2.0,
            }
            self.model.set_dynamic_syntax(data)
        self.model.read_bdf(bdf_filename)
        cc = self.model.caseControlDeck
        #print cc.subcases
        analysisCases = []
        for (isub, subcase) in sorted(cc.subcases.iteritems()):
            if subcase.has_parameter('LOAD'):
                analysisCases.append(subcase)

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
        self.f06_file.close()
        if self.op2_file is not None:
            self.op2_file.close()
            self.op2_pack_file.close()

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

            if case.has_parameter('LABEL'):
                (self.label, options) = case.get_parameter('LABEL')
            else:
                self.label = ''
            self.iSubcaseNameMap[isubcase] = [self.Subtitle, self.label]

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
        if model.grid.n:
            cd = set(model.grid.cd)
            if cd == 0:
                raise NotImplementedError('Set cd=0; cd=%s is not supported.')

        for ni in xrange(model.grid.n):  # GRIDs
            nid = model.grid.node_id[ni]
            ps = model.grid.ps[ni]

            while ps > 0:
                psi = ps % 10
                self.iUsg.append(i + psi - 1)
                self.Usg.append(0.0)
                ps = ps // 10

            nidComponentToID[(nid, 1)] = i
            nidComponentToID[(nid, 2)] = i + 1
            nidComponentToID[(nid, 3)] = i + 2
            nidComponentToID[(nid, 4)] = i + 3
            nidComponentToID[(nid, 5)] = i + 4
            nidComponentToID[(nid, 6)] = i + 5
            i += 6
            #print('iUsg[%i] = %s' % (nid, self.iUsg))

        spoint = model.spoint
        if spoint.n:
            for nid in sorted(model.spoint.spoint):  # SPOINTS
                nidComponentToID[(nid, 1)] = i
                i += 1
        assert i > 0, 'no DOFs'

        #: starting index for MPC cards
        self.mp_index = model.grid.n + spoint.n

        return(nidComponentToID, i)

    def run_sol_101(self, model, case):
        #print("case = ", case)
        assert model.sol == 101, 'model.sol=%s is not 101' % (model.sol)

        if 'WTMASS' in self.params:
            wtmass = self.params['WTMASS'].value1
        else:
            wtmass = 1.0

        if 'COUPMASS' in self.params:
            coupmass = self.params['COUPMASS'].value1
        else:
            coupmass = -1

        if case.has_parameter('FMETHOD'):
            iflutter = self.subcase.get_parameter('FMETHOD')
            self.flutter[iflutter]
            self.flfact[iflutter]

        if case.has_parameter('METHOD'):
            imethod = self.subcase.get_parameter('METHOD')
            self.eigb[imethod]
            self.eigc[imethod]
            self.eigr[imethod]
            self.eigrl[imethod]

        if "GRDPNT" in model.params and model.params["GRDPNT"] >= 0:
            g0 = model.params["GRDPNT"]
            reference_point = None
            if g0 in model.nodes:
                reference_point = model.nodes[g0].Position()
            #(mass, cg, I) = model.mass_properties(reference_point, sym_axis=None, num_cpus=1)
            #mass *= wtmass
            #I *= wtmass
            # calculate mass - really should use the Mass matrix...

        ## define IDs of grid point components in matrices
        if 1:
            # analysis
            (Kgg, Fg, n) = self.setup_sol_101(model, case)
            self.build_dof_sets()
            print("------------------------\n")
            print("solving...")
            Ua = self.solve_sol_101(Kgg, Fg)
            print("Ua = ", Ua)
            print("Us = ", self.Us)

            dofsAll = set([i for i in xrange(n)])
            #dofsA = remove_dofs(remove_dofs(dofsAll, self.iUs), self.iUm))
            dofsA = remove_dofs(dofsAll, self.iUs)
            dofsA.sort()
            U = zeros(n, 'float64')
            print("U   = ", U)
            print("iUs = ", self.iUs)
            #print("iUm = ", self.iUm)

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
        #nnodes = model.grid.n


        # bars
        #ncbars = len(cbars)   # not tested

        # not implemented - solids
        #nctetra4s = model.elements_solid.ctetra4.n
        #ncpenta6s = model.elements_solid.cpenta6.n
        #nchexa8s  = model.elements_solid.chexa8.n

        #=========================
        if self.is_stress or self.is_strain or self.is_force:
            # SPRINGS
            nsprings = 0
            elementTypes = [
                model.elements_spring.celas1,
                model.elements_spring.celas2,
                model.elements_spring.celas3,
                model.elements_spring.celas4,
            ]
            for elementType in elementTypes:
                nsprings += elementType.n

            if nsprings:
                o1 = zeros(nsprings, 'float64')
                e1 = zeros(nsprings, 'float64')
                f1 = zeros(nsprings, 'float64')

                ispring = 0
                for elementType in elementTypes:
                    n = elementType.n
                    #print("n%s = %s" % (elementType.type, n))
                    if n:
                        elementType.displacement_stress(model, self.positions, q, self.nidComponentToID,
                            ispring, o1, e1, f1)
                        eids = elementType.element_id
                        print("eids =", eids)
                    ispring += n
                #if self.is_strain:
                self._store_spring_oes(model, eids, e1, case, elementType.type, Type='strain')
                #if self.is_stress:
                self._store_spring_oes(model, eids, o1, case, elementType.type, Type='stress')
                #if self.is_force:
                self._store_spring_oef(model, eids, f1, case, elementType.type)
                del e1
                del o1
                del f1
            #del elementType model.elements_springs

            # RODS
            elementTypes = [model.crod, model.conrod]  # model.ctube

            for elementType in elementTypes:
                n = elementType.n
                if n:
                    #margin_1 =
                    #margin_2 =
                    #margin_12 =
                    (e1, e4,
                     o1, o4,
                     f1, f4) = elementType.displacement_stress(model, self.positions, q, self.nidComponentToID)
                    eids = elementType.element_id
                    if self.is_strain:
                        self._store_rod_oes(model, eids, e1, e4, case, elementType.type, Type='strain')
                    del e1, e4
                    if self.is_stress:
                        self._store_rod_oes(model, eids, o1, o4, case, elementType.type, Type='stress')
                    del o1, o4
                    if self.is_force:
                        self._store_rod_oef(model, eids, f1, f4, case, elementType.type)
                    del f1, f4
                del elementType

            #=========================
            # CHSEAR
            #ncshears = model.elements_shell.cshear.n # half implemented
            ncshears = 0
            if ncshears:
                #stress = zeros((ncshears, 3), 'float64')
                #strain = zeros((ncshears, 3), 'float64')
                #force  = zeros((ncshears, 16), 'float64')

                stress, strain, force = elementType.displacement_stress(model, self.positions, q, self.nidComponentToID)
                #if self.is_strain:
                self._store_cshear_oes(model, cshears, strain, case, 'CSHEAR', Type='strain')
                #if self.is_stress:
                self._store_cshear_oes(model, cshears, stress, case, 'CSHEAR', Type='stress')
                #if self.is_force:
                self._store_cshear_oef(model, cshears, force, case, 'CSHEAR')
                del stress
                del strain
                del force

            #=========================
            # BARS
            ncbars = 0
            if ncbars:
                print("ncbars", ncbars)
                o1 = zeros(ncbars, 'float64')
                e1 = zeros(ncbars, 'float64')
                f1 = zeros(ncbars, 'float64')
                for i, eid in enumerate(cbars):
                    element = elements[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID)
                    o1[i] = oxi
                    e1[i] = exi
                    f1[i] = fxi
                #if self.is_strain:
                self._store_bar_oes(model, cbars, e1, case, Type='strain')
                #if self.is_stress:
                self._store_bar_oes(model, cbars, o1, case, Type='stress')
                #if self.is_force:
                self._store_bar_oef(model, cbars, f1, case)
                del e1
                del o1
                del f1

            #=========================
            # BEAMS
            ncbeams = 0
            if ncbeams:
                print("ncbeams", ncbeams)
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
            #ncquad4s = model.elements_shell.cquad4.n
            #nctria3s = model.elements_shell.ctria3.n
            nctria3s = 0
            ncquad4s = 0
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
        return
        results = [self.displacements]
        header = None
        pageStamp = None

        self._op2_header(f, packing=packing)

        for result in results:
            for subcase, case in sorted(result.iteritems()):
                case.write_op2(header, pageStamp, f, is_mag_phase=False, packing=packing)
                if not packing:
                    f.write('\n')
        marker1 = [4, 0, 4]
        marker2 = [4, 0, 4]
        marker3 = [4, 0, 4]
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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}
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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}

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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}
        return data_code

    def _store_cshear_oef(self, model, eids, force, case, elementType):
        if len(eids) == 0:
            return
        data_code = self._OEF_f06_header(case, elementType)
        is_sort1 = True
        isubcase = case.id

        if elementType == 'CSHEAR':
            forces = RealCShearForce(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(elementType)

        data = []
        i = 0
        #(elementID,
        #   f12, f14, tau1,
        #   ...) = line
        for (eid, Fi) in zip(eids, force):
            line = [eid] + list(Fi)
            data.append(line)

        dt = None
        forces.add_f06_data(data, dt)

        if elementType == 'CSHEAR':
            self.shearForces[isubcase] = forces
        else:
            raise NotImplementedError(elementType)
        #stress.dt = None

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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}

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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}

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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}
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
                    'nonlinear_factor': None, 'dataNames':['lsdvmn']}
        disp = DisplacementObject(data_code, is_sort1, isubcase, dt=None)

        data = []

        i = 0
        #(nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
        for ni in xrange(model.grid.n):
            nid = model.grid.nid[ni]
            line = [nid, 'G']
            xyz = U[i:i + 6]  # 1,2,3,4,5,6
            print("nid=%s txyz,rxyz=%s" % ( nid, xyz))
            line += xyz
            i += 6
            data.append(line)
        disp.add_f06_data(data, dt)
        self.displacements[isubcase] = disp
        self.iSubcases.append(isubcase)

    def setup_sol_101(self, model, case):
        # the (GridID,componentID) -> internalID
        (self.nidComponentToID, i) = self.build_nid_component_to_id(model)
        self.apply_SPCs(model, case, self.nidComponentToID)
        self.apply_MPCs(model, case, self.nidComponentToID)

        #spcDOFs = self.iUs
        #mpcDOFs = self.iUm

        #Mgg = zeros((i, i), 'float64')
        Kgg, Kgg_sparse = self.assemble_global_stiffness(model, i, self.nidComponentToID)
        Fg = self.assemble_forces(model, i, case, self.nidComponentToID)
        return Kgg, Fg, i

    def solve_sol_101(self, Kgg, Fg):
        for (i, j, a) in zip(self.iUm, self.jUm, self.Um):
            print("Kgg[%s, %s] = %s" % (i, j, a))
            Kgg[i, j] = a

        (self.IDtoNidComponents) = reverse_dict(self.nidComponentToID)
        print("IDtoNidComponents = ", self.IDtoNidComponents)
        print("Kgg =\n" + print_annotated_matrix(Kgg, self.IDtoNidComponents, self.IDtoNidComponents))
        #print("Kgg = \n", Kgg)
        #print("iSize = ", i)

        #(Kaa, Fa) = self.Partition(Kgg)
        #sys.exit('verify Kgg')

        print("Kgg/%g = \n%s" % (self.knorm, print_matrix(Kgg / self.knorm)))
        Kaa, dofs2 = partition_dense_symmetric(Kgg, self.iUs)
        print("Kaa/%g = \n%s" % (self.knorm, print_matrix(Kaa / self.knorm)))
        #print("Kaa.shape = ",Kaa.shape)

        #sys.exit('verify Kaa')
        Fa, _dofs2 = partition_dense_vector(Fg, self.iUs)
        #print("Kaa = \n%s" % (print_matrix(Kaa)))

        print("Fg/%g = %s" % (self.fnorm, Fg/self.fnorm))
        print("Fa/%g = %s" % (self.fnorm, Fa/self.fnorm))
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

    def add_stiffness(self, K, dofs, nijv):
        Kgg = self.Kgg
        #print(type(Kgg))
        for i, dof1 in enumerate(dofs):
            #dof1i = self.nidComponentToID[dof1]
            for j, dof2 in enumerate(dofs):
                if abs(K[i, j]) > 0.0:
                    print('i=%s j=%s dof1=%s dof2=%s Ke[i,j]=%s' % (i, j, dof1, dof2, K[i,j]/self.knorm))
                    #dof2i = self.nidComponentToID[dof2]
                    #assert isinstance(dof1i, int), dof1i
                    #assert isinstance(dof2i, int), dof2i
                    Kgg[dof1, dof2] += K[i, j]
                    #Kgg[dof1i, dof2i] += K[i, j]
                    #print('Kgg[%i,%i]=%d' % (dof1i, dof2i, Kgg[dof1i, dof2i]) )
        print('Kggi =\n', Kgg)

    def assemble_global_stiffness(self, model, i, Dofs):
        self.Kgg = zeros((i, i), 'float64')
        print("Kgg.shape", self.Kgg.shape)

        dof_mapper = []

        nnodes = model.grid.n
        nspoints = model.spoint.n
        assert nnodes > 0
        print("nnodes =", nnodes)
        #ndofs = 6 * nnodes + nspoints

        i = 0
        nids = model.grid.node_id
        #nids.sort()
        print('nids =', nids, i)
        #dof_ids = arange(0, 6*nnodes, 6)

        #dofs_0 = [nid=2, 1] -> searchsorted(nids, nid)[0]

        #for nid in sorted(self.nodes.iterkeys()):
            #nid_dof_mapper[]

        self.positions = {}
        index0s = {}
        for i in xrange(model.grid.n):
            nid = model.grid.node_id[i]
            self.positions[nid] = model.grid.xyz[i]
            index0s[nid] = 6 * i

        # spring
        for i in xrange(model.elements_spring.celas1.n):
            K, dofs, nijv = model.elements_spring.celas1.get_stiffness(i, model, self.positions, index0s)
            print("Kcelas1 =\n", K)
            self.add_stiffness(K, dofs, nijv)

        for i in xrange(model.elements_spring.celas2.n):
            K, dofs, nijv = model.elements_spring.celas2.get_stiffness(i, model, self.positions, index0s)
            self.add_stiffness(K, dofs, nijv)

        if 0:
            for i in xrange(model.elements_spring.celas3.n):
                K, dofs, nijv = model.elements_spring.celas3.get_stiffness(i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

            for i in xrange(model.elements_spring.celas4.n):
                K, dofs, nijv = model.elements_spring.celas4.get_stiffness(i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        # conrod
        for i in xrange(model.conrod.n):
            K, dofs, nijv = model.conrod.get_stiffness(i, model, self.positions, index0s)
            self.add_stiffness(K, dofs, nijv)

        # crod
        for i in xrange(model.crod.n):
            K, dofs, nijv = model.crod.get_stiffness(i, model, self.positions, index0s)
            self.add_stiffness(K, dofs, nijv)

        # ctube

        # shells
        for i in xrange(model.elements_shell.ctria3.n):
            K, dofs, nijv = model.ctria3.get_stiffness(i, model, self.positions, index0s)
            self.add_stiffness(K, dofs, nijv)

        for i in xrange(model.elements_shell.cquad4.n):
            K, dofs, nijv = model.cquad4.get_stiffness(i, model, self.positions, index0s)
            self.add_stiffness(K, dofs, nijv)

        #Kgg_sparse = coo_matrix((entries, (rows, cols)), shape=(i, i))
        Kgg_sparse = None
        Kgg = self.Kgg
        return Kgg, Kgg_sparse

    def apply_SPCs(self, model, case, nidComponentToID):
        has_spcs = False
        if not case.has_parameter('SPC'):
            spc_ids = self.model.get_SPCx_ids(exclude_spcadd=True)
            has_spcs = True
            ## todo:  is this correct???
        else:
            # get the value, 1 is the options (SPC has no options)
            spc_ids = [ case.get_parameter('SPC')[0] ]

        if case.has_parameter('SPC') or has_spcs:
            for spc_id in spc_ids:
                self.log.debug('applying SPC=%i' % spc_id)
                SpcSet = model.SPC(spc_id)

                #print("spc_set =", SpcSet)
                for spc in SpcSet:
                    if spc.type == 'SPC1':
                        for dof, node_ids in spc.components.iteritems():
                            #print("dof =", dof)
                            for dofi in dof:
                                dofi = int(dofi)
                                for nid in node_ids:
                                    key = (nid, dofi)
                                    i = nidComponentToID[key]
                                    print("i=%s Us=%s" % (i, 0.0))
                                    if i not in self.iUsb:
                                        self.iUsb.append(i)
                                        self.Usb.append(0.0)
                                    #else:
                                        #raise RuntimeError('duplicate ')
                    elif spc.type == 'SPC':
                        for dof, node_id in spc.components:
                            key = (node_id, dof)
                            i = nidComponentToID[key]
                            if i not in self.iUsb:
                                self.iUsb.append(i)
                                self.Usb.append(0.0)
                    else:
                        raise NotImplementedError(spc.type)

    def apply_MPCs(self, model, case, nidComponentToID):
        isMPC = False
        mp_index = self.mp_index
        if case.has_parameter('MPC'):
            # get the value, 1 is the options (MPC has no options)
            mpc_id = case.get_parameter('MPC')[0]
            mpcs = model.MPC(mpc_id)

            iconstraint = mp_index
            for mpc in mpcs:
                if mpc.type == 'MPC':
                    for constraints in mpc.constraints:
                        i = mp_index + iconstraint
                        for (G, C, A) in constraints:
                            key = (G, C)
                            j = nidComponentToID[key]
                            iconstraint += 1

                            self.Ump.append(A)
                            self.iUmp.append(i)
                            self.jUmp.append(j)
                else:
                    raise NotImplementedError(mpc.type)

    def build_dof_sets(self):
        # s = sb + sg
        self.Us  = self.Usb  + self.Usg
        self.iUs = self.iUsb + self.iUsg

        # l = b + c + lm
        self.Ul  = self.Uc  + self.Ulm
        self.iUl = self.iUc + self.iUlm

        # t = l + r
        self.Ut  = self.Ul  + self.Ur
        self.iUt = self.iUl + self.iUr

        # a = t + q
        self.Ua  = self.Ut  + self.Uq
        self.iUa = self.iUt + self.iUq

        # d = a + e
        self.Ud  = self.Ua  + self.Ue
        self.iUd = self.iUa + self.iUe

        # f = a + o
        self.Uf  = self.Ua  + self.Uo
        self.iUf = self.iUa + self.iUo

        # fe = f + e
        self.Ufe  = self.Uf  + self.Ue
        self.iUfe = self.iUf + self.iUe

        # n = f + s
        self.Un  = self.Uf  + self.Us
        self.iUn = self.iUf + self.iUs

        # ne = n + e
        self.Une  = self.Un  + self.Ue
        self.iUne = self.iUn + self.iUe

        # m = mp + mr
        self.Um  = self.Ump  + self.Umr
        self.iUm = self.iUmp + self.iUmr
        self.jUm = self.jUmp + self.jUmr

        # g = n + m
        self.Ug  = self.Un  + self.Um
        self.iUg = self.iUn + self.iUm

        # p = g + e
        self.Up  = self.Ug  + self.Ue
        self.iUp = self.iUg + self.iUe

        # ks = k + sa
        self.Uks  = self.Uk  + self.Usa
        self.iUks = self.iUk + self.iUsa

        # js = j + sa
        self.Ujs  = self.Uj  + self.Usa
        self.iUjs = self.iUj + self.iUsa

        # fr = o + l = f - q - r
        self.Ufr  = self.Uo  + self.Ul
        self.iUfr = self.iUo + self.iUl

        # v = o + c + r
        self.Uv  = self.Uo  + self.Uc  + self.Ur
        self.iUv = self.iUo + self.iUc + self.iUr
        return

    def assemble_forces(self, model, i, case, Dofs):
        """very similar to writeCodeAster loads"""
        Fg = zeros(i, 'float64')
        #print(model.loads)
        (loadID, junk) = model.caseControlDeck.get_subcase_parameter(case.id, 'LOAD')
        print("loadID = ", loadID)
        loads = model.loadcase.resolve(int(loadID))

        for load in loads:
            print(load)
            if load.type in ['FORCE', 'MOMENT']:
                if load.type in ['FORCE']:
                    ni = 0
                elif  load.type in ['MOMENT']:
                    ni = 3
                else:
                    raise NotImplementedError(load.type)

                for i in xrange(load.n):
                    nid = load.node_id[i]
                    x, y, z = load.mag[i] * load.xyz[i]

                    if abs(x) > 0.:
                        Fg[Dofs[(nid, ni + 1)]] += x
                    if abs(y) > 0.:
                        Fg[Dofs[(nid, ni + 2)]] += y
                    if abs(z) > 0.:
                        Fg[Dofs[(nid, ni + 3)]] += z
            else:
                raise NotImplementedError(load.type)
        return Fg

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

        assert case.has_parameter('SPCFORCES') == True
        if case.has_parameter('SPCFORCES'):
            (value, options) = case.get_parameter('SPCFORCES')
            if options is not []:
                if value != 'NONE':
                    SPCForces = Ksa * Ua + Kss * Us
                    if isMPC:
                        SPCForces += Ksm * Um

                    result = SPCForcesObject(data_code, transient)
                    result.add_f06_data()

                    flag = 0
                    if 'PRINT' in options:
                        f06.write(result.write_f06(header, pageStamp, pageNum))
                        flag += 1
                    if 'PLOT' in options:
                        op2.write(result.write_op2(Title, Subtitle))
                        flag += 1
                    if not flag:
                        f06.write(result.write_f06(header, pageStamp, pageNum))

        if case.has_parameter('MPCFORCES'):
            if options is not []:
                (value, options) = case.get_parameter('MPCFORCES')
                if value != 'NONE':
                    MPCForces = Kma * Ua + Kmm * Um
                    if isSPC:
                        MPCForces += Kms * Us

                    result = MPCForcesObject(data_code, transient)
                    result.add_f06_data()
                    flag = 0
                    if 'PRINT' in options:
                        f06.write(result.write_f06(header, pageStamp, pageNum))
                        flag += 1
                    if 'PLOT' in options:
                        f06.write(result.write_op2(Title, Subtitle))
                        flag += 1
                    if not flag:
                        f06.write(result.write_f06(header, pageStamp, pageNum))

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
                      'MPC','MPCADD',

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

def main():
    from pyNastran.bdf2.solver.solver_args import run_arg_parse
    fargs = run_arg_parse()

    s = Solver(fargs)
    s.run_solver()

    #if os.path.exists(s.op2_name):
        #op2_name = s.op2_name
        #op2 = OP2(op2FileName=op2_name, make_geom=False, debug=True, log=None)
    #else:
        #op2_name = 'solid_bending.op2'
        #op2 = OP2(op2_name)
        #op2.make_op2_debug = True
        #op2.read_op2()

if __name__ == '__main__':
    main()