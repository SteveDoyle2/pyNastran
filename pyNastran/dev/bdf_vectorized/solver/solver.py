# -*- coding: utf-8 -*-
# pylint#: disable=C0103
"""
http://slideplayer.com/slide/3330177/
"""
import os
from datetime import date
from struct import pack

import numpy as np
from numpy import (array, zeros, ones, arange,
                   searchsorted, diag)
from numpy.linalg import solve  # type: ignore

import scipy
from scipy.sparse import dok_matrix  # type: ignore
from cpylog import get_logger2

# pyNastran
import pyNastran.bdf.bdf_interface.dev.matrices
from pyNastran.bdf.bdf_interface.dev.mass import make_gpwg
from pyNastran.dev.bdf_vectorized.solver.utils import (
    reverse_dict, partition_dense_symmetric, partition_dense_vector, remove_dofs)
#from pyNastran.f06.f06_writer import sorted_bulk_data_header
from pyNastran.utils.dev import list_print
from pyNastran.utils.mathematics import print_matrix, print_annotated_matrix
from pyNastran.dev.bdf_vectorized.bdf import BDF #, SPC, SPC1
#from pyNastran.f06.f06_writer import F06Writer
from pyNastran.op2.op2 import OP2

# Tables
#from pyNastran.op2.tables.opg_appliedLoads.opg_objects import (
    #RealAppliedLoadsVectorArray, AppliedLoadsVectorArray)

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
#from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject
#from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject
#from pyNastran.f06.dev.tables.oload_resultant import Resultant

# springs
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (
    RealSpringStressArray, RealSpringStrainArray)
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealSpringForceArray

# rods
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import (
    RealRodStressArray, RealRodStrainArray)
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealRodForceArray

# shear
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import (
    RealShearStressArray, RealShearStrainArray)
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealCShearForceArray

# beams
from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (
    RealBeamStressArray, RealBeamStrainArray)
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealCBeamForceArray


from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import (
    RealLoadVectorArray,
)

def partition_dense_matrix(a, b, c=None):
    raise NotImplementedError('partition_dense_matrix a=%s b=%s c=%s' % (str(a), str(b), str(c)))


class Resultant:
    def __init__(self, table_name: str, total_load, subcase_id: int):
        pass
        #Resultant('OLOAD', total_load, self.subcase_id)

class Solver(OP2):
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
      .. todo:: SETx not supported for output requests; ALL/NONE

    Bulk Data:
      GRID,CORDx
       - Position
       - ps constraint
       .. todo:: analysis & output coordinate system != 0

      CONROD, CROD/PROD, CTUBE/PTUBE
       - .. todo:: CTUBE not tested

      CELAS1, CELAS2, CELAS3, CELAS4, PELAS
       - ge not supported (it's damping, so not an issue)
       - .. todo:: test CELAS3/CELAS4

      CSHEAR/PSHEAR
       - .. todo:: what loads can be applied?
       - .. todo:: non-constant a, b not tested
       - .. todo:: f1, f2 not supported and not checked
         .. todo:: support PLOAD2
      MAT1

      LOAD, FORCE, MOMENT
       - coord 0
       - .. todo:: coordx not tested
      PLOAD1
       - .. todo:: validate...
       - distributed load; forces/moments
       - LE: x1=0.0*L; x2=1.0*L
       - FR: x1=0.0; x2=1.0
       - LE_PR:  .. todo:: add this...
       - FR_PR:  .. todo:: add this...
       - .. todo:: support alternate coordinate system
       - .. todo:: static load at x1=0.5

      SPC, SPC1
       - .. todo:: constraints in alternate coordinate system (specified by GRID cards)

      MPC/RBE2/RBE3/RBAR
       - .. todo:: not supported

     CQUAD4/CTRIA3
       - .. todo:: not supported
       - .. todo:: PLOAD, PLOAD2, PLOAD4

     CHEXA, CPENTA, CTETRA
       - .. todo:: not supported
       - .. todo:: PLOAD, PLOAD4

    Results:
      CONROD, CROD/PROD
       - FORCE
       - STRESS/STRAIN
         - margins not supported
         - ge not supported (it's damping, so not an issue)

      CELAS1, CELAS2, CELAS3, CELAS4, PELAS
       - STRESS/STRAIN/FORCE
       - .. todo:: test CELAS3/CELAS4

      CBEAM
          - PLOAD1
        - .. todo:: not done... but close
      CSHEAR
        - FORCE/STRESS/STRAIN
          - not calculated
          - no tables created
    """
    def __init__(self, fargs, log=None):
        """
        fargs : dict[key] : value
            'BDFNAME' : str
                the bdf filename
            '--k' : float > 1.0
                scale factor for stiffness matrix
            '--m' : float > 1.0
                scale factor for mass matrix
            '--f' : float > 1.0
                scale factor for force array
        """
        #F06Writer.__init_data__(self)
        OP2.__init__(self, debug=False, log=None, debug_file=None) # make_geom=False,
        debug = fargs['--debug']
        self.log = get_logger2(log, debug)

        self.page_num = 1
        self.fargs = fargs

        # normalization of stiffness matrix
        self.knorm = fargs['--k']
        # normalization of load vector
        self.fnorm = fargs['--f']
        # normalization of mass matrix
        self.mnorm = fargs['--m']

        self.isubcases = []
        self.nU = 0
        self.nUs = 0
        self.nUm = 0

        #==============================
        # files
        self.f06_filename = None
        self.op2_filename = None
        self.op2_pack_filename = None

        self.f06_file = None
        self.op2_file = None
        self.op2_pack_file = None

        #==============================
        # op2 objects
        self.oload_resultant = None

        #==============================
        # minor data members
        self.Subtitle = None
        self.IDtoNidComponents = None
        self.nidComponentToID = None
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

        # g-set before elimination of any degrees of freedom
        # n-set after elimination of multipoint constraints
        # f-set after elimination of automatic constraints and SPC's
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

        #==============================
        self.is_displacement_result = True
        self.is_stress_result = True
        self.is_strain_result = True
        self.is_force_result = True
        self.positions = {}
        self.subtitle = None
        self.model = None
        self.gravLoad = None
        self.label = None
        self.mp_index = None

        #self.Ujs = None
        #self.Kgg = None
        #self.iUv = None
        #self.iUjs = None
        #self.iUfr = None
        #self.iUne = None
        #self.iUp = None
        #self.IDtoNidComponents = None
        #self.iUf = None
        #self.iUg = None
        #self.iUd = None
        #self.Ua = None
        #self.Uf = None
        #self.Ud = None
        #self.Une = None
        #self.Up = None
        #self.Ufe = None

        #self.Uks = None
        #self.iUks = None
        self.Kgg = None
        self.Mgg = None
        self.Mgg_sparse = None
        #------------------------------
        self.Us = None
        self.iUs = None

        # l = b + c + lm
        self.Ul = None
        self.iUl = None

        # t = l + r
        self.Ut = None
        self.iUt = None

        # a = t + q
        self.Ua = None
        self.iUa = None

        # d = a + e
        self.Ud = None
        self.iUd = None

        # f = a + o
        self.Uf = None
        self.iUf = None

        # fe = f + e
        self.Ufe = None
        self.iUfe = None

        # n = f + s
        self.Un = None
        self.iUn = None

        # ne = n + e
        self.Une = None
        self.iUne = None

        # m = mp + mr
        self.Um = None
        self.iUm = None
        self.jUm = None

        # g = n + m
        self.Ug = None
        self.iUg = None

        # p = g + e
        self.Up = None
        self.iUp = None

        # ks = k + sa
        self.Uks = None
        self.iUks = None

        # js = j + sa
        self.Ujs = None
        self.iUjs = None

        # fr = o + l = f - q - r
        self.Ufr = None
        self.iUfr = None

        # v = o + c + r
        self.Uv = None
        self.iUv = None

    def write_summary(self, f06_file, card_count):
        """dummy function"""
        pass

    def _solve(self, K, F, dofs):  # can be overwritten
        r"""solves \f$ [K]{x} = {F}\f$ for \f${x}\f$"""
        self.log.info("--------------")
        self.log.info("Kaa_norm / %s = \n" % self.knorm + list_print(K / self.knorm))
        self.log.info("--------------")
        self.log.info("Fa/%g = %s" % (self.fnorm, F / self.fnorm))
        if F[0] == 0.0:
            assert max(F) != min(F), 'no load is applied...'
        self.log.info("--------------")

        try:
            U = scipy.sparse.linalg.spsolve(K, F)
            #U = solve(K, F) # numpy
        except Exception:
            failed = []
            faileds = []
            for i, iu in enumerate(dofs):
                #absF = abs(F)
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
                    ilist = list(set(range(len(dofs))).difference(set(faileds)))
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
        self.f06_file.write('-' * 80 + '\n')
        self.f06_file.flush()
        return U

    def run_solver(self):
        """the main interface to the solver"""
        fargs = self.fargs
        bdf_filename = os.path.abspath(fargs['BDFNAME'])
        bdf_base = os.path.abspath(fargs['BDFBASE'])

        #bdf_base, ext = os.path.splitext(bdfName)
        self.f06_filename = bdf_base + '.f06'
        self.op2_filename = bdf_base + '.op2'
        self.op2_pack_filename = bdf_base + '_pack.op2'

        #self.op2_file = open(self.op2_name, 'wb')
        #self.op2_pack_file = open(self.op2_pack_name, 'w')
        self.op2_file = None
        self.op2_pack_file = None

        d = date.today()
        self.date = (d.month, d.day, d.year)
        self.title = 'pyNastran Job'

        #------------------------------------------
        # read the deck
        self.model = BDF(log=self.log, debug=False)
        self.model.cards_to_read = get_solver_cards()
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
        #------------------------------------------
        self.f06_file = open(self.f06_filename, 'w') # , encoding=self.model._encoding

        self.f06_file.write(self.make_f06_header())
        #self.f06_file.write(sorted_bulk_data_header())

        #------------------------------------------
        # start of analysis

        self.log.info(self.model.case_control_deck)
        self.log.info('--------------------')
        for subcase_id, subcase in sorted(self.model.subcases.items()):
            self.log.info(subcase)
            if 'TITLE' in subcase:
                self.title = subcase.get_parameter('TITLE')[0]
                break

        # TODO: there is a bug with the Case Control Deck reading
        #       the TITLE doesn't get into SUBCASE 0 because the line
        #       is missing from case_control_lines
        #
        #assert self.title != 'pyNastran Job'

        page_stamp = self.make_stamp(self.title, self.date)
        self.page_stamp = page_stamp

        cc = self.model.case_control_deck
        #print(cc.subcases)
        analysis_cases = []
        for (isub, subcase) in sorted(cc.subcases.items()):
            self.subcase_key[isub] = [isub]
            self.log.info(subcase)
            if 'LOAD' in subcase:
                analysis_cases.append(subcase)
                #print('analyzing subcase = \n%s' % subcase)
            #else:
                #raise RuntimeError('A LOAD card was not set')

        self.write_summary(self.f06_file, card_count=self.model.card_count)

        #print(analysis_cases)
        for case in analysis_cases:
            self.log.info(case)
            #(value, options) = case.get_parameter('STRESS')
            #print("STRESS value   = %s" % value)
            #print("STRESS options = %s" % options)

            if 'TEMPERATURE(INITIAL)' in case:
                (value, options) = case.get_parameter('TEMPERATURE(INITIAL)')
                self.log.info('value   = %s' % value)
                self.log.info('options = %s' % options)
                raise NotImplementedError('TEMPERATURE(INITIAL) not supported')
                #integrate(B.T*E*alpha*dt*Ads)
            #sys.exit('starting case')
            self.log.info('starting case')
            self.run_case(self.model, case)

        self.f06_file.write('Kgg / %s =\n%s\n\n' % (self.knorm, list_print(self.Kgg / self.knorm)))
        self.f06_file.write('Fg =\n%s\n\n' % list_print(self.Fg))

        self.f06_file.write('Kaa / %s =\n%s\n\n' % (self.knorm, list_print(self.Kaa / self.knorm)))
        self.f06_file.write('Fa =\n%s\n\n' % list_print(self.Fa))

        self.f06_file.close()
        if self.op2_file is not None:
            self.op2_file.close()
            self.op2_pack_file.close()

    def run_case(self, model, case):
        sols = {
            101: self.run_sol_101,
            #103: self.run_sol_103,
        }

        isubcase = case.id
        if model.sol in sols:
            if 'SUBTITLE' in case:
                (self.Subtitle, options) = case.get_parameter('SUBTITLE')
            else:
                self.Subtitle = 'DEFAULT'

            if 'LABEL' in case:
                (self.label, options) = case.get_parameter('LABEL')
            else:
                self.label = ''
            self.isubcase_name_map[isubcase] = [self.Subtitle, self.label]

            # really should be is_f06_stress, is_op2_stress, etc.
            # also should have SET support
            #print(case.params)
            #if 'DISPLACEMENT' in case:
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

            #if 'STRESS' in case:
                #value, options = case.get_parameter('STRESS')
                #if value == 'ALL':
                    #self.is_stress = True
                #elif value == 'NONE':
                    #self.is_stress = False
                #else:
                    #raise NotImplementedError('STRESS = %r is not supported' % value)

            #if 'STRAIN' in case:
                #value, options = case.get_parameter('STRAIN')
                #if value == 'ALL':
                    #self.is_strain = True
                #elif value == 'NONE':
                    #self.is_strain = False
                #else:
                    #raise NotImplementedError('STRAIN = %r is not supported' % value)

            #if 'FORCE' in case:
                #value, options = case.get_parameter('FORCE')
                #if value == 'ALL':
                    #self.is_force = True
                #elif value == 'NONE':
                    #self.is_force = False
                #else:
                    #raise NotImplementedError('FORCE = %r is not supported' % value)

            self.is_displacement_result = True
            self.is_stress_result = True
            self.is_strain_result = True
            self.is_force_result = True
            if not self.case_result_flags['DISPLACEMENT'][0]:
                self.is_displacement_result = False
            if not self.case_result_flags['STRESS'][0]:
                self.is_stress_result = False
            if not self.case_result_flags['STRAIN'][0]:
                self.is_strain_result = False
            if not self.case_result_flags['FORCE'][0]:
                self.is_force_result = False

            if not(self.is_displacement_result or
                   self.is_stress_result or self.is_strain_result or
                   self.is_force_result):
                msg = 'No results selected...'
                raise RuntimeError(msg)
            if model.sol not in sols:
                raise NotImplementedError('sol=%r is not supported' % model.sol)

            self.subcase_id = case.id
            sol = sols[model.sol]
            self.log.info('starting SOL %s' % model.sol)
            sol(model, case)
        else:
            raise NotImplementedError('model.sol=%s not in %s' % (model.sol, sols.keys()))

    def get_case_parameter(self, case, param_name):
        if param_name in case:
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
        nid_component_to_id_map = {}
        if model.grid.n:
            cd = set(model.grid.cd)
            if cd == 0:
                raise NotImplementedError('Set cd=0; cd=%s is not supported.')

        for ni in range(model.grid.n):  # GRIDs
            nid = model.grid.node_id[ni]
            ps = model.grid.ps[ni]

            if ps > -1:
                self.log.info('node_id=%s ps=%s' % (nid, ps))
                iUsg = []
                Usg = []
                j = None
                while ps > 0:
                    psi = ps % 10
                    j = i + psi - 1
                    iUsg.append(j)
                    Usg.append(0.0)
                    ps //= 10
                self.log.info("  Usg=%s" % Usg)
                self.iUsg += iUsg
                self.Usg += Usg

            nid_component_to_id_map[(nid, 1)] = i
            nid_component_to_id_map[(nid, 2)] = i + 1
            nid_component_to_id_map[(nid, 3)] = i + 2
            nid_component_to_id_map[(nid, 4)] = i + 3
            nid_component_to_id_map[(nid, 5)] = i + 4
            nid_component_to_id_map[(nid, 6)] = i + 5
            i += 6
        self.log.info(f'iUsg = {self.iUsg}')

        spoint = model.spoint
        if spoint.n:
            for nid in sorted(model.spoint.spoint):  # SPOINTS
                nid_component_to_id_map[(nid, 1)] = i
                i += 1
        assert i > 0, 'no DOFs'

        #: starting index for MPC cards
        self.mp_index = model.grid.n + spoint.n

        return nid_component_to_id_map, i

    def get_Mgg(self, model, ndofs, force_calcs=False):
        Mgg = None
        if force_calcs:
            Mgg, Mgg_sparse = self.assemble_global_mass_matrix(model, ndofs, self.nidComponentToID)
            model.params['GRDPNT'] = 0

        if 'GRDPNT' in model.params:
            grid_point = model.params['GRDPNT']
            if grid_point == -1:
                self.log.debug('No mass requested')
            else:
                self.log.debug('PARAM,GRIDPNT = %s' % grid_point)
                self.make_gpwg(grid_point, Mgg)
                #sys.exit('check the mass...')
        return Mgg

    def run_sol_103(self, model, case):
        self.end_options['SEMR'] = True
        raise NotImplementedError()
        #"""
        #ug = un+um All structural DOF including scalar DOF
        #um DOF eliminated by multipoint constraints and rigid elements
        #un = uf+us All structural DOF not constrained by multipoint constraints
        #us DOF eliminated by single-point constraints
        #uf = ua+uo Unconstrained (free) structural DOF
        #uo DOF omitted by static condensation (Guyan Reduction)
        #ua = ur+ul DOF used in real eigenvalue analysis
        #ur DOF to which determinate reactions are applied for the solution of free body models
        #ul The remaining structural DOF used in static analysis (points left over)
        #ue Extra DOF introduced in dynamic analysis
        #ud = ua+ue DOF used in dynamic analysis by the direct method
        #up = ug+ue The g-set plus EXTRA points for dynamic analysis
        #uz DOF representing modal coordinates
        #uh = ue+uz DOF used in dynamic analysis by the modal method

        #SPC Forces
        #qs = Kfs * Uf + Kss * Ys
        #qs = Kfs * Uf + Kss * (Ys + Ys_dot*dt + Ys_dotdot*dt*dt)  # where *dt means integration
        #"""
        #assert model.sol == 103, 'model.sol=%s is not 101' % (model.sol)
        #if 'WTMASS' in model.params:
            ## converts weight units to mass
            #wtmass = model.params['WTMASS'].value1
        #else:
            #wtmass = 1.0

        #if 'COUPMASS' in model.params:
            #coupled_mass = model.params['COUPMASS'].value1
            #lumped_mass = False
        #else:
            #coupled_mass = -1
            #lumped_mass = True

        #if 'METHOD' in case:
            #imethod = model.subcase.get_parameter('METHOD')
            #raise NotImplementedError('METHOD')
            ##self.eigb[imethod]
            ##self.eigc[imethod]
            ##self.eigr[imethod]
            ##self.eigrl[imethod]

        #Mgg = self.get_Mgg(model, ndofs, force_calcs=True)
        ## analysis
        #(Kgg, Fg, n) = self.setup_sol_101(model, case)

        #self.build_dof_sets()
        #Lambda, Ua = self.solve_sol_103(Kgg, Mgg)

        #dofsAll = {i for i in range(n)}
        #dofsA = remove_dofs(dofsAll, self.iUs)
        #dofsA.sort()
        #U = zeros(n, 'float64')

        ## TODO handle MPCs
        #for (i, iu) in enumerate(self.iUs):
            #U[iu] = self.Us[i]
        #for (i, iu) in enumerate(dofsA):
            #U[iu] = Ua[i]

        #if self.is_displacement_result:
            #self._store_displacements(model, U, case)
        #q = U


    def run_sol_101(self, model, case):
        #print("case = ", case)
        assert model.sol == 101, 'model.sol=%s is not 101' % model.sol
        self.log.info('case = %s' % case)
        self.log.info('subcase_id = %s' % self.subcase_id)

        #if 'WTMASS' in model.params:
            #wtmass = model.params['WTMASS'].value1
        #else:
            #wtmass = 1.0

        #if 'COUPMASS' in model.params:
            #coupmass = model.params['COUPMASS'].value1
        #else:
            #coupmass = -1

        #if 'FMETHOD' in case:
            #iflutter = model.subcase.get_parameter('FMETHOD')

        #if 'METHOD' in case:
            #imethod = model.subcase.get_parameter('METHOD')
            #self.eigb[imethod]
            #self.eigc[imethod]
            #self.eigr[imethod]
            #self.eigrl[imethod]

        model.params = {'GRDPNT': 0,}
        if "GRDPNT" in model.params and model.params["GRDPNT"] >= 0:
            g0 = model.params["GRDPNT"]
            reference_point = None

            if g0 in model.grid.node_id:
                reference_point = model.grid.get_position_by_node_id(g0)
                #i = model.grid.get_node_index_by_node_id(g0)
                #reference_point = model.grid[g0].get_position()
            #(mass, cg, I) = model.mass_properties(reference_point, sym_axis=None, num_cpus=1)
            #mass *= wtmass
            #I *= wtmass
            # calculate mass - really should use the Mass matrix...

        ## define IDs of grid point components in matrices
        if 1:
            # analysis
            self.log.info('setup_sol_101')
            (Kgg, Fg, n) = self.setup_sol_101(model, case)
            self.build_dof_sets()
            self.log.info("------------------------\n")
            self.log.info("solving...")
            Ua = self.solve_sol_101(Kgg, Fg)
            self.log.info("Ua =\n%s" % Ua)
            self.log.info("Us =\n%s" % self.Us)

            dofsAll = {i for i in range(n)}
            #dofsA = remove_dofs(remove_dofs(dofsAll, self.iUs), self.iUm))
            dofsA = remove_dofs(dofsAll, self.iUs)
            dofsA.sort()
            U = zeros(n, 'float64')
            self.log.info("U   =\n%s" % U)
            self.log.info("iUs =\n%s" % self.iUs)
            #print("iUm = ", self.iUm)

            # TODO handle MPCs
            for (i, iu) in enumerate(self.iUs):
                U[iu] = self.Us[i]
            for (i, iu) in enumerate(dofsA):
                U[iu] = Ua[i]

            self.log.info("*U = \n%s" % U)
            self.log.info("dofsA = %s" % dofsA)

            if self.is_displacement_result:
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
            try:
                unused_val, options = case.get_parameter('OLOAD')
            except KeyError:
                self.log.warning('No OLOAD...')
                #self.log.warning(case)
                #raise
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
        if self.is_stress_result or self.is_strain_result or self.is_force_result:
            # SPRINGS
            nsprings = 0
            element_types = [
                model.celas1,
                model.celas2,
                model.celas3,
                model.celas4,
            ]
            for element_type in element_types:
                nsprings = element_type.n

                if nsprings == 0:
                    continue
                o1 = zeros(nsprings, 'float64')
                e1 = zeros(nsprings, 'float64')
                f1 = zeros(nsprings, 'float64')

                ispring = 0
                for element_type in element_types:
                    n = element_type.n
                    #print("n%s = %s" % (element_type.type, n))
                    if n:
                        element_type.displacement_stress(
                            model, self.positions, q, self.nidComponentToID,
                            ispring, o1, e1, f1)
                        eids = element_type.element_id
                        self.log.info("eids = %s" % eids)

                        self._store_spring_oes(model, eids, e1, case, element_type.type,
                                               Type='strain')
                        #if self.is_stress_result:
                        self._store_spring_oes(model, eids, o1, case, element_type.type,
                                               Type='stress')
                        #if self.is_force_result:
                        self._store_spring_oef(model, eids, f1, case, element_type.type)
                        del e1
                        del o1
                        del f1
                    ispring += n
                    #if self.is_strain:
            #del element_type model.elements_springs

            # RODS
            element_types = [model.crod, model.conrod]  # model.ctube

            self.is_stress_result = True
            self.is_strain_result = True
            self.is_force_result = True
            for element_type in element_types:
                n = element_type.n
                self.log.info('Type=%s n=%s displacement_stress' % (element_type.type, n))
                if n:
                    #margin_1 =
                    #margin_2 =
                    #margin_12 =
                    (e1, e4,
                     o1, o4,
                     f1, f4) = element_type.displacement_stress(
                         model, self.positions, q, self.nidComponentToID)
                    eids = element_type.element_id
                    if self.is_strain_result:
                        self._store_rod_oes(
                            model, eids, e1, e4, case, element_type.type, Type='strain')
                    del e1, e4
                    if self.is_stress_result:
                        self._store_rod_oes(
                            model, eids, o1, o4, case, element_type.type, Type='stress')
                    del o1, o4
                    if self.is_force_result:
                        self._store_rod_oef(model, eids, f1, f4, case, element_type.type)
                    del f1, f4
                del element_type

            #=========================
            # CHSEAR
            #ncshears = model.elements_shell.cshear.n # half implemented
            ncshears = 0
            if ncshears:
                cshears = model.cshear
                #stress = zeros((ncshears, 3), 'float64')
                #strain = zeros((ncshears, 3), 'float64')
                #force  = zeros((ncshears, 16), 'float64')

                stress, strain, force = element_type.displacement_stress(
                    model, self.positions, q, self.nidComponentToID)
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
                cbars = model.cbar
                self.log.info("ncbars = %s" % ncbars)
                o1 = zeros(ncbars, 'float64')
                e1 = zeros(ncbars, 'float64')
                f1 = zeros(ncbars, 'float64')
                for i, eid in enumerate(cbars):
                    element = cbars[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID)
                    o1[i] = oxi
                    e1[i] = exi
                    f1[i] = fxi
                #if self.is_strain_result:
                self._store_bar_oes(model, cbars, e1, case, Type='strain')
                #if self.is_stress_result:
                self._store_bar_oes(model, cbars, o1, case, Type='stress')
                #if self.is_force_result:
                self._store_bar_oef(model, cbars, f1, case)
                del e1
                del o1
                del f1

            #=========================
            # BEAMS
            ncbeams = 0
            if ncbeams:
                cbeams = model.cbeam
                self.log.info("ncbeams = %s" % ncbeams)
                o1 = zeros(ncbeams, 'float64')
                e1 = zeros(ncbeams, 'float64')
                f1 = zeros(ncbeams, 'float64')
                for i, eid in enumerate(cbeams):
                    element = cbeams[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID)
                    o1[i] = oxi
                    e1[i] = exi
                    f1[i] = fxi
                #if self.is_strain_result:
                self._store_beam_oes(model, cbeams, e1, case, Type='strain')
                #if self.is_stress_result:
                self._store_beam_oes(model, cbeams, o1, case, Type='stress')
                #if self.is_force_result:
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
                self.log.info("nctria3 = %s" % nctria3s)
                stress = zeros((nctria3s+ncquad4s, 3), 'float64')
                strain = zeros((nctria3s+ncquad4s, 3), 'float64')
                force = zeros((nctria3s+ncquad4s, 3), 'float64')

            i0 = i
            if nctria3s:
                ctria3s = model.ctria3
                for i, eid in enumerate(ctria3s):
                    element = elements[eid]
                    (stressi, straini, forcei) = element.displacement_stress(
                        model, q, self.nidComponentToID)
                    stress[i, :] = stressi
                    strain[i, :] = straini
                    force[i, :] = forcei
                i0 = i

            if ncquad4s:
                cquad4s = model.cquad4
                for i, eid in enumerate(ncquad4s):
                    element = cquad4s[eid]
                    (stressi, straini, forcei) = element.displacement_stress(
                        model, q, self.nidComponentToID)
                    stress[i0+i, :] = stressi
                    strain[i0+i, :] = straini
                    force[i0+i, :] = forcei

            if nctria3s or ncquad4s:
                #if self.is_strain_result:
                self._store_plate_oes(model, cbeams, stress, case, Type='strain')
                #if self.is_stress_result:
                self._store_plate_oes(model, cbeams, strain, case, Type='stress')
                #if self.is_force_result:
                self._store_plate_oef(model, cbeams, force, case)
                del stress, strain, force

            # SOLIDS
        #=========================
        self.log.debug(self.displacements[1].data)
        #self.log.debug(self.displacements[1])
        self.log.debug(self.conrod_force)
        self.log.debug(self.conrod_stress)
        self.log.debug(self.conrod_strain)
        self.write_f06(self.f06_file, end_flag=True, quiet=True, close=False)
        self.write_op2(self.op2_file, packing=True)
        self.write_op2(self.op2_pack_file, packing=False)
        self.log.info('finished SOL 101')

    def _op2_header(self, op2_file, packing=True):
        data = [
            4, 3, 4,
            1, 28, 12,
            4, 7, 4,   # 7
            'NASTRAN FORT TAPE ID CODE - ',   # 28 = 7*4
            4, 2, 4,   # 7
            4, -1, 4,
            4, 0, 4,
            4, 2, 4,
            4, 0, 4,
            4, 2, 4]
        if packing:
            op2_file.write(pack('9i28s18i', *data))
        if not packing:
            op2_file.write(str(data) + '\n')

    def write_op2(self, op2_file, packing=False):
        return
        #results = [self.displacements]
        #header = None
        #page_stamp = None

        #self._op2_header(op2_file, packing=packing)

        #for result in results:
            #for subcase, case in sorted(result.items()):
                #case.write_op2(header, page_stamp, op2_file, is_mag_phase=False, packing=packing)
                #if not packing:
                    #op2_file.write('\n')
        #marker1 = [4, 0, 4]
        #marker2 = [4, 0, 4]
        #marker3 = [4, 0, 4]
        #marker = marker1 + marker2 + marker3
        #if packing:
            #nmarker = len(marker)
            #p = pack('%ii' % nmarker, *marker)
            #op2_file.write(p)
        #else:
            #op2_file.write(str(marker)+'\n')
        #op2_file.close()

    def _store_beam_oes(self, model, eids, axial, case, element_type='CBEAM', Type='strain'):
        #print('eids =', eids)
        if len(eids) == 0:
            return
        analysis_code = 1
        #transient = False
        isubcase = case.id
        is_sort1 = False
        dt = None
        format_code = 1  # ???
        s_code = 1

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 1, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
            'element_name': element_type, 'format_code':format_code,
            's_code': s_code,
            'nonlinear_factor': None, 'data_names':['lsdvmn']}
        if Type == 'stress':
            if element_type == 'CBEAM':
                stress = RealBeamStressArray(data_code, is_sort1, isubcase, dt=False)
        elif Type == 'strain':
            if element_type == 'CBEAM':
                stress = RealBeamStrainArray(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(Type)

        data = []
        for (eid, axiali) in zip(eids, axial):
            element = model.Element(eid)
            n1, n2 = element.node_ids
            self.log.info(n1, n2)
            #      (eid, grid, sd,  sxc,   sxd, sxe, sxf,  smax, smin, mst, msc) = out
            line = [eid, n1, 0.0, axiali, 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            data.append(line)

            line = [eid, n2, 1.0, axiali, 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            data.append(line)
        stress.add_f06_data(data, dt)

        if element_type == 'CBEAM' and Type == 'stress':
            self.cbeam_stress[isubcase] = stress
        elif element_type == 'CBEAM' and Type == 'strain':
            self.cbeam_strain[isubcase] = stress
        else:
            raise NotImplementedError('element_type=%r Type=%r' % (element_type, Type))
        stress.dt = None


    def _store_beam_oef(self, model, eids, fx, case, element_type='CBEAM'):
        #print('eids =', eids)
        if len(eids) == 0:
            return
        analysis_code = 1
        #transient = False
        isubcase = case.id
        is_sort1 = False
        dt = None
        format_code = 1  # ???
        #s_code = None

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 1, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OEF',
            'element_name': element_type, 'format_code':format_code,
            #'s_code': s_code,
            'nonlinear_factor': None, 'data_names':['lsdvmn']}

        if element_type == 'CBEAM':
            forces = RealCBeamForceArray(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(element_type)

        data = []
        for (eid, fxi) in zip(eids, fx):
            element = model.Element(eid)
            n1, n2 = element.node_ids
            self.log.info('***(*', n1, n2)
            #      [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
            line = [eid, n1, 0.0, 0., 0., 0., 0., 0., 0., 0.0]
            data.append(line)
            line = [eid, n1, 1.0, 0., 0., 0., 0., 0., 0., 0.0]
            data.append(line)
            line = [eid, n2, 0.0, 0., 0., 0., 0., 0., 0., 0.0]
            data.append(line)
            line = [eid, n2, 1.0, 0., 0., 0., 0., 0., 0., 0.0]
            #data.append(line)
        self.log.info(data)
        forces.add_f06_data(data, dt)

        if element_type == 'CBEAM':
            self.cbeam_force[isubcase] = forces
        else:
            raise NotImplementedError(element_type)
        #stress.dt = None

    def _OEF_f06_header(self, case, element_name):
        analysis_code = 1
        #transient = False
        #is_sort1 = False
        #dt = None
        format_code = 1  # ???
        s_code = None

        if element_name == 'CELAS1':
            element_type = 11
        elif element_name == 'CELAS2':
            element_type = 12
        elif element_name == 'CELAS3':
            element_type = 13
        elif element_name == 'CELAS4':
            element_type = 14
        elif element_name == 'CROD':
            element_type = None
        elif element_name == 'CONROD':
            element_type = None
        else:
            raise NotImplementedError(element_name)

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'sort_code': 0,
            'tCode' : 1, 'table_code':1,  # are these the same?
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OEF',
            'element_name': element_name, 'element_type': element_type,
            'format_code':format_code,
            #'s_code': s_code,
            'nonlinear_factor': None, 'data_names':['lsdvmn']}
        return data_code

    def _store_cshear_oef(self, model, eids, force, case, element_name):
        if len(eids) == 0:
            return
        data_code = self._OEF_f06_header(case, element_name)
        is_sort1 = True
        isubcase = case.id

        if element_name == 'CSHEAR':
            forces = RealCShearForceArray(data_code, is_sort1, isubcase, dt=False)
        else:
            raise NotImplementedError(element_name)

        data = []
        #(elementID,
        #   f12, f14, tau1,
        #   ...) = line
        for (eid, Fi) in zip(eids, force):
            line = [eid] + list(Fi)
            data.append(line)

        dt = None
        forces.add_f06_data(data, dt)

        if element_name == 'CSHEAR':
            self.cshear_force[isubcase] = forces
        else:
            raise NotImplementedError(element_name)
        #stress.dt = None

    def _store_bar_oes(self, model, cbars, e1, case, Type='strain'):
        raise NotImplementedError()

    def _store_bar_oef(self, model, cbars, f1, case):
        raise NotImplementedError()

    def _store_plate_oes(self, model, cbeams, stress, case, Type='strain'):
        raise NotImplementedError()

    def _store_plate_oef(self, model, cbeams, force, case):
        raise NotImplementedError()

    def _store_cshear_oes(self, model, eids, results, case, element_type, Type='strain'):
        """
        self.cshear_forces = {}
        self.cshear_stress = {}
        self.cshear_strain = {}

                    S T R E S S E S   I N   S H E A R   P A N E L S   ( C S H E A R )
      ELEMENT       MAX            AVG        SAFETY    ELEMENT    MAX      AVG    SAFETY
        ID.        SHEAR          SHEAR       MARGIN      ID.     SHEAR    SHEAR   MARGIN
          328   1.721350E+03   1.570314E+03   7.2E+01
        """
        if len(eids) == 0:
            return
        analysis_code = 1
        #transient = False
        isubcase = case.id
        is_sort1 = True
        dt = None
        format_code = 1  # ???
        s_code = None

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 1, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
            'element_name': element_type, 'format_code':format_code,
            's_code': s_code,
            'nonlinear_factor': None, 'data_names':['lsdvmn']}

        if Type == 'stress':
            #if element_type == 'CELAS2':
            stress = RealShearStressArray(data_code, is_sort1, isubcase, dt=None)
            #else:
                #raise NotImplementedError(element_type)
        elif Type == 'strain':
            #if element_type == 'CELAS2':
            stress = RealShearStrainArray(data_code, is_sort1, isubcase, dt=None)
            #else:
                #raise NotImplementedError(element_type)
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
            self.cshear_stress[isubcase] = stress
        elif Type == 'strain':
            self.cshear_strain[isubcase] = stress
        else:
            raise NotImplementedError('element_type=%r Type=%r' % (element_type, Type))
        #stress.dt = None

    def _store_spring_oes(self, model, eids, axial, case, element_name, Type='stress'):
        if len(eids) == 0:
            return
        analysis_code = 1
        #transient = False
        isubcase = case.id
        is_sort1 = True
        #dt = None
        format_code = 1  # ???
        s_code = None

        stress_code = 0
        if Type == 'strain':
            stress_code = 1

        if element_name == 'CELAS1':
            element_type = 11
        elif element_name == 'CELAS2':
            element_type = 12
        elif element_name == 'CELAS3':
            element_type = 13
        elif element_name == 'CELAS4':
            element_type = 14
        else:
            raise NotImplementedError(element_name)

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'sort_code': 0,
            'tCode' : 1, 'table_code':1,  # are these the same?
            'sort_bits': [0, 0, 0],
            'num_wide': 8, 'table_name': 'OES',
            'element_name': element_name, 'element_type': element_type,
            'format_code':format_code,
            's_code': s_code,
            'nonlinear_factor': None, 'data_names':['lsdvmn'],
            'stress_bits' : [None, stress_code, None, stress_code, None],
        }

        if Type == 'stress':
            result = RealSpringStressArray(data_code, is_sort1, isubcase, dt=None)
        elif Type == 'strain':
            result = RealSpringStrainArray(data_code, is_sort1, isubcase, dt=None)
        else:
            raise NotImplementedError(Type)
        #result.element_type = element_type

        #print('axial =', axial)
        ntimes = 1
        nelements = eids.size
        dtype = 'float32'
        result.build_data(ntimes, nelements, dtype)
        result.data[0, :, 0] = axial
        result.element = eids

        #print('axial %s = %s' % (Type, axial))
        if Type == 'stress':
            if element_name == 'CELAS1':
                self.celas1_stress[isubcase] = result
            elif element_name == 'CELAS2':
                self.celas2_stress[isubcase] = result
            elif element_name == 'CELAS3':
                self.celas3_stress[isubcase] = result
            elif element_name == 'CELAS4':
                self.celas4_stress[isubcase] = result
            else:
                raise NotImplementedError('element_name=%r Type=%r' % (element_name, Type))
        elif Type == 'strain':
            if element_name == 'CELAS1':
                self.celas1_strain[isubcase] = result
            elif element_name == 'CELAS2':
                self.celas2_strain[isubcase] = result
            elif element_name == 'CELAS3':
                self.celas3_strain[isubcase] = result
            elif element_name == 'CELAS4':
                self.celas4_strain[isubcase] = result
            else:
                raise NotImplementedError('element_name=%r Type=%r' % (element_name, Type))
        else:
            raise NotImplementedError('element_name=%r Type=%r' % (element_name, Type))
        #stress.dt = None

    def _store_spring_oef(self, model, eids, axial, case, element_name):
        if len(eids) == 0:
            return
        #print('storing spring force...')
        data_code = self._OEF_f06_header(case, element_name)

        is_sort1 = True
        isubcase = case.id
        #dt = None
        forces = RealSpringForceArray(data_code, is_sort1, isubcase, dt=None)

        ntimes = 1
        nelements = eids.size
        dtype = 'float32'
        forces.build_data(ntimes, nelements, dtype)
        forces.data[0, :, 0] = axial
        forces.element = eids

        if element_name == 'CELAS1':
            self.celas1_force[isubcase] = forces
        elif element_name == 'CELAS2':
            self.celas2_force[isubcase] = forces
        elif element_name == 'CELAS3':
            self.celas3_force[isubcase] = forces
        elif element_name == 'CELAS4':
            self.celas4_force[isubcase] = forces
        else:
            raise NotImplementedError(element_name)
        #stress.dt = None

    def _store_rod_oef(self, model, eids, axial, torsion, case, element_type):
        if len(eids) == 0:
            return
        data_code = self._OEF_f06_header(case, element_type)
        is_sort1 = True
        isubcase = case.id

        forces = RealRodForceArray(data_code, is_sort1, isubcase, dt=False)

        #data = []
        #i = 0
        #(elementID, axial, torsion) = line
        #for (eid, axiali, torsioni) in zip(eids, axial, torsion):
            #line = [eid, axiali, torsioni]
            #data.append(line)

        #dt = None
        #forces.add_f06_data(data, dt)
        ntimes = 1
        nelements = eids.size
        forces.build_data(ntimes, nelements, float_fmt='float32')
        forces.data[0, :, 0] = axial
        forces.data[0, :, 1] = torsion
        forces.element = eids

        if element_type == 'CROD':
            self.crod_force[isubcase] = forces
        elif element_type == 'CONROD':
            self.conrod_force[isubcase] = forces
        elif element_type == 'CTUBE':
            self.ctube_force[isubcase] = forces
        else:
            raise NotImplementedError(element_type)
        #stress.dt = None

    def _store_rod_oes(self, model, eids, axial, torsion, case, element_type, Type='stress'):
        #if len(eids) == 0:
            #return
        analysis_code = 1
        #transient = False
        isubcase = case.id
        is_sort1 = True
        #dt = None
        format_code = 1  # ???
        s_code = None

        stress_code = 0
        if Type == 'strain':
            stress_code = 1

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'sort_code': 0,
            'tCode' : 1, 'table_code':1,  # are these the same?
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OES',
            'element_name': element_type, 'format_code':format_code,
            's_code': s_code,
            'nonlinear_factor': None, 'data_names':['lsdvmn'],
            'stress_bits' : [None, stress_code, None, stress_code, None],
        }

        nelements = len(eids)
        data = np.zeros((nelements, 4), dtype='float32')

        #(elementID, axial, torsion, margin_axial, margin_torsion) = line
        data[:, 0] = axial
        data[:, 2] = torsion

        ntimes = 1
        nelements = eids.size
        dtype = 'float32'
        if Type == 'stress':
            result = RealRodStressArray(data_code, is_sort1, isubcase, dt=False)
            stress = self.op2_results.stress
            if element_type == 'CROD':
                stress.crod_stress[isubcase] = result
            elif element_type == 'CONROD':
                stress.conrod_stress[isubcase] = result
            elif element_type == 'CTUBE':
                stress.ctube_stress[isubcase] = result
            else:
                raise NotImplementedError('element_type=%r Type=%r' % (element_type, Type))
        elif Type == 'strain':
            result = RealRodStrainArray(data_code, is_sort1, isubcase, dt=False)
            strain = self.op2_results.strain
            if element_type == 'CROD':
                strain.crod_strain[isubcase] = result
            elif element_type == 'CONROD':
                strain.conrod_strain[isubcase] = result
            elif element_type == 'CTUBE':
                strain.ctube_strain[isubcase] = result
            else:
                raise NotImplementedError('element_type=%r Type=%r' % (element_type, Type))
        else:
            raise NotImplementedError('element_type=%r Type=%r' % (element_type, Type))
        result.build_data(ntimes, nelements, dtype)
        result.data[0, :, :] = data
        result.element = eids

        #stress.write_f06(self.header, pageStamp)
        #stress.dt = None

    def _store_displacements(self, model, U, case):
        """
        fills the displacement object
        """
        self.isubcases = []
        analysis_code = 1
        #transient = False
        isubcase = case.id
        is_sort1 = True
        #dt = None
        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OUGV1',
            'nonlinear_factor': None, 'data_names':['lsdvmn'],
            'tCode' : 1, 'table_code':1,  # are these the same?
        }
        ntimes = 1
        nnodes = model.grid.n
        disp = RealDisplacementArray(data_code, is_sort1, isubcase, dt=None)

        ntotal = nnodes
        disp.build_data(ntimes, nnodes, ntotal,
                        ntimes, nnodes, float_fmt='float32')
        #data = []

        #i = 0
        #(nodeID, gridType, t1, t2, t3, r1, r2, r3) = line

        #grid_type = 0  # SECTOR/HARMONIC/RING POINT; H
        #grid_type = 1  # GRID; G
        #grid_type = 2  # SPOINT; S
        #grid_type = 7  # RIGID POINT (e.g. RBE3); L
        disp.node_gridtype[:, 0] = model.grid.node_id
        disp.node_gridtype[:, 1] = 1 # GRID (TODO: no SPOINTs)
        disp.data[0, :, :] = U.reshape(1, nnodes, 6)
        self.displacements[isubcase] = disp
        self.isubcases.append(isubcase)

    def setup_sol_101(self, model, case):
        # the (GridID,componentID) -> internalID
        (self.nidComponentToID, ndofs) = self.build_nid_component_to_id(model)
        #self.log.info('apply SPCs')
        #self.apply_SPCs(model, case, self.nidComponentToID)
        #self.log.info('apply MPCs')
        #self.apply_MPCs(model, case, self.nidComponentToID)

        #spc_dofs = self.iUs
        #mpc_dofs = self.iUm

        #self.log.info('building Mgg')
        #Mgg = self.get_Mgg(model, ndofs, force_calcs=True)

        ngrids = model.grid.n
        xyz_cid0 = np.zeros((ngrids, 3), dtype='float64')
        xyz = model.grid.xyz  ## TODO: update for coords
        for i in range(ngrids):
            xyz_cid0[i, :] = xyz[i, :]

        self.log.info('building Fg')
        Fg = self.assemble_forces(model, ndofs, case, self.nidComponentToID, xyz_cid0)

        self.log.info('building Kgg')
        Kgg, Kgg_sparse = self.assemble_global_stiffness_matrix(model, ndofs, self.nidComponentToID)

        self.log.info('ready to run...')
        return Kgg, Fg, ndofs

    def make_gpwg(self, grid_point, Mgg):
        """
        Parameters
        ----------
        grid_point : int
            0->origin, x>0, that grid point
        Mgg : (N, N) matrix
            the mass matrix
        """
        #return
        example = False
        if example:
            # nid, xb, yb, zb
            node1 = [0., 0., 0.]
            node2 = [1., 0., 0.]
            node3 = [0.5, 1., 0.]
            node4 = [0.5, 0., 1.]
            xyz_cid0 = array([node1, node2, node3, node4],
                             dtype='float32')
            cg = xyz_cid0.mean(axis=0)
            self.log.info('cg = %s' % cg)
        else:
            xyz_cid0 = zeros((self.model.grid.n, 3), dtype='float32')
            for i, (key, xyz) in enumerate(sorted(self.positions.items())):
                xyz_cid0[i, :] = xyz

            #xyz_cid0 = self.model.grid.get_xyz(cid=0)

        grid_point = 2
        if grid_point in [-1, 0]:
            # -1 is don't compute a mass matrix, but since we're in this function,
            # we use the origin.
            ref_point = zeros(3, dtype='float32')
        else:
            # find mass/inertia about point G
            ref_point = self.positions[grid_point]

        grid_cps = self.model.grid.cp
        coords = self.model.coords
        Mo, S, mass, cg, II, IQ, Q = make_gpwg(
            Mgg, ref_point, xyz_cid0, grid_cps, coords, self.log)

        self.grid_point_weight.reference_point = grid_point
        self.grid_point_weight.MO = Mo
        self.grid_point_weight.S = S
        self.grid_point_weight.mass = mass
        self.grid_point_weight.cg = cg
        self.grid_point_weight.IS = II
        self.grid_point_weight.IQ = diag(IQ)
        self.grid_point_weight.Q = Q



    def _save_applied_load(self, Fg):
        data_code = {
            'nonlinear_factor':None,
            'sort_bits':[0, 0, 0], 'analysis_code':0,
            'is_msc':True, 'format_code':0, 'sort_code': 0,
            'data_names':['lsdvmn'],
            'num_wide': 8, 'table_name': 'OPG',
            'tCode' : 1, 'table_code':1,  # are these the same?

            #'log': self.log, 'analysis_code': analysis_code,
            #'device_code': 1, 'table_code': 1, 'sort_code': 0,
            #'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OUG',
            #'nonlinear_factor': None, 'data_names':['lsdvmn']
        }
        isubcase = self.subcase_id
        dt = None
        is_sort1 = True
        applied_loads = RealLoadVectorArray(
            data_code,
            is_sort1,
            isubcase,
            dt)
        ntimes = 1
        nnodes = Fg.size // 6
        nx = ntimes
        ny = nnodes
        float_fmt = 'float32'
        ntotal = nnodes
        applied_loads.build_data(ntimes, nnodes, ntotal, nx, ny, float_fmt)
        applied_loads.node_gridtype[:, 0] = self.model.grid.node_id
        applied_loads.node_gridtype[:, 1] = 1 # G
        self.load_vectors[self.subcase_id] = applied_loads
        applied_loads.data[0, :, :] = Fg.reshape((nnodes, 6))

    def solve_sol_101(self, Kgg, Fg):
        self.Kgg = Kgg
        self.Fg = Fg

        self._save_applied_load(Fg)
        for (i, j, a) in zip(self.iUm, self.jUm, self.Um):
            self.log.info("Kgg[%s, %s] = %s" % (i, j, a))
            Kgg[i, j] = a

        self.IDtoNidComponents = reverse_dict(self.nidComponentToID)
        self.log.info("IDtoNidComponents = %s" % self.IDtoNidComponents)
        self.log.info("Kgg =\n" + print_annotated_matrix(Kgg, self.IDtoNidComponents,
                                                         self.IDtoNidComponents))
        #print("Kgg = \n", Kgg)
        #print("iSize = ", i)

        #(Kaa, Fa) = self.Partition(Kgg)
        #sys.exit('verify Kgg')

        self.log.info("Kgg/%g = \n%s" % (self.knorm, print_matrix(Kgg / self.knorm)))
        Kaa, dofs2 = partition_dense_symmetric(Kgg, self.iUs)
        self.log.info("Kaa/%g = \n%s" % (self.knorm, print_matrix(Kaa / self.knorm)))
        #print("Kaa.shape = ",Kaa.shape)

        #sys.exit('verify Kaa')
        Fa, _dofs2 = partition_dense_vector(Fg, self.iUs)
        #print("Kaa = \n%s" % (print_matrix(Kaa)))

        self.log.info("Fg/%g = \n%s" % (self.fnorm, Fg/self.fnorm))
        self.log.info("Fa/%g = \n%s" % (self.fnorm, Fa/self.fnorm))
        #print("Us = ", self.Us)

        #self.Us = array(self.Us, 'float64')  # SPC
        #self.Um = array(self.Um, 'float64')  # MPC

        #is_spc = False
        #is_mpc = False
        #if 0:
            #zero = array([])
            #MPCgg = zero
            #Ksa = Kas = Kss = Cam = Cma = Kma = Kam = Kmm = Kms = Ksm = zero
            #Kaa1 = Kaa2 = zero
            #Fm = Fs = zero

            ##Kaa = partition_dense_matrix(Kgg,iFree)
            #Kaa0 = Kaa
            #Fa0 = Fa

            #if is_spc:
               ##Fs  = partition_dense_vector(Fg, self.iUs)
                #Ksa = partition_dense_matrix(Kgg, self.iUs, iFree)
                #Kas = Ksa.transpose()
                #Kss = partition_dense_matrix(Kgg, self.iUs)

            #if is_mpc:
                #Fm = partition_dense_vector(Fg, self.iUm)
                #Cam = partition_dense_matrix(MPCgg, iFree)
                #Cma = partition_dense_matrix(MPCgg, self.iUm)

                #Kaa1 = Cam * Kmm * Cma
                #Kaa2 = Cam * Kma + Kam * Cma
                #assert Cam.transpose() == Cma

                #Kma = partition_dense_matrix(Kgg, self.iUm, iFree)
                #Kam = Kma.transpose()
                #Kmm = partition_dense_matrix(Kgg, self.iUm)
                #if is_spc:
                    #Kms = partition_dense_matrix(Kgg, self.iUm, self.iUs)
                    #Ksm = Kms.transpose()

            #Fa = Fa0  # + Cam*Fm
            #Kaa = Kaa0  # +Kaa1+Kaa2

        self.Kaa = Kaa
        self.Fa = Fa
        Ua = self._solve(Kaa, Fa, dofs2)
        #self.Um = Kma*Ua
        self.Ua = Ua
        return Ua

    def element_dof_start(self, elem, nids):
        node_ids = elem.node_ids
        index0s = searchsorted(nids, node_ids)
        index0s *= 6
        return node_ids, index0s

    def add_stiffness(self, K, dofs, nijv):
        #print('self.IDtoNidComponents = %s' % self.IDtoNidComponents) # None
        self.log.debug('self.nidComponentToID')
        for k, value in sorted(self.nidComponentToID.items()):
            self.log.debug('  %s : %s' % (k, value))
        #nidComponentToID

        Kgg = self.Kgg
        kabs = np.abs(K)
        #print(type(Kgg))
        self.log.debug('Ki =\n\n%s' % K)

        for i, dof1 in enumerate(dofs):
            if isinstance(dof1, tuple):
                dof1i = self.nidComponentToID[dof1]
            else:
                dof1i = dof1

            for j, dof2 in enumerate(dofs):
                if kabs[i, j] > 0.0:
                    self.log.debug('  i=%s j=%s dof1=%s dof2=%s Ke[i,j]=%s' % (
                        i, j, dof1, dof2, K[i, j] / self.knorm))

                    if isinstance(dof2, tuple):
                        dof2i = self.nidComponentToID[dof2]
                    else:
                        dof2i = dof2
                    #assert isinstance(dof1i, int), dof1i
                    #assert isinstance(dof2i, int), dof2i
                    Kgg[dof1i, dof2i] += K[i, j]
                    #Kgg[dof1i, dof2i] += K[i, j]
                    #print('Kgg[%s,%s]=%s' % (dof1, dof2, str(Kgg[dof1i, dof2i])) )
        self.log.info('Kggi =\n%s' % Kgg)

    def add_mass(self, M, dofs, nijv):
        Mgg = self.Mgg
        Mgg_sparse = self.Mgg_sparse
        for i, dof1 in enumerate(dofs):
            dof1i = self.nidComponentToID[dof1]
            for j, dof2 in enumerate(dofs):
                if abs(M[i, j]) > 0.0:
                    dof2i = self.nidComponentToID[dof2]
                    Mgg[dof1i, dof2i] += M[i, j]
                    Mgg_sparse[dof1i, dof2i] += M[i, j]

    def assemble_global_stiffness_matrix(self, model, i, Dofs):
        self.Kgg = zeros((i, i), 'float64')
        self.log.info("Kgg.shape = %s" % str(self.Kgg.shape))

        #dof_mapper = []

        nnodes = model.grid.n
        #nspoints = model.spoint.n
        assert nnodes > 0
        self.log.info("nnodes = %s" % nnodes)
        #ndofs = 6 * nnodes + nspoints

        i = 0
        nids = model.grid.node_id
        #nids.sort()
        self.log.info('nids = %s %s' %(nids, i))
        #dof_ids = arange(0, 6*nnodes, 6)

        #dofs_0 = [nid=2, 1] -> searchsorted(nids, nid)[0]

        #for nid in sorted(self.nodes.keys()):
            #nid_dof_mapper[]

        self.log.info('start calculating xyz_cid0')
        self.positions = {}
        index0s = {}
        for i in range(model.grid.n):
            nid = model.grid.node_id[i]
            self.positions[nid] = model.grid.xyz[i]
            index0s[nid] = 6 * i
        self.log.info('end calculating xyz_cid0')

        # spring
        if model.celas1.n:
            self.log.info('start calculating Kcelas1')
            for i in range(model.celas1.n):
                K, dofs, nijv = model.celas1.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.log.info("Kcelas1 =\n%s" % K)
                self.add_stiffness(K, dofs, nijv)

        if model.celas2.n:
            self.log.info('start calculating Kcelas2')
            for i in range(model.celas2.n):
                K, dofs, nijv = model.celas2.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        if model.celas3.n:
            self.log.info('start calculating Kcelas3')
            for i in range(model.celas3.n):
                K, dofs, nijv = model.celas3.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        if model.celas4.n:
            self.log.info('start calculating Kcelas4')
            for i in range(model.celas4.n):
                K, dofs, nijv = model.celas4.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        # conrod
        if model.conrod.n:
            self.log.info('start calculating Kconrod')
            for i in range(model.conrod.n):
                #self.conrod_stress = {}
                K, dofs, nijv = model.conrod.get_stiffness_matrix(i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        # crod
        if model.crod.n:
            self.log.info('start calculating Kcrod')
            for i in range(model.crod.n):
                K, dofs, nijv = model.crod.get_stiffness_matrix(i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        # ctube
        if model.ctube.n:
            self.log.info('start calculating Kctube')
            for i in range(model.crod.n):
                K, dofs, nijv = model.ctube.get_stiffness_matrix(i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        # shells
        if model.ctria3.n:
            for i in range(model.ctria3.n):
                K, dofs, nijv = model.ctria3.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        if model.cquad4.n:
            for i in range(model.cquad4.n):
                K, dofs, nijv = model.cquad4.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)

        # solids
        if model.ctetra4.n:
            #K, dofs, nijv = model.ctetra4.get_stiffness_matrices(
            #    model, xyz_cid0, index0s)
            for i in range(model.ctetra4.n):
                K, dofs, nijv = model.ctetra4.get_stiffness_matrix(
                    i, model, self.positions, index0s)
                self.add_stiffness(K, dofs, nijv)


        #Kgg_sparse = coo_matrix((entries, (rows, cols)), shape=(i, i))
        Kgg_sparse = None
        Kgg = self.Kgg
        return Kgg, Kgg_sparse

    #def assemble_global_damping_matrix(self, model, i, Dofs):

    def assemble_global_mass_matrix(self, model, ndofs, Dofs):
        self.Mgg = zeros((ndofs, ndofs), 'float64')
        self.Mgg_sparse = dok_matrix((ndofs, ndofs), dtype='float64')

        #dof_mapper = []

        nnodes = model.grid.n
        #nspoints = model.spoint.n
        assert nnodes > 0, nnodes

        i = 0
        #nids = model.grid.node_id

        self.positions = {}
        index0s = {}
        for i in range(model.grid.n):
            nid = model.grid.node_id[i]
            self.positions[nid] = model.grid.xyz[i]
            index0s[nid] = 6 * i

        # mass
        conm1 = model.mass.conm1
        for i in range(conm1.n):
            M = conm1.get_mass_matrix(i)
            i0 = index[conm1.node_id[i]]
            coord_id = conm1.coord_id[i]
            if coord_id != 0:
                msg = 'CONM1 doesnt support coord_id != 0 for element %i; coord_id=%i' % (
                    model.conm1.element_id[i], coord_id)
                raise RuntimeError(msg)
            # CONM1 doesn't consider coord ID
            self.Mgg[i0:i0+6, i0:i0+6] += M[:, :]
            for ii1, i01 in arange(i0, i0+6):
                for ii2, i02 in arange(i0, i0+6):
                    self.Mgg_sparse[i01, i02] = M[ii1, ii2]
        for i in range(model.mass.conm2.n):
            M, dofs, nijv = model.conm1.get_mass_matrix(i, model, self.positions, index0s)
            self.add_mass(M, dofs, nijv)

        #if 0:
            ## cmass
            #for i in range(model.cmass1.n):
                #M, dofs, nijv = model.cmass1.get_mass_matrix(i, model, self.positions, index0s)
                #self.add_mass(M, dofs, nijv)
            #for i in range(model.cmass2.n):
                #M, dofs, nijv = model.cmass2.get_mass_matrix(i, model, self.positions, index0s)
                #self.add_mass(M, dofs, nijv)
            #for i in range(model.cmass3.n):
                #M, dofs, nijv = model.cmass3.get_mass_matrix(i, model, self.positions, index0s)
                #self.add_mass(M, dofs, nijv)
            #for i in range(model.cmass4.n):
                #M, dofs, nijv = model.cmass4.get_mass_matrix(i, model, self.positions, index0s)
                #self.add_mass(M, dofs, nijv)

        # conrod
        for i in range(model.conrod.n):
            M, dofs, nijv = model.conrod.get_mass_matrix(i, model, self.positions, index0s)
            self.add_mass(M, dofs, nijv)
        # crod
        for i in range(model.crod.n):
            M, dofs, nijv = model.crod.get_mass_matrix(i, model, self.positions, index0s)
            self.add_mass(M, dofs, nijv)
        #if 0:
            ## ctube
            #for i in range(model.ctube.n):
                #M, dofs, nijv = model.ctube.get_mass_matrix(i, model, self.positions, index0s)
                #self.add_mass(M, dofs, nijv)

        # shells
        for i in range(model.elements_shell.ctria3.n):
            M, dofs, nijv = model.elements_shell.ctria3.get_mass_matrix(
                i, model, self.positions, index0s)
            self.add_mass(M, dofs, nijv)
        for i in range(model.elements_shell.cquad4.n):
            M, dofs, nijv = model.elements_shell.cquad4.get_mass_matrix(
                i, model, self.positions, index0s)
            self.add_mass(M, dofs, nijv)

        # solids
        nsolids = model.elements_solid.n
        if nsolids:
            solid = model.elements_solid
            for i in range(solid.ctetra4.n):
                M, dofs, nijv = solid.ctetra4.get_mass_matrix(i, model, self.positions, index0s)
                self.add_mass(M, dofs, nijv)
            for i in range(solid.cpenta6.n):
                M, dofs, nijv = solid.ctetra4.get_mass_matrix(i, model, self.positions, index0s)
                self.add_mass(M, dofs, nijv)
            for i in range(solid.chexa8.n):
                M, dofs, nijv = solid.ctetra4.get_mass_matrix(i, model, self.positions, index0s)
                self.add_mass(M, dofs, nijv)
        # ctetra10
        # cpenta15
        # chexa20

        #Mgg_sparse = coo_matrix((entries, (rows, cols)), shape=(i, i))
        Mgg_sparse = None
        Mgg = self.Mgg
        self.log.info('returning Mgg')
        return Mgg, Mgg_sparse

    def apply_SPCs(self, model, case, nidComponentToID):
        has_spcs = False
        if not 'SPC' in case:
            spc_ids = self.model.get_SPCx_ids(exclude_spcadd=True)
            has_spcs = True
            ## todo:  is this correct???
        else:
            # get the value, 1 is the options (SPC has no options)
            spc_ids = [case.get_parameter('SPC')[0]]

        if case.has_parameter('SPC') or has_spcs:
            for spc_id in spc_ids:
                self.log.debug('applying SPC=%i' % spc_id)
                spc_set = model.SPC(spc_id)

                #print("spc_set =", spc_set)
                for spc in spc_set:
                    if spc.type == 'SPC1':
                        for dof, node_ids in spc.components.items():
                            #print("dof =", dof)
                            for dofi in dof:
                                dofi = int(dofi)
                                for nid in node_ids:
                                    key = (nid, dofi)
                                    i = nidComponentToID[key]
                                    self.log.info("  i=%s Us=%s" % (i, 0.0))
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
        return has_spcs

    def apply_MPCs(self, model, case, nidComponentToID):
        is_mpc = False
        mp_index = self.mp_index
        if 'MPC' in case:
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
        return is_mpc

    def build_dof_sets(self):
        # s = sb + sg
        self.Us = self.Usb + self.Usg
        self.iUs = self.iUsb + self.iUsg

        # l = b + c + lm
        self.Ul = self.Uc + self.Ulm
        self.iUl = self.iUc + self.iUlm

        # t = l + r
        self.Ut = self.Ul + self.Ur
        self.iUt = self.iUl + self.iUr

        # a = t + q
        self.Ua = self.Ut + self.Uq
        self.iUa = self.iUt + self.iUq

        # d = a + e
        self.Ud = self.Ua + self.Ue
        self.iUd = self.iUa + self.iUe

        # f = a + o
        self.Uf = self.Ua + self.Uo
        self.iUf = self.iUa + self.iUo

        # fe = f + e
        self.Ufe = self.Uf + self.Ue
        self.iUfe = self.iUf + self.iUe

        # n = f + s
        self.Un = self.Uf + self.Us
        self.iUn = self.iUf + self.iUs

        # ne = n + e
        self.Une = self.Un + self.Ue
        self.iUne = self.iUn + self.iUe

        # m = mp + mr
        self.Um = self.Ump + self.Umr
        self.iUm = self.iUmp + self.iUmr
        self.jUm = self.jUmp + self.jUmr

        # g = n + m
        self.Ug = self.Un + self.Um
        self.iUg = self.iUn + self.iUm

        # p = g + e
        self.Up = self.Ug + self.Ue
        self.iUp = self.iUg + self.iUe

        # ks = k + sa
        self.Uks = self.Uk + self.Usa
        self.iUks = self.iUk + self.iUsa

        # js = j + sa
        self.Ujs = self.Uj + self.Usa
        self.iUjs = self.iUj + self.iUsa

        # fr = o + l = f - q - r
        self.Ufr = self.Uo + self.Ul
        self.iUfr = self.iUo + self.iUl

        # v = o + c + r
        self.Uv = self.Uo + self.Uc + self.Ur
        self.iUv = self.iUo + self.iUc + self.iUr
        return

    def write_oload_resultant(self, Fg, xyz_cid0):
        nrows = Fg.size // 6
        Fxyz = Fg.reshape(nrows, 6)
        Fxyz_sum = Fxyz.sum(axis=0)
        forces = Fxyz[:, :3]

        #print('Fxyz =', Fxyz)
        #print('Fxyz_sum =', Fxyz_sum)
        #print('forces =', forces)
        #print('Fg =\n', Fg)
        #print('xyz_cid0 =\n', xyz_cid0)
        #print(xyz_cid0.shape)
        #print(forces.shape)
        moments = np.cross(xyz_cid0, forces)
        #print('moments =', moments)
        self.log.debug(Fg)
        moments_sum = moments.sum(axis=0)
        #print('moments_sum =', moments_sum)
        total_load = Fxyz_sum
        total_load[:3] += moments_sum
        self.oload_resultant = Resultant('OLOAD', total_load, self.subcase_id)
        page_stamp = self.page_stamp
        page_num = self.page_num
        self.oload_resultant.write_f06(self.f06_file, page_stamp, page_num)
        #print(self.oload_resultant)
        #self.f06_file.flush()


    def assemble_forces(self, model, ndofs, case, Dofs, xyz_cid0):
        """builds loads"""
        self.log.info('assemble forces')
        Fg = zeros(ndofs, 'float64')
        #print(model.loads)
        load_id = model.case_control_deck.get_subcase_parameter(case.id, 'LOAD')[0]
        self.log.info("load_id = %s" % load_id)
        loads = model.loadcase.resolve(int(load_id))

        for load in loads:
            self.log.info(load)
            if load.type in ['FORCE', 'MOMENT']:
                if load.type in ['FORCE']:
                    ni = 0
                elif  load.type in ['MOMENT']:
                    ni = 3
                else:
                    raise NotImplementedError(load.type)

                for i in range(load.n):
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
        #print('Fg = %s' % Fg)
        self.write_oload_resultant(Fg, xyz_cid0=xyz_cid0)
        return Fg

        #self.grav_load = array([0., 0., 0.])
        #for load_set in LoadSet:
            #self.log.info("load_set = %r" % str(load_set))
            ##print("type", type(load_set))
            #(types_found, force_loads, moment_loads,
             #force_constraints, moment_constraints,
             #gravity_load) = load_set.organize_loads(model)

            #self.log.info('types_found = %s' % types_found)
            #if not (isinstance(types_found, list) or  isinstance(types_found, set)):
                #raise RuntimeError(type(types_found))
            #assert isinstance(force_loads, dict), type(force_loads)
            #assert isinstance(moment_loads, dict), type(moment_loads)
            #assert isinstance(force_constraints, dict), type(force_constraints)
            #assert isinstance(moment_constraints, dict), type(moment_constraints)
            #assert isinstance(gravity_load, list), type(gravity_load)

            #nids = set()
            #for nid in force_loads:
                #nids.add(nid)
            #for nid in moment_loads:
                #nids.add(nid)

            #if gravity_load != []:
                #self.log.info("gravity_load = %s" % gravity_load)
                #self.grav_load += gravity_load

            #for nid in nids:
                #self.log.info("nid = %s" % nid)

                #if nid in force_loads:
                    #force = force_loads[nid]
                    ##print("force", force)
                    #if abs(force[0]) > 0.:
                        #Fg[Dofs[(nid, 1)]] += force[0]
                    #if abs(force[1]) > 0.:
                        #Fg[Dofs[(nid, 2)]] += force[1]
                    #if abs(force[2]) > 0.:
                        #Fg[Dofs[(nid, 3)]] += force[2]

                #if nid in moment_loads:
                    #moment = moment_loads[nid]
                    ##print("moment", moment)
                    #if abs(moment[0]) > 0.:
                        #Fg[Dofs[(nid, 4)]] += moment[0]
                    #if abs(moment[1]) > 0.:
                        #Fg[Dofs[(nid, 5)]] += moment[1]
                    #if abs(moment[2]) > 0.:
                        #Fg[Dofs[(nid, 6)]] += moment[2]

        #if sum(abs(self.grav_load)) > 0.0:
            #for (eid, elem) in sorted(model.elements.items()):  # CROD, CONROD
                #self.log.info("----------------------------")
                #node_ids, index0s = self.element_dof_start(elem, nids)
                #self.log.info("node_ids=%s index0s=%s" % (node_ids, index0s))

                ## n_ijv is the position of the values of K in the dof
                #fnorm = 10.
                #(Fgi, nGrav) = elem.Fg(model, self.grav_load, fnorm)
                #for (fg, dof) in zip(Fgi, nGrav):
                    ##print("dof = ",dof)
                    #if dof in Dofs:
                        #Fg[Dofs[dof]] += fg
        #self.log.info("Fg  = %s" % Fg)
        #return Fg

    #def write_results(self, case):
        #Us = self.Us  # SPC set
        #Um = self.Um  # MPC set
        #Ua = self.Ua  # constrained set
        ##Ug - global set

        #iUs = self.iUs
        #iUm = self.iUm
        #iUa = self.iUa
        #page_num = 1

        #if case.has_parameter('DISPLACEMENT'):
            #(value, options) = case.get_parameter('DISPLACEMENT')
            #if options is not []:
                #Ug_separate = [[Ua, iUa], [Us, iUs], [Um, iUm]]
                #Ug = departition_dense_vector(Ug_separate)

                #result = RealDisplacementArray(data_code, transient)
                #result.add_f06_data()

                #if 'PRINT' in options:
                    #f06.write(result.write_f06(header, pageStamp, page_num))
                #if 'PLOT' in options:
                    #op2.write(result.write_op2(self.title, self.Subtitle))

        #assert case.has_parameter('SPCFORCES') == True
        #if case.has_parameter('SPCFORCES'):
            #(value, options) = case.get_parameter('SPCFORCES')
            #if options is not []:
                #if value != 'NONE':
                    #SPCForces = Ksa * Ua + Kss * Us
                    #if isMPC:
                        #SPCForces += Ksm * Um

                    #result = RealSPCForces(data_code, transient)
                    #result.add_f06_data()

                    #flag = 0
                    #if 'PRINT' in options:
                        #f06.write(result.write_f06(header, pageStamp, page_num))
                        #flag += 1
                    #if 'PLOT' in options:
                        #op2.write(result.write_op2(Title, Subtitle))
                        #flag += 1
                    #if not flag:
                        #f06.write(result.write_f06(header, pageStamp, page_num))

        #if case.has_parameter('MPCFORCES'):
            #if options is not []:
                #(value, options) = case.get_parameter('MPCFORCES')
                #if value != 'NONE':
                    #MPCForces = Kma * Ua + Kmm * Um
                    #if is_spc:
                        #MPCForces += Kms * Us

                    #result = RealMPCForces(data_code, transient)
                    #result.add_f06_data()
                    #flag = 0
                    #if 'PRINT' in options:
                        #f06.write(result.write_f06(header, pageStamp, page_num))
                        #flag += 1
                    #if 'PLOT' in options:
                        #f06.write(result.write_op2(Title, Subtitle))
                        #flag += 1
                    #if not flag:
                        #f06.write(result.write_f06(header, pageStamp, page_num))

        #if case.has_parameter('GPFORCE'):
            #if options is not []:
                #(value, options) = case.get_parameter('GPFORCE')
                #AppliedLoads = Kaa * Ua
                #if is_spc:
                    #AppliedLoads += Kas * Us
                #if is_mpc:
                    #AppliedLoads += Kam * Um

                #result = AppliedLoadsObject(data_code, transient)
                #result.add_f06_data()
                #if 'PRINT' in options:
                    #f06.write(result.write_f06(header, pageStamp, page_num))
                #if 'PLOT' in options:
                    #op2.write(result.write_op2(Title, Subtitle))

        #if case.has_parameter('STRAIN'):
            #if options is not []:
                #(value, options) = case.get_parameter('STRAIN')

                #for (eid, elem) in sorted(model.elements()):
                    #pass

                #result = xxxObject(data_code, transient)
                #result.add_f06_data()
                #if 'PRINT' in options:
                    #f06.write(result.write_f06(header, pageStamp, page_num))
                #if 'PLOT' in options:
                    #op2.write(result.write_op2(Title, Subtitle))


def get_solver_cards():
    cards_to_read = set([
        'PARAM',
        'GRID', 'GRDSET',

        # elements
        'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',

        'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CONROD',  #'CBEND',
        'CTRIA3', 'CTRIA6',
        'CQUAD4', 'CQUAD8',
        'CTETRA', 'CPENTA', 'CHEXA',
        'CSHEAR',

        # rigid elements - represent as MPCs???
        #'RBAR','RBAR1','RBE1','RBE2','RBE3',

        # properties
        'PELAS',
        'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PTUBE',
        #'PBEND',
        'PSHELL', 'PCOMP', 'PSHEAR',  # 'PCOMPG',
        'PSOLID',

        # materials
        'MAT1', 'MAT2', 'MAT8',

        # spc/mpc constraints
        'SPC', 'SPC1', 'SPCADD',
        'MPC', 'MPCADD',

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
    return cards_to_read

def main():
    from pyNastran.dev.bdf_vectorized.solver.solver_args import run_arg_parse
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

def test_mass():
    fargs = {
        '--k': 1.0,
        '--f': 1.0,
        '--m': 1.0,
    }
    s = Solver(fargs)
    Mgg = None
    grid_point = 0
    s.make_gpwg(grid_point, Mgg)

if __name__ == '__main__':  # pragma: no cover
    #test_mass()
    main()
