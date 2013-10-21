# pylint: disable=E0602,C0103
from __future__ import print_function
import os
import sys
from itertools import izip

# 3rd party
from numpy import array, zeros, ones
from numpy import searchsorted
from numpy.linalg import solve

from scipy.sparse import coo_matrix

# pyNastran
from pyNastran.utils.mathematics import print_matrix, print_annotated_matrix
from pyNastran.bdf.bdf import BDF, SPC, SPC1
from pyNastran.f06.f06 import F06
from pyNastran.op2.op2 import OP2

# Tables
from pyNastran.op2.tables.oug.oug_displacements import DisplacementObject
#from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import SPCForcesObject
#from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import MPCForcesObject

# rods
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RodStressObject, RodStrainObject, ConrodStressObject, ConrodStrainObject, CtubeStressObject, CtubeStrainObject
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealRodForce, RealConrodForce, RealCtubeForce

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

    def solve(self, K, F, dofs):  # can be overwritten
        r"""solves \f$ [K]{x} = {F}\f$ for \f${x}\f$"""
        print("--------------")
        print("Kaa_norm = \n" + str(K / self.fnorm))
        print("--------------")
        print("Fa = ", F)
        if F[0] == 0.0:
            assert max(F) != min(F), 'no load is applied...'
        print("--------------")

        #asdf
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
                    faileds.append(i)
                    failed.append([nid, dof])
            msg = self.make_grid_point_singularity_table(failed)
            self.f06_file.write(msg)

            if 'AUTOSPC' in self.model.params:
                autospc = self.model.params['AUTOSPC']
                value = autospc.values[0]

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
                    U[failed] = U2
            else:
                self.f06_file.close()
                raise

        return U

    def run_solver(self, bdfName):
        bdfName = os.path.abspath(bdfName)
        bdf_base, ext = os.path.splitext(bdfName)
        self.f06_name = bdf_base + '.f06'
        self.op2_name = bdf_base + '.op2'

        print("**is3D =", self.is3D)

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
        Kgg, Kgg_sparse = self.assemble_global_stiffness(model, i, Kgg, nidComponentToID, Fg)
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
        # bars
        crods = []
        conrods = []
        ctubes = []

        # bars
        cbeams = []
        cbars = []

        # shells
        cquad4s = []
        ctria3s = []

        # solids
        ctetra4s = []
        cpenta5s = []
        chexa8s = []

        type_map = {
            'CONROD' : conrods,
            'CROD'   : crods,
            'CTUBE'  : ctubes,
            'CBEAM'  : cbeams,
            'CBAR'   : cbars,

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

        # rods
        ncrods = len(crods)
        nconrods = len(conrods)
        nctubes = len(ctubes)

        # bars
        ncbars = len(cbars)
        ncbeams = len(cbeams) # half implemented

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

        # shells
        ctria3s = array(ctria3s)
        cquad4s = array(cquad4s)

        #=========================
        if self.is_stress or self.is_strain or self.is_force:
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
                         f1i, f4i) = element.displacement_stress(model, q, self.nidComponentToID, is3D=self.is3D)
                        o1[i] = o1i
                        e1[i] = e1i
                        f1[i] = f1i

                        o4[i] = o4i
                        e4[i] = e4i
                        f4[i] = f4i
                    if self.is_strain:
                        self.store_rod_oes(model, eids, e1, e4, case, elementType, Type='strain')
                    del e1, e4
                    if self.is_stress:
                        self.store_rod_oes(model, eids, o1, o4, case, elementType, Type='stress')
                    del o1, o4
                    if self.is_force:
                        self.store_rod_oef(model, eids, f1, f4, case, elementType)
                    del f1, f4

            #=========================
            # BARS / BEAMS
            print("ncbeams", ncbeams)
            print("cbeams", cbeams)
            if ncbeams:
                o1 = zeros(ncbeams, 'float64')
                e1 = zeros(ncbeams, 'float64')
                f1 = zeros(ncbeams, 'float64')
                for i, eid in enumerate(cbeams):
                    element = elements[eid]
                    (exi, oxi, fxi) = element.displacement_stress(model, q, self.nidComponentToID, is3D=self.is3D)
                    o1[i] = oxi
                    e1[i] = exi
                    f1[i] = fxi
                #if self.is_strain:
                self.store_beam_oes(model, cbeams, e1, case, Type='strain')
                #if self.is_stress:
                self.store_beam_oes(model, cbeams, o1, case, Type='stress')
                #if self.is_force:
                self.store_beam_oef(model, cbeams, f1, case)
                del e1
                del o1
                del f1

            #=========================
            # SHELLS
            print("nctria3", nctria3s)
            if nctria3s or ncquad4s:
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
                self.store_plate_oes(model, cbeams, stress, case, Type='strain')
                #if self.is_stress:
                self.store_plate_oes(model, cbeams, strain, case, Type='stress')
                #if self.is_force:
                self.store_plate_oef(model, cbeams, force, case)
                del stress, strain, force

            # SOLIDS
        #=========================
        self.write_f06(self.f06_file, end_flag=True)


    def store_beam_oes(self, model, eids, axial, case, elementType='CBEAM', Type='strain'):
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


    def store_beam_oef(self, model, eids, fx, case, elementType='CBEAM'):
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

    def store_rod_oef(self, model, eids, axial, torsion, case, elementType):
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
                    'element_name': elementType, 'format_code':format_code,
                    #'s_code': s_code,
                    'nonlinear_factor': None}

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

    def store_rod_oes(self, model, eids, axial, torsion, case, elementType, Type='stress'):
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
        Ua = self.solve_sol_101(Kaa, Fa, dofs2)
        return Ua, i

    def solve_sol_101(self, Kaa, Fa, dofs2):
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

        Ua = self.solve(Kaa, Fa, dofs2)
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
                (K, dofs, nIJV, Fg, nGrav) = elem.Stiffness(model, node_ids, index0s, self.gravLoad, self.is3D, self.fnorm)
            except TypeError as e:
                msg = 'elem %s must take:\n>>>Stiffness(model, node_ids, index0s, self.gravLoad, self.is3D, self.fnorm)' % elem.type
                print(msg)
                raise
            print("K[%s] = \n%s" % (eid, K / self.fnorm))
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
            if not (isinstance(typesFound, list) or  isinstance(typesFound, set)):
                raise RuntimeError(type(typesFound))
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
                    print(moment)
                    if abs(moment[0]) > 0.:
                        Fg[Dofs[(nid, 4)]] += moment[0]
                    if abs(moment[1]) > 0.:
                        Fg[Dofs[(nid, 5)]] += moment[1]
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