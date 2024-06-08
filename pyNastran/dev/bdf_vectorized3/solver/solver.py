"""
References:
https://www.oasys-software.com/help/gsa/9.0/GSA_Theory.pdf
https://download.strusoft.com/FEM-Design/inst110x/theory.pdf
https://people.duke.edu/~hpgavin/StructuralDynamics/StructuralElements.pdf

Newton Rhapson
--------------
https://www.youtube.com/watch?v=uXvoN4OleeE&ab_channel=Dr.ClaytonPettit
[K][x] = {F0}
[Kt][dx] = {F0} - {F_xi}
{x} += {dx}
{F0}: is the starting force
{Kt}: tangent stiffness

force convergence:
|({F_xi} - {F0})| / |{F0}| < Ftol; Ftol=0.05%

displacement convergence:
|{dx}| / |{dxi}| < utol; utol=

"""
import os
import copy
from datetime import date
from collections import defaultdict
from itertools import count
from typing import Union, TextIO, Any

import numpy as np
import scipy as sp
import scipy.sparse as sci_sparse
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.dok import dok_matrix

from cpylog import SimpleLogger
import pyNastran
from pyNastran.nptyping_interface import (
    NDArrayNbool, NDArrayNint, NDArrayN2int, NDArrayNfloat, NDArrayNNfloat)
from pyNastran.utils.numpy_utils import integer_types # , float_types
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase

from pyNastran.f06.f06_writer import make_end
from pyNastran.f06.tables.oload_resultant import Resultant

from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_interface.op2_classes import (
    RealDisplacementArray, RealSPCForcesArray, RealLoadVectorArray,
    RealEigenvalues, RealEigenvectorArray)
from pyNastran.op2.result_objects.grid_point_weight import make_grid_point_weight
#from pyNastran.bdf.mesh_utils.loads import get_ndof

from pyNastran.dev.bdf_vectorized3.solver.recover.static_force import recover_force_101
#from .recover.static_stress import recover_stress_101
#from .recover.static_strain import recover_strain_101
from .recover.strain_energy import recover_strain_energy_101
from .recover.utils import get_plot_request
from .build_stiffness import build_Kgg, DOF_MAP, Kbb_to_Kgg
from pyNastran.dev.bdf_vectorized3.bdf_interface.breakdowns import NO_MASS # , NO_VOLUME


class Solver:
    """defines the Nastran knockoff class"""
    def __init__(self, model: BDF):
        self.model = model
        model.setup(run_geom_check=True)

        self.superelement_id = 0
        self.op2 = OP2(log=model.log, mode='nx')
        self.log = model.log

        d = date.today()
        # the date stamp used in the F06
        self.op2.date = (d.month, d.day, d.year)

        # the "as solved" a-set displacement (after AUTOSPC is applied)
        self.xa_ = None
        self.Kaa = None
        self.Kgg = None
        self.aset = None
        self.sset = None

        base_name = os.path.splitext(model.bdf_filename)[0]
        self._bdf_filename = base_name + '.solver.bdf'
        self.f06_filename = base_name + '.solver.f06'
        self.op2_filename = base_name + '.solver.op2'

    def run(self):
        page_num = 1
        model = self.model
        model.write_bdf(self._bdf_filename)
        sol = model.sol
        solmap = {
            101 : self.run_sol_101_statics,  # static
            103 : self.run_sol_103_modes,  # modes
            111 : self.run_sol_111_harmonic,  # modal frequency response
        }
        model.cross_reference()
        self._update_card_count()

        title = ''
        title = f'pyNastran {pyNastran.__version__}'
        for subcase in model.subcases.values():
            if 'TITLE' in subcase:
                title = subcase.get_parameter('TITLE')
                break

        today = None
        page_stamp = self.op2.make_stamp(title, today) # + '\n'
        with open(self.f06_filename, 'w') as f06_file:
            f06_file.write(self.op2.make_f06_header())
            self.op2._write_summary(f06_file, card_count=model.card_count)
            f06_file.write('\n')
            if sol in [101, 103, 105, 107, 109, 111, 112]:
                for subcase_id, subcase in sorted(model.subcases.items()):
                    if subcase_id == 0:
                        continue
                    self.log.debug(f'subcase_id={subcase_id}')
                    #isubcase = subcase.id
                    subtitle = f'SUBCASE {subcase_id}'
                    label = ''
                    if 'SUBTITLE' in subcase:
                        subtitle = subcase.get_parameter('SUBTITLE')
                    if 'LABEL' in subcase:
                        label = subcase.get_parameter('LABEL')

                    runner = solmap[sol]
                    end_options = runner(
                        subcase, f06_file, page_stamp,
                        title=title, subtitle=subtitle, label=label,
                        page_num=page_num,
                        idtype='int32', fdtype='float64')
            else:
                raise NotImplementedError(sol)
            end_flag = True
            f06_file.write(make_end(end_flag, end_options))

    def _update_card_count(self) -> None:
        for card_type, values in self.model._type_to_id_map.items():
            self.model.card_count[card_type] = len(values)

    def build_Fb(self, xg: NDArrayNfloat, sset_b,
                 dof_map: DOF_MAP,
                 ndof: int, subcase: Subcase) -> NDArrayNfloat:
        model = self.model
        model.log.info('starting build_Fb')

        Fb = np.zeros(ndof, dtype='float32')
        if 'LOAD' not in subcase:
            return Fb

        load_id, unused_options = subcase['LOAD']
        #print('load_id =', load_id)
        reduced_loads = model.get_reduced_static_load()

        #loads, scales, is_grav = model.get_reduced_loads(
            #load_id, scale=1., consider_load_combinations=True,
            #skip_scale_factor0=False, stop_on_failure=True, msg='')
        loads = reduced_loads[load_id]
        #loads : list[loads]
            #a series of load objects
        #scale_factors : list[float]
            #the associated scale factors
        #is_grav : bool
            #is there a gravity card
        for (scale, load) in loads:
            if load.type == 'SLOAD':
                #print(load.get_stats())
                for mag, nid in zip(load.mags, load.nodes):
                    i = dof_map[(nid, 0)]  # TODO: wrong...
                    Fb[i] = mag * scale
            elif load.type == 'FORCE':
                # TODO: wrong because it doesn't handle SPOINTs
                fxyz_myz = load.sum_forces_moments()
                fxyz = fxyz_myz[:, :3]
                nids = load.node_id
                self.log.debug(f'  FORCE nids={nids} Fxyz={fxyz}')
                for fxyzi, nid in zip(fxyz, nids):
                    assert len(fxyzi) == 3
                    fi = dof_map[(nid, 1)]
                    Fb[fi:fi+3] = fxyzi

            elif load.type == 'MOMENT':
                fxyz_myz = load.sum_forces_moments()
                mxyz = fxyz_myz[:, 3:]
                nids = load.node_id
                self.log.debug(f'  MOMENT nid={nids} Mxyz={mxyz}')
                for mxyzi, nid in zip(mxyz, nids):
                    fi = dof_map[(nid, 4)]
                    assert len(mxyzi) == 3
                    Fb[fi:fi+3] = mxyzi
            elif load.type == 'SPCD':
                for nid, components, enforced in zip(load.nodes, load.components, load.enforced):
                    for component in str(components):
                        dof = int(component)
                        fi = dof_map[(nid, dof)]
                        #print(f'(nid, dof) = ({nid}, {dof}) => {fi}')
                        xg[fi] = enforced
                        Fb[fi] = np.nan
                        sset_b[fi] = True
            else:
                print(load.get_stats())
                raise NotImplementedError(load)
        #print(subcase)
        model.log.info('end of build_Fb')
        return Fb

    def get_mpc_constraints(self, subcase: Subcase,
                            dof_map: DOF_MAP,
                            fdtype: str='float64') -> Any:
        ndof = len(dof_map)
        model = self.model
        #Fb = np.zeros(ndof, dtype='float32')
        constraints = []
        if 'MPC' not in subcase:
            return constraints

        log = model.log
        mpc_id, unused_options = subcase['MPC']

        #print('mpc_id =', mpc_id)

        ieq = 0
        ieqs_list = []
        dofs_list = []
        coefficients_list = []
        dependents_list = []
        independents_list = []
        independents_eq = defaultdict(list)
        rigid_element_cards = (card for card in model.rigid_element_cards
                               if card.n > 0)
        for element in rigid_element_cards:
            #etype = element.type
            #ieq += 1
            raise NotImplementedError(element.get_stats())

        mpc = model.mpc
        nequations = len(mpc)
        Cmpc = dok_matrix((nequations, ndof), dtype=fdtype)
        # TODO: assumes a single set (no MPCADD)

        for j, mpc_idi, (idim0, idim1) in zip(count(), mpc.mpc_id, mpc.idim):
            #coefficients : array([ 1., -1.])
            #components : array([1, 3])
            #id     : array([], dtype=int32)
            #idim   : array([[0, 2]])
            #max_id : 10
            #mpc_id : array([10])
            #n      : 1
            #nnode  : array([2])
            #node_id : array([2, 3])
            coefficients = mpc.coefficients[idim0:idim1]
            components = mpc.components[idim0:idim1]
            nodes = mpc.node_id[idim0:idim1]

            # The first degree-of-freedom (G1, C1) in the sequence is defined to
            # be the dependent DOF.
            # A dependent DOF assigned by one MPC entry cannot be assigned
            #  dependent by another MPC entry or by a rigid element.

            #: Component number. (Any one of the Integers 1 through 6 for grid
            #: points; blank or zero for scalar points.)
            #self.components = components

            #: Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
            #self.coefficients = coefficients

            for i, nid, component, coeff in zip(count(), nodes, components, coefficients):
                self.log.debug(f'ieq={ieq} (g,c)={(nid, component)} coeff={coeff}')
                assert isinstance(component, integer_types), component
                idof = dof_map[(nid, component)]
                ieqs_list.append(i)
                Cmpc[j, idof] = coeff

                dofs_list.append(idof)
                coefficients_list.append(coeff)
                if i == 0:
                    dependents_list.append(idof)
                else:
                    independents_list.append(idof)
                    independents_eq[ieq].append(idof)
            ieq += 1

        #print(f'Cmpc         = {Cmpc}')
        dofs = np.array(dofs_list, dtype='int32')
        ieqs = np.array(ieqs_list, dtype='int32')
        independents = np.array(independents_list, dtype='int32')
        dependents = np.array(dependents_list, dtype='int32')
        independents_dependents = np.hstack([independents, dependents])
        coefficients = np.array(coefficients_list, dtype='float64')
        #print(f'ieqs         = {ieqs}')
        #print(f'dofs         = {dofs}')
        #print(f'coefficients = {coefficients}')
        #print(f'dependents   = {dependents}')
        #print(f'independents = {independents}')

        #remove_dependents()
        udependents = np.unique(dependents)
        assert len(udependents) == len(dependents)
        uindependents = np.unique(independents)
        dofs = independents + dependents
        udofs, idofs = np.unique(independents_dependents, return_counts=True)

        while 1:
            for udof in uindependents:
                n_dof = np.where(udof == udofs)[0][0]
                log.debug(f'dof={udof} N={n_dof}')
                if n_dof == 1:
                    break
            else:
                raise RuntimeError('cannot find a DOF to remove')
            #np.where(condition)

            break
        #for udof, in udofs:

        assert nequations > 0, 'No MPC equations despite MPCs being found'

        CmpcI = dok_matrix((nequations, ndof*2), dtype=fdtype)
        CmpcI[:, :ndof] = Cmpc
        for row in range(nequations):
            CmpcI[row, ndof+row] = 1
        log.info(f'CmpcI = {CmpcI}')

        log.info(f'Cmpc = {Cmpc}')
        print('rows:')
        for row in range(nequations):
            dense_row = Cmpc[row, :].toarray()
            abs_row = abs_row = np.abs(dense_row)
            max_val = abs_row.max()
            imax = np.where(dense_row == max_val)[0][0]
            if max_val == 0.:
                continue
            dense_row = CmpcI[row, :].toarray()
            for row2 in range(row + 1, nequations):
                val = CmpcI[row2, imax]
                print(row, row2, imax, val)
                CmpcI[row2, :] -= dense_row * (val / max_val)

            print(dense_row)
            print(f'max={max_val} imax={imax}')
            print(CmpcI)
        #self.log.info(constraints)

        #uindependents, nicount = np.unique(independents, return_counts=True)
        #print(udofs, idofs)
        return constraints

    def run_sol_101_statics(self, subcase: Subcase,
                            f06_file: TextIO,
                            page_stamp: str, title: str='', subtitle: str='', label: str='',
                            page_num: int=1,
                            idtype: str='int32',
                            fdtype: str='float64') -> Any:
        """
        Runs a SOL 101

        Analysis (ASET): This set contains all boundary DOFs of the superelement.
                         It is considered fixed by default.
        Fixed Boundary (BSET): This subset of the ASET contains all fixed boundary DOFs.
        Free Boundary  (CSET): This subset of the ASET contains all free boundary DOFs.

        SOL 101 Sets
        ------------
        b = DOFs fixed during component mode analysis or dynamic reduction.
        c = DOFs that are free during component mode synthesis or dynamic reduction.
        lm = Lagrange multiplier DOFs created by the rigid elements
        r = Reference DOFs used to determine free body motion.
        l = b + c + lm
        t = l + r
        q = Generalized DOFs assigned to component modes and residual vectors
        a = t + q
        #-------------------------
        {s} = {sb} + {sg}  # done
        f =
        f = a + o        unconstrained (free) structural DOFs
        n = f + s        all DOFs not constrained by multipoint constraints
        ne = n + e       all DOFs not constrained by multipoint constraints plus extra DOFs
        m = mp + mr      all DOFs eliminated by multipoint constraints
        g = n + m        all DOFs including scalar DOFs

        Per Mystran
        -----------
        n = g - m        all DOFs not constrained by multipoint constraints
        f = n - s        unconstrained (free structural DOFs)
        a = f - o        ananlysis??? set
        L = a - r

        My Simplification
        -----------------
        g = n + (mp + mr)
        g = f + s + (mp + mr)
        g = f + (sb + sg) + (mp + mr)
        g = (a + o) + (sb + sg) + (mp + mr)
        a + o = g - (sb + sg) + (mp + mr)
        a = g - o - (sb + sg) + (mp + mr)

        l = b + c   # left over dofs
        t = l + r   # total set of DOFs for super
        a = t + q
        d = a + e   # dynamic analysis
        f = a + o
        n = f + s
        g = n + m
        p = g + e
        ps = p + sa
        pa = ps + k
        fe = f + e
        ne = n + e
        fr = f - q - r
        v = o + c + r

        New eqs
        f = n - s
        a = f - o
        t = a - q
        l = t - r
        b = l - c ???
        c = l - b ???
        """
        self.log.debug(f'run_sol_101')
        end_options = [
            'SEMG', # STIFFNESS AND MASS MATRIX GENERATION STEP
            'SEMR', # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            'SEKR', # STIFFNESS MATRIX REDUCTION STEP
            'SELG', # LOAD MATRIX GENERATION STEP
            'SELR', # LOAD MATRIX REDUCTION STEP
        ]
        #basic
        itime = 0
        ntimes = 1  # static
        isubcase = subcase.id
        #-----------------------------------------------------------------------
        model = self.model
        model.setup(run_geom_check=True)

        log = model.log

        dof_map, ps = _get_dof_map(model)

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)

        gset_b = ps_to_sg_set(ndof, ps)
        Kgg = build_Kgg(model, dof_map,
                         ndof, ngrid, ndof_per_grid,
                         idtype='int32', fdtype=fdtype)
        Mbb = build_Mbb(model, subcase, dof_map, ndof, fdtype=fdtype)
        #print(self.op2.grid_point_weight)
        reference_point, MO = grid_point_weight(model, Mbb, dof_map, ndof)
        weight = make_grid_point_weight(
            reference_point, MO,
            approach_code=1, table_code=13,
            title=title, subtitle=subtitle, label=label,
            superelement_adaptivity_index='')
        self.op2.grid_point_weight[label] = weight
        page_num = weight.write_f06(f06_file, page_stamp, page_num)

        self.get_mpc_constraints(subcase, dof_map)

        Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid, inplace=False)
        del Mbb

        gset = np.arange(ndof, dtype=idtype)
        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)
        #sset_b = np.union1d(sset_b, sset_g)
        #print(sset_g)
        #print(sset_b)
        #print(sset)
        Fb = self.build_Fb(xg, sset_b, dof_map, ndof, subcase)
        Fg = Fb
        # ----------------------MPC reduction----------------------
        if 0:
            mset_b = np.zeros(ndof, dtype='bool')
            nelement_mpc = len(model.mpc) + len(model.rbe1) + len(model.rbe2) + len(model.rbe3)
            if nelement_mpc == 0:
                Knn = Kgg
                Fn = Fg
                xn = xg
            else:  # pragma: no cover
                raise RuntimeError()

            nset_b = gset_b & ~mset_b  # n = g-m
            del Kgg, Mgg, Fg
            # ----------------------SPC reduction NEW----------------------

            fset_b = nset_b & sset_b
            fset = gset[fset_b]
            K = partition_matrix(Knn, [('f', fset), ('s', sset), ('0', set0)])
            Ff, Fs = partition_vector2(Fn, [['f', fset], ['s', sset]])
            xf, xs, x0 = partition_vector3(xn, [['f', fset], ['s', sset], ['0', set0]])

            Kff = K['ff']
            Kfs = K['fs']
            Kss = K['ss']

            # ----------------------OMIT reduction----------------------
            #a = f - o
            omit_b = np.zeros(ndof, dtype='bool')
            if len(model.omit) == 0:
                Kaa = Kff
                Fa = Ff
                xa = xf
            else:  # pragma: no cover
                raise RuntimeError()

            aset_b = fset_b & omit_b
            aset = gset[aset_b]

            # q=0 for statics/modes, so a=t
            Ktt = Kaa
            xt = xa
            Ft = Fa
            tset_b = aset_b
            tset = aset

            # ----------------------OMIT reduction----------------------
            # L = t - r
            # rset
            rset_b = np.zeros(ndof, dtype='bool')
            if len(model.suport) == 0:
                Kll = Ktt
            else:  # pragma: no cover
                raise RuntimeError()
            lset_b = tset_b & rset_b
            lset = gset[lset_b]

            self.mset = gset[mset_b]
            self.rset = gset[rset_b]
            self.aset = aset
            self.lset = lset

        # ----------------------SPC reduction OLD----------------------
        # Constrained set

        # sset = sb_set | sg_set

        # free structural DOFs
        #fset = ~sset; # & ~mset
        fset_b = gset_b & sset_b # g-b
        #self.log.debug(f'gset_b = {gset_b}')
        #self.log.debug(f'sset_b = {sset_b}')
        #self.log.debug(f'fset_b = {fset_b}')
        aset, tset = get_residual_structure(model, dof_map, fset_b)

        asetmap = get_aset(model)
        if asetmap:
            aset = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype)
        else:
            aset = np.setdiff1d(gset, sset) # a = g-s
            #fset = np.setdiff1d(nset, sset) # f = n-s
        naset = aset.sum()
        nsset = sset.sum()
        if naset == 0 and nsset == 0:
            raise RuntimeError('no residual structure found')
        # The a-set and o-set are created in the following ways:
        # 1. If only OMITi entries are present, then the o-set consists
        #    of DOFs listed explicitly on OMITi entries. The remaining
        #    f-set DOFs are placed in the b-set, which is a subset of
        #    the a-set.
        # 2. If ASETi or QSETi entries are present, then the a-set consists
        #    of all DOFs listed on ASETi entries and any entries listing its
        #    subsets, such as QSETi, SUPORTi, CSETi, and BSETi entries.  Any
        #    OMITi entries are redundant. The remaining f-set DOFs are placed
        #    in the o-set.
        # 3. If there are no ASETi, QSETi, or OMITi entries present but there
        #    are SUPORTi, BSETi, or CSETi entries present, then the entire
        #    f-set is placed in the a-set and the o-set is not created.

        page_num = write_oload(Fg, dof_map, isubcase, ngrid, ndof_per_grid,
                               f06_file, page_stamp, page_num, log)

        # aset - analysis set
        # sset - SPC set
        #print('aset = ', aset)
        #print('sset = ', sset)

        # u1 = Kaa^-1 * (F1k - Kas*u2k)
        abs_xg = np.abs(xg)
        finite_xg = np.any(np.isfinite(abs_xg))
        if finite_xg and np.nanmax(abs_xg) > 0.:
            self.log.info(f'SPCD found')
            self.log.info(f'  xg = {xg}')
            self.log.info(f'  Fg = {Fg}')
            set0 = (xg == 0.)
            set0_ = np.where(set0)

            self.log.info(f'  aset = {aset}')
            if np.any(set0_):
                self.log.info(f'  removing set0_ from sset')
                self.log.info(f'    set0_ = {set0_}')
                self.log.info(f'    sset = {sset}')
                sset = np.setdiff1d(sset, set0_)
                self.log.info(f'    -> sset = {sset}')
            else:
                self.log.info(f'  sset = {sset}')

            #self.log.info(f'  Fa_solve = {Fa_solve}')
            #Fa = Fa.reshape(len(Fa), 1)
            #xa = Fa.reshape(len(xa), 1)
            #print(f'Fa.shape = {Fa.shape}')
            #print(f'xa.shape = {xa.shape}')
            #print(f'Kas.shape = {Kas.shape}')
            #Kas_xa = Ksa @ xa
            #print(f'Kas @ xa.shape = {Kas_xa.shape}')
            #Fa_solve = Fa - Kas @ xa
            #Kspc = Kas[ixa, :][:, ixa]
            #print(Kspc.shape)
        else:
            set0 = sset
            sset = []
            #self.sset = sset
        self.sset = sset
        self.set0 = set0
        self.aset = aset

        Fa, Fs = partition_vector2(Fb, [['a', aset], ['s', sset]])
        del Fb, Fs

        xa, xs, x0 = partition_vector3(xg, [['a', aset], ['s', sset], ['0', set0]])
        #self.log.info(f'xg = {xg}')
        del xg
        #self.log.info(f'xa = {xa}')
        #self.log.info(f'xs = {xs}')
        #self.log.info(f'x0 = {x0}')
        del x0

        #print(Kgg)
        self.Kgg: np.ndarray = Kgg.toarray()
        K = partition_matrix(Kgg, [('a', aset), ('s', sset), ('0', set0)])
        Kaa = K['aa']
        Kss = K['ss']
        #Kas = K['as']
        Ksa = K['sa']
        #Ks = partition_matrix(Kggs, [['a', aset], ['s', sset], ['0', set0]])
        #Kaas = Ks['aa']
        #assert Kaa.shape == Kaas.shape

        self.Kaa = Kaa

        #M = partition_matrix(Mgg, [['a', aset], ['s', sset]])
        #Maa = M['aa']
        #Kss = K['ss']
        Kas = K['as']
        Ksa = K['sa']
        K0a = K['0a']
        K0s = K['0s']
        #self.Maa = Maa
        #[Kaa]{xa} + [Kas]{xs} = {Fa}
        #[Ksa]{xa} + [Kss]{xs} = {Fs}

        #{xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
        #{Fs} = [Ksa]{xa} + [Kss]{xs}

        #print(Kaa)
        #print(Kas)
        #print(Kss)
        Fa_solve: np.ndarray = Fa
        is_sset = len(sset)
        is_set0 = len(set0)
        is_aset = len(aset)
        if is_sset:
            Fa_solve = Fa - Kas @ xs
            self.log.info(f'  Fa_solve = {Fa_solve}')

        Fg_oload = Fg.copy()
        if is_aset:
            xa_, ipositive, inegative = solve(Kaa, Fa_solve, aset, log, idtype=idtype)
            Fa_ = Fa[ipositive]

            log.info(f'aset_ = {ipositive}')
            log.info(f'xa_ = {xa_}')
            log.info(f'Fa_ = {Fa_}')

            xa[ipositive] = xa_
            xa[inegative] = 0.
            self.xa_ = xa_
            self.Fa_ = Fa_
        else:
            self.log.warning('A-set is empty; all DOFs are constrained')
            self.xa_ = []
            self.Fa_ = []


        xg = np.full(ndof, np.nan, dtype=fdtype)
        #print('aset =', aset)
        #print('sset =', sset)
        #print('set0 =', set0)
        xg[aset] = xa
        xg[sset] = xs
        xg[set0] = 0.
        Fg[aset] = Fa
        #print(xg)
        self.xg = xg
        self.Fg = Fg

        fspc = np.full(ndof, 0., dtype=fdtype)

        if is_sset:
            log.info(f'fspc_s recovery')
            log.debug(f'  aset = {aset}')
            log.debug(f'  sset = {sset}')
            log.debug(f'  xa = {xa}')
            log.debug(f'  xs = {xs}')
            fspc_s = Ksa @ xa + Kss @ xs
            log.debug(f'  fspc_s = {fspc_s}')
            log.debug(f'  Ksa @ xa = {Ksa @ xa}')
            log.debug(f'  Kss @ xs = {Kss @ xs}')
            log.debug(f'  sum = {Ksa @ xa + Kss @ xs}')
            Fg[sset] = fspc_s
            fspc[sset] = fspc_s
        if is_set0:
            log.info(f'fspc_0 recovery')
            fspc_0 = K0a @ xa + K0s @ xs
            log.debug(f'  fspc_0 = {fspc_0}')
            Fg[set0] = fspc_0
            fspc[set0] = fspc_0

        log.debug(f'xa = {xa}')
        log.debug(f'xs = {xs}')
        log.debug(f'fspc = {fspc}')

        xb = xg_to_xb(model, xg, ngrid, ndof_per_grid)
        Fb = xg_to_xb(model, Fg, ngrid, ndof_per_grid)

        #log.debug(f'Fs = {Fs}')
        #log.debug(f'Fb = {Fb}')
        #log.debug(f'xb = {xb}')

        self._save_displacment(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, xg,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)

        self._save_applied_load(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, Fg_oload,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)

        self._save_spc_forces(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, fspc,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)

        op2 = self.op2
        page_stamp += '\n'
        if 'FORCE' in subcase:
            recover_force_101(f06_file, op2, self.model, dof_map, subcase, xb,
                              title=title, subtitle=subtitle, label=label,
                              page_stamp=page_stamp)

        #if 'STRAIN' in subcase:
            #recover_strain_101(f06_file, op2, self.model, dof_map, subcase, xb,
                               #title=title, subtitle=subtitle, label=label,
                               #page_stamp=page_stamp)
        #if 'STRESS' in subcase:
            #recover_stress_101(f06_file, op2, self.model, dof_map, subcase, xb,
                               #title=title, subtitle=subtitle, label=label,
                               #page_stamp=page_stamp)
        if 'ESE' in subcase:
            recover_strain_energy_101(f06_file, op2, self.model, dof_map, subcase, xb,
                                      title=title, subtitle=subtitle, label=label,
                                      page_stamp=page_stamp)
        #Fg[sz_set] = -1
        #xg[sz_set] = -1
        op2.write_op2(self.op2_filename, post=-1, endian=b'<', skips=None, nastran_format='nx')
        self.log.info('finished')
        return end_options

    def _save_displacment(self, f06_file: TextIO,
                          subcase: Subcase, itime: int, ntimes: int,
                          node_gridtype: NDArrayN2int,
                          xg: NDArrayNfloat,
                          ngrid: int, ndof_per_grid: int,
                          title: str='', subtitle: str='', label: str='',
                          fdtype: str='float32',
                          page_num: int=1, page_stamp: str='PAGE %s') -> int:
        f06_request_name = 'DISPLACEMENT'
        table_name = 'OUGV1'
        #self.log.debug(f'xg = {xg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, xg,
            RealDisplacementArray, f06_request_name, table_name, self.op2.displacements,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

    def _save_spc_forces(self, f06_file: TextIO,
                         subcase: Subcase, itime: int, ntimes: int,
                         node_gridtype: NDArrayN2int, fspc: NDArrayNfloat,
                         ngrid: int, ndof_per_grid: int,
                         title: str='', subtitle: str='', label: str='',
                         fdtype: str='float32',
                         page_num: int=1, page_stamp: str='PAGE %s') -> int:
        f06_request_name = 'SPCFORCES'
        table_name = 'OQG1'
        #self.log.debug(f'Fg = {Fg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, fspc,
            RealSPCForcesArray, f06_request_name, table_name, self.op2.spc_forces,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

    def _save_applied_load(self, f06_file: TextIO,
                           subcase: Subcase, itime: int, ntimes: int,
                           node_gridtype: NDArrayN2int,
                           Fg: NDArrayNfloat,
                           ngrid: int, ndof_per_grid: int,
                           title: str='', subtitle: str='', label: str='',
                           fdtype: str='float32',
                           page_num: int=1, page_stamp: str='PAGE %s') -> int:
        f06_request_name = 'OLOAD'
        table_name = 'OPG1'
        #self.log.debug(f'Fg = {Fg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, Fg,
            RealLoadVectorArray, f06_request_name, table_name, self.op2.load_vectors,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

    def _recast_data(self, idtype: str, fdtype: str):
        idtype = 'int32'
        fdtype = 'float32'
        return idtype, fdtype

    def _save_static_table(self, f06_file: TextIO,
                           subcase: Subcase, itime: int, ntimes: int,
                           node_gridtype: NDArrayN2int, Fg: NDArrayNfloat,
                           obj: Union[RealSPCForcesArray],
                           f06_request_name: str,
                           table_name: str, slot: dict[Any, RealSPCForcesArray],
                           ngrid: int, ndof_per_grid: int,
                           title: str='', subtitle: str='', label: str='',
                           idtype: str='int32', fdtype: str='float32',
                           page_num: int=1, page_stamp: str='PAGE %s') -> int:
        idtype, fdtype = self._recast_data(idtype, fdtype)
        isubcase = subcase.id
        #self.log.debug(f'saving {f06_request_name} -> {table_name}')
        unused_nids_write, write_f06, write_op2, quick_return = get_plot_request(
            subcase, f06_request_name)
        if quick_return:
            return page_num
        nnodes = node_gridtype.shape[0]
        data = np.zeros((ntimes, nnodes, 6), dtype=fdtype)
        ngrid_dofs = ngrid * ndof_per_grid
        if ndof_per_grid == 6:
            _fgi = Fg[:ngrid_dofs].reshape(ngrid, ndof_per_grid)
            data[itime, :ngrid, :] = _fgi
        else:
            raise NotImplementedError(ndof_per_grid)
        data[itime, ngrid:, 0] = Fg[ngrid_dofs:]

        spc_forces = obj.add_static_case(
            table_name, node_gridtype, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
        if write_f06:
            page_num = spc_forces.write_f06(
                f06_file, header=None,
                page_stamp=page_stamp, page_num=page_num,
                is_mag_phase=False, is_sort1=True)
            f06_file.write('\n')
        if write_op2:
            slot[isubcase] = spc_forces
        return page_num

    def run_sol_103_modes(self, subcase: Subcase,
                          f06_file: TextIO,
                          page_stamp: str,
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1,
                          idtype: str='int32', fdtype: str='float64'):
        """
        [M]{xdd} + [C]{xd} + [K]{x} = {F}
        [M]{xdd} + [K]{x} = {F}
        -[M]{xdd}λ^2 + [K]{x} = {0}
        {X}(λ^2 - [M]^-1[K]) = {0}
        λ^2 - [M]^-1[K] = {0}
        λ^2 = [M]^-1[K]
        [A][X] = [X]λ^2
        """
        self.log.debug(f'run_sol_103 (modes)')
        model = self.model
        assert len(model.methods), 'SOL 103 (modes) requires a METHOD and a EIGR/EIGRL card'
        end_options = [
            'SEMR',  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            'SEKR',  # STIFFNESS MATRIX REDUCTION STEP
        ]
        op2 = self.op2
        #log = model.log
        #write_f06 = True

        node_gridtype = _get_node_gridtype(model, idtype=idtype)

        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)
        out = _run_modes(model, subcase,
                         ngrid, ndof_per_grid, ndof,
                         node_gridtype, dof_map,
                         idtype=idtype, fdtype=fdtype)
        xg_out = out['xg_out']
        eigenvalues = out['eigenvalues']
        nmodes = len(eigenvalues)

        isubcase = subcase.id
        mode_cycles = eigenvalues
        eigenvalues_obj = RealEigenvalues(title, 'LAMA', nmodes=nmodes)

        cycle = np.sqrt(np.abs(eigenvalues)) / (2. * np.pi)
        radians = np.sqrt(np.abs(eigenvalues))

        eigenvalues_obj.mode = np.arange(nmodes, dtype='int32') + 1
        eigenvalues_obj.extraction_order = np.arange(nmodes, dtype='int32') + 1
        eigenvalues_obj.eigenvalues = eigenvalues
        eigenvalues_obj.radians = radians
        eigenvalues_obj.cycles = cycle
        eigenvalues_obj.generalized_mass = np.ones(nmodes, dtype='float64')
        eigenvalues_obj.generalized_stiffness = np.ones(nmodes, dtype='float64')

        op2.eigenvalues[title] = eigenvalues_obj
        #op2.eigenvalues[isubcase] = eigenvalues_obj

        nnodes = node_gridtype.shape[0]
        data = xg_out.reshape((nmodes, nnodes, 6)).astype('float32')
        eigenvalues = eigenvalues.astype('float32')
        table_name = 'OUGV1'
        modes = np.arange(1, nmodes + 1, dtype=idtype)
        eigenvectors_obj = RealEigenvectorArray.add_modal_case(
            table_name, node_gridtype, data, isubcase, modes, eigenvalues, mode_cycles,
            is_sort1=True, is_random=False, is_msc=True, random_code=0,
            title=title, subtitle=subtitle, label=label)
        eigenvectors_obj.nonlinear_factor = 1
        op2.eigenvectors[isubcase] = eigenvectors_obj

        write_f06 = True
        if write_f06:
            str(page_num)
            header = []
            page_num = eigenvalues_obj.write_f06(
                f06_file, header=header,
                page_stamp=page_stamp, page_num=page_num)
            f06_file.write('\n')

            page_num = eigenvectors_obj.write_f06(
                f06_file, header=None,
                page_stamp=page_stamp, page_num=page_num,
                is_mag_phase=False, is_sort1=True)
            f06_file.write('\n')
        #fspc = Ksa @ xa + Kss @ xs
        #Fs[ipositive] = Fsi

        #Fg[aset] = Fa
        #Fg[sset] = fspc
        #log.debug(f'xa = {xa}')
        #log.debug(f'Fs = {Fs}')
        #log.debug(f'xg = {xg}')
        #log.debug(f'Fg = {Fg}')

        str(f06_file)
        str(page_stamp)
        op2.write_op2(self.op2_filename, post=-1, endian=b'<', skips=None, nastran_format='nx')
        return end_options
        #raise NotImplementedError(subcase)

    def run_sol_111_harmonic(self, subcase: Subcase, f06_file: TextIO,
                             page_stamp: str, title: str='', subtitle: str='', label: str='',
                             page_num: int=1,
                             idtype: str='int32', fdtype: str='float64'):
        """
        frequency
        [M]{xdd} + [C]{xd} + [K]{x} = {F}
        {φ} ([M]s^2 + [C]{s} + [K]) = {F}
        {φ} ([M]s^2 + [C]{s} + [K]) = {F}
        {φ} ([M]s^2 + [K]) = {F}

        {d}*sin(wt) ([M]s^2 + [K]) = {0}
        [M]{xdd} + [K]{x} = {0}
        x = A*sin(wt)
        xd = dx/dt = A*w*cos(wt)
        xdd = -A*w^2 * sin(wt)
        -w^2[M]{A}*sin(wt) + [K]{A}*sin(wt) = {0}
        (-w^2[M] + [K]){A} = {0}
        -w^2[M] + [K] = 0
        w^2 = [M]^-1 * [K]

        -w^2[M]{A}*sin(wt) + [K]{A}*sin(wt) = {F}sin(wt)
        (-w^2[M] + [K]){A} = {F}
        {A} = (-w^2[M] + [K])^-1 * {F}
        {A} = ([K] - w^2[M])^-1 * {F}

        {φ}^T[M]{φ}{qdd} + {φ}^T[B]{φ}{qd} {φ}^T[K]{φ}{q} = {φ}^T{P(t)}
        {qdd} + 2*zeta*w_n{qd} w_n^2{q} = 1/M_n {φ}^T{P(t)}
        {u_total} = {u_dynamic} + {u_static} = [PHI]{q} + {u_static}

        """
        str(subcase)
        str(f06_file)
        str(fdtype)
        str(idtype)
        # https://pyfrf.readthedocs.io/en/latest/Showcase.html
        #ndof = K.shape[0]

        #end_options = [
            #'SEMR',  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            #'SEKR',  # STIFFNESS MATRIX REDUCTION STEP
        #]
        model = self.model
        #op2 = self.op2
        #---------------------------------------------------
        #title = ''
        #today = None
        #page_stamp = self.op2.make_stamp(title, today) # + '\n'

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)

        out = _run_modes(
            model, subcase,
            ngrid, ndof_per_grid, ndof,
            node_gridtype, dof_map,
            idtype=idtype, fdtype=fdtype)

        aset = out['aset']
        Kaa = out['Kaa']
        Maa = out['Maa']
        eigenvalues = out['eigenvalues']
        omegas = np.sqrt(np.abs(eigenvalues))

        freqs = get_frequencies(model, subcase, omegas)
        nfreq = len(freqs)
        omegai = 2 * np.pi * freqs

        ndof_a = Kaa.shape[0]
        Caa = np.eye(ndof_a)
        np.fill_diagonal(Caa, 0.01)

        #nfreq = len(omegas)
        shape = (ndof, ndof, nfreq)
        tf_matrix = np.zeros(shape, dtype='complex128')  # full system 3x3 TF matrix
        tf_mask = np.ones((ndof, ndof, nfreq), dtype='bool')
        tf_mask[~aset, :, :] = 0
        tf_mask[:, ~aset, :] = 0
        tf_matrix_a = np.ones((ndof_a, ndof_a, nfreq), dtype='complex128')
        #tf_matrix_a = np.ma.masked_array(tf_matrix, tf_mask)
        for i, omegai in enumerate(omegas):
            mat = -omegai**2*Maa + 1j*omegai*Caa + Kaa
            res = sp.linalg.inv(mat)
            tf_matrix_a[:,:,i] = res

        tf_matrix = np.zeros(shape, dtype='complex128')  # full system 3x3 TF matrix
        tf_matrix.ravel()[tf_mask.ravel()] = tf_matrix_a.ravel()

        # TODO: save tf_matrix to ACCEL response
        return tf_matrix


def partition_matrix(matrix, sets: list[tuple[str, np.ndarray]]) -> dict[tuple[str, str], NDArrayNNfloat]:
    """partitions a matrix"""
    matrices = {}
    for aname, aset in sets:
        for bname, bset in sets:
            matrices[aname + bname] = matrix[aset, :][:, bset]
    return matrices

def partition_vector(vector, sets, fdtype: str='float64') -> list[NDArrayNfloat]:  # pragma: no cover
    """partitions a vector"""
    vectors = []
    for unused_aname, aset in sets:
        if len(aset) == 0:
            vectors.append(np.array([], dtype=fdtype))
            continue
        vectori = vector[aset]
        vectors.append(vectori)
    return vectors

def partition_vector2(vector, sets,
                      fdtype: str='float64') -> tuple[NDArrayNfloat, NDArrayNfloat]:
    """partitions a vector"""
    assert len(sets) == 2, sets
    #vectors = partition_vector(vector, sets, fdtype=fdtype)
    #return tuple(vectors)
    (unused_name0, set0), (unused_name1, set1) = sets
    vectors = (vector[set0], vector[set1])
    return vectors

def partition_vector3(vector, sets,
                      fdtype: str='float64') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """partitions a vector"""
    assert len(sets) == 3, sets
    #vectors = partition_vector(vector, sets, fdtype=fdtype)
    (unused_name0, set0), (unused_name1, set1), (unused_name2, set2) = sets
    vectors = (vector[set0], vector[set1], vector[set2])
    return tuple(vectors)

def remove_rows(Kgg: NDArrayNNfloat, aset: NDArrayNint,
                idtype: str='int32') -> NDArrayNNfloat:
    """
    Applies AUTOSPC to the model (sz)

    mp  DOFs eliminated by multipoint constraints.
    mr  DOFs eliminated by multipoint constraints created by the rigid
        elements using the LGELIM method on the Case Control command RIGID.
    sb* DOFs eliminated by single-point constraints that are included
        in boundary condition changes and by the AUTOSPC feature.
        (See the sz set)
    sg* DOFs eliminated by single-point constraints that are specified
        on the PS field on GRID Bulk Data entries.
    sz  DOFs eliminated by the AUTOSPC feature.
    o   DOFs omitted by structural matrix partitioning.
    q   Generalized DOFs assigned to component modes and residual vectors.
    r   Reference DOFs used to determine free body motion.
    c   DOFs that are free during component mode synthesis or dynamic reduction.
    b   DOFs fixed during component mode analysis or dynamic reduction.
    lm  Lagrange multiplier DOFs created by the rigid elements
        using the LAGR method on the Case Control command, RIGID.
    e   Extra DOFs introduced in dynamic analysis.
    sa  Permanently constrained aerodynamic DOFs.
    k   Aerodynamic mesh point set for forces and displacements on the aero mesh.
    j   Aerodynamic mesh collocation point set (exact physical
        interpretation is dependent on the aerodynamic theory).


    Input Sets
    ----------
    mp  DOFs eliminated by multipoint constraints.
    mr  DOFs eliminated by multipoint constraints created by the rigid
        elements using the LGELIM method on the Case Control command RIGID.
    g   all DOFs including scalar DOFs
    o   Omit set
    s   SPC sets (sb + sg)

    Source of Sets
    --------------
    g     6-DOFs per GRID, 1-DOF SPOINTs,
          Nmodes-DOF EPOINTs (Extra Dynamic/Modal Point; ),
          nRBE-DOF LPOINTs ("Lagrange" Point; e.g., RBE3),
          nHarmonic-DOF HPOINTs (Harmonic Point; e.g., RINGFL)
    m     MPC, MPCADD, RBAR/RBAR1, RBE1/RBE2/RBE3, RROD, RSPLINE,
          RJOINT, RSSCON, RTRPLT/RTRPLT1, GMBC, GMSPC*
    r     SUPORT/SUPORT1, SUPAX
    s     SPC/SPC1, SPCADD, BNDGRID, AUTOSPC
    sb    SPC/SPC1, SPCADD, SPCAX, FLSYM, GMSPC*, BNDGRID, PARAM,AUTOSPC,YES
    sz    PARAM,AUTOSPC,YES
    sg    GRID, GRIDB, GRDSET (PS field)
    o     OMIT/OMIT1, GRID (SEID field), SESET
    q     QSET/QSET1
    sa    CAEROi
    a     ASET/ASET1, CSUPEXT, superelement exterior DOFs
    c     CSET/CSET1, BNDREE/BNDFRE1
    b     BSET/BSET1, BNDFIX/BNDFIX1
    ap    ACCSPT
    rv    RVDOF/RVDOF1
    u1-u6 USET/USET1, SEUSET/SEUSET1

    * placed in the "m" set only if constraints are not specified in
      the basic coordinate system


    Set Name         Description
    --------         -----------
    s = sb + sg      all DOFs eliminated by single point constraints
    l = b + c + lm   the DOFs remaining after the reference DOFs are removed (DOF left over)
    t = l + r        the total set of physical boundary DOF for superelements
    a = t + q        the analysis set used in eigensolution
    d = a + e        the set used in dynamic analysis by the direct method
    f = a + o        unconstrained (free) structural DOFs
    fe = f + e       free DOFs plus extra DOFs
    n = f + s        all DOFs not constrained by multipoint constraints
    ne = n + e       all DOFs not constrained by multipoint constraints plus extra DOFs
    m = mp + mr      all DOFs eliminated by multipoint constraints
    g = n + m        all DOFs including scalar DOFs
    p = g + e        all physical DOFs including extra point DOFs
    ks = k + sa      the union of k and the re-used s-set (6 dof per grid)
    js = j + sa      the union of j and the re-used s-set (6 dof per grid)
    fr = o + l       statically independent set minus the statically determinate supports (fr = f – q – r)
    v = o + c + r    the set free to vibrate in dynamic reduction and component mode synthesis
    al = a – lm      a-set  without Lagrange multiplier DOFs
    dl = d – lm      d-set  without Lagrange multiplier DOFs
    gl = g – lm      g-set  without Lagrange multiplier DOFs
    ll = l – lm      l-set  without Lagrange multiplier DOFs
    nf = ne – lm     ne-set without Lagrange multiplier DOFs
    pl = p – lm      p-set  without Lagrange multiplier DOFs
    tl = t – lm      t-set  without Lagrange multiplier DOFs
    nl = n – lm      n-set  without Lagrange multiplier DOFs
    fl = f – lm      f-set  without Lagrange multiplier DOFs
    ff = fe – lm     fe-set without Lagrange multiplier DOFs


    Per MSC Nastran
    ---------------
    js = j + sa  the union of j and the re-used s-set (6 dof per grid)
    ks = k + sa  the union of k and the re-used s-set (6 dof per grid)

    [K]{x} = {F}
    a - active
    s - SPC
    [Kaa Kas]{xa} = {Fa}
    [Ksa Kss]{xs}   {Fs}
    """
    abs_kgg = np.abs(Kgg)
    col_kgg = abs_kgg.max(axis=0)
    row_kgg = abs_kgg.max(axis=1)
    #print(abs_kgg)
    #print(col_kgg)
    #print(row_kgg)
    #izero = np.where((col_kgg == 0.) & (row_kgg == 0))[0]
    if isinstance(Kgg, np.ndarray):
        ipositive = np.where((col_kgg > 0.) | (row_kgg > 0))[0]
    elif isinstance(Kgg, (csc_matrix, lil_matrix)):
        ipositive1 = col_kgg.todense().nonzero()[1]
        ipositive2 = row_kgg.todense().T.nonzero()[1]
        ipositive = np.union1d(ipositive1, ipositive2)
        if isinstance(Kgg, lil_matrix):
            Kgg = Kgg.tocsc()
    else:
        raise NotImplementedError(type(Kgg))
    all_rows = np.arange(Kgg.shape[0], dtype=idtype)
    inegative = np.setdiff1d(all_rows, ipositive)
    apositive = aset[ipositive]
    sz_set = np.setdiff1d(aset, apositive)
    Kaa = Kgg[ipositive, :][:, ipositive]
    return Kaa, ipositive, inegative, sz_set

def guyan_reduction(matrix, set1, set2):
    """
    https://en.wikipedia.org/wiki/Guyan_reduction

    [A11 A12]{d1} = {F1}
    [A21 A22]{d2} = {0}
    [A11]{d1} + [A12]{d2} = {F}
    [A21]{d1} + [A22]{d2} = {0}
    {d2} = -[A22]^-1[A21]{d1}
    [A11]{d1} + [A12](-[A22]^-1[A21]{d1}) = {F}
    ([A11] + ([A12]-[A22]^-1[A21]){d1} = {F}

    (AB)^-1 = B^-1 A^-1
    (cA)^-1 = A^-1/c
    (A)^-n = (A^-1)^n
    (A^-1)^T = (A^T)^-1

    https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
    (A + B)^-1 = A^-1 - 1/(1+g)A^-1 B A^-1; g=trace(B A^-1) where g != -1
               = A^-1 - A^-1*B*(A+B^-1)
    (A + B)^-1 = A^-1 + X
    (A^-1 + X)(A + B) = I
    A^-1 A + X A + A^-1 B + X B = I
    X(A + B) = A^-1 B
    X = -A^-1 B (A^-1 + B)^-1
    X = -A^-1 B (A^-1 + X)
    (I + A^-1 B) X = A^-1 B A^-1
    X = -(I + A^-1 B)^-1 (A^-1 B A^-1)
    """
    sets = [
        ['1', set1],
        ['2', set2],
    ]
    A = partition_matrix(matrix, sets)
    A11 = A['11']
    A22 = A['22']
    #A12 = A['12']
    A21 = A['21']
    nA11 = A11.shape[0]

    # ([A11] + ([A12]-[A22]^-1[A21]){d1} = {F}
    A22m1 = np.linalg.inv(A22)
    #A12_A22m1_A21 = np.linalg.multi_dot([A11, A22m1, A21])
    T = np.vstack([np.eye(nA11), A22m1 @ A21])
    return np.linalg.multi_dot([T.T, A, T])
    #return A11 + A12_A22m1_A21

def _get_node_gridtype(model: BDF, idtype: str='int32') -> NDArrayN2int:
    """
    Helper method for results post-processing

    Point type (per NX 10; OUG table; p.5-663):
    =1, GRID Point
    =2, Scalar Point
    =3, Extra Point
    =4, Modal
    =5, p-elements, 0-DOF
    -6, p-elements, number of DOF

    """
    node_gridtype = []
    #fdtype = 'float64'
    if model.grid.n:
        nodes = model.grid.node_id
        ones = np.ones(model.grid.n, dtype=idtype)
        node_gridtype.append(np.column_stack([nodes, ones]))

    if model.spoint.n:
        #raise NotImplementedError()
        spoint_id = model.spoint.spoint_id
        spoint_id.sort()
        twos = np.ones(model.spoint.n, dtype=idtype) * 2
        node_gridtype.append(np.column_stack([spoint_id, twos]))

    assert len(node_gridtype) > 0
    #print(node_gridtype)
    node_gridtype_array = np.vstack(node_gridtype)
    return node_gridtype_array

def get_aset(model: BDF) -> set[tuple[int, int]]:
    return model.aset.set_map

def get_bset(model: BDF) -> set[tuple[int, int]]:
    """creates the b-set"""
    return model.bset.set_map

def get_cset(model: BDF) -> set[tuple[int, int]]:
    """creates the c-set"""
    return model.cset.set_map

def get_omit_set(model: BDF) -> set[tuple[int, int]]:
    """creates the o-set"""
    return model.omit.set_map

def get_rset(model: BDF) -> set[tuple[int, int]]:
    """creates the r-set"""
    rset_map = set()
    return rset_map
    for rset in model.suport:
        for nid, comp in zip(rset.ids, rset.components):
            for compi in comp:
                rset_map.add((nid, int(compi)))

    for suport in model.suport1:
        comp = suport.components
        for nid in suport.ids:
            for compi in comp:
                rset_map.add((nid, int(compi)))
    return rset_map

def get_qset(model: BDF) -> set[tuple[int, int]]:
    """creates the q-set"""
    return model.qset.set_map

def get_residual_structure(model: BDF,
                           dof_map: DOF_MAP,
                           fset: NDArrayNbool, idtype: str='int32') -> NDArrayNbool:
    """gets the residual structure dofs"""
    asetmap = get_aset(model)
    bsetmap = get_bset(model)
    csetmap = get_cset(model)
    rsetmap = get_rset(model)
    qsetmap = get_qset(model)
    osetmap = get_omit_set(model)
    aset = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype, use_ints=False)
    bset = apply_dof_map_to_set(bsetmap, dof_map, idtype=idtype, use_ints=False)
    cset = apply_dof_map_to_set(csetmap, dof_map, idtype=idtype, use_ints=False)
    rset = apply_dof_map_to_set(rsetmap, dof_map, idtype=idtype, use_ints=False)
    oset = apply_dof_map_to_set(osetmap, dof_map, idtype=idtype, use_ints=False)
    qset = apply_dof_map_to_set(qsetmap, dof_map, idtype=idtype, use_ints=False)

    #ndof = len(dof_map)
    #print('aset =', aset)
    #print('qset =', qset)
    #print('fset =', fset)
    #print('oset =', oset)

    aqo_set = np.vstack([aset, qset, oset])
    bcr_set = np.vstack([bset, cset, rset])
    # The a-set and o-set are created in the following ways:
    #    1. If only OMITi entries are present, then the o-set consists
    #       of degrees-of-freedom listed explicitly on OMITi entries.
    #       The remaining f-set degrees-of-freedom are placed in the
    #       b-set, which is a subset of the a-set.
    if np.any(oset) and ~np.any([aset, bset, cset, qset, rset]):
        # b = f - o
        bset = np.setdiff1d(fset, oset)
        aset = bset
        # 2. If ASETi or QSETi entries are present, then the a-set
        #    consists of all degrees-of-freedom listed on ASETi entries
        #    and any entries listing its subsets, such as QSETi, SUPORTi
        #    CSETi, and BSETi entries. Any OMITi entries are redundant.
        #    The remaining f-set degrees-of-freedom are placed in the
        #    o-set.
    elif np.any([aset, qset]):
        abcqr_set = np.vstack([aset, bset, cset, qset, rset])
        #print(abcqr_set)
        aset = np.sum(abcqr_set, axis=0).astype('bool') # dim=2 -> row
        if np.any(oset): # check no overlap with O set
            ao_set = aset | oset
            if np.any(ao_set):
                raise RuntimeError('OMITi entries cannot overlap with ASETi entries '
                                   'or any ASET subsets, such as QSETi, '
                                   'SUPORTi, CSETi, and BSETi entries.')

        oset = fset & ~aset # assign remaining to O set
        # 3. If there are no ASETi, QSETi, or OMITi entries present but
        #    there are SUPORTi, BSETi, or CSETi entries present, then
        #    the entire f-set is placed in the a-set and the o-set is
        #    not created.
    elif not np.any(aqo_set) and np.any(bcr_set):
        aset = fset
        # 4. There must be at least one explicit ASETi, QSETi, or OMITi
        #    entry for the o-set to exist, even if the ASETi, QSETi, or
        #    OMITi entry is redundant. (related to item 3)
    else:
        # No model reduction - same as previous option
        aset = fset
    # Add ASET to residual structure TSET
    tset = aset & ~qset

    # Exclusive Degrees-of-freedom sets
    # ---------------------------------
    # m  # ([nGdof,1] logical) Degrees-of-freedom eliminated by multiple constraints
    # sb # ([nGdof,numSID] logical) Degrees-of-freedom eliminated by single-point constraints that are included in boundary conditions
    # sg # ([nGdof,1] logical) Degrees-of-freedom eliminated by single-point constraints that are specified on the PS field on node entries
    # o  # ([nGdof,1] logical) Degrees-of-freedom omitted by structural matrix partitioning
    # q  # ([nGdof,1] logical) Generalized degrees-of-freedom for dynamic reduction or component mode synthesis
    # r  # ([nGdof,1] logical) Reference degrees-of-freedom used to determine free body motion
    # c  # ([nGdof,1] logical) Degrees-of-freedom that are free during component mode synthesis or dynamic reduction
    # b  # ([nGdof,1] logical) Degrees-of-freedom fixed during component mode analysis or dynamic reduction
    # e  # ([nGdof,1] logical) extra degrees-of-freedom introduced in dynamic analysis
    # sa # Permanently constrained aerodynamic degrees-of-freedom
    # k  # Aerodynamic degrees-of-freedom

    # Nonexclusive Degrees-of-freedom sets
    # ------------------------------------
    # s  # ([nGdof,numSID] logical) [sb + sg] Degrees-of-freedom eliminated by single point constraints
    # l  # ([nGdof,1] logical) [b + c] Structural degrees-of-freedom remaining after the reference degrees-of-freedom are removed (degrees-of-freedom left over)
    #    # ([nGdof,1] logical) [l + r] Total set of physical boundary degrees-of-freedom for superelements
    #    # ([nGdof,1] logical) [t + q] Set assembled in superelement analysis
    # d  # ([nGdof,1] logical) [a + e] Set used in dynamic analysis by the direct method
    #    # ([nGdof,1] logical) [a + o] Unconstrained (free) structural degrees-of-freedom
    # fe # ([nGdof,1] logical) [f + e] Free structural degrees-of-freedom plus extra degrees-of-freedom
    #    # ([nGdof,1] logical) [f + s] Degrees-of-freedom not constrained by multipoint constraints
    # ne % ([nGdof,1] logical) [n + e] Structural degrees-of-freedom not constrained by multipoint constraints plus extra degrees-of-freedom
    # g = true(nGdof,1) [n + m] All structural degrees-of-freedom including scalar degrees-of-freedom
    # p = [g + e] Physical degrees-of-freedom
    # ps = [p + sa] Physical and constrained (SPCi) aerodynamic degrees-of-freedom
    # pa = [ps + k] Physical set for aerodynamics
    # fr = [f ? q ? r] Statically independent set minus the statically determinate supports
    # v = [o + c + r] Set free to vibrate in dynamic reduction and component mode synthesis
    return aset, tset

def apply_dof_map_to_set(set_map,
                         dof_map: DOF_MAP,
                         idtype: str='int32',
                         use_ints: bool=True) -> NDArrayNbool:
    """changes a set defined in terms of (nid, comp) into an array of integers"""
    if use_ints:
        ndof = len(set_map)
        aset = np.full(ndof, 0, dtype=idtype)
        for i, dofi in enumerate(set_map):
            aset[i] = dof_map[dofi]
        aset.sort()
    else:
        ndof = len(dof_map)
        aset = np.full(ndof, 0, dtype='bool')
        for dofi in set_map:
            i = dof_map[dofi]
            aset[i] = True
    return aset

def xg_to_xb(model, xg: NDArrayNfloat, ngrid: int, ndof_per_grid: int,
             inplace: bool=True) -> NDArrayNfloat:
    assert isinstance(xg, np.ndarray)
    str(ngrid)

    xb = xg
    if not inplace:
        xb = copy.deepcopy(xb)

    nids = model._type_to_id_map['GRID']
    for i, nid in enumerate(nids):
        node = model.nodes[nid]
        if node.cd:
            model.log.debug(f'node {nid} has a CD={node.cd}')
            cd_ref = node.cd_ref
            T = cd_ref.beta_n(n=2)
            i1 = i * ndof_per_grid
            i2 = (i+1) * ndof_per_grid
            xi = xg[i1:i2]
            xb[i1:i2] = xi @ T  # TODO: verify the transform; I think it's right
    return xb

def write_oload(Fb: NDArrayNfloat,
                dof_map: DOF_MAP,
                isubcase: int, ngrid: int, ndof_per_grid: int,
                f06_file: TextIO, page_stamp: str, page_num: int,
                log: SimpleLogger) -> int:
    """writes the OLOAD RESULTANT table"""
    str(ngrid)
    str(dof_map)
    fxyz_mxyz = Fb[:ngrid*ndof_per_grid].reshape(ngrid, ndof_per_grid)
    fxyz_mxyz_sum = fxyz_mxyz.sum(axis=0)

    f06_file.write(
        ' *** USER INFORMATION MESSAGE 7310 (VECPRN)\n'
        '     ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM WILL BE USED AS REFERENCE LOCATION.\n'
        '     RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN '
        'SUPERELEMENT BASIC SYSTEM COORDINATES.\n'
    )
    oload = Resultant('OLOAD', fxyz_mxyz_sum, isubcase)
    log.info(f'OLOAD {fxyz_mxyz_sum}')
    page_num = oload.write_f06(f06_file, page_stamp, page_num)
    return page_num + 1

def solve(Kaa: lil_matrix,
          Fa_solve: np.ndarray,
          aset: np.ndarray,
          log: SimpleLogger,
          idtype: str='int32'):
    """solves [K]{u} = {F}"""
    log.info('starting solve')
    Kaa_, ipositive, inegative, unused_sz_set = remove_rows(Kaa, aset, idtype=idtype)
    if len(ipositive) == 0:
        log.error('empty Kaa')
    #if np.linalg.det(Kaa) == 0.:
        #log.error('singular Kaa')

    isolve = len(ipositive)
    if isolve == 0:
        log.error(f'  ipositive = {ipositive}')
        log.error(f'  Kaa_ = {Kaa_}')
        #print(Kaa_)
        raise RuntimeError('no residual structure found')

    #Maa_ = Maa[ipositive, :][:, ipositive]
    #print(f'Fg = {Fg}')
    #print(f'Fa = {Fa}')
    #print(f'Fs = {Fs}')

    Fa_ = Fa_solve[ipositive]
    # [A]{x} = {b}
    # [Kaa]{x} = {F}
    # {x} = [Kaa][F]
    #print(f'Kaa:\n{Kaa}')
    #print(f'Fa: {Fa}')

    Kaa_dense = Kaa_.toarray()
    log.debug(f'  Kaas_:\n{Kaa_dense}')
    log.debug(f'  Kaa_:\n{Kaa_}')
    log.debug(f'  Fa_: {Fa_}')
    xa_ = np.linalg.solve(Kaa_dense, Fa_)
    xas_ = sci_sparse.linalg.spsolve(Kaa_, Fa_)
    #xas_ = np.linalg.solve(Kaas_.toarray(), Fa_)
    sparse_error = np.linalg.norm(xa_ - xas_)
    if sparse_error > 1e-12:
        log.warning(f'  sparse_error = {sparse_error}')
    log.info('finished solve')
    return xas_, ipositive, inegative

def build_Mbb(model: BDF,
              subcase: Subcase,
              dof_map: DOF_MAP,
              ndof: int, fdtype='float64') -> NDArrayNNfloat:
    """builds the mass matrix in the basic frame, [Mbb]"""
    log = model.log
    log.info('starting build_Mbb')
    wtmass = 1.
    #wtmass = model.get_param('WTMASS', 1.0)
    #Mbb = np.eye(ndof, dtype=fdtype)
    Mbb = np.zeros((ndof, ndof), dtype=fdtype)
    #Mbb = np.eye(ndof, dtype='int32')
    str(model)
    str(subcase)
    #no_mass = {
        #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
    #}

    mass_rod_2x2 = np.array([
        [2, 0, 1, 0],
        [0, 2, 0, 1],
        [1, 0, 2, 0],
        [0, 1, 0, 2],
    ], dtype='float64') / 3.
    mass_tri = np.array([
        [4, 0, 2, 0, 1, 0],
        [0, 4, 0, 2, 0, 1],
        [2, 0, 4, 0, 2, 0],
        [0, 2, 0, 4, 0, 2],
        [1, 0, 2, 0, 4, 0],
        [0, 1, 0, 2, 0, 4],
    ], dtype='float64') / 6.
    #mass_quad_1x1 = np.array([
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
        #[1, 0, 1, 0, 1, 0, 1, 0],
        #[0, 1, 0, 1, 0, 1, 0, 1],
    #], dtype='float64') / 36.

    mass_quad_2x2 = np.array([
        [4, 0, 2, 0, 1, 0, 2, 0],
        [0, 4, 0, 2, 0, 1, 0, 2],
        [2, 0, 4, 0, 2, 0, 1, 0],
        [0, 2, 0, 4, 0, 2, 0, 1],
        [1, 0, 2, 0, 4, 0, 2, 0],
        [0, 1, 0, 2, 0, 4, 0, 2],
        [2, 0, 1, 0, 2, 0, 4, 0],
        [0, 2, 0, 1, 0, 2, 0, 4],
    ], dtype='float64') / 36.

    mass_total = 0.
    if model.conm1:
        nid = elem.nid
        nid_ref = elem.nid_ref
        if nid_ref.type == 'GRID':
            i1 = dof_map[(nid, 1)]

            # TODO: support CID
            if nid_ref.cd != elem.cid:
                log.warning(f'  CONM1 eid={eid} nid={nid} CD={nid_ref.cd} to cid={elem.cid} is not supported')
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)
        Mbb[i1:i1+6, i1:i1+6] += elem.mass_matrix

    # has possibility of mass
    has_mass = False

    for elem in model.element_cards:
        etype = elem.type
        if etype in NO_MASS:
            continue
        if elem.n == 0:
            continue

        has_mass = True
        if etype == 'CONM2':
            mass_total = conm2_fill_Mbb(model, mass_total, Mbb, dof_map)

        elif etype in ['CROD', 'CONROD', 'CTUBE']:
            # verified
            mass = elem.mass()
            if mass == 0.0:
                log.warning(f'  no mass for {etype} eid={eid}')
                continue

            nids1 = elem.nodes[:, 0]
            nids2 = elem.nodes[:, 1]
            for nid1, nid2 in zip(nids1, nids2):
                i1 = dof_map[(nid1, 1)]
                j1 = dof_map[(nid2, 1)]
                ii = [i1, i1 + 1,
                      j1, j1 + 1]
                #Mbb[i1, i1] = Mbb[i1+1, i1+1] = \
                #Mbb[j1, j1] = Mbb[j1+1, j1+1] = mass / 3

                #Mbb[i1, j1] = Mbb[j1, i1] = \
                #Mbb[i1+1, j1+1] = Mbb[j1+1, i1+1] = mass / 6
                iii, jjj = np.meshgrid(ii, ii)
                Mbb[iii, jjj] += mass_rod_2x2 * mass
                #ii = [i1, i1 + 1, j1, j1 + 1]
                #print(Mbb[ii, :][:, ii])
        elif etype in {'CBAR', 'CBEAM'}:
            # TODO: verify
            # TODO: add rotary inertia
            mass = elem.mass()
            if mass.sum() == 0.0:
                log.warning(f'  no mass for {etype} eid={elem.element_id}')
                continue
            for (nid1, nid2), massi in zip(elem.nodes, mass):
                i1 = dof_map[(nid1, 1)]
                j1 = dof_map[(nid2, 1)]

                #Mbb[i1, i1] = Mbb[i1+1, i1+1] = \
                #Mbb[j1, j1] = Mbb[j1+1, j1+1] = mass / 3

                #Mbb[i1, j1] = Mbb[j1, i1] = \
                #Mbb[i1+1, j1+1] = Mbb[j1+1, i1+1] = mass / 6
                ii = [i1, i1 + 1,
                      j1, j1 + 1]
                iii, jjj = np.meshgrid(ii, ii)
                Mbb[iii, jjj] += mass_rod_2x2 * massi
        elif etype == 'CTRIA3':
            # TODO: verify
            # TODO: add rotary inertia
            mass = elem.Mass()
            nid1, nid2, nid3 = elem.nodes
            i1 = dof_map[(nid1, 1)]
            i2 = dof_map[(nid2, 1)]
            i3 = dof_map[(nid3, 1)]
            ii = [
                i1, i1 + 1,
                i2, i2 + 1,
                i3, i3 + 1,
            ]
            iii, jjj = np.meshgrid(ii, ii)
            Mbb[iii, jjj] += mass_tri * mass
            #Mbb[i1, i1] = Mbb[i1+1, i1+1] = Mbb[i1+2, i1+2] = \
            #Mbb[i2, i2] = Mbb[i2+1, i2+1] = Mbb[i2+2, i2+2] = \
            #Mbb[i3, i3] = Mbb[i3+1, i3+1] = Mbb[i3+2, i3+2] = mass / 3
        elif etype == 'CQUAD4':
            # TODO: verify
            # TODO: add rotary inertia
            try:
                mass = elem.Mass()
            except AttributeError:
                pid_ref = elem.pid_ref
                if pid_ref.mid1 is None and pid_ref.mid2 is None:
                    log.warning(f'  no mass for CQUAD4 eid={eid}')
                    continue
                #raise
                #mid_ref = elem.mid_ref
                #rho = mid_ref.Rho()
            nid1, nid2, nid3, nid4 = elem.nodes
            i1 = dof_map[(nid1, 1)]
            i2 = dof_map[(nid2, 1)]
            i3 = dof_map[(nid3, 1)]
            i4 = dof_map[(nid4, 1)]
            if mass == 0.:
                pid_ref = elem.pid_ref
                ptype = pid_ref.type
                if ptype == 'PSHELL':
                    mid_ref = pid_ref.mid_ref
                    rho = mid_ref.rho
                    log.warning(f'  no mass for CQUAD4 eid={eid} ptype={ptype} rho={rho}')
                else:
                    log.warning(f'  no mass for CQUAD4 eid={eid} ptype={ptype}')
            else:
                pid_ref = elem.pid_ref
                ptype = pid_ref.type
                log.info(f'  mass={mass} for CQUAD4 eid={eid} ptype={ptype}')

            ii = [
                i1, i1 + 1,
                i2, i2 + 1,
                i3, i3 + 1,
                i4, i4 + 1,
            ]
            iii, jjj = np.meshgrid(ii, ii)
            Mbb[iii, jjj] += mass_quad_2x2 * mass
            #if 0:  # pragma: no cover
                #mass4 = mass / 9. # 4/36
                #mass2 = mass / 18. # 2/36
                #mass1 = mass / 36.
                #print(mass1, mass2, mass4)
                #Mbb[i1, i1] += mass4
                #Mbb[i1+1, i1+1] += mass4
                #Mbb[i2, i2] += mass4
                #Mbb[i2+1, i2+1] += mass4
                #Mbb[i3, i3] += mass4
                #Mbb[i3+1, i3+1] += mass4
                #Mbb[i4, i4] += mass4
                #Mbb[i4+1, i4+1] += mass4

                #Mbb[i1, i3] += mass1
                #Mbb[i1, i2] += mass2
                #Mbb[i1, i4] += mass2
                #Mbb[i3, i1] += mass1
                #Mbb[i2, i1] += mass2
                #Mbb[i4, i1] += mass2

                #Mbb[i1+1, i3+1] += mass1
                #Mbb[i1+1, i2+1] += mass1
                #Mbb[i1+1, i4+1] += mass2
                #Mbb[i3+1, i1+1] += mass1
                #Mbb[i2+1, i1+1] += mass2
                #Mbb[i4+1, i1+1] += mass2

                #Mbb[i2, i4] += mass1
                #Mbb[i2+1, i4+1] += mass1
                #Mbb[i2, i3] += mass2
                #Mbb[i2+1, i3+1] += mass2

                #Mbb[i4, i2] += mass1
                #Mbb[i4+1, i2+1] += mass1
                #Mbb[i3, i2] += mass2
                #Mbb[i3+1, i2+1] += mass2

                #Mbb[i3, i4] += mass2
                #Mbb[i3+1, i4+1] = mass2
                #Mbb[i4, i3] += mass2
                #Mbb[i4+1, i3+1] += mass2

                #Mbb[i2+1, i3+1] = 1
                #Mbb[i2+1, i2+1] = Mbb[i1+1, i4+1] = 2

                #Mbb[i2, i1+1] = Mbb[i3, i1+1] = Mbb[i4, i1+1] = 1
                #Mbb[i1, i2+1] = Mbb[i3, i2+1] = Mbb[i4, i2+1] = 1
                #Mbb[i1, i3+1] = Mbb[i2, i3+1] = Mbb[i4, i3+1] = 1
                #Mbb[i1, i4+1] = Mbb[i2, i4+1] = Mbb[i3, i4+1] = 1
                #print(Mbb)
            #print(Mbb[ii, :][:, ii])
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)

    if wtmass != 1.0:
        Mbb *= wtmass

    # TODO: working on a hack to NOTt fake the mass when there are no SPOINTs
    #       as mass can be easily added to a spring model
    #
    # TODO: non-zero mass case doesn't handle SPOINTs
    #
    has_special_points = 'SPOINT' in model.card_count or 'EPOINT' in model.card_count
    is_all_grids = not has_special_points
    unused_can_dof_slice = is_all_grids and not has_mass
    if Mbb.sum() != 0.0:
    #if Mbb.sum() != 0.0 or can_dof_slice:
        #print(f'is_all_grids={is_all_grids} has_mass={has_mass}; can_dof_slice={can_dof_slice} Mbb.shape={Mbb.shape}')
        i = np.arange(0, ndof).reshape(ndof//6, 6)[:, :3].ravel()
        #print(Mbb[i, i])
        massi = Mbb[i, i].sum()
        log.info(f'finished build_Mbb; M={massi:.6g}; mass_total={mass_total:.6g}')
    else:
        Mbb = np.eye(ndof, dtype=fdtype)
        log.error(f'finished build_Mbb; faking mass; M={Mbb.sum()} ndof={ndof}')
    return Mbb

def grid_point_weight(model: BDF, Mbb, dof_map: DOF_MAP, ndof: int):
    str(dof_map)
    str(ndof)
    z = np.zeros((3, 3), dtype='float64')
    nnodes = model.grid.n
    nspoints = model.spoint.n
    #nepoints = len(model.epoints)
    ndof_grid = 6 * nnodes
    ndof2 = 6 * nnodes + nspoints
    #print(f'ndof={ndof} ndof2={ndof2}')
    D = np.zeros((ndof2, 6), dtype='float64')
    D[ndof_grid:, 0] = 1
    #print('D.shape ', D.shape)
    inid = 0
    #print(model.nodes)
    #print(f'D.shape = {D.shape}')
    reference_point = 0 #  model.get_param('GRDPNT', 0)
    if reference_point == 0:
        dxyz = np.zeros(3, dtype=Mbb.dtype)
    else:
        dxyz = model.nodes[reference_point].get_position()

    coord = model.coord
    xyz_cid0 = model.grid.xyz_cid0()
    cds = model.grid.cd
    cids = coord.coord_id
    #icds = coord.index(cds, assume_sorted=False, inverse=False)
    beta = coord.j
    for xyz0, cd in zip(xyz_cid0, cds):
        #coord = model.coord.s
        #print(f'nid={nid}')

        # TODO: what about otuput coordinate frames; specifically cylindrical and spherical frames?
        xi, yi, zi = xyz0 - dxyz
        Tr = np.array([
            [0, zi, -yi],
            [-zi, 0, xi],
            [yi, -xi, 0],
        ], dtype='float64')
        #cd_ref = node.cd_ref
        beta = coord.xyz_to_global_transform[cd]
        Ti = beta
        TiT_Tr = Ti.T @ Tr
        d = np.block([
            [Ti.T, TiT_Tr],
            [z, Ti.T],
        ])
        #print(f'd.shape={d.T.shape}')
        #print(f'd={d}')
        #print(f'D[:, {inid*6}:{(inid+1)*6}]')
        D[inid*6:(inid+1)*6, :] = d.T
        inid += 1

    # Location Basic System (CP Fields) Mass Global System (CD Fields)
    # Grid ID xb yb zb xCD yCD zCD
    # 1       0   0 0  2   3   5
    # 2       1   0 0  2   3   5
    # 3       0.5 1 0  2   3   5
    # 4       0.5 0 1  2   3   5
    #print(f'Dt.shape = {D.T.shape}')
    #print(f'Mbb.shape = {Mbb.shape}')
    #print(f'D.shape = {D.shape}')
    #print(f'D.T =\n{D.T}')
    M0 = D.T @ Mbb @ D
    return reference_point, M0

def dof_map_to_tr_set(dof_map, ndof: int) -> tuple[NDArrayNbool, NDArrayNbool]:
    """creates the translation/rotation sets"""
    trans_set = np.zeros(ndof, dtype='bool')
    rot_set = np.zeros(ndof, dtype='bool')
    for (unused_nid, dof), idof in dof_map.items():
        if dof in [0, 1, 2, 3]:
            trans_set[idof] = True
        else:
            rot_set[idof] = True
    return trans_set, rot_set

def ps_to_sg_set(ndof: int, ps: list[int]):
    """creates the SPC on the GRID (PS-field) set, {sg}"""
    # all DOFs are initially assumed to be active
    sg_set = np.ones(ndof, dtype='bool')

    # False means it's constrained
    sg_set[ps] = False
    return sg_set


def _get_dof_map(model: BDF) -> tuple[DOF_MAP, list[int]]:
    """helper method for ``get_static_force_vector_from_subcase_id``"""
    i = 0
    dof_map = {}
    spoints: list[int] = []
    ps: list[int] = []
    for nid, psi in zip(model.grid.node_id, model.grid.ps):
        for dof in range(1, 7):
            dof_map[(nid, dof)] = i
            i += 1
        if psi != 0:
            asf
            for psii in str(psi):
                nid_dof = (nid, int(psii))
                j = dof_map[nid_dof]
                ps.append(j)

    # we want the GRID points to be first
    assert len(spoints) == 0, spoints

    if model.spoint.n:
        for nid in model.spoint.spoint_id:
            key = (nid, 0)
            if key not in dof_map:
                dof_map[key] = i
                i += 1
    assert len(dof_map) > 0
    return dof_map, ps

def get_ndof(model: BDF, subcase: Subcase) -> tuple[int, int, int]:
    ndof_per_grid = 6
    if 'HEAT' in subcase:
        ndof_per_grid = 1
    ngrid = model.grid.n
    nspoint = model.spoint.n
    #nepoint = model.epoint.n
    nepoint = 0
    ndof = ngrid * ndof_per_grid + nspoint + nepoint
    #print(f'ngrid={ngrid} nspoint={nspoint}')
    assert ndof > 0, model.card_count
    return ngrid, ndof_per_grid, ndof

def conm2_fill_Mbb(model: BDF,
                   mass_total: float,
                   Mbb,
                   dof_map: dict[tuple[int, int], int]) -> float:
    eye3 = np.eye(3, dtype='float64')
    conm2 = model.conm2
    log = model.log
    #mass = elem.Mass()
    #nid = elem.nid
    #nid_ref = elem.nid_ref
    inid = model.grid.index(conm2.node_id)
    cds = model.grid.cd[inid]
    for eid, nid, cd, cid, mass, elem_x, elem_i in zip(
            conm2.element_id, conm2.node_id, cds, conm2.coord_id, conm2.mass(),
            conm2.xyz_offset, conm2.inertia):
        i1 = dof_map[(nid, 1)]
        if cd != cid:
            log.warning(f'  CONM2 eid={eid} nid={nid} CD={cd} to CONM2 cid={cid} is not supported')
        #Mbb[i1, i1] = mass
        #Mbb[i1+1, i1+1] = mass
        #Mbb[i1+2, i1+2] = mass
        # TODO: support CID
        I11, I21, I22, I31, I32, I33 = elem_i
        x1, x2, x3 = elem_x
        mxx = np.array([
            [x1 * x1, -x1 * x2, -x1 * x3],
            [-x2 * x1, x2 * x2, -x2 * x3],
            [-x3 * x1, x3 * x2, x3 * x3],
        ]) * mass
        Tr = np.array([
            [0, x3, -x2],
            [-x3, 0, x1],
            [x2, -x1, 0],
        ], dtype='float64')
        mx = Tr * mass
        I = np.array([
            [I11, -I21, I31],
            [-I21, I22, -I32],
            [-I31, -I32, I33],
        ]) + mxx

        #[mass, 01, 02, 03, mass * X3, -mass * X2]
        #[10, mass, 12, -mass * X3, 14, mass * X1]
        #[20, 21, mass, mass * X2, -mass * X1, 25]
        #[30, -mass * X3, mass * X2,        I11 + mass * X2 * X2 + mass * X3 * X3, -I21 - mass * X2 * X1,                  -I31 - mass * X3 * X1]
        #[mass * X3, 41, -mass * X1,       -I21 - mass * X2 * X1,                   I22 + mass * X1 * X1 + mass * X3 * X3, -I32 - mass * X3 * X2]
        #[-mass * X2, mass * X1, 52,       -I31 - mass * X3 * X1,                  -I32 - mass * X3 * X2,                   I33 + mass * X2 * X2 + mass * X1 * X1]

        Mbb[i1:i1+3, i1:i1+3] += eye3 * mass
        Mbb[i1:i1+3, i1+3:i1+6] += mx
        Mbb[i1+3:i1+6, i1:i1+3] += mx.T
        mass_total += mass
        Mbb[i1+3:i1+6, i1+3:i1+6] += I
    return mass_total

def _run_modes(model: BDF,
               subcase: Subcase,
               ngrid: int,
               ndof_per_grid: int,
               ndof: int,
               node_gridtype: np.ndarray,
               dof_map: dict,
               idtype: str = 'int32',
               fdtype: str = 'floaat64'):
    out = {}
    Kgg = build_Kgg(model, dof_map,
                    ndof, ngrid,
                    ndof_per_grid,
                    idtype='int32', fdtype=fdtype)
    out['Kgg'] = Kgg

    Mbb = build_Mbb(model, subcase, dof_map, ndof, fdtype=fdtype)

    Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid)
    out['Mgg'] = Mgg
    del Mbb

    #gset = np.arange(ndof, dtype=idtype)
    gset = np.ones(ndof, dtype='bool')
    sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)
    #aset = np.setdiff1d(gset, sset) # a = g-s
    aset = gset & ~sset_b # a = g-s
    out['aset'] = aset

    # aset - analysis set
    # sset - SPC set
    xa, xs = partition_vector2(xg, [['a', aset], ['s', sset_b]])
    out['xg'] = xg

    #del xg
    #print(f'xa = {xa}')
    #print(f'xs = {xs}')
    #print(Kgg)
    M = partition_matrix(Mgg, [['a', aset], ['s', sset_b]])
    Maa = M['aa']
    out['Maa'] = Maa

    K = partition_matrix(Kgg, [['a', aset], ['s', sset_b]])
    Kaa = K['aa']
    out['Kaa'] = Kaa
    #Kss = K['ss']
    #Kas = K['as']
    #Ksa = K['sa']

    #[Kaa]{xa} + [Kas]{xs} = {Fa}
    #[Ksa]{xa} + [Kss]{xs} = {Fs}

    #{xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
    #{Fs} = [Ksa]{xa} + [Kss]{xs}

    # TODO: apply AUTOSPCs correctly
    #print(Kaa)
    #print(Kas)
    #print(Kss)
    #Maa_, ipositive, inegative, unused_sz_set = remove_rows(Maa, aset)
    Kaa_, ipositive, unused_inegative, unused_sz_set = remove_rows(Kaa, aset)
    Maa_ = Maa[ipositive, :][:, ipositive]
    #Fs = np.zeros(ndof, dtype=fdtype)
    #print(f'Fg = {Fg}')
    #print(f'Fa = {Fa}')
    #print(f'Fs = {Fs}')
    #Fa_ = Fa[ipositive]
    # [A]{x} = {b}
    # [Kaa]{x} = {F}
    # {x} = [Kaa][F]
    #print(f'Kaa:\n{Kaa}')
    #print(f'Fa: {Fa}')

    #print(f'Kaa_:\n{Kaa_}')
    #print(f'Fa_: {Fa_}')
    #na = Kaa_.shape[0]
    ndof_ = Kaa_.shape[0]

    method_id, options = subcase['METHOD']
    assert isinstance(method_id, int), method_id
    method = model.methods[method_id]
    if method.type == 'EIGRL':
        neigenvalues = method.nd  # nroots
    #neigenvalues = 10

    eigenvalues, xa_ = solve_eigenvector(Kaa_, Maa_, ndof, neigenvalues)
    nmodes = len(eigenvalues)
    model.log.debug(f'eigenvalues = {eigenvalues}')
    out['eigenvalues'] = eigenvalues
    #print(f'xa_ = {xa_} {xa_.shape}')
    #xa2 = xa_.reshape(nmodes, na, na)

    ndof_a = len(xa)
    xg_out = np.full((nmodes, ndof), 0., dtype=fdtype)
    #xa_out = xg_out[:, aset]
    #assert xa_out.shape == (nmodes, ndof_a)
    #xa_out = np.full((nmodes, ndof_a), np.nan, dtype=fdtype)
    xa_out = np.zeros((nmodes, ndof_a), dtype=fdtype)
    #xa[ipositive] = xa_
    for imode in range(nmodes):
        xa_out[imode, ipositive] = xa_[:, imode]
    #xa_out[:, ipositive] = xa_
    #xg = np.arange(ndof, dtype=fdtype)
    #xg[aset] = xa
    #xg[sset] = xs
    xg_out[:, aset] = xa_out
    xg_out[:, sset] = xs

    assert np.all(np.isfinite(xa_out))
    assert np.all(np.isfinite(xg_out))
    out['xg_out'] = xg_out
    return out

def solve_eigenvector(Kaa: csc_matrix, Maa: np.ndarray,
                      ndof: int, neigenvalues: int) -> tuple[np.ndarray, np.ndarray]:
    Kaa_dense = Kaa.toarray()
    assert isinstance(Maa, np.ndarray), type(Maa)
    #Maa_dense = Maa.toarray()

    # reduce the size of the matrix going into the solver
    # by removing empty rows/columns
    Kaai = Kaa_dense
    is_kaa = (np.abs(Kaai) > 0.)
    is_maa = (np.abs(Maa) > 0.)
    is_kaa_maa = (is_kaa & is_maa)
    rows = np.any(is_kaa_maa, axis=0)
    cols = np.any(is_kaa_maa, axis=1)
    assert len(rows) == len(cols)
    is_modes = rows & cols
    no_modes = ~is_modes
    Kaa2 = Kaa[is_modes, :][:, is_modes]
    Maa2 = Maa[is_modes, :][:, is_modes]

    ndof2 = Maa.shape[0]
    if ndof2 <= neigenvalues:
        Kaa2_dense = Kaa2.todense()
        eigenvalues, xa = sp.linalg.eigh(Kaa2_dense, Maa2)
    else:
        #If M is specified, solves ``A * x[i] = w[i] * M * x[i]``
        # which = SM (smallest magnitude eigenvalues)
        eigenvalues, xa_ = sp.sparse.linalg.eigsh(
            Kaa2, k=neigenvalues, M=Maa2,
            sigma=None, which='SM', v0=None, ncv=None, maxiter=None, tol=0,
            return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
        xa = np.full((ndof, neigenvalues), np.nan, dtype=xa_.dtype)
        xa[no_modes, :] = 0
        xa[is_modes, :] = xa_
    return eigenvalues, xa

def _build_xg(model: BDF, dof_map: DOF_MAP,
              ndof: int, subcase: Subcase) -> NDArrayNfloat:
    """
    Builds the {xg} vector, which has all SPCs in the analysis (cd) frame
    (called global g by NASTRAN)

    {s} = {sb} + {sg}
    {sb} = SPC set on SPC/SPC1/SPCADD cards (boundary)
    {sg} = SPCs from PS field on GRID card (grid)

    """
    log = model.log
    #model = self.model
    xspc = np.full(ndof, np.nan, dtype='float64')
    #get_parameter(self, param_name, msg='', obj=False)
    spc_id, unused_options = subcase['SPC']
    if 'SPC' not in subcase:
        log.warning(f'no spcs...{spc_id}')
        log.warning(str(subcase))
        return xspc
    spc_id, unused_options = subcase['SPC']
    spcs = []
    #for spc in model.spcs:
        #if spc.n == 0:
            #continue
        #spci = spc.slice_card_by_id(spc_id)
        #spcs.append(spci)
    spc_cards = [spc for spc in model.spc_cards if spc.n > 0]
    spcs = [spc.slice_card_by_id(spc_id, sort_ids=True) for spc in spc_cards]
    model.spc1
    #spcs = model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)

    spc_set = []
    sset = np.zeros(ndof, dtype='bool')
    for spc in spcs:
        if spc.type == 'SPC1':
            #print(spc.get_stats())
            dofs_missed = []
            for comp, (inode0, inode1) in zip(spc.components, spc.inode):
                for dof in str(comp):
                    dofi = int(dof)
                    spc_nodes = spc.node_id[inode0:inode1]
                    for nid in spc_nodes:
                        try:
                            idof = dof_map[(nid, dofi)]
                        except Exception:
                            dofs_missed.append((nid, dofi))
                            #print('dof_map =', dof_map)
                            #print((nid, dofi))
                            continue
                        sset[idof] = True
                        spc_set.append(idof)
                        xspc[idof] = 0.
            if dofs_missed:
                dof_str = ', '.join(str(dofi) for dofi in dofs_missed)
                log.warning(f'Missing (nid,dof) pairs:')
                log.warning(f'  {dof_str}\n{spc.rstrip()}')
                #=({nid},{dofi}) doesn't exist...skipping
            del dofs_missed
        elif spc.type == 'SPC':
            for nid, components, enforcedi in zip(spc.node_id, spc.components, spc.enforced):
                for component in str(components):
                    dofi = int(component)
                    idof = dof_map[(nid, dofi)]
                    sset[idof] = True
                    spc_set.append(idof)
                    xspc[idof] = enforcedi
        else:
            raise NotImplementedError(spc)
    spc_set = np.array(spc_set, dtype='int32')
    #print('spc_set =', spc_set, xspc)
    return spc_set, sset, xspc


def get_frequencies(model: BDF, subcase: Subcase,
                    omega_ns: np.ndarray) -> np.ndarray:

    is_freqs = (len(model.frequencies) == 0)
    is_frequency = ('FREQUENCY' not in subcase)
    natural_freq = np.unique(omega_ns) / (2 * np.pi)
    del omega_ns

    if is_freqs or is_frequency:
        nfreq = 1001
        #freq_max = 2 * np.pi * omega_max
        fmax_default = 1.5 * natural_freq.max()

        #freq = np.linspace(0., 1.5*freq_max, num=nfreq)
        frequencies = np.linspace(1., fmax_default, num=nfreq)
    else:
        freq_id, unused_options = subcase['FREQUENCY']
        frequencies_list = []
        for freq in model.frequencies[freq_id]:
            if freq.type in {'FREQ', 'FREQ1', 'FREQ2'}:
                frequencies_list.append(freq.freqs)
            elif freq.type in {'FREQ3', 'FREQ4', 'FREQ5'}:
                freqi = freq.get_frequencies(natural_freq)
                frequencies_list.append(freqi)
            else:  # pragma: no cover
                raise RuntimeError(freq)
        frequencies = np.unique(np.hstack(frequencies_list))
    return frequencies
