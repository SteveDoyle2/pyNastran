import copy
from datetime import date
from collections import defaultdict
from itertools import count
from typing import List, Dict, Tuple, Union, Any

import numpy as np
import scipy.sparse as sci_sparse

import pyNastran
from pyNastran.bdf.bdf import BDF, Subcase
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_interface.hdf5_interface import (
    RealDisplacementArray, RealSPCForcesArray, RealLoadVectorArray)
from pyNastran.bdf.mesh_utils.loads import _get_dof_map, get_ndof

from pyNastran.dev.solver.recover_static_strains import recover_strain_101
from pyNastran.dev.solver.build_stiffness import build_Kbb

class Solver:
    """defines the Nastran knockoff class"""
    def __init__(self, model: BDF):
        self.model = model
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
        self.f06_filename = 'junk.f06'

    def run(self):
        model = self.model
        sol = model.sol
        solmap = {
            101 : self.run_sol_101,  # static
            103 : self.run_sol_103,  # modes
        }
        model.cross_reference()
        self._update_card_count()

        title = ''
        today = None
        page_stamp = self.op2.make_stamp(title, today)
        with open(self.f06_filename, 'w') as f06_file:
            f06_file.write(self.op2.make_f06_header())
            self.op2._write_summary(f06_file, card_count=model.card_count)
            f06_file.write('\n')
            if sol in [101, 103, 105, 107, 109, 111, 112]:
                for subcase_id, subcase in sorted(model.subcases.items()):
                    if subcase_id == 0:
                        continue
                    self.log.debug(f'subcase_id={subcase_id}')
                    runner = solmap[sol]
                    runner(subcase, f06_file, page_stamp,
                           idtype='int32', fdtype='float64')
            else:
                raise NotImplementedError(sol)

    def _update_card_count(self):
        for card_type, values in self.model._type_to_id_map.items():
            self.model.card_count[card_type] = len(values)

    def build_xg(self, dof_map: Dict[Any, int], ndof: int, subcase: Subcase) -> np.ndarray:
        """
        Builds the {xg} vector, which has all SPCs in the analysis (cd) frame
        (called global g by NASTRAN)

        {s} = {sb} + {sg}
        {sb} = SPC set on SPC/SPC1/SPCADD cards (boundary)
        {sg} = SPCs from PS field on GRID card (grid)

        """
        model = self.model
        xspc = np.full(ndof, np.nan, dtype='float64')
        #get_parameter(self, param_name, msg='', obj=False)
        spc_id, unused_options = subcase['SPC']
        if 'SPC' not in subcase:
            model.log.warning(f'no spcs...{spc_id}')
            model.log.warning(str(subcase))
            return xspc
        spc_id, unused_options = subcase['SPC']
        spcs = model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)

        spc_set = []
        sset = np.zeros(ndof, dtype='bool')
        for spc in spcs:
            if spc.type == 'SPC1':
                #print(spc.get_stats())
                for dofi in spc.components:
                    dofi = int(dofi)
                    for nid in spc.nodes:
                        try:
                            idof = dof_map[(nid, dofi)]
                        except:
                            print(spc)
                            print('dof_map =', dof_map)
                            print((nid, dofi))
                        sset[idof] = True
                        spc_set.append(idof)
                        xspc[idof] = 0.
        spc_set = np.array(spc_set, dtype='int32')
        #print('spc_set =', spc_set, xspc)
        return spc_set, sset, xspc


    def build_Fb(self, dof_map: Dict[Any, int], ndof: int, subcase: Subcase) -> np.array:
        model = self.model
        Fb = np.zeros(ndof, dtype='float32')
        if 'LOAD' not in subcase:
            return Fb

        load_id, unused_options = subcase['LOAD']
        #print('load_id =', load_id)
        loads, scales, is_grav = model.get_reduced_loads(
            load_id, scale=1., consider_load_combinations=True,
            skip_scale_factor0=False, stop_on_failure=True, msg='')
        #loads : List[loads]
            #a series of load objects
        #scale_factors : List[float]
            #the associated scale factors
        #is_grav : bool
            #is there a gravity card
        for load, scale in zip(loads, scales):
            if load.type == 'SLOAD':
                #print(load.get_stats())
                for mag, nid in zip(load.mags, load.nodes):
                    i = dof_map[(nid, 1)]  # TODO: wrong...
                    Fb[i] = mag * scale
            elif load.type == 'FORCE':
                fxyz = load.to_global()
                nid = load.node
                self.log.debug(f'FORCE nid={nid} Fxyz={fxyz}')
                for i, dof in enumerate([1, 2, 3]):
                    # TODO: wrong because it doesn't handle SPOINTs
                    fi = dof_map[(nid, dof)]
                    Fb[fi] = fxyz[i]
            elif load.type == 'MOMENT':
                fxyz = load.to_global()
                nid = load.node
                self.log.debug(f'MOMENT nid={nid} Fxyz={fxyz}')
                for i, dof in enumerate([4, 5, 6]):
                    # TODO: wrong because it doesn't handle SPOINTs
                    fi = dof_map[(nid, dof)]
                    Fb[fi] = fxyz[i]
            else:
                print(load.get_stats())
                raise NotImplementedError(load)
        #print(subcase)
        return Fb

    def xg_to_xb(self, xg: np.ndarray, ngrid: int, ndof_per_grid: int, inplace=True):
        assert isinstance(xg, np.ndarray)
        model = self.model

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

    def Kbb_to_Kgg(self, Kbb: np.ndarray, ngrid: int, ndof_per_grid: int, inplace=True) -> np.ndarray:
        """does an in-place transformation"""
        assert isinstance(Kbb, (np.ndarray, sci_sparse.csc.csc_matrix)), type(Kbb)
        if not isinstance(Kbb, np.ndarray):
            Kbb = Kbb.tolil()
        model = self.model
        assert ngrid > 0, model.card_count
        nids = model._type_to_id_map['GRID']

        Kgg = Kbb
        if not inplace:
            Kgg = copy.deepcopy(Kgg)

        for i, nid in enumerate(nids):
            node = model.nodes[nid]
            if node.cd:
                model.log.debug(f'node {nid} has a CD={node.cd}')
                cd_ref = node.cd_ref
                T = cd_ref.beta_n(n=2)
                i1 = i * ndof_per_grid
                i2 = (i+1) * ndof_per_grid
                Ki = Kbb[i1:i2, i1:i2]
                Kgg[i1:i2, i1:i2] = T.T @ Ki @ T
        return Kgg

    def get_mpc_constraints(self, subcase: Subcase, dof_map):
        model = self.model
        #Fb = np.zeros(ndof, dtype='float32')
        constraints = []
        if 'MPC' not in subcase:
            return constraints

        mpc_id, unused_options = subcase['MPC']
        mpcs = model.get_reduced_mpcs(mpc_id, consider_mpcadd=True, stop_on_failure=True)
        #print('mpc_id =', mpc_id)

        ieq = 0
        ieqs = []
        dofs = []
        coefficients = []
        dependents = []
        independents = []
        independents_eq = defaultdict(list)
        for element in model.rigid_elements.values():
            #etype = element.type
            #ieq += 1
            raise NotImplementedError(element.get_stats())

        for mpc in mpcs:
            mpc_type = mpc.type
            if mpc_type == 'MPC':
                # The first degree-of-freedom (G1, C1) in the sequence is defined to
                # be the dependent DOF.
                # A dependent DOF assigned by one MPC entry cannot be assigned
                #  dependent by another MPC entry or by a rigid element.

                #: Component number. (Any one of the Integers 1 through 6 for grid
                #: points; blank or zero for scalar points.)
                #self.components = components

                #: Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
                #self.coefficients = coefficients
                #mpc.
                for i, nid, component, coeff in zip(count(), mpc.nodes, mpc.components, mpc.coefficients):
                    self.log.debug(f'ieq={ieq} (g,c)={(nid, component)} coeff={coeff}')
                    assert isinstance(component, int), component
                    idof = dof_map[(nid, component)]
                    ieqs.append(i)
                    dofs.append(idof)
                    coefficients.append(coeff)
                    if i == 0:
                        dependents.append(idof)
                    else:
                        independents.append(idof)
                        independents_eq[ieq].append(idof)
            else:
                raise NotImplementedError(mpc.get_stats())
            ieq += 1

        ieqs = np.array(ieqs, dtype='int32')
        independents = np.array(independents, dtype='int32')
        dependents = np.array(dependents, dtype='int32')
        independents_dependents = np.hstack([independents, dependents])
        coefficients = np.array(coefficients, dtype='float64')
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
                self.log.debug(f'dof={udof} N={n_dof}')
                if n_dof == 1:
                    break
            else:
                raise RuntimeError('cannot find a DOF to remove')
            #np.where(condition)

            break
        #for udof, in udofs:

        #uindependents, nicount = np.unique(independents, return_counts=True)
        #print(udofs, idofs)
        return constraints

    def run_sol_101(self, subcase: Subcase, f06_file,
                    page_stamp, idtype='int32', fdtype='float64'):
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
        d = a + e
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
        #basic
        itime = 0
        ntimes = 1  # static
        isubcase = subcase.id
        title = f'pyNastran {pyNastran.__version__}'
        subtitle = f'SUBCASE {isubcase}'
        label = ''
        if 'TITLE' in subcase:
            title = subcase.get_parameter('TITLE')
        if 'SUBTITLE' in subcase:
            subtitle = subcase.get_parameter('SUBTITLE')
        if 'LABEL' in subcase:
            label = subcase.get_parameter('LABEL')
        #-----------------------------------------------------------------------
        model = self.model
        log = model.log

        dof_map = _get_dof_map(model)
        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)
        Kbb, Kbbs = build_Kbb(model, subcase, dof_map, ndof)
        self.get_mpc_constraints(subcase, dof_map)

        Kgg = self.Kbb_to_Kgg(Kbb, ngrid, ndof_per_grid, inplace=False)
        Kggs = self.Kbb_to_Kgg(Kbbs, ngrid, ndof_per_grid)
        del Kbb, Kbbs

        gset = np.arange(ndof, dtype=idtype)
        gset_b = np.ones(ndof, dtype='bool')
        sset, sset_b, xg = self.build_xg(dof_map, ndof, subcase)


        # Constrained set
        # sset = sb_set | sg_set

        # free structural DOFs
        #fset = ~sset; # & ~mset
        fset_b = gset_b & sset_b # g-b
        self.log.debug(f'gset_b = {gset_b}')
        self.log.debug(f'sset_b = {sset_b}')
        self.log.debug(f'fset_b = {fset_b}')
        aset, tset = get_residual_structure(model, dof_map, fset_b)

        asetmap = get_aset(model)
        if asetmap:
            aset = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype)
        else:
            aset = np.setdiff1d(gset, sset) # a = g-s

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
        self.sset = sset
        self.aset = aset

        Fb = self.build_Fb(dof_map, ndof, subcase)
        Fg = Fb
        Fa, Fs = partition_vector(Fb, [['a', aset], ['s', sset]])
        del Fb
        # Mgg = self.build_Mgg(subcase)

        # aset - analysis set
        # sset - SPC set
        #print('aset = ', aset)
        #print('sset = ', sset)
        xa, xs = partition_vector(xg, [['a', aset], ['s', sset]])
        del xg
        #print(f'xa = {xa}')
        #print(f'xs = {xs}')
        #print(Kgg)
        self.Kgg = Kgg
        K = partition_matrix(Kgg, [['a', aset], ['s', sset]])
        Kaa = K['aa']
        Kss = K['ss']
        #Kas = K['as']
        Ksa = K['sa']
        Ks = partition_matrix(Kggs, [['a', aset], ['s', sset]])
        Kaas = Ks['aa']

        self.Kaa = Kaa
        #[Kaa]{xa} + [Kas]{xs} = {Fa}
        #[Ksa]{xa} + [Kss]{xs} = {Fs}

        #{xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
        #{Fs} = [Ksa]{xa} + [Kss]{xs}

        #print(Kaa)
        #print(Kas)
        #print(Kss)
        Kaas_, ipositive, inegative, sz_set = remove_rows(Kaas, aset, idtype=idtype)
        Kaa_, ipositive, inegative, sz_set = remove_rows(Kaa, aset, idtype=idtype)
        Fs = np.zeros(ndof, dtype=fdtype)
        #print(f'Fg = {Fg}')
        #print(f'Fa = {Fa}')
        #print(f'Fs = {Fs}')
        Fa_ = Fa[ipositive]
        # [A]{x} = {b}
        # [Kaa]{x} = {F}
        # {x} = [Kaa][F]
        #print(f'Kaa:\n{Kaa}')
        #print(f'Fa: {Fa}')

        log.debug(f'Kaa_:\n{Kaa_}')
        log.debug(f'Fa_: {Fa_}')
        xa_ = np.linalg.solve(Kaa_, Fa_)
        xas_ = sci_sparse.linalg.spsolve(Kaas_, Fa_)
        sparse_error = np.linalg.norm(xa_ - xas_)
        if sparse_error > 1e-12:
            log.warning(f'sparse_error = {sparse_error}')

        self.xa_ = xa
        log.info(f'xa_ = {xa_}')

        xa[ipositive] = xa_
        xa[inegative] = 0.

        xg = np.full(ndof, np.nan, dtype=fdtype)
        xg[aset] = xa
        xg[sset] = xs
        Fg[aset] = Fa

        fspc = np.full(ndof, 0., dtype=fdtype)
        fspci = Ksa @ xa + Kss @ xs
        Fg[sset] = fspci
        fspc[sset] = fspci

        xb = self.xg_to_xb(xg, ngrid, ndof_per_grid)
        Fb = self.xg_to_xb(Fg, ngrid, ndof_per_grid)

        log.debug(f'Fs = {Fs}')
        page_num = 1

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
            node_gridtype, Fg,
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
        if 'STRAIN' in subcase:
            recover_strain_101(f06_file, op2, self.model, dof_map, isubcase, xb,
                               title=title, subtitle=subtitle, label=label,
                               page_stamp=page_stamp)
        #if 'STRESS' in subcase:
            #recover_stress_101(self.model, xb, dof_map)
        #Fg[sz_set] = -1
        #xg[sz_set] = -1
        log.debug(f'xa = {xa}')
        op2.write_op2('junk.op2', post=-1, endian=b'<', skips=None, nastran_format='nx')
        return xa_

    def _save_displacment(self, f06_file,
                          subcase: Subcase, itime: int, ntimes: int,
                          node_gridtype, xg,
                          ngrid: int, ndof_per_grid: int,
                          title='', subtitle='', label='',
                          fdtype='float32', page_num=1, page_stamp='PAGE %s') -> int:
        f06_request_name = 'DISPLACEMENT'
        table_name = 'OUGV1'
        #self.log.debug(f'xg = {xg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, xg,
            RealDisplacementArray, f06_request_name, table_name,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

    def _save_spc_forces(self, f06_file,
                         subcase: Subcase, itime: int, ntimes: int,
                         node_gridtype, fspc,
                         ngrid: int, ndof_per_grid: int,
                         title='', subtitle='', label='',
                         fdtype='float32', page_num=1, page_stamp='PAGE %s') -> int:
        f06_request_name = 'SPCFORCES'
        table_name = 'OQG1'
        #self.log.debug(f'Fg = {Fg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, fspc,
            RealSPCForcesArray, f06_request_name, table_name,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

    def _save_applied_load(self, f06_file,
                           subcase: Subcase, itime: int, ntimes: int,
                           node_gridtype, Fg,
                           ngrid: int, ndof_per_grid: int,
                           title='', subtitle='', label='',
                           fdtype='float32', page_num=1, page_stamp='PAGE %s') -> int:
        f06_request_name = 'OLOAD'
        table_name = 'OPG1'
        #self.log.debug(f'Fg = {Fg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, Fg,
            RealLoadVectorArray, f06_request_name, table_name,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

    def _save_static_table(self, f06_file,
                           subcase: Subcase, itime: int, ntimes: int,
                           node_gridtype, Fg,
                           obj: Union[RealSPCForcesArray], f06_request_name, table_name,
                           ngrid: int, ndof_per_grid: int,
                           title='', subtitle='', label='',
                           fdtype='float32', page_num=1,
                           page_stamp='PAGE %s') -> int:
        if f06_request_name not in subcase:
            return
        isubcase = subcase.id
        self.log.debug(f'saving {f06_request_name} -> {table_name}')
        unused_nids_write, write_f06, write_op2 = get_plot_request(
            subcase, f06_request_name)
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
        self.op2.spc_forces[isubcase] = spc_forces
        return page_num

    def run_sol_103(self, subcase: Subcase, f06_file,
                    page_stamp,
                    idtype='int32', fdtype='float64'):
        """
        [M]{xdd} + [C]{xd} + [K]{x} = {F}
        [M]{xdd} + [K]{x} = {F}
        [M]{xdd}λ^2 + [K]{x} = {0}
        {X}(λ^2 + [M]^-1[K]) = {0}
        λ^2 + [M]^-1[K] = {0}
        λ^2 = -[M]^-1[K]
        [A][X] = [X]λ^2
        """
        log = self.model.log

        dof_map = _get_dof_map(self.model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)
        Kbb, Kbbs = build_Kbb(self.model, subcase, dof_map, ndof)

        Mbb = np.eye(Kbb.shape[0], dtype=fdtype)

        Kgg = self.Kbb_to_Kgg(Kbb, ngrid, ndof_per_grid)
        Mgg = self.Kbb_to_Kgg(Mbb, ngrid, ndof_per_grid)
        Kggs = self.Kbb_to_Kgg(Kbbs, ngrid, ndof_per_grid)
        del Kbb, Mgg

        gset = np.arange(ndof, dtype=idtype)
        sset, xg = self.build_xg(dof_map, ndof, subcase)
        aset = np.setdiff1d(gset, sset) # a = g-s

        # Mgg = self.build_Mgg(subcase)

        # aset - analysis set
        # sset - SPC set
        xa, xs = partition_vector(xg, [['a', aset], ['s', sset]])
        del xg
        #print(f'xa = {xa}')
        #print(f'xs = {xs}')
        #print(Kgg)
        K = partition_matrix(Kgg, [['a', aset], ['s', sset]])
        Kaa = K['aa']
        #Kss = K['ss']
        #Kas = K['as']
        #Ksa = K['sa']
        Ks = partition_matrix(Kggs, [['a', aset], ['s', sset]])
        Kaas = Ks['aa']

            #[Kaa]{xa} + [Kas]{xs} = {Fa}
            #[Ksa]{xa} + [Kss]{xs} = {Fs}

        #{xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
        #{Fs} = [Ksa]{xa} + [Kss]{xs}

        # TODO: apply AUTOSPCs correctly
        #print(Kaa)
        #print(Kas)
        #print(Kss)
        Kaa_, ipositive, inegative, unused_sz_set = remove_rows(Kaa, aset)
        Kaas_, ipositive, inegative, unused_sz_set = remove_rows(Kaas, aset)
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
        xa_ = np.linalg.eigh(Kaa_)
        #print(f'xa_ = {xa_}')

        xa[ipositive] = xa_
        xg = np.arange(ndof, dtype=fdtype)
        xg[aset] = xa
        xg[sset] = xs
        #fspc = Ksa @ xa + Kss @ xs
        #Fs[ipositive] = Fsi

        #Fg[aset] = Fa
        #Fg[sset] = fspc
        log.debug(f'xa = {xa}')
        #log.debug(f'Fs = {Fs}')
        log.debug(f'xg = {xg}')
        #log.debug(f'Fg = {Fg}')

        f06_file
        page_stamp
        return xa_
        #raise NotImplementedError(subcase)

    def run_sol_111(self, subcase: Subcase, f06_file, idtype='int32', fdtype='float64'):
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


def partition_matrix(matrix, sets) -> Dict[Tuple[str, str], np.ndarray]:
    """partitions a matrix"""
    matrices = {}
    for aname, aset in sets:
        for bname, bset in sets:
            matrices[aname + bname] = matrix[aset, :][:, bset]
    return matrices

def partition_vector(vector, sets) -> List[np.ndarray]:
    """partitions a vector"""
    vectors = []
    for unused_aname, aset in sets:
        vectori = vector[aset]
        vectors.append(vectori)
    return vectors


def remove_rows(Kgg: np.ndarray, aset: np.ndarray, idtype='int32') -> np.ndarray:
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
    elif isinstance(Kgg, (sci_sparse.csc.csc_matrix, sci_sparse.lil.lil_matrix)):
        ipositive1 = np.atleast_1d(col_kgg.todense()).nonzero()
        ipositive2 = np.atleast_1d(row_kgg.todense()).nonzero()
        ipositive = np.union1d(ipositive1, ipositive2)
        if isinstance(Kgg, sci_sparse.lil.lil_matrix):
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

def _get_node_gridtype(model: BDF, idtype='int32') -> np.ndarray:
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
    for nid, node_ref in sorted(model.nodes.items()):
        if node_ref.type == 'GRID':
            node_gridtype.append((nid, 1))
        elif node_ref.type == 'SPOINT':
            node_gridtype.append((nid, 2))
        else:
            raise NotImplementedError(node_ref)
    assert len(node_gridtype) > 0
    #print(node_gridtype)
    node_gridtype_array = np.array(node_gridtype, dtype=idtype)
    return node_gridtype_array

def get_plot_request(subcase: Subcase, request: str) -> Tuple[str, bool, bool]:
    """
    request = 'SPCFORCES'
    """
    value, options = subcase.get_parameter(request)
    write_f06 = False
    write_f06 = True
    if 'PRINT' in options:
        write_f06 = True
    if 'PLOT' in options:
        write_op2 = True
    if not(write_f06 or write_op2):
        write_op2 = True
    nids_write = value
    return nids_write, write_f06, write_op2

def get_aset(model: BDF) -> Any:
    aset_map = set()
    for aset in model.asets:
        if aset.type == 'ASET1':
            comp = aset.components
            for nid in aset.ids:
                for compi in comp:
                    aset_map.add((nid, int(compi)))
        elif aset.type == 'ASET':
            for nid, comp in zip(aset.ids, aset.components):
                for compi in comp:
                    aset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(aset)
    return aset_map

def get_bset(model: BDF) -> Any:
    """creates the b-set"""
    bset_map = set()
    for bset in model.bsets:
        if bset.type == 'BSET1':
            comp = bset.components
            for nid in bset.ids:
                for compi in comp:
                    bset_map.add((nid, int(compi)))
        elif bset.type == 'BSET':
            for nid, comp in zip(bset.ids, bset.components):
                for compi in comp:
                    bset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(bset)
    return bset_map

def get_cset(model: BDF) -> Any:
    """creates the c-set"""
    cset_map = set()
    for cset in model.csets:
        if cset.type == 'CSET1':
            comp = cset.components
            for nid in cset.ids:
                for compi in comp:
                    cset_map.add((nid, int(compi)))
        elif cset.type == 'CSET':
            for nid, comp in zip(cset.ids, cset.components):
                for compi in comp:
                    cset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(cset)
    return cset_map

def get_omit_set(model: BDF) -> Any:
    """creates the o-set"""
    omit_set_map = set()
    for omit in model.omits:
        if omit.type == 'OMIT1':
            comp = omit.components
            for nid in omit.ids:
                for compi in comp:
                    omit_set_map.add((nid, int(compi)))
        elif omit.type == 'OMIT':
            for nid, comp in zip(omit.ids, omit.components):
                for compi in comp:
                    omit_set_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(omit)
    return omit_set_map

def get_rset(model: BDF) -> Any:
    """creates the r-set"""
    rset_map = set()
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

def get_qset(model: BDF) -> Any:
    """creates the r-set"""
    qset_map = set()
    for qset in model.qsets:
        if qset.type == 'QSET1':
            comp = qset.components
            for nid in qset.ids:
                for compi in comp:
                    qset_map.add((nid, int(compi)))
        elif qset.type == 'QSET':
            for nid, comp in zip(qset.ids, qset.components):
                for compi in comp:
                    qset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(qset)
    return qset_map

def get_residual_structure(model: BDF, dof_map, fset, idtype='int32'):
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
    elif (not np.any(aqo_set) and np.any(bcr_set)):
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

def apply_dof_map_to_set(set_map, dof_map, idtype='int32', use_ints=True):
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
