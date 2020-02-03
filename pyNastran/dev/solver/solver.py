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
    RealDisplacementArray, RealSPCForcesArray)
from pyNastran.bdf.mesh_utils.loads import _get_dof_map, get_ndof


class Solver:
    """defines the Nastran knockoff class"""
    def __init__(self, model: BDF):
        self.model = model
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

    def build_Kbb(self, subcase: Subcase, dof_map, ndof, dtype='float32') -> Tuple[np.array, Any]:
        """[K] = d{P}/dx"""
        model = self.model

        Kbb = np.zeros((ndof, ndof), dtype=dtype)
        Kbbs = sci_sparse.dok_matrix((ndof, ndof), dtype=dtype)
        #print(dof_map)

        #crods = model._type_to_id_map['CROD']
        #ctubes = model._type_to_id_map['CTUBE']
        #print('celas1s =', celas1s)
        #_get_loadid_ndof(model, subcase_id)
        nelements = 0
        nelements += _build_kbb_celas1(model, Kbb, Kbbs, dof_map)
        nelements += _build_kbb_celas2(model, Kbb, Kbbs, dof_map)
        nelements += _build_kbb_conrod(model, Kbb, Kbbs, dof_map)
        nelements += _build_kbb_crod(model, Kbb, Kbbs, dof_map)
        nelements += _build_kbb_ctube(model, Kbb, Kbbs, dof_map)
        nelements += _build_kbb_cbar(model, Kbb, Kbbs, dof_map)
        assert nelements > 0, nelements
        Kbbs2 = Kbbs.tocsc()
        Kbb2 = Kbbs2.toarray()
        error = np.linalg.norm(Kbb - Kbb2)
        if error > 1e-12:
            self.log.warning(f'error = {error}')
        return Kbb, Kbbs2

    def _recover_strain_101(self, xg, dof_map, fdtype='float32'):
        """recovers the strains"""
        eids = 'ALL'
        nelements = 0
        nelements += _recover_strain_celas1(self.model, dof_map, xg, eids)
        nelements += _recover_strain_celas2(self.model, dof_map, xg, eids)
        #assert nelements > 0, nelements

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
                        spc_set.append(idof)
                        xspc[idof] = 0.
        spc_set = np.array(spc_set, dtype='int32')
        #print('spc_set =', spc_set, xspc)
        return spc_set, xspc


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

        log = self.model.log

        dof_map = _get_dof_map(self.model)
        node_gridtype = _get_node_gridtype(self.model, idtype=idtype)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)
        Kbb, Kbbs = self.build_Kbb(subcase, dof_map, ndof)
        self.get_mpc_constraints(subcase, dof_map)

        Kgg = self.Kbb_to_Kgg(Kbb, ngrid, ndof_per_grid, inplace=False)
        Kggs = self.Kbb_to_Kgg(Kbbs, ngrid, ndof_per_grid)
        del Kbb, Kbbs

        gset = np.arange(ndof, dtype='int32')
        sset, xg = self.build_xg(dof_map, ndof, subcase)
        aset = np.setdiff1d(gset, sset) # a = g-s

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
        fspc = Ksa @ xa + Kss @ xs
        Fg[sset] = fspc

        xb = self.xg_to_xb(xg, ngrid, ndof_per_grid)
        Fb = self.xg_to_xb(Fg, ngrid, ndof_per_grid)

        log.debug(f'Fs = {Fs}')
        page_num = 1
        self._save_spc_forces(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, Fg,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        self._save_displacment(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, xg,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)

        if 'STRAIN' in subcase:
            self._recover_strain_101(xb, dof_map)
        #if 'STRESS' in subcase:
        #Fg[sz_set] = -1
        #xg[sz_set] = -1
        log.debug(f'xa = {xa}')
        return xa_

    def _save_spc_forces(self, f06_file,
                         subcase: Subcase, itime: int, ntimes: int,
                         node_gridtype, Fg,
                         ngrid: int, ndof_per_grid: int,
                         title='', subtitle='', label='',
                         fdtype='float32', page_num=1, page_stamp='PAGE %s') -> int:
        f06_request_name = 'SPCFORCES'
        table_name = 'OQG1'
        #self.log.debug(f'Fg = {Fg}')
        page_num = self._save_static_table(
            f06_file,
            subcase, itime, ntimes,
            node_gridtype, Fg,
            RealSPCForcesArray, f06_request_name, table_name,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            fdtype=fdtype, page_num=page_num, page_stamp=page_stamp)
        return page_num

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
            page_num = spc_forces.write_f06(f06_file, header=None,
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
        fdtype = 'float64'
        log = self.model.log

        dof_map = _get_dof_map(self.model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)
        Kbb, Kbbs = self.build_Kbb(subcase, dof_map, ndof)

        Mbb = np.eye(Kbb.shape[0], dtype=fdtype)

        Kgg = self.Kbb_to_Kgg(Kbb, ngrid, ndof_per_grid)
        Mgg = self.Kbb_to_Kgg(Mbb, ngrid, ndof_per_grid)
        Kggs = self.Kbb_to_Kgg(Kbbs, ngrid, ndof_per_grid)
        del Kbb, Mgg

        gset = np.arange(ndof, dtype='int32')
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
        Kaa_, ipositive, inegative, sz_set = remove_rows(Kaa, aset)
        Kaas_, ipositive, inegative, sz_set = remove_rows(Kaas, aset)
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
        pass

def get_ieids_eids(model: BDF, etype, eids_str, idtype='int32', fdtype='float32'):
    """helper for the stress/strain/displacment recovery"""
    eids = model._type_to_id_map[etype]
    if eids_str == 'ALL':
        ieids = np.arange(len(eids), dtype=idtype)
    else:
        ieids = np.searchsorted(eids_str, eids)
    empty_array = np.full(len(ieids), np.nan, dtype=fdtype)
    return ieids, eids, empty_array

def _recover_strain_celas1(model: BDF, dof_map, xg, eids, fdtype='float32'):
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    ielas, celas1s, strains = get_ieids_eids(model, 'CELAS1', eids, fdtype=fdtype)
    for ieid, eid in zip(ielas, celas1s):
        elem = model.elements[eid]
        strain = _recover_straini_celas12(xg, dof_map, elem)
        strains[ielas] = strain
    return len(celas1s)

def _recover_strain_celas2(model: BDF, dof_map, xg, eids, fdtype='float32'):
    """
    recovers static strain
    TODO: write the OP2/F06

    """
    ielas, celas2s, strains = get_ieids_eids(model, 'CELAS2', eids, fdtype=fdtype)
    #celas3s = model._type_to_id_map['CELAS3']
    #celas4s = model._type_to_id_map['CELAS4']
    for ieid, eid in zip(ielas, celas2s):
        elem = model.elements[eid]
        strain = _recover_straini_celas12(xg, dof_map, elem)
        strains[ielas] = strain
    return len(celas2s)

def _build_kbb_celas1(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CELAS1 Kbb matrix"""
    celas1s = model._type_to_id_map['CELAS1']
    for eid in celas1s:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    return len(celas1s)

def _build_kbb_celas2(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CELAS2 Kbb matrix"""
    celas2s = model._type_to_id_map['CELAS2']
    #celas3s = model._type_to_id_map['CELAS3']
    #celas4s = model._type_to_id_map['CELAS4']
    for eid in celas2s:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    return len(celas2s)

def _recover_straini_celas12(xg, dof_map, elem):
    """get the static strain"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[i] - xg[j]  # TODO: check the sign
    return strain

def _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki):
    """fill the CELASx Kbb matrix"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    k = ki * np.array([[1, -1,],
                       [-1, 1]])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += k[ie1, ie2]
            Kbbs[ib1, ib2] += k[ie1, ie2]
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2

def _build_kbb_cbar(model, Kbb, Kbbs, dof_map):
    """fill the CBAR Kbb matrix"""
    cbars = model._type_to_id_map['CBAR']
    for eid in cbars:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(cbars)

def _build_kbb_crod(model, Kbb, Kbbs, dof_map):
    """fill the CROD Kbb matrix"""
    crods = model._type_to_id_map['CROD']
    for eid in crods:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(crods)

def _build_kbb_ctube(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CTUBE Kbb matrix"""
    ctubes = model._type_to_id_map['CTUBE']
    for eid in ctubes:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(ctubes)

def _build_kbb_conrod(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CONROD Kbb matrix"""
    conrods = model._type_to_id_map['CONROD']
    for eid in conrods:
        elem = model.elements[eid]
        mat = elem.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(conrods)

def _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat, fdtype='float64'):
    """fill the ith rod Kbb matrix"""
    nid1, nid2 = elem.nodes
    #mat = elem.mid_ref
    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    L = np.linalg.norm(dxyz12)
    E = mat.E
    G = mat.G()
    J = elem.J()
    A = elem.Area()
    E = elem.E()
    #L = elem.Length()
    k_axial = A * E / L
    k_torsion = G * J / L

    assert isinstance(k_axial, float), k_axial
    assert isinstance(k_torsion, float), k_torsion
    #Kbb[i, i] += ki[0, 0]
    #Kbb[i, j] += ki[0, 1]
    #Kbb[j, i] = ki[1, 0]
    #Kbb[j, j] = ki[1, 1]
    k = np.array([[1., -1.],
                  [-1., 1.]])  # 1D rod
    Lambda = _lambda1d(dxyz12, debug=False)
    K = Lambda.T @ k @ Lambda
    #i11 = dof_map[(n1, 1)]
    #i12 = dof_map[(n1, 2)]

    #i21 = dof_map[(n2, 1)]
    #i22 = dof_map[(n2, 2)]

    nki, nkj = K.shape
    K2 = np.zeros((nki*2, nkj*2), dtype=fdtype)

    i1 = 0
    i2 = 3 # dof_map[(n1, 2)]
    if k_axial == 0.0 and k_torsion == 0.0:
        dofs = []
        n_ijv = []
        K2 = []
    elif k_torsion == 0.0: # axial; 2D or 3D
        K2 = K * k_axial
        #dofs = np.array([
            #i1, i1+1, i1+2,
            #i2, i2+1, i2+2,
        #], 'int32')
        #n_ijv = [
            ## axial
            #(nid1, 1), (nid1, 2), (nid1, 3),
            #(nid2, 1), (nid2, 2), (nid2, 3),
        #]
    elif k_axial == 0.0: # torsion; assume 3D
        K2 = K * k_torsion
        #dofs = np.array([
            #i1+3, i1+4, i1+5,
            #i2+3, i2+4, i2+5,
        #], 'int32')
        #n_ijv = [
            ## torsion
            #(nid1, 4), (nid1, 5), (nid1, 6),
            #(nid2, 4), (nid2, 5), (nid2, 6),
        #]

    else:  # axial + torsion; assume 3D
        # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
        K2[:nki, :nki] = K * k_axial

        # u1mx, u1my, u1mz, u2mx, u2my, u2mz
        K2[nki:, nki:] = K * k_torsion

        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,

            i1+3, i1+4, i1+5,
            i2+3, i2+4, i2+5,
        ], 'int32')
        n_ijv = [
            # axial
            (nid1, 1), (nid1, 2), (nid1, 3),
            (nid2, 1), (nid2, 2), (nid2, 3),

            # torsion
            (nid1, 4), (nid1, 5), (nid1, 6),
            (nid2, 4), (nid2, 5), (nid2, 6),
        ]
    for dof1, nij1 in zip(dofs, n_ijv):
        i1 = dof_map[nij1]
        for dof2, nij2 in zip(dofs, n_ijv):
            i2 = dof_map[nij2]
            #ki = K2[i1, i2]  #old
            ki = K2[dof1, dof2]  # new
            if abs(ki) > 0.:
                #print(nij1, nij2, f'({i1}, {i2});', (dof1, dof2), ki)
                #Kbb[dof1, dof2] = ki #  old
                Kbb[i1, i2] = ki # new
                Kbbs[i1, i2] = ki
        #print(K2)
    #print(Kbb)
    return

def _lambda1d(dxyz, debug=True):
    """
    ::
      3d  [l,m,n,0,0,0]  2x6
          [0,0,0,l,m,n]
    """
    #R = self.Rmatrix(model,is3D)

    #xyz1 = model.Node(n1).get_position()
    #xyz2 = model.Node(n2).get_position()
    #v1 = xyz2 - xyz1
    if debug:
        print("v1=%s" % dxyz)
    n = np.linalg.norm(dxyz)
    if n == 0:
        raise ZeroDivisionError(dxyz)

    (l, m, n) = dxyz / n
    #l = 1
    #m = 2
    #n = 3
    Lambda = np.zeros((2, 6), 'd')
    Lambda[0, 0] = Lambda[1, 3] = l
    Lambda[0, 1] = Lambda[1, 4] = m
    Lambda[0, 2] = Lambda[1, 5] = n

    #print("R = \n",R)
    #debug = True
    if debug:
        print("Lambda = \n" + str(Lambda))
    return Lambda

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

