from typing import List, Dict, Tuple, Any
import numpy as np

from pyNastran.bdf.bdf import BDF, Subcase
from pyNastran.bdf.mesh_utils.loads import _get_dof_map, _get_loadid_ndof, get_ndof

class Solver:
    def __init__(self, model: BDF):
        self.model = model
        self.log = model.log

    def run(self):
        sol = self.model.sol
        solmap = {
            101 : self.run_sol_101,
            103 : self.run_sol_103,
        }
        self.model.cross_reference()
        self._update_card_count()

        if sol in [101, 103, 105, 107, 109, 111, 112]:
            for subcase_id, subcase in sorted(self.model.subcases.items()):
                if subcase_id == 0:
                    continue
                self.log.debug(f'subcase_id={subcase_id}')
                runner = solmap[sol]
                runner(subcase)
        else:
            raise NotImplementedError(sol)

    def _update_card_count(self):
        for card_type, values in self.model._type_to_id_map.items():
            self.model.card_count[card_type] = len(values)

    def build_Kbb(self, subcase: Subcase) -> np.array:
        model = self.model
        unused_ndof_per_grid, ndof = get_ndof(model, subcase)

        Kbb = np.zeros((ndof, ndof), dtype='float32')
        dof_map = _get_dof_map(model)
        #print(dof_map)

        #celas2s = model._type_to_id_map['CELAS2']
        #celas3s = model._type_to_id_map['CELAS3']
        #celas4s = model._type_to_id_map['CELAS4']

        #crods = model._type_to_id_map['CROD']
        #ctubes = model._type_to_id_map['CTUBE']
        #print('celas1s =', celas1s)
        #_get_loadid_ndof(model, subcase_id)
        _build_kbb_celas1(model, Kbb, dof_map)
        _build_kbb_conrod(model, Kbb, dof_map)

        #print(Kbb)
        return Kbb, dof_map, ndof

    def build_xb(self, dof_map: Dict[Any, int], ndof: int, subcase: Subcase) -> np.ndarray:
        model = self.model
        xspc = np.full(ndof, np.nan, dtype='float64')

        spc_id, unused_options = subcase['SPC']
        if 'SPC' not in subcase:
            model.log.warning(f'no spcs...{spc_id}')
            model.log.warning(str(subcase))
            spc_set = np.array([], dtype='int32')
            return spc_set, xspc
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
                    fi = dof_map[(nid, dof)]  # TODO: wrong...
                    Fb[fi] = fxyz[i]
            else:
                print(load.get_stats())
                raise NotImplementedError(load)
        #print(subcase)
        return Fb

    def Kbb_to_Kgg(self, Kbb: np.ndarray) -> np.ndarray:
        """TODO: transform"""
        Kgg = Kbb
        return Kgg

    def remove_rows(self, Kgg: np.ndarray) -> np.ndarray:
        """
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
        q   Generalized degrees-of-freedom assigned to component modes
            and residual vectors.
        r   reference degrees-of-freedom used to determine free body motion.
        c   DOFs that are free during component mode synthesis or dynamic reduction.
        b   DOFs fixed during component mode analysis or dynamic reduction.
        lm  Lagrange multiplier degrees-of-freedom created by the rigid
            elements using the LAGR method on the Case Control command, RIGID.
        e   Extra degrees-of-freedom introduced in dynamic analysis.
        sa  Permanently constrained aerodynamic degrees-of-freedom.
        k   Aerodynamic mesh point set for forces and displacements on the aero mesh.
        j   Aerodynamic mesh collocation point set (exact physical
            interpretation is dependent on the aerodynamic theory).
        s = sb + sg      all DOFs eliminated by single point constraints
        l = b + c + lm   the DOFs remaining after the reference DOFs are removed (DOF left over)
        t = l + r        the total set of physical boundary DOF for superelements
        a = t + q        the analysis set used in eigensolution
        d = a + e        the set used in dynamic analysis by the direct method
        f = a + o        unconstrained (free) structural DOFs
        fe = f + e       free DOFs plus extra DOFs
        n = f + s        all DOFs not constrained by multipoint constraints
        ne = n + e       all DOFs not constrained by multipoint constraints plus
                         extra degrees-off reedom
        m = mp + mr      all DOFs eliminated by multipoint constraints
        g = n + m        all DOFs including scalar DOFs
        p = g + e        all physical DOFs including extra point DOFs
        ks = k + sa      the union of k and the re-used s-set (6 dof per grid)
        js = j + sa      the union of j and the re-used s-set (6 dof per grid)

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
        ipositive = np.where((col_kgg > 0.) | (row_kgg > 0))[0]
        Kaa = Kgg[ipositive, :][:, ipositive]
        return Kaa, ipositive

    def run_sol_101(self, subcase: Subcase):
        log = self.model.log
        Kbb, dof_map, ndof = self.build_Kbb(subcase)

        gset = np.arange(ndof, dtype='int32')
        sset, xb = self.build_xb(dof_map, ndof, subcase)
        aset = np.setdiff1d(gset, sset) # a = g-s

        Fb = self.build_Fb(dof_map, ndof, subcase)
        Kgg = self.Kbb_to_Kgg(Kbb)
        Fg = Fb
        Fa, Fs = partition_vector(Fb, [['a', aset], ['s', sset]])
        del Kbb, Fb
        # Mgg = self.build_Mgg(subcase)

        # aset - analysis set
        # sset - SPC set
        xa, xs = partition_vector(xb, [['a', aset], ['s', sset]])
        #print(f'xa = {xa}')
        #print(f'xs = {xs}')
        #print(Kgg)
        K = partition_matrix(Kgg, [['a', aset], ['s', sset]])
        Kaa = K['aa']
        Kss = K['ss']
        Kas = K['as']
        Ksa = K['sa']
            #[Kaa]{xa} + [Kas]{xs} = {Fa}
            #[Ksa]{xa} + [Kss]{xs} = {Fs}

        #{xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
        #{Fs} = [Ksa]{xa} + [Kss]{xs}

        # TODO: apply SPCs
        #print(Kaa)
        #print(Kas)
        #print(Kss)
        Kaa_, ipositive = self.remove_rows(Kaa)
        Fs = np.zeros(ndof, dtype='float64')
        #print(f'Fg = {Fg}')
        #print(f'Fa = {Fa}')
        #print(f'Fs = {Fs}')
        Fa_ = Fa[ipositive]
        # [A]{x} = {b}
        # [Kaa]{x} = {F}
        # {x} = [Kaa][F]
        #print(f'Kaa:\n{Kaa}')
        #print(f'Fa: {Fa}')

        #print(f'Kaa_:\n{Kaa_}')
        #print(f'Fa_: {Fa_}')
        xa_ = np.linalg.solve(Kaa_, Fa_)
        #print(f'xa_ = {xa_}')

        xa[ipositive] = xa_
        Fsi = Ksa @ xa + Kss @ xs
        #Fs[ipositive] = Fsi
        log.debug(f'xa = {xa}')
        log.debug(f'Fs = {Fs}')
        return xa_

    def run_sol_103(self, subcase: Subcase):
        raise NotImplementedError(subcase)

def _build_kbb_celas1(model, Kbb, dof_map):
    celas1s = model._type_to_id_map['CELAS1']

    for eid in celas1s:
        elem = model.elements[eid]
        ki = elem.K()

        #print(elem, ki)
        #print(elem.get_stats())
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
        #Kbb[j, i] += ki
        #Kbb[i, j] += ki
        del i, j, ki, nid1, nid2, c1, c2

def _build_kbb_conrod(model, Kbb, dof_map):
    conrods = model._type_to_id_map['CONROD']
    for eid in conrods:
        elem = model.elements[eid]
        nid1, nid2 = elem.nodes
        mat = elem.mid_ref
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
        K2 = np.zeros((nki*2, nkj*2), 'float64')
        if k_axial == 0.0 and k_torsion == 0.0:
            dofs = []
            n_ijv = []
            K2 = []
        elif k_torsion == 0.0: # axial; 2D or 3D
            K2 = K * k_axial
            dofs = np.array([
                i1, i1+1, i1+2,
                i2, i2+1, i2+2,
            ], 'int32')
            n_ijv = [
                # axial
                (nid1, 1), (nid1, 2), (nid1, 3),
                (nid2, 1), (nid2, 2), (nid2, 3),
            ]
        elif k_axial == 0.0: # torsion; assume 3D
            K2 = K * k_torsion
            dofs = np.array([
                i1+3, i1+4, i1+5,
                i2+3, i2+4, i2+5,
            ], 'int32')
            n_ijv = [
                # torsion
                (nid1, 4), (nid1, 5), (nid1, 6),
                (nid2, 4), (nid2, 5), (nid2, 6),
            ]

        else:  # axial + torsion; assume 3D
            # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
            K2[:nki, :nki] = K * k_axial

            # u1mx, u1my, u1mz, u2mx, u2my, u2mz
            K2[nki:, nki:] = K * k_torsion

            i1 = 0
            i2 = 3 # dof_map[(n1, 2)]
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
                    Kbb[dof1, dof2] = K2[i1, i2]
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
    Lambda = np.zeros((2, 6), dtype='float64')
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
