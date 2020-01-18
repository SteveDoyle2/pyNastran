from typing import List, Dict, Any
import numpy as np

from pyNastran.bdf.bdf import BDF, CaseControlDeck, Subcase
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

        celas1s = model._type_to_id_map['CELAS1']
        #celas2s = model._type_to_id_map['CELAS2']
        #celas3s = model._type_to_id_map['CELAS3']
        #celas4s = model._type_to_id_map['CELAS4']

        #crods = model._type_to_id_map['CROD']
        conrods = model._type_to_id_map['CONROD']
        #ctubes = model._type_to_id_map['CTUBE']
        #print('celas1s =', celas1s)
        #_get_loadid_ndof(model, subcase_id)

        for eid in celas1s:
            elem = model.elements[eid]
            ki = elem.K()

            #print(elem, ki)
            #print(elem.get_stats())
            n1, n2 = elem.nodes
            c1, c2 = elem.c1, elem.c2
            i = dof_map[(n1, c1)]
            j = dof_map[(n2, c2)]
            Kbb[j, i] += ki
            Kbb[i, j] += ki
            del i, j, ki, n1, n2, c1, c2

        for eid in conrods:
            elem = model.elements[eid]
            n1, n2 = elem.nodes
            i11 = dof_map[(n1, 1)]
            i12 = dof_map[(n1, 2)]

            i21 = dof_map[(n2, 1)]
            i22 = dof_map[(n2, 2)]
            mat = elem.mid_ref
            n1 = elem.nodes_ref[0].get_position()
            n2 = elem.nodes_ref[1].get_position()
            v = n2 - n1
            L = np.linalg.norm(v)
            E = mat.E
            G = mat.G()
            J = elem.J()
            A = elem.Area()
            E = elem.E()
            #L = elem.Length()
            print(A, E, L)
            ka = A * E / L
            #print(G, J, L)
            kt = G * J / L
            assert isinstance(ka, float), ka
            assert isinstance(kt, float), kt
            #Kbb[i, i] += ki[0, 0]
            #Kbb[i, j] += ki[0, 1]
            #Kbb[j, i] = ki[1, 0]
            #Kbb[j, j] = ki[1, 1]

        #print(Kbb)
        return Kbb, dof_map, ndof

    def build_Fb(self, dof_map: Dict[Any, int], ndof: int, subcase: Subcase) -> np.array:
        model = self.model
        Fb = np.zeros(ndof, dtype='float32')
        if 'LOAD' in subcase:
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
        ipositive = np.where((col_kgg > 0.) | (row_kgg > 0))[0]
        Kaa = Kgg[ipositive, :][:, ipositive]
        return Kaa, ipositive

    def run_sol_101(self, subcase: Subcase):
        Kbb, dof_map, ndof = self.build_Kbb(subcase)
        Fb = self.build_Fb(dof_map, ndof, subcase)
        Kgg = self.Kbb_to_Kgg(Kbb)
        # Mgg = self.build_Mgg(subcase)

        # TODO: apply SPCs
        Kaa, ipositive = self.remove_rows(Kgg)
        Fg = Fb[ipositive]
        # [A]{x} = {b}
        # [Kaa]{x} = {F}
        # {x} = [Kaa][F]
        x = np.linalg.solve(Kaa, Fg)
        #print(x)
        return x

    def run_sol_103(self, subcase: Subcase):
        raise NotImplementedError(subcase)
