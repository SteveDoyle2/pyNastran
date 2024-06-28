from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, SUPORT

class SetMethods:
    def __init__(self, model: BDF):
        self.model = model

    def get_sset(self, suport_id: int=0) -> np.ndarray:
        model = self.model
        nid_dof_list = []
        if self.model.grid.has_ps():
            ps = self.model.grid.ps
            ips = np.where(ps > 0)[0]
            if len(ips):
                nids = self.model.grid.node_id[ips]
                dof = ps[ips]
                nid_dofi = np.column_stack([nids, dof])
                nid_dof_list.append(nid_dofi)

        nid_dof = _nid_dof(nid_dof_list)
        return nid_dof

    def get_rset(self, suport_id: int=0) -> np.ndarray:
        model = self.model
        assert suport_id >= 0, model.suport.get_stats()
        suport: SUPORT = model.suport

        suport_ids = suport.suport_id
        nids = suport.node_id
        dofs = suport.component
        assert len(nids) == len(suport_ids)
        assert len(nids) == len(dofs)

        nid_dof_list = []
        i0 = np.where(suport_ids == 0)[0]
        if len(i0):
            nid_dofi = np.column_stack([nids[i0], dofs[i0]])
            nid_dof_list.append(nid_dofi)

        if suport_id > 0:
            i1 = np.where(suport_ids == suport_id)[0]
            assert len(i1) > 0, f'cant find SUPORT1={suport_id:d}'
            nid_dofi = np.column_stack([nids[i1], dofs[i1]])
            nid_dof_list.append(nid_dofi)

        nid_dof = _nid_dof(nid_dof_list)
        return nid_dof

    def get_mset(self, mpc_id: int=0, include_rbe: bool=True) -> np.ndarray:
        """gets the dependent degrees of freedom"""
        model = self.model
        log = model.log
        rigid_cards = [card for card in model.rigid_element_cards if len(card)]
        nid_dof = np.zeros((0, 2), dtype='int32')
        if mpc_id == 0 and len(rigid_cards) == 0:
            return nid_dof

        nid_dof_list = []
        if include_rbe:
            for elem in rigid_cards:
                if elem.type == 'RBE2':
                    for dof, (i0, i1) in zip(elem.independent_dof, elem.idim):
                        dep_nodes = elem.dependent_nodes[i0:i1]
                        dofs = np.ones(len(dep_nodes), dtype='int32') * dof
                        nid_dofi = np.column_stack([dep_nodes, dofs])
                        nid_dof_list.append(nid_dofi)
                elif elem.type == 'RBE3':
                    print(elem.get_stats())
                    rbe3
                elif elem.type in {'RBAR', 'RBAR1', 'RROD', 'RBE1'}:
                    log.warning(f'skipping {elem.type}...')
                else:  # pragma: no cover
                    raise RuntimeError(elem.type)
            if len(rigid_cards) == 0:
                assert len(nid_dof_list) > 0, nid_dof_list

        if mpc_id in model.mpcadd.mpc_id:
            aaa

        mpc_cards = [card for card in model.mpc_cards if len(card)]
        for mpc in mpc_cards:
            if mpc_id in mpc.mpc_id:
                card = mpc.slice_card_by_id(mpc_id)
                if card.type == 'MPC':
                    #print(card.get_stats())
                    assert len(card.components) == len(card.node_id), card.get_stats()
                    nid_dofi = np.column_stack([card.node_id, card.components])
                else:
                    print(card.get_stats())
                    raise RuntimeError(card.type)
                nid_dof_list.append(nid_dofi)

        nid_dof = _nid_dof(nid_dof_list)
        return nid_dof

    def get_spline_nodes(self) -> np.ndarray:
        model = self.model
        log = model.log
        set_ids_type = _get_setx_type(model)
        set_ids = set_ids_type[:, 0]
        set_types = set_ids_type[:, 1]

        spline_nodes_list = []
        spline_cards = [card for card in model.aero_spline_cards if len(card)]
        for spline in spline_cards:
            if spline.type in {'SPLINE1', 'SPLINE2', 'SPLINE4'}:
                missing_sets = np.setdiff1d(spline.set_id, set_ids)
                assert len(missing_sets) == 0, missing_sets
                iset = np.searchsorted(set_ids, spline.set_id)
                set_type = set_types[iset]
                _get_spline_set_nodes(model, spline_nodes_list, set_ids, set_type)
            else:  # pragma: no cover
                raise RuntimeError(spline.get_stats())
        assert len(spline_nodes_list) > 0, spline_nodes_list
        node_ids = np.unique(np.hstack(spline_nodes_list))
        return node_ids

def _get_spline_set_nodes(model: BDF,
                          spline_nodes_list: list[np.ndarray],
                          set_ids: np.ndarray,
                          set_type: np.ndarray) -> None:
    for set_idi, set_typei in zip(set_ids, set_type):
        if set_typei == 1:
            set_card = model.set1
            card = set_card.slice_card_by_id(set_idi)
            ids = card.ids
            spline_nodes_list.append(ids)
        else:  # pragma: no cover
            raise RuntimeError(f'cant find SET{set_typei}')


def _get_setx_type(model: BDF) -> np.ndarray:
    """build a data structure to access the SET1-4 cards"""
    set_ids_list = []
    for set in model.setx_cards:
        nset = len(set)
        if nset == 0:
            continue
        iset = int(set.type[-1])
        isets = np.ones(nset) * iset
        # print(set.get_stats())
        set_idsi = np.column_stack([set.set_id, isets])
        set_ids_list.append(set_idsi)
    assert len(set_ids_list) > 0, set_ids_list
    set_ids_type = np.vstack(set_ids_list)
    set_ids = np.unique(set_ids_type[:, 0])
    isort = np.argsort(set_ids_type[:, 0])
    set_ids_type = set_ids_type[isort, :]
    assert len(set_ids) == len(set_ids_type), 'duplicate set_ids'
    return set_ids_type


def _nid_dof(nid_dof_list: list[np.ndarray]) -> np.ndarray:
    if len(nid_dof_list):
        nid_dof = np.vstack(nid_dof_list)
    else:
        nid_dof = np.zeros((0, 2), dtype='int32')
    return nid_dof
