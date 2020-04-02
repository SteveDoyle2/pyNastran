"""
defines various methods to access high level BDF data:
 - GetCard()
   - get_card_ids_by_card_types(self, card_types=None, reset_type_to_slot_map=False,
                                stop_on_missing_card=False, combine=False)
   - get_rslot_map(self, reset_type_to_slot_map=False)
   - get_cards_by_card_types(self, card_types, reset_type_to_slot_map=False,
                             stop_on_missing_card=False)
   - get_SPCx_node_ids(self, spc_id, stop_on_failure=True)
   - get_SPCx_node_ids_c1( spc_id, stop_on_failure=True)
   - get_reduced_loads(self, load_id, scale=1., skip_scale_factor0=True, msg='')
   - get_reduced_dloads(self, dload_id, scale=1., skip_scale_factor0=True, msg='')
   - get_node_ids_with_elements(self, eids, msg='')
   - get_elements_nodes_by_property_type(self, dtype='int32',
                                         save_element_types=False)
   - get_elements_properties_nodes_by_element_type(self, dtype='int32', solids=None)
   - get_element_ids_list_with_pids(self, pids=None)
   - get_pid_to_node_ids_and_elements_array(self, pids=None, etypes=None, idtype='int32')
   - get_element_ids_dict_with_pids(self, pids=None, stop_if_no_eids=True)
   - get_node_id_to_element_ids_map(self)
   - get_node_id_to_elements_map(self)
   - get_property_id_to_element_ids_map(self)
   - get_material_id_to_property_ids_map(self)
   - get_reduced_mpcs(self, mpc_id)
   - get_reduced_spcs(self, spc_id)
   - get_spcs(self, spc_id, consider_nodes=False)

"""
# pylint: disable=C0103
from copy import deepcopy
from collections import defaultdict
from typing import List, Dict, Set, Tuple, Optional, Union, Any

import numpy as np

from pyNastran.bdf.bdf_interface.get_methods import GetMethods
from pyNastran.utils.numpy_utils import integer_types

from pyNastran.bdf.mesh_utils.dvxrel import get_dvprel_ndarrays
#from pyNastran.bdf.mesh_utils.forces_moments import (
    # get_load_arrays, get_forces_moments_array,
    #get_pressure_array, get_temperatures_array,
#)
from pyNastran.bdf.mesh_utils.mpc_dependency import (
    get_mpc_node_ids, get_mpc_node_ids_c1,
    get_rigid_elements_with_node_ids, get_dependent_nid_to_components,
    get_lines_rigid, get_mpcs)


class GetCard(GetMethods):
    """defines various methods to access high level BDF data"""
    def __init__(self) -> None:
        self._type_to_slot_map = {}
        GetMethods.__init__(self)

    def get_card_ids_by_card_types(self, card_types: List[str]=None,
                                   reset_type_to_slot_map: bool=False,
                                   stop_on_missing_card: bool=False,
                                   combine: bool=False) -> Dict[str, List[int]]:
        """
        Parameters
        ----------
        card_types : str / List[str] / default=None
            the list of keys to consider (list of strings; string)
            None : all cards
        reset_type_to_slot_map : bool
            should the mapping dictionary be rebuilt (default=False);
            set to True if you added cards
        stop_on_missing_card : bool
            crashes if you request a card and it doesn't exist
        combine : bool; default=False
            change out_dict into out_list
            combine the list of cards

        Returns
        -------
        out_dict: dict[str]=List[ids]
            the key=card_type, value=the ID of the card object
        out_list: List[ids]
            value=the ID of the card object
            useful

        Examples
        ---------
        >>> out_dict = model.get_card_ids_by_card_types(
            card_types=['GRID', 'CTRIA3', 'CQUAD4'], combine=False)
        >>> out_dict = {
            'GRID' : [1, 2, 10, 42, 1000],
            'CTRIA3' : [1, 2, 3, 5],
            'CQUAD4' : [4],
        }

        **shell elements**

          >>> out_dict = model.get_card_ids_by_card_types(
              card_types=['CTRIA3', 'CQUAD4'], combine=True)
          >>> out_dict = {
              [1, 2, 3, 4, 5],
          }

        """
        if card_types is None:
            card_types = list(self.cards_to_read)

        if isinstance(card_types, str):
            card_types = [card_types]
        elif not isinstance(card_types, (list, tuple)):
            raise TypeError('card_types must be a list/tuple; type=%s' % type(card_types))

        #if reset_type_to_slot_map or self._type_to_slot_map is None:
            #self._type_to_slot_map = rslot_map
        if reset_type_to_slot_map:
            self._reset_type_to_slot_map()
        #out_dict = {
            #(key) : (self._type_to_id_map[key] if key in self.card_count else [])
            #for key in card_types
        #}
        out_dict = {}
        for key in card_types:
            if key in self.card_count:
                out_dict[key] = sorted(self._type_to_id_map[key])
            else:
                if stop_on_missing_card:
                    raise RuntimeError('%r is not in the card_count; keys=%s' %
                                       str(sorted(self.card_count.keys())))
                out_dict[key] = []
        if combine:
            out_list = []
            for key, value in sorted(out_dict.items()):
                out_list += value
            return out_list
        return out_dict

    def _reset_type_to_slot_map(self) -> Dict[str, str]:
        """resets self._type_to_slot_map"""
        rslot_map = defaultdict(list)
        for dict_name, card_names in self._slot_to_type_map.items():
            #print('card_names=%s dict_name=%s' % (card_names, dict_name))
            card_name0 = card_names[0]
            if card_name0 in ['DTABLE', 'GRDSET', 'SESUP', 'DOPTPRM', 'MONPNT1', 'SUPORT',
                              'MKAERO1', 'MATHP']:
                pass
            else:
                adict = getattr(self, dict_name)
                if isinstance(adict, dict):
                    for key, card in adict.items():
                        if isinstance(card, list):
                            alist = card
                            for cardi in alist:
                                rslot_map[cardi.type].append(key)
                            #msg = '%s; names=%s \ncard=%s' % (type(card), card_names, card)
                            #raise NotImplementedError(msg)
                        else:
                            rslot_map[card.type].append(key)
                elif isinstance(adict, list):
                    alist = adict
                    for value in alist:
                        if isinstance(value, list):
                            msg = '%s; names=%s value=%s' % (type(value), card_names, value)
                            raise NotImplementedError(msg)

                        if value.type in ['CSET1', 'CSET']:
                            pass
                            #rslot_map[value.type] = value.
                        else:
                            raise NotImplementedError('list; names=%s' % card_names)
                else:
                    raise NotImplementedError('%s; names=%s' % (type(adict), card_names))
        return rslot_map

    def get_rslot_map(self, reset_type_to_slot_map=False) -> Dict[str, str]:
        """gets the rslot_map"""
        if (reset_type_to_slot_map or self._type_to_slot_map is None or
                len(self._type_to_slot_map) == 0):
            self.reset_rslot_map()

        rslot_map = self._type_to_slot_map
        assert 'GRID' in rslot_map
        return rslot_map

    def reset_rslot_map(self) -> None:
        """helper method for get_rslot_map"""
        rslot_map = {}
        for key, values in self._slot_to_type_map.items():
            for value in values:
                rslot_map[value] = key
        self._type_to_slot_map = rslot_map

    @property
    def nid_map(self) -> Dict[int, int]:
        """
        Gets the GRID/SPOINT/EPOINT ids to a sorted order.

        Parameters
        ----------
        sort_ids : bool; default=True
            sort the ids

        Returns
        -------
        nid_map : Dict[nid] : i
            nid : int
                the GRID/SPOINT/EPOINT id
            i : int
                the index

        ..note ::  GRIDs, SPOINTs, & EPOINTs are stored in separate slots,
                   so they are unorganized.
        ..note :: see ``self.get_nid_map(sort_ids=False)`` for the unsorted version

        """
        return self.get_nid_map(sort_ids=True)

    def get_nid_map(self, sort_ids: bool=True) -> Dict[int, int]:
        """
        Maps the GRID/SPOINT/EPOINT ids to a sorted/unsorted order.

        Parameters
        ----------
        sort_ids : bool; default=True
            sort the ids

        Returns
        -------
        nid_map : Dict[nid] : i
            nid : int
                the GRID/SPOINT/EPOINT id
            i : int
                the index

        ..note ::  GRIDs, SPOINTs, & EPOINTs are stored in separate slots,
                   so they are unorganized.

        """
        nids = []
        index_nids = []
        i = 0
        for nid in self.nodes:
            nids.append(nid)
            index_nids.append(i)
            i += 1
        for nid in self.spoints:
            nids.append(nid)
            index_nids.append(i)
            i += 1
        for nid in self.epoints:
            nids.append(nid)
            index_nids.append(i)
            i += 1

        if sort_ids:
            inids = np.argsort(nids)
            nids = np.sort(nids)
            index_nids = np.array(index_nids)[inids]

        nid_map = {}
        for nid, i in zip(nids, index_nids):
            nid_map[nid] = i
        return nid_map

    def get_cards_by_card_types(self, card_types: List[str],
                                reset_type_to_slot_map: bool=False,
                                stop_on_missing_card: bool=False) -> None:
        """
        Parameters
        ----------
        card_types : List[str]
            the list of keys to consider
        reset_type_to_slot_map : bool
            should the mapping dictionary be rebuilt (default=False);
            set to True if you added cards
        stop_on_missing_card : bool
            crashes if you request a card and it doesn't exist

        Returns
        -------
        out_dict : dict[str] = List[BDFCard()]
            the key=card_type, value=the card object

        """
        if not isinstance(card_types, (list, tuple)):
            raise TypeError(f'card_types={card_types!r} must be a list/tuple; '
                            f'type={type(card_types)}')

        #self._type_to_id_map = {
        #    'CQUAD4' : [1, 2, 3]
        #}
        #self._slot_to_type_map = {'elements' : [CQUAD4, CTRIA3]}
        rslot_map = self.get_rslot_map(reset_type_to_slot_map=False)

        out = {}
        for card_type in card_types:
            if card_type not in self.card_count:
                if stop_on_missing_card:
                    keys = str(sorted(self.card_count.keys()))
                    raise RuntimeError(f'{card_type} is not in the card_count; keys={keys}')
                out[card_type] = []
                continue

            #print('card_type=%r' % card_type)
            try:
                key = rslot_map[card_type]  # update attributes.py ~line 740
            except:
                print(rslot_map.keys())
                self.log.error("card_type=%r' hasn't been added to "
                               "self._slot_to_type_map...check for typos")
                raise
            try:
                slot = getattr(self, key)
            except AttributeError:
                if hasattr(self.zona, key):
                    slot = getattr(self.zona, key)
                else:
                    raise
            ids = self._type_to_id_map[card_type]
            cards = []
            if isinstance(ids, bool):
                continue

            for idi in ids:
                try:
                    card = slot[idi]
                except KeyError:
                    print(slot)
                    msg = 'key=%r id=%r cannot be found\n' % (key, idi)
                    msg += 'id=%s not found.  Allowed=%s' % (
                        key, np.unique(ids))
                    #print(msg)
                    raise KeyError(msg)
                except TypeError:
                    msg = 'key=%s id=%s cannot be found' % (key, idi)
                    #print(msg)
                    raise TypeError(msg)

                if isinstance(card, list):
                    for cardi in card:  # loads/spc/mpc
                        if cardi.type == card_type: # loads
                            cards.append(cardi)
                else:
                    cards.append(card)
                #for card in cards:
                    #print('%s' % str(card).split('\n')[0])
            out[card_type] = cards
        return out

    def get_SPCx_node_ids(self, spc_id: int,
                          consider_spcadd: bool=True,
                          stop_on_failure: bool=True) -> List[int]:
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        ----------
        spc_id : int
            the SPC id
        stop_on_failure : bool; default=True
            errors if parsing something new

        Returns
        -------
        node_ids : List[int]
            the constrained associated node ids

        """
        spcs = self.get_reduced_spcs(
            spc_id, consider_spcadd=consider_spcadd, stop_on_failure=stop_on_failure)

        warnings = ''
        node_ids = []
        for card in spcs:
            if card.type == 'SPC':
                nids = card.node_ids
            elif card.type == 'SPC1':
                nids = card.node_ids
            elif card.type in ['GMSPC', 'SPCAX']:
                warnings += str(card)
                continue
            else:
                warnings += str(card)
                continue
            node_ids += nids
        if warnings:
            self.log.warning("get_SPCx_node_ids doesn't consider:\n%s" % warnings.rstrip('\n'))
        return node_ids

    def get_SPCx_node_ids_c1(self, spc_id: int, stop_on_failure: bool=True) -> Dict[str, List[int]]:
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        ----------
        spc_id : int
            the SPC id
        stop_on_failure : bool; default=True
            errors if parsing something new

        Returns
        -------
        node_ids_c1 : Dict[component] = node_ids
            component : str
                the DOF to constrain
            node_ids : List[int]
                the constrained node ids

        """
        spcs = self.get_reduced_spcs(spc_id, stop_on_failure=stop_on_failure)

        node_ids_c1 = defaultdict(str)
        #print('spcs = ', spcs)
        warnings = ''
        for card in spcs:
            if card.type == 'SPC':
                for nid, c1 in zip(card.node_ids, card.components):
                    assert nid is not None, card.node_ids
                    node_ids_c1[nid] += c1
            elif card.type == 'SPC1':
                nids = card.node_ids
                c1 = card.components
                for nid in nids:
                    node_ids_c1[nid] += c1
            elif card.type in ['GMSPC', 'SPCAX']:
                warnings += str(card)
            else:
                msg = 'get_SPCx_node_ids_c1 doesnt supprt %r' % card.type
                if stop_on_failure:
                    raise RuntimeError(msg)
                self.log.warning(msg)

        if warnings:
            self.log.warning("get_SPCx_node_ids_c1 doesn't consider:\n%s" % warnings.rstrip('\n'))
        return node_ids_c1

    def get_MPCx_node_ids(self, mpc_id: int,
                          consider_mpcadd: bool=True,
                          stop_on_failure: bool=True) -> List[List[int]]:
        """see ``pyNastran.bdf.mesh_utils.mpc_dependency.get_mpc_node_ids(...)``"""
        lines = get_mpc_node_ids(
            self, mpc_id,
            consider_mpcadd=consider_mpcadd,
            stop_on_failure=stop_on_failure)
        return lines

    def get_MPCx_node_ids_c1(self, mpc_id: int,
                             consider_mpcadd: bool=True,
                             stop_on_failure: bool=True) -> Tuple[Dict[str, List[int]],
                                                                  Dict[str, List[int]]]:
        """see ``pyNastran.bdf.mesh_utils.mpc_dependency.get_mpc_node_ids_c1(...)``"""
        independent_node_ids_c1, dependent_node_ids_c1 = get_mpc_node_ids_c1(
            self, mpc_id,
            consider_mpcadd=consider_mpcadd,
            stop_on_failure=stop_on_failure)
        return independent_node_ids_c1, dependent_node_ids_c1

    def get_mpcs(self, mpc_id: int, consider_mpcadd: bool=True,
                 stop_on_failure: bool=True) -> Tuple[List[int], List[str]]:
        """see ``pyNastran.bdf.mesh_utils.mpc_dependency.get_mpcs(...)``"""
        nids, comps = get_mpcs(self, mpc_id, consider_mpcadd=consider_mpcadd,
                 stop_on_failure=stop_on_failure)
        return nids, comps

    def get_rigid_elements_with_node_ids(self, node_ids):
        """see ``pyNastran.bdf.mesh_utils.mpc_dependency.get_rigid_elements_with_node_ids(...)``"""
        rbes = get_rigid_elements_with_node_ids(self, node_ids)
        return rbes

    def get_dependent_nid_to_components(self, mpc_id=None, stop_on_failure=True):
        """see ``pyNastran.bdf.mesh_utils.mpc_dependency.get_dependent_nid_to_components(...)``"""
        dependent_nid_to_components = get_dependent_nid_to_components(
            self, mpc_id=mpc_id, stop_on_failure=stop_on_failure)
        return dependent_nid_to_components

    def _get_rigid(self) -> Any:
        """see ``pyNastran.bdf.mesh_utils.mpc_dependency.get_lines_rigid(...)``"""
        lines_rigid = get_lines_rigid(self)
        return lines_rigid

    def _get_dvprel_ndarrays(self, nelements: int, pids: np.ndarray,
                             fdtype='float32', idtype='int32'):
        """see ``pyNastran.bdf.mesh_utils.dvxrel.get_dvprel_ndarrays(...)``"""
        dvprel_dict = get_dvprel_ndarrays(
            self, nelements, pids, fdtype=fdtype, idtype=idtype)
        return dvprel_dict

    #def get_load_arrays(self, subcase_id, eid_map, node_ids, normals, nid_map=None,
                        #stop_on_failure=True):
        #"""see ``pyNastran.bdf.mesh_utils.forces_moments.get_load_arrays(...)``"""
        #out = get_load_arrays(self, subcase_id, eid_map, node_ids, normals,
                              #nid_map=nid_map, stop_on_failure=stop_on_failure)
        #return out

    #def _get_forces_moments_array(self, p0, load_case_id,
                                  #eid_map, nnodes, normals, dependents_nodes,
                                  #nid_map=None, include_grav=False):
        #"""see ``pyNastran.bdf.mesh_utils.forces_moments.get_forces_moments_array(...)``"""
        #centroidal_pressures, forces, spcd = get_forces_moments_array(
            #self, p0, load_case_id,
            #eid_map, nnodes, normals, dependents_nodes,
            #nid_map=nid_map, include_grav=include_grav)
        #return centroidal_pressures, forces, spcd

    #def get_pressure_array(self, load_case_id, eids, stop_on_failure=True):
        #"""see ``pyNastran.bdf.mesh_utils.forces_moments.get_pressure_array(...)``"""
        #return get_pressure_array(
            #self, load_case_id, eids, stop_on_failure=stop_on_failure)

    #def _get_temperatures_array(self, load_case_id, nid_map=None, dtype='float32'):
        #"""see ``pyNastran.bdf.mesh_utils.forces_moments.get_temperatures_array(...)``"""
        #return get_temperatures_array(self, load_case_id, nid_map=nid_map, dtype=dtype)

    def get_reduced_loads(self, load_case_id, scale=1.,
                          consider_load_combinations=True,
                          skip_scale_factor0=False,
                          stop_on_failure=True, msg=''):
        """
        Accounts for scale factors.

        Parameters
        ----------
        load_case_id : int
            the desired LOAD id
        consider_load_combinations : bool; default=True
            look at the LOAD card
        scale : float; default=1.0
            additional scale factor on top of the existing LOADs
        skip_scale_factor0 : bool; default=False
            Skip loads with scale factor=0.0.
            Nastran does not do this.
            Nastran will fail if referenced loads do not exist.
        stop_on_failure : bool; default=True
            errors if parsing something new
        msg : str
            debug message

        Returns
        -------
        loads : List[loads]
            a series of load objects
        scale_factors : List[float]
            the associated scale factors
        is_grav : bool
            is there a gravity card

        .. warning:: assumes xref=True

        """
        if not isinstance(load_case_id, integer_types):
            msg = 'load_case_id must be an integer; type=%s, load_case_id:\n%r' % (
                type(load_case_id), load_case_id)
            raise TypeError(msg)

        try:
            load_case = self.Load(
                load_case_id, consider_load_combinations=consider_load_combinations, msg=msg)
        except KeyError:
            if stop_on_failure:
                raise
            self.log.error("could not find expected LOAD/LOADSET id=%s" % load_case_id)
            return []

        loads, scale_factors, is_grav = self._reduce_load_case(load_case, scale=scale)
        assert len(loads) == len(scale_factors)
        return loads, scale_factors, is_grav

    def _reduce_load_case(self, load_case, scale=1., consider_load_combinations=True,
                          unallowed_load_ids=None, msg=''):
        """reduces a load case"""
        scale_factors_out = []
        loads_out = []
        is_grav_out = False
        if unallowed_load_ids is None:
            unallowed_load_ids = []

        for load in load_case:
            if load.type == 'LOAD':
                load_ids = load.get_load_ids()
                load_scale = load.scale * scale
                scale_factors = load.scale_factors
                assert len(load_ids) == len(scale_factors), str(load)
                scale_factors_temp = [load_scale * scalei for scalei in scale_factors]
                for load_idi, scalei in zip(load_ids, scale_factors_temp):
                    # prevents recursion
                    if load_idi in unallowed_load_ids:
                        msg = 'There is a recursion error.  LOAD trace=%s; load_id=%s' % (
                            unallowed_load_ids, load_idi)
                        raise RuntimeError(msg)
                    unallowed_load_ids2 = deepcopy(unallowed_load_ids)
                    unallowed_load_ids2.append(load_idi)

                    load_casei = self.Load(
                        load_idi, consider_load_combinations=consider_load_combinations, msg=msg)
                    loadsi, scale_factorsi, is_gravi = self._reduce_load_case(
                        load_casei, scale=scalei,
                        consider_load_combinations=consider_load_combinations,
                        unallowed_load_ids=unallowed_load_ids2)
                    if is_gravi:
                        is_grav_out = True
                    scale_factors_out += scale_factorsi
                    loads_out += loadsi
            elif load.type in 'GRAV':
                scale_factors_out.append(scale)
                loads_out.append(load)
                is_grav_out = True
            else:
                scale_factors_out.append(scale)
                loads_out.append(load)
        return loads_out, scale_factors_out, is_grav_out

    def get_reduced_dloads(self, dload_id, scale=1., consider_dload_combinations=True,
                           skip_scale_factor0=False, msg=''):
        """
        Accounts for scale factors.

        Parameters
        ----------
        dload_id : int
            the desired DLOAD id
        consider_dload_combinations : bool; default=True
            look at the DLOAD card
        scale : float; default=1.0
            additional scale factor on top of the existing LOADs
        skip_scale_factor0 : bool; default=False
            Skip loads with scale factor=0.0.
            Nastran does not do this.
            Nastran will fail if referenced loads do not exist.
        msg : str
            debug message

        Returns
        -------
        dloads : List[loads]
            a series of dload objects
        scale_factors : List[float]
            the associated scale factors

        .. warning:: assumes xref=True

        """
        dload_case = self.DLoad(
            dload_id,
            consider_dload_combinations=consider_dload_combinations,
            msg=msg)
        dloads, scale_factors = self._reduce_dload_case(
            dload_case, scale=scale, skip_scale_factor0=skip_scale_factor0,
            msg=msg)
        return dloads, scale_factors

    def _reduce_dload_case(self, dload_case, scale=1., unallowed_dload_ids=None,
                           skip_scale_factor0=False, msg=''):
        """
        Reduces a dload case

        Parameters
        ----------
        dload_case : List[???]
            a series of DLOAD cards
        scale : float; default=1.0
            additional scale factor on top of the existing LOADs
        unallowed_dload_ids : List[int]; default=None
            helper to prevent recursion
        skip_scale_factor0 : bool; default=False
            Skip loads with scale factor=0.0.
            Nastran does not do this.
            Nastran will fail if referenced loads do not exist.
        msg : str
            debug message

        Returns
        -------
        dloads : List[loads]
            a series of dload objects
        scale_factors : List[float]
            the associated scale factors

        """
        scale_factors_out = []
        dloads_out = []
        if unallowed_dload_ids is None:
            unallowed_dload_ids = []

        assert isinstance(dload_case, list), dload_case
        for dload in dload_case:
            if dload.type == 'DLOAD':
                dload_ids = dload.get_load_ids()
                load_scale = dload.scale * scale
                scale_factors = dload.scale_factors
                if len(dload_ids) != len(scale_factors):
                    msg = 'dload_ids=%s scale_factors=%s\n%s' % (
                        dload_ids, scale_factors, str(dload))
                    raise ValueError(msg)

                scale_factors_temp = [load_scale * scalei for scalei in scale_factors]
                for dload_idi, scalei in zip(dload_ids, scale_factors_temp):
                    # prevents recursion
                    if dload_idi in unallowed_dload_ids:
                        msg = 'There is a recursion error.  DLOAD trace=%s; dload_id=%s' % (
                            unallowed_dload_ids, dload_idi)
                        raise RuntimeError(msg)
                    unallowed_dload_ids2 = deepcopy(unallowed_dload_ids)
                    unallowed_dload_ids2.append(dload_idi)

                    dload_casei = self.DLoad(dload_idi, msg=msg)
                    dloadsi, scale_factorsi = self._reduce_dload_case(
                        dload_casei, scale=scalei, unallowed_dload_ids=unallowed_dload_ids2, )
                    scale_factors_out += scale_factorsi
                    dloads_out += dloadsi
            else:
                scale_factors_out.append(scale)
                dloads_out.append(dload)
        return dloads_out, scale_factors_out

    def _get_maps(self, eids: Optional[List[int]]=None,
                  map_names: Optional[List[str]]=None,
                  consider_0d: bool=True,
                  consider_0d_rigid: bool=True,
                  consider_1d: bool=True,
                  consider_2d: bool=True,
                  consider_3d: bool=True) -> Any:
        """
        Gets a series of mappings (e.g. node_id to element_id)

        eids : List[int]
            the element ids to consider
        map_names : List[str]; default=None -> all
            'edge_to_eid_map', 'eid_to_edge_map', 'nid_to_edge_map', 'nid_to_eid_map'
        consider_0d : bool; default=True
            considers CELASx, CDAMPx, CFAST
        consider_0d_rigid : bool; default=True
            considers MPC, RBAR, RBE2, RBE3, RSPLINE elements
        consider_1d : bool; default=True
            considers CONROD, CROD, CBAR, CBEAM elements
        consider_2d : bool; default=True
            considers CQUAD4, CQUAD8, CQUADR, CQUAD,
            CTRIA3, CTRIA6, CTRIAX, CTRIAX6, CSHEAR elements
        consider_3d : bool; default=True
            considers CTETRA, CPENTA, CPYRAM, CHEXA elements

        .. todo:: consider_0d support
        .. todo:: consider_0d_rigid support

        """
        allowed_maps = [
            'edge_to_eid_map',
            'eid_to_edge_map',
            'nid_to_edge_map',
            #'edge_to_nid_map',  # unnecessary
            #'eid_to_eid_map',  # not added yet
            'nid_to_eid_map',
            #'face_to_edge_map',  # what does this look like?
        ]
        if map_names is None:
            map_names = allowed_maps
        else:
            if isinstance(map_names, str):
                map_names = [map_names]
            if not isinstance(map_names, (list, tuple)):
                msg = 'map_names=%s must be a list or tuple; not %s' % (
                    map_names, type(map_names))
                raise TypeError(msg)
            for name in map_names:
                if name not in allowed_maps:
                    msg = 'name=%r; allowed=%s' % (name, sorted(allowed_maps.keys()))
                    raise RuntimeError(msg)

        eid_to_edge_map = {}
        eid_to_nid_map = {}

        edge_to_eid_map = defaultdict(set)
        nid_to_edge_map = defaultdict(set)  #set() ???
        nid_to_eid_map = defaultdict(set)

        if eids is None:
            eids = self.elements.keys()

        types_to_consider = []
        if consider_0d:
            types_to_consider += []
        if consider_0d_rigid:
            types_to_consider += []
        if consider_1d:
            types_to_consider += ['CROD', 'CONROD', 'CBAR', 'CBEAM', 'CBEAM3']
        if consider_2d:
            types_to_consider += ['CTRIA3', 'CTRIAX', 'CTRIA6', 'CTRIAX6',
                                  'CQUAD4', 'CQUAD', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUADX8',
                                  'CSHEAR']
        if consider_3d:
            types_to_consider += ['CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA']

        for eid in eids:
            elem = self.elements[eid]
            if elem.type not in types_to_consider:
                continue
            node_ids = elem.node_ids
            edges = elem.get_edge_ids()
            eid_to_edge_map[eid] = edges
            eid_to_nid_map[eid] = node_ids

            for nid in node_ids:
                nid_to_eid_map[nid].add(eid)
            for edge in edges:
                assert not isinstance(edge, integer_types), 'edge=%s elem=\n%s' % (edge, elem)
                assert edge[0] < edge[1], 'edge=%s elem=\n%s' % (edge, elem)
                try:
                    edge_to_eid_map[edge].add(eid)
                except TypeError:
                    print(elem)
                    raise
                for nid in edge:
                    nid_to_edge_map[nid].add(tuple(edge))

        out = {}
        allowed_maps = [
            #'edge_to_eid_map',
            #'eid_to_edge_map',
            #'nid_to_edge_map',
            #'edge_to_nid_map', # unnecessary
            #'nid_to_eid_map',
        ]

        for key in map_names:
            if key == 'edge_to_eid_map':
                out[key] = edge_to_eid_map
            elif key == 'eid_to_edge_map':
                out[key] = eid_to_edge_map
            elif key == 'nid_to_edge_map':
                out[key] = nid_to_edge_map
            elif key == 'nid_to_eid_map':
                out[key] = nid_to_eid_map
            #elif key == 'eid_to_eid_map': # not added yet
                #out[key] = eid_to_eid_map
            else:
                self.log.error('missing map %r' % key)
        return out

    def get_node_ids_with_elements(self, eids: List[int], msg: str='') -> Set[int]:
        """
        Get the node IDs associated with a list of element IDs

        Parameters
        ----------
        eids : List[int]
            list of element ID
        msg : str
            An additional message to print out if an element is not found

        Returns
        -------
        node_ids : Set[int]
            set of node IDs

        For example::

          eids = [1, 2, 3]  # list of elements with pid=1
          msg = ' which are required for pid=1'
          node_ids = bdf.get_node_ids_with_elements(eids, msg=msg)

        """
        if isinstance(eids, integer_types):
            eids = [eids]

        nids2 = set()
        for eid in eids:
            element = self.Element(eid, msg=msg)
            self.log.debug("element.pid = %s" % (element.pid))
            nids = set(element.node_ids)
            nids2.update(nids)
        return nids2

    def get_elements_nodes_by_property_type(
            self,
            dtype: str='int32',
            save_element_types: bool=False) -> Tuple[Dict[Tuple[str, int], Tuple[List[int], List[int]]]]:
        """
        Gets a dictionary of (etype, pid) to [eids, node_ids]

        Parameters
        ----------
        dtype : str; default='int32'
            the type of the integers
        save_element_types : bool; default=False
            adds the etype_to_eids_pids_nids output

        Returns
        -------
        etype_pid_to_eids_nids : dict[(etype, pid)] : [eids, nids]
            etype : str
                the element type
            pid : int
                the property id
                CONRODS have a pid of 0
            eids : (neids, ) int ndarray
                the elements with the property id of pid
            nids : (neids, nnodes/element) int ndarray
                the nodes corresponding to the element
        etype_to_eids_pids_nids : dict[etype] : [eids, pids, nids]
            Enabled by save_element_types; default=None
            etype : str
                the element type
            eids : (neids, ) int ndarray
                the elements with the property id of pid
            pids : (neids, ) int ndarray
                the property ids
                CONRODS have a pid of 0
            nids : (neids, nnodes/element) int ndarray
                the nodes corresponding to the element

        """
        etype_to_eids_pids_nids = self.get_elements_properties_nodes_by_element_type(dtype=dtype)
        output = {}
        for etype, (eids, pids, nids) in etype_to_eids_pids_nids.items():
            upids = np.unique(pids)
            for upid in upids:
                ipid = np.where(pids == upid)[0]
                output[(etype, upid)] = [eids[ipid], nids[ipid, :]]

        if save_element_types:
            return output, etype_to_eids_pids_nids
        return output, None

    def get_element_nodes_by_element_type(self, dtype='int32', solids=None):
        # type: (str, bool) -> Any
        """see ``get_elements_properties_nodes_by_element_type``"""
        self.deprecated('get_element_nodes_by_element_type',
                        'get_elements_properties_nodes_by_element_type', '1.2')
        return self.get_elements_properties_nodes_by_element_type(
            dtype=dtype, solids=solids)

    def _upcast_int_dtype(self, dtype: str) -> str:
        """helper for 64-bit integers"""
        if dtype == 'int32' and len(self.nodes) and max(self.nodes) > 2147483647:
            # or max(self.elements) > 2147483647):
            dtype = 'int64'
        return dtype

    def get_elements_properties_nodes_by_element_type(self,
                                                      dtype: str='int32',
                                                      solids: Optional[Dict[str, Any]]=None,
                                                      stop_if_no_eids: bool=True) -> np.array:
        """
        Gets a dictionary of element type to [eids, pids, node_ids]

        Parameters
        ----------
        dtype : str; default='int32'
            the type of the integers
        solids : dict[etype] : value
            etype : str
                the element type
                should only be CTETRA, CHEXA, CPENTA, CPYRAM
            value : varies
                (nnodes_min, nnodes_max) : Tuple(int, int)
                    the min/max number of nodes for the element
                (nnodes, ) : Tuple(int, )
                    the number of nodes
                    useful if you only have CTETRA4s or only want CTETRA10s
                    fails if you're wrong (and too low)

        Returns
        -------
        etype_to_eids_pids_nids : dict[etype] : [eids, pids, nids]
            etype : str
                the element type
            eids : (neids, ) int ndarray
                the elements with the property id of pid
            pids : (neids, ) int ndarray
                the property ids
                CONRODS have a pid of 0
            nids : (neids, nnodes/element) int ndarray
                the nodes corresponding to the element

        """
        dtype = self._upcast_int_dtype(dtype)

        etypes_no_pids = {
            'CELAS4', 'CDAMP4', 'CHBDYG', 'GENEL',
        }

        etypes = {
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CROD', 'CONROD', 'CTUBE',
            'CBAR', 'CBEAM', 'CBEND', 'CBEAM3',
            'CSHEAR', 'CVISC',
            'CTRIA3', 'CTRIA6', 'CTRIAR',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
            'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
            #'CPLSTS3', 'CPLSTS6', 'CPLSTS4', 'CPLSTS8',
            'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6',
            'CQUADX', 'CQUADX4', 'CQUADX8',
            'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
            'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CFAST', 'CGAP',

            # not supported
            'GENEL', 'CHBDYG',
        }
        output = {}

        if solids is None:
            solids = {
                'CTETRA' : (4, 10),
                #'CTETRA' : (10, ),
                'CHEXA' : (8, 20),
                'CPENTA' : (6, 15),
                'CPYRAM' : (5, 13),
            }

        etypes_found = []
        for etype in etypes:
            if etype not in self._type_to_id_map:
                continue
            eids_list = self._type_to_id_map[etype]
            if not eids_list:
                continue
            etypes_found.append(etype)
            eids = np.array(eids_list, dtype=dtype)
            neids = len(eids)
            eid0 = eids[0]

            elem0 = self.elements[eid0]
            nnodes = len(elem0.nodes)

            if etype not in solids or len(solids[etype]) == 1:
                pids = np.zeros(neids, dtype=dtype)
                nids = np.zeros((neids, nnodes), dtype=dtype)
                for i, eid in enumerate(eids):
                    elem = self.elements[eid]
                    if elem.type in etypes_no_pids:
                        pid = 0
                    else:
                        pid = elem.Pid()
                    assert pid is not None, elem
                    pids[i] = pid
                    nidsi = elem.node_ids
                    try:
                        nids[i, :] = nidsi
                    except TypeError:
                        #print(elem)
                        #print('nidsi =', nidsi)
                        nidsi2 = [nid if nid is not None else 0
                                  for nid in nidsi]
                        try:
                            nids[i, :] = nidsi2
                        except:
                            print(elem)
                            print(nidsi)
                            print(nidsi2)
                            raise
                output[etype] = [eids, pids, nids]
            else:
                # SOLID elements can be variable length
                nnodes_min = min(solids[etype])
                nnodes_max = max(solids[etype])
                pids = np.zeros(neids, dtype='int32')
                nids = np.zeros((neids, nnodes_max), dtype=dtype)
                ieids_max = []
                ieids_min = []
                for i, eid in enumerate(eids):
                    elem = self.elements[eid]
                    pid = elem.Pid()
                    assert pid is not None, elem
                    pids[i] = pid
                    nidsi = elem.node_ids
                    nnodesi = len(nidsi)
                    if nnodesi == nnodes_max:
                        ieids_max.append(i)
                    else:
                        ieids_min.append(i)
                    #self.log.info(str(elem))
                    try:
                        nids[i, :nnodesi] = nidsi
                    except TypeError:
                        #print(elem)
                        #print('nidsi =', nidsi)
                        nidsi2 = [nid  if nid is not None else 0
                                  for nid in nidsi]
                        try:
                            nids[i, :] = nidsi2
                        except:
                            raise
                if len(ieids_max):
                    etype_max = elem.type + str(nnodes_max)
                    ieids_max = np.array(ieids_max, dtype=dtype)
                    output[etype_max] = [eids[ieids_max], pids[ieids_max], nids[ieids_max, :]]
                if len(ieids_min):
                    etype_min = elem.type + str(nnodes_min)
                    ieids_min = np.array(ieids_min, dtype=dtype)
                    output[etype_min] = [eids[ieids_min], pids[ieids_min], nids[ieids_min, :nnodes_min]]
        if stop_if_no_eids:
            msg = (
                'get_elements_properties_nodes_by_element_type output is empty; '
                'nelements=%s; etypes_found=%s' % (
                    len(self.elements), etypes_found)) # etypes_found
            self.log.warning(msg)
        else:
            assert len(output), 'get_elements_properties_nodes_by_element_type output is empty...'
        return output

    #--------------------
    # ELEMENT CARDS

    def get_element_ids_list_with_pids(self, pids: Optional[List[int]]=None) -> List[int]:
        """
        Gets all the element IDs with a specific property ID.

        Parameters
        ----------
        pids : List[int]; default=None -> all
            list of property ID

        Returns
        -------
        element_ids : List[int]
            the element ids

        For example, we want to get all the element ids with ``pids=[1, 2, 3]``

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename)
           pids = [1, 2, 3]
           eids_list = model.get_element_ids_list_with_pids(pids)
           >>> eids_list
           [10, 11, 20, 21, 30, 31]

        """
        etypes_no_pids = [
            'CELAS4', 'CDAMP4', 'CHBDYG', 'GENEL',
        ]
        if pids is None:
            pids = set(list(self.properties.keys()))
        elif isinstance(pids, integer_types):
            pids = {pids}
        else:
            assert isinstance(pids, (list, tuple)), 'pids=%s type=%s' % (pids, type(pids))

        eids2 = []
        for eid, element in sorted(self.elements.items()):
            if element.type in etypes_no_pids:
                pid = 0
            else:
                pid = element.Pid()
            if pid in pids:
                eids2.append(eid)
        return eids2

    def get_pid_to_node_ids_and_elements_array(self,
                                               pids: Union[List[int], int, None]=None,
                                               etypes: Optional[List[str]]=None,
                                               idtype: str='int32',
                                               msg: str='') -> Tuple[Dict[int, str], np.ndarray]:
        """
        a work in progress

        Parameters
        ----------
        pids : List[int]
            list of property ID
        etypes : List[str]
            element types to consider

        Returns
        -------
        pid_to_eids_ieids_map : dict[(pid, etype)] = eid_ieid
            eid_ieid : (Nelements, 2) int ndarray
                eid is the element id
                ieid is the index in the node_ids array
        node_ids : (nelements, nnodes) int ndarray
            nelements : int
                the number of elements in the property type
            nnodes : int
                varies based on the element type

        """
        if pids is None:
            pids = list(self.properties)
        elif isinstance(pids, integer_types):
            pids = [int]

        assert isinstance(pids, (list, tuple, np.ndarray)), 'pids=%s type=%s' % (pids, type(pids))
        pid_to_eids_map = {}
        for pid in pids:
            pid_to_eids_map[pid] = []

        skip_elements = ['CONROD']
        etypes_ = self._slot_to_type_map['elements']
        etype_to_nids_map = {}
        pid_to_eids_ieids_map = defaultdict(list)

        etypes_no_pids = [
            'CELAS4', 'CDAMP4', 'CHBDYG', 'GENEL',
        ]
        etypes_none_nodes = [
            'CELAS1', 'CELAS2', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP4', 'CDAMP5',
            'CBUSH', 'CBUSH1D', 'CFAST',
            'CTRIAX', 'CQUADX', 'CTRIAX6',
            'CTRIA6', 'CQUAD8', 'CQUAD',
            'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
            'CRAC2D', 'CRAC3D', 'CHBDYP', #'CHBDYG',
            'CHACAB', 'CAABSF',
        ]

        if etypes is None:
            etypes = etypes_

        try:
            for etype in etypes_:
                if etype not in etypes_:
                    continue
                eids = self._type_to_id_map[etype]
                if len(eids) == 0:
                    continue
                if etype in skip_elements:
                    self.log.warning('skipping etype=%s because there are no properties%s' % (
                        etype, msg))
                    continue

                # get the number of nodes of the first element
                eid = eids[0]
                element0 = self.elements[eid]
                nnodes = len(element0.node_ids)

                neids = len(eids)
                node_ids = np.zeros((neids, nnodes), dtype=idtype)
                if etype in etypes_none_nodes:
                    for ieid, eid in enumerate(eids):
                        element = self.elements[eid]
                        try:
                            node_ids[ieid, :] = [nid  if nid is not None else 0
                                                 for nid in element.node_ids]
                        except:
                            self.log.error('This error can occur when you have '
                                           'linear and quadratic solid elements '
                                           'within the same model\n%s' % element)
                            raise
                        if etype in etypes_no_pids:
                            pid = 0
                        else:
                            pid = element.Pid()
                        #nids_to_pids_map[]
                        pid_to_eids_ieids_map[(pid, etype)].append((eid, ieid))
                else:
                    try:
                        for ieid, eid in enumerate(eids):
                            element = self.elements[eid]
                            try:
                                node_ids[ieid, :] = element.node_ids
                            except TypeError:
                                print(element)
                                raise
                            if etype in etypes_no_pids:
                                pid = 0
                            else:
                                pid = element.Pid()
                            #nids_to_pids_map[]
                            pid_to_eids_ieids_map[(pid, etype)].append((eid, ieid))
                    except TypeError:
                        print(etype)
                        print(element)
                        raise
                etype_to_nids_map[etype] = node_ids
            for key, value in pid_to_eids_ieids_map.items():
                pid_to_eids_ieids_map[key] = np.array(value, dtype=idtype)
        except OverflowError:
            assert idtype == 'int32', 'idtype=%r while overflowing...' % idtype
            pid_to_eids_ieids_map = self.get_pid_to_node_ids_and_elements_array(
                pids=pids, etypes=etypes, idtype='int64')
        return pid_to_eids_ieids_map

    def get_element_ids_dict_with_pids(self,
                                       pids: Union[List[int], int, None]=None,
                                       stop_if_no_eids: bool=True, msg: str='') -> Dict[int, List[int]]:
        """
        Gets all the element IDs with a specific property ID.

        Parameters
        ----------
        pids : List[int] / int
            list of property ID
        stop_if_no_eids : bool; default=True
            prevents crashing if there are no elements
            setting this to False really doesn't make sense for non-DMIG models

        Returns
        -------
        element_ids : dict[pid] = List[eid]
            dictionary of lists by property
            pid : int
                property id
            eid : int
                element id

        For example, we want all the elements with ``pids=[4, 5, 6]``,
        but we want them in separate groups

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename)
           pids = [4, 5, 6]
           eids_dict = model.get_element_ids_dict_with_pids(pids)
           >>> eids_dict
           {
               4 : [40, 41],
               5 : [50, 51],
               6 : [60, 61],
           }

           # consider all properties
           eids_dict = model.get_element_ids_dict_with_pids()
           {
               1 : [1, 2, 3],
               4 : [40, 41],
               5 : [50, 51],
               6 : [60, 61],
           }

        Notes
        -----
        What happens with CONRODs?

        """
        if pids is None:
            pids = list(self.properties)
        elif isinstance(pids, integer_types):
            pids = [pids]

        assert isinstance(pids, (list, tuple, np.ndarray)), 'pids=%s type=%s' % (pids, type(pids))
        pid_to_eids_map = {}
        for pid in pids:
            pid_to_eids_map[pid] = []

        elem_count = 0
        #element_type_to_dmap_id = {
            #'CONROD' : -10,
            #'CELAS2' : -12,
            #'CELAS4' : -14,
            #'CDAMP2' : -21,
            #'CDAMP4' : -23,
            #'CHBDYG' : -108,
        #}
        elements_without_properties = {
            'CONROD', 'CELAS2', 'CELAS4', 'CDAMP2', 'CDAMP4',
            'CHBDYG', 'GENEL'}
        for eid, element in self.elements.items():
            try:
                pid = element.Pid()
            except AttributeError:
                if element.type in elements_without_properties:
                    continue
                print(element)
                raise
            if pid in pids:
                pid_to_eids_map[pid].append(eid)
            elem_count += 1

        if elem_count == 0 and stop_if_no_eids:
            raise RuntimeError('no elements with properties found%s\ncard_count=%s' % (
                msg, str(self.card_count)))
        elif elem_count == 0:
            self.log.warning('no elements with properties found%s' % msg)
        return pid_to_eids_map

    def get_node_id_to_element_ids_map(self) -> Dict[int, List[int]]:
        """
        Returns a dictionary that maps node IDs to a list of elemnent IDs

        .. todo:: support 0d or 1d elements
        .. todo:: support elements with missing nodes
                  (e.g. CQUAD8 with missing nodes)

        """
        nid_to_eids_map = {}
        for nid in self.nodes:  # initalize the mapper
            nid_to_eids_map[nid] = []

        if self.spoints:  # SPOINTs
            for nid in sorted(self.spoints):  # SPOINTs
                nid_to_eids_map[nid] = []

        skip_cards = {
            'CCONEAX',
        }
        for (eid, element) in self.elements.items():  # load the mapper
            if element.type in skip_cards:
                continue
            try:
                # not supported for 0-D and 1-D elements
                nids = element.node_ids
            except AttributeError:
                print(element.type)
            else:
                for nid in nids:  # (e.g. CQUAD8 with missing node)
                    if nid:
                        nid_to_eids_map[nid].append(eid)

        return nid_to_eids_map

    def get_node_id_to_elements_map(self) -> Dict[int, List[int]]:
        """
        Returns a dictionary that maps node IDs to a list of elements.

        Returns
        -------
        nid_to_elements_map : Dict[nid]=List[eid]
            node id to a list of elements

        .. todo:: support 0d or 1d elements
        .. todo:: support elements with missing nodes
                  (e.g. CQUAD8 with missing nodes)

        """
        nid_to_elements_map = {}
        for nid in self.nodes:  # initalize the mapper
            nid_to_elements_map[nid] = []

        for nid in self.spoints:
            nid_to_elements_map[nid] = []
        for nid in self.epoints:
            nid_to_elements_map[nid] = []

        skip_cards = {
            'CCONEAX',
        }
        for element in self.elements.values():  # load the mapper
            if element.type in skip_cards:
                continue

            try:
                # not supported for 0-D and 1-D elements
                nids = element.node_ids
            except AttributeError:
                print(element.type)
            else:
                for nid in nids:  # (e.g. CQUAD8 with missing node)
                    if nid:
                        nid_to_elements_map[nid].append(element)

        return nid_to_elements_map

    def get_property_id_to_element_ids_map(self, msg: str='') -> Dict[int, List[int]]:
        """
        Returns a dictionary that maps a property ID to a list of elements.

        Returns
        -------
        pid_to_eids_map : Dict[pid]=List[eid]
            property id to a list of elements
        msg : str; default=''
            a message added to the error message

        """
        pid_to_eids_map = {}
        pids = self.property_ids
        for pid in pids:
            pid_to_eids_map[pid] = []
        #for pid in self.phbdys.keys():
            #assert pid not in pid_to_eids_map, 'pid=%s is already used and must be used by PHBDY' % pid
            #pid_to_eids_map[pid] = []

        elements_without_properties = {
            'CONROD', 'CONM2', 'CELAS2', 'CELAS4', 'CDAMP2', 'CDAMP4',
            'GENEL', 'CHACAB', 'CAABSF', }
        thermal_elements = {'CHBDYP'}
        elements_without_properties.update(thermal_elements)
        skip_elements = elements_without_properties

        for eid in self.element_ids:
            element = self.Element(eid)
            element_type = element.type
            if element_type in skip_elements:
                continue
            if hasattr(element, 'pid'):
                pid = element.Pid()
                if pid < 0: # CTRIAX6
                    continue
                try:
                    pid_to_eids_map[pid].append(eid)
                except KeyError:
                    print(element)
                    raise KeyError('pid=%s is invalid for card%s=\n%s' % (pid, msg, str(element)))
        return pid_to_eids_map

    def get_material_id_to_property_ids_map(self, msg: str='') -> Dict[int, List[int]]:
        """
        Returns a dictionary that maps a material ID to a list of properties

        Returns
        -------
        mid_to_pids_map : dict[int] = int
            the mapping
        msg : str; default=''
            a message added to the error message

        .. code-block:: python

           >>> mid_to_pid_map = get_material_id_to_property_ids_map()
           >>> mid = 1
           >>> pids = get_material_id_to_property_ids_map[mid]
           >>> pids
           [1, 2, 3]

        .. note:: all properties require an mid to be counted (except for
                PCOMP, which has multiple mids)

        """
        mid_to_pids_map = {}
        mids = self.get_material_ids()
        for mid in mids:
            mid_to_pids_map[mid] = []

        properties_without_materials = {
            'PGAP', 'PELAS', 'PVISC', 'PBUSH', 'PDAMP', 'PFAST', 'PBUSH1D',
            'PACABS', 'PAABSF', 'PACBAR',
        }

        for pid in self.property_ids:
            prop = self.Property(pid)
            prop_type = prop.type
            if prop_type in ['PCOMP', 'PCOMPG']:
                mids = prop.Mids()

                for mid in mids:
                    if pid not in mid_to_pids_map[mid]:
                        mid_to_pids_map[mid].append(pid)
                    else:  # PCOMP
                        if hasattr(prop, 'mid') and prop.Mid() in mids:
                            if pid not in mid_to_pids_map[mid]:
                                mid_to_pids_map[mid].append(pid)
            elif prop_type in properties_without_materials:
                pass
            elif prop_type in ['PSHELL']:
                mids = prop.material_ids
                for i, mid in enumerate(mids):
                    if mid is None or mid == 0:
                        continue
                    try:
                        mid_to_pids_map[mid].append(pid)
                    except KeyError:
                        print(prop)
                        raise KeyError('i=%s mid=%s is invalid for card%s=\n%s' % (
                            i, mid, msg, str(prop)))
            else:
                mid = prop.Mid()
                try:
                    mid_to_pids_map[mid].append(pid)
                except KeyError:
                    print(prop)
                    raise KeyError('mid=%s is invalid for card %s=\n%s' % (mid, msg, str(prop)))
        return mid_to_pids_map

    def get_reduced_nsms(self, nsm_id: int,
                         consider_nsmadd: bool=True,
                         stop_on_failure: bool=True) -> List[Any]:
        """
        Get all traced NSMs that are part of a set

        Parameters
        ----------
        nsm_id : int
            the NSM id
        consider_nsmadd : bool
            NSMADDs should not be considered when referenced from an NSMADD
            from a case control, True should be used.
        stop_on_failure : bool; default=True
            errors if parsing something new

        Returns
        -------
        mpcs : List[NSM]
            the various NSMs

        """
        if not isinstance(nsm_id, integer_types):
            msg = 'nsm_id must be an integer; type=%s, nsm_id=\n%r' % (type(nsm_id), nsm_id)
            raise TypeError(msg)

        try:
            nsms = self.NSM(nsm_id, consider_nsmadd=consider_nsmadd)
        except KeyError:
            if stop_on_failure:
                raise
            self.log.error("could not find expected NSM id=%s" % nsm_id)
            return []

        nsms2 = []
        for nsm in nsms:
            if nsm.type == 'NSMADD':
                for nsmi in nsm.nsm_ids:
                    if isinstance(nsmi, list):
                        for nsmii in nsmi:
                            if isinstance(nsmii, integer_types):
                                nsmiii = nsmii
                            else:
                                nsmiii = nsmii.conid
                            nsms2i = self.get_reduced_nsms(
                                nsmiii, consider_nsmadd=False,
                                stop_on_failure=stop_on_failure)
                            nsms2 += nsms2i
                    else:
                        assert isinstance(nsmi, integer_types), nsmi
                        nsms2i = self.get_reduced_nsms(
                            nsmi, consider_nsmadd=False,
                            stop_on_failure=stop_on_failure)
                        nsms2 += nsms2i
            else:
                nsms2.append(nsm)
        return nsms2

    def get_reduced_mpcs(self, mpc_id: int,
                         consider_mpcadd: bool=False,
                         stop_on_failure: bool=True) -> List[Any]:
        """
        Get all traced MPCs that are part of a set

        Parameters
        ----------
        mpc_id : int
            the MPC id
        consider_mpcadd : bool
            MPCADDs should not be considered when referenced from an MPCADD
            from a case control, True should be used.
        stop_on_failure : bool; default=True
            errors if parsing something new

        Returns
        -------
        mpcs : List[MPC]
            the various MPCs

        """
        if not isinstance(mpc_id, integer_types):
            msg = 'mpc_id must be an integer; type=%s, mpc_id=\n%r' % (type(mpc_id), mpc_id)
            raise TypeError(msg)

        try:
            mpcs = self.MPC(mpc_id, consider_mpcadd=consider_mpcadd)
        except KeyError:
            if stop_on_failure:
                raise
            self.log.error("could not find expected MPC id=%s" % mpc_id)
            return []

        mpcs2 = []
        for mpc in mpcs:
            if mpc.type == 'MPCADD':
                for mpci in mpc.mpc_ids:
                    if isinstance(mpci, list):
                        for mpcii in mpci:
                            if isinstance(mpcii, integer_types):
                                mpciii = mpcii
                            else:
                                mpciii = mpcii.conid
                            mpcs2i = self.get_reduced_mpcs(
                                mpciii, consider_mpcadd=False,
                                stop_on_failure=stop_on_failure)
                            mpcs2 += mpcs2i
                    else:
                        assert isinstance(mpci, integer_types), mpci
                        mpcs2i = self.get_reduced_mpcs(
                            mpci, consider_mpcadd=False,
                            stop_on_failure=stop_on_failure)
                        mpcs2 += mpcs2i
            else:
                mpcs2.append(mpc)
        return mpcs2

    def get_reduced_spcs(self, spc_id: int,
                         consider_spcadd: bool=True,
                         stop_on_failure: bool=True) -> List[Any]:
        """
        Get all traced SPCs that are part of a set

        Parameters
        ----------
        spc_id : int
            the SPC id
        consider_spcadd : bool
            SPCADDs should not be considered when referenced from an SPCADD
            from a case control, True should be used.
        stop_on_failure : bool; default=True
            errors if parsing something new

        Returns
        -------
        spcs : List[SPC]
            the various SPCs

        """
        if not isinstance(spc_id, integer_types):
            msg = 'spc_id must be an integer; type=%s, spc_id=\n%r' % (type(spc_id), spc_id)
            raise TypeError(msg)

        try:
            spcs = self.SPC(spc_id, consider_spcadd=consider_spcadd)
        except KeyError:
            if stop_on_failure:
                raise
            self.log.error("could not find expected SPC id=%s" % spc_id)
            return []

        spcs2 = []
        for spc in spcs:
            if spc.type == 'SPCADD':
                for spci in spc.sets:
                    if isinstance(spci, list):
                        for spcii in spci:
                            if isinstance(spcii, integer_types):
                                spciii = spcii
                            else:
                                spciii = spcii.conid
                            spcs2i = self.get_reduced_spcs(
                                spciii,
                                consider_spcadd=False,
                                stop_on_failure=stop_on_failure)
                            spcs2 += spcs2i
                    else:
                        assert isinstance(spci, integer_types), spci
                        spcs2i = self.get_reduced_spcs(
                            spci,
                            consider_spcadd=False,
                            stop_on_failure=stop_on_failure)
                        spcs2 += spcs2i
            else:
                spcs2.append(spc)
        return spcs2

    def get_spcs(self, spc_id: int,
                 consider_nodes: bool=False,
                 stop_on_failure: bool=True) -> Tuple[List[int], List[str]]:
        """
        Gets the SPCs in a semi-usable form.

        Parameters
        ----------
        spc_id : int
            the desired SPC ID
        consider_nodes : bool; default=False
            True : consider the GRID card PS field
            False: consider the GRID card PS field

        Returns
        -------
        nids : List[int]
            the constrained nodes
        comps : List[str]
            the components that are constrained on each node

        Considers:
          - SPC
          - SPC1
          - SPCADD
          - GRID

        Doesn't consider:
          - non-zero enforced value on SPC
          - GMSPC

        """
        warnings = ''
        spcs = self.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=stop_on_failure)
        nids = []
        comps = []
        for spc in spcs:
            if spc.type == 'SPC1':
                nodes = spc.nodes
                nnodes = len(nodes)
                nids += nodes
                comps += [str(spc.components)] * nnodes
            elif spc.type == 'SPC':
                for nid, comp, unused_enforced in zip(spc.nodes, spc.components, spc.enforced):
                    nids.append(nid)
                    comps.append(comp)
            else:
                warnings += str(spc)
                #raise NotImplementedError(spc.type)
        if warnings:
            self.log.warning("get_spcs doesn't consider:\n%s" % warnings.rstrip('\n'))

        if consider_nodes:
            for nid, node in self.nodes.items():
                if node.ps:
                    nids.append(nid)
                    comps.append(node.ps)
        return nids, comps

    def get_mklist(self) -> np.ndarray:
        """gets the MKLIST vector from MKAERO1/MKAERO2"""
        mklist = []
        mkarray = np.array([])
        for mkaero in self.mkaeros:
            mklist += mkaero.mklist()
        if mklist:
            mkarray = np.hstack([mklist])
            #new_array = [tuple(row) for row in mkarray]
            #unique_pairs = np.lib.arraysetops.unique(new_array, axis=0).tolist()
        return mkarray
