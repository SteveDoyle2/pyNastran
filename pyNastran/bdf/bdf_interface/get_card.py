# pylint: disable=E1101,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from collections import defaultdict
from six import string_types, iteritems, iterkeys, itervalues

import numpy as np

from pyNastran.bdf.bdf_interface.get_methods import GetMethods
#from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.utils import integer_types


class GetCard(GetMethods):
    def __init__(self):
        self._type_to_slot_map = {}
        GetMethods.__init__(self)

    def get_card_ids_by_card_types(self, card_types=None, reset_type_to_slot_map=False,
                                   stop_on_missing_card=False, combine=False):
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

        Example 1
        ---------
        out_dict = model.get_card_ids_by_card_types(
            card_types=['GRID', 'CTRIA3', 'CQUAD4'], combine=False)
        out_dict = {
            'GRID' : [1, 2, 10, 42, 1000],
            'CTRIA3' : [1, 2, 3, 5],
            'CQUAD4' : [4],
        }

        Example 2 - Shell Elements
        --------------------------
        out_dict = model.get_card_ids_by_card_types(
            card_types=['CTRIA3', 'CQUAD4'], combine=True)
        out_dict = {
            [1, 2, 3, 4, 5],
        }
        """
        if card_types is None:
            card_types = list(self.cards_to_read)

        if isinstance(card_types, string_types):
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
            for key, value in sorted(iteritems(out_dict)):
                out_list += value
            return out_list
        return out_dict

    def _reset_type_to_slot_map(self):
        rslot_map = defaultdict(list)
        for dict_name, card_names in iteritems(self._slot_to_type_map):
            print('card_names=%s dict_name=%s' % (card_names, dict_name))
            card_name0 = card_names[0]
            if card_name0 in ['DTABLE', 'GRDSET', 'SESUP', 'DOPTPRM', 'MONPNT1', 'SUPORT',
                              'MKAERO1', 'MATHP']:
                pass
            else:
                adict = getattr(self, dict_name)
                if isinstance(adict, dict):
                    for key, card in iteritems(adict):
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
                        else:
                            if value.type in ['CSET1', 'CSET']:
                                pass
                                #rslot_map[value.type] = value.
                            else:
                                raise NotImplementedError('list; names=%s' % card_names)
                else:
                    raise NotImplementedError('%s; names=%s' % (type(adict), card_names))
        return rslot_map

    def get_rslot_map(self, reset_type_to_slot_map=False):
        if (reset_type_to_slot_map or self._type_to_slot_map is None or
                len(self._type_to_slot_map) == 0):
            rslot_map = {}
            for key, values in iteritems(self._slot_to_type_map):
                for value in values:
                    rslot_map[value] = key
            self._type_to_slot_map = rslot_map
        else:
            rslot_map = self._type_to_slot_map
        assert 'GRID' in rslot_map
        return rslot_map

    def get_cards_by_card_types(self, card_types, reset_type_to_slot_map=False,
                                stop_on_missing_card=False):
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
            raise TypeError('card_types must be a list/tuple; type=%s' % type(card_types))

        #self._type_to_id_map = {
        #    'CQUAD4' : [1, 2, 3]
        #}
        #self._slot_to_type_map = {'elements' : [CQUAD4, CTRIA3]}
        rslot_map = self.get_rslot_map(reset_type_to_slot_map=False)

        out = {}
        for card_type in card_types:
            if card_type not in self.card_count:
                if stop_on_missing_card:
                    raise RuntimeError('%r is not in the card_count; keys=%s'
                                       % str(sorted(self.card_count.keys())))
                out[card_type] = []
                continue

            #print('card_type=%r' % card_type)
            try:
                key = rslot_map[card_type]  # update attributes.py ~line 520
            except:
                self.log.error("card_type=%r' hasn't been added to "
                               "self._slot_to_type_map...check for typos")
                raise
            slot = getattr(self, key)
            ids = self._type_to_id_map[card_type]
            cards = []
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

    def get_SPCx_node_ids(self, spc_id, exclude_spcadd=False, stop_on_failure=True):
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        ----------
        spc_id : int
            the SPC id
        exclude_spcadd : bool
            you can exclude SPCADD if you just want a list of all the
            SPCs in the model.  For example, apply all the SPCs when
            there is no SPC=N in the case control deck, but you don't
            need to apply SPCADD=N twice.
        stop_on_failure : bool; default=True
            errors if parsing something new
        """
        try:
            spcs = self.spcs[spc_id]
        except KeyError:
            self.log.warning('spc_id=%s not found' % spc_id)
            return []
        except TypeError:
            print('spc_id=%r' % spc_id)
            raise


        node_ids = []
        for card in spcs:
            if card.type == 'SPC':
                nids = card.node_ids
            elif card.type == 'SPC1':
                nids = card.node_ids
            elif card.type == 'SPCADD':
                nids = []
                for new_spc_id in card.spc_ids:
                    nidsi = self.get_SPCx_node_ids(new_spc_id, exclude_spcadd=False,
                                                   stop_on_failure=stop_on_failure)
                    nids += nidsi
            else:
                self.log.warning('get_SPCx_node_ids doesnt supprt %r' % card.type)
                continue
            node_ids += nids
        return node_ids

    def get_SPCx_node_ids_c1(self, spc_id, exclude_spcadd=False, stop_on_failure=True):
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        ----------
        spc_id : int
            the SPC id
        exclude_spcadd : bool
            you can exclude SPCADD if you just want a list of all the
            SPCs in the model.  For example, apply all the SPCs when
            there is no SPC=N in the case control deck, but you don't
            need to apply SPCADD=N twice.
        stop_on_failure : bool; default=True
            errors if parsing something new
        """
        try:
            spcs = self.spcs[spc_id]
        except KeyError:
            self.log.warning('spc_id=%s not found' % spc_id)
            return {}

        node_ids_c1 = defaultdict(str)
        #print('spcs = ', spcs)
        for card in spcs:  # used to be sorted(spcs)
            if card.type == 'SPC':
                for nid, c1 in zip(card.node_ids, card.components):
                    assert nid is not None, card.node_ids
                    node_ids_c1[nid] += c1
            elif card.type == 'SPC1':
                nids = card.node_ids
                c1 = card.components
                for nid in nids:
                    node_ids_c1[nid] += c1
            elif card.type == 'SPCADD':
                nids = []
                for new_spc_id in card.spc_ids:
                    nids_c1i = self.get_SPCx_node_ids_c1(new_spc_id, exclude_spcadd=False,
                                                         stop_on_failure=stop_on_failure)
                    for nid, c1 in iteritems(nids_c1i):
                        node_ids_c1[nid] += c1
            elif card.type in ['GMSPC', 'SPCAX']:
                self.log.warning('get_SPCx_node_ids doesnt supprt %r' % card.type)
            else:
                msg = 'get_SPCx_node_ids_c1 doesnt supprt %r' % card.type
                if stop_on_failure:
                    raise RuntimeError(msg)
                else:
                    self.log.warning(msg)
        return node_ids_c1

    def get_MPCx_node_ids_c1(self, mpc_id, exclude_mpcadd=False, stop_on_failure=True):
        r"""
        Get the MPC/MPCADD IDs.

        Parameters
        ----------
        mpc_id : int
            the MPC id
        exclude_mpcadd : bool
            you can exclude MPCADD if you just want a list of all the
            MPCs in the model.  For example, apply all the MPCs when
            there is no MPC=N in the case control deck, but you don't
            need to apply MPCADD=N twice.
            TODO: not used
        stop_on_failure : bool; default=True
            errors if parsing something new

        I      I
          \   /
        I---D---I
        """
        lines = []
        if not isinstance(mpc_id, integer_types):
            msg = 'mpc_id must be an integer; type=%s, mpc_id=\n%r' % (type(mpc_id), mpc_id)
            raise TypeError(msg)

        try:
            mpcs = self.mpcs[mpc_id]
        except KeyError:
            self.log.warning('mpc_id=%s not found' % mpc_id)
            return []

        # dependent, independent
        for card in mpcs:
            if card.type == 'MPC':
                nids = card.node_ids
                nid0 = nids[0]
                #constraint0 = card.constraints[0]
                #enforced0 = card.enforced[0]
                #card.constraints[1:]
                for nid, enforced in zip(nids[1:], card.enforced[1:]):
                    if enforced != 0.0:
                        lines.append([nid0, nid])
            elif card.type == 'MPCADD':
                nids = []
                for new_mpc_id in card.mpc_ids:
                    linesi = self.get_MPCx_node_ids_c1(new_mpc_id, exclude_mpcadd=False)
                    lines += linesi
            else:
                msg = 'get_MPCx_node_ids_c1 doesnt supprt %r' % card.type
                if stop_on_failure:
                    raise RuntimeError(msg)
                else:
                    self.log.warning(msg)
        return lines

    def get_load_arrays(self, subcase_id, nid_map, eid_map, node_ids, normals):
        """
        Gets the following load arrays for the GUI

        Loads include:
         - Temperature
         - Pressure (Centroidal)
         - Forces
         - SPCD

        Parameters
        ----------
        model : BDF()
            the BDF object
        subcase_id : int
            the subcase id

        Returns
        -------
        found_load : bool
            a flag that indicates if load data was found
        found_temperature : bool
            a flag that indicates if temperature data was found
        temperature_data : tuple(temperature_key, temperatures)
            temperature_key : str
                One of the following:
                  TEMPERATURE(MATERIAL)
                  TEMPERATURE(INITIAL)
                  TEMPERATURE(LOAD)
                  TEMPERATURE(BOTH)
            temperatures : (nnodes, 1) float ndarray
                the temperatures
        load_data : tuple(centroidal_pressures, forces, spcd)
            centroidal_pressures : (nelements, 1) float ndarray
                the pressure
            forces : (nnodes, 3) float ndarray
                the pressure
            spcd : (nnodes, 3) float ndarray
                the SPCD load application
        """
        subcase = self.subcases[subcase_id]
        is_loads = False
        is_temperatures = False

        load_keys = (
            'LOAD', 'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
            'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')
        temperature_keys = (
            'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
            'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')

        centroidal_pressures = None
        forces = None
        spcd = None
        temperature_key = None
        temperatures = None
        for key in load_keys:
            try:
                load_case_id = subcase.get_parameter(key)[0]
            except KeyError:
                # print('no %s for isubcase=%s' % (key, subcase_id))
                continue
            try:
                load_case = self.loads[load_case_id]
            except KeyError:
                self.log.warning('LOAD=%s not found' % load_case_id)
                continue

            if key == 'LOAD':
                p0 = np.array([0., 0., 0.], dtype='float32')
                centroidal_pressures, forces, spcd = self._get_forces_moments_array(
                    p0, load_case_id,
                    nid_map=nid_map,
                    eid_map=eid_map,
                    node_ids=node_ids,
                    normals=normals,
                    dependents_nodes=self.node_ids,
                    include_grav=False)
                if centroidal_pressures is not None: # or any of the others
                    is_loads = True
            elif key in temperature_keys:
                is_temperatures, temperatures = self._get_temperatures_array(load_case_id)
                temperature_key = key
            else:
                raise NotImplementedError(key)
        temperature_data = (temperature_key, temperatures)
        load_data = (centroidal_pressures, forces, spcd)
        return is_loads, is_temperatures, temperature_data, load_data

    def _get_dvprel_ndarrays(self, nelements, pids):
        """creates arrays for dvprel results"""
        dvprel_t_init = np.zeros(nelements, dtype='float32')
        dvprel_t_min = np.zeros(nelements, dtype='float32')
        dvprel_t_max = np.zeros(nelements, dtype='float32')
        design_region = np.zeros(nelements, dtype='int32')

        for key, dvprel in iteritems(self.dvprels):
            if dvprel.type == 'DVPREL1':
                prop_type = dvprel.prop_type
                desvars = dvprel.dvids
                coeffs = dvprel.coeffs
                if hasattr(dvprel, 'pid_ref'):
                    pid = dvprel.pid_ref.pid
                else:
                    pid = dvprel.pid
                var_to_change = dvprel.pname_fid
                assert len(desvars) == 1, len(desvars)

                if prop_type == 'PSHELL':
                    i = np.where(pids == pid)
                    design_region[i] = dvprel.oid
                    assert len(i) > 0, i
                    if var_to_change == 'T':
                        #value = 0.
                        lower_bound = 0.
                        upper_bound = 0.
                        for desvar, coeff in zip(desvars, coeffs):
                            if isinstance(desvar, integer_types):
                                desvar_ref = self.desvars[desvar]
                            else:
                                desvar_ref = desvar.desvar_ref
                            xiniti = desvar_ref.xinit
                            if desvar_ref.xlb != -1e20:
                                xiniti = max(xiniti, desvar_ref.xlb)
                                lower_bound = desvar_ref.xlb
                            if desvar_ref.xub != 1e20:
                                xiniti = min(xiniti, desvar_ref.xub)
                                upper_bound = desvar_ref.xub

                            # code validation
                            if desvar_ref.delx is not None and desvar_ref.delx != 1e20:
                                pass

                            # TODO: haven't quite decided what to do
                            if desvar_ref.ddval is not None:
                                msg = 'DESVAR id=%s DDVAL is not None\n%s' % str(desvar_ref)
                            assert desvar_ref.ddval is None, desvar_ref
                            xinit = coeff * xiniti
                        dvprel_t_init[i] = xinit
                        dvprel_t_min[i] = lower_bound
                        dvprel_t_max[i] = upper_bound
                    else:
                        msg = 'var_to_change=%r; dvprel=\n%s' % (var_to_change, str(dvprel))
                        raise NotImplementedError(msg)
                else:
                    msg = 'prop_type=%r; dvprel=\n%s' % (prop_type, str(dvprel))
                    raise NotImplementedError(msg)
            else:
                msg = 'dvprel.type=%r; dvprel=\n%s' % (dvprel.type, str(dvprel))
                raise NotImplementedError(msg)

            # TODO: haven't quite decided what to do
            if dvprel.p_max != 1e20:
                dvprel.p_max

            # TODO: haven't quite decided what to do
            if dvprel.p_min is not None:
                dvprel.p_min
        return dvprel_t_init, dvprel_t_min, dvprel_t_max, design_region

    def _get_forces_moments_array(self, p0, load_case_id,
                                  nid_map, eid_map, node_ids, normals, dependents_nodes,
                                  include_grav=False):
        """
        Gets the forces/moments on the nodes

        Parameters
        ----------
        p0 : (3, ) float ndarray
            the reference location
        load_case_id : int
            the load id
        nidmap : ???
            ???
        eidmap : ???
            ???
        node_ids : ???
            ???
        normals : (nelements, 3) float ndarray
            the normal vectors for the shells
        dependents_nodes : ???
            ???
        include_grav : bool; default=False
            is the mass of the elements considered; unused

        Returns
        -------
        temperature_data : tuple(temperature_key, temperatures)
            temperature_key : str
                One of the following:
                  TEMPERATURE(MATERIAL)
                  TEMPERATURE(INITIAL)
                  TEMPERATURE(LOAD)
                  TEMPERATURE(BOTH)
            temperatures : (nnodes, 1) float ndarray
                the temperatures
        load_data : tuple(centroidal_pressures, forces, spcd)
            centroidal_pressures : (nelements, 1) float ndarray
                the pressure
            forces : (nnodes, 3) float ndarray
                the pressure
            spcd : (nnodes, 3) float ndarray
                the SPCD load application

        Considers
        ---------
        FORCE
        PLOAD2 - CTRIA3, CQUAD4, CSHEAR
        PLOAD4 - CTRIA3, CTRIA6, CTRIAR
                 CQUAD4, CQUAD8, CQUAD, CQUADR, CSHEAR
                 CTETRA, CPENTA, CHEXA
        SPCD
        """
        if not any(['FORCE' in self.card_count, 'PLOAD2' in self.card_count,
                    'PLOAD4' in self.card_count, 'SPCD' in self.card_count]):
            return None, None, None
        nids = sorted(self.nodes.keys())
        nnodes = len(nids)

        load_case = self.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)

        #eids = sorted(self.elements.keys())
        centroidal_pressures = np.zeros(len(self.elements), dtype='float32')
        nodal_pressures = np.zeros(len(self.node_ids), dtype='float32')

        forces = np.zeros((nnodes, 3), dtype='float32')
        spcd = np.zeros((nnodes, 3), dtype='float32')
        # loop thru scaled loads and plot the pressure
        cards_ignored = {}

        assert normals is not None
        fail_nids = set()
        fail_count = 0
        fail_count_max = 3
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'FORCE':
                scale2 = load.mag * scale  # does this need a magnitude?
                nid = load.node
                if nid in dependents_nodes:
                    fail_nids.add(nid)
                    fail_count += 1
                    if fail_count < fail_count_max:
                        print('    nid=%s is a dependent node and has a FORCE applied\n%s' % (
                            nid, str(load)))
                forces[nid_map[nid]] += load.xyz * scale2

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        node_ids = elem.node_ids
                        nnodes = len(node_ids)
                        ie = eid_map[eid]
                        normal = normals[ie, :]

                        area = elem.Area()
                        forcei = pressure * normal * area / nnodes
                        # r = elem.Centroid() - p0
                        # m = cross(r, f)
                        for nid in node_ids:
                            if nid in dependents_nodes:
                                fail_nids.add(nid)
                                fail_count += 1
                                if fail_count < fail_count_max:
                                    print('    nid=%s is a dependent node and has a PLOAD2 applied\n'
                                          '%s' % (nid, str(load)))
                            forces[nid_map[nid]] += forcei
                        forces += forcei
                        # F += f
                        # M += m
                    else:
                        self.log.debug('    case=%s etype=%r loadtype=%r not supported' % (
                            load_case_id, elem.type, load.type))

            elif load.type == 'PLOAD4':
                # single element per PLOAD
                #eid = elem.eid
                #pressures[eids.index(eid)] = p
                #pressure = load.pressures[0] * scale

                # multiple elements
                for elem in load.eids:
                    ie = eid_map[elem.eid]
                    normal = normals[ie, :]
                    # pressures[eids.index(elem.eid)] += p
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR']:
                        area = elem.get_area()
                        elem_node_ids = elem.node_ids
                        nface = len(elem_node_ids)

                        if load.surf_or_line == 'SURF':
                            if np.linalg.norm(load.nvector) != 0.0 or load.Cid() != 0:
                                normal = load.nvector / np.linalg.norm(load.nvector)
                                assert load.Cid() == 0, 'cid=%r on a PLOAD4 is not supported\n%s' % (load.Cid(), str(load))
                        else:
                            msg = 'surf_or_line=%r on PLOAD4 is not supported\n%s' % (
                                load.surf_or_line, str(load))
                            self.log.debug(msg)
                            continue

                        pressures = load.pressures[:nface]
                        if min(pressures) != max(pressures):
                            pressure = np.mean(pressures)
                        else:
                            pressure = pressures[0]

                        forcei = pressure * area * normal / nface
                        for nid in elem_node_ids:
                            if nid in dependents_nodes:
                                fail_nids.add(nid)
                                fail_count += 1
                                if fail_count < fail_count_max:
                                    print('    nid=%s is a dependent node and has a'
                                          ' PLOAD4 applied\n%s' % (nid, str(load)))
                            #forces[nids.index(nid)] += F
                            i = nid_map[nid]
                            try:
                                forces[i, :] += forcei
                            except IndexError:
                                print('i = %s' % i)
                                print('normals.shape = %s' %  str(normals.shape))
                                print('forces.shape = %s' % str(forces.shape))
                                print('normal = ', normal)
                                print('forces[i, :] = ', forces[i, :])
                                raise
                        #nface = 3
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        area = elem.get_area()
                        elem_node_ids = elem.node_ids
                        nface = len(elem_node_ids)

                        if load.surf_or_line == 'SURF':
                            if np.linalg.norm(load.nvector) != 0.0 or load.Cid() != 0:
                                normal = load.nvector / np.linalg.norm(load.nvector)
                                assert load.Cid() == 0, 'cid=%r on a PLOAD4 is not supported\n%s' % (load.Cid(), str(load))
                        else:
                            msg = 'surf_or_line=%r on PLOAD4 is not supported\n%s' % (
                                load.surf_or_line, str(load))
                            self.log.debug(msg)
                            continue

                        pressures = load.pressures[:nface]
                        if min(pressures) != max(pressures):
                            pressure = np.mean(pressures)
                        else:
                            pressure = pressures[0]

                        forcei = pressure * area * normal / nface

                        for nid in elem_node_ids:
                            if nid in dependents_nodes:
                                fail_nids.add(nid)
                                fail_count += 1
                                if fail_count < fail_count_max:
                                    print('    nid=%s is a dependent node and has a'
                                          ' PLOAD4 applied\n%s' % (nid, str(load)))
                            #forces[nids.index(nid)] += F
                            i = nid_map[nid]
                            try:
                                forces[i, :] += forcei
                            except IndexError:
                                print('i = %s' % i)
                                print('normals.shape = %s' %  str(normals.shape))
                                print('forces.shape = %s' % str(forces.shape))
                                print('normal = ', normal)
                                print('forces[i, :] = ', forces[i, :])
                                raise
                            nface = 4
                    else:
                        elem_node_ids = elem.node_ids
                        if elem.type == 'CTETRA':
                            #face1 = elem.get_face(load.g1.nid, load.g34.nid)
                            face, area, centroid, normal = elem.get_face_area_centroid_normal(
                                load.g1.nid, load.g34.nid)
                            #assert face == face1
                            nface = 3
                        elif elem.type == 'CHEXA':
                            #face1 = elem.get_face(load.g34.nid, load.g1.nid)
                            face, area, centroid, normal = elem.get_face_area_centroid_normal(
                                load.g34.nid, load.g1.nid)
                            #assert face == face1
                            nface = 4
                        elif elem.type == 'CPENTA':
                            g1 = load.g1.nid
                            if load.g34 is None:
                                #face1 = elem.get_face(g1)
                                face, area, centroid, normal = elem.get_face_area_centroid_normal(g1)
                                nface = 3
                            else:
                                #face1 = elem.get_face(g1, load.g34.nid)
                                face, area, centroid, normal = elem.get_face_area_centroid_normal(
                                    g1, load.g34.nid)
                                nface = 4
                            #assert face == face1
                        else:
                            msg = ('case=%s eid=%s etype=%r loadtype=%r not supported'
                                   % (load_case_id, eid, elem.type, load.type))
                            self.log.debug(msg)
                            continue

                        pressures = load.pressures[:nface]
                        assert len(pressures) == nface
                        if min(pressures) != max(pressures):
                            pressure = np.mean(pressures)
                            #msg = ('%s%s\npressure.min=%s != pressure.max=%s using average'
                                   #' of %%s; load=%s eid=%%s'  % (
                                       #str(load), str(elem), min(pressures), max(pressures),
                                       #load.sid)
                            #print(msg % (pressure, eid))
                        else:
                            pressure = pressures[0]
                        #centroidal_pressures

                        if  load.surf_or_line == 'SURF':
                            if np.linalg.norm(load.nvector) != 0.0 or load.Cid() != 0:
                                normal = load.nvector / np.linalg.norm(load.nvector)
                                assert load.Cid() == 0, 'cid=%r on a PLOAD4 is not supported\n%s' % (load.Cid(), str(load))
                        else:
                            msg = 'surf_or_line=%r on PLOAD4 is not supported\n%s' % (
                                load.surf_or_line, str(load))
                            self.log.debug(msg)
                            continue

                        f = pressure * area * normal * scale
                        for inid in face:
                            inidi = nid_map[elem_node_ids[inid]]
                            nodal_pressures[inid] += pressure * scale / nface
                            forces[inidi, :] += f / nface
                        centroidal_pressures[ie] += pressure

                        #r = centroid - p
                        #load.cid.transformToGlobal()
                        #m = cross(r, f)
                        #M += m

            #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
            elif load.type == 'SPCD':
                #self.gids = [integer(card, 2, 'G1'),]
                #self.constraints = [components_or_blank(card, 3, 'C1', 0)]
                #self.enforced = [double_or_blank(card, 4, 'D1', 0.0)]
                for nid, c1, d1 in zip(load.node_ids, load.constraints, load.enforced):
                    if nid in dependents_nodes:
                        fail_nids.add(nid)
                        fail_count += 1
                        if fail_count < fail_count_max:
                            self.log.warning('    nid=%s is a dependent node and has an'
                                             ' SPCD applied\n%s' % (nid, str(load)))
                    c1 = int(c1)
                    assert c1 in [1, 2, 3, 4, 5, 6], c1
                    if c1 < 4:
                        spcd[nid_map[nid], c1 - 1] = d1
            elif load.type in ['FORCE1', 'MOMENT1']:
                pass
            else:
                print(load)
                if load.type not in cards_ignored:
                    cards_ignored[load.type] = True
                    self.log.warning('  _get_forces_moments_array - unsupported '
                                     'load.type = %s' % load.type)
        if fail_count:
            fail_nids_list = list(fail_nids)
            fail_nids_list.sort()
            self.log.warning('fail_nids = %s' % np.array(fail_nids_list))
        return centroidal_pressures, forces, spcd

    def get_pressure_array(self, load_case, eids, normals):
        """
        Gets the pressures for a load case

        Parameters
        ----------
        load_case : ???
        eids : ???
        normals : (nelements, 3) float ndarray
            the element normals

        Returns
        -------
        is_pressure : bool
            the pressure data
        pressures : (nelements, 1) float ndarray
            the centroidal pressures
        """
        if 'PLOAD4' not in self.card_count:
            return False, None

        # account for scale factors
        loads2 = []
        scale_factors2 = []
        for load in load_case:
            if load.type == 'LOAD':
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        pressures = np.zeros(len(self.elements), dtype='float32')

        iload = 0
        nloads = len(loads2)
        show_nloads = nloads > 5000
        # loop thru scaled loads and plot the pressure
        for load, scale in zip(loads2, scale_factors2):
            if show_nloads and iload % 5000 == 0:
                self.log.debug('  NastranIOv iload=%s/%s' % (iload, nloads))
            if load.type == 'PLOAD4':
                #print(load.object_attributes())
                for elem in load.eids:
                    #elem = self.elements[eid]
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
                                     'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        pressure = load.pressures[0] * scale

                        # single element per PLOAD
                        #eid = elem.eid
                        #pressures[eids.index(eid)] = pressure

                        # multiple elements
                        #for elem in load.eids:
                        ie = np.searchsorted(eids, elem.eid)
                        #pressures[ie] += p  # correct; we can't assume model orientation
                        normal = normals[ie, 2]  # considers normal of shell
                        pressures[ie] += pressure * normal

                    #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        #A, centroid, normal = elem.get_face_area_centroid_normal(
                            #load.g34.nid, load.g1.nid)
                        #r = centroid - p
            iload += 1
        return True, pressures

    def _get_temperatures_array(self, load_case_id, dtype='float32'):
        """
        Builds the temperature array based on thermal cards.
        Used by the GUI.

        Parameters
        ----------
        load_case_id : int
            the load id
        dtype : str; default='float32'
            the type of the temperature array

        Returns
        -------
        is_temperatures : bool
            is there temperature data
        temperatures : (nnodes, ) float ndarray
            the temperatures
        """
        if 'TEMP' not in self.card_count:
            return False, None
        is_temperatures = True
        nids = sorted(self.nodes.keys())

        load_case = self.loads[load_case_id]
        loads2, scale_factors2 = self._get_loads_and_scale_factors(load_case)
        tempd = self.tempds[load_case_id].temperature if load_case_id in self.tempds else 0.
        temperatures = np.ones(len(self.nodes), dtype=dtype) * tempd
        for load, scale in zip(loads2, scale_factors2):
            if load.type == 'TEMP':
                temps_dict = load.temperatures
                for nid, val in iteritems(temps_dict):
                    nidi = nids.index(nid)
                    temperatures[nidi] = val
            else:
                self.log.debug(load.type)
        return is_temperatures, temperatures

    def _get_rigid(self):
        """
        GUI helper function

        dependent = (lines[:, 0])
        independent = np.unique(lines[:, 1])
        """
        lines_rigid = []
        for eid, elem in iteritems(self.rigid_elements):
            if elem.type == 'RBE3':
                if elem.Gmi != []:
                    msg = 'UM is not supported; RBE3 eid=%s Gmi=%s' % (elem.eid, elem.Gmi)
                    raise RuntimeError(msg)
                #list_fields = ['RBE3', elem.eid, None, elem.ref_grid_id, elem.refc]
                n1 = elem.ref_grid_id
                assert isinstance(n1, integer_types), 'RBE3 eid=%s ref_grid_id=%s' % (elem.eid, n1)
                for (_weight, ci, Gij) in elem.WtCG_groups:
                    Giji = elem._nodeIDs(nodes=Gij, allow_empty_nodes=True)
                    # list_fields += [wt, ci] + Giji
                    for n2 in Giji:
                        assert isinstance(n2, integer_types), 'RBE3 eid=%s Giji=%s' % (elem.eid, Giji)
                        lines_rigid.append([n1, n2])
            elif elem.type == 'RBE2':
                #list_fields = ['RBE2', elem.eid, elem.Gn(), elem.cm
                               #] + elem.Gmi_node_ids + [elem.alpha]
                n2 = elem.Gn() # independent
                nids1 = elem.Gmi_node_ids # dependent
                for n1 in nids1:
                    lines_rigid.append([n1, n2])
            elif elem.type in ['RBAR', 'RBAR1', 'RROD']: ## TODO: these aren't quite right
                dependent = elem.Ga()
                independent = elem.Gb()
                lines_rigid.append([dependent, independent])
            else:
                print(str(elem))
        return lines_rigid

    def _get_loads_and_scale_factors(self, load_case):
        """account for scale factors"""
        loads2 = []
        scale_factors2 = []
        for load in load_case:
            if load.type == 'LOAD':
                scale_factors, loads = load.get_reduced_loads()
                for scale_factor, loadi in zip(scale_factors, loads):
                    if scale_factor == 0.0:
                        continue
                    scale_factors2.append(scale_factor)
                    loads2.append(loadi)
                #else:
                    #scale_factors2 += scale_factors
                    #loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)
        return loads2, scale_factors2

    def get_rigid_elements_with_node_ids(self, node_ids):
        """
        Gets the series of rigid elements that use specific nodes

        Parameters
        ----------
        node_ids : List[int]
            the node ids to check

        Returns
        -------
        rbes : List[int]
            the set of self.rigid_elements
        """
        try:
            nids = set(node_ids)
        except TypeError:
            print(node_ids)
            raise
        rbes = []
        for eid, rigid_element in iteritems(self.rigid_elements):
            if rigid_element.type in ['RBE3', 'RBE2', 'RBE1', 'RBAR', 'RSPLINE']:
                independent_nodes = set(rigid_element.independent_nodes)
                dependent_nodes = set(rigid_element.dependent_nodes)
                rbe_nids = independent_nodes | dependent_nodes
                if nids.intersection(rbe_nids):
                    rbes.append(eid)
            else:
                raise RuntimeError(rigid_element.type)
        return rbes

    def get_dependent_nid_to_components(self, mpc_id=None):
        """
        Gets a dictionary of the dependent node/components.

        Parameters
        ----------
        mpc_id : int; default=None -> no MPCs are checked
            TODO: add

        Returns
        -------
        dependent_nid_to_components : dict[node_id] : components
            node_id : int
                the node_id
            components : str
                the DOFs that are linked

        Nastran can either define a load/motion at a given node.
        SPCs define constraints that may not have loads/motions.

        MPCs and rigid elements define independent and dependent nodes on
        specific DOFs.
          - independent nodes : loads/motions may be defined
          - dependent nodes : loads/motions may not be defined
        """
        if mpc_id is not None:
            raise NotImplementedError('MPCs')
            #mpcs = self.get_mpcs(mpc_id)

        dependent_nid_to_components = {}
        for eid, rigid_element in iteritems(self.rigid_elements):
            if rigid_element.type == 'RBE2':
                dependent_nodes = set(rigid_element.dependent_nodes)
                components = rigid_element.cm
                for nid in dependent_nodes:
                    dependent_nid_to_components[nid] = components
            elif rigid_element.type == 'RBE3':
                dependent_nid_to_components[rigid_element.ref_grid_id] = rigid_element.refc
                for gmi, cmi in zip(rigid_element.Gmi_node_ids, rigid_element.Cmi):
                    dependent_nid_to_components[gmi] = cmi
            #if rigid_element.type in ['RBE3', 'RBE2', 'RBE1', 'RBAR']:
                ##independent_nodes = set(rigid_element.independent_nodes)
                #dependent_nodes = set(rigid_element.dependent_nodes)
                #rbe_nids = independent_nodes | dependent_nodes
                #if nids.intersection(rbe_nids):
                    #rbes.append(eid)
            #elif rigid_element == 'RSPLINE':
            else:
                raise RuntimeError(rigid_element.type)
        return dependent_nid_to_components

    #def get_node_ids_with_element(self, eid, msg=''):
        #self.deprecated('get_node_ids_with_element(eid)', 'get_node_ids_with_elements(eid)', '0.9')
        #return self.get_node_ids_with_elements(eid, msg=msg)

    def _get_maps(self, eids=None, map_names=None,
                  consider_0d=True, consider_0d_rigid=True,
                  consider_1d=True, consider_2d=True, consider_3d=True):
        """
        Gets a series of mappings (e.g. node_id to element_id)

        eids : List[int]
            the element ids to consider
        map_names : List[str]; default=None -> []
            does nothing
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

        .. todo:: map_names support
        .. todo:: consider_0d support
        .. todo:: consider_0d_rigid support
        """
        if map_names is None:
            map_names = []
        allowed_maps = [
            'edge_to_eid_map',
            'eid_to_edge_map',
            'nid_to_edge_map',
            'edge_to_nid_map',
            'eid_to_eid_map',
        ]

        for name in map_names:
            if name not in allowed_maps:
                msg = 'name=%s; allowed=%s' % (name, sorted(allowed_maps.keys()))
                raise RuntimeError(msg)

        eid_to_edge_map = {}
        eid_to_nid_map = {}

        edge_to_eid_map = defaultdict(set)
        #edge_to_nid_map = {}  # unnecessary

        nid_to_edge_map = defaultdict(set)  #set([]) ???
        nid_to_eid_map = defaultdict(list)


        if eids is None:
            eids = iterkeys(self.elements)

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
                nid_to_eid_map[nid].append(eid)
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
        out = (
            edge_to_eid_map,
            eid_to_edge_map,
            nid_to_edge_map,
        )
        return out

    def get_node_ids_with_elements(self, eids, msg=''):
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

        For example

        .. code-block:: python

          >>> eids = [1, 2, 3]  # list of elements with pid=1
          >>> msg = ' which are required for pid=1'
          >>> node_ids = bdf.get_node_ids_with_elements(eids, msg=msg)
        """
        if isinstance(eids, integer_types):
            eids = [eids]

        nids2 = set([])
        for eid in eids:
            element = self.Element(eid, msg=msg)
            self.log.debug("element.pid = %s" % (element.pid))
            nids = set(element.node_ids)
            nids2.update(nids)
        return nids2

    def get_elements_nodes_by_property_type(self, dtype='int32',
                                            save_element_types=False):
        """
        Gets a dictionary of (etype, pid) to [eids, node_ids]

        Parameters
        ----------
        dtype : str; default='int32'
            the type of the integers
        save_element_types : bool; default=False
            adds the etypes output

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
        etypes = self.get_element_nodes_by_element_type(dtype=dtype)
        output = {}
        for etype, (eids, pids, nids) in iteritems(etypes):
            upids = np.unique(pids)
            for upid in upids:
                ipid = np.where(pids == upid)[0]
                output[(etype, upid)] = [eids[ipid], nids[ipid, :]]
        if save_element_types:
            return output, None
        else:
            return output, etypes

    def get_element_nodes_by_element_type(self, dtype='int32', solids=None):
        """
        Gets a dictionary of element type to [eids, pids, node_ids]

        Parameters
        ----------
        dtype : str; default='int32'
            the type of the integers

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
                    fails if you're wrong
        """
        etypes = [
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
        ]
        output = {}

        if solids is None:
            solids = {
                'CTETRA' : (4, 10),
                #'CTETRA' : (10, ),
                'CHEXA' : (8, 20),
                'CPENTA' : (6, 15),
                'CPYRAM' : (5, 13),
            }
        for etype in etypes:
            if etype not in self._type_to_id_map:
                continue
            eids = np.array(self._type_to_id_map[etype], dtype=dtype)
            neids = len(eids)
            eid0 = eids[0]

            elem0 = self.elements[eid0]
            nnodes = len(elem0.nodes)

            if etype not in solids or len(solids[etype]) == 1:
                pids = np.zeros(neids, dtype=dtype)
                nids = np.zeros((neids, nnodes), dtype=dtype)
                for i, eid in enumerate(eids):
                    elem = self.elements[eid]
                    pid = elem.Pid()
                    assert pid is not None, elem
                    nidsi = elem.node_ids
                    #self.log.info(str(elem))
                    try:
                        nids[i, :] = nidsi
                    except TypeError:
                        #print(elem)
                        #print('nidsi =', nidsi)
                        nidsi2 = [nid  if nid is not None else 0
                                 for nid in nidsi]
                        try:
                            nids[i, :] = nidsi2
                        except:
                            print(elem)
                            print(nidsi)
                            print(nidsi2)
                            raise
                pids[i] = pid
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
                pids[i] = pid
                if len(ieids_max):
                    etype_max = elem.type + str(nnodes_max)
                    ieids_max = np.array(ieids_max, dtype=dtype)
                    output[etype_max] = [eids[ieids_max], pids[ieids_max], nids[ieids_max, :]]
                if len(ieids_min):
                    etype_min = elem.type + str(nnodes_min)
                    ieids_min = np.array(ieids_min, dtype=dtype)
                    output[etype_min] = [eids[ieids_min], pids[ieids_min], nids[ieids_min, :nnodes_min]]
        assert len(output), 'output is empty...'
        return output

    #--------------------
    # ELEMENT CARDS

    def get_element_ids_list_with_pids(self, pids=None):
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
        """
        if pids is None:
            pids = iterkeys(self.properties)
        elif isinstance(pids, integer_types):
            pids = [int]
        assert isinstance(pids, (list, tuple)), 'pids=%s type=%s' % (pids, type(pids))
        eids2 = []
        for eid, element in sorted(iteritems(self.elements)):
            pid = element.Pid()
            if pid in pids:
                eids2.append(eid)
        return eids2

    def get_pid_to_node_ids_and_elements_array(self, pids=None, etypes=None):
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
        node_ids : (nelements, nnodes)
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
        pid_to_eids_ieids_map = {}
        if etypes is None:
            etypes = etypes_
        for etype in etypes_:
            if etype not in etypes_:
                continue
            if etype in skip_elements:
                self.log.warning('skipping etype=%s' % etype)
                continue
            eids = self._type_to_id_map[etype]
            element0 = self.elements[eids[0]]
            nnodes = len(element0.node_ids)
            neids = len(eids)
            node_ids = np.zeros((neids, nnodes), dtype='int32')
            for ieid, eid in enumerate(eids):
                element = self.elements[eid]
                node_ids[ieid, :] = element.node_ids
                pid = element.Pid()
                #nids_to_pids_map[]
                pid_to_eids_ieids_map[(pid, etype)].append((eid, ieid))
            etype_to_nids_map[etype] = node_ids
        for key, value in iteritems(pid_to_eids_ieids_map):
            pid_to_eids_ieids_map[key] = np.array(value, dtype='int32')
        return pid_to_eids_ieids_map

    def get_element_ids_dict_with_pids(self, pids=None, stop_if_no_eids=True):
        """
        Gets all the element IDs with a specific property ID.

        Parameters
        ----------
        pids : List[int]
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

          # consider all properties
          eids_dict = model.get_element_ids_dict_with_pids()

        Note
        ----
        What happens with CONRODs?
        """
        if pids is None:
            pids = list(self.properties)
        elif isinstance(pids, integer_types):
            pids = [int]

        assert isinstance(pids, (list, tuple, np.ndarray)), 'pids=%s type=%s' % (pids, type(pids))
        pid_to_eids_map = {}
        for pid in pids:
            pid_to_eids_map[pid] = []

        for eid, element in iteritems(self.elements):
            try:
                pid = element.Pid()
                #if element.type == 'CONROD':
                    #raise RuntimeError('CONROD pid=%r' % pid)
                if pid in pids:
                    pid_to_eids_map[pid].append(eid)
            except AttributeError:
                #eids2[0].append(eid)
                pass

        if stop_if_no_eids:
            for eids in itervalues(pid_to_eids_map):
                if len(eids):
                    return pid_to_eids_map
            raise RuntimeError('no elements found')
        return pid_to_eids_map

    def get_node_id_to_element_ids_map(self):
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
            for nid in sorted(self.spoints.points):  # SPOINTs
                nid_to_eids_map[nid] = []

        for (eid, element) in iteritems(self.elements):  # load the mapper
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

    def get_node_id_to_elements_map(self):
        """
        Returns a dictionary that maps node IDs to a list of elemnents

        .. todo:: support 0d or 1d elements
        .. todo:: support elements with missing nodes
                  (e.g. CQUAD8 with missing nodes)
        """
        nid_to_elements_map = {}
        for nid in self.nodes:  # initalize the mapper
            nid_to_elements_map[nid] = []

        if self.spoints:  # SPOINTs
            for nid in sorted(self.spoints.spoints):  # SPOINTs
                nid_to_elements_map[nid] = []

        for (eid, element) in iteritems(self.elements):  # load the mapper
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

    def get_property_id_to_element_ids_map(self):
        """
        Returns a dictionary that maps a property ID to a list of elemnents
        """
        pid_to_eids_map = {}
        pids = self.property_ids
        for pid in pids:
            pid_to_eids_map[pid] = []

        for eid in self.element_ids:
            element = self.Element(eid)
            element_type = element.type
            if element_type in ['CONROD', 'CONM2', 'CELAS2', 'CELAS4', 'CDAMP2', 'CDAMP4']:
                continue
            if hasattr(element, 'pid'):
                pid = element.Pid()
                try:
                    pid_to_eids_map[pid].append(eid)
                except KeyError:
                    print(element)
                    raise KeyError('pid=%s is invalid for card=\n%s' % (pid, str(element)))
        return pid_to_eids_map

    def get_material_id_to_property_ids_map(self):
        """
        Returns a dictionary that maps a material ID to a list of properties

        Returns
        -------
        mid_to_pids_map : dict[int] = int
            the mapping

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
        for mid in self.get_material_ids():
            mid_to_pids_map[mid] = []

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
            elif prop_type in ['PGAP', 'PELAS', 'PVISC', 'PBUSH', 'PDAMP', 'PFAST']:
                pass
            elif prop_type in ['PSHELL']:
                mids = prop.material_ids
                for i, mid in enumerate(mids):
                    if mid is None:
                        continue
                    try:
                        mid_to_pids_map[mid].append(pid)
                    except KeyError:
                        print(prop)
                        raise KeyError('i=%s mid=%s is invalid for card=\n%s' % (i, mid, str(prop)))
            else:
                mid = prop.Mid()
                try:
                    mid_to_pids_map[mid].append(pid)
                except KeyError:
                    print(prop)
                    raise KeyError('mid=%s is invalid for card=\n%s' % (mid, str(prop)))
        return mid_to_pids_map

    def get_reduced_mpcs(self, mpc_id):
        """get all MPCs that are part of a set"""
        mpcs = self.MPC(mpc_id)
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
                            mpcs2i = self.get_reduced_mpcs(mpciii)
                            mpcs2 += mpcs2i
                    else:
                        assert isinstance(mpci, integer_types), mpci
                        mpcs2i = self.get_reduced_mpcs(mpci)
                        mpcs2 += mpcs2i
            else:
                mpcs2.append(mpc)
        return mpcs2

    def get_reduced_spcs(self, spc_id):
        """get all SPCs that are part of a set"""
        spcs = self.SPC(spc_id)
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
                            spcs2i = self.get_reduced_spcs(spciii)
                            spcs2 += spcs2i
                    else:
                        # print('spci =', spci)
                        # print(spci.object_attributes())
                        assert isinstance(spci, integer_types), spci
                        spcs2i = self.get_reduced_spcs(spci)
                        spcs2 += spcs2i
            else:
                spcs2.append(spc)
        return spcs2

    def get_spcs(self, spc_id, consider_nodes=False):
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

        Doesn't consider:
          - non-zero enforced value on SPC
        """
        spcs = self.get_reduced_spcs(spc_id)
        #self.spcs[key] = [constraint]
        nids = []
        comps = []
        for spc in spcs:
            if spc.type == 'SPC1':
                nodes = spc.nodes
                nnodes = len(nodes)
                nids += nodes
                comps += [str(spc.constraints)] * nnodes
            elif spc.type == 'SPC':
                for nid, comp, enforced in zip(spc.gids, spc.constraints, spc.enforced):
                    nids.append(nid)
                    comps.append(comp)
            else:
                self.log.warning('not considering:\n%s' % str(spc))
                #raise NotImplementedError(spc.type)

        if consider_nodes:
            for nid, node in iteritems(self.nodes):
                if node.ps:
                    nids.append(nid)
                    comps.append(node.ps)
        return nids, comps

    def get_mpcs(self, mpc_id):
        """
        Gets the MPCs in a semi-usable form.

        Parameters
        ----------
        mpc_id : int
            the desired MPC ID

        Returns
        -------
        nids : List[int]
            the constrained nodes
        comps : List[str]
            the components that are constrained on each node

        Considers:
          - MPC
          - MPCADD
        """
        mpcs = self.get_reduced_mpcs(mpc_id)
        nids = []
        comps = []
        for mpc in mpcs:
            if mpc.type == 'MPC':
                for nid, comp, enforced in zip(mpc.gids, mpc.constraints, mpc.enforced):
                    nids.append(nid)
                    comps.append(comp)
            else:
                self.log.warning('not considering:\n%s' % str(mpc))
                #raise NotImplementedError(mpc.type)
        return nids, comps
