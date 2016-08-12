# pylint: disable=E1101,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, iteritems, iterkeys
from collections import defaultdict

import numpy as np
from pyNastran.utils import integer_types
from pyNastran.bdf.deprecated import GetMethodsDeprecated
from pyNastran.bdf.cards.nodes import SPOINT, EPOINT
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes


class GetMethods(GetMethodsDeprecated, BDFAttributes):
    def __init__(self):
        self._type_to_slot_map = {}
        BDFAttributes.__init__(self)

    def get_card_ids_by_card_types(self, card_types=None, reset_type_to_slot_map=False,
                                   stop_on_missing_card=False):
        """
        Parameters
        ----------
        card_types : List[str]
            the list of keys to consider (list of strings; string)
        reset_type_to_slot_map : bool
            should the mapping dictionary be rebuilt (default=False);
            set to True if you added cards
        stop_on_missing_card : bool
            crashes if you request a card and it doesn't exist

        Returns
        -------
        out_dict: dict[str]=List[ids]
            the key=card_type, value=the ID of the card object
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
        #out = {
            #(key) : (self._type_to_id_map[key] if key in self.card_count else [])
            #for key in card_types
        #}
        out = {}
        for key in card_types:
            if key in self.card_count:
                out[key] = sorted(self._type_to_id_map[key])
            else:
                if stop_on_missing_card:
                    raise RuntimeError('%r is not in the card_count; keys=%s' %
                                       str(sorted(self.card_count.keys())))
                out[key] = []
        return out

    def _reset_type_to_slot_map(self):
        rslot_map = defaultdict(list)
        for dict_name, card_names in iteritems(self._slot_to_type_map):
            print('card_names=%s dict_name=%s' % (card_names, dict_name))
            card_name0 = card_names[0]
            if card_name0 in ['DTABLE', 'GRDSET', 'SESUP', 'DOPTPRM', 'MONPNT1', 'SUPORT', 'MKAERO1',
                              'MATHP']:
                pass
            else:
                adict = getattr(self, dict_name)
                if isinstance(adict, dict):
                    for key, card in iteritems(adict):
                        if isinstance(card, list):
                            alist = card
                            for cardi in alist:
                                rslot_map[cardi.type].append(key)
                            #raise NotImplementedError('%s; names=%s \ncard=%s' % (type(card), card_names, card))
                        else:
                            rslot_map[card.type].append(key)
                elif isinstance(adict, list):
                    alist = adict
                    for value in alist:
                        if isinstance(value, list):
                            raise NotImplementedError('%s; names=%s value=%s' % (type(value), card_names, value))
                        else:
                            if value.type in ['CSET1', 'CSET']:
                                pass
                                #rslot_map[value.type] = value.
                            else:
                                raise NotImplementedError('list; names=%s' % card_names)
                else:
                    raise NotImplementedError('%s; names=%s' % (type(adict), card_names))
        rslot_map
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
        if reset_type_to_slot_map or self._type_to_slot_map is None:
            rslot_map = {}
            for key, values in iteritems(self._slot_to_type_map):
                for value in values:
                    rslot_map[value] = key
            self._type_to_slot_map = rslot_map
        else:
            rslot_map = self._type_to_slot_map

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
                key = rslot_map[card_type]
            except:
                self.log.error("card_type=%r' hasn't been added to self._slot_to_type_map...check for typos")
                raise
            slot = getattr(self, key)
            ids = self._type_to_id_map[card_type]
            cards = []
            for idi in ids:
                try:
                    card = slot[idi]
                except KeyError:
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

    def get_rigid_elements_with_node_ids(self, node_ids):
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

    def get_node_ids_with_element(self, eid, msg=''):
        return self.get_node_ids_with_elements([eid], msg=msg)

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
        consider_2d : bool; default=True
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
            assert name in allowed_maps, 'name=%s; allowed=%s' % (name, sorted(allowed_maps.keys()))

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
                assert not isinstance(edge, int), 'edge=%s elem=\n%s' % (edge, elem)
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
        nids2 = set([])
        for eid in eids:
            element = self.Element(eid, msg=msg)
            self.log.debug("element.pid = %s" % (element.pid))
            nids = set(element.node_ids)
            nids2.update(nids)
        return nids2

    def Node(self, nid, allow_empty_nodes=False, msg=''):
        if (nid == 0 or nid is None) and allow_empty_nodes:
            return None
        elif nid in self.nodes:
            return self.nodes[nid]
        elif self.spoints and nid in self.spoints.points:
            return SPOINT(nid)
        elif self.epoints and nid in self.epoints.points:
            return EPOINT(nid)
        else:
            assert isinstance(nid, integer_types), 'nid should be an integer; not %s' % type(nid)
            nid_list = np.unique(list(self.nodes.keys()))
            raise RuntimeError('nid=%s is not a GRID, SPOINT, or EPOINT%s\n%s' % (nid, msg, nid_list))

    def Nodes(self, nids, allow_empty_nodes=False, msg=''):
        """
        Returns a series of node objects given a list of node IDs
        """
        nodes = []
        for nid in nids:
            nodes.append(self.Node(nid, allow_empty_nodes=allow_empty_nodes, msg=msg))
        return nodes

    #--------------------
    # ELEMENT CARDS

    def get_element_ids_list_with_pids(self, pids):
        """
        Gets all the element IDs with a specific property ID.

        Parameters
        ----------
        pids : List[int]
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
        assert isinstance(pids, (list, tuple)), pids
        eids2 = []
        for eid, element in sorted(iteritems(self.elements)):
            pid = element.Pid()
            if pid in pids:
                eids2.append(eid)
        return eids2

    def get_element_ids_dict_with_pids(self, pids=None):
        """
        Gets all the element IDs with a specific property ID.

        Parameters
        ----------
        pids : List[int]
            list of property ID

        Returns
        -------
        element_ids : dict[int] = List[int]
            as a dictionary of lists by property

        For example, we want all the elements with ``pids=[4, 5, 6]``,
        but we want them in separate groups

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          pids = [4, 5, 6]
          eids_dict = model.get_element_ids_dict_with_pids(pids)

          # consider all properties
          eids_dict = model.get_element_ids_dict_with_pids()
        """
        if pids is None:
            pids = self.properties.keys()
        elif isinstance(pids, int):
            pids = [int]
        assert isinstance(pids, (list, tuple)), pids
        eids2 = {}
        for pid in pids:
            eids2[pid] = []
        for eid, element in iteritems(self.elements):
            try:
                pid = element.Pid()
                if pid in pids:
                    eids2[pid].append(eid)
            except AttributeError:
                #eids2[0].append(eid)
                pass
        return eids2

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
            for nid in sorted(self.spoints.spoints):  # SPOINTs
                nid_to_eids_map[nid] = []

        for (eid, element) in iteritems(self.elements):  # load the mapper
            try:
                # not supported for 0-D and 1-D elements
                nids = element.node_ids
            except AttributeError:
                print(element.type)
            else:
                for nid in nids:  # (e.g. CQUAD8 with missing node)
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

    def Element(self, eid, msg=''):
        """gets an element (not rigid (RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE) or mass (CMASS1, CONM2))"""
        try:
            return self.elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed elements=%s'
                           % (eid, msg, np.unique(list(self.elements.keys()))))

    def Elements(self, eids, msg=''):
        elements = []
        for eid in eids:
            elements.append(self.Element(eid, msg))
        return elements

    def Mass(self, eid, msg=''):
        """gets a mass element (CMASS1, CONM2)"""
        try:
            return self.masses[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed masses=%s'
                           % (eid, msg, np.unique(list(self.masses.keys()))))

    def RigidElement(self, eid, msg=''):
        """gets a rigid element (RBAR, RBE2, RBE3, RBAR, RBAR1, RSPLINE)"""
        try:
            return self.rigid_elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed rigid_elements=%s'
                           % (eid, msg, np.unique(list(self.rigid_elements.keys()))))

    #--------------------
    # PROPERTY CARDS

    def Property(self, pid, msg=''):
        """
        gets an elemental property (e.g. PSOLID, PLSOLID, PCOMP, PSHELL, PSHEAR);
        not mass property (PMASS)
        """
        try:
            return self.properties[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Pids=%s'
                           % (pid, msg, self.property_ids))

    def Properties(self, pids, msg=''):
        properties = []
        for pid in pids:
            properties.append(self.Property(pid, msg))
        return properties

    def PropertyMass(self, pid, msg=''):
        """
        gets a mass property (PMASS)
        """
        try:
            return self.properties_mass[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Mass Pids=%s'
                           % (pid, msg, np.unique(list(self.mass_property.keys()))))

    def Phbdy(self, pid, msg=''):
        """gets a PHBDY"""
        try:
            return self.phbdys[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PHBDY Pids=%s'
                           % (pid, msg, np.unique(list(self.phbdys.keys()))))

    #--------------------
    # MATERIAL CARDS

    def get_structural_material_ids(self):
        return self.materials.keys()

    def get_material_ids(self):
        return self.materials.keys() + self.thermal_materials.keys() + self.hyperelastic_materials.keys()

    def get_thermal_material_ids(self):
        return self.thermal_materials.keys()

    def Material(self, mid, msg=''):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermal_materials:
            return self.thermal_materials[mid]
        else:
            msg = '\n' + msg
            raise KeyError('Invalid Material ID:  mid=%s%s' % (mid, msg))

    def StructuralMaterial(self, mid, msg=''):
        try:
            mat = self.materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Structural Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def ThermalMaterial(self, mid, msg=''):
        try:
            mat = self.thermal_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Thermal Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def HyperelasticMaterial(self, mid, msg=''):
        try:
            mat = self.hyperelastic_materials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Hyperelastic Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def Materials(self, mids, msg=''):
        materials = []
        for mid in mids:
            materials.append(self.Material(mid, msg))
        return materials

    #--------------------
    # LOADS

    def Load(self, sid, msg=''):
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid
        if sid in self.loads:
            load = self.loads[sid]
        else:
            raise KeyError('cannot find LoadID=%r%s.\nLoadIDs=%s\n' % (
                sid, msg, np.unique(list(self.loads.keys()))))
        return load

    def DLoad(self, sid, msg=''):
        """gets a DLOAD"""
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid
        if sid in self.dloads:
            load = self.dloads[sid]
        else:
            raise KeyError('cannot find DLoadID=%r%s.\nDLoadIDs=%s\n' % (
                sid, msg, np.unique(list(self.dloads.keys()))))
        return load

    def get_dload_entries(self, sid, msg=''):
        assert isinstance(sid, integer_types), 'sid=%s is not an integer\n' % sid
        if sid in self.dload_entries:
            load = self.dload_entries[sid]
        else:
            raise KeyError('cannot find DLoad Entry ID=%r%s.\nDLoadEntryIDs=%s\n' % (
                sid, msg, np.unique(list(self.dload_entries.keys()))))
        return load

    def DELAY(self, delay_id, msg=''):
        """gets a DELAY"""
        try:
            return self.delays[delay_id]
        except KeyError:
            raise KeyError('delay_id=%s not found%s.  Allowed DELAY=%s'
                           % (delay_id, msg, np.unique(list(self.delays.keys()))))

    def DPHASE(self, dphase_id, msg=''):
        """gets a DPHASE"""
        try:
            return self.dphases[dphase_id]
        except KeyError:
            raise KeyError('dphase_id=%s not found%s.  Allowed DPHASE=%s'
                           % (dphase_id, msg, np.unique(list(self.dphases.keys()))))

    #--------------------
    def MPC(self, mpc_id, msg=''):
        """gets an MPC"""
        assert isinstance(mpc_id, integer_types), 'mpc_id=%s is not an integer\n' % mpc_id
        if mpc_id in self.mpcs:
            constraint = self.mpcs[mpc_id]
        else:
            raise KeyError('cannot find MPC ID=%r%s.\nAllowed MPCs=%s' % (
                mpc_id, msg, np.unique(list(self.mpcs.keys()))))
        return constraint

    def SPC(self, spc_id, msg=''):
        """gets an SPC"""
        assert isinstance(spc_id, integer_types), 'spc_id=%s is not an integer\n' % spc_id
        if spc_id in self.spcs:
            constraint = self.spcs[spc_id]
        else:
            raise KeyError('cannot find SPC ID=%r%s.\nAllowed SPCs=%s' % (
                spc_id, msg, np.unique(list(self.spcs.keys()))))
        return constraint

    def get_reduced_mpcs(self, mpc_id):
        """get all MPCs that are part of a set"""
        mpcs = self.MPC(mpc_id)
        mpcs2 = []
        for mpc in mpcs:
            if mpc.type == 'MPCADD':
                for mpci in mpc.sets:
                    if isinstance(mpci, list):
                        for mpcii in mpci:
                            if isinstance(mpcii, int):
                                mpciii = mpcii
                            else:
                                mpciii = mpcii.conid
                            mpcs2i = self.get_reduced_mpcs(mpciii)
                            mpcs2 += mpcs2i
                    else:
                        assert isinstance(mpci, int), mpci
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
                            if isinstance(spcii, int):
                                spciii = spcii
                            else:
                                spciii = spcii.conid
                            spcs2i = self.get_reduced_spcs(spciii)
                            spcs2 += spcs2i
                    else:
                        # print('spci =', spci)
                        # print(spci.object_attributes())
                        assert isinstance(spci, int), spci
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

        if consider_nodes and final_flag:
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

    #--------------------
    # Sets
    def SET1(self, set_id, msg=''):
        """gets a SET1"""
        assert isinstance(set_id, integer_types), 'set_id=%s is not an integer\n' % set_id
        if set_id in self.sets:
            set1 = self.sets[set_id]
        else:
            raise KeyError('cannot find SET1 ID=%r.\n%s' % (set_id, msg))
        return set1

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid, msg=''):
        """gets an COORDx"""
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coord_ids))

    #--------------------
    # AERO CARDS

    def AEList(self, aelist, msg=''):
        """gets an AELIST"""
        try:
            return self.aelists[aelist]
        except KeyError:
            raise KeyError('aelist=%s not found%s.  Allowed AELIST=%s'
                           % (aelist, msg, np.unique(list(self.aelists.keys()))))

    def AEFact(self, aefact, msg=''):
        """gets an AEFACT"""
        try:
            return self.aefacts[aefact]
        except KeyError:
            raise KeyError('aefact=%s not found%s.  Allowed AEFACT=%s'
                           % (aefact, msg, np.unique(list(self.aefacts.keys()))))

    def Acsid(self, msg=''):
        """gets the aerodynamic system coordinate"""
        if self.aero is not None:
            acsid_aero = self.aero.Acsid()
        if self.aeros is not None:
            acsid_aeros = self.aeros.Acsid()

        if self.aero is not None:
            if self.aeros is not None:
                assert acsid_aero == acsid_aeros, 'AERO acsid=%s, AEROS acsid=%s' % (acsid_aero,
                                                                                     acsid_aeros)
            cid = self.Coord(acsid_aero, msg=msg)
        elif self.aeros is not None:
            cid = self.Coord(acsid_aeros, msg=msg)
        else:
            msg = 'neither AERO nor AEROS cards exist.'
            raise RuntimeError(msg)
        return cid

    def Aero(self, msg=''):
        """gets the AERO"""
        if self.aero is not None:
            return self.aero
        else:
            raise RuntimeError('no AERO card found%s.' % (msg))

    def Aeros(self, msg=''):
        """gets the AEROS"""
        if self.aeros is not None:
            return self.aeros
        else:
            raise RuntimeError('no AEROS card found%s.' % (msg))

    def Spline(self, eid, msg=''):
        """gets a SPLINEx"""
        try:
            return self.splines[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, np.unique(list(self.splines.keys()))))

    def CAero(self, eid, msg=''):
        """gets an CAEROx"""
        try:
            return self.caeros[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, np.unique(list(self.caero_ids))))

    def PAero(self, pid, msg=''):
        """gets a PAEROx"""
        try:
            return self.paeros[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, np.unique(list(self.paeros.keys()))))

    def Gust(self, sid, msg=''):
        """gets a GUST"""
        try:
            return self.gusts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, np.unique(list(self.gusts.keys()))))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid, msg=''):
        """gets an AESTAT"""
        try:
            return self.aestats[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, np.unique(list(self.aestats.keys()))))

    def AELIST(self, aid, msg=''):
        """gets an AELIST"""
        try:
            return self.aelists[aid]
        except KeyError:
            raise KeyError('id=%s not found%s.  Allowed AELISTs=%s'
                           % (aid, msg, np.unique(list(self.aelists.keys()))))

    def AELink(self, link_id, msg=''):
        """gets an AELINK"""
        try:
            return self.aelinks[link_id]
        except KeyError:
            raise KeyError('link_id=%s not found%s.  Allowed AELINKs=%s'
                           % (link_id, msg, np.unique(list(self.aelinks.keys()))))

    def AEParam(self, aid, msg=''):
        """gets an AEPARM"""
        try:
            return self.aeparams[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, np.unique(list(self.aeparams.keys()))))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid, msg=''):
        """gets an FLFACT"""
        try:
            return self.flfacts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, np.unique(list(self.flfacts.keys()))))

    def Flutter(self, fid, msg=''):
        """gets a FLUTTER"""
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, np.unique(list(self.flutters.keys()))))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid, msg=''):
        """gets a DCONSTR"""
        try:
            return self.dconstrs[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, np.unique(list(self.dconstrs.keys()))))

    def DResp(self, dresp_id, msg=''):
        """gets a DRESPx"""
        try:
            return self.dresps[dresp_id]
        except KeyError:
            raise KeyError('dresp_id=%s not found%s.  Allowed DRESPx=%s'
                           % (dresp_id, msg, np.unique(list(self.dresps.keys()))))

    def Desvar(self, desvar_id, msg=''):
        """gets a DESVAR"""
        try:
            return self.desvars[desvar_id]
        except KeyError:
            raise KeyError('desvar_id=%s not found%s.  Allowed DESVARs=%s'
                           % (desvar_id, msg, np.unique(list(self.desvars.keys()))))

    def DDVal(self, oid, msg=''):
        """gets a DDVAL"""
        try:
            return self.ddvals[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, np.unique(list(self.ddvals.keys()))))

    def DVcrel(self, dv_id, msg=''):
        """gets a DVCREL1/DVCREL2"""
        try:
            return self.dvcrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVCRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvcrels.keys()))))

    def DVmrel(self, dv_id, msg=''):
        """gets a DVMREL1/DVMREL2"""
        try:
            return self.dvmrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVMRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvmrels.keys()))))

    def DVprel(self, dv_id, msg=''):
        """gets a DVPREL1/DVPREL2"""
        try:
            return self.dvprels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVPRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvprels.keys()))))

    #--------------------
    # SET CARDS

    def Set(self, sid, msg=''):
        try:
            return self.sets[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, np.unique(list(self.sets.keys()))))

    def SetSuper(self, seid, msg=''):
        try:
            return self.setsSuper[seid]
        except KeyError:
            raise KeyError('seid=%s not found%s.  Allowed SETx=%s'
                           % (seid, msg, np.unique(list(self.setsSuper.keys()))))

    #--------------------
    # METHOD CARDS
    def Method(self, sid, msg=''):
        """gets a METHOD (EIGR, EIGRL)"""
        try:
            return self.methods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, np.unique(list(self.methods.keys()))))

    def CMethod(self, sid, msg=''):
        """gets a METHOD (EIGC)"""
        try:
            return self.cmethods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, np.unique(list(self.cmethods.keys()))))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg=''):
        """gets a TABLEx (TABLED1, TABLED2, TABLD3)"""
        try:
            return self.tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, np.unique(list(self.tables.keys()))))

    def RandomTable(self, tid, msg=''):
        try:
            return self.random_tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, np.unique(list(self.random_tables.keys()))))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid, msg=''):
        """gets an NLPARM"""
        try:
            return self.nlparms[nid]
        except KeyError:
            raise KeyError('nid=%s not found%s.  Allowed NLPARMs=%s'
                           % (nid, msg, np.unique(list(self.nlparms.keys()))))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname, msg=''):
        """gets a DMIG"""
        try:
            return self.dmig[dname]
        except KeyError:
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, np.unique(list(self.dmig.keys()))))

    def DEQATN(self, equation_id, msg=''):
        """gets a DEQATN"""
        try:
            return self.dequations[equation_id]
        except KeyError:
            raise KeyError('equation_id=%s not found%s.  Allowed DMIGs=%s'
                           % (equation_id, msg, np.unique(list(self.dequations.keys()))))

