# pylint: disable=E1101,C0103,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, iteritems
#import sys
from numpy import ndarray

from pyNastran.bdf.deprecated import GetMethodsDeprecated
from pyNastran.bdf.cards.nodes import SPOINT


class GetMethods(GetMethodsDeprecated):
    def __init__(self):
        pass

    def get_node_ids_with_element(self, eid, msg=''):
        return self.get_node_ids_with_elements([eid], msg=msg)

    def get_node_ids_with_elements(self, eids, msg=''):
        """
        Get the node IDs associated with a list of element IDs

        :param self: the BDF object
        :param eids: list of element ID
        :param msg:  a additional message to print out if an element is not found
        :returns node_ids: set of node IDs

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
            nids = set(element.nodeIDs())
            nids2.update(nids)
        return nids2

    def Node(self, nid, allowEmptyNodes=False, msg=''):
        if (nid == 0 or nid is None) and allowEmptyNodes:
            return None
        elif nid in self.nodes:
            return self.nodes[nid]
        elif self.spoints and nid in self.spoints.spoints:
            return SPOINT(nid)
        else:
            assert isinstance(nid, int), 'nid should be an integer; not %s' % type(nid)
            raise RuntimeError('nid=%s is not a GRID or SPOINT%s' % (nid, msg))

    def Nodes(self, nids, allowEmptyNodes=False, msg=''):
        """
        Returns a series of node objects given a list of node IDs
        """
        nodes = []
        for nid in nids:
            nodes.append(self.Node(nid, allowEmptyNodes, msg))
        return nodes

    #--------------------
    # ELEMENT CARDS

    def get_element_ids_list_with_pids(self, pids):
        """
        Gets all the element IDs with a specific property ID.

        :param self: the BDF object
        :param pids: list of property ID
        :returns element_ids: as a list

        For example, we want to get all the element ids with ``pids=[1, 2, 3]``

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          pids = [1, 2, 3]
          eids_list = model.get_element_ids_with_pids(pids, mode='list')
        """
        assert isinstance(pids, list), pids
        eids2 = []
        for eid, element in sorted(iteritems(self.elements)):
            pid = element.Pid()
            if pid in pids:
                eids2.append(eid)
        return eids2

    def get_element_ids_dict_with_pids(self, pids):
        """
        Gets all the element IDs with a specific property ID.

        :param self: the BDF object
        :param pids: list of property ID
        :returns element_ids: as a dictionary of lists by property

        For example, we want all the elements with ``pids=[4, 5, 6]``,
        but we want them in separate groups

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          pids = [4, 5, 6]
          eids_dict = model.get_element_ids_with_pids(pids, mode='dict')
        """
        assert isinstance(pids, list), pids
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
        Returns a dictionary that maps a node ID to a list of elemnents

        .. todo:: support 0d or 1d elements
        .. todo:: support elements with missing nodes
                  (e.g. CQUAD8 with missing nodes)
        """
        nidToElementsMap = {}
        for nid in self.nodes:  # initalize the mapper
            nidToElementsMap[nid] = []

        if self.spoints:  # SPOINTs
            for nid in sorted(self.spoints.spoints):  # SPOINTs
                nidToElementsMap[nid] = []

        for (eid, element) in iteritems(self.elements):  # load the mapper
            try:
                # not supported for 0-D and 1-D elements
                nids = element.nodeIDs()
                for nid in nids:  # (e.g. CQUAD8 with missing node)
                    nidToElementsMap[nid].append(eid)
            except:
                pass
        return nidToElementsMap

    def get_property_id_to_element_ids_map(self):
        """
        Returns a dictionary that maps a property ID to a list of elemnents
        """
        pidToEidsMap = {}
        pids = self.property_ids
        for pid in pids:
            pidToEidsMap[pid] = []

        for eid in self.element_ids:
            element = self.Element(eid)

            if hasattr(element, 'pid'):
                pid = element.Pid()
                if pid == 0:  # CONM2
                    continue
                pidToEidsMap[pid].append(eid)
        return pidToEidsMap

    def get_material_id_to_property_ids_map(self):
        """
        Returns a dictionary that maps a material ID to a list of properties

        :returns mid_to_pids_map: the mapping

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
            if prop.type in ['PCOMP', 'PCOMPG']:
                mids = prop.Mids()

                for mid in mids:
                    if pid not in mid_to_pids_map[mid]:
                        mid_to_pids_map[mid].append(pid)
                    else:  # PCOMP
                        if hasattr(prop, 'mid') and prop.Mid() in mids:
                            if pid not in mid_to_pids_map[mid]:
                                mid_to_pids_map[mid].append(pid)
            else:
                if hasattr(prop, 'Mids'):
                    mids = prop.Mids()
                else:
                    mids = [prop.Mid()]

                for mid in mids:
                    mid_to_pids_map[mid].append(pid)
        return mid_to_pids_map

    def Element(self, eid, msg=''):
        try:
            return self.elements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed elements=%s'
                           % (eid, msg, self.elements.keys()))

    def Elements(self, eids, msg=''):
        elements = []
        for eid in eids:
            elements.append(self.Element(eid, msg))
        return elements

    def Mass(self, eid, msg=''):
        try:
            return self.masses[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed masses=%s'
                           % (eid, msg, self.masses.keys()))

    def RigidElement(self, eid, msg=''):
        try:
            return self.rigidElements[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed rigidElements=%s'
                           % (eid, msg, self.rigidElements.keys()))

    #--------------------
    # PROPERTY CARDS

    def Property(self, pid, msg=''):
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
        try:
            return self.properties_mass[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed Mass Pids=%s'
                           % (pid, msg, self.mass_property.keys()))

    def Phbdy(self, pid, msg=''):
        try:
            return self.phbdys[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PHBDY Pids=%s'
                           % (pid, msg, self.phbdys.keys()))

    #--------------------
    # MATERIAL CARDS

    def get_structural_material_ids(self):
        return self.materials.keys()

    def get_material_ids(self):
        return self.materials.keys() + self.thermalMaterials.keys()

    def get_thermal_material_ids(self):
        return self.thermalMaterials.keys()

    def Material(self, mid, msg=''):
        if mid in self.materials:
            return self.materials[mid]
        elif mid in self.thermalMaterials:
            return self.thermalMaterials[mid]
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
            mat = self.thermalMaterials[mid]
        except KeyError:
            msg = '\n' + msg
            raise KeyError('Invalid Thermal Material ID:  mid=%s%s' % (mid, msg))
        return mat

    def HyperelasticMaterial(self, mid, msg=''):
        try:
            mat = self.hyperelasticMaterials[mid]
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
        assert isinstance(sid, int), 'sid=%s is not an integer\n' % sid
        if sid in self.loads:
            load = self.loads[sid]
        else:
            raise KeyError('cannot find LoadID=%r%s.\nLoadIDs=%s\n' % (sid, msg, sorted(self.loads.keys())))
        return load

    def DLoad(self, sid, msg=''):
        assert isinstance(sid, int), 'sid=%s is not an integer\n' % sid
        if sid in self.dloads:
            load = self.dloads[sid]
        else:
            raise KeyError('cannot find DLoadID=%r%s.\nDLoadIDs=%s\n' % (sid, msg, sorted(self.dloads.keys())))
        return load

    def get_dload_entries(self, sid, msg=''):
        assert isinstance(sid, int), 'sid=%s is not an integer\n' % sid
        if sid in self.dload_entries:
            load = self.dload_entries[sid]
        else:
            raise KeyError('cannot find DLoad Entry ID=%r%s.\nDLoadEntryIDs=%s\n' % (sid, msg, sorted(self.dload_entries.keys())))
        return load

    #--------------------
    # SPCs

    def SPC(self, conid, msg=''):
        assert isinstance(conid, int), 'conid=%s is not an integer\n' % conid
        if conid in self.spcs:
            constraint = self.spcs[conid]
        else:
            raise KeyError('cannot find ConstraintID=%r.\n%s' % (conid, msg))
        return constraint

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid, msg=''):
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coord_ids))

    #--------------------
    # AERO CARDS

    def AEList(self, aelist, msg=''):
        try:
            return self.aelists[aelist]
        except KeyError:
            raise KeyError('aelist=%s not found%s.  Allowed AELIST=%s'
                           % (aelist, msg, self.aelists.keys()))

    def AEFact(self, aefact, msg=''):
        try:
            return self.aefacts[aefact]
        except KeyError:
            raise KeyError('aefact=%s not found%s.  Allowed AEFACT=%s'
                           % (aefact, msg, self.aefacts.keys()))

    def Aero(self, acsid, msg=''):
        try:
            return self.aero[acsid]
        except KeyError:
            raise KeyError('acsid=%s not found%s.  Allowed AERO=%s'
                           % (acsid, msg, self.aero.keys()))

    def Aeros(self, acsid, msg=''):
        try:
            return self.aeros[acsid]
        except KeyError:
            raise KeyError('acsid=%s not found%s.  Allowed AEROS=%s'
                           % (acsid, msg, self.aeros.keys()))

    def Spline(self, eid, msg=''):
        try:
            return self.splines[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, self.splines.keys()))

    def CAero(self, eid, msg=''):
        try:
            return self.caeros[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, self.caero_ids))

    def PAero(self, pid, msg=''):
        try:
            return self.paeros[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, self.paeros.keys()))

    def Gust(self, sid, msg=''):
        try:
            return self.gusts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, self.gusts.keys()))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid, msg=''):
        try:
            return self.aestats[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, self.aestats.keys()))

    def AELink(self, linkID, msg=''):
        try:
            return self.aelinks[linkID]
        except KeyError:
            raise KeyError('linkID=%s not found%s.  Allowed AELINKs=%s'
                           % (linkID, msg, self.aelinks.keys()))

    def AEParam(self, aid, msg=''):
        try:
            return self.aeparams[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, self.aeparams.keys()))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid, msg=''):
        try:
            return self.flfacts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, self.flfacts.keys()))

    def Flutter(self, fid, msg=''):
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, self.flutters.keys()))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid, msg=''):
        try:
            return self.dconstrs[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, self.dconstrs.keys()))

    def Desvar(self, oid, msg=''):
        try:
            return self.desvars[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DESVARs=%s'
                           % (oid, msg, self.desvars.keys()))

    def DDVal(self, oid, msg=''):
        try:
            return self.ddvals[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, self.ddvals.keys()))

    #--------------------
    # SET CARDS

    def Set(self, sid, msg=''):
        try:
            return self.sets[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, self.sets.keys()))

    def SetSuper(self, seid, msg=''):
        try:
            return self.setsSuper[seid]
        except KeyError:
            raise KeyError('seid=%s not found%s.  Allowed SETx=%s'
                           % (seid, msg, self.setsSuper.keys()))

    #--------------------
    # METHOD CARDS
    def Method(self, sid, msg=''):
        try:
            return self.methods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, self.methods.keys()))

    def CMethod(self, sid, msg=''):
        try:
            return self.cmethods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, self.cmethods.keys()))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg=''):
        try:
            return self.tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, self.tables.keys()))

    def RandomTable(self, tid, msg=''):
        try:
            return self.randomTables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, self.randomTables.keys()))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid, msg=''):
        try:
            return self.nlparms[nid]
        except KeyError:
            raise KeyError('nid=%s not found%s.  Allowed NLPARMs=%s'
                           % (nid, msg, self.nlparms.keys()))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname, msg=''):
        try:
            return self.dmig[dname]
        except KeyError:
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, self.dmig.keys()))

    def get_x_associated_with_y(self, xdict, xkeys, ykeys, stop_on_failure=True):
        """
        Get the range of sub-properties of a card.

        .. note::
           Assumes you're taking a single path through the cards.
           You could probably explicitly code these queries faster, but
           this method has a lot of flexibility with very little user code.

        :param self:
            the BDF object
        :param xdict:
            the BDF attribute that should be querried
            (e.g. self.elements)
        :param xkeys:
            the list of object keys that should be stepped through
            associated with xdict
            (e.g. eids=[1, 2, 3])
        :param ykeys:
            the list of response keys that should be stepped through
            (e.g. ['pid', 'mid', 'mid'])
        :param stop_on_failure:
            Should an error be raised if there is an invalid key?
            For example, get all material used by elements, but don't crash
            on CONRODs.
        :returns results:
            The set of all values used

        .. code-block:: python

          # Get nodes associated with eid=[1, 2, 3]
          nodes = self.get_x_associated_with_y(
              self.elements, [1, 2, 3], ['nodes'])

          # Get node IDs associated with eid=[1, 2, 3]
          nodesIDs = self.get_x_associated_with_y(
              self.elements, [1, 2, 3], ['nodes', 'nid'])

          # Get coord IDs associated with eid=[1, 2, 3]
          coordIDs = self.get_x_associated_with_y(
              self.elements, [1, 2, 3], ['nodes', 'cp', 'cid'])

          # Get properties associated with eid=[1, 2, 3]
          properties = self.get_x_associated_with_y(
              self.elements, [1, 2, 3], ['pid'])

          # Get materials associated with eid=[1, 2, 3]
          materials = self.get_x_associated_with_y(
              self.elements, [1, 2, 3], ['pid', 'mid'])

          # Get material IDs associated with eid=[1, 2, 3]
          materialIDs = self.get_x_associated_with_y(
              self.elements, [1, 2, 3], ['pid', 'mid', 'mid'])

          # Get the values for Young's Modulus
          E = self.get_x_associated_with_y(
              self.elements, None, ['pid', 'mid', 'e'])
        """
        assert isinstance(xdict, dict), type(xdict)
        assert self._xref == True, self._xref

        if isinstance(xkeys, list) or isinstance(xkeys, tuple):
            pass
        elif isinstance(xkeys, ndarray):
            assert len(xkeys.shape) == 1, xkeys.shape
            assert isinstance(xkeys[0], int), type(xkeys[0])
        elif xkeys is None:
            xkeys = xdict.iterkeys()
        else:
            raise RuntimeError('invalid type; type(xkeys)=%r' % type(xkeys))

        assert isinstance(ykeys[0], string_types), ykeys

        out_set = set([])
        for xkey in xkeys:
            xi = xdict[xkey]
            _getattr(out_set, xi, ykeys, len(ykeys)-1, stop_on_failure)
        return out_set

    def __test_method(self):
        # getElementsAssociatedWithMaterialIDs(self):
        #getElementIDsAssociatedWithPropertyIDs
        #getPropertyIDsAssociatedWithMaterialIDs
        #get_x_associated_with_y(self, mids, 'pid', yname='mid')

        #CTETRA   1       1       8       13      67      33
        #CTETRA   2       1       8       7       62      59
        #CTETRA   3       1       8       45      58      66
        #nids 7 8 13 33 45 58 59 62 66 67

        cps = self.get_x_associated_with_y(
            self.elements, [1, 2, 3], ['nodes', 'cp', 'cast'], stop_on_failure=False)
        #print('*cps', cps)

        nids2 = self.get_x_associated_with_y(self.elements, [1, 2, 3], ['nodes', 'nid'])
        nids2 = list(nids2)
        nids2.sort()
        #print('*nids2', nids2)

        mids = self.get_x_associated_with_y(
            self.elements, [1, 2, 3, 4, 5], ['pid', 'mid'])
        #print('*mids', mids)

        mids2 = self.get_x_associated_with_y(
            self.elements, None, ['pid', 'mid', 'mid'])

        E = self.get_x_associated_with_y(
            self.elements, None, ['pid', 'mid', 'e'])
        #print('*E', E)

        #get_nodes_associated_with_elements
        #nids = get_x_associated_with_y(self, self.elements, 'nodes', 'nid')
        #pids = get_x_associated_with_y(self, self.elements, 'pid', 'pid')
        #mids = get_x_associated_with_y(self, self.properties, 'mid', 'mid')


def _getattr(out_set, xi, xkeys, nlevels_left=0, stop_on_failure=True):
    """
    Recursive method to help get_x_associated_with_y get the value of xkeys
    for the given variable xi

    :param out_set:
        the SET of all outputs that will be filled and implicitly returned
    :type out_set:
        SET

    :param xi:
        the current variable being iterated over
    :type xi:
        BDF BaseCard

    :param xkeys:
        the variables left to iterate over
    :type xkeys:
        LIST of STRINGS

    :param nlevels_left:
        the number of levels still to iterate over
    :type nlevels_left:
        INT >= 0

    :param stop_on_failure:
        should the code crash if a certain xkey cannot be found
    :type stop_on_failure:
        bool
    """
    try:
        y = getattr(xi, xkeys[0])
    except AttributeError:
        if stop_on_failure:
            raise
        return

    if nlevels_left:
        if isinstance(y, list):
            # handles nodes
            for yi in y:
                _getattr(out_set, yi, xkeys[1:], nlevels_left-1, stop_on_failure=stop_on_failure)
        else:
            _getattr(out_set, y, xkeys[1:], nlevels_left-1, stop_on_failure=stop_on_failure)
    else:
        # handles nodes
        if isinstance(y, list):
            for yi in y:
                out_set.add(yi)
        else:
            out_set.add(y)
