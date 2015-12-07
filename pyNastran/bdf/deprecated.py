# pylint: disable=E1101,C0103
import warnings
from numpy import array


class GetMethodsDeprecated(object):
    """defines deprecated methods for GetMethods"""
    def getElementIDsWithPID(self, pid):
        """
        Gets all the element IDs with a specific property ID

        Parameters
        ----------
        pid : int
             property ID

        Returns
        -------
        elementIDs : List[int]
            as a list

        .. deprecated:: will be removed in version 0.8

        The same functionality may be used by calling
          >>> self.getElementIDsWithPIDs([pid], mode='list')
        """
        self.deprecated('getElementIDsWithPID(pid)', 'get_element_ids_list_with_pids([pid])', '0.8')
        return self.getElementIDsWithPIDs([pid], mode='list')

    def Flfact(self, sid, msg):
        """
        .. deprecated:: will be removed in version 0.8
        """
        self.deprecated('Flfact(sid)', 'FLFACT(sid)', '0.8')
        return self.FLFACT(sid, msg)

    def getNodes(self):
        """deprecated"""
        self.deprecated('getNodes()', 'get_nodes()', '0.8')
        return self.get_nodes()

    #def getNodeIDsWithElement(self, eid):
        #"""deprecated"""
        #self.deprecated('getNodeIDsWithElement(eid)', 'get_node_ids_with_element(eid)', '0.8')
        #return self.get_node_ids_with_element(eid)

    #def getNodeIDsWithElements(self, eids, msg=''):
        #"""deprecated"""
        #self.deprecated('getNodeIDsWithElements(eids)', 'get_node_ids_with_elements(eids)', '0.7')
        #return self.get_node_ids_with_elements(eids, msg=msg)

    #------------------------------
    #def nNodes(self):
        #"""deprecated"""
        #self.deprecated('nNodes()', 'nnodes', '0.7')
        #return self.nnodes

    #def get_nnodes(self):
        #"""deprecated"""
        #self.deprecated('get_nnodes()', 'nnodes', '0.7')
        #return self.nnodes

    #def nodeIDs(self):
        #"""deprecated"""
        #self.deprecated('nodeIDs()', 'node_ids', '0.7')
        #return self.node_ids

    #def get_node_ids(self):
        #"""deprecated"""
        #self.deprecated('get_node_ids()', 'node_ids', '0.7')
        #return self.node_ids

    #------------------------------
    #def nElements(self):
        #"""deprecated"""
        #self.deprecated('nElements()', 'nelements', '0.7')
        #return self.nelements

    #def elementIDs(self):
        #"""deprecated"""
        #self.deprecated('elementIDs()', 'element_ids', '0.7')
        #return self.element_ids

    #def get_nelements(self):
        #"""deprecated"""
        #self.deprecated('get_nelements()', 'nelements', '0.7')
        #return self.element_ids

    #def get_element_ids(self):
        #"""deprecated"""
        #self.deprecated('get_element_ids()', 'element_ids', '0.7')
        #return self.element_ids

    #------------------------------
    #def nProperties(self):
        #"""deprecated"""
        #return self.nproperties

    #def get_nproperties(self):
        #"""deprecated"""
        #return len(self.properties)

    #def get_property_ids(self):
        #"""deprecated"""
        #self.deprecated('get_property_ids()', 'property_ids', '0.7')
        #return self.property_ids

    #def propertyIDs(self):
        #"""deprecated"""
        #self.deprecated('propertyIDs()', 'property_ids', '0.7')
        #return self.property_ids

    ##------------------------------

    #def get_coord_ids(self):
        #self.deprecated('get_coord_ids()', 'coord_ids', '0.7')
        #return self.coord_ids

    #def coordIDs(self):
        #"""deprecated"""
        #self.deprecated('coordIDs()', 'coord_ids', '0.7')
        #return self.coord_ids

    #------------------------------

    #def nCAeros(self):
        #"""deprecated"""
        #self.deprecated('nCAeros()', 'ncaeros', '0.7')
        #return self.get_ncaeros()

    #def get_ncaeros(self):
        #self.deprecated('get_ncaeros()', 'ncaeros', '0.7')
        #return self.ncaeros

    #==================================================

    #def getElementIDsWithPIDs(self, pids, mode='list'):
        #"""deprecated"""
        #self.deprecated('getElementIDsWithPIDs()', 'get_element_ids_with_pids(pids, mode)', '0.7')
        #return self.get_element_ids_with_pids(pids, mode=mode)


    def _get_element_ids_with_pids(self, pids, mode='list'):
        if mode not in ['list', 'dict']:
            msg = "mode=%r is not supported.  Use 'list' or 'dict'\n" % mode
            raise ValueError(msg)
        if mode == 'list':
            self.deprecated('get_element_ids_with_pids(pids, "list")', 'get_element_ids_list_with_pids(pids)', '0.8')
            return self.get_element_ids_list_with_pids(pids)
        else:
            self.deprecated('get_element_ids_with_pids(pids, "dict")', 'get_element_ids_dict_with_pids(pids)', '0.8')
            return self.get_element_ids_dict_with_pids(pids)

    #def getNodeIDToElementIDsMap(self):
        #"""deprecated"""
        #self.deprecated('getNodeIDToElementIDsMap()', 'get_element_ids_with_pids()', '0.7')
        #return self.get_node_id_to_element_ids_map()

    #def getPropertyIDToElementIDsMap(self):
        #"""deprecated"""
        #self.deprecated('getPropertyIDToElementIDsMap()', 'get_property_id_to_element_ids_map()', '0.7')
        #return self.get_property_id_to_element_ids_map()

    #def getMaterialIDToPropertyIDsMap(self):
        #"""deprecated"""
        #self.deprecated('getMaterialIDToPropertyIDsMap()', 'get_material_id_to_property_ids_map()', '0.7')
        #return self.get_material_id_to_property_ids_map()

    #def structuralMaterialIDs(self):
        #"""deprecated"""
        #self.deprecated('structuralMaterialIDs()', 'get_structural_material_ids()', '0.7')
        #return self.get_structural_material_ids()

    #def thermalMaterialIDs(self):
        #"""deprecated"""
        #self.deprecated('thermalMaterialIDs()', 'get_thermal_material_ids()', '0.7')
        #return self.get_thermal_material_ids()

    #def materialIDs(self):
        #"""deprecated"""
        #self.deprecated('materialIDs()', 'get_material_ids()', '0.7')
        #return self.get_material_ids()


class DeprecatedCompositeShellProperty(object):
    """
    To be deprecated in:
      - Version 0.7
    To be removed in:
      - Version 0.8
    """
    def MassPerArea(self, iply='all', method='nplies'):
        #self.deprecated('MassPerArea(iply, method)', 'get_mass_per_area(iply, method)', '0.8')
        return self.get_mass_per_area(iply, method)

    def Thickness(self, iply='all'):
        return self.get_thickness(iply)

    def nPlies(self):
        return self.get_nplies()

    def Nsm(self):
        return self.get_nonstructural_mass()

    def isSymmetrical(self):
        return self.is_symmetrical()

    def Rho(self, iply):
        return self.get_density(iply)

    def Theta(self, iply):
        return self.get_theta(iply)

    def sout(self, iply):
        return self.get_sout(iply)
