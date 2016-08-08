# pylint: disable=E1101,C0103


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

        .. deprecated:: 0.7
           Use :func:`get_element_ids_list_with_pids([pid], mode='list')`
        """
        self.deprecated('getElementIDsWithPID(pid)', 'get_element_ids_list_with_pids([pid])', '0.8')
        return self.get_element_ids_list_with_pids([pid], mode='list')

    def Flfact(self, sid, msg):
        """
        .. deprecated:: 0.7
            Use :func:`FLFACT(sid, msg)`
        """
        self.deprecated('Flfact(sid)', 'FLFACT(sid)', '0.8')
        return self.FLFACT(sid, msg)

    def getNodes(self):
        """
        .. deprecated:: 0.7
            Use :func:`get_nodes()`
        """
        self.deprecated('getNodes()', 'get_nodes()', '0.8')
        return self.get_nodes()

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

    def get_nplies(self):
        return self.nplies

    def get_material_ids(self):
        return self.material_ids

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
