from numpy import zeros, searchsorted, where, asarray

from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class ShellElement(Element):
    def __init__(self, model):
        Element.__init__(self, model)

    def get_index_by_element_id(self, element_id=None, msg=''):
        i = self._get_sorted_index(self.element_id, element_id, self.n, 'element_id in %s%s' % (self.type, msg), check=True)
        return i

    def get_property_index_from_property_id(self, property_id=None, i=None):
        """Find all the j-indicies where seid=seidi for some given subset of i-indicies"""
        return self._get_index_from_param('property_id', self.property_id, property_id, i)

    def __getitem__(self, element_ids):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.element_id, element_ids)
        return self.slice_by_index(i)

    def get_mass(self, element_id=None, total=False, node_ids=None, xyz_cid0=None):
        return self.get_mass_by_element_id(element_id, total, node_ids, xyz_cid0)

    def get_mass_by_element_id(self, element_id=None, total=False, node_ids=None, xyz_cid0=None):
        """
        Gets the mass of the CQUAD4s on a total or per element basis.

        :param self: the CQUAD4 object
        :param element_id: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        mass, _area, _normal = self._mass_area_normal(element_id=element_id,
            xyz_cid0=xyz_cid0,
            calculate_mass=True, calculate_area=False,
            calculate_normal=False)

        if total:
            return mass.sum()
        else:
            #print('mass.shape = %s' % mass.shape)
            return mass

    def get_normal(self, element_id=None, xyz_cid0=None):
        return self.get_normal_by_element_id(element_id, xyz_cid0)

    def get_normal(self, element_id=None, xyz_cid0=None):
        """
        Gets the normals of the CQUAD4s on per element basis.

        :param self: the CQUAD4 object
        :param element_id: the elements to consider (default=None -> all)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, normal = self._mass_area_normal(element_id=element_id,
            xyz_cid0=xyz_cid0,
            calculate_mass=False, calculate_area=False,
            calculate_normal=True)
        return normal


    def get_area(self, element_id=None, total=False, xyz_cid0=None):
        return self.get_area_by_element_id(element_id, total, xyz_cid0)

    def get_area_by_element_id(self, element_id=None, total=False, xyz_cid0=None):
        """
        Gets the area of the CQUAD4s on a total or per element basis.

        :param self: the CQUAD4 object
        :param element_id: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        ..note:: If node_ids is None, the positions of all the GRID cards
                 must be calculated
        """
        _mass, area, _normal = self._mass_area_normal(element_id=element_id,
            xyz_cid0=xyz_cid0,
            calculate_mass=False, calculate_area=True,
            calculate_normal=False)
        if total:
            return area.sum()
        else:
            return area

    def get_thickness(self, element_id=None):
        return self.get_thickness_by_element_id(element_id)

    def get_thickness_by_element_id(self, element_id=None):
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]
        #print 'element_id =', element_id
        #print 'property_id =', property_id
        t = self.model.properties_shell.get_thickness(property_id)
        return t

    def get_nonstructural_mass(self, element_id=None):
        return self.get_nonstructural_mass_by_element_id(element_id)

    def get_nonstructural_mass_by_element_id(self, element_id=None):
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]
        nsm = self.model.properties_shell.get_nonstructural_mass(property_id)
        return nsm

    def get_density(self, element_id=None):
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]
        #print('density - element_id = %s' % element_id)
        density = self.model.properties_shell.get_density(property_id)
        #print('density_out = %s' % density)
        return density

    def slice_by_index(self, i):
        i = asarray(i)
        #name = self.__class__.__name__
        #print('name = %r' % name)
        #obj = CQUAD4(self.model)
        #obj_class = type(name, (ShellElement, ), {})
        obj_class = self.__class__#.__class__
        obj = obj_class(self.model)
        #print(type(obj))

        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        obj.zoffset = self.zoffset[i]
        obj.t_flag = self.t_flag[i]
        obj.thickness = self.thickness[i, :]
        return obj

