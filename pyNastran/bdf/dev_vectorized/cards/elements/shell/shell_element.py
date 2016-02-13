from numpy import searchsorted, asarray, isnan

from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class ShellElement(Element):
    def __init__(self, model):
        Element.__init__(self, model)

    def get_element_id_by_element_index(self, i=None):
        if i is None:
            element_id = self.element_id
        else:
            element_id = self.element_id[i]
        return element_id

    def get_element_index_by_element_id(self, element_id=None, msg=''):
        i = self._get_sorted_index(self.element_id, element_id, 'element_id', 'element_id in %s%s' % (self.type, msg), check=True)
        return i

    def get_property_id_by_element_id(self, element_id=None):
        i = self.get_element_index_by_element_id(element_id)
        return self.get_property_id_by_element_index(i)

    def get_property_id_by_element_index(self, i=None):
        if i is None:
            property_id = self.property_id
        else:
            property_id = self.property_id[i]
        return property_id

    def get_property_index_by_property_id(self, property_id=None, i=None):
        """Find all the j-indicies where seid=seidi for some given subset of i-indicies"""
        return self._get_index_by_param('property_id', self.property_id, property_id, i)

    def __getitem__(self, i):
        return self.slice_by_index(i)

    def slice_by_element_id(self, element_id):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.element_id, element_id)
        return self.slice_by_index(i)

    def get_mass_by_element_id(self, element_id=None, total=False, node_ids=None, xyz_cid0=None):
        """
        Gets the mass of the CQUAD4s on a total or per element basis.

        :param element_id: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
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

    def get_mass_per_area_by_element_id(self, element_id=None, node_ids=None, xyz_cid0=None):
        """
        Gets the mass per area of the CQUAD4s on a total or per element basis.

        :param element_id: the elements to consider (default=None -> all)
        :param total: should the mass be summed (default=False)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        mass, _area, _normal = self._mass_area_normal(element_id=element_id,
            xyz_cid0=xyz_cid0,
            calculate_mass=True, calculate_area=True,
            calculate_normal=False)
        mpa = mass / _area
        return mpa

    def get_normal_by_element_id(self, element_id=None, xyz_cid0=None):
        """
        Gets the normals of the CQUAD4s on per element basis.

        :param element_id: the elements to consider (default=None -> all)

        :param xyz_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
                  must be calculated
        """
        _mass, area, normal = self._mass_area_normal(element_id=element_id,
                                                     xyz_cid0=xyz_cid0,
                                                     calculate_mass=False, calculate_area=False,
                                                     calculate_normal=True)
        return normal

    def get_area_by_element_id(self, element_id=None, total=False, xyz_cid0=None):
        """
        Gets the area of the CQUAD4s on a total or per element basis.

        :param element_id: the elements to consider (default=None -> all)
        :param total: should the area be summed (default=False)

        :param node_ids:   the GRIDs as an (N, )  NDARRAY (or None)
        :param grids_cid0: the GRIDs as an (N, 3) NDARRAY in CORD2R=0 (or None)

        .. note:: If node_ids is None, the positions of all the GRID cards
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

    def get_thickness_by_element_id(self, element_id=None):
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]
        #print 'element_ids =', element_ids
        #print 'property_ids =', property_id
        t = self.model.properties_shell.get_thickness_by_property_id(property_id)
        return t

    def get_nonstructural_mass_by_element_id(self, element_id=None):
        if element_id is None:
            element_ids = self.element_id
            property_id = self.property_id
            i = None
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]
        nsm = self.model.properties_shell.get_nonstructural_mass_by_property_id(property_id)
        return nsm

    def get_density_by_element_id(self, element_id=None):
        if element_id is None:
            element_id = self.element_id
            property_id = self.property_id
            i = None
            n = self.n
        else:
            i = searchsorted(self.element_id, element_id)
            property_id = self.property_id[i]
            n = len(element_id)

        #print('density - element_ids = %s' % element_id)
        density = self.model.properties_shell.get_density_by_property_id(property_id)
        #self.model.log.debug('density_out = %s' % density)
        assert density.shape == (n, ), density.shape
        assert isnan(density) == False, density
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

