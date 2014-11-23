from six.moves import StringIO
from numpy import zeros, searchsorted, where, asarray, array

from pyNastran.bdf.dev_vectorized.cards.elements.element import Element

class SolidElement(Element):
    def __init__(self, model):
        Element.__init__(self, model)

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

    def __repr__(self):
        f = StringIO()
        f.write('<%s object> n=%s\n' % (self.type, self.n))
        self.write_bdf(f)
        return f.getvalue()

    def _get_node_locations_by_element_id(self, element_id=None, xyz_cid0=None):
        i = self._get_sorted_index(self.element_id, element_id, self.n,
                                   'element_id', 'element_id in %s' % self.type, check=True)
        self.model.log.debug('ielem = %s' % i)
        if xyz_cid0 is None:
            xyz_cid0 = self.model.grid.get_position_by_index()
        return self._get_node_locations_by_index(i, xyz_cid0)

    def allocate(self, ncards):
        self.n = ncards
        float_fmt = self.model.float
        self.element_id = zeros(ncards, 'int32')
        self.property_id = zeros(ncards, 'int32')
        self.node_ids = zeros((ncards, self.nnodes), 'int32')
        #self._comments.append(comment)

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def get_mass_by_element_id(self, element_id=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more SolidElement elements.

        :param element_id: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_id is None:
            element_id = self.element_id
        V = self.get_volume_by_eLement_id(element_id, xyz_cid0)

        mid = self.model.properties_solid.get_material_id_by_property_id(self.property_id)
        rho = self.model.materials.get_density_by_material_id(mid)

        rho = self.model.properties_solid.psolid.get_density_by_property_id(self.property_id)
        #rho = self.model.materials.get_density_by_material_id(mid)

        try:
            mass = V * rho
        except ValueError:
            msg = 'element_id = %s; n=%s\n' % (element_id, len(element_id))
            msg += 'mid=%s\n' % mid
            msg += 'rho=%s\n' % rho
            msg += 'V.shape = %s\n' % str(V.shape)
            msg += 'rho.shape = %s' % str(rho.shape)
            print(msg)
            raise

        n = len(element_id)
        assert mass.shape == (n, ), mass.shape
        if total:
            mass = mass.sum()
        return mass

    def get_mass_centroid_inertia_by_eLement_id(self, p=None, element_id=None, xyz_cid0=None, total=False):
        """
        Calculates the mass, centroid, and (3, 3) moment of interia
        matrix.  Considers position, but not the (hopefully) small
        elemental term.

        :param p: the point to take the moment of inertia about (default=None -> origin)

        a  = integral(mu * (y^2 + z^2), dV)
        b  = integral(mu * (x^2 + z^2), dV)
        c  = integral(mu * (y^2 + y^2), dV)
        a' = integral(mu * (yz), dV)
        b' = integral(mu * (xz), dV)
        c' = integral(mu * (xy), dV)

        I = [ a  -b', -c']
            [-b'  b   -a']
            [-c' -a'   c ]

        Exact MOI for tetrahedron
        http://www.thescipub.com/abstract/?doi=jmssp.2005.8.11
        """
        if p is None:
            p = zeros(3, self.model.float)

        r = centroid - p  # 2D array - 1D array
        I = mass * r**2 # column vector * 2D array
        return mass, centroid, I

    def get_density_by_element_id(self, element_id=None):
        if element_id is None:
            element_id = self.element_id

        n = len(element_id)
        rho = zeros(n, dtype='float64')
        i = where(element_id == self.element_id)[0]
        for pid in self.property_id[i]:
            j = where(pid == self.property_id[i])[0]
            rhoi = self.model.properties_solid.psolid.get_density_by_property_id(pid)
            rho[j] = rhoi
        assert rho.shape == (n, ), rho.shape
        return rho

    def slice_by_index(self, i):
        i = asarray(i)
        #name = self.__class__.__name__
        #obj_class = type(name, (SolidElement, ), {})
        obj_class = self.__class__#.__class__
        obj = obj_class(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        return obj
