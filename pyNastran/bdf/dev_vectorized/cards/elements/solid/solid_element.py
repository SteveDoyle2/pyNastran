from numpy import zeros, searchsorted, where

class SolidElement(object):
    def __init__(self, model):
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []
        #self.comments = {}

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

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def get_mass(self, element_ids=None, xyz_cid0=None, total=False):
        """
        Gets the mass for one or more SolidElement elements.

        :param element_ids: the elements to consider (default=None -> all)
        :param xyz_cid0: the positions of the GRIDs in CID=0 (default=None)
        :param total: should the centroid be summed (default=False)
        """
        if element_ids is None:
            element_ids = self.element_id
        V = self.get_volume(element_ids, xyz_cid0)
        mid = self.model.properties_solid.get_mid(self.property_id)
        rho = self.model.materials.get_rho(mid)

        mass = V * rho
        if total:
            mass = mass.sum()
        return mass

    def get_mass_centroid_inertia(self, p=None, element_ids=None, xyz_cid0=None, total=False):
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

    def get_density(self, element_ids=None):
        if element_ids is None:
            element_ids = self.element_id

        rho = []
        i = where(element_ids == self.element_id)[0]
        for pid in self.property_id[i]:
            rhoi = self.model.properties_solid.psolid.get_density(pid)
            rho += rhoi
        return rho
