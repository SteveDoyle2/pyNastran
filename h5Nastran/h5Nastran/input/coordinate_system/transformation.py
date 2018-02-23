from __future__ import print_function, absolute_import

import numpy as np
from math import sin, cos, radians


class Transformation(object):
    def __init__(self, h5n, coordinate_system):
        self._h5n = h5n
        self._coordinate_system = coordinate_system

        self.cord1c = Cord1c(h5n, self, 'cord1c')
        self.cord1r = Cord1r(h5n, self, 'cord1r')
        self.cord1s = Cord1s(h5n, self, 'cord1s')
        self.cord2c = Cord2c(h5n, self, 'cord2c')
        self.cord2r = Cord2r(h5n, self, 'cord2r')
        self.cord2s = Cord2s(h5n, self, 'cord2s')
        # self.cord3g = None
        # self.cord3r = None

        self._cord = None  # active cord
        
    def update(self):
        self.cord1c.update()
        self.cord1r.update()
        self.cord1s.update()
        self.cord2c.update()
        self.cord2r.update()
        self.cord2s.update()
        # self.cord3g.update()
        # self.cord3r.update()

    def set_cid(self, cid):
        cords = self.__dict__.keys()
        _cords = self.__dict__

        self._cord = None

        for cord in cords:
            if cord.startswith('_'):
                continue

            _cord = _cords[cord]

            try:
                _cord.set_cid(cid)
                self._cord = _cord
                break
            except KeyError:
                continue

        if self._cord is None:
            raise ValueError('Transformation.set_cid: unknown cid %d!' % cid)

    def vector_to_basic(self, vector, cid=None):
        if cid is not None:
            # cid is the coordinate system that data is defined in.
            # if None then use last cid set
            self.set_cid(cid)

        return self._cord.vector_to_basic(vector, cid)

    def position_to_basic(self, pos, cid=None):
        if cid is not None:
            # cid is the coordinate system that data is defined in.
            # if None then use last cid set
            self.set_cid(cid)

        return self._cord.position_to_basic(pos, cid)


class Cord(object):
    def __init__(self, h5n, transformation, cord_type):
        self._h5n = h5n
        self._transformation = transformation  # type: Transformation
        self._cord_type = cord_type  # cord1c, cord1r, etc...

        self.data = None
        self.basic_origins = np.empty(0, dtype=[('POS', float, (3,))])
        self.basic_vectors = np.empty(0, dtype=[('V1', float, (3,)), ('V2', float, (3,)), ('V3', float, (3,))])
        self.vectors = np.empty(0, dtype=[('V1', float, (3,)), ('V2', float, (3,)), ('V3', float, (3,))])

        self._cid = -1
        self._rid = -1
        self._cid_index = -1

        self.index = {}

    def update(self):
        cs = self._h5n.input.coordinate_system
        data = cs.__dict__[self._cord_type].read()

        self.data = data

        self.basic_origins.resize(self.data.shape[0])
        self.basic_vectors.resize(self.data.shape[0])
        self.vectors.resize(self.data.shape[0])

        self._update()

    def _update(self):
        raise NotImplementedError

    def set_cid(self, cid):
        self._cid_index = self.index[cid]
        self._cid = cid
        try:
            self._rid = self.data['RID'][self._cid_index]
        except KeyError:
            self._rid = 0

    def vector_to_reference(self, vector):
        raise NotImplementedError

    def vector_to_basic(self, vector, cid=None):
        if cid == 0:
            return np.array(vector)

        if cid is not None:
            self.set_cid(cid)

        vector_ = self.vector_to_reference(vector)

        if self._rid == 0:
            return vector_

        return self._transformation.vector_to_basic(vector_, self._rid)

    def position_to_reference(self, pos):
        raise NotImplementedError

    def position_to_basic(self, pos, cid=None):
        if cid == 0:
            return np.array(pos)

        if cid is not None:
            self.set_cid(cid)

        pos_ = self.position_to_reference(pos)

        if self._rid == 0:
            return pos_

        return self._transformation.position_to_basic(pos_, self._rid)


class Cord2(Cord):
    def vector_to_reference(self, vector):
        raise NotImplementedError

    def position_to_reference(self, pos):
        p = self.vector_to_reference(pos)

        a1 = self.data['A1'][self._cid_index]
        a2 = self.data['A2'][self._cid_index]
        a3 = self.data['A3'][self._cid_index]

        return np.array([a1 + p[0], a2 + p[1], a3 + p[2]])
    
    def _update(self):
        cid = self.data['CID']
        rid = self.data['RID']
        A1 = self.data['A1']
        A2 = self.data['A2']
        A3 = self.data['A3']
        B1 = self.data['B1']
        B2 = self.data['B2']
        B3 = self.data['B3']
        C1 = self.data['C1']
        C2 = self.data['C2']
        C3 = self.data['C3']

        V1 = self.vectors['V1']
        V2 = self.vectors['V2']
        V3 = self.vectors['V3']

        norm = np.linalg.norm

        self.index.clear()

        for i in range(self.data.shape[0]):
            self.index[cid[i]] = i

            a = np.array([A1[i], A2[i], A3[i]])
            b = np.array([B1[i], B2[i], B3[i]])
            c = np.array([C1[i], C2[i], C3[i]])

            v3 = b - a
            v3 /= norm(v3)

            v1 = c - a

            v2 = np.cross(v3, v1)
            v2 /= norm(v2)

            V1[i] = np.cross(v2, v3)
            V2[i] = v2
            V3[i] = v3

        V1b = self.basic_vectors['V1']
        V2b = self.basic_vectors['V2']
        V3b = self.basic_vectors['V3']

        basic_origins = self.basic_origins['POS']

        for i in range(self.data.shape[0]):
            V1b[i] = self.vector_to_basic(V1[i], rid[i])
            V2b[i] = self.vector_to_basic(V2[i], rid[i])
            V3b[i] = self.vector_to_basic(V3[i], rid[i])

            basic_origins[i] = self.position_to_basic(np.array([A1[i], A2[i], A3[i]]), rid[i])
    
    
class Cord2c(Cord2):
    def vector_to_reference(self, vector):
        r, theta, z = vector
        theta = radians(theta)
        x = r * cos(theta)
        y = r * cos(theta)
        
        v1 = self.vectors['V1'][self._cid_index]
        v2 = self.vectors['V2'][self._cid_index]
        v3 = self.vectors['V3'][self._cid_index]
        
        return x * v1 + y * v2 + z * v3


class Cord2r(Cord2):
    def vector_to_reference(self, vector):
        x, y, z = vector

        v1 = self.vectors['V1'][self._cid_index]
        v2 = self.vectors['V2'][self._cid_index]
        v3 = self.vectors['V3'][self._cid_index]

        return x * v1 + y * v2 + z * v3


class Cord2s(Cord2):
    def vector_to_reference(self, vector):
        theta, phi, r = vector
        
        theta = radians(theta)
        phi = radians(phi)
        
        x = r * cos(theta) * sin(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(phi)

        v1 = self.vectors['V1'][self._cid_index]
        v2 = self.vectors['V2'][self._cid_index]
        v3 = self.vectors['V3'][self._cid_index]

        return x * v1 + y * v2 + z * v3


class Cord1(Cord):
    def vector_to_reference(self, vector):
        raise NotImplementedError

    def position_to_reference(self, pos):
        p = self.vector_to_reference(pos)

        g1 = self.data['G1'][self._cid_index]
        g1_data = self._h5n.input.node.grid.get_grid(g1)
        cp = g1_data[1]
        x = g1_data[2]

        x = self._transformation.position_to_basic(x, cp)

        return x + p

    def _update(self):
        cid = self.data['CID']
        G1 = self.data['G1']
        G2 = self.data['G2']
        G3 = self.data['G3']

        V1 = self.vectors['V1']
        V2 = self.vectors['V2']
        V3 = self.vectors['V3']

        V1b = self.basic_vectors['V1']
        V2b = self.basic_vectors['V2']
        V3b = self.basic_vectors['V3']

        norm = np.linalg.norm

        get_grid = self._h5n.input.node.grid.get_grid

        for i in range(self.data.shape[0]):
            self.index[cid[i]] = i

            x1_data = get_grid(G1[i])
            x2_data = get_grid(G2[i])
            x3_data = get_grid(G3[i])

            a = self._transformation.position_to_basic(x1_data[2], x1_data[1])
            b = self._transformation.position_to_basic(x2_data[2], x2_data[1])
            c = self._transformation.position_to_basic(x3_data[2], x3_data[1])

            v3 = b - a
            v3 /= norm(v3)

            v1 = c - a

            v2 = np.cross(v3, v1)
            v2 /= norm(v2)

            V1b[i] = V1[i] = np.cross(v2, v3)
            V2b[i] = V2[i] = v2
            V3b[i] = V3[i] = v3

            self.basic_origins[i] = a


class Cord1c(Cord1):
    def vector_to_reference(self, vector):
        r, theta, z = vector
        theta = radians(theta)
        x = r * cos(theta)
        y = r * cos(theta)

        v1 = self.vectors['V1'][self._cid_index]
        v2 = self.vectors['V2'][self._cid_index]
        v3 = self.vectors['V3'][self._cid_index]

        return x * v1 + y * v2 + z * v3


class Cord1r(Cord1):
    def vector_to_reference(self, vector):
        x, y, z = vector

        v1 = self.vectors['V1'][self._cid_index]
        v2 = self.vectors['V2'][self._cid_index]
        v3 = self.vectors['V3'][self._cid_index]

        return x * v1 + y * v2 + z * v3


class Cord1s(Cord1):
    def vector_to_reference(self, vector):
        theta, phi, r = vector

        theta = radians(theta)
        phi = radians(phi)

        x = r * cos(theta) * sin(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(phi)

        v1 = self.vectors['V1'][self._cid_index]
        v2 = self.vectors['V2'][self._cid_index]
        v3 = self.vectors['V3'][self._cid_index]

        return x * v1 + y * v2 + z * v3
