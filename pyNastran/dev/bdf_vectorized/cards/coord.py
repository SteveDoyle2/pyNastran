from numpy import (array, searchsorted, zeros, full,
                   nan, where, vstack, dot, cross, degrees, radians, arctan2,
                   cos, sin, arccos, hstack, eye, ndarray, sqrt, unique,
                   transpose, asarray, isnan, array_equal)
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.cards.coordinate_systems import (
    CORD1R, CORD1C, CORD1S,
    CORD2R, CORD2C, CORD2S)

#from pyNastran.bdf.field_writer_8 import print_card_8
#from pyNastran.bdf.field_writer_16 import print_card_16
#from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard


def normalize(vector):
    #print('-----------')
    #print('normalize(v); v=%s ' % v)
    ni = norm(vector, axis=0)
    #print('isnan', isnan(ni))
    if ni.min() == 0.0:
        msg = 'vector=%s\n' % vector
        msg = 'ni=%s\n' % ni
        raise RuntimeError(msg)
    return vector / ni

# .. todo:: incomplete
class Coord(VectorizedCard):
    def get_global_position_by_node_id(self, node_id, cp):
        i = self.model.grid.get_node_index_by_node_id(node_id)
        self.model.log.info('i = %s; type=%s' % (i, type(i)))
        xyz = self.model.grid.xyz[i, :]
        xyz = xyz.reshape(len(i), 3)
        return self.get_global_position_by_xyz(xyz, cp)

    def get_global_position_by_xyz(self, xyz, cp):
        assert isinstance(cp, int), cp
        assert isinstance(xyz, ndarray), xyz
        coord = self.coords[cp]

        T = coord.beta()
        assert coord.is_resolved, coord
        origin = coord.origin
        assert origin is not None, 'origin is None...\n%s' % str(coord)
        #self.model.log.info('originA = %s' % origin)

        origin = origin.reshape(3, 1)
        #self.model.log.info('origin = %s' % origin)
        #self.model.log.info('T.shape=%s' % str(T.shape))
        #self.model.log.info('T=\n%s' % T)
        #self.model.log.info('xyz=%s' % xyz)
        #self.model.log.info('xyz.shape=%s' % str(xyz.shape))
        xyz2 = (dot(T.T, xyz.T)  + origin).T  # TODO: should this have T.T?
        #print('xyz = %s' % xyz.T)
        #print('xyz2 = %s' % xyz2)
        assert xyz.shape == xyz2.shape, "xyz.shape=%s xyz2.shape=%s" % (xyz.shape, xyz2.shape)
        return xyz2

    def write_card(self, bdf_file, size, is_double, coord_id=None):
        assert size in [8, 16], size
        assert is_double in [True, False], is_double

        if self.n:
            #if coord_id is None:
            #i = arange(self.n)
            #else:
                #assert len(unique(coord_id))==len(coord_id), unique(coord_id)
                #i = searchsorted(self.coord_id, coord_id)

            #if size == 8:
                #print_cardi = print_card_8
            #elif is_double:
                #print_cardi = print_card_double
            #else:
                #print_cardi = print_card_16

            if len(self.coords) > 1:
                bdf_file.write('$COORDs\n')

            if max(self.coords) > self.model.max_int:
                size = 16
            for cid, coord in self.coords.items():
                if cid > 0:
                    #if cid in self._comments:
                        #bdf_file.write(self._comments[cid])
                    #if coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
                        #card = [coord.type, cid, coord.g1, coord.g2, coord.g3]
                    #else:
                        #card = [coord.type, cid, coord.rid] + list(coord.e1) + list(coord.e2) + list(coord.e3)
                    #bdf_file.write(print_cardi(card))
                    bdf_file.write(coord.write_card(size=size))

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        msg = '<Coord object>; Instead call:\n'
        msg += '>>> coords[coord_id]'
        return msg

    def __getitem__(self, i):
        self.log.debug('i**** = %s' % i)
        #i = self.get_coord_index_for_coord_id()
        #return self.coords[value]
        return self.slice_by_index(i)

    def get_coord_index_by_coord_id(self, coord_id=None, msg=''):
        if coord_id == 0:
            return self.coords[0]
        i = self._get_sorted_index(self.coord_id, coord_id, 'coord_id', 'coord_id in %s%s' % (self.type, msg), check=True)
        return i

    def slice_by_coord_id(self, coord_id=None):
        self.log.debug('coord_id = %s' % coord_id)
        #self.log.debug('self.coord_id =', self.coord_id)
        i = self.get_coord_index_by_coord_id(coord_id)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        self.log.debug('i**** = %s' % i)
        assert i.max() < self.n, 'i=%s n=%s i.shape=%s' % (i, self.n, i.shape)
        obj = Coord(self.model)
        obj.n = len(i)
        obj.coord_id = self.coord_id[i]
        for coord_id in obj.coord_id:
            obj.coords[coord_id] = self.coords[coord_id]
        obj.Type = self.Type[i]
        obj.origin = self.origin[i, :]
        obj.is_resolved = self.is_resolved[i]
        obj.T = self.T[i, :, :]
        return obj

    def get_cid_by_coord_id(self, coord_id=None):
        i = self.get_coord_index_by_coord_id(coord_id, msg='')
        return self.get_cid_by_coord_index(i)

    def get_cid_by_coord_index(self, i=None):
        return self.coord_id[i]

    def get_coord_id_by_coord_index(self, i):
        return self.coord_id[i]

    def get_rid_by_coord_id(self, coord_id=None):
        #i = self.get_coord_index_by_coord_id(coord_id, msg='')
        #return self.get_rid_by_coord_index(i)
        if coord_id is None:
            coord_id = self.coord_id
        coord_id = self._set_as_array(coord_id)
        rids = []
        for cid in coord_id:
            rid = self.coords[cid].Rid()
            rids.append(rid)
        return array(rids, dtype='int32')

    def get_rid_by_coord_index(self, i=None):
        coord_id = self.get_coord_id_by_coord_index(i)
        return self.get_rid_by_coord_id(coord_id)

    def allocate(self, ncards=None, card_count=None):
        float_fmt = self.model.float_fmt
        assert ncards is not None or card_count is not None
        if ncards is None:
            ncards = array([card_count[name]
                            for name in ['CORD1R', 'CORD1C', 'CORD1S',
                                         'CORD2R', 'CORD2C', 'CORD2S']
                            if name in card_count], dtype='int32').sum() + 1
        #ncards += 1
        self.model.log.debug('nCOORDcards = %s' % ncards)
        #print('ncards coord = %s' % ncards)
        self.coord_id = zeros(ncards, dtype='int32')
        self.Type = full(ncards, nan, dtype='|U1')  # R-CORD2R, S-CORD2S, C-CORD2C
        self.T = full((ncards, 3, 3), nan, dtype=float_fmt)
        self.origin = full((ncards, 3), nan, dtype=float_fmt)
        self.is_resolved = full(ncards, False, dtype='bool')

        self.Type[0] = 'R'
        self.T[0, :, :] = eye(3)
        self.origin[0, :] = [0., 0., 0.]
        self.is_resolved[0] = True
        self.i = 1

    def __init__(self, model):
        """
        Defines the ShellProperties object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        VectorizedCard.__init__(self, model)
        float_fmt = self.model.float_fmt

        self.n = 1
        ncards = 1
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        self.coords = {0: CORD2R(0, origin, zaxis, xzplane),}
        self.coord_id = zeros(ncards, dtype='int32')
        self.Type = full(ncards, 'R', dtype='|U1')  # R-CORD2R, S-CORD2S, C-CORD2C

        self.T = full((ncards, 3, 3), nan, dtype=float_fmt)
        self.T[0, :, :] = eye(3)
        self.origin = zeros((ncards, 3), dtype=float_fmt)
        self.is_resolved = full(ncards, True, dtype='bool')

        #self.cord2r = CORD2R()
        #self.cord2c = CORD2C()
        #self.cord2s = CORD2S()

        #cid = concatenate(pshell.cid, pcomp.cid)
        #unique_cids = unique(cid)
        #if unique_cids != len(cid):
        #    raise RuntimeError('There are duplicate PSHELL/PCOMP IDs...')

    def build(self, coord_id=None):
        self.model.log.debug('----------building COORDx-------------')
        cids_to_resolve = []
        ncoords = len(self.coords.keys())
        self.allocate(ncards=ncoords)
        #print('coord_ids = %s' % self.coords.keys())
        #print('T = \n%s' % self.T)
        for i, (cid, coord) in enumerate(sorted(self.coords.items())):
            #self.model.log.debug('i=%s cid=%s' % (i, cid))
            self.coord_id[i] = cid
            self.Type[i] = coord.Type
            if coord.is_resolved:
                self.is_resolved[i] = True
                self.origin[i, :] = coord.origin
                self.T[i, :, :] = vstack([
                    coord.i[:],
                    coord.j[:],
                    coord.k[:]])
            else:
                self.model.log.debug('need to resolve cid=%i rid=%i Type=%s' % (cid, coord.rid, coord.Type))
                cids_to_resolve.append(cid)
        self.model.log.debug('coord_id = %s' % self.coord_id)
        #self.model.log.debug('T =\n%s' % self.T)
        self.resolve_coords(cids_to_resolve)

    def update(self, maps):
        """
        maps = {
            'node' : nid_map,
            'coord' : cid_map,
        }
        """
        nid_map = maps['node']
        cid_map = maps['coord']
        coords2 = {}
        for i, cid in enumerate(self.coord_id):
            cid2 = cid_map[cid]
            self.coord_id[i] = cid2
            coord = self.coords[cid]
            coord.update(nid_map, cid_map)
            coords2[cid2] = coord

    def resolve_coords(self, cids_to_resolve):
        while cids_to_resolve:
            cids_to_resolve2 = []
            #print('need to resolve cids=%s' % cids_to_resolve)
            #print('is_resolved=%s' % self.is_resolved)
            #i = where(self.is_resolved == False)[0]
            #cids_to_resolve = self.coord_id[i]
            for cid in cids_to_resolve:
                i = where(self.coord_id == cid)[0]
                coord = self.coords[cid]
                n = coord.type[-2]
                if n == '2':
                    rid = coord.rid
                    #print('need to resolve rid=%s...' % rid)
                    j = where(self.coord_id == rid)[0]
                    is_resolved = self.is_resolved[j]
                    #print('  is rid=%i resolved -> %s' % (rid, is_resolved))
                    if is_resolved:
                        ref_coord = self.coords[rid]
                        ref_coord_type = self.Type[j][0]
                        self.resolve_coord2(i, coord,
                                            j, ref_coord_type, ref_coord)
                elif n == '1':
                    g123 = array([coord.g1, coord.g2, coord.g3])
                    grids = self.model.grid
                    #print('COORD1 grids = %s' % g123)
                    inode = grids.get_node_index_by_node_id(g123)
                    node_cids = grids.cp[inode]
                    xyz = grids.xyz[inode, :]
                    #print('need to resolve cid=%s node_cps=%s...' % (cid, node_cids))
                    j = searchsorted(self.coord_id, node_cids)
                    is_resolved = self.is_resolved[j].all()
                    if is_resolved:
                        #print('  resolved -> True is_resolved=%s' % is_resolved)
                        ref_coord = [
                            self.coords[node_cids[0]],
                            self.coords[node_cids[1]],
                            self.coords[node_cids[2]], ]
                        ref_coord_type = [
                            self.Type[j[0]],
                            self.Type[j[1]],
                            self.Type[j[2]], ]
                        #print('COORD1x j=%s' % j)
                        self.resolve_coord1(i, coord, xyz,
                                            j, ref_coord_type, ref_coord)
                else:
                    #print('rid=%s is not resolved' % rid)
                    cids_to_resolve2.append(cid)
            if len(cids_to_resolve) >= cids_to_resolve2:
                msg = 'circular reference...\ncoord_ids=%s\n' % (cids_to_resolve2)
                raise RuntimeError(msg)
            cids_to_resolve = cids_to_resolve2
            #break
        #print('is_resolved = %s' % self.is_resolved)
        #print('origin = \n%s' % self.origin)
        #print('T = \n%s' % self.T)

    def resolve_coord1(self, i, coord, xyz,
                       r, ref_coord_type, ref_coord):
        #print('  resolving cid=%s' % coord.cid)
        # get the reference transform
        T1 = self.T[r[0], :, :].reshape(3, 3)
        T2 = self.T[r[1], :, :].reshape(3, 3)
        T3 = self.T[r[2], :, :].reshape(3, 3)

        e = []
        T = [T1, T2, T3]
        for Ti, xyzi, ri, ref_coord_typei, ref_coordi in zip(T, xyz, r, ref_coord_type, ref_coord):
            #print("  xyzi = %s" % xyzi)
            #print("  Ti = %s" % Ti)
            #print("  ri = %s" % ri)
            #print("  refi = %s" % ref_coordi)
            if ref_coord_typei == 'R':  # Rectangular
                ei = dot(Ti, xyzi)
            elif ref_coord_typei == 'C':  # Cylindrical
                ei = dot(T, self.cylindrical_to_rectangular(xyzi))
                 #+ self.origin[ri]
            elif ref_coord_typei == 'S':  # Spherical
                ei = dot(Ti, self.spherical_to_rectangular(xyzi))
            else:
                raise NotImplementedError(ref_coord_typei)
            e.append(ei)
        #print("  e = %s" % e)
        e13 = e[2] - e[0]
        e12 = e[1] - e[0]

        coord.k = normalize(e12)
        coord.j = normalize(cross(coord.k, e13))
        coord.i = cross(coord.j, coord.k)

        #print("  e13 = %s" % e13)
        #print("  e12 = %s" % e12)
        #print("  coord.i = %s" % coord.i)
        #print("  coord.j = %s" % coord.j)
        #print("  coord.k = %s" % coord.k)
        coord.is_resolved = True
        self.is_resolved[i] = True
        self.origin[i] = e[0]
        #print('i** = %s' % i)
        #print('coord = \n%s' % coord)
        T = vstack([coord.i, coord.j, coord.k])
        #print('T.shape = %s' % str(T.shape))
        self.T[i, :, :] = T

        coord.origin = e[0]

    def resolve_coord2(self, i, coord,
                       r, ref_coord_type, ref_coord):
        #print('  icoord = %s' % i)
        #print('  ref_coord = %s' % ref_coord)
        #print('  ref_coord_type = %s' % ref_coord_type)
        #print('  ref_coord.i = %s' % ref_coord.i)
        #print('  ref_coord.j = %s' % ref_coord.j)
        #print('  ref_coord.k = %s' % ref_coord.k)
        #print('  ref_coord.o = %s' % ref_coord.origin)

        #e1 = coord.origin

        # get the reference transform
        T = self.T[r, :, :].reshape(3, 3)
        #e1 = coord.origin

        #print('  ref_coord_type = %s' % ref_coord_type)
        #print('  coord.e1 = %s' % coord.e1)
        #print('  coord.e2 = %s' % coord.e2)
        #print('  coord.e3 = %s' % coord.e3)
        #print("  T[%i] = %s" % (i, T))
        # transform coord.e2/e3 to the global
        if ref_coord_type == 'R':  # Rectangular
            e1 = dot(T, coord.e1)
            # proving math is right...
            pts = vstack([coord.e2, coord.e3])
            e123 = dot(T, pts.T).T
            e2 = e123[0, :]
            e3 = e123[1, :]
        elif ref_coord_type == 'C':  # Cylindrical
            _e1 = dot(T, self.cylindrical_to_rectangular(coord.e1)).T
            _e2 = dot(T, self.cylindrical_to_rectangular(coord.e2)).T
            _e3 = dot(T, self.cylindrical_to_rectangular(coord.e3)).T
            _e123 = vstack([_e1, _e2, _e3])

            pts = vstack([coord.e1, coord.e2, coord.e3])
            e123 = dot(T, self.cylindrical_to_rectangular(pts).T)
            assert array_equal(e123, _e123), "e123=\n%s\n\n_e123=\n%s" % (e123, _e123)
            e1 = e123[0, :]
            e2 = e123[1, :]
            e3 = e123[2, :]

        elif ref_coord_type == 'S':  # Spherical
            #print([coord.e1, coord.e2, coord.e3])
            # column
            #_e1 = dot(T, self.spherical_to_rectangular(coord.e1))
            #_e2 = dot(T, self.spherical_to_rectangular(coord.e2))
            #_e3 = dot(T, self.spherical_to_rectangular(coord.e3))
            #print("_e1", _e1
            #_e123 = vstack([_e1, _e2, _e3])

            pts = vstack([coord.e1, coord.e2, coord.e3])
            #print("pts", pts, pts.shape)
            if 1:
                _e1 = self.spherical_to_rectangular(coord.e1)  # good are _
                _e2 = self.spherical_to_rectangular(coord.e2)
                _e3 = self.spherical_to_rectangular(coord.e3)
                _e123 = vstack([_e1, _e2, _e3])
                e123 = _e123
                #t123 = self.spherical_to_rectangular(pts)
                #assert array_equal(t123, _t123), "t123=\n%s\n\n_t123=\n%s" % (t123, _t123)

                #e123 = dot(T, t123).T
                #e123 = dot(T, self.spherical_to_rectangular(pts).T).T
                #assert array_equal(e123, _e123), "e123=\n%s\n\n_e123=\n%s" % (e123, _e123)
            else:
                e123 = dot(T, self.spherical_to_rectangular(pts).T).T
            e1 = e123[0, :]
            e2 = e123[1, :]
            e3 = e123[2, :]

            #pts = vstack([coord.e2, coord.e3])
            #print("pts2", pts, pts.shape)
            #e123 = dot(T, self.spherical_to_rectangular(coord, pts).T).T
            #e2 = e123[0, :]
            #e3 = e123[1, :]
            #assert array_equal(e1, _e1), "e1=%s\n_e1=%s" % (e1, _e1)
            #assert allclose(e2, _e2), "\n e2=%s\n_e2=%s\ndiff=%s" % (e2, _e2, e2 - _e2)
            #assert allclose(e3, _e3), "\n e3=%s\n_e3=%s\ndiff=%s" % (e3, _e3, e3 - _e3)

            #e1 = dot(T, self.spherical_to_rectangular(coord, coord.e1))
            #e2 = dot(T, self.spherical_to_rectangular(coord, coord.e2))
            #e3 = dot(T, self.spherical_to_rectangular(coord, coord.e3))
            #raise NotImplementedError(ref_coord)
        else:
            raise RuntimeError(ref_coord_type)
        #e2, rid_matrix = self.rid.transformToGlobal(coord.e2)
        #e3, rid_matrix = self.rid.transformToGlobal(coord.e3)
        #print("  e1 = %s" % list(e1))
        #print("  e2 = %s" % list(e2))
        #print("  e3 = %s" % list(e3))
        e13 = e3 - e1
        e12 = e2 - e1
        #print('e13 = %s' % list(e13))
        #print('e12 = %s' % list(e12))
        coord.k = normalize(e12)
        #print('k = %s' % list(coord.k))
        coord.j = normalize(cross(coord.k, e13))
        #print('j = %s' % list(coord.j))
        coord.i = cross(coord.j, coord.k)

        #print("  e13 = %s" % e13)
        #print("  e12 = %s" % e12)
        #print("  coord.i = %s" % coord.i)
        #print("  coord.j = %s" % coord.j)
        #print("  coord.k = %s" % coord.k)
        coord.is_resolved = True
        self.is_resolved[i] = True
        self.origin[i] = e1
        self.T[i, :, :] = vstack([coord.i, coord.j, coord.k])
        coord.origin = e1

    def cylindrical_to_rectangular(self, r_theta_z):
        if len(r_theta_z.shape) == 2:
            R = r_theta_z[:, 0]
            theta = radians(r_theta_z[:, 1])
            z = r_theta_z[:, 2]
        else:
            R = r_theta_z[0]
            theta = radians(r_theta_z[1])
            z = r_theta_z[2]
        x = R * cos(theta)
        y = R * sin(theta)
        xyz = vstack([x, y, z])
        return xyz

    def rectangular_to_cylindrical(self, xyz):
        if len(xyz.shape) == 2:
            x = xyz[:, 0]
            y = xyz[:, 1]
            z = xyz[:, 2]
        else:
            x = xyz[0]
            y = xyz[1]
            z = xyz[2]
        theta = degrees(arctan2(y, x))
        R = sqrt(x * x + y * y)
        r_theta_z = vstack([R, theta, z])
        return r_theta_z

    def cylindrical_to_spherical(self, r_theta_z):
        xyz = self.cylindrical_to_rectangular(r_theta_z)
        r_theta_phi = self.rectangular_to_spherical(xyz)
        return r_theta_phi

    def spherical_to_cylindrical(self, r_theta_phi):
        xyz = self.spherical_to_rectangular(r_theta_phi)
        r_theta_z = self.rectangular_to_cylindrical(xyz)
        return r_theta_z

    def spherical_to_rectangular(self, r_theta_phi):
        if len(r_theta_phi.shape) == 2:
            #print('2d', r_theta_phi)
            R = r_theta_phi[:, 0]
            #print('  R = %s' % R)
            theta = radians(r_theta_phi[:, 1])
            phi = radians(r_theta_phi[:, 2])
            x = R * sin(theta) * cos(phi)
            y = R * sin(theta) * sin(phi)
            z = R * cos(theta)
            xyz = vstack([x, y, z]).T
            #print('  x = %s' % x)
            #print('  xyz = %s' % xyz)
        else:
            #dddf
            #print('1d')
            #print('  r_theta_phi = %s' % r_theta_phi)
            #print('  shape = %s' % str(r_theta_phi.shape))
            R = r_theta_phi[0]
            theta = radians(r_theta_phi[1])
            phi = radians(r_theta_phi[2])
            x = R * sin(theta) * cos(phi)
            y = R * sin(theta) * sin(phi)
            z = R * cos(theta)
            xyz = hstack([x, y, z])
        return xyz

    def rectangular_to_spherical(self, xyz):
        if len(xyz.shape) == 2:
            x = xyz[:, 0]
            y = xyz[:, 1]
            z = xyz[:, 2]
        else:
            x = xyz[0]
            y = xyz[1]
            z = xyz[2]
        R = sqrt(x * x + y * y + z * z)
        phi = degrees(arctan2(y, x))
        theta = degrees(arccos(z / R))
        r_theta_phi = vstack([R, theta, phi])
        return r_theta_phi

    def add_cord1r(self, card, comment=''):
        """adds a CORD1R card"""
        coord = CORD1R.add_card(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord1c(self, card, comment=''):
        """adds a CORD1C card"""
        coord = CORD1C.add_card(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord1s(self, card, comment=''):
        """adds a CORD1S card"""
        coord = CORD1S.add_card(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord2r(self, card, comment=''):
        """adds a CORD2R card"""
        coord = CORD2R.add_card(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1
        #print('adding cord2r; cids=%s' % self.coords.keys())

    def add_cord2x(self, card_name, card, comment=''):
        """adds a CORD2x card"""
        if card_name == 'CORD2R':
            self.add_cord2r(card, comment=comment)
        elif card_name == 'CORD2C':
            self.add_cord2c(card, comment=comment)
        elif card_name == 'CORD2S':
            self.add_cord2s(card, comment=comment)
        else:
            raise NotImplementedError(card_name)

    def add_cord2c(self, card, comment=''):
        """adds a CORD2C card"""
        coord = CORD2C.add_card(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord2s(self, card, comment=''):
        """adds a CORD2S card"""
        coord = CORD2S.add_card(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def transform_xyz_to_global_by_coord_id(self, xyz, cp):
        xyz = _check_xyz_shape(xyz)
        icp = self.get_coord_index_by_coord_id(cp)
        T = self.T[icp, :, :].reshape(3, 3)
        origin = self.origin[icp, :]
        xyz_global = dot(xyz - origin, transpose(T))
        return xyz_global

    def transform_node_id_to_local_by_coord_id(self, node_id, cp_local):
        xyz_global = self.transform_node_id_to_global_xyz(node_id)
        icp = self.get_coord_index_by_coord_id(cp_local)
        T = self.T[icp, :, :].reshape(3, 3)
        origin = self.origin[icp, :]
        xyz_local = dot(xyz_global - origin, transpose(T))
        return xyz_local

    def transform_node_id_to_global_xyz(self, node_id):
        coord = self.coords[0]

        grids = self.model.grid
        inode = grids.get_node_index_by_node_id(node_id)
        cp = grids.cp[inode]
        ucp = unique(cp)
        node_id = self._set_as_array(node_id)
        nnodes = len(node_id)
        #print('node_id = ', node_id)
        #print('inode = ', inode)

        xyz = zeros((nnodes, 3), dtype='float64')
        #print('cp = %s' % cp)
        #print('********')
        #print('ucp = %s' % ucp)
        for cpi in ucp:
            #print('-----------')
            #print('cpi = %s' % cpi)
            inodei = where(cpi == cp)[0]
            xyzi = grids.xyz[inodei, :]
            if cpi == 0:
                xyz[inodei, :] = xyzi
            else:
                icp = self.get_coord_index_by_coord_id(cpi)
                #print('icp = %s' % icp)
                T = self.T[icp, :, :].reshape(3, 3)
                origin = self.origin[icp, :]
                # dot(p - self.origin, transpose(beta))
                xyz[inodei, :] = dot(xyzi - origin, transpose(T))
        return xyz
        #return self.transform_vector_to_global(p) + self.origin

    def transform(self, cp):
        #print('cp = %s' % cp)
        return self.coords[cp].beta()

    def origin(self, cp):
        return self.coords[cp].origin

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        for i in self.n:
            yield i

    def values(self):
        for i in range(self.n):
            yield self.__getitem__(i)

    def items(self):
        for i in range(self.n):
            yield i, self.__getitem__(i)

def _check_xyz_shape(xyz):
    xyz = asarray(xyz)
    ndim = len(xyz.shape)
    if ndim == 1:
        xyz = xyz.reshape((1, 3))
    return xyz
