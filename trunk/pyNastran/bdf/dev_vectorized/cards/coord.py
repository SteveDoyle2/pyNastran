from numpy import (array, concatenate, searchsorted, unique, zeros, array, full,
                   nan, where, vstack, dot, cross, degrees, radians, arctan2,
                   cos, sin, hstack, array_equal, allclose, eye)
from numpy.linalg import norm

from pyNastran.bdf.cards.coordinateSystems import (
    CORD1R, CORD1C, CORD1S,
    CORD2R, CORD2C, CORD2S)

def normalize(v):
    print(v)
    return v / norm(v, axis=0)

# ..todo:: incomplete

class Coord(object):
    def get_global_position(self, xyz, cp):
        assert isinstance(cp, int), cp
        coord = self.coords[cp]

        T = coord.beta()
        origin = coord.origin.reshape(3, 1)
        print('origin = %s' % origin)
        print('T.shape=%s' % str(T.shape))
        print('xyz.shape=%s' % str(xyz.shape))
        xyz2 = (dot(T, xyz.T)  + origin).T
        #print('xyz = %s' % xyz.T)
        #print('xyz2 = %s' % xyz2)
        assert xyz.shape == xyz2.shape, "xyz.shape=%s xyz2.shape=%s" % (xyz.shape, xyz2.shape)
        return xyz2

    def __repr__(self):
        print('dummy')

    def __str__(self):
        return 'dummy'

    def __getitem__(self, value):
        return self.coords[value]

    def allocate(self, ncards=None, card_count=None):
        float_fmt = self.model.float
        assert ncards is not None or card_count is not None
        if ncards is None:
            ncards = array([card_count[name]
                      for name in ['CORD1R', 'CORD1C', 'CORD1S',
                                   'CORD2R', 'CORD2C', 'CORD2S']
                      if name in card_count], dtype='int32').sum() + 1
        print('nCOORDcards = %s' % ncards)
        #print('ncards coord = %s' % ncards)
        self.coord_id = zeros(ncards, dtype='int32')
        self.Type = full(ncards, nan, dtype='|S1')  # R-CORD2R, S-CORD2S, C-CORD2C
        self.T = full((ncards, 3, 3), nan, dtype=float_fmt)
        self.origin = full((ncards, 3), nan, dtype=float_fmt)
        self.is_resolved = full(ncards, False, dtype='bool')

    def __init__(self, model):
        """
        Defines the ShellProperties object.

        :param self: the ShellProperties object
        :param model: the BDF object
        :param pshells: the list of PSHELL cards
        :param pcomps: the list of PCOMP cards
        :param pshears: the list of PSHEAR cards
        """
        self.model = model
        float_fmt = self.model.float

        self.n = 1
        ncards = 1
        self.coords = {0: CORD2R(),}
        self.coord_id = zeros(ncards, dtype='int32')
        self.Type = full(ncards, 'R', dtype='|S1')  # R-CORD2R, S-CORD2S, C-CORD2C

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
        print('----------building COORDx-------------')
        cids_to_resolve = []
        ncoords = len(self.coords.keys())
        self.allocate(ncards=ncoords)
        #print('coord_ids = %s' % self.coords.keys())
        #print('T = \n%s' % self.T)
        for i, (cid, coord) in enumerate(sorted(self.coords.iteritems())):
            self.coord_id[i] = cid
            self.Type[i] = coord.Type
            if coord.isResolved:
                self.is_resolved[i] = True
                self.origin[i, :] = coord.origin
                self.T[i, :, :] = vstack([coord.i[:],
                                          coord.j[:],
                                          coord.k[:]] )
            else:
                print('need to resolve cid=%i rid=%i Type=%s' % (cid, coord.rid, coord.Type))
                cids_to_resolve.append(cid)
        print('coord_id = %s' % self.coord_id)
        #print('T =\n%s' % self.T)
        self.resolve_coords(cids_to_resolve)

    def resolve_coords(self, cids_to_resolve):
        while cids_to_resolve:
            cids_to_resolve2 = []
            print('need to resolve cids=%s' % cids_to_resolve)
            print('is_resolved=%s' % self.is_resolved)
            #i = where(self.is_resolved == False)[0]
            #cids_to_resolve = self.coord_id[i]
            for cid in cids_to_resolve:
                i = where(self.coord_id == cid)[0]
                coord = self.coords[cid]
                n = coord.type[-2]
                if n == '2':
                    rid = coord.rid
                    print('need to resolve rid=%s...' % rid)
                    j = where(self.coord_id == rid)[0]
                    is_resolved = self.is_resolved[j]
                    print('  is rid=%i resolved -> %s' % (rid, is_resolved))
                    if is_resolved:
                        ref_coord = self.coords[rid]
                        ref_coord_type = self.Type[j][0]
                        self.resolve_coord2(i, coord,
                                            j, ref_coord_type, ref_coord)
                elif n == '1':
                    g123 = array([coord.g1, coord.g2, coord.g3])
                    inode = searchsorted(self.model.grid.node_id, g123)
                    node_cids = self.model.grid.cp[inode]
                    xyz = self.model.grid.xyz[inode, :]
                    print('need to resolve cid=%s node_cps=%s...' % (cid, node_cids))
                    j = searchsorted(self.coord_id, node_cids)
                    is_resolved = self.is_resolved[j].all()
                    if is_resolved:
                        print('  resolved -> True is_resolved=%s' % is_resolved)
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
                    print('rid=%s is not resolved' % rid)
                    cids_to_resolve2.append(cid)
            assert len(cids_to_resolve) < cids_to_resolve2, 'circular reference...\ncoord_ids=%s\n' % (cids_to_resolve2)
            cids_to_resolve = cids_to_resolve2
            #break
        #print('is_resolved = %s' % self.is_resolved)
        #print('origin = \n%s' % self.origin)
        #print('T = \n%s' % self.T)
        #aaa

    def resolve_coord1(self, i, coord, xyz,
                       r, ref_coord_type, ref_coord):
        print('  resolving cid=%s' % coord.cid)
        # get the reference transform
        T1 = self.T[r[0], :, :].reshape(3, 3)
        T2 = self.T[r[1], :, :].reshape(3, 3)
        T3 = self.T[r[2], :, :].reshape(3, 3)

        e2 = []
        T = [T1, T2, T3]
        for Ti, xyzi, ri, ref_coord_typei, ref_coordi in zip(T, xyz, r, ref_coord_type, ref_coord):
            #print("  xyzi = %s" % xyzi)
            #print("  Ti = %s" % Ti)
            #print("  ri = %s" % ri)
            #print("  refi = %s" % ref_coordi)
            if ref_coord_typei == 'R':  # Rectangular
                ei2 = dot(Ti, xyzi)
            elif ref_coord_typei == 'C':  # Cylindrical
                ei2 = dot(T, self.cylindrical_to_rectangular(xyzi))
                 #+ self.origin[ri]
            elif ref_coord_typei == 'S':  # Spherical
                ei2 = dot(Ti, self.spherical_to_rectangular(xyzi))
            else:
                raise NotImplementedError(ref_coord_typei)
            e2.append(ei2)
        print("  e2 = %s" % e2)
        e13 = e2[2] - e2[0]
        e12 = e2[1] - e2[0]

        coord.k = normalize(e12)
        coord.j = normalize(cross(coord.k, e13))
        coord.i = cross(coord.j, coord.k)

        print("  e13 = %s" % e13)
        print("  e12 = %s" % e12)
        print("  coord.i = %s" % coord.i)
        print("  coord.j = %s" % coord.j)
        print("  coord.k = %s" % coord.k)
        coord.isResolved = True
        self.is_resolved[i] = True
        self.origin[i] = e2[0]
        print('i** = %s' % i)
        print('coord = \n%s' % coord)
        T = vstack([coord.i, coord.j, coord.k])
        print('T.shape = %s' % str(T.shape))
        self.T[i, :, :] = T

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

        #print("  coord.e1 = %s" % coord.e1)
        #print("  coord.e2 = %s" % coord.e2)
        #print("  coord.e3 = %s" % coord.e3)
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
            e1 = dot(T, self.cylindrical_to_rectangular(coord.e1))
            e2 = dot(T, self.cylindrical_to_rectangular(coord.e2))
            e3 = dot(T, self.cylindrical_to_rectangular(coord.e3))
        elif ref_coord_type == 'S':  # Spherical
            #print [coord.e1, coord.e2, coord.e3]
            # column
            #_e1 = dot(T, self.spherical_to_rectangular(coord.e1))
            #_e2 = dot(T, self.spherical_to_rectangular(coord.e2))
            #_e3 = dot(T, self.spherical_to_rectangular(coord.e3))
            #print "_e1", _e1
            #_e123 = vstack([_e1, _e2, _e3])

            pts = vstack([coord.e1, coord.e2, coord.e3])
            #print "pts", pts, pts.shape
            #_t1 = self.spherical_to_rectangular(coord.e1)  # good are _
            #_t2 = self.spherical_to_rectangular(coord.e2)
            #_t3 = self.spherical_to_rectangular(coord.e3)
            #_t123 = vstack([_t1, _t2, _t3])
            #t123 = self.spherical_to_rectangular(pts)
            #assert array_equal(t123, _t123), "t123=\n%s\n\n_t123=\n%s" % (t123, _t123)

            #e123 = dot(T, t123).T
            e123 = dot(T, self.spherical_to_rectangular(pts).T).T
            #assert array_equal(e123, _e123), "e123=\n%s\n\n_e123=\n%s" % (e123, _e123)
            e1 = e123[0, :]
            e2 = e123[1, :]
            e3 = e123[2, :]

            #pts = vstack([coord.e2, coord.e3])
            #print "pts2", pts, pts.shape
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
        #e2, rid_matrix = self.rid.transformToGlobal(coord.e2)
        #e3, rid_matrix = self.rid.transformToGlobal(coord.e3)
        #print("  e2 = %s" % e2)
        #print("  e3 = %s" % e3)
        e13 = e3 - e1
        e12 = e2 - e1
        coord.k = normalize(e12)
        coord.j = normalize(cross(coord.k, e13))
        coord.i = cross(coord.j, coord.k)

        #print("  e13 = %s" % e13)
        #print("  e12 = %s" % e12)
        #print("  coord.i = %s" % coord.i)
        #print("  coord.j = %s" % coord.j)
        #print("  coord.k = %s" % coord.k)
        coord.isResolved = True
        self.is_resolved[i] = True
        self.origin[i] = e1
        self.T[i, :, :] = vstack([coord.i, coord.j, coord.k])

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
        theta = degrees(atan2(y, x))
        R = sqrt(x * x + y * y)
        r_theta_z = vstack([R, theta, z])
        return r_theta_z

    def spherical_to_rectangular(self, r_theta_phi):
        if len(r_theta_phi.shape) == 2:
            #print '2d', r_theta_phi
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
        phi = degrees(atan2(y, x))
        theta = degrees(acos(z / R))
        r_theta_phi = vstack([R, theta, phi])
        return r_theta_phi

    def add_cord1r(self, card, comment=''):
        coord = CORD1R(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord1c(self, card, comment=''):
        coord = CORD1C(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord1s(self, card, comment=''):
        coord = CORD1S(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord2r(self, card, comment=''):
        coord = CORD2R(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1
        #print('adding cord2r; cids=%s' % self.coords.keys())

    def add_cord2c(self, card, comment=''):
        coord = CORD2C(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def add_cord2s(self, card, comment=''):
        coord = CORD2S(card, comment=comment)
        self.coords[coord.cid] = coord
        self.n += 1

    def transform(self, cp):
        #print('cp = %s' % cp)
        return self.coords[cp].beta()

    def origin(self, cp):
        return self.coords[cp].origin

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        for cid in self.coords:
            yield cid

    def values(self):
        for coord in self.coords.itervalues():
            yield coord

    def items(self):
        for cid, coord in self.coords.iteritems():
            yield cid, coord

    def iterkeys(self):
        return self.keys()
    def itervalues(self):
        return self.values()
    def iteritems(self):
        return self.items()