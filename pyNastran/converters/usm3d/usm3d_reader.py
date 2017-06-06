"""
Defines:
  - Usm3d(log=None, debug=False)
     - read_cogsg(cogsg_file, stop_after_header=False)
     - read_usm3d(self, basename, dimension_flag, read_loads=True)
     - write_usm3d(basename)
     - read_mapbc(mapbc_filename)
     - read_bc(self, bc_filename, stop_after_header=False, get_lbouf=False)
     - read_flo(self, flo_filename, n=None, node_ids=None)
"""
from __future__ import print_function
import os
from struct import pack, unpack

from six.moves import range
from numpy import array, zeros, where
import numpy as np
from pyNastran.utils.log import get_logger2


class Usm3d(object):
    bcmap_to_bc_name = {
        0 : 'Supersonic Inflow',
        1 : 'Reflection plane',
        2 : 'Supersonic Outflow',
        3 : 'Subsonic Outer Boundaries',
        4 : 'Viscous Surfaces',
        5 : 'Inviscid aerodynamic surface',
        44 : 'Blunt base',
        55 : 'Thick Trailing Edges',

        103 : 'Engine-exhaust (Fan)',
        102 : 'Engine-exhaust (Jet Core)',
        101 : 'Engine-intake',
        203 : 'Engine-exhaust (Fan)',
        202 : 'Engine-exhaust (Jet Core)',
        201 : 'Engine-intake',

        1001 : 'Special inflow',
        1002 : 'Special Outflow (Fixed Pressure)',

        #'0 - Freestream - Supersonic Inflow (Bounding Box)
        #'2 - Extrapolation - Supersonic Outflow (Bounding Box)

        #'1 - Reflection Plane - Tangent Flow - (Symmetry Plane)
        #'3 - Characteristic Inflow - Subsonic Inflow/Outflow/Sideflow (Bounding Box)

        #'4 - Inviscid Surface (Physical)
        #'5', Viscous Surface (Physical)
    }

    def __init__(self, log=None, debug=None):
        self.nodes = None
        self.tris = None
        self.tets = None
        self.bcs = None
        self.precision = 'double'
        self.log = get_logger2(log, debug=debug)


    def read_mapbc(self, mapbc_filename):
        """
        0 - Supersonic Inflow
        1 - Reflection plane
        2 - Supersonic Outflow
        3 - Subsonic Outer Boundaries
        4 - Viscous Surfaces
        5 - Inviscid aerodynamic surface
        44 - Blunt base
        55 - Thick Trailing Edges
        n*100+3 - Engine-exhaust (Fan)
        n*100+2 - Engine-exhaust (Jet Core)
        n*100+1 - Engine-intake
        1001 - Special inflow
        1002 - Special Outflow (Fixed Pressure)

        0 - Freestream - Supersonic Inflow (Bounding Box)
        2 - Extrapolation - Supersonic Outflow (Bounding Box)

        1 - Reflection Plane - Tangent Flow - (Symmetry Plane)
        3 - Characteristic Inflow - Subsonic Inflow/Outflow/Sideflow (Bounding Box)

        4 - Inviscid Surface (Physical)
        5 - Viscous Surface (Physical)
        #==============================

        #Thu Dec 19 11:46:03 2013
        #bc.map
        Patch #    BC   Family   #surf   surfIDs      Family
        #--------------------------------------------------------
        1          44     44         0         0      Base  -> Blunt base
        2          4       4         0         0      Bay   -> Viscous Surfaces
        3          0       0         0         0      Inflow -> Supersonic Inflow
        4          2       2         0         0      Outflow -> Supersonic Outflow
        5          3       3         0         0      Sideflow -> Characteristic Inflow/Outflow
        7          1       1         0         0      Symmetry -> Reflection plane
        """
        with open(mapbc_filename, 'r') as mapbc_file:
            lines = mapbc_file.readlines()

        lines2 = []
        for line in lines:
            if len(line.strip().split('#')[0]) > 0:
                lines2.append(line)
        lines = lines2

        mapbc = {}
        for line in lines[1:]:
            sline = line.split()
            #self.log.info(sline)
            patch_id, bc, family, surf, surf_ids = sline[:5]
            mapbc[int(patch_id)] = [int(bc), int(family), int(surf), surf_ids]
        return mapbc

    def read_usm3d(self, basename, dimension_flag, read_loads=True):
        """
        Parameters
        ----------
        basename : str
            the root path to the *.cogsg, *.bc, *.mapbc, *.face, *.front files
        dimension_flag : int
            ???
            2?/3
        read_loads : bool; default=True
            ???

        Returns
        -------
        nodes : (nnodes, 3) float ndarray
           the xyz coordinates
        tris_tets : ???
           ???
        tris : (ntri, 3) int ndarray
           0-based node indices
        bcs : (ntris, ) int ndarray
           1-based patch/"property" IDs
        mapbc : dict[patch_id] : map_line
           patch_id : int
           map_line : List[bc, family, surf, surf_ids]
               bc : int
                   boundary condition number
               family : int
                   family number
               surf : int
                   surface number
               surf_ids : str
                   ???
        loads : dict[]
            ???
        flo_filename : str
            the latest result filename
            None : no *.flo file could be found
            static : *.flo
            transient : *_xxx.flo
        """
        cogsg_filename = basename + '.cogsg'
        bc_filename = basename + '.bc'
        face_filename = basename + '.face'
        front_filename = basename + '.front'
        mapbc_filename = basename + '.mapbc'
        flo_filename = None

        if 1:
            # pick the highest N value or use "basename.flo"
            dirname = os.path.dirname(basename)
            if dirname == '':
                dirname = os.getcwd()
            flo_filenames = os.listdir(dirname)

            # get the max N value
            nmax = -1
            for flo_filename in flo_filenames:
                base, ext = os.path.splitext(flo_filename)
                if ext == '.flo':
                    n = base.split('_')[-1]
                    try: # get the incrementation index
                        n = int(n)
                        nmax = max(n, nmax)
                    except ValueError: # don't bother incrementing
                        pass
            # determine .flo file name
            if nmax > 0:
                flo_filename = basename + '_%s.flo' % (nmax)
            else:
                flo_filename = basename + '.flo'
        else:
            # hardcoded flo file
            flo_filename = basename + '_160.flo'

        nodes, elements = self.read_cogsg(cogsg_filename)
        try:
            header, tris, bcs = self.read_bc(bc_filename)
        except IOError:
            tris = None
            bcs = None
            self.log.error('Cannot find %r...skipping; required for geometry' % bc_filename)
        try:
            mapbc = self.read_mapbc(mapbc_filename)
        except IOError:
            mapbc = {}
            self.log.warning('Cannot find %r...skipping' % mapbc_filename)
        self.tris = tris
        self.bcs = bcs
        self.mapbc = mapbc

        loads = {}
        if read_loads and os.path.exists(flo_filename):
            npoints = nodes.shape[0]
            try:
                node_ids_volume, loads = self.read_flo(flo_filename, n=npoints)
            except:
                self.log.error('Had trouble reading %r...' % flo_filename)
                raise
        else:
            self.log.warning('Cannot find %r...skipping' % flo_filename)
        self.loads = loads

        return nodes, elements, tris, bcs, mapbc, loads, flo_filename
        #self.read_front(front_file)
        #self.read_face(face_file)


    def write_usm3d(self, basename):
        write_usm3d_volume(self, basename)

    def read_bc(self, bc_filename, stop_after_header=False, get_lbouf=False):
        self.log.info("bc_filename = %r" % bc_filename)
        with open(bc_filename, 'r') as bc_file:
            lines = bc_file.readlines()

        #mbouf,dum1,dum1,igrid
        header = lines[0].strip().split()

        (nbouf, nbou1, npatch, igrid) = header
       #(ntris, nbou1, npatch, igrid) = header
        header[0] = int(nbouf)
        header[1] = int(nbou1)
        header[2] = int(npatch)
        header[3] = int(igrid)
        ntris = int(nbouf)

        if stop_after_header:
            return header, None, None

        if get_lbouf:
            lbouf = zeros((ntris, 4), dtype='int32')
            for i in range(ntris):
                line = lines[i+2].strip()
                #print('%r' % line)
                (n, isurf, n1, n2, n3) = line.split()
                lbouf[i, :] = [isurf, n1, n2, n3]
            return header, lbouf
        else:
            tris = zeros((ntris, 3), dtype='int32')
            bcs = zeros(ntris, dtype='int32')

            for i in range(ntris):
                (n, isurf, n1, n2, n3) = lines[i+2].split()
                tris[i] = [n1, n2, n3]
                bcs[i] = isurf
            tris = tris - 1
            #self.bcs = [tris, bcs]
            return header, tris, bcs

    def read_cogsg(self, cogsg_filename, stop_after_header=False):
        """
        Reads the *.cogsg file

        Returns
        -------
        nodes : (N, 3) float ndarray
           the nodes
        tet_elements : ???
           ???
        """
        cogsg_file = open(cogsg_filename, 'rb')

        # nelements * 4 * 4 + 32 ???
        dummy = cogsg_file.read(4)  # 1022848
        dummy_int, = unpack('>i', dummy)
        #assert dummy_int == 1022848, 'dummy_int = %s' % dummy_int

        # file header
        if self.precision == 'single':
            Format = '>6if'
            nbytes = 6 * 4 + 4
        elif self.precision == 'double':
            Format = '>6id'
            nbytes = 6 * 4 + 8
        else:
            raise RuntimeError('invalid precision format')
        data = cogsg_file.read(nbytes)

        (inew, ne, np, nb, npv, nev, tc) = unpack(Format, data)
        self.header = {
            'dummy'    : dummy_int,
            'inew'     : inew, # dummy int
            'nElements': ne,  # nc;  number of tets
            'nPoints'  : np,  # npo; number of grid points including nbn
            'nBoundPts': nb,  # nbn; number of boundary points including nbc
            'nViscPts' : npv, # npv; number of viscous points (=0 for Euler)
            'nViscElem': nev, # ncv; number of viscous cells (=0 for Euler)
            'tc'       : tc,  # dummy double
                              # nbc
        }
        if stop_after_header:
            return self.header
        self.log.info(self.header)

        # nbn nodes
        #
        #del ne, np

        if 1:
            nodes, tets = self._read_cogsg_volume(cogsg_file)
            return nodes, tets
        #else:
        #----------------------------------------------------------------------
        # elements
        # face elements
        nnodes_per_face = 3
        nfaces = ne

        if nfaces > 0:
            data_length = nnodes_per_face * nfaces
            Format = '>' + 'i' * data_length
            data = cogsg_file.read(4 * data_length)

            faces = unpack(Format, data)
            faces = array(faces)
            faces = faces.reshape((nfaces, 3))
        else:
            faces = None
        #----------------------------------------------------------------------
        # nodes
        nbound_pts = nb
        nnodes = nbound_pts

        #data_length = nnodes
        if self.precision == 'double':
            data_length = 8 * nnodes
        elif self.precision == 'single':
            data_length = 4 * nnodes
        else:
            raise RuntimeError('precision = %r' % self.precision)

        skip_nodes = False
        if skip_nodes:
            t = cogsg_file.tell()
            cogsg_file._goto(t + data_length * 3)
            nodes = None
        else:
            if self.precision == 'double':
                Format = '>%sd' % nnodes
                node_array_format = 'float64'
            elif self.precision == 'single':
                Format = '>%sd' % nnodes
                node_array_format = 'float32'
            else:
                raise RuntimeError('precision = %r' % self.precision)

            if 0:
                data = cogsg_file.read(data_length)
                X = unpack(Format, data)
                data = cogsg_file.read(data_length)
                Y = unpack(Format, data)
                data = cogsg_file.read(data_length)
                Z = unpack(Format, data)
                nodes = np.array([X, Y, Z])

                nodes = np.zeros((nnodes, 3), node_array_format)
                nodes[:, 0] = X
                nodes[:, 1] = Y
                nodes[:, 2] = Z
                del X, Y, Z
            else:
                data = cogsg_file.read(3 * data_length)
                assert self.precision == 'single', self.precision
                nodes = np.fromstring(data, '>4f').reshape(3, nnodes).T
                #nodes = np.fromstring(data, '>4f').reshape(nnodes, 3)


        cogsg_file.read(nnodes * 3 * 8)  # 3 -> xyz, 8 -> double precision ???

        #----------------------------------------------------------------------
        # elements
        # boundary layer elements
        nnodes_per_tet = 4
        ntets = nev

        if ntets:
            data_length = nnodes_per_tet * ntets
            Format = '>' + 'i' * data_length
            data = cogsg_file.read(4 * data_length)

            tets = unpack(Format, data)
            tets = np.array(tets)
            tets = tets.reshape((tets, 4))

        #----------------------------------------------------------------------
        # volume points
        nnodes = npv

        Format = '>%si' % nnodes
        data = cogsg_file.read(4 * nnodes)

        nodes_vol = unpack(Format, data)
        nodes_vol = np.array(nodes_vol)
        nodes_vol = nodes_vol.reshape((tets, 3))

    def _read_cogsg_volume(self, cogsg_file):
        # volume cells
        self.log.debug('tell volume = %s' % cogsg_file.tell())
        # surface + volume cells ???
        nelements = self.header['nElements']
        Format = '>%si' % nelements


        self.log.debug("fv.tell = %s" % cogsg_file.tell())
        use_fromstring = True

        if use_fromstring:
            ndata = 4 * (4 * nelements)
            data = cogsg_file.read(ndata)

            # the 4 means that we make a (nelements, 4) array?
            elements = np.fromstring(data, dtype='>4f') - 1
        else:
            elements = np.zeros((nelements, 4), 'int32')
            for i in range(4): #  tets
                data = cogsg_file.read(4 * nelements)
                elements[:, i] = unpack(Format, data)
            elements -= 1

        assert elements.shape == (nelements, 4), elements.shape

        dummy2 = cogsg_file.read(4)
        self.log.debug("dummy2 = %s %s" % (unpack('>i', dummy2), unpack('>f', dummy2)))
        dummy_int2, = unpack('>i', dummy2)

        # 32 = dummy_int2 - 4 * nelements * 4
        assert self.header['dummy'] == dummy_int2

        #-----------------------------------
        # nodes
        nnodes = self.header['nPoints']
        Format = '>%sd' % nnodes

        dummy3 = cogsg_file.read(4)  # nnodes * 3 * 8
        dummy3_int, = unpack('>i', dummy3)
        #assert dummy3_int == 298560
        self.log.debug("dummy3 = %i" % unpack('>i', dummy3)) #, unpack('>f', dummy3)

        if use_fromstring:
            data_length = 8 * nnodes
            data = cogsg_file.read(3 * data_length)
            assert self.precision == 'double', self.precision
            nodes = np.fromstring(data, '>d').reshape(3, nnodes).T

            # the ravel creates a copy that we can then use to put in
            # a contigous order
            nodes = np.asarray(nodes.ravel(), dtype='<d').reshape(nnodes, 3)
        else:
            nodes = zeros((nnodes, 3), 'float64')
            for i in range(3): #  x, y, z
                data = cogsg_file.read(8 * nnodes)
                assert len(data) == (8 * nnodes)
                nodes[:, i] = unpack(Format, data)


        dummy4 = cogsg_file.read(4) # nnodes * 3 * 8
        dummy4_int, = unpack('>i', dummy4)
        #print("dummy4 = ", unpack('>i', dummy4), unpack('>f', dummy4))

        assert dummy3_int == dummy4_int
        self.nodes = nodes
        self.tets = elements
        return nodes, elements

    def read_flo(self, flo_filename, n=None, node_ids=None):
        """
        ipltqn is a format code where:
         - ipltqn = 0  (no printout)
         - ipltqn = 1  (unformatted)
         - ipltqn = 2  (formatted) - default

        Parameters
        ----------
        flo_filename : str
            the name of the file to read
        n : int; default=None
            the number of points to read (initializes the array)
            n is typically the number of points, but is not required to be
            this lets you read nodes 1...n, but not greater than n+1.
            node_ids must be set to None.
        node_ids : List[int]; default=None
            the specific points to read (n must be set to None).

        nvars = 5
          - (nodeID,rho,rhoU,rhoV,rhoW) = sline
            (e) = line

        nvars = 6
          - (nodeID,rho,rhoU,rhoV,rhoW,e) = line

        Also, stupid Nastran-esque float formatting is sometimes used,
        so 5.0-100 exists, which is 5.0E-100.  We just assume it's 0.

        Returns
        -------
        node_id : (nnodes, ) int ndarray
            the node ids for the selected loads
        loads : dict[result_name] : data
            result_name : str
                the name of the result
            data : (nnodes, ) float ndarray
                the data corresponding to the result_name
        """
        result_names = ['Mach', 'U', 'V', 'W', 'T', 'rho', 'rhoU', 'rhoV', 'rhoW', 'p', 'Cp']

        is_sparse = None
        if n is None:
            assert node_ids is not None, node_ids
            assert len(node_ids) > 0, node_ids
            n = len(node_ids)
            is_sparse = True
        else:
            assert node_ids is None, node_ids
            is_sparse = False

        #formatCode = 2
        node_id = np.zeros(n, 'int32')
        rho = np.zeros(n, 'float32')
        rhoU = np.zeros(n, 'float32')
        rhoV = np.zeros(n, 'float32')
        rhoW = np.zeros(n, 'float32')
        e = np.zeros(n, 'float32')

        with open(flo_filename, 'r') as flo_file:
            line = flo_file.readline().strip()
            try:
                #file is messsed up
                mach = float(line)
            except:
                raise
                #loads['Cp'] = e  # it's 0 anyways...
                #return node_id, loads


            # determine the number of variables on each line
            sline1 = flo_file.readline().strip().split()
            nvars = None
            if len(sline1) == 6:
                nvars = 6
                rhoi, rhoui, rhovi, rhowi, ei = Float(sline1[1:], 5)
            else:
                nvars = 5
                rhoi, rhoui, rhovi, rhowi = Float(sline1[1:], 4)
                sline2 = flo_file.readline().strip().split()
                ei = Float(sline2, 1)[0]

            # set the i=0 values
            if not is_sparse:
                nmax = n
                i = 0
                node_id[i] = sline1[0]
                rho[i] = rhoi
                rhoU[i] = rhoui
                rhoV[i] = rhovi
                rhoW[i] = rhowi
                e[i] = ei
            else:
                ni = 0
                node_ids_minus_1 = array(node_ids) - 1
                nmax = node_ids_minus_1.max() + 1
                if 0 in node_ids_minus_1:
                    i = 0
                    node_id[i] = sline1[0]
                    rho[i] = rhoi
                    rhoU[i] = rhoui
                    rhoV[i] = rhovi
                    rhoW[i] = rhowi
                    e[i] = ei
                    ni += 1

            # loop over the rest of the data in the flo file
            if node_ids is None:
                ni = n
                # extract nodes 1, 2, ... 10, but not 11+
                if nvars == 6:  # sequential nvars=6
                    for i in range(1, n):
                        sline1 = flo_file.readline().strip().split()
                        rhoi, rhoui, rhovi, rhowi, ei = Float(sline1[1:], 5)
                        node_id[i] = sline1[0]
                        rho[i] = rhoi
                        rhoU[i] = rhoui
                        rhoV[i] = rhovi
                        rhoW[i] = rhowi
                        e[i] = ei
                        assert len(sline1) == 6, 'len(sline1)=%s' % len(sline1)
                else:  # sequential nvars=5
                    for i in range(1, n):
                        sline1 = flo_file.readline().strip().split()
                        rhoi, rhoui, rhovi, rhowi = Float(sline1[1:], 4)
                        assert len(sline1) == 5, 'len(sline1)=%s' % len(sline1)

                        sline2 = flo_file.readline().strip().split()
                        ei = Float(sline2, 1)[0]

                        node_id[i] = sline1[0]
                        rho[i] = rhoi
                        rhoU[i] = rhoui
                        rhoV[i] = rhovi
                        rhoW[i] = rhowi
                        e[i] = ei
                        assert len(sline2) == 1, 'len(sline2)=%s' % len(sline2)
            else:
                # extract node 1, 2, and 10
                if nvars == 6:  # dynamic nvars=6
                    for i in range(1, nmax):
                        if i in node_ids_minus_1:
                            sline1 = flo_file.readline().strip().split()
                            rhoi, rhoui, rhovi, rhowi, ei = Float(sline1[1:], 5)

                            node_id[ni] = sline1[0]
                            rho[ni] = rhoi
                            rhoU[ni] = rhoui
                            rhoV[ni] = rhovi
                            rhoW[ni] = rhowi
                            e[ni] = ei
                            assert len(sline1) == 6, 'len(sline1)=%s' % len(sline1)
                            ni += 1
                        else:
                            line1 = flo_file.readline()
                else:  # dynamic nvars=5
                    for i in range(1, nmax):
                        if i in node_ids_minus_1:
                            sline1 = flo_file.readline().strip().split()
                            rhoi, rhoui, rhovi, rhowi = Float(sline1[1:], 4)
                            assert len(sline1) == 5, 'len(sline1)=%s' % len(sline1)

                            sline2 = flo_file.readline().strip().split()
                            ei = Float(sline2, 1)[0]

                            node_id[ni] = sline1[0]
                            rho[ni] = rhoi
                            rhoU[ni] = rhoui
                            rhoV[ni] = rhovi
                            rhoW[ni] = rhowi
                            e[ni] = ei
                            assert len(sline2) == 1, 'len(sline2)=%s' % len(sline2)
                            ni += 1
                        else:
                            line1 = flo_file.readline()
                            line2 = flo_file.readline()

        assert len(rho) == ni

        # llimit the minimum density (to prevent division errors)
        rho_min = 0.001
        irho_zero = where(rho < rho_min)[0]
        rho[irho_zero] = rho_min

        loads = {}

        if '.aux.' in flo_filename:
            # the names (rho, e, rhoU, etc.) aren't correct, but that's OK
            # the load names are correct
            loads['inst vor'] = rho
            loads['timeavg vor'] = rhoU
            loads['inst visc'] = rhoV
            loads['timeavg visc'] = rhoW
            loads['local CFL'] = e
            return node_id, loads

        # standard outputs
        gamma = 1.4
        two_over_Mach2 = 2.0 / mach ** 2
        one_over_gamma = 1.0 / gamma

        gm1 = gamma - 1

        # node_id, rhoi, rhoui, rhovi, rhowi, ei
        rhoVV = (rhoU**2 + rhoV**2 + rhoW**2) / rho
        if 'p' in result_names or 'Mach' in result_names or 'Cp' in result_names:
            pND = gm1*(e - rhoVV/2.)
            if 'p' in result_names:
                loads['p'] = pND
        if 'Mach' in result_names:
            Mach = (rhoVV / (gamma * pND))**0.5
            loads['Mach'] = Mach
        if 'Cp' in result_names:
            Cp = two_over_Mach2 * (pND - one_over_gamma)
            loads['Cp'] = Cp

        T = gamma * pND / rho # =a^2 as well
        #a = T.sqrt()
        if 'T' in result_names:
            loads['T'] = T

        if 'rho' in result_names:
            loads['rho'] = rho

        if 'rhoU' in result_names:
            loads['rhoU'] = rhoU
        if 'rhoV' in result_names:
            loads['rhoV'] = rhoV
        if 'rhoW' in result_names:
            loads['rhoW'] = rhoW

        if 'U' in result_names:
            loads['U'] = rhoU / rho
        if 'V' in result_names:
            loads['V'] = rhoV / rho
        if 'W' in result_names:
            loads['W'] = rhoW / rho
        return node_id, loads




def Float(sline, n):
    """floats a value"""
    vals = []
    for val in sline:
        try:
            vals.append(float(val))
        except:
            vals.append(0.0)
    return vals

def write_usm3d_volume(model, basename):
    """
    writes a *.cogsg, *.front, *.face file
    """
    cogsg_file = basename + '.cogsg'
    #face_file = basename + '.face'
    #front_file = basename + '.front'

    write_cogsg_volume(model, cogsg_file)
    #write_front(model, front_file)
    #write_face(model, face_file)


#def write_front(model):
    #pass

#def write_face(model):
    #pass

def write_cogsg_volume(model, cogsg_file):
    """
    writes a *.cogsg file
    """
    #n = 0
    self = model
    nnodes = self.nodes.shape[0]
    ntets = self.tets.shape[0]

    with open(cogsg_file, 'wb') as outfile:
        # file header
        values = [32 + ntets * 4 * 4,]
        block_size = pack('>i', *values)
        outfile.write(block_size)

        header = self.header
        values_header = [
            header['inew'],
            header['nElements'],
            header['nPoints'],
            header['nBoundPts'],
            header['nViscPts'],
            header['nViscElem'],
            header['tc'], # d
            ]
        #n = 36
        Format = '>6id'
        model.log.debug("Format = %r" % Format)
        data = pack(Format, *values_header)
        outfile.write(data)
        #print("outfile.tell = %s" % outfile.tell())
        #outfile.write(block_size)
        #--------------------------------------------------------------------------
        tets = self.tets

        # tet header
        #values = [ntets * 4 * 4]  # n1, n2, n3, n4 -> 4; 4 -> int
        #block_size = pack('>i', *values)
        #outfile.write(block_size)

        # tets
        Format = '>%si' % ntets
        n0 = tets[:, 0] + 1
        n1 = tets[:, 1] + 1
        n2 = tets[:, 2] + 1
        n3 = tets[:, 3] + 1
        #print("n0 = %s" % n0)
        outfile.write(pack(Format, *n0))
        outfile.write(pack(Format, *n1))
        outfile.write(pack(Format, *n2))
        outfile.write(pack(Format, *n3))
        #n += 4 * 8 * ntets
        #print("outfile.tell 2 = ", outfile.tell(), n)

        # tet footer
        outfile.write(block_size)
        #--------------------------------------------------------------------------

        # nodes header
        values = [nnodes * 3 * 8]  # xyz -> 3; 8 -> double precision (???)
        block_size = pack('>i', *values)
        outfile.write(block_size)

        # nodes
        #npoints = header['nPoints']
        Format = '>%sd' % nnodes
        nodes = self.nodes
        n0 = nodes[:, 0]
        n1 = nodes[:, 1]
        n2 = nodes[:, 2]
        outfile.write(pack(Format, *n0))
        outfile.write(pack(Format, *n1))
        outfile.write(pack(Format, *n2))

        # nodes footer
        outfile.write(block_size)


def main():  # pragma: no cover
    """test problem"""
    model = Usm3d()
    if 1:
        #basename = 'HSCT_inviscid'
        #basename = 'box'
        basename = 'new2'
        model.read_usm3d(basename, 3)
        model.write_usm3d(basename + '_2')

        model.read_usm3d(basename + '_2', 3)
    else:
        basename = 'new2'

        #cogsg_filename = basename + '.cogsg'
        #bc_filename = basename + '.bc'
        #face_filename = basename + '.face'
        #front_filename = basename + '.front'
        #mapbc_filename = basename + '.mapbc'
        flo_filename = basename + '.flo'

        #model.read_usm3d(basename, 3)

        node_ids, loads = model.read_flo(flo_filename, node_ids=[10])


if __name__ == '__main__':  # pragma: no cover
    main()
