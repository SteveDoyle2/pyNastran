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
import os
from struct import pack, unpack
from collections import OrderedDict

import numpy as np
from pyNastran.utils import check_path
from cpylog import get_logger2


def read_usm3d(basename, log=None, debug=None):
    """reads a usm3d file"""
    model = Usm3d(log=log, debug=debug)
    #model.read_cogsg(cogsg_filename, stop_after_header=False)
    unused_dimension_flag = None
    model.read_usm3d(basename, unused_dimension_flag, read_loads=True)
    return model

class Usm3d:
    """Usm3d interface class"""
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
        """
        Initializes the Usm3d object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        """
        self.nodes = None
        self.tris = None
        self.tets = None
        self.bcs = None

        self.header = None
        self.loads = None
        self.mapbc = None
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

    def read_usm3d(self, basename, unused_dimension_flag, read_loads=True):
        """
        Parameters
        ----------
        basename : str
            the root path to the *.cogsg, *.bc, *.mapbc, *.face, *.front files
        unused_dimension_flag : int; unused
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
        unused_face_filename = basename + '.face'
        unused_front_filename = basename + '.front'
        mapbc_filename = basename + '.mapbc'
        flo_filename = None

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

        nodes, elements = self.read_cogsg(cogsg_filename)

        try:
            unused_header, tris, bcs = self.read_bc(bc_filename)
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
                unused_node_ids_volume, loads = self.read_flo(flo_filename, n=npoints)
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
        """
        writes a *.cogsg, *.front, *.face file
        """
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
            lbouf = np.zeros((ntris, 4), dtype='int32')
            for i in range(ntris):
                line = lines[i+2].strip()
                #print('%r' % line)
                (unused_n, isurf, n1, n2, n3) = line.split()
                lbouf[i, :] = [isurf, n1, n2, n3]
            return header, lbouf

        tris = np.zeros((ntris, 3), dtype='int32')
        bcs = np.zeros(ntris, dtype='int32')
        for i in range(ntris):
            (unused_n, isurf, n1, n2, n3) = lines[i+2].split()
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
        check_path(cogsg_filename, 'cogsg file')
        with open(cogsg_filename, 'rb') as cogsg_file:
            # nelements * 4 * 4 + 32 ???
            dummy = cogsg_file.read(4)  # 1022848
            dummy_int, = unpack('>i', dummy)
            #assert dummy_int == 1022848, 'dummy_int = %s' % dummy_int

            # file header
            if self.precision == 'single':
                sformat = '>6if'
                nbytes = 6 * 4 + 4
            elif self.precision == 'double':
                sformat = '>6id'
                nbytes = 6 * 4 + 8
            else:
                raise RuntimeError('invalid precision format')
            data = cogsg_file.read(nbytes)

            (inew, ne, npoints, nb, npv, nev, tc) = unpack(sformat, data)
            self.header = {
                'dummy'    : dummy_int,
                'inew'     : inew, # dummy int
                'nElements': ne,  # nc;  number of tets
                'nPoints'  : npoints,  # npo; number of grid points including nbn
                'nBoundPts': nb,  # nbn; number of boundary points including nbc
                'nViscPts' : npv, # npv; number of viscous points (=0 for Euler)
                'nViscElem': nev, # ncv; number of viscous cells (=0 for Euler)
                'tc'       : tc,  # dummy double
                                  # nbc
            }
            if stop_after_header:
                return self.header
            self.log.info(str(self.header))

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
                str_format = '>' + 'i' * data_length
                data = cogsg_file.read(4 * data_length)

                faces = unpack(str_format, data)
                faces = np.array(faces)
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
                    str_format = '>%sd' % nnodes
                    unused_node_array_format = 'float64'
                elif self.precision == 'single':
                    str_format = '>%sd' % nnodes
                    unused_node_array_format = 'float32'
                else:
                    raise RuntimeError('precision = %r' % self.precision)

                data = cogsg_file.read(3 * data_length)
                assert self.precision == 'single', self.precision
                nodes = np.frombuffer(data, '>4f').reshape(3, nnodes).T.copy()
                #nodes = np.frombuffer(data, '>4f').reshape(nnodes, 3).copy()


            cogsg_file.read(nnodes * 3 * 8)  # 3 -> xyz, 8 -> double precision ???

            #----------------------------------------------------------------------
            # elements
            # boundary layer elements
            nnodes_per_tet = 4
            ntets = nev

            if ntets:
                data_length = nnodes_per_tet * ntets
                str_format = '>' + 'i' * data_length
                data = cogsg_file.read(4 * data_length)

                tets = unpack(str_format, data)
                tets = np.array(tets)
                tets = tets.reshape((tets, 4))

            #----------------------------------------------------------------------
            # volume points
            nnodes = npv

            str_format = '>%si' % nnodes
            data = cogsg_file.read(4 * nnodes)

        nodes_vol = unpack(str_format, data)
        nodes_vol = np.array(nodes_vol)
        nodes_vol = nodes_vol.reshape((tets, 3))

    def _read_cogsg_volume(self, cogsg_file):
        # volume cells
        self.log.debug('tell volume = %s' % cogsg_file.tell())
        # surface + volume cells ???
        nelements = self.header['nElements']
        #str_format = '>%si' % nelements


        self.log.debug("fv.tell = %s" % cogsg_file.tell())
        ndata = 4 * (4 * nelements)
        data = cogsg_file.read(ndata)

        # the 4 means that we make a (nelements, 4) array?
        elements = np.frombuffer(data, dtype='>i').copy() - 1
        elements = elements.reshape((4, nelements)).T
        assert elements.shape == (nelements, 4), elements.shape

        dummy2 = cogsg_file.read(4)
        self.log.debug("dummy2 = %s %s" % (unpack('>i', dummy2), unpack('>f', dummy2)))
        dummy_int2, = unpack('>i', dummy2)

        # 32 = dummy_int2 - 4 * nelements * 4
        assert self.header['dummy'] == dummy_int2

        #-----------------------------------
        # nodes
        nnodes = self.header['nPoints']
        #str_format = '>%sd' % nnodes

        dummy3 = cogsg_file.read(4)  # nnodes * 3 * 8
        dummy3_int, = unpack('>i', dummy3)
        #assert dummy3_int == 298560
        self.log.debug("dummy3 = %i" % unpack('>i', dummy3)) #, unpack('>f', dummy3)

        data_length = 8 * nnodes
        data = cogsg_file.read(3 * data_length)
        assert self.precision == 'double', self.precision
        nodes = np.frombuffer(data, '>d').reshape(3, nnodes).T

        # the ravel creates a copy that we can then use to put in
        # a contigous order
        nodes = np.asarray(nodes.ravel(), dtype='<d').reshape(nnodes, 3)

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
          - (nodeID, rho, rhoU, rhoV, rhoW) = sline
            (e) = line

        nvars = 6
          - (nodeID, rho, rhoU, rhoV, rhoW, e) = line

        Also, Nastran-esque float formatting is sometimes used,
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
        node_id, loads = read_flo(flo_filename, n=n, node_ids=node_ids)
        return node_id, loads

def read_flo(flo_filename, n=None, node_ids=None):
    """reads a *.flo file"""
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
            rhoi = parse_float(sline1[1])
            rhoui = parse_float(sline1[2])
            rhovi = parse_float(sline1[3])
            rhowi = parse_float(sline1[4])
            ei = parse_float(sline1[5])
        else:
            nvars = 5
            rhoi = parse_float(sline1[1])
            rhoui = parse_float(sline1[2])
            rhovi = parse_float(sline1[3])
            rhowi = parse_float(sline1[4])

            sline2 = flo_file.readline().strip().split()
            ei = parse_float(sline2)

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
            node_ids_minus_1 = np.array(node_ids) - 1
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
                    rhoi = parse_float(sline1[1])
                    rhoui = parse_float(sline1[2])
                    rhovi = parse_float(sline1[3])
                    rhowi = parse_float(sline1[4])
                    ei = parse_float(sline1[5])

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
                    rhoi = parse_float(sline1[1])
                    rhoui = parse_float(sline1[2])
                    rhovi = parse_float(sline1[3])
                    rhowi = parse_float(sline1[4])
                    assert len(sline1) == 5, 'len(sline1)=%s' % len(sline1)

                    sline2 = flo_file.readline().strip().split()
                    ei = parse_float(sline2)

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
                        rhoi = parse_float(sline1[1])
                        rhoui = parse_float(sline1[2])
                        rhovi = parse_float(sline1[3])
                        rhowi = parse_float(sline1[4])
                        ei = parse_float(sline1[5])

                        node_id[ni] = sline1[0]
                        rho[ni] = rhoi
                        rhoU[ni] = rhoui
                        rhoV[ni] = rhovi
                        rhoW[ni] = rhowi
                        e[ni] = ei
                        assert len(sline1) == 6, 'len(sline1)=%s' % len(sline1)
                        ni += 1
                    else:
                        unused_line1 = flo_file.readline()
            else:  # dynamic nvars=5
                for i in range(1, nmax):
                    if i in node_ids_minus_1:
                        sline1 = flo_file.readline().strip().split()
                        rhoi = parse_float(sline1[1])
                        rhoui = parse_float(sline1[2])
                        rhovi = parse_float(sline1[3])
                        rhowi = parse_float(sline1[4])
                        assert len(sline1) == 5, 'len(sline1)=%s' % len(sline1)

                        sline2 = flo_file.readline().strip().split()
                        ei = parse_float(sline2[1])

                        node_id[ni] = sline1[0]
                        rho[ni] = rhoi
                        rhoU[ni] = rhoui
                        rhoV[ni] = rhovi
                        rhoW[ni] = rhowi
                        e[ni] = ei
                        assert len(sline2) == 1, 'len(sline2)=%s' % len(sline2)
                        ni += 1
                    else:
                        unused_line1 = flo_file.readline()
                        unused_line2 = flo_file.readline()

    assert len(rho) == ni

    # limit the minimum density (to prevent division errors)
    rho_min = 0.001
    irho_zero = np.where(rho < rho_min)[0]
    rho[irho_zero] = rho_min

    loads = OrderedDict()

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
    two_over_mach2 = 2.0 / mach ** 2
    one_over_gamma = 1.0 / gamma

    gm1 = gamma - 1

    # node_id, rhoi, rhoui, rhovi, rhowi, ei
    rhoVV = (rhoU ** 2 + rhoV ** 2 + rhoW ** 2) / rho
    if 'p' in result_names or 'Mach' in result_names or 'Cp' in result_names:
        pND = gm1 * (e - rhoVV / 2.)
        if 'p' in result_names:
            loads['p'] = pND
    if 'Mach' in result_names:
        pabs = np.abs(pND)
        Mach = np.full(n, np.nan, dtype='float32')

        ipwhere = np.where(pabs > 0.0)[0]
        if len(ipwhere):
            inner = rhoVV[ipwhere] / (gamma * pabs[ipwhere])
            inwhere = ipwhere[np.where(inner >= 0.0)[0]]
            if len(inwhere):
                Mach[inwhere] = np.sqrt(rhoVV[inwhere] / (gamma * pabs[inwhere]))
        loads['Mach'] = Mach
    if 'Cp' in result_names:
        Cp = two_over_mach2 * (pND - one_over_gamma)
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

def parse_float(svalue):
    """floats a value"""
    try:
        val = float(svalue)
    except TypeError:
        val = 0.0
    return val

def write_usm3d_volume(model, basename):
    """
    writes a *.cogsg, *.front, *.face file
    """
    cogsg_filename = basename + '.cogsg'
    #face_fileame = basename + '.face'
    #front_fileame = basename + '.front'

    write_cogsg_volume(model, cogsg_filename)
    #write_front(model, front_fileame)
    #write_face(model, face_fileame)


#def write_front(model):
    #pass

#def write_face(model):
    #pass

def write_cogsg_volume(model, cogsg_fileame):
    """
    writes a *.cogsg file
    """
    #n = 0
    self = model
    nnodes = self.nodes.shape[0]
    ntets = self.tets.shape[0]

    with open(cogsg_fileame, 'wb') as outfile:
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
        str_format = '>6id'
        model.log.debug("str_format = %r" % str_format)
        data = pack(str_format, *values_header)
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
        str_format = '>%si' % ntets
        n0 = tets[:, 0] + 1
        n1 = tets[:, 1] + 1
        n2 = tets[:, 2] + 1
        n3 = tets[:, 3] + 1
        #print("n0 = %s" % n0)
        outfile.write(pack(str_format, *n0))
        outfile.write(pack(str_format, *n1))
        outfile.write(pack(str_format, *n2))
        outfile.write(pack(str_format, *n3))
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
        str_format = '>%sd' % nnodes
        nodes = self.nodes
        n0 = nodes[:, 0]
        n1 = nodes[:, 1]
        n2 = nodes[:, 2]
        outfile.write(pack(str_format, *n0))
        outfile.write(pack(str_format, *n1))
        outfile.write(pack(str_format, *n2))

        # nodes footer
        outfile.write(block_size)


def main():  # pragma: no cover
    """test problem"""
    model = Usm3d()
    if 1:
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
        unused_node_ids, unused_loads = model.read_flo(flo_filename, node_ids=[10])


if __name__ == '__main__':  # pragma: no cover
    main()
