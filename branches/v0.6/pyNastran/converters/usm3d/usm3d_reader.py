import os
from struct import pack, unpack
from numpy import array, transpose, zeros, where
from pyNastran.utils.log import get_logger

def write_front(model):
    pass

def write_face(model):
    pass

class Usm3dReader(object):
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
        self.precision = 'double'
        self.log = get_logger(log, 'debug' if debug else 'info')

    def write_usm3d(self, basename):
        if 1:
            write_usm3d_volume(self, basename)
        else:
            asdf

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
        Patch #        BC             Family   #surf   surfIDs         Family
        #----------------------------------------------------------------------
        1              44             44             0            0        Base  -> Blunt base
        2              4              4              0            0        Bay   -> Viscous Surfaces
        3              0              0              0            0        Inflow -> Supersonic Inflow
        4              2              2              0            0        Outflow -> Supersonic Outflow
        5              3              3              0            0        Sideflow -> Characteristic Inflow/Outflow
        7              1              1              0            0        Symmetry -> Reflection plane
        """
        mapbc_file = open(mapbc_filename, 'r')
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
            patch_id, bc, family, surf, surf_ids  = sline[:5]
            mapbc[int(patch_id)] = int(bc)
        mapbc_file.close()
        return mapbc

    def read_usm3d(self, basename, dimension_flag):
        cogsg_filename = basename + '.cogsg'
        bc_filename = basename + '.bc'
        face_filename = basename + '.face'
        front_filename = basename + '.front'
        mapbc_filename = basename + '.mapbc'

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
                    n = int(n)
                    nmax = max(n, nmax)

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
            tris, bcs = self.read_bc(bc_filename)
        except IOError:
            tris = None
            bcs = None
            self.log.warning('Cannot find %r...skipping' % bc_filename)
        try:
            mapbc = self.read_mapbc(mapbc_filename)
        except IOError:
            mapbc = {}
            self.log.warning('Cannot find %r...skipping' % mapbc_filename)
        self.tris = tris
        self.bcs = bcs
        self.mapbc = mapbc

        loads = {}
        if os.path.exists(flo_filename):
            npoints, three = nodes.shape
            loads = self.read_flo(flo_filename, n=npoints)
        else:
            self.log.warning('Cannot find %r...skipping' % flo_filename)
        self.loads = loads

        return nodes, elements, tris, bcs, mapbc, loads
        #self.read_front(front_file)
        #self.read_face(face_file)

    def read_bc(self, bc_file):
        f = open(bc_file, 'r')
        lines = f.readlines()
        f.close()
        (ntris, dunnoA, dunnoB, dunnoC) = lines[0].strip().split()
        ntris = int(ntris)
        tris = zeros((ntris, 3), dtype='int32')
        bcs  = zeros(ntris, dtype='int32')

        for i in xrange(ntris):
            (n, isurf, n1, n2, n3) = lines[i+2].split()
            tris[i] = [n1, n2, n3]
            bcs[i] = isurf
        tris = tris - 1
        return tris, bcs
        #self.bcs = [tris, bcs]


    def read_cogsg(self, cogsg_file):
        f = open(cogsg_file, 'rb')

        # nelements * 4 * 4 + 32 ???
        dummy = f.read(4)  # 1022848
        dummy_int, = unpack('>i', dummy)
        #assert dummy_int == 1022848, 'dummy_int = %s' % dummy_int

        # file header
        if self.precision=='single':
            Format = '>6if'
            nbytes = 6 * 4 + 4
            self.log.debug("***single precision")
        elif self.precision=='double':
            Format = '>6id'
            nbytes = 6 * 4 + 8
        else:
            raise Exception('invalid precision format')
        data = f.read(nbytes)

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
        self.log.info(self.header)

        # nbn nodes
        #
        #del ne, np

        if 1:
            return self._read_cogsg_volume(f)
        #else:
        #----------------------------------------------------------------------
        # elements
        # face elements
        nnodes_per_face = 3
        nfaces = ne

        if nfaces > 0:
            data_length = nnodes_per_face * nfaces
            Format = '>' + 'i' * data_length
            data = f.read(4 * data_length)

            faces = unpack(Format, data)
            faces = array(faces)
            faces = faces.reshape((nfaces, 3))
        else:
            faces = None
        #----------------------------------------------------------------------
        # nodes
        nBoundPts = nb
        nnodes = nBoundPts

        #data_length = nnodes
        if self.precision == 'double':
            data_length = 8 * nnodes
        elif self.precision == 'single':
            data_length = 4 * nnodes
        else:
            raise RuntimeError('precision = %r' % self.precision)

        skip_nodes = False
        if skip_nodes == True:
            t = self.tell()
            f.goto(t + data_length * 3)
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

            data = f.read(data_length)
            X = unpack(Format, data)

            data = f.read(data_length)
            Y = unpack(Format, data)

            data = f.read(data_length)
            Z = unpack(Format, data)

            nodes = array([X, Y, Z])

            nodes = zeros((nnodes, 3), node_array_format)
            nodes[:, 0] = X
            nodes[:, 1] = Y
            nodes[:, 2] = Z
            del X, Y, Z


        f.read(nnodes * 3 * 8)  # 3 -> xyz, 8 -> double precision ???

        #----------------------------------------------------------------------
        # elements
        # boundary layer elements
        nnodes_per_tet = 4
        ntets = nev

        if ntets:
            data_length = nnodes_per_tet * ntets
            Format = '>' + 'i' * data_length
            data = f.read(4 * data_length)

            tets = unpack(Format, data)
            tets = array(tets)
            tets = tets.reshape((tets, 4))

        #----------------------------------------------------------------------
        # volume points
        nnodes = npv

        Format = '>%si' % nnodes
        data = f.read(4 * nnodes)

        nodes_vol = unpack(Format, data)
        nodes_vol = array(nodes_vol)
        nodes_vol = nodes_vol.reshape((tets, 3))



    def _read_cogsg_volume(self, f):
        # volume cells
        self.log.debug('tell volume = %s' % f.tell())
        # surface + volume cells ???
        nelements = self.header['nElements']
        Format = '>%si' % nelements

        elements = zeros((nelements, 4), 'int32')

        self.log.debug("fv.tell = %s" % f.tell())
        for i in range(4): #  tets
            data = f.read(4 * nelements)
            elements[:, i] = unpack(Format, data)
            #print "elements[:, %s] =" %i, elements[:, i]
        elements -= 1

        #print 'tell volume2 =', f.tell()
        dummy2 = f.read(4)
        self.log.debug("dummy2 = %s %s" % (unpack('>i', dummy2), unpack('>f', dummy2)))
        dummy_int2, = unpack('>i', dummy2)

        # 32 = dummy_int2 - 4 * nelements * 4
        assert self.header['dummy'] == dummy_int2

        #-----------------------------------
        # nodes
        nnodes = self.header['nPoints']
        Format = '>%sd' % nnodes

        dummy3 = f.read(4)  # nnodes * 3 * 8
        dummy3_int, = unpack('>i', dummy3)
        #assert dummy3_int == 298560
        print "dummy3 = ", unpack('>i', dummy3), unpack('>f', dummy3)

        nodes = zeros((nnodes, 3), 'float64')
        for i in range(3): #  x, y, z
            data = f.read(8 * nnodes)
            assert len(data) == (8 * nnodes)
            nodes[:, i] = unpack(Format, data)

        dummy4 = f.read(4) # nnodes * 3 * 8
        dummy4_int, = unpack('>i', dummy4)
        #print "dummy4 = ", unpack('>i', dummy4), unpack('>f', dummy4)

        assert dummy3_int == dummy4_int
        self.nodes = nodes
        self.tets = elements
        return nodes, elements

    def read_flo(self, flo_filename, n=None):
        """
        ipltqn is a format code where:
         - ipltqn = 0  (no printout)
         - ipltqn = 1  (unformatted)
         - ipltqn = 2  (formatted) - default

        :param flo_filename: the name of the file to read
        :param n:            the number of points to read (initializes the array)

        nvars = 5
          - (nodeID,rho,rhoU,rhoV,rhoW) = sline
            (e) = line

        nvars = 6
          - (nodeID,rho,rhoU,rhoV,rhoW,e) = line
        """
        formatCode = 2

        flo_file = open(flo_filename, 'r')
        mach = float(flo_file.readline().strip())

        node_id = zeros(n, 'int32')
        rho = zeros(n, 'float32')
        rhoU = zeros(n, 'float32')
        rhoV = zeros(n, 'float32')
        rhoW = zeros(n, 'float32')
        e = zeros(n, 'float32')

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
        i = 0
        node_id[i] = sline1[0]
        rho[i] = rhoi
        rhoU[i] = rhoui
        rhoV[i] = rhovi
        rhoW[i] = rhowi
        e[i] = ei

        # loop over the rest of the data in the flo file
        if n:
            if nvars == 6:
                for i in xrange(1, n):
                    sline1 = flo_file.readline().strip().split()
                    rhoi, rhoui, rhovi, rhowi, ei = Float(sline1[1:], 5)
                    node_id[i] = sline1[0]
                    rho[i] = rhoi
                    rhoU[i] = rhoui
                    rhoV[i] = rhovi
                    rhoW[i] = rhowi
                    e[i] = ei
                    assert len(sline1) == 6, 'len(sline1)=%s' % len(sline1)
            else:  # nvars = 5
                for i in xrange(1, n):
                    sline1 = flo_file.readline().strip().split()
                    node_id[i] = sline1[0]
                    rhoi, rhoui, rhovi, rhowi = Float(sline1[1:], 4)

                    assert len(sline1) == 5, 'len(sline1)=%s' % len(sline1)

                    sline2 = flo_file.readline().strip().split()
                    ei = Float(sline2, 1)[0]

                    rho[i] = rhoi
                    rhoU[i] = rhoui
                    rhoV[i] = rhovi
                    rhoW[i] = rhowi
                    e[i] = ei
                    assert len(sline2) == 1, 'len(sline2)=%s' % len(sline2)
        else:
            raise RuntimeError('nnodes is not defined...')
        flo_file.close()

        # llimit the minimum density (to prevent division errors)
        rho_min = 0.001
        irho_zero = where(rho < rho_min)[0]
        rho[irho_zero] = rho_min

        result_names = ['Mach', 'U', 'V', 'W', 'T', 'rhoU', 'p', 'Cp']
        loads = {}

        gamma = 1.4
        two_over_Mach2 = 2.0 / mach ** 2
        one_over_gamma = 1.0 / gamma

        gm1 = gamma - 1
        rhoVV = (rhoU**2 + rhoV**2 + rhoW**2) / rho
        if 'p' in result_names or 'Mach' in result_names or 'Cp' in result_names:
            pND = gm1*(e - rhoVV/2. )
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
        return loads




def Float(sline, n):
    vals = []
    for val in sline:
        try:
            vals.append(float(val))
        except:
            vals.append(0.0)
    return vals

def write_usm3d_volume(model, basename):
    cogsg_file = basename + '.cogsg'
    face_file = basename + '.face'
    front_file = basename + '.front'

    write_cogsg_volume(model, cogsg_file)
    #write_front(model, front_file)
    #write_face(model, face_file)

def write_cogsg_volume(model, cogsg_file):
    #n = 0
    self = model
    nnodes, three = self.nodes.shape
    ntets, four = self.tets.shape

    outfile = open(cogsg_file, 'wb')

    # file header

    values = [32 + ntets * 4 * 4,]
    block_size = pack('>i', *values)
    outfile.write(block_size)

    header = self.header
    values_header = [header['inew'],
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
    #print "outfile.tell = ", outfile.tell()
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
    #print "n0 =", n0
    outfile.write( pack(Format, *n0) )
    outfile.write( pack(Format, *n1) )
    outfile.write( pack(Format, *n2) )
    outfile.write( pack(Format, *n3) )
    #n += 4 * 8 * ntets
    #print "outfile.tell 2 = ", outfile.tell(), n

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
    outfile.write( pack(Format, *n0) )
    outfile.write( pack(Format, *n1) )
    outfile.write( pack(Format, *n2) )

    # nodes footer
    outfile.write(block_size)

    #--------------------------------------------------------------------------
    outfile.close()


def main():
    model = Usm3dReader()
    #basename = 'HSCT_inviscid'
    basename = 'box'
    model.read_usm3d(basename, 3)
    model.write_usm3d(basename + '_2')

    model.read_usm3d(basename + '_2', 3)


if __name__ == '__main__':
    main()