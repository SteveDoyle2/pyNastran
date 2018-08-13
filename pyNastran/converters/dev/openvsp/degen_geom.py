from __future__ import print_function
from copy import deepcopy
from collections import defaultdict
from six import iteritems

import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.converters.panair.panair_grid import PanairGrid, PanairPatch


class Geom(object):
    """defines a lifting surfaces"""
    def __init__(self, name, lifting_surface_xyz,
                 lifting_surface_nx, lifting_surface_ny):
        self.name = name
        self.xyz = lifting_surface_xyz
        self.nx = lifting_surface_nx
        self.ny = lifting_surface_ny

    def write_bdf_file_obj(self, bdf_file, nid0=1, eid=1, pid=1):
        """exports a GEOM to Nastran formatted file"""
        nx = self.nx
        ny = self.ny
        nxy = nx * ny
        cp = None

        for ni in range(nxy):
            x, y, z = self.xyz[ni, :]
            card = ['GRID', nid0 + ni, cp, x, y, z]
            bdf_file.write(print_card_8(card))
        ni = nxy - 1

        #print('ni', ni, self.xyz.shape)
        eidi = 0
        for i in range(nx - 1):
            for j in range(ny - 1):
                g1 = nid0 + i*ny + j
                g2 = nid0 + i*ny + j + 1
                g3 = nid0 + (i+1) * ny + j + 1
                g4 = nid0 + (i+1) * ny + j
                card = ['CQUAD4', eid + eidi, pid, g1, g2, g3, g4]
                bdf_file.write(print_card_8(card))
                eidi += 1
        nid0 += ni + 1
        eid += eidi + 1
        return nid0, eid, pid

    @property
    def elements(self):
        nid0 = 1
        #eidi = 0
        k = 0

        nx = self.nx
        ny = self.ny
        nxy = nx * ny
        elements = np.zeros((nxy, 4), dtype='int32')
        for i in range(nx - 1):
            for j in range(ny - 1):
                g1 = nid0 + i*ny + j
                g2 = nid0 + i*ny + j + 1
                g3 = nid0 + (i+1) * ny + j + 1
                g4 = nid0 + (i+1) * ny + j
                card = [g1, g2, g3, g4]
                elements[k, :] = card
                k += 1
        return elements

    def __repr__(self):
        msg = ('Geom(name=%s, lifting_surface_xyz, '
               'lifting_surface_nx, lifting_surface_ny)' % (self.name))
        return msg


class DegenGeom(object):
    def __init__(self, log=None, debug=False):
        self.log = log
        self.debug = debug
        self.components = defaultdict(list)

    def write_bdf(self, bdf_filename):
        bdf_file = open(bdf_filename, 'wb')
        bdf_file.write('$pyNastran: VERSION=NX\n')
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')

        nid = 1
        eid = 1
        pid = 1

        mid = 1
        t = 0.1
        E = 3.0e7
        G = None
        nu = 0.3
        card = ['MAT1', mid, E, G, nu]
        bdf_file.write(print_card_8(card))
        for name, comps in sorted(iteritems(self.components)):
            bdf_file.write('$ name = %r\n' % name)
            for comp in comps:
                card = ['PSHELL', pid, mid, t]
                bdf_file.write(print_card_8(card))
                nid, eid, pid = comp.write_bdf_file_obj(bdf_file, nid, eid, pid)
                pid += 1

    def write_panair(self, panair_filename, panair_case_filename):  # pragma: no cover
        pan = PanairGrid()
        pan.mach = 0.5
        pan.is_end = True
        pan.ncases = 2
        pan.alphas = [0., 5.]


        i = 0
        pan.nNetworks = 1
        kt = 1
        cp_norm = 1
        for name, comps in sorted(iteritems(self.components)):
            #panair_file.write('$ name = %r\n' % name)
            for comp in comps:
                namei = name + str(i)
                x = deepcopy(comp.lifting_surface_xyz[:, 0])
                y = deepcopy(comp.lifting_surface_xyz[:, 1])
                z = deepcopy(comp.lifting_surface_xyz[:, 2])
                x = x.reshape((comp.lifting_surface_nx, comp.lifting_surface_ny))
                y = y.reshape((comp.lifting_surface_nx, comp.lifting_surface_ny))
                z = z.reshape((comp.lifting_surface_nx, comp.lifting_surface_ny))

                xyz = np.dstack([x, y, z])
                assert xyz.shape[2] == 3
                patch = PanairPatch(pan.nNetworks, namei, kt, cp_norm, xyz, self.log)
                pan.patches[i] = patch
                pan.nNetworks += 1
                i += 1

                if 'wing' in name.lower():  # make a wing cap
                    namei = 'cap%i' % i
                    #assert comp.lifting_surface_nx == 6, comp.lifting_surface_nx
                    assert comp.lifting_surface_ny == 33, comp.lifting_surface_ny
                    #print(x.shape)
                    xend = deepcopy(x[-1, :])
                    print(xend)
                    yend = deepcopy(y[-1, :])
                    zend = deepcopy(z[-1, :])
                    imid = comp.lifting_surface_ny // 2
                    x = np.zeros((imid+1, 2), dtype='float32')
                    y = np.zeros((imid+1, 2), dtype='float32')
                    z = np.zeros((imid+1, 2), dtype='float32')
                    print(imid, xend[imid], xend.min())
                    xflip = list(xend[0:imid+1])
                    yflip = list(yend[0:imid+1])
                    zflip = list(zend[0:imid+1])
                    x[:, 0] = xflip[::-1]
                    y[:, 0] = yflip[::-1]
                    z[:, 0] = zflip[::-1]
                    x[:, 1] = xend[imid:]
                    y[:, 1] = yend[imid:]
                    z[:, 1] = zend[imid:]
                    print(x)

                    #x = xend[0:imid:-1].extend(x[imid:])
                    #y = yend[0:imid:-1].extend(y[imid:])
                    #z = zend[0:imid:-1].extend(z[imid:])
                    #print(x)
                    x = x.reshape((2, imid+1))
                    y = y.reshape((2, imid+1))
                    z = z.reshape((2, imid+1))

                    #print(xend)
                    xyz = np.dstack([x, y, z])
                    assert xyz.shape[2] == 3
                    patch = PanairPatch(pan.nNetworks, namei, kt, cp_norm,
                                        xyz, self.log)
                    pan.patches[i] = patch
                    pan.nNetworks += 1
                    i += 1
                #i += 1
        pan.write_panair(panair_filename)
        #self.nNetworks = i

    def read_degen_geom(self, degen_geom_csv):
        with open(degen_geom_csv, 'r') as degen_geom_file:
            lines = degen_geom_file.readlines()

        degen_geom_file = open(degen_geom_csv, 'r')
        for i in range(4):
            line = degen_geom_file.readline()

        iline = 3
        line = lines[iline]

        ncomponents = int(line)
        self.log.info('ncomponents=%s' % ncomponents)
        degen_geom_file.readline()
        iline += 1

        for icomp in range(ncomponents):
            sline = degen_geom_file.readline().strip().split(',')
            iline += 1

            sline2 = lines[iline].strip().split(',')
            assert sline == sline2, 'iline=%s \nsline1=%s \nsline2=%s' % (iline, sline, sline2)
            surface_type = sline[0]
            self.log.info('surface_type=%r' % surface_type)

            if surface_type == 'LIFTING_SURFACE':
                lifting_surface, name, dunno_zero = sline
                if int(dunno_zero) not in [0, 1]:
                    raise RuntimeError('unexpected line %s; should be 0 or 1 '
                                       'for the 3rd value:\n%s' % (iline, sline))

                self.log.info('lifting_surface name=%r' % name)
                (iline, nelements, normals, area, lifting_surface_xyz,
                 lifting_surface_nx, lifting_surface_ny) = read_surface_node(
                     degen_geom_file, lines, iline, self.log)

                self.log.debug('iline=%s' % iline)

                iline = read_surface_face(
                    degen_geom_file, lines, iline, nelements, self.log,
                    normals, area)

                iline = read_plate(degen_geom_file, lines, iline, self.log)
                iline = read_stick_node(degen_geom_file, lines, iline, self.log)
                iline = read_stick_face(degen_geom_file, lines, iline, self.log)

                # DegenGeom Type
                # POINT
                #(vol, volWet, area, areaWet,
                 #Ishellxx, Ishellyy, Ishellzz, Ishellxy, Ishellxz, Ishellyz,
                 #Isolidxx, Isolidyy, Isolidzz, Isolidxy, Isolidxz, Isolidyz,
                 #cgShellx, cgShelly, cgShellz, cgSolidx, cgSolidy, cgSolidz)
                degen_geom_file.readline()
                iline += 1

                point = degen_geom_file.readline().strip()
                iline += 1
                assert point == 'POINT', point

                degen_geom_file.readline()
                iline += 1

                cg_line = degen_geom_file.readline().split()
                assert len(cg_line) == 22, len(cg_line)
                sline2 = lines[iline].strip().split(',')
                assert len(sline2) == 22, len(sline2)
                iline += 1

                degen_geom_file.readline()
                iline += 1
                component = Geom(name, lifting_surface_xyz,
                                 lifting_surface_nx, lifting_surface_ny)
            elif surface_type == 'BODY':
                lifting_surface, name = sline

                # DegenGeom Type,nXsecs, nPnts/Xsec
                #line = degen_geom_file.readline()
                #iline += 1

                #SURFACE_NODE,21,17
                #line = degen_geom_file.readline()
                #iline += 1
                #print(line)

                (iline, nelements, normals, area, lifting_surface_xyz,
                 lifting_surface_nx, lifting_surface_ny) = read_surface_node(
                     degen_geom_file, lines, iline, self.log)

                iline = read_surface_face(
                    degen_geom_file, lines, iline, nelements, self.log,
                    normals, area)
                iline = read_plate(degen_geom_file, lines, iline, self.log)
                iline = read_plate(degen_geom_file, lines, iline, self.log)

                iline = read_stick_node(degen_geom_file, lines, iline, self.log)
                iline = read_stick_face(degen_geom_file, lines, iline, self.log)

                #iline = read_stick_node(degen_geom_file, lines, iline, self.log)
                #iline = read_stick_face(degen_geom_file, lines, iline, self.log)

                line = degen_geom_file.readline()
                iline += 1
                #print('**', line)
                component = Geom(name, lifting_surface_xyz,
                                 lifting_surface_nx, lifting_surface_ny)

            else:
                raise RuntimeError(sline)
            self.components[name].append(component)

def read_surface_node(degen_geom_file, lines, iline, log):
    """
    degenGeom, Type, nxsections, npoints/xsection
    SURFACE_NODE,6,33
    # nnodes -> 6*33=198
    # nelements -> 160
    """
    line = degen_geom_file.readline()
    #print(line)
    iline += 1

    line = degen_geom_file.readline().strip()
    iline += 1

    sline = line.split(',')
    log.info(str(sline))
    surface_node, lifting_surface_nx, lifting_surface_ny = sline
    assert surface_node == 'SURFACE_NODE', surface_node
    lifting_surface_nx = int(lifting_surface_nx)
    lifting_surface_ny = int(lifting_surface_ny)
    npoints = lifting_surface_nx * lifting_surface_ny
    nelements = (lifting_surface_nx - 1) * (lifting_surface_ny - 1)
    log.info('npoints = %r' % npoints)
    lifting_surface_xyz = np.zeros((npoints, 3), dtype='float64')
    normals = np.zeros((nelements, 3), dtype='float64')
    area = np.zeros(nelements, dtype='float64')

    # x, y, z, u, v
    degen_geom_file.readline()
    iline += 1
    for ipoint in range(npoints):
        sline = degen_geom_file.readline().strip().split(',')
        iline += 1
        x, y, z, u, v = sline
        sline2 = lines[iline].strip().split(',')
        assert sline == sline2, 'iline=%s \nsline1=%s \nsline2=%s' % (iline, sline, sline2)
        #log.debug('%s: %s' % (iline, sline2))
        lifting_surface_xyz[ipoint, :] = [x, y, z]
    iline += 1
    return iline, nelements, normals, area, lifting_surface_xyz, lifting_surface_nx, lifting_surface_ny

def read_surface_face(degen_geom_file, lines, iline, nelements, log,
                      normals, area):
    """SURFACE_FACE,5,32"""
    sline = degen_geom_file.readline().strip().split(',')
    surface_node, nx, ny = sline
    assert surface_node == 'SURFACE_FACE', 'surface_node=%s, nx+%s, ny=%s' % (surface_node, nx, ny)
    sline2 = lines[iline].strip().split(',')
    assert sline == sline2, 'iline=%s \nsline1=%s \nsline2=%s' % (iline, sline, sline2)

    line = degen_geom_file.readline() # nx,ny,nz,area
    iline += 1
    log.debug(line)
    for ielem in range(nelements):
        line = degen_geom_file.readline()
        iline += 1
        nx, ny, nz, areai = line.split(',')
        normals[ielem, :] = [nx, ny, nz]
        area[ielem] = areai
    iline += 1
    return iline

def read_plate(degen_geom_file, lines, iline, log):
    """
    DegenGeom Type,nXsecs,nPnts/Xsec
    PLATE,6,17
    nx,ny,nz

    TODO: the plate is very unclear...it's 6 lines with 3 normals on each line
          but 6*17?
    """
    degen_geom_file.readline()
    iline += 1

    sline = degen_geom_file.readline().strip().split(',')
    sline2 = lines[iline].strip().split(',')
    assert sline == sline2, 'iline=%s \nsline1=%s \nsline2=%s' % (iline, sline, sline2)
    plate, nx, ny = sline

    iline += 1
    nx = int(nx)
    ny = int(ny)
    nxy = nx * ny
    degen_geom_file.readline()
    iline += 1
    for i in range(nx):
        degen_geom_file.readline()
        iline += 1

    # x,y,z,zCamber,t,nCamberx,nCambery,nCamberz,u,wTop,wBot
    line = degen_geom_file.readline()
    iline += 1
    for i in range(nxy):
        degen_geom_file.readline()
        iline += 1
    return iline

def read_stick_node(degen_geom_file, lines, iline, log):
    """
    DegenGeom Type, nXsecs
    STICK_NODE, 6
    (lex, ley, lez, tex, tey, tez, cgShellx, cgShelly, cgShellz,
     cgSolidx, cgSolidy, cgSolidz, toc, tLoc, chord,
     Ishell11, Ishell22, Ishell12, Isolid11, Isolid22, Isolid12,
     sectArea, sectNormalx, sectNormaly, sectNormalz,
     perimTop, perimBot, u,
     t00, t01, t02, t03, t10, t11, t12, t13,
     t20, t21, t22, t23, t30, t31, t32, t33,
     it00, it01, it02, it03, it10, it11, it12, it13,
     it20, it21, it22, it23, it30, it31, it32, it33)
    """
    degen_geom_file.readline()
    iline += 1
    stick_node, nx = degen_geom_file.readline().split(',')
    iline += 1
    assert stick_node == 'STICK_NODE', stick_node
    nx = int(nx)
    degen_geom_file.readline()
    iline += 1
    for i in range(nx):
        degen_geom_file.readline()
        iline += 1
    return iline

def read_stick_face(degen_geom_file, lines, iline, log):
    """
    DegenGeom Type, nXsecs
    STICK_FACE, 5
    sweeple,sweepte,areaTop,areaBot
    """
    degen_geom_file.readline()
    iline += 1

    sline = degen_geom_file.readline().strip().split(',')
    sline2 = lines[iline].strip().split(',')
    assert sline == sline2, 'iline=%s \nsline1=%s \nsline2=%s' % (iline, sline, sline2)
    stick_face, nx = sline
    iline += 1

    assert stick_face == 'STICK_FACE', stick_face
    nx = int(nx)
    degen_geom_file.readline()
    for i in range(nx):
        degen_geom_file.readline()
        iline += 1
    return iline

def main():  # pragma: no cover
    degen_geom_csv = 'model_DegenGeom.csv'
    d = DegenGeom()
    d.read_degen_geom(degen_geom_csv)
    d.write_bdf('model.bdf')

    panair_filename = 'panair.inp'
    panair_case_filename = 'model_DegenGeom.vspaero'
    d.write_panair(panair_filename, panair_case_filename)


if __name__ == '__main__':  # pragma: no cover
    main()
