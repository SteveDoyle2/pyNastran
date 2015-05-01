from six import iteritems
from collections import defaultdict
from numpy import zeros
from pyNastran.bdf.field_writer_8 import print_card_8

class Geom(object):
    def __init__(self, name, lifting_surface_xyz, lifting_surface_nx, lifting_surface_ny):
        self.name = name
        self.lifting_surface_xyz = lifting_surface_xyz
        self.lifting_surface_nx = lifting_surface_nx
        self.lifting_surface_ny = lifting_surface_ny

    def write_bdf_file_obj(self, bdf_file, nid=1, eid=1, pid=1):
        nx = self.lifting_surface_nx
        ny = self.lifting_surface_ny

        cp = None
        for ni, (x, y, z) in enumerate(self.lifting_surface_xyz):
            card = ['GRID', nid + ni, cp, x, y, z]
            bdf_file.write(print_card_8(card))

        eidi = 0
        for j in range(ny - 1):
            for i in range(nx - 1):
                # this is wrong...
                g1 = nid + i
                g2 = nid + i + 1
                g3 = nid + nx + i + 1
                g4 = nid + nx + i
                card = ['CQUAD4', eid + eidi, pid, g1, g2, g3, g4]
                bdf_file.write(print_card_8(card))
                eidi += 1
            nid += nx
        print(nx, ny)
        print(eidi)
        nid += ni + 1
        eid += eidi + 1
        return nid, eid, pid

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
            card = ['PSHELL', pid, mid, t]
            bdf_file.write(print_card_8(card))
            for comp in comps:
                nid, eid, pid = comp.write_bdf_file_obj(bdf_file, nid, eid, pid)

    def read_degen_geom(self, degen_geom_csv):
        f = open(degen_geom_csv)
        for i in range(4):
            line = f.readline()
        ncomponents = int(line)
        f.readline()

        for i in range(ncomponents):
            line = f.readline().strip().split(',')
            lifting_surface, name = line
            if lifting_surface == 'LIFTING_SURFACE':
                #degenGeom, Type, nxsections, npoints/xsection
                # SURFACE_NODE,6,33
                # nnodes -> 6*33=198
                # nelements -> 160
                f.readline()
                surface_node, lifting_surface_nx, lifting_surface_ny = f.readline().strip().split(',')
                assert surface_node == 'SURFACE_NODE', surface_node
                lifting_surface_nx = int(lifting_surface_nx)
                lifting_surface_ny = int(lifting_surface_ny)
                npoints = lifting_surface_nx * lifting_surface_ny
                nelements = (lifting_surface_nx - 1) * (lifting_surface_ny - 1)
                print('npoints = %r' % npoints)
                lifting_surface_xyz = zeros((npoints, 3), dtype='float64')
                normals = zeros((nelements, 3), dtype='float64')
                area = zeros(nelements, dtype='float64')

                # x, y, z, u, v
                f.readline()
                for i in range(npoints):
                    line = f.readline()
                    x, y, z, u, v = line.split(',')
                    lifting_surface_xyz[i, :] = [x, y, z]

                # SURFACE_FACE,5,32
                surface_node, nx, ny = f.readline().strip().split(',')
                line = f.readline() # nx,ny,nz,area
                print(line)
                for i in range(nelements):
                    line = f.readline()
                    nx, ny, nz, areai = line.split(',')
                    normals[i, :] = [nx, ny, nz]
                    area[i] = areai


                # TODO: the plate is very unclear...it's 6 lines with 3 normals on each line
                #       but 6*17?
                # DegenGeom Type,nXsecs,nPnts/Xsec
                # PLATE,6,17
                # nx,ny,nz
                f.readline()
                plate, nx, ny = f.readline().strip().split(',')
                nx = int(nx)
                ny = int(ny)
                nxy = nx * ny
                f.readline()
                for i in range(nx):
                    f.readline()

                # x,y,z,zCamber,t,nCamberx,nCambery,nCamberz,u,wTop,wBot
                line = f.readline()
                for i in range(nxy):
                    f.readline()

                # DegenGeom Type, nXsecs
                # STICK_NODE, 6
                # lex,ley,lez,tex,tey,tez,cgShellx,cgShelly,cgShellz,cgSolidx,cgSolidy,cgSolidz,toc,tLoc,chord,Ishell11,Ishell22,Ishell12,Isolid11,Isolid22,Isolid12,sectArea,sectNormalx,sectNormaly,sectNormalz,perimTop,perimBot,u,t00,t01,t02,t03,t10,t11,t12,t13,t20,t21,t22,t23,t30,t31,t32,t33,it00,it01,it02,it03,it10,it11,it12,it13,it20,it21,it22,it23,it30,it31,it32,it33,
                f.readline()
                stick_node, nx = f.readline().split(',')
                assert stick_node == 'STICK_NODE', stick_node
                nx = int(nx)
                f.readline()
                for i in range(nx):
                    f.readline()

                # DegenGeom Type, nXsecs
                # STICK_FACE, 5
                # sweeple,sweepte,areaTop,areaBot
                f.readline()
                stick_face, nx = f.readline().split(',')
                assert stick_face == 'STICK_FACE', stick_face
                nx = int(nx)
                f.readline()
                for i in range(nx):
                    f.readline()

                # DegenGeom Type
                # POINT
                # vol,volWet,area,areaWet,Ishellxx,Ishellyy,Ishellzz,Ishellxy,Ishellxz,Ishellyz,Isolidxx,Isolidyy,Isolidzz,Isolidxy,Isolidxz,Isolidyz,cgShellx,cgShelly,cgShellz,cgSolidx,cgSolidy,cgSolidz
                f.readline()
                point = f.readline().strip()
                assert point == 'POINT', point
                f.readline()
                cg_line = f.readline()
                f.readline()

            else:
                raise RuntimeError(line)
            component = Geom(name, lifting_surface_xyz, lifting_surface_nx, lifting_surface_ny)
            self.components[name].append(component)



def main():
    degen_geom_csv = 'model_DegenGeom.csv'
    d = DegenGeom()
    d.read_degen_geom(degen_geom_csv)
    d.write_bdf('model.bdf')


if __name__ == '__main__':
    main()
