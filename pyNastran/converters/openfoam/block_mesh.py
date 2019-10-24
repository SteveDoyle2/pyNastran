"""
defines:
 - model = read_block_mesh(block_mesh_dict_filename, log=None, debug=False)

"""
from collections import defaultdict
from itertools import count

import numpy as np
from numpy.linalg import norm  # type: ignore

from pyNastran.converters.openfoam.openfoam_parser import FoamFile, convert_to_dict #, write_dict
from pyNastran.bdf.field_writer_8 import print_card_8
from cpylog import get_logger2


def read_block_mesh(block_mesh_dict_filename, log=None, debug=False):
    """functional interface to BlockMesh"""
    model = BlockMesh(log=log, debug=debug)
    model.read_openfoam(block_mesh_dict_filename)
    return model

def area_centroid(n1, n2, n3, n4):
    """calculates the area and centroid of a quad"""
    centroid = (n1 + n2 + n3 + n4) / 4.
    n = np.cross(n3 - n1, n4 - n2)
    area = 0.5 * norm(n)
    return area, centroid

def volume8(n1, n2, n3, n4, n5, n6, n7, n8):
    """calculates the volume of a hex"""
    (A1, c1) = area_centroid(n1, n2, n3, n4)
    (A2, c2) = area_centroid(n5, n6, n7, n8)
    V = (A1 + A2) / 2. * norm(c1 - c2)
    return V

def get_not_indexes(a, indices):
    """only works for 1D"""
    ia = np.indices(a.shape)
    not_indices = np.setxor1d(ia, indices)
    return not_indices


class BlockMesh:
    """defines BlockMesh"""
    def __init__(self, log=None, debug=False):
        """creates BlockMesh"""
        debug = False
        #log = None
        self.log = get_logger2(log, debug=debug)

        # arrays
        self.nodes = None
        self.quads = None
        self.hexas = None
        self.npoints = None
        self.grading = None

        self.iname_to_quads = None
        self.inames = None
        self.bcs = None
        self.iname_to_name = None
        self.iname_to_type = None

    def make_hex_bar(self, unused_bias, ncells):
        """
        bias = Llast/Lfirst
        """
        #k = bias = ** 1/N  # close to this-ish
        k = 1.
        assert isinstance(ncells, int), ncells
        ipoints = np.arange(ncells + 1) # ipoint
        # xmax = d0 * k**n
        kn = k**ipoints
        L0 = 1.0 / sum(kn)
        x = kn * L0
        return ipoints, x

    def make_hex_bars(self, unused_hex_id):
        #self.hexas = hexas
        #self.npoints = npoints
        #self.grading = grading
        points = []
        line_pairs = [
            # a set is 2 values
            #set1  set2 idir,face
            (0, 1, 4, 5, 0, 2),
            (3, 2, 7, 6, 0, 2),
            (0, 3, 1, 2, 1, 1),
            (4, 7, 5, 6, 1, 1),
            (0, 4, 1, 5, 2, 0),
            (3, 7, 2, 6, 2, 0),
        ]

        points = []
        for grading, unused_hexa in zip(self.grading, self.hexas):
            # x = 0-1 to 4-5
            #     3-2 to 7-6
            # y = 0-3 to 1-2
            #     4-7 to 5-6
            # z = 0-4 to 1-5
            #     3-7 to 2-6
            nx, ny, nz, unused_method = grading
            bias = (nx, ny, nz)
            for line_pair in line_pairs:
                i1, i2, i3, i4, idir, iface = line_pair
                n1 = self.nodes[i1, :]
                n2 = self.nodes[i2, :]
                n3 = self.nodes[i3, :]
                n4 = self.nodes[i4, :]
                ncells_x = grading[idir]
                ncells_y = grading[iface]
                bias_x = bias[idir]
                bias_y = bias[iface]

                unused_npx, x = self.make_hex_bar(bias_x, ncells_x)
                unused_npy, y = self.make_hex_bar(bias_y, ncells_y)

                da = n1 - n2
                db = n3 - n4
                La = norm(da)
                Lb = norm(db)
                xa = La * x
                xb = Lb * x

                ncells = min(ncells_x, ncells_y)  ## TODO: what should this be?

                for i in range(1, ncells):
                    #print(n2)
                    #print(xa[i])
                    p1 = n2 + xa[i]
                    p2 = n4 + xb[i]
                    dp = p2 - p1
                    L = norm(dp)
                    unused_yout = L * y

                    #unused_new_points = p1 + dp * y
                    points.append(1)
                    i / L
                    #p1 =

    def read_openfoam(self, block_mesh_name='blockMeshDict'):
        """reads the BlockMesh file"""
        self.log.info('block_mesh_name = %r' % block_mesh_name)
        foam_file = FoamFile(block_mesh_name)
        foam_lines = foam_file.read_foam_file()

        foam_file_dict = convert_to_dict(self, foam_lines, debug=True)
        #print(write_dict(foam_file_dict))
        #unused_keys = foam_file_dict.keys()

        vertices = foam_file_dict['vertices']
        blocks = foam_file_dict['blocks']
        boundaries = foam_file_dict['boundary']
        #self.log.info(boundaries)

        nodes = []
        for unused_ivertex, vertex in vertices.items():
            x, y, z = vertex.strip('() ').split()
            nodes.append([x, y, z])

        hexas = []
        npoints = []
        grading = []
        for key, block in blocks.items():
            hexa, npointsi, gradingi, blank = block.split(')')
            assert blank == '', '%r' % blank

            hexa = hexa.replace('hex (', '').split()
            npointsi = npointsi.strip('( ').split()
            gradingi = gradingi.split('(')[1].split()
            #self.log.info('%s %s' % (npointsi, gradingi))
            #self.log.info(hexa)
            hexa = [int(i) for i in hexa]
            #print('hexa', key, hexa)
            hexas.append(hexa)
            npoints.append(npointsi)
            grading.append(gradingi)

        bc_map = {
            'wall' : 0,
            'patch' : 1,
            'symmetry' : 2,
        }
        iname_to_quads = defaultdict(list)
        inames = []
        bcs = []
        iname = 1
        iname_to_name = {}
        iname_to_type = {}
        quads = []
        for key, boundary in boundaries.items():
            #print key, boundary.keys()
            Type = boundary['type']
            bc = bc_map[Type]
            #bcs.append(bc)
            faces = boundary['faces']
            iname_to_name[iname] = key
            iname_to_type[iname] = Type
            for unused_iface, face in faces.items():
                quad = face.strip('() ').split()
                quad = [int(i) for i in quad]
                quads.append(quad)
                iname_to_quads[iname].append(quad)
                inames.append(iname)
                bcs.append(bc)
            self.log.info('iname=%s -> %s; bc=%s -> %s' % (iname, key, bc, Type))
            iname += 1

        nodes = np.array(nodes, dtype='float32')
        self.nodes = nodes

        hexas = np.array(hexas, dtype='int32')
        npoints = np.array(npoints, dtype='int32')
        grading = np.array(grading, dtype='float32')
        self.hexas = hexas
        self.npoints = npoints
        self.grading = grading

        quads = np.array(quads, dtype='int32')
        inames = np.array(inames, dtype='int32')
        bcs = np.array(bcs, dtype='int32')

        self.quads = quads
        self.iname_to_quads = iname_to_quads
        self.inames = inames
        self.bcs = bcs
        self.iname_to_name = iname_to_name
        self.iname_to_type = iname_to_type

        #bdf_filename = 'blockMesh.bdf'
        #self.write_bdf(bdf_filename, nodes, hexas)
        return nodes, hexas, quads, inames, bcs

    def write_bdf(self, bdf_filename, nodes, hexas):
        """writes the BlockMesh as a Nastran BDF"""
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            for inode, node in enumerate(nodes):
                (x, y, z) = node
                fields = ['GRID', inode + 1, None, float(x), float(y), float(z)]
                bdf_file.write(print_card_8(fields))

            pid = 1
            for ielement, hexa in enumerate(hexas):
                bdf_file.write(print_card_8(['CHEXA', ielement + 1, pid,] + list(hexa)))

            bdf_file.write('PSOLID, 1, 1\n')
            bdf_file.write('MAT1, 1, 1.0,,0.3\n')
            bdf_file.write('ENDDATA\n')
        #print(write_dict(d, baseword='BlockMesh'))

    def adjust_nodes_to_symmetry(self):
        # find where nodes have -y value
        unused_y_locations = self.nodes[:, 1]

        #i = where(y_locations < -0.0508)[0]
        #self.nodes[i, 1] = -0.0508

        # set -y to 0
        #i = where(y_locations < 0.0)[0]
        #self.nodes[i, 1] = 0.0
        neq = 1
        ieq = 0
        while neq > 0:
            neq = self.equivalence_nodes(0.0001)
            #self.write_block_mesh('blockMeshDict_n' + str(neq), make_symmetry=True)
            ieq += 1
            if ieq == 10:
                asdf
                break
            self.log.debug('neq = %s' % neq)
        self.log.debug('ieq = %s' % ieq)

    def equivalence_nodes(self, Rtol=0.1):
        """equivalences the nodes of a BlockMeshDict"""
        #same_location = []
        inode_map = {}
        neq = 0
        nodes = self.nodes
        nnodes = self.nodes.shape[0]
        self.log.info('nnodes = %s' % nnodes)
        for inode, nodei in enumerate(self.nodes):
            for jnode, nodej in enumerate(self.nodes):
                if inode < jnode:
                    delta = norm(nodei-nodej)
                    if delta < Rtol:
                        #print('delta=%6g %2s %2s' % (delta, inode, jnode))
                        #same_location.append([inode, jnode])
                        inode_map[jnode] = inode
                        neq += 1
        #same_location = array(same_location, dtype='int32')
        #print('same_location = \n%s' % same_location)
        #nleft = nnodes - neq
        self.log.info('neq = %s' % neq)
        if neq > 15:
            asdf

        # get the node mapping as a dict
        new_ids_map = {}
        new_ids_map2 = {}
        for i in range(nnodes):
            if i in inode_map:  # nodes to adjust
                new_ids_map[i] = inode_map[i]
            else:
                new_ids_map[i] = i

        stack_nodes = set(self.hexas.flatten())
        self.log.info(stack_nodes)
        for iname, faces in self.iname_to_quads.items():
            for face in faces:
                stack_nodes.update(set(list(face)))

        #if 0:
            #print(stack_nodes, len(stack_nodes))
            #for key, value in sorted(new_ids_map.items()):
                #if key not in stack_nodes:
                    #print(' #k=%s v=%s' % (key, value))
                #if key != value:
                    #print(' *k=%s v=%s' % (key, value))
                #else:
                    #print('  k=%s v=%s' % (key, value))

        # map them to 0...n
        i = 0
        j = 0
        j0 = 0
        self.log.info('-------------------')
        for key, value in sorted(new_ids_map.items()):  # the dict of collapsed nodes
            if key not in stack_nodes:
                j += 1
                continue

            if key == value:
                new_ids_map2[i+j] = value - j
                i += 1
                j0 = j
            else:
                new_ids_map2[i+j] = value - j0
                j += 1

        #for key, value in sorted(new_ids_map2.items()):
            #self.log.info('  k=%s v=%s' % (key, value))
        new_ids_map = new_ids_map2

        # get new array of nodes
        #not_indexes = get_not_indexes(self.nodes, same_location[:, 1])
        #nodes2 = self.nodes[not_indexes, :]

        #a[not_indices] = 888
        #print('nodes2 =', nodes2)
        nodes2 = []
        nodes2_written = []
        for inode, jnode in sorted(new_ids_map2.items()):
            if jnode not in nodes2_written:
                node2 = nodes[inode, :]
                nodes2.append(node2)
                nodes2_written.append(jnode)

        nodes2 = np.array(nodes2, dtype='float32')
        nnodes2 = nodes2.shape[0]
        assert nnodes2 <= nnodes, 'nnodes=%s nnodes2=%s' % (nnodes, nnodes2)

        # update the hexas
        hexas2 = []
        npoints2 = []  # don't change
        grading2 = []  # don't change
        for unused_ihexa, hexa, npointsi, gradingi in zip(count(), self.hexas,
                                                          self.npoints, self.grading):
            hexa2 = []
            for j in hexa:
                i = new_ids_map[j]
                hexa2.append(i)
            #print hexa2
            #asdf
            if len(np.unique(hexa2)) == 8:

                hexas2.append(hexa2)
                npoints2.append(npointsi)
                grading2.append(gradingi)
        hexas2 = np.array(hexas2, dtype='int32')
        npoints2 = np.array(npoints2, dtype='int32')
        grading2 = np.array(grading2, dtype='float32')


        # update the faces
        iname_to_quads = {}
        for iname, faces in self.iname_to_quads.items():
            faces2 = []
            for face in faces:
                face2 = []
                for j in face:
                    i = new_ids_map[j]
                    face2.append(i)
                faces2.append(face2)
            iname_to_quads[iname] = np.array(faces2, dtype='int32')

        #self.log.info(nodes.shape)
        #self.log.info(nodes2.shape)
        self.nodes = nodes2
        self.hexas = hexas2
        self.grading = grading2
        self.npoints = npoints2
        self.iname_to_quads = iname_to_quads

        return neq

    def write_block_mesh(self, block_mesh_name_out='blockMeshDict.out', make_symmetry=False):
        """filename interface to ``_write_block_mesh``"""
        with open(block_mesh_name_out, 'w') as block_mesh_file:
            self.log.info('writing %s' % block_mesh_name_out)
            self._write_block_mesh(block_mesh_file, make_symmetry)

    def _write_block_mesh(self, block_mesh_file, make_symmetry):
        """writes a BlockMeshDict with a file object"""
        nodes = self.nodes
        hexas = self.hexas

        header = (
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\    /   O peration     | Version:  2.2.2                                 |\n'
            '|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n'
            '|    \\/     M anipulation  |                                                 |\n'
            r'\*---------------------------------------------------------------------------*/' +
            '\n' +
            'FoamFile\n'
            '{\n'
            '    version     0.0508;\n'
            '    format      ascii;\n'
            '    class       dictionary;\n'
            '    object      blockMeshDict;\n'
            '}\n'
            '\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n'
        )
        block_mesh_file.write(header)

        block_mesh_file.write('convertToMeters 1.0;\n')
        block_mesh_file.write('\n')
        block_mesh_file.write('vertices\n')
        block_mesh_file.write('(\n')

        unique_x = np.unique(nodes[:, 0])
        unique_y = np.unique(nodes[:, 1])
        unique_z = np.unique(nodes[:, 2])
        unique_x.sort()
        unique_y.sort()
        unique_z.sort()
        #print('unique_x = %s' % unique_x)
        self.log.info('unique_y = %s' % unique_y)
        #print('unique_z = %s' % unique_z)

        nnodes = nodes.shape[0]
        fmt_node = '%%%si' % len(str(nnodes))
        #print('fmt-node = %r' % fmt_node)
        for inode, (x, y, z) in enumerate(nodes):
            block_mesh_file.write('    (%7s %7s %7s) // %s\n' % (x, y, z, inode))
            if (inode + 1) % 5 == 0:
                block_mesh_file.write('\n')
        block_mesh_file.write(');\n\n')

        block_mesh_file.write('blocks\n')
        block_mesh_file.write('(\n')
        block_mesh_file.write('      // hex                        npoints in each dir;    stretching\n')
        #print "nodes = ", nodes.shape
        hexai_fmt = 'hex (%s %s %s %s %s %s %s %s)' % tuple([fmt_node] * 8)

        m_to_inch = 39.3701
        for ihexa, (hexa, npointsi, gradingi) in enumerate(zip(hexas, self.npoints, self.grading)):
            save_element = False
            for inode in hexa:
                #print('inode = %s' % inode)
                node = nodes[inode, :]
                #print('nodes[%s] = %s' % (inode, str(node)))
                if node[1] >= 0.0:
                    save_element = True

            #npointsi = self.npoints[ihexa]
            #gradingi = self.grading[ihexa]
            save_element = True
            if save_element:
                voli = volume8(
                    nodes[hexa[0], :], nodes[hexa[1], :], nodes[hexa[2], :], nodes[hexa[3], :],
                    nodes[hexa[4], :], nodes[hexa[5], :], nodes[hexa[6], :], nodes[hexa[7], :],)
                if voli > 0:
                    shexai = hexai_fmt % tuple(hexa)
                    snpointsi = '(%s %s %s)' % tuple(npointsi)
                    sgradingi = 'simpleGrading (%g %g %g)' % tuple(gradingi)
                    svoli = 'vol = %g [in]' % (voli * m_to_inch**3)
                    block_mesh_file.write('      %s %s %s  // %s\n' % (
                        shexai, snpointsi, sgradingi, svoli))

            if (ihexa + 1) % 4 == 0:
                block_mesh_file.write('\n')
        block_mesh_file.write(');\n\n')

        block_mesh_file.write('edges\n')
        block_mesh_file.write('(\n')
        block_mesh_file.write(');\n')

        block_mesh_file.write('\n')
        block_mesh_file.write('boundary\n')
        block_mesh_file.write('(\n')

        #iname_quads = {}

        symmetry_faces = []
        facei_fmt = '           (%s %s %s %s)' % tuple([fmt_node] * 4)
        for iname in sorted(self.iname_to_quads):
            name = self.iname_to_name[iname]
            Type = self.iname_to_type[iname]
            block_mesh_file.write('    %s\n' % name)
            block_mesh_file.write('    {\n')
            block_mesh_file.write('        type %s;\n' % Type)
            block_mesh_file.write('        faces\n')
            block_mesh_file.write('        (\n')
            faces = self.iname_to_quads[iname]
            for face in faces:
                n1 = nodes[face[0], :]
                n2 = nodes[face[1], :]
                n3 = nodes[face[2], :]
                n4 = nodes[face[3], :]
                area, centroid = area_centroid(n1, n2, n3, n4)
                centroid *= m_to_inch
                if area > 0.:
                    area *= m_to_inch ** 2
                    y_centroid = centroid[1]
                    if make_symmetry and np.allclose(y_centroid, 0.0):
                        symmetry_faces.append(face)
                    else:
                        if y_centroid <= 0.0:
                            block_mesh_file.write(
                                facei_fmt % tuple(face) +
                                ' // centroid=(%3g, %3g, %3g) [in]  Area=%.2f [in^2]\n' % (
                                    centroid[0], centroid[1], centroid[2], area))
                        else:
                            block_mesh_file.write(facei_fmt % tuple(face) + '\n') #+ ' // c=(%7s, %7s, %7s)\n' % tuple(centroid)
            block_mesh_file.write('        );\n')
            block_mesh_file.write('    }\n')
            #save_face = False

        if make_symmetry:
            #name = self.iname_to_name[iname]
            name = 'symmetry'
            foam_type = 'symmetryPlane'
            #foam_type = self.iname_to_type[iname]

            block_mesh_file.write('    %s\n' % name)
            block_mesh_file.write('    {\n')
            block_mesh_file.write('        type %s;\n' % foam_type)
            block_mesh_file.write('        faces\n')
            block_mesh_file.write('        (\n')
            for face in symmetry_faces:
                n1 = nodes[face[0], :]
                n2 = nodes[face[1], :]
                n3 = nodes[face[2], :]
                n4 = nodes[face[3], :]
                area, centroid = area_centroid(n1, n2, n3, n4)
                centroid *= m_to_inch
                area *= m_to_inch ** 2

                y_centroid = centroid[1]
                if y_centroid <= 0.0:
                    block_mesh_file.write(
                        facei_fmt % tuple(face) +
                        ' // centroid=(%3g, %3g, %3g) [in]  Area=%.2f [in^2]\n' % (
                            centroid[0], centroid[1], centroid[2], area))
                else:
                    block_mesh_file.write(facei_fmt % tuple(face) + '\n')
                    #+ ' // c=(%7s, %7s, %7s)\n' % tuple(centroid)

            block_mesh_file.write('        );\n')
            block_mesh_file.write('    }\n')

        block_mesh_file.write(');\n\n')
        block_mesh_file.write('mergeMatchPairs\n')
        block_mesh_file.write('(\n')
        block_mesh_file.write(');\n\n')
        block_mesh_file.write('// ************************************************************************* //\n')

def mirror_block_mesh(block_mesh_name, block_mesh_name_out, log=None, debug=True):
    """mirrors a blockMeshDict"""
    make_symmetry = True
    block_mesh_model = read_block_mesh(block_mesh_name, log=log, debug=debug)
    #out = block_mesh_model.read_openfoam(block_mesh_name)
    #unused_nodes, unused_hexas, unused_quads, unused_names, unused_bcs = out
    if make_symmetry:
        block_mesh_model.adjust_nodes_to_symmetry()
    block_mesh_model.write_block_mesh(block_mesh_name_out, make_symmetry=make_symmetry)

def main():  # pragma: no cover
    import sys
    block_mesh_name = sys.argv[1]
    block_mesh_name_out = sys.argv[2]
    mirror_block_mesh(block_mesh_name, block_mesh_name_out)

if __name__ == '__main__':  # pragma: no cover
    main()
