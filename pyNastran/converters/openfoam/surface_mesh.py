from numpy import zeros
from pyNastran.converters.openfoam.openfoam_parser import (
    #write_dict,
    FoamFile, convert_to_dict)
from pyNastran.bdf.field_writer import print_card_8


class Points:
    def __init__(self):
        foam_points = FoamFile('points')
        lines = foam_points.read_foam_file()
        self.points_dict = convert_to_dict(self, lines)

    def read_points(self):
        #print write_dict(d, baseword='Points = ')
        keys = self.points_dict.keys()

        ifoam = keys.index('FoamFile')
        keys.pop(ifoam)
        assert len(keys) == 1, keys
        npoints = keys[0]
        print("npoints = ", npoints)

        nodes = self.points_dict[npoints]
        #f = open('points2.bdf', 'wb')
        #f.write('CEND\n')
        #f.write('BEGIN BULK\n')
        nnodes = len(nodes)
        nodes_array = zeros((nnodes, 3), dtype='float32')
        for inode, node in enumerate(nodes):
            x, y, z = node.strip('() ').split()
            nodes_array[inode, :] = [x, y, z]
            #f.write(print_card_8(['GRID', inode + 1, None, float(x), float(y), float(z)]))
            #f.write(print_card_8(['CONM2', inode + 1, inode + 1, 0, 0., 0., 0.]))
        #f.write(print_card_8(['CQUAD4', inode + 10, inode + 10, 1, 2, 3, 4]))
        #f.write('PSHELL, 1, 1, 0.1\n')
        #f.write('MAT1, 1, 1.0,,0.3\n')
        #f.write('ENDDATA\n')
        #f.close()
        return nodes_array

class Faces:
    def __init__(self):
        foam_file = FoamFile('faces')
        lines = foam_file.read_foam_file()
        self.foam_dict = convert_to_dict(self, lines)

    def read_faces(self):
        keys = self.foam_dict.keys()

        ifoam = keys.index('FoamFile')
        keys.pop(ifoam)
        assert len(keys) == 1, keys
        nfaces = keys[0]
        print("nfaces = ", nfaces)

        faces = self.foam_dict[nfaces]
        nfaces = len(faces)
        faces_array = zeros((nfaces, 4), dtype='int32')
        for iface, face in enumerate(faces):
            #4(11 132 133 12)
            n, face2 = face.strip(') ').split('(')
            n = n.strip()
            #print face
            #print face2
            #print
            sface = face2.split()
            #print sface
            assert n == '4', n
            faces_array[iface, :] = sface
        return faces_array

def main():
    foam_points = Points()
    nodes = foam_points.read_points()
    unused_nnodes = nodes.shape[0]

    foam_faces = Faces()
    quads = foam_faces.read_faces()
    quads += 1

    with open('points_faces.bdf', 'wb') as bdf_file:
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')

        for inode, node in enumerate(nodes):
            card = ['GRID', inode + 1, None] + list(node)
            bdf_file.write(print_card_8(card))

        ielement = 1
        pid = 1
        for quad in quads:
            card = ['CQUAD4', ielement, pid] + list(quad)
            bdf_file.write(print_card_8(card))
            ielement += 1
        bdf_file.write('PSHELL, 1, 1, 0.1\n')
        bdf_file.write('MAT1, 1, 1.0,,0.3\n')
        bdf_file.write('ENDDATA\n')

if __name__ == '__main__':   # pragma: no cover
    main()
