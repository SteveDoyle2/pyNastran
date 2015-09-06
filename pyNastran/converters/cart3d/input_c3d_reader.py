from __future__ import print_function
from six import iteritems
from collections import defaultdict
from numpy import array, linspace, hstack, vstack
from pyNastran.bdf.cards.aero import points_elements_from_quad_points


class InputC3dReader(object):
    def __init__(self):
        pass

    def read_input_c3d(self, input_c3d_filename, stack=True):
        with open(input_c3d_filename, 'r') as c3d_filename:
            lines = c3d_filename.readlines()

        packs = defaultdict(list)
        is_comment = False
        for line in lines:
            line = line.strip()
            #print(is_comment, line)
            if '"' in line:
                is_comment = not is_comment
                continue
            if is_comment or line.startswith('+'):
                continue

            line = line.split('#')[0].strip()
            if not line:
                continue
            #if line[0].isalpha():
                #print('# %s' % line)
                #continue

            print('%r' % line)
            if len(line) > 2 and line[0].isdigit() and line[1] == '.':
                pack_key = int(line.split('.')[0])
            else:
                pack_values = line.split()
                packs[pack_key].append(pack_values)

        for key, values in sorted(iteritems(packs)):
            print(key)
            for value in values:
                print('    %s' % value)
        print(packs[1])
        geometry_filename = packs[1][0]
        xyz_size = [float(val) for val in packs[2][1]]
        xyz_dim = [int(val) for val in packs[3][0]]
        return self.create_nodes_elements(xyz_size, xyz_dim, stack=stack)

    def create_nodes_elements(self, xyz_size, xyz_dim, stack=True):
        xmin, xmax, ymin, ymax, zmin, zmax = xyz_size
        xdim, ydim, zdim = xyz_dim
        x = linspace(0., 1., num=xdim)
        y = linspace(0., 1., num=ydim)
        z = linspace(0., 1., num=zdim)
        print(x)

        # organized in x, y, z planes order
        planes = (
            # yz
            # xmin, xmax
            ((xmin, ymin, zmin), (xmin, ymin, zmax), (xmin, ymax, zmax), (xmin, ymax, zmin), y, z),
            ((xmax, ymin, zmin), (xmax, ymin, zmax), (xmax, ymax, zmax), (xmax, ymax, zmin), y, z), # wrong

            ## xz
            # ymin, ymax
            ((xmin, ymin, zmin), (xmin, ymin, zmax), (xmax, ymin, zmax), (xmax, ymin, zmin), x, z),
            ((xmin, ymax, zmin), (xmin, ymax, zmax), (xmax, ymax, zmax), (xmax, ymax, zmin), x, z),

            # xy
            # zmin, zmax
            ((xmin, ymin, zmin), (xmin, ymax, zmin), (xmax, ymax, zmin), (xmax, ymin, zmin), x, y),
            ((xmin, ymin, zmax), (xmin, ymax, zmax), (xmax, ymax, zmax), (xmax, ymin, zmax), x, y),
        )
        points = []
        elements = []
        nnodes = 0
        for plane in planes:
            p1, p2, p3, p4, xi, yi = plane
            p1 = array(p1, dtype='float32')
            p2 = array(p2, dtype='float32')
            p3 = array(p3, dtype='float32')
            p4 = array(p4, dtype='float32')
            print(plane[:1])

            pointsi, elementsi = points_elements_from_quad_points(p1, p2, p3, p4, xi, yi)
            points.append(pointsi)
            if stack:
                elements.append(elementsi + nnodes)
                nnodes += pointsi.shape[0]
            else:
                elements.append(elementsi)
            #break
            #print(points)

        if stack:
            points = vstack(points)
            elements = vstack(elements)
        return points, elements


def main():
    input_c3d_filename = r'F:\work\pyNastran\pyNastran\master2\pyNastran\converters\cart3d\models\bJet\input.c3d'
    c3d = InputC3dReader()
    nodes, elements = c3d.read_input_c3d(input_c3d_filename)

if __name__ == '__main__':
    main()