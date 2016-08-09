from __future__ import print_function
from numpy import array, unique, hstack, zeros

class OBJ(object):
    def __init__(self):
        pass

    def read_obj(self, obj_filename):
        f = open(obj_filename, 'r')

        # v -0.0817245 0.000635 0.00421862
        # v -0.0817245 0.000580371 0.00421862
        # v -0.0817245 -0.000635 0.00421862
        # l 1 2
        # l 2 3

        nodes = []
        lines = []
        #faces = []
        for line in f:
            sline = line.strip().split()
            #print(sline)
            Type = sline[0]
            if Type == 'v':  # vertex
                nodes.append(sline[1:])
            elif Type == 'l':  # line
                lines.append(sline[1:])
            #elif Type == 'vt':  # texture coordinate
                #lines.append(sline[1:])
            #elif Type == 'vn':  # normal vector (not unit vector)
                #lines.append(sline[1:])
            #elif Type == 'vp':  # parameter space vertex
                #lines.append(sline[1:])
            else:
                raise NotImplementedErrror(sline)
        self.nodes = array(nodes, dtype='float64')

        # make it 0-based instead of 1 based
        self.lines = array(lines, dtype='int32') - 1

        self.make_elements()

    def make_elements(self):
        print(self.nodes.shape)
        unodes, indicies = unique_rows(self.nodes, return_inverse=True)

        print(unodes)
        #print(list(indicies))
        print(unodes.shape)
        print(indicies.shape)

        n1 = self.lines[:, 0]
        n2 = self.lines[:, 1]
        i1 = indicies[n1]
        i2 = indicies[n2]
        nrows = len(i1)
        #self.lines = hstack([i1, i2], dtype='int32')
        self.lines = hstack([i1, i2])
        lines2 = zeros((nrows, 2), dtype='int32')
        lines2[:, 0] = i1
        lines2[:, 1] = i2
        self.lines = lines2
        print(self.lines.shape)
        self.nodes = unodes

    def write_obj(self, obj_filename):
        f = open(obj_filename, 'wb')

        float_fmt = '8.6f'
        int_fmt = 'i'
        node_fmt = 'v %%%s %%%s %%%s\n'  % (float_fmt, float_fmt, float_fmt)
        line_fmt = 'l %%%s %%%s\n'  % (int_fmt, int_fmt)
        print(node_fmt)
        for node in self.nodes:
            f.write(node_fmt % tuple(node))
        for line in self.lines + 1:
            f.write(line_fmt % tuple(line))



def unique_rows(data, return_inverse=False):
    ncols = data.shape[1]
    dtype = data.dtype.descr * ncols
    struct = data.view(dtype)

    uniq, indicies = unique(struct, return_inverse=return_inverse)
    uniq = uniq.view(data.dtype).reshape(-1, ncols)
    return uniq, indicies

def main():
    obj_filename = '6.5e-06_edges.txt'
    obj = OBJ()
    obj.read_obj(obj_filename)
    obj.write_obj('b.txt')

if __name__ == '__main__':
    main()
