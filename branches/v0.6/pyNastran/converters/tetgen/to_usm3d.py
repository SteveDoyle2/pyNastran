from struct import pack, unpack
from numpy import array, transpose, zeros

from pyNastran.utils.log import get_logger

from pyNastran.converters.tetgen.tetgen_reader import TetgenReader
from pyNastran.converters.usm3d.usm3d_reader import write_usm3d


def main():
    m = TetgenReader()
    base = 'tetgen_test_flipped.1'
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')
    write_usm3d(base)

    #model = Usm3dReader()
    #basename = 'HSCT_inviscid'
    #basename = 'box'
    #model.read_usm3d(basename)
    #model.write_usm3d(basename + '_2')


if __name__ == '__main__':
    main()