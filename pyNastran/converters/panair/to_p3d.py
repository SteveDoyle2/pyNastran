import os
import pyNastran
from pyNastran.converters.panair.panair_grid import PanairGrid
PKG_PATH = pyNastran.__path__[0]

def run():
    infilename = os.path.join(PKG_PATH, 'converters', 'panair', 'M100', 'M100.inp')
    model = PanairGrid(log=None, debug=True)
    model.read_panair(infilename)

    p3d_name = 'M100.plt'
    model.write_plot3d(p3d_name)

if __name__ == '__main__':  # pragma: no cover
    run()
