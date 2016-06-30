import os
import pyNastran
from pyNastran.converters.panair.panair_grid import PanairGrid
pkg_path = pyNastran.__path__[0]

def run():
    infileName = os.path.join(pkg_path, 'converters', 'panair', 'M100', 'M100.inp')
    model = PanairGrid(log=None, debug=True)
    model.read_panair(infileName)

    p3d_name = 'M100.plt'
    model.write_plot3d(p3d_name)

if __name__ == '__main__':  # pragma: no cover
    run()
