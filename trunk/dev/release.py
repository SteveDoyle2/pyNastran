import pyNastran
pkg_name = pyNastran.__path__[0]
from shutil import rmtree

rmtree(pkg_name, 'nastranwrapper-using_pynastran')
rmtree(pkg_name, 'pch')
rmtree(pkg_name, 'latex')
rmtree(pkg_name, 'converters', 'calculix')
rmtree(pkg_name, 'converters', 'cgns')
rmtree(pkg_name, 'converters', 'iges')
rmtree(pkg_name, 'converters', 'plot3d')
rmtree(pkg_name, 'converters', 'tecplot'
os.remove(pkg_name, 'converters', 'toCodeAster.py')
os.remove(pkg_name, 'converters', 'toSfepy.py')
os.remove(pkg_name, 'converters', 'nastran', 'bdf_to_p3d.py')

rmtree(pkg_name, 'applications')
#rmtree(pkg_name, 'applications', 'solver')
#rmtree(pkg_name, 'applications', 'hyper')
#rmtree(pkg_name, 'applications', 'contact')
#rmtree(pkg_name, 'applications', 'aero_tet_mesh')
