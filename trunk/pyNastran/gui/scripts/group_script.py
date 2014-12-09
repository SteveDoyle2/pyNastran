
import pyNastran
pkg_path = pyNastran.__path__[0]

self.save_data = True
if 0:
    bdf_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf')
    op2_filename = os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.op2')


    self.on_load_geometry(infile_name=bdf_filename, geometry_format='nastran')

    element_id = range(1, 105)
    #self.create_group('nastran', '0', element_id=element_id)  # panel ID
    self.create_group('nastran', 'A-Nastran', element_id=element_id, coord_id=3, show=False)
    self.post_groups(['A-Nastran'])

if 1:
    cart3d_filename = os.path.join(pkg_path, 'converters', 'cart3d', 'threePlugs.a.tri')
    self.on_load_geometry(infile_name=cart3d_filename, geometry_format='cart3d')

    element_id = range(1, 6000)
    self.create_group('cart3d', 'A-Cart3d', element_id=element_id, show=False)
    self.post_groups(['A-Cart3d'])

