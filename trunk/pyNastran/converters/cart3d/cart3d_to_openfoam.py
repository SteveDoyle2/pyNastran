from pyNastran.converters.cart3d.cart3d_reader import Cart3DInpt
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename

def cart3d_to_openfoam(cart3d_filename, inpt_filename):
    cart3d_to_stl_filename(cart3d_filename, 'cart3d.stl')
    create_openfoam_inputs(inpt_filename)

def create_openfoam_inputs(inpt_filename):
    inpt = read_cart3d_inpt(int_filename)
    inpt.write_openfoam()

if __name__ == '__main__':  # pragma: no cover
    #bdf_filename = 'g278.bdf'
    cart3d_geo_filename = 'g278.tri'
    cart3d_inpt_filename = 'g278.inpt'
    cart3d_to_openfoam(cartcart3d_geo_filename, cart3d_inpt_filename)