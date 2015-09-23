import os

from numpy import array
# from pyNastran.converters.cart3d.cart3d_reader import Cart3DInpt
from pyNastran.converters.cart3d.cart3d import Cart3d
from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename

def get_header():
    header = ''
    header += '/*------------OpenFOAM------------*\\n'
    header += ' |                                |\n'
    header += ' |       this is a header         |\n'
    header += ' |                                |\n'
    header += '/*------------OpenFOAM------------*\\n'
    return header

def write_foam_file(version, obj, location=None):
    msg = ''
    msg += 'FoamFile\n'
    msg += '{\n'
    msg += '    version     %s;\n' % version
    msg += '    format      ascii;\n'
    msg += '    class       dictionary;\n'
    if location is not None:
        msg += '    location    "%s";\n' % location
    msg += '    version     %s;\n' % obj
    msg += '}\n'
    msg += 'convertToMeters 1.0;\n'
    msg += '\n'
    return msg

def cart3d_to_openfoam(cart3d_filename, inpt_filename, basepath):
    """
    https://openfoamwiki.net/images/f/f0/Final-AndrewJacksonSlidesOFW7.pdf
    """
    stl_filename = 'cart3d.stl'
    cart3d_to_stl_filename(cart3d_filename, stl_filename)

    cart3d = Cart3d()
    cart3d.read_cart3d(cart3d_filename)
    points = cart3d.points
    elements = cart3d.elements

    # convert to 0 based
    elements -= 1

    cart3d.void = (1., 0., 0.)

    #create_openfoam_inputs(inpt_filename)
    inpt = Inpt()
    inpt.read_cart3d_inpt(inpt_filename)

    basepath = 'hsct'
    os.makedirs(basepath)
    os.makedirs(os.path.join(basepath, 'constant', 'polyMesh'))
    os.mkdir(os.path.join(basepath, 'system'))

    fsnappy = open(os.path.join(basepath, 'constant', 'snappyHexMeshDict'))
    fblock = open(os.path.join(basepath, 'constant', 'polyMesh', 'blockMeshDict'))
    fcontrol = open(os.path.join(basepath, 'system', 'controlMeshDict'))
    fsnappy.write(get_header())
    fblock.write(get_header())
    fcontrol.write(get_header())

    snappy = ''
    block = ''
    control = ''

    block += write_foam_file(data, 'blockMeshDict', 'constant/polyMesh')
    snappy += write_foam_file(data, 'autoMeshDict')
    control += write_foam_file(data, 'controlMeshDict', 'constant/polyMesh')

    maxR = inpt.maxR
    xmax, ymax, zmax = points.max(axis=1)
    xmin, ymin, zmin = points.min(axis=1)

    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    xmin -= dx * maxR
    ymin -= dy * maxR
    zmin -= dz * maxR

    xmax += dx * maxR
    ymax += dy * maxR
    zmax += dz * maxR

    p = array([
        [xmin, ymin, zmin],
        [xmin, ymin, zmax],
        [xmin, ymax, zmax],
        [xmin, ymax, zmin],

        [xmax, ymin, zmin],
        [xmax, ymin, zmax],
        [xmax, ymax, zmax],
        [xmax, ymax, zmin],
    ])

    block += 'vertices\n'
    block += '{\n'
    for point in p:
        block += '    (%s, %s, %s)\n' % tuple(point)
    block += '};\n'

    block += 'elements\n'
    block += '{\n'
    block += '    hex (0 1 2 3 4 5 6 7) (%s %s %s) (1.0 1.0 1.0)\n' % inpt.xyz_start
    block += '};\n\n'

    snappy += 'castellatedMesh true;\n'
    snappy += 'snap            true;\n'
    snappy += 'addLayers       false;\n'

    block += 'patches\n'
    block += '{\n'
    block += '    inflow\n'
    block += '    (\n'
    block += '       4 (0 1 2 3)\n'
    block += '    outflow\n'
    block += '       4 (4 5 6 7)\n'
    block += '    )\n'
    block += '};\n\n'
    block += 'mergePatchParis();\n'

    assert os.path.exists(stl_filename), stl_filename
    stl_base, ext = os.path.splitext(stl_filename)

    snappy += 'geometry;\n'
    snappy += '{\n'
    snappy += '    %s // STL filename;\n' % os.path.basename(stl_filename)
    snappy += '    {\n'
    snappy += '        type triSurfaceMesh;;\n'
    snappy += '        name %s;\n' % stl_base
    snappy += '    }\n'
    snappy += '};\n'

    snappy += 'snapControls\n'
    snappy += '{\n'
    snappy += '    nSmoothPatch     3;\n'
    snappy += '    tolerance        1.0;\n'
    snappy += '    nSolveIter       300;\n'
    snappy += '    nRelaxIter       5;\n'
    snappy += '    nFeatureSnapIter 10;\n'
    snappy += '};\n'

    inpt.max_cells = 2
    inpt.min_refinement_cells = 1
    inpt.ncells_between_levels = 2
    castel += '{'
    castel += '    maxGlobalCells      %i\n' % inpt.max_cells * 1e6
    castel += '    minRefinementCells  %s\n' % inpt.min_refinement_cells
    castel += '    nCellsBetweenLevels %s\n\n' % inpt.ncells_between_levels

    castel += '    features();\n'
    castel += '}\n\n'

    castel += '// refinementSurfaces'
    castel += 'resolveFeatureAngle        %s\n' % inpt.feature_angle
    castel += 'locationInMesh             (%s %s %s);\n' % tuple(cart3d.void)
    castel += 'allowFreeStandingZoneFaces true;  // ?????\n'
    castel += '}'

    fsnappy.write(snappy)
    fblock.write(block)
    fcontrol.write(control)
    fcastel.write(castel)

    fsnappy.close()
    fblock.close()
    fcontrol.close()
    fcastle.close()


class Inpt(object):
    def __init__(self):
        # base mesh
        self.xyz_start = (5, 5, 5)
        self.maxR = 30

        # refine mesh
        self.angle_refine = 25.
        self.dir_threshold = 10.
        self.sharp_refine = True

        # castle mesh
        self.keep_xyz = (1., 0., 0.)

        # layer mesh
        self.nbuffer_layers = 5

        # aero
        self.alpha = 10.
        self.beta = 0.
        self.Mach = 2.4
        self.S = 100.      # ft^2
        self.c = 50.       # ft
        self.ainf = 1000.  # ft/s

        self.pinf = 14.7 # psi
        self.pinf * 144.  # psf
        self.gamma = 1.4
        self.qinf = self.gamma / 2 * self.pinf * self.Mach**2  # psf
        self.Vinf = self.Mach * self.qinf


    def read_cart3d_inpt(self, inpt_filename):
        pass


def create_openfoam_inputs(inpt_filename):
    inpt = read_cart3d_inpt(inpt_filename)


def main():
    #bdf_filename = 'g278.bdf'
    cart3d_geo_filename = 'g278.tri'
    cart3d_inpt_filename = 'g278.inpt'
    cart3d_to_openfoam(cartcart3d_geo_filename, cart3d_inpt_filename)


if __name__ == '__main__':  # pragma: no cover
    main()
