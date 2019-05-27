from copy import deepcopy
from docopt import docopt

import pyNastran
from pyNastran.converters.stl.stl import read_stl

def main():
    # --ascii <fmt> doesn't work for binary properly, so we set:
    # --ascii False
    msg = ''
    msg += 'Usage:\n'
    msg += '  stl_reshape [--ascii <fmt>] <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--ascii <fmt>] [--xy | --yz | --xz] [--scale S] <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--ascii <fmt>] [--xscale X] [--yscale Y] [--zscale Z] <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--ascii <fmt>] --delta <xshift> <yshift> <zshift> <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--ascii <fmt>] --mirror <xyz> <tol> <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--ascii <fmt>] --flip_normals <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--ascii <fmt>] --stats <in_stl_filename>\n'
    #msg += '  stl_reshape [--ascii <fmt>] <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape -h | --help\n'
    msg += '  stl_reshape -v | --version\n'
    msg += '\n'
    msg += 'Options:\n'
    msg += '  -h, --help        show this help message and exit\n'
    msg += "  -v, --version     show program's version number and exit\n"
    msg += '  --ascii           writes the model is ascii\n'
    msg += '  fmt               the sscii format; False -> binary\n'

    msg += '  --scale S         scale the xyz values by S\n'
    msg += '  --xy              flip the x and y axes (pick one)\n'
    msg += '  --yz              flip the y and z axes (pick one)\n'
    msg += '  --xz              flip the x and z axes (pick one)\n'

    msg += '  --xscale X        scale the x values by X (use \'-1.0\' if < 0)\n'
    msg += '  --yscale Y        scale the y values by Y (use \'-1.0\' if < 0)\n'
    msg += '  --zscale Z        scale the z values by Z (use \'-1.0\' if < 0)\n'

    msg += '  --delta           shift the model some <xshift,yshift,shift>\n'
    msg += '  xshift            shift the model by xshift (use \'-1.0\' if < 0)\n'
    msg += '  yshift            shift the model by yshift (use \'-1.0\' if < 0)\n'
    msg += '  zshift            shift the model by zshift (use \'-1.0\' if < 0)\n'

    msg += '  --mirror          create a mirror model\n'
    msg += '  xyz               the x, y, or z direction to mirror about\n'
    msg += '  tol               the tolerance\n'

    msg += '  --flip_normals    flip the normals of the elements\n'

    msg += '  --stats           print the min/max locations\n'

    msg += '  in_stl_filename   the input filename\n'
    msg += '  out_stl_filename  the output filename\n'

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    print(data)
    stl_reshape(data)


def stl_reshape(data):
    if '--xy' not in data:
        data['--xy'] = False
    if '--yz' not in data:
        data['--yz'] = False
    if '--xz' not in data:
        data['--xz'] = False
    if '--scale' not in data:
        data['--scale'] = None

    if '--xscale' not in data:
        data['--xscale'] = None
    if '--yscale' not in data:
        data['--yscale'] = None
    if '--zscale' not in data:
        data['--zscale'] = None

    if '<xshift>' not in data:
        data['<xshift>'] = None
    if '<yshift>' not in data:
        data['<yshift>'] = None
    if '<zshift>' not in data:
        data['<zshift>'] = None

    if '--stats' not in data:
        data['--stats'] = None
    if '--mirror' not in data:
        data['--mirror'] = None
    if '--flip_normals' not in data:
        data['--flip_normals'] = None

    in_stl_filename = data['<in_stl_filename>']
    out_stl_filename = data['<out_stl_filename>']
    assert in_stl_filename != out_stl_filename

    stl = read_stl(in_stl_filename)

    if data['<fmt>'] in ['False', False]:
        is_binary = True
        fmt = None
    else:
        fmt = data['<fmt>']
        is_binary = False
    print('is_binary=%s' % is_binary)

    if data['--xy'] or data['--yz'] or data['--xz']:
        scale = 1.
        if data['--scale'] is not None:
            scale = float(data['--scale'])

        if data['--xy']:
            assert data['--yz'] is False
            assert data['--xz'] is False
            axes = 'xy'
        elif data['--yz']:
            assert data['--xy'] is False
            assert data['--xz'] is False
            axes = 'yz'
        elif data['--xz']:
            assert data['--xy'] is False
            assert data['--yz'] is False
            axes = 'xz'
        #print('flip_axes = %r' % axes)
        #print(data)
        stl.flip_axes(axes, scale)

    elif data['--xscale'] or data['--yscale'] or data['--zscale']:
        xscale = 1.
        yscale = 1.
        zscale = 1.
        if data['--xscale'] is not None:
            xscale = float(data['--xscale'].strip("'"))
        if data['--yscale'] is not None:
            yscale = float(data['--yscale'].strip("'"))
        if data['--zscale'] is not None:
            zscale = float(data['--zscale'].strip("'"))
        x = deepcopy(stl.nodes[:, 0])
        y = deepcopy(stl.nodes[:, 1])
        z = deepcopy(stl.nodes[:, 2])
        stl.nodes[:, 0] = x * xscale
        stl.nodes[:, 1] = y * yscale
        stl.nodes[:, 2] = z * zscale
    elif data['<xshift>'] or data['<yshift>'] or data['<zshift>']:
        xshift = 1.
        yshift = 1.
        zshift = 1.
        if data['<xshift>'] is not None:
            if isinstance(xshift, str):
                xshift = float(data['<xshift>'].strip("'"))
            else:
                xshift = float(data['<xshift>'])

        if data['<yshift>'] is not None:
            if isinstance(xshift, str):
                yshift = float(data['<yshift>'].strip("'"))
            else:
                yshift = float(data['<yshift>'])

        if data['<zshift>'] is not None:
            if isinstance(xshift, str):
                zshift = float(data['<zshift>'].strip("'"))
            else:
                zshift = float(data['<zshift>'])

        print('delta = (%s, %s, %s)' % (xshift, yshift, zshift))
        stl.nodes[:, 0] += xshift
        stl.nodes[:, 1] += yshift
        stl.nodes[:, 2] += zshift
    elif data['--scale']:
        scale = float(data['--scale'])
        stl.nodes *= scale
    elif data['--stats']:
        xmax, ymax, zmax = stl.nodes.max(axis=0)
        xmin, ymin, zmin = stl.nodes.min(axis=0)
        print('xyz_max = (%g, %g, %g)' % (xmax, ymax, zmax))
        print('xyz_min = (%g, %g, %g)' % (xmin, ymin, zmin))
        return
    elif data['--mirror']:
        #plane = data['plane']
        #assert plane in ['xy', 'yz', 'xz'], 'plane=%r' % plane
        xyz = data['<xyz>']
        tol = float(data['<tol>'])
        stl.create_mirror_model(xyz, tol)
    elif data['--flip_normals']:
        stl.flip_normals()
    else:
        raise RuntimeError('unsupported reshape...data=%s' % data)

    stl.write_stl(out_stl_filename, is_binary=is_binary, float_fmt=fmt)

if __name__ == '__main__':  # pragma: no cover
    main()
