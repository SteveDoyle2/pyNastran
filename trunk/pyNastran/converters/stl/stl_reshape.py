from docopt import docopt
from copy import deepcopy

import pyNastran
from pyNastran.converters.stl.stl_reader import STLReader

def main():
    msg = ''
    msg += 'Usage:\n'
    msg += '  stl_reshape <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--xy | --yz | --xz] [--scale S] <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape [--xscale X] [--yscale Y] [--zscale Z] <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape --delta <xshift> <yshift> <zshift> <in_stl_filename> <out_stl_filename>\n'
    msg += '  stl_reshape -h | --help\n'
    msg += '  stl_reshape -v | --version\n'
    msg += '\n'
    msg += 'Options:\n'
    msg += '  -h, --help     show this help message and exit\n'
    msg += '  --xy            flip the x and y axes (pick one)\n'
    msg += '  --yz            flip the y and z axes (pick one)\n'
    msg += '  --xz            flip the x and z axes (pick one)\n'

    msg += '  --scale S        scale the xyz values by S\n'
    msg += '  --xscale X       scale the x values by X\n'
    msg += '  --yscale Y       scale the y values by Y\n'
    msg += '  --zscale Z       scale the z values by Z\n'

    msg += '  --delta          shift the model some <xshift,yshift,shift>'
    msg += '  xshift       shift the model by xshift\n'
    msg += '  yshift       shift the model by yshift\n'
    msg += '  zshift       shift the model by zshift\n'

    msg += "  -v, --version  show program's version number and exit\n"
    msg += '  in_stl_filename  the input filename\n'
    msg += '  out_stl_filename  the output filename\n'

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)

    print data
    in_stl_filename = data['<in_stl_filename>']
    out_stl_filename = data['<out_stl_filename>']


    stl = STLReader()
    stl.read_stl(in_stl_filename)
    
    if data['--xy'] or data['--yz'] or data['--xz']:
        scale = 1.
        if data['--scale'] is not None:
            scale = float(data['--scale'])

        if data['--xy']:
            assert data['--yz'] is False
            assert data['--xz'] is False
            x = deepcopy(stl.nodes[:, 0])
            y = deepcopy(stl.nodes[:, 1])
            stl.nodes[:, 0] = y * scale
            stl.nodes[:, 1] = x * scale
            #stl.scale_nodes(xscale, yscale, zscale)
        elif data['--yz']:
            assert data['--xy'] is False
            assert data['--xz'] is False
            y = deepcopy(stl.nodes[:, 1])
            z = deepcopy(stl.nodes[:, 2])
            stl.nodes[:, 1] = z * scale
            stl.nodes[:, 2] = y * scale
        elif data['--xz']:
            assert data['--xy'] is False
            assert data['--yz']  is False
            x = deepcopy(stl.nodes[:, 0])
            z = deepcopy(stl.nodes[:, 2])
            stl.nodes[:, 0] = z * scale
            stl.nodes[:, 2] = x * scale
    elif data['--xscale'] or data['--yscale'] or data['--zscale']:
        xscale = 1.
        yscale = 1.
        zscale = 1.
        if data['--xscale'] is not None:
            xscale = float(data['--xscale'])
        if data['--yscale'] is not None:
            yscale = float(data['--yscale'])
        if data['--zscale'] is not None:
            zscale = float(data['--zscale'])
        x = deepcopy(stl.nodes[:, 0])
        y = deepcopy(stl.nodes[:, 1])
        z = deepcopy(stl.nodes[:, 2])
        stl.nodes[:, 0] = x * xscale
        stl.nodes[:, 1] = y * yscale
        stl.nodes[:, 2] = z * zscale
    elif data['--xshift'] or data['--yshift'] or data['--zshift']:
        xshift = 1.
        yshift = 1.
        zshift = 1.
        if data['--xshift'] is not None:
            xshift = float(data['--xshift'])
        if data['--yshift'] is not None:
            yshift = float(data['--yshift'])
        if data['--zshift'] is not None:
            zshift = float(data['--zshift'])
        x = deepcopy(stl.nodes[:, 0])
        y = deepcopy(stl.nodes[:, 1])
        z = deepcopy(stl.nodes[:, 2])
        stl.nodes[:, 0] = x + xshift
        stl.nodes[:, 1] = y + yshift
        stl.nodes[:, 2] = z + zshift


    stl.write_stl(out_stl_filename, is_binary=False, float_fmt='%6.7f')

if __name__ == '__main__':
    main()