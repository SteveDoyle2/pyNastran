mapper = {
    's' : 1,
    'i' : 4,
    'f' : 4,

    'l' : 8,
    'd' : 8,
}

import struct

def iformat(Format, precision=2):
    if precision == 2:  # double
        f = Format.replace('i', 'l').replace('f', 'd')
    elif precision == 1: # single
        f = Format.replace('l', 'i').replace('d', 'f')
    else:
        raise NotImplementedError(precision)


    packs = []
    pack = ['', '']
    for fi in f:
        if fi.isdigit():
            pack[0] += fi
        else:
            pack[1] = fi
            packs.append(pack)
            pack = ['', '']
    if pack != ['', '']:
        packs.append(pack)

    f2 = ''
    ndata = 0
    for pack in packs:
        if pack[0] == '':
            pack[0] = 1
        else:
            pack[0] = int(pack[0])

        ndata += pack[0] * mapper[pack[1]]
        f2 += str(pack[0]) + pack[1]
    print packs
    return f2, ndata

print iformat('4si3f', 2)
#print iformat('4si3f')