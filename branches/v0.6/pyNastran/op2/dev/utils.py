from struct import calcsize

def iformat(Format, precision=2):
    if precision == 2:  # double
        f = Format.replace('i', 'l').replace('f', 'd')
    elif precision == 1: # single
        f = Format.replace('l', 'i').replace('d', 'f')
    else:
        raise NotImplementedError(precision)

    ndata = calcsize(f)
    return f, ndata

print iformat('4si3f', 2)
#print iformat('4si3f')