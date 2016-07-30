from numpy import angle, float32
from pyNastran.utils import object_attributes


def writeFloats10E(vals):
    vals2 = []
    is_all_zeros = True
    for v in vals:
        v2 = '%10.3E' % v
        if v2 in (' 0.000E+00', '-0.000E+00'):
            v2 = ' 0.0'
        else:
            is_all_zeros = False
        vals2.append(v2)
    return vals2, is_all_zeros


def writeFloats12E(vals):
    vals2 = []
    is_all_zeros = True
    for v in vals:
        v2 = '%12.5E' % v
        if v2 in (' 0.00000E+00', '-0.00000E+00'):
            v2 = ' 0.0'
        else:
            is_all_zeros = False
        vals2.append(v2)
    return vals2, is_all_zeros


def write_float_12E(val):
    vals2 = []
    v2 = '%12.5E' % val
    if v2 in (' 0.00000E+00', '-0.00000E+00'):
        v2 = ' 0.0'
    return v2


def write_float_13e(val):
    val2 = '%13.6E' % val
    if val2 in (' 0.000000E+00', '-0.000000E+00'):
        val2 = ' 0.0'
    return val2


#def writeFloats13E(vals):
    #vals2 = []
    #is_all_zeros = True
    #for v in vals:
        #v2 = '%13.6E' % v
        #if v2 in (' 0.000000E+00', '-0.000000E+00'):
            #v2 = ' 0.0'
        #else:
            #is_all_zeros = False
        #vals2.append(v2)
    #return vals2, is_all_zeros

def write_floats_13e(vals):
    vals2 = []
    for v in vals:
        v2 = '%13.6E' % v
        if v2 in (' 0.000000E+00', '-0.000000E+00'):
            v2 = ' 0.0'
        vals2.append(v2)
    return vals2


def write_imag_floats_13e(vals, is_mag_phase):
    vals2 = []

    if is_mag_phase:
        for v in vals:
            v2 = '%13.6E' % abs(v)
            if v2 in (' 0.000000E+00', '-0.000000E+00'):
                v2 = ' 0.0'
            vals2.append(v2)

        # phase
        for v in vals:
            v2 = angle(v, deg=True)

            v3 = '%-13.4f' % v2 if v2 >= 0.0 else '%-13.4f' % (v2 + 360.)
            if v3 == '0.0000       ':
                v3 = '   0.0'
            vals2.append(v3)
    else:
        for v in vals:
            v2 = '%13.6E' % v.real
            if v2 in (' 0.000000E+00', '-0.000000E+00'):
                v2 = ' 0.0'
            vals2.append(v2)

        for v in vals:
            v3 = '%13.6E' % v.imag
            if v3 in (' 0.000000E+00', '-0.000000E+00'):
                v3 = ' 0.0'
            vals2.append(v3)
    return vals2


def writeFloats8p4F(vals):
    vals2 = []
    is_all_zeros = True
    for v in vals:
        if v >= 1000.0 or v <= -100.0:
            raise RuntimeError(v)
        v2 = '%8.4f' % v
        if v2 in ('  0.0000', ' -0.0000'):
            v2 = '  0.0   '
        else:
            is_all_zeros = False
        vals2.append(v2)
    return vals2, is_all_zeros


def _eigenvalue_header(obj, header, itime, ntimes, dt):
    if obj.nonlinear_factor is not None:
        name = obj.data_code['name']
        if isinstance(dt, int):
            dt_line = ' %14s = %i\n' % (name.upper(), dt)
        elif isinstance(dt, (float, float32)):
            dt_line = ' %14s = %12.5E\n' % (name, dt)
        else:
            dt_line = ' %14s = %12.5E %12.5Ej\n' % (name, dt.real, dt.imag)
        header[1] = dt_line
        codes = getattr(obj, name + 's')
        if not len(codes) == ntimes:
            msg = '%ss in %s the wrong size; ntimes=%s; %ss=%s\n' % (name,
                obj.__class__.__name__, ntimes, name, codes)
            atts = object_attributes(obj)
            msg += 'names=%s\n' % atts
            msg += 'data_names=%s\n' % obj.data_names
            raise IndexError(msg)

        if hasattr(obj, 'eigr'):
            try:
                eigenvalue_real = obj.eigrs[itime]
            except IndexError:
                msg = 'eigrs[%s] not found; ntimes=%s; eigrs=%s' % (itime, ntimes, obj.eigrs)
                msg += 'names=%s' % object_attributes(obj)
                raise IndexError(msg)
            eigr_line = ' %14s = %12.6E\n' % ('EIGENVALUE', eigenvalue_real)
            header[2] = eigr_line
    return header


def get_key0_compare(adict):
    """Gets the "first" key in a dictionary

    The entry is kind of irrelevant.
    """
    keys = list(adict.keys())
    return keys[0]


def get_key0(adict):
    """Gets the "first" key in a dictionary

    The entry is kind of irrelevant.
    """
    keys = list(adict.keys())
    return keys[0]
