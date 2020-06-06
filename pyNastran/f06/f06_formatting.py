from typing import List
import numpy as np
from pyNastran.utils import object_attributes


def write_float_11e(val: float) -> str:
    """writes a Nastran formatted 11.4 float"""
    v2 = '%11.4E' % val
    if v2 in (' 0.0000E+00', '-0.0000E+00'):
        v2 = ' 0.0'
    return v2

def write_float_12e(val: float) -> str:
    """writes a Nastran formatted 12.5 float"""
    v2 = '%12.5E' % val
    if v2 in (' 0.00000E+00', '-0.00000E+00'):
        v2 = ' 0.0'
    return v2


def write_float_13e(val: float) -> str:
    """writes a Nastran formatted 13.6 float"""
    val2 = '%13.6E' % val
    if val2 in (' 0.000000E+00', '-0.000000E+00'):
        val2 = ' 0.0'
    return val2


def write_floats_10e(vals: List[float]) -> List[str]:
    """writes a series of Nastran formatted 10.3 floats"""
    vals2 = []
    for v in vals:
        v2 = '%10.3E' % v
        if v2 in (' 0.000E+00', '-0.000E+00'):
            v2 = ' 0.0'
        vals2.append(v2)
    return vals2


def write_floats_12e(vals: List[float]) -> List[str]:
    """writes a series of Nastran formatted 12.5 floats"""
    vals2 = []
    for v in vals:
        v2 = '%12.5E' % v
        if v2 in (' 0.00000E+00', '-0.00000E+00'):
            v2 = ' 0.0'
        vals2.append(v2)
    return vals2


def write_floats_13e(vals: List[float]) -> List[str]:
    """writes a series of Nastran formatted 13.6 floats"""
    vals2 = []
    for v in vals:
        v2 = '%13.6E' % v
        if v2 in (' 0.000000E+00', '-0.000000E+00'):
            v2 = ' 0.0'
        vals2.append(v2)
    return vals2


def write_imag_floats_13e(vals: List[float], is_mag_phase: bool) -> str:
    vals2 = []

    if is_mag_phase:
        for v in vals:
            v2 = '%13.6E' % abs(v)
            if v2 in (' 0.000000E+00', '-0.000000E+00'):
                v2 = ' 0.0'
            vals2.append(v2)

        # phase
        for v in vals:
            v2 = np.angle(v, deg=True)

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


def write_floats_8p4f(vals: List[float]) -> List[str]:
    """writes an 8.4F formatted number"""
    vals2 = []
    for val in vals:
        if val >= 1000.0 or val <= -100.0:
            raise RuntimeError(val)
        val2 = '%8.4f' % val
        if val2 in ('  0.0000', ' -0.0000'):
            val2 = '  0.0   '
        vals2.append(val2)
    return vals2

def write_floats_8p1e(vals: List[float]) -> List[str]:
    """writes an 8.1E formatted number"""
    vals2 = []
    for val in vals:
        if abs(val) < 1e-10:
            vals2.append('       ')
        else:
            vals2.append('%8.1E' % val)
    return vals2

def _eigenvalue_header(obj, header, itime: int, ntimes: int, dt):
    if obj.nonlinear_factor not in (None, np.nan):
        name = obj.data_code['name']
        if isinstance(dt, (int, np.int32, np.int64)):
            dt_line = ' %14s = %i\n' % (name.upper(), dt)
        elif isinstance(dt, (float, np.float32, np.float64)):
            dt_line = ' %14s = %12.5E\n' % (name, dt)
        else:
            dt_line = ' %14s = %12.5E %12.5Ej\n' % (name, dt.real, dt.imag)
        header[1] = dt_line
        codes = getattr(obj, name + 's')
        if not len(codes) == ntimes:
            msg = (f'{name}s in {obj.__class__.__name__} the wrong size; '
                   f'ntimes={ntimes}; {name}s={codes}\n')
            atts = object_attributes(obj)
            msg += f'names={atts}\n'
            msg += f'data_names={obj.data_names}\n'
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
