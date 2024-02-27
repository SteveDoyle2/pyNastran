# -*- coding: utf-8 -*-
from scipy.interpolate import splrep, splev  # type: ignore
from scipy.integrate import quad  # type: ignore


def build_spline(x, y):
    """
    Builds a cubic spline or 1st order spline if there are less than 3 terms

    Parameters
    ----------
    x : list[float]
        the independent variable
    y : list[float]
        the dependent variable

    Returns
    -------
    splrep : splrep object
        linear or cubic spline depending on the length of x

    .. note:: a 1st order spline is the same as linear interpolation

    """
    #return splrep(x, y, k=1) if len(x) < 3 else splrep(x, y)

    if len(x) == 2:
        # build a linearly interpolated spline
        k = 1
    elif len(x) == 3:
        # build a quadratic spline
        k = 2
    else:
        # build a cubic spline
        k = 3
    return splrep(x, y, k=k)


def integrate_positive_unit_line(x, y, min_value=0.):
    """
    Integrates a line of length 1.0 by spline interpolation

    Parameters
    ----------
    x : list[float]
        the independent variable
    y : list[float]
        the dependent variable
    min_value : float; default=0.0
        ???

    Returns
    -------
    integrated_value : float
        the area under the curve
    """
    if len(set(y)) == 1:
        return y[0]  # (x1-x0 = 1., so yBar*1 = yBar)
    try:
        assert len(x) == len(y), 'x=%s y=%s' % (x, y)
        # now integrate the area
        eval_posit_spline = lambda x, spl, min_val: max(splev([x], spl), min_val)
        out = quad(eval_posit_spline, 0., 1., args=(build_spline(x, y), min_value))
    except Exception:
        raise RuntimeError('spline Error x=%s y=%s' % (x, y))
    return out[0]
