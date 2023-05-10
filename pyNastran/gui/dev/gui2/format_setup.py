from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.gui.dev.gui2 import IS_TESTING, IS_OFFICIAL_RELEASE
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.gui.dev.gui2.gui2 import MainWindow2

CLASS_MAP = {}
try:
    from pyNastran.converters.cart3d.cart3d_io import Cart3dIO
    CLASS_MAP['cart3d'] = Cart3dIO
except ImportError:  # pragma: no cover
    pass

try:
    from pyNastran.converters.stl.stl_io import STL_IO
    CLASS_MAP['stl'] = STL_IO
except ImportError:  # pragma: no cover
    pass

def build_fmts(gui: MainWindow2,
               format_class_map,
               fmt_order: list[str],
               log: SimpleLogger,
               stop_on_failure: bool=False,
               ) -> Any:
    """populates the formats that will be supported"""
    fmts = []
    for fmt in fmt_order:
        geom_results_funcs = 'get_%s_wildcard_geometry_results_functions' % fmt

        if fmt in format_class_map:
            cls = format_class_map[fmt](gui)
            data = getattr(cls, geom_results_funcs)()
        #elif hasattr(self, geom_results_funcs):
            #data = getattr(self, geom_results_funcs)()
        else:
            msg = 'get_%s_wildcard_geometry_results_functions does not exist' % fmt
            if stop_on_failure:
                raise RuntimeError(msg)
            if not IS_OFFICIAL_RELEASE:
                if log is None:
                    print('***', msg)
                else:
                    gui.log_error(msg)
        _add_fmt(fmts, fmt, geom_results_funcs, data)

    if len(fmts) == 0:
        RuntimeError('No formats...expected=%s' % fmt_order)
    #self.fmts = fmts
    #print("fmts =", fmts)
    supported_formats = [fmt[0] for fmt in fmts]
    if not IS_TESTING:  # pragma: no cover
        print('supported_formats = %s' % supported_formats)
    #assert 'cart3d' in self.supported_formats, self.supported_formats
    if len(fmts) == 0:
        print('supported_formats = %s' % supported_formats)
        raise RuntimeError('no modules were loaded...')
    return fmts, supported_formats

def _add_fmt(fmts: list[str], fmt: str, geom_results_funcs, data):
    """
    Adds a format

    Parameters
    ----------
    fmts : list[formats]
        format : list[fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func]
        macro_name : ???
            ???
        geo_fmt : ???
            ???
        geo_func : ???
            ???
        res_fmt : ???
            ???
        res_func : ???
            ???
    fmt : str
        nastran, cart3d, etc.
    geom_results_funcs : str
        'get_nastran_wildcard_geometry_results_functions'
        'get_cart3d_wildcard_geometry_results_functions'
    data : function
        the outputs from ``get_nastran_wildcard_geometry_results_functions()``
        so 1 or more formats (macro_name, geo_fmt, geo_func, res_fmt, res_func)

    """
    msg = 'macro_name, geo_fmt, geo_func, res_fmt, res_func = data\n'
    msg += 'data = %s'
    if isinstance(data, tuple):
        assert len(data) == 5, msg % str(data)
        macro_name, geo_fmt, geo_func, res_fmt, res_func = data
        fmts.append((fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func))
    elif isinstance(data, list):
        for datai in data:
            assert len(datai) == 5, msg % str(datai)
            macro_name, geo_fmt, geo_func, res_fmt, res_func = datai
            fmts.append((fmt, macro_name, geo_fmt, geo_func, res_fmt, res_func))
    else:
        raise TypeError(data)
