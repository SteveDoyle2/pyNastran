from copy import deepcopy
import numpy as np

from pyNastran.gui.gui_objects.table import Table


class LayeredTableResults(Table):
    def __init__(self, subcase_id, headers, eids, eid_max, scalars,
                 methods,
                 data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 set_max_min=False, uname='LayeredTableResults'):
        """this is a centroidal result

        Parameters
        ----------
        headers : list[str]
            the sidebar word
        titles : list[str]
            the legend title

        """
        location = 'centroid'
        titles = None
        Table.__init__(
            self, subcase_id, location, titles, headers, scalars,
            data_formats=data_formats, nlabels=nlabels,
            labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, set_max_min=set_max_min,
            uname=uname)
        self.methods = methods
        self.eids = eids
        self.eid_max = eid_max
        self.form_names = []

    def finalize(self):
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)

    def get_methods(self, i: int, name: str) -> list[str]:
        return self.methods

    def has_coord_transform(self, i: int, name: str) -> bool:
        return True
    def has_derivation_transform(self, i: int, resname: str) -> bool:  # min/max/avg
        return True
    def has_nodal_combine_transform(self, i: int, resname: str) -> bool:  # elemental -> nodal
        return True

    def deflects(self, unused_i, unused_res_name):
        return False

    def get_default_legend_title(self, i, name):
        """legend title"""
        (itime, ilayer, imethod, unused_header) = name
        return self.methods[imethod]

    def get_legend_title(self, i, name):
        """legend title"""
        (itime, ilayer, imethod, unused_header) = name
        return self.methods[imethod]

    def get_header(self, i, name):
        """a header shows up in the text"""
        (itime, ilayer, imethod, header) = name
        return self.methods[imethod] + ': ' + header

    def get_data_format(self, i, name):
        return '%.3f'  # TODO: update
    def get_default_data_format(self, i, name):
        return '%.3f'  # TODO: update
    def get_scale(self, i, name):
        return None
    def get_default_scale(self, i, name):
        return None

    def get_scalar(self, i, name, method):
        return self.get_result(i, name, method)

    def get_magnitude(self, i, name, method):
        scalar = self.get_scalar(i, name, method)  # TODO: update
        mag = scalar
        if scalar.dtype.name in ['complex64']:
            mag = np.sqrt(scalar.real ** 2 + scalar.imag ** 2)
        return mag

    def get_min_max(self, i, name, method=''):
        mag = self.get_magnitude(i, name, method)
        if np.any(np.isfinite(mag)):
            return np.nanmin(mag), np.nanmax(mag)
        return np.nan, np.nan

    def get_default_min_max(self, i: int, name: str, method: str) -> tuple[float, float]:
        mag = self.get_magnitude(i, name, method)
        if np.any(np.isfinite(mag)):
            return np.nanmin(mag), np.nanmax(mag)
        return np.nan, np.nan

    def get_result(self, i: int, name: str, method: str):
        (itime, ilayer, imethod, unused_header) = name
        scalars = self.scalars[itime, :, ilayer, imethod]

        if len(scalars) == self.eid_max:
            return scalars
        data = np.full(self.eid_max, np.nan, dtype=scalars.dtype)
        #print(f'data.shape={data.shape} eids.shape={self.eids.shape} scalars.shape={scalars.shape}')
        #print(self.methods)
        data[self.eids] = scalars
        return data

    def __repr__(self):
        """defines str(self)"""
        msg = f'LayeredTableResults:\n'
        msg += f'    title={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    methods={self.methods}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg
