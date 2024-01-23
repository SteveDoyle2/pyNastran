from copy import deepcopy
import numpy as np

from pyNastran.gui.gui_objects.table import Table


class SimpleTableResults(Table):
    def __init__(self, subcase_id, headers, eids, eid_max, scalars,
                 methods,
                 data_format=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 location='centroid',
                 set_max_min=False, uname='Geometry'):
        """this is a centroidal result

        Parameters
        ----------
        headers : list[str]
            the sidebar word
        titles : list[str]
            the legend title

        """
        titles = None
        assert data_format is not None, data_format
        ntimes = scalars.shape[0]
        data_formats = [data_format] * len(methods) * ntimes
        Table.__init__(
            self, subcase_id, location, titles, headers, scalars,
            data_formats=data_formats, nlabels=nlabels,
            labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, set_max_min=set_max_min,
            uname=uname)
        self.methods = methods
        self.eids = eids
        self.eid_max = eid_max
        assert len(eids) == scalars.shape[1], f'len(eids)={len(eids)} scalars.shape={scalars.shape}'
        if self.is_complex:
            self.phases = np.zeros(len(methods) * ntimes)

    def finalize(self):
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)

    def get_methods(self, i: int, resname: str):
        return self.methods

    def get_default_legend_title(self, i, name):
        """legend title"""
        (itime, imethod, unused_header) = name
        return self.methods[imethod]

    def get_legend_title(self, i, name):
        """legend title"""
        (itime, imethod, unused_header) = name
        return self.methods[imethod]

    def get_annotation(self, i, name):
        """a header shows up in the text"""
        (itime, imethod, header) = name
        return self.methods[imethod] + ': ' + header

    #def get_data_format(self, i, name):
        #return '%.3f'  # TODO: update
    #def get_default_data_format(self, i, name):
        #return '%.3f'  # TODO: update
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
        return np.nanmin(mag), np.nanmax(mag)

    def get_default_min_max(self, i, name) -> tuple[float, float]:
        mag = self.get_magnitude(i, name, method='')
        return np.nanmin(mag), np.nanmax(mag)

    def get_phase(self, i, name):
        if self.is_real:
            return None
        j = self._get_j(i, name)
        return self.phases[j]

    def get_result(self, i, name, method: str):
        #print(i, name)
        (itime, imethod, unused_header) = name
        scalars = self.scalars[itime, :, imethod]

        if len(scalars) == self.eid_max:
            return scalars
        data = np.full(self.eid_max, np.nan, dtype=scalars.dtype)
        #print(f'data.shape={data.shape} eids.shape={self.eids.shape} scalars.shape={scalars.shape}')
        #print(self.methods)
        try:
            data[self.eids] = scalars
        except IndexError:
            raise RuntimeError(f'{self.uname!r} eids.max()={self.eids.max()} scalars.shape={scalars.shape}')
        return data

    def _get_j(self, i, name):
        (itime, imethod, unused_header) = name
        ntimes = self.scalars.shape[0]
        j = ntimes * imethod + itime
        return j

    def has_coord_transform(self, i: int, name: str) -> tuple[bool, list[str]]:
        return True, ['Material']
    def has_derivation_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """min/max/avg"""
        out = {'derivation': ['Absolute Max'], }
        return True, out
    def has_nodal_combine_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return True, ['Centroid']
    #def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool]:
        #is_enabled_fringe = False
        #is_checked_fringe = True
        #is_enabled_disp = False
        #is_checked_disp = False
        #is_enabled_vector = False
        #is_checked_vector = False
        #out = (
            #is_enabled_fringe, is_checked_fringe,
            #is_enabled_disp, is_checked_disp,
            #is_enabled_vector, is_checked_vector)
        #return out

    def get_data_format(self, i, name):
        j = self._get_j(i, name)
        try:
            return self.data_formats[j]
        except IndexError:
            print(f'data_formats = {self.data_formats}')
            print(str(self))
            print("ires =", i)
            print(name)
            raise

    def __repr__(self):
        """defines str(self)"""
        msg = f'SimpleTableResults:\n'
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
