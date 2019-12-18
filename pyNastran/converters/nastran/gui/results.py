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
        headers : List[str]
            the sidebar word
        titles : List[str]
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

    def finalize(self):
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)

    def get_methods(self, i):
        return self.methods

    def deflects(self, unused_i, unused_res_name):
        return False

    def get_default_title(self, i, name):
        """legend title"""
        (itime, ilayer, imethod, unused_header) = name
        return self.methods[imethod]

    def get_title(self, i, name):
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

    def get_scalar(self, i, name):
        return self.get_result(i, name)

    def get_magnitude(self, i, name):
        scalar = self.get_scalar(i, name)  # TODO: update
        mag = scalar
        if scalar.dtype.name in ['complex64']:
            mag = np.sqrt(scalar.real ** 2 + scalar.imag ** 2)
        return mag

    def get_min_max(self, i, name):
        mag = self.get_magnitude(i, name)
        if np.any(np.isfinite(mag)):
            return np.nanmin(mag), np.nanmax(mag)
        return np.nan, np.nan

    def get_default_min_max(self, i, name):
        mag = self.get_magnitude(i, name)
        if np.any(np.isfinite(mag)):
            return np.nanmin(mag), np.nanmax(mag)
        return np.nan, np.nan

    def get_result(self, i, name):
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
        headers : List[str]
            the sidebar word
        titles : List[str]
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

    def get_methods(self, i):
        return self.methods

    #def deflects(self, unused_i, unused_res_name):
        #return False

    def get_default_title(self, i, name):
        """legend title"""
        (itime, imethod, unused_header) = name
        return self.methods[imethod]

    def get_title(self, i, name):
        """legend title"""
        (itime, imethod, unused_header) = name
        return self.methods[imethod]

    def get_header(self, i, name):
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

    def get_scalar(self, i, name):
        return self.get_result(i, name)

    def get_magnitude(self, i, name):
        scalar = self.get_scalar(i, name)  # TODO: update
        mag = scalar
        if scalar.dtype.name in ['complex64']:
            mag = np.sqrt(scalar.real ** 2 + scalar.imag ** 2)
        return mag

    def get_min_max(self, i, name):
        mag = self.get_magnitude(i, name)
        return np.nanmin(mag), np.nanmax(mag)

    def get_default_min_max(self, i, name):
        mag = self.get_magnitude(i, name)
        return np.nanmin(mag), np.nanmax(mag)

    def get_phase(self, i, name):
        if self.is_real:
            return None
        j = self._get_j(i, name)
        return self.phases[j]

    def get_result(self, i, name):
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
