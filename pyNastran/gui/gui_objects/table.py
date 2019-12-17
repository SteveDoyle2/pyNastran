from copy import deepcopy
from typing import List, Optional
import numpy as np
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


class Table(GuiResultCommon):
    def __init__(self, subcase_id: int, location: str,
                 titles: List[str], headers: List[str],
                 scalars,
                 data_formats: Optional[str]=None,
                 nlabels: Optional[int]=None,
                 labelsize: Optional[int]=None,
                 ncolors: Optional[int]=None,
                 colormap: str='jet',
                 set_max_min: bool=False,
                 uname: str='Geometry'):
        """
        Defines a Displacement/Eigenvector result

        Parameters
        ----------
        subcase_id : int
            the flag that points to self.subcases for a message
        headers : List[str]
            the sidebar word
        titles : List[str]
            the legend title
        #xyz : (nnodes, 3)
            #the nominal xyz locations
        #scalars : (nnodes,n) float ndarray
            ##the data to make a contour plot with
            #does nothing
        dxyz : (nnodes, 3)
            the delta xyz values
        linked_scale_factor : bool
            is the displacement scale factor linked
            displacements/loads steps should be
            force/eigenvectors should not be
        scales : List[float]
            the table (e.g., deflection, SPC Forces) scale factors
            nominally, this starts as an empty list and is filled later
        data_formats : List[str]
            the type of data result (e.g. '%i', '%.2f', '%.3f')
        ncolors : int; default=None
            sets the default for reverting the legend ncolors
        set_max_min : bool; default=False
            set default_mins and default_maxs

        Unused
        ------
        #deflects : bool; default=True
            #seems to be an unused parameter...
        uname : str
            some unique name for ...
        """
        GuiResultCommon.__init__(self)

        self.subcase_id = subcase_id
        self.location = location
        assert location in ['node', 'centroid'], 'location=%r' % location
        #assert self.subcase_id > 0, self.subcase_id

        #self.dxyz = dxyz
        #self.dim = len(self.dxyz.shape)


        self.uname = uname
        #self.dxyz_norm = norm(dxyz, axis=1)

        #self.deflects = deflects
        self.titles = titles
        self.headers = headers
        self.subcase_id = subcase_id
        self.scalars = scalars
        self.data_type = self.scalars.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #print('self.data_type = %r' % self.data_type)
        self.is_complex = not self.is_real
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

        if data_formats is not None:
            for data_format in data_formats:
                assert '%' in data_format, data_formats

        self.data_formats = data_formats
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)
        self.data_formats_default = deepcopy(self.data_formats)
        ntimes = scalars.shape[0]
        if not self.is_real:
            #: stored in degrees
            self.phases = np.zeros(ntimes)

        if set_max_min:
            self.min_values = deepcopy(self.default_mins)
            self.max_values = deepcopy(self.default_maxs)
        else:
            self.max_values = None
            self.min_values = None

    def save_defaults(self):
        self.data_formats = ['%g'] * len(self.titles)
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)
        self.scales_default = deepcopy(self.scales)
        self.data_formats_default = deepcopy(self.data_formats)

        size = self.default_maxs.size
        assert size == len(self.titles), 'len(maxs)=%s len(titles)=%s' % (
            size, len(self.titles))
        self.default_mins = self.default_mins.reshape(size)
        self.default_maxs = self.default_maxs.reshape(size)

        self.min_values = deepcopy(self.default_mins)
        self.max_values = deepcopy(self.default_maxs)

    # getters

    def get_location(self, i, unused_name):
        return self.location

    def get_header(self, i, unused_name):
        #j = self.titles_default.index(name)
        #return self.titles[j]
        return self.headers[i]

    def get_phase(self, i, name):
        if self.is_real:
            return None
        return self.phases[i]

    def get_data_format(self, i, name):
        try:
            return self.data_formats[i]
        except IndexError:
            print(f'data_formats = {self.data_formats}')
            print(str(self))
            print("ires =", i)
            print(name)
            raise

    def get_scale(self, i, name):
        return self.scales[i]

    def get_title(self, i, name):
        return self.titles[i]

    def get_min_max(self, i, name):
        return self.min_values[i], self.max_values[i]

    #-------------------------------------
    # setters

    def set_data_format(self, i, name, data_format):
        self.data_formats[i] = data_format

    #def set_scale(self, i, name, scale):
        ##j = self.titles_default.index(name)
        #if self.linked_scale_factor:
            #self.scales[:] = scale
        #else:
            #self.scales[i] = scale

    def set_phase(self, i, name, phase):
        if self.is_real:
            return
        #j = self.titles_default.index(name)
        self.phases[i] = phase

    def set_title(self, i, name, title):
        self.titles[i] = title

    def set_min_max(self, i, name, min_value, max_value):
        self.min_values[i] = min_value
        self.max_values[i] = max_value

    #-------------------------------------
    # default getters
    def get_default_data_format(self, i, name):
        (itime, imethod, unused_header) = name
        ntimes = self.scalars.shape[0]
        j = ntimes * imethod + itime
        return self.data_formats_default[j]

    def get_default_min_max(self, i, name):
        return self.default_mins[i], self.default_maxs[i]

    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #def get_default_min_max(self, i, name):
        #return self.min_default[i], self.max_default[i]

    def get_default_scale(self, i, name):
        return self.scales_default[i]

    def get_default_phase(self, i, name):
        if self.is_real:
            return None
        return 0.0

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        # TODO: do this right
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    def get_default_title(self, i, name):
        return self.titles_default[i]

    #-------------------------------------
    # unmodifyable getters

    def get_data_type(self, unused_i, unused_name):
        """the precision of the data"""
        return self.data_type

    #def get_vector_size(self, i, unused_name):
        #"""the result size"""
        #print(i)
        #j = self.titles_default.index(name)
        #return 3

    def get_plot_value(self, i, unused_name):
        asdf
        return self
        #if self.is_real:
            #if self.dim == 2:
                #dxyz = self.dxyz
            #elif self.dim == 3:
                #dxyz = self.dxyz[i, :]
            #else:
                #raise NotImplementedError('dim=%s' % self.dim)
        #else:
            #dxyz = self._get_complex_displacements(i)

        #assert len(dxyz.shape) == 2, dxyz.shape
        #return dxyz

    #def _get_complex_displacements_by_phase(self, i, phase=0.):
        #"""
        #Get displacements for a complex eigenvector result.
        #"""
        #theta = np.radians(phase)
        #dxyz = self.dxyz[i, :].real * np.cos(theta) + self.dxyz[i, :].imag * np.sin(theta)
        #return dxyz

    #def _get_complex_displacements(self, i):
        #"""see ``_get_complex_displacements_by_phase``"""
        #dxyz = self._get_complex_displacements_by_phase(i, self.phases[i])
        #return dxyz

    #def get_result(self, i, name):
        #if self.is_real:
            #if self.dim == 2:
                ## single result
                #dxyz = self.dxyz
            #elif self.dim == 3:
                ## multiple results
                ## .0006 -> 0.0
                ## .057 -> 0.0123
                ## min
                #dxyz = self.dxyz[i, :]
            #else:
                #raise NotImplementedError('dim=%s' % self.dim)
        #else:
            #dxyz = self._get_complex_displacements(i)

        #assert len(dxyz.shape) == 2, dxyz.shape
        #return dxyz

    def get_vector_result(self, i, name):
        bbb
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        #if self.is_real:
            #xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                #i, name, self.scales[i])
        #else:
            #xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                #i, name, self.scales[i], self.phases[i])
        #return xyz, deflected_xyz

    def __repr__(self):
        """defines str(self)"""
        #self.subcase_id = subcase_id
        #self.location = location
        #assert location in ['node', 'centroid'], 'location=%r' % location
        #self.linked_scale_factor = linked_scale_factor

        #self.dxyz = dxyz
        #self.dim = len(self.dxyz.shape)


        #self.uname = uname
        #self.dxyz_norm = norm(dxyz, axis=1)

        #self.deflects = deflects
        #self.titles = titles
        #self.headers = headers
        #self.scales = scales
        #self.subcase_id = subcase_id
        #self.data_type = self.dxyz.dtype.str # '<c8', '<f4'
        #self.nlabels = nlabels
        #self.labelsize = labelsize
        #self.ncolors = ncolors
        #self.colormap = colormap
        #self.data_formats = data_formats
        #self.titles_default = deepcopy(self.titles)
        #self.headers_default = deepcopy(self.headers)
        #self.scales_default = deepcopy(self.scales)
        #self.data_formats_default = deepcopy(self.data_formats)
        ##elif self.dim == 3:
        #self.default_mins[itime] = normi.min().real
        #self.default_maxs[itime] = normi.max().real

        #self.phases = np.zeros(ntimes)
        #self.min_values = deepcopy(self.default_mins)
        #self.max_values = deepcopy(self.default_maxs)

        msg = f'Table ({self.__class__.__name__}):\n'
        msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg
