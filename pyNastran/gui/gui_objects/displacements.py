from copy import deepcopy
from typing import List, Optional, Any

from numpy import zeros
import numpy as np
from numpy.linalg import norm  # type: ignore
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


class VectorTable(GuiResultCommon):
    def __init__(self, subcase_id: int, location: str,
                 titles: List[str], headers: List[str],
                 dxyz: Any, linked_scale_factor: bool, #xyz, scalar,
                 scales,
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
        self.linked_scale_factor = linked_scale_factor
        #assert self.subcase_id > 0, self.subcase_id

        self.dxyz = dxyz
        self.dim = len(self.dxyz.shape)


        self.uname = uname
        #self.dxyz_norm = norm(dxyz, axis=1)

        #self.deflects = deflects
        self.titles = titles
        self.headers = headers
        self.scales = scales
        self.subcase_id = subcase_id
        self.data_type = self.dxyz.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<f4', '<f8'] else False
        #print('self.data_type = %r' % self.data_type)
        self.is_complex = not self.is_real
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

        self.data_formats = data_formats
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)
        self.scales_default = deepcopy(self.scales)
        self.data_formats_default = deepcopy(self.data_formats)
        if self.dim == 2:
            ntimes = 1
            self.default_mins = zeros(1, dtype=self.dxyz.dtype)
            self.default_maxs = zeros(1, dtype=self.dxyz.dtype)
            normi = norm(self.dxyz, axis=1)
            self.default_mins[0] = normi.min().real
            self.default_maxs[0] = normi.max().real
        elif self.dim == 3:
            ntimes = self.dxyz.shape[0]
            self.default_mins = zeros(ntimes)
            self.default_maxs = zeros(ntimes)
            for itime in range(ntimes):
                normi = norm(self.dxyz[itime, :, :], axis=1)
                self.default_mins[itime] = normi.min().real
                self.default_maxs[itime] = normi.max().real

            if not self.is_real:
                #: stored in degrees
                self.phases = np.zeros(ntimes)
        else:
            raise NotImplementedError('dim=%s' % self.dim)

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
        return self.data_formats[i]

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

    def set_scale(self, i, name, scale):
        #j = self.titles_default.index(name)
        if self.linked_scale_factor:
            self.scales[:] = scale
        else:
            self.scales[i] = scale

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
        return self.data_formats_default[i]

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

    def get_vector_size(self, i, unused_name):
        """the result size"""
        #print(i)
        #j = self.titles_default.index(name)
        return 3

    def get_plot_value(self, i, unused_name):
        if self.is_real:
            if self.dim == 2:
                dxyz = self.dxyz
            elif self.dim == 3:
                dxyz = self.dxyz[i, :]
            else:
                raise NotImplementedError('dim=%s' % self.dim)
        else:
            dxyz = self._get_complex_displacements(i)

        assert len(dxyz.shape) == 2, dxyz.shape
        return dxyz

    def _get_complex_displacements_by_phase(self, i, phase=0.):
        """
        Get displacements for a complex eigenvector result.
        """
        theta = np.radians(phase)
        dxyz = self.dxyz[i, :].real * np.cos(theta) + self.dxyz[i, :].imag * np.sin(theta)
        return dxyz

    def _get_complex_displacements(self, i):
        """see ``_get_complex_displacements_by_phase``"""
        dxyz = self._get_complex_displacements_by_phase(i, self.phases[i])
        return dxyz

    def get_result(self, i, name):
        if self.is_real:
            if self.dim == 2:
                # single result
                dxyz = self.dxyz
            elif self.dim == 3:
                # multiple results
                # .0006 -> 0.0
                # .057 -> 0.0123
                # min
                dxyz = self.dxyz[i, :]
            else:
                raise NotImplementedError('dim=%s' % self.dim)
        else:
            dxyz = self._get_complex_displacements(i)

        assert len(dxyz.shape) == 2, dxyz.shape
        return dxyz

    def get_vector_result(self, i, name):
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                i, name, self.scales[i])
        else:
            xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                i, name, self.scales[i], self.phases[i])
        return xyz, deflected_xyz

    def __repr__(self):
        """defines str(self)"""
        #self.subcase_id = subcase_id
        #self.location = location
        #assert location in ['node', 'centroid'], 'location=%r' % location
        #self.linked_scale_factor = linked_scale_factor

        #self.dxyz = dxyz
        #self.dim = len(self.dxyz.shape)


        self.uname = uname
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

        msg = 'VectorTable:\n'
        msg += f'    title={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg


class ElementalTableResults(VectorTable):
    def __init__(self, subcase_id, titles, headers, dxyz, unused_scalar,
                 scales, data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 set_max_min=False, uname='Geometry'):
        """this is a centroidal result"""
        linked_scale_factor = False
        location = 'centroid'
        VectorTable.__init__(
            self, subcase_id, location, titles, headers, dxyz, linked_scale_factor, scales,
            data_formats=data_formats, nlabels=nlabels,
            labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, set_max_min=set_max_min,
            uname=uname)

    def get_methods(self, i):
        if self.is_real:
            return ['magnitude', 'x', 'y', 'z']
        else:
            raise NotImplementedError('self.is_real=%s' % self.is_real)

    def deflects(self, unused_i, unused_res_name):
        return False

    def get_vector_result_by_scale_phase(self, i, unused_name, unused_scale, phase=0.):
        xyz = None
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            if self.dim == 2:
                # single result
                deflected_xyz = self.dxyz
            elif self.dim == 3:
                deflected_xyz = self.dxyz[i, :]
            else:
                raise NotImplementedError('dim=%s' % self.dim)
        else:
            deflected_xyz = self._get_complex_displacements_by_phase(i, phase)
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return xyz, deflected_xyz


class ForceTableResults(VectorTable):
    def __init__(self, subcase_id, titles, headers, dxyz, unused_scalar,
                 scales, data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 set_max_min=False, uname='ForceTableResults'):
        # type: (int, List[str], List[str], Any, Any, List[float], Any, Optional[int], Optional[int], Optional[int], str, bool, str) -> None
        """this is a nodal force result"""
        linked_scale_factor = False
        location = 'node'
        VectorTable.__init__(
            self, subcase_id, location, titles, headers, dxyz, linked_scale_factor, scales,
            data_formats=data_formats, nlabels=nlabels,
            labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, set_max_min=set_max_min,
            uname=uname)

    def get_methods(self, i):
        if self.is_real:
            return ['magnitude', 'tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        else:
            return ['node real', 'node imag', 'node magnitude', 'node phase']

    #def get_location(self, i, name):
        #"""the result type"""
        #return 'node'

    def deflects(self, unused_i, unused_res_name):
        return False

    def get_vector_result_by_scale_phase(self, i, unused_name, unused_scale, phase=0.):
        """
        Gets the real/complex deflection result

        Parameters
        ----------
        i : int
            mode/time/loadstep number
        name : str
            unused; useful for debugging
        scale : float
            deflection scale factor; true scale
            unused because it's a force
        phase : float; default=0.0
            phase angle (degrees); unused for real results

        Returns
        -------
        xyz : (nnodes, 3) float ndarray
            the nominal state
        deflected_xyz : (nnodes, 3) float ndarray
            the deflected state
        """
        xyz = None
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            if self.dim == 2:
                # single result
                deflected_xyz = self.dxyz
            elif self.dim == 3:
                deflected_xyz = self.dxyz[i, :]
            else:
                raise NotImplementedError('dim=%s' % self.dim)
        else:
            deflected_xyz = self._get_complex_displacements_by_phase(i, phase)
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return xyz, deflected_xyz

    def __repr__(self):
        """defines str(self)"""
        msg = 'ForceTableResults\n'
        msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg


class DisplacementResults(VectorTable):
    def __init__(self, subcase_id, titles, headers, xyz, dxyz, unused_scalar,
                 scales, data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 set_max_min=False, uname='DisplacementResults'):
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
        xyz : (nnodes, 3)
            the nominal xyz locations
        dxyz : (nnodes, 3)
            the delta xyz values
        scalars : (nnodes,n) float ndarray
            #the data to make a contour plot with
            does nothing
        scales : List[float]
            the deflection scale factors
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
        linked_scale_factor = False
        location = 'node'
        VectorTable.__init__(
            self, subcase_id, location, titles, headers, dxyz, linked_scale_factor,
            scales,
            data_formats=data_formats, nlabels=nlabels,
            labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, set_max_min=set_max_min,
            uname=uname)
        ##assert self.subcase_id > 0, self.subcase_id
        self.xyz = xyz

    #-------------------------------------
    # unmodifyable getters

    def deflects(self, unused_i, unused_res_name):
        return True

    def get_location(self, unused_i, unused_name):
        """the result type"""
        return 'node'

    #-------------------------------------

    def get_methods(self, i):
        if self.is_real:
            return ['magnitude', 'tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        else:
            return ['node real', 'node imag', 'node magnitude', 'node phase']

    #@property
    #def scalar(self):
        #return self.dxyz_norm

    #def get_scalar(self, i, name):
        #print(self.dxyz_norm)
        #return self.dxyz_norm

    def get_vector_result_by_scale_phase(self, i, unused_name, scale, phase=0.):
        """
        Gets the real/complex deflection result

        Parameters
        ----------
        i : int
            mode/time/loadstep number
        name : str
            unused; useful for debugging
        scale : float
            deflection scale factor; true scale
        phase : float; default=0.0
            phase angle (degrees); unused for real results

        Returns
        -------
        xyz : (nnodes, 3) float ndarray
            the nominal state
        deflected_xyz : (nnodes, 3) float ndarray
            the deflected state
        """
        assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            if self.dim == 2:
                # single result
                deflected_xyz = self.xyz + scale * self.dxyz
            elif self.dim == 3:
                deflected_xyz = self.xyz + scale * self.dxyz[i, :]
            else:
                raise NotImplementedError('dim=%s' % self.dim)
        else:
            assert isinstance(i, int), (i, phase)
            assert isinstance(phase, float), (i, phase)
            dxyz = self._get_complex_displacements_by_phase(i, phase)
            deflected_xyz = self.xyz + scale * dxyz
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return self.xyz, deflected_xyz

    def __repr__(self):
        """defines str(self)"""
        msg = 'DisplacementResults\n'
        msg += f'    titles={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    header={self.headers!r}\n'
        msg += f'    data_format={self.data_formats!r}\n'
        msg += f'    uname={self.uname!r}\n'
        return msg
