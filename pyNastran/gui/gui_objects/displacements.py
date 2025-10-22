from copy import deepcopy
from typing import Optional, Any

from numpy import zeros
import numpy as np
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon
from pyNastran.femutils.utils import safe_norm


class VectorTable(GuiResultCommon):
    def __init__(self, subcase_id: int, location: str,
                 titles: list[str], headers: list[str],
                 dxyz: np.ndarray,
                 linked_scale_factor: bool,  # xyz, scalar,
                 scales: list[float],
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
        headers : list[str]
            the sidebar word
        titles : list[str]
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
        scales : list[float]
            the table (e.g., deflection, SPC Forces) scale factors
            nominally, this starts as an empty list and is filled later
        data_formats : list[str]
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
        self.arrow_scales = scales
        self.subcase_id = subcase_id
        self.data_type = self.dxyz.dtype.str  # '<c8', '<f4'
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
            normi = safe_norm(self.dxyz, axis=1)
            self.default_mins[0] = normi.min().real
            self.default_maxs[0] = normi.max().real
        elif self.dim == 3:
            ntimes = self.dxyz.shape[0]
            self.default_mins = zeros(ntimes)
            self.default_maxs = zeros(ntimes)
            for itime in range(ntimes):
                normi = safe_norm(self.dxyz[itime, :, :], axis=1)
                self.default_mins[itime] = normi.min().real
                self.default_maxs[itime] = normi.max().real

            if not self.is_real:
                #: stored in degrees
                self.phases = np.zeros(ntimes)
        else:  # pragma: no cover
            raise NotImplementedError('dim=%s' % self.dim)

        if set_max_min:
            self.min_values = deepcopy(self.default_mins)
            self.max_values = deepcopy(self.default_maxs)
        else:
            self.max_values = None
            self.min_values = None

    def save_defaults(self) -> None:
        ntitles = len(self.titles)
        self.data_formats = ['%g'] * ntitles
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

    def validate(self) -> None:
        ntitles = len(self.titles)
        nheaders = len(self.headers)
        ndata_formats = len(self.data_formats)
        nmin_values = len(self.min_values)
        nmax_values = len(self.max_values)
        assert ntitles > 0, self.titles
        assert nheaders > 0, self.headers
        assert ndata_formats == ntitles, self.data_formats
        for title in self.titles:
            assert isinstance(title, str), title
        for header in self.headers:
            assert isinstance(header, str), header
        for data_format in self.data_formats:
            assert isinstance(data_format, str), data_format
        return

    # getters
    def has_coord_transform(self, i: int,
                            resname: str) -> tuple[bool, list[str]]:
        return True, ['Global']

    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        return False, {}

    def has_nodal_combine_transform(self, i: int,
                                    resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []

    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        is_enabled_fringe = True
        is_checked_fringe = True
        is_enabled_disp = True
        is_checked_disp = True
        is_enabled_vector = False
        is_checked_vector = False
        out = (
            is_enabled_fringe, is_checked_fringe,
            is_enabled_disp, is_checked_disp,
            is_enabled_vector, is_checked_vector)
        return out

    def get_location(self, i: int, unused_name: str) -> str:
        return self.location

    def get_annotation(self, i: int, unused_name: str) -> str:
        #j = self.titles_default.index(name)
        #return self.titles[j]
        return self.headers[i]

    def get_phase(self, i: int, name: str) -> Optional[float]:
        if self.is_real:
            return None
        return self.phases[i]

    def get_data_format(self, i: int, name: str) -> str:
        return self.data_formats[i]

    def get_scale(self, i: int, name: str) -> float:
        return self.scales[i]

    def get_arrow_scale(self, i: int, name: str) -> float:
        return self.arrow_scales[i]

    def get_legend_title(self, i: int, name: str) -> str:
        return self.titles[i]

    def get_min_max(self, i: int, name: str) -> tuple[float, float]:
        return self.min_values[i], self.max_values[i]

    def get_imin_imax(self, i: int, name: str) -> tuple[None, None]:
        return None, None

    #-------------------------------------
    # setters

    def set_data_format(self, i: int, name: str, data_format: str) -> None:
        self.data_formats[i] = data_format

    def set_scale(self, i: int, name: str, scale: float) -> None:
        #j = self.titles_default.index(name)
        if self.linked_scale_factor:
            self.scales[:] = scale
        else:
            self.scales[i] = scale

    def set_arrow_scale(self, i: int, name: str, scale: float) -> None:
        #j = self.titles_default.index(name)
        if self.linked_scale_factor:
            self.arrow_scales[:] = scale
        else:
            self.arrow_scales[i] = scale

    def set_phase(self, i: int, name: str, phase: float) -> None:
        if self.is_real:
            return
        #j = self.titles_default.index(name)
        self.phases[i] = phase

    def set_legend_title(self, i: int, name: str, title: str) -> None:
        self.titles[i] = title

    def set_min_max(self, i: int, name: str,
                    min_value: float, max_value: float) -> None:
        self.min_values[i] = min_value
        self.max_values[i] = max_value

    #-------------------------------------
    # default getters
    def get_default_data_format(self, i: int, name: str):
        return self.data_formats_default[i]

    def get_default_min_max(self, i: int, name: str) -> tuple[float, float]:
        return self.default_mins[i], self.default_maxs[i]

    def get_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def set_nlabels_labelsize_ncolors_colormap(self, i: int, name: str,
                                               nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #def get_default_min_max(self, i: int, name: str) -> tuple[float, float]:
        #return self.min_default[i], self.max_default[i]

    def get_default_scale(self, i: int, name: str) -> float:
        return self.scales_default[i]

    def get_default_arrow_scale(self, i: int, name: str) -> float:
        return self.scales_default[i]

    def get_default_phase(self, i: int, name: str) -> Optional[float]:
        if self.is_real:
            return None
        return 0.0

    def get_default_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):
        # TODO: do this right
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    def get_default_legend_title(self, i: int, resname: str) -> str:
        return self.titles_default[i]

    #-------------------------------------
    # unmodifyable getters

    def get_data_type(self, unused_i, resname: str) -> str:
        """the precision of the data"""
        return self.data_type

    def get_vector_size(self, i: int, resname: str) -> int:
        """vector_size=1 is the default and displacement has 3 components"""
        #print(i)
        #j = self.titles_default.index(name)
        return 3

    #def get_plot_value(self, i:int, resname: str) -> np.ndarray:
        #"""plot returns the displacement..."""
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

    def _get_complex_displacements_by_phase(self, i: int, name: str,
                                            phase: float=0.) -> np.ndarray:
        """
        Get displacements for a complex eigenvector result.
        """
        theta = np.radians(phase)
        dxyz = self.dxyz[i, :].real * np.cos(theta) + self.dxyz[i, :].imag * np.sin(theta)
        return dxyz

    def _get_complex_displacements(self, i: int, name: str) -> np.ndarray:
        """see ``_get_complex_displacements_by_phase``"""
        dxyz = self._get_complex_displacements_by_phase(i, name, self.phases[i])
        return dxyz

    def get_fringe_result(self, i: int, name: str) -> np.ndarray:
        fringe, vector = self.get_fringe_vector_result(i, name)
        return fringe

    def get_fringe_vector_result(self, i: int, name: str) -> tuple[np.ndarray, np.ndarray]:
        """gets the 'typical' result which is a vector"""
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
            dxyz = self._get_complex_displacements(i, name)

        assert len(dxyz.shape) == 2, dxyz.shape
        normi = safe_norm(dxyz, axis=1)
        return normi, dxyz

    def get_force_vector_result(self, i: int, name: str,
                                ) -> tuple[np.ndarray, np.ndarray]:
        return self.get_vector_result(i, name)

    def get_vector_result(self, i: int, name: str,
                          ) -> tuple[np.ndarray, np.ndarray]:
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

        msg = 'VectorTable:\n'
        msg += f'    title={self.titles!r}\n'
        msg += f'    subcase_id={self.subcase_id}\n'
        msg += f'    data_type={self.data_type!r}\n'
        msg += f'    is_real={self.is_real} is_complex={self.is_complex}\n'
        msg += f'    location={self.location!r}\n'
        msg += f'    headers={self.headers!r}\n'
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

    def get_methods(self, i: int) -> list[str]:
        if self.is_real:
            return ['magnitude', 'x', 'y', 'z']
        raise NotImplementedError('self.is_real=%s' % self.is_real)

    def get_vector_result_by_scale_phase(self, i: int, name: str,
                                         unused_scale: float,
                                         phase: float=0.,
                                         ) -> tuple[Optional[np.ndarray], np.ndarray]:
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
            deflected_xyz = self._get_complex_displacements_by_phase(
                i, name, phase)
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return xyz, deflected_xyz


forces = ['Magnitude', 'Fx', 'Fy', 'Fz']
moments = ['Magnitude', 'Mx', 'My', 'Mz']


class ForceTableResults(VectorTable):
    def __init__(self, subcase_id: int, titles: list[str], headers: list[str],
                 dxyz: Any,
                 unused_scalar: Any,
                 scales: list[float],
                 data_formats: Optional[str]=None,
                 nlabels: Optional[int]=None,
                 labelsize: Optional[int]=None,
                 ncolors: Optional[int]=None,
                 colormap: str='jet',
                 set_max_min: bool=False,
                 uname: str='ForceTableResults'):
        """this is a nodal force result"""
        linked_scale_factor = False
        location = 'node'
        VectorTable.__init__(
            self, subcase_id, location, titles, headers, dxyz, linked_scale_factor, scales,
            data_formats=data_formats, nlabels=nlabels,
            labelsize=labelsize, ncolors=ncolors,
            colormap=colormap, set_max_min=set_max_min,
            uname=uname)

    def get_location(self, i: int, name):
        """the result type"""
        return 'node'

    def get_methods(self, i: int, name) -> list[str]:
        return ['Magnitude']

    def get_vector_result_by_scale_phase(self, i: int, name: str,
                                         unused_scale: float,
                                         phase: float=0.,
                                         ) -> tuple[Optional[np.ndarray], np.ndarray]:
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
            deflected_xyz = self._get_complex_displacements_by_phase(
                i, name, phase)
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return xyz, deflected_xyz

    def __repr__(self) -> str:
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


translation = ['Magnitude', 'tx', 'ty', 'tz']
rotation = ['Magnitude', 'rx', 'ry', 'rz']


class DisplacementResults(VectorTable):
    def __init__(self, subcase_id: int,
                 titles: list[str],
                 headers: list[str],
                 xyz: np.ndarray,
                 dxyz: np.ndarray,
                 unused_scalar: np.ndarray,
                 scales: np.ndarray,
                 data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None,
                 colormap: str='jet',
                 set_max_min: bool=False,
                 uname: str='DisplacementResults'):
        """
        Defines a Displacement/Eigenvector result

        Parameters
        ----------
        subcase_id : int
            the flag that points to self.subcases for a message
        headers : list[str]
            the sidebar word
        titles : list[str]
            the legend title
        xyz : (nnodes, 3)
            the nominal xyz locations
        dxyz : (nnodes, 3)
            the delta xyz values
        scalars : (nnodes,n) float ndarray
            #the data to make a contour plot with
            does nothing
        scales : list[float]
            the deflection scale factors
            nominally, this starts as an empty list and is filled later
        data_formats : list[str]
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

    def deflects(self, unused_i: int, unused_res_name: str) -> bool:
        """deflection is opt-in"""
        return True

    def get_location(self, unused_i: int, unused_name: str) -> str:
        """the result type"""
        return 'node'

    #-------------------------------------

    def get_methods(self, i: int, name: str) -> list[str]:
        return ['Magnitude']

    #@property
    #def scalar(self):
        #return self.dxyz_norm

    #def get_scalar(self, i:int, name):
        #print(self.dxyz_norm)
        #return self.dxyz_norm

    def get_vector_result_by_scale_phase(self, i: int, name: str,
                                         scale: float,
                                         phase: float=0.) -> tuple[np.ndarray, np.ndarray]:
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
            dxyz = self._get_complex_displacements_by_phase(
                i, name, phase)
            deflected_xyz = self.xyz + scale * dxyz
        assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        return self.xyz, deflected_xyz

    def __repr__(self) -> str:
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
