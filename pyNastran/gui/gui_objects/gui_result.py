"""
defines:
 - GuiResultCommon
 - GuiResult

"""
from abc import abstractmethod
from typing import Any, Optional
import numpy as np
in1d = np.in1d
from pyNastran.utils.numpy_utils import integer_float_types

REAL_TYPES = ['<i4', '<i8', '<f4', '<f8',
              '|i1', # this is a boolean
              #'|b1', # this is officially a boolean...
              '>i4', '>i8', '>f4', '>f8']
#COMPLEX_TYPES = ['<c8']
INT_TYPES = ['<i4', '<i8', '|i1',
             '>i4', '>i8']

class GuiResultCommon:
    def __init__(self):
        self.class_name = self.__class__.__name__
        self.is_real = False
        self.is_complex = False
        self.has_methods_table(0, 'test')
        self.has_coord_transform(0, 'test')
        self.has_derivation_transform(0, 'test')
        self.has_nodal_combine_transform(0, 'test')
        self.has_output_checks(0, 'test')

    #def get_data_type(self, i: int, name: str):  # pragma: no cover
        #raise NotImplementedError(f'{self.class_name}.get_data_type')

    #def get_annotation(self, i: int, name: str):  # pragma: no cover
        #raise NotImplementedError(f'{self.class_name}.get_annotation')

    def set_sidebar_args(self, i: int, name: str, **kwargs) -> None:
        pass
    def has_methods_table(self, i: int, res_name: str) -> bool:
        return False
    def has_coord_transform(self, i: int, name: str) -> tuple[bool, list[str]]:  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.has_coord_transform')
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:  # pragma: no cover
        """min/max/avg"""
        raise NotImplementedError(f'{self.class_name}.has_derivation_transform')
    def has_nodal_combine_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        raise NotImplementedError(f'{self.class_name}.has_nodal_combine_transform')
    def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool,
                                                               bool, bool, bool]:
        is_enabled_fringe = False
        is_checked_fringe = True
        is_enabled_disp = False
        is_checked_disp = False
        is_enabled_vector = False
        is_checked_vector = False
        out = (
            is_enabled_fringe, is_checked_fringe,
            is_enabled_disp, is_checked_disp,
            is_enabled_vector, is_checked_vector)
        return out

    def deflects(self, unused_i: int, unused_res_name: str) -> bool:
        """deflection is opt-in"""
        return False
    def is_normal_result(self, i: int, name: str) -> bool:
        """normal result is opt-in"""
        return False

    def get_data_format(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_data_format')

    def get_location(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_location')

    def get_legend_title(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_legend_title')

    def get_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_nlabels_labelsize_ncolors_colormap')

    @abstractmethod
    def get_min_max(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_min_max')

    #@abstractmethod
    #def get_scalar(self, i: int, name: str):  # pragma: no cover
        #raise NotImplementedError(f'{self.class_name}.get_scalar')

    @abstractmethod
    def get_methods(self, i: int, name: str) -> list[str]:  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_methods')

    @abstractmethod
    def get_fringe_vector_result(self, i: int, name: str) -> tuple[Any, Any]:  # pragma: no cover
        raise NotImplementedError(f'{self.class_name}.get_fringe_vector_result')

    def get_vector_size(self, i: int, name: str) -> int:
        """vector_size=1 is the default"""
        return 1

    def get_scale(self, i: int, name: str) -> float:
        return 0.0

    def get_phase(self, i: int, name: str) -> Optional[float]:
        return None

    #------------
    # setters

    def set_data_format(self, i: int, name: str, data_format):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def set_min_max(self, i: int, name: str, min_value, max_value):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def set_legend_title(self, i: int, name: str, title):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def set_nlabels_labelsize_ncolors_colormap(self, i: int, name: str, nlabels, labelsize,
                                               ncolors, colormap):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def set_scale(self, i: int, name: str, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')

    def set_phase(self, i: int, name: str, phase) -> None:
        pass

    #------------
    # default getters
    def get_default_data_format(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def get_default_min_max(self, i: int, name: str) -> tuple[float, float]:  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def get_default_legend_title(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def get_default_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):  # pragma: no cover
        raise NotImplementedError(self.class_name)

    def get_default_scale(self, i: int, name: str):
        return 0.

    def get_default_phase(self, i: int, name: str):
        return None

    def _get_complex_displacements_by_phase(self, i: int, phase):
        raise NotImplementedError(self.class_name)


class NullResult(GuiResultCommon):
    def __init__(self):
        super(NullResult, self).__init__()

    #def get_scalar(self, i: int, name: str) -> None:
        #return None

    def __repr__(self) -> str:
        msg = '<NormalResult>'
        return msg

class GridPointForceResult(GuiResultCommon):
    def __init__(self, subcase_id: int,
                 header: str, title: str,
                 gpforce_array,
                 uname: str='GridPointForceResult'):
        """
        Parameters
        ----------
        subcase_id : int
            the flag that points to self.subcases for a message
        header : str
            the sidebar word
        title : str
            the legend title
        uname : str
            some unique name for ...
        """
        GuiResultCommon.__init__(self)

        self.subcase_id = subcase_id
        self.title = title
        self.header = header
        #self.scale = scale
        #self.location = location
        self.location = 'node'
        self.subcase_id = subcase_id
        self.uname = uname
        super(GridPointForceResult, self).__init__()
        self.gpforce_array = gpforce_array

    def has_coord_transform(self, i: int, name: str) -> tuple[bool, list[str]]:
        return False, []
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        return False, {}
    def has_nodal_combine_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []

    #def get_scalar(self, i: int, name: str) -> None:
        #return None

    def get_fringe_vector_result(self, i: int, name: str) -> tuple[None, None]:
        return None, None
    def get_legend_title(self, i: int, name: str) -> str:
        return self.title
    def get_location(self, i: int, name: str) -> str:
        return self.location
    def get_annotation(self, i: int, name: str) -> str:
        return self.header
    def get_methods(self, i: int, name: str) -> list[None]:
        return [None]
    def get_data_format(self, i: int, name: str) -> None:
        return None
    def get_nlabels_labelsize_ncolors_colormap(self, i: int, name: str) -> tuple[Any, Any, Any, Any]:
        return None, None, None, None
    def get_default_data_format(self, i: int, name: str) -> None:
        return None
    def get_default_min_max(self, i: int, name: str) -> tuple[Optional[float], Optional[float]]:
        return None, None
    def get_default_legend_title(self, i: int, name: str) -> str:
        return self.title
    def get_default_nlabels_labelsize_ncolors_colormap(self, i: int, name: str) -> tuple[Any, Any, Any, Any]:
        return None, None, None, None
    def set_nlabels_labelsize_ncolors_colormap(self, i: int, name: str,
                                               nlabels, labelsize, ncolors, colormap) -> None:
        return
    def get_min_max(self, i: int, name: str) -> tuple[Any, Any]:
        return None, None

    #def save_vtk_result(self, used_titles: set[str]) -> None:

    def __repr__(self) -> str:
        msg = '<GridPointForceResult>'
        return msg

class NormalResult(GuiResultCommon):
    def __init__(self, subcase_id: int,
                 header: str, title: str,
                 nlabels=2, labelsize=5, ncolors=2,
                 colormap: str='jet',
                 data_format: str='%.1f',
                 uname: str='NormalResult'):
        """
        subcase_id : int
            the flag that points to self.subcases for a message
        header : str
            the sidebar word
        title : str
            the legend title
        #location : str
            #node, centroid
        #scalar : (n,) ndarray
            #the data to make a contour plot with
        data_format : str
            the type of data result (e.g. '%i', '%.2f', '%.3f')
        uname : str
            some unique name for ...
        """
        GuiResultCommon.__init__(self)
        self.scalar = None
        self.subcase_id = subcase_id
        self.title = title
        self.header = header
        self.data_format = data_format

        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

        self.title_default = self.title
        self.header_default = self.header
        self.data_format_default = self.data_format

        self.min_value = -1.
        self.max_value = 1.
        self.min_default = -1.
        self.max_default = 1.
        self.uname = uname

    def has_coord_transform(self, i: int, name: str) -> tuple[bool, list[str]]:
        return False, []
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        return False, {}
    def has_nodal_combine_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []

    def get_data_type(self, i: int, name: str):
        #print('Aname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_type

    def get_data_format(self, i: int, name: str):
        #print('Bname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_format

    def get_location(self, i: int, name: str):
        #print('Cname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return None

    def get_annotation(self, i: int, name: str):
        return self.header

    def get_legend_title(self, i: int, name: str):
        return self.title

    def get_phase(self, i: int, name: str):
        return None

    def get_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):
        #self.ncolors = 1000
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def get_min_max(self, i: int, name: str):
        return self.min_value, self.max_value

    #def get_scalar(self, i: int, name: str) -> np.ndarray:
        #return self.scalar

    #------------
    # setters

    def set_data_format(self, i: int, name: str, data_format):
        self.data_format = data_format

    def set_min_max(self, i: int, name: str, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def set_scale(self, i: int, name: str, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')

    def set_legend_title(self, i: int, res_name: str, title: str) -> None:
        self.title = title

    #def set_phase(self, i: int, name: str, phase):
        #pass

    def set_nlabels_labelsize_ncolors_colormap(self, i: int, name: str, nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #------------
    # default getters
    def get_default_data_format(self, i: int, name: str):
        return self.data_format_default

    def get_default_min_max(self, i: int, name: str) -> tuple[float, float]:
        return self.min_default, self.max_default

    def get_default_scale(self, i: int, name: str):
        return 0.

    def get_default_phase(self, i: int, name: str):
        return None

    def get_default_legend_title(self, i: int, name: str) -> str:
        return self.title_default

    def get_default_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    #------------
    # unmodifyable getters
    def get_scale(self, i: int, name: str) -> float:
        return 0.

    def get_methods(self, i: int, name: str) -> list[str]:
        return ['centroid']

    def get_fringe_vector_result(self, i: int, name: str) -> tuple[None, None]:
        return None, None

    def is_normal_result(self, i: int, name: str) -> bool:
        """normal result is opt-in"""
        return True

    def __repr__(self):
        msg = 'NormalResult\n'
        msg += '    uname=%r\n' % self.uname
        return msg


class GuiResult(GuiResultCommon):
    deflects = False
    def __init__(self, subcase_id: int, header: str, title: str, location: str,
                 scalar: np.ndarray,
                 mask_value: Optional[int]=None, nlabels: Optional[int]=None,
                 labelsize: Optional[int]=None, ncolors: Optional[int]=None,
                 colormap: str='jet',
                 data_map: Any=None,
                 data_format: Optional[str]=None,
                 uname: str='GuiResult'):
        """
        Parameters
        ----------
        subcase_id : int
            the flag that points to self.subcases for a message
        header : str
            the sidebar word that goes in the lower left
        title : str
            the legend title
        location : str
            node, centroid
        scalar : (n,) int/float ndarray
            the data to make a contour plot with
        mask_value : int; default=None
            the NaN marker when scalars are ints
        data_format : str
            the type of data result (e.g. '%i', '%.2f', '%.3f')
        data_map : dict[???, ???]
            ???
        uname : str
            some unique name for ...

        """
        GuiResultCommon.__init__(self)

        self.data_map = data_map
        self.subcase_id = subcase_id
        #assert self.subcase_id > 0, self.subcase_id

        self.title = title
        self.header = header
        #self.scale = scale
        self.location = location
        assert location in ('node', 'centroid'), location
        self.subcase_id = subcase_id
        self.uname = uname

        if scalar is None:
            raise RuntimeError('title=%r scalar is None...' % self.title)
        assert scalar.shape[0] == scalar.size, 'shape=%s size=%s' % (str(scalar.shape), scalar.size)
        self.scalar = scalar
        #self.data_type = self.dxyz.dtype.str # '<c8', '<f4'
        self.data_type = self.scalar.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in REAL_TYPES else False
        self.is_complex = not self.is_real
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

        #print('title=%r data_type=%r' % (self.title, self.data_type))
        if self.data_type in INT_TYPES:
            self.data_format = '%i'
        elif data_format is None:
            self.data_format = '%.2f'
        else:
            self.data_format = data_format

        self.title_default = self.title
        self.header_default = self.header
        self.data_format_default = self.data_format

        self.min_default = np.nanmin(self.scalar)
        self.max_default = np.nanmax(self.scalar)
        if self.data_type in INT_TYPES:
            # turns out you can't have a NaN/inf with an integer array
            # we need to recast it
            if mask_value is not None:
                inan_short = np.where(self.scalar == mask_value)[0]
                if len(inan_short):
                    # overly complicated way to allow us to use ~inan to invert the array
                    inan = in1d(np.arange(len(self.scalar)), inan_short)
                    inan_remaining = self.scalar[~inan]

                    self.scalar = np.asarray(self.scalar, 'f')
                    self.data_type = self.scalar.dtype.str
                    self.data_format = '%.0f'
                    self.scalar[inan] = np.nan
                    try:
                        self.min_default = inan_remaining.min()
                    except ValueError:  # pragma: no cover
                        print('inan_remaining =', inan_remaining)
                        raise
                    self.max_default = inan_remaining.max()
        else:
            # handling VTK NaN oddinty
            # filtering the inf values and replacing them with NaN
            # 1.#R = inf
            # 1.#J = nan
            ifinite = np.isfinite(self.scalar)
            if not np.all(ifinite):
                self.scalar[~ifinite] = np.nan
                try:
                    self.min_default = self.scalar[ifinite].min()
                except ValueError:
                    print(self.title)
                    print(self.scalar)
                    raise
                self.max_default = self.scalar[ifinite].max()
        self.min_value = self.min_default
        self.max_value = self.max_default

    #------------
    def _validate(self, new):
        subcase_id = self.subcase_id
        header = self.header
        title = self.title
        location = self.location

        if isinstance(new, integer_float_types):
            return subcase_id, header, title, location

        assert self.location == new.location, f'location={self.location} new.location={new.location}'
        assert self.scalar.shape == new.scalar.shape, f'scalar.shape={self.scalar.shape} new.scalar.shape={new.scalar.shape}'
        return subcase_id, header, title, location

    def has_coord_transform(self, i: int, name: str) -> tuple[bool, list[str]]:
        return False, []
    def has_derivation_transform(self, i: int, res_name: str,
                                 ) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        return False, {}
    def has_nodal_combine_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []

    def __neg__(self):
        return self.__mul__(-1)

    def __pos__(self):
        return self.__mul__(1)

    def __radd__(self, new):
        self.__add__(new)

    def __add__(self, new):
        subcase_id, header, title, location = self._validate(new)
        if isinstance(new, integer_float_types):
            scalar = self.scalar + new
        else:
            scalar = self.scalar + new.scalar
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    def __rsub__(self, new):
        subcase_id, header, title, location = self._validate(new)
        if isinstance(new, integer_float_types):
            scalar = new - self.scalar
        else:
            scalar = new.scalar - self.scalar
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    def __sub__(self, new):
        return self.__add__(-new)

    def __rmul__(self, new):
        return self.__mul__(new)

    def __mul__(self, new):
        subcase_id, header, title, location = self._validate(new)
        if isinstance(new, integer_float_types):
            scalar = self.scalar * new
        elif isinstance(new, GuiResult):
            scalar = self.scalar * new.scalar
        else:
            raise NotImplementedError(new)
        #subcase_id, header, title, location = self._validate(new)
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    def __mod__(self, new):
        subcase_id, header, title, location = self._validate(new)
        if isinstance(new, integer_float_types):
            scalar = self.scalar % new
        elif isinstance(new, GuiResult):
            scalar = self.scalar % new.scalar
        else:
            raise NotImplementedError(new)
        #subcase_id, header, title, location = self._validate(new)
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    def __abs__(self):
        subcase_id, header, title, location = self._validate(self)
        scalar = np.abs(self.scalar)
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    def __rtruediv__(self, new):
        subcase_id, header, title, location = self._validate(new)
        if isinstance(new, integer_float_types):
            scalar = new / self.scalar
        else:
            raise NotImplementedError(new)
        #subcase_id, header, title, location = self._validate(new)
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    def __truediv__(self, new):
        return self.__mul__(1. / new)

    def __pow__(self, new):
        subcase_id, header, title, location = self._validate(self)
        if isinstance(new, integer_float_types):
            scalar = self.scalar ** new
        elif isinstance(new, GuiResult):
            scalar = self.scalar ** new.scalar
        else:
            raise NotImplementedError(new)
        return GuiResult(subcase_id, header, title, location, scalar,
                         mask_value=self.max_value,
                         nlabels=None, labelsize=None, ncolors=None,
                         colormap='jet', data_map=None, data_format=None, uname='GuiResult')

    #def __div__(self, new):
        #return self.__mul__(1. / new)

    #------------
    # getters
    #def export_hdf5_file(self, hdf5_file, exporter):
        #asd
    def get_data_type(self, i: int, name: str) -> str:
        #print('Aname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_type

    def get_data_format(self, i: int, name: str) -> str:
        #print('Bname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_format

    def get_location(self, i: int, name: str):
        #print('Cname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.location

    def get_annotation(self, i: int, name: str) -> str:
        return self.header

    def get_legend_title(self, i: int, name: str) -> str:
        """a title is the lagend label"""
        return self.title

    def get_phase(self, i: int, name: str) -> None:
        return None

    def get_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):
        #self.ncolors = 1000
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def get_min_max(self, i: int, name: str):
        return self.min_value, self.max_value

    #def get_scalar(self, i: int, name: str):
        #return self.scalar

    #------------
    # setters

    def set_data_format(self, i: int, name: str, data_format):
        self.data_format = data_format

    def set_min_max(self, i: int, name: str, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def set_scale(self, i: int, name: str, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')

    def set_title(self, i: int, name: str, title):
        self.title = title

    #def set_phase(self, i: int, name: str, phase):
        #pass

    def set_nlabels_labelsize_ncolors_colormap(self, i: int, name: str, nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #------------
    # default getters
    def get_default_data_format(self, i: int, name: str):
        return self.data_format_default

    def get_default_min_max(self, i: int, name: str) -> tuple[float, float]:
        return self.min_default, self.max_default

    def get_default_scale(self, i: int, name: str):
        return 0.

    def get_default_phase(self, i: int, name: str):
        return None

    def get_default_legend_title(self, i: int, name: str):
        return self.title_default

    def get_default_nlabels_labelsize_ncolors_colormap(self, i: int, name: str):
        # TODO: do this right
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    #------------
    # unmodifyable getters
    def get_scale(self, i: int, name: str) -> int:
        return 0.

    def get_methods(self, i: int, name: str) -> list[str]:
        if self.is_real:
            return [self.location]
        else:
            raise NotImplementedError('title=%s is not real; fmt=%s' % (self.title, self.data_type))
            #return ['node real', 'node imag', 'node magnitude', 'node phase']

    #def get_plot_value(self, i: int, name: str):
        #if self.is_real:
            #return self.dxyz[i, :]

    def get_vector_array_by_phase(self, i: int, unused_name, phase=0.):
        raise RuntimeError()
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            # e(i*theta) = cos(theta) + i*sin(theta)
            if self.dim == 2:
                # single result
                dxyz = self.dxyz
            elif self.dim == 3:
                dxyz = self.dxyz[i, :]
            else:
                raise NotImplementedError('dim=%s' % self.dim)
        else:
            dxyz = self._get_complex_displacements_by_phase(i, phase)
        assert len(dxyz.shape) == 2, dxyz.shape
        return xyz, dxyz

    def get_fringe_vector_result(self, i: int, name: str) -> tuple[np.ndarray, None]:
        if self.is_real:
            #return self.dxyz[i, :]
            return self.scalar, None
        else:
            raise NotImplementedError(f'title={self.title!r} is not real; fmt={self.data_type}')

    #def get_vector_result(self, i: int, name: str):
        #if self.is_real:
            #xyz = self.xyz + self.scales[i] * self.dxyz[i, :]
        #else:
            ## z = x + i*y
            ##
            #phase = self.phase
        #return self.xyz, xyz

    def save_vtk_result(self, used_titles: set[str]):
        titlei = self.title
        if self.subcase_id > 0:
            titlei = f'{self.title}_subcase={self.subcase_id:d}'

        from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk
        vtk_array = numpy_to_vtk(self.scalar, deep=0, array_type=None)

        title_out = titlei
        i = 1
        while title_out in used_titles:
            title_out = f'{titlei}_{i}'
            i += 1

        #if i != 1:
            #log.warning(f'duplicate GuiResult {titlei} because it is already used -> {title_out}')
        check_title(title_out, used_titles)
        vtk_array.SetName(title_out)
        return vtk_array

    def __repr__(self) -> str:
        msg = 'GuiResult\n'
        msg += '    title=%r\n' % self.title
        msg += '    data_type=%r\n' % self.data_type
        msg += '    uname=%r\n' % self.uname
        return msg


class GuiResultIDs(GuiResult):
    def __init__(self, subcase_id: int, header: str, title: str, location: str,
                 ids: Any, scalar: Any,
                 mask_value: Optional[int]=None, nlabels: Optional[int]=None,
                 labelsize: Optional[int]=None, ncolors: Optional[int]=None,
                 colormap: str='jet', data_map: Any=None,
                 data_format: Optional[str]=None, uname: str='GuiResult'):
        """a GuiResult with node/element ids"""
        self.ids = ids
        super().__init__(
            subcase_id, header, title, location, scalar,
            mask_value, nlabels, labelsize, ncolors, colormap, data_map,
            data_format, uname)


def check_title(title: str, used_titles: set[str]) -> None:
    assert title not in used_titles, title
    used_titles.add(title)
