"""
defines:
 - GuiResultCommon
 - GuiResult

"""
from typing import Any, Optional
import numpy as np
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

    #def get_data_type(self, i, name):
        #raise NotImplementedError(self.class_name)

    #def get_header(self, i, name):
        #raise NotImplementedError(self.class_name)

    def is_normal_result(self, i, name):
        return False

    def get_data_format(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_location(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_title(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_min_max(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_scalar(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_methods(self, i):
        raise NotImplementedError(self.class_name)

    def get_result(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_vector_size(self, i, name):
        return 1

    def get_scale(self, i, name):
        return 0.

    def get_phase(self, i, name):
        return None

    #------------
    # setters

    def set_data_format(self, i, name, data_format):
        raise NotImplementedError(self.class_name)

    def set_min_max(self, i, name, min_value, max_value):
        raise NotImplementedError(self.class_name)

    def set_title(self, i, name, title):
        raise NotImplementedError(self.class_name)

    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize,
                                               ncolors, colormap):
        raise NotImplementedError(self.class_name)

    def set_scale(self, i, name, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')

    def set_phase(self, i, name, phase):
        pass

    #------------
    # default getters
    def get_default_data_format(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_default_min_max(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_default_title(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        raise NotImplementedError(self.class_name)

    def get_default_scale(self, i, name):
        return 0.

    def get_default_phase(self, i, name):
        return None

    def _get_complex_displacements_by_phase(self, i, phase):
        raise NotImplementedError(self.class_name)


class NullResult(GuiResultCommon):
    def __init__(self):
        super(NullResult, self).__init__()

    def get_scalar(self, i, name):
        return None

    def __repr__(self):
        msg = '<NormalResult>'
        return msg

class GridPointForceResult(GuiResultCommon):
    def __init__(self, subcase_id, header, title, gpforce_array, uname='GridPointForceResult'):
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

    def get_scalar(self, i, name):
        return None

    def get_result(self, i, name):
        return None
    def get_title(self, i, name):
        return self.title
    def get_location(self, i, name):
        return self.location
    def get_header(self, i, name):
        return self.header
    def get_methods(self, i):
        return None
    def get_data_format(self, i, name):
        return None
    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        return None, None, None, None
    def get_default_data_format(self, i, name):
        return None
    def get_default_min_max(self, i, name):
        return None, None
    def get_default_title(self, i, name):
        return self.title
    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        return None, None, None, None
    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize, ncolors, colormap):
        return
    def get_min_max(self, i, name):
        return None, None

    #def get_default_min_max(self, i, name):
        #return None, None

    def __repr__(self):
        msg = '<GridPointForceResult>'
        return msg

class NormalResult(GuiResultCommon):
    def __init__(self, subcase_id, header, title,
                 nlabels=2, labelsize=5, ncolors=2, colormap='jet',
                 data_format='%.1f', uname='NormalResult'):
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

    def get_data_type(self, i, name):
        #print('Aname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_type

    def get_data_format(self, i, name):
        #print('Bname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_format

    def get_location(self, i, name):
        #print('Cname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return None

    def get_header(self, i, name):
        return self.header

    def get_title(self, i, name):
        return self.title

    def get_phase(self, i, name):
        return None

    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        #self.ncolors = 1000
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def get_min_max(self, i, name):
        return self.min_value, self.max_value

    def get_scalar(self, i, name):
        return self.scalar

    #------------
    # setters

    def set_data_format(self, i, name, data_format):
        self.data_format = data_format

    def set_min_max(self, i, name, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def set_scale(self, i, name, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')

    def set_title(self, i, name, title):
        self.title = title

    #def set_phase(self, i, name, phase):
        #pass

    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #------------
    # default getters
    def get_default_data_format(self, i, name):
        return self.data_format_default

    def get_default_min_max(self, i, name):
        return self.min_default, self.max_default

    def get_default_scale(self, i, name):
        return 0.

    def get_default_phase(self, i, name):
        return None

    def get_default_title(self, i, name):
        return self.title_default

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    #------------
    # unmodifyable getters
    def get_scale(self, i, name):
        return 0.

    def get_vector_size(self, i, name):
        return 1

    def get_methods(self, i):
        return None

    def get_result(self, i, name):
        return None

    def is_normal_result(self, i, name):
        return True

    def __repr__(self):
        msg = 'NormalResult\n'
        msg += '    uname=%r\n' % self.uname
        return msg


class GuiResult(GuiResultCommon):
    deflects = False
    def __init__(self, subcase_id: int, header: str, title: str, location: str, scalar: Any,
                 mask_value: Optional[int]=None, nlabels: Optional[int]=None,
                 labelsize: Optional[int]=None, ncolors: Optional[int]=None,
                 colormap: str='jet', data_map: Any=None,
                 data_format: Optional[str]=None, uname: str='GuiResult'):
        """
        Parameters
        ----------
        subcase_id : int
            the flag that points to self.subcases for a message
        header : str
            the sidebar word
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
        assert location in ['node', 'centroid'], location
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
                    inan = np.in1d(np.arange(len(self.scalar)), inan_short)
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
    def get_data_type(self, i, name):
        #print('Aname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_type

    def get_data_format(self, i, name):
        #print('Bname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.data_format

    def get_location(self, i, name):
        #print('Cname=%r data_type=%s fmt=%s' % (self.title, self.data_type, self.data_format))
        return self.location

    def get_header(self, i, name):
        return self.header

    def get_title(self, i, name):
        return self.title

    def get_phase(self, i, name):
        return None

    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        #self.ncolors = 1000
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def get_min_max(self, i, name):
        return self.min_value, self.max_value

    def get_scalar(self, i, name):
        return self.scalar

    #------------
    # setters

    def set_data_format(self, i, name, data_format):
        self.data_format = data_format

    def set_min_max(self, i, name, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def set_scale(self, i, name, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')

    def set_title(self, i, name, title):
        self.title = title

    #def set_phase(self, i, name, phase):
        #pass

    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #------------
    # default getters
    def get_default_data_format(self, i, name):
        return self.data_format_default

    def get_default_min_max(self, i, name):
        return self.min_default, self.max_default

    def get_default_scale(self, i, name):
        return 0.

    def get_default_phase(self, i, name):
        return None

    def get_default_title(self, i, name):
        return self.title_default

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        # TODO: do this right
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    #------------
    # unmodifyable getters
    def get_scale(self, i, name):
        return 0.

    def get_vector_size(self, i, name):
        return 1

    def get_methods(self, i):
        if self.is_real:
            return self.location
        else:
            raise NotImplementedError('title=%s is not real; fmt=%s' % (self.title, self.data_type))
            #return ['node real', 'node imag', 'node magnitude', 'node phase']

    #def get_plot_value(self, i, name):
        #if self.is_real:
            #return self.dxyz[i, :]

    def get_vector_array_by_phase(self, i, unused_name, phase=0.):
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

    def get_result(self, i, name):
        if self.is_real:
            #return self.dxyz[i, :]
            return self.scalar
        else:
            raise NotImplementedError('title=%s is not real; fmt=%s' % (self.title, self.data_type))

    #def get_vector_result(self, i, name):
        #if self.is_real:
            #xyz = self.xyz + self.scales[i] * self.dxyz[i, :]
        #else:
            ## z = x + i*y
            ##
            #phase = self.phase
        #return self.xyz, xyz

    def __repr__(self):
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
