"""
defines:
 - GuiResultCommon
 - GuiResult
"""
from __future__ import print_function
import numpy as np

REAL_TYPES = ['<i4', '<i8', '<f4', '<f8',
              '|i1', # this is a boolean
              '>i4', '>i8', '>f4', '>f8']
INT_TYPES = ['<i4', '<i8', '|i1',
             '>i4', '>i8']

class GuiResultCommon(object):
    def __init__(self):
        self.class_name = self.__class__.__name__

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

class NormalResult(GuiResultCommon):
    def __init__(self, subcase_id, header, title,
                 nlabels=2, labelsize=5, ncolors=2, colormap='jet',
                 data_format='%.1f', uname='GuiResult'):
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
    def __init__(self, subcase_id, header, title, location, scalar,
                 mask_value=None, nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 data_format=None, uname='GuiResult'):
        """
        subcase_id : int
            the flag that points to self.subcases for a message
        header : str
            the sidebar word
        title : str
            the legend title
        location : str
            node, centroid
        scalar : (n,) ndarray
            the data to make a contour plot with
        data_format : str
            the type of data result (e.g. '%i', '%.2f', '%.3f')
        uname : str
            some unique name for ...
        """
        GuiResultCommon.__init__(self)

        self.subcase_id = subcase_id
        #assert self.subcase_id > 0, self.subcase_id

        self.title = title
        self.header = header
        #self.scale = scale
        self.location = location
        assert location in ['node', 'centroid'], location
        self.subcase_id = subcase_id
        self.uname = uname

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

        self.min_default = self.scalar.min()
        self.max_default = self.scalar.max()
        if self.data_type in INT_TYPES:
            # turns out you can't have a NaN/inf with an integer array
            # we need to recast it
            if mask_value is not None:
                inan_short = np.where(self.scalar == mask_value)[0]
                if len(inan_short):
                    # overly complicated way to allow us to use ~inan to invert the array
                    inan = np.in1d(np.arange(len(self.scalar)), inan_short)

                    self.scalar = np.asarray(self.scalar, 'f')
                    self.data_type = self.scalar.dtype.str
                    self.data_format = '%.0f'
                    self.scalar[inan] = np.nan
                    self.min_default = self.scalar[~inan].min()
                    self.max_default = self.scalar[~inan].max()
        else:
            # handling VTK NaN oddinty
            # filtering the inf values and replacing them with NaN
            # 1.#R = inf
            # 1.#J = nan
            ifinite = np.isfinite(self.scalar)
            if not np.all(ifinite):
                self.scalar[~ifinite] = np.nan
                self.min_default = self.scalar[ifinite].min()
                self.max_default = self.scalar[ifinite].max()
        self.min_value = self.min_default
        self.max_value = self.max_default

    #------------
    # getters
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
        msg += '    uname=%r\n' % self.uname
        return msg
