from __future__ import print_function

from copy import deepcopy
from numpy import zeros
import numpy as np
from numpy.linalg import norm


class NastranComplexDisplacementResults(object):
    def __init__(self, subcase_id, titles, xyz, dxyz, scalar,
                 default_scale=40., uname='NastranGeometry'):
        self.subcase_id = subcase_id
        self.data_formats = ['%g'] * len(titles)
        self.xyz = xyz
        self.dxyz = dxyz
        self.dxyz_norm = norm(dxyz, axis=1)
        self.titles = titles
        self.scale = default_scale

        # displacement results can change the scale, so we need this
        # for defaulting the result and to not locking out the
        # displacement
        self.default_scale = default_scale

        # titles point to the scalar bar and thus can change
        self.titles_default = deepcopy(titles)

        # data formats are modified on the legend
        self.data_formats_default = deepcopy(self.data_formats)
        #self.default_scale = default_scales

        #theta = (2*np.pi * i/frame) % (2 * pi)
        theta = 0.0

        # calculate deflections
        eigvs = model.eigenvectors[1000].data[6, :, :]
        scale = 1.0
        defl = scale * (np.real(eigvs[:, :3]) * np.cos(theta) +
                        np.imag(eigvs[:, :3]) * np.sin(theta))


class NastranDisplacementResults(object):
    def __init__(self, subcase_id, titles, headers, xyz, dxyz, scalar,
                 scales, nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 deflects=True, uname='NastranGeometry'):
        self.subcase_id = subcase_id
        assert self.subcase_id > 0, self.subcase_id

        self.xyz = xyz
        self.dxyz = dxyz
        self.uname = uname
        #self.dxyz_norm = norm(dxyz, axis=1)

        self.deflects = deflects
        self.titles = titles
        self.headers = headers
        self.scales = scales
        self.subcase_id = subcase_id
        self.data_type = self.dxyz.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type == '<f4' else False
        self.is_complex = not self.is_real
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

        self.data_formats = None
        self.titles_default = None
        self.headers_default = None
        self.scales_default = None
        self.data_formats_default = None

        ntimes = self.dxyz.shape[0]
        self.default_mins = zeros(ntimes, dtype=self.dxyz.dtype)
        self.default_maxs = zeros(ntimes, dtype=self.dxyz.dtype)
        for itime in range(ntimes):
            normi = norm(self.dxyz[itime, :, :], axis=1)
            self.default_mins[itime] = normi.min()
            self.default_maxs[itime] = normi.max()

        self.max_values = None
        self.min_values = None

    def save_defaults(self):
        self.data_formats = ['%g'] * len(self.titles)
        self.titles_default = deepcopy(self.titles)
        self.headers_default = deepcopy(self.headers)
        self.data_formats_default = deepcopy(self.data_formats)
        self.scales_default = deepcopy(self.scales)

        size = self.default_maxs.size
        assert size == len(self.titles), 'len(maxs)=%s len(titles)=%s' % (
            size, len(self.titles))
        self.default_mins = self.default_mins.reshape(size)
        self.default_maxs = self.default_maxs.reshape(size)

        self.min_values = deepcopy(self.default_mins)
        self.max_values = deepcopy(self.default_maxs)
    #-------------------------------------
    # getters

    def get_header(self, i, name):
        #j = self.titles_default.index(name)
        #return self.titles[j]
        return self.headers[i]

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
        j = self.titles_default.index(name)
        self.scales[i] = scale

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

    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize, ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #def get_default_min_max(self, i, name):
        #return self.min_default[i], self.max_default[i]

    def get_default_scale(self, i, name):
        return self.scales_default[i]

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        # TODO: do this right
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    def get_default_title(self, i, name):
        return self.titles_default[i]

    #-------------------------------------
    # unmodifyable getters

    def get_data_type(self, i, name):
        return self.data_type

    def get_location(self, i, name):
        return 'node'

    def get_vector_size(self, i, name):
        #print(i)
        #j = self.titles_default.index(name)
        return 3

    #-------------------------------------

    def get_methods(self, i):
        if self.is_real:
            return ['node']
        else:
            return ['node real', 'node imag', 'node magnitude', 'node phase']

    def get_plot_value(self, i, name):
        if self.is_real:
            return self.dxyz[i, :]
        #if method == 'real':
            #return self.dxyz[i, :].real
        #elif method == 'imag':
            #return self.dxyz[i, :].imag
        #elif method == 'magnitude':
            #return abs(self.dxyz[i, :])
        #elif method == 'phase':
            #return angle(self.dxyz[i, :], deg=True)
        #else:
            #raise RuntimeError(method)

    def get_result(self, i, name):
        if self.is_real:
            # .0006 -> 0.0
            # .057 -> 0.0123
            # min
            return self.dxyz[i, :]
        else:
            raise NotImplementedError(self.is_real)
            #if method == 'real':
                #return self.dxyz[i, :].real
            #elif method == 'imag':
                #return self.dxyz[i, :].imag
            #elif method == 'magnitude':
                #return abs(self.dxyz[i, :])
            #elif method == 'phase':
                #return angle(self.dxyz[i, :], deg=True)
            #else:
                #raise RuntimeError(method)

    #@property
    #def scalar(self):
        #return self.dxyz_norm

    #def get_scalar(self, i, name):
        ##print(self.dxyz_norm)
        #return self.dxyz_norm

    def get_vector_result(self, i, name):
        if self.is_real:
            xyz = self.xyz + self.scales[i] * self.dxyz[i, :]
        else:
            # z = x + i*y
            #
            phase = self.phase
        return self.xyz, xyz

    def set_phase(self, i, name, phase):
        j = self.titles_default.index(name)
        self.phase[i] = phase

    def __repr__(self):
        msg = 'NastranDisplacementResults\n'
        msg += '    uname=%r\n' % self.uname
        return msg
