from copy import deepcopy
from numpy import zeros
from numpy.linalg import norm  # type: ignore
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


class TransientElementResults:
    deflects = False
    def __init__(self, subcase_id, titles, headers, unused_scalars,
                 scales, data_formats=None,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                 set_max_min=False, uname='NastranGeometry'):
        """
        subcase_id : int
            the flag that points to self.subcases for a message
        headers : List[str]
            the sidebar word
        titles : List[str]
            the legend title
        scalars : (nnodes,n) float ndarray
            #the data to make a contour plot with
            does nothing
        scales : ???
            the deflection scale factors
        data_formats : List[str]
            the type of data result (e.g. '%i', '%.2f', '%.3f')
        uname : str
            some unique name for ...

        """
        self.subcase_id = subcase_id
        #assert self.subcase_id > 0, self.subcase_id

        self.xyz = xyz
        self.dxyz = dxyz
        self.dim = len(self.dxyz.shape)

        self.uname = uname
        #self.dxyz_norm = norm(dxyz, axis=1)

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
                raise NotImplementedError('result must be real')
                #: stored in degrees
                #self.phases = np.zeros(ntimes)
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
        self.data_formats_default = deepcopy(self.data_formats)

        size = self.default_maxs.size
        assert size == len(self.titles), 'len(maxs)=%s len(titles)=%s' % (
            size, len(self.titles))
        self.default_mins = self.default_mins.reshape(size)
        self.default_maxs = self.default_maxs.reshape(size)

        self.min_values = deepcopy(self.default_mins)
        self.max_values = deepcopy(self.default_maxs)
    #-------------------------------------
    # getters

    def get_header(self, i, unused_name):
        return self.headers[i]

    #def get_phase(self, i, name):
        #if self.is_real:
            #return None
        #return self.phases[i]

    def get_data_format(self, i, unused_name):
        return self.data_formats[i]

    def get_title(self, i, unused_name):
        return self.titles[i]

    def get_min_max(self, i, unused_name):
        return self.min_values[i], self.max_values[i]

    #-------------------------------------
    # setters

    def set_data_format(self, i, unused_name, data_format):
        self.data_formats[i] = data_format

    #def set_scale(self, i, name, scale):
        #self.scales[i] = scale

    #def set_phase(self, i, name, phase):
        #if self.is_real:
            #return
        #self.phases[i] = phase

    def set_title(self, i, unused_name, title):
        self.titles[i] = title

    def set_min_max(self, i, unused_name, min_value, max_value):
        self.min_values[i] = min_value
        self.max_values[i] = max_value

    #-------------------------------------
    # default getters
    def get_default_data_format(self, i, unused_name):
        return self.data_formats_default[i]

    def get_default_min_max(self, i, unused_name):
        return self.default_mins[i], self.default_maxs[i]

    def get_nlabels_labelsize_ncolors_colormap(self, unused_i, unused_name):
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def set_nlabels_labelsize_ncolors_colormap(self, unused_i, unused_name, nlabels, labelsize,
                                               ncolors, colormap):
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

    #def get_default_min_max(self, i, name):
        #return self.min_default[i], self.max_default[i]

    #def get_default_scale(self, i, name):
        #return self.scales_default[i]

    #def get_default_phase(self, i, name):
        #if self.is_real:
            #return None
        #return 0.0

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    def get_default_title(self, i, unused_sname):
        return self.titles_default[i]

    #-------------------------------------
    # unmodifyable getters

    def get_data_type(self, unused_i, unused_name):
        """the precision of the data"""
        return self.data_type

    def get_location(self, unused_i, unused_name):
        """the result type"""
        return 'centroid'

    def get_vector_size(self, unused_i, unused_name):
        """the result size"""
        return 3

    #-------------------------------------

    def get_methods(self, unused_i):
        if self.is_real:
            return ['magnitude', 'tx', 'ty', 'tz', 'rx', 'ry', 'rz']
        else:
            return ['node real', 'node imag', 'node magnitude', 'node phase']

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

    def get_result(self, i, unused_name):
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

    #@property
    #def scalar(self):
        #return self.dxyz_norm

    #def get_scalar(self, i, name):
        ##print(self.dxyz_norm)
        #return self.dxyz_norm

    def get_vector_result(self, i, name):
        assert len(self.xyz.shape) == 2, self.xyz.shape
        if self.is_real:
            xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                i, name, self.scales[i])
        else:
            xyz, deflected_xyz = self.get_vector_result_by_scale_phase(
                i, name, self.scales[i], self.phases[i])
        return xyz, deflected_xyz

    #def get_vector_result_by_scale_phase(self, i, name, scale, phase=0.):
        #"""
        #Gets the real/complex deflection result

        #Parameters
        #----------
        #i : int
            #mode/time/loadstep number
        #name : str
            #unused; useful for debugging
        #scale : float
            #deflection scale factor; true scale
        #phase : float; default=0.0
            #phase angle (degrees); unused for real results

        #Returns
        #-------
        #xyz : (nnodes, 3) float ndarray
            #the nominal state
        #deflected_xyz : (nnodes, 3) float ndarray
            #the deflected state
        #"""
        #assert len(self.xyz.shape) == 2, self.xyz.shape
        #if self.is_real:
            #if self.dim == 2:
                ## single result
                #deflected_xyz = self.xyz + scale * self.dxyz
            #elif self.dim == 3:
                #deflected_xyz = self.xyz + scale * self.dxyz[i, :]
            #else:
                #raise NotImplementedError('dim=%s' % self.dim)
        #else:
            #dxyz = self._get_complex_displacements_by_phase(i, phase)
            #deflected_xyz = self.xyz + scale * dxyz
        #assert len(deflected_xyz.shape) == 2, deflected_xyz.shape
        #return self.xyz, deflected_xyz

    def __repr__(self):
        """defines str(self)"""
        msg = 'TransientElementResults\n'
        msg += '    uname=%r\n' % self.uname
        return msg
