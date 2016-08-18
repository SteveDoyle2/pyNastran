from copy import deepcopy


class GuiResult(object):
    def __init__(self, subcase_id, header, title, location, scalar,
                 nlabels=None, labelsize=None, ncolors=None, colormap='jet',
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
        self.subcase_id = subcase_id
        #assert self.subcase_id > 0, self.subcase_id

        self.deflects = False
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
        self.is_real = True if self.data_type in ['<i4', '<i8', '<f4', '<f8', '|i1'] else False
        self.is_complex = not self.is_real
        self.nlabels = nlabels
        self.labelsize = labelsize
        self.ncolors = ncolors
        self.colormap = colormap

        #print('title=%r data_type=%r' % (self.title, self.data_type))
        if self.data_type in ['<i4', '<i8', '|i1']:
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

    def set_title(self, i, name):
        return self.title

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
