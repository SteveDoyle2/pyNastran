from copy import deepcopy

class GuiResult(object):
    def __init__(self, subcase_id, header, title, location, scalar, data_format=None, uname='GuiResult'):
        self.subcase_id = subcase_id
        assert self.subcase_id > 0, self.subcase_id

        self.deflects = False
        self.title = title
        self.header = header
        #self.scale = scale
        self.location = location
        assert location in ['node', 'centroid'], location
        self.subcase_id = subcase_id
        self.uname = uname

        self.scalar = scalar
        #self.data_type = self.dxyz.dtype.str # '<c8', '<f4'
        self.data_type = self.scalar.dtype.str # '<c8', '<f4'
        self.is_real = True if self.data_type in ['<i4', '<f4', '<f8'] else False
        self.is_complex = not self.is_real


        if self.data_type in ['<i4', '<f4']:
            self.data_format = '%i'
        else:
            self.data_format = '%g'

        self.title_default = deepcopy(self.title)
        self.header_default = deepcopy(self.header)
        self.data_format_default = deepcopy(self.data_format)
        #self.scale_default = deepcopy(self.scale)

    def get_data_type(self, i, name):
        return self.data_type

    def get_location(self, i, name):
        return self.location

    def get_title(self, i, name):
        return self.title

    def get_header(self, i, name):
        return self.header

    def get_default_title(self, i, name):
        return self.title_default

    def get_data_format(self, i, name):
        return self.data_format

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

    #def get_scalar(self, i, name):
        #return self.dxyz_norm

    #def get_vector_result(self, i, name):
        #if self.is_real:
            #xyz = self.xyz + self.scales[i] * self.dxyz[i, :]
        #else:
            ## z = x + i*y
            ##
            #phase = self.phase
        #return self.xyz, xyz

    def get_scale(self, i, name):
        return 0.
        #return self.scale

    def get_default_scale(self, i, name):
        return 0.
        #return self.scale_default

    def set_scale(self, i, name, scale):
        raise RuntimeError('This object cannot set a displacement scale factor.')
        #self.scale = scale

    #def set_phase(self, i, name, phase):
        #self.phase = phase

    def __repr__(self):
        msg = 'GuiResult\n'
        msg += '    uname=%r\n' % self.uname
        return msg
