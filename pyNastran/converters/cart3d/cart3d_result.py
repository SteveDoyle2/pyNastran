"""
defines:
 - Cart3dGeometry
"""
from __future__ import print_function
from copy import deepcopy
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


class Cart3dGeometry(GuiResultCommon):
    """
    Stores the cart3d results.
    """
    def __init__(self, subcase_id, labels,
                 nodes, elements, regions, area, cnormals, colormap='jet',
                 uname='Cart3dGeometry'):
        GuiResultCommon.__init__(self)
        self.colormap_default = colormap
        self.uname = uname
        self.n = 2
        self.nodes = nodes # 0
        self.elements = elements # 1
        self.regions = regions # 2
        self.area = area
        self.centroid_normals = cnormals
        self.labels = labels
        self.data_formats = ['%i', '%i', '%i', '%.3f',
                             '%.3f', '%.3f', '%.3f']
        self.titles = ['NodeID', 'ElementID', 'Region', 'Area',
                       'NormalX', 'NormalY', 'NormalZ', ]
        self.result_types = ['NodeID', 'ElementID', 'Region', 'Area',
                             'NormalX', 'NormalY', 'NormalZ', ]
        self.subcase_id = subcase_id

        self.min_default = [
            nodes.min(), elements.min(), regions.min(), area.min(),
            cnormals[:, 0].min(), cnormals[:, 1].min(), cnormals[:, 2].min()]
        self.max_default = [
            nodes.max(), elements.max(), regions.max(), area.max(),
            cnormals[:, 0].max(), cnormals[:, 1].max(), cnormals[:, 2].max()]

        self.min_value = deepcopy(self.min_default)
        self.max_value = deepcopy(self.max_default)

        self.title_default = deepcopy(self.titles)
        self.data_format_default = deepcopy(self.data_formats)

        ntitles = len(self.titles)
        self.nlabels = [None] * ntitles
        self.labelsize = [None] * ntitles
        self.ncolors = [None] * ntitles
        self.colormap = [self.colormap_default] * ntitles

    def get_header(self, i, name):
        """
        header : str
            the sidebar word
        """
        return name

    def set_nlabels_labelsize_ncolors_colormap(self, i, name, nlabels, labelsize,
                                               ncolors, colormap):
        j = self.titles.index(name)
        self.nlabels[j] = nlabels
        self.labelsize[j] = labelsize
        self.ncolors[j] = ncolors
        self.colormap[j] = colormap

    def get_location(self, i, name):
        j = self.titles.index(name)
        if name == 'NodeID':
            return 'node'
        elif name == 'ElementID':
            return 'centroid'
        elif name == 'Region':
            return 'centroid'
        elif name == 'Area':
            return 'centroid'
        elif name in ['NormalX', 'NormalY', 'NormalZ',]:
            return 'centroid'
        raise NotImplementedError('i=%s' % str(i))

    #def get_scale(self, i, name):
        #j = self.titles.index(name)
        #return 0.0

    def get_vector_size(self, i, name):
        #j = self.titles.index(name)
        return 1

    def get_methods(self, i):
        if i == 1:
            return ['centroid']
        return ['node']

    def get_scalar(self, i, name):
        return self.get_result(i, name)

    def get_result(self, i, name):
        if name == 'NodeID':
            res = self.nodes
        elif name == 'ElementID':
            res = self.elements
        elif name == 'Region':
            res = self.regions
        elif name == 'Area':
            res = self.area
        elif name == 'NormalX':
            res = self.centroid_normals[:, 0]
        elif name == 'NormalY':
            res = self.centroid_normals[:, 1]
        elif name == 'NormalZ':
            res = self.centroid_normals[:, 2]
        else:
            raise NotImplementedError('i=%s' % str(i))
        return res

    #----------------------------------------------------
    # colormap
    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        try:
            j = self.titles.index(name)
        except ValueError:
            print('name=%s' % str(name))
            print('self.titles=%r' % self.titles)
            raise
        return self.nlabels[j], self.labelsize[j], self.ncolors[j], self.colormap[j]

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    #----------------------------------------------------
    # title
    def get_title(self, i, name):
        j = self.titles.index(name)
        return self.result_types[j]

    def set_title(self, i, name, title):
        j = self.titles.index(name)
        self.result_types[j] = title

    def get_default_title(self, i, name):
        j = self.titles.index(name)
        return self.title_default[j]

    #----------------------------------------------------
    # data format
    def get_data_format(self, i, name):
        j = self.titles.index(name)
        return self.data_formats[j]

    def set_data_format(self, i, name, data_format):
        j = self.titles.index(name)
        self.data_formats[j] = data_format

    def get_default_data_format(self, i, name):
        j = self.titles.index(name)
        return self.data_format_default[j]

    #----------------------------------------------------
    # min/max
    def get_min_max(self, i, name):
        j = self.titles.index(name)
        return self.min_value[j], self.max_value[j]

    def set_min_max(self, i, name, min_value, max_value):
        j = self.titles.index(name)
        self.min_value[j] = min_value
        self.max_value[j] = max_value

    def get_default_min_max(self, i, name):
        j = self.titles.index(name)
        return self.min_default[j], self.max_default[j]

    #----------------------------------------------------
    #def get_default_scale(self, i, name):
        #return 0.

    #def get_default_phase(self, i, name):
        #return None



    def __repr__(self):
        msg = 'Cart3dGeometry\n'
        msg += '    uname=%r\n' % self.uname
        msg += '    n=%r' % self.n
        return msg


class Cart3dResult(GuiResultCommon):
    """is this used???"""
    def __init__(self, subcase_id, labels, loads, colormap='jet', uname='Cart3dResult'):
        GuiResultCommon.__init__(self)
        self.uname = uname
        self.labels = labels
        self.n = 5
        self.loads = loads
        self.data_formats = ['%.3g', '%.3g', '%.3g', '%.3g', '%.3g']
        self.titles = ['rho', 'rhoU' 'rhoV', 'rhoW', 'rhoE']
        self.labels = labels

        self.min_default = loads.min(axis=1).tolist()
        assert len(self.min_default) == 5, len(self.min_default)
        self.max_default = loads.max(axis=1).tolist()

        self.min_value = deepcopy(self.min_default)
        self.max_value = deepcopy(self.max_default)
        self.data_format_default = deepcopy(self.data_formats)

        ntitles = len(self.titles)
        self.nlabels = [None] * ntitles
        self.labelsize = [None] * ntitles
        self.ncolors = [None] * ntitles
        self.colormap = [colormap] * ntitles

    def get_location(self, i, name):
        return 'centroid'

    def get_scalar(self, i, method):
        return self.get_result(i, method)

    def get_result(self, i, method):
        #print('method = %r' % method)
        ii, name = i
        return self.loads[name]
        #if name == 'rho':
            #return self.rho
        #elif name == 'rhoU':
            #return self.rhoU
        #elif name == 'rhoV':
            #return self.rhoV
        #elif name == 'rhoW':
            #return self.rhoW
        #elif name == 'rhoE':
            #return self.rhoE
        #raise NotImplementedError('i=%s' % i)

    #def get_data_fmt(self, i, name):
        #asdf

    def get_vector_size(self, i, name):
        j = self.titles.index(name)
        return 1

    def get_methods(self, i):
        if i == 1:
            return ['centroid']
        return ['node']

    #----------------------------------------------------
    # colormap
    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        j = self.titles.index(name)
        return self.nlabels[j], self.labelsize[j], self.ncolors[j], self.colormap[j]

    def get_default_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.get_nlabels_labelsize_ncolors_colormap(i, name)

    #----------------------------------------------------
    # title
    def get_title(self, i, name):
        #j = self.titles.index(name)
        return self.titles[i]

    def set_title(self, i, name, title):
        #j = self.titles.index(name)
        self.titles[i] = title

    def get_default_title(self, i, name):
        #j = self.titles.index(name)
        return self.titles[i]

    #----------------------------------------------------
    # data format
    def get_data_format(self, i, name):
        j = self.titles.index(name)
        return self.data_formats[j]

    def set_data_format(self, i, name, data_format):
        j = self.titles.index(name)
        self.data_formats[j] = data_format

    def get_default_data_format(self, i, name):
        j = self.titles.index(name)
        return self.data_format_default[j]

    #----------------------------------------------------
    # min/max
    def get_min_max(self, i, name):
        j = self.titles.index(name)
        return self.min_value[j], self.max_value[j]

    def set_min_max(self, i, name, min_value, max_value):
        j = self.titles.index(name)
        self.min_value[j] = min_value
        self.max_value[j] = max_value

    def get_default_min_max(self, i, name):
        j = self.titles.index(name)
        return self.min_default[j], self.max_default[j]

    #----------------------------------------------------

    def __repr__(self):
        msg = 'Cart3dResult\n'
        msg += '  i | Title\n'
        msg += '----+----------\n'
        #msg += '  0 | NodeID\n'
        #msg += '  1 | ElementID\n'
        msg += '  0 | Rho\n'
        msg += '  1 | RhoU\n'
        msg += '  2 | RhoV\n'
        msg += '  3 | RhoW\n'
        msg += '  4 | RhoE\n'
        msg += '----+----------\n'
        return msg
