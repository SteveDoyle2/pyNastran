class Cart3dGeometry(object):
    def __init__(self, subcase_id, labels,
                 nodes, elements, regions, area, cnormals,
                 uname='Cart3dGeometry'):
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

        self.min_value = self.min_default
        self.max_value = self.max_default


        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = 'jet'

    def get_header(self, i, name):
        """
        header : str
            the sidebar word
        """
        return name

    def get_min_max(self, i, name):
        j = self.titles.index(name)
        return self.min_value[j], self.max_value[j]

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

    #def get_data(self, i, name):
        #j = self.titles.index(name)
        #if name == 'NodeID':
            #return self.nodes
        #elif name == 'ElementID':
            #return self.elements
        #elif name == 'Region':
            #return self.regions
        #raise NotImplementedError('i=%s' % str(i))

    def get_scale(self, i, name):
        j = self.titles.index(name)
        return 0.0

    def get_title(self, i, name):
        j = self.titles.index(name)
        return self.result_types[j]

    #def get_data_fmt(self, i, name):
        #asdf

    def get_data_format(self, i, name):
        j = self.titles.index(name)
        return self.data_formats[j]

    def get_vector_size(self, i, name):
        j = self.titles.index(name)
        return 1

    def get_methods(self, i):
        if i == 1:
            return ['centroid']
        return ['node']

    def get_result(self, i, name):
        if name == 'NodeID':
            return self.nodes
        elif name == 'ElementID':
            return self.elements
        elif name == 'Region':
            return self.regions
        elif name == 'Area':
            return self.area
        elif name == 'NormalX':
            return self.centroid_normals[:, 0]
        elif name == 'NormalY':
            return self.centroid_normals[:, 1]
        elif name == 'NormalZ':
            return self.centroid_normals[:, 2]
        raise NotImplementedError('i=%s' % str(i))

    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

    def __repr__(self):
        msg = 'Cart3dGeometry\n'
        msg += '    uname=%r\n' % self.uname
        msg += '    n=%r' % self.n
        return msg

    def get_phase(self, i, name):
        return None

    def set_phase(self, i, name):
        pass


class Cart3dResult(object):
    def __init__(self, subcase_id, labels, loads, uname='Cart3dResult'):
        self.uname = uname
        self.labels = labels
        self.n = 5
        self.loads = loads
        self.data_formats = ['%.3g', '%.3g', '%.3g', '%.3g', '%.3g']
        self.titles = ['rho', 'rhoU' 'rhoV', 'rhoW', 'rhoE']
        self.labels = labels

        self.nlabels = None
        self.labelsize = None
        self.ncolors = None
        self.colormap = 'jet'

    def get_result(self, i, method):
        print('method = %r' % method)
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

    def get_title(self, i, name):
        j = self.titles.index(name)
        return self.result_types[j]

    def get_data_fmt(self, i, name):
        asdf

    def get_data_format(self, i, name):
        j = self.titles.index(name)
        return self.data_formats[j]

    def get_vector_size(self, i, name):
        j = self.titles.index(name)
        return 1

    def get_methods(self, i):
        if i == 1:
            return ['centroid']
        return ['node']

    def get_nlabels_labelsize_ncolors_colormap(self, i, name):
        return self.nlabels, self.labelsize, self.ncolors, self.colormap

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

    def get_phase(self, i, name):
        return None

    def set_phase(self, i, name):
        pass
