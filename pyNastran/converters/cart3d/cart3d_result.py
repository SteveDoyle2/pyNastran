"""
defines:
 - Cart3dGeometry

"""
from typing import Any
from copy import deepcopy
import numpy as np
from pyNastran.gui.gui_objects.gui_result import GuiResultCommon


class Cart3dGeometry(GuiResultCommon):
    """Stores the cart3d results."""
    def __init__(self, subcase_id, labels,
                 nodes, elements, regions, area, cnormals,
                 colormap: str='jet',
                 uname: str='Cart3dGeometry'):
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

    def get_annotation(self, i, name):
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
        #j = self.titles.index(name)
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

    def has_coord_transform(self, i: int, name: str) -> tuple[bool, list[str]]:
        return False, []
    def has_derivation_transform(self, i: int, resname: str) -> tuple[bool, dict[str, Any]]:
        """min/max/avg"""
        return False, {}
    def has_nodal_combine_transform(self, i: int, resname: str) -> tuple[bool, list[str]]:
        """elemental -> nodal"""
        return False, []
    #def has_output_checks(self, i: int, resname: str) -> tuple[bool, bool, bool]:
        #is_enabled_fringe = False
        #is_checked_fringe = True
        #is_enabled_disp = False
        #is_checked_disp = False
        #is_enabled_vector = False
        #is_checked_vector = False
        #out = (
            #is_enabled_fringe, is_checked_fringe,
            #is_enabled_disp, is_checked_disp,
            #is_enabled_vector, is_checked_vector)
        #return out

    def get_methods(self, i: str, name: str):
        if i == 1:
            return ['centroid']
        return ['node']

    def get_fringe_vector_result(self, i, name) -> tuple[np.ndarray, None]:
        fringe = self.get_fringe_result(i, name)
        return fringe, None

    def get_fringe_result(self, i, name) -> np.ndarray:
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
        else:  # pragma: no cover
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
    def get_legend_title(self, i, name):
        j = self.titles.index(name)
        return self.result_types[j]

    def set_legend_title(self, i, name, title):
        j = self.titles.index(name)
        self.result_types[j] = title

    def get_default_legend_title(self, i, name):
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
    def get_imin_imax(self, i, name) -> tuple[None, None]:
        return None, None
    def get_min_max(self, i, name) -> tuple[float, float]:
        j = self.titles.index(name)
        return self.min_value[j], self.max_value[j]

    def set_min_max(self, i, name, min_value, max_value):
        j = self.titles.index(name)
        self.min_value[j] = min_value
        self.max_value[j] = max_value

    def get_default_min_max(self, i, name) -> tuple[float, float]:
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
