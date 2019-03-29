# coding: utf-8
from __future__ import division, unicode_literals, print_function

from collections import OrderedDict

from six import string_types
import numpy as np

import vtk

from pyNastran.gui.qt_version import qt_version
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.qt_files.gui_qt_common import GuiQtCommon
from pyNastran.gui.qt_files.mark_actions import create_annotation
from pyNastran.gui.menus.menus import Group


#class Interactor(vtk.vtkGenericRenderWindowInteractor):
    #def __init__(self):
        ##vtk.vtkGenericRenderWindowInteractor()
        #pass

    #def HighlightProp(self):
        #print('highlight')


#class PyNastranRenderWindowInteractor(QVTKRenderWindowInteractor):
    #def __init__(self, parent=None):

        #render_window = vtk.vtkRenderWindow()
        #iren = Interactor()
        #iren.SetRenderWindow(render_window)
        #kwargs = {
            #'iren' : iren,
            #'rw' : render_window,
        #}
        #QVTKRenderWindowInteractor.__init__(self, parent=parent,
                                            #iren=iren, rw=render_window)
        #self.Highlight

class GuiVTKCommon(GuiQtCommon):
    """this class has VTK functionality, but no interactive/menu capability"""
    def __init__(self, **kwds):
        if qt_version == 'pyqt4':
            GuiQtCommon.__init__(self, **kwds)
        elif qt_version == 'pyqt5':
            super(GuiVTKCommon, self).__init__(**kwds)
        elif qt_version in ['pyside', 'pyside2']:
            GuiQtCommon.__init__(self, **kwds)
        else:  #: pragma: no cover
            raise NotImplementedError(qt_version)

    def fake_init(self):
        """vtk setup"""
        self.create_vtk_actors(create_rend=False)
        #self._create_vtk_objects()
        #self.build_vtk_frame()
    #---------------------------------------------------------------------------
    # basic init
    def create_vtk_actors(self, create_rend=True):
        """creates the vtk actors used by the GUI"""
        if create_rend:
            self.rend = vtk.vtkRenderer()

        xyz = [0., 0., 0.]
        self.min_max_actors.append(create_annotation(self, 'Min', *xyz))
        self.min_max_actors.append(create_annotation(self, 'Max', *xyz))
        for actor in self.min_max_actors:
            actor.SetVisibility(False)

        # vtk actors
        self.grid = vtk.vtkUnstructuredGrid()

        # edges
        self.edge_actor = vtk.vtkLODActor()
        self.edge_actor.DragableOff()
        self.edge_mapper = vtk.vtkPolyDataMapper()

        self.create_cell_picker()

    def create_cell_picker(self):
        """creates the vtk picker objects"""
        self.cell_picker = vtk.vtkCellPicker()
        self.node_picker = vtk.vtkPointPicker()

        self.area_picker = vtk.vtkAreaPicker()  # vtkRenderedAreaPicker?
        self.rubber_band_style = vtk.vtkInteractorStyleRubberBandPick()
        #vtk.vtkInteractorStyleRubberBand2D
        #vtk.vtkInteractorStyleRubberBand3D
        #vtk.vtkInteractorStyleRubberBandZoom
        #vtk.vtkInteractorStyleAreaSelectHover
        #vtk.vtkInteractorStyleDrawPolygon

        #vtk.vtkAngleWidget
        #vtk.vtkAngleRepresentation2D
        #vtk.vtkAngleRepresentation3D
        #vtk.vtkAnnotation
        #vtk.vtkArrowSource
        #vtk.vtkGlyph2D
        #vtk.vtkGlyph3D
        #vtk.vtkHedgeHog
        #vtk.vtkLegendBoxActor
        #vtk.vtkLegendScaleActor
        #vtk.vtkLabelPlacer

        self.cell_picker.SetTolerance(0.001)
        self.node_picker.SetTolerance(0.001)

    def _build_vtk_frame_post(self):
        self.build_lookup_table()

        text_size = self.settings.set_text_size # was 14
        dtext_size = text_size + 1
        self.create_text([5, 5 + 3 * dtext_size], 'Max  ', text_size)  # text actor 0
        self.create_text([5, 5 + 2 * dtext_size], 'Min  ', text_size)  # text actor 1
        self.create_text([5, 5 + 1 * dtext_size], 'Word1', text_size)  # text actor 2
        self.create_text([5, 5], 'Word2', text_size)  # text actor 3

        self.get_edges()
        if self.is_edges:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()

    #---------------------------------------------------------------------------
    # properties

    @property
    def logo(self):
        """Gets the pyNastran icon path, which can be overwritten"""
        return self._logo

    @logo.setter
    def logo(self, logo):
        """Sets the pyNastran icon path, which can be overwritten"""
        self._logo = logo

    @property
    def legend_shown(self):
        """determines if the legend is shown"""
        return self.scalar_bar.is_shown

    @property
    def scalar_bar_actor(self):
        """gets the scalar bar actor"""
        return self.scalar_bar.scalar_bar

    @property
    def color_function(self):
        """gets the scalar bar's color function"""
        return self.scalar_bar.color_function

    @property
    def result_name(self):
        """
        creates the self.result_name variable
        """
        # case_key = (1, 'ElementID', 1, 'centroid', '%.0f')
        case_key = self.case_keys[self.icase]
        assert isinstance(case_key, integer_types), case_key
        unused_obj, (unused_i, res_name) = self.result_cases[case_key]
        return res_name

    #---------------------------------------------------------------------------
    # basic interaction
    def on_show_debug(self):
        """sets a flag for showing/hiding DEBUG messages"""
        self.settings.show_debug = not self.settings.show_debug

    def on_show_info(self):
        """sets a flag for showing/hiding INFO messages"""
        self.settings.show_info = not self.settings.show_info

    def on_show_command(self):
        """sets a flag for showing/hiding COMMAND messages"""
        self.settings.show_command = not self.settings.show_command

    def on_show_warning(self):
        """sets a flag for showing/hiding WARNING messages"""
        self.settings.show_warning = not self.settings.show_warning

    def on_show_error(self):
        """sets a flag for showing/hiding ERROR messages"""
        self.settings.show_error = not self.settings.show_error

    def hide_axes(self, cids=None):
        """
        Show a set of coordinate systems

        ..todo :: fix the coords
        """
        if cids is None:
            cids = self.axes.keys()
        for cid in cids:
            axis = self.axes[cid]
            axis.VisibilityOff()
        self.corner_axis.EnabledOff()

    def show_axes(self, cids=None):
        """
        Show a set of coordinate systems

        ..todo :: fix the coords
        """
        if cids is None:
            cids = self.axes.keys()
        for cid in cids:
            axis = self.axes[cid]
            axis.VisibilityOn()
        self.corner_axis.EnabledOn()

    def delete_actor(self, name):
        """deletes an actor and associated properties"""
        if name != 'main':
            if name in self.geometry_actors:
                actor = self.geometry_actors[name]
                self.rend.RemoveActor(actor)
                del self.geometry_actors[name]

            if name in self.geometry_properties:
                unused_prop = self.geometry_properties[name]
                del self.geometry_properties[name]
            self.Render()

    def show_only(self, names):
        """
        Show these actors only

        names : str, List[str]
            names to show
            If they're hidden, show them.
            If they're shown and shouldn't be, hide them.

        ..todo :: update the GeomeryProperties
        """
        raise NotImplementedError('show_only')

    def hide_actors(self, except_names=None):
        """
        Hide all the actors

        except_names : str, List[str], None
            list of names to exclude
            None : hide all

        ..note :: If an actor is hidden and in the except_names, it will still be hidden.
        ..todo :: update the GeomeryProperties
        """
        if except_names is None:
            except_names = []
        elif isinstance(except_names, string_types):
            except_names = [except_names]

        # hide everything but the main grid
        for key, actor in self.geometry_actors.items():
            if key not in except_names:
                actor.VisibilityOff()

        self.hide_axes()
        self.hide_legend()
        #self.settings.set_background_color_to_white()

    #---------------------------------------------------------------------------
    def Render(self):
        """Renders the GUI"""
        #self.vtk_interactor.Render()
        self.vtk_interactor.GetRenderWindow().Render()

    def cell_centroid(self, cell_id):
        """gets the cell centroid"""
        cell = self.grid_selected.GetCell(cell_id)
        nnodes = cell.GetNumberOfPoints()
        points = cell.GetPoints()
        centroid = np.zeros(3, dtype='float32')
        for ipoint in range(nnodes):
            point = np.array(points.GetPoint(ipoint), dtype='float32')
            centroid += point
        centroid /= nnodes
        return centroid

    def build_lookup_table(self):
        """build the nominal lookup table for the scalar bar"""
        scalar_range = self.grid_selected.GetScalarRange()
        self.grid_mapper.SetScalarRange(scalar_range)
        self.grid_mapper.SetLookupTable(self.color_function)
        self.rend.AddActor(self.scalar_bar_actor)

    def _set_results(self, form, cases):
        assert len(cases) > 0, cases
        if isinstance(cases, OrderedDict):
            self.case_keys = list(cases.keys())
        else:
            self.case_keys = sorted(cases.keys())
            assert isinstance(cases, dict), type(cases)

        self.result_cases = cases

        if len(self.case_keys) > 1:
            self.icase = -1
            self.ncases = len(self.result_cases)  # number of keys in dictionary
        elif len(self.case_keys) == 1:
            self.icase = -1
            self.ncases = 1
        else:
            self.icase = -1
            self.ncases = 0
        self.icase_disp = None
        self.icase_vector = None
        self.icase_fringe = None
        self.set_form(form)

    #---------------------------------------------------------------------------
    # groups
    def create_group_with_name(self, name, eids, model_name='main'):
        """
        Creates a group from the root model (model_name)

        Parameters
        ----------
        name : str
            the name of the group
        eids : (neids, ) int ndarray
            the elements in the group
        model_name : str; default='main'
            the name of the parent model
        """
        elements_pound = self.groups[model_name].elements_pound
        element_str = ''
        group = Group(
            name, element_str, elements_pound,
            editable=True)

        # TODO: make sure all the eids exist
        group.element_ids = eids
        self.log_command('create_group_with_name(%r, %r)' % (name, eids))
        self.groups[name] = group
