# coding: utf-8
from typing import Optional
import numpy as np

from pyNastran.gui.vtk_rendering_core import (
    vtkRenderer, # vtkRenderWindow, vtkRenderWindowInteractor,
    #vtkActor, vtkCamera,
    #vtkDataSetMapper,
    vtkPolyDataMapper)
#from vtk import (vtkLODActor,
                 #vtkCellPicker, vtkPointPicker, vtkAreaPicker, vtkDataSetMapper,
                 #vtkInteractorStyleRubberBandPick,
                 #vtkArrowSource,
                 #vtkGlyph3D, vtkExtractEdges,
#)
from vtkmodules.vtkRenderingLOD import vtkLODActor
from vtkmodules.vtkRenderingCore import vtkCellPicker, vtkPointPicker, vtkAreaPicker, vtkDataSetMapper
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleRubberBandPick

from vtkmodules.vtkFiltersCore import vtkGlyph3D
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

from pyNastran.gui.qt_version import qt_version
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.qt_files.gui_qt_common import GuiQtCommon
from pyNastran.gui.qt_files.mark_actions import create_annotation
from pyNastran.gui.utils.vtk.vtk_vector import build_glyph
from pyNastran.gui.utils.vtk.vtk_edges import create_edges_from_grid
from pyNastran.gui.utils.vtk.vtk_utils import map_element_centroid_to_node_fringe_result
from pyNastran.gui.menus.menus import Group


#class Interactor(vtkGenericRenderWindowInteractor):
    #def __init__(self):
        ##vtkGenericRenderWindowInteractor()
        #pass

    #def HighlightProp(self):
        #print('highlight')


#class PyNastranRenderWindowInteractor(QVTKRenderWindowInteractor):
    #def __init__(self, parent=None):

        #render_window = vtkRenderWindow()
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
        if qt_version in {'pyqt5', 'pyqt6'}:
            super(GuiVTKCommon, self).__init__(**kwds)
        elif qt_version in {'pyside2', 'pyside6'}:
            GuiQtCommon.__init__(self, **kwds)
        else:  #: pragma: no cover
            raise NotImplementedError(qt_version)

    def fake_init(self) -> None:
        """vtk setup"""
        self.create_vtk_actors(create_rend=False)
        #self._create_vtk_objects()

        #self.build_vtk_frame()
        #self.add_geometry()
        self._build_vtk_frame_post(build_lookup_table=False)
    #---------------------------------------------------------------------------
    # basic init
    def create_vtk_actors(self, create_rend: bool=True) -> False:
        """creates the vtk actors used by the GUI"""
        if create_rend:
            self.rend = vtkRenderer()

        xyz = [0., 0., 0.]
        self.min_max_actors.append(create_annotation(self, 'Min', *xyz))
        self.min_max_actors.append(create_annotation(self, 'Max', *xyz))
        for actor in self.min_max_actors:
            actor.SetVisibility(False)

        # vtk actors
        self.grid = vtkUnstructuredGrid()

        # edges
        self.edge_actor = vtkLODActor()
        self.edge_actor.DragableOff()
        self.edge_mapper = vtkPolyDataMapper()

        self.create_cell_picker()

    def create_cell_picker(self) -> None:
        """creates the vtk picker objects"""
        self.cell_picker = vtkCellPicker()
        self.node_picker = vtkPointPicker()

        self.area_picker = vtkAreaPicker()  # vtkRenderedAreaPicker?
        self.rubber_band_style = vtkInteractorStyleRubberBandPick()
        #vtkInteractorStyleRubberBand2D
        #vtkInteractorStyleRubberBand3D
        #vtkInteractorStyleRubberBandZoom
        #vtkInteractorStyleAreaSelectHover
        #vtkInteractorStyleDrawPolygon

        #vtkAngleWidget
        #vtkAngleRepresentation2D
        #vtkAngleRepresentation3D
        #vtkAnnotation
        #vtkArrowSource
        #vtkGlyph2D
        #vtkGlyph3D
        #vtkHedgeHog
        #vtkLegendBoxActor
        #vtkLegendScaleActor
        #vtkLabelPlacer

        self.cell_picker.SetTolerance(0.001)
        self.node_picker.SetTolerance(0.001)

    def add_geometry(self) -> None:
        """
        #(N,)  for stress, x-disp
        #(N,3) for warp vectors/glyphs
        grid_result = vtkFloatArray()

        point_data: vtkPointData  = self.grid.GetPointData()
        cell_data: vtkCellData = self.grid.GetCellData()

        self.grid.GetCellData().SetScalars(grid_result)
        self.grid.GetPointData().SetScalars(grid_result)


        self.grid_mapper   <-input-> self.grid
        vtkDataSetMapper() <-input-> vtkUnstructuredGrid()

        self.grid_mapper   <--map--> self.geom_actor <-add-> self.rend
        vtkDataSetMapper() <--map--> vtkActor()      <-add-> vtkRenderer()
        """
        if self.is_groups:
            # solid_bending: eids 1-182
            self._setup_element_mask()
            #eids = np.arange(172)
            #eids = arange(172)
            #self.update_element_mask(eids)
        else:
            self.grid_selected = self.grid
        #print('grid_selected =', self.grid_selected)

        self.grid_mapper = vtkDataSetMapper()
        self.grid_mapper.SetInputData(self.grid_selected)

        self._make_contour_filter()

        #if 0:
            #self.warp_filter = vtkWarpVector()
            #self.warp_filter.SetScaleFactor(50.0)
            #self.warp_filter.SetInput(self.grid_mapper.GetUnstructuredGridOutput())

            #self.geom_filter = vtkGeometryFilter()
            #self.geom_filter.SetInput(self.warp_filter.GetUnstructuredGridOutput())

            #self.geom_mapper = vtkPolyDataMapper()
            #self.geom_actor.setMapper(self.geom_mapper)

        #if 0:
            #from vtk.numpy_interface import algorithms
            #arrow = vtkArrowSource()
            #arrow.PickableOff()

            #self.glyph_transform = vtkTransform()
            #self.glyph_transform_filter = vtkTransformPolyDataFilter()
            #self.glyph_transform_filter.SetInputConnection(arrow.GetOutputPort())
            #self.glyph_transform_filter.SetTransform(self.glyph_transform)

            #self.glyph = vtkGlyph3D()
            #self.glyph.setInput(xxx)
            #self.glyph.SetSource(self.glyph_transform_filter.GetOutput())

            #self.glyph.SetVectorModeToUseVector()
            #self.glyph.SetColorModeToColorByVector()
            #self.glyph.SetScaleModeToScaleByVector()
            #self.glyph.SetScaleFactor(1.0)

            #self.append_filter = vtkAppendFilter()
            #self.append_filter.AddInputConnection(self.grid.GetOutput())


        #self.warpVector = vtkWarpVector()
        #self.warpVector.SetInput(self.grid_mapper.GetUnstructuredGridOutput())
        #grid_mapper.SetInput(Filter.GetOutput())

        self.geom_actor = vtkLODActor()
        self.geom_actor.DragableOff()
        self.geom_actor.SetMapper(self.grid_mapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        #self.geom_actor.GetProperty().SetDiffuseColor(1, 0, 0)  # red

        #if 0:
            #id_filter = vtkIdFilter()

            #ids = np.array([1, 2, 3], dtype='int32')
            #id_array = numpy_to_vtk(
                #num_array=ids,
                #deep=True,
                #array_type=VTK_INT,
            #)

            #id_filter.SetCellIds(id_array.GetOutputPort())
            #id_filter.CellIdsOff()
            #self.grid_mapper.SetInputConnection(id_filter.GetOutputPort())
        self.rend.AddActor(self.geom_actor)
        self.build_glyph()

    def build_glyph(self) -> None:
        """builds the glyph actor"""
        glyph_source, glyphs, glyph_mapper, arrow_actor = build_glyph(self.grid)
        self.rend.AddActor(arrow_actor)

        self.glyph_source = glyph_source
        self.glyphs = glyphs
        self.glyph_mapper = glyph_mapper
        self.arrow_actor = arrow_actor
        #-----------------------------------------
        glyphs_centroid = vtkGlyph3D()
        glyphs_centroid.SetVectorModeToUseVector()
        glyphs_centroid.SetScaleModeToScaleByVector()
        glyphs_centroid.SetColorModeToColorByScale()
        glyphs_centroid.ScalingOn()
        glyphs_centroid.ClampingOn()
        glyphs_centroid.SetSourceConnection(glyph_source.GetOutputPort())

        glyph_mapper_centroid = vtkPolyDataMapper()
        glyph_mapper_centroid.SetInputConnection(glyphs_centroid.GetOutputPort())
        glyph_mapper_centroid.ScalarVisibilityOff()

        arrow_actor_centroid = vtkLODActor()
        arrow_actor_centroid.SetMapper(glyph_mapper_centroid)

        self.glyphs_centroid = glyphs_centroid
        self.glyph_mapper_centroid = glyph_mapper_centroid
        self.arrow_actor_centroid = arrow_actor_centroid

    def _build_vtk_frame_post(self, build_lookup_table: bool=True) -> None:
        if build_lookup_table:
            self.build_lookup_table()

        corner_text_size = self.settings.corner_text_size # was 14
        dtext_size = corner_text_size + 1

        # we build these in reverse order
        create_text = self.tool_actions.create_text
        create_text([5, 5 + 3 * dtext_size], 'Max  ', text_size=corner_text_size)  # text actor 0
        create_text([5, 5 + 2 * dtext_size], 'Min  ', text_size=corner_text_size)  # text actor 1
        create_text([5, 5 + 1 * dtext_size], 'Word1', text_size=corner_text_size)  # text actor 2
        create_text([5, 5], 'Word2', text_size=corner_text_size)  # text actor 3

        self.get_edges()
        if self.settings.is_edges_visible:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOn()
        else:
            prop = self.edge_actor.GetProperty()
            prop.EdgeVisibilityOff()

    def get_edges(self) -> None:
        """Create the edge actor"""
        edge_mapper = self.edge_mapper
        edge_actor = self.edge_actor
        color_function = (
            self.color_function_black if self.settings.is_edges_black
            else self.color_function)
        create_edges_from_grid(
            self.grid_selected, edge_mapper, edge_actor,
            color_function,
            is_edges_visible=self.settings.is_edges_visible)
        self.rend.AddActor(edge_actor)

    #---------------------------------------------------------------------------
    # properties

    @property
    def logo(self) -> str:
        """Gets the pyNastran icon path, which can be overwritten"""
        return self._logo

    @logo.setter
    def logo(self, logo: str) -> None:
        """Sets the pyNastran icon path, which can be overwritten"""
        self._logo = logo

    @property
    def legend_shown(self) -> bool:
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
    def result_name(self) -> str:
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
    def on_show_debug(self) -> bool:
        """sets a flag for showing/hiding DEBUG messages"""
        self.settings.show_debug = not self.settings.show_debug

    def on_show_info(self) -> bool:
        """sets a flag for showing/hiding INFO messages"""
        self.settings.show_info = not self.settings.show_info

    def on_show_command(self) -> bool:
        """sets a flag for showing/hiding COMMAND messages"""
        self.settings.show_command = not self.settings.show_command

    def on_show_warning(self) -> bool:
        """sets a flag for showing/hiding WARNING messages"""
        self.settings.show_warning = not self.settings.show_warning

    def on_show_error(self) -> bool:
        """sets a flag for showing/hiding ERROR messages"""
        self.settings.show_error = not self.settings.show_error

    def hide_axes(self, cids=None) -> None:
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

    def show_axes(self, cids=None) -> None:
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

    def delete_actor(self, name: str) -> None:
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

    def show_only(self, names: list[str]) -> None:
        """
        Show these actors only

        names : str, list[str]
            names to show
            If they're hidden, show them.
            If they're shown and shouldn't be, hide them.

        ..todo :: update the GeomeryProperties
        """
        raise NotImplementedError('show_only')

    def hide_actors(self, except_names=None, hide_legend: bool=True) -> None:
        """
        Hide all the actors

        except_names : str, list[str], None
            list of names to exclude
            None : hide all

        ..note :: If an actor is hidden and in the except_names, it will still be hidden.
        ..todo :: update the GeomeryProperties
        """
        if except_names is None:
            except_names = []
        elif isinstance(except_names, str):
            except_names = [except_names]

        # hide everything but the main grid
        for key, actor in self.geometry_actors.items():
            if key not in except_names:
                actor.VisibilityOff()

        self.hide_axes()
        if hide_legend:
            self.hide_legend()
        #self.settings.set_background_color_to_white()

    #---------------------------------------------------------------------------
    def Render(self) -> None:
        """Renders the GUI"""
        #self.vtk_interactor.Render()
        self.vtk_interactor.GetRenderWindow().Render()

    def cell_centroid_grid(self, grid, cell_id: int,
                           dtype: int='float32') -> Optional[np.ndarray]:
        """gets the cell centroid"""
        if not self.is_gui:
            centroid = np.zeros(3, dtype=dtype)
            return centroid
        cell = grid.GetCell(cell_id)
        ncells = grid.GetNumberOfCells()
        if cell_id > ncells:
            return None

        try:
            nnodes = cell.GetNumberOfPoints()
        except AttributeError:
            return None
        points = cell.GetPoints()
        assert nnodes > 0, f'nnodes={nnodes:d} cell_id={cell_id:d} cell={cell}'
        centroid = np.zeros(3, dtype=dtype)
        for ipoint in range(nnodes):
            point = np.array(points.GetPoint(ipoint), dtype=dtype)
            centroid += point
        centroid /= nnodes
        return centroid

    def cell_centroid(self, cell_id: int, dtype: str='float32') -> None:
        """gets the cell centroid"""
        return self.cell_centroid_grid(self.grid_selected, cell_id, dtype=dtype)

    def build_lookup_table(self) -> None:
        """build the nominal lookup table for the scalar bar"""
        scalar_range = self.grid_selected.GetScalarRange()
        self.grid_mapper.SetScalarRange(scalar_range)
        self.grid_mapper.SetLookupTable(self.color_function)
        self.rend.AddActor(self.scalar_bar_actor)

    def _set_results(self, form, cases) -> None:
        assert len(cases) > 0, cases
        self.case_keys = list(cases.keys())
        #self.case_keys = sorted(cases.keys())
        assert isinstance(cases, dict), type(cases)

        self.model_data.result_cases = cases

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

    #===========================================================================
    # groups
    # fake functions
    def show_eids(self, eids) -> None:
        pass

    # fake functions
    #---------------
    # real functions
    def find_result_by_name(self, desired_name: str,
                            restype: str='either') -> np.ndarray:
        for icase in range(self.ncases):
            name, result = self.get_name_result_data(icase, restype=restype)
            if name == desired_name:
                return result
        raise RuntimeError('cannot find name=%r' % desired_name)

    def post_group_by_name(self, name: str) -> None:
        """posts a group with a specific name"""
        assert isinstance(name, str), name
        group = self.groups[name]
        self.post_group(group, update_groups=False)

    def post_group(self, group, update_groups=False):
        """posts a group object"""
        name = group.name
        if update_groups and name not in self.groups:
            self.groups[name] = group
        eids = group.element_ids
        self.show_eids(eids)
        assert isinstance(name, str), name
        self.group_active = name

    def create_group_with_name(self, name: str,
                               eids: np.ndarray,
                               model_name: str='main') -> Group:
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

        Returns
        -------
        group : Group
            the corresponding group

        """
        elements_pound = self.groups[model_name].elements_pound
        element_str = ''
        group = Group(
            name, element_str, elements_pound,
            editable=True)

        # TODO: make sure all the eids exist
        group.element_ids = eids
        self.log_command('self.create_group_with_name(%r, %r)' % (name, eids))
        self.groups[name] = group
        return group

    def map_element_centroid_to_node_fringe_result(self, update_limits: bool=True,
                                                   show_msg: bool=True) -> bool:
        """
        Maps elemental fringe results to nodal fringe results.

        If you have a 5-noded CQUAD4 (e.g., centroid + 4 nodes), only the
        centroidal value will be mapped, even though you could map the
        average the nodal values instead.  It's not wrong to do it this
        way, but it could be more accurate.

        If you have CQUAD4 with centroidal only or something like strain
        energy, this will map properly.

        >>> is_passed = self.map_element_centroid_to_node_fringe_result(update_limits=True)
        """
        is_passed = False
        icase_fringe = self.icase_fringe
        if icase_fringe is None:
            self.log.error('No fringe result shown')
            return is_passed

        (obj, (i, name)) = self.result_cases[icase_fringe]
        location = obj.get_location(i, name)
        is_passed, out_data = map_element_centroid_to_node_fringe_result(
            self.grid, location, self.log)
        if not is_passed:
            return is_passed
        imin, imax, min_value_actual, max_value_actual = out_data

        legend_title = obj.get_legend_title(i, name)
        min_value, max_value = obj.get_min_max(i, name)
        if update_limits:
            min_value = min_value_actual
            max_value = max_value_actual

        data_format = obj.get_data_format(i, name)
        nlabels, labelsize, ncolors, colormap = obj.get_nlabels_labelsize_ncolors_colormap(i, name)
        is_legend_shown = self.scalar_bar.is_shown

        self.update_scalar_bar(legend_title, min_value, max_value,
                               data_format,
                               nlabels=nlabels, labelsize=labelsize,
                               ncolors=ncolors, colormap=colormap,
                               is_shown=is_legend_shown)

        location_nodal = 'node'
        self._update_min_max_actors(location_nodal, icase_fringe,
                                    imin, min_value_actual,
                                    imax, max_value_actual)

        subcase_id = obj.subcase_id
        subtitle, label = self.get_subtitle_label(subcase_id)
        label2 = obj.get_annotation(i, name)
        if label2:
            label += '; ' + label2

        self.tool_actions.update_corner_text_actors(
            location=location_nodal,
            subcase_id=subcase_id,
            subtitle=subtitle,
            label=label,
            imin=imin, min_value=min_value,
            imax=imax, max_value=max_value,
        )
        self.vtk_interactor.Render()
        self.log_command(f'self.map_element_centroid_to_node_fringe_result('
                         f'update_limits={update_limits}, show_msg={show_msg})')
        return is_passed
