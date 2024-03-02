"""
defines:
 - ShearMomentTorqueObject

"""
from __future__ import annotations
import os
from typing import cast, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    get_element_centroids)

from pyNastran.gui.menus.cutting_plane.shear_moment_torque import ShearMomentTorqueWindow, ResultsDialog
from pyNastran.gui.qt_files.colors import PURPLE_FLOAT
from pyNastran.gui.qt_files.base_gui import BaseGui

from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.op2.tables.ogf_gridPointForces.smt import (
    get_nid_cd_xyz_cid0, setup_coord_from_plane,
    get_xyz_stations, write_smt_to_csv, plot_smt,)

#from pyNastran.gui.utils.vtk.vtk_vector import build_glyph
from pyNastran.gui.typing import ColorFloat
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.gui_objects.settings import Settings
    from vtk import vtkActor

LINE_NAME = 'smt_vector'
POINT_NAME = 'smt_points'
ARROW_NAME = 'smt_arrow'

class ShearMomentTorqueObject(BaseGui):
    """wrapper around ShearMomentTorqueWindow"""
    def __init__(self, gui: MainWindow):
        #self.gui = gui
        super().__init__(gui)
        self._smt_shown = False
        self._smt_window = None
        self._smt_window_shown = False

    def set_font_size(self, font_size):
        """sets the font size for the preferences window"""
        if self._smt_window_shown:
            self._smt_window.set_font_size(font_size)

    def set_shear_moment_torque_menu(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        """
        #if not hasattr(self, 'case_keys'):  # TODO: maybe include...
            #self.log_error('No model has been loaded.')
            #return

        #camera = self.gui.GetCamera()
        #min_clip, max_clip = camera.GetClippingRange()
        gui = self.gui
        settings: Settings = gui.settings
        model_name = gui.name
        model = gui.models[model_name]
        if hasattr(model, 'coords'):
            cids = list(model.coords.keys())
            cids.sort()
        else:
            cids = [0]

        icase = self.gui.icase
        if icase == -1 or len(gui.result_cases) == 0:
            gui.log.error('Select a Grid Point Forces result.')
            return

        (obj, (unused_i, unused_name)) = gui.result_cases[icase]
        if not hasattr(obj, 'gpforce_array'):
            gui.log.error('Select a Grid Point Forces result.')
            return

        gpforce: RealGridPointForcesArray = obj.gpforce_array
        data = {
            'font_size' : settings.font_size,
            'cids' : cids,
            'plane_color' : (1., 0., 1.), # purple
            'plane_opacity' : 0.9,

            'gpforce' : gpforce,
            'model_name' : model_name,
            'clicked_ok' : False,
            'close' : False,
        }
        if not self._smt_window_shown:
            self._smt_window = ShearMomentTorqueWindow(data, win_parent=gui)
            self._smt_window.show()
            self._smt_window_shown = True
            self._smt_window.exec_()
        else:
            self._smt_window.activateWindow()

        if data['close']:
            #if not self._smte_window._updated_preference:
            #    settings.on_set_font_size(data['font_size'])
            del self._smt_window
            self.on_clear_plane_actors(render=True)
            self._smt_window_shown = False
        else:
            self._smt_window.activateWindow()

    def set_plane_properties(self, opacity: float, color: ColorFloat) -> None:
        gui: MainWindow = self.gui
        if hasattr(gui, 'plane_actor'):
            plane_actor: vtkActor = gui.plane_actor
            self.set_plane_opacity_color(plane_actor, opacity, color)
            if plane_actor.GetVisibility() == 0:
                plane_actor.VisibilityOn()
                plane_actor.Modified()
        self.render()

    def on_clear_plane_actors(self, render: bool=True) -> None:
        gui: MainWindow = self.gui
        #if hasattr(gui, 'plane_actor'):
            #del gui.plane_actor
        #gui.clear_actor('smt_plane')
        names = [
            'plane_actor', 'smt_plane',
            LINE_NAME, POINT_NAME, ARROW_NAME,
        ]
        clear_actors(gui, names)

        # clear the coordinate system
        cid = 'smt_axes'
        if cid in gui.axes:
            actor = gui.axes[cid]
            gui.rend.RemoveActor(actor)
            del gui.axes[cid]
        #gui.geometry_properties[label] = CoordProperties(
            #label, coord_type, is_visible, coord_scale)
        #gui.geometry_actors[label] = axes

        if render:
            self.render()

    def make_plane_from_data(self, data) -> None:
        """the callback to validate the plane/coordinate system"""
        model_name = data['model_name']
        #gpforce: RealGridPointForcesArray = data['gpforce']
        nplanes = data['nplanes']
        #model = self.models[model_name]

        cid_p1, p1 = data['p1'] # start
        cid_p2, p2 = data['p2'] # xzplane
        cid_p3, p3 = data['p3'] # end
        cid_zaxis, zaxis = data['zaxis']
        plane_color = data['plane_color']
        plane_opacity = data['plane_opacity']
        method = data['method']

        self.plot_cutting_plane(
            model_name,
            p1, p2, p3, zaxis,
            method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            plane_color=plane_color, plane_opacity=plane_opacity, nplanes=nplanes,
        )
        return

    def make_smt_from_data(self, data, show: bool=True) -> None:
        """the callback for the shear, moment, torque plotter"""
        #self.out_data['method'] = method
        #self.out_data['p1'] = [p1_cid, p1]
        #self.out_data['p2'] = [p2_cid, p2]
        #self.out_data['p3'] = [p2_cid, p3]
        #self.out_data['zaxis'] = [zaxis_method, zaxis_cid, zaxis]
        #self.out_data['ytol'] = ytol
        #self.out_data['zero_tol'] = zero_tol
        #self.out_data['plane_color'] = self.plane_color_float
        #self.out_data['plane_opacity'] = 0.6
        #self.out_data['csv_filename'] = csv_filename
        #self.out_data['clicked_ok'] = True

        model_name = data['model_name']
        gpforce: RealGridPointForcesArray = data['gpforce']
        nplanes = data['nplanes']
        #model = self.models[model_name]

        cid_p1, p1 = data['p1'] # start
        cid_p2, p2 = data['p2'] # xzplane
        cid_p3, p3 = data['p3'] # end
        cid_zaxis, zaxis = data['zaxis']
        plane_color = data['plane_color']
        plane_opacity = data['plane_opacity']
        method = data['method']
        station_location = data['station_location']
        force_scale, force_unit = data['force']
        moment_scale, moment_unit = data['moment']

        csv_filename = None
        if 'csv_filename' in data:
            csv_filename = data['csv_filename']
        self.plot_shear_moment_torque(
            model_name, gpforce,
            p1, p2, p3, zaxis,
            method=method,
            station_location=station_location,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes, plane_color=plane_color, plane_opacity=plane_opacity,
            csv_filename=csv_filename,
            force_scale=force_scale, force_unit=force_unit,
            moment_scale=moment_scale, moment_unit=moment_unit,
            show=show)
        return

    def plot_cutting_plane(self,
                           model_name: str,
                           p1: np.ndarray,
                           p2: np.ndarray,
                           p3: np.ndarray,
                           zaxis: np.ndarray,
                           method: str='Vector',
                           cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                           nplanes: int=20, plane_color: Optional[ColorFloat]=None,
                           plane_opacity: float=0.5,
                           stop_on_failure: bool=False) -> None:
        """Plane Actor is drawn in the i-k plane"""
        if plane_color is None:
            plane_color = PURPLE_FLOAT
        assert len(plane_color) == 3, plane_color

        model = self.gui.models[model_name]
        nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
        #xyz1, xyz2, xyz3, i, k, origin, xzplane, dim_max, stations

        #nplanes = 1
        unused_out = self.plot_plane(
            model, xyz_cid0, plane_color,
            p1, p2, p3,
            zaxis,
            method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes, plane_opacity=plane_opacity,
            stop_on_failure=stop_on_failure)
        return

    def plot_shear_moment_torque(self,
                                 model_name: str,
                                 gpforce: RealGridPointForcesArray,
                                 p1: np.ndarray,
                                 p2: np.ndarray,
                                 p3: np.ndarray,
                                 zaxis: np.ndarray,
                                 method: str='Vector',
                                 station_location: str='End-Origin',
                                 cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                                 nplanes: int=20, plane_color: Optional[ColorFloat]=None,
                                 plane_opacity: float=0.5,
                                 csv_filename=None,
                                 force_scale: float=1.0, force_unit: str='',
                                 moment_scale: float=1.0, moment_unit: str='',
                                 show: bool=True,
                                 stop_on_failure: bool=False,
                                 ) -> tuple[np.ndarray, np.ndarray]:
        """
        Creates a shear moment torque plot for the active plot result
        Plane Actor is drawn in the i-k plane

        Parameters
        ----------
        model_name : str
            the name of the model
        p1: (3,) float ndarray
            defines the starting point for the shear, moment, torque plot
        p3: (3,) float ndarray
            defines the end point for the shear, moment, torque plot
        p2: (3,) float ndarray
            defines the XZ plane for the shears/moments
        zaxis: (3,) float ndarray
            the direction of the z-axis
        station_location : str; default='End-Origin'
            'End-Origin': magnitude along p3-p1
            'X' : Global X location
            'Y' : Global Y location
            'Z' : Global Z location
        cid_p1 / cid_p2 / cid_p3; default=0
            the coordinate systems for p1, p2, and p3
        method : str
            'CORD2R':
               zaxis:  point on the z-axis
               p2:     point on the xz-plane
            'Vector':
               zaxis:  k vector
               p2:     xz-plane vector
             'Z-Axis Projection':
               zaxis:  point on the z-axis
               p2:     p2 is a point on the xz-plane
        show: bool; default=True
            shows the plots

        Returns
        -------
        force_sum / moment_sum : (nstations, 3) float ndarray
            the forces/moments at the station

        """
        if plane_color is None:
            plane_color = PURPLE_FLOAT
        plane_color = cast(ColorFloat, plane_color)
        assert len(plane_color) == 3, plane_color

        gui: MainWindow = self.gui
        log = gui.log

        model = gui.models[model_name]
        nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)

        is_failed, stations, coord, iaxis_march = self.plot_plane(
            model, xyz_cid0, plane_color,
            p1, p2, p3,
            zaxis, method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes, plane_opacity=plane_opacity,
            stop_on_failure=stop_on_failure)
        if is_failed:
            force_sum = np.zeros((0, 3))
            moment_sum = np.zeros((0, 3))
            return force_sum, moment_sum

        coord = cast(CORD2R, coord)
        eids, element_centroids_cid0 = get_element_centroids(model)
        force_sum, moment_sum, new_coords, nelems, nnodes = gpforce.shear_moment_diagram(
            nids, xyz_cid0, nid_cd, icd_transform,
            eids, element_centroids_cid0,
            stations, model.coords, coord,
            iaxis_march=iaxis_march,
            itime=0, idir=0,
            nodes_tol=None, debug=False, log=log)
        self.render()

        force_sum *= force_scale
        moment_sum *= moment_scale
        root_filename = os.path.join(gui.last_dir, 'shear_moment_torque')

        origins, xyz_stations, xlabel = get_xyz_stations(
            stations, new_coords,
            station_location=station_location)

        length_unit = ''
        length_unit2 = f' ({length_unit})' if force_unit else ''
        force_unit2 = f' ({force_unit})' if force_unit else ''
        moment_unit2 = f' ({moment_unit})' if moment_unit else ''
        labels = [
            f'Station{length_unit2}',
            f'X{length_unit2}', f'Y{length_unit2}', f'Z{length_unit2}',
            f'Fx{force_unit2}', f'Fy{force_unit2}', f'Fz{force_unit2}',
            f'Mx{moment_unit2}', f'My{moment_unit2}', f'Mz{moment_unit2}'
        ]
        data = np.column_stack([xyz_stations, origins, force_sum, moment_sum])
        dlg = ResultsDialog(self._smt_window.p1_label, data, labels,
                            title='Shear, Moment, Torque Results')

        if csv_filename:
            root_filename = os.path.splitext(csv_filename)[0]
            write_smt_to_csv(
                csv_filename,
                stations, nelems, nnodes, new_coords,
                force_sum, moment_sum,
                force_unit=force_unit, moment_unit=moment_unit)
        plot_smt(
            xyz_stations,
            force_sum, moment_sum,
            nelems, nnodes, show=show,
            xtitle=xlabel, xlabel=xlabel,
            length_unit=length_unit,
            force_unit=force_unit, moment_unit=moment_unit,
            root_filename=root_filename,
            plot_force_components=False,
            plot_moment_components=False,
        )
        return force_sum, moment_sum

    def plot_plane(self,
                   model: BDF,
                   xyz_cid0: np.ndarray,
                   plane_color: ColorFloat,
                   p1: np.ndarray,
                   p2: np.ndarray,
                   p3: np.ndarray,
                   zaxis: np.ndarray,
                   method: str='Vector',
                   cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                   nplanes: int=20, plane_opacity: float=0.5,
                   show: bool=True,
                   stop_on_failure: bool=False,
                   ) -> tuple[bool, np.ndarray, CORD2R, np.ndarray]:
        """
        Parameters
        ----------
        xyz_cid0 : (nnode, 3) float array
           the nodes in the global (basic) frame
        p1 : (3,) float array
            the starting point
        p3 : (3,) float array
            the ending point
        p2 : (3,) float array
            Meaning varies depending on method
            CORD2R:
             - origin + xz-plane
            Vector:
             - xz-plane
            Z-Axis Projection:
             - see method
        method : str; default='Vector'
            'CORD2R':
               zaxis: point on the z-axis
               p2:     point on the xz-plane
            'Coord ID':
               zaxis: point on the z-axis
               p2:     coord id defines the axes
            'Vector':
               zaxis:  k vector
               p2:     xz-plane vector
             'Z-Axis Projection':
               zaxis:  point on the z-axis
               p2:     p2 is a point on the xz-plane

        Returns
        -------
        is_failed : bool
           flag to indicate success/failure
        coord_out : Coord
            the generated coordinate system where the x-axis defines
            the direction to be marched
        stations : (nplanes,) float ndarray
            the coordinates in the x-axis that will be marched down
        iaxis_march : (3,) np.ndarray
            the normalized vector of p3-p1

        """
        self.on_clear_plane_actors(render=False)
        is_failed = True
        try:
            xyz1, xyz2, xyz3, i, k, coord_out, iaxis_march, dim_max, stations = setup_coord_from_plane(
                model, xyz_cid0,
                p1, p2, p3, zaxis,
                method=method,
                cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
                nplanes=nplanes,
            )
        except Exception as verror:
            self.render()
            model.log.error(str(verror))
            if stop_on_failure:
                raise
            stations = np.array([])
            coord = CORD2R.init_from_empty()
            iaxis_march = np.array([0., 0., 0.])
            return is_failed, stations, coord, iaxis_march
        coord = coord_out

        #debug = True
        #if debug or 1:
            #print('coord_out:')
            #print(f'  origin: {coord_out.origin}')
            #print(f'  zaxis:  {coord_out.e2}')
            #print(f'  xzplane: {coord_out.e3}')
            #print(f'  method: {method!r}')

        # the plane actor defines the plane of the output results,
        # not the plane of the march direction
        # xyz1: origin
        # xyz2: xzplane
        unused_plane_actor = self.create_plane_actor(
            xyz1, xyz2,
            coord, i, k, dim_max,
            plane_color, plane_opacity,
        )

        gui = self.gui
        nodes = xyz1[np.newaxis, :] + iaxis_march[np.newaxis, :] * stations[:, np.newaxis]
        irange = np.arange(len(stations)-1, dtype='int32')
        elements = np.column_stack([irange, irange+1])
        color = plane_color

        opacity = plane_opacity
        vtk_actor_actions = gui.vtk_actor_actions

        settings: Settings = self.gui.settings
        line_width = settings.highlight_line_thickness
        point_size = settings.highlight_point_size
        line_width: int = gui._get_geometry_property_items(
            LINE_NAME,
            'line_width', 5)[0]
        point_size: int = gui._get_geometry_property_items(
            POINT_NAME,
            'point_size', 5)[0]
        #if LINE_NAME in gui.geometry_properties:
            #prop = gui.geometry_properties[LINE_NAME]
            #line_width = prop.line_width

        #point_size = 5
        #if POINT_NAME in gui.geometry_properties:
            #prop = gui.geometry_properties[POINT_NAME]
            #point_size = prop.point_size

        unused_ugrid1 = vtk_actor_actions.set_line_grid(
            LINE_NAME, nodes, elements,
            color, line_width=line_width, opacity=opacity,
            representation='wire', add=True)
        unused_ugrid2 = vtk_actor_actions.set_line_grid(
            POINT_NAME, nodes, elements,
            color, point_size=point_size, opacity=opacity,
            representation='point', add=True)

        #if ugrid2 is not None:
            #glyph_source, glyphs, glyph_mapper, arrow_actor = build_glyph(ugrid2)
            #gui.alt_grids[ARROW_NAME] = arrow_actor
            #gui.rend.AddActor(arrow_actor)

        is_failed = False
        self.render()
        return is_failed, stations, coord, iaxis_march

    def render(self) -> None:
        self.gui.rend.Render()
        self.gui.Render()

    def create_plane_actor(self,
                           xyz1: np.ndarray,
                           xyz2: np.ndarray,
                           coord: CORD2R,
                           i: np.ndarray,
                           k: np.ndarray,
                           dim_max: float,
                           plane_color: ColorFloat,
                           plane_opacity: float,
                           ) -> vtkActor:
        """
        The plane is defined in the jk plane

        Parameters
        ----------
        xyz1 / xyz2 / xyz3 : (3,) float ndarray
            the 1=starting 2=ending, 3=normal coordinates of the
            coordinate frames to create in the cid=0 frame
        i / k : (3,) float ndarray
            the i and k vectors of the coordinate system

        """
        origin = coord.origin
        beta = coord.beta().T

        cid = 'smt_axes'
        gui: MainWindow = self.gui
        vtk_actor_actions = gui.vtk_actor_actions

        # TODO: make coord show up in EditGeometryProperties
        #       label !='' makes it show up, but it messes up the coord text
        gui.tool_actions.create_coordinate_system(
            cid, dim_max, label='', origin=origin,
            matrix_3x3=beta, coord_type='xyz')

        j = np.cross(k, i)
        #center = (xyz1 + xyz2) / 2.
        plane_actor = vtk_actor_actions.create_plane_actor_from_points(
            xyz1, j, k, dim_max,
            representation='surface',
            color=plane_color,  # floats
            opacity=plane_opacity, # 0=transparent, 1=solid
            actor_name='smt_plane')
        self.set_plane_opacity_color(plane_actor, plane_opacity, plane_color)
        plane_actor.VisibilityOn()
        return plane_actor

    def set_plane_opacity_color(self, plane_actor: vtkActor,
                                opacity: float,
                                color: ColorFloat) -> None:
        """sets the opacity and color of the plane"""
        props = self.gui.geometry_properties['smt_plane']
        props.set_color(color)
        props.opacity = opacity
        prop = plane_actor.GetProperty()
        prop.SetColor(*color)  # floats
        prop.SetDiffuseColor(*color)
        prop.SetOpacity(opacity) # 0=transparent, 1=solid
        prop.Modified()

def clear_actors(gui: MainWindow, names: list[str]) -> None:
    """clears the gui actors"""
    for name in names:
        if hasattr(gui, name):
            # del gui.plane_actor
            delattr(gui, name)
        gui.clear_actor(name)
        if hasattr(gui.geometry_properties, name):
            del gui.geometry_properties[name]
