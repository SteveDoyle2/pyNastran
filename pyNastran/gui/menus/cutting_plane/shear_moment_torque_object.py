"""
defines:
 - ShearMomentTorqueObject

"""
from __future__ import annotations
import os
from typing import cast, Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    get_element_centroids)

from pyNastran.gui.menus.groups_modify.groups import Group
from pyNastran.gui.menus.cutting_plane.shear_moment_torque import ShearMomentTorqueWindow, ResultsDialog
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
    from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
    from pyNastran.gui.qt_files.vtk_actor_actions import VtkActorActions
    from vtk import vtkActor, vtkProperty

PLANE_NAME = 'smt_plane'
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

        self.model_name = ''
        self.model_data = {}

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the preferences window"""
        if self._smt_window_shown:
            self._smt_window.set_font_size(font_size)

    def get_gpforce_from_icase(self, icase: int) -> tuple[bool, Optional[RealGridPointForcesArray]]:
        is_failed = True
        gui = self.gui
        if icase == -1 or len(gui.result_cases) == 0:
            gui.log.error('Select a Grid Point Forces result.')
            return is_failed, None

        (obj, (unused_i, unused_name)) = gui.result_cases[icase]
        if not hasattr(obj, 'gpforce_array'):
            gui.log.error('Select a Grid Point Forces result.')
            return is_failed, None
        is_failed = False
        gpforce: RealGridPointForcesArray = obj.gpforce_array
        return is_failed, gpforce

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
        model = self.setup_model_data(model_name)
        if hasattr(model, 'coords'):
            cids = list(model.coords.keys())
            cids.sort()
        else:
            cids = [0]

        icase = self.gui.icase
        is_failed, gpforce = self.get_gpforce_from_icase(icase)
        if is_failed:
            return

        group: Group = self.gui.groups[gui.group_active]
        elements_pound = group.elements_pound
        data = {
            'font_size' : settings.font_size,
            'icase': icase,
            'cids' : cids,
            'elements_pound': elements_pound,
            'plane_color' : settings.shear_moment_torque_color,
            'plane_opacity' : settings.shear_moment_torque_opacity,
            'vector_line_width' : settings.shear_moment_torque_line_width,
            'vector_point_size' : settings.shear_moment_torque_point_size,

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

    def get_model(self, model_name: str) -> BDF:
        model = self.gui.models[model_name]
        return model

    def setup_model_data(self, model_name: str) -> BDF:
        model = self.get_model(model_name)

        if model_name != self.model_name:
            nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
            eids, element_centroids_cid0 = get_element_centroids(model)
            self.model_data = {
                'nids': nids,
                'nid_cd': nid_cd,
                'icd_transform': icd_transform,
                'xyz_cid0': xyz_cid0,
                'eids': eids,
                'element_centroids_cid0': element_centroids_cid0,
            }
            self.model_name = model_name

    def set_plane_properties(self) -> None:
        gui: MainWindow = self.gui
        settings: Settings = gui.settings

        opacity = settings.shear_moment_torque_opacity
        line_width = settings.shear_moment_torque_line_width
        color = settings.shear_moment_torque_color
        point_size = settings.shear_moment_torque_point_size

        if PLANE_NAME in gui.geometry_actors:
            actor: vtkActor = gui.geometry_actors[PLANE_NAME]
            alt_geometry: AltGeometry = gui.geometry_properties[PLANE_NAME]
            set_plane_opacity_color(
                alt_geometry, actor,
                opacity=opacity, color=color)

        if LINE_NAME in gui.geometry_actors:
            actor: vtkActor = gui.geometry_actors[LINE_NAME]
            alt_geometry: AltGeometry = gui.geometry_properties[LINE_NAME]
            set_plane_opacity_color(
                alt_geometry, actor, color=color,
                opacity=opacity, line_width=line_width)

        if POINT_NAME in gui.geometry_actors:
            actor: vtkActor = gui.geometry_actors[POINT_NAME]
            alt_geometry = gui.geometry_properties[POINT_NAME]
            set_plane_opacity_color(
                alt_geometry, actor, color=color,
                opacity=opacity, point_size=point_size)

        self.render()

    def on_clear_plane_actors(self, render: bool=True) -> None:
        gui: MainWindow = self.gui
        names = [
            PLANE_NAME, LINE_NAME, POINT_NAME, ARROW_NAME,
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

    def make_plane_from_data(self, data: dict[str, Any]) -> None:
        """the callback to validate the plane/coordinate system"""
        #gpforce: RealGridPointForcesArray = data['gpforce']
        nplanes = data['nplanes']
        #model = self.models[model_name]

        cid_p1, p1 = data['p1'] # start
        cid_p2, p2 = data['p2'] # xzplane
        cid_p3, p3 = data['p3'] # end
        cid_zaxis, zaxis = data['zaxis']
        method = data['method']

        self.plot_cutting_plane(
            p1, p2, p3, zaxis,
            method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes,
        )
        return

    def make_smt_from_data(self, data: dict[str, Any], show: bool=True) -> None:
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

        #gpforce: RealGridPointForcesArray = data['gpforce']
        icase: int = data['icase']
        nplanes = data['nplanes']

        cid_p1, p1 = data['p1'] # start
        cid_p2, p2 = data['p2'] # xzplane
        cid_p3, p3 = data['p3'] # end
        cid_zaxis, zaxis = data['zaxis']
        method = data['method']
        station_location = data['station_location']
        element_ids = data['element_ids']

        length_scale, length_unit = data['length']
        force_scale, force_unit = data['force']
        moment_scale, moment_unit = data['moment']
        csv_filename = None
        if 'csv_filename' in data:
            csv_filename = data['csv_filename']
        self.plot_shear_moment_torque(
            icase,
            p1, p2, p3, zaxis,
            method=method,
            station_location=station_location,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes,
            element_ids=element_ids,
            csv_filename=csv_filename,
            length_scale=length_scale, length_unit=length_unit,
            force_scale=force_scale, force_unit=force_unit,
            moment_scale=moment_scale, moment_unit=moment_unit,
            show=show)
        return

    def plot_cutting_plane(self,
                           p1: np.ndarray,
                           p2: np.ndarray,
                           p3: np.ndarray,
                           zaxis: np.ndarray,
                           method: str='Vector',
                           cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                           nplanes: int=20,
                           stop_on_failure: bool=False) -> None:
        """Plane Actor is drawn in the i-k plane"""
        #model = self.gui.models[self.model_name]
        #nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
        #self.model_data = {
            #'nids': nids,
            #'nid_cd': nid_cd,
            #'icd_transform': icd_transform,
            #'xyz_cid0': xyz_cid0,
        #}
        model = self.get_model(self.model_name)
        xyz_cid0 = self.model_data['xyz_cid0']

        #nplanes = 1
        unused_out = self.plot_plane(
            model, xyz_cid0,
            p1, p2, p3,
            zaxis,
            method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes,
            stop_on_failure=stop_on_failure)
        return

    def plot_shear_moment_torque(self,
                                 icase: int,
                                 p1: np.ndarray,
                                 p2: np.ndarray,
                                 p3: np.ndarray,
                                 zaxis: np.ndarray,
                                 method: str='Vector',
                                 station_location: str='End-Origin',
                                 cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                                 nplanes: int=20,
                                 element_ids: Optional[np.ndarray]=None,
                                 csv_filename=None,
                                 length_scale: float=1.0, length_unit: str='',
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
        icase : int
            a result case with grid point forces
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
            'coord id':
               zaxis:  N/A
               p2:     use the coord id from p2 as the output frame
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
            this will be empty if there no case is found

        """
        is_failed, gpforce = self.get_gpforce_from_icase(icase)
        if is_failed:
            return np.array([]), np.array([])

        gui: MainWindow = self.gui
        log = gui.log

        model = self.get_model(self.model_name)
        nids = self.model_data['nids']
        nid_cd = self.model_data['nid_cd']
        icd_transform = self.model_data['icd_transform']
        xyz_cid0 = self.model_data['xyz_cid0']

        eids = self.model_data['eids']
        element_centroids_cid0 = self.model_data['element_centroids_cid0']

        if element_ids is not None:
            element_ids = np.unique(element_ids)
            ieids = np.searchsorted(eids, element_ids)
            assert np.array_equal(eids[ieids], element_ids)

            eids = element_ids
            element_centroids_cid0 = element_centroids_cid0[ieids, :]

        is_failed, stations, coord, iaxis_march = self.plot_plane(
            model, xyz_cid0,
            p1, p2, p3,
            zaxis, method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            nplanes=nplanes,
            stop_on_failure=stop_on_failure)
        if is_failed:
            force_sum = np.zeros((0, 3))
            moment_sum = np.zeros((0, 3))
            return force_sum, moment_sum

        coord = cast(CORD2R, coord)
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

        cids, origins, xyz_stations, xlabel = get_xyz_stations(
            stations*length_scale, new_coords,
            station_location=station_location)
        origins *= length_scale

        length_unit2 = f' ({length_unit})' if length_unit else ''
        force_unit2 = f' ({force_unit})' if force_unit else ''
        moment_unit2 = f' ({moment_unit})' if moment_unit else ''
        labels = [
            f'Station{length_unit2}',
            f'X{length_unit2}', f'Y{length_unit2}', f'Z{length_unit2}',
            f'Fx{force_unit2}', f'Fy{force_unit2}', f'Fz{force_unit2}',
            f'Mx{moment_unit2}', f'My{moment_unit2}', f'Mz{moment_unit2}'
        ]
        data = np.column_stack([xyz_stations, origins, force_sum, moment_sum])

        if self._smt_window is not None:
            # automatically detroy when the window is closed by referencing
            # an object in the window
            parent = self._smt_window.p1_label
            dlg = ResultsDialog(parent, data, labels,
                                title='Shear, Moment, Torque Results')

        if csv_filename:
            root_filename = os.path.splitext(csv_filename)[0]
            write_smt_to_csv(
                csv_filename,
                stations, nelems, nnodes, cids, origins,
                force_sum, moment_sum,
                length_unit=length_unit,
                force_unit=force_unit,
                moment_unit=moment_unit)
        plot_smt(
            xyz_stations,
            force_sum, moment_sum,
            nelems, nnodes, show=show,
            xtitle=xlabel, xlabel=xlabel,
            length_unit=length_unit,
            force_unit=force_unit,
            moment_unit=moment_unit,
            root_filename=root_filename,
            plot_force_components=False,
            plot_moment_components=False,
        )
        msg = (
            f'model_name = {self.model_name!r}\n'
            f'model = self.shear_moment_torque_obj.setup_model_data(model_name)\n'
            f'icase = {icase:d}\n'
            f'p1 = np.array([{p1[0]}, {p1[1]}, {p1[2]}])\n'
            f'p2 = np.array([{p2[0]}, {p2[1]}, {p2[2]}])\n'
            f'p3 = np.array([{p3[0]}, {p3[1]}, {p3[2]}])\n'
            f'zaxis = np.array([{zaxis[0]}, {zaxis[1]}, {zaxis[2]}])\n'
            f'force_sum, moment_sum = self.shear_moment_torque_obj.plot_shear_moment_torque(\n'
            f'    icase, p1, p2, p3, zaxis, method={method!r},\n'
            f'    station_location={station_location!r},\n'
            f'    cid_p1={cid_p1:d}, cid_p2={cid_p2:d}, cid_p3={cid_p3:d}, cid_zaxis={cid_zaxis:d},\n'
            f'    nplanes={nplanes:d}, csv_filename={csv_filename!r},\n'
            f'    length_scale={length_scale}, length_unit={length_unit!r},\n'
            f'    force_scale={force_scale}, force_unit={force_unit!r},\n'
            f'    moment_scale={moment_scale}, moment_unit={moment_unit!r},\n'
            f'    show={show}, stop_on_failure={stop_on_failure})'
        )
        gui.log_command(msg)

        return force_sum, moment_sum

    def plot_plane(self,
                   model: BDF,
                   xyz_cid0: np.ndarray,
                   p1: np.ndarray,
                   p2: np.ndarray,
                   p3: np.ndarray,
                   zaxis: np.ndarray,
                   method: str='Vector',
                   cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                   nplanes: int=20,
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
             - origin + xz-plane
        method : str; default='Vector'
            'Vector':
               zaxis:  k vector
               p2:     xz-plane vector
            'CORD2R':
               zaxis: point on the z-axis
               p2:     point on the xz-plane
            'Coord ID':
               zaxis: point on the z-axis
               p2:    coord id defines the axes
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

        gui: MainWindow = self.gui
        settings: Settings = gui.settings
        plane_color = settings.shear_moment_torque_color
        plane_opacity = settings.shear_moment_torque_opacity
        point_size = settings.shear_moment_torque_point_size
        line_width = settings.shear_moment_torque_line_width

        # the plane actor defines the plane of the output results,
        # not the plane of the march direction
        # xyz1: origin
        # xyz2: xzplane
        unused_plane_actor = self._create_plane_actor(
            xyz1, xyz2,
            coord, i, k, dim_max,
            plane_color, plane_opacity,
        )
        self._create_vector_points(
            xyz1, iaxis_march, stations,
            plane_color, plane_opacity,
            point_size, line_width,
        )

        is_failed = False
        self.render()
        return is_failed, stations, coord, iaxis_march

    def _create_vector_points(self, xyz1: np.ndarray,
                              iaxis_march: np.ndarray,
                              stations: np.ndarray,
                              color: ColorFloat,
                              opacity: float,
                              point_size: float,
                              line_width: float) -> None:
        """
        Creates additional actors to represent the location of
        the cutting planes
        """
        gui = self.gui
        nodes = xyz1[np.newaxis, :] + iaxis_march[np.newaxis, :] * stations[:, np.newaxis]
        irange = np.arange(len(stations)-1, dtype='int32')
        elements = np.column_stack([irange, irange+1])

        vtk_actor_actions: VtkActorActions = gui.vtk_actor_actions

        # update the settings
        for actor_name in (LINE_NAME, POINT_NAME):
            if actor_name in gui.geometry_properties:
                alt_geometry: AltGeometry = gui.geometry_properties[actor_name]
                alt_geometry.color = color
                alt_geometry.opacity = opacity
                alt_geometry.line_width = line_width
                alt_geometry.point_size = point_size

        # plot line and points to indicate plane locations
        unused_ugrid1 = vtk_actor_actions.set_line_grid(
            LINE_NAME, nodes, elements,
            color, line_width=line_width, opacity=opacity,
            representation='wire', add=True,
            visible_in_geometry_properties=False)
        unused_ugrid2 = vtk_actor_actions.set_line_grid(
            POINT_NAME, nodes, elements,
            color, point_size=point_size, opacity=opacity,
            representation='point', add=True,
            visible_in_geometry_properties=False)

        self.set_plane_properties()
        gui.geometry_actors[LINE_NAME]
        gui.geometry_actors[POINT_NAME]
        #if ugrid2 is not None:
            #glyph_source, glyphs, glyph_mapper, arrow_actor = build_glyph(ugrid2)
            #gui.alt_grids[ARROW_NAME] = arrow_actor
            #gui.rend.AddActor(arrow_actor)

    def render(self) -> None:
        self.gui.rend.Render()
        self.gui.Render()

    def _create_plane_actor(self,
                            xyz1: np.ndarray,
                            xyz2: np.ndarray,
                            coord: CORD2R,
                            i: np.ndarray,
                            k: np.ndarray,
                            dim_max: float,
                            color: ColorFloat,
                            opacity: float,
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
        vtk_actor_actions: VtkActorActions = gui.vtk_actor_actions

        # TODO: make coord show up in EditGeometryProperties
        #       label !='' makes it show up, but it messes up the coord text
        gui.tool_actions.create_coordinate_system(
            cid, dim_max, label='', origin=origin,
            matrix_3x3=beta, coord_type='xyz')

        j = np.cross(k, i)
        #center = (xyz1 + xyz2) / 2.

        actor_name = PLANE_NAME
        if actor_name in gui.geometry_properties:
            alt_geometry: AltGeometry = gui.geometry_properties[actor_name]
            alt_geometry.color = color
            alt_geometry.opacity = opacity
            #alt_geometry.color = color
            #alt_geometry.color = color

        plane_actor = vtk_actor_actions.create_plane_actor_from_points(
            xyz1, j, k, dim_max,
            representation='surface',
            color=color,  # floats
            opacity=opacity, # 0=transparent, 1=solid
            actor_name=actor_name,
            visible_in_geometry_properties=False,
        )
        #alt_geometry = gui.geometry_properties[PLANE_NAME]
        #set_plane_opacity_color(alt_geometry, plane_actor,
                                #opacity=opacity, color=color)
        #plane_actor.VisibilityOn()
        return plane_actor

def set_plane_opacity_color(alt_geom: AltGeometry,
                            actor: vtkActor,
                            color: ColorFloat,
                            opacity: Optional[float]=None,
                            point_size: Optional[float]=None,
                            line_width: Optional[float]=None,
                            ) -> None:
    """sets the opacity and color of the plane"""
    alt_geom.set_color(color)
    prop: vtkProperty = actor.GetProperty()
    prop.SetColor(*color)  # floats
    prop.SetDiffuseColor(*color)
    if opacity is not None:
        alt_geom.opacity = opacity
        prop.SetOpacity(opacity) # 0=transparent, 1=solid
    if line_width is not None:
        alt_geom.line_width = line_width
        prop.SetLineWidth(line_width)
    if point_size is not None:
        alt_geom.point_size = point_size
        prop.SetPointSize(point_size)
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
