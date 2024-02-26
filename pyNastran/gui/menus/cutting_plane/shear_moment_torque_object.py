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

from pyNastran.gui.menus.cutting_plane.shear_moment_torque import ShearMomentTorqueWindow
from pyNastran.gui.qt_files.colors import PURPLE_FLOAT
from pyNastran.gui.qt_files.base_gui import BaseGui
from pyNastran.gui.typing import Color

from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.op2.tables.ogf_gridPointForces.smt import (
    get_nid_cd_xyz_cid0, plot_smt, setup_coord_from_plane)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.typing import Color
    from vtk import vtkActor


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
        settings = gui.settings
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

    def on_clear_plane_actors(self, render: bool=True) -> None:
        gui: MainWindow = self.gui
        #if hasattr(gui, 'plane_actor'):
            #del gui.plane_actor
        #gui.clear_actor('smt_plane')
        clear_actors(gui, ['plane_actor', 'smt_plane'])

        # clear the coordinate system
        cid = ''
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
        #nplanes = data['nplanes']
        #model = self.models[model_name]

        cid_p1, p1 = data['p1'] # start
        cid_p2, p2 = data['p2'] # xzplane
        cid_p3, p3 = data['p3'] # end
        cid_zaxis, zaxis = data['zaxis']
        plane_color = data['plane_color']
        plane_opacity = data['plane_opacity']
        method = data['method']
        csv_filename = data['csv_filename']

        self.plot_cutting_plane(
            model_name,
            p1, p2, p3, zaxis,
            method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
            plane_color=plane_color, plane_opacity=plane_opacity,
            csv_filename=csv_filename)
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
        force_scale, force_unit = data['force']
        moment_scale, moment_unit = data['moment']

        csv_filename = None
        if 'csv_filename' in data:
            csv_filename = data['csv_filename']
        self.plot_shear_moment_torque(
            model_name, gpforce,
            p1, p2, p3, zaxis,
            method=method,
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
                           plane_color: Optional[Color]=None, plane_opacity: float=0.5,
                           csv_filename=None,
                           stop_on_failure: bool=False) -> None:
        """Plane Actor is drawn in the i-k plane"""
        if plane_color is None:
            plane_color = PURPLE_FLOAT
        assert len(plane_color) == 3, plane_color

        model = self.gui.models[model_name]
        nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
        #xyz1, xyz2, xyz3, i, k, origin, xzplane, dim_max, stations

        nplanes = 1
        unused_out = self.plot_plane(
            model, xyz_cid0, plane_color,
            p1, p2, p3,
            zaxis, method=method,
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
                                 cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                                 nplanes: int=20, plane_color: Optional[Color]=None, plane_opacity: float=0.5,
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
        cid_p1 / cid_p2 / cid_p3
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
        plane_color = cast(Color, plane_color)
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
        if csv_filename:
            root_filename = os.path.splitext(csv_filename)[0]

            force_label = f'({force_unit})' if force_unit else ''
            moment_label = f'({moment_unit})' if moment_unit else ''
            with open(csv_filename, 'w') as csv_file:
                header = (
                    'Station,nelements,nnodes,coord_id,origin_x,origin_y,origin_z,'
                    f'Fx{force_label},Fy{force_label},Fz{force_label},'
                    f'Mx{moment_label},My{moment_label},Mz{moment_label}\n')
                csv_file.write(header)
                for station, nelem, nnode, coord_id, force_sumi, moment_sumi in zip(
                    stations, nelems, nnodes, new_coords, force_sum, moment_sum):
                    coord: CORD2R = new_coords[coord_id]
                    origin: np.ndarray = coord.origin
                    csv_file.write(
                        f'{station},{nelem:d},{nnode:d},{coord_id:d},'
                        f'{origin[0]},{origin[1]},{origin[2]},'
                        f'{force_sumi[0]},{force_sumi[1]},{force_sumi[2]},'
                        f'{moment_sumi[0]},{moment_sumi[1]},{moment_sumi[2]}\n')
        plot_smt(
            stations,
            force_sum, moment_sum,
            nelems, nnodes, show=show,
            xtitle='i Station', xlabel='i Station',
            force_unit=force_unit, moment_unit=moment_unit,
            root_filename=root_filename,
        )
        return force_sum, moment_sum

    def plot_plane(self,
                   model: BDF,
                   xyz_cid0: np.ndarray,
                   plane_color: Color,
                   p1: np.ndarray,
                   p2: np.ndarray,
                   p3: np.ndarray,
                   zaxis: np.ndarray,
                   method: str='Vector',
                   cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                   nplanes: int=20, plane_opacity: float=0.5,
                   show: bool=True,
                   stop_on_failure: bool=False,
                   ) -> tuple[bool, np.ndarray, np.ndarray, np.ndarray]:
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
            Neaning varies depending on method
            CORD2R:
             - origin + xz-plane
            Z-Axis Projection:
             - see method
        method : str; default='Vector'
            'CORD2R' :
               zaxis: point on the z-axis
               p2:     point on the xz-plane
            'Vector'
               zaxis:  k vector
               p2:     xz-plane vector
             'Z-Axis Projection'
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
        #try:
            ## i/j/k vector is nan
            #coord = CORD2R(1, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                           #comment='')
        #except Exception:
            #log.error('The coordinate system is invalid; check your cutting plane.')
            #if stop_on_failure:
                #raise
            #return None, None

        # the plane actor defines the plane of the output results,
        # not the plane of the march direction
        # xyz1: origin
        # xyz2: xzplane
        #j = np.cross(k, i)
        unused_plane_actor = self.create_plane_actor(
            xyz1, xyz2,
            coord, i, k, dim_max,
            plane_color, plane_opacity,
        )
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
                           plane_color: Color,
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

        cid = ''
        gui: MainWindow = self.gui
        gui.tool_actions.create_coordinate_system(
            cid, dim_max, label='%s' % cid, origin=origin,
            matrix_3x3=beta, coord_type='xyz')

        j = np.cross(k, i)
        plane_actor = gui._create_plane_actor_from_points(
            xyz1, xyz2, j, k, dim_max,
            representation='surface',
            color=plane_color,  # floats
            opacity=plane_opacity, # 0=transparent, 1=solid
            actor_name='smt_plane')
        props = self.gui.geometry_properties['smt_plane']
        props.set_color(plane_color)
        props.opacity = plane_opacity
        prop = plane_actor.GetProperty()
        prop.SetColor(*plane_color)
        prop.SetOpacity(plane_opacity) # 0=transparent, 1=solid
        plane_actor.VisibilityOn()
        return plane_actor

def clear_actors(gui: MainWindow, names: list[str]) -> None:
    """clears the gui actors"""
    for name in names:
        # del gui.plane_actor
        if hasattr(gui, name):
            delattr(gui, name)
        gui.clear_actor(name)
