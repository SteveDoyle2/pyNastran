"""
defines:
 - CuttingPlaneObject

"""
from __future__ import annotations
import os
from typing import cast, TYPE_CHECKING
import numpy as np

import matplotlib
from pyNastran.gui.qt_version import qt_version
if qt_version == 'pyside2':
    matplotlib.rcParams["backend"] = 'Qt5Agg'
#else:
    #valid strings are ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo',
                       #'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf',
                       #'ps', 'svg', 'template']
    #raise NotImplementedError(qt_version)

from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.mesh_utils.cut_model_by_plane import _p1_p2_zaxis_to_cord2r
from pyNastran.bdf.mesh_utils.cutting_plane_plotter import cut_and_plot_model

from pyNastran.gui.menus.cutting_plane.cutting_plane import CuttingPlaneWindow
from pyNastran.gui.qt_files.colors import PURPLE_FLOAT
from pyNastran.gui.qt_files.base_gui import BaseGui
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.typing import ColorFloat


class CuttingPlaneObject(BaseGui):
    def __init__(self, gui: MainWindow):
        self.gui: MainWindow = gui
        self._cutting_plane_window_shown = False
        self._cutting_plane_window = None

    def set_font_size(self, font_size):
        """sets the font size for the preferences window"""
        if self._cutting_plane_window_shown:
            self._cutting_plane_window.set_font_size(font_size)

    def set_cutting_plane_menu(self):
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

        data = {
            'font_size' : settings.font_size,
            'cids' : cids,
            'plane_color' : (1., 0., 1.), # purple
            'name' : model_name,

            'model_name' : model_name,
            'clicked_ok' : False,
            'close' : False,
        }
        if not self._cutting_plane_window_shown:
            self._cutting_plane_window = CuttingPlaneWindow(data, win_parent=self.gui)
            self._cutting_plane_window.show()
            self._cutting_plane_window = True
            self._cutting_plane_window.exec_()
        else:
            self._cutting_plane_window.activateWindow()

        if data['close']:
            #if not self._cutting_plane_window._updated_preference:
            #    settings.on_set_font_size(data['font_size'])
            del self._cutting_plane_window
            gui: MainWindow = self.gui
            if hasattr(gui, 'plane_actor'):
                del gui.plane_actor
            gui.clear_actor('smt_plane')
            self._cutting_plane_window_shown = False
        else:
            self._cutting_plane_window.activateWindow()

    def make_cutting_plane_from_data(self, data):
        """Creates a cutting plane of the aero_model for the active plot result"""
        model_name = data['model_name']
        #model = self.models[model_name]

        cid_p1, p1 = data['p1']
        cid_p2, p2 = data['p2']
        cid_zaxis, zaxis = data['zaxis']
        ytol = data['ytol']
        zero_tol = data['zero_tol']
        plane_color = data['plane_color']
        method = data['method']
        plane_opacity = data['plane_opacity']
        csv_filename = None
        if 'csv_filename' in data:
            csv_filename = data['csv_filename']
        self.gui.make_cutting_plane(
            model_name,
            p1, p2, zaxis,
            method=method,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_zaxis=cid_zaxis,
            ytol=ytol, plane_atol=zero_tol,
            plane_color=plane_color,
            plane_opacity=plane_opacity, csv_filename=csv_filename)

    def make_cutting_plane(self,
                           model_name: str,
                           p1: np.ndarray,
                           p2: np.ndarray,
                           zaxis: np.ndarray,
                           method: str='Z-Axis Projection',
                           cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                           ytol=1., plane_atol=1e-5,
                           plane_color=None, plane_opacity: float=0.5,
                           csv_filename=None,
                           show=True, stop_on_failure=False):
        """Creates a cutting plane of the aero_model for the active plot result

        Plane Actor is drawn in the i-k plane
        """
        if plane_color is None:
            plane_color = PURPLE_FLOAT

        plane_color = cast(ColorFloat, plane_color)
        assert len(plane_color) == 3, plane_color
        gui: MainWindow = self.gui
        log = gui.log

        model = gui.models[model_name]
        class_name = model.__class__.__name__
        if class_name in ['BDF', 'OP2Geom']:
            out = model.get_displacement_index_xyz_cp_cd()
            unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            nids = nid_cp_cd[:, 0]
            xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
                xyz_cp, nids, icp_transform,
                cid=0)
        else:
            msg = f'{class_name!r} is not supported'
            log.error(msg)
            if stop_on_failure:
                raise RuntimeError(msg)
            return

        #xyz_min, xyz_max = model.xyz_limits
        xyz_min = xyz_cid0.min(axis=0)
        xyz_max = xyz_cid0.max(axis=0)
        assert len(xyz_min) == 3, xyz_min

        dxyz = np.abs(xyz_max - xyz_min)
        dim_max = dxyz.max()
        izero = np.where(dxyz == 0)
        dxyz[izero] = dim_max

        xyz1, xyz2, unused_z_global, i, k, origin, zaxis, xzplane = _p1_p2_zaxis_to_cord2r(
            model, p1, p2, zaxis,
            cid_p1=cid_p1, cid_p2=cid_p2, cid_zaxis=cid_zaxis,
            method=method)

        plane_actor = gui._create_plane_actor_from_points(
            xyz1, xyz2, i, k, dim_max,
            actor_name='plane')
        prop = plane_actor.GetProperty()
        prop.SetColor(*plane_color)
        prop.SetOpacity(plane_opacity) # 0=transparent, 1=solid

        gui.rend.Render()

        try:
            # i/j/k vector is nan
            coord = CORD2R(1, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                           comment='')
        except Exception:
            msg = 'The coordinate system is invalid; check your cutting plane.'
            log.error(msg)
            if stop_on_failure:
                raise RuntimeError(msg)
            return
        #print(coord)
        origin = coord.origin
        beta = coord.beta().T

        cid = ''
        gui.tool_actions.create_coordinate_system(
            cid, dim_max, label='%s' % cid, origin=origin,
            matrix_3x3=beta, coord_type='xyz')

        nnodes = xyz_cid0.shape[0]

        #case = self.result_cases[self.icase_aero]
        (obj, (i, name)) = gui.result_cases[gui.icase_fringe]
        fringe, vector = obj.get_fringe_vector_result(i, name)
        res_scalars = vector if vector is not None else fringe
        location = obj.get_location(i, name)

        if res_scalars is None:
            if plane_actor is not None:
                plane_actor.VisibilityOn()
                #print(gui.plane_actor)
                #print(dir(gui.plane_actor))
                #plane_actor.VisibilityOff()
            msg = 'No result is selected.'
            log.error(msg)
            if stop_on_failure:
                raise RuntimeError(msg)
            return

        try:
            if hasattr(obj, 'titles'):
                title = obj.titles[0]
            elif hasattr(obj, 'title'):
                title = obj.title
            else:
                raise NotImplementedError(obj)
        except Exception: # pragma: no cover
            print(f'icase={gui.icase_fringe}\n{obj}')
            raise

        plane_actor.VisibilityOn()
        if location == 'centroid':
            if hasattr(model, 'map_centroidal_result'):
                nodal_result = model.map_centroidal_result(res_scalars)
            else:
                msg = 'Centroidal results are not supported.'
                log.error(msg)
                if stop_on_failure:
                    raise RuntimeError(msg)
                return
            #np.savetxt('Cp_centroid.csv', res_scalars, header='# Cp')
            #np.savetxt('Cp_nodal.csv', nodal_result, header='# Cp')
        else:
            nodal_result = res_scalars
        assert len(nodal_result) == nnodes

        invert_yaxis = False
        if title in ['Cp']:
            invert_yaxis = True

        cut_and_plot_model(title, p1, p2, zaxis,
                           model, coord, nodal_result,
                           log,
                           ytol,
                           plane_atol=plane_atol,
                           csv_filename=csv_filename,
                           invert_yaxis=invert_yaxis,
                           cut_type='edge', show=show)

    #def create_plane_actor(self,
                           #xyz1: np.ndarray,
                           #xyz2: np.ndarray,
                           #coord: CORD2R,
                           #i: np.ndarray,
                           #k: np.ndarray,
                           #dim_max: float,
                           #plane_color: tuple[float, float, float],
                           #plane_opacity: float,
                           #) -> tuple[vtkActor, vtkProperty]:
        #"""
        #The plane is defined in the ik plane

        #Parameters
        #----------
        #xyz1 / xyz2 / xyz3 : (3,) float ndarray
            #the 1=starting 2=ending, 3=normal coordinates of the
            #coordinate frames to create in the cid=0 frame
        #i / k : (3,) float ndarray
            #the i and k vectors of the coordinate system

        #"""
        #gui: MainWindow = self.gui
        #plane_actor = gui._create_plane_actor_from_points(
            #xyz1, xyz2, i, k, dim_max,
            #actor_name='plane')
        #prop = plane_actor.GetProperty()
        #prop.SetColor(*plane_color)
        #prop.SetOpacity(plane_opacity) # 0=transparent, 1=solid
        #plane_actor.VisibilityOn()
        #gui.rend.Render()
        #return plane_actor, prop
