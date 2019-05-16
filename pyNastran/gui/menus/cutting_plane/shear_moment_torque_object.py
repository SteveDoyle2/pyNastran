"""
defines:
 - ShearMomentTorqueObject
"""
from __future__ import print_function
import numpy as np

from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.mesh_utils.cut_model_by_plane import _p1_p2_zaxis_to_cord2r

from pyNastran.gui.menus.cutting_plane.shear_moment_torque import ShearMomentTorqueWindow
from pyNastran.gui.qt_files.colors import PURPLE_FLOAT


class ShearMomentTorqueObject(object):
    """wrapper around ShearMomentTorqueWindow"""
    def __init__(self, gui):
        self.gui = gui
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
        (obj, (unused_i, unused_name)) = gui.result_cases[icase]
        if not hasattr(obj, 'gpforce_array'):
            gui.log.error('Select a Grid Point Forces result.')
            return

        gpforce = obj.gpforce_array
        data = {
            'font_size' : settings.font_size,
            'cids' : cids,
            'plane_color' : (1., 0., 1.), # purple
            'plane_opacity' : 0.9,

            'model_name' : model_name,
            'gpforce' : gpforce,
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
            if hasattr(gui, 'plane_actor'):
                del gui.plane_actor
            gui.clear_actor('smt_plane')
            self._smt_window_shown = False
        else:
            self._smt_window.activateWindow()

    def make_smt_from_data(self, data, show=True):
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
        gpforce = data['gpforce']
        nplanes = data['nplanes']
        #model = self.models[model_name]

        cid_p1, p1 = data['p1']
        cid_p2, p2 = data['p2']
        cid_p3, p3 = data['p3']
        cid_zaxis, zaxis = data['zaxis']
        plane_color = data['plane_color']
        plane_opacity = data['plane_opacity']
        method = data['method']
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
            show=show)

    def plot_shear_moment_torque(self, model_name, gpforce,
                                 p1, p2, p3, zaxis,
                                 method='Z-Axis Projection',
                                 cid_p1=0, cid_p2=0, cid_p3=0, cid_zaxis=0,
                                 nplanes=20, plane_color=None, plane_opacity=0.5,
                                 csv_filename=None, show=True):
        """Creates a shear moment torque plot for the active plot result"""
        log = self.gui.log
        if plane_color is None:
            plane_color = PURPLE_FLOAT
        assert len(plane_color) == 3, plane_color

        model = self.gui.models[model_name]
        out = model.get_displacement_index_xyz_cp_cd()
        icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        nids = nid_cp_cd[:, 0]
        nid_cd = nid_cp_cd[:, [0, 2]]
        xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nids, icp_transform,
            cid=0)

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
        xyz3 = model.coords[cid_p3].transform_node_to_global(p3)
        dx = xyz3[0] - xyz1[0]
        stations = np.linspace(0., dx, num=nplanes, endpoint=True)
        x = stations

        try:
            # i/j/k vector is nan
            coord = CORD2R(1, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                           comment='')
        except:
            log.error('The coordinate system is invalid; check your cutting plane.')
            return None, None
        origin = coord.origin
        beta = coord.beta().T

        cid = ''
        self.gui.create_coordinate_system(
            cid, dim_max, label='%s' % cid, origin=origin,
            matrix_3x3=beta, coord_type='xyz')

        plane_actor = self.gui._create_plane_actor_from_points(
            xyz1, xyz2, i, k, dim_max,
            actor_name='smt_plane')
        props = self.gui.geometry_properties['smt_plane']
        props.set_color(plane_color)
        props.opacity = plane_opacity
        prop = plane_actor.GetProperty()
        prop.SetColor(*plane_color)
        prop.SetOpacity(plane_opacity) # 0=transparent, 1=solid
        plane_actor.VisibilityOn()

        if 0:
            point_actor = self.gui.create_point_actor_from_points(
                [xyz1, xyz3], point_size=8, actor_name='smt_points')
            prop = point_actor.GetProperty()
            prop.SetColor(*plane_color)
            prop.SetOpacity(plane_opacity) # 0=transparent, 1=solid
            point_actor.VisibilityOn()
        self.gui.rend.Render()
        self.gui.Render()

        eids, element_centroids_cid0 = get_element_centriods(model)
        force_sum, moment_sum = gpforce.shear_moment_diagram(
            xyz_cid0, eids, nids, icd_transform,
            element_centroids_cid0,
            model.coords, nid_cd, stations, coord,
            idir=0, itime=0, debug=False, logger=log)
        plot_smt(x, force_sum, moment_sum, show=show)
        return force_sum, moment_sum


def plot_smt(x, force_sum, moment_sum, show=True):
    """plots the shear, moment, torque plots"""
    import matplotlib.pyplot as plt
    plt.close()
    #f, ax = plt.subplots()
    # ax = fig.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 0], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Axial')
    ax.grid(True)

    fig = plt.figure(2)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 1], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Shear Y')
    ax.grid(True)

    fig = plt.figure(3)
    ax = fig.gca()
    ax.plot(x, force_sum[:, 2], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Shear Z')
    ax.grid(True)

    fig = plt.figure(4)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 0], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Torque')
    ax.grid(True)

    fig = plt.figure(5)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 1], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Moment Y')
    ax.grid(True)

    fig = plt.figure(6)
    ax = fig.gca()
    ax.plot(x, moment_sum[:, 2], '-*')
    ax.set_xlabel('X')
    ax.set_ylabel('Moment Z')
    ax.grid(True)

    if show:
        plt.show()


def get_element_centriods(model):
    """gets the element ids and their centroids"""
    eids = []
    element_centroids_cid0 = []
    for eid, elem in sorted(model.elements.items()):
        eids.append(eid)
        element_centroids_cid0.append(elem.Centroid())

    eids = np.array(eids, dtype='int32')
    element_centroids_cid0 = np.array(element_centroids_cid0, dtype='float64')
    return eids, element_centroids_cid0
