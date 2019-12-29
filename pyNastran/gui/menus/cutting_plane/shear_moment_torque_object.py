"""
defines:
 - ShearMomentTorqueObject

"""
#from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    get_element_centroids)

from pyNastran.gui.menus.cutting_plane.shear_moment_torque import ShearMomentTorqueWindow
from pyNastran.gui.qt_files.colors import PURPLE_FLOAT
from pyNastran.op2.tables.ogf_gridPointForces.smt import (
    get_nid_cd_xyz_cid0, plot_smt, setup_coord_from_plane)

class ShearMomentTorqueObject:
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
                                 method: str='Z-Axis Projection',
                                 cid_p1: int=0, cid_p2: int=0, cid_p3: int=0, cid_zaxis: int=0,
                                 nplanes: int=20, plane_color=None, plane_opacity=0.5,
                                 csv_filename=None, show=True, stop_on_failure=False):
        """
        Creates a shear moment torque plot for the active plot result

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
           'CORD2R' : typical CORD2R
            'Z-Axis Projection' : project p2 on the z-axis
        """
        log = self.gui.log
        if plane_color is None:
            plane_color = PURPLE_FLOAT
        assert len(plane_color) == 3, plane_color

        model = self.gui.models[model_name]
        nids, nid_cd, icd_transform, xyz_cid0 = get_nid_cd_xyz_cid0(model)
        #xyz1, xyz2, xyz3, i, k, origin, xzplane, dim_max, stations
        try:
            xyz1, xyz2, xyz3, i, k, coord_out, dim_max, stations = setup_coord_from_plane(
                model, xyz_cid0,
                p1, p2, p3, zaxis,
                method=method,
                cid_p1=cid_p1, cid_p2=cid_p2, cid_p3=cid_p3, cid_zaxis=cid_zaxis,
                nplanes=nplanes,
            )
        except ValueError as verror:
            model.log.error(str(verror))
            if stop_on_failure:
                raise
            return
        coord = coord_out
        #try:
            ## i/j/k vector is nan
            #coord = CORD2R(1, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                           #comment='')
        #except:
            #log.error('The coordinate system is invalid; check your cutting plane.')
            #if stop_on_failure:
                #raise
            #return None, None
        unused_plane_actor, unused_prop = self.create_plane_actor(
            xyz1, xyz2,
            coord, i, k, dim_max,
            plane_color, plane_opacity,
        )
        #if 0:
            #point_actor = self.gui.create_point_actor_from_points(
                #[xyz1, xyz3], point_size=8, actor_name='smt_points')
            #prop = point_actor.GetProperty()
            #prop.SetColor(*plane_color)
            #prop.SetOpacity(plane_opacity) # 0=transparent, 1=solid
            #point_actor.VisibilityOn()
        self.gui.rend.Render()
        self.gui.Render()

        eids, element_centroids_cid0 = get_element_centroids(model)
        force_sum, moment_sum = gpforce.shear_moment_diagram(
            xyz_cid0, eids, nids, icd_transform,
            element_centroids_cid0,
            model.coords, nid_cd, stations, coord,
            idir=0, itime=0, debug=False, log=log)
        plot_smt(stations, force_sum, moment_sum, show=show)
        return force_sum, moment_sum

    def create_plane_actor(self, xyz1, xyz2, coord, i, k, dim_max: float,
                           plane_color, plane_opacity: float):
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
        return plane_actor, prop
