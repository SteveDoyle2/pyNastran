"""
defines:
 - EditGeometryPropertiesObject
"""
from copy import deepcopy
import numpy as np
import vtk
from pyNastran.gui.menus.edit_geometry_properties.manage_actors import EditGeometryProperties
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry


# EDIT ACTOR PROPERTIES
class EditGeometryPropertiesObject:
    """defines EditGeometryPropertiesObject"""
    def __init__(self, gui):
        """creates EditGeometryPropertiesObject"""
        self.gui = gui
        self._edit_geometry_properties_window_shown = False
        self._edit_geometry_properties = None

    def set_font_size(self, font_size):
        """sets the font size for the edit geometry properties window"""
        if self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties.set_font_size(font_size)

    def edit_geometry_properties(self):
        """
        Opens a dialog box to set:

        +--------+----------+
        |  Name  |  String  |
        +--------+----------+
        |  Min   |  Float   |
        +--------+----------+
        |  Max   |  Float   |
        +--------+----------+
        | Format | pyString |
        +--------+----------+
        """
        if not hasattr(self.gui, 'case_keys'):
            self.gui.log_error('No model has been loaded.')
            return
        if not len(self.gui.geometry_properties):
            self.gui.log_error('No secondary geometries to edit.')
            return
        #print('geometry_properties.keys() =', self.geometry_properties.keys())
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        data = deepcopy(self.gui.geometry_properties)
        data['font_size'] = self.gui.settings.font_size
        if not self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties = EditGeometryProperties(data, win_parent=self.gui)
            self._edit_geometry_properties.show()
            self._edit_geometry_properties_window_shown = True
            self._edit_geometry_properties.exec_()
        else:
            self._edit_geometry_properties.activateWindow()

        if 'clicked_ok' not in data:
            self._edit_geometry_properties.activateWindow()
            return

        if data['clicked_ok']:
            self.on_update_geometry_properties(data)
            self._save_geometry_properties(data)
            del self._edit_geometry_properties
            self._edit_geometry_properties_window_shown = False
        elif data['clicked_cancel']:
            self.on_update_geometry_properties(self.gui.geometry_properties)
            del self._edit_geometry_properties
            self._edit_geometry_properties_window_shown = False

    def _save_geometry_properties(self, out_data):
        for name, group in out_data.items():
            if name in ['clicked_ok', 'clicked_cancel']:
                continue

            if name not in self.gui.geometry_properties:
                # we've deleted the actor
                continue

            geom_prop = self.gui.geometry_properties[name]
            if isinstance(geom_prop, CoordProperties):
                pass
            elif isinstance(geom_prop, AltGeometry):
                geom_prop.color = group.color
                geom_prop.line_width = group.line_width
                geom_prop.opacity = group.opacity
                geom_prop.point_size = group.point_size
            else:
                raise NotImplementedError(geom_prop)

    def on_update_geometry_properties_override_dialog(self, geometry_properties):
        """
        Update the goemetry properties and overwite the options in the
        edit geometry properties dialog if it is open.

        Parameters
        -----------
        geometry_properties : dict {str : CoordProperties or AltGeometry}
            Dictionary from name to properties object. Only the names included in
            ``geometry_properties`` are modified.
        """
        if self._edit_geometry_properties_window_shown:
            # Override the output state in the edit geometry properties diaglog
            # if the button is pushed while the dialog is open. This prevent the
            # case where you close the dialog and the state reverts back to
            # before you hit the button.
            for name, prop in geometry_properties.items():
                self._edit_geometry_properties.out_data[name] = prop
                if self._edit_geometry_properties.active_key == name:
                    index = self._edit_geometry_properties.table.currentIndex()
                    self._edit_geometry_properties.update_active_key(index)
        self.on_update_geometry_properties(geometry_properties)

    def on_update_geometry_properties_window(self, geometry_properties):
        """updates the EditGeometryProperties window"""
        if self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties.on_update_geometry_properties_window(
                geometry_properties)

    def on_update_geometry_properties(self, out_data, name=None, write_log=True):
        """
        Applies the changed properties to the different actors if
        something changed.

        Note that some of the values are limited.  This prevents
        points/lines from being shrunk to 0 and also the actor being
        actually "hidden" at the same time.  This prevents confusion
        when you try to show the actor and it's not visible.
        """
        lines = []
        if name is None:
            for namei, group in out_data.items():
                if namei in ['clicked_ok', 'clicked_cancel']:
                    continue
                self._update_ith_geometry_properties(namei, group, lines, render=False)
        else:
            group = out_data[name]
            self._update_ith_geometry_properties(name, group, lines, render=False)

        self.gui.vtk_interactor.Render()
        if write_log and lines:
            msg = 'out_data = {\n'
            msg += ''.join(lines)
            msg += '}\n'
            msg += 'self.on_update_geometry_properties(out_data)'
            self.gui.log_command(msg)

    def _update_ith_geometry_properties(self, namei, group, lines, render=True):
        """updates a geometry"""
        if namei not in self.gui.geometry_actors:
            #print('cant find %r' % namei)
            # we've deleted the actor
            return

        actor = self.gui.geometry_actors[namei]
        if isinstance(actor, vtk.vtkActor):
            alt_prop = self.gui.geometry_properties[namei]
            label_actors = alt_prop.label_actors
            lines += self._update_geometry_properties_actor(namei, group, actor, label_actors)
        elif isinstance(actor, vtk.vtkAxesActor):
            changed = False
            is_visible1 = bool(actor.GetVisibility())
            is_visible2 = group.is_visible
            if is_visible1 != is_visible2:
                actor.SetVisibility(is_visible2)
                alt_prop = self.gui.geometry_properties[namei]
                alt_prop.is_visible = is_visible2
                actor.Modified()
                changed = True

            if changed:
                lines.append('    %r : CoordProperties(is_visible=%s),\n' % (
                    namei, is_visible2))
        else:
            raise NotImplementedError(actor)
        if render:
            self.gui.vtk_interactor.Render()

    def _update_geometry_properties_actor(self, name, group, actor, label_actors):
        """
        Updates an actor

        Parameters
        ----------
        name : str
            the geometry proprety to update
        group : AltGeometry()
            a storage container for all the actor's properties
        actor : vtkActor()
            the actor where the properties will be applied
        label_actors : ???
            ???
        linewidth1 : int
            the active linewidth; unused???
        linewidth2 : int
            the new linewidth; unused???
        """
        lines = []
        changed = False
        #mapper = actor.GetMapper()
        prop = actor.GetProperty()
        backface_prop = actor.GetBackfaceProperty()

        if name == 'main' and backface_prop is None:
            # don't edit these
            # we're lying about the colors to make sure the
            # colors aren't reset for the Normals
            color1 = prop.GetDiffuseColor()
            color2 = color1
            assert color1[1] <= 1.0, color1
        else:
            color1 = prop.GetDiffuseColor()
            assert color1[1] <= 1.0, color1
            color2 = group.color_float
            #print('line2646 - name=%s color1=%s color2=%s' % (name, str(color1), str(color2)))
            #color2 = group.color

        opacity1 = prop.GetOpacity()
        opacity2 = group.opacity
        opacity2 = max(0.1, opacity2)

        line_width1 = prop.GetLineWidth()
        line_width2 = group.line_width
        line_width2 = max(1, line_width2)

        point_size1 = prop.GetPointSize()
        point_size2 = group.point_size
        point_size2 = max(1, point_size2)

        representation = group.representation
        alt_prop = self.gui.geometry_properties[name]
        #representation = alt_prop.representation
        #is_visible1 = alt_prop.is_visible
        is_visible1 = bool(actor.GetVisibility())
        is_visible2 = group.is_visible
        #print('is_visible1=%s is_visible2=%s'  % (is_visible1, is_visible2))

        bar_scale1 = alt_prop.bar_scale
        bar_scale2 = group.bar_scale
        # bar_scale2 = max(0.0, bar_scale2)

        #print('name=%s color1=%s color2=%s' % (name, str(color1), str(color2)))
        if color1 != color2:
            #print('color_2662[%s] = %s' % (name, str(color1)))
            assert isinstance(color1[0], float), color1
            prop.SetDiffuseColor(color2)
            changed = True
        if line_width1 != line_width2:
            line_width2 = max(1, line_width2)
            prop.SetLineWidth(line_width2)
            changed = True
        if opacity1 != opacity2:
            #if backface_prop is not None:
                #backface_prop.SetOpacity(opacity2)
            prop.SetOpacity(opacity2)
            changed = True
        if point_size1 != point_size2:
            prop.SetPointSize(point_size2)
            changed = True
        if representation == 'bar' and bar_scale1 != bar_scale2:
            #print('name=%s rep=%r bar_scale1=%s bar_scale2=%s' % (
                #name, representation, bar_scale1, bar_scale2))
            self.set_bar_scale(name, bar_scale2)

        if is_visible1 != is_visible2:
            actor.SetVisibility(is_visible2)
            alt_prop.is_visible = is_visible2
            #prop.SetViPointSize(is_visible2)
            actor.Modified()
            for label_actor in label_actors:
                label_actor.SetVisibility(is_visible2)
                label_actor.Modified()
            changed = True

        if changed:
            lines.append('    %r : AltGeometry(self, %r, color=(%s, %s, %s), '
                         'line_width=%s, opacity=%s, point_size=%s, bar_scale=%s, '
                         'representation=%r, is_visible=%s),\n' % (
                             name, name, color2[0], color2[1], color2[2], line_width2,
                             opacity2, point_size2, bar_scale2, representation, is_visible2))
            prop.Modified()
        return lines


    def set_bar_scale(self, name, bar_scale):
        """
        Sets the bar scale

        Parameters
        ----------
        name : str
           the parameter to scale (e.g. TUBE_y, TUBE_z)
        bar_scale : float
           the scaling factor
        """
        #print('set_bar_scale - GuiCommon2; name=%s bar_scale=%s' % (name, bar_scale))
        if bar_scale <= 0.0:
            return
        assert bar_scale > 0.0, 'bar_scale=%r' % bar_scale

        # bar_y : (nbars, 6) float ndarray
        #     the xyz coordinates for (node1, node2) of the y/z axis of the bar
        #     xyz1 is the centroid
        #     xyz2 is the end point of the axis with a length_xyz with a bar_scale of 1.0
        bar_y = self.gui.bar_lines[name]

        #dy = c - yaxis
        #dz = c - zaxis
        #print('bary:\n%s' % bar_y)
        xyz1 = bar_y[:, :3]
        xyz2 = bar_y[:, 3:]
        dxyz = xyz2 - xyz1

        # vectorized version of L = sqrt(dx^2 + dy^2 + dz^2)
        length_xyz = np.linalg.norm(dxyz, axis=1)
        izero = np.where(length_xyz == 0.0)[0]
        if len(izero):
            bad_eids = self.gui.bar_eids[name][izero]
            self.gui.log.error('The following elements have zero length...%s' % bad_eids)

        # v = dxyz / length_xyz *  bar_scale
        # xyz2 = xyz1 + v

        nnodes = len(length_xyz)
        grid = self.gui.alt_grids[name]
        points = grid.GetPoints()
        for i in range(nnodes):
            #unused_point = points.GetPoint(2*i+1)
            #print(unused_point)
            node = xyz1[i, :] + length_xyz[i] * bar_scale * dxyz[i, :]
            #print(unused_point, node)
            points.SetPoint(2 * i + 1, *node)

        grid.Modified()
        #print('update2...')
