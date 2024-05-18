"""
defines:
 - EditGeometryPropertiesObject
"""
from __future__ import annotations
from copy import deepcopy
from typing import Union, Any, TYPE_CHECKING

import numpy as np

from pyNastran.gui.vtk_rendering_core import vtkActor, vtkAxesActor, vtkProperty
from pyNastran.gui.menus.edit_geometry_properties.manage_actors import EditGeometryProperties
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.gui.qt_files.base_gui import BaseGui
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui import MainWindow

# EDIT ACTOR PROPERTIES
class EditGeometryPropertiesObject(BaseGui):
    """defines EditGeometryPropertiesObject"""
    def __init__(self, gui: MainWindow):
        """creates EditGeometryPropertiesObject"""
        #self.gui = gui
        super().__init__(gui)
        self._edit_geometry_properties_window_shown = False
        self._edit_geometry_properties = None

    def set_font_size(self, font_size: int) -> None:
        """sets the font size for the edit geometry properties window"""
        if self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties.set_font_size(font_size)

    def edit_geometry_properties(self) -> None:
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
        gui = self.gui
        if not hasattr(gui, 'case_keys'):
            gui.log_error('No model has been loaded.')
            return
        if not len(gui.geometry_properties):
            gui.log_error('No secondary geometries to edit.')
            return
        #print('geometry_properties.keys() =', self.geometry_properties.keys())
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        # remove smt_plane, etc. from visible list, so the menu handles it
        data = {}
        for key, value in gui.geometry_properties.items():
            if isinstance(value, AltGeometry) and not value.visible_in_geometry_properties:
                continue
            data[key] = deepcopy(value)

        data['font_size'] = gui.settings.font_size
        if not self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties = EditGeometryProperties(data, win_parent=gui)
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
            self.on_update_geometry_properties(gui.geometry_properties)
            del self._edit_geometry_properties
            self._edit_geometry_properties_window_shown = False

    def _save_geometry_properties(self, out_data):
        for name, group in out_data.items():
            if name in ['clicked_ok', 'clicked_cancel']:
                continue
            gui = self.gui
            if name not in gui.geometry_properties:
                # we've deleted the actor
                continue

            geom_prop = gui.geometry_properties[name]
            if isinstance(geom_prop, CoordProperties):
                pass
            elif isinstance(geom_prop, AltGeometry):
                geom_prop.color = group.color
                geom_prop.line_width = group.line_width
                geom_prop.opacity = group.opacity
                geom_prop.point_size = group.point_size
            else:  # pragma: no cover
                raise NotImplementedError(geom_prop)

    def on_update_geometry_properties_override_dialog(
        self, geometry_properties: dict[str, AltGeometry]) -> None:
        """
        Update the goemetry properties and overwrite the options in the
        edit geometry properties dialog if it is open.

        Parameters
        -----------
        geometry_properties : dict {str : CoordProperties or AltGeometry}
            Dictionary from name to properties object. Only the names included
            in ``geometry_properties`` are modified.
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

    def on_update_geometry_properties_window(
        self, geometry_properties: dict[str, AltGeometry]) -> None:
        """updates the EditGeometryProperties window"""
        if self._edit_geometry_properties_window_shown:
            self._edit_geometry_properties.on_update_geometry_properties_window(
                geometry_properties)

    def on_update_geometry_properties(self, out_data, name=None,
                                      write_log: bool=True) -> None:
        """
        Applies the changed properties to the different actors if
        something changed.

        Note that some of the values are limited.  This prevents
        points/lines from being shrunk to 0 and also the actor being
        actually "hidden" at the same time.  This prevents confusion
        when you try to show the actor and it's not visible.
        """
        lines1 = []
        lines2 = []
        gui = self.gui
        geometry_properties = gui.geometry_properties
        if name is None:
            for namei, group1 in out_data.items():
                if namei in ['clicked_ok', 'clicked_cancel']:
                    continue
                group2 = geometry_properties.get(namei, None)
                update_group2 = map_group1_results_to_group2(group1, group2)
                self._update_ith_geometry_properties(namei, group1, lines1, render=False)
                if update_group2:
                    self._update_ith_geometry_properties(namei, group2, lines2, render=False)
        else:
            group1 = out_data[name]
            group2 = geometry_properties.get(name, None)
            update_group2 = map_group1_results_to_group2(group1, group2)
            self._update_ith_geometry_properties(name, group1, lines1, render=False)
            if update_group2:
                self._update_ith_geometry_properties(name, group2, lines2, render=False)

        gui.vtk_interactor.Render()
        if write_log and lines1:
            msg = 'out_data = {\n'
            msg += ''.join(lines1)
            msg += '}\n'
            msg += 'self.on_update_geometry_properties(out_data)'
            gui.log_command(msg)

    def _update_ith_geometry_properties(
        self, namei: str,
        group: Union[AltGeometry, CoordProperties],
        lines: list[str],
        render: bool=True) -> None:
        """updates a geometry"""
        gui = self.gui
        if namei not in gui.geometry_actors:
            #print('cant find %r' % namei)
            # we've deleted the actor
            return

        actor: vtkActor = gui.geometry_actors[namei]
        if isinstance(actor, vtkActor):
            alt_prop = gui.geometry_properties[namei]
            #if alt_prop.is_v
            label_actors = alt_prop.label_actors
            lines += self._update_geometry_properties_actor(namei, group, actor, label_actors)
        elif isinstance(actor, vtkAxesActor):
            changed = False
            is_visible1 = bool(actor.GetVisibility())
            is_visible2 = group.is_visible
            if is_visible1 != is_visible2:
                actor.SetVisibility(is_visible2)
                alt_prop = gui.geometry_properties[namei]
                alt_prop.is_visible = is_visible2
                actor.Modified()
                changed = True

            if changed:
                lines.append(f'    {namei!r} : CoordProperties(is_visible={is_visible2}),\n')
        else:  # pragma: no cover
            raise NotImplementedError(actor)
        if render:
            gui.vtk_interactor.Render()

    def _update_geometry_properties_actor(self, name: str,
                                          group: AltGeometry,
                                          actor: vtkActor,
                                          label_actors: list[Any]) -> list[str]:
        """
        Applies limits to the variables.  Then, checks to see if
        something in the group has changed.  If it has, updates
        the actor.

        Parameters
        ----------
        name : str
            the geometry proprety to update
        group : AltGeometry()
            a storage container for all the actor's properties
        actor : vtkActor()
            the actor where the properties will be applied
        label_actors : list[???]
            ???
        linewidth1 : int
            the active linewidth; unused???
        linewidth2 : int
            the new linewidth; unused???

        """
        lines = []
        changed = False
        #mapper = actor.GetMapper()
        prop: vtkProperty = actor.GetProperty()
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
            lines.append(
                f'    {name!r} : AltGeometry(self, {name!r}, '
                f'color=({color2[0]}, {color2[1]}, {color2[2]}), '
                f'line_width={line_width2}, opacity={opacity2}, '
                f'point_size={point_size2}, bar_scale={bar_scale2}, '
                f'representation={representation!r}, is_visible={is_visible2}),\n')
            prop.Modified()
        return lines

    def set_bar_scale(self, name: str, bar_scale: float) -> None:
        """
        Sets the bar scale

        Parameters
        ----------
        name : str
           the parameter to scale (e.g. TUBE_y, TUBE_z)
        bar_scale : float
           the scaling factor

        """
        print(f'set_bar_scale - GuiCommon2; name={name!r} bar_scale={bar_scale}')
        if bar_scale <= 0.0:
            return
        assert bar_scale > 0.0, 'bar_scale=%r' % bar_scale
        gui: MainWindow = self.gui

        # bar_y : (nbars, 6) float ndarray
        #     the xyz coordinates for (node1, node2) of the y/z axis of the bar
        #     xyz1 is the centroid
        #     xyz2 is the end point of the axis with a length_xyz with a bar_scale of 1.0
        bar_y = gui.bar_lines[name]

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
            bad_eids = gui.bar_eids[name][izero]
            gui.log.error('The following elements have zero length...%s' % bad_eids)

        # v = dxyz / length_xyz *  bar_scale
        # xyz2 = xyz1 + v

        nnodes = len(length_xyz)
        grid = gui.alt_grids[name]
        points = grid.GetPoints()
        for i in range(nnodes):
            #unused_point = points.GetPoint(2*i+1)
            #print(unused_point)
            node = xyz1[i, :] + length_xyz[i] * bar_scale * dxyz[i, :]
            #print(unused_point, node)
            points.SetPoint(2 * i + 1, *node)

        grid.Modified()
        #print('update2...')

def map_group1_results_to_group2(group1: Union[CoordProperties, AltGeometry],
                                 group2: Union[CoordProperties, AltGeometry, None]) -> bool:
    """helper to also update main data instead of just the actors"""
    update_group2 = False
    if group2 is None:
        return update_group2
    elif isinstance(group1, CoordProperties):
        #label='Global XYZ', coord_type'xyz', is_visible=False, scale=9.699980926513673
        keys = ('is_visible', )
    else:
        assert isinstance(group1, AltGeometry), group1
        keys = ('point_size', 'opacity',
                '_color', # 'color_float',
                'is_visible', 'line_width', 'representation', 'bar_scale')

    for key in keys:
        if not hasattr(group1, key):
            continue
        value1 = getattr(group1, key)
        value2 = getattr(group2, key)
        if value1 != value2:
            setattr(group2, key, value1)
            update_group2 = True
    return update_group2
