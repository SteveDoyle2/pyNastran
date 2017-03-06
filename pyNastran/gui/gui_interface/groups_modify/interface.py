from __future__ import print_function
from six import integer_types
from pyNastran.gui.gui_interface.groups_modify.groups_modify import GroupsModify, Group

def on_modify_group(self):
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
    if not len(self.groups):  # no 'main' group
        self.log_error('No main group to create.')
        return
    print('groups.keys() =', self.groups.keys())

    data = {0 : self.groups['main']}

    i = 1
    for name, group in sorted(iteritems(self.groups)):
        if name == 'main':
            continue
        data[i] = group
        i += 1
    #data = deepcopy(self.groups)

    if not self._modify_groups_window_shown:
        self._modify_groups = GroupsModify(
            data, win_parent=self, group_active=self.group_active)
        self._modify_groups.show()
        self._modify_groups_window_shown = True
        self._modify_groups.exec_()
    else:
        self._modify_groups.activateWindow()

    if 'clicked_ok' not in data:
        self._modify_groups.activateWindow()
        return

    if data['clicked_ok']:
        self.on_update_modify_group(data)
        imain = self._modify_groups.imain
        name = self._modify_groups.keys[imain]
        self.post_group_by_name(name)
        #name =
        #self._save_geometry_properties(data)
        del self._modify_groups
        self._modify_groups_window_shown = False
    elif data['clicked_cancel']:
        self.on_update_modify_group(data)
        del self._modify_groups
        self._modify_groups_window_shown = False
