from pyNastran.gui.menus.groups_modify.groups_modify import GroupsModify

def on_set_modify_groups(self):
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
    print('groups.keys() = %s' % self.groups.keys())

    group_active = self.group_active
    assert isinstance(group_active, str), group_active
    data = {
        'font_size' : self.settings.font_size,
        0 : self.groups['main'],
        'clicked_ok' : False,
        'close' : False,
    }

    i = 1
    for name, group in sorted(self.groups.items()):
        if name == 'main':
            continue
        data[i] = group
        i += 1

    if not self._modify_groups_window_shown:
        self._modify_groups_window = GroupsModify(
            data, win_parent=self, group_active=group_active)
        self._modify_groups_window.show()
        self._modify_groups_window_shown = True
        self._modify_groups_window.exec_()
    else:
        self._modify_groups_window.activateWindow()

    if data['close']:
        if not self._modify_groups_window._updated_groups:
            self._apply_modify_groups(data)
        self._modify_groups_window_shown = False
        del self._modify_groups_window
    else:
        self._modify_groups_window.activateWindow()
