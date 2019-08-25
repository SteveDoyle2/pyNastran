import os

#from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication, QLabel, QPushButton, QLineEdit,
    QGridLayout, QHBoxLayout, QVBoxLayout,
    QComboBox, )
#from qtpy.compat import getexistingdirectory
from qtpy.QtWidgets import QLabel, QGridLayout, QVBoxLayout, QLineEdit

from pyNastran.gui.utils.qt.qelement_edit import QNodeEdit, QElementEdit
from pyNastran.gui.utils.qt.pydialog import PyDialog  # check_int, check_float
from pyNastran.gui.menus.results_sidebar import ResultsWindow
#from pyNastran.gui.menus.results_sidebar_utils import (
    #get_cases_from_tree, #build_pruned_tree)

field_map = {
    'CAERO1' : [
        'eid',
        'pid',
        ('lspan', 'nspan'),
        ('lchord', 'nchord'),
        ('p1', 'x12'),
        ('p4', 'x43'),
    ],
}
print(field_map)

class ModifyWindow(PyDialog):
    def __init__(self, data, win_parent=None):
        card = data['card']

class ElementWindow(PyDialog):
    """
    +------------------------+
    |                        |
    |- + controls            |
    |  |                     |
    |  +-- name              |
    |  |                     |
    |  +-- name              |
    |                        |
    | Create/Edit/Delete     |
    | Type             SKIN3 |
    | eid              ____  |
    | pid              ____  |
    | n1               ____  |
    | n2               ____  |
    | n3               ____  |
    |                        |
    |                 Apply  |
    |                        |
    +------------------------+

    TODO: handle losing focus not still picking, which can cause the property_id
          to get set to the element_id/node_id
    """
    def __init__(self, data, controls, win_parent=None):
        """create a cONTROL surface"""
        PyDialog.__init__(self, data, win_parent)
        self._updated_menu = False

        self.controls = controls

        self.comment_label = QLabel('Comment')
        self.comment_edit = QLineEdit()

        self.eid_label = QLabel('Element ID')
        self.eid_edit = QElementEdit(self, str(''), pick_style='single', tab_to_next=False)

        self.pid_label = QLabel('Property ID')
        self.pid_edit = QLineEdit()

        self.mid_label = QLabel('Material ID')
        self.mid_edit = QLineEdit()

        self.n1_label = QLabel('Node 1')
        self.n2_label = QLabel('Node 2')
        self.n3_label = QLabel('Node 3')
        self.n4_label = QLabel('Node 4')
        self.n5_label = QLabel('Node 5')
        self.n6_label = QLabel('Node 6')
        self.n7_label = QLabel('Node 7')
        self.n8_label = QLabel('Node 8')
        self.n9_label = QLabel('Node 9')
        self.n10_label = QLabel('Node 10')
        self.n1_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n2_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n3_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n4_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n5_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n6_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n7_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n8_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n9_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)
        self.n10_edit = QNodeEdit(self, str(''), pick_style='single', tab_to_next=True)

        for inode in range(3, 10+1): # 3-10
            inode_label = 'n%i_label' % inode
            inode_edit = 'n%i_edit' % inode
            getattr(self, inode_label).setVisible(False)
            getattr(self, inode_edit).setVisible(False)

        self.mcsid_label = QLabel('Material Coord')
        self.mcsid_pulldown = QComboBox()

        self.element_type_label = QLabel('Element Type')
        self.element_type_pulldown = QComboBox()

        ELEMENT_TYPES = ['CROD', 'CONROD', 'CTRIA3', 'CQUAD4']
        self.element_types = ELEMENT_TYPES
        for element_type in ELEMENT_TYPES:
            self.element_type_pulldown.addItem(element_type)
        self.element_type = ELEMENT_TYPES[0]


        self.method_type_label = QLabel('Method Type')
        self.method_type_pulldown = QComboBox()
        self.method_types = ['Create', 'Edit', 'Delete']
        METHOD_TYPES = ['Create', 'Edit', 'Delete']
        for method_type in METHOD_TYPES:
            self.method_type_pulldown.addItem(method_type)
        self.method_type = METHOD_TYPES[0]

        #cases = get_cases_from_tree(self.controls)
        #parent = self
        #name = 'main'
        #data = self.controls
        #choices = cases

        #self.results_widget_label = QLabel('Results:')
        #self.results_widget = ResultsWindow(
            #parent, name, data, choices,
            #include_clear=False, include_delete=True,
            #include_results=False)


        self.add_button = QPushButton('Create')
        self.delete_button = QPushButton('Delete')
        self.apply_button = QPushButton('Apply')
        self.setup_layout()
        self.setup_connections()

    def setup_connections(self):
        self.element_type_pulldown.currentIndexChanged.connect(self.on_element_type)
        self.method_type_pulldown.currentIndexChanged.connect(self.on_method_type)

    def on_element_type(self, value):
        """
        element type pulldown (BEAM, SKIN3, SKIN4)

        Parameters
        ----------
        value : int
            index in element_types
        """
        #animation_types = [BEAM, SKIN3, SKIN4]
        element_type = self.element_types[value]
        #func_map = {
            #'CROD': self.on_show_crod,
            #'CONROD': self.on_show_conrod,
            #'CTRIA3': self.on_show_shell,
            #'CQUAD4': self.on_show_shell,
        #}
        #func = func_map[element_type]
        param_map = {
            'CROD': ['pid'],
            'CONROD': ['mid'],
            'CTRIA3': ['pid', 'n3', 'mcsid'],
            'CQUAD4': ['pid', 'n3', 'n4', 'mcsid'],
        }
        params = param_map[element_type]
        self._update_func(element_type, params)

    def _update_func(self, element_type, params):
        """updates the menu"""
        is_pid = 'pid' in params
        is_mid = 'mid' in params
        is_mcsid = 'mcsid' in params
        self.pid_label.setVisible(is_pid)
        self.pid_edit.setVisible(is_pid)

        self.mid_label.setVisible(is_mid)
        self.mid_edit.setVisible(is_mid)

        for inode in range(3, 10+1): # 3-10
            is_ni = 'n%i' % inode in params
            inode_label = 'n%i_label' % inode
            inode_edit = 'n%i_edit' % inode
            getattr(self, inode_label).setVisible(is_ni)
            getattr(self, inode_edit).setVisible(is_ni)

        self.mcsid_label.setVisible(is_mcsid)
        self.mcsid_pulldown.setVisible(is_mcsid)

    def on_method_type(self, value):
        """
        element type pulldown (Create, Edit, Delete)

        Parameters
        ----------
        value : int
            index in element_types
        """
        #animation_types = [BEAM, SKIN3, SKIN4]
        method_type = self.method_types[value]
        enable = True
        if method_type == 'Create':
            pass
            #self.on_show_beam()
        elif method_type == 'Edit':
            pass
            #self.on_show_skin3()
        elif method_type == 'Delete':
            enable = False
            #self.on_show_skin4()
        else:
            raise NotImplementedError('value = ', value)
        self.pid_edit.setEnabled(enable)
        self.mid_edit.setEnabled(enable)
        self.n1_edit.setEnabled(enable)
        self.n2_edit.setEnabled(enable)

        for inode in range(3, 10+1): # 3-10
            inode_label = 'n%i_label' % inode
            inode_edit = 'n%i_edit' % inode
            getattr(self, inode_label).setEnabled(enable)
            getattr(self, inode_edit).setEnabled(enable)
        #self.n3_edit.setEnabled(enable)
        #self.n4_edit.setEnabled(enable)
        self.mcsid_pulldown.setEnabled(enable)
        self.comment_edit.setEnabled(enable)

    def setup_layout(self):
        irow = 0
        grid = QGridLayout()

        grid.addWidget(self.method_type_label, irow, 1)
        grid.addWidget(self.method_type_pulldown, irow, 2)
        irow += 1

        grid.addWidget(self.element_type_label, irow, 1)
        grid.addWidget(self.element_type_pulldown, irow, 2)
        irow += 1

        grid.addWidget(self.eid_label, irow, 1)
        grid.addWidget(self.eid_edit, irow, 2)
        irow += 1

        grid.addWidget(self.pid_label, irow, 1)
        grid.addWidget(self.pid_edit, irow, 2)
        irow += 1

        grid.addWidget(self.mid_label, irow, 1)
        grid.addWidget(self.mid_edit, irow, 2)
        irow += 1

        for inode in range(1, 10+1): # 1-10
            inode_label = 'n%i_label' % inode
            inode_edit = 'n%i_edit' % inode
            ni_label = getattr(self, inode_label)
            ni_edit = getattr(self, inode_edit)

            grid.addWidget(ni_label, irow, 1)
            grid.addWidget(ni_edit, irow, 2)
            irow += 1

        grid.addWidget(self.mcsid_label, irow, 1)
        grid.addWidget(self.mcsid_pulldown, irow, 2)
        irow += 1

        grid.addWidget(self.comment_label, irow, 1)
        grid.addWidget(self.comment_edit, irow, 2)
        irow += 1


        ok_cancel_box = QHBoxLayout()
        ok_cancel_box.addWidget(self.apply_button)
        ok_cancel_box.addWidget(self.add_button)
        ok_cancel_box.addWidget(self.delete_button)

        vbox = QVBoxLayout()
        #vbox.addWidget(self.results_widget_label)
        #vbox.addWidget(self.results_widget)
        vbox.addLayout(grid)
        vbox.addStretch()
        vbox.addLayout(ok_cancel_box)

        self.setLayout(vbox)
        self.setWindowTitle('Define Element')

        self.mid_label.setVisible(False)
        self.mid_edit.setVisible(False)

        self.n3_label.setVisible(False)
        self.n3_edit.setVisible(False)

        self.n4_label.setVisible(False)
        self.n4_edit.setVisible(False)

        self.mcsid_label.setVisible(False)
        self.mcsid_pulldown.setVisible(False)


    #def on_add(self):
        #self.validate()
    #def validate(self):
        #name = self.name_edit.text().strip()
        #points = self.points_edit.text().strip().split()
        #skins = self.skins_edit.text().strip().split()
        ##self.controls
        #return True

def main(): # pragma: no cover
    """test example for AnimationWindow"""
    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)


    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window

    #from pyNastran.gui.menus.legend.animation import AnimationWindow
    data2 = {
        'font_size' : 8,
        'icase_fringe' : 1,
        'icase_disp' : 2,
        'icase_vector' : 3,

        'name' : 'cat',
        'time' : 2,
        'frames/sec' : 30,
        'resolution' : 1,
        'iframe' : 0,
        'is_scale' : False,
        'dirname' : os.getcwd(),
        'scale' : 2.0,
        'default_scale' : 10,

        'arrow_scale' : 3.0,
        'default_arrow_scale' : 30,

        #'phase' : 0.,
        'phase' : None,
        'default_phase' : 120.,
        #'default_phase' : None,

        #'start_time' : 0.,
        #'end_time' : 0.,
        'default_time' : 0.,
        'icase_start' : 10,
        'icase_delta' : 3,
        'stress_min' : 0.,
        'stress_max' : 1000.,
    }
    data2['phase'] = 0.  # uncomment for phase

    form = [
        [u'Controls', None, [
            (u'Elevator', 0, []),
            (u'Flap', 1, []),
        ]],
    ]
    #[0, 1, 2, 3, 4, 5, 6, 7, 8]
    main_window = ElementWindow(data2, controls=form)
    main_window.show()
    # Enter the main loop
    app.exec_()

if __name__ == "__main__": # pragma: no cover
    main()
