"""
defines:
 - QElementEdit
 - QNodeEdit
"""
from __future__ import print_function

from qtpy import QtCore
from qtpy.QtGui import QFocusEvent
from qtpy.QtWidgets import QLineEdit

from pyNastran.bdf.utils import write_patran_syntax_dict


class QElementEdit(QLineEdit):
    """creates a QLineEdit that can pick element ids"""
    def __init__(self, win_parent, name, parent=None, *args, **kwargs):
        super(QElementEdit, self).__init__(parent, *args, **kwargs)
        self.win_parent = win_parent
        self.name = name
        #self.focusInEvent.connect(self.on_focus)

    def focusInEvent(self, event):
        self.on_focus()
        QLineEdit.focusInEvent(self, QFocusEvent(QtCore.QEvent.FocusIn))

    def on_focus_callback(self, eids, nids, name):
        """the callback method for ``on_focus``"""
        if eids is None:
            return
        eids_str = write_patran_syntax_dict({'' : eids})
        self.setText(eids_str)

    def on_focus(self):
        """called when the QElementEdit is activated"""
        gui = self.win_parent.win_parent
        gui.mouse_actions.on_area_pick(is_eids=True, is_nids=False,
                                       name=self.name,
                                       callback=self.on_focus_callback,
                                       force=True)


class QNodeEdit(QLineEdit):
    """creates a QLineEdit that can pick node ids"""
    def __init__(self, win_parent, name, parent=None, *args, **kwargs):
        super(QNodeEdit, self).__init__(parent, *args, **kwargs)
        self.win_parent = win_parent
        self.name = name
        #self.focusInEvent.connect(self.on_focus)

    def focusInEvent(self, event):
        self.on_focus()
        QLineEdit.focusInEvent(self, QFocusEvent(QtCore.QEvent.FocusIn))

    def on_focus_callback(self, eids, nids, name):
        """the callback method for ``on_focus``"""
        if nids is None:
            return
        nids_str = write_patran_syntax_dict({'' : nids})
        self.setText(nids_str)

    def on_focus(self):
        """called when the QElementEdit is activated"""
        gui = self.win_parent.win_parent
        gui.mouse_actions.on_area_pick(is_eids=True, is_nids=True,
                                       name=self.name,
                                       callback=self.on_focus_callback,
                                       force=True)
