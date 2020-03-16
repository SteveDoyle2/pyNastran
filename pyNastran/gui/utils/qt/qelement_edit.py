"""
defines:
 - QElementEdit
 - QNodeEdit
 - QNodeElementEdit

"""
from qtpy import QtCore
from qtpy.QtGui import QFocusEvent
from qtpy.QtWidgets import QLineEdit

from pyNastran.bdf.utils import write_patran_syntax_dict


class QNodeElementEdit(QLineEdit):
    """creates a QLineEdit that can pick node/element ids"""
    def __init__(self, win_parent, name, parent=None, is_eids=True, is_nids=True,
                 representation='wire', pick_style='area', cleanup=True, *args, **kwargs):
        """
        A node/element picker

        Parameters
        ----------
        win_parent : QDialog
            set this to self
        name : str
            the name of the model
        parent : the gui
            this is the gui, it could be win_parent.parent, win_parent.parent.gui, etc.
        is_eids : bool; default=True
            should elements be picked
        is_nids : bool; default=True
            should nodes picked
        representation : str; default='wire'
            the way that the highlighted actor should be displayed
            {wire, points, surface}
        pick_style : str; default='area'
            the style of the picker
            {area, single}

        """
        super(QNodeElementEdit, self).__init__(parent=parent, *args, **kwargs)
        self.win_parent = win_parent
        self.name = name
        self.is_eids = is_eids
        self.is_nids = is_nids
        self.representation = representation
        self.pick_style = pick_style
        self.cleanup = cleanup
        assert self.representation in ['wire', 'points', 'surface', 'points+wire', 'points+surface'], 'representation=%r' % representation
        assert self.pick_style in ['area', 'single'], 'pick_style=%r' % pick_style
        self.style = None
        self._is_updated = False
        #self.focusInEvent.connect(self.on_focus)

    def focusInEvent(self, event):
        self.on_focus()
        QLineEdit.focusInEvent(self, QFocusEvent(QtCore.QEvent.FocusIn))

    def on_focus_callback(self, eids, nids, name):
        """the callback method for ``on_focus``"""
        raise NotImplementedError('write_patran_syntax_dict callback')
        #eids_str = write_patran_syntax_dict({'' : eids})
        #self.setText(eids_str)

    def on_focus(self):
        """called when the QNodeElementEdit is activated"""
        win_parent = self.win_parent
        if win_parent is None:
            return
        gui = win_parent.win_parent
        if gui is None:
            return

        self._is_updated = False
        if self.pick_style == 'area':
            style = gui.mouse_actions.on_area_pick(
                is_eids=self.is_eids, is_nids=self.is_nids,
                representation=self.representation,
                name=self.name,
                callback=self.on_focus_callback,
                cleanup=self.cleanup,
                force=True)
        elif self.pick_style == 'single':
            style = gui.mouse_actions.on_highlight(
                is_eids=self.is_eids, is_nids=self.is_nids,
                representation=self.representation,
                name=self.name,
                callback=self.on_focus_callback,
                cleanup=self.cleanup,
                force=True)
        else:
            raise NotImplementedError('pick_style=%r and must be [area, single]' % self.pick_style)

        # if style has updated the gui
        if self._is_updated and 0:
            if self.style is not None:
                self.style.remove_actors()
                gui.vtk_interactor.Render()

        #self.style = style

class QElementEdit(QNodeElementEdit):
    """creates a QLineEdit that can pick element ids"""
    def __init__(self, win_parent, name, parent=None, pick_style='area', tab_to_next=False,
                 cleanup=True, *args, **kwargs):
        """
        An element picker

        Parameters
        ----------
        win_parent : QDialog
            set this to self
        name : str
            the name of the model
        parent : the gui
            this is the gui, it could be win_parent.parent, win_parent.parent.gui, etc.
        #representation : str; default='wire'
            #the way that the highlighted actor should be displayed
            #{wire, points, surface}
        pick_style : str; default='area'
            the style of the picker
            {area, single}
        tab_to_next : bool; default=False
            ???
        """
        self.tab_to_next = tab_to_next
        super(QElementEdit, self).__init__(win_parent, name, parent=parent,
                                           is_eids=True, is_nids=False,
                                           representation='wire', pick_style=pick_style,
                                           cleanup=cleanup, *args, **kwargs)

    def on_focus_callback(self, eids, nids, name):
        """the callback method for ``on_focus``"""
        if eids is None:
            return
        if len(eids) == 1:
            self.setText(str(eids[0]))
            return
        eids_str = write_patran_syntax_dict({'' : eids})
        self.setText(eids_str)

        self._is_updated = True
        #if self.tab_to_next:
            #self.keyPressEvent(QtCore.Qt.Key_Tab)


class QNodeEdit(QNodeElementEdit):
    """creates a QLineEdit that can pick node ids"""
    def __init__(self, win_parent, name, parent=None, pick_style='area', tab_to_next=False,
                 cleanup=True, *args, **kwargs):
        """
        A node picker

        Parameters
        ----------
        win_parent : QDialog
            set this to self
        name : str
            the name of the model
        parent : the gui
            this is the gui, it could be win_parent.parent, win_parent.parent.gui, etc.
        pick_style : str; default='area'
            the style of the picker
            {area, single}
        tab_to_next : bool; default=False
            ???
        """
        self.tab_to_next = tab_to_next
        super(QNodeEdit, self).__init__(win_parent, name, parent=parent,
                                        is_eids=False, is_nids=True,
                                        representation='points', pick_style=pick_style,
                                        cleanup=cleanup, *args, **kwargs)

    def on_focus_callback(self, eids, nids, name):
        """the callback method for ``on_focus``"""
        if nids is None:
            return
        if len(nids) == 1:
            self.setText(str(nids[0]))
            return
        nids_str = write_patran_syntax_dict({'' : nids})
        self.setText(nids_str)

        self._is_updated = True
        #if self.tab_to_next:
            #self.keyPressEvent(QtCore.Qt.Key_Tab)
