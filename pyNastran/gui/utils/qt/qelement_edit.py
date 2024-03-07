"""
defines:
 - QElementLineEdit
 - QNodeLineEdit
 - QNodeElementLineEdit

 - QElementTextEdit
 - QNodeTextEdit
 - QNodeElementTextEdit

"""
import numpy as np
from qtpy import QtCore
from qtpy.QtGui import QFocusEvent
from qtpy.QtWidgets import QLineEdit, QTextEdit, QApplication

from pyNastran.bdf.utils import write_patran_syntax_dict


class QNodeElementLineEdit(QLineEdit):
    """creates a QLineEdit that can pick node/element ids"""
    def __init__(self, win_parent, name: str, parent=None,
                 is_eids: bool=True, is_nids: bool=True,
                 representation: str='wire', pick_style: str='area',
                 max_length: int=32767, cleanup: bool=True, *args, **kwargs):
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
        super(QNodeElementLineEdit, self).__init__(parent=parent, *args, **kwargs)
        self.win_parent = win_parent
        self.name = name
        self.is_eids = is_eids
        self.is_nids = is_nids
        self.representation = representation
        self.pick_style = pick_style
        self.cleanup = cleanup
        assert self.representation in {'wire', 'points', 'surface', 'points+wire', 'points+surface'}, f'representation={representation!r}'
        assert self.pick_style in {'area', 'single'}, f'pick_style={pick_style!r}'
        if max_length != 32767:
            self.setMaxLength(max_length) # default is 32767

        self.style = None
        self._is_updated = False
        self.clipboard_data = ''
        #self.focusInEvent.connect(self.on_focus)

    def focusInEvent(self, event):
        self.on_focus()
        QLineEdit.focusInEvent(self, QFocusEvent(QtCore.QEvent.FocusIn))

    #def inputRejected(self, event):
        #clipboard = QApplication.clipboard()
        #self.clipboard_data = clipboard
        #x = 1

    def on_focus_callback(self, eids: np.ndarray, nids: np.ndarray, name: str):
        """the callback method for ``on_focus``"""
        raise NotImplementedError('write_patran_syntax_dict callback')
        #eids_str = write_patran_syntax_dict({'' : eids})
        #self.setText(eids_str)

    def on_focus(self):
        """called when the QNodeElementEdit is activated"""
        win_parent = self.win_parent
        if win_parent is None:
            return
        if not hasattr(win_parent, 'win_parent'):
            print('expected gui')
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
        else:  # pragma: no cover
            raise NotImplementedError(f'pick_style={self.pick_style!r} and must be [area, single]')
        del style

        # if style has updated the gui
        if self._is_updated and 0:
            if self.style is not None:
                self.style.remove_actors()
                gui.vtk_interactor.Render()

        #self.style = style


class QElementLineEdit(QNodeElementLineEdit):
    """creates a QLineEdit that can pick element ids"""
    def __init__(self, win_parent, name: str, parent=None,
                 pick_style: str='area', tab_to_next: bool=False,
                 cleanup: bool=True, max_length: int=32767, *args, **kwargs):
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
        super(QElementLineEdit, self).__init__(
            win_parent, name, parent=parent,
            is_eids=True, is_nids=False,
            representation='wire', pick_style=pick_style,
            max_length=max_length, cleanup=cleanup, *args, **kwargs)

    def on_focus_callback(self, eids: np.ndarray, nids: np.ndarray,
                          name: str) -> None:
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


class QNodeLineEdit(QNodeElementLineEdit):
    """creates a QLineEdit that can pick node ids"""
    def __init__(self, win_parent, name: str, parent=None,
                 pick_style: str='area', tab_to_next: bool=False,
                 cleanup: bool=True, max_length: int=32767, *args, **kwargs):
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
        super(QNodeLineEdit, self).__init__(
            win_parent, name, parent=parent,
            is_eids=False, is_nids=True,
            representation='points', pick_style=pick_style,
            cleanup=cleanup, max_length=max_length, *args, **kwargs)

    def on_focus_callback(self, eids: np.ndarray, nids: np.ndarray,
                          name: str) -> None:
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



class QNodeElementTextEdit(QTextEdit):
    """creates a QTextEdit that can pick node/element ids"""
    def __init__(self, win_parent, name: str, parent=None,
                 is_eids: bool=True, is_nids: bool=True,
                 representation: str='wire', pick_style: str='area',
                 max_length: int=32767, cleanup: bool=True, *args, **kwargs):
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
        super(QNodeElementTextEdit, self).__init__(parent=parent, *args, **kwargs)
        self.win_parent = win_parent
        self.name = name
        self.is_eids = is_eids
        self.is_nids = is_nids
        self.representation = representation
        self.pick_style = pick_style
        self.cleanup = cleanup
        assert self.representation in {'wire', 'points', 'surface', 'points+wire', 'points+surface'}, f'representation={representation!r}'
        assert self.pick_style in {'area', 'single'}, f'pick_style={pick_style!r}'
        #if max_length != 32767:
            #self.setMaxLength(max_length) # default is 32767

        self.style = None
        self._is_updated = False
        #self.clipboard_data = ''
        #self.focusInEvent.connect(self.on_focus)
        if 0:  # pragma: no cover
            from qtpy.QtGui import QValidator, QDoubleValidator, QRegExpValidator
            from qtpy.QtCore import QRegularExpression, QRegExp

            # any number of:
            #  - 0-9
            #  - spaces
            #  - colon
            regex = QRegExp('([0-9]* :)*')

            #m_Validator = QRegExpValidator(regex, self.toPlainText()) # doesn't work?
            m_Validator = QRegExpValidator(regex, self)
            #"^(?=.{1,12}$)^[a-zA-Z][a-zA-Z0-9]*$"

            #m_Validator.setRegularExpression(regex)
            self.setValidator(m_Validator)
            self.connect


    def focusInEvent(self, event):
        self.on_focus()
        QTextEdit.focusInEvent(self, QFocusEvent(QtCore.QEvent.FocusIn))

    #def inputRejected(self, event):
        #clipboard = QApplication.clipboard()
        #self.clipboard_data = clipboard
        #x = 1

    def on_focus_callback(self, eids: np.ndarray, nids: np.ndarray, name: str):
        """the callback method for ``on_focus``"""
        raise NotImplementedError('write_patran_syntax_dict callback')
        #eids_str = write_patran_syntax_dict({'' : eids})
        #self.setText(eids_str)

    def on_focus(self):
        """called when the QNodeElementEdit is activated"""
        win_parent = self.win_parent
        if win_parent is None:
            return
        if not hasattr(win_parent, 'win_parent'):
            print('expected gui')
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
        else:  # pragma: no cover
            raise NotImplementedError(f'pick_style={self.pick_style!r} and must be [area, single]')
        del style

        # if style has updated the gui
        if self._is_updated and 0:  # pragma: no cover
            if self.style is not None:
                self.style.remove_actors()
                gui.vtk_interactor.Render()

        #self.style = style

    #def textChanged(self, event):
        #pass

class QElementTextEdit(QNodeElementTextEdit):
    """creates a QTextEdit that can pick element ids"""
    def __init__(self, win_parent, name: str, parent=None,
                 pick_style: str='area', tab_to_next: bool=False,
                 cleanup: bool=True, *args, **kwargs):
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
        super(QElementTextEdit, self).__init__(
            win_parent, name, parent=parent,
            is_eids=True, is_nids=False,
            representation='wire', pick_style=pick_style,
            cleanup=cleanup, *args, **kwargs)

    def on_focus_callback(self, eids: np.ndarray, nids: np.ndarray,
                          name: str) -> None:
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

class QNodeTextEdit(QNodeElementTextEdit):
    """creates a QLineEdit that can pick node ids"""
    def __init__(self, win_parent, name: str, parent=None,
                 pick_style: str='area', tab_to_next: bool=False,
                 cleanup: bool=True, max_length: int=32767, *args, **kwargs):
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
        super(QNodeTextEdit, self).__init__(
            win_parent, name, parent=parent,
            is_eids=False, is_nids=True,
            representation='points', pick_style=pick_style,
            cleanup=cleanup, max_length=max_length, *args, **kwargs)

    def on_focus_callback(self, eids: np.ndarray, nids: np.ndarray,
                          name: str) -> None:
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


def main():  # pragma: no cover
    from qtpy.QtWidgets import QVBoxLayout, QApplication, QDialog, QPushButton
    class TestWindow(QDialog):
        def __init__(self):
            super(QDialog, self).__init__()

            name = 'test'
            self.element_edit = QElementTextEdit(
                self, name, parent=None, pick_style='area',
                tab_to_next=False, cleanup=True, max_length=5)
            apply_button = QPushButton('Apply')

            vbox = QVBoxLayout()
            vbox.addWidget(self.element_edit)
            vbox.addWidget(apply_button)
            self.setLayout(vbox)

            apply_button.clicked.connect(self.on_apply)

        def on_apply(self):
            pass

    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    #1:1000000
    import sys
    # Someone is launching this directly
    # Create the QApplication
    app = QApplication(sys.argv)
    #The Main window

    main_window = TestWindow()
    main_window.show()
    app.exec_()

if __name__ == '__main__':  # pragma: no cover
    main()
