from qtpy import QtCore
from qtpy.QtWidgets import QMainWindow

from cpylog import SimpleLogger
from cpylog.html_utils import str_to_html
from pyNastran.gui.menus.application_log import ApplicationLogWidget


class LoggableGui(QMainWindow):
    def __init__(self, html_logging: bool):
        self.html_logging = html_logging
        super().__init__()

    def _logg_msg(self, log_type: str, filename: str,
                  lineno: int, msg: str) -> None:
        """
        Add message to log widget trying to choose right color for it.

        Parameters
        ----------
        log_type : str
            {DEBUG, INFO, ERROR, COMMAND, WARNING} or prepend 'GUI '
        filename : str
            the active file
        lineno : int
            line number
        msg : str
            message to be displayed
        """
        if not self.html_logging:
            # standard logger
            name = '%-8s' % (log_type + ':')
            filename_n = '%s:%s' % (filename, lineno)
            msg2 = ' %-28s %s\n' % (filename_n, msg)
            print(name, msg2)
            return

        # if 'DEBUG' in log_type and not self.settings.show_debug:
        #     return
        # elif 'INFO' in log_type and not self.settings.show_info:
        #     return
        # elif 'COMMAND' in log_type and not self.settings.show_command:
        #     return
        # elif 'WARNING' in log_type and not self.settings.show_warning:
        #     return
        # elif 'ERROR' in log_type and not self.settings.show_error:
        #     return

        if log_type in ['GUI ERROR', 'GUI COMMAND', 'GUI DEBUG', 'GUI INFO', 'GUI WARNING']:
            log_type = log_type[4:]  # drop the GUI

        html_msg = str_to_html(log_type, filename, lineno, msg)
        self._log_msg(html_msg)

    def _log_msg(self, msg: str) -> None:
        """prints an HTML log message"""
        self.log_mutex.lockForWrite()
        text_cursor = self.log_widget.textCursor()
        end = text_cursor.End
        text_cursor.movePosition(end)
        text_cursor.insertHtml(msg)
        self.log_widget.ensureCursorVisible() # new message will be visible
        self.log_mutex.unlock()

    def setup_logging(self):
        self.log = SimpleLogger(
            level='debug', encoding='utf-8',
            log_func=lambda w, x, y, z: self._logg_msg(w, x, y, z))
        if not self.html_logging:
            return

        # logging needs synchronizing, so the messages from different
        # threads would not be interleave
        self.log_mutex = QtCore.QReadWriteLock()

        self.log_dock_widget = ApplicationLogWidget(self)
        self.log_widget = self.log_dock_widget.log_widget
        #self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.log_dock_widget)
        return self.log_widget
