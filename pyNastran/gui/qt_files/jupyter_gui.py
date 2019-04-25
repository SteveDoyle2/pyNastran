import sys
import pyNastran
from pyNastran.gui.arg_handling import get_inputs
from pyNastran.gui.gui import QApplication, MainWindow
from PyQt4.QtCore import QCoreApplication


def main(argv=None):
    if argv is None:
        argv = []
        print('jupyter.amain; argv was None -> %s' % argv)

    app_created = False
    app = QCoreApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
        app_created = True

    QApplication.setOrganizationName("pyNastran")
    QApplication.setOrganizationDomain(pyNastran.__website__)
    QApplication.setApplicationName("pyNastran")
    QApplication.setApplicationVersion(pyNastran.__version__)

    inputs = get_inputs(argv)
    window = MainWindow(inputs)
    try:
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        app.exec_()
        #sys.exit(app.exec_())
