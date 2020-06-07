from qtpy import QtCore
from qtpy.QtWidgets import (
    QVBoxLayout, QWidget, QDockWidget, QComboBox, QGridLayout, QSizePolicy, QLabel, QDialog,
    QPushButton,
)
#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib.figure import Figure
#from pyNastran.bdf.cards.elements.beam_connectivity import box_faces
from pyNastran.bdf.cards.elements.beam_line_connectivity import box_lines
from pyNastran.gui.qt_files.qt_matplotlib_interface import matplotlib_backend
FigureCanvas = matplotlib_backend.FigureCanvasQTAgg
NavigationToolbar = matplotlib_backend.NavigationToolbar2QT


class MatplotlibDockWdiget(QDockWidget):
    def __init__(self, title, parent):
        QDockWidget.__init__(self, title, parent)
        self.widget = QWidget(self)
        self.figure = Figure() # Figure(figsize=(4., 5), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self.widget)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.hide()
        FigureCanvas.setSizePolicy(self.canvas,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self.canvas)
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.toolbar)
        self.toolbar.show()
        self.widget.setLayout(vbox)
        self.setWidget(self.widget)


class PropertyPlotter(QDialog):
    def __init__(self, model, parent):
        super(PropertyPlotter, self).__init__(parent)
        self.parent = parent
        self.listeners = []
        self.model = model
        self.property_ids = []
        self.setWindowTitle('Property Plotter')
        #
        # property  v-------
        # layer
        #
        property_label = QLabel('Property')
        self.property_pulldown_edit = QComboBox()
        for pid, prop in sorted(self.model.properties.items()):
            name = 'Property %i: %s' % (pid, prop.type)
            self.property_ids.append(pid)
            self.property_pulldown_edit.addItem(name)

        grid = QGridLayout()
        grid.addWidget(property_label, 0, 0)
        grid.addWidget(self.property_pulldown_edit, 0, 1)
        self.property_pulldown_edit.currentIndexChanged.connect(self.on_plot)
        #self.setLayout(grid)

        plotter_widget = QWidget(self)
        self.figure = Figure(dpi=100)
        self.subplot = self.figure.add_subplot(111)
        self.subplot.grid(True)

        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(plotter_widget)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setMinimumHeight(100)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        #self.toolbar = NavigationToolbar(self.canvas, self)
        FigureCanvas.setSizePolicy(self.canvas,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self.canvas)

        vbox2 = QVBoxLayout()
        vbox2.addWidget(self.canvas)
        plotter_widget.setLayout(vbox2)

        close = QPushButton('Close')
        vbox = QVBoxLayout()
        vbox.addLayout(grid)
        vbox.addWidget(plotter_widget)
        #vbox2.addWidget(self.toolbar)
        #vbox.addStretch()
        vbox.addWidget(close)

        close.clicked.connect(self.on_close)

        self.setLayout(vbox)
        self.on_plot()
        self.show()

    def on_close(self, event):
        self.close()

    def on_plot(self, value=0):

        pid = self.property_ids[value]
        prop = self.model.properties[pid]
        plot_property(self.subplot, prop)
        self.canvas.draw()

def plot_property(subplot, prop):
    zmean = []
    thickness = []
    mids = []
    thetas = []
    title = None
    xlabel = None
    ylabel = 'Z'

    lines = []
    xs = [0., 1.]
    if prop.type in ['PCOMP', 'PCOMPG']:
        zs = prop.get_z_locations()
        zmean = (zs[1:] + zs[:-1]) / 2.
        zmin = zs.min()
        zmax = zs.max()
        for zi in zs:
            ys = [zi, zi]
            lines.append((xs, ys))
        lines.append(([0., 0.], [zmin, zmax]))
        lines.append(([1., 1.], [zmin, zmax]))
        thickness = prop.get_thicknesses()
        mids = prop.get_material_ids()
        thetas = prop.get_thetas()
    elif prop.type == 'PSHELL':
        t = prop.Thickness()
        zs = prop.get_z_locations()
        zmean = (zs[1:] + zs[:-1]) / 2.
        zmin = zs.min()
        zmax = zs.max()
        for zi in zs:
            ys = [zi, zi]
            lines.append((xs, ys))
        thickness = [t]
        mids = [prop.Mid1()]
        thetas = [None]
        lines.append(([0., 0.], [zmin, zmax]))
        lines.append(([1., 1.], [zmin, zmax]))
    elif prop.type in ['PBARL']:
        dim = prop.dim
        title = 'Mid = %s' % prop.mid
        xlabel = 'Z'
        ylabel = 'Y'
        if prop.Type == 'BOX':
            ilines, points = box_lines(dim)
            for line in ilines:
                xs = points[line, 1]
                ys = points[line, 0]
                lines.append((xs, ys))
        else:
            xs = [0., 1., 2., 3.]
            ys = [1., 0., 1., 0.]
            lines.append((xs, ys))
    else:
        xs = [0., 1., 2., 3.]
        ys = [1., 0., 1., 0.]
        zmean = []
        thickness = []
        mids = []
        thetas = []
        lines.append((xs, ys))

    #print(lines)
    subplot.clear()
    subplot.grid(True)
    for (xs, ys) in lines:
        subplot.plot(xs, ys, 'r')

    #print('zmean = %s' % zmean)
    for zmeani, thicknessi, mid, theta in zip(zmean, thickness, mids, thetas):
        stheta = ' theta=%s' % theta if theta is not None else ''
        subplot.text(0.5, zmeani, 't=%s mid=%s%s' % (thicknessi, mid, stheta),
                     horizontalalignment='center')
    if title:
        subplot.set_title(title)

    if xlabel:
        subplot.set_xlabel(xlabel)
        subplot.set_aspect('equal', 'box')
    else:
        # the x-axis is free; 1D
        subplot.set_aspect('auto', None)

    if ylabel:
        subplot.set_ylabel(ylabel)


def main():
    from pyNastran.bdf.bdf import BDF
    model = BDF()

    pid = 1
    model.add_pshell(pid, mid1=1, t=0.1, mid2=None, twelveIt3=1.0,
                     mid3=None, tst=0.833333, nsm=0.0,
                     z1=None, z2=None, mid4=None,
                     comment='')

    pid = 2
    mids = [1, 2]
    thicknesses = [0.1, 0.2]
    thetas = None
    model.add_pcomp(pid, mids, thicknesses, thetas=thetas, souts=None, nsm=0.,
                    sb=0., ft=None, tref=0., ge=0.,
                    lam=None, z0=None, comment='')

    pid = 3
    mids = [1, 2]
    thicknesses = [0.1, 0.2]
    thetas = None
    model.add_pcomp(pid, mids, thicknesses, thetas=thetas, souts=None, nsm=0.,
                    sb=0., ft=None, tref=0., ge=0.,
                    lam='SYM', z0=0., comment='')

    pid = 5
    mid = 1
    Type = 'BOX'
    dim = [2., 1., 0.2, 0.1]
    model.add_pbarl(pid, mid, Type, dim)
    model.validate()

    pid = 6
    mid = 1
    Type = 'BAR'
    dim = [2., 1.]
    model.add_pbarl(pid, mid, Type, dim)
    model.validate()

    import sys
    from qtpy.QtWidgets import QApplication
    app = QApplication(sys.argv)
    window = PropertyPlotter(model, parent=None)

    # Enter the main loop
    app.exec_()


if __name__ == '__main__':   # pragma: no cover
    main()
