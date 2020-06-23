import os
import sys
#from functools import partial

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from qtpy.QtWidgets import QApplication, QPushButton

from pyNastran.gui.utils.qt.dialogs import save_file_dialog
from pyNastran.gui.utils.qt.results_window import ResultsWindow
from ..wildcards import GEOM_BDF_SAVE

from .modify_map import MODIFY_MAP, UPDATE_MAP
from .modify_menu import  ModifyMenu
from .model_sidebar import Sidebar


def build_form_from_model(model):
    form = []
    objs = []

    i = 0
    nodes_form = []
    properties_form = []
    materials_form = []
    caeros_form = []
    #splines_form = []
    import numpy as np

    nids = list(model.nodes.keys())
    spoints = list(model.spoints.keys())
    ngrids = len(nids)
    nspoints = len(spoints)
    nnodes = ngrids + nspoints
    idtype = 'int32'
    max_nid = 0
    if ngrids:
        max_nid = max(nids)
    if nspoints:
        max_nid = max(max_nid, max(spoints))
    if max_nid > 100000000:
        idtype = 'int64'

    nid_type = np.zeros((nnodes, 2), dtype=idtype)
    nid_type[:ngrids, 0] = nids
    nid_type[:ngrids, 1] = 0 # grid
    nid_type[ngrids:ngrids+nspoints, 0] = spoints
    nid_type[ngrids:ngrids+nspoints, 1] = 1 # spoint

    #for nid, nid_typei in nid_type:
        #if nid_typei == 0:
            #nodes_form.append(('GRID %i' % nid, i, []))
        #elif nid_typei == 1:
            #nodes_form.append(('SPOINT %i' % nid, i, []))
        #else:  # pragma: no cover
            #raise RuntimeError(nid_type)
        #i += 1

    elements_form = []
    mass_form = []
    rigid_elements_form = []
    for eid, elem in sorted(model.elements.items()):
        elements_form.append(('%s %i' % (elem.type, eid), i, []))
        objs.append(elem)
        i += 1
    for eid, elem in sorted(model.masses.items()):
        mass_form.append(('%s %i' % (elem.type, eid), i, []))
        objs.append(elem)
        i += 1
    for eid, elem in sorted(model.rigid_elements.items()):
        rigid_elements_form.append(('%s %i' % (elem.type, eid), i, []))
        objs.append(elem)
        i += 1

    elem_form = []
    if elements_form and mass_form and rigid_elements_form:
        elem_form = [
            ('Elastic', None, elements_form),
            ('Masses', None, mass_form),
            ('Rigid', None, rigid_elements_form),
        ]
    elif elements_form and mass_form:
        elem_form = [
            ('Elastic', None, elements_form),
            ('Masses', None, mass_form),
        ]
    elif elements_form and rigid_elements_form:
        elem_form = [
            ('Elastic', None, elements_form),
            ('Rigid', None, rigid_elements_form),
        ]
    elif mass_form and rigid_elements_form:
        elem_form = [
            ('Masses', None, mass_form),
            ('Rigid', None, rigid_elements_form),
        ]
    elif elements_form:
        elem_form = [
            ('Elastic', None, elements_form),
        ]
    elif mass_form:
        elem_form = [
            ('Masses', None, mass_form),
        ]
    elif rigid_elements_form:
        elem_form = [
            ('Rigid', None, rigid_elements_form),
        ]

    for pid, prop in sorted(model.properties.items()):
        properties_form.append(('%s %i' % (prop.type, pid), i, []))
        objs.append(prop)
        i += 1
    for mid, mat in sorted(model.materials.items()):
        materials_form.append(('%s %i' % (mat.type, mid), i, []))
        objs.append(mat)
        i += 1
    for caeroid, caero in sorted(model.caeros.items()):
        caeros_form.append(('%s %i' % (caero.type, caeroid), i, []))
        objs.append(caero)
        i += 1
    #for splineid, spline in sorted(model.splines.items()):
        #splines_form.append(('%s %i' % (spline.type, splineid), i, []))
        #i += 1

    #coords_form = []
    #for cid, coord in sorted(model.coords.items()):
        #coords_form.append(('%s %i' % (coord.type, cid), i, []))
        #i += 1

    #nodes = ('Nodes', None, nodes_form)
    elements = ('Elements', None, elem_form)
    properties = ('Properties', None, properties_form)
    materials = ('Materials', None, materials_form)
    caeros = ('CAEROs', None, caeros_form)
    #splines = ('SPLINEs', None, splines_form)
    #coords = ('Coords', None, coords_form)

    #if nodes_form:
        #form.append(nodes)
    if elem_form:
        form.append(elements)
    if properties_form:
        form.append(properties)
    if materials_form:
        form.append(materials)
    if caeros_form:
        form.append(caeros)
    #if splines_form:
        #form.append(splines)
    #if coords_form:
        #form.append(coords)
    #print(form)
    assert len(objs) > 0, objs
    return form, objs

class ModelSidebar(Sidebar):
    def __init__(self, parent, nastran_io=None):
        self.bdf_filename = None
        self.model = None
        self.objs = None
        self.parent2 = parent
        self.nastran_io = nastran_io

        right_click_actions = [
            ('Print...', self.on_print, True),
            ('Stats...', self.on_stats, True),
            ('Modify...', self.on_modify, True),
        ]
        export_button = QPushButton('Export')
        export_button.clicked.connect(self.on_export)
        setup_dict = {
            4: export_button,
        }

        form = []
        parent = self.parent2
        super(ModelSidebar, self).__init__(parent, data=form, actions=right_click_actions,
                         results_window_title='Model',
                         clear_data=False, setup_dict=setup_dict, debug=True)
        self.apply_button.setVisible(False)
        self.show()

    def on_print(self, icase):
        """prints the selected card"""
        obj = self.objs[icase]
        print(obj)
        self.model.log.info('\n' + str(obj))

    def on_stats(self, icase):
        """prints the stats for the selected card"""
        stats = self.objs[icase].get_stats()
        print(stats)
        self.model.log.info(str(stats))

    def on_modify(self, icase):
        """opens a menu to modify a card"""
        obj = self.objs[icase]

        update_function_name = None
        if obj.type in UPDATE_MAP:
            update_function_name = UPDATE_MAP[obj.type]

        try:
            variables = MODIFY_MAP[obj.type]
        except KeyError:
            self.model.log.warning(f'{obj.type} does not support Modify...')
            keys = list(MODIFY_MAP.keys())
            keys.sort()
            print('keys=', keys)
            print(obj)
            return
        self.load_menu(self.model, obj, variables, update_function_name, win_parent=None)

    def load_menu(self, model, obj, variables, update_function_name, win_parent=None):
        data = {
            'font_size': 9,
            'model' : model,
            'obj' : obj,
            'variables': variables,
            'update_function_name': update_function_name,
        }
        parent = self.parent
        self.menu = ModifyMenu(data, nastran_io=self.nastran_io, win_parent=None)
        self.menu.show()
        self.menu.exec_()

    def set_model(self, model):
        self.model = model
        self.bdf_filename = model.bdf_filename
        form, objs = build_form_from_model(self.model)
        self.objs = objs
        self.result_case_window.update_data(form)

    def on_export(self):
        """exports a modified bdf model"""
        bdf_filename_base, ext = os.path.splitext(self.bdf_filename)
        default_name = bdf_filename_base + '.modified' + ext
        default_dirname = os.path.dirname(default_name)
        fname, unused_flt = save_file_dialog(self, 'Save As...', default_dirname, GEOM_BDF_SAVE)
        if not fname:
            return
        self.model.write_bdf(fname)

    def __repr__(self):
        return '<nastran_menu.ModelSidebar(...)>'

def main():  # pragma: no cover
    app = QApplication(sys.argv)

    import pyNastran
    PKG_PATH = pyNastran.__path__[0]
    MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
    bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
    #bdf_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.bdf')

    from pyNastran.bdf.bdf import read_bdf
    model = read_bdf(bdf_filename)
    print(model.get_bdf_stats())
    unused_name = 'name'
    #res_widget.update_results(form, name)
    #--------------------------------------------

    m = ModelSidebar(app)
    m.set_model(model)
    sys.exit(app.exec_())


if __name__ == "__main__":  # pragma: no cover
    main()
