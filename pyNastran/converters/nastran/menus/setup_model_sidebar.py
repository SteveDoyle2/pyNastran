import sys

# kills the program when you hit Cntl+C from the command line
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from qtpy.QtWidgets import QApplication, QPushButton

from pyNastran.gui.utils.qt.dialogs import save_file_dialog
from pyNastran.gui.utils.qt.results_window import ResultsWindow
from pyNastran.converters.nastran.menus.model_sidebar import Sidebar
from pyNastran.converters.nastran.menus.modify_menu import  ModifyMenu
from pyNastran.converters.nastran.wildcards import GEOM_BDF_SAVE

class Var:
    def __init__(self, name, var, vartype='lineedit', pulldown_objs=None, enabled=True):
        self.name = name
        self.var = var
        self.vartype = vartype
        self.pulldown_objs = pulldown_objs
        self.enabled = enabled
        assert vartype in ['lineedit', 'pulldown', 'spinner'], vartype

MODIFY_MAP = {
    'CAERO1': [
        Var('Element ID', 'eid', enabled=False),
        Var('Property ID', 'pid', vartype='pulldown', pulldown_objs='paeros'),
        Var('iGroup', 'igroup'),
        Var('nSpan Boxes', 'nspan', vartype='spinner'),
        Var('nChord Boxes', 'nchord', vartype='spinner'),
        Var('AEFACT Span', 'lspan', vartype='pulldown', pulldown_objs='aefacts'),
        Var('AEFACT Chord', 'lchord', vartype='pulldown', pulldown_objs='aefacts'),
        Var('Point 1', 'p1'),
        Var('Distance 12', 'x12'),
        Var('Point 4', 'p4'),
        Var('Distance 43', 'x43'),
    ],
}

def build_form_from_model(model):
    form = []
    objs = []

    i = 0
    #nodes_form = []
    #properties_form = []
    #materials_form = []
    caeros_form = []
    #splines_form = []
    import numpy as np

    nids = list(model.nodes.keys())
    spoints = list(model.spoints.keys())
    ngrids = len(nids)
    nspoints = len(spoints)
    nnodes = ngrids + nspoints
    nid_type = np.zeros((nnodes, 2), dtype='int32')
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

    #elements_form = []
    #mass_form = []
    #rigid_elements_form = []
    #for eid, elem in sorted(model.elements.items()):
        #elements_form.append(('%s %i' % (elem.type, eid), i, []))
        #i += 1
    #for eid, elem in sorted(model.masses.items()):
        #mass_form.append(('%s %i' % (elem.type, eid), i, []))
        #i += 1
    #for eid, elem in sorted(model.rigid_elements.items()):
        #rigid_elements_form.append(('%s %i' % (elem.type, eid), i, []))
        #i += 1

    #elem_form = []
    #if elements_form and mass_form and rigid_elements_form:
        #elem_form = [
            #('Elastic', None, elements_form),
            #('Masses', None, mass_form),
            #('Rigid', None, rigid_elements_form),
        #]
    #elif elements_form and mass_form:
        #elem_form = [
            #('Elastic', None, elements_form),
            #('Masses', None, mass_form),
        #]
    #elif elements_form and rigid_elements_form:
        #elem_form = [
            #('Elastic', None, elements_form),
            #('Rigid', None, rigid_elements_form),
        #]
    #elif mass_form and rigid_elements_form:
        #elem_form = [
            #('Masses', None, mass_form),
            #('Rigid', None, rigid_elements_form),
        #]
    #elif elements_form:
        #elem_form = [
            #('Elastic', None, elements_form),
        #]
    #elif mass_form:
        #elem_form = [
            #('Masses', None, mass_form),
        #]
    #elif rigid_elements_form:
        #elem_form = [
            #('Rigid', None, rigid_elements_form),
        #]

    #for pid, prop in sorted(model.properties.items()):
        #properties_form.append(('%s %i' % (prop.type, pid), i, []))
        #i += 1
    #for mid, mat in sorted(model.materials.items()):
        #materials_form.append(('%s %i' % (mat.type, mid), i, []))
        #i += 1
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
    #elements = ('Elements', None, elem_form)
    #properties = ('Properties', None, properties_form)
    #materials = ('Materials', None, materials_form)
    caeros = ('CAEROs', None, caeros_form)
    #splines = ('SPLINEs', None, splines_form)
    #coords = ('Coords', None, coords_form)

    #if nodes_form:
        #form.append(nodes)
    #if elem_form:
        #form.append(elements)
    #if properties_form:
        #form.append(properties)
    #if materials_form:
        #form.append(materials)
    if caeros_form:
        form.append(caeros)
    #if splines_form:
        #form.append(splines)
    #if coords_form:
        #form.append(coords)
    #print(form)
    return form, objs

def main():  # pragma: no cover
    app = QApplication(sys.argv)

    import os
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

    class ModelSidebar(Sidebar):
        def __init__(self, model):
            form, objs = build_form_from_model(model)
            def on_print(icase):
                """prints the selected card"""
                print(objs[icase])
            def on_stats(icase):
                """prints the stats for the selected card"""
                print(objs[icase].get_stats())

            def on_modify(icase):
                """opens a menu to modify a card"""
                obj = objs[icase]
                try:
                    variables = MODIFY_MAP[obj.type]
                except KeyError:
                    print(obj)
                    return
                load_menu(model, obj, variables, win_parent=None)

            def load_menu(model, obj, variables, win_parent=None):
                data = {
                    'font_size': 9,
                    'model' : model,
                    'obj' : obj,
                    'variables': variables,
                }
                menu = ModifyMenu(data, win_parent=None)
                menu.show()
                menu.exec_()

            right_click_actions = [
                ('Print...', on_print, True),
                ('Stats...', on_stats, True),
                ('Modify...', on_modify, True),
            ]
            export_button = QPushButton('Export')
            export_button.clicked.connect(self.on_export)
            setup_dict = {
                4: export_button,
            }
            super().__init__(self, data=form, actions=right_click_actions,
                             results_window_title='Model',
                             clear_data=False, setup_dict=setup_dict, debug=True)
            self.apply_button.setVisible(False)
            self.show()

        def on_export(self):
            """exports a modified bdf model"""
            bdf_filename_base, ext = os.path.splitext(bdf_filename)
            default_name = bdf_filename_base + '.modified' + ext
            default_dirname = os.path.dirname(default_name)
            fname, flt = save_file_dialog(self, 'Save As...', default_dirname, GEOM_BDF_SAVE)
            model.write_bdf(fname)

    m = ModelSidebar(model)
    sys.exit(app.exec_())


if __name__ == "__main__":  # pragma: no cover
    main()
