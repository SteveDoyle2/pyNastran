from pyNastran.dev.h5.h5_nastran2 import get_gui_nastran_ugrid
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.qt_files.colors import (
    RED_FLOAT, BLUE_FLOAT, GREEN_FLOAT, LIGHT_GREEN_FLOAT, PINK_FLOAT, PURPLE_FLOAT,
    YELLOW_FLOAT, ORANGE_FLOAT)

class H5NastranIO():
    def __init__(self, gui):
        self.gui = gui

    def get_h5nastran_wildcard_geometry_results_functions(self):
        data = ('VSPAero',
                'H5Nastran (*.h5)', self.load_h5nastran_geometry,
                None, None
               )
        return data

    def load_h5nastran_geometry(self, hdf5_filename, name='main', plot=True, **kwargs):
        out = get_gui_nastran_ugrid(
            hdf5_filename,
            self.gui.grid,
            add_property_info=True,
            add_material_info=True,
            subcases=None,  # default=None -> all
            modes=None, # default=None -> all
            results=None, # default=None -> all,
        )
        model, ugrid, root, alt_grids, node_ids, element_ids, form, cases = out
        self.node_ids = node_ids

        ugrid = alt_grids['main']
        del alt_grids['main']
        for name, ugrid in alt_grids.items():
            self.gui.create_alternate_vtk_grid(
                name, color=ORANGE_FLOAT, line_width=5, opacity=1., point_size=4,
                representation='point', follower_function=None)
        if self.gui.alt_grids:
            self.gui._add_alt_actors(self.gui.alt_grids)

        model_name = name
        self.gui.isubcase_name_map = {1: ['OpenVSP', '']}
        ID = 1
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids

        #form, cases = self._fill_cases(cases, ID, model)
        self.gui._finish_results_io2(model_name, form, cases)
        #self.gui.grid = ugrid

    def _fill_cases(cases, ID, model):
        pass
