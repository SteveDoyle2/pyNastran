"""Defines the GUI IO file for OBJ."""
import numpy as np

from pyNastran.gui.vtk_interface import vtkTriangle, vtkQuad
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points
from pyNastran.converters.dev.obj.obj import read_obj


class ObjIO:
    """Defines the GUI class for OBJ."""
    def __init__(self, gui):
        self.gui = gui

    def get_obj_wildcard_geometry_results_functions(self):
        """gets the OBJ wildcard loader used in the file load menu"""
        data = ('obj',
                'OBJ (*.obj)', self.load_obj_geometry,
                None, None)
        return data

    def _remove_old_obj_geometry(self, filename):
        #return self._remove_old_geometry(filename)

        self.gui.eid_map = {}
        self.gui.nid_map = {}
        if filename is None:
            self.gui.scalar_bar_actor.VisibilityOff()
            skip_reading = True
        else:
            self.gui.turn_corner_text_off()
            self.gui.grid.Reset()

            self.gui.result_cases = {}
            self.gui.ncases = 0
            try:
                del self.gui.case_keys
                del self.gui.icase
                del self.gui.isubcase_name_map
            except Exception:
                # print("cant delete geo")
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()
        return skip_reading

    def load_obj_geometry(self, obj_filename, name='main', plot=True):
        """
        The entry point for OBJ geometry loading.

        Parameters
        ----------
        obj_filename : str
            the obj filename to load
        name : str
            the name of the "main" actor for the GUI
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded
        """
        model_name = name
        skip_reading = self._remove_old_obj_geometry(obj_filename)
        if skip_reading:
            return

        log = self.gui.log
        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}
        model = read_obj(obj_filename, log=log, debug=False)
        self.model_type = 'obj'
        nodes = model.nodes
        nelements = model.nelements

        self.gui.nnodes = model.nnodes
        self.gui.nelements = nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        assert nodes is not None
        #nnodes = nodes.shape[0]

        mmax = nodes.max(axis=0)
        mmin = nodes.min(axis=0)
        dim_max = (mmax - mmin).max()
        xmax, ymax, zmax = mmax
        xmin, ymin, zmin = mmin
        log.info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        log.info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        log.info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        self.gui.create_global_axes(dim_max)
        points = numpy_to_vtk_points(nodes)

        #assert elements.min() == 0, elements.min()

        tri_etype = 5 # vtkTriangle().GetCellType()
        #self.create_vtk_cells_of_constant_element_type(grid, elements, etype)
        quad_etype = 9 # vtkQuad().GetCellType()

        all_tris = [tri_faces for tri_faces in model.tri_faces.values()]
        all_quads = [quad_faces for quad_faces in model.quad_faces.values()]
        if len(all_tris):
            tris = np.vstack(all_tris)
            for eid, element in enumerate(tris):
                elem = vtkTriangle()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                grid.InsertNextCell(tri_etype, elem.GetPointIds())
        if len(all_quads):
            quads = np.vstack(all_quads)
            for eid, element in enumerate(quads):
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                elem.GetPointIds().SetId(3, element[3])
                grid.InsertNextCell(quad_etype, elem.GetPointIds())
        grid.SetPoints(points)
        grid.Modified()

        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['OBJ', '']}
        cases = {}
        ID = 1
        form, cases, icase, node_ids, element_ids = self._fill_obj_geometry_objects(
            cases, ID, nodes, nelements, model)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    def clear_obj(self):
        pass

    def _fill_obj_geometry_objects(self, cases, ID, nodes, nelements, model):
        nnodes = nodes.shape[0]

        eids = np.arange(1, nelements + 1)
        nids = np.arange(1, nnodes + 1)

        subcase_id = 0
        nid_res = GuiResult(subcase_id, 'NodeID', 'NodeID', 'node', nids)
        eid_res = GuiResult(subcase_id, 'ElementID', 'ElementID', 'centroid', eids)


        cases[0] = (nid_res, (0, 'NodeID'))
        cases[1] = (eid_res, (0, 'ElementID'))
            #2 : (area_res, (0, 'Area')),
            #4 : (cart3d_geo, (0, 'NormalX')),
            #5 : (cart3d_geo, (0, 'NormalY')),
            #6 : (cart3d_geo, (0, 'NormalZ')),

        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            #('Area', 2, []),
            #('Normal X', 4, []),
            #('Normal Y', 5, []),
            #('Normal Z', 6, []),
        ]

        form = [
            ('Geometry', None, geometry_form),
        ]
        icase = 2
        return form, cases, icase, nids, eids
