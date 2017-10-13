"""
Defines the GUI IO file for OBJ.
"""
from __future__ import print_function
import os
from six import iteritems
from six.moves import range

from numpy import arange, mean, vstack, unique, where, sqrt
import numpy as np

from pyNastran.utils import integer_types
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.dev.obj.obj import read_obj


class ObjIO(object):
    """
    Defines the GUI class for OBJ.
    """
    def __init__(self):
        pass

    def get_obj_wildcard_geometry_results_functions(self):
        """
        gets the OBJ wildcard loader used in the file load menu
        """
        data = ('OBJ',
                'OBJ (*.obj)', self.load_obj_geometry,
                None)
        return data

    def _remove_old_obj_geometry(self, filename):
        #return self._remove_old_geometry(filename)

        self.eid_map = {}
        self.nid_map = {}
        if filename is None:
            #self.emptyResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.scalarBar.VisibilityOff()
            skip_reading = True
        else:
            self.turn_text_off()
            self.grid.Reset()
            #self.gridResult.Reset()
            #self.gridResult.Modified()

            self.result_cases = {}
            self.ncases = 0
            try:
                del self.case_keys
                del self.icase
                del self.isubcase_name_map
            except:
                # print("cant delete geo")
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
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
        skip_reading = self._remove_old_obj_geometry(obj_filename)
        if skip_reading:
            return

        self.eid_maps[name] = {}
        self.nid_maps[name] = {}
        model = read_obj(obj_filename, log=self.log, debug=False)
        self.model_type = 'obj'
        nodes = model.nodes
        elements = model.tri_faces #elements

        self.nnodes = model.nnodes
        self.nelements = model.nelements

        grid = self.grid
        grid.Allocate(self.nelements, 1000)

        #if 0:
            #fraction = 1. / self.nnodes  # so you can color the nodes by ID
            #for nid, node in sorted(iteritems(nodes)):
                #self.grid_result.InsertNextValue(nid * fraction)

        assert nodes is not None
        #nnodes = nodes.shape[0]

        mmax = nodes.max(axis=0)
        mmin = nodes.min(axis=0)
        dim_max = (mmax - mmin).max()
        xmax, ymax, zmax = mmax
        xmin, ymin, zmin = mmin
        self.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        self.create_global_axes(dim_max)
        points = self.numpy_to_vtk_points(nodes)

        #assert elements.min() == 0, elements.min()

        etype = 5 # vtkTriangle().GetCellType()
        self.create_vtk_cells_of_constant_element_type(grid, elements, etype)

        grid.SetPoints(points)
        grid.Modified()
        if hasattr(grid, 'Update'):
            grid.Update()
        self._create_obj_free_edges(model, nodes, elements)


        # loadCart3dResults - regions/loads
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.isubcase_name_map = {1: ['OBJ', '']}
        cases = {}
        ID = 1
        form, cases, icase = self._fill_obj_geometry_objects(
            cases, ID, nodes, elements, model)
        self._finish_results_io2(form, cases)

    def _create_obj_free_edges(self, model, nodes, elements):
        """creates the free edges to help identify unclosed models"""
        return
        free_edges_array = model.get_free_edges(elements)
        nfree_edges = len(free_edges_array)

        if nfree_edges:
            # yellow = (1., 1., 0.)
            pink = (0.98, 0.4, 0.93)
            npoints = 2 * nfree_edges
            if 'free_edges' not in self.alt_grids:
                self.create_alternate_vtk_grid(
                    'free_edges', color=pink, line_width=3, opacity=1.0,
                    representation='surface')

            alt_grid = self.alt_grids['free_edges']
            etype = 3  # vtk.vtkLine().GetCellType()
            elements2 = np.arange(0, npoints, dtype='int32').reshape(nfree_edges, 2)
            self.create_vtk_cells_of_constant_element_type(alt_grid, elements2, etype)

            #alt_grid.Allocate(nfree_edges, 1000)
            free_edge_nodes = nodes[free_edges_array.ravel(), :]
            points = self.numpy_to_vtk_points(free_edge_nodes)
            alt_grid.SetPoints(points)

        else:
            # TODO: clear free edges
            pass

        if 'free_edges' in self.alt_grids:
            self._add_alt_actors(self.alt_grids)
            self.geometry_actors['free_edges'].Modified()
            if hasattr(self.geometry_actors['free_edges'], 'Update'):
                self.geometry_actors['free_edges'].Update()

    def clear_obj(self):
        pass

    def _fill_obj_geometry_objects(self, cases, ID, nodes, elements, model):
        nelements = elements.shape[0]
        nnodes = nodes.shape[0]

        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)
        #area = model.get_area()
        #cnormals = model.get_normals()
        #cnnodes = cnormals.shape[0]
        #assert cnnodes == nelements, len(cnnodes)

        #print('nnodes =', nnodes)
        #print('nelements =', nelements)
        #print('regions.shape =', regions.shape)
        subcase_id = 0
        labels = ['NodeID', 'ElementID', 'Area',
                  'Normal X', 'Normal Y', 'Normal Z']
        nid_res = GuiResult(subcase_id, 'NodeID', 'NodeID', 'node', nids)
        eid_res = GuiResult(subcase_id, 'ElementID', 'ElementID', 'centroid', eids)
        #area_res = GuiResult(subcase_id, 'NodeID', 'NodeID', 'centroid', area)

        cases = {
            0 : (nid_res, (0, 'NodeID')),
            1 : (eid_res, (0, 'ElementID')),
            #2 : (area_res, (0, 'Area')),
            #4 : (cart3d_geo, (0, 'NormalX')),
            #5 : (cart3d_geo, (0, 'NormalY')),
            #6 : (cart3d_geo, (0, 'NormalZ')),
        }
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
        return form, cases, icase

