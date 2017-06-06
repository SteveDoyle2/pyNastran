"""
Defines the GUI IO file for Cart3d.
"""
from __future__ import print_function
import os
from six import iteritems
from six.moves import range

from numpy import arange, mean, vstack, zeros, unique, where, sqrt
import numpy as np

import vtk
from vtk import vtkTriangle

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.cart3d.cart3d import read_cart3d
from pyNastran.converters.cart3d.cart3d_result import Cart3dGeometry #, Cart3dResult

from pyNastran.converters.cart3d.input_c3d_reader import read_input_c3d
from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl


class Cart3dIO(object):
    def __init__(self):
        pass

    def get_cart3d_wildcard_geometry_results_functions(self):
        data = ('Cart3d',
                'Cart3d (*.tri; *.triq)', self.load_cart3d_geometry,
                'Cart3d (*.triq)', self.load_cart3d_results)
        return data

    def _remove_old_cart3d_geometry(self, filename):
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
                del self.iSubcaseNameMap
            except:
                # print("cant delete geo")
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalarBar.VisibilityOff()
        self.scalarBar.Modified()
        return skip_reading

    def load_cart3d_geometry(self, cart3d_filename, dirname, name='main', plot=True):
        skip_reading = self._remove_old_cart3d_geometry(cart3d_filename)
        if skip_reading:
            return

        self.eid_map = {}
        self.nid_map = {}
        model = read_cart3d(cart3d_filename, log=self.log, debug=False)
        self.model_type = 'cart3d'
        nodes = model.nodes
        elements = model.elements
        regions = model.regions
        loads = model.loads

        self.nNodes = model.npoints
        self.nElements = model.nelements

        grid = self.grid
        grid.Allocate(self.nElements, 1000)

        self.nid_map = {}
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                self.grid_result.InsertNextValue(nid * fraction)

        assert nodes is not None
        nnodes = nodes.shape[0]

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
        self._create_cart3d_free_edges(model, nodes, elements)


        # loadCart3dResults - regions/loads
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        assert loads is not None
        if 'Mach' in loads:
            avg_mach = mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''
        self.iSubcaseNameMap = {1: ['Cart3d%s' % note, '']}
        cases = {}
        ID = 1
        form, cases, icase = self._fill_cart3d_case2(cases, ID, nodes, elements, regions, model)
        mach, alpha, beta = self._create_box(cart3d_filename, ID, form, cases, icase, regions)
        #mach = None
        self._fill_cart3d_results(cases, form, icase, ID, loads, model, mach)
        self._finish_results_io2(form, cases)

    def _create_box(self, cart3d_filename, ID, form, cases, icase, regions):
        stack = True
        dirname = os.path.dirname(os.path.abspath(cart3d_filename))
        input_c3d_filename = os.path.join(dirname, 'input.c3d')
        input_cntl_filename = os.path.join(dirname, 'input.cntl')
        mach = None
        alpha = None
        beta = None
        if os.path.exists(input_cntl_filename):
            cntl = read_input_cntl(input_cntl_filename, log=self.log, debug=self.debug)
            mach, alpha, beta = cntl.get_flow_conditions()
            bcs = cntl.get_boundary_conditions()
            bc_xmin, bc_xmax, bc_ymin, bc_ymax, bc_xmin, bc_xmax, surfbcs = bcs
            #stack = False

            if surfbcs:
                bc_form = [
                    ('Rho', icase, []),
                    ('xVelocity', icase + 1, []),
                    ('yVelocity', icase + 2, []),
                    ('zVelocity', icase + 3, []),
                    ('Mach', icase + 4, []),
                    ('Pressure', icase + 5, []),
                ]
                icase += 5
                nelements = self.nElements
                rho = zeros(nelements, dtype='float32')
                xvel = zeros(nelements, dtype='float32')
                yvel = zeros(nelements, dtype='float32')
                zvel = zeros(nelements, dtype='float32')
                vel = zeros(nelements, dtype='float32')
                pressure = zeros(nelements, dtype='float32')

                uregions = set(unique(regions))
                surf_bc_regions = set(surfbcs.keys())
                invalid_regions = surf_bc_regions - uregions
                if len(invalid_regions) != 0:
                    assert len(invalid_regions) == 0, invalid_regions

                for bc_id, bc_values in sorted(iteritems(surfbcs)):
                    rhoi, xveli, yveli, zveli, pressi = bc_values
                    i = where(regions == bc_id)[0]
                    rho[i] = rhoi
                    xvel[i] = xveli
                    yvel[i] = yveli
                    zvel[i] = zveli
                    pressure[i] = pressi

                mach = sqrt(xvel ** 2 + yvel ** 2 + zvel ** 2)

                rho_res = GuiResult(ID, header='Rho', title='Rho',
                                    location='centroid', scalar=rho)
                xvel_res = GuiResult(ID, header='xVelocity', title='xVelocity',
                                     location='centroid', scalar=xvel)
                yvel_res = GuiResult(ID, header='yVelocity', title='yVelocity',
                                     location='centroid', scalar=yvel)
                zvel_res = GuiResult(ID, header='zVelocity', title='zVelocity',
                                     location='centroid', scalar=zvel)
                mach_res = GuiResult(ID, header='Mach', title='Mach',
                                     location='centroid', scalar=mach)
                pressure_res = GuiResult(ID, header='Pressure', title='Pressure',
                                         location='centroid', scalar=pressure)

                cases[icase] = (rho_res, (ID, 'Rho'))
                cases[icase + 1] = (xvel_res, (ID, 'xVelocity'))
                cases[icase + 2] = (yvel_res, (ID, 'yVelocity'))
                cases[icase + 3] = (zvel_res, (ID, 'zVelocity'))
                cases[icase + 4] = (mach_res, (ID, 'Mach'))
                cases[icase + 5] = (pressure_res, (ID, 'Pressure'))
                form.append(('Boundary Conditions', None, bc_form))
        else:
            self.log.warning('input_cntl_filename doesnt exist = %s' % input_cntl_filename)


        if os.path.exists(input_c3d_filename):
            nodes, elements = read_input_c3d(input_c3d_filename, stack=stack,
                                             log=self.log, debug=self.debug)

            # Planes
            # ----------
            # xmin, xmax
            # ymin, ymax
            # zmin, zmax

            if stack:
                red = (1., 0., 0.)
                color = red
                self.set_quad_grid('box', nodes, elements, color, line_width=1, opacity=1.)
            else:
                red = (1., 0., 0.)
                inflow_nodes = []
                inflow_elements = []

                green = (0., 1., 0.)
                symmetry_nodes = []
                symmetry_elements = []

                colori = (1., 1., 0.)
                outflow_nodes = []
                outflow_elements = []

                blue = (0., 0., 1.)
                farfield_nodes = []
                farfield_elements = []

                ifarfield = 0
                isymmetry = 0
                iinflow = 0
                ioutflow = 0

                nfarfield_nodes = 0
                nsymmetry_nodes = 0
                ninflow_nodes = 0
                noutflow_nodes = 0
                for bcsi, nodesi, elementsi in zip(bcs, nodes, elements):
                    # 0 = FAR FIELD
                    # 1 = SYMMETRY
                    # 2 = INFLOW  (specify all)
                    # 3 = OUTFLOW (simple extrap)
                    self.log.info('bcsi = %s' % bcsi)
                    nnodes = nodesi.shape[0]
                    bc = bcsi
                    if isinstance(bc, int):
                        if bc == 0:
                            farfield_nodes.append(nodesi)
                            farfield_elements.append(elementsi + nfarfield_nodes)
                            nfarfield_nodes += nnodes
                            ifarfield += 1
                        elif bc == 1:
                            symmetry_nodes.append(nodesi)
                            symmetry_elements.append(elementsi + nsymmetry_nodes)
                            nsymmetry_nodes += nnodes
                            isymmetry += 1
                        elif bc == 2:
                            inflow_nodes.append(nodesi)
                            inflow_elements.append(elementsi + ninflow_nodes)
                            ninflow_nodes += nnodes
                            iinflow += 1
                        elif bc == 3:
                            outflow_nodes.append(nodesi)
                            outflow_elements.append(elementsi + noutflow_nodes)
                            noutflow_nodes += nnodes
                            ioutflow += 1
                        else:
                            msg = 'bc=%s' % str(bc)
                            raise NotImplementedError(msg)
                    elif isinstance(bc, dict):
                        continue
                    else:
                        msg = 'bc=%s' % str(bc)
                        raise NotImplementedError(msg)

                if ifarfield:
                    color = blue
                    nodes = vstack(farfield_nodes)
                    elements = vstack(farfield_elements)
                    self.set_quad_grid('farfield', nodes, elements, color, line_width=1, opacity=1.)

                if isymmetry:
                    color = green
                    nodes = vstack(symmetry_nodes)
                    elements = vstack(symmetry_elements)
                    self.set_quad_grid('symmetry', nodes, elements, color, line_width=1, opacity=1.)

                if iinflow:
                    color = red
                    nodes = vstack(inflow_nodes)
                    elements = vstack(inflow_elements)
                    self.set_quad_grid('inflow', nodes, elements, color, line_width=1, opacity=1.)

                if ioutflow:
                    color = colori
                    nodes = vstack(outflow_nodes)
                    elements = vstack(outflow_elements)
                    self.set_quad_grid('outflow', nodes, elements, color, line_width=1, opacity=1.)

                #i = 0
                #for nodesi, elementsi in zip(nodes, elements):
                    #self.set_quad_grid('box_%i' % i, nodesi, elementsi, color,
                                       #line_width=1, opacity=1.)
                    #i += 1
        else:
            self.log.warning('input_c3d_filename doesnt exist = %s' % input_c3d_filename)
        return mach, alpha, beta

    def _create_cart3d_free_edges(self, model, nodes, elements):
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
            elements2 = np.arange(0, nfree_edges * 2, dtype='int32').reshape(nfree_edges, 2)
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

    def clear_cart3d(self):
        pass

    def load_cart3d_results(self, cart3d_filename, dirname):
        self.load_cart3d_geometry(cart3d_filename, dirname)

    def _fill_cart3d_case2(self, cases, ID, nodes, elements, regions, model):
        nelements = elements.shape[0]
        nnodes = nodes.shape[0]

        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)
        area = model.get_area()
        cnormals = model.get_normals()
        cnnodes = cnormals.shape[0]
        assert cnnodes == nelements, len(cnnodes)

        #print('nnodes =', nnodes)
        #print('nelements =', nelements)
        #print('regions.shape =', regions.shape)
        subcase_id = 0
        labels = ['NodeID', 'ElementID', 'Region', 'Area',
                  'Normal X', 'Normal Y', 'Normal Z']
        cart3d_geo = Cart3dGeometry(subcase_id, labels,
                                    nids, eids, regions, area, cnormals,
                                    uname='Cart3dGeometry')

        cases = {
            0 : (cart3d_geo, (0, 'NodeID')),
            1 : (cart3d_geo, (0, 'ElementID')),
            2 : (cart3d_geo, (0, 'Region')),
            3 : (cart3d_geo, (0, 'Area')),
            4 : (cart3d_geo, (0, 'NormalX')),
            5 : (cart3d_geo, (0, 'NormalY')),
            6 : (cart3d_geo, (0, 'NormalZ')),
        }
        geometry_form = [
            ('NodeID', 0, []),
            ('ElementID', 1, []),
            ('Region', 2, []),
            ('Area', 3, []),
            ('Normal X', 4, []),
            ('Normal Y', 5, []),
            ('Normal Z', 6, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]
        icase = 7
        return form, cases, icase

        #cnormals = model.get_normals(nodes, elements)
        #nnormals = model.get_normals_at_nodes(nodes, elements, cnormals)

        #cases_new[i] = (ID, nnormals[:, 0], 'Normal X', 'node', '%.3f')
        #cases_new[i + 1] = (ID, nnormals[:, 1], 'Normal Y', 'node', '%.3f')
        #cases_new[i + 2] = (ID, nnormals[:, 2], 'Normal Z', 'node', '%.3f')
        #i += 3

    def _fill_cart3d_results(self, cases, form, icase, ID, loads, model, mach):
        new = False
        results_form = []
        cases_new = []
        result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                        'rhoU', 'rhoV', 'rhoW', 'rhoE', 'a', 'T', 'q', 'Pressure']

        for result_name in result_names:
            #print('result_name = %r' % result_name)
            if result_name in loads:
                nodal_data = loads[result_name]

                rho_res = GuiResult(ID, header=result_name, title=result_name,
                                    location='node', scalar=nodal_data)
                cases[icase] = (rho_res, (0, result_name))
                results_form.append((result_name, icase, []))
                icase += 1

        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases, icase
