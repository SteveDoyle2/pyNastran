"""Defines the GUI IO file for Cart3d."""
import os
from collections import OrderedDict
import collections

from numpy import arange, mean, vstack, unique, where, sqrt
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.qt_files.colors import (
    RED_FLOAT, BLUE_FLOAT, GREEN_FLOAT, PINK_FLOAT, YELLOW_FLOAT)
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)
from pyNastran.converters.cart3d.cart3d import read_cart3d
from pyNastran.converters.cart3d.cart3d_result import Cart3dGeometry

from pyNastran.converters.cart3d.input_c3d_reader import read_input_c3d
from pyNastran.converters.cart3d.input_cntl_reader import read_input_cntl


class Cart3dIO:
    """Defines the GUI class for Cart3d."""
    def __init__(self, gui):
        self.gui = gui
        self.data_map = None

    def get_cart3d_wildcard_geometry_results_functions(self):
        """
        gets the Cart3d wildcard loader used in the file load menu
        """
        data = ('Cart3d',
                'Cart3d (*.tri; *.triq)', self.load_cart3d_geometry,
                'Cart3d (*.triq)', self.load_cart3d_results)
        return data

    def _remove_old_cart3d_geometry(self, filename, model_name):
        #return
        #if model_name != 'main':
            #return
        #return self._remove_old_geometry(filename)

        self.gui.eid_map = {}
        self.gui.nid_map = {}
        if filename is None:
            self.gui.scalar_bar_actor.VisibilityOff()
            skip_reading = True
        else:
            self.gui.turn_text_off()
            self.gui.grid.Reset()

            self.gui.result_cases = OrderedDict()
            self.gui.ncases = 0
            try:
                del self.gui.case_keys
                del self.gui.icase
                del self.gui.isubcase_name_map
            except:
                # print("cant delete geo")
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()
        return skip_reading

    def load_cart3d_geometry(self, cart3d_filename, name='main', plot=True):
        """
        The entry point for Cart3d geometry loading.

        Parameters
        ----------
        cart3d_filename : str
            the cart3d filename to load
        name : str
            the name of the "main" actor for the GUI
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded
        """
        model_name = name
        skip_reading = self._remove_old_cart3d_geometry(cart3d_filename, model_name)
        if skip_reading:
            return

        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}
        model = read_cart3d(cart3d_filename, log=self.gui.log, debug=False)
        self.model = model
        self.gui.model_type = 'cart3d'
        nodes = model.nodes
        elements = model.elements
        regions = model.regions
        loads = model.loads

        self.gui.nnodes = model.npoints
        self.gui.nelements = model.nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        #if 0:
            #fraction = 1. / self.nnodes  # so you can color the nodes by ID
            #for nid, node in sorted(nodes.items()):
                #self.grid_result.InsertNextValue(nid * fraction)

        assert nodes is not None
        #nnodes = nodes.shape[0]

        mmax = nodes.max(axis=0)
        mmin = nodes.min(axis=0)
        dim_max = (mmax - mmin).max()
        xmax, ymax, zmax = mmax
        xmin, ymin, zmin = mmin
        self.gui.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        self.gui.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        self.gui.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        self.gui.create_global_axes(dim_max)
        points = numpy_to_vtk_points(nodes)

        #assert elements.min() == 0, elements.min()

        etype = 5 # vtkTriangle().GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, etype)

        grid.SetPoints(points)
        grid.Modified()
        self._create_cart3d_free_edges(model, nodes, elements)


        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        assert loads is not None
        if 'Mach' in loads:
            avg_mach = mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''
        self.gui.isubcase_name_map = {1: ['Cart3d%s' % note, '']}
        cases = OrderedDict()
        ID = 1

        icase = self.gui.get_new_icase()
        form, cases, icase, node_ids, element_ids, data_map_dict = _fill_cart3d_geometry_objects(
            cases, ID, nodes, elements, regions, model, model_name, icase=icase)
        if model_name != 'main':
            #print('icase1a =', len(self.gui.result_cases))
            #print('icase1b =', len(cases))
            self.gui.update_result_cases(cases)
            self.gui.result_cases.update(cases)
            form = self.gui.get_form() + form

        self.data_map = data_map_dict

        mach, unused_alpha, unused_beta = self._create_box(
            cart3d_filename, ID, form, cases, icase, regions, model_name)
        #mach = None
        _fill_cart3d_results(cases, form, icase, ID, loads, model, mach)

        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        build_map_centroidal_result(model)
        self.gui._finish_results_io2(model_name, form, cases)

    def _create_box(self, cart3d_filename, ID, form, cases, icase, regions, model_name):
        """creates the bounding box for boundary conditions"""
        dirname = os.path.dirname(os.path.abspath(cart3d_filename))
        input_c3d_filename = os.path.join(dirname, 'input.c3d')
        input_cntl_filename = os.path.join(dirname, 'input.cntl')
        mach = None
        alpha = None
        beta = None
        unused_gamma = None

        bcs = None
        if os.path.exists(input_cntl_filename):
            cntl = read_input_cntl(input_cntl_filename,
                                   log=self.gui.log, debug=self.gui.debug)
            mach, alpha, beta, unused_gamma = cntl.get_flow_conditions()
            (unused_bc_xmin, unused_bc_xmax,
             unused_bc_ymin, unused_bc_ymax,
             unused_bc_xmin, unused_bc_xmax, surfbcs) = cntl.get_boundary_conditions()
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
                nelements = self.gui.nelements
                rho = np.full(nelements, np.nan, dtype='float32')
                xvel = np.full(nelements, np.nan, dtype='float32')
                yvel = np.full(nelements, np.nan, dtype='float32')
                zvel = np.full(nelements, np.nan, dtype='float32')
                #vel = np.full(nelements, np.nan, dtype='float32')
                pressure = np.full(nelements, np.nan, dtype='float32')

                uregions = set(unique(regions))
                surf_bc_regions = set(surfbcs.keys())
                invalid_regions = surf_bc_regions - uregions
                if len(invalid_regions) != 0:
                    assert len(invalid_regions) == 0, invalid_regions

                for bc_id, bc_values in sorted(surfbcs.items()):
                    rhoi, xveli, yveli, zveli, pressi = bc_values
                    i = where(regions == bc_id)[0]
                    rho[i] = rhoi
                    xvel[i] = xveli
                    yvel[i] = yveli
                    zvel[i] = zveli
                    pressure[i] = pressi

                inan = np.where(rho == 0.0)
                rho[inan] = np.nan
                xvel[inan] = np.nan
                yvel[inan] = np.nan
                zvel[inan] = np.nan
                #vel[inan] = np.nan
                pressure[inan] = np.nan

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
                form.append((get_name('Boundary Conditions', model_name), None, bc_form))
        else:
            self.gui.log.warning('input_cntl_filename doesnt exist = %s' % input_cntl_filename)


        if os.path.exists(input_c3d_filename):
            # put in one group

            # Planes
            # ----------
            # xmin, xmax
            # ymin, ymax
            # zmin, zmax
            nodes, elements = read_input_c3d(input_c3d_filename, stack=True,
                                             log=self.gui.log, debug=self.gui.debug)

            #print('name=%r model_name=%r' % (self.gui.name, model_name))
            box_name = get_name('box', self.gui.name)
            self.gui.set_quad_grid(box_name, nodes, elements,
                                   color=RED_FLOAT, line_width=1, opacity=1.)

            #-------------------------------------------------------------------
            # put in multiple groups
            nodes, elements = read_input_c3d(input_c3d_filename, stack=False,
                                             log=self.gui.log, debug=self.gui.debug)

            inflow_nodes = []
            inflow_elements = []

            symmetry_nodes = []
            symmetry_elements = []

            outflow_nodes = []
            outflow_elements = []

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
            if bcs is None:
                bcs = [None] * len(nodes)
            for bcsi, nodesi, elementsi in zip(bcs, nodes, elements):
                # 0 = FAR FIELD
                # 1 = SYMMETRY
                # 2 = INFLOW  (specify all)
                # 3 = OUTFLOW (simple extrap)
                self.gui.log.info('bcsi = %s' % bcsi)
                nnodes = nodesi.shape[0]
                bc = bcsi
                if bc is None:  # fake case
                    continue
                elif isinstance(bc, integer_types):
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
                elif isinstance(bc, dict): # ???
                    if len(bc) == 0:
                        continue
                    # bc = {
                    #    2: [2.0, 3.0, 0.0, 0.0, 5.0],
                    #    3: [1.0, 1.5, 0.0, 0.0, 0.714285]
                    # }
                    continue
                    #msg = 'bc=%s' % str(bc)
                    #raise NotImplementedError(msg)
                else:
                    msg = 'bc=%s' % str(bc)
                    raise NotImplementedError(msg)

            if ifarfield:
                nodes = vstack(farfield_nodes)
                elements = vstack(farfield_elements)
                name = get_name('farfield', self.gui.name)
                self.gui.set_quad_grid(name, nodes, elements, color=BLUE_FLOAT,
                                       line_width=1, opacity=1.)

            if isymmetry:
                nodes = vstack(symmetry_nodes)
                elements = vstack(symmetry_elements)
                name = get_name('symmetry', self.gui.name)
                self.gui.set_quad_grid(name, nodes, elements, color=GREEN_FLOAT,
                                       line_width=1, opacity=1.)

            if iinflow:
                nodes = vstack(inflow_nodes)
                elements = vstack(inflow_elements)
                name = get_name('inflow', self.gui.name)
                self.gui.set_quad_grid(name, nodes, elements, color=RED_FLOAT,
                                       line_width=1, opacity=1.)

            if ioutflow:
                nodes = vstack(outflow_nodes)
                elements = vstack(outflow_elements)
                name = get_name('outflow', self.gui.name)
                self.gui.set_quad_grid(name, nodes, elements, color=YELLOW_FLOAT,
                                       line_width=1, opacity=1.)

            #i = 0
            #for nodesi, elementsi in zip(nodes, elements):
                #self.set_quad_grid('box_%i' % i, nodesi, elementsi, color,
                                   #line_width=1, opacity=1.)
                #i += 1
        else:
            self.gui.log.warning('input_c3d_filename doesnt exist = %s' % input_c3d_filename)
        return mach, alpha, beta

    def _create_cart3d_free_edges(self, model, nodes, elements):
        """creates the free edges to help identify unclosed models"""
        free_edges_array = model.get_free_edges(elements)
        nfree_edges = len(free_edges_array)

        name = get_name('free_edges', self.gui.name)
        if nfree_edges:
            npoints = 2 * nfree_edges
            if name not in self.gui.alt_grids:
                self.gui.create_alternate_vtk_grid(
                    name, color=PINK_FLOAT, line_width=3, opacity=1.0,
                    representation='surface')

            alt_grid = self.gui.alt_grids[name]
            etype = 3  # vtk.vtkLine().GetCellType()
            elements2 = np.arange(0, npoints, dtype='int32').reshape(nfree_edges, 2)
            create_vtk_cells_of_constant_element_type(alt_grid, elements2, etype)

            #alt_grid.Allocate(nfree_edges, 1000)
            free_edge_nodes = nodes[free_edges_array.ravel(), :]
            points = numpy_to_vtk_points(free_edge_nodes)
            alt_grid.SetPoints(points)

        else:
            # TODO: clear free edges
            pass

        if name in self.gui.alt_grids:
            self.gui._add_alt_actors(self.gui.alt_grids)
            self.gui.geometry_actors[name].Modified()
            if hasattr(self.gui.geometry_actors[name], 'Update'):
                self.gui.geometry_actors[name].Update()

    def clear_cart3d(self):
        pass

    def load_cart3d_results(self, cart3d_filename):
        """Loads the Cart3d results into the GUI"""
        self.load_cart3d_geometry(cart3d_filename)

def _fill_cart3d_geometry_objects(cases, unused_id, nodes, elements, regions, model,
                                  model_name, icase=0):
    """Creates the results form for Cart3d Geometry"""
    nelements = elements.shape[0]
    nnodes = nodes.shape[0]

    eids = arange(1, nelements + 1)
    nids = arange(1, nnodes + 1)
    area = model.get_area()
    cnormals = model.get_normals()
    cnnodes = cnormals.shape[0]
    assert cnnodes == nelements, len(cnnodes)

    inv_counter = _node_inverse_counter(model, nnodes)
    def data_map_func(data):
        res = np.zeros(nnodes, dtype='float32')
        for elem, datai in zip(model.elements, data):
            res[elem] += datai
        #print(res)

        #res = np.zeros(nnodes, dtype='float32')
        #res[model.elements] = data
        #print(res)
        #print('----')
        #print(inv_counter)
        ## neids * ??? = (3,nnodes)
        ##(3,) * ??? = (5,)
        ## eids * map -> node_ids
        ##results = np.zeros(nnodes, dtype='float32')
        #print(model.elements)
        #node_results = data[model.elements.ravel()]
        #assert node_results.shape == model.elements.shape
        #node_results_sum = node_results.sum(axis=1)
        #assert node_results_sum.shape == nnodes
        #return node_results_sum * inv_counter
        return res * inv_counter

    data_map_dict = {('centroid', 'Node') : data_map_func}

    #print('nnodes =', nnodes)
    #print('nelements =', nelements)
    #print('regions.shape =', regions.shape)
    subcase_id = 0
    labels = ['NodeID', 'ElementID', 'Region', 'Area',
              'Normal X', 'Normal Y', 'Normal Z']
    cart3d_geo = Cart3dGeometry(subcase_id, labels,
                                nids, eids, regions, area, cnormals,
                                uname='Cart3dGeometry')

    #normal_z = cnormals[:, 2]
    #node_normal = data_map_func(normal_z)
    #result_name = 'Normal Z-nodal'
    #node_res = GuiResult(subcase_id, header=result_name, title=result_name,
                         #location='node', scalar=node_normal)

    cases = OrderedDict()
    cases[icase + 0] = (cart3d_geo, (0, 'NodeID'))
    cases[icase + 1] = (cart3d_geo, (0, 'ElementID'))
    cases[icase + 2] = (cart3d_geo, (0, 'Region'))
    cases[icase + 3] = (cart3d_geo, (0, 'Area'))
    cases[icase + 4] = (cart3d_geo, (0, 'NormalX'))
    cases[icase + 5] = (cart3d_geo, (0, 'NormalY'))
    cases[icase + 6] = (cart3d_geo, (0, 'NormalZ'))
    #cases[icase + 7] = (node_res, (0, 'NormalZ-nodal'))

    geometry_form = [
        ('NodeID', icase + 0, []),
        ('ElementID', icase + 1, []),
        ('Region', icase + 2, []),
        ('Area', icase + 3, []),
        ('Normal X', icase + 4, []),
        ('Normal Y', icase + 5, []),
        ('Normal Z', icase + 6, []),
        #('Normal Z-nodal', icase + 7, []),
    ]
    form = [
        (get_name('Geometry', model_name), None, geometry_form),
    ]
    icase += len(geometry_form)
    return form, cases, icase, nids, eids, data_map_dict
    #cnormals = model.get_normals(nodes, elements)
    #nnormals = model.get_normals_at_nodes(nodes, elements, cnormals)


def _node_inverse_counter(model, nnodes):
    node_ids = model.elements.ravel()
    #max_nid = node_ids.max()
    node_count = collections.Counter(node_ids)
    #Counter({0: 7, 1: 4, 3: 2, 2: 1, 4: 1})

    # we're going to multiply by 1/node_count
    # if we have no nodes, then we have no results, so we have a sum of 0
    # so we make it 1 to prevent division by 0
    counter = np.ones(nnodes)
    for nid, count in node_count.items():
        counter[nid] = count
    inv_counter = 1. / counter
    return inv_counter

def _fill_cart3d_results(cases, form, icase, ID, loads, unused_model, unused_mach):
    """Creates the results form for Cart3d Results"""
    results_form = []
    unused_cases_new = []
    result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                    'rhoU', 'rhoV', 'rhoW', 'rhoE', 'a', 'T', 'q', 'Pressure']

    inan = None
    if 'rho' in loads:
        rho = loads['rho']
        inan = np.where(rho == 0.)
    for result_name in result_names:
        #print('result_name = %r' % result_name)
        if result_name in loads:
            nodal_data = loads[result_name]
            if inan is not None:
                nodal_data[inan] = np.nan
            rho_res = GuiResult(ID, header=result_name, title=result_name,
                                location='node', scalar=nodal_data)
            cases[icase] = (rho_res, (0, result_name))
            results_form.append((result_name, icase, []))
            icase += 1

    if len(results_form):
        form.append(('Results', None, results_form))
    return form, cases, icase

def build_map_centroidal_result(model):
    """
    Sets up map_centroidal_result.  Used for:
     - cutting plane
     - nodal Cp
    """
    if hasattr(model, 'map_centroidal_result'):
        return
    #mapped_node_ids = []
    nnodes = model.nnodes

    elem_ravel = model.elements.flatten()
    unused_node_count_ids, node_count = np.unique(elem_ravel, return_counts=True)
    #print('node_count_ids =', node_count_ids)
    #print('node_count =', node_count)

    # tehre are some nodes that are used 118 times...
    # this can be at a converging tip
    #ibig = np.where(node_count > 7)
    #print(node_count_ids[ibig])
    #print(node_count[ibig])

    #print(str(node_count.tolist()).replace('L', ''))

    # calculate inv_node_count
    inv_node_count = 1. / node_count

    # build the centroidal mapper
    def map_centroidal_result(centroidal_data):
        """maps centroidal data onto nodal data"""
        nodal_data = np.zeros(nnodes, dtype=centroidal_data.dtype)
        for datai, node_ids in zip(centroidal_data, model.elements):
            for nid in node_ids:
                nodal_data[nid] += datai
        return nodal_data * inv_node_count
    model.map_centroidal_result = map_centroidal_result

def get_name(base, gui_name):
    return base if gui_name == 'main' else '%s - %s' % (base, gui_name)
