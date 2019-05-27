"""Defines the GUI IO file for Panair."""
import os
from collections import OrderedDict

import numpy as np
from numpy import zeros, amax, amin, arange

import vtk
from vtk import vtkQuad

from pyNastran.converters.panair.panair_grid import PanairGrid
from pyNastran.converters.panair.agps import AGPS
from pyNastran.converters.panair.panair_out import read_panair_out

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


class PanairIO:
    """Defines the GUI class for Panair."""
    def __init__(self, gui):
        self.gui = gui
        self.colormap = 'viridis'
        self.elements = None

    def get_panair_wildcard_geometry_results_functions(self):
        data = ('Panair',
                'Panair (*.inp)', self.load_panair_geometry,
                #'Panair (*.agps);;Panair (*.out)',  self.load_panair_results)
                'Panair (*agps);;All files (*)', self.load_panair_results)
        return data

    def load_panair_geometry(self, panair_filename, name='main', plot=True):
        model_name = name
        self.gui.nid_map = {}
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self.gui._remove_old_geometry(panair_filename)
        if skip_reading:
            return

        model = PanairGrid(log=self.gui.log, debug=self.gui.debug)
        self.gui.model_type = model.model_type
        model.read_panair(panair_filename)
        self.gui.geom_model = model
        # get_wakes=True shows explicit wakes
        #
        # TODO: bad for results...what about just not adding it to the patches/bcs?
        nodes, elements, regions, kt, cp_norm = model.get_points_elements_regions(
            get_wakes=True)

        self.gui.nnodes = len(nodes)
        self.gui.nelements = len(elements)

        #print("nnodes = ",self.nnodes)
        #print("nelements = ", self.nelements)

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.gui.nnodes)

        assert len(nodes) > 0
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)
        points = numpy_to_vtk_points(nodes)

        assert len(elements) > 0
        elem = vtkQuad()
        quad_type = elem.GetCellType()
        create_vtk_cells_of_constant_element_type(grid, elements, quad_type)

        grid.SetPoints(points)
        grid.Modified()

        # loadPanairResults - regions/loads
        if plot:
            self.gui.scalar_bar_actor.VisibilityOn()
            self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['Panair', '']}
        cases = OrderedDict()
        ID = 1

        loads = []
        form, cases, node_ids, element_ids = self._fill_panair_geometry_case(
            cases, ID, nodes, elements, regions, kt, cp_norm, loads)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids

        #if plot:
        self.gui._finish_results_io2(model_name, form, cases)
        #else:
        #self._set_results([form], cases)

    def clear_panair(self):
        del self.elements

    def _fill_panair_geometry_case(self, cases, ID, nodes, elements, regions,
                                   kt, cp_norm, unused_loads):
        self.elements = elements
        colormap = 'jet' # self.colormap
        nnids = nodes.shape[0]
        neids = elements.shape[0]
        nids = arange(0., nnids, dtype='int32') + 1
        eids = arange(0., neids, dtype='int32') + 1

        # centroidal
        #nelements = len(elements)
        #print('nelements = ', nelements)
        #print('nnodes = ', nodes.shape)

        p1 = nodes[elements[:, 0]]
        p2 = nodes[elements[:, 1]]
        p3 = nodes[elements[:, 2]]
        p4 = nodes[elements[:, 3]]
        xyz_centroid = (p1 + p2 + p3 + p4) / 4.

        n = np.cross(p3 - p1, p4 - p2)
        n_norm = np.linalg.norm(n, axis=1)
        area = n_norm / 2.
        normal = n / n_norm[:, np.newaxis]

        itime = 0
        # ID, header, title, location, values, format, uname
        region_res = GuiResult(ID, 'Network', 'Network', 'centroid', regions,
                               data_format=None, colormap=colormap, uname='Network')
        eid_res = GuiResult(ID, 'ElementID', 'ElementID', 'centroid', eids,
                            data_format=None, colormap=colormap, uname='ElementID')
        nid_res = GuiResult(ID, 'NodeID', 'NodeID', 'node', nids,
                            data_format=None, colormap=colormap, uname='NodeID')
        area_res = GuiResult(ID, 'Area', 'Area', 'centroid', area,
                             data_format=None, colormap=colormap, uname='Area')

        nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                nlabels=2, labelsize=5, ncolors=2,
                                #colormap=colormap,
                                data_format='%.1f',
                                uname='NormalResult')
        nx_res = GuiResult(ID, 'normal_x', 'NormalX', 'centroid', normal[:, 0],
                           data_format='%.3f', colormap=colormap, uname='NormalX')
        ny_res = GuiResult(ID, 'normal_y', 'NormalY', 'centroid', normal[:, 1],
                           data_format='%.3f', colormap=colormap, uname='NormalY')
        nz_res = GuiResult(ID, 'normal_z', 'NormalZ', 'centroid', normal[:, 2],
                           data_format='%.3f', colormap=colormap, uname='NormalZ')
        cenx_res = GuiResult(ID, 'centroid_x', 'CentroidX', 'centroid', xyz_centroid[:, 0],
                             data_format=None, colormap=colormap, uname='CentroidX')
        ceny_res = GuiResult(ID, 'centroid_y', 'CentroidY', 'centroid', xyz_centroid[:, 1],
                             data_format=None, colormap=colormap, uname='CentroidY')
        cenz_res = GuiResult(ID, 'centroid_z', 'CentroidZ', 'centroid', xyz_centroid[:, 2],
                             data_format=None, colormap=colormap, uname='CentroidZ')

        kt_res = GuiResult(ID, 'Kt', 'Kt', 'centroid', kt,
                           data_format=None, colormap=colormap, uname='Kt')
        cp_res = GuiResult(ID, 'CpNorm', 'CpNorm', 'centroid', cp_norm,
                           data_format=None, colormap=colormap, uname='CpNorm')

        # nodal
        cenx_res = GuiResult(ID, 'node_x', 'NodeX', 'node', nodes[:, 0],
                             data_format='%.2f', colormap=colormap, uname='node_x')
        ceny_res = GuiResult(ID, 'node_y', 'NodeY', 'node', nodes[:, 1],
                             data_format='%.2f', colormap=colormap, uname='node_y')
        cenz_res = GuiResult(ID, 'node_z', 'NodeZ', 'node', nodes[:, 2],
                             data_format='%.2f', colormap=colormap, uname='node_z')

        icase = 0
        location_form = [
            ('Normal', icase + 4, []),
            ('normal_x', icase + 5, []),
            ('normal_y', icase + 6, []),
            ('normal_z', icase + 7, []),

            ('centroid_x', icase + 8, []),
            ('centroid_y', icase + 9, []),
            ('centroid_z', icase + 10, []),

            ('node_x', icase + 11, []),
            ('node_y', icase + 12, []),
            ('node_z', icase + 13, []),
        ]

        geometry_form = [
            ('Patch', icase, []),
            ('ElementID', icase + 1, []),
            ('NodeID', icase + 2, []),
            ('Area', icase + 3, []),
            ('Location', None, location_form),
            ('Kt', icase + 14, []),
            ('CpNorm', icase + 15, []),
        ]
        form = [
            ('Geometry', None, geometry_form),
        ]

        cases[icase] = (region_res, (itime, 'Patch'))
        cases[icase + 1] = (eid_res, (itime, 'ElementID'))
        cases[icase + 2] = (nid_res, (itime, 'NodeID'))
        cases[icase + 3] = (area_res, (itime, 'Area'))

        # location_form
        cases[icase + 4] = (nxyz_res, (itime, 'Normal'))
        cases[icase + 5] = (nx_res, (itime, 'NormalX'))
        cases[icase + 6] = (ny_res, (itime, 'NormalY'))
        cases[icase + 7] = (nz_res, (itime, 'NormalZ'))
        #---
        cases[icase + 8] = (cenx_res, (itime, 'CentroidX'))
        cases[icase + 9] = (ceny_res, (itime, 'CentroidY'))
        cases[icase + 10] = (cenz_res, (itime, 'CentroidZ'))
        #---
        cases[icase + 11] = (cenx_res, (itime, 'node_x'))
        cases[icase + 12] = (ceny_res, (itime, 'node_y'))
        cases[icase + 13] = (cenz_res, (itime, 'node_z'))

        cases[icase + 14] = (kt_res, (itime, 'Kt'))
        cases[icase + 15] = (cp_res, (itime, 'CpNorm'))

        return form, cases, nids, eids

    def load_panair_results(self, panair_filename):
        model_name = 'main'
        #colormap = self.colormap
        colormap = 'jet'
        cases = self.gui.result_cases
        model = AGPS(log=self.gui.log, debug=self.gui.debug)
        model.read_agps(panair_filename)

        icase = len(self.gui.result_cases)
        results_form = []

        # get the Cp on the nodes
        cp0 = model.pressures[0]
        ncp = cp0.shape[0]

        nnodes = self.gui.nnodes
        nelements = self.gui.nelements

        geom_model = self.gui.geom_model
        mach = geom_model.mach
        #ncases = geom_model.ncases

        ID = 1
        is_beta0 = min(geom_model.betas) == max(geom_model.betas) and geom_model.betas[0] == 0.
        for icp in range(ncp):
            alpha = geom_model.alphas[icp]
            beta = geom_model.betas[icp]
            case_name = 'alpha=%s' % alpha
            if not is_beta0:
                case_name += ' beta=%s' %  beta

            imin = 0
            Cp_array = zeros(nnodes, dtype='float32')
            for unused_ipatch, Cp in sorted(model.pressures.items()):
                Cpv = Cp[icp].ravel()
                nCp = len(Cpv)
                try:
                    Cp_array[imin:imin + nCp] = Cpv
                except ValueError:
                    # agps stores implicit and explicit wakes
                    # we're skipping all wakes
                    pass
                imin += nCp

                #Cp_array2 = (Cp_array[self.elements[:, 0]] +
                             #Cp_array[self.elements[:, 1]] +
                             #Cp_array[self.elements[:, 2]] +
                             #Cp_array[self.elements[:, 3]]) / 4.
            results_form += [
                ('Cp Nodal - %s' % case_name, icase, []),
                #('Cp Centroidal', icase + 1, [],),
            ]
            Cpn_res = GuiResult(ID, 'Cp Nodal - %s' % case_name, 'Cp', 'node', Cp_array,
                                data_format='%.3f', colormap=colormap, uname='Cp_nodal')
            #Cpc_res = GuiResult(ID, 'Cp Centroidal', 'Cp', 'centroid', Cp_array2,
                                #data_format='%.3f', uname='Cp_centroidal')
            cases[icase] = (Cpn_res, (0, 'Cp_nodal'))
            icase += 1


        form = self.gui.get_form()
        assert len(results_form) > 0, results_form
        form.append(('Results: Mach=%s' % mach, None, results_form))

        dirname = os.path.dirname(panair_filename)
        panair_out_filename = os.path.join(dirname, 'panair.out')
        ft13_filename = os.path.join(dirname, 'ft13')

        if os.path.exists(panair_out_filename):
            out = read_panair_out(panair_out_filename, log=self.gui.log)
            unused_alphas = geom_model.alphas
            unused_betas = geom_model.betas
            icase, out_form = add_networks(out.networks, out.headers, is_beta0,
                                           ID, icase, cases, geom_model, nelements,
                                           self.gui.log, colormap=colormap)
            assert len(out_form) > 0, out_form
            form.append(('Out: Mach=%s' % mach, None, out_form))

            if os.path.exists(ft13_filename):
                out.read_ft13(ft13_filename)
                icase, ft13_form = add_networks(out.networks_ft13, out.headers_ft13, is_beta0,
                                                ID, icase, cases, geom_model, nelements,
                                                self.gui.log, colormap=colormap)
                assert len(ft13_form) > 0, ft13_form
                form.append(('Ft13: Mach=%s' % mach, None, ft13_form))


        #cases[icase + 1] = (Cpc_res, (0, 'Cp_centroidal'))
        self.gui._finish_results_io2(model_name, form, cases)

def add_networks(out_networks, out_headers, is_beta0,
                 ID, icase, cases, geom_model, nelements, log, colormap='jet'):
    out_form = []
    unused_nsolutions = len(out_networks)
    for isolution, networks in sorted(out_networks.items()):
        if networks == {}:
            log.info('skipping isolution=%s' % (isolution))
            continue
        print('isolution = ', isolution, geom_model.alphas)
        alpha = geom_model.alphas[isolution-1]
        beta = geom_model.betas[isolution-1]
        case_name = 'alpha=%s' % alpha
        if not is_beta0:
            case_name += ' beta=%s' %  beta

        out_form2 = []
        #print('----------')
        for iheader, title in enumerate(out_headers):
            #header = '%s - %s' % (title, case_name)
            header = title
            if title in ['x', 'y', 'z']:
                if isolution > 1:
                    continue
                header = title
                #break

            #print('%i %r' % (iheader, title))
            icell = 0
            icell0 = 0
            data_array = np.full(nelements, np.nan, dtype='float32')
            #print('----------')
            #print('***keys =', networks.keys())
            for inetwork, network in sorted(networks.items()):
                patch = geom_model.patches[inetwork-1]

                if patch.is_wake():  # wakes don't have results
                    continue

                datai = network.data[:, iheader]
                nrows = patch.nrows
                ncols = patch.ncols

                #print(patch.get_header())
                #print(inetwork, patch.inetwork, patch.network_name)
                #print(len(datai), datai.shape, nrows, ncols, nrows*ncols, (nrows-1)*(ncols-1))

                datai2 = datai.reshape(ncols-1, nrows-1)
                npatch_cells = len(datai)
                icell += npatch_cells
                data_array[icell0:icell] = datai2.ravel()
                icell0 += npatch_cells

            if header in ['x', 'y', 'z']:
                out_form += [
                    (header, icase, []),
                    #('Cp Centroidal', icase + 1, [],),
                ]
            else:
                out_form2 += [
                    (header, icase, []),
                    #('Cp Centroidal', icase + 1, [],),
                ]

            Cpn_res = GuiResult(ID, header, title, 'centroid', data_array,
                                data_format='%.3f', colormap=colormap, uname=header)
            #Cpc_res = GuiResult(ID, 'Cp Centroidal', 'Cp', 'centroid', Cp_array2,
                                #data_format='%.3f', uname='Cp_centroidal')
            cases[icase] = (Cpn_res, (0, header))
            icase += 1
        out_formi = (case_name, None, out_form2)
        out_form.append(out_formi)
    return icase, out_form
