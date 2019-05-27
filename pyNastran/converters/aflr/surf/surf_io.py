import os
from collections import OrderedDict

from numpy import vstack, amax, amin, arange, ones, zeros, where

#VTK_TRIANGLE = 5
import vtk

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.aflr.surf.surf_reader import SurfReader, TagReader
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_types, numpy_to_vtk_points)
from pyNastran.gui.qt_files.colors import YELLOW_FLOAT


class SurfIO:
    def __init__(self, gui):
        self.gui = gui

    def get_surf_wildcard_geometry_results_functions(self):
        data = (
            'AFLR3 Surf',
            'AFLR3 Surf (*.surf)', self.load_surf_geometry,
            None, None)
        return data

    def load_surf_geometry(self, surf_filename, name=None, plot=True):
        model_name = name
        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skip_reading:
        #    return

        model = SurfReader()

        self.gui.model_type = 'surf'
        self.gui.log.debug('surf_filename = %s' % surf_filename)

        model.read_surf(surf_filename)
        nnodes = model.nodes.shape[0]
        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        nelements = ntris + nquads

        nodes = model.nodes
        self.gui.nelements = nelements
        self.gui.nnodes = nnodes

        #print("nNodes = %s" % self.nnodes)
        #print("nElements = %s" % self.nelements)
        assert nelements > 0, nelements


        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)
        self.gui.log.info('max = %s' % mmax)
        self.gui.log.info('min = %s' % mmin)

        points = numpy_to_vtk_points(nodes)
        tris = model.tris - 1
        quads = model.quads - 1

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        elements = []
        etypes = []
        if ntris:
            elements.append(tris)
            etypes.append(5) # vtkTriangle().GetCellType()
        if nquads:
            elements.append(quads)
            etypes.append(9) # vtkQuad().GetCellType()
        create_vtk_cells_of_constant_element_types(grid, elements, etypes)


        model.read_surf_failnode(surf_filename)
        if len(model.nodes_failed):
            if 'failed_nodes' not in self.gui.alt_grids:
                self.gui.create_alternate_vtk_grid('failed_nodes', color=YELLOW_FLOAT,
                                                   line_width=3, opacity=1.0)

            ifailed = where(model.nodes_failed == 1)[0]
            nfailed = len(ifailed)
            failed_grid = self.gui.alt_grids['failed_nodes']
            failed_grid.Allocate(nfailed, 1000)
            points2 = vtk.vtkPoints()
            points2.SetNumberOfPoints(nfailed)

            for j, nid in enumerate(model.nodes_failed):
                elem = vtk.vtkVertex()
                c = nodes[nid - 1, :]
                self.gui.log.debug('%s %s' % (nid, c))
                points2.InsertPoint(j, *c)
                elem.GetPointIds().SetId(0, j)
                failed_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            failed_grid.SetPoints(points2)
            self.gui._add_alt_actors(self.gui.alt_grids)

            actor = self.gui.geometry_actors['failed_nodes']
            actor.Modified()
            prop = actor.GetProperty()
            prop.SetRepresentationToPoints()
            prop.SetPointSize(10)

        self.gui.nelements = nelements
        grid.SetPoints(points)
        grid.Modified()

        # loadSurfResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['AFLR Surface', '']}
        cases = OrderedDict()
        ID = 1

        form, cases, node_ids, element_ids = self._fill_surf_case(
            surf_filename, cases, ID, nnodes, nelements, model)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        if plot:
            self.gui._finish_results_io2(model_name, form, cases)

    def clear_surf(self):
        pass

    #def _load_surf_results(self, openfoam_filename):
        #raise NotImplementedError()

    def _fill_surf_case(self, surf_filename, cases, unused_ID, nnodes, nelements, model):
        """builds the results for the *.surf AFLR3 input file"""
        base, ext = os.path.splitext(surf_filename)
        assert ext == '.surf', surf_filename

        tag_filename = base + '.tags'

        unused_cases_new = []
        has_tag_data = False
        results_form = []
        geometry_form = [
            #('Region', 0, []),
            ('ElementID', 0, []),
            ('NodeID', 1, []),
            ('SurfaceID', 2, []),
            ('ReconFlag', 3, []),
            ('GridBC', 4, []),

            ('NormalX', 5, []),
            ('NormalY', 6, []),
            ('NormalZ', 7, []),

            ('normSpacing', 8, []),
            ('BL_thick', 9, []),
        ]
        nids = arange(1, nnodes + 1)
        norm_spacing = model.node_props[:, 0]
        bl_thickness = model.node_props[:, 1]

        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        #nelements = ntris + nquads
        eids = arange(1, nelements + 1)

        if ntris and nquads:
            element_props = vstack([model.tri_props, model.quad_props])
        elif ntris:
            element_props = model.tri_props
        elif nquads:
            element_props = model.quad_props

        surf_ids = element_props[:, 0]
        recon_flags = element_props[:, 1]
        grid_bcs = element_props[:, 2]
        #print(unique(grid_bcs))

        normals = model.get_normals()
        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        surface_res = GuiResult(0, header='SurfaceID', title='SurfaceID',
                                location='centroid', scalar=surf_ids)

        recon_res = GuiResult(0, header='ReconFlag', title='ReconFlag',
                              location='centroid', scalar=recon_flags)
        gridbc_res = GuiResult(0, header='GridBC', title='GridBC',
                               location='centroid', scalar=grid_bcs)
        normalx_res = GuiResult(0, header='NormalX', title='NormalX',
                                location='centroid', scalar=normals[:, 0])
        normaly_res = GuiResult(0, header='NormalY', title='NormalY',
                                location='centroid', scalar=normals[:, 1])
        normalz_res = GuiResult(0, header='NormalZ', title='NormalZ',
                                location='centroid', scalar=normals[:, 2])

        normspacing_res = GuiResult(0, header='NormSpacing', title='NormSpacing',
                                    location='node', scalar=norm_spacing)
        blthick_res = GuiResult(0, header='BL_thick', title='BL_thick',
                                location='node', scalar=bl_thickness)

        icase = 0
        cases[icase] = (eid_res, (0, 'ElementID'))
        cases[icase + 1] = (nid_res, (0, 'NodeID'))
        cases[icase + 2] = (surface_res, (0, 'SurfaceID'))

        cases[icase + 3] = (recon_res, (0, 'ReconFlag'))
        cases[icase + 4] = (gridbc_res, (0, 'GridBC'))
        cases[icase + 5] = (normalx_res, (0, 'NormalX'))
        cases[icase + 6] = (normaly_res, (0, 'NormalY'))
        cases[icase + 7] = (normalz_res, (0, 'NormalZ'))

        cases[icase + 8] = (normspacing_res, (0, 'NormSpacing'))
        cases[icase + 9] = (blthick_res, (0, 'BL_thick'))
        icase += 10

        if os.path.exists(tag_filename):
            tagger = TagReader()
            data = tagger.read_tag_filename(tag_filename)

            int_data = ones((nelements, 8), dtype='int32') * -10.
            float_data = zeros((nelements, 2), dtype='float64')
            for key, datai in sorted(data.items()):
                #self.log.info(datai)
                [name, is_visc, is_recon, is_rebuild, is_fixed, is_source,
                 is_trans, is_delete, bl_spacing, bl_thickness, nlayers] = datai
                i = where(surf_ids == key)[0]
                int_data[i, :] = [is_visc, is_recon, is_rebuild, is_fixed,
                                  is_source, is_trans, is_delete, nlayers]
                float_data[i, :] = [bl_spacing, bl_thickness]
                self.gui.log.info('data[%i] = %s' % (key, name))

            has_tag_data = True
            tag_form = []
            tag_form.append(('is_visc', icase, []))
            tag_form.append(('is_recon', icase + 1, []))
            tag_form.append(('is_rebuild', icase + 2, []))
            tag_form.append(('is_fixed', icase + 3, []))
            tag_form.append(('is_source', icase + 4, []))
            tag_form.append(('is_trans', icase + 5, []))
            tag_form.append(('is_delete', icase + 6, []))
            tag_form.append(('nlayers', icase + 7, []))
            tag_form.append(('bl_spacing', icase + 8, []))
            tag_form.append(('bl_thickness', icase + 9, []))

            visc_res = GuiResult(0, header='is_visc', title='is_visc',
                                 location='centroid', scalar=int_data[:, 0])
            recon_res = GuiResult(0, header='is_recon', title='is_recon',
                                  location='centroid', scalar=int_data[:, 1])
            rebuild_res = GuiResult(0, header='is_rebuild', title='is_rebuild',
                                    location='centroid', scalar=int_data[:, 2])
            fixed_res = GuiResult(0, header='is_fixed', title='is_fixed',
                                  location='centroid', scalar=int_data[:, 3])
            source_res = GuiResult(0, header='is_source', title='is_source',
                                   location='centroid', scalar=int_data[:, 4])
            trans_res = GuiResult(0, header='is_trans', title='is_trans',
                                  location='centroid', scalar=int_data[:, 5])
            delete_res = GuiResult(0, header='is_delete', title='is_delete',
                                   location='centroid', scalar=int_data[:, 6])
            nlayers_res = GuiResult(0, header='nlayers', title='nlayers',
                                    location='centroid', scalar=int_data[:, 7])

            spacing_res = GuiResult(0, header='bl_spacing', title='bl_spacing',
                                    location='centroid', scalar=float_data[:, 0])
            blthickness_res = GuiResult(0, header='bl_thickness', title='bl_thickness',
                                        location='centroid', scalar=float_data[:, 1])

            cases[icase] = (visc_res, (0, 'is_visc'))
            cases[icase + 1] = (recon_res, (0, 'is_recon'))
            cases[icase + 2] = (rebuild_res, (0, 'is_rebuild'))
            cases[icase + 3] = (fixed_res, (0, 'is_fixed'))
            cases[icase + 4] = (source_res, (0, 'is_source'))
            cases[icase + 5] = (trans_res, (0, 'is_trans'))
            cases[icase + 6] = (delete_res, (0, 'is_delete'))
            cases[icase + 7] = (nlayers_res, (0, 'nlayers'))
            cases[icase + 8] = (spacing_res, (0, 'bl_spacing'))
            cases[icase + 9] = (blthickness_res, (0, 'bl_thickness'))

        form = [
            ('Geometry', None, geometry_form),
        ]
        if has_tag_data:
            form.append(('Tag Data', None, tag_form),)

        results_form = []
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases, nids, eids
