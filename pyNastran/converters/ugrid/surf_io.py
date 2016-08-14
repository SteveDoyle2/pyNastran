from __future__ import print_function
import os

from six import iteritems
from six.moves import range
#from numpy import unique
from numpy import vstack, amax, amin, arange, ones, zeros, where
#from numpy import zeros, arange, mean, amax, amin, array

#VTK_TRIANGLE = 5
import vtk
from vtk import vtkTriangle, vtkQuad

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.ugrid.surf_reader import SurfReader, TagReader


class SurfIO(object):
    def __init__(self):
        pass

    def get_surf_wildcard_geometry_results_functions(self):
        data = (
            'AFLR3 Surf',
            'AFLR3 Surf (*.surf)', self.load_surf_geometry,
            None, None)
        return data

    def load_surf_geometry(self, surf_filename, dirname, plot=True):
        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skip_reading:
        #    return

        model = SurfReader()

        self.model_type = 'surf'
        print('surf_filename = %s' % surf_filename)

        model.read_surf(surf_filename)
        nnodes = model.nodes.shape[0]
        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        nelements = ntris + nquads

        nodes = model.nodes
        self.nElements = nelements
        self.nNodes = nnodes

        print("nNodes = %s" % self.nNodes)
        print("nElements = %s" % self.nElements)
        assert nelements > 0, nelements

        self.grid.Allocate(self.nElements, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.create_global_axes(dim_max)
        self.log.info('max = %s' % mmax)
        self.log.info('min = %s' % mmin)

        for inode, node in enumerate(nodes):
            points.InsertPoint(inode, node)

        tris = model.tris - 1
        quads = model.quads - 1


        if ntris:
            for eid, element in enumerate(tris):
                elem = vtkTriangle()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())
        if nquads:
            for eid, element in enumerate(quads):
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, element[0])
                elem.GetPointIds().SetId(1, element[1])
                elem.GetPointIds().SetId(2, element[2])
                elem.GetPointIds().SetId(3, element[3])
                self.grid.InsertNextCell(elem.GetCellType(),
                                         elem.GetPointIds())

        model.read_surf_failnode(surf_filename)
        if len(model.nodes_failed):
            if 'failed_nodes' not in self.alt_grids:
                yellow = (1., 1., 0.)
                self.create_alternate_vtk_grid('failed_nodes', color=yellow,
                                               line_width=3, opacity=1.0)

            ifailed = where(model.nodes_failed == 1)[0]
            nfailed = len(ifailed)
            self.alt_grids['failed_nodes'].Allocate(nfailed, 1000)
            grid2 = self.alt_grids['failed_nodes']
            points2 = vtk.vtkPoints()
            points2.SetNumberOfPoints(nfailed)

            for j, nid in enumerate(model.nodes_failed):
                elem = vtk.vtkVertex()
                c = nodes[nid - 1, :]
                print(nid, c)
                points2.InsertPoint(j, *c)
                elem.GetPointIds().SetId(0, j)
                self.alt_grids['failed_nodes'].InsertNextCell(elem.GetCellType(),
                                                              elem.GetPointIds())
            self.alt_grids['failed_nodes'].SetPoints(points2)
            self._add_alt_actors(self.alt_grids)

            actor = self.geometry_actors['failed_nodes']
            actor.Modified()
            prop = actor.GetProperty()
            prop.SetRepresentationToPoints()
            prop.SetPointSize(10)

        self.nElements = nelements
        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
        #print("updated grid")

        # loadCart3dResults - regions/loads
        self. turn_text_on()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['AFLR Surface', '']}
        cases = {}
        ID = 1

        form, cases = self._fill_surf_case(surf_filename, cases, ID, nnodes, nelements, model)
        if plot:
            self._finish_results_io2(form, cases)

    def clear_surf(self):
        pass

    def _load_surf_results(self, openfoam_filename, dirname):
        raise NotImplementedError()

    def _fill_surf_case(self, surf_filename, cases, ID, nnodes, nelements, model):
        base, ext = os.path.splitext(surf_filename)
        assert ext == '.surf', surf_filename

        tag_filename = base + '.tags'

        cases_new = []
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

        icase = 0
        cases[icase] = (eid_res, (0, 'ElementID'))
        cases[icase + 1] = (nid_res, (0, 'NodeID'))
        cases[icase + 2] = (surface_res, (0, 'SurfaceID'))

        cases[(ID, 3, 'ReconFlag', 1, 'centroid', '%i', '')] = recon_flags
        cases[(ID, 4, 'GridBC', 1, 'centroid', '%i', '')] = grid_bcs
        cases[(ID, 5, 'NormalX', 1, 'centroid', '%.3f', '')] = normals[:, 0]
        cases[(ID, 6, 'NormalY', 1, 'centroid', '%.3f', '')] = normals[:, 1]
        cases[(ID, 7, 'NormalZ', 1, 'centroid', '%.3f', '')] = normals[:, 2]
        cases[(ID, 8, 'normSpacing', 1, 'node', '%.3e', '')] = norm_spacing
        cases[(ID, 9, 'BL_thick', 1, 'node', '%.3e', '')] = bl_thickness

        if os.path.exists(tag_filename):
            tagger = TagReader()
            data = tagger.read_tag_filename(tag_filename)

            int_data = ones((nelements, 8), dtype='int32') * -10.
            float_data = zeros((nelements, 2), dtype='float64')
            for key, datai in sorted(iteritems(data)):
                #self.log.info(datai)
                [name, is_visc, is_recon, is_rebuild, is_fixed, is_source,
                 is_trans, is_delete, bl_spacing, bl_thickness, nlayers] = datai
                i = where(surf_ids == key)[0]
                int_data[i, :] = [is_visc, is_recon, is_rebuild, is_fixed,
                                  is_source, is_trans, is_delete, nlayers]
                float_data[i, :] = [bl_spacing, bl_thickness]
                self.log.info('data[%i] = %s' % (key, name))

            has_tag_data = True
            tag_form = []
            icase = 10
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
                                 location='node', scalar=int_data[:, 0])
            recon_res = GuiResult(0, header='is_recon', title='is_recon',
                                  location='node', scalar=int_data[:, 1])
            rebuild_res = GuiResult(0, header='is_rebuild', title='is_rebuild',
                                    location='node', scalar=int_data[:, 2])
            fixed_res = GuiResult(0, header='is_fixed', title='is_fixed',
                                  location='node', scalar=int_data[:, 3])
            source_res = GuiResult(0, header='is_source', title='is_source',
                                   location='node', scalar=int_data[:, 4])
            trans_res = GuiResult(0, header='is_trans', title='is_trans',
                                  location='node', scalar=int_data[:, 5])
            delete_res = GuiResult(0, header='is_delete', title='is_delete',
                                   location='node', scalar=int_data[:, 6])
            nlayers_res = GuiResult(0, header='nlayers', title='nlayers',
                                    location='node', scalar=int_data[:, 7])

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
        return form, cases
