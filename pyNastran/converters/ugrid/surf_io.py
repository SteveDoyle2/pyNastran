from numpy import vstack, amax, amin, arange, ones, zeros, where
from pyNastran.converters.ugrid.surf_reader import SurfReader, TagReader

#VTK_TRIANGLE = 5
from six import iteritems
from six.moves import range
import os
from numpy import unique
#from numpy import zeros, arange, mean, amax, amin, array

import vtk
from vtk import vtkTriangle, vtkQuad
from pyNastran.utils import print_bad_path


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
        self.grid2.Allocate(1, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.update_axes_length(dim_max)
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
                self.create_alternate_vtk_grid('failed_nodes', color=yellow, line_width=3, opacity=1.0)

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
                self.alt_grids['failed_nodes'].InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            self.alt_grids['failed_nodes'].SetPoints(points2)
            self._add_alt_actors(self.alt_grids)

            actor = self.geometry_actors['failed_nodes']
            actor.Modified()
            prop = actor.GetProperty()
            prop.SetRepresentationToPoints()
            prop.SetPointSize(10)


            # self.

        self.nElements = nelements
        self.grid.SetPoints(points)
        self.grid.Modified()
        #self.grid2.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            #self.grid2.Update()
        #print("updated grid")

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
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
        model = Cart3DReader(log=self.log, debug=False)
        #self.model_type = model.model_type
        #(nodes, elements, regions, loads) = model.read_cart3d(cart3dFileName)

        model.infilename = cart3d_filename
        if is_binary(infilename):
            model.infile = open(cart3d_filename, 'rb')
            (model.nPoints, model.nElements) = self.read_header_binary()
            points = model.read_points_binary(self.nPoints)
            elements = model.read_elements_binary(self.nElements)
            regions = model.read_regions_binary(self.nElements)
            #loads = {}
        else:
            model.infile = open(cart3d_filename, 'r')
            model.read_header_ascii()
            points = model.read_points_ascii(bypass=True)
            elements = model.read_elements_ascii(bypass=True)
            regions = model.read_regions_ascii(bypass=True)
            loads = model.read_results_ascii(0, model.infile, result_names=result_names)
        self.load_cart3d_geometry(cart3d_filename, dirname)


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
        print(unique(grid_bcs))

        normals = model.get_normals()
        cases[(ID, 0, 'ElementID', 1, 'centroid', '%i', '')] = eids
        cases[(ID, 1, 'NodeID',    1, 'node', '%i', '')] = nids
        cases[(ID, 2, 'SurfaceID', 1, 'centroid', '%i', '')] = surf_ids
        cases[(ID, 3, 'ReconFlag', 1, 'centroid', '%i', '')] = recon_flags
        cases[(ID, 4, 'GridBC',    1, 'centroid', '%i', '')] = grid_bcs
        cases[(ID, 5, 'NormalX',   1, 'centroid', '%.3f', '')] = normals[:, 0]
        cases[(ID, 6, 'NormalY',   1, 'centroid', '%.3f', '')] = normals[:, 1]
        cases[(ID, 7, 'NormalZ',   1, 'centroid', '%.3f', '')] = normals[:, 2]
        cases[(ID, 8, 'normSpacing', 1, 'node', '%.3e', '')] = norm_spacing
        cases[(ID, 9, 'BL_thick',    1, 'node', '%.3e', '')] = bl_thickness

        if os.path.exists(tag_filename):
            tagger = TagReader()
            data = tagger.read_tag_filename(tag_filename)

            int_data = ones((nelements, 8), dtype='int32') * -10.
            float_data = zeros((nelements, 2), dtype='float64')
            for key, datai in sorted(iteritems(data)):
                #self.log.info(datai)
                [name, is_visc, is_recon, is_rebuild, is_fixed, is_source, is_trans, is_delete, bl_spacing, bl_thickness, nlayers] = datai
                i = where(surf_ids == key)[0]
                int_data[i, :] = [is_visc, is_recon, is_rebuild, is_fixed, is_source, is_trans, is_delete, nlayers]
                float_data[i, :] = [bl_spacing, bl_thickness]
                self.log.info('data[%i] = %s' % (key, name))

            has_tag_data = True
            tag_form = []
            i = 10
            tag_form.append( ('is_visc',    i, []) )
            tag_form.append( ('is_recon',   i + 1, []) )
            tag_form.append( ('is_rebuild', i + 2, []) )
            tag_form.append( ('is_fixed',   i + 3, []) )
            tag_form.append( ('is_source',  i + 4, []) )
            tag_form.append( ('is_trans',   i + 5, []) )
            tag_form.append( ('is_delete',  i + 6, []) )
            tag_form.append( ('nlayers',    i + 7, []) )
            tag_form.append( ('bl_spacing',  i + 8, []) )
            tag_form.append( ('bl_thickness', i + 9, []) )

            cases[(ID, i, 'is_visc',    1, 'centroid', '%i', '')] = int_data[:, 0]
            cases[(ID, i + 1, 'is_recon',   1, 'centroid', '%i', '')] = int_data[:, 1]
            cases[(ID, i + 2, 'is_rebuild', 1, 'centroid', '%i', '')] = int_data[:, 2]
            cases[(ID, i + 3, 'is_fixed',   1, 'centroid', '%i', '')] = int_data[:, 3]
            cases[(ID, i + 4, 'is_source',  1, 'centroid', '%i', '')] = int_data[:, 4]
            cases[(ID, i + 5, 'is_trans',   1, 'centroid', '%i', '')] = int_data[:, 5]
            cases[(ID, i + 6, 'is_delete', 1, 'centroid', '%i', '')] = int_data[:, 6]
            cases[(ID, i + 7, 'nlayers',   1, 'centroid', '%i', '')] = int_data[:, 7]

            cases[(ID, i + 8, 'bl_spacing', 1, 'centroid', '%.3e', '')] = float_data[:, 0]
            cases[(ID, i + 9, 'bl_thickness', 1, 'centroid', '%.3e', '')] = float_data[:, 1]

        form = [
            ('Geometry', None, geometry_form),
        ]
        if has_tag_data:
            form.append(('Tag Data', None, tag_form),)

        results_form = []
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
