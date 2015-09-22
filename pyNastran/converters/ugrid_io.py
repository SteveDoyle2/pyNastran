from numpy import vstack, amax, amin, arange, ones, zeros, where
from ulst.formats.aflr3.surf_reader import TagReader
from ulst.formats.aflr3.ugrid_reader import UGRID
from ulst.formats.aflr2.ugrid2d_reader import UGRID2D_Reader

#VTK_TRIANGLE = 5
from six import iteritems
from six.moves import range
import os
from numpy import zeros, unique, where
#from numpy import zeros, arange, mean, amax, amin, array, where

import vtk
from vtk import vtkTriangle, vtkQuad
from pyNastran.utils import print_bad_path, is_binary_file


class UGRID_IO(object):
    def __init__(self):
        pass

    def get_ugrid_wildcard_geometry_results_functions(self):
        data = (
            'AFLR3 Ugrid',
            'AFLR3 Ugrid (*.ugrid)', self.load_ugrid_geometry,
             None, None)
        return data

    def load_ugrid_geometry(self, ugrid_filename, dirname, plot=True):
        #skipReading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skipReading:
        #    return
        if is_binary_file(ugrid_filename):
            model = UGRID(log=self.log, debug=True)
            base, fmt, ext = os.path.basename(ugrid_filename).split('.')
            is_2d = False
        else:
            base, ext = os.path.basename(ugrid_filename).split('.')
            model = UGRID2D_Reader(log=self.log, debug=True)
            is_2d = True

        self.modelType = 'ugrid'
        print('ugrid_filename = %s' % ugrid_filename)


        assert ext == 'ugrid', ugrid_filename
        model.read_ugrid(ugrid_filename)

        if is_2d:
            tris = model.tris
            quads = model.quads
        else:
            tris = model.tris - 1
            quads = model.quads - 1

        #self.nodes = nodes
        #self.tris  = tris
        #self.quads = quads
        #self.pids = pids

        #self.tets = tets
        #self.penta5s = penta5s
        #self.penta6s = penta6s
        #self.hexas = hexas

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

        self.nElements = nelements
        self.grid.SetPoints(points)
        self.grid.Modified()
        #self.grid2.Modified()
        print('update...')
        if hasattr(self.grid, 'Update'):
            self.grid.Update()
            #self.grid2.Update()
        #print("updated grid")

        # loadCart3dResults - regions/loads
        self.TurnTextOn()
        self.scalarBar.VisibilityOn()
        self.scalarBar.Modified()

        self.iSubcaseNameMap = {1: ['AFLR UGRID Surface', '']}
        cases = {}
        ID = 1

        if hasattr(model, 'pids'):
            form, cases = self._fill_ugrid3d_case(ugrid_filename, cases, ID, nnodes, nelements, model)
        else:
            form, cases = self._fill_ugrid2d_case(ugrid_filename, cases, ID, nnodes, nelements, model)

        if plot:
            self._finish_results_io2(form, cases)

    def clear_surf(self):
        pass

    #def _load_ugrid_results(self, openfoam_filename, dirname):
        #model = Cart3DReader(log=self.log, debug=False)
        ##self.modelType = model.modelType
        ##(nodes, elements, regions, loads) = model.read_cart3d(cart3dFileName)

        #model.infilename = cart3d_filename
        #if is_binary(infilename):
            #model.infile = open(cart3d_filename, 'rb')
            #(model.nPoints, model.nElements) = self.read_header_binary()
            #points = model.read_points_binary(self.nPoints)
            #elements = model.read_elements_binary(self.nElements)
            #regions = model.read_regions_binary(self.nElements)
            ##loads = {}
        #else:
            #model.infile = open(cart3d_filename, 'r')
            #model.read_header_ascii()
            #points = model.read_points_ascii(bypass=True)
            #elements = model.read_elements_ascii(bypass=True)
            #regions = model.read_regions_ascii(bypass=True)
            #loads = model.read_results_ascii(0, model.infile, result_names=result_names)
        #self.load_cart3d_geometry(cart3d_filename, dirname)


    def _fill_ugrid2d_case(self, base, cases, ID, nnodes, nelements, model):
        cases_new = []
        results_form = []

        geometry_form = [
            ('ElementID', 0, []),
            ('NodeID', 1, []),
        ]

        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        cases[(ID, 0, 'ElementID', 1, 'centroid', '%i')] = eids
        cases[(ID, 1, 'NodeID', 1, 'node', '%i')] = nids

        form = [
            ('Geometry', None, geometry_form),
        ]
        #if has_tag_data:
            #form.append(('Tag Data', None, tag_form),)

        #results_form = []
        #if len(results_form):
            #form.append(('Results', None, results_form))
        return form, cases

    def _fill_ugrid3d_case(self, base, cases, ID, nnodes, nelements, model):
        tag_filename = base + '.tags'
        mapbc_filename = base.split('.')[0] + '.mapbc'
        print('mapbc_filename = %r' % mapbc_filename)

        cases_new = []
        has_tag_data = False
        has_mapbc_data = False
        results_form = []
        mapbc_form = []

        geometry_form = [
            #('Region', 0, []),
            ('ElementID', 0, []),
            ('NodeID', 1, []),
            ('SurfaceID', 2, []),
            #('normSpacing', 3, []),
            #('BL_thick', 4, []),
            #('ReconFlag', 5, []),
            #('GridBC', 6, []),
        ]

        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        #nelements = ntris + nquads
        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        #grid_bcs = element_props[:, 2]

        #npids = len(model.pids)
        pids = model.pids
        cases[(ID, 0, 'ElementID', 1, 'centroid', '%i')] = eids
        cases[(ID, 1, 'NodeID', 1, 'node', '%i')] = nids
        cases[(ID, 2, 'SurfaceID', 1, 'centroid', '%i')] = pids

        n = 3
        if os.path.exists(tag_filename):
            #surf_ids = element_props[:, 0]
            #recon_flags = element_props[:, 1]
            #cases[(ID, 2, 'ReconFlag', 1, 'centroid', '%i')] = recon_flags
            #cases[(ID, 3, 'GridBC',    1, 'centroid', '%i')] = grid_bcs

            tagger = TagReader()
            data = tagger.read_tag_filename(tag_filename)

            int_data = ones((nelements, 8), dtype='int32') * -10.
            float_data = zeros((nelements, 2), dtype='float64')
            for key, datai in sorted(iteritems(data)):
                #self.log.info(datai)
                [name, is_visc, is_recon, is_rebuild, is_fixed, is_source, is_trans, is_delete, bl_spacing, bl_thickness, nlayers] = datai
                i = where(pids == key)[0]
                int_data[i, :] = [is_visc, is_recon, is_rebuild, is_fixed, is_source, is_trans, is_delete, nlayers]
                float_data[i, :] = [bl_spacing, bl_thickness]
                self.log.info('data[%i] = %s' % (key, name))

            has_tag_data = True
            tag_form = []
            tag_form.append( ('is_visc',      n, []) )
            tag_form.append( ('is_recon',     n+1, []) )
            tag_form.append( ('is_rebuild',   n+2, []) )
            tag_form.append( ('is_fixed',     n+3, []) )
            tag_form.append( ('is_source',    n+4, []) )
            tag_form.append( ('is_trans',     n+5, []) )
            tag_form.append( ('is_delete',    n+6, []) )
            tag_form.append( ('nlayers',      n+7, []) )
            tag_form.append( ('bl_spacing',   n+8, []) )
            tag_form.append( ('bl_thickness', n+9, []) )

            cases[(ID, n, 'is_visc',      1, 'centroid', '%i')] = int_data[:, 0]
            cases[(ID, n + 1, 'is_recon',   1, 'centroid', '%i')] = int_data[:, 1]
            cases[(ID, n + 2, 'is_rebuild', 1, 'centroid', '%i')] = int_data[:, 2]
            cases[(ID, n + 3, 'is_fixed',   1, 'centroid', '%i')] = int_data[:, 3]
            cases[(ID, n + 4, 'is_source',  1, 'centroid', '%i')] = int_data[:, 4]
            cases[(ID, n + 5, 'is_trans',   1, 'centroid', '%i')] = int_data[:, 5]
            cases[(ID, n + 6, 'is_delete',  1, 'centroid', '%i')] = int_data[:, 6]
            cases[(ID, n + 7, 'nlayers',    1, 'centroid', '%i')] = int_data[:, 7]

            cases[(ID, n + 8, 'bl_spacing',   1, 'centroid', '%.3e')] = float_data[:, 0]
            cases[(ID, n + 9, 'bl_thickness', 1, 'centroid', '%.3e')] = float_data[:, 1]
            n += 10
        else:
            self.log_info('tag_filename=%r could not be found' % tag_filename)

        if os.path.exists(mapbc_filename):
            has_mapbc_data = True
            mapbc = open(mapbc_filename, 'r')
            lines = mapbc.readlines()
            lines = [line.strip() for line in lines
                     if not line.strip().startswith('#') and line.strip()]
            npatches = int(lines[0])
            mapbcs = zeros(pids.shape, dtype='int32')
            for ipatch in range(npatches):
                line = lines[ipatch + 1]
                iline, bc_num, name = line.split()
                iline = int(iline)
                bc_num = int(bc_num)
                assert ipatch + 1 == iline, 'line=%r; ipatch=%s iline=%s' % (line, ipatch + 1, iline)
                islot = where(pids == ipatch + 1)[0]
                if len(islot) == 0:
                    upids = unique(pids)
                    msg = 'ipatch=%s not found in pids=%s' % (ipatch + 1, upids)
                    raise RuntimeError(msg)
                mapbcs[islot] = bc_num
                print(line)
            mapbc_form.append(('Map BC', n, []))
            cases[(ID, n, 'Map BC', 1, 'centroid', '%i')] = mapbcs
        else:
            self.log_info('mapbc_filename=%r could not be found' % mapbc_filename)


        #norm_spacing = model.node_props[:, 0]
        #bl_thickness = model.node_props[:, 1]
        #cases[(ID, 1, 'normSpacing', 1, 'node', '%.3e')] = norm_spacing
        #cases[(ID, 2, 'BL_thick',    1, 'node', '%.3e')] = bl_thickness

        form = [
            ('Geometry', None, geometry_form),
        ]
        if has_tag_data:
            form.append(('Tag Data', None, tag_form),)
        if has_mapbc_data:
            form.append(('Map BC Data', None, mapbc_form),)

        results_form = []
        if len(results_form):
            form.append(('Results', None, results_form))
        print(form)
        return form, cases
