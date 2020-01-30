import os
from collections import OrderedDict

from numpy import amax, amin, arange, ones, zeros, where, unique

#VTK_TRIANGLE = 5
import vtk
#from vtk import vtkTriangle, vtkQuad

from pyNastran.converters.aflr.surf.surf_reader import TagReader
from pyNastran.converters.aflr.ugrid.ugrid_reader import UGRID
from pyNastran.converters.aflr.ugrid.ugrid2d_reader import UGRID2D_Reader
from pyNastran.utils import is_binary_file
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_types, numpy_to_vtk_points)
from pyNastran.gui.qt_files.colors import RED_FLOAT


class UGRID_IO:
    def __init__(self, gui):
        self.gui = gui

    def get_ugrid_wildcard_geometry_results_functions(self):
        data = (
            'ugrid',
            'AFLR2/AFLR3 UGrid2D (*.ugrid)', self.load_ugrid_geometry,  # 2d
            None, None)
        return data

    def get_ugrid3d_wildcard_geometry_results_functions(self):
        data = (
            'ugrid3d',
            'AFLR3 Ugrid3D (*.ugrid)', self.load_ugrid3d_geometry,
            None, None)
        return data

    #def load_ugrid_geometry_2d(self, ugrid_filename, name='main', plot=True):
        #"""Loads a UGRID3D as a 2D file"""
        #self._load_ugrid_geometry(ugrid_filename, read_solids=False, name=name, plot=plot)

    def load_ugrid3d_geometry(self, ugrid_filename, name='main', plot=True):
        """Loads a UGRID3D as a 3D file"""
        self.load_ugrid_geometry(ugrid_filename, read_solids=True, name=name, plot=plot)


    def load_ugrid_geometry(self, ugrid_filename, read_solids=False, name='main', plot=True):
        """
        The entry point for UGRID geometry loading.

        Parameters
        ----------
        ugrid_filename : str
            the ugrid filename to load
        read_solids : bool
            True : load the tets/pentas/hexas from the UGRID3D model
            False : UGRID2D or limits the UGRID3D model to tris/quads
        name : str
            the name of the "main" actor for the GUI
        plot : bool; default=True
            should the model be generated or should we wait until
            after the results are loaded
        """
        model_name = name
        #skip_reading = self.remove_old_openfoam_geometry(openfoam_filename)
        #if skip_reading:
        #    return
        model, is_2d, is_3d = get_ugrid_model(ugrid_filename, read_solids, self.gui.log)
        self.gui.model_type = 'ugrid'
        self.gui.log.debug('ugrid_filename = %s' % ugrid_filename)

        model.read_ugrid(ugrid_filename)
        self.gui.model = model

        nnodes = model.nodes.shape[0]
        ntris = model.tris.shape[0]
        nquads = model.quads.shape[0]
        ntets = 0
        npenta5s = 0
        npenta6s = 0
        nhexas = 0
        if is_2d:
            tris = model.tris
            quads = model.quads
            nelements = ntris + nquads
        else:
            if read_solids:
                ntets = model.tets.shape[0]
                npenta5s = model.penta5s.shape[0]
                npenta6s = model.penta6s.shape[0]
                nhexas = model.hexas.shape[0]
                tets = model.tets - 1
                penta5s = model.penta5s - 1
                penta6s = model.penta6s - 1
                hexas = model.hexas - 1
                nelements = ntets + npenta5s + npenta6s + nhexas
            else:
                tris = model.tris - 1
                quads = model.quads - 1
                nelements = ntris + nquads

        #self.nodes = nodes
        #self.tris  = tris
        #self.quads = quads
        #self.pids = pids

        #self.tets = tets
        #self.penta5s = penta5s
        #self.penta6s = penta6s
        #self.hexas = hexas


        nodes = model.nodes
        self.gui.nelements = nelements
        self.gui.nnodes = nnodes

        self.gui.log.info("nnodes=%s nelements=%s" % (self.gui.nnodes, self.gui.nelements))
        assert nelements > 0, nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)
        self.gui.log.info('max = %s' % mmax)
        self.gui.log.info('min = %s' % mmin)

        if is_3d and read_solids:
            diff_node_ids = model.check_hanging_nodes(stop_on_diff=False)
            if len(diff_node_ids):
                self.gui.create_alternate_vtk_grid(
                    'hanging_nodes', color=RED_FLOAT, line_width=5, opacity=1.,
                    point_size=10, representation='point')
                self._add_ugrid_nodes_to_grid('hanging_nodes', diff_node_ids, nodes)
                self.gui._add_alt_actors(self.gui.alt_grids)

        points = numpy_to_vtk_points(nodes)

        elements = []
        etypes = []
        if is_2d or not read_solids:
            if ntris:
                elements.append(tris)
                etypes.append(5) # vtkTriangle().GetCellType()
            if nquads:
                elements.append(quads)
                etypes.append(9) # vtkQuad().GetCellType()
        elif is_3d:
            if ntets:
                elements.append(tets)
                etypes.append(10) # VTK_TETRA().GetCellType()
            if npenta5s:
                elements.append(penta5s)
                etypes.append(14) # vtk.vtkPyramid().GetCellType()
            if npenta6s:
                elements.append(penta6s)
                etypes.append(13) # VTK_WEDGE().GetCellType()
            if nhexas:
                elements.append(hexas)
                etypes.append(12) # VTK_HEXAHEDRON().GetCellType()

        self.gui.model.elements = elements
        self.gui.model.etypes = etypes
        create_vtk_cells_of_constant_element_types(grid, elements, etypes)

        self.gui.nelements = nelements
        grid.SetPoints(points)
        grid.Modified()

        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['AFLR UGRID Surface', '']}
        cases = OrderedDict()
        ID = 1

        if hasattr(model, 'pids'):
            form, cases, node_ids, element_ids = self._fill_ugrid3d_case(
                ugrid_filename, cases, ID, nnodes, nelements, model, read_solids)
        else:
            form, cases, node_ids, element_ids = self._fill_ugrid2d_case(
                cases, ID, nnodes, nelements)

        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        if plot:
            self.gui._finish_results_io2(model_name, form, cases)

    def _add_ugrid_nodes_to_grid(self, name, diff_node_ids, nodes):
        """
        based on:
          _add_nastran_nodes_to_grid
        """
        nnodes = nodes.shape[0]
        assert nnodes > 0, nnodes
        # if nnodes == 0:
            # return
        nnodes = len(diff_node_ids)
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nnodes)

        alt_grid = self.gui.alt_grids[name]
        for nid in diff_node_ids:
            node = nodes[nid, :]
            self.gui.log.info('nid=%s node=%s' % (nid, node))
            points.InsertPoint(nid, *node)

            #if 1:
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, nid)
            #else:
                #elem = vtk.vtkSphere()
                #sphere_size = self._get_sphere_size(dim_max)
                #elem.SetRadius(sphere_size)
                #elem.SetCenter(points.GetPoint(nid))

            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        alt_grid.SetPoints(points)

    def clear_surf(self):
        pass

    # def _load_ugrid_results(self, openfoam_filename, dirname):
        # pass

    def _fill_ugrid2d_case(self, cases, unused_id, nnodes, nelements):
        #cases_new = []
        #results_form = []
        colormap = self.gui.settings.colormap
        geometry_form = [
            ('ElementID', 0, []),
            ('NodeID', 1, []),
        ]

        #ntris = model.tris.shape[0]
        #nquads = model.quads.shape[0]
        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids, colormap=colormap)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nids, colormap=colormap)

        icase = 0
        cases[icase] = (eid_res, (0, 'ElementID'))
        cases[icase + 1] = (nid_res, (0, 'NodeID'))

        form = [
            ('Geometry', None, geometry_form),
        ]
        #if has_tag_data:
            #form.append(('Tag Data', None, tag_form),)

        #results_form = []
        #if len(results_form):
            #form.append(('Results', None, results_form))
        return form, cases, nids, eids

    def _fill_ugrid3d_case(self, base, cases, ID, nnodes, nelements, model, read_solids):
        if os.path.exists(base):
            # base = 'C:/data/'
            # tag_filename = 'C:/data/.tags'
            self.gui.log.info('mapbc_filename does not exist')
            self.gui.log.info('tag_filename does not exist')
            tag_filename = None
            mapbc_filename = None
        else:
            tag_filename = base + '.tags'
            mapbc_filename = base.split('.')[0] + '.mapbc'
            self.gui.log.info('mapbc_filename = %r' % mapbc_filename)

        colormap = self.gui.settings.colormap

        #cases_new = []
        has_tag_data = False
        has_mapbc_data = False
        results_form = []
        mapbc_form = []

        geometry_form = [
            #('Region', 0, []),
            ('ElementID', 0, []),
            ('NodeID', 1, []),
            ('Normals', 2, []),
            #('normSpacing', 3, []),
            #('BL_thick', 4, []),
            #('ReconFlag', 5, []),
            #('GridBC', 6, []),
        ]
        if not read_solids:
            geometry_form.append(('SurfaceID', 3, []))


        #ntris = model.tris.shape[0]
        #nquads = model.quads.shape[0]
        #nelements = ntris + nquads
        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        #grid_bcs = element_props[:, 2]

        #npids = len(model.pids)
        pids = model.pids
        eid_res = GuiResult(0, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids, colormap=colormap)
        nid_res = GuiResult(0, header='NodeID', title='NodeID',
                            location='node', scalar=nids, colormap=colormap)
        nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                nlabels=2, labelsize=5, ncolors=2,
                                colormap=colormap, data_format='%.1f',
                                uname='NormalResult')

        icase = 0
        cases[icase] = (eid_res, (0, 'ElementID'))
        cases[icase + 1] = (nid_res, (0, 'NodeID'))
        cases[icase + 2] = (nxyz_res, (0, 'Normals'))
        icase += 3
        if not read_solids:
            surface_res = GuiResult(0, header='SurfaceID', title='SurfaceID',
                                    location='centroid', scalar=pids, colormap=colormap)
            cases[icase] = (surface_res, (0, 'SurfaceID'))
            icase += 1

        if tag_filename is not None and os.path.exists(tag_filename) and not read_solids:
            #surf_ids = element_props[:, 0]
            #recon_flags = element_props[:, 1]
            #cases[(ID, 2, 'ReconFlag', 1, 'centroid', '%i')] = recon_flags
            #cases[(ID, 3, 'GridBC',    1, 'centroid', '%i')] = grid_bcs

            tagger = TagReader()
            data = tagger.read_tag_filename(tag_filename)

            int_data = ones((nelements, 8), dtype='int32') * -10.
            float_data = zeros((nelements, 2), dtype='float64')
            for key, datai in sorted(data.items()):
                #self.gui.log.info(datai)
                [name, is_visc, is_recon, is_rebuild, is_fixed, is_source,
                 is_trans, is_delete, bl_spacing, bl_thickness, nlayers] = datai
                i = where(pids == key)[0]
                int_data[i, :] = [is_visc, is_recon, is_rebuild, is_fixed,
                                  is_source, is_trans, is_delete, nlayers]
                float_data[i, :] = [bl_spacing, bl_thickness]
                self.gui.log.info('data[%i] = %s' % (key, name))

            has_tag_data = True
            tag_form = []
            tag_form.append(('is_visc', icase, []))
            tag_form.append(('is_recon', icase+1, []))
            tag_form.append(('is_rebuild', icase+2, []))
            tag_form.append(('is_fixed', icase+3, []))
            tag_form.append(('is_source', icase+4, []))
            tag_form.append(('is_trans', icase+5, []))
            tag_form.append(('is_delete', icase+6, []))
            tag_form.append(('nlayers', icase+7, []))
            tag_form.append(('bl_spacing', icase+8, []))
            tag_form.append(('bl_thickness', icase+9, []))

            visc_res = GuiResult(0, header='is_visc', title='is_visc',
                                 location='node', scalar=int_data[:, 0], colormap=colormap)
            recon_res = GuiResult(0, header='is_recon', title='is_recon',
                                  location='node', scalar=int_data[:, 1], colormap=colormap)
            rebuild_res = GuiResult(0, header='is_rebuild', title='is_rebuild',
                                    location='node', scalar=int_data[:, 2], colormap=colormap)
            fixed_res = GuiResult(0, header='is_fixed', title='is_fixed',
                                  location='node', scalar=int_data[:, 3], colormap=colormap)
            source_res = GuiResult(0, header='is_source', title='is_source',
                                   location='node', scalar=int_data[:, 4], colormap=colormap)
            trans_res = GuiResult(0, header='is_trans', title='is_trans',
                                  location='node', scalar=int_data[:, 5], colormap=colormap)
            delete_res = GuiResult(0, header='is_delete', title='is_delete',
                                   location='node', scalar=int_data[:, 6], colormap=colormap)
            nlayers_res = GuiResult(0, header='nlayers', title='nlayers',
                                    location='node', scalar=int_data[:, 7], colormap=colormap)

            spacing_res = GuiResult(0, header='bl_spacing', title='bl_spacing',
                                    location='centroid', scalar=float_data[:, 0],
                                    colormap=colormap)
            blthickness_res = GuiResult(0, header='bl_thickness', title='bl_thickness',
                                        location='centroid', scalar=float_data[:, 1],
                                        colormap=colormap)

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

            icase += 10
        elif tag_filename is not None:
            self.gui.log.warning('tag_filename=%r could not be found' % tag_filename)

        if mapbc_filename is not None and os.path.exists(mapbc_filename) and not read_solids:
            has_mapbc_data = True
            with open(mapbc_filename, 'r') as mapbc:
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
                self.gui.log.info(line)
            mapbc_form.append(('Map BC', icase, []))

            mapbc_res = GuiResult(0, header='Map BC', title='Map BC',
                                  location='centroid', scalar=mapbcs, colormap=colormap)
            cases[icase + 9] = (mapbc_res, (0, 'Map BC'))
        elif mapbc_filename is not None:
            self.gui.log.warning('mapbc_filename=%r could not be found' % mapbc_filename)


        #norm_spacing = model.node_props[:, 0]
        #bl_thickness = model.node_props[:, 1]
        #cases[(ID, 1, 'normSpacing', 1, 'node', '%.3e', '')] = norm_spacing
        #cases[(ID, 2, 'BL_thick',    1, 'node', '%.3e', '')] = bl_thickness

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
        self.gui.log.info(form)
        return form, cases, nids, eids

def get_ugrid_model(ugrid_filename, read_solids, log):
    """helper method for UGRID_IO"""
    if read_solids or is_binary_file(ugrid_filename):
        model = UGRID(log=log, debug=True, read_solids=read_solids)
        ext = os.path.basename(ugrid_filename).split('.')[2] # base, fmt, ext
        is_2d = False
    else:
        ext = os.path.basename(ugrid_filename).split('.')[1] # base, ext
        model = UGRID2D_Reader(log=log, debug=True)
        is_2d = True
    is_3d = not is_2d

    assert ext == 'ugrid', ugrid_filename
    return model, is_2d, is_3d
