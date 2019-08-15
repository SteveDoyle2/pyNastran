"""Defines the GUI IO file for Tecplot."""
from collections import OrderedDict

import numpy as np
#from numpy import arange, mean, amax, amin, array
from vtk import vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra

from pyNastran.converters.tecplot.tecplot import read_tecplot
#from pyNastran.converters.tecplot.utils import merge_tecplot_files
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points


class TecplotIO:
    def __init__(self, gui):
        self.gui = gui

    def _remove_old_cart3d_geometry(self, tecplot_filename):
        pass

    def get_tecplot_wildcard_geometry_results_functions(self):
        data = ('Tecplot Binary FEBlock',
                'Tecplot Binary FEBlock (*.dat; *.plt; *.tec)', self.load_tecplot_geometry,
                None, None)
        return data

    def load_tecplot_geometry(self, tecplot_filename, name='main', plot=True):
        model_name = name
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_cart3d_geometry(tecplot_filename)
        #skip_reading = False
        if skip_reading:
            return

        #if 0:
            #fnames = os.listdir('time20000')
            #fnames = [os.path.join('time20000', fname) for fname in fnames]
            #model = merge_tecplot_files(fnames, tecplot_filename_out=None, log=self.log)
        #else:
        model = read_tecplot(tecplot_filename, log=self.gui.log, debug=False)

        self.gui.model_type = 'tecplot'
        self.gui.nnodes = sum([zone.nnodes for zone in model.zones])
        variables = None
        for zone in model.zones:
            variables = zone.variables
            break

        #self._make_tecplot_geometry(model, self.nnodes, quads_only=True) # cart3d
        is_surface = self._make_tecplot_geometry(model, quads_only=False)

        #self._create_cart3d_free_edegs(model, nodes, elements)


        # loadTecplotResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        loads = []
        assert loads is not None
        if 'Mach' in loads:
            avg_mach = np.mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''
        self.gui.isubcase_name_map = {1: ['Tecplot%s' % note, '']}
        cases = OrderedDict()
        ID = 1

        form, cases, node_ids, element_ids = self._fill_tecplot_case(cases, ID, model, variables, is_surface)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

        #if 0:
            # http://www.vtk.org/Wiki/VTK/Examples/Cxx/Filtering/AppendFilter
            #points = vtkAppendFilter()
            #if VTK_MAJOR_VERSION <= 5:
                #appendFilter.AddInput(polydata)
                #appendFilter.AddInput(ug)
            #else:
            #appendFilter.AddInputData(polydata)
            #appendFilter.AddInputData()
            #appendFilter.Update()

    def _make_tecplot_geometry(self, model, quads_only=False):
        """
        Returns
        -------
        is_surface : bool
            the model is made up of only shells (not 100%)
        """
        grid = self.gui.grid
        nnodes = self.gui.nnodes
        nodes = np.zeros((nnodes, 3), dtype='float32')
        #print('nnodes=', nnodes)

        inode = 0

        quads = []
        tris = []
        tets = []
        hexas = []
        for zone in model.zones:
            nodes2d = zone.xy
            nodes3d = zone.xyz
            nnodes2d = nodes2d.shape[0]
            nnodes3d = nodes3d.shape[0]

            # elements
            if 'I' in zone.headers_dict:
                i = zone.headers_dict['I']
                if 'J' in zone.headers_dict:
                    j = zone.headers_dict['J']
                    if 'K' in zone.headers_dict:
                        k = zone.headers_dict['K']
                        nnodes = i * j * k
                        elements = np.arange(0, nnodes).reshape(k, j, i)
                        #print(elements[0, :, :])
                        #print(elements[1:, :, :])
                        n1 = elements[:-1, :-1, :-1].ravel()
                        n2 = elements[:-1, :-1, 1:].ravel()
                        n3 = elements[:-1, 1:, 1:].ravel()
                        n4 = elements[:-1, 1:, :-1].ravel()

                        n5 = elements[1:, :-1, :-1].ravel()
                        n6 = elements[1:, :-1, 1:].ravel()
                        n7 = elements[1:, 1:, 1:].ravel()
                        n8 = elements[1:, 1:, :-1].ravel()
                        #nhexas = (i - 1) * (j - 1) * (k - 1)
                        quadsi = []
                        hexasi = np.vstack([n1, n2, n3, n4, n5, n6, n7, n8]).T
                    else:
                        nnodes = i * j
                        elements = np.arange(0, nnodes).reshape(j, i)
                        #print('elements:')
                        #print(elements)
                        n1 = elements[:-1, :-1].ravel()
                        n2 = elements[:-1, 1:].ravel()
                        n3 = elements[1:, 1:].ravel()
                        n4 = elements[1:, :-1].ravel()
                        #nquads = (i - 1) * (j - 1)
                        quadsi = np.vstack([n1, n2, n3, n4]).T
                        hexasi = []
                    #trisi = None
                    #tetsi = None
                    #hexasi = None
                    trisi = []
                    tetsi = []
            else:
                quadsi = zone.quad_elements
                hexasi = zone.hexa_elements
                tetsi = zone.tet_elements
                trisi = zone.tri_elements

            if len(quadsi):
                quads.append(inode + quadsi)
            if len(trisi):
                tris.append(inode + trisi)
            if len(tetsi):
                tets.append(inode + tetsi)
            if len(hexasi):
                hexas.append(inode + hexasi)

            # nodes
            if nnodes2d and nnodes3d:
                raise RuntimeError('2d and 3d nodes is not supported')
            elif nnodes2d:
                nodes[inode:inode+nnodes2d, :2] = nodes2d
                inode += nnodes2d
            elif nnodes3d:
                nodes[inode:inode+nnodes3d] = nodes3d
                inode += nnodes3d
            else:  # pragma: no cover
                raise RuntimeError('failed to find 2d/3d nodes')
            #print('inode', inode)
            #print(nodes)
            #print('-------------')
        #print('stack', len(quads))
        quads = stack(quads)
        tris = stack(tris)
        tets = stack(tets)
        hexas = stack(hexas)

        nquads = len(quads)
        ntris = len(tris)
        nhexas = len(hexas)
        ntets = len(tets)

        nshells = nquads + ntris
        nsolids = ntets + nhexas
        if nshells:
            is_surface = True
            grid = self.gui.grid
            _create_tecplot_shells(grid, nquads, quads, ntris, tris)
            self.gui.nelements = nshells

        elif nsolids:
            #if 0:
                #tris, quads = model.skin_elements()
                #is_tris = bool(len(tris))
                #is_quads = bool(len(quads))
                #self._create_tecplot_shells(is_quads, quads, is_tris, tris)
            #else:
            is_surface = False
            grid = self.gui.grid
            nelements = _create_tecplot_solids(grid, zone, nsolids, ntets, tets, nhexas, hexas,
                                               is_surface=is_surface)
            self.gui.nelements = nelements
        else:
            raise NotImplementedError()

        #----------------------------------------------
        #print('nnodes', nnodes, inode)
        mmax = np.amax(nodes, axis=0)
        mmin = np.amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)

        #print('nodes', nodes)
        #nnodes = nodes.shape[0]
        points = numpy_to_vtk_points(nodes)

        grid.SetPoints(points)
        grid.Modified()
        return is_surface

    def clear_tecplot(self):
        pass

    #def load_tecplot_results(self, cart3d_filename):
        #model = Cart3D(log=self.log, debug=False)
        #self.load_cart3d_geometry(cart3d_filename)

    def _fill_tecplot_case(self, cases, ID, model, variables, is_surface):
        #'x', 'y', 'z',
        #result_names = ['rho', 'U', 'V', 'W', 'p']
        #result_names = zone.variables[3:]
        #nelements = elements.shape[0]
        nelements = self.gui.nelements
        nnodes = self.gui.nnodes
        #nnodes = zone.nnodes
        #nnodes = nodes.shape[0]

        #is_results = False
        is_results = True
        results_form = []

        if is_surface:
            element_id = 'FaceID'
        else:
            element_id = 'ElementID'

        geometry_form = [
            ('NodeID', 0, []),
            (element_id, 1, []),
        ]
        assert isinstance(nnodes, int), 'nnodes=%s type=%s' % (nnodes, type(nnodes))
        assert isinstance(nelements, int), 'nelements=%s type=%s' % (nelements, type(nelements))
        assert nnodes > 0, nnodes
        assert nelements > 0, nelements

        nids = np.arange(1, nnodes + 1, dtype='int32')
        eids = np.arange(1, nelements + 1, dtype='int32')

        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        eid_res = GuiResult(ID, header=element_id, title=element_id,
                            location='centroid', scalar=eids)

        icase = 0
        cases[icase] = (nid_res, (0, 'NodeID'))
        cases[icase + 1] = (eid_res, (0, element_id))
        icase += 2

        nvars = len(variables)
        if model.nzones == 1:
            zone0 = model.zones[0]
            results = zone0.nodal_results
        else:
            results = []
            for zonei in model.zones:
                results.append(zonei.nodal_results)
            results = np.vstack(results)
            assert results.shape == (nnodes, nvars), results.shape

        if is_results and nvars:
            for iresult, result_name in enumerate(variables):
                #if results.shape[1] == 1:
                    #nodal_data = results
                    #assert len(variables) == 1, variables
                #else:
                nodal_data = results[:, iresult]

                node_res = GuiResult(ID, header=result_name, title=result_name,
                                     location='node', scalar=nodal_data)
                cases[icase] = (node_res, (0, result_name))

                results_form.append((result_name, icase, []))
                icase += 1
        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases, nids, eids

def _create_tecplot_shells(grid, is_quads, quads, is_tris, tris):
    if is_quads:
        #elem.GetCellType() = 9  # vtkQuad
        for face in quads:
            elem = vtkQuad()
            epoints = elem.GetPointIds()
            epoints.SetId(0, face[0])
            epoints.SetId(1, face[1])
            epoints.SetId(2, face[2])
            epoints.SetId(3, face[3])
            grid.InsertNextCell(9, epoints)

    if is_tris:
        #elem.GetCellType() = 5  # vtkTriangle
        for face in tris:
            elem = vtkTriangle()
            epoints = elem.GetPointIds()
            epoints.SetId(0, face[0])
            epoints.SetId(1, face[1])
            epoints.SetId(2, face[2])
            grid.InsertNextCell(5, epoints)

def _create_tecplot_solids(grid, model, nsolids, ntets, tets, nhexas, hexas, is_surface=True):
    """
    add a model with solid elements

    Parameters
    ----------
    is_surface : bool; default=True
        True : skin the model (good for large models, but doesn't load everything)
        False : load the model normally
    """
    if is_surface:
        if nhexas:
            free_faces = np.array(model.get_free_faces(), dtype='int32')# + 1
            nfaces = len(free_faces)
            nelements = nfaces
            unused_elements = free_faces
            grid.Allocate(nfaces, 1000)

            #elem.GetCellType() = 9  # vtkQuad
            for face in free_faces:
                elem = vtkQuad()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0])
                epoints.SetId(1, face[1])
                epoints.SetId(2, face[2])
                epoints.SetId(3, face[3])
                grid.InsertNextCell(9, epoints)
    else:
        # is_volume
        grid.Allocate(nsolids, 1000)
        nelements = nsolids
        if ntets:
            for node_ids in tets:
                elem = vtkTetra()
                epoints = elem.GetPointIds()
                epoints.SetId(0, node_ids[0])
                epoints.SetId(1, node_ids[1])
                epoints.SetId(2, node_ids[2])
                epoints.SetId(3, node_ids[3])
                #elem.GetCellType() = 5  # vtkTriangle
                grid.InsertNextCell(elem.GetCellType(), epoints)


        if nhexas:
            for node_ids in hexas:
                elem = vtkHexahedron()
                epoints = elem.GetPointIds()
                epoints.SetId(0, node_ids[0])
                epoints.SetId(1, node_ids[1])
                epoints.SetId(2, node_ids[2])
                epoints.SetId(3, node_ids[3])
                epoints.SetId(4, node_ids[4])
                epoints.SetId(5, node_ids[5])
                epoints.SetId(6, node_ids[6])
                epoints.SetId(7, node_ids[7])
                #elem.GetCellType() = 5  # vtkTriangle
                grid.InsertNextCell(elem.GetCellType(), epoints)
    assert nelements > 0, nelements
    return nelements


def stack(elements):
    if len(elements) == 0:
        pass
    elif len(elements) == 1:
        elements = elements[0]
        #print(elements)
    else:
        #print('----stack------')
        #for elementsi in elements:
            #print(elementsi)
        #print('----stack------')

        elements = np.vstack(elements)
    return elements
