"""
Defines the GUI IO file for Tecplot.
"""
from __future__ import print_function
#import os
from six import integer_types
from collections import OrderedDict

from numpy import arange, mean, amax, amin, array
from vtk import vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra

from pyNastran.converters.tecplot.tecplot import read_tecplot
#from pyNastran.converters.tecplot.utils import merge_tecplot_files
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points


class TecplotIO(object):
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
        self.gui.nnodes = model.nnodes

        #self._make_tecplot_geometry(model, self.nnodes, quads_only=True) # cart3d
        is_surface = self._make_tecplot_geometry(model, quads_only=False)

        #self._create_cart3d_free_edegs(model, nodes, elements)


        # loadTecplotResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        loads = []
        assert loads is not None
        if 'Mach' in loads:
            avg_mach = mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''
        self.gui.isubcase_name_map = {1: ['Tecplot%s' % note, '']}
        cases = OrderedDict()
        ID = 1

        form, cases, node_ids, element_ids = self._fill_tecplot_case(cases, ID, model, is_surface)
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
        nodes = model.xyz
        unused_nnodes = self.gui.nnodes
        grid = self.gui.grid

        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)

        points = numpy_to_vtk_points(nodes)

        #elements = model.elements
        quads = model.quad_elements
        hexas = model.hexa_elements
        tets = model.tet_elements
        tris = model.tri_elements

        nquads = len(quads)
        ntris = len(tris)
        nhexas = len(hexas)
        ntets = len(tets)

        nshells = nquads + ntris
        nsolids = ntets + nhexas
        if nshells:
            is_surface = True
            self._create_tecplot_shells(nquads, quads, ntris, tris)
            self.gui.nelements = nshells

        elif nsolids:
            #if 0:
                #tris, quads = model.skin_elements()
                #is_tris = bool(len(tris))
                #is_quads = bool(len(quads))
                #self._create_tecplot_shells(is_quads, quads, is_tris, tris)
            #else:
            is_surface = False
            if is_surface:
                if nhexas:
                    free_faces = array(model.get_free_faces(), dtype='int32')# + 1
                    nfaces = len(free_faces)
                    self.gui.nelements = nfaces
                    unused_elements = free_faces
                    grid.Allocate(nfaces, 1000)

                    for face in free_faces:
                        elem = vtkQuad()
                        epoints = elem.GetPointIds()
                        epoints.SetId(0, face[0])
                        epoints.SetId(1, face[1])
                        epoints.SetId(2, face[2])
                        epoints.SetId(3, face[3])
                        #elem.GetCellType() = 5  # vtkTriangle
                        grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            else:
                # is_volume
                grid.Allocate(nsolids, 1000)
                self.gui.nelements = nsolids
                if ntets:
                    for node_ids in tets:
                        elem = vtkTetra()
                        epoints = elem.GetPointIds()
                        epoints.SetId(0, node_ids[0])
                        epoints.SetId(1, node_ids[1])
                        epoints.SetId(2, node_ids[2])
                        epoints.SetId(3, node_ids[3])
                        #elem.GetCellType() = 5  # vtkTriangle
                        grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())


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
                        grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        else:
            raise NotImplementedError()

        grid.SetPoints(points)
        grid.Modified()
        return is_surface

    def _create_tecplot_shells(self, is_quads, quads, is_tris, tris):
        grid = self.gui.grid
        if is_quads:
            for face in quads:
                elem = vtkQuad()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0])
                epoints.SetId(1, face[1])
                epoints.SetId(2, face[2])
                epoints.SetId(3, face[3])
                #elem.GetCellType() = 5  # vtkTriangle
                grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

        if is_tris:
            for face in tris:
                elem = vtkTriangle()
                epoints = elem.GetPointIds()
                epoints.SetId(0, face[0])
                epoints.SetId(1, face[1])
                epoints.SetId(2, face[2])
                #elem.GetCellType() = 5  # vtkTriangle
                grid.InsertNextCell(5, elem.GetPointIds())

    def clear_tecplot(self):
        pass

    #def load_tecplot_results(self, cart3d_filename):
        #model = Cart3D(log=self.log, debug=False)
        #self.load_cart3d_geometry(cart3d_filename)

    def _fill_tecplot_case(self, cases, ID, model, is_surface):
        #'x', 'y', 'z',
        #result_names = ['rho', 'U', 'V', 'W', 'p']
        result_names = model.variables[3:]
        #nelements = elements.shape[0]
        nelements = self.gui.nelements
        nnodes = model.nnodes
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
        assert isinstance(nnodes, integer_types), 'nnodes=%s type=%s' % (nnodes, type(nnodes))
        assert isinstance(nelements, integer_types), 'nelements=%s type=%s' % (nelements, type(nelements))
        assert nnodes > 0, nnodes
        assert nelements > 0, nelements

        nids = arange(1, nnodes + 1, dtype='int32')
        eids = arange(1, nelements + 1, dtype='int32')

        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        eid_res = GuiResult(ID, header=element_id, title=element_id,
                            location='centroid', scalar=eids)

        icase = 0
        cases[icase] = (nid_res, (0, 'NodeID'))
        cases[icase + 1] = (eid_res, (0, element_id))
        icase += 2

        results = model.results
        if is_results and len(results):
            for iresult, result_name in enumerate(result_names):
                if results.shape[1] == 1:
                    nodal_data = results
                    assert len(result_names) == 1, result_names
                else:
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
