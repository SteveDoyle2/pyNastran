from numpy import arange, mean, amax, amin, array

from pyNastran.gui.vtk_interface import vtkHexahedron, vtkQuad, vtkTriangle, vtkTetra

from pyNastran.converters.dev.avus.avus_grid import AvusGrid
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points

class AvusIO:
    def __init__(self, gui):
        self.gui = gui

    def get_avus_wildcard_geometry_results_functions(self):
        data = ('Avus',
                'Avus (*.grd)', self.load_avus_geometry,
                None, None)
        return data

    #def removeOldGeometry(self, filename):
        #self._remove_old_cart3d_geometry(filename)

    def _remove_old_cart3d_geometry(self, filename):
        #return self._remove_old_geometry(filename)

        self.gui.eid_map = {}
        self.gui.nid_map = {}
        if filename is None:
            self.gui.scalar_bar_actor.VisibilityOff()
            skip_reading = True
        else:
            self.gui.turn_corner_text_off()
            self.gui.grid.Reset()

            self.gui.result_cases = {}
            self.gui.ncases = 0
            try:
                del self.gui.case_keys
                del self.gui.icase
                del self.gui.isubcase_name_map
            except Exception:
                # print("cant delete geo")
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()
        return skip_reading

    def load_avus_geometry(self, avus_filename, name='main', plot=True):
        model_name = name
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_cart3d_geometry(avus_filename)
        if skip_reading:
            return

        model = AvusGrid(log=self.gui.log, debug=False)
        model.read_avus_grid(avus_filename)

        self.model_type = 'avus'
        #self.model_type = model.model_type
        self.nnodes = model.nnodes

        #self._make_tecplot_geometry(model, self.nnodes, quads_only=True) # cart3d
        is_surface = self._make_avus_geometry(model, quads_only=False)

        #self._create_cart3d_free_edegs(model, nodes, elements)


        # loadAvusResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        loads = []
        assert loads is not None
        if 'Mach' in loads:
            avg_mach = mean(loads['Mach'])
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''
        self.gui.isubcase_name_map = {1: ['Avus%s' % note, '']}
        cases = {}
        ID = 1

        form, cases, node_ids, element_ids = self._fill_avus_case(cases, ID, model, is_surface)
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


    def _make_avus_geometry(self, model, quads_only=False):
        nodes = model.nodes
        #nnodes = self.nnodes

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

        is_shells = nquads + ntris
        is_solids = ntets + nhexas
        if is_shells:
            is_surface = True
            if nquads:
                elements = quads
                for face in quads:
                    elem = vtkQuad()
                    epoints = elem.GetPointIds()
                    epoints.SetId(0, face[0])
                    epoints.SetId(1, face[1])
                    epoints.SetId(2, face[2])
                    epoints.SetId(3, face[3])
                    #elem.GetCellType() = 5  # vtkTriangle
                    grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            if ntris:
                elements = tris
                for face in tris:
                    elem = vtkTriangle()
                    epoints = elem.GetPointIds()
                    epoints.SetId(0, face[0])
                    epoints.SetId(1, face[1])
                    epoints.SetId(2, face[2])
                    #elem.GetCellType() = 5  # vtkTriangle
                    grid.InsertNextCell(5, elem.GetPointIds())

        elif is_solids:
            if ntets:
                elements = tets
                is_surface = False
                self.nelements = model.nelements
                grid.Allocate(self.nelements, 1000)

                nelements = elements.shape[0]
                for node_ids in elements:
                    elem = vtkTetra()
                    epoints = elem.GetPointIds()
                    epoints.SetId(0, node_ids[0])
                    epoints.SetId(1, node_ids[1])
                    epoints.SetId(2, node_ids[2])
                    epoints.SetId(3, node_ids[3])
                    #elem.GetCellType() = 5  # vtkTriangle
                    grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())


            if nhexas:
                elements = hexas
                is_surface = True
                # is_surface = False
                is_volume = not is_surface

                if is_surface:
                    self.nelements = model.nelements
                    free_faces = array(model.get_free_faces(), dtype='int32')# + 1
                    nfaces = len(free_faces)
                    elements = free_faces
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

                elif is_volume:
                    self.nelements = model.nelements
                    grid.Allocate(self.nelements, 1000)

                    nelements = elements.shape[0]
                    for eid in range(nelements):
                        elem = vtkHexahedron()
                        node_ids = elements[eid, :]
                        epoints = elem.GetPointIds()
                        epoints.SetId(0, node_ids[0])
                        epoints.SetId(1, node_ids[1])
                        epoints.SetId(2, node_ids[2])
                        epoints.SetId(3, node_ids[3])
                        epoints.SetId(4, node_ids[4])
                        epoints.SetId(5, node_ids[5])
                        epoints.SetId(6, node_ids[6])
                        epoints.SetId(7, node_ids[7])
                        #elem.GetCellType() = 8  # vtkHexa
                        grid.InsertNextCell(8, elem.GetPointIds())
        else:
            raise NotImplementedError()

        grid.SetPoints(points)
        grid.Modified()
        return is_surface

    def clear_avus(self):
        pass

    #def load_tecplot_results(self, cart3d_filename):
        #model = Cart3D(log=self.log, debug=False)
        #self.load_cart3d_geometry(cart3d_filename)

    def _fill_avus_case(self, cases, ID, model, is_surface):
        #'x', 'y', 'z',
        #result_names = ['rho', 'U', 'V', 'W', 'p']
        #result_names = []
        # result_names = model.variables[3:]
        #nelements = elements.shape[0]
        nelements = model.nelements
        nnodes = model.nnodes
        #nnodes = nodes.shape[0]


        #cases_new = []
        #new = False

        #is_results = False
        #is_results = True
        #results_form = []

        if is_surface:
            element_id = 'FaceID'
        else:
            element_id = 'ElementID'

        geometry_form = [
            ('NodeID', 0, []),
            (element_id, 1, []),
            #('Region', 2, []),
        ]

        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes + 1)

        eid_res = GuiResult(ID, header=element_id, title=element_id,
                            location='centroid', scalar=eids)
        nid_res = GuiResult(ID, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        icase = 0
        itime = 0
        cases[icase] = (eid_res, (itime, element_id))
        cases[icase + 1] = (nid_res, (itime, 'NodeID'))

        #cases[(ID, 2, 'Region', 1, 'centroid', '%i')] = regions

        return geometry_form, cases, nids, eids

        #results = model.results
        #if is_results and len(results):
            #i = 2
            #for iresult, result_name in enumerate(result_names):
                #if results.shape[1] == 1:
                    #nodal_data = results
                    #assert len(result_names) == 1, result_names
                #else:
                    #nodal_data = results[:, iresult]

                #if new:
                    #cases_new[i] = (result, i, result_name, 1, 'node', '%.3f', '')
                #else:
                    #cases[(ID, i, result_name, 1, 'node', '%.3f', '')] = nodal_data
                #results_form.append((result_name, i, []))
                #i += 1
        #form = [
            #('Geometry', None, geometry_form),
        #]
        #if len(results_form):
            #form.append(('Results', None, results_form))
        #return form, cases
