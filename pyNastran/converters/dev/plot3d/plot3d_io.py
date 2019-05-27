from collections import OrderedDict
import vtk
from vtk import vtkQuad

from numpy import zeros, array, cross
from numpy.linalg import norm  # type: ignore

from pyNastran.converters.dev.plot3d.plot3d import Plot3d
from pyNastran.gui.gui_objects.gui_result import GuiResult

raise NotImplementedError()

class Plot3d_io:  # pragma: no cover
    def __init__(self, gui):
        self.gui = gui

    def get_plot3d_wildcard_geometry_results_functions(self):
        data = ('Plot3D',
                'Plot3D (*.p3d; *.p3da)', self.load_plot3d_geometry,
                None, None)
        return data

    def load_plot3d_geometry(self, p3d_filename, name='main'):
        print("load_plot3d_geometry")
        self.nid_map = {}

        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_geometry(p3d_filename)
        if skip_reading:
            return

        model = Plot3d(log=self.log, debug=self.debug)
        #self.model_type = model.model_type
        model.read_plot3d(p3d_filename)

        npoints = 0
        nelements = 0
        for iblock, shape in sorted(model.block_shapes.items()):
            npoints += shape[0] * shape[1] * shape[2]
            nelements += (shape[0] - 1)  * (shape[1] - 1) * (shape[2] - 1)

        nblocks = len(model.block_shapes)
        self.nnodes = npoints
        self.nelements = nelements


        #nodes, elements, regions = model.get_points_elements_regions()
        #for nid,node in enumerate(nodes):
            #print "node[%s] = %s" % (nid, str(node))

        self.grid.Allocate(self.nelements, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nnodes)

        nid = 0
        nid_base = 0

        eid_base = 0

        elem = vtkQuad()
        quad_type = elem.GetCellType()
        #nblocks = len(model.x)
        for iblock in range(nblocks):
            print("iblock =", iblock)
            nid_base = nid
            x = model.x[iblock]
            y = model.y[iblock]
            z = model.z[iblock]
            print("x.shape[%s] =" % iblock, x.shape)
            print(x)
            (ni, nj, nk) = x.shape
            assert nk == 1
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        points.InsertPoint(
                            nid, x[i, j, 0], y[i, j, 0], z[i, j, 0])
                        nid += 1

            for j in range(nj - 1):
                jstart = nid_base + j * ni
                for i in range(ni - 1):
                    elem = vtkQuad()

                    p1 = jstart + (i)
                    p2 = jstart + (i + 1)
                    p3 = jstart + (ni) + (i + 1)
                    p4 = jstart + (ni) + (i)

                    elem.GetPointIds().SetId(0, p1)
                    elem.GetPointIds().SetId(1, p2)
                    elem.GetPointIds().SetId(2, p3)
                    elem.GetPointIds().SetId(3, p4)
                    element = [p1, p2, p3, p4]
                    self.grid.InsertNextCell(quad_type, elem.GetPointIds())
                    print(element)
                #jstart += ni

            #nid_base += ni * nj * nk
            eid_base += (ni-1) * (nj-1) * (nk-1)
            break

        #print("eid = ", eid)
        grid = self.gui.grid
        grid.SetPoints(points)
        grid.Modified()
        grid.Update()
        self.log_info("updated grid")

        #return

        # loadPlot3dResults - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        self.gui.isubcase_name_map = {1: ['Plot3d', '']}
        cases = OrderedDict()
        ID = 1

        #cases = self._fill_stl_case(cases, ID, elements)
        self.result_cases = cases
        self.case_keys = sorted(cases.keys())
        #print "case_keys = ",self.case_keys
        #print "type(case_keys) = ",type(self.case_keys)
        self.ncases = min(0, len(self.result_cases) - 1)  # number of keys in dictionary
        self.icase = 0 if self.ncases == 0 else -1
        self.cycle_results()  # start at ncase=0

    def fill_plot3d_geometry_case(self, cases, ID, nodes, elements, regions, loads):
        #print "regions**** = ",regions
        #nNodes = self.nnodes
        #nElements = self.nelements

        #result_names = ['Cp', 'Mach', 'U', 'V', 'W', 'E', 'rho',
                                      #'rhoU', 'rhoV', 'rhoW', 'rhoE']
        #nelements, three = elements.shape
        #print regions
        icase = 0
        itime = 0

        region_res = GuiResult(ID, header='Region', title='Region',
                               location='centroid', scalar=regions)
        cases[icase] = region_res(nid_res, (itime, 'Region'))
        icase += 1

        # centroidal
        Xc = zeros(len(elements), 'float64')
        Yc = zeros(len(elements), 'float64')
        Zc = zeros(len(elements), 'float64')
        area = zeros(len(elements), 'float64')

        Xn = zeros(len(nodes), 'float64')
        Yn = zeros(len(nodes), 'float64')
        Zn = zeros(len(nodes), 'float64')

        for i, element in enumerate(elements):
            p1, p2, p3, p4 = element
            P1 = array(nodes[p1])
            P2 = array(nodes[p2])
            P3 = array(nodes[p3])
            P4 = array(nodes[p4])
            a = P3 - P1
            b = P4 - P2
            A = 0.5 * norm(cross(a, b))
            x, y, z = (P1 + P2 + P3 + P4) / 4.0
            Xc[i] = x
            Yc[i] = y
            Zc[i] = z
            area[i] = A
        for i, node in enumerate(nodes):
            Xn[i] = node[0]
            Yn[i] = node[1]
            Zn[i] = node[2]

        cases[(ID, icase, 'node_x', 1, 'node', '%.2f', '')] = Xn
        cases[(ID, icase + 1, 'node_y', 1, 'node', '%.2f', '')] = Yn
        cases[(ID, icase + 2, 'node_z', 1, 'node', '%.2f', '')] = Zn

        cases[(ID, icase + 3, 'centroid_x', 1, 'centroid', '%.2f', '')] = Xc
        cases[(ID, icase + 4, 'centroid_y', 1, 'centroid', '%.2f', '')] = Yc
        cases[(ID, icase + 5, 'centroid_z', 1, 'centroid', '%.2f', '')] = Zc

        cases[(ID, icase + 6, 'Area', 1, 'centroid', '%.2f', '')] = area
        icase += 7

        #print("load.keys() = ", loads.keys())
        #for key in result_names:
            #if key in loads:
                #nodal_data = loads[key]
                #cases[(ID, key, 1, 'node', '%.3f')] = nodal_data
        return cases
