from numpy import arange, amax, amin, zeros, hstack

from pyNastran.gui.vtk_common_core import vtkPoints
from pyNastran.gui.vtk_interface import vtkLine, vtkTriangle

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.dev.openvsp.adb_reader import ADB_Reader


class ADB_IO:  # pragma: no cover
    def __init__(self, gui):
        self.gui = gui

    def get_openvsp_wildcard_geometry_results_functions(self):  # pragma: no cover
        data = ('VSPAero',
                'VSPAero (*.adb)', self.load_vsp_aero_geometry,
                #'Cart3d (*.triq)', self.load_cart3d_results,
                None, None
               )
        return data

    def _remove_old_adb_geometry(self, adb_filename):  # pragma: no cover
        self.gui.eid_map = {}
        self.gui.nid_map = {}
        if adb_filename is None:
            self.gui.scalar_bar_actor.VisibilityOff()
            skip_reading = True
        else:
            self.gui.turn_text_off()
            self.gui.grid.Reset()

            self.gui.result_cases = {}
            self.gui.ncases = 0
            try:
                del self.gui.case_keys
                del self.gui.icase
                del self.gui.isubcase_name_map
            except Exception:
                # print('cant delete geo')
                pass

            #print(dir(self))
            skip_reading = False
        #self.scalar_bar_actor.VisibilityOff()
        self.gui.scalar_bar_actor.Modified()
        return skip_reading

    def load_vsp_aero_geometry(self, adb_filename,
                               name='main', plot=True):  # pragma: no cover
        model_name = name
        #key = self.case_keys[self.icase]
        #case = self.result_cases[key]

        skip_reading = self._remove_old_adb_geometry(adb_filename)
        if skip_reading:
            return

        if 0:
            plot_wakes = False
        else:
            plot_wakes = True # doesn't work perfectly

        log = self.gui.log
        model = ADB_Reader(log=log, debug=False)
        self.gui.model_type = 'vspaero'
        #self.model_type = model.model_type
        (nodes, elements) = model.read_adb(adb_filename)
        nxyz_nodes = nodes.shape[0]
        nwake_nodes = model.wake_xyz.shape[0]

        nxyz_elements = model.tris.shape[0]
        nwake_elements = max(0, model.wake_xyz.shape[0] - 1)
        if plot_wakes:
            nnodes = nxyz_nodes + nwake_nodes
            nelements = nxyz_elements + nwake_elements
        else:
            nnodes = nxyz_nodes
            nelements = nxyz_elements
        self.gui.nnodes = nnodes
        self.gui.nelements = nelements

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = vtkPoints()
        points.SetNumberOfPoints(self.gui.nnodes)
        #vectorReselt.SetNumberOfComponents(3)
        self.gui.nid_map = {}

        assert nodes is not None

        nid = 0
        #print('nxyz_nodes=%s' % nxyz_nodes)
        mmax = amax(nodes, axis=0)
        mmin = amin(nodes, axis=0)
        dim_max = (mmax - mmin).max()
        self.gui.create_global_axes(dim_max)
        for i in range(nxyz_nodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1
        log.info('nxyz_nodes=%s nwake_nodes=%s total=%s' % (
            nxyz_nodes, nwake_nodes, nxyz_nodes + nwake_nodes))
        log.info('nxyz_elements=%s nwake_elements=%s total=%s' % (
            nxyz_elements, nwake_elements, nxyz_elements + nwake_elements))

        elements -= 1
        for eid in range(nxyz_elements):
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem.GetPointIds().SetId(0, node_ids[0])
            elem.GetPointIds().SetId(1, node_ids[1])
            elem.GetPointIds().SetId(2, node_ids[2])
            #elem.GetCellType() = 5  # vtkTriangle
            grid.InsertNextCell(5, elem.GetPointIds())

        if plot_wakes:
            for j in range(nwake_nodes):
                node = model.wake_xyz[j, :]
                points.InsertPoint(nid, node)
                nid += 1

            for wake in model.wake_elements:
                i0, i1 = wake
                #node_ids = range(i0, i1 + 1)
                #elem.SetPoints(points[i0:i1 + 1])
                #print(i0, i1)
                for ii in range(i0+1, i1):
                    elem = vtkLine()
                    elem.GetPointIds().SetId(0, ii - 1)
                    elem.GetPointIds().SetId(1, ii)
                    grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #eid += 1
                #break
            #assert i1 == nelements, 'ii=%s nelements=%s' % (ii, nelements)

        grid.SetPoints(points)
        grid.Modified()
        #self.log_info("updated grid")

        # load results - regions/loads
        self.gui.scalar_bar_actor.VisibilityOn()
        self.gui.scalar_bar_actor.Modified()

        mach = model.machs[0]
        alpha = model.alphas[0]
        beta = model.betas[0]
        note = ':  Mach=%.2f, alpha=%.1f, beta=%.1f' % (mach, alpha, beta)
        self.gui.isubcase_name_map = {1: ['OpenVSP%s' % note, '']}
        cases = {}
        ID = 1

        form, cases = self._fill_adb_case(cases, ID, model, plot_wakes)
        self.gui._finish_results_io2(model_name, form, cases)

    #def clear_adb(self):
        #pass

    #def load_adb_results(self, cart3d_filename):
        #raise NotImplementedError()


    def _fill_adb_case(self, cases, ID, model, plot_wakes=False):  # pragma: no cover
        nxyz_nodes = model.nodes.shape[0]
        nxyz_elements = model.tris.shape[0]
        nwake_nodes = model.wake_xyz.shape[0]
        nwake_elements = max(0, model.wake_xyz.shape[0] - 1)

        if plot_wakes:
            nnodes = nxyz_nodes + nwake_nodes
            nelements = nxyz_elements + nwake_elements
            #assert nnodes == 408 + 1584, nnodes
            #assert nelements == 2223, nelements
            is_normals = False
        else:
            nnodes = nxyz_nodes
            nelements = nxyz_elements
            nodes = model.nodes
            elements = model.tris
            is_normals = True

        Cp = model.Cp
        surf_id = model.surf_id
        area = model.area

        results_form = []
        geometry_form = [
            ('Region', 0, []),
            ('NodeID', 1, []),
            ('ElementID', 2, []),
            ('Area', 3, []),
        ]

        eids = arange(1, nelements + 1)
        nids = arange(1, nnodes+1)

        if nwake_elements and plot_wakes:
            fzero_pad = zeros(nwake_elements, dtype='float32')
            izero_pad = zeros(nwake_elements, dtype='int32')
            surf_id = hstack([surf_id, izero_pad])
            area = hstack([area, fzero_pad])

        assert len(eids) == nelements, len(eids)
        assert len(surf_id) == nelements, len(surf_id)
        assert len(area) == nelements, len(area)

        region_res = GuiResult(ID, 'Region', 'Region', 'centroid', surf_id)
        nid_res = GuiResult(ID, 'NodeID', 'NodeID', 'centroid', nids)
        eid_res = GuiResult(ID, 'ElementID', 'ElementID', 'centroid', eids)
        area_res = GuiResult(ID, 'Area', 'Area', 'centroid', area,
                             data_format='%.3f')

        i = 0
        cases[i] = (region_res, (0, 'Region'))
        cases[i + 1] = (nid_res, (0, 'NodeID'))
        cases[i + 2] = (eid_res, (0, 'ElementID'))
        cases[i + 3] = (area_res, (0, 'Area'))

        if is_normals:
            geometry_form.append(('Normal X', i, []))
            geometry_form.append(('Normal Y', i + 1, []))
            geometry_form.append(('Normal Z', i + 2, []))

            cnormals = model.get_normals(nodes, elements, shift_nodes=False)
            #nnormals = model.get_normals_at_nodes(nodes, elements, cnormals, shift_nodes=False)
            cnnodes = cnormals.shape[0]
            assert cnnodes == nelements, len(cnnodes)

            nx_res = GuiResult(ID, 'Normal X', 'Normal X', 'centroid', cnormals[:, 0],
                               data_format='%.3f')
            ny_res = GuiResult(ID, 'Normal Y', 'Normal Y', 'centroid', cnormals[:, 0],
                               data_format='%.3f')
            nz_res = GuiResult(ID, 'Normal Z', 'Normal Z', 'centroid', cnormals[:, 0],
                               data_format='%.3f')
            cases[i] = (nx_res, (0, 'Normal X'))
            cases[i + 1] = (ny_res, (0, 'Normal Y'))
            cases[i + 2] = (nz_res, (0, 'Normal Z'))
            i += 3

        cp_res = GuiResult(ID, 'Cp', 'Cp', 'centroid', Cp, data_format='%.3f')
        cases[i] = (cp_res, (0, 'Cp'))
        results_form.append(('Cp', i, []))
        i += 1

        form = [
            ('Geometry', None, geometry_form),
        ]
        if len(results_form):
            form.append(('Results', None, results_form))
        return form, cases
