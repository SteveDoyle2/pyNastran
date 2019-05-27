"""Defines the GUI IO file for Fast."""
import os
from collections import OrderedDict

import numpy as np
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.fast.fgrid_reader import read_fgrid
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points
from pyNastran.gui.gui_objects.gui_result import GuiResult


class FastIO:
    def __init__(self, gui):
        self.gui = gui

    def get_fast_wildcard_geometry_results_functions(self):
        data = ('FAST',
                'FAST (*.fgrid)', self.load_fast_geometry,
                None, None)
        return data

    def load_fast_geometry(self, fgrid_filename, name='main', plot=True):
        model_name = name
        skip_reading = self.gui._remove_old_geometry(fgrid_filename)
        if skip_reading:
            return

        ext = os.path.splitext(fgrid_filename)[1]
        ext = ext.lower()
        if ext == '.fgrid':
            dimension_flag = 3
        else:
            msg = 'unsupported extension=%r.  Use "fgrid".' % ext
            raise RuntimeError(msg)

        #read_loads = True
        model = read_fgrid(
            fgrid_filename,
            dimension_flag,
            log=self.gui.log, debug=False,
        )

        nodes = model.nodes

        ntris = 0
        ntets = 0
        if model.tets is None:
            dimension_flag = 2
            if model.tris is not None:
                tris = model.tris - 1
                ntris = tris.shape[0]
        else:
            dimension_flag = 3
            tets = model.tets - 1
            ntets = tets.shape[0]

        nnodes = nodes.shape[0]
        model.log.info('ntris=%s ntets=%s' % (ntris, ntets))

        #print('node0 = %s' % str(nodes[0, :]))
        #print('node%i = %s' % (1, str(nodes[1, :])))
        #print('node%i = %s' % (2, str(nodes[2, :])))
        #print('node%i = %s' % (nnodes, str(nodes[-1, :])))
        #mapbc = model.mapbc
        #loads = model.loads

        self.gui.nnodes = nnodes
        nelements = ntris + ntets
        self.gui.nelements = nelements

        #print("nnodes = %i" % self.nnodes)
        #print("nelements = %i" % self.nelements)

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
        self.gui.nid_map = {}

        assert nodes is not None
        nnodes = nodes.shape[0]

        if dimension_flag == 2:
            for (n0, n1, n2) in tris:
                elem = vtkTriangle()
                elem.GetPointIds().SetId(0, n0)
                elem.GetPointIds().SetId(1, n1)
                elem.GetPointIds().SetId(2, n2)
                #elem.GetCellType() = 5  # vtkTriangle
                grid.InsertNextCell(5, elem.GetPointIds())
        elif dimension_flag == 3:
            for (n0, n1, n2, n3) in tets:
                elem = vtkTetra()
                elem.GetPointIds().SetId(0, n0)
                elem.GetPointIds().SetId(1, n1)
                elem.GetPointIds().SetId(2, n2)
                elem.GetPointIds().SetId(3, n3)
                #elem.GetCellType() = 10  # vtkTetra
                grid.InsertNextCell(10, elem.GetPointIds())
        else:
            raise RuntimeError('dimension_flag=%r' % dimension_flag)

        grid.SetPoints(points)
        grid.Modified()

        # regions/loads
        self.gui.scalar_bar_actor.Modified()

        cases = OrderedDict()
        #cases = self.result_cases
        form = []
        node_ids, element_ids = self._fill_fast_results(
            form, cases, model,
            nnodes, nelements, dimension_flag,
            results=True)
        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids
        self.gui._finish_results_io2(model_name, form, cases)

    def clear_fast(self):
        pass

    def _fill_fast_results(self, form, cases, model,
                           nnodes, nelements, dimension_flag,
                           results=False):
        note = ''
        self.gui.isubcase_name_map = {
            1: ['Fast%s' % note, ''],
            #2: ['Fast%s' % note, ''],
        }
        cases, node_ids, element_ids = self._fill_fast_case(
            form, cases, model,
            nnodes, nelements, dimension_flag,
            results=results)
        return node_ids, element_ids

    def _fill_fast_case(self, form, cases, model,
                        nnodes, nelements, dimension_flag,
                        results=False):
        self.gui.scalar_bar_actor.VisibilityOff()

        icase = 0
        geometry_form = [
            ('NodeID', icase + 0, []),
            ('ElementID', icase + 1, []),
        ]
        res_id = 1
        nids = np.arange(1, nnodes + 1)
        eids = np.arange(1, nelements + 1)

        nid_res = GuiResult(res_id, header='NodeID', title='NodeID',
                            location='node', scalar=nids)
        cases[icase] = (nid_res, (0, 'NodeID'))
        icase += 1

        eid_res = GuiResult(res_id, header='ElementID', title='ElementID',
                            location='centroid', scalar=eids)
        cases[icase] = (eid_res, (0, 'ElementID'))
        icase += 1

        if dimension_flag == 2:
            geometry_form.append(('BC', icase + 2, []))
            bc_res = GuiResult(res_id, header='BoundaryCondition', title='BC',
                               location='centroid', scalar=model.bcs)
            cases[icase] = (bc_res, (0, 'BoundaryCondition'))
            icase += 1

        #if results:
            #res_id = 1
            #if bcs is not None:
                #cases[(res_id, icase, 'Region', 1, 'centroid', '%i', '')] = bcs
                #icase += 1

                #mapbc_print = defaultdict(list)
                #for region, bcnum in sorted(mapbc.items()):
                    #mapbc_print[bcnum].append(region)
                    #try:
                        #name = bcmap_to_bc_name[bcnum]
                    #except KeyError:
                        #name = '???'
                    ##self.log.info('Region=%i BC=%s name=%r' % (region, bcnum, name))

                #for bcnum, regions in sorted(mapbc_print.items()):
                    #try:
                        #name = bcmap_to_bc_name[bcnum]
                    #except KeyError:
                        #name = '???'
                    #self.log.info('BC=%s Regions=%s name=%r' % (bcnum, regions, name))
                #self.scalar_bar_actor.VisibilityOn()

            ##==============================
            #res_id = 2
            #if len(loads):
                #for key, load in loads.items():
                    #cases[(res_id, icase, key, 1, 'node', '%.3f')] = load
                    #icase += 1
                #self.scalar_bar_actor.VisibilityOn()

        form.append(('Geometry', None, geometry_form))
        return cases, nids, eids
