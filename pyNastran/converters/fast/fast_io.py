"""
Defines the GUI IO file for Fast.
"""
from __future__ import print_function
from six import iteritems
from six.moves import range
import os
from collections import defaultdict

import vtk
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.fast.fgrid_reader import FGridReader


class FastIO(object):
    def __init__(self):
        pass

    def get_fast_wildcard_geometry_results_functions(self):
        data = ('FAST',
                'FAST (*.fgrid)', self.load_fast_geometry,
                None, None)
        return data

    def load_fast_geometry(self, fgrid_filename, dirname, name='main', plot=True):
        skip_reading = self._remove_old_geometry(fgrid_filename)
        if skip_reading:
            return

        model = FGridReader(log=self.log, debug=False)

        base_filename, ext = os.path.splitext(fgrid_filename)
        if '.fgrid' == ext:
            dimension_flag = 3
        #elif '.ele' == ext:
            #dimension_flag = 3
        else:
            raise RuntimeError('unsupported extension.  Use "cogsg" or "front".')

        read_loads = True
        model.read_fgrid(fgrid_filename, dimension_flag)

        dimension_flag = 3
        nodes = model.nodes
        tris = model.tris - 1
        tets = model.tets - 1

        nnodes = nodes.shape[0]
        ntris = tris.shape[0]
        ntets = tets.shape[0]

        #print('node0 = %s' % str(nodes[0, :]))
        #print('node%i = %s' % (1, str(nodes[1, :])))
        #print('node%i = %s' % (2, str(nodes[2, :])))
        #print('node%i = %s' % (nnodes, str(nodes[-1, :])))
        #print('tris.max/min = ', tris.max(), tris.min())
        #print('tets.max/min = ', tets.max(), tets.min())
        #bcs = model.bcs
        #mapbc = model.mapbc
        #loads = model.loads

        self.nNodes = nnodes
        self.nElements = ntris + ntets

        #print("nNodes = %i" % self.nNodes)
        #print("nElements = %i" % self.nElements)

        grid = self.grid
        grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = self.numpy_to_vtk_points(nodes)
        self.nid_map = {}
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                self.gridResult.InsertNextValue(nid * fraction)

        assert nodes is not None
        nnodes = nodes.shape[0]

        if dimension_flag == 2:
            for (n0, n1, n2) in tris:
                elem = vtkTriangle()
                #node_ids = elements[eid, :]
                elem.GetPointIds().SetId(0, n0)
                elem.GetPointIds().SetId(1, n1)
                elem.GetPointIds().SetId(2, n2)
                grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        elif dimension_flag == 3:
            if ntets:
                for (n0, n1, n2, n3) in tets:
                    elem = vtkTetra()
                    elem.GetPointIds().SetId(0, n0)
                    elem.GetPointIds().SetId(1, n1)
                    elem.GetPointIds().SetId(2, n2)
                    elem.GetPointIds().SetId(3, n3)
                    grid.InsertNextCell(10, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        else:
            raise RuntimeError('dimension_flag=%r' % dimension_flag)

        grid.SetPoints(points)
        grid.Modified()
        if hasattr(grid, 'Update'):
            grid.Update()

        # regions/loads
        self.scalarBar.Modified()

        cases = {}
        #cases = self.result_cases
        self._fill_fast_results(cases, model, results=False)
        self._finish_results_io(cases)

    def clear_fast(self):
        pass

    def _fill_fast_results(self, cases, model, results=False):
        note = ''
        if 0:
            if 'Mach' in loads:
                avg_mach = loads['Mach'].mean()
                note = ':  avg(Mach)=%g' % avg_mach
            else:
                note = ''

        self.iSubcaseNameMap = {
            1: ['Fast%s' % note, ''],
            2: ['Fast%s' % note, ''],
        }

        #ID = 1
        cases = self._fill_fast_case(cases, model, results=results)

    def _fill_fast_case(self, cases, model, results=False):
        self.scalarBar.VisibilityOff()

        icase = 0
        geometry_form = [
            ('ElementID', icase, [])
        ]
        if results:
            ID = 1
            if bcs is not None:
                cases[(ID, icase, 'Region', 1, 'centroid', '%i', '')] = bcs
                icase += 1

                mapbc_print = defaultdict(list)
                for region, bcnum in sorted(iteritems(mapbc)):
                    mapbc_print[bcnum].append(region)
                    try:
                        name = bcmap_to_bc_name[bcnum]
                    except KeyError:
                        name = '???'
                    #self.log.info('Region=%i BC=%s name=%r' % (region, bcnum, name))

                for bcnum, regions in sorted(iteritems(mapbc_print)):
                    try:
                        name = bcmap_to_bc_name[bcnum]
                    except KeyError:
                        name = '???'
                    self.log.info('BC=%s Regions=%s name=%r' % (bcnum, regions, name))
                self.scalarBar.VisibilityOn()

            #==============================
            ID = 2
            if len(loads):
                for key, load in iteritems(loads):
                    cases[(ID, icase, key, 1, 'node', '%.3f')] = load
                    icase += 1
                self.scalarBar.VisibilityOn()
        return cases
