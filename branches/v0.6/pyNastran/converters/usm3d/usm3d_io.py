import os
from collections import defaultdict

import vtk
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.usm3d.usm3d_reader import Usm3dReader


class Usm3dIO(object):
    def __init__(self):
        pass

    def load_usm3d_results(self, flo_filename, dirname):
        model = Usm3dReader(log=self.log, debug=False)
        #self.resultCases = {}
        npoints = self.nElements
        node_ids_volume, loads = model.read_flo(flo_filename, n=npoints)

        cases = self.resultCases
        bcs = None
        mapbc = None
        bcmap_to_bc_name = None
        self._fill_usm3d_results(cases, bcs, mapbc, bcmap_to_bc_name, loads)


    def load_usm3d_geometry(self, cogsg_filename, dirname):
        print "load_usm3d_geometry..."
        skipReading = self.removeOldGeometry(cogsg_filename)
        if skipReading:
            return

        model = Usm3dReader(log=self.log, debug=False)

        base_filename, ext = os.path.splitext(cogsg_filename)
        #node_filename = base_filename + '.node'
        #ele_filename = base_filename + '.ele'
        if '.cogsg' == ext:
            dimension_flag = 3
        #elif '.ele' == ext:
            #dimension_flag = 3
        else:
            raise RuntimeError('unsupported extension.  Use "cogsg" or "front".')
        model.read_usm3d(base_filename, dimension_flag)
        nodes = model.nodes
        tris = model.tris
        tets = model.tets
        bcs = model.bcs
        mapbc = model.mapbc
        loads = model.loads
        bcmap_to_bc_name = model.bcmap_to_bc_name

        self.nNodes, three = nodes.shape
        ntris = 0
        ntets = 0
        if tris is not None:
            ntris, three = tris.shape

        if dimension_flag == 2:
            pass
        elif dimension_flag == 3:
            ntets, four = tets.shape
            ntets = 0
        else:
            raise RuntimeError()
        self.nElements = ntris + ntets

        print("nNodes = ",self.nNodes)
        print("nElements = ", self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)
        self.grid2.Allocate(1, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(nodes.iteritems()):
                points.InsertPoint(nid - 1, *node)
                self.gridResult.InsertNextValue(nid * fraction)
                #print str(element)

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes, three = nodes.shape

        nid = 0
        print "nnodes=%s" % nnodes
        for i in xrange(nnodes):
            points.InsertPoint(nid, nodes[i, :])
            nid += 1

        #elements -= 1
        if ntris:
            for (n0, n1, n2) in tris:
                elem = vtkTriangle()
                #node_ids = elements[eid, :]
                elem.GetPointIds().SetId(0, n0)
                elem.GetPointIds().SetId(1, n1)
                elem.GetPointIds().SetId(2, n2)
                self.grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle

        if dimension_flag == 2:
            pass
        elif dimension_flag == 3:
            if ntets:
                for (n0, n1, n2, n3) in tets:
                    elem = vtkTetra()
                    #assert elem.GetCellType() == 10, elem.GetCellType()
                    elem.GetPointIds().SetId(0, n0)
                    elem.GetPointIds().SetId(1, n1)
                    elem.GetPointIds().SetId(2, n2)
                    elem.GetPointIds().SetId(3, n3)
                    self.grid.InsertNextCell(10, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        else:
            raise RuntimeError()

        self.grid.SetPoints(points)
        self.grid.Modified()
        self.grid.Update()
        print("updated grid")

        # regions/loads
        self.TurnTextOn()
        self.scalarBar.Modified()

        cases = {}
        #cases = self.resultCases
        self._fill_usm3d_results(cases, bcs, mapbc, bcmap_to_bc_name, loads)

    def _fill_usm3d_results(self, cases, bcs, mapbc, bcmap_to_bc_name, loads):
        if 'Mach' in loads:
            avgMach = loads['Mach'].mean()
            note = ':  avg(Mach)=%g' % avgMach
        else:
            note = ''

        self.iSubcaseNameMap = {
            1: ['Usm3d%s' % note, ''],
            2: ['Usm3d%s' % note, ''],
        }

        #ID = 1
        cases = self._fill_usm3d_case(cases, bcs, mapbc, bcmap_to_bc_name, loads)

        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print("caseKeys = ",self.caseKeys)
        #print "type(caseKeys) = ",type(self.caseKeys)
        if len(self.resultCases) == 0:
            self.nCases = 1
        elif len(self.resultCases) == 1:
            self.nCases = 1
        else:
            self.nCases = len(self.resultCases) - 1  # number of keys in dictionary

        self.iCase = 0 if self.nCases == 0 else -1
        self.cycleResults()  # start at nCase=0

    def _fill_usm3d_case(self, cases, bcs, mapbc, bcmap_to_bc_name, loads):
        self.scalarBar.VisibilityOff()

        ID = 1
        if bcs is not None and self.is_centroidal:
            cases[(ID, 'Region', 1, 'centroid', '%.0f')] = bcs

            mapbc_print = defaultdict(list)
            for region, bcnum in sorted(mapbc.iteritems()):
                mapbc_print[bcnum].append(region)
                try:
                    name = bcmap_to_bc_name[bcnum]
                except KeyError:
                    name = '???'
                #self.log.info('Region=%i BC=%s name=%r' % (region, bcnum, name))

            for bcnum, regions in sorted(mapbc_print.iteritems()):
                try:
                    name = bcmap_to_bc_name[bcnum]
                except KeyError:
                    name = '???'
                self.log.info('BC=%s Regions=%s name=%r' % (bcnum, regions, name))
            self.scalarBar.VisibilityOn()

        #==============================
        ID = 2
        if self.is_nodal and len(loads):
            for key, load in loads.iteritems():
                cases[(ID, key, 1, 'nodal', '%.3f')] = load
            self.scalarBar.VisibilityOn()
        return cases