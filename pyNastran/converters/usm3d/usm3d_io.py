from six import iteritems
from six.moves import range
import os
from collections import defaultdict

import vtk
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.usm3d.usm3d_reader import Usm3d
from pyNastran.converters.usm3d.time_accurate_results import get_n_list


class Usm3dIO(object):
    def __init__(self):
        pass

    def get_usm3d_wildcard_geometry_results_functions(self):
        data = ('Usm3D',
                'USM3D (*.cogsg; *.front)', self.load_usm3d_geometry,
                'Usm3d (*.flo)', self.load_usm3d_results)
        return data

    def step_results_usm3d(self):
        # minimum is 1
        nstep = 100

        assert self.out_filename is not None, self.out_filename
        flo_filename = self.out_filename
        dirname = os.path.dirname(flo_filename)
        if dirname == '':
            dirname = os.getcwd()
        basename = os.path.basename(flo_filename)
        base = os.path.splitext(basename)[0]

        if '_' in base:
            model_name, n = base.rsplit('_', 1)
            #print("model_name=%r n=%r" % (model_name, n))
            n = int(n)
            n_list = get_n_list(dirname, model_name)
            inn = n_list.index(n)
            if inn+nstep < len(n_list):
                nnew = n_list[inn+nstep]
            else:
                nnew = max(n_list)
                if nnew == n:
                    raise RuntimeError('%r is the last file' % self.out_filename)
            #print("inn=%r nnew=%r" % (inn, nnew))
            flo_filename = model_name + '_%s.flo' % nnew
        else:
            raise RuntimeError('The current file is must have the format of xxx_%%i.flo, not %r' % self.out_filename)
        #print("loading %r" % flo_filename)
        self.load_usm3d_results(flo_filename, dirname)
        self.out_filename = os.path.join(dirname, flo_filename)

        print("done stepping...")

    #def _get_next_n(self, base):
        #n = int(n)
        ## get the max N value
        #nmax = -1
        #for flo_filename in flo_filenames:
            #base, ext = os.path.splitext(flo_filename)
            #if ext == '.flo':
                #n = base.split('_')[-1]
                #try: # get the incrementation index
                    #n = int(n)
                    #if n > nold:
                        #return n
                #except:
                    #raise NotImplementedError()
        #return None

    def load_usm3d_results(self, flo_filename, dirname):
        model = Usm3d(log=self.log, debug=False)
        #self.result_cases = {}
        npoints = self.nNodes
        node_ids_volume, loads = model.read_flo(flo_filename, n=npoints)

        cases = self.result_cases
        bcs = None
        mapbc = None
        bcmap_to_bc_name = None
        self._fill_usm3d_results(cases, bcs, mapbc, bcmap_to_bc_name, loads)

    def load_usm3d_geometry(self, cogsg_filename, dirname, name='main', plot=True):
        #print("load_usm3d_geometry...")
        skip_reading = self._remove_old_geometry(cogsg_filename)
        if skip_reading:
            return

        model = Usm3d(log=self.log, debug=False)

        base_filename, ext = os.path.splitext(cogsg_filename)
        #node_filename = base_filename + '.node'
        #ele_filename = base_filename + '.ele'
        if '.cogsg' == ext:
            dimension_flag = 3
        #elif '.ele' == ext:
            #dimension_flag = 3
        else:
            raise RuntimeError('unsupported extension.  Use "cogsg" or "front".')

        read_loads = True
        nodes, tris_tets, tris, bcs, mapbc, loads, flo_filename = model.read_usm3d(base_filename, dimension_flag, read_loads=read_loads)
        del tris_tets
        nodes = model.nodes
        tris = model.tris
        tets = model.tets
        bcs = model.bcs
        mapbc = model.mapbc
        loads = model.loads

        self.out_filename = None
        if flo_filename is not None:
            self.out_filename = flo_filename

        bcmap_to_bc_name = model.bcmap_to_bc_name

        self.nNodes = nodes.shape[0]
        ntris = 0
        ntets = 0
        if tris is not None:
            ntris = tris.shape[0]

        if dimension_flag == 2:
            pass
        elif dimension_flag == 3:
            ntets = tets.shape[0]
            ntets = 0
        else:
            raise RuntimeError()
        self.nElements = ntris + ntets

        print("nNodes = %i" % self.nNodes)
        print("nElements = %i" % self.nElements)

        self.grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        #self.gridResult.Allocate(self.nNodes, 1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nid_map = {}
        #elem.SetNumberOfPoints(nNodes)
        if 0:
            fraction = 1. / self.nNodes  # so you can color the nodes by ID
            for nid, node in sorted(iteritems(nodes)):
                points.InsertPoint(nid - 1, *node)
                self.gridResult.InsertNextValue(nid * fraction)
                #print(str(element))

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes = nodes.shape[0]

        nid = 0
        #print("nnodes=%s" % nnodes)
        for i in range(nnodes):
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
            raise RuntimeError('dimension_flag=%r' % dimension_flag)

        self.grid.SetPoints(points)
        self.grid.Modified()
        if hasattr(self.grid, 'Update'):
            self.grid.Update()

        # regions/loads
        self. turn_text_on()
        self.scalarBar.Modified()

        cases = {}
        #cases = self.result_cases
        self._fill_usm3d_results(cases, bcs, mapbc, bcmap_to_bc_name, loads)
        self._finish_results_io(cases)

    def clear_usm3d(self):
        pass

    def _fill_usm3d_results(self, cases, bcs, mapbc, bcmap_to_bc_name, loads):
        if 'Mach' in loads:
            avg_mach = loads['Mach'].mean()
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''

        self.iSubcaseNameMap = {
            1: ['Usm3d%s' % note, ''],
            2: ['Usm3d%s' % note, ''],
        }

        #ID = 1
        cases = self._fill_usm3d_case(cases, bcs, mapbc, bcmap_to_bc_name, loads)

    def _fill_usm3d_case(self, cases, bcs, mapbc, bcmap_to_bc_name, loads):
        self.scalarBar.VisibilityOff()

        ID = 1
        if bcs is not None:
            cases[(ID, 'Region', 1, 'centroid', '%i')] = bcs

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

        ID = 2
        if len(loads):
            for key, load in iteritems(loads):
                cases[(ID, key, 1, 'node', '%.3f', '')] = load
            self.scalarBar.VisibilityOn()
        return cases
