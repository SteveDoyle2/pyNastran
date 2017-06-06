from six import iteritems
from six.moves import range
import os
from collections import defaultdict

import numpy as np
import vtk
from vtk import vtkTriangle, vtkTetra

from pyNastran.converters.usm3d.usm3d_reader import Usm3d
from pyNastran.converters.usm3d.time_accurate_results import get_n_list

from pyNastran.gui.gui_objects.gui_result import GuiResult


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
            msg = 'The current file is must have the format of xxx_%%i.flo, not %r' % self.out_filename
            raise RuntimeError(msg)
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
        form = self.get_form()
        bcs = None
        mapbc = None
        bcmap_to_bc_name = None

        self._fill_usm3d_results(cases, form,
                                 bcs, mapbc, bcmap_to_bc_name, loads,
                                 is_geometry=False)

    def load_usm3d_geometry(self, cogsg_filename, dirname, name='main', plot=True):
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
        nodes, tris_tets, tris, bcs, mapbc, loads, flo_filename = model.read_usm3d(
            base_filename, dimension_flag, read_loads=read_loads)
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

        self.log.debug("nNodes = %i" % self.nNodes)
        self.log.debug("nElements = %i" % self.nElements)

        grid = self.grid
        grid.Allocate(self.nElements, 1000)
        #self.gridResult.SetNumberOfComponents(self.nElements)

        self.nid_map = {}
        self.eid_map = {}
        #elem.SetNumberOfPoints(nNodes)
        #if 0:
            #fraction = 1. / self.nNodes  # so you can color the nodes by ID
            #for nid, node in sorted(iteritems(nodes)):
                #points.InsertPoint(nid - 1, *node)
                #self.gridResult.InsertNextValue(nid * fraction)
                #print(str(element))

                #elem = vtk.vtkVertex()
                #elem.GetPointIds().SetId(0, i)
                #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

        assert nodes is not None
        nnodes = nodes.shape[0]

        points = self.numpy_to_vtk_points(nodes)
        if ntris:
            self.element_ids = np.arange(1, ntris + 1, dtype='int32')
            etype = 5  # vtkTriangle().GetCellType()
            self.create_vtk_cells_of_constant_element_type(grid, tris, etype)
            #for (n0, n1, n2) in tris:
                #elem = vtkTriangle()
                #elem.GetPointIds().SetId(0, n0)
                #elem.GetPointIds().SetId(1, n1)
                #elem.GetPointIds().SetId(2, n2)
                #grid.InsertNextCell(5, elem.GetPointIds())  #elem.GetCellType() = 5  # vtkTriangle
        else:
            ntets = tets.shape[0]
            self.element_ids = np.arange(1, ntets + 1, dtype='int32')

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
        form = []
        form, cases = self._fill_usm3d_results(cases, form,
                                               bcs, mapbc, bcmap_to_bc_name, loads,
                                               is_geometry=True)
        self._finish_results_io2(form, cases)

    def clear_usm3d(self):
        """dummy function"""
        pass

    def _fill_usm3d_results(self, cases, form,
                            bcs, mapbc, bcmap_to_bc_name, loads,
                            is_geometry=True):
        """sets up usm3d results"""
        if 'Mach' in loads:
            avg_mach = loads['Mach'].mean()
            note = ':  avg(Mach)=%g' % avg_mach
        else:
            note = ''

        self.iSubcaseNameMap = {
            1: ['Usm3d%s' % note, ''],
            2: ['Usm3d%s' % note, ''],
        }

        form, cases = self._fill_usm3d_case(
            cases, form,
            bcs, mapbc, bcmap_to_bc_name, loads,
            is_geometry=is_geometry)
        return form, cases

    def _fill_usm3d_case(self, cases, form,
                         bcs, mapbc, bcmap_to_bc_name, loads, is_geometry=True):
        """actually fills the sidebar"""
        self.scalarBar.VisibilityOff()

        subcasemap_id = 1
        icase = len(cases)
        itime = 0
        if is_geometry:
            assert self.element_ids is not None, self.element_ids
            assert len(self.element_ids) > 0, self.element_ids
            eid_res = GuiResult(
                subcasemap_id, 'ElementID', 'ElementID', 'centroid', self.element_ids,
                nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                data_format='%i', uname='GuiResult')

            region_res = GuiResult(
                subcasemap_id, 'Patch', 'Patch', 'centroid', bcs,  # patch_id
                nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                data_format='%i', uname='GuiResult')
            cases[icase] = (eid_res, (itime, 'ElementID'))
            cases[icase + 1] = (region_res, (itime, 'Patch'))
            form.append(('ElementID', icase, []))
            form.append(('Patch', icase + 1, []))
            icase += 2

        if bcs is not None:
            patch_id = bcs

            form += [
                ('BC', icase, []),
                ('Family', icase + 1, []),
            ]
            bc_value = np.zeros(bcs.shape, dtype='int32')
            family = np.zeros(bcs.shape, dtype='int32')
            mapbc_print = defaultdict(list)
            for region, mapi in sorted(iteritems(mapbc)):
                bcnum = mapi[0]
                familyi = mapi[1]
                mapbc_print[bcnum].append(region)
                try:
                    name = bcmap_to_bc_name[bcnum]
                except KeyError:
                    name = '???'
                #self.log.info('Region=%i BC=%s name=%r' % (region, bcnum, name))
                ipatch = np.where(patch_id == region)[0]
                bc_value[ipatch] = bcnum
                family[ipatch] = familyi

            bc_res = GuiResult(subcasemap_id, 'BC', 'BC', 'centroid', bc_value,
                               nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                               data_format='%i', uname='GuiResult')
            family_res = GuiResult(subcasemap_id, 'Family', 'Family', 'centroid', family,
                                   nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                                   data_format='%i', uname='GuiResult')
            cases[icase] = (bc_res, (itime, 'BC'))
            cases[icase + 1] = (family_res, (itime, 'Family'))
            icase += 2


            for bcnum, regions in sorted(iteritems(mapbc_print)):
                try:
                    name = bcmap_to_bc_name[bcnum]
                except KeyError:
                    name = '???'
                self.log.info('BC=%s Regions=%s name=%r' % (bcnum, regions, name))

            self.scalarBar.VisibilityOn()

        subcasemap_id = 2
        if len(loads):
            form0 = []
            for key, load in iteritems(loads):
                load_res = GuiResult(subcasemap_id, key, key, 'node', load,
                                     nlabels=None, labelsize=None, ncolors=None, colormap='jet',
                                     data_format='%.3f', uname='GuiResult')
                cases[icase] = (load_res, (itime, key))
                formi = (key, icase, [])
                form0.append(formi)
                icase += 1

            if form0:
                form.append(('Results', None, form0))
        self.scalarBar.VisibilityOn()
        return form, cases
