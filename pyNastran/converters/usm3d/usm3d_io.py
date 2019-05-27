"""Defines the GUI IO file for Usm3d."""
import os
from collections import defaultdict, OrderedDict

import numpy as np

from pyNastran.utils import object_attributes
from pyNastran.utils.numpy_utils import integer_float_types
from pyNastran.converters.usm3d.usm3d_reader import Usm3d
from pyNastran.converters.usm3d.time_accurate_results import get_n_list

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type, numpy_to_vtk_points)


class Usm3dIO:
    def __repr__(self):
        return '<Usm3dIO class>'

    def __init__(self, gui):
        self.gui = gui
        assert gui is not None

    def get_usm3d_wildcard_geometry_results_functions(self):
        data = ('Usm3D',
                'USM3D (*.cogsg; *.front)', self.load_usm3d_geometry,
                'Usm3d (*.flo)', self.load_usm3d_results)
        return data

    def on_reload_usm3d(self):
        """
        For USM3D, we dynamically load the latest CFD results time step,
        hich is really handy when you're running a job.
        """
        # minimum is 1
        nstep = 100

        if self.gui.out_filename is None:
            msg = 'usm3d_filename=%r must not be None\n' % self.gui.out_filename
            dir_gui = []
            for key in object_attributes(self.gui):
                try:
                    value = getattr(self.gui, key)
                except KeyError:
                    # self.edge_actor is a
                    if key not in ['edge_actor']:
                        self.gui.log.warning('key=%s is undefined...' % key)

                if isinstance(value, (integer_float_types, str)):
                    dir_gui.append(key)
            dir_gui.sort()
            msg += 'dir(gui) = [%s]' % ', '.join(dir_gui)
            raise RuntimeError(msg)
        flo_filename = self.gui.out_filename
        dirname = os.path.dirname(flo_filename)
        if dirname == '':
            dirname = os.getcwd()
        basename = os.path.basename(flo_filename)
        base = os.path.splitext(basename)[0]


        # box.flo -> box_100.flo
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
                    raise RuntimeError('%r is the last file' % self.gui.out_filename)
            #print("inn=%r nnew=%r" % (inn, nnew))
            flo_filename = model_name + '_%s.flo' % nnew
        else:
            flo_filename = self.gui.out_filename
            #msg = (
                #'The current file is must have the format of '
                #'xxx_%%i.flo, not %r' % self.out_filename)
            #raise RuntimeError(msg)
        #print("loading %r" % flo_filename)
        self.load_usm3d_results(flo_filename)
        self.gui.out_filename = os.path.join(dirname, flo_filename)

        #print("done stepping...")

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

    def load_usm3d_results(self, flo_filename):
        model = Usm3d(log=self.gui.log, debug=False)
        npoints = self.gui.nnodes
        unused_node_ids_volume, loads = model.read_flo(flo_filename, n=npoints)

        cases = self.gui.result_cases
        form = self.gui.get_form()
        bcs = None
        mapbc = None
        bcmap_to_bc_name = None

        self._fill_usm3d_results(cases, form,
                                 bcs, mapbc, bcmap_to_bc_name, loads,
                                 is_geometry=False)

    def load_usm3d_geometry(self, cogsg_filename, name='main', plot=True):
        model_name = name
        skip_reading = self.gui._remove_old_geometry(cogsg_filename)
        if skip_reading:
            return

        self.gui.eid_maps[name] = {}
        self.gui.nid_maps[name] = {}
        model = Usm3d(log=self.gui.log, debug=False)

        base_filename, ext = os.path.splitext(cogsg_filename)
        #node_filename = base_filename + '.node'
        #ele_filename = base_filename + '.ele'
        if ext == '.cogsg':
            dimension_flag = 3
        #elif ext == '.ele':
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

        self.gui.out_filename = None
        if flo_filename is not None:
            self.gui.out_filename = flo_filename

        bcmap_to_bc_name = model.bcmap_to_bc_name

        self.gui.nnodes = nodes.shape[0]
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
        self.gui.nelements = ntris + ntets

        self.gui.log.debug("nnodes = %i" % self.gui.nnodes)
        self.gui.log.debug("nelements = %i" % self.gui.nelements)

        grid = self.gui.grid
        grid.Allocate(self.gui.nelements, 1000)

        self.gui.nid_map = {}
        self.gui.eid_map = {}

        assert nodes is not None
        nnodes = nodes.shape[0]
        node_ids = np.arange(1, nnodes + 1, dtype='int32')

        points = numpy_to_vtk_points(nodes)
        if ntris:
            element_ids = np.arange(1, ntris + 1, dtype='int32')
            etype = 5  # vtkTriangle().GetCellType()
            create_vtk_cells_of_constant_element_type(grid, tris, etype)
        else:
            ntets = tets.shape[0]
            element_ids = np.arange(1, ntets + 1, dtype='int32')

        if dimension_flag == 2:
            pass
        elif dimension_flag == 3:
            if ntets:
                etype = 10 # vtkTetra().GetCellType()
                assert tets.max() > 0, tets.min()
                create_vtk_cells_of_constant_element_type(grid, tets, etype)
        else:
            raise RuntimeError('dimension_flag=%r' % dimension_flag)

        grid.SetPoints(points)
        grid.Modified()

        self.gui.node_ids = node_ids
        self.gui.element_ids = element_ids

        # regions/loads
        self.gui.scalar_bar_actor.Modified()

        cases = OrderedDict()
        form = []
        form, cases = self._fill_usm3d_results(cases, form,
                                               bcs, mapbc, bcmap_to_bc_name, loads,
                                               is_geometry=True)
        self.gui._finish_results_io2(model_name, form, cases)

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

        self.gui.isubcase_name_map = {
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
        self.gui.scalar_bar_actor.VisibilityOff()
        colormap = self.gui.settings.colormap

        subcasemap_id = 1
        icase = len(cases)
        itime = 0
        if is_geometry:
            assert self.gui.element_ids is not None, self.gui.element_ids
            assert len(self.gui.element_ids) > 0, self.gui.element_ids
            eid_res = GuiResult(
                subcasemap_id, 'ElementID', 'ElementID', 'centroid', self.gui.element_ids,
                nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
                data_format='%i', uname='GuiResult')

            cases[icase] = (eid_res, (itime, 'ElementID'))
            form.append(('ElementID', icase, []))
            icase += 1
            if bcs is not None:
                region_res = GuiResult(
                    subcasemap_id, 'Patch', 'Patch', 'centroid', bcs,  # patch_id
                    nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
                    data_format='%i', uname='GuiResult')
                cases[icase] = (region_res, (itime, 'Patch'))
                form.append(('Patch', icase, []))
                icase += 1

        if bcs is not None:
            patch_id = bcs

            form += [
                ('BC', icase, []),
                ('Family', icase + 1, []),
            ]
            bc_value = np.zeros(bcs.shape, dtype='int32')
            family = np.zeros(bcs.shape, dtype='int32')
            mapbc_print = defaultdict(list)
            for region, mapi in sorted(mapbc.items()):
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
                               nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
                               data_format='%i', uname='GuiResult')
            family_res = GuiResult(subcasemap_id, 'Family', 'Family', 'centroid', family,
                                   nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
                                   data_format='%i', uname='GuiResult')
            cases[icase] = (bc_res, (itime, 'BC'))
            cases[icase + 1] = (family_res, (itime, 'Family'))
            icase += 2


            for bcnum, regions in sorted(mapbc_print.items()):
                try:
                    name = bcmap_to_bc_name[bcnum]
                except KeyError:
                    name = '???'
                self.gui.log.info('BC=%s Regions=%s name=%r' % (bcnum, regions, name))

            self.gui.scalar_bar_actor.VisibilityOn()

        subcasemap_id = 2
        if len(loads):
            form0 = []
            for key, load in loads.items():
                load_res = GuiResult(subcasemap_id, key, key, 'node', load,
                                     nlabels=None, labelsize=None, ncolors=None, colormap=colormap,
                                     data_format='%.3f', uname='GuiResult')
                cases[icase] = (load_res, (itime, key))
                formi = (key, icase, [])
                form0.append(formi)
                icase += 1

            if form0:
                form.append(('Results', None, form0))
        self.gui.scalar_bar_actor.VisibilityOn()
        return form, cases
