import os
import sys
import traceback
import time as time_module
from typing import Optional, Any, TYPE_CHECKING

import numpy as np
from qtpy.compat import getopenfilename
#from qtpy.QtWidgets import QFileDialog
from pyNastran.bdf.patran_utils.read_patran_custom_results import load_patran_nod
from pyNastran.utils import print_bad_path

from pyNastran.gui.qt_files.base_gui import BaseGui
from pyNastran.gui.utils.load_results import (
    load_csv, load_deflection_csv, create_res_obj,
    load_user_geom,
)
from pyNastran.gui.qt_files.user_geometry import add_user_geometry

if TYPE_CHECKING:
    #from pyNastran.gui.menus.results_sidebar import Sidebar
    from pyNastran.gui.gui_objects.settings import Settings
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.qt_files.tool_actions import ToolActions
IS_TESTING = 'test' in sys.argv[0]


class LoadActions(BaseGui):
    """performance mode should be handled in the main gui to minimize flipping"""

    #@property
    #def log(self):
        #"""links the the GUI's log"""
        #return self.gui.log

    def on_load_geometry(self, infile_name=None, geometry_format=None,
                         name: str='main',
                         plot: bool=True,
                         stop_on_failure: bool=False):
        """
        Loads a baseline geometry

        Parameters
        ----------
        infile_name : str; default=None -> popup
            path to the filename
        geometry_format : str; default=None
            the geometry format for programmatic loading
        name : str; default='main'
            the name of the actor; don't use this
        plot : bool; default=True
            Should the baseline geometry have results created and plotted/rendered?
            If you're calling the on_load_results method immediately after, set it to False
        stop_on_failure : bool; default=True
            stop the code if True

        """
        assert isinstance(name, str), 'name=%r type=%s' % (name, type(name))
        is_failed, out = self._load_geometry_filename(
            geometry_format, infile_name)
        print("is_failed =", is_failed)
        if is_failed:
            return

        if not hasattr(self, 'model_objs'):
            self.model_objs = {}

        gui = self.gui
        has_results = False
        infile_name, load_function, filter_index, formats, geometry_format2 = out
        if load_function is not None:
            #settings: Settings = gui.settings
            #if settings.use_startup_directory and os.path.exists(settings.startup_directory):
                #last_dir = settings.startup_directory
            #else:
            self._set_last_dir(infile_name)

            if gui.name == '':
                name = 'main'
            else:
                print('name = %r' % name)

            if name != gui.name:
                #scalar_range = self.grid_selected.GetScalarRange()
                #self.grid_mapper.SetScalarRange(scalar_range)
                gui.grid_mapper.ScalarVisibilityOff()
                #self.grid_mapper.SetLookupTable(self.color_function)
            gui.name = name
            gui._reset_model(name)

            # reset alt grids
            names = gui.alt_grids.keys()
            for name in names:
                gui.alt_grids[name].Reset()
                gui.alt_grids[name].Modified()

            if not os.path.exists(infile_name) and geometry_format:
                msg = f'input file={infile_name!r} does not exist'
                gui.log_error(msg)
                gui.log_error(print_bad_path(infile_name))
                return

            # clear out old data
            if gui.model_type is not None:
                clear_name = 'clear_' + gui.model_type
                try:
                    dy_method = getattr(self, clear_name)  # 'self.clear_nastran()'
                    dy_method()
                except Exception:
                    self.gui.log_error("method %r does not exist" % clear_name)
            gui.log_info(f'reading {geometry_format} file {infile_name!r}')

            try:
                time0 = time_module.time()

                if geometry_format2 in gui.format_class_map:
                    # initialize the class
                    #print('geometry_format=%r geometry_format2=%s' % (geometry_format, geometry_format2))

                    # TODO: was geometry_format going into this...
                    cls = gui.format_class_map[geometry_format2](gui)

                    function_name2 = f'load_{geometry_format2}_geometry'
                    load_function2 = getattr(cls, function_name2)
                    has_results = load_function2(infile_name, name=name, plot=plot)
                    self.model_objs[name] = cls
                else:
                    #  nastran
                    assert 'nastran' in geometry_format2
                    has_results = load_function(infile_name, name=name, plot=plot) # self.last_dir,

                dt = time_module.time() - time0
                print('dt_load = %.2f sec = %.2f min' % (dt, dt / 60.))
                #else:
                    #name = load_function.__name__
                    #self.log_error(str(args))
                    #self.log_error("'plot' needs to be added to %r; "
                                   #"args[-1]=%r" % (name, args[-1]))
                    #has_results = load_function(infile_name) # , self.last_dir
                    #form, cases = load_function(infile_name) # , self.last_dir
            except Exception as error:
                #raise
                msg = traceback.format_exc()
                gui.log_error(msg)
                if stop_on_failure or gui.dev:
                    raise
                #return
            #self.vtk_panel.Update()
            gui.rend.ResetCamera()

        # the model has been loaded, so we enable load_results
        if filter_index >= 0:
            geometry_format_out = formats[filter_index].lower()
            gui.format = geometry_format_out
            unused_enable = has_results
            #self.load_results.Enable(enable)
        else: # no file specified
            return
        #print("on_load_geometry(infile_name=%r, geometry_format=None)" % infile_name)
        gui.infile_name = infile_name
        gui.out_filename = None
        #if self.out_filename is not None:
            #msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)

        if name == 'main':
            msg = '%s - %s' % (geometry_format_out, gui.infile_name)
            gui.window_title = msg
            gui.update_menu_bar()
            main_str = ''
        else:
            main_str = ', name=%r' % name
        self._add_recent_file(infile_name, geometry_format)

        gui.settings.recent_files
        gui.log_command("on_load_geometry(infile_name=%r, geometry_format=%r%s)" % (
            infile_name, geometry_format_out, main_str))

    def _set_last_dir(self, infile_name: str) -> None:
        """
        Update the directory after finding a file not after completing
        sucessfully.  If the user explicitly loaded a file as part of
        the script, that doesn't count.
        """
        gui: MainWindow = self.gui
        last_dir = os.path.split(infile_name)[0]
        gui.last_dir = last_dir
        gui.settings.startup_directory = last_dir

    def _add_recent_file(self, infile_name: str,
                         geometry_format: str) -> None:
        gui: MainWindow = self.gui
        arg = (infile_name, geometry_format)
        if arg in gui.settings.recent_files:
            gui.settings.recent_files.remove(arg)
        gui.settings.recent_files.insert(0, arg)

    def _load_geometry_filename(self, geometry_format: str, infile_name: str) -> tuple[bool, Any]:
        """gets the filename and format"""
        wildcard = ''
        is_failed = False

        if geometry_format and geometry_format.lower() not in self.gui.supported_formats:
            is_failed = True
            msg = f'The import for the {geometry_format!r} module failed.\n'
            self.gui.log_error(msg)
            if IS_TESTING:  # pragma: no cover
                raise RuntimeError(msg)
            return is_failed, None

        if infile_name:
            if geometry_format is None:
                is_failed = True
                msg = 'infile_name=%r and geometry_format=%r; both must be specified\n' % (
                    infile_name, geometry_format)
                self.gui.log_error(msg)
                return is_failed, None

            geometry_format = geometry_format.lower()

            for fmt in self.gui.fmts:
                fmt_name, _major_name, _geom_wildcard, geom_func, res_wildcard, _resfunc = fmt
                if geometry_format == fmt_name:
                    load_function = geom_func
                    if res_wildcard is None:
                        unused_has_results = False
                    else:
                        unused_has_results = True
                    break
            else:
                self.gui.log_error('---invalid format=%r' % geometry_format)
                is_failed = True
                return is_failed, None
            formats = [geometry_format]
            filter_index = 0
        else:
            # load a pyqt window
            formats = []
            load_functions = []
            has_results_list = []
            wildcard_list = []

            # setup the selectable formats
            for fmt in self.gui.fmts:
                fmt_name, _major_name, geom_wildcard, geom_func, res_wildcard, _res_func = fmt
                formats.append(_major_name)
                wildcard_list.append(geom_wildcard)
                load_functions.append(geom_func)

                if res_wildcard is None:
                    has_results_list.append(False)
                else:
                    has_results_list.append(True)

            # the list of formats that will be selectable in some odd syntax
            # that pyqt uses
            wildcard = ';;'.join(wildcard_list)

            # get the filter index and filename
            if infile_name is not None and geometry_format is not None:
                filter_index = formats.index(geometry_format)
            else:
                title = 'set_legend_title File to Load'
                wildcard_index, infile_name = self.create_load_file_dialog(wildcard, title)
                if not infile_name:
                    # user clicked cancel
                    is_failed = True
                    return is_failed, None
                filter_index = wildcard_list.index(wildcard_index)

            geometry_format = formats[filter_index]
            load_function = load_functions[filter_index]
            unused_has_results = has_results_list[filter_index]
        return is_failed, (infile_name, load_function, filter_index, formats, geometry_format)

    def on_load_results(self, out_filename=None):
        """
        Loads a results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        out_filename : str / None
            the path to the results file

        """
        geometry_format = self.gui.format
        if self.gui.format is None:
            msg = 'on_load_results failed:  You need to load a file first...'
            self.gui.log_error(msg)
            return
            #raise RuntimeError(msg)

        if out_filename in [None, False]:
            title = 'Select a Results File for %s' % self.gui.format
            wildcard = None
            load_function = None

            for fmt in self.gui.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                if geometry_format == fmt_name:
                    wildcard = _reswild
                    load_function = _resfunc
                    break
            else:
                msg = f'format={geometry_format!r} is not supported'
                self.gui.log_error(msg)
                raise RuntimeError(msg)

            if wildcard is None:
                msg = f'format={geometry_format!r} has no method to load results'
                self.gui.log_error(msg)
                return
            out_filename = self.create_load_file_dialog(wildcard, title)[1]
        else:

            for fmt in self.gui.fmts:
                fmt_name, _major_name, _geowild, _geofunc, _reswild, _resfunc = fmt
                #print('fmt_name=%r geometry_format=%r' % (fmt_name, geometry_format))
                if fmt_name == geometry_format:
                    load_function = _resfunc
                    break
            else:
                msg = ('format=%r is not supported.  '
                       'Did you load a geometry model?' % geometry_format)
                self.gui.log_error(msg)
                raise RuntimeError(msg)

        if out_filename == '':
            return
        if isinstance(out_filename, str):
            out_filename = [out_filename]
        for out_filenamei in out_filename:
            if not os.path.exists(out_filenamei):
                msg = 'result file=%r does not exist' % out_filenamei
                self.gui.log_error(msg)
                return
                #raise IOError(msg)
            self.gui.last_dir = os.path.split(out_filenamei)[0]

            name = 'main'
            if name in self.model_objs:
                model = self.model_objs[name]
                if fmt[0] == 'nastran3':
                    fmti, geo1, geo2, res1, res2 = model.get_nastran3_wildcard_geometry_results_functions()
                    load_function = res2

            try:
                load_function(out_filenamei)
            except Exception: #  as e
                msg = traceback.format_exc()
                self.gui.log_error(msg)
                print(msg)
                return
                #raise

            self.gui.out_filename = out_filenamei
            msg = '%s - %s - %s' % (self.gui.format, self.gui.infile_name, out_filenamei)
            self.gui.window_title = msg
            print("self.load_actions.on_load_results(%r)" % out_filenamei)
            self.gui.out_filename = out_filenamei
            self.gui.log_command("on_load_results(%r)" % out_filenamei)

    #---------------------------------------------------------------------------
    def on_load_custom_results(self, out_filename=None, restype=None,
                               stop_on_failure: bool=False) -> bool:
        """will be a more generalized results reader"""
        is_failed, out_filename, iwildcard = self._on_load_custom_results_load_filename(
            out_filename=out_filename, restype=restype)

        if is_failed:
            if stop_on_failure:  # pragma: no cover
                raise RuntimeError(f'failed getting filename={out_filename!r}')
            return is_failed
        if out_filename == '':
            is_failed = True
            return is_failed

        is_failed = True
        if not os.path.exists(out_filename):
            msg = 'result file=%r does not exist' % out_filename
            self.gui.log_error(msg)
            if stop_on_failure:  # pragma: no cover
                raise RuntimeError(msg)
            return is_failed

        self._set_last_dir(out_filename)
        try:
            if iwildcard == 0:
                self._on_load_nodal_elemental_results(
                    'Nodal', out_filename, stop_on_failure=stop_on_failure)
                restype = 'Node'
            elif iwildcard == 1:
                self._on_load_nodal_elemental_results(
                    'Elemental', out_filename, stop_on_failure=stop_on_failure)
                restype = 'Element'
            elif iwildcard == 2:
                self._load_deflection(out_filename)
                restype = 'Deflection'
            elif iwildcard == 3:
                self._load_force(out_filename)
                restype = 'Force'
            elif iwildcard == 4:
                self.load_patran_nod(out_filename)
                restype = 'Patran_nod'
            else:  # pramga: no cover
                raise NotImplementedError('iwildcard = %s' % iwildcard)
        except Exception:
            msg = traceback.format_exc()
            self.gui.log_error(msg)
            if stop_on_failure:  # pragma: no cover
                raise RuntimeError(msg)
            return is_failed
        self.gui.log_command('self.load_actions.on_load_custom_results('
                             f'{out_filename!r}, restype={restype!r})')
        is_failed = False
        return is_failed

    def _on_load_custom_results_load_filename(self, out_filename=None, restype=None):
        is_failed = True
        #unused_geometry_format = self.format
        if self.gui.format is None:
            msg = 'on_load_results failed:  You need to load a file first...'
            self.gui.log_error(msg)
            return is_failed, None, None

        if out_filename in [None, False]:
            title = 'Select a Custom Results File for %s' % (self.gui.format)

            #print('wildcard_level =', wildcard_level)
            #self.wildcard_delimited = 'Delimited Text (*.txt; *.dat; *.csv)'
            fmts = [
                'Node - Delimited Text (*.txt; *.dat; *.csv)',
                'Element - Delimited Text (*.txt; *.dat; *.csv)',
                'Nodal Deflection - Delimited Text (*.txt; *.dat; *.csv)',
                'Nodal Force - Delimited Text (*.txt; *.dat; *.csv)',
                'Patran nod (*.nod)',
            ]
            fmt = ';;'.join(fmts)
            wildcard_level, out_filename = self.create_load_file_dialog(fmt, title)
            if not out_filename:
                return is_failed, None, None # user clicked cancel
            iwildcard = fmts.index(wildcard_level)
        else:
            fmts = [
                'node', 'element', 'deflection', 'force', 'patran_nod',
            ]
            iwildcard = fmts.index(restype.lower())
        is_failed = False
        return is_failed, out_filename, iwildcard

    def _load_deflection(self, out_filename):
        """loads a deflection file"""
        self._load_deflection_force(out_filename, is_deflection=True, is_force=False)

    def _load_force(self, out_filename):
        """loads a force file"""
        self._load_deflection_force(out_filename, is_deflection=False, is_force=True)

    def _load_deflection_force(self, out_filename, is_deflection=False, is_force=False):
        out_filename_short = os.path.basename(out_filename)
        A, nids_index, fmt_dict, headers = load_deflection_csv(out_filename)
        #nrows, ncols, fmts
        header0 = headers[0]
        result0 = A[header0]
        nrows = result0.shape[0]

        nnodes = self.gui.nnodes
        if nrows != nnodes:
            #'nrows=%s nnodes=%s' % (nrows, self.gui.nnodes)
            self.log.warning(f'The deflection CSV has {nrows:d} rows, but '
                             f'there are {nnodes:d} nodes in the model.  '
                             'Verify that the result is for the correct '
                             "model and that it's not an elemental result.")
            A = _resize_array(A, nids_index, self.gui.node_ids, nrows, nnodes)

        result_type = 'node'
        self._add_cases_to_form(A, fmt_dict, headers, result_type,
                                out_filename_short, update=True, is_scalar=False,
                                is_deflection=is_deflection, is_force=is_force)

    def _on_load_nodal_elemental_results(self, result_type: str,
                                         out_filename=None,
                                         stop_on_failure: bool=False):
        """
        Loads a CSV/TXT results file.  Must have called on_load_geometry first.

        Parameters
        ----------
        result_type : str
            'Nodal', 'Elemental'
        out_filename : str / None
            the path to the results file

        """
        try:
            self._load_csv(result_type, out_filename, stop_on_failure=stop_on_failure)
        except Exception:
            msg = traceback.format_exc()
            self.gui.log_error(msg)
            if stop_on_failure:  # pragma: no cover
                raise
            #return
            raise

        #if 0:
            #self.out_filename = out_filename
            #msg = '%s - %s - %s' % (self.format, self.infile_name, out_filename)
            #self.window_title = msg
            #self.out_filename = out_filename

    #---------------------------------------------------------------------------
    def on_load_user_geom(self, csv_filename: Optional[str]=None,
                          name: Optional[str]=None,
                          color: Optional[list[float]]=None) -> None:
        """
        Loads a User Geometry CSV File of the form:

        #    id  x    y    z
        GRID, 1, 0.2, 0.3, 0.3
        GRID, 2, 1.2, 0.3, 0.3
        GRID, 3, 2.2, 0.3, 0.3
        GRID, 4, 5.2, 0.3, 0.3
        grid, 5, 5.2, 1.3, 2.3  # case insensitive

        #    ID, nodes
        BAR,  1, 1, 2
        TRI,  2, 1, 2, 3
        # this is a comment

        QUAD, 3, 1, 5, 3, 4
        QUAD, 4, 1, 2, 3, 4  # this is after a blank line

        #RESULT,4,CENTROID,AREA(%f),PROPERTY_ID(%i)
        # in element id sorted order: value1, value2
        #1.0, 2.0 # bar
        #1.0, 2.0 # tri
        #1.0, 2.0 # quad
        #1.0, 2.0 # quad

        #RESULT,NODE,NODEX(%f),NODEY(%f),NODEZ(%f)
        # same difference

        #RESULT,VECTOR3,GEOM,DXYZ
        # 3xN

        Parameters
        ----------
        csv_filename : str (default=None -> load a dialog)
            the path to the user geometry CSV file
        name : str (default=None -> extract from fname)
            the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0

        """
        gui = self.gui
        if csv_filename in [None, False]:
            title = 'Load User Geometry'
            csv_filename = gui.load_actions.create_load_file_dialog(
                #gui.wildcard_delimited + ';;STL (*.stl)', title)[1]
                gui.wildcard_delimited, title)[1]
            if not csv_filename:
                return
            self._set_last_dir(csv_filename)
        assert isinstance(csv_filename, str), csv_filename

        if color is None:
            # we mod the num_user_points so we don't go outside the range
            icolor = gui.num_user_points % len(gui.color_order)
            color = gui.color_order[icolor]
        if name is None:
            name = os.path.basename(csv_filename).rsplit('.', 1)[0]

        self._add_user_geometry(csv_filename, name, color)
        gui.log_command(f'self.on_load_user_geom({csv_filename}, name={name!r}, color={str(color)})')

    def _add_user_geometry(self, csv_filename: str,
                           name: str,
                           color: list[float]) -> None:
        """
        helper method for ``on_load_user_geom``

        A custom geometry can be the pyNastran custom form or an STL

        """
        gui: MainWindow = self.gui
        tool_actions: ToolActions = gui.tool_actions
        if name in gui.geometry_actors:
            msg = f'Name: {name!r} is already in geometry_actors\nChoose a different name.'
            raise ValueError(msg)
        if len(name) == 0:
            msg = f'Invalid Name: name={name!r}'
            raise ValueError(msg)

        point_name = name + '_point'
        geom_name = name + '_geom'

        grid_ids, xyz, bars, tris, quads = load_user_geom(
            csv_filename, self.gui.log, encoding='latin1')
        nbars = len(bars)
        ntris = len(tris)
        nquads = len(quads)
        nelements = nbars + ntris + nquads
        gui.create_alternate_vtk_grid(point_name, color=color, opacity=1.0,
                                      point_size=5, representation='point')

        if nelements > 0:
            nid_map = {}
            i = 0
            for nid in grid_ids:
                nid_map[nid] = i
                i += 1
            gui.create_alternate_vtk_grid(geom_name, color=color, opacity=1.0,
                                          line_width=5, representation='toggle')

        # allocate
        nnodes = len(grid_ids)
        #self.alt_grids[point_name].Allocate(npoints, 1000)
        #if nelements > 0:
            #self.alt_grids[geom_name].Allocate(npoints, 1000)

        alt_grid = gui.alt_grids[point_name]
        geom_grid = gui.alt_grids[geom_name]

        add_user_geometry(
            alt_grid, geom_grid,
            xyz, nid_map, nnodes,
            bars, tris, quads,
            nelements, nbars, ntris, nquads)

        # create actor/mapper
        tool_actions.add_alt_geometry(alt_grid, point_name)
        if nelements > 0:
            tool_actions.add_alt_geometry(geom_grid, geom_name)

        # set representation to points
        #self.geometry_properties[point_name].representation = 'point'
        #self.geometry_properties[geom_name].representation = 'toggle'
        #actor = self.geometry_actors[name]
        #prop = actor.GetProperty()
        #prop.SetRepresentationToPoints()
        #prop.SetPointSize(4)

    #---------------------------------------------------------------------------
    def load_patran_nod(self, nod_filename: str) -> None:
        """reads a Patran formatted *.nod file"""
        A, fmt_dict, headers = load_patran_nod(nod_filename, self.gui.node_ids)

        out_filename_short = os.path.relpath(nod_filename)
        result_type = 'node'
        self._add_cases_to_form(A, fmt_dict, headers, result_type,
                                out_filename_short, update=True,
                                is_scalar=True)

    def _load_csv(self, result_type: str,
                  out_filename: str,
                  stop_on_failure: bool=False) -> None:
        """
        common method between:
          - on_add_nodal_results(filename)
          - on_add_elemental_results(filename)

        Parameters
        ----------
        result_type : str
            ???
        out_filename : str
            the CSV filename to load

        """
        out_filename_short = os.path.relpath(out_filename)
        A, fmt_dict, headers = load_csv(out_filename)
        #nrows, ncols, fmts
        header0 = headers[0]
        result0 = A[header0]
        nrows = result0.size

        if result_type == 'Nodal':
            nnodes = self.gui.nnodes
            if nrows != nnodes:
                self.log.warning('The fringe CSV has %i rows, but there are %i nodes in the '
                                 'model.  Verify that the result is for the correct model and '
                                 "that it's not an elemental result." % (nrows, nnodes))
                A = _resize_array(A, A['index'], self.gui.node_ids, nrows, nnodes)
            result_type2 = 'node'
        elif result_type == 'Elemental':
            nelements = self.gui.nelements
            if nrows != nelements:
                self.log.warning('The fringe CSV has %i rows, but there are %i elements in the '
                                 'model.  Verify that the result is for the correct model and '
                                 "that it's not a nodal result." % (nrows, nelements))
                A = _resize_array(A, A['index'], self.gui.element_ids, nrows, nelements)
            result_type2 = 'centroid'
        else:  # pragma: no cover
            raise NotImplementedError(f'result_type={result_type!r}')

        #num_ids = len(ids)
        #if num_ids != nrows:
            #A2 = {}
            #for key, matrix in A.items():
                #fmt = fmt_dict[key]
                #assert fmt not in ['%i'], 'fmt=%r' % fmt
                #if len(matrix.shape) == 1:
                    #matrix2 = np.full(num_ids, dtype=matrix.dtype)
                    #iids = np.searchsorted(ids, )
            #A = A2
        self._add_cases_to_form(A, fmt_dict, headers, result_type2,
                                out_filename_short, update=True, is_scalar=True)

    def _add_cases_to_form(self, A: dict[str, np.ndarray],
                           fmt_dict,
                           headers: list[str],
                           result_type: str,
                           out_filename_short: str,
                           update: bool=True,
                           is_scalar: bool=True,
                           is_deflection: bool=False,
                           is_force: bool=False):
        """
        common method between:
          - _load_csv
          - _load_deflection_csv

        Parameters
        ----------
        A : dict[key] = (n, m) array
            the numpy arrays
            key : str
                the name
            n : int
                number of nodes/elements
            m : int
                secondary dimension
                N/A : 1D array
                3 : deflection
        fmt_dict : dict[header] = fmt
            the format of the arrays
            header : str
                the name
            fmt : str
                '%i', '%f'
        headers : list[str]???
            the titles???
        result_type : str
            'node', 'centroid'
        out_filename_short : str
            the display name
        update : bool; default=True
            update the res_widget

        # A = np.loadtxt('loadtxt_spike.txt', dtype=('float,int'))
        # dtype=[('f0', '<f8'), ('f1', '<i4')])
        # A['f0']
        # A['f1']
        """
        #print('A =', A)
        formi: list[tuple[str, Any, Any]]= []
        form = self.gui.get_form()
        icase = len(self.gui.case_keys)
        islot = 0
        for case_key in self.gui.case_keys:
            if isinstance(case_key, tuple):
                islot = case_key[0]
                break

        #assert len(headers) > 0, 'headers=%s' % (headers)
        #assert len(headers) < 50, 'headers=%s' % (headers)
        settings: Settings = self.gui.settings
        for header in headers:
            if is_scalar:
                out = create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                                     colormap=settings.colormap)
            else:
                out = create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                                     is_deflection=is_deflection, is_force=is_force,
                                     dim_max=self.gui.settings.dim_max, xyz_cid0=self.gui.xyz_cid0,
                                     colormap=settings.colormap)
            res_obj, title = out

            #cases[icase] = (stress_res, (subcase_id, 'Stress - isElementOn'))
            #form_dict[(key, itime)].append(('Stress - IsElementOn', icase, []))
            #key = (res_obj, (0, title))
            self.gui.case_keys.append(icase)
            self.gui.result_cases[icase] = (res_obj, (islot, title))
            formi.append((header, icase, []))

            # TODO: double check this should be a string instead of an int
            self.gui.label_actors[icase] = []
            self.gui.label_ids[icase] = set()
            icase += 1
        form.append((out_filename_short, None, formi))

        self.gui.ncases += len(headers)
        #cases[(ID, 2, 'Region', 1, 'centroid', '%i')] = regions
        if update:
            self.gui.res_widget.update_results(form, 'main')

    #---------------------------------------------------------------------------
    def create_load_file_dialog(self, qt_wildcard: str, title: str,
                                default_filename: Optional[str]=None) -> tuple[str, str]:
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        #fname, flt = QFileDialog.getOpenFileName(
            #self, title, default_filename, file_types, options=options)
        #flt = str(filt).strip()
        #return fname, flt

        use_startup_directory = self.settings.use_startup_directory
        startup_directory = self.settings.startup_directory.strip()
        last_dir = self.gui.last_dir

        basedir = ''
        if use_startup_directory and os.path.exists(startup_directory):
            basedir = startup_directory
        else:
            basedir = last_dir

        if not os.path.exists(basedir):
            basedir = ''

        #if startup_directory == '' and default_filename is None and os.path.exists(last_dir):
        #elif os.path.exists(startup_directory):
            #basedir = startup_directory

        fname, wildcard_level = getopenfilename(
            parent=self.gui, caption=title,
            basedir=basedir, filters=qt_wildcard,
            selectedfilter='', options=None)
        return wildcard_level, fname

    #def create_load_file_dialog2(self, qt_wildcard, title):
        ## getOpenFileName return QString and we want Python string
        ##title = 'Load a Tecplot Geometry/Results File'
        #last_dir = ''
        ##qt_wildcard = ['Tecplot Hex Binary (*.tec; *.dat)']
        #dialog = MultiFileDialog()
        #dialog.setWindowTitle(title)
        #dialog.setDirectory(self.last_dir)
        #dialog.setFilters(qt_wildcard.split(';;'))
        #if dialog.exec_() == QtGui.QDialog.Accepted:
            #outfiles = dialog.selectedFiles()
            #wildcard_level = dialog.selectedFilter()
            #return str(wildcard_level), str(fname)
        #return None, None
    # --------------------------------------------------------------------------

    def on_run_script(self, python_file=False) -> bool:
        """pulldown for running a python script"""
        is_passed = False
        gui: MainWindow = self.gui
        if python_file in [None, False]:
            title = 'Choose a Python Script to Run'
            wildcard = 'Python (*.py)'
            infile_name = self.create_load_file_dialog(
                wildcard, title, gui._default_python_file)[1]
            if not infile_name:
                return is_passed # user clicked cancel

            #python_file = os.path.join(script_path, infile_name)
            python_file = os.path.join(infile_name)
            self._set_last_dir(python_file)

        if not os.path.exists(python_file):
            msg = 'python_file = %r does not exist' % python_file
            gui.log_error(msg)
            return is_passed

        with open(python_file, 'r') as python_file_obj:
            txt = python_file_obj.read()
        is_passed = gui._execute_python_code(txt, show_msg=False)
        if not is_passed:
            return is_passed
        gui._default_python_file = python_file
        gui.log_command(f'self.on_run_script({python_file!r})')
        print(f'self.on_run_script({python_file!r})')
        return is_passed


def _resize_array(A: dict[str, np.ndarray],
                  nids_index: np.ndarray,
                  node_ids: np.ndarray,
                  nrows: int, nnodes: int) -> dict[str, np.ndarray]:
    """
    Resizes an array to be the right size.

    Let's say we have 5 nodes in the output csv that aren't in the model.
    We need to filter them out.  Alternatively, we may have extra nodes.

    We need to get the array to have results at all the node ids.

    Parameters
    ----------
    A : np.ndarray
        the dictionary-like results array
    node_ids : int np.ndarray
        the node ids in the model
    nrows : int
        the number of rows in the csv; same length as nids_index
    nnodes : int
        the number of nodes in the real model; same length as node_ids

    Returns
    -------
    A2 : np.ndarray
        the properly sized dictionary-like results array

    """
    # we might have extra nodes or missing nodes, so
    # find the list of valid indices
    inids = np.searchsorted(node_ids, nids_index)
    iexist = np.where(inids < nnodes)[0]

    A2 = {}
    for key, Ai in A.items():
        #print('Ai.shape', Ai.shape, len(iexist))
        if key == 'index':
            A2[key] = node_ids
            continue

        if len(Ai.shape) == 1:
            if isinstance(Ai[0], (np.float32, np.float64)):
                new_array = np.full((nnodes, ), np.nan, dtype=Ai.dtype)
            elif isinstance(Ai[0], (np.int32, np.int64)):
                new_array = np.full((nnodes, ), -1, dtype=Ai.dtype)
            else:
                raise NotImplementedError(Ai[0].dtype)
            #print('iexist', iexist.shape, Ai.shape)
            new_array[iexist] = Ai[iexist]
            A2[key] = new_array

        elif len(Ai.shape) == 2:
            ncols = Ai.shape[1]
            if isinstance(Ai[0, 0], (np.float32, np.float64)):
                new_array = np.full((nnodes, ncols), np.nan, dtype=Ai.dtype)
            elif isinstance(Ai[0, 0], (np.int32, np.int64)):
                new_array = np.full((nnodes, ncols), -1, dtype=Ai.dtype)
            else:
                raise NotImplementedError(Ai[0].dtype)
            #print('iexist', iexist.shape, Ai.shape)
            new_array[iexist] = Ai[iexist]
            A2[key] = new_array
        else:
            raise  NotImplementedError(Ai.shape)

        #A2[key] = Ai[iexist]
        #print('A2[%s].shape = %s' % (key, A2[key].shape))
    #print()
    return A2
