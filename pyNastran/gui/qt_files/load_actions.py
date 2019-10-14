import os
import sys
import traceback
import time as time_module

import numpy as np
from qtpy.compat import getopenfilename
#from qtpy.QtWidgets import QFileDialog
from pyNastran.bdf.patran_utils.read_patran_custom_results import load_patran_nod
from pyNastran.utils import print_bad_path

from pyNastran.gui.utils.load_results import load_csv, load_deflection_csv
from pyNastran.gui.utils.load_results import create_res_obj
IS_TESTING = 'test' in sys.argv[0]


class LoadActions:
    """performance mode should be handled in the main gui to minimize flipping"""
    def __init__(self, gui):
        self.gui = gui

    @property
    def log(self):
        """links the the GUI's log"""
        return self.gui.log

    def on_load_geometry(self, infile_name=None, geometry_format=None, name='main',
                         plot=True, raise_error=False):
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
        raise_error : bool; default=True
            stop the code if True

        """
        assert isinstance(name, str), 'name=%r type=%s' % (name, type(name))
        is_failed, out = self._load_geometry_filename(
            geometry_format, infile_name)
        print("is_failed =", is_failed)
        if is_failed:
            return

        has_results = False
        infile_name, load_function, filter_index, formats, geometry_format2 = out
        if load_function is not None:
            self.gui.last_dir = os.path.split(infile_name)[0]

            if self.gui.name == '':
                name = 'main'
            else:
                print('name = %r' % name)

            if name != self.gui.name:
                #scalar_range = self.grid_selected.GetScalarRange()
                #self.grid_mapper.SetScalarRange(scalar_range)
                self.gui.grid_mapper.ScalarVisibilityOff()
                #self.grid_mapper.SetLookupTable(self.color_function)
            self.gui.name = name
            self.gui._reset_model(name)

            # reset alt grids
            names = self.gui.alt_grids.keys()
            for name in names:
                self.gui.alt_grids[name].Reset()
                self.gui.alt_grids[name].Modified()

            if not os.path.exists(infile_name) and geometry_format:
                msg = 'input file=%r does not exist' % infile_name
                self.gui.log_error(msg)
                self.gui.log_error(print_bad_path(infile_name))
                return

            # clear out old data
            if self.gui.model_type is not None:
                clear_name = 'clear_' + self.gui.model_type
                try:
                    dy_method = getattr(self, clear_name)  # 'self.clear_nastran()'
                    dy_method()
                except:
                    print("method %r does not exist" % clear_name)
            self.gui.log_info("reading %s file %r" % (geometry_format, infile_name))

            try:
                time0 = time_module.time()

                if geometry_format2 in self.gui.format_class_map:
                    # intialize the class
                    #print('geometry_format=%r geometry_format2=%s' % (geometry_format, geometry_format2))

                    # TODO: was geometry_format going into this...
                    cls = self.gui.format_class_map[geometry_format2](self.gui)

                    function_name2 = 'load_%s_geometry' % geometry_format2
                    load_function2 = getattr(cls, function_name2)
                    has_results = load_function2(infile_name, name=name, plot=plot)
                else:
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
                self.gui.log_error(msg)
                if raise_error or self.gui.dev:
                    raise
                #return
            #self.vtk_panel.Update()
            self.gui.rend.ResetCamera()

        # the model has been loaded, so we enable load_results
        if filter_index >= 0:
            self.gui.format = formats[filter_index].lower()
            unused_enable = has_results
            #self.load_results.Enable(enable)
        else: # no file specified
            return
        #print("on_load_geometry(infile_name=%r, geometry_format=None)" % infile_name)
        self.gui.infile_name = infile_name
        self.gui.out_filename = None
        #if self.out_filename is not None:
            #msg = '%s - %s - %s' % (self.format, self.infile_name, self.out_filename)

        if name == 'main':
            msg = '%s - %s' % (self.gui.format, self.gui.infile_name)
            self.gui.window_title = msg
            self.gui.update_menu_bar()
            main_str = ''
        else:
            main_str = ', name=%r' % name

        self.gui.log_command("on_load_geometry(infile_name=%r, geometry_format=%r%s)" % (
            infile_name, self.gui.format, main_str))

    def _load_geometry_filename(self, geometry_format, infile_name):
        """gets the filename and format"""
        wildcard = ''
        is_failed = False

        if geometry_format and geometry_format.lower() not in self.gui.supported_formats:
            is_failed = True
            msg = 'The import for the %r module failed.\n' % geometry_format
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
                title = 'Choose a Geometry File to Load'
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
                msg = 'format=%r is not supported' % geometry_format
                self.gui.log_error(msg)
                raise RuntimeError(msg)

            if wildcard is None:
                msg = 'format=%r has no method to load results' % geometry_format
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

            try:
                load_function(out_filenamei)
            except: #  as e
                msg = traceback.format_exc()
                self.gui.log_error(msg)
                print(msg)
                #return
                raise

            self.gui.out_filename = out_filenamei
            msg = '%s - %s - %s' % (self.gui.format, self.gui.infile_name, out_filenamei)
            self.gui.window_title = msg
            print("on_load_results(%r)" % out_filenamei)
            self.gui.out_filename = out_filenamei
            self.gui.log_command("on_load_results(%r)" % out_filenamei)

    def on_load_custom_results(self, out_filename=None, restype=None, stop_on_failure=False):
        """will be a more generalized results reader"""
        is_failed, out_filename, iwildcard = self._on_load_custom_results_load_filename(
            out_filename=out_filename, restype=restype)

        if is_failed:
            if stop_on_failure:  # pragma: no cover
                raise RuntimeError('failed getting filename')
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

        try:
            if iwildcard == 0:
                self._on_load_nodal_elemental_results('Nodal', out_filename, stop_on_failure=stop_on_failure)
                restype = 'Node'
            elif iwildcard == 1:
                self._on_load_nodal_elemental_results('Elemental', out_filename, stop_on_failure=stop_on_failure)
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
            else:
                raise NotImplementedError('iwildcard = %s' % iwildcard)
        except:
            msg = traceback.format_exc()
            self.gui.log_error(msg)
            if stop_on_failure:  # pragma: no cover
                raise RuntimeError(msg)
            return is_failed
        self.gui.log_command("on_load_custom_results(%r, restype=%r)" % (out_filename, restype))
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
            self.log.warning('The deflection CSV has %i rows, but there are %i nodes in the model.'
                             "  Verify that the result is for the correct model and that it's "
                             'not an elemental result.' % (nrows, nnodes))
            A = _resize_array(A, nids_index, self.gui.node_ids, nrows, nnodes)

        result_type = 'node'
        self._add_cases_to_form(A, fmt_dict, headers, result_type,
                                out_filename_short, update=True, is_scalar=False,
                                is_deflection=is_deflection, is_force=is_force)

    def _on_load_nodal_elemental_results(self, result_type, out_filename=None, stop_on_failure=False):
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
        except:
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

    def load_patran_nod(self, nod_filename):
        """reads a Patran formatted *.nod file"""
        A, fmt_dict, headers = load_patran_nod(nod_filename, self.gui.node_ids)

        out_filename_short = os.path.relpath(nod_filename)
        result_type = 'node'
        self._add_cases_to_form(A, fmt_dict, headers, result_type,
                                out_filename_short, update=True,
                                is_scalar=True)

    def _load_csv(self, result_type, out_filename, stop_on_failure=False):
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
        else:
            raise NotImplementedError('result_type=%r' % result_type)

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

    def _add_cases_to_form(self, A, fmt_dict, headers, result_type,
                           out_filename_short, update=True, is_scalar=True,
                           is_deflection=False, is_force=False):
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
        headers : List[str]???
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
        formi = []
        form = self.gui.get_form()
        icase = len(self.gui.case_keys)
        islot = 0
        for case_key in self.gui.case_keys:
            if isinstance(case_key, tuple):
                islot = case_key[0]
                break

        #assert len(headers) > 0, 'headers=%s' % (headers)
        #assert len(headers) < 50, 'headers=%s' % (headers)
        for header in headers:
            if is_scalar:
                out = create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                                     colormap='jet')
            else:
                out = create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                                     is_deflection=is_deflection, is_force=is_force,
                                     dim_max=self.gui.settings.dim_max, xyz_cid0=self.gui.xyz_cid0,
                                     colormap='jet')
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

    def create_load_file_dialog(self, qt_wildcard, title, default_filename=None):
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        #fname, flt = QFileDialog.getOpenFileName(
            #self, title, default_filename, file_types, options=options)
        #flt = str(filt).strip()
        #return fname, flt

        if default_filename is None:
            default_filename = self.gui.last_dir
        fname, wildcard_level = getopenfilename(
            parent=self.gui, caption=title,
            basedir=default_filename, filters=qt_wildcard,
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


def _resize_array(A, nids_index, node_ids, nrows, nnodes):
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
