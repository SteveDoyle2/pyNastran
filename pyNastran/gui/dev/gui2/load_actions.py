from __future__ import annotations
import os
from typing import Optional, TYPE_CHECKING
import time as time_module
import traceback

from qtpy.compat import getopenfilename
from pyNastran.utils import print_bad_path
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui2 import MainWindow2
    from pyNastran.gui.gui_objects.settings import Settings
IS_TESTING = False

class LoadActions:
    """performance mode should be handled in the main gui to minimize flipping"""
    def __init__(self, gui: MainWindow2):
        self.gui = gui

    @property
    def log(self):
        """links the the GUI's log"""
        return self.gui.log

    def on_load_geometry(self,
                         infile_name: Optional[str]=None,
                         geometry_format: Optional[str]=None,
                         name: str='main',
                         plot: bool=True,
                         raise_error: bool=False):
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
        log = self.log
        infile_name, load_function, filter_index, formats, geometry_format2 = out
        if load_function is not None:
            last_dir = os.path.split(infile_name)[0]
            settings: Settings = self.gui.settings
            self.gui.last_dir = last_dir
            settings.startup_directory = last_dir

            if self.gui.name == '':
                name = 'main'
            else:
                print('name = %r' % name)

            active_name = self.gui.name
            alt_grids = self.gui.alt_grids
            grid_mappers = self.gui.grid_mappers
            if name != active_name and active_name in grid_mappers:
                #scalar_range = self.grid_selected.GetScalarRange()
                #self.grid_mapper.SetScalarRange(scalar_range)
                grid_mappers[active_name].ScalarVisibilityOff()
                #self.grid_mapper.SetLookupTable(self.color_function)
            self.gui.name = name
            self.gui._reset_model(name)

            # reset alt grids
            names = self.gui.alt_grids.keys()
            for name in names:
                alt_grid = alt_grids[name]
                alt_grid.Reset()
                alt_grid.Modified()

            if not os.path.exists(infile_name) and geometry_format:
                msg = 'input file=%r does not exist' % infile_name
                log.error(msg)
                log.error(print_bad_path(infile_name))
                return

            # clear out old data
            if self.gui.model_type is not None:
                clear_name = 'clear_' + self.gui.model_type
                try:
                    dy_method = getattr(self, clear_name)  # 'self.clear_nastran()'
                    dy_method()
                except Exception:
                    self.gui.log.error("method %r does not exist" % clear_name)
            log.info("reading %s file %r" % (geometry_format, infile_name))

            try:
                time0 = time_module.time()

                if geometry_format2 in self.gui.format_class_map:
                    # initialize the class
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
                log.error(msg)
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
            assert len(load_functions) > 0, load_functions

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

    def create_load_file_dialog(self, qt_wildcard: str, title: str,
                                default_filename: Optional[str]=None) -> tuple[str, str]:
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
