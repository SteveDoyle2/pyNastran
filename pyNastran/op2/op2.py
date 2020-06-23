#pylint: disable=W0201,W0223,R0901,R0902,R0904
"""
Defines the main OP2 class.  Defines:

 - read_op2(op2_filename=None, combine=True, subcases=None,
            exclude_results=None, include_results=None,
            log=None, debug=True, debug_file=None, build_dataframe=None,
            skip_undefined_matrices=True, mode='msc', encoding=None)

 - OP2(debug=True, log=None, debug_file=None, mode='msc')
   - build_dataframe()
   - combine_results(combine=True)
   - create_objects_from_matrices()
   - object_attributes(mode='public', keys_to_skip=None, filter_properties=False)
   - object_methods(mode='public', keys_to_skip=None)
   - print_subcase_key()
   - read_op2(op2_filename=None, combine=True, build_dataframe=None,
              skip_undefined_matrices=False, encoding=None)
   - set_mode(mode)
   - transform_displacements_to_global(i_transform, coords, xyz_cid0=None, debug=False)
   - transform_gpforce_to_global(nids_all, nids_transform, i_transform, coords, xyz_cid0=None)

"""
from __future__ import annotations
import sys
from collections import defaultdict
from pickle import load, dump, dumps
from typing import List, Dict, Tuple, Optional, Any, TYPE_CHECKING

import numpy as np

#import pyNastran
from pyNastran.utils import (
    object_attributes, object_methods, ipython_info)
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.monpnt import MONPNT1, MONPNT3

from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import (SortCodeError, DeviceCodeError,
                                  FortranMarkerError,  SixtyFourBitError,
                                  OverwriteTableError)
from pyNastran.op2.writer.op2_writer import OP2Writer
#from pyNastran.op2.op2_interface.op2_f06_common import Op2F06Attributes
from pyNastran.op2.op2_interface.op2_scalar import OP2_Scalar
from pyNastran.op2.op2_interface.transforms import (
    transform_displacement_to_global, transform_gpforce_to_globali)
from pyNastran.utils import check_path
if TYPE_CHECKING:  # pragma: no cover
    from h5py import File as H5File


class OP2(OP2_Scalar, OP2Writer):
    _properties = ['is_real', 'is_complex', 'is_random',
                   '_sort_method', 'is_sort1', 'is_sort2',
                   'matrix_tables', 'table_name_str']

    def __init__(self,
                 debug: bool=True,
                 log: Any=None,
                 debug_file: Optional[str]=None,
                 mode: Optional[str]=None) -> None:
        """
        Initializes the OP2 object

        Parameters
        ----------
        debug : bool; default=False
            enables the debug log and sets the debug in the logger
        log : Log()
            a logging object to write debug messages to
         (.. seealso:: import logging)
        debug_file : str; default=None (No debug)
            sets the filename that will be written to
        mode : str; default=None -> 'msc'
            {msc, nx}

        """
        self.encoding = None
        self.mode = mode
        if mode is not None:
            self.set_mode(mode)
        make_geom = False
        assert make_geom is False, make_geom
        OP2_Scalar.__init__(self, debug=debug, log=log, debug_file=debug_file)
        self.ask = False
        self.post = None
        self.table_count = defaultdict(int)

    def __del__(self) -> None:
        if hasattr(self, 'h5_file') and self.h5_file is not None:
            self.h5_file.close()

    def object_attributes(self, mode: str='public', keys_to_skip: Optional[List[str]]=None,
                          filter_properties: bool=False) -> List[str]:
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : List[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong

        """
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public',
                       keys_to_skip: Optional[List[str]]=None) -> List[str]:
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : List[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong

        """
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def __eq__(self, op2_model) -> bool:
        """
        Diffs the current op2 model vs. another op2 model.
        Crashes if they're not equal.
        """
        try:
            is_equal = self.assert_op2_equal(op2_model,
                                             stop_on_failure=True, debug=False)
        except (AssertionError, ValueError):
            is_equal = False
            raise
        return is_equal

    def assert_op2_equal(self, op2_model, skip_results: Optional[List[str]]=None,
                         stop_on_failure: bool=True, debug: bool=False) -> None:
        """
        Diffs the current op2 model vs. another op2 model.

        Parameters
        ----------
        op2_model : OP2()
            the model to compare to
        skip_results : List[str]; default=None -> []
            results that shouldn't be compred
        stop_on_failure : bool; default=True
            True : Crashes if they're not equal
            False : Go to the next object
        debug : bool; default=False
            give slightly more debugging messages

        Returns
        -------
        is_equal : bool
            are the objects equal?

        Raises
        ------
        AssertionError/ValueError : stop_on_failure=True and and error occurred
        NotImplementedError : this is a sign of an unsupported object

        """
        if skip_results is None:
            skip_results = set()
        else:
            skip_results = set(skip_results)

        skip_results.add('gpdt')
        skip_results.add('bgpdt')
        skip_results.add('eqexin')
        skip_results.add('psds')

        if not self.read_mode == op2_model.read_mode:
            self.log.warning('self.read_mode=%s op2_model.read_mode=%s ... assume True' % (
                self.read_mode, op2_model.read_mode))
            return True

        table_types = self.get_table_types()
        for table_type in table_types:
            if table_type in skip_results or table_type.startswith('responses.'):
                continue
            # model.displacements
            adict = self.get_result(table_type)
            bdict = op2_model.get_result(table_type)
            if adict is None and bdict is None:
                continue

            # check number of subcases
            if len(adict) != len(bdict):
                self.log.warning('len(self.%s)=%s len(op2_model.%s)=%s' % (
                    table_type, len(adict), table_type, len(bdict)))
                if stop_on_failure:
                    return False
                continue

            # loop over each DisplacementArray
            for key, avalue in adict.items():
                if debug:
                    self.log.debug('working on %r subcase=%s' % (table_type, str(key)))

                # get the displacement for model B
                bvalue = bdict[key]
                is_equal = self._is_op2_case_equal(table_type, key, avalue, bvalue,
                                                   stop_on_failure=stop_on_failure, debug=debug)
                if not is_equal and stop_on_failure:
                    return is_equal
        return True

    def _is_op2_case_equal(self, table_type: str,
                           key, a_obj, b_obj,
                           stop_on_failure: bool=True, debug: bool=False) -> bool:
        """
        Helper method for ``assert_op2_equal``

        Parameters
        ----------
        table_type : str
            the type of table (e.g., ``displacements``)
        key : subcase_id / tuple_obj
            subcase_id : int
                the subcase_id
            tuple_obj : Tuple(???, ???, ...)
                the fancy tuple thingy that you see in single subcase buckling...

                subcase_id : int
                    the subcase_id
                sort_code : int
                    1 : SORT1
                    2 : SORT2
                title??? : str
                    the case title
                subtitle??? : str
                    the case subtitle
                superelement_id : str???
                    the superelement
                other terms???
                TODO: document better
        a_obj : Op2Object()
            a RealDisplacementArray, ComplexDisplacementArray, RealSolidStressArray, etc.
            for the self model
        b_obj : Op2Object()
            a RealDisplacementArray, ComplexDisplacementArray, RealSolidStressArray, etc.
            for the comparison model
        stop_on_failure : bool; default=True
            True : Crashes if they're not equal
            False : Go to the next object
        debug : bool; default=False
            give slightly more debugging messages

        Returns
        -------
        is_equal : bool
            are the objects equal?

        Raises
        ------
        AssertionError/ValueError : stop_on_failure=True and and error occurred
        NotImplementedError : this is a sign of an unsupported object

        """
        # check the name (e.g., RealDisplacementArray vs. ComplexDisplacementArray)
        aname = a_obj.__class__.__name__
        bname = b_obj.__class__.__name__
        if not aname == bname:
            self.log.warning('type(a)=%s type(b)=%s' % (aname, bname))
            return False

        if aname == 'PARAM': # TODO: update this
            return True

        # does this ever hit?
        if not any(word in aname for word in ['Array', 'Eigenvalues', 'GridPointWeight']):
            msg = '%s is not an Array ... assume equal' % aname
            self.log.warning(msg)
            raise NotImplementedError('%s __eq__' % aname)
            #continue

        # use the array methods to check for equality
        # TODO: this can crash
        try:
            is_not_equal = a_obj != b_obj
        except ValueError:
            if stop_on_failure:
                raise

        if is_not_equal:
            self.log.warning(f'key={key} table_type={table_type!r} are not equal; class_name={aname!r}')
            return False
        return True

    def set_mode(self, mode: str) -> None:
        """
        Sets the mode as 'msc', 'nx', 'autodesk', 'nasa95', or 'optistruct'
        """
        if mode.lower() == 'msc':
            self.set_as_msc()
        elif mode.lower() == 'nx':
            self.set_as_nx()
        elif mode.lower() == 'autodesk':
            self.set_as_autodesk()
        elif mode.lower() == 'nasa95':
            self.set_as_nasa95()
        elif mode.lower() == 'optistruct': # radioss,
            self.set_as_optistruct()
        else:
            raise RuntimeError(f'mode={mode!r} and must be in [msc, nx, '
                               f'autodesk, nasa95, optistruct]')

    def to_nx(self) -> None:
        if self.is_msc:
            self.log.warning('switching to NX')
            self.set_as_nx()
            self.set_table_type()

    def to_msc(self) -> None:
        if self.is_nx:
            self.log.warning('switching to MSC')
            self.set_as_msc()

    def include_exclude_results(self,
                                exclude_results: Optional[List[str]]=None,
                                include_results: Optional[List[str]]=None) -> None:
        """
        Sets results to include/exclude

        Parameters
        ----------
        exclude_results / include_results : List[str] / str; default=None
            a list of result types to exclude/include
            one of these must be None

        """
        if exclude_results and include_results:
            msg = (
                'exclude_results or include_results must be None\n'
                'exclude_results=%r\n'
                'include_results=%r\n' % (exclude_results, include_results)
            )
            raise RuntimeError(msg)

        if exclude_results:
            self.remove_results(exclude_results)
        elif include_results:
            self.set_results(include_results)

    def saves(self) -> str:
        """Saves a pickled string"""
        return dumps(self)

    def __getstate__(self):
        """clears out a few variables in order to pickle the object"""
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['log']
        if hasattr(self, 'results') and hasattr(self._results, 'log'):
            del state['_results'].log
        #if hasattr(self, '_card_parser_b'):
            #del state['_card_parser_b']
        #if hasattr(self, '_card_parser_prepare'):
            #del state['_card_parser_prepare']

        # this block let's us identify the objects that are problematic
        # we just play with the value of i to delete all objects past
        # some threshold.  Once we find where the code breaks, we dig
        # into the objects further
        if 0:  # pragma: no cover
            i = 0
            for key, value in sorted(state.items()):
                if isinstance(value, dict) and len(value) == 0:
                    continue
                #if not isinstance(value, (str, int, float)):
                if i > 100:
                    #print('deleting', key)
                    del state[key]
                else:
                    #print('***', key, value)
                    i += 1
                #else:
                    #print(key, type(value), value)
                    #break
                #i += 1
        return state

    def save(self, obj_filename: str='model.obj', unxref: bool=True) -> None:
        """Saves a pickleable object"""
        #del self.log
        #del self._card_parser, self._card_parser_prepare
        if hasattr(self, 'generalized_tables'):
            del self.generalized_tables
        if hasattr(self, 'op2_reader'):
            del self.op2_reader

        #print(object_attributes(self, mode="all", keys_to_skip=[]))
        with open(obj_filename, 'wb') as obj_file:
            dump(self, obj_file)

    def load(self, obj_filename: str='model.obj') -> None:
        """Loads a pickleable object"""
        with open(obj_filename, 'rb') as obj_file:
            obj = load(obj_file)

        keys_to_skip = [
            'ask',
            'binary_debug',
            '_close_op2',
            '_count',
            '_results',
            '_table_mapper',
            'additional_matrices',
            'apply_symmetry',
            'debug_file',
            'expected_times',
            'f',
            'generalized_tables',
            'is_all_subcases',
            'is_debug_file',
            'is_geometry',
            'is_vectorized',
            'isubcase',
            'log',
            'matrix_tables',
            'mode',
            'n',
            'ntotal',
            'num_wide',
            'op2_reader',
            'table_name',
            'use_vector',
            'words',
        ]
        keys = object_attributes(self, mode="all", keys_to_skip=keys_to_skip,
                                 filter_properties=True)
        for key in keys:
            if key.startswith('__') and key.endswith('__'):
                continue

            try:
                val = getattr(obj, key)
            except NameError:
                self.log.warning('key=%r val=%s' % (key, val))
                continue
            #print(key)
            #if isinstance(val, types.FunctionType):
                #continue
            try:
                setattr(self, key, val)
            except AttributeError:
                print('key=%r val=%s' % (key, val))
                raise

        #self.case_control_deck = CaseControlDeck(self.case_control_lines, log=self.log)
        self.log.debug('done loading!')

    @property
    def is_geometry(self) -> bool:
        return False

    def read_op2(self, op2_filename: Optional[str]=None,
                 combine: bool=True,
                 build_dataframe: Optional[bool]=None,
                 skip_undefined_matrices: bool=False,
                 encoding: Optional[str]=None) -> None:
        """
        Starts the OP2 file reading

        Parameters
        ----------
        op2_filename : str (default=None -> popup)
            the op2_filename
        combine : bool; default=True
            True : objects are isubcase based
            False : objects are (isubcase, subtitle) based;
                    will be used for superelements regardless of the option
        #load_as_h5 : default=False
            #loads the op2 as an h5 file to save memory
            #stores the result.element/data attributes in h5 format
        build_dataframe : bool (default=None -> True if in iPython, False otherwise)
            builds a pandas DataFrame for op2 objects
        skip_undefined_matrices : bool; default=False
             True : prevents matrix reading crashes
        encoding : str
            the unicode encoding (default=None; system default)

        """
        if op2_filename:
            check_path(op2_filename, name='op2_filename')
        mode = self.mode
        if build_dataframe is None:
            build_dataframe = False
            if ipython_info():
                build_dataframe = True

        if encoding is None:
            encoding = sys.getdefaultencoding()
        self.encoding = encoding

        self.skip_undefined_matrices = skip_undefined_matrices
        assert self.ask in [True, False], self.ask
        self.is_vectorized = True
        self.log.debug('combine=%s' % combine)
        self.log.debug('-------- reading op2 with read_mode=1 (array sizing) --------')
        self.read_mode = 1
        self._close_op2 = False

        load_as_h5 = False
        if hasattr(self, 'load_as_h5'):
            load_as_h5 = self.load_as_h5

        try:
            # get GUI object names, build objects, but don't read data
            table_names = OP2_Scalar.read_op2(self, op2_filename=op2_filename,
                                              load_as_h5=load_as_h5, mode=mode)
            self.table_names = table_names

            # TODO: stuff to figure out objects
            # TODO: stuff to show gui of table names
            # TODO: clear out objects the user doesn't want
            self.read_mode = 2
            self._close_op2 = True
            self.log.debug('-------- reading op2 with read_mode=2 (array filling) --------')
            _create_hdf5_info(self.op2_reader.h5_file, self)
            OP2_Scalar.read_op2(self, op2_filename=self.op2_filename, mode=mode)
        except FileNotFoundError:
            raise
        except:
            OP2_Scalar.close_op2(self, force=True)
            raise
        self._finalize()
        if build_dataframe:
            self.build_dataframe()
        self.create_objects_from_matrices()
        self.combine_results(combine=combine)
        self.log.debug('finished reading op2')
        str(self.op2_results)

    def create_objects_from_matrices(self) -> None:
        """
        creates the following objects:
          - monitor3 : MONPNT3 object from the MP3F matrix
          - monitor1 : MONPNT1 object from the PMRF, PERF, PFRF, AGRF, PGRF, AFRF matrices

        """
        #assert len(self._frequencies) > 0, self._frequencies
        if 'MP3F' in self.matrices:
            self.monitor3 = MONPNT3(self._frequencies, self.matrices['MP3F'])

        # these are totally wrong...it doesn't go by component;
        # it goes by inertial, external, flexibility, etc.
        if 'PERF' in self.matrices:
            #self.monitor1 = MONPNT1(
                #self._frequencies, self.matrices, [
                # :)       ?       :)      :)      ?       ?
                #'PMRF', 'AFRF', 'PFRF', 'PGRF', 'AGRF', 'PERF', ])

            self.monitor1 = MONPNT1(
                self._frequencies, self.matrices,
                #  :)       ?       :)      :)2     ?       ?
                ['PMRF', 'PERF', 'PFRF', 'AGRF', 'PGRF', 'AFRF', ])

    def _finalize(self) -> None:
        """internal method"""
        if hasattr(self, 'subcase'):
            del self.subcase

        result_types = self.get_table_types()
        for result_type in result_types:
            if result_type in ['params', 'gpdt', 'bgpdt', 'eqexin', 'psds'] or result_type.startswith('responses.'):
                continue
            result = self.get_result(result_type)
            try:
                values = result.values()
            except AttributeError:
                self.log.error('result_type = %s' % result_type)
                raise

            for obj in values:
                if hasattr(obj, 'finalize'):
                    obj.finalize()
                elif hasattr(obj, 'tCode') and not obj.is_sort1:
                    raise RuntimeError('object has not implemented finalize\n%s' % (
                        ''.join(obj.get_stats())))
        self.del_structs()

    def build_dataframe(self) -> None:
        """
        Converts the OP2 objects into pandas DataFrames

        .. todo:: fix issues with:
         - RealDisplacementArray
         - RealPlateStressArray (???)
         - RealPlateStrainArray (???)
         - RealCompositePlateStrainArray (???)

        """
        # TODO: sorter = uniques.argsort()
        #C:\Anaconda\lib\site-packages\pandas\core\algorithms.py:198:
        #    DeprecationWarning: unorderable dtypes;
        #    returning scalar but in the future this will be an error

        no_sort2_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
        result_types = self.get_table_types()

        if len(self.matrices):
            for key, matrix in sorted(self.matrices.items()):
                if hasattr(matrix, 'build_dataframe'):
                    matrix.build_dataframe()
                else:
                    self.log.warning('pandas: build_dataframe is not supported for key=%s type=%s' % (
                        key, str(type(matrix))))
                    raise NotImplementedError()
                    #continue

        skip_pandas = ['params', 'gpdt', 'bgpdt', 'eqexin', 'grid_point_weight', 'psds']
        for result_type in result_types:
            if result_type in skip_pandas or result_type.startswith('responses.'):
                #self.log.debug('skipping %s' % result_type)
                continue

            result = self.get_result(result_type)
            for obj in result.values():
                class_name = obj.__class__.__name__
                #print('working on %s' % class_name)
                obj.object_attributes()
                obj.object_methods()
                str(obj)
                obj.get_stats()

                if class_name in no_sort2_classes:
                    try:
                        obj.build_dataframe()
                        assert obj.data_frame is not None
                    except MemoryError:
                        raise
                    except:
                        self.log.error(obj)
                        self.log.error('build_dataframe is broken for %s' % class_name)
                        raise
                    continue
                if obj.is_sort2:
                    #self.log.warning(obj)
                    self.log.warning('build_dataframe is not supported for %s - SORT2' % class_name)
                    continue

                # SORT1
                try:
                    obj.build_dataframe()
                #except TypeError:
                    #self.log.error(obj)
                    #self.log.error('build_dataframe is broken with a TypeError for %s' % class_name)
                except MemoryError:
                    raise
                except NotImplementedError:
                    self.log.warning(obj)
                    self.log.warning('build_dataframe is broken for %s' % class_name)
                    raise
                except:
                    self.log.error(obj)
                    self.log.error('build_dataframe is broken for %s' % class_name)
                    raise

    def load_hdf5(self, hdf5_filename: str, combine: bool=True) -> None:  # praga: no cover
        """Loads an h5 file into an OP2 object"""
        self.deprecated('load_hdf5', 'load_hdf5_filename', '1.2')
        return self.load_hdf5_filename(hdf5_filename, combine=True)

    def load_hdf5_filename(self, hdf5_filename: str, combine: bool=True) -> None:
        """
        Loads an h5 file into an OP2 object

        Parameters
        ----------
        hdf5_filename : str
            the path to the an hdf5 file
        combine : bool; default=True
            runs the combine routine

        """
        check_path(hdf5_filename, 'hdf5_filename')
        from pyNastran.op2.op2_interface.hdf5_interface import load_op2_from_hdf5_file
        import h5py
        self.op2_filename = hdf5_filename

        self.log.info('hdf5_op2_filename = %r' % hdf5_filename)
        debug = False
        with h5py.File(hdf5_filename, 'r') as h5_file:
            load_op2_from_hdf5_file(self, h5_file, self.log, debug=debug)
        self.combine_results(combine=combine)

    def load_hdf5_file(self, h5_file: H5File, combine: bool=True) -> None:
        """
        Loads an h5 file object into an OP2 object

        Parameters
        ----------
        h5_file : H5File()
            an h5py file object
        combine : bool; default=True
            runs the combine routine

        """
        from pyNastran.op2.op2_interface.hdf5_interface import load_op2_from_hdf5_file
        #self.op2_filename = hdf5_filename
        #self.log.info('hdf5_op2_filename = %r' % hdf5_filename)
        debug = False
        load_op2_from_hdf5_file(self, h5_file, self.log, debug=debug)
        self.combine_results(combine=combine)

    def export_hdf5(self, hdf5_filename: str) -> None:
        """Converts the OP2 objects into hdf5 object"""
        self.deprecated('export_hdf5', 'export_hdf5_filename', '1.2')
        return self.export_hdf5_filename(hdf5_filename)

    def export_hdf5_filename(self, hdf5_filename: str) -> None:
        """
        Converts the OP2 objects into hdf5 object

        TODO: doesn't support:
          - BucklingEigenvalues

        """
        from pyNastran.op2.op2_interface.hdf5_interface import export_op2_to_hdf5_filename
        export_op2_to_hdf5_filename(hdf5_filename, self)

    def export_hdf5_file(self, hdf5_file: H5File, exporter=None) -> None:
        """
        Converts the OP2 objects into hdf5 object

        Parameters
        ----------
        hdf5_file : H5File()
            an h5py object
        exporter : HDF5Exporter; default=None
            unused

        TODO: doesn't support:
          - BucklingEigenvalues

        """
        ## type (file, Any) -> None
        from pyNastran.op2.op2_interface.hdf5_interface import export_op2_to_hdf5_file
        export_op2_to_hdf5_file(hdf5_file, self)

    def combine_results(self, combine: str=True) -> None:
        """
        we want the data to be in the same format and grouped by subcase, so
        we take

        .. code-block:: python

           stress = {
               # isubcase, analysis_code, sort_method, count, superelement_adaptivity_index, pval_step
               (1, 2, 1, 0, 'SUPERELEMENT 0', '') : result1,
               (1, 2, 1, 0, 'SUPERELEMENT 10', '') : result2,
               (1, 2, 1, 0, 'SUPERELEMENT 20', '') : result3,
               (2, 2, 1, 0, 'SUPERELEMENT 0', '') : result4,

           code = (isubcase, analysis_code, sort_method, count, ogs,
                   superelement_adaptivity_index, pval_step)
             }
        and convert it to:

        .. code-block:: python

           stress = {
               1 : result1 + result2 + results3,
               2 : result4,
           }

        """
        self.combine = combine
        result_types = self.get_table_types()
        results_to_skip = ['bgpdt', 'gpdt', 'eqexin', 'grid_point_weight', 'psds']

        # set subcase_key
        for result_type in result_types:
            if result_type in results_to_skip or result_type.startswith('responses.'):
                continue
            result = self.get_result(result_type)
            try:
                case_keys = sorted(result.keys())
            except AttributeError:
                self.log.error(f'result_type = {result_type}')
                raise
            # unique_isubcases = []  # List[int]
            for case_key in case_keys:
                #print('case_key =', case_key)
                if isinstance(case_key, tuple):
                    isubcasei, analysis_codei, sort_methodi, counti, isuperelmemnt_adaptivity_index, pval_step, ogs = case_key
                    #isubcasei, analysis_codei, sort_methodi, counti, isuperelmemnt_adaptivity_index, table_name = case_key
                    if ogs == 0:
                        value = (analysis_codei, sort_methodi, counti,
                                 isuperelmemnt_adaptivity_index, pval_step)
                    else:
                        value = (analysis_codei, sort_methodi, counti,
                                 isuperelmemnt_adaptivity_index, pval_step, ogs)

                    if value not in self.subcase_key[isubcasei]:
                        #print('isubcase=%s value=%s' % (isubcasei, value))
                        self.subcase_key[isubcasei].append(value)
                else:
                    #print('combine - case_key =', case_keys)
                    break

        if not combine:
            subcase_key2 = {}
            for key, value in self.subcase_key.items():
                subcase_key2[key] = value
            self.subcase_key = subcase_key2
            #print('self.subcase_key =', self.subcase_key)
            #print('skipping combine results')
            return
        del result, case_keys

        isubcases = np.unique(list(self.subcase_key.keys()))
        unique_isubcases = np.unique(isubcases)

        self.log.debug('combine_results')
        for result_type in result_types:
            if result_type in results_to_skip or result_type.startswith('responses.'):
                continue
            result = self.get_result(result_type)
            if len(result) == 0:
                continue
            for isubcase in unique_isubcases:
                try:
                    keys = self.subcase_key[isubcase]
                except TypeError:
                    print('isubcase =', isubcase)
                    print('isubcases =', isubcases)
                    print('self.subcase_key =', self.subcase_key)
                    raise

                #print('keys = %s' % keys)
                key0 = tuple([isubcase] + list(keys[0]))

                if len(key0) == 5:
                    # ogs is optional
                    isubcase, analysis_code, unused_sort_code, count, isuperelmemnt_adaptivity_index, pval_step = key0
                    key1 = (isubcase, analysis_code, 1, count, isuperelmemnt_adaptivity_index, pval_step)
                    key2 = (isubcase, analysis_code, 2, count, isuperelmemnt_adaptivity_index, pval_step)
                else:
                    isubcase, analysis_code, unused_sort_code, count, isuperelmemnt_adaptivity_index, pval_step, ogs = key0
                    key1 = (isubcase, analysis_code, 1, count, isuperelmemnt_adaptivity_index, pval_step, ogs)
                    key2 = (isubcase, analysis_code, 2, count, isuperelmemnt_adaptivity_index, pval_step, ogs)

                #isubcase, analysis_code, sort_code, count, isuperelmemnt_adaptivity_index, table_name = key0
                #key1 = (isubcase, analysis_code, 1, count, isuperelmemnt_adaptivity_index, table_name)
                #key2 = (isubcase, analysis_code, 2, count, isuperelmemnt_adaptivity_index, table_name)
                if len(keys) == 1:
                    if key0 not in result:
                        continue
                    # rename the case since we have only one tuple for the result
                    # key0 = tuple([isubcase] + list(key0))
                    result[isubcase] = result[key0]
                    del result[key0]
                elif len(keys) == 2 and key1 in keys and key2 in keys:
                    # continue
                    #print('key0 =', result_type, key0)
                    # res0 = result[key0]

                    isubcase, analysis_code, unused_sort_code, count, isuperelmemnt_adaptivity_index = key0
                    #isubcase, analysis_code, sort_code, count, isuperelmemnt_adaptivity_index, table_name = key0
                    if not (key1 in result and key2 in result):
                        if key1 in result:
                            res1 = result[key1]
                            self.log.info("res=%s has a single case; trivial" %
                                          res1.__class__.__name__)
                            result[isubcase] = result[key1]
                            #print('del key1=%s' % str(key1))
                            del result[key1]
                        elif key2 in result:
                            res2 = result[key2]
                            self.log.info("res=%s has a single case; trivial" %
                                          res2.__class__.__name__)
                            result[isubcase] = result[key2]
                            #print('del key2=%s' % str(key2))
                            del result[key2]
                        continue

                    res1 = result[key1]
                    if not hasattr(res1, 'combine'):
                        self.log.info("res=%s has no method combine" % res1.__class__.__name__)
                        continue

                    self.log.info("res=%s has combine" % res1.__class__.__name__)
                    res2 = result[key2]
                    del result[key1]
                    del result[key2]
                    res1.combine(res2)
                    result[isubcase] = res1
                     #print('r[isubcase] =', result[isubcase])
                else:
                    #self.log.info("continue")
                    continue
            setattr(self, result_type, result)
        #print('subcase_key =', self.subcase_key)

        subcase_key2 = {}
        not_results = ['eigenvalues', 'eigenvalues_fluid', 'params', 'gpdt', 'bgpdt',
                       'eqexin', 'desvars', 'grid_point_weight', 'psds']
        for result_type in result_types:
            if result_type in not_results or result_type.startswith('responses.'):
                continue
            result = self.get_result(result_type)
            try:
                case_keys = list(result.keys())
            except AttributeError:
                self.log.error('result_type = %s' % result_type)
                raise

            try:
                case_keys = sorted(case_keys)  # TODO: causes DeprecationWarning
            except TypeError:
                self.log.error('result.keys() = %s' % case_keys)

            if len(result) == 0:
                continue
            for isubcase in unique_isubcases:
                if isubcase not in subcase_key2:
                    subcase_key2[isubcase] = []

            for isubcase in unique_isubcases:
                for case_key in case_keys:
                    #print('isubcase=%s case_key=%s' % (isubcase, case_key))
                    assert not isinstance(case_key, str), result_type
                    if isinstance(case_key, integer_types):
                        if isubcase == case_key and case_key not in subcase_key2[isubcase]:
                            subcase_key2[isubcase] = [isubcase]
                    else:
                        try:
                            subcasei = case_key[0]
                        except IndexError:
                            msg = 'case_key=%s; type(case_key)=%s; case_key[0] is not the subcase id' % (
                                case_key, type(case_key))
                            raise IndexError(msg)
                        #if not subcasei == isubcase:
                            #continue
                        if case_key not in subcase_key2[subcasei]:
                            subcase_key2[isubcase].append(case_key)
        self.subcase_key = subcase_key2
        #print('subcase_key = %s' % self.subcase_key)

    def get_key_order(self) -> List[Tuple[int, int, int, int, int, str]]:
        """
        Returns
        -------
        keys3 : List[int, int, int, int, int, str]
            the keys in order
        """
        keys = []
        table_types = self.get_table_types()
        skip_tables = ['gpdt', 'bgpdt', 'eqexin', 'grid_point_weight', 'psds']
        for table_type in sorted(table_types):
            if table_type in skip_tables or table_type.startswith('responses.'):
                continue
            result_type_dict = self.get_result(table_type)
            #if result_type_dict is None: # gpdt, eqexin
                #continue
            if len(result_type_dict) == 0:
                continue
            for key in result_type_dict:
                if isinstance(key, str):
                    if table_type not in ['eigenvalues', 'eigenvalues_fluid', 'params']:
                        self.log.warning('table_type = %s' % table_type)
                    continue
                if key not in keys:
                    keys.append(key)

        #print(self.get_op2_stats())
        #keys_order = []

        # subcase_ids = self.subcase_key.keys()

        #self.isubcase_name_map[self.isubcase] = [self.subtitle, self.analysis_code, self.label]
        #subcase_ids = list(self.isubcase_name_map.keys())
        #subcase_ids.sort()
        #print('subcase_ids =', subcase_ids)


        # isubcase, analysis_code, sort_method, count, ogs, superelement_adaptivity_index, pval_step
        #(1, 2, 1, 0, 0, 'SUPERELEMENT 0')  : result1

        isubcases = set()
        analysis_codes = set()
        sort_methods = set()
        counts = set()
        ogss = set()
        superelement_adaptivity_indexs = set()
        pval_steps = set()

        for key in keys:
            #print('key = %s' % str(key))
            if len(key) == 6:
                isubcase, analysis_code, sort_method, count, superelement_adaptivity_index, pval_step = key
                ogs = 0
            elif len(key) == 7:
                isubcase, analysis_code, sort_method, count, ogs, superelement_adaptivity_index, pval_step = key
            else:
                print('  %s' % str(key))
                raise RuntimeError(key)
                #isubcase, analysis_code, sort_method, count, ogs, superelement_adaptivity_index = key
            isubcases.add(isubcase)
            analysis_codes.add(analysis_code)
            sort_methods.add(sort_method)
            counts.add(count)
            ogss.add(ogs)
            superelement_adaptivity_indexs.add(superelement_adaptivity_index)
            pval_steps.add(pval_step)

        isubcase_list = list(isubcases)
        analysis_code_list = list(analysis_codes)
        sort_method_list = list(sort_methods)
        count_list = list(counts)
        ogs_list = list(ogss)
        superelement_adaptivity_index_list = list(superelement_adaptivity_indexs)
        pval_step_list = list(pval_steps)

        isubcase_list.sort()
        analysis_code_list.sort()
        sort_method_list.sort()
        count_list.sort()
        ogs_list.sort()
        superelement_adaptivity_index_list.sort()
        pval_step_list.sort()

        keys3 = []  # type: List[Tuple[int, int, int, int, int, int, str]]
        for isubcase in isubcase_list:
            for count in count_list:
                for analysis_code in analysis_code_list:
                    for superelement_adaptivity_index in superelement_adaptivity_index_list:
                        for pval_step in pval_step_list:
                            for sort_method in sort_method_list:
                                for ogs in ogs_list:
                                    key = (isubcase, analysis_code, sort_method,
                                           count, ogs, superelement_adaptivity_index, pval_step)
                                    if key not in keys3:
                                        #print('adding ', key)
                                        keys3.append(key)
        if len(keys3) == 0:
            self.log.warning('No results...\n' + self.get_op2_stats(short=True))
        #assert len(keys3) > 0, keys3
        return keys3

    def print_subcase_key(self) -> None:
        self.log.info('---self.subcase_key---')
        for isubcase, keys in sorted(self.subcase_key.items()):
            if len(keys) == 1:
                self.log.info('subcase_id=%s : keys=%s' % (isubcase, keys))
            else:
                self.log.info('subcase_id=%s' % isubcase)
                for key in keys:
                    self.log.info('  %s' % str(key))
        #self.log.info('subcase_key = %s' % self.subcase_key)

    def transform_displacements_to_global(self, icd_transform: Any, coords: Dict[int, Any], xyz_cid0: Any=None,
                                          debug: bool=False) -> None:
        """
        Transforms the ``data`` of displacement-like results into the
        global coordinate system for those nodes with different output
        coordinate systems. Takes indicies and transformation matricies
        for nodes with their output in coordinate systems other than the
        global.

        Used in combination with  ``BDF.get_displacement_index``

        Parameters
        ----------
        icd_transform : dict{int cid : int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``BDF.point_ids`` that their output (`CD`) in that
            coordinate system.
        coords : dict{int cid :Coord()}
            Dictionary of coordinate id to the coordinate object
            Use this if CD is only rectangular
            Use this if CD is not rectangular
        xyz_cid0 : (nnodes+nspoints, 3) float ndarray
            the nodes in the global frame
            Don't use this if CD is only rectangular
            Use this if CD is not rectangular
        debug : bool; default=False
            developer debug

        .. warning:: only works if all nodes are included...
                     ``test_pynastrangui isat_tran.dat isat_tran.op2 -f nastran``
        .. note:: Nastran has this concept of a basic (cid=0) and global (cid=cd)
                  coordinate system.  They occur at the same time.  Basic is for
                  positions/properties, while global is for result outputs.

                  pyNastran's OP2 interface uses:
                    - cd=0 for global frames
                    - cd>0 are local frames
                  pyNastran's BDF interface uses:
                    - cp=0 for global frames
                    - cp>0 are local frames

        """
        #output = {}
        ato = self.op2_results.ato
        crm = self.op2_results.crm
        psd = self.op2_results.psd
        rms = self.op2_results.rms
        #no = self.op2_results.no
        disp_like_dicts = [
            # should NO results be transformed?
            #no.displacements, no.velocities, no.accelerations,
            #no.spc_forces, no.mpc_forces,

            self.displacements,
            ato.displacements, crm.displacements, psd.displacements, rms.displacements,
            self.displacements_scaled,
            self.displacement_scaled_response_spectra_abs,
            self.displacement_scaled_response_spectra_nrl,

            self.velocities,
            ato.velocities, crm.velocities, psd.velocities, rms.velocities,
            self.velocity_scaled_response_spectra_abs,

            self.accelerations,
            ato.accelerations, crm.accelerations, psd.accelerations, rms.accelerations,
            self.acceleration_scaled_response_spectra_abs,
            self.acceleration_scaled_response_spectra_nrl,

            self.eigenvectors,
            self.op2_results.RADCONS.eigenvectors, self.op2_results.RADEFFM.eigenvectors,
            self.op2_results.RADEATC.eigenvectors, self.op2_results.ROUGV1.eigenvectors,

            self.spc_forces, ato.spc_forces, crm.spc_forces, psd.spc_forces, rms.spc_forces,
            self.mpc_forces, ato.mpc_forces, crm.mpc_forces, psd.mpc_forces, rms.mpc_forces,

            self.applied_loads,
            self.load_vectors,
        ]
        for disp_like_dict in disp_like_dicts:
            if not disp_like_dict:
                continue
            #print('-----------')
            for subcase, result in disp_like_dict.items():
                if result.table_name in ['BOUGV1', 'BOPHIG', 'TOUGV1']:
                    continue
                self.log.debug("transforming %s" % result.table_name)
                transform_displacement_to_global(subcase, result, icd_transform, coords, xyz_cid0,
                                                 self.log, debug=debug)

    def transform_gpforce_to_global(self, nids_all, nids_transform, icd_transform, coords, xyz_cid0=None):
        """
        Transforms the ``data`` of GPFORCE results into the
        global coordinate system for those nodes with different output
        coordinate systems. Takes indicies and transformation matricies
        for nodes with their output in coordinate systems other than the
        global.

        Used in combination with ``BDF.get_displacement_index``

        Parameters
        ----------
        nids_all : ???
            ???
        nids_transform : dict{int cid : int ndarray nds}
            Dictionary from coordinate id to corresponding node ids.
        icd_transform : dict{int cid : int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``BDF.point_ids`` that their output (`CD`) in that
            coordinate system.
        coords : dict{int cid :Coord()}
            Dictionary of coordinate id to the coordinate object
            Use this if CD is only rectangular
            Use this if CD is not rectangular
        xyz_cid0 : ???
            required for cylindrical/spherical coordinate systems

        """
        disp_like_dicts = [
            # TODO: causes test_op2_solid_shell_bar_01_gpforce_xyz to fail
            #       even though it should be uncommented
            self.grid_point_forces,
        ]
        for disp_like_dict in disp_like_dicts:
            if not disp_like_dict:
                continue
            self.log.debug('-----------')
            for subcase, result in disp_like_dict.items():
                transform_gpforce_to_globali(subcase, result,
                                             nids_all, nids_transform,
                                             icd_transform, coords, xyz_cid0, self.log)
        self.log.debug('-----------')


def read_op2(op2_filename: Optional[str]=None,
             combine: bool=True,
             subcases: Optional[List[int]]=None,
             exclude_results: Optional[List[str]]=None,
             include_results: Optional[List[str]]=None,
             log: Any=None,
             debug: bool=True,
             debug_file: Optional[str]=None,
             build_dataframe: Optional[bool]=None,
             skip_undefined_matrices: bool=True,
             mode: Optional[str]=None,
             encoding: Optional[str]=None) -> OP2:
    """
    Creates the OP2 object without calling the OP2 class.

    Parameters
    ----------
    op2_filename : str (default=None -> popup)
        the op2_filename
    combine : bool; default=True
        True : objects are isubcase based
        False : objects are (isubcase, subtitle) based;
                will be used for superelements regardless of the option
    subcases : List[int, ...] / int; default=None->all subcases
        list of [subcase1_ID,subcase2_ID]
    exclude_results / include_results : List[str] / str; default=None
        a list of result types to exclude/include
        one of these must be None
    build_dataframe : bool (default=None -> True if in iPython, False otherwise)
        builds a pandas DataFrame for op2 objects
    skip_undefined_matrices : bool; default=False
         True : prevents matrix reading crashes
    debug : bool; default=False
        enables the debug log and sets the debug in the logger
    log : Log()
        a logging object to write debug messages to
        (.. seealso:: import logging)
    mode : str; default=None -> 'msc'
        the version of the Nastran you're using
        {nx, msc, autodesk, optistruct}
    debug_file : str; default=None (No debug)
        sets the filename that will be written to
    encoding : str
        the unicode encoding (default=None; system default)

    Returns
    -------
    model : OP2()
        an OP2 object

    .. todo:: creates the OP2 object without all the read methods

    .. note :: this method will change in order to return an object that
               does not have so many methods

    """
    if op2_filename:
        check_path(op2_filename, name='op2_filename')
    model = OP2(log=log, debug=debug, debug_file=debug_file, mode=mode)
    model.set_subcases(subcases)
    model.include_exclude_results(exclude_results=exclude_results,
                                  include_results=include_results)

    model.read_op2(op2_filename=op2_filename, build_dataframe=build_dataframe,
                   skip_undefined_matrices=skip_undefined_matrices, combine=combine,
                   encoding=encoding)
    ## TODO: this will go away when OP2 is refactored
    ## TODO: many methods will be missing, but it's a start...
    ## doesn't support F06 writer
    #obj = Op2F06Attributes()
    #attr_names = object_attributes(obj, mode="public", keys_to_skip=None)
    #for attr_name in attr_names:
        #attr = getattr(model, attr_name)
        #setattr(obj, attr_name, attr)
    #obj.get_op2_stats()
    return model


def _create_hdf5_info(h5_file, op2_model):
    """exports the h5 info group"""
    load_as_h5 = False
    if hasattr(op2_model, 'load_as_h5'):
        load_as_h5 = op2_model.load_as_h5
    if not load_as_h5:
        return
    from pyNastran.op2.op2_interface.hdf5_interface import create_info_group
    create_info_group(h5_file, op2_model)
