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
   - object_attributes(mode='public', keys_to_skip=None)
   - object_methods(mode='public', keys_to_skip=None)
   - print_subcase_key()
   - read_op2(op2_filename=None, combine=True, build_dataframe=None,
              skip_undefined_matrices=False, encoding=None)
   - set_mode(mode)
   - transform_displacements_to_global(i_transform, coords, xyz_cid0=None, debug=False)
   - transform_gpforce_to_global(nids_all, nids_transform, i_transform, coords, xyz_cid0=None)

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
from copy import deepcopy
from six import iterkeys, iteritems, string_types, itervalues, b
from six.moves.cPickle import load, dump

import numpy as np

import pyNastran
from pyNastran.utils import (
    object_attributes, object_methods, integer_types, ipython_info)
from pyNastran.op2.tables.monpnt import MONPNT1, MONPNT3

from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import SortCodeError, DeviceCodeError, FortranMarkerError
#from pyNastran.op2.op2_interface.op2_writer import OP2Writer
from pyNastran.op2.op2_interface.op2_f06_common import Op2F06Attributes
from pyNastran.op2.op2_interface.op2_scalar import OP2_Scalar


def read_op2(op2_filename=None, combine=True, subcases=None,
             exclude_results=None, include_results=None,
             log=None, debug=True, debug_file=None, build_dataframe=None,
             skip_undefined_matrices=True, mode='msc', encoding=None):
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
    mode : str; default='msc'
        the version of the Nastran you're using
        {nx, msc, optistruct}
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
    model = OP2(log=log, debug=debug, debug_file=debug_file, mode=mode)
    model.set_subcases(subcases)
    if exclude_results and include_results:
        msg = (
            'exclude_results or include_results must be None\n'
            'exclude_results=%r\n'
            'include_results=%r\n' % (exclude_results, include_results)
        )
        raise RuntimeError(msg)
    elif exclude_results:
        model.remove_results(exclude_results)
    elif include_results:
        model.set_results(include_results)

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


#class OP2(OP2_Scalar, OP2Writer):
class OP2(OP2_Scalar):

    def __init__(self,
                 debug=True, log=None,
                 debug_file=None, mode='msc'):
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
        mode : str; default='msc'
            {msc, nx}

        """
        self.encoding = None
        self.set_mode(mode)
        make_geom = False
        assert make_geom == False, make_geom
        OP2_Scalar.__init__(self, debug=debug, log=log, debug_file=debug_file)
        self.ask = False

    def object_attributes(self, mode='public', keys_to_skip=None):
        # type: (str, Optional[List[str]]) -> List[str]
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
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        # type: (str, Optional[List[str]]) -> List[str]
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

    def __eq__(self, op2_model):
        """
        Diffs the current op2 model vs. another op2 model.
        Crashes if they're not equal.
        """
        try:
            is_equal = self.assert_op2_equal(op2_model,
                                             stop_on_failure=True, debug=False)
        except (AssertionError, ValueError):
            is_equal = False
        return is_equal

    def assert_op2_equal(self, op2_model, stop_on_failure=True, debug=False):
        """
        Diffs the current op2 model vs. another op2 model.

        Parameters
        ----------
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
        if not self.read_mode == op2_model.read_mode:
            self.log.warning('self.read_mode=%s op2_model.read_mode=%s ... assume True' % (
                self.read_mode, op2_model.read_mode))
            return True

        table_types = self.get_table_types()
        for table_type in table_types:
            # model.displacements
            adict = getattr(self, table_type)
            bdict = getattr(op2_model, table_type)

            # check number of subcases
            if len(adict) != len(bdict):
                self.log.warning('len(self.%s)=%s len(op2_model.%s)=%s' % (
                    table_type, len(adict), table_type, len(bdict)))
                if stop_on_failure:
                    return False
                continue

            # loop over each DisplacementArray
            for key, avalue in iteritems(adict):
                if debug:
                    self.log.debug('working on %r subcase=%s' % (table_type, str(key)))

                # get the displacement for model B
                bvalue = bdict[key]
                is_equal = self._is_op2_case_equal(table_type, key, avalue, bvalue,
                                                   stop_on_failure=stop_on_failure, debug=debug)
                if not is_equal and stop_on_failure:
                    return is_equal
        return True

    def _is_op2_case_equal(self, table_type, key, a_obj, b_obj, stop_on_failure=True, debug=False):
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

        # does this ever hit?
        if not any(word in aname for word in ['Array', 'Eigenvalues']):
            msg = '%s is not an Array ... assume equal' % aname
            self.log.warning(msg)
            raise NotImplementedError('%s __eq__' % aname)
            #continue

        # use the array methods to check for equality
        # TODO: this can crash
        if a_obj != b_obj:
            self.log.warning('key=%s table_type=%r are not equal; class_name=%r' % (
                key, table_type, aname))
            return False
        return True

    def set_mode(self, mode):
        """
        Sets the mode as 'msc' or 'nx'
        """
        if mode.lower() == 'msc':
            self.set_as_msc()
        elif mode.lower() == 'nx':
            self.set_as_nx()
        else:
            raise RuntimeError("mode=%r and must be 'msc' or 'nx'")

    def saves(self):
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

    def save(self, obj_filename='model.obj', unxref=True):
        # type: (str, bool) -> None
        """Saves a pickleable object"""
        #del self.log
        #del self._card_parser, self._card_parser_prepare

        #print(object_attributes(self, mode="all", keys_to_skip=[]))
        with open(obj_filename, 'wb') as obj_file:
            dump(self, obj_file)

    def load(self, obj_filename='model.obj'):
        # type: (str) -> None
        """Loads a pickleable object"""
        with open(obj_filename, 'rb') as obj_file:
            obj = load(obj_file)

        keys_to_skip = [
            'total_effective_mass_matrix',
            'effective_mass_matrix',
            'rigid_body_mass_matrix',
            'modal_effective_mass_fraction',
            'modal_participation_factors',
            'modal_effective_mass',
            'modal_effective_weight',
        ]
        for key in object_attributes(self, mode="all", keys_to_skip=keys_to_skip):
            if key.startswith('__') and key.endswith('__'):
                continue

            val = getattr(obj, key)
            print(key)
            #if isinstance(val, types.FunctionType):
                #continue
            try:
                setattr(self, key, val)
            except AttributeError:
                print('key=%r val=%s' % (key, val))
                raise

        #self.case_control_deck = CaseControlDeck(self.case_control_lines, log=self.log)
        self.log.debug('done loading!')

    #def _set_ask_vectorized(self, ask=False):
        #"""
        #Enables vectorization

        #The code will degenerate to dictionary based results when
        #a result does not support vectorization.

        #Vectorization is always True here.

        #Parameters
        #----------
        #ask: bool
            #Do you want to see a GUI of result types.

        #+--------+---------------+---------+------------+
        #| Case # | Vectorization |  Ask    | Read Modes |
        #+========+===============+=========+============+
        #|    1   | True          |  True   |  1, 2      |
        #+--------+---------------+---------+------------+
        #|    2   | True          |  False  |  1, 2      |
        #+--------+---------------+---------+------------+
        #|    3   | False         |  True   |  1, 2      |
        #+--------+---------------+---------+------------+
        #|    4   | False         |  False  |  0         |
        #+--------+---------------+---------+------------+

        #Definitions
        #===========
          #Vectorization - A storage structure that allows for faster read/access
                          #speeds and better memory usage, but comes with a more
                          #difficult to use data structure.

                          #It limits the node IDs to all be integers (e.g. element
                          #centroid).  Composite plate elements (even for just CTRIA3s)
                          #with an inconsistent number of layers will have a more
                          #difficult data structure.
          #Scanning   - a quick check used to figure out how many results to process
                      #that takes almost no time
          #Reading    - process the op2 data
          #Build      - call the __init__ on a results object (e.g. RealDisplacementArray)
          #Start Over - Go to the start of the op2 file
          #Ask        - launch a GUI dialog to let the user click which results to load

        #Read Mode Definitions
        #=====================
          #0.   The default OP2 dictionary based-approach with no asking GUI (removed)
          #1.   The first read of a result to get the shape of the data
          #2.   The second read of a result to get the results

        #Cases
        #======
          #1.   Scan the block to get the size, build the object (read_mode=1),
               #ask the user, start over, fill the objects (read_mode=2).
               #Degenerate to read_mode=0 when read_mode=2 cannot be used based
               #upon the value of ask.
          #2.   Same as case #1, but don't ask the user.
               #Scan the block to get the size, build the object (read_mode=1),
               #start over, fill the objects (read_mode=2).
          #3.   Scan the block to get the object types (read_mode=1), ask the user,
               #build the object & fill it (read_mode=2)
          #4.   Read the block to get the size, build the object & fill it (read_mode=0; removed)
        #"""
        #self.ask = ask

    def read_op2(self, op2_filename=None, combine=True, build_dataframe=None,
                 skip_undefined_matrices=False, encoding=None):
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
        build_dataframe : bool (default=None -> True if in iPython, False otherwise)
            builds a pandas DataFrame for op2 objects
        skip_undefined_matrices : bool; default=False
             True : prevents matrix reading crashes
        encoding : str
            the unicode encoding (default=None; system default)

        """
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

        # get GUI object names, build objects, but don't read data
        OP2_Scalar.read_op2(self, op2_filename=op2_filename)

        # TODO: stuff to figure out objects
        # TODO: stuff to show gui of table names
        # TODO: clear out objects the user doesn't want
        self.read_mode = 2
        self._close_op2 = True
        self.log.debug('-------- reading op2 with read_mode=2 (array filling) --------')
        OP2_Scalar.read_op2(self, op2_filename=self.op2_filename)

        self._finalize()
        if build_dataframe:
            self.build_dataframe()
        self.create_objects_from_matrices()
        self.combine_results(combine=combine)
        self.log.debug('finished reading op2')

    def create_objects_from_matrices(self):
        """
        creates the following objects:
          - sonitor3 : MONPNT3 object from the MP3F matrix
          - monitor1 : MONPNT1 object from the PMRF, PERF, PFRF, AGRF, PGRF, AFRF matrices

        """
        if 'MP3F' in self.matrices:
            self.monitor3 = MONPNT3(self._frequencies, self.matrices['MP3F'])

        # these are totally wrong...it doesn't go by component; it goes by inertial, external, flexibility, etc.
        if 'PERF' in self.matrices:
            #                                                           :)       ?       :)      :)      ?       ?
            #self.monitor1 = MONPNT1(self.frequencies, self.matrices, ['PMRF', 'AFRF', 'PFRF', 'PGRF', 'AGRF', 'PERF', ])

            #                                                           :)       ?       :)      :)2     ?       ?
            self.monitor1 = MONPNT1(self._frequencies, self.matrices, ['PMRF', 'PERF', 'PFRF', 'AGRF', 'PGRF', 'AFRF', ])

    def _finalize(self):
        """internal method"""
        result_types = self.get_table_types()
        for result_type in result_types:
            result = getattr(self, result_type)
            for obj in itervalues(result):
                if hasattr(obj, 'finalize'):
                    obj.finalize()
        self.del_structs()

    def build_dataframe(self):
        """
        Converts the OP2 objects into pandas DataFrames

        .. todo:: fix issues with:
         - RealDisplacementArray
         - RealPlateStressArray (???)
         - RealPlateStrainArray (???)
         - RealCompositePlateStrainArray (???)

        """
        # TODO: sorter = uniques.argsort()
        #C:\Anaconda\lib\site-packages\pandas\core\algorithms.py:198: DeprecationWarning: unorderable dtypes;
            #returning scalar but in the future this will be an error

        no_sort2_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
        result_types = self.get_table_types()

        if len(self.matrices):
            for key, matrix in sorted(iteritems(self.matrices)):
                if hasattr(matrix, 'build_dataframe'):
                    matrix.build_dataframe()
                else:
                    self.log.warning('pandas: build_dataframe is not supported for key=%s type=%s' % (key, str(type(matrix))))
                    raise NotImplementedError()
                    #continue

        for result_type in result_types:
            result = getattr(self, result_type)
            for obj in itervalues(result):
                class_name = obj.__class__.__name__
                #print('working on %s' % class_name)
                obj.object_attributes()
                obj.object_methods()

                if class_name in no_sort2_classes:
                    try:
                        obj.build_dataframe()
                        obj.object_methods()
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
                try:
                    obj.build_dataframe()
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

    def load_hdf5(self, hdf5_filename, combine=True):
        """loads an h5 file into an OP2 object"""
        assert os.path.exists(hdf5_filename), print_bad_path(hdf5_filename)
        from pyNastran.op2.op2_interface.hdf5_interface import load_op2_from_hdf5_file
        import h5py
        self.op2_filename = hdf5_filename

        self.log.info('hdf5_op2_filename = %r' % hdf5_filename)
        debug = False
        with h5py.File(hdf5_filename, 'r') as h5_file:
            load_op2_from_hdf5_file(self, h5_file, self.log, debug=debug)
        self.combine_results(combine=combine)

    def export_to_hdf5(self, hdf5_filename):
        """
        Converts the OP2 objects into hdf5 object

        TODO: doesn't support:
          - matrices
          - BucklingEigenvalues

        """
        from pyNastran.op2.op2_interface.hdf5_interface import export_op2_to_hdf5_file
        export_op2_to_hdf5_file(hdf5_filename, self)

    def combine_results(self, combine=True):
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

        # set subcase_key
        for result_type in result_types:
            result = getattr(self, result_type)
            case_keys = sorted(result.keys())
            unique_isubcases = []
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
            for key, value in iteritems(self.subcase_key):
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
            result = getattr(self, result_type)
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
                    isubcase, analysis_code, sort_code, count, isuperelmemnt_adaptivity_index, pval_step = key0
                    key1 = (isubcase, analysis_code, 1, count, isuperelmemnt_adaptivity_index, pval_step)
                    key2 = (isubcase, analysis_code, 2, count, isuperelmemnt_adaptivity_index, pval_step)
                else:
                    isubcase, analysis_code, sort_code, count, isuperelmemnt_adaptivity_index, pval_step, ogs = key0
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

                    isubcase, analysis_code, sort_code, count, isuperelmemnt_adaptivity_index = key0
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
                    else:
                        self.log.info("res=%s has combine" % res1.__class__.__name__)
                    res2 = result[key2]
                    del result[key1]
                    del result[key2]
                    res1.combine(res2)
                    result[isubcase] = res1
                    # print('r[isubcase] =', result[isubcase])
                else:
                    #self.log.info("continue")
                    continue
            setattr(self, result_type, result)
        #print('subcase_key =', self.subcase_key)

        subcase_key2 = {}
        for result_type in result_types:
            if result_type == 'eigenvalues':
                continue
            result = getattr(self, result_type)
            case_keys = list(result.keys())
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
                    assert not isinstance(case_key, string_types), result_type
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

    def get_key_order(self):
        keys = []
        table_types = self.get_table_types()
        for table_type in sorted(table_types):
            result_type_dict = getattr(self, table_type)
            if len(result_type_dict) == 0:
                continue
            for key in result_type_dict:
                if isinstance(key, string_types):
                    if table_type not in ['eigenvalues']:
                        print(table_type)
                    continue
                if key not in keys:
                    keys.append(key)

        #print(self.get_op2_stats())
        keys_order = []

        # subcase_ids = self.subcase_key.keys()

        #self.isubcase_name_map[self.isubcase] = [self.subtitle, self.analysis_code, self.label]
        #subcase_ids = list(self.isubcase_name_map.keys())
        #subcase_ids.sort()
        #print('subcase_ids =', subcase_ids)


        # isubcase, analysis_code, sort_method, count, ogs, superelement_adaptivity_index, pval_step
        #(1, 2, 1, 0, 0, 'SUPERELEMENT 0')  : result1

        isubcases = set([])
        analysis_codes = set([])
        sort_methods = set([])
        counts = set([])
        ogss = set([])
        superelement_adaptivity_indexs = set([])
        pval_steps = set([])

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

        isubcases = list(isubcases)
        analysis_codes = list(analysis_codes)
        sort_methods = list(sort_methods)
        counts = list(counts)
        ogss = list(ogss)
        superelement_adaptivity_indexs = list(superelement_adaptivity_indexs)
        pval_steps = list(pval_steps)

        isubcases.sort()
        analysis_codes.sort()
        sort_methods.sort()
        counts.sort()
        ogss.sort()
        superelement_adaptivity_indexs.sort()
        pval_steps.sort()

        keys3 = []
        for isubcase in isubcases:
            for count in counts:
                for analysis_code in analysis_codes:
                    for superelement_adaptivity_index in superelement_adaptivity_indexs:
                        for pval_step in pval_steps:
                            for sort_method in sort_methods:
                                for ogs in ogss:
                                    key = (isubcase, analysis_code, sort_method,
                                           count, ogs, superelement_adaptivity_index, pval_step)
                                    if key not in keys3:
                                        #print('adding ', key)
                                        keys3.append(key)
        if len(keys3) == 0:
            self.log.warning('No results...\n' + self.get_op2_stats(short=True))
        #assert len(keys3) > 0, keys3
        return keys3

    def print_subcase_key(self):
        self.log.info('---self.subcase_key---')
        for isubcase, keys in sorted(iteritems(self.subcase_key)):
            if len(keys) == 1:
                self.log.info('subcase_id=%s : keys=%s' % (isubcase, keys))
            else:
                self.log.info('subcase_id=%s' % isubcase)
                for key in keys:
                    self.log.info('  %s' % str(key))
        #self.log.info('subcase_key = %s' % self.subcase_key)

    def transform_displacements_to_global(self, i_transform, coords, xyz_cid0=None, debug=False):
        """
        Transforms the ``data`` of displacement-like results into the
        global coordinate system for those nodes with different output
        coordinate systems. Takes indicies and transformation matricies
        for nodes with their output in coordinate systems other than the
        global.

        Used in combination with  ``BDF.get_displacement_index``

        Parameters
        ----------
        i_transform : dict{int cid : int ndarray}
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
        disp_like_dicts = [
            # should NO results be transformed?
            #self.displacements_NO, self.velocities_NO, self.accelerations_NO,
            #self.spc_forces_NO, self.mpc_forces_NO,

            self.displacements,
            self.displacements_ATO, self.displacements_CRM, self.displacements_PSD, self.displacements_RMS,
            self.displacements_scaled,
            self.displacement_scaled_response_spectra_ABS,
            self.displacement_scaled_response_spectra_NRL,

            self.velocities,
            self.velocities_ATO, self.velocities_CRM, self.velocities_PSD, self.velocities_RMS,
            self.velocity_scaled_response_spectra_ABS,

            self.accelerations,
            self.accelerations_ATO, self.accelerations_CRM, self.accelerations_PSD, self.accelerations_RMS,
            self.acceleration_scaled_response_spectra_ABS,
            self.acceleration_scaled_response_spectra_NRL,

            self.eigenvectors,
            self.eigenvectors_RADCONS, self.eigenvectors_RADEFFM,
            self.eigenvectors_RADEATC, self.eigenvectors_ROUGV1,

            self.spc_forces, self.spc_forces_ATO, self.spc_forces_CRM, self.spc_forces_PSD, self.spc_forces_RMS,
            self.mpc_forces, self.mpc_forces_ATO, self.mpc_forces_CRM, self.mpc_forces_PSD, self.mpc_forces_RMS,

            self.applied_loads, self.load_vectors,
        ]
        for disp_like_dict in disp_like_dicts:
            if not disp_like_dict:
                continue
            #print('-----------')
            for subcase, result in iteritems(disp_like_dict):
                #print('result.name = ', result.class_name)
                data = result.data
                for cid, inode in iteritems(i_transform):
                    if cid in [-1, 0]:
                        continue
                    coord = coords[cid]
                    coord_type = coord.type
                    cid_transform = coord.beta()

                    # a global coordinate system has 1s along the main diagonal
                    is_global_cid = False
                    if np.array_equal([1., 1., 1.], np.diagonal(cid_transform)):
                        is_global_cid = True

                    if not is_global_cid and debug:
                        self.log.debug('coord\n%s' % coord)
                        self.log.debug(cid_transform)
                        self.log.debug('inode = %s' % [str(val).rstrip('L') for val in inode.tolist()])
                        self.log.debug('data.shape = %s' % str(data.shape))
                        self.log.debug('len(inode) = %s' % len(inode))
                        assert np.array_equal(inode, np.unique(inode))
                    if coord_type in ['CORD2R', 'CORD1R']:
                        if is_global_cid:
                            #self.log.debug('is_global_cid')
                            continue
                        #self.log.debug('rectangular')


                        # isat_tran.op2
                        #  - nspoint = 4
                        #  - ngrid = 5379
                        #
                        #  - ntotal = 5383
                        #  data.shape = (101, 8, 6)
                        translation = data[:, inode, :3]
                        rotation = data[:, inode, 3:]
                        data[:, inode, :3] = translation.dot(cid_transform)
                        data[:, inode, 3:] = rotation.dot(cid_transform)
                    elif coord_type in ['CORD2C', 'CORD1C']:
                        #self.log.debug('cylindrical')
                        if xyz_cid0 is None:
                            msg = ('xyz_cid is required for cylindrical '
                                   'coordinate transforms')
                            raise RuntimeError(msg)
                        xyzi = xyz_cid0[inode, :]
                        rtz_cid = coord.xyz_to_coord_array(xyzi)
                        theta = rtz_cid[:, 1]
                        #self.log.debug('theta = %s' % list(theta))
                        for itime in range(data.shape[0]):
                            translation = data[itime, inode, :3]
                            rotation = data[itime, inode, 3:]
                            #theta_max1 = translation[:, 1].max()
                            #theta_max2 = rotation[:, 1].max()
                            #print('theta_max = ', max(theta_max1, theta_max2))

                            if 0:
                                # actually somewhat close
                                translation[:, 1] += theta
                                rotation[:, 1] += theta
                                translation = coord.coord_to_xyz_array(data[itime, inode, :3])
                                rotation = coord.coord_to_xyz_array(data[itime, inode, 3:])
                            elif 0:
                                # bad for spc (global_radial_cd)
                                translation[:, 1] += theta
                                rotation[:, 1] += theta
                                translation = coord.coord_to_xyz_array(translation)
                                rotation = coord.coord_to_xyz_array(rotation)
                            elif 0:
                                # very wrong...
                                #print(translation.shape)
                                translation[:, 1] += theta
                                rotation[:, 1] += theta
                            elif 1:
                                pass
                            else:
                                raise RuntimeError('no option selected...')

                            #if is_global_cid:
                                #data[itime, inode, :3] = translation
                                #data[itime, inode, 3:] = rotation
                                #print('is_global_cid')
                                #continue
                            data[itime, inode, :3] = translation.dot(cid_transform)
                            data[itime, inode, 3:] = rotation.dot(cid_transform)
                    elif coord_type in ['CORD2S', 'CORD1S']:
                        #print('spherical')
                        if xyz_cid0 is None:
                            msg = ('xyz_cid is required for spherical '
                                   'coordinate transforms')
                            raise RuntimeError(msg)
                        xyzi = xyz_cid0[inode, :]
                        rtp_cid = coord.xyz_to_coord_array(xyzi)
                        theta = rtp_cid[:, 1]
                        phi = rtp_cid[:, 2]
                        for itime in range(data.shape[0]):
                            translation = data[itime, inode, :3]
                            rotation = data[itime, inode, 3:]
                            if 0:
                                translation[:, 1] += theta
                                translation[:, 2] += phi
                                rotation[:, 1] += theta
                                rotation[:, 2] += phi
                            translation = coord.coord_to_xyz_array(translation)
                            rotation = coord.coord_to_xyz_array(rotation)
                            if is_global_cid:
                                data[itime, inode, :3] = translation
                                data[itime, inode, 3:] = rotation
                                continue
                            data[itime, inode, :3] = translation.dot(cid_transform)
                            data[itime, inode, 3:] = rotation.dot(cid_transform)
                    else:
                        raise RuntimeError(coord)

    def transform_gpforce_to_global(self, nids_all, nids_transform, i_transform, coords, xyz_cid0=None):
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
        i_transform : dict{int cid : int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``BDF.point_ids`` that their output (`CD`) in that
            coordinate system.
        coords : dict{int cid :Coord()}
            Dictionary of coordinate id to the coordinate object
            Use this if CD is only rectangular
            Use this if CD is not rectangular
        xyz_cid0 : ???
            ???

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
            for subcase, result in iteritems(disp_like_dict):
                self.log.debug('result.name = %s' % result.class_name)
                data = result.data

                # inode_xyz :
                #    the indices of the nodes in the model grid point list
                for cid, inode_xyz in iteritems(i_transform):
                    self.log.debug('cid = %s' % cid)
                    if cid in [-1, 0]:
                        continue
                    coord = coords[cid]
                    coord_type = coord.type
                    cid_transform = coord.beta()

                    # a global coordinate system has 1s along the main diagonal
                    is_global_cid = False
                    if np.array_equal([1., 1., 1.], np.diagonal(cid_transform)):
                        is_global_cid = True
                    nids = np.array(nids_transform[cid])

                    #from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
                    #result = RealGridPointForcesArray()
                    if result.is_unique: # TODO: doesn't support preload
                        #self.node_element = zeros((self.ntimes, self.ntotal, 2), dtype='int32')
                        nids_all_gp = result.node_element[0, :, 0]

                        self.log.debug('nids_all_gp = %s' % list(nids_all_gp))
                        self.log.debug('nids = %s' % list(nids))

                        # the indices of the grid points that we're transforming
                        inode_gp = np.where(np.in1d(nids_all_gp, nids))[0]
                        self.log.debug('inode_gp = %s' % list(inode_gp))
                        nids_gp = nids_all_gp[inode_gp]
                        self.log.debug('nids_gp = %s' % list(nids_gp))
                        assert np.array_equal(np.unique(nids_gp), np.unique(nids)), 'nids_gp=%s nids=%s' % (nids_gp, nids)
                        self.log.debug('---------')
                        # the transformation index to go from xyz to the grid point forces
                        inode_gp_xyz = np.searchsorted(nids_all, nids_all_gp)
                        self.log.debug('nids_all = %s' % list(nids_all))
                        self.log.debug('nids_all_gp = %s' % list(nids_all_gp))
                        self.log.debug('inode_xyz = %s' % list(inode_xyz))
                        self.log.debug('inode_gp_xyz = %s' % list(inode_gp_xyz))
                        sys.stdout.flush()
                        assert len(inode_gp_xyz) == len(nids_all_gp), len(nids_all_gp)

                    else:
                        raise NotImplementedError(result)

                    if coord_type in ['CORD2R', 'CORD1R']:
                        if is_global_cid:
                            continue
                        #print('coord\n', coord)
                        #print(cid_transform)
                        #print('inode = %s' % inode)
                        #print('rectangular')
                        translation = data[:, inode_gp, :3]
                        rotation = data[:, inode_gp, 3:]
                        data[:, inode_gp, :3] = translation.dot(cid_transform)
                        data[:, inode_gp, 3:] = rotation.dot(cid_transform)
                    elif coord_type in ['CORD2C', 'CORD1C']:
                        #print('cylindrical')
                        if xyz_cid0 is None:
                            msg = ('xyz_cid is required for cylindrical '
                                   'coordinate transforms')
                            raise RuntimeError(msg)
                        xyzi = xyz_cid0[inode_xyz, :]
                        rtz_cid = coord.xyz_to_coord_array(xyzi)
                        theta_xyz = rtz_cid[inode_xyz, 1]
                        theta = rtz_cid[inode_gp_xyz, 1]
                        self.log.debug('coord\n%s' % coord)
                        self.log.debug(cid_transform)
                        #print('theta_xyz = %s' % list(theta_xyz))
                        #print('theta     = %s' % list(theta))

                        for itime in range(data.shape[0]):
                            #inode = np.where(np.in1d(nids_all, 2))[0]
                            #print('start', data[itime, inode_gp, :3])
                            translation = data[itime, inode_gp, :3]
                            rotation = data[itime, inode_gp, 3:]
                            if 0:
                                # horrible
                                translation = coord.coord_to_xyz_array(data[itime, inode_gp, :3])
                                rotation = coord.coord_to_xyz_array(data[itime, inode_gp, 3:])
                                data[itime, inode_gp, :3] = translation.dot(cid_transform)
                                data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
                            elif 0:
                                # expectedly horrible
                                translation[:, 1] += theta
                                rotation[:, 1] += theta
                                translation = coord.coord_to_xyz_array(data[itime, inode_gp, :3])
                                rotation = coord.coord_to_xyz_array(data[itime, inode_gp, 3:])
                                data[itime, inode_gp, :3] = translation.dot(cid_transform)
                                data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
                            elif 0:
                                # actually not bad...
                                translation = coord.xyz_to_coord_array(translation)
                                rotation = coord.xyz_to_coord_array(rotation)
                                translation[:, 1] += theta
                                rotation[:, 1] += theta
                                translation = coord.coord_to_xyz_array(translation)
                                rotation = coord.coord_to_xyz_array(rotation)
                                data[itime, inode_gp, :3] = translation.dot(cid_transform)
                                data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
                            elif 0:
                                # not that bad, worse than previous
                                translation = coord.xyz_to_coord_array(translation)
                                rotation = coord.xyz_to_coord_array(rotation)
                                translation[:, 1] += theta
                                rotation[:, 1] += theta
                                translation = coord.coord_to_xyz_array(translation)
                                rotation = coord.coord_to_xyz_array(rotation)
                                data[itime, inode_gp, :3] = translation.dot(cid_transform.T)
                                data[itime, inode_gp, 3:] = rotation.dot(cid_transform.T)
                            elif 1:
                                # doesn't work...actually pretty close
                                data[itime, inode_gp, :3] = translation.dot(cid_transform)
                                data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
                            elif 0:
                                # very, very close
                                data[itime, inode_gp, :3] = translation.dot(cid_transform.T)
                                data[itime, inode_gp, 3:] = rotation.dot(cid_transform.T)
                            elif 0:
                                # is this just the same as one of the previous?
                                data[itime, inode_gp, :3] = cid_transform.T.dot(translation.T).T
                                data[itime, inode_gp, 3:] = cid_transform.T.dot(rotation.T).T
                            else:
                                raise RuntimeError('no option selected...')

                            #if is_global_cid:
                                #data[itime, inode, :3] = translation
                                #data[itime, inode, 3:] = rotation
                                #continue
                            #data[itime, inode_gp, :3] = translation.dot(cid_transform)
                            #data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
                            #print('end', data[itime, inode_gp, :3])

                    elif coord_type in ['CORD2S', 'CORD1S']:
                        self.log.debug('spherical')
                        if xyz_cid0 is None:
                            msg = ('xyz_cid is required for spherical '
                                   'coordinate transforms')
                            raise RuntimeError(msg)
                        xyzi = xyz_cid0[inode, :]
                        rtp_cid = coord.xyz_to_coord_array(xyzi)
                        theta = rtp_cid[:, 1]
                        phi = rtp_cid[:, 2]
                        for itime in range(data.shape[0]):
                            translation = data[itime, inode, :3]
                            rotation = data[itime, inode, 3:]
                            #if 0:
                                #translation[:, 1] += theta
                                #translation[:, 2] += phi
                                #rotation[:, 1] += theta
                                #rotation[:, 2] += phi
                            translation = coord.coord_to_xyz_array(data[itime, inode, :3])
                            rotation = coord.coord_to_xyz_array(data[itime, inode, 3:])
                            #if is_global_cid:
                                #data[itime, inode, :3] = translation
                                #data[itime, inode, 3:] = rotation
                                #continue
                            data[itime, inode, :3] = translation.dot(cid_transform)
                            data[itime, inode, 3:] = rotation.dot(cid_transform)
                    else:
                        raise RuntimeError(coord)
        self.log.debug('-----------')
        return

def main():  # pragma: no cover
    """testing new ideas"""
    import pyNastran
    pkg_path = pyNastran.__path__[0]

    # we don't want the variable name to get picked up by the class
    _op2_filename = os.path.join(pkg_path, '..', 'models',
                                 'sol_101_elements', 'solid_shell_bar.op2')

    model = OP2()
    model.read_op2(_op2_filename)
    isubcase = 1

    # ============displacement================
    # same for velocity/acceleration/mpc forces/spc forces/applied load
    # maybe not temperature because it's only 1 result...
    displacement = model.displacements[isubcase]
    data = displacement.data
    grid_type = model.displacements[isubcase].grid_type

    # t1, t2, t3, r1, r2, r3
    assert displacement.ntimes == 1, displacement.ntimes

    # get all the nodes for element 1
    inode1 = data.getNodeIndex([1])    # [itransient, node, t1/t2]
    datai = data[0, inode1, :]
    grid_typei = grid_type[inode1]

    # ============solid stress=================
    # same for solid strain
    solid_stress = model.ctetra_stress[isubcase]
    data = solid_stress.data
    cid = model.ctetra_stress[isubcase].cid

    # oxx, oyy, ozz, txy, tyz, txz
    assert solid_stress.ntimes == 1, solid_stress.ntimes

    # get the indexs for cid, element 1
    ielem1 = solid_stress.getElementPropertiesIndex([1])  # element-specific properties
    datai = cid[ielem1]

    # get all the nodes for element 1
    ielem1 = solid_stress.getElementIndex([1])
    # [itransient, elem*node, oxx/oyy, etc.]
    datai = data[0, ielem1, :]

    # get all the nodes for element 1
    ielem1 = solid_stress.getElementIndex([1])
    # [itransient, elem*node, oxx/oyy, etc.]
    datai = data[0, ielem1, :]

    # get all the nodes for element 4 and 5
    ielem45 = solid_stress.getElementIndex([[1, 4, 5]])
    datai = data[0, ielem45, :]

    # get the index for element 1, centroid
    ielem1_centroid = solid_stress.getElementNodeIndex([[1, 0]])
    datai = data[0, ielem1_centroid, :]

    # get the index for element 1, node 1
    ielem1_node1 = solid_stress.getElementNodeIndex([[1, 1]])
    datai = data[0, ielem1_node1, :]

    # get the index for element 1, node 1 and element 1, node 2
    ielem1_node12 = solid_stress.getElementNodeIndex([[1, 1],
                                                      [1, 2],])
    datai = data[0, ielem1_node12, :]

    # ============plate stress=================
    # same for plate strain
    plate_stress = model.cquad4_stress[isubcase]
    data = plate_stress.data

    # oxx, oyy, ozz, txy, tyz, txz
    assert plate_stress.ntimes == 1, plate_stress.ntimes

    # get all the nodes for element 1
    ielem1 = plate_stress.getElementIndex([1])
    # [itransient, elem*node, oxx/oyy, etc.]
    datai = plate_stress[0, ielem1, :]

    # ========composite plate stress============
    # same for plate strain
    comp_plate_stress = model.cquad4_composite_stress[isubcase]
    data = comp_plate_stress.data
    cid = model.cquad4_stress[isubcase].cid

    # get the indexs for cid, element 1
    ielem1 = comp_plate_stress.getElementPropertiesIndex([1])  # element-specific properties
    datai = cid[ielem1]

    # oxx, oyy, ozz, txy, tyz, txz
    assert comp_plate_stress.ntimes == 1, comp_plate_stress.ntimes

    # get all the nodes/layers for element 1
    ielem1 = comp_plate_stress.getElementIndex([1])
    # [itransient, elem*node*layer, oxx/oyy, etc.]
    datai = data[0, ielem1, :]

    # get all the layers for element 1, centroid, and all the layers
    ielem1_centroid = solid_stress.getElementNodeIndex([[1, 0]])
    datai = data[0, ielem1_centroid, :]

    # get the index for element 1, centroid, layer 0
    ielem1_centroid_layer = solid_stress.getElementNodeLayerIndex([[1, 0, 0]])
    datai = data[0, ielem1_centroid_layer, :]

    # get the index for element 1, layer 0, and all the nodes
    ielem1_layer = solid_stress.getElementLayerIndex([[1, 0]])
    datai = data[0, ielem1_layer, :]


if __name__ == '__main__':  # pragma: no cover
    main()
