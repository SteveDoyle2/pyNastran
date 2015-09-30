#pylint: disable=W0201,W0223,R0901,R0902,R0904
"""
Main OP2 class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os

from numpy import unique, vstack
from pyNastran.op2.op2_scalar import OP2_Scalar

from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import SortCodeError, DeviceCodeError, FortranMarkerError


#from pyNastran.op2.op2_writer import OP2Writer

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
        """
        self.set_mode(mode)
        make_geom = False
        assert make_geom == False, make_geom
        OP2_Scalar.__init__(self,
                     debug=debug, log=log, debug_file=debug_file)
        self.ask = False

    def set_mode(self, mode):
        if mode.lower() == 'msc':
            self.set_as_msc()
        elif mode.lower() == 'nx':
            self.set_as_nx()
        else:
            raise RuntimeError("mode=%r and must be 'msc' or 'nx'")

    def set_as_vectorized(self, ask=False):
        """
        Enables vectorization

        The code will degenerate to dictionary based results when
        a result does not support vectorization.

        Vectorization is always True here.
        :param ask:  Do you want to see a GUI of result types.

        +--------+---------------+---------+------------+
        | Case # | Vectorization |  Ask    | Read Modes |
        +========+===============+=========+============+
        |    1   | True          |  True   |  1, 2      |
        +--------+---------------+---------+------------+
        |    2   | True          |  False  |  1, 2      |
        +--------+---------------+---------+------------+
        |    3   | False         |  True   |  1, 2      |
        +--------+---------------+---------+------------+
        |    4   | False         |  False  |  0         |
        +--------+---------------+---------+------------+

        Definitions
        ===========
          Vectorization - A storage structure that allows for faster read/access
                          speeds and better memory usage, but comes with a more
                          difficult to use data structure.

                          It limits the node IDs to all be integers (e.g. element
                          centroid).  Composite plate elements (even for just CTRIA3s)
                          with an inconsistent number of layers will have a more
                          difficult data structure.  Models with solid elements of
                          mixed type will also be more complicated (or potentially
                          split up).
          Scanning   - a quick check used to figure out how many results to process
                      that takes almost no time
          Reading    - process the op2 data
          Build      - call the __init__ on a results object (e.g. DisplacementObject)
          Start Over - Go to the start of the op2 file
          Ask        - launch a GUI dialog to let the user click which results to load

        Read Mode Definitions
        =====================
          0.   The default OP2 dictionary based-approach with no asking GUI
          1.   The first read of a result to get the shape of the data
               (and result types if vectorization=True or result types
               if vectorization=False)
          2.   The second read of a result to get the results

        Cases
        ======
          1.   Scan the block to get the size, build the object (read_mode=1),
               ask the user, start over, fill the objects (read_mode=2).
               Degenerate to read_mode=0 when read_mode=2 cannot be used based
               upon the value of ask.
          2.   Same as case #1, but don't ask the user.
               Scan the block to get the size, build the object (read_mode=1),
               start over, fill the objects (read_mode=2).
          3.   Scan the block to get the object types (read_mode=1), ask the user,
               build the object & fill it (read_mode=2)
          4.   Read the block to get the size, build the object & fill it (read_mode=0)
        """
        self.ask = ask

    def read_op2(self, op2_filename=None, vectorized=True, combine=True):
        """
        Starts the OP2 file reading

        Parameters
        ----------
        op2_filename : str (default=None -> popup)
            the op2_filename
        vectorized : bool; default=True
            should the vectorized objects be used
        combine : bool; default=True
            True : objects are isubcase based
            False : objects are (isubcase, subtitle) based;
                    will be used for superelements regardless of the option
        """
        self.log.debug('vectorized=%s combine=%s' % (vectorized, combine))
        assert self.ask in [True, False], self.ask
        self.is_vectorized = vectorized
        if self.is_vectorized:
            self.log.debug('-------- reading op2 with read_mode=1 --------')
            self.read_mode = 1
            self._close_op2 = False

            # get GUI object names, build objects, but don't read data
            OP2_Scalar.read_op2(self, op2_filename=op2_filename)

            # TODO: stuff to figure out objects
            # TODO: stuff to show gui of table names
            # TODO: clear out objects the user doesn't want
            self.read_mode = 2
            self._close_op2 = True
            self.log.debug('-------- reading op2 with read_mode=2 --------')
            OP2_Scalar.read_op2(self, op2_filename=self.op2_filename)
        else:
            self.read_mode = 0
            self._close_op2 = True
            OP2_Scalar.read_op2(self, op2_filename=op2_filename)

        self.combine_results(combine=combine)
        self.log.info('finished reading op2')

    def combine_results(self, combine=True):
        """
        we want the data to be in the same format and grouped by subcase, so
        we take

        .. code-block:: python

          stress = {
              # isubcase, analysis_code, sort_method, count, subtitle
              (1, 2, 1, 0, 'SUPERELEMENT 0') : result1,
              (1, 2, 1, 0, 'SUPERELEMENT 10') : result2,
              (1, 2, 1, 0, 'SUPERELEMENT 20') : result3,
              (2, 2, 1, 0, 'SUPERELEMENT 0') : result4,
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

        for result_type in result_types:
            result = getattr(self, result_type)
            case_keys = sorted(result.keys())

            unique_isubcases = []
            for case_key in case_keys:
                if isinstance(case_key, tuple):
                    isubcasei, analysis_codei, sort_methodi, counti, isubtitle = case_key
                    value = (analysis_codei, sort_methodi, counti, isubtitle)
                    if value not in self.subcase_key[isubcasei]:
                        self.subcase_key[isubcasei].append(value)
                else:
                    break
        if not combine:
            return
        del result, case_keys
        isubcases = unique(self.subcase_key.keys())
        unique_isubcases = unique(isubcases)

        self.log.info('combine_results')
        for result_type in result_types:
            result = getattr(self, result_type)
            if len(result) == 0:
                continue
            for isubcase in unique_isubcases:
                keys = self.subcase_key[isubcase]
                key0 = tuple([isubcase] + list(keys[0]))
                if len(keys) == 1:
                    if key0 not in result:
                        continue
                    # rename the case since we have only one tuple for the result
                    # key0 = tuple([isubcase] + list(key0))
                    result[isubcase] = result[key0]
                    del result[key0]
                elif len(keys) == 2:
                    # continue
                    print('key0 =', result_type, key0)
                    # res0 = result[key0]

                    isubcase, analysis_code, sort_code, count, subtitle = key0
                    key1 = (isubcase, analysis_code, 1, count, subtitle)
                    key2 = (isubcase, analysis_code, 2, count, subtitle)
                    if not (key1 in result and key2 in result):
                        if key1 in result:
                            res1 = result[key1]
                            self.log.info("res=%s has a single case; trivial" % res1.__class__.__name__)
                            result[isubcase] = result[key1]
                            del result[key1]
                        elif key2 in result:
                            res2 = result[key2]
                            self.log.info("res=%s has a single case; trivial" % res2.__class__.__name__)
                            result[isubcase] = result[key2]
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
                    self.log.info("continue")
                    continue
            setattr(self, result_type, result)

        subcase_key2 = {}
        for result_type in result_types:
            result = getattr(self, result_type)
            case_keys = sorted(result.keys())
            if len(result) == 0:
                continue
            for isubcase in unique_isubcases:
                for case_key in case_keys:
                    if isubcase in subcase_key2 and case_key not in subcase_key2[isubcase]:
                        subcase_key2[isubcase].append(case_key)
        self.subcase_key = subcase_key2


def main():
    import pyNastran
    pkg_path = pyNastran.__path__[0]

    # we don't want the variable name to get picked up by the class
    _op2_filename = os.path.join(pkg_path, '..', 'models', 'sol_101_elements', 'solid_shell_bar.op2')

    model = OP2()
    model.set_as_vectorized(ask=False)
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
