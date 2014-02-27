#pylint: disable=C0103,W0201,W0223,R0901,R0902,R0904
"""
Main OP2 class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
#import warnings
from numpy import array, unique, where

#from pyNastran.utils import is_binary
#from pyNastran.utils.gui_io import load_file_dialog
from pyNastran.op2.dev_explicit.op2 import OP2

class OP2_Vectorized(OP2):

    def __init__(self, make_geom=False, save_skipped_cards=False,
                 debug=True, log=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param save_skipped_cards: creates the skippedCards.out file (default=False)
        :param debug: prints data about how the OP2 was parsed (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        """
        debug = False
        OP2.__init__(self, make_geom=make_geom, save_skipped_cards=save_skipped_cards,
                     debug=debug, log=log)
        #print(self.binary_debug)
        #self.binary_debug = None
        self.ask = False

    def set_as_vectorized(self, ask=False):
        """
        Enables vectorization (currently subject to change and break frequently)

        The code will degenerate to dictionary based results when
        a result does not support vectorization.

        Vectorization is always True here.
        :param ask:  Do you want to see a GUI of result types.

        +--------+---------------+---------+------------+
        | Case # | Vectorization |  Ask    | Read Modes |
        +========+===============+=========+============+
        |    1   | True          |  True   |  1, 2      |
        |    2   | True          |  False  |  1, 2      |
        |    3   | False         |  True   |  1, 2      |
        |    4   | False         |  False  |  0         | <------ < v0.8
        +--------+---------------+---------+------------|

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

    def read_op2(self, op2_filename=None):
        """
        Starts the OP2 file reading
        """
        assert self.ask in [True, False], self.ask
        self.is_vectorized = True
        if self.is_vectorized:
            self.log.info('-------- reading the op2 with read_mode=1 --------')
            self.read_mode = 1
            self._close_op2 = False

            # get GUI object names, build objects, but don't read data
            OP2.read_op2(self, op2_filename=op2_filename)

            # TODO: stuff to figure out objects
            # TODO: stuff to show gui of table names
            # TODO: clear out objects the user doesn't want
            self.read_mode = 2
            self._close_op2 = True
            self.log.info('-------- reading the op2 with read_mode=2 --------')
            OP2.read_op2(self, op2_filename=op2_filename)
        else:
            #self.read_mode = 0
            OP2.read_op2(self, op2_filename=op2_filename)
            return
            #raise NotImplementedError()
        self.f.close()
        self.skippedCardsFile.close()
        self.combine_results()
        self.log.info('finished reading op2')

    def combine_results(self, combine=True):
        """
        we want the data to be in the same format and grouped by subcase, so
        we take

        stress = {
            (1, 'SUPERELEMENT 0') : result1,
            (1, 'SUPERELEMENT 10') : result2,
            (1, 'SUPERELEMENT 20') : result3,
            (2, 'SUPERELEMENT 0') : result4,
        }
        and convert it to:

        stress = {
            1 : result1 + result2 + results3,
            2 : result4,
        }
        """
        if not combine:
            return
        self.log.info('compress_results')
        result_types = [
            'displacements', 'eigenvectors', 'spcForces', 'mpcForces',

            'ctetra_stress', 'chexa_stress', 'cpenta_stress',
            'ctetra_strain', 'chexa_strain', 'cpenta_strain',
            'crod_stress', 'conrod_stress', 'ctube_stress',
            'crod_strain', 'conrod_strain', 'ctube_strain',

            'plateStress', 'plateStrain',
            'barStress', 'beamStress',
            'barStrain', 'beamStrain',

            'rodForces',
            'barForces',
            'beamForces',
            'plateForces', 'plateForces2',
            'strainEnergy',
            ]

        for result_type in result_types:
            result = getattr(self, result_type)
            case_keys = sorted(result.keys())
            isubcases = [isubcase for (isubcase, name) in case_keys]
            unique_isubcases = unique(isubcases)

            #print(case_keys)
            #print(isubcases)
            #print(unique_isubcases)
            for isubcase in unique_isubcases:
                #wherei = where(isubcases==isubcase)[0]
                wherei = [i for i, isubcase in enumerate(isubcases) if isubcase == 1]
                keys = [case_keys[i] for i in wherei]
                if len(wherei) == 1:
                    # rename the case since we have only one tuple for the result
                    result[isubcase] = result[keys[0]]
                    del result[keys[0]]
                else:
                    # multiple results to combine
                    res1 = result[keys[0]]
                    if not hasattr(res1, 'combine'):
                        print("res1=%s has no method combine" % res1.__class__.__name__)
                        continue

                    #raise NotImplementedError('multiple results to combine')
                    del result[keys[0]]

                    results_list = []
                    for key in keys[1:]:
                        results_list.append(result[key])
                        del result[key]

                    #res1.data = hstack(res1.data, res2.data)
                    # combination method depends on result, so stresses are
                    # combined differently than displacements
                    res1.combine(results_list)
                    result[isubcase] = res1

                #print(keys)
            setattr(self, result_type, result)


if __name__ == '__main__':
    import pyNastran
    pkg_path = pyNastran.__path__[0]

    # we don't want the variable name to get picked up by the class
    _op2_filename = os.path.join(pkg_path, '..', 'models', 'sol_101_elements', 'solid_shell_bar.op2')

    model = OP2_Vectorized()
    model.set_as_vectorized(ask=False)
    model.read_op2(_op2_filename)
    isubcase = 1

    # ============displacement================
    # same for velocity/acceleration/mpc forces/spc forces/applied load
    # maybe not temperature because it's only 1 result...
    displacement = model.displacement[isubcase]
    data = displacement.data
    grid_type = model.displacement[isubcase].grid_type

    # t1, t2, t3, r1, r2, r3
    assert displacement.ntimes == 1, displacement.ntimes

    # get all the nodes for element 1
    inode1 = data.getNodeIndex( [1] )
    # [itransient, node, t1/t2]
    datai = data[0, inode1, :]
    grid_typei = grid_type[inode]

    # ============solid stress=================
    # same for solid strain
    solid_stress = model.solidStress[isubcase]
    data = solid_stress.data
    cid = model.solidStress[isubcase].cid

    # oxx, oyy, ozz, txy, tyz, txz
    assert solid_stress.ntimes == 1, solid_stress.ntimes

    # get the indexs for cid, element 1
    ielem1 = solid_stress.getElementPropertiesIndex( [1] )  # element-specific properties
    datai = cid[ielem1]

    # get all the nodes for element 1
    ielem1 = solid_stress.getElementIndex( [1] )
    # [itransient, elem*node, oxx/oyy, etc.]
    datai = data[0, ielem1, :]

    # get all the nodes for element 1
    ielem1 = solid_stress.getElementIndex( [1] )
    # [itransient, elem*node, oxx/oyy, etc.]
    datai = data[0, ielem1, :]

    # get all the nodes for element 4 and 5
    ielem45 = solid_stress.getElementIndex( [[1, 4, 5]] )
    datai = data[0, ielem45, :]

    # get the index for element 1, centroid
    ielem1_centroid = solid_stress.getElementNodeIndex( [[1, 0]] )
    datai = data[0, ielem1_centroid, :]

    # get the index for element 1, node 1
    ielem1_node1 = solid_stress.getElementNodeIndex( [[1, 1]] )
    datai = data[0, ielem1_node1, :]

    # get the index for element 1, node 1 and element 1, node 2
    ielem1_node12 = solid_stress.getElementNodeIndex( [[1, 1],
                                                       [1, 2],])
    datai = data[0, ielem1_node12, :]

    # ============plate stress=================
    # same for plate strain
    plate_stress = model.plateStress[isubcase]
    data = plate_stress.data

    # oxx, oyy, ozz, txy, tyz, txz
    assert plate_stress.ntimes == 1, plate_stress.ntimes

    # get all the nodes for element 1
    ielem1 = plate_stress.getElementIndex( [1] )
    # [itransient, elem*node, oxx/oyy, etc.]
    datai = plate_stress[0, ielem1, :]

    # ========composite plate stress============
    # same for plate strain
    comp_plate_stress = model.compositePlateStress[isubcase]
    data = comp_plate_stress.data
    cid = model.plateStress[isubcase].cid

    # get the indexs for cid, element 1
    ielem1 = comp_plate_stress.getElementPropertiesIndex( [1] )  # element-specific properties
    datai = cid[ielem1]

    # oxx, oyy, ozz, txy, tyz, txz
    assert comp_plate_stress.ntimes == 1, comp_plate_stress.ntimes

    # get all the nodes/layers for element 1
    ielem1 = comp_plate_stress.getElementIndex( [1] )
    # [itransient, elem*node*layer, oxx/oyy, etc.]
    datai = data[0, ielem1, :]

    # get all the layers for element 1, centroid, and all the layers
    ielem1_centroid = solid_stress.getElementNodeIndex( [[1, 0]] )
    datai = data[0, ielem1_centroid, :]

    # get the index for element 1, centroid, layer 0
    ielem1_centroid_layer = solid_stress.getElementNodeLayerIndex( [[1, 0, 0]] )
    datai = data[0, ielem1_centroid_layer, :]

    # get the index for element 1, layer 0, and all the nodes
    ielem1_layer = solid_stress.getElementLayerIndex( [[1, 0]] )
    datai = data[0, ielem1_layer, :]
