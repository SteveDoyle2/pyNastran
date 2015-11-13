from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import unittest
#from numpy import allclose, array
#from numpy.linalg import norm

import pyNastran
from pyNastran.bdf.bdf import BDF

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'bdf', 'test')
#print("test_path = %r" % test_path)


class BDFUpdater(BDF):
    def __init__(self):
        BDF.__init__(self)
        self._type_map = {
            'GRID': self.nodes,
            #'RINGAX': self.nodes,

            'MAT1': self.materials,
            'MAT2': self.materials,
            'MAT3': self.materials,
            'MAT4': self.materials,
            'MAT5': self.materials,
            'MAT8': self.materials,
            'MAT10': self.materials,

            # shell
            'CTRIAX': self.elements,
            'CTRIA3': self.elements,
            'CTRIA6': self.elements,
            'CQUAD': self.elements,
            'CQUAD4': self.elements,
            'CQUAD8': self.elements,
            'CQUADX': self.elements,
            'PCOMP': self.properties,
            'PCOMPG': self.properties,
            'PSHELL': self.properties,

            # shear,
            'CSHEAR': self.elements,
            'PSHEAR': self.properties,

            # plane
            #'PLPLANE' : self.elements,

            # solid
            'CTETRA' : self.elements,
            'CPENTA' : self.elements,
            'CHEXA'  : self.elements,
            'PSOLID' : self.properties,
            'PLSOLID': self.properties,

            # rod
            'CONROD': self.elements,
            'CROD'  : self.elements,
            'PROD'  : self.properties,

            # tube
            'CTUBE' : self.elements,
            'PTUBE' : self.properties,

            # beam
            'CBEAM'  : self.elements,
            'PBEAM'  : self.properties,
            'PBEAML' : self.properties,

            # bar
            'CBAR'  : self.elements,
            'PBAR'  : self.properties,
            'PBARL' : self.properties,

            # bend
            'CBEND': self.elements,
            'PBEND': self.properties,
            # bush
            'CBUSH'  : self.elements,
            'CBUSH1D': self.elements,
            'CBUSH2D': self.elements,
            'PBUSH'  : self.properties,
            'PBUSH1D': self.properties,
            #'PBUSH2D': self.properties,
            #'PBUSHT' : self.properties,

            # spring
            'CELAS1': self.elements,
            'CELAS2': self.elements,
            'CELAS3': self.elements,
            'CELAS4': self.elements,
            'PELAS' : self.properties,

            # dampers
            'CFAST': self.elements,
            'PFAST': self.properties,

            # dampers
            'CDAMP1': self.elements,
            'CDAMP2': self.elements,
            'CDAMP3': self.elements,
            'CDAMP4': self.elements,
            'CDAMP5': self.elements,
            'PDAMP' : self.properties,
            'PDAMPT': self.properties,
            'PDAMP5': self.properties,
            'PVISC' : self.properties,

            #'CRAC2D' : self.elements,
            #'CRAC3D' : self.elements,
            #'PRAC2D' : self.properties,
            #'PRAC3D' : self.properties,

            # mass
            #'CONM1': self.elements,
            'CONM2' : self.elements,
            'CMASS1': self.elements,
            'CMASS2': self.elements,
            'CMASS3': self.elements,
            'CMASS4': self.elements,
            'PMASS' : self.properties,
            #'NSM'   : self.elements, # ???

            # rigid elements
            'RBAR'  : self.rigid_elements,
            'RBAR1' : self.rigid_elements,
            'RBE1'  : self.rigid_elements,
            'RBE2'  : self.rigid_elements,
            'RBE3'  : self.rigid_elements,

            # methods
            'EIGB' : self.methods,
            'EIGC' : self.methods,
            'EIGP' : self.methods,
            'EIGR' : self.methods,
            'EIGRL': self.methods,

            # frequencies
            #'FREQ' : self.freqs,
            #'FREQ1' : self.freqs,
            #'FREQ2' : self.freqs,

            # param
            'PARAM' : self.params,

            # ========= remove these =============
            # loads
            'FORCE'   : self.loads,
            'FORCE1'  : self.loads,
            'FORCE2'  : self.loads,
            'MOMENT'  : self.loads,
            'MOMENT1' : self.loads,
            'MOMENT2' : self.loads,
            'RFORCE'  : self.loads,
            'LOAD'    : self.loads,
            'DLOAD'   : self.loads,
            'SLOAD'   : self.loads,
            'TLOAD1'  : self.loads,
            'TLOAD2'  : self.loads,
            'RLOAD1'  : self.loads,
            'RLOAD2'  : self.loads,
            'DAREA'   : self.loads,
            'GRAV'    : self.loads,
            #'ACCEL'  : self.loads,
            'ACCEL1'  : self.loads,

            # constraints
            'SPC1'   : self.spcObject,
            'SPC'    : self.spcObject,
            'SPCADD' : self.spcObject,
            'SPCAX'  : self.spcObject,

            'MPC'    : self.mpcObject,
            'MPCADD' : self.mpcObject,

        }
        for key in self._type_map:
            assert key in self.cards_to_read, 'key=%r is not a valid card' % key

    def update_card(self, Type, iType, ifield, value):
        for key in self.card_count:
            if key not in self._type_map:
                raise RuntimeError('add %r to self._type_map' % key)

        # get the storage object
        try:
            objs = self._type_map[Type]
        except KeyError:
            raise KeyError('Updating card Type=%r is not supported' % Type)

        # get the specific card
        try:
            card = objs[iType]
        except KeyError:
            msg = 'Could not find %s ID=%r' % (Type, iType)
            raise KeyError(msg)

        # update the card
        card.update_field(ifield, value)
        return card

class TestOpenMDAO(unittest.TestCase):

    def test_openmdao_good_1(self):
        updates = [
            #['MAT1', 3, 10.0],  # 3 is E -> set to 10.0
            #['MAT1', 4, 10.0],  # 3 is G -> set to 10.0
            ['GRID', 1, 3, 10.0],  # 3 is x1 -> set to 10.0
            ['GRID', 1, 4, 20.0],  # 4 is x2 -> set to 20.0
            ['CPENTA', 9, 2, 10],  # 2 is property_id -> set to 10
            ['CPENTA', 9, 3, 20],  # 3 is node1 -> set to 20
            ['PSOLID', 4, 1, 2],   # 1 is material_id
            ['PARAM', 'WTMASS', 1, 'WTMASs'],  # key
            ['PARAM', 'WTMASS', 2, 0.0025],  # value1
            ['PCOMP', 1, 2, 1.],
            ['PCOMP', 1, 3, 2.],
            ['PCOMP', 1, 12, 'YES_A'],
            ['PCOMP', 1, 16, 'YES_B'],
            ['PCOMP', 1, 20, 'YES_C'],
            ['CTETRA', 8, 3, 1], # nid[0]
            ['CTETRA', 8, 4, 2], # nid[1]
            ['CTETRA', 8, 5, 3], # nid[2]
            ['CTETRA', 8, 6, 4], # nid[3]
        ]
        #GRID           1       0      0.      0.      0.       0
        #GRID           2       0      1.      0.      0.       0
        #GRID           3       0      1.      1.      0.       0
        #GRID           4       0      0.      1.      0.       0
        #CPENTA         9       4      21      22      23      24      25      26
        #PSOLID   4       1       0
        #CTETRA         8       4      11      12      13      15

        bdf_filename = os.path.join(test_path, 'unit', 'test_mass.dat')

        model = BDFUpdater()
        model.read_bdf(bdf_filename)

        for iupdate in updates:
            Type, iType, ifield, value = iupdate
            card = model.update_card(Type, iType, ifield, value)
            #print(card)

    def test_openmdao_bad_1(self):
        updates = [  # KeyError
            ['GRID', 1, 0, 20.0],
            ['GRID', 1, 10, 20.0],
            ['CHEXA', 100, 3, 20],
            ['FAKECARD', 100, 3, 20],
        ]
        bdf_filename = os.path.join(test_path, 'unit', 'test_mass.dat')

        model = BDFUpdater()
        model.read_bdf(bdf_filename)

        for iupdate in updates:
            Type, iType, ifield, value = iupdate
            self.assertRaises(KeyError, model.update_card, Type, iType, ifield, value)
            #print("tried to apply %r=%s" % (Type, iType))

    def test_openmdao_bad_2(self):
        updates = [  # IndexError
            ['PARAM', 'WTMASS', 3, 0.005],  # value2; invalid b/c WTMASS
            ['PCOMP', 1, 24, 'YES_D'],  # too many plies
        ]
        bdf_filename = os.path.join(test_path, 'unit', 'test_mass.dat')

        model = BDFUpdater()
        model.read_bdf(bdf_filename)

        for iupdate in updates:
            Type, iType, ifield, value = iupdate
            self.assertRaises(IndexError, model.update_card, Type, iType, ifield, value)
            #print("tried to apply %r=%s" % (Type, iType))

    def test_openmaod_bad_3(self):
        params_bad = {1 : '10'}
        params_good = {'cat' : '10'}
        model = BDF(debug=False)

        with self.assertRaises(TypeError):
            model.set_dynamic_syntax(params_bad)

        model.set_dynamic_syntax(params_good)
        val = model._parse_dynamic_syntax('cat')
        self.assertEqual('10', val)
        #with self.assertRaises(KeyError):
            #self._parse_dynamic_syntax(key)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
