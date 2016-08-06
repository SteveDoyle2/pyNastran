from __future__ import print_function
import os
import unittest
from six import iteritems
from pyNastran.bdf.bdf import read_bdf

from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris

class TestMeshUtils(unittest.TestCase):

    def test_fix_bad_quads(self):
        msg = [
            'SOL 101',
            'CEND',
            'BEGIN BULK',
            'GRID,1,,0.,0.,0.',
            'GRID,2,,1.,0.,0.',
            'GRID,3,,1.,1.,0.',
            'GRID,4,,0.,1.,0.',
            'GRID,5,,1.,1.,0.00001',
            'GRID,6,,0.,0.,0.00001',
            'CQUAD4,1, 2, 1,2,3,4',
            'CQUAD4,2, 2, 1,2,3,5',
            'CQUAD4,3, 2, 1,6,3,5',
        ]
        bdf_filename = 'fix_bad_quads.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('\n'.join(msg))
        model = read_bdf(bdf_filename, xref=False, debug=False)
        model.cross_reference(xref=True, xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        convert_bad_quads_to_tris(model, tol=0.01)
        for eid, elem in sorted(iteritems(model.elements)):
            print(elem)
        os.remove(bdf_filename)
