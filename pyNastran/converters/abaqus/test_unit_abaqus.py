import os
import unittest
from io import StringIO
from cpylog import get_logger

import pyNastran
from pyNastran.converters.abaqus.abaqus import read_abaqus
from pyNastran.converters.abaqus.abaqus_to_nastran import (
    nastran_to_abaqus_filename, abaqus_to_nastran_filename)
from pyNastran.converters.format_converter import cmd_line_format_converter

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, 'converters', 'abaqus', 'models')
NASTRAN_MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestAbaqus(unittest.TestCase):
    def test_abaqus_to_nastran_1(self):
        """plate conversion"""
        nastran_filename = os.path.join(NASTRAN_MODEL_PATH, 'plate', 'plate.bdf')
        abaqus_inp_filename = os.path.join(MODEL_PATH, 'plate_out.inp')
        nastran_to_abaqus_filename(nastran_filename, abaqus_inp_filename)

        #model = read_abaqus(abaqus_filename, debug=True)
        nastran_filename_out = os.path.join(MODEL_PATH, 'plate2_out.bdf')
        abaqus_to_nastran_filename(abaqus_inp_filename, nastran_filename_out)

    def test_abaqus_to_nastran_2(self):
        """plate conversion"""
        #nastran_filename = os.path.join(MODEL_PATH, 'plate.inp')
        abaqus_inp_filename = os.path.join(MODEL_PATH, 'in.inp')
        #nastran_to_abaqus_filename(nastran_filename, abaqus_inp_filename)

        #model = read_abaqus(abaqus_filename, debug=True)
        nastran_filename_out = os.path.join(MODEL_PATH, 'out.bdf')
        abaqus_to_nastran_filename(abaqus_inp_filename, nastran_filename_out)

    def test_abaqus_to_nastran_3(self):
        """ctetra4 conversion"""
        #nastran_filename = os.path.join(MODEL_PATH, 'plate.inp')
        abaqus_inp_filename = os.path.join(MODEL_PATH, 'test_bracket.inp')
        #nastran_to_abaqus_filename(nastran_filename, abaqus_inp_filename)

        #model = read_abaqus(abaqus_filename, debug=True)
        nastran_filename_out = os.path.join(MODEL_PATH, 'test_bracket_out.bdf')
        abaqus_to_nastran_filename(abaqus_inp_filename, nastran_filename_out)

    def test_abaqus_to_nastran_4(self):
        """ctetra4 conversion"""
        #nastran_filename = os.path.join(MODEL_PATH, 'plate.inp')
        abaqus_inp_filename = os.path.join(MODEL_PATH, 'test_bracket_separate.inp')
        #nastran_to_abaqus_filename(nastran_filename, abaqus_inp_filename)

        #model = read_abaqus(abaqus_filename, debug=True)
        nastran_filename_out = os.path.join(MODEL_PATH, 'test_bracket_separate_out.bdf')
        abaqus_to_nastran_filename(abaqus_inp_filename, nastran_filename_out)

    def test_abaqus_1(self):
        """simple test"""
        lines = make_model()
        log = get_logger(level='warning', encoding='utf-8')
        model = read_abaqus(lines, log=log, debug=False)
        str(model)
        model.write('spike.inp')
        os.remove('spike.inp')

        abaqus_filename = os.path.join(MODEL_PATH, 'abaqus_out.inp')
        with open(abaqus_filename, 'w') as abaqus_file:
            abaqus_file.writelines('\n'.join(lines))

    def test_abaqus_2(self):
        """two hex blocks with duplicate node ids"""
        abaqus_filename = os.path.join(MODEL_PATH, 'single_block.inp')
        log = get_logger(level='error', encoding='utf-8')

        model = read_abaqus(abaqus_filename, log=log, debug=False)
        str(model)
        model.write('spike.inp')
        os.remove('spike.inp')

    def _test_abaqus_3(self):  # pragma: no cover
        log = get_logger(level='debug', encoding='utf-8')
        lines = [
            '*part, name=test\n'
            '*NODE, NSET=NALL\n'
            '1,1.000000e+00,1.000000e+01,0.000000e+00\n'
            '2,1.000000e+00,1.000000e+01,1.000000e+02\n'
            '3,1.000000e+00,0.000000e+00,0.000000e+00\n'
            '*ELEMENT,TYPE=C3D10,ELSET=C3D10\n'
            '    82543,   13582,   49204,   12456,   13462,  182815,  182816,  100651,  100652,  182817,  100262,\n'
            '    82544,   50978,   16525,   15771,   15591,  182819,  109464,  182818,  182820,  109463,  107203,\n'
            '*NSET,NSET=FIX\n'
            '    7,\n'
            '    3,\n'
            '    1,\n'
            '    5,\n'
            '*NSET,NSET=LOAD\n'
            '    6,\n'
            '    2,\n'
            '    4,\n'
            '    8,\n'
            '*end part\n'
            '*BOUNDARY\n'
            'FIX,1,3\n'
            '*MATERIAL,NAME=EL\n'
            '*ELASTIC\n'
            '210000.,.3\n'
            '*SOLID SECTION,ELSET=C3D10,MATERIAL=EL\n'
            '*STEP\n'
            '*STATIC\n'
            '*CLOAD\n'
            'LOAD,2,-25\n'
            '*EL FILE\n'
            'U,S\n'
            '*NODE PRINT, NSET=LOAD\n'
            'U\n'
            '*END STEP\n'
        ]
        abaqus_file = StringIO()
        abaqus_file.writelines(lines)
        abaqus_file.seek(0)
        model = read_abaqus(abaqus_file, log=log, debug=True)
        str(model)
        del model

    def test_abaqus_to_nastran_5(self):
        """convert to nastran small field"""
        abaqus_filename = os.path.join(MODEL_PATH, 'solid2.inp')
        log = get_logger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'solid2.bdf')

        argv = ['format_converter', 'abaqus', abaqus_filename,
                'nastran', bdf_filename, '--encoding', 'utf-8-sig']
        cmd_line_format_converter(argv=argv, quiet=True, log=log)
        os.remove(bdf_filename)

    def test_abaqus_to_nastran_6(self):
        """convert to nastran small field"""
        abaqus_filename = os.path.join(MODEL_PATH, 'test_xform.inp')
        log = get_logger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'test_xform.bdf')

        argv = ['format_converter', 'abaqus', abaqus_filename,
                'nastran', bdf_filename, '--encoding', 'utf-8-sig']
        cmd_line_format_converter(argv=argv, quiet=True, log=log)
        os.remove(bdf_filename)

def make_model():
    """makes a test model"""
    dummy_part1 = _make_part('dummy1')
    dummy_part2 = _make_part('dummy2')
    assembly_name = 'combo'
    lines = dummy_part1 + dummy_part2 + [
        '*material,elastic,name=steel',
        #'*elastic',
        #'*plastic',
        '*user material,constants=3',
        '1,2',
        '*assembly, name=%s' % assembly_name,
        '*instance',
        '1',
        '*end instance',
        '*node',
        '1,0.,0.,0.',
        '2,1.,0.,0.',
        '3,1.,1.,0.',
        '4,0.,1.,0.',
        #'*element, type=cpe3',
        #'1,1,2,3',
        '*nset,instance=dummy1',
        '1,2,3,4,5',
        '*elset,instance=dummy1',
        '1',
        '*elset,instance=dummy2,generate',
        '1,10,3',
        '*end assembly',

        '*step',
        '*static',
        '1,3,4,2.0',
        '*end step',
    ]
    return lines

def _make_part(part_name):
    """makes a test part"""
    part = [
        '*part, name=%s' % part_name,
        '*node',
        '1,0.,0.,0.',
        '2,1.,0.,0.',
        '3,1.,1.,0.',
        '4,0.,1.,0.',
        '*element, type=cpe3',
        '1,1,2,3',
        '*element, type=cpe4',
        '2,1,2,3,4',
        '*element, type=cpe4r',
        '3,1,2,3,4',
        '*element, type=r2d2',
        '4,1,2',
        '*element, type=cax3',
        '5,1,2,3',
        '*element, type=cpe4r',
        '6,1,2,3,4',
        '*element, type=cax4r',
        '7,1,2,3,4',
        '*element, type=cohax4',
        '8,1,2,3,4',
        '*element, type=coh2d4',
        '9,1,2,3,4',
        '*mass',
        'mass_str',
        '*rotary inertia',
        'rotary_inertia_str',
        '*elset,elset=eset',
        '1',
        '*mass',
        'mass_str',
        '*solid section,material=steel',
        '42',
        #''
        '*end part',
        #'*material, name=steel',
        #'42',
    ]
    return part

if __name__ == '__main__':  #  pragma: no cover
    unittest.main()
