"""
ResultData is the target format when writing results to the database.
Punch results are converted to this format on the fly.
F06 results would be too, but not yet implemented.
OP2 results could simply target this format in the beginning.
"""
from __future__ import print_function, absolute_import


class ResultData(object):
    def __init__(self):
        self._is_elemental = False
        self._is_nodal = False
        self.element_code = -1  # e.g. linear cquad4 = 33
        self._is_real = False
        self._is_complex = False
        self._is_random = False
        
        self._is_elemental_forces = False
        self._is_elemental_stresses = False
        self._is_elemental_strains = False
        
        self._is_displacements = False
        self._is_grid_forces = False
        self._is_oloads = False
        self._is_mpcf = False
        self._is_spcf = False

        # data holds the actual data in a dict-like object (preferably an ndarray)
        # all values must have the same length
        # keys should be the same as what is defined in the MSC spec for the table columns
        self.data = None

        # options are mostly for punch/f06 results... for example if shell loading is in material system
        # options should contain the value 'MATERIAL'... the loading will be converted to the element system before
        # writing
        # that's the only option implemented for now
        self.options = set()
        # where does the data come from?  punch, f06, or op2
        self.source = None
        
        self.subcase_id = None
        
    def set_result_type(self, result_type):
        if 'ELEMENT' in result_type:           
            if 'STRESS' in result_type:
                self.elemental_stresses()
            elif 'STRAIN' in result_type:
                self.elemental_strains()
            elif 'FORCE' in result_type:
                self.elemental_forces()
                
            tmp = result_type.split()
            
            elem_code = -1
            
            for i in tmp:
                try:
                    elem_code = int(i)
                except ValueError:
                    pass
                
            if elem_code == -1:
                raise Exception('Unknown result type! %s' % result_type)
            
            self.element_code = elem_code
            
        else:
            if 'DISPLACEMENTS' in result_type:
                self.displacements()
            elif 'GRID POINT FORCE BALANCE' in result_type:
                self.grid_forces()
            elif 'OLOADS':
                self.oloads()
            elif 'MPCF' in result_type:
                self.mpcf()
            elif 'SPCF' in result_type:
                self.spcf()
            else:
                raise Exception('Unknown result type! %s' % result_type)
        
        if 'REAL' in result_type:
            self.real()
        elif 'COMPLEX' in result_type:
            self.complex()
        elif 'RANDOM' in result_type:
            self.random()
        else:
            raise Exception('Unknown result type! %s' % result_type)

    @property
    def result_type(self):
        result_type = []
        if self._is_elemental:
            result_type.append('ELEMENT')
        elif self._is_nodal:
            if self._is_displacements:
                result_type.append('DISPLACEMENTS')
            elif self._is_grid_forces:
                result_type.append('GRID POINT FORCE BALANCE')
            elif self._is_oloads:
                result_type.append('OLOADS')
            elif self._is_mpcf:
                result_type.append('MPCF')
            elif self._is_spcf:
                result_type.append('SPCF')
        
        if len(result_type) == 0:
            raise Exception('Unknown result type!')
        
        if self._is_elemental:
            if self._is_elemental_forces:
                result_type.append('FORCES')
            elif self._is_elemental_stresses:
                result_type.append('STRESSES')
            elif self._is_elemental_strains:
                result_type.append('STRAINS')
            else:
                raise Exception('Unknown result type!')
            
            if self.element_code <= 0:
                raise Exception('Unknown element code! %d' % self.element_code)
            
            result_type.append(str(self.element_code))
            
        if self._is_real:
            result_type.append('REAL')
        elif self._is_complex:
            result_type.append('COMPLEX')
        elif self._is_random:
            result_type.append('RANDOM')
        else:
            raise Exception('Unknown result type!')
        
        return ' '.join(result_type)
        
    def displacements(self):
        self.nodal()
        self._is_displacements = True
        self._is_grid_forces = False
        self._is_oloads = False
        self._is_mpcf = False
        self._is_spcf = False
        
    def grid_forces(self):
        self.nodal()
        self._is_displacements = False
        self._is_grid_forces = True
        self._is_oloads = False
        self._is_mpcf = False
        self._is_spcf = False
        
    def oloads(self):
        self.nodal()
        self._is_displacements = False
        self._is_grid_forces = False
        self._is_oloads = True
        self._is_mpcf = False
        self._is_spcf = False
        
    def mpcf(self):
        self.nodal()
        self._is_displacements = False
        self._is_grid_forces = False
        self._is_oloads = False
        self._is_mpcf = True
        self._is_spcf = False
        
    def spcf(self):
        self.nodal()
        self._is_displacements = False
        self._is_grid_forces = False
        self._is_oloads = False
        self._is_mpcf = False
        self._is_spcf = True
        
    def elemental(self):
        self._is_elemental = True
        self._is_nodal = False
        
    def nodal(self):
        self._is_nodal = True
        self._is_elemental = False
        
    def real(self):
        self._is_real = True
        self._is_complex = False
        self._is_random = False
        
    def complex(self):
        self._is_complex = True
        self._is_real = False
        self._is_random = False
        
    def random(self):
        self._is_random = True
        self._is_real = False
        self._is_complex = False
        
    def elemental_forces(self):
        self.elemental()
        self._is_elemental_forces = True
        self._is_elemental_stresses = False
        self._is_elemental_strains = False
        
    def elemental_stresses(self):
        self.elemental()
        self._is_elemental_stresses = True
        self._is_elemental_forces = False
        self._is_elemental_strains = False
        
    def elemental_strains(self):
        self.elemental()
        self._is_elemental_strains = True
        self._is_elemental_stresses = False
        self._is_elemental_forces = False
        
    def punch_results(self):
        self.source = 'punch'
        
    def f06_results(self):
        self.source = 'f06'
        
    def op2_results(self):
        self.source = 'op2'
