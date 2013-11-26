#from .mat1 import MAT1
#from .mats1 import MATS1
#from .mat4 import MAT4
#from .mat10 import MAT10

from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT4, MAT5, MAT8,
    MAT10, MAT11, MATS1)
from pyNastran.bdf.bdfInterface.assign_type import integer

class Materials(object):
    def __init__(self, model):
        self.model = model

        #self.mat1 = MAT1(model)
        #self.mats1 = MATS1(model)

        self.mat1 = {}
        self.mats1 = {}
        #==================
        self.mat2 = {}
        self.mat4 = {}
        self.mat8 = {}
        self.mat10 = {}
        #self.mat2 = MAT1(model)
        #self.mat4 = MAT1(model)
        #self.mat8 = MAT1(model)
        #self.mat10 = MAT1(model)

    def add_mat1(self, card, comment):
        #self.mat1.add(card, comment)
        mid = integer(card, 1, 'materialid')
        mat = MAT1(card=card, comment=comment)
        #mat.add(card=card, comment=comment)
        self.mat1[mid] = mat

    def add_mats1(self, card, comment):
        #self.mats1.add(card, comment)
        mid = integer(card, 1, 'materialid')
        mat = MATS1(card=card, comment=comment)
        #mat.add()
        self.mats1[mid] = mat

    def add_mat2(self, card, comment):
        self.mat2.add(card, comment)

    def add_mat4(self, card, comment):
        self.mat4.add(card, comment)

    def add_mat8(self, card, comment):
        self.mat8.add(card, comment)

    def add_mat10(self, card, comment):
        self.mat10.add(card, comment)

    def build(self):
        pass
        #types = self._get_types()
        #for mat in types:
        #    mat.build()
        #self.mat1.build()
        #self.mat1s.build()
        #==================
        #self.mat2.build()
        #self.mat4.build()
        #self.mat8.build()
        #self.mat10.build()

    def get_stats(self):
        msg = []
        types = self._get_types()
        names = self._get_type_names()
        for mat, name in zip(types, names):
            msg.append('  %-8s: %i' % (name, len(mat)))
            #nmat = mat.n
            #if nmat:
                #msg.append('  %-8s: %i' % (mat.type, nmat))
        return msg

    def get_shell_material(self, material_id):
        if material_id in self.mats1:
            return self.mats1[material_id]
        elif material_id in self.mat1:
            return self.mat1[material_id]
        elif material_id in self.mat2:
            return self.mat2[material_id]
        elif material_id in self.mat8:
            return self.mat2[material_id]

    def get_solid_material(self, material_id):
        if material_id in self.mats1:
            return self.mats1[material_id]
        elif material_id in self.mat11:
            return self.mat11[material_id]

        elif material_id in self.mat1:
            return self.mat1[material_id]

        elif material_id in self.mat9:
            return self.mat9[material_id]
        elif material_id in self.mat10:
            return self.mat10[material_id]
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def get_axial_material(self, material_id):
        if material_id in self.mat3:
            return self.mat3[material_id]
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def get_thermal_material(self, material_id):
        if material_id in self.mat4:
            return self.mat4[material_id]
        elif material_id in self.mat5:
            return self.mat5[material_id]
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def __getitem__(self, material_id):
        assert isinstance(material_id, int), 'material_id=%r' % material_id
        types = self._get_types()
        for mat in types:
            if material_id in mat.material_id:
                return mat
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def __len__(self):
        types = self._get_types()
        n = 0
        for mat in types:
            #n += mat.n
            n += len(mat)
        return n

    def __iter__(self):
        types = self._get_types()
        for materials in types:
            for mid, mat in materials.iteritems():
                yield mat
            #if mat.type == 'MATS1':
                #yield mat.material_id
            #else:
                #if mat.n:
                    #for mid in mat.material_id:
                        #yield mid

    def _get_types(self):
        return [self.mat1, self.mat2, self.mat4, self.mat8, self.mat10, 
                self.mats1,]

    def _get_type_names(self):
        return ['MAT1', 'MAT2', 'MAT4', 'MAT8', 'MAT10',
                'MATS1']

    def _verify(self, xref=True):
        for mid, material in sorted(self.mat1.iteritems()):
            material._verify(xref)

    def write_bdf(self, f, size=8, material_ids=None):
        f.write('$MATERIALS\n')
        types = self._get_types()
        for materials in types:
            for mid, mat in sorted(materials.iteritems()):
                if material_ids is None or mid in material_ids:
                    #mat.write_bdf(f, size=size, material_ids=material_ids)
                    f.write(str(mat))