from numpy import zeros, where, array, unique, concatenate

#from pyNastran.dev.bdf_vectorized.cards.materials.mat1 import MAT1
#from pyNastran.dev.bdf_vectorized.cards.materials.mats1 import MATS1
#from pyNastran.dev.bdf_vectorized.cards.materials.mat4 import MAT4
#from pyNastran.dev.bdf_vectorized.cards.materials.mat10 import MAT10

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.dev.bdf_vectorized.utils import slice_to_iter

#from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT4, MAT5, MAT8,
    #MAT10, MAT11) #, MATS1)


class Materials:
    def __init__(self, model):
        self.model = model
        self.n = 0
        self.is_built = False

        self.mat1 = model.mat1
        self.mats1 = model.mats1
        #self.mat2 = model.mat2
        #self.mat4 = model.mat4
        #self.mat5 = model.mat5
        self.mat8 = model.mat8
        #self.mat10 = model.mat10
        #self.mat11 = model.mat11

    #def add_mat1(self, card, comment):
        #self.model.log.debug('adding mat1')
        #self.mat1.add(card=card, comment=comment)

    #def add_mats1(self, card, comment):
        #self.mats1.add(card, comment)

    #def add_mat2(self, card, comment):
        #self.mat2.add(card, comment)

    #def add_mat4(self, card, comment):
        #self.mat4.add(card, comment)

    #def add_mat5(self, card, comment):
        #self.mat5.add(card, comment)

    #def add_mat8(self, card, comment):
        #self.mat8.add(card, comment)

    #def add_mat10(self, card, comment):
        #self.mat10.add(card, comment)

    #def add_mat11(self, card, comment):
        #self.mat11.add(card, comment)

    def allocate(self, card_count):
        self.model.log.info('allocate Materials')
        n = 0
        types = self._get_types()
        names = self._get_type_names()
        for name, mat in zip(names, types):
            if name in card_count:
                self.model.log.info('    allocate %s->%s' % (name, card_count[name]))
                n = card_count[name]
                mat.allocate(n)
                self.n += n
        self.n = n
        #aaa

    def build(self):
        self.is_built = True
        n = 0
        types = self._get_types()
        names = self._get_type_names()
        self.model.log.debug('building Materials')
        for name, mat in zip(names, types):
            ni = len(mat)
            n += ni
            mat.build()
            self.model.log.debug('    building %s; n=%s' % (mat.__class__.__name__, ni))
        self.n = n

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

    def get_material_id_by_property_id(self, property_id):
        n = len(property_id)
        types = self._get_types(nlimit=True)
        _material_id = concatenate([ptype.material_id for ptype in types])
        self.model.log.debug(_material_id)
        assert _material_id.shape == (n, ), _material_id.shape
        return _material_id

    #def get_density_by_material_id(self, material_id):
        #rho = zeros(len(material_id), dtype='float64')
        #for i, mid in enumerate(material_id):
            #mat = self.get_structural_material(mid)
            #rho[i] = mat.rho
        #return rho

    def get_density_by_material_id(self, material_id):
        int_flag = True if isinstance(material_id, integer_types) else False
        if isinstance(material_id, integer_types):
            material_id = array([material_id], dtype='int32')

        n = len(material_id)
        self.model.log.debug('material_id=%s' % material_id)
        umids = unique(material_id)
        umids_found = zeros(len(umids), dtype='int32')

        density = zeros(n, dtype='float64')
        mat_types = [
            self.mat1,
            #self.mat2,
            self.mat8,
        ]
        for mat_type in mat_types:
            for iumid, mid in enumerate(umids):
                if mid in mat_type.material_id:
                    rho = mat_type.get_density_by_material_id([mid])
                    i = where(mid == material_id)[0]
                    density[i] = rho
                    umids_found[iumid] = 1
                #else:
                    #print('mid=%s not found in Type=%s materials=%s' % (mid, mat_type.type, mat_type.material_id))
        assert density.shape == (n, ), density.shape

        if umids_found.min() == 0:
            izero = where(umids_found == 0)[0]
            msg = 'mids=%s not found' % umids[izero]
            raise RuntimeError(msg)
        return density[0] if int_flag else density

    def get_density_E_by_material_id(self, material_id):
        assert self.is_built, self.is_built
        int_flag = True if isinstance(material_id, integer_types) else False
        n = len(material_id)
        density = zeros(n, dtype='float64')
        E = zeros(n, dtype='float64')

        for i, mid in enumerate(material_id):
            mat = self.get_structural_material(mid)
            density[i] = mat.rho
            E[i] = mat.E()
        assert density.shape == (n, ), density.shape
        assert E.shape == (n, ), E.shape
        return density[0], E[0] if int_flag else density, E

    #def _get_param_by_material_id(self, material_id, name):
        #n = len(material_id)
        #nsm = zeros(n)
        #for i, mid in enumerate(material_id):
            #mat = self.get_structural_material(mid)
            #nsm[i] = getattr(mat, name)
        #assert nsm.shape == (n, ), nsm.shape
        #return nsm

    def get_nonstructural_mass_by_material_id(self, material_id):
        #return _get_param_by_material_id(material_id, 'nsm')
        n = len(material_id)
        nsm = zeros(n)
        for i, mid in enumerate(material_id):
            mat = self.get_structural_material(mid)
            nsm[i] = mat.nsm
        assert nsm.shape == (n, ), nsm.shape
        return nsm

    def get_E_by_material_id(self, material_id):
        #return _get_param_by_material_id(material_id, 'E')
        n = len(material_id)
        E = zeros(n)
        for i, mid in enumerate(material_id):
            mat = self.get_shell_material(mid)
            E[i] = mat.get_E_by_material_id(mid)
        assert E.shape == (n, ), E.shape
        return E

    def get_G_by_material_id(self, material_id):
        n = len(material_id)
        G = zeros(n)
        for i, mid in enumerate(material_id):
            mat = self.get_shell_material(mid)
            G[i] = mat.get_G_by_material_id(mid)
        assert G.shape == (n, ), G.shape
        return G

    #def get_material_index_by_material_id(self, property_id):
        #if property_id is None:
            #return arange(self.n)
        #return searchsorted(property_id, self.property_id)

    def get_structural_material(self, material_id):
        if material_id in self.mat1:
            return self.mat1.slice_by_material_id(material_id)
        elif material_id in self.mat8:
            return self.mat2.slice_by_material_id(material_id)
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def get_shell_material(self, material_id):
        #if material_id in self.mats1.material_id:
            #return self.mats1.slice_by_material_id(material_id)
        if material_id in self.mat1.material_id:
            return self.mat1.slice_by_material_id(material_id)
        #elif material_id in self.mat2.material_id:
            #return self.mat2.slice_by_material_id(material_id)
        elif material_id in self.mat8.material_id:
            return self.mat8.slice_by_material_id(material_id)
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def get_solid_material(self, material_id):
        #if material_id in self.mats1:
            #return self.mats1.slice_by_material_id(material_id)
        #if material_id in self.mat11:
            #return self.mat11.slice_by_material_id(material_id)
        if material_id in self.mat1.material_id:
            return self.mat1.slice_by_material_id(material_id)
        elif material_id in self.mat9.material_id:
            return self.mat9.slice_by_material_id(material_id)
        elif material_id in self.mat10.material_id:
            return self.mat10.slice_by_material_id(material_id)
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def get_axial_material(self, material_id):
        if material_id in self.mat3.material_id:
            return self.mat3.slice_by_material_id(material_id)
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def get_thermal_material(self, material_id):
        if material_id in self.mat4.material_id:
            return self.mat4.slice_by_material_id(material_id)
        elif material_id in self.mat5.material_id:
            return self.mat5.slice_by_material_id(material_id)
        raise RuntimeError('Could not find material_id=%r' % material_id)

    def __getitem__(self, material_id):
        TypeMap = {
            'MAT1'  : (self.mat1, self.mat1.material_id),
            #'MAT2'  : (self.mat2, self.mat2.material_id),
            #'MAT4'  : (self.mat4, self.mat4.material_id),
            #'MAT5'  : (self.mat5, self.mat5.material_id),
            'MAT8'  : (self.mat8, self.mat8.material_id),
            #'MAT10' : (self.mat10, self.mat10.material_id),

            #'MATS1' : (self.mats1, self.mats1.material_id),
        }
        material_id, int_flag = slice_to_iter(material_id)
        out = []
        self.model.log.debug('material_ids = %s' % material_id)
        for mid in material_id:
            obj = None
            for Type, (mdict, mids) in TypeMap.items():
                if mdict.n:
                    self.model.log.debug('Type=%s mid=%s mids=%s' % (Type, mid, mids))
                    if mid in mids:
                        self.model.log.debug(' *mid=%i Type=%s was found' % (mid, Type))
                        i = where(mid == mids)[0]
                        obj = mdict.slice_by_material_id([mid])
                        break
                    else:
                        self.model.log.debug('  mid=%i Type=%s was not found' % (mid, Type))
            out.append(obj)
        return out[0] if int_flag else out

    #def __getitem2__(self, material_id):
        #assert isinstance(material_id, integer_types), 'material_id=%r' % material_id
        #types = self._get_types()
        #for mat in types:
            #if material_id in mat.material_id:
                #return mat
        #raise RuntimeError('Could not find material_id=%r' % material_id)

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
            if materials.n:
                for mid, mat in materials.items():
                    yield mat

    def _get_types(self):
        return [self.mat1,
                #self.mat2, self.mat4,
                self.mat8,
                #self.mat10,
                self.mats1,
                ]

    def _get_type_names(self):
        return ['MAT1',
                #'MAT2', 'MAT4',
                'MAT8',
                #'MAT10',
                'MATS1',
                'MATHP',
                ]

    def _verify(self, xref=True):
        for mid, material in sorted(self.mat1.items()):
            material._verify(xref)

    def write_card(self, bdf_file, size=8, is_double=False, material_id=None):
        self.model.log.debug('Materials.write_card.n=%s' % self.n)
        if self.n:
            bdf_file.write('$MATERIALS\n')
            types = [
                self.mat1,
                #self.mat2, self.mat4,
                self.mat8,
                #self.mat10,
                self.mats1,
            ]
            for materials in types:
                if materials.n:
                    self.model.log.debug('    writing %s' % materials.type)
                    materials.write_card(bdf_file, size=size, material_id=material_id)
