from .mat1 import MAT1
#from .mat4 import MAT4
#from .mat10 import MAT10

class Materials(object):
    def __init__(self, model):
        self.model = model

        self.mat1 = MAT1(model)
        #==================
        
        self.mat2 = MAT1(model)
        self.mat4 = MAT1(model)
        self.mat8 = MAT1(model)
        self.mat10 = MAT1(model)
    
    def add_mat1(self, card, comment):
        self.mat1.add(card, comment)

    def add_mat2(self, card, comment):
        self.mat2.add(card, comment)

    def add_mat4(self, card, comment):
        self.mat4.add(card, comment)

    def add_mat8(self, card, comment):
        self.mat8.add(card, comment)

    def add_mat10(self, card, comment):
        self.mat10.add(card, comment)

    def build(self):
        self.mat1.build()
        #==================
        #self.mat1s.build()
        self.mat2.build()
        self.mat4.build()
        self.mat8.build()
        self.mat10.build()

    def get_stats(self):
        msg = []
        types = self._get_types()
        for mat in types:
            nmat = mat.n
            if nmat:
                msg.append('  %-8s: %i' % (mat.type, nmat))
        return msg

    def __getitem__(self, mid):
        assert isinstance(mid, int), 'mid=%r' % mid
        types = self._get_types()
        for mat in types:
            if mid in mat.material_id:
                return mat
        raise RuntimeError('mid=%r does not exist' % mid)

    def __len__(self):
        types = self._get_types()
        n = 0
        for mat in types:
            n += mat.n
        return n

    def __iter__(self):
        types = self._get_types()
        for mat in types:
            if mat.n:
                for mid in mat.material_id:
                    yield mid

    def _get_types(self):
        return [self.mat1, self.mat2, self.mat4, self.mat8, self.mat10]

    def _verify(self, xref=True):
        self.mat1._verify(xref)

    def write_bdf(self, f, size=8, mids=None):
        f.write('$MATERIALS\n')
        types = self._get_types()
        for mat in types:
            if mat.n:
                self.mat1.write_bdf(f, size=size, mids=mids)
