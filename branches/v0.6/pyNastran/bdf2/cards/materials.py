from .mat1 import MAT1

class Materials(object):
    def __init__(self, model):
        self.model = model

        self._mat1 = []
        self._mat1_comments = []
        self.mat1 = MAT1(model)
    
    def add_mat1(self, card, comment):
        self._mat1.append(card)
        self._mat1_comments.append(comment)

    def build(self):
        self.mat1.build(self._mat1)
        #self.mat1s.build(self._mat1s)
        #self.mat2.build(self._mat2)
        #self.mat8.build(self._mat8)

    def get_stats(self):
        msg = []
        types = self._get_types()
        for mat in types:
            nmat = len(mat.mid)
            if nmat:
                msg.append('  %-8s: %i' % (mat.type, nmat))
        return msg

    def __getitem__(self, mid):
        assert isinstance(mid, int), 'mid=%r' % mid
        types = self._get_types()
        for mat in types:
            if mid in mat.mid:
                return mat
        raise RuntimeError('mid=%r does not exist' % mid)

    def __len__(self):
        types = self._get_types()
        n = 0
        for mat in types:
            n += len(mat.mid)
        return n

    def __iter__(self):
        types = self._get_types()
        for mat in types:
            for mid in mat.mid:
                yield mid

    def _get_types(self):
        return [self.mat1]

    def _verify(self, xref=True):
        self.mat1._verify(xref)

    def write_bdf(self, f, size=8):
        f.write('$MATERIALS\n')
        self.mat1.write_bdf(f, size=size)
        #self.mat1s.write_bdf(f, size=size)
        #self.mat2.write_bdf(f, size=size)
        #self.mat8.write_bdf(f, size=size)