class EQEXIN:
    def __init__(self, nid, dof, doftype):
        self.nid = nid
        self.dof = dof
        self.doftype = doftype

    def get_stats(self, short=True):
        return str(self)

    def __repr__(self):
        return 'EQEXIN(nid, ndof, doftype); nnodes=%s\n' % len(self.nid)

