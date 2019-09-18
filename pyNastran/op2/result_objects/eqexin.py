class EQEXIN:
    def __init__(self, nid, dof, doftype):
        self.nid = nid
        self.dof = dof
        self.doftype = doftype

    def write_csv(self, file_obj):
        file_obj.write('# nnodes=%s\n' % len(self.nid))
        file_obj.write('# nid, dof, dof_type\n')
        for nid, dof, dof_type in zip(self.nid, self.dof, self.doftype):
            file_obj.write('%s, %s, %s\n' % (nid, dof, dof_type))

    def get_stats(self, short=True):
        return str(self)

    def __repr__(self):
        return 'EQEXIN(nid, ndof, doftype); nnodes=%s\n' % len(self.nid)

