class GPDT:
    def __init__(self, nid_cp_cd_ps, xyz):
        self.nid_cp_cd_ps = nid_cp_cd_ps
        self.xyz = xyz

    def get_stats(self, short=True):
        return str(self)

    def __repr__(self):
        return 'GPDT(nid_cp_cd_ps, xyz); nnodes=%s\n' % (self.xyz.shape[0])


class BGPDT:
    def __init__(self, cd, xyz):
        self.nid_cp_cd_ps = cd
        self.xyz = xyz

    def get_stats(self, short=True):
        return str(self)

    def __repr__(self):
        return 'BGPDT(cd, xyz); nnodes=%s\n' % (self.xyz.shape[0])
