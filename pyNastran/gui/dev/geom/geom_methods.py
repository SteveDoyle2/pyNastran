from pyNastran.bdf.field_writer_8 import print_card_8
class GeomMethods:
    def add_hexa(name, ):
        """
        CBRICK, eid, pid, nu, nv, nw,
        G111, G222, G333, ..., G(nu, nv, nw)

        CBRICK, 1, 1, 4, 2, 3
        111, 211, 311, 411, 121, 221, 321, 421,
        112, 212, 312, 412, 122, 222, 322, 422,
        113, 213, 313, 413, 123, 223, 323, 423,
        """
        pass

    def create_point_between_points(self, p1, p2, x1=0.5):
        p3 = p1 * x1 + p2 (1 - x1)
        return p3

    def add_tri_rod(name, p1, p2, p3):
        self.trirods[name] = [p1, p2, p3]
    def trirod_to_brick(self, p1, p2, p3, width, height):
        nu = 3
        nv = 2
        nw = 2
        nids = write_cbrick(eid, pid, nids, nu, nv, nw)

def write_cbrick(eid, pid, nids, nu, nv, nw):
    """
    nids = [1, 2, 3, 4, 5, 6, 7, 8]
    """
    list_fields = ['CBRICK', eid, pid, nu, nv, nw, None, None, None]

    nids = []
    for iw in range(nw):
        for iv in range(nv):
            for iu in range(nu):
                nids.append(1 * (iu + 1) + 10 * (iv + 1) + 100 * (iw + 1))
    list_fields += nids
    return nids
    #return print_card_8(list_fields)

print(write_cbrick(1, 2, [], 3, 2, 2))
print(write_cbrick(1, 2, [], 2, 2, 2))
