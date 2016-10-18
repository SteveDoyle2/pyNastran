from pyNastran.bdf.dev_vectorized.bdf_interface.attributes2 import BDFAttributes


class WriteMesh(BDFAttributes):
    """defines methods for writing cards"""
    def __init___(self):
        """creates methods for writing cards"""
        BDFAttributes.__init__(self)

    def write_bdf(bdf_filename, interspersed=True):
        """writes a BDF"""
        pass
