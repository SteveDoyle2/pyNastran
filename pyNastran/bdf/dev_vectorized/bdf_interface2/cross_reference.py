from __future__ import print_function

from pyNastran.bdf.dev_vectorized.bdf_interface2.attributes import BDFAttributes


class CrossReference(BDFAttributes):
    """defines methods for writing cards"""

    def __init___(self):
        """creates methods for writing cards"""
        BDFAttributes.__init__(self)

    def cross_reference(self, xref=True):
        pass

    def uncross_reference(self):
        pass
