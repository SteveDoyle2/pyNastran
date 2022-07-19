"""Defines the MatrixDict class"""
#import numpy as np
#from pyNastran.op2.op2_interface.write_utils import export_to_hdf5
#from pyNastran.utils import object_attributes, object_methods
from pyNastran.op2.op2_interface.op2_codes import MSC_ELEMENTS


class MatrixDict:
    """storage object for KDICT, MDICT, BDICT, etc. is op2.matdicts"""
    def __init__(self, name):
        self.name = name
        self.element_types = []
        self.numwides = []
        self.numgrids = []
        self.dof_per_grids = []

        self.eids = []
        self.ge = []
        self.address = []
        self.forms = []
        self.sils = []
        self.xforms = []

    def add(self, eltype, numwids, numgrid, dof_per_grid, form,
            eids, ge, address, sil, xform=None):
        """Sets the next set of the KDICT"""
        self.element_types.append(eltype)
        self.numwides.append(numwids)
        self.numgrids.append(numgrid)
        self.dof_per_grids.append(dof_per_grid)
        self.forms.append(form)

        self.eids.append(eids)
        self.ge.append(ge)
        self.address.append(address)
        self.sils.append(sil)
        self.xforms.append(xform)

    #@property
    #def nodes(self):
        #return [sil // 10 for sil in self.sils]

    #@property
    #def dofs(self):
        #return [sil % 10 for sil in self.sils]

    @property
    def nelements(self):
        return sum([len(eids) for eids in self.eids])

    @property
    def element_names(self):
        return [MSC_ELEMENTS[etype] for etype in self.element_types]

    def __repr__(self):
        msg = 'MatrixDict(name=%r, nelements=%s element_types=%s, element_names=[%s])' % (
            self.name, self.nelements, self.element_types, ', '.join(self.element_names))
        return msg
