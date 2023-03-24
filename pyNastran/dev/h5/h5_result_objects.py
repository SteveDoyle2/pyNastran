from typing import Any
import numpy as np
import vtk
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk # , numpy_to_vtkIdTypeArray


class Table:
    def __init__(self,
                 name: str,
                 itime: int, iresult: int, iresults: list[int],
                 domain: int, position: int, length: int):
        # unique name
        self.name = name
        # what is the current time step
        self.itime = itime
        # what is the current result index
        self.iresult = iresult
        # get the other result cases
        self.iresults = iresults

        self.domain = domain
        self.position = position
        self.length = length



class RealVectorTable(Table):
    def __init__(self,
                 name: str,
                 itime: int, iresult: int, iresults: list[int],
                 domain: int, position: int, length: int,
                 ID, TX, TY, TZ, RX, RY, RZ, DOMAIN,
                 location: str):
        Table.__init__(self, name, itime, iresult, iresults,
                       domain, position, length)
        # location of results (node/element)
        self.location = location

        self.ID = ID
        self.TX = TX
        self.TY = TY
        self.TZ = TZ

        self.RX = RX
        self.RY = RY
        self.RZ = RZ
        self.DOMAIN = DOMAIN

    def headers(self) -> list[str]:
        return ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']

    def get_results(self) -> tuple[list[str], list[Any]]:
        txyz = np.stack([self.TX, self.TY, self.TZ], axis=1)
        txyz_array = numpy_to_vtk(txyz, deep=1, array_type=None)
        return [self.name], [txyz_array]

    def __repr__(self) -> str:
        msg = (
            'HDF5Table:\n'
            f'  itime: {self.itime}\n'
            f'  iresults: {self.iresults}\n'
            f'  domain: {self.domain}\n'
            f'  location: {self.location}\n'
            #'  itime: {itime}\n'
        )
        return msg

class RealVectorTableOptistruct(Table):
    def __init__(self,
                 name: str,
                 itime: int, iresult: int, iresults: list[int],
                 domain: int, position: int, length: int,
                 ID, VALUE, DOMAIN: int,
                 location: str):
        Table.__init__(self, name, itime, iresult, iresults,
                       domain, position, length)
        # location of results (node/element)
        self.location = location
        assert VALUE.shape

        self.ID = ID
        #self.VALUE = VALUE # np.ascontiguousarray(VALUE)
        self.VALUE = VALUE[:, :3].copy()
        self.DOMAIN = DOMAIN
        assert self.location == 'node', self.location

    def get_results(self) -> list[Any]:
        print(self)
        #txyz = self.VALUE[:, :3].copy()
        txyz = self.VALUE
        assert txyz.shape[1] == 3, txyz.shape
        txyz_array = numpy_to_vtk(txyz, deep=1, array_type=None)
        return [self.name], [txyz_array]

    def headers(self) -> list[str]:
        return ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']

    def __repr__(self) -> str:
        msg = (
            'RealVectorTableOptistruct:\n'
            f'  itime: {self.itime}\n'
            f'  iresults: {self.iresults}\n'
            f'  domain: {self.domain}\n'
            f'  location: {self.location}\n'
            #'  itime: {itime}\n'
        )
        return msg

class RealStrainEnergyOptistruct(Table):
    def __init__(self,
                 name: str,
                 itime: int, iresult: int, iresults: list[int],
                 domain: int, position: int, length: int,
                 EID, ENERGY, PERCENT, DENSITY, DOMAIN: int,
                 location: str):
        #assert name == 'STRAIN_ELEM', name
        #name = 'Strain Energy'
        Table.__init__(self, name, itime, iresult, iresults,
                       domain, position, length)
        # location of results (node/element)
        self.location = location
        self.eids = None

        # TODO: is this contiguous?
        #self.VALUE = VALUE # np.ascontiguousarray(VALUE)
        self.EID = EID
        self.ENERGY = ENERGY
        self.PERCENT = PERCENT
        self.DENSITY = DENSITY
        self.DOMAIN = DOMAIN
        assert self.location == 'element', self.location

    def get_results(self) -> list[Any]:
        print(self)
        isort = self.eids.argsort()
        eids_sorted = self.eids[isort]

        isort_EID = self.EID.argsort()
        i = isort[isort_EID]

        #txyz = self.VALUE[:, :3].copy()
        names = [f'{self.name}: {header}' for header in self.headers()]
        results = [numpy_to_vtk(arrayi, deep=1, array_type=None)
                   for arrayi in (self.ENERGY, self.PERCENT, self.DENSITY)]
        return names, results

    def headers(self) -> list[str]:
        return ['Energy', 'Percent', 'Density']

    def __repr__(self) -> str:
        msg = (
            'RealStrainEnergyOptistruct:\n'
            f'  itime: {self.itime}\n'
            f'  iresults: {self.iresults}\n'
            f'  domain: {self.domain}\n'
            f'  location: {self.location}\n'
            #'  itime: {itime}\n'
        )
        return msg


class StressTensor3D:
    """
    >>> M = StressTensor(SX, SY, SZ, TXY, TYZ, TZX)
    >>> M.SX
    >>> Mmax, Mmin, Mmid = M.principal()
    >>> Mmax[i]

    # slice first for faster access
    >>> M[i].principal()
    """
    def __init__(self, SX, SY, SZ, TXY, TYZ, TZX):
        self.SX = SX
        self.SY = SY
        self.SZ = SZ
        self.TXY = TXY
        self.TYZ = TYZ
        self.TZX = TZX
    def __getitem__(self, i):
        return StressTensor3D(self.SX[i], self.SY[i], self.SZ[i],
                              self.TXY[i], self.TYZ[i], self.TZX[i])

class ShellStressTable:
    def __init__(self,
                 itime: int, iresult: int, iresults: list[int],
                 domain: int, position: int, length: int,
                 EID,
                 FD1, X1, Y1, XY1,
                 FD2, X2, Y2, XY2,
                 DOMAIN,
                 etype: str):
        Table.__init__(self, itime, iresult, iresults,
                       domain, position, length)
