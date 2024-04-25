"""
defines readers for BDF objects in the OP2 AXIC table
"""
#pylint: disable=C0103,R0914
from __future__ import annotations
from typing import TYPE_CHECKING
#from struct import unpack, Struct

#import numpy as np

#from pyNastran import is_release
#from pyNastran.bdf.cards.properties.mass import PMASS, NSM, NSML
#from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBEND
#from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP
#from pyNastran.bdf.cards.properties.bush import PBUSH
#from pyNastran.bdf.cards.properties.damper import PDAMP, PVISC
#from pyNastran.bdf.cards.properties.properties import PFAST, PGAP
#from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
#from pyNastran.bdf.cards.properties.shell import PSHEAR, PSHELL, PCOMP
#from pyNastran.bdf.cards.properties.solid import PSOLID
#from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

#from pyNastran.bdf.cards.thermal.thermal import PCONV, PHBDY
# PCOMPG, PBUSH1D, PBEAML, PBEAM3
#from pyNastran.bdf.cards.axisymmetric.axisymmetric import (
    #AXIF, RINGFL,
    #AXIC, RINGAX, POINTAX, CCONEAX, PCONEAX, PRESAX, TEMPAX,)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2_geom import OP2Geom


class AXIC:
    """defines methods for reading op2 properties"""

    def read_axic_4(self, data: bytes, ndata: int):
        """reads the AXIC table"""
        return self.op2._read_geom_4(self.axic_map, data, ndata)

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def read_fake(self, data: bytes, n: int) -> int:
        return self.op2.read_fake(data, n)

    def __init__(self, op2: OP2Geom):
        self.op2 = op2
        self.axic_map = { # per NX
            (515, 5, 144) : ['AXIC', self.read_fake], # per NX
            (2115, 21, 156) : ['FORCEAX', self.read_fake], # per NX
            (2315, 23, 146) : ['CCONEAX', self.read_fake], # per NX
            (5615, 56, 145) : ['RINGAX', self.read_fake], # per NX
        }
