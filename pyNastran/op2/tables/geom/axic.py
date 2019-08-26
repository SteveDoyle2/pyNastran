"""
defines readers for BDF objects in the OP2 AXIC table
"""
#pylint: disable=C0103,R0914
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
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class AXIC(GeomCommon):
    """defines methods for reading op2 properties"""

    def _read_axic_4(self, data, ndata):
        """reads the AXIC table"""
        return self._read_geom_4(self._axic_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._axic_map = { # per NX
            (515, 5, 144) : ['AXIC', self._read_fake], # per NX
            (2115, 21, 156) : ['FORCEAX', self._read_fake], # per NX
            (2315, 23, 146) : ['CCONEAX', self._read_fake], # per NX
            (5615, 56, 145) : ['RINGAX', self._read_fake], # per NX
        }
