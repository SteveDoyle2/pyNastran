from __future__ import annotations
import sys
from itertools import count
from collections import defaultdict
from typing import Optional, Any, cast, TYPE_CHECKING

import numpy as np

from pyNastran.dev.bdf_vectorized3.bdf_interface.add_methods import AddMethods
from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import BDFAttributes
from pyNastran.dev.bdf_vectorized3.cards.elements.bar import BAROR
from pyNastran.dev.bdf_vectorized3.cards.deqatn import DEQATN

from pyNastran.bdf.cards.dmig import DMIG, DMIG_UACCEL, DMI, DMIJ, DMIJI, DMIK, DMIAX, DTI, DTI_UNITS
from pyNastran.bdf.cards.materials import NXSTRAT
from pyNastran.bdf.cards.methods import EIGRL
from pyNastran.bdf.cards.dynamic import (
    FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
    TSTEP, TSTEP1, TSTEPNL, NLPARM, NLPCI, ROTORG, ROTORD)
from pyNastran.bdf.cards.optimization import DOPTPRM
from pyNastran.bdf.cards.bdf_tables import TABLEH1, TABLEHT
from pyNastran.bdf.field_writer_8 import print_card_8


from pyNastran.bdf.cards.aero.dynamic_loads import AERO, MKAERO1, MKAERO2
from pyNastran.bdf.cards.aero.static_loads import AEROS
from pyNastran.utils.numpy_utils import integer_types, integer_string_types
from pyNastran.bdf.cards.params import MDLPRM, PARAM
from pyNastran.dev.bdf_vectorized3.bdf import DTI_UNITS
from pyNastran.bdf.cards.contact import BCTPARA, BCTPARM
from pyNastran.bdf.cards.bdf_tables import (
    TABLED1, TABLED2, TABLED3, TABLED4,
    TABLEM1, TABLEM2, TABLEM3, TABLEM4,
    DTABLE,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArray3float, NDArray66float
    from pyNastran.dev.bdf_vectorized3.bdf import PARAM # BDF,
    #from pyNastran.dev.bdf_vectorized3.cards.grid import GRID
    #from pyNastran.dev.bdf_vectorized3.cards.coord import COORD # CORD1R, CORD1C, CORD1S, CORD2R, CORD2C, CORD2S
    #from pyNastran.dev.bdf_vectorized3.cards.loads.static_loads import LOAD, FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2, LOADSET
    #from pyNastran.dev.bdf_vectorized3.cards.loads.static_pressure_loads import PLOAD, PLOAD1, PLOAD2, PLOAD4 # , PLOADX1
    #from pyNastran.dev.bdf_vectorized3.cards.loads.dynamic_loads import (
        #DAREA, DELAY, DLOAD, DPHASE, LSEQ, QVECT, RANDPS,
        #TIC, RLOAD1, RLOAD2, TLOAD1, TLOAD2)
    #from pyNastran.dev.bdf_vectorized3.cards.elements.mass import CONM1, CONM2
    #from pyNastran.dev.bdf_vectorized3.cards.elements.plot import PLOTEL
    #from pyNastran.dev.bdf_vectorized3.cards.elements.shear import CSHEAR, PSHEAR
    #from pyNastran.dev.bdf_vectorized3.cards.elements.shell import (
        #CTRIA3, CTRIA6, CTRIAR,
        #CQUAD4, CQUAD8, CQUAD, CQUADR)
    #from pyNastran.dev.bdf_vectorized3.cards.elements.shell_axi import (
        #CTRIAX, CTRIAX6,
        #CQUADX, CQUADX4, CQUADX8)
    #from pyNastran.dev.bdf_vectorized3.cards.aero.aero import (
        #CAERO1, CAERO2, CAERO3, CAERO4, CAERO5, CAERO7,
        #PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
        #SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5,
        #MONPNT1, MONPNT2, MONPNT3,
        #AECOMP, AECOMPL, AEFACT, AELIST, AEPARM, AELINK,
        #AESURF, AESURFS,
        #CSSCHD, FLFACT, GUST)
    #from pyNastran.dev.bdf_vectorized3.cards.elements.mass import CONM1, CONM2
    #from pyNastran.dev.bdf_vectorized3.cards.elements.spring import CELAS1, CELAS2, CELAS3, CELAS4, PELAS, PELAST
    #from pyNastran.dev.bdf_vectorized3.cards.elements.damper import CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, PDAMP, PDAMPT, CGAP, CVISC, PGAP, PVISC
    #from pyNastran.dev.bdf_vectorized3.cards.elements.nsm import NSM, NSM1, NSML, NSML1, NSMADD
    #from pyNastran.dev.bdf_vectorized3.cards.elements.rod import CROD, CTUBE, CONROD, PROD, PTUBE
    #from pyNastran.dev.bdf_vectorized3.cards.elements.solid import CTETRA, CPYRAM, CPENTA, CHEXA, PSOLID, PLSOLID # , PCOMPS, PCOMPLS
    #from pyNastran.dev.bdf_vectorized3.cards.materials import MAT1, MAT2, MAT8, MAT9 # MAT3, MAT4, MAT5, , MAT10, MAT10C, MAT11


class AddCoords(BDFAttributes):
    def add_cord2r(self, cid: int,
                   origin: Optional[list[float] | NDArray3float],
                   zaxis: Optional[list[float] | NDArray3float],
                   xzplane: Optional[list[float] | NDArray3float],
                   rid: int=0, setup: bool=True, comment: str='') -> int:
        """
        Creates the CORD2R card, which defines a rectangular coordinate
        system using 3 vectors.

        Parameters
        ----------
        cid : int
            coordinate system id
        rid : int; default=0
            the referenced coordinate system that defines the system the
            vectors
        origin : list[float, float, float]
            the origin of the coordinate system
        zaxis : list[float, float, float]
            the z-axis of the coordinate system
        xzplane : list[float, float, float]
            a point on the xz plane
        comment : str; default=''
            a comment for the card

        """
        coord = self.coord.add_cord2r(cid, origin, zaxis, xzplane, rid=rid,
                                      setup=setup, comment=comment)
        return coord

    def add_cord2c(self, cid: int,
                   origin: Optional[list[float] | NDArray3float],
                   zaxis: Optional[list[float] | NDArray3float],
                   xzplane: Optional[list[float] | NDArray3float],
                   rid: int=0, setup: bool=True, comment: str='') -> int:
        """
        Creates the CORD2C card, which defines a cylindrical coordinate
        system using 3 vectors.

        Parameters
        ----------
        cid : int
            coordinate system id
        rid : int; default=0
            the referenced coordinate system that defines the system the
            vectors
        origin : list[float, float, float]
            the origin of the coordinate system
        zaxis : list[float, float, float]
            the z-axis of the coordinate system
        xzplane : list[float, float, float]
            a point on the xz plane
        comment : str; default=''
            a comment for the card

        """
        coord = self.coord.add_cord2c(
            cid, rid=rid, origin=origin, zaxis=zaxis, xzplane=xzplane,
            setup=setup, comment=comment)
        return coord

    def add_cord2s(self, cid: int,
                   origin: Optional[list[float] | NDArray3float],
                   zaxis: Optional[list[float] | NDArray3float],
                   xzplane: Optional[list[float] | NDArray3float],
                   rid: int=0, setup: bool=True, comment: str='') -> int:
        """
        Creates the CORD2C card, which defines a spherical coordinate
        system using 3 vectors.

        Parameters
        ----------
        cid : int
            coordinate system id
        origin : list[float, float, float]
            the origin of the coordinate system
        zaxis : list[float, float, float]
            the z-axis of the coordinate system
        xzplane : list[float, float, float]
            a point on the xz plane
        rid : int; default=0
            the referenced coordinate system that defines the system the
            vectors
        comment : str; default=''
            a comment for the card

        """
        coord = self.coord.add_cord2s(
            cid, rid=rid, origin=origin, zaxis=zaxis, xzplane=xzplane,
            setup=setup, comment=comment)
        return coord

    def add_cord1r(self, cid: int, g1: int, g2: int, g3: int, comment: str='') -> int:
        """
        Creates the CORD1R card, which defines a rectangular coordinate
        system using 3 GRID points.

        Parameters
        ----------
        cid : int
            the coordinate id
        g1 : int
            grid point 1
        g2 : int
            grid point 2
        g3 : int
            grid point 3
        comment : str; default=''
            a comment for the card

        """
        coord = self.coord.add_cord1r(cid, g1, g2, g3, comment=comment)
        return coord

    def add_cord1c(self, cid: int, g1: int, g2: int, g3: int, comment: str='') -> int:
        """
        Creates the CORD1C card, which defines a cylindrical coordinate
        system using 3 GRID points.

        Parameters
        ----------
        cid : int
            the coordinate id
        g1 : int
            grid point 1
        g2 : int
            grid point 2
        g3 : int
            grid point 3
        comment : str; default=''
            a comment for the card

        """
        coord = self.coord.add_cord1c(cid, g1, g2, g3, comment=comment)
        return coord

    def add_cord1s(self, cid: int, g1: int, g2: int, g3: int, comment: str='') -> int:
        """
        Creates the CORD1S card, which defines a spherical coordinate
        system using 3 GRID points.

        Parameters
        ----------
        cid : int
            the coordinate id
        g1 : int
            grid point 1
        g2 : int
            grid point 2
        g3 : int
            grid point 3
        comment : str; default=''
            a comment for the card

        """
        coord = self.coord.add_cord1s(cid, g1, g2, g3, comment=comment)
        return coord

    def add_cord3g(self, cid: int,
                   method_es, method_int, form,
                   thetas: list[int],
                   rid: int,
                   comment: str='') -> CORD3G:
        """Creates a CORD3G card"""
        coord = CORD3G(cid, method_es, method_int, form, thetas, rid, comment=comment)
        return coord

    #def add_gmcord(self, cid, entity, gm_ids, comment: str='') -> GMCORD:
        #"""Creates a GMCORD coordinate card"""
        #coord = GMCORD(cid, entity, gm_ids, comment=comment)
        #return coord


class AddBolts(BDFAttributes):
    def add_bolt_nx(self, bolt_id: int,
                    element_type: int,
                    eids: Optional[list]=None,  # element_type=1
                    nids: Optional[list]=None,  # element_type=2
                    csid=None,  # element_type=2
                    idir=None,  # element_type=2
                    comment: str='') -> int:
        bolt = self.bolt.add_nx(bolt_id, element_type,
                                eids=eids, nids=nids,
                                csid=csid, idir=idir, comment=comment)
        return bolt

    def add_boltseq_nx(self, sid: int,
                       s_nos: list[int],
                       b_ids: list[int],
                       n_incs: Optional[list[int]]=None,
                       comment: str='') -> int:
        bolt = self.boltseq.add(sid, s_nos, b_ids, n_incs=n_incs, comment=comment)
        return bolt

    def add_boltfor_nx(self, sid: int, load_value: float, bolt_ids: list[int],
                       comment: str='') -> int:
        boltfor = self.boltfor.add_nx(sid, load_value, bolt_ids, comment=comment)
        return boltfor


class Add0dElements(BDFAttributes):
    def add_pelas(self, pid: int, k: float, ge: float=0., s: float=0.,
                  comment: str='') -> int:
        """
        Creates a PELAS card

        Parameters
        ----------
        pid : int
            property id
        k : float
            spring stiffness
        ge : int; default=0.0
            damping coefficient
        s : float; default=0.0
            stress coefficient
        comment : str; default=''
            a comment for the card

        """
        prop = self.pelas.add(pid, k, ge, s, comment=comment)
        return prop

    def add_celas1(self, eid: int, pid: int, nids: list[int],
                   c1: int=0, c2: int=0, comment: str='') -> int:
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        elem = self.celas1.add(eid, pid, nids, c1, c2, comment=comment)
        return elem

    def add_celas2(self, eid: int, k: float, nids: list[int],
                   c1: int=0, c2: int=0, ge: float=0., s: float=0.,
                   comment: str='') -> int:
        """
        Creates a CELAS2 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : list[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        ge : int; default=0.0
            damping coefficient
        s : float; default=0.0
            stress coefficient
        comment : str; default=''
            a comment for the card

        """
        elem = self.celas2.add(eid, k, nids, c1=c1, c2=c2, ge=ge, s=s, comment=comment)
        return elem

    def add_celas3(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CELAS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.celas3.add(eid, pid, nids, comment=comment)
        return elem

    def add_celas4(self, eid: int, k: float, nids: list[int], comment: str='') -> int:
        """
        Creates a CELAS4 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.celas4.add(eid, k, nids, comment=comment)
        return elem

    def add_cdamp1(self, eid: int, pid: int, nids: list[int], c1: int=0, c2: int=0,
                   comment: str='') -> int:
        """
        Creates a CDAMP1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        elem = self.cdamp1.add(eid, pid, nids, c1=c1, c2=c2, comment=comment)
        return elem

    def add_cdamp2(self, eid: int, b: float, nids: list[int],
                   c1: int=0, c2: int=0, comment: str='') -> int:
        """
        Creates a CDAMP2 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : list[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        elem = self.cdamp2.add(eid, b, nids, c1=c1, c2=c2, comment=comment)
        return elem

    def add_cdamp3(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CDAMP3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.cdamp3.add(eid, pid, nids, comment=comment)
        return elem

    def add_cdamp4(self, eid: int, b: float, nids: list[int],
                   comment: str='') -> int:
        """
        Creates a CDAMP4 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.cdamp4.add(eid, b, nids, comment=comment)
        return elem

    def add_cdamp5(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CDAMP5 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP5)
        nids : list[int, int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.cdamp5.add(eid, pid, nids, comment=comment)
        return elem

    def add_pdamp(self, pid: int, b: float, comment: str='') -> int:
        """
        Creates a PDAMP card

        Parameters
        ----------
        pid : int
            property id
        b : float
            viscous damping
        comment : str; default=''
            a comment for the card

        """
        prop = self.pdamp.add(pid, b, comment=comment)
        return prop

    def add_pdampt(self, pid: int, tbid: int, comment: str='') -> int:
        """
        Creates a PDAMPT card

        Parameters
        ----------
        pid : int
            property id
        tbid : int
            TABLED1? id
        comment : str; default=''
            a comment for the card

        """
        prop = self.pdampt.add(pid, tbid, comment=comment)
        return prop

    def add_pdamp5(self, pid: int, mid: int, b: float, comment: str='') -> int:
        """Creates a PDAMP5 card"""
        prop = self.pdamp5.add(pid, mid, b, comment=comment)
        return prop

    def add_cvisc(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CVISC card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PVISC)
        nids : list[int, int]
            GRID ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.cvisc.add(eid, pid, nids, comment=comment)
        return elem

    def add_pvisc(self, pid: int, ce: float, cr: float, comment: str='') -> int:
        """
        Creates a PVISC card

        Parameters
        ----------
        pid : int
            property id for a CVISC
        ce : float
            Viscous damping values for extension in units of force per unit velocity
        cr : float
            Viscous damping values for rotation in units of moment per unit velocity.
        comment : str; default=''
            a comment for the card

        """
        prop = self.pvisc.add(pid, ce, cr, comment=comment)
        return prop

    def add_cgap(self, eid: int, pid: int, nids: list[int],
                 x: Optional[list[int]], g0: Optional[int],
                 cid: Optional[int]=None, comment: str='') -> int:
        """
        Creates a CGAP card

        Parameters
        ----------
        eid : int
            Element ID
        pid : int
            Property ID (PGAP)
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        x : list[float, float, float]
            Components of the orientation vector,
            from GA, in the displacement coordinate system at GA
        g0 : int
            GO Alternate method to supply the orientation vector using
            grid point GO. Direction of is from GA to GO
        cid : int; default=None
            Element coordinate system identification number.
            CID must be specified if GA and GB are coincident
            (distance from GA to GB < 10^-4)
        comment : str; default=''
            a comment for the card

        """
        elem = self.cgap.add(eid, pid, nids, x, g0, cid=cid, comment=comment)
        return elem

    def add_pgap(self, pid: int, u0: float=0., f0: float=0.,
                 ka: float=1.e8, kb: Optional[float]=None, mu1: float=0.,
                 kt: Optional[float]=None, mu2: Optional[float]=None,
                 tmax: float=0., mar: float=100., trmin: float=0.001,
                 comment: str='') -> int:
        """
        Defines the properties of the gap element (CGAP entry).

        Parameters
        ----------
        pid : int
            property id for a CGAP
        u0 : float; default=0.
            Initial gap opening
        f0 : float; default=0.
            Preload
        ka : float; default=1.e8
            Axial stiffness for the closed gap
        kb : float; default=None -> 1e-14 * ka
            Axial stiffness for the open gap
        mu1 : float; default=0.
            Coefficient of static friction for the adaptive gap element
            or coefficient of friction in the y transverse direction
            for the nonadaptive gap element
        kt : float; default=None -> mu1*ka
            Transverse stiffness when the gap is closed
        mu2 : float; default=None -> mu1
            Coefficient of kinetic friction for the adaptive gap element
            or coefficient of friction in the z transverse direction
            for the nonadaptive gap element
        tmax : float; default=0.
            Maximum allowable penetration used in the adjustment of
            penalty values. The positive value activates the penalty
            value adjustment
        mar : float; default=100.
            Maximum allowable adjustment ratio for adaptive penalty
            values KA and KT
        trmin : float; default=0.001
            Fraction of TMAX defining the lower bound for the allowable
            penetration
        comment : str; default=''
            a comment for the card

        """
        prop = self.pgap.add(pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin,
                             comment=comment)
        return prop

    def add_cfast(self, eid: int, pid: int, fast_type: str,
                  ida: int, idb: int,
                  gs=None, ga=None, gb=None,
                  xs=None, ys=None, zs=None, comment: str='') -> int:
        """Creates a CFAST card"""
        elem = self.cfast.add(eid, pid, fast_type, [ida, idb],
                              gs=gs, ga=ga, gb=gb,
                              xs=xs, ys=ys, zs=zs, comment=comment)
        return elem

    def add_pfast(self, pid: int, d: float,
                  kt1: float, kt2: float, kt3: float,
                  mcid: int=-1, mflag: int=0,
                  kr1: float=0., kr2: float=0., kr3: float=0.,
                  mass: float=0., ge: float=0.,
                  comment: str='') -> int:
        """
        Creates a PAST card

        Parameters
        ----------
        pid : int
            property id
        d : int
            diameter of the fastener
        kt1, kt2, kt3 : float
            stiffness values in directions 1-3
        mcid : int; default=-1
            specifies the element stiffness coordinate system
        mflag : int; default=0
            0-absolute; 1-relative
        kr1, kr2, kr3 : float; default=0.0
            rotational stiffness values in directions 1-3
        mass : float; default=0.0
            lumped mass of the fastener
        ge : float; default=0.0
            structural damping
        comment : str; default=''
            a comment for the card

        """
        prop = self.pfast.add(pid, d, kt1, kt2, kt3, mcid=mcid, mflag=mflag,
                              kr1=kr1, kr2=kr2, kr3=kr3, mass=mass, ge=ge,
                              comment=comment)
        return prop

    def add_cbush(self, eid: int, pid: int, nids,
                  x: Optional[list[float]], g0: Optional[int], cid=None,
                  s: float=0.5, ocid: int=-1,
                  si: Optional[list[float]]=None, comment: str='') -> int:
        """
        Creates a CBUSH card

        Parameters
        ----------
        eid : int
            Element id
        pid : int
            Property id (PBUSH)
        nids : list[int, int]
            node ids; connected grid points at ends A and B
            The nodes may be coincident, but then cid is required.
        x : list[float, float, float]; None
            list : the directional vector used to define the stiffnesses
                   or damping from the PBUSH card
            None : use g0
        g0 : int/None
            int : the directional vector used to define the stiffnesses
                  or damping from the PBUSH card
            None : use x
        cid : int; default=None
            Element coordinate system identification. A 0 means the basic
            coordinate system. If CID is blank, then the element coordinate
            system is determined from GO or Xi.
        s: float; default=0.5
            Location of spring damper (0 <= s <= 1.0)
        ocid : int; default=-1
            Coordinate system identification of spring-damper offset.
            (Integer > -1; Default = -1, which means the offset
            point lies on the line between GA and GB)
        si : list[float, float, float]; default=None
            Components of spring-damper offset in the OCID coordinate system
            if OCID > 0.
            None : [None, None, None]
        comment : str; default=''
            a comment for the card

        """
        elem = self.cbush.add(eid, pid, nids, x, g0, cid=cid, s=s, ocid=ocid, si=si,
                              comment=comment)
        return elem

    def add_pbush(self, pid: int, k: list[float], b: list[float], ge: list[float],
                  rcv: Optional[list[float]]=None, mass: Optional[float]=None,
                  alpha: float=0., tref: float=0., coincident_length=None,
                  comment: str='') -> int:
        """
        Creates a PBUSH card, which defines a property for a PBUSH

        Parameters
        ----------
        pid : int
            property id
        k : list[float]
            Nominal stiffness values in directions 1 through 6.
            len(k) = 6
        b : list[float]
            Nominal damping coefficients in direction 1 through 6 in units of
            force per unit velocity
            len(b) = 6
        ge : list[float]
            Nominal structural damping constant in directions 1 through 6.
            len(ge) = 6
        rcv : list[float]; default=None -> (None, None, None, None)
            [sa, st, ea, et] = rcv
            length(rcv) = 4
        mass : float; default=None
            lumped mass of the CBUSH
            This is an MSC only parameter.
        comment : str; default=''
            a comment for the card

        """
        prop = self.pbush.add(pid, k, b, ge, rcv=rcv, mass=mass, comment=comment)
        return prop

    def add_cbush1d(self, eid: int, pid: int, nids: list[int], cid: Optional[int]=None,
                    comment: str='') -> int:
        """Creates a CBUSH1D card"""
        elem = self.cbush1d.add(eid, pid, nids, cid=cid, comment=comment)
        return elem

    def add_cbush2d(self, eid: int, pid: int, nids: list[int], cid: int=0,
                    plane: str='XY', sptid: Optional[int]=None, comment: str='') -> int:
        """Creates a CBUSH2D card"""
        elem = self.cbush2d.add(eid, pid, nids, cid=cid, plane=plane, sptid=sptid, comment=comment)
        return elem

    def add_pbush1d(self, pid: int,
                    k: float=0., c: float=0., m: float=0.,
                    sa: float=0., se: float=0., optional_vars=None,
                    comment: str='') -> int:
        """Creates a PBUSH1D card"""
        prop = self.pbush1d.add(pid, k=k, c=c, m=m, sa=sa, se=se,
                                optional_vars=optional_vars, comment=comment)
        return prop

    #def add_pbush2d(self, pid, k, c, m, sa, se, optional_vars, comment: str='') -> int:
        #"""
        #Creates a PBUSH2D card
        #"""
        #prop = PBUSH2D(pid. comment=comment)
        #self.add_property_object(prop)
        #return prop

    def add_pbusht(self, pid: int,
                   k_tables: Optional[list[int]]=None,
                   b_tables: Optional[list[int]]=None,
                   ge_tables: Optional[list[int]]=None,
                   kn_tables: Optional[list[int]]=None, comment: str='') -> int:
        """Creates a PBUSHT card"""
        prop = self.pbusht.add(pid, k_tables, b_tables, ge_tables, kn_tables,
                               comment=comment)
        return prop

    def add_pelast(self, pid: int, tkid: int=0, tgeid: int=0, tknid: int=0,
                   comment: str='') -> int:
        """
        Creates a PELAST card

        Parameters
        ----------
        pid : int
            property id
        tkid : float
            TABLEDx that defines k vs. frequency
        tgeid : int; default=0
            TABLEDx that defines ge vs. frequency
        s : float; default=0.
            TABLEDx that defines force vs. displacement
        comment : str; default=''
            a comment for the card

        """
        prop = self.pelast.add(pid, tkid, tgeid, tknid, comment=comment)
        return prop

    def add_conm1(self, eid: int, nid: int, mass_matrix: NDArray66float,
                  cid: int=0, comment: str='') -> int:
        """
        Creates a CONM1 card

        Parameters
        ----------
        eid : int
            element id
        nid : int
            the node to put the mass matrix
        mass_matrix : (6, 6) float ndarray
            the 6x6 mass matrix, M
        cid : int; default=0
            the coordinate system for the mass matrix
        comment : str; default=''
            a comment for the card

        ::

          [M] = [M11 M21 M31 M41 M51 M61]
                [    M22 M32 M42 M52 M62]
                [        M33 M43 M53 M63]
                [            M44 M54 M64]
                [    Sym         M55 M65]
                [                    M66]

        """
        mass = self.conm1.add(eid, nid, mass_matrix, cid=cid, comment=comment)
        return mass

    def add_conm2(self, eid: int, nid: int, mass: float, cid: int=0,
                  X: Optional[list[float]]=None, I: Optional[list[float]]=None,
                  comment: str='') -> int:
        """
        Creates a CONM2 card

        Parameters
        ----------
        eid : int
           element id
        nid : int
           node id
        mass : float
           the mass of the CONM2
        cid : int; default=0
           coordinate frame of the offset (-1=absolute coordinates)
        X : (3, ) list[float]; default=None -> [0., 0., 0.]
            xyz offset vector relative to nid
        I : (6, ) list[float]; default=None -> [0., 0., 0., 0., 0., 0.]
            mass moment of inertia matrix about the CG
            I11, I21, I22, I31, I32, I33 = I
        comment : str; default=''
            a comment for the card

        """
        mass_obj = self.conm2.add(eid, nid, mass, cid=cid, X=X, I=I, comment=comment)
        return mass_obj

    def add_pmass(self, pid: int, mass: float, comment: str='') -> int:
        """
        Creates an PMASS card, which defines a mass applied to a single DOF

        Parameters
        ----------
        pid : int
            Property id used by a CMASS1/CMASS3 card
        mass : float
            the mass to apply
        comment : str; default=''
            a comment for the card

        """
        prop = self.pmass.add(pid, mass, comment=comment)
        return prop

    def add_cmass1(self, eid: int, pid: int, nids: list[int],
                   c1: int=0, c2: int=0, comment: str='') -> int:
        """
        Creates a CMASS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        mass_obj = self.cmass1.add(eid, pid, nids, c1, c2, comment=comment)
        return mass_obj

    def add_cmass2(self, eid: int, mass: float, nids: list[int],
                   c1: int, c2: int, comment: str='') -> int:
        """
        Creates a CMASS2 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            mass
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        mass_obj = self.cmass2.add(eid, mass, nids, c1, c2, comment=comment)
        return mass_obj

    def add_cmass3(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CMASS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        mass = self.cmass3.add(eid, pid, nids, comment=comment)
        return mass

    def add_cmass4(self, eid: int, mass: float, nids: list[int], comment: str='') -> int:
        """
        Creates a CMASS4 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            SPOINT mass
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        mass_obj = self.cmass4.add(eid, mass, nids, comment=comment)
        return mass_obj


class Add1dElements(BDFAttributes):
    def add_conrod(self, eid: int, mid: int, nids: list[int],
                   A: float=0.0, j: float=0.0, c: float=0.0, nsm: float=0.0,
                   comment: str='') -> int:
        """
        Creates a CONROD card

        Parameters
        ----------
        eid : int
            element id
        mid : int
            material id
        nids : list[int, int]
            node ids
        A : float; default=0.
            area
        j : float; default=0.
            polar moment of inertia
        c : float; default=0.
            stress factor
        nsm : float; default=0.
            non-structural mass per unit length
        comment : str; default=''
            a comment for the card

        """
        elem = self.conrod.add(eid, mid, nids, A=A, j=j, c=c, nsm=nsm, comment=comment)
        return elem

    def add_crod(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CROD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PROD)
        nids : list[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.crod.add(eid, pid, nids, comment=comment)
        return elem

    def add_prod(self, pid: int, mid: int, A: float,
                 j: float=0., c: float=0., nsm: float=0., comment: str='') -> int:
        """
        Creates a PROD card

        Parameters
        ----------
        pid : int
           property id
        mid : int
           material id
        A : float
           area
        J : float; default=0.
           polar moment of inertia
        c : float; default=0.
           stress factor
        nsm : float; default=0.
           nonstructural mass per unit length
        comment : str; default=''
            a comment for the card

        """
        prop = self.prod.add(pid, mid, A, j=j, c=c, nsm=nsm, comment=comment)
        return prop

    def add_ctube(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CTUBE card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        nids : list[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.ctube.add(eid, pid, nids, comment=comment)
        return elem

    def add_ptube(self, pid: int, mid: int, OD1: float, t: Optional[float]=None,
                  nsm: float=0., OD2: Optional[float]=None, comment: str='') -> int:
        """
        Adds a PTUBE card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        OD1 : float
            outer diameter at End A
        t : float; default=None -> OD1/2.
            thickness
        nsm : float; default=0.
            non-structural mass per unit length
        OD2 : float; default=None -> OD1
            outer diameter at End B
        comment : str; default=''
            a comment for the card

        """
        prop = self.ptube.add(pid, mid, OD1, t=t, nsm=nsm, OD2=OD2, comment=comment)
        return prop

    def add_baror(self, pid: int, is_g0, g0, x, offt: str='GGG', comment: str='') -> BAROR:
        baror = BAROR(pid, g0, x, offt=offt, comment=comment)
        assert self.baror is None
        self.baror = baror
        return self.baror

    def add_cbarao(self, eid: int, scale: str, x: list[float], comment: str='') -> CBARAO:
        """
        Creates a CBARAO card, which defines additional output locations
        for the CBAR card.

        It also changes the OP2 element type from a CBAR-34 to a CBAR-100.
        However, it is ignored if there are no PLOAD1s in the model.
        Furthermore, the type is changed for the whole deck, regardless of
        whether there are PLOAD1s in the other load cases.

        Parameters
        ----------
        eid : int
            element id
        scale : str
            defines what x means
            LE : x is in absolute coordinates along the bar
            FR : x is in fractional
        x : list[float]
            the additional output locations
            len(x) <= 6
        comment : str; default=''
            a comment for the card

        Notes
        -----
        MSC only

        """
        elem_flag = self.cbarao.add(eid, scale, x, comment=comment)
        return self.cbarao

    def add_cbar(self, eid: int, pid: int, nids: list[int],
                 x: Optional[list[float]], g0: Optional[int],
                 offt: str='GGG', pa: int=0, pb: int=0,
                 wa: Optional[list[float]]=None, wb: Optional[list[float]]=None,
                 comment: str='', validate: bool=False) -> int:
        """
        Adds a CBAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        x : list[float, float, float]
            Components of orientation vector, from GA, in the displacement
            coordinate system at GA (default), or in the basic coordinate system
        g0 : int
            Alternate method to supply the orientation vector using grid
            point G0. Direction of is from GA to G0. is then transferred
            to End A
        offt : str; default='GGG'
            Offset vector interpretation flag
        pa / pb : int; default=0
            Pin Flag at End A/B.  Releases the specified DOFs
        wa / wb : list[float, float, float]
            Components of offset vectors from the grid points to the end
            points of the axis of the shear center
        comment : str; default=''
            a comment for the card

        """
        if validate:
            for nid in nids:
                assert nid in self.nodes, f'nid={nid!r} does not exist'
        elem = self.cbar.add(eid, pid, nids, x, g0, offt=offt, pa=pa, pb=pb,
                             wa=wa, wb=wb, comment=comment)
        return elem

    #def add_pbarl_dvprel1(self, pid: int, mid: int,
                          #Type: str, dim: list[float], dim_constraints: list[Any],
                          #group: str='MSCBML0', nsm: float=0.,
                          #comment: str='') -> tuple[PBARL, list[DESVAR], list[DVPREL1]]:
        #"""
        #dim = [0.1, 0.2, 0.3, 0.4]
        #dim_constraints = [
            #None,
            #[0.01, 1.0],
            #[None, 1.0],
            #None,
        #]"""
        #assert len(dim) == len(dim_constraints), f'len(dim)={len(dim)} len(dim_constraints)={len(dim_constraints)}'
        #pbarl = self.add_pbarl(
            #pid, mid, Type, dim, group=group, nsm=nsm,
            #comment: str='')
        #prop_type = 'PBAR'
        #desvar_id = max(self.desvars) + 1
        #oid = max(self.dvprels) + 1
        #desvars = []
        #dvprels = []
        #for i, dim, dim_min_max in zip(count(), dim, dim_constraints):
            #if dim_min_max is None:
                #continue
            #xinit = dim
            #dim_min, dim_max = dim_min_max
            #pname_fid = f'DIM{i+1:d}'
            #label = pname_fid
            #dvids = [desvar_id]
            #coeffs = [1.]
            #xlb = -1e20 if dim_min is None else dim_min
            #xub = 1e20 if dim_max is None else dim_max
            #desvar = self.add_desvar(
                #desvar_id, label, xinit,
                #xlb=xlb, xub=xub,
                #delx=None, ddval=None,
                #comment: str='')
            #dvprel1 = self.add_dvprel1(
                #oid, prop_type, pid, pname_fid, dvids, coeffs,
                #p_min=None, p_max=1e20, c0=0.0,
                #validate=True, comment: str='')
            #desvars.append(desvar)
            #dvprels.append(dvprel1)
            #desvar_id += 1
            #oid += 1
        #return pbarl, desvars, dvprels

    def add_pbar(self, pid: int, mid: int, A: float=0.,
                 i1: float=0., i2: float=0., i12: float=0., j: float=0.,
                 nsm: float=0.,
                 c1: float=0., c2: float=0.,
                 d1: float=0., d2: float=0.,
                 e1: float=0., e2: float=0.,
                 f1: float=0., f2: float=0.,
                 k1: float=1.e8, k2: float=1.e8, comment: str='') -> int:
        """
        Creates a PBAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        area : float
            area
        i1, i2, i12, j : float
            moments of inertia
        nsm : float; default=0.
            nonstructural mass per unit length
        c1/c2, d1/d2, e1/e2, f1/f2 : float
           the y/z locations of the stress recovery points
           c1 - point C.y
           c2 - point C.z

        k1 / k2 : float; default=1.e8
            Shear stiffness factor K in K*A*G for plane 1/2.
        comment : str; default=''
            a comment for the card

        """
        prop = self.pbar.add(pid, mid, A=A, i1=i1, i2=i2, i12=i12, j=j, nsm=nsm,
                             c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2,
                             f1=f1, f2=f2, k1=k1, k2=k2, comment=comment)
        return prop

    def add_pbarl(self, pid: int, mid: int, bar_type: str, dim: list[float],
                  group: str='MSCBML0', nsm: float=0., comment: str='') -> int:
        """
        Creates a PBARL card, which defines A, I1, I2, I12, and J using
        dimensions rather than explicit values.

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        bar_type : str
            type of the bar
            valid_types = {
                ROD, TUBE, TUBE2,
                I, CHAN, T, BOX, BAR, CROSS, H, T1,
                I1, CHAN1, Z, CHAN2, T2, BOX1, HEXA, HAT, HAT1, DBOX
            }
        dim : list[float]
            dimensions for cross-section corresponding to Type;
            the length varies
        group : str default='MSCBML0'
            this parameter can lead to a very broken deck with a very
            bad error message; don't touch it!
        nsm : float; default=0.
           non-structural mass
        comment : str; default=''
            a comment for the card

        The shear center and neutral axis do not coincide when:
           - Type = I and dim2 != dim3
           - Type = CHAN, CHAN1, CHAN2
           - Type = T
           - Type = T1, T2
           - Type = BOX1
           - Type = HAT, HAT1
           - Type = DBOX

        """
        prop = self.pbarl.add(pid, mid, bar_type, dim,
                              group=group, nsm=nsm, comment=comment)
        return prop

    def add_cbeam(self, eid: int, pid: int, nids: list[int],
                  x: Optional[list[float]], g0: Optional[int],
                  offt: str='GGG', bit=None,
                  pa: int=0, pb: int=0,
                  wa=None, wb=None,
                  sa: int=0, sb: int=0, comment: str='') -> int:
        """
        Adds a CBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        x : list[float, float, float]
            Components of orientation vector, from GA, in the displacement
            coordinate system at GA (default), or in the basic coordinate system
        g0 : int
            Alternate method to supply the orientation vector using grid
            point G0. Direction of is from GA to G0. is then transferred
            to End A
        offt : str; default='GGG'
            Offset vector interpretation flag
            None : bit is active
        bit : float; default=None
            Built-in twist of the cross-sectional axes about the beam axis
            at end B relative to end A.
            For beam p-elements ONLY!
            None : offt is active
        pa / pb : int; default=0
            Pin Flag at End A/B.  Releases the specified DOFs
        wa / wb : list[float, float, float]
            Components of offset vectors from the grid points to the end
            points of the axis of the shear center
        sa / sb : int; default=0
            Scalar or grid point identification numbers for the ends A and B,
            respectively. The degrees-of-freedom at these points are the
            warping variables . SA and SB cannot be specified for
            beam p-elements
        comment : str; default=''
            a comment for the card

        Notes
        -----
        offt/bit are MSC specific fields

        """
        elem = self.cbeam.add(eid, pid, nids, x, g0, offt=offt, bit=bit,
                              pa=pa, pb=pb, wa=wa, wb=wb, sa=sa, sb=sb, comment=comment)
        return elem

    def add_pbeam(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                  c1=None, c2=None, d1=None, d2=None,
                  e1=None, e2=None, f1=None, f2=None,
                  k1=1., k2=1., s1=0., s2=0.,
                  nsia=0., nsib=None, cwa=0., cwb=None,
                  m1a=0., m2a=0., m1b=None, m2b=None,
                  n1a=0., n2a=0., n1b=None, n2b=None,
                  comment: str='') -> int:
        """
        .. todo:: fix 0th entry of self.so, self.xxb

        Creates a PBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        xxb : list[float]
            The percentage locations along the beam [0., ..., 1.]
        so : list[str]
            YES, YESA, NO
        area : list[float]
            area
        i1, i2, i12, j : list[float]
            moments of inertia
        nsm : list[float]; default=None -> [0.]*nxxb
            nonstructural mass per unit length
        c1/c2, d1/d2, e1/e2, f1/f2 : list[float]; default=None -> [0.]*nxxb
           the y/z locations of the stress recovery points
           c1 - point C.y
           c2 - point C.z
        k1 / k2 : float; default=1.
            Shear stiffness factor K in K*A*G for plane 1/2.
        s1 / s2 : float; default=0.
            Shear relief coefficient due to taper for plane 1/2.
        nsia / nsia : float; default=0. / nsia
            non structural mass moment of inertia per unit length
            about nsm center of gravity at Point A/B.
        cwa / cwb : float; default=0. / cwa
            warping coefficient for end A/B.
        m1a / m2a : float; default=0. / 0.
            y/z coordinate of center of gravity of
            nonstructural mass for end A.
        m1b / m2b : float; default=m1a / m2a
            y/z coordinate of center of gravity of
            nonstructural mass for end B.
        n1a / n2a : float; default=0. / 0.
            y/z coordinate of neutral axis for end A.
        n1b / n2b : float; default=n1a / n2a
            y/z coordinate of neutral axis for end B.
        comment : str; default=''
            a comment for the card

        """
        prop = self.pbeam.add(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                              c1, c2, d1, d2, e1, e2, f1, f2,
                              k1=k1, k2=k2, s1=s1, s2=s2,
                              nsia=nsia, nsib=nsib, cwa=cwa, cwb=cwb,
                              m1a=m1a, m2a=m2a, m1b=m1b, m2b=m2b,
                              n1a=n1a, n2a=n2a, n1b=n1b, n2b=n2b, comment=comment)
        return prop

    def add_pbcomp(self, pid, mid, y, z, c, mids,
                   area=0.0, i1=0.0, i2=0.0, i12=0.0, j=0.0, nsm=0.0,
                   k1=1.0, k2=1.0, m1=0.0, m2=0.0, n1=0.0, n2=0.0,
                   symopt=0, comment: str='') -> int:
        """
        Creates a PBCOMP card

        Parameters
        ----------
        pid : int
            Property ID
        mid : int
            Material ID
        mids : list[int]
            Material ID for the i-th integration point
        y / z : list[float]
            The (y,z) coordinates of the lumped areas in the element
            coordinate system
        c : list[float]; default=0.0
            Fraction of the total area for the i-th lumped area
            default not supported...
        area : float
            Area of beam cross section
        i1 / i2 : float; default=0.0
            Area moment of inertia about plane 1/2 about the neutral axis
        i12 : float; default=0.0
           area product of inertia
        j : float; default=0.0
            Torsional moment of interia
        nsm : float; default=0.0
            Nonstructural mass per unit length
        k1 / k2 : float; default=1.0
            Shear stiffness factor K in K*A*G for plane 1/2
        m1 / m2 : float; default=0.0
            The (y,z) coordinates of center of gravity of nonstructural mass
        n1 / n2 : float; default=0.0
            The (y,z) coordinates of neutral axis
        symopt : int; default=0
            Symmetry option to input lumped areas for the beam cross section
            0 < Integer < 5
        comment : str; default=''
            a comment for the card

        """
        prop = self.pbcomp.add(pid, mid, y, z, c, mids,
                               area, i1, i2, i12, j, nsm,
                               k1, k2, m1, m2, n1, n2,
                               symopt, comment=comment)
        return prop

    def add_pbmsect(self, pid, mid, form, options, comment: str='') -> int:
        """Creates a PBMSECT card"""
        prop = self.pbmsect.add(pid, mid, form, options, comment=comment)
        return prop

    def add_pbrsect(self, pid, mid, form, options, comment: str='') -> int:
        """Creates a PBRSECT card"""
        prop = self.pbrsect.add(pid, mid, form, options, comment=comment)
        return prop

    def add_pbeam3(self, pid, mid, A, iz, iy, iyz=0., j=None, nsm=0.,
                   cy=0., cz=0., dy=0., dz=0., ey=0., ez=0., fy=0., fz=0.,
                   comment: str='') -> int:
        """Creates a PBEAM3 card"""
        prop = self.pbeam3.add(pid, mid, A, iz, iy, iyz=iyz, j=j, nsm=nsm,
                               cy=cy, cz=cz, dy=dy, dz=dz, ey=ey, ez=ez, fy=fy, fz=fz,
                               comment=comment)
        return prop

    def add_pbeaml(self, pid: int, mid: int, beam_type: str,
                   xxb, dims, so=None, nsm=None,
                   group: str='MSCBML0', comment: str='') -> int:
        """
        Creates a PBEAML card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        beam_type : str
            the section profile
        xxb : list[float]
            The percentage locations along the beam [0., ..., 1.]
        dims : list[dim]
            dim : list[float]
                The dimensions for each section
        so : list[str]; default=None
            YES, YESA, NO
            None : [0.] * len(xxb)
        nsm : list[float]; default=None
            nonstructural mass per unit length
            None : [0.] * len(xxb)
        group : str; default='MSCBML0'
            this parameter can lead to a very broken deck with a very
            bad error message; don't touch it!
        comment : str; default=''
            a comment for the card

        """
        prop = self.pbeaml.add(pid, mid, beam_type, xxb, dims,
                               group=group, so=so, nsm=nsm, comment=comment)
        return prop

    def add_pbeaml_dvprel1(self, pid: int, mid: int, beam_type: str,
                           xxb, dims, dim_constraints,
                           so=None, nsm=None,
                           #static_stress_constraints=None,
                           #static_strain_constraints=None,
                           #static_force_constraints=None,
                           group: str='MSCBML0',
                           comment: str='') -> tuple[int, list[int], list[int]]:
        """
        Returns
        -------
        PBARL id : int
        DESVAR ids : list[int]
        DVPREL1 ids : list[int]

        dim = [0.1, 0.2, 0.3, 0.4]
        dim_constraints = [
            None,
            [0.01, 1.0],
            [None, 1.0],
            None,
        ]"""
        dim_station0 = dims[0]
        assert len(dim_station0) == len(dim_constraints), f'len(dim_station0)={len(dim_station0)} len(dim_constraints)={len(dim_constraints)}'
        pbeaml = self.add_pbeaml(pid, mid, beam_type, xxb, dims,
                                 group=group, so=so, nsm=nsm, comment=comment)

        prop_type = 'PBEAML'
        desvar_id = 1 if len(self.desvars) == 0 else max(self.desvars) + 1
        oid = 1 if len(self.dvprels) == 0 else max(self.dvprels) + 1
        #dresp_id = 1 if len(self.dresps) == 0 else max(self.dresps) + 1
        #dconstr_id = 1 # if len(self.dconstrs) == 0 else max(self.dconstrs) + 1
        desvars = []
        dvprels = []
        comment = f'PBEAML-pid={pid:d}'

        desvar_id0 = desvar_id
        for i, dim, dim_min_max in zip(count(), dim_station0, dim_constraints):
            if dim_min_max is None:
                continue
            pname_fid = f'DIM{i+1:d}(A)'  # DIM1(A)
            desvar_label = f'DIM{i+1:d}_{pid}'
            if len(desvar_label) > 8:
                desvar_label = f'D{i+1:d}_{pid}'
            if len(desvar_label) > 8:
                desvar_label = f'D{i+1:d}{pid}'
            assert len(desvar_label) <= 8, desvar_label
            #label = pname_fid

            if beam_type == 'TUBE':  # t
                # desvar defines Ri, t
                # dvprel defines Ro, Ri
                desvar_Ri = desvar_id0
                desvar_t = desvar_id0 + 1
                Ro_init = dim_station0[0]
                Ri_init = dim_station0[1]
                Ro_min, Ro_max = dim_constraints[0]
                Ri_min, Ri_max = dim_constraints[1]
                t_init = Ro_init - Ri_init
                t_min = Ro_min - Ri_min
                t_max = Ro_max - Ri_max
                assert t_max > t_min, f'pid={pid} t_min={t_min} t_max={t_max}'
                if i == 0:
                    Ri_min_max = [Ri_min, Ri_max]
                    dim_min_max = Ri_min_max
                    xinit = Ri_init
                else:
                    t_min_max = [t_min, t_max]
                    dim_min_max = t_min_max
                    xinit = t_init
            else:
                xinit = dim

            dim_min, dim_max = dim_min_max
            xlb = -1e20 if dim_min is None else dim_min
            xub = 1e20 if dim_max is None else dim_max

            desvar = self.add_desvar(
                desvar_id, desvar_label, xinit,
                xlb=xlb, xub=xub,
                delx=None, ddval=None,
                comment=comment)

            if beam_type == 'TUBE':  # Ri
                # desvar defines Ri, t
                # dvprel defines Ro, Ri
                dvids = [desvar_Ri, desvar_t]  # Ro
                #Ri = 0*Ro + 1*Ri
                if i == 0:
                    coeffs = [1., 1.]  # Ro
                else:
                    coeffs = [1., 0.]
            else:
                dvids = [desvar_id]
                coeffs = [1.]

            dvprel1 = self.dvprel1.add(
                oid, prop_type, pid, pname_fid, dvids, coeffs,
                p_min=None, p_max=1e20, c0=0.0,
                validate=True, comment=comment)
            desvars.append(desvar)
            dvprels.append(dvprel1)
            desvar_id += 1
            oid += 1
            comment = ''
        #dconstrs = []
        #if static_stress_constraints:
            #label = f'o_resp{pid:d}'
            #response_type = 'STRESS'
            #dresp_id = add_beam_stress_strain_constraints(self, pid, label, response_type,
                                                          #static_stress_constraints,
                                                          #dresp_id, dconstr_id,
                                                          #dconstrs)
        #if static_strain_constraints:
            #label = f'o_resp{pid:d}'
            #response_type = 'STRAIN'
            #dresp_id = add_beam_stress_strain_constraints(self, pid, label, response_type,
                                                          #static_stress_constraints, dresp_id,
                                                          #dconstrs)
        return pbeaml, desvars, dvprels

    def add_cbend(self, eid: int, pid: int, nids,
                  g0: Optional[int],
                  x: Optional[list[float]],
                  geom: int, comment: str='') -> int:
        """Creates a CBEND card"""
        elem = self.cbend.add(eid, pid, nids, g0, x, geom, comment=comment)
        return elem

    def add_pbend(self, pid: int, mid: int, beam_type, A, i1, i2, j,
                  c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                  nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b,
                  comment: str='') -> int:
        """Creates a PBEND card"""
        prop = self.pbend.add(pid, mid, beam_type, A, i1, i2, j,
                     c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b,
                     comment=comment)
        return prop

    def add_cbeam3(self, eid, pid, nids, x=None, g0=None,
                   wa=None, wb=None, wc=None, tw=None,
                   sa: int=0, sb: int=0, sc: int=0,
                   comment: str='') -> int:
        """Creates a CBEAM3 card"""
        assert isinstance(sa, integer_types), sa
        assert isinstance(sb, integer_types), sb
        assert isinstance(sc, integer_types), sc
        assert x is not None or g0 is not None, (x, g0)
        elem = self.cbeam3.add(eid, pid, nids, x=x, g0=g0,
                               wa=wa, wb=wb, wc=wc, tw=tw,
                               sa=sa, sb=sb, sc=sc,
                               comment=comment)
        return elem


class Add2dElements(BDFAttributes):
    def add_cshear(self, eid: int, pid: int, nids: list[int],
                   comment: str='') -> int:
        """
        Creates a CSHEAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHEAR)
        nids : list[int, int, int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        elem = self.cshear.add(eid, pid, nids, comment=comment)
        return elem

    def add_pshear(self, pid: int, mid: int, t: float, nsm: float=0.,
                   f1: float=0., f2: float=0., comment: str='') -> int:
        """
        Creates a PSHEAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        t : float
            shear panel thickness
        nsm : float; default=0.
            nonstructural mass per unit length
        f1 : float; default=0.0
            Effectiveness factor for extensional stiffness along edges 1-2 and 3-4
        f2 : float; default=0.0
            Effectiveness factor for extensional stiffness along edges 2-3 and 1-4
        comment : str; default=''
            a comment for the card

        """
        prop = self.pshear.add(pid, mid, t, nsm=nsm, f1=f1, f2=f2, comment=comment)
        return prop

    def add_ctria3(self, eid: int, pid: int, nids: list[int],
                   zoffset: float=0., theta_mcid: float=0.0,
                   tflag: int=0, T1=None, T2=None, T3=None,
                   comment: str='') -> int:
        """
        Creates a CTRIA3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        elem = self.ctria3.add(
            eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
            tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        return elem

    def add_cquad4(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: int | float=0.0, zoffset: float=None,
                   tflag: int=0, T1=None, T2=None, T3=None, T4=None,
                   comment: str='') -> int:
        """
        Creates a CQUAD4 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        elem = self.cquad4.add(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
            tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        return elem

    def add_ctria6(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: float=0., zoffset: float=0.,
                   tflag: int=0, T1=None, T2=None, T3=None, comment: str='') -> int:
        """
        Creates a CTRIA6 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        elem = self.ctria6.add(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
            tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        return elem

    def add_cquad8(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: int | float=0., zoffset: float=0.,
                   tflag: int=0, T1=None, T2=None, T3=None, T4=None, comment: str='') -> int:
        """
        Creates a CQUAD8 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int, int/None, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        elem = self.cquad8.add(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
            tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        return elem

    def add_cquad(self, eid: int, pid: int, nids: list[int],
                  theta_mcid: int | float=0., comment: str='') -> int:
        """
        Creates a CQUAD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int, int/None, int/None,
                    int/None, int/None, int/None]
            node ids
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        comment : str; default=''
            a comment for the card

        """
        elem = self.cquad.add(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)
        return elem

    def add_ctriar(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: int | float=0.0, zoffset: float=0.0,
                   tflag: int=0, T1=None, T2=None, T3=None, comment: str='') -> int:
        """
        Creates a CTRIAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        elem = self.ctriar.add(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
            tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        return elem

    def add_cquadr(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: int | float=0.0, zoffset: float=0., tflag: int=0,
                   T1=None, T2=None, T3=None, T4=None, comment: str='') -> int:
        """
        Creates a CQUADR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        elem = self.cquadr.add(
            eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
            tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        return elem

    def add_snorm(self, nid: int, normal: list[float], cid: int=0, comment: str='') -> int:
        snorm = self.snorm.add(nid, normal, cid=cid, comment=comment)
        return snorm

    def add_pshell(self, pid: int, mid1: int=None, t: float=None,
                   mid2: int=None, twelveIt3: float=1.0,
                   mid3: int=None, tst: float=0.833333, nsm: float=0.0,
                   z1: float=None, z2: float=None, mid4: int=None,
                   comment: str='') -> int:
        """
        Creates a PSHELL card

        Parameters
        ----------
        pid : int
            property id
        mid1 : int; default=None
            defines membrane material
            defines element density (unless blank)
        mid2 : int; default=None
            defines bending material
            defines element density if mid1=None
        mid3 : int; default=None
            defines transverse shear material
        mid4 : int; default=None
            defines membrane-bending coupling material
        twelveIt3 : float; default=1.0
            Bending moment of inertia ratio, 12I/T^3. Ratio of the actual
            bending moment inertia of the shell, I, to the bending
            moment of inertia of a homogeneous shell, T^3/12. The default
            value is for a homogeneous shell.
        nsm : float; default=0.0
            non-structural mass per unit area
        z1 / z2 : float; default=None
            fiber distance location 1/2 for stress/strain calculations
            z1 default : -t/2 if thickness is defined
            z2 default : t/2 if thickness is defined
        comment : str; default=''
            a comment for the card

        """
        prop = self.pshell.add(pid, mid1=mid1, t=t, mid2=mid2, twelveIt3=twelveIt3,
                               mid3=mid3, tst=tst, nsm=nsm,
                               z1=z1, z2=z2, mid4=mid4,
                               comment=comment)
        return prop

    def add_pcomp(self, pid, mids, thicknesses, thetas=None, souts=None,
                  nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None,
                  z0=None, comment: str='') -> int:
        """
        Creates a PCOMP card

        Parameters
        ----------
        pid : int
            property id
        mids : list[int, ..., int]
            material ids for each ply
        thicknesses : list[float, ..., float]
            thicknesses for each ply
        thetas : list[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : list[str, ..., str]; default=None
            should the stress? be printed; {YES, NO}
            None : [NO] * nplies
        nsm : float; default=0.
            nonstructural mass per unit area
        sb : float; default=0.
            Allowable shear stress of the bonding material.
            Used by the failure theory
        ft : str; default=None
            failure theory; {HILL, HOFF, TSAI, STRN, None}
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            structural damping
        lam : str; default=None
            symmetric flag; {SYM, MEM, BEND, SMEAR, SMCORE, None}
            None : not symmmetric
        z0 : float; default=None
            Distance from the reference plane to the bottom surface
            None : -1/2 * total_thickness
        comment : str; default=''
            a comment for the card

        """
        prop = self.pcomp.add(pid, mids, thicknesses, thetas, souts,
                              nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam,
                              z0=z0, comment=comment)
        return prop

    def add_pcompg(self, pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                   nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0, lam=None, z0=None,
                   comment: str='') -> int:
        """
        Creates a PCOMPG card

        Parameters
        ----------
        pid : int
            property id
        global_ply_ids : list[int]
            the ply id
        mids : list[int, ..., int]
            material ids for each ply
        thicknesses : list[float, ..., float]
            thicknesses for each ply
        thetas : list[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : list[str, ..., str]; default=None
            should the stress? be printed; {YES, NO}
            None : [NO] * nplies
        nsm : float; default=0.
            nonstructural mass per unit area
        sb : float; default=0.
            Allowable shear stress of the bonding material.
            Used by the failure theory
        ft : str; default=None
            failure theory; {HILL, HOFF, TSAI, STRN, None}
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            structural damping
        lam : str; default=None
            symmetric flag; {SYM, MEM, BEND, SMEAR, SMCORE, None}
            None : not symmmetric
        z0 : float; default=None
            Distance from the reference plane to the bottom surface
            None : -1/2 * total_thickness
        comment : str; default=''
            a comment for the card

        """
        prop = self.pcompg.add(
            pid, global_ply_ids, mids, thicknesses, thetas=thetas, souts=souts,
            nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam, z0=z0,
            comment=comment)
        return prop


class Add3dElements(BDFAttributes):
    def add_ctetra(self, eid: int, pid: int, nids: list[int],
                   comment: str='') -> int:
        """
        Creates a CTETRA4/CTETRA10

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : list[int]
            node ids; n=4 or 10
        comment : str; default=''
            a comment for the card

        """
        elem = self.ctetra.add(eid, pid, nids, comment=comment)
        return elem

    def add_cpyram(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CPYRAM5/CPYRAM13

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : list[int]
            node ids; n=5 or 13
        comment : str; default=''
            a comment for the card

        """
        elem = self.cpyram.add(eid, pid, nids, comment=comment)
        return elem

    def add_cpenta(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CPENTA6/CPENTA15

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : list[int]
            node ids; n=6 or 15
        comment : str; default=''
            a comment for the card

        """
        elem = self.cpenta.add(eid, pid, nids, comment=comment)
        return elem

    def add_chexa(self, eid: int, pid: int,
                  nids: list[int], comment: str='') -> int:
        """
        Creates a CHEXA8/CHEXA20

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : list[int]
            node ids; n=8 or 20
        comment : str; default=''
            a comment for the card

        """
        elem = self.chexa.add(eid, pid, nids, comment=comment)
        return elem

    def add_psolid(self, pid: int, mid: int, cordm: int=0,
                   integ: int=None, stress: str=None, isop=None,
                   fctn: str='SMECH', comment: str='') -> int:
        """
        Creates a PSOLID card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        cordm : int; default=0
            material coordinate system
        integ : int; default=None
            None-varies depending on element type
            0, 'BUBBLE'
            1, 'GAUSS'
            2, 'TWO'
            3, 'THREE'
            REDUCED
            FULL
        stress : int/str; default=None
            None/GRID, 1-GAUSS
        isop : int/str; default=None
            0-REDUCED
            1-FULL
        fctn : str; default='SMECH'
            PFLUID/SMECH
        comment : str; default=''
            a comment for the card

        """
        prop = self.psolid.add(
            pid, mid, cordm=cordm, integ=integ, stress=stress, isop=isop,
            fctn=fctn, comment=comment)
        return prop

    def add_plsolid(self, pid: int, mid: int, stress_strain: str='GRID',
                    ge: float=0., comment: str='') -> int:
        """
        Creates a PLSOLID card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        stress_strain : str
            Location of stress and strain output
            valid types = {GRID, GAUSS}
        ge : float; default=0.
            damping coefficient
        comment : str; default=''
            a comment for the card

        """
        prop = self.plsolid.add(pid, mid, stress_strain=stress_strain,
                                ge=ge, comment=comment)
        return prop

    def add_pcomps(self, pid, global_ply_ids, mids, thicknesses, thetas,
                   cordm=0, psdir=13, sb=None, nb=None, tref=0.0, ge=0.0,
                   failure_theories=None, interlaminar_failure_theories=None,
                   souts=None, comment: str='') -> int:
        """Creates a PCOMPS card"""
        prop = self.pcomps.add(
            pid, global_ply_ids, mids, thicknesses, thetas,
            cordm, psdir, sb, nb, tref, ge,
            failure_theories, interlaminar_failure_theories, souts,
            comment=comment)
        return prop


class AddRigidElements(BDFAttributes):
    def add_rrod(self, eid: int, nids: list[int],
                 cma: str='', cmb: str='',
                 alpha: float=0.0, comment: str='') -> int:
        """
        Creates a RROD element

        Parameters
        ----------
        eid : int
            element id
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        cma / cmb : str; default=''
            dependent DOFs
        alpha : float; default=0.0
            coefficient of thermal expansion
        comment : str; default=''
            a comment for the card

        """
        elem = self.rrod.add(eid, nids, cma=cma, cmb=cmb, alpha=alpha, comment=comment)
        return elem

    def add_rbe1(self, eid: int,
                 Gni: list[int], Cni: list[int],
                 Gmi: list[int], Cmi: list[int],
                 alpha: float=0., comment: str='') -> int:
        """
        Creates an RBE1 element

        Parameters
        ----------
        eid : int
            element id
        Gni : list[int]
            independent node ids
        Cni : list[str]
            the independent components (e.g., '123456')
        Gmi : list[int]
            dependent node ids
        Cmi : list[str]
            the dependent components (e.g., '123456')
        alpha : float; default=0.
            thermal expansion coefficient
        comment : str; default=''
            a comment for the card

        """
        elem = self.rbe1.add(eid, Gni, Cni, Gmi, Cmi, alpha=alpha, comment=comment)
        return elem

    def add_rbe2(self, eid: int,
                 gn: int, # independent
                 cm: str, Gmi: list[int],  # dependent
                 alpha: float=0.0, tref: float=0.0, comment: str='',
                 validate: bool=False) -> int:
        """
        Creates an RBE2 element

        Parameters
        ----------
        eid : int
            element id
        gn : int
           Identification number of grid point to which all six independent
           degrees-of-freedom for the element are assigned.
        cm : str
            Component numbers of the dependent degrees-of-freedom in the
            global coordinate system at grid points GMi.
        Gmi : list[int]
            dependent nodes
        alpha : float; default=0.0
            thermal expansion coefficient
        tref : float; default=0.0
            reference temperature

        """
        elem = self.rbe2.add(eid, gn, cm, Gmi, alpha=alpha, tref=tref, comment=comment)
        if validate:
            assert gn in self.nodes, f'gn={gn!r} does not exist'
            for nid in elem.Gmi:
                assert nid in self.nodes, f'Gm={nid!r} does not exist'
        return elem

    def add_rbe3(self, eid: int, refgrid: int, refc: str,
                 weights: list[float], comps: list[str], Gijs: list[list[int]],
                 Gmi=None, Cmi=None,
                 alpha: float=0.0, tref: float=0.0,
                 validate: bool=True,
                 comment: str='') -> int:
        """
        Creates an RBE3 element

        Parameters
        ----------
        eid : int
            element id
        refgrid : int
            dependent node
        refc - str
            dependent components for refgrid???
        weights : list[float]
            independent weights for the importance of the DOF
        comps : list[str]
            independent components
            len(comps) = len(weights)
        GiJs : varies
            independent nodes
            list[list[int]]:
                allows for different nodes for the different weights
                len(GiJs) = len(weights)
            list[int, ..., int]:
                intended for a single weight
                This will be expanded into list[list[int]]
        Gmi : list[int, ..., int]; default=None -> []
            dependent nodes / UM Set
        Cmi : list[str, ..., str]; default=None -> []
            dependent components / UM Set
        alpha : float; default=0.0
            thermal expansion coefficient
        tref : float; default=0.0
            reference temperature
        comment : str; default=''
            a comment for the card

        """
        #weights: list[float], comps: list[str], Gijs: list[int],
        if isinstance(Gijs[0], integer_types):
            Gijs2 = []
            for Gij in Gijs:
                assert isinstance(Gij, integer_types), 'Gij=%s type=%s' % (Gij, type(Gij))
                Gijs2.append([Gij])
            Gijs = Gijs2

        if validate:
            assert len(weights) == len(comps)
            assert len(weights) == len(Gijs)
            for Gij in Gijs:
                try:
                    len(Gij)
                except TypeError:
                    msg = f'weights={weights} comps={comps} Gijs={Gijs}'
                    raise TypeError(f'RBE3 eid={eid} ' + msg)
            #independent_nodes.extend(Gij)
            #ngrid_per_weight.append(len(Gij))

        elem = self.rbe3.add(eid, refgrid, refc, weights, comps, Gijs,
                             Gmi=Gmi, Cmi=Cmi, alpha=alpha, tref=tref, comment=comment)
        return elem

    def add_rbar(self, eid: int, nids: list[int],
                 cna: str, cnb: str,
                 cma: str, cmb: str,
                 alpha: float=0.,
                 tref: float=0.,
                 comment: str='') -> int:
        """
        Creates a RBAR element

        Parameters
        ----------
        eid : int
            element id
        nids : list[int, int]
            node ids; connected grid points at ends A and B
        cna / cnb : str
            independent DOFs in '123456'
        cma / cmb : str
            dependent DOFs in '123456'
        alpha : float; default=0.0
            coefficient of thermal expansion
        tref : float; default=0.0
            reference temperature
        comment : str; default=''
            a comment for the card

        """
        elem = self.rbar.add(eid, nids, cna, cnb, cma, cmb,
                             alpha=alpha, tref=tref, comment=comment)
        return elem

    def add_rbar1(self, eid: int, nids: list[int], cb: str,
                  alpha: float=0., comment: str='') -> int:
        """Creates an RBAR1 element"""
        elem = self.rbar1.add(eid, nids, cb, alpha=alpha, comment=comment)
        return elem

    def add_rspline(self, eid: int, independent_nid: int,
                    dependent_nids: list[int], dependent_components: list[int],
                    diameter_ratio: float=0.1, comment: str='') -> int:
        """
        Creates an RSPLINE card, which uses multipoint constraints for the
        interpolation of displacements at grid points

        Parameters
        ----------
        eid : int
            element id
        independent_nid : int
            the independent node id
        dependent_nids : list[int]
            the dependent node ids
        dependent_components : list[str]
            Components to be constrained
        diameter_ratio : float; default=0.1
            Ratio of the diameter of the elastic tube to the sum of the
            lengths of all segments
        comment : str; default=''
            a comment for the card

        """
        elem = self.rspline.add(eid, independent_nid, dependent_nids, dependent_components,
                                diameter_ratio=diameter_ratio, comment=comment)
        return elem

    def add_rsscon(self, eid: int, rigid_type: str,
                   shell_eid=None, solid_eid=None,
                   a_solid_grids=None, b_solid_grids=None, shell_grids=None,
                   comment: str='') -> int:
        """
        Creates an RSSCON card, which defines multipoint constraints to
        model clamped connections of shell-to-solid elements.

        Parameters
        ----------
        eid : int
            element id
        rigid_type : str
            GRID/ELEM
        shell/solid_eid : int; default=None
            the shell/solid element id (if rigid_type=ELEM)
        shell/solid_grids : list[int, int]; default=None
            the shell/solid node ids (if rigid_type=GRID)
        comment : str; default=''
            a comment for the card

        """
        elem = self.rsscon.add(
            eid, rigid_type,
            shell_eid=shell_eid, solid_eid=solid_eid,
            a_solid_grids=a_solid_grids, b_solid_grids=b_solid_grids,
            shell_grids=shell_grids,
            comment=comment)
        return elem


class AddAero(BDFAttributes):
    def add_aeros(self, cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0,
                  comment: str='') -> AEROS:
        """
        Creates an AEROS card

        Parameters
        ----------
        cref : float
            the aerodynamic chord
        bref : float
            the wing span
        sref : float
            the wing area
        acsid : int; default=0
            aerodyanmic coordinate system
        rcsid : int; default=0
            coordinate system for rigid body motions
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        aeros = AEROS(cref, bref, sref, acsid=acsid, rcsid=rcsid,
                      sym_xz=sym_xz, sym_xy=sym_xy, comment=comment)
        self._add_methods.add_aeros_object(aeros)
        return aeros

    def add_aero(self, velocity: float, cref: float, rho_ref: float,
                 acsid: int=0, sym_xz: int=0, sym_xy: int=0,
                 comment: str='') -> AERO:
        """
        Creates an AERO card

        Parameters
        ----------
        velocity : float
            the airspeed
        cref : float
            the aerodynamic chord
        rho_ref : float
            FLFACT density scaling factor
        acsid : int; default=0
            aerodyanmic coordinate system
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=sym_xz, sym_xy=sym_xy,
                    comment=comment)
        self._add_methods.add_aero_object(aero)
        return aero

    def add_caero1(self, eid: int, pid: int, igroup: int,
                   p1: NDArray3float, x12: float,
                   p4: NDArray3float, x43: float,
                   cp: int=0,
                   nspan: int=0, lspan: int=0,
                   nchord: int=0, lchord: int=0, comment: str='') -> int:
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        pid : int
            int : PAERO1 ID
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int; default=0
            int : coordinate system
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
        comment : str; default=''
             a comment for the card

        """
        caero = self.caero1.add(
            eid, pid, igroup, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment=comment)
        return caero

    def add_caero2(self, eid: int, pid: int, igroup: int,
                   p1: list[float], x12: float,
                   cp: int=0,
                   nsb: int=0, nint: int=0,
                   lsb: int=0, lint: int=0, comment: str='') -> int:
        """
        Defines a CAERO2 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO2
            int : PAERO2 ID
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (forward position)
        x12 : float
            length of the CAERO2
        cp : int; default=0
            int : coordinate system
        nsb : int; default=0
            Number of slender body elements
        lsb : int; default=0
            AEFACT id for defining the location of the slender body elements
        nint : int; default=0
            Number of interference elements
        lint : int; default=0
            AEFACT id for defining the location of interference elements
        comment : str; default=''
            a comment for the card

        """
        caero = self.caero2.add(eid, pid, igroup, p1, x12, cp=cp, nsb=nsb, nint=nint, lsb=lsb,
                                lint=lint, comment=comment)
        return caero

    def add_caero3(self, eid: int, pid: int,
                   p1: np.ndarray, x12: float,
                   p4: np.ndarray, x43: float,
                   cp: int=0, list_w: int=0, list_c1=None, list_c2=None, comment: str='') -> int:
        """Creates a CAERO3 card"""
        caero = self.caero3.add(
            eid, pid, p1, x12, p4, x43,
            cp=cp, list_w=list_w, list_c1=list_c1, list_c2=list_c2, comment=comment)
        return caero

    def add_caero4(self, eid: int, pid: int,
                   p1: np.ndarray, x12: float,
                   p4: np.ndarray, x43: float,
                   cp: int=0, nspan: int=0, lspan: int=0, comment: str='') -> int:
        """
        Defines a CAERO4 card, which defines a strip theory surface.

        Parameters
        ----------
        eid : int
            element id
        pid : int
            int : PAERO4 ID
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2
            (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3
            (typically x, tip chord)
        cp : int; default=0
            int : coordinate system
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
        comment : str; default=''
             a comment for the card

        """
        caero = self.caero4.add(
            eid, pid, p1, x12, p4, x43,
            cp=cp, nspan=nspan, lspan=lspan, comment=comment)
        return caero

    def add_caero5(self, eid: int, pid: int,
                   p1: list[float], x12: float,
                   p4: list[float], x43: float,
                   cp: int=0,
                   nspan: int=0, lspan: int=0,
                   ntheory: int=0,
                   nthick: int=0, comment: str='') -> int:
        """Creates a CAERO5 card"""
        caero = self.caero5.add(
            eid, pid, p1, x12, p4, x43, cp=cp, nspan=nspan, lspan=lspan,
            ntheory=ntheory, nthick=nthick, comment=comment)
        return caero

    def add_caero7(self, eid: int, label: str,
                   p1: np.ndarray, x12: float,
                   p4: np.ndarray, x43: float,
                   cp: int=0,
                   nspan: int=0, lspan: int=0,
                   nchord: int=0,
                   p_airfoil: int=0, ztaic: int=0, comment: str='') -> int:
        caero = self.caero7.add(
            eid, label, p1, x12, p4, x43,
            cp=cp, nspan=nspan,
            nchord=nchord, lspan=lspan,
            p_airfoil=p_airfoil, ztaic=ztaic, comment=comment)
        return caero

    def add_paero1(self, pid: int, caero_body_ids: Optional[list[int]]=None, comment: str='') -> int:
        """
        Creates a PAERO1 card, which defines associated bodies for the
        panels in the Doublet-Lattice method.

        Parameters
        ----------
        pid : int
            PAERO1 id
        caero_body_ids : list[int]; default=None
            CAERO2 ids that are within the same IGID group
        comment : str; default=''
            a comment for the card

        """
        paero = self.paero1.add(pid, caero_body_ids=caero_body_ids, comment=comment)
        return paero

    def add_paero2(self, pid: int, orient: str, width: float, AR: float,
                   thi: list[int], thn: list[int],
                   lrsb: Optional[int]=None,
                   lrib: Optional[int]=None,
                   lth: Optional[int]=None,
                   comment: str='') -> int:
        """
        Creates a PAERO2 card, which defines additional cross-sectional
        properties for the CAERO2 geometry.

        Parameters
        ----------
        pid : int
            PAERO1 id
        orient : str
            Orientation flag. Type of motion allowed for bodies. Refers
            to the aerodynamic coordinate system of ACSID. See AERO entry.
            valid_orientations = {Z, Y, ZY}
        width : float
            Reference half-width of body and the width of the constant
            width interference tube
        AR : float
            Aspect ratio of the interference tube (height/width)
        thi / thn : list[int]
            The first (thi) and last (thn) interference element of a body
            to use the theta1/theta2 array
        lrsb : int; default=None
            int : AEFACT id containing a list of slender body half-widths
                  at the end points of the slender body elements
            None : use width
        lrib : int; default=None
            int : AEFACT id containing a list of interference body
                  half-widths at the end points of the interference elements
            None : use width
        lth : list[int, int]; default=None
            AEFACT ids for defining theta arrays for interference calculations
            for theta1/theta2; length=2
        comment : str; default=''
            a comment for the card

        """
        paero = self.paero2.add(
            pid, orient, width, AR, thi, thn, lrsb=lrsb, lrib=lrib,
            lth=lth, comment=comment)
        return paero

    def add_paero3(self, pid: int, nbox: int, ncontrol_surfaces: int,
                   x: list[float], y: list[float], comment: str='') -> int:
        """
        Creates a PAERO3 card, which defines the number of Mach boxes
        in the flow direction and the location of cranks and control
        surfaces of a Mach box lifting surface.

        Parameters
        ----------
        pid : int
            PAERO1 id
        nbox : int
            Number of Mach boxes in the flow direction; 0 < nbox < 50
        ncontrol_surfaces : int
            Number of control surfaces. (0, 1, or 2)
        x / y : list[float, None]
            float : locations of points 5 through 12, which are in the
            aerodynamic coordinate system, to define the cranks and
            control surface geometry.
        comment : str; default=''
            a comment for the card

        """
        paero = self.paero3.add(
            pid, nbox, ncontrol_surfaces, x, y, comment=comment)
        return paero

    def add_paero4(self, pid: int,
                   docs: list[float], caocs: list[float], gapocs: list[float],
                   cla: int=0, lcla: int=0,
                   circ: int=0, lcirc: int=0, comment: str='') -> int:
        """
        Parameters
        ----------
        PID : int
            Property identification number. (Integer > 0)
        CLA : int; default=0
            Select Prandtl-Glauert correction. (Integer = -1, 0, 1)
            -1 Compressibility correction made to lift curve slope data for a reference Mach number.
            0  No correction and no list needed. (Default)
            +1 No correction and lift curve slope provided by a list as a
               function of strip location and Mach number.
        LCLA : int
            ID number of the AEFACT entry that lists the lift curve slope
            on all strips for each Mach number on the MKAEROi entry. See
            Remark 2 below. (Integer = 0 if CLA = 0, > 0 if CLA ≠ 0)
        CIRC : int; default=0
            Select Theodorsen’s function C(k) or the number of exponential
            coefficients used to approximate C(k).
            (Integer = 0, 1, 2, 3; Must be zero if CLA ≠ 0.)
            0 Theodorsen function.
            1, 2, 3 Approximate function with b0, b1, β1, ..., bn, βn n = 1, 2, 3.
        LCIRC : int
            Identification number of the AEFACT entry that lists the b, β values
            for each Mach number. See Remark 3, 4, and 5 below; variable b’s
            and β’s for each mi on the MKAEROi entry.
            (Integer = 0 if CIRC = 0, > 0 if CIRC ≠ 0)
        DOCi : list[float]
            d/c = distance of the control surface hinge aft of the quarter-chord
            divided by the strip chord (Real ≥ 0.0)
        CAOCi : list[float]
            ca/c = control surface chord divided by the strip chord. (Real ≥ 0.0)
        GAPOCi : list[float]
            g/c = control surface gap divided by the strip chord. (Real ≥ 0.0)

        """
        paero = self.paero4.add(
            pid, docs, caocs, gapocs, cla=cla, lcla=lcla,
            circ=circ, lcirc=lcirc, comment=comment)
        return paero

    def add_paero5(self, pid: int, caoci: list[float],
                   nalpha: int=0, lalpha: int=0,
                   nxis: int=0, lxis: int=0,
                   ntaus: int=0, ltaus: int=0,
                   comment: str='') -> int:
        """Creates a PAERO5 card"""
        paero = self.paero5.add(
            pid, caoci, nalpha=nalpha, lalpha=lalpha, nxis=nxis, lxis=lxis,
            ntaus=ntaus, ltaus=ltaus, comment=comment)
        return paero

    def add_spline1(self, eid: int, caero: int, box1: int, box2: int, setg: int,
                    dz: float=0., method: str='IPS',
                    usage: str='BOTH', nelements: int=10,
                    melements: int=10, comment: str='') -> int:
        """
        Creates a SPLINE1, which defines a surface spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements : int; default=10
            The number of FE elements along the local spline x-axis if
            using the FPS option
        melements : int; default=10
            The number of FE elements along the local spline y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        spline = self.spline1.add(
            eid, caero, box1, box2, setg, dz=dz, method=method,
            usage=usage, nelements=nelements, melements=melements,
            comment=comment)
        return spline

    def add_spline2(self, eid: int, caero: int,
                    box1: int, box2: int, setg: int,
                    dz: float=0.0, dtor: float=1.0,
                    cid: int=0,
                    dthx: float=0.0, dthy: float=0.0,
                    usage: str='BOTH',
                    comment: str='') -> int:
        """
        Creates a SPLINE2 card, which defines a beam spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box/body id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the beam spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        dtor : float; default=1.0
            Torsional flexibility ratio (EI/GJ).
            Use 1.0 for bodies (CAERO2).
        cid : int; default=0
            Rectangular coordinate system for which the y-axis defines the
            axis of the spline. Not used for bodies, CAERO2
        dthx : float; default=None
            Rotational attachment flexibility.
            DTHX : Used for rotation about the spline's x-axis (in-plane
                   bending rotations).  It is not used for bodies (CAERO2).
            DTHY : Used for rotation about the spline's y-axis (torsion).
                   It is used for slope of bodies.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        spline = self.spline2.add(
            eid, caero, box1, box2, setg, dz=dz, dtor=dtor, cid=cid,
            dthx=dthx, dthy=dthy, usage=usage, comment=comment)
        return spline

    def add_spline3(self, eid: int, caero: int, box_id: int,
                    components: int,
                    nodes: list[int],
                    displacement_components: list[int],
                    coeffs: list[float],
                    usage: str='BOTH', comment: str='') -> int:
        """
        Creates a SPLINE3 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box_id : int
           Identification number of the aerodynamic box number.
        components : int
           The component of motion to be interpolated.
           3, 5          (CAERO1)
           2, 3, 5, 6    (CAERO2)
           3             (CAERO3)
           3, 5, 6       (CAERO4)
           3, 5, 6       (CAERO5)
           1, 2, 3, 5, 6 (3D Geometry)
           2-lateral displacement
           3-transverse displacement
           5-pitch angle
           6-relative control angle for CAERO4/5; yaw angle for CAERO2

        nodes : list[int]
           Grid point identification number of the independent grid point.
        displacement_components : list[int]
           Component numbers in the displacement coordinate system.
           1-6 (GRIDs)
           0 (SPOINTs)
        coeffs : list[float]
           Coefficient of the constraint relationship.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        spline = self.spline3.add(
            eid, caero, box_id, components, nodes,
            displacement_components, coeffs, usage=usage,
            comment=comment)
        return spline

    def add_spline4(self, eid: int, caero: int, aelist: int, setg: int,
                    dz: float=0., method: str='IPS', usage: str='BOTH',
                    nelements: int=10, melements: int=10,
                    ftype: Optional[str]='WF2',
                    rcore: Optional[float]=None,
                    comment: str='') -> int:
        """
        Creates a SPLINE4 card, which defines a curved Infinite Plate,
        Thin Plate, or Finite Plate Spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements / melements : int; default=10
            The number of FE elements along the local spline x/y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        spline = self.spline4.add(
            eid, caero, aelist, setg, dz, method, usage,
            nelements, melements, ftype=ftype, rcore=rcore, comment=comment)
        return spline

    def add_spline5(self, eid: int, caero: int, aelist: int, setg: int, thx, thy,
                    dz: float=0.0, dtor: float=1.0, cid: int=0,
                    usage: str='BOTH', method: str='BEAM',
                    ftype: str='WF2', rcore=None, comment: str='') -> int:
        """Creates a SPLINE5 card"""
        assert isinstance(cid, int), cid
        spline = self.spline5.add(
            eid, caero, aelist, setg, thx, thy, dz=dz, dtor=dtor, cid=cid,
            usage=usage, method=method, ftype=ftype, rcore=rcore, comment=comment)
        return spline

    def add_trim(self, sid: int, mach: float, q: float,
                 labels: list[str], uxs: list[float], aeqr: float=1.0,
                 trim_type: int=1, comment: str='') -> int:
        """
        Creates a TRIM/TRIM2 card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        mach : float
            the mach number
        q : float
            dynamic pressure
        labels : list[str]
            names of the fixed variables
        uxs : list[float]
            values corresponding to labels
        aeqr : float
            0.0 : rigid trim analysis
            1.0 : elastic trim analysis (default)
        trim_type : int
            1 : creates a TRIM
            2 : creates a TRIM2
        comment : str; default=''
            a comment for the card

        """
        if trim_type == 1:
            trim = self.trim.add(sid, mach, q, labels, uxs, aeqr=aeqr, comment=comment)
            #trim = self.trim
        elif trim_type == 2:
            trim = self.trim2.add(sid, mach, q, labels, uxs, aeqr=aeqr, comment=comment)
            #trim = self.trim2
        else:  # pragma: no cover
            raise ValueError(trim_type)
        return trim

    def add_mkaero1(self, machs: list[float], reduced_freqs: list[float],
                    comment: str='') -> MKAERO1:
        """
        Creates an MKAERO1 card, which defines a set of mach and
        reduced frequencies.

        Parameters
        ----------
        machs : list[float]
            series of Mach numbers
        reduced_freqs : list[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card

        """
        mkaero = MKAERO1(machs, reduced_freqs, comment=comment)
        self._add_methods.add_mkaero_object(mkaero)
        return mkaero

    def add_mkaero2(self, machs, reduced_freqs, comment: str='') -> MKAERO2:
        """
        Creates an MKAERO2 card, which defines a set of mach and
        reduced frequency pairs.

        Parameters
        ----------
        machs : list[float]
            series of Mach numbers
        reduced_freqs : list[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card

        """
        mkaero = MKAERO2(machs, reduced_freqs, comment=comment)
        self._add_methods.add_mkaero_object(mkaero)
        return mkaero

    def add_gust(self, sid: int, dload: int, wg: float, x0: float, V: Optional[float]=None,
                 comment: str='') -> int:
        """
        Creates a GUST card, which defines a stationary vertical gust
        for use in aeroelastic response analysis.

        Parameters
        ----------
        sid : int
            gust load id
        dload : int
            TLOADx or RLOADx entry that defines the time/frequency
            dependence
        wg : float
            Scale factor (gust velocity/forward velocity) for gust
            velocity
        x0 : float
            Streamwise location in the aerodynamic coordinate system of
            the gust reference point.
        V : float; default=None
            float : velocity of the vehicle (must be the same as the
                    velocity on the AERO card)
            None : ???
        comment : str; default=''
            a comment for the card

        """
        gust = self.gust.add(sid, dload, wg, x0, V=V, comment=comment)
        return gust

    def add_flutter(self, sid: int, method: str,
                    density: int, mach: int, reduced_freq_velocity: int,
                    imethod: str='L',
                    nvalue: Optional[int]=None,
                    omax: Optional[float]=None,
                    epsilon: float=1.0e-3,
                    comment: str='', validate: bool=False) -> int:
        """
        Creates a FLUTTER card, which is required for a flutter (SOL 145)
        analysis.

        Parameters
        ----------
        sid : int
            flutter id
        method : str
            valid methods = [K, KE,
                             PKS, PKNLS, PKNL, PKE]
        density : int
            defines a series of air densities in units of mass/volume
            PARAM,WTMASS does not affect this
            AERO affects this
            references an FLFACT id
        mach : int
            defines a series of the mach numbers
            references an FLFACT id
        reduced_freq_velocity : int
            Defines a series of either:
               1) reduced frequencies - K, KE
               2) velocities - PK, PKNL, PKS, PKNLS
            depending on the method chosen.
            references an FLFACT id
        imethod : str; default='L'
            Choice of interpolation method for aerodynamic matrix interpolation.
            imethods :
               1) L - linear
               2) S - surface
               3) TCUB - termwise cubic
        nvalue : int
            Number of eigenvalues beginning with the first eigenvalue for
            output and plots
        omax : float
            For the PKS and PKNLS methods, OMAX specifies the maximum
            frequency (Hz), to be used in he flutter sweep.
            MSC only.
        epsilon : float; default=1.0e-3
            Convergence parameter for k. Used in the PK and PKNL methods only
        comment : str; default=''
            a comment for the card

        """
        flutter = self.flutter.add(
            sid, method, density, mach, reduced_freq_velocity,
            imethod=imethod, nvalue=nvalue,
            omax=omax, epsilon=epsilon,
            comment=comment, validate=validate)
        return flutter

    def add_flfact(self, sid: int, factors: list[float], comment: str='') -> int:
        """
        Creates an FLFACT card, which defines factors used for flutter
        analysis.  These factors define either:
         - density
         - mach
         - velocity
         - reduced frequency
        depending on the FLUTTER method chosen (e.g., PK, PKNL, PKNLS)

        Parameters
        ----------
        sid : int
            the id of a density, reduced_frequency, mach, or velocity table
            the FLUTTER card defines the meaning
        factors : varies
            values : list[float, ..., float]
                list of factors
            list[f1, THRU, fnf, nf, fmid]
                f1 : float
                    first value
                THRU : str
                    the word THRU
                fnf : float
                    second value
                nf : int
                    number of values
                fmid : float; default=(f1 + fnf) / 2.
                    the mid point to bias the array
                TODO: does f1 need be be greater than f2/fnf???
        comment : str; default=''
            a comment for the card

        """
        flfact = self.flfact.add(sid, factors, comment=comment)
        return flfact

    def add_aecomp(self, name: str, list_type: list[str], lists: int | list[int],
                   comment: str='') -> int:
        """
        Creates an AECOMP card

        Parameters
        ----------
        name : str
            the name of the component
        list_type : str
            One of CAERO, AELIST or CMPID for aerodynamic components and
            SET1 for structural components. Aerodynamic components are
            defined on the aerodynamic ks-set mesh while the structural
            components are defined on the g-set mesh.
        lists : list[int, int, ...]; int
            The identification number of either SET1, AELIST or CAEROi
            entries that define the set of grid points that comprise
            the component
        comment : str; default=''
            a comment for the card

        """
        aecomp = self.aecomp.add(name, list_type, lists, comment=comment)
        return aecomp

    def add_aecompl(self, name: str, labels: list[str], comment: str='') -> AECOMPL:
        """
        Creates an AECOMPL card

        Parameters
        ----------
        name : str
            the name of the component
        labels : list[str, str, ...]; str
            A string of 8 characters referring to the names of other components
            defined by either AECOMP or other AECOMPL entries.
        comment : str; default=''
            a comment for the card

        """
        aecompl = AECOMPL(name, labels, comment=comment)
        self._add_methods.add_aecomp_object(aecompl)
        return aecompl

    def add_aestat(self, aestat_id: int, label: str, comment: str='') -> int:
        """
        Creates an AESTAT card, which is a variable to be used in a TRIM analysis

        Parameters
        ----------
        aestat_id : int
            unique id
        label : str
            name for the id
        comment : str; default=''
            a comment for the card

        """
        aestat = self.aestat.add(aestat_id, label, comment=comment)
        return aestat

    def add_aelink(self, aelink_id: int, label: str,
                   independent_labels: list[str], linking_coefficients: list[float],
                   comment: str='') -> int:
        """
        Creates an AELINK card, which defines an equation linking
        AESTAT and AESURF cards

        Parameters
        ----------
        aelink_id : int
            unique id
        label : str
            name of the dependent AESURF card
        independent_labels : list[str, ..., str]
            name for the independent variables (AESTATs)
        linking_coefficients : list[float]
            linking coefficients
        comment : str; default=''
            a comment for the card

        """
        aelink = self.aelink.add(aelink_id, label, independent_labels,
                                 linking_coefficients, comment=comment)
        return aelink

    def add_aelist(self, sid: int, elements: list[int], comment: str='') -> int:
        """
        Creates an AELIST card, which defines the aero boxes for
        an AESURF/SPLINEx.

        Parameters
        ----------
        sid : int
            unique id
        elements : list[int, ..., int]
            list of box ids
        comment : str; default=''
            a comment for the card

        """
        aelist = self.aelist.add(sid, elements, comment=comment)
        return aelist

    def add_aefact(self, sid: int, fractions: list[float], comment: str='') -> int:
        """
        Creates an AEFACT card, which is used by the CAEROx / PAEROx card
        to adjust the spacing of the sub-paneleing (and grid point
        paneling in the case of the CAERO3).

        Parameters
        ----------
        sid : int
            unique id
        fractions : list[float, ..., float]
            list of percentages
        comment : str; default=''
            a comment for the card

        """
        aefact = self.aefact.add(sid, fractions, comment=comment)
        return aefact

    def add_diverg(self, sid: int, nroots: int, machs: list[float], comment: str='') -> int:
        """
        Creates an DIVERG card, which is used in divergence
        analysis (SOL 144).

        Parameters
        ----------
        sid : int
            The name
        nroots : int
            the number of roots
        machs : list[float, ..., float]
            list of Mach numbers
        comment : str; default=''
            a comment for the card

        """
        diverg = self.diverg.add(sid, nroots, machs, comment=comment)
        return diverg

    def add_csschd(self, sid: int, aesurf_id: int,
                   lschd: int, lalpha: int=None, lmach: int=None,  # aefact
                   comment: str='') -> int:
        """
        Creates an CSSCHD card, which defines a specified control surface
        deflection as a function of Mach and alpha (used in SOL 144/146).

        Parameters
        ----------
        sid : int
            the unique id
        aesurf_id/aesid : int
            the control surface (AESURF) id
        lalpha : int; default=None
            the angle of attack profile (AEFACT) id
        lmach : int; default=None
            the mach profile (AEFACT) id
        lschd : int; default=None
            the control surface deflection profile (AEFACT) id
        comment : str; default=''
            a comment for the card

        """
        csschd = self.csschd.add(sid, aesurf_id, lschd, lalpha=lalpha, lmach=lmach,
                                 comment=comment)
        return csschd

    def add_aesurf(self, aesid: int, label: str,
                   cid1: int, aelist_id1: int,
                   cid2: Optional[int]=None, aelist_id2: Optional[int]=None,
                   eff: float=1.0, ldw: str='LDW',
                   crefc: float=1.0, crefs: float=1.0,
                   pllim: float=-np.pi/2., pulim: float=np.pi/2.,
                   hmllim=None, hmulim=None, # hinge moment lower/upper limits
                   tqllim=None, tqulim=None, # TABLEDi deflection limits vs. dynamic pressure
                   comment: str='') -> int:
        """
        Creates an AESURF card, which defines a control surface

        Parameters
        ----------
        aesid : int
            controller number
        label : str
            controller name
        cid1 / cid2 : int / None
            coordinate system id for primary/secondary control surface
        aelist_id1 / aelist_id2 : int / None
            AELIST id for primary/secondary control surface
        eff : float; default=1.0
            Control surface effectiveness
        ldw : str; default='LDW'
            Linear downwash flag;  ['LDW', 'NODLW']
        crefc : float; default=1.0
            reference chord for the control surface
        crefs : float; default=1.0
            reference area for the control surface
        pllim / pulim : float; default=-pi/2 / pi/2
            Lower/Upper deflection limits for the control surface in radians
        hmllim / hmulim : float; default=None
            Lower/Upper hinge moment limits for the control surface in
            force-length units
        tqllim / tqulim : int; default=None
            Set identification numbers of TABLEDi entries that provide the
            lower/upper deflection limits for the control surface as a
            function of the dynamic pressure
        comment : str; default=''
            a comment for the card

        """
        aesurf = self.aesurf.add(aesid, label, cid1, aelist_id1, cid2=cid2, aelist_id2=aelist_id2,
                                 eff=eff, ldw=ldw, crefc=crefc, crefs=crefs,
                                 pllim=pllim, pulim=pulim,
                                 hmllim=hmllim, hmulim=hmulim,
                                 tqllim=tqllim, tqulim=tqulim, comment=comment)
        return aesurf

    def add_aesurfs(self, aesurfs_id: int, label: str,
                    list1: int, list2: int, comment: str='') -> int:
        """
        Creates an AESURFS card

        Parameters
        ----------
        aesid : int
            the unique id
        label : str
            the AESURF name
        list1 / list2 : int / None
            the list (SET1) of node ids for the primary/secondary
            control surface(s) on the AESURF card
        comment : str; default=''
            a comment for the card

        """
        aesurfs = self.aesurfs.add(aesurfs_id, label, list1, list2, comment=comment)
        return aesurfs

    def add_aeparm(self, aeparm_id: int, label: str, units: str, comment: str='') -> int:
        """
        Creates an AEPARM card, which defines a new trim variable.

        Parameters
        ----------
        aeparm_id : int
            the unique id
        label : str
            the variable name
        units : str
            unused by Nastran
        comment : str; default=''
            a comment for the card

        """
        aeparm = self.aeparm.add(aeparm_id, label, units, comment=comment)
        return aeparm

    def add_aepress(self, mach, sym_xz: str, sym_xy: str, ux_id: int, dmij: str, dmiji: str):
        #AEPRESS MACH SYMXZ SYMXY UXID DMIJ DMIJI
        """adds an AEPRESS card"""
        assert isinstance(sym_xz, str), sym_xz
        assert isinstance(sym_xy, str), sym_xy
        assert isinstance(dmij, str), dmij
        assert isinstance(dmiji, str), dmiji
        fields = ['AEPRESS', mach, sym_xz, sym_xy, ux_id, dmij, dmiji]
        self.reject_card_lines('AEPRESS', print_card_(fields).split('\n'), show_log=False)

    def add_aeforce(self, mach: float, sym_xz: str, sym_xy: str, ux_id: int,
                    mesh: str, force: int, dmik: str, perq: str) -> None:
        """adds an AEPRESS card"""
        assert isinstance(mesh, str), mesh
        assert isinstance(sym_xz, str), sym_xz
        assert isinstance(sym_xy, str), sym_xy
        assert isinstance(dmik, str), dmik
        assert isinstance(perq, str), perq
        fields = ['AEFORCE', mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq]
        self.reject_card_lines('AEPRESS', print_card_(fields).split('\n'), show_log=False)


class AddOptimization(BDFAttributes):
    def add_deqatn(self, equation_id: int, eqs: list[str], comment: str='') -> DEQATN:
        """
        Creates a DEQATN card

        Parameters
        ----------
        equation_id : int
            the id of the equation
        eqs : list[str]
            the equations, which may overbound the field
            split them by a semicolon (;)
        comment : str; default=''
            a comment for the card

        DEQATN  41      F1(A,B,C,D,R) = A+B *C-(D**3 + 10.0) + sin(PI(1) * R)
                        + A**2 / (B - C); F = A + B - F1 * D

        def F1(A, B, C, D, R):
            F1 = A+B *C-(D**3 + 10.0) + sin(PI(1) * R) + A**2 / (B - C)
            F = A + B - F1 * D
            return F

        eqs = [
            'F1(A,B,C,D,R) = A+B *C-(D**3 + 10.0) + sin(PI(1) * R) + A**2 / (B - C)',
            'F = A + B - F1 * D',
        ]
        >>> deqatn = model.add_deqatn(41, eqs, comment: str='')

        """
        deqatn = DEQATN(equation_id, eqs, ifile=0, comment=comment)
        self._add_methods.add_deqatn_object(deqatn)
        return deqatn

    def add_desvar(self, desvar_id: int, label: str, xinit: float,
                   xlb: float=-1e20, xub: float=1e20,
                   delx=None, ddval: Optional[int]=None,
                   comment: str='') -> int:
        """
        Creates a DESVAR card

        Parameters
        ----------
        desvar_id : int
            design variable id
        label : str
            name of the design variable
        xinit : float
            the starting point value for the variable
        xlb : float; default=-1.e20
            the lower bound
        xub : float; default=1.e20
            the lower bound
        delx : float; default=1.e20
            fractional change allowed for design variables during
            approximate optimization
            NX  if blank : take from DOPTPRM; otherwise 1.0
            MSC if blank : take from DOPTPRM; otherwise 0.5
        ddval : int; default=None
            int : DDVAL id
                  allows you to set discrete values
            None : continuous
        comment : str; default=''
            a comment for the card

        """
        desvar = self.desvar.add(desvar_id, label, xinit, xlb, xub, delx=delx,
                                 ddval=ddval, comment=comment)
        return desvar

    def add_topvar(self, opt_id, label, ptype, xinit, pid, xlb=0.001, delxv=0.2,
                   power=3.0) -> int:
        """adds a TOPVAR"""
        topvar = TOPVAR(opt_id, label, ptype, xinit, pid, xlb=xlb, delxv=delxv, power=power)
        self._add_methods.add_topvar_object(topvar)
        return topvar

    def add_dresp1(self, dresp_id: int, label: str,
                   response_type: str, property_type: str, region: str,
                   atta: Optional[int | float | str],
                   attb: Optional[int | float | str],
                   atti: list[int | float | str],
                   validate: bool=True, comment: str='') -> int:
        """
        Creates a DRESP1 card.

        A DRESP1 is used to define a "simple" output result that may be
        optimized on.  A simple result is a result like stress, strain,
        force, displacement, eigenvalue, etc. for a node/element that
        may be found in a non-optimization case.

        Parameters
        ----------
        dresp_id : int
            response id
        label : str
            Name of the response
        response_type : str
            Response type
        property_type : str
            Element flag (PTYPE = 'ELEM'), or property entry name, or panel
            flag for ERP responses (PTYPE = 'PANEL' - See Remark 34), or
            RANDPS ID. Blank for grid point responses. 'ELEM' or property
            name used only with element type responses (stress, strain,
            force, etc.) to identify the relevant element IDs, or the property
            type and relevant property IDs.

            Must be {ELEM, PBAR, PSHELL, PCOMP, PANEL, etc.)
            PTYPE = RANDPS ID when RTYPE=PSDDISP, PSDVELO, or PSDACCL.
        region : str
            Region identifier for constraint screening
        atta : int / float / str / blank
            Response attribute
        attb : int / float / str / blank
            Response attribute
        atti : list[int / float / str]
            the response values to pull from
            list[int]:
                list of grid ids
                list of property ids
            list[str]
                'ALL'
        comment : str; default=''
            a comment for the card
        validate : bool; default=True
            should the card be validated when it's created

        Examples
        --------
        **Stress/PSHELL**

        >>> dresp_id = 103
        >>> label = 'resp1'
        >>> response_type = 'STRESS'
        >>> property_type = 'PSHELL'
        >>> pid = 3
        >>> atta = 9 # von mises upper surface stress
        >>> region = None
        >>> attb = None
        >>> atti = [pid]
        >>> DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)


        **Stress/PCOMP**

        >>> dresp_id = 104
        >>> label = 'resp2'
        >>> response_type = 'STRESS'
        >>> property_type = 'PCOMP'
        >>> pid = 3
        >>> layer = 4
        >>> atta = 9 # von mises upper surface stress
        >>> region = None
        >>> attb = layer
        >>> atti = [pid]
        >>> DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        """
        assert len(label) <= 8, label
        dresp = self.dresp1.add(
            dresp_id, label, response_type, property_type, region,
            atta, attb, atti, validate=validate, comment=comment)
        return dresp

    def add_dresp2(self, dresp_id: int, label: str, dequation: int, region: int,
                   params: dict[tuple[int, str], list[int]],
                   method: str='MIN',
                   c1: float=1., c2: float=0.005, c3: float=10.,
                   validate: bool=True, comment: str='') -> int:
        """
        Creates a DRESP2 card.

        A DRESP2 is used to define a "complex" output result that may be
        optimized on.  A complex result is a result that uses:
          - simple (DRESP1) results
          - complex (DRESP2) results
          - default values (DTABLE)
          - DVCRELx values
          - DVMRELx values
          - DVPRELx values
          - DESVAR values
        Then, an equation (DEQATN) is used to formulate an output response.

        Parameters
        ----------
        dresp_id : int
            response id
        label : str
            Name of the response
        dequation : int
            DEQATN id
        region : str
            Region identifier for constraint screening
        params : dict[(index, card_type)] = values
            the storage table for the response function
            index : int
                a counter
            card_type : str
                the type of card to pull from
                DESVAR, DVPREL1, DRESP2, etc.
            values : list[int]
                the values for this response
        method : str; default=MIN
            flag used for FUNC=BETA/MATCH
            FUNC = BETA
                valid options are {MIN, MAX}
            FUNC = MATCH
                valid options are {LS, BETA}
        c1 / c2 / c3 : float; default=1. / 0.005 / 10.0
            constants for FUNC=BETA or FUNC=MATCH
        comment : str; default=''
            a comment for the card
        validate : bool; default=False
            should the card be validated when it's created

        params = {
           (0, 'DRESP1') = [10, 20],
           (1, 'DESVAR') = [30],
           (2, 'DRESP1') = [40],
        }

        """
        assert len(label) <= 8, label
        dresp = self.dresp2.add(
            dresp_id, label, dequation, region, params,
            method=method, c1=c1, c2=c2, c3=c3, comment=comment,
            validate=validate)
        return dresp

    def add_dresp3(self, dresp_id, label, group, Type, region, params,
                   validate=True, comment: str='') -> DRESP3:
        """Creates a DRESP3 card"""
        dresp = DRESP3(dresp_id, label, group, Type, region, params,
                       validate=validate, comment=comment)
        self._add_methods.add_dresp_object(dresp)
        return dresp

    def add_dvcrel1(self, dvcrel_id: int, element_type: str, eid: int, cp_name: str,
                    desvar_ids: list[int], coeffs: list[float],
                    cp_min=None, cp_max: float=1e20, c0: float=0.,
                    validate: bool=True, comment: str='') -> int:
        """
        Creates a DVCREL1 card

        Parameters
        ----------
        element_type : int
            optimization id
        element_type : str
            element card name (e.g., CONM2)
        eid : int
            element id
        cp_name : str/int
            optimization parameter as an element connectivity name (e.g., X1)
        desvar_ids : list[int]
            DESVAR ids
        coeffs : list[float]
            scale factors for DESVAR ids
        cp_min : float; default=None
            minimum value
        cp_max : float; default=1e20
            maximum value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        """
        dvcrel = self.dvcrel1.add(
            dvcrel_id, element_type, eid, cp_name, desvar_ids, coeffs,
            cp_min=cp_min, cp_max=cp_max, c0=c0,
            validate=validate, comment=comment)
        return dvcrel

    def add_dvcrel2(self, dvcrel_id: int, element_type: str, eid: int, cp_name: str,
                    deqation, desvar_ids: list[int], labels: list[float],
                    cp_min=None, cp_max: float=1e20,
                    validate: bool=True, comment: str='') -> int:
        """Creates a DVCREL2 card"""
        dvcrel = self.dvcrel2.add(dvcrel_id, element_type, eid, cp_name,
                                  deqation, desvar_ids, labels,
                                  cp_min=cp_min, cp_max=cp_max,
                                  validate=validate, comment=comment)
        return dvcrel

    def add_dvprel1(self, dvprel_id: int, prop_type: str, pid: int, pname_fid: int | str,
                    desvar_ids: list[int],
                    coeffs: list[float],
                    p_min=None, p_max: float=1e20, c0: float=0.0,
                    validate: bool=True, comment: str='') -> int:
        """
        Creates a DVPREL1 card

        Parameters
        ----------
        dvprel_id : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        dvids : list[int]
            DESVAR ids
        coeffs : list[float]
            scale factors for DESVAR ids
        p_min : float; default=None
            minimum property value
        p_max : float; default=1e20
            maximum property value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        """
        dvprel = self.dvprel1.add(dvprel_id, prop_type, pid, pname_fid, desvar_ids, coeffs,
                                  p_min=p_min, p_max=p_max,
                                  c0=c0, validate=validate, comment=comment)
        return dvprel

    def add_dvprel2(self, dvprel_id: int, prop_type: str, pid: int,
                    pname_fid: int | str, deqation: int,
                    dvids: list[int]=None,
                    labels: list[str]=None,
                    p_min: Optional[float]=None, p_max: float=1.0e20,
                    validate: bool=True, comment: str='') -> int:
        """
        Creates a DVPREL2 card

        Parameters
        ----------
        dvprel_id : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        deqation : int
            DEQATN id
        dvids : list[int]; default=None
            DESVAR ids
        labels : list[str]; default=None
            DTABLE names
        p_min : float; default=None
            minimum property value
        p_max : float; default=1e20
            maximum property value
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        Notes
        -----
        either dvids or labels is required

        """
        dvprel = self.dvprel2.add(
            dvprel_id, prop_type, pid, pname_fid, deqation, dvids, labels,
            p_min=p_min, p_max=p_max, validate=validate, comment=comment)
        return dvprel

    def add_dvmrel1(self, oid: int, mat_type: str, mid: int, mp_name: str,
                    dvids: list[int], coeffs: list[float],
                    mp_min: Optional[float]=None, mp_max: float=1e20, c0: float=0.,
                    validate: bool=True, comment: str='') -> int:
        """
        Creates a DVMREL1 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            material card name (e.g., MAT1)
        mid : int
            material id
        mp_name : str
            optimization parameter as a pname (material name; E)
        dvids : list[int]
            DESVAR ids
        coeffs : list[float]
            scale factors for DESVAR ids
        mp_min : float; default=None
            minimum material property value
        mp_max : float; default=1e20
            maximum material property value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        """
        dvmrel = self.dvmrel1.add(
            oid, mat_type, mid, mp_name, dvids, coeffs,
            mp_min, mp_max, c0=c0, validate=validate, comment=comment)
        return dvmrel

    def add_dvmrel2(self, dvmrel_id: int, mat_type: str, mid: int, mp_name: str,
                    deqatn_id: int, desvar_ids: list[int], labels: list[str],
                    mp_min: Optional[float]=None, mp_max: float=1e20,
                    validate: bool=True, comment: str='') -> int:
        """
        Creates a DVMREL2 card

        Parameters
        ----------
        dvmrel_id : int
            optimization id
        mat_type : str
            material card name (e.g., MAT1)
        mid : int
            material id
        mp_name : str
            optimization parameter as a pname (material name; E)
        deqatn_id : int
            DEQATN id
        desvar_ids : list[int]; default=None
            DESVAR ids
        labels : list[str]; default=None
            DTABLE names
        mp_min : float; default=None
            minimum material property value
        mp_max : float; default=1e20
            maximum material property value
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        Notes
        -----
        either dvids or labels is required

        """
        dvmrel = self.dvmrel2.add(
            dvmrel_id, mat_type, mid, mp_name, deqatn_id, desvar_ids, labels,
            mp_min=mp_min, mp_max=mp_max,
            validate=validate, comment=comment)
        return dvmrel

    def add_dvgrid(self, dvid: int, nid: int, dxyz: NDArray3float,
                   cid: int=0, coeff: float=1.0, comment: str='') -> int:
        """
        Creates a DVGRID card

        Parameters
        ----------
        dvid : int
            DESVAR id
        nid : int
            GRID/POINT id
        dxyz : (3, ) float ndarray
            the amount to move the grid point
        cid : int; default=0
            Coordinate system for dxyz
        coeff : float; default=1.0
            the dxyz scale factor
        comment : str; default=''
            a comment for the card

        """
        dvgrid = self.dvgrid.add(dvid, nid, dxyz, cid=cid, coeff=coeff, comment=comment)
        return dvgrid

    def add_ddval(self, ddval_id: int, ddvals: list[int], comment: str='') -> int:
        """Creates a DDVAL card"""
        ddval = self.ddval.add(ddval_id, ddvals, comment=comment)
        return ddval

    def add_dlink(self, oid: int, dependent_desvar: int,
                  independent_desvars: list[int],
                  coeffs: list[float],
                  c0: float=0.0, cmult: float=1.0, comment: str='') -> int:
        """
        Creates a DLINK card, which creates a variable that is a lienar
        ccombination of other design variables

        Parameters
        ----------
        oid : int
            optimization id
        dependent_desvar : int
            the DESVAR to link
        independent_desvars : list[int]
            the DESVARs to combine
        coeffs : list[int]
            the linear combination coefficients
        c0 : float; default=0.0
            an offset
        cmult : float; default=1.0
            an scale factor
        comment : str; default=''
            a comment for the card

        """
        dlink = self.dlink.add(
            oid, dependent_desvar,
            independent_desvars, coeffs, c0=c0, cmult=cmult, comment=comment)
        return dlink

    def add_dconstr(self, dconstr_id: int, dresp_id: int,
                    lid: float=-1.e20, uid: float=1.e20,
                    lowfq: float=0.0, highfq: float=1.e20, comment: str='') -> int:
        """
        Creates a DCONSTR card

        Parameters
        ----------
        dconstr_id : int
            unique optimization id
        dresp_id : int
            DRESP1/2 id
        lid / uid=-1.e20 / 1.e20
            lower/upper bound
        lowfq / highfq : float; default=0. / 1.e20
            lower/upper end of the frequency range
        comment : str; default=''
            a comment for the card

        """
        dconstr = self.dconstr.add(dconstr_id, dresp_id, lid=lid, uid=uid, lowfq=lowfq,
                                   highfq=highfq, comment=comment)
        return dconstr

    def add_dconadd(self, dconstr_id: int, dconstrs: list[int], comment: str='') -> int:
        """Creates a DCONADD card"""
        dconadd = self.dconadd.add(dconstr_id, dconstrs, comment=comment)
        return dconadd

    def add_doptprm(self, params: dict[str, int | float], comment: str='') -> DOPTPRM:
        """Creates a DOPTPRM card"""
        doptprm = DOPTPRM(params, comment=comment)
        self._add_methods.add_doptprm_object(doptprm)
        return doptprm

    def add_dscreen(self, response_type: str, trs: float=-0.5, nstr: int=20,
                    comment: str='') -> int:
        """
        Creates a DSCREEN object

        Parameters
        ----------
        response_type : str
            Response type for which the screening criteria apply
        trs : float
            Truncation threshold
        nstr : int
            Maximum number of constraints to be retained per region per
            load case
        comment : str; default=''
            a comment for the card

        """
        dscreen = self.dscreen.add(response_type, trs=trs, nstr=nstr, comment=comment)
        return dscreen
    # ------------------------------------------
    # nx optimization
    def add_dvtrel1(self, dvtrel_id: int, label: str, group_id: int,
                    state: str='ACTIVE', dsv_flag: int=0, dvid1: int=0,
                    validate=False, comment: str='') -> int:
        dvtrel = DVTREL1(dvtrel_id, label, group_id, state=state,
                         dsv_flag=dsv_flag, dvid1=dvid1, comment=comment)
        self._add_methods.add_dvtrel_object(dvtrel)
        return dvtrel

    def add_dmncon(self, constraint_id: int, constraint_type: str,
                   xyz=None, normal=None, size=None,
                   m=None, d=None, nsections=None,
                   angle=None, mind=None, off_flag=None,
                   comment: str='') -> DMNCON:
        dmncon = DMNCON(constraint_id, constraint_type, xyz=xyz, normal=normal,
                        size=size, m=m, d=d, nsections=nsections,
                        angle=angle, mind=mind,
                        off_flag=off_flag, comment=comment)
        self._add_methods.add_dmncon_object(dmncon)
        return dmncon

    def add_dscons(self, dscid: int, label: str, constraint_type: str,
                   nid_eid: int, comp: int,
                   limit: float=0.0, opt: str='MAX', layer_id: int=1) -> None:
        """
        Design Constraint
        Defines a design constraint in design sensitivity analysis
        (original DSA). See the MSC.Nastran Reference Manual, Chapter 15.

        | DSCONS | DSCID |  LABEL | TYPE | ID | COMP | LIMIT | OPT | LAMNO |
        | DSCONS |   21  | G101DX | DISP |  4 |   1  |  0.06 | MAX |   6   |
        """
        assert opt in ['MIN', 'MAX'], opt
        assert isinstance(limit, float), limit
        assert isinstance(layer_id, int), layer_id
        fields = ['DSCONS', dscid, label, nid_eid, comp, limit, opt, layer_id]
        self.reject_card_lines('DSCONS', print_card_(fields).split('\n'), show_log=False)

    def add_dvset(self, vid: int, dv_type: str, field: int, pref: float, pids: list[float],
                  alpha: float=1.0) -> None:
        """
        Design Variable Set Property

        Defines a set of element properties that vary in a fixed relation
        to a design variable for design sensitivity analysis (original DSA).
        See the MSC.Nastran Reference Manual, Chapter 15.

        | DVSET | VID  |  TYPE  | FIELD | PREF | ALPHA | PIDl | PID2 | PID3 |
        |       | PID4 |  PID5  |  etc. |      |       |      |      |      |
        | DVSET |  21  | PSHELL |   4   | 0.20 |  1.0  |  99  |  101 |  110 |
        |       |  111 |   122  |       |      |       |      |      |      |

        | DVSET |  VID |  TYPE  | FIELD | PREF | ALPHA | PID  | THRU | PID2 |
        | DVSET |  21  | PSHELL |   4   | 0.20 |  1.0  |  101 | THRU |  105 |
        | DVSET |  VID |  TYPE  | FIELD | MIDV |       | PIDl | PID2 | PID3 |
        | DVSET |  21  | PSHELL |   3   | 134  |       |  87  |  101 |      |
        MSC 2001
        """
        fields = ['DVSET', vid, dv_type, field, pref, alpha] + pids
        self.reject_card_lines('DVSET', print_card_(fields).split('\n'), show_log=False)

    def add_dvar(self, bid: int, label: str, vids: list[float],
                 deltab: float=0.02) -> None:
        """
        Design Variable

        Defines a design variable for design sensitivity analysis
        (original DSA) described in the MSC.Nastran Reference Manual,
        Chapter 15.

        | DVAR |  BID |  LABEL | DELTAB | VID1 | VID2 | VID3 | VID4 | VID5 |
        |      | VID6 |  etc.  |        |      |      |      |      |      |
        | DVAR |  10  | LFDOOR |  0.01  |   2  |   4  |   5  |   6  |   9  |
        |      |  10  |        |        |      |      |      |      |      |
        MSC 2001

        Parameters
        ----------
        BID : int
            Design variable identification number. Must be unique for all DVAR.
        LABEL : str
            Label used to describe variable in output.
        DELTAB : float; default=0.02
            The change in the dimensionless design variable B to be used in the
            calculation of the design sensitivity coefficients.
        VIDi : list[int]
            Identification number of DVSET entry.

        """
        fields = ['DVAR', bid, label, deltab] + vids
        self.reject_card_lines('DVAR', print_card_(fields).split('\n'), show_log=False)


class AddMaterial(BDFAttributes):
    def add_creep(self, mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                  a, b, c, d, e, f, g, comment: str='') -> int:
        """Creates a CREEP card"""
        mat = CREEP(mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                    a, b, c, d, e, f, g, comment=comment)
        self._add_methods.add_creep_material_object(mat)
        return mat

    def add_mat1(self, mid: int,
                 E: Optional[float],
                 G: Optional[float],
                 nu: Optional[float],
                 rho: float=0.0, alpha: float=0.0, tref: float=0.0, ge: float=0.0,
                 St: float=0.0, Sc: float=0.0, Ss: float=0.0,
                 mcsid: int=0, comment: str='') -> int:
        """
        Creates a MAT1 card

        Parameters
        ----------
        mid : int
            material id
        E : float / None
            Young's modulus
        G : float / None
            Shear modulus
        nu : float / None
            Poisson's ratio
        rho : float; default=0.
            density
        a : float; default=0.
            coefficient of thermal expansion
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            damping coefficient
        St / Sc / Ss : float; default=0.
            tensile / compression / shear allowable
        mcsid : int; default=0
            material coordinate system id
            used by PARAM,CURV
        comment : str; default=''
            a comment for the card

        If E, G, or nu is None (only 1), it will be calculated

        """
        mat = self.mat1.add(
            mid, E, G, nu, rho=rho, alpha=alpha, tref=tref, ge=ge, St=St,
            Sc=Sc, Ss=Ss, mcsid=mcsid, comment=comment)
        return mat

    def add_mat2(self, mid: float, G11: float, G12: float, G13: float,
                 G22: float, G23: float, G33: float, rho: float=0.,
                 a1: Optional[float]=None, a2: Optional[float]=None, a3: Optional[float]=None,
                 tref: float=0., ge: float=0.,
                 St: Optional[float]=None, Sc: Optional[float]=None, Ss: Optional[float]=None,
                 mcsid: Optional[int]=None, comment: str='') -> int:
        """Creates an MAT2 card"""
        mat = self.mat2.add(mid, G11, G12, G13, G22, G23, G33,
                            rho, a1, a2, a3,
                            tref=tref, ge=ge, St=St, Sc=Sc,
                            Ss=Ss, mcsid=mcsid, comment=comment)
        return mat

    def add_mat3(self, mid: int, ex: float, eth: float, ez: float,
                 nuxth: float, nuthz: float, nuzx: float,
                 rho: float=0.0, gzx: Optional[float]=None,
                 ax: float=0., ath: float=0., az: float=0.,
                 tref: float=0., ge: float=0.,
                 comment: str='') -> int:
        """Creates a MAT3 card"""
        mat = self.mat3.add(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho=rho, gzx=gzx,
                            ax=ax, ath=ath, az=az, tref=tref, ge=ge,
                            comment=comment)
        return mat

    def add_mat4(self, mid: int, k: float, cp: float=0.0, rho: float=1.0,
                 H: Optional[float]=None, mu: Optional[float]=None, hgen: float=1.0,
                 ref_enthalpy: Optional[float]=None, tch: Optional[float]=None, tdelta: Optional[float]=None,
                 qlat: Optional[float]=None, comment: str='') -> int:
        """Creates a MAT4 card"""
        mat = self.mat4.add(mid, k, cp=cp, rho=rho, H=H, mu=mu, hgen=hgen,
                            ref_enthalpy=ref_enthalpy, tch=tch, tdelta=tdelta,
                            qlat=qlat, comment=comment)
        return mat

    def add_mat5(self, mid: int, kxx: float=0., kxy: float=0., kxz: float=0.,
                 kyy: float=0., kyz: float=0., kzz: float=0., cp: float=0.,
                 rho: float=1., hgen: float=1., comment: str='') -> int:
        """Creates a MAT5 card"""
        mat = self.mat5.add(mid, kxx=kxx, kxy=kxy, kxz=kxz, kyy=kyy,
                            kyz=kyz, kzz=kzz, cp=cp,
                            rho=rho, hgen=hgen, comment=comment)
        return mat

    def add_mat8(self, mid: int, e11: float, e22: float, nu12: float,
                 g12: float=0.0, g1z: float=1e8, g2z: float=1e8,
                 rho: float=0., a1: float=0., a2: float=0., tref: float=0.,
                 Xt: float=0., Xc: Optional[float]=None,
                 Yt: float=0., Yc: Optional[float]=None,
                 S: float=0., ge: float=0., F12: float=0., strn: float=0.,
                 comment: str='') -> int:
        """Creates a MAT8 card"""
        mat = self.mat8.add(mid, e11, e22, nu12, g12, g1z, g2z, rho=rho, a1=a1, a2=a2,
                            tref=tref, xt=Xt, xc=Xc, yt=Yt, yc=Yc,
                            s=S, ge=ge, f12=F12, strn=strn, comment=comment)
        return mat

    def add_mat9(self, mid: int,
                 G11=0., G12=0., G13=0., G14=0., G15=0., G16=0.,
                 G22=0., G23=0., G24=0., G25=0., G26=0.,
                 G33=0., G34=0., G35=0., G36=0.,
                 G44=0., G45=0., G46=0.,
                 G55=0., G56=0., G66=0.,
                 rho=0., A=None, tref=0., ge=0., comment: str='') -> int:
        """Creates a MAT9 card"""
        mat = self.mat9.add(mid, G11, G12, G13, G14, G15, G16,
                            G22, G23, G24, G25, G26,
                            G33, G34, G35, G36, G44, G45, G46, G55,
                            G56, G66, rho, A, tref, ge, comment=comment)
        return mat

    def add_mat10(self, mid: int, bulk: float, rho: float, c: float,
                  ge: float=0.0,
                  gamma: Optional[float]=None,
                  table_bulk: Optional[int]=None,
                  table_rho: Optional[int]=None,
                  table_ge: Optional[int]=None,
                  table_gamma: Optional[int]=None,
                  comment: str='') -> int:
        """
        Creates a MAT10 card

        Parameters
        ----------
        mid : int
            material id
        bulk : float; default=None
            Bulk modulus
        rho : float; default=None
            Density
        c : float; default=None
            Speed of sound
        ge : float; default=0.
            Damping
        gamma : float; default=None
            NX : ratio of imaginary bulk modulus to real bulk modulus; default=0.0
            MSC : normalized admittance coefficient for porous material
        table_bulk : int; default=None
            TABLEDx entry defining bulk modulus vs. frequency
            None for MSC Nastran
        table_rho : int; default=None
            TABLEDx entry defining rho vs. frequency
            None for MSC Nastran
        table_ge : int; default=None
            TABLEDx entry defining ge vs. frequency
            None for MSC Nastran
        table_gamma : int; default=None
            TABLEDx entry defining gamma vs. frequency
            None for MSC Nastran
        comment : str; default=''
            a comment for the card

        """
        mat = self.mat10.add(mid, bulk, rho, c, ge=ge, gamma=gamma,
                             table_bulk=table_bulk, table_rho=table_rho,
                             table_ge=table_ge, table_gamma=table_gamma,
                             comment=comment)
        return mat

    def add_mat11(self, mid: int, e1: float, e2: float, e3: float,
                  nu12: float, nu13: float, nu23: float,
                  g12: float, g13: float, g23: float,
                  rho: float=0.0,
                  a1: float=0.0, a2: float=0.0, a3: float=0.0,
                  tref: float=0.0, ge: float=0.0, comment: str='') -> int:
        """Creates a MAT11 card"""
        mat = self.mat11.add(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=rho,
                             a1=a1, a2=a2, a3=a3, tref=tref, ge=ge, comment=comment)
        return mat

    def add_mat10c(self, mid: int, form: str='REAL',
                   rho_real: float=0.0, rho_imag: float=0.0,
                   c_real: float=0.0, c_imag: float=0.0,
                   comment: str='') -> int:
        mat = self.mat10c.add(mid, form=form, rho_real=rho_real, rho_imag=rho_imag,
                              c_real=c_real, c_imag=c_imag, comment=comment)
        return mat

    def add_matort(self, mid: int, E1: float, E2: float, E3: float,
                   nu12: float, nu23: float, nu31: float,
                   G12: float, G23: float, G31: float,
                   rho: float=0.0,
                   alpha1: float=0.0, alpha2: float=0.0, alpha3: float=0.0,
                   tref: float=0.0, ge: float=0.0,
                   comment: str=''):
        mat = self.matort.add(
            mid, E1, E2, E3, nu12, nu23, nu31, G12, G23, G31,
            rho=rho, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3,
            tref=tref, ge=ge,
            comment=comment)
        return mat

    def add_mat3d(self, mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=0.0,
                  comment: str='') -> int:
        """
        This is a VABS specific card that is almost identical to the MAT11.
        """
        mat = MAT3D(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=rho,
                    comment=comment)
        self._add_methods.add_structural_material_object(mat)
        return mat

    def add_matg(self, mid, idmem, behav, tabld, tablu, yprs, epl, gpl, gap=0.,
                 tab_yprs=None, tab_epl=None,
                 tab_gpl=None, tab_gap=None, comment: str='') -> int:
        """Creates a MATG card"""
        mat = MATG(mid, idmem, behav, tabld, tablu, yprs, epl, gpl, gap=gap,
                   tab_yprs=tab_yprs, tab_epl=tab_epl,
                   tab_gpl=tab_gpl, tab_gap=tab_gap, comment=comment)
        self._add_methods.add_structural_material_object(mat)
        return mat

    def add_mathe(self, mid: int, model: str, bulk: float,
                  mus, alphas, betas,
                  mooney, sussbat, aboyce, gent,
                  rho: float=0., texp: float=0., tref: float=0., ge: float=0.,
                  comment: str='') -> int:
        """Creates a MATHE card, which is a hyperelastic material"""
        if model == 'MOONEY':
            (c10, c01,
             c20, c11, c02,
             c30, c21, c12, c03, ) = mooney
            mat = self.mathe.add_mooney(
                mid, bulk, rho=rho, texp=texp,
                c10=c10, c01=c01,
                c20=c20, c11=c11, c02=c02,
                c30=c30, c21=c21, c12=c12, c03=c03,
                comment=comment)
        else:  # pragma: no cover
            mat = self.mathe.add(mid, model, bulk, mus, alphas, betas,
                                 mooney, sussbat, aboyce, gent,
                                 rho=rho, texp=texp, tref=tref, ge=ge,
                                 comment=comment)
        return mat

    def add_mathp(self, mid, a10=0., a01=0., d1=None, rho=0., av=0., tref=0., ge=0., na=1, nd=1,
                  a20=0., a11=0., a02=0., d2=0.,
                  a30=0., a21=0., a12=0., a03=0., d3=0.,
                  a40=0., a31=0., a22=0., a13=0., a04=0., d4=0.,
                  a50=0., a41=0., a32=0., a23=0., a14=0., a05=0., d5=0.,
                  tab1=None, tab2=None, tab3=None, tab4=None, tabd=None, comment: str='') -> int:
        """Creates a MATHP card"""
        mat = self.mathp.add(mid, a10, a01, d1, rho, av, tref, ge, na, nd,
                             a20, a11, a02, d2,
                             a30, a21, a12, a03, d3,
                             a40, a31, a22, a13, a04,
                             d4, a50, a41, a32, a23, a14, a05, d5, tab1, tab2, tab3,
                             tab4, tabd, comment=comment)
        return mat

    def add_mats1(self, mid: int, Type: str,
                  h: float, hr: float, yf: float, limit1: float, limit2: float,
                  tid: int=0, comment: str='') -> int:
        """Creates a MATS1 card"""
        mat = self.mats1.add(mid, Type, h, hr, yf, limit1, limit2,
                             tid=tid, comment=comment)
        return mat

    def add_matt1(self, mid: int, e_table=None, g_table=None, nu_table=None, rho_table=None,
                  a_table=None, ge_table=None, st_table=None, sc_table=None, ss_table=None,
                  comment: str='') -> int:
        """Creates a MATT1 card"""
        mat = self.matt1.add(
            mid, e_table, g_table, nu_table, rho_table, a_table,
            ge_table, st_table, sc_table, ss_table,
            comment=comment)
        return mat

    def add_matt2(self, mid, g11_table=None, g12_table=None, g13_table=None, g22_table=None,
                  g23_table=None, g33_table=None, rho_table=None,
                  a1_table=None, a2_table=None, a3_table=None, ge_table=None,
                  st_table=None, sc_table=None, ss_table=None,
                  comment: str='') -> int:
        """Creates a MATT2 card"""
        mat = self.matt2.add(mid,
                             g11_table, g22_table, g33_table,
                             g12_table, g13_table, g23_table,
                             rho_table,
                             a1_table, a2_table, a3_table,
                             st_table, sc_table, ss_table, ge_table,
                             comment=comment)
        return mat

    def add_matt3(self, mid: int,
                  ex_table: int=0, eth_table: int=0, ez_table: int=0,
                  nuth_table: int=0, nuxz_table: int=0, rho_table: int=0,
                  gzx_table: int=0,
                  ax_table: int=0, ath_table: int=0, az_table: int=0,
                  ge_table: int=0, comment: str='') -> int:
        """Creates a MATT3 card"""
        mat = self.matt3.add(mid, ex_table, eth_table, ez_table,
                             nuth_table, nuxz_table, rho_table, gzx_table,
                             ax_table, ath_table, az_table, ge_table, comment=comment)
        return mat

    def add_matt4(self, mid: int,
                  k_table: int=0, cp_table: int=0, h_table: int=0,
                  mu_table: int=0, hgen_table: int=0, comment: str='') -> int:
        """Creates a MATT4 card"""
        mat = self.matt4.add(mid, k_table, cp_table, h_table, mu_table, hgen_table,
                             comment=comment)
        return mat

    def add_matt5(self, mid, kxx_table: int=0, kxy_table: int=0, kxz_table: int=0,
                  kyy_table: int=0, kyz_table: int=0, kzz_table: int=0,
                  cp_table: int=0, hgen_table: int=0, comment: str='') -> int:
        """Creates a MATT5 card"""
        mat = self.matt5.add(
            mid, kxx_table, kxy_table, kxz_table, kyy_table,
            kyz_table, kzz_table, cp_table,
            hgen_table, comment=comment)
        return mat

    def add_matt8(self, mid: int, e1_table=None, e2_table=None, nu12_table=None,
                  g12_table=None, g1z_table=None, g2z_table=None, rho_table=None,
                  a1_table=None, a2_table=None,
                  xt_table=None, xc_table=None, yt_table=None, yc_table=None,
                  s_table=None, ge_table=None, f12_table=None, comment: str='') -> int:
        """Creates a MATT8 card"""
        mat = self.matt8.add(mid, e1_table=e1_table, e2_table=e2_table, nu12_table=nu12_table,
                             g12_table=g12_table, g1z_table=g1z_table, g2z_table=g2z_table,
                             rho_table=rho_table, a1_table=a1_table, a2_table=a2_table,
                             xt_table=xt_table, xc_table=xc_table, yt_table=yt_table, yc_table=yc_table,
                             s_table=s_table, ge_table=ge_table, f12_table=f12_table, comment=comment)
        return mat

    def add_matt9(self, mid: int,
                  g11_table=None, g12_table=None, g13_table=None, g14_table=None,
                  g15_table=None, g16_table=None, g22_table=None, g23_table=None,
                  g24_table=None, g25_table=None, g26_table=None, g33_table=None,
                  g34_table=None, g35_table=None, g36_table=None, g44_table=None,
                  g45_table=None, g46_table=None, g55_table=None, g56_table=None,
                  g66_table=None, rho_table=None,
                  a1_table=None, a2_table=None, a3_table=None,
                  a4_table=None, a5_table=None, a6_table=None,
                  ge_table=None, comment: str='') -> int:
        mat = self.matt9.add(
            mid,
            g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
            g22_table, g23_table, g24_table, g25_table, g26_table,
            g33_table, g34_table, g35_table, g36_table,
            g44_table, g45_table, g46_table,
            g55_table, g56_table,
            g66_table, rho_table,
            a1_table, a2_table, a3_table, a4_table, a5_table, a6_table,
            ge_table, comment=comment)
        return mat


class AddContact(BDFAttributes):
    def add_bsurf(self, glue_id: int, eids: list[int], comment: str='') -> int:
        """Creates a BSURF card, which defines a contact/glue source/target for shells"""
        bsurf = self.bsurf.add(glue_id, eids, comment=comment)
        return bsurf

    def add_bsurfs(self, glue_id: int, eids: list[int], comment: str='') -> int:
        """Creates a BSURFS card, which defines a contact/glue source/target for solids"""
        bsurfs = self.bsurfs.add(glue_id, eids, comment=comment)
        return bsurfs

    def add_bctset(self, contact_id: int, source_ids: list[int],
                   target_ids: list[int],
                   frictions: list[float],
                   min_distances: list[float],
                   max_distances: list[float],
                   desc_id: int=0,
                   comment: str='') -> int:
        """Creates a BCTSET card, which defines a generalized contact set"""
        bctset = self.bctset.add(contact_id, source_ids, target_ids,
                                 frictions, min_distances,
                                 max_distances, desc_id=desc_id, comment=comment)
        return bctset

    def add_bctadd(self, contact_id: int, contact_sets: list[int], comment: str='') -> int:
        """Creates a BCTADD card"""
        bctadd = self.bctadd.add(contact_id, contact_sets, comment=comment)
        return bctadd

    def add_bctpara(self, csid, params: dict[str, Any], comment: str='') -> BCTPARA:
        """Creates a BCTPARA card"""
        bctpara = BCTPARA(csid, params, comment=comment)
        self._add_methods.add_bctpara_object(bctpara)
        return bctpara

    def add_bctparm(self, contact_id, params: dict[str, Any], comment: str='') -> BCTPARM:
        """Creates a BCTPARA card"""
        bctpara = BCTPARM(contact_id, params, comment=comment)
        self._add_methods.add_bctparm_object(bctpara)
        return bctpara

    def add_blseg(self, line_id: int, nodes: list[int], comment: str='') -> int:
        """Creates a BLSEG card"""
        blseg = self.blseg.add(line_id, nodes, comment=comment)
        return blseg

    def add_bconp(self, contact_id: int, slave: int, master: int,
                  fric_id: int, sfac: float=1.0, ptype: int=1, cid: int=0,
                  comment: str='') -> int:
        """Creates a BCONP card"""
        bconp = self.bconp.add(contact_id, slave, master, sfac, fric_id, ptype,
                               cid, comment=comment)
        return bconp

    def add_bfric(self, friction_id: int, mu1: float, fstiff: float=None, comment: str=''):
        bfric = self.bfric.add(friction_id, mu1, fstiff=fstiff, comment=comment)
        return bfric

    def add_bcrpara(self, contact_region_id: int,
                    surf: str='TOP',
                    offset: Optional[float]=None,
                    surface_type: str='FLEX',
                    grid_point: int=0,
                    comment: str='') -> int:
        """
        Creates a BCRPARA card

        Parameters
        ----------
        crid : int
            CRID Contact region ID.
        offset : float; default=None
            Offset distance for the contact region (Real > 0.0).
            None : OFFSET value in BCTPARA entry
        surf : str; default='TOP'
            SURF Indicates the contact side. See Remark 1.  {'TOP', 'BOT'; )
        surface_type : str; default='FLEX'
            Indicates whether a contact region is a rigid surface if it
            is used as a target region. {'RIGID', 'FLEX'}.
            This is not supported for SOL 101.
        grid_point : int; default=0
            Control grid point for a target contact region with TYPE=RIGID
            or when the rigid-target algorithm is used.  The grid point
            may be used to control the motion of a rigid surface.
            (Integer > 0).  This is not supported for SOL 101.
        comment : str; default=''
            a comment for the card

        """
        bcrpara = self.bcrpara.add(
            contact_region_id, surf=surf, offset=offset,
            surface_type=surface_type, grid_point=grid_point,
            comment=comment)
        return bcrpara

    def add_bedge(self, bedge_id: int,
            eids: list[int],
            grids: list[tuple[int, int]],
            comment: str='') -> int:
        """Creates a BEDGE card"""
        bcrpara = self.bedge.add(bedge_id, eids, grids, comment=comment)
        return bcrpara

class AddSuperelements(BDFAttributes):
    def add_sebset(self, seid: int, ids: list[int], components: int | list[int], comment: str='') -> int:
        """Creates an SEBSET/SEBSET1 card"""
        ifile = 0
        if isinstance(components, integer_string_types):
            sebset = self.sebset.add_set1(seid, ids, components, ifile=ifile, comment=comment)
        else:
            sebset = self.sebset.add_set(seid, ids, components, ifile=ifile, comment=comment)
        return sebset

    def add_sebset1(self, seid, ids, components, comment: str='') -> int:
        """.. seealso:: ``add_secset``"""
        return self.add_sebset(seid, ids, components, comment=comment)

    def add_secset(self, seid: int, ids: list[int], components: int | list[int], comment: str='') -> int:
        """Creates an SECSET/SECSET1 card"""
        ifile = 0
        if isinstance(components, integer_string_types):
            secset = self.sebset.add_set1(seid, ids, components, ifile=ifile, comment=comment)
        else:
            secset = self.secset.add_set(seid, ids, components, ifile=ifile, comment=comment)
        return secset

    def add_secset1(self, seid: int, ids: list[int], components: int, comment: str='') -> int:
        """.. seealso:: ``add_secset``"""
        return self.add_secset(seid, ids, components, comment=comment)

    def add_seqset(self, seid, ids, components, comment: str='') -> int:
        """Creates an SEQSET card"""
        ifile = 0
        if isinstance(components, integer_string_types):
            seqset = self.seqset.add_set1(seid, ids, components, ifile=ifile, comment=comment)
        else:
            seqset = self.seqset.add_set(seid, ids, components, ifile=ifile, comment=comment)
        return seqset

    def add_seqset1(self, seid, ids, components, comment: str='') -> int:
        """.. seealso:: ``add_secset``"""
        return self.add_seqset(seid, ids, components, comment=comment)

    def add_seset(self, seid: int, node_ids: list[int], comment: str='') -> int:
        """Creates an SESET card"""
        seset = self.seset.add(seid, node_ids, comment=comment)
        return seset

    def add_sesup(self, nodes, Cs, comment: str='') -> SESUP:
        """Creates an SESUP card"""
        se_suport = SESUP(nodes, Cs, comment=comment)
        self._add_methods.add_sesuport_object(se_suport)
        return se_suport

    def add_release(self, seid, comp: int, nids: list[int], comment: str='') -> int:
        release = self.release.add(seid, comp, nids, comment=comment)
        return release

    def add_sebndry(self, seid_a, seid_b, ids, comment: str='') -> SEBNDRY:
        sebndry = SEBNDRY(seid_a, seid_b, ids, comment=comment)
        self._add_methods.add_sebndry_object(sebndry)
        return sebndry

    def add_sebulk(self, seid, superelement_type, rseid,
                   method='AUTO', tol=1e-5, loc='YES', unitno=None, comment: str='') -> SEBULK:
        sebulk = SEBULK(seid, superelement_type, rseid,
                        method=method, tol=tol, loc=loc,
                        unitno=unitno, comment=comment)
        self._add_methods.add_sebulk_object(sebulk)
        return sebulk

    def add_seconct(self, seid_a: int, seid_b: int, tol: float, loc: str,
                    nodes_a: list[int], nodes_b: list[int], comment: str='') -> SECONCT:
        seconct = SECONCT(seid_a, seid_b, tol, loc, nodes_a, nodes_b,
                          comment=comment)
        self._add_methods.add_seconct_object(seconct)
        return seconct

    def add_seelt(self, seid, ids, comment: str='') -> SEELT:
        seelt = SEELT(seid, ids, comment=comment)
        self._add_methods.add_seelt_object(seelt)
        return seelt

    def add_seexcld(self, seid_a, seid_b, nodes, comment: str='') -> SEEXCLD:
        seexcld = SEEXCLD(seid_a, seid_b, nodes, comment=comment)
        self._add_methods.add_seexcld_object(seexcld)
        return seexcld

    def add_selabel(self, seid, label, comment: str='') -> SELABEL:
        selabel = SELABEL(seid, label, comment=comment)
        self._add_methods.add_selabel_object(selabel)
        return selabel

    def add_seloc(self, seid, nodes_seid, nodes0, comment: str='') -> SELOC:
        """
        Creates an SELOC card, which transforms the superelement SEID
        from PA to PB.  Basically, define two CORD1Rs.

        Parameters
        ----------
        seid : int
            the superelement to transform
        nodes_seid : list[int, int, int]
            the nodes in the superelement than define the resulting coordinate system
        nodes0 : list[int, int, int]
            the nodes in the superelement than define the starting coordinate system
        comment : str; default=''
            a comment for the card

        """
        seloc = SELOC(seid, nodes_seid, nodes0, comment=comment)
        self._add_methods.add_seloc_object(seloc)
        return seloc

    def add_seload(self, lid_s0, seid, lid_se, comment: str='') -> SELOAD:
        seload = SELOAD(lid_s0, seid, lid_se, comment=comment)
        self._add_methods.add_seload_object(seload)
        return seload

    def add_sempln(self, seid, p1, p2, p3, comment: str='') -> SEMPLN:
        sempln = SEMPLN(seid, p1, p2, p3, comment=comment)
        self._add_methods.add_sempln_object(sempln)
        return sempln

    def add_setree(self, seid, seids, comment: str='') -> SETREE:
        setree = SETREE(seid, seids, comment=comment)
        self._add_methods.add_setree_object(setree)
        return setree

    def add_csuper(self, seid, psid, nodes, comment: str='') -> CSUPER:
        csuper = CSUPER(seid, psid, nodes, comment=comment)
        self._add_methods.add_csuper_object(csuper)
        return csuper

    def add_csupext(self, seid, nodes, comment: str='') -> CSUPEXT:
        csupext = CSUPEXT(seid, nodes, comment=comment)
        self._add_methods.add_csupext_object(csupext)
        return csupext

    def add_senqset(self, set_id, n, comment: str='') -> SENQSET:
        senqset = SENQSET(set_id, n, comment=comment)
        self._add_methods.add_senqset_object(senqset)
        return senqset

    def add_extrn(self, nids: list[int], comps: list[str]) -> None:
        fields = ['EXTRN']
        for nid, comp in zip(nids, comps):
            fields.extend([nid, comp])
        self.reject_card_lines('EXTRN', print_card_(fields).split('\n'), show_log=False)


class AddThermal(BDFAttributes):
    def add_temp(self, load_id: int, temperature_dict: dict[int, float],
                 comment: str='') -> int:
        """
        Creates a TEMP card

        Parameters
        ----------
        load_id : int
            Load set identification number
        temperature_dict : dict[nid] : temperature
            nid : int
                node id
            temperature : float
                the nodal temperature
        comment : str; default=''
            a comment for the card

        """
        temp = self.temp.add(load_id, temperature_dict, comment=comment)
        return temp

    #def add_tempp1(self) -> TEMPP1:
        #temp = TEMPP1()
        #self._add_thermal_load_object(temp)
        #return temp

    def add_tempd(self, sid: int, temperature: float, comment: str='') -> int:
        """
        Creates a TEMPD card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        temperature : float
            default temperature
        comment : str; default=''
            a comment for the card

        """
        tempd = self.tempd.add(sid, temperature, comment=comment)
        return tempd

    def add_qhbdy(self, sid: int, flag: str, q0: float, grids: list[int],
                  area_factor: Optional[float]=None, comment: str='') -> int:
        """
        Creates a QHBDY card

        Parameters
        ----------
        sid : int
            load id
        flag : str
            valid_flags = {POINT, LINE, REV, AREA3, AREA4, AREA6, AREA8}
        q0 : float
            Magnitude of thermal flux into face. Q0 is positive for heat
            into the surface
        area_factor : float; default=None
            Area factor depends on type
        grids : list[int]
            Grid point identification of connected grid points
        comment : str; default=''
            a comment for the card

        """
        load = self.qhbdy.add(sid, flag, q0, grids, area_factor=area_factor, comment=comment)
        return load

    def add_qbdy1(self, sid: int, qflux: float, eids: list[int],
                  comment: str='') -> int:
        """Creates a QBDY1 card"""
        load = self.qbdy1.add(sid, qflux, eids, comment=comment)
        return load

    def add_qbdy2(self, sid: int, eid: int, qfluxs: list[float], comment: str='') -> int:
        """Creates a QBDY1 card"""
        load = self.qbdy2.add(sid, eid, qfluxs, comment=comment)
        return load

    def add_qbdy3(self, sid: int, q0, eids: list[int],
                  cntrlnd: int=0, comment: str='') -> int:
        """
        Creates a QBDY3 card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        q0 : float; default=None
            Magnitude of thermal flux vector into face
        control_id : int; default=0
            Control point
        eids : list[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        load = self.qbdy3.add(sid, q0, eids, cntrlnd=cntrlnd, comment=comment)
        return load

    def add_qvol(self, sid: int, qvol: float,
                 control_point: int, elements: list[int],
                 comment: str='') -> int:
        """Creates a QVOL card"""
        load = self.qvol.add(sid, qvol, control_point, elements, comment=comment)
        return load

    def add_qvect(self, sid: int, q0: float, eids: list[int],
                  t_source: float=None,
                  ce: int=0,
                  vector_tableds: list[int | float]=0.0,
                  control_id: int=0, comment: str='') -> int:
        """
        Creates a QVECT card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        q0 : float; default=None
            Magnitude of thermal flux vector into face
        t_source : float; default=None
            Temperature of the radiant source
        ce : int; default=0
            Coordinate system identification number for thermal vector flux
        vector_tableds : list[int/float, int/float, int/float]
            vector : float; default=0.0
                directional cosines in coordinate system CE) of
                the thermal vector flux
            tabled : int
                TABLEDi entry identification numbers defining the
                components as a function of time
        control_id : int; default=0
            Control point
        eids : list[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        load = self.qvect.add(
            sid, q0, eids, t_source=t_source, ce=ce,
            vector_tableds=vector_tableds, control_id=control_id,
            comment=comment)
        return load

    def add_chbdyg(self, eid: int, surface_type: str, nodes: list[int],
            iview_front: int=0, iview_back: int=0,
            rad_mid_front: int=0, rad_mid_back: int=0,
            comment: str='') -> int:
        """Creates a CHBDYG card"""
        elem = self.chbdyg.add(
            eid, surface_type, nodes,
            iview_front=iview_front, iview_back=iview_back,
            rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
            comment=comment)
        return elem

    def add_chbdyp(self, eid: int, pid: int, surface_type: str,
                   g1: int, g2: int,
                   g0: int=0, gmid: int=0,
                   ce: int=0,
                   iview_front: int=0, iview_back: int=0,
                   rad_mid_front: int=0, rad_mid_back: int=0,
                   e1=None, e2=None, e3=None,
                   comment: str='') -> int:
        """
        Creates a CHBDYP card

        Parameters
        ----------
        eid : int
            Surface element ID
        pid : int
            PHBDY property entry identification numbers. (Integer > 0)
        surface_type : str
            Surface type
            Must be {POINT, LINE, ELCYL, FTUBE, TUBE}
        iview_front : int; default=0
            A VIEW entry identification number for the front face.
        iview_back : int; default=0
            A VIEW entry identification number for the back face.
        g1 / g2 : int
            Grid point identification numbers of grids bounding the surface
        g0 : int; default=0
            Orientation grid point
        rad_mid_front : int
            RADM identification number for front face of surface element
        rad_mid_back : int
            RADM identification number for back face of surface element.
        gmid : int
            Grid point identification number of a midside node if it is used
            with the line type surface element.
        ce : int; default=0
            Coordinate system for defining orientation vector
        e1 / e2 / e3 : float; default=None
            Components of the orientation vector in coordinate system CE.
            The origin of the orientation vector is grid point G1.
        comment : str; default=''
            a comment for the card

        """
        elem = self.chbdyp.add(
            eid, pid, surface_type, g1, g2,
            g0=g0, gmid=gmid, ce=ce,
            iview_front=iview_front, iview_back=iview_back,
            rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
            e1=e1, e2=e2, e3=e3,
            comment=comment)
        return elem

    def add_chbdye(self, eid: int, eid2: int, side: int,
                   iview_front: int=0, iview_back: int=0,
                   rad_mid_front: int=0, rad_mid_back: int=0,
                   comment: str='') -> int:
        """
        Creates a CHBDYE card

        Parameters
        ----------
        eid : int
            surface element ID number for a side of an element
        eid2: int
            a heat conduction element identification
        side: int
            a consistent element side identification number (1-6)
        iview_front: int; default=0
            a VIEW entry identification number for the front face
        iview_back: int; default=0
            a VIEW entry identification number for the back face
        rad_mid_front: int; default=0
            RADM identification number for front face of surface element
        rad_mid_back: int; default=0
            RADM identification number for back face of surface element
        comment : str; default=''
            a comment for the card

        """
        elem = self.chbdye.add(
            eid, eid2, side,
            iview_front=iview_front, iview_back=iview_back,
            rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
            comment=comment)
        return elem

    def add_phbdy(self, pid: int,
                  area_factor: float=None,
                  d1=None, d2=None,
                  comment: str='') -> int:
        """
        Creates a PHBDY card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        area_factor: float
            Area factor of the surface used only for CHBDYP element
            Must be {POINT, LINE, TUBE, ELCYL}
            TUBE : constant thickness of hollow tube
        d1, d2 : float; default=None
            Diameters associated with the surface
            Used with CHBDYP [ELCYL, TUBE, FTUBE] surface elements
        comment : str; default=''
            a comment for the card

        """
        prop = self.phbdy.add(pid, area_factor=area_factor,
                              d1=d1, d2=d2, comment=comment)
        return prop

    def add_conv(self, eid: int, pconid: int,
                 ta: list[int], film_node: int=0, cntrlnd: int=0,
                 comment: str='') -> int:
        """
        Creates a CONV card

        Parameters
        ----------
        eid : int
            element id
        pconid : int
            Convection property ID
        mid : int
            Material ID
        ta : list[int]
            Ambient points used for convection 0's are allowed for TA2
            and higher
        film_node : int; default=0
            Point for film convection fluid property temperature
        cntrlnd : int; default=0
            Control point for free convection boundary condition
        comment : str; default=''
            a comment for the card

        """
        boundary_condition = self.conv.add(
            eid, pconid, ta,
            film_node=film_node, cntrlnd=cntrlnd,
            comment=comment)
        return boundary_condition

    def add_convm(self, eid: int, pconvm: int, ta1: int,
                  film_node: int=0, cntmdot: int=0,
                  ta2: int=0, mdot: float=1.0,
                  comment: str='') -> int:
        """
        Creates a CONVM card

        Parameters
        ----------
        eid : int
            element id (CHBDYP)
        pconid : int
            property ID (PCONVM)
        mid : int
            Material ID
        ta1 : int
            ambient point for convection
        ta2 : int; default=None
            None : ta1
            ambient point for convection
        film_node : int; default=0
        cntmdot : int; default=0
            control point used for controlling mass flow
            0/blank is only allowed when mdot > 0
        mdot : float; default=1.0
            a multiplier for the mass flow rate in case there is no
            point associated with the CNTRLND field
            required if cntmdot = 0
        comment : str; default=''
            a comment for the card

        """
        boundary_condition = self.convm.add(
            eid, pconvm, ta1,
            film_node=film_node, cntmdot=cntmdot,
            ta2=ta2, mdot=mdot,
            comment=comment)
        return boundary_condition

    def add_tempbc(self, sid: int, bc_type: str,
                   nodes: list[int], temps: list[float],
                   comment: str='') -> int:
        tempbc = self.tempbc.add(sid, bc_type, nodes, temps, comment=comment)
        return tempbc

    def add_radm(self, radmid, absorb, emissivity, comment: str='') -> int:
        """Creates a RADM card"""
        boundary_condition = self.radm.add(radmid, absorb, emissivity, comment=comment)
        return boundary_condition

    def add_radbc(self, nodamb, famb, cntrlnd, eids, comment: str='') -> int:
        """Creates a RADBC card"""
        boundary_condition = self.radbc.add(nodamb, famb, cntrlnd, eids, comment=comment)
        return boundary_condition

    def add_view(self, iview: int, icavity: int, shade: str='BOTH',
                 nbeta: int=1, ngamma: int=1, dislin: float=0.0, comment: str='') -> int:
        """Creates a VIEW card"""
        view = self.view.add(iview, icavity, shade=shade, nbeta=nbeta, ngamma=ngamma,
                             dislin=dislin, comment=comment)
        return view

    def add_view3d(self, icavity: int,
                   gitb: int=4, gips: int=4, cier: int=4,
                   error_tol: float=0.1, zero_tol: float=1e-10, warp_tol: float=0.01,
                   rad_check: int=3, comment: str='') -> int:
        """Creates a VIEW3D card"""
        view3d = self.view3d.add(icavity, gitb=gitb, gips=gips, cier=cier,
                                 error_tol=error_tol, zero_tol=zero_tol, warp_tol=warp_tol,
                                 rad_check=rad_check, comment=comment)
        return view3d

    def add_pconv(self, pconid: int, mid: int=None,
                  form: int=0,
                  exponent_free_convection=0.0,
                  free_convection_type: int=0,
                  table_id: int=0,
                  chlen: float=None, gidin: int=None, ce: int=0,
                  e1: float=None, e2: float=None, e3: float=None,
                  comment: str='') -> int:
        """
        Creates a PCONV card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid : int; default=None
            Material ID
        form : int; default=0
            Type of formula used for free convection
            Must be {0, 1, 10, 11, 20, or 21}
        expf : float; default=0.0
            Free convection exponent as implemented within the context
            of the particular form that is chosen
        ftype : int; default=0
            Formula type for various configurations of free convection
        tid : int; default=None
            Identification number of a TABLEHT entry that specifies the
            two variable tabular function of the free convection heat
            transfer coefficient
        chlen : float; default=None
            Characteristic length
        gidin : int; default=None
            Grid ID of the referenced inlet point
        ce : int; default=0
            Coordinate system for defining orientation vector.
        e1 / e2 / e3 : list[float]; default=None
            Components of the orientation vector in coordinate system CE.
            The origin of the orientation vector is grid point G1
        comment : str; default=''
            a comment for the card

        """
        prop = self.pconv.add(
            pconid, mid=mid,
            form=form, exponent_free_convection=exponent_free_convection,
            free_convection_type=free_convection_type,
            table_id=table_id, chlen=chlen, gidin=gidin,
            ce=ce, e1=e1, e2=e2, e3=e3, comment=comment)
        return prop

    def add_pconvm(self, pconid: int, mid: int, coeff: float,
                   form: int=0, flag: int=0,
                   expr: float=0.0,
                   exppi: float=0.0, exppo: float=0.0,
                   comment: str='') -> int:
        """
        Creates a PCONVM card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid: int
            Material ID
        coeff: float
            Constant coefficient used for forced convection
        form: int; default=0
            Type of formula used for free convection
            Must be {0, 1, 10, 11, 20, or 21}
        flag: int; default=0
            Flag for mass flow convection
        expr: float; default=0.0
            Reynolds number convection exponent
        exppi: float; default=0.0
            Prandtl number convection exponent for heat transfer into
            the working fluid
        exppo: float; default=0.0
            Prandtl number convection exponent for heat transfer out of
            the working fluid
        comment : str; default=''
            a comment for the card

        """
        prop = self.pconvm.add(
            pconid, mid, coeff, form=form, flag=flag,
            expr=expr, exppi=exppi, exppo=exppo, comment=comment)
        return prop

    def add_radset(self, cavity_ids: list[int], comment: str='') -> int:
        """Creates a RADSET card"""
        set_obj = self.radset.add(cavity_ids, comment=comment)
        return set_obj


class AddAcoustic(BDFAttributes):
    def add_acmodl(self, infor, fset, sset,
                   normal=0.5, olvpang=60., search_unit='REL', intol=0.2,
                   area_op=0, ctype='STRONG', method='BW',
                   sk_neps=0.5, dsk_neps=0.75, all_set='NO', inter='DIFF',
                   nastran_version='nx', comment: str='') -> ACMODL:
        acmodl = ACMODL(infor, fset, sset,
                        normal=normal, olvpang=olvpang, search_unit=search_unit,
                        intol=intol, area_op=area_op, ctype=ctype, method=method,
                        sk_neps=sk_neps, dsk_neps=dsk_neps, all_set=all_set,
                        inter=inter, nastran_version=nastran_version, comment=comment)
        self._add_methods.add_acmodl_object(acmodl)
        return acmodl

    def add_chacab(self, eid, pid, nodes, comment: str='') -> int:
        chacab = self.chacab.add(eid, pid, nodes, comment=comment)
        return chacab

    def add_caabsf(self, eid: int, pid: int, nodes: list[int],
                   comment: str='') -> int:
        caabsf = self.caabsf.add(eid, pid, nodes, comment=comment)
        return caabsf

    def add_chacbr(self, eid, pid, nodes, comment: str='') -> int:
        chacbr = self.chacbr.add(eid, pid, nodes, comment=comment)
        return chacbr

    def add_paabsf(self, pid: int,
                   table_reactance_real: int=0,
                   table_reactance_imag: int=0,
                   s: float=1.0, a: float=1.0, b: float=0.0,
                   k: float=0.0, rhoc: float=1.0,
                   comment: str=''):
        paabsf = self.paabsf.add(
            pid, table_reactance_real=table_reactance_real,
            table_reactance_imag=table_reactance_imag,
            s=s, a=a, b=b, k=k, rhoc=rhoc, comment=comment)
        return paabsf

    def add_pacbar(self, pid: int, mback: float, mseptm: float,
                   freson: float, kreson: float, comment: str='') -> int:
        """
        Creates a PACBAR card

        Parameters
        ----------
        pid : int
            Property identification number. (Integer > 0)
        mback : float
            Mass per unit area of the backing material
        mseptm : float
            Mass per unit area of the septum material
        freson : float; default=None
            Resonant frequency of the sandwich construction in hertz.
        kreson : float; default=None
            Resonant stiffness of the sandwich construction.

        """
        pacbar = self.pacbar.add(pid, mback, mseptm, freson, kreson, comment=comment)
        return pacbar

    def add_pacabs(self, pid, cutfr, b, k, m,
                   synth=True, tid_resistance=None, tid_reactance=None,
                   tid_weight=None, comment: str='') -> PACBAR:
        """
        Creates a PACABS card

        Parameters
        ----------
        pid : int
            Property identification number.
        synth : bool; default=True
            Request the calculation of B, K, and M from the tables TIDi below
        tid_resistance : int; default=None
            Identification of the TABLEDi entry that defines the resistance.
        tid_reactance : int; default=None
            Identification of the TABLEDi entry that defines the reactance.
        tid_weight : int; default=None
            Identification of the TABLEDi entry that defines the weighting function.
        cutfr : float
            Cutoff frequency for tables referenced above. (Real > 0.0)
        B, K, M : float
            Equivalent damping, stiffness and mass values per unit area. (Real > 0.0)

        ..note:: tables are defined as a function of frequency in cycles/time
        """
        pacabs = PACABS(pid, cutfr, b, k, m,
                        synth=synth, tid_resistance=tid_resistance,
                        tid_reactance=tid_reactance, tid_weight=tid_weight,
                        comment=comment)
        self._add_methods.add_acoustic_property_object(pacabs)
        return pacabs


class AddCards(AddCoords, Add0dElements, Add1dElements, Add2dElements, Add3dElements,
               AddRigidElements, AddAero, AddOptimization, AddMaterial, AddContact,
               AddSuperelements, AddThermal, AddBolts, AddAcoustic):
    def __init__(self):
        #BDFAttributes.__init__(self)
        super().__init__()
        self._add_methods = AddMethods(self)
        self.is_superelements = False
        self._is_cards_dict = True
        self._solmap_to_value = {}
        self.rsolmap_to_str = {}
        self._type_to_id_map = defaultdict(list)

        self._iparse_errors = 0
        self._stored_parse_errors = []
        self.values_to_skip = None
        self._duplicate_elements = None
        self._duplicate_properties = None
        self._duplicate_masses = None
        self._duplicate_materials = None
        self._duplicate_coords = None
        self._duplicate_thermal_materials = None
        self._stop_on_parsing_error = True
        self._stop_on_xref_error = True

    def parse_cards(self) -> None:
        cards = self._cards_to_setup
        log = self.log
        for card in cards:
            try:
                class_name = card.type
            except:  # pragma: no cover
                print(card)
                raise
            #print(class_name)
            if not (class_name == 'COORD' or
                    class_name in self._card_parser_prepare):
                msg = (
                    f'{class_name} has not been added to the list of cards to load,\n'
                    'but is included in the set of cards to process\n'
                    'add it to bdf.py into _card_parse_prepare \n'
                )
                if card not in self.cards_to_read and card not in {'ASET'}:
                    msg += 'add it to self.cards_to_read'
                raise RuntimeError(msg)
            assert isinstance(card.n, int), f'{card.type} n={card.n} type={type(card.n)}'
            #print(class_name, card.n)
            if card.n > 0:
                try:
                    card.parse_cards()
                    #log.info(class_name)
                except Exception as error:
                    log.error(str(error))
                    log.error(class_name)
                    raise
                ncard = len(card)
                icard = np.arange(ncard)
                card.slice_card_by_index(icard)

    def setup(self, run_setup: bool=True,
              run_geom_check: bool=True) -> None:
        #if not run_setup:
            #return
        self.parse_cards()
        self._check_duplicates()
        if not run_geom_check:
            return
        #missing = {
            #'node_id': np.array([], dtype='int32'),
            #'coord_id': np.array([], dtype='int32'),
            #'property_id': np.array([], dtype='int32'),
            #'material_id': np.array([], dtype='int32'),
        #}

        self._geom_check()

    def _check_duplicates(self) -> None:
        self.spring_element_ids
        self.damper_element_ids
        self.shell_element_ids
        self.solid_element_ids
        self.bar_property_ids
        self.beam_property_ids
        self.shell_property_ids
        self.solid_property_ids

    def _geom_check(self) -> None:
        log = self.log
        NO_GEOM_CHECK = {'GRID', 'COORD', 'SPOINT', 'EPOINT'}
        cards = [card for card in self._cards_to_setup if card.n > 0]

        for card in cards:
            #missing = {
                #'node_id': np.array([], dtype='int32'),
                #'coord_id': np.array([], dtype='int32'),
                #'property_id': np.array([], dtype='int32'),
                #'material_id': np.array([], dtype='int32'),
            #}
            class_name = card.type
            if class_name not in NO_GEOM_CHECK:
                missing = {
                    #'node_id': np.array([], dtype='int32'),
                    #'coord_id': np.array([], dtype='int32'),
                    #'property_id': np.array([], dtype='int32'),
                    #'material_id': np.array([], dtype='int32'),
                }
                if not hasattr(card, 'geom_check'):
                    self.log.error(f"{class_name!r} has no method 'geom_check'")
                    continue
                try:
                    card.geom_check(missing)
                #except (AttributeError, ValueError) as e:
                    #self.log.error(str(e))
                except (AttributeError, ValueError) as e:
                    log.error(str(e))
                    raise
                for key, value in missing.items():
                    if len(value):
                        log.warning(f'{card.type}: {key}={value} referenced by and is missing')

        #for k, v in missing.items():
            #if len(v):
                #log.warning(f'{k}={v} referenced by and is missing')

    @property
    def subcases(self) -> dict[int, Optional[Any]]:
        """gets the subcases"""
        if self.case_control_deck is None:
            return {}
        return self.case_control_deck.subcases

    def _get_maps(self, *args, **kwargs):
        self.log.warning('no _get_maps')

    def get_encoding(self, encoding: Optional[str]=None) -> str:
        """gets the file encoding"""
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
            if encoding is None:
                encoding = sys.getdefaultencoding()
            elif isinstance(encoding, bytes):
                # needed for hdf5 loader for some reason...
                encoding = self._encoding.decode('latin1')
                self._encoding = encoding
        encoding = cast(str, encoding)  # just for typing
        assert isinstance(encoding, str), encoding
        return encoding

    def cross_reference(self, run_setup: bool=True,
                        run_geom_check: bool=True):
        self.setup(run_setup=run_setup,
                   run_geom_check=run_geom_check)
        #self.log.warning('no xref')

    def add_param(self, key: str, values: list[int | float | str],
                  comment: str='') -> PARAM:
        return self._add_param_nastran(key, values, comment=comment)

    def _add_param_nastran(self, key: str, values: list[int | float | str],
                           comment: str='') -> PARAM:
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/list
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        param = PARAM(key, values, comment=comment)
        self._add_methods.add_param_object(param)
        return param

    # -----------------------------------------------------------------------------------
    def add_grid(self, nid: int, xyz: None | list[float] | NDArray3float,
                 cp: int=0, cd: int=0, ps: int=0,
                 seid: int=0, comment: str='') -> int:
        """
        Creates the GRID card

        Parameters
        ----------
        nid : int
            node id
        cp : int; default=0
            the xyz coordinate frame
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=0
            Additional SPCs in the analysis coordinate frame (e.g. 123).
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card

        """
        grid = self.grid.add(nid, xyz, cp=cp, cd=cd, ps=ps, seid=seid, comment=comment)
        return grid

    def add_grdset(self, cp: int, cd: int, ps: str, seid: int, comment: str='') -> GRDSET:
        """
        Creates the GRDSET card

        Parameters
        ----------
        cp : int; default=0
            the xyz coordinate frame
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card

        """
        grdset = GRDSET(cp, cd, ps, seid, comment=comment)
        self.grdset = grdset
        return grdset

    def add_seqgp(self, nids: list[int], seqids: list[int | float],
                  comment: str='') -> SEQGP:
        """
        Creates the SEQGP card

        Parameters
        ----------
        nids : int
           the node ids
        seqid : int/float
           the superelement id
        comment : str; default=''
            a comment for the card

        """
        seqgp = SEQGP(nids, seqids, comment=comment)
        self._add_methods.add_seqgp_object(seqgp)
        return seqgp

    def add_spoint(self, ids: int | list[int], comment: str='') -> int:
        """
        Creates the SPOINTs card that contains many SPOINTs

        Parameters
        ----------
        ids : list[int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        spoint = self.spoint.add(ids, comment=comment)
        return spoint

    def add_epoint(self, ids: int | list[int], comment: str='') -> int:
        """
        Creates the EPOINTs card that contains many EPOINTs

        Parameters
        ----------
        ids : list[int]
            EPOINT ids
        comment : str; default=''
            a comment for the card

        """
        epoint = self.epoint.add(ids, comment=comment)
        return epoint

    def add_point(self, nid: int, xyz: Any, cp: int=0, comment: str='') -> int:
        """
        Creates the POINT card

        Parameters
        ----------
        nid : int
            node id
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cp : int; default=0
            coordinate system for the xyz location
        comment : str; default=''
            a comment for the card

        """
        point = self.point.add(nid, xyz, cp=cp, comment=comment)
        return point

    #def add_cgen(self, Type, field_eid, pid, field_id, th_geom_opt,
                 #eidl, eidh, t_abcd=None, direction='L', comment: str='') -> CGEN:
        #"""Creates a CGEN element card"""
        #elem = CGEN(Type, field_eid, pid, field_id, th_geom_opt,
                    #eidl, eidh, t_abcd=t_abcd, direction=direction, comment=comment)
        #self.add_element_object(elem)
        #return elem

    #def add_gmload(self, sid, normal, entity, entity_id, method, load_magnitudes,
                   #cid=0, comment: str='') -> GMLOAD:
        #"""Creates a GMLOAD load card"""
        #load = GMLOAD(sid, normal, entity, entity_id, method, load_magnitudes,
                      #cid=cid, comment=comment)
        #self.add_load_object(load)
        #return load

    def add_param(self, key: str, values: list[int | float | str],
                  comment: str='') -> PARAM:
        return self._add_param_nastran(key, values, comment=comment)

    def _add_param_nastran(self, key: str, values: list[int | float | str],
                           comment: str='') -> PARAM:
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/list
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        param = PARAM(key, values, comment=comment)
        self._add_methods.add_param_object(param)
        return param

    def add_mdlprm(self, mdlprm_dict: dict[str, int | float],
                   comment: str='') -> MDLPRM:
        """
        Creates a MDLPRM card

        Parameters
        ----------
        mdlprm_dict : dict[name, value]
            name : str
                the name of the MDLPRM
            value: int/float
                varies depending on the type of MDLPRM
        comment : str; default=''
            a comment for the card

        """
        mdlprm = MDLPRM(mdlprm_dict, comment=comment)
        self._add_methods.add_mdlprm_object(mdlprm)
        return mdlprm


    def _add_param_mystran(self, key: str, values: list[int | float | str],
                           comment: str='') -> PARAM_MYSTRAN:
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/list
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        param = PARAM_MYSTRAN(key, values, comment=comment)
        self._add_methods.add_param_object(param)
        return param

    def add_plotel(self, eid: int, nodes: list[int], comment: str='') -> int:
        """
        Adds a PLOTEL card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int]
            Unique GRID point IDs
        comment : str; default=''
            a comment for the card

        """
        elem = self.plotel.add(eid, nodes, comment=comment)
        return elem

    def add_nsm(self, sid: int, nsm_type: str, pid_eid: int, value: float,
                comment: str='') -> int:
        """
        Creates an NSM card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        pid_eid : list[int]; int
            property id or element id depending on nsm_type
        value : list[float]; float
            the non-structural pass per unit length/area
            same length as pid_eid
        comment : str; default=''
            a comment for the card

        """
        if isinstance(pid_eid, int):
            pid_eid = [pid_eid]
        if isinstance(value, float):
            value = [value]
        assert isinstance(pid_eid, list), pid_eid
        assert isinstance(value, list), value
        assert len(pid_eid) == len(value), 'len(pid_eid)=%s len(value)=%s' % (len(pid_eid), len(value))

        #nsms = []
        nsm = self.nsm.add(sid, nsm_type, pid_eid, value, comment=comment)
        #for pid_eidi, valuei in zip(pid_eid, value):
            #nsm = NSM(sid, nsm_type, pid_eidi, valuei, comment=comment)
            #self._add_methods.add_nsm_object(nsm)
            #nsms.append(nsm)
        return nsm

    def add_nsm1(self, sid: int, nsm_type: str, value: float, ids: list[int],
                 comment: str='') -> int:
        """
        Creates an NSM1 card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        value : float
            the non-structural pass per unit length/area
        ids : list[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card

        """
        assert isinstance(value, (float, str)), f'value={value} type={type(value)}'
        nsm = self.nsm1.add(sid, nsm_type, value, ids, comment=comment)
        return nsm

    def add_nsml(self, sid: int, nsm_type: str, pid_eid: int, value: float,
                 comment: str='') -> int:
        """
        Creates an NSML card, which defines lumped non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        pid_eid : list[int]; int
            property id or element id depending on nsm_type
        value : list[float]; float
            the non-structural pass per unit length/area
            same length as pid_eid
        comment : str; default=''
            a comment for the card

        """
        if isinstance(pid_eid, int):
            pid_eid = [pid_eid]
        if isinstance(value, float):
            value = [value]
        assert isinstance(pid_eid, list), pid_eid
        assert isinstance(value, list), value
        assert len(pid_eid) == len(value), 'len(pid_eid)=%s len(value)=%s' % (len(pid_eid), len(value))

        nsm = self.nsml.add(sid, nsm_type, pid_eid, value, comment=comment)
        #nsms = []
        #for pid_eidi, valuei in zip(pid_eid, value):
            #nsm = self.nsml.add(sid, nsm_type, pid_eidi, valuei, comment=comment)
            #nsms.append(nsm)
        return nsm

    def add_nsml1(self, sid: int, nsm_type: str, value: float, ids: list[int],
                  comment: str='') -> int:
        """
        Creates an NSML1 card, which defines lumped non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        value : float
            the non-structural pass per unit length/area
        ids : list[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card

        """
        nsm = self.nsml1.add(sid, nsm_type, value, ids, comment=comment)
        return nsm

    def add_nsmadd(self, sid: int, sets: list[int], comment: str='') -> int:
        """
        Creates an NSMADD card, which sum NSM sets

        Parameters
        ----------
        sid : int
            the NSM Case Control value
        sets : list[int]
            the NSM, NSM1, NSML, NSML1 values
        comment : str; default=''
            a comment for the card

        """
        nsmadd = self.nsmadd.add(sid, sets, comment=comment)
        return nsmadd

    def add_plplane(self, pid, mid, cid=0, stress_strain_output_location='GRID',
                    comment: str='') -> int:
        """Creates a PLPLANE card"""
        prop = self.plplane.add(
            pid, mid, cid=cid,
            stress_strain_output_location=stress_strain_output_location,
            comment=comment)
        return prop

    def add_pplane(self, pid: int, mid: int, t: float=0.0, nsm: float=0.0,
                   formulation_option: int=0, comment: str='') -> int:
        """Creates a PPLANE card"""
        prop = self.pplane.add(
            pid, mid, t=t, nsm=nsm,
            formulation_option=formulation_option, comment=comment)
        return prop

    def add_pshln1(self, pid: int, mid1: int=0, mid2: int=0, analysis: str='ISH',
                   behx=None, integration=None, behxh=None, integration_h=None,
                   comment: str='') -> int:
        prop = self.pshln1.add(
            pid, mid1=mid1, mid2=mid2, analysis=analysis,
            behx=behx, integration=integration,
            behxh=behxh, integration_h=integration_h, comment='')
        return prop

    def add_pshln2(self, pid: int, mid: int, direct: int=1, thickness: float=1.0,
                   analysis='ISH',
                   behx: list[str]=None,
                   integration: list[str]=None,
                   behxh: list[str]=None,
                   integration_h: list[str]=None,
                   comment: str='') -> int:
        prop = self.pshln2.add(
            pid, mid, direct=1, thickness=1.0, analysis='ISH',
            behx=None, integration=None, behxh=None, integration_h=None, comment='')
        return prop

    def add_cplstn3(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    comment: str='') -> int:
        """Creates a CPLSTN4 card"""
        elem = self.cplstn3.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_cplstn4(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    comment: str='') -> int:
        """Creates a CPLSTN4 card"""
        elem = self.cplstn4.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_cplstn6(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    comment: str='') -> int:
        """Creates a CPLSTN6 card"""
        elem = self.cplstn6.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_cplstn8(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    comment: str='') -> int:
        """Creates a CPLSTN8 card"""
        elem = self.cplstn8.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_cplsts3(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    tflag: int=0, T1=None, T2=None, T3=None, comment: str='') -> int:
        """Creates a CPLSTS3 card"""
        elem = self.cplsts3.add(
            eid, pid, nids, theta=theta,
            tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        return elem

    def add_cplsts4(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    tflag: int=0, T1=None, T2=None, T3=None, T4=None,
                    comment: str='') -> int:
        """Creates a CPLSTS4 card"""
        elem = self.cplsts4.add(
            eid, pid, nids, theta=theta,
            tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        return elem

    def add_cplsts6(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    tflag=0, thickness=None,
                    comment: str='') -> int:
        """Creates a CPLSTS6 card"""
        elem = self.cplsts6.add(eid, pid, nids, theta=theta,
                                tflag=tflag, thickness=thickness,
                                comment=comment)
        return elem

    def add_cplsts8(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
                    tflag: int=0,
                    thickness=None, comment: str='') -> int:
        """Creates a CPLSTS8 card"""
        elem = self.cplsts8.add(eid, pid, nids, theta=theta,
                                tflag=tflag, thickness=thickness,
                                comment=comment)
        return elem

    def add_crac2d(self, eid, pid, nids, comment: str='') -> int:
        """Creates a PRAC2D card"""
        elem = CRAC2D(eid, pid, nids, comment=comment)
        return elem

    def add_prac2d(self, pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180.,
                   comment: str='') -> int:
        """Creates a PRAC2D card"""
        prop = PRAC2D(pid, mid, thick, iplane, nsm=nsm, gamma=gamma, phi=phi,
                      comment=comment)
        return prop

    def add_crac3d(self, eid, pid, nids, comment: str='') -> int:
        """Creates a CRAC3D card"""
        elem = CRAC3D(eid, pid, nids, comment=comment)
        return elem

    def add_prac3d(self, pid, mid, gamma=0.5, phi=180., comment: str='') -> int:
        """Creates a PRAC3D card"""
        prop = PRAC3D(pid, mid, gamma=gamma, phi=phi, comment=comment)
        return prop

    def add_genel_stiffness(self, eid, ul, ud, k, s=None) -> int:
        """creates a GENEL card using the stiffness (K) approach"""
        assert k is not None
        genel = self.genel.add(eid, ul, ud, k, None, s)
        return genel

    def add_genel_flexibility(self, eid, ul, ud, z, s=None) -> int:
        """creates a GENEL card using the flexiblity (Z) approach"""
        assert z is not None
        genel = self.genel.add(eid, ul, ud, None, z, s)
        return genel

    #def add_axic(self, nharmonics, comment: str='') -> AXIC:
        #"""Creates a AXIC card"""
        #axic = AXIC(nharmonics, comment=comment)
        #self._add_methods.add_axic_object(axic)
        #return axic

    #def add_pointax(self, nid, ringax, phi, comment: str='') -> POINTAX:
        #"""Creates a POINTAX card"""
        #node = POINTAX(nid, ringax, phi, comment=comment)
        #self._add_methods.add_ringax_object(node)
        #return node

    #def add_ringax(self, nid, R, z, ps=None, comment: str='') -> RINGAX:
        #"""Creates a RINGAX card"""
        #node = RINGAX(nid, R, z, ps=ps, comment=comment)
        #self._add_methods.add_ringax_object(node)
        #return node

    #def add_cconeax(self, eid, pid, rings, comment: str='') -> CCONEAX:
        #"""Creates a CCONEAX card"""
        #elem = CCONEAX(eid, pid, rings, comment=comment)
        #self._add_methods.add_element_object(elem)
        #return elem

    #def add_pconeax(self, pid, mid1, t1=None, mid2=0, i=None, mid3=None, t2=None,
                    #nsm=0., z1=None, z2=None, phi=None, comment: str='') -> PCONEAX:
        #"""Creates a PCONEAX card"""
        #prop = PCONEAX(pid, mid1, t1, mid2, i, mid3, t2, nsm, z1, z2, phi,
                       #comment=comment)
        #self._add_methods.add_property_object(prop)
        #return prop

    def add_tempax(self, sid, ring, phi, temperature, comment: str='') -> int:
        """Creates a TEMPAX card"""
        load = self.tempax.add(sid, ring, phi, temperature, comment=comment)
        return load

    def add_presax(self, sid, pressure, rid1, rid2, phi1=0., phi2=360., comment: str='') -> int:
        """Creates a PRESAX card"""
        load = self.presax.add(sid, pressure, rid1, rid2, phi1, phi2, comment=comment)
        return load

    def add_forceax(self, sid, ring_id, hid, scale, f_rtz, comment: str='') -> int:
        forceax = self.forceax.add(sid, ring_id, hid, scale, f_rtz, comment=comment)
        return forceax

    def add_ctrax3(self, eid, pid, nids, theta=0., comment: str='') -> int:
        """Creates a CTRAX3 card"""
        elem = self.ctrax3.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_ctrax6(self, eid, pid, nids, theta=0., comment: str='') -> int:
        """Creates a CTRAX6 card"""
        elem = self.ctrax6.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_ctriax(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: int|float=0., comment: str='') -> int:
        """Creates a CTRIAX card"""
        elem = self.ctriax.add(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)
        return elem

    def add_ctriax6(self, eid: int, mid: int, nids: list[int], theta: float=0.,
                    comment: str='') -> int:
        """Creates a CTRIAX6 card"""
        elem = self.ctriax6.add(eid, mid, nids, theta=theta, comment=comment)
        return elem

    def add_cquadx(self, eid: int, pid: int, nids: list[int],
                   theta_mcid: int | float=0., comment: str='') -> int:
        """Creates a CQUADX card"""
        elem = self.cquadx.add(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)
        return elem

    def add_cquadx4(self, eid: int, pid: int, nids: list[int],
                    theta: float=0., comment: str='') -> int:
        """Creates a CQUADX4 card"""
        elem = self.cquadx4.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_cquadx8(self, eid: int, pid: int, nids: list[int],
                    theta: float=0., comment: str='') -> int:
        """Creates a CQUADX8 card"""
        elem = self.cquadx8.add(eid, pid, nids, theta=theta, comment=comment)
        return elem

    def add_nxstrat(self, sid, params, comment: str='') -> NXSTRAT:
        """Creates an NXSTRAT card"""
        nxstrat = NXSTRAT(sid, params, comment=comment)
        self._add_methods.add_nxstrat_object(nxstrat)
        return nxstrat

    def add_load(self, sid: int, scale: float,
                 scale_factors: list[float],
                 load_ids: list[int], comment: str='') -> int:
        """
        Creates a LOAD card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            overall scale factor
        scale_factors : list[float]
            individual scale factors (corresponds to load_ids)
        load_ids : list[int]
            individual load_ids (corresponds to scale_factors)
        comment : str; default=''
            a comment for the card

        """
        load = self.load.add(sid, scale, scale_factors, load_ids, comment=comment)
        return load

    def add_cload(self, sid: int, scale: float,
                  scale_factors: list[float], load_ids: list[int],
                  comment: str='') -> int:
        """
        Creates a CLOAD card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            overall scale factor
        scale_factors : list[float]
            individual scale factors (corresponds to load_ids)
        load_ids : list[int]
            individual load_ids (corresponds to scale_factors)
        comment : str; default=''
            a comment for the card

        """
        load = self.cload.add(sid, scale, scale_factors, load_ids, comment=comment)
        return load

    def add_lseq(self, sid: int, excite_id: int,
                 lid: int=0, tid: int=0,
                 comment: str='') -> int:
        """
        Creates a LSEQ card

        Parameters
        ----------
        sid : int
            loadset id; LOADSET points to this
        excite_id : int
            set id assigned to this static load vector
        lid : int
            load set id of a set of static load entries;
            LOAD in the Case Control
        tid : int; default=None
            temperature set id of a set of thermal load entries;
            TEMP(LOAD) in the Case Control
        comment : str; default=''
            a comment for the card

        """
        load = self.lseq.add(sid, excite_id, lid, tid=tid, comment=comment)
        return load

    def add_sload(self, sid: int, nids: list[int], mags: list[float],
                  comment: str='') -> int:
        """
        Creates an SLOAD (GRID/SPOINT load)

        Parameters
        ----------
        sid : int
            load id
        nids : int; list[int]
            the GRID/SPOINT ids
        mags : float; list[float]
            the load magnitude
        comment : str; default=''
            a comment for the card

        """
        sload = self.sload.add(sid, nids, mags, comment=comment)
        return sload

    def add_dload(self, sid: int, scale: float,
                  scale_factors: list[float], load_ids: list[int],
                  comment: str='') -> int:
        """
        Creates a DLOAD card

        Parameters
        ----------
        sid : int
            Load set identification number. See Remarks 1. and 4. (Integer > 0)
        scale : float
            Scale factor. See Remarks 2. and 8. (Real)
        Si : list[float]
            Scale factors. See Remarks 2., 7. and 8. (Real)
        load_ids : list[int]
            Load set identification numbers of RLOAD1, RLOAD2, TLOAD1,
            TLOAD2, and ACSRCE entries. See Remarks 3. and 7. (Integer > 0)
        comment : str; default=''
            a comment for the card

        """
        dload = self.dload.add(sid, scale, scale_factors, load_ids, comment=comment)
        return dload

    def add_darea(self, sid: int, nid: int, component: str, scale: float,
                  comment: str='') -> int:
        """
        Creates a DAREA card

        Parameters
        ----------
        sid : int
            darea id
        nid : int
            GRID, EPOINT, SPOINT id
        component : str
            Component number. (0-6; 0-EPOINT/SPOINT; 1-6 GRID)
        scale : float
            Scale (area) factor

        """
        darea = self.darea.add(sid, nid, component, scale, comment=comment)
        return darea

    def add_tload1(self, sid: int, excite_id: int, tabled_id: int, delay: int=0,
                   load_type: str='LOAD', us0: float=0.0, vs0: foat=0.0,
                   comment: str='') -> int:
        """
        Creates a TLOAD1 card, which defines a load based on a table

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        tabled_id : int
            TABLEDi id that defines F(t) for all degrees of freedom in
            EXCITEID entry
            float : MSC not supported
        delay : int/float; default=0
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        Type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        us0 : float; default=0.
            Factor for initial displacements of the enforced degrees-of-freedom
            MSC only
        vs0 : float; default=0.
            Factor for initial velocities of the enforced degrees-of-freedom
            MSC only
        comment : str; default=''
            a comment for the card

        """
        load = self.tload1.add(sid, excite_id, tabled_id, delay=delay,
                               load_type=load_type, us0=us0, vs0=vs0,
                               comment=comment)
        return load

    def add_tload2(self, sid: int, excite_id: int, delay: int=0,
                   load_type: str='LOAD',
                   T1: float=0., T2: Optional[float]=None,
                   frequency: float=0., phase: float=0.,
                   c: float=0., b: float=0.,
                   us0: float=0., vs0: float=0.,
                   comment: str='') -> int:
        """
        Creates a TLOAD2 card, which defines a exponential time load

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        load_type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        T1 : float; default=0.
            time constant (t1 > 0.0)
            times below this are ignored
        T2 : float; default=None
            time constant (t2 > t1)
            times above this are ignored
        frequency : float; default=0.
            Frequency in cycles per unit time.
        phase : float; default=0.
            Phase angle in degrees.
        c : float; default=0.
            Exponential coefficient.
        b : float; default=0.
            Growth coefficient.
        us0 : float; default=0.
            Factor for initial displacements of the enforced degrees-of-freedom
            MSC only
        vs0 : float; default=0.
            Factor for initial velocities of the enforced degrees-of-freedom
            MSC only
        comment : str; default=''
            a comment for the card

        """
        load = self.tload2.add(
            sid, excite_id, delay=delay, load_type=load_type, T1=T1, T2=T2,
            frequency=frequency, phase=phase, c=c, b=b,
            us0=us0, vs0=vs0, comment=comment)
        return load

    def add_rload1(self, sid: int, excite_id: int,
                   delay: int | float=0,
                   dphase: int | float=0,
                   tc: int | float=0,
                   td: int | float=0,
                   load_type='LOAD', comment: str='') -> int:
        """
        Creates an RLOAD1 card, which defines a frequency-dependent load
        based on TABLEDs.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        dphase : int/float; default=None
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        tc : int/float; default=0
            TABLEDi id that defines C(f) for all degrees of freedom in
            EXCITEID entry
        td : int/float; default=0
            TABLEDi id that defines D(f) for all degrees of freedom in
            EXCITEID entry
        Type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        comment : str; default=''
            a comment for the card

        """
        load = self.rload1.add(
            sid, excite_id, delay=delay, dphase=dphase, tc=tc, td=td,
            load_type=load_type, comment=comment)
        return load

    def add_rload2(self, sid: int, excite_id: int,
                   delay: int | float=0,
                   dphase: int | float=0,
                   tb: int | float=0,
                   tphi: int | float=0,
                   load_type: str='LOAD', comment: str='') -> int:
        """
        Creates an RLOAD2 card, which defines a frequency-dependent load
        based on TABLEDs.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        dphase : int/float; default=None
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        tb : int/float; default=0
            TABLEDi id that defines B(f) for all degrees of freedom in
            EXCITEID entry
        tphi : int/float; default=0
            TABLEDi id that defines phi(f) for all degrees of freedom in
            EXCITEID entry
        load_type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        comment : str; default=''
            a comment for the card

        """
        load = self.rload2.add(
            sid, excite_id, delay=delay, dphase=dphase, tb=tb, tphi=tphi,
            load_type=load_type, comment=comment)
        return load

    def add_rforce(self, sid: int, nid: int, scale: float, r123: list[float],
                   cid: int=0, method: int=1, racc: float=0.,
                   main_bulk: int=0, idrf: int=0, comment: str='') -> int:
        """Creates an RFORCE card"""
        load = self.rforce.add(sid, nid, scale, r123, cid=cid, method=method, racc=racc,
                               main_bulk=main_bulk, idrf=idrf, comment=comment)
        return load

    def add_rforce1(self, sid: int, nid: int, scale: float,
                    group_id: int, cid: int=0, r123: Optional[list[float]]=None,
                    racc: float=0., main_bulk: int=0, method: int=2,
                    comment: str='') -> int:
        """
        Creates an RFORCE1 card

        Parameters
        ----------
        sid : int
            load set id
        nid : int
            grid point through which the rotation vector acts
        scale : float
            scale factor of the angular velocity in revolutions/time
        r123 : list[float, float, float] / (3, ) float ndarray
            rectangular components of the rotation vector R that passes
            through point G
        racc : int; default=0.0
            ???
        main_bulk : int; default=0
            Indicates whether the CID coordinate system is defined in the main
            Bulk Data Section (MB = -1) or the partitioned superelement Bulk
            Data Section (MB = 0). Coordinate systems referenced in the main
            Bulk Data Section are considered stationary with respect to the
            assembly basic coordinate system.
        group_id : int
            Group identification number. The GROUP entry referenced in the
            GROUPID field selects the grid points to which the load is applied.
        cid : int; default=0
            Coordinate system defining the components of the rotation vector.
        method : int; default=2
            Method used to compute centrifugal forces due to angular velocity.
        comment : str; default=''
            a comment for the card

        """
        load = self.rforce1.add(
            sid, nid, scale, group_id, cid=cid, r123=r123, racc=racc,
            main_bulk=main_bulk, method=method, comment=comment)
        return load

    def add_randps(self, sid, j, k, x=0., y=0., tid=0, comment: str='') -> int:
        """
        Creates a RANDPS card

        Parameters
        ----------
        sid : int
            random analysis set id
            defined by RANDOM in the case control deck
        j : int
            Subcase id of the excited load set
        k : int
            Subcase id of the applied load set
            k > j
        x / y : float; default=0.0
            Components of the complex number
        tid : int; default=0
            TABRNDi id that defines G(F)
        comment : str; default=''
            a comment for the card

        """
        randps = self.randps.add(sid, j, k, x=x, y=y, tid=tid, comment=comment)
        return randps

    def add_randt1(self, sid, n, t0, tmax, comment: str='') -> int:
        """
        Creates a RANDT1 card

        Parameters
        ----------
        sid : int
            random analysis set id
            defined by RANDOM in the case control deck
        n : int
            ???
        t0 : int
            ???
        tmax : float
            ???
        comment : str; default=''
            a comment for the card

        """
        dload = self.randt1.add(sid, n, t0, tmax, comment=comment)
        return dload

    def add_acsrce(self, sid: int, excite_id: int, rho: float, b: float,
                   delay: int | float=0,
                   dphase: int | float=0,
                   power: int | float=0,
                   comment: str='') -> int:
        """
        Creates an ACSRCE card

        Parameters
        ----------
        sid : int
            load set id number (referenced by DLOAD)
        excite_id : int
            Identification number of a DAREA or SLOAD entry that lists
            each degree of freedom to apply the excitation and the
            corresponding scale factor, A, for the excitation
        rho : float
            Density of the fluid
        b : float
            Bulk modulus of the fluid
        delay : int; default=0
            Time delay, τ.
        dphase : int / float; default=0
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        power : int; default=0
            Power as a function of frequency, P(f).
            float : value of P(f) used over all frequencies for all
                    degrees of freedom in EXCITEID entry.
            int : TABLEDi entry that defines P(f) for all degrees of
                  freedom in EXCITEID entry.
        comment : str; default=''
            a comment for the card

        """
        load = self.acsrce.add(
            sid, excite_id, rho, b,
            delay=delay, dphase=dphase, power=power,
            comment=comment)
        return load

    def add_loadcyn(self, sid, scale, segment_id, scales, load_ids,
                    segment_type=None, comment: str='') -> LOADCYN:
        """Creates a LOADCYN card"""
        load = LOADCYN(sid, scale, segment_id, scales, load_ids,
                       segment_type=segment_type, comment=comment)
        self._add_methods.add_load_object(load)
        return load

    def add_loadcyh(self, sid, scale, hid, htype, scales, load_ids,
                    comment: str='') -> LOADCYH:
        """Creates a LOADCYH card"""
        load = LOADCYH(sid, scale, hid, htype, scales, load_ids, comment=comment)
        self._add_methods.add_load_object(load)
        return load

    def add_force(self, sid: int, node: int, mag: float, xyz: np.ndarray, cid: int=0,
                  comment: str='') -> int:
        """
        Creates a FORCE card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        xyz : (3, ) float ndarray
            the load direction in the cid frame
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card

        """
        load = self.force.add(sid, node, mag, xyz, cid=cid, comment=comment)
        return load

    def add_force1(self, sid: int, node: int, mag: float,
                   g1: int, g2: int, comment: str='') -> int:
        """
        Creates a FORCE1 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        n1 / n2 : int / int
            defines the load direction
            n = n2 - n1
        comment : str; default=''
            a comment for the card

        """
        load = self.force1.add(sid, node, mag, g1, g2, comment=comment)
        return load

    def add_force2(self, sid: int, node: int, mag: float,
                   g1: int, g2: int, g3: int, g4: int, comment: str='') -> int:
        """
        Creates a FORCE2 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        g1 / g2 / g3 / g4 : int / int / int / int
            defines the load direction
            n = (g2 - g1) x (g4 - g3)
        comment : str; default=''
            a comment for the card

        """
        load = self.force2.add(sid, node, mag, g1, g2, g3, g4, comment=comment)
        return load

    def add_moment(self, sid, node, mag, xyz, cid=0, comment: str='') -> int:
        """
        Creates a MOMENT card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        cid : int; default=0
            the coordinate system for the load
        xyz : (3, ) float ndarray; default
            the load direction in the cid frame
        comment : str; default=''
            a comment for the card

        """
        load = self.moment.add(sid, node, mag, xyz, cid=cid, comment=comment)
        return load

    def add_moment1(self, sid, node, mag, g1, g2, comment: str='') -> int:
        """
        Creates a MOMENT1 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        n1 / n2 : int / int
            defines the load direction
            n = n2 - n1
        comment : str; default=''
            a comment for the card

        """
        load = self.moment1.add(sid, node, mag, g1, g2, comment=comment)
        return load

    def add_moment2(self, sid, node, mag, g1, g2, g3, g4, comment: str='') -> int:
        """
        Creates a MOMENT2 card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        g1 / g2 / g3 / g4 : int / int / int / int
            defines the load direction
            n = (g2 - g1) x (g4 - g3)
        comment : str; default=''
            a comment for the card

        """
        load = self.moment2.add(sid, node, mag, g1, g2, g3, g4, comment=comment)
        return load

    def add_deform(self, sid: int, eid: int, deformation: float, comment: str='') -> int:
        """
        Creates an DEFORM card, which defines applied deformation on
        a 1D element.  Links to the DEFORM card in the case control
        deck.

        Parameters
        ----------
        sid : int
            load id
        eid : int
            CTUBE/CROD/CONROD/CBAR/CBEAM element id
        deformation : float
            the applied deformation
        comment : str; default=''
            a comment for the card

        """
        load = self.deform.add(sid, eid, deformation, comment=comment)
        return load

    def add_accel(self, sid: int, N: list[float], direction: str,
                  locs: list[float], vals: list[float], cid: int=0,
                  comment: str='') -> int:
        """
        Creates an ACCEL card

        Parameters
        ----------
        sid : int
            load id
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        direction : str
            Component direction of acceleration variation
            {X, Y, Z}
        locs : list[float]
            Location along direction DIR in coordinate system CID for
            specification of a load scale factor.
        vals : list[float]
            The load scale factor associated with location LOCi
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card

        """
        load = self.accel.add(sid, N, direction, locs, vals, cid=cid,
                              comment=comment)
        return load

    def add_accel1(self, sid: int, scale: float,
                   N: list[float], nodes: list[int],
                   cid: int=0, comment: str='') -> int:
        """
        Creates an ACCEL1 card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            scale factor for load
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        nodes : list[int]
            the nodes to apply acceleration to
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card

        """
        load = self.accel1.add(sid, scale, N, nodes, cid=cid, comment=comment)
        return load

    def add_grav(self, sid: int, scale: float, N: list[float],
                 cid: int=0, mb: int=0, comment: str='') -> int:
        """
        Creates an GRAV card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            scale factor for load
        N : (3, ) float ndarray
            the acceleration vector in the cid frame
        cid : int; default=0
            the coordinate system for the load
        mb : int; default=0
            ???
        comment : str; default=''
            a comment for the card

        """
        grav = self.grav.add(sid, scale, N, cid=cid, mb=mb, comment=comment)
        return grav

    def add_pload(self, sid: int, pressure: float, nodes: list[int],
                  comment: str='') -> int:
        """
        Creates a PLOAD card, which defines a uniform pressure load on a
        shell/solid face or arbitrarily defined quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply
        nodes : list[int]
            The nodes that are used to define the normal are defined
            using the same method as the CTRIA3/CQUAD4 normal.
            n = 3 or 4
        comment : str; default=''
            a comment for the card

        """
        load = self.pload.add(sid, pressure, nodes, comment=comment)
        return load

    def add_pload1(self, sid: int, eid: int, load_type: str, scale: float,
                   x1: float, p1: float,
                   x2: Optional[float]=None, p2: Optional[float]=None,
                   comment: str='') -> int:
        """
        Creates a PLOAD1 card, which may be applied to a CBAR/CBEAM

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element to apply the load to
        load_type : str
            type of load that's applied
            valid_types = {FX, FY, FZ, FXE, FYE, FZE,
                           MX, MY, MZ, MXE, MYE, MZE}
        scale : str
            Determines scale factor for X1, X2.
            {LE, FR, LEPR, FRPR}
        x1 / x2 : float / float
            the starting/end position for the load application
            the default for x2 is x1
        p1 / p2 : float / float
            the magnitude of the load at x1 and x2
            the default for p2 is p1
        comment : str; default=''
            a comment for the card

        Point Load       : x1 == x2
        Distributed Load : x1 != x2

        """
        load = self.pload1.add(sid, eid, load_type, scale, x1, p1, x2=x2, p2=p2,
                               comment=comment)
        return load

    def add_pload2(self, sid: int, pressure: float, eids: list[int],
                   comment: str='') -> int:
        """
        Creates a PLOAD2 card, which defines an applied load normal
        to the quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply to the elements
        eids : list[int]
            the elements to apply pressure to
            n < 6 or a continouus monotonic list of elements (e.g., [1, 2, ..., 1000])
        comment : str; default=''
            a comment for the card

        """
        load = self.pload2.add(sid, pressure, eids, comment=comment)
        return load

    def add_pload4(self,
                   sid: int,
                   eids: list[int],
                   pressures: list[float],
                   g1=-1, g34=-1, cid: int=0,
                   nvector=None,
                   surf_or_line: str='SURF', line_load_dir: str='NORM',
                   comment: str='') -> int:
        """
        Creates a PLOAD4 card

        Parameters
        ----------
        sid : int
            the load id
        eids : list[int, ...]
            shells : the range of element ids; must be sequential
            solids : must be length 1
        pressures : list[float, float, float, float]
            tri : must be length 4 (the last value should be the same as the 0th value)
            quad : must be length 4
        g1 : int/None
            only used for solid elements
        g34 : int / None
            only used for solid elements
        cid : int; default=0
            the coordinate system for ???
        nvector : (3, ) float ndarray
           blank : load acts normal to the face
           the local pressure vector
        surf_or_line : str; default='SURF'
           SURF : surface load
           LINE : line load (only defined for QUADR, TRIAR)
           not supported
        line_load_dir : str; default='NORM'
           direction of the line load (see surf_or_line); {X, Y, Z, TANG, NORM}
           not supported
        comment : str; default=''
            a comment for the card

        TODO: fix the way "pressures" works

        """
        load = self.pload4.add(sid, eids, pressures, g1=g1, g34=g34, cid=cid,
                               nvector=nvector, surf_or_line=surf_or_line,
                               line_load_dir=line_load_dir, comment=comment)
        return load

    def add_ploadx1(self, sid: int, eid: int, pa: float,
                    nids: list[int], pb: Optional[float]=None, theta: float=0.,
                    comment: str='') -> int:
        """
        Creates a PLOADX1 card, which defines surface traction for
        axisymmetric elements.

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element id (CQUADX, CTRIAX, or CTRIAX6)
        nids : list[int, int]
            Corner grid points.
            GA and GB are any two adjacent corner grid points of the element
        pa / pb : float / None
            Surface traction at grid point GA or GB
            pb : default is None -> pa
        theta : float; default=0.0
            Angle between surface traction and inward normal to the line
            segment.
        comment : str; default=''
            a comment for the card

        """
        load = self.ploadx1.add(sid, eid, pa, nids, pb=pb, theta=theta,
                                comment=comment)
        return load

    def add_spc(self, conid :int, nodes: list[int], components: list[str],
                enforced: list[float], comment: str='') -> int:
        """
        Creates an SPC card, which defines the degree of freedoms to be
        constrained

        Parameters
        ----------
        conid : int
            constraint id
        nodes : list[int]
            GRID/SPOINT ids
        components : list[int]
            the degree of freedoms to constrain (e.g., 1, 123)
        enforced : list[float]
            the constrained value for the given node (typically 0.0)
        comment : str; default=''
            a comment for the card

        Notes
        -----
        len(nodes) == len(components) == len(enforced)

        .. warning:: non-zero enforced deflection requires an SPCD as well

        """
        spc = self.spc.add(conid, nodes, components, enforced, comment=comment)
        return spc

    def add_spc1(self, spc_id: int, components: int, nodes: list[int],
                 comment: str='') -> int:
        """
        Creates an SPC1 card, which defines the degree of freedoms to be
        constrained to a value of 0.0

        Parameters
        ----------
        spc_id : int
            constraint id
        components : int
            the degree of freedoms to constrain (e.g., 1, 123)
        nodes : list[int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        spc = self.spc1.add(spc_id, components, nodes, comment=comment)
        return spc

    def add_spcd(self, spc_id: int, nodes: list[int],
                 components: list[int], enforced: list[float],
                 comment: str='') -> int:
        """
        Creates an SPCD card, which defines the degree of freedoms to be
        set during enforced motion

        Parameters
        ----------
        spc_id : int
            constraint id
        nodes : list[int]
            GRID/SPOINT ids
        components : list[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        enforced : list[float]
            the constrained value for the given node (typically 0.0)
        comment : str; default=''
            a comment for the card

        Notes
        -----
        len(nodes) == len(components) == len(enforced)

        .. warning:: Non-zero enforced deflection requires an SPC/SPC1 as well.
                     Yes, you really want to constrain the deflection to 0.0
                     with an SPC1 card and then reset the deflection using an
                     SPCD card.

        """
        spc = self.spcd.add(spc_id, nodes, components, enforced, comment=comment)
        return spc

    def add_spcadd(self, conid, sets, comment: str='') -> int:
        """Creates a SPCADD card"""
        spcadd = self.spcadd.add(conid, sets, comment=comment)
        return spcadd

    def add_spcax(self, conid, ringax, hid, component, enforced,
                  comment: str='') -> int:
        """Creates an SPCAX card"""
        spcax = self.spcax.add(conid, ringax, hid, component, enforced, comment=comment)
        return spcax

    #def add_gmspc(self, conid, component, entity, entity_id, comment: str='') -> GMSPC:
        #"""Creates a GMSPC card"""
        #spc = GMSPC(conid, component, entity, entity_id, comment=comment)
        #self._add_methods.add_constraint_spc_object(spc)
        #return spc

    def add_mpc(self, mpc_id: int,
                nodes: list[int], components: list[str], coefficients: list[float],
                comment: str='') -> int:
        """
        Creates an MPC card

        Parameters
        ----------
        conid : int
            Case Control MPC id
        nodes : list[int]
            GRID/SPOINT ids
        components : list[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        coefficients : list[float]
            the scaling coefficients

        """
        mpc = self.mpc.add(mpc_id, nodes, components, coefficients, comment=comment)
        return mpc

    def add_mpcadd(self, mpc_id, sets, comment: str='') -> int:
        """Creates an MPCADD card"""
        mpcadd = self.mpcadd.add(mpc_id, sets, comment=comment)
        return mpcadd

    def add_suport(self, nodes: list[int], components: list[int],
                   comment: str='') -> int:
        """
        Creates a SUPORT card, which defines free-body reaction points.
        This is always active.

        Parameters
        ----------
        nodes : list[int]
            the nodes to release
        components : list[int]
            components to support at each node (1, 13)
        comment : str; default=''
            a comment for the card

        """
        self.suport.add_set(nodes, components, comment=comment)
        #suport = SUPORT(nodes, Cs, comment=comment)
        #self._add_methods.add_suport_object(suport)
        #return suport
        return self.suport

    def add_suport1(self, suport_id, nodes: list[int], components: list[int],
                    comment: str='') -> int:
        """
        Creates a SUPORT card, which defines free-body reaction points.

        Parameters
        ----------
        conid : int
            Case Control SUPORT id
        nodes : list[int]
            the nodes to release
        components : list[int]
            components to support at each node (1, 13)
        comment : str; default=''
            a comment for the card

        """
        suport1 = self.suport.add_set1(suport_id, nodes, components, comment=comment)
        return suport1

    def add_eigr(self, sid, method='LAN', f1=None, f2=None, ne=None, nd=None,
                 norm='MASS', G=None, C=None, comment: str='') -> EIGR:
        """
        Adds a EIGR card

        Parameters
        ----------
        sid : int
            method id
        method : str; default='LAN'
            eigenvalue method
            recommended: {LAN, AHOU}
            obsolete : {INV, SINV, GIV, MGIV, HOU, MHOU, AGIV}
        f1 / f2 : float; default=None
            lower/upper bound eigenvalue
        f2 : float; default=None
            upper bound eigenvalue
        ne : int; default=None
            estimate of number of roots (used for INV)
        nd : int; default=None
            desired number of roots
        msglvl : int; default=0
            debug level; 0-4
        maxset : int; default=None
            Number of vectors in block or set
        shfscl : float; default=None
            estimate of first flexible mode natural frequency
        norm : str; default=None
            {MAX, MASS, AF, POINT}
            default=MASS (NX)
        G : int; default=None
            node id for normalization; only for POINT
        C : int; default=None
            component for normalization (1-6); only for POINT
        comment : str; default=''
            a comment for the card

        """
        method = EIGR(sid, method=method, f1=f1, f2=f2, ne=ne, nd=nd,
                      norm=norm, G=G, C=C, comment=comment)
        self._add_methods.add_method_object(method)
        return method

    def add_eigrl(self, sid, v1=None, v2=None, nd=None, msglvl=0, maxset=None, shfscl=None,
                  norm=None, options=None, values=None, comment: str='') -> EIGRL:
        """
        Adds an EIGRL card

        Parameters
        ----------
        sid : int
            method id
        v1 : float; default=None
            lower bound eigenvalue
        v2 : float; default=None
            upper bound eigenvalue
        nd : int
            number of roots
        msglvl : int; default=0
            debug level; 0-4
        maxset : int; default=None
            Number of vectors in block or set
        shfscl : float; default=None
            estimate of first flexible mode natural frequency
        norm : str; default=None
            {MAX, MASS}
        options : ???; default=None -> []
            ???
        values : ???; default=None -> []
            ???
        comment : str; default=''
            a comment for the card

        """
        method = EIGRL(sid, v1, v2, nd, msglvl, maxset, shfscl, norm, options,
                       values, comment=comment)
        self._add_methods.add_method_object(method)
        return method

    def add_eigb(self, sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                 comment: str='') -> EIGB:
        """Creates an EIGB card"""
        method = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                      comment=comment)
        self._add_methods.add_method_object(method)
        return method

    def add_eigc(self, sid, method, grid, component, epsilon, neigenvalues,
                 norm='MAX', # common
                 mblkszs=None, iblkszs=None, ksteps=None, NJIs=None,
                 alphaAjs=None, omegaAjs=None,
                 alphaBjs=None, omegaBjs=None,
                 LJs=None, NEJs=None, NDJs=None,
                 shift_r1=None, shift_i1=None, isrr_flag=None,
                 nd1=None, comment: str='') -> EIGC:
        """Creates an EIGC card"""
        method = EIGC(sid, method, grid, component, epsilon, neigenvalues, norm=norm,
                      mblkszs=mblkszs, iblkszs=iblkszs, ksteps=ksteps, NJIs=NJIs,

                      alphaAjs=alphaAjs, omegaAjs=omegaAjs,
                      alphaBjs=alphaBjs, omegaBjs=omegaBjs,

                      LJs=LJs, NEJs=NEJs, NDJs=NDJs,
                      shift_r1=shift_r1, shift_i1=shift_i1, isrr_flag=isrr_flag,
                      nd1=nd1, comment=comment)
        self._add_methods.add_cmethod_object(method)
        return method

    def add_eigp(self, sid, alpha1, omega1, m1, alpha2, omega2, m2,
                 comment: str='') -> EIGP:
        """Creates an EIGP card"""
        method = EIGP(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=comment)
        self._add_methods.add_cmethod_object(method)
        return method

    def add_set1(self, sid: int, ids: list[int], is_skin: bool=False,
                 comment: str='') -> int:
        """
        Creates a SET1 card, which defines a list of structural grid
        points or element identification numbers.

        Parameters
        ----------
        sid : int
            set id
        ids : list[int, str]
            AECOMP, SPLINEx, PANEL : all grid points must exist
            XYOUTPUT : missing grid points are ignored
            The only valid string is THRU
            ``ids = [1, 3, 5, THRU, 10]``
        is_skin : bool; default=False
            if is_skin is used; ids must be empty
        comment : str; default=''
            a comment for the card

        """
        set_obj = self.set1.add(sid, ids, is_skin=is_skin, comment=comment)
        return set_obj

    def add_set2(self, sid: int, macro: int,
                 sp1: float, sp2: float,
                 ch1: float, ch2: float,
                 zmax: float=0.0, zmin: float=0.0,
                 comment: str='') -> int:
        """
        Creates a SET2 card, which defines a list of structural grid
        points in terms of aerodynamic macro elements.

        Parameters
        ----------
        sid : int
            set id
        macro : int
            the aerodynamic macro element id
        sp1 / sp2 : float
            lower/higher span division point defining the prism containing the set
        ch1 / ch2 : float
            lower/higher chord division point defining the prism containing the set
        zmax / zmin : float; default=0.0/0.0
            z-coordinate of top/bottom of the prism containing the set
            a zero value implies a value of infinity
        comment : str; default=''
            a comment for the card

        """
        set_index = self.set2.add(sid, macro, sp1, sp2, ch1, ch2,
                                  zmax=zmax, zmin=zmin, comment=comment)
        return set_index

    def add_set3(self, sid: int, desc: str, ids: list[int],
                 comment: str='') -> int:
        """Creates a SET3 card"""
        set_obj = self.set3.add(sid, desc, ids, comment=comment)
        return set_obj

    def add_aset(self, ids: list[int],
                 components: str | list[str],
                 comment: str='') -> int:
        """
        Creates an ASET/ASET1 card, which defines the degree of freedoms
        that will be retained during an ASET modal reduction.

        Parameters
        ----------
        ids : list[int]
            the GRID/SPOINT ids
        components : list[str]; str
            the degree of freedoms to be retained (e.g., '1', '123')
            if a list is passed in, a ASET is made
            if a str is passed in, a ASET1 is made
        comment : str; default=''
            a comment for the card

        Notes
        -----
        the length of components and ids must be the same
        """
        ifile = 0
        if isinstance(components, list):
            aset = self.aset.add_set(ids, components, ifile=ifile, comment=comment)
        else:
            assert isinstance(components, (str, integer_types)), components
            aset = self.aset.add_set1(ids, components, ifile=ifile, comment=comment)
        return aset

    def add_aset1(self, ids: list[int],
                  components: str | list[str],
                  comment: str='') -> int:
        """.. .. seealso:: ``add_aset``"""
        return self.add_aset(ids, components, comment=comment)

    def add_bset(self, ids: list[int],
                 components: str | list[str],
                 comment: str='') -> int:
        """
        Creates an BSET/BSET1 card, which defines the degree of freedoms
        that will be fixed during a generalized dynamic reduction or
        component model synthesis calculation.

        Parameters
        ----------
        ids : list[int]
            the GRID/SPOINT ids
        components : list[str]; str
            the degree of freedoms to be fixed (e.g., '1', '123')
            if a list is passed in, a ASET is made
            if a str is passed in, a ASET1 is made
        comment : str; default=''
            a comment for the card

        Notes
        -----
        the length of components and ids must be the same
        """
        ifile = 0
        if isinstance(components, list):
            bset = self.bset.add_set(ids, components, ifile=ifile, comment=comment)
        else:
            assert isinstance(components, (str, integer_types)), components
            bset = self.bset.add_set1(ids, components, ifile=ifile, comment=comment)
        return bset

    def add_bset1(self, ids: list[int],
                  components: str | list[str],
                  comment: str='') -> int:
        """.. .. seealso:: ``add_bset``"""
        return self.add_bset(ids, components, comment=comment)

    def add_cset(self, ids: list[int],
                 components: str | list[str],
                 comment: str='') -> int:
        """
        Creates an CSET/CSET1 card, which defines the degree of freedoms
        that will be free during a generalized dynamic reduction or
        component model synthesis calculation.

        Parameters
        ----------
        ids : list[int]
            the GRID/SPOINT ids
        components : list[str]; str
            the degree of freedoms to be free (e.g., '1', '123')
            if a list is passed in, a CSET is made
            if a str is passed in, a CSET1 is made
        comment : str; default=''
            a comment for the card

        Notes
        -----
        the length of components and ids must be the same

        """
        ifile = 0
        if isinstance(components, list):
            cset = self.cset.add_set(ids, components, ifile=ifile, comment=comment)
        else:
            assert isinstance(components, (str, integer_types)), components
            cset = self.cset.add_set1(ids, components, ifile=ifile, comment=comment)
        return cset

    def add_cset1(self, ids: list[int],
                  components: str | list[str],
                  comment: str='') -> int:
        """.. seealso:: ``add_cset``"""
        return self.add_cset(ids, components, comment=comment)

    def add_omit1(self, ids: list[int],
                  components: str | list[str],
                  comment: str='') -> int:
        """.. seealso:: ``add_omit``"""
        return self.add_omit(ids, components, comment=comment)

    def add_omit(self, ids: list[int],
                  components: str | list[str],
                  comment: str='') -> int:
        """
        Creates an OMIT1 card, which defines the degree of freedoms that
        will be excluded (o-set) from the analysis set (a-set).

        Parameters
        ----------
        ids : list[int]
            the GRID/SPOINT ids
        components : list[str]; str
            the degree of freedoms to be retained (e.g., '1', '123')
            if a list is passed in, a OMIT is made
            if a str is passed in, a OMIT1 is made
        comment : str; default=''
            a comment for the card

        """
        ifile = 0
        if isinstance(components, list):
            omit = self.omit.add_set(ids, components, ifile=ifile, comment=comment)
        else:
            assert isinstance(components, (str, integer_types)), components
            omit = self.omit.add_set1(ids, components, ifile=ifile, comment=comment)
        return omit

    def add_qset(self, ids: list[int],
                 components: str | list[str],
                 comment: str='') -> int:
        """
        Creates a QSET/QSET1 card, which defines generalized degrees of
        freedom (q-set) to be used for dynamic reduction or component
        mode synthesis.

        Parameters
        ----------
        ids : list[int]
            the GRID/SPOINT ids
        components : list[str]; str
            the degree of freedoms to be created (e.g., '1', '123')
            if a list is passed in, a QSET is made
            if a str is passed in, a QSET1 is made
        comment : str; default=''
            a comment for the card

        """
        ifile = 0
        if isinstance(components, list):
            qset = self.qset.add_set(ids, components, ifile=ifile, comment=comment)
        else:
            assert isinstance(components, (str, integer_types)), components
            qset = self.qset.add_set1(ids, components, ifile=ifile, comment=comment)
        return qset

    def add_qset1(self, ids: list[int],
                  components: str | list[str],
                  comment: str='') -> int:
        """.. seealso:: ``add_qset``"""
        return self.add_qset(ids, components, comment=comment)

    def add_uset(self, name: str, ids: list[int], components: list[int],
                 comment: str='') -> int:
        """
        Creates a USET card, which defines a degrees-of-freedom set.

        Parameters
        ----------
        name : str
            SNAME Set name. (One to four characters or the word 'ZERO'
            followed by the set name.)
        ids : list[int]
            the GRID/SPOINT ids
        components : list[int]
            the degree of freedoms (e.g., 1, 123)
            if a list is passed in, a USET is made
            if an int is passed in, a USET1 is made
        comment : str; default=''
            a comment for the card

        """
        ifile = 0
        if isinstance(components, list):
            uset = self.uset.add_set(name, ids, components, ifile=ifile, comment=comment)
        else:
            assert isinstance(components, (str, integer_types)), components
            uset = self.uset.add_set1(name, ids, components, ifile=ifile, comment=comment)
        assert uset is not None
        return uset

    def add_uset1(self, name, ids, components, comment: str='') -> int:
        """.. seealso:: ``add_uset``"""
        return self.add_uset(name, ids, components, comment=comment)

    def add_dtable(self, default_values: dict[str, float],
                   comment: str='') -> DTABLE:
        """
        Creates a DTABLE card

        Parameters
        ----------
        default_values : dict
            key : str
                the parameter name
            value : float
                the value
        comment : str; default=''
            a comment for the card

        """
        dtable = DTABLE(default_values, comment=comment)
        self._add_methods.add_dtable_object(dtable)
        return dtable

    def add_tabled1(self, tid: int,
                    x: np.ndarray, y: np.ndarray,
                    xaxis: str='LINEAR', yaxis: str='LINEAR', extrap: int=0,
                    comment: str='') -> TABLED1:
        """
        Creates a TABLED1, which is a dynamic load card that is applied
        by the DAREA card

        Parameters
        ----------
        tid : int
            table id
        x : list[float]
            nvalues
        y : list[float]
            nvalues
        xaxis : str
            LINEAR, LOG
        yaxis : str
            LINEAR, LOG
        extrap : int; default=0
            Extrapolation method:
                0 : linear
                1 : constant
            .. note:: this is NX specific
        comment : str; default=''
            a comment for the card

        """
        table = TABLED1(tid, x, y, xaxis=xaxis, yaxis=yaxis,
                        extrap=extrap, comment=comment)
        self._add_methods.add_tabled_object(table)
        return table

    def add_tabled2(self, tid: int, x1: float,
                    x: np.ndarray, y: np.ndarray,
                    extrap: int=0, comment: str='') -> TABLED2:
        """Creates a TABLED2 card"""
        table = TABLED2(tid, x1, x, y, extrap=extrap, comment=comment)
        self._add_methods.add_tabled_object(table)
        return table

    def add_tabled3(self, tid: int, x1: float, x2: float,
                    x: np.ndarray, y: np.ndarray,
                    extrap: int=0, comment: str='') -> TABLED3:
        """Creates a TABLED3 card"""
        table = TABLED3(tid, x1, x2, x, y, extrap=extrap, comment=comment)
        self._add_methods.add_tabled_object(table)
        return table

    def add_tabled4(self, tid: int,
                    x1: float, x2: float, x3: float, x4: float,
                    a: list[float], comment: str='') -> TABLED4:
        """Creates a TABLED4 card"""
        table = TABLED4(tid, x1, x2, x3, x4, a, comment=comment)
        self._add_methods.add_tabled_object(table)
        return table

    def add_tablem1(self, tid: int, x: np.ndarray, y: np.ndarray,
                    xaxis: str='LINEAR', yaxis: str='LINEAR',
                    extrap: int=0, comment: str='') -> TABLEM1:
        """Creates a TABLEM1 card"""
        table = TABLEM1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)
        self._add_methods.add_tablem_object(table)
        return table

    def add_tablem2(self, tid: int, x1: float,
                    x: np.ndarray, y: np.ndarray,
                    extrap: int=0, comment: str='') -> TABLEM2:
        """Creates a TABLEM2 card"""
        table = TABLEM2(tid, x1, x, y, extrap=extrap, comment=comment)
        self._add_methods.add_tablem_object(table)
        return table

    def add_tablem3(self, tid: int, x1: float, x2: float,
                    x: np.ndarray, y: np.ndarray,
                    extrap: int=0, comment: str='') -> TABLEM3:
        """Creates a TABLEM3 card"""
        table = TABLEM3(tid, x1, x2, x, y, extrap=extrap, comment=comment)
        self._add_methods.add_tablem_object(table)
        return table

    def add_tablem4(self, tid: int,
                    x1: float, x2: float, x3: float, x4: float,
                    a: list[float], comment: str='') -> TABLEM4:
        """Creates a TABLEM4 card"""
        table = TABLEM4(tid, x1, x2, x3, x4, a, comment=comment)
        self._add_methods.add_tablem_object(table)
        return table

    def add_tables1(self, tid: int, x: np.ndarray, y: np.ndarray,
                    Type: int=1, comment: str='') -> TABLES1:
        """
        Adds a TABLES1 card, which defines a stress dependent material

        Parameters
        ----------
        tid : int
            Table ID
        Type : int; default=1
            Type of stress-strain curve (1 or 2)
            1 - Cauchy (true) stress vs. total true strain
            2 - Cauchy (true) stress vs. plastic true strain (MSC only)
        x, y : list[float]
            table values
        comment : str; default=''
            a comment for the card
        """
        table = TABLES1(tid, x, y, Type=Type, comment=comment)
        self._add_methods.add_table_object(table)
        return table

    def add_tablest(self, tid, x, y, comment: str='') -> TABLEST:
        """Creates an TABLEST card"""
        table = TABLEST(tid, x, y, comment=comment)
        self._add_methods.add_table_object(table)
        return table

    def add_tabrnd1(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment: str='') -> TABRND1:
        """Creates an TABRND1 card"""
        table = TABRND1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)
        self._add_methods.add_random_table_object(table)
        return table

    def add_tabrndg(self, tid: int, Type: int, LU: float, WG: float,
                    comment: str='') -> TABRNDG:
        """
        Creates a TABRNDG card

        Parameters
        ----------
        tid : int
            table id
        Type : int
           PSD type
           1 : von Karman
           2 : Dryden
        LU : float
            Scale of turbulence divided by velocity (units of time)
        WG : float
            Root-mean-square gust velocity
        comment : str; default=''
            a comment for the card
        """
        table = TABRNDG(tid, Type, LU, WG, comment=comment)
        self._add_methods.add_random_table_object(table)
        return table

    def add_tabdmp1(self, tid: int, x, y, Type: str='G', comment: str='') -> TABDMP1:
        """Creates a TABDMP1 card"""
        table = TABDMP1(tid, x, y, Type=Type, comment=comment)
        self._add_methods.add_table_sdamping_object(table)
        return table

    def add_tableht(self, tid: int, x, y, comment: str='') -> TABLEHT:
        table = TABLEHT(tid, x, y, comment=comment)
        self._add_methods.add_table_object(table)
        return table

    def add_tableh1(self, tid: int, x, y, comment: str='') -> TABLEH1:
        table = TABLEH1(tid, x, y, comment=comment)
        self._add_methods.add_table_object(table)
        return table

    def add_freq(self, sid, freqs, comment: str='') -> FREQ:
        """
        Creates a FREQ card

        Parameters
        ----------
        sid : int
            set id referenced by case control FREQUENCY
        freqs : list[float]
            the frequencies for a FREQx object
        comment : str; default=''
            a comment for the card

        """
        freq = FREQ(sid, freqs, comment=comment)
        self._add_methods.add_freq_object(freq)
        return freq

    def add_freq1(self, sid: int, f1: float, df: float, ndf: int=1,
                  comment: str='') -> FREQ1:
        """
        Creates a FREQ1 card

        Parameters
        ----------
        sid : int
            set id referenced by case control FREQUENCY
        f1 : float
            first frequency
        df : float
            frequency increment
        ndf : int; default=1
            number of frequency increments
        comment : str; default=''
            a comment for the card

        """
        freq = FREQ1(sid, f1, df, ndf=ndf, comment=comment)
        self._add_methods.add_freq_object(freq)
        return freq

    def add_freq2(self, sid: int, f1: float, f2: float, nf: int=1,
                  comment: str='') -> FREQ2:
        """
        Creates a FREQ2 card

        Parameters
        ----------
        sid : int
            set id referenced by case control FREQUENCY
        f1 : float
            first frequency
        f2 : float
            last frequency
        nf : int; default=1
            number of logorithmic intervals
        comment : str; default=''
            a comment for the card

        """
        freq = FREQ2(sid, f1, f2, nf=nf, comment=comment)
        self._add_methods.add_freq_object(freq)
        return freq

    def add_freq3(self, sid: int, f1: float, f2=None, Type: str='LINEAR',
                  nef: int=10, cluster: float=1.0,
                  comment: str='') -> FREQ3:
        """Creates a FREQ3 card"""
        freq = FREQ3(sid, f1, f2, Type, nef, cluster, comment)
        self._add_methods.add_freq_object(freq)
        return freq

    def add_freq4(self, sid: int, f1: float=0., f2: float=1e20,
                  fspread: float=0.1, nfm: int=3, comment: str='') -> FREQ4:
        """
        Creates a FREQ4 card

        Parameters
        ----------
        sid : int
            set id referenced by case control FREQUENCY
        f1 : float; default=0.0
            Lower bound of frequency range in cycles per unit time.
        f2 : float; default=1E20
            Upper bound of frequency range in cycles per unit time.
        nfm : int; default=3
            Number of evenly spaced frequencies per 'spread' mode.
        comment : str; default=''
            a comment for the card

        """
        freq = FREQ4(sid, f1, f2, fspread, nfm, comment=comment)
        self._add_methods.add_freq_object(freq)
        return freq

    def add_freq5(self, sid: int, fractions, f1: float=0., f2: float=1e20,
                  comment: str='') -> FREQ5:
        """
        Creates a FREQ5 card

        Parameters
        ----------
        sid : int
            set id referenced by case control FREQUENCY
        f1 : float; default=0.0
            Lower bound of frequency range in cycles per unit time.
        f2 : float; default=1e20
            Upper bound of frequency range in cycles per unit time.
        fractions : list[float]
            Fractions of the natural frequencies in the range F1 to F2.
        comment : str; default=''
            a comment for the card

        Notes
        -----
        FREQ5 is only valid in modal frequency-response solutions
        (SOLs 111, 146, and 200) and is ignored in direct frequency
        response solutions.

        """
        freq = FREQ5(sid, fractions, f1=f1, f2=f2, comment=comment)
        self._add_methods.add_freq_object(freq)
        return freq

    def add_tf(self, tf_id: int,
               nid0: int, component: int,
               b0: float, b1: float, b2: float,
               nids: list[int], components: list[int],
               a: list[tuple[float, float, float]],
               comment: str='') -> int:
        """Creates a TF card"""
        tf = self.tf.add(tf_id, nid0, component, b0, b1, b2,
                         nids, components, a, comment=comment)
        return tf

    def add_group(self, group_id: int, nodes, elements, properties, comment: str='') -> GROUP:
        self.log.warning('skipping GROUP')
        #group = GROUP(group_id, nodes, elements, properties, comment=comment)
        #self._add_methods.add_group_object(group)
        #return group

    # ------------------------------------------
    def add_modtrak(self, modtrak_id: int, low_range: int, high_range: int,
                    mt_filter: float, comment: str='') -> int:
        modtrak = self.modtrak.add(modtrak_id, low_range, high_range, mt_filter, comment=comment)
        return modtrak

    def add_mondsp1(self, name, label, axes, aecomp_name, xyz,
                          cp=0, cd=None, ind_dof='123', comment: str='') -> MONDSP1:
        """
        Creates a MONDSP1 card

        Parameters
        ----------
        name : str
            Character string of up to 8 characters identifying the
            monitor point
        label : str
            A string comprising no more than 56 characters
            that identifies and labels the monitor point.
        axes : str
            components {1,2,3,4,5,6}
        aecomp_name : str
            name of the AECOMP/AECOMPL entry
        xyz : list[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, Coord; default=0
           coordinate system of XYZ
        cd : int; default=None -> cp
            the coordinate system for load outputs
        ind_dof : str; default='123'
            the dofs to map
        comment : str; default=''
            a comment for the card

        Notes
        -----
        MSC specific card

        """
        mondsp1 = MONDSP1(name, label, axes, aecomp_name, xyz,
                          cp=cp, cd=cd, ind_dof=ind_dof, comment=comment)
        self._add_methods.add_monpnt_object(mondsp1)
        return mondsp1

    def add_monpnt1(self, name: str, label: str, axes: str, aecomp_name: str,
                    xyz: list[float],
                    cp: int=0,
                    cd: Optional[int]=None, comment: str='') -> int:
        """
        Creates a MONPNT1 card

        Parameters
        ----------
        name : str
            Character string of up to 8 characters identifying the
            monitor point
        label : str
            A string comprising no more than 56 characters
            that identifies and labels the monitor point.
        axes : str
            components {1,2,3,4,5,6}
        aecomp_name : str
            name of the AECOMP/AECOMPL entry
        xyz : list[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, Coord; default=0
           int : coordinate system
        cd : int; default=None -> cp
            the coordinate system for load outputs
        comment : str; default=''
            a comment for the card

        Notes
        -----
        CD - MSC specific field

        """
        monitor_point = self.monpnt1.add(
            name, label, axes, aecomp_name, xyz, cp=cp, cd=cd,
            comment=comment)
        return self.monpnt1

    def add_monpnt2(self, name: str, label: str, table: int, element_type: str,
                    nddl_item: int, eid: int, comment: str='') -> int:
        """Creates a MONPNT2 card"""
        monitor_point = self.monpnt2.add(name, label, table, element_type, nddl_item, eid,
                                         comment=comment)
        return monitor_point

    def add_monpnt3(self, name: str, label: str, axes: str,
                    grid_set: int, xyz: list[float],
                    elem_set: int=0,
                    cp: int=0, cd: Optional[int]=None,
                    xflag=None, comment: str='') -> int:
        """Creates a MONPNT3 card"""
        monitor_point = self.monpnt3.add(name, label, axes, grid_set, xyz,
                         elem_set=elem_set,
                         cp=cp, cd=cd, xflag=xflag, comment=comment)
        return monitor_point

    def add_tic(self, sid: int, nodes: list[int], components: list[int],
                u0: float=0., v0: float=0., comment: str='') -> int:
        """
        Creates a TIC card

        Parameters
        ----------
        sid : int
            Case Control IC id
        nodes : int / list[int]
            the nodes to which apply the initial conditions
        components : int / list[int]
            the DOFs to which apply the initial conditions
        u0 : float / list[float]
            Initial displacement.
        v0 : float / list[float]
            Initial velocity.
        comment : str; default=''
            a comment for the card

        """
        tic = self.tic.add(sid, nodes, components, u0=u0, v0=v0, comment=comment)
        return tic

    def add_tstep1(self, sid, tend, ninc, nout, comment: str='') -> TSTEP1:
        """Creates a TSTEP1 card"""
        tstep1 = TSTEP1(sid, tend, ninc, nout, comment=comment)
        self._add_methods.add_tstep_object(tstep1)
        return tstep1

    def add_tstep(self, sid, N, DT, NO, comment: str='') -> TSTEP:
        """Creates a TSTEP card"""
        tstep = TSTEP(sid, N, DT, NO, comment=comment)
        self._add_methods.add_tstep_object(tstep)
        return tstep

    def add_tstepnl(self, sid, ndt, dt, no, method='ADAPT', kstep=None,
                    max_iter=10, conv='PW', eps_u=1.e-2, eps_p=1.e-3,
                    eps_w=1.e-6, max_div=2, max_qn=10, max_ls=2,
                    fstress=0.2, max_bisect=5, adjust=5, mstep=None,
                    rb=0.6, max_r=32., utol=0.1, rtol_b=20.,
                    min_iter=None, comment: str='') -> TSTEPNL:
        """Creates a TSTEPNL card"""
        tstepnl = TSTEPNL(sid, ndt, dt, no, method=method, kstep=kstep,
                          max_iter=max_iter, conv=conv, eps_u=eps_u, eps_p=eps_p,
                          eps_w=eps_w, max_div=max_div, max_qn=max_qn, max_ls=max_ls,
                          fstress=fstress, max_bisect=max_bisect, adjust=adjust,
                          mstep=mstep, rb=rb, max_r=max_r, utol=utol, rtol_b=rtol_b,
                          min_iter=min_iter, comment=comment)
        self._add_methods.add_tstepnl_object(tstepnl)
        return tstepnl

    def add_nlparm(self, nlparm_id, ninc=None, dt=0.0, kmethod='AUTO', kstep=5,
                   max_iter=25, conv='PW', int_out='NO', eps_u=0.01,
                   eps_p=0.01, eps_w=0.01, max_div=3, max_qn=None, max_ls=4,
                   fstress=0.2, ls_tol=0.5, max_bisect=5, max_r=20., rtol_b=20.,
                   comment: str='') -> NLPARM:
        """Creates an NLPARM card"""
        nlparm = NLPARM(nlparm_id, ninc=ninc, dt=dt, kmethod=kmethod, kstep=kstep,
                        max_iter=max_iter, conv=conv, int_out=int_out, eps_u=eps_u,
                        eps_p=eps_p, eps_w=eps_w, max_div=max_div, max_qn=max_qn,
                        max_ls=max_ls, fstress=fstress, ls_tol=ls_tol,
                        max_bisect=max_bisect, max_r=max_r, rtol_b=rtol_b, comment=comment)
        self._add_methods.add_nlparm_object(nlparm)
        return nlparm

    def add_nlpci(self, nlpci_id, Type='CRIS', minalr=0.25, maxalr=4.,
                  scale=0., desiter=12, mxinc=20, comment: str='') -> NLPCI:
        """Creates an NLPCI card"""
        nlpci = NLPCI(nlpci_id, Type=Type, minalr=minalr, maxalr=maxalr,
                      scale=scale, desiter=desiter, mxinc=mxinc, comment=comment)
        self._add_methods.add_nlpci_object(nlpci)
        return nlpci

    def add_delay(self, delay_id: int, node: int, component: int, delay: float,
                  comment: str='') -> int:
        """
        Creates a DELAY card

        Parameters
        ----------
        delay_id : int
            DELAY id that is referenced by a TLOADx, RLOADx or ACSRCE card
        nodes : int
            node that see the delay
        component : int
            the component corresponding to the node that see the delay
        delays : list[float]
            Time delay (tau) for designated point Pi and component Ci
        comment : str; default=''
            a comment for the card

        """
        delay = self.delay.add(delay_id, node, component, delay, comment=comment)
        return delay

    def add_dphase(self, dphase_id: int, node: int, component: int, phase_lead: float,
                   comment: str='') -> int:
        """
        Creates a DPHASE card

        Parameters
        ----------
        dphase_id : int
            DPHASE id that is referenced by a RLOADx or ACSRCE card
        node : int
            node that see the delay
        component : int
            the component corresponding to the node that see the delay
        phase_lead : float
            Phase lead θ in degrees.
        comment : str; default=''
            a comment for the card

        """
        dphase = self.dphase.add(dphase_id, node, component, phase_lead, comment=comment)
        return dphase

    def add_rotorg(self, rotor_id: int, nids: list[int], comment: str='') -> int:
        """Creates a ROTORG card"""
        rotor = self.rotorg.add(rotor_id, nids, comment=comment)
        return rotor

    def add_rotord(self, sid, rstart, rstep, numstep,
                   rids, rsets, rspeeds, rcords, w3s, w4s, rforces, brgsets,
                   refsys='ROT', cmout=0.0, runit='RPM', funit='RPM',
                   zstein='NO', orbeps=1.e-6, roprt=0, sync=1, etype=1,
                   eorder=1.0, threshold=0.02, maxiter=10, comment: str='') -> ROTORD:
        """Creates a ROTORD card"""
        rotor = ROTORD(sid, rstart, rstep, numstep, rids, rsets, rspeeds,
                       rcords, w3s, w4s, rforces,
                       brgsets, refsys=refsys, cmout=cmout,
                       runit=runit, funit=funit, zstein=zstein, orbeps=orbeps,
                       roprt=roprt, sync=sync, etype=etype, eorder=eorder,
                       threshold=threshold, maxiter=maxiter, comment=comment)
        self._add_methods.add_rotor_object(rotor)
        return rotor

    def add_bulk_lines(self, lines: list[str]) -> None:
        if isinstance(lines, str):
            lines = [lines]
        self.reject_lines.append(lines)
        #self.reject_card_lines('dummy', lines, show_log=True)

    def add_include_file(self, include_filename: str,
                         is_windows: Optional[bool]=None) -> None:
        lines = write_include(include_filename, is_windows=is_windows).rstrip().split('\n')
        self.reject_lines.append(lines)
        #self.reject_card_lines('INCLUDE', lines, show_log=True)

    def add_panel(self, names: list[str], set_ids: list[int]) -> None:
        fields = ['PANEL']
        for name, set_id in zip(names, set_ids):
            fields.extend([name, set_id])
        self.reject_card_lines('PANEL', print_card_8(fields).split('\n'), show_log=False)

    def add_cmfree(self, eid, s, s2, y, n) -> None:
        fields = ['CMFREE', eid, s, s2, y, n]
        self.reject_card_lines('CMFREE', print_card_8(fields).split('\n'), show_log=False)

    def add_cfluid2(self, eid, ringfls, rho, b, harmonic) -> None:
        fields = ['CFLUID2', eid] + ringfls + [rho, b, harmonic]
        self.reject_card_lines('CFLUID2', print_card_8(fields).split('\n'), show_log=False)

    def add_cfluid3(self, eid, ringfls, rho, b, harmonic) -> None:
        fields = ['CFLUID3', eid] + ringfls + [rho, b, harmonic]
        self.reject_card_lines('CFLUID3', print_card_8(fields).split('\n'), show_log=False)

    def add_cfluid4(self, eid, ringfls, rho, b, harmonic) -> None:
        fields = ['CFLUID4', eid] + ringfls + [rho, b, harmonic]
        self.reject_card_lines('CFLUID4', print_card_8(fields).split('\n'), show_log=False)

    def add_rgyro(self, sid, asynci, refrot, unit, speed_low, speed_high, speed,
                  comment: str='') -> None:
        """Creates an RGYRO card"""
        fields = ['RGYRO', sid, asynci, refrot, unit, speed_low, speed_high, speed]
        self.reject_card_lines('RGYRO', print_card_8(fields).split('\n'), show_log=False)

    def add_rspint(self, rid, grida, gridb, gr, unit, table_id, comment: str='') -> None:
        """Creates an RSPINT card"""
        fields = ['RSPINT', rid, grida, gridb, gr, unit, table_id]
        self.reject_card_lines('RSPINT', print_card_8(fields).split('\n'), show_log=False)

    def add_dti(self, name, fields, comment: str='') -> DTI | DTI_UNITS:
        """Creates a DTI card"""
        if name == 'UNITS':
            dti = DTI_UNITS(name, fields, comment=comment)
        else:
            dti = DTI(name, fields, comment=comment)
        self._add_methods.add_dti_object(dti)
        return dti

    def add_dmig_uaccel(self, tin, ncol, load_sequences, comment: str='') -> DMIG_UACCEL:
        """Creates a DMIG,UACCEL card"""
        dmig = DMIG_UACCEL(tin, ncol, load_sequences, comment=comment)
        self._add_methods.add_dmig_object(dmig)
        return dmig

    def add_dmig(self, name, ifo, tin, tout, polar, ncols, GCj, GCi,
                 Real, Complex=None, comment: str='') -> DMIG:
        """
        Creates a DMIG card

        Parameters
        ----------
        name : str
            the name of the matrix
        ifo : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : list[(node, dof)]
            the [jnode, jDOFs]
        GCi  : list[(node, dof)]
            the inode, iDOFs
        Real : list[float]
            The real values
        Complex : list[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        dmig = DMIG(name, ifo, tin, tout, polar, ncols, GCj, GCi,
                    Real, Complex, comment=comment)
        self._add_methods.add_dmig_object(dmig)
        return dmig

    def add_dmi(self, name, form, tin, tout, nrows, ncols, GCj, GCi,
                Real, Complex=None, comment: str='') -> DMI:
        """Creates a DMI card"""
        dmi = DMI(name, form, tin, tout, nrows, ncols, GCj, GCi, Real,
                  Complex, comment=comment)
        self._add_methods.add_dmi_object(dmi)
        return dmi

    def add_dmiax(self, name, matrix_form, tin, tout, ncols,
                  GCNj, GCNi, Real, Complex=None, comment: str='') -> DMIAX:
        """Creates a DMIAX card"""
        dmiax = DMIAX(name, matrix_form, tin, tout, ncols,
                      GCNj, GCNi, Real, Complex=Complex, comment=comment)
        self._add_methods.add_dmiax_object(dmiax)
        return dmiax

    def add_dmij(self, name, form, tin, tout, nrows, ncols, GCj, GCi,
                 Real, Complex=None, comment: str='') -> DMIJ:
        """Creates a DMIJ card"""
        dmij = DMIJ(name, form, tin, tout, nrows, ncols, GCj, GCi,
                    Real, Complex, comment=comment)
        self._add_methods.add_dmij_object(dmij)
        return dmij

    def add_dmiji(self, name, ifo, tin, tout, nrows, ncols, GCj, GCi,
                  Real, Complex=None, comment: str='') -> DMIJI:
        """
        | DMIJI | NAME | 0 | IFO | TIN | TOUT POLAR | | NCOL |
        """
        dmiji = DMIJI(name, ifo, tin, tout, nrows, ncols, GCj, GCi,
                      Real, Complex, comment=comment)
        self._add_methods.add_dmiji_object(dmiji)
        return dmiji

    def add_dmik(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment: str='') -> DMIK:
        """
        Creates a DMIK card

        Parameters
        ----------
        name : str
            the name of the matrix
        ifo : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : list[(node, dof)]
            the jnode, jDOFs
        GCi  : list[(node, dof)]
            the inode, iDOFs
        Real : list[float]
            The real values
        Complex : list[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        dmik = DMIK(name, ifo, tin, tout, polar, ncols,
                    GCj, GCi, Real, Complex, comment=comment)
        self._add_methods.add_dmik_object(dmik)
        return dmik

    #---------------------------------------------------------------------

    def add_nlrsfd(self, sid, ga, gb, plane, bdia, blen, bclr, soln,
                   visco, pvapco, nport,
                   pres1, theta1, pres2, theat2, npnt,
                   offset1, offset2) -> None:
        fields = [
            'NLRSFD', sid, ga, gb, plane, bdia, blen, bclr, soln,
            visco, pvapco, nport,
            pres1, theta1, pres2, theat2, npnt,
            offset1, offset2]
        self.reject_card_lines('NLRSFD', print_card_8(fields).split('\n'), show_log=False)

    def add_nolin1(self, sid, gi, ci, s, gj, cj, t) -> None:
        fields = ['NOLIN1', sid, gi, ci, s, gj, cj, t]
        self.reject_card_lines('NOLIN1', print_card_8(fields).split('\n'), show_log=False)

    def add_nolin2(self, sid, g, ci, s, gj, cj, gk, ck) -> None:
        fields = ['NOLIN2', sid, g, ci, s, gj, cj, gk, ck]
        self.reject_card_lines('NOLIN2', print_card_8(fields).split('\n'), show_log=False)

    def add_nolin3(self, sid, gi, ci, s, gj, cj, a) -> None:
        fields = ['NOLIN3', sid, gi, ci, s, gj, cj, a]
        self.reject_card_lines('NOLIN3', print_card_8(fields).split('\n'), show_log=False)

    def add_nolin4(self, sid, gi, ci, s, gj, cj, a) -> None:
        fields = ['NOLIN4', sid, gi, ci, s, gj, cj, a]
        self.reject_card_lines('NOLIN4', print_card_8(fields).split('\n'), show_log=False)

    def add_dynred(self, sid, fmax, nirv, nit, idir, nqdes) -> None:
        fields = ['DYNRED', sid, fmax, nirv, nit, idir, nqdes]
        self.reject_card_lines('DYNRED', print_card_8(fields).split('\n'), show_log=False)

    def add_rcross(self, sid, rtype1, id1, comp1, rtype2, id2, comp2, curid) -> None:
        fields = ['RCROSS', sid, rtype1, id1, comp1, rtype2, id2, comp2, curid]
        self.reject_card_lines('RCROSS', print_card_8(fields).split('\n'), show_log=False)

    def add_uxvec(self, idi: int, labels: list[str], uxs: list[float], comment: str='') -> UXVEC:
        uxvec = UXVEC(idi, labels, uxs, comment=comment)
        self._add_methods.add_uxvec_object(uxvec)
        return uxvec

    #----------------------------------------------------------------------------------
    # parametric
    #def add_pset(self, idi, poly1, poly2, poly3, cid, typei, typeids, comment: str='') -> PSET:
        #"""PSET ID POLY1 POLY2 POLY3 CID SETTYP ID"""
        #pset = PSET(idi, poly1, poly2, poly3, cid, typei, typeids, comment=comment)
        #self.pset[idi] = pset
        #return pset

    #def add_pval(self, idi, poly1, poly2, poly3, cid, typei, typeids, comment: str='') -> PVAL:
        #"""PVAL ID POLY1 POLY2 POLY3 CID SETTYP ID"""
        #pval = PVAL(idi, poly1, poly2, poly3, cid, typei, typeids, comment=comment)
        #self._add_methods.add_pval(pval, allow_overwrites=False)
        #return pval

    #def add_gmcurv(self, curve_id, group, data, cid_in=0, cid_bc=0, comment: str='') -> GMCURV:
        #curve = GMCURV(curve_id, group, data, cid_in=cid_in, cid_bc=cid_bc,
                       #comment=comment)
        #self._add_methods.add_gmcurv(curve, allow_overwrites=False)
        #return curve

    #def add_gmsurf(self, curve_id, group, data, cid_in=0, cid_bc=0, comment: str='') -> GMSURF:
        #surf = GMSURF(curve_id, group, data, cid_in=cid_in, cid_bc=cid_bc, comment=comment)
        #self._add_methods.add_gmsurf(surf, allow_overwrites=False)
        #return surf

    #def add_feedge(self, edge_id, nids, cid, geom_ids, geomin='POINT', comment: str='') -> FEEDGE:
        #edge = FEEDGE(edge_id, nids, cid, geom_ids, geomin=geomin, comment=comment)
        #self._add_methods.add_feedge(edge, allow_overwrites=False)
        #return edge

    #def add_feface(self, face_id, nids, cid, surf_ids, comment: str='') -> FEFACE:
        #face = FEFACE(face_id, nids, cid, surf_ids, comment=comment)
        #self._add_methods.add_feface(face, allow_overwrites=False)
        #return face
    #----------------------------------------------------------------------------------
    # cyclic
    def add_cyax(self, nids, comment: str='') -> CYAX:
        cyax = CYAX(nids, comment=comment)
        self._add_methods.add_cyax_object(cyax)
        return cyax

    def add_cyjoin(self, side, coord, nids, comment: str='') -> CYJOIN:
        cyjoin = CYJOIN(side, coord, nids, comment=comment)
        self._add_methods.add_cyjoin_object(cyjoin)
        return cyjoin

    #def add_cysym(self, side, coord, nids, comment: str='') -> CYSYM:
        #cysym = CYSYM(side, coord, nids, comment=comment)
        #self._add_cysym_object(cysym)
        #return cysym
