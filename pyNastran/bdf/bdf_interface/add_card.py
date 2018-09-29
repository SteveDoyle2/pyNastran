# coding: utf-8
# pylint: disable=R0913, R0914, C0103
"""
Defines a method to add a card that is faster than add_card
and far less error prone for a user

That said, there are still a few bugs.
"""
from __future__ import print_function

from typing import Any, Optional, List, Union
from six import string_types
import numpy as np

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.add_methods import AddMethods

from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP, PRAC2D, PRAC3D
from pyNastran.bdf.cards.properties.solid import PLSOLID, PSOLID, PIHEX, PCOMPS
from pyNastran.bdf.cards.msgmesh import CGEN

from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.elements.solid import (
    #CTETRA, CPYRAM, CPENTA, CHEXA,
    CIHEX1, CIHEX2,
    CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
    CTETRA10, CPYRAM13, CPENTA15, CHEXA20,
)
from pyNastran.bdf.cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE, RSSCON

from pyNastran.bdf.cards.axisymmetric.axisymmetric import (
    AXIC, RINGAX, POINTAX, CCONEAX, PCONEAX, PRESAX, TEMPAX,)
from pyNastran.bdf.cards.elements.axisymmetric_shells import (
    CTRAX3, CTRAX6, CTRIAX, CTRIAX6, CQUADX, CQUADX4, CQUADX8)
from pyNastran.bdf.cards.elements.shell import (
    CQUAD, CQUAD4, CQUAD8, CQUADR, CSHEAR,
    CTRIA3, CTRIA6, CTRIAR,
    CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
    CPLSTS3, #CPLSTS4,
    SNORM,
)
from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE, PPLANE
from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D, PBUSHT #PBUSH2D
from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                                 CDAMP5)
from pyNastran.bdf.cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE
from pyNastran.bdf.cards.elements.bars import CBAR, CBARAO, CBEAM3, CBEND#, BAROR
from pyNastran.bdf.cards.elements.beam import CBEAM#, BEAMOR
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBRSECT, PBEND, PBEAM3
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP, PBMSECT
# CMASS5
from pyNastran.bdf.cards.elements.mass import CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4
from pyNastran.bdf.cards.properties.mass import PMASS, NSM, NSM1, NSML, NSML1, NSMADD
from pyNastran.bdf.cards.constraints import (SPC, SPCADD, SPCAX, SPC1, SPCOFF, SPCOFF1,
                                             MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                             GMSPC)
from pyNastran.bdf.cards.coordinate_systems import (CORD1R, CORD1C, CORD1S,
                                                    CORD2R, CORD2C, CORD2S, #CORD3G,
                                                    GMCORD)
from pyNastran.bdf.cards.deqatn import DEQATN
from pyNastran.bdf.cards.dynamic import (
    DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
    TSTEP, TSTEP1, TSTEPNL, NLPARM, NLPCI, TF, ROTORG, ROTORD, TIC)
from pyNastran.bdf.cards.loads.loads import (
    LSEQ, SLOAD, DAREA, RFORCE, RFORCE1, SPCD, LOADCYN, DEFORM)
from pyNastran.bdf.cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
from pyNastran.bdf.cards.loads.static_loads import (LOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                    FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                    PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                                    GMLOAD)
from pyNastran.bdf.cards.loads.random_loads import (
    RANDPS, RANDT1)

from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11, MAT3D,
                                           MATG, MATHE, MATHP, CREEP, NXSTRAT, EQUIV)
from pyNastran.bdf.cards.material_deps import (
    MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, MATS1)

from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL
from pyNastran.bdf.cards.nodes import GRID, GRDSET, SPOINTs, EPOINTs, POINT, SEQGP
from pyNastran.bdf.cards.aero.aero import (
    AECOMP, AEFACT, AELINK, AELIST, AEPARM, AESURF, AESURFS,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    MONPNT1, MONPNT2, MONPNT3,
    SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5)
from pyNastran.bdf.cards.aero.static_loads import AESTAT, AEROS, CSSCHD, TRIM, TRIM2, DIVERG
from pyNastran.bdf.cards.aero.dynamic_loads import AERO, FLFACT, FLUTTER, GUST, MKAERO1, MKAERO2
from pyNastran.bdf.cards.optimization import (
    DCONADD, DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
    DRESP1, DRESP2, DRESP3,
    DVCREL1, DVCREL2,
    DVMREL1, DVMREL2,
    DVPREL1, DVPREL2,
    DVGRID, DSCREEN)
from pyNastran.bdf.cards.bdf_sets import (
    ASET, BSET, CSET, QSET, USET,
    ASET1, BSET1, CSET1, QSET1, USET1,
    OMIT1,
    SET1, SET3, #RADSET,
    SEBSET, SECSET, SEQSET, # SEUSET
    SEBSET1, SECSET1, SEQSET1, # SEUSET1
    SESET, #SEQSEP,
    RADSET,
)
from pyNastran.bdf.cards.params import PARAM
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL, DTI
from pyNastran.bdf.cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD,
                                               QVOL, QVECT)
from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                 PHBDY, CONV, CONVM, TEMPBC)
from pyNastran.bdf.cards.thermal.radiation import RADM, RADBC, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                            TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                            TABLES1, TABDMP1, TABLEST, TABRND1, TABRNDG,
                                            DTABLE)
from pyNastran.bdf.cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA, BCONP, BLSEG
from pyNastran.utils.numpy_utils import integer_string_types
from pyNastran.bdf.cards.superelements import (
    CSUPER, CSUPEXT, SEBNDRY, SECONCT, SEELT, SEEXCLD, SELABEL, SELOAD,
    SELOC, SEMPLN, SENQSET, SETREE)

class AddCards(AddMethods):
    """defines the add_cardname functions that use the object inits"""
    def __init__(self):
        AddMethods.__init__(self)

    def add_grid(self, nid, xyz, cp=0, cd=0, ps='', seid=0, comment=''):
        # type: (int, Union[None, List[float], np.ndarray], int, int, str, int, str) -> None
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
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card

        """
        grid = GRID(nid, xyz, cp=cp, cd=cd, ps=ps, seid=seid, comment=comment)
        self._add_node_object(grid)
        return grid

    def add_grdset(self, cp, cd, ps, seid, comment=''):
        # type: (int, int, str, int, str) -> GRDSET
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

    def add_seqgp(self, nids, seqids, comment=''):
        # type: (List[int], List[Union[int, float]], str) -> SEQGP
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
        self._add_seqgp_object(seqgp)
        return seqgp

    def add_spoint(self, ids, comment=''):
        # type: (Union[int, List[int]], str) -> SPOINTs
        """
        Creates the SPOINTs card that contains many SPOINTs

        Parameters
        ----------
        ids : List[int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        spoint = SPOINTs(ids, comment=comment)
        self._add_spoint_object(spoint)
        return spoint

    def add_epoint(self, ids, comment=''):
        # type: (Union[int, List[int]], str) -> EPOINTs
        """
        Creates the EPOINTs card that contains many EPOINTs

        Parameters
        ----------
        ids : List[int]
            EPOINT ids
        comment : str; default=''
            a comment for the card

        """
        epoint = EPOINTs(ids, comment=comment)
        self._add_epoint_object(epoint)
        return epoint

    def add_point(self, nid, xyz, cp=0, comment=''):
        # type: (int, int, Any, str) -> POINT
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
        point = POINT(nid, xyz, cp=cp, comment=comment)
        self._add_point_object(point)
        return point

    def add_cord2r(self, cid, origin, zaxis, xzplane, rid=0, comment=''):
        # type: (int, Optional[Union[List[float], np.ndarray]], Any, Any, int, str) -> CORD2R
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
        origin : List[float, float, float]
            the origin of the coordinate system
        zaxis : List[float, float, float]
            the z-axis of the coordinate system
        xzplane : List[float, float, float]
            a point on the xz plane
        comment : str; default=''
            a comment for the card

        """
        coord = CORD2R(cid, origin, zaxis, xzplane, rid=rid, comment=comment)
        self._add_coord_object(coord)
        return coord

    def add_cord2c(self, cid, origin, zaxis, xzplane, rid=0, comment=''):
        # type: (int, Optional[Union[List[float], np.ndarray]], Any, Any, int, str) -> CORD2C
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
        origin : List[float, float, float]
            the origin of the coordinate system
        zaxis : List[float, float, float]
            the z-axis of the coordinate system
        xzplane : List[float, float, float]
            a point on the xz plane
        comment : str; default=''
            a comment for the card

        """
        coord = CORD2C(cid, origin, zaxis, xzplane, rid=rid, comment=comment)
        self._add_coord_object(coord)
        return coord

    def add_cord2s(self, cid, origin, zaxis, xzplane, rid=0, comment=''):
        # type: (int, Optional[Union[List[float], np.ndarray]], Any, Any, int, str) -> CORD2S
        """
        Creates the CORD2C card, which defines a spherical coordinate
        system using 3 vectors.

        Parameters
        ----------
        cid : int
            coordinate system id
        origin : List[float, float, float]
            the origin of the coordinate system
        zaxis : List[float, float, float]
            the z-axis of the coordinate system
        xzplane : List[float, float, float]
            a point on the xz plane
        rid : int; default=0
            the referenced coordinate system that defines the system the
            vectors
        comment : str; default=''
            a comment for the card

        """
        coord = CORD2S(cid, rid=rid, origin=origin, zaxis=zaxis, xzplane=xzplane,
                       comment=comment)
        self._add_coord_object(coord)
        return coord

    def add_cord1r(self, cid, g1, g2, g3, comment=''):
        # type: (int, int, int, int, str) -> CORD1R
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
        coord = CORD1R(cid, g1, g2, g3, comment=comment)
        self._add_coord_object(coord)
        return coord

    def add_cord1c(self, cid, g1, g2, g3, comment=''):
        # type: (int, int, int, int, str) -> CORD1C
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
        coord = CORD1C(cid, g1, g2, g3, comment=comment)
        self._add_coord_object(coord)
        return coord

    def add_cord1s(self, cid, g1, g2, g3, comment=''):
        # type: (int, int, int, int, str) -> CORD1S
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
        coord = CORD1S(cid, g1, g2, g3, comment=comment)
        self._add_coord_object(coord)
        return coord

    def add_param(self, key, values, comment=''):
        # type: (str, List[Union[int, float, str]], str) -> PARAM
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/List
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        param = PARAM(key, values, comment=comment)
        self._add_param_object(param)
        return param

    def add_plotel(self, eid, nodes, comment=''):
        # type: (int, List[int], str) -> PLOTEL
        """
        Adds a PLOTEL card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : List[int, int]
            Unique GRID point IDs
        comment : str; default=''
            a comment for the card

        """
        elem = PLOTEL(eid, nodes, comment=comment)
        self._add_plotel_object(elem)
        return elem

    def add_conm1(self, eid, nid, mass_matrix, cid=0, comment=''):
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
        mass = CONM1(eid, nid, mass_matrix, cid=cid, comment=comment)
        self._add_mass_object(mass)
        return mass

    def add_conm2(self, eid, nid, mass, cid=0, X=None, I=None, comment=''):
        # type: (int, int, float, int, Any, Any, str) -> CONM2
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
        X : (3, ) List[float]; default=None -> [0., 0., 0.]
            xyz offset vector relative to nid
        I : (6, ) List[float]; default=None -> [0., 0., 0., 0., 0., 0.]
            mass moment of inertia matrix about the CG
            I11, I21, I22, I31, I32, I33 = I
        comment : str; default=''
            a comment for the card

        """
        mass_obj = CONM2(eid, nid, mass, cid=cid, X=X, I=I, comment=comment)
        self._add_mass_object(mass_obj)
        return mass_obj

    def add_nsm(self, sid, nsm_type, pid_eid, value, comment=''):
        # type: (int, str, int, float, str) -> NSM
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
        pid_eid : List[int]; int
            property id or element id depending on nsm_type
        value : List[float]; float
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

        nsms = []
        for pid_eidi, valuei in zip(pid_eid, value):
            nsm = NSM(sid, nsm_type, pid_eidi, valuei, comment=comment)
            self._add_nsm_object(nsm)
            nsms.append(nsm)
        return nsms

    def add_nsm1(self, sid, nsm_type, value, ids, comment=''):
        # type: (int, str, float, List[int], str) -> NSM1
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
        ids : List[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card

        """
        assert isinstance(value, float)
        nsm = NSM1(sid, nsm_type, value, ids, comment=comment)
        self._add_nsm_object(nsm)
        return nsm

    def add_nsml(self, sid, nsm_type, pid_eid, value, comment=''):
        # type: (int, str, int, float, str) -> NSML
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
        pid_eid : List[int]; int
            property id or element id depending on nsm_type
        value : List[float]; float
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

        nsms = []
        for pid_eidi, valuei in zip(pid_eid, value):
            nsm = NSML(sid, nsm_type, pid_eidi, valuei, comment=comment)
            self._add_nsm_object(nsm)
            nsms.append(nsm)
        return nsms

    def add_nsml1(self, sid, nsm_type, value, ids, comment=''):
        # type: (int, str, float, List[int], str) -> NSML1
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
        ids : List[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card

        """
        nsm = NSML1(sid, nsm_type, value, ids, comment=comment)
        self._add_nsm_object(nsm)
        return nsm

    def add_nsmadd(self, sid, sets, comment=''):
        # type: (int, List[int], str) -> NSMADD
        """
        Creates an NSMADD card, which sum NSM sets

        Parameters
        ----------
        sid : int
            the NSM Case Control value
        sets : List[int]
            the NSM, NSM1, NSML, NSML1 values
        comment : str; default=''
            a comment for the card

        """
        nsmadd = NSMADD(sid, sets, comment=comment)
        self._add_nsmadd_object(nsmadd)
        return nsmadd

    def add_pmass(self, pid, mass, comment=''):
        # type: (int, float, str) -> PMASS
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
        prop = PMASS(pid, mass, comment=comment)
        self._add_property_mass_object(prop)
        return prop

    def add_cmass1(self, eid, pid, nids, c1=0, c2=0, comment=''):
        # type: (int, int, [int, int], int, int, str) -> CMASS1
        """
        Creates a CMASS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        mass_obj = CMASS1(eid, pid, nids, c1, c2, comment=comment)
        self._add_mass_object(mass_obj)
        return mass_obj

    def add_cmass2(self, eid, mass, nids, c1, c2, comment=''):
        # type: (int, float, [int, int], int, int, str) -> CMASS2
        """
        Creates a CMASS2 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            mass
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        mass_obj = CMASS2(eid, mass, nids, c1, c2, comment=comment)
        self._add_mass_object(mass_obj)
        return mass_obj

    def add_cmass3(self, eid, pid, nids, comment=''):
        # type: (int, int, [int, int], str) -> CMASS3
        """
        Creates a CMASS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        mass = CMASS3(eid, pid, nids, comment=comment)
        self._add_mass_object(mass)
        return mass

    def add_cmass4(self, eid, mass, nids, comment=''):
        # type: (int, float, [int, int], str) -> CMASS4
        """
        Creates a CMASS4 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            SPOINT mass
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        mass_obj = CMASS4(eid, mass, nids, comment=comment)
        self._add_mass_object(mass_obj)
        return mass_obj

    def add_pelas(self, pid, k, ge=0., s=0., comment=''):
        # type: (int, float, float, float, str) -> PELAS
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
        prop = PELAS(pid, k, ge, s, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_celas1(self, eid, pid, nids, c1=0, c2=0, comment=''):
        # type: (int, int, List[int], int, int, str) -> CELAS1
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        elem = CELAS1(eid, pid, nids, c1, c2, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_celas2(self, eid, k, nids, c1=0, c2=0, ge=0., s=0., comment=''):
        # type: (int, float, List[int], int, int, float, float, str) -> CELAS2
        """
        Creates a CELAS2 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : List[int, int]
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
        elem = CELAS2(eid, k, nids, c1=c1, c2=c2, ge=ge, s=s, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_celas3(self, eid, pid, nids, comment=''):
        # type: (int, int, List[int], str) -> CELAS3
        """
        Creates a CELAS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = CELAS3(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_celas4(self, eid, k, nids, comment=''):
        # type: (int, float, List[int], str) -> CELAS4
        """
        Creates a CELAS4 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = CELAS4(eid, k, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cdamp1(self, eid, pid, nids, c1=0, c2=0, comment=''):
        # type: (int, int, List[int], int, int, str) -> CDAMP1
        """
        Creates a CDAMP1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        elem = CDAMP1(eid, pid, nids, c1=c1, c2=c2, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cdamp2(self, eid, b, nids, c1=0, c2=0, comment=''):
        # type: (int, float, List[int], int, int, str) -> CDAMP2
        """
        Creates a CDAMP2 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : List[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        elem = CDAMP2(eid, b, nids, c1=c1, c2=c2, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cdamp3(self, eid, pid, nids, comment=''):
        # type: (int, int, List[int], str) -> CDAMP3
        """
        Creates a CDAMP3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = CDAMP3(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cdamp4(self, eid, b, nids, comment=''):
        # type: (int, float, List[int], str) -> CDAMP4
        """
        Creates a CDAMP4 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = CDAMP4(eid, b, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cdamp5(self, eid, pid, nids, comment=''):
        # type: (int, float, List[int], str) -> CDAMP5
        """
        Creates a CDAMP5 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP5)
        nids : List[int, int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        elem = CDAMP5(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pdamp(self, pid, b, comment=''):
        # type: (int, float, str) -> PDAMP
        """Creates a PDAMP card"""
        prop = PDAMP(pid, b, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pdampt(self, pid, tbid, comment=''):
        # type: (int, int, str) -> PDAMPT
        """Creates a PDAMPT card"""
        prop = PDAMPT(pid, tbid, comment=comment)
        self._add_pdampt_object(prop)
        return prop

    def add_pdamp5(self, pid, mid, b, comment=''):
        # type: (int, int, float, str) -> PDAMP5
        """Creates a PDAMP5 card"""
        prop = PDAMP5(pid, mid, b, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cvisc(self, eid, pid, nids, comment=''):
        # type: (int, int, List[int], str) -> CVISC
        """
        Creates a CVISC card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PVISC)
        nids : List[int, int]
            GRID ids
        comment : str; default=''
            a comment for the card

        """
        elem = CVISC(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pvisc(self, pid, ce, cr, comment=''):
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
        # type: (int, float, float, str) -> PVISC
        prop = PVISC(pid, ce, cr, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cgap(self, eid, pid, nids, x, g0, cid=None, comment=''):
        # type: (int, int, int, List[int], Any, Optional[int], str) -> CGAP
        """
        Creates a CGAP card

        Parameters
        ----------
        eid : int
            Element ID
        pid : int
            Property ID (PGAP)
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        x : List[float, float, float]
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
        elem = CGAP(eid, pid, nids, x, g0, cid=cid, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pgap(self, pid, u0=0., f0=0., ka=1.e8, kb=None, mu1=0.,
                 kt=None, mu2=None, tmax=0., mar=100., trmin=0.001,
                 comment=''):
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
        # type: (int, float, float, float, Any, float, Any, Any, float, float, float, str) -> PGAP
        prop = PGAP(pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin,
                    comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cfast(self, eid, pid, Type, ida, idb, gs=None, ga=None, gb=None,
                  xs=None, ys=None, zs=None, comment=''):
        # type: (int, int, str, int, int, Any, Any, Any, Any, Any, Any, str) -> CFAST
        """Creates a CFAST card"""
        elem = CFAST(eid, pid, Type, ida, idb, gs=gs, ga=ga, gb=gb,
                     xs=xs, ys=ys, zs=zs, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pfast(self, pid, d, kt1, kt2, kt3, mcid=-1, mflag=0,
                  kr1=0., kr2=0., kr3=0., mass=0., ge=0., comment=''):
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
        mcid : int; default=01
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
        prop = PFAST(pid, d, kt1, kt2, kt3, mcid=mcid, mflag=mflag,
                     kr1=kr1, kr2=kr2, kr3=kr3, mass=mass, ge=ge, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cbush(self, eid, pid, nids, x, g0, cid=None,
                  s=0.5, ocid=-1, si=None, comment=''):
        # type(int, int, List[int], Any, Optional[int], Optional[int], float, int, Optional[List[float]], str) -> CBUSH
        """
        Creates a CBUSH card

        Parameters
        ----------
        eid : int
            Element id
        pid : int
            Property id (PBUSH)
        nids : List[int, int]
            node ids; connected grid points at ends A and B
            The nodes may be coincident, but then cid is required.
        x : List[float, float, float]; None
            List : the directional vector used to define the stiffnesses
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
        si : List[float, float, float]; default=None
            Components of spring-damper offset in the OCID coordinate system
            if OCID > 0.
            None : [None, None, None]
        comment : str; default=''
            a comment for the card

        """
        elem = CBUSH(eid, pid, nids, x, g0, cid=cid, s=s, ocid=ocid, si=si, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pbush(self, pid, k, b, ge, rcv=None, mass=None, comment=''):
        """
        Creates a PBUSH card, which defines a property for a CBUSH

        Parameters
        ----------
        pid : int
            property id
        k : List[float]
            Nominal stiffness values in directions 1 through 6.
            len(k) = 6
        b : List[float]
            Nominal damping coefficients in direction 1 through 6 in units of
            force per unit velocity
            len(b) = 6
        ge : List[float]
            Nominal structural damping constant in directions 1 through 6.
            len(ge) = 6
        rcv : List[float]; default=None -> (None, None, None, None)
            [sa, st, ea, et] = rcv
            length(rcv) = 4
        mass : float; default=None
            lumped mass of the CBUSH
            This is an MSC only parameter.
        comment : str; default=''
            a comment for the card

        """
        prop = PBUSH(pid, k, b, ge, rcv=rcv, mass=mass, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cbush1d(self, eid, pid, nids, cid=None, comment=''):
        # type: (int, int, List[int], int, str) -> CBUSH1D
        """Creates a CBUSH1D card"""
        elem = CBUSH1D(eid, pid, nids, cid=cid, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cbush2d(self, eid, pid, nids, cid=0, plane='XY', sptid=None, comment=''):
        # type: (int, int, List[int], int, str, int, str) -> CBUSH2D
        """Creates a CBUSH2D card"""
        elem = CBUSH2D(eid, pid, nids, cid=cid, plane=plane, sptid=sptid, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pbush1d(self, pid, k=0., c=0., m=0., sa=0., se=0.,
                    optional_vars=None, comment=''):
        """Creates a PBUSH1D card"""
        prop = PBUSH1D(pid, k=k, c=c, m=m, sa=sa, se=se,
                       optional_vars=optional_vars, comment=comment)
        self._add_property_object(prop)
        return prop

    #def add_pbush2d(self, pid, k, c, m, sa, se, optional_vars, comment=''):
        #"""
        #Creates a PBUSH2D card
        #"""
        #prop = PBUSH2D(pid. comment=comment)
        #self._add_property_object(prop)
        #return prop

    def add_pbusht(self, pid, k_tables, b_tables, ge_tables, kn_tables,
                   comment=''):
        """Creates a PBUSHT card"""
        prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pelast(self, pid, tkid=0, tgeid=0, tknid=0, comment=''):
        # type: (int, int, int, int, str) -> PELAST
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
        prop = PELAST(pid, tkid, tgeid, tknid, comment=comment)
        self._add_pelast_object(prop)
        return prop

    def add_conrod(self, eid, mid, nids, A=0.0, j=0.0, c=0.0, nsm=0.0, comment=''):
        # type: (int, int, List[int], float, float, float, float, str) -> CONROD
        """
        Creates a CONROD card

        Parameters
        ----------
        eid : int
            element id
        mid : int
            material id
        nids : List[int, int]
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
        elem = CONROD(eid, mid, nids, A=A, j=j, c=c, nsm=nsm, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_crod(self, eid, pid, nids, comment=''):
        # type: (int, int, List[int], str) -> CROD
        """
        Creates a CROD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PROD)
        nids : List[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        elem = CROD(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_prod(self, pid, mid, A, j=0., c=0., nsm=0., comment=''):
        # type: (int, int, float, float, float, float, str) -> PROD
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
        prop = PROD(pid, mid, A, j=j, c=c, nsm=nsm, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_ctube(self, eid, pid, nids, comment=''):
        # type: (int, int, List[int], str) -> CTUBE
        """
        Creates a CTUBE card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        nids : List[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        elem = CTUBE(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ptube(self, pid, mid, OD1, t=None, nsm=0., OD2=None, comment=''):
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
        prop = PTUBE(pid, mid, OD1, t=t, nsm=nsm, OD2=OD2, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cbarao(self, eid, scale, x, comment=''):
        # type: (int, str, List[float], str) -> CBARAO
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
        x : List[float]
            the additional output locations
            len(x) <= 6
        comment : str; default=''
            a comment for the card

        Notes
        -----
        MSC only

        """
        elem_flag = CBARAO(eid, scale, x, comment=comment)
        self._add_ao_object(elem_flag, allow_overwrites=False)
        return elem_flag

    def add_cbar(self, eid, pid, nids, x, g0, offt='GGG', pa=0, pb=0,
                 wa=None, wb=None, comment=''):
        # type: (int, int, List[int], Any, int, str, int, int, Any, Any, str) -> CBAR
        """
        Adds a CBAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        x : List[float, float, float]
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
        wa / wb : List[float, float, float]
            Components of offset vectors from the grid points to the end
            points of the axis of the shear center
        comment : str; default=''
            a comment for the card

        """
        elem = CBAR(eid, pid, nids, x, g0, offt=offt, pa=pa, pb=pb,
                    wa=wa, wb=wb, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pbar(self, pid, mid, A=0., i1=0., i2=0., i12=0., j=0., nsm=0.,
                 c1=0., c2=0., d1=0., d2=0., e1=0., e2=0.,
                 f1=0., f2=0., k1=1.e8, k2=1.e8, comment=''):
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
        prop = PBAR(pid, mid, A=A, i1=i1, i2=i2, i12=i12, j=j, nsm=nsm,
                    c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2,
                    f1=f1, f2=f2, k1=k1, k2=k2, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pbarl(self, pid, mid, Type, dim, group='MSCBML0', nsm=0., comment=''):
        """
        Creates a PBARL card, which defines A, I1, I2, I12, and J using
        dimensions rather than explicit values.

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        Type : str
            type of the bar
            valid_types = {
                ROD, TUBE, I, CHAN, T, BOX, BAR, CROSS, H, T1,
                I1, CHAN1, Z, CHAN2, T2, BOX1, HEXA, HAT, HAT1, DBOX
            }
        dim : List[float]
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
        prop = PBARL(pid, mid, Type, dim, group=group, nsm=nsm, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cbeam(self, eid, pid, nids, x, g0, offt='GGG', bit=None,
                  pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment=''):
        """
        Adds a CBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        x : List[float, float, float]
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
        wa / wb : List[float, float, float]
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
        elem = CBEAM(eid, pid, nids, x, g0, offt=offt, bit=bit,
                     pa=pa, pb=pb, wa=wa, wb=wb, sa=sa, sb=sb, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pbeam(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                  c1=None, c2=None, d1=None, d2=None,
                  e1=None, e2=None, f1=None, f2=None,
                  k1=1., k2=1., s1=0., s2=0.,
                  nsia=0., nsib=None, cwa=0., cwb=None,
                  m1a=0., m2a=None, m1b=0., m2b=None,
                  n1a=0., n2a=None, n1b=0., n2b=None,
                  comment=''):
        """
        .. todo:: fix 0th entry of self.so, self.xxb

        Creates a PBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        xxb : List[float]
            The percentage locations along the beam [0., ..., 1.]
        so : List[str]
            YES, YESA, NO
        area : List[float]
            area
        i1, i2, i12, j : List[float]
            moments of inertia
        nsm : List[float]; default=None -> [0.]*nxxb
            nonstructural mass per unit length
        c1/c2, d1/d2, e1/e2, f1/f2 : List[float]; default=None -> [0.]*nxxb
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
        m1a / m2a : float; default=0. / m1a
            y/z coordinate of center of gravity of
            nonstructural mass for end A.
        m1b / m2b : float; default=0. / m1b
            y/z coordinate of center of gravity of
            nonstructural mass for end B.
        n1a / n2a : float; default=0. / n1a
            y/z coordinate of neutral axis for end A.
        n1b / n2b : float; default=0. / n1b
            y/z coordinate of neutral axis for end B.
        comment : str; default=''
            a comment for the card

        """
        prop = PBEAM(pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                     c1, c2, d1, d2, e1, e2, f1, f2,
                     k1=k1, k2=k2, s1=s1, s2=s2,
                     nsia=nsia, nsib=nsib, cwa=cwa, cwb=cwb,
                     m1a=m1a, m2a=m2a, m1b=m1b,
                     m2b=m2b, n1a=n1a, n2a=n2a, n1b=n1b, n2b=n2b, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pbcomp(self, pid, mid, y, z, c, mids,
                   area=0.0, i1=0.0, i2=0.0, i12=0.0, j=0.0, nsm=0.0,
                   k1=1.0, k2=1.0, m1=0.0, m2=0.0, n1=0.0, n2=0.0,
                   symopt=0, comment=''):
        """
        Creates a PBCOMP card

        Parameters
        ---------
        pid : int
            Property ID
        mid : int
            Material ID
        mids : List[int]
            Material ID for the i-th integration point
        y / z : List[float]
            The (y,z) coordinates of the lumped areas in the element
            coordinate system
        c : List[float]; default=0.0
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
        prop = PBCOMP(pid, mid, y, z, c, mids,
                      area, i1, i2, i12, j, nsm,
                      k1, k2, m1, m2, n1, n2,
                      symopt, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pbmsect(self, pid, mid, form, options, comment=''):
        """Creates a PBMSECT card"""
        prop = PBMSECT(pid, mid, form, options, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pbrsect(self, pid, mid, form, options, comment=''):
        """Creates a PBRSECT card"""
        prop = PBRSECT(pid, mid, form, options, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pbeam3(self, pid, mid, A, iz, iy, iyz, j, nsm=0.,
                   cy=0., cz=0., dy=0., dz=0., ey=0., ez=0., fy=0., fz=0.,
                   comment=''):
        """Creates a PBEAM3 card"""
        prop = PBEAM3(pid, mid, A, iz, iy, iyz, j, nsm=nsm,
                      cy=cy, cz=cz, dy=dy, dz=dz, ey=ey, ez=ez, fy=fy, fz=fz,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pbeaml(self, pid, mid, beam_type, xxb, dims, so=None, nsm=None,
                   group='MSCBML0', comment=''):
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
        xxb : List[float]
            The percentage locations along the beam [0., ..., 1.]
        dims : List[dim]
            dim : List[float]
                The dimensions for each section
        so : List[str]; default=None
            YES, YESA, NO
            None : [0.] * len(xxb)
        nsm : List[float]; default=None
            nonstructural mass per unit length
            None : [0.] * len(xxb)
        group : str; default='MSCBML0'
            this parameter can lead to a very broken deck with a very
            bad error message; don't touch it!
        comment : str; default=''
            a comment for the card

        """
        prop = PBEAML(pid, mid, beam_type, xxb, dims,
                      group=group, so=so, nsm=nsm, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cbend(self, eid, pid, nids, g0, x, geom, comment=''):
        """Creates a CBEND card"""
        elem = CBEND(eid, pid, nids, g0, x, geom, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pbend(self, pid, mid, beam_type, A, i1, i2, j,
                  c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                  nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b,
                  comment=''):
        """Creates a PBEND card"""
        prop = PBEND(pid, mid, beam_type, A, i1, i2, j,
                     c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b,
                     comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cbeam3(self, eid, pid, nids, x, g0, wa, wb, wc, tw, s,
                   comment=''):
        """Creates a CBEAM3 card"""
        elem = CBEAM3(eid, pid, nids, x, g0, wa, wb, wc, tw, s,
                      comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cshear(self, eid, pid, nids, comment=''):
        """
        Creates a CSHEAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHEAR)
        nids : List[int, int, int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        elem = CSHEAR(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pshear(self, pid, mid, t, nsm=0., f1=0., f2=0., comment=''):
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
        prop = PSHEAR(pid, mid, t, nsm=nsm, f1=f1, f2=f2, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_ctria3(self, eid, pid, nids, zoffset=0., theta_mcid=0.0,
                   tflag=0, T1=None, T2=None, T3=None, comment=''):
        """
        Creates a CTRIA3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int]
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
        elem = CTRIA3(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquad4(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                   tflag=0, T1=None, T2=None, T3=None, T4=None,
                   comment=''):
        """
        Creates a CQUAD4 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int]
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
        elem = CQUAD4(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ctria6(self, eid, pid, nids, theta_mcid=0., zoffset=0.,
                   tflag=0, T1=None, T2=None, T3=None, comment=''):
        """
        Creates a CTRIA6 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int/None, int/None, int/None]
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
        elem = CTRIA6(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquad8(self, eid, pid, nids, theta_mcid=0., zoffset=0.,
                   tflag=0, T1=None, T2=None, T3=None, T4=None, comment=''):
        """
        Creates a CQUAD8 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int, int/None, int/None, int/None, int/None]
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
        elem = CQUAD8(eid, pid, nids,
                      theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquad(self, eid, pid, nids, theta_mcid=0., comment=''):
        """
        Creates a CQUAD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int, int/None, int/None,
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
        elem = CQUAD(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ctriar(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.0,
                   tflag=0, T1=None, T2=None, T3=None, comment=''):
        """
        Creates a CTRIAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int]
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
        elem = CTRIAR(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquadr(self, eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                   T1=None, T2=None, T3=None, T4=None, comment=''):
        """
        Creates a CQUADR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int]
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
        elem = CQUADR(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pshell(self, pid, mid1=None, t=None, mid2=None, twelveIt3=1.0,
                   mid3=None, tst=0.833333, nsm=0.0,
                   z1=None, z2=None, mid4=None,
                   comment=''):
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
        prop = PSHELL(pid, mid1=mid1, t=t, mid2=mid2, twelveIt3=twelveIt3,
                      mid3=mid3, tst=tst, nsm=nsm,
                      z1=z1, z2=z2, mid4=mid4,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pcomp(self, pid, mids, thicknesses, thetas=None, souts=None,
                  nsm=0., sb=0., ft=None, tref=0., ge=0., lam=None,
                  z0=None, comment=''):
        """
        Creates a PCOMP card

        Parameters
        ----------
        pid : int
            property id
        mids : List[int, ..., int]
            material ids for each ply
        thicknesses : List[float, ..., float]
            thicknesses for each ply
        thetas : List[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : List[str, ..., str]; default=None
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
        prop = PCOMP(pid, mids, thicknesses, thetas, souts,
                     nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam,
                     z0=z0, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pcompg(self, pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                   nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0, lam=None, z0=None, comment=''):
        """
        Creates a PCOMPG card

        Parameters
        ----------
        pid : int
            property id
        global_ply_ids : List[int]
            the ply id
        mids : List[int, ..., int]
            material ids for each ply
        thicknesses : List[float, ..., float]
            thicknesses for each ply
        thetas : List[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : List[str, ..., str]; default=None
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
        prop = PCOMPG(pid, global_ply_ids, mids, thicknesses, thetas=thetas, souts=souts,
                      nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam, z0=z0,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pcomps(self, pid, global_ply_ids, mids, thicknesses, thetas,
                   cordm=0, psdir=13, sb=None, nb=None, tref=0.0, ge=0.0,
                   failure_theories=None, interlaminar_failure_theories=None,
                   souts=None, comment=''):
        """Creates a PCOMPS card"""
        prop = PCOMPS(pid, global_ply_ids, mids, thicknesses, thetas,
                      cordm, psdir, sb, nb, tref, ge,
                      failure_theories, interlaminar_failure_theories, souts,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_plplane(self, pid, mid, cid=0, stress_strain_output_location='GRID',
                    comment=''):
        """Creates a PLPLANE card"""
        prop = PLPLANE(pid, mid, cid=cid,
                       stress_strain_output_location=stress_strain_output_location,
                       comment=comment)
        self._add_property_object(prop)
        return prop

    def add_pplane(self, pid, mid, t=0., nsm=0., formulation_option=0,
                   comment=''):
        """Creates a PPLANE card"""
        prop = PPLANE(pid, mid, t=t, nsm=nsm, formulation_option=formulation_option,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_cplstn3(self, eid, pid, nids, theta=0.0, comment=''):
        """Creates a CPLSTN4 card"""
        elem = CPLSTN3(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cplstn4(self, eid, pid, nids, theta=0.0, comment=''):
        """Creates a CPLSTN4 card"""
        elem = CPLSTN4(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cplstn6(self, eid, pid, nids, theta=0.0, comment=''):
        """Creates a CPLSTN6 card"""
        elem = CPLSTN6(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cplstn8(self, eid, pid, nids, theta=0.0, comment=''):
        """Creates a CPLSTN8 card"""
        elem = CPLSTN8(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cplsts3(self, eid, pid, nids, theta=0.0, comment=''):
        """Creates a CPLSTS3 card"""
        elem = CPLSTS3(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ctetra(self, eid, pid, nids, comment=''):
        """
        Creates a CTETRA4/CTETRA10

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=4 or 10
        comment : str; default=''
            a comment for the card

        """
        if len(nids) == 4:
            elem = CTETRA4(eid, pid, nids, comment=comment)
        else:
            elem = CTETRA10(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cpyram(self, eid, pid, nids, comment=''):
        """
        Creates a CPYRAM5/CPYRAM13

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=5 or 13
        comment : str; default=''
            a comment for the card

        """
        if len(nids) == 5:
            elem = CPYRAM5(eid, pid, nids, comment=comment)
        else:
            elem = CPYRAM13(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cpenta(self, eid, pid, nids, comment=''):
        """
        Creates a CPENTA6/CPENTA15

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=6 or 15
        comment : str; default=''
            a comment for the card

        """
        if len(nids) == 6:
            elem = CPENTA6(eid, pid, nids, comment=comment)
        else:
            elem = CPENTA15(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_chexa(self, eid, pid, nids, comment=''):
        """
        Creates a CHEXA8/CHEXA20

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=8 or 20
        comment : str; default=''
            a comment for the card

        """
        if len(nids) == 8:
            elem = CHEXA8(eid, pid, nids, comment=comment)
        else:
            elem = CHEXA20(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_psolid(self, pid, mid, cordm=0, integ=None, stress=None, isop=None,
                   fctn='SMECH', comment=''):
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
        prop = PSOLID(pid, mid, cordm=cordm, integ=integ, stress=stress, isop=isop,
                      fctn=fctn, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_plsolid(self, pid, mid, stress_strain='GRID', ge=0., comment=''):
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
        prop = PLSOLID(pid, mid, stress_strain=stress_strain, ge=ge, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_crac2d(self, eid, pid, nids, comment=''):
        """Creates a PRAC2D card"""
        elem = CRAC2D(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_prac2d(self, pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180.,
                   comment=''):
        """Creates a PRAC2D card"""
        prop = PRAC2D(pid, mid, thick, iplane, nsm=nsm, gamma=gamma, phi=phi,
                      comment=comment)
        self._add_property_object(prop)
        return prop

    def add_crac3d(self, eid, pid, nids, comment=''):
        """Creates a CRAC3D card"""
        elem = CRAC3D(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_prac3d(self, pid, mid, gamma=0.5, phi=180., comment=''):
        """Creates a PRAC3D card"""
        prop = PRAC3D(pid, mid, gamma=gamma, phi=phi, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_axic(self, nharmonics, comment=''):
        """Creates a AXIC card"""
        axic = AXIC(nharmonics, comment=comment)
        self._add_axic_object(axic)
        return axic

    def add_pointax(self, nid, ringax, phi, comment=''):
        """Creates a POINTAX card"""
        node = POINTAX(nid, ringax, phi, comment=comment)
        self._add_ringax_object(node)
        return node

    def add_ringax(self, nid, R, z, ps=None, comment=''):
        """Creates a RINGAX card"""
        node = RINGAX(nid, R, z, ps=ps, comment=comment)
        self._add_ringax_object(node)
        return node

    def add_cconeax(self, eid, pid, rings, comment=''):
        """Creates a CCONEAX card"""
        elem = CCONEAX(eid, pid, rings, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pconeax(self, pid, mid1, t1=None, mid2=0, i=None, mid3=None, t2=None,
                    nsm=None, z1=None, z2=None, phi=None, comment=''):
        """Creates a PCONEAX card"""
        prop = PCONEAX(pid, mid1, t1, mid2, i, mid3, t2, nsm, z1, z2, phi,
                       comment=comment)
        self._add_property_object(prop)
        return prop

    def add_tempax(self, sid, ring, phi, temperature, comment=''):
        """Creates a TEMPAX card"""
        load = TEMPAX(sid, ring, phi, temperature, comment=comment)
        self._add_load_object(load)
        return load

    def add_presax(self, sid, pressure, rid1, rid2, phi1=0., phi2=360., comment=''):
        """Creates a PRESAX card"""
        load = PRESAX(sid, pressure, rid1, rid2, phi1, phi2, comment=comment)
        self._add_load_object(load)
        return load

    def add_ctrax3(self, eid, pid, nids, theta=0., comment=''):
        """Creates a CTRAX3 card"""
        elem = CTRAX3(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ctrax6(self, eid, pid, nids, theta=0., comment=''):
        """Creates a CTRAX6 card"""
        elem = CTRAX6(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ctriax(self, eid, pid, nids, theta_mcid=0., comment=''):
        """Creates a CTRIAX card"""
        elem = CTRIAX(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_ctriax6(self, eid, mid, nids, theta=0., comment=''):
        """Creates a CTRIAX6 card"""
        elem = CTRIAX6(eid, mid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquadx(self, eid, pid, nids, theta_mcid=0., comment=''):
        """Creates a CQUADX card"""
        elem = CQUADX(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquadx4(self, eid, pid, nids, theta=0., comment=''):
        """Creates a CQUADX4 card"""
        elem = CQUADX4(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cquadx8(self, eid, pid, nids, theta=0., comment=''):
        """Creates a CQUADX8 card"""
        elem = CQUADX8(eid, pid, nids, theta=theta, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cihex1(self, eid, pid, nids, comment=''):
        """see CHEXA"""
        elem = CIHEX1(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_cihex2(self, eid, pid, nids, comment=''):
        """see CHEXA"""
        elem = CIHEX2(eid, pid, nids, comment=comment)
        self._add_element_object(elem)
        return elem

    def add_pihex(self, pid, mid, cordm=0, integ=None, stress=None, isop=None,
                  fctn='SMECH', comment=''):
        """.. seealso:: PSOLID"""
        prop = PIHEX(pid, mid, cordm=cordm, integ=integ, stress=stress, isop=isop,
                     fctn=fctn, comment=comment)
        self._add_property_object(prop)
        return prop

    def add_creep(self, mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                  a, b, c, d, e, f, g, comment=''):
        """Creates a CREEP card"""
        mat = CREEP(mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                    a, b, c, d, e, f, g, comment=comment)
        self._add_creep_material_object(mat)
        return mat

    def add_mat1(self, mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0, St=0.0,
                 Sc=0.0, Ss=0.0, mcsid=0, comment=''):
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
        mat = MAT1(mid, E, G, nu, rho=rho, a=a, tref=tref, ge=ge, St=St,
                   Sc=Sc, Ss=Ss, mcsid=mcsid, comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat2(self, mid, G11, G12, G13, G22, G23, G33, rho=0.,
                 a1=None, a2=None, a3=None,
                 tref=0., ge=0., St=None, Sc=None,
                 Ss=None, mcsid=None, comment=''):
        """Creates an MAT2 card"""
        mat = MAT2(mid, G11, G12, G13, G22, G23, G33,
                   rho, a1, a2, a3,
                   tref=tref, ge=ge, St=St, Sc=Sc,
                   Ss=Ss, mcsid=mcsid, comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat3(self, mid, ex, eth, ez, nuxth, nuthz, nuzx, rho=0.0, gzx=None,
                 ax=0., ath=0., az=0., tref=0., ge=0.,
                 comment=''):
        """Creates a MAT3 card"""
        mat = MAT3(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho=rho, gzx=gzx,
                   ax=ax, ath=ath, az=az, tref=tref, ge=ge,
                   comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat4(self, mid, k, cp=0.0, rho=1.0, H=None, mu=None, hgen=1.0,
                 ref_enthalpy=None, tch=None, tdelta=None,
                 qlat=None, comment=''):
        """Creates a MAT4 card"""
        mat = MAT4(mid, k, cp=cp, rho=rho, H=H, mu=mu, hgen=hgen,
                   ref_enthalpy=ref_enthalpy, tch=tch, tdelta=tdelta,
                   qlat=qlat, comment=comment)
        self._add_thermal_material_object(mat)
        return mat

    def add_mat5(self, mid, kxx=0., kxy=0., kxz=0., kyy=0., kyz=0., kzz=0., cp=0.,
                 rho=1., hgen=1., comment=''):
        """Creates a MAT5 card"""
        mat = MAT5(mid, kxx=kxx, kxy=kxy, kxz=kxz, kyy=kyy, kyz=kyz, kzz=kzz, cp=cp,
                   rho=rho, hgen=hgen, comment=comment)
        self._add_thermal_material_object(mat)
        return mat

    def add_mat8(self, mid, e11, e22, nu12, g12=0.0, g1z=1e8, g2z=1e8,
                 rho=0., a1=0., a2=0.,
                 tref=0., Xt=0., Xc=None, Yt=0., Yc=None,
                 S=0., ge=0., F12=0., strn=0., comment=''):
        """Creates a MAT8 card"""
        mat = MAT8(mid, e11, e22, nu12, g12, g1z, g2z, rho=rho, a1=a1, a2=a2,
                   tref=tref, Xt=Xt, Xc=Xc, Yt=Yt, Yc=Yc,
                   S=S, ge=ge, F12=F12, strn=strn, comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat9(self, mid,
                 G11=0., G12=0., G13=0., G14=0., G15=0., G16=0.,
                 G22=0., G23=0., G24=0., G25=0., G26=0.,
                 G33=0., G34=0., G35=0., G36=0.,
                 G44=0., G45=0., G46=0.,
                 G55=0., G56=0., G66=0.,
                 rho=0., A=None, tref=0., ge=0., comment=''):
        """Creates a MAT9 card"""
        mat = MAT9(mid, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26,
                   G33, G34, G35, G36, G44, G45, G46, G55,
                   G56, G66, rho, A, tref, ge, comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat10(self, mid, bulk, rho, c, ge=0.0, gamma=None,
                  table_bulk=None, table_rho=None, table_ge=None, table_gamma=None,
                  comment=''):
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
        mat = MAT10(mid, bulk, rho, c, ge=ge, gamma=gamma,
                    table_bulk=table_bulk, table_rho=table_rho,
                    table_ge=table_ge, table_gamma=table_gamma,
                    comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat11(self, mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=0.0,
                  a1=0.0, a2=0.0, a3=0.0, tref=0.0, ge=0.0, comment=''):
        """Creates a MAT11 card"""
        mat = MAT11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=rho,
                    a1=a1, a2=a2, a3=a3, tref=tref, ge=ge, comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mat3d(self, mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=0.0,
                  comment=''):
        """
        This is a VABS specific card that is almost identical to the MAT11.
        """
        mat = MAT3D(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho=rho,
                    comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_matg(self, mid, idmem, behav, tabld, tablu, yprs, epl, gpl, gap=0.,
                 tab_yprs=None, tab_epl=None,
                 tab_gpl=None, tab_gap=None, comment=''):
        """Creates a MATG card"""
        mat = MATG(mid, idmem, behav, tabld, tablu, yprs, epl, gpl, gap=gap,
                   tab_yprs=tab_yprs, tab_epl=tab_epl,
                   tab_gpl=tab_gpl, tab_gap=tab_gap, comment=comment)
        self._add_structural_material_object(mat)
        return mat

    def add_mathe(self, mid, model, bulk, rho, texp, mus, alphas, betas, mooney,
                  sussbat, aboyce, comment=''):
        """Creates a MATHE card"""
        mat = MATHE(mid, model, bulk, rho, texp, mus, alphas, betas, mooney,
                    sussbat, aboyce, comment=comment)
        self._add_hyperelastic_material_object(mat)
        return mat

    def add_mathp(self, mid, a10=0., a01=0., d1=None, rho=0., av=0., tref=0., ge=0., na=1, nd=1,
                  a20=0., a11=0., a02=0., d2=0.,
                  a30=0., a21=0., a12=0., a03=0., d3=0.,
                  a40=0., a31=0., a22=0., a13=0., a04=0., d4=0.,
                  a50=0., a41=0., a32=0., a23=0., a14=0., a05=0., d5=0.,
                  tab1=None, tab2=None, tab3=None, tab4=None, tabd=None, comment=''):
        """Creates a MATHP card"""
        mat = MATHP(mid, a10, a01, d1, rho, av, tref, ge, na, nd,
                    a20, a11, a02, d2,
                    a30, a21, a12, a03, d3,
                    a40, a31, a22, a13, a04,
                    d4, a50, a41, a32, a23, a14, a05, d5, tab1, tab2, tab3,
                    tab4, tabd, comment=comment)
        self._add_hyperelastic_material_object(mat)
        return mat

    def add_mats1(self, mid, tid, Type, h, hr, yf, limit1, limit2, comment=''):
        """Creates a MATS1 card"""
        mat = MATS1(mid, tid, Type, h, hr, yf, limit1, limit2, comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_matt1(self, mid, e_table=None, g_table=None, nu_table=None, rho_table=None,
                  a_table=None, ge_table=None, st_table=None, sc_table=None, ss_table=None,
                  comment=''):
        """Creates a MATT1 card"""
        mat = MATT1(mid, e_table, g_table, nu_table, rho_table, a_table,
                    ge_table, st_table, sc_table, ss_table,
                    comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_matt2(self, mid, g11_table=None, g12_table=None, g13_table=None, g22_table=None,
                  g23_table=None, g33_table=None, rho_table=None,
                  a1_table=None, a2_table=None, a3_table=None, ge_table=None,
                  st_table=None, sc_table=None, ss_table=None,
                  comment=''):
        """Creates a MATT2 card"""
        mat = MATT2(mid, g11_table, g12_table, g13_table, g22_table,
                    g23_table, g33_table, rho_table,
                    a1_table, a2_table, a3_table, ge_table,
                    st_table, sc_table, ss_table,
                    comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_matt3(self, mid, ex_table=None, eth_table=None, ez_table=None,
                  nuth_table=None, nuxz_table=None, rho_table=None,
                  gzx_table=None, ax_table=None, ath_table=None, az_table=None,
                  ge_table=None, comment=''):
        """Creates a MATT3 card"""
        mat = MATT3(mid, ex_table, eth_table, ez_table,
                    nuth_table, nuxz_table, rho_table, gzx_table,
                    ax_table, ath_table, az_table, ge_table, comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_matt4(self, mid, k_table=None, cp_table=None, h_table=None,
                  mu_table=None, hgen_table=None, comment=''):
        """Creates a MATT4 card"""
        mat = MATT4(mid, k_table, cp_table, h_table, mu_table, hgen_table,
                    comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_matt5(self, mid, kxx_table=None, kxy_table=None, kxz_table=None,
                  kyy_table=None, kyz_table=None, kzz_table=None, cp_table=None,
                  hgen_table=None, comment=''):
        """Creates a MATT5 card"""
        mat = MATT5(mid, kxx_table, kxy_table, kxz_table, kyy_table,
                    kyz_table, kzz_table, cp_table,
                    hgen_table, comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_matt8(self, mid, e1_table=None, e2_table=None, nu12_table=None,
                  g12_table=None, g1z_table=None, g2z_table=None, rho_table=None,
                  a1_table=None, a2_table=None,
                  xt_table=None, xc_table=None, yt_table=None, yc_table=None,
                  s_table=None, ge_table=None, f12_table=None, comment=''):
        """Creates a MATT8 card"""
        mat = MATT8(mid, e1_table=e1_table, e2_table=e2_table, nu12_table=nu12_table,
                    g12_table=g12_table, g1z_table=g1z_table, g2z_table=g2z_table,
                    rho_table=rho_table, a1_table=a1_table, a2_table=a2_table,
                    xt_table=xt_table, xc_table=xc_table, yt_table=yt_table, yc_table=yc_table,
                    s_table=s_table, ge_table=ge_table, f12_table=f12_table, comment=comment)
        self._add_material_dependence_object(mat)
        return mat

    def add_nxstrat(self, sid, params, comment=''):
        """Creates an NXSTRAT card"""
        nxstrat = NXSTRAT(sid, params, comment=comment)
        self._add_nxstrat_object(nxstrat)
        return nxstrat

    def add_load(self, sid, scale, scale_factors, load_ids, comment=''):
        """
        Creates a LOAD card

        Parameters
        ----------
        sid : int
            load id
        scale : float
            overall scale factor
        scale_factors : List[float]
            individual scale factors (corresponds to load_ids)
        load_ids : List[int]
            individual load_ids (corresponds to scale_factors)
        comment : str; default=''
            a comment for the card

        """
        load = LOAD(sid, scale, scale_factors, load_ids, comment=comment)
        self._add_load_combination_object(load)
        return load

    def add_lseq(self, sid, excite_id, lid, tid=None, comment=''):
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
        load = LSEQ(sid, excite_id, lid, tid=tid, comment=comment)
        self._add_lseq_object(load)
        return load

    def add_sload(self, sid, nids, mags, comment=''):
        """
        Creates an SLOAD (GRID/SPOINT load)

        Parameters
        ----------
        sid : int
            load id
        nids : int; List[int]
            the GRID/SPOINT ids
        mags : float; List[float]
            the load magnitude
        comment : str; default=''
            a comment for the card

        """
        load = SLOAD(sid, nids, mags, comment=comment)
        self._add_load_object(load)
        return load

    def add_dload(self, sid, scale, scale_factors, load_ids, comment=''):
        """
        Creates a DLOAD card

        Parameters
        ----------
        sid : int
            Load set identification number. See Remarks 1. and 4. (Integer > 0)
        scale : float
            Scale factor. See Remarks 2. and 8. (Real)
        Si : List[float]
            Scale factors. See Remarks 2., 7. and 8. (Real)
        load_ids : List[int]
            Load set identification numbers of RLOAD1, RLOAD2, TLOAD1,
            TLOAD2, and ACSRCE entries. See Remarks 3. and 7. (Integer > 0)
        comment : str; default=''
            a comment for the card

        """
        load = DLOAD(sid, scale, scale_factors, load_ids, comment=comment)
        self._add_dload_object(load)
        return load

    def add_darea(self, sid, p, c, scale, comment=''):
        """
        Creates a DAREA card

        Parameters
        ----------
        sid : int
            darea id
        Pi : int
            GRID, EPOINT, SPOINT id
        c : str
            Component number. (0-6; 0-EPOINT/SPOINT; 1-6 GRID)
        scale : float
            Scale (area) factor

        """
        darea = DAREA(sid, p, c, scale, comment=comment)
        self._add_darea_object(darea)

    def add_tload1(self, sid, excite_id, tid, delay=0, Type='LOAD', us0=0.0, vs0=0.0,
                   comment=''):
        """
        Creates a TLOAD1 card, which defienes a load based on a table

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        tid : int
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
        load = TLOAD1(sid, excite_id, tid, delay=delay, Type=Type, us0=us0, vs0=vs0,
                      comment=comment)
        self._add_dload_entry(load)
        return load

    def add_tload2(self, sid, excite_id, delay=0, Type='LOAD', T1=0., T2=None,
                   frequency=0., phase=0., c=0., b=0.,
                   us0=0., vs0=0., comment=''):
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
        Type : int/str; default='LOAD'
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
        load = TLOAD2(sid, excite_id, delay=delay, Type=Type, T1=T1, T2=T2,
                      frequency=frequency, phase=phase, c=c, b=b,
                      us0=us0, vs0=vs0, comment=comment)
        self._add_dload_entry(load)
        return load

    def add_rload1(self, sid, excite_id, delay=0, dphase=0, tc=0, td=0,
                   Type='LOAD', comment=''):
        """
        Creates an RLOAD1 card, which defienes a frequency-dependent load
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
        load = RLOAD1(sid, excite_id, delay=delay, dphase=dphase, tc=tc, td=td,
                      Type=Type, comment=comment)
        self._add_dload_entry(load)
        return load

    def add_rload2(self, sid, excite_id, delay=0, dphase=0, tb=0, tp=0,
                   Type='LOAD', comment=''):
        """
        Creates a nRLOAD2 card, which defienes a frequency-dependent load
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
        tc : int/float; default=0
            TABLEDi id that defines C(f) for all degrees of freedom in
            EXCITEID entry
        td : int/float; default=0
            TABLEDi id that defines D(f) for all degrees of freedom in
            EXCITEID entry
        tp : int/float; default=0
            TABLEDi id that defines phi(f) for all degrees of freedom in
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
        load = RLOAD2(sid, excite_id, delay=delay, dphase=dphase, tb=tb, tp=tp,
                      Type=Type, comment=comment)
        self._add_dload_entry(load)
        return load

    def add_rforce(self, sid, nid, scale, r123, cid=0, method=1, racc=0.,
                   mb=0, idrf=0, comment=''):
        """Creates an RFORCE card"""
        load = RFORCE(sid, nid, scale, r123, cid=cid, method=method, racc=racc,
                      mb=mb, idrf=idrf, comment=comment)
        self._add_load_object(load)
        return load

    def add_rforce1(self, sid, nid, scale, group_id, cid=0, r123=None, racc=0.,
                    mb=0, method=2, comment=''):
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
        r123 : List[float, float, float] / (3, ) float ndarray
            rectangular components of the rotation vector R that passes
            through point G
        racc : int; default=0.0
            ???
        mb : int; default=0
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
        load = RFORCE1(sid, nid, scale, group_id, cid=cid, r123=r123, racc=racc,
                       mb=mb, method=method, comment=comment)
        self._add_load_object(load)
        return load

    def add_randps(self, sid, j, k, x=0., y=0., tid=0, comment=''):
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
        randps = RANDPS(sid, j, k, x=x, y=y, tid=tid, comment=comment)
        self._add_dload_entry(randps)
        return randps

    def add_randt1(self, sid, n, t0, tmax, comment=''):
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
        dload = RANDT1(sid, n, t0, tmax, comment=comment)
        self._add_dload_entry(dload)
        return dload

    def add_acsrce(self, sid, excite_id, rho, b, delay=0, dphase=0, power=0,
                   comment=''):
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
        load = ACSRCE(sid, excite_id, rho, b, delay=delay, dphase=dphase, power=power,
                      comment=comment)
        self._add_dload_entry(load)
        return load

    def add_loadcyn(self, sid, scale, segment_id, scales, load_ids,
                    segment_type=None, comment=''):
        """Creates a LOADCYN card"""
        load = LOADCYN(sid, scale, segment_id, scales, load_ids,
                       segment_type=segment_type, comment=comment)
        self._add_load_object(load)
        return load

    def add_force(self, sid, node, mag, xyz, cid=0, comment=''):
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
        load = FORCE(sid, node, mag, xyz, cid=cid, comment=comment)
        self._add_load_object(load)
        return load

    def add_force1(self, sid, node, mag, g1, g2, comment=''):
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
        load = FORCE1(sid, node, mag, g1, g2, comment=comment)
        self._add_load_object(load)
        return load

    def add_force2(self, sid, node, mag, g1, g2, g3, g4, comment=''):
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
        load = FORCE2(sid, node, mag, g1, g2, g3, g4, comment=comment)
        self._add_load_object(load)
        return load

    def add_moment(self, sid, node, mag, xyz, cid=0, comment=''):
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
        load = MOMENT(sid, node, mag, xyz, cid=cid, comment=comment)
        self._add_load_object(load)
        return load

    def add_moment1(self, sid, node, mag, g1, g2, comment=''):
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
        load = MOMENT1(sid, node, mag, g1, g2, comment=comment)
        self._add_load_object(load)
        return load

    def add_moment2(self, sid, node, mag, g1, g2, g3, g4, comment=''):
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
        load = MOMENT2(sid, node, mag, g1, g2, g3, g4, comment=comment)
        self._add_load_object(load)
        return load

    def add_deform(self, sid, eid, deformation, comment=''):
        """
        Creates an DEFORM card, which defines applied deformation on
        a 1D elemment.  Links to the DEFORM card in the case control
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
        load = DEFORM(sid, eid, deformation, comment=comment)
        self._add_load_object(load)
        return load

    def add_accel(self, sid, N, direction, locs, vals, cid=0, comment=''):
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
        locs : List[float]
            Location along direction DIR in coordinate system CID for
            specification of a load scale factor.
        vals : List[float]
            The load scale factor associated with location LOCi
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card

        """
        load = ACCEL(sid, N, direction, locs, vals, cid=cid, comment=comment)
        self._add_load_object(load)
        return load

    def add_accel1(self, sid, scale, N, nodes, cid=0, comment=''):
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
        direction : str
            Component direction of acceleration variation
            {X, Y, Z}
        nodes : List[int]
            the nodes to apply acceleration to
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card

        """
        load = ACCEL1(sid, scale, N, nodes, cid=cid, comment=comment)
        self._add_load_object(load)
        return load

    def add_grav(self, sid, scale, N, cid=0, mb=0, comment=''):
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
        grav = GRAV(sid, scale, N, cid=cid, mb=mb, comment=comment)
        self._add_load_object(grav)
        return grav

    def add_pload(self, sid, pressure, nodes, comment=''):
        """
        Creates a PLOAD card, which defines a uniform pressure load on a
        shell/solid face or arbitrarily defined quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply
        nodes : List[int]
            The nodes that are used to define the normal are defined
            using the same method as the CTRIA3/CQUAD4 normal.
            n = 3 or 4
        comment : str; default=''
            a comment for the card

        """
        load = PLOAD(sid, pressure, nodes, comment=comment)
        self._add_load_object(load)
        return load

    def add_pload1(self, sid, eid, load_type, scale, x1, p1, x2=None, p2=None, comment=''):
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
        load = PLOAD1(sid, eid, load_type, scale, x1, p1, x2=x2, p2=p2, comment=comment)
        self._add_load_object(load)
        return load

    def add_pload2(self, sid, pressure, eids, comment=''):
        """
        Creates a PLOAD2 card, which defines an applied load normal to the quad/tri face

        Parameters
        ----------
        sid : int
            load id
        pressure : float
            the pressure to apply to the elements
        eids : List[int]
            the elements to apply pressure to
            n < 6 or a continouus monotonic list of elements (e.g., [1, 2, ..., 1000])
        comment : str; default=''
            a comment for the card

        """
        load = PLOAD2(sid, pressure, eids, comment=comment)
        self._add_load_object(load)
        return load

    def add_pload4(self, sid, eids, pressures, g1=None, g34=None, cid=0,
                   nvector=None, surf_or_line='SURF', line_load_dir='NORM', comment=''):
        """
        Creates a PLOAD4 card

        Parameters
        ----------
        sid : int
            the load id
        eids : List[int, ...]
            shells : the range of element ids; must be sequential
            solids : must be length 1
        pressures : List[float, float, float, float]
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
        load = PLOAD4(sid, eids, pressures, g1=g1, g34=g34, cid=cid,
                      nvector=nvector, surf_or_line=surf_or_line,
                      line_load_dir=line_load_dir, comment=comment)
        self._add_load_object(load)
        return load

    def add_ploadx1(self, sid, eid, pa, nids, pb=None, theta=0., comment=''):
        """
        Creates a PLOADX1 card, which defines surface traction for
        axisymmetric elements.

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element id (CQUADX, CTRIAX, or CTRIAX6)
        nids : List[int, int]
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
        load = PLOADX1(sid, eid, pa, nids, pb=pb, theta=theta, comment=comment)
        self._add_load_object(load)
        return load

    def add_gmload(self, sid, normal, entity, entity_id, method, load_magnitudes,
                   cid=0, comment=''):
        """Creates a GMLOAD object"""
        load = GMLOAD(sid, normal, entity, entity_id, method, load_magnitudes,
                      cid=cid, comment=comment)
        self._add_load_object(load)
        return load

    def add_spc(self, conid, nodes, components, enforced, comment=''):
        """
        Creates an SPC card, which defines the degree of freedoms to be
        constrained

        Parameters
        ----------
        conid : int
            constraint id
        nodes : List[int]
            GRID/SPOINT ids
        components : List[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        enforced : List[float]
            the constrained value for the given node (typically 0.0)
        comment : str; default=''
            a comment for the card

        Notes
        -----
        len(nodes) == len(components) == len(enforced)

        .. warning:: non-zero enforced deflection requires an SPCD as well

        """
        spc = SPC(conid, nodes, components, enforced, comment=comment)
        self._add_constraint_spc_object(spc)
        return spc

    def add_spc1(self, conid, components, nodes, comment=''):
        """
        Creates an SPC1 card, which defines the degree of freedoms to be
        constrained to a value of 0.0

        Parameters
        ----------
        conid : int
            constraint id
        components : str
            the degree of freedoms to constrain (e.g., '1', '123')
        nodes : List[int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        spc = SPC1(conid, components, nodes, comment=comment)
        self._add_constraint_spc_object(spc)
        return spc

    def add_spcd(self, sid, nodes, components, enforced, comment=''):
        """
        Creates an SPCD card, which defines the degree of freedoms to be
        set during enforced motion

        Parameters
        ----------
        conid : int
            constraint id
        nodes : List[int]
            GRID/SPOINT ids
        components : List[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        enforced : List[float]
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
        spc = SPCD(sid, nodes, components, enforced, comment=comment)
        self._add_load_object(spc)
        return spc

    def add_spcadd(self, conid, sets, comment=''):
        """Creates a SPCADD card"""
        spcadd = SPCADD(conid, sets, comment=comment)
        self._add_constraint_spcadd_object(spcadd)
        return spcadd

    def add_spcax(self, conid, ringax, hid, component, enforced, comment=''):
        """Creates an SPCAX card"""
        spcax = SPCAX(conid, ringax, hid, component, enforced, comment=comment)
        self._add_constraint_spc_object(spcax)
        return spcax

    def add_gmspc(self, conid, component, entity, entity_id, comment=''):
        """Creates a GMSPC card"""
        spc = GMSPC(conid, component, entity, entity_id, comment=comment)
        self._add_constraint_spc_object(spc)
        return spc

    def add_mpc(self, conid, nodes, components, coefficients, comment=''):
        """
        Creates an MPC card

        Parameters
        ----------
        conid : int
            Case Control MPC id
        nodes : List[int]
            GRID/SPOINT ids
        components : List[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        coefficients : List[float]
            the scaling coefficients

        """
        mpc = MPC(conid, nodes, components, coefficients, comment=comment)
        self._add_constraint_mpc_object(mpc)
        return mpc

    def add_mpcadd(self, conid, sets, comment=''):
        """Creates an MPCADD card"""
        mpcadd = MPCADD(conid, sets, comment=comment)
        self._add_constraint_mpcadd_object(mpcadd)
        return mpcadd

    def add_suport(self, nodes, Cs, comment=''):
        """
        Creates a SUPORT card, which defines free-body reaction points.
        This is always active.

        Parameters
        ----------
        nodes : List[int]
            the nodes to release
        Cs : List[str]
            compoents to support at each node
        comment : str; default=''
            a comment for the card

        """
        suport = SUPORT(nodes, Cs, comment=comment)
        self._add_suport_object(suport)
        return suport

    def add_suport1(self, conid, nodes, Cs, comment=''):
        """
        Creates a SUPORT card, which defines free-body reaction points.

        Parameters
        ----------
        conid : int
            Case Control SUPORT id
        nodes : List[int]
            the nodes to release
        Cs : List[str]
            compoents to support at each node
        comment : str; default=''
            a comment for the card

        """
        suport1 = SUPORT1(conid, nodes, Cs, comment=comment)
        self._add_suport1_object(suport1)
        return suport1

    def add_aeros(self, cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0,
                  comment=''):
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
        self._add_aeros_object(aeros)
        return aeros

    def add_aero(self, velocity, cref, rho_ref, acsid=0, sym_xz=0, sym_xy=0,
                 comment=''):
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
        self._add_aero_object(aero)
        return aero

    def add_caero1(self, eid, pid, igroup, p1, x12, p4, x43,
                   cp=0, nspan=0, lspan=0, nchord=0, lchord=0, comment=''):
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO1
            int : PAERO1 ID
            PAERO1 : PAERO1 object (xref)
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
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card

        """
        caero = CAERO1(eid, pid, igroup, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment=comment)
        self._add_caero_object(caero)
        return caero

    def add_caero2(self, eid, pid, igroup, p1, x12, cp=0, nsb=0, nint=0, lsb=0,
                   lint=0, comment=''):
        """
        Defines a CAERO2 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO2
            int : PAERO2 ID
            PAERO2 : PAERO2 object (xref)
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (forward position)
        x12 : float
            length of the CAERO2
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
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
        caero = CAERO2(eid, pid, igroup, p1, x12, cp=cp, nsb=nsb, nint=nint, lsb=lsb,
                       lint=lint, comment=comment)
        self._add_caero_object(caero)
        return caero

    def add_caero3(self, eid, pid, list_w, p1, x12, p4, x43,
                   cp=0, list_c1=None, list_c2=None, comment=''):
        """Creates a CAERO3 card"""
        caero = CAERO3(eid, pid, list_w, p1, x12, p4, x43,
                       cp=cp, list_c1=list_c1, list_c2=list_c2, comment=comment)
        self._add_caero_object(caero)
        return caero

    def add_caero4(self, eid, pid, p1, x12, p4, x43,
                   cp=0, nspan=0, lspan=0, comment=''):
        """
        Defines a CAERO4 card, which defines a strip theory surface.

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO4
            int : PAERO4 ID
            PAERO4 : PAERO4 object (xref)
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
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card

        """
        caero = CAERO4(eid, pid, p1, x12, p4, x43,
                       cp=cp, nspan=nspan, lspan=lspan, comment=comment)
        self._add_caero_object(caero)
        return caero

    def add_caero5(self, eid, pid, p1, x12, p4, x43, cp=0, nspan=0, lspan=0,
                   ntheory=0, nthick=0, comment=''):
        """Creates a CAERO5 card"""
        caero = CAERO5(eid, pid, p1, x12, p4, x43, cp=cp, nspan=nspan, lspan=lspan,
                       ntheory=ntheory, nthick=nthick, comment=comment)
        self._add_caero_object(caero)
        return caero

    def add_paero1(self, pid, caero_body_ids=None, comment=''):
        """
        Creates a PAERO1 card, which defines associated bodies for the
        panels in the Doublet-Lattice method.

        Parameters
        ----------
        pid : int
            PAERO1 id
        caero_body_ids : List[int]; default=None
            CAERO2 ids that are within the same IGID group
        comment : str; default=''
            a comment for the card

        """
        paero = PAERO1(pid, caero_body_ids=caero_body_ids, comment=comment)
        self._add_paero_object(paero)
        return paero

    def add_paero2(self, pid, orient, width, AR, thi, thn,
                   lrsb=None, lrib=None, lth=None, comment=''):
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
        thi / thn : List[int]
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
        lth : List[int, int]; default=None
            AEFACT ids for defining theta arrays for interference calculations
            for theta1/theta2; length=2
        comment : str; default=''
            a comment for the card

        """
        #paero = PAERO2(pid, orient, width, AR, thi, thn,
                       #lrsb, lrib, lth1, lth2, comment=comment)
        paero = PAERO2(pid, orient, width, AR, thi, thn, lrsb=lrsb, lrib=lrib,
                       lth=lth, comment=comment)
        self._add_paero_object(paero)
        return paero

    def add_paero3(self, pid, nbox, ncontrol_surfaces, x, y, comment=''):
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
        x / y : List[float, None]
            float : locations of points 5 through 12, which are in the
            aerodynamic coordinate system, to define the cranks and
            control surface geometry.
        comment : str; default=''
            a comment for the card

        """
        paero = PAERO3(pid, nbox, ncontrol_surfaces, x, y, comment=comment)
        self._add_paero_object(paero)
        return paero

    def add_paero4(self, pid, docs, caocs, gapocs, cla=0, lcla=0,
                   circ=0, lcirc=0, comment=''):
        """Creates a PAERO4 card"""
        paero = PAERO4(pid, docs, caocs, gapocs, cla=cla, lcla=lcla,
                       circ=circ, lcirc=lcirc, comment=comment)
        self._add_paero_object(paero)
        return paero

    def add_paero5(self, pid, caoci, nalpha=0, lalpha=0, nxis=0, lxis=0,
                   ntaus=0, ltaus=0, comment=''):
        """Creates a PAERO5 card"""
        paero = PAERO5(pid, caoci, nalpha=nalpha, lalpha=lalpha, nxis=nxis, lxis=lxis,
                       ntaus=ntaus, ltaus=ltaus, comment=comment)
        self._add_paero_object(paero)
        return paero

    def add_spline1(self, eid, caero, box1, box2, setg, dz=0., method='IPS',
                    usage='BOTH', nelements=10,
                    melements=10, comment=''):
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
        spline = SPLINE1(eid, caero, box1, box2, setg, dz=dz, method=method,
                         usage=usage, nelements=nelements, melements=melements,
                         comment=comment)
        self._add_spline_object(spline)
        return spline

    def add_spline2(self, eid, caero, id1, id2, setg, dz=0.0, dtor=1.0, cid=0,
                    dthx=0., dthy=0., usage='BOTH', comment=''):
        """
        Creates a SPLINE2 card, which defines a beam spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        id1 / id2 : int
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
        spline = SPLINE2(eid, caero, id1, id2, setg, dz=dz, dtor=dtor, cid=cid,
                         dthx=dthx, dthy=dthy, usage=usage, comment=comment)
        self._add_spline_object(spline)
        return spline

    def add_spline3(self, eid, caero, box_id, components, nodes,
                    displacement_components,
                    coeffs, usage='BOTH', comment=''):
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

        nodes : List[int]
           Grid point identification number of the independent grid point.
        displacement_components :  : List[int]
           Component numbers in the displacement coordinate system.
           1-6 (GRIDs)
           0 (SPOINTs)
        coeffs :  : List[float]
           Coefficient of the constraint relationship.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        spline = SPLINE3(eid, caero, box_id, components, nodes,
                         displacement_components,
                         coeffs, usage=usage,
                         comment=comment)
        self._add_spline_object(spline)
        return spline

    def add_spline4(self, eid, caero, aelist, setg, dz, method, usage,
                    nelements, melements, comment=''):
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
        nelements : int; default=10
            The number of FE elements along the local spline x-axis if
            using the FPS option
        melements : int; default=10
            The number of FE elements along the local spline y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        spline = SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                         nelements, melements, comment=comment)
        self._add_spline_object(spline)
        return spline

    def add_spline5(self, eid, caero, aelist, setg, thx, thy, dz=0., dtor=1.0, cid=0,
                    usage='BOTH', method='BEAM', ftype='WF2', rcore=None, comment=''):
        """Creates a SPLINE5 card"""
        assert isinstance(cid, int), cid
        spline = SPLINE5(eid, caero, aelist, setg, thx, thy, dz=dz, dtor=dtor, cid=cid,
                         usage=usage, method=method, ftype=ftype, rcore=rcore, comment=comment)
        self._add_spline_object(spline)
        return spline

    def add_trim(self, sid, mach, q, labels, uxs, aeqr=1.0, trim_type=1, comment=''):
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
        labels : List[str]
            names of the fixed variables
        uxs : List[float]
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
            trim = TRIM(sid, mach, q, labels, uxs, aeqr=aeqr, comment=comment)
        elif trim_type == 2:
            trim = TRIM2(sid, mach, q, labels, uxs, aeqr=aeqr, comment=comment)
        else:  # pragma: no cover
            raise ValueError(trim_type)
        self._add_trim_object(trim)
        return trim

    def add_mkaero1(self, machs, reduced_freqs, comment=''):
        """
        Creates an MKAERO1 card, which defines a set of mach and
        reduced frequencies.

        Parameters
        ----------
        machs : List[float]
            series of Mach numbers
        reduced_freqs : List[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card

        """
        mkaero = MKAERO1(machs, reduced_freqs, comment=comment)
        self._add_mkaero_object(mkaero)
        return mkaero

    def add_mkaero2(self, machs, reduced_freqs, comment=''):
        """
        Creates an MKAERO2 card, which defines a set of mach and
        reduced frequency pairs.

        Parameters
        ----------
        machs : List[float]
            series of Mach numbers
        reduced_freqs : List[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card

        """
        mkaero = MKAERO2(machs, reduced_freqs, comment=comment)
        self._add_mkaero_object(mkaero)
        return mkaero

    def add_gust(self, sid, dload, wg, x0, V=None, comment=''):
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
        gust = GUST(sid, dload, wg, x0, V=V, comment=comment)
        self._add_gust_object(gust)
        return gust

    def add_eigr(self, sid, method='LAN', f1=None, f2=None, ne=None, nd=None,
                 norm='MASS', G=None, C=None, comment=''):
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
        self._add_method_object(method)
        return method

    def add_eigrl(self, sid, v1=None, v2=None, nd=None, msglvl=0, maxset=None, shfscl=None,
                  norm=None, options=None, values=None, comment=''):
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
        self._add_method_object(method)
        return method

    def add_eigb(self, sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                 comment=''):
        """Creates an EIGB card"""
        method = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                      comment=comment)
        self._add_method_object(method)
        return method

    def add_eigc(self, sid, method, grid, component, epsilon, neigenvalues,
                 norm='MAX', # common
                 mblkszs=None, iblkszs=None, ksteps=None, NJIs=None,
                 alphaAjs=None, omegaAjs=None,
                 alphaBjs=None, omegaBjs=None,
                 LJs=None, NEJs=None, NDJs=None,
                 shift_r1=None, shift_i1=None, isrr_flag=None,
                 nd1=None, comment=''):
        """Creates an EIGC card"""
        method = EIGC(sid, method, grid, component, epsilon, neigenvalues, norm=norm,
                      mblkszs=mblkszs, iblkszs=iblkszs, ksteps=ksteps, NJIs=NJIs,

                      alphaAjs=alphaAjs, omegaAjs=omegaAjs,
                      alphaBjs=alphaBjs, omegaBjs=omegaBjs,

                      LJs=LJs, NEJs=NEJs, NDJs=NDJs,
                      shift_r1=shift_r1, shift_i1=shift_i1, isrr_flag=isrr_flag,
                      nd1=nd1, comment=comment)
        self._add_cmethod_object(method)
        return method

    def add_eigp(self, sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=''):
        """Creates an EIGP card"""
        method = EIGP(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=comment)
        self._add_cmethod_object(method)
        return method

    def add_set1(self, sid, ids, is_skin=False, comment=''):
        """
        Creates a SET1 card, which defines a list of structural grid
        points or element identification numbers.

        Parameters
        ----------
        sid : int
            set id
        ids : List[int, str]
            AECOMP, SPLINEx, PANEL : all grid points must exist
            XYOUTPUT : missing grid points are ignored
            The only valid string is THRU
            ``ids = [1, 3, 5, THRU, 10]``
        is_skin : bool; default=False
            if is_skin is used; ids must be empty
        comment : str; default=''
            a comment for the card

        """
        set_obj = SET1(sid, ids, is_skin=is_skin, comment=comment)
        self._add_set_object(set_obj)
        return set_obj

    def add_set3(self, sid, desc, ids, comment=''):
        """Creates a SET3 card"""
        set_obj = SET3(sid, desc, ids, comment=comment)
        self._add_set_object(set_obj)
        return set_obj

    def add_aset(self, ids, components, comment=''):
        """
        Creates an ASET/ASET1 card, which defines the degree of freedoms
        that will be retained during an ASET modal reduction.

        Parameters
        ----------
        ids : List[int]
            the GRID/SPOINT ids
        components : List[str]; str
            the degree of freedoms to be retained (e.g., '1', '123')
            if a list is passed in, a ASET is made
            if a str is passed in, a ASET1 is made
        comment : str; default=''
            a comment for the card

        Notes
        -----
        the length of components and ids must be the same
        """
        if isinstance(components, string_types):
            aset = ASET1(ids, components, comment=comment)
        else:
            aset = ASET(ids, components, comment=comment)
        self._add_aset_object(aset)
        return aset

    def add_aset1(self, ids, components, comment=''):
        """.. .. seealso:: ``add_aset``"""
        return self.add_aset(ids, components, comment=comment)

    def add_bset(self, ids, components, comment=''):
        """
        Creates an BSET/BSET1 card, which defines the degree of freedoms
        that will be fixed during a generalized dynamic reduction or
        component model synthesis calculation.

        Parameters
        ----------
        ids : List[int]
            the GRID/SPOINT ids
        components : List[str]; str
            the degree of freedoms to be fixed (e.g., '1', '123')
            if a list is passed in, a ASET is made
            if a str is passed in, a ASET1 is made
        comment : str; default=''
            a comment for the card

        Notes
        -----
        the length of components and ids must be the same
        """
        if isinstance(components, string_types):
            bset = BSET1(ids, components, comment=comment)
        else:
            bset = BSET(ids, components, comment=comment)
        self._add_bset_object(bset)
        return bset

    def add_bset1(self, ids, components, comment=''):
        """.. .. seealso:: ``add_bset``"""
        return self.add_bset(ids, components, comment=comment)

    def add_cset(self, ids, components, comment=''):
        """
        Creates an CSET/CSET1 card, which defines the degree of freedoms
        that will be free during a generalized dynamic reduction or
        component model synthesis calculation.

        Parameters
        ----------
        ids : List[int]
            the GRID/SPOINT ids
        components : List[str]; str
            the degree of freedoms to be free (e.g., '1', '123')
            if a list is passed in, a CSET is made
            if a str is passed in, a CSET1 is made
        comment : str; default=''
            a comment for the card

        Notes
        -----
        the length of components and ids must be the same

        """
        if isinstance(components, string_types):
            cset = CSET1(ids, components, comment=comment)
        else:
            cset = CSET(ids, components, comment=comment)
        self._add_cset_object(cset)
        return cset

    def add_cset1(self, ids, components, comment=''):
        """.. seealso:: ``add_cset``"""
        return self.add_cset(ids, components, comment=comment)

    #def add_omit1(self, ids, components, comment=''):
        #""".. seealso:: ``add_omit``"""
        #return self.add_omit(ids, components, comment=comment)

    #def add_omit(self, ids, components, comment=''):
        #"""
        #Creates an OMIT1 card, which defines the degree of freedoms that
        #will be excluded (o-set) from the analysis set (a-set).

        #Parameters
        #----------
        #ids : List[int]
            #the GRID/SPOINT ids
        #components : List[str]; str
            #the degree of freedoms to be retained (e.g., '1', '123')
            #if a list is passed in, a OMIT is made
            #if a str is passed in, a OMIT1 is made
        #comment : str; default=''
            #a comment for the card
        #"""
        #if isinstance(components, string_types):
            #omit = OMIT1(ids, components, comment=comment)
        #else:
            #omit = OMIT(ids, components, comment=comment)
        #self._add_omit_object(omit)
        #return omit

    def add_omit1(self, ids, components, comment=''):
        """
        Creates an OMIT1 card, which defines the degree of freedoms that
        will be excluded (o-set) from the analysis set (a-set).

        Parameters
        ----------
        ids : List[int]
            the GRID/SPOINT ids
        components : str
            the degree of freedoms to be omitted (e.g., '1', '123')
        comment : str; default=''
            a comment for the card

        """
        omit1 = OMIT1(ids, components, comment=comment)
        self._add_omit_object(omit1)
        return omit1

    def add_qset(self, ids, components, comment=''):
        """
        Creates a QSET/QSET1 card, which defines generalized degrees of
        freedom (q-set) to be used for dynamic reduction or component
        mode synthesis.

        Parameters
        ----------
        ids : List[int]
            the GRID/SPOINT ids
        components : List[str]; str
            the degree of freedoms to be created (e.g., '1', '123')
            if a list is passed in, a QSET is made
            if a str is passed in, a QSET1 is made
        comment : str; default=''
            a comment for the card

        """
        if isinstance(components, string_types):
            qset = QSET1(ids, components, comment=comment)
        else:
            qset = QSET(ids, components, comment=comment)
        self._add_qset_object(qset)
        return qset

    def add_qset1(self, ids, components, comment=''):
        """.. seealso:: ``add_qset``"""
        return self.add_qset(ids, components, comment=comment)

    def add_uset(self, name, ids, components, comment=''):
        """
        Creates a USET card, which defines a degrees-of-freedom set.

        Parameters
        ----------
        name : str
            SNAME Set name. (One to four characters or the word 'ZERO'
            followed by the set name.)
        ids : List[int]
            the GRID/SPOINT ids
        components : List[str]
            the degree of freedoms (e.g., '1', '123')
            if a list is passed in, a USET is made
            if a str is passed in, a USET1 is made
        comment : str; default=''
            a comment for the card

        """
        if isinstance(components, integer_string_types):
            uset = USET1(name, ids, components, comment=comment)
        else:
            uset = USET(name, ids, components, comment=comment)
        self._add_uset_object(uset)
        return uset

    def add_uset1(self, name, ids, components, comment=''):
        """.. seealso:: ``add_uset``"""
        return self.add_uset(name, ids, components, comment=comment)

    def add_sebset(self, seid, ids, components, comment=''):
        """Creates an SEBSET/SEBSET1 card"""
        if isinstance(components, integer_string_types):
            sebset = SEBSET1(seid, ids, components, comment=comment)
        else:
            sebset = SEBSET(seid, ids, components, comment=comment)
        self._add_sebset_object(sebset)
        return sebset

    def add_sebset1(self, seid, ids, components, comment=''):
        """.. seealso:: ``add_secset``"""
        return self.add_sebset(seid, ids, components, comment=comment)

    def add_secset(self, seid, ids, components, comment=''):
        """Creates an SECSET/SECSET1 card"""
        if isinstance(components, integer_string_types):
            secset = SECSET1(seid, ids, components, comment=comment)
        else:
            secset = SECSET(seid, ids, components, comment=comment)
        self._add_secset_object(secset)
        return secset

    def add_secset1(self, seid, ids, components, comment=''):
        """.. seealso:: ``add_secset``"""
        return self.add_secset(seid, ids, components, comment=comment)

    def add_seqset(self, seid, ids, components, comment=''):
        """Creates an SEQSET card"""
        if isinstance(components, integer_string_types):
            seqset = SEQSET1(seid, ids, components, comment=comment)
        else:
            seqset = SEQSET(seid, ids, components, comment=comment)
        self._add_seqset_object(seqset)
        return seqset

    def add_seqset1(self, seid, ids, components, comment=''):
        """.. seealso:: ``add_secset``"""
        return self.add_seqset(seid, ids, components, comment=comment)

    def add_seset(self, seid, ids, comment=''):
        """Creates an SEUSET card"""
        seset = SESET(seid, ids, comment=comment)
        self._add_seset_object(seset)
        return seset

    def add_sesup(self, nodes, Cs, comment=''):
        """Creates an SESUP card"""
        se_suport = SESUP(nodes, Cs, comment=comment)
        self._add_sesuport_object(se_suport)
        return se_suport

    def add_flutter(self, sid, method, density, mach, reduced_freq_velocity,
                    imethod='L', nvalue=None,
                    omax=None, epsilon=1.0e-3, comment=''):
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
            For the PKS and PKNLS methods, OMAX specifies the maximum frequency, in
            Hz., to be used in he flutter sweep.
            MSC only.
        epsilon : float; default=1.0e-3
            Convergence parameter for k. Used in the PK and PKNL methods only
        comment : str; default=''
            a comment for the card

        """
        flutter = FLUTTER(sid, method, density, mach, reduced_freq_velocity,
                          imethod=imethod, nvalue=nvalue,
                          omax=omax, epsilon=epsilon,
                          comment=comment)
        self._add_flutter_object(flutter)
        return flutter

    def add_flfact(self, sid, factors, comment=''):
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
            values : List[float, ..., float]
                list of factors
            List[f1, THRU, fnf, nf, fmid]
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
        flfact = FLFACT(sid, factors, comment=comment)
        self._add_flfact_object(flfact)
        return flfact

    def add_aecomp(self, name, list_type, lists, comment=''):
        # type: (str, List[str], Union[int, List[int]], str) -> None
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
        lists : List[int, int, ...]; int
            The identification number of either SET1, AELIST or CAEROi
            entries that define the set of grid points that comprise
            the component
        comment : str; default=''
            a comment for the card

        """
        aecomp = AECOMP(name, list_type, lists, comment=comment)
        self._add_aecomp_object(aecomp)
        return aecomp

    def add_aestat(self, aestat_id, label, comment=''):
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
        aestat = AESTAT(aestat_id, label, comment=comment)
        self._add_aestat_object(aestat)

    def add_aelink(self, aelink_id, label, independent_labels, linking_coefficents, comment=''):
        """
        Creates an AELINK card, which defines an equation linking
        AESTAT and AESURF cards

        Parameters
        ----------
        aelink_id : int
            unique id
        label : str
            name of the dependent AESURF card
        independent_labels : List[str, ..., str]
            name for the independent variables (AESTATs)
        linking_coefficents : List[float]
            linking coefficients
        comment : str; default=''
            a comment for the card

        """
        aelink = AELINK(aelink_id, label, independent_labels, linking_coefficents, comment=comment)
        self._add_aelink_object(aelink)
        return aelink

    def add_aelist(self, sid, elements, comment=''):
        """
        Creates an AELIST card, which defines the aero boxes for
        an AESURF/SPLINEx.

        Parameters
        ----------
        sid : int
            unique id
        elements : List[int, ..., int]
            list of box ids
        comment : str; default=''
            a comment for the card

        """
        aelist = AELIST(sid, elements, comment=comment)
        self._add_aelist_object(aelist)
        return aelist

    def add_aefact(self, sid, fractions, comment=''):
        """
        Creates an AEFACT card, which is used by the CAEROx / PAEROx card
        to adjust the spacing of the sub-paneleing (and grid point
        paneling in the case of the CAERO3).

        Parameters
        ----------
        sid : int
            unique id
        fractions : List[float, ..., float]
            list of percentages
        comment : str; default=''
            a comment for the card

        """
        aefact = AEFACT(sid, fractions, comment=comment)
        self._add_aefact_object(aefact)
        return aefact

    def add_diverg(self, sid, nroots, machs, comment=''):
        """
        Creates an DIVERG card, which is used in divergence
        analysis (SOL 144).

        Parameters
        ----------
        sid : int
            The name
        nroots : int
            the number of roots
        machs : List[float, ..., float]
            list of Mach numbers
        comment : str; default=''
            a comment for the card

        """
        diverg = DIVERG(sid, nroots, machs, comment=comment)
        self._add_diverg_object(diverg)
        return diverg

    def add_csschd(self, sid, aesid, lschd, lalpha=None, lmach=None,
                   comment=''):
        """
        Creates an CSSCHD card, which defines a specified control surface
        deflection as a function of Mach and alpha (used in SOL 144/146).

        Parameters
        ----------
        sid : int
            the unique id
        aesid : int
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
        csschd = CSSCHD(sid, aesid, lschd, lalpha=lalpha, lmach=lmach,
                        comment=comment)
        self._add_csschd_object(csschd)
        return csschd

    def add_aesurf(self, aesid, label, cid1, alid1, cid2=None, alid2=None,
                   eff=1.0, ldw='LDW', crefc=1.0, crefs=1.0,
                   pllim=-np.pi/2., pulim=np.pi/2.,
                   hmllim=None, hmulim=None, # hinge moment lower/upper limits
                   tqllim=None, tqulim=None, # TABLEDi deflection limits vs. dynamic pressure
                   comment=''):
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
        alid1 / alid2 : int / None
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
        aesurf = AESURF(aesid, label, cid1, alid1, cid2=cid2, alid2=alid2,
                        eff=eff, ldw=ldw, crefc=crefc, crefs=crefs,
                        pllim=pllim, pulim=pulim,
                        hmllim=hmllim, hmulim=hmulim,
                        tqllim=tqllim, tqulim=tqulim, comment=comment)
        self._add_aesurf_object(aesurf)
        return aesurf

    def add_aesurfs(self, aesid, label, list1, list2, comment=''):
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
        aesurfs = AESURFS(aesid, label, list1, list2, comment=comment)
        self._add_aesurfs_object(aesurfs)
        return aesurfs

    def add_aeparm(self, aeparm_id, label, units, comment=''):
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
        aeparm = AEPARM(aeparm_id, label, units, comment=comment)
        self._add_aeparm_object(aeparm)
        return aeparm

    def add_dtable(self, default_values, comment=''):
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
        self._add_dtable_object(dtable)
        return dtable

    def add_tabled1(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment=''):
        """
        Creates a TABLED1, which is a dynamic load card that is applied
        by the DAREA card

        Parameters
        ----------
        tid : int
            table id
        x : List[float]
            nvalues
        y : List[float]
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
        table = TABLED1(tid, x, y, xaxis=xaxis, yaxis=yaxis, extrap=extrap, comment=comment)
        self._add_tabled_object(table)
        return table

    def add_tabled2(self, tid, x1, x, y, comment=''):
        """Creates a TABLED2 card"""
        table = TABLED2(tid, x1, x, y, comment=comment)
        self._add_tabled_object(table)
        return table

    def add_tabled3(self, tid, x1, x2, x, y, comment=''):
        """Creates a TABLED3 card"""
        table = TABLED3(tid, x1, x2, x, y, comment=comment)
        self._add_tabled_object(table)
        return table

    def add_tabled4(self, tid, x1, x2, x3, x4, a, comment=''):
        """Creates a TABLED4 card"""
        table = TABLED4(tid, x1, x2, x3, x4, a, comment=comment)
        self._add_tabled_object(table)
        return table

    def add_tablem1(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        """Creates a TABLEM1 card"""
        table = TABLEM1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)
        self._add_tablem_object(table)
        return table

    def add_tablem2(self, tid, x1, x, y, extrap=0, comment=''):
        """Creates a TABLEM2 card"""
        table = TABLEM2(tid, x1, x, y, extrap=extrap, comment=comment)
        self._add_tablem_object(table)
        return table

    def add_tablem3(self, tid, x1, x2, x, y, extrap=0, comment=''):
        """Creates a TABLEM3 card"""
        table = TABLEM3(tid, x1, x2, x, y, extrap=extrap, comment=comment)
        self._add_tablem_object(table)
        return table

    def add_tablem4(self, tid, x1, x2, x3, x4, a, comment=''):
        """Creates a TABLEM4 card"""
        table = TABLEM4(tid, x1, x2, x3, x4, a, comment=comment)
        self._add_tablem_object(table)
        return table

    def add_tables1(self, tid, x, y, Type=1, comment=''):
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
        x, y : List[float]
            table values
        comment : str; default=''
            a comment for the card
        """
        table = TABLES1(tid, x, y, Type=Type, comment=comment)
        self._add_table_object(table)
        return table

    def add_tablest(self, tid, x, y, comment=''):
        """Creates an TABLEST card"""
        table = TABLEST(tid, x, y, comment=comment)
        self._add_table_object(table)
        return table

    def add_tabrnd1(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        """Creates an TABRND1 card"""
        table = TABRND1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)
        self._add_random_table_object(table)
        return table

    def add_tabrndg(self, tid, Type, LU, WG, comment=''):
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
        self._add_random_table_object(table)
        return table

    def add_tabdmp1(self, tid, x, y, Type='G', comment=''):
        """Creates a TABDMP1 card"""
        table = TABDMP1(tid, x, y, Type=Type, comment=comment)
        self._add_table_sdamping_object(table)
        return table

    def add_freq(self, sid, freqs, comment=''):
        """
        Creates a FREQ card

        Parameters
        ----------
        sid : int
            set id referenced by case control FREQUENCY
        freqs : List[float]
            the frequencies for a FREQx object
        comment : str; default=''
            a comment for the card

        """
        freq = FREQ(sid, freqs, comment=comment)
        self._add_freq_object(freq)
        return freq

    def add_freq1(self, sid, f1, df, ndf=1, comment=''):
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
        self._add_freq_object(freq)
        return freq

    def add_freq2(self, sid, f1, f2, nf=1, comment=''):
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
        self._add_freq_object(freq)
        return freq

    def add_freq3(self, sid, f1, f2=None, Type='LINEAR', nef=10, cluster=1.0,
                  comment=''):
        """Creates a FREQ3 card"""
        freq = FREQ3(sid, f1, f2, Type, nef, cluster, comment)
        self._add_freq_object(freq)
        return freq

    def add_freq4(self, sid, f1=0., f2=1e20, fspread=0.1, nfm=3, comment=''):
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
        self._add_freq_object(freq)
        return freq

    def add_freq5(self, sid, fractions, f1=0., f2=1e20, comment=''):
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
        fractions : List[float]
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
        self._add_freq_object(freq)
        return freq

    def add_rrod(self, eid, nids, cma='', cmb='', alpha=0.0, comment=''):
        """
        Creates a RROD element

        Parameters
        ----------
        eid : int
            element id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        cma / cmb : str; default=''
            dependent DOFs
        alpha : float; default=0.0
            coefficient of thermal expansion
        comment : str; default=''
            a comment for the card

        """
        elem = RROD(eid, nids, cma=cma, cmb=cmb, alpha=alpha, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rbe1(self, eid, Gni, Cni, Gmi, Cmi, alpha=0., comment=''):
        """
        Creates an RBE1 element

        Parameters
        ----------
        eid : int
            element id
        Gni : List[int]
            independent node ids
        Cni : List[str]
            the independent components (e.g., '123456')
        Gmi : List[int]
            dependent node ids
        Cmi : List[str]
            the dependent components (e.g., '123456')
        alpha : float; default=0.
            thermal expansion coefficient
        comment : str; default=''
            a comment for the card

        """
        elem = RBE1(eid, Gni, Cni, Gmi, Cmi, alpha=alpha, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rbe2(self, eid, gn, cm, Gmi, alpha=0.0, comment=''):
        """Creates an RBE2 element"""
        elem = RBE2(eid, gn, cm, Gmi, alpha=alpha, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rbe3(self, eid, refgrid, refc, weights, comps, Gijs,
                 Gmi=None, Cmi=None, alpha=0.0, comment=''):
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
        GiJs : List[int, ..., int]
            independent nodes
        comps : List[str, ..., str]
            independent components
        weights : List[float, ..., float]
            independent weights for the importance of the DOF
        Gmi : List[int, ..., int]; default=None -> []
            dependent nodes / UM Set
        Cmi : List[str, ..., str]; default=None -> []
            dependent components / UM Set
        alpha : float; default=0.0
            thermal expansion coefficient
        comment : str; default=''
            a comment for the card

        """
        elem = RBE3(eid, refgrid, refc, weights, comps, Gijs,
                    Gmi=Gmi, Cmi=Cmi, alpha=alpha, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rbar(self, eid, nids, cna, cnb, cma, cmb, alpha=0., comment=''):
        """
        Creates a RBAR element

        Parameters
        ----------
        eid : int
            element id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        cna / cnb : str
            independent DOFs in '123456'
        cma / cmb : str
            dependent DOFs in '123456'
        alpha : float; default=0.0
            coefficient of thermal expansion
        comment : str; default=''
            a comment for the card

        """
        elem = RBAR(eid, nids, cna, cnb, cma, cmb, alpha=alpha, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rbar1(self, eid, nids, cb, alpha=0., comment=''):
        """Creates an RBAR1 element"""
        elem = RBAR1(eid, nids, cb, alpha=alpha, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rspline(self, eid, independent_nid, dependent_nids, dependent_components,
                    diameter_ratio=0.1, comment=''):
        """
        Creates an RSPLINE card, which uses multipoint constraints for the
        interpolation of displacements at grid points

        Parameters
        ----------
        eid : int
            element id
        independent_nid : int
            the independent node id
        dependent_nids : List[int]
            the dependent node ids
        dependent_components : List[str]
            Components to be constrained
        diameter_ratio : float; default=0.1
            Ratio of the diameter of the elastic tube to the sum of the
            lengths of all segments
        comment : str; default=''
            a comment for the card

        """
        elem = RSPLINE(eid, independent_nid, dependent_nids, dependent_components,
                       diameter_ratio=diameter_ratio, comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_rsscon(self, eid, rigid_type,
                   shell_eid=None, solid_eid=None,
                   a_solid_grids=None, b_solid_grids=None, shell_grids=None,
                   comment=''):
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
        shell/solid_grids : List[int, int]; default=None
            the shell/solid node ids (if rigid_type=GRID)
        comment : str; default=''
            a comment for the card

        """
        elem = RSSCON(eid, rigid_type,
                      shell_eid=shell_eid, solid_eid=solid_eid,
                      a_solid_grids=a_solid_grids, b_solid_grids=b_solid_grids,
                      shell_grids=shell_grids,
                      comment=comment)
        self._add_rigid_element_object(elem)
        return elem

    def add_tf(self, sid, nid0, c, b0, b1, b2, nids, components, a, comment=''):
        """Creates a TF card"""
        tf = TF(sid, nid0, c, b0, b1, b2, nids, components, a, comment=comment)
        self._add_tf_object(tf)
        return tf

    def add_deqatn(self, equation_id, eqs, comment=''):
        """
        Creates a DEQATN card

        Parameters
        ----------
        equation_id : int
            the id of the equation
        eqs : List[str]
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
        >>> deqatn = model.add_deqatn(41, eqs, comment='')

        """
        deqatn = DEQATN(equation_id, eqs, comment=comment)
        self._add_deqatn_object(deqatn)
        return deqatn

    def add_desvar(self, desvar_id, label, xinit, xlb=-1e20, xub=1e20,
                   delx=None, ddval=None, comment=''):
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
        desvar = DESVAR(desvar_id, label, xinit, xlb, xub, delx=delx,
                        ddval=ddval, comment=comment)
        self._add_desvar_object(desvar)
        return desvar

    def add_dresp1(self, dresp_id, label, response_type, property_type, region,
                   atta, attb, atti, validate=True, comment=''):
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
        lable : str
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
        atti : List[int / float / str]
            the response values to pull from
            List[int]:
                list of grid ids
                list of property ids
            List[str]
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
        dresp = DRESP1(dresp_id, label, response_type, property_type, region,
                       atta, attb, atti, validate=validate, comment=comment)
        self._add_dresp_object(dresp)
        return dresp

    def add_dresp2(self, dresp_id, label, dequation, region, params,
                   method='MIN', c1=1., c2=0.005, c3=10.,
                   validate=True, comment=''):
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
            values : List[int]
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
        dresp = DRESP2(dresp_id, label, dequation, region, params,
                       method=method, c1=c1, c2=c2, c3=c3, comment=comment,
                       validate=validate)
        self._add_dresp_object(dresp)
        return dresp

    def add_dresp3(self, dresp_id, label, group, Type, region, params,
                   validate=True, comment=''):
        """Creates a DRESP3 card"""
        dresp = DRESP3(dresp_id, label, group, Type, region, params,
                       validate=validate, comment=comment)
        self._add_dresp_object(dresp)
        return dresp

    def add_dvcrel1(self, oid, Type, eid, cp_name, dvids, coeffs,
                    cp_min=None, cp_max=1e20, c0=0., validate=True, comment=''):
        """
        Creates a DVCREL1 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        EID : int
            element id
        cp_name : str/int
            optimization parameter as an element connectivity name (e.g., X1)
        dvids : List[int]
            DESVAR ids
        coeffs : List[float]
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
        dvcrel = DVCREL1(oid, Type, eid, cp_name, dvids, coeffs,
                         cp_min=cp_min, cp_max=cp_max, c0=c0,
                         validate=validate, comment=comment)
        self._add_dvcrel_object(dvcrel)
        return dvcrel

    def add_dvcrel2(self, oid, Type, eid, cp_name, deqation, dvids, labels,
                    cp_min=None, cp_max=1e20, validate=True, comment=''):
        """Creates a DVCREL2 card"""
        dvcrel = DVCREL2(oid, Type, eid, cp_name, deqation, dvids, labels,
                         cp_min=cp_min, cp_max=cp_max,
                         validate=validate, comment=comment)
        self._add_dvcrel_object(dvcrel)
        return dvcrel

    def add_dvprel1(self, oid, prop_type, pid, pname_fid, dvids, coeffs,
                    p_min=None, p_max=1e20, c0=0.0, validate=True, comment=''):
        """
        Creates a DVPREL1 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        dvids : List[int]
            DESVAR ids
        coeffs : List[float]
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
        dvprel = DVPREL1(oid, prop_type, pid, pname_fid, dvids, coeffs,
                         p_min=p_min, p_max=p_max,
                         c0=c0, validate=validate, comment=comment)
        self._add_dvprel_object(dvprel)
        return dvprel

    def add_dvprel2(self, oid, prop_type, pid, pname_fid, deqation,
                    dvids=None, labels=None, p_min=None, p_max=1e20,
                    validate=True, comment=''):
        """
        Creates a DVPREL2 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        deqation : int
            DEQATN id
        dvids : List[int]; default=None
            DESVAR ids
        labels : List[str]; default=None
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
        dvprel = DVPREL2(oid, prop_type, pid, pname_fid, deqation, dvids, labels,
                         p_min=p_min, p_max=p_max, validate=validate, comment=comment)
        self._add_dvprel_object(dvprel)
        return dvprel

    def add_dvmrel1(self, oid, mat_type, mid, mp_name, dvids, coeffs,
                    mp_min=None, mp_max=1e20, c0=0., validate=True, comment=''):
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
        dvids : List[int]
            DESVAR ids
        coeffs : List[float]
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
        dvmrel = DVMREL1(oid, mat_type, mid, mp_name, dvids, coeffs,
                         mp_min, mp_max, c0=c0, validate=validate, comment=comment)
        self._add_dvmrel_object(dvmrel)
        return dvmrel

    def add_dvmrel2(self, oid, mat_type, mid, mp_name, deqation, dvids, labels,
                    mp_min=None, mp_max=1e20, validate=True, comment=''):
        """
        Creates a DVMREL2 card

        Parameters
        ----------
        oid : int
            optimization id
        mat_type : str
            material card name (e.g., MAT1)
        mid : int
            material id
        mp_name : str
            optimization parameter as a pname (material name; E)
        deqation : int
            DEQATN id
        dvids : List[int]; default=None
            DESVAR ids
        labels : List[str]; default=None
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
        dvmrel = DVMREL2(oid, mat_type, mid, mp_name, deqation, dvids, labels,
                         mp_min=mp_min, mp_max=mp_max,
                         validate=validate, comment=comment)
        self._add_dvmrel_object(dvmrel)
        return dvmrel

    def add_dvgrid(self, dvid, nid, dxyz, cid=0, coeff=1.0, comment=''):
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
        dvgrid = DVGRID(dvid, nid, dxyz, cid=cid, coeff=coeff, comment=comment)
        self._add_dvgrid_object(dvgrid)
        return dvgrid

    def add_ddval(self, oid, ddvals, comment=''):
        """Creates a DDVAL card"""
        ddval = DDVAL(oid, ddvals, comment=comment)
        self._add_ddval_object(ddval)
        return ddval

    def add_dlink(self, oid, dependent_desvar,
                  independent_desvars, coeffs, c0=0., cmult=1., comment=''):
        """
        Creates a DLINK card, which creates a variable that is a lienar
        ccombination of other design variables

        Parameters
        ----------
        oid : int
            optimization id
        dependent_desvar : int
            the DESVAR to link
        independent_desvars : List[int]
            the DESVARs to combine
        coeffs : List[int]
            the linear combination coefficients
        c0 : float; default=0.0
            an offset
        cmult : float; default=1.0
            an scale factor
        comment : str; default=''
            a comment for the card

        """
        dlink = DLINK(oid, dependent_desvar,
                      independent_desvars, coeffs, c0=c0, cmult=cmult, comment=comment)
        self._add_dlink_object(dlink)
        return dlink

    def add_dconstr(self, oid, dresp_id, lid=-1.e20, uid=1.e20, lowfq=0.,
                    highfq=1.e20, comment=''):
        """
        Creates a DCONSTR card

        Parameters
        ----------
        oid : int
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
        dconstr = DCONSTR(oid, dresp_id, lid=lid, uid=uid, lowfq=lowfq,
                          highfq=highfq, comment=comment)
        self._add_dconstr_object(dconstr)
        return dconstr

    def add_dconadd(self, oid, dconstrs, comment=''):
        """Creates a DCONADD card"""
        dconadd = DCONADD(oid, dconstrs, comment=comment)
        self._add_dconstr_object(dconadd)
        return dconadd

    def add_doptprm(self, params, comment=''):
        """Creates a DOPTPRM card"""
        doptprm = DOPTPRM(params, comment=comment)
        self._add_doptprm_object(doptprm)
        return doptprm

    def add_dscreen(self, rtype, trs=-0.5, nstr=20, comment=''):
        """
        Creates a DSCREEN object

        Parameters
        ----------
        rtype : str
            Response type for which the screening criteria apply
        trs : float
            Truncation threshold
        nstr : int
            Maximum number of constraints to be retained per region per
            load case
        comment : str; default=''
            a comment for the card

        """
        dscreen = DSCREEN(rtype, trs=trs, nstr=nstr, comment=comment)
        self._add_dscreen_object(dscreen)
        return dscreen

    def add_monpnt1(self, name, label, axes, aecomp_name, xyz, cp=0, cd=None,
                    comment=''):
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
        xyz : List[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, CORDx; default=0
           int : coordinate system
        cd : int; default=None -> cp
            the coordinate system for load outputs
        comment : str; default=''
            a comment for the card

        Notes
        -----
        CD - MSC specific field

        """
        monitor_point = MONPNT1(name, label, axes, aecomp_name, xyz, cp=cp, cd=cd,
                                comment=comment)
        self._add_monpnt_object(monitor_point)
        return monitor_point

    def add_monpnt2(self, name, label, table, Type, nddl_item, eid,
                    comment=''):
        """Creates a MONPNT2 card"""
        monitor_point = MONPNT2(name, label, table, Type, nddl_item, eid,
                                comment=comment)
        self._add_monpnt_object(monitor_point)
        return monitor_point

    def add_monpnt3(self, name, label, axes, grid_set, elem_set, xyz,
                    cp=0, cd=None, xflag=None, comment=''):
        """Creates a MONPNT3 card"""
        monitor_point = MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                                cp=cp, cd=cd, xflag=xflag, comment=comment)
        self._add_monpnt_object(monitor_point)
        return monitor_point

    def add_bsurfs(self, id, eids, g1s, g2s, g3s, comment=''):
        """Creates a BSURFS card"""
        bsurfs = BSURFS(id, eids, g1s, g2s, g3s, comment=comment)
        self._add_bsurfs_object(bsurfs)
        return bsurfs

    def add_bsurf(self, sid, eids, comment=''):
        """Creates a BSURF card"""
        bsurf = BSURF(sid, eids, comment=comment)
        self._add_bsurf_object(bsurf)
        return bsurf

    def add_bctset(self, csid, sids, tids, frictions, min_distances,
                   max_distances, comment='', sol=101):
        """Creates a BCTSET card"""
        bctset = BCTSET(csid, sids, tids, frictions, min_distances,
                        max_distances, comment=comment, sol=sol)
        self._add_bctset_object(bctset)
        return bctset

    def add_bctadd(self, csid, contact_sets, comment=''):
        """Creates a BCTADD card"""
        bctadd = BCTADD(csid, contact_sets, comment=comment)
        self._add_bctadd_object(bctadd)
        return bctadd

    def add_bctpara(self, csid, params, comment=''):
        """Creates a BCTPARA card"""
        bctpara = BCTPARA(csid, params, comment=comment)
        self._add_bctpara_object(bctpara)
        return bctpara

    def add_bcrpara(self, crid, surf='TOP', offset=None, Type='FLEX',
                    grid_point=0, comment=''):
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
        Type : str; default='FLEX'
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
        bcrpara = BCRPARA(crid, surf=surf, offset=offset, Type=Type,
                          grid_point=grid_point, comment=comment)
        self._add_bcrpara_object(bcrpara)
        return bcrpara

    def add_tic(self, sid, nodes, components, u0=0., v0=0., comment=''):
        """
        Creates a TIC card

        Parameters
        ----------
        sid : int
            Case Control IC id
        nodes : int / List[int]
            the nodes to which apply the initial conditions
        components : int / List[int]
            the DOFs to which apply the initial conditions
        u0 : float / List[float]
            Initial displacement.
        v0 : float / List[float]
            Initial velocity.
        comment : str; default=''
            a comment for the card

        """
        tic = TIC(sid, nodes, components, u0=u0, v0=v0, comment=comment)
        self._add_tic_object(tic)
        return tic

    def add_tstep1(self, sid, tend, ninc, nout, comment=''):
        """Creates a TSTEP1 card"""
        tstep1 = TSTEP1(sid, tend, ninc, nout, comment=comment)
        self._add_tstep_object(tstep1)
        return tstep1

    def add_tstep(self, sid, N, DT, NO, comment=''):
        """Creates a TSTEP card"""
        tstep = TSTEP(sid, N, DT, NO, comment=comment)
        self._add_tstep_object(tstep)
        return tstep

    def add_tstepnl(self, sid, ndt, dt, no, method='ADAPT', kstep=None,
                    max_iter=10, conv='PW', eps_u=1.e-2, eps_p=1.e-3,
                    eps_w=1.e-6, max_div=2, max_qn=10, max_ls=2,
                    fstress=0.2, max_bisect=5, adjust=5, mstep=None,
                    rb=0.6, max_r=32., utol=0.1, rtol_b=20.,
                    min_iter=None, comment=''):
        """Creates a TSTEPNL card"""
        tstepnl = TSTEPNL(sid, ndt, dt, no, method=method, kstep=kstep,
                          max_iter=max_iter, conv=conv, eps_u=eps_u, eps_p=eps_p,
                          eps_w=eps_w, max_div=max_div, max_qn=max_qn, max_ls=max_ls,
                          fstress=fstress, max_bisect=max_bisect, adjust=adjust,
                          mstep=mstep, rb=rb, max_r=max_r, utol=utol, rtol_b=rtol_b,
                          min_iter=min_iter, comment=comment)
        self._add_tstepnl_object(tstepnl)
        return tstepnl

    def add_nlparm(self, nlparm_id, ninc=None, dt=0.0, kmethod='AUTO', kstep=5,
                   max_iter=25, conv='PW', int_out='NO', eps_u=0.01,
                   eps_p=0.01, eps_w=0.01, max_div=3, max_qn=None, max_ls=4,
                   fstress=0.2, ls_tol=0.5, max_bisect=5, max_r=20., rtol_b=20., comment=''):
        """Creates an NLPARM card"""
        nlparm = NLPARM(nlparm_id, ninc=ninc, dt=dt, kmethod=kmethod, kstep=kstep,
                        max_iter=max_iter, conv=conv, int_out=int_out, eps_u=eps_u,
                        eps_p=eps_p, eps_w=eps_w, max_div=max_div, max_qn=max_qn,
                        max_ls=max_ls, fstress=fstress, ls_tol=ls_tol,
                        max_bisect=max_bisect, max_r=max_r, rtol_b=rtol_b, comment=comment)
        self._add_nlparm_object(nlparm)
        return nlparm

    def add_nlpci(self, nlpci_id, Type='CRIS', minalr=0.25, maxalr=4.,
                  scale=0., desiter=12, mxinc=20, comment=''):
        """Creates an NLPCI card"""
        nlpci = NLPCI(nlpci_id, Type=Type, minalr=minalr, maxalr=maxalr,
                      scale=scale, desiter=desiter, mxinc=mxinc, comment=comment)
        self._add_nlpci_object(nlpci)
        return nlpci

    def add_delay(self, sid, nodes, components, delays, comment=''):
        """
        Creates a DELAY card

        Parameters
        ----------
        sid : int
            DELAY id that is referenced by a TLOADx, RLOADx or ACSRCE card
        nodes : List[int]
            list of nodes that see the delay
            len(nodes) = 1 or 2
        components : List[int]
            the components corresponding to the nodes that see the delay
            len(nodes) = len(components)
        delays : List[float]
            Time delay (tau) for designated point Pi and component Ci
            len(nodes) = len(delays)
        comment : str; default=''
            a comment for the card

        """
        delay = DELAY(sid, nodes, components, delays, comment=comment)
        self._add_delay_object(delay)
        return delay

    def add_dphase(self, sid, nodes, components, phase_leads, comment=''):
        """
        Creates a DPHASE card

        Parameters
        ----------
        sid : int
            DPHASE id that is referenced by a RLOADx or ACSRCE card
        nodes : List[int]
            list of nodes that see the delay
            len(nodes) = 1 or 2
        components : List[int]
            the components corresponding to the nodes that see the delay
            len(nodes) = len(components)
        phase_leads : List[float]
            Phase lead θ in degrees.
            len(nodes) = len(delays)
        comment : str; default=''
            a comment for the card

        """
        dphase = DPHASE(sid, nodes, components, phase_leads, comment=comment)
        self._add_dphase_object(dphase)
        return dphase

    def add_rotorg(self, sid, nids, comment=''):
        """Creates a ROTORG card"""
        rotor = ROTORG(sid, nids, comment=comment)
        self._add_rotor_object(rotor)
        return rotor

    def add_rotord(self, sid, rstart, rstep, numstep,
                   rids, rsets, rspeeds, rcords, w3s, w4s, rforces, brgsets,
                   refsys='ROT', cmout=0.0, runit='RPM', funit='RPM',
                   zstein='NO', orbeps=1.e-6, roprt=0, sync=1, etype=1,
                   eorder=1.0, threshold=0.02, maxiter=10, comment=''):
        """Creates a ROTORD card"""
        rotor = ROTORD(sid, rstart, rstep, numstep, rids, rsets, rspeeds,
                       rcords, w3s, w4s, rforces,
                       brgsets, refsys=refsys, cmout=cmout,
                       runit=runit, funit=funit, zstein=zstein, orbeps=orbeps,
                       roprt=roprt, sync=sync, etype=etype, eorder=eorder,
                       threshold=threshold, maxiter=maxiter, comment=comment)
        self._add_rotor_object(rotor)
        return rotor

    def add_cmfree(self, eid, s, s2, y, n):
        fields = ['CMFREE', eid, s, s2, y, n]
        self.reject_card_lines('CMFREE', print_card_8(fields).split('\n'))

    def add_cfluid2(self, eid, ringfls, rho, b, harmonic):
        fields = ['CFLUID2', eid] + ringfls + [rho, b, harmonic]
        self.reject_card_lines('CFLUID2', print_card_8(fields).split('\n'))

    def add_cfluid3(self, eid, ringfls, rho, b, harmonic):
        fields = ['CFLUID3', eid] + ringfls + [rho, b, harmonic]
        self.reject_card_lines('CFLUID3', print_card_8(fields).split('\n'))

    def add_cfluid4(self, eid, ringfls, rho, b, harmonic):
        fields = ['CFLUID4', eid] + ringfls + [rho, b, harmonic]
        self.reject_card_lines('CFLUID4', print_card_8(fields).split('\n'))

    def add_rgyro(self, sid, asynci, refrot, unit, speed_low, speed_high, speed, comment=''):
        """Creates an RGYRO card"""
        fields = ['RGYRO', sid, asynci, refrot, unit, speed_low, speed_high, speed]
        self.reject_card_lines('RGYRO', print_card_8(fields).split('\n'))

    def add_rspint(self, rid, grida, gridb, gr, unit, table_id, comment=''):
        """Creates an RSPINT card"""
        fields = ['RSPINT', rid, grida, gridb, gr, unit, table_id]
        self.reject_card_lines('RSPINT', print_card_8(fields).split('\n'))

    def add_temp(self, sid, temperatures, comment=''):
        """
        Creates a TEMP card

        Parameters
        ----------
        sid : int
            Load set identification number
        temperatures : dict[nid] : temperature
            nid : int
                node id
            temperature : float
                the nodal temperature
        comment : str; default=''
            a comment for the card

        """
        temp = TEMP(sid, temperatures, comment=comment)
        self._add_thermal_load_object(temp)
        return temp

    #def add_tempp1(self):
        #temp = TEMPP1()
        #self._add_thermal_load_object(temp)
        #return temp

    def add_tempd(self, sid, temperature, comment=''):
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
        tempd = TEMPD(sid, temperature, comment=comment)
        self._add_tempd_object(tempd)
        return tempd

    def add_qhbdy(self, sid, flag, q0, grids, af=None, comment=''):
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
        af : float; default=None
            Area factor depends on type
        grids : List[int]
            Grid point identification of connected grid points
        comment : str; default=''
            a comment for the card

        """
        load = QHBDY(sid, flag, q0, grids, af=af, comment=comment)
        self._add_thermal_load_object(load)
        return load

    def add_qbdy1(self, sid, qflux, eids, comment=''):
        """Creates a QBDY1 card"""
        load = QBDY1(sid, qflux, eids, comment=comment)
        self._add_thermal_load_object(load)
        return load

    def add_qbdy2(self, sid, eid, qfluxs, comment=''):
        """Creates a QBDY1 card"""
        load = QBDY2(sid, eid, qfluxs, comment=comment)
        self._add_thermal_load_object(load)
        return load

    def add_qbdy3(self, sid, Q0, cntrlnd, eids, comment=''):
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
        eids : List[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        load = QBDY3(sid, Q0, cntrlnd, eids, comment=comment)
        self._add_thermal_load_object(load)
        return load

    def add_qvol(self, sid, qvol, control_point, elements, comment=''):
        """Creates a QVOL card"""
        load = QVOL(sid, qvol, control_point, elements, comment=comment)
        self._add_load_object(load)
        return load

    def add_qvect(self, sid, q0, eids, t_source=None,
                  ce=0, vector_tableds=None, control_id=0, comment=''):
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
        vector_tableds : List[int/float, int/float, int/float]
            vector : float; default=0.0
                directional cosines in coordinate system CE) of
                the thermal vector flux
            tabled : int
                TABLEDi entry identification numbers defining the
                components as a function of time
        control_id : int; default=0
            Control point
        eids : List[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """
        load = QVECT(sid, q0, eids, t_source=t_source, ce=ce,
                     vector_tableds=vector_tableds, control_id=control_id,
                     comment=comment)
        self._add_dload_entry(load)
        return load

    def add_chbdyg(self, eid, surface_type, nodes,
                   iview_front=0, ivew_back=0,
                   rad_mid_front=0, rad_mid_back=0, comment=''):
        """Creates a CHBDYG card"""
        elem = CHBDYG(eid, surface_type, nodes,
                      iview_front=iview_front, ivew_back=ivew_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      comment=comment)
        self._add_thermal_element_object(elem)
        return elem

    def add_chbdyp(self, eid, pid, surface_type, g1, g2,
                   g0=0, gmid=None, ce=0,
                   iview_front=0, ivew_back=0,
                   rad_mid_front=0, rad_mid_back=0,
                   e1=None, e2=None, e3=None,
                   comment=''):
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
        ivew_back : int; default=0
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
        elem = CHBDYP(eid, pid, surface_type, g1, g2,
                      g0=g0, gmid=gmid, ce=ce,
                      iview_front=iview_front, ivew_back=ivew_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      e1=e1, e2=e2, e3=e3,
                      comment=comment)
        self._add_thermal_element_object(elem)
        return elem

    def add_chbdye(self, eid, eid2, side,
                   iview_front=0, ivew_back=0,
                   rad_mid_front=0, rad_mid_back=0,
                   comment=''):
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
        ivew_back: int; default=0
            a VIEW entry identification number for the back face
        rad_mid_front: int; default=0
            RADM identification number for front face of surface element
        rad_mid_back: int; default=0
            RADM identification number for back face of surface element
        comment : str; default=''
            a comment for the card

        """
        elem = CHBDYE(eid, eid2, side,
                      iview_front=iview_front, ivew_back=ivew_back,
                      rad_mid_front=rad_mid_front, rad_mid_back=rad_mid_back,
                      comment=comment)
        self._add_thermal_element_object(elem)
        return elem

    def add_phbdy(self, pid, af=None, d1=None, d2=None, comment=''):
        """
        Creates a PHBDY card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        af : int
            Area factor of the surface used only for CHBDYP element
            Must be {POINT, LINE, TUBE, ELCYL}
            TUBE : constant thickness of hollow tube
        d1, d2 : float; default=None
            Diameters associated with the surface
            Used with CHBDYP [ELCYL, TUBE, FTUBE] surface elements
        comment : str; default=''
            a comment for the card

        """
        prop = PHBDY(pid, af=af, d1=d1, d2=d2, comment=comment)
        self._add_phbdy_object(prop)
        return prop

    def add_conv(self, eid, pconid, ta, film_node=0, cntrlnd=0, comment=''):
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
        ta : List[int]
            Ambient points used for convection 0's are allowed for TA2
            and higher
        film_node : int; default=0
            Point for film convection fluid property temperature
        cntrlnd : int; default=0
            Control point for free convection boundary condition
        comment : str; default=''
            a comment for the card

        """
        boundary_condition = CONV(eid, pconid, ta,
                                  film_node=film_node, cntrlnd=cntrlnd,
                                  comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def add_convm(self, eid, pconvm, ta1, film_node=0, cntmdot=0,
                  ta2=None, mdot=1.0,
                  comment=''):
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
        boundary_condition = CONVM(eid, pconvm, ta1,
                                   film_node=film_node, cntmdot=cntmdot,
                                   ta2=ta2, mdot=mdot,
                                   comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def add_radm(self, radmid, absorb, emissivity, comment=''):
        """Creates a RADM card"""
        boundary_condition = RADM(radmid, absorb, emissivity, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.radmid)
        return boundary_condition

    def add_radbc(self, nodamb, famb, cntrlnd, eids, comment=''):
        """Creates a RADBC card"""
        boundary_condition = RADBC(nodamb, famb, cntrlnd, eids, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.nodamb)
        return boundary_condition

    def add_pconv(self, pconid, mid, form=0, expf=0.0, ftype=0, tid=None,
                  chlen=None, gidin=None, ce=0,
                  e1=None, e2=None, e3=None,
                  comment=''):
        """
        Creates a PCONV card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid : int
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
        e1 / e2 / e3 : List[float]; default=None
            Components of the orientation vector in coordinate system CE.
            The origin of the orientation vector is grid point G1
        comment : str; default=''
            a comment for the card

        """
        prop = PCONV(pconid, mid,
                     form=form, expf=expf, ftype=ftype,
                     tid=tid, chlen=chlen, gidin=gidin,
                     ce=ce, e1=e1, e2=e2, e3=e3, comment=comment)
        self._add_convection_property_object(prop)
        return prop

    def add_pconvm(self, pconid, mid, coef, form=0, flag=0,
                   expr=0.0, exppi=0.0, exppo=0.0, comment=''):
        """
        Creates a PCONVM card

        Parameters
        ----------
        pconid : int
            Convection property ID
        mid: int
            Material ID
        coef: float
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
        prop = PCONVM(pconid, mid, coef, form=form, flag=flag,
                      expr=expr, exppi=exppi, exppo=exppo, comment=comment)
        self._add_convection_property_object(prop)
        return prop

    def add_dti(self, name, fields, comment=''):
        """Creates a DTI card"""
        if name.upper() == 'UNITS':
            dti = DTI(name, fields, comment=comment)
            self._add_dti_object(dti)
        else:
            if comment:
                self.reject_lines.append([_format_comment(comment)])
                msg = "DTI only supports name='UNITS'; name=%r fields=%s" % (name, str(fields))
            raise NotImplementedError(msg)
            #self.reject_cards.append(card_obj)
            #self._write_reject_message(card_name, card_obj, comment=comment)
        return dti

    def add_dmig_uaccel(self, tin, ncol, load_sequences, comment=''):
        """Creates a DMIG,UACCEL card"""
        dmig = DMIG_UACCEL(tin, ncol, load_sequences, comment=comment)
        self._add_dmig_object(dmig)
        return dmig

    def add_dmig(self, name, ifo, tin, tout, polar, ncols, GCj, GCi,
                 Real, Complex=None, comment=''):
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
        GCj  : List[(node, dof)]
            the [jnode, jDOFs]
        GCi  : List[(node, dof)]
            the inode, iDOFs
        Real : List[float]
            The real values
        Complex : List[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        dmig = DMIG(name, ifo, tin, tout, polar, ncols, GCj, GCi,
                    Real, Complex, comment=comment)
        self._add_dmig_object(dmig)
        return dmig

    def add_dmi(self, name, form, tin, tout, nrows, ncols, GCj, GCi,
                Real, Complex=None, comment=''):
        """Creates a DMI card"""
        dmi = DMI(name, form, tin, tout, nrows, ncols, GCj, GCi, Real,
                  Complex, comment=comment)
        self._add_dmi_object(dmi)
        return dmi

    def add_dmij(self, name, form, tin, tout, nrows, ncols, GCj, GCi,
                 Real, Complex=None, comment=''):
        """Creates a DMIJ card"""
        dmij = DMIJ(name, form, tin, tout, nrows, ncols, GCj, GCi,
                    Real, Complex, comment=comment)
        self._add_dmij_object(dmij)
        return dmij

    def add_dmiji(self, name, ifo, tin, tout, nrows, ncols, GCj, GCi,
                  Real, Complex=None, comment=''):
        """
        | DMIJI | NAME | 0 | IFO | TIN | TOUT POLAR | | NCOL |
        """
        dmiji = DMIJI(name, ifo, tin, tout, nrows, ncols, GCj, GCi,
                      Real, Complex, comment=comment)
        self._add_dmiji_object(dmiji)
        return dmiji

    def add_dmik(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
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
        GCj  : List[(node, dof)]
            the jnode, jDOFs
        GCi  : List[(node, dof)]
            the inode, iDOFs
        Real : List[float]
            The real values
        Complex : List[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        dmik = DMIK(name, ifo, tin, tout, polar, ncols,
                    GCj, GCi, Real, Complex, comment=comment)
        self._add_dmik_object(dmik)
        return dmik

    def add_cgen(self, Type, field_eid, pid, field_id, th_geom_opt,
                 eidl, eidh, t_abcd=None, direction='L', comment=''):
        """Creates a CGEN card"""
        elem = CGEN(Type, field_eid, pid, field_id, th_geom_opt,
                    eidl, eidh, t_abcd=t_abcd, direction=direction, comment=comment)
        self._add_element_object(elem)
        return elem
