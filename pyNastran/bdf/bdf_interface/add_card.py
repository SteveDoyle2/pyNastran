"""defines various methods to add cards"""
# pylint: disable=R0913, R0914, C0103
from __future__ import print_function

import numpy as np

from pyNastran.bdf.bdf_interface.add_methods import AddMethods
from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP, PRAC2D, PRAC3D, PCONEAX
from pyNastran.bdf.cards.properties.solid import PLSOLID, PSOLID, PIHEX, PCOMPS

from pyNastran.bdf.cards.elements.springs import (CELAS1, CELAS2, CELAS3, CELAS4,)
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.elements.solid import (CTETRA, CPYRAM, CPENTA, CHEXA, CIHEX1)
from pyNastran.bdf.cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE

from pyNastran.bdf.cards.elements.shell import (CQUAD, CQUAD4, CQUAD8, CQUADR, CQUADX,
                                                CSHEAR, CTRIA3, CTRIA6, CTRIAX,
                                                CTRIAX6, CTRIAR,
                                                CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8)
from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE, PPLANE
from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D, PBUSHT
from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                                 CDAMP5)
from pyNastran.bdf.cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE
from pyNastran.bdf.cards.elements.bars import CBAR, CBEAM3, CBEND
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBRSECT  # PBEND
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP, PBMSECT
# CMASS5
from pyNastran.bdf.cards.elements.mass import CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4
from pyNastran.bdf.cards.properties.mass import PMASS#, NSM
from pyNastran.bdf.cards.constraints import (SPC, SPCADD, SPCAX, SPC1,
                                             MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                             GMSPC)
from pyNastran.bdf.cards.coordinate_systems import (CORD1R, CORD1C, CORD1S,
                                                    CORD2R, CORD2C, CORD2S, #CORD3G,
                                                    GMCORD)
from pyNastran.bdf.cards.deqatn import DEQATN
from pyNastran.bdf.cards.dynamic import (
    DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ4,
    TSTEP, TSTEPNL, NLPARM, NLPCI, TF, ROTORG, ROTORD)
from pyNastran.bdf.cards.loads.loads import (
    LSEQ, SLOAD, DAREA, RANDPS, RFORCE, RFORCE1, SPCD, LOADCYN)
from pyNastran.bdf.cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
from pyNastran.bdf.cards.loads.static_loads import (LOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                    FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                    PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                                    GMLOAD)

from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11,
                                           MATHE, MATHP, CREEP, EQUIV)
from pyNastran.bdf.cards.material_deps import MATT1, MATT2, MATT4, MATT5, MATS1

from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL
from pyNastran.bdf.cards.nodes import GRID, GRDSET, SPOINTs, EPOINTs, POINT
from pyNastran.bdf.cards.aero import (
    AECOMP, AEFACT, AELINK, AELIST, AEPARM, AESTAT,
    AESURF, AESURFS, AERO, AEROS, CSSCHD,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    MONPNT1, MONPNT2, MONPNT3,
    FLFACT, FLUTTER, GUST, MKAERO1,
    MKAERO2, SPLINE1, SPLINE2, SPLINE3, SPLINE4,
    SPLINE5, TRIM, DIVERG)
from pyNastran.bdf.cards.optimization import (
    DCONADD, DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
    DRESP1, DRESP2, DRESP3,
    DVCREL1, DVCREL2,
    DVMREL1, DVMREL2,
    DVPREL1, DVPREL2,
    DVGRID)
from pyNastran.bdf.cards.bdf_sets import (
    ASET, BSET, CSET, QSET, USET,
    ASET1, BSET1, CSET1, QSET1, USET1,
    SET1, SET3, #RADSET,
    SEBSET, SECSET, SEQSET, # SEUSET
    SEBSET1, SECSET1, SEQSET1, # SEUSET1
    SESET, #SEQSEP
)
from pyNastran.bdf.cards.params import PARAM
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL
from pyNastran.bdf.cards.thermal.loads import QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, QVOL
from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                 PHBDY, CONV, CONVM, RADM, RADBC)
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                            TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                            TABLES1, TABDMP1, TABLEST, TABRND1, TABRNDG, #TIC,
                                            DTABLE)
from pyNastran.bdf.cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA


class AddCards(AddMethods):
    """defines the add_cardname functions that use the object inits"""
    def __init__(self):
        AddMethods.__init__(self)

    def add_grid(self, nid, cp=0, xyz=None, cd=0, ps='', seid=0, comment=''):
        """
        Creates the GRID card in a functional way

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
        """
        grid = GRID(nid, cp=cp, xyz=xyz, cd=cd, ps=ps, seid=seid, comment=comment)
        self.add_node_object(grid)

    def add_grdset(self, cp, cd, ps, seid, comment):
        grdset = GRDSET(cp, cd, ps, seid, comment=comment)
        self.gridSet = grdset

    def add_spoint(self, ids, comment=''):
        """
        Creates the SPOINTs card that contains many SPOINTs

        Parameters
        ----------
        ids : List[int]
            SPOINT ids
        comment : str
            a comment for the card
        """
        spoint = SPOINTs(ids, comment=comment)
        self.add_spoint_object(spoint)

    def add_epoint(self, ids, comment=''):
        """
        Creates the EPOINTs card that contains many EPOINTs

        Parameters
        ----------
        ids : List[int]
            EPOINT ids
        comment : str
            a comment for the card
        """
        epoint = EPOINTs(ids, comment=comment)
        self.add_epoint_object(epoint)

    def add_point(self, nid, cp, xyz, comment=''):
        """
        Creates the POINT card
        """
        point = POINT(nid, cp, xyz, comment=comment)
        self.add_point_object(point)

    def add_cord2r(self, cid, rid=0, origin=None, zaxis=None, xzplane=None,
                   comment=''):
        coord = CORD2R(cid, rid=rid, origin=origin, zaxis=zaxis, xzplane=xzplane,
                       comment=comment)
        self.add_coord_object(coord)

    def add_cord2c(self, cid, rid=0, origin=None, zaxis=None, xzplane=None,
                   comment=''):
        coord = CORD2C(cid, rid=rid, origin=origin, zaxis=zaxis, xzplane=xzplane,
                       comment=comment)
        self.add_coord_object(coord)

    def add_cord2s(self, cid, rid=0, origin=None, zaxis=None, xzplane=None,
                   comment=''):
        coord = CORD2S(cid, rid=rid, origin=origin, zaxis=zaxis, xzplane=xzplane,
                       comment=comment)
        self.add_coord_object(coord)

    def add_cord1r(self, cid, g1, g2, g3, comment=''):
        coord = CORD1R(cid, g1, g2, g3, comment=comment)
        self.add_coord_object(coord)

    def add_cord1c(self, cid, g1, g2, g3, comment=''):
        coord = CORD1C(cid, g1, g2, g3, comment=comment)
        self.add_coord_object(coord)

    def add_cord1s(self, cid, g1, g2, g3, comment=''):
        coord = CORD1S(cid, g1, g2, g3, comment=comment)
        self.add_coord_object(coord)

    def add_param(self, key, values, comment=''):
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/List
            varies depending on the type of PARAM
        comment : str; default=''
            optional string
        """
        param = PARAM(key, values, comment=comment)
        self.add_param_object(param)

    def add_plotel(self, eid, nodes, comment=''):
        elem = PLOTEL(eid, nodes, comment=comment)
        self.add_plotel_object(elem)

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
        ::

          [M] = [M11 M21 M31 M41 M51 M61]
                [    M22 M32 M42 M52 M62]
                [        M33 M43 M53 M63]
                [            M44 M54 M64]
                [    Sym         M55 M65]
                [                    M66]
        """
        mass = CONM1(eid, nid, mass_matrix, cid=cid, comment=comment)
        self.add_mass_object(mass)

    def add_conm2(self, eid, nid, cid, mass, X=None, I=None, comment=''):
        """
        Creates a CONM2 card

        Parameters
        ----------
        eid : int
           element ID
        nid : int
           node ID
        mass : float
           the mass of the CONM2
        cid : int; default=0
           coordinate frame of the offset (-1=absolute coordinates)
        X : (3, ) List[float]; default=None -> [0., 0., 0.]
            xyz offset vector relative to nid
        I : (6, ) List[float]; default=None -> [0., 0., 0., 0., 0., 0.]
            mass moment of inertia matrix about the CG
            I11, I21, I22, I31, I32, I33 = I
        """
        mass = CONM2(eid, nid, mass, cid=cid, X=X, I=I, comment=comment)
        self.add_mass_object(mass)

    def add_pmass(self, pid, mass, comment=''):
        prop = PMASS(pid, mass, comment=comment)
        self.add_property_mass_object(prop)

    def add_cmass1(self, eid, pid, g1, c1, g2, c2, comment=''):
        mass = CMASS1(eid, pid, g1, c1, g2, c2, comment=comment)
        self.add_mass_object(mass)

    def add_cmass2(self, eid, mass, g1, c1, g2, c2, comment=''):
        mass = CMASS2(eid, mass, g1, c1, g2, c2, comment=comment)
        self.add_mass_object(mass)

    def add_cmass3(self, eid, pid, s1, s2, comment=''):
        mass = CMASS3(eid, pid, s1, s2, comment=comment)
        self.add_mass_object(mass)

    def add_cmass4(self, eid, mass, s1, s2=0, comment=''):
        mass = CMASS4(eid, mass, s1, s2=s2, comment=comment)
        self.add_mass_object(mass)

    def add_pelas(self, pid, k, ge, s, comment=''):
        prop = PELAS(pid, k, ge, s, comment=comment)
        self.add_property_object(prop)

    def add_celas1(self, eid, pid, nids, c1, c2, comment=''):
        elem = CELAS1(eid, pid, nids, c1, c2, comment=comment)
        self.add_element_object(elem)

    def add_celas2(self, eid, k, nids, c1=0, c2=0, ge=0., s=0., comment=''):
        elem = CELAS2(eid, k, nids, c1=c1, c2=c2, ge=ge, s=s, comment=comment)
        self.add_element_object(elem)

    def add_celas3(self, eid, pid, s1, s2, comment=''):
        elem = CELAS3(eid, pid, s1, s2, comment=comment)
        self.add_element_object(elem)

    def add_celas4(self, eid, k, s1, s2, comment=''):
        elem = CELAS4(eid, k, s1, s2, comment=comment)
        self.add_element_object(elem)

    def add_cdamp1(self, eid, pid, nids, c1, c2, comment=''):
        elem = CDAMP1(eid, pid, nids, c1, c2, comment=comment)
        self.add_element_object(elem)

    def add_cdamp2(self, eid, b, nids, c1, c2, comment=''):
        elem = CDAMP2(eid, b, nids, c1, c2, comment=comment)
        self.add_element_object(elem)

    def add_cdamp3(self, eid, pid, nids, comment=''):
        elem = CDAMP3(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_cdamp4(self, eid, b, nids, comment=''):
        elem = CDAMP4(eid, b, nids, comment=comment)
        self.add_element_object(elem)

    def add_cdamp5(self, eid, pid, nids, comment=''):
        elem = CDAMP5(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_pdamp(self, pid, b, comment=''):
        prop = PDAMP(pid, b, comment=comment)
        self.add_property_object(prop)

    def add_pdampt(self, pid, tbid, comment=''):
        prop = PDAMPT(pid, tbid, comment=comment)
        self.add_pdampt_object(prop)

    def add_pdamp5(self, pid, mid, b, comment=''):
        prop = PDAMP5(pid, mid, b, comment=comment)
        self.add_property_object(prop)

    def add_cvisc(self, eid, pid, nids, comment=''):
        elem = CVISC(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_pvisc(self, pid, ce, cr, comment=''):
        prop = PVISC(pid, ce, cr, comment=comment)
        self.add_property_object(prop)

    def add_cgap(self, eid, pid, ga, gb, x, g0, cid, comment=''):
        elem = CGAP(eid, pid, ga, gb, x, g0, cid, comment=comment)
        self.add_element_object(elem)

    def add_pgap(self, pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin,
                 comment=''):
        prop = PGAP(pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin,
                    comment=comment)
        self.add_property_object(prop)

    def add_cfast(self, eid, pid, Type, ida, idb, gs, ga, gb, xs, ys, zs,
                  comment=''):
        elem = CFAST(eid, pid, Type, ida, idb, gs, ga, gb, xs, ys, zs,
                     comment=comment)
        self.add_element_object(elem)

    def add_pfast(self, pid, d, mcid, mflag, kt1, kt2, kt3, kr1, kr2, kr3, mass,
                  ge, comment=''):
        prop = PFAST(pid, d, mcid, mflag, kt1, kt2, kt3, kr1, kr2, kr3, mass,
                     ge, comment=comment)
        self.add_property_object(prop)

    def add_cbush(self, eid, pid, ga, gb, x, g0, cid, s, ocid, si, comment=''):
        elem = CBUSH(eid, pid, ga, gb, x, g0, cid, s, ocid, si, comment=comment)
        self.add_element_object(elem)

    def add_pbush(self, pid, k, b, ge, rcv, comment=''):
        prop = PBUSH(pid, k, b, ge, rcv, comment=comment)
        self.add_property_object(prop)

    def add_cbush1d(self, eid, pid, ga, gb, cid, comment=''):
        elem = CBUSH1D(eid, pid, ga, gb, cid, comment=comment)
        self.add_element_object(elem)

    def add_cbush2d(self, eid, pid, ga, gb, cid, plane, sptid, comment=''):
        elem = CBUSH2D(eid, pid, ga, gb, cid, plane, sptid, comment=comment)
        self.add_element_object(elem)

    def add_pbush1d(self, pid, k, c, m, sa, se, optional_vars, comment=''):
        prop = PBUSH1D(pid, k, c, m, sa, se, optional_vars, comment=comment)
        self.add_property_object(prop)

    def add_pbusht(self, pid, k_tables, b_tables, ge_tables, kn_tables,
                   comment=''):
        prop = PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables,
                      comment=comment)
        self.add_property_object(prop)

    def add_pelast(self, pid, tkid, tgeid, tknid, comment=''):
        prop = PELAST(pid, tkid, tgeid, tknid, comment=comment)
        self.add_property_object(prop)

    def add_conrod(self, eid, mid, nids, A, j=0.0, c=0.0, nsm=0.0, comment=''):
        elem = CONROD(eid, mid, nids, A, j=j, c=c, nsm=nsm, comment=comment)
        self.add_element_object(elem)

    def add_crod(self, eid, pid, nids, comment=''):
        elem = CROD(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_prod(self, pid, mid, A, j=0., c=0., nsm=0., comment=''):
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
        self.add_property_object(prop)

    def add_ctube(self, eid, pid, nids, comment=''):
        elem = CTUBE(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_ptube(self, pid, mid, OD1, t, nsm, OD2, comment=''):
        prop = PTUBE(pid, mid, OD1, t, nsm, OD2, comment=comment)
        self.add_property_object(prop)

    def add_cbar(self, eid, pid, ga, gb, x, g0, offt='GGG', pa=0, pb=0,
                 wa=None, wb=None, comment=''):
        elem = CBAR(eid, pid, ga, gb, x, g0, offt=offt, pa=pa, pb=pb,
                    wa=wa, wb=wb, comment=comment)
        self.add_element_object(elem)

    def add_pbar(self, pid, mid, A=0., i1=0., i2=0., i12=0., j=0., nsm=0.,
                 c1=0., c2=0., d1=0., d2=0., e1=0., e2=0.,
                 f1=0., f2=0., k1=1.e8, k2=1.e8, comment=''):
        prop = PBAR(pid, mid, A=A, i1=i1, i2=i2, i12=i12, j=j, nsm=nsm,
                    c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2,
                    f1=f1, f2=f2, k1=k1, k2=k2,
                    comment=comment)
        self.add_property_object(prop)

    def add_pbarl(self, pid, mid, group, Type, dim, nsm=0., comment=''):
        prop = PBARL(pid, mid, group, Type, dim, nsm=nsm, comment=comment)
        self.add_property_object(prop)

    def add_cbeam(self, eid, pid, ga, gb, x, g0, is_offt, offt, bit, pa, pb, wa,
                  wb, sa, sb, comment=''):
        elem = CBEAM(eid, pid, ga, gb, x, g0, is_offt, offt, bit, pa, pb, wa,
                     wb, sa, sb, comment=comment)
        self.add_element_object(elem)

    def add_pbeam(self, pid, mid, xxb, so, area, i1, i2, i12, j, nsm,
                  c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, s1, s2,
                  nsia, nsib, cwa, cwb, m1a, m2a, m1b,
                  m2b, n1a, n2a, n1b, n2b, comment=''):
        prop = PBEAM(pid, mid, xxb, so, area, i1, i2, i12, j, nsm, c1, c2, d1,
                     d2, e1, e2, f1, f2, k1, k2, s1, s2,
                     nsia, nsib, cwa, cwb, m1a, m2a, m1b,
                     m2b, n1a, n2a, n1b, n2b, comment=comment)
        self.add_property_object(prop)

    def add_pbeaml(self, pid, mid, group, Type, xxb, so, dims, nsm, comment=''):
        prop = PBEAML(pid, mid, group, Type, xxb, so, dims, nsm, comment=comment)
        self.add_property_object(prop)

    def add_cshear(self, eid, pid, nids, comment=''):
        elem = CSHEAR(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_pshear(self, pid, t, mid, nsm, f1, f2, comment=''):
        prop = PSHEAR(pid, t, mid, nsm, f1, f2, comment=comment)
        self.add_property_object(prop)

    def add_ctria3(self, eid, pid, nids, zOffset, thetaMcid=0.0, TFlag=0,
                   T1=1.0, T2=1.0, T3=1.0, comment=''):
        elem = CTRIA3(eid, pid, nids, zOffset, thetaMcid=thetaMcid, TFlag=TFlag,
                      T1=T1, T2=T2, T3=T3, comment=comment)
        self.add_element_object(elem)

    def add_cquad4(self, eid, pid, nids, thetaMcid=0.0, zOffset=0., TFlag=0,
                   T1=1.0, T2=1.0, T3=1.0, T4=1.0,
                   comment=''):
        elem = CQUAD4(eid, pid, nids, thetaMcid=thetaMcid, zOffset=zOffset, TFlag=TFlag,
                      T1=T1, T2=T2, T3=T3, T4=T4,
                      comment=comment)
        self.add_element_object(elem)

    def add_ctria6(self, eid, pid, nids, theta_mcid, zOffset, TFlag, T1, T2, T3,
                   comment=''):
        elem = CTRIA6(eid, pid, nids, theta_mcid, zOffset, TFlag, T1, T2, T3,
                      comment=comment)
        self.add_element_object(elem)

    def add_cquad8(self, eid, pid, nids, T1, T2, T3, T4, thetaMcid, zOffset,
                   TFlag, comment=''):
        elem = CQUAD8(eid, pid, nids, T1, T2, T3, T4, thetaMcid, zOffset,
                      TFlag, comment=comment)
        self.add_element_object(elem)

    def add_cquad(self, eid, pid, nids, comment=''):
        elem = CQUAD(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_pshell(self, pid, mid1=None, t=None, mid2=None, twelveIt3=1.0,
                   mid3=None, tst=0.833333, nsm=0.0,
                   z1=None, z2=None, mid4=None,
                   comment=''):
        prop = PSHELL(pid, mid1=mid1, t=t, mid2=mid2, twelveIt3=twelveIt3,
                      mid3=mid3, tst=tst, nsm=nsm,
                      z1=z1, z2=z2, mid4=mid4,
                      comment=comment)
        self.add_property_object(prop)

    def add_pcomp(self, pid, mids, thicknesses, thetas, souts, nsm=0., sb=0.,
                  ft=None, TRef=0., ge=0., lam=None,
                  z0=None, comment=''):
        prop = PCOMP(pid, mids, thicknesses, thetas, souts, nsm=nsm, sb=sb,
                     ft=ft, TRef=TRef, ge=ge, lam=lam,
                     z0=z0, comment=comment)
        self.add_property_object(prop)

    def add_pcompg(self, pid, global_ply_ids, mids, thicknesses, thetas, souts,
                   nsm, sb, ft, TRef, ge, lam, z0,
                   comment=''):
        prop = PCOMPG(pid, global_ply_ids, mids, thicknesses, thetas, souts,
                      nsm, sb, ft, TRef, ge, lam, z0,
                      comment=comment)
        self.add_property_object(prop)

    def add_pcomps(self, pid, cordm, psdir, sb, nb, tref, ge, global_ply_ids,
                   mids, thicknesses, thetas,
                   failure_theories,
                   interlaminar_failure_theories,
                   souts, comment=''):
        prop = PCOMPS(pid, cordm, psdir, sb, nb, tref, ge, global_ply_ids,
                      mids, thicknesses, thetas,
                      failure_theories,
                      interlaminar_failure_theories,
                      souts, comment=comment)
        self.add_property_object(prop)

    def add_plplane(self, pid, mid, cid=0, stress_strain_output_location='GRID',
                    comment=''):
        prop = PLPLANE(pid, mid, cid=cid,
                       stress_strain_output_location=stress_strain_output_location,
                       comment=comment)
        self.add_property_object(prop)

    def add_pplane(self, pid, mid, t=0., nsm=0., formulation_option=0,
                   comment=''):
        prop = PPLANE(pid, mid, t=t, nsm=nsm, formulation_option=formulation_option,
                      comment=comment)
        self.add_property_object(prop)

    def add_ctetra(self, eid, pid, nids, comment=''):
        elem = CTETRA(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_cpenta(self, eid, pid, nids, comment=''):
        elem = CPENTA(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_chexa(self, eid, pid, nids, comment=''):
        elem = CHEXA(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_pyram(self, eid, pid, nids, comment=''):
        elem = CPYRAM(eid, pid, nids, comment=comment)
        self.add_element_object(elem)

    def add_psolid(self, pid, mid, cordm=0, integ=None, stress=None, isop=None,
                   fctn='SMECH', comment=''):
        prop = PSOLID(pid, mid, cordm=cordm, integ=integ, stress=stress, isop=isop,
                      fctn=fctn, comment=comment)
        self.add_property_object(prop)

    def add_plsolid(self, pid, mid, stress_strain='GRID', ge=0., comment=''):
        prop = PLSOLID(pid, mid, stress_strain=stress_strain, ge=ge, comment=comment)
        self.add_property_object(prop)

    def add_mat1(self, mid, E, G, nu, rho=0.0, a=0.0, TRef=0.0, ge=0.0, St=0.0,
                 Sc=0.0, Ss=0.0, Mcsid=0, comment=''):
        mat = MAT1(mid, E, G, nu, rho=rho, a=a, TRef=TRef, ge=ge, St=St,
                   Sc=Sc, Ss=Ss, Mcsid=Mcsid, comment=comment)
        self.add_structural_material_object(mat)

    def add_mat2(self, mid, G11, G12, G13, G22, G23, G33, rho, a1, a2, a3,
                 TRef=0., ge=0., St=None, Sc=None,
                 Ss=None, Mcsid=None, comment=''):
        mat = MAT2(mid, G11, G12, G13, G22, G23, G33, rho, a1, a2, a3,
                   TRef=TRef, ge=ge, St=St, Sc=Sc,
                   Ss=Ss, Mcsid=Mcsid, comment=comment)
        self.add_structural_material_object(mat)

    def add_mat3(self, mid, ex, eth, ez, nuxth, nuthz, nuzx, rho=0.0, gzx=None,
                 ax=0., ath=0., az=0., TRef=0., ge=0.,
                 comment=''):
        mat = MAT3(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho=rho, gzx=gzx,
                   ax=ax, ath=ath, az=az, TRef=TRef, ge=ge,
                   comment=comment)
        self.add_structural_material_object(mat)

    def add_mat4(self, mid, k, cp=0.0, rho=1.0, H=None, mu=None, hgen=1.0,
                 refEnthalpy=None, tch=None, tdelta=None,
                 qlat=None, comment=''):
        mat = MAT4(mid, k, cp=cp, rho=rho, H=H, mu=mu, hgen=hgen,
                   refEnthalpy=refEnthalpy, tch=tch, tdelta=tdelta,
                   qlat=qlat, comment=comment)
        self.add_thermal_material_object(mat)

    def add_mat5(self, mid, kxx=0., kxy=0., kxz=0., kyy=0., kyz=0., kzz=0., cp=0.,
                 rho=1., hgen=1., comment=''):
        mat = MAT5(mid, kxx=kxx, kxy=kxy, kxz=kxz, kyy=kyy, kyz=kyz, kzz=kzz, cp=cp,
                   rho=rho, hgen=hgen, comment=comment)
        self.add_thermal_material_object(mat)

    def add_mat8(self, mid, e11, e22, nu12, g12, g1z, g2z, rho=0., a1=0., a2=0.,
                 TRef=0., Xt=0., Xc=None, Yt=0., Yc=None,
                 S=0., ge=0., F12=0., strn=0., comment=''):
        mat = MAT8(mid, e11, e22, nu12, g12, g1z, g2z, rho=rho, a1=a1, a2=a2,
                   TRef=TRef, Xt=Xt, Xc=Xc, Yt=Yt, Yc=Yc,
                   S=S, ge=ge, F12=F12, strn=strn, comment=comment)
        self.add_structural_material_object(mat)

    def add_mat9(self, mid, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26,
                 G33, G34, G35, G36, G44, G45, G46, G55,
                 G56, G66, rho, A, TRef, ge, comment=''):
        mat = MAT9(mid, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26,
                   G33, G34, G35, G36, G44, G45, G46, G55,
                   G56, G66, rho, A, TRef, ge, comment=comment)
        self.add_structural_material_object(mat)

    def add_mat10(self, mid, bulk, rho, c, ge, comment=''):
        mat = MAT10(mid, bulk, rho, c, ge, comment=comment)
        self.add_structural_material_object(mat)

    def add_mat11(self, mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho,
                  a1, a2, a3, TRef, ge, comment=''):
        mat = MAT11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho,
                    a1, a2, a3, TRef, ge, comment=comment)
        self.add_structural_material_object(mat)

    def add_load(self, sid, scale, scale_factors, load_ids, comment=''):
        load = LOAD(sid, scale, scale_factors, load_ids, comment=comment)
        self.add_load_object(load)

    def add_force(self, sid, node, mag, cid=0, xyz=None, comment=''):
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
        cid : int; default=0
            the coordinate system for the load
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the load direction in the cid frame
        comment : str; default=''
            a comment for the card
        """
        load = FORCE(sid, node, mag, cid=cid, xyz=xyz, comment=comment)
        self.add_load_object(load)

    def add_force1(self, sid, node, mag, g1, g2, comment=''):
        load = FORCE1(sid, node, mag, g1, g2, comment=comment)
        self.add_load_object(load)

    def add_force2(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        load = FORCE2(sid, node, mag, g1, g2, g3, g4, comment=comment)
        self.add_load_object(load)

    def add_moment(self, sid, node, cid, mag, xyz, comment=''):
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
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the load direction in the cid frame
        comment : str; default=''
            a comment for the card
        """
        load = MOMENT(sid, node, cid, mag, xyz, comment=comment)
        self.add_load_object(load)

    def add_moment1(self, sid, node, mag, g1, g2, comment=''):
        load = MOMENT1(sid, node, mag, g1, g2, comment=comment)
        self.add_load_object(load)

    def add_moment2(self, sid, node, mag, g1, g2, g3, g4, comment=''):
        load = MOMENT2(sid, node, mag, g1, g2, g3, g4, xyz=None, comment=comment)
        self.add_load_object(load)

    def add_accel(self, sid, cid, N, direction, locs, vals, comment=''):
        load = ACCEL(sid, cid, N, direction, locs, vals, comment=comment)
        self.add_load_object(load)

    def add_accel1(self, sid, cid, scale, N, nodes, comment=''):
        load = ACCEL1(sid, cid, scale, N, nodes, comment=comment)
        self.add_load_object(load)

    def add_grav(self, sid, scale, N, cid=0, mb=0, comment=''):
        load = GRAV(sid, scale, N, cid=cid, mb=mb, comment=comment)
        self.add_load_object(load)

    def add_pload(self, sid, p, nodes, comment=''):
        load = PLOAD(sid, p, nodes, comment=comment)
        self.add_load_object(load)

    def add_pload1(self, sid, eid, Type, scale, x1, p1, x2, p2, comment=''):
        load = PLOAD1(sid, eid, Type, scale, x1, p1, x2, p2, comment=comment)
        self.add_load_object(load)

    def add_pload2(self, sid, pressure, eids, comment=''):
        load = PLOAD2(sid, pressure, eids, comment=comment)
        self.add_load_object(load)

    def add_pload4(self, sid, eids, pressures, g1=None, g34=None, cid=0,
                   NVector=None, sorl='SURF', ldir='NORM', comment=''):
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
        NVector : (3, ) float ndarray
           the local pressure vector
        sorl : str; default='SURF'
           ???
        ldir : str; default='NORM'
           ???
        comment : str; default=''
            a comment for the card

        TODO: fix the way "pressures" works
        """
        load = PLOAD4(sid, eids, pressures, g1=g1, g34=g34, cid=cid,
                      NVector=NVector, sorl=sorl,
                      ldir=ldir, comment=comment)
        self.add_load_object(load)

    def add_spc(self, conid, gids, constraints, enforced, comment=''):
        spc = SPC(conid, gids, constraints, enforced, comment=comment)
        self.add_constraint_spc_object(spc)

    def add_spc1(self, conid, constraints, nodes, comment=''):
        spc = SPC1(conid, constraints, nodes, comment=comment)
        self.add_constraint_spc_object(spc)

    def add_spcd(self, sid, gids, constraints, enforced, comment=''):
        spc = SPCD(sid, gids, constraints, enforced, comment=comment)
        self.add_load_object(spc)

    def add_spcadd(self, conid, sets, comment=''):
        spcadd = SPCADD(conid, sets, comment=comment)
        self.add_constraint_spc_object(spcadd)

    def add_mpc(self, conid, gids, constraints, enforced, comment=''):
        mpc = MPC(conid, gids, constraints, enforced, comment=comment)
        self.add_constraint_mpc_object(mpc)

    def add_mpcadd(self, conid, sets, comment=''):
        mpcadd = MPCADD(conid, sets, comment=comment)
        self.add_constraint_mpc_object(mpcadd)

    def add_suport(self, ids, Cs, comment=''):
        suport = SUPORT(ids, Cs, comment=comment)
        self.add_suport_object(suport)

    def add_suport1(self, conid, ids, Cs, comment=''):
        suport1 = SUPORT1(conid, ids, Cs, comment=comment)
        self.add_suport1_object(suport1)

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
        self.add_aeros_object(aeros)

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
        self.add_aero_object(aero)

    def add_caero1(self, eid, pid, cp, nspan, lspan, nchord, lchord, igid, p1,
                   x12, p4, x43, comment=''):
        caero = CAERO1(eid, pid, cp, nspan, lspan, nchord, lchord, igid, p1,
                       x12, p4, x43, comment=comment)
        self.add_caero_object(caero)

    def add_caero2(self, eid, pid, igid, p1, x12, cp=0, nsb=0, nint=0, lsb=0,
                   lint=0, comment=''):
        caero = CAERO2(eid, pid, igid, p1, x12, cp=cp, nsb=nsb, nint=nint, lsb=lsb,
                       lint=lint, comment=comment)
        self.add_caero_object(caero)

    def add_caero3(self, eid, pid, cp, list_w, list_c1, list_c2, p1, x12, p4,
                   x43, comment=''):
        caero = CAERO3(eid, pid, cp, list_w, list_c1, list_c2, p1, x12, p4,
                       x43, comment=comment)
        self.add_caero_object(caero)

    def add_caero4(self, eid, pid, cp, nspan, lspan, p1, x12, p4, x43,
                   comment=''):
        caero = CAERO4(eid, pid, cp, nspan, lspan, p1, x12, p4, x43,
                       comment=comment)
        self.add_caero_object(caero)

    def add_caero5(self, eid, pid, p1, x12, p4, x43, cp=0, nspan=0, lspan=0,
                   ntheory=0, nthick=0, comment=''):
        caero = CAERO5(eid, pid, p1, x12, p4, x43, cp=cp, nspan=nspan, lspan=lspan,
                       ntheory=ntheory, nthick=nthick, comment=comment)
        self.add_caero_object(caero)

    def add_paero1(self, pid, Bi=None, comment=''):
        paero = PAERO1(pid, Bi=Bi, comment=comment)
        self.add_paero_object(paero)

    def add_paero2(self, pid, orient, width, AR, lrsb, lrib, lth1, lth2, thi,
                   thn, comment=''):
        paero = PAERO2(pid, orient, width, AR, lrsb, lrib, lth1, lth2, thi,
                       thn, comment=comment)
        self.add_paero_object(paero)

    def add_paero3(self, pid, nbox, ncontrol_surfaces, x, y, comment=''):
        paero = PAERO3(pid, nbox, ncontrol_surfaces, x, y, comment=comment)
        self.add_paero_object(paero)

    def add_paero4(self, pid, docs, caocs, gapocs, cla=0, lcla=0,
                   circ=0, lcirc=0, comment=''):
        paero = PAERO4(pid, docs, caocs, gapocs, cla=cla, lcla=lcla,
                       circ=circ, lcirc=lcirc, comment=comment)
        self.add_paero_object(paero)

    def add_paero5(self, pid, caoci, nalpha=0, lalpha=0, nxis=0, lxis=0,
                   ntaus=0, ltaus=0, comment=''):
        paero = PAERO5(pid, caoci, nalpha=nalpha, lalpha=lalpha, nxis=nxis, lxis=lxis,
                       ntaus=ntaus, ltaus=ltaus, comment=comment)
        self.add_paero_object(paero)

    def add_spline1(self, eid, caero, box1, box2, setg, dz=0., method='IPS',
                    usage='BOTH', nelements=10,
                    melements=10, comment=''):
        spline = SPLINE1(eid, caero, box1, box2, setg, dz=dz, method=method,
                         usage=usage, nelements=nelements, melements=melements,
                         comment=comment)
        self.add_spline_object(spline)

    def add_spline2(self, eid, caero, id1, id2, setg, dz, dtor, cid, dthx,
                    dthy, usage, comment=''):
        spline = SPLINE2(eid, caero, id1, id2, setg, dz, dtor, cid, dthx,
                         dthy, usage, comment=comment)
        self.add_spline_object(spline)

    def add_spline3(self, eid, caero, box_id, components, nids,
                    displacement_components,
                    coeffs, usage='BOTH', comment=''):
        spline = SPLINE3(eid, caero, box_id, components, nids,
                         displacement_components,
                         coeffs, usage=usage,
                         comment=comment)
        self.add_spline_object(spline)

    def add_spline4(self, eid, caero, aelist, setg, dz, method, usage,
                    nelements, melements, comment=''):
        spline = SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                         nelements, melements, comment=comment)
        self.add_spline_object(spline)

    def add_spline5(self, eid, caero, aelist, setg, dz, dtor, cid, thx, thy,
                    usage, method, ftype, rcore, comment=''):
        spline = SPLINE5(eid, caero, aelist, setg, dz, dtor, cid, thx, thy,
                         usage, method, ftype, rcore, comment=comment)
        self.add_spline_object(spline)

    def add_trim(self, sid, mach, q, labels, uxs, aeqr=0.0, comment=''):
        trim = TRIM(sid, mach, q, labels, uxs, aeqr=aeqr, comment=comment)
        self.add_trim_object(trim)

    def add_mkaero1(self, machs, reduced_freqs, comment=''):
        mkaero = MKAERO1(machs, reduced_freqs, comment=comment)
        self.add_mkaero_object(mkaero)

    def add_mkaero2(self, machs, reduced_freqs, comment=''):
        mkaero = MKAERO2(machs, reduced_freqs, comment=comment)
        self.add_mkaero_object(mkaero)

    def add_gust(self, sid, dload, wg, x0, V, comment=''):
        gust = GUST(sid, dload, wg, x0, V, comment=comment)
        self.add_gust_object(gust)

    def add_eigr(self, sid, method, f1, f2, ne, nd, norm, G, C, comment=''):
        method = EIGR(sid, method, f1, f2, ne, nd, norm, G, C, comment=comment)
        self.add_method_object(method)

    def add_eigrl(self, sid, v1, v2, nd, msglvl, maxset, shfscl, norm, options,
                  values, comment=''):
        method = EIGRL(sid, v1, v2, nd, msglvl, maxset, shfscl, norm, options,
                       values, comment=comment)
        self.add_method_object(method)

    def add_eigb(self, sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                 comment=''):
        method = EIGB(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                      comment=comment)
        self.add_method_object(method)

    def add_eigc(self, sid, method, norm, G, C, E, ndo, mblkszs=None,
                 iblkszs=None, ksteps=None,
                 NJIs=None, alphaAjs=None,
                 omegaAjs=None, alphaBjs=None,
                 omegaBjs=None, LJs=None, NEJs=None,
                 NDJs=None, shift_r1=None,
                 shift_i1=None, isrr_flag=None,
                 nd1=None, comment=''):
        method = EIGC(sid, method, norm, G, C, E, ndo, mblkszs=mblkszs,
                      iblkszs=iblkszs, ksteps=ksteps,
                      NJIs=NJIs, alphaAjs=alphaAjs,
                      omegaAjs=omegaAjs, alphaBjs=alphaBjs,
                      omegaBjs=omegaBjs, LJs=LJs, NEJs=NEJs,
                      NDJs=NDJs, shift_r1=shift_r1,
                      shift_i1=shift_i1, isrr_flag=isrr_flag,
                      nd1=nd1, comment=comment)
        self.add_cmethod_object(method)

    def add_eigp(self, sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=''):
        method = EIGP(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment=comment)
        self.add_cmethod_object(method)

    def add_set1(self, sid, ids, is_skin=False, comment=''):
        set_obj = SET1(sid, ids, is_skin=is_skin, comment=comment)
        self.add_set_object(set_obj)

    def add_set3(self, sid, desc, ids, comment=''):
        set_obj = SET3(sid, desc, ids, comment=comment)
        self.add_set_object(set_obj)

    def add_aset(self, ids, components, comment=''):
        aset = ASET(ids, components, comment=comment)
        self.add_aset_object(aset)

    def add_aset1(self, components, ids, comment=''):
        aset = ASET1(components, ids, comment=comment)
        self.add_aset_object(aset)

    def add_bset(self, ids, components, comment=''):
        bset = BSET(ids, components, comment=comment)
        self.add_bset_object(bset)

    def add_bset1(self, components, ids, comment=''):
        bset = BSET1(components, ids, comment=comment)
        self.add_bset_object(bset)

    def add_cset(self, ids, components, comment=''):
        cset = CSET(ids, components, comment=comment)
        self.add_cset_object(cset)

    def add_cset1(self, ids, components, comment=''):
        cset = CSET1(ids, components, comment=comment)
        self.add_cset_object(cset)

    def add_qset(self, ids, components, comment=''):
        qset = QSET(ids, components, comment=comment)
        self.add_qset_object(qset)

    def add_qset1(self, components, ids, comment=''):
        qset = QSET1(components, ids, comment=comment)
        self.add_qset_object(qset)

    def add_uset(self, name, components, ids, comment=''):
        uset = USET(name, components, ids, comment=comment)
        self.add_uset_object(uset)

    def add_uset1(self, name, components, ids, comment=''):
        uset = USET1(name, components, ids, comment=comment)
        self.add_uset_object(uset)

    def add_flutter(self, sid, method, density, mach, reduced_freq_velocity,
                    imethod='L', nvalue=None,
                    omax=None, epsilon=None,
                    comment=''):
        flutter = FLUTTER(sid, method, density, mach, reduced_freq_velocity,
                          imethod=imethod, nvalue=nvalue,
                          omax=omax, epsilon=epsilon,
                          comment=comment)
        self.add_flutter_object(flutter)

    def add_flfact(self, sid, factors, comment=''):
        flfact = FLFACT(sid, factors, comment=comment)
        self.add_flfact_object(flfact)

    def add_aecomp(self, name, list_type, lists, comment=''):
        aecomp = AECOMP(name, list_type, lists, comment=comment)
        self.add_aecomp_object(aecomp)

    def add_aestat(self, id, label, comment=''):
        aestat = AESTAT(id, label, comment=comment)
        self.add_aestat_object(aestat)

    def add_aelink(self, id, label, independent_labels, Cis, comment=''):
        aelink = AELINK(id, label, independent_labels, Cis, comment=comment)
        self.add_aelink_object(aelink)

    def add_aelist(self, sid, elements, comment):
        aelist = AELIST(sid, elements, comment=comment)
        self.add_aelist_object(aelist)

    def add_aefact(self, sid, Di, comment=''):
        aefact = AEFACT(sid, Di, comment=comment)
        self.add_aefact_object(aefact)

    def add_diverg(self, sid, nroots, machs, comment=''):
        diverg = DIVERG(sid, nroots, machs, comment=comment)
        self.add_diverg_object(diverg)

    def add_csschd(self, sid, aesid, lschd, lalpha=None, lmach=None,
                   comment=''):
        csschd = CSSCHD(sid, aesid, lschd, lalpha=lalpha, lmach=lmach,
                        comment=comment)
        self.add_csschd_object(csschd)

    def add_aesurf(self, aesid, label, cid1, alid1, cid2=None, alid2=None,
                   eff=1.0, ldw='LDW', crefc=1.0,
                   crefs=1.0, pllim=-np.pi/2.,
                   pulim=np.pi/2., hmllim=None,
                   hmulim=None, tqllim=None,
                   tqulim=None, comment=''):
        """
        Creates an AESURF card

        Parameters
        ----------
        aesid : int
            controller number
        label : str
            controller name
        cid1 : int
            coordinate system id for primary control surface
        alid1 : int
            AELIST id for primary control surface
        cid2 : int; default=None
            coordinate system id for secondary control surface
        alid2 : int; default=None
            AELIST id for secondary control surface
        eff : float; default=1.0
            Control surface effectiveness
        ldw : str; default='LDW'
            Linear downwash flag;  ['LDW', 'NODLW']
        crefc : float; default=1.0
            reference chord for the control surface
        crefs : float; default=1.0
            reference area for the control surface
        pllim : float; default=-pi/2
            Lower deflection limits for the control surface in radians
        pulim : float; default=pi/2
            Upper deflection limits for the control surface in radians
        hmllim : float; default=None
            Lower hinge moment limits for the control surface in
            force-length units
        hmulim : float; default=None
            Upper hinge moment limits for the control surface in
            force-length units
        tqllim : int; default=None
            Set identification numbers of TABLEDi entries that provide the
            lower deflection limits for the control surface as a function
            of the dynamic pressure
        tqulim : int; default=None
            Set identification numbers of TABLEDi entries that provide the
            upper deflection limits for the control surface as a function
            of the dynamic pressure
        comment : str; default=''
            a comment for the card
        """
        aesurf = AESURF(aesid, label, cid1, alid1, cid2=cid2, alid2=alid2,
                        eff=eff, ldw=ldw, crefc=crefc,
                        crefs=crefs, pllim=pllim,
                        pulim=pulim, hmllim=hmllim,
                        hmulim=hmulim, tqllim=tqllim,
                        tqulim=tqulim, comment=comment)
        self.add_aesurf_object(aesurf)

    def add_aesurfs(self, aesid, label, list1, list2, comment=''):
        aesurfs = AESURFS(aesid, label, list1, list2, comment=comment)
        self.add_aesurfs_object(aesurfs)

    def add_aeparm(self, id, label, units, comment=''):
        aeparm = AEPARM(id, label, units, comment=comment)
        self.add_aeparm_object(aeparm)

    def add_tabled1(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        table = TABLED1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)
        self.add_table_object(table)

    def add_tabled2(self, tid, x1, x, y, comment=''):
        table = TABLED2(tid, x1, x, y, comment=comment)
        self.add_table_object(table)

    def add_tabled3(self, tid, x1, x2, x, y, comment=''):
        table = TABLED3(tid, x1, x2, x, y, comment=comment)
        self.add_table_object(table)

    def add_tabled4(self, tid, x1, x2, x3, x4, a, comment=''):
        table = TABLED4(tid, x1, x2, x3, x4, a, comment=comment)
        self.add_table_object(table)

    def add_tablem1(self, tid, x, y, comment=''):
        table = TABLEM1(tid, x, y, comment=comment)
        self.add_table_object(table)

    def add_tablem2(self, tid, x1, x, y, comment=''):
        table = TABLEM2(tid, x1, x, y, comment=comment)
        self.add_table_object(table)

    def add_tablem3(self, tid, x1, x2, x, y, comment=''):
        table = TABLEM3(tid, x1, x2, x, y, comment=comment)
        self.add_table_object(table)

    def add_tablem4(self, tid, x1, x2, x3, x4, a, comment=''):
        table = TABLEM4(tid, x1, x2, x3, x4, a, comment=comment)
        self.add_table_object(table)

    def add_tables1(self, tid, Type, x, y, comment=''):
        table = TABLES1(tid, Type, x, y, comment=comment)
        self.add_table_object(table)

    def add_tablest(self, tid, x, y, comment=''):
        table = TABLEST(tid, x, y, comment=comment)
        self.add_table_object(table)

    def add_freq(self, sid, freqs, comment=''):
        freq = FREQ(sid, freqs, comment=comment)
        self.add_freq_object(freq)

    def add_freq1(self, sid, f1, df, ndf, comment=''):
        freq = FREQ1(sid, f1, df, ndf, comment=comment)
        self.add_freq_object(freq)

    def add_freq2(self, sid, f1, f2, ndf=1, comment=''):
        freq = FREQ2(sid, f1, f2, ndf=ndf, comment=comment)
        self.add_freq_object(freq)

    def add_freq4(self, sid, f1, f2, fspread, nfm, comment=''):
        freq = FREQ4(sid, f1, f2, fspread, nfm, comment=comment)
        self.add_freq_object(freq)

    def add_rrod(self, eid, ga, gb, cma, cmb, alpha, comment=''):
        elem = RROD(eid, ga, gb, cma, cmb, alpha, comment=comment)
        self.add_rigid_element_object(elem)

    def add_rbe1(self, eid, Gni, Cni, Gmi, Cmi, alpha, comment=''):
        elem = RBE1(eid, Gni, Cni, Gmi, Cmi, alpha, comment=comment)
        self.add_rigid_element_object(elem)

    def add_rbe2(self, eid, gn, cm, Gmi, alpha=0.0, comment=''):
        elem = RBE2(eid, gn, cm, Gmi, alpha=alpha, comment=comment)
        self.add_rigid_element_object(elem)

    def add_rbe3(self, eid, refgrid, refc, weights, comps, Gijs, Gmi=None,
                 Cmi=None, alpha=0.0, comment=''):
        elem = RBE3(eid, refgrid, refc, weights, comps, Gijs, Gmi=Gmi,
                    Cmi=Cmi, alpha=alpha, comment=comment)
        self.add_rigid_element_object(elem)

    def add_rbar(self, eid, ga, gb, cna, cnb, cma, cmb, alpha=0., comment=''):
        elem = RBAR(eid, ga, gb, cna, cnb, cma, cmb, alpha=alpha, comment=comment)
        self.add_rigid_element_object(elem)

    def add_rbar1(self, eid, ga, gb, cb, alpha, comment=''):
        elem = RBAR1(eid, ga, gb, cb, alpha, comment=comment)
        self.add_rigid_element_object(elem)

    def add_rspline(self, eid, nids, components, diameter_ratio=0.1, comment=''):
        elem = RSPLINE(eid, nids, components, diameter_ratio=diameter_ratio, comment=comment)
        self.add_rigid_element_object(elem)

    def add_tf(self, sid, nid0, c, b0, b1, b2, nids, components, a, comment=''):
        tf = TF(sid, nid0, c, b0, b1, b2, nids, components, a, comment=comment)
        self.add_tf_object(tf)

    def add_desvar(self, desvar_id, label, xinit, xlb, xub, delx=None,
                   ddval=None, comment=''):
        desvar = DESVAR(desvar_id, label, xinit, xlb, xub, delx=delx,
                        ddval=ddval, comment=comment)
        self.add_desvar_object(desvar)

    def add_dresp1(self, dresp_id, label, response_type, property_type, region,
                   atta, attb, atti, validate=True, comment=''):
        dresp = DRESP1(dresp_id, label, response_type, property_type, region,
                       atta, attb, atti, validate=validate, comment=comment)
        self.add_dresp_object(dresp)

    def add_dresp2(self, dresp_id, label, dequation, region, method, c1, c2, c3,
                   params, comment=''):
        dresp = DRESP2(dresp_id, label, dequation, region, method, c1, c2, c3,
                       params, comment=comment)
        self.add_dresp_object(dresp)

    def add_dresp3(self, dresp_id, label, group, Type, region, params,
                   comment=''):
        dresp = DRESP3(dresp_id, label, group, Type, region, params,
                       comment=comment)
        self.add_dresp_object(dresp)

    def add_dvcrel1(self, oid, Type, eid, cp_name, cp_min, cp_max, dvids,
                    coeffs, c0=0., validate=True, comment=''):
        dvcrel = DVCREL1(oid, Type, eid, cp_name, cp_min, cp_max, dvids,
                         coeffs, c0=c0, validate=validate, comment=comment)
        self.add_dvcrel_object(dvcrel)

    def add_dvcrel2(self, oid, Type, eid, cp_name, cp_min, cp_max, deqation,
                    dvids, labels,
                    validate=False, comment=''):
        dvcrel = DVCREL2(oid, Type, eid, cp_name, cp_min, cp_max, deqation,
                         dvids, labels, validate=validate, comment=comment)
        self.add_dvcrel_object(dvcrel)

    def add_dvprel1(self, oid, Type, pid, pname_fid, p_min, p_max, dvids,
                    coeffs, c0=0.0, validate=False, comment=''):
        dvprel = DVPREL1(oid, Type, pid, pname_fid, p_min, p_max, dvids,
                         coeffs, c0=c0, validate=validate, comment=comment)
        self.add_dvprel_object(dvprel)

    def add_dvprel2(self, oid, Type, pid, pname_fid, p_min, p_max, deqation,
                    dvids, labels, validate=True, comment=''):
        dvprel = DVPREL2(oid, Type, pid, pname_fid, p_min, p_max, deqation,
                         dvids, labels, validate=validate, comment=comment)
        self.add_dvprel_object(dvprel)

    def add_dvmrel1(self, oid, Type, mid, mp_name, mp_min, mp_max, dvids,
                    coeffs, c0=0., validate=True, comment=''):
        dvmrel = DVMREL1(oid, Type, mid, mp_name, mp_min, mp_max, dvids,
                         coeffs, c0=c0, validate=validate, comment=comment)
        self.add_dvmrel_object(dvmrel)

    def add_dvmrel2(self, oid, Type, mid, mp_name, mp_min, mp_max, deqation,
                    dvids, labels, validate=True, comment=''):
        dvmrel = DVMREL2(oid, Type, mid, mp_name, mp_min, mp_max, deqation,
                         dvids, labels, validate=validate, comment=comment)
        self.add_dvmrel_object(dvmrel)

    def add_dvgrid(self, dvid, nid, dxyz, cid=0, coeff=1.0, comment=''):
        dvgrid = DVGRID(dvid, nid, dxyz, cid=cid, coeff=coeff, comment=comment)
        self.add_dvgrid_object(dvgrid)

    def add_ddval(self, oid, ddvals, comment=''):
        ddval = DDVAL(oid, ddvals, comment=comment)
        self.add_ddval_object(ddval)

    def add_dlink(self, oid, ddvid, c0, cmult, IDv, Ci, comment=''):
        dlink = DLINK(oid, ddvid, c0, cmult, IDv, Ci, comment=comment)
        self.add_dlink_object(dlink)

    def add_dconstr(self, oid, dresp_id, lid=1.e20, uid=1.e20, lowfq=0.,
                    highfq=1.e20, comment=''):
        dconstr = DCONSTR(oid, dresp_id, lid=lid, uid=uid, lowfq=lowfq,
                          highfq=highfq, comment=comment)
        self.add_dconstr_object(dconstr)

    def add_dconadd(self, oid, dconstrs, comment=''):
        dconadd = DCONADD(oid, dconstrs, comment=comment)
        self.add_dconstr_object(dconadd)

    def add_doptprm(self, params, comment=''):
        doptprm = DOPTPRM(params, comment=comment)
        self.add_doptprm_object(doptprm)

    def add_monpnt1(self, name, label, axes, comp, xyz, cp=0, cd=None,
                    comment=''):
        monitor_point = MONPNT1(name, label, axes, comp, xyz, cp=cp, cd=cd,
                                comment=comment)
        self.add_monpnt_object(monitor_point)

    def add_monpnt2(self, name, label, table, Type, nddl_item, eid,
                    comment=''):
        monitor_point = MONPNT2(name, label, table, Type, nddl_item, eid,
                                comment=comment)
        self.add_monpnt_object(monitor_point)

    def add_monpnt3(self, name, label, axes, grid_set, elem_set, xyz,
                    cp=0, cd=None, xflag=None, comment=''):
        monitor_point = MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                                cp=cp, cd=cd, xflag=xflag, comment=comment)
        self.add_monpnt_object(monitor_point)
