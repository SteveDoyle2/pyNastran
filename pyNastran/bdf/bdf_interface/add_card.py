"""defines various methods to add cards"""
# pylint: disable=R0913, R0914, C0103
from __future__ import print_function

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
        self.add_node(grid)

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

    def add_param(self, key, values, comment=''):
        """
        Creates a PARAM card.

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
        comment : str
            a comment for the card
        """
        prop = PROD(pid, mid, A, j=j, c=c, nsm=nsm, comment=comment)
        self.add_property_object(prop)

    #def add_ctube(self):
    #def add_ptube(self):

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

    def add_pshell(self, pid, mid1=None, t=None, mid2=None, twelveIt3=1.0,
                   mid3=None, tst=0.833333, nsm=0.0,
                   z1=None, z2=None, mid4=None,
                   comment=''):
        prop = PSHELL(pid, mid1=mid1, t=t, mid2=mid2, twelveIt3=twelveIt3,
                      mid3=mid3, tst=tst, nsm=nsm,
                      z1=z1, z2=z2, mid4=mid4,
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

    def add_mat1(self, mid, E, G, nu, rho=0.0, a=0.0, TRef=0.0, ge=0.0, St=0.0,
                 Sc=0.0, Ss=0.0, Mcsid=0, comment=''):
        mat = MAT1(mid, E, G, nu, rho=rho, a=a, TRef=TRef, ge=ge, St=St,
                   Sc=Sc, Ss=Ss, Mcsid=Mcsid, comment=comment)
        self.add_structural_material(mat)

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
        """
        load = FORCE(sid, node, mag, cid=cid, xyz=xyz, comment=comment)
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
        """
        load = MOMENT(sid, node, cid, mag, xyz, comment=comment)
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

        TODO: fix the way "pressures" works
        """
        load = PLOAD4(sid, eids, pressures, g1=g1, g34=g34, cid=cid,
                      NVector=NVector, sorl=sorl,
                      ldir=ldir, comment=comment)
        self.add_load_object(load)

    def add_spc(self, conid, gids, constraints, enforced, comment=''):
        spc = SPC(conid, gids, constraints, enforced, comment='')
        self.add_constraint_spc(spc)

    def add_spc1(self, conid, constraints, nodes, comment=''):
        spc = SPC1(conid, constraints, nodes, comment=comment)
        self.add_constraint_spc(spc)

    def add_spcadd(self, conid, sets, comment=''):
        spcadd = SPCADD(conid, sets, comment=comment)
        self.add_constraint_spc(spcadd)

    def add_mpc(self, conid, gids, constraints, enforced, comment=''):
        mpc = MPC(conid, gids, constraints, enforced, comment=comment)
        self.add_constraint_mpc(mpc)

    def add_mpcadd(self, conid, sets, comment=''):
        mpcadd = MPCADD(conid, sets, comment=comment)
        self.add_constraint_mpc(mpcadd)

    def add_aeros(self, cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0,
                  comment=''):
        aeros = AEROS(cref, bref, sref, acsid=acsid, rcsid=rcsid, sym_xz=sym_xz, sym_xy=sym_xy,
                      comment=comment)
        self.add_aeros_object(aeros)

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

