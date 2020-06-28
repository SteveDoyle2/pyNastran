# pylint: disable=R0902,R0904,R0914
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.cards.nodes import SPOINT, EPOINT
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import (
        CYAX, CYJOIN, AXIF,
        TOPVAR, MPCAX, CORD3G,
        SESUPORT, SEUSET, SEUSET1,
    )
    from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL, GENEL
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
        AXIF, RINGFL,
        AXIC, RINGAX, POINTAX, CCONEAX, PCONEAX, PRESAX, TEMPAX,)
    from pyNastran.bdf.cards.elements.axisymmetric_shells import (
        CTRAX3, CTRAX6, CTRIAX, CTRIAX6, CQUADX, CQUADX4, CQUADX8)
    from pyNastran.bdf.cards.elements.shell import (
        CQUAD, CQUAD4, CQUAD8, CQUADR, CSHEAR,
        CTRIA3, CTRIA6, CTRIAR,
        CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
        CPLSTS3, CPLSTS4, CPLSTS6, CPLSTS8,
        SNORM,
    )
    from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE, PPLANE
    from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
    from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D, PBUSHT #PBUSH2D
    from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                                     CDAMP5)
    from pyNastran.bdf.cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
    from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE
    from pyNastran.bdf.cards.elements.bars import CBAR, CBARAO, CBEAM3, CBEND, BAROR
    from pyNastran.bdf.cards.elements.beam import CBEAM, BEAMOR
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
        LSEQ, SLOAD, DAREA, RFORCE, RFORCE1, SPCD, DEFORM, LOADCYN)
    from pyNastran.bdf.cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
    from pyNastran.bdf.cards.loads.static_loads import (LOAD, CLOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                        FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                        PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                                        GMLOAD)
    from pyNastran.bdf.cards.loads.random_loads import RANDPS, RANDT1

    from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                               MAT8, MAT9, MAT10, MAT11, MAT3D,
                                               MATG, MATHE, MATHP, CREEP, EQUIV,
                                               NXSTRAT)
    from pyNastran.bdf.cards.material_deps import (
        MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, MATS1)

    from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL, MODTRAK
    from pyNastran.bdf.cards.nodes import GRID, GRDSET, SPOINTs, EPOINTs, POINT, SEQGP, GRIDB

    from pyNastran.bdf.cards.aero.aero import (
        AECOMP, AECOMPL, AEFACT, AELINK, AELIST, AEPARM, AESURF, AESURFS,
        CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
        PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
        MONPNT1, MONPNT2, MONPNT3,
        SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5)
    from pyNastran.bdf.cards.aero.static_loads import AESTAT, AEROS, CSSCHD, TRIM, TRIM2, DIVERG
    from pyNastran.bdf.cards.aero.dynamic_loads import AERO, FLFACT, FLUTTER, GUST, MKAERO1, MKAERO2
    #from pyNastran.bdf.cards.aero.zona import (
        #ACOORD, AEROZ, AESURFZ, BODY7, CAERO7, MKAEROZ, PAFOIL7, PANLST1, PANLST3,
        #SEGMESH, SPLINE1_ZONA, SPLINE2_ZONA, SPLINE3_ZONA, TRIMLNK, TRIMVAR, TRIM_ZONA,
        #ZONA)

    from pyNastran.bdf.cards.optimization import (
        DCONADD, DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
        DRESP1, DRESP2, DRESP3,
        DVCREL1, DVCREL2,
        DVMREL1, DVMREL2,
        DVPREL1, DVPREL2,
        DVGRID, DSCREEN)
    from pyNastran.bdf.cards.superelements import (
        RELEASE, SEBNDRY, SEBULK, SECONCT, SEELT, SEEXCLD,
        SELABEL, SELOAD, SELOC, SEMPLN, SENQSET, SETREE,
        CSUPER, CSUPEXT,
    )
    from pyNastran.bdf.cards.bdf_sets import (
        ASET, BSET, CSET, QSET, USET, OMIT,
        ASET1, BSET1, CSET1, QSET1, USET1, OMIT1,
        SET1, SET2, SET3,
        SEBSET, SECSET, SEQSET,  # SEUSET
        SEBSET1, SECSET1, SEQSET1,  # SEUSET1
        SESET,  # SEQSEP,
        RADSET
    )
    from pyNastran.bdf.cards.params import PARAM
    from pyNastran.bdf.cards.dmig import DMIG, DMIAX, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL, DTI
    from pyNastran.bdf.cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, TEMPB3,
                                                   QVOL, QVECT)
    from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                     PHBDY, CONV, CONVM, TEMPBC)
    from pyNastran.bdf.cards.thermal.radiation import RADM, RADBC, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
    from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                                TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                                TABLES1, TABDMP1, TABLEST, TABLEHT, TABLEH1,
                                                TABRND1, TABRNDG,
                                                DTABLE)
    from pyNastran.bdf.cards.contact import (
        BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA, BCONP, BLSEG,
        BFRIC)
    from pyNastran.bdf.cards.parametric.geometry import PSET, PVAL, FEEDGE, FEFACE, GMCURV, GMSURF

class AddMethods(BDFAttributes):
    """defines methods to add card objects to the BDF"""
    def __init__(self) -> None:
        BDFAttributes.__init__(self)

    def _add_dmi_object(self, dmi: DMI, allow_overwrites: bool=False) -> None:
        """adds a DMI object"""
        name = dmi.name
        self.dmi[name] = dmi
        self._type_to_id_map[dmi.type].append(name)

    def _add_dmig_object(self, dmig: DMIG, allow_overwrites: bool=False) -> None:
        """adds a DMIG object"""
        name = dmig.name
        self.dmig[name] = dmig
        self._type_to_id_map[dmig.type].append(name)

    def _add_dmiax_object(self, dmiax: DMIAX, allow_overwrites: bool=False) -> None:
        """adds a DMI object"""
        name = dmiax.name
        self.dmiax[name] = dmiax
        self._type_to_id_map[dmiax.type].append(name)

    def _add_dmij_object(self, dmij: DMIJ, allow_overwrites: bool=False) -> None:
        """adds a DMIJ object"""
        name = dmij.name
        self.dmij[name] = dmij
        self._type_to_id_map[dmij.type].append(name)

    def _add_dmiji_object(self, dmiji: DMIJI, allow_overwrites: bool=False) -> None:
        """adds a DMIJI object"""
        name = dmiji.name
        self.dmiji[name] = dmiji
        self._type_to_id_map[dmiji.type].append(name)

    def _add_dmik_object(self, dmik: DMIK, allow_overwrites: bool=False) -> None:
        """adds a DMIK object"""
        name = dmik.name
        self.dmik[name] = dmik
        self._type_to_id_map[dmik.type].append(name)

    def _add_dti_object(self, dti: DTI, allow_overwrites: bool=False) -> None:
        """adds an DTI object"""
        name = dti.name
        if name == 'UNITS' or name not in self.dti:
            self.dti[name] = dti
            self._type_to_id_map[dti.type].append(name)
        else:
            old_dti = self.dti[name]
            key = list(dti.fields.keys())[0]
            assert key not in old_dti.fields, 'key=%i old_fields=%s fields=%s' % (key, old_dti.fields, dti.fields)
            old_dti.fields[key] = dti.fields[key]

    def _add_param_object(self, param: PARAM, allow_overwrites: bool=False) -> None:
        """adds a PARAM object"""
        key = param.key
        if key in self.params and not allow_overwrites:
            if not param == self.params[key]:
                #if param.key in self.params:
                    #msg = 'key=%s param=%s old_param=%s' % (key, param, self.params[key])
                    #raise KeyError(msg)
                self.log.warning('key=%s param=%s old_param=%s' %
                                 (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
            self._type_to_id_map[param.type].append(key)

    def _add_node_object(self, node: Union[GRID], allow_overwrites: bool=False) -> None:
        """adds a GRID card"""
        key = node.nid
        if key in self.nodes and not allow_overwrites:
            if not node == self.nodes[key]:
                assert node.nid not in self.nodes, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, self.nodes[key], node)
            else:
                #print('GRID was duplicated...nid=%s; node=\n%s' % (key, node))
                pass
        else:
            assert key > 0, 'nid=%s node=%s' % (key, node)
            self.nodes[key] = node
            self._type_to_id_map[node.type].append(key)

    def _add_gridb_object(self, node: GRIDB, allow_overwrites: bool=False) -> None:
        """adds a GRIDB card"""
        key = node.nid
        assert key > 0, 'eid=%s node=%s' % (key, node)
        if key in self.gridb and not allow_overwrites:
            assert node.nid not in self.gridb, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, self.gridb[key], node)
        self.gridb[key] = node
        self._type_to_id_map[node.type].append(key)
        self._is_axis_symmetric = True

    def _add_ringfl_object(self, ringfl: RINGFL, allow_overwrites: bool=False) -> None:
        """adds a RINGFL card"""
        key = ringfl.ringfl
        assert key > 0, 'eid=%s ringfl=%s' % (key, ringfl)
        if key in self.ringfl and not allow_overwrites:
            assert ringfl.ringfl not in self.ringfl, 'ringfl=%s\nold_ringfl=\n%snew_ringfl=\n%s' % (ringfl.ringfl, self.ringfl[key], ringfl)
        self.ringfl[key] = ringfl
        self._type_to_id_map[ringfl.type].append(key)
        self._is_axis_symmetric = True

    def _add_ringax_object(self, ringax: Union[RINGAX, POINTAX], allow_overwrites: bool=False) -> None:
        """adds a RINGAX card"""
        key = ringax.nid
        if key in self.ringaxs and not allow_overwrites:
            if not ringax == self.ringaxs[key]:
                assert ringax.nid not in self.ringaxs, 'nid=%s\nold_ringax=\n%snew_ringax=\n%s' % (ringax.nid, self.ringaxs[key], ringax)
            else:
                #print('RINGAX was duplicated...nid=%s; ringax=\n%s' % (key, ringax))
                pass
        else:
            assert key > 0, 'nid=%s ringax=%s' % (key, ringax)
            self.ringaxs[key] = ringax
            self._type_to_id_map[ringax.type].append(key)
        self._is_axis_symmetric = True

    def _add_seqgp_object(self, seqgp: SEQGP) -> None:
        """adds an SEQGP card"""
        if self.seqgp is None:
            self.seqgp = seqgp
        else:
            self.seqgp.append(seqgp)

    def _add_point_object(self, point: POINT, allow_overwrites: bool=False) -> None:
        """adds a POINT card"""
        key = point.nid
        if key in self.points and not allow_overwrites:
            if not point == self.points[key]:
                assert point.nid not in self.points, 'nid=%s\nold_point=\n%snew_point=\n%s' % (point.nid, self.points[key], point)
            else:
                #print('POINT was duplicated...nid=%s; point=\n%s' % (key, point))
                pass
        else:
            assert key > 0, 'nid=%s point=%s' % (key, point)
            self.points[key] = point
            self._type_to_id_map[point.type].append(key)

    def _add_spoint_object(self, spoints: SPOINTs) -> None:
        """adds an SPOINT card"""
        comment = spoints.comment
        if hasattr(spoints, 'ifile'):
            ifile = spoints.ifile
            for nid in spoints.points:
                if nid in self.spoints:
                    continue
                spoint = SPOINT(nid, comment=comment)
                spoint.ifile = ifile
                comment = ''
                self.spoints[nid] = spoint
        else:
            for nid in spoints.points:
                if nid in self.spoints:
                    continue
                spoint = SPOINT(nid, comment=comment)
                comment = ''
                self.spoints[nid] = spoint

    def _add_epoint_object(self, epoints: EPOINTs) -> None:
        """adds an EPOINT card"""
        comment = epoints.comment
        for nid in epoints.points:
            if nid in self.epoints:
                continue
            epoint = EPOINT(nid, comment=comment)
            comment = ''
            self.epoints[nid] = epoint

    def _add_setree_object(self, setree: SETREE) -> None:
        key = setree.seid
        self.setree[key] = setree
        self._type_to_id_map[setree.type].append(key)
    def _add_senqset_object(self, senqset: SENQSET) -> None:
        key = senqset.set_id
        self.senqset[key] = senqset
        self._type_to_id_map[senqset.type].append(key)
    def _add_sebulk_object(self, sebulk: SEBULK) -> None:
        key = sebulk.seid
        self.sebulk[key] = sebulk
        self._type_to_id_map[sebulk.type].append(key)
    def _add_release_object(self, release: RELEASE) -> None:
        key = release.seid
        self.release[key] = release
        self._type_to_id_map[release.type].append(key)
    def _add_sebndry_object(self, sebndry: SEBNDRY) -> None:
        key = (sebndry.seid_a, sebndry.seid_b)
        self.sebndry[key] = sebndry
    def _add_seloc_object(self, seloc: SELOC) -> None:
        key = seloc.seid
        self.seloc[key] = seloc
        self._type_to_id_map[seloc.type].append(key)
    def _add_sempln_object(self, sempln: SEMPLN) -> None:
        key = sempln.seid
        self.sempln[key] = sempln
        self._type_to_id_map[sempln.type].append(key)

    def _add_seconct_object(self, seconct: SECONCT) -> None:
        key = (seconct.seid_a, seconct.seid_b)
        self.seconct[key] = seconct
        self._type_to_id_map[seconct.type].append(key)
    def _add_selabel_object(self, selabel: SELABEL) -> None:
        key = selabel.seid
        self.selabel[key] = selabel
        self._type_to_id_map[selabel.type].append(key)
    def _add_seexcld_object(self, seexcld: SEEXCLD) -> None:
        key = (seexcld.seid_a, seexcld.seid_b)
        self.seexcld[key] = seexcld
        self._type_to_id_map[seexcld.type].append(key)

    def _add_seelt_object(self, seelt: SEELT) -> None:
        #self.seelt.append(seelt)
        key = seelt.seid
        self.seelt[key] = seelt
        self._type_to_id_map[seelt.type].append(key)
    def _add_seload_object(self, seload: SELOAD) -> None:
        key = seload.seid
        self.seload[key] = seload
        self._type_to_id_map[seload.type].append(key)
    def _add_csuper_object(self, csuper: CSUPER) -> None:
        key = csuper.seid
        self.csuper[key] = csuper
        self._type_to_id_map[csuper.type].append(key)
    def _add_csupext_object(self, csupext: CSUPEXT) -> None:
        key = csupext.seid
        self.csupext[key] = csupext
        self._type_to_id_map[csupext.type].append(key)

    def _add_plotel_object(self, elem: PLOTEL, allow_overwrites: bool=False) -> None:
        """adds an PLOTEL object"""
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if not allow_overwrites:
            if key in self.elements:
                if elem == self.elements[key]:
                    self._duplicate_elements.append(elem)
                    if self._stop_on_duplicate_error:
                        self.pop_parse_errors()
            elif key in self.plotels:
                if not elem == self.plotels[key]:
                    assert elem.eid not in self.plotels, 'eid=%s\nold_element=\n%snew_element=\n%s' % (elem.eid, self.plotels[elem.eid], elem)
        self.plotels[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def _add_element_object(self, elem: Union[CELAS1, CELAS2, CELAS3, CELAS4,
                                              CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5,
                                              CVISC, CBUSH, CBUSH1D, CBUSH2D, CFAST, #CWELD
                                              CGAP, GENEL, CCONEAX,
                                              CROD, CTUBE, CONROD,
                                              CBAR, CBEAM, CBEAM3, CBEND, CSHEAR,
                                              CTRIA3, CTRIA6, CTRIAR,
                                              CQUAD4, CQUAD8, CQUADR, CQUAD,
                                              CTRIAX, CTRIAX6,
                                              CQUADX, CQUADX4, CQUADX8,
                                              CRAC2D, CRAC3D,
                                              CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
                                              CPLSTS3, #CPLSTS4, CPLSTS6, CPLSTS8,
                                              CTETRA4, CTETRA10, CPENTA6, CPENTA15,
                                              CHEXA8, CHEXA20, CPYRAM5, CPYRAM13,
                                              CTRAX3, CTRAX6,
                                              # thermal
                                              CHBDYE, CHBDYG, CHBDYP,
                                              # Nastran 95
                                              CIHEX1, CIHEX2,],
                            allow_overwrites: bool=False) -> None:
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if key in self.elements and not allow_overwrites:
            if not elem == self.elements[key]:
                self._duplicate_elements.append(elem)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.elements[key] = elem
            self._type_to_id_map[elem.type].append(key)

    def _add_ao_object(self, elem_flag: CBARAO, allow_overwrites: bool=False) -> None:
        """adds a CBARAO"""
        key = elem_flag.eid
        assert key > 0, 'eid=%s must be positive; elem_flag=\n%s' % (key, elem_flag)
        if key in self.ao_element_flags and not allow_overwrites:
            if not elem_flag == self.ao_element_flags[key]:
                #self._duplicate_elements.append(elem_flag)
                #if self._stop_on_duplicate_error:
                    #self.pop_parse_errors()
                assert elem_flag.eid not in self.ao_element_flags, 'eid=%s\nold_ao_element=\n%snew_ao_element=\n%s' % (
                    elem_flag.eid, self.ao_element_flags[elem_flag.eid], elem_flag)
        else:
            self.ao_element_flags[key] = elem_flag
            self._type_to_id_map[elem_flag.type].append(key)

    def _add_doptprm_object(self, doptprm: DOPTPRM) -> None:
        """adds a DOPTPRM"""
        self.doptprm = doptprm

    def _add_nsm_object(self, nsm: Union[NSM, NSM1, NSML, NSML1], allow_overwrites: bool=False) -> None:
        """adds a nsm object to a nsm set"""
        key = nsm.sid
        assert key > 0, 'sid=%s must be positive; nsm=\n%s' % (key, nsm)
        if key in self.nsms:
            self.nsms[key].append(nsm)
        else:
            self.nsms[key] = [nsm]
            self._type_to_id_map[nsm.type].append(key)

    def _add_nsmadd_object(self, nsmadd: NSMADD, allow_overwrites: bool=False) -> None:
        """adds a nsmadd object to a nsm set"""
        key = nsmadd.sid
        assert key > 0, 'sid=%s must be positive; nsmadd=\n%s' % (key, nsmadd)
        if key in self.nsmadds:
            self.nsmadds[key].append(nsmadd)
        else:
            self.nsmadds[key] = [nsmadd]
            self._type_to_id_map[nsmadd.type].append(key)

    def _add_mass_object(self, mass: Union[CMASS1, CMASS2, CMASS3, CMASS4,
                                           CONM1, CONM2], allow_overwrites: bool=False) -> None:
        key = mass.eid
        assert key > 0, 'eid=%s must be positive; mass=\n%s' % (key, mass)
        if key in self.masses and not allow_overwrites:
            if not mass == self.masses[key]:
                self._duplicate_masses.append(mass)
        else:
            self.masses[key] = mass
            self._type_to_id_map[mass.type].append(key)

    def _add_damper_object(self, elem, allow_overwrites: bool=False) -> None:
        """.. warning:: can dampers have the same ID as a standard element?"""
        return self._add_element_object(elem, allow_overwrites)

    def _add_rigid_element_object(self, elem: Union[RBAR, RBAR1,
                                                    RBE1, RBE2, RBE3,
                                                    RROD, RSPLINE, RSSCON],
                                  allow_overwrites: bool=False) -> None:
        key = elem.eid
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        if key in self.rigid_elements and not allow_overwrites:
            assert elem.eid not in self.rigid_elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, self.rigid_elements[elem.eid], elem)
        self.rigid_elements[key] = elem
        self._type_to_id_map[elem.type].append(key)

    def _add_thermal_element_object(self, elem: Union[CHBDYE, CHBDYG, CHBDYP]) -> None:
        """same as add_element at the moment..."""
        self._add_element_object(elem)

    def _add_deqatn_object(self, deqatn: DEQATN, allow_overwrites: bool=False) -> None:
        """adds an DEQATN object"""
        key = deqatn.equation_id
        assert key > 0, 'ID=%s deqatn\n%s' % (key, deqatn)
        if key in self.dequations and not allow_overwrites:
            if not deqatn.write_card() == self.dequations[key].write_card():
                assert key not in self.dequations, 'id=%s old_eq=\n%snew_eq=\n%s' % (
                    key, self.dequations[key], deqatn)
        self.dequations[key] = deqatn
        self._type_to_id_map[deqatn.type].append(key)

    def _add_acoustic_property_object(self, prop: PACABS) -> None:
        self._add_property_object(prop)

    def _add_property_object(self, prop: Union[PELAS, PBUSH, PBUSH1D, PDAMP, PDAMP5, # PBUSH2D,
                                               PFAST, PVISC, PGAP, PRAC2D, PRAC3D, # PWELD
                                               PROD, PTUBE,
                                               PBAR, PBARL, PBRSECT, PCONEAX,
                                               PBEAM, PBEAML, PBCOMP, PBMSECT,
                                               PBEND, PBEAM3, PPLANE, PSHEAR,
                                               PSHELL, PCOMP, PCOMPG, PLPLANE,
                                               PSOLID, PLSOLID, PIHEX, PCOMPS],
                             allow_overwrites: bool=False) -> None:
        """
        adds one of the following objects:
          PELAS, PBUSH, PBUSH1D, PBUSH2D, PDAMP,
          PROD, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
          PSHELL, PCOMP, PCOMPG,
          PSOLID, PLSOLID
        """
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        if key in self.properties and not allow_overwrites:
            if not prop == self.properties[key]:
                self._duplicate_properties.append(prop)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.properties[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_property_mass_object(self, prop: PMASS, allow_overwrites: bool=False) -> None:
        """adds an PMASS object"""
        key = prop.pid
        if key in self.properties_mass and not allow_overwrites:
            if not prop == self.properties_mass[key]:
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key,self.properties_mass[key],prop))
                assert key not in self.properties_mass, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, self.properties_mass[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.properties_mass[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_dtable_object(self, dtable: DTABLE, allow_overwrites: bool=False) -> None:
        """adds an DTABLE object"""
        if self.dtable is not None:
            if not dtable == self.dtable:
                raise RuntimeError('DTABLE cannot be overwritten\nold:\n%s\nnew:\n%s',
                                   self.dtable, dtable)
        else:
            self.dtable = dtable
            #self._type_to_id_map[dtable.type].append(1)

    def _add_bcrpara_object(self, card: BCRPARA, allow_overwrites: bool=False) -> None:
        """adds an BCRPARA object"""
        key = card.crid
        self.bcrparas[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bctadd_object(self, card: BCTADD, allow_overwrites: bool=False) -> None:
        """adds an BCTADD object"""
        key = card.csid
        self.bctadds[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bctpara_object(self, card: BCTPARA, allow_overwrites: bool=False) -> None:
        """adds an BCTPARA object"""
        key = card.csid
        self.bctparas[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bctset_object(self, card: BCTSET, allow_overwrites: bool=False) -> None:
        """adds an BCTSET object"""
        key = card.csid
        self.bctsets[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bconp_object(self, bconp: BCONP) -> None:
        """adds an BCONP object"""
        key = bconp.contact_id
        self.bconp[key] = bconp
        self._type_to_id_map[bconp.type].append(key)

    def _add_blseg_object(self, blseg: BLSEG) -> None:
        """adds an BLSEG object"""
        key = blseg.line_id
        self.blseg[key] = blseg
        self._type_to_id_map[blseg.type].append(key)

    def _add_bfric_object(self, bfric: BFRIC) -> None:
        """adds an BFRIC object"""
        key = bfric.friction_id
        self.bfric[key] = bfric
        self._type_to_id_map[bfric.type].append(key)

    def _add_bsurf_object(self, card: BSURF, allow_overwrites: bool=False) -> None:
        """adds an BSURF object"""
        key = card.sid
        self.bsurf[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_bsurfs_object(self, card: BSURFS, allow_overwrites: bool=False) -> None:
        """adds an BSURFS object"""
        key = card.id
        self.bsurfs[key] = card
        self._type_to_id_map[card.type].append(key)

    def _add_radcav_object(self, radcav: RADCAV, allow_overwrites: bool=False) -> None:
        """adds an RADCAV object"""
        key = radcav.icavity
        if key in self.radcavs and not allow_overwrites:
            if not radcav == self.radcavs[key]:
                assert key not in self.radcavs, 'pid=%s old RADCAV=\n%snew RADCAV=\n%s' % (key, self.radcavs[key], radcav)
        else:
            assert key > 0, 'pid=%s radcav=%s' % (key, radcav)
            self.radcavs[key] = radcav
            self._type_to_id_map[radcav.type].append(key)

    def _add_radmtx_object(self, radmtx: RADMTX, allow_overwrites: bool=False) -> None:
        """adds an RADMTX object"""
        key = radmtx.icavity
        if key in self.radmtx and not allow_overwrites:
            if not radmtx == self.radmtx[key]:
                assert key not in self.radmtx, 'pid=%s old RADMTX=\n%snew RADMTX=\n%s' % (key, self.radmtx[key], radmtx)
        else:
            assert key > 0, 'pid=%s radmtx=%s' % (key, radmtx)
            self.radmtx[key] = radmtx
            self._type_to_id_map[radmtx.type].append(key)

    def _add_tempd_object(self, tempd: TEMPD, allow_overwrites: bool=False) -> None:
        """adds an TEMPD object"""
        key = tempd.sid
        if key in self.tempds and not allow_overwrites:
            if not tempd == self.tempds[key]:
                assert key not in self.tempds, 'TEMPD.sid=%s old=\n%snew=\n%s' % (
                    key, self.tempds[key], tempd)
        else:
            assert key > 0, 'sid=%s tempd=%s' % (key, tempd)
            self.tempds[key] = tempd
            self._type_to_id_map[tempd.type].append(key)

    def _add_pbusht_object(self, prop: PBUSHT, allow_overwrites: bool=False) -> None:
        """adds an PBUSHT object"""
        key = prop.pid
        if key in self.pbusht and not allow_overwrites:
            if not prop == self.pbusht[key]:
                assert key not in self.pbusht, 'PBUSHT.pid=%s old=\n%snew=\n%s' % (
                    key, self.pbusht[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pbusht[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_pdampt_object(self, prop: PDAMPT, allow_overwrites: bool=False) -> None:
        """adds an PDAMPT object"""
        key = prop.pid
        if key in self.pdampt and not allow_overwrites:
            if not prop == self.pdampt[key]:
                assert key not in self.pdampt, 'PDAMPT.pid=%s old=\n%snew=\n%s' % (
                    key, self.pdampt[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            self.pdampt[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_pelast_object(self, prop: PELAST, allow_overwrites: bool=False) -> None:
        """adds an PELAST object"""
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        if key in self.pelast and not allow_overwrites:
            if not prop == self.pelast[key]:
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' % (key, self.pelast[key],prop))
                assert key not in self.pelast, 'PELAST.pid=%s old=\n%snew=\n%s' % (
                    key, self.pelast[key], prop)
        else:
            self.pelast[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_tf_object(self, tf: TF, allow_overwrites: bool=False) -> None:
        """adds an TF (transfer function) object"""
        key = tf.sid
        assert key > 0, 'sid=%s tf=%s' % (key, tf)
        if key in self.transfer_functions:
            self.transfer_functions[key].append(tf)
        else:
            self.transfer_functions[key] = [tf]
            self._type_to_id_map[tf.type].append(key)

    def _add_structural_material_object(self, material: Union[MAT1, MAT2, MAT3, MAT8, MAT9,
                                                              MAT10, MAT11, MAT3D, MATG],
                                        allow_overwrites: bool=False) -> None:
        """adds an MAT1, MAT2, MAT8 object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.materials and not allow_overwrites:
            if not material == self.materials[key]:
                self._duplicate_materials.append(material)
        else:
            self.materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_thermal_material_object(self, material: Union[MAT4, MAT5],
                                     allow_overwrites: bool=False) -> None:
        """adds an MAT4, MAT5 object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.thermal_materials and not allow_overwrites:
            if not material == self.thermal_materials[key]:
                self._duplicate_thermal_materials.append(material)
        else:
            self.thermal_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_hyperelastic_material_object(self, material: Union[MATHE, MATHP],
                                          allow_overwrites: bool=False) -> None:
        """adds an MATHP, MATHE object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.hyperelastic_materials and not allow_overwrites:
            if not material == self.hyperelastic_materials[key]:
                assert key not in self.hyperelastic_materials, 'mid=%s\nold=\n%snew=\n%s' % (key, self.hyperelastic_materials[key], material)
        else:
            self.hyperelastic_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_material_dependence_object(self, material: Union[MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9,
                                                              MATS1], # MATS3, MATS8
                                        allow_overwrites: bool=False) -> None:
        """
        adds the following objects:
            MATS1, MATS3, MATS8,
            MATT1, MATT2, MATT3,
            MATT4, MATT5, MATT8, MATT9
        """
        Type = material.type
        key = material.mid
        mapper = {
            'MATS1' : self.MATS1,
            'MATS3' : self.MATS3,
            'MATS8' : self.MATS8,

            'MATT1' : self.MATT1,
            'MATT2' : self.MATT2,
            'MATT3' : self.MATT3,
            'MATT4' : self.MATT4,
            'MATT5' : self.MATT5,
            'MATT8' : self.MATT8,
            'MATT9' : self.MATT9,
        }
        slot = mapper[Type]
        if key in slot and not allow_overwrites:
            if not material == slot[key]:
                assert key not in slot, 'dMATx.mid=%s Type=%r\nold=\n%snew=\n%s' % (key, Type, slot[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            slot[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_creep_material_object(self, material: CREEP, allow_overwrites: bool=False) -> None:
        """
        Adds a CREEP material

        Notes
        -----
        May be removed in the future.  Are CREEP cards materials?
        They have an MID, but reference structural materials.

        """
        key = material.mid
        if key in self.thermal_materials and not allow_overwrites:
            if not material == self.creep_materials[key]:
                assert key not in self.creep_materials, 'Material.mid=%s\nold=\n%snew=\n%s' % (key, self.creep_materials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.creep_materials[key] = material
            self._type_to_id_map[material.type].append(key)

    def _add_coord_object(self, coord: Union[CORD1R, CORD1C, CORD1S,
                                             CORD2R, CORD2C, CORD2S], # CORD3G
                          allow_overwrites: bool=False) -> None:
        """adds a CORDx object"""
        key = coord.cid
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        if key in self.coords:
            #if not allow_overwrites:
            if not coord == self.coords[key]:
                self._duplicate_coords.append(coord)
        else:
            self.coords[key] = coord
            self._type_to_id_map[coord.type].append(key)

    def _add_load_combination_object(self, load: Union[LOAD, CLOAD]) -> None:
        """adds a load object to a load case"""
        key = load.sid
        if key in self.load_combinations:
            self.load_combinations[key].append(load)
        else:
            self.load_combinations[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_load_object(self, load: Union[FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                           PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                           GRAV, ACCEL, ACCEL1, SPCD, SLOAD,
                                           QBDY1, QBDY2, QBDY3, QVOL, TEMPAX, PRESAX,
                                           RFORCE, RFORCE1, LOADCYN, LOADCYH, DEFORM,
                                           GMLOAD]) -> None:
        """adds a load object to a load case"""
        key = load.sid
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_dload_object(self, load: DLOAD) -> None:
        """adds a dload object to a load case"""
        key = load.sid
        if key in self.dloads:
            self.dloads[key].append(load)
        else:
            self.dloads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_dload_entry(self, dload: Union[ACSRCE, RANDPS, RANDT1,
                                            TLOAD1, TLOAD2, RLOAD1, RLOAD2,
                                            QVECT]) -> None:
        """adds a sub-dload object to a load case"""
        key = dload.sid
        if key in self.dload_entries:
            self.dload_entries[key].append(dload)
        else:
            self.dload_entries[key] = [dload]
            self._type_to_id_map[dload.type].append(key)

    def _add_lseq_object(self, load: LSEQ) -> None:
        """adds a LSEQ object to a load case"""
        key = load.sid
        if key in self.load_combinations:
            self.load_combinations[key].append(load)
        else:
            self.load_combinations[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_thermal_load_object(self, load: Union[TEMP, TEMPB3, QHBDY, QBDY1, QBDY2, QBDY3]) -> None:
        # same function at the moment...
        key = load.sid
        assert key > 0, 'key=%s; load=%s\n' % (key, load)
        if key in self.loads:
            self.loads[key].append(load)
        else:
            self.loads[key] = [load]
            self._type_to_id_map[load.type].append(key)

    def _add_phbdy_object(self, prop: PHBDY) -> None:
        key = prop.pid
        if key in self.phbdys:
            if not prop == self.phbdys[key]:
                assert key not in self.phbdys, 'PHBDY.pid=%s\nold=\n%snew=\n%s' % (
                    key, self.phbdys[key], prop)
        else:
            assert key > 0, 'pid=%s prop=\n%s' % (key, prop)
            self.phbdys[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def _add_view_object(self, view: VIEW) -> None:
        """adds a VIEW object"""
        key = view.iview
        assert key > 0, 'key=%s; view=%s\n' % (key, view)
        if key in self.views:
            if not view == self.views[key]:
                assert key not in self.views, 'VIEW.iview=%s\nold=\n%snew=\n%s' % (
                    key, self.views[key], view)
        else:
            assert key > 0, 'iview=%s view=\n%s' % (key, view)
            self.views[key] = view
            self._type_to_id_map[view.type].append(key)

    def _add_view3d_object(self, view3d: VIEW3D) -> None:
        """adds a VIEW3D object"""
        key = view3d.icavity
        assert key > 0, 'key=%s; view3d=%s\n' % (key, view3d)
        if key in self.view3ds:
            if not view3d == self.view3ds[key]:
                assert key not in self.view3ds, 'VIEW3D.icavity=%s\nold=\n%snew=\n%s' % (
                    key, self.view3ds[key], view3d)
        else:
            assert key > 0, 'icavity=%s view3d=\n%s' % (key, view3d)
            self.view3ds[key] = view3d
            self._type_to_id_map[view3d.type].append(key)

    def _add_normal_object(self, snorm: SNORM) -> None:
        """adds an SNORM object"""
        key = snorm.nid
        assert key > 0, 'key=%s; snorm=%s\n' % (key, snorm)
        if key in self.normals:
            if not snorm == self.normals[key]:
                assert key not in self.normals, 'VIEW.iview=%s\nold=\n%snew=\n%s' % (
                    key, self.normals[key], snorm)
        else:
            assert key > 0, 'pid=%s SNORM=\n%s' % (key, snorm)
            self.normals[key] = snorm
            self._type_to_id_map[snorm.type].append(key)

    def _add_convection_property_object(self, prop: Union[PCONV, PCONVM]) -> None:
        key = prop.pconid
        assert key > 0, key
        assert key not in self.convection_properties, key
        self.convection_properties[key] = prop
        self._type_to_id_map[prop.type].append(key)

    def _add_thermal_bc_object(self, bc: Union[CONV, CONVM, RADM, TEMPBC], key) -> None:
        assert key > 0
        if key in self.bcs:
            self.bcs[key].append(bc)
        else:
            self.bcs[key] = [bc]
            self._type_to_id_map[bc.type].append(key)

    def _add_constraint_mpc_object(self, constraint: MPC) -> None: # MPCAX
        key = constraint.conid
        if key in self.mpcs:
            self.mpcs[key].append(constraint)
        else:
            self.mpcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_constraint_mpcadd_object(self, constraint: MPCADD) -> None:
        key = constraint.conid
        if key in self.mpcadds:
            self.mpcadds[key].append(constraint)
        else:
            self.mpcadds[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spc_object(self, constraint: Union[SPC, SPC1, SPCAX, GMSPC]) -> None:
        key = constraint.conid
        if key in self.spcs:
            self.spcs[key].append(constraint)
        else:
            self.spcs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spcadd_object(self, constraint: SPCADD) -> None:
        key = constraint.conid
        if key in self.spcadds:
            self.spcadds[key].append(constraint)
        else:
            self.spcadds[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spcoff_object(self, constraint: Union[SPCOFF, SPCOFF1]) -> None:
        """dumb key, but good enough..."""
        key = constraint.type
        if key in self.spcoffs:
            self.spcoffs[key].append(constraint)
        else:
            self.spcoffs[key] = [constraint]
            self._type_to_id_map[constraint.type].append(key)

    def _add_sesuport_object(self, se_suport: Union[SESUP, SESUPORT]) -> None:
        """adds an SESUPORT"""
        self._type_to_id_map[se_suport.type].append(len(self.se_suport))
        self.se_suport.append(se_suport)

    def _add_suport_object(self, suport: SUPORT) -> None:
        """adds a SUPORT"""
        self._type_to_id_map[suport.type].append(len(self.suport))
        self.suport.append(suport)

    def _add_suport1_object(self, suport1: SUPORT1) -> None:
        """adds a SUPORT1"""
        key = suport1.conid
        if key in self.suport1:
            self.suport1[key].add_suport1_to_set(suport1)
        else:
            assert suport1.conid > 0
            self.suport1[key] = suport1
            self._type_to_id_map[suport1.type].append(key)

    def _add_tic_object(self, tic: TIC, allow_overwrites: bool=False) -> None:
        """adds a TIC object"""
        key = tic.sid
        if key in self.tics:
            self.tics[key].add(tic)
        else:
            assert tic.sid > 0
            self.tics[key] = tic
            self._type_to_id_map[tic.type].append(key)

    def _add_darea_object(self, darea: DAREA, allow_overwrites: bool=False) -> None:
        """adds a DAREA object"""
        #key = (darea.sid, darea.p, darea.c)
        key = darea.sid
        if key in self.dareas:
            self.dareas[key].add(darea)
        else:
            assert darea.sid > 0
            self.dareas[key] = darea
            self._type_to_id_map[darea.type].append(key)

    def _add_dphase_object(self, dphase: DPHASE, allow_overwrites: bool=False) -> None:
        """adds a DPHASE object"""
        #key = (dphase.sid, dphase.nid, dphase.component) # dphase.phase_lead,
        key = dphase.sid
        if key in self.dphases:
            self.dphases[key].add(dphase)
        else:
            assert dphase.sid > 0, key
            self.dphases[key] = dphase
            self._type_to_id_map[dphase.type].append(key)

    def _add_delay_object(self, delay: DELAY, allow_overwrites: bool=False) -> None:
        """adds an DELAY object"""
        #key = (delay.sid, delay.nid, delay.component)
        key = delay.sid
        assert key > 0, 'sid=%s delay=%s' % (key, delay)
        if key in self.delays:
            self.delays[key].add(delay)
        else:
            self.delays[key] = delay
            self._type_to_id_map[delay.type].append(key)

    def _add_aero_object(self, aero: AERO) -> None:
        """adds an AERO object"""
        # only one AERO card allowed
        assert self.aero is None, '\naero=\n%s old=\n%s' % (aero, self.aero)
        self.aero = aero
        #self._type_to_id_map[aero.type].append(key)

    def _add_aeros_object(self, aeros: AEROS) -> None:
        """adds an AEROS object"""
        # only one AEROS card allowed
        assert self.aeros is None, '\naeros=\n%s old=\n%s' % (aeros, self.aeros)
        self.aeros = aeros
        #self._type_to_id_map[aeros.type].append(key)

    #def _add_aeroz_object(self, aeroz: AEROZ) -> None:
        #"""adds an AEROZ object"""
        #key = aeroz.sid
        #if key in self.aeroz and not allow_overwrites:
            #if not aeroz == self.zona.aeroz[key]:
                #assert key not in self.aeroz, 'AEROZ.sid=%s\nold=\n%snew=\n%s' % (key, self.aeroz[key], aeroz)
        #else:
            #assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            #self.aeroz[key] = aeroz
            #self._type_to_id_map[aeroz.type].append(key)

    def _add_baror_object(self, baror: BAROR) -> None:
        """adds an BAROR object"""
        # only one BAROR card allowed
        assert self.baror is None, '\nBAROR=\n%s old=\n%s' % (baror, self.baror)
        if self.baror is None:
            self.baror = baror

    def _add_beamor_object(self, beamor: BEAMOR) -> None:
        """adds an BEAMOR object"""
        # only one BAROR card allowed
        assert self.beamor is None, '\nBEAMOR=\n%s old=\n%s' % (beamor, self.beamor)
        if self.beamor is None:
            self.beamor = beamor

    def _add_axic_object(self, axic: AXIC) -> None:
        """adds an AXIC object"""
        # only one AXIC card allowed
        assert self.axic is None, '\naxic=\n%s old=\n%s' % (axic, self.axic)
        self.axic = axic

    def _add_axif_object(self, axif: AXIF) -> None:
        """adds an AXIF object"""
        # only one AXIC card allowed
        assert self.axif is None, '\naxif=\n%s old=\n%s' % (axif, self.axif)
        self.axif = axif

    def _add_acmodl_object(self, acmodl) -> None:
        """adds a ACMODL object"""
        assert self.acmodl is None, self.acmodl
        self.acmodl = acmodl

    def _add_cyax_object(self, cyax: CYAX) -> None:
        """adds an CYAX object"""
        # only one CYAX card allowed
        assert self.cyax is None, '\ncyax=\n%s old=\n%s' % (cyax, self.cyax)
        self.cyax = cyax

    def _add_cyjoin_object(self, cyjoin: CYJOIN) -> None:
        """adds an CYJOIN object"""
        key = cyjoin.side
        assert key not in self.cyjoin, 'CYJOIN.side=%s\nold=\n%snew=\n%s' % (key, self.cyjoin[key], cyjoin)
        assert key >= 0
        self.cyjoin[key] = cyjoin
        self._type_to_id_map[cyjoin.type].append(key)

    def _add_modtrak_object(self, modtrak: MODTRAK) -> None:
        """adds an MODTRAK object"""
        # only one CYAX card allowed
        assert self.modtrak is None, '\nmodtrak=\n%s old=\n%s' % (modtrak, self.modtrak)
        self.modtrak = modtrak

    def _add_aefact_object(self, aefact: AEFACT, allow_overwrites: bool=False) -> None:
        """adds an AEFACT object"""
        key = aefact.sid
        if key in self.aefacts and not allow_overwrites:
            if not aefact == self.aefacts[key]:
                assert key not in self.aefacts, 'AEFACT.sid=%s\nold=\n%snew=\n%s' % (key, self.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.aefacts[key] = aefact
            self._type_to_id_map[aefact.type].append(key)

    def _add_aelist_object(self, aelist: AELIST) -> None:
        """adds an AELIST object"""
        key = aelist.sid
        assert key not in self.aelists, 'AELIST.sid=%s\nold=\n%snew=\n%s' % (key, self.aelists[key], aelist)
        assert key >= 0
        self.aelists[key] = aelist
        self._type_to_id_map[aelist.type].append(key)

    def _add_aelink_object(self, aelink: AELINK) -> None:
        """adds an AELINK object"""
        key = aelink.aelink_id
        assert key >= 0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        self._type_to_id_map[aelink.type].append(key)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.aestats[key])

    def _add_aecomp_object(self, aecomp: Union[AECOMP, AECOMPL]) -> None:
        """adds an AECOMP object"""
        key = aecomp.name
        assert key not in self.aecomps, '\naecomp=\n%s oldAECOMP=\n%s' % (aecomp, self.aecomps[key])
        self.aecomps[key] = aecomp
        self._type_to_id_map[aecomp.type].append(key)

    def _add_aeparm_object(self, aeparam: AEPARM) -> None:
        """adds an AEPARM object"""
        key = aeparam.aeparm_id
        assert key not in self.aeparams, '\naeparam=\n%s oldAEPARM=\n%s' % (aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam
        self._type_to_id_map[aeparam.type].append(key)

    def _add_aestat_object(self, aestat: AESTAT) -> None:
        """adds an AESTAT object"""
        key = aestat.aestat_id
        assert key not in self.aestats, '\naestat=\n%s old=\n%s' % (
            aestat, self.aestats[key])
        assert key >= 0
        self.aestats[key] = aestat
        self._type_to_id_map[aestat.type].append(key)

    def _add_aesurf_object(self, aesurf: AESURF) -> None:
        """adds an AESURF object"""
        key = aesurf.aesid
        assert key not in self.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, self.aesurf[key])
        assert key >= 0
        self.aesurf[key] = aesurf
        self._type_to_id_map[aesurf.type].append(key)

    def _add_aesurfs_object(self, aesurfs: AESURFS) -> None:
        """adds an AESURFS object"""
        key = aesurfs.aesid
        assert key not in self.aesurfs, '\naesurfs=\n%s old=\n%s' % (
            aesurfs, self.aesurfs[key])
        assert key >= 0
        self.aesurfs[key] = aesurfs
        self._type_to_id_map[aesurfs.type].append(key)

    def _add_csschd_object(self, csschd: CSSCHD) -> None:
        """adds an CSSCHD object"""
        key = csschd.sid
        assert key not in self.csschds, '\nCSSCHD=\n%s old=\n%s' % (csschd, self.csschds[key])
        assert key >= 0
        self.csschds[key] = csschd
        self._type_to_id_map[csschd.type].append(key)

    def _add_caero_object(self, caero: Union[CAERO1, CAERO2, CAERO3, CAERO4, CAERO5]) -> None:
        """adds an CAERO1/CAERO2/CAERO3/CAERO4/CAERO5 object"""
        key = caero.eid
        assert key not in self.caeros, '\nkey=%s; caero=\n%r old_caero=\n%r' % (
            key, caero, self.caeros[key])
        assert key > 0
        self.caeros[key] = caero
        self._type_to_id_map[caero.type].append(key)

    def _add_paero_object(self, paero: Union[PAERO1, PAERO2, PAERO3, PAERO4, PAERO5]) -> None:
        """adds an PAERO1/PAERO2/PAERO3/PAERO4/PAERO5 object"""
        key = paero.pid
        assert key not in self.paeros, '\npaero=\n%r old_paero=\n%r' % (
            paero, self.paeros[key])
        assert key > 0, 'paero.pid = %r' % (key)
        self.paeros[key] = paero
        self._type_to_id_map[paero.type].append(key)

    def _add_monpnt_object(self, monitor_point: Union[MONPNT1, MONPNT2, MONPNT3]) -> None:
        """adds an MONPNT object"""
        key = monitor_point.name
        assert key not in self.monitor_points, '\nmonitor_point=\n%soldMOTPNT=\n%s' % (
            monitor_point, self.monitor_points[key])
        self.monitor_points.append(monitor_point)
        self._type_to_id_map[monitor_point.type].append(len(self.monitor_points) - 1)

    def _add_spline_object(self, spline: Union[SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5]) -> None:
        """adds an SPLINE1/SPLINE2/SPLINE3/SPLINE4/SPLINE5 object"""
        key = spline.eid
        assert spline.eid not in self.splines, f'\nspline={spline}\n%sold_spline=\n{self.splines[key]}'
        assert spline.eid > 0, spline
        self.splines[key] = spline
        self._type_to_id_map[spline.type].append(key)

    def _add_gust_object(self, gust: GUST) -> None:
        """adds an GUST object"""
        key = gust.sid
        assert key not in self.gusts
        assert key > 0
        self.gusts[key] = gust
        self._type_to_id_map[gust.type].append(key)

    def _add_trim_object(self, trim: Union[TRIM, TRIM2], allow_overwrites: bool=False) -> None:
        """adds an TRIM object"""
        key = trim.sid
        if not allow_overwrites:
            assert key not in self.trims, 'TRIM=%s  old=\n%snew=\n%s' % (key, self.trims[key], trim)
        assert key > 0, 'key=%r trim=\n%s' % (key, trim)
        self.trims[key] = trim
        self._type_to_id_map[trim.type].append(key)

    def _add_diverg_object(self, diverg: DIVERG, allow_overwrites: bool=False) -> None:
        """adds an DIVERG object"""
        key = diverg.sid
        if not allow_overwrites:
            assert key not in self.divergs, 'DIVERG=%s  old=\n%snew=\n%s' % (key, self.divergs[key], diverg)
        assert key > 0, 'key=%r diverg=\n%s' % (key, diverg)
        self.divergs[key] = diverg
        self._type_to_id_map[diverg.type].append(key)

    def _add_flutter_object(self, flutter: FLUTTER) -> None:
        """adds an FLUTTER object"""
        key = flutter.sid
        assert key not in self.flutters, 'FLUTTER=%s old=\n%snew=\n%s' % (key, self.flutters[key], flutter)
        assert key > 0
        self.flutters[key] = flutter
        self._type_to_id_map[flutter.type].append(key)

    def _add_flfact_object(self, flfact: FLFACT) -> None:
        """adds an FLFACT object"""
        key = flfact.sid
        #assert key not in self.flfacts
        assert key > 0
        self.flfacts[key] = flfact  # set id...
        self._type_to_id_map[flfact.type].append(key)

    def _add_dconstr_object(self, dconstr: [DCONSTR, DCONADD]) -> None:
        """adds a DCONSTR object"""
        #key = (dconstr.oid, dconstr.rid)
        key = dconstr.oid
        #assert key not in self.dconstrs, 'key=%r DCONSTR/DCONADD=\n%s' % (key, dconstr)
        assert dconstr.oid > 0
        #assert dconstr.dresp_id > 0
        if key in self.dconstrs:
            self.dconstrs[key].append(dconstr)
        else:
            self.dconstrs[key] = [dconstr]
        self._type_to_id_map[dconstr.type].append(key)

    #def add_DCONADD(self, dconadd, allow_overwrites: bool=False) -> None:
        #key = dconadd.oid
        #if key in self.dconstrs and not allow_overwrites:
            #if not dconadd == self.dconstrs[key]:
                #assert key not in self.dconstrs, 'DCONADD=%s old=\n%snew=\n%s' % (
                    #key, self.dconstrs[key], dconadd)
        #else:
            #assert key > 0, 'dcid=%s dconadd=%s' % (key, dconadd)
            #self.dconstrs[key] = dconadd
            #self._type_to_id_map[dconadd.type].append(key)

    def _add_desvar_object(self, desvar: DESVAR) -> None:
        """adds a DESVAR object"""
        key = desvar.desvar_id
        assert key not in self.desvars, 'DESVAR=%s old=\n%snew=\n%s' % (
            key, self.desvars[key], desvar)
        assert key > 0
        self.desvars[key] = desvar
        self._type_to_id_map[desvar.type].append(key)

    def _add_topvar_object(self, topvar: TOPVAR) -> None:
        """adds a TOPVAR object"""
        key = topvar.opt_id
        assert key not in self.topvar, 'TOPVAR=%s old=\n%snew=\n%s' % (
            key, self.topvar[key], topvar)
        assert key > 0
        self.topvar[key] = topvar
        self._type_to_id_map[topvar.type].append(key)

    def _add_ddval_object(self, ddval: DDVAL) -> None:
        """adds a DDVAL object"""
        key = ddval.oid
        assert key not in self.ddvals, 'DDVAL=%s old=\n%snew=\n%s' % (
            key, self.ddvals[key], ddval)
        assert key > 0
        self.ddvals[key] = ddval
        self._type_to_id_map[ddval.type].append(key)

    def _add_dlink_object(self, dlink: DLINK) -> None:
        """adds a DLINK object"""
        key = dlink.oid
        assert key not in self.dlinks, 'DLINK=%s old=\n%snew=\n%s' % (
            key, self.dlinks[key], dlink)
        assert key > 0
        self.dlinks[key] = dlink
        self._type_to_id_map[dlink.type].append(key)

    def _add_dscreen_object(self, dscreen: DSCREEN) -> None:
        """adds a DSCREEN object"""
        key = dscreen.rtype
        assert key not in self.dscreen, 'DSCREEN=%s old=\n%snew=\n%s' % (
            key, self.dscreen[key], dscreen)
        assert len(key) > 0, 'key=%r' % key
        self.dscreen[key] = dscreen
        self._type_to_id_map[dscreen.type].append(key)

    def _add_dresp_object(self, dresp: Union[DRESP1, DRESP2, DRESP3]) -> None:
        """adds a DRESP1/DRESP2/DRESP3 object"""
        key = dresp.dresp_id
        assert key not in self.dresps, 'DRESPx=%s old=\n%snew=\n%s' % (
            key, self.dresps[key], dresp)
        assert key > 0
        self.dresps[key] = dresp
        self._type_to_id_map[dresp.type].append(key)

    def _add_dvcrel_object(self, dvcrel: Union[DVCREL1, DVCREL2]) -> None:
        """adds a DVCREL1/DVCREL2 object"""
        key = dvcrel.oid
        assert key not in self.dvcrels, 'DVCRELx=%s old\n%snew=\n%s' % (
            key, self.dvcrels[key], dvcrel)
        assert key > 0
        self.dvcrels[key] = dvcrel
        self._type_to_id_map[dvcrel.type].append(key)

    def _add_dvmrel_object(self, dvmrel: Union[DVMREL1, DVMREL2]) -> None:
        """adds a DVMREL1/DVMREL2 object"""
        key = dvmrel.oid
        assert key not in self.dvmrels, 'DVMRELx=%s old=\n%snew=\n%s' % (
            key, self.dvmrels[key], dvmrel)
        assert key not in self.dvmrels
        assert key > 0
        self.dvmrels[key] = dvmrel
        self._type_to_id_map[dvmrel.type].append(key)

    def _add_dvprel_object(self, dvprel: Union[DVPREL1, DVPREL2]) -> None:
        """adds a DVPREL1/DVPREL2 object"""
        key = dvprel.oid
        assert key not in self.dvprels, 'DVPRELx=%s old\n%snew=\n%s' % (
            key, self.dvprels[key], dvprel)
        assert key > 0
        self.dvprels[key] = dvprel
        self._type_to_id_map[dvprel.type].append(key)

    def _add_dvgrid_object(self, dvgrid: DVGRID) -> None:
        """adds a DVGRID object"""
        key = dvgrid.dvid
        assert key > 0
        if key not in self.dvgrids:
            self.dvgrids[key] = []
            self._type_to_id_map[dvgrid.type].append(key)
        self.dvgrids[key].append(dvgrid)

    def _add_nlparm_object(self, nlparm: NLPARM) -> None:
        """adds a NLPARM object"""
        key = nlparm.nlparm_id
        assert key not in self.nlparms
        assert key > 0, 'key=%s; nlparm=%s\n' % (key, nlparm)
        self.nlparms[key] = nlparm
        self._type_to_id_map[nlparm.type].append(key)

    def _add_rotor_object(self, rotor: Union[ROTORD, ROTORG]) -> None:
        """adds a ROTORD/ROTORG object"""
        key = rotor.sid
        assert key > 0, 'key=%s; rotor=%s\n' % (key, rotor)
        if key in self.rotors:
            rotor_old = self.rotors[key]
            assert rotor.type == rotor_old.type
            self.rotors[key].nids += rotor.nids
        else:
            self.rotors[key] = rotor
        self._type_to_id_map[rotor.type].append(key)

    def _add_nlpci_object(self, nlpci: NLPCI) -> None:
        """adds a NLPCI object"""
        key = nlpci.nlpci_id
        assert key not in self.nlpcis
        assert key > 0
        self.nlpcis[key] = nlpci
        self._type_to_id_map[nlpci.type].append(key)

    def _add_nxstrat_object(self, nxstrat: NXSTRAT) -> None:
        key = nxstrat.sid
        assert key not in self.nxstrats, 'nxstrats=%s nxstrat=%s' % (self.nxstrats, nxstrat)
        assert key > 0
        self.nxstrats[key] = nxstrat
        self._type_to_id_map[nxstrat.type].append(key)

    def _add_tstep_object(self, tstep: Union[TSTEP, TSTEP1],
                          allow_overwrites: bool=False) -> None:
        """adds a TSTEP object"""
        key = tstep.sid
        if key in self.tsteps and not allow_overwrites:
            if not tstep == self.tsteps[key]:
                assert key not in self.tsteps, 'TSTEP=%s\nold=\n%snew=\n%s' % (key, self.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.tsteps[key] = tstep
            self._type_to_id_map[tstep.type].append(key)

    def _add_tstepnl_object(self, tstepnl: TSTEPNL,
                            allow_overwrites: bool=False) -> None:
        """adds a TSTEPNL object"""
        key = tstepnl.sid
        if key in self.tstepnls and not allow_overwrites:
            if not tstepnl == self.tstepnls[key]:
                assert key not in self.tstepnls, 'TSTEPNL=%s\nold=\n%snew=\n%s' % (key, self.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.tstepnls[key] = tstepnl
            self._type_to_id_map[tstepnl.type].append(key)

    def _add_freq_object(self, freq: Union[FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5]) -> None:
        key = freq.sid
        assert key > 0
        if key in self.frequencies:
            freq0 = self.frequencies[key][0]
            if freq0.type == 'FREQ' and freq.type == 'FREQ':
                freq0.add_frequency_object(freq)
            else:
                self.frequencies[key].append(freq)
        else:
            self.frequencies[key] = [freq]
            self._type_to_id_map[freq.type].append(key)

    def _add_set_object(self, set_obj: Union[SET1, SET2, SET3]) -> None:
        """adds an SET1/SET3 object"""
        key = set_obj.sid
        assert key >= 0
        if key in self.sets:
            self.sets[key].add_set(set_obj)
        else:
            self.sets[key] = set_obj
            self._type_to_id_map[set_obj.type].append(key)

    def _add_radset_object(self, set_obj: RADSET) -> None:
        """adds an RADSET object"""
        if self.radset:
            self.radset.add_set(set_obj)
        else:
            self.radset = set_obj
            #self._type_to_id_map[set_obj.type].append(key)

    def _add_aset_object(self, set_obj: Union[ASET, ASET1]) -> None:
        """adds an ASET/ASET1 object"""
        self.asets.append(set_obj)
        n = len(self._type_to_id_map['ASET'])
        self._type_to_id_map['ASET'].append(n)

    def _add_omit_object(self, set_obj: Union[OMIT, OMIT1]) -> None:
        """adds an OMIT/OMIT1 object"""
        self.omits.append(set_obj)
        n = len(self._type_to_id_map['OMIT'])
        self._type_to_id_map['OMIT'].append(n)

    def _add_bset_object(self, set_obj: Union[BSET, BSET1]) -> None:
        """adds an BSET/BSET1 object"""
        self.bsets.append(set_obj)
        n = len(self._type_to_id_map['BSET'])
        self._type_to_id_map['BSET'].append(n)

    def _add_cset_object(self, set_obj: Union[CSET, CSET1]) -> None:
        """adds an CSET/USET1 object"""
        self.csets.append(set_obj)
        n = len(self._type_to_id_map['CSET'])
        self._type_to_id_map['CSET'].append(n)

    def _add_qset_object(self, set_obj: Union[QSET, QSET1]) -> None:
        """adds an QSET/QSET1 object"""
        self.qsets.append(set_obj)
        n = len(self._type_to_id_map['QSET'])
        self._type_to_id_map['QSET'].append(n)

    def _add_uset_object(self, set_obj: Union[USET, USET1]) -> None:
        """adds an USET/USET1 object"""
        key = set_obj.name
        if key in self.usets:
            self.usets[key].append(set_obj)
        else:
            self.usets[key] = [set_obj]
        self._type_to_id_map[set_obj.type].append(key)

    def _add_sebset_object(self, set_obj: Union[SEBSET, SEBSET1]) -> None:
        """adds an SEBSET/SEBSET1 object"""
        self.se_bsets.append(set_obj)

    def _add_secset_object(self, set_obj: Union[SECSET, SECSET1]) -> None:
        """adds an SECSET/SECSTE1 object"""
        self.se_csets.append(set_obj)

    def _add_seqset_object(self, set_obj: Union[SEQSET, SEQSET1]) -> None:
        """adds an SEQSET/SEQSET1 object"""
        self.se_qsets.append(set_obj)

    def _add_seuset_object(self, set_obj: Union[SEUSET, SEUSET1]) -> None:
        """adds an SEUSET/SEUSET1 object"""
        key = set_obj.name
        if key in self.se_usets:
            self.se_usets[key].append(set_obj)
        else:
            self.se_usets[key] = [set_obj]
        self._type_to_id_map[set_obj.type].append(key)

    def _add_seset_object(self, set_obj: SESET) -> None:
        """adds an SESET object"""
        key = set_obj.seid
        assert key >= 0
        if key in self.se_sets:
            old_set = self.se_sets[key]
            set_obj.add_seset(old_set)
        self.se_sets[key] = set_obj
        self._type_to_id_map[set_obj.type].append(key)

    def _add_table_object(self, table: Union[TABLEH1, TABLEHT, TABLES1, TABLEST]) -> None:
        """adds a TABLES1, TABLEST object"""
        key = table.tid
        if key in self.tables:
            if not table == self.tables[key]:
                assert key not in self.tables, '\ntable=\n%s old_table=\n%s' % (
                    table, self.tables[key])
        assert key > 0
        self.tables[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_tabled_object(self, table: Union[TABLED1, TABLED2, TABLED3, TABLED4]) -> None:
        """adds a TABLED1, TABLED2, TABLED3, TABLED4 object"""
        key = table.tid
        assert key not in self.tables_d, '\ntabled=\n%s old_tabled=\n%s' % (
            table, self.tables_d[key])
        #assert key > 0; yes you can have negative tables...
        self.tables_d[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_tablem_object(self, table: Union[TABLEM1, TABLEM2, TABLEM3, TABLEM4]) -> None:
        """adds a TABLEM1, TABLEM2, TABLEM3, TABLEM4 object"""
        key = table.tid
        assert key not in self.tables_m, '\ntablem=\n%s old_tablem=\n%s' % (
            table, self.tables_m[key])
        #assert key > 0; yes you can have negative tables...
        self.tables_m[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_table_sdamping_object(self, table: TABDMP1) -> None:
        """adds a TABDMP1 object"""
        key = table.tid
        assert key not in self.tables_sdamping, '\nTable=\n%s oldTable=\n%s' % (
            table, self.tables_sdamping[key])
        #assert key > 0; yes you can have negative tables...
        self.tables_sdamping[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_random_table_object(self, table: Union[TABRND1, TABRNDG]) -> None:
        """adds a TABRND1, TABRNDG object"""
        key = table.tid
        assert key not in self.random_tables, '\nTable=\n%s old=\n%s' % (
            table, self.random_tables[key])
        assert key > 0
        self.random_tables[key] = table
        self._type_to_id_map[table.type].append(key)

    def _add_method_object(self, method: Union[EIGR, EIGRL, EIGB],
                           allow_overwrites: bool=False) -> None:
        """adds a EIGR/EIGRL object"""
        key = method.sid
        if key in self.methods and not allow_overwrites:
            if not method == self.methods[key]:
                assert key not in self.methods, 'sid=%s\nold_method=\n%snew_method=\n%s' % (key, self.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.methods[key] = method
            self._type_to_id_map[method.type].append(key)

    def _add_cmethod_object(self, method: Union[EIGC, EIGP],
                            allow_overwrites: bool=False) -> None:
        """adds a EIGB/EIGC object"""
        key = method.sid
        if key in self.cMethods and not allow_overwrites:
            if not method == self.cMethods[key]:
                assert key not in self.cMethods, 'sid=%s\nold_cmethod=\n%snew_cmethod=\n%s' % (key, self.cMethods[key], method)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, method)
            self.cMethods[key] = method
            self._type_to_id_map[method.type].append(key)

    def _add_mkaero_object(self, mkaero: Union[MKAERO1, MKAERO2]) -> None:
        """adds an MKAERO1/MKAERO2 object"""
        self.mkaeros.append(mkaero)

    #---------------------------------------------------------------------------
    # parametric
    def _add_pset(self, pset: PSET, allow_overwrites: bool=False) -> None:
        assert pset.idi not in self.pset, pset
        self.pset[pset.idi] = pset

    def _add_pval(self, pval: PVAL, allow_overwrites: bool=False) -> None:
        if pval.idi not in self.pval:
            self.pval[pval.idi] = []
        self.pval[pval.idi].append(pval)

    def _add_gmcurv(self, curve: GMCURV, allow_overwrites: bool=False) -> None:
        assert curve.curve_id not in self.gmcurv, curve
        self.gmcurv[curve.curve_id] = curve

    def _add_gmsurf(self, surf: GMSURF, allow_overwrites: bool=False) -> None:
        assert surf.surf_id not in self.gmsurf, surf
        self.gmsurf[surf.surf_id] = surf

    def _add_feface(self, face: FEFACE, allow_overwrites: bool=False) -> None:
        key = face.face_id
        if key in self.feface and not allow_overwrites:
            if not face == self.feface[key]:
                raise RuntimeError(f'feface is duplicated\n{face}\nold:\n{self.feface[key]}')
        else:
            self.feface[face.face_id] = face
            self._type_to_id_map[face.type].append(key)
        #assert face.face_id not in self.feface, face
        #self.feface[face.face_id] = face

    def _add_feedge(self, edge: FEEDGE, allow_overwrites: bool=False) -> None:
        key = edge.edge_id
        if key in self.feedge and not allow_overwrites:
            if not edge == self.feedge[key]:
                raise RuntimeError(f'feedge is duplicated\n{edge}\nold:\n{self.feedge[key]}')
        else:
            self.feedge[edge.edge_id] = edge
            self._type_to_id_map[edge.type].append(key)

    #---------------------------------------------------------------------------
