# pylint: disable=R0902,R0904,R0914
from __future__ import annotations
from typing import Union, TYPE_CHECKING
from pyNastran.bdf.cards.nodes import SPOINT, EPOINT
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import (
        BDF,
        CYAX, CYJOIN, AXIF,
        TOPVAR, MPCAX, CORD3G,
        SESUPORT, SEUSET, SEUSET1,
    )
    from pyNastran.bdf.cards.bolt import BOLT, BOLTFOR, BOLTSEQ, BOLTLD
    from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL, GENEL
    from pyNastran.bdf.cards.properties.properties import PFAST, PGAP, PRAC2D, PRAC3D
    from pyNastran.bdf.cards.properties.solid import PLSOLID, PSOLID, PIHEX, PCOMPS, PCOMPLS
    #from pyNastran.bdf.cards.msgmesh import CGEN, GMCORD

    from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
    from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

    from pyNastran.bdf.cards.elements.solid import (
        #CTETRA, CPYRAM, CPENTA, CHEXA,
        CIHEX1, CIHEX2, CHEXA1, CHEXA2,
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
                                                        MATCID,
                                                        )
    from pyNastran.bdf.cards.deqatn import DEQATN
    from pyNastran.bdf.cards.dynamic import (
        DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
        TSTEP, TSTEP1, TSTEPNL, NLPARM, NLPCI, TF, ROTORG, ROTORD, TIC)
    from pyNastran.bdf.cards.loads.loads import (
        LSEQ, SLOAD, DAREA, RFORCE, RFORCE1, SPCD, DEFORM,
        LOADCYN, LOADCYH)
    from pyNastran.bdf.cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
    from pyNastran.bdf.cards.loads.static_loads import (LOAD, CLOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                        FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                        PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                                        GMLOAD)
    from pyNastran.bdf.cards.loads.random_loads import RANDPS, RANDT1

    from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                               MAT8, MAT9, MAT10, MAT11, MAT3D,
                                               MATG, MATHE, MATHP, CREEP, EQUIV,
                                               MATDMG,
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
        DVTREL1,
        DVGRID, DSCREEN)
    from pyNastran.bdf.cards.optimization_nx import (
        GROUP, DMNCON, DMNCON,
    )
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
    from pyNastran.bdf.cards.params import PARAM, MDLPRM
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
        BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCPARA, BCTPARA, BCONP, BLSEG,
        BFRIC, BCTPARM, BGADD, BGSET, BCBODY)
    from pyNastran.bdf.cards.parametric.geometry import PSET, PVAL, FEEDGE, FEFACE, GMCURV, GMSURF
    from pyNastran.bdf.cards.elements.acoustic import PACABS

class AddMethods:
    """defines methods to add card objects to the BDF"""
    def __init__(self, model: BDF) -> None:
        #BDFAttributes.__init__(self)
        self.model = model

    #@property
    #def log(self) -> SimpleLogger:
        #return self.model.log

    #@property
    #def _type_to_id_map(self) -> dict[str, Any]:
        #return self.model._type_to_id_map

    def _add_dmi_object(self, dmi: DMI, allow_overwrites: bool=False) -> None:
        """adds a DMI object"""
        name = dmi.name
        self.model.dmi[name] = dmi
        self.model._type_to_id_map[dmi.type].append(name)

    def _add_dmig_object(self, dmig: DMIG, allow_overwrites: bool=False) -> None:
        """adds a DMIG object"""
        name = dmig.name
        self.model.dmig[name] = dmig
        self.model._type_to_id_map[dmig.type].append(name)

    def _add_dmiax_object(self, dmiax: DMIAX, allow_overwrites: bool=False) -> None:
        """adds a DMI object"""
        name = dmiax.name
        self.model.dmiax[name] = dmiax
        self.model._type_to_id_map[dmiax.type].append(name)

    def _add_dmij_object(self, dmij: DMIJ, allow_overwrites: bool=False) -> None:
        """adds a DMIJ object"""
        name = dmij.name
        self.model.dmij[name] = dmij
        self.model._type_to_id_map[dmij.type].append(name)

    def _add_dmiji_object(self, dmiji: DMIJI, allow_overwrites: bool=False) -> None:
        """adds a DMIJI object"""
        name = dmiji.name
        self.model.dmiji[name] = dmiji
        self.model._type_to_id_map[dmiji.type].append(name)

    def _add_dmik_object(self, dmik: DMIK, allow_overwrites: bool=False) -> None:
        """adds a DMIK object"""
        name = dmik.name
        self.model.dmik[name] = dmik
        self.model._type_to_id_map[dmik.type].append(name)

    def _add_dti_object(self, dti: DTI, allow_overwrites: bool=False) -> None:
        """adds an DTI object"""
        name = dti.name
        model = self.model
        if name == 'UNITS' or name not in model.dti:
            model.dti[name] = dti
            model._type_to_id_map[dti.type].append(name)
        else:
            old_dti = model.dti[name]
            key = list(dti.fields.keys())[0]
            assert key not in old_dti.fields, 'key=%i old_fields=%s fields=%s' % (key, old_dti.fields, dti.fields)
            old_dti.fields[key] = dti.fields[key]

    def _add_param_object(self, param: PARAM, allow_overwrites: bool=False) -> None:
        """adds a PARAM object"""
        key = param.key
        model = self.model
        if key in model.params and not allow_overwrites:
            if not param == model.params[key]:
                #if param.key in self.params:
                    #msg = 'key=%s param=%s old_param=%s' % (key, param, self.params[key])
                    #raise KeyError(msg)
                model.log.warning('key=%s param=%s old_param=%s' %
                                  (key, param, model.params[key]))
                model.params[key] = param
        else:
            model.params[key] = param
            model._type_to_id_map[param.type].append(key)

    def _add_mdlprm_object(self, mdlprm: MDLPRM, allow_overwrites: bool=False) -> None:
        """adds a MDLPRM object"""
        if self.model.mdlprm is None:
            self.model.mdlprm = mdlprm
        else:
            model_mdlprm_dict = self.model.mdlprm.mdlprm_dict
            for key, value in mdlprm.mdlprm_dict.items():
                if key in model_mdlprm_dict:
                    assert self.model.mdlprm is None, self.model.mdlprm
                else:
                    model_mdlprm_dict[key] = value
        #model._type_to_id_map[param.type].append(key)

    def _add_node_object(self, node: GRID, allow_overwrites: bool=False) -> None:
        """adds a GRID card"""
        key = node.nid
        model = self.model
        if key in model.nodes and not allow_overwrites:
            if not node == model.nodes[key]:
                assert node.nid not in model.nodes, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, model.nodes[key], node)
            else:
                #print('GRID was duplicated...nid=%s; node=\n%s' % (key, node))
                pass
        else:
            assert key > 0, 'nid=%s node=%s' % (key, node)
            model.nodes[key] = node
            model._type_to_id_map[node.type].append(key)

    def _add_gridb_object(self, node: GRIDB, allow_overwrites: bool=False) -> None:
        """adds a GRIDB card"""
        key = node.nid
        model = self.model
        assert key > 0, 'eid=%s node=%s' % (key, node)
        if key in self.model.gridb and not allow_overwrites:
            assert node.nid not in model.gridb, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, model.gridb[key], node)
        model.gridb[key] = node
        model._type_to_id_map[node.type].append(key)
        model._is_axis_symmetric = True

    def _add_ringfl_object(self, ringfl: RINGFL, allow_overwrites: bool=False) -> None:
        """adds a RINGFL card"""
        key = ringfl.ringfl
        assert key > 0, 'eid=%s ringfl=%s' % (key, ringfl)
        if key in self.model.ringfl and not allow_overwrites:
            assert ringfl.ringfl not in self.model.ringfl, 'ringfl=%s\nold_ringfl=\n%snew_ringfl=\n%s' % (ringfl.ringfl, self.model.ringfl[key], ringfl)
        self.model.ringfl[key] = ringfl
        self.model._type_to_id_map[ringfl.type].append(key)
        self.model._is_axis_symmetric = True

    def _add_ringax_object(self, ringax: Union[RINGAX, POINTAX], allow_overwrites: bool=False) -> None:
        """adds a RINGAX card"""
        key = ringax.nid
        model = self.model
        if key in self.model.ringaxs and not allow_overwrites:
            if not ringax == model.ringaxs[key]:
                assert ringax.nid not in model.ringaxs, 'nid=%s\nold_ringax=\n%snew_ringax=\n%s' % (ringax.nid, model.ringaxs[key], ringax)
            else:
                #print('RINGAX was duplicated...nid=%s; ringax=\n%s' % (key, ringax))
                pass
        else:
            assert key > 0, 'nid=%s ringax=%s' % (key, ringax)
            model.ringaxs[key] = ringax
            model._type_to_id_map[ringax.type].append(key)
        model._is_axis_symmetric = True

    def _add_seqgp_object(self, seqgp: SEQGP) -> None:
        """adds an SEQGP card"""
        if self.model.seqgp is None:
            self.model.seqgp = seqgp
        else:
            self.model.seqgp.append(seqgp)

    def _add_point_object(self, point: POINT, allow_overwrites: bool=False) -> None:
        """adds a POINT card"""
        key = point.nid
        model = self.model
        if key in model.points and not allow_overwrites:
            if not point == model.points[key]:
                assert point.nid not in model.points, 'nid=%s\nold_point=\n%snew_point=\n%s' % (point.nid, model.points[key], point)
            else:
                #print('POINT was duplicated...nid=%s; point=\n%s' % (key, point))
                pass
        else:
            assert key > 0, 'nid=%s point=%s' % (key, point)
            model.points[key] = point
            model._type_to_id_map[point.type].append(key)

    def _add_spoint_object(self, spoints: SPOINTs) -> None:
        """adds an SPOINT card"""
        comment = spoints.comment
        if hasattr(spoints, 'ifile'):
            ifile = spoints.ifile
            for nid in spoints.points:
                if nid in self.model.spoints:
                    continue
                spoint = SPOINT(nid, comment=comment)
                spoint.ifile = ifile
                comment = ''
                self.model.spoints[nid] = spoint
        else:
            for nid in spoints.points:
                if nid in self.model.spoints:
                    continue
                spoint = SPOINT(nid, comment=comment)
                comment = ''
                self.model.spoints[nid] = spoint

    def _add_epoint_object(self, epoints: EPOINTs) -> None:
        """adds an EPOINT card"""
        comment = epoints.comment
        for nid in epoints.points:
            if nid in self.model.epoints:
                continue
            epoint = EPOINT(nid, comment=comment)
            comment = ''
            self.model.epoints[nid] = epoint

    def _add_setree_object(self, setree: SETREE) -> None:
        key = setree.seid
        self.model.setree[key] = setree
        self.model._type_to_id_map[setree.type].append(key)
    def _add_senqset_object(self, senqset: SENQSET) -> None:
        key = senqset.set_id
        self.model.senqset[key] = senqset
        self.model._type_to_id_map[senqset.type].append(key)
    def _add_sebulk_object(self, sebulk: SEBULK) -> None:
        key = sebulk.seid
        self.model.sebulk[key] = sebulk
        self.model._type_to_id_map[sebulk.type].append(key)
    def _add_release_object(self, release: RELEASE) -> None:
        key = release.seid
        self.model.release[key] = release
        self.model._type_to_id_map[release.type].append(key)
    def _add_sebndry_object(self, sebndry: SEBNDRY) -> None:
        key = (sebndry.seid_a, sebndry.seid_b)
        self.model.sebndry[key] = sebndry
    def _add_seloc_object(self, seloc: SELOC) -> None:
        key = seloc.seid
        self.model.seloc[key] = seloc
        self.model._type_to_id_map[seloc.type].append(key)
    def _add_sempln_object(self, sempln: SEMPLN) -> None:
        key = sempln.seid
        self.model.sempln[key] = sempln
        self.model._type_to_id_map[sempln.type].append(key)

    def _add_seconct_object(self, seconct: SECONCT) -> None:
        key = (seconct.seid_a, seconct.seid_b)
        self.model.seconct[key] = seconct
        self.model._type_to_id_map[seconct.type].append(key)
    def _add_selabel_object(self, selabel: SELABEL) -> None:
        key = selabel.seid
        self.model.selabel[key] = selabel
        self.model._type_to_id_map[selabel.type].append(key)
    def _add_seexcld_object(self, seexcld: SEEXCLD) -> None:
        key = (seexcld.seid_a, seexcld.seid_b)
        self.model.seexcld[key] = seexcld
        self.model._type_to_id_map[seexcld.type].append(key)

    def _add_seelt_object(self, seelt: SEELT) -> None:
        #self.model.seelt.append(seelt)
        key = seelt.seid
        self.model.seelt[key] = seelt
        self.model._type_to_id_map[seelt.type].append(key)
    def _add_seload_object(self, seload: SELOAD) -> None:
        key = seload.seid
        self.model.seload[key] = seload
        self.model._type_to_id_map[seload.type].append(key)
    def _add_csuper_object(self, csuper: CSUPER) -> None:
        key = csuper.seid
        self.model.csuper[key] = csuper
        self.model._type_to_id_map[csuper.type].append(key)
    def _add_csupext_object(self, csupext: CSUPEXT) -> None:
        key = csupext.seid
        self.model.csupext[key] = csupext
        self.model._type_to_id_map[csupext.type].append(key)

    def _add_plotel_object(self, elem: PLOTEL, allow_overwrites: bool=False) -> None:
        """adds an PLOTEL object"""
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if not allow_overwrites:
            if key in self.model.elements:
                if elem == self.model.elements[key]:
                    self.model._duplicate_elements.append(elem)
                    if self.model._stop_on_duplicate_error:
                        self.model.pop_parse_errors()
            elif key in self.model.plotels:
                if not elem == self.model.plotels[key]:
                    assert elem.eid not in self.model.plotels, 'eid=%s\nold_element=\n%snew_element=\n%s' % (elem.eid, self.model.plotels[elem.eid], elem)
        self.model.plotels[key] = elem
        self.model._type_to_id_map[elem.type].append(key)

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
                                              CIHEX1, CIHEX2,
                                              CHEXA1, CHEXA2,
                                              ],
                            allow_overwrites: bool=False) -> None:
        key = elem.eid
        model = self.model
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        if key in model.elements and not allow_overwrites:
            if not elem == model.elements[key]:
                model._duplicate_elements.append(elem)
                if model._stop_on_duplicate_error:
                    model.pop_parse_errors()
        else:
            model.elements[key] = elem
            model._type_to_id_map[elem.type].append(key)

    def _add_ao_object(self, elem_flag: CBARAO, allow_overwrites: bool=False) -> None:
        """adds a CBARAO"""
        key = elem_flag.eid
        model = self.model
        assert key > 0, 'eid=%s must be positive; elem_flag=\n%s' % (key, elem_flag)
        if key in model.ao_element_flags and not allow_overwrites:
            if not elem_flag == model.ao_element_flags[key]:
                #self.model._duplicate_elements.append(elem_flag)
                #if self.model._stop_on_duplicate_error:
                    #self.model.pop_parse_errors()
                assert elem_flag.eid not in self.model.ao_element_flags, 'eid=%s\nold_ao_element=\n%snew_ao_element=\n%s' % (
                    elem_flag.eid, model.ao_element_flags[elem_flag.eid], elem_flag)
        else:
            model.ao_element_flags[key] = elem_flag
            model._type_to_id_map[elem_flag.type].append(key)

    def _add_doptprm_object(self, doptprm: DOPTPRM) -> None:
        """adds a DOPTPRM"""
        self.model.doptprm = doptprm

    def _add_nsm_object(self, nsm: Union[NSM, NSM1, NSML, NSML1], allow_overwrites: bool=False) -> None:
        """adds a nsm object to a nsm set"""
        key = nsm.sid
        assert key > 0, 'sid=%s must be positive; nsm=\n%s' % (key, nsm)
        if key in self.model.nsms:
            self.model.nsms[key].append(nsm)
        else:
            self.model.nsms[key] = [nsm]
            self.model._type_to_id_map[nsm.type].append(key)

    def _add_nsmadd_object(self, nsmadd: NSMADD, allow_overwrites: bool=False) -> None:
        """adds a nsmadd object to a nsm set"""
        key = nsmadd.sid
        assert key > 0, 'sid=%s must be positive; nsmadd=\n%s' % (key, nsmadd)
        if key in self.model.nsmadds:
            self.model.nsmadds[key].append(nsmadd)
        else:
            self.model.nsmadds[key] = [nsmadd]
            self.model._type_to_id_map[nsmadd.type].append(key)

    def _add_mass_object(self, mass: Union[CMASS1, CMASS2, CMASS3, CMASS4,
                                           CONM1, CONM2], allow_overwrites: bool=False) -> None:
        key = mass.eid
        model = self.model
        assert key > 0, 'eid=%s must be positive; mass=\n%s' % (key, mass)
        if key in model.masses and not allow_overwrites:
            if not mass == self.model.masses[key]:
                model._duplicate_masses.append(mass)
        else:
            model.masses[key] = mass
            model._type_to_id_map[mass.type].append(key)

    def _add_damper_object(self, elem, allow_overwrites: bool=False) -> None:
        """.. warning:: can dampers have the same ID as a standard element?"""
        return self._add_element_object(elem, allow_overwrites)

    def _add_rigid_element_object(self, elem: Union[RBAR, RBAR1,
                                                    RBE1, RBE2, RBE3,
                                                    RROD, RSPLINE, RSSCON],
                                  allow_overwrites: bool=False) -> None:
        key = elem.eid
        model = self.model
        assert key > 0, 'eid=%s elem=%s' % (key, elem)
        if key in model.rigid_elements and not allow_overwrites:
            assert elem.eid not in model.rigid_elements, 'eid=%s\noldElement=\n%snewElement=\n%s' % (elem.eid, model.rigid_elements[elem.eid], elem)
        model.rigid_elements[key] = elem
        model._type_to_id_map[elem.type].append(key)

    def _add_thermal_element_object(self, elem: Union[CHBDYE, CHBDYG, CHBDYP]) -> None:
        """same as add_element at the moment..."""
        self._add_element_object(elem)

    def _add_deqatn_object(self, deqatn: DEQATN, allow_overwrites: bool=False) -> None:
        """adds an DEQATN object"""
        key = deqatn.equation_id
        assert key > 0, 'ID=%s deqatn\n%s' % (key, deqatn)
        model = self.model
        if key in model.dequations and not allow_overwrites:
            if not deqatn.write_card() == model.dequations[key].write_card():
                assert key not in model.dequations, 'id=%s old_eq=\n%snew_eq=\n%s' % (
                    key, model.dequations[key], deqatn)
        model.dequations[key] = deqatn
        model._type_to_id_map[deqatn.type].append(key)

    def _add_acoustic_property_object(self, prop: PACABS) -> None:
        self._add_property_object(prop)

    def _add_property_object(self, prop: Union[PELAS, PBUSH, PBUSH1D, PDAMP, PDAMP5, # PBUSH2D,
                                               PFAST, PVISC, PGAP, PRAC2D, PRAC3D, # PWELD
                                               PROD, PTUBE,
                                               PBAR, PBARL, PBRSECT, PCONEAX,
                                               PBEAM, PBEAML, PBCOMP, PBMSECT,
                                               PBEND, PBEAM3, PPLANE, PSHEAR,
                                               PSHELL, PCOMP, PCOMPG, PLPLANE,
                                               PSOLID, PLSOLID, PIHEX, PCOMPS, PCOMPLS],
                             allow_overwrites: bool=False) -> None:
        """
        adds one of the following objects:
          PELAS, PBUSH, PBUSH1D, PBUSH2D, PDAMP,
          PROD, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
          PSHELL, PCOMP, PCOMPG,
          PSOLID, PLSOLID, PIHEX, PCOMPS, PCOMPLS
        """
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        model = self.model
        if key in model.properties and not allow_overwrites:
            if not prop == model.properties[key]:
                model._duplicate_properties.append(prop)
                if model._stop_on_duplicate_error:
                    model.pop_parse_errors()
        else:
            model.properties[key] = prop
            model._type_to_id_map[prop.type].append(key)

    def _add_property_mass_object(self, prop: PMASS, allow_overwrites: bool=False) -> None:
        """adds an PMASS object"""
        key = prop.pid
        model = self.model
        if key in model.properties_mass and not allow_overwrites:
            if not prop == self.model.properties_mass[key]:
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' %(key, model.properties_mass[key],prop))
                assert key not in model.properties_mass, 'pid=%s oldProperty=\n%snewProperty=\n%s' % (key, model.properties_mass[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            model.properties_mass[key] = prop
            model._type_to_id_map[prop.type].append(key)

    def _add_dtable_object(self, dtable: DTABLE, allow_overwrites: bool=False) -> None:
        """adds an DTABLE object"""
        model = self.model
        if model.dtable is not None:
            if not dtable == model.dtable:
                raise RuntimeError('DTABLE cannot be overwritten\nold:\n%s\nnew:\n%s',
                                   model.dtable, dtable)
        else:
            self.model.dtable = dtable
            #self.model._type_to_id_map[dtable.type].append(1)

    def _add_bcrpara_object(self, card: BCRPARA, allow_overwrites: bool=False) -> None:
        """adds an BCRPARA object"""
        key = card.crid
        self.model.bcrparas[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bgadd_object(self, card: BGADD, allow_overwrites: bool=False) -> None:
        """adds an BGADD object"""
        key = card.glue_id
        self.model.bgadds[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bctadd_object(self, card: BCTADD, allow_overwrites: bool=False) -> None:
        """adds an BCTADD object"""
        key = card.csid
        self.model.bctadds[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bcpara_object(self, card: BCPARA, allow_overwrites: bool=False) -> None:
        """adds an BCPARA object"""
        key = card.csid
        self.model.bcparas[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bctpara_object(self, card: BCTPARA, allow_overwrites: bool=False) -> None:
        """adds an BCTPARA object"""
        key = card.csid
        self.model.bctparas[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bctparam_object(self, card: BCTPARM, allow_overwrites: bool=False) -> None:
        """adds an BCTPARM object"""
        key = card.csid
        self.model.bctparms[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bctset_object(self, card: BCTSET, allow_overwrites: bool=False) -> None:
        """adds an BCTSET object"""
        key = card.csid
        self.model.bctsets[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bgset_object(self, card: BGSET, allow_overwrites: bool=False) -> None:
        """adds an BGSET object"""
        key = card.glue_id
        self.model.bgsets[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bconp_object(self, bconp: BCONP) -> None:
        """adds an BCONP object"""
        key = bconp.contact_id
        self.model.bconp[key] = bconp
        self.model._type_to_id_map[bconp.type].append(key)

    def _add_bcbody_object(self, bcbody: BCBODY) -> None:
        """adds an BCBODY object"""
        key = bcbody.contact_id
        self.model.bcbodys[key] = bcbody
        self.model._type_to_id_map[bcbody.type].append(key)

    def _add_blseg_object(self, blseg: BLSEG) -> None:
        """adds an BLSEG object"""
        key = blseg.line_id
        self.model.blseg[key] = blseg
        self.model._type_to_id_map[blseg.type].append(key)

    def _add_bfric_object(self, bfric: BFRIC) -> None:
        """adds an BFRIC object"""
        key = bfric.friction_id
        self.model.bfric[key] = bfric
        self.model._type_to_id_map[bfric.type].append(key)

    def _add_bsurf_object(self, card: BSURF, allow_overwrites: bool=False) -> None:
        """adds an BSURF object"""
        key = card.sid
        self.model.bsurf[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_bsurfs_object(self, card: BSURFS, allow_overwrites: bool=False) -> None:
        """adds an BSURFS object"""
        key = card.id
        self.model.bsurfs[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def _add_radcav_object(self, radcav: RADCAV, allow_overwrites: bool=False) -> None:
        """adds an RADCAV object"""
        key = radcav.icavity
        model = self.model
        if key in self.model.radcavs and not allow_overwrites:
            if not radcav == model.radcavs[key]:
                assert key not in model.radcavs, 'pid=%s old RADCAV=\n%snew RADCAV=\n%s' % (key, model.radcavs[key], radcav)
        else:
            assert key > 0, 'pid=%s radcav=%s' % (key, radcav)
            model.radcavs[key] = radcav
            model._type_to_id_map[radcav.type].append(key)

    def _add_radmtx_object(self, radmtx: RADMTX, allow_overwrites: bool=False) -> None:
        """adds an RADMTX object"""
        key = radmtx.icavity
        model = self.model
        if key in model.radmtx and not allow_overwrites:
            if not radmtx == model.radmtx[key]:
                assert key not in model.radmtx, 'pid=%s old RADMTX=\n%snew RADMTX=\n%s' % (key, model.radmtx[key], radmtx)
        else:
            assert key > 0, 'pid=%s radmtx=%s' % (key, radmtx)
            model.radmtx[key] = radmtx
            model._type_to_id_map[radmtx.type].append(key)

    def _add_tempd_object(self, tempd: TEMPD, allow_overwrites: bool=False) -> None:
        """adds an TEMPD object"""
        key = tempd.sid
        model = self.model
        if key in model.tempds and not allow_overwrites:
            if not tempd == model.tempds[key]:
                assert key not in model.tempds, 'TEMPD.sid=%s old=\n%snew=\n%s' % (
                    key, model.tempds[key], tempd)
        else:
            assert key > 0, 'sid=%s tempd=%s' % (key, tempd)
            model.tempds[key] = tempd
            model._type_to_id_map[tempd.type].append(key)

    def _add_pbusht_object(self, prop: PBUSHT, allow_overwrites: bool=False) -> None:
        """adds an PBUSHT object"""
        key = prop.pid
        model = self.model
        if key in model.pbusht and not allow_overwrites:
            if not prop == model.pbusht[key]:
                assert key not in model.pbusht, 'PBUSHT.pid=%s old=\n%snew=\n%s' % (
                    key, model.pbusht[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            model.pbusht[key] = prop
            model._type_to_id_map[prop.type].append(key)

    def _add_pdampt_object(self, prop: PDAMPT, allow_overwrites: bool=False) -> None:
        """adds an PDAMPT object"""
        key = prop.pid
        model = self.model
        if key in model.pdampt and not allow_overwrites:
            if not prop == model.pdampt[key]:
                assert key not in model.pdampt, 'PDAMPT.pid=%s old=\n%snew=\n%s' % (
                    key, model.pdampt[key], prop)
        else:
            assert key > 0, 'pid=%s prop=%s' % (key, prop)
            model.pdampt[key] = prop
            model._type_to_id_map[prop.type].append(key)

    def _add_pelast_object(self, prop: PELAST, allow_overwrites: bool=False) -> None:
        """adds an PELAST object"""
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        model = self.model
        if key in model.pelast and not allow_overwrites:
            if not prop == model.pelast[key]:
                #print('pid=%s\noldProperty=\n%snewProperty=\n%s' % (key, model.pelast[key],prop))
                assert key not in model.pelast, 'PELAST.pid=%s old=\n%snew=\n%s' % (
                    key, model.pelast[key], prop)
        else:
            model.pelast[key] = prop
            model._type_to_id_map[prop.type].append(key)

    def _add_tf_object(self, tf: TF, allow_overwrites: bool=False) -> None:
        """adds an TF (transfer function) object"""
        key = tf.sid
        assert key > 0, 'sid=%s tf=%s' % (key, tf)
        model = self.model
        if key in model.transfer_functions:
            model.transfer_functions[key].append(tf)
        else:
            model.transfer_functions[key] = [tf]
            model._type_to_id_map[tf.type].append(key)

    def _add_structural_material_object(self, material: Union[MAT1, MAT2, MAT3, MAT8, MAT9,
                                                              MAT10, MAT11, MAT3D, MATG],
                                        allow_overwrites: bool=False) -> None:
        """adds an MAT1, MAT2, MAT8 object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        model = self.model
        if key in model.materials and not allow_overwrites:
            if not material == model.materials[key]:
                model._duplicate_materials.append(material)
        else:
            model.materials[key] = material
            model._type_to_id_map[material.type].append(key)

    def _add_thermal_material_object(self, material: Union[MAT4, MAT5],
                                     allow_overwrites: bool=False) -> None:
        """adds an MAT4, MAT5 object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        model = self.model
        if key in model.thermal_materials and not allow_overwrites:
            if not material == model.thermal_materials[key]:
                model._duplicate_thermal_materials.append(material)
        else:
            model.thermal_materials[key] = material
            model._type_to_id_map[material.type].append(key)

    def _add_hyperelastic_material_object(self, material: Union[MATHE, MATHP],
                                          allow_overwrites: bool=False) -> None:
        """adds an MATHP, MATHE object"""
        key = material.mid
        assert key > 0, 'mid=%s material=\n%s' % (key, material)
        if key in self.model.hyperelastic_materials and not allow_overwrites:
            if not material == self.model.hyperelastic_materials[key]:
                assert key not in self.model.hyperelastic_materials, 'mid=%s\nold=\n%snew=\n%s' % (key, self.model.hyperelastic_materials[key], material)
        else:
            self.model.hyperelastic_materials[key] = material
            self.model._type_to_id_map[material.type].append(key)

    def _add_material_dependence_object(self, material: Union[MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9,
                                                              MATS1, MATDMG], # MATS3, MATS8
                                        allow_overwrites: bool=False) -> None:
        """
        adds the following objects:
            MATS1, MATS3, MATS8,
            MATT1, MATT2, MATT3,
            MATT4, MATT5, MATT8, MATT9,
            MATDMG
        """
        Type = material.type
        key = material.mid
        mapper = {
            'MATS1' : self.model.MATS1,
            'MATS3' : self.model.MATS3,
            'MATS8' : self.model.MATS8,

            'MATT1' : self.model.MATT1,
            'MATT2' : self.model.MATT2,
            'MATT3' : self.model.MATT3,
            'MATT4' : self.model.MATT4,
            'MATT5' : self.model.MATT5,
            'MATT8' : self.model.MATT8,
            'MATT9' : self.model.MATT9,
            'MATDMG': self.model.MATDMG,
        }
        slot = mapper[Type]
        if key in slot and not allow_overwrites:
            if not material == slot[key]:
                assert key not in slot, 'dMATx.mid=%s Type=%r\nold=\n%snew=\n%s' % (key, Type, slot[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            slot[key] = material
            self.model._type_to_id_map[material.type].append(key)

    def _add_creep_material_object(self, material: CREEP, allow_overwrites: bool=False) -> None:
        """
        Adds a CREEP material

        Notes
        -----
        May be removed in the future.  Are CREEP cards materials?
        They have an MID, but reference structural materials.

        """
        key = material.mid
        if key in self.model.thermal_materials and not allow_overwrites:
            if not material == self.model.creep_materials[key]:
                assert key not in self.model.creep_materials, 'Material.mid=%s\nold=\n%snew=\n%s' % (key, self.model.creep_materials[key], material)
        else:
            assert key > 0, 'mid=%s material=\n%s' % (key, material)
            self.model.creep_materials[key] = material
            self.model._type_to_id_map[material.type].append(key)

    def _add_coord_object(self, coord: Union[CORD1R, CORD1C, CORD1S,
                                             CORD2R, CORD2C, CORD2S], # CORD3G
                          allow_overwrites: bool=False) -> None:
        """adds a CORDx object"""
        key = coord.cid
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        if key in self.model.coords:
            #if not allow_overwrites:
            if not coord == self.model.coords[key]:
                self.model._duplicate_coords.append(coord)
        else:
            self.model.coords[key] = coord
            self.model._type_to_id_map[coord.type].append(key)

    def _add_matcid_object(self, matcid: Union[MATCID]) -> None:
        """adds a MATCID object"""
        key = matcid.cid
        assert matcid.cid > -1, 'cid=%s coord=\n%s' % (key, matcid)

        # Multiple MATCIDs can share the same CID
        if key in self.model.MATCID:
            self.model.MATCID[key].append(matcid)
        else:
            self.model.MATCID[key] = [matcid]
            self.model._type_to_id_map[matcid.type].append(key)

    def _add_load_combination_object(self, load: Union[LOAD, CLOAD]) -> None:
        """adds a load object to a load case"""
        key = load.sid
        if key in self.model.load_combinations:
            self.model.load_combinations[key].append(load)
        else:
            self.model.load_combinations[key] = [load]
            self.model._type_to_id_map[load.type].append(key)

    def _add_load_object(self, load: Union[FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                           PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                           GRAV, ACCEL, ACCEL1, SPCD, SLOAD,
                                           QBDY1, QBDY2, QBDY3, QVOL, TEMPAX, PRESAX,
                                           RFORCE, RFORCE1, LOADCYN, LOADCYH, DEFORM,
                                           GMLOAD]) -> None:
        """adds a load object to a load case"""
        key = load.sid
        if key in self.model.loads:
            self.model.loads[key].append(load)
        else:
            self.model.loads[key] = [load]
            self.model._type_to_id_map[load.type].append(key)

    def _add_dload_object(self, load: DLOAD) -> None:
        """adds a dload object to a load case"""
        key = load.sid
        if key in self.model.dloads:
            self.model.dloads[key].append(load)
        else:
            self.model.dloads[key] = [load]
            self.model._type_to_id_map[load.type].append(key)

    def _add_dload_entry(self, dload: Union[ACSRCE, RANDPS, RANDT1,
                                            TLOAD1, TLOAD2, RLOAD1, RLOAD2,
                                            QVECT]) -> None:
        """adds a sub-dload object to a load case"""
        key = dload.sid
        if key in self.model.dload_entries:
            self.model.dload_entries[key].append(dload)
        else:
            self.model.dload_entries[key] = [dload]
            self.model._type_to_id_map[dload.type].append(key)

    def _add_lseq_object(self, load: LSEQ) -> None:
        """adds a LSEQ object to a load case"""
        key = load.sid
        if key in self.model.load_combinations:
            self.model.load_combinations[key].append(load)
        else:
            self.model.load_combinations[key] = [load]
            self.model._type_to_id_map[load.type].append(key)

    def _add_thermal_load_object(self, load: Union[TEMP, TEMPB3, QHBDY, QBDY1, QBDY2, QBDY3]) -> None:
        # same function at the moment...
        key = load.sid
        assert key > 0, 'key=%s; load=%s\n' % (key, load)
        if key in self.model.loads:
            self.model.loads[key].append(load)
        else:
            self.model.loads[key] = [load]
            self.model._type_to_id_map[load.type].append(key)

    def _add_phbdy_object(self, prop: PHBDY) -> None:
        key = prop.pid
        if key in self.model.phbdys:
            if not prop == self.model.phbdys[key]:
                assert key not in self.model.phbdys, 'PHBDY.pid=%s\nold=\n%snew=\n%s' % (
                    key, self.model.phbdys[key], prop)
        else:
            assert key > 0, 'pid=%s prop=\n%s' % (key, prop)
            self.model.phbdys[key] = prop
            self.model._type_to_id_map[prop.type].append(key)

    def _add_view_object(self, view: VIEW) -> None:
        """adds a VIEW object"""
        key = view.iview
        assert key > 0, 'key=%s; view=%s\n' % (key, view)
        if key in self.model.views:
            if not view == self.model.views[key]:
                assert key not in self.model.views, 'VIEW.iview=%s\nold=\n%snew=\n%s' % (
                    key, self.model.views[key], view)
        else:
            assert key > 0, 'iview=%s view=\n%s' % (key, view)
            self.model.views[key] = view
            self.model._type_to_id_map[view.type].append(key)

    def _add_view3d_object(self, view3d: VIEW3D) -> None:
        """adds a VIEW3D object"""
        key = view3d.icavity
        assert key > 0, 'key=%s; view3d=%s\n' % (key, view3d)
        if key in self.model.view3ds:
            if not view3d == self.model.view3ds[key]:
                assert key not in self.model.view3ds, 'VIEW3D.icavity=%s\nold=\n%snew=\n%s' % (
                    key, self.model.view3ds[key], view3d)
        else:
            assert key > 0, 'icavity=%s view3d=\n%s' % (key, view3d)
            self.model.view3ds[key] = view3d
            self.model._type_to_id_map[view3d.type].append(key)

    def _add_normal_object(self, snorm: SNORM) -> None:
        """adds an SNORM object"""
        key = snorm.nid
        assert key > 0, 'key=%s; snorm=%s\n' % (key, snorm)
        if key in self.model.normals:
            if not snorm == self.model.normals[key]:
                assert key not in self.model.normals, 'VIEW.iview=%s\nold=\n%snew=\n%s' % (
                    key, self.model.normals[key], snorm)
        else:
            assert key > 0, 'pid=%s SNORM=\n%s' % (key, snorm)
            self.model.normals[key] = snorm
            self.model._type_to_id_map[snorm.type].append(key)

    def _add_convection_property_object(self, prop: Union[PCONV, PCONVM]) -> None:
        key = prop.pconid
        if key in self.model.convection_properties:
            if not prop == self.model.convection_properties[key]:
                assert key not in self.model.convection_properties, 'PCONV/PCONVM.pconid=%s\nold=\n%snew=\n%s' % (
                    key, self.model.convection_properties[key], prop)
        else:
            assert key > 0, 'pid=%s PCONV/PCONVM=\n%s' % (key, prop)
            #assert key > 0, key
            #assert key not in self.model.convection_properties, key
            self.model.convection_properties[key] = prop
            self.model._type_to_id_map[prop.type].append(key)

    def _add_thermal_bc_object(self, bc: Union[CONV, CONVM, RADM, TEMPBC], key) -> None:
        assert key > 0
        if key in self.model.bcs:
            self.model.bcs[key].append(bc)
        else:
            self.model.bcs[key] = [bc]
            self.model._type_to_id_map[bc.type].append(key)

    def _add_constraint_mpc_object(self, constraint: MPC) -> None: # MPCAX
        key = constraint.conid
        if key in self.model.mpcs:
            self.model.mpcs[key].append(constraint)
        else:
            self.model.mpcs[key] = [constraint]
            self.model._type_to_id_map[constraint.type].append(key)

    def _add_constraint_mpcadd_object(self, constraint: MPCADD) -> None:
        key = constraint.conid
        if key in self.model.mpcadds:
            self.model.mpcadds[key].append(constraint)
        else:
            self.model.mpcadds[key] = [constraint]
            self.model._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spc_object(self, constraint: Union[SPC, SPC1, SPCAX, GMSPC]) -> None:
        key = constraint.conid
        if key in self.model.spcs:
            self.model.spcs[key].append(constraint)
        else:
            self.model.spcs[key] = [constraint]
            self.model._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spcadd_object(self, constraint: SPCADD) -> None:
        key = constraint.conid
        if key in self.model.spcadds:
            self.model.spcadds[key].append(constraint)
        else:
            self.model.spcadds[key] = [constraint]
            self.model._type_to_id_map[constraint.type].append(key)

    def _add_constraint_spcoff_object(self, constraint: Union[SPCOFF, SPCOFF1]) -> None:
        """dumb key, but good enough..."""
        key = constraint.type
        if key in self.model.spcoffs:
            self.model.spcoffs[key].append(constraint)
        else:
            self.model.spcoffs[key] = [constraint]
            self.model._type_to_id_map[constraint.type].append(key)

    def _add_sesuport_object(self, se_suport: Union[SESUP, SESUPORT]) -> None:
        """adds an SESUPORT"""
        self.model._type_to_id_map[se_suport.type].append(len(self.model.se_suport))
        self.model.se_suport.append(se_suport)

    def _add_suport_object(self, suport: SUPORT) -> None:
        """adds a SUPORT"""
        self.model._type_to_id_map[suport.type].append(len(self.model.suport))
        self.model.suport.append(suport)

    def _add_suport1_object(self, suport1: SUPORT1) -> None:
        """adds a SUPORT1"""
        key = suport1.conid
        if key in self.model.suport1:
            self.model.suport1[key].add_suport1_to_set(suport1)
        else:
            assert suport1.conid > 0
            self.model.suport1[key] = suport1
            self.model._type_to_id_map[suport1.type].append(key)

    def _add_tic_object(self, tic: TIC, allow_overwrites: bool=False) -> None:
        """adds a TIC object"""
        key = tic.sid
        if key in self.model.tics:
            self.model.tics[key].add(tic)
        else:
            assert tic.sid > 0
            self.model.tics[key] = tic
            self.model._type_to_id_map[tic.type].append(key)

    def _add_darea_object(self, darea: DAREA, allow_overwrites: bool=False) -> None:
        """adds a DAREA object"""
        #key = (darea.sid, darea.p, darea.c)
        key = darea.sid
        if key in self.model.dareas:
            self.model.dareas[key].add(darea)
        else:
            assert darea.sid > 0
            self.model.dareas[key] = darea
            self.model._type_to_id_map[darea.type].append(key)

    def _add_dphase_object(self, dphase: DPHASE, allow_overwrites: bool=False) -> None:
        """adds a DPHASE object"""
        #key = (dphase.sid, dphase.nid, dphase.component) # dphase.phase_lead,
        key = dphase.sid
        if key in self.model.dphases:
            self.model.dphases[key].add(dphase)
        else:
            assert dphase.sid > 0, key
            self.model.dphases[key] = dphase
            self.model._type_to_id_map[dphase.type].append(key)

    def _add_delay_object(self, delay: DELAY, allow_overwrites: bool=False) -> None:
        """adds an DELAY object"""
        #key = (delay.sid, delay.nid, delay.component)
        key = delay.sid
        assert key > 0, 'sid=%s delay=%s' % (key, delay)
        if key in self.model.delays:
            self.model.delays[key].add(delay)
        else:
            self.model.delays[key] = delay
            self.model._type_to_id_map[delay.type].append(key)

    def _add_aero_object(self, aero: AERO) -> None:
        """adds an AERO object"""
        # only one AERO card allowed
        assert self.model.aero is None, '\naero=\n%s old=\n%s' % (aero, self.model.aero)
        self.model.aero = aero
        #self.model._type_to_id_map[aero.type].append(key)

    def _add_aeros_object(self, aeros: AEROS) -> None:
        """adds an AEROS object"""
        # only one AEROS card allowed
        assert self.model.aeros is None, '\naeros=\n%s old=\n%s' % (aeros, self.model.aeros)
        self.model.aeros = aeros
        #self.model._type_to_id_map[aeros.type].append(key)

    #def _add_aeroz_object(self, aeroz: AEROZ) -> None:
        #"""adds an AEROZ object"""
        #key = aeroz.sid
        #if key in self.model.aeroz and not allow_overwrites:
            #if not aeroz == self.model.zona.aeroz[key]:
                #assert key not in self.model.aeroz, 'AEROZ.sid=%s\nold=\n%snew=\n%s' % (key, self.model.aeroz[key], aeroz)
        #else:
            #assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            #self.model.aeroz[key] = aeroz
            #self.model._type_to_id_map[aeroz.type].append(key)

    def _add_baror_object(self, baror: BAROR) -> None:
        """adds an BAROR object"""
        # only one BAROR card allowed
        assert self.model.baror is None, '\nBAROR=\n%s old=\n%s' % (baror, self.model.baror)
        if self.model.baror is None:
            self.model.baror = baror

    def _add_beamor_object(self, beamor: BEAMOR) -> None:
        """adds an BEAMOR object"""
        # only one BAROR card allowed
        assert self.model.beamor is None, '\nBEAMOR=\n%s old=\n%s' % (beamor, self.model.beamor)
        if self.model.beamor is None:
            self.model.beamor = beamor

    def _add_axic_object(self, axic: AXIC) -> None:
        """adds an AXIC object"""
        # only one AXIC card allowed
        assert self.model.axic is None, '\naxic=\n%s old=\n%s' % (axic, self.model.axic)
        self.model.axic = axic

    def _add_axif_object(self, axif: AXIF) -> None:
        """adds an AXIF object"""
        # only one AXIC card allowed
        assert self.model.axif is None, '\naxif=\n%s old=\n%s' % (axif, self.model.axif)
        self.model.axif = axif

    def _add_acmodl_object(self, acmodl) -> None:
        """adds a ACMODL object"""
        assert self.model.acmodl is None, self.model.acmodl
        self.model.acmodl = acmodl

    def _add_cyax_object(self, cyax: CYAX) -> None:
        """adds an CYAX object"""
        # only one CYAX card allowed
        assert self.model.cyax is None, '\ncyax=\n%s old=\n%s' % (cyax, self.model.cyax)
        self.model.cyax = cyax

    def _add_cyjoin_object(self, cyjoin: CYJOIN) -> None:
        """adds an CYJOIN object"""
        key = cyjoin.side
        assert key not in self.model.cyjoin, 'CYJOIN.side=%s\nold=\n%snew=\n%s' % (key, self.model.cyjoin[key], cyjoin)
        assert key >= 0
        self.model.cyjoin[key] = cyjoin
        self.model._type_to_id_map[cyjoin.type].append(key)

    def _add_modtrak_object(self, modtrak: MODTRAK) -> None:
        """adds an MODTRAK object"""
        # only one CYAX card allowed
        assert self.model.modtrak is None, '\nmodtrak=\n%s old=\n%s' % (modtrak, self.model.modtrak)
        self.model.modtrak = modtrak

    def _add_aefact_object(self, aefact: AEFACT, allow_overwrites: bool=False) -> None:
        """adds an AEFACT object"""
        key = aefact.sid
        if key in self.model.aefacts and not allow_overwrites:
            if not aefact == self.model.aefacts[key]:
                assert key not in self.model.aefacts, 'AEFACT.sid=%s\nold=\n%snew=\n%s' % (key, self.model.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.model.aefacts[key] = aefact
            self.model._type_to_id_map[aefact.type].append(key)

    def _add_aelist_object(self, aelist: AELIST) -> None:
        """adds an AELIST object"""
        key = aelist.sid
        assert key not in self.model.aelists, 'AELIST.sid=%s\nold=\n%snew=\n%s' % (key, self.model.aelists[key], aelist)
        assert key >= 0
        self.model.aelists[key] = aelist
        self.model._type_to_id_map[aelist.type].append(key)

    def _add_aelink_object(self, aelink: AELINK) -> None:
        """adds an AELINK object"""
        key = aelink.aelink_id
        assert key >= 0
        if key not in self.model.aelinks:
            self.model.aelinks[key] = []
        self.model.aelinks[key].append(aelink)
        self.model._type_to_id_map[aelink.type].append(key)
        #assert key not in self.model.aestats,'\naestat=%s oldAESTAT=\n%s' %(aestat,self.model.aestats[key])

    def _add_aecomp_object(self, aecomp: Union[AECOMP, AECOMPL]) -> None:
        """adds an AECOMP object"""
        key = aecomp.name
        assert key not in self.model.aecomps, '\naecomp=\n%s oldAECOMP=\n%s' % (aecomp, self.model.aecomps[key])
        self.model.aecomps[key] = aecomp
        self.model._type_to_id_map[aecomp.type].append(key)

    def _add_aeparm_object(self, aeparam: AEPARM) -> None:
        """adds an AEPARM object"""
        key = aeparam.aeparm_id
        assert key not in self.model.aeparams, '\naeparam=\n%s oldAEPARM=\n%s' % (aeparam, self.model.aeparams[key])
        assert key >= 0
        self.model.aeparams[key] = aeparam
        self.model._type_to_id_map[aeparam.type].append(key)

    def _add_aestat_object(self, aestat: AESTAT) -> None:
        """adds an AESTAT object"""
        key = aestat.aestat_id
        assert key not in self.model.aestats, '\naestat=\n%s old=\n%s' % (
            aestat, self.model.aestats[key])
        assert key >= 0
        self.model.aestats[key] = aestat
        self.model._type_to_id_map[aestat.type].append(key)

    def _add_aesurf_object(self, aesurf: AESURF) -> None:
        """adds an AESURF object"""
        key = aesurf.aesid
        assert key not in self.model.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, self.model.aesurf[key])
        assert key >= 0
        self.model.aesurf[key] = aesurf
        self.model._type_to_id_map[aesurf.type].append(key)

    def _add_aesurfs_object(self, aesurfs: AESURFS) -> None:
        """adds an AESURFS object"""
        key = aesurfs.aesid
        assert key not in self.model.aesurfs, '\naesurfs=\n%s old=\n%s' % (
            aesurfs, self.model.aesurfs[key])
        assert key >= 0
        self.model.aesurfs[key] = aesurfs
        self.model._type_to_id_map[aesurfs.type].append(key)

    def _add_csschd_object(self, csschd: CSSCHD) -> None:
        """adds an CSSCHD object"""
        key = csschd.sid
        assert key not in self.model.csschds, '\nCSSCHD=\n%s old=\n%s' % (csschd, self.model.csschds[key])
        assert key >= 0
        self.model.csschds[key] = csschd
        self.model._type_to_id_map[csschd.type].append(key)

    def _add_caero_object(self, caero: Union[CAERO1, CAERO2, CAERO3, CAERO4, CAERO5]) -> None:
        """adds an CAERO1/CAERO2/CAERO3/CAERO4/CAERO5 object"""
        key = caero.eid
        assert key not in self.model.caeros, '\nkey=%s; caero=\n%r old_caero=\n%r' % (
            key, caero, self.model.caeros[key])
        assert key > 0
        self.model.caeros[key] = caero
        self.model._type_to_id_map[caero.type].append(key)

    def _add_paero_object(self, paero: Union[PAERO1, PAERO2, PAERO3, PAERO4, PAERO5]) -> None:
        """adds an PAERO1/PAERO2/PAERO3/PAERO4/PAERO5 object"""
        key = paero.pid
        assert key not in self.model.paeros, '\npaero=\n%r old_paero=\n%r' % (
            paero, self.model.paeros[key])
        assert key > 0, 'paero.pid = %r' % (key)
        self.model.paeros[key] = paero
        self.model._type_to_id_map[paero.type].append(key)

    def _add_monpnt_object(self, monitor_point: Union[MONPNT1, MONPNT2, MONPNT3]) -> None:
        """adds an MONPNT object"""
        key = monitor_point.name
        assert key not in self.model.monitor_points, f'\nmonitor_point:\n{monitor_point}oldMOTPNT:\n{self.model.monitor_points[key]}'
        self.model.monitor_points.append(monitor_point)
        self.model._type_to_id_map[monitor_point.type].append(len(self.model.monitor_points) - 1)

    def _add_spline_object(self, spline: Union[SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5]) -> None:
        """adds an SPLINE1/SPLINE2/SPLINE3/SPLINE4/SPLINE5 object"""
        key = spline.eid
        assert spline.eid not in self.model.splines, f'\nspline:\n{spline}\nold_spline:\n{self.model.splines[key]}'
        assert spline.eid > 0, spline
        self.model.splines[key] = spline
        self.model._type_to_id_map[spline.type].append(key)

    def _add_gust_object(self, gust: GUST) -> None:
        """adds an GUST object"""
        key = gust.sid
        assert key not in self.model.gusts
        assert key > 0
        self.model.gusts[key] = gust
        self.model._type_to_id_map[gust.type].append(key)

    def _add_trim_object(self, trim: Union[TRIM, TRIM2], allow_overwrites: bool=False) -> None:
        """adds an TRIM object"""
        key = trim.sid
        if not allow_overwrites:
            assert key not in self.model.trims, 'TRIM=%s  old=\n%snew=\n%s' % (key, self.model.trims[key], trim)
        assert key > 0, 'key=%r trim=\n%s' % (key, trim)
        self.model.trims[key] = trim
        self.model._type_to_id_map[trim.type].append(key)

    def _add_diverg_object(self, diverg: DIVERG, allow_overwrites: bool=False) -> None:
        """adds an DIVERG object"""
        key = diverg.sid
        if not allow_overwrites:
            assert key not in self.model.divergs, 'DIVERG=%s  old=\n%snew=\n%s' % (key, self.model.divergs[key], diverg)
        assert key > 0, 'key=%r diverg=\n%s' % (key, diverg)
        self.model.divergs[key] = diverg
        self.model._type_to_id_map[diverg.type].append(key)

    def _add_flutter_object(self, flutter: FLUTTER) -> None:
        """adds an FLUTTER object"""
        key = flutter.sid
        assert key not in self.model.flutters, 'FLUTTER=%s old=\n%snew=\n%s' % (key, self.model.flutters[key], flutter)
        assert key > 0
        self.model.flutters[key] = flutter
        self.model._type_to_id_map[flutter.type].append(key)

    def _add_flfact_object(self, flfact: FLFACT) -> None:
        """adds an FLFACT object"""
        key = flfact.sid
        #assert key not in self.model.flfacts
        assert key > 0
        self.model.flfacts[key] = flfact  # set id...
        self.model._type_to_id_map[flfact.type].append(key)

    def _add_dconstr_object(self, dconstr: Union[DCONSTR, DCONADD]) -> None:
        """adds a DCONSTR object"""
        #key = (dconstr.oid, dconstr.rid)
        key = dconstr.oid
        #assert key not in self.model.dconstrs, 'key=%r DCONSTR/DCONADD=\n%s' % (key, dconstr)
        assert dconstr.oid > 0
        #assert dconstr.dresp_id > 0
        if key in self.model.dconstrs:
            self.model.dconstrs[key].append(dconstr)
        else:
            self.model.dconstrs[key] = [dconstr]
        self.model._type_to_id_map[dconstr.type].append(key)

    #def add_DCONADD(self, dconadd, allow_overwrites: bool=False) -> None:
        #key = dconadd.oid
        #if key in self.model.dconstrs and not allow_overwrites:
            #if not dconadd == self.model.dconstrs[key]:
                #assert key not in self.model.dconstrs, 'DCONADD=%s old=\n%snew=\n%s' % (
                    #key, self.model.dconstrs[key], dconadd)
        #else:
            #assert key > 0, 'dcid=%s dconadd=%s' % (key, dconadd)
            #self.model.dconstrs[key] = dconadd
            #self.model._type_to_id_map[dconadd.type].append(key)

    def _add_desvar_object(self, desvar: DESVAR) -> None:
        """adds a DESVAR object"""
        key = desvar.desvar_id
        assert key not in self.model.desvars, 'DESVAR=%s old=\n%snew=\n%s' % (
            key, self.model.desvars[key], desvar)
        assert key > 0
        self.model.desvars[key] = desvar
        self.model._type_to_id_map[desvar.type].append(key)

    def _add_topvar_object(self, topvar: TOPVAR) -> None:
        """adds a TOPVAR object"""
        key = topvar.opt_id
        assert key not in self.model.topvar, 'TOPVAR=%s old=\n%snew=\n%s' % (
            key, self.model.topvar[key], topvar)
        assert key > 0
        self.model.topvar[key] = topvar
        self.model._type_to_id_map[topvar.type].append(key)

    def _add_ddval_object(self, ddval: DDVAL) -> None:
        """adds a DDVAL object"""
        key = ddval.oid
        assert key not in self.model.ddvals, 'DDVAL=%s old=\n%snew=\n%s' % (
            key, self.model.ddvals[key], ddval)
        assert key > 0
        self.model.ddvals[key] = ddval
        self.model._type_to_id_map[ddval.type].append(key)

    def _add_dlink_object(self, dlink: DLINK) -> None:
        """adds a DLINK object"""
        key = dlink.oid
        assert key not in self.model.dlinks, 'DLINK=%s old=\n%snew=\n%s' % (
            key, self.model.dlinks[key], dlink)
        assert key > 0
        self.model.dlinks[key] = dlink
        self.model._type_to_id_map[dlink.type].append(key)

    def _add_dscreen_object(self, dscreen: DSCREEN) -> None:
        """adds a DSCREEN object"""
        key = dscreen.rtype
        assert key not in self.model.dscreen, 'DSCREEN=%s old=\n%snew=\n%s' % (
            key, self.model.dscreen[key], dscreen)
        assert len(key) > 0, 'key=%r' % key
        self.model.dscreen[key] = dscreen
        self.model._type_to_id_map[dscreen.type].append(key)

    def _add_dresp_object(self, dresp: Union[DRESP1, DRESP2, DRESP3]) -> None:
        """adds a DRESP1/DRESP2/DRESP3 object"""
        key = dresp.dresp_id
        assert key not in self.model.dresps, 'DRESPx=%s old=\n%snew=\n%s' % (
            key, self.model.dresps[key], dresp)
        assert key > 0
        self.model.dresps[key] = dresp
        self.model._type_to_id_map[dresp.type].append(key)

    def _add_dvcrel_object(self, dvcrel: Union[DVCREL1, DVCREL2]) -> None:
        """adds a DVCREL1/DVCREL2 object"""
        key = dvcrel.oid
        assert key not in self.model.dvcrels, 'DVCRELx=%s old\n%snew=\n%s' % (
            key, self.model.dvcrels[key], dvcrel)
        assert key > 0
        self.model.dvcrels[key] = dvcrel
        self.model._type_to_id_map[dvcrel.type].append(key)

    def _add_dvmrel_object(self, dvmrel: Union[DVMREL1, DVMREL2]) -> None:
        """adds a DVMREL1/DVMREL2 object"""
        key = dvmrel.oid
        assert key not in self.model.dvmrels, 'DVMRELx=%s old=\n%snew=\n%s' % (
            key, self.model.dvmrels[key], dvmrel)
        assert key not in self.model.dvmrels
        assert key > 0
        self.model.dvmrels[key] = dvmrel
        self.model._type_to_id_map[dvmrel.type].append(key)

    def _add_dvprel_object(self, dvprel: Union[DVPREL1, DVPREL2]) -> None:
        """adds a DVPREL1/DVPREL2 object"""
        key = dvprel.oid
        assert key not in self.model.dvprels, 'DVPRELx=%s old\n%snew=\n%s' % (
            key, self.model.dvprels[key], dvprel)
        assert key > 0
        self.model.dvprels[key] = dvprel
        self.model._type_to_id_map[dvprel.type].append(key)

    def _add_dvgrid_object(self, dvgrid: DVGRID) -> None:
        """adds a DVGRID object"""
        key = dvgrid.dvid
        assert key > 0
        if key not in self.model.dvgrids:
            self.model.dvgrids[key] = []
            self.model._type_to_id_map[dvgrid.type].append(key)
        self.model.dvgrids[key].append(dvgrid)

    #-----------------------------------------------------------
    # nx optimization
    def _add_dvtrel_object(self, dvtrel: DVTREL1) -> None: # DVTREL2
        """adds a DVTREL1/DVTREL2 object"""
        key = dvtrel.dvtrel_id
        assert key not in self.model.dvtrels, 'DVTRELx=%s old\n%snew=\n%s' % (
            key, self.model.dvtrels[key], dvtrel)
        assert key > 0
        self.model.dvtrels[key] = dvtrel
        self.model._type_to_id_map[dvtrel.type].append(key)

    def _add_group_object(self, group: GROUP) -> None:
        """adds a GROUP object"""
        key = group.group_id
        assert key > 0
        self.model.group[key] = group
        self.model._type_to_id_map[group.type].append(key)

    def _add_dmncon_object(self, dmncon: DMNCON) -> None:
        """adds a DMNCON object"""
        key = dmncon.constraint_id
        assert key > 0
        self.model.dmncon[key] = dmncon
        self.model._type_to_id_map[dmncon.type].append(key)

    #-----------------------------------------------------------

    def _add_nlparm_object(self, nlparm: NLPARM) -> None:
        """adds a NLPARM object"""
        key = nlparm.nlparm_id
        assert key not in self.model.nlparms
        assert key > 0, 'key=%s; nlparm=%s\n' % (key, nlparm)
        self.model.nlparms[key] = nlparm
        self.model._type_to_id_map[nlparm.type].append(key)

    def _add_rotor_object(self, rotor: Union[ROTORD, ROTORG]) -> None:
        """adds a ROTORD/ROTORG object"""
        key = rotor.sid
        assert key > 0, 'key=%s; rotor=%s\n' % (key, rotor)
        if key in self.model.rotors:
            rotor_old = self.model.rotors[key]
            assert rotor.type == rotor_old.type
            self.model.rotors[key].nids += rotor.nids
        else:
            self.model.rotors[key] = rotor
        self.model._type_to_id_map[rotor.type].append(key)

    def _add_nlpci_object(self, nlpci: NLPCI) -> None:
        """adds a NLPCI object"""
        key = nlpci.nlpci_id
        assert key not in self.model.nlpcis
        assert key > 0
        self.model.nlpcis[key] = nlpci
        self.model._type_to_id_map[nlpci.type].append(key)

    def _add_nxstrat_object(self, nxstrat: NXSTRAT) -> None:
        key = nxstrat.sid
        assert key not in self.model.nxstrats, 'nxstrats=%s nxstrat=%s' % (self.model.nxstrats, nxstrat)
        assert key > 0
        self.model.nxstrats[key] = nxstrat
        self.model._type_to_id_map[nxstrat.type].append(key)

    def _add_tstep_object(self, tstep: Union[TSTEP, TSTEP1],
                          allow_overwrites: bool=False) -> None:
        """adds a TSTEP object"""
        key = tstep.sid
        if key in self.model.tsteps and not allow_overwrites:
            if not tstep == self.model.tsteps[key]:
                assert key not in self.model.tsteps, 'TSTEP=%s\nold=\n%snew=\n%s' % (key, self.model.tsteps[key], tstep)
        else:
            assert key > 0, 'sid=%s tstep=\n%s' % (key, tstep)
            self.model.tsteps[key] = tstep
            self.model._type_to_id_map[tstep.type].append(key)

    def _add_tstepnl_object(self, tstepnl: TSTEPNL,
                            allow_overwrites: bool=False) -> None:
        """adds a TSTEPNL object"""
        key = tstepnl.sid
        if key in self.model.tstepnls and not allow_overwrites:
            if not tstepnl == self.model.tstepnls[key]:
                assert key not in self.model.tstepnls, 'TSTEPNL=%s\nold=\n%snew=\n%s' % (key, self.model.tstepnls[key], tstepnl)
        else:
            assert key > 0, 'sid=%s tstepnl=\n%s' % (key, tstepnl)
            self.model.tstepnls[key] = tstepnl
            self.model._type_to_id_map[tstepnl.type].append(key)

    def _add_freq_object(self, freq: Union[FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5]) -> None:
        key = freq.sid
        assert key > 0
        if key in self.model.frequencies:
            freq0 = self.model.frequencies[key][0]
            if freq0.type == 'FREQ' and freq.type == 'FREQ':
                freq0.add_frequency_object(freq)
            else:
                self.model.frequencies[key].append(freq)
        else:
            self.model.frequencies[key] = [freq]
            self.model._type_to_id_map[freq.type].append(key)

    def _add_set_object(self, set_obj: Union[SET1, SET2, SET3]) -> None:
        """adds an SET1/SET3 object"""
        key = set_obj.sid
        assert key >= 0
        if key in self.model.sets:
            self.model.sets[key].add_set(set_obj)
        else:
            self.model.sets[key] = set_obj
            self.model._type_to_id_map[set_obj.type].append(key)

    def _add_radset_object(self, set_obj: RADSET) -> None:
        """adds an RADSET object"""
        if self.model.radset:
            self.model.radset.add_set(set_obj)
        else:
            self.model.radset = set_obj
            #self.model._type_to_id_map[set_obj.type].append(key)

    def _add_aset_object(self, set_obj: Union[ASET, ASET1]) -> None:
        """adds an ASET/ASET1 object"""
        self.model.asets.append(set_obj)
        n = len(self.model._type_to_id_map['ASET'])
        self.model._type_to_id_map['ASET'].append(n)

    def _add_omit_object(self, set_obj: Union[OMIT, OMIT1]) -> None:
        """adds an OMIT/OMIT1 object"""
        self.model.omits.append(set_obj)
        n = len(self.model._type_to_id_map['OMIT'])
        self.model._type_to_id_map['OMIT'].append(n)

    def _add_bset_object(self, set_obj: Union[BSET, BSET1]) -> None:
        """adds an BSET/BSET1 object"""
        self.model.bsets.append(set_obj)
        n = len(self.model._type_to_id_map['BSET'])
        self.model._type_to_id_map['BSET'].append(n)

    def _add_cset_object(self, set_obj: Union[CSET, CSET1]) -> None:
        """adds an CSET/USET1 object"""
        self.model.csets.append(set_obj)
        n = len(self.model._type_to_id_map['CSET'])
        self.model._type_to_id_map['CSET'].append(n)

    def _add_qset_object(self, set_obj: Union[QSET, QSET1]) -> None:
        """adds an QSET/QSET1 object"""
        self.model.qsets.append(set_obj)
        n = len(self.model._type_to_id_map['QSET'])
        self.model._type_to_id_map['QSET'].append(n)

    def _add_uset_object(self, set_obj: Union[USET, USET1]) -> None:
        """adds an USET/USET1 object"""
        key = set_obj.name
        if key in self.model.usets:
            self.model.usets[key].append(set_obj)
        else:
            self.model.usets[key] = [set_obj]
        self.model._type_to_id_map[set_obj.type].append(key)

    def _add_sebset_object(self, set_obj: Union[SEBSET, SEBSET1]) -> None:
        """adds an SEBSET/SEBSET1 object"""
        self.model.se_bsets.append(set_obj)

    def _add_secset_object(self, set_obj: Union[SECSET, SECSET1]) -> None:
        """adds an SECSET/SECSTE1 object"""
        self.model.se_csets.append(set_obj)

    def _add_seqset_object(self, set_obj: Union[SEQSET, SEQSET1]) -> None:
        """adds an SEQSET/SEQSET1 object"""
        self.model.se_qsets.append(set_obj)

    def _add_seuset_object(self, set_obj: Union[SEUSET, SEUSET1]) -> None:
        """adds an SEUSET/SEUSET1 object"""
        key = set_obj.name
        if key in self.model.se_usets:
            self.model.se_usets[key].append(set_obj)
        else:
            self.model.se_usets[key] = [set_obj]
        self.model._type_to_id_map[set_obj.type].append(key)

    def _add_seset_object(self, set_obj: SESET) -> None:
        """adds an SESET object"""
        key = set_obj.seid
        assert key >= 0
        if key in self.model.se_sets:
            old_set = self.model.se_sets[key]
            set_obj.add_seset(old_set)
        self.model.se_sets[key] = set_obj
        self.model._type_to_id_map[set_obj.type].append(key)

    def _add_table_object(self, table: Union[TABLEH1, TABLEHT, TABLES1, TABLEST]) -> None:
        """adds a TABLES1, TABLEST object"""
        key = table.tid
        if key in self.model.tables:
            if not table == self.model.tables[key]:
                assert key not in self.model.tables, '\ntable=\n%s old_table=\n%s' % (
                    table, self.model.tables[key])
        assert key > 0
        self.model.tables[key] = table
        self.model._type_to_id_map[table.type].append(key)

    def _add_tabled_object(self, table: Union[TABLED1, TABLED2, TABLED3, TABLED4]) -> None:
        """adds a TABLED1, TABLED2, TABLED3, TABLED4 object"""
        key = table.tid
        assert key not in self.model.tables_d, '\ntabled=\n%s old_tabled=\n%s' % (
            table, self.model.tables_d[key])
        #assert key > 0; yes you can have negative tables...
        self.model.tables_d[key] = table
        self.model._type_to_id_map[table.type].append(key)

    def _add_tablem_object(self, table: Union[TABLEM1, TABLEM2, TABLEM3, TABLEM4]) -> None:
        """adds a TABLEM1, TABLEM2, TABLEM3, TABLEM4 object"""
        key = table.tid
        assert key not in self.model.tables_m, '\ntablem=\n%s old_tablem=\n%s' % (
            table, self.model.tables_m[key])
        #assert key > 0; yes you can have negative tables...
        self.model.tables_m[key] = table
        self.model._type_to_id_map[table.type].append(key)

    def _add_table_sdamping_object(self, table: TABDMP1) -> None:
        """adds a TABDMP1 object"""
        key = table.tid
        assert key not in self.model.tables_sdamping, '\nTable=\n%s oldTable=\n%s' % (
            table, self.model.tables_sdamping[key])
        #assert key > 0; yes you can have negative tables...
        self.model.tables_sdamping[key] = table
        self.model._type_to_id_map[table.type].append(key)

    def _add_random_table_object(self, table: Union[TABRND1, TABRNDG]) -> None:
        """adds a TABRND1, TABRNDG object"""
        key = table.tid
        assert key not in self.model.random_tables, '\nTable=\n%s old=\n%s' % (
            table, self.model.random_tables[key])
        assert key > 0
        self.model.random_tables[key] = table
        self.model._type_to_id_map[table.type].append(key)

    def _add_method_object(self, method: Union[EIGR, EIGRL, EIGB],
                           allow_overwrites: bool=False) -> None:
        """adds a EIGR/EIGRL object"""
        key = method.sid
        if key in self.model.methods and not allow_overwrites:
            if not method == self.model.methods[key]:
                assert key not in self.model.methods, 'sid=%s\nold_method=\n%snew_method=\n%s' % (key, self.model.methods[key], method)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, method)
            self.model.methods[key] = method
            self.model._type_to_id_map[method.type].append(key)

    def _add_cmethod_object(self, method: Union[EIGC, EIGP],
                            allow_overwrites: bool=False) -> None:
        """adds a EIGB/EIGC object"""
        key = method.sid
        if key in self.model.cMethods and not allow_overwrites:
            if not method == self.model.cMethods[key]:
                assert key not in self.model.cMethods, 'sid=%s\nold_cmethod=\n%snew_cmethod=\n%s' % (key, self.model.cMethods[key], method)
        else:
            assert key > 0, 'sid=%s cMethod=\n%s' % (key, method)
            self.model.cMethods[key] = method
            self.model._type_to_id_map[method.type].append(key)

    def _add_mkaero_object(self, mkaero: Union[MKAERO1, MKAERO2]) -> None:
        """adds an MKAERO1/MKAERO2 object"""
        self.model.mkaeros.append(mkaero)

    #---------------------------------------------------------------------------
    # parametric
    def _add_pset(self, pset: PSET, allow_overwrites: bool=False) -> None:
        assert pset.idi not in self.model.pset, pset
        self.model.pset[pset.idi] = pset

    def _add_pval(self, pval: PVAL, allow_overwrites: bool=False) -> None:
        if pval.idi not in self.model.pval:
            self.model.pval[pval.idi] = []
        self.model.pval[pval.idi].append(pval)

    def _add_gmcurv(self, curve: GMCURV, allow_overwrites: bool=False) -> None:
        assert curve.curve_id not in self.model.gmcurv, curve
        self.model.gmcurv[curve.curve_id] = curve

    def _add_gmsurf(self, surf: GMSURF, allow_overwrites: bool=False) -> None:
        assert surf.surf_id not in self.model.gmsurf, surf
        self.model.gmsurf[surf.surf_id] = surf

    def _add_feface(self, face: FEFACE, allow_overwrites: bool=False) -> None:
        key = face.face_id
        if key in self.model.feface and not allow_overwrites:
            if not face == self.model.feface[key]:
                raise RuntimeError(f'feface is duplicated\n{face}\nold:\n{self.model.feface[key]}')
        else:
            self.model.feface[face.face_id] = face
            self.model._type_to_id_map[face.type].append(key)
        #assert face.face_id not in self.model.feface, face
        #self.model.feface[face.face_id] = face

    def _add_feedge(self, edge: FEEDGE, allow_overwrites: bool=False) -> None:
        key = edge.edge_id
        if key in self.model.feedge and not allow_overwrites:
            if not edge == self.model.feedge[key]:
                raise RuntimeError(f'feedge is duplicated\n{edge}\nold:\n{self.model.feedge[key]}')
        else:
            self.model.feedge[edge.edge_id] = edge
            self.model._type_to_id_map[edge.type].append(key)

    #---------------------------------------------------------------------------
    # nx bolts
    def _add_bolt_object(self, bolt: BOLT, allow_overwrites: bool=False) -> None:
        key = bolt.bolt_id
        if key in self.model.bolt and not allow_overwrites:
            if not bolt == self.model.bolt[key]:
                raise RuntimeError(f'bolt is duplicated\n{bolt}\nold:\n{self.model.bolt[key]}')
        else:
            self.model.bolt[bolt.bolt_id] = bolt
            self.model._type_to_id_map[bolt.type].append(key)

    def _add_boltseq_object(self, boltseq: BOLTSEQ, allow_overwrites: bool=False) -> None:
        key = boltseq.sid
        if key in self.model.boltseq and not allow_overwrites:
            if not boltseq == self.model.boltseq[key]:
                raise RuntimeError(f'boltseq is duplicated\n{boltseq}\nold:\n{self.model.boltseq[key]}')
        else:
            self.model.boltseq[boltseq.sid] = boltseq
            self.model._type_to_id_map[boltseq.type].append(key)

    def _add_boltfor_object(self, boltfor: BOLTFOR, allow_overwrites: bool=False) -> None:
        key = boltfor.sid
        if key in self.model.boltfor and not allow_overwrites:
            if not boltfor == self.model.boltfor[key]:
                raise RuntimeError(f'boltfor is duplicated\n{boltfor}\nold:\n{self.model.boltfor[key]}')
        else:
            self.model.boltfor[boltfor.sid] = boltfor
            self.model._type_to_id_map[boltfor.type].append(key)
