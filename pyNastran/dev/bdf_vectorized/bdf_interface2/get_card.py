# pylint: disable=E1101,C0103
import numpy as np
from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.deprecated import GetMethodsDeprecated
#from pyNastran.bdf.cards.nodes import SPOINT, EPOINT
from pyNastran.dev.bdf_vectorized.bdf_interface2.attributes import BDFAttributes


class GetMethods(BDFAttributes):
    def __init__(self):
        self._type_to_slot_map = {}
        BDFAttributes.__init__(self)

    #--------------------
    # SET CARDS
    def Set(self, sid, msg=''):
        try:
            return self.sets[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed SETx=%s'
                           % (sid, msg, np.unique(list(self.sets.keys()))))

    def SET1(self, set_id, msg=''):
        """gets a SET1"""
        assert isinstance(set_id, integer_types), 'set_id=%s is not an integer\n' % set_id
        if set_id in self.sets:
            set1 = self.sets[set_id]
        else:
            raise KeyError('cannot find SET1 ID=%r.\n%s' % (set_id, msg))
        return set1

    #--------------------
    # COORDINATES CARDS
    def Coord(self, cid, msg=''):
        """gets an COORDx"""
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coord_ids))

    #--------------------
    # AERO CARDS

    def AEList(self, aelist, msg=''):
        """gets an AELIST"""
        try:
            return self.aelists[aelist]
        except KeyError:
            raise KeyError('aelist=%s not found%s.  Allowed AELIST=%s'
                           % (aelist, msg, np.unique(list(self.aelists.keys()))))

    def AEFact(self, aefact, msg=''):
        """gets an AEFACT"""
        try:
            return self.aefacts[aefact]
        except KeyError:
            raise KeyError('aefact=%s not found%s.  Allowed AEFACT=%s'
                           % (aefact, msg, np.unique(list(self.aefacts.keys()))))

    def Acsid(self, msg=''):
        """gets the aerodynamic system coordinate"""
        if self.aero is not None:
            acsid_aero = self.aero.Acsid()
        if self.aeros is not None:
            acsid_aeros = self.aeros.Acsid()

        if self.aero is not None:
            if self.aeros is not None:
                assert acsid_aero == acsid_aeros, 'AERO acsid=%s, AEROS acsid=%s' % (acsid_aero,
                                                                                     acsid_aeros)
            cid = self.Coord(acsid_aero, msg=msg)
        elif self.aeros is not None:
            cid = self.Coord(acsid_aeros, msg=msg)
        else:
            msg = 'neither AERO nor AEROS cards exist.'
            raise RuntimeError(msg)
        return cid

    def Aero(self, msg=''):
        """gets the AERO"""
        if self.aero is not None:
            return self.aero
        else:
            raise RuntimeError(f'no AERO card found{msg}.')

    def Aeros(self, msg=''):
        """gets the AEROS"""
        if self.aeros is not None:
            return self.aeros
        else:
            raise RuntimeError(f'no AEROS card found{msg}.')

    def Spline(self, eid, msg=''):
        """gets a SPLINEx"""
        try:
            return self.splines[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed SPLINEx=%s'
                           % (eid, msg, np.unique(list(self.splines.keys()))))

    def CAero(self, eid, msg=''):
        """gets an CAEROx"""
        try:
            return self.caeros[eid]
        except KeyError:
            raise KeyError('eid=%s not found%s.  Allowed CAEROx=%s'
                           % (eid, msg, np.unique(list(self.caeros.keys()))))

    def PAero(self, pid, msg=''):
        """gets a PAEROx"""
        try:
            return self.paeros[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed PAEROx=%s'
                           % (pid, msg, np.unique(list(self.paeros.keys()))))

    def Gust(self, sid, msg=''):
        """gets a GUST"""
        try:
            return self.gusts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed GUSTs=%s'
                           % (sid, msg, np.unique(list(self.gusts.keys()))))

    #--------------------
    # AERO CONTROL SURFACE CARDS
    def AEStat(self, aid, msg=''):
        """gets an AESTAT"""
        try:
            return self.aestats[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AESTATs=%s'
                           % (aid, msg, np.unique(list(self.aestats.keys()))))

    def AELIST(self, aid, msg=''):
        """gets an AELIST"""
        try:
            return self.aelists[aid]
        except KeyError:
            raise KeyError('id=%s not found%s.  Allowed AELISTs=%s'
                           % (aid, msg, np.unique(list(self.aelists.keys()))))

    def AELink(self, link_id, msg=''):
        """gets an AELINK"""
        try:
            return self.aelinks[link_id]
        except KeyError:
            raise KeyError('link_id=%s not found%s.  Allowed AELINKs=%s'
                           % (link_id, msg, np.unique(list(self.aelinks.keys()))))

    def AEParam(self, aid, msg=''):
        """gets an AEPARM"""
        try:
            return self.aeparams[aid]
        except KeyError:
            raise KeyError('aid=%s not found%s.  Allowed AEPARMs=%s'
                           % (aid, msg, np.unique(list(self.aeparams.keys()))))

    #--------------------
    # FLUTTER CARDS

    def FLFACT(self, sid, msg=''):
        """gets an FLFACT"""
        try:
            return self.flfacts[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           % (sid, msg, np.unique(list(self.flfacts.keys()))))

    def Flutter(self, fid, msg=''):
        """gets a FLUTTER"""
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, np.unique(list(self.flutters.keys()))))

    #--------------------
    # OPTIMIZATION CARDS

    def DConstr(self, oid, msg=''):
        """gets a DCONSTR"""
        try:
            return self.dconstrs[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DCONSTRs=%s'
                           % (oid, msg, np.unique(list(self.dconstrs.keys()))))

    def DResp(self, dresp_id, msg=''):
        """gets a DRESPx"""
        try:
            return self.dresps[dresp_id]
        except KeyError:
            raise KeyError('dresp_id=%s not found%s.  Allowed DRESPx=%s'
                           % (dresp_id, msg, np.unique(list(self.dresps.keys()))))

    def Desvar(self, desvar_id, msg=''):
        """gets a DESVAR"""
        try:
            return self.desvars[desvar_id]
        except KeyError:
            raise KeyError('desvar_id=%s not found%s.  Allowed DESVARs=%s'
                           % (desvar_id, msg, np.unique(list(self.desvars.keys()))))

    def DDVal(self, oid, msg=''):
        """gets a DDVAL"""
        try:
            return self.ddvals[oid]
        except KeyError:
            raise KeyError('oid=%s not found%s.  Allowed DDVALs=%s'
                           % (oid, msg, np.unique(list(self.ddvals.keys()))))

    def DVcrel(self, dv_id, msg=''):
        """gets a DVCREL1/DVCREL2"""
        try:
            return self.dvcrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVCRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvcrels.keys()))))

    def DVmrel(self, dv_id, msg=''):
        """gets a DVMREL1/DVMREL2"""
        try:
            return self.dvmrels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVMRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvmrels.keys()))))

    def DVprel(self, dv_id, msg=''):
        """gets a DVPREL1/DVPREL2"""
        try:
            return self.dvprels[dv_id]
        except KeyError:
            raise KeyError('dv_id=%s not found%s.  Allowed DVPRELx=%s'
                           % (dv_id, msg, np.unique(list(self.dvprels.keys()))))

    #--------------------
    # METHOD CARDS
    def Method(self, sid, msg=''):
        """gets a METHOD (EIGR, EIGRL)"""
        try:
            return self.methods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed METHODs=%s'
                           % (sid, msg, np.unique(list(self.methods.keys()))))

    def CMethod(self, sid, msg=''):
        """gets a METHOD (EIGC)"""
        try:
            return self.cmethods[sid]
        except KeyError:
            raise KeyError('sid=%s not found%s.  Allowed CMETHODs=%s'
                           % (sid, msg, np.unique(list(self.cmethods.keys()))))

    #--------------------
    # TABLE CARDS
    def Table(self, tid, msg=''):
        """gets a TABLEx (TABLED1, TABLED2, TABLD3)"""
        try:
            return self.tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, np.unique(list(self.tables.keys()))))

    def RandomTable(self, tid, msg=''):
        try:
            return self.random_tables[tid]
        except KeyError:
            raise KeyError('tid=%s not found%s.  Allowed TABLEs=%s'
                           % (tid, msg, np.unique(list(self.random_tables.keys()))))

    #--------------------
    # NONLINEAR CARDS

    def NLParm(self, nid, msg=''):
        """gets an NLPARM"""
        try:
            return self.nlparms[nid]
        except KeyError:
            raise KeyError('nid=%s not found%s.  Allowed NLPARMs=%s'
                           % (nid, msg, np.unique(list(self.nlparms.keys()))))

    #--------------------
    # MATRIX ENTRY CARDS
    def DMIG(self, dname, msg=''):
        """gets a DMIG"""
        try:
            return self.dmig[dname]
        except KeyError:
            raise KeyError('dname=%s not found%s.  Allowed DMIGs=%s'
                           % (dname, msg, np.unique(list(self.dmig.keys()))))

    def DEQATN(self, equation_id, msg=''):
        """gets a DEQATN"""
        try:
            return self.dequations[equation_id]
        except KeyError:
            raise KeyError('equation_id=%s not found%s.  Allowed DMIGs=%s'
                           % (equation_id, msg, np.unique(list(self.dequations.keys()))))
