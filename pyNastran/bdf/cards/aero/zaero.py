# coding: utf-8
# pylint: disable=W0212,C0103
from __future__ import annotations

from collections import defaultdict
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils import object_attributes, object_methods

from typing import Any
from pyNastran.bdf.bdf_interface.utils import sorteddict
from pyNastran.bdf.bdf_interface.utils import _prep_comment
from pyNastran.bdf.bdf_interface.pybdf import _clean_comment

from .zaero_cards.zaero_sets import SETADD
from pyNastran.bdf.cards.aero.zaero_cards.atm import ATMOS, FIXMATM, FIXHATM, FIXMACH, FIXMDEN
from pyNastran.bdf.cards.aero.zaero_cards.spline import (
    SPLINE1_ZAERO,
    SPLINE2_ZAERO,
    SPLINE3_ZAERO,
    SPLINEM,
)
from pyNastran.bdf.cards.aero.zaero_cards.geometry import (
    PANLST1,
    PANLST2,
    PANLST3,
    SEGMESH,
    CAERO7,
    BODY7,
    PAFOIL7,
    PAFOIL8,
    AESURFZ,
    AESLINK,
)
from pyNastran.bdf.cards.aero.zaero_cards.plot import (
    PLTAERO,
    PLTMODE,
    PLTVG,
    PLTFLUT,
    PLTTIME,
    PLTCP,
    PLTMIST,
    PLTSURF,
    PLTBODE,
    PLTTRIM,
)
from pyNastran.bdf.cards.aero.zaero_cards.flutter import FLUTTER_ZAERO, MKAEROZ
from pyNastran.bdf.cards.aero.zaero_cards.trim import (
    TRIM_ZAERO,
    TRIMVAR,
    TRIMLNK,
)
from pyNastran.bdf.cards.aero.zaero_cards.manuever import (
    MLOADS,
    EXTINP,
    EXTOUT,
    TRIMFNC,
    ACTU,
    LOADMOD,
    RBRED,
)
from pyNastran.bdf.cards.aero.zaero_cards.gust import GLOADS, DGUST, CGUST, MFTGUST
from pyNastran.bdf.cards.aero.zaero_cards.ase import (
    ASE,
    ASECONT,
    ASESNSR,
    ASESNS1,
    CJUNCT,
    CONCT,
    TFSET,
    MIMOSS,
    SISOTF,
    SENSET,
    SURFSET,
    CNCTSET,
    ASEGAIN,
    GAINSET,
    AEROLAG,
    MIMOTF,
    SENSR,
    GAIN,
    SUMBLK,
    DEADBN,
    DELAY_ZAERO,
    FILTFL,
    LIMTR,
)
from pyNastran.bdf.cards.aero.zaero_cards.bdf_tables import TABLED1_ZAERO, TABDMP1_ZAERO
from pyNastran.bdf.cards.aero.zaero_cards.dmi import DMIL
from pyNastran.bdf.cards.aero.zaero_cards.cards import (
    MLDPRNT,
    MLDSTAT,
    MINSTAT,
    MLDTRIM,
    MLDCOMD,
    MLDTIME,
    AEROZ,
    ACOORD,
    ATTACH,
    EXTFILE,
    CONMLST,
    CPFACT,
    APCONST,
    SPLINE0,
    PBODY7,
    TRIMFLT,
    ASEOUT,
    CMARGIN,
    OUTPUT4,
    CROSPSD,
    FOILSEC,
    GENGUST,
    MFTGUST,
    DMIS,
    CELLWNG,
    CELLBOX,
    INPCFD,
    OMITCFD,
    WT1AJJ,
    WT1FRC,
    WT2AJJ,
    WTUCP,
    TRIMOBJ,
    TRIMCON,
    APCNSND,
    APCNSCP,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

ZAERO_CARDS = [
    # atmosphere
    "ATMOS", "FIXHATM", "FIXMATM", "FIXMACH", "FIXMDEN",
    # already added
    "AEFACT", "CORD2R",
    # geometry
    "CAERO7", "AEROZ", "AESURFZ", "AESLINK", "ATTACH",
    "PANLST1", "PANLST2", "PANLST3",
    "PAFOIL7", "PAFOIL8",
    "BODY7", "PBODY7", "SEGMESH",
    "ACOORD", "MKAEROZ",
    "SPLINEM",
    # -------------
    # trim
    "PLTCP", "TRIMVAR", "TRIMLNK",
    "TRIMFNC",  # optimization
    "PLTTRIM",
    # -------------
    # flutter
    "FLUTTER", "FIXHATM", "FIXMATM",  # 'FIXMDEN',
    "PLTVG", "PLTFLUT", "PLTSURF",
    # -------------
    # plotting
    "PLTMODE", "PLTAERO",  #'PLTTIME',
    "PLTMIST",
    # -------------
    # mloads
    "MLOADS", "CJUNCT", "CONCT", "TFSET",
    "MLDSTAT", "MLDTRIM", "MLDPRNT", "MLDCOMD",
    "EXTINP", "EXTOUT",
    "LOADMOD", "PLTTIME",
    # -------------
    # gust
    "GLOADS", "DGUST", "CGUST",
    "GENGUST", "MFTGUST",
    # -------------
    # ase
    "ASE", "ASECONT", "ASESNSR", "ASESNS1",
    "SENSET", "ACTU", "MIMOSS", "SISOTF", "ASEGAIN",
    "GAINSET", "MIMOTF", "SENSR", "GAIN",
    "SUMBLK", "DEADBN", "DELAY", "FILTFL",
    "LIMTR", "PLTBODE", "AEROLAG",
    # -------------
    # other
    "SETADD", "DMIL", "EXTFILE", "MLDTIME", "MLDCOMD",
    "MINSTAT", "APCONST", "CONMLST", "CPFACT", "SPLINE0",
    "TRIMFLT", "ASEOUT", "CMARGIN", "OUTPUT4", "CROSPSD",
    "FOILSEC", "DMIS",
    "CELLWNG", "CELLBOX", "INPCFD", "OMITCFD",
    "WT1AJJ", "WT1FRC", "WT2AJJ", "WTUCP",
    "TRIMOBJ", "TRIMCON",
    "APCNSND", "APCNSCP",
    "RBRED",
    "CNCTSET", "SURFSET",
]


class AddMethods:
    def __init__(self, model: BDF):
        self.model = model

    @property
    def zona(self):
        return self.model.zaero

    def add_mldprnt_object(self, mldprnt: MLDPRNT) -> None:
        """adds an MLDPRNT object"""
        key = mldprnt.mldprnt_id
        if key in self.model.zaero.mldprnt:
            self.model.log.warning(f"skipping MLDPRNT\n{str(mldprnt)}")
            return
        assert key not in self.model.zaero.mldprnt, key
        assert key > 0, key
        self.model.zaero.mldprnt[key] = mldprnt
        self.model._type_to_id_map[mldprnt.type].append(key)

    def add_mldcomd_object(self, mldcomd: MLDCOMD) -> None:
        """adds an MLDTRIM object"""
        key = mldcomd.mldcomd_id
        assert key not in self.model.zaero.mldcomd, key
        assert key > 0, key
        self.model.zaero.mldcomd[key] = mldcomd
        self.model._type_to_id_map[mldcomd.type].append(key)

    def add_extfile_object(self, extfile: EXTFILE) -> None:
        """adds an EXTFILE object"""
        key = extfile.extfile_id
        assert key not in self.model.zaero.extfile, key
        assert key > 0, key
        self.model.zaero.extfile[key] = extfile
        self.model._type_to_id_map[extfile.type].append(key)

    def add_dmil_object(self, dmil: DMIL) -> None:
        """adds an DMIL object"""
        name = dmil.name
        key = (name, dmil.row, dmil.col)
        assert key not in self.model.zaero.dmil, key
        self.model.zaero.dmil[key] = dmil
        self.model._type_to_id_map[dmil.type].append(name)

    def add_mldtime_object(self, mldtime: MLDTIME) -> None:
        """adds an MLDTRIM object"""
        key = mldtime.mldtime_id
        assert key not in self.model.zaero.mldtime, key
        assert key > 0, key
        self.model.zaero.mldtime[key] = mldtime
        self.model._type_to_id_map[mldtime.type].append(key)

    def add_mldtrim_object(self, mldtrim: MLDTRIM) -> None:
        """adds an MLDTRIM object"""
        key = mldtrim.mldtrim_id
        assert key not in self.model.zaero.mldtrim, key
        assert key > 0, key
        self.model.zaero.mldtrim[key] = mldtrim
        self.model._type_to_id_map[mldtrim.type].append(key)

    def add_mldstat_object(self, mldstat: MLDSTAT) -> None:
        """adds an MLDSTAT object"""
        key = mldstat.mldstat_id
        assert key > 0, key
        zaero = self.model.zaero
        assert key not in zaero.mldstat, key
        zaero.mldstat[key] = mldstat
        self.model._type_to_id_map[mldstat.type].append(key)

    def add_minstat_object(self, minstat: MINSTAT) -> None:
        """adds an MINSTAT object"""
        key = minstat.minstat_id
        assert key not in self.model.zaero.minstat, key
        assert key > 0, key
        self.model.zaero.minstat[key] = minstat
        self.model._type_to_id_map[minstat.type].append(key)

    def add_conmlst_object(self, conmlst: CONMLST) -> None:
        """adds a CONMLST object"""
        key = conmlst.conmlst_id
        assert key > 0, key
        assert key not in self.model.zaero.conmlst, key
        self.model.zaero.conmlst[key] = conmlst
        self.model._type_to_id_map[conmlst.type].append(key)

    def add_cpfact_object(self, cpfact: CPFACT) -> None:
        """adds a CPFACT object"""
        key = cpfact.cpfact_id
        assert key > 0, key
        assert key not in self.model.zaero.cpfact, key
        self.model.zaero.cpfact[key] = cpfact
        self.model._type_to_id_map[cpfact.type].append(key)

    def add_apconst_object(self, apconst: APCONST) -> None:
        """adds an APCONST object"""
        key = apconst.sid
        assert key > 0, key
        assert key not in self.model.zaero.apconst, key
        self.model.zaero.apconst[key] = apconst
        self.model._type_to_id_map[apconst.type].append(key)

    def add_spline0_object(self, spline0: SPLINE0) -> None:
        """adds a SPLINE0 object"""
        key = spline0.eid
        assert key > 0, key
        assert key not in self.model.zaero.spline0, key
        self.model.zaero.spline0[key] = spline0
        self.model._type_to_id_map[spline0.type].append(key)

    def add_pbody7_object(self, pbody7: PBODY7) -> None:
        """adds a PBODY7 object"""
        key = pbody7.pid
        assert key > 0, key
        assert key not in self.model.zaero.pbody7, key
        self.model.zaero.pbody7[key] = pbody7
        self.model._type_to_id_map[pbody7.type].append(key)

    def add_trimflt_object(self, trimflt: TRIMFLT) -> None:
        """adds a TRIMFLT object"""
        key = trimflt.trimflt_id
        assert key > 0, key
        assert key not in self.model.zaero.trimflt, key
        self.model.zaero.trimflt[key] = trimflt
        self.model._type_to_id_map[trimflt.type].append(key)

    def _add_generic_card(self, card, storage_attr: str) -> None:
        """Add a generic ZAERO card to the specified storage dict."""
        key = card.card_id
        storage = getattr(self.model.zaero, storage_attr)
        if key in storage:
            # Allow duplicate ID=0 by using a counter
            key = len(storage) + 1
        storage[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def add_aseout_object(self, card) -> None:
        self._add_generic_card(card, "aseout")

    def add_cmargin_object(self, card) -> None:
        """adds a CMARGIN object"""
        key = card.cmargin_id
        assert key > 0, key
        assert key not in self.model.zaero.cmargin, key
        self.model.zaero.cmargin[key] = card
        self.model._type_to_id_map[card.type].append(key)

    def add_output4_object(self, card) -> None:
        self._add_generic_card(card, "output4")

    def add_crospsd_object(self, card) -> None:
        self._add_generic_card(card, "crospsd")

    def add_foilsec_object(self, card) -> None:
        self._add_generic_card(card, "foilsec")

    def add_gengust_card_object(self, card) -> None:
        self._add_generic_card(card, "gengust_card")

    def add_mftgust_object(self, card) -> None:
        self._add_generic_card(card, "mftgust")

    def add_dmis_object(self, card) -> None:
        self._add_generic_card(card, "dmis")

    def add_cellwng_object(self, card) -> None:
        self._add_generic_card(card, "cellwng")

    def add_cellbox_object(self, card) -> None:
        self._add_generic_card(card, "cellbox")

    def add_inpcfd_object(self, card) -> None:
        self._add_generic_card(card, "inpcfd")

    def add_omitcfd_object(self, card) -> None:
        self._add_generic_card(card, "omitcfd")

    def add_wt1ajj_object(self, card) -> None:
        self._add_generic_card(card, "wt1ajj")

    def add_wt1frc_object(self, card) -> None:
        self._add_generic_card(card, "wt1frc")

    def add_wt2ajj_object(self, card) -> None:
        self._add_generic_card(card, "wt2ajj")

    def add_wtucp_object(self, card) -> None:
        self._add_generic_card(card, "wtucp")

    def add_trimobj_object(self, card) -> None:
        self._add_generic_card(card, "trimobj_card")

    def add_trimcon_object(self, card) -> None:
        self._add_generic_card(card, "trimcon_card")

    def add_apcnsnd_object(self, card) -> None:
        self._add_generic_card(card, "apcnsnd")

    def add_apcnscp_object(self, card) -> None:
        self._add_generic_card(card, "apcnscp")

    def add_extinp_object(self, extinp: EXTINP) -> None:
        """adds an EXTINP object"""
        key = extinp.extinp_id
        assert key not in self.model.zaero.extinp, key
        assert key > 0, key
        self.model.zaero.extinp[key] = extinp
        self.model._type_to_id_map[extinp.type].append(key)

    def add_extout_object(self, extout: EXTOUT) -> None:
        """adds an EXTOUT object"""
        key = extout.extout_id
        assert key not in self.model.zaero.extout, key
        assert key > 0, key
        self.model.zaero.extout[key] = extout
        self.model._type_to_id_map[extout.type].append(key)

    def add_splinem_object(self, splinem: SPLINEM) -> None:
        """adds an SPLINEM object"""
        assert self.model.zaero.splinem is None, self.model.zaero.splinem
        self.model.zaero.splinem = splinem
        self.model._type_to_id_map[splinem.type].append(1)

    def add_trimfnc_object(self, trimfnc: TRIMFNC) -> None:
        """adds an TRIMFNC object"""
        key = trimfnc.trimfnc_id
        assert key not in self.model.zaero.trimfnc, key
        assert key > 0, key
        self.model.zaero.trimfnc[key] = trimfnc
        self.model._type_to_id_map[trimfnc.type].append(key)

    def add_mimoss_object(self, mimoss: MIMOSS) -> None:
        """adds an MIMOSS object"""
        key = mimoss.mimoss_id
        assert key not in self.model.zaero.mimoss, key
        assert key > 0, key
        self.model.zaero.mimoss[key] = mimoss
        self.model._type_to_id_map[mimoss.type].append(key)

    def add_mimotf_object(self, mimotf: MIMOTF) -> None:
        """adds a MIMOTF object"""
        key = mimotf.mimotf_id
        assert key not in self.model.zaero.mimotf, key
        assert key > 0, key
        self.model.zaero.mimotf[key] = mimotf
        self.model._type_to_id_map[mimotf.type].append(key)

    def add_sensr_object(self, sensr: SENSR) -> None:
        """adds a SENSR object"""
        key = sensr.sensr_id
        assert key not in self.model.zaero.sensr, key
        assert key > 0, key
        self.model.zaero.sensr[key] = sensr
        self.model._type_to_id_map[sensr.type].append(key)

    def add_gain_object(self, gain: GAIN) -> None:
        """adds a GAIN object"""
        key = gain.gain_id
        assert key not in self.model.zaero.gain, key
        assert key > 0, key
        self.model.zaero.gain[key] = gain
        self.model._type_to_id_map[gain.type].append(key)

    def add_sumblk_object(self, sumblk: SUMBLK) -> None:
        """adds a SUMBLK object"""
        key = sumblk.sumblk_id
        assert key not in self.model.zaero.sumblk, key
        assert key > 0, key
        self.model.zaero.sumblk[key] = sumblk
        self.model._type_to_id_map[sumblk.type].append(key)

    def add_deadbn_object(self, deadbn: DEADBN) -> None:
        """adds a DEADBN object"""
        key = deadbn.deadbn_id
        assert key not in self.model.zaero.deadbn, key
        assert key > 0, key
        self.model.zaero.deadbn[key] = deadbn
        self.model._type_to_id_map[deadbn.type].append(key)

    def add_delay_zaero_object(self, delay: DELAY_ZAERO) -> None:
        """adds a DELAY (ZAERO) object"""
        key = delay.delay_id
        assert key not in self.model.zaero.delay_zaero, key
        assert key > 0, key
        self.model.zaero.delay_zaero[key] = delay
        self.model._type_to_id_map[delay.type].append(key)

    def add_filtfl_object(self, filtfl: FILTFL) -> None:
        """adds a FILTFL object"""
        key = filtfl.filtfl_id
        assert key not in self.model.zaero.filtfl, key
        assert key > 0, key
        self.model.zaero.filtfl[key] = filtfl
        self.model._type_to_id_map[filtfl.type].append(key)

    def add_limtr_object(self, limtr: LIMTR) -> None:
        """adds a LIMTR object"""
        key = limtr.limtr_id
        assert key not in self.model.zaero.limtr, key
        assert key > 0, key
        self.model.zaero.limtr[key] = limtr
        self.model._type_to_id_map[limtr.type].append(key)

    def add_sisotf_object(self, sisotf: SISOTF) -> None:
        """adds an SISOTF object"""
        key = sisotf.sisotf_id
        assert key not in self.model.zaero.sisotf, key
        assert key > 0, key
        self.model.zaero.sisotf[key] = sisotf
        self.model._type_to_id_map[sisotf.type].append(key)

    def add_tfset_object(self, tfset: TFSET) -> None:
        """adds an TFSET object"""
        key = tfset.tfset_id
        assert key not in self.model.zaero.tfset, key
        assert key > 0, key
        self.model.zaero.tfset[key] = tfset
        self.model._type_to_id_map[tfset.type].append(key)

    def add_setadd_object(self, setadd: SETADD) -> None:
        """adds an SETADD object"""
        key = setadd.setadd_id
        assert key not in self.model.zaero.setadd, key
        assert key > 0, key
        self.model.zaero.setadd[key] = setadd
        self.model._type_to_id_map[setadd.type].append(key)

    def add_senset_object(self, senset: SENSET) -> None:
        """adds an SENSET object"""
        key = senset.senset_id
        assert key not in self.model.zaero.senset, key
        assert key > 0, key
        self.model.zaero.senset[key] = senset
        self.model._type_to_id_map[senset.type].append(key)

    def add_cnctset_object(self, cnctset: CNCTSET) -> None:
        """adds an CNCTSET object"""
        key = cnctset.cnctset_id
        assert key not in self.model.zaero.cnctset, key
        assert key > 0, key
        self.model.zaero.cnctset[key] = cnctset
        self.model._type_to_id_map[cnctset.type].append(key)

    def add_surfset_object(self, surfset: SURFSET) -> None:
        """adds an SURFSET object"""
        key = surfset.surfset_id
        assert key not in self.model.zaero.surfset, key
        assert key > 0, key
        self.model.zaero.surfset[key] = surfset
        self.model._type_to_id_map[surfset.type].append(key)

    def add_loadmod_object(self, loadmod: LOADMOD) -> None:
        """adds an LOADMOD object"""
        key = loadmod.loadmod_id
        assert key not in self.model.zaero.loadmod, key
        assert key > 0, key
        self.model.zaero.loadmod[key] = loadmod
        self.model._type_to_id_map[loadmod.type].append(key)

    def add_rbred_object(self, rbred: RBRED) -> None:
        """adds an LOADMOD object"""
        key = rbred.sid
        assert key not in self.model.zaero.rbred, key
        assert key > 0, key
        self.model.zaero.rbred[key] = rbred
        self.model._type_to_id_map[rbred.type].append(key)

    def add_panlst_object(self, panlst: PANLST1 | PANLST2 | PANLST3) -> None:
        """adds an PANLST1/PANLST2/PANLST3 object"""
        key = panlst.eid
        assert key > 0, key
        zaero = self.model.zaero
        # assert key not in zaero.panlsts, '\npanlst=\n%s old=\n%s' % (
        #     panlst, zaero.panlsts[key])
        if key not in zaero.panlsts:
            zaero.panlsts[key] = []
        zaero.panlsts[key].append(panlst)
        self.model._type_to_id_map[panlst.type].append(key)

    def add_pafoil_object(self, pafoil: PAFOIL7 | PAFOIL8) -> None:
        """adds an PAFOIL7/PAFOIL8 object"""
        key = pafoil.pid
        assert pafoil.pid > 0
        zaero = self.model.zaero
        assert key not in zaero.pafoil, "\npafoil7=\n%s old=\n%s" % (pafoil, zaero.pafoil[key])
        zaero.pafoil[key] = pafoil
        self.model._type_to_id_map[pafoil.type].append(key)

    def add_aesurfz_object(self, aesurf: AESURFZ) -> None:
        """adds an AESURFZ object"""
        key = aesurf.aesid
        model = self.model
        assert key not in model.aesurf, "\naesurf=\n%s old=\n%s" % (aesurf, model.aesurf[key])
        model.aesurf[key] = aesurf
        model._type_to_id_map[aesurf.type].append(key)

    def add_aeslink_object(self, aeslink: AESLINK) -> None:
        """adds an AESLINK object"""
        key = aeslink.label
        model = self.model
        zaero = model.zaero
        assert key not in zaero.aeslink, "\naeslink=\n%s old=\n%s" % (aeslink, zaero.aeslink[key])
        zaero.aeslink[key] = aeslink
        model._type_to_id_map[aeslink.type].append(key)

    def add_mloads_object(self, mloads: MLOADS) -> None:
        """adds an MLOADS object"""
        key = mloads.mloads_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.mloads, "\nmloads=\n%s old=\n%s" % (mloads, zaero.mloads[key])
        zaero.mloads[key] = mloads
        model._type_to_id_map[mloads.type].append(key)

    def add_gloads_object(self, gloads: GLOADS) -> None:
        """adds an GLOADS object"""
        key = gloads.gloads_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.gloads, "\ngloads=\n%s old=\n%s" % (gloads, model.gloads[key])
        zaero.gloads[key] = gloads
        model._type_to_id_map[gloads.type].append(key)

    def add_dgust_object(self, dgust: DGUST) -> None:
        """adds an DGUST object"""
        key = dgust.dgust_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.dgust, "\ndgust=\n%s old=\n%s" % (dgust, zaero.dgust[key])
        zaero.dgust[key] = dgust
        model._type_to_id_map[dgust.type].append(key)

    def add_cgust_object(self, cgust: CGUST) -> None:
        """adds an CGUST object"""
        key = cgust.cgust_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.cgust, "\ncgust=\n%s old=\n%s" % (cgust, zaero.cgust[key])
        zaero.cgust[key] = cgust
        model._type_to_id_map[cgust.type].append(key)

    def add_ase_object(self, ase: ASE) -> None:
        """adds an ASE object"""
        key = ase.ase_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.ase, "\nase=\n%s old=\n%s" % (ase, zaero.ase[key])
        zaero.ase[key] = ase
        model._type_to_id_map[ase.type].append(key)

    def add_asecont_object(self, asecont: ASECONT) -> None:
        """adds an ASECONT object"""
        key = asecont.asecont_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asecont, "\nasecont=\n%s old=\n%s" % (asecont, zaero.asecont[key])
        zaero.asecont[key] = asecont
        model._type_to_id_map[asecont.type].append(key)

    def add_asegain_object(self, asegain: ASEGAIN) -> None:
        """adds an ASEGAIN object"""
        key = asegain.asegain_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asegain, "\nasegain=\n%s old=\n%s" % (asegain, zaero.asegain[key])
        zaero.asegain[key] = asegain
        model._type_to_id_map[asegain.type].append(key)

    def add_asesnsr_object(self, asesnsr: ASESNSR) -> None:
        """adds an ASESNSR object"""
        key = asesnsr.asesnsr_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asesnsr, "\nasesnsr=\n%s old=\n%s" % (asesnsr, zaero.asesnsr[key])
        zaero.asesnsr[key] = asesnsr
        model._type_to_id_map[asesnsr.type].append(key)

    def add_asesns1_object(self, asesns1: ASESNS1) -> None:
        """adds an ASESNS1 object"""
        key = asesns1.asesns1_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asesns1, "\nasesns1=\n%s old=\n%s" % (asesns1, zaero.asesns1[key])
        zaero.asesns1[key] = asesns1
        model._type_to_id_map[asesns1.type].append(key)

    def add_gainset_object(self, gainset: GAINSET) -> None:
        """adds an GAINSET object"""
        key = gainset.gainset_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.gainset, "\ngainset=\n%s old=\n%s" % (gainset, zaero.gainset[key])
        zaero.gainset[key] = gainset
        model._type_to_id_map[gainset.type].append(key)

    def add_pltbode_object(self, pltbode: PLTBODE) -> None:
        """adds an PLTBODE object"""
        key = pltbode.set_id
        model = self.model
        zaero = model.zaero
        if key in zaero.pltbode:
            model.log.warning(f"skipping duplicate PLTBODE\n{str(pltbode)}")
            return
        assert key not in zaero.pltbode, "\npltbode=\n%s old=\n%s" % (pltbode, zaero.pltbode[key])
        zaero.pltbode[key] = pltbode
        model._type_to_id_map[pltbode.type].append(key)

    def add_cjunct_object(self, cjunct: CJUNCT) -> None:
        """adds an CJUNCT object"""
        key = cjunct.cjunct_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.cjunct, "\ncjunct=\n%s old=\n%s" % (cjunct, zaero.cjunct[key])
        zaero.cjunct[key] = cjunct
        model._type_to_id_map[cjunct.type].append(key)

    def add_conct_object(self, conct: CONCT) -> None:
        """adds an CONCT object"""
        key = conct.conct_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.conct, "\nconct=\n%s old=\n%s" % (conct, zaero.conct[key])
        zaero.conct[key] = conct
        model._type_to_id_map[conct.type].append(key)

    def add_aerolag_object(self, aerolag: AEROLAG) -> None:
        """adds an CONCT object"""
        key = aerolag.aerolag_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.conct, "\naerolag=\n%s old=\n%s" % (aerolag, zaero.aerolag[key])
        zaero.aerolag[key] = aerolag
        model._type_to_id_map[aerolag.type].append(key)

    def add_actu_object(self, actu: ACTU) -> None:
        """adds an AESURFZ object"""
        key = actu.actu_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.actu, "\nactu=\n%s old=\n%s" % (actu, zaero.actu[key])
        zaero.actu[key] = actu
        model._type_to_id_map[actu.type].append(key)

    def add_mkaeroz_object(self, mkaeroz: MKAEROZ) -> None:
        """adds an MKAEROZ object"""
        assert mkaeroz.sid not in self.model.zaero.mkaeroz
        assert mkaeroz.sid > 0
        key = mkaeroz.sid
        self.model.zaero.mkaeroz[key] = mkaeroz
        self.model._type_to_id_map[mkaeroz.type].append(key)

    def add_trimvar_object(self, trimvar: TRIMVAR) -> None:
        """adds an TRIMVAR object"""
        key = trimvar.var_id
        assert trimvar.var_id not in self.model.zaero.trimvar, "\ntrimvar=\n%s old=\n%s" % (
            trimvar,
            self.model.zaero.trimvar[key],
        )
        assert trimvar.var_id > 0
        self.model.zaero.trimvar[key] = trimvar
        self.model._type_to_id_map[trimvar.type].append(key)

    def add_trimlnk_object(self, trimlnk: TRIMLNK) -> None:
        """adds an TRIMLNK object"""
        assert trimlnk.link_id not in self.model.zaero.trimlnk
        assert trimlnk.link_id > 0
        key = trimlnk.link_id
        self.model.zaero.trimlnk[key] = trimlnk
        self.model._type_to_id_map[trimlnk.type].append(key)

    def add_attach_object(self, attach: ATTACH) -> None:
        """adds an ATTACH object"""
        assert attach.attach_id not in self.model.zaero.attach
        assert attach.attach_id > 0
        key = attach.attach_id
        self.model.zaero.attach[key] = attach
        self.model._type_to_id_map[attach.type].append(key)

    def add_pltmode_object(self, plot: PLTMODE) -> None:
        """adds an PLTMODE object"""
        assert plot.set_id not in self.model.zaero.pltmode, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zaero.pltmode[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_pltaero_object(self, plot: PLTAERO) -> None:
        """adds an PLTAERO object"""
        assert plot.set_id not in self.model.zaero.pltaero, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zaero.pltaero[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_pltvg_object(self, plot: PLTVG) -> None:
        """adds an PLTVG object"""
        assert plot.set_id not in self.model.zaero.pltvg, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zaero.pltvg[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_pltsurf_object(self, plot: PLTSURF) -> None:
        """adds an PLTSURF object"""
        assert plot.set_id not in self.model.zaero.pltsurf, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zaero.pltsurf[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_pltcp_object(self, plot: PLTCP) -> None:
        """adds an PLTCP object"""
        # assert plot.set_id not in self.model.zaero.pltcp, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        if key not in self.model.zaero.pltcp:
            self.model.zaero.pltcp[key] = []
        self.model.zaero.pltcp[key].append(plot)
        self.model._type_to_id_map[plot.type].append(key)

    def add_plttrim_object(self, plot: PLTTRIM) -> None:
        """adds an PLTTRIM object"""
        # assert plot.set_id not in self.model.zaero.plttrim, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        if key not in self.model.zaero.plttrim:
            self.model.zaero.plttrim[key] = []
        self.model.zaero.plttrim[key].append(plot)
        self.model._type_to_id_map[plot.type].append(key)

    def add_plttime_object(self, plot: PLTTIME) -> None:
        """adds an PLTTIME object"""
        # assert plot.set_id not in self.model.zaero.pltcp, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        if key not in self.model.zaero.plttime:
            self.model.zaero.plttime[key] = []
        self.model.zaero.plttime[key].append(plot)
        self.model._type_to_id_map[plot.type].append(key)

    def add_pltflut_object(self, plot: PLTFLUT) -> None:
        """adds an PLTFLUT object"""
        # assert plot.set_id not in self.model.zaero.pltcp, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        if key not in self.model.zaero.pltflut:
            self.model.zaero.pltflut[key] = []
        self.model.zaero.pltflut[key].append(plot)
        self.model._type_to_id_map[plot.type].append(key)

    def add_pltmist_object(self, plot: PLTMIST) -> None:
        """adds an PLTMIST object"""
        assert plot.set_id not in self.model.zaero.pltmist, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zaero.pltmist[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_flutter_table_object(
        self, flutter_table: FIXHATM | FIXMATM | FIXMDEN | FIXMACH
    ) -> None:
        """adds an FIXMATM object"""
        key = flutter_table.sid
        model = self.model
        zaero = model.zaero
        assert key not in zaero.flutter_table, "\nflutter_table=\n%s old=\n%s" % (
            flutter_table,
            zaero.flutter_table[key],
        )
        zaero.flutter_table[key] = flutter_table
        model._type_to_id_map[flutter_table.type].append(key)

    def add_atmos_object(self, atmos: ATMOS) -> None:
        """adds an ATMOS object"""
        key = atmos.atmos_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.atmos, "\natmos=\n%s old=\n%s" % (atmos, zaero.atm[key])
        zaero.atmos[key] = atmos
        model._type_to_id_map[atmos.type].append(key)


class ZAERO:
    def __init__(self, model):
        self.model = model
        self.caero_to_name_map = {}
        self._add_methods = AddMethods(model)

        # singletons
        self.splinem: Optional[SPLINEM] = None

        # aero models
        self.atmos: dict[int, ATMOS] = {}
        self.flutter_table: dict[int, FIXHATM | FIXMATM | FIXMDEN | FIXMACH] = {}

        #: store PANLST1,PANLST2,PANLST3
        self.pltsurf: dict[int, PLTSURF] = {}
        self.panlsts: dict[int, PANLST1 | PANLST2 | PANLST3] = {}
        self.attach: dict[int, PLTAERO] = {}
        self.pltaero: dict[int, PLTAERO] = {}
        self.pltmode: dict[int, PLTMODE] = {}
        self.pltmist: dict[int, PLTMIST] = {}

        #: store PAFOIL7/PAFOIL8
        self.pafoil: dict[int, PAFOIL7 | PAFOIL8] = {}
        self.mloads: dict[int, MLOADS] = {}

        # TODO: add me
        self.dse: dict[int, DSE] = {}  # add me
        self.extfile: dict[int, EXTFILE] = {}
        self.conmlst: dict[int, CONMLST] = {}
        self.cpfact: dict[int, CPFACT] = {}
        self.apconst: dict[int, APCONST] = {}
        self.spline0: dict[int, SPLINE0] = {}
        self.pbody7: dict[int, PBODY7] = {}
        self.trimflt: dict[int, TRIMFLT] = {}
        self.output4: dict[int, OUTPUT4] = {}
        self.crospsd: dict[int, CROSPSD] = {}
        self.foilsec: dict[int, FOILSEC] = {}
        self.gengust_card: dict[int, GENGUST] = {}
        self.mftgust: dict[int, MFTGUST] = {}
        self.dmis: dict[int, DMIS] = {}
        self.cellwng: dict[int, CELLWNG] = {}
        self.cellbox: dict[int, CELLBOX] = {}
        self.inpcfd: dict[int, INPCFD] = {}
        self.omitcfd: dict[int, OMITCFD] = {}
        self.wt1ajj: dict[int, WT1AJJ] = {}
        self.wt1frc: dict[int, WT1FRC] = {}
        self.wt2ajj: dict[int, WT2AJJ] = {}
        self.wtucp: dict[int, WTUCP] = {}
        self.trimobj_card: dict[int, TRIMOBJ] = {}
        self.trimcon_card: dict[int, TRIMCON] = {}
        # FOILSEC
        # TRIMFLT

        # transient
        self.eloads: dict[int, ELOADS] = {}  # add me
        self.mldcomd: dict[int, MLDCOMD] = {}
        self.plttime: dict[int, PLTTIME] = {}

        # trim
        self.pltcp: dict[int, PLTCP] = {}
        self.plttrim: dict[int, PLTTRIM] = {}
        self.aeslink: dict[int, AESLINK] = {}
        self.trimvar: dict[int, TRIMVAR] = {}
        self.trimlnk: dict[int, TRIMLNK] = {}
        self.mldtrim: dict[int, MLDTRIM] = {}
        self.trimfnc: dict[int, TRIMFNC] = {}
        self.trimobj: dict[int, TRIMOBJ] = {}  # add me
        self.trimcon: dict[int, TRIMCON] = {}  # add me
        self.trimadd: dict[int, TRIMADD] = {}  # add me
        self.trimflt: dict[int, TRIMFLT] = {}  # add me

        # flutter
        self.pltvg: dict[int, PLTVG] = {}
        self.mkaeroz: dict[int, MKAEROZ] = {}
        self.nlfltr: dict[int, NLFLTR] = {}  # add me
        self.pltflut: dict[int, PLTFLUT] = {}

        # gust
        self.gloads: dict[int, GLOADS] = {}
        self.dgust: dict[int, DGUST] = {}
        self.cgust: dict[int, CGUST] = {}
        self.gengust: dict[int, GENGUST] = {}  # add me
        self.gustinp: dict[int, GUSTINP] = {}  # add me
        self.dfs: dict[int, Any] = {}

        # ase
        self.ase: dict[int, ASE] = {}
        self.asecont: dict[int, ASECONT] = {}
        self.asesnsr: dict[int, ASESNSR] = {}
        self.asesns1: dict[int, ASESNS1] = {}
        self.asegain: dict[int, ASEGAIN] = {}
        self.gainset: dict[int, GAINSET] = {}
        self.pltbode: dict[int, PLTBODE] = {}
        self.aseout: dict[int, ASEOUT] = {}  # add me
        self.apcnsnd: dict[int, APCNSND] = {}  # add me
        self.apcnscp: dict[int, APCNSCP] = {}  # add me
        self.mimoss: dict[int, MIMOSS] = {}
        self.mimotf: dict[int, MIMOTF] = {}
        self.sisotf: dict[int, SISOTF] = {}
        self.cmargin: dict[int, CMARGIN] = {}
        self.aerolag: dict[int, AEROLAG] = {}
        self.sensr: dict[int, SENSR] = {}
        self.gain: dict[int, GAIN] = {}
        self.sumblk: dict[int, SUMBLK] = {}
        self.deadbn: dict[int, DEADBN] = {}
        self.delay_zaero: dict[int, DELAY_ZAERO] = {}
        self.filtfl: dict[int, FILTFL] = {}
        self.limtr: dict[int, LIMTR] = {}

        # other
        self.extinp: dict[int, EXTINP] = {}
        self.extout: dict[int, EXTOUT] = {}
        self.loadmod: dict[int, LOADMOD] = {}
        self.rbred: dict[int, RBRED] = {}

        self.mldstat: dict[int, MLDSTAT] = {}
        self.minstat: dict[int, MINSTAT] = {}
        self.mldprnt: dict[int, MLDPRNT] = {}
        self.mldtime: dict[int, MLDTIME] = {}
        self.dmil: dict[tuple[str, int, int], DMIL] = {}
        self.actu: dict[int, ACTU] = {}
        self.cjunct: dict[int, CJUNCT] = {}
        self.conct: dict[int, CONCT] = {}
        self.tfset: dict[int, TFSET] = {}
        self.senset: dict[int, SENSET] = {}
        self.cnctset: dict[int, CNCTSET] = {}
        self.surfset: dict[int, SURFSET] = {}
        self.setadd: dict[int, SETADD] = {}
        self.thermo: dict[int, THERMO] = {}  # add me

        # TODO:
        self.extfile: dict[int, EXTFILE] = {}

    @classmethod
    def _init_from_self(cls, model: BDF):
        """helper method for dict_to_h5py"""
        return cls(model)

    def clear(self):
        """clears out the ZAERO object"""
        self.panlsts = {}
        self.mkaeroz = {}
        self.trimvar = {}
        self.trimlnk = {}
        self.pafoil = {}
        # self.aeroz = {}

    def object_attributes(
        self,
        mode: str = "public",
        keys_to_skip: Optional[list[str]] = None,
        filter_properties: bool = False,
    ):
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : list[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong
        """
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            "log",
            "model",
        ]
        return object_attributes(
            self,
            mode=mode,
            keys_to_skip=keys_to_skip + my_keys_to_skip,
            filter_properties=filter_properties,
        )

    def object_methods(
        self, mode: str = "public", keys_to_skip: Optional[list[str]] = None
    ) -> list[str]:
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : list[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong
        """
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip: list[str] = []

        my_keys_to_skip = [
            "log",
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip + my_keys_to_skip)

    def add_segmesh(self,
                    segmesh_id: int,
                    naxial: int,
                    nradial: int,
                    nose_radius: float,
                    iaxis: int,
                    itypes: list[int],
                    xs: list[float],
                    cambers: list[float],
                    ys: list[float],
                    zs: list[float],
                    idys: list[int],
                    idzs: list[int],
                    comment: str = "",
                    ) -> SEGMESH:
        card = SEGMESH(
            segmesh_id, naxial, nradial, nose_radius,
            iaxis, itypes, xs, cambers,
            ys, zs,
            idys, idzs,
            comment=comment,)
        self.model._add_methods.add_paero_object(card)
        return card

    def add_atmos(
        self,
        atm_id: int,
        mass_unit: str,
        length_unit: str,
        temperature_unit: str,
        atmosphere_table: list[float],
        comment: str = "") -> ATMOS:
        # alt: np.ndarray, sos: np.ndarray,
        # density: np.ndarray, temperature: np.ndarray) -> ATMOS:
        card = ATMOS(
            atm_id, mass_unit, length_unit, temperature_unit, atmosphere_table, comment=comment
        )
        self._add_methods.add_atmos_object(card)
        return card

    # ---------------------------------------------------------------
    # atmosphere / flutter table
    # ---------------------------------------------------------------
    def add_fixhatm(
        self,
        sid: int,
        alt: float,
        atm_id: int,
        mass_unit: str,
        length_unit: str,
        mkaeroz_ids: list[int],
        vref: float=1.0,
        fluttf_id: int=0,
        print_flag: int=0,
        comment: str = "",) -> FIXHATM:
        card = FIXHATM(
            sid, alt, atm_id, mass_unit, length_unit, mkaeroz_ids,
            fluttf_id=fluttf_id, print_flag=print_flag,
            vref=vref, comment=comment,
        )
        self._add_methods.add_flutter_table_object(card)
        return card

    def add_fixmatm(
        self,
        sid: int,
        mkaeroz_id: int,
        mass_unit: str,
        length_unit: str,
        alts: list[float],
        atm_id: int = 0, fluttf_id: int=0,
        print_flag: int=0, vref: float = 1.0, comment: str = "",) -> FIXMATM:
        card = FIXMATM(
            sid, mkaeroz_id,  mass_unit, length_unit, alts,
            atm_id=atm_id,
            fluttf_id=fluttf_id, print_flag=print_flag,
            vref=vref, comment=comment,)
        self._add_methods.add_flutter_table_object(card)
        return card

    def add_fixmach(
        self,
        sid: int,
        mkaeroz_id: int,
        mass_unit: str,
        length_unit: str,
        velocity: list[float],
        rho: list[float],
        fluttf_id: int = 0,
        print_flag: int = 0,
        vref: float = 1.0,
        comment: str = "",) -> FIXMACH:
        card = FIXMACH(
            sid, mkaeroz_id,
            mass_unit, length_unit,
            velocity, rho,
            fluttf_id=fluttf_id,
            print_flag=print_flag,
            vref=vref,
            comment=comment,
        )
        self._add_methods.add_flutter_table_object(card)
        return card

    def add_fixmden(
        self,
        sid: int,
        mkaeroz_id: int,
        rho: float,
        mass_unit: str,
        length_unit: str,
        velocity: list[float],
        fluttf_id: int = 0,
        print_flag: int = 0,
        vref: float = 1.0,
        comment: str = "",
    ) -> FIXMDEN:
        card = FIXMDEN(
            sid,
            mkaeroz_id,
            rho,
            mass_unit, length_unit,
            velocity,
            fluttf_id=fluttf_id,
            print_flag=print_flag,
            vref=vref,
            comment=comment,
        )
        self._add_methods.add_flutter_table_object(card)
        return card

    # ---------------------------------------------------------------
    # geometry
    # ---------------------------------------------------------------
    def add_splinem(self, save_flag: str, filename, comment: str = "") -> SPLINEM:
        card = SPLINEM(save_flag, filename, comment=comment)
        self._add_methods.add_splinem_object(card)
        return card

    def add_panlst1(
        self, eid: int, macro_id: int, box1: int, box2: int, comment: str = ""
    ) -> PANLST1:
        card = PANLST1(eid, macro_id, box1, box2, comment=comment)
        self._add_methods.add_panlst_object(card)
        return card

    def add_panlst2(self, eid: int, macro_id: int, boxes: list[int], comment: str = "") -> PANLST2:
        card = PANLST2(eid, macro_id, boxes, comment=comment)
        self._add_methods.add_panlst_object(card)
        return card

    def add_panlst3(self, eid: int, panel_groups, comment: str = "") -> PANLST3:
        card = PANLST3(eid, panel_groups, comment=comment)
        self._add_methods.add_panlst_object(card)
        return card

    def add_pafoil7(
        self,
        pid: int,
        i_axial: int,
        i_thickness_root: int,
        i_camber_root: int,
        le_radius_root: float,
        i_thickness_tip: int,
        i_camber_tip: int,
        le_radius_tip: float,
        comment: str = "",
    ) -> PAFOIL7:
        card = PAFOIL7(
            pid,
            i_axial,
            i_thickness_root,
            i_camber_root,
            le_radius_root,
            i_thickness_tip,
            i_camber_tip,
            le_radius_tip,
            comment=comment,
        )
        self._add_methods.add_pafoil_object(card)
        return card

    def add_pafoil8(
        self,
        pid: int,
        i_axial: int,
        i_thickness_root: int,
        i_camber_root: int,
        le_radius_root: float,
        i_thickness_tip: int,
        comment: str = "",
    ) -> PAFOIL8:
        card = PAFOIL8(
            pid,
            i_axial,
            i_thickness_root,
            i_camber_root,
            le_radius_root,
            i_thickness_tip,
            comment=comment,
        )
        self._add_methods.add_pafoil_object(card)
        return card

    def add_aesurfz(
        self,
        label: str,
        surface_type: str,
        cid: int,
        panlst: int,
        setg: int,
        actuator_tf: int,
        comment: str = "",
    ) -> AESURFZ:
        card = AESURFZ(label, surface_type, cid, panlst, setg, actuator_tf, comment=comment)
        self._add_methods.add_aesurfz_object(card)
        return card

    def add_aeslink(
        self,
        label,
        link_type: str,
        actu_id: int,
        independent_labels: list[str],
        linking_coefficients: list[float],
        comment: str = "",
    ) -> AESLINK:
        card = AESLINK(
            label, link_type, actu_id, independent_labels, linking_coefficients, comment=comment
        )
        self._add_methods.add_aeslink_object(card)
        return card

    def add_attach(
        self, attach_id: int, model_name: str, setk: int, refgrid: int, comment: str = ""
    ) -> ATTACH:
        card = ATTACH(attach_id, model_name, setk, refgrid, comment=comment)
        self._add_methods.add_attach_object(card)
        return card

    # ---------------------------------------------------------------
    # plot
    # ---------------------------------------------------------------
    def add_pltmode(
        self,
        set_id: int,
        symmetry: str,
        mode: int,
        output_format: str,
        filename: str,
        max_disp: float=1.0,
        comment: str = "",
    ) -> PLTMODE:
        card = PLTMODE(set_id, symmetry, mode, output_format, filename,
                       max_disp=max_disp, comment=comment)
        self._add_methods.add_pltmode_object(card)
        return card

    def add_pltaero(
        self,
        set_id: int,
        femgrid: str,
        offset: int,
        out_format: str,
        filename: str,
        cell: str = "NO",
        vct: str = "NO",
        comment: str = "",
    ) -> PLTAERO:
        card = PLTAERO(
            set_id, femgrid, offset, out_format, filename, cell=cell, vct=vct, comment=comment
        )
        self._add_methods.add_pltaero_object(card)
        return card

    def add_pltvg(
        self,
        set_id: int,
        flutter_id: int,
        xaxis: str,
        filename,
        nmode: int = 0,
        output_format: str = "TABLE",
        rho_ref: float = 1.0,
        comment: str = "",
    ) -> PLTVG:
        card = PLTVG(
            set_id,
            flutter_id,
            xaxis,
            filename,
            nmode=nmode,
            output_format=output_format,
            rho_ref=rho_ref,
            comment=comment,
        )
        self._add_methods.add_pltvg_object(card)
        return card

    def add_pltcp(
        self,
        set_id: int,
        sym_flag: str,
        mkaeroz_id: int,
        ik: int,
        mode: int,
        out_format: str,
        filename,
        aero_filename="",
        comment: str = "",
    ) -> PLTCP:
        card = PLTCP(
            set_id,
            sym_flag,
            mkaeroz_id,
            ik,
            mode,
            out_format,
            filename,
            aero_filename=aero_filename,
            comment=comment,
        )
        self._add_methods.add_pltcp_object(card)
        return card

    def add_pltmist(
        self,
        set_id: int,
        ase_id: int,
        irow: int,
        icol: int,
        klist: int,
        out_format: str,
        filename: str,
        comment: str = "",
    ) -> PLTMIST:
        card = PLTMIST(set_id, ase_id, irow, icol, klist, out_format, filename, comment=comment)
        self._add_methods.add_pltmist_object(card)
        return card

    def add_pltsurf(
        self,
        set_id: int,
        label: str,
        out_format: str,
        filename: str,
        scale_factor: float = 1.0,
        comment: str = "",
    ) -> PLTSURF:
        card = PLTSURF(
            set_id, label, out_format, filename, scale_factor=scale_factor, comment=comment
        )
        self._add_methods.add_pltsurf_object(card)
        return card

    def add_pltflut(
        self,
        set_id: int,
        out_format: str,
        filename: str,
        scale_factor: float = 1.0,
        comment: str = "",
    ) -> PLTFLUT:
        card = PLTFLUT(set_id, out_format, filename, scale_factor=scale_factor, comment=comment)
        self._add_methods.add_pltflut_object(card)
        return card

    def add_pltbode(
        self,
        set_id: int,
        cmargin_id: int,
        fmin: float,
        fmax: float,
        nf: int,
        filename: str,
        log_scale_flag: int = 0,
        draw_flag: int = 0,
        comment: str = "",
    ) -> PLTBODE:
        card = PLTBODE(
            set_id,
            cmargin_id,
            fmin,
            fmax,
            nf,
            filename,
            log_scale_flag=log_scale_flag,
            draw_flag=draw_flag,
            comment=comment,
        )
        self._add_methods.add_pltbode_object(card)
        return card

    def add_plttime(
        self,
        set_id: int,
        mloads_id: int,
        tstart: float,
        tend: float,
        ndt: int,
        out_type: str,
        filename: str,
        aero_filename: str,
        output_format: str = "TECPLOT",
        scale_factor: float = 1.0,
        comment: str = "",
    ) -> PLTTIME:
        card = PLTTIME(
            set_id,
            mloads_id,
            tstart,
            tend,
            ndt,
            out_type,
            filename,
            aero_filename,
            output_format=output_format,
            scale_factor=scale_factor,
            comment=comment,
        )
        self._add_methods.add_plttime_object(card)
        return card

    def add_plttrim(
        self,
        set_id: int,
        trim_id: int,
        out_type: str,
        filename: str,
        aero_filename: str,
        flex: str = "FLEX",
        output_format: str = "TECPLOT",
        scale_factor: float = 1.0,
        comment: str = "",
    ) -> PLTTRIM:
        card = PLTTRIM(
            set_id,
            trim_id,
            out_type,
            filename,
            aero_filename,
            flex=flex,
            output_format=output_format,
            scale_factor=scale_factor,
            comment=comment,
        )
        self._add_methods.add_plttrim_object(card)
        return card

    # ---------------------------------------------------------------
    # flutter
    # ---------------------------------------------------------------
    def add_mkaeroz(
        self,
        sid: int,
        mach: float,
        freqs: list[float],
        method: int = 0,
        flt_id: int = 0,
        print_flag: int = 0,
        filename: str = '',
        save: str = "SAVE",
        comment: str = "",
    ) -> MKAEROZ:
        card = MKAEROZ(
            sid,
            mach,
            freqs,
            method=method,
            flt_id=flt_id,
            filename=filename,
            print_flag=print_flag,
            save=save,
            comment=comment,
        )
        self._add_methods.add_mkaeroz_object(card)
        return card

    # ---------------------------------------------------------------
    # trim
    # ---------------------------------------------------------------
    def add_trimvar(
        self,
        var_id: int,
        label: str,
        lower: float,
        upper: float,
        trimlnk_id: int,
        dmi,
        sym: int,
        initial: Optional[float] = None,
        dcd="NONE",
        dcy="NONE",
        dcl="NONE",
        dcr="NONE",
        dcm="NONE",
        dcn="NONE",
        comment: str = "",
    ) -> TRIMVAR:
        card = TRIMVAR(
            var_id,
            label,
            lower,
            upper,
            trimlnk_id,
            dmi,
            sym,
            initial=initial,
            dcd=dcd,
            dcy=dcy,
            dcl=dcl,
            dcr=dcr,
            dcm=dcm,
            dcn=dcn,
            comment=comment,
        )
        self._add_methods.add_trimvar_object(card)
        return card

    def add_trimlnk(
        self, link_id: int, sym: str, coeffs: list[float], var_ids: list[int], comment: str = ""
    ) -> TRIMLNK:
        card = TRIMLNK(link_id, sym, coeffs, var_ids, comment=comment)
        self._add_methods.add_trimlnk_object(card)
        return card

    def add_trimfnc(
        self,
        trimfnc_id: int,
        fcn_type: str,
        label: str,
        rhs_flag: str,
        is_set,
        ia_set,
        remark: str,
        comment: str = "",
    ) -> TRIMFNC:
        card = TRIMFNC(
            trimfnc_id, fcn_type, label, rhs_flag, is_set, ia_set, remark, comment=comment
        )
        self._add_methods.add_trimfnc_object(card)
        return card

    # ---------------------------------------------------------------
    # maneuver loads
    # ---------------------------------------------------------------
    def add_mloads(
        self,
        mloads_id: int,
        asecont_id: int,
        flutter_id: int,
        minstat_id: int,
        mldstat_id: int,
        mldcomd_id: int,
        mldtime_id: int,
        mldprnt_id: int,
        fmax: float,
        save_freq: str,
        df: float = 0.01,
        filename: str = "",
        comment: str = "",
    ) -> MLOADS:
        card = MLOADS(
            mloads_id,
            asecont_id,
            flutter_id,
            minstat_id,
            mldstat_id,
            mldcomd_id,
            mldtime_id,
            mldprnt_id,
            fmax,
            save_freq,
            df=df,
            filename=filename,
            comment=comment,
        )
        self._add_methods.add_mloads_object(card)
        return card

    def add_extinp(
        self,
        extinp_id: int,
        input_type: int,
        itf_id: int,
        itf_component: int,
        label: str,
        comment: str = "",
    ) -> EXTINP:
        card = EXTINP(extinp_id, input_type, itf_id, itf_component, label, comment=comment)
        self._add_methods.add_extinp_object(card)
        return card

    def add_extout(
        self,
        extout_id: int,
        input_type: int,
        itf_id: int,
        itf_component: int,
        label: str,
        comment: str = "",
    ) -> EXTOUT:
        card = EXTOUT(extout_id, input_type, itf_id, itf_component, label, comment=comment)
        self._add_methods.add_extout_object(card)
        return card

    def add_actu(self, actu_id: int, a0: float, a1: float, a2: float, comment: str = "") -> ACTU:
        card = ACTU(actu_id, a0, a1, a2, comment=comment)
        self._add_methods.add_actu_object(card)
        return card

    def add_loadmod(
        self, loadmod_id: int, label, cid: int, set_k: int, set_g: int, comment: str = ""
    ) -> LOADMOD:
        card = LOADMOD(loadmod_id, label, cid, set_k, set_g, comment=comment)
        self._add_methods.add_loadmod_object(card)
        return card

    def add_rbred(
        self, sid: int, id_ase: int, component: str, node_id: int, phugoid0: str, comment: str = ""
    ) -> RBRED:
        card = RBRED(sid, id_ase, component, node_id, phugoid0, comment=comment)
        self._add_methods.add_rbred_object(card)
        return card

    # ---------------------------------------------------------------
    # gust
    # ---------------------------------------------------------------
    def add_gloads(
        self,
        gloads_id: int,
        asecont_id: int,
        flutter_id: int,
        minstat_id: int,
        mldstat_id: int,
        mldcomd_id: int,
        mldtime_id: int,
        mldprnt_id: int,
        save_flag: str,
        form: str,
        filename: str,
        save_freq: str,
        filename_freq: str,
        comment: str = "",
    ) -> GLOADS:
        card = GLOADS(
            gloads_id,
            asecont_id,
            flutter_id,
            minstat_id,
            mldstat_id,
            mldcomd_id,
            mldtime_id,
            mldprnt_id,
            save_flag,
            form,
            filename,
            save_freq,
            filename_freq,
            comment=comment,
        )
        self._add_methods.add_gloads_object(card)
        return card

    def add_dgust(
        self,
        dgust_id: int,
        gust_type: str,
        length_gust,
        gust_velocity: float,
        x0: float,
        fmax: float = 0.0,
        df: float = 0.01,
        nap: int = 0,
        comment: str = "",
    ) -> DGUST:
        card = DGUST(
            dgust_id,
            gust_type,
            length_gust,
            gust_velocity,
            x0,
            fmax=fmax,
            df=df,
            nap=nap,
            comment=comment,
        )
        self._add_methods.add_dgust_object(card)
        return card

    def add_cgust(
        self,
        cgust_id: int,
        gust_type: str,
        one_over_velocity: float,
        x0: float,
        rms_velocity: float,
        length_gust: float = 2500.0,
        fmax: float = 0.0,
        df: float = 0.01,
        comment: str = "",
    ) -> CGUST:
        card = CGUST(
            cgust_id,
            gust_type,
            one_over_velocity,
            x0,
            rms_velocity,
            length_gust=length_gust,
            fmax=fmax,
            df=df,
            comment=comment,
        )
        self._add_methods.add_cgust_object(card)
        return card

    # ---------------------------------------------------------------
    # ASE (aeroservoelastic)
    # ---------------------------------------------------------------
    def add_ase(
        self,
        ase_id: int,
        asecont_id: int,
        flutter_id: int,
        mldstat_id: int,
        minstat_id: int,
        cmargin_id: int,
        comment: str = "",
    ) -> ASE:
        card = ASE(
            ase_id, asecont_id, flutter_id, mldstat_id, minstat_id, cmargin_id, comment=comment
        )
        self._add_methods.add_ase_object(card)
        return card

    def add_asecont(
        self,
        asecont_id: int,
        surf_id: int,
        sens_id: int,
        tf_id: int,
        gain_id: int,
        conct_id: int,
        extinp_set_id: int,
        extout_set_id: int,
        comment: str = "",
    ) -> ASECONT:
        card = ASECONT(
            asecont_id,
            surf_id,
            sens_id,
            tf_id,
            gain_id,
            conct_id,
            extinp_set_id,
            extout_set_id,
            comment=comment,
        )
        self._add_methods.add_asecont_object(card)
        return card

    def add_asesnsr(
        self,
        asesnsr_id: int,
        sensor_type: int,
        sgid: int,
        component: int,
        factor: float,
        sum_method: str,
        comment: str = "",
    ) -> ASESNSR:
        card = ASESNSR(
            asesnsr_id, sensor_type, sgid, component, factor, sum_method, comment=comment
        )
        self._add_methods.add_asesnsr_object(card)
        return card

    def add_asesns1(
        self,
        asesns1_id: int,
        label: str,
        ikey,
        factor: float,
        sum_method: str = "NO",
        comment: str = "",
    ) -> ASESNS1:
        card = ASESNS1(asesns1_id, label, ikey, factor, sum_method=sum_method, comment=comment)
        self._add_methods.add_asesns1_object(card)
        return card

    def add_asegain(
        self,
        asegain_id: int,
        otf_id: int,
        c_out: int,
        itf_id: int,
        c_in: int,
        gain: float,
        gain_type: str,
        comment: str = "",
    ) -> ASEGAIN:
        card = ASEGAIN(asegain_id, otf_id, c_out, itf_id, c_in, gain, gain_type, comment=comment)
        self._add_methods.add_asegain_object(card)
        return card

    def add_gainset(self, gainset_id: int, ids: list[int], comment: str = "") -> GAINSET:
        card = GAINSET(gainset_id, ids, comment=comment)
        self._add_methods.add_gainset_object(card)
        return card

    def add_cjunct(
        self, cjunct_id: int, nu: int, ny: int, values: list[float], comment: str = ""
    ) -> CJUNCT:
        card = CJUNCT(cjunct_id, nu, ny, values, comment=comment)
        self._add_methods.add_cjunct_object(card)
        return card

    def add_conct(
        self,
        conct_id: int,
        output_tf_id: int,
        output_component: int,
        input_tf_id: int,
        input_component: int,
        comment: str = "",
    ) -> CONCT:
        card = CONCT(
            conct_id, output_tf_id, output_component, input_tf_id, input_component, comment=comment
        )
        self._add_methods.add_conct_object(card)
        return card

    def add_tfset(self, tfset_id: int, ids: list[int], comment: str = "") -> TFSET:
        card = TFSET(tfset_id, ids, comment=comment)
        self._add_methods.add_tfset_object(card)
        return card

    def add_mimoss(
        self,
        mimoss_id: int,
        ntf: int,
        nu: int,
        ny: int,
        dmi_label: str,
        mimoss_type: str,
        print_flag: int,
        values: list[float],
        labels: list[str],
        comment: str = "",
    ) -> MIMOSS:
        card = MIMOSS(
            mimoss_id,
            ntf,
            nu,
            ny,
            dmi_label,
            mimoss_type,
            print_flag,
            values,
            labels,
            comment=comment,
        )
        self._add_methods.add_mimoss_object(card)
        return card

    def add_mimotf(
        self, mimotf_id: int, n_input: int, n_output: int, tf_ids: list[int], comment: str = ""
    ) -> MIMOTF:
        card = MIMOTF(mimotf_id, n_input, n_output, tf_ids, comment=comment)
        self._add_methods.add_mimotf_object(card)
        return card

    def add_sisotf(
        self,
        sisotf_id: int,
        nnumerator: int,
        ndenominator: int,
        b: list[float],
        a: list[float],
        comment: str = "",
    ) -> SISOTF:
        card = SISOTF(sisotf_id, nnumerator, ndenominator, b, a, comment=comment)
        self._add_methods.add_sisotf_object(card)
        return card

    def add_aerolag(
        self, aerolag_id: int, nlag: int, lag_values: list[float], comment: str = ""
    ) -> AEROLAG:
        card = AEROLAG(aerolag_id, nlag, lag_values, comment=comment)
        self._add_methods.add_aerolag_object(card)
        return card

    def add_sensr(
        self,
        sensr_id: int,
        sensor_type: int,
        sgid: int,
        component: int,
        factor: float,
        sum_method: str,
        comment: str = "",
    ) -> SENSR:
        card = SENSR(sensr_id, sensor_type, sgid, component, factor, sum_method, comment=comment)
        self._add_methods.add_sensr_object(card)
        return card

    def add_gain(self, gain_id: int, k: float, comment: str = "") -> GAIN:
        card = GAIN(gain_id, k, comment=comment)
        self._add_methods.add_gain_object(card)
        return card

    def add_sumblk(
        self, sumblk_id: int, nsignal: int, signs: list[float], comment: str = ""
    ) -> SUMBLK:
        card = SUMBLK(sumblk_id, nsignal, signs, comment=comment)
        self._add_methods.add_sumblk_object(card)
        return card

    def add_deadbn(self, deadbn_id: int, threshold: float, comment: str = "") -> DEADBN:
        card = DEADBN(deadbn_id, threshold, comment=comment)
        self._add_methods.add_deadbn_object(card)
        return card

    def add_delay_zaero(
        self, delay_id: int, tau: float, order: int, comment: str = ""
    ) -> DELAY_ZAERO:
        card = DELAY_ZAERO(delay_id, tau, order, comment=comment)
        self._add_methods.add_delay_zaero_object(card)
        return card

    def add_filtfl(
        self,
        filtfl_id: int,
        filter_type: int,
        freq: float,
        order: int,
        zeta: float,
        comment: str = "",
    ) -> FILTFL:
        card = FILTFL(filtfl_id, filter_type, freq, order, zeta, comment=comment)
        self._add_methods.add_filtfl_object(card)
        return card

    def add_limtr(self, limtr_id: int, lower: float, upper: float, comment: str = "") -> LIMTR:
        card = LIMTR(limtr_id, lower, upper, comment=comment)
        self._add_methods.add_limtr_object(card)
        return card

    # ---------------------------------------------------------------
    # sets
    # ---------------------------------------------------------------
    def add_senset(self, senset_id: int, ids: list[int], comment: str = "") -> SENSET:
        card = SENSET(senset_id, ids, comment=comment)
        self._add_methods.add_senset_object(card)
        return card

    def add_surfset(self, surfset_id: int, ids: list[int], comment: str = "") -> SURFSET:
        card = SURFSET(surfset_id, ids, comment=comment)
        self._add_methods.add_surfset_object(card)
        return card

    def add_cnctset(self, cnctset_id: int, ids: list[int], comment: str = "") -> CNCTSET:
        card = CNCTSET(cnctset_id, ids, comment=comment)
        self._add_methods.add_cnctset_object(card)
        return card

    def add_setadd(self, setadd_id: int, ids: list[int], comment: str = "") -> SETADD:
        card = SETADD(setadd_id, ids, comment=comment)
        self._add_methods.add_setadd_object(card)
        return card

    # ---------------------------------------------------------------
    # MLD / transient
    # ---------------------------------------------------------------
    def add_mldprnt(
        self,
        mldprnt_id: int,
        filename: str,
        form: str,
        tspnt: float,
        tepnt: float,
        labels: list[str],
        ikeys: list[int],
        psd_time: str = "TIME",
        sof="NO",
        comment: str = "",
    ) -> MLDPRNT:
        card = MLDPRNT(
            mldprnt_id,
            filename,
            form,
            tspnt,
            tepnt,
            labels,
            ikeys,
            psd_time=psd_time,
            sof=sof,
            comment=comment,
        )
        self._add_methods.add_mldprnt_object(card)
        return card

    def add_mldstat(
        self,
        mldstat_id: int,
        mldtrim_id: int,
        transform,
        filename: str,
        states: list[str],
        values: list[float],
        dx_tox: str = "YES",
        state_space_arr: str = "",
        state_space_brr: str = "",
        comment: str = "",
    ) -> MLDSTAT:
        card = MLDSTAT(
            mldstat_id,
            mldtrim_id,
            transform,
            filename,
            states,
            values,
            dx_tox=dx_tox,
            state_space_arr=state_space_arr,
            state_space_brr=state_space_brr,
            comment=comment,
        )
        self._add_methods.add_mldstat_object(card)
        return card

    def add_minstat(
        self,
        minstat_id,
        aerolag_id,
        itmax,
        apcid,
        pweight_id,
        dinit_id,
        klist_id,
        min_inp,
        print_flag,
        save_flag,
        filename,
        msmod,
        aerolag_gust_id=0,
        dinit_gust_id=0,
        comment: str = "",
    ) -> MINSTAT:
        card = MINSTAT(
            minstat_id,
            aerolag_id,
            itmax,
            apcid,
            pweight_id,
            dinit_id,
            klist_id,
            min_inp,
            print_flag,
            save_flag,
            filename,
            msmod,
            aerolag_gust_id=aerolag_gust_id,
            dinit_gust_id=dinit_gust_id,
            comment=comment,
        )
        self._add_methods.add_minstat_object(card)
        return card

    def add_mldtrim(
        self,
        mldtrim_id,
        gravity: float,
        nz: float,
        thkcam: str,
        trim_vars: list[str],
        values: list[float],
        modal_dmi: str = "",
        comment: str = "",
    ) -> MLDTRIM:
        card = MLDTRIM(
            mldtrim_id, gravity, nz, thkcam, trim_vars, values, modal_dmi=modal_dmi, comment=comment
        )
        self._add_methods.add_mldtrim_object(card)
        return card

    def add_mldcomd(
        self, mldcomd_id: int, extinp_ids: list[int], table_ids: list[int], comment: str = ""
    ) -> MLDCOMD:
        card = MLDCOMD(mldcomd_id, extinp_ids, table_ids, comment=comment)
        self._add_methods.add_mldcomd_object(card)
        return card

    def add_mldtime(
        self,
        mldtime_id: int,
        tstart: float,
        tend: float,
        dt: float,
        out_dt: int,
        print_flag: int,
        method: int | None,
        comment: str = "",
    ) -> MLDTIME:
        card = MLDTIME(mldtime_id, tstart, tend, dt, out_dt, print_flag, method, comment=comment)
        self._add_methods.add_mldtime_object(card)
        return card

    # ---------------------------------------------------------------
    # other
    # ---------------------------------------------------------------
    def add_extfile(self, extfile_id: int, filename: str, comment: str = "") -> EXTFILE:
        card = EXTFILE(extfile_id, filename, comment=comment)
        self._add_methods.add_extfile_object(card)
        return card

    def add_dmil(
        self, name: str, col: int, row: int, values: list[float], comment: str = ""
    ) -> DMIL:
        card = DMIL(name, col, row, values, comment=comment)
        self._add_methods.add_dmil_object(card)
        return card

    def add_conmlst(
        self, conmlst_id: int, factors: list[float], conm_ids: list[int], comment: str = ""
    ) -> CONMLST:
        card = CONMLST(conmlst_id, factors, conm_ids, comment=comment)
        self._add_methods.add_conmlst_object(card)
        return card

    def add_cpfact(
        self,
        cpfact_id: int,
        idmk: int,
        sym: str,
        comp: str,
        cptype: str,
        panlst: str,
        factor1: float,
        factor2: float,
        strips: list[int],
        comment: str = "",
    ) -> CPFACT:
        card = CPFACT(
            cpfact_id, idmk, sym, comp, cptype, panlst, factor1, factor2, strips, comment=comment
        )
        self._add_methods.add_cpfact_object(card)
        return card

    def add_apconst(
        self,
        sid: int,
        da0: int,
        da1: int,
        da2: int,
        nrp: int,
        ncp: int,
        fr_values: list[float],
        fc_values: list[float],
        comment: str = "",
    ) -> APCONST:
        card = APCONST(sid, da0, da1, da2, nrp, ncp, fr_values, fc_values, comment=comment)
        self._add_methods.add_apconst_object(card)
        return card

    def add_spline0(
        self, eid: int, model_name: str, cp: int, setk: int, comment: str = ""
    ) -> SPLINE0:
        card = SPLINE0(eid, model_name, cp, setk, comment=comment)
        self._add_methods.add_spline0_object(card)
        return card

    def add_body7(self, eid: int, label: str, pid: int,
        nseg: int, idmeshes: list[int], acoord: int=0, comment: str="",) -> BODY7:
        card = BODY7(eid, label, pid, nseg, idmeshes, acoord=acoord, comment=comment)
        self.model._add_methods.add_caero_object(card)
        return card

    def add_pbody7(self, pid: int, wake: int, inlet: int,
                 idps: list[int], flowrates: list[float],
                 cp_base: float=-0.2, xs_wake: float=1.3,
                 xd_wake: float=1.1, yoffset: float=0.0, zoffset: float=0.0,
                 comment: str = "") -> PBODY7:
        card = PBODY7(pid, wake, inlet, idps, flowrates, cp_base=cp_base,
                      xs_wake=xs_wake, xd_wake=xd_wake, yoffset=yoffset,
                      zoffset=zoffset, comment=comment)
        self._add_methods.add_pbody7_object(card)
        return card

    def add_trimflt(
        self, trimflt_id: int, title: str, alpha: float, fields: list, comment: str = ""
    ) -> TRIMFLT:
        card = TRIMFLT(trimflt_id, title, alpha, fields, comment=comment)
        self._add_methods.add_trimflt_object(card)
        return card

    def add_cmargin(
        self,
        cmargin_id: int,
        gm_high: float,
        gm_low: float,
        pm_high: float,
        pm_low: float,
        df: float = 1e-4,
        print_flag: int = 0,
        nroot: int = 0,
        comment: str = "",
    ) -> CMARGIN:
        card = CMARGIN(
            cmargin_id,
            gm_high,
            gm_low,
            pm_high,
            pm_low,
            df=df,
            print_flag=print_flag,
            nroot=nroot,
            comment=comment,
        )
        self._add_methods.add_cmargin_object(card)
        return card

    def verify(self, xref):
        if self.model.nastran_format not in {"zona", "zaero"}:
            return
        for panlsts in self.panlsts.values():
            for panlst in panlsts:
                panlst._verify(xref)
        for mkaeroz in self.mkaeroz.values():
            mkaeroz._verify(xref)
        for trimvar in self.trimvar.values():
            trimvar._verify(xref)
        for trimlnk in self.trimlnk.values():
            trimlnk._verify(xref)
        for pafoil in self.pafoil.values():
            pafoil._verify(xref)
        for attach in self.attach.values():
            attach._verify(xref)

    def validate(self):
        if self.model.nastran_format not in {"zona", "zaero"}:
            return
        for panlst in self.panlsts.values():
            for panlsti in panlst:
                panlsti.validate()

        singletons, dicts, dicts_list = get_dicts(self, "write")
        for items in dicts_list:
            for item in items.values():
                for itemi in item:
                    itemi.validate()

        for items in dicts:
            for item in items.values():
                if isinstance(item, list):
                    if len(item):
                        print(item)
                        asdf
                    continue
                # self.model.log.info(f'xref {item.type}')
                item.validate()

    def update_for_zaero(self):
        """updates for zaero"""
        card_parser = self.model._card_parser
        add_methods = self.model._add_methods
        zaero_add = self._add_methods
        card_parser2 = {
            # aero models
            "AEROZ": (AEROZ, add_methods.add_aeros_object),
            "ATMOS": (ATMOS, zaero_add.add_atmos_object),
            "FIXMATM": (FIXMATM, zaero_add.add_flutter_table_object),
            "FIXHATM": (FIXHATM, zaero_add.add_flutter_table_object),
            "FIXMACH": (FIXMACH, zaero_add.add_flutter_table_object),
            "FIXMDEN": (FIXMDEN, zaero_add.add_flutter_table_object),
            # trim
            "TRIM": (TRIM_ZAERO, add_methods.add_trim_object),
            "TABLED1": (TABLED1_ZAERO, add_methods.add_tabled_object),
            "TABDMP1": (TABDMP1_ZAERO, add_methods.add_table_sdamping_object),
            "TRIMFNC": (TRIMFNC, zaero_add.add_trimfnc_object),
            "PLTTIME": (PLTTIME, zaero_add.add_plttime_object),
            "TRIMVAR": (TRIMVAR, zaero_add.add_trimvar_object),
            "TRIMLNK": (TRIMLNK, zaero_add.add_trimlnk_object),
            # geometry
            "SPLINE1": (SPLINE1_ZAERO, add_methods.add_spline_object),
            "SPLINE2": (SPLINE2_ZAERO, add_methods.add_spline_object),
            "SPLINE3": (SPLINE3_ZAERO, add_methods.add_spline_object),
            "SPLINEM": (SPLINEM, zaero_add.add_splinem_object),
            "PANLST1": (PANLST1, zaero_add.add_panlst_object),
            "PANLST2": (PANLST2, zaero_add.add_panlst_object),
            "PANLST3": (PANLST3, zaero_add.add_panlst_object),
            "PAFOIL7": (PAFOIL7, zaero_add.add_pafoil_object),
            "PAFOIL8": (PAFOIL8, zaero_add.add_pafoil_object),
            "MKAEROZ": (MKAEROZ, zaero_add.add_mkaeroz_object),
            "SEGMESH": (SEGMESH, add_methods.add_paero_object),
            "BODY7": (BODY7, add_methods.add_caero_object),
            "CAERO7": (CAERO7, add_methods.add_caero_object),
            "ACOORD": (ACOORD, add_methods.add_coord_object),
            "AESURFZ": (AESURFZ, zaero_add.add_aesurfz_object),
            "AESLINK": (AESLINK, zaero_add.add_aeslink_object),
            # flutter
            "FLUTTER": (FLUTTER_ZAERO, add_methods.add_flutter_object),
            "PLTVG": (PLTVG, zaero_add.add_pltvg_object),
            "PLTFLUT": (PLTFLUT, zaero_add.add_pltflut_object),
            # mloads
            "MLOADS": (MLOADS, zaero_add.add_mloads_object),
            # gust
            "GLOADS": (GLOADS, zaero_add.add_gloads_object),
            "DGUST": (DGUST, zaero_add.add_dgust_object),
            "CGUST": (CGUST, zaero_add.add_cgust_object),
            # ase
            "ASE": (ASE, zaero_add.add_ase_object),
            "ASECONT": (ASECONT, zaero_add.add_asecont_object),
            "ASESNSR": (ASESNSR, zaero_add.add_asesnsr_object),
            "ASESNS1": (ASESNS1, zaero_add.add_asesns1_object),
            "ASEGAIN": (ASEGAIN, zaero_add.add_asegain_object),
            "GAINSET": (GAINSET, zaero_add.add_gainset_object),
            "PLTBODE": (PLTBODE, zaero_add.add_pltbode_object),
            "MIMOSS": (MIMOSS, zaero_add.add_mimoss_object),
            "SISOTF": (SISOTF, zaero_add.add_sisotf_object),
            "CJUNCT": (CJUNCT, zaero_add.add_cjunct_object),
            "CONCT": (CONCT, zaero_add.add_conct_object),
            "AEROLAG": (AEROLAG, zaero_add.add_aerolag_object),
            "MIMOTF": (MIMOTF, zaero_add.add_mimotf_object),
            "SENSR": (SENSR, zaero_add.add_sensr_object),
            "GAIN": (GAIN, zaero_add.add_gain_object),
            "SUMBLK": (SUMBLK, zaero_add.add_sumblk_object),
            "DEADBN": (DEADBN, zaero_add.add_deadbn_object),
            "DELAY": (DELAY_ZAERO, zaero_add.add_delay_zaero_object),
            "FILTFL": (FILTFL, zaero_add.add_filtfl_object),
            "LIMTR": (LIMTR, zaero_add.add_limtr_object),
            # other
            "SETADD": (SETADD, zaero_add.add_setadd_object),
            "SENSET": (SENSET, zaero_add.add_senset_object),
            "CNCTSET": (CNCTSET, zaero_add.add_cnctset_object),
            "SURFSET": (SURFSET, zaero_add.add_surfset_object),
            "ACTU": (ACTU, zaero_add.add_actu_object),
            "LOADMOD": (LOADMOD, zaero_add.add_loadmod_object),
            "RBRED": (RBRED, zaero_add.add_rbred_object),
            "ATTACH": (ATTACH, zaero_add.add_attach_object),
            "PLTMODE": (PLTMODE, zaero_add.add_pltmode_object),
            "PLTAERO": (PLTAERO, zaero_add.add_pltaero_object),
            "PLTCP": (PLTCP, zaero_add.add_pltcp_object),
            "PLTTRIM": (PLTTRIM, zaero_add.add_plttrim_object),
            "PLTSURF": (PLTSURF, zaero_add.add_pltsurf_object),
            "PLTMIST": (PLTMIST, zaero_add.add_pltmist_object),
            "EXTINP": (EXTINP, zaero_add.add_extinp_object),
            "EXTOUT": (EXTOUT, zaero_add.add_extout_object),
            "TFSET": (TFSET, zaero_add.add_tfset_object),
            "MLDSTAT": (MLDSTAT, zaero_add.add_mldstat_object),
            "MINSTAT": (MINSTAT, zaero_add.add_minstat_object),
            "CONMLST": (CONMLST, zaero_add.add_conmlst_object),
            "CPFACT": (CPFACT, zaero_add.add_cpfact_object),
            "APCONST": (APCONST, zaero_add.add_apconst_object),
            "SPLINE0": (SPLINE0, zaero_add.add_spline0_object),
            "PBODY7": (PBODY7, zaero_add.add_pbody7_object),
            "TRIMFLT": (TRIMFLT, zaero_add.add_trimflt_object),
            "ASEOUT": (ASEOUT, zaero_add.add_aseout_object),
            "CMARGIN": (CMARGIN, zaero_add.add_cmargin_object),
            "OUTPUT4": (OUTPUT4, zaero_add.add_output4_object),
            "CROSPSD": (CROSPSD, zaero_add.add_crospsd_object),
            "FOILSEC": (FOILSEC, zaero_add.add_foilsec_object),
            "GENGUST": (GENGUST, zaero_add.add_gengust_card_object),
            "MFTGUST": (MFTGUST, zaero_add.add_mftgust_object),
            "DMIS": (DMIS, zaero_add.add_dmis_object),
            "CELLWNG": (CELLWNG, zaero_add.add_cellwng_object),
            "CELLBOX": (CELLBOX, zaero_add.add_cellbox_object),
            "INPCFD": (INPCFD, zaero_add.add_inpcfd_object),
            "OMITCFD": (OMITCFD, zaero_add.add_omitcfd_object),
            "WT1AJJ": (WT1AJJ, zaero_add.add_wt1ajj_object),
            "WT1FRC": (WT1FRC, zaero_add.add_wt1frc_object),
            "WT2AJJ": (WT2AJJ, zaero_add.add_wt2ajj_object),
            "WTUCP": (WTUCP, zaero_add.add_wtucp_object),
            "TRIMOBJ": (TRIMOBJ, zaero_add.add_trimobj_object),
            "TRIMCON": (TRIMCON, zaero_add.add_trimcon_object),
            "APCNSND": (APCNSND, zaero_add.add_apcnsnd_object),
            "APCNSCP": (APCNSCP, zaero_add.add_apcnscp_object),
            "MLDTRIM": (MLDTRIM, zaero_add.add_mldtrim_object),
            "MLDCOMD": (MLDCOMD, zaero_add.add_mldcomd_object),
            "MLDTIME": (MLDTIME, zaero_add.add_mldtime_object),
            "DMIL": (DMIL, zaero_add.add_dmil_object),
            "EXTFILE": (EXTFILE, zaero_add.add_extfile_object),
            "MLDPRNT": (MLDPRNT, zaero_add.add_mldprnt_object),
        }
        skip_keys = [
            "TRIM",
            "TABLED1",
            "TABDMP1",
            "SPLINE1",
            "SPLINE2",
            "SPLINE3",
        ]
        for key in card_parser2:
            assert key in ZAERO_CARDS or key in skip_keys, f"add key={key!r} to zaero card_parser2"
        skip_cards = ["AEFACT", "CORD2R", "SET1"]
        for key in ZAERO_CARDS:
            assert key in card_parser2 or key in skip_cards, f"add key={key!r} to card_parser2"
        card_parser.update(card_parser2)
        self.model.cards_to_read.update(set(ZAERO_CARDS))
        # print('update for zona!!!!!!!!!!!')

    def cross_reference(self):
        model = self.model
        if model.nastran_format not in {"zona", "zaero"}:
            return

        # these will be xref'd twice
        for caero in model.caeros.values():
            caero.cross_reference(model)
            self.caero_to_name_map[caero.label] = caero.eid

        singletons, dicts, dicts_list = get_dicts(self, "xref")
        for items in dicts_list:
            for item in items.values():
                for itemi in item:
                    itemi.cross_reference(model)

        for items in dicts:
            for item in items.values():
                # self.model.log.info(f'xref {item.type}')
                item.cross_reference(model)

        for unused_id, panlst in self.panlsts.items():
            for panlsti in panlst:
                panlsti.cross_reference(model)
        self._checks()

    def _checks(self):
        self.view_block_diagram()
        # self._check_tfset_cjunct()
        # self._check_cntcset_conct()

    def view_block_diagram(self, subcase_id: int = -1) -> None:
        from .zaero_interface.graphviz_interface import view_block_diagram
        view_block_diagram(self.model, subcase_id=subcase_id)

    def _check_tfset_cjunct(self):  # pragma: no cover
        assert len(self.tfset) == 1, self.tfset
        sisotfs = set(list(self.sisotf))
        cjunct_ids = set(list(self.cjunct))
        expected_tfs = cjunct_ids.union(sisotfs)
        for tfset_id, tfset in self.tfset.items():
            tfset_ids = set(tfset.ids)

        extra = expected_tfs - tfset_ids
        missing = tfset_ids - expected_tfs
        assert len(extra) == 0, f"There are more CJUNCTs than values in TFSET; extra={extra}"
        assert len(missing) == 0, f"There are fewer CJUNCTs than values in TFSET; missing={missing}"

    def _check_cntcset_conct(self):  # pragma: no cover
        # assert len(self.cnctset) in [0, 1], len(self.cnctset)
        cntcset_ids = set()
        for idi, cnctset in self.cnctset.items():
            cntcset_ids.update(cnctset.ids)
        contc_ids = set(list(self.conct))
        extra = contc_ids - cntcset_ids
        missing = cntcset_ids - contc_ids
        assert len(extra) == 0, f"There are more CONTCs than values in CNTCADD; extra={extra}"
        assert len(missing) == 0, (
            f"There are fewer CONTCs than values in CNTCADD; missing={missing}"
        )
        all_blocks = []
        # assert len(self.mimoss) == 0, self.mimoss
        for idi, card in self.sisotf.items():
            all_blocks.append(f"{card.type}={idi}-1")
        for idi, card in self.mimoss.items():
            all_blocks.append(f"{card.type}={idi}-1")
        for idi, card in self.actu.items():
            all_blocks.append(f"{card.type}={idi}-1")
        for idi, card in self.cjunct.items():
            # i = 0
            for i in range(1, card.nu + 1):  # inputs
                all_blocks.append(f"{card.type}={idi}-{i}")
            # for j in range(i, card.ny+1):
            #     all_blocks.append(f'{card.type}={idi}-{j}')
            # nu: 2
            # ny: 1

        log = self.model.log
        print(f"all_blocks = {all_blocks}")
        for idi, card in self.conct.items():
            # 'SISOTF=31004-1', 'ACTU=21001-1', 'ACTU=21002-1
            # print(card)
            input_ref = card.input_ref
            if input_ref is None:
                log.warning(f"missing input-type for:\n{str(card)}")
                continue

            input_name = f"{input_ref.type}={card.input_tf_id}-{card.input_component}"
            log.debug(f"found {input_name}")
            assert input_name in all_blocks, f"input={input_name!r} not in all_blocks\n{str(card)}"

            output_ref = card.output_ref
            if output_ref is None:
                log.warning(f"missing output-type for:\n{str(card)}")
                continue
            output_name = f"{output_ref.type}={card.output_tf_id}-{card.output_component}"
            if output_name not in all_blocks:
                log.warning(f"output={output_name!r} not in all_blocks\n{str(card)}")
                continue
            log.debug(f"found {output_name}")
            # asdf
        # for
        # asdf

    def safe_cross_reference(self, xref_errors=None):
        model = self.model
        if model.nastran_format not in {"zona", "zaero"}:
            return
        if xref_errors is None:
            xref_errors = defaultdict(list)

        for caero in model.caeros.values():
            caero.safe_cross_reference(model, xref_errors)
            self.caero_to_name_map[caero.label] = caero.eid

        singletons, dicts, dicts_list = get_dicts(self, "xref")
        for items in dicts_list:
            for item in items.values():
                for itemi in item:
                    itemi.safe_cross_reference(model, xref_errors)

        for items in dicts:
            for item in items.values():
                # self.model.log.info(f'xref {item.type}')
                item.safe_cross_reference(model, xref_errors)

        for unused_id, panlst in self.panlsts.items():
            for panlsti in panlst:
                panlsti.safe_cross_reference(model, xref_errors)
        self._checks()

    def uncross_reference(zaero: ZAERO):
        singletons, dicts, dicts_list = get_dicts(zaero, "write")
        for panlsts in zaero.panlsts.values():
            for panlst in panlsts:
                panlst.uncross_reference()

        for items in dicts_list:
            for item in items.values():
                for itemi in item:
                    itemi.uncross_reference()

        for dicti in dicts:
            if isinstance(dicti, list):
                if len(dicti):
                    print(dicti)
                    asdf
                continue
            for value in dicti.values():
                value.uncross_reference()

    def write_bdf(
        self, bdf_file: TextIO, size: int = 8, is_double: bool = False, sort_cards: bool = True
    ):
        # if self.model.nastran_format != 'zona':
        # return
        for unused_id, panlst in self.panlsts.items():
            for panlsti in panlst:
                bdf_file.write(panlsti.write_card(size=size, is_double=is_double))
        # for unused_id, mkaeroz in self.mkaeroz.items():
        #     bdf_file.write(mkaeroz.write_card(size=size, is_double=is_double))
        # for unused_id, trimvar in self.trimvar.items():
        #     bdf_file.write(trimvar.write_card(size=size, is_double=is_double))
        # for unused_id, trimlnk in self.trimlnk.items():
        #     bdf_file.write(trimlnk.write_card(size=size, is_double=is_double))
        # for unused_id, pafoil in self.pafoil.items():
        #     bdf_file.write(pafoil.write_card(size=size, is_double=is_double))
        # for unused_id, attach in self.attach.items():
        #     bdf_file.write(attach.write_card(size=size, is_double=is_double))

        singletons, dicts, dicts_list = get_dicts(self, "write")
        for item in singletons:
            if item is not None:
                bdf_file.write(item.write_card(size=size, is_double=is_double))

        for items in dicts_list:
            for item in items.values():
                for itemi in item:
                    bdf_file.write(itemi.write_card(size=size, is_double=is_double))

        for items in dicts:
            for key, value in sorteddict(items, sort_cards):
                bdf_file.write(value.write_card(size=size, is_double=is_double))

    def convert_to_nastran(self, save: bool = True):
        """Converts a ZAERO model to Nastran"""
        from pyNastran.bdf.cards.aero.zaero_interface.zaero_to_nastran import zaero_to_nastran

        return zaero_to_nastran(self, save=save)

    def add_caero2s(self, caero2s, add=False):
        """Converts ZAERO BODY7 to CAERO2/PAERO2/AEFACT"""
        model = self.model
        add_methods = model._add_methods
        caero_body_ids = []
        for caero2 in caero2s:
            caero_id = caero2.eid
            out = caero2.convert_to_nastran(model)
            caero_new, paero2, aefact_xs, aefact_width, aefact_theta1, aefact_theta2 = out
            caero_body_ids.append(caero_id)
            if add:
                add_methods.add_aefact_object(aefact_xs)
                add_methods.add_aefact_object(aefact_width)
                add_methods.add_aefact_object(aefact_theta1)
                add_methods.add_aefact_object(aefact_theta2)
                add_methods.add_paero_object(paero2)
                add_methods.add_caero_object(caero_new)
        return

    def __repr__(self):
        msg = "<ZAERO>; nPANLSTs=%s nmkaeroz=%s" % (
            len(self.panlsts),
            len(self.mkaeroz),
        )
        return msg

    def get_bdf_cards(
        self,
        bulk_data_lines: list[str],
        bulk_data_ilines: Optional[np.ndarray] = None,
        use_dict: bool = False,
    ) -> tuple[Any, Any, Any]:
        """Parses the BDF lines into a list of card_lines"""
        # self.log.warning('get_bdf_cards')
        use_dict = False

        dict_cards = {}
        model = self.model
        allow_tabs = False
        if bulk_data_ilines is None:
            bulk_data_ilines = np.zeros((len(bulk_data_lines), 2), dtype="int32")

        cards_list: list[Any] = []
        cards_dict: dict[str, list[Any]] = defaultdict(list)
        # cards = defaultdict(list)
        card_count: dict[str, int] = defaultdict(int)
        full_comment = ""
        card_lines = []
        card_ilines = []
        old_ifile_iline = None
        old_card_name = None
        backup_comment = ""
        nlines = len(bulk_data_lines)

        # self.echo = True
        # self.force_echo_off = False

        log = self.model.log
        for iline_bulk, line in enumerate(bulk_data_lines):
            ifile_iline = bulk_data_ilines[iline_bulk, :]
            # print(iline_bulk, ifile_iline)
            # print(iline_bulk, ifile_iline, line)
            # print('    backup={backup_comment!r}')
            comment = ""
            if "$" in line and (line.lstrip().startswith("$") or line.index("$") >= 72):
                line, comment = line.split("$", 1)
                comment = comment.rstrip()
                # if line.strip():
                #     print(line)
                strip_comment = comment.strip()
                if strip_comment.lower().startswith("group:"):
                    continue_flag = model._store_group(strip_comment)
                    if continue_flag:
                        continue
            # if not self.allow_tabs and '\t' in line:
            # raise RuntimeError(f'There are tabs in:\n{line}')
            # self.log.warning(f'There are tabs in:\n{line}')

            card_name = line.split(",", 1)[0].split("\t", 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ["+", "*"]:
                if old_card_name:
                    # multiline card is finished
                    if card_name not in ZAERO_CARDS:
                        pass
                    elif not allow_tabs and "\t" in (joined_lines_n := "\n".join(card_lines)):
                        joined_lines_n2 = "\n".join((f"{line!r}" for line in card_lines))
                        log.warning(f"There are tabs in:\n{joined_lines_n2}")
                        # raise RuntimeError(f'There are tabs in:\n{joined_lines_n2}')

                    if model.echo and not model.force_echo_off:
                        model.log.info(
                            "Reading %s:\n" % old_card_name + full_comment + "".join(card_lines)
                        )

                    # if full_comment:
                    # print('full_comment = ', full_comment)
                    cards_list.append(
                        [old_card_name, _prep_comment(full_comment), card_lines, card_ilines[-1]]
                    )

                    card_count[old_card_name] += 1
                    card_lines = []
                    card_ilines = []
                    full_comment = ""

                    if old_card_name == "ECHOON":
                        self.echo = True
                    elif old_card_name == "ECHOOFF":
                        self.echo = False
                old_ifile_iline = ifile_iline
                old_card_name = card_name.rstrip(" *")

                if old_card_name == "ENDDATA":
                    model.card_count["ENDDATA"] = 1
                    if nlines - iline_bulk > 1:
                        nleftover = nlines - iline_bulk - 1
                        msg = "exiting due to ENDDATA found with %i lines left" % nleftover
                        model.log.debug(msg)
                    cards_list, cards_dict, card_count = fix_card_list(
                        cards_list, cards_dict, card_count
                    )
                    return cards_list, cards_dict, card_count
                # print("card_name = %s" % card_name)

            comment = _clean_comment(comment)

            # TODO: these additional \n need to be there for rejected cards
            #      but not parsed cards
            if line.rstrip():
                card_lines.append(line)
                card_ilines.append(ifile_iline)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + comment + "\n"
                    else:
                        full_comment += backup_comment
                    backup_comment = ""
                elif comment:
                    full_comment += comment + "\n"
                    backup_comment = ""

            elif comment:
                backup_comment += comment + "\n"
            # elif comment:
            # backup_comment += '$' + comment + '\n'

        if card_lines:
            if not allow_tabs and "\t" in (joined_lines_n := "\n".join(card_lines)):
                log.error(f"There are tabs in:\n{joined_lines_n}")
                # raise RuntimeError(f'There are tabs in:\n{joined_lines_n}')

            if model.echo and not model.force_echo_off:
                model.log.info("Reading %s:\n" % old_card_name + full_comment + "".join(card_lines))
            # print('end_add %s' % card_lines)

            # old dictionary version
            # cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            # if backup_comment + full_comment:
            # print('backup_comment + full_comment = ', backup_comment + full_comment)
            if old_card_name in dict_cards:
                cards_dict[old_card_name].append(
                    [_prep_comment(backup_comment + full_comment), card_lines, ifile_iline]
                )
            else:
                # cards_list.append([old_card_name, _prep_comment(
                #     backup_comment + full_comment), card_lines, ifile_iline])
                cards_list.append(
                    [
                        old_card_name,
                        _prep_comment(backup_comment + full_comment),
                        card_lines,
                        card_ilines[-1],
                    ]
                )
            card_count[old_card_name] += 1
        self.echo = False

        cards_list, cards_dict, card_count = fix_card_list(cards_list, cards_dict, card_count)
        return cards_list, cards_dict, card_count


def fix_card_list(cards_list, cards_dict, card_count):
    assert len(cards_dict) == 0, cards_dict
    card_count = defaultdict(int)
    # skip_cards = {'SPLINE1', 'SPLINE2', 'SPLINE3', 'AEFACT', 'CONM2',
    #               'GENEL', 'DMI', 'DMIG', 'TABLED1', 'TABLED2',
    #               'EIGR', 'EIGRL', 'DAREA'}
    include_cards = {
        "GRID",
        "SPOINT",
        "EPOINT",
        "CORD2R",
        "CORD2C",
        "CORD2S",
        "CORD1R",
        "CORD1C",
        "CORD1S",
        "CROD",
        "PROD",
        "CONROD",
        "CBUSH",
        "PBUSH",
        "CBUSH1D",
        "PBUSH1D",
        "CBAR",
        "PBAR",
        "PBARL",
        "CBEAM",
        "PBEAM",
        "PBEAML",
        "CTETRA",
        "CHEXA",
        "CPENTA",
        "CPYRAM",
        "PSOLID",
        "PLSOLID",
        "CQUAD4",
        "CTRIA3",
        "CQUAD8",
        "CTRIA6",
        "CTRIAR",
        "CQUADR",
        "CQUAD",
        "PSHELL",
        "PCOMP",
        "PCOMPG",
        "MAT1",
        "MAT2",
        "MAT8",
        "MAT9",
        "MAT10",
    }
    cards_list2 = []
    for card_name, comment, card_lines, ifile_iline in cards_list:
        # (ifile, iline) = ifile_iline
        if ifile_iline[0] == 1000 and card_name not in include_cards:  # f06/prt
            continue
        card_count[card_name] += 1
        cards_list2.append((card_name, comment, card_lines, ifile_iline))
    return cards_list2, cards_dict, dict(card_count)


def get_dicts(zaero: ZAERO, method: str) -> tuple[list, dict[int, list], list[dict]]:
    assert method in ["xref", "write"], f"method={method!r}"
    dicts = [
        # --------------general-------------
        # zaero.aeroz,
        zaero.atmos, zaero.flutter_table,
        # -------------geometry-------------
        # zaero.panlsts,  # special-list
        zaero.pafoil, zaero.attach,
        zaero.pltsurf, zaero.pltmode, zaero.pltmist, zaero.pltbode,
        # -------------transient------------
        zaero.mloads, zaero.eloads,
        # --------------flutter-------------
        zaero.nlfltr, zaero.mkaeroz,
        # ---------------trim---------------
        # zaero.trim,  # part of the main BDF
        zaero.aeslink, zaero.trimvar, zaero.trimlnk, zaero.trimfnc,
        zaero.trimobj, zaero.trimcon, zaero.actu,
        # ---------------ase---------------
        zaero.cjunct, zaero.conct, zaero.tfset, zaero.cnctset, zaero.ase,
        zaero.asecont, zaero.asesnsr, zaero.asesns1,
        zaero.asegain, zaero.gainset,
        zaero.mimoss, zaero.mimotf, zaero.sisotf,
        zaero.sensr, zaero.gain, zaero.sumblk, zaero.deadbn,
        zaero.delay_zaero, zaero.filtfl, zaero.limtr,
        #
        zaero.senset, zaero.surfset, zaero.mldtrim, zaero.mldstat,
        zaero.minstat, zaero.mldprnt, zaero.mldcomd, zaero.mldtime,
        # zaero.extinp, zaero.extout,
        zaero.loadmod,
        zaero.rbred,
        zaero.aerolag,
        # ---------------gust---------------
        zaero.gloads, zaero.dgust, zaero.cgust,
        # ---------------other--------------
        zaero.extfile,
        zaero.dse, zaero.dmil,
        # plotting
        zaero.pltvg, zaero.pltbode, zaero.pltaero,
    ]
    dict_lists = [
        zaero.pltcp, zaero.pltflut, zaero.plttime, zaero.plttrim,
    ]
    if method == "write":
        # these are xref'd by their parent
        dicts.extend([zaero.setadd, zaero.extinp, zaero.extout])
    singletons = [zaero.splinem]
    return singletons, dicts, dict_lists
