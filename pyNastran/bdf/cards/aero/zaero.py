# coding: utf-8
# pylint: disable=W0212,C0103
from __future__ import annotations

from collections import defaultdict
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils import (
    object_attributes, object_methods, PathLike)

from typing import Any
from pyNastran.bdf.bdf_interface.utils import sorteddict
from pyNastran.bdf.bdf_interface.utils import _prep_comment
from pyNastran.bdf.bdf_interface.pybdf import _clean_comment

from .zaero_cards.zaero_sets import (
    SETADD)
from pyNastran.bdf.cards.aero.zaero_cards.atm import (
    ATMOS, FIXMATM, FIXHATM, FIXMACH, FIXMDEN)
from pyNastran.bdf.cards.aero.zaero_cards.spline import (
    SPLINE1_ZAERO, SPLINE2_ZAERO, SPLINE3_ZAERO, SPLINEM)
from pyNastran.bdf.cards.aero.zaero_cards.geometry import (
    PANLST1, PANLST2, PANLST3, SEGMESH,
    CAERO7, BODY7, PAFOIL7, PAFOIL8, AESURFZ, AESLINK)
from pyNastran.bdf.cards.aero.zaero_cards.plot import (
    PLTAERO, PLTMODE, PLTVG, PLTFLUT, PLTTIME,
    PLTCP, PLTMIST, PLTSURF, PLTBODE, PLTTRIM)
from pyNastran.bdf.cards.aero.zaero_cards.flutter import (
    FLUTTER_ZAERO, MKAEROZ)
from pyNastran.bdf.cards.aero.zaero_cards.trim import (
    TRIM_ZAERO, TRIMVAR, TRIMLNK,)
from pyNastran.bdf.cards.aero.zaero_cards.manuever import (
    MLOADS, EXTINP, EXTOUT, TRIMFNC, ACTU, LOADMOD, RBRED,)
from pyNastran.bdf.cards.aero.zaero_cards.gust import (
    GLOADS, DGUST, CGUST, MFTGUST)
from pyNastran.bdf.cards.aero.zaero_cards.ase import (
    ASE, ASECONT, ASESNSR, ASESNS1,
    CJUNCT, CONCT, TFSET, MIMOSS, SISOTF,
    SENSET, SURFSET, CNCTSET,
    ASEGAIN, GAINSET,
    AEROLAG,
)
from pyNastran.bdf.cards.aero.zaero_cards.bdf_tables import (
    TABLED1_ZAERO, TABDMP1_ZAERO)
from pyNastran.bdf.cards.aero.zaero_cards.dmi import DMIL
from pyNastran.bdf.cards.aero.zaero_cards.cards import (
    MLDPRNT, MLDSTAT, MINSTAT, MLDTRIM, MLDCOMD, MLDTIME,
    AEROZ, ACOORD, ATTACH, EXTFILE,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

ZAERO_CARDS = [
    # atmosphere
    'ATMOS',
    'FIXHATM', 'FIXMATM', 'FIXMACH', 'FIXMDEN',

    # already added
    'AEFACT', 'CORD2R',

    # geometry
    'CAERO7', 'AEROZ', 'AESURFZ', 'AESLINK',
    'ATTACH',
    'PANLST1', 'PANLST2', 'PANLST3',
    'PAFOIL7', 'PAFOIL8',
    'SEGMESH', 'BODY7', 'ACOORD', 'MKAEROZ',
    'SPLINEM',
    # -------------
    # trim
    'PLTCP',
    'TRIMVAR', 'TRIMLNK',
    'TRIMFNC', # optimization
    'PLTTRIM',
    # -------------
    # flutter
    'FLUTTER',
    'FIXHATM', 'FIXMATM', # 'FIXMDEN',
    'PLTVG', 'PLTFLUT', 'PLTSURF',
    # -------------
    # plotting
    'PLTMODE', 'PLTAERO', #'PLTTIME',
    'PLTMIST',
    # -------------
    # mloads
    'MLOADS',
    'CJUNCT', 'CONCT', 'TFSET',
    'MLDSTAT', 'MLDTRIM', 'MLDPRNT', 'MLDCOMD',
    'EXTINP', 'EXTOUT', 'LOADMOD',
    'PLTTIME',
    # -------------
    # gust
    'GLOADS', 'DGUST', 'CGUST',
    # -------------
    # ase
    'ASE', 'ASECONT',
    'ASESNSR', 'ASESNS1', 'SENSET',
    'ACTU',
    'MIMOSS', 'SISOTF',
    'ASEGAIN', 'GAINSET',
    'PLTBODE', 'AEROLAG',
    # -------------
    # other
    'SETADD',
    # 'DMIL',
    'EXTFILE',
    'MLDTIME', 'MLDCOMD',
    'MINSTAT', #'APCONST',
    'RBRED',
    # 'SPLINE0', 'PBODY7',
    'CNCTSET', 'SURFSET',
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
        assert key not in self.mldprnt
        assert key > 0, key
        self.mldprnt[key] = mldprnt
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

    def add_mldprnt_object(self, mldprnt: MLDPRNT) -> None:
        """adds an MLDPRNT object"""
        key = mldprnt.mldprnt_id
        if key in self.model.zaero.mldprnt:
            self.model.log.warning(f'skipping MLDPRNT\n{str(mldprnt)}')
            return
        assert key not in self.model.zaero.mldprnt, key
        assert key > 0, key
        self.model.zaero.mldprnt[key] = mldprnt
        self.model._type_to_id_map[mldprnt.type].append(key)

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
        assert key not in zaero.pafoil, '\npafoil7=\n%s old=\n%s' % (
            pafoil, zaero.pafoil[key])
        zaero.pafoil[key] = pafoil
        self.model._type_to_id_map[pafoil.type].append(key)

    def add_aesurfz_object(self, aesurf: AESURFZ) -> None:
        """adds an AESURFZ object"""
        key = aesurf.aesid
        model = self.model
        assert key not in model.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, model.aesurf[key])
        model.aesurf[key] = aesurf
        model._type_to_id_map[aesurf.type].append(key)

    def add_aeslink_object(self, aeslink: AESLINK) -> None:
        """adds an AESLINK object"""
        key = aeslink.label
        model = self.model
        zaero = model.zaero
        assert key not in zaero.aeslink, '\naeslink=\n%s old=\n%s' % (
            aeslink, zaero.aeslink[key])
        zaero.aeslink[key] = aeslink
        model._type_to_id_map[aeslink.type].append(key)

    def add_mloads_object(self, mloads: MLOADS) -> None:
        """adds an MLOADS object"""
        key = mloads.mloads_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.mloads, '\nmloads=\n%s old=\n%s' % (
            mloads, zaero.mloads[key])
        zaero.mloads[key] = mloads
        model._type_to_id_map[mloads.type].append(key)

    def add_gloads_object(self, gloads: GLOADS) -> None:
        """adds an GLOADS object"""
        key = gloads.gloads_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.gloads, '\ngloads=\n%s old=\n%s' % (
            gloads, model.gloads[key])
        zaero.gloads[key] = gloads
        model._type_to_id_map[gloads.type].append(key)

    def add_dgust_object(self, dgust: DGUST) -> None:
        """adds an DGUST object"""
        key = dgust.dgust_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.dgust, '\ndgust=\n%s old=\n%s' % (
            dgust, zaero.dgust[key])
        zaero.dgust[key] = dgust
        model._type_to_id_map[dgust.type].append(key)

    def add_cgust_object(self, cgust: CGUST) -> None:
        """adds an CGUST object"""
        key = cgust.cgust_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.cgust, '\ncgust=\n%s old=\n%s' % (
            cgust, zaero.cgust[key])
        zaero.cgust[key] = cgust
        model._type_to_id_map[cgust.type].append(key)

    def add_ase_object(self, ase: ASE) -> None:
        """adds an ASE object"""
        key = ase.ase_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.ase, '\nase=\n%s old=\n%s' % (
            ase, zaero.ase[key])
        zaero.ase[key] = ase
        model._type_to_id_map[ase.type].append(key)

    def add_asecont_object(self, asecont: ASECONT) -> None:
        """adds an ASECONT object"""
        key = asecont.asecont_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asecont, '\nasecont=\n%s old=\n%s' % (
            asecont, zaero.asecont[key])
        zaero.asecont[key] = asecont
        model._type_to_id_map[asecont.type].append(key)

    def add_asegain_object(self, asegain: ASEGAIN) -> None:
        """adds an ASEGAIN object"""
        key = asegain.asegain_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asegain, '\nasegain=\n%s old=\n%s' % (
            asegain, zaero.asegain[key])
        zaero.asegain[key] = asegain
        model._type_to_id_map[asegain.type].append(key)

    def add_asesnsr_object(self, asesnsr: ASESNSR) -> None:
        """adds an ASESNSR object"""
        key = asesnsr.asesnsr_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asesnsr, '\nasesnsr=\n%s old=\n%s' % (
            asesnsr, zaero.asesnsr[key])
        zaero.asesnsr[key] = asesnsr
        model._type_to_id_map[asesnsr.type].append(key)

    def add_asesns1_object(self, asesns1: ASESNS1) -> None:
        """adds an ASESNS1 object"""
        key = asesns1.asesns1_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.asesns1, '\nasesns1=\n%s old=\n%s' % (
            asesns1, zaero.asesns1[key])
        zaero.asesns1[key] = asesns1
        model._type_to_id_map[asesns1.type].append(key)

    def add_gainset_object(self, gainset: GAINSET) -> None:
        """adds an GAINSET object"""
        key = gainset.gainset_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.gainset, '\ngainset=\n%s old=\n%s' % (
            gainset, zaero.gainset[key])
        zaero.gainset[key] = gainset
        model._type_to_id_map[gainset.type].append(key)

    def add_pltbode_object(self, pltbode: PLTBODE) -> None:
        """adds an PLTBODE object"""
        key = pltbode.set_id
        model = self.model
        zaero = model.zaero
        if key in zaero.pltbode:
            model.log.warning(f'skipping duplicate PLTBODE\n{str(pltbode)}')
            return
        assert key not in zaero.pltbode, '\npltbode=\n%s old=\n%s' % (
            pltbode, zaero.pltbode[key])
        zaero.pltbode[key] = pltbode
        model._type_to_id_map[pltbode.type].append(key)

    def add_cjunct_object(self, cjunct: CJUNCT) -> None:
        """adds an CJUNCT object"""
        key = cjunct.cjunct_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.cjunct, '\ncjunct=\n%s old=\n%s' % (
            cjunct, zaero.cjunct[key])
        zaero.cjunct[key] = cjunct
        model._type_to_id_map[cjunct.type].append(key)

    def add_conct_object(self, conct: CONCT) -> None:
        """adds an CONCT object"""
        key = conct.conct_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.conct, '\nconct=\n%s old=\n%s' % (
            conct, zaero.conct[key])
        zaero.conct[key] = conct
        model._type_to_id_map[conct.type].append(key)

    def add_aerolag_object(self, aerolag: AEROLAG) -> None:
        """adds an CONCT object"""
        key = aerolag.aerolag_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.conct, '\naerolag=\n%s old=\n%s' % (
            aerolag, zaero.aerolag[key])
        zaero.aerolag[key] = aerolag
        model._type_to_id_map[aerolag.type].append(key)

    def add_actu_object(self, actu: ACTU) -> None:
        """adds an AESURFZ object"""
        key = actu.actu_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.actu, '\nactu=\n%s old=\n%s' % (
            actu, zaero.actu[key])
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
        assert trimvar.var_id not in self.model.zaero.trimvar, '\ntrimvar=\n%s old=\n%s' % (
            trimvar, self.model.zaero.trimvar[key])
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

    def add_flutter_table_object(self, flutter_table: FIXHATM | FIXMATM | FIXMDEN | FIXMACH) -> None:
        """adds an FIXMATM object"""
        key = flutter_table.sid
        model = self.model
        zaero = model.zaero
        assert key not in zaero.flutter_table, '\nflutter_table=\n%s old=\n%s' % (
            flutter_table, zaero.flutter_table[key])
        zaero.flutter_table[key] = flutter_table
        model._type_to_id_map[flutter_table.type].append(key)

    def add_atmos_object(self, atmos: ATMOS) -> None:
        """adds an ATMOS object"""
        key = atmos.atmos_id
        model = self.model
        zaero = model.zaero
        assert key not in zaero.atmos, '\natmos=\n%s old=\n%s' % (
            atmos, zaero.atm[key])
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
        self.dse: dict[int, DSE] = {}
        self.extfile: dict[int, EXTFILE] = {}
        # FOILSEC
        # CPFACT
        # TRIMFLT

        # transient
        self.eloads: dict[int, ELOADS] = {}
        self.mldcomd: dict[int, MLDCOMD] = {}
        self.plttime: dict[int, PLTTIME] = {}

        # trim
        self.pltcp: dict[int, PLTCP] = {}
        self.plttrim: dict[int, PLTTRIM] = {}
        self.aeslink: dict[int, AESLINK] = {}
        self.trimvar: dict[int, TRIMVAR] = {}
        self.trimlnk: dict[int, TRIMLNK] = {}
        self.mldtrim: dict[int, MLDTRIM] = {}

        # flutter
        self.pltvg: dict[int, PLTVG] = {}
        self.mkaeroz: dict[int, MKAEROZ] = {}
        self.nlfltr: dict[int, NLFLTR] = {}
        self.pltflut: dict[int, PLTFLUT] = {}

        # trim
        self.trimfnc: dict[int, TRIMFNC] = {}
        self.trimobj: dict[int, TRIMOBJ] = {}
        self.trimcon: dict[int, TRIMCON] = {}

        # gust
        self.gloads: dict[int, GLOADS] = {}
        self.dgust: dict[int, DGUST] = {}
        self.cgust: dict[int, CGUST] = {}
        self.gengust: dict[int, GENGUST] = {}
        self.gustinp: dict[int, GUSTINP] = {}

        # ase
        self.ase: dict[int, ASE] = {}
        self.asecont: dict[int, ASECONT] = {}
        self.asesnsr: dict[int, ASESNSR] = {}
        self.asesns1: dict[int, ASESNS1] = {}
        self.asegain: dict[int, ASEGAIN] = {}
        self.gainset: dict[int, GAINSET] = {}
        self.pltbode: dict[int, PLTBODE] = {}
        self.aseout: dict[int, ASEOUT] = {}
        self.apcnsnd: dict[int, APCNSND] = {}
        self.apcnscp: dict[int, APCNSCP] = {}
        self.mimoss: dict[int, MIMOSS] = {}
        self.sisotf: dict[int, SISOTF] = {}
        self.cmargin: dict[int, CMARGIN] = {}
        self.aerolag: dict[int, AEROLAG] = {}

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
        #self.aeroz = {}

    def object_attributes(self, mode: str='public', keys_to_skip: Optional[list[str]]=None,
                          filter_properties: bool=False):
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
            'log', 'model',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public',
                       keys_to_skip: Optional[list[str]]=None) -> list[str]:
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

        my_keys_to_skip = ['log',]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def verify(self, xref):
        if self.model.nastran_format not in {'zona', 'zaero'}:
            return
        for panlst in self.panlsts.values():
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
        if self.model.nastran_format not in {'zona', 'zaero'}:
            return
        for panlst in self.panlsts.values():
            for panlsti in panlst:
                panlsti.validate()

        singletons, dicts, dicts_list = get_dicts(self, 'write')
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

    def PAFOIL(self, pid: int, msg: str=''):
        """gets a pafoil profile (PAFOIL7/PAFOIL8)"""
        try:
            return self.pafoil[pid]
        except KeyError:
            pafoils = np.unique(list(self.pafoil.keys()))
            raise KeyError(f'pid={pid} not found{msg}.  Allowed pafoils={pafoils}')

    def update_for_zaero(self):
        """updates for zaero"""
        card_parser = self.model._card_parser
        add_methods = self.model._add_methods
        zaero_add = self._add_methods
        card_parser2 = {
            # aero models
            'AEROZ': (AEROZ, add_methods.add_aeros_object),
            'ATMOS': (ATMOS, zaero_add.add_atmos_object),
            'FIXMATM': (FIXMATM, zaero_add.add_flutter_table_object),
            'FIXHATM': (FIXHATM, zaero_add.add_flutter_table_object),
            'FIXMACH': (FIXMACH, zaero_add.add_flutter_table_object),
            'FIXMDEN': (FIXMDEN, zaero_add.add_flutter_table_object),
            # trim
            'TRIM': (TRIM_ZAERO, add_methods.add_trim_object),
            'TABLED1': (TABLED1_ZAERO, add_methods.add_tabled_object),
            'TABDMP1': (TABDMP1_ZAERO, add_methods.add_table_sdamping_object),
            'TRIMFNC': (TRIMFNC, zaero_add.add_trimfnc_object),
            'PLTTIME': (PLTTIME, zaero_add.add_plttime_object),
            'TRIMVAR': (TRIMVAR, zaero_add.add_trimvar_object),
            'TRIMLNK': (TRIMLNK, zaero_add.add_trimlnk_object),
            # geometry
            'SPLINE1': (SPLINE1_ZAERO, add_methods.add_spline_object),
            'SPLINE2': (SPLINE2_ZAERO, add_methods.add_spline_object),
            'SPLINE3': (SPLINE3_ZAERO, add_methods.add_spline_object),
            'SPLINEM': (SPLINEM, zaero_add.add_splinem_object),
            'PANLST1': (PANLST1, zaero_add.add_panlst_object),
            'PANLST2': (PANLST2, zaero_add.add_panlst_object),
            'PANLST3': (PANLST3, zaero_add.add_panlst_object),
            'PAFOIL7': (PAFOIL7, zaero_add.add_pafoil_object),
            'PAFOIL8': (PAFOIL8, zaero_add.add_pafoil_object),
            'MKAEROZ': (MKAEROZ, zaero_add.add_mkaeroz_object),
            'SEGMESH': (SEGMESH, add_methods.add_paero_object),
            'BODY7': (BODY7, add_methods.add_caero_object),
            'CAERO7': (CAERO7, add_methods.add_caero_object),
            'ACOORD': (ACOORD, add_methods.add_coord_object),
            'AESURFZ': (AESURFZ, zaero_add.add_aesurfz_object),
            'AESLINK': (AESLINK, zaero_add.add_aeslink_object),
            # flutter
            'FLUTTER': (FLUTTER_ZAERO, add_methods.add_flutter_object),
            'PLTVG': (PLTVG, zaero_add.add_pltvg_object),
            'PLTFLUT': (PLTFLUT, zaero_add.add_pltflut_object),
            # mloads
            'MLOADS': (MLOADS, zaero_add.add_mloads_object),
            # gust
            'GLOADS': (GLOADS, zaero_add.add_gloads_object),
            'DGUST': (DGUST, zaero_add.add_dgust_object),
            'CGUST': (CGUST, zaero_add.add_cgust_object),
            # ase
            'ASE': (ASE, zaero_add.add_ase_object),
            'ASECONT': (ASECONT, zaero_add.add_asecont_object),
            'ASESNSR': (ASESNSR, zaero_add.add_asesnsr_object),
            'ASESNS1': (ASESNS1, zaero_add.add_asesns1_object),
            'ASEGAIN': (ASEGAIN, zaero_add.add_asegain_object),
            'GAINSET': (GAINSET, zaero_add.add_gainset_object),
            'PLTBODE': (PLTBODE, zaero_add.add_pltbode_object),
            'MIMOSS': (MIMOSS, zaero_add.add_mimoss_object),
            'SISOTF': (SISOTF, zaero_add.add_sisotf_object),
            'CJUNCT': (CJUNCT, zaero_add.add_cjunct_object),
            'CONCT': (CONCT, zaero_add.add_conct_object),
            'AEROLAG': (AEROLAG, zaero_add.add_aerolag_object),
            # other
            'SETADD': (SETADD, zaero_add.add_setadd_object),
            'SENSET': (SENSET, zaero_add.add_senset_object),
            'CNCTSET': (CNCTSET, zaero_add.add_cnctset_object),
            'SURFSET': (SURFSET, zaero_add.add_surfset_object),
            'ACTU': (ACTU, zaero_add.add_actu_object),
            'LOADMOD': (LOADMOD, zaero_add.add_loadmod_object),
            'RBRED': (RBRED, zaero_add.add_rbred_object),
            'ATTACH': (ATTACH, zaero_add.add_attach_object),
            'PLTMODE': (PLTMODE, zaero_add.add_pltmode_object),
            'PLTAERO': (PLTAERO, zaero_add.add_pltaero_object),
            'PLTCP': (PLTCP, zaero_add.add_pltcp_object),
            'PLTTRIM': (PLTTRIM, zaero_add.add_plttrim_object),
            'PLTSURF': (PLTSURF, zaero_add.add_pltsurf_object),
            'PLTMIST': (PLTMIST, zaero_add.add_pltmist_object),
            'EXTINP': (EXTINP, zaero_add.add_extinp_object),
            'EXTOUT': (EXTOUT, zaero_add.add_extout_object),
            'TFSET': (TFSET, zaero_add.add_tfset_object),
            'MLDSTAT': (MLDSTAT, zaero_add.add_mldstat_object),
            'MINSTAT': (MINSTAT, zaero_add.add_minstat_object),
            'MLDTRIM': (MLDTRIM, zaero_add.add_mldtrim_object),
            'MLDCOMD': (MLDCOMD, zaero_add.add_mldcomd_object),
            'MLDTIME': (MLDTIME, zaero_add.add_mldtime_object),
            #'DMIL': (DMIL, zaero_add.add_dmil_object),
            'EXTFILE': (EXTFILE, zaero_add.add_extfile_object),
            'MLDPRNT': (MLDPRNT, zaero_add.add_mldprnt_object),
        }
        skip_keys = [
            'TRIM', 'TABLED1', 'TABDMP1',
            'SPLINE1', 'SPLINE2', 'SPLINE3',
        ]
        for key in card_parser2:
            assert key in ZAERO_CARDS or key in skip_keys, f'add key={key!r} to zaero card_parser2'
        skip_cards = ['AEFACT', 'CORD2R', 'SET1']
        for key in ZAERO_CARDS:
            assert key in card_parser2 or key in skip_cards, f'add key={key!r} to card_parser2'
        card_parser.update(card_parser2)
        self.model.cards_to_read.update(set(ZAERO_CARDS))
        # print('update for zona!!!!!!!!!!!')

    def cross_reference(self):
        model = self.model
        if model.nastran_format not in {'zona', 'zaero'}:
            return

        # these will be xref'd twice
        for caero in model.caeros.values():
            caero.cross_reference(model)
            self.caero_to_name_map[caero.label] = caero.eid

        singletons, dicts, dicts_list = get_dicts(self, 'xref')
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
        self.build_block()
        # self._check_tfset_cjunct()
        # self._check_cntcset_conct()

    def build_block(self):
        # g = graphviz.Graph('G', filename='process2.gv', engine='sfdp')
        # g = graphviz.Diagram('G', filename='process2.gv', engine='sfdp')
        try:
            import graphviz
            from graphviz import Digraph, ExecutableNotFound
        except ImportError:
            return

        if not isinstance(self.model.bdf_filename, PathLike):
            return
        # g = graphviz.Digraph('G', filename='hello2.gv')
        # g.edge('Hello', 'World')
        # g.view()
        # asdf

        filename = str(self.model.bdf_filename) + '_ase'
        g = Digraph('G', filename=filename)
        # g.attr('node', shape='circle')

        mloads_id = 3
        # mloads_id = 100
        ase_id = 0
        # asecont_id = 100001
        if mloads_id == 0:
            subcases = self.model.subcases
            assert len(subcases) > 0, subcases
            print(subcases)
            subcase0 = subcases.pop(0)
            if 'MLOADS' in subcase0:
                mloads_id = subcase0['MLOADS'][0]
            elif 'ASE' in subcase0:
                ase_id = subcase0['ASE'][0]
        print(f'mloads_id = {mloads_id}')
        print(f'ase_id = {ase_id}')

        asecont = None
        if mloads_id in self.mloads:
            mloads = self.mloads[mloads_id]
            # print(mloads)
            mldcomd_id = mloads.mldcomd_id
            asecont = mloads.asecont_ref
        if ase_id in self.ase:
            ase = self.ase[ase_id]
            # print(ase)
            asecont = ase.asecont_ref

        #if mloads_id in self.mloads:
        if asecont is not None:
            #asecont_id = mloads.asecont_id
            tfset_id = asecont.tf_id
            gainset_id = asecont.gain_id
            cnctset_id = asecont.conct_id
            senset_id = asecont.sens_id
        else:
            tfset_id = 1000001
            cnctset_id = 1000001
            gainset_id = 1000001
            senset_id = 1000001
            mldcomd_id = 0
            # raise RuntimeError('mloads')
            # mldcomd_id = 401
        # tfset_id = 0
        # cnctset_id = 0
        # senset_id =0

        if mldcomd_id and mldcomd_id in self.mldcomd:
            g.attr('node', shape='box')
            mldcomd = self.mldcomd[mldcomd_id]
            assert mldcomd.extinps_ref is not None, mldcomd
            for extinp_ref in mldcomd.extinps_ref:
                # input_type: '2'
                # itf_component: 3
                # itf_id: 400006
                itf_ref = extinp_ref.itf_ref
                if itf_ref.type == 'CJUNCT':
                    # output_name = f'{itf_ref.type}={extinp_ref.itf_component}'
                    output_name = f'{itf_ref.type}={itf_ref.cjunct_id} (in={itf_ref.nu}, out={itf_ref.ny})'
                    assert itf_ref.ny == 1, itf_ref.get_stats()
                    valuei = itf_ref.values[extinp_ref.itf_component-1, 0]
                    # ki = itf_ref.
                    # nu: 3
                    # ny: 1
                    # values: [0.0, 0.0, 1.0]
                    output_comment = clean_comment(itf_ref.comment)
                else:
                    raise RuntimeError(itf_ref)

                input_name = f'EXTINP={extinp_ref.extinp_id}'
                input_comment = clean_comment(extinp_ref.comment)
                tag = f'({valuei}*I={extinp_ref.itf_component})'
                g.edge(input_name + input_comment,
                       output_name + output_comment,
                       label=f'MLDCOMD={mldcomd_id} {tag}')
                # ITFID
                # CI

        #---------------------------------------------
        tfset_ids = []
        conct_ids = []
        asegain_ids = []
        asesnsr_ids = []

        if tfset_id in self.tfset: # SISOTF/CJUNCT/MIMOSS
            tfset = self.tfset[tfset_id]
            tfset_ids = set(tfset.ids)

        if cnctset_id in self.cnctset: # CONCT
            cnctset = self.cnctset[cnctset_id]
            conct_ids = cnctset.ids

        if gainset_id in self.gainset: # ASEGAIN
            gainset = self.gainset[gainset_id]
            asegain_ids = gainset.ids

        if senset_id in self.senset: # ASESNSR
            senset = self.senset[senset_id]
            asesnsr_ids = senset.ids

        print(f'tfset_ids = {tfset_ids}')
        print(f'all_conct_ids = {conct_ids}')
        print(f'asegain_ids = {asegain_ids}')
        print(f'asesnsr_ids = {asesnsr_ids}')

        log = self.model.log
        
        # draw ASESNSRs
        for idi, card in self.asegain.items():
            if idi not in asegain_ids:
                continue
            output_ref = card.output_ref
            output_name = f'{output_ref.type}={output_ref.asesnsr_id} ({output_ref.name})'
            # print(output_ref.get_stats())

            # if output_ref.type in tfset_ids:
            output_comment = clean_comment(output_ref.comment)
            if output_ref.asesnsr_id in asesnsr_ids:
                g.attr('node', shape='ellipse')
                g.node(output_name+output_comment)
            else:
                g.attr('node', shape='box')
                g.node(output_name+output_comment)

        for actu_id, card in self.actu.items():
            name = f'{card.type}={actu_id}'
            comment = clean_comment(card.comment)
            g.node(name+comment)

        # g.attr('node', shape='diamond')
        # for cjunct_id, card in self.cjunct.items():
        #     name = f'{card.type}={cjunct_id} (in={card.nu}, out={card.ny})'
        #     comment = clean_comment(card.comment)
        #     g.node(name+comment)

        g.attr('node', shape='box')
        for idi, card in self.sisotf.items():
            if idi not in tfset_ids:
                continue
            name = f'{card.type}={card.sisotf_id}'
            comment = clean_comment(card.comment)
            g.node(name+comment)

        g.attr('node', shape='box')
        for idi, card in self.asegain.items():
            if idi not in asegain_ids:
                continue
            # c_in: 1
            # c_out: 1
            # gain: 4001
            # gain_type: 'Q'
            # itf_id: 400001
            # otf_id: 1096004
            input_ref = card.input_ref
            # if input_ref is None:
            #     log.warning(f'missing input-type for:\n{str(card)}')
            # elif input_ref.type != 'ASEGAIN':
            #     input_name = f'{input_ref.type}={input_ref.input_tf_id}'
            if input_ref.type == 'SISOTF':
                itag = '' if input_ref.sisotf_id in tfset_ids else 'x'
                assert itag == '', input_ref
                input_name = f'{itag}{input_ref.type}={input_ref.sisotf_id}'
                input_comment = clean_comment(input_ref.comment)
            else:
                raise RuntimeError(input_ref)

            output_ref = card.output_ref
            # if output_ref is None:
            #     log.warning(f'missing output-type for:\n{str(card)}')
                # asdf
            if output_ref.type == 'ASESNSR':
                otag = '' if output_ref.asesnsr_id in asesnsr_ids else 'x'
                output_name = f'{otag}{output_ref.type}={output_ref.asesnsr_id} ({output_ref.name})'
            # else:
            #     output_name = f'{output_ref.type}={output_ref.input_tf_id}'
            else:
                raise RuntimeError(output_ref)
            output_comment = clean_comment(output_ref.comment)
            # print(f'{output_comment!r}')
            # output_comment = ''

            tag = f'\n(I={card.c_in}, O={card.c_out})'

            # e.node('name1', label='name')
            g.edge(output_name+output_comment,
                   input_name+input_comment,
                   label=f'ASEGAIN={idi} {tag}')

        log = self.model.log
        g.attr('node', shape='box')
        for idi, card in self.conct.items():
            if idi not in conct_ids:
                continue

            # this is a CONCT
            # 'SISOTF=31004-1', 'ACTU=21001-1', 'ACTU=21002-1
            sivalue = ''
            input_ref = card.input_ref

            # CJUNCT, MIMOSS, SISOTF or ACTU
            if input_ref is None:
                log.warning(f'missing input-type for:\n{str(card)}')
                input_name = f'Input CONCT={idi}'
                input_comment = '\n???'
            elif input_ref.type == 'CJUNCT':
                sivalue = f'{input_ref.values[card.input_component-1,0]}*'
                input_name = f'{input_ref.type}={card.input_tf_id} (in={input_ref.nu}, out={input_ref.ny})'
                input_comment = clean_comment(input_ref.comment)
            elif input_ref.type == 'MIMOSS':
                itag = '' if card.input_tf_id in tfset_ids else 'x'
                input_name = f'{itag}{input_ref.type}={card.input_tf_id} (in={input_ref.nu}, out={input_ref.ny})'
                input_comment = clean_comment(input_ref.comment)

            elif input_ref.type == 'SISOTF':
                itag = '' if card.input_tf_id in tfset_ids else 'x'
                input_name = f'{itag}{input_ref.type}={card.input_tf_id}'
                input_comment = clean_comment(input_ref.comment)
            elif input_ref.type == 'ACTU':
                input_name = f'{input_ref.type}={card.input_tf_id}'
                input_comment = clean_comment(input_ref.comment)
            else:  # pragma: no cover
                raise RuntimeError(input_ref)

            log.debug(f'found {input_name}')
            # assert input_name in all_blocks, f'input={input_name!r} not in all_blocks\n{str(card)}'

            # this is a CONCT
            sovalue = ''
            output_ref = card.output_ref
            # CJUNCT, MIMOSS, SISOTF, ASESNSR, or ASESNS1
            if output_ref is None:
                log.warning(f'missing output-type for:\n{str(card)}')
                output_name = f'Output CONCT={idi}'
                output_comment = '\n???'
            elif output_ref.type == 'SISOTF':
                output_name = f'{output_ref.type}={card.output_tf_id}'
                output_comment = clean_comment(output_ref.comment)
            elif output_ref.type == 'CJUNCT':
                # print('output_ref.values', output_ref.values)
                outputs = output_ref.values[:, card.output_component-1].tolist()
                # if len(outputs) == 1:
                #     sovalue = f'{outputs[0]}*'
                # else:
                #     sovalue = f'{outputs}*'
                output_name = f'{output_ref.type}={card.output_tf_id} (in={output_ref.nu}, out={output_ref.ny})'
                output_comment = clean_comment(output_ref.comment)
            elif output_ref.type == 'MIMOSS':
                # print('output_ref.values', output_ref.values)
                #outputs = output_ref.values[:, card.output_component - 1].tolist()
                # if len(outputs) == 1:
                #     sovalue = f'{outputs[0]}*'
                # else:
                #     sovalue = f'{outputs}*'
                output_name = f'{output_ref.type}={card.output_tf_id} (in={output_ref.nu}, out={output_ref.ny})'
                output_comment = clean_comment(output_ref.comment)
            elif output_ref.type == 'ASESNSR':
                output_name = f'{output_ref.type}={card.output_tf_id}'
                output_comment = clean_comment(output_ref.comment)
            else:  # pragma: no cover
                raise RuntimeError(output_ref)
                output_name = f'{output_ref.type}={card.output_tf_id}'
                output_comment = '' if output_ref.type != 'ACTU' else clean_comment(output_ref.comment)
            log.debug(f'found {output_name}')

            tag = f'({sivalue}I={card.input_component}, {sovalue}O={card.output_component})'
            g.edge(output_name+output_comment,
                   input_name+input_comment,
                   label=f'CONCT={idi}\n{tag}')

        #-----------------
        # aeslinks
        for aeslink_label, aeslink in self.aeslink.items():
            actu_ref = aeslink.actu_ref
            output_name = f'{actu_ref.type}={actu_ref.actu_id}'
            output_comment = clean_comment(actu_ref.comment)

            # actu_id: int, independent_labels: list[str],
            # linking_coefficients: list[float],

            for label, label_ref, coeff in zip(aeslink.independent_labels,
                                               aeslink.independent_labels_ref,
                                               aeslink.linking_coefficients):
                # print(label_ref.get_stats())
                input_name = f'{label_ref.type}={label}'
                input_comment = clean_comment(label_ref.comment)
                g.edge(output_name+output_comment,
                       input_name+input_comment,
                       label=f'AESLINK={coeff}*{aeslink_label}')

        #-----------------
        try:
            g.view()
        except ExecutableNotFound:
            return

    def _check_tfset_cjunct(self):  # pragma: no cover
        assert len(self.tfset) == 1, self.tfset
        sisotfs = set(list(self.sisotf))
        cjunct_ids = set(list(self.cjunct))
        expected_tfs = cjunct_ids.union(sisotfs)
        for tfset_id, tfset in self.tfset.items():
            tfset_ids = set(tfset.ids)

        extra = expected_tfs - tfset_ids
        missing = tfset_ids - expected_tfs
        assert len(extra) == 0, f'There are more CJUNCTs than values in TFSET; extra={extra}'
        assert len(missing) == 0, f'There are fewer CJUNCTs than values in TFSET; missing={missing}'

    def _check_cntcset_conct(self):  # pragma: no cover
        # assert len(self.cnctset) in [0, 1], len(self.cnctset)
        cntcset_ids = set()
        for idi, cnctset in self.cnctset.items():
            cntcset_ids.update(cnctset.ids)
        contc_ids = set(list(self.conct))
        extra = contc_ids - cntcset_ids
        missing = cntcset_ids - contc_ids
        assert len(extra) == 0, f'There are more CONTCs than values in CNTCADD; extra={extra}'
        assert len(missing) == 0, f'There are fewer CONTCs than values in CNTCADD; missing={missing}'
        all_blocks = []
        # assert len(self.mimoss) == 0, self.mimoss
        for idi, card in self.sisotf.items():
            all_blocks.append(f'{card.type}={idi}-1')
        for idi, card in self.mimoss.items():
            all_blocks.append(f'{card.type}={idi}-1')
        for idi, card in self.actu.items():
            all_blocks.append(f'{card.type}={idi}-1')
        for idi, card in self.cjunct.items():
            # i = 0
            for i in range(1, card.nu+1):  # inputs
                all_blocks.append(f'{card.type}={idi}-{i}')
            # for j in range(i, card.ny+1):
            #     all_blocks.append(f'{card.type}={idi}-{j}')
            # nu: 2
            # ny: 1

        log = self.model.log
        print(f'all_blocks = {all_blocks}')
        for idi, card in self.conct.items():
            # 'SISOTF=31004-1', 'ACTU=21001-1', 'ACTU=21002-1
            # print(card)
            input_ref = card.input_ref
            if input_ref is None:
                log.warning(f'missing input-type for:\n{str(card)}')
                continue

            input_name = f'{input_ref.type}={card.input_tf_id}-{card.input_component}'
            log.debug(f'found {input_name}')
            assert input_name in all_blocks, f'input={input_name!r} not in all_blocks\n{str(card)}'

            output_ref = card.output_ref
            if output_ref is None:
                log.warning(f'missing output-type for:\n{str(card)}')
                continue
            output_name = f'{output_ref.type}={card.output_tf_id}-{card.output_component}'
            if output_name not in all_blocks:
                log.warning(f'output={output_name!r} not in all_blocks\n{str(card)}')
                continue
            log.debug(f'found {output_name}')
            # asdf
        # for
        # asdf

    def safe_cross_reference(self, xref_errors=None):
        model = self.model
        if model.nastran_format not in {'zona', 'zaero'}:
            return
        if xref_errors is None:
            xref_errors = defaultdict(list)

        for caero in model.caeros.values():
            caero.safe_cross_reference(model, xref_errors)
            self.caero_to_name_map[caero.label] = caero.eid

        singletons, dicts, dicts_list = get_dicts(self, 'xref')
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
        singletons, dicts, dicts_list = get_dicts(zaero, 'write')
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

    def write_bdf(self, bdf_file: TextIO, size: int=8,
                  is_double: bool=False,
                  sort_cards: bool=True):
        #if self.model.nastran_format != 'zona':
            #return
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

        singletons, dicts, dicts_list = get_dicts(self, 'write')
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

    def convert_to_nastran(self, save: bool=True):
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
        msg = '<ZAERO>; nPANLSTs=%s nmkaeroz=%s' % (
            len(self.panlsts), len(self.mkaeroz),
        )
        return msg

    def get_bdf_cards(self, bulk_data_lines: list[str],
                      bulk_data_ilines: Optional[np.ndarray]=None,
                      use_dict: bool=False) -> tuple[Any, Any, Any]:
        """Parses the BDF lines into a list of card_lines"""
        #self.log.warning('get_bdf_cards')
        use_dict = False

        dict_cards = {}
        model = self.model
        allow_tabs = False
        if bulk_data_ilines is None:
            bulk_data_ilines = np.zeros((len(bulk_data_lines), 2), dtype='int32')

        cards_list: list[Any] = []
        cards_dict: dict[str, list[Any]] = defaultdict(list)
        #cards = defaultdict(list)
        card_count: dict[str, int] = defaultdict(int)
        full_comment = ''
        card_lines = []
        card_ilines = []
        old_ifile_iline = None
        old_card_name = None
        backup_comment = ''
        nlines = len(bulk_data_lines)

        # self.echo = True
        # self.force_echo_off = False

        log = self.model.log
        for iline_bulk, line in enumerate(bulk_data_lines):
            ifile_iline = bulk_data_ilines[iline_bulk, :]
            # print(iline_bulk, ifile_iline)
            # print(iline_bulk, ifile_iline, line)
            # print('    backup={backup_comment!r}')
            comment = ''
            if '$' in line and (line.lstrip().startswith('$') or line.index('$') >= 72):
                line, comment = line.split('$', 1)
                # if line.strip():
                #     print(line)
                strip_comment = comment.strip()
                if strip_comment.lower().startswith('group:'):
                    continue_flag = model._store_group(strip_comment)
                    if continue_flag:
                        continue
            #if not self.allow_tabs and '\t' in line:
                #raise RuntimeError(f'There are tabs in:\n{line}')
                #self.log.warning(f'There are tabs in:\n{line}')

            card_name = line.split(',', 1)[0].split('\t', 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ['+', '*']:
                if old_card_name:
                    # multiline card is finished
                    if card_name not in ZAERO_CARDS:
                        pass
                    elif not allow_tabs and '\t' in (joined_lines_n := '\n'.join(card_lines)):
                        joined_lines_n2 = '\n'.join((f'{line!r}' for line in card_lines))
                        log.warning(f'There are tabs in:\n{joined_lines_n2}')
                        # raise RuntimeError(f'There are tabs in:\n{joined_lines_n2}')

                    if model.echo and not model.force_echo_off:
                        model.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    #if full_comment:
                        #print('full_comment = ', full_comment)
                    cards_list.append([old_card_name, _prep_comment(full_comment),
                                       card_lines, card_ilines[-1]])

                    card_count[old_card_name] += 1
                    card_lines = []
                    card_ilines = []
                    full_comment = ''

                    if old_card_name == 'ECHOON':
                        self.echo = True
                    elif old_card_name == 'ECHOOFF':
                        self.echo = False
                old_ifile_iline = ifile_iline
                old_card_name = card_name.rstrip(' *')

                if old_card_name == 'ENDDATA':
                    model.card_count['ENDDATA'] = 1
                    if nlines - iline_bulk > 1:
                        nleftover = nlines - iline_bulk - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        model.log.debug(msg)
                    cards_list, cards_dict, card_count = fix_card_list(
                        cards_list, cards_dict, card_count)
                    return cards_list, cards_dict, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)

            #TODO: these additional \n need to be there for rejected cards
            #      but not parsed cards
            if line.rstrip():
                card_lines.append(line)
                card_ilines.append(ifile_iline)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += comment + '\n'
            #elif comment:
                #backup_comment += '$' + comment + '\n'

        if card_lines:
            if not allow_tabs and '\t' in (joined_lines_n := '\n'.join(card_lines)):
                log.error(f'There are tabs in:\n{joined_lines_n}')
                # raise RuntimeError(f'There are tabs in:\n{joined_lines_n}')

            if model.echo and not model.force_echo_off:
                model.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))
            #print('end_add %s' % card_lines)

            # old dictionary version
            #cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            #if backup_comment + full_comment:
                #print('backup_comment + full_comment = ', backup_comment + full_comment)
            if old_card_name in dict_cards:
                cards_dict[old_card_name].append([_prep_comment(
                    backup_comment + full_comment), card_lines, ifile_iline])
            else:
                # cards_list.append([old_card_name, _prep_comment(
                #     backup_comment + full_comment), card_lines, ifile_iline])
                cards_list.append([old_card_name, _prep_comment(
                    backup_comment + full_comment), card_lines, card_ilines[-1]])
            card_count[old_card_name] += 1
        self.echo = False

        cards_list, cards_dict, card_count = fix_card_list(
            cards_list, cards_dict, card_count)
        return cards_list, cards_dict, card_count


def fix_card_list(cards_list, cards_dict, card_count):
    assert len(cards_dict) == 0, cards_dict
    card_count = defaultdict(int)
    # skip_cards = {'SPLINE1', 'SPLINE2', 'SPLINE3', 'AEFACT', 'CONM2',
    #               'GENEL', 'DMI', 'DMIG', 'TABLED1', 'TABLED2',
    #               'EIGR', 'EIGRL', 'DAREA'}
    include_cards = {
        'GRID', 'SPOINT', 'EPOINT',
        'CORD2R', 'CORD2C', 'CORD2S',
        'CORD1R', 'CORD1C', 'CORD1S',
        'CROD', 'PROD', 'CONROD',
        'CBUSH', 'PBUSH',
        'CBUSH1D', 'PBUSH1D',
        'CBAR', 'PBAR', 'PBARL',
        'CBEAM', 'PBEAM', 'PBEAML',
        'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
        'PSOLID', 'PLSOLID',
        'CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6',
        'CTRIAR', 'CQUADR', 'CQUAD',
        'PSHELL', 'PCOMP', 'PCOMPG',
        'MAT1', 'MAT2', 'MAT8', 'MAT9', 'MAT10',
    }
    cards_list2 = []
    for card_name, comment, card_lines, ifile_iline in cards_list:
        # (ifile, iline) = ifile_iline
        if ifile_iline[0] == 1000 and card_name not in include_cards: # f06/prt
            continue
        card_count[card_name] += 1
        cards_list2.append((card_name, comment, card_lines, ifile_iline))
    return cards_list2, cards_dict, dict(card_count)


def get_dicts(zaero: ZAERO, method: str) -> tuple[list,
                                                  dict[int, list],
                                                  list[dict]]:
    assert method in ['xref', 'write'], f'method={method!r}'
    dicts = [
        # --------------general-------------
        # zaero.aeroz,
        zaero.atmos, zaero.flutter_table,
        # -------------geometry-------------
        # zaero.panlsts,  # special-list
        zaero.pafoil,
        zaero.attach,
        zaero.pltsurf, zaero.pltmode, zaero.pltmist,
        zaero.pltbode,
        # -------------transient------------
        zaero.mloads, zaero.eloads,
        # --------------flutter-------------
        zaero.nlfltr, zaero.mkaeroz,
        # ---------------trim---------------
        # zaero.trim,  # part of the main BDF
        zaero.aeslink, zaero.trimvar, zaero.trimlnk,
        zaero.trimfnc, zaero.trimobj, zaero.trimcon,
        zaero.actu,
        # ---------------ase---------------
        zaero.cjunct, zaero.conct, zaero.tfset, zaero.cnctset,
        zaero.ase, zaero.asecont, zaero.asesnsr, zaero.asesns1,
        zaero.asegain, zaero.gainset,
        zaero.mimoss, zaero.sisotf,
        #
        zaero.senset, zaero.surfset,
        zaero.mldtrim, zaero.mldstat, zaero.minstat, zaero.mldprnt,
        zaero.mldcomd, zaero.mldtime,
        # zaero.extinp, zaero.extout,
        zaero.loadmod, zaero.rbred,
        zaero.aerolag,
        # ---------------gust---------------
        zaero.gloads, zaero.dgust, zaero.cgust,
        # ---------------other--------------
        zaero.extfile,
        zaero.dse,
        zaero.dmil,
        # plotting
        zaero.pltvg, zaero.pltbode,
        zaero.pltaero,
    ]
    dict_lists = [
        zaero.pltcp, zaero.pltflut,
        zaero.plttime, zaero.plttrim,
    ]
    if method == 'write':
        # these are xref'd by their parent
        dicts.extend([zaero.setadd, zaero.extinp, zaero.extout])
    singletons = [zaero.splinem]
    return singletons, dicts, dict_lists


def clean_comment(comment: str) -> str:
    lines = comment.split('\n')
    lines2 = []
    for line in lines:
        commenti = line.strip('$ -')
        if commenti.startswith('#'):
            continue
        if commenti:
            lines2.append(commenti)
    comment2 = '\n'.join(lines2).replace('\r', '').replace(':', '-')
    if comment2:
        return '\n' + comment2
    return comment2
