# coding: utf-8
# pylint: disable=W0212,C0103
from __future__ import annotations
from collections import defaultdict
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils import object_attributes, object_methods

from typing import Any
from pyNastran.bdf.bdf_interface.utils import _prep_comment
from pyNastran.bdf.bdf_interface.pybdf import _clean_comment

from pyNastran.bdf.cards.aero.zona_cards.atm import (
    ATMOS, FIXMATM, FIXHATM)
from pyNastran.bdf.cards.aero.zona_cards.spline import (
    SPLINE1_ZONA, SPLINE2_ZONA, SPLINE3_ZONA,)
from pyNastran.bdf.cards.aero.zona_cards.geometry import (
    PANLST1, PANLST2, PANLST3, SEGMESH,
    CAERO7, BODY7, PAFOIL7, AESURFZ, AESLINK)
from pyNastran.bdf.cards.aero.zona_cards.plot import (
    PLTAERO, PLTMODE, )
from pyNastran.bdf.cards.aero.zona_cards.flutter import (
    FLUTTER_ZONA, MKAEROZ)
from pyNastran.bdf.cards.aero.zona_cards.trim import (
    TRIM_ZONA, TRIMVAR, TRIMLNK,)
from pyNastran.bdf.cards.aero.zona_cards.manuever import (
    MLOADS, EXTINP, EXTOUT, TRIMFNC, ACTU, LOADMOD, RBRED,)
from pyNastran.bdf.cards.aero.zona_cards.gust import (
    GLOADS, DGUST, CGUST, MFTGUST)
from pyNastran.bdf.cards.aero.zona_cards.ase import (
    ASE, ASECONT, ASESNSR, ASESNS1,
    CJUNCT, CONCT, TFSET, MIMOSS, SISOTF,
    SENSET, SURFSET, CNCTSET,
    ASEGAIN, GAINSET,
)
from pyNastran.bdf.cards.aero.zona_cards.bdf_tables import (
    TABLED1_ZONA, TABDMP1_ZONA)
from pyNastran.bdf.cards.aero.zona_cards.dmi import DMIL
from pyNastran.bdf.cards.aero.zona_cards.cards import (
    MLDPRNT, MLDSTAT, MINSTAT, MLDTRIM, MLDCOMD, MLDTIME,
    AEROZ, ACOORD, ATTACH, EXTFILE,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

ZONA_CARDS = [
    # atmosphere
    'FIXHATM', 'FIXMATM',  # flutter table
    'ATMOS',

    # geometry
    'CAERO7', 'AEROZ', 'AESURFZ', 'AESLINK',
    'ATTACH',
    'PANLST1', 'PANLST2', 'PANLST3', 'PAFOIL7',
    'SEGMESH', 'BODY7', 'ACOORD', 'MKAEROZ',
    'TRIMVAR', 'TRIMLNK',
    # -------------
    # flutter
    'FLUTTER',
    'FIXMDEN', 'FIXHATM', 'FIXMACH',
    # -------------
    # plotting
    'PLTMODE', 'PLTAERO', 'PLTTIME',
    # -------------
    # mloads
    'MLOADS',
    'CJUNCT', 'CONCT', 'TFSET',
    'MLDSTAT', 'MLDTRIM', 'MLDPRNT', 'MLDCOMD',
    'EXTINP', 'EXTOUT', 'TRIMFNC', 'LOADMOD',
    # -------------
    # gust
    'GLOADS', 'DGUST', 'CGUST',
    # -------------
    # ase
    'ASE', 'ASECONT',
    'ASESNSR', 'ASESNS1', 'SENSET',
    'AEROLAG', 'ACTU',
    'MIMOSS', 'SISOTF',
    'ASEGAIN', 'GAINSET',
    # -------------
    # other
    'DMIL', 'EXTFILE',
    'MLDTIME', 'MLDCOMD',
    'MINSTAT', 'APCONST',
    # 'PLTFLUT', 'PLTVG', 'PLTCP', 'PLTMIST',
    'RBRED',
    # 'SPLINE0', 'PBODY7',
    'CNCTSET', 'SURFSET',
]


class AddMethods:
    def __init__(self, model: BDF):
        self.model = model

    @property
    def zona(self):
        return self.model.zona

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
        assert key not in self.model.zona.mldcomd, key
        assert key > 0, key
        self.model.zona.mldcomd[key] = mldcomd
        self.model._type_to_id_map[mldcomd.type].append(key)

    def add_extfile_object(self, extfile: EXTFILE) -> None:
        """adds an EXTFILE object"""
        key = extfile.extfile_id
        assert key not in self.model.zona.extfile, key
        assert key > 0, key
        self.model.zona.extfile[key] = extfile
        self.model._type_to_id_map[extfile.type].append(key)

    def add_dmil_object(self, dmil: DMIL) -> None:
        """adds an DMIL object"""
        name = dmil.name
        key = (name, dmil.row, dmil.col)
        assert key not in self.model.zona.dmil, key
        self.model.zona.dmil[key] = dmil
        self.model._type_to_id_map[dmil.type].append(name)

    def add_mldtime_object(self, mldtime: MLDTIME) -> None:
        """adds an MLDTRIM object"""
        key = mldtime.mldtime_id
        assert key not in self.model.zona.mldtime, key
        assert key > 0, key
        self.model.zona.mldtime[key] = mldtime
        self.model._type_to_id_map[mldtime.type].append(key)

    def add_mldprnt_object(self, mldprnt: MLDPRNT) -> None:
        """adds an MLDPRNT object"""
        key = mldprnt.mldtime_id
        assert key not in self.model.zona.mldtime, key
        assert key > 0, key
        self.model.zona.mldprnt[key] = mldprnt
        self.model._type_to_id_map[mldprnt.type].append(key)

    def add_mldtrim_object(self, mldtrim: MLDTRIM) -> None:
        """adds an MLDTRIM object"""
        key = mldtrim.mldtrim_id
        assert key not in self.model.zona.mldtrim, key
        assert key > 0, key
        self.model.zona.mldtrim[key] = mldtrim
        self.model._type_to_id_map[mldtrim.type].append(key)

    def add_mldstat_object(self, mldstat: MLDSTAT) -> None:
        """adds an MLDSTAT object"""
        key = mldstat.mldstat_id
        assert key > 0, key
        zona = self.model.zona
        assert key not in zona.mldstat, key
        zona.mldstat[key] = mldstat
        self.model._type_to_id_map[mldstat.type].append(key)

    def add_minstat_object(self, minstat: MINSTAT) -> None:
        """adds an MINSTAT object"""
        key = minstat.minstat_id
        assert key not in self.model.zona.minstat, key
        assert key > 0, key
        self.model.zona.minstat[key] = minstat
        self.model._type_to_id_map[minstat.type].append(key)

    def add_extinp_object(self, extinp: EXTINP) -> None:
        """adds an EXTINP object"""
        key = extinp.extinp_id
        assert key not in self.model.zona.extinp, key
        assert key > 0, key
        self.model.zona.extinp[key] = extinp
        self.model._type_to_id_map[extinp.type].append(key)

    def add_extout_object(self, extout: EXTOUT) -> None:
        """adds an EXTOUT object"""
        key = extout.extout_id
        assert key not in self.model.zona.extout, key
        assert key > 0, key
        self.model.zona.extout[key] = extout
        self.model._type_to_id_map[extout.type].append(key)

    def add_trimfnc_object(self, trimfnc: TRIMFNC) -> None:
        """adds an TRIMFNC object"""
        key = trimfnc.trimfnc_id
        assert key not in self.model.zona.trimfnc, key
        assert key > 0, key
        self.model.zona.trimfnc[key] = trimfnc
        self.model._type_to_id_map[trimfnc.type].append(key)

    def add_mimoss_object(self, mimoss: MIMOSS) -> None:
        """adds an MIMOSS object"""
        key = mimoss.mimoss_id
        assert key not in self.model.zona.mimoss, key
        assert key > 0, key
        self.model.zona.mimoss[key] = mimoss
        self.model._type_to_id_map[mimoss.type].append(key)


    def add_sisotf_object(self, sisotf: SISOTF) -> None:
        """adds an SISOTF object"""
        key = sisotf.sisotf_id
        assert key not in self.model.zona.sisotf, key
        assert key > 0, key
        self.model.zona.sisotf[key] = sisotf
        self.model._type_to_id_map[sisotf.type].append(key)

    def add_tfset_object(self, tfset: TFSET) -> None:
        """adds an TFSET object"""
        key = tfset.tfset_id
        assert key not in self.model.zona.tfset, key
        assert key > 0, key
        self.model.zona.tfset[key] = tfset
        self.model._type_to_id_map[tfset.type].append(key)

    def add_senset_object(self, senset: SENSET) -> None:
        """adds an SENSET object"""
        key = senset.senset_id
        assert key not in self.model.zona.senset, key
        assert key > 0, key
        self.model.zona.senset[key] = senset
        self.model._type_to_id_map[senset.type].append(key)

    def add_cnctset_object(self, cnctset: CNCTSET) -> None:
        """adds an CNCTSET object"""
        key = cnctset.cnctset_id
        assert key not in self.model.zona.cnctset, key
        assert key > 0, key
        self.model.zona.cnctset[key] = cnctset
        self.model._type_to_id_map[cnctset.type].append(key)

    def add_surfset_object(self, surfset: SURFSET) -> None:
        """adds an SURFSET object"""
        key = surfset.surfset_id
        assert key not in self.model.zona.surfset, key
        assert key > 0, key
        self.model.zona.surfset[key] = surfset
        self.model._type_to_id_map[surfset.type].append(key)

    def add_loadmod_object(self, loadmod: LOADMOD) -> None:
        """adds an LOADMOD object"""
        key = loadmod.loadmod_id
        assert key not in self.model.zona.loadmod, key
        assert key > 0, key
        self.model.zona.loadmod[key] = loadmod
        self.model._type_to_id_map[loadmod.type].append(key)

    def add_rbred_object(self, rbred: RBRED) -> None:
        """adds an LOADMOD object"""
        key = rbred.sid
        assert key not in self.model.zona.rbred, key
        assert key > 0, key
        self.model.zona.rbred[key] = rbred
        self.model._type_to_id_map[rbred.type].append(key)

    def add_panlst_object(self, panlst: PANLST1 | PANLST2 | PANLST3) -> None:
        """adds an PANLST1/PANLST2/PANLST3 object"""
        key = panlst.eid
        assert key > 0, key
        zona = self.model.zona
        # assert key not in zona.panlsts, '\npanlst=\n%s old=\n%s' % (
        #     panlst, zona.panlsts[key])
        if key not in zona.panlsts:
            zona.panlsts[key] = []
        zona.panlsts[key].append(panlst)
        self.model._type_to_id_map[panlst.type].append(key)

    def add_pafoil_object(self, pafoil: PAFOIL7) -> None:
        """adds an PAFOIL7/PAFOIL8 object"""
        key = pafoil.pid
        assert pafoil.pid > 0
        zona = self.model.zona
        assert key not in zona.pafoil, '\npafoil7=\n%s old=\n%s' % (
            pafoil, zona.pafoil[key])
        zona.pafoil[key] = pafoil
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
        zona = model.zona
        assert key not in zona.aeslink, '\naeslink=\n%s old=\n%s' % (
            aeslink, zona.aeslink[key])
        zona.aeslink[key] = aeslink
        model._type_to_id_map[aeslink.type].append(key)

    def add_mloads_object(self, mloads: MLOADS) -> None:
        """adds an MLOADS object"""
        key = mloads.mloads_id
        model = self.model
        zona = model.zona
        assert key not in zona.mloads, '\nmloads=\n%s old=\n%s' % (
            mloads, zona.mloads[key])
        zona.mloads[key] = mloads
        model._type_to_id_map[mloads.type].append(key)

    def add_gloads_object(self, gloads: GLOADS) -> None:
        """adds an GLOADS object"""
        key = gloads.gloads_id
        model = self.model
        zona = model.zona
        assert key not in zona.gloads, '\ngloads=\n%s old=\n%s' % (
            gloads, model.gloads[key])
        zona.gloads[key] = gloads
        model._type_to_id_map[gloads.type].append(key)

    def add_dgust_object(self, dgust: DGUST) -> None:
        """adds an DGUST object"""
        key = dgust.dgust_id
        model = self.model
        zona = model.zona
        assert key not in zona.dgust, '\ndgust=\n%s old=\n%s' % (
            dgust, zona.dgust[key])
        zona.dgust[key] = dgust
        model._type_to_id_map[dgust.type].append(key)

    def add_cgust_object(self, cgust: CGUST) -> None:
        """adds an CGUST object"""
        key = cgust.cgust_id
        model = self.model
        zona = model.zona
        assert key not in zona.cgust, '\ncgust=\n%s old=\n%s' % (
            cgust, zona.cgust[key])
        zona.cgust[key] = cgust
        model._type_to_id_map[cgust.type].append(key)

    def add_ase_object(self, ase: ASE) -> None:
        """adds an ASE object"""
        key = ase.ase_id
        model = self.model
        zona = model.zona
        assert key not in zona.ase, '\nase=\n%s old=\n%s' % (
            ase, zona.ase[key])
        zona.ase[key] = ase
        model._type_to_id_map[ase.type].append(key)

    def add_asecont_object(self, asecont: ASECONT) -> None:
        """adds an ASECONT object"""
        key = asecont.asecont_id
        model = self.model
        zona = model.zona
        assert key not in zona.asecont, '\nasecont=\n%s old=\n%s' % (
            asecont, zona.asecont[key])
        zona.asecont[key] = asecont
        model._type_to_id_map[asecont.type].append(key)

    def add_asegain_object(self, asegain: ASEGAIN) -> None:
        """adds an ASEGAIN object"""
        key = asegain.asegain_id
        model = self.model
        zona = model.zona
        assert key not in zona.asegain, '\nasegain=\n%s old=\n%s' % (
            asegain, zona.asegain[key])
        zona.asecont[key] = asegain
        model._type_to_id_map[asegain.type].append(key)

    def add_asesnsr_object(self, asesnsr: ASESNSR) -> None:
        """adds an ASESNSR object"""
        key = asesnsr.asesnsr_id
        model = self.model
        zona = model.zona
        assert key not in zona.asesnsr, '\nasesnsr=\n%s old=\n%s' % (
            asesnsr, zona.asesnsr[key])
        zona.asesnsr[key] = asesnsr
        model._type_to_id_map[asesnsr.type].append(key)

    def add_asesns1_object(self, asesns1: ASESNS1) -> None:
        """adds an ASESNS1 object"""
        key = asesns1.asesns1_id
        model = self.model
        zona = model.zona
        assert key not in zona.asesns1, '\nasesns1=\n%s old=\n%s' % (
            asesns1, zona.asesns1[key])
        zona.asesns1[key] = asesns1
        model._type_to_id_map[asesns1.type].append(key)

    def add_gainset_object(self, gainset: GAINSET) -> None:
        """adds an ASEGAIN object"""
        key = gainset.asegain_id
        model = self.model
        zona = model.zona
        assert key not in zona.gainset, '\ngainset=\n%s old=\n%s' % (
            gainset, zona.gainset[key])
        zona.asecont[key] = gainset
        model._type_to_id_map[gainset.type].append(key)

    def add_cjunct_object(self, cjunct: CJUNCT) -> None:
        """adds an CJUNCT object"""
        key = cjunct.cjunct_id
        model = self.model
        zona = model.zona
        assert key not in zona.cjunct, '\ncjunct=\n%s old=\n%s' % (
            cjunct, zona.cjunct[key])
        zona.cjunct[key] = cjunct
        model._type_to_id_map[cjunct.type].append(key)

    def add_conct_object(self, conct: CONCT) -> None:
        """adds an CONCT object"""
        key = conct.cjunct_id
        model = self.model
        zona = model.zona
        assert key not in zona.conct, '\nconct=\n%s old=\n%s' % (
            conct, zona.conct[key])
        zona.conct[key] = conct
        model._type_to_id_map[conct.type].append(key)

    def add_actu_object(self, actu: ACTU) -> None:
        """adds an AESURFZ object"""
        key = actu.actu_id
        model = self.model
        zona = model.zona
        assert key not in zona.actu, '\nactu=\n%s old=\n%s' % (
            actu, zona.actu[key])
        zona.actu[key] = actu
        model._type_to_id_map[actu.type].append(key)

    def add_mkaeroz_object(self, mkaeroz: MKAEROZ) -> None:
        """adds an MKAEROZ object"""
        assert mkaeroz.sid not in self.model.zona.mkaeroz
        assert mkaeroz.sid > 0
        key = mkaeroz.sid
        self.model.zona.mkaeroz[key] = mkaeroz
        self.model._type_to_id_map[mkaeroz.type].append(key)

    def add_trimvar_object(self, trimvar: TRIMVAR) -> None:
        """adds an TRIMVAR object"""
        assert trimvar.var_id not in self.model.zona.trimvar
        assert trimvar.var_id > 0
        key = trimvar.var_id
        self.model.zona.trimvar[key] = trimvar
        self.model._type_to_id_map[trimvar.type].append(key)

    def add_trimlnk_object(self, trimlnk: TRIMLNK) -> None:
        """adds an TRIMLNK object"""
        assert trimlnk.link_id not in self.model.zona.trimlnk
        assert trimlnk.link_id > 0
        key = trimlnk.link_id
        self.model.zona.trimlnk[key] = trimlnk
        self.model._type_to_id_map[trimlnk.type].append(key)

    def add_attach_object(self, attach: ATTACH) -> None:
        """adds an ATTACH object"""
        assert attach.attach_id not in self.model.zona.attach
        assert attach.attach_id > 0
        key = attach.attach_id
        self.model.zona.attach[key] = attach
        self.model._type_to_id_map[attach.type].append(key)

    def add_plotmode_object(self, plot: PLTMODE) -> None:
        """adds an PLTMODE object"""
        assert plot.set_id not in self.model.zona.plotmode, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zona.plotmode[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_plotaero_object(self, plot: PLTAERO) -> None:
        """adds an PLTAERO object"""
        assert plot.set_id not in self.model.zona.plotaero, str(plot)
        assert plot.set_id > 0
        key = plot.set_id
        self.model.zona.plotaero[key] = plot
        self.model._type_to_id_map[plot.type].append(key)

    def add_flutter_table_object(self, flutter_table: FIXHATM | FIXMATM) -> None:
        """adds an FIXMATM object"""
        key = flutter_table.sid
        model = self.model
        zona = model.zona
        assert key not in zona.flutter_table, '\nflutter_table=\n%s old=\n%s' % (
            flutter_table, zona.flutter_table[key])
        zona.flutter_table[key] = flutter_table
        model._type_to_id_map[flutter_table.type].append(key)

    def add_atmos_object(self, atmos: ATMOS) -> None:
        """adds an ATMOS object"""
        key = atmos.atmos_id
        model = self.model
        zona = model.zona
        assert key not in zona.atmos, '\natmos=\n%s old=\n%s' % (
            atmos, zona.atm[key])
        zona.atmos[key] = atmos
        model._type_to_id_map[atmos.type].append(key)

class ZONA:
    def __init__(self, model):
        self.model = model
        self.caero_to_name_map = {}
        self._add_methods = AddMethods(model)

        # aero models
        self.atmos: dict[int, ATMOS] = {}
        self.flutter_table: dict[int, FIXHATM | FIXMATM] = {}

        #: store PANLST1,PANLST2,PANLST3
        self.panlsts: dict[int, PANLST1 | PANLST2 | PANLST3] = {}
        self.attach: dict[int, PLTAERO] = {}
        self.plotaero: dict[int, PLTAERO] = {}
        self.plotmode: dict[int, PLTMODE] = {}

        #: store PAFOIL7/PAFOIL8
        self.pafoil: dict[int, PAFOIL7] = {}
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

        # trim
        self.aeslink: dict[int, AESLINK] = {}
        self.trimvar: dict[int, TRIMVAR] = {}
        self.trimlnk: dict[int, TRIMLNK] = {}
        self.mldtrim: dict[int, MLDTRIM] = {}

        # flutter
        self.mkaeroz: dict[int, MKAEROZ] = {}
        self.nlfltr: dict[int, NLFLTR] = {}

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
        self.aseout: dict[int, ASEOUT] = {}
        self.apcnsnd: dict[int, APCNSND] = {}
        self.apcnscp: dict[int, APCNSCP] = {}
        self.mimoss: dict[int, MIMOSS] = {}
        self.sisotf: dict[int, SISOTF] = {}
        self.cmargin: dict[int, CMARGIN] = {}

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
        """clears out the ZONA object"""
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
        if self.model.nastran_format != 'zona':
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
        if self.model.nastran_format != 'zona':
            return
        for panlst in self.panlsts.values():
            for panlsti in panlst:
                panlsti.validate()

        dicts = get_dicts(self, 'write')
        for items in dicts:
            for item in items.values():
                # self.model.log.info(f'xref {item.type}')
                item.validate()

    def PAFOIL(self, pid, msg=''):
        """gets a pafoil profile (PAFOIL7/PAFOIL8)"""
        try:
            return self.pafoil[pid]
        except KeyError:
            pafoils = np.unique(list(self.pafoil.keys()))
            raise KeyError(f'pid={pid} not found{msg}.  Allowed pafoils={pafoils}')

    def update_for_zona(self):
        """updates for zona"""
        card_parser = self.model._card_parser
        add_methods = self.model._add_methods
        zona_add = self._add_methods
        card_parser2 = {
            # aero models
            'ATMOS': (ATMOS, zona_add.add_atmos_object),
            'FIXMATM': (FIXMATM, zona_add.add_flutter_table_object),
            'FIXHATM': (FIXHATM, zona_add.add_flutter_table_object),
            # trim
            'TRIM': (TRIM_ZONA, add_methods.add_trim_object),
            'TABLED1': (TABLED1_ZONA, add_methods.add_tabled_object),
            'TABDMP1': (TABDMP1_ZONA, add_methods.add_table_sdamping_object),
            'AEROZ': (AEROZ, add_methods.add_aeros_object),
            # geometry
            'SPLINE1': (SPLINE1_ZONA, add_methods.add_spline_object),
            'SPLINE2': (SPLINE2_ZONA, add_methods.add_spline_object),
            'SPLINE3': (SPLINE3_ZONA, add_methods.add_spline_object),
            'PANLST1': (PANLST1, zona_add.add_panlst_object),
            'PANLST2': (PANLST2, zona_add.add_panlst_object),
            'PANLST3': (PANLST3, zona_add.add_panlst_object),
            'PAFOIL7': (PAFOIL7, zona_add.add_pafoil_object),
            'MKAEROZ': (MKAEROZ, zona_add.add_mkaeroz_object),
            'SEGMESH': (SEGMESH, add_methods.add_paero_object),
            'BODY7': (BODY7, add_methods.add_caero_object),
            'CAERO7': (CAERO7, add_methods.add_caero_object),
            'ACOORD': (ACOORD, add_methods.add_coord_object),
            'AESURFZ': (AESURFZ, zona_add.add_aesurfz_object),
            'AESLINK': (AESLINK, zona_add.add_aeslink_object),
            # flutter
            'FLUTTER': (FLUTTER_ZONA, add_methods.add_flutter_object),
            # mloads
            'MLOADS': (MLOADS, zona_add.add_mloads_object),
            # gust
            'GLOADS': (GLOADS, zona_add.add_gloads_object),
            'DGUST': (DGUST, zona_add.add_dgust_object),
            'CGUST': (CGUST, zona_add.add_cgust_object),
            # ase
            'ASE': (ASE, zona_add.add_ase_object),
            'ASECONT': (ASECONT, zona_add.add_asecont_object),
            'ASESNSR': (ASESNSR, zona_add.add_asesnsr_object),
            'ASESNS1': (ASESNS1, zona_add.add_asesns1_object),
            'MIMOSS': (MIMOSS, zona_add.add_mimoss_object),
            'SISOTF': (SISOTF, zona_add.add_sisotf_object),
            'CJUNCT': (CJUNCT, zona_add.add_cjunct_object),
            'CONCT': (CONCT, zona_add.add_conct_object),
            # other
            'SENSET': (SENSET, zona_add.add_senset_object),
            'CNCTSET': (CNCTSET, zona_add.add_cnctset_object),
            'SURFSET': (SURFSET, zona_add.add_surfset_object),
            'ACTU': (ACTU, zona_add.add_actu_object),
            'LOADMOD': (LOADMOD, zona_add.add_loadmod_object),
            'RBRED': (RBRED, zona_add.add_rbred_object),
            'TRIMVAR': (TRIMVAR, zona_add.add_trimvar_object),
            'TRIMLNK': (TRIMLNK, zona_add.add_trimlnk_object),
            'ATTACH': (ATTACH, zona_add.add_attach_object),
            'PLTMODE': (PLTMODE, zona_add.add_plotmode_object),
            'PLTAERO': (PLTAERO, zona_add.add_plotaero_object),
            'EXTINP': (EXTINP, zona_add.add_extinp_object),
            'EXTOUT': (EXTOUT, zona_add.add_extout_object),
            'TFSET': (TFSET, zona_add.add_tfset_object),
            'TRIMFNC': (TRIMFNC, zona_add.add_trimfnc_object),
            'MLDSTAT': (MLDSTAT, zona_add.add_mldstat_object),
            'MINSTAT': (MINSTAT, zona_add.add_minstat_object),
            'MLDTRIM': (MLDTRIM, zona_add.add_mldtrim_object),
            'MLDCOMD': (MLDCOMD, zona_add.add_mldcomd_object),
            'MLDTIME': (MLDTIME, zona_add.add_mldtime_object),
            'DMIL': (DMIL, zona_add.add_dmil_object),
            'EXTFILE': (EXTFILE, zona_add.add_extfile_object),
            'MLDPRNT': (MLDPRNT, zona_add.add_mldprnt_object),
            # PLTCP
            # PLTFLUT
            # PLTVG
            # PLTMIST
        }
        card_parser.update(card_parser2)
        self.model.cards_to_read.update(set(ZONA_CARDS))

    def cross_reference(self):
        if self.model.nastran_format != 'zona':
            return
        dicts = get_dicts(self, 'xref')
        for items in dicts:
            for item in items.values():
                # self.model.log.info(f'xref {item.type}')
                item.cross_reference(self.model)

        for caero in self.model.caeros.values():
            #print('%s uses CAERO eid=%s' % (caero.label, caero.eid))
            self.caero_to_name_map[caero.label] = caero.eid

    def safe_cross_reference(self, xref_errors=None):
        if self.model.nastran_format != 'zona':
            return
        if xref_errors is None:
            xref_errors = defaultdict(list)
        dicts = get_dicts(self, 'xref')
        for items in dicts:
            for item in items.values():
                # self.model.log.info(f'xref {item.type}')
                item.safe_cross_reference(self.model, {})

        for caero in self.model.caeros.values():
            #print('%s uses CAERO eid=%s' % (caero.label, caero.eid))
            self.caero_to_name_map[caero.label] = caero.eid

    def write_bdf(self, bdf_file: TextIO, size: int=8,
                  is_double: bool=False):
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

        dicts = get_dicts(self, 'write')
        for items in dicts:
            for key, value in items.items():
                bdf_file.write(value.write_card(size=size, is_double=is_double))

    def convert_to_nastran(self, save=True):
        """Converts a ZONA model to Nastran"""
        if self.model.nastran_format != 'zona':
            caeros = {}
            caero2s = []
            make_paero1 = False
            return caeros, caero2s, make_paero1

        caeros, caero2s, make_paero1 = self._convert_caeros()
        splines = self._convert_splines()
        aesurf, aelists = self._convert_aesurf_aelist()

        trims = self._convert_trim()
        aeros, aero = self.model.aeros.convert_to_zona(self.model)

        aelinks = self._convert_trimlnk()

        if save:
            self.clear()
            self.model.splines = splines
            self.model.aesurf = aesurf
            self.model.aelists = aelists
            self.model.aelinks = aelinks
            self.model.trims = trims
            self.model.aeros = aeros
            self.model.aero = aero
        return caeros, caero2s, make_paero1

    def _convert_caeros(self):
        """Converts ZONA CAERO7/BODY7 to CAERO1/CAERO2"""
        model = self.model
        caeros = {}
        caero2s = []
        make_paero1 = False
        for caero_id, caero in sorted(model.caeros.items()):
            if caero.type == 'CAERO7':
                caero_new = caero.convert_to_nastran()
                make_paero1 = True
            elif caero.type == 'BODY7':
                caero2s.append(caero)
                continue
            else:
                raise NotImplementedError(caero)
            caeros[caero_id] = caero_new

        self.add_caero2s(caero2s, add=False)
        return caeros, caero2s, make_paero1

    def add_caero2s(self, caero2s, add=False):
        """Converts ZONA BODY7 to CAERO2/PAERO2/AEFACT"""
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

    def _convert_splines(self):
        """Converts ZONA splines to splines"""
        splines = {}
        for unused_spline_id, spline in self.model.splines.items():
            #print(spline)
            if spline.type == 'SPLINE1_ZONA':
                splines_new = spline.convert_to_nastran(self.model)
            elif spline.type == 'SPLINE3_ZONA':
                splines_new = spline.convert_to_nastran(self.model)
            else:
                raise NotImplementedError(spline)
            for spline_new in splines_new:
                splines[spline.eid] = spline_new
        return splines

    def _convert_aesurf_aelist(self):
        """
        Converts ZONA AESURFZ to AESURF/AELIST

        +---------+--------+-------+-------+-------+--------+--------+
        |    1    |   2    |   3   |   4   |   5   |   6    |    7   |
        +=========+========+=======+=======+=======+========+========+
        | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |
        +---------+--------+-------+-------+-------+--------+--------+
        | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |
        +---------+--------+-------+-------+-------+--------+--------+
        """
        model = self.model
        aelist_id = max(model.aelists) + 1 if model.aelists else 1
        aesurf_id = aelist_id
        aesurf = {}
        aelists = {}
        for unused_aesurf_name, aesurfi in sorted(model.aesurf.items()):
            aelist, aesurfi2 = aesurfi.convert_to_nastran(model, aesurf_id, aelist_id)
            aelists[aelist.sid] = aelist
            aesurf[aesurfi2.aesurf_id] = aesurfi2
            aesurf_id += 1
            aelist_id += 1
        return aesurf, aelists

    def _convert_trim(self):
        """Converts ZONA TRIM to TRIM"""
        trims = {}
        model = self.model
        for trim_id, trim in sorted(model.trims.items()):
            trim_new = trim.convert_to_nastran(model)
            trims[trim_id] = trim_new
        return trims

    def _convert_trimlnk(self):
        """Converts ZONA TRIMLNK to AELINK"""
        model = self.model
        assert isinstance(model.aelinks, dict), model.aelinks
        aelinks = {}
        for trim_id, trimlnk in sorted(self.trimlnk.items()):
            aelink = trimlnk.convert_to_nastran(model)
            aelinks[trim_id] = aelink
        return aelinks

    def __repr__(self):
        msg = '<ZONA>; nPANLSTs=%s nmkaeroz=%s' % (
            len(self.panlsts), len(self.mkaeroz),
        )
        return msg

    def get_bdf_cards(self, bulk_data_lines: list[str],
                      bulk_data_ilines: Optional[Any]=None,
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
                if line.strip():
                    print(line)
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
                    if card_name not in ZONA_CARDS:
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
                raise RuntimeError(f'There are tabs in:\n{joined_lines_n}')

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
        return cards_list, cards_dict, card_count


def get_dicts(zona: ZONA, method: str) -> list[dict]:
    assert method in ['xref', 'write'], f'method={method!r}'
    dicts = [
        # --------------general-------------
        # zona.aeroz,
        zona.atmos, zona.flutter_table,
        # -------------geometry-------------
        # zona.panlsts,  # special-list
        zona.pafoil,
        zona.attach,
        # -------------transient------------
        zona.mloads, zona.eloads,
        # --------------flutter-------------
        zona.nlfltr, zona.mkaeroz,
        # ---------------trim---------------
        # zona.trim,  # part of the main BDF
        zona.aeslink, zona.trimvar, zona.trimlnk,
        zona.trimfnc, zona.trimobj, zona.trimcon,
        # ---------------ase---------------
        zona.cjunct, zona.conct, zona.tfset, zona.cnctset,
        zona.ase, zona.asecont, zona.asesnsr, zona.asesns1,
        zona.asegain, zona.gainset,
        zona.mimoss, zona.sisotf,
        #
        zona.senset, zona.surfset,
        zona.mldtrim, zona.mldstat, zona.minstat, zona.mldprnt,
        zona.mldcomd, zona.mldtime,
        # zona.extinp, zona.extout,
        zona.loadmod, zona.rbred,
        # ---------------gust---------------
        zona.gloads, zona.dgust, zona.cgust,
        # ---------------other--------------
        zona.extfile,
        zona.dse,
        zona.dmil,
    ]
    if method == 'write':
        # these are xref'd by their parent
        dicts.extend([zona.extinp, zona.extout])
    return dicts

