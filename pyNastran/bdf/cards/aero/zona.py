# coding: utf-8
# pylint: disable=W0212,C0103
from __future__ import annotations
from typing import TextIO, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils import object_attributes, object_methods

from pyNastran.bdf.cards.aero.zona_cards.spline import (
    SPLINE1_ZONA, SPLINE2_ZONA, SPLINE3_ZONA,)
from pyNastran.bdf.cards.aero.zona_cards.geometry import (
    PANLST1, PANLST2, PANLST3, SEGMESH,
    CAERO7, BODY7, PAFOIL7, AESURFZ)
from pyNastran.bdf.cards.aero.zona_cards.plot import (
    PLTAERO, PLTMODE, )
from pyNastran.bdf.cards.aero.zona_cards.flutter import (
    FLUTTER_ZONA, MKAEROZ)
from pyNastran.bdf.cards.aero.zona_cards.trim import (
    TRIM_ZONA, TRIMVAR, TRIMLNK,)
from pyNastran.bdf.cards.aero.zona_cards.manuever import (
    MLOADS, EXTINP, EXTOUT, TRIMFNC, ACTU, LOADMOD, RBRED,)
from pyNastran.bdf.cards.aero.zona_cards.ase import (
    ASE, ASECONT, ASESNSR,
    CJUNCT, CONCT, TFSET, MIMOSS,
    SENSET, SURFSET, CNCTSET,
)
from pyNastran.bdf.cards.aero.zona_cards.bdf_tables import (
    TABLED1_ZONA, TABDMP1_ZONA)
from pyNastran.bdf.cards.aero.zona_cards.dmi import DMIL
from pyNastran.bdf.cards.aero.zona_cards.cards import (
    MLDPRNT, MLDSTAT, MINSTAT, MLDTRIM, MLDCOMD, MLDTIME,
    AEROZ, ACOORD, ATTACH,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

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
        assert key not in self.model.zona.mldstat, key
        assert key > 0, key
        self.model.zona.mldstat[key] = mldstat
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
        assert panlst.eid not in self.model.zona.panlsts, key
        assert key > 0, key
        self.model.zona.panlsts[key] = panlst
        self.model._type_to_id_map[panlst.type].append(key)

    def add_pafoil_object(self, pafoil: PAFOIL7) -> None:
        """adds an PAFOIL7/PAFOIL8 object"""
        key = pafoil.pid
        assert pafoil.pid > 0
        assert key not in self.model.zona.pafoil, '\npafoil7=\n%s old=\n%s' % (
            pafoil, self.model.zona.pafoil[key])
        self.model.zona.pafoil[key] = pafoil
        self.model._type_to_id_map[pafoil.type].append(key)

    def add_aesurfz_object(self, aesurf: AESURFZ) -> None:
        """adds an AESURFZ object"""
        key = aesurf.aesid
        model = self.model
        assert key not in model.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, model.aesurf[key])
        model.aesurf[key] = aesurf
        model._type_to_id_map[aesurf.type].append(key)

    def add_mloads_object(self, mloads: MLOADS) -> None:
        """adds an MLOADS object"""
        key = mloads.mloads_id
        model = self.model
        assert key not in model.aesurf, '\nmloads=\n%s old=\n%s' % (
            mloads, model.mloads[key])
        self.model.zona.mloads[key] = mloads
        model._type_to_id_map[mloads.type].append(key)

    def add_ase_object(self, ase: ASE) -> None:
        """adds an ASE object"""
        key = ase.ase_id
        model = self.model
        assert key not in model.aesurf, '\nase=\n%s old=\n%s' % (
            ase, model.ase[key])
        self.model.zona.asecont[key] = ase
        model._type_to_id_map[ase.type].append(key)

    def add_asecont_object(self, asecont: ASECONT) -> None:
        """adds an ASECONT object"""
        key = asecont.asecont_id
        model = self.model
        assert key not in model.aesurf, '\nasecont=\n%s old=\n%s' % (
            asecont, model.asecont[key])
        self.model.zona.asecont[key] = asecont
        model._type_to_id_map[asecont.type].append(key)

    def add_asesnsr_object(self, asesnsr: ASESNSR) -> None:
        """adds an ASESNSR object"""
        key = asesnsr.asesnsr_id
        model = self.model
        assert key not in model.aesurf, '\nasesnsr=\n%s old=\n%s' % (
            asesnsr, model.asesnsr[key])
        self.model.zona.asesnsr[key] = asesnsr
        model._type_to_id_map[asesnsr.type].append(key)

    def add_cjunct_object(self, cjunct: CJUNCT) -> None:
        """adds an CJUNCT object"""
        key = cjunct.cjunct_id
        model = self.model
        assert key not in model.aesurf, '\ncjunct=\n%s old=\n%s' % (
            cjunct, model.cjunct[key])
        self.model.zona.cjunct[key] = cjunct
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
        assert key not in model.aesurf, '\nactu=\n%s old=\n%s' % (
            actu, model.actu[key])
        self.model.zona.actu[key] = actu
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


class ZONA:
    def __init__(self, model):
        self.model = model
        self.caero_to_name_map = {}
        self._add_methods = AddMethods(model)

        #: store PANLST1,PANLST2,PANLST3
        self.panlsts: dict[int, PANLST1 | PANLST2 | PANLST3] = {}
        self.mkaeroz: dict[int, MKAEROZ] = {}
        self.trimvar: dict[int, TRIMVAR] = {}
        self.trimlnk: dict[int, TRIMLNK] = {}
        self.attach: dict[int, PLTAERO] = {}
        self.plotaero: dict[int, PLTAERO] = {}
        self.plotmode: dict[int, PLTMODE] = {}
        #: store PAFOIL7/PAFOIL8
        self.pafoil: dict[int, PAFOIL7] = {}
        self.mloads: dict[int, MLOADS] = {}
        # TODO: add me
        self.eloads: dict[int, MLOADS] = {}
        self.gloads: dict[int, MLOADS] = {}
        self.nlfltr: dict[int, MLOADS] = {}
        self.dse: dict[int, MLOADS] = {}

        self.extinp: dict[int, EXTINP] = {}
        self.extout: dict[int, EXTOUT] = {}
        self.trimfnc: dict[int, TRIMFNC] = {}
        self.loadmod: dict[int, LOADMOD] = {}
        self.rbred: dict[int, RBRED] = {}

        self.cmargin: dict[int, CMARGIN] = {}
        self.mldtrim: dict[int, MLDTRIM] = {}
        self.mldstat: dict[int, MLDSTAT] = {}
        self.minstat: dict[int, MINSTAT] = {}
        self.mldprnt: dict[int, MLDPRNT] = {}
        self.mldcomd: dict[int, MLDCOMD] = {}
        self.mldtime: dict[int, MLDTIME] = {}
        self.dmil: dict[tuple[str, int, int], DMIL] = {}
        self.actu: dict[int, ACTU] = {}
        self.ase: dict[int, ASE] = {}
        self.asecont: dict[int, ASECONT] = {}
        self.asesnsr: dict[int, ASESNSR] = {}
        self.cjunct: dict[int, CJUNCT] = {}
        self.conct: dict[int, CONCT] = {}
        self.tfset: dict[int, TFSET] = {}
        self.mimoss: dict[int, MIMOSS] = {}
        self.senset: dict[int, SENSET] = {}
        self.cnctset: dict[int, CNCTSET] = {}
        self.surfset: dict[int, SURFSET] = {}

        # TODO:
        self.sisotf: dict[int, SISOTF] = {}
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
            panlst.validate()
        for mkaeroz in self.mkaeroz.values():
            mkaeroz.validate()
        for trimvar in self.trimvar.values():
            trimvar.validate()
        for trimlnk in self.trimlnk.values():
            trimlnk.validate()
        for pafoil in self.pafoil.values():
            pafoil.validate()
        for attach in self.attach.values():
            attach.validate()

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
            'TRIM': (TRIM_ZONA, add_methods.add_trim_object),
            'TABLED1': (TABLED1_ZONA, add_methods.add_tabled_object),
            'TABDMP1': (TABDMP1_ZONA, add_methods.add_table_sdamping_object),
            'CAERO7': (CAERO7, add_methods.add_caero_object),
            'AEROZ': (AEROZ, add_methods.add_aeros_object),
            'AESURFZ': (AESURFZ, zona_add.add_aesurfz_object),
            'FLUTTER': (FLUTTER_ZONA, add_methods.add_flutter_object),
            'MLOADS': (MLOADS, zona_add.add_mloads_object),
            'ASE': (ASE, zona_add.add_ase_object),
            'ASECONT': (ASECONT, zona_add.add_asecont_object),
            'ASESNSR': (ASESNSR, zona_add.add_asesnsr_object),
            'SENSET': (SENSET, zona_add.add_senset_object),
            'CNCTSET': (CNCTSET, zona_add.add_cnctset_object),
            'SURFSET': (SURFSET, zona_add.add_surfset_object),
            'CJUNCT': (CJUNCT, zona_add.add_cjunct_object),
            'CONCT': (CONCT, zona_add.add_conct_object),
            'ACTU': (ACTU, zona_add.add_actu_object),
            'LOADMOD': (LOADMOD, zona_add.add_loadmod_object),
            'RBRED': (RBRED, zona_add.add_rbred_object),
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
            'MIMOSS': (MIMOSS, zona_add.add_mimoss_object),
            'DMIL': (DMIL, zona_add.add_dmil_object),
            #'MLDPRNT': (MLDPRNT, self.add_mldprnt_object),
            # PLTCP
            # PLTFLUT
            # PLTVG
            # PLTMIST
        }
        card_parser.update(card_parser2)
        cards = [
            'CAERO7', 'AEROZ', 'AESURFZ', 'ATTACH',
            'PANLST1', 'PANLST2', 'PANLST3', 'PAFOIL7',
            'SEGMESH', 'BODY7', 'ACOORD', 'MKAEROZ',
            'TRIMVAR', 'TRIMLNK', 'FLUTTER',
            'PLTMODE', 'PLTAERO',
            # mloads
            'MLOADS',
            'ASE', 'ASECONT', 'ASESNSR',
            'CJUNCT', 'CONCT', 'TFSET',
            'MLDSTAT', 'MLDTRIM', 'MLDPRNT', 'MLDCOMD',
            'EXTINP', 'EXTOUT', 'TRIMFNC', 'LOADMOD',
            'DMIL',
            'MLDTIME', 'MLDCOMD', 'ACTU',
            'PLTTIME', 'MINSTAT', 'AEROLAG', 'APCONST',
            # 'PLTFLUT', 'PLTVG', 'PLTCP', 'PLTMIST',
            'FIXMDEN', 'FIXHATM', 'FIXMACH',
            'RBRED',
            # 'SPLINE0', 'PBODY7',
            'MIMOSS',
            'SENSET', 'CNCTSET', 'SURFSET',
        ]

        self.model.cards_to_read.update(set(cards))

    def cross_reference(self):
        if self.model.nastran_format != 'zona':
            return

        dicts = [
            self.mkaeroz, self.trimvar, self.trimlnk,
            self.attach, self.pafoil,
            self.mloads, self.eloads, self.gloads, self.dse, self.nlfltr,
            self.cjunct, self.conct, self.tfset,
            self.ase, self.asecont, self.asesnsr,
            self.senset, self.cnctset, self.surfset,
            self.mimoss,
            self.mldtrim, self.mldstat, self.mldprnt,
            self.mldcomd, self.mldtime, self.pafoil,
            # self.extinp, self.extout,
            self.trimfnc,
            self.loadmod, self.rbred,
        ]

        for items in dicts:
            for item in items.values():
                # self.model.log.info(f'xref {item.type}')
                item.cross_reference(self.model)
        # for mkaeroz in self.mkaeroz.values():
        #     mkaeroz.cross_reference(self.model)
        # for trimvar in self.trimvar.values():
        #     trimvar.cross_reference(self.model)
        # for trimlnk in self.trimlnk.values():
        #     trimlnk.cross_reference(self.model)
        # for unused_id, pafoil in self.pafoil.items():
        #     pafoil.cross_reference(self.model)
        # for unused_id, attach in self.attach.items():
        #     attach.cross_reference(self.model)
        #for aeroz in self.aeroz.values():
            #aeroz.cross_reference(self.model)

        for caero in self.model.caeros.values():
            #print('%s uses CAERO eid=%s' % (caero.label, caero.eid))
            self.caero_to_name_map[caero.label] = caero.eid

    def safe_cross_reference(self):
        self.cross_reference()

    def write_bdf(self, bdf_file: TextIO, size: int=8,
                  is_double: bool=False):
        #if self.model.nastran_format != 'zona':
            #return
        for unused_id, panlst in self.panlsts.items():
            bdf_file.write(panlst.write_card(size=size, is_double=is_double))
        for unused_id, mkaeroz in self.mkaeroz.items():
            bdf_file.write(mkaeroz.write_card(size=size, is_double=is_double))
        for unused_id, trimvar in self.trimvar.items():
            bdf_file.write(trimvar.write_card(size=size, is_double=is_double))
        for unused_id, trimlnk in self.trimlnk.items():
            bdf_file.write(trimlnk.write_card(size=size, is_double=is_double))
        for unused_id, pafoil in self.pafoil.items():
            bdf_file.write(pafoil.write_card(size=size, is_double=is_double))
        for unused_id, attach in self.attach.items():
            bdf_file.write(attach.write_card(size=size, is_double=is_double))

        dicts = [
            self.mloads, self.eloads, self.gloads, self.dse, self.nlfltr,
            self.cjunct, self.conct, self.tfset,
            self.ase, self.asecont, self.asesnsr,
            self.senset, self.cnctset, self.surfset,
            self.mimoss,
            self.mldtrim, self.mldstat, self.minstat, self.mldprnt,
            self.mldcomd, self.mldtime,
            self.extinp, self.extout, self.trimfnc,
            self.loadmod, self.rbred,
        ]
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
