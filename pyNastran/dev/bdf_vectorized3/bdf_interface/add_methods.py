from __future__ import annotations
from typing import Union, TYPE_CHECKING

#import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import BDFAttributes
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, PARAM
    from pyNastran.bdf.bdf import (
        MDLPRM,
        EIGRL, EIGR, EIGB, EIGP, EIGC,
        DOPTPRM,
        TSTEP, TSTEP1, TSTEPNL,
        AERO, AEROS,
        NLPARM, NLPCI,
        FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
        CAERO2, CAERO3, CAERO4, CAERO5,
        PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
        SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5,
        AESURF, AESURFS, TRIM, TRIM2, DIVERG, AESTAT,
        FLUTTER, MKAERO1, MKAERO2, GUST,
        DMIG, DMI, DMIAX, DMIJ, DMIJI, DMIK, DTI,
        #PACABS,
        SUPORT, SUPORT1,)

class AddMethods():
    def __init__(self, model: BDF):
        self.model = model

    # general values
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


    def _add_nxstrat_object(self, nxstrat: NXSTRAT) -> None:
        key = nxstrat.sid
        assert key not in self.model.nxstrats, 'nxstrats=%s nxstrat=%s' % (self.model.nxstrats, nxstrat)
        assert key > 0
        self.model.nxstrats[key] = nxstrat
        self.model._type_to_id_map[nxstrat.type].append(key)

    # SOL 101

    # SOL 103
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

    # SOL 107? - complex eigenvectors
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

    # SOL xxx - frequency
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

    # SOL xxx - transient
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

    # nonlinear
    def _add_nlparm_object(self, nlparm: NLPARM) -> None:
        """adds an NLPARM object"""
        key = nlparm.nlparm_id
        assert key not in self.model.nlparms
        assert key > 0, 'key=%s; nlparm=%s\n' % (key, nlparm)
        self.model.nlparms[key] = nlparm
        self.model._type_to_id_map[nlparm.type].append(key)

    def _add_nlpci_object(self, nlpci: NLPCI) -> None:
        """adds an NLPCI object"""
        key = nlpci.nlpci_id
        assert key not in self.model.nlpcis
        assert key > 0
        self.model.nlpcis[key] = nlpci
        self.model._type_to_id_map[nlpci.type].append(key)

    #---------------------------------------------------------------------------
    # SOL 144/145/146 - general aero
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

    # SOL 144 - static aero
    def _add_diverg_object(self, diverg: DIVERG, allow_overwrites: bool=False) -> None:
        """adds an DIVERG object"""
        key = diverg.sid
        if not allow_overwrites:
            assert key not in self.model.divergs, 'DIVERG=%s  old=\n%snew=\n%s' % (key, self.model.divergs[key], diverg)
        assert key > 0, 'key=%r diverg=\n%s' % (key, diverg)
        self.model.divergs[key] = diverg
        self.model._type_to_id_map[diverg.type].append(key)

    # SOL 145 - flutter
    def _add_flutter_object(self, flutter: FLUTTER) -> None:
        """adds an FLUTTER object"""
        key = flutter.sid
        assert key not in self.model.flutters, 'FLUTTER=%s old=\n%snew=\n%s' % (key, self.model.flutters[key], flutter)
        assert key > 0
        self.model.flutters[key] = flutter
        self.model._type_to_id_map[flutter.type].append(key)

    def _add_mkaero_object(self, mkaero: Union[MKAERO1, MKAERO2]) -> None:
        """adds an MKAERO1/MKAERO2 object"""
        self.model.mkaeros.append(mkaero)

        # SOL 200
    def _add_doptprm_object(self, doptprm: DOPTPRM) -> None:
        """adds a DOPTPRM"""
        self.model.doptprm = doptprm

    #  matricies
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
