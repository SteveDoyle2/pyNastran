from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import cast_encoding_strip, cast_encoding_strip0, get_group_name
from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.case_control_deck import CaseControlDeck
if TYPE_CHECKING: # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def _load_h5_parameter_pvt(model: BDF, group: Group, encoding: str):
    log = model.log
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            group_name = get_group_name(group)
            raise RuntimeError(group_name)
            #msg += f'{indent}Group: {group_name}\n'
            #indent += '    '
            #print_group(group, indent=indent)
            #continue
        elif class_id == 'TABLE':
            group_name = get_group_name(group)
            #log.warning(f'PARAM table name={group_name}')
        else: # pragma: no cover
            raise RuntimeError(class_id)

        if group_name == 'PVT':
            datatype = h5_element.name
            data = h5_element.read()
            name_bytes = data['NAME']
            name = cast_encoding_strip(name_bytes, encoding)
            if datatype == 'INT':
                # NAME, VALUE
                value = data['VALUE']
                for key, valuei in zip(name, value):
                    model.add_param(key, [valuei])
            elif datatype == 'CHAR':
                # NAME, VALUE
                value_bytes = data['VALUE']
                value = cast_encoding_strip(value_bytes, encoding)
                for key, valuei in zip(name, value):
                    model.add_param(key, [valuei])
            elif datatype in {'DOUBLE'}:  # TODO: does FLOAT exist?
                # NAME, VALUE
                value = data['VALUE']
                for key, valuei in zip(name, value):
                    model.add_param(key, [valuei])
            else: # pragma: no cover
                raise NotImplementedError(datatype)
        else: # pragma: no cover
            raise NotImplementedError(group_name)
        x = 1

def _load_h5_parameter_group(model: BDF, group: Group, encoding: str):
    log = model.log
    assert log is not None, log
    group_name = get_group_name(group)
    skip_names = {'TSTEP', } # , 'PVT',
    if group_name in skip_names:
        log.warning(f'skipping {group_name} in _load_h5_parameter_group')
    elif group_name == 'CASECC':
        subcase_data = group['SUBCASE'].read()
        #print(f'subcase_data = {subcase_data}')
        #print(f'  dtype subcase_data = {subcase_data.dtype}')
        #for h5_elementi in group._f_iter_nodes():
            #print(h5_elementi)
            #print(h5_elementi.dtype.name)
        nsubcase = len(subcase_data)
        subcase_id_dict = {sid: {} for sid in subcase_data['SID']}
        cc = CaseControlDeck([], log=log)
        #cc.log = log
        for sid in subcase_id_dict:
            subcase = cc.create_new_subcase(sid)
            subcase.log = log
        model.case_control_deck = cc
        subcases = cc.subcases

        # 'SYMFLG'
        basic_ints = {
            'LOADSET', 
            'DESOBJ', 'DESSUB',
            'TIC',
            'TRIM', 'DIVERG', 'GUST',
            'RANDOM',
            'CYCLIC', 'CLOAD',
            'DYNRED',
            'BSQUEAL', }

        integer_map = {
            'SPCSET': 'SPC',
            'MPCSET': 'MPC',
            'TFSET': 'TF',
            'TFLSET': 'TFL',
            'MONSET': 'MONITOR',
            'FLUTTER': 'FMETHOD',  # FLUTTER
            'DAMPTBL': 'SDAMPING', # TABDMP1
            #'MODLIST': '', # MODTRAK
            #THMATSET, 
        }
        _matrices = [
            'M2GG', 'B2GG', 'K2PP', 'B2PP', 'K2GG', 'P2G',
            'K42GG', 'A2GG', 'M2PP', 'BCONTACT',
        ]
        strings = ['SUBTITLE', 'LABEL'] + _matrices
        bools = [
            'AESYMXY', 'AESYMXZ',
            'ExtDROut', 'ExtDROmd', 'ExtDROun',
            'ExtDRIn', 'ExtDRImd', 'ExtDRIun',
        ]

        kvo3 = {  # key-value-options-3
            #'SPCFORCE':
            'MPCFORCE': ('MPCFSET', 'MPCMEDIA', 'MPCFFMT'),
            'DISP':     ('DPLPTSET', 'DPLMEDIA', 'DPLFMT'),
            'VELO':     ('VELPTSET', 'VELMEDIA', 'VELFMT'),
            'ACCEL':    ('ACCPTSET', 'ACCMEDIA', 'ACCFMT'),
            'STRAIN':   ('STNSET', 'STNMEDIA', 'STNFMT'),
            #'GPKE':     ('GPKESET', 'GPKEMEDIA', 'GPKEFMT'),  # broken
            'ESE':      ('ESESET', 'ESEMEDIA', 'ESEFMT'),
            'GPFORCE':  ('GPFSET', 'GPFMEDIA', 'GPFFMT'),
            #'STRESS:   VONMISES,
            #'EKE':      ('EKEPTSET', 'EKEMEDIA', 'EKEFMT', 'EKETHRSH', ???
            #'FORCE':
        }
        kvo3_skip = []
        for values in kvo3.values():
            kvo3_skip.extend(list(values))

        sids = subcase_data['SID']
        all_names = list(subcase_data.dtype.names)
        for key in basic_ints:
            data = subcase_data[key]
            for (sid, valuei) in zip(sids, data):
                assert isinstance(valuei, integer_types), (key, valuei)
                #assert valuei > 0, valuei
                if valuei == 0:
                    continue
                subcases[sid].add_integer_type(key, valuei)
            log.debug(f'basic int {key}: {data}')

        for key, mapped_name in integer_map.items():
            data = subcase_data[key]
            for (sid, valuei) in zip(sids, data):
                assert isinstance(valuei, integer_types), (key, mapped_name, valuei)
                if valuei == 0:
                    continue
                subcases[sid].add_integer_type(mapped_name, valuei)
            log.debug(f'int {key}: {data}')

        for key in strings:  # comes in as bytes
            data = subcase_data[key]
            for (sid, valuei) in zip(sids, data):
                svalue = valuei.decode('latin1').strip()
                if valuei == '':
                    continue
                subcases[sid].add_string_type(key, svalue)
                #log.debug(f'{key}: {data}; value={value!r} type={type(value)}')
            log.debug(f'str {key}: {data}')

        for key, (kset, kmedia, kfmt) in kvo3.items():
            set_value = subcase_data[kset]
            iset = (set_value != 0)
            #if set_value.max() == 0 and set_value.min() == 0:
            #    continue
            media = subcase_data[kmedia][iset]
            fmt = subcase_data[kfmt][iset]
            for (sid, seti, mediai, fmti) in zip(sids[iset], set_value[iset], media, fmt):
                assert isinstance(seti, integer_types), seti
                assert isinstance(mediai, integer_types), mediai
                assert isinstance(fmti, integer_types), fmti
                seti = 'ALL' if seti == -1 else seti
                value = seti
                smedia = 'PLOT'
                sort = 'SORT1' if fmti == 1 else 'SORT2'
                options = [smedia, sort]  # media
                    
                #assert mediai == 0, (sid, key, seti, mediai, fmti)
                assert fmti in [1, 3], (sid, key, seti, mediai, fmti)
                subcases[sid].add(key, value, options, 'STRESS-type')

        #----------------------------------
        for key in subcase_data.dtype.names:
            if key in basic_ints or key in integer_map:
                continue
            if key in strings or key in kvo3_skip:
                continue
            data = subcase_data[key]
            data2 = data.flatten()
            #log.info(f'{key}: data={data} data2={data2}')

            value = data2[0]
            if isinstance(value, integer_types) and data.max() == 0 and data.min() == 0:
                continue
            if isinstance(value, float_types) and data.max() == 0. and data.min() == 0.:
                continue
            elif isinstance(value, str):
                data = np.array([value.strip() for value in data])
                value = data[0]
                if np.all(data == ''):
                    continue
            elif isinstance(value, bytes):
                data = np.array([value.strip() for value in data])
                value = data[0]
                if np.all(data == b''):
                    continue

            if key == 'TITLE':  # bytes?
                assert len(data) == 1, data
                valuei = data[0]
                svalue = valuei.decode('latin1').strip()
                if svalue == '':
                    continue
                #raise RuntimeError((key, svalue, type(sval))
                subcases[0].add_string_type(key, svalue)
                continue
            #log.debug(f'{key}: {data}; value={value!r} type={type(value)}')

        # has a bug...
        #log.debug(str(cc))
        #
        #
        # bytes:
        # int: SYMFLG,
        #
        # LDSPTSET, LDSMEDIA, LDSFMT,
        # NONPTSET, NONMEDIA, NONFMT,
        # STSPTSET, STSMEDIA, STSFMT,
        # FCEPTSET, FCEMEDIA, FCEFMT,
        # FOCPTSET, FOCMEDIA, FOCFMT,
        # SSDSET,   SSDMEDIA, SSDFMT,
        # SSVSET,   SSVMEDIA, SSVFMT,
        # SSASET,   SSAMEDIA, SSAFMT,
        # ARFPTSET, ARFMEDIA, ARFFMT,
        # SVSET,    SVMEDIA,  SVFMT,
        # FLUPTSET, FLUMEDIA, FLUFMT,
        # GPSPTSET, GPSMEDIA, GPSFMT,
        # STFSET,   STFMEDIA, STFFMT,
        # ESDPTSET, ESDMEDIA, ESDFMT,
        # GSDPTSET, GSDMEDIA, GSDFMT,
        # CNTSET,   CNTMEDIA, CNTFMT,
        # DATSET,   DATMEDIA, DATFMT,
        # VUGSET,   VUGMEDIA, VUGFMT,
        # NLSSET,   NLSMEDIA, NLSFMT,
        # GPEPTSET, GPEMEDIA, GPEFMT,
        # EDEPTSET, EDEMEDIA, EDEFMT,
        # GVSET,    GVMEDIA,  GVFMT, 
        # MDLSSET,  MDLSMEDIA, MDLSFMT, MDLSESRT, MDLSTHRE, MDLSTFVL,
        # MDLKSET,  MDLKMEDIA, MDLKFMT, MDLKESRT, MDLKTHRE, MDLKTFVL
        # VITSOSET, VITSRSET, VITSSSET, VITSTHRS, VITSBITS
        # NVELOSET, NVELFDEF, NVELFDTA, NVELTHRS, NVELBITS, 
        # ICFSET,   ICFMED,   ICFFMT,
        # ACPOWSET, ACPOWMED, ACPOWFMT,
        # MXMNGSET, MXMNGMDA, MXMNGFMT,
        # MXMNESET, MXMNEMDA, MXMNEFMT,
        # PFMSSID,  PFMSMED,  PFMSFMT,
        # PFMFSID,  PFMFMED,  PFMFFMT,
        # PFPSID,   PFPMED,   PFPFMT,
        # PFGSID,   PFGMED,   PFGFMT,
        # ACFPMSET, ACFPMMED, ACFPMFMT,
        # INTENSET, INTENMED, INTENFMT,
        # ERPSID,   ERPMED,   ERPFMT,
        # PFPSSID,  PFPSMED,  PFPSFMT,
        # RCRSET,             RCRFMT,
        # HISTSET,  HISTTYPE, HISTFMT,
        # FATIGUE,  FTGMED,   FTGFMT,
        # ICFGENST, ICFGENNM, ICFUSEST, ICFUSENM,
        # PACOSET,  PACSSET, PACBITS,
        # WTSOSET,  WTSRSET, WTSSSET, WTSTHRS, WTSBITS,
        # APRESS,
        # ELSOSET, ELSRSET, ELSSSET, ELSTHRS, ELSBITS,
        # PFPSSSET, PFPSFLTR, PFPSOPTS, PFPSPSET,
        # NLOPCTRL, NLOPDBG, NLOPPOST, NLOPMPCH,
        # TFOSET,  TFUNIT,  TFOPTS,
        # ERPSSET, ERPFLTR, ERPOPTS, ERPCSV,
        # ASMOUTFL, ASMOUTNM,
        # SEMG, SEKR, SELG, SELR, SEDR, SEMR, SEDV, SERE, SERS, SEEX,
        # LSEM, SYM_LEN, SYM_POS, 
        # {'names': [ESLSET, REESET, ELDSET, THLDSET, DYMLDSET, FEQRESET, TSTEPTRN, STPLTFLG, AXSYMSET, NOHARMON, TSTRU, SETKP, FLAGKP, SETMP, FLAGMP, SETBP, FLAGBP, OUTRESPV, FLDBNDY, CEESET, NONLINLD, PARTIT, NONPARAM, LCC, SEID, LCN, SEFINAL, SETKG, FLAGKG, SETMG, FLAGMG, SETBG, FLAGBG, HOUT, NOUT, SETPG, FLAGPG, SECMDFLG, SET2ID, DSAPRT, DSASTORE, DSAOUTPT, REESETF, OUTRCV, STATSUBP, DFT1, DFT2, ADAPT, SUBSPAN, DESGLB, ANALYSIS, GPQSTRS, GPQFORC, GPQSTRN, SUPORT1, STATSUBB, BCID, AUXMODEL, ADACT, REUESET, DAMPTBLF, ITERMETH, MODTRKID, DSAFORM, DSAEXPO, DSABEGIN, DSAINTVL, DSAFINAL, DSASETID, SORTFLG, RANDBIT, AECONFIG, OCIDREQ, TEMPMAT, AECSSSET, EDETHRSH, DFT3, DFR1, DFR2, DFR3, SETK4G, FLAGK4G, NK42GG, NA2GG, EFFMASET, EFFMAGID, EFFMATHR, EQUILMED, EQUILGRD, AEUXREF, GCHK, GCHKOUT, GCHKSET, GCHKGID, GCHKTHR, GCHKRTHR, GCHKDREC, ASPCMED, ASPCEPS, ASPCPRT, ASPCPCH, NK2PP, NM2PP, NB2PP, NK2GG, NM2GG, NB2GG, NP2G, GEODSET, GEODMXMN, GEODOCID, GEODNUMB, GEOLSET, GEOLMXMN, GEOLOCID, GEOLNUMB, GEOSSET, GEOSMXMN, GEOSOCID, GEOSNUMB, GEOMSET, GEOMMXMN, GEOMOCID, GEOMNUMB, GEOASET, GEOAMXMN, GEOAOCID, GEOANUMB, GEOVSET, GEOVMXMN, GEOVOCID, GEOVNUMB, NTFL, SEDAMP, WCHK, WCHKOUT, WCHKSET, WCHKGID, WCHKCGI, WCHKWM, EXSEOUT, EXSEMED, EXSEUNIT, EXSERES1, EXSERES2, FK2PP, FM2PP, FB2PP, FK2GG, FM2GG, FB2GG, TICTYPE, FK42GG, FA2GG, SUBSTEP, STEPID, NSMID, ROUTDISP, ROUTVELO, ROUTACCE, ROUTLOAD, ROUTSPCF, ROUTSTRS, ROUTFORC, ROUTSTRN, ROUTMSCF, ACPOWCSV, NLOUT, SEEFMNO, SEEFMHV, SEEFMDLF, SEEFMBND, SEEFMSDP, CONNECTOR, ESETHRSH, POSTUNIT, POSTOPT1, POSTOPT2, TICDIFF, DSAESEID, MCFRSET, MCFRSOLN, MCFRFILT, MCFROPT, ELSUMID, ELSUMOPT, ELSUMDUM, RGYRO, CMSESET, CMSEMDIA, CMSEOPTS, CMSETHRE, CMSETOPN, GPRSORT, MASSSET, AESOLN, POSTO2NM, RANDVAR, RSVCRQTS, RSVCOPTS, RSVCSTBS, RSVCRQTC, RSVCOPTC, RSVCSTBC, DESVAR, BCONTACTI, BCONTACTC, MODSELS1, MODSELS2, MODSELS3, MODSELS4, MODSELS5, MODSELS6, MODSELS7, MODSELS8, MODSELS9, MODSELF1, MODSELF2, MODSELF3, MODSELF4, MODSELF5, MODSELF6, MODSELF7, MODSELF8, MODSELF9, FTNURN, SUFNAM1, SUFNAM2, ENVELOP1, ENVELOP2, GPFLXSET, GPFLXMED, CAMPBELL, SPLINOUT, MONITOR, FBODYLD, STOCHAST, EXPTLDID, EXPTLDNM, EXPTLDSI, AERCONFIG, NLHARM, PFMSSSET, PFMSFLTR, PFMSOPTS, PFMSSTMP, PFMFSSET, PFMFPSET, PFMFFLTR, PFMFOPTS, PFMFFLMP, PFMFSTMP, PFPSSET, PFPFLTR, PFPOPTS, PFPPSET, PFGSSET, PFGGSET, NLICCASE, NLICSTEP, NLICLFAC, FRFFLAG, FRFCMPID, FRFCONST, FRFUNTNO, FRFCMPNM, DESMOD, RSDAMPST, RSDAMPFL, VCCT, FP2G, FRQVAR, HADAPT, BCHANGE, BCMOVE, UNGLUE, HSUBCASE, HSTEP, HTIME, TESTTHRR, TESTTHRI, ICFUNTNO, NLOOPH, NLSTEP, NMODES, RCPARM, NLICSSTP, DEACTEL, ACTIVAT, NLICTOLR, ACTISET, IRLOAD, NLOPDELI, NLOPGRID, DASAVE, ERMPF, THLDVER, THMATVER, NPEAK, NEAR, LFREQ, HFREQ, RTYPE, PSCALE, NSAMP, SEED, OFFD, FSORT2, MFREQ, SOLNID, FemCheck, DOMAIN_ID]
        #datatype = group.name
        #model.log.warning(f'skipping {group_name}={datatype} in _load_h5_parameter_group')
        #data = group.read()
        #log.warning(type(group))
    elif group_name == 'PVT':
        _load_h5_parameter_pvt(model, group, encoding)
    elif group_name == 'AECOMP':
        #('NAME', 'LISTTYPE', 'LISTIDS_POS', 'LISTIDS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        name = identity['NAME']
        list_type = cast_encoding_strip(identity['LISTTYPE'], encoding)
        nlists = identity['LISTIDS_LEN']
        lists = group['LISTIDS'].read()
        model.aecomp._save(name, list_type, nlists, lists)
    elif group_name == 'AECOMPL':
        #('NAME', 'LABELS_POS', 'LABELS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        name = identity['NAME']
        nlabels = identity['LABELS_LEN']
        labels_bytes = group['LABELS'].read()
        labels = cast_encoding_strip0(labels_bytes, encoding)
        model.aecompl._save(name, nlabels, labels)
    elif group_name == 'AEFACT':
        #('SID', 'DIS_POS', 'DIS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        aefact_id = identity['SID']
        nfractions = identity['DIS_LEN']
        fractions = group['DIS'].read()
        model.aefact._save(aefact_id, nfractions, fractions)
    else:
        raise NotImplementedError(group_name)

def load_h5_parameter(model: BDF, group: Group, encoding: str='latin1'):
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            _load_h5_parameter_group(model, h5_element, encoding)
            continue
        elif class_id == 'TABLE':
            pass
        else: # pragma: no cover
            raise RuntimeError(class_id)

        name = h5_element.name
        data = h5_element.read()
        #if name in {'PVT'}:
            #model.log.warning(f'skipping {element_name} in load_h5_parameter')
        if name == 'TSTEPNL':
            model.log.warning(f'skipping {name} in load_h5_parameter')
        elif name == 'AEROS':
            #('ACSID', 'RCSID', 'REFC', 'REFB', 'REFS', 'SYMXZ', 'SYMXY', 'DOMAIN_ID')
            acsid = data['ACSID'][0]
            rcsid = data['RCSID'][0]
            bref = data['REFB'][0]
            cref = data['REFC'][0]
            sref = data['REFS'][0]
            sym_xz = data['SYMXZ'][0]
            sym_xy = data['SYMXY'][0]
            aeros = model.add_aeros(cref, bref, sref, acsid=acsid,
                                    rcsid=rcsid, sym_xz=sym_xz, sym_xy=sym_xy, comment='')
            str(aeros)
        elif name == 'AESTAT':
            #('ID', 'LABEL', 'DOMAIN_ID')
            aestat_ids = data['ID']
            labels = data['LABEL']
            for aestat_id, label in zip(aestat_ids, labels):
                label_str = label.strip().decode('latin1')
                model.add_aestat(aestat_id, label_str, comment='')
        elif name == 'MDLPRM':
            #('NAME', 'VALUE', 'DOMAIN_ID')
            names = data['NAME']
            values = data['VALUE']
            mdlprm_dict = {}
            for name, value in zip(names, values):
                name_str = name.strip().decode('latin1')
                mdlprm_dict[name_str] = value
            assert len(mdlprm_dict) == 1, mdlprm_dict
            assert 'HDF5' in mdlprm_dict, mdlprm_dict
            model.add_mdlprm(mdlprm_dict, comment='')
        elif name == 'MONPNT1':
            name = data['NAME']
            label = data['LABEL']
            axes = data['AXES']
            comp = data['COMP']
            cp = data['CP']
            xyz = np.stack(
                [data['X'], data['Y'], data['Z']],
                axis=1)
            assert xyz.shape[1] == 3, xyz.shape
            cd = data['CD']
            model.monpnt1._save(
                name, label, comp, axes, cp, xyz, cd)
        else: # pragma: no cover
            raise NotImplementedError(name)
