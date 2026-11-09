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

def _decode_media_fmt(mediai: int, fmti: int) -> list[str]:
    """Decode the MEDIA and FMT bit fields from CASECC into option strings.

    MEDIA: bit 0=PRINT, bit 1=PLOT, bit 2=PUNCH
    FMT:   bit 0: 1=SORT1 (else SORT2), bit 1: 1=PHASE (else REAL)
           0 means unset/default — omit sort/format options
    """
    options: list[str] = []
    if mediai & 1:
        options.append('PRINT')
    if mediai & 2:
        options.append('PLOT')
    if mediai & 4:
        options.append('PUNCH')
    if fmti == 0:
        return options
    if not (fmti & 1):
        options.append('SORT2')
    if fmti & 2:
        options.append('PHASE')
    return options


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

        # --- H5 CASECC field name → case control command ---
        #
        # Integer scalars: h5_name maps to case_control_name (or same name)
        # SET/MEDIA/FMT triplets: output request commands
        #
        # MEDIA bit field: 1=PRINT, 2=PLOT, 3=PRINT+PLOT, 4=PUNCH
        # FMT bit field:   bit 0: 0=SORT2, 1=SORT1
        #                  bit 1: 0=REAL,  1=PHASE
        #                  so 0=default, 1=SORT1+REAL, 2=SORT2+REAL, 3=SORT1+PHASE
        # SET values:      -1=ALL, 0=not requested, >0=SET ID

        basic_ints = {
            'LOADSET',
            'DESOBJ', 'DESSUB',
            'TIC',  # TODO: thermal solutions use IC, not TIC
            'TRIM', 'DIVERG', 'GUST',
            'RANDOM',
            'CYCLIC', 'CLOAD',
            'DYNRED',
            'BSQUEAL',
            'ADAPT',
            'SUPORT1',
            'NLSTEP', 'NMODES',
            'APRESS',
        }

        integer_map = {
            'SPCSET': 'SPC',
            'MPCSET': 'MPC',
            'ESLSET': 'LOAD',
            'REESET': 'METHOD',
            'CEESET': 'CMETHOD',
            'DYMLDSET': 'DLOAD',
            'FEQRESET': 'FREQUENCY',
            'TSTEPTRN': 'TSTEP',  # TODO: nonlinear solutions use TSTEPNL, not TSTEP
            'ELDSET': 'DEFORM',
            'THLDSET': 'TEMP(LOAD)',
            'THMATSET': 'TEMP(MAT)',
            'TFSET': 'TF',
            'TFLSET': 'TFL',
            'MONSET': 'MONITOR',
            'FLUTTER': 'FMETHOD',
            'DAMPTBL': 'SDAMPING',
            'NONLINLD': 'NONLINEAR',
            'STATSUBP': 'STATSUB',
            'MODTRKID': 'MODTRAK',
            'NSMID': 'NSM',
            'BCID': 'BC',
            'STATSUBB': 'STATSUB(BUCKLING)',
        }
        _matrices = [
            'M2GG', 'B2GG', 'K2PP', 'B2PP', 'K2GG', 'P2G',
            'K42GG', 'A2GG', 'M2PP', 'BCONTACT',
        ]
        strings = ['SUBTITLE', 'LABEL', 'ANALYSIS'] + _matrices

        # output request triplets: case_control_name → (SET, MEDIA, FMT)
        kvo3 = {
            'DISPLACEMENT': ('DPLPTSET', 'DPLMEDIA', 'DPLFMT'),
            'VELOCITY':     ('VELPTSET', 'VELMEDIA', 'VELFMT'),
            'ACCELERATION': ('ACCPTSET', 'ACCMEDIA', 'ACCFMT'),
            'STRESS':       ('STSPTSET', 'STSMEDIA', 'STSFMT'),
            'STRAIN':       ('STNSET',   'STNMEDIA', 'STNFMT'),
            'FORCE':        ('FOCPTSET', 'FOCMEDIA', 'FOCFMT'),
            'SPCFORCES':    ('FCEPTSET', 'FCEMEDIA', 'FCEFMT'),
            'MPCFORCES':    ('MPCFSET',  'MPCMEDIA', 'MPCFFMT'),
            'OLOAD':        ('LDSPTSET', 'LDSMEDIA', 'LDSFMT'),
            'GPFORCE':      ('GPFSET',   'GPFMEDIA', 'GPFFMT'),
            'ESE':          ('ESESET',   'ESEMEDIA', 'ESEFMT'),
            'GPKE':         ('GPKESET',  'GPKEMEDI', 'GPKEFMT'),
            'EKE':          ('EKEPTSET', 'EKEMEDIA', 'EKEFMT'),
            'EDE':          ('EDEPTSET', 'EDEMEDIA', 'EDEFMT'),
            'NLSTRESS':     ('NLSSET',   'NLSMEDIA', 'NLSFMT'),
            'NONLINEAR':    ('NONPTSET', 'NONMEDIA', 'NONFMT'),
            'SDISPLACEMENT': ('SSDSET',  'SSDMEDIA', 'SSDFMT'),
            'SVELOCITY':    ('SSVSET',   'SSVMEDIA', 'SSVFMT'),
            'SACCELERATION': ('SSASET',  'SSAMEDIA', 'SSAFMT'),
            'SVECTOR':      ('SVSET',    'SVMEDIA',  'SVFMT'),
            'FLUX':         ('FLUPTSET', 'FLUMEDIA', 'FLUFMT'),
            'GPSTRESS':     ('GPSPTSET', 'GPSMEDIA', 'GPSFMT'),
            'GPSTRAIN':     ('STFSET',   'STFMEDIA', 'STFFMT'),
            'BOUTPUT':      ('CNTSET',   'CNTMEDIA', 'CNTFMT'),
            'ELSDCON':      ('ESDPTSET', 'ESDMEDIA', 'ESDFMT'),
            'GPSDCON':      ('GSDPTSET', 'GSDMEDIA', 'GSDFMT'),
            'AEROF':        ('ARFPTSET', 'ARFMEDIA', 'ARFFMT'),
            'GPEPSILON':    ('GPEPTSET', 'GPEMEDIA', 'GPEFMT'),
            'DATAREC':      ('DATSET',   'DATMEDIA', 'DATFMT'),
        }
        kvo3_skip = set()
        for values in kvo3.values():
            kvo3_skip.update(values)

        sids = subcase_data['SID']
        for key in basic_ints:
            # if key not in subcase_data.dtype.names:
            #     continue
            data = subcase_data[key]
            for (sid, valuei) in zip(sids, data):
                if valuei == 0:
                    continue
                subcases[sid].add_integer_type(key, valuei)
            log.debug(f'basic int {key}: {data}')

        for key, mapped_name in integer_map.items():
            # if key not in subcase_data.dtype.names:
            #     continue
            data = subcase_data[key]
            for (sid, valuei) in zip(sids, data):
                if valuei == 0:
                    continue
                subcases[sid].add_integer_type(mapped_name, valuei)
            log.debug(f'int {key} -> {mapped_name}: {data}')

        for key in strings:
            # if key not in subcase_data.dtype.names:
            #     continue
            data = subcase_data[key]
            for (sid, valuei) in zip(sids, data):
                svalue = valuei.decode('latin1').strip()
                if svalue == '':
                    continue
                subcases[sid].add_string_type(key, svalue)
            log.debug(f'str {key}: {data}')

        for key, (kset, kmedia, kfmt) in kvo3.items():
            # if kset not in subcase_data.dtype.names:
            #     continue
            set_value = subcase_data[kset]
            iset = (set_value != 0)
            if not iset.any():
                continue
            media = subcase_data[kmedia][iset]
            fmt = subcase_data[kfmt][iset]
            for (sid, seti, mediai, fmti) in zip(
                    sids[iset], set_value[iset], media, fmt):
                value = 'ALL' if seti == -1 else seti
                options = _decode_media_fmt(mediai, fmti)
                subcases[sid].add(key, value, options, 'STRESS-type')

        # VONMISES flag applies to STRESS/STRAIN output options
        # if 'VONMISES' in subcase_data.dtype.names:
        vonmises = subcase_data['VONMISES']
        ivm = (vonmises != 0)
        for sid, vm in zip(sids[ivm], vonmises[ivm]):
            if vm == 0:
                continue
            sc = subcases[sid]
            for result_key in ('STRESS', 'STRAIN'):
                if result_key not in sc.params:
                    continue
                value, options, param_type = sc.params[result_key]
                if vm & 1:
                    if 'VONMISES' not in options:
                        options.append('VONMISES')
                if vm & 2:
                    if 'BILIN' not in options:
                        options.append('BILIN')

        # Log any CASECC fields that aren't handled yet
        handled_fields = (
            {'SID', 'TITLE', 'VONMISES', 'DOMAIN_ID'}
            | basic_ints
            | set(integer_map.keys())
            | set(strings)
            | kvo3_skip
        )
        for field_name in subcase_data.dtype.names:
            if field_name in handled_fields:
                continue
            data = subcase_data[field_name]
            if data.dtype.kind in ('i', 'f'):
                if (data != 0).any():
                    log.warning(f'unhandled CASECC field {field_name}: {data}')
            elif data.dtype.kind in ('S', 'U'):
                for val in data:
                    sval = val.decode('latin1').strip() if isinstance(val, bytes) else val.strip()
                    if sval:
                        log.warning(f'unhandled CASECC field {field_name}: {data}')
                        break

        # TITLE goes to the global subcase (subcase 0)
        if 'TITLE' in subcase_data.dtype.names:
            title_data = subcase_data['TITLE']
            if len(title_data) > 0:
                svalue = title_data[0].decode('latin1').strip()
                if svalue:
                    subcases[0].add_string_type('TITLE', svalue)

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
