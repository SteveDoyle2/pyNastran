from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.dev.bdf_vectorized3.bdf import BDF as BDFv
from pyNastran.utils import object_attributes, object_methods
from pyNastran.bdf.bdf_interface.utils import (
    to_fields, _get_card_name,
)

OBJ_NAME_DICTS = [
    'elements', 'rigid_elements', 'materials',
    'masses', 'properties_mass',
    'dmig', 'dmi', 'dmij', 'dmik', 'dmiji', 'dti', 'dmiax',

    'caeros', 'paeros', 'trims', 'splines', 'csschds', 'divergs',
    'flutters', 'flfacts', 'gusts', 
    'aefacts', 'aelists', 'aeparams',
    'aestats', 'aesurf', 'aesurfs', 
    'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9', 'MATT11',
    'MATS1', 'MATS3', 'MATS8',
    'nlparms', 'nlpcis', 'normals', 'plotels', 
    'aecomps',
    'desvars', 'dvprels', 'dvmrels', 'dvtrels',
    'pbusht', 'pdampt', 'pelast',
    'sets',
]
OBJ_NAME_DICT_LISTS = [
    'loads', 'load_combinations',
    'dloads', 'dload_entries',
    'nsmadds', 'nsms',
    'spcs', 'mpcs',
    'spcadds', 'mpcadds',
    'dresps', 'dconstrs',
    'bctadds', 'bctparas',
    'bgadds', 'bgsets', 
    'dconadds',
    'aelinks',
    'se_sets', 'se_usets', 'spcoffs',
]
OBJ_NAME_LISTS = [
    'asets', 'bsets', 'csets', 'omits', 'qsets',
    'se_bsets', 'se_csets', 'se_qsets',
    'se_suport',
    'suport',
    'units',
    'mkaeros', 
]

COMMON_SINGLETON_CARDS = ['aero', 'aeros', 'dtable',]
COMMON_CARDS = [
    #'mkaeros',
    #-----------------------------------
    'params', #'rotors',
    'tables', 'tables_d', #'tables_h'
    'tables_m',
    'tables_sdamping',
    'methods', 'cMethods',
]
COMMON_KEYS = [
    'superelement_models', 'log',
    'system_command_lines',
    'executive_control_lines',
    'case_control_deck',
    'punch', 'dumplines', 
    'reject_cards', 'reject_lines',
    'active_filenames',
    # bool
    'echo', 'debug', 'read_includes',
    # str
    'include_dir',
]

REMOVE_PARAMS = [
    'card_count', 
    'active_filename', 
    'allow_duplicate_element_rbe_mass',
    'writer', 'xref_obj', 'zaero', #'zona',
    'is_msc', 'is_mystran', 'is_nx', 'is_optistruct',
    'is_strict_card_parser', 'is_superelements', 'is_zaero',
    'allow_overwrites_set', 'allow_tabs',
    'rsolmap_to_str', 'save_file_structure',
    'use_new_deck_parser', 'values_to_skip',
    'force_echo_off', 'dict_of_vars',
    'cards_to_read', 'replace_includes',
]
#ADD_CARDS = [
#    'suport1',
#]
SKIP_FINAL = [
    # TODO: remove these
    'aedw', 'aeforce', 'aepress',
    'feedge', 'feface', 'gmcurv', 'gmsurf', 'pset', 'pval',
    'sol', 'sol_iline', 

    # parameters
    'type_slot_str',
    'is_long_ids',
    'cards_to_read',
    'coord_ids', 'ncoords',
    'caero_ids', 'ncaeros',
    'element_ids', 'nelements',
    'nastran_format',
    # not in BDFv
    'bdf_filename', 
    'case_control_lines',
    'special_cards',
    'include_filenames',
    # ?
    'app',
]
SKIP = [
    'MATDMG', 'acmodl', 'acplnw', 'amlreg', 'ao_element_flags',
    'baror', 'bcbodys', 'bconp', 'bcparas', 'bcrparas', 'bcs', 'bctparms', 'bctsets', 'beamor', 'bfric', 'blseg', 'bolt', 'boltfor', 'boltfrc', 'boltld', 'boltseq', 'bsurf', 'bsurfs', 'convection_properties', 'creep_materials', 'csuper', 'csupext', 'cyax', 'cyjoin', 'dareas', 'ddvals', 'delays', 'dequations', 'dlinks', 'dmncon', 'doptprm', 'dphases', 'dscreen', 'dvcrels', 'dvgrids', 'epoints', 'frequencies', 'grdset', 'group', 'hyperelastic_materials', 'matcid', 'mdlprm', 'micpnt', 'model_groups', 'modtrak', 'nxstrats', 'phbdys', 'points', 'radcavs', 'radmtx', 'radset', 'random_tables', 'reject_count', 'release', 'rotors', 'sebndry', 'sebulk', 'seconct', 'seelt', 'seexcld', 'selabel', 'seload', 'seloc', 'sempln', 'senqset', 'seqgp', 'setree', 'sol_method', 'spoints', 'tempds', 'thermal_materials', 'tics', 'topvar', 'transfer_functions', 'tstepnls', 'tsteps', 'type_slot_str',  'usets', 'uxvec', 'view3ds', 'views', ]

def bdf_to_bdf_vectorized(model: BDF,
                          clean_parent: bool=True) -> BDFv:
    modelv = BDFv(log=model.log)
    log = modelv.log
    for nid, node in model.nodes.items():
        assert node.type == 'GRID', node
        modelv.add_grid(
            nid, node.xyz, cp=node.cp,
            cd=node.cd, ps=node.ps, seid=node.seid)

    for cid, coord in model.coords.items():
        if cid == 0:
            continue
        fields = coord.raw_fields()
        card_name = fields[0]
        card = BDFCard(fields)
        modelv.add_card(card, card_name)

    for monpnt in model.monitor_points:
        if monpnt.type == 'MONPNT1':
            modelv.add_monpnt1(
                monpnt.name, monpnt.label,
                monpnt.axes, monpnt.comp,
                monpnt.xyz,
                cp=monpnt.cp,
                cd=monpnt.cd,
                comment=monpnt.comment)
        else:
            raise NotImplementedError(monpnt)

    assert isinstance(model.suport1, dict), model.suport1
    for sid, suport in model.suport1.items():
        assert not isinstance(suport, list), suport
        card = suport.raw_fields()
        #print(card)
        modelv.add_card(card, 'SUPORT1')

    for pid, prop in model.properties.items():
        prop_type = prop.type
        obj = getattr(modelv, prop_type.lower())
        if prop_type == 'PSHELL':
            modelv.add_pshell(
                pid, mid1=prop.mid1, t=prop.t,
                mid2=prop.mid2, twelveIt3=prop.twelveIt3,
                mid3=prop.mid3, mid4=prop.mid4, tst=prop.tst)
        else:
            fields = prop.raw_fields()
            card = BDFCard(fields)
            obj.add_card(card, ifile=0)

    #for mid, mat in model.materials.items():
    #    mat_type = mat.type
    #    fields = mat.raw_fields()
    #    card = BDFCard(fields)
    #    obj.add_card(card)

    for obj_name in SKIP:
        obj = getattr(model, obj_name)
        if isinstance(obj, list):
            log.warning(f'{obj_name} is a list')
    #return
    for obj_name in OBJ_NAME_DICTS:
        obj_dict = getattr(model, obj_name)
        assert isinstance(obj_dict, dict), (obj_name, obj_dict)
        for eid, elem in obj_dict.items():
            elem_type = elem.type
            objv = getattr(modelv, elem_type.lower())
            fields = elem.raw_fields()
            card = BDFCard(fields)
            objv.add_card(card, ifile=0)

    for obj_name in OBJ_NAME_DICT_LISTS:
        loads_dict = getattr(model, obj_name)
        assert isinstance(loads_dict, dict), (obj_name, loads_dict)
        for sid, loads in loads_dict.items():
            assert isinstance(loads, list), (obj_name, loads)
            for load in loads:
                load_type = load.type
                objv = getattr(modelv, load_type.lower())
                fields = load.raw_fields()
                card = BDFCard(fields)
                objv.add_card(card, ifile=0)

    for obj_name in OBJ_NAME_LISTS:
        lst = getattr(model, obj_name)
        assert isinstance(lst, list), (obj_name, lst)
        for obj in lst:
            obj_type = obj.type
            objv = getattr(modelv, obj_type.lower())
            fields = obj.raw_fields()
            card = BDFCard(fields)
            objv.add_card(card, ifile=0)

    _copy_common(model, modelv, clean_parent=clean_parent)

    if clean_parent:
        del model.nodes
        del model.coords
        del model.properties
        del model.monitor_points
        del model.suport1

        for obj_name in OBJ_NAME_DICTS:
            delattr(model, obj_name)
        for obj_name in OBJ_NAME_DICT_LISTS:
            delattr(model, obj_name)
        for obj_name in OBJ_NAME_LISTS:
            delattr(model, obj_name)
        for obj_name in SKIP:
            obj = getattr(model, obj_name)
            if obj is None:
                delattr(model, obj_name)
            elif isinstance(obj, (dict, list)) and len(obj) == 0:
                delattr(model, obj_name)

        for key in REMOVE_PARAMS:
            delattr(model, key)

        attrs = object_attributes(model)
        set_attrs = set(attrs)
        set_attrs -= set(SKIP_FINAL)
        attrs = list(set_attrs)
        attrs.sort()
        if len(attrs):
            log.warning(str(attrs))
            for key in attrs:
                obj = getattr(model, key)
                if hasattr(modelv, key):
                    objv = getattr(modelv, key)
                    if isinstance(obj, list) and isinstance(objv, list):
                        log.warning(f' {key} is a common list')
                    elif isinstance(obj, dict) and isinstance(objv, dict):
                        log.warning(f' {key} is a common dict')
                    else:
                        raise TypeError(obj)
                else:
                    if isinstance(obj, list):
                        log.warning(f' {key} is a list')
                    elif isinstance(obj, dict):
                        log.warning(f' {key} is a dict')
                    else:
                        raise TypeError(f'{key} is invalid')
 
        assert len(attrs) == 0, attrs
    modelv.setup()
    return modelv


REMOVE_ATTRSV = [
    # bool
    'replace_includes', 'force_echo_off',
    'is_strict_card_parser',
    'save_file_structure',
    'use_new_deck_parser',
    'is_superelements',
    'filter_midside_nodes',
    'run_testing_checks',
    'allow_empty_sets', 'allow_nan_area',
    'allow_nan_density', 'allow_nan_length',
    'allow_nan_mass', 'allow_nan_thickness',
    'allow_nan_volume', 'allow_overwrites_set',
    'allow_tabs',
    'apply_clip_to_desvar_range',
    # str
    'idtype', 'fdtype',
    #'active_filename',
    # bool-version
    'is_msc', 'is_mystran', 'is_nx', 'is_optistruct', 'is_zaero',
    # list[str]
    'cards_to_read',
    # list[str]?
    'include_filenames', 
    # dict
    'rsolmap_to_str',
    'reject_count',
    'card_count', # dict[str, int]
    # other
    'writer', 'zaero',
    'values_to_skip', 'special_cards',
    'sol_iline', 'sol_method',
]
IGNORE_ATTRSV = [
    # non-cards
    'nastran_format', 'active_filename', 'app',
    
    # cards
    'optimization_cards', 'mpc_cards', 'spc_cards',
    'bolt_cards',

    # elements/properties
    'element_cards', 'property_cards',
    'plot_element_cards',
    'rigid_element_cards',
    'spring_element_cards',
    'line_property_cards',
    'bar_property_cards',
    'beam_property_cards',
    'shell_element_cards', 'shell_property_cards',
    'solid_element_cards', 'solid_property_cards', 
    'nonlinear_shell_property_cards',

    # loads
    'load_cards', 
    'static_load_cards',
    'dynamic_cards', 'dynamic_load_cards',
    
    # aero
    'aero_element_cards', 'aero_load_cards', 'aero_objects',
    'aero_property_cards',  'aero_spline_cards',

    'hyperelastic_material_cards',
    'nsm_cards', 'nonstructural_mass_cards',
    'monitor_point_cards',
    'seset_cards',
    'set_cards', 'setx_cards',
    'contact_cards',
    'damper_element_cards',

    'thermal_material_dep_cards',
    'thermal_material_cards',
    'thermal_boundary_condition_cards', 'thermal_element_cards',
    # ids
    'spring_element_ids',
    'damper_element_ids',
    'bar_property_ids',
    'beam_property_ids', 
    'shell_element_ids', 'shell_property_ids',
    'solid_element_ids', 'solid_property_ids',
    'dresp_ids',
]
OBJ_NAME_DICTSV = [
    'dmig', 'dmi', 'dmij', 'dmik', 'dmiji', 'dti', 'dmiax',
    'nlparms', 'nlpcis',
    'tsteps', 'tstepnls', 'nxstrats',
    'random_tables', 'frequencies', 'dequations',
]
OBJ_NAME_LISTSV = ['mkaeros',]

class CardStorageWriter:
    def __init__(self, debug: bool=False):
        self.cards = []
        self.debug = debug
    def write(self, msg: str) -> None:
        if self.debug:
            print(msg)
            print('------------------')
        self.cards.append(msg)
    def writelines(self, lines: list[str]) -> None:
        self.cards.extend(lines)


def _v2b_grid(modelv: BDFv, model: BDF) -> None:
    grid = modelv.grid
    coord = modelv.coord

    grid = modelv.grid
    if len(grid):
        for nid, xyz, cp, cd, ps, seid in zip(
            grid.node_id.tolist(), grid.xyz, grid.cp.tolist(),
            grid.cd.tolist(), grid.ps.tolist(), grid.seid.tolist()):
            ps = '' if ps == 0 else str(ps)
            #seid = None if seid == 0 else seid
            #cd = None if cd == 0 else cd
            model.add_grid(nid, xyz, cp=cp, cd=cd, ps=ps, seid=seid)
        model.card_count['GRID'] = len(grid)

    # no card count needed for coords because add_card does it
    file_obj = CardStorageWriter()
    modelv.coord.write_file(file_obj, size=16)
    for card in file_obj.cards:
        card_lines = card.split('\n')
        card_name = _get_card_name(card_lines, '')
        fields = to_fields(card_lines, card_name)
        model.add_card(fields, card_name)
    return ['grid', 'coord']

def _v2b_prop(modelv: BDFv, model: BDF) -> None:
    pcomp = modelv.pcomp
    if len(pcomp):
        for pid, nsm, sb, ft, ge, tref, lam, z0, \
            nlayer, ilayer in zip(
                pcomp.property_id.tolist(),
                pcomp.nsm.tolist(),
                pcomp.shear_bonding.tolist(),
                pcomp.failure_theory.tolist(),
                pcomp.tref.tolist(), pcomp.ge.tolist(),
                pcomp.lam.tolist(), pcomp.z0.tolist(),
                pcomp.nlayer, pcomp.ilayer):

            ilayer0, ilayer1 = ilayer
            material_ids = pcomp.material_id[ilayer0 : ilayer1].tolist()
            thetas = pcomp.theta[ilayer0 : ilayer1].tolist()
            thicknesses = pcomp.thickness[ilayer0 : ilayer1].tolist()
            souts = pcomp.sout[ilayer0 : ilayer1].tolist()
            ft = ft if ft != '' else None
            lam = lam if lam != '' else None
            model.add_pcomp(
                pid, mids=material_ids, thetas=thetas, thicknesses=thicknesses,
                souts=souts,
                nsm=nsm, sb=sb, ft=ft,
                ge=ge, tref=tref, lam=lam, z0=z0, comment='')
        model.card_count['PCOMP'] = len(pcomp)
    return ['pcomp']

def _v2b_mat(modelv: BDFv, model: BDF) -> None:
    mat1 = modelv.mat1
    if len(mat1):
        for mid, e, g, nu, rho, alpha, tref, st, sc, ss, mcsid in zip(
            mat1.material_id, mat1.E, mat1.G, mat1.nu, mat1.rho,
            mat1.alpha, mat1.tref,
            mat1.St, mat1.Sc, mat1.Ss, mat1.mcsid):
            mat1i = model.add_mat1(
                mid, e, g, nu, rho=rho, a=alpha, tref=tref,
                St=st, Sc=sc, Ss=ss, mcsid=mcsid,
                comment='')
        model.card_count['MAT1'] = len(mat1)

    #mat2 = modelv.mat2
    #mat4 = modelv.mat4
    #mat5 = modelv.mat5
    return ['mat1']

def _v2b_rigid(modelv: BDFv, model: BDF) -> None:
    rbe2 = modelv.rbe2
    if len(rbe2):
        for eid, ind_node, ind_dof, alpha, tref, idim in zip(
            rbe2.element_id.tolist(),
            rbe2.independent_node.tolist(),
            rbe2.independent_dof.tolist(),
            rbe2.alpha.tolist(), rbe2.tref.tolist(), rbe2.idim):
            idim0, idim1 = idim
            dep = rbe2.dependent_nodes[idim0:idim1].tolist()
            gn = ind_node
            cm = str(ind_dof)
            Gmi = dep
            elem = model.add_rbe2(eid,
                           gn, # independent
                           cm, Gmi,  # dependent
                           alpha=alpha, tref=tref, comment='',
                           validate=True)
        model.card_count['RBE2'] = len(rbe2)
    return ['rbe2']

def _v2b_aero(modelv: BDFv, model: BDF) -> None:
    monpnt1 = modelv.monpnt1

    if len(monpnt1):
        for name, label, axes, comp, xyz, cp, cd in zip(
            monpnt1.name, monpnt1.label, monpnt1.axes, monpnt1.comp,
            monpnt1.xyz, monpnt1.cp, monpnt1.cd):
            model.add_monpnt1(
                name, label, axes, comp, xyz,
                cp=cp, cd=cd, comment='')
        model.card_count['MONPNT1'] = len(monpnt1)

    return ['monpnt1']

def bdf_vectorized_to_bdf(modelv: BDFv,
                          clean_parent: bool=False) -> BDF:
    model = BDF(log=modelv.log)
    log = model.log

    EXPLICITV = []

    EXPLICITV += _v2b_grid(modelv, model)
    EXPLICITV += _v2b_prop(modelv, model)
    EXPLICITV += _v2b_mat(modelv, model)
    EXPLICITV += _v2b_rigid(modelv, model)
    EXPLICITV += _v2b_aero(modelv, model)

    for obj_name in OBJ_NAME_DICTSV:
        objv_dict = getattr(modelv, obj_name)
        assert isinstance(objv_dict, dict), (obj_name, objv_dict)
        for eid, elem in objv_dict.items():
            elem_type = elem.type
            obj = getattr(model, elem_type.lower())
            fields = elem.raw_fields()
            card = BDFCard(fields)
            model.add_card(card, ifile=0)

    for obj_name in OBJ_NAME_LISTSV:
        lst = getattr(modelv, obj_name)
        assert isinstance(lst, list), (obj_name, lst)
        for obji in lst:
            obj_type = obji.type
            obj = getattr(model, obj_type.lower())  # mkaeros
            fields = obj.raw_fields()
            card = BDFCard(fields)
            model.add_card(card, ifile=0)

    _copy_common(modelv, model, clean_parent=clean_parent)

    remove_attrs = (
        IGNORE_ATTRSV + OBJ_NAME_DICTSV + 
        REMOVE_ATTRSV + OBJ_NAME_LISTSV + EXPLICITV)

    # short & slow way of writing the cards
    # relies on using writelines
    all_attrs = object_attributes(modelv)
    all_attrs = list(set(all_attrs) - set(remove_attrs))
    all_attrs.sort()
    converted_attrs = []
    for attr in all_attrs:
        obj = getattr(modelv, attr)
        if not hasattr(obj, 'n'):
            continue
        if obj.n == 0:
            continue
        
        #debug = (attr in {'pcomp'})
        debug = False
        file_obj = CardStorageWriter(debug=debug)
        size = 8 if debug else 16

        obj.write_file(file_obj, size=size)
        for card in file_obj.cards:
            card_lines = card.split('\n')
            card_name = _get_card_name(card_lines, '')
            fields = to_fields(card_lines, card_name)
            model.add_card(fields, card_name)
        converted_attrs.append(attr)
    
    if clean_parent:
        for obj_name in EXPLICITV:
            delattr(modelv, obj_name)
        for obj_name in OBJ_NAME_DICTSV:
            delattr(modelv, obj_name)
        for obj_name in REMOVE_ATTRSV:
            delattr(modelv, obj_name)
        for obj_name in OBJ_NAME_LISTSV:
            delattr(modelv, obj_name)
        for obj_name in converted_attrs:
            delattr(modelv, obj_name)

        set_attrs = set(object_attributes(modelv)) - set(IGNORE_ATTRSV)
        attrs = list(set_attrs)
        attrs.sort()

        attrs2 = []
        for attr in attrs:
            obj = getattr(modelv, attr)
            if obj is None:
                #log.warning(f' {attr} is None')
                continue
            if isinstance(obj, dict) and len(obj) == 0:
                continue
            attrs2.append(attr)
        attrs = attrs2

        if len(attrs):
            #log.warning(str(attrs))
            for attr in attrs:
                obj = getattr(modelv, attr)

                try:
                    obj.n
                except AttributeError:
                    log.warning(f' {attr} is not a card; {str(obj.__class__.__name__)}')
                    log.warning(str(obj))
    return model


def _copy_common(model: BDF | BDFv,
                 model_out: BDF | BDFv,
                 clean_parent: bool=False) -> None:
    #model_out.case_control_deck = model.case_control_deck
    COMMON = COMMON_KEYS+COMMON_SINGLETON_CARDS+COMMON_CARDS
    for key in COMMON:
        obj = getattr(model, key)
        objv = getattr(model_out, key)

        # verify they're the same base type
        # e.g., params is a dict
        if objv is not None:
            assert obj is not None, key
            assert isinstance(obj, type(objv)), (key, type(obj), type(objv))

        if key in COMMON_SINGLETON_CARDS and obj is not None:
            assert obj is not None, obj
            assert objv is None, objv
            model_out.card_count[key.upper()] = 1
        elif key in COMMON_CARDS:
            assert isinstance(obj, dict), (key, obj)
            for idi, obji in obj.items():
                assert hasattr(obji, 'type'), obji
                cls_type = obji.type
                if cls_type in model_out.card_count:
                    model_out.card_count[cls_type] += 1
                else:
                    model_out.card_count[cls_type] = 1
        #print(f'saving {key}')
        setattr(model_out, key, obj)
    if 'ENDDATA' in model.card_count:
        model_out.card_count['ENDDATA'] = model.card_count['ENDDATA']
    #print(model_out)

    if isinstance(model, BDF):
        func_to_model = bdf_to_bdf_vectorized
    else:
        func_to_model = bdf_vectorized_to_bdf

    for modeli in model.initial_superelement_models:
        modelvi = func_to_model(
            model, clean_parent=clean_parent)
        model.initial_superelement_models.append(modelvi)

    if clean_parent:
        for key in COMMON:
            delattr(model, key)
        del model.initial_superelement_models
