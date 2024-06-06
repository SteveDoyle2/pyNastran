import numpy as np

from pyNastran.dev.bdf_vectorized.bdf import BDF
from pyNastran.utils.numpy_utils import integer_types


def bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                 starting_id_dict=None, round_ids=False, cards_to_skip=None):
    """
    Renumbers a BDF

    Parameters
    ----------
    bdf_filename : str
        a bdf_filename (string; supported) or a BDF model (BDF)
        that has been cross referenced and is fully valid (a equivalenced deck is not valid)
    bdf_filename_out : str
        a bdf_filename to write
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    starting_id_dict : dict, None (default=None)
        None : renumber everything starting from 1
        dict : {key : starting_id}
            key : str
                the key (e.g. eid, nid, cid, ...)
            starting_id : int, None
                int : the value to start from
                None : don't renumber this key
    round_ids : bool; default=False
        Should a rounding up be applied for each variable?
        This makes it easier to read a deck and verify that it's been renumbered properly.
        This only really applies when starting_id_dict is None
    cards_to_skip : list[str]; (default=None -> don't skip any cards)
        There are edge cases (e.g. FLUTTER analysis) where things can break due to
        uncross-referenced cards.  You need to disable entire classes of cards in
        that case (e.g. all aero cards).

    .. todo:: bdf_model option for bdf_filename hasn't been tested
    .. todo:: add support for subsets (e.g. renumber only a subset of nodes/elements)
    .. todo:: doesn't support partial renumbering
    .. todo:: doesn't support element material coordinate systems

    ..warning :: spoints might be problematic...check
    ..warning :: still in development, but it usually brutally crashes if it's not supported
    ..warning :: be careful of card unsupported cards (e.g. ones not read in)

    Supports
     - GRIDs
       - no superelements
     - COORDx
     - elements
        - CELASx/CONROD/CBAR/CBEAM/CQUAD4/CTRIA3/CTETRA/CPENTA/CHEXA
        - RBAR/RBAR1/RBE1/RBE2/RBE3/RSPLINE
     - properties
        - PSHELL/PCOMP/PCOMPG/PSOLID/PSHEAR/PBAR/PBARL
          PROD/PTUBE/PBEAM
     - mass
        - CMASSx/CONMx/PMASS
     - aero
       - FLFACT
       - SPLINEx
       - FLUTTER
     - partial case control
       - METHOD/CMETHOD/FREQENCY
       - LOAD/DLOAD/LSEQ/LOADSET...LOADSET/LSEQ is iffy
       - SET cards
         - nodes
         - elements
       - SPC/MPC/FLUTTER/FLFACT
    - constraints
       - SPC/SPCADD/SPCAX/SPCD
       - MPC/MPCADD
       - SUPORT/SUPORT1
    - solution control/methods
       - TSTEP/TSTEPNL
       - NLPARM
       - EIGB/EIGC/EIGRL/EIGR
    - sets
       - USET
    - other
      - tables
      - materials
      - loads/dloads

    Not Done
     - SPOINT
     - any cards with SPOINTs
       - DMIG/DMI/DMIJ/DMIJI/DMIK/etc.
       - CELASx
       - CDAMPx
     - superelements
     - aero cards
       - CAEROx
       - PAEROx
     - thermal cards?
     - optimization cards
     - SETx
     - PARAM,GRDPNT,x; where x>0
     - GRID SEID
     - case control
       - STATSUB
       - SUBCASE
       - global SET cards won't be renumbered properly

    Examples
    --------
    # Renumber Everything; Start from 1
    >>> bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     round_ids=False)

    # Renumber Everything; Start Material IDs from 100
    >>> starting_id_dict = {
        'mid' : 100,
    }
    >>> bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     starting_ids_dict=starting_ids_dict, round_ids=False)

    # Only Renumber Material IDs
    >>> starting_id_dict = {
        'cid' : None,
        'nid' : None,
        'eid' : None,
        'pid' : None,
        'mid' : 1,
        'spc_id' : None,
        'mpc_id' : None,
        'load_id' : None,
        'dload_id' : None,

        'method_id' : None,
        'cmethod_id' : None,
        'spline_id' : None,
        'table_id' : None,
        'flfact_id' : None,
        'flutter_id' : None,
        'freq_id' : None,
        'tstep_id' : None,
        'tstepnl_id' : None,
        'suport_id' : None,
        'suport1_id' : None,
        'tf_id' : None,
    }
    >>> bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     starting_ids_dict=starting_ids_dict, round_ids=False)

    """
    starting_id_dict_default = {
        'cid' : 1,
        'nid' : 1,
        'eid' : 1,
        'pid' : 1,
        'mid' : 1,
        'spc_id' : 1,
        'mpc_id' : 1,
        'load_id' : 1,
        'dload_id' : 1,

        'method_id' : 1,
        'cmethod_id' : 1,
        'spline_id' : 1,
        'table_id' : 1,
        'flfact_id' : 1,
        'flutter_id' : 1,
        'freq_id' : 1,
        'tstep_id' : 1,
        'tstepnl_id' : 1,
        'suport_id' : 1,
        'suport1_id' : 1,
        'tf_id' : 1,
    }
    # fill up starting_id_dict
    if starting_id_dict is None:
        starting_id_dict = starting_id_dict_default
    else:
        for key, value in starting_id_dict_default.items():
            if key not in starting_id_dict:
                starting_id_dict[key] = value

    nid = None
    cid = None
    eid = None
    pid = None
    mid = None
    spc_id = None
    mpc_id = None
    load_id = None
    dload_id = None
    method_id = None
    cmethod_id = None
    spline_id = None
    table_id = None
    flfact_id = None
    flutter_id = None
    freq_id = None
    tstep_id = None
    tstepnl_id = None
    suport_id = None
    suport1_id = None
    tf_id = None

    # turn them into variables
    for key, value in sorted(starting_id_dict.items()):
        #assert isinstance(key, str), key
        assert key in starting_id_dict_default, f'key={key!r} is invalid'
        if value is None:
            pass
        else:
            if not isinstance(value, integer_types):
                msg = 'key=%r value=%r must be an integer; type(value)=%s' % (
                    key, value, type(value))
                raise TypeError(msg)

        if key == 'nid':
            nid = int(value)
        elif key == 'cid':
            if value is None:
                cid = None
            else:
                cid = int(value)
        elif key == 'eid':
            eid = int(value)
        elif key == 'pid':
            if value is None:
                pid = None
            else:
                pid = int(value)
        elif key == 'mid':
            if value is None:
                mid = None
            else:
                mid = int(value)
        elif key == 'spc_id':
            spc_id = int(value)
        elif key == 'mpc_id':
            mpc_id = int(value)

        elif key == 'load_id':
            load_id = int(value)
        elif key == 'dload_id':
            dload_id = int(value)
        elif key == 'method_id':
            method_id = int(value)
        elif key == 'cmethod_id':
            cmethod_id = int(value)
        elif key == 'spline_id':
            spline_id = int(value)
        elif key == 'table_id':
            table_id = int(value)
        elif key == 'flfact_id':
            flfact_id = int(value)
        elif key == 'flutter_id':
            flutter_id = int(value)
        elif key == 'freq_id':
            freq_id = int(value)
        elif key == 'tstep_id':
            tstep_id = int(value)
        elif key == 'tstepnl_id':
            tstepnl_id = int(value)
        elif key == 'suport_id':
            suport_id = int(value)
        elif key == 'suport1_id':
            suport1_id = int(value)
        elif key == 'tf_id':
            tf_id = int(value)
        else:
            raise NotImplementedError('key=%r' % key)

    # build the maps
    #eid_map = {}
    #nid_map = {}
    #reverse_nid_map = {}
    #pid_map = {}
    #mid_map = {}
    #mpc_map = {}
    #spc_map = {}
    dload_map = {}
    load_map = {}

    cmethod_map = {}
    method_map = {}
    #flfact_map = {}
    #flutter_map = {}
    #aefact_map = {}
    #freq_map = {}
    tstep_map = {}
    tstepnl_map = {}
    suport_map = {}
    suport1_map = {}

    if isinstance(bdf_filename, str):
        model = BDF(debug=False)
        model.disable_cards(cards_to_skip)
        model.read_bdf(bdf_filename)
    else:
        model = bdf_filename

    elements = [
        model.celas1,
        model.celas2,
        model.celas3,
        model.celas4,
        #model.cdamp1,
        #model.cdamp2,
        #model.cdamp3,
        #model.cdamp4,
        model.conrod,
        model.crod,
        model.ctube,
        model.cbar,
        model.cbeam,
        model.cshear,

        model.cquad4,
        model.ctria3,
        model.cquad8,
        model.ctria6,
        #model.ctriax,
        model.ctriax6,
        model.ctetra4, model.ctetra10,
        model.cpenta6, model.cpenta15,
        #model.cpyram5, model.cpyram13,
        model.chexa8, model.chexa20,
    ]
    props = [
        model.pelas,
        #model.pdamp,
        model.ptube,
        model.ptube,

        model.pbar,
        model.pbarl,
        model.pbeam,
        model.pbeaml,

        model.pshear,
        model.pshell,
        model.pcomp,
        #model.pcompg,
        model.psolid,
        model.plsolid,
    ]
    materials = [
        model.mat1,
        #model.mat2,
        #model.mat3,
        #model.mat4,
        #model.mat5,
        model.mat8,
        #model.mat9,
        #model.mat10,
        #model.mat11,
    ]
    loads = [
        model.force,
        model.force1,
        model.force2,
        model.moment,
        model.moment1,
        model.moment2,
        model.pload,
        model.pload1,
        model.pload2,
        model.pload4,
        #model.rforce,
        #model.dload,
        #model.load,
        #model.sload,
    ]


    nid_map = {}
    #reverse_nid_map = {}
    if 'nid' in starting_id_dict and nid is not None:
        #i = nid
        #nid_map = {}
        ngrids = model.grid.n
        if ngrids:
            nids = model.grid.node_id - 1
            nid_map.update(
                {(nid + nids[i]) : (i+1) for i in range(ngrids)})
            #reverse_nid_map.update(
                #{(i+1) : (nid + nids[i]) for i in range(ngrids)}) # index to nid
        #print(min(nid_map))
        # TODO: SPOINTs
        # TODO: EPOINTs
        #if model.spoints.points:

    cid_map = {}
    if 'cid' in starting_id_dict and cid is not None:
        cids = list(model.coords.coords.keys())
        #print(cids)
        ncoords = len(cids)
        cids.sort()
        #print('cids =', cids)
        cid_map.update(
            {(cid + cids[i] - 1) : i for i in range(ncoords)})
        #print('cid_map =', cid_map)

    eid_map = {}
    if 'eid' in starting_id_dict and pid is not None:

        eids = [element.element_id for element in elements if element.n]
        eids.append(list(model.rigid_elements.keys()))
        eids = np.hstack(eids)
        neids = len(eids)
        eids.sort()
        eid_map.update(
            {eids[i] : (i + eid) for i in range(neids)})

    pid_map = {}
    if 'pid' in starting_id_dict and pid is not None:
        pids = np.hstack([prop.property_id for prop in props if prop.n])
        npids = len(pids)
        pids.sort()
        pid_map.update(
            {pids[i] : (i + pid) for i in range(npids)})

    mid_map = {}
    if 'mid' in starting_id_dict and mid is not None:
        mids = np.hstack([mat.material_id for mat in materials])
        nmids = len(mids)
        mids.sort()
        mid_map.update(
            {mids[i] : (mid + i) for i in range(nmids)})

    load_map = {}
    if 'load_id' in starting_id_dict and load_id is not None:
        loadsi = [load.load_id for load in loads if load.n]
        if loadsi:
            #print(loadsi)
            load_ids = np.unique(np.hstack(loadsi))
            del loadsi
            nload_ids = len(load_ids)
            load_ids.sort()
            load_map.update(
                {load_ids[i] : (load_id + i) for i in range(nload_ids)})

    spc_map = {}
    if 'spc_id' in starting_id_dict and spc_map is not None:
        spcs = [model.spc, model.spc1, model.spcadd]
        spcsi = []
        for spc_obj in spcs:
            #print(spc_obj)
            #spc_ids = [spc.spc_id for spc_id, spc in spc_obj.items() if spc.n]
            spc_ids = spc_obj.keys()
            if spc_ids:
                spcsi.append(spc_ids)
        #spcsi = [spc_id for spc_id in spcs]
        #print('spcsi =', spcsi)
        #asdf
        if spcsi:
            #print(loadsi)
            spc_ids = np.unique(np.hstack(spcsi))
            #del spcsi
            nspc_ids = len(spc_ids)
            spc_ids.sort()
            spc_map.update(
                {spc_ids[i] : (spc_id + i) for i in range(nspc_ids)})
            #print('spc_ids =', spc_ids)

    caero_map = {}
    #if model.caeros:
        #caeros = model.caeros.keys()
        #caeros.sort()
        #ncaeros = len(caeros)
        #caero_map = {caeros[i] : (i+1) for i in range(ncaeros)}

    paero_map = {}
    if model.paeros:
        paeros = model.paeros.keys()
        paeros.sort()
        npaeros = len(paeros)
        paero_map = {paeros[i] : (i+1) for i in range(npaeros)}

    aefact_map = {}
    if model.aefacts:
        aefacts = model.aefacts.keys()
        aefacts.sort()
        naefacts = len(aefacts)
        aefact_map = {aefacts[i] : (i+1) for i in range(naefacts)}

    spline_map = {}
    #if model.splines:
        #splines = model.splines.keys()
        #splines.sort()
        #nsplines = len(splines)
        #spline_map = {splines[i] : (i+1) for i in range(nsplines)}

    set_map = {}
    if model.sets:
        sets = model.sets.keys()
        sets.sort()
        nsets = len(sets)
        set_map = {sets[i] : (i+1) for i in range(nsets)}

    #load_map = {}
    maps = {
        'node' : nid_map,
        'coord' : cid_map,
        'property' : pid_map,
        'element' : eid_map,
        'load' : load_map,
        'material' : mid_map,

        'caero' : caero_map, # ???
        'paero' : paero_map, # PAEROx
        'aefact' : aefact_map, # ???
        'spline' : spline_map, # SPLINE1-SPLINE5
        'set' : set_map,
    }

    model.grid.update(maps)
    model.coords.update(maps)
    for elem in elements:
        elem.update(maps)

    rigid_elements2 = {}
    for eid, elem in model.rigid_elements.values():
        eid2 = eid_map[eid]
        rigid_elements2[eid2] = eid_map[eid]
        elem.update(maps)

    for prop in props:
        prop.update(maps)
    for mat in materials:
        mat.update(maps)
    for spc_dict in spcs:
        for spc_id, spc in spc_dict.values():
            spc.update(maps)

    if model.aero is not None:
        model.aero.update(maps)
    if model.aeros is not None:
        model.aeros.update(maps)

    for caero in model.caeros.values():
        caero.update(maps)
    for spline in model.splines.values():
        spline.update(model, maps)
    for flutter in model.flutters.values():
        flutter.update(maps)
    for flfact in model.flfacts.values():
        flfact.update(maps)
    for flutter in model.flutters.values():
        flutter.update(maps)

    for desvar in model.desvars.values():
        desvar.update(maps)
    for dconstr in model.dconstrsvalues():
        dconstr.update(maps)
    for dresp in model.dresps.values():
        dresp.update(maps)
    for dconadd in model.dconadds.values():
        dconadd.update(maps)
    for dvgrid in model.dvgrids.values():
        dvgrid.update(maps)
    for dvcrel in model.dvcrels.values():
        dvcrel.update(maps)
    for dvmrel in model.dvmrels.values():
        dvmrel.update(maps)
    for dvprel in model.dvprels.values():
        dvprel.update(maps)

    model.darea.update(maps)
    model.dphase.update(maps)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double,
                        interspersed=False)
    return model
