from pyNastran.bdf.bdf import BDF, PLOAD2

def remove_missing_loads(model: BDF) -> None:
    """
    removes loads where:
     - PLOAD2s missing elements
       - eids=1,THRU,10 with 5-10 missing becomes eids=1,THRU,4 (on separate cards)
    """
    all_eids = set(list(model.elements))
    loads2_dict = {}
    for sid, loads in model.loads.items():
        loads2 = []
        pload2s = []
        for iload, load in enumerate(loads):
            if load.type == 'PLOAD2':
                load.eids_ref = None
                pload2s.append(load)
            else:
                loads2.append(load)

        # fix the PLOAD2s dropping elements
        pload2s_new = []
        for pload2 in pload2s:
            pressure = pload2.pressure
            eids_to_add = [eid for eid in pload2.eids
                           if eid in all_eids]
            for eid in eids_to_add:
                pload2_new = PLOAD2(sid, pressure, [eid])
                pload2s_new.append(pload2_new)
        loads2 = loads2 + pload2s_new
        loads2_dict[sid] = loads2
    model.loads = loads2_dict
