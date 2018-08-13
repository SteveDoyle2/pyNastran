def get_eid_dt_from_eid_device(eid_device, dt, sort_method):
    """common unvectorized method for transforming SORT2 data into SORT1"""
    if sort_method == 1:
        eid = eid_device // 10
        #print("SORT1 dt=%s eid_device=%s eid=%s" % (dt, eid_device, eid))
    else:
        eid, dt = dt, eid_device
    return eid, dt
