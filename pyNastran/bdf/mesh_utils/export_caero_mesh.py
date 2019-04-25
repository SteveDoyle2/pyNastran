"""
defines:
 - export_caero_mesh(model, caero_bdf_filename='caero.bdf', is_subpanel_model=True)

"""
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.field_writer_8 import print_card_8

def export_caero_mesh(model, caero_bdf_filename='caero.bdf', is_subpanel_model=True):
    # type: (BDF, str, bool) -> None
    """write the CAERO cards as CQUAD4s that can be visualized"""
    export_caero_model_aesurf(model, caero_bdf_filename='caero.bdf', is_subpanel_model=is_subpanel_model)

#def export_caero_model_base(model, caero_bdf_filename='caero.bdf', is_subpanel_model=True):
    ## type: (BDF, str, bool) -> None
    #inid = 1
    #mid = 1
    #model.log.debug('---starting export_caero_model of %s---' % caero_bdf_filename)
    #with open(caero_bdf_filename, 'w') as bdf_file:
        #bdf_file.write('CEND\n')
        #bdf_file.write('BEGIN BULK\n')

        #for caero_eid, caero in sorted(model.caeros.items()):
            #if caero.type == 'CAERO2':
                #continue

            ##assert caero_eid != 1, 'CAERO eid=1 is reserved for non-flaps'
            #scaero = str(caero).rstrip().split('\n')
            #bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
            #points, elements = caero.panel_points_elements()
            #npoints = points.shape[0]
            ##nelements = elements.shape[0]
            #for ipoint, point in enumerate(points):
                #x, y, z = point
                #bdf_file.write(print_card_8(['GRID', inid+ipoint, None, x, y, z]))

def export_caero_model_aesurf(model, caero_bdf_filename='caero.bdf', is_subpanel_model=True):
    # type: (BDF, str, bool) -> None
    inid = 1
    mid = 1
    model.log.debug('---starting export_caero_model of %s---' % caero_bdf_filename)
    with open(caero_bdf_filename, 'w') as bdf_file:
        #bdf_file.write('$ pyNastran: punch=True\n')
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')
        #if is_subpanel_model:
        for aesurf_id, unused_aesurf in model.aesurf.items():
            #cid = aesurf.cid1

            #aesurf_mid = aesurf_id
            aesurf_mid = 1
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (aesurf_id, aesurf_mid))
            #print(cid)
            #ax, ay, az = cid.i
            #bx, by, bz = cid.j
            #cx, cy, cz = cid.k
            #bdf_file.write('CORD2R,%s,,%s,%s,%s,%s,%s,%s\n' % (
                #cid, ax, ay, az, bx, by, bz))
            #bdf_file.write(',%s,%s,%s\n' % (cx, cy, cz))
            #print(cid)
            #aesurf.elements

        for caero_eid, caero in sorted(model.caeros.items()):
            #assert caero_eid != 1, 'CAERO eid=1 is reserved for non-flaps'
            scaero = str(caero).rstrip().split('\n')
            if is_subpanel_model:
                if caero.type == 'CAERO2':
                    continue

                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
                points, elements = caero.panel_points_elements()
                npoints = points.shape[0]
                #nelements = elements.shape[0]
                for ipoint, point in enumerate(points):
                    x, y, z = point
                    bdf_file.write(print_card_8(['GRID', inid+ipoint, None, x, y, z]))

                #pid = caero_eid
                #mid = caero_eid
                jeid = 0
                for elem in elements + inid:
                    p1, p2, p3, p4 = elem
                    eid2 = jeid + caero_eid
                    pidi = _get_subpanel_property(model, eid2)
                    fields = ['CQUAD4', eid2, pidi, p1, p2, p3, p4]
                    bdf_file.write(print_card_8(fields))
                    jeid += 1
            else:
                if caero.type == 'CAERO2':
                    continue
                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
                points = caero.get_points()
                npoints = 4
                for ipoint, point in enumerate(points):
                    x, y, z = point
                    bdf_file.write(print_card_8(['GRID', inid+ipoint, None, x, y, z]))

                pid = _get_subpanel_property(model, caero_eid)
                p1 = inid
                p2 = inid + 1
                p3 = inid + 2
                p4 = inid + 3
                bdf_file.write(print_card_8(['CQUAD4', caero_eid, pid, p1, p2, p3, p4]))
            inid += npoints
        bdf_file.write('PSHELL,%s,%s,0.1\n' % (1, 1))
        bdf_file.write('MAT1,%s,3.0E7,,0.3\n' % mid)
        bdf_file.write('ENDDATA\n')

def _get_subpanel_property(model, eid):
    """gets the property id for the subpanel"""
    pidi = None
    for aesurf_id, aesurf in model.aesurf.items():
        aelist_id = aesurf.aelist_id1()
        aelist = model.aelists[aelist_id]
        if eid in aelist.elements:
            pidi = aesurf_id
            break
    if pidi is None:
        #pidi = pid
        pidi = 1
    return pidi
