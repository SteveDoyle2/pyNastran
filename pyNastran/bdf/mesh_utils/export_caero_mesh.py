"""
defines:
 - export_caero_mesh(model, caero_bdf_filename='caero.bdf', is_subpanel_model=True)

"""
from __future__ import annotations
import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.field_writer_8 import print_card_8

def export_caero_mesh(model: BDF, caero_bdf_filename: str='caero.bdf',
                      is_subpanel_model: bool=True,
                      pid_method: str='aesurf',
                      write_panel_xyz: bool=True) -> None:
    """
    Write the CAERO cards as CQUAD4s that can be visualized

    model: BDF
        a valid geometry
    caero_bdf_filename : str
        the file to write
    is_subpanel_model : bool; default=True
        True : write the subpanels as CQUAD4s
        False : write the macro elements as CQUAD4s
    pid_method : str; default='aesurf'
        'aesurf' : write the referenced AESURF as the property ID
                   main structure will be pid=1
        'caero' : write the CAERO1 as the property id
        'paero' : write the PAERO1 as the property id
    write_panel_xyz : bool; default=True
        write the following table...
        $$  CAEROID      EID       XLE       YLE       ZLE     CHORD      SPAN   XLE+C/4
        $$        1        1    0.0000    0.2500    0.0000    0.0988    0.5000    0.0247
        $$        1        2    0.0988    0.2500    0.0000    0.0988    0.5000    0.1234

    """
    if not pid_method in {'aesurf', 'caero', 'paero'}:
        raise RuntimeError(f'pid_method={pid_method!r} is not [aesurf, caero, paero]')

    inid = 1
    mid = 1
    model.log.debug('---starting export_caero_model of %s---' % caero_bdf_filename)
    with open(caero_bdf_filename, 'w') as bdf_file:
        #bdf_file.write('$ pyNastran: punch=True\n')
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')

        _write_properties(model, bdf_file, pid_method=pid_method)
        for caero_eid, caero in sorted(model.caeros.items()):
            #assert caero_eid != 1, 'CAERO eid=1 is reserved for non-flaps'
            scaero = str(caero).rstrip().split('\n')
            if is_subpanel_model:
                if caero.type == 'CAERO2':
                    continue

                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')

                #bdf_file.write("$   CAEROID       ID       XLE      YLE      ZLE     CHORD      SPAN\n")
                points, elements = caero.panel_points_elements()
                if write_panel_xyz:
                    _write_subpanel_strips(bdf_file, model, caero_eid, points, elements)

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
                    pidi = _get_subpanel_property(
                        model, caero_eid, eid2, pid_method=pid_method)
                    fields = ['CQUAD4', eid2, pidi, p1, p2, p3, p4]
                    bdf_file.write(print_card_8(fields))
                    jeid += 1
            else:
                # macro model
                if caero.type == 'CAERO2':
                    continue
                bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
                points = caero.get_points()
                npoints = 4
                for ipoint, point in enumerate(points):
                    x, y, z = point
                    bdf_file.write(print_card_8(['GRID', inid+ipoint, None, x, y, z]))

                pid = _get_subpanel_property(
                    model, caero_eid, caero_eid, pid_method=pid_method)
                p1 = inid
                p2 = inid + 1
                p3 = inid + 2
                p4 = inid + 3
                bdf_file.write(print_card_8(['CQUAD4', caero_eid, pid, p1, p2, p3, p4]))
            inid += npoints

        # aluminum
        E = 350e9 # 350 GPa
        #G = None
        nu = 0.3
        rho = 2700. # 2700 kg/m^3
        bdf_file.write(f'MAT1,%s,{E},,{nu},{rho}\n' % mid)
        bdf_file.write('ENDDATA\n')


def _write_subpanel_strips(bdf_file, model, caero_eid, points, elements):
    """writes the strips for the subpanels"""
    #bdf_file.write("$   CAEROID       ID       XLE      YLE      ZLE     CHORD      SPAN\n")
    bdf_file.write('$$\n$$ XYZ_LE is taken at the center of the leading edge; (p1+p4)/2\n$$\n')
    bdf_file.write('$$ %8s %8s %9s %9s %9s %9s %9s %9s\n' % (
        'CAEROID', 'EID', 'XLE', 'YLE', 'ZLE', 'CHORD', 'SPAN', 'XLE+C/4'))

    for i in range(elements.shape[0]):
        # The point numbers here are consistent with the CAERO1
        p1 = points[elements[i, 0], :]
        p4 = points[elements[i, 1], :]
        p2 = points[elements[i, 2], :]
        p3 = points[elements[i, 3], :]
        le = (p1 + p4)*0.5
        te = (p2 + p3)*0.5
        dy = (p4 - p1)[1]
        dz = (p4 - p1)[2]
        span = math.sqrt(dy**2 + dz**2)
        chord = te[0] - le[0]
        xqc = le + chord / 4.
        bdf_file.write("$$ %8d %8d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n" % (
            caero_eid, caero_eid+i, le[0], le[1], le[2], chord, span, xqc[0]))

def _get_subpanel_property(model: BDF, caero_id: int, eid: int, pid_method: str='aesurf') -> int:
    """gets the property id for the subpanel"""
    pid = None
    if pid_method == 'aesurf':
        for aesurf_id, aesurf in model.aesurf.items():
            aelist_id = aesurf.aelist_id1()
            aelist = model.aelists[aelist_id]
            if eid in aelist.elements:
                pid = aesurf_id
                break
    elif pid_method == 'caero':
        pid = caero_id
    elif pid_method == 'paero':
        caero = model.caeros[caero_id]
        pid = caero.pid
    else:  # pragma: no cover
        raise RuntimeError(f'pid_method={pid_method!r} is not [aesurf, caero, paero]')

    if pid is None:
        pid = 1
    return pid

def _write_properties(model: BDF, bdf_file, pid_method: str='aesurf') -> None:
    if pid_method == 'aesurf':
        _write_aesurf_properties(model, bdf_file)
    elif pid_method == 'paero':
        for paero_id, paero in sorted(model.paeros.items()):
            spaero = str(paero).rstrip().split('\n')
            bdf_file.write('$ ' + '\n$ '.join(spaero) + '\n')
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (paero_id, 1))
    elif pid_method == 'caero':
        for caero_eid, caero in sorted(model.caeros.items()):
            scaero = str(caero).rstrip().split('\n')
            bdf_file.write('$ ' + '\n$ '.join(scaero) + '\n')
            bdf_file.write('PSHELL,%s,%s,0.1\n' % (caero_eid, 1))
    else:  # pragma: no cover
        raise RuntimeError(f'pid_method={repr(pid_method)} is not [aesurf, caero, paero]')


def _write_aesurf_properties(model: BDF, bdf_file):
    aesurf_mid = 1
    for aesurf_id, aesurf in model.aesurf.items():
        #cid = aesurf.cid1

        #aesurf_mid = aesurf_id
        saesurf = str(aesurf).rstrip().split('\n')
        bdf_file.write('$ ' + '\n$ '.join(saesurf) + '\n')
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
    # dummy property
    bdf_file.write('PSHELL,%s,%s,0.1\n' % (1, 1))
    return
