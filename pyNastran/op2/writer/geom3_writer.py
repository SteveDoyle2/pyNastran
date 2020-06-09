from struct import pack, Struct
from collections import defaultdict

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header, write_header_nvalues

def write_geom3(op2, op2_ascii, obj, endian=b'<', nastran_format='nx'):
    if not hasattr(obj, 'loads') and not hasattr(obj, 'load_combinations'):
        return
    loads_by_type = defaultdict(list)
    for unused_load_id, loads in obj.loads.items():
        for load in loads:
            loads_by_type[load.type].append(load)
    for unused_load_id, loads in obj.load_combinations.items():
        for load in loads:
            loads_by_type[load.type].append(load)
    for unused_load_id, load in obj.tempds.items():
        loads_by_type[load.type].append(load)
    # pedge, pface

    # return if no supported cards are found
    wrong_tables = [
        # see EDT
        'DEFORM', 'CLOAD', # these are be in the
    ]
    cards_to_skip = [
        'TEMPRB',
    ]
    supported_cards = [
        'FORCE', 'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1', 'MOMENT2',
        'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4', 'PLOADX1',
        'GRAV', 'SLOAD', 'RFORCE',
        #'ACCEL',
        'ACCEL1',
        'TEMP', 'TEMPP1', 'QBDY1', 'QBDY2', 'QBDY3', 'QVOL',
        'LOAD',
        'TEMPD',
        'QHBDY',
    ]
    for key in wrong_tables:
        if key not in loads_by_type:
            continue
        del loads_by_type[key]

    is_loads = False
    for load_type in sorted(loads_by_type.keys()):
        if load_type in supported_cards:
            is_loads = True
            continue
            #break
        elif load_type in cards_to_skip:
            obj.log.warning('skipping GEOM3-%s' % load_type)
            continue
    #else:
        #return

    if not is_loads:
        return
    write_geom_header(b'GEOM3', op2, op2_ascii)

    itable = -3
    for load_type, loads in sorted(loads_by_type.items()):
        if load_type in ['SPCD']: # not a GEOM3 load
            continue
        elif load_type in cards_to_skip:
            obj.log.warning('skipping GEOM3-%s' % load_type)
            continue

        #elif load_type not in supported_cards:
            #continue
        #print('GEOM3', itable, load_type)
        try:
            nbytes = write_card(op2, op2_ascii, load_type, loads, endian, obj.log,
                                nastran_format=nastran_format)
        except:  # pragma: no cover
            obj.log.error('failed GEOM3-%s' % load_type)
            raise
        op2.write(pack('i', nbytes))

        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)

    #-------------------------------------

def write_card(op2, op2_ascii, load_type, loads, endian, log,
               nastran_format: str='nx'):
    nloads = len(loads)
    if load_type == 'FORCE':
        key = (4201, 42, 18)
        nfields = 7
        spack = Struct(endian + b'3i 4f')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            data = [load.sid, load.node_id, load.Cid(), load.mag] + list(load.xyz)
            op2_ascii.write('  FORCE data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    elif load_type == 'FORCE1':
        key = (4001, 40, 20)
        nfields = 5
        spack = Struct(endian + b'iifii')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            #(sid, node, mag, n1, n2) = out
            data = [load.sid, load.node_id, load.mag, load.G1(), load.G2()]
            op2_ascii.write('  FORCE1 data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    elif load_type == 'FORCE2':
        key = (4101, 41, 22)
        nfields = 7
        spack = Struct(endian + b'iif4i')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            #(sid, node_id, mag, n1, n2, n3, n4) = out
            data = [load.sid, load.node_id, load.mag,
                    load.G1(), load.G2(), load.G3(), load.G4()]
            op2_ascii.write('  FORCE2 data=%s\n' % str(data))
            op2.write(spack.pack(*data))

    elif load_type == 'MOMENT':
        key = (4801, 48, 19)
        nfields = 7
        spack = Struct(endian + b'3i 4f')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            data = [load.sid, load.node_id, load.Cid(), load.mag] + list(load.xyz)
            op2_ascii.write('  MOMENT data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    elif load_type == 'MOMENT1':
        key = (4601, 46, 21)
        nfields = 5
        spack = Struct(endian + b'iifii')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            #(sid, node, mag, n1, n2) = out
            data = [load.sid, load.node_id, load.mag, load.G1(), load.G2()]
            op2_ascii.write('  MOMENT1 data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    elif load_type == 'MOMENT2':
        key = (4701, 47, 23)
        nfields = 7
        spack = Struct(endian + b'iif4i')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            #(sid, node_id, mag, n1, n2, n3, n4) = out
            data = [load.sid, load.node_id, load.mag,
                    load.G1(), load.G2(), load.G3(), load.G4()]
            op2_ascii.write('  MOMENT2 data=%s\n' % str(data))
            op2.write(spack.pack(*data))

    elif load_type == 'GRAV':
        nbytes = _write_grav(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'PLOAD':
        nbytes = _write_pload(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'PLOAD1':
        nbytes = _write_pload1(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'PLOAD2':
        nbytes = _write_pload2(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'PLOAD4': # msc
        nbytes = _write_pload4(load_type, loads, op2, op2_ascii,
                               endian, nastran_format=nastran_format)

    elif load_type == 'PLOADX1':
        key = (7309, 73, 351)
        nfields = 7
        spack = Struct(endian + b'2i2f iif')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            #(sid, eid, pa, pb, ga, gb, theta) = out
            data = [load.sid, load.eid, load.pa, load.pb, load.ga, load.gb, load.theta]
            op2_ascii.write('  PLOADX1 data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    elif load_type == 'SLOAD':
        nbytes = _write_sload(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'RFORCE':
        nbytes = _write_rforce(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'TEMP':
        nbytes = _write_temp(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'QVOL':
        nbytes = _write_qvol(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'QBDY1':
        nbytes = _write_qbdy1(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'QBDY2':
        nbytes = _write_qbdy2(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'QBDY3':
        nbytes = _write_qbdy3(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'QHBDY':
        nbytes = _write_qhbdy(load_type, loads, nloads, op2, op2_ascii, endian, log)

    elif load_type == 'TEMPP1':
        nbytes = _write_tempp1(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'TEMPD':
        nbytes = _write_tempd(load_type, loads, nloads, op2, op2_ascii, endian)


    #elif load_type == 'ACCEL1':
        #key = (7401,74,601)
        #nfields = 3
        #spack = Struct(endian + b'2i f')
        #nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        #for load in loads:
            #1 SID    I Load set identification number
            #2 CID    I Coordinate system identification number
            #3 A     RS Acceleration vector scale factor
            #4 N(3)  RS Components of a vector coordinate system defined by CID
            #7 GRIDID I Grid ID or THRU or BY code
            #Words 7 repeats until (-1) occurs.
            #for nid, mag in zip(load.nodes, load.mags):
                #data = [load.sid, nid, mag]
                #op2_ascii.write('  SLOAD data=%s\n' % str(data))
                #op2.write(spack.pack(*data))
    elif load_type == 'LOAD':
        nbytes = _write_load(load_type, loads, op2, op2_ascii, endian)
    elif load_type == 'LSEQ':
        nbytes = _write_lseq(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'ACCEL':
        nbytes = _write_accel(load_type, loads, nloads, op2, op2_ascii, endian)
    elif load_type == 'ACCEL1':
        nbytes = _write_accel1(load_type, loads, op2, op2_ascii, endian)

    else:  # pragma: no cover
        load0 = loads[0]
        raise NotImplementedError(load0)
    return nbytes

def _write_pload4(load_type, loads, op2, op2_ascii, endian, nastran_format='nx'):
    """writes the PLOAD4s"""
    key = (7209, 72, 299)

    nloads = 0
    for load in loads:
        nloads += len(load.eids)

    if nastran_format == 'msc':
        nfields = 16
        spack = Struct(endian + b'2i 4f 3i 3f 8s 8s')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            #surf_or_line = surf_or_line.rstrip().decode('latin1')
            #line_load_dir = line_load_dir.rstrip().decode('latin1')
            #if line_load_dir == '':
                ## TODO: not 100%
                #line_load_dir = 'NORM'

            ## forces NX pload4 function to get called if it should be
            #assert surf_or_line in ['SURF', 'LINE']
            #assert line_load_dir in ['LINE', 'X', 'Y', 'Z', 'TANG', 'NORM'], 'line_load_dir=%r' % line_load_dir

            #self.pressures = np.asarray(pressures, dtype='float64')
            #self.nvector = nvector
            #self.surf_or_line = surf_or_line
            #self.line_load_dir = line_load_dir
            pressures = list(load.pressures)
            g1 = load.g1 if load.g1 is not None else 0
            g34 = load.g34 if load.g34 is not None else 0
            cid = load.cid if load.cid is not None else 0
            nids_cid = [g1, g34, cid]
            nvector = list(load.nvector)
            assert len(load.pressures) == 4, load.pressures
            assert None not in nids_cid, nids_cid

            pnn = pressures + nids_cid + nvector
            for eid in load.eids:
                #(sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3, surf_or_line, line_load_dir) = out
                surf_or_line = load.surf_or_line.encode('ascii')
                line_load_dir = load.line_load_dir.encode('ascii')
                data = [load.sid, eid] + pnn + [surf_or_line, line_load_dir]
                assert None not in data, data
                op2_ascii.write('  PLOAD4 data=%s\n' % str(data))
                op2.write(spack.pack(*data))
    elif nastran_format == 'nx':
        #Word Name Type Description
        #1 SID          I Load set identification number
        #2 EID          I Element identification number
        #3 P(4)        RS Pressures
        #7 G1           I Grid point identification number at a corner of the face
        #8 G34          I Grid point identification number at a diagonal from G1 or CTETRA corner
        #9  CID         I Coordinate system identification number
        #10 N(3)       RS Components of a vector coordinate system defined by CID
        nfields = 12
        spack = Struct(endian + b'2i 4f 3i 3f')
        nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
        for load in loads:
            pressures = list(load.pressures)
            g1 = load.g1 if load.g1 is not None else 0
            g34 = load.g34 if load.g34 is not None else 0
            cid = load.cid if load.cid is not None else 0
            nids_cid = [g1, g34, cid]
            nvector = list(load.nvector)
            assert len(load.pressures) == 4, load.pressures
            assert None not in nids_cid, nids_cid

            pnn = pressures + nids_cid + nvector
            for eid in load.eids:
                #(sid, eid, p1, p2, p3, p4, g1, g34, cid, n1, n2, n3) = out
                data = [load.sid, eid] + pnn

                assert None not in data, data
                op2_ascii.write('  PLOAD4 data=%s\n' % str(data))
                op2.write(spack.pack(*data))
    return nbytes

def _get_nloads_from_eids(loads):
    nloads = 0
    for load in loads:
        nloads += len(load.eids)
    return nloads

def _get_nloads_from_elements(loads):
    nloads = 0
    for load in loads:
        nloads += len(load.elements)
    return nloads

def _write_accel(load_type, loads, nloads, op2, op2_ascii, endian):
    """
    ACCEL(7401, 74, 601) - NX  (uses CHAR1)
    ACCEL(11302,113,600) - MSC (uses CHAR4)

    Word Name Type Description
    1 SID     I Load set identification number
    2 CID     I Coordinate system identification number
    3 N(3)   RS Components of a vector coordinate system defined by CID
    6 DIR CHAR1 Component direction of acceleration variation
    7 LOCi   RS Location along direction DIR in coordinate system
    8 VALi   RS The load scale factor associated with location LOCi
    Words 7 through 8 repeat until (-1,-1) occurs.
    """
    raise NotImplementedError(load_type)
    #key = (7401, 74, 601)
    #fmt = endian
    #data = []
    #for load in loads:
        #nids = load.node_ids
        #fmt += b'iif%ii' % (len(nids) + 1)
        #datai = [load.sid, load.Cid(), load.scale] + nids + [-1]
        #op2_ascii.write('  LOAD data=%s\n' % str(datai))
        #data += datai
    #nfields = len(data)
    #nbytes = write_header_nvalues(load_type, nfields, key, op2, op2_ascii)
    #op2.write(pack(fmt, *data))
    #return nbytes

def _write_accel1(load_type, loads, op2, op2_ascii, endian):
    """
    ACCEL1(7501, 75, 602) - NX
    ACCEL1(11402,114,601) - MSC

    1 SID    I Load set identification number
    2 CID    I Coordinate system identification number
    3 A     RS Acceleration vector scale factor
    4 N(3)  RS Components of a vector coordinate system defined by CID
    7 GRIDID I Grid ID or THRU or BY code
    Words 7 repeats until (-1) occurs.
    NX/MSC
    """
    key = (7501, 75, 602)
    fmt = endian
    data = []
    for load in loads:
        nids = load.node_ids
        fmt += b'iif%ii' % (len(nids) + 1)
        datai = [load.sid, load.Cid(), load.scale] + nids + [-1]
        op2_ascii.write('  LOAD data=%s\n' % str(datai))
        data += datai
    nfields = len(data)
    nbytes = write_header_nvalues(load_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *data))
    return nbytes

def _write_grav(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the GRAVs"""
    key = (4401, 44, 26)
    nfields = 7
    #(sid, cid, a, n1, n2, n3, mb) = out
    spack = Struct(endian + b'ii4fi')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        data = [load.sid, load.Cid(), load.scale] + list(load.N) + [load.mb]
        op2_ascii.write('  GRAV data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_load(load_type, loads, op2, op2_ascii, endian):
    """writes the LOADs"""
    key = (4551, 61, 84)
    fmt = endian
    data = []
    for load in loads:
        nscales = len(load.scale_factors)
        fmt += b'if' + b'fi' * nscales + b'ii'
        datai = [load.sid, load.scale, ]
        for scale, load_id in zip(load.scale_factors, load.load_ids):
            datai.append(scale)
            datai.append(load_id)
        datai += [-1, -1]
        op2_ascii.write('  LOAD data=%s\n' % str(datai))
        data += datai
    nfields = len(data)
    nbytes = write_header_nvalues(load_type, nfields, key, op2, op2_ascii)
    op2.write(pack(fmt, *data))
    return nbytes

def _write_lseq(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the LSEQs"""
    key = (3609, 36, 188)
    nfields = 5
    spack = Struct(endian + b'5i')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, darea, load_id, temperature_id, undef) = out
        tid = 0 if load.tid is None else load.tid
        datai = [load.sid, load.excite_id, load.lid, tid, 0]
        op2_ascii.write('  LSEQ data=%s\n' % str(datai))
        op2.write(spack.pack(*datai))
    return nbytes

def _write_pload(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the PLOADs"""
    key = (5101, 51, 24)
    nfields = 6
    spack = Struct(endian + b'i f 4i')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, pressure, n1, n2, n3, n4) = out
        nids = list(load.node_ids)
        if len(nids) == 3:
            nids.append(0)
        data = [load.sid, load.pressure] + nids
        op2_ascii.write('  PLOAD data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_pload1(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the PLOAD1s"""
    key = (6909, 69, 198)
    nfields = 8
    spack = Struct(endian + b'4i4f')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, eid, load_type, scale, x1, p1, x2, p2) = out
        load_typei = load.valid_types.index(load.load_type) + 1 # 1-based
        scale = load.valid_scales.index(load.scale) + 1
        data = [load.sid, load.eid, load_typei, scale,
                load.x1, load.p1, load.x2, load.p2]
        op2_ascii.write('  PLOAD1 data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_pload2(load_type, loads, op2, op2_ascii, endian):
    """writes the PLOAD2s"""
    key = (6802, 68, 199)
    nfields = 3
    spack = Struct(endian + b'ifi')

    nloads = _get_nloads_from_eids(loads)
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, p, eid) = out
        for eid in load.eids:
            data = [load.sid, load.pressure, eid]
            op2_ascii.write('  PLOAD2 data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    return nbytes

def _write_qbdy1(load_type, loads, op2, op2_ascii, endian):
    """writes the QBDY1s"""
    key = (4509, 45, 239)
    nfields = 3
    spack = Struct(endian + b'ifi')

    nloads = _get_nloads_from_eids(loads)
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, q0, eid) = out
        for eid in load.eids:
            data = [load.sid, load.qflux, eid]
            op2_ascii.write('  QBDY1 data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    return nbytes

def _write_qbdy2(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the QBDY2s"""
    key = (4909, 49, 240)
    nfields = 10
    spack = Struct(endian + b'ii8f')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, eid, q1, q2, q3, q4, q5, q6, q7, q8) = out
        qflux = list(load.qfluxs)
        nflux = len(qflux)
        if nflux < 8:
            qflux = qflux + [0.] * (8 - nflux)
        data = [load.sid, load.eid] + qflux
        op2_ascii.write('  QBDY2 data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_qbdy3(load_type, loads, op2, op2_ascii, endian):
    """writes the QBDY3s"""
    key = (2109, 21, 414)
    nfields = 4
    spack = Struct(endian + b'ifii')

    nloads = _get_nloads_from_eids(loads)
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, q0, cntrlnd, eid) = out
        for eid in load.eids:
            data = [load.sid, load.q0, load.cntrlnd, eid]
            op2_ascii.write('  QBDY3 data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    return nbytes

def _write_qhbdy(load_type, loads, nloads, op2, op2_ascii, endian, log):
    """writes the QHBDYs"""
    key = (4309, 43, 233)
    nfields = 12
    spack = Struct(endian + b'2i 2f 8i')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #1 SID   I Load set identification number
        #2 FLAG  I Face type
        #3 Q0   RS Magnitude of thermal flux into face
        #4 AF   RS Area factor
        #5 G(8)  I Grid point identification numbers
        #print(load.get_stats())

        nids = [0] * 8
        #if flag == 1:
            #flag = 'POINT'
        #elif flag == 2:
            #flag = 'LINE'
        #elif flag == 5:
            #flag = 'AREA4'
        #elif flag == 9:
            #flag = 'AREA8'

        if load.flag == 'POINT':
            # C:\NASA\m4\formats\git\examples\move_tpl\ex8a.op2 - flag=1
            flag = 1 # 0?
            nnodes = 1
        elif load.flag == 'LINE':
            flag = 2
            nnodes = 2
        elif load.flag == 'REV':
            flag = 3
            nnodes = 2
        elif load.flag == 'AREA4':
            flag = 5
            nnodes = 4
        elif load.flag == 'AREA8':
            flag = 9
            nnodes = 8
        else:  # pragma: no cover
            print(load.get_stats())
            raise NotImplementedError(load.get_stats())

        grids = load.grids
        for i in range(nnodes):
            nids[i] = grids[i]
            if nids[i] <= 0:
                log.warning(f'QHBDY: nids[{i}]={nids[i]} nids={nids}')
            #assert nids[i] > 0, f'QHBDY: nids[{i}]={nids[i]} nids={nids}'

        data = [
            load.sid, flag, load.q0,
            0.0 if load.af is None else load.af,
            ] + nids
        op2_ascii.write('  QHBDY data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_qvol(load_type, loads, op2, op2_ascii, endian):
    """writes the QVOLs"""
    key = (2309, 23, 416)
    nfields = 4
    spack = Struct(endian + b'if2i')

    nloads = _get_nloads_from_elements(loads)
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #(sid, qvol, cntrlnd, eid) = out
        for eid in load.elements:
            data = [load.sid, load.qvol, load.control_point, eid]
            op2_ascii.write('  QVOL data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    return nbytes

def _get_nloads_from_temperatures(loads):
    nloads = 0
    for load in loads:
        nloads += len(load.temperatures)
    return nloads

def _write_sload(load_type, loads, op2, op2_ascii, endian):
    """writes the SLOADs"""
    key = (5401, 54, 25)
    data = []
    for load in loads:
        for nid, mag in zip(load.nodes, load.mags):
            #(sid, nid, scale_factor) = out
            datai = [load.sid, nid, mag]
            op2_ascii.write('  SLOAD data=%s\n' % str(datai))
            data += datai
    nfields = len(data)
    nloads = nfields // 3
    spack = Struct(endian + b'iif' * nloads)
    nbytes = write_header_nvalues(load_type, nfields, key, op2, op2_ascii)
    op2.write(spack.pack(*data))
    return nbytes

def _write_rforce(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the RFORCEs"""
    key = (5509, 55, 190)
    nfields = 10
    spack = Struct(endian + b'3i 4f ifi')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #self.idrf = idrf
        #sid, nid, cid, a, r1, r2, r3, method, racc, mb = data
        #scale = 1.0
        data = [load.sid, load.nid, load.cid, load.scale] + load.r123 + [
            load.method, load.racc, load.mb]
        assert load.idrf == 0, load
        op2_ascii.write('  RFORCE data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_temp(load_type, loads, op2, op2_ascii, endian):
    """writes the TEMPs"""
    key = (5701, 57, 27)
    nfields = 3
    spack = Struct(endian + b'iif')

    nloads = _get_nloads_from_temperatures(loads)
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        for nid, temp in sorted(load.temperatures.items()):
            #(sid, g, T) = out
            data = [load.sid, nid, temp]
            op2_ascii.write('  TEMP data=%s\n' % str(data))
            op2.write(spack.pack(*data))
    return nbytes

def _write_tempd(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the TEMPDs"""
    key = (5641, 65, 98)
    nfields = 2
    spack = Struct(endian + b'if')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #sid, T = data
        data = [load.sid, load.temperature]
        op2_ascii.write('  TEMPD data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes

def _write_tempp1(load_type, loads, nloads, op2, op2_ascii, endian):
    """writes the TEMPP1s"""
    key = (8109, 81, 201)
    nfields = 6
    spack = Struct(endian + b'2i 4f')
    nbytes = write_header(load_type, nfields, nloads, key, op2, op2_ascii)
    for load in loads:
        #sid, eid, t, tprime, ts1, ts2 = data
        data = [load.sid, load.eid, load.tbar, load.tprime] + load.t_stress
        op2_ascii.write('  TEMPP1 data=%s\n' % str(data))
        op2.write(spack.pack(*data))
    return nbytes
