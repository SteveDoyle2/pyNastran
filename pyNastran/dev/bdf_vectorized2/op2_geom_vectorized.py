import struct
import numpy as np
from pyNastran.op2.op2_geom import OP2GeomCommon
from pyNastran.dev.bdf_vectorized2.bdf_vectorized import BDF

class OP2Geom(OP2GeomCommon, BDF):
    """creates an vectorized interface for the OP2 and BDF classes"""
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None, mode='msc'):
        """
        Initializes the OP2 object

        Parameters
        ----------
        make_geom : bool; default=False
            reads the BDF tables
        debug : bool; default=False
            enables the debug log and sets the debug in the logger
        log: log()
            a logging object to write debug messages to
            (.. seealso:: import logging)
        debug_file : default=None -> no debug
            sets the filename that will be written to
        mode : str; default='msc'
            {msc, nx}

        """
        BDF.__init__(self, debug=debug, log=log, mode=mode)
        OP2GeomCommon.__init__(self, debug=debug, log=log, debug_file=debug_file, mode=mode)

    def read_op2(self, op2_filename=None, combine=True, build_dataframe=False,
                 skip_undefined_matrices=False, encoding=None):
        OP2GeomCommon.read_op2(
            self, op2_filename=op2_filename, combine=combine, build_dataframe=build_dataframe,
            skip_undefined_matrices=skip_undefined_matrices,
            encoding=encoding)
        #assert len(self.elements) == 0, self.card_count
        #assert len(self.elements) == 0, self.elements

    def _read_grid(self, data, n):
        """(4501,45,1) - the marker for Record 17"""
        ntotal = 32
        nentries = (len(data) - n) // ntotal
        #datan = data[n:]
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #  0    1   2   3   4   5   6   7
        #(nid, cp, x1, x2, x3, cd, ps, seid) = out
        #nid, cp, cd, seid = ints[:, [0, 1, 5, 7]]
        nid = ints[:, 0]
        cp = ints[:, 1]
        cd = ints[:, 5]
        seid = ints[:, 7]
        xyz = floats[:, [2, 3, 4]]
        ps = np.full(nentries, '', dtype='|U8') # TODO: fake

        grid = self.grid
        grid.nid = nid
        grid.cp = cp
        grid.cd = cd
        grid.xyz = xyz
        grid.ps = ps
        grid.seid = seid

        nfailed = 0
        #isave = np.where(nid < 10000000)
        self.increase_card_count('GRID', nentries - nfailed)
        return len(data)

    def _read_ctetra(self, data, n):
        """
        CTETRA(5508,55,217)      - the marker for Record 88
        CTETPR(7609,76,9993)     - the marker for Record 89
        CTETR10F(16600,166,9999) - the marker for Record 90
        CTETR4FD(16100,161,9999) - the marker for Record 91
        """
        ntotal = 48
        nentries = (len(data) - n) // ntotal
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        extra_nodes = nodes[:, 4:]
        assert extra_nodes.shape[1] == 6, extra_nodes.shape
        max_nid = extra_nodes.min(axis=1)
        izero = np.where(max_nid == 0)[0]
        inonzero = np.where(max_nid > 0)[0]
        if len(izero):
            self.ctetra4.eid = eid[izero]
            self.ctetra4.pid = pid[izero]
            self.ctetra4.nids = nodes[izero, :4]
        if len(inonzero):
            self.ctetra10.eid = eid[inonzero]
            self.ctetra10.pid = pid[inonzero]
            self.ctetra10.nids = nodes[inonzero, :]
        self.increase_card_count('CTETRA', nentries)
        return n

    def _read_chexa(self, data, n):
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        ntotal = 88
        nentries = (len(data) - n) // ntotal
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        extra_nodes = nodes[:, 8:]
        assert extra_nodes.shape[1] == 12, extra_nodes.shape
        max_nid = extra_nodes.min(axis=1)
        izero = np.where(max_nid == 0)[0]
        inonzero = np.where(max_nid > 0)[0]
        if len(izero):
            self.chexa8.eid = eid[izero]
            self.chexa8.pid = pid[izero]
            self.chexa8.nids = nodes[izero, :8]
        if len(inonzero):
            self.chexa20.eid = eid[inonzero]
            self.chexa20.pid = pid[inonzero]
            self.chexa20.nids = nodes[inonzero, :]
        self.increase_card_count('CHEXA', nentries)
        return n

    def _read_cpenta(self, data, n):
        """
        CPENTA(4108,41,280)      - the marker for Record 63
        CPENPR(7509,75,9992)     - the marker for Record 64
        CPENT15F(16500,165,9999) - the marker for Record 65
        CPENT6FD(16000,160,9999) - the marker for Record 66
        """
        ntotal = 68
        nentries = (len(data) - n) // ntotal
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).frombuffer()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        extra_nodes = nodes[:, 6:]
        assert extra_nodes.shape[1] == 9, extra_nodes.shape
        max_nid = extra_nodes.min(axis=1)
        izero = np.where(max_nid == 0)[0]
        inonzero = np.where(max_nid > 0)[0]
        if len(izero):
            self.cpenta6.eid = eid[izero]
            self.cpenta6.pid = pid[izero]
            self.cpenta6.nids = nodes[izero, :6]
        if len(inonzero):
            self.cpenta15.eid = eid[inonzero]
            self.cpenta15.pid = pid[inonzero]
            self.cpenta15.nids = nodes[inonzero, :]
        self.increase_card_count('CPENTA', nentries)
        return n

    def _read_ctria3(self, data, n):
        """
        CTRIA3(5959,59,282)    - the marker for Record 94
        """
        element_type = 'CTRIA3'
        return self._read_ctria3_ctriar_helper(data, n, element_type)

    def _read_ctria3_ctriar_helper(self, data, n, element_type):
        """used by ``read_ctria3`` and ``read_ctriar``"""
        ntotal = 52
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:5]
        assert nodes.shape[1] == 3, nodes.shape
        # 0    1    2   3   4   5       6     7        8       9
        #(eid, pid, n1, n2, n3, theta, zoffs, blank1,  blank2, tflag,
        #            10   11  12   13
                     #t1, t2, t3) = out

        theta = floats[:, 5]
        zoffset = floats[:, 6]
        thickness_flag = ints[:, 9]
        thickness = floats[:, 10:]

        if element_type == 'CTRIA3':
            elem = self.ctria3
        elif element_type == 'CTRIAR':
            elem = self.ctriar
        else:
            raise NotImplementedError(element_type)
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        elem.theta = theta
        elem.zoffset = zoffset
        elem.thickness_flag = thickness_flag
        elem.thickness = thickness
        self.card_count[element_type] = nentries
        return n

    def _read_ctriar(self, data, n):
        """
        CTRIAR(9200,92,385)    - the marker for Record 99
        """
        element_type = 'CTRIAR'
        return self._read_ctria3_ctriar_helper(data, n, element_type)

    def _read_cquad4(self, data, n):
        """
        common method for CQUAD4, CQUADR
        CQUAD4(2958,51,177)    - the marker for Record 70
        CQUAD4(13900,139,9989) - the marker for Record 71
        256,480 /

        Word Name Type Description
        1 EID    I Element identification number
        2 PID    I Property identification number
        3 G(4)   I Grid point identification numbers of connection points
        7 THETA RS Material property orientation angle or coordinate system ID
        8 ZOFFS RS Offset from the surface of grid points reference plane
        9 UNDEF None
        10 TFLAG I Alternate thickness flag
        11 T(4) RS Membrane thickness of element at grid points
        """
        element_type = 'CQUAD4'
        return self._read_cquad4_cquadr_helper(data, n, element_type)

    def _read_cquadr(self, data, n):
        """see ``read_cquad4``"""
        element_type = 'CQUADR'
        return self._read_cquad4_cquadr_helper(data, n, element_type)

    def _read_cquad4_cquadr_helper(self, data, n, element_type):
        """used by ``read_cquad4`` and ``read_cquadr``"""
        ntotal = 56
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:6]
        # 0    1    2   3   4   5   6       7     8         9
        #(eid, pid, n1, n2, n3, n4, theta, zoffs, blank, tflag,
        #            10   11  12   13
                     #t1, t2, t3, t4) = out
        theta = floats[:, 7]
        zoffset = floats[:, 8]
        thickness_flag = ints[:, 9]
        thickness = floats[:, 10:]

        if element_type == 'CQUAD4':
            elem = self.cquad4
        elif element_type == 'CQUAD':
            elem = self.cquad
        elif element_type == 'CQUADR':
            elem = self.cquadr
        else:
            raise NotImplementedError(element_type)

        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        elem.theta = theta
        elem.zoffset = zoffset
        elem.thickness_flag = thickness_flag
        elem.thickness = thickness
        self.card_count[element_type] = nentries
        return n

    def _read_celas1(self, data, n):
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        ntotal = 24
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()

        #(eid, pid, g1, g2, c1, c2) = out
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:4]
        dofs = ints[:, 4:]
        elem = self.celas1
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        elem.dofs = dofs
        self.card_count['CELAS1'] = nentries
        return n

    def _read_celas2(self, data, n):
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        ntotal = 32
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #(eid, k, g1, g2, c1, c2, ge, s) = out
        eid = ints[:, 0]
        k = floats[:, 1]
        nodes = ints[:, 2:4]
        dofs = ints[:, 4:6]
        ge = floats[:, 6]
        s = floats[:, 7]
        elem = self.celas2
        elem.eid = eid
        elem.k = k
        elem.nids = nodes
        elem.dofs = dofs
        elem.ge = ge
        elem.s = s
        self.card_count['CELAS2'] = nentries
        return n


    def _read_celas3(self, data, n):
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #(eid, pid, s1, s2) = out
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        elem = self.celas3
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        self.card_count['CELAS3'] = nentries
        return n

    def _read_celas4(self, data, n):
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #(eid, pid, s1, s2) = out
        eid = ints[:, 0]
        k = floats[:, 1]
        nodes = ints[:, 2:]
        elem = self.celas3
        elem.eid = eid
        elem.k = k
        elem.nids = nodes
        self.card_count['CELAS4'] = nentries
        return n

    def _read_cdamp1(self, data, n):
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        ntotal = 24
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()

        #(eid, pid, g1, g2, c1, c2) = out
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:4]
        dofs = ints[:, 4:]
        elem = self.cdamp1
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        elem.dofs = dofs
        self.card_count['CDAMP1'] = nentries
        return n

    def _read_cdamp2(self, data, n):
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        ntotal = 24
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #(eid, bdamp, g1, g2, c1, c2) = out
        eid = ints[:, 0]
        bdamp = floats[:, 1]
        nodes = ints[:, 2:4]
        dofs = ints[:, 4:]
        elem = self.cdamp2
        elem.eid = eid
        elem.b = bdamp
        elem.nids = nodes
        self.card_count['CDAMP2'] = nentries
        return n

    def _read_cdamp3(self, data, n):
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()

        #(eid, pid, s1, s2) = out
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        elem = self.cdamp3
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        self.card_count['CDAMP3'] = nentries
        return n

    def _read_cdamp4(self, data, n):
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #(eid, pid, s1, s2) = out
        eid = ints[:, 0]
        bdamp = floats[:, 1]
        nodes = ints[:, 2:]
        elem = self.cdamp4
        elem.eid = eid
        elem.b = bdamp
        elem.nids = nodes
        self.card_count['CDAMP4'] = nentries
        return n

    def _read_cvisc(self, data, n):
        """CVISC(3901,39,50) - the marker for Record 105"""
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        #floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()

        #(eid, pid, n1, n2) = out
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        elem = self.cvisc
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        self.card_count['CVISC'] = nentries
        return n

    def _read_conrod(self, data, n):
        """
        CONROD(1601,16,47) - the marker for Record 58
        """
        ntotal = 32
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        floats = np.frombuffer(data[n:], self.fdtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        nodes = ints[:, 1:3]
        mid = ints[:, 3]
        a = floats[:, 4]
        j = floats[:, 5]
        c = floats[:, 6]
        nsm = floats[:, 7]
        elem = self.conrod
        elem.eid = eid
        elem.nids = nodes
        elem.mid = mid
        elem.A = a
        elem.j = j
        elem.c = c
        elem.nsm = nsm
        self.card_count['CONROD'] = nentries
        return n

        #ntotal = 32  # 8*4
        #s = Struct(b(self._endian + '4i4f'))
        #nelements = (len(data) - n) // ntotal
        #for i in range(nelements):
            #edata = data[n:n+32]
            #out = s.unpack(edata)
            #if self.is_debug_file:
                #self.binary_debug.write('  CONROD=%s\n' % str(out))
            #(eid, n1, n2, mid, a, j, c, nsm) = out
            #elem = CONROD.add_op2_data(out)
            #self.add_op2_element(elem)
            #n += ntotal
        #self.card_count['CONROD'] = nelements
        #return n

    def _read_crod(self, data, n):
        """
        CROD(3001,30,48)    - the marker for Record 81
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        elem = self.crod
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        self.card_count['CROD'] = nentries
        return n

    def _read_ctube(self, data, n):
        """
        CTUBE(3701,37,49) - the marker for Record 104
        """
        ntotal = 16
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        elem = self.ctube
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        self.card_count['CTUBE'] = nentries
        return n

    def _read_cbar(self, data, n):
        """
        CBAR(2408,24,180) - the marker for Record 8

        MSC/NX
        Word Name Type Description
        1 EID    I  Element identification number
        2 PID    I  Property identification number
        3 GA     I  Grid point identification number at end A
        4 GB     I  Grid point identification number at end B

        F=0* XYZ option -- basic coordinate system
           5 X1 RS  T1 component of orientation vector from GA
           6 X2 RS  T2 component of orientation vector from GA
           7 X3 RS  T3 component of orientation vector from GA
           8 FE  I  Orientation vector flag (encoded)
        F=1* XYZ option -- global coordinate system
           5 X1 RS  T1 component of orientation vector from GA
           6 X2 RS  T2 component of orientation vector from GA
           7 X3 RS  T3 component of orientation vector from GA
           8 FE  I   Orientation vector flag (encoded)
        F=2* Grid option
           5 GO I Grid point identification number at end of orientation vector
           6 UNDEF(2) none Not used
           8 FE I Orientation vector flag (encoded)
        *F = FE bit-wise AND with 3
        End F

        9  PA   I Pin flags for end A
        10 PB   I Pin flags for end B
        11 W1A RS T1 component of offset vector from GA
        12 W2A RS T2 component of offset vector from GA
        13 W3A RS T3 component of offset vector from GA
        14 W1B RS T1 component of offset vector from GB
        15 W2B RS T2 component of offset vector from GB
        16 W3B RS T3 component of offset vector from GB
        F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sebload1.op2
        """
        nelements = (len(data) - n) // 64
        for i in range(nelements):
            edata = data[n:n + 64]  # 16*4
            fe, = self.struct_i.unpack(edata[28:32])
            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            x = None
            g0 = None
            if f == 0:
                out = struct.unpack(self._endian + b'4i3f3i6f', edata)
                (eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
                x = [x1, x2, x3]
            elif f == 1:
                out = struct.unpack(self._endian + b'4i3f3i6f', edata)
                (eid, pid, ga, gb, x1, x2, x3, _f, pa, pb,
                 w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
                x = [x1, x2, x3]
            elif f == 2:
                out = struct.unpack(self._endian + b'7ii2i6f', edata)
                (eid, pid, ga, gb, g0, junk, junk, _f, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                data_in = [[eid, pid, ga, gb, pa, pb, w1a,
                            w2a, w3a, w1b, w2b, w3b], [f, g0]]
            else:
                raise RuntimeError(f'invalid f value...f={f}')

            nids = [ga, gb]
            elem = self.cbar.add(eid, pid, nids, x, g0, offt='GGG',
                                 pin_flags=[pa, pb],
                                 wa=[w1a, w2a, w3a],
                                 wb=[w1b, w2b, w3b], comment='')

            assert f == fe, 'f=%s type(f)=%s fe=%s\n%s' % (f, type(f), fe, elem)
            n += 64
        self.card_count['CBAR'] = nelements
        return n

    def _read_cbeam(self, data, n):
        """
        CBEAM(5408,54,261) - the marker for Record 10
        """
        nelements = (len(data) - n) // 72
        for i in range(nelements):
            edata = data[n:n + 72]  # 18*4
            fe, = self.struct_i.unpack(edata[40:44])

            # per DMAP: F = FE bit-wise AND with 3
            f = fe & 3
            x = None
            g0 = None
            if f == 0:  # basic cid
                out = struct.unpack(self._endian + b'6i3f3i6f', edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; basic cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
                x = [x1, x2, x3]
            elif f == 1:  # global cid
                out = struct.unpack(self._endian + b'6i3f3i6f', edata)
                (eid, pid, ga, gb, sa, sb, x1, x2, x3, fe, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; global cid' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, x1, x2, x3]]
                x = [x1, x2, x3]
            elif f == 2:  # grid option
                out = struct.unpack(self._endian + b'12i6f', edata)
                (eid, pid, ga, gb, sa, sb, g0, xx, xx, fe, pa,
                 pb, w1a, w2a, w3a, w1b, w2b, w3b) = out
                #self.log.info('CBEAM: eid=%s fe=%s f=%s; grid option' % (eid, fe, f))
                data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                           [f, g0]]
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            elem = self.cbeam.add(eid, pid, [ga, gb], x, g0,
                                  pin_flags=[pa, pb],
                                  wa=[w1a, w2a, w3a],
                                  wb=[w1b, w2b, w3b],
                                  sa=sa, sb=sb,
                                  comment='')
            #self.add_op2_element(elem)
            n += 72
        self.card_count['CBEAM'] = nelements
        return n

    def _read_cbush(self, data, n):
        """
        CBUSH(2608,26,60) - the marker for Record 13
        """
        nelements = (len(data) - n) // 56
        struct_obj1 = struct.Struct(self._endian + b'4i iii i ifi3f')
        struct_obj2 = struct.Struct(self._endian + b'4i fff i ifi3f')
        for i in range(nelements):
            edata = data[n:n + 56]  # 14*4
            out = struct_obj1.unpack(edata)
            eid, pid, ga, gb, five, sixi, seven, f, cid, s, ocid, s1, s2, s3 = out
            si = [s1, s2, s3]
            if f == -1: # Use Element CID below for orientation
                x = [None, None, None]
                g0 = None
            elif f in [0, 1]:
                # 5, 6, 7, f
                #0:4 4:8 8:12 12:16 16:20 20:24 24:28 28:32 32:36
                #0   1   2    3     4     5     6     7     8
                #x1, x2, x3, f2 = unpack('3f i', edata[20:36])
                out = struct_obj2.unpack(edata)
                eid, pid, ga, gb, x1, x2, x3, f2, cid, s, ocid, s1, s2, s3 = out

                assert f == f2, 'f=%s f2=%s' % (f, f2)
                x = [x1, x2, x3]
                g0 = None
            elif f == 2:
                x = [None, None, None]
                g0 = five
            else:
                raise RuntimeError('invalid f value...f=%r' % f)
            if cid == -1:
                cid = None
            data_in = [[eid, pid, ga, gb, cid, s, ocid, si], x, g0]

            #elem = CBUSH.add_op2_data(data_in, f)
            self.cbush.add(eid, pid, [ga, gb], x, g0, cid=cid, s=s, ocid=ocid,
                           si=si, comment='')
            n += 56
        self.card_count['CBUSH'] = nelements
        return n

    def _read_cshear(self, data, n):
        """
        CSHEAR(3101,31,61)    - the marker for Record 84
        """
        ntotal = 24
        nentries = (len(data) - n) // ntotal
        assert ntotal % 4 == 0, 'ntotal=%s nentries=%s ndata=%s n=%s' % (ntotal, nentries, len(data), n)
        nwords = ntotal // 4
        ints = np.frombuffer(data[n:], self.idtype).reshape(nentries, nwords).copy()

        #(eid, pid, n1, n2, n3, n4) = out
        eid = ints[:, 0]
        pid = ints[:, 1]
        nodes = ints[:, 2:]
        elem = self.cshear
        elem.eid = eid
        elem.pid = pid
        elem.nids = nodes
        self.card_count['CSHEAR'] = nentries
        return n

    def _add_load_object(self, load):
        pass
    def _add_constraint_spc_object(self, constraint):
        pass
    def _add_constraint_mpc_object(self, constraint):
        pass
    def _add_rigid_element_object(self, constraint):
        pass
    def _add_suport_object(self, constraint):
        pass
    def _add_uset_object(self, obj):
        pass
    def _add_thermal_load_object(self, load):
        pass
    def _add_nlparm_object(self, card, allow_overwrites=True):
        pass
    def _add_thermal_material_object(self, material, allow_overwrites=True):
        pass
    def _add_tstepnl_object(self, card, allow_overwrites=True):
        pass

    def _add_coord_object(self, coord, allow_overwrites: bool=False) -> None:
        """adds a Coord object"""
        key = coord.cid
        assert coord.cid > -1, 'cid=%s coord=\n%s' % (key, coord)
        if key in self.coords:
            #if not allow_overwrites:
            if not coord == self.coords[key]:
                self._duplicate_coords.append(coord)
        else:
            self.coords[key] = coord
            self._type_to_id_map[coord.type].append(key)

    #def _add_element_object(self, elem, allow_overwrites=False):
        #key = elem.eid
        #assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        #if key in self.elements and not allow_overwrites:
            #if not elem == self.elements[key]:
                #self._duplicate_elements.append(elem)
                #if self._stop_on_duplicate_error:
                    #self.pop_parse_errors()
        #else:
            #self.elements[key] = elem
            #self._type_to_id_map[elem.type].append(key)

    def _add_property_object(self, prop, allow_overwrites=False):
        """
        adds one of the following objects:
          PELAS, PBUSH, PBUSH1D, PBUSH2D, PDAMP,
          PROD, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
          PSHELL, PCOMP, PCOMPG,
          PSOLID, PLSOLID
        """
        key = prop.pid
        assert key > 0, 'pid=%s prop=%s' % (key, prop)
        if key in self.properties and not allow_overwrites:
            if not prop == self.properties[key]:
                self._duplicate_properties.append(prop)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.properties[key] = prop
            self._type_to_id_map[prop.type].append(key)

    def increase_card_count(self, name, count_num=1):
        self.card_count[name] = count_num
