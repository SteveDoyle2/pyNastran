"""
Defines:
  - Cart3D(log=None, debug=False)
     - read_cart3d(self, infilename, result_names=None)
     - write_cart3d(self, outfilename, is_binary=False, float_fmt='%6.7f')

     - flip_model()
     - make_mirror_model(self, nodes, elements, regions, loads, axis='y', tol=0.000001)
     - make_half_model(self, axis='y', remap_nodes=True)
     - get_free_edges(self, elements)
     - get_area(self)
     - get_normals(self)
     - get_normals_at_nodes(self, cnormals)

  - comp2tri(in_filenames, out_filename,
             is_binary=False, float_fmt='%6.7f')

"""
import sys
from struct import pack, unpack
from typing import TextIO, BinaryIO, Optional

import numpy as np
from cpylog import get_logger2, SimpleLogger
from pyNastran.utils import _filename

class Cart3dReaderWriter:
    """Cart3d IO class"""
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)
        self._endian = b''
        self._encoding = 'latin1'
        self.n = 0
        self.infilename = None
        self.points = np.zeros((0, 3), dtype='float64')
        self.elements = np.zeros((0, 3), dtype='int32')
        self.regions = np.zeros(0, dtype='int32')
        self.loads = {}

    @property
    def nresults(self) -> int:
        """get the number of results"""
        return len(self.loads)

    @property
    def nnodes(self) -> int:
        """alternate way to access number of points"""
        return self.npoints

    @property
    def npoints(self) -> int:
        """get the number of points"""
        return self.points.shape[0]

    @property
    def nodes(self) -> np.ndarray:
        """alternate way to access the points"""
        return self.points

    @nodes.setter
    def nodes(self, points) -> None:
        """alternate way to access the points"""
        self.points = points

    @property
    def nelements(self) -> int:
        """get the number of elements"""
        return self.elements.shape[0]

    def _read_elements_ascii(self, infile: TextIO, npoints: int, nelements: int) -> np.ndarray:
        """
        An element is defined by n1,n2,n3 and the ID is the location in elements.

        """
        assert nelements > 0, 'npoints=%s nelements=%s' % (npoints, nelements)
        elements = np.zeros((nelements, 3), dtype='int32')

        ieid = 0
        data = []
        while ieid < nelements:
            data += infile.readline().strip().split()
            while len(data) > 2:
                n1 = int(data.pop(0))
                n2 = int(data.pop(0))
                n3 = int(data.pop(0))
                elements[ieid] = [n1, n2, n3]
                ieid += 1

        nid_min = elements.min()
        if nid_min != 1:
            nid_max = elements.max()
            nnodes = self.nodes.shape[0]
            if nid_max == nnodes:
                msg = (
                    'Possible Cart3d error due to unused nodes\n'
                    'min(nids)=%s; expected 1; nid_max=%s nnodes=%s' % (
                        nid_min, nid_max, nnodes))
                self.log.warning(msg)
            else:
                msg = 'elements:\n%s\nmin(nids)=%s; expected 1; nid_max=%s nnodes=%s' % (
                    elements, nid_min, nid_max, nnodes, )
                #raise RuntimeError(msg)
                self.log.warning(msg)
            #assert elements.min() == 1, elements.min()
        return elements - 1

    def _read_cart3d_ascii(self, cart3d_filename: str, encoding: str,
                           result_names=None) -> None:
        log = self.log
        cart3d_filename = _filename(cart3d_filename)
        with open(cart3d_filename, 'r', encoding=self._encoding) as infile:
            try:
                npoints, nelements, nresults = _read_header_ascii(infile)
                self.points = _read_points_ascii(infile, npoints)
                self.elements = self._read_elements_ascii(infile, npoints, nelements)
                self.regions = _read_regions_ascii(infile, nelements)
                results, result_names = _read_results_ascii(0, infile, npoints, nresults, log,
                                                            result_names=result_names)
                if results is not None:
                    self.loads = _calculate_results(result_names, results, self.log)
            except Exception:
                msg = f'failed reading {cart3d_filename!r}'
                log.error(msg)
                raise

    def _read_cart3d_binary(self, cart3d_filename: str, endian: bytes) -> None:
        self.n = 0
        with open(cart3d_filename, 'rb') as infile:
            try:
                npoints, nelements, nresults, endian = self._read_header_binary(infile)
                self.points = self._read_points_binary(infile, npoints, endian)
                self.elements = self._read_elements_binary(infile, nelements, endian)
                self.regions = self._read_regions_binary(infile, nelements, endian)
                # TODO: loads
            except Exception:
                msg = f'failed reading {cart3d_filename!r}'
                self.log.error(msg)
                raise
            assert self.n == infile.tell(), 'n=%s tell=%s' % (self.n, infile.tell())

    def _read_header_binary(self, infile: BinaryIO) -> tuple[int, int, int, bytes]:
        """
        Reads the header::

          npoints nelements          # geometry
          npoints nelements nresults # results

        """
        log = self.log

        data = infile.read(4)
        size_little, = unpack(b'<i', data)
        size_big, = unpack(b'>i', data)

        if size_big in [12, 8]:
            endian = b'>'
            size = size_big
        elif size_little in [8, 12]:
            endian = b'<'
            size = size_little
        else:
            self._rewind(infile)
            self.show(infile, 100)
            raise RuntimeError('unknown endian')
        self._endian = endian
        self.n += 4

        data = infile.read(size)
        self.n += size

        so4 = size // 4  # size over 4
        if so4 == 3:
            (npoints, nelements, nresults) = unpack(endian + b'iii', data)
            log.info(f'npoints={npoints:d} nelements={nelements:d} nresults={nresults:d}')
        elif so4 == 2:
            (npoints, nelements) = unpack(self._endian + b'ii', data)
            nresults = 0
            log.info(f'npoints={npoints:d} nelements={nelements:d}')
        else:
            self._rewind(infile)
            self.show(infile, 100, endian=endian)
            raise RuntimeError(f'in the wrong spot...endian...size/4={so4}')
        infile.read(8)  # end of first block, start of second block
        self.n += 8
        return npoints, nelements, nresults, endian

    def _read_points_binary(self, infile: BinaryIO, npoints: int, endian: bytes) -> np.ndarray:
        """reads the xyz points"""
        size = npoints * 12  # 12=3*4 all the points
        data = infile.read(size)
        self.n += size

        dtype = np.dtype(endian + b'f4')
        points = np.frombuffer(data, dtype=dtype).reshape((npoints, 3)).copy()

        infile.read(8)  # end of second block, start of third block
        self.n += 8
        return points

    def _read_elements_binary(self, infile: BinaryIO, nelements: int, endian: bytes) -> np.ndarray:
        """reads the triangles"""
        size = nelements * 12  # 12=3*4 all the elements
        data = infile.read(size)
        self.n += size

        dtype = np.dtype(endian + b'i4')
        elements = np.frombuffer(data, dtype=dtype).reshape((nelements, 3)).copy()

        infile.read(8)  # end of third (element) block, start of regions (fourth) block
        self.n += 8
        assert elements.min() == 1, elements.min()
        return elements - 1

    def _read_regions_binary(self, infile: BinaryIO, nelements: int, endian: bytes) -> np.ndarray:
        """reads the regions"""
        size = nelements * 4  # 12=3*4 all the elements
        data = infile.read(size)
        self.n += size

        regions = np.zeros(nelements, dtype='int32')
        dtype = endian + b'i'
        regions = np.frombuffer(data, dtype=dtype).copy()

        infile.read(4)  # end of regions (fourth) block
        self.n += 4
        return regions

    def _read_results_binary(self, i: int, infile, result_names=None):
        """binary results are not supported"""
        pass

    def _rewind(self, infile):  # pragma: no cover
        """go back to the beginning of the file"""
        self.n = 0
        infile.seek(self.n)

    def show(self, infile: BinaryIO, n: int, types='ifs', endian=None):  # pragma: no cover
        assert self.n == infile.tell(), 'n=%s tell=%s' % (self.n, infile.tell())
        #nints = n // 4
        data = infile.read(4 * n)
        strings, ints, floats = self.show_data(data, types=types, endian=endian)
        infile.seek(self.n)
        return strings, ints, floats

    def show_data(self, data: bytes, types='ifs', endian=None):  # pragma: no cover
        if endian is None:
            endian = self._endian
        return show_data(sys.stdout, data, endian, types=types)

    def show_ndata(self, infile: BinaryIO, n: int, types: str='ifs'):  # pragma: no cover
        return self._write_ndata(infile, sys.stdout, n, types=types)

    def _write_ndata(self, infile: BinaryIO, outfile, n: int, types: str='ifs'):  # pragma: no cover
        """Useful function for seeing what's going on locally when debugging."""
        nold = self.n
        data = infile.read(n)
        self.n = nold
        infile.seek(self.n)
        return _write_data(outfile, data, types=types, endian=self._endian)

def _write_cart3d_ascii(outfilename: str,
                        points: np.ndarray, elements: np.ndarray, regions: np.ndarray,
                        loads, float_fmt: str) -> None:
    npoints = points.shape[0]
    nelements = elements.shape[0]
    is_loads = len(loads) > 0
    with open(outfilename, 'w') as outfile:
        int_fmt = _write_header_ascii(outfile, npoints, nelements, is_loads)
        _write_points_ascii(outfile, points, float_fmt)
        _write_elements_ascii(outfile, elements, int_fmt)
        _write_regions_ascii(outfile, regions)

        if is_loads:
            _write_loads_ascii(outfile, loads, npoints, float_fmt='%6.6f')

def _write_cart3d_binary(outfilename: str,
                         points: np.ndarray, elements: np.ndarray, regions: np.ndarray,
                         loads, endian):
    npoints = points.shape[0]
    nelements = elements.shape[0]

    is_loads = len(loads) > 0
    with open(outfilename, 'wb') as outfile:
        _write_header_binary(outfile, npoints, nelements, is_loads, endian)
        _write_points_binary(outfile, points, endian)
        _write_elements_binary(outfile, elements, endian)
        _write_regions_binary(outfile, regions, endian)

        if is_loads:
            raise NotImplementedError('loads writing is not supported; set model.loads={}')
            #self._write_loads(outfile, self.loads, is_binary, float_fmt)


def convert_to_float(svalues: list[str]) -> list[float]:
    """Takes a list of strings and converts them to floats."""
    values = []
    for value in svalues:
        values.append(float(value))
    return values

def _get_list(sline: list[str]) -> list[float]:
    """Takes a list of strings and converts them to floats."""
    try:
        sline2 = convert_to_float(sline)
    except ValueError:
        print("sline = %s" % sline)
        raise SyntaxError('cannot parse %s' % sline)
    return sline2

def b(mystr: str) -> bytes:
    """reimplementation of six.b(...) to work in Python 2"""
    return mystr.encode('ascii')

def show_data(outfile, data: bytes, endian: bytes, types: str='ifs'):  # pragma: no cover
    return _write_data(outfile, data, endian, types=types)

def _write_data(outfile: TextIO, data: bytes, endian: bytes, types: str='ifs'):  # pragma: no cover
    """Useful function for seeing what's going on locally when debugging."""
    n = len(data)
    nints = n // 4
    ndoubles = n // 8
    strings = None
    ints = None
    floats = None
    longs = None

    if endian is None:
        raise RuntimeError(endian)
        #endian = self._endian

    endian_str = endian.decode('ascii')
    if 's' in types:
        strings = unpack('%s%is' % (endian_str, n), data)
        outfile.write("strings = %s\n" % str(strings))
    if 'i' in types:
        ints = unpack('%s%ii' % (endian_str, nints), data)
        outfile.write("ints    = %s\n" % str(ints))
    if 'f' in types:
        floats = unpack('%s%if' % (endian_str, nints), data)
        outfile.write("floats  = %s\n" % str(floats))

    if 'l' in types:
        longs = unpack('%s%il' % (endian_str, nints), data)
        outfile.write("long  = %s\n" % str(longs))
    if 'I' in types:
        ints2 = unpack('%s%iI' % (endian_str, nints), data)
        outfile.write("unsigned int = %s\n" % str(ints2))
    if 'L' in types:
        longs2 = unpack('%s%iL' % (endian_str, nints), data)
        outfile.write("unsigned long = %s\n" % str(longs2))
    if 'q' in types:
        longs = unpack('%s%iq' % (endian_str, ndoubles), data[:ndoubles*8])
        outfile.write("long long = %s\n" % str(longs))
    return strings, ints, floats

def _read_header_ascii(infile: TextIO) -> tuple[int, int, int]:
    """
    Reads the header::

      npoints nelements          # geometry
      npoints nelements nresults # results

    """
    line = infile.readline()
    sline = line.strip().split()
    if len(sline) == 2:
        npoints, nelements = int(sline[0]), int(sline[1])
        nresults = 0
    elif len(sline) == 3:
        npoints = int(sline[0])
        nelements = int(sline[1])
        nresults = int(sline[2])
    else:
        raise ValueError('invalid result type')
    return npoints, nelements, nresults

def _read_points_ascii(infile: TextIO, npoints: int) -> np.ndarray:
    """A point is defined by x,y,z and the ID is the location in points."""
    p = 0
    data = []
    assert npoints > 0, 'npoints=%s' % npoints
    points = np.zeros((npoints, 3), dtype='float32')
    while p < npoints:
        data += infile.readline().strip().split()
        while len(data) > 2:
            x = data.pop(0)
            y = data.pop(0)
            z = data.pop(0)
            try:
                points[p] = [x, y, z]
            except:
                msg = f'failed reading ipoint={p}\nxyz=[{x},{y},{z}]'
                x2 = _bad_parse_float(x)
                y2 = _bad_parse_float(y)
                z2 = _bad_parse_float(z)
                raise RuntimeError(msg + f'; update=[{x2!r},{y2!r},{z2!r}]')
            p += 1
    return points

def _bad_parse_float(in_value: str) -> str:
    """hack to help xyz parsing"""
    try:
        out = float(in_value)
        out = in_value.strip()
    except ValueError:
        out = '???'
    return out

def _read_regions_ascii(infile: TextIO, nelements: int) -> np.ndarray:
    """reads the region section"""
    regions = np.zeros(nelements, dtype='int32')
    iregion = 0
    data = []
    while iregion < nelements:
        data = infile.readline().strip().split()
        ndata = len(data)
        regions[iregion : iregion + ndata] = data
        iregion += ndata
    return regions


def _read_results_ascii(i: int, infile: TextIO,
                        npoints: int, nresults: int, log: SimpleLogger,
                        result_names: Optional[list[str]]=None) -> tuple[Optional[np.ndarray],
                                                                         list[str]]:
    """
    Reads the Cp results.
    Results are read on a nodal basis from the following table:
      Cp
      rho,rhoU,rhoV,rhoW,rhoE

    With the following definitions:
      Cp = (p - 1/gamma) / (0.5*M_inf*M_inf)
      rhoVel^2 = rhoU^2+rhoV^2+rhoW^2
      M^2 = rhoVel^2/rho^2

    Thus:
      p = (gamma-1)*(e- (rhoU**2+rhoV**2+rhoW**2)/(2.*rho))
      p_dimensional = qInf * Cp + pInf

    # ???
    rho,rhoU,rhoV,rhoW,rhoE

    Parameters
    ----------
    result_names : list[str]; default=None (All)
        result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                        'Mach', 'U', 'V', 'W', 'E']

    """
    if nresults == 0:
        return None, []
    if result_names is None:
        result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                        'Mach', 'U', 'V', 'W', 'E', 'a', 'T', 'Pressure', 'q']
    log.debug('---starting read_results---')

    if nresults == 6:
        results = np.zeros((npoints, 6), dtype='float32')

        i0 = infile.tell()
        try:
            i2 = i
            all_slines = []
            for ipoint in range(npoints):
                sline = [infile.readline().strip()]
                sline += infile.readline().strip().split()
                i2 += 2
                all_slines.append(sline)
            results = np.array(all_slines, dtype='float64')
            i = i2
        except ValueError:
            infile.seek(i0)
            for ipoint in range(npoints):
                # Cp
                # rho       rhoU      rhoV      rhoW      E
                # 0.416594
                # 1.095611  0.435676  0.003920  0.011579  0.856058
                sline = [infile.readline().strip()]
                sline += infile.readline().strip().split()
                i += 2
                values = _get_list(sline)
                results[ipoint, :] = values
    else:
        raise RuntimeError('only nresults=6 is supported')
        #p=0
        #cp = sline[0]
        #rho = float(sline[1])
        #if(rho > abs(0.000001)):
            #rhoU = float(sline[2])
            #rhoV = float(sline[3])
            #rhoW = float(sline[4])
            #rhoE = float(sline[5])
            #mach2 = (rhoU) ** 2 + (rhoV) ** 2 + (rhoW) ** 2 / rho ** 2
            #mach = sqrt(mach2)
            #if mach > 10:
                #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (
                    #pointNum, cp, mach, rho, rhoU, rhoV, rhoW))
        #print("pt=%s i=%s Cp=%s p=%s" %(pointNum,i,sline[0],p))
    return results, result_names


def _calculate_results(result_names: list[str], results: np.ndarray,
                       log: SimpleLogger, loads=None):
    """
    Takes the Cart3d variables and calculates additional variables

    Parameters
    ----------
    result_names : list[str]
        the variables to calculate
    results : (n,6) ndarray
        the non-dimensional primitive flow variables
    loads : dict; default=None -> {}
        key : ???
        value : ???

    """
    if loads is None:
        loads = {}
    Cp = results[:, 0]
    rho = results[:, 1]
    rho_u = results[:, 2]
    rho_v = results[:, 3]
    rho_w = results[:, 4]
    E = results[:, 5]

    ibad = np.where(rho <= 0.000001)[0]
    if len(ibad) > 0:

        if 'Mach' in result_names:
            Mach = np.sqrt(rho_u**2 + rho_v**2 + rho_w**2)# / rho
            Mach[ibad] = 0.0
        if 'U' in result_names:
            U = rho_u / rho
            U[ibad] = 0.0
        if 'U' in result_names:
            V = rho_v / rho
            V[ibad] = 0.0
        if 'W' in result_names:
            W = rho_w / rho
            W[ibad] = 0.0
        #if 'rhoE' in result_names:
            #rho_e = rhoE / rho
            #e[ibad] = 0.0

        is_bad = True
        #n = 0
        #for i in ibad:
            #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (
                #i, Cp[i], Mach[i], rho[i], rho_u[i], rho_v[i], rho_w[i]))
            #Mach[i] = 0.0
            #n += 1
            #if n > 10:
            #    break
    else:
        is_bad = False


    #loc = locals()
    if 'Cp' in result_names:
        loads['Cp'] = Cp
    if 'rhoU' in result_names:
        loads['rhoU'] = rho_u
    if 'rhoV' in result_names:
        loads['rhoV'] = rho_v
    if 'rhoW' in result_names:
        loads['rhoW'] = rho_w
    #if 'rhoE' in result_names:
        #loads['rhoE'] = rho_e

    if 'rho' in result_names:
        loads['rho'] = rho

    is_mach = False
    if 'Mach' in result_names:
        if not is_bad:
            #Mach = np.sqrt(rho_u**2 + rho_v**2 + rho_w**2) / rho
            Mach = np.sqrt(rho_u**2 + rho_v**2 + rho_w**2)
        loads['Mach'] = Mach
        is_mach = True

    if 'U' in result_names:
        if not is_bad:
            U = rho_u / rho
        loads['U'] = U
    if 'V' in result_names:
        if not is_bad:
            V = rho_v / rho
        loads['V'] = V
    if 'W' in result_names:
        if not is_bad:
            W = rho_w / rho
        loads['W'] = W
    if 'E' in result_names:
        #if not is_bad:
            #E = rhoE / rho
        loads['E'] = E

    gamma = 1.4
    qinf = 1.0
    pinf = 1. / gamma
    Tinf = 1.0
    #Cp = (p - pinf) / qinf
    p = Cp * qinf + pinf

    #Cp_neg = -1.196898
    #rho_neg = 0.874096
    #p_neg = Cp_neg * qinf + pinf
    #T_neg = (Tinf * gamma * p_neg / rho_neg)

    T = (Tinf * gamma) * p / rho

    if 'a' in result_names:
        # T.min() = -0.773
        # is nondimensional qinf wrong?
        # -> more likely the formula for T is wrong
        #print('T: min=%s max=%s' % (T.min(), T.max()))
        loads['a'] = np.sqrt(T)
    if 'T' in result_names:
        loads['T'] = T

    if 'Pressure' in result_names:
        loads['Pressure'] = p
    if 'q' in result_names:
        q = 0.5 * rho * Mach ** 2
        loads['q'] = q
    # dynamic pressure
    # speed of sound
    # total pressure = p0/rhoi*ainf**2
    # total density
    # entropy
    # kinetic energy
    # enthalpy
    # energy, E
    # total energy
    # total enthalpy

    #i = where(Mach == max(Mach))[0][0]
    #log.info("i=%s Cp=%s rho=%s rho_u=%s rho_v=%s rho_w=%s Mach=%s" % (
        #i, Cp[i], rho[i], rho_u[i], rho_v[i], rho_w[i], Mach[i]))
    log.debug('---finished read_results---')
    return loads

def _write_header_binary(outfile: BinaryIO, npoints: int, nelements: int,
                         is_loads: bool, endian: bytes) -> None:
    """writes the cart3d header"""
    if is_loads:
        fmt = endian + b'iiiii'
        msg = pack(fmt, 3*4, npoints, nelements, 6, 4)
    else:
        fmt = endian + b'iiii'
        msg = pack(fmt, 2*4, npoints, nelements, 4)
    outfile.write(msg)

def _write_header_ascii(outfile: TextIO, npoints: int, nelements: int,
                        is_loads: bool) -> str:
    """
    Writes the header::

      npoints nelements          # geometry
      npoints nelements nresults # results

    """
    if is_loads:
        msg = '%i %i 6\n' % (npoints, nelements)
    else:
        msg = '%i %i\n' % (npoints, nelements)
    outfile.write(msg)

    # take the max value, string it, and length it
    # so 123,456 is length 6
    int_fmt = '%%%si' % len(str(nelements))
    return int_fmt

def _write_points_binary(outfile: BinaryIO, points: np.ndarray, endian: bytes) -> None:
    """writes the points"""
    four = pack(endian + b'i', 4)
    outfile.write(four)

    npoints = points.shape[0]
    fmt = endian + b'%if' % (npoints * 3)
    floats = pack(fmt, *np.ravel(points))

    outfile.write(floats)
    outfile.write(four)

def _write_points_ascii(outfile: TextIO, points: np.ndarray, float_fmt='%6.6f') -> None:
    """writes the points"""
    # if isinstance(float_fmt, bytes):
    #     fmt_ascii = float_fmt
    # else:
    #     fmt_ascii = float_fmt.encode('latin1')
    np.savetxt(outfile, points, fmt=float_fmt)

def _write_elements_binary(outfile: BinaryIO, elements: np.ndarray, endian: bytes) -> None:
    """writes the triangles"""
    fmt = endian + b'i'
    four = pack(fmt, 4)
    outfile.write(four)
    nelements = elements.shape[0]
    fmt = endian + b'%ii' % (nelements * 3)
    ints = pack(fmt, *np.ravel(elements+1))

    outfile.write(ints)
    outfile.write(four)

def _write_elements_ascii(outfile: BinaryIO, elements: np.ndarray, int_fmt: str) -> None:
    """writes the triangles"""
    #fmt_ascii = int_fmt.encode('latin1')
    np.savetxt(outfile, elements+1, fmt=int_fmt)

def _write_regions_binary(outfile: BinaryIO, regions: np.ndarray, endian: bytes) -> None:
    """writes the regions"""
    fmt = endian + b'i'
    four = pack(fmt, 4)
    outfile.write(four)

    nregions = len(regions)
    fmt = endian + b'%ii' % nregions
    ints = pack(fmt, *regions)
    outfile.write(ints)

    outfile.write(four)

def _write_regions_ascii(outfile: TextIO, regions: np.ndarray) -> None:
    """writes the regions"""
    fmt = '%i'
    np.savetxt(outfile, regions, fmt=fmt)

def _write_loads_ascii(outfile, loads, npoints: int, float_fmt='%6.6f'):
    """writes the *.triq loads

    Results are read on a nodal basis from the following table:
      Cp
      rho,rhoU,rhoV,rhoW,rhoE

    """
    Cp = loads['Cp']
    rho = loads['rho']
    rhoU = loads['rhoU']
    rhoV = loads['rhoV']
    rhoW = loads['rhoW']
    E = loads['E']
    assert len(Cp) == npoints, 'len(Cp)=%s npoints=%s' % (len(Cp), npoints)

    fmt = '%s\n%s %s %s %s %s\n' % (float_fmt, float_fmt, float_fmt,
                                    float_fmt, float_fmt, float_fmt)
    for (cpi, rhoi, rhou, rhov, rhoe, e) in zip(Cp, rho, rhoU, rhoV, rhoW, E):
        outfile.write(fmt % (cpi, rhoi, rhou, rhov, rhoe, e))
