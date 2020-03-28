# pylint: disable=R0902,R0904,R0914
from math import sin, cos, radians, atan2, sqrt, degrees
from itertools import count
from typing import Tuple # , TYPE_CHECKING

import numpy as np
from numpy import array, zeros
from scipy.sparse import coo_matrix  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string, string_or_blank,
    parse_components, interpret_value, integer_double_string_or_blank)


class DTI(BaseCard):
    """
    +-----+-------+-----+------+-------+--------+------+-------------+
    |  1  |   2   |  3  |   4  |   5   |    6   |   7  |       8     |
    +=====+=======+=====+======+=======+========+======+=============+
    | DTI | UNITS | "1" | MASS | FORCE | LENGTH | TIME |   STRESS    |
    +-----+-------+-----+------+-------+--------+------+-------------+

    MSC

    +-----+-------+-----+------+-------+--------+------+-------------+
    |  1  |   2   |  3  |   4  |   5   |    6   |   7  |       8     |
    +=====+=======+=====+======+=======+========+======+=============+
    | DTI | UNITS | "1" | MASS | FORCE | LENGTH | TIME | TEMPERATURE |
    +-----+-------+-----+------+-------+--------+------+-------------+

    NX
    """
    type = 'DTI'
    #_properties = ['shape', 'ifo', 'is_real', 'is_complex', 'is_polar', 'matrix_type', 'tin_dtype', 'tout_dtype']

    @classmethod
    def _init_from_empty(cls):
        name = 'name'
        fields = []
        return DTI(name, fields, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        keys, values = self.fields

        # nan != nan
        values = [value if value == value else None for value in values]
        values_str = [value.decode(encoding) if isinstance(value, bytes) else value
                      for value in values]
        #values = [valuei.decode(encoding) if isinstance(valuei, bytes) else (
        #    None if np.isnan(valuei) else valuei)
        #          for valuei in values]
        self.fields = {key : value for key, value in zip(keys, values_str)}

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        """exports the elements in a vectorized way"""
        from pyNastran.bdf.bdf_interface.hdf5_exporter import _export_list
        for name, dti in sorted(model.dti.items()):
            if name == 'UNITS':
                i = 0
                for key, value in sorted(dti.fields.items()):
                    #print(key, value)
                    h5_group = h5_file.create_group(str(key))
                    if value is None:
                        h5_group.create_dataset(str(i), data=np.nan)
                    else:
                        h5_group.create_dataset(str(i), data=value)
                    i += 1
                #fields = {
                    #'mass' : mass,
                    #'force' : force,
                    #'length' : length,
                    #'time' : time,
                    #'temp_stress' : temp_stress
                #}
            else:
                for irecord, fields in sorted(dti.fields.items()):
                    #h5_group = h5_file.create_group(str(irecord))
                    attr = 'irecord=%s' % irecord
                    namei = str(irecord)
                    values = fields
                    _export_list(h5_file, attr, namei, values, encoding)
                    #print(h5_group)
                    #print(irecord, fields)

    def __init__(self, name, fields, comment=''):
        """
        Creates a DTI card

        Parameters
        ----------
        name : str
            UNITS
        fields : List[varies]
            the fields
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.name = name
        self.fields = fields
        assert len(fields) > 0, fields

    @classmethod
    def add_card(cls, card, comment):
        """
        Adds a DTI card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        if name == 'UNITS':
            integer(card, 2, '1')
            mass = string_or_blank(card, 3, 'mass')
            force = string_or_blank(card, 4, 'force')
            length = string_or_blank(card, 5, 'length')
            time = string_or_blank(card, 6, 'time')
            temp_stress = string_or_blank(card, 7, 'stress/temperature')
            fields = {
                'mass' : mass,
                'force' : force,
                'length' : length,
                'time' : time,
                'temp_stress' : temp_stress
            }
        else:
            fields = []
            #field2 = card[2]

            list_fields = []
            irecord = integer(card, 2, 'record')
            if irecord == 0:
                for i in range(3, len(card)):
                    val = integer_double_string_or_blank(
                        card, i, 'T%i' % (i-1), default=32767)
                    list_fields.append(val)
            else:
                for i in range(3, len(card)):
                    val = integer_double_string_or_blank(
                        card, i, 'T%i' % (i-1), default=None)
                    list_fields.append(val)
            fields = {irecord: list_fields,}
        return DTI(name, fields, comment=comment)

    def raw_fields(self):
        if self.name == 'UNITS':
            mass = self.fields['mass']
            force = self.fields['force']
            length = self.fields['length']
            time = self.fields['time']
            temp_stress = self.fields['temp_stress']
            list_fields = ['DTI', self.name, '1', mass, force, length, time, temp_stress]
        else:
            list_fields = []
            for irecord, fields in sorted(self.fields.items()):
                nfields = len(fields)
                list_fields += ['DTI', self.name] + fields
                nleftover = nfields % 8
                if nleftover:
                    list_fields += [None] * nleftover
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        if self.name == 'UNITS':
            card = self.repr_fields()
            return self.comment + print_card_8(card)

        msg = self.comment
        for irecord, fields in sorted(self.fields.items()):
            list_fields = ['DTI', self.name, irecord, ] + fields
            msg += print_card_8(list_fields)
        return msg


class NastranMatrix(BaseCard):
    """
    Base class for the DMIG, DMIJ, DMIJI, DMIK matrices
    """
    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        self.finalize()

    def __init__(self, name, matrix_form, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment='', finalize=True):
        """
        Creates a NastranMatrix

        Parameters
        ----------
        name : str
            the name of the matrix
        matrix_form : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : List[(node, dof)]
            the jnode, jDOFs
        GCi  : List[(node, dof)]
            the inode, iDOFs
        Real : List[float]
            The real values
        Complex : List[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        if Complex is None:
            Complex = []
        if tout is None:
            tout = 0

        polar = _set_polar(polar)

        if matrix_form not in [1, 2, 4, 5, 6, 8, 9]:
            msg = (
                'matrix_form=%r must be [1, 2, 4, 5, 6, 8, 9]\n'
                '  1: Square\n'
                '  2: Rectangular\n'
                #'  4: Lower Triangular\n'
                #'  5: Upper Triangular\n'
                '  6: Symmetric\n'
                #'  8: Identity (m=nRows, n=m)\n'
                '  9: Rectangular\n' % matrix_form)
            raise ValueError(msg)
        self.name = name

        #: 4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.matrix_form = matrix_form

        #: 1-Real, Single Precision; 2=Real,Double Precision;
        #  3=Complex, Single; 4=Complex, Double
        self.tin = tin

        #: 0-Set by cell precision
        self.tout = tout

        #: Input format of Ai, Bi. (Integer=blank or 0 indicates real, imaginary format;
        #: Integer > 0 indicates amplitude, phase format.)
        self.polar = polar

        self.ncols = ncols
        self.GCj = GCj
        self.GCi = GCi

        self.Real = Real
        if len(Complex) or self.is_complex:
            self.Complex = Complex
            assert self.tin in [3, 4], 'tin=%r and must 3 or 4 to be complex' % self.tin
            assert self.tout in [0, 3, 4], 'tin=%r and must 0, 3 or 4 to be complex' % self.tout
        assert isinstance(matrix_form, integer_types), 'matrix_form=%r type=%s' % (matrix_form, type(matrix_form))
        assert not isinstance(matrix_form, bool), 'matrix_form=%r type=%s' % (matrix_form, type(matrix_form))
        if finalize:
            self.finalize()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a NastranMatrix (DMIG, DMIJ, DMIK, DMIJI) card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        #zero

        matrix_form = integer(card, 3, 'ifo')
        tin = integer(card, 4, 'tin')
        tout = integer_or_blank(card, 5, 'tout', 0)
        polar = integer_or_blank(card, 6, 'polar', 0)
        if matrix_form == 1: # square
            ncols = integer_or_blank(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        elif matrix_form == 6: # symmetric
            ncols = integer_or_blank(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        elif matrix_form in [2, 9]: # rectangular
            ncols = integer(card, 8, 'matrix_form=%s; ncol' % (matrix_form))
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'matrix_form=%s; ncol' % self.matrix_form)

            msg = (
                '%s name=%r matrix_form=%r is not supported.  Valid forms:\n'
                '  4=Lower Triangular\n'
                '  5=Upper Triangular\n'
                '  6=Symmetric\n'
                '  8=Identity (m=nRows, n=m)\n' % (cls.type, name, matrix_form)
            )
            raise NotImplementedError(msg)

        GCj = []
        GCi = []
        Real = []
        Complex = []
        return cls(name, matrix_form, tin, tout, polar, ncols,
                   GCj, GCi, Real, Complex, comment=comment, finalize=False)

    @property
    def matrix_type(self):
        """gets the matrix type"""
        if not isinstance(self.matrix_form, integer_types):
            msg = 'ifo must be an integer; matrix_form=%r type=%s name=%s' % (
                self.matrix_form, type(self.matrix_form), self.name)
            raise TypeError(msg)
        if isinstance(self.matrix_form, bool):
            msg = 'matrix_form must not be a boolean; matrix_form=%r type=%s name=%s' % (
                self.matrix_form, type(self.matrix_form), self.name)
            raise TypeError(msg)

        if self.matrix_form == 1:
            matrix_type = 'square'
        elif self.matrix_form == 6:
            matrix_type = 'symmetric'
        elif self.matrix_form in [2, 9]:
            matrix_type = 'rectangular'
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'matrix_form=%s; ncol' % self.matrix_form)
            raise NotImplementedError('%s matrix_form=%r is not supported' % (
                self.type, self.matrix_form))
        return matrix_type

    def finalize(self):
        """converts the lists into numpy arrays"""
        self.GCi = np.asarray(self.GCi)
        self.GCj = np.asarray(self.GCj)
        self.Real = np.asarray(self.Real)
        if self.is_complex:
            self.Complex = np.asarray(self.Complex)

    @property
    def shape(self):
        """gets the matrix shape"""
        if self.matrix_form in [1, 6]: # square, symmetric
            if self.ncols is not None:
                shape = (self.ncols, self.ncols)
            else:
                nrows, ncols = get_row_col_map(
                    self, self.GCi, self.GCj, self.matrix_form)[:2]
                shape = (nrows, ncols)
        elif self.matrix_form in [2, 9]:
            raise NotImplementedError('need to pull the nrows after reading in everything')
            #shape = (self.ncols, self.ncols)
        else:
            raise NotImplementedError('matrix_form=%s' % self.matrix_form)
        return shape

    def _add_column(self, card, comment=''):
        """adds an additional column entry to the matrix"""
        if comment:
            if hasattr(self, '_comment'):
                self.comment += comment
            else:
                self.comment = comment

        name = string(card, 1, 'name')
        if name == 'UACCEL':
            self._add_column_uaccel()
            return

        Gj = integer(card, 2, 'Gj')
        # Cj = integer(card, 3, 'Cj')
        Cj = integer_or_blank(card, 3, 'Cj', 0)
        #Cj = parse_components(card, 3, 'Cj')
        assert 0 <= Cj <= 6, 'C%i must be between [0, 6]; Cj=%s' % (0, Cj)

        nfields = len(card)
        #print("nfields = %i" % nfields)
        #print("card[5:] =", card[5:])
        #print("(nfields - 5) %% 4 = %i" % ((nfields - 5) % 4))

        nloops = (nfields - 5) // 4
        if (nfields - 5) % 4 in [2, 3]:  # real/complex
            nloops += 1
        #assert nfields <= 8,'nfields=%s' % nfields
        #print("nloops = %i" % nloops)
        assert nloops > 0, 'nloops=%s' % nloops

        for i in range(nloops):
            self.GCj.append((Gj, Cj))

        if self.is_complex:
            if self.is_polar:
                for i in range(nloops):
                    n = 5 + 4 * i
                    Gi = integer(card, n, 'Gi')
                    # Ci = integer(card, n + 1, 'Ci')
                    Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                    #Ci = parse_components(card, n + 1, 'Ci')
                    assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                    self.GCi.append((Gi, Ci))
                    magi = double(card, n + 2, 'ai')
                    phasei = double(card, n + 3, 'bi')
                    reali = magi * cos(radians(phasei))
                    complexi = magi * sin(radians(phasei))
                    self.Real.append(reali)
                    self.Complex.append(complexi)
            else:
                for i in range(nloops):
                    n = 5 + 4 * i
                    Gi = integer(card, n, 'Gi')
                    # Ci = integer(card, n + 1, 'Ci')
                    Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                    #Ci = parse_components(card, n + 1, 'Ci')
                    assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                    self.GCi.append((Gi, Ci))
                    reali = double(card, n + 2, 'real')
                    complexi = double(card, n + 3, 'complex')
                    self.Real.append(reali)
                    self.Complex.append(complexi)
        else:
            # real
            for i in range(nloops):
                n = 5 + 4 * i
                Gi = integer(card, n, 'Gi')
                # Ci = integer(card, n + 1, 'Ci')
                Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                #Ci = parse_components(card, n + 1, 'Ci')
                assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                reali = double(card, n + 2, 'real')
                self.GCi.append((Gi, Ci))
                self.Real.append(reali)
                #print("GC=%s,%s real=%s" % (Gi, Ci, reali))

        msg = '(len(GCj)=%s len(GCi)=%s' % (len(self.GCj), len(self.GCi))
        assert len(self.GCj) == len(self.GCi), msg
        #if self.is_complex:
            #self.Complex(double(card, v, 'complex')

    def get_matrix(self, is_sparse=False, apply_symmetry=True):
        """
        Builds the Matrix

        Parameters
        ----------
        is_sparse : bool; default=False
            should the matrix be returned as a sparse matrix.
            Slower for dense matrices.
        apply_symmetry : bool; default=True
            If the matrix is symmetric (ifo=6), returns a symmetric matrix.
            Supported as there are symmetric matrix routines.

        Returns
        -------
        M : numpy.ndarray or scipy.coomatrix
            the matrix
        rows : dict[int] = [int, int]
            dictionary of keys=rowID, values=(Grid,Component) for the matrix
        cols: dict[int] = [int, int]
            dictionary of keys=columnID, values=(Grid,Component) for the matrix

        .. warning:: is_sparse=True WILL fail

        """
        return get_matrix(self, is_sparse=is_sparse, apply_symmetry=apply_symmetry)

    @property
    def is_real(self):
        """real vs. complex attribute"""
        return not self.is_complex

    @property
    def is_complex(self):
        """real vs. complex attribute"""
        if self.tin in [1, 2]: # real
            return False
        elif self.tin in [3, 4]: # complex
            return True
        msg = ('Matrix %r must have a value of TIN = [1, 2, 3, 4].\n'
               'TIN defines the type (real, complex) '
               'of the matrix.  TIN=%r.\n'
               '  TIN=1,2 -> real\n'
               '  TIN=3,4 -> complex' % (self.name, self.tin))
        raise ValueError(msg)

    @property
    def is_polar(self):
        """
        Used by:
          - DMIG
          - DMIJ
          - DMIJI
          - DMIK

        Not used by:
          - DMI
          - DMIAX
          - DMIG, UACCEL
          - DMIGOUT
          - DMIGROT

        """
        if self.polar == 0: # real, imag
            return False
        elif self.polar == 1: # mag, phase
            return True
        elif self.polar is None:
            return False
        msg = ('Matrix %r must have a value of POLAR = [0, 1].\n'
               'POLAR defines the type (real/imag or mag/phase) complex) '
               'of the matrix.  POLAR=%r.' % (self.name, self.polar))
        raise ValueError(msg)

    @property
    def tin_dtype(self):
        """gets the input dtype"""
        return _get_dtype(self.is_complex, self.tin)

    @property
    def tout_dtype(self):
        """gets the output dtype"""
        return _get_dtype(self.is_complex, self.tout)

    def __repr__(self):
        return self.write_card(size=8, is_double=False)

    def fill_in_default_components(self, model):
        for i, (Gi, Ci) in enumerate(self.GCi):
            if Ci is None:
                node = model.nodes[Gi]
                if node.type == 'GRID':
                    msg = ('Ci on DMIG card must be 1, 2, 3, 4, 5, or 6; '
                           'Node=%i (GRID); Ci=%s' % (Gi, Ci))
                    raise RuntimeError(msg)
                elif node.type in ['SPOINT', 'EPOINT']:
                    Ci = 0
                else:
                    raise NotImplementedError(node)
                self.GCi[i] = [Gi, Ci]

        for i, (Gj, Cj) in enumerate(self.GCj):
            if Cj is None:
                node = model.nodes[Gj]
                if node.type == 'GRID':
                    msg = ('Cj on DMIG card must be 1, 2, 3, 4, 5, or 6; '
                           'Node=%i (GRID); Cj=%s' % (Gj, Cj))
                    raise RuntimeError(msg)
                elif node.type in ['SPOINT', 'EPOINT']:
                    Cj = 0
                else:
                    raise NotImplementedError(node)
                self.GCj[i] = [Gj, Cj]
        return

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        if self.tin in [1, 3]:
            is_double = False
        elif self.tin in [2, 4]:
            is_double = True
            size = 16
        else:
            raise RuntimeError('tin=%r must be 1, 2, 3, or 4' % self.tin)

        assert isinstance(self.GCi, (list, np.ndarray)), 'type(GCi)=%s' % type(self.GCi)
        assert isinstance(self.GCj, (list, np.ndarray)), 'type(GCj)=%s' % type(self.GCj)
        assert isinstance(self.Real, (list, np.ndarray)), 'type(Real)=%s' % type(self.Real)
        #assert isinstance(self.GCi[0], (list, np.ndarray)), 'type(GCi[0])=%s' % type(self.GCi[0])
        #assert isinstance(self.GCj[0], (list, np.ndarray)), 'type(GCj[0])=%s' % type(self.GCj[0])

        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % (self.type, self.name)
        list_fields = [self.type, self.name, 0, self.matrix_form, self.tin,
                       self.tout, self.polar, None, self.ncols]
        if size == 8:
            msg += print_card_8(list_fields)
        else:
            msg += print_card_16(list_fields)

        if self.is_complex:
            if self.is_polar:
                for (GCi, GCj, reali, complexi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                    magi = sqrt(reali**2 + complexi**2)
                    if reali == 0.0:
                        phasei = 0.0
                    else:
                        phasei = degrees(atan2(complexi, reali))
                    list_fields = [self.type, self.name, GCj[0], GCj[1],
                                   None, GCi[0], GCi[1], magi, phasei]
                    if size == 8:
                        msg += print_card_8(list_fields)
                    elif is_double:
                        msg += print_card_double(list_fields)
                    else:
                        msg += print_card_16(list_fields)
            else:
                for (GCi, GCj, reali, complexi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                    list_fields = [self.type, self.name, GCj[0], GCj[1],
                                   None, GCi[0], GCi[1], reali, complexi]
                    if size == 8:
                        msg += print_card_8(list_fields)
                    elif is_double:
                        msg += print_card_double(list_fields)
                    else:
                        msg += print_card_16(list_fields)
        else:
            for (GCi, GCj, reali) in zip(self.GCi, self.GCj, self.Real):
                list_fields = [self.type, self.name, GCj[0], GCj[1],
                               None, GCi[0], GCi[1], reali, None]
                if size == 8:
                    msg += print_card_8(list_fields)
                elif is_double:
                    msg += print_card_double(list_fields)
                else:
                    msg += print_card_16(list_fields)

        #msg += '\n\nGCi[0]=%s\n' % self.GCi[0]
        #msg += 'GCj[0]=%s\n' % self.GCj[0]
        #msg += 'Real[0]=%s\n' % self.Real[0]
        #assert isinstance(self.GCi[0], (list, np.ndarray)), msg
        #assert isinstance(self.GCj[0], (list, np.ndarray)), msg
        #assert isinstance(self.Real[0], (list, np.ndarray)), msg

        return msg


class DMIG_UACCEL(BaseCard):
    """
    Direct Matrix Input of Enforced Static Acceleration
    Defines rigid body accelerations in the basic coordinate system.

    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    |   1  |   2    |  3  |  4  |  5  |  6  |  7  |   8   |       |
    +======+========+=====+=====+=====+=====+=====+=======+=======+
    | DMIG | UACCEL | "0" | "9" | TIN |     |     |       | NCOL  |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  L  |     |     |  G1 | C1  |  X1   |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    |      |   G2   |  C2 | X2  |     |  G3 | C3  |  X3   |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+

    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  0  |  9  |  1  |     |     |       |   4   |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  2  |     |     |  2  |  3  | 386.4 |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  3  |     |     |  2  |  4  |  3.0  |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  4  |     |     |  2  |  6  |  1.0  |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    """
    type = 'DMIG'
    name = 'UACCEL'
    def __init__(self, tin, ncol, load_sequences, comment=''):
        if comment:
            self.comment = comment
        self.tin = tin
        self.ncol = ncol
        self.load_sequences = load_sequences
        #print(str(self))

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmig_to_hdf5(h5_file, model, model.dmig, encoding)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DMIG,UACCEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        tin = integer(card, 4, 'tin')
        ncol = integer_or_blank(card, 8, 'ncol')
        return DMIG_UACCEL(tin, ncol, load_sequences={}, comment=comment)

    def _add_column(self, card, comment=''):
        if comment:
            if hasattr(self, '_comment'):
                self.comment += comment
            else:
                self.comment = comment
        load_seq = integer(card, 2, 'load_seq')

        i = 0
        ifield = 5
        self.load_sequences[load_seq] = []
        assert len(card) >= 8, 'len=%s card=%s' % (len(card), card)
        while ifield < len(card):
            g1 = integer(card, ifield, 'nid%d' % i)
            c1 = parse_components(card, ifield+1, 'c%d' % i)
            x1 = double(card, ifield+2, 'x%d' % i)
            #assert len(card) <= 8, 'len=%s card=%s' % (len(card), card)
            gcx = [g1, c1, x1]
            self.load_sequences[load_seq].append(gcx)
            ifield += 4
            i += 1


    @staticmethod
    def finalize():
        """a passer method"""
        pass

    def raw_fields(self):
        list_fields = [
            'DMIG', 'UACCEL', 0, 9, self.tin, None, None, None, self.ncol
        ]
        for lseq, ncx in sorted(self.load_sequences.items()):
            list_fields += [lseq, None, None]
            for ncxi in ncx:
                list_fields += ncxi
           #for (nid, comp, xi) in ncx:
        #print('list_fields= %s' % list_fields)
        self.write_card()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        if self.tin in [1, 3]:
            is_double = False
            msg = self.write_card_8()
        elif self.tin in [2, 4]:
            is_double = True
            size = 16
            msg = self.write_card_16()
        else:
            raise RuntimeError('tin=%r must be 1, 2, 3, or 4' % self.tin)
        return msg

    def write_card_8(self):
        """writes the card in small field format"""
        return self._write_card(print_card_8)

    def write_card_16(self):
        """writes the card in small large format"""
        return self._write_card(print_card_16)

    def _write_card(self, func):
        """writes the card"""
        msg = '\n$' + '-' * 80
        msg += '\n$ DMIG Matrix UACCEL\n'
        list_fields = [
            'DMIG', 'UACCEL', 0, 9, self.tin, None, None, None, self.ncol,
        ]
        msg += func(list_fields)

        for lseq, ncx in sorted(self.load_sequences.items()):
            list_fields = ['DMIG', 'UACCEL']
            list_fields += [lseq, None, None]
            for ncxi in ncx:
                list_fields += ncxi + [None]
            list_fields.pop()
            msg += func(list_fields)
        #print(msg)
        #if self.is_complex:
            #msg += self._get_complex_fields(func)
        #else:
            #msg += self._get_real_fields(func)
        return msg

    def __repr__(self):
        return self.write_card(size=8)

class DMIG(NastranMatrix):
    """
    Defines direct input matrices related to grid, extra, and/or scalar points.
    The matrix is defined by a single header entry and one or more column
    entries. A column entry is required for each column with nonzero elements.

    +------+------+----+-----+-----+------+-------+----+------+
    |   1  |  2   | 3  |  4  |  5  |   6  |   7   | 8  |  9   |
    +======+======+====+=====+=====+======+=======+====+======+
    | DMIG | NAME | 0  | IFO | TIN | TOUT | POLAR |    | NCOL |
    +------+------+----+-----+-----+------+-------+----+------+
    | DMIG | NAME | GJ | CJ  |     |  G1  |  C1   | A1 |  B1  |
    +------+------+----+-----+-----+------+-------+----+------+
    |      |  G2  | C2 | A2  |  B2 |      |       |    |      |
    +------+------+----+-----+-----+------+-------+----+------+
    """
    type = 'DMIG'
    _properties = ['is_real', 'is_complex', 'is_polar', 'matrix_type', 'shape',
                   'tin_dtype', 'tout_dtype']

    #@classmethod
    #def _init_from_empty(cls):
        #name = 'name'
        #ifo = 1
        #tin = 1
        #tout = 1
        #polar = 0
        #ncols = 1
        #GCj = []
        #GCi = []
        #Real = []
        #return DMIG(name, ifo, tin, tout, polar, ncols, GCj, GCi, Real,
                    #Complex=None, comment='', finalize=True)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmig_to_hdf5(h5_file, model, model.dmig, encoding)

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment='', finalize=True):
        """
        Creates a DMIG card

        Parameters
        ----------
        name : str
            the name of the matrix
        ifo : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : List[(node, dof)]
            the [jnode, jDOFs]
        GCi  : List[(node, dof)]
            the inode, iDOFs
        Real : List[float]
            The real values
        Complex : List[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment,
                               finalize=finalize)


class DMIAX(BaseCard):
    """
    Direct Matrix Input for Axisymmetric Analysis

    Defines axisymmetric (fluid or structure) related direct input matrix
    terms.  The matrix is defined by a single header entry and one or
    more column entries. Only one header entry is required. A column
    entry is required for each column with nonzero elements.

    +-------+------+----+--------+------+--------+-------+----+------+
    |   1   |  2   | 3  |    4   |   5  |   6    |   7   | 8  |  9   |
    +=======+======+====+========+======+========+=======+====+======+
    | DMIAX | NAME | 0  |   IFO  | TIN  |  TOUT  |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+

    +-------+------+----+--------+------+--------+-------+----+------+
    |   1   |  2   | 3  |    4   |   5  |   6    |   7   | 8  |  9   |
    +=======+======+====+========+======+========+=======+====+======+
    | DMIAX | NAME | GJ |   CJ   |  NJ  |        |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+
    |       |  G1  | C1 |   N1   |  A1  |   B1   |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+
    |       |  G2  | C2 |  etc.  |      |        |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+

    +-------+------+----+--------+------+--------+-------+----+------+
    |   1   |  2   | 3  |    4   |   5  |   6    |   7   | 8  |  9   |
    +=======+======+====+========+======+========+=======+====+======+
    | DMIAX | B2PP | 0  |   1    |  3   |        |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+
    | DMIAX | B2PP | 32 |        |      |        |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+
    |       | 1027 | 3  | 4.25+6 |      | 2.27+3 |       |    |      |
    +-------+------+----+--------+------+--------+-------+----+------+

    """
    type = 'DMIAX'

    def __init__(self, name, matrix_form, tin, tout, ncols,
                 GCNj, GCNi, Real, Complex=None, comment=''):
        """
        Creates a DMIAX card

        Parameters
        ----------
        name : str
            the name of the matrix
        matrix_form : int
            matrix shape
            1=Square
            2=General Rectangular
            6=Symmetric
        tin : int
            matrix input precision
            1=Real, Single Precision
            3=Complex, Single Precision
        tout : int
            matrix output precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        GCNj  : List[(node, dof, harmonic_number)]???
            the jnode, jDOFs
        GCNi  : List[(node, dof, harmonic_number)]???
            the inode, iDOFs
        Real : List[float]???
            The real values
        Complex : List[float]???; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        ncols = None

        if comment:
            self.comment = comment

        if Complex is None:
            Complex = []

        if tout is None:
            tout = 0

        self.name = name

        #: ifo/4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.matrix_form = matrix_form

        #: 1-Real, Single Precision; 2=Real,Double Precision;
        #  3=Complex, Single; 4=Complex, Double
        self.tin = tin

        #: 0-Set by cell precision
        self.tout = tout

        self.ncols = ncols
        self.GCNj = GCNj
        self.GCNi = GCNi

        self.Real = Real
        if len(Complex) or self.is_complex:
            self.Complex = Complex
            if matrix_form not in [1]:  #4, 5, 6, 8
                msg = (
                    f'{self.type} name={name!r} matrix_form={matrix_form!r} '
                    'must be [1, 2, 6]\n'
                    '  1: Square\n'
                    '  2: General Rectangular\n'
                    '  4: Lower Triangular\n'
                    '  5: Upper Triangular\n'
                    '  6: Symmetric\n'
                    '  8: Identity (m=nRows, n=m)\n')
                raise ValueError(msg)

        assert isinstance(matrix_form, integer_types), 'matrix_form=%r type=%s' % (matrix_form, type(matrix_form))
        assert not isinstance(matrix_form, bool), 'matrix_form=%r type=%s' % (matrix_form, type(matrix_form))

    def finalize(self):
        """converts the lists into numpy arrays"""
        return
        #self.GCi = np.asarray(self.GCi)
        #self.GCj = np.asarray(self.GCj)
        self.Real = np.asarray(self.Real)
        if self.is_complex:
            self.Complex = np.asarray(self.Complex)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmiax_to_hdf5(h5_file, model, model.dmiax, encoding)

    @property
    def is_real(self):
        """is the matrix real?"""
        if self.tin in [1, 2]:
            return True
        return False

    @property
    def is_complex(self):
        """is the matrix complex"""
        return not self.is_real

    @property
    def is_polar(self):
        """is the matrix polar (vs real/imag)?"""
        return False

    @property
    def tin_dtype(self):
        """gets the input dtype"""
        return _get_dtype(self.is_complex, self.tin)

    @property
    def tout_dtype(self):
        """gets the output dtype"""
        return _get_dtype(self.is_complex, self.tout)

    @property
    def matrix_type(self):
        """gets the matrix type"""
        if not isinstance(self.matrix_form, integer_types):
            msg = 'ifo must be an integer; matrix_form=%r type=%s name=%s' % (
                self.matrix_form, type(self.matrix_form), self.name)
            raise TypeError(msg)
        if isinstance(self.matrix_form, bool):
            msg = 'matrix_form must not be a boolean; matrix_form=%r type=%s name=%s' % (
                self.matrix_form, type(self.matrix_form), self.name)
            raise TypeError(msg)

        if self.matrix_form == 1:
            matrix_type = 'square'
        #elif self.matrix_form == 6:
            #matrix_type = 'symmetric'
        #elif self.matrix_form in [2, 9]:
            #matrix_type = 'rectangular'
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'matrix_form=%s; ncol' % self.matrix_form)
            raise NotImplementedError(f'{self.type} matrix_form={self.matrix_form} '
                                      'is not supported')
        return matrix_type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a NastranMatrix (DMIAX) card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        #zero

        matrix_form = integer(card, 3, 'ifo')
        tin = integer(card, 4, 'tin')
        tout = integer_or_blank(card, 5, 'tout', 0)
        if matrix_form == 1: # square
            ncols = integer_or_blank(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        elif matrix_form == 6: # symmetric
            ncols = integer_or_blank(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        elif matrix_form in [2, 9]: # rectangular
            ncols = integer(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'matrix_form=%s; ncol' % self.matrix_form)
            raise NotImplementedError('matrix_form=%s is not supported' % matrix_form)

        GCj = []
        GCi = []
        Real = []
        Complex = []
        return DMIAX(name, matrix_form, tin, tout, ncols,
                     GCj, GCi, Real, Complex, comment=comment)

    def _add_column(self, card, comment=''):
        if comment:
            if hasattr(self, '_comment'):
                self.comment += comment
            else:
                self.comment = comment

        unused_name = string(card, 1, 'name')

        Gj = integer(card, 2, 'Gj')
        # Cj = integer(card, 3, 'Cj')
        Cj = integer_or_blank(card, 3, 'Cj', 0)
        #Cj = parse_components(card, 3, 'Cj')
        Nj = integer_or_blank(card, 4, 'Nj')

        assert 0 <= Cj <= 6, 'C%i must be between [0, 6]; Cj=%s' % (0, Cj)

        nfields = len(card)
        #print("nfields = %i" % nfields)
        #print("card[5:] =", card[5:])
        #print("(nfields - 5) %% 4 = %i" % ((nfields - 5) % 4))

        nloops = (nfields - 8) // 8
        if nfields - 8 % 8:
            nloops += 1
        #assert nfields <= 8,'nfields=%s' % nfields
        #print("nloops = %i" % nloops)
        assert nloops > 0, 'nloops=%s' % nloops

        self.GCNj.append((Gj, Cj, Nj))
        GCNi = []
        self.GCNi.append(GCNi)
        if self.is_complex:
            for i in range(nloops):
                #print(dir(card))
                n = 9 + 8 * i
                Gi = integer(card, n, f'Gi{i}')
                # Ci = integer(card, n + 1, 'Ci')
                Ci = integer_or_blank(card, n + 1, f'Ci{i}', 0)
                #Ci = parse_components(card, n + 1, 'Ci')
                Ni = integer_or_blank(card, n + 2, f'Ni{i}')

                assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                GCNi.append((Gi, Ci, Ni))
                reali = double(card, n + 3, 'real')
                complexi = double(card, n + 4, 'complex')
                self.Real.append(reali)
                self.Complex.append(complexi)
        else:
            # real
            for i in range(nloops):
                n = 9 + 9 * i
                Gi = integer(card, n, 'Gi')
                # Ci = integer(card, n + 1, 'Ci')
                Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                #Ci = parse_components(card, n + 1, 'Ci')
                Ni = integer(card, n + 2, 'Ni')

                assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                reali = double(card, n + 3, 'real')
                GCNi.append((Gi, Ci, Ni))
                self.Real.append(reali)
                #print("GC=%s,%s real=%s" % (Gi, Ci, reali))

        msg = '(len(GCNj)=%s len(GCNi)=%s' % (len(self.GCNj), len(self.GCNi))
        assert len(self.GCNj) == len(self.GCNi), msg
        #if self.is_complex:
            #self.Complex(double(card, v, 'complex')

    def raw_fields(self):
        list_fields = [
            'DMIAX', self.name, 0, self.matrix_form, self.tin, None, None, None, self.ncols,
        ]
        k = 0
        if self.is_real:
            for i, GCNj in enumerate(self.GCNj):
                gj, cj, nj = GCNj
                list_fields += ['DMIAX', self.name, gj, cj, nj, None, None, None, None]
                for unused_j, GCNi in enumerate(self.GCNi[i]):
                    gi, ci, ni = GCNi
                    reali = self.Real[k]
                    list_fields += [gi, ci, ni, reali, None, None, None, None]
                    k += 1
        else:
            for i, GCNj in enumerate(self.GCNj):
                gj, cj, nj = GCNj
                list_fields += ['DMIAX', self.name, gj, cj, nj, None, None, None, None]
                for unused_j, GCNi in enumerate(self.GCNi[i]):
                    gi, ci, ni = GCNi
                    reali = self.Real[k]
                    imagi = self.Complex[k]
                    list_fields += [gi, ci, ni, reali, imagi, None, None, None, None]
                    k += 1

        self.write_card()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        if self.tin in [1, 3]:
            is_double = False
            msg = self.write_card_8()
        elif self.tin in [2, 4]:
            is_double = True
            size = 16
            msg = self.write_card_16()
        else:
            raise RuntimeError('tin=%r must be 1, 2, 3, or 4' % self.tin)
        return msg

    def write_card_8(self):
        """writes the card in small field format"""
        return self._write_card(print_card_8)

    def write_card_16(self):
        """writes the card in small large format"""
        return self._write_card(print_card_16)

    def _write_card(self, func):
        """writes the card"""
        msg = '\n$' + '-' * 80
        msg += f'\n$ DMIAX Matrix {self.name}\n'
        list_fields = [
            'DMIAX', self.name, 0, self.matrix_form, self.tin, None, None, None, self.ncols,
        ]
        msg += func(list_fields)
        k = 0
        assert len(self.GCNj) > 0, self.get_stats()
        assert len(self.GCNi) > 0, self.get_stats()
        if self.is_real:
            for i, GCNj in enumerate(self.GCNj):
                gj, cj, nj = GCNj
                list_fields = ['DMIAX', self.name, gj, cj, nj, None, None, None, None]
                for unused_j, GCNi in enumerate(self.GCNi[i]):
                    gi, ci, ni = GCNi
                    reali = self.Real[k]
                    list_fields += [gi, ci, ni, reali, None, None, None, None]
                    k += 1
                msg += func(list_fields)
        else:
            for i, GCNj in enumerate(self.GCNj):
                gj, cj, nj = GCNj
                list_fields = ['DMIAX', self.name, gj, cj, nj, None, None, None, None]
                for unused_j, GCNi in enumerate(self.GCNi[i]):
                    gi, ci, ni = GCNi
                    reali = self.Real[k]
                    imagi = self.Complex[k]
                    list_fields += [gi, ci, ni, reali, imagi, None, None, None]
                    k += 1
                msg += func(list_fields)
        return msg

    def __repr__(self):
        return self.write_card(size=8)

class DMIJ(NastranMatrix):
    """
    Direct Matrix Input at js-Set of the Aerodynamic Mesh
    Defines direct input matrices related to collation degrees-of-freedom
    (js-set) of aerodynamic mesh points for CAERO1, CAERO3, CAERO4 and CAERO5
    and for the slender body elements of CAERO2. These include W2GJ, FA2J and
    input pressures and downwashes associated with AEPRESS and AEDW entries.
    The matrix is described by a single header entry and one or more column
    entries. A column entry is required for each column with nonzero elements.
    For entering data for the interference elements of a CAERO2, use DMIJI
    or DMI.

    """
    type = 'DMIJ'
    _properties = ['shape', 'ifo', 'is_real', 'is_complex', 'is_polar', 'matrix_type',
                   'tin_dtype', 'tout_dtype']

    @classmethod
    def _init_from_empty(cls):
        name = 'name'
        ifo = 1
        tin = 1
        tout = 1
        polar = 0
        ncols = 1
        GCj = []
        GCi = []
        Real = []
        return DMIJ(name, ifo, tin, tout, polar, ncols, GCj, GCi, Real,
                    Complex=None, comment='', finalize=True)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmig_to_hdf5(h5_file, model, model.dmij, encoding)

    def __init__(self, name, matrix_form, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment='',
                 finalize=True):
        """
        Creates a DMIJ card

        Parameters
        ----------
        name : str
            the name of the matrix
        matrix_form : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : List[(node, dof)]???
            the jnode, jDOFs
        GCi  : List[(node, dof)]???
            the inode, iDOFs
        Real : List[float]???
            The real values
        Complex : List[float]???; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        NastranMatrix.__init__(self, name, matrix_form, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment,
                               finalize=finalize)


class DMIJI(NastranMatrix):
    """
    Direct Matrix Input at js-Set of the Interference Body
    Defines direct input matrices related to collation degrees-of-freedom
    (js-set) of aerodynamic mesh points for the interference elements of CAERO2.
    These include W2GJ, FA2J and input pressures and downwashes associated with
    AEPRESS and AEDW entries. The matrix is described by a single header entry
    and one or more column entries. A column entry is required for each column
    with nonzero elements.  For entering data for the slender elements of a
    CAERO2, or a CAERO1, 3, 4 or 5 use DMIJ or DMI.

    """
    type = 'DMIJI'
    _properties = ['shape', 'ifo', 'is_real', 'is_complex', 'is_polar', 'matrix_type',
                   'tin_dtype', 'tout_dtype']

    #@classmethod
    #def _init_from_empty(cls):
        #name = 'name'
        #ifo = 1
        #tin = 1
        #tout = 1
        #polar = 0
        #ncols = 1
        #GCj = []
        #GCi = []
        #Real = []
        #return DMIJI(name, ifo, tin, tout, polar, ncols, GCj, GCi, Real,
                     #Complex=None, comment='', finalize=True)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmig_to_hdf5(h5_file, model, model.dmiji, encoding)

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment='', finalize=True):
        """
        Creates a DMIJI card

        Parameters
        ----------
        name : str
            the name of the matrix
        ifo : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : List[(node, dof)]???
            the jnode, jDOFs
        GCi  : List[(node, dof)]???
            the inode, iDOFs
        Real : List[float]???
            The real values
        Complex : List[float]???; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment,
                               finalize=finalize)


class DMIK(NastranMatrix):
    """
    Direct Matrix Input at ks-Set of the Aerodynamic Mesh
    Defines direct input matrices related to physical (displacement)
    degrees-of-freedom (ks-set) of aerodynamic grid points. These include WKK,
    WTFACT and input forces associated with AEFORCE entries. The matrix is
    described by a single header entry and one or more column entries. A column
    entry is required for each column with nonzero elements.

    +------+-------+----+-----+-----+------+-------+----+------+
    |   1  |   2   | 3  |  4  |  5  |   6  |   7   | 8  |  9   |
    +======+=======+====+=====+=====+======+=======+====+======+
    | DMIK | NAME  | 0  | IFO | TIN | TOUT | POLAR |    | NCOL |
    +------+-------+----+-----+-----+------+-------+----+------+
    | DMIK | NAME  | GJ | CJ  |     |  G1  |  C1   | A1 |  B1  |
    +------+-------+----+-----+-----+------+-------+----+------+
    |      |  G2   | C2 | A2  |  B2 |      |       |    |      |
    +------+-------+----+-----+-----+------+-------+----+------+
    | DMIK | ALPH1 | 0  |  9  |  2  |  0   |   1   |    |      |
    +------+-------+----+-----+-----+------+-------+----+------+
    | DMIK | ALPH1 | 1  |  1  |  1  |  1   |  1.0  |    |      |
    +------+-------+----+-----+-----+------+-------+----+------+
    |      |   2   | 1  | 1.0 |     |      |       |    |      |
    +------+-------+----+-----+-----+------+-------+----+------+
    """
    type = 'DMIK'
    _properties = ['shape', 'ifo', 'is_real', 'is_complex', 'is_polar', 'matrix_type',
                   'tin_dtype', 'tout_dtype']

    #@classmethod
    #def _init_from_empty(cls):
        #name = 'name'
        #ifo = 1
        #tin = 1
        #tout = 1
        #polar = 0
        #ncols = 1
        #GCj = []
        #GCi = []
        #Real = []
        #return DMIK(name, ifo, tin, tout, polar, ncols, GCj, GCi, Real,
                    #Complex=None, comment='', finalize=True)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmig_to_hdf5(h5_file, model, model.dmik, encoding)

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment='', finalize=True):
        """
        Creates a DMIK card

        Parameters
        ----------
        name : str
            the name of the matrix
        ifo : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : List[(node, dof)]
            the jnode, jDOFs
        GCi  : List[(node, dof)]
            the inode, iDOFs
        Real : List[float]
            The real values
        Complex : List[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment,
                               finalize=finalize)


class DMI(NastranMatrix):
    """
    +------+-------+------+------+---------+----------+-----------+-----------+------+
    |  1   |   2   |  3   |   4  |    5    |    6     |     7     | 8         |  9   |
    +======+=======+======+======+=========+==========+===========+===========+======+
    | DMI  |  NAME |  0   | FORM |   TIN   |   TOUT   |           |     M     |  N   |
    +------+-------+------+------+---------+----------+-----------+-----------+------+
    | DMI  |  NAME |  J   |  I1  | A(I1,J) |  A(I1,J) | A(I1+1,J) | A(I1+2,J) | etc. |
    +------+-------+------+------+---------+----------+-----------+-----------+------+
    |      |  I2   | etc. |      |         |          |           |           |      |
    +------+-------+------+------+---------+----------+-----------+-----------+------+
    """
    type = 'DMI'
    _properties = ['shape', 'ifo', 'is_real', 'is_complex', 'is_polar', 'matrix_type',
                   'tin_dtype', 'tout_dtype']

    @classmethod
    def _init_from_empty(cls):
        name = 'name'
        matrix_form = 8
        tin = 1
        tout = 1
        nrows = 5
        ncols = 5
        GCj = []
        GCi = []
        Real = []
        return DMI(name, matrix_form, tin, tout, nrows, ncols, GCj, GCi, Real,
                   Complex=None, comment='', finalize=False)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        _export_dmig_to_hdf5(h5_file, model, model.dmi, encoding)

    def __init__(self, name, matrix_form, tin, tout, nrows, ncols,
                 GCj, GCi, Real, Complex=None, comment='', finalize=True):
        #NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               #GCj, GCi, Real, Complex, comment='')
        if comment:
            self.comment = comment

        if Complex is None:
            Complex = []

        if tout is None:
            tout = 0

        if matrix_form not in [1, 2, 3, 4, 5, 6, 8]:
            msg = (
                '%s name=%r matrix_form=%r must be [1, 2, 3, 4, 5, 6, 8]\n'
                '  1: Square\n'
                '  2: Rectangular\n'
                '  3: Diagonal matrix (M=number of rows, N=1)\n'
                '  4: Lower Triangular\n'
                '  5: Upper Triangular\n'
                '  6: Symmetric\n'
                '  8: Identity (m=nRows, n=m)\n'
                #'  9: Rectangular\n'
                % (self.type, name, matrix_form))
            raise ValueError(msg)

        self.name = name
        self.matrix_form = matrix_form
        self.tin = tin
        self.tout = tout
        self.nrows = nrows
        self.ncols = ncols
        self.GCi = GCi
        self.GCj = GCj
        self.Real = Real
        if len(Complex) or self.is_complex:
            self.Complex = Complex
        if finalize:
            self.finalize()

    #@property
    #def form(self):
        #"""gets the matrix_form"""
        #self.deprecated('form', 'matrix_form', '1.1')
        #return self.matrix_form

    #@form.setter
    #def form(self, matrix_form):
        #"""sets the matrix_form"""
        #self.deprecated('form', 'matrix_form', '1.1')
        #self.matrix_form = matrix_form

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DMI card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        #zero

        #: Form of the matrix:  1=Square (not symmetric); 2=Rectangular;
        #: 3=Diagonal (m=nRows,n=1);  4=Lower Triangular; 5=Upper Triangular;
        #: 6=Symmetric; 8=Identity (m=nRows, n=m)
        matrix_form = integer(card, 3, 'form')

        #: 1-Real, Single Precision; 2=Real,Double Precision;
        #: 3=Complex, Single; 4=Complex, Double
        tin = integer(card, 4, 'tin')

        #: 0-Set by cell precision
        tout = integer_or_blank(card, 5, 'tout', 0)

        nrows = integer(card, 7, 'nrows')
        ncols = integer(card, 8, 'ncols')

        assert len(card) == 9, 'len(DMI card) = %i\ncard=%s' % (len(card), card)

        GCj = []
        GCi = []
        Real = []
        Complex = []
        return DMI(name, matrix_form, tin, tout, nrows, ncols,
                   GCj, GCi, Real, Complex, comment=comment, finalize=False)

    def finalize(self):
        self.GCi = np.asarray(self.GCi)
        self.GCj = np.asarray(self.GCj)
        self.Real = np.asarray(self.Real)
        if self.is_complex:
            self.Complex = np.asarray(self.Complex)

    @property
    def matrix_type(self):
        """
        gets the matrix type

        1 Square matrix (not symmetric)
        2 General rectangular matrix
        3 Diagonal matrix (M=number of rows, N = 1)
        #4 Lower triangular factor
        #5 Upper triangular factor
        6 Symmetric matrix
        8 Identity matrix (M=number of rows, N = M)
        """
        if not isinstance(self.matrix_form, integer_types):
            msg = 'ifo must be an integer; matrix_form=%r type=%s name=%s' % (
                self.matrix_form, type(self.matrix_form), self.name)
            raise TypeError(msg)
        if isinstance(self.matrix_form, bool):
            msg = 'matrix_form must not be a boolean; matrix_form=%r type=%s name=%s' % (
                self.matrix_form, type(self.matrix_form), self.name)
            raise TypeError(msg)

        if self.matrix_form == 1:
            matrix_type = 'square'
        elif self.matrix_form == 2: # 9 ???
            matrix_type = 'rectangular'
        elif self.matrix_form == 3:
            matrix_type = 'diagonal'
        elif self.matrix_form == 6:
            matrix_type = 'symmetric'
        elif self.matrix_form == 9:
            matrix_type = 'identity'
        else:
            raise NotImplementedError('%s matrix_form=%r is not supported' % (
                self.type, self.matrix_form))
        return matrix_type

    @property
    def is_polar(self):
        if self.tin in [1, 2]:
            is_polar = False
        elif self.tin in [3, 4]:
            is_polar = False # TODO: could be wrong...
        else:
            raise NotImplementedError('nrows=%s ncols=%s' % (self.nrows, self.ncols))
        return is_polar

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def ifo(self):
        """
        ifo
        #: 4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)

        #: Form of the matrix:  1=Square (not symmetric); 2=Rectangular;
        #: 3=Diagonal (m=nRows,n=1);  4=Lower Triangular; 5=Upper Triangular;
        #: 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.matrix_form = integer(card, 3, 'matrix_form')

        """
        return self.matrix_form
        #if self.nrows == self.ncols:
            ## symmetric
            #ifo = 6
        ##elif self.nrows > 1 and self.ncols > 1:
            ##ifo = 2
        #else:
            #raise NotImplementedError('matrix_form=%r nrows=%s ncols=%s' % (
                #self.matrix_form, self.nrows, self.ncols))
        #return ifo

    def _add_column(self, card, comment=''):
        """
        .. todo:: support comment
        """
        if self.is_complex:
            self._read_complex(card)
        else:
            self._read_real(card)

    def _read_real(self, card):
        """reads a real DMI column"""
        # column number
        j = integer(card, 2, 'icol')

        # counter
        i = 0
        fields = [interpret_value(field, card) for field in card[3:]]

        # Real, starts at A(i1,j), goes to A(i2,j) in a column
        while i < len(fields):
            i1 = fields[i]
            if isinstance(i1, integer_types):
                i += 1
                is_done_reading_floats = False
                while not is_done_reading_floats and i < len(fields):
                    real_value = fields[i]
                    if isinstance(real_value, integer_types):
                        is_done_reading_floats = True
                    elif isinstance(real_value, float):
                        #print('adding j=%s i1=%s val=%s' % (j, i1, real_value))
                        self.GCj.append(j)
                        self.GCi.append(i1)
                        self.Real.append(real_value)
                        i += 1
                        i1 += 1
                    else:
                        real_value = self.Real[-1]
                        end_i = fields[i + 1]
                        for ii in range(i1, end_i + 1):
                            #print('adding j=%s i1=%s val=%s' % (j, ii, real_value))
                            self.GCj.append(j)
                            self.GCi.append(ii)
                            self.Real.append(real_value)
                        i += 1
                        is_done_reading_floats = True

    def _read_complex(self, card):
        """reads a complex DMI column"""
        #msg = 'complex matrices not supported in the DMI reader...'
        #raise NotImplementedError(msg)
        # column number
        j = integer(card, 2, 'icol')
        # counter
        i = 0
        fields = [interpret_value(field, card) for field in card[3:]]
        # Complex, starts at A(i1,j)+imag*A(i1,j), goes to A(i2,j) in a column
        if 0: # pragma: no cover
            is_real = True
            gci = None
            for field in fields:
                if isinstance(field, integer_types):
                    gci = field
                elif isinstance(field, float):
                    if is_real:
                        real = field
                    else:
                        self.GCj.append(j)
                        self.GCi.append(gci)
                        self.Real.append(real)
                        self.Complex.append(field)
                    is_real = not is_real

        while i < len(fields):
            i1 = fields[i]
            assert isinstance(i1, int), card
            i += 1
            is_done_reading_floats = False
            while not is_done_reading_floats and i < len(fields):
                value = fields[i]
                #print("i=%s len(fields)=%s value=%s" % (
                    #i, len(fields), value))
                if isinstance(value, integer_types):
                    is_done_reading_floats = True
                elif isinstance(value, float):
                    complex_value = fields[i + 1]
                    assert isinstance(complex_value, float), card
                    self.GCj.append(j)
                    self.GCi.append(i1)
                    self.Real.append(value)
                    self.Complex.append(complex_value)
                    i += 2
                else:
                    raise NotImplementedError()

    @property
    def is_real(self):
        """real vs. complex attribute"""
        return not self.is_complex

    @property
    def is_complex(self):
        """real vs. complex attribute"""
        if self.tin in [3, 4]:
            return True
        return False

    def raw_fields(self):
        """
        .. warning:: All the writers are bad because Nastran insists on
                      making columns a single DMI card.  This makes
                      writing a card much harder, so there are a lot of
                      NotImplementedErrors floating about.

                      This is an invalid method, but is not disabled
                      because it's currently needed for checking results

        """
        list_fields = ['DMI', self.name, 0, self.matrix_form, self.tin,
                       self.tout, None, self.nrows, self.ncols]

        if self.is_complex:
            for (gci, gcj, reali, imagi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                list_fields += ['DMI', self.name, gcj, gci, reali, imagi]
        else:
            for (gci, gcj, reali) in zip(self.GCi, self.GCj, self.Real):
                list_fields += ['DMI', self.name, gcj, gci, reali]
        return list_fields

    def write_card_8(self):
        """writes the card in single precision"""
        return self._write_card(print_card_8)

    def _get_real_fields(self, func):
        msg = ''
        uGCj = np.unique(self.GCj)
        for gcj in uGCj:
            i = np.where(gcj == self.GCj)[0]
            gcis = self.GCi[i]
            reals = self.Real[i]
            isort = np.argsort(gcis)
            list_fields = ['DMI', self.name, gcj]

            # will always write the first one
            gci_last = -1
            for gci, real in zip(gcis[isort], reals[isort]):
                if gci == gci_last + 1:
                    pass
                else:
                    list_fields.append(gci)
                list_fields.append(real)
                gci_last = gci
            msg += func(list_fields)
        return msg

    def _get_complex_fields(self, func):
        msg = ''
        uGCj = np.unique(self.GCj)
        for gcj in uGCj:
            i = np.where(gcj == self.GCj)[0]
            gcis = self.GCi[i]
            reals = self.Real[i]
            complexs = self.Complex[i]
            isort = np.argsort(gcis)
            list_fields = ['DMI', self.name, gcj]

            # will always write the first one
            gci_last = -10
            #print('gcis=%s \nreals=%s \ncomplexs=%s' % (
                #gcis[isort], reals[isort], complexs[isort]))
            if max(gcis) == min(gcis):
                list_fields += [gcis[0]]
                for reali, complexi in zip(reals, complexs):
                    list_fields.extend([reali, complexi])
                msg += func(list_fields)
            else:
                #print(f'list_fields0 = {list_fields}')
                for i, gci, reali, complexi in zip(count(), gcis[isort], reals[isort], complexs[isort]):
                    #print('B', gci, reali, complexi, gci_last)
                    if gci != gci_last + 1 and i != 0:
                        pass
                    else:
                        list_fields.append(gci)
                    list_fields.append(reali)
                    list_fields.append(complexi)
                    gci_last = gci
                #print(f'list_fields = {list_fields}')
                msg += func(list_fields)
        return msg

    def get_matrix(self, is_sparse=False, apply_symmetry=True):
        """
        Builds the Matrix

        Parameters
        ----------
        is_sparse : bool; default=False
            should the matrix be returned as a sparse matrix.
            Slower for dense matrices.
        apply_symmetry : bool; default=True
            If the matrix is symmetric (ifo=6), returns a symmetric matrix.
            Supported as there are symmetric matrix routines.

        Returns
        -------
        M : numpy.ndarray or scipy.coomatrix
            the matrix
        rows : dict[int] = [int, int]
            dictionary of keys=rowID, values=(Grid,Component) for the matrix
        cols: dict[int] = [int, int]
            dictionary of keys=columnID, values=(Grid,Component) for the matrix

        .. warning:: is_sparse=True WILL fail

        """
        return get_dmi_matrix(self, is_sparse=is_sparse, apply_symmetry=apply_symmetry)

    def write_card_16(self):
        """writes the card in single precision"""
        return self._write_card(print_card_16)

    def write_card_double(self):
        """writes the card in double precision"""
        return self._write_card(print_card_16)

    def _write_card(self, func):
        """writes the card in single/double precision"""
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % ('DMI', self.name)
        list_fields = ['DMI', self.name, 0, self.matrix_form, self.tin,
                       self.tout, None, self.nrows, self.ncols]
        msg += print_card_8(list_fields)

        if self.is_complex:
            msg += self._get_complex_fields(func)
        else:
            msg += self._get_real_fields(func)
        return msg

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        if size == 8:
            return self.write_card_8()
        if is_double:
            return self.write_card_double()
        return self.write_card_16()

    def __repr__(self):
        """
        .. todo:: support shortened output format.  There's a very low 1000
                  DMI cap, I assume this is entries and not matrices.

        """
        return self.write_card(size=8, is_double=False)


def get_row_col_map(matrix, GCi, GCj, ifo):
    ndim = len(GCi.shape)
    #print('ndim=%s' % ndim)
    #print('GCj=%s' % GCj)
    #print('GCi=%s' % GCi)
    if ndim == 1:
        rows, cols, rows_reversed, cols_reversed = _get_row_col_map_1d(matrix, GCi, GCj, ifo)
    else:
        rows, cols, rows_reversed, cols_reversed = _get_row_col_map_2d(matrix, GCi, GCj, ifo)

    nrows = len(rows)
    ncols = len(cols)
    assert nrows > 0, 'nrows=%s' % nrows
    assert ncols > 0, 'ncols=%s' % ncols
    return nrows, ncols, ndim, rows, cols, rows_reversed, cols_reversed

def _get_row_col_map_1d(matrix, GCi, GCj, ifo):
    """helper for ``get_row_col_map``"""
    rows = {}
    rows_reversed = {}

    cols = {}
    cols_reversed = {}
    i = 0
    #nrows = np.unique(GCi)
    #ncols = np.unique(GCj)
    for gci in GCi:
        if gci not in rows:
            rows[gci] = i
            rows_reversed[i] = gci
            i += 1

    if ifo == 6:
        # symmetric
        #print(GCj)
        for gcj in GCj:
            if gcj not in rows:
                #print('row.gcj = %s' % str(gcj))
                rows[gcj] = i
                rows_reversed[i] = gcj
                i += 1
        cols = rows
        cols_reversed = rows_reversed
    else:
        j = 0
        for gcj in GCj:
            if gcj not in cols:
                cols[gcj] = j
                cols_reversed[j] = gcj
                j += 1
    return rows, cols, rows_reversed, cols_reversed

def _get_row_col_map_2d(matrix, GCi, GCj, ifo):
    """helper for ``get_row_col_map``"""
    rows = {}
    rows_reversed = {}

    cols = {}
    cols_reversed = {}
    #print('i0=%s j0=%s' % (i, j))
    #nrows = len(GCi)
    #ncols = len(GCj)
    #rows_array = np.zeros((nrows, 2), dtype='int32')
    #cols_array = np.zeros((ncols, 2), dtype='int32')
    #for i, (nid, comp) in enumerate(GCi):
        ##print('i=%s nid=%s comp=%s nrows=%s rows_array.shape=%s' % (
            ##i, nid, comp, nrows, str(rows_array.shape)))
        #rows_array[i, :] = [nid, comp]
    #print('rows_array = \n%s' % rows_array)

    #for j, (nid, comp) in enumerate(GCj):
        #cols_array[j, :] = [nid, comp]
    #print('cols_array = \n%s' % cols_array)

    #print(GCi)
    #print(GCj)
    i = 0
    for (nid, comp) in GCi:
        gci = (nid, comp)
        if gci not in rows:
            #print('row.gci = %s' % str(gci))
            rows[gci] = i
            rows_reversed[i] = gci
            i += 1
    if ifo == 6:
        # symmetric
        for (nid, comp) in GCj:
            gcj = (nid, comp)
            if gcj not in rows:
                #print('row.gcj = %s' % str(gcj))
                rows[gcj] = i
                rows_reversed[i] = gcj
                i += 1
        cols = rows
        cols_reversed = rows_reversed
    else:
        j = 0
        for (nid, comp) in GCj:
            gcj = (nid, comp)
            if gcj not in cols:
                #print('col.gcj = %s' % str(gcj))
                cols[gcj] = j
                cols_reversed[j] = gcj
                j += 1
    return rows, cols, rows_reversed, cols_reversed

def _fill_sparse_matrix(matrix, nrows, ncols):
    """helper method for get_matrix"""
    GCj = array(matrix.GCj, dtype='int32') - 1
    GCi = array(matrix.GCi, dtype='int32') - 1
    reals = array(matrix.Real, dtype='float32')

    # TODO: matrix size:  is this correct?
    nrows = max(GCi) + 1
    ncols = max(GCj) + 1

    dtype = _get_dtype(matrix.is_complex, matrix.tin)
    # TODO: no check for symmetry
    # TODO: no check for dtype
    if matrix.is_complex:
        complexs = array(matrix.Complex, dtype='float32')
        data = array([reals, complexs]).astype(complex)
    else:
        data = reals

    if matrix.matrix_form in [1, 6]:
        nrows = max(nrows, ncols)
        ncols = nrows

    #A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype=dtype) # test
    sparse_matrix = coo_matrix((data, (matrix.GCi, matrix.GCj)),
                               shape=(nrows, ncols), dtype=dtype)
    #sparse_matrix = coo_matrix( (data,(matrix.GCi,matrix.GCj)),shape=(i,j)) # old
    #sparse_matrix = coo_matrix( (data,(matrix.GCi,matrix.GCj)),shape=(nrows,ncols))
    #print(sparse_matrix.toarray())
    #print(sparse_matrix)
    return sparse_matrix


def _fill_dense_rectangular_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry):
    """helper method for get_matrix"""
    is_sparse = False
    if self.is_complex:
        dense_mat = zeros((nrows, ncols), dtype='complex128')
        if self.matrix_form == 6 and apply_symmetry:  # symmetric
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[(gci[0], gci[1])]
                j = cols[(gcj[0], gcj[1])]
                dense_mat[i, j] = complex(reali, complexi)
                dense_mat[j, i] = complex(reali, complexi)
        else:
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[(gci[0], gci[1])]
                j = cols[(gcj[0], gcj[1])]
                dense_mat[i, j] = complex(reali, complexi)
    else:
        dense_mat = zeros((nrows, ncols), dtype='float64')
        if self.matrix_form == 6 and apply_symmetry:  # symmetric
            try:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[(gci[0], gci[1])]
                    j = cols[(gcj[0], gcj[1])]
                    dense_mat[i, j] = reali
                    dense_mat[j, i] = reali
            except IndexError:
                msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                       'is_polar=%s ncols=%s M.shape=%s\n' % (
                           self.name, ndim, i, j, self.matrix_type,
                           self.is_polar, self.ncols, dense_mat.shape))
                msg += 'Rows:\n'
                for i, row in enumerate(rows):
                    msg += 'i=%s row=%s\n' % (i, row)
                raise RuntimeError(msg)
        else:
            try:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[(gci[0], gci[1])]
                    j = cols[(gcj[0], gcj[1])]
                    dense_mat[i, j] = reali
            except KeyError:
                msg = ('name=%s ndim=%s gci=%s gcj=%s matrix_type=%s '
                       'is_polar=%s is_sparse=%s ncols=%s M.shape=%s\n\n' % (
                           self.name, ndim, str(gci), str(gcj), self.matrix_type,
                           self.is_polar, is_sparse, self.ncols, dense_mat.shape))

                gci2 = (gci[0], gci[1])
                gcj2 = (gcj[0], gcj[1])
                if gci2 in rows:
                    msg += 'gci/row_key=%s found\n' % str(gci2)
                else:
                    msg += 'gci/row_key=%s not found\n' % str(gci2)
                    msg += 'Rows:\n'
                    for i, row in enumerate(rows):
                        msg += '  i=%s row=%s\n' % (i, row)

                if gcj2 in cols:
                    msg += '\ngcj/col_key=%s found\n' % str(gcj2)
                else:
                    msg += '\ngcj/col_key=%s not found\n' % str(gcj2)
                    msg += 'Cols:\n'
                    for j, col in enumerate(cols):
                        msg += '  j=%s row=%s\n' % (j, col)

                msg += '\n'
                print(msg)
                raise KeyError(msg)
            except IndexError:
                msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                       'is_polar=%s is_sparse=%s ncols=%s M.shape=%s\n' % (
                           self.name, ndim, i, j, self.matrix_type,
                           self.is_polar, is_sparse, self.ncols, dense_mat.shape))
                msg += 'Rows:\n'
                for i, row in enumerate(rows):
                    msg += '  i=%s row=%s\n' % (i, row)

                msg += '\nCols:\n'
                for j, row in enumerate(cols):
                    msg += '  j=%s row=%s\n' % (j, col)
                raise RuntimeError(msg)
    return dense_mat


def _fill_dense_column_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry):
    """helper method for get_matrix"""
    is_sparse = False
    if self.is_complex:
        dense_mat = zeros((nrows, ncols), dtype='complex128')
        if self.matrix_form == 6 and apply_symmetry:  # symmetric
            assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[gci]
                j = cols[gcj]
                dense_mat[i, j] = complex(reali, complexi)
                dense_mat[j, i] = complex(reali, complexi)
        elif self.matrix_form == 2:  # rectangular
            assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[gci]
                j = cols[gcj]
                dense_mat[i, j] = complex(reali, complexi)
        else:
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[gci]
                j = cols[gcj]
    else:
        #print('nrows=%s ncols=%s' % (nrows, ncols))
        dense_mat = zeros((nrows, ncols), dtype='float64')
        if self.matrix_form == 6 and apply_symmetry:  # symmetric
            assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
            for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                i = rows[gci]
                j = cols[gcj]
                dense_mat[i, j] = reali
                dense_mat[j, i] = reali
        else:
            try:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[gci]
                    j = cols[gcj]
                    dense_mat[i, j] = reali
            except IndexError:
                msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                       'is_polar=%s is_sparse=%s ncols=%s M.shape=%s\n' % (
                           self.name, ndim, i, j, self.matrix_type,
                           self.is_polar, is_sparse, self.ncols, dense_mat.shape))
                msg += 'Rows:\n'
                for i, row in enumerate(rows):
                    msg += '  i=%s row=%s\n' % (i, row)
                raise RuntimeError(msg)
    return dense_mat

def get_dmi_matrix(matrix: DMI, is_sparse: bool=False,
                   apply_symmetry: bool=True) -> Tuple[np.array, None, None]:
    """
    Builds the Matrix

    Parameters
    ----------
    is_sparse : bool
        should the matrix be returned as a sparse matrix (default=True).
        Slower for dense matrices.
    apply_symmetry: bool
        If the matrix is symmetric (matrix_form=6), returns a symmetric matrix.
        Supported as there are symmetric matrix routines.
        TODO: unused...

    Returns
    -------
    M : ndarray
        the matrix
    rows : None
        unused
    cols : None
        unused

    .. warning:: is_sparse=True WILL fail

    """
    ifo = matrix.ifo
    GCj = array(matrix.GCj, dtype='int32') - 1
    GCi = array(matrix.GCi, dtype='int32') - 1

    dtype = matrix.tin_dtype

    if matrix.is_complex:
        data = matrix.Real + matrix.Complex * 1j
    else:
        data = matrix.Real

    if ifo == 2:
        # rectangular
        nrows = matrix.nrows
        ncols = matrix.ncols

        M = coo_matrix((data, (GCi, GCj)),
                       shape=(nrows, ncols), dtype=dtype)
        if not is_sparse:
            M = M.toarray()

    else:
        nrows = matrix.nrows
        ncols = matrix.ncols
        if ifo == 6:
            nrows = max(nrows, ncols)
            ncols = nrows
        M = coo_matrix((data, (GCi, GCj)),
                       shape=(nrows, ncols), dtype=dtype)
        if not is_sparse:
            M = M.toarray()
    #else:
        #ifo : int
        #    matrix shape
        #    4=Lower Triangular
        #    5=Upper Triangular
        #    6=Symmetric
        #    8=Identity (m=nRows, n=m)
        #raise RuntimeError(matrix.get_stats())
    return M, None, None

def get_matrix(self, is_sparse=False, apply_symmetry=True):
    """
    Builds the Matrix

    Parameters
    ----------
    is_sparse : bool
        should the matrix be returned as a sparse matrix (default=True).
        Slower for dense matrices.
    apply_symmetry: bool
        If the matrix is symmetric (matrix_form=6), returns a symmetric matrix.
        Supported as there are symmetric matrix routines.
        TODO: unused...

    Returns
    -------
    M : ndarray
        the matrix
    rows : Dict[(nid, nid)] = float
        dictionary of keys=rowID,    values=(Grid,Component) for the matrix
    cols : Dict[](int, int)] = float
        dictionary of keys=columnID, values=(Grid,Component) for the matrix

    .. warning:: is_sparse=True WILL fail

    """
    nrows, ncols, ndim, rows, cols, rows_reversed, cols_reversed = get_row_col_map(
        self, self.GCi, self.GCj, self.matrix_form)
    #print('rows = ', rows)
    #print('cols = ', cols)
    #print('i=%s j=%s' % (i, j))
    #nrows = len(rows2)
    #ncols = len(cols2)

    #A = ss.lil_matrix((3,3), dtype='d') # double precision
    #rows = []
    #cols = []
    #data = []
    #for i in range(3):
       #for j in range(3):
           #k = float((i+1)*(j+1))
           #rows.append(i)
           #cols.append(j)
           #data.append(k)
           #A[i,j] = k

    #is_sparse = False
    if is_sparse:
        M = _fill_sparse_matrix(self, nrows, ncols)
        return M
    else:
        if ndim == 1:
            M = _fill_dense_column_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry)
        else:
            M = _fill_dense_rectangular_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry)

        #print(M)
        return (M, rows_reversed, cols_reversed)


def _export_dmig_to_hdf5(h5_file, model, dict_obj, encoding):
    """export dmigs, dmij, dmiji, dmik, dmi"""
    for name, dmig in dict_obj.items():
        dmig_group = h5_file.create_group(name)
        dmig_group.create_dataset('tin', data=dmig.tin)

        if hasattr(dmig, 'tout'):
            dmig_group.create_dataset('tout', data=dmig.tout)

        if dmig.type == 'DMIG' and name == 'UACCEL':
            if dmig.ncol is not None:
                dmig_group.create_dataset('ncol', data=dmig.ncol)
            #load_seq_group = dmig_group.create_group('load_sequences')

            nids = []
            dofs = []
            values = []
            for lseq, ncx in sorted(dmig.load_sequences.items()):
                lseq_group = dmig_group.create_group(str(lseq))
                #list_fields += [lseq, None, None]
                for (nid, dof, value) in ncx:
                    nids.append(nid)
                    dofs.append(int(dof))
                    values.append(value)

                    #print('nids =', nids)
                    #print('dofs =', dofs)
                    #print('values =', values)
                lseq_group.create_dataset('nids', data=nids)
                lseq_group.create_dataset('dofs', data=dofs)
                lseq_group.create_dataset('values', data=values)
        else:
            if hasattr(dmig, 'nrows') and dmig.nrows is not None:
                dmig_group.create_dataset('nrows', data=dmig.nrows)
            if dmig.ncols is not None:
                dmig_group.create_dataset('ncols', data=dmig.ncols)
            if hasattr(dmig, 'polar'):
                dmig_group.create_dataset('polar', data=dmig.polar)

            dmig_group.create_dataset('matrix_form', data=dmig.matrix_form)
            dmig_group.create_dataset('tin_dtype', data=dmig.tin_dtype)
            dmig_group.create_dataset('tout_dtype', data=dmig.tout_dtype)

            dmig_group.create_dataset('matrix_type', data=dmig.matrix_type)
            dmig_group.create_dataset('is_complex', data=dmig.is_complex)
            dmig_group.create_dataset('is_real', data=dmig.is_real)
            dmig_group.create_dataset('is_polar', data=dmig.is_polar)

            dmig_group.create_dataset('GCi', data=dmig.GCi)
            dmig_group.create_dataset('GCj', data=dmig.GCj)
            dmig_group.create_dataset('Real', data=dmig.Real)
            if hasattr(dmig, 'Complex') and dmig.Complex is not None:
                dmig_group.create_dataset('Complex', data=dmig.Complex)


def _export_dmiax_to_hdf5(h5_file, model, dict_obj, encoding):
    """export dmiax"""
    for name, dmiax in dict_obj.items():
        #print(f'exporting {dmiax.type} name={name!r}')
        dmiax_group = h5_file.create_group(name)
        dmiax_group.create_dataset('tin', data=dmiax.tin)

        if hasattr(dmiax, 'tout'):
            dmiax_group.create_dataset('tout', data=dmiax.tout)

        if hasattr(dmiax, 'nrows') and dmiax.nrows is not None:
            dmiax_group.create_dataset('nrows', data=dmiax.nrows)
        if dmiax.ncols is not None:
            dmiax_group.create_dataset('ncols', data=dmiax.ncols)
        if hasattr(dmiax, 'polar'):
            dmiax_group.create_dataset('polar', data=dmiax.polar)

        dmiax_group.create_dataset('matrix_form', data=dmiax.matrix_form)
        dmiax_group.create_dataset('tin_dtype', data=dmiax.tin_dtype)
        dmiax_group.create_dataset('tout_dtype', data=dmiax.tout_dtype)

        dmiax_group.create_dataset('matrix_type', data=dmiax.matrix_type)
        dmiax_group.create_dataset('is_complex', data=dmiax.is_complex)
        dmiax_group.create_dataset('is_real', data=dmiax.is_real)
        dmiax_group.create_dataset('is_polar', data=dmiax.is_polar)

        gcnj = []
        j_none_flags = []

        gcni = []
        i_none_flags = []
        for j, GCNj in enumerate(dmiax.GCNj):
            gj, cj, nj = GCNj
            is_none_flag_j = False
            if nj is None:
                nj = 0
                is_none_flag_j = True
            j_none_flags.append(is_none_flag_j)
            gcnj.append((gj, cj, nj))
            for unused_i, GCNi in enumerate(dmiax.GCNi[j]):
                gi, ci, ni = GCNi
                is_none_flag_i = False
                if ni is None:
                    ni = 0
                    is_none_flag_i = True
                i_none_flags.append(is_none_flag_i)
                gcni.append((gi, ci, ni, j))

        dmiax_group.create_dataset('GCNi_j', data=gcni)
        dmiax_group.create_dataset('GCNj', data=gcnj)
        dmiax_group.create_dataset('i_none_flags', data=i_none_flags)
        dmiax_group.create_dataset('j_none_flags', data=j_none_flags)

        dmiax_group.create_dataset('Real', data=dmiax.Real)
        if hasattr(dmiax, 'Complex') and dmiax.Complex is not None:
            dmiax_group.create_dataset('Complex', data=dmiax.Complex)

def _set_polar(polar):
    if polar in [None, 0, False]:
        polar = 0
    elif polar in [1, True]:
        polar = 1
    else:
        raise ValueError('polar=%r and must be 0 or 1' % polar)
    return polar

def _get_dtype(is_complex, type_flag):
    if type_flag == 1:
        dtype = 'float32'
    elif type_flag == 2:
        dtype = 'float64'
    elif type_flag == 3:
        dtype = 'complex64'
    elif type_flag == 4:
        dtype = 'complex128'
    elif type_flag == 0:
        if is_complex:
            dtype = 'complex128'
        else:
            dtype = 'float64'
    else:
        raise RuntimeError("invalid option for matrix format")
    return dtype

