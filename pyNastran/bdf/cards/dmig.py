# pylint: disable=R0902,R0904,R0914
from __future__ import annotations
from copy import deepcopy
from math import sin, cos, radians, atan2, sqrt, degrees
from itertools import count
import warnings
from typing import Any, Callable, TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.femutils.utils import unique2d
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs, collapse_thru_ipacks
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, blank, integer_or_blank,
    double, string, string_or_blank,
    parse_components, interpret_value, integer_double_string_or_blank)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.bdf.bdf import BDF


class DTI_UNITS(BaseCard):
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
        name = 'UNITS'
        fields = []
        return DTI_UNITS(name, fields, comment='')

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
        #from pyNastran.bdf.bdf_interface.hdf5_exporter import _export_list
        for name, dti in sorted(model.dti.items()):
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

    def __init__(self, name, fields, comment=''):
        """
        Creates a DTI,UNITS card

        Parameters
        ----------
        name : str
            UNITS
        fields : list[varies]
            the fields
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.name = name
        self.fields = fields
        #print(fields)
        assert len(fields) > 0, fields
        for key, fieldsi in fields.items():
            assert fieldsi is None or isinstance(fieldsi, str), fields

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
            'temp_stress' : temp_stress,
        }
        return DTI_UNITS(name, fields, comment=comment)

    def raw_fields(self):
        mass = self.fields['mass']
        force = self.fields['force']
        length = self.fields['length']
        time = self.fields['time']
        temp_stress = self.fields['temp_stress']
        list_fields = ['DTI', 'UNITS', '1', mass, force, length, time, temp_stress]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

# 0 means the same as tin
# we won't handle it with this
TOUT_DTYPE_MAP = {
    1: 'float32',
    2: 'float64',
    3: 'complex64',
    4: 'complex128',
}
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

    def _finalize_hdf5(self, encoding: str) -> None:
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
        print('values_str', values_str)
        print('fields', self.fields)
        asdf

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding: str):
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
                h5_group = h5_file.create_group(str(name))
                for irecord, fields in sorted(dti.fields.items()):
                    #h5_group = h5_file.create_group(str(irecord))
                    attr = 'irecord=%s' % irecord
                    namei = str(irecord)
                    values = fields
                    _export_list(h5_group, attr, namei, values, encoding)
                    #print(h5_group)
                    #print(irecord, fields)

    def __deepcopy__(self, memo: dict[str, Any]):
        """performs a deepcopy"""
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        bad_attrs = []
        for key, value in self.__dict__.items():
            try:
                memo2 = deepcopy(value, memo)
            except SyntaxError:
                if isinstance(value, dict):
                    for keyi, valuei in value.items():
                        self.log.warn(valuei)
                        #print(valuei.object_attributes())
                        #break
                #raise
                bad_attrs.append(key)
            setattr(result, key, memo2)
        if bad_attrs:
            raise RuntimeError(f'failed copying {bad_attrs}')
        return result

    def __init__(self, name: str, fields: dict[int, list],
                 comment: str=''):
        """
        Creates a DTI card

        Parameters
        ----------
        name : str
            UNITS
        fields : dict[int, list[Any]]
            the fields
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.name = name
        self.fields = fields
        for key, fieldsi in fields.items():
            assert fieldsi is not None, fields
            assert isinstance(fieldsi, list), fieldsi
            for fieldi in fieldsi:
                assert not isinstance(fieldi, bytes), fieldsi
                #assert not isinstance(fieldi, np.ndarray), fieldsi
        assert len(fields) > 0, fields
        assert name != 'UNITS', name

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        assert name != 'UNITS', name

        #fields = []
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
        dict_fields = {irecord: list_fields,}
        return DTI(name, dict_fields, comment=comment)

    def raw_fields(self):
        list_fields = []
        for irecord, fields in sorted(self.fields.items()):
            nfields = len(fields)
            list_fields += ['DTI', self.name] + fields
            print('dti-tree', fields)

            nleftover = nfields % 8
            if nleftover:
                list_fields += [None] * nleftover
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        msg = self.comment
        for irecord, fields in sorted(self.fields.items()):
            assert isinstance(fields, list), fields
            list_fields = ['DTI', self.name, irecord, ] #+ fields
            for fieldsi in fields:
                if isinstance(fieldsi, (list, tuple, np.ndarray)):
                    list_fields += list(fieldsi)
                else:
                    list_fields.append(fieldsi)
            msg += print_card_8(list_fields)
        return msg

    def __repr__(self) -> str:
        """
        Prints a card in the simplest way possible
        (default values are left blank).

        """
        comment = self.comment
        try:
            return self.write_card(size=8)
        except Exception:
            try:
                return self.write_card(size=16)
            except Exception:
                print('problem printing %s card' % self.type)
                #print("list_fields = ", list_fields)
                raise


class NastranMatrix(BaseCard):
    """
    Base class for the DMIG, DMIJ, DMIJI, DMIK matrices
    """
    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        self.finalize()

    def __init__(self, name: str, matrix_form: int,
                 tin: int, tout: int, polar: int, ncols: int,
                 GCj: list[tuple[int, int]],
                 GCi: list[tuple[int, int]],
                 Real: list[float], Complex=None, comment: str='', finalize: bool=True):
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
        GCj  : list[(node, dof)]
            the jnode, jDOFs
        GCi  : list[(node, dof)]
            the inode, iDOFs
        Real : list[float]
            The real values
        Complex : list[float]; default=None
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

        if matrix_form not in {1, 2, 4, 5, 6, 8, 9}:
            msg = (
                f'matrix_form={matrix_form!r} must be [1, 2, 4, 5, 6, 8, 9]\n'
                '  1: Square\n'
                '  2: Rectangular\n'
                #'  4: Lower Triangular\n'
                #'  5: Upper Triangular\n'
                '  6: Symmetric\n'
                #'  8: Identity (m=nRows, n=m)\n'
                '  9: Rectangular\n')
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
            assert self.tin in [3, 4], f'tin={self.tin!r} and must 3 or 4 to be complex'
            assert self.tout in [0, 3, 4], f'tin={self.tout!r} and must 0, 3 or 4 to be complex'
        assert isinstance(matrix_form, integer_types), 'matrix_form=%r type=%s' % (matrix_form, type(matrix_form))
        assert not isinstance(matrix_form, bool), 'matrix_form=%r type=%s' % (matrix_form, type(matrix_form))
        if finalize:
            self.finalize()

    def __deepcopy__(self, memo):
        """doesn't copy the label_actors to speed things up?"""
        #keys = ['name', '_color', 'display', 'line_width', 'point_size', '_opacity',
                #'_representation', 'is_visible', 'bar_scale', 'is_pickable']
        cls = self.__class__
        result = cls.__new__(cls)
        idi = id(self)
        memo[idi] = result
        for key in list(self.__dict__.keys()):
            value = self.__dict__[key]
            setattr(result, key, deepcopy(value, memo))
        return result

    @property
    def matrix_form_str(self) -> str:
        return DMI_MATRIX_MAP[self.matrix_form]

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

        matrix_form = integer(card, 3, 'ifo/matrix_form')
        tin = integer(card, 4, 'tin')
        tout = integer_or_blank(card, 5, 'tout', default=0)
        polar = integer_or_blank(card, 6, 'polar', default=0)
        if matrix_form == 1: # square
            ncols = integer_or_blank(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        elif matrix_form == 6: # symmetric
            ncols = integer_or_blank(card, 8, 'matrix_form=%s; ncol' % matrix_form)
        elif matrix_form in {2, 9}: # rectangular
            # If NCOL is not used for rectangular matrices:
            # - IFO=9, GJ and CJ will determine the sorted sequence, but will
            #   otherwise be ignored; a rectangular matrix will be generated
            #   with the columns submitted being in the 1 to N positions, where
            #   N is the number of logical entries submitted (not counting
            #   the header entry).
            # - IFO=2, the number of columns of the rectangular matrix will be
            #   equal to the index of the highest numbered non-null column (in
            #   internal sort). Trailing null columns of the g- or p-size matrix
            #   will be truncated.
            ncols = integer_or_blank(card, 8, f'matrix_form={matrix_form}; ncol')
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'matrix_form=%s; ncol' % self.matrix_form)

            msg = (
                f'{cls.type} name={name!r} matrix_form={matrix_form!r} is not supported.  Valid forms:\n'
                '  4=Lower Triangular\n'
                '  5=Upper Triangular\n'
                '  6=Symmetric\n'
                '  8=Identity (m=nRows, n=m)\n'
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
        elif self.matrix_form in {2, 9}:
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
        if self.matrix_form in {1, 6}: # square, symmetric
            if self.ncols is not None:
                shape = (self.ncols, self.ncols)
            else:
                nrows, ncols = get_row_col_map(
                    self, self.GCi, self.GCj, self.matrix_form)[:2]
                shape = (nrows, ncols)
        elif self.matrix_form in {2, 9}:
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

    def get_matrix(self, is_sparse: bool=False,
                   apply_symmetry: bool=True) -> tuple[np.ndarray | scipy.coomatrix,
                                                       dict[int, tuple[int, int]],
                                                       dict[int, tuple[int, int]],
                                                       ]:
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
    def is_real(self) -> bool:
        """real vs. complex attribute"""
        return not self.is_complex

    @property
    def is_complex(self) -> bool:
        """real vs. complex attribute"""
        if self.tin in {1, 2}: # real
            return False
        elif self.tin in {3, 4}: # complex
            return True
        msg = ('Matrix %r must have a value of TIN = [1, 2, 3, 4].\n'
               'TIN defines the type (real, complex) '
               'of the matrix.  TIN=%r.\n'
               '  TIN=1,2 -> real\n'
               '  TIN=3,4 -> complex' % (self.name, self.tin))
        raise ValueError(msg)

    @property
    def is_polar(self) -> bool:
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
        msg = (f'Matrix {self.name!r} must have a value of POLAR = [0, 1].\n'
               'POLAR defines the type (real/imag or mag/phase) complex) '
               f'of the matrix.  POLAR={self.polar!r}.')
        raise ValueError(msg)

    @property
    def tin_dtype(self) -> str:
        """gets the input dtype"""
        return _get_dtype(self.is_complex, self.tin)

    @property
    def tout_dtype(self) -> str:
        """gets the output dtype"""
        return _get_dtype(self.is_complex, self.tout)

    def __repr__(self) -> str:
        return self.write_card(size=8, is_double=False)

    def fill_in_default_components(self, model: BDF) -> None:
        for i, (Gi, Ci) in enumerate(self.GCi):
            if Ci is None:
                node = model.nodes[Gi]
                if node.type == 'GRID':
                    msg = ('Ci on DMIG card must be 1, 2, 3, 4, 5, or 6; '
                           'Node=%i (GRID); Ci=%s' % (Gi, Ci))
                    raise RuntimeError(msg)
                elif node.type in {'SPOINT', 'EPOINT'}:
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
                elif node.type in {'SPOINT', 'EPOINT'}:
                    Cj = 0
                else:
                    raise NotImplementedError(node)
                self.GCj[i] = [Gj, Cj]
        return

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        size, is_double = _determine_size_double_from_tin(
            self.tin, size, is_double)

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

        if size == 8 and len(self.GCi):
            Gi = np.array(self.GCi)[:, 0]
            Gj = np.array(self.GCj)[:, 0]
            if max(Gi.max(), Gj.max()) >= 100000000:
                size = 16
            del Gi, Gj

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

def _determine_size_double_from_tin(tin: int,
                                    size: int, is_double: bool) -> tuple[int, bool]:
    """
    we ignore the requested is_double flag because otherwise Nastran
    can't read in the matrix
    """
    if tin in {1, 3}:
        is_double = False
    elif tin in {2, 4}:
        is_double = True
        size = 16
    else:
        raise RuntimeError('tin=%r must be 1, 2, 3, or 4' % tin)
    return size, is_double

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
        if self.tin in {1, 3}:
            is_double = False
            msg = self.write_card_8()
        elif self.tin in {2, 4}:
            is_double = True
            size = 16
            msg = self.write_card_16()
        else:
            raise RuntimeError('tin=%r must be 1, 2, 3, or 4' % self.tin)
        return msg

    def write_card_8(self) -> str:
        """writes the card in small field format"""
        return self._write_card(print_card_8)

    def write_card_16(self) -> str:
        """writes the card in small large format"""
        return self._write_card(print_card_16)

    def _write_card(self, func) -> str:
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
            2/9=Rectangular
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
        GCj  : list[(node, dof)]
            the [jnode, jDOFs] columns
        GCi  : list[(node, dof)]
            the [inode, iDOFs] rows
        Real : list[float]
            The real values
        Complex : list[float]; default=None
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
        GCNj  : list[(node, dof, harmonic_number)]???
            the jnode, jDOFs
        GCNi  : list[(node, dof, harmonic_number)]???
            the inode, iDOFs
        Real : list[float]???
            The real values
        Complex : list[float]???; default=None
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
    def is_real(self) -> bool:
        """is the matrix real?"""
        if self.tin in [1, 2]:
            return True
        return False

    @property
    def is_complex(self) -> bool:
        """is the matrix complex"""
        return not self.is_real

    @property
    def is_polar(self) -> bool:
        """is the matrix polar (vs real/imag)?"""
        return False

    @property
    def tin_dtype(self) -> str:
        """gets the input dtype"""
        return _get_dtype(self.is_complex, self.tin)

    @property
    def tout_dtype(self) -> str:
        """gets the output dtype"""
        return _get_dtype(self.is_complex, self.tout)

    @property
    def matrix_type(self) -> str:
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
        #elif self.matrix_form in {2, 9}:
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
        elif matrix_form in {2, 9}: # rectangular
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
        if self.tin in {1, 3} and size == 8:
            is_double = False
            msg = self.write_card_8()
        elif self.tin in {1, 3} and size == 16:
            is_double = False
            msg = self.write_card_16()
        elif self.tin in [2, 4]:
            is_double = True
            size = 16
            msg = self.write_card_double()
        else:
            raise RuntimeError('tin=%r must be 1, 2, 3, or 4' % self.tin)
        return msg

    def write_card_8(self):
        """writes the card in small field format"""
        return self._write_card(print_card_8)

    def write_card_16(self):
        """writes the card in small large format"""
        return self._write_card(print_card_16)

    def write_card_double(self):
        """writes the card in small large format"""
        return self._write_card(print_card_double)

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
        GCj  : list[(node, dof)]???
            the jnode, jDOFs
        GCi  : list[(node, dof)]???
            the inode, iDOFs
        Real : list[float]???
            The real values
        Complex : list[float]???; default=None
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

    def __init__(self, name: str, ifo: int,
                 tin: int, tout: int, polar: int,
                 ncols: int,
                 GCj, GCi, Real, Complex=None,
                 comment: str='', finalize: bool=True):
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
        GCj  : list[(node, dof)]???
            the jnode, jDOFs
        GCi  : list[(node, dof)]???
            the inode, iDOFs
        Real : list[float]???
            The real values
        Complex : list[float]???; default=None
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
        GCj  : list[(node, dof)]
            the jnode, jDOFs
        GCi  : list[(node, dof)]
            the inode, iDOFs
        Real : list[float]
            The real values
        Complex : list[float]; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card

        """
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment,
                               finalize=finalize)

DMI_MATRIX_MAP = {
    1: 'square',
    2: 'rectangular',  # 9 ???
    3: 'diagonal',
    6: 'symmetric',
    9: 'identity',
}
REVERSE_DMI_MAP = {value: key for key, value in DMI_MATRIX_MAP.items()}

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

    def __init__(self, name: str, matrix_form: int, tin: int, tout: int,
                 nrows: int, ncols: int,
                 GCj, GCi, Real, Complex=None, comment='', finalize=True):
        """

        Parameters
        ----------
        name : str
            The name of the matrix
        matrix_form : str | int
            The shape of the matrix
        tin / tout : int
            The matrix input precision
        nrows / ncols : int
            The number of rows and columns
        GCj : list[int]
            list of column ids
        GCi : list[int]
            list of row ids
        Real : list[float]
            The real values
        Complex : Optional[list[float]]
            The complex values (if the matrix is complex)
        comment : str; default=''
            The comment for the card
        finalize : bool; default=True
            Finish creating the card (set to True if user)

        """
        #NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               #GCj, GCi, Real, Complex, comment='')
        if comment:
            self.comment = comment

        if Complex is None:
            Complex = []

        #-------------------------------------------------------------------------------------
        if isinstance(tin, str):
            reverse_tout_map = {value: key for key, value in TOUT_DTYPE_MAP.items()}
            tin2 = tin.lower().strip()
            try:
                tin = reverse_tout_map[tin2]
            except:
                keys = list(TOUT_DTYPE_MAP) + list(reverse_tout_map)
                raise SyntaxError(f'tin={tin!r} is not in allowed={keys}')

        if tout is None:
            tout = 0
        if isinstance(tout, str):
            reverse_tout_map = {value: key for key, value in TOUT_DTYPE_MAP.items()}
            tout2 = tout.lower().strip()
            try:
                tout = reverse_tout_map[tout2]
            except:
                keys = list(TOUT_DTYPE_MAP) + list(reverse_tout_map)
                raise SyntaxError(f'tout={tout!r} is not in allowed={keys}')

        if isinstance(matrix_form, str):
            matrix_form2 = matrix_form.lower().strip()
            try:
                matrix_form = REVERSE_DMI_MAP[matrix_form2]
            except KeyError:
                keys = list(DMI_MATRIX_MAP) + list(REVERSE_DMI_MAP)
                raise SyntaxError(f'matrix_form={matrix_form!r} is not in allowed={keys}')

        #-------------------------------------------------------------------------------------

        if tout not in {0, 1, 2, 3, 4}:
            raise SyntaxError(f'tout={tout!r} is not in allowed; [1, 2, 3, 4]')

        if matrix_form not in {1, 2, 3, 4, 5, 6, 8}:
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

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        tout = integer_or_blank(card, 5, 'tout', default=0)
        #blank(card, 6, 'blank')
        nrows = integer(card, 7, 'nrows')
        ncols = integer(card, 8, 'ncols')

        assert len(card) == 9, f'len(DMI card) = {len(card):d}\ncard={card}'

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

        try:
            matrix_type = DMI_MATRIX_MAP[self.matrix_form]
        except KeyError:
            raise NotImplementedError('%s matrix_form=%r is not supported' % (
                self.type, self.matrix_form))
        return matrix_type

    @property
    def is_polar(self) -> bool:
        if self.tin in {1, 2}:
            is_polar = False
        elif self.tin in {3, 4}:
            is_polar = False # TODO: could be wrong...
        else:
            raise NotImplementedError(f'nrows={self.nrows} ncols={self.ncols}; tin={self.tin} not [1,2,3,4]')
        return is_polar

    @property
    def shape(self) -> tuple[int, int]:
        return self.nrows, self.ncols

    @property
    def ifo(self) -> int:
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

    def _add_column(self, card: BDFCard, comment: str=''):
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
                        assert real_value == 'THRU', real_value
                        real_value = self.Real[-1]
                        end_i = fields[i + 1]
                        for ii in range(i1, end_i + 1):
                            #print('THRU adding j=%s i1=%s val=%s' % (j, ii, real_value))
                            self.GCj.append(j)
                            self.GCi.append(ii)
                            self.Real.append(real_value)
                        i += 1
                        is_done_reading_floats = True
        #print(self.GCi)
        #print(self.GCj)
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
    def is_real(self) -> bool:
        """real vs. complex attribute"""
        return not self.is_complex

    @property
    def is_complex(self) -> bool:
        """real vs. complex attribute"""
        if self.tin in {3, 4}:
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
        msg = _get_real_matrix_columns(
            self.name, self.GCi, self.GCj, self.Real, func)
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

            #if reals.max() == 0. and reals.min() == 0. and complexs.max() == 0. and complexs.min() == 0.:
                #continue

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

    def get_matrix(self,
                   is_sparse: bool=False,
                   apply_symmetry: bool=True) -> tuple[np.array, None, None]:
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
        mat, rows, cols = get_dmi_matrix(
            self, is_sparse=is_sparse, apply_symmetry=apply_symmetry)
        return mat, rows, cols

    def write_card_16(self):
        """writes the card in single precision"""
        return self._write_card(print_card_16)

    def write_card_double(self):
        """writes the card in double precision"""
        return self._write_card(print_card_double)

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
        size, is_double = _determine_size_double_from_tin(
            self.tin, size, is_double)

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


def _get_real_matrix_columns(name: str, GCi, GCj, Real,
                             func: Callable[float, str]) -> str:
    msg = ''
    uGCj = np.unique(GCj)
    #print(f'uGCj={uGCj}')
    for gcj in uGCj:
        # get a single column and filter 0 values
        i = np.where((gcj == GCj) & (Real != 0.))[0]
        if len(i) == 0:
            list_fields = ['DMI', name, gcj, 1, 0.0]
            msg += func(list_fields)
            continue
        assert len(i) > 0, i
        singles, doubles = collapse_thru_ipacks(i, GCi[i].tolist())
        assert len(singles) + len(doubles) > 0, (singles, doubles)
        # if singles:
        #     gcis = GCi[singles]
        #     reals = Real[singles]
        #     list_fields = ['DMI', name, gcj]
        #     for gci, real in zip(gcis, reals):
        #         list_fields.extend([gci, real])
        #     msg += func(list_fields)

        list_fields = ['DMI', name, gcj]
        for (start, thru, end) in doubles:
            i2 = slice(start, end+1)
            gcis2 = GCi[i2]
            reals2 = Real[i2]
            if reals2.max() == reals2.min():
                gci1 = gcis2[0]
                gci2 = gcis2[-1]
                real = reals2[0]
                list_fields.extend([gci1, real, 'THRU', gci2])
            else:
                for gci, real in zip(gcis2, reals2):
                    list_fields.extend([gci, real])
                #list_fields.extend([gci, 'THRU', real])
                #print(double)
            assert len(list_fields) > 3, list_fields
        #print(list_fields)
        if singles:
            gcis1 = GCi[singles]
            reals1 = Real[singles]
            for gci, real in zip(gcis1, reals1):
                list_fields.extend([gci, real])
        assert len(list_fields) > 3, list_fields
        msg += func(list_fields)
    return msg

def _get_real_matrix_columns2(name: str, GCi, GCj, Real,
                              func: Callable[float, str]):  # pramga: no cover
    msg = ''
    uGCj = np.unique(GCj)
    #print(f'uGCj={uGCj}')
    for gcj in uGCj:
        i = np.where((gcj == GCj) & (Real != 0.))[0]
        gcis = GCi[i].copy()
        reals = Real[i].copy()
        isort = np.argsort(gcis)
        gcis = gcis[isort]
        reals = reals[isort]

        gc1 = gcis[0]
        nreals = len(reals)
        # print(f'gcj={gcj} gcis={gcis}')
        # print(f'reals={reals}')
        list_fields = ['DMI', name, gcj]
        max_value = reals.max()

        if nreals == 1:
            # 1 value
            list_fields.extend([gc1, max_value])
            msg += func(list_fields)
            return
        # print(f'gcis = {gcis}')
        assert len(gcis) == len(np.unique(gcis)), gcis
        gc_end = gcis[-1]
        gc_range = np.arange(gc1, gc_end + 1)

        # dgci = gc_end - gc1 #+ 1
        # print(f'dgci={dgci}')
        # print(f'gc_range={gc_range}')
        if len(gc_range) == len(i):
            # dense
            # print(f'dense0, dgci={dgci}')
            if reals.max() == reals.min():
                # DMI     WKK     1       1       1.0     THRU    112
                list_fields.extend([gc1, max_value, 'THRU', gc_end])
                msg += func(list_fields)
                return
            else:
                list_fields.append(gc1)
                list_fields.extend(reals)
            msg += func(list_fields)
            return

        real_range = np.zeros(len(gc_range), dtype=reals.dtype)
        igc_range = np.searchsorted(gc_range, gcis)
        real_range[igc_range] = reals
        gcis = gc_range
        reals = real_range
        isort = np.argsort(gcis)

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


def get_row_col_map(matrix: DMIG,
                    GCi: np.ndarray, GCj: np.ndarray,
                    ifo: int) -> tuple[int, int, int,
                                       np.ndarray, np.ndarray,
                                       dict[int, Any],
                                       dict[int, Any]]:
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

def _get_row_col_map_1d(matrix, GCi, GCj, ifo: int):
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

    i = 0
    for gc in GCi:
        tgc = tuple(gc)
        if tgc not in rows:
            rows[tgc] = i
            rows_reversed[i] = tgc
            i += 1

    if ifo == 6:
        # symmetric
        for gc in GCj:
            tgc = tuple(gc)
            if tgc not in rows:
                rows[tgc] = i
                rows_reversed[i] = tgc
                i += 1
        cols = rows
        cols_reversed = rows_reversed
    else:
        j = 0
        for gc in GCj:
            tgc = tuple(gc)
            if tgc not in cols:
                cols[tgc] = j
                cols_reversed[j] = tgc
                j += 1
    return rows, cols, rows_reversed, cols_reversed

def _fill_sparse_matrix(matrix: DMIG, nrows: int, ncols: int,
                        apply_symmetry: bool) -> coo_matrix:
    """helper method for ``get_matrix``"""
    if matrix.GCi.ndim == 1:
        assert matrix.GCj.ndim == 1, matrix.GCj.ndim
        rows = matrix.GCi
        cols = matrix.GCj
        GCj = np.array(matrix.GCj, dtype='int32') - 1
        GCi = np.array(matrix.GCi, dtype='int32') - 1
        # TODO: matrix size:  is this correct?
        nrows = max(GCi) + 1
        ncols = max(GCj) + 1
    else:
        assert matrix.GCi.ndim == 2, matrix.GCi.ndim
        assert matrix.GCj.ndim == 2, matrix.GCj.ndim
        GCi = matrix.GCi
        GCj = matrix.GCj

        # TODO: matrix size:  is this correct?
        #Gi = GCi[:, 0]
        #Gj = GCj[:, 0]
        #GCij = np.vstack([GCi, GCj])
        #uGCij, idx = unique2d(GCi, return_index=True)
        if matrix.matrix_form == 6:  # symmetric
            ngc = GCi.shape[0]
            GCij = np.vstack([GCi, GCj])
            uGCij, idx = unique2d(GCij, return_index=True)
            rows_cols, nrows = gc_to_index(GCij)
            ncols = nrows

            rows = rows_cols[:ngc]
            cols = rows_cols[ngc:]
        else:
            # symmetric matrices will be wrong if the values are unordered...
            rows, nrows = gc_to_index(GCi)
            cols, ncols = gc_to_index(GCj)
        #rows = unique2d(GCi)
        #cols = unique2d(GCj)
        #nrows = unique2d(rows).shape[0]
        #ncols = unique2d(cols).shape[0]

    float_dtype = _get_real_dtype(matrix.tin)
    reals = np.array(matrix.Real, dtype=float_dtype)

    dtype = _get_dtype(matrix.is_complex, matrix.tin)

    if matrix.is_complex:
        complexs = np.array(matrix.Complex, dtype=float_dtype)
        data = reals + 1j * complexs
    else:
        data = reals

    if matrix.matrix_form in {1, 6}:
        nrows = max(nrows, ncols)
        ncols = nrows

    assert len(rows) == len(cols)
    assert len(data) == len(rows)
    if matrix.matrix_form == 6 and apply_symmetry:
        is_diagonal, not_diagonal = _get_diagonal_symmetric(matrix)
        if np.any(not_diagonal):
            rows2 = np.hstack([rows, cols[not_diagonal]])
            cols2 = np.hstack([cols, rows[not_diagonal]])
            data = np.hstack([data, data[not_diagonal]])
            assert len(rows2) == len(cols2)
            assert len(data) == len(rows2)
            rows = rows2
            cols = cols2

    sparse_matrix = coo_matrix(
        (data, (rows, cols)),
        shape=(nrows, ncols), dtype=dtype)
    return sparse_matrix

def _build_gc_map(GC: np.ndarray) -> dict[tuple[int, int], int]:
    """helper method for ``gc_to_index``"""
    i = 0
    gc_map = {}
    for gc in GC:
        tgc = tuple(gc)
        if tgc in gc_map:
            continue
        gc_map[tgc] = i
        i += 1
    return gc_map

def gc_to_index(GC: np.ndarray) -> tuple[np.ndarray, int]:
    """helper method for ``_fill_sparse_matrix``"""
    gc_map = _build_gc_map(GC)
    ngrid_map = len(gc_map)
    ngrid = GC.shape[0]
    index = np.zeros(ngrid, dtype='int32')
    for i, gc in enumerate(GC):
        tgc = tuple(gc)
        j = gc_map[tgc]
        index[i] = j
    return index, ngrid_map

def _fill_dense_rectangular_matrix(matrix: DMIG,
                                   nrows: int, ncols: int, ndim: int,
                                   rows: dict[Any, int], cols: dict[Any, int],
                                   apply_symmetry: bool) -> Any:
    """helper method for ``get_matrix``"""
    if matrix.is_complex:
        dense_mat = _fill_dense_rectangular_matrix_complex(
            matrix, nrows, ncols, ndim, rows, cols, apply_symmetry)
    else:
        dense_mat = _fill_dense_rectangular_matrix_real(
            matrix, nrows, ncols, ndim, rows, cols, apply_symmetry)
        assert isinstance(dense_mat, np.ndarray), type(dense_mat)
    return dense_mat

def _fill_dense_rectangular_matrix_complex(matrix: DMIG,
                                           nrows: int, ncols: int, ndim: int,
                                           rows: dict[Any, int], cols: dict[Any, int],
                                           apply_symmetry: bool) -> np.ndarray:
    """helper method for ``_fill_dense_rectangular_matrix``"""
    dense_mat = np.zeros((nrows, ncols), dtype=matrix.tin_dtype)
    real_imag = matrix.Real + 1j * matrix.Complex
    if matrix.matrix_form == 6 and apply_symmetry:  # symmetric
        is_diagonal, not_diagonal = _get_diagonal_symmetric(matrix)
        for (gcj, real_imagi) in zip(matrix.GCj[is_diagonal], real_imag[is_diagonal]):
            j = cols[(gcj[0], gcj[1])]
            dense_mat[j, j] += real_imagi

        for (gcj, gci, real_imagi) in zip(matrix.GCj[not_diagonal], matrix.GCi[not_diagonal],
                                          real_imag[not_diagonal]):
            i = rows[(gci[0], gci[1])]
            j = cols[(gcj[0], gcj[1])]
            dense_mat[i, j] += real_imagi
            dense_mat[j, i] += real_imagi
    else:
        for (gcj, gci, real_imagi) in zip(matrix.GCj, matrix.GCi, real_imag):
            i = rows[(gci[0], gci[1])]
            j = cols[(gcj[0], gcj[1])]
            dense_mat[i, j] += real_imagi
    return dense_mat

def _get_diagonal_symmetric(matrix: DMIG) -> tuple[np.ndarray, np.ndarray]:
    """helper for ``apply_symmetry``"""
    assert matrix.GCi.ndim == 2, matrix.GCi.ndim
    assert matrix.GCj.ndim == 2, matrix.GCj.ndim
    dij = matrix.GCi - matrix.GCj
    dij[:, 0] == dij[:, 1]
    is_diagonal = (dij[:, 0] == 0) & (dij[:, 1] == 0)
    not_diagonal = ~is_diagonal
    return is_diagonal, not_diagonal

def _fill_dense_rectangular_matrix_real(matrix: DMIG,
                                        nrows: int, ncols: int, ndim: int,
                                        rows: dict[Any, int], cols: dict[Any, int],
                                        apply_symmetry: bool) -> np.ndarray:
    """helper method for ``_fill_dense_rectangular_matrix``"""
    dense_mat = np.zeros((nrows, ncols), dtype=matrix.tin_dtype)
    if matrix.matrix_form == 6 and apply_symmetry:  # symmetric
        is_diagonal, not_diagonal = _get_diagonal_symmetric(matrix)
        try:
            for (gcj, reali) in zip(matrix.GCj[is_diagonal], matrix.Real[is_diagonal]):
                i = rows[tuple(gcj)]
                dense_mat[i, i] += reali

            for (gcj, gci, reali) in zip(matrix.GCj[not_diagonal], matrix.GCi[not_diagonal], matrix.Real[not_diagonal]):
                i = rows[tuple(gci)]
                j = cols[tuple(gcj)]
                dense_mat[i, j] += reali
                dense_mat[j, i] += reali
        except IndexError:
            msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                   'is_polar=%s ncols=%s M.shape=%s\n' % (
                       matrix.name, ndim, i, j, matrix.matrix_type,
                       matrix.is_polar, matrix.ncols, dense_mat.shape))
            msg += 'Rows:\n'
            for i, row in enumerate(rows):
                msg += 'i=%s row=%s\n' % (i, row)
            raise RuntimeError(msg)
    else:
        try:
            for (gcj, gci, reali) in zip(matrix.GCj, matrix.GCi, matrix.Real):
                tgci = tuple(gci)
                tgcj = tuple(gcj)
                i = rows[tgci]
                j = cols[tgcj]
                dense_mat[i, j] += reali
        except KeyError:
            msg = ('name=%s ndim=%s gci=%s gcj=%s matrix_type=%s '
                   'is_polar=%s ncols=%s M.shape=%s\n\n' % (
                       matrix.name, ndim, str(gci), str(gcj), matrix.matrix_type,
                       matrix.is_polar, matrix.ncols, dense_mat.shape))

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
                   'is_polar=%s ncols=%s M.shape=%s\n' % (
                       matrix.name, ndim, i, j, matrix.matrix_type,
                       matrix.is_polar, matrix.ncols, dense_mat.shape))
            msg += 'Rows:\n'
            for i, row in enumerate(rows):
                msg += '  i=%s row=%s\n' % (i, row)

            msg += '\nCols:\n'
            for j, row in enumerate(cols):
                msg += '  j=%s row=%s\n' % (j, row)
            raise RuntimeError(msg)
    return dense_mat


def _fill_dense_column_matrix(matrix: DMIG,
                              nrows: int, ncols: int, ndim: int,
                              rows: dict[Any, int], cols: dict[Any, int],
                              apply_symmetry: bool) -> np.ndarray:
    """helper method for ``get_matrix``"""
    if matrix.is_complex:
        dense_mat = _fill_dense_column_matrix_complex(
            matrix, nrows, ncols, ndim, rows, cols, apply_symmetry)
    else:
        dense_mat = _fill_dense_column_matrix_real(
            matrix, nrows, ncols, ndim, rows, cols, apply_symmetry)
    return dense_mat

def _fill_dense_column_matrix_real(matrix: DMIG,
                                   nrows: int, ncols: int, ndim: int,
                                   rows: dict[Any, int], cols: dict[Any, int],
                                   apply_symmetry: bool) -> np.ndarray:
    """helper method for ``_fill_dense_column_matrix``

    What does symmetry mean for a column matrix?!!!
    """
    #print('nrows=%s ncols=%s' % (nrows, ncols))
    dense_mat = np.zeros((nrows, ncols), dtype=matrix.tin_dtype)
    if matrix.matrix_form == 6 and apply_symmetry:  # symmetric
        assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
        raise RuntimeError('What does symmetry mean for a column matrix?!!!')
        #for (gcj, gci, reali) in zip(matrix.GCj, matrix.GCi, matrix.Real):
            #i = rows[gci]
            #j = cols[gcj]
            #dense_mat[i, j] += reali
            #dense_mat[j, i] += reali
    else:
        try:
            for (gcj, gci, reali) in zip(matrix.GCj, matrix.GCi, matrix.Real):
                i = rows[gci]
                j = cols[gcj]
                dense_mat[i, j] += reali
        except IndexError:
            msg = ('name=%s ndim=%s gci=%s gcj=%s matrix_type=%s '
                   'is_polar=%s ncols=%s M.shape=%s\n' % (
                       matrix.name, ndim, gci, gcj, matrix.matrix_type,
                       matrix.is_polar, matrix.ncols, dense_mat.shape))
            msg += 'Rows:\n'
            for i, row in enumerate(rows):
                msg += '  i=%s row=%s\n' % (i, row)
            raise RuntimeError(msg)
    return dense_mat

def _fill_dense_column_matrix_complex(matrix: DMIG,
                                      nrows: int, ncols: int, ndim: int,
                                      rows: dict[Any, int], cols: dict[Any, int],
                                      apply_symmetry: bool) -> np.ndarray:
    """
    helper method for ``_fill_dense_column_matrix``

    What does symmetry mean for a column matrix?!!!
    """
    dense_mat = np.zeros((nrows, ncols), dtype=matrix.tin_dtype)
    if matrix.matrix_form == 6 and apply_symmetry:  # symmetric
        assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
        raise RuntimeError('What does symmetry mean for a column matrix?!!!')
        #for (gcj, gci, reali, complexi) in zip(matrix.GCj, matrix.GCi,
                                               #matrix.Real, matrix.Complex):
            #i = rows[gci]
            #j = cols[gcj]
            #dense_mat[i, j] += complex(reali, complexi)
            #dense_mat[j, i] += complex(reali, complexi)
    elif matrix.matrix_form == 2:  # rectangular
        assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
        for (gcj, gci, reali, complexi) in zip(matrix.GCj, matrix.GCi,
                                               matrix.Real, matrix.Complex):
            i = rows[gci]
            j = cols[gcj]
            dense_mat[i, j] += complex(reali, complexi)
    else:
        for (gcj, gci, reali, complexi) in zip(matrix.GCj, matrix.GCi,
                                               matrix.Real, matrix.Complex):
            i = rows[gci]
            j = cols[gcj]
            dense_mat[i, j] += complex(reali, complexi)
    return dense_mat


def get_dmi_matrix(matrix: DMI,
                   is_sparse: bool=False,
                   apply_symmetry: bool=True) -> tuple[np.array, None, None]:
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
    if isinstance(matrix.GCi, np.ndarray):
        assert matrix.GCi.ndim == 1, matrix.GCi.ndim
        assert matrix.GCj.ndim == 1, matrix.GCj.ndim
    else:
        # TestAero.test_zona_2
        warnings.warn(f'matrix={matrix.name!r} GCi is not a numpy array...type={type(matrix.GCi)}')
        #print(matrix)
        #print('GCi =', matrix.GCi)
        #print('GCj =', matrix.GCj)
    GCj = np.array(matrix.GCj, dtype='int32') - 1
    GCi = np.array(matrix.GCi, dtype='int32') - 1

    dtype = matrix.tin_dtype

    if matrix.is_complex:
        data = matrix.Real + matrix.Complex * 1j
    else:
        data = matrix.Real

    matrix_form_str = matrix.matrix_form_str
    if matrix_form_str in {'rectangular', 'identity'}:  # 2, 9
        # rectangular
        nrows = matrix.nrows
        ncols = matrix.ncols
        M = _set_matrix(nrows, ncols,
                        data, GCi, GCj,
                        dtype)
    elif matrix_form_str == 'diagonal':
        nrows = max(matrix.nrows, matrix.ncols)
        assert matrix.ncols == 1, (matrix.nrows, matrix.ncols)
        GCj = np.zeros(len(GCi), dtype=GCi.dtype)
        ncols = 1
        M = _set_matrix(nrows, ncols,
                        data, GCi, GCj,
                        dtype)
    elif matrix_form_str == 'square':
        nrows = matrix.nrows
        ncols = matrix.ncols
        assert nrows == ncols, (nrows, ncols)
        M = _set_matrix(nrows, ncols,
                        data, GCi, GCj,
                        dtype)
    else:
        raise RuntimeError(matrix_form_str)
        nrows = matrix.nrows
        ncols = matrix.ncols
        if matrix_form_str == 'symmetric':
            nrows = max(nrows, ncols)
            ncols = nrows
            #matrix_form_str = 'square'
        M = _set_matrix(nrows, ncols,
                        data, GCi, GCj,
                        dtype)

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

def get_matrix(self: DMIG,
               is_sparse: bool=False,
               apply_symmetry: bool=False) -> tuple[Any,
                                                   dict[int, Any],
                                                   dict[int, Any]]:
    """
    Builds the Matrix

    Parameters
    ----------
    is_sparse : bool; default=False
        should the matrix be returned as a sparse matrix.
        Slower for dense matrices.
    apply_symmetry: bool; default=False
        If the matrix is symmetric (matrix_form=6), returns a symmetric matrix.
        Supported as there are symmetric matrix routines.

    Returns
    -------
    M : ndarray
        the matrix
    rows : dict[(nid, nid)] = float
        dictionary of keys=rowID,    values=(Grid,Component) for the matrix
    cols : dict[(int, int)] = float
        dictionary of keys=columnID, values=(Grid,Component) for the matrix

    """
    nrows, ncols, ndim, rows, cols, rows_reversed, cols_reversed = get_row_col_map(
        self, self.GCi, self.GCj, self.matrix_form)

    #is_sparse = False
    if is_sparse:
        #assert isinstance(self, (DMIG, DMIK, DMIJI)), type(self)
        M = _fill_sparse_matrix(self, nrows, ncols, apply_symmetry)
    else:
        if ndim == 1:
            #assert isinstance(self, int), type(self)
            M = _fill_dense_column_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry)
            assert isinstance(M, np.ndarray), type(M)
        else:
            #assert isinstance(self, (DMIG, DMIK, DMIJ, DMIJI)), type(self)
            M = _fill_dense_rectangular_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry)
            assert isinstance(M, np.ndarray), type(M)
    return M, rows_reversed, cols_reversed


def _set_matrix(nrows: int, ncols: int,
                data: np.ndarray,
                GCi: np.ndarray, GCj: np.ndarray,
                dtype: str) -> coo_matrix:
    try:
        matrixi = coo_matrix((data, (GCi, GCj)),
                             shape=(nrows, ncols), dtype=dtype)
    except ValueError:
        print(f'nrows, cols = ({nrows}, {ncols})')
        print('data = ', data)
        print('GCi = ', GCi)
        print('GCj = ', GCj)
        raise
    return matrixi

def _export_dmig_to_hdf5(h5_file, model: BDF, dict_obj, encoding: str) -> None:
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


def _export_dmiax_to_hdf5(h5_file, model, dict_obj, encoding: str) -> None:
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

def _set_polar(polar) -> int:
    if polar in {None, 0, False}:
        polar = 0
    elif polar in {1, True}:
        polar = 1
    else:  # pragma: no cover
        raise ValueError(f'polar={polar!r} and must be 0 or 1')
    return polar

def _get_dtype(is_complex: bool, type_flag: int) -> str:
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
    else:  # pragma: no cover
        raise RuntimeError(f'invalid option for matrix format {type_flag}')
    return dtype

def _get_real_dtype(type_flag: int) -> str:
    """A complex64 array is made up of two float32 arrays."""
    if type_flag in {1, 3}:
        dtype = 'float32'
    elif type_flag in {0, 2, 4}:
        dtype = 'float64'
    else:  # pragma: no cover
        raise RuntimeError(f'invalid option for matrix format {type_flag}')
    return dtype

def dtype_to_tin_tout_str(myarray: np.ndarray) -> str:
    tin_real = myarray.real.dtype.itemsize
    tin_total = myarray.dtype.itemsize
    if tin_real == 8 and tin_total == 16:
        tin = 'complex128'
    elif tin_real == 8 and tin_total == 8:
        tin = 'float64'
    elif tin_real == 4 and tin_total == 8:
        tin = 'complex64'
    elif tin_real == 4 and tin_total == 4:
        tin = 'float32'
    else:
        raise NotImplementedError('dtype_to_tin_tout')
    return tin
