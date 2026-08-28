"""
Defines methods for the op2 & hdf5 writer
"""
from struct import Struct, pack
from typing import BinaryIO, TextIO

import numpy as np
#from pyNastran.utils.numpy_utils import integer_float_types


def set_table3_field(str_fields, ifield: int, value):
    """
    ifield is 1 based
    """
    return str_fields[:ifield-1] + value + str_fields[ifield:]

def _write_markers(op2_file: BinaryIO, fascii: TextIO, markers):
    """
    writes pairs of markers

    Parameters
    ----------
    op2_file : file
        the op2 file object
    markers : list[int]
        a set of 3 markers such as [-3, 1, 0] will write as
        [4, -3, 4, 4, 1, 4, 4, 0, 4]
    """
    out = []
    n = 0
    for marker in markers:
        out += [4, marker, 4]
        n += 3
        fascii.write(f'marker = [4, {marker:d}, 4]\n')
    op2_file.write(pack(b'<%ii' % n, *out))


def write_table_header(op2_file: BinaryIO, fascii: TextIO, table_name: str):
    """
    Writes the beginning of an op2 table

    Parameters
    ----------
    op2_file : file
        the op2 file object
    table_name : str
        the table name to write
    """
    table0 = [
        4, 2, 4,
        8, table_name.encode('ascii'), 8,
        #4, 0, 4,
    ]
    assert len(table_name) == 8, table_name
    table0_format = '<4i 8s i'
    struct_table = Struct(table0_format)
    op2_file.write(struct_table.pack(*table0))
    fascii.write('write_table_header: %s header0 = %s\n' % (table_name, table0))


def to_column_bytes(data_list: list[np.ndarray], dtype_out: str,
                    debug: bool=False) -> np.ndarray:
    """
    Takes an stackable numpy array of mixed types (e.g., ints/strings)
    and casts them to the appropriate output datatype
    (typically float32/float64).

    An array is stackable if it's the same shape (e.g., ints/floats).  This
    requirement is a bit looser for strings (4 characters per 32-bit float)
    """
    #shape = data_list[0].shape
    for i, datai in enumerate(data_list):
        #if isinstance(datai, bytes):
            ##print('bytes')
            #data_list[i] = np.frombuffer(datai, dtype=dtype_out)
        if datai.dtype != dtype_out:
            #print(datai.dtype, dtype_out)
            data_list[i] = view_dtype(datai, dtype_out)  # TODO: is this faster/correct?
            #data_list[i] = datai.view(dtype_out)  # TODO: is this faster/correct?
            #data_list[i] = np.frombuffer(datai.tobytes(), dtype=dtype_out)
        elif debug:
            #print('floats...')
            print(datai.shape)
        if debug:
            print(data_list[i].shape)
    try:
        out = np.column_stack(data_list)
    except ValueError:
        for i, datai in enumerate(data_list):
            print(i, datai.shape)
        raise
    return out

def get_complex_fdtype(dtype):
    """complex64 -> float32; complex128 -> float64"""
    if dtype.itemsize == 8:
        return np.float32(1).dtype
    return np.float64(1).dtype # 8

def view_idtype_as_fdtype(int_array: np.ndarray, fdtype: str) -> np.ndarray:
    """
    If we're downcasting from int64 to float32, we can't directly go to float32.
    We need to first go to int32, then to float32.
    """
    if int_array.dtype == np.int64:
        int_array = view_dtype(int_array.astype('int32'), fdtype)
    else:
        #print(f'array_obj.dtype.itemsize={nodedevice_gridtype.dtype.itemsize} dtype.itemsize={fdtype.itemsize}')
        int_array = view_dtype(int_array, fdtype)
    return int_array

def view_dtype(array_obj: np.ndarray, dtype) -> np.ndarray:
    """handles downcasting data"""
    if array_obj.dtype.itemsize == dtype.itemsize:
        return array_obj.view(dtype)
    return array_obj.astype(dtype)


def get_title_subtitle_label(title: str,
                             subtitle: str,
                             label: str,
                             superelement_adaptivity_index: str='',
                             ) -> tuple[bytes, bytes, bytes]:
    """
    TODO: subtitle is missing
      superelement_adaptivity_index
    """
    title_out = b'%-128s' % title.encode('ascii')
    #subtitle_out = b'%-128s' % subtitle.encode('ascii')
    subtitle_out = _write_subtitle_adaptivity_index(
        subtitle, superelement=superelement_adaptivity_index,
        adaptivity_index='')
    label_out = b'%-128s' % label.encode('ascii')

    assert len(title_out) == 128, len(title_out)
    assert len(subtitle_out) == 128, len(subtitle_out)
    assert len(label_out) == 128, len(label_out)
    return title_out, subtitle_out, label_out


def _write_subtitle_adaptivity_index(
        subtitle: str | bytes,
        superelement: int | str='',
        adaptivity_index: int | str='') -> bytes:
    if isinstance(subtitle, str):
        subtitle_prefix = b'%-67s' % subtitle[:67].encode('ascii')
    else:
        assert isinstance(subtitle, bytes), subtitle
        subtitle_prefix = b'%-67s' % subtitle[:67]
    assert len(subtitle_prefix) == 67, (len(subtitle_prefix), subtitle_prefix)

    if superelement:
        # if isinstance(superelement, bytes):
        #     super_adapt_bytes = b'SUPERELEMENT %b' % superelement
        if isinstance(superelement, int):
            super_adapt_bytes = b'SUPERELEMENT %d' % superelement
        else:
            assert isinstance(superelement, str), superelement
            if superelement.startswith('SUPERELEMENT '):
                superelement = superelement[13:]
            super_adapt_bytes = b'SUPERELEMENT %s' % superelement.encode('ascii')
            # title + 'SUPERELEMENT 0, 1'
            # title + 'SUPERELEMENT 0 (id=2000)'
    elif adaptivity_index:
        raise RuntimeError(adaptivity_index)
        super_adapt_bytes = b'ADAPTIVITY INDEX=%7d' % int(adaptivity_index)
    else:
        super_adapt_bytes = b''
    subtitle_out = b'%67s%61s' % (subtitle_prefix, super_adapt_bytes)
    if len(subtitle_out) > 128:
        # subtitle (n=52) = 'TESTING OF ROTOR DYNAMICS MATH AND DMAP CAPABILITIES'
        # subtitle_OG     = 'TESTING OF ROTOR DYNAMICS MATH AND DMAP CAPABILITIES'
        raise RuntimeError(
            f'Too long (n={len(subtitle_out)}\n'
            f'subtitle    = {subtitle!r}\n'
            f'subtitle_prefix   (n={len(subtitle_prefix):d} = {subtitle_prefix!r}\n'
            f'super_adapt_bytes (n={len(super_adapt_bytes):d} = {super_adapt_bytes!r}\n'
            f'out         = {subtitle_out!r}')
    return subtitle_out
