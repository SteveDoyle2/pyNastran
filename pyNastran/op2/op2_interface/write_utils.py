"""
Defines methods for the op2 & hdf5 writer
"""
from struct import Struct, pack
from typing import BinaryIO, TextIO

import numpy as np
import scipy
import scipy.sparse as sp
#from pyNastran.utils.numpy_utils import integer_float_types
from pyNastran.utils import int_version

SCIPY_VERSION = int_version('scipy', scipy.__version__)[:2]

# address scipy.sparse.coo -> scipy.sparse._coo
IS_NEW_SCIPY = (SCIPY_VERSION >= [1, 8])
IS_OLD_SCIPY = not IS_NEW_SCIPY


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

def export_to_hdf5(self, group, log):
    """exports the object to HDF5 format"""
    #headers = self.get_headers()

    # for some reason we can't just not write the properties...
    names = self.object_attributes(filter_properties=False)
    dynamic_string = [
        'headers', 'data_names', 'words', 'gridtype_str', 'element_data_type', 'location',
        'failure_theory',
    ]

    for name in names:
        if name in ['data_code', 'dataframe', 'data_frame', 'element_mapper', 'h5_file']:
            continue
        value = getattr(self, name)
        if value is None:
            continue
        elif isinstance(value, dict):
            log.warning(f'HDF5: skipping name={name!r} value={value:d}')
            continue
        #elif isinstance(value, (integer_float_types, str, bytes, np.ndarray, list, h5py._hl.dataset.Dataset)):
            #pass
        elif ((IS_NEW_SCIPY and isinstance(value, sp._coo.coo_matrix)) or
              (IS_OLD_SCIPY and isinstance(value, sp.coo.coo_matrix))):
            # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_bsh111svd2.op2
            #
            # https://stackoverflow.com/questions/43390038/storing-scipy-sparse-matrix-as-hdf5
            #g = group.create_group('Mcoo')
            group.create_dataset('data', data=value.data)
            group.create_dataset('row', data=value.row)
            group.create_dataset('col', data=value.col)
            group.attrs['shape'] = value.shape
            continue
        #else:  #pragma, no cover
            #msg = f'type={type(value)} value={value}'
            #raise TypeError(msg)

        #if name in ['dt', 'nonlinear_factor', 'element'] and value is None:
            #continue

        # h5py doesn't support unicode, so we have to turn the data into ASCII.
        # All these are fine, but this routine will probably fail at some point
        # on the subtitle/label being actual unicode.
        if name in ['element_names']:  # grid point forces
            value = np.asarray(value, dtype='|S8').tolist()
        elif name in dynamic_string:
            svalue = [str(valuei) for valuei in value]

            # the size of the array is the |S8 or |S12 or whatever
            max_len = max([(len(valuei)) for valuei in svalue])
            dtype = '|S%i' % max_len
            value = np.array(svalue, dtype=dtype)

        elif name in ['element', 'element_type'] and isinstance(value, np.ndarray):
            if value.dtype is np.dtype(np.int32):
                pass
            else:
                # unicode
                #value = value.tolist()
                value = np.asarray(value, dtype='|S8').tolist()

        #if hasattr(value, 'export_to_hdf5'):
            #msg = 'sub-object export_to_hdf5 not supported\nkey=%s value=%s' % (key, value)
            #raise NotImplementedError(msg)
        if isinstance(value, np.ndarray) and value.dtype.name.startswith('str'):
             # str256
            n = value.dtype.name[3:]
            value = np.asarray(value, dtype='|S'+n)

        try:
            group.create_dataset(name, data=value)
        except TypeError:
            print('name = %r; type=%s' % (name, type(value)))
            print(value)
            print('------------------')
            raise
            #continue
        #print('done')
