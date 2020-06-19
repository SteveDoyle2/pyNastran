"""
Defines methods for the op2 & hdf5 writer
"""
from struct import Struct, pack
import numpy as np
import scipy.sparse as sp

def set_table3_field(str_fields, ifield, value):
    """
    ifield is 1 based
    """
    return str_fields[:ifield-1] + value + str_fields[ifield:]

def _write_markers(op2_file, fascii, markers):
    """
    writes pairs of markers

    Parameters
    ----------
    op2_file : file
        the op2 file object
    markers : List[int]
        a set of 3 markers such as [-3, 1, 0] will write as
        [4, -3, 4, 4, 1, 4, 4, 0, 4]
    """
    out = []
    n = 0
    for marker in markers:
        out += [4, marker, 4]
        n += 3
        fascii.write('marker = [4, %i, 4]\n' % marker)
    op2_file.write(pack(b'<%ii' % n, *out))


def write_table_header(op2_file, fascii, table_name):
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
    st = Struct(table0_format)
    op2_file.write(st.pack(*table0))
    fascii.write('%s header0 = %s\n' % (table_name, table0))

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
        elif isinstance(value, sp.coo.coo_matrix):
            # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_bsh111svd2.op2
            #
            # https://stackoverflow.com/questions/43390038/storing-scipy-sparse-matrix-as-hdf5
            #g = group.create_group('Mcoo')
            group.create_dataset('data', data=value.data)
            group.create_dataset('row', data=value.row)
            group.create_dataset('col', data=value.col)
            group.attrs['shape'] = value.shape
            continue

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
        try:
            group.create_dataset(name, data=value)
        except TypeError:
            print('name = %r; type=%s' % (name, type(value)))
            print(value)
            print('------------------')
            raise
            #continue
        #print('done')
