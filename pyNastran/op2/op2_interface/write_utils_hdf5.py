import numpy as np
import scipy.sparse as sp
from pyNastran.utils.scipy_sparse_utils import IS_NEW_SCIPY, IS_OLD_SCIPY


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
