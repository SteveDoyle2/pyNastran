from __future__ import annotations
import warnings
import numpy as np
from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.op2.result_objects.utils import (
    real_modes_to_omega_freq, complex_damping_frequency)
from typing import TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd


def build_dataframe_transient_header(self: ScalarObject) -> tuple[list[str], list[np.ndarray]]:
    """builds the header for the Pandas DataFrame/table"""
    assert isinstance(self.name, (str, bytes)), 'name=%s type=%s' % (self.name, type(self.name))
    # print('self.name = %r' % self.name)
    # name = self.name #data_code['name']
    times = self._times
    utimes = np.unique(times)
    if not len(times) == len(utimes):
        msg = 'WARNING : %s - times=%s unique_times=%s...assuming new values...new_times=%s' % (
            self.__class__.__name__, times, utimes, np.arange(len(times)))
        warnings.warn(msg)
        # raise RuntimeError(msg)
        times = np.arange(len(times))
    column_names = []
    column_values = []

    data_names = self.data_code['data_names']
    # print(f'data_names = {data_names}')
    for name in data_names:
        # if name == primary_name:
        # times = self.da
        times = np.array(getattr(self, name + 's'))
        if name == 'mode':
            column_names.append('Mode')
            column_values.append(times)

            # if freq not in data_names:
            # if name == 'freq':
            # #if hasattr(self, 'freqs'):
            #     column_names.append('Freq')
            #     column_values.append(self.freqs)
            # elif name == 'eigr':
            #     column_names.append('eigenvalue_real')
            #     column_values.append(self.eigrs)
            # elif hasattr(self, 'eigrs') and 0:
            #     try:
            #         abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
            #     except FloatingPointError:
            #         msg = 'Cant analyze freq = sqrt(eig)/(2*pi)\neigr=%s\n' % (self.eigrs)
            #         abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
            #         msg += 'freq = sqrt(abs(self.eigrs)) / (2 * np.pi)=%s' % abs_freqs
            #         raise FloatingPointError(msg)
            #     column_names.append('Freq')
            #     column_values.append(abs_freqs)
            # else:
            #     pass

            # Convert eigenvalues to frequencies
            # TODO: add damping header
        elif name in ['eign']:
            omega_radians, abs_freqs = real_modes_to_omega_freq(self.eigns)
            column_names.append('Freq')
            column_values.append(abs_freqs)
            column_names.append('Eigenvalue')
            column_values.append(times)
            column_names.append('Radians')
            column_values.append(omega_radians)

        elif name in ['eigr']:
            column_names.append('EigenvalueReal')
            column_values.append(times)

        elif name in ['eigi']:
            column_names.append('EigenvalueImag')
            column_values.append(times)
            eigr = np.array(self.eigrs)
            eigi = np.array(self.eigis)
            damping, frequncy = complex_damping_frequency(eigr, eigi)
            column_names.append('Damping')
            column_values.append(damping)
        elif name in ['mode_cycle']:
            continue
            # column_names.append('mode_cycle(Freq?)')
            # column_values.append(times)
        elif name in ['mode2']:
            continue
            # column_names.append('mode2(Freq?)')
            # column_values.append(times)
        elif name in ['cycle']:
            continue
            # column_names.append('Freq (Cycles/s)')
            # column_values.append(times)

        elif name in ['freq']:
            column_names.append('Freq')
            column_values.append(times)
        elif name in ['dt', 'time']:
            column_names.append('Time')
            column_values.append(times)
        elif name in ['lftsfq', 'lsdvmn', 'load_step', 'loadID', 'load_factor', 'loadIDs']:
            column_names.append('LoadStep')
            column_values.append(times)
        elif name == 'node_id':
            column_names.append('NodeID')
            column_values.append(times)
        elif name == 'element_id':
            column_names.append('ElementID')
            column_values.append(times)
        else:  # pragma: no cover
            msg = 'build_dataframe; name=%r' % name
            print(msg)
            raise NotImplementedError(msg)
    assert len(column_names) > 0, column_names
    assert len(column_names) == len(column_values), 'names=%s values=%s' % (column_names, column_values)
    assert len(self.get_headers()) == self.data.shape[-1], 'headers=%s; n=%s\ndata.headers=%s' % (
    self.get_headers(), len(self.get_headers()), self.data.shape[-1])
    return column_names, column_values


def build_pandas_transient_elements(self: ScalarObject,
                                    column_values: np.ndarray, column_names: list[str],
                                    headers: list[str],
                                    element: np.ndarray, data: np.ndarray):
    """common method to build a transient dataframe"""
    import pandas as pd
    columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

    eid_item = []
    for eid in element:
        for header in headers:
            eid_item.append([eid, header])
    ntimes, nelements = data.shape[:2]

    nheaders = len(headers)
    if self.table_name in ['OEFPSD1']:
        ifilter = ~np.isnan(np.max(data[-1, :, :], axis=1))
        if ifilter.sum():
            warnings.warn(f'filtering NaNs from {self.class_name}')
            data2 = data[:, ifilter, :]
            nelements = data2.shape[1]
            A = data2.reshape(ntimes, nelements*nheaders).T
        else:
            A = data.reshape(ntimes, nelements*nheaders).T
    else:
        A = data.reshape(ntimes, nelements*nheaders).T

    names = ['ElementID', 'Item']
    index = pd.MultiIndex.from_tuples(eid_item, names=names)
    try:
        data_frame = pd.DataFrame(A, columns=columns, index=index)
    except ValueError:
        print('A.shape =', A.shape)
        print('len(element) =', len(element))
        print('columns =', columns)
        raise
    # old
    #data_frame = pd.Panel(data, items=column_values, major_axis=element, minor_axis=headers).to_frame()
    #data_frame.columns.names = column_names
    #data_frame.index.names = ['ElementID', 'Item']
    #print(data_frame)
    return data_frame


def build_pandas_transient_from_dict(self: ScalarObject,
                                     column_values, column_names, headers: list[str],
                                     mydict: dict[str, np.ndarray],
                                     data: np.ndarray,
                                     names: list[str]) -> pd.DataFrame:  # pragma: no cover
    """
    common method to build a transient dataframe
    handles mixed type arrays for element_node
    TODO: not done...
    """
    # Freq                  0.00001  10.00000 20.00000 30.00000                 40.00000 50.00000 60.00000
    # ElementID NodeID Item
    # 1         0      oxx        0j       0j       0j       0j    (3200.0806+6017.714j)       0j       0j
    #                  oyy        0j       0j       0j       0j    (410.68146+772.2816j)       0j       0j
    #                  ozz        0j       0j       0j       0j    (0.306115+0.5756457j)       0j       0j
    #                  txy        0j       0j       0j       0j  (-120.69606-226.96753j)       0j       0j
    #                  tyz        0j       0j       0j       0j  (0.70554054+1.3267606j)       0j       0j
    #                  txz        0j       0j       0j       0j     (5193.834+9766.943j)       0j       0j
    # 2                oxx        0j       0j       0j       0j    (8423.371+15840.051j)       0j       0j
    #                  oyy        0j       0j       0j       0j    (-3364.359-6326.637j)       0j       0j
    #                  ozz        0j       0j       0j       0j  (-74931.664-140908.11j)       0j       0j
    #                  txy        0j       0j       0j       0j  (-261.20972-491.20178j)       0j       0j
    #                  tyz        0j       0j       0j       0j   (121.57285+228.61633j)       0j       0j
    #                  txz        0j       0j       0j       0j     (5072.678+9539.112j)       0j       0j
    import pandas as pd
    columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

    #print(data.shape)
    ntimes, nelements = data.shape[:2]
    nheaders = len(headers)
    try:
        A = data.reshape(ntimes, nelements*nheaders).T
    except ValueError:  # pragma: no cover
        ntotal = ntimes * nelements * nheaders
        print(f'data.shape={data.shape}; ntimes={ntimes} nelements={nelements} nheaders={nheaders}; ntotal={ntotal}')
        raise
    assert ntimes == len(column_values[0]), (ntimes, column_values[0])

    if names is None:
        names = ['ElementID', 'NodeID', 'Item']

    element_node = {}
    nvars = len(element_node)
    assert len(names) == nvars + 1, f'names={names} element_node={element_node} (n={len(element_node)})'
    neid = len(element_node[0])

    data_list = []
    for key, datai in mydict.items():
        irange = np.arange(len(datai))

    if 1:
        all_headers = headers * nelements
        # print(all_headers)
    if 0:
        all_headers = headers * neid
        print(all_headers)
    if 0:
        ndata = len(element_node[0])
        neid = len(np.unique(element_node[0]))
    if 0:
        names_list = []
        for header in headers:
            namei = np.full(nelements, header)
            names_list.append(namei)
        all_headers = np.hstack(names_list)
        print(all_headers)
    eid_nid_item.append(all_headers)
    # print('nheaders = ', len(all_headers))

    index = pd.MultiIndex.from_arrays(eid_nid_item, names=names)
    data_frame = pd.DataFrame(A, columns=columns, index=index)

    #element_node = [element_node[:, 0], element_node[:, 1]]
    #data_frame = pd.Panel(data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
    #data_frame.columns.names = column_names
    #data_frame.index.names = ['ElementID', 'NodeID', 'Item']
    #print(data_frame)
    return data_frame


def build_pandas_transient_element_node(self: ScalarObject,
                                        column_values, column_names,
                                        headers: list[str],
                                        element_node: np.ndarray, data: np.ndarray,
                                        names=None,
                                        from_tuples: bool=True,
                                        from_array: bool=False):
    """
    common method to build a transient dataframe
    TODO: doesn't handle mixed type arrays for element_node
    """
    # Freq                  0.00001  10.00000 20.00000 30.00000                 40.00000 50.00000 60.00000
    # ElementID NodeID Item
    # 1         0      oxx        0j       0j       0j       0j    (3200.0806+6017.714j)       0j       0j
    #                  oyy        0j       0j       0j       0j    (410.68146+772.2816j)       0j       0j
    #                  ozz        0j       0j       0j       0j    (0.306115+0.5756457j)       0j       0j
    #                  txy        0j       0j       0j       0j  (-120.69606-226.96753j)       0j       0j
    #                  tyz        0j       0j       0j       0j  (0.70554054+1.3267606j)       0j       0j
    #                  txz        0j       0j       0j       0j     (5193.834+9766.943j)       0j       0j
    # 2                oxx        0j       0j       0j       0j    (8423.371+15840.051j)       0j       0j
    #                  oyy        0j       0j       0j       0j    (-3364.359-6326.637j)       0j       0j
    #                  ozz        0j       0j       0j       0j  (-74931.664-140908.11j)       0j       0j
    #                  txy        0j       0j       0j       0j  (-261.20972-491.20178j)       0j       0j
    #                  tyz        0j       0j       0j       0j   (121.57285+228.61633j)       0j       0j
    #                  txz        0j       0j       0j       0j     (5072.678+9539.112j)       0j       0j
    import pandas as pd
    columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

    #print(data.shape)
    ntimes, nelements = data.shape[:2]
    nheaders = len(headers)
    try:
        A = data.reshape(ntimes, nelements*nheaders).T
    except ValueError:  # pragma: no cover
        ntotal = ntimes * nelements * nheaders
        print(f'data.shape={data.shape}; ntimes={ntimes} nelements={nelements} nheaders={nheaders}; ntotal={ntotal}')
        raise
    assert ntimes == len(column_values[0]), (ntimes, column_values[0])

    if names is None:
        names = ['ElementID', 'NodeID', 'Item']

    assert not(from_tuples and from_array)
    if from_tuples:
        nvars = element_node.shape[1]
        assert len(names) == nvars + 1, f'names={names} nvars+1={nvars+1}; element_node={element_node} {element_node.shape}'
        eid_nid_item = []
        for eid, nid in element_node:
            for header in headers:
                eid_nid_item.append([eid, nid, header])
        # print('neid_item =', len(eid_nid_item))
        index = pd.MultiIndex.from_tuples(eid_nid_item, names=names)
    elif from_array:
        nvars = len(element_node)
        assert len(names) == nvars + 1, f'names={names} nvars+1={nvars+1:d}; element_node={element_node} (n={len(element_node)})'
        eid_nid_item = []
        neid = len(element_node[0])
        for eid in element_node:
            eidi = np.vstack([eid]*nheaders)
            eid_nid_item.append(eidi.ravel())
        if 1:
            all_headers = headers * nelements
            # print(all_headers)
        if 0:
            all_headers = headers * neid
            print(all_headers)
        if 0:
            ndata = len(element_node[0])
            neid = len(np.unique(element_node[0]))
        if 0:
            names_list = []
            for header in headers:
                namei = np.full(nelements, header)
                names_list.append(namei)
            all_headers = np.hstack(names_list)
            print(all_headers)
        eid_nid_item.append(all_headers)
        # print('nheaders = ', len(all_headers))

        index = pd.MultiIndex.from_arrays(eid_nid_item, names=names)
    else:  # pragma: no cover
        raise RuntimeError('from_tuple, from_array')
    data_frame = pd.DataFrame(A, columns=columns, index=index)

    #element_node = [element_node[:, 0], element_node[:, 1]]
    #data_frame = pd.Panel(data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
    #data_frame.columns.names = column_names
    #data_frame.index.names = ['ElementID', 'NodeID', 'Item']
    #print(data_frame)
    return data_frame
