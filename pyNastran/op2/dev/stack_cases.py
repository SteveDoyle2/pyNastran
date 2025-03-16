import numpy as np
from cpylog import SimpleLogger
from pyNastran.op2.op2 import OP2, read_op2

def find_critical_cases(op2_filenames: list[str]):
    log = SimpleLogger(level='debug')
    strain_results = (
        'cquad4_strain', 'ctria3_strain',
        'ctria3_composite_strain', 'cquad4_composite_strain',
        'cbar_strain',
        'chexa_strain', 'cpenta_strain', 'cpyram_strain', 'ctetra_strain',
    )
    data_mapper = [
        # desired results
        ('strain', strain_results)
    ]
    beam_strain_results = ['emaxa', 'emina']
    bar_strain_results = ['emaxa', 'emina']
    shell_strain_results = []
    composite_strain_results = ['e11', 'e22', 'e12']
    solid_strain_results = []
    types_to_result_types_to_pull = {
        'cquad4_strain': ('element_node', shell_strain_results),
        'ctria3_strain': ('element_node', shell_strain_results),
        'ctria3_composite_strain': ('element_layer', composite_strain_results),
        'cquad4_composite_strain': ('element_layer', composite_strain_results),
        #'cbar_strain': ('element', bar_strain_results),
        #'chexa_strain': ('element_node', solid_strain_results),
        #'cpenta_strain': ('element_node', solid_strain_results),
        #'cpyram_strain': ('element_node', solid_strain_results),
        #'ctetra_strain': ('element_node', solid_strain_results),
    }

    #--------------------------------------------------------
    data_storage = _setup_data_storage(
        log, data_mapper, types_to_result_types_to_pull)
    #--------------------------------------------------------
    #out = {}
    # o11, subcase, file

    # pull data
    icase = 0
    case_map = {}
    for ifile, op2_filename in enumerate(op2_filenames):
        model_results = read_op2(op2_filename, debug=False)
        op2_results = model_results.op2_results
        for grouping, results_strings in data_mapper:
            results_objs = getattr(op2_results, grouping)  # Strain()
            for result_str in results_strings:
                if result_str not in types_to_result_types_to_pull:
                    #log.warning(f'no results for {result_str} found')
                    continue

                # [e11, e22]
                id_str, result_types_to_pull = types_to_result_types_to_pull[result_str]
                if len(result_types_to_pull) == 0:
                    continue

                data_dict = getattr(results_objs, result_str)  # strain.ctria3_strain
                for subcase, res in data_dict.items():
                    # [e11, e22, e12, e1z, e2z, angle, major, minor, max_shear]
                    headers = res.get_headers()

                    # element_layer
                    try:
                        id_data = getattr(res, id_str)
                    except AttributeError:
                        log.error('\n' + ''.join(res.get_stats()))
                        raise

                    # data; assuming static
                    data = res.data
                    assert data.shape[0] == 1, data.shape

                    for header in result_types_to_pull:
                        key = (grouping, result_str, header)
                        iheader = headers.index(header)

                        # (strain, ctria3_strain, e11)
                        #data_storage[(grouping, result_str, header)] = {
                            #'id': [], 'case_min': [], min': [], 'case_max': [], 'max': []}
                        values = data[0, :, iheader]  # assuming static
                        nvalues = len(values)
                        cases = np.ones(nvalues, dtype='int32') * icase
                        #print(data_storage)
                        #data_storage[key]['case'].append(cases)
                        data_storage[key]['id'].append(id_data)
                        data_storage[key]['case_min'].append(cases)
                        data_storage[key]['case_max'].append(cases)
                        data_storage[key]['min'].append(values)
                        data_storage[key]['max'].append(values)

                        ids_out, case_min, mins_out, case_max, maxs_out = squash(
                            data_storage[key]['id'],
                            data_storage[key]['case_min'],
                            data_storage[key]['min'],
                            data_storage[key]['case_max'],
                            data_storage[key]['max'])
                        data_storage[key] = {
                            'id': ids_out,
                            'case_min': case_min,
                            'min': mins_out,
                            'case_max': case_max,
                            'max': maxs_out,
                        }
                    case_map[icase] = (ifile, op2_filename, subcase)
                    icase += 1
    for icase, (ifile, pathi, subcase) in case_map.items():
        # ifile, op2_filename, subcase
        print(f'{icase}: ifile={ifile} subcase={subcase} path={str(pathi)}')
    print(key)
    print(data_storage[key])
    return data_storage, case_map

def _setup_data_storage(log: SimpleLogger,
                        data_mapper,
                        types_to_result_types_to_pull) -> dict[tuple[str, str, str],
                                                               dict[str, np.ndarray]]:
    # size the output data
    data_storage = {}
    for grouping, results_strings in data_mapper:
        for result_str in results_strings:
            if result_str not in types_to_result_types_to_pull:
                log.warning(f'no results for {result_str} found')
                continue
            id_str, result_types_to_pull = types_to_result_types_to_pull[result_str]
            if len(result_types_to_pull) == 0:
                continue
            for header in result_types_to_pull:
                for result_str in results_strings:
                    # (strain, ctria3_strain, e11)
                    key = (grouping, result_str, header)
                    #print(key)
                    data_storage[key] = {
                        'id': [],
                        'min': [], 'case_min': [],
                        'max': [], 'case_max': []
                    }
    return data_storage

def squash(ids: list[np.ndarray],
           cases_min: list[np.ndarray],
           mins: list[np.ndarray],
           cases_max: list[np.ndarray],
           maxs: list[np.ndarray]):
    # validation
    assert len(cases_min), f'{len(cases_min)} cases found'
    assert len(cases_max), f'{len(cases_max)} cases found'
    assert len(cases_min) == len(ids), f'cases={cases_min} ids={ids}'
    assert len(cases_min) == len(mins)
    assert len(cases_min) == len(maxs)
    assert len(cases_min) == len(cases_max), f'cases={cases_min} cases_max={cases_max}'

    if len(ids) == 1:
        # trivial case
        return ids, cases_min, mins, cases_max, maxs

    # simplified for testing
    if isinstance(cases_min, list) and isinstance(cases_min[0], int):
        cases_min = [np.ones(len(idi), dtype='int32')*case for case, idi in zip(cases_min, ids)]
    if isinstance(cases_max, list) and isinstance(cases_max[0], int):
        cases_max = [np.ones(len(idi), dtype='int32')*case for case, idi in zip(cases_max, ids)]

    for idi, case_min, mini, case_max, maxi in zip(ids, cases_min, mins, cases_max, maxs):
        assert len(idi) == len(case_min), (ids, case_min)
        assert len(idi) == len(case_max), (ids, case_max)
        assert len(idi) == len(mini), (ids, mini)
        assert len(idi) == len(maxi), (ids, maxi)
    #------------------------------------------------------------------

    # these are in order
    ids0 = ids[0]
    ids1 = ids[1]
    allow_nan = True

    if np.array_equal(ids0, ids1):
        ids_out = ids0
        #for mini in mins:
            #print(mini.shape)

        min_array = np.column_stack(mins)
        max_array = np.column_stack(maxs)
        cases_array_min = np.column_stack(cases_min)
        cases_array_max = np.column_stack(cases_max)
        allow_nan = False
    else:
        # not sure they're the same ids
        if ids0.ndim == 1:
            ids_out, cases_array_min, min_array, cases_array_max, max_array = _squash1d(
                ids, cases_min, mins, cases_max, maxs)
        elif ids0.ndim == 2:
            ids_out, cases_array_min, min_array, cases_array_max, max_array = _squash2d(
                ids, cases_min, mins, cases_max, maxs)
        else:
            raise RuntimeError((ids0.shape, ids1.shape))

    #if allow_nan:
    mins2 = np.nanmin(min_array, axis=1)
    maxs2 = np.nanmax(max_array, axis=1)
    imin_out = np.nanargmin(min_array, axis=1, keepdims=True)
    imax_out = np.nanargmax(max_array, axis=1, keepdims=True)
    # else:
    #     mins2 = min_array.min(axis=1)
    #     maxs2 = max_array.max(axis=1)
    #     # find the index of the critical case
    #     imin_out = np.argmin(min_array, axis=1)
    #     imax_out = np.argmax(max_array, axis=1)

    # map the index to the case array
    if 0:  # pragma: no cover
        case_min = cases_array_min[imin_out]
        case_max = cases_array_max[imax_out]
    else:
        # requires keepdims=True
        case_min = np.take_along_axis(cases_array_min, imin_out, axis=1).flatten()
        case_max = np.take_along_axis(cases_array_max, imax_out, axis=1).flatten()

    #print('min_array', min_array)
    #print('cases_array_min', cases_array_min)
    #print('imin_out', imin_out)
    #print('case_min', case_min)
    assert case_min.shape == mins2.shape
    # print('ids_out', ids_out)
    # print('mins2', mins2)
    # print('maxs2', maxs2, '\n')
    assert len(ids_out) == len(mins2)
    return [ids_out], [case_min], [mins2], [case_max], [maxs2]

def _squash1d(case_array: list[np.ndarray],
              ids: list[np.ndarray],
              cases_min: list[np.ndarray],
              mins: list[np.ndarray],
              cases_max: list[np.ndarray],
              maxs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # these are in order
    ids0 = ids[0]
    ids1 = ids[1]

    min0 = mins[0]
    min1 = mins[1]

    max0 = maxs[0]
    max1 = maxs[1]

    # 1d arrays
    ids_out = np.unique(np.hstack(ids))
    n_ids = len(ids_out)
    cases_array_min = np.full((n_ids, 2), -1, dtype='int32')
    cases_array_max = np.full((n_ids, 2), -1, dtype='int32')
    min_array = np.full((n_ids, 2), np.nan, min0.dtype)
    max_array = np.full((n_ids, 2), np.nan, max0.dtype)
    i0 = np.searchsorted(ids_out, ids0)
    i1 = np.searchsorted(ids_out, ids1)

    cases_array_min[i0, 0] = cases_min[0]
    cases_array_min[i1, 1] = cases_min[1]
    cases_array_max[i0, 0] = cases_max[0]
    cases_array_max[i1, 1] = cases_max[1]

    min_array[i0, 0] = min0
    min_array[i1, 1] = min1
    max_array[i0, 0] = max0
    max_array[i1, 1] = max1
    return ids_out, cases_array_min, min_array, cases_array_max, max_array

def _squash2d(case_array: list[np.ndarray],
              ids: list[np.ndarray],
              cases_min: list[np.ndarray],
              mins: list[np.ndarray],
              cases_max: list[np.ndarray],
              maxs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # these are in order
    ids0 = ids[0]
    ids1 = ids[1]

    min0 = mins[0]
    min1 = mins[1]

    max0 = maxs[0]
    max1 = maxs[1]

    ids0_str = [f'{idi[0]}_{idi[1]}' for idi in ids0]
    ids1_str = [f'{idi[0]}_{idi[1]}' for idi in ids1]
    uids_str = np.unique(np.hstack([ids0_str, ids1_str]))
    ids_out = np.array([val.split('_') for val in uids_str], dtype=ids0.dtype)

    n_ids = len(ids_out)
    cases_array_min = np.full((n_ids, 2), -1, dtype='int32')
    cases_array_max = np.full((n_ids, 2), -1, dtype='int32')
    min_array = np.full((n_ids, 2), np.nan, min0.dtype)
    max_array = np.full((n_ids, 2), np.nan, max0.dtype)
    #print(f'uids_str = {uids_str}')
    #print(f'ids0_str = {ids0_str}')
    i0 = np.searchsorted(uids_str, ids0_str)
    i1 = np.searchsorted(uids_str, ids1_str)
    #print('min_array = ', min_array)
    #print('i0 = ', i0)
    #print('min0 = ', min0)

    cases_array_min[i0, 0] = cases_min[0]
    cases_array_min[i1, 1] = cases_min[1]
    cases_array_max[i0, 0] = cases_max[0]
    cases_array_max[i1, 1] = cases_max[1]
    min_array[i0, 0] = min0
    min_array[i1, 1] = min1
    max_array[i0, 0] = max0
    max_array[i1, 1] = max1
    return ids_out, cases_array_min, min_array, cases_array_max, max_array


def squash_test():
    cases = [1001, 1010]

    # 2d example - different length
    eid_layer1 = np.array([[1, 1]])
    eid_layer2 = np.array([[1, 2],
                           [1, 3]])
    ids = [eid_layer1, eid_layer2]
    mins = [np.array([1.]), np.array([0., 1.])]
    maxs = [np.array([1.]), np.array([3., 4])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

    # 2d example - different ids - overlap
    eid_layer1 = np.array([[1, 1],
                           [1, 2]])
    eid_layer2 = np.array([[1, 1],
                           [1, 3]])
    ids = [eid_layer1, eid_layer2]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

    # 2d example - different ids - unique
    eid_layer1 = np.array([[1, 1],
                           [1, 2]])
    eid_layer2 = np.array([[1, 4],
                           [1, 5]])
    ids = [eid_layer1, eid_layer2]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

    # 2d example - same ids
    eid_layer_same = np.array([[1, 1],
                               [1, 2]])
    ids = [eid_layer_same, eid_layer_same]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

    # 1d-same length
    ids_same = np.array([1, 2])
    ids = [ids_same, ids_same]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

    # 1d-different_length
    ids = [np.array([1, 2]), np.array([3, 4, 5])]
    mins = [np.array([1., 2.]), np.array([4., 5., 6])]
    maxs = [np.array([1., 2.]), np.array([4., 5., 6])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

    # 1d-dissimilar ids
    ids = [np.array([1, 2, 3]), np.array([3, 4, 5])]
    mins = [np.array([1., 2., 3.]), np.array([4., 5., 6])]
    maxs = [np.array([1., 2., 3.]), np.array([4., 5., 6])]
    ids2, imin2, mins2, imax2, maxs2 = squash(cases, ids, mins, maxs)

def main():
    from pathlib import Path
    import pyNastran
    pkg_path = Path(pyNastran.__path__[0])
    model_dir = pkg_path / '..' / 'models' / 'sol_101_elements'
    op2_filenames = [model_dir / 'static_solid_shell_bar.op2',
                     model_dir / 'static_solid_shell_bar_xyz.op2']
    find_critical_cases(op2_filenames)

    #squash_test()

if __name__ == '__main__':
    main()
