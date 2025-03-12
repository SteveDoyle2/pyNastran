import numpy as np
from pyNastran.op2.op2 import OP2, read_op2

def main2():

    data_mapper = [
        ('strain', ('cquad4_strain', 'ctria3_strain',
                    'ctria3_composite_strain', 'cquad4_composite_strain',
                    'cbar_strain',

                    'chexa_strain',
                    'cpenta_strain',
                    'cpyram_strain',
                    'ctetra_strain',
                    ),
         )
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
        #'cbar_strain': ('element_node', bar_strain_results),
        #'chexa_strain': ('element_node', solid_strain_results),
        #'cpenta_strain': ('element_node', solid_strain_results),
        #'cpyram_strain': ('element_node', solid_strain_results),
        #'ctetra_strain': ('element_node', solid_strain_results),
    }

    #--------------------------------------------------------
    data_storage = {}
    for grouping, results_strings in data_mapper:
        for result_str in results_strings:
            result_types_to_pull = types_to_result_types_to_pull[result_str]
            if len(result_types_to_pull) == 0:
                continue
            data_storage[(grouping, result_str)] = {
                'id': [], 'min': [], 'max': []}

    #--------------------------------------------------------


    out = {}
    # o11, subcase, file
    op2_filenames = ['a.op2', 'b.op2', 'c.op2']

    for ifile, op2_filename in enumerate(op2_filenames):
        model_results = read_op2(op2_filename)
        for grouping, results_strings in data_mapper:
            results_objs = getattr(model_results, grouping)  # Strain()
            for result_str in results_strings:
                data_dict = getattr(results_objs, result_str)  # strain.ctria3_strain

                # [e11, e22]
                id_str, result_types_to_pull = types_to_result_types_to_pull[result_str]
                if len(result_types_to_pull) == 0:
                    continue

                for subcase, res in data_dict.items():
                    # [e11, e22, e12, e1z, e2z, angle, major, minor, max_shear]
                    headers = res.get_headers()

                    # element_layer
                    id_data = getattr(res, id_str)

                    # data; assuming static
                    data = res.data
                    assert data.shape[0] == 1, data.shape

                    for header in result_types_to_pull:
                        iheader = headers.index(header)

                        data_storage[(grouping, result_str)] = {
                            'id': [], 'min': [], 'max': []}
                        values = data[0, :, iheader]  # assuming static
                        data_storage['id'].append(id_data)
                        data_storage['min'].append(values)
                        data_storage['max'].append(values)



def squash(ids: list[np.ndarray],
           mins: list[np.ndarray],
           maxs: list[np.ndarray]):
    for idi, mini in zip(ids, mins):
        assert len(idi) == len(mini), (ids, mini)
    # these are in order
    ids0 = ids[0]
    ids1 = ids[1]

    min0 = mins[0]
    min1 = mins[1]

    max0 = maxs[0]
    max1 = maxs[1]
    allow_nan = True

    if np.array_equal(ids0, ids1):
        ids_out = ids0
        min_array = np.column_stack(mins)
        max_array = np.column_stack(maxs)
        allow_nan = False
    else:
        # not sure they're the same ids
        if ids0.ndim == 1:
            ids_out, min_array, max_array = _squash1d(ids, mins, maxs)
        elif ids0.ndim == 2:
            ids_out, min_array, max_array = _squash2d(ids, mins, maxs)
        else:
            raise RuntimeError((ids0.shape, ids1.shape))

    if allow_nan:
        mins2 = np.nanmin(min_array, axis=1)
        maxs2 = np.nanmax(max_array, axis=1)
    else:
        mins2 = min_array.min(axis=1)
        maxs2 = max_array.max(axis=1)
    # print('ids_out', ids_out)
    # print('mins2', mins2)
    # print('maxs2', maxs2, '\n')
    assert len(ids_out) == len(mins2)
    return ids_out, mins2, maxs2

def _squash1d(ids: list[np.ndarray],
              mins: list[np.ndarray],
              maxs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    min_array = np.full((n_ids, 2), np.nan, min0.dtype)
    max_array = np.full((n_ids, 2), np.nan, max0.dtype)
    i0 = np.searchsorted(ids_out, ids0)
    i1 = np.searchsorted(ids_out, ids1)
    min_array[i0, 0] = min0
    min_array[i1, 1] = min1
    max_array[i0, 0] = max0
    max_array[i1, 1] = max1
    return ids_out, min_array, max_array

def _squash2d(ids: list[np.ndarray],
              mins: list[np.ndarray],
              maxs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    min_array = np.full((n_ids, 2), np.nan, min0.dtype)
    max_array = np.full((n_ids, 2), np.nan, max0.dtype)
    #print(f'uids_str = {uids_str}')
    #print(f'ids0_str = {ids0_str}')
    i0 = np.searchsorted(uids_str, ids0_str)
    i1 = np.searchsorted(uids_str, ids1_str)
    #print('min_array = ', min_array)
    #print('i0 = ', i0)
    #print('min0 = ', min0)
    min_array[i0, 0] = min0
    min_array[i1, 1] = min1
    max_array[i0, 0] = max0
    max_array[i1, 1] = max1
    return ids_out, min_array, max_array


def main():
    # 2d example - different length
    eid_layer1 = np.array([[1, 1]])
    eid_layer2 = np.array([[1, 2],
                           [1, 3]])
    ids = [eid_layer1, eid_layer2]
    mins = [np.array([1.]), np.array([0., 1.])]
    maxs = [np.array([1.]), np.array([3., 4])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 2d example - different ids - overlap
    eid_layer1 = np.array([[1, 1],
                           [1, 2]])
    eid_layer2 = np.array([[1, 1],
                           [1, 3]])
    ids = [eid_layer1, eid_layer2]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 2d example - different ids - unique
    eid_layer1 = np.array([[1, 1],
                           [1, 2]])
    eid_layer2 = np.array([[1, 4],
                           [1, 5]])
    ids = [eid_layer1, eid_layer2]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 2d example - same ids
    eid_layer_same = np.array([[1, 1],
                               [1, 2]])
    ids = [eid_layer_same, eid_layer_same]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 1d-same length
    ids_same = np.array([1, 2])
    ids = [ids_same, ids_same]
    mins = [np.array([1., 2.]), np.array([0., 1.])]
    maxs = [np.array([1., 2.]), np.array([3., 4])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 1d-different_length
    ids = [np.array([1, 2]), np.array([3, 4, 5])]
    mins = [np.array([1., 2.]), np.array([4., 5., 6])]
    maxs = [np.array([1., 2.]), np.array([4., 5., 6])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 1d-dissimilar ids
    ids = [np.array([1, 2, 3]), np.array([3, 4, 5])]
    mins = [np.array([1., 2., 3.]), np.array([4., 5., 6])]
    maxs = [np.array([1., 2., 3.]), np.array([4., 5., 6])]
    ids2, mins2, maxs2 = squash(ids, mins, maxs)

    # 2d ids example (element_layer, element_node)
    #


if __name__ == '__main__':
    main()
