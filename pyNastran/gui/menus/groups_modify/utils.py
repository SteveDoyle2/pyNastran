from typing import List
import numpy as np

def get_groups_sorted_by_name(data) -> List[str]:
    keys = []
    values = []
    for key, group in data.items():
        if isinstance(key, int):
            keys.append(key)
            values.append(group.name)
    #print(values)
    isort = np.argsort(keys)
    group_names = list(np.array(values, dtype='object')[isort])
    #print(group_names)
    return group_names
