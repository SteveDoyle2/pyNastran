from typing import Any
import numpy as np
Group = Any

def get_groups_sorted_by_name(data: dict[str, Group]) -> list[str]:
    return [group.name for key, group in data.keys() if isinstance(key, int)]

    keys: list[str] = []
    values: list[str] = []
    for key, group in data.items():
        if isinstance(key, int):
            keys.append(key)
            values.append(group.name)
    #print(values)
    isort = np.argsort(keys)
    group_names = list(np.array(values, dtype='object')[isort])
    #print(group_names)
    return group_names
