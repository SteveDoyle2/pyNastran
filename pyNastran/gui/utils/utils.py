"""
Simple pure python utilities used by the GUI
"""
from typing import List


def find_next_value_in_sorted_list(lst: List[int], old: int, new: int):
    """
    iold=1 and inew=2, but there is no value of 2, so we skip to 3
    """
    if new in lst:
        return new
    iold = lst.index(old)

    if new > old:
        list_high = lst[iold+1:]
        list_low = lst[:iold+1]

        # increasing icase
        #  skipping an internal value
        if list_high:
            value = list_high[0]
        else:
            # list_low has to exist because we've wrapped around
            # just return the first value
            value = list_low[0]
    else:
        # decreasing icase
        list_high = lst[iold:]
        list_low = lst[:iold]
        if list_low:
            # take the next value
            # [1, 2, *, 4, 5]
            #  4 -> 2
            #  2 -> 1
            value = list_low[-1]
        else:
            # wrapped around
            # [1, 2, *, 4, 5]
            # 1 -> 5
            value = list_high[-1]
    return value
