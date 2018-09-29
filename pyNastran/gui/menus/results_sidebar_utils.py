"""helper for results_sidebar.py"""
from __future__ import print_function
from six import string_types, integer_types


def build_pruned_tree(tree, cases):
    """
    Build new tree with only the specifiied cases

    Parameters
    ----------
    tree : List[...]
        the sidebar tree
    cases : List[int]
        the cases to keep

    Returns
    -------
    tree_final : List[...]
        the updated sidebar tree

    Examples
    --------
    form = [
        [u'Geometry', None, [
            (u'NodeID', 0, []),
            (u'ElementID', 1, []),
            (u'PropertyID', 2, []),
            (u'MaterialID', 3, []),
            (u'E', 4, []),
            (u'Element Checks', None, [
                (u'ElementDim', 5, []),
                (u'Min Edge Length', 6, []),
                (u'Min Interior Angle', 7, []),
                (u'Max Interior Angle', 8, [])],
            ),],
        ],
    ]
    cases = [1, 2, 5]
    tree2 = prune_tree(form, cases)
    >>> tree2
    [
        [u'Geometry', None, [
            (u'ElementID', 1, []),
            (u'PropertyID', 2, []),
            (u'Element Checks', None, [
                (u'ElementDim', 5, []),
            ),],
        ],
    ]
    """
    tree_final = []
    is_results, tree_final = _build_pruned_tree(tree, cases, tree_final)
    return tree_final

def _build_pruned_tree(tree, cases, tree2):
    """helper method for ``build_pruned_tree``"""
    is_results = False
    if isinstance(tree[0], string_types):
        try:
            name, icase, cases2 = tree
        except ValueError:
            print(tree)
            raise

        if isinstance(icase, integer_types):
            assert cases2 == [], tree
            if icase in cases:
                is_results = True
                tree2.append(tree)
        else:
            assert icase is None, tree
            tree3 = []
            for case in cases2:
                is_resultsi, tree3 = _build_pruned_tree(case, cases, tree3)
                if is_resultsi:
                    is_results = True
            if is_results:
                tree2.append((name, icase, tree3))
    else:
        tree3 = []
        for case in tree:
            is_resultsi, tree3 = _build_pruned_tree(case, cases, tree3)
            if is_resultsi:
                is_results = True
        tree2 = tree3
    return is_results, tree2


def get_cases_from_tree(tree):
    """
    Get the cases found in the tree

    Parameters
    ----------
    tree : List[...]
        the sidebar tree

    Returns
    -------
    cases : List[int]
        the cases in the tree

    Examples
    --------
    form = [
        [u'Geometry', None, [
            (u'NodeID', 0, []),
            (u'ElementID', 1, []),
            (u'PropertyID', 2, []),
            (u'MaterialID', 3, []),
            (u'E', 4, []),
            (u'Element Checks', None, [
                (u'ElementDim', 5, []),
                (u'Min Edge Length', 6, []),
                (u'Min Interior Angle', 7, []),
                (u'Max Interior Angle', 8, [])],
            ),],
        ],
    ]
    cases = get_cases_from_tree(form)
    >>> cases
    [0, 1, 2, 3, 4, 5, 6, 7, 8]
    """
    cases = []
    if isinstance(tree[0], string_types):
        try:
            name, icase, cases2 = tree
        except ValueError:  # pragma: no cover
            print(tree)
            raise
        if isinstance(icase, integer_types):
            assert cases2 == [], tree
            cases.append(icase)
        else:
            assert icase is None or icase == '', tree
            icase = None
            for case in cases2:
                cases += get_cases_from_tree(case)
    else:
        for case in tree:
            cases += get_cases_from_tree(case)
    return cases
