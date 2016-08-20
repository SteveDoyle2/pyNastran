
from __future__ import print_function
from six import iteritems
from pyNastran.bdf.bdf import read_bdf, CTRIA3

def split_elements(bdf_filename):
    model = read_bdf(bdf_filename, xref=True)
    for eid, elem in iteritems(model.elements):
        if elem.type == 'CTRIA3':
            #
            #        3
            #       /|\
            #      / | \
            #     /  |  \
            #    /   4   \
            #   /  /   \  \
            #  / /       \ \
            # 1-------------2
            #
            p1, p2, p3 = elem.get_node_positions()
            centroid = (p1 + p2 + p3) / 3.

            #
            #      3
            #     /|\
            #    / | \
            #   /  |  \
            #  /   |   \
            # 1----4----2
            #
        elif elem.type == 'CQUAD4':
            #
            #
            # 4---------3
            # | \     / |
            # |   \  /  |
            # |    5    |
            # |  /   \  |
            # |/       \|
            # 1---------2
            #
            # the same thing shown in a rotated view
            #           4
            #          /| \
            #       /   |   \
            #     /     |     \
            #   /       |       \
            # 1---------5---------3
            #   \       |       /
            #     \     |     /
            #       \   |   /
            #         \ | /
            #           2
            #
            # max_area, taper_ratio, area_ratio
            # 4----7----3
            # |    |    |
            # |    |    |
            # 8----9----6
            # |    |    |
            # |    |    |
            # 1----4----2
            #
            # max_interior_angle
            #      4---------3
            #     / \       /
            #    /   \     /
            #   /     \   /
            #  /       \ /
            # 1---------2
            #
            # taper_ratio
            #     4--6--3
            #    /   |   \
            #   /    |    \
            #  /     |     \
            # 1------5------2
            #
            # taper_ratio
            #     4------3
            #    / \    / \
            #   /   \  /   \
            #  /     \/     \
            # 1-------5------2
            #
            # taper_ratio
            #     4------3
            #    / \      \
            #   /   \      \
            #  /     \      \
            # 1-------5------2
            pass
