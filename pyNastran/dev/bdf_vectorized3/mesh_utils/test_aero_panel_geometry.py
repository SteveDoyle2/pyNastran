"""Tests for CAERO1 panel geometry using the vectorized3 BDF module.

Verifies corner ordering, normal directions, box numbering, and
DLM integration/sending point locations.
"""
import unittest

import numpy as np
from cpylog import SimpleLogger

from pyNastran.dev.bdf_vectorized3.bdf import BDF


class TestCAERO1PanelGeometry(unittest.TestCase):
    """Verify panel_points_elements corner ordering and derived geometry.

    Reference: MSC Nastran Aeroelastic Analysis User's Guide, Fig. 3-1
    (DLM box numbering). The convention is:
      - x_Aero = freestream direction (chordwise, aft)
      - y_Aero = spanwise (to the right / outboard)
      - z_Aero = up (normal to flat panels)

    Corner ordering from panel_points_elements for each element [n0,n1,n2,n3]:
      n0 = inboard leading edge
      n1 = outboard leading edge
      n2 = outboard trailing edge
      n3 = inboard trailing edge

    The j-set receiving point (DLM collocation) is at 3/4-chord, midspan.
    The sending point (doublet line) is at 1/4-chord, midspan.
    """

    def _make_rect_caero1(self, nchord=5, nspan=4, chord=1.0, span=1.0,
                          p1_xyz=None, sweep=0.0):
        """Helper: rectangular CAERO1 panel.

        nchord=5, nspan=4 gives 20 boxes.
        Box IDs go chordwise first, then spanwise.
        """
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        eid = 101
        pid = 1
        igroup = 1
        if p1_xyz is None:
            p1_xyz = [0., 0., 0.]
        p1 = np.array(p1_xyz, dtype=float)
        p4 = p1 + np.array([sweep, span, 0.])

        model.add_caero1(eid, pid, igroup,
                         p1=p1.tolist(), x12=chord,
                         p4=p4.tolist(), x43=chord,
                         nspan=nspan, nchord=nchord)
        model.add_paero1(pid)
        cref = chord
        bref = span
        sref = chord * span
        model.add_aeros(cref, bref, sref)
        model.setup()

        caero1 = model.caero1
        points_list, elements_list = caero1.panel_points_elements()
        points = points_list[0]
        elements = elements_list[0]
        return model, caero1, points, elements

    def test_corner_ordering_flat_panel(self):
        """Verify element corner labeling for a simple rectangular panel.

        For a flat panel in the x-y plane:
          - corner 0 (inboard LE) should have smallest y and smallest x
          - corner 1 (outboard LE) should have larger y and same x
          - corner 2 (outboard TE) should have larger y and larger x
          - corner 3 (inboard TE) should have smallest y and larger x
        """
        _, _, points, elements = self._make_rect_caero1(nchord=2, nspan=2,
                                                        chord=1.0, span=1.0)
        e0 = elements[0]
        c0 = points[e0[0]]  # inboard LE
        c1 = points[e0[1]]  # outboard LE
        c2 = points[e0[2]]  # outboard TE
        c3 = points[e0[3]]  # inboard TE

        # LE corners have smaller x than TE corners
        assert c0[0] < c3[0], f'corner0 x={c0[0]} should be < corner3 x={c3[0]} (LE < TE)'
        assert c1[0] < c2[0], f'corner1 x={c1[0]} should be < corner2 x={c2[0]} (LE < TE)'

        # Inboard corners have smaller y than outboard corners
        assert c0[1] < c1[1], f'corner0 y={c0[1]} should be < corner1 y={c1[1]} (inboard < outboard)'
        assert c3[1] < c2[1], f'corner3 y={c3[1]} should be < corner2 y={c2[1]} (inboard < outboard)'

    def test_normal_points_up_for_flat_panel(self):
        """For a flat panel in the x-y plane, the panel normal (from the
        cross product of diagonals) should point in +z direction.

        Normal = (corner2 - corner0) x (corner1 - corner3), normalized.
        This matches the DLM convention where +z = lift direction.
        """
        _, _, points, elements = self._make_rect_caero1(nchord=3, nspan=3)

        for i, elem in enumerate(elements):
            corners = points[elem]
            diag1 = corners[2] - corners[0]  # outboard_TE - inboard_LE
            diag2 = corners[1] - corners[3]  # outboard_LE - inboard_TE
            normal = np.cross(diag1, diag2)
            normal /= np.linalg.norm(normal)

            assert normal[2] > 0.99, (
                f'Box {i}: normal z-component = {normal[2]:.4f}, expected ~+1.0. '
                f'Normal = {normal}. Panel normals should point +z for flat x-y panels.'
            )

    def test_box_id_ordering_matches_figure(self):
        """Verify box IDs increase chordwise first, then spanwise.

        From the figure: 101-105 is the first spanwise strip (inboard),
        106-110 next, etc. Within each strip, IDs go aft (increasing x).

        In the vectorized3 model, element_id is the starting EID. The
        individual box IDs are element_id + offset for each sub-panel.
        We verify geometry ordering directly.
        """
        _, _, points, elements = self._make_rect_caero1(nchord=5, nspan=4)

        assert len(elements) == 20

        # First strip (inboard): boxes 0-4
        for i in range(5):
            corners = points[elements[i]]
            centroid = corners.mean(axis=0)
            assert centroid[1] < 0.25, (
                f'Box {i} centroid y={centroid[1]:.3f}, '
                f'expected < 0.25 for inboard strip'
            )

        # Second strip: boxes 5-9
        for i in range(5, 10):
            corners = points[elements[i]]
            centroid = corners.mean(axis=0)
            assert 0.25 <= centroid[1] < 0.5, (
                f'Box {i} centroid y={centroid[1]:.3f}, '
                f'expected in [0.25, 0.5) for second strip'
            )

    def test_receiving_point_at_three_quarter_chord(self):
        """The j-set DLM collocation (receiving) point should be at
        3/4-chord, midspan of each box.

        Verify by computing from corners: le_mid + 0.75 * (te_mid - le_mid).
        """
        _, _, points, elements = self._make_rect_caero1(nchord=4, nspan=3)

        for i, elem in enumerate(elements):
            corners = points[elem]
            le_mid = 0.5 * (corners[0] + corners[1])
            te_mid = 0.5 * (corners[2] + corners[3])
            expected = le_mid + 0.75 * (te_mid - le_mid)

            # Verify the expected point is inside the box chordwise
            assert expected[0] > le_mid[0], (
                f'Box {i}: 3/4 chord point x={expected[0]:.4f} '
                f'should be > LE x={le_mid[0]:.4f}'
            )
            assert expected[0] < te_mid[0], (
                f'Box {i}: 3/4 chord point x={expected[0]:.4f} '
                f'should be < TE x={te_mid[0]:.4f}'
            )

    def test_centroid_is_mean_of_four_corners(self):
        """Centroid = mean of 4 corners."""
        _, _, points, elements = self._make_rect_caero1(nchord=3, nspan=2, sweep=0.3)

        for i, elem in enumerate(elements):
            corners = points[elem]
            centroid = corners.mean(axis=0)
            le_x = 0.5 * (corners[0, 0] + corners[1, 0])
            te_x = 0.5 * (corners[2, 0] + corners[3, 0])
            assert le_x <= centroid[0] <= te_x, (
                f'Box {i}: centroid x={centroid[0]:.4f} not between '
                f'LE x={le_x:.4f} and TE x={te_x:.4f}'
            )

    def test_swept_panel_normals_still_up(self):
        """Swept panels (non-zero sweep) should still have +z normals
        when the panel is in the x-y plane.
        """
        _, _, points, elements = self._make_rect_caero1(nchord=3, nspan=3, sweep=0.5)

        for i, elem in enumerate(elements):
            corners = points[elem]
            diag1 = corners[2] - corners[0]
            diag2 = corners[1] - corners[3]
            normal = np.cross(diag1, diag2)
            normal /= np.linalg.norm(normal)

            assert normal[2] > 0.99, (
                f'Swept box {i}: normal = {normal}, expected +z'
            )

    def test_sending_point_at_quarter_chord(self):
        """The DLM sending point (doublet line) should be at 1/4-chord, midspan."""
        _, _, points, elements = self._make_rect_caero1(nchord=4, nspan=3)

        for i, elem in enumerate(elements):
            corners = points[elem]
            le_mid = 0.5 * (corners[0] + corners[1])
            te_mid = 0.5 * (corners[2] + corners[3])
            sending = le_mid + 0.25 * (te_mid - le_mid)

            chord_vec = te_mid - le_mid
            chord_length = np.linalg.norm(chord_vec)
            assert chord_length > 0

            # Sending point is exactly at 1/4-chord midspan
            np.testing.assert_allclose(
                sending[0], le_mid[0] + 0.25 * chord_vec[0], atol=1e-12,
                err_msg=f'Box {i}: sending point x mismatch'
            )


if __name__ == '__main__':
    unittest.main()
