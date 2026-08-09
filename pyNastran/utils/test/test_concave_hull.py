import unittest
import numpy as np
from pyNastran.utils.concave_hull import (
    concave_hull_2d, order_hull_by_angle,)


class TestConcaveHull(unittest.TestCase):
    """Tests for concave_hull_2d and order_hull_by_angle.

    Tolerances
    ----------
    hull area : atol=1e-10 for exact geometric shapes
    hull closure : first == last point, exact
    point count : exact integer match
    """

    def test_convex_hull_circle(self) -> None:
        """A circle's convex hull should contain all points and close properly."""
        npts = 24
        theta = np.linspace(0, 2 * np.pi, npts, endpoint=False)
        pts = np.column_stack([np.cos(theta), np.sin(theta)])

        hull = concave_hull_2d(pts, alpha=0.0)

        assert hull.shape[1] == 2
        np.testing.assert_array_equal(hull[0], hull[-1])
        assert len(hull) == npts + 1

    def test_convex_hull_square(self) -> None:
        """Four corners of a unit square — hull must be exactly 4 boundary points."""
        pts = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ])
        hull = concave_hull_2d(pts, alpha=0.0)
        np.testing.assert_array_equal(hull[0], hull[-1])
        assert len(hull) == 5

    def test_convex_hull_square_with_interior_point(self) -> None:
        """Interior point should not appear on the convex hull boundary."""
        pts = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ])
        hull = concave_hull_2d(pts, alpha=0.0)
        np.testing.assert_array_equal(hull[0], hull[-1])
        assert len(hull) == 5

    def test_concave_hull_removes_large_triangles(self) -> None:
        """With a high alpha, the hull should be tighter (fewer or equal boundary points)."""
        npts = 40
        theta = np.linspace(0, 2 * np.pi, npts, endpoint=False)
        r = 1.0 + 0.3 * np.cos(3 * theta)
        pts = np.column_stack([r * np.cos(theta), r * np.sin(theta)])

        hull_convex = concave_hull_2d(pts, alpha=0.0)
        hull_concave = concave_hull_2d(pts, alpha=1.5)

        assert len(hull_convex) >= 3
        assert len(hull_concave) >= 3
        np.testing.assert_array_equal(hull_convex[0], hull_convex[-1])
        np.testing.assert_array_equal(hull_concave[0], hull_concave[-1])

    def test_hull_too_few_points(self) -> None:
        """Should raise ValueError with fewer than 3 points."""
        pts = np.array([[0.0, 0.0], [1.0, 1.0]])
        with self.assertRaises(ValueError):
            concave_hull_2d(pts, alpha=0.0)

    def test_order_hull_by_angle_ccw(self) -> None:
        """After ordering, angles should be monotonically increasing."""
        pts = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [-1.0, 0.0],
            [0.0, -1.0],
            [1.0, 0.0],
        ])
        ordered = order_hull_by_angle(pts)
        np.testing.assert_array_equal(ordered[0], ordered[-1])

        open_pts = ordered[:-1]
        centroid = open_pts.mean(axis=0)
        angles = np.arctan2(
            open_pts[:, 1] - centroid[1],
            open_pts[:, 0] - centroid[0],
        )
        diffs = np.diff(angles)
        assert np.all(diffs >= 0), f"angles not monotonic: {angles}"

    def test_convex_hull_ellipse_area(self) -> None:
        """Convex hull of a dense ellipse should approximate pi*a*b.

        Tolerances: 1% relative for npts=200.
        """
        npts = 200
        a, b = 3.0, 1.5
        theta = np.linspace(0, 2 * np.pi, npts, endpoint=False)
        pts = np.column_stack([a * np.cos(theta), b * np.sin(theta)])

        hull = concave_hull_2d(pts, alpha=0.0)
        hull = order_hull_by_angle(hull)

        hull_open = hull[:-1]
        area = 0.5 * np.abs(
            np.sum(
                hull_open[:, 0] * np.roll(hull_open[:, 1], -1)
                - np.roll(hull_open[:, 0], -1) * hull_open[:, 1]
            )
        )
        expected_area = np.pi * a * b
        np.testing.assert_allclose(area, expected_area, rtol=0.01)


if __name__ == '__main':
    unittest.main()
