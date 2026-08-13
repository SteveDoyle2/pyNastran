import numpy as np
from scipy.spatial import Delaunay


def concave_hull_2d(points_2d: np.ndarray,
                    alpha: float = 0.0) -> np.ndarray:
    """Compute the concave hull (alpha shape) of a 2D point set.

    Uses the Delaunay triangulation and removes triangles whose
    circumradius exceeds 1/alpha.  When alpha=0 the convex hull
    is returned.

    Parameters
    ----------
    points_2d : (n, 2) float ndarray
        2D point cloud (e.g. y-z of a cross-section)
    alpha : float
        Controls concavity.  Larger values produce tighter hulls.
        0 gives the convex hull.

    Returns
    -------
    hull_points : (m, 2) float ndarray
        Ordered boundary points of the concave hull, closed
        (first == last).
    """
    if len(points_2d) < 3:
        raise ValueError(f"need >= 3 points for a hull, got {len(points_2d)}")

    tri = Delaunay(points_2d)
    triangles = tri.simplices

    if alpha > 0.0:
        keep = np.ones(len(triangles), dtype=bool)
        for i, simplex in enumerate(triangles):
            pa = points_2d[simplex[0]]
            pb = points_2d[simplex[1]]
            pc = points_2d[simplex[2]]
            a = np.linalg.norm(pb - pa)
            b = np.linalg.norm(pc - pb)
            c = np.linalg.norm(pa - pc)
            s = (a + b + c) / 2.0
            area = np.sqrt(max(s * (s - a) * (s - b) * (s - c), 0.0))
            if area == 0.0:
                keep[i] = False
                continue
            circum_r = (a * b * c) / (4.0 * area)
            if circum_r > 1.0 / alpha:
                keep[i] = False
        triangles = triangles[keep]

    edges: dict[tuple[int, int], int] = {}
    for simplex in triangles:
        for j in range(3):
            edge = (simplex[j], simplex[(j + 1) % 3])
            edge_key = (min(edge), max(edge))
            edges[edge_key] = edges.get(edge_key, 0) + 1

    boundary_edges = [e for e, cnt in edges.items() if cnt == 1]

    if not boundary_edges:
        raise ValueError("no boundary edges found — check alpha value")

    adj: dict[int, list[int]] = {}
    for a_idx, b_idx in boundary_edges:
        adj.setdefault(a_idx, []).append(b_idx)
        adj.setdefault(b_idx, []).append(a_idx)

    start = boundary_edges[0][0]
    ordered = [start]
    prev = -1
    current = start
    for _ in range(len(boundary_edges) + 1):
        neighbors = adj[current]
        next_node = neighbors[0] if neighbors[0] != prev else neighbors[1]
        if next_node == start:
            break
        ordered.append(next_node)
        prev = current
        current = next_node

    ordered.append(ordered[0])
    return points_2d[ordered]


def order_hull_by_angle(hull_yz: np.ndarray) -> np.ndarray:
    """Re-order hull points by angle around the centroid.

    Parameters
    ----------
    hull_yz : (m, 2) float ndarray
        closed hull points (first == last)

    Returns
    -------
    ordered : (m, 2) float ndarray
        points ordered CCW by angle, closed
    """
    pts = hull_yz[:-1]
    centroid = pts.mean(axis=0)
    angles = np.arctan2(pts[:, 1] - centroid[1], pts[:, 0] - centroid[0])
    order = np.argsort(angles)
    ordered = pts[order]
    return np.vstack([ordered, ordered[0:1]])

