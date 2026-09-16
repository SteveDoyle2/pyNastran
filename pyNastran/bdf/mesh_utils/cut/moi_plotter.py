"""
Section-stiffness extraction along a structural span.

Cuts a shell (and optionally bar/beam) FEM at user-specified stations,
integrates the section properties at each cut, and produces:

- area, second moments of area, and polar moment (A, I, J)
- modulus-weighted stiffnesses (EA, EI, GJ) with composite laminate support
- Bredt-Batho torsion and transverse-shear-flow shear center
- an equivalent CBEAM stick model that reproduces every stiffness exactly
- span-wise CSV and PNG plots of all quantities

The primary entry point is ``cut_and_plot_moi``.
"""
from __future__ import annotations
import os
import copy
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING
from itertools import count

import numpy as np
try:
    import matplotlib.pyplot as plt  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.coordinate_systems import (
    CORD2R, Coord,
    xyz_to_rtz_array, rtz_to_xyz_array)
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.cut.torsion import bredt_batho_gj
from pyNastran.bdf.mesh_utils.cut.shear_center import shear_center
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    cut_face_model_by_coord, _setup_faces,
    # is_element_cut,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.cards.elements.shell import CTRIA3, CQUAD4
Rods = tuple[np.ndarray, np.ndarray, np.ndarray]


def cut_and_plot_moi(bdf_filename: PathLike | BDF,
                     normal_plane: np.ndarray,
                     log: SimpleLogger,
                     stations: list[float] | np.ndarray,
                     coords: list[CORD2R],
                     x_vector: list[float],
                     include_lines: bool=False,
                     include_solids: bool=False,
                     face_data: Optional[Any]=None,
                     dirname: PathLike='',
                     ifig: int=1,
                     debug_vectorize: bool=True,
                     debug_v3: bool=False,
                     rho: float=1.0,
                     xyz_round: int | None=None,
                     area_round: int | None=None,
                     inertia_round: int | None=None,
                     beam_grid_xyz: Optional[np.ndarray]=None,
                     beam_grid_ids: Optional[np.ndarray]=None,
                     beam_id0: int=1,
                     stop_on_failure: bool=False,
                     cut_data_span_filename: PathLike='cut_data_vs_span.csv',
                     beam_model_bdf_filename: PathLike='equivalent_beam_model.bdf',
                     thetas_csv_filename: PathLike='thetas.csv',
                     normalized_inertia_png_filename: PathLike='normalized_inertia_vs_span.png',
                     area_span_png_filename: PathLike='area_vs_span.png',
                     amoi_span_png_filename: PathLike='amoi_vs_span.png',
                     e_amoi_span_png_filename: PathLike='e_amoi_vs_span.png',
                     centroid_span_png_filename: PathLike='centroid_vs_span.png',
                     plot: bool=True,
                     show: bool=False) -> tuple[dict[str, np.ndarray],       # y, L, A, I, J, ExI, EyI, GJ, avg_centroid,
                                                list[str], list[str], int]:  # plane_bdf_filenames1, plane_bdf_filenames2, ifig
    """
    Cut a shell (and optionally bar/beam) FEM at prescribed stations and
    return the section stiffness distribution along the span.

    The cutting plane marches along the coord's local y-axis.  At each
    station, a cut in the local xz-plane intersects every shell face that
    straddles the plane, producing a ring of wall segments whose areas,
    centroids and moduli are integrated into the six independent second
    moments, Bredt-Batho torsion constant, and shear center.

    The primary second moments in the cut-plane frame are::

        Ixx = sum(A * x^2)      second moment about the z-axis
        Izz = sum(A * z^2)      second moment about the x-axis
        Ixz = sum(A * x * z)    product of inertia
        J   = Ixx + Izz         polar moment

    Terms involving the out-of-plane coordinate y are zero by construction.

    Composite laminates are handled via the ``[A]^{-1}`` equivalent
    modulus: each element's Ex (normal to the cut) and Ey (tangential)
    are resolved into the cut frame before multiplying by the element's
    area, so ``ExI`` and ``EyI`` carry the real bending stiffness even
    for multi-material sections.

    When ``include_lines=True``, CBAR/CBEAM elements that straddle the
    cut plane contribute concentrated area, bending stiffness (EA, EI),
    and torsion stiffness (GJ).  Their own bending inertia is rotated
    from the element's local axes into the cut-plane frame.  They do
    NOT participate in the thin-walled Bredt-Batho or shear-center
    solves, which remain shell-only.

    Parameters
    ----------
    bdf_filename : PathLike | BDF
        path to a bulk-data file, or an already-loaded cross-referenced
        ``BDF`` object
    normal_plane : (3,) float ndarray
        unit normal of the cutting plane in the basic frame; for a wing
        cut along +y this is ``[0, 1, 0]``
    log : SimpleLogger
        logging object
    stations : list[float] | (nstation,) ndarray
        the y-coordinates (in the basic frame) at which to cut
    coords : list[CORD2R]
        one coordinate system per station, centred on that station;
        its local y-axis is the march direction (normal to the cut
        plane) and its local xz-plane defines the in-plane axes for
        the section integrals
    x_vector : list[float]
        the CBEAM orientation vector *v*, in the basic frame (the emitted
        elements use ``offt='GGG'``).  This is not cosmetic: it sets the
        element y/z axes, and therefore which way round I1 and I2 come out
        and what sign I12 takes.  With a beam along +y, ``[1,0,0]`` puts I1
        on the chordwise moment while ``[0,0,1]`` puts I1 on the flapwise
        one.  Must not be parallel to the beam axis.
    include_lines : bool; default=False
        when True, CBAR and CBEAM elements that straddle each cut plane
        are included in the EA, EI and GJ totals.  Their own bending
        inertia (I1, I2, I12 from the PBAR/PBEAM) is rotated into the
        cut-plane frame and added on top of the parallel-axis A*d^2 term.
        Torsion uses the element's actual J, not a polar-moment
        approximation.  Bredt-Batho and shear-center remain shell-only.
    include_solids : bool; default=False
        unused, reserved for future solid-element support
    face_data : tuple | None; default=None
        pre-computed face topology from ``_setup_faces``; if None it is
        built automatically from the model.  Structure::

            (nids, xyz_cid0, elements)

        where *elements* is a dict keyed by ``'tri3'`` (etc.) whose
        values are ``(eids, node_ids, zoffset)`` tuples
    dirname : PathLike; default=''
        directory for all output files (CSV, BDF, PNG)
    ifig : int; default=1
        starting matplotlib figure number; useful when making multiple
        cuts so the plots do not overwrite each other
    debug_vectorize : bool; default=True
        use the faster vectorized cutting-plane method
    debug_v3 : bool; default=False
        use the experimental v3 cutting path
    rho : float; default=1.0
        density written on the equivalent beam model's MAT1; the PBEAM A
        field is the real geometric area, so ``rho * A`` gives a
        meaningful mass
    xyz_round : int | None; default=None
        decimal places to round GRID coordinates to
    area_round : int | None; default=None
        decimal places to round area values to
    inertia_round : int | None; default=None
        decimal places to round I / J values to
    beam_grid_xyz : (nstation, 3) float ndarray | None; default=None
        where to put the equivalent beam model's GRIDs.  By default a
        GRID is dropped on each cut's shear center.  Pass an explicit
        set of points (load control points, an existing loads-model
        grid, a straight reference axis) and the GRIDs go there instead,
        with the difference carried on the CBEAM WA/WB offsets so the
        elastic axis still runs through the real shear centers.  One row
        per station, in the same order as *stations*; rows for stations
        that fail to cut are dropped.
    beam_grid_ids : (nstation,) int ndarray | None; default=None
        GRID ids to use with *beam_grid_xyz*.  Reusing the ids from the
        source deck lets parts that share a point merge into a connected
        model.
    beam_id0 : int; default=1
        first CBEAM/PBEAM id (and the MAT1 id) in the equivalent beam
        model.  Only matters when several beam models are merged without
        renumbering.
    stop_on_failure : bool; default=False
        if True, raise on any station that cannot be cut; if False,
        skip it and continue to the next station
    cut_data_span_filename : PathLike; default='cut_data_vs_span.csv'
        CSV file with all section data vs. span; set to ``''`` to skip
    beam_model_bdf_filename : PathLike; default='equivalent_beam_model.bdf'
        punch-format BDF of the equivalent CBEAM model; set to ``''``
        to skip
    thetas_csv_filename : PathLike; default='thetas.csv'
        per-element material-angle diagnostic file
    normalized_inertia_png_filename : PathLike; default='normalized_inertia_vs_span.png'
        filename for the normalised inertia plot
    area_span_png_filename : PathLike; default='area_vs_span.png'
        filename for the area-vs-span plot
    amoi_span_png_filename : PathLike; default='amoi_vs_span.png'
        filename for the area-MOI-vs-span plot
    e_amoi_span_png_filename : PathLike; default='e_amoi_vs_span.png'
        filename for the E*I-vs-span plot
    centroid_span_png_filename : PathLike; default='centroid_vs_span.png'
        filename for the centroid-vs-span plot
    plot : bool; default=True
        generate matplotlib PNG plots
    show : bool; default=False
        call ``plt.show()`` after plotting (blocks until closed)

    Returns
    -------
    out_dict : dict[str, ndarray]
        section data keyed by name; the nine core entries are
        ``'stations'``, ``'L'``, ``'A'``, ``'I'``, ``'J'``,
        ``'ExI'``, ``'EyI'``, ``'GJ'``, ``'avg_centroid'``.
        Additional keys: ``'neutral_axis'``, ``'shear_center'``
        (both in the basic frame).
    plane_bdf_filenames1 : list[str]
        paths to the cut-plane face BDFs (side 1)
    plane_bdf_filenames2 : list[str]
        paths to the cut-plane face BDFs (side 2)
    ifig : int
        the next unused figure number
    """
    assert isinstance(x_vector, list), x_vector
    assert len(x_vector) == 3, x_vector

    if isinstance(dirname, str):
        dirname = Path(dirname)
    if isinstance(bdf_filename, PathLike):
        model = read_bdf(bdf_filename, log=log)
        model_static = copy.deepcopy(model)
    else:
        model = bdf_filename
        model_static = bdf_filename

    out = _get_station_data(
        model, model_static,
        stations, coords, normal_plane,
        dirname, face_data=face_data,
        include_lines=include_lines, include_solids=include_solids,
        debug_vectorize=debug_vectorize,
        debug_v3=debug_v3,
        stop_on_failure=stop_on_failure,
    )
    (thetas, stations, dx, dz, L, A, I, J, ExI, EyI, GJ, avg_centroid,
     plane_bdf_filenames, plane_bdf_filenames2, ExA, EyA, GA,
     avg_centroid_global, neutral_axis_global, shear_center_global,
     neutral_axis_offset) = out

    assert len(stations) > 0, stations
    thetas_csv_filename = dirname / thetas_csv_filename
    with open(thetas_csv_filename, 'w') as csv_filename:
        csv_filename.write('# eid(%d),theta,Ex,Ey,Gxy\n')
        for eid, (theta, ex, ey, gxy) in sorted(thetas.items()):
            csv_filename.write(f'{eid:d},{theta},{ex},{ey},{gxy}\n')

    avg_centroid[:, 1] = stations

    #   0    1    2    3    4    5
    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    Ix = I[:, 0]
    # Iy = I[:, 1]
    Iz = I[:, 2]
    Ixz = I[:, 5]

    ExIx = ExI[:, 0]
    # ExIy = ExI[:, 1]
    ExIz = ExI[:, 2]
    ExIxz = ExI[:, 5]
    # Ex = ExIx / Ix
    # Ey = ExIz / Iz
    J = Ix + Iz
    G = GJ / J
    #i1, i2, i12 = Ix, Iy, Ixy

    if beam_model_bdf_filename:
        beam_model_bdf_filename = dirname / beam_model_bdf_filename
        # Ex* rather than Ey* because Ex is the modulus along the beam axis
        # (normal to the cut plane); see the note in the docstring.
        #
        # *_global, NOT avg_centroid: the latter is in the cut coord's local
        # frame with column 1 overwritten by the station.  That happens to be
        # the basic frame for a wing cut (the coord is built so the local axes
        # coincide with the global ones), but for a fuselage cut the coord is
        # rotated and the GRIDs would come out permuted.
        #
        # ExIx/ExIz/ExIxz are second moments about the *cut coord's* local x
        # and z axes, so _write_beam_model needs to know where those axes
        # point in the basic frame to rotate them into the element frame.
        plane_i = np.array([coord.i for coord in coords], dtype='float64')
        plane_k = np.array([coord.k for coord in coords], dtype='float64')
        _write_beam_model(
            neutral_axis_global, shear_center_global, neutral_axis_offset,
            A, ExA, GA,
            ExIx, ExIz, ExIxz, GJ,
            x_vector=x_vector,
            plane_i=plane_i, plane_k=plane_k,
            bdf_filename=beam_model_bdf_filename,
            rho=rho, xyz_round=xyz_round, area_round=area_round,
            inertia_round=inertia_round,
            beam_grid_xyz=beam_grid_xyz, beam_grid_ids=beam_grid_ids,
            beam_id0=beam_id0, log=log)

    if cut_data_span_filename:
        cut_data_span_filename = dirname / cut_data_span_filename
        inotnan = np.isfinite(stations)
        X = np.column_stack([stations, dx, dz, A,
                             I, J, #Ex, Ey, G,
                             ExI,
                             EyI,
                             GJ,
                             avg_centroid])[inotnan, :]

        header = (
            'station,dx,dz,A,'
            'Ix,Iy,Iz,Ixy,Ixz,Iyz,J,'
            'Ex*Ix,Ex*Iy,Ex*Iz,Ex*Ixy,Ex*Ixz,Ex*Iyz,'
            'Ey*Ix,Ey*Iy,Ey*Iz,Ey*Ixy,Ey*Ixz,Ey*Iyz,'
            'GJ,'
            'xcentroid,ycentroid,zcentroid')
        np.savetxt(cut_data_span_filename, X, header=header, delimiter=',')

    if plot:
        ifig = plot_inertia(
            log, stations, A, I, J,
            ExI, EyI, GJ, avg_centroid, show=show,
            dirname=dirname, ifig=ifig,
            normalized_inertia_png_filename=normalized_inertia_png_filename,
            amoi_span_png_filename=amoi_span_png_filename,
            e_amoi_span_png_filename=e_amoi_span_png_filename,
            centroid_span_png_filename=centroid_span_png_filename,
        )

    out_dict = {
        'stations': stations, 'L': L, 'A': A, 'I': I, 'J': J,
        'ExI': ExI, 'EyI': EyI, 'GJ': GJ, 'avg_centroid': avg_centroid,
        # basic frame; 'avg_centroid' is area-weighted and in the cut frame
        'neutral_axis': neutral_axis_global,
        'shear_center': shear_center_global,
    }
    return out_dict, plane_bdf_filenames, plane_bdf_filenames2, ifig


def load_moi_data(csv_filename: PathLike) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                                                   np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a ``cut_data_vs_span.csv`` written by ``cut_and_plot_moi`` and
    return the same arrays it produces in memory.

    Parameters
    ----------
    csv_filename : PathLike
        path to the CSV file

    Returns
    -------
    y : (nstation,) float ndarray
        span stations
    A : (nstation,) float ndarray
        section area at each station
    I : (nstation, 6) float ndarray
        second moments ``[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]``
    J : (nstation,) float ndarray
        polar moment ``Ixx + Izz``
    ExI : (nstation, 6) float ndarray
        modulus-weighted second moments (Ex * I)
    EyI : (nstation, 6) float ndarray
        modulus-weighted second moments (Ey * I)
    GJ : (nstation,) float ndarray
        torsion stiffness
    avg_centroid : (nstation, 3) float ndarray
        section centroid ``(x, y, z)``
    """
    import pandas as pd
    # 'station, dx, dz, A, '
    # 'Ix, Iy, Iz, Ixy, Ixz, Iyz, J, '
    # 'Ex*Ix, Ex*Iy, Ex*Iz, Ex*Ixy, Ex*Ixz, Ix*Iyz, '
    # 'Ey*Ix, Ey*Iy, Ey*Iz, Ey*Ixy, Ey*Ixz, Iy*Iyz, '
    # 'GJ,'
    # 'xcentroid, ycentroid, zcentroid')
    df = pd.read_csv(csv_filename)
    df.columns = df.columns.str.strip(' #')
    # print(df.columns)
    y = df['station'].to_numpy()
    A = df['A'].to_numpy()

    # The column *names* are local-frame labels and carry no global meaning --
    # on a swept cut there is no right answer for which in-plane product is
    # "Ixz" vs "Iyz", so they are not worth arguing about.  What does matter is
    # that this reads back in the same order `cut_and_plot_moi` wrote, so that a
    # CSV round-trip is the identity and `I[:, 5]` still means what it meant in
    # memory.  These lists must therefore match the header exactly.
    #
    # They used to end '..., Iyz, Ixz' against a header ending '..., Ixz, Iyz',
    # which permuted the last two columns on load.  `plot_compare_inertia` takes
    # `Ixz = I[:, 5]` straight off this, so it plotted the ~1e-14 out-of-plane
    # term instead of the real one (2.4e3 on the wing) -- a flat-zero curve.
    I = df[['Ix', 'Iy', 'Iz', 'Ixy', 'Ixz', 'Iyz']].to_numpy()
    ExI = df[['Ex*Ix', 'Ex*Iy', 'Ex*Iz', 'Ex*Ixy', 'Ex*Ixz', 'Ex*Iyz']].to_numpy()
    EyI = df[['Ey*Ix', 'Ey*Iy', 'Ey*Iz', 'Ey*Ixy', 'Ey*Ixz', 'Ey*Iyz']].to_numpy()
    GJ = df['GJ'].to_numpy()
    J = df['J'].to_numpy()
    avg_centroid = df[['xcentroid', 'ycentroid', 'zcentroid']].to_numpy()
    return y, A, I, J, ExI, EyI, GJ, avg_centroid


def _element_triad(xyz_a: np.ndarray,
                   xyz_b: np.ndarray,
                   x_vector: np.ndarray,
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    The MSC CBEAM element triad for ``offt='GGG'`` (orientation vector in the
    basic frame)::

        x_e  along GA+WA -> GB+WB
        y_e  the part of v normal to x_e   (plane 1 contains x_e and v)
        z_e  = x_e cross y_e

    Raises if v is parallel to the element axis, which leaves plane 1
    undefined -- Nastran would reject the CBEAM for the same reason.
    """
    dxyz = xyz_b - xyz_a
    norm = np.linalg.norm(dxyz)
    if norm == 0.0:
        raise ValueError(f'coincident beam ends {xyz_a} and {xyz_b}')
    x_e = dxyz / norm

    y_e = x_vector - x_vector.dot(x_e) * x_e
    norm_y = np.linalg.norm(y_e)
    if norm_y < 1e-8:
        raise ValueError(
            f'x_vector={x_vector} is parallel to the beam axis {x_e}; '
            'the CBEAM orientation plane is undefined')
    y_e /= norm_y
    z_e = np.cross(x_e, y_e)
    return x_e, y_e, z_e


def _project_section(ixx: float, izz: float, ixz: float,
                     y_e: np.ndarray, z_e: np.ndarray,
                     plane_i: np.ndarray, plane_k: np.ndarray,
                     ) -> tuple[float, float, float]:
    """
    Rotates the cut-plane second moments onto the element axes.

    ``ixx``, ``izz``, ``ixz`` are ``int(x_loc^2)``, ``int(z_loc^2)`` and
    ``int(x_loc*z_loc)`` about the cut coord's local axes ``plane_i`` and
    ``plane_k``.  They are the three independent entries of the symmetric
    2x2 tensor ``M`` in that basis, so for any in-plane direction pair::

        I1  = a.M.a      I2 = b.M.b      I12 = a.M.b

    with ``a`` and ``b`` the components of ``y_e`` and ``z_e`` along
    ``(plane_i, plane_k)``.  The element axes are normal to ``x_e``, not to
    the cut plane, so for a swept beam they poke slightly out of plane;
    the in-plane parts are renormalized, which is the closest the cut-plane
    integrals can get.

    Returns ``(I1, I2, I12)``.
    """
    def _in_plane(vec: np.ndarray) -> tuple[float, float]:
        c = vec.dot(plane_i)
        s = vec.dot(plane_k)
        norm = np.hypot(c, s)
        if norm < 1e-8:
            raise ValueError(
                f'element axis {vec} is normal to the cut plane spanned by '
                f'{plane_i} and {plane_k}; the section cannot be projected')
        return c / norm, s / norm

    cy, sy = _in_plane(y_e)
    cz, sz = _in_plane(z_e)

    i1 = cy * cy * ixx + 2. * cy * sy * ixz + sy * sy * izz
    i2 = cz * cz * ixx + 2. * cz * sz * ixz + sz * sz * izz
    i12 = cy * cz * ixx + (cy * sz + sy * cz) * ixz + sy * sz * izz
    return i1, i2, i12


def _write_beam_model(neutral_axis: np.ndarray,
                      shear_center: np.ndarray,
                      neutral_axis_offset: np.ndarray,
                      A: np.ndarray,
                      ExA: np.ndarray,
                      GA: np.ndarray,
                      ExIxx: np.ndarray,
                      ExIzz: np.ndarray,
                      ExIxz: np.ndarray,
                      GJ: np.ndarray,
                      x_vector: list[float],
                      plane_i: np.ndarray,
                      plane_k: np.ndarray,
                      bdf_filename: PathLike='',
                      rho: float=0.1,
                      xyz_round: int | None=None,
                      area_round: int | None=None,
                      inertia_round: int | None=None,
                      beam_grid_xyz: Optional[np.ndarray]=None,
                      beam_grid_ids: Optional[np.ndarray]=None,
                      beam_id0: int=1,
                      log=None,
                      ):
    """
    Assume y is down the axis of the beam
    The beam cross section is defined in the x-z plane.
    z
    |
    |
    |
    ----------> x

    Parameters
    ----------
    neutral_axis : (nstation, 3) float ndarray
        the modulus-weighted centroid of each cut, in the basic frame
    shear_center : (nstation, 3) float ndarray
        where a transverse shear produces no twist, in the basic frame;
        NaN for any station where it could not be found
    neutral_axis_offset : (nstation, 2) float ndarray
        neutral axis minus *area* centroid, in the CUT frame's in-plane axes
        ``(plane_i, plane_k)``.  ``ExIxx``/``ExIzz``/``ExIxz`` are reported
        about the area centroid, so this is what moves them to the neutral
        axis; it is identically zero for a homogeneous section.

    The section stiffnesses that come out of the cut are modulus-weighted
    integrals (EA, E*I, G*J), so they have to be split into a material part
    and a property part before they can be written.  Rather than write a
    unit-modulus MAT1 and fold the moduli into the PBEAM (which couples EA to
    the transverse shear term K*G*A, because both read the same A field), the
    MAT1 carries real reference moduli::

        E_ref = sum(Ex_i*dA_i) / sum(dA_i)      area-weighted axial modulus
        G_ref = sum(Gxy_i*dA_i) / sum(dA_i)     area-weighted shear modulus
        nu    = E_ref/(2*G_ref) - 1             consistent by construction

    and the PBEAM carries effective section properties referenced to them::

        A   = EA / E_ref        -> E_ref*A   == EA      (exact)
        I1  = E*Iyy_e / E_ref   -> E_ref*I1  == E*Iyy_e (exact)
        I2  = E*Izz_e / E_ref   -> E_ref*I2  == E*Izz_e (exact)
        I12 = E*Iyz_e / E_ref
        J   = GJ / G_ref        -> G_ref*J   == GJ      (exact)
        K1  = K2 = GA / (G_ref*A)                       (exact)

    Every stiffness is then reproduced independently.  For a homogeneous
    section this degenerates to the physical answer: E_ref/G_ref are the real
    moduli, A is the real area (so ``rho`` gives a meaningful mass), and
    K1 = K2 = 1.

    **The element axis rides the shear center, not the centroid.**  MSC's
    Figure 16-142 lays out three distinct lines through the section::

        GA --w_a--> element origin (0,0,0)  ...on the SHEAR CENTER line
                    + (N1(A), N2(A))        ...to the NEUTRAL AXIS line
                    + (M1(A), M2(A))        ...to the nonstructural mass

    so ``GA+WA`` is the point about which a transverse shear produces no
    twist, and ``N1``/``N2`` say where the modulus-weighted centroid sits
    relative to it, in the *element* frame (``N1`` along ``y_e``, ``N2``
    along ``z_e``).  I1/I2/I12/J are taken about the neutral axis.

    Putting the axis on the centroid instead and writing ``N1=N2=0`` asserts
    centroid == shear center, which is true for a closed doubly symmetric
    section (a fuselage barrel, a boom) and wrong for anything with a single
    plane of symmetry or none -- a wing box, a fin, a stabilizer.  There the
    error is not small: the shear center of an airfoil box sits well forward
    of the centroid, and collapsing them silently couples bending into
    torsion.  ``shear_center`` solves the transverse shear flow for the real
    location; where it cannot (a degenerate or collinear cut) it returns NaN
    and the axis falls back to the neutral axis with ``N1=N2=0``, which is
    the old behavior.

    By default a GRID is written on each cut's shear center.  Passing
    ``beam_grid_xyz`` moves the GRIDs onto a prescribed set of points (load
    control points, for instance) and puts the difference on the CBEAM WA/WB
    offset vectors::

        WA = shear_center[i]   - beam_grid_xyz[i]
        WB = shear_center[i+1] - beam_grid_xyz[i+1]

    A CBEAM's elastic axis runs from GA+WA to GB+WB, so the offsets restore
    exactly the geometry the un-offset model would have had -- the section
    properties are unchanged and still referenced to the neutral axis, only
    the nodes have moved.  ``offt='GGG'`` puts the offsets in the global
    frame, which is what the section points are already in.

    **The second moments are re-referenced to the neutral axis.**  The cutter
    reports ``ExI`` about the *area* centroid; for a heterogeneous section
    the neutral axis is somewhere else, and the parallel-axis theorem moves
    it::

        ExIxx_na = ExIxx - ExA*du**2       du, dv = neutral_axis_offset
        ExIzz_na = ExIzz - ExA*dv**2
        ExIxz_na = ExIxz - ExA*du*dv

    The shift is toward the centroid of the *stiffness*, so the corrections
    are subtractive and vanish identically for a homogeneous section.

    **Section properties are rotated into the element frame.**  The cutter
    reports second moments about the *cut coord's* in-plane axes::

        Ixx = int(x_loc^2 dA)   Izz = int(z_loc^2 dA)   Ixz = int(x_loc*z_loc dA)

    where ``x_loc``/``z_loc`` are the coord's local x and z (``plane_i`` and
    ``plane_k`` give those axes in the basic frame).  MSC wants them about the
    *element* axes instead::

        I1 = int(y_e^2 dA)   I2 = int(z_e^2 dA)   I12 = int(y_e*z_e dA)

    and the element triad follows from the CBEAM orientation vector::

        x_e = (GB+WB - GA-WA) / |...|
        y_e = (v - (v.x_e) x_e) / |...|        v = ``x_vector``, offt='GGG'
        z_e = x_e cross y_e

    Those two bases differ by a rotation in the cut plane, so the three
    second moments transform as a 2x2 tensor rather than mapping one-to-one.
    Writing ``I1 = Ixx`` and ``I12 = +Ixz`` (what this used to do) is only
    right when ``y_e`` happens to land on ``+x_loc``; with ``v = [0, 0, 1]``
    and a spanwise beam it lands on ``z_loc`` instead, which swapped I1/I2
    and flipped the sign of I12.  ``_project_section`` does the rotation.

    TODO: the section is integrated in the cut plane, which is normal to the
          coord, not normal to the element axis.  For a swept/canted beam
          (vtail) those differ and the properties are a cosine off; a warning
          is emitted when the misalignment exceeds ~15 deg.
    """
    if isinstance(bdf_filename, str) and len(bdf_filename) == 0:
        return

    nstation = len(A)
    x_vector = np.asarray(x_vector, dtype='float64')
    plane_i = np.asarray(plane_i, dtype='float64')
    plane_k = np.asarray(plane_k, dtype='float64')
    neutral_axis = np.asarray(neutral_axis, dtype='float64')
    shear_center = np.asarray(shear_center, dtype='float64')
    neutral_axis_offset = np.asarray(neutral_axis_offset, dtype='float64')
    for nm, arr in (('plane_i', plane_i), ('plane_k', plane_k),
                    ('neutral_axis', neutral_axis),
                    ('shear_center', shear_center)):
        if arr.shape != (nstation, 3):
            raise ValueError(
                f'{nm} must be ({nstation:d}, 3) to match the stations; '
                f'got {arr.shape}')
    if neutral_axis_offset.shape != (nstation, 2):
        raise ValueError(
            f'neutral_axis_offset must be ({nstation:d}, 2) to match the '
            f'stations; got {neutral_axis_offset.shape}')
    if beam_grid_xyz is not None:
        beam_grid_xyz = np.asarray(beam_grid_xyz, dtype='float64')
        if beam_grid_xyz.shape != (nstation, 3):
            raise ValueError(
                f'beam_grid_xyz must be ({nstation:d}, 3) to match the '
                f'stations; got {beam_grid_xyz.shape}')
        if not np.isfinite(beam_grid_xyz).all():
            raise ValueError('beam_grid_xyz contains NaN/inf')
    if beam_grid_ids is not None:
        if beam_grid_xyz is None:
            raise ValueError('beam_grid_ids requires beam_grid_xyz')
        beam_grid_ids = np.asarray(beam_grid_ids, dtype='int64')
        if beam_grid_ids.shape != (nstation,):
            raise ValueError(
                f'beam_grid_ids must be ({nstation:d},) to match the '
                f'stations; got {beam_grid_ids.shape}')
        if len(np.unique(beam_grid_ids)) != nstation:
            raise ValueError('beam_grid_ids are not unique')

    # stations where no cut was found come back as NaN; writing them produces
    # blank GRID/PBEAM fields and an unreadable deck
    ivalid = np.where(
        np.isfinite(A) & np.isfinite(ExA) & np.isfinite(GA) &
        np.isfinite(neutral_axis).all(axis=1))[0]
    if len(ivalid) < 2:
        raise RuntimeError(
            f'cannot write an equivalent beam model; only {len(ivalid):d} '
            'valid station(s) were cut (2 are needed to make a CBEAM)')
    if beam_grid_xyz is not None and len(ivalid) != nstation and log is not None:
        # a prescribed grid point is usually there because something else
        # attaches to it, so silently dropping one is worth a shout
        idropped = np.setdiff1d(np.arange(nstation), ivalid)
        dropped = (beam_grid_ids[idropped].tolist()
                   if beam_grid_ids is not None else idropped.tolist())
        log.warning(
            f'{len(idropped):d} of {nstation:d} prescribed beam grid points '
            f'had no cut and were dropped: {dropped}. If any of them is an '
            'interface point, the merged model will not connect there.')

    neutral_axis = neutral_axis[ivalid, :]
    shear_center = shear_center[ivalid, :]
    neutral_axis_offset = neutral_axis_offset[ivalid, :]
    if beam_grid_xyz is not None:
        beam_grid_xyz = beam_grid_xyz[ivalid, :]
    if beam_grid_ids is not None:
        beam_grid_ids = beam_grid_ids[ivalid]
    A = A[ivalid]
    ExA = ExA[ivalid]
    GA = GA[ivalid]
    ExIxx = ExIxx[ivalid]
    ExIzz = ExIzz[ivalid]
    ExIxz = ExIxz[ivalid]
    GJ = GJ[ivalid]
    plane_i = plane_i[ivalid, :]
    plane_k = plane_k[ivalid, :]

    # The element axis goes on the shear center.  Where the shear flow solve
    # could not place one, fall back to the neutral axis -- that is only right
    # for a doubly symmetric section, so say so rather than quietly writing a
    # beam whose torsion is referenced to the wrong line.
    nofound = ~np.isfinite(shear_center).all(axis=1)
    axis_xyz = np.where(nofound[:, np.newaxis], neutral_axis, shear_center)
    if nofound.any() and log is not None:
        log.warning(
            f'{int(nofound.sum()):d} of {len(ivalid):d} stations have no shear '
            'center; the beam axis there falls back to the neutral axis and '
            'N1/N2 are written as 0.')

    # Parallel axis: ExI comes back about the AREA centroid, but I1/I2/I12 are
    # defined about the neutral axis.  The two coincide for a homogeneous
    # section, so this is exactly zero there.
    du = neutral_axis_offset[:, 0]
    dv = neutral_axis_offset[:, 1]
    ExIxx = ExIxx - ExA * du ** 2
    ExIzz = ExIzz - ExA * dv ** 2
    ExIxz = ExIxz - ExA * du * dv

    # area-weighted reference moduli; a single MAT1 for the whole beam
    Atotal = A.sum()
    E_ref = ExA.sum() / Atotal
    G_ref = GA.sum() / Atotal
    nu = E_ref / (2. * G_ref) - 1.

    # effective section properties referenced to E_ref/G_ref.  These are still
    # in the cut-plane basis; the rotation into the element frame happens per
    # element, because that is where the element axis is known.
    area_eff = ExA / E_ref
    ixx_eff = ExIxx / E_ref
    izz_eff = ExIzz / E_ref
    ixz_eff = ExIxz / E_ref
    j_eff = GJ / G_ref
    # shear correction factor; 1.0 when E/G is uniform over the section
    k_eff = GA / (G_ref * area_eff)

    mid = beam_id0
    beam_model = BDF(debug=False)
    beam_model.add_mat1(mid=mid, E=E_ref, G=G_ref, nu=nu, rho=rho)

    # where the GRIDs go, and how far that is from the shear center
    if beam_grid_xyz is None:
        grid_xyz = axis_xyz
        offset = None
        if xyz_round is not None:
            grid_xyz = grid_xyz.round(xyz_round)
    else:
        # prescribed points are left exactly as supplied -- rounding them
        # would break the coincidence with whatever deck they came from --
        # so xyz_round is applied to the offsets instead
        grid_xyz = beam_grid_xyz
        offset = axis_xyz - beam_grid_xyz
        if xyz_round is not None:
            offset = offset.round(xyz_round)

    nid = (np.arange(1, len(grid_xyz) + 1) if beam_grid_ids is None
           else beam_grid_ids)
    if area_round is not None:
        area_eff = area_eff.round(area_round)
    if inertia_round is not None:
        j_eff = j_eff.round(inertia_round)

    for nidi, xyz in zip(nid, grid_xyz):
        beam_model.add_grid(int(nidi), xyz)

    # where the element axis really ends up once WA/WB have been rounded; N1/N2
    # are measured from that, not from the unrounded shear center
    axis_written = grid_xyz if offset is None else grid_xyz + offset

    # the element axis is only normal to the cut plane for an unswept beam;
    # warn once rather than per element
    skewed = 0
    for ielem in range(1, len(A)):
        eid = pid = beam_id0 + ielem - 1
        nids = [int(nid[ielem-1]), int(nid[ielem])]
        g0 = None
        if offset is None:
            wa = wb = None
        else:
            # GA+WA -> GB+WB is the elastic axis, so this puts the element
            # back on the shear center line no matter where the GRIDs sit
            wa = offset[ielem-1, :].tolist()
            wb = offset[ielem, :].tolist()
        beam_model.add_cbeam(eid, pid, nids, x_vector.tolist(), g0,
                             offt='GGG', bit=None,
                             pa=0, pb=0, wa=wa, wb=wb, sa=0, sb=0, comment='')
        so = ['YES', 'YES']
        xxb = [0., 1.]
        area = [area_eff[ielem-1], area_eff[ielem]]

        # rotate each end's cut-plane second moments onto this element's axes
        xyz_a = axis_written[ielem-1, :]
        xyz_b = axis_written[ielem, :]
        x_e, y_e, z_e = _element_triad(xyz_a, xyz_b, x_vector)
        i1, i2, i12, n1, n2 = [], [], [], [], []
        for iend in (ielem-1, ielem):
            normal = np.cross(plane_i[iend, :], plane_k[iend, :])
            if abs(x_e.dot(normal)) < 0.966:  # ~15 deg
                skewed += 1
            i1i, i2i, i12i = _project_section(
                ixx_eff[iend], izz_eff[iend], ixz_eff[iend],
                y_e, z_e, plane_i[iend, :], plane_k[iend, :])
            if inertia_round is not None:
                i1i = round(i1i, inertia_round)
                i2i = round(i2i, inertia_round)
                i12i = round(i12i, inertia_round)
            i1.append(i1i)
            i2.append(i2i)
            i12.append(i12i)

            # N1/N2 locate the neutral axis relative to the element axis (which
            # is on the shear center), resolved on the element's own y_e/z_e.
            # The offset is in the cut plane and the element axis is not quite
            # normal to it on a swept beam, so a sliver of it falls along x_e
            # and is simply dropped -- there is no PBEAM field for it.
            dxyz = neutral_axis[iend, :] - axis_written[iend, :]
            n1i = float(dxyz.dot(y_e))
            n2i = float(dxyz.dot(z_e))
            if xyz_round is not None:
                n1i = round(n1i, xyz_round)
                n2i = round(n2i, xyz_round)
            n1.append(n1i)
            n2.append(n2i)
        j = [j_eff[ielem-1], j_eff[ielem]]
        # K is constant over the element; average the two ends
        k1 = k2 = 0.5 * (k_eff[ielem-1] + k_eff[ielem])
        beam_model.add_pbeam(pid, mid, xxb, so, area, i1, i2, i12, j, nsm=None,
                             c1=None, c2=None, d1=None, d2=None, e1=None, e2=None, f1=None, f2=None,
                             k1=k1, k2=k2, s1=0., s2=0., nsia=0., nsib=None, cwa=0., cwb=None,
                             m1a=0., m2a=0., m1b=None, m2b=None,
                             n1a=n1[0], n2a=n2[0], n1b=n1[1], n2b=n2[1],
                             comment='')
    if skewed and log is not None:
        log.warning(
            f'{skewed:d} of {2*(len(A)-1):d} beam ends have the element axis '
            'more than ~15 deg off the cut-plane normal (swept/canted beam). '
            'The section was integrated in the cut plane, so I1/I2/I12/A are '
            'overstated by roughly 1/cos(angle).')
    beam_model.write_bdf(bdf_filename)


def _get_station_data(model: BDF,
                      model_static: BDF,
                      dys: list[float],
                      coords: list[CORD2R],
                      normal_plane: np.ndarray,
                      dirname: Path,
                      plane_atol: float=1e-5,
                      include_lines: bool=False,
                      include_solids: bool=False,
                      debug_vectorize: bool=True,
                      debug_v3: bool=False,
                      stop_on_failure: bool=False,
                      face_data=None) -> tuple[
                         dict[int, tuple[float, float, float, float]],  # thetas
                         #y, dx, dz,
                         Any, Any, Any,
                         #A, I, J,
                         Any, Any, Any,
                         #ExI, EyI, GJ,
                         Any, Any, Any,
                         #avg_centroid, plane_bdf_filenames, plane_bdf_filenames2,
                         Any, list[str], list[str]]:
    """
    Loop over stations, cut the model at each one, and accumulate all
    section properties into span-length arrays.

    This is the inner workhorse of ``cut_and_plot_moi``; it handles the
    station loop, the optional bar/beam crossing detection, and the call
    to ``calculate_area_moi`` for each successful cut.

    Parameters
    ----------
    model : BDF
        cross-referenced model used for element property look-ups and
        bar/beam crossing detection
    model_static : BDF
        a (possibly deep-copied) model used by the face-cutting
        geometry routines, so that adding temporary coords to *model*
        does not mutate the caller's object
    dys : list[float]
        y-stations (in the basic frame) at which to cut
    coords : list[CORD2R]
        one coordinate system per station; cutting plane is its local
        xz-plane (y_local = 0)
    normal_plane : (3,) float ndarray
        unit normal of the cutting planes in the basic frame
    dirname : Path
        base directory for intermediate BDF / CSV files
    plane_atol : float; default=1e-5
        absolute tolerance for the cutting-plane intersection
    include_lines : bool; default=False
        find CBAR/CBEAM elements straddling each cut and pass them to
        ``calculate_area_moi`` for inclusion in EA, EI and GJ
    include_solids : bool; default=False
        unused
    debug_vectorize : bool; default=True
        use the faster vectorized cutting-plane method
    debug_v3 : bool; default=False
        use the experimental v3 cutting path
    stop_on_failure : bool; default=False
        if True, raise when a station cannot be cut; if False, skip it
    face_data : tuple | None; default=None
        pre-computed face topology; built automatically when None

    Returns
    -------
    tuple
        ``(thetas, y, dx, dz, L, A, I, J, ExI, EyI, GJ,
        avg_centroid, plane_bdf_filenames1, plane_bdf_filenames2,
        ExA, EyA, GA, avg_centroid_global, neutral_axis_global,
        shear_center_global, neutral_axis_offset)``
    """
    log = model.log

    # initialize theta
    thetas = {}
    for eid in model.elements:
        #  theta, Ex, Ey, Gxy
        thetas[eid] = (0., 0., 0., 0.)

    if face_data is None:
        # TODO: could filter out unused nodes
        _log, *face_data = _setup_faces(
            model,
            include_lines=include_lines, include_solids=include_solids)
    nodes, xyz_cid0, elements = face_data
    tri_eids, tri_nodes, _zoffset = elements['tri3']
    # nnode = len(nodes)
    # ntri = len(tri_eids)

    #p1 = np.array([466.78845, 735.9053, 0.0])
    #p2 = np.array([624.91345, 639.68896, -0.99763656])
    #dx = p2 - p1
    plane_bdf_filenames1 = []
    plane_bdf_filenames2 = []

    ny = len(dys)
    assert ny > 0, dys
    y = np.full(ny, np.nan, dtype='float64')
    dx = np.full(ny, np.nan, dtype='float64')
    dz = np.full(ny, np.nan, dtype='float64')
    length = np.full(ny, np.nan, dtype='float64')
    area = np.full(ny, np.nan, dtype='float64')
    inertia = np.full((ny, 6), np.nan, dtype='float64')
    J = np.full(ny, np.nan, dtype='float64')
    ExI = np.full((ny, 6), np.nan, dtype='float64')
    EyI = np.full((ny, 6), np.nan, dtype='float64')
    GJ = np.full(ny, np.nan, dtype='float64')
    avg_centroid = np.full((ny, 3), np.nan, dtype='float64')
    ExA = np.full(ny, np.nan, dtype='float64')
    EyA = np.full(ny, np.nan, dtype='float64')
    GA = np.full(ny, np.nan, dtype='float64')
    # avg_centroid is reported in the cut coord's LOCAL frame (the plots and
    # the csv want in-plane coordinates), so the beam model needs its own copy
    # transformed back to the basic frame or the GRIDs come out rotated.
    avg_centroid_global = np.full((ny, 3), np.nan, dtype='float64')
    # the neutral axis (modulus-weighted centroid) and the shear center, both
    # in the basic frame.  These are what the CBEAM actually needs: the element
    # axis rides the shear center and N1/N2 point from it to the neutral axis.
    neutral_axis_global = np.full((ny, 3), np.nan, dtype='float64')
    shear_center_global = np.full((ny, 3), np.nan, dtype='float64')
    # ...and the in-plane offset from the area centroid to the neutral axis,
    # kept in the CUT frame because that is the frame ExI is reported in
    neutral_axis_offset = np.zeros((ny, 2), dtype='float64')

    log.debug(f'dys={dys}; n={len(dys):d}')
    assert len(dys) == len(coords), (len(dys), len(coords))

    # Lazily-populated cache for shell material properties.  The expensive
    # material_coordinate_system / get_Ainv_equivalent_pshell calls depend
    # only on the element, not the station, so each element is computed at
    # most once.  An empty dict signals "use and populate the cache"; the
    # get_element_inertias fast path fills it on first encounter.
    shell_prop_cache: dict[int, tuple[float, float, float, float, float]] = {}

    ncuts_found = 0
    for icut, dy, coord in zip(count(), dys, coords):
        # itri_nodes = np.searchsorted(nodes, tri_nodes)
        # xyz_cid = coord.transform_node_to_local_array(xyz_cid0)
        # y_cid = xyz_cid[:, 1]
        # is_tri_cut = fis_tri_cut(y_cid, itri_nodes, ntri)

        model.coords[coord.cid] = coord
        plane_bdf_filename1 = dirname / f'plane_face1_{icut:d}.bdf'
        plane_bdf_filename2 = dirname / f'plane_face2_{icut:d}.bdf'
        cut_face_filename = dirname / f'cut_face_{icut:d}.csv'
        if os.path.exists(cut_face_filename):
            os.remove(cut_face_filename)

        found_cut, rods = _get_station_datai(
            model, model_static,
            dy, coord, plane_atol=plane_atol,
            debug_vectorize=debug_vectorize,
            debug_v3=debug_v3,
            stop_on_failure=stop_on_failure,
            plane_bdf_filename1=plane_bdf_filename1,
            plane_bdf_filename2=plane_bdf_filename2,
            face_data=face_data, log=log)

        # if not os.path.exists(plane_bdf_filename1) or len(rods) == 0:
        if not found_cut:
            log.debug(coord)
            log.debug(f'skipping calculate_area_moi {icut:d} (station={dy:g})')
            continue
            # break
        plane_bdf_filenames1.append(plane_bdf_filename1)
        plane_bdf_filenames2.append(plane_bdf_filename2)
        # eid, nid, inid1, inid2
        #print(unique_geometry_array)
        #moi_filename = 'amoi_%i.bdf' % i
        moi_filename = None
        log.info(f'calculate_area_moi {icut:d} (station={dy})')
        bar_data = None
        if include_lines:
            bar_data = _find_bar_beam_crossings(model, coord, log=log)
        (dxi, dzi, lengthi, areai,
         inertiai, Ji,
         ExIi, EyIi, GJi, avg_centroidi,
         ExAi, EyAi, GAi,
         neutral_axisi, shear_centeri) = calculate_area_moi(
            model, rods, normal_plane, thetas,
            moi_filename=moi_filename,
            bar_data=bar_data,
            shell_prop_cache=shell_prop_cache)

        #print(out)
        y[icut] = dy
        dx[icut] = dxi  # length
        dz[icut] = dzi  # height
        length[icut] = lengthi
        area[icut] = areai
        inertia[icut, :] = inertiai
        # print(Ji, EIi, GJi)
        # print(len(Ji), len(EIi), len(GJi))
        J[icut] = Ji
        ExI[icut, :] = ExIi
        EyI[icut, :] = EyIi
        GJ[icut] = GJi
        avg_centroid[icut, :] = avg_centroidi
        # the cut plane passes through the coord origin, so the local
        # out-of-plane component is ~0 and this lands on the real section
        # centroid in basic coordinates
        avg_centroid_global[icut, :] = coord.transform_node_to_global(avg_centroidi)
        neutral_axis_global[icut, :] = coord.transform_node_to_global(neutral_axisi)
        neutral_axis_offset[icut, :] = (neutral_axisi[[0, 2]] -
                                        avg_centroidi[[0, 2]])
        if np.isfinite(shear_centeri).all():
            shear_center_global[icut, :] = coord.transform_node_to_global(
                shear_centeri)
        ExA[icut] = ExAi
        EyA[icut] = EyAi
        GA[icut] = GAi
        ncuts_found += 1
        #break
    if ncuts_found == 0:
        raise RuntimeError('no cuts found...')

    out = (
        thetas, y, dx, dz,
        length, area, inertia, J,
        ExI, EyI, GJ,
        avg_centroid, plane_bdf_filenames1, plane_bdf_filenames2,
        ExA, EyA, GA, avg_centroid_global,
        neutral_axis_global, shear_center_global, neutral_axis_offset,
    )
    return out


def _get_station_datai(model: BDF,
                       model_static: BDF,
                       dy: float,
                       coord: CORD2R,
                       plane_atol: float=1e-5,
                       debug_vectorize: bool=True,
                       debug_v3: bool=False,
                       stop_on_failure: bool=False,
                       plane_bdf_filename1: PathLike='',
                       plane_bdf_filename2: PathLike='',
                       face_data=None,
                       log=None):
    """
    Cut the model at a single station and return the rod connectivity.

    Wraps ``cut_face_model_by_coord`` with error handling: when
    *stop_on_failure* is False, a station that cannot be cut returns
    ``(False, [])`` instead of raising.

    Returns
    -------
    found_cut : bool
        whether the cutting plane intersected any element faces
    rods : tuple
        ``(rod_eid_nodes, rod_nids, rod_xyzs)`` describing the wall
        segments that make up the cut ring; empty when *found_cut*
        is False
    """
    nodal_result = None
    try:
        out = cut_face_model_by_coord(
            model_static, coord,
            nodal_result, plane_atol=plane_atol,
            skip_cleanup=True,
            # csv_filename=cut_face_filename,
            csv_filename='',
            # plane_bdf_filename='')
            plane_bdf_filename1=plane_bdf_filename1,
            plane_bdf_filename2=plane_bdf_filename2,
            plane_y_offset=dy, face_data=face_data,
            debug_vectorize=debug_vectorize,
            debug_v3=debug_v3,
            stop_on_failure=stop_on_failure,
        )
    except PermissionError:
        print(f'failed to delete {plane_bdf_filename1}')
        raise
        # continue
    except RuntimeError as error:
        # incorrect ivalues=[0, 1, 2]; dy=771. for CRM
        #
        # A station that lands off the end of the structure (or exactly on the
        # last ring of nodes) legitimately has nothing to cut, and the cutter
        # signals that by raising rather than returning found_cut=False.  That
        # used to abort the whole run, which is the opposite of what
        # stop_on_failure=False asks for -- and it makes prescribed stations
        # (LCPs, which often sit at or just past a tip) unusable.  Honor the
        # flag: re-raise when the caller said the cut must succeed, otherwise
        # report the station as empty and march on.
        if stop_on_failure:
            raise
        if log is not None:
            log.warning(f'no cut at station={dy:g} (coord {coord.cid:d}): '
                        f'{error}')
        return False, []
    found_cut, unused_unique_geometry_array, unused_unique_results_array, rods = out
    return found_cut, rods


def plot_inertia(log: SimpleLogger,
                 station: np.ndarray, A: np.ndarray,
                 I: np.ndarray, J: np.ndarray,
                 ExI: np.ndarray, EyI: np.ndarray, GJ: np.ndarray,
                 avg_centroid: np.ndarray,
                 linestyle: str='-',
                 dirname: PathLike='',
                 x: str='x',
                 y: str='y',
                 z: str='z',
                 station_word='Span',
                 save: bool=True,
                 show: bool=True,
                 ifig: int=1,
                 tag: str='',
                 normalized_inertia_png_filename: PathLike='normalized_inertia_vs_span.png',
                 amoi_span_png_filename: PathLike='amoi_vs_span.png',
                 e_amoi_span_png_filename: PathLike='e_amoi_vs_span.png',
                 centroid_span_png_filename: PathLike='centroid_vs_span.png') -> int:
    """
    Generate matplotlib span-wise plots of section properties.

    Produces four figures: normalised inertia, raw area-MOI, E*I, and
    centroid position, all versus station.

    Parameters
    ----------
    log : SimpleLogger
        logging object
    station : (nstation,) float ndarray
        span positions
    A : (nstation,) float ndarray
        section area
    I : (nstation, 6) float ndarray
        second moments ``[Ixx, Iyy, Izz, Ixy, Iyz, Ixz]``
    J : (nstation,) float ndarray
        polar moment
    ExI : (nstation, 6) float ndarray
        Ex-weighted second moments
    EyI : (nstation, 6) float ndarray
        Ey-weighted second moments
    GJ : (nstation,) float ndarray
        torsion stiffness
    avg_centroid : (nstation, 3) float ndarray
        section centroid at each station
    linestyle : str; default='-'
        matplotlib line-style string
    dirname : PathLike; default=''
        directory for saving PNG files
    x, y, z : str; default='x', 'y', 'z'
        axis labels used in legends and axis titles
    station_word : str; default='Span'
        label for the x-axis (e.g. ``'Span'`` or ``'Buttline'``)
    save : bool; default=True
        write PNG files to *dirname*
    show : bool; default=True
        call ``plt.show()`` after the last figure
    ifig : int; default=1
        starting figure number
    tag : str; default=''
        prefix for legend labels (useful when overlaying multiple cuts)
    normalized_inertia_png_filename : PathLike
        filename for the normalized-inertia figure
    amoi_span_png_filename : PathLike
        filename for the area-MOI figure
    e_amoi_span_png_filename : PathLike
        filename for the E*I figure
    centroid_span_png_filename : PathLike
        filename for the centroid figure

    Returns
    -------
    ifig : int
        the next unused figure number
    """
    absI = np.abs(I)
    absExI = np.abs(ExI)
    # absGJ = np.abs(GJ)

    assert isinstance(ifig, int), ifig
    fig = plt.figure(ifig)
    ax = fig.gca()
    ai_max = absI[:, :3].max(axis=0)
    aei_max = absExI[:, :3].max(axis=0)
    ai_max[ai_max == 0] = 1.
    aei_max[aei_max == 0] = 1.
    assert len(ai_max) == 3, (ai_max.shape, absI)
    xx = f'{x}{x}'
    zz = f'{z}{z}'
    xz = f'{x}{z}'
    log.info(f'ai_max={ai_max}')
    log.info(f'aei_max={aei_max}')
    ax.plot(station, I[:, 0] / ai_max[0], 'ro-', label=f'I{xx}')
    ax.plot(station, I[:, 1] / ai_max[1], 'bo-', label=f'I{zz}')
    ax.plot(station, I[:, 2] / ai_max[2], 'go-', label=f'I{xz}')

    ax.plot(station, ExI[:, 0] / aei_max[0], 'ro', label=f'E{x}I{xx}', linestyle='--')
    ax.plot(station, ExI[:, 1] / aei_max[1], 'bo', label=f'E{x}I{zz}', linestyle='--')
    ax.plot(station, ExI[:, 2] / aei_max[2], 'go', label=f'E{x}I{xz}', linestyle='--')
    #ax.plot(station, GJ / aGJ.max(), 'go-', label='GJ', linestyle='--')

    ax.grid(True)
    ax.set_xlabel(f'{station_word}, {y}')
    ax.set_ylabel('Normalized Area MOI, I')
    ax.legend()
    png_filename = os.path.join(dirname, normalized_inertia_png_filename)
    if save:
        log.info(f'saving {png_filename}')
        fig.savefig(png_filename)
    #-------------------------------------------------------

    fig = plt.figure(ifig + 2)
    ax = fig.gca()
    ax.plot(station, I[:, 0], 'ro', linestyle=linestyle, label='Ixx')
    ax.plot(station, I[:, 1], 'bo', linestyle=linestyle, label='Iyy')
    ax.plot(station, I[:, 2], 'go', linestyle=linestyle, label='Izz')
    ax.grid(True)
    ax.set_xlabel(f'{station_word}, {y}')
    ax.set_ylabel('Area MOI, I')
    ax.legend()
    png_filename = os.path.join(dirname, amoi_span_png_filename)
    if save:
        log.info(f'saving {png_filename}')
        fig.savefig(png_filename)
    #-------------------------------------------------------


    fig = plt.figure(ifig + 3)
    ax = fig.gca()
    ax.plot(station, ExI[:, 0], 'ro', linestyle=linestyle, label=f'EI{xx}')
    #ax.plot(station, I[:, 0], 'bo-', label='Ixx')
    ax.grid(True)
    ax.set_xlabel(f'{station_word}, {y}')
    ax.set_ylabel(f'E{xx}*Area MOI, E{xx}*I')
    ax.legend()
    png_filename = os.path.join(dirname, e_amoi_span_png_filename)
    if save:
        log.info(f'saving {png_filename}')
        fig.savefig(png_filename)
    #-------------------------------------------------------

    fig = plt.figure(ifig + 4)
    ax = fig.gca()
    ax.plot(station, avg_centroid[:, 0], 'ro', linestyle=linestyle, label=f'{x}cg')
    ax.plot(station, avg_centroid[:, 2], 'bo', linestyle=linestyle, label=f'{z}cg')
    ax.grid(True)
    ax.set_xlabel(f'{station_word}, {y}')
    ax.set_ylabel('Centroid')
    ax.legend()
    png_filename = os.path.join(dirname, centroid_span_png_filename)
    if save:
        log.info(f'saving {png_filename}')
        fig.savefig(png_filename)
    #-------------------------------------------------------

    if show:
        plt.show()
    ifig += 4
    return ifig


def calculate_area_moi(model: BDF,
                       rods: Rods,
                       normal_plane: np.ndarray,
                       thetas: dict[int, tuple[float, float, float, float]],
                       moi_filename: PathLike='',
                       eid_filename: PathLike='eid_file.csv',
                       use_bredt_batho: bool=True,
                       use_shear_center: bool=True,
                       bar_data=None,
                       shell_prop_cache: dict[int, tuple[float, float, float, float, float]] | None = None,
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray,               # dxi, dyi, total_area,
                                  np.ndarray, np.ndarray,                           # Isum, Jsum,
                                  np.ndarray, np.ndarray, np.ndarray, np.ndarray]:  # ExIsum, EyIsum, GJsum, avg_centroid
    """
    Integrate section properties at a single cut.

    Each shell wall segment contributes ``A = t * L`` of area at its
    centroid.  The second moments are the parallel-axis (A * d^2) terms;
    the elements' own bending inertia is negligible for thin shells.
    When bar/beam data is supplied, each bar adds a concentrated area at
    its crossing point plus its own bending inertia rotated into the
    cut-plane frame.

    Parameters
    ----------
    model : BDF
        cross-referenced model (used for element look-ups and logging)
    rods : tuple[ndarray, ndarray, ndarray]
        ``(rod_eid_nodes, rod_nids, rod_xyzs)`` from the face cutter,
        describing the wall segments that make up the cut ring
    normal_plane : (3,) float ndarray
        unit normal of the cutting plane in the basic frame
    thetas : dict[int, tuple[float, float, float, float]]
        mutable mapping ``{eid: (theta_deg, Ex, Ey, Gxy)}``; updated
        in-place with the material-angle data for every shell element
        that participates in this cut
    moi_filename : PathLike; default=''
        when non-empty, write a diagnostic BDF/CSV of the cut geometry
    eid_filename : PathLike; default='eid_file.csv'
        companion CSV for *moi_filename*
    use_bredt_batho : bool; default=True
        compute the torsion constant from the closed-cell topology
        (Bredt-Batho); when False, fall back to the polar-moment
        approximation ``GJ = G * (Ix + Iz)``
    use_shear_center : bool; default=True
        solve the transverse shear flow for the shear center; when
        False, the area centroid is reported for all three reference
        points (centroid, neutral axis, shear center)
    bar_data : tuple | None; default=None
        output of ``_find_bar_beam_crossings``; when not None, the five
        arrays ``(centroids, areas, own_I, own_J, E_arr)`` are merged
        into the shell data before integration.  Bar elements do NOT
        enter the Bredt-Batho or shear-center solves.
    shell_prop_cache : dict | None; default=None
        pre-computed ``{eid: (thickness, theta_deg, Ex, Ey, Gxy)}``
        from ``_precompute_shell_props``.  Eliminates the expensive
        per-element ``material_coordinate_system`` and
        ``get_Ainv_equivalent_pshell`` calls inside the station loop.

    Returns
    -------
    dxi : float
        section width (max x - min x after rotation)
    dyi : float
        section height
    total_length : float
        sum of wall-segment arc lengths (perimeter)
    total_area : float
        sum of all wall-segment and bar areas
    Isum : (6,) float ndarray
        ``[Ixx, Iyy, Izz, Ixy, Iyz, Ixz]`` about the area centroid
    Jsum : float
        polar moment ``Ixx + Izz``
    ExIsum : (6,) float ndarray
        modulus-weighted ``Ex * I`` about the area centroid
    EyIsum : (6,) float ndarray
        modulus-weighted ``Ey * I``
    GJsum : float
        torsion stiffness (Bredt-Batho for shells + G*J for bars)
    avg_centroid : (3,) float ndarray
        area-weighted centroid in the cut coord's local frame
    ExAsum : float
        section axial stiffness ``sum(Ex * dA)``
    EyAsum : float
        section transverse stiffness ``sum(Ey * dA)``
    GAsum : float
        section shear stiffness ``sum(G * dA)``
    xyz_neutral_axis : (3,) float ndarray
        modulus-weighted centroid (neutral axis) in the local frame;
        identical to *avg_centroid* for a homogeneous section
    xyz_shear_center : (3,) float ndarray
        where a transverse shear produces no twist; NaN when it could
        not be determined.  This is where the CBEAM element axis belongs.
    """
    assert isinstance(rods, tuple), type(rods)
    assert isinstance(thetas, dict), type(thetas)
    rod_eid_nodes, rod_nids, rod_xyzs = rods
    assert isinstance(rod_eid_nodes, np.ndarray), type(rod_eid_nodes)
    assert isinstance(rod_nids, np.ndarray), type(rod_nids)
    assert isinstance(rod_xyzs, np.ndarray), type(rod_xyzs)

    if 0:
        print(f'rod_eid_nodes:\n{rod_eid_nodes}')
        print(f'rod_nids:\n{rod_nids}')
        print(f'rod_xyzs:\n{rod_xyzs}')

    eids = np.abs(rod_eid_nodes[:, 0])
    neids = len(eids)
    all_nids = rod_nids
    n1 = rod_eid_nodes[:, 1]
    n2 = rod_eid_nodes[:, 2]
    inid1 = np.searchsorted(all_nids, n1)
    inid2 = np.searchsorted(all_nids, n2)
    xyz1 = rod_xyzs[inid1, :]
    xyz2 = rod_xyzs[inid2, :]
    centroid = (xyz1 + xyz2) / 2.
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    assert len(length) == neids

    centroid, length, area, thickness, E = get_element_inertias(
        model, normal_plane, thetas,
        eids, length, centroid,
        shell_prop_cache=shell_prop_cache)

    # ---- bar / beam contributions (lumped areas) ----
    # Bars add concentrated area at their crossing point but do NOT
    # participate in the thin-walled Bredt-Batho or shear-center solves,
    # so the shell-only wall arrays (xyz1, xyz2, thickness, length) are
    # kept untouched.
    nshell = len(area)
    bar_own_I_arr = None
    bar_own_J_arr = None
    if bar_data is not None:
        bar_centroids, bar_areas, bar_own_I, bar_own_J, bar_E = bar_data
        if len(bar_areas) > 0:
            centroid = np.vstack([centroid, bar_centroids])
            area = np.concatenate([area, bar_areas])
            E = np.vstack([E, bar_E])
            bar_own_I_arr = bar_own_I
            bar_own_J_arr = bar_own_J
            model.log.debug(f'  {len(bar_areas):d} bar/beam element(s) cross '
                           f'this cut (total bar area = {bar_areas.sum():g})')

    # [Ixx, Iyy, Izz, Ixy, Iyz, Ixz]
    inertia: np.ndarray = np.zeros((len(area), 6), dtype='float64')

    # (Ex, Ey, Gxy)
    ex = E[:, 0]
    ey = E[:, 1]
    gxy = E[:, 2]

    total_length = length.sum()
    total_area = area.sum()
    if total_area == 0.0:
        avg_centroid = centroid.mean(axis=0)
    else:
        avg_centroid = (centroid * area[:, np.newaxis]) .sum(axis=0) / total_area
    assert len(avg_centroid) == 3, len(avg_centroid)
    # y corresponds to the station in the plane of the coordinate system
    # and is 0. because we're in the local plane
    x = centroid[:, 0] - avg_centroid[0]
    y = centroid[:, 1] - avg_centroid[1]
    z = centroid[:, 2] - avg_centroid[2]

    xmin = x.min()
    xmax = x.max()
    ixmins = np.where(x == xmin)[0]
    ixmaxs = np.where(x == xmax)[0]
    try:
        ixmin = ixmins[0]
        ixmax = ixmaxs[0]
    except IndexError:
        print(f'x = {x}')
        print(f'ixmins = {ixmins}')
        print(f'ixmaxs = {ixmaxs}')
        raise RuntimeError('bad cut?')

    xyz_min = centroid[ixmin, :]
    xyz_max = centroid[ixmax, :]
    d = xyz_max - xyz_min
    dx = d[0]
    dz = d[2]
    theta = np.arctan2(dx, dz)
    thetad = np.degrees(theta)

    nnodes = len(x)
    delta = np.zeros((nnodes, 3))
    delta[:, 1] = thetad

    xyz = np.zeros((nnodes, 3))
    # we're swapping what axes we have to make the transform easier
    xyz[:, 0] = x
    xyz[:, 1] = z
    xyz[:, 2] = 0.
    rtz = xyz_to_rtz_array(xyz)
    rtz2 = rtz + delta
    # NOTE: this used to be called xyz2, which shadowed the wall end
    # coordinate of the same name computed above; _write_moi_file was being
    # handed this rotation scratch array instead of the real wall end
    xyz_rot = rtz_to_xyz_array(rtz2)
    x2 = xyz_rot[:, 0]
    y2 = xyz_rot[:, 1]
    #z2 = xyz_rot[:, 2]

    #origin = d
    #zaxis = np.array([0., 1., 0.])
    #xzplane = d
    dxi = x2.max() - x2.min()
    dyi = y2.max() - y2.min()
    #dzi = z2.max() - z2.min()  # zero by definition

    inertia[:, 0] = area * (x * x)  # Izz
    inertia[:, 1] = area * (y * y)  # just 0
    inertia[:, 2] = area * (z * z)  # Ixx
    inertia[:, 3] = area * (x * y)  # just 0
    inertia[:, 4] = area * (y * z)  # just 0
    inertia[:, 5] = area * (x * z)  # Ixz

    # Add bar/beam elements' own bending inertia (about their centroids,
    # already rotated into the cut-plane frame).  The A*d^2 parallel-axis
    # terms are already included above; these are the self-inertia remainder.
    if bar_own_I_arr is not None and len(bar_own_I_arr) > 0:
        inertia[nshell:, 0] += bar_own_I_arr[:, 0]  # Ixx_own
        inertia[nshell:, 2] += bar_own_I_arr[:, 1]  # Izz_own
        inertia[nshell:, 5] += bar_own_I_arr[:, 2]  # Ixz_own

    # cut is in xz plane
    ix = inertia[:, 0]
    iz = inertia[:, 2]
    J = ix + iz

    Isum = inertia.sum(axis=0)
    Jsum = J.sum()
    ExIsum = (ex[:, np.newaxis] * inertia).sum(axis=0)
    EyIsum = (ey[:, np.newaxis] * inertia).sum(axis=0)
    # Shell GJ uses the polar-moment approximation as a default; Bredt-Batho
    # may replace it below.  Bar/beam GJ is always the element's own torsion
    # constant G*J_bar (added after Bredt-Batho), NOT its polar second moment.
    GJsum = (gxy[:nshell] * J[:nshell]).sum()
    assert len(Isum) == 6, len(Isum)

    # Modulus-weighted areas.  ExA = sum(Ex_i*dA_i) is the section axial
    # stiffness EA; GA = sum(Gxy_i*dA_i) is the transverse shear stiffness.
    # Neither is recoverable from ExI/I on a heterogeneous section, because
    # ExI is weighted by x^2 while EA is weighted by area, so they have to be
    # accumulated here.
    ExAsum = (ex * area).sum()
    EyAsum = (ey * area).sum()
    GAsum = (gxy * area).sum()

    # Torsion.  G*(Ix+Iz) is the polar moment, which is the torsion constant
    # only for a closed circular section; a wing box is far stiffer in the
    # polar measure than it really is in torsion.  Recover the cell topology
    # from the wall connectivity and solve Bredt-Batho instead.
    if use_bredt_batho and len(xyz1) == len(thickness):
        gj_bredt, torsion_method, ncells = bredt_batho_gj(
            xyz1, xyz2, length, thickness, gxy[:nshell], log=model.log)
        if torsion_method == 'none':
            model.log.warning(
                'no usable walls for the torsion calculation; '
                'falling back to GJ = G*(Ix+Iz)')
        else:
            if torsion_method == 'open':
                # orders of magnitude softer than a closed cell, so this must
                # not pass silently; the usual cause is a station landing on a
                # rib/bulkhead plane, where the in-plane element filter drops
                # the coincident shells and breaks the loop
                model.log.warning(
                    'no closed cell found in this cut; using the open-section '
                    f'GJ={gj_bredt:g}, which is far softer than a closed '
                    'section. If the section really is closed, move the '
                    'station off the rib/bulkhead plane.')
            else:
                model.log.debug(
                    f'GJ: {torsion_method} section, {ncells:d} cell(s); '
                    f'Bredt-Batho GJ={gj_bredt:g} vs polar G*Ip={GJsum:g}')
            GJsum = gj_bredt

    # Shear center.  MSC puts the CBEAM element axis on it (Figure 16-142):
    # GA+WA rides the shear-center line and N1/N2 then point at the neutral
    # axis.  For a doubly symmetric section all three points coincide and none
    # of this matters; for an airfoil box they are inches apart, and putting
    # the axis on the centroid instead silently couples bending into torsion.
    xyz_shear_center = np.full(3, np.nan, dtype='float64')
    xyz_neutral_axis = avg_centroid.copy()
    if use_shear_center and len(xyz1) == len(thickness):
        result = shear_center(xyz1, xyz2, length, thickness,
                              ex[:nshell], gxy[:nshell], log=model.log)
        # the section lives in the cut plane, so the out-of-plane coordinate is
        # whatever the centroid has; carrying it keeps the point on the plane
        # when the caller transforms back to the basic frame
        xyz_neutral_axis[0] = result.xy_neutral_axis[0]
        xyz_neutral_axis[2] = result.xy_neutral_axis[1]
        if result.method == 'none':
            model.log.warning(
                'no shear center for this cut; the beam axis will fall back '
                'to the neutral axis, which is only right for a doubly '
                'symmetric section')
        else:
            xyz_shear_center[:] = xyz_neutral_axis
            xyz_shear_center[0] = result.xy_shear_center[0]
            xyz_shear_center[2] = result.xy_shear_center[1]
            model.log.debug(
                f'shear center: {result.method} section, {result.ncells:d} '
                f'cell(s), offset from the neutral axis '
                f'({result.xy_shear_center - result.xy_neutral_axis}), '
                f'force residual {result.force_error:g}')

    # Add bar/beam torsion stiffness (G * J_bar) on top of the shell GJ.
    # This is the element's true torsion constant, not a polar-moment
    # approximation -- it is correct even for open cross-sections.
    if bar_own_J_arr is not None and len(bar_own_J_arr) > 0:
        GJsum += (gxy[nshell:] * bar_own_J_arr).sum()

    if moi_filename is not None:
        dirname = os.path.dirname(moi_filename)
        eid_filename = os.path.join(dirname, eid_filename)
        _write_moi_file(
            moi_filename, eid_filename,
            eids, n1, n2, xyz1, xyz2, length, thickness, area,
            centroid, avg_centroid, inertia, E,
        )
    out = (
        dxi, dyi, total_length, total_area,
        Isum, Jsum,
        ExIsum, EyIsum, GJsum, avg_centroid,
        ExAsum, EyAsum, GAsum,
        xyz_neutral_axis, xyz_shear_center,
    )
    return out


def _write_moi_file(moi_filename: PathLike,
                    eid_filename: PathLike,
                    eids, n1, n2, xyz1, xyz2,
                    length, thickness, area,
                    centroid, avg_centroid, I, E) -> None:
    """
    Write a diagnostic BDF of the cut geometry as CONROD elements and a
    companion CSV of per-element section data.

    The BDF can be loaded in a viewer to visually verify that the cut
    ring looks correct; the CONROD areas carry the element's cut-area
    contribution so the ring's total area matches the integrated value.
    """
    eidi = 1
    mid = 1
    nid0 = max(n1.max(), n2.max()) + 1
    with open(moi_filename, 'w') as bdf_file, open(eid_filename, 'w') as eid_file:
        bdf_file.write('$ pyNastran: punch=True\n')
        bdf_file.write('MAT1,1,3.0e7,,0.3\n')
        grid = ['GRID', nid0, 0, avg_centroid[0], avg_centroid[2], 0.]
        bdf_file.write(print_card_8(grid))
        bdf_file.write(f'CONM2   {1:8d}{nid0:8d}\n')

        fmt = ('%s,' * 7)[:-1] + '\n'
        eid_file.write('# eid(%i),pid(%i),area,thickness,Ixx,Izz,Ixz\n')
        for eid, n1i, n2i, xyz1i, xyz2i, lengthi, thicknessi, areai, centroidi, Ii, Ei in zip(
                eids, n1, n2, xyz1, xyz2, length, thickness, area, centroid, I, E):
            actual_eid = abs(eid)

            assert nid0 not in [n1i, n2i], (n1i, n2i)
            pidi = actual_eid
            #pid = eidi
            grid1 = ['GRID', n1i, None] + xyz1i.tolist()
            grid2 = ['GRID', n2i, None] + xyz2i.tolist()
            #crod = ['CROD', eidi, pid, n1i, n2i]
            A, J, nsm = Ii
            #prod = ['PROD', pid, mid, A, J, 0., nsm]
            assert eidi > 0, eidi
            conrod = ['CONROD', eidi, n1i, n2i, mid, A, J, 0., nsm]
            bdf_file.write(print_card_8(grid1))
            bdf_file.write(print_card_8(grid2))
            #bdf_file.write(print_card_8(crod))
            #bdf_file.write(print_card_8(prod))
            bdf_file.write(print_card_8(conrod))
            eidi += 1
            #PID | MID |  A  |  J  |  C  | NSM
            eid_file.write(fmt % (eidi, pidi, areai, thicknessi, Ii[0], Ii[1], Ii[2]))


def get_element_inertias(model: BDF,
                         normal_plane: np.ndarray,
                         thetas: dict[int, tuple[float, float, float, float]],
                         eids: list[int],
                         length: list[float],
                         centroid: list[np.ndarray],
                         shell_prop_cache: dict[int, tuple[float, float, float, float, float]] | None = None,
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract per-element thickness, area, and equivalent moduli for every
    shell element in the cut ring.

    Iterates over the element ids produced by the face cutter, resolves
    each element's material orientation angle relative to the cut plane,
    and calls ``get_Ainv_equivalent_pshell`` to obtain ``(Ex, Ey, Gxy)``
    in the cut frame.  Elements whose normals are nearly in-plane with
    the cut (> ~25 deg off the cut-plane normal) are zeroed out so they
    do not contribute stiffness (they are skins seen edge-on, not
    structural walls).

    Parameters
    ----------
    model : BDF
        cross-referenced model
    normal_plane : (3,) float ndarray
        unit normal of the cutting plane
    thetas : dict
        mutable ``{eid: (theta, Ex, Ey, Gxy)}``; updated in-place
    eids : (n,) int ndarray
        element ids from the face cutter
    length : (n,) float ndarray
        segment lengths of each wall piece
    centroid : (n, 3) float ndarray
        segment centroids
    shell_prop_cache : dict or None
        pre-computed ``{eid: (thickness, theta_deg, Ex, Ey, Gxy)}``
        from ``_precompute_shell_props``.  When provided the expensive
        per-element ``material_coordinate_system`` and
        ``get_Ainv_equivalent_pshell`` calls are skipped.

    Returns
    -------
    centroid : (n, 3) float ndarray
        (may be filtered from the input)
    length : (n,) float ndarray
    area : (n,) float ndarray
        ``thickness * length`` for each segment
    thickness : (n,) float ndarray
    E : (n, 3) float ndarray
        ``(Ex, Ey, Gxy)`` per element
    """
    cg_list: list[np.ndarray] = []
    area_list: list[float] = []
    length_list: list[float] = []
    thickness_list: list[float] = []
    E_list: list[tuple[float, float, float]] = []

    log = model.log
    normal_plane_vector = normal_plane.copy().reshape((3, 1))
    for eid, lengthi, centroidi in zip(eids, length, centroid):
        element = model.elements[eid]
        if element.type not in _SHELL_TYPES:
            log.warning(element)
            continue

        # look up the cache first; on a miss compute once and store
        if shell_prop_cache is not None and eid in shell_prop_cache:
            thicknessi, thetad, Ex, Ey, Gxy = shell_prop_cache[eid]
        else:
            thicknessi, _areai, thetad, Ex, Ey, Gxy, nu_xy = _get_shell_inertia(
                element, normal_plane, normal_plane_vector, lengthi)
            if shell_prop_cache is not None:
                shell_prop_cache[eid] = (thicknessi, thetad, Ex, Ey, Gxy)

        areai = thicknessi * lengthi
        thetas[eid] = (thetad, Ex, Ey, Gxy)
        thickness_list.append(thicknessi)
        length_list.append(lengthi)
        area_list.append(areai)
        cg_list.append(centroidi)
        E_list.append((Ex, Ey, Gxy))

    centroid = np.array(cg_list, dtype='float64')
    length2 = np.array(length_list, dtype='float64')
    assert np.allclose(length, length2)
    area = np.array(area_list, dtype='float64')
    thickness = np.array(thickness_list, dtype='float64')
    E = np.array(E_list, dtype='float64')
    return centroid, length, area, thickness, E


# element types handled by the section-property extraction
_SHELL_TYPES = frozenset(['CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8'])


def _precompute_shell_props(
        model: BDF,
        normal_plane: np.ndarray,
        candidate_eids: set[int] | None = None,
        ) -> dict[int, tuple[float, float, float, float, float]]:
    """
    Pre-compute the thickness and equivalent moduli for shell elements.

    The result is a cache ``{eid: (thickness, theta_deg, Ex, Ey, Gxy)}``
    that is valid for all stations sharing the same *normal_plane*.
    Elements whose normals are nearly in-plane with the cut
    (``|cos θ| > 0.9``) are stored with zeroed properties so they
    contribute no stiffness, matching the behaviour of
    ``_get_shell_inertia``.

    When *candidate_eids* is supplied, only those elements are
    processed (typically the set of element ids that appear in the
    triangulated face mesh).  This avoids computing expensive material
    properties for elements that can never be intersected by any
    station.

    This is the single biggest performance win in the module: the
    ``material_coordinate_system`` and ``get_Ainv_equivalent_pshell``
    calls are expensive (~0.3 ms each), and the per-element material
    properties do not depend on *which* station the cut is at — only
    the segment length changes.
    """
    normal_plane_vector = normal_plane.copy().reshape((3, 1))
    cache: dict[int, tuple[float, float, float, float, float]] = {}

    for eid, element in model.elements.items():
        if element.type not in _SHELL_TYPES:
            continue
        if candidate_eids is not None and eid not in candidate_eids:
            continue

        pid_ref = element.pid_ref
        thicknessi = element.Thickness()

        dxyz, centroid_unused, imat, unused_jmat, element_normal = \
            element.material_coordinate_system()

        n1, n2, n3 = element_normal
        R1 = np.array([
            [0., -n3, n2],
            [n3, 0., -n1],
            [-n2, n1, 0.],
        ], dtype='float64')
        R2 = np.array([
            [1 - n1 ** 2, -n1 * n2, -n1 * n3],
            [-n1 * n2, 1 - n2 ** 2, -n2 * n3],
            [-n1 * n3, -n2 * n3, 1 - n3 ** 2],
        ])
        imat_col = imat.reshape(3, 1)
        b = np.linalg.multi_dot([normal_plane_vector.T, R1, imat_col])
        c = np.linalg.multi_dot([normal_plane_vector.T, R2, imat_col])
        imat_rotation_angle = np.arctan2(b, c).item()
        imat_rotation_angle_deg = np.degrees(imat_rotation_angle)
        if imat_rotation_angle_deg <= -90.:
            imat_rotation_angle_deg += 180.
        elif imat_rotation_angle_deg > 90.:
            imat_rotation_angle_deg -= 180.

        abs_cos_theta = abs(normal_plane @ element_normal)
        if abs_cos_theta > 0.9:  # <25.8 degrees off the cut → in-plane element
            cache[eid] = (0., 0., 0., 0., 0.)
        else:
            Ex, Ey, Gxy, nu_xy = pid_ref.get_Ainv_equivalent_pshell(
                imat_rotation_angle_deg, thicknessi)
            cache[eid] = (thicknessi, imat_rotation_angle_deg, Ex, Ey, Gxy)
    return cache


def _get_shell_inertia(element: CTRIA3 | CQUAD4,
                       normal_plane: np.ndarray,
                       normal_plane_vector: np.ndarray,
                       lengthi: float,) -> tuple[float, float, float,
                                                 float, float, float, float]:
    """
    Parameters
    ----------
    element : CTRIA3 / CQUAD4
        the object to cut
    normal_plane : (3,) float ndarray
        the normal vector of the cutting plane (should be roughly normal to the element face)
    normal_plane_vector : (3,1) float ndarray
        the normal vector of the cutting plane (should be roughly normal to the element face)
    lengthi : float
        the length the cutting plane makes with the element

    Returns
    -------
    thicknessi : float
        the total thickness of the element
    areai : float
        the cut area of the element
    imat_rotation_angle_deg : float
        the angle between the cutting plane and the normal_plane / normal_plane_vector
        this is NOT the angle of the fiber
    Ex : float
        the moduli normal to the cut plane
    Ey : float
        the moduli parallel to the cut plane (normal to Ex)
    Gxy : float
        the inplane shear moduli
    nu_xy : float
        the correlary to in-plane nu12

    """
    pid_ref = element.pid_ref
    thicknessi = element.Thickness()
    dxyz, centroid, imat, unused_jmat, element_normal = element.material_coordinate_system()
    #print('imat = ', imat)
    #print('normal = ', normal)
    n1, n2, n3 = element_normal
    n12 = n1 * n2
    n13 = n1 * n3
    n23 = n2 * n3
    # http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
    # expanding eq 20 into
    # R(n,theta) = R0 + R1*sin(theta) + R2*cos(theta)
    # R0 = np.array([
    #     [n1 ** 2, n12, n13],
    #     [n12, n2 ** 2, n23],
    #     [n13, n23, n3 ** 2],
    # ], dtype='float64')
    R1 = np.array([
        [0., -n3, n2],
        [n3, 0., -n1],
        [-n2, n1, 0.],
    ], dtype='float64')
    R2 = np.array([
        [1 - n1 ** 2, -n12, -n13],
        [-n12, 1 - n2 ** 2, -n23],
        [-n13, -n23, 1 - n3 ** 2],
    ])
    imat = imat.reshape(3, 1)
    #print(normal_plane.shape, R1.shape, imat.shape)
    #a = np.linalg.multi_dot([normal_plane.T, R0, imat])
    b = np.linalg.multi_dot([normal_plane_vector.T, R1, imat])
    c = np.linalg.multi_dot([normal_plane_vector.T, R2, imat])

    #  maximize m' dot p = p.T dot m
    # m' = R dot m
    #    = a + b*sin(theta) + c*cos(theta)
    #  d/d(theta) = b*cos(theta)*sin(theta) = 0
    #
    #  d/d(theta) = b*cos(theta) - c*sin(theta) = 0
    #  b*cos(theta) = c*sin(theta)
    #  tan(theta) = b/c
    #
    # the theta to rotate by in order to orient imat with the normal
    #print(b, c)
    imat_rotation_angle = np.arctan2(b, c).item()
    imat_rotation_angle_deg = np.degrees(imat_rotation_angle)
    if imat_rotation_angle_deg <= -90.:
        imat_rotation_angle_deg += 180.
    elif imat_rotation_angle_deg > 90.:
        imat_rotation_angle_deg -= 180.

    #element_normal = element.Normal()
    # cos(theta) = a o b / (|a| * |b|)
    # |a| = length of normal vector = 1.0
    # |b| = length of normal_plane vector = 1.0
    #
    # cos(theta) = a o b
    # then we take the absolute value because we don't care if the element is +/- theta off

    abs_cos_theta = abs(normal_plane @ element_normal)
    assert isinstance(imat_rotation_angle, float), imat_rotation_angle
    if abs_cos_theta > 0.9:  # <25.8 degrees
        # filter out elements that are in-plane
        thicknessi = 0.
        areai = 0.
        Ex = 0.
        Ey = 0.
        Gxy = 0.
        nu_xy = 0.
    else:
        Ex, Ey, Gxy, nu_xy = pid_ref.get_Ainv_equivalent_pshell(
            imat_rotation_angle_deg, thicknessi, # degrees=True,
        )

        #thicknessi = prop.Thickness()
        areai = thicknessi * lengthi

    # pid = pid_ref.pid
    # if pid == 10:
    #     import copy
    #     pid_ref45 = copy.deepcopy(pid_ref)
    #     pid_ref45.mids_ref = [copy.deepcopy(pid_ref.mids_ref[0])]
    #     pid_ref45.thetas = [copy.deepcopy(pid_ref.thetas[0])]
    #     pid_ref45.thicknesses = [copy.deepcopy(pid_ref.thicknesses[0])]
    #     pid_ref45.mids = [copy.deepcopy(pid_ref.mids[0])]
    #     pid_ref45.get_thetas()
    #     Ex45, Ey45, Gxy45, nu_xy45 = pid_ref45.get_Ainv_equivalent_pshell(
    #         imat_rotation_angle_deg, thicknessi)
    #
    #     pid_ref0 = copy.deepcopy(pid_ref)
    #     pid_ref0.mids_ref = [copy.deepcopy(pid_ref.mids_ref[1])]
    #     pid_ref0.thetas = [copy.deepcopy(pid_ref.thetas[1])]
    #     pid_ref0.thicknesses = [copy.deepcopy(pid_ref.thicknesses[1])]
    #     pid_ref0.mids = [copy.deepcopy(pid_ref.mids[1])]
    #     pid_ref0.get_thetas()
    #     Ex0, Ey0, Gxy0, nu_xy0 = pid_ref0.get_Ainv_equivalent_pshell(
    #         imat_rotation_angle_deg, thicknessi)
    return thicknessi, areai, imat_rotation_angle_deg, Ex, Ey, Gxy, nu_xy


def _get_bar_section_props(prop, t: float=0.5,
                           ) -> tuple[float, float, float, float, float]:
    """
    Extract cross-section constants from a bar/beam property card.

    For PBEAM (tapered), properties are linearly interpolated at fraction
    *t* along the element's ``xxb`` stations (0 = end A, 1 = end B).
    For PBAR, PBARL and PBEAML the section is constant and *t* is
    ignored.

    Parameters
    ----------
    prop : PBAR | PBEAM | PBARL | PBEAML
        cross-referenced property object
    t : float; default=0.5
        fractional position along the element (used only for PBEAM)

    Returns
    -------
    A : float
        cross-section area
    i1 : float
        second moment of area about the element y-axis (I1)
    i2 : float
        second moment of area about the element z-axis (I2)
    i12 : float
        product of inertia (I12); zero for symmetric standard shapes
    j : float
        torsion constant (Saint-Venant J)
    """
    # PBEAM: tapered, properties vary along xxb
    if hasattr(prop, 'xxb') and prop.xxb is not None and len(prop.xxb) > 1:
        xxb = np.asarray(prop.xxb, dtype='float64')
        A = float(np.interp(t, xxb, prop.A))
        i1 = float(np.interp(t, xxb, prop.i1))
        i2 = float(np.interp(t, xxb, prop.i2))
        i12 = float(np.interp(t, xxb, prop.i12))
        j = float(np.interp(t, xxb, prop.j))
        return A, i1, i2, i12, j

    # PBAR, PBARL, PBEAML: constant section
    A = prop.Area()
    i1 = prop.I11()
    i2 = prop.I22()
    j = prop.J()
    # I12: available as attribute on PBAR; 0 for standard shapes (PBARL/PBEAML)
    i12 = getattr(prop, 'i12', 0.0)
    if i12 is None:
        i12 = 0.0
    return A, i1, i2, i12, j


def _bar_own_I_in_cut_frame(i1: float, i2: float, i12: float,
                             ga_global: np.ndarray,
                             gb_global: np.ndarray,
                             elem,
                             plane_i: np.ndarray,
                             plane_k: np.ndarray,
                             model: BDF,
                             ) -> tuple[float, float, float]:
    """
    Rotate a bar/beam element's own bending inertias from its element
    y_e / z_e axes into the cut-plane ``(plane_i, plane_k)`` system.

    The element triad is built from the node positions and the CBAR/CBEAM
    orientation vector (``elem.x`` or ``elem.g0``), following the same
    convention as MSC Nastran with ``offt='GGG'``.  The element's
    ``y_e`` and ``z_e`` directions are then projected onto the cut-plane
    axes and the 2x2 inertia tensor is rotated::

        M_cut = R^T  M_elem  R

    where ``R = [[cy, sy], [cz, sz]]`` are the in-plane direction
    cosines of ``y_e`` and ``z_e`` relative to ``plane_i`` and
    ``plane_k``.

    Parameters
    ----------
    i1, i2, i12 : float
        second moments and product of inertia about the element's own
        y_e and z_e axes (from the PBAR/PBEAM card)
    ga_global, gb_global : (3,) float ndarray
        end-point positions of the bar in the basic frame
    elem : CBAR | CBEAM
        the element object (used to read the orientation vector)
    plane_i, plane_k : (3,) float ndarray
        unit vectors of the cut coord's local x and z axes
    model : BDF | None
        used only when ``elem.g0`` is set (to look up the G0 grid);
        may be None when the orientation is given by ``elem.x``

    Returns
    -------
    ixx : float
        ``int(x_loc^2 dA)`` — second moment about the cut-plane z-axis
    izz : float
        ``int(z_loc^2 dA)`` — second moment about the cut-plane x-axis
    ixz : float
        ``int(x_loc * z_loc dA)`` — product of inertia in the cut plane
    """
    dxyz = gb_global - ga_global
    L = np.linalg.norm(dxyz)
    if L < 1e-12:
        return 0., 0., 0.
    x_e = dxyz / L

    # orientation vector  (global frame, same convention as CBEAM offt='GGG')
    if elem.g0 is not None:
        g0_xyz = elem.g0_ref.get_position()
        v = g0_xyz - ga_global
    else:
        v = np.asarray(elem.x, dtype='float64')

    y_e = v - v.dot(x_e) * x_e
    norm_y = np.linalg.norm(y_e)
    if norm_y < 1e-8:
        return 0., 0., 0.
    y_e /= norm_y
    z_e = np.cross(x_e, y_e)

    # project element axes onto the cut plane
    def _components(vec):
        c = vec.dot(plane_i)
        s = vec.dot(plane_k)
        n = np.hypot(c, s)
        if n < 1e-8:
            return 0., 0.
        return c / n, s / n

    cy, sy = _components(y_e)
    cz, sz = _components(z_e)

    # 2-D tensor rotation:  M_cut = R^T  M_elem  R
    # with R = [[cy, sy], [cz, sz]]
    ixx = cy * cy * i1 + 2. * cy * cz * i12 + cz * cz * i2
    izz = sy * sy * i1 + 2. * sy * sz * i12 + sz * sz * i2
    ixz = cy * sy * i1 + (cy * sz + sy * cz) * i12 + cz * sz * i2
    return ixx, izz, ixz


def _find_bar_beam_crossings(model: BDF,
                              coord: CORD2R,
                              log=None,
                              ) -> tuple[np.ndarray, np.ndarray,
                                         np.ndarray, np.ndarray, np.ndarray]:
    """
    Find every CBAR / CBEAM that straddles the cut plane and return
    its section contribution.

    The cut plane is defined as ``y_local = 0`` in *coord*.  A bar
    "straddles" it when its two end nodes have opposite signs of
    y_local (strictly: ``ya * yb <= 0`` with ``|ya - yb| > 0``).
    The crossing point is found by linear interpolation along the
    element, and the section properties are evaluated at that fraction
    (constant for PBAR, interpolated for tapered PBEAM).

    The bar's own bending inertia ``(I1, I2, I12)`` is rotated from
    the element's local axes into the cut-plane frame via
    ``_bar_own_I_in_cut_frame``, so that a rectangular cross-section
    oriented tangentially to the skin contributes the correct amounts
    to the section's Ixx and Izz.

    Parameters
    ----------
    model : BDF
        cross-referenced model
    coord : CORD2R
        the cutting coordinate system; the cut plane is its local
        xz-plane (y_local = 0)
    log : SimpleLogger | None
        logger for debug messages

    Returns
    -------
    centroids : (nbar, 3) float ndarray
        crossing-point coordinates in *coord*'s local frame
    areas : (nbar,) float ndarray
        cross-section area at the cut
    own_I : (nbar, 3) float ndarray
        ``(Ixx, Izz, Ixz)`` of each bar's own bending inertia about
        its centroid, already rotated into the cut-plane frame
    own_J : (nbar,) float ndarray
        torsion constant (Saint-Venant J) at the cut
    E_arr : (nbar, 3) float ndarray
        ``(E, E, G)`` per bar — axial modulus (twice, for Ex/Ey
        compatibility with the shell convention) and shear modulus
    """
    centroids: list[np.ndarray] = []
    areas: list[float] = []
    own_I: list[tuple[float, float, float]] = []
    own_J: list[float] = []
    E_arr: list[tuple[float, float, float]] = []

    plane_i = coord.i
    plane_k = coord.k

    for eid, elem in model.elements.items():
        if elem.type not in ('CBAR', 'CBEAM'):
            continue

        nids = elem.node_ids
        ga = model.nodes[nids[0]].get_position()
        gb = model.nodes[nids[1]].get_position()

        ga_local = coord.transform_node_to_local(ga)
        gb_local = coord.transform_node_to_local(gb)

        ya = ga_local[1]
        yb = gb_local[1]

        # must straddle (or touch) y = 0
        if ya * yb > 0.:
            continue
        dy = yb - ya
        if abs(dy) < 1e-12:
            continue

        t = -ya / dy  # fraction from A to B at y = 0
        xyz_cross = ga_local + t * (gb_local - ga_local)

        # section properties at the crossing fraction
        prop = elem.pid_ref
        mat = prop.mid_ref
        E_val = mat.E()
        G_val = mat.G()

        A_val, i1, i2, i12_val, j_val = _get_bar_section_props(prop, t)
        if A_val <= 0.:
            continue

        ixx, izz, ixz = _bar_own_I_in_cut_frame(
            i1, i2, i12_val, ga, gb, elem, plane_i, plane_k, model)

        centroids.append(xyz_cross)
        areas.append(A_val)
        own_I.append((ixx, izz, ixz))
        own_J.append(j_val)
        E_arr.append((E_val, E_val, G_val))

    if len(centroids) == 0:
        return (np.empty((0, 3), dtype='float64'),
                np.empty(0, dtype='float64'),
                np.empty((0, 3), dtype='float64'),
                np.empty(0, dtype='float64'),
                np.empty((0, 3), dtype='float64'))

    return (np.array(centroids, dtype='float64'),
            np.array(areas, dtype='float64'),
            np.array(own_I, dtype='float64'),
            np.array(own_J, dtype='float64'),
            np.array(E_arr, dtype='float64'))


def plot_compare_inertia(log: SimpleLogger,
                         csv_filenames: list[tuple[Path, str, str]],
                         # ifig: int=1,
                         x: str='x',
                         y: str='y',
                         z: str='z',
                         yrange=None,
                         span_label = 'Span, y (in)',
                         dirname: Path='', save: bool=True,
                         ylim_GJ_ratio=None,
                         ylim_EyIzz_ratio=None,
                         show: bool=True) -> int:
    """
    Overlay section-property span plots from multiple CSV files for
    comparison (e.g. shell-only vs. shell+bars, or two mesh densities).

    Each CSV is loaded via ``load_moi_data`` and plotted on shared axes
    with the caller-supplied line style and tag.

    Parameters
    ----------
    log : SimpleLogger
        logging object
    csv_filenames : list[tuple[Path, str, str]]
        ``[(csv_path, legend_tag, linestyle), ...]`` for each data set
    x, y, z : str
        axis-label characters (default ``'x'``, ``'y'``, ``'z'``)
    yrange : tuple | None
        ``(ymin, ymax)`` to clip the span axis; None for auto
    span_label : str
        x-axis label, e.g. ``'Span, y (in)'``
    dirname : Path
        directory for saving PNG files
    save : bool; default=True
        write PNG files
    ylim_GJ_ratio : tuple | None
        y-axis limits for the GJ-ratio subplot
    ylim_EyIzz_ratio : tuple | None
        y-axis limits for the EyIzz-ratio subplot
    show : bool; default=True
        call ``plt.show()`` after the last figure

    Returns
    -------
    ifig : int
        the next unused figure number
    """
    xx = f'{x}{x}'
    xy = f'{x}{y}'
    xz = f'{x}{z}'
    # yz = f'{y}{z}'
    zz = f'{z}{z}'
    # yy = f'{y}{y}'
    xy = ''.join(sorted(xy))
    xz = ''.join(sorted(xz))

    marker = ''
    # absI = np.abs(I)
    # absExI = np.abs(ExI)
    # absGJ = np.abs(GJ)

    save0 = save
    save = False
    ilast_file = len(csv_filenames) - 1
    data = []
    for ifile, (cut_data_span_filename, tag, linestyle) in enumerate(csv_filenames):
        is_last_file = (ifile == ilast_file)
        if is_last_file:
            # last file
            save = save0

        log.info(f'loading {str(cut_data_span_filename)}')
        datai = load_moi_data(cut_data_span_filename)
        data.append(datai)
        station, A, I, J, ExI, EyI, GJ, avg_centroid = datai
        # assert station.max() < 2000., station.max()

        # the y-terms are all 0.0
        Ixx = I[:, 0]
        # Iyy = I[:, 1]
        Izz = I[:, 2]
        # Ixy = I[:, 3]
        # Iyz = I[:, 4]
        Ixz = I[:, 5]

        # These are average moduli
        Ex = ExI[:, 0] / Ixx
        Ey = EyI[:, 0] / Ixx
        G = GJ / J
        EyA =  Ex * A

        # assert isinstance(ifig, int), ifig
        #-------------------------------------------------------
        ifig = 1
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, A, 'r', marker=marker, label=f'{tag}Area', linestyle=linestyle)
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel('Area, A ($in^2$)')
        ax.legend()

        png_filename = dirname / f'Area{y}_vs_span.png'
        if save:
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1
        #-------------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.semilogy(station, Ixx, 'C0', linestyle=linestyle, marker=marker, label=f'{tag}I{xx}')
        # ax.semilogy(station, Iyy, 'C1', linestyle=linestyle, marker=marker, label='Iyy')
        ax.semilogy(station, Izz, 'C2', linestyle=linestyle, marker=marker, label=f'{tag}I{zz}')
        # ax.semilogy(station, Ixy, 'C3', linestyle=linestyle, marker=marker, label='Ixy')
        # ax.semilogy(station, Iyz, 'C4', linestyle=linestyle, marker=marker, label='Iyz')
        ax.semilogy(station, np.abs(Ixz), 'C5', linestyle=linestyle, marker=marker, label=f'{tag}I{xz}')
        ax.semilogy(station, J, 'k', linestyle=linestyle, marker=marker, label=f'{tag}J')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel('Area MOI, I ($in^4$)')
        ax.legend()
        png_filename = dirname / f'I{xx}_I{zz}_I{xz}_J_log.png'
        if save:
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1

        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, Ixx, 'C0', linestyle=linestyle, marker=marker, label=f'{tag}I{xx}')
        # ax.plot(station, Iyy, 'C1', linestyle=linestyle, marker=marker, label='Iyy')
        ax.plot(station, Izz, 'C2', linestyle=linestyle, marker=marker, label=f'{tag}I{zz}')
        # ax.plot(station, Ixy, 'C3', linestyle=linestyle, marker=marker, label='Ixy')
        # ax.plot(station, Iyz, 'C4', linestyle=linestyle, marker=marker, label='Iyz')
        ax.plot(station, Ixz, 'C5', linestyle=linestyle, marker=marker, label=f'{tag}I{xz}')
        ax.plot(station, J, 'k', linestyle=linestyle, marker=marker, label=f'{tag}J')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel('Area MOI, I ($in^4$)')
        ax.legend()
        png_filename = dirname / f'I{xx}_I{zz}_I{xz}_J.png'
        if save:
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1

        #-------------------------------------------------------
        # ifig_EyIzz = ifig
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, ExI[:, 2], color='r', linestyle=linestyle, marker=marker, label=f'{tag}E{y}*I{zz}')  # Ey*Izz
        #ax.plot(station, I[:, 0], 'bo-', label='Ixx')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness: E{y}*I{zz}')
        ax.legend()
        png_filename = dirname / f'stiffness_E{y}I{zz}.png'
        if save:
            log.debug(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1

        #---------------------------------------------------
        # ifig_GJ = ifig
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, GJ, color='k', linestyle=linestyle, marker=marker, label=f'{tag}G{xy}*J{xz}')
        #ax.plot(station, I[:, 0], 'b-', marker=marker, label='Ixx')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness: G{xy}*J{xz}')
        ax.legend()
        png_filename = dirname / f'stiffness_G{xy}J{xz}.png'
        if save:
            log.info(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1

        #---------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, GJ, 'k', linestyle=linestyle, marker=marker, label=f'{tag}G{xy}*J{xz}')
        ax.plot(station, ExI[:, 0], 'r', linestyle=linestyle, marker=marker, label=f'{tag}E{y}*I{xx}')
        ax.plot(station, EyA, 'b', linestyle=linestyle, marker=marker, label=f'{tag}E{y}*A{y}')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness: G{xy}*J{xz}')
        ax.legend()
        png_filename = dirname / f'stiffness_E{x}I{xx}_G{xy}J{xz}.png'
        if save:
            log.info(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1
        #---------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, G, color='k', linestyle=linestyle, marker=marker, label=f'G{xy}')
        ax.plot(station, Ex, color='r', linestyle=linestyle, marker=marker, label=f'E{y}')  # this is really flipped
        ax.plot(station, Ey, color='b', linestyle=linestyle, marker=marker, label=f'E{x}')  # this is really flipped
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Effective Modulus: E{x}, E{z}, G{xy}')
        ax.legend()
        png_filename = dirname / f'stiffness_E{x}_E{z}_G{xy}.png'
        if save:
            log.info(f'saving {png_filename}')
            fig.savefig(png_filename)
        ifig += 1

        #-------------------------------------------------------
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.plot(station, avg_centroid[:, 0], color='r', linestyle=linestyle, marker=marker, label=f'{tag}xcg')
        ax.plot(station, avg_centroid[:, 2], color='b', linestyle=linestyle, marker=marker, label=f'{tag}zcg')
        ax.grid(True)
        ax.set_xlabel(span_label)
        ax.set_ylabel('Centroid (in)')
        ax.legend()
        png_filename = dirname / 'centroid_vs_span.png'
        if save:
            log.info(f'saving {png_filename}')
            fig.savefig(png_filename)
        #-------------------------------------------------------
        ifig += 4

    #----------------------------------------------------------
    ifig = 1
    if len(data) == 2:
        tag1 = csv_filenames[0][1]
        tag2 = csv_filenames[1][1]
        station1, A1, I1, J1, ExI1, EyI1, GJ1, avg_centroid1 = data[0]
        station2, A2, I2, J2, ExI2, EyI2, GJ2, avg_centroid2 = data[1]
        ustation = np.unique(np.hstack([station1, station2]))
        common_station = np.intersect1d(station1, station2)
        istation1 = np.searchsorted(station1, common_station)
        istation2 = np.searchsorted(station2, common_station)

        # assert np.allclose(station1, station2)
        # Ixx1 = I1[:, 0]
        # Izz1 = I1[:, 2]
        # Ixz1 = I1[:, 5]
        Ex1 = ExI1[:, 0]/I1[:,0]
        # Ey1 = EyI1[:, 0]/I1[:,0]
        # G1 = GJ1/J1
        ExA1 =  Ex1 * A1

        # Ixx2 = I2[:, 0]
        # Izz2 = I2[:, 2]
        # Ixz2 = I2[:, 5]
        Ex2 = ExI2[:, 0]/I2[:,0]
        # Ey2 = EyI2[:, 0]/I2[:,0]
        # G = GJ2/J2
        ExA2 =  Ex2 * A2

        # fig = plt.figure(ifig)
        # ax = fig.gca()

        ax_yscale = 'linear'
        ifig = 10
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.set_yscale(ax_yscale)
        ax2 = ax.twinx()
        # 1 / 2 will give us less
        ExA_ratio = ExA1[istation1] / ExA2[istation2] - 1
        ax.plot(station1, ExA1, color='r', linestyle='-', marker=marker, label=f'{tag1}E{y}*A{y}')  # this is really Ey*Ay; just some bad names upstream
        ax.plot(station2, ExA2, color='b', linestyle='--', marker=marker, label=f'{tag2}E{y}*A{y}')
        ax2.plot(common_station, ExA_ratio*100, color='k', linestyle='-', marker=marker, label='ratio')  # Ey*A
        ax.grid()
        ax.legend()
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness ($lb_f*in^2$): E{y}*A{y}')
        ax2.set_ylabel('%Difference')
        fig.savefig(dirname / f'compare_E{y}A{y}.png')

        ifig += 1
        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.set_yscale(ax_yscale)
        ax2 = ax.twinx()
        # 1 / 2 will give us less
        ax.grid()
        ExI_ratio = ExI1[istation1, 2] / ExI2[istation2, 2] - 1
        ax.plot(station1, ExI1[:, 2], color='r', linestyle='-', marker=marker, label=f'{tag1}E{y}*I{xx}')  # this is really Ey*Ixx; just some bad names upstream
        ax.plot(station2, ExI2[:, 2], color='b', linestyle='--', marker=marker, label=f'{tag2}E{y}*I{xx}')
        ax2.plot(common_station, ExI_ratio*100, color='k', linestyle='-', marker=marker, label='ratio')
        if ylim_EyIzz_ratio:
            ax2.set_ylim(ylim_EyIzz_ratio)
        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness ($lb_f*in^2$): E{y}*I{xx}')
        ax2.set_ylabel('%Difference')
        ax.legend()
        fig.savefig(dirname / f'compare_E{y}I{xx}.png')
        ifig += 1

        fig = plt.figure(ifig)
        ax = fig.gca()
        ax.set_yscale(ax_yscale)
        ax2 = ax.twinx()
        GJ_ratio = GJ1[istation1] / GJ2[istation2] - 1
        ax.plot(station1, GJ1, color='r', linestyle='-', marker=marker, label=f'{tag1}G{xy}*J{xz}')
        ax.plot(station2, GJ2, color='b', linestyle='--', marker=marker, label=f'{tag2}G{xy}*J{xz}')
        ax2.plot(common_station, GJ_ratio*100, color='k', linestyle='-', marker=marker, label='ratio')
        if ylim_GJ_ratio:
            ax2.set_ylim(ylim_GJ_ratio)

        ax.set_xlabel(span_label)
        ax.set_ylabel(f'Stiffness ($lb_f*in^2$): G{xy}*J{xz}')
        ax2.set_ylabel('%Difference')
        ax.grid()
        ax.legend()
        fig.savefig(dirname / f'compare_G{xy}_J{xz}.png')
        plt.show()
    ifig += 1

    if show:
        plt.show()
    return ifig
