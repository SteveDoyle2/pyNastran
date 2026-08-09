"""Tests for cut_body7: concave hull and ZAERO BODY7 generation."""

import os
import unittest
from pathlib import Path

import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cpylog import SimpleLogger
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import BDF, CORD2R

from pyNastran.bdf.mesh_utils.cut.cut_body7 import (
    get_cut_points,
    cut_and_generate_body7,
)

TEST_PATH = Path(__file__).parent
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'
assert PKG_PATH.exists(), print_bad_path(PKG_PATH)
assert MODEL_PATH.exists(), print_bad_path(MODEL_PATH)


class TestGetCutPoints(unittest.TestCase):
    """Tests for get_cut_points."""

    def test_extracts_unique_points(self) -> None:
        """Duplicate xyz values between rod endpoints should be collapsed."""
        rod_eid_nodes = np.array([
            [10, 1, 2],
            [20, 2, 3],
        ], dtype="int32",)
        rod_nids = np.array([1, 2, 3], dtype="int32")
        rod_xyzs = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ])
        rods = (rod_eid_nodes, rod_nids, rod_xyzs)
        pts = get_cut_points(rods)
        assert pts.shape == (3, 3)
        np.testing.assert_array_equal(pts[0], [0.0, 0.0, 0.0])
        np.testing.assert_array_equal(pts[1], [1.0, 0.0, 0.0])
        np.testing.assert_array_equal(pts[2], [2.0, 0.0, 0.0])


class TestCutAndGenerateBody7(unittest.TestCase):
    """Integration tests for cut_and_generate_body7.

    Tolerances
    ----------
    station values : exact float match (no interpolation involved)
    nradial : exact integer match in AEFACT card count
    BODY7/SEGMESH text : substring match for card names
    """
    def test_tube_body7_generation(self) -> None:
        """Cut a cylinder at 3 stations and verify BODY7/SEGMESH/AEFACT output."""
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = build_tube_model(log, radius=5.0, length=20.0, ncircum=16, nspan=5)

        stations = [4.0, 10.0, 16.0]
        coords = []
        for i, dy in enumerate(stations):
            coord = CORD2R(
                100 + i,
                rid=0,
                origin=np.array([0.0, dy, 0.0]),
                zaxis=np.array([0.0, dy, 1.0]),
                xzplane=np.array([1.0, dy, 0.0]),
            )
            coords.append(coord)

        normal_plane = np.array([0.0, 1.0, 0.0])

        result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            stations,
            coords,
            body_id=10,
            label="TUBE",
            nradial=12,
            alpha=0.0,
            segmesh_id_start=200,
            aefact_id_start=5000,
        )

        assert len(result["stations_found"]) == 3
        np.testing.assert_array_equal(result["stations_found"], stations)
        assert len(result["hulls"]) == 3
        assert len(result["hull_yz_resampled"]) == 3
        assert result["centroids"].shape == (3, 3)

        for yz in result["hull_yz_resampled"]:
            assert yz.shape == (12, 2)

        model = result["model"]
        assert len(model.caeros) == 1 #"BODY7" in cards_text
        assert len(model.paeros) == 1 #"SEGMESH" in cards_text
        assert len(model.aefacts) == 6 # "AEFACT" in cards_text
        #assert "TUBE" in cards_text
        #assert cards_text.count("AEFACT") == 6

        plot_cross_sections(
            [result],
            ["R=5 tube"],
            ["blue"],
            stations,
            "test_tube_body7_generation: y-axis tube, 3 stations",
            TEST_PATH / "tmp_tube_body7_generation.png",
            nominal_radii=[5.0],
        )

    def test_tube_body7_with_pbody7(self) -> None:
        """Verify PBODY7 card is generated when pbody7_id > 0."""
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = build_tube_model(
            log, radius=3.0, length=10.0, ncircum=12, nspan=4)

        stations = [2.5, 5.0, 7.5]
        coords = []
        for i, dy in enumerate(stations):
            coord = CORD2R(
                200 + i,
                rid=0,
                origin=np.array([0.0, dy, 0.0]),
                zaxis=np.array([0.0, dy, 1.0]),
                xzplane=np.array([1.0, dy, 0.0]),
            )
            coords.append(coord)

        normal_plane = np.array([0.0, 1.0, 0.0])

        result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            stations,
            coords,
            body_id=20,
            label="TUBE2",
            nradial=8,
            pbody7_id=5,
        )

        zaero_model = result['model']
        assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
        assert len(zaero_model.paeros) == 1, len(zaero_model.paeros)
        assert len(zaero_model.aefacts) == 6, len(zaero_model.aefacts)
        #cards_text = result["cards_text"]
        #assert "PBODY7" in cards_text
        #assert "BODY7" in cards_text

        plot_cross_sections(
            [result],
            ["R=3 tube"],
            ["green"],
            stations,
            "test_tube_body7_with_pbody7: y-axis tube + PBODY7",
            TEST_PATH / "tmp_tube_body7_with_pbody7.png",
            nominal_radii=[3.0],
        )

    def test_output_file_written(self) -> None:
        """Verify that output_filename writes the cards to disk."""
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = build_tube_model(
            log, radius=2.0, length=8.0, ncircum=12, nspan=3)

        stations = [2.0, 4.0, 6.0]
        coords = []
        for i, dy in enumerate(stations):
            coord = CORD2R(
                300 + i,
                rid=0,
                origin=np.array([0.0, dy, 0.0]),
                zaxis=np.array([0.0, dy, 1.0]),
                xzplane=np.array([1.0, dy, 0.0]),
            )
            coords.append(coord)

        normal_plane = np.array([0.0, 1.0, 0.0])
        out_file = TEST_PATH / "tmp_body7_output.dat"

        try:
            result = cut_and_generate_body7(
                model,
                normal_plane,
                log,
                stations,
                coords,
                body_id=30,
                label="FUSE",
                nradial=8,
                output_filename=out_file,
            )

            assert out_file.exists()
            zaero_model = result['model']
            assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
            assert len(zaero_model.paeros) == 1, len(zaero_model.paeros)
            assert len(zaero_model.aefacts) == 6, len(zaero_model.aefacts)
            #contents = out_file.read_text()
            #assert "BODY7" in contents
            #assert "SEGMESH" in contents
            #assert contents == result["cards_text"]

            plot_cross_sections(
                [result],
                ["R=2 tube"],
                ["red"],
                stations,
                "test_output_file_written: y-axis tube, file output",
                TEST_PATH / "tmp_output_file_written.png",
                nominal_radii=[2.0],
            )
        finally:
            if out_file.exists():
                os.remove(out_file)

    def test_no_cuts_found_raises(self) -> None:
        """Stations placed outside the model should raise RuntimeError."""
        log = SimpleLogger(level="error", encoding="utf-8")
        model = build_tube_model(
            log, radius=2.0, length=8.0, ncircum=12, nspan=3)

        stations = [100.0, 200.0]
        coords = []
        for i, dy in enumerate(stations):
            coord = CORD2R(
                400 + i,
                rid=0,
                origin=np.array([0.0, dy, 0.0]),
                zaxis=np.array([0.0, dy, 1.0]),
                xzplane=np.array([1.0, dy, 0.0]),
            )
            coords.append(coord)

        normal_plane = np.array([0.0, 1.0, 0.0])
        with self.assertRaises(RuntimeError):
            cut_and_generate_body7(
                model,
                normal_plane,
                log,
                stations,
                coords,
                body_id=40,
                label="MISS",
                nradial=8,
            )

    def test_single_station(self) -> None:
        """A single cut station should still produce valid output."""
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = build_tube_model(
            log, radius=4.0, length=12.0, ncircum=16, nspan=4)

        stations = [6.0]
        coords = [
            CORD2R(
                500,
                rid=0,
                origin=np.array([0.0, 6.0, 0.0]),
                zaxis=np.array([0.0, 6.0, 1.0]),
                xzplane=np.array([1.0, 6.0, 0.0]),
            )
        ]

        normal_plane = np.array([0.0, 1.0, 0.0])
        result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            stations,
            coords,
            body_id=50,
            label="ONESTN",
            nradial=8,
        )
        zaero_model = result['model']
        assert len(result["stations_found"]) == 1
        assert len(result["hulls"]) == 1
        #print(zaero_model.get_bdf_stats())
        assert len(zaero_model.caeros) == 1, len(zaero_model.caeros) # "BODY7" in result["cards_text"]
        assert result["hull_yz_resampled"][0].shape == (8, 2)

        plot_cross_sections(
            [result],
            ["R=4 tube"],
            ["purple"],
            stations,
            "test_single_station: 1 cut on y-axis tube",
            TEST_PATH / "tmp_single_station.png",
            nominal_radii=[4.0],
        )


class TestTwoCylindersXAxis(unittest.TestCase):
    """Two x-axis cylinders: a fuselage and an offset nacelle.

    Tests the cutting convention with marching along x instead of y.
    Each cylinder gets its own BODY7.

    Tolerances
    ----------
    station values : exact float match
    nradial : exact integer match in AEFACT card count
    hull shape : cross-section radius within 10% of nominal (mesh discretization)
    """

    def test_two_x_cylinders(self) -> None:
        """Cut a fuselage + nacelle model and generate two BODY7 cards.

        The x-ranges are disjoint so each body's cut stations only
        intersect that body's elements (cut_face_model_by_coord cuts
        ALL elements in the model at the cut plane).

        Fuselage: radius=3, y=0 z=0, x=[0..20], 8 bays
        Nacelle:  radius=1, y=5 z=0, x=[25..35], 4 bays
        Stations:  fuselage at x=[2.5, 7.5, 12.5, 17.5]
                   nacelle  at x=[27, 30, 33]
        """
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = BDF(log=log, debug=False)
        model.add_mat1(1, E=1.0e7, G=None, nu=0.3)
        model.add_pshell(1, 1, t=0.1)

        next_nid, next_eid = add_x_cylinder(
            model,
            nid_start=1,
            eid_start=1,
            pid=1,
            center_y=0.0,
            center_z=0.0,
            radius=3.0,
            x_start=0.0,
            x_end=20.0,
            ncircum=16,
            nspan=8,
        )
        add_x_cylinder(
            model,
            nid_start=next_nid,
            eid_start=next_eid,
            pid=1,
            center_y=5.0,
            center_z=0.0,
            radius=1.0,
            x_start=25.0,
            x_end=35.0,
            ncircum=12,
            nspan=4,
        )
        model.cross_reference()

        normal_plane = np.array([1.0, 0.0, 0.0])

        # --- fuselage BODY7 ---
        fuse_stations = [2.5, 7.5, 12.5, 17.5]
        fuse_coords = [
            make_x_cut_coord(100 + i, dx, center_y=0.0, center_z=0.0)
            for i, dx in enumerate(fuse_stations)
        ]
        fuse_result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            fuse_stations,
            fuse_coords,
            body_id=1,
            label="FUSELAG",
            nradial=16,
            alpha=0.0,
            segmesh_id_start=100,
            aefact_id_start=1000,
        )

        assert len(fuse_result["stations_found"]) == 4
        np.testing.assert_array_equal(fuse_result["stations_found"], fuse_stations)
        for yz in fuse_result["hull_yz_resampled"]:
            assert yz.shape == (16, 2)
            radii = np.sqrt(yz[:, 0] ** 2 + yz[:, 1] ** 2)
            np.testing.assert_allclose(radii, 3.0, atol=0.5)

        zaero_model = fuse_result['model']
        assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
        assert len(zaero_model.paeros) == 1, len(zaero_model.paeros)
        assert len(zaero_model.aefacts) == 8, len(zaero_model.aefacts)
        #assert "BODY7" in fuse_result["cards_text"]
        #assert "FUSELAG" in fuse_result["cards_text"]
        #assert "SEGMESH" in fuse_result["cards_text"]
        #assert fuse_result["cards_text"].count("AEFACT") == 8

        # --- nacelle BODY7 ---
        nac_stations = [27.0, 30.0, 33.0]
        nac_coords = [
            make_x_cut_coord(200 + i, dx, center_y=5.0, center_z=0.0)
            for i, dx in enumerate(nac_stations)
        ]
        nac_result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            nac_stations,
            nac_coords,
            body_id=2,
            label="NACELLE",
            nradial=12,
            alpha=0.0,
            segmesh_id_start=200,
            aefact_id_start=2000,
        )

        assert len(nac_result["stations_found"]) == 3
        np.testing.assert_array_equal(nac_result["stations_found"], nac_stations)
        for yz in nac_result["hull_yz_resampled"]:
            assert yz.shape == (12, 2)
            radii = np.sqrt(yz[:, 0] ** 2 + yz[:, 1] ** 2)
            np.testing.assert_allclose(radii, 1.0, atol=0.25)

        zaero_model = fuse_result['model']
        assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
        assert len(zaero_model.paeros) == 1, len(zaero_model.paeros)
        assert len(zaero_model.aefacts) == 8, len(zaero_model.aefacts)
        #assert "BODY7" in nac_result["cards_text"]
        #assert "NACELLE" in nac_result["cards_text"]
        #assert nac_result["cards_text"].count("AEFACT") == 6

        plot_cross_sections(
            [fuse_result],
            ["fuselage R=3"],
            ["blue"],
            fuse_stations,
            "test_two_x_cylinders: fuselage (x-axis)",
            TEST_PATH / "tmp_two_x_cyl_fuselage.png",
            nominal_radii=[3.0],
        )
        plot_cross_sections(
            [nac_result],
            ["nacelle R=1"],
            ["red"],
            nac_stations,
            "test_two_x_cylinders: nacelle (x-axis)",
            TEST_PATH / "tmp_two_x_cyl_nacelle.png",
            nominal_radii=[1.0],
        )

    def test_two_x_cylinders_write_output(self) -> None:
        """Write both BODY7 card sets to a single file."""
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = BDF(log=log, debug=False)
        model.add_mat1(1, E=1.0e7, G=None, nu=0.3)
        model.add_pshell(1, 1, t=0.1)

        next_nid, next_eid = add_x_cylinder(
            model,
            nid_start=1,
            eid_start=1,
            pid=1,
            center_y=0.0,
            center_z=0.0,
            radius=3.0,
            x_start=0.0,
            x_end=20.0,
            ncircum=16,
            nspan=8,
        )
        add_x_cylinder(
            model,
            nid_start=next_nid,
            eid_start=next_eid,
            pid=1,
            center_y=5.0,
            center_z=0.0,
            radius=1.0,
            x_start=25.0,
            x_end=35.0,
            ncircum=12,
            nspan=4,
        )
        model.cross_reference()

        normal_plane = np.array([1.0, 0.0, 0.0])
        out_file = TEST_PATH / "tmp_two_cylinders.dat"

        try:
            fuse_stations = [5.0, 10.0, 15.0]
            fuse_coords = [
                make_x_cut_coord(300 + i, dx) for i, dx in enumerate(fuse_stations)
            ]
            fuse_result = cut_and_generate_body7(
                model,
                normal_plane,
                log,
                fuse_stations,
                fuse_coords,
                body_id=1,
                label="FUSELAG",
                nradial=12,
                segmesh_id_start=100,
                aefact_id_start=1000,
                output_filename=out_file,
            )
            zaero_model = fuse_result['model']
            assert out_file.exists()
            contents = out_file.read_text()
            assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
            #assert "FUSELAG" in contents
            #assert contents == fuse_result["cards_text"]

            nac_stations = [27.5, 31.25]
            nac_coords = [
                make_x_cut_coord(400 + i, dx, center_y=5.0)
                for i, dx in enumerate(nac_stations)
            ]
            nac_result = cut_and_generate_body7(
                model,
                normal_plane,
                log,
                nac_stations,
                nac_coords,
                body_id=2,
                label="NACELLE",
                nradial=10,
                segmesh_id_start=200,
                aefact_id_start=2000,
                zaero_model = zaero_model,
            )

            #with open(out_file, "a") as f:
            #    f.write(nac_result["cards_text"])
            #model = nac_result['model']

            combined = out_file.read_text()
            assert len(zaero_model.caeros) == 2
            assert len(zaero_model.paeros) == 2
            assert "FUSELAG" in combined
            assert "NACELLE" in combined
            body7_lines = [ln for ln in combined.splitlines() if ln.startswith("BODY7")]
            segmesh_lines = [ln for ln in combined.splitlines() if ln.startswith("SEGMESH")]
            assert len(body7_lines) == 2
            assert len(segmesh_lines) == 2

            plot_cross_sections(
                [fuse_result],
                ["fuselage R=3"],
                ["blue"],
                fuse_stations,
                "test_write_output: fuselage",
                TEST_PATH / "tmp_write_fuse.png",
                nominal_radii=[3.0],
            )
            plot_cross_sections(
                [nac_result],
                ["nacelle R=1"],
                ["red"],
                nac_stations,
                "test_write_output: nacelle",
                TEST_PATH / "tmp_write_nac.png",
                nominal_radii=[1.0],
            )
        finally:
            if out_file.exists():
                os.remove(out_file)

    def test_two_cylinders_four_stations_plot(self) -> None:
        """Two concentric-ish cylinders in one model, cut as a single body.

        The convex hull at each station wraps both cross-sections into
        one outer boundary, discarding interior points.

        Cylinder 1: R=5,  center y=0  z=0, x=[-11..101], ncircum=20, nspan=1
        Cylinder 2: R=10, center y=0  z=0, x=[-11..101], ncircum=20, nspan=1
        Stations: x = 0, 25, 50, 100
        Coord origin at (x_station, 0, 0) so both cylinders are centered.

        The hull should match the outer cylinder (R=10) since it
        encloses the inner one (R=5).

        Tolerances
        ----------
        station count : 4 exact
        hull radii : within 5% of R=10 (outer cylinder dominates)
        """
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = BDF(log=log, debug=False)
        model.add_mat1(1, E=1.0e7, G=None, nu=0.3)
        model.add_pshell(1, 1, t=0.1)

        r1 = 5.0
        r2 = 10.0

        next_nid, next_eid = add_x_cylinder(
            model,
            nid_start=1,
            eid_start=1,
            pid=1,
            center_y=0.0,
            center_z=0.0,
            radius=r1,
            x_start=-11.0,
            x_end=101.0,
            ncircum=20,
            nspan=1,
        )
        add_x_cylinder(
            model,
            nid_start=next_nid,
            eid_start=next_eid,
            pid=1,
            center_y=0.0,
            center_z=0.0,
            radius=r2,
            x_start=-11.0,
            x_end=101.0,
            ncircum=20,
            nspan=1,
        )
        model.cross_reference()

        normal_plane = np.array([1.0, 0.0, 0.0])
        stations = [0.0, 25.0, 50.0, 100.0]
        coords = [
            make_x_cut_coord(100 + i, dx, center_y=0.0, center_z=0.0)
            for i, dx in enumerate(stations)
        ]

        result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            stations,
            coords,
            body_id=1,
            label="COMBINED",
            nradial=20,
            alpha=0.0,
            segmesh_id_start=100,
            aefact_id_start=1000,
        )

        assert len(result["stations_found"]) == 4
        np.testing.assert_array_equal(result["stations_found"], stations)

        for yz in result["hull_yz_resampled"]:
            assert yz.shape == (20, 2)
            radii = np.sqrt(yz[:, 0] ** 2 + yz[:, 1] ** 2)
            np.testing.assert_allclose(radii, r2, rtol=0.05)

        zaero_model = result['model']
        assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
        assert len(zaero_model.paeros) == 1, len(zaero_model.paeros)
        assert len(zaero_model.aefacts) == 8, len(zaero_model.aefacts)
        #assert "BODY7" in result["cards_text"]
        #assert "COMBINED" in result["cards_text"]
        #assert result["cards_text"].count("AEFACT") == 8

        plot_cross_sections(
            [result],
            ["combined hull"],
            ["red"],
            stations,
            "Hull wraps outer cylinder (R=10), inner points (R=5) discarded",
            TEST_PATH / "tmp_two_cylinders_4stations.png",
            nominal_radii=[r1, r2],
        )

    def test_bwb(self) -> None:
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        normal_plane = None
        log = SimpleLogger(level="warning", encoding="utf-8")
        stations = np.linspace(-100., 1600., num=101)
        coords = [
            make_x_cut_coord(100 + i, dx, center_y=0.0, center_z=0.0)
            for i, dx in enumerate(stations)
        ]

        result = cut_and_generate_body7(
            bdf_filename,
            normal_plane,
            log,
            stations,
            coords,
            body_id=1,
            label="BWB",
            nradial=20,
            alpha=0.0,
            segmesh_id_start=100,
            aefact_id_start=1000,
        )

        zaero_model = result["model"]
        # body7 = zaero_model.zaero.body7[1]
        body7 = zaero_model.caeros[1]
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        body7.plot(ax)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("BWB BODY7 3D")
        fig.tight_layout()
        plt.show()

        plot_cross_sections(
            [result],
            ["BWB hull"],
            ["blue"],
            stations,
            "test_bwb: BWB cross-sections",
            TEST_PATH / "tmp_bwb_xsec.png",
        )

    def test_sears_haack_body(self) -> None:
        """Sears-Haack body of revolution, 10 cut stations.

        R(x) = R_max * [4*x*(1-x)]^(3/4), x in [0,1].
        Minimum wave drag for given length and volume.

        Body: L=100, R_max=5, x-axis, 24 circumferential, 40 axial.
        10 evenly spaced interior stations from x=5 to x=95.

        Tolerances
        ----------
        hull radii : within 6% of analytic R(x) (24-gon inscribed polygon)
        station count : 10 exact
        """
        log = SimpleLogger(level="warning", encoding="utf-8")
        model = BDF(log=log, debug=False)
        model.add_mat1(1, E=1.0e7, G=None, nu=0.3)
        model.add_pshell(1, 1, t=0.1)

        body_length = 100.0
        r_max = 5.0
        ncircum = 24
        nspan = 40

        theta_arr = np.linspace(0, 2 * np.pi, ncircum, endpoint=False)
        x_arr = np.linspace(0.0, body_length, nspan + 1)

        nid = 1
        for xval in x_arr:
            x_norm = xval / body_length
            x_norm = np.clip(x_norm, 1e-12, 1.0 - 1e-12)
            r = sears_haack_radius(r_max, x_norm)
            for th in theta_arr:
                y = r * np.cos(th)
                z = r * np.sin(th)
                model.add_grid(nid, [xval, y, z])
                nid += 1

        eid = 1
        for ix in range(nspan):
            for it in range(ncircum):
                n1 = 1 + ix * ncircum + it
                n2 = 1 + ix * ncircum + (it + 1) % ncircum
                n3 = 1 + (ix + 1) * ncircum + (it + 1) % ncircum
                n4 = 1 + (ix + 1) * ncircum + it
                model.add_cquad4(eid, 1, [n1, n2, n3, n4])
                eid += 1
        model.cross_reference()

        normal_plane = np.array([1.0, 0.0, 0.0])
        stations = np.linspace(5.0, 95.0, 10).tolist()
        coords = [
            make_x_cut_coord(100 + i, dx, center_y=0.0, center_z=0.0)
            for i, dx in enumerate(stations)]

        result = cut_and_generate_body7(
            model,
            normal_plane,
            log,
            stations,
            coords,
            body_id=1,
            label="SEARSHK",
            nradial=24,
            alpha=0.0,
            segmesh_id_start=100,
            aefact_id_start=1000,
        )

        assert len(result["stations_found"]) == 10
        np.testing.assert_array_equal(result["stations_found"], stations)

        analytic_radii = []
        for dx in stations:
            x_norm = dx / body_length
            r_expected = sears_haack_radius(r_max, x_norm)
            analytic_radii.append(r_expected)

        for i, yz in enumerate(result["hull_yz_resampled"]):
            assert yz.shape == (24, 2)
            radii = np.sqrt(yz[:, 0] ** 2 + yz[:, 1] ** 2)
            np.testing.assert_allclose(
                radii,
                analytic_radii[i],
                rtol=0.06,
                err_msg=f"station {stations[i]}: expected R={analytic_radii[i]:.3f}",
            )

        zaero_model = result['model']
        assert len(zaero_model.caeros) == 1, len(zaero_model.caeros)
        assert len(zaero_model.paeros) == 1, len(zaero_model.paeros)
        assert len(zaero_model.aefacts) == 20, len(zaero_model.aefacts)
        #assert "BODY7" in result["cards_text"]
        #assert "SEARSHK" in result["cards_text"]
        #assert result["cards_text"].count("AEFACT") == 20

        plot_cross_sections(
            [result],
            ["Sears-Haack hull"],
            ["blue"],
            stations,
            f"Sears-Haack body: L={body_length}, R_max={r_max}, 10 stations",
            TEST_PATH / "tmp_sears_haack_xsec.png",
            nominal_radii=[r_max],
        )

        # --- side-view profile: R(x) analytic vs. hull mean radius ---
        hull_mean_radii = []
        for yz in result["hull_yz_resampled"]:
            radii = np.sqrt(yz[:, 0] ** 2 + yz[:, 1] ** 2)
            hull_mean_radii.append(radii.mean())

        x_fine = np.linspace(0, body_length, 500)
        r_fine = [sears_haack_radius(x / body_length) for x in x_fine]

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.fill_between(x_fine, r_fine, -np.array(r_fine), color="lightblue", alpha=0.3)
        ax.plot(x_fine, r_fine, "b-", linewidth=1.5, label="analytic R(x)")
        ax.plot(x_fine, -np.array(r_fine), "b-", linewidth=1.5)

        ax.plot(
            stations,
            hull_mean_radii,
            "ro",
            markersize=6,
            label="hull mean radius",
        )
        ax.plot(stations, [-r for r in hull_mean_radii], "ro", markersize=6)

        for i, dx in enumerate(stations):
            ax.plot(
                [dx, dx],
                [-hull_mean_radii[i], hull_mean_radii[i]],
                "r-",
                linewidth=0.8,
                alpha=0.5,
            )

        ax.set_xlabel("x")
        ax.set_ylabel("R(x)")
        ax.set_title(
            f"Sears-Haack body profile: L={body_length}, R_max={r_max}, "
            f"R(x) = R_max * [4x(1-x)]^{{3/4}}"
        )
        ax.set_aspect("equal")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        profile_path = TEST_PATH / "tmp_sears_haack_profile.png"
        fig.savefig(profile_path, dpi=150)
        plt.show()
        # plt.close(fig)
        assert profile_path.exists()
        assert profile_path.stat().st_size > 0
        os.remove(profile_path)


def sears_haack_radius(r_max: float, x_norm: float) -> float:
    """R(x) = R_max * [4*x*(1-x)]^(3/4), x in [0,1]."""
    return r_max * (4.0 * x_norm * (1.0 - x_norm)) ** 0.75


def make_x_cut_coord(
    cid: int,
    x_station: float,
    center_y: float = 0.0,
    center_z: float = 0.0,) -> CORD2R:
    """Build a CORD2R whose local y-axis = global x.

    The cut happens at local y = 0, which corresponds to
    global x = x_station.  Local x = -global y, local z = global z.

    Parameters
    ----------
    cid : int
        coordinate system ID
    x_station : float
        where to cut along global x
    center_y, center_z : float
        cross-section center (shifts the coord origin so the hull
        is centered on the body axis)
    """
    origin = np.array([x_station, center_y, center_z])
    zaxis = origin + np.array([0.0, 0.0, 1.0])
    xzplane = origin + np.array([0.0, -1.0, 0.0])
    return CORD2R(cid, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)


def build_tube_model(log: SimpleLogger,
                     radius: float = 5.0,
                     length: float = 20.0,
                     ncircum: int = 16,
                     nspan: int = 5,) -> BDF:
    """Build a cylindrical tube shell model aligned along y-axis.

    Grid layout: ncircum points around circumference (xz-plane)
    at each of nspan+1 stations from y=0 to y=length. CQUAD4
    elements connect adjacent rings. The cutting convention
    marches along y.
    """
    model = BDF(log=log, debug=False)

    mid = 1
    pid = 1
    t = 0.1
    E = 1.0e7
    model.add_mat1(mid, E=E, G=None, nu=0.3)
    model.add_pshell(pid, mid, t=t)

    nid = 1
    theta_arr = np.linspace(0, 2 * np.pi, ncircum, endpoint=False)
    y_arr = np.linspace(0.0, length, nspan + 1)

    for iy, yval in enumerate(y_arr):
        for it, th in enumerate(theta_arr):
            x = radius * np.cos(th)
            z = radius * np.sin(th)
            model.add_grid(nid, [x, yval, z])
            nid += 1

    eid = 1
    for iy in range(nspan):
        for it in range(ncircum):
            n1 = iy * ncircum + it + 1
            n2 = iy * ncircum + (it + 1) % ncircum + 1
            n3 = (iy + 1) * ncircum + (it + 1) % ncircum + 1
            n4 = (iy + 1) * ncircum + it + 1
            model.add_cquad4(eid, pid, [n1, n2, n3, n4])
            eid += 1

    model.cross_reference()
    return model


def add_x_cylinder(model: BDF,
                   nid_start: int,
                   eid_start: int,
                   pid: int,
                   center_y: float,
                   center_z: float,
                   radius: float,
                   x_start: float,
                   x_end: float,
                   ncircum: int = 16,
                   nspan: int = 8,) -> tuple[int, int]:
    """
    Add a cylinder along x-axis to an existing BDF.

    Parameters
    ----------
    model : BDF
        model to add grids/elements to (already has MAT1/PSHELL)
    nid_start, eid_start : int
        starting node/element IDs
    pid : int
        PSHELL property ID (must already exist in model)
    center_y, center_z : float
        cross-section center offset from origin
    radius : float
        cylinder radius
    x_start, x_end : float
        axial extent
    ncircum, nspan : int
        mesh density

    Returns
    -------
    next_nid, next_eid : int
        next available node/element IDs after this cylinder
    """
    theta_arr = np.linspace(0, 2 * np.pi, ncircum, endpoint=False)
    x_arr = np.linspace(x_start, x_end, nspan + 1)

    nid = nid_start
    for ix, xval in enumerate(x_arr):
        for it, th in enumerate(theta_arr):
            y = center_y + radius * np.cos(th)
            z = center_z + radius * np.sin(th)
            model.add_grid(nid, [xval, y, z])
            nid += 1

    eid = eid_start
    for ix in range(nspan):
        for it in range(ncircum):
            n1 = nid_start + ix * ncircum + it
            n2 = nid_start + ix * ncircum + (it + 1) % ncircum
            n3 = nid_start + (ix + 1) * ncircum + (it + 1) % ncircum
            n4 = nid_start + (ix + 1) * ncircum + it
            model.add_cquad4(eid, pid, [n1, n2, n3, n4])
            eid += 1
    return nid, eid

def plot_cross_sections(
    results: list[dict],
    labels: list[str],
    colors: list[str],
    stations: list[float],
    title: str,
    plot_path: Path,
    nominal_radii: list[float] | None = None,) -> None:
    """Save a cross-section plot for one or more cut results.

    Parameters
    ----------
    results : list[dict]
        each entry is a return dict from cut_and_generate_body7
    labels : list[str]
        legend label per result
    colors : list[str]
        matplotlib color per result
    stations : list[float]
        x-axis station values (subplot titles)
    title : str
        figure suptitle
    plot_path : Path
        where to save the .png
    nominal_radii : list[float] or None
        if given, overlay dashed reference circles at these radii
    """
    nstations = len(stations)
    fig, axes = plt.subplots(1, nstations, figsize=(4 * nstations, 4))
    if nstations == 1:
        axes = [axes]

    for istation, ax in enumerate(axes):
        for ires, (res, lbl, clr) in enumerate(zip(results, labels, colors)):
            if istation >= len(res["hulls"]):
                continue
            raw = res["cut_points_2d"][istation]
            hull = res["hulls"][istation]
            yz = res["hull_yz_resampled"][istation]

            ax.plot(raw[:, 0], raw[:, 1], ".", color=clr, markersize=2, alpha=0.4)
            ax.plot(hull[:, 0], hull[:, 1], "-", color=clr, linewidth=0.6, alpha=0.5)
            ax.plot(yz[:, 0], yz[:, 1], "o-", color=clr, markersize=3, label=lbl)

        if nominal_radii:
            theta = np.linspace(0, 2 * np.pi, 100)
            for r in nominal_radii:
                ax.plot(
                    r * np.cos(theta),
                    r * np.sin(theta),
                    "--",
                    color="gray",
                    linewidth=0.7,
                    alpha=0.5,
                    label=f"R={r:.1f}",
                )

        ax.set_title(f"station = {stations[istation]:.1f}")
        ax.set_xlabel("local y")
        ax.set_ylabel("local z")
        ax.set_aspect("equal")
        ax.legend(fontsize=6, loc="upper right")
        ax.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    fig.show()
    # plt.close(fig)
    assert plot_path.exists()
    assert plot_path.stat().st_size > 0
    os.remove(plot_path)


if __name__ == "__main__":
    unittest.main()
