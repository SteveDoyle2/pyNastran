from pathlib import Path
import unittest
import numpy as np

from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent

from pyNastran.dev.tools.pressure_map import (
    pressure_map, pressure_filename_to_fa2j, pressure_filename_to_wkk_diag)
from pyNastran.dev.tools.pressure_map_aero_setup import (
    get_aero_model, get_aero_pressure_centroid)

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_DIR = PKG_PATH / '..' / 'models'
DIRNAME = Path(__file__).parent


class TestPressureMap(unittest.TestCase):
    def test_pressure_map_cart3d(self):
        aero_format = 'cart3d'
        cart3d_filename = PKG_PATH / 'converters' / 'cart3d' / 'models' / 'threePlugs.bin.tri'
        bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.bdf'
        caero_bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.caero.bdf'

        skip_cards = ['CBAR']
        bdf_model = read_bdf(bdf_filename, skip_cards=skip_cards)
        if not caero_bdf_filename.exists():  # pragma: no cover
            export_caero_mesh(
                bdf_filename, caero_bdf_filename,
                is_aerobox_model=True,
                write_panel_xyz=False)

        aero_model, variables = get_aero_model(
            cart3d_filename, aero_format,
                   aero_xyz_scale=1.0,
                   xyz_units='in',
                   stop_on_failure=True)

        neids = len(aero_model.elements)
        eids = np.arange(neids)
        Cp = np.sin(eids/neids)
        aero_model.loads['Cp'] = Cp

        pressure_map(
            aero_model,
            bdf_model,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='pressure',
            method='full_model',
            xyz_units='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'cart3d_pressure_fullmodel_1.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None)

        pressure_map(
            aero_model,
            bdf_model,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force',
            method='full_model',
            xyz_units='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'cart3d_force_fullmodel_2.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None)

        with self.assertRaises(RuntimeError):
            pressure_map(
                aero_model,
                bdf_model,
                # eids_structure=np.array([]),
                # eid_csv_filename='',
                eid_load_id=-1,
                aero_format=aero_format,
                map_type='force_moment',
                method='full_model',
                xyz_units='in',
                pressure_units='psi',
                pressure_sid=1,
                force_sid=2,
                moment_sid=3,
                idtype='int32', fdtype='float64',
                pressure_filename=DIRNAME/'cart3d_forcemoment_fullmodel_3.bdf',
                aero_xyz_scale=1.0, qinf=1.0,
                sref=1.0, cref=1.0, bref=1.0,
                reference_point=None,
                regions_to_include=None,
                regions_to_remove=None)

        pressure_filename = DIRNAME / 'cart3d_forcemoment_panelmodel_4.bdf'
        pressure_map(
            aero_model, #cart3d_filename,
            caero_bdf_filename,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force_moment',
            method='panel_model',
            xyz_units='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=pressure_filename,
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None)

        fa2j_filename = DIRNAME/'cart3d_fa2j_5.bdf'
        pressure_filename_to_fa2j(pressure_filename, fa2j_filename, sid=1)

        wkk_filename = DIRNAME/'cart3d_wkk_6.bdf'
        pressure_filename1 = pressure_filename
        pressure_filename2 = pressure_filename
        pressure_filename_to_wkk_diag(
            pressure_filename1,
            pressure_filename2,
            np.zeros(2),
            np.ones(2),
            wkk_filename,
            force_sid1=2, moment_sid1=3,
            force_sid2=1, moment_sid2=3)

    def test_pressure_map_fluent(self):
        aero_format = 'fluent'
        # map_type = 'pressure'
        bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.bdf'
        vrt_filename = MODEL_DIR / 'bwb' / 'bwb-saero.vrt'
        caero_bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.caero.bdf'
        if not caero_bdf_filename.exists():  # pragma: no cover
            export_caero_mesh(
                bdf_filename, caero_bdf_filename,
                is_aerobox_model=True,
                write_panel_xyz=False)

        log = SimpleLogger(level='info')
        if not vrt_filename.exists():  # pragma: no cover
            nastran_to_fluent(bdf_filename, vrt_filename, log=log)

        aero_model, variables = get_aero_model(
            vrt_filename, aero_format,
                   aero_xyz_scale=1.0,
                   xyz_units='in',
                   stop_on_failure=True)
        aero_model.titles = ['ElementID', 'Pressure Coefficient']
        # get_aero_pressure_centroid(
        #     aero_model, aero_format,
        #     map_type, variable='Cp',
        #     regions_to_include=None,
        #     regions_to_remove=None)
        pressure_map(
            aero_model, #cart3d_filename,
            caero_bdf_filename,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force_moment',
            method='panel_model',
            xyz_units='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'fluent_forcemoment_panelmodel_1.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None)

        pressure_map(
            aero_model, #cart3d_filename,
            caero_bdf_filename,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force',
            method='full_model',
            xyz_units='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'fluent_force_fullmodel_2.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None)


if __name__ == '__main__':
    unittest.main()
