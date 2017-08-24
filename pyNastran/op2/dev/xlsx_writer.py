#pylint: disable=W0201,C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import copy
from datetime import date
from collections import defaultdict
from typing import Any
from six import iteritems

import pyNastran
from pyNastran.op2.op2_interface.op2_f06_common import OP2_F06_Common
from pyNastran.op2.op2_interface.write_utils import _write_markers

def make_stamp(title, today=None):
    # (str, Any) -> str
    if title is None:
        title = ''

    #lenghts = [7, 8, 5, 5, 3, 4, 4, 6, 9, 7, 8, 8]
    months = [' January', 'February', 'March', 'April', 'May', 'June',
              'July', 'August', 'September', 'October', 'November', 'December']
    if today is None:
        today = date.today()
        str_month = months[today.month - 1].upper()
        str_today = '%-9s %2s, %4s' % (str_month, today.day, today.year)
    else:
        (month, day, year) = today
        str_month = months[month - 1].upper()
        str_today = '%-9s %2s, %4s' % (str_month, day, year)
    str_today = str_today  #.strip()

    release_date = '02/08/12'  # pyNastran.__releaseDate__
    release_date = ''
    build = 'pyNastran v%s %s' % (pyNastran.__version__, release_date)
    if title is None:
        title = ''
    out = '1    %-67s   %-19s %-22s PAGE %%5i\n' % (title.strip(), str_today, build)
    return out


class XlsxWriter(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)
        self.card_count = {}

    def make_f06_header(self):
        """If this class is inherited, the F06 Header may be overwritten"""
        return make_f06_header()

    def make_stamp(self, title, today):
        # (str, Any) -> str
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return make_stamp(title, today)

    def write_xlsx(self, xlsx_filename, is_mag_phase=False,
                  delete_objects=True):
        # type: (str, bool, bool) -> None
        """
        Writes an XLSX file based on the data we have stored in the object

        Parameters
        ----------
        xlsx_filename:  str
            the name of the XLSX file to write
        is_mag_phase : bool;  default=False
            True : write complex using Magnitude/Phase
            False : write complex using Real/Imaginary
            Real objects don't use this parameter.
        delete_objects : bool; default=True
            should objects be deleted after they're written to reduce
            memory
        """
        print('writing %s' % xlsx_filename)

        if isinstance(xlsx_filename, str):
            from xlwings import Workbook, Sheet  # type: ignore
            from pywintypes import com_error  # type: ignore

            workbook = Workbook()  # Creates a connection with a new workbook
            try:
                workbook.save(xlsx_filename)
            except com_error:
                raise RuntimeError('Close %s' % xlsx_filename)
        else:
            raise NotImplementedError(type(xlsx_filename))
            # assert isinstance(xlsx_filename, file), 'type(xlsx_filename)= %s' % xlsx_filename
            # op2 = xlsx_filename
            # xlsx_filename = op2.name
            # print('xlsx_filename =', xlsx_filename)

        #op2_ascii.write('writing [3, 7, 0] header\n')
        #if markers == [3,]:  # PARAM, POST, -1
            #self.read_markers([3])
            #data = self.read_block()

            #self.read_markers([7])
            #data = self.read_block()
            #data = self._read_record()
            #self.read_markers([-1, 0])
        #elif markers == [2,]:  # PARAM, POST, -2
        isheet = 1
        #if 0:
        ##_write_markers(op2, op2_ascii, [3, 0, 7])
            ##tape_code = b'NASTRAN FORT TAPE ID CODE - '
            #Sheet(isheet).name = sheet_name
            #sheet = Sheet(isheet)
            #sheet['A1'].value = 'NASTRAN FORT TAPE ID CODE'
            #sheet['B1'].value = tape_code

            #nastran_version = b'NX8.5   ' if self.is_nx else b'XXXXXXXX'
            #sheet['A2'].value = 'nastran_version'
            #sheet['B2'].value = nastran_version
            #isheet =+ 1

        if self.grid_point_weight.reference_point is not None and 0:
            if hasattr(result, 'write_xlsx'):
                self.grid_point_weight.write_xlsx(workbook, page_stamp, self.page_num)
            else:
                print("*op2 - grid_point_weight not written")


        #is_mag_phase = False

        # eigenvalues are written first
        for ikey, result in sorted(iteritems(self.eigenvalues)):
            # header
            #print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
            if hasattr(result, 'write_xlsx'):
                result.write_xlsx(xlsx)
                if delete_objects:
                    del result
            else:
                print("*xlsx - %s not written" % result.__class__.__name__)
                # asdf

        # then eigenvectors
        # has a special header
        for isubcase, result in sorted(iteritems(self.eigenvectors)):
            (subtitle, label) = self.isubcase_name_map[isubcase]

            if hasattr(result, 'write_xlsx'):
                print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
                result.write_xlsx(workbook, is_mag_phase=is_mag_phase)
                if delete_objects:
                    del result
            else:
                print("*xlsx - %s not written" % result.__class__.__name__)
                # asdf

        # finally, we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order


        #if markers == [3,]:  # PARAM, POST, -2
            #self.read_markers([3])
            #data = self.read_block()
            #self.read_markers([7])
            #data = self.read_block()
            #data = self._read_record()


        oef = [
            # OEF - forces
            # alphabetical order...
            self.cbar_force,

            # beam
            #self.cbeam_forces,
            #self.cbar100_forces,
            #self.cbend_forces,
            self.cbeam_force,
            self.cbush_force,

            self.celas1_force,
            self.celas2_force,
            self.celas3_force,
            self.celas4_force,

            self.cdamp1_force,
            self.cdamp2_force,
            self.cdamp3_force,
            self.cdamp4_force,

            #self.plateForces,   # centroidal elements
            #self.plateForces2,  # bilinear elements

            self.conrod_force,
            self.cquad4_force,
            self.cquad8_force,
            self.crod_force,
            self.cshear_force,
            self.ctria3_force,
            self.ctria6_force,
            self.ctube_force,

            # other
            #self.cgap_force, #self.solidPressureForces,
        ]

        oes1x1 = [
            # cbars/cbeams
            self.cbar_stress, self.cbeam_stress,

            # bush
            self.cbush_stress, self.cbush1d_stress_strain,

            # rods
            #self.nonlinearRodStress,

            self.celas1_stress, self.celas2_stress, self.celas3_stress, self.celas4_stress,

            self.chexa_stress,
            self.conrod_stress,
            self.cpenta_stress,
            self.cquad4_stress, self.cquad8_stress, self.cquadr_stress,
            self.crod_stress,
            self.cshear_stress,
            self.ctetra_stress,
            self.ctria3_stress, self.ctria6_stress, self.ctriar_stress,
            self.ctube_stress,
        ]
        oes1c = [
            self.cquad4_composite_stress, self.cquad8_composite_stress, self.cquadr_composite_stress,
            self.ctria3_composite_stress, self.ctria6_composite_stress, self.ctriar_composite_stress,

            #self.nonlinearPlateStress,
            self.ctriax_stress, #self.hyperelasticPlateStrain,
        ]
        #stress = oes1x1 + oes1c

        strain = [
            # bars/beams
            self.cbar_strain, self.cbeam_strain,

            # springs,
            self.celas1_strain, self.celas2_strain, self.celas3_strain, self.celas4_strain,

            # plates
            self.ctria3_strain, self.cquad4_strain,
            self.cshear_strain,
            self.cquad4_composite_strain, self.cquad8_composite_strain, self.cquadr_composite_strain,
            self.ctria3_composite_strain, self.ctria6_composite_strain, self.ctriar_composite_strain,

            #self.nonlinearPlateStrain,
            #self.ctriax_strain, self.hyperelasticPlateStress,

            # solids
            self.ctetra_strain,
            self.cpenta_strain,
            self.chexa_strain,

            # rods
            #self.nonlinearRodStrain,  # non-vectorized

            self.chexa_strain,
            self.conrod_strain,
            self.cpenta_strain,
            self.cquad4_strain,
            self.cquad8_strain,
            self.cquadr_strain,
            self.crod_strain,
            self.cshear_strain,
            self.ctetra_strain,
            self.ctria3_strain,
            self.ctria6_strain,
            self.ctriar_strain,
            self.ctube_strain,

            # bush
            self.cbush_strain,
        ]

        oug = [
            self.accelerations,
            self.displacements, self.displacements_PSD, self.displacements_ATO, self.displacements_RMS,
            #self.scaled_displacements,  # ???
            self.temperatures,
            self.velocities, self.eigenvectors,
        ]
        oqg_mpc = [self.mpc_forces,]
        oqg_spc = [self.spc_forces,]
        ogs = [self.grid_point_stresses, self.grid_point_volume_stresses, ]
        ogp = [self.grid_point_forces,]
        other = [
            #self.forceVectors,
            #self.loadVectors,
            self.thermal_load_vectors,
        ]
        isubcases = sorted(self.isubcase_name_map.keys())
        #title = self.title

        res_categories = [
            ('OUG', oug),
            ('OQG_SPC', oqg_spc),
            ('OQG_MPC', oqg_mpc),
            ('OEF', oef),
            ('OES1X', oes1x1),
            ('OES1C', oes1c),
            ('OSTR', strain),
            ('OGS', ogs),
            ('OGP', ogp),
            ('other', other)
        ]
        #res_outs = {}

        # TODO: add a summary sheet

        # TODO: this may need to be reworked such that all of subcase 1
        #is printed before subcase 2
        for res_category_name, res_category in res_categories:
            #print("res_category_name = %s" % res_category_name)
            for res_type in res_category:
                res_keys = isubcases
                for res_key in res_keys:
                    isubcase = res_key
                    if isubcase in res_type:
                        #(subtitle, label) = self.isubcase_name_map[isubcase]
                        result = res_type[isubcase]
                        element_name = ''
                        if hasattr(result, 'element_name'):
                            element_name = ' - ' + result.element_name
                        if hasattr(result, '_write_xlsx'):
                            class_name = result.__class__.__name__
                            sheet_name = '%s_%s' % (class_name, isubcase)
                            assert len(sheet_name) < 31, sheet_name
                            Sheet(isheet).name = sheet_name
                            sheet = Sheet(isheet)

                            print(' %s - isubcase=%i%s' % (result.__class__.__name__, isubcase, element_name))
                            result._write_xlsx(sheet, is_mag_phase=False)
                            isheet += 1
                        else:
                            print("  *xlsx - %s not written" % result.__class__.__name__)

        #footer = [4, 0, 4]
        workbook.save()
