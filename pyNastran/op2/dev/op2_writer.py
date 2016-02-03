#pylint: disable=W0201,C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
import sys
import copy
from datetime import date
from collections import defaultdict
from struct import pack

import pyNastran
from pyNastran.op2.op2_f06_common import OP2_F06_Common
from pyNastran.op2.write_utils import _write_markers

def make_stamp(Title, today=None):
    if 'Title' is None:
        Title = ''

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
    if Title is None:
        Title = ''
    out = '1    %-67s   %-19s %-22s PAGE %%5i\n' % (Title.strip(), str_today, build)
    return out


class OP2Writer(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)
        self.card_count = {}

    def make_f06_header(self):
        """If this class is inherited, the F06 Header may be overwritten"""
        return make_f06_header()

    def make_stamp(self, Title, today):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return make_stamp(Title, today)

    def write_op2(self, op2_outname, is_mag_phase=False,
                  delete_objects=True):
        """
        Writes an OP2 file based on the data we have stored in the object

        :param op2_outname:  the name of the F06 file to write
        :param is_mag_phase: should complex data be written using Magnitude/Phase
                         instead of Real/Imaginary (default=False; Real/Imag)
                         Real objects don't use this parameter.
        :param delete_objects: should objects be deleted after they're written
                         to reduce memory (default=True)
        """
        assert op2_outname != 'ctria3.op2'
        print('writing %s' % op2_outname)
        if isinstance(op2_outname, str):
            op2 = open(op2_outname, 'wb')
            op2_ascii = open(op2_outname + '.txt', 'wb')
        else:
            assert isinstance(op2_outname, file), 'type(op2_outname)= %s' % op2_outname
            op2 = op2_outname
            op2_outname = op2.name
            print('op2_outname =', op2_outname)

        #op2_ascii.write('writing [3, 7, 0] header\n')
        #if markers == [3,]:  # PARAM, POST, -1
            #self.read_markers([3])
            #data = self.read_block()

            #self.read_markers([7])
            #data = self.read_block()
            #data = self._read_record()
            #self.read_markers([-1, 0])
        #elif markers == [2,]:  # PARAM, POST, -2
        if 1:
        #_write_markers(op2, op2_ascii, [3, 0, 7])
            op2.write(pack('3i', *[4, 3, 4,]))
            tape_code = b'NASTRAN FORT TAPE ID CODE - '
            op2.write(pack('7i 28s i', *[4, 1, 4,
                                         4, 7, 4,
                                         28, tape_code, 28]))

            nastran_version = b'NX8.5   ' if self.is_nx else b'XXXXXXXX'
            op2.write(pack(b'4i 8s i', *[4, 2, 4,
                                         #4, 2, 4,
                                         #4, 1, 4,
                                         #4, 8, 4,
                                         8, nastran_version, 8]))
            op2.write(pack(b'6i', *[4, -1, 4,
                                    4, 0, 4,]))
        else:
            _write_markers(op2, op2_ascii, [2, 4])

        if self.grid_point_weight.reference_point is not None:
            if hasattr(result, 'write_op2'):
                print("grid_point_weight")
                self.grid_point_weight.write_op2(op2, page_stamp, self.page_num)
            else:
                print("*op2 - grid_point_weight not written")


        #is_mag_phase = False

        # eigenvalues are written first
        for ikey, result in sorted(iteritems(self.eigenvalues)):
            header
            #print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
            if hasattr(result, 'write_op2'):
                result.write_op2(op2, op2_ascii)
                if delete_objects:
                    del result
            else:
                print("*op2 - %s not written" % result.__class__.__name__)
                write_op2

        # then eigenvectors
        # has a special header
        for isubcase, result in sorted(iteritems(self.eigenvectors)):
            (subtitle, label) = self.iSubcaseNameMap[isubcase]

            if hasattr(result, 'write_op2'):
                print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
                result.write_op2(op2, op2_ascii, is_mag_phase=is_mag_phase)
                if delete_objects:
                    del result
            else:
                print("*op2 - %s not written" % result.__class__.__name__)
                write_op2

        # finally, we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order


        #if markers == [3,]:  # PARAM, POST, -2
            #self.read_markers([3])
            #data = self.read_block()
            #self.read_markers([7])
            #data = self.read_block()
            ##self.show(100)
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
            self.displacements, self.displacementsPSD, self.displacementsATO, self.displacementsRMS,
            self.scaledDisplacements,  # ???
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
        isubcases = sorted(self.iSubcaseNameMap.keys())
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
        res_outs = {}

        # TODO: this may need to be reworked such that all of subcase 1
        #is printed before subcase 2
        for res_name, res_category in res_categories:
            print("res_name = %s" % res_name)
            for res_type in res_category:
                res_keys = isubcases
                itable = -1
                for res_key in res_keys:
                    isubcase = res_key
                    if isubcase in res_type:
                        #(subtitle, label) = self.iSubcaseNameMap[isubcase]
                        result = res_type[isubcase]
                        element_name = ''
                        if hasattr(result, 'element_name'):
                            element_name = ' - ' + result.element_name
                        if hasattr(result, 'write_op2'):
                            print(' %s - isubcase=%i%s' % (result.__class__.__name__, isubcase, element_name))
                            result.write_op2(op2, op2_ascii, itable, self.date, is_mag_phase=False)
                        else:
                            print("  *op2 - %s not written" % result.__class__.__name__)

            footer = [4, 0, 4]
            op2.write(pack(b'3i', *footer))
        #footer = [4, 0, 4]
        op2.write(pack(b'3i', *footer))
        op2.close()
