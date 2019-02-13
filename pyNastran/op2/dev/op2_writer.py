#pylint: disable=W0201,C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
#import copy
from datetime import date
#from collections import defaultdict
from struct import pack, Struct

import pyNastran
from pyNastran.op2.op2_interface.op2_f06_common import OP2_F06_Common
from pyNastran.op2.op2_interface.write_utils import _write_markers
from .writer.geom1 import write_geom1
from .writer.geom2 import write_geom2
from .writer.ept import write_ept
from .writer.mpt import write_mpt


def make_stamp(title, today=None):
    if 'Title' is None:
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


class OP2Writer(OP2_F06_Common):
    def __init__(self, op2):
        self.log = op2.log
        OP2_F06_Common.__init__(self)
        self.card_count = {}

    #def make_f06_header(self):
        #"""If this class is inherited, the F06 Header may be overwritten"""
        #return make_f06_header()

    def make_stamp(self, title, today):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return make_stamp(title, today)

    def write_op2(self, op2_outname, obj=None, is_mag_phase=False,
                  delete_objects=True, post=-1, endian=b'<'):
        """
        Writes an OP2 file based on the data we have stored in the object

        Parameters
        ----------
        op2_outname : str
            the name of the F06 file to write
        obj : OP2(); default=None -> self
            the OP2 object if you didn't inherit the class
        is_mag_phase : bool; default=False
            should complex data be written using Magnitude/Phase
            instead of Real/Imaginary (default=False; Real/Imag)
           Real objects don't use this parameter.
        delete_objects : bool; default=True
            should objects be deleted after they're written
            to reduce memory (default=True)
        """
        assert op2_outname != 'ctria3.op2'
        print('writing %s' % op2_outname)
        struct_3i = Struct(endian + '3i')

        if obj is None:
            obj = self
        if isinstance(op2_outname, str):
            fop2 = open(op2_outname, 'wb')
            fop2_ascii = open(op2_outname + '.txt', 'w')
            print('op2 out = %r' % op2_outname)
        else:
            assert isinstance(op2_outname, file), 'type(op2_outname)= %s' % op2_outname
            fop2 = op2_outname
            op2_outname = op2.name
            print('op2_outname =', op2_outname)

        #op2_ascii.write('writing [3, 7, 0] header\n')
        #if markers == [3,]:  # PARAM, POST, -1
            #self.op2_reader.read_markers([3])
            #data = self.read_block()

            #self.op2_reader.read_markers([7])
            #data = self.read_block()
            #data = self._read_record()
            #self.op2_reader.read_markers([-1, 0])
        #elif markers == [2,]:  # PARAM, POST, -2
        if post == -1:
        #_write_markers(op2, op2_ascii, [3, 0, 7])
            fop2.write(struct_3i.pack(*[4, 3, 4,]))
            tape_code = b'NASTRAN FORT TAPE ID CODE - '
            fop2.write(pack('7i 28s i', *[4, 1, 4,
                                          4, 7, 4,
                                          28, tape_code, 28]))

            nastran_version = b'NX8.5   ' if obj.is_nx else b'XXXXXXXX'
            fop2.write(pack(b'4i 8s i', *[4, 2, 4,
                                         #4, 2, 4,
                                         #4, 1, 4,
                                         #4, 8, 4,
                                         8, nastran_version, 8]))
            fop2.write(pack(b'6i', *[4, -1, 4,
                                     4, 0, 4,]))
        elif post == -2:
            _write_markers(fop2, fop2_ascii, [2, 4])
        else:
            raise RuntimeError('post = %r; use -1 or -2' % post)

        write_geom1(fop2, fop2_ascii, obj)
        write_geom2(fop2, fop2_ascii, obj)
        write_ept(fop2, fop2_ascii, obj)
        write_mpt(fop2, fop2_ascii, obj)
        if obj.grid_point_weight.reference_point is not None:
            if hasattr(result, 'write_op2'):
                print("grid_point_weight")
                obj.grid_point_weight.write_op2(fop2, endian=endian)
            else:
                print("*op2 - grid_point_weight not written")


        #is_mag_phase = False

        # eigenvalues are written first
        for ikey, result in sorted(obj.eigenvalues.items()):
            header
            #print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
            if hasattr(result, 'write_op2'):
                result.write_op2(fop2, fop2_ascii, endian=endian)
                #if delete_objects:
                    #del result
            else:
                print("*op2 - %s not written" % result.__class__.__name__)
                write_op2

        # then eigenvectors
        # has a special header
        for isubcase, result in sorted(obj.eigenvectors.items()):
            (subtitle, label) = obj.isubcase_name_map[isubcase]

            if hasattr(result, 'write_op2'):
                print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
                result.write_op2(fop2, fop2_ascii, is_mag_phase=is_mag_phase, endian=endian)
                #if delete_objects:
                    #del result
            else:
                print("*op2 - %s not written" % result.__class__.__name__)
                write_op2

        # finally, we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order


        #if markers == [3,]:  # PARAM, POST, -2
            #self.op2_reader.read_markers([3])
            #data = self.read_block()
            #self.op2_reader.read_markers([7])
            #data = self.read_block()
            ##self.show(100)
            #data = self._read_record()


        oef = [
            # OEF - forces
            # alphabetical order...
            obj.cbar_force,

            # beam
            #obj.cbeam_forces,
            #obj.cbar100_forces,
            #obj.cbend_forces,
            obj.cbeam_force,
            obj.cbush_force,

            obj.celas1_force,
            obj.celas2_force,
            obj.celas3_force,
            obj.celas4_force,

            obj.cdamp1_force,
            obj.cdamp2_force,
            obj.cdamp3_force,
            obj.cdamp4_force,

            #obj.plateForces,   # centroidal elements
            #obj.plateForces2,  # bilinear elements

            obj.conrod_force,
            obj.cquad4_force,
            obj.cquad8_force,
            obj.crod_force,
            obj.cshear_force,
            obj.ctria3_force,
            obj.ctria6_force,
            obj.ctube_force,

            # other
            #obj.cgap_force, #obj.solidPressureForces,
        ]

        oes1x1 = [
            # cbars/cbeams
            obj.cbar_stress, obj.cbeam_stress,

            # bush
            obj.cbush_stress, obj.cbush1d_stress_strain,

            # rods
            #obj.nonlinearRodStress,

            obj.celas1_stress, obj.celas2_stress, obj.celas3_stress, obj.celas4_stress,

            obj.chexa_stress,
            obj.conrod_stress,
            obj.cpenta_stress,
            obj.cquad4_stress, obj.cquad8_stress, obj.cquadr_stress,
            obj.crod_stress,
            obj.cshear_stress,
            obj.ctetra_stress,
            obj.ctria3_stress, obj.ctria6_stress, obj.ctriar_stress,
            obj.ctube_stress,
        ]
        oes1c = [
            obj.cquad4_composite_stress, obj.cquad8_composite_stress, obj.cquadr_composite_stress,
            obj.ctria3_composite_stress, obj.ctria6_composite_stress, obj.ctriar_composite_stress,

            #obj.nonlinearPlateStress,
            obj.ctriax_stress, #obj.hyperelasticPlateStrain,
        ]
        #stress = oes1x1 + oes1c

        strain = [
            # bars/beams
            obj.cbar_strain, obj.cbeam_strain,

            # springs,
            obj.celas1_strain, obj.celas2_strain, obj.celas3_strain, obj.celas4_strain,

            # plates
            obj.ctria3_strain, obj.cquad4_strain,
            obj.cshear_strain,

            #obj.nonlinearPlateStrain,
            #obj.ctriax_strain, obj.hyperelasticPlateStress,

            # solids
            obj.ctetra_strain,
            obj.cpenta_strain,
            obj.chexa_strain,

            # rods
            #obj.nonlinearRodStrain,  # non-vectorized

            obj.chexa_strain,
            obj.conrod_strain,
            obj.cpenta_strain,
            obj.cquad4_strain,
            obj.cquad8_strain,
            obj.cquadr_strain,
            obj.crod_strain,
            obj.cshear_strain,
            obj.ctetra_strain,
            obj.ctria3_strain,
            obj.ctria6_strain,
            obj.ctriar_strain,
            obj.ctube_strain,

            # bush
            obj.cbush_strain,
        ]
        ostr1c = [
            obj.cquad4_composite_strain, obj.cquad8_composite_strain, obj.cquadr_composite_strain,
            obj.ctria3_composite_strain, obj.ctria6_composite_strain, obj.ctriar_composite_strain,
        ]

        oug = [
            obj.accelerations,
            obj.displacements, #obj.displacementsPSD, obj.displacementsATO, obj.displacementsRMS,
            #obj.scaledDisplacements,  # ???
            obj.temperatures,
            obj.velocities, obj.eigenvectors,
        ]
        oqg_mpc = [obj.mpc_forces]
        oqg_spc = [obj.spc_forces]
        ogs = [
            #obj.grid_point_stresses,
            #obj.grid_point_volume_stresses,
        ]
        opg = [obj.load_vectors]
        ogp = [obj.grid_point_forces]
        other = [
            #obj.forceVectors,
            #obj.loadVectors,
            obj.thermal_load_vectors,
        ]
        isubcases = sorted(obj.isubcase_name_map.keys())
        #title = obj.title

        res_categories = [
            ('OUG', oug),
            ('OQG_SPC', oqg_spc),
            ('OQG_MPC', oqg_mpc),
            ('OPG1', opg),
            ('OEF', oef),
            ('OES1X', oes1x1),
            ('OES1C', oes1c),
            ('OSTR', strain),
            ('OSTR1C', ostr1c),
            ('OGS', ogs),
            ('OGP', ogp),
            ('other', other)
        ]
        res_outs = {}
        self._write_categories(obj, res_categories, isubcases,
                               fop2, fop2_ascii,
                               struct_3i, endian)

    def _write_categories(self, obj, res_categories, isubcases,
                          fop2, fop2_ascii,
                          struct_3i, endian=b'<'):
        # TODO: this may need to be reworked such that all of subcase 1
        #is printed before subcase 2
        for res_category_name, res_category in res_categories:
            case_count = self._write_category(
                obj, res_category_name, res_category, isubcases,
                fop2, fop2_ascii,
                struct_3i, endian=b'<')

            if case_count:
                print('res_category_name=%s case_count=%s'  % (res_category_name, case_count))
                # close off the result
                footer = [4, 0, 4]
                fop2.write(struct_3i.pack(*footer))
                fop2_ascii.write('close_a = %s\n' % footer)

        # close off the op2
        footer = [4, 0, 4]
        fop2.write(struct_3i.pack(*footer))
        fop2_ascii.write('close_b = %s\n' % footer)
        fop2.close()
        fop2_ascii.close()

    def _write_category(self, obj, res_category_name, res_category, isubcases,
                        fop2, fop2_ascii,
                        struct_3i, endian=b'<'):
        case_count = 0
        itable = -1
        print_msg = True
        for ires_type, res_type in enumerate(res_category):
            res_keys = isubcases
            for res_key in res_keys:
                isubcase = res_key
                if isubcase in res_type:
                    if print_msg:
                        print("res_category_name = %s" % res_category_name)
                        print_msg = False
                    #(subtitle, label) = obj.isubcase_name_map[isubcase]
                    result = res_type[isubcase]
                    element_name = ''
                    if hasattr(result, 'element_name'):
                        element_name = ' - ' + result.element_name

                    if hasattr(result, 'write_op2'):
                        #if hasattr(result, 'is_bilinear') and result.is_bilinear():
                            #obj.log.warning("  *op2 - %s (%s) not written" % (
                                #result.__class__.__name__, result.element_name))
                            #continue
                        #print(' %s - isubcase=%i%s' % (result.__class__.__name__, isubcase, element_name))
                        itable = result.write_op2(fop2, fop2_ascii, itable, obj.date, is_mag_phase=False, endian=endian)
                    else:
                        obj.log.warning("  *op2 - %s not written" % result.__class__.__name__)
                        continue

                    case_count += 1
                    header = [
                        4, itable, 4,
                        4, 1, 4,
                        4, 0, 4,
                    ]
                    #print('writing itable=%s' % itable)
                    assert itable is not None, '%s itable is None' % result.__class__.__name__
                    fop2.write(pack(b'9i', *header))
                    fop2_ascii.write('footer2 = %s\n' % header)
                    #print('bailing...')
                    #return case_count
        return case_count
