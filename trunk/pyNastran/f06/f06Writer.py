#pylint: disable=W0201,C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import copy
from datetime import date
import warnings

import pyNastran
from pyNastran.f06.tables.grid_point_weight import GridPointWeight


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
        #print "day=%s month=%s year=%s" % (day, month, year)
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


def make_f06_header():
    spaces = ''
    lines1 = [
        spaces + '/* -------------------------------------------------------------------  */\n',
        spaces + '/*                              PYNASTRAN                               */\n',
        spaces + '/*                      - NASTRAN FILE INTERFACE -                      */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*              A Python reader/editor/writer for the various           */\n',
        spaces + '/*                        NASTRAN file formats.                         */\n',
        spaces + '/*                       Copyright (C) 2011-2013                        */\n',
        spaces + '/*               Steven Doyle, Al Danial, Marcin Garrozik               */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*    This program is free software; you can redistribute it and/or     */\n',
        spaces + '/*    modify it under the terms of the GNU Lesser General Public        */\n',
        spaces + '/*    License as published by the Free Software Foundation;             */\n',
        spaces + '/*    version 3 of the License.                                         */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*    This program is distributed in the hope that it will be useful,   */\n',
        spaces + '/*    but WITHOUT ANY WARRANTY; without even the implied warranty of    */\n',
        spaces + '/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */\n',
        spaces + '/*    GNU Lesser General Public License for more details.               */\n',
        spaces + '/*                                                                      */\n',
        spaces + '/*    You should have received a copy of the GNU Lesser General Public  */\n',
        spaces + '/*    License along with this program; if not, write to the             */\n',
        spaces + '/*    Free Software Foundation, Inc.,                                   */\n',
        spaces + '/*    675 Mass Ave, Cambridge, MA 02139, USA.                           */\n',
        spaces + '/* -------------------------------------------------------------------  */\n',
        '\n']

    spaces = 46 * ' '
    version = 'Version %8s' % pyNastran.__version__
    lines2 = [
        spaces + '* * * * * * * * * * * * * * * * * * * *\n',
        spaces + '* * * * * * * * * * * * * * * * * * * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *            pyNastran            * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *%s* *\n' % version.center(33),
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *          %15s        * *\n' % pyNastran.__releaseDate2__,
        spaces + '* *                                 * *\n',
        spaces + '* *            Questions            * *\n',
        spaces + '* *        mesheb82@gmail.com       * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* *                                 * *\n',
        spaces + '* * * * * * * * * * * * * * * * * * * *\n',
        spaces + '* * * * * * * * * * * * * * * * * * * *\n\n\n']
    return ''.join(lines1 + lines2)


def sorted_bulk_data_header():
    msg = '0                                                 S O R T E D   B U L K   D A T A   E C H O                                         \n'
    msg += '                 ENTRY                                                                                                              \n'
    msg += '                 COUNT        .   1  ..   2  ..   3  ..   4  ..   5  ..   6  ..   7  ..   8  ..   9  ..  10  .                      \n'
    return msg

def make_end(end_flag=False):
    lines = []
    lines2 = []
    if end_flag:
        lines = ['', '',
        '0                                   * * * *  A N A L Y S I S  S U M M A R Y  T A B L E  * * * *',
        '0 SEID  PEID PROJ VERS APRCH      SEMG SEMR SEKR SELG SELR MODES DYNRED SOLLIN PVALID SOLNL LOOPID DESIGN CYCLE SENSITIVITY',
        ' --------------------------------------------------------------------------------------------------------------------------']
    #     0     0    1    1 '        '    T    T    T    T    T     F      F      T      0     F     -1            0           F

        seid = 0
        peid = 0
        proj = 1
        vers = 1
        approach = '        '

        SELG = 'T'
        SEMG = 'T'
        SEMR = 'T'
        SEKR = 'T'
        SELR = 'T'
        MODES = 'F'
        DYNRED = 'F'

        SOLLIN = 'T'
        PVALID = 0
        SOLNL = 'F'
        LOOPID = -1
        CYCLE = 0
        SENSITIVITY = 'F'

        msg = '     %s     %s    %s    %s %8r    %s    %s    %s    %s    %s     %s      %s      %s      %s     %s     %s            %s           %s' % (
            seid, peid, proj, vers, approach, SEMG, SEMR, SEKR, SELG, SELR, MODES, DYNRED, SOLLIN, PVALID, SOLNL,
            LOOPID, CYCLE, SENSITIVITY)
        lines.append(msg)

        lines2 = ['0SEID = SUPERELEMENT ID.',
        ' PEID = PRIMARY SUPERELEMENT ID OF IMAGE SUPERELEMENT.',
        ' PROJ = PROJECT ID NUMBER.',
        ' VERS = VERSION ID.',
        ' APRCH = BLANK FOR STRUCTURAL ANALYSIS.  HEAT FOR HEAT TRANSFER ANALYSIS.',
        ' SEMG = STIFFNESS AND MASS MATRIX GENERATION STEP.',
        ' SEMR = MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES).',
        ' SEKR = STIFFNESS MATRIX REDUCTION STEP.',
        ' SELG = LOAD MATRIX GENERATION STEP.',
        ' SELR = LOAD MATRIX REDUCTION STEP. ',
        ' MODES = T (TRUE) IF NORMAL MODES OR BUCKLING MODES CALCULATED.',
        ' DYNRED = T (TRUE) MEANS GENERALIZED DYNAMIC AND/OR COMPONENT MODE REDUCTION PERFORMED.',
        ' SOLLIN = T (TRUE) IF LINEAR SOLUTION EXISTS IN DATABASE.',
        ' PVALID = P-DISTRIBUTION ID OF P-VALUE FOR P-ELEMENTS',
        ' LOOPID = THE LAST LOOPID VALUE USED IN THE NONLINEAR ANALYSIS.  USEFUL FOR RESTARTS.',
        ' SOLNL = T (TRUE) IF NONLINEAR SOLUTION EXISTS IN DATABASE.',
        ' DESIGN CYCLE = THE LAST DESIGN CYCLE (ONLY VALID IN OPTIMIZATION).',
        ' SENSITIVITY = SENSITIVITY MATRIX GENERATION FLAG.',
        ' ',
        ' No PARAM values were set in the Control File.']

    lines3 = [' ',
             '1                                        * * * END OF JOB * * *',
             ' ',
             ' ']
    return '\n'.join(lines+lines2+lines3)


#class F06WriterDeprecated(object):
#    def __init__(self):
#        pass

class F06Writer(object):
#class F06Writer(F06WriterDeprecated):
    def __init__(self):
        #F06WriterDeprecated.__init__(self)
        self.card_count = {}

        #: BDF Title
        self.Title = None

        #: a dictionary that maps an integer of the subcaseName to the
        #: subcaseID
        self.iSubcaseNameMap = {}

        self.page_num = 1

        self.iSubcases = []
        self.__objects_vector_init__()
        self.__objects_init__()


    def __objects_vector_init__(self):
        """
        All OUG table is simple to vectorize, so we declere it in __objects_init__
        On the other hand, the rodForces object contains CROD/CTUBE/CONROD elements.
        It is difficult to handle initializing the CRODs/CONRODs given a
        mixed type case, so we split out the elements.
        """
        #======================================================================
        # rods
        self.crod_force = {}
        self.conrod_force = {}
        self.ctube_force = {}

        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}

        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}

        #======================================================================
        # springs
        self.celas1_force = {}
        self.celas2_force = {}
        self.celas3_force = {}
        self.celas4_force = {}

        self.celas1_stress = {}
        self.celas2_stress = {}
        self.celas3_stress = {}
        self.celas4_stress = {}

        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        #======================================================================
        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}
        #======================================================================

    def __objects_init__(self):
        """More variable declarations"""
        #: the date the job was run on
        self.date = None

        #: Grid Point Weight Table
        #: create with:
        #:   PARAM   GRDPNT    0  (required for F06/OP2)
        #:   PARAM   POSTEXT YES  (required for OP2)
        self.grid_point_weight = GridPointWeight()

        #: ESE
        self.eigenvalues = {}

        #: OUG - displacement
        self.displacements = {}           # tCode=1 thermal=0
        self.displacementsPSD = {}        # random
        self.displacementsATO = {}        # random
        self.displacementsRMS = {}        # random
        self.displacementsCRM = {}        # random
        self.displacementsNO = {}        # random
        self.scaledDisplacements = {}     # tCode=1 thermal=8

        #: OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        #: OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0

        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        #: OUG - acceleration
        self.accelerations = {}           # tCode=11 thermal=0

        # OEF - Forces - tCode=4 thermal=0
        # rods
        self.rodForces = {}

        self.barForces = {}
        self.bar100Forces = {}
        self.beamForces = {}
        self.bendForces = {}
        self.bushForces = {}
        self.coneAxForces = {}
        self.damperForces = {}
        self.gapForces = {}

        self.plateForces = {}
        self.ctria3_force = {}
        self.cquad4_force = {}

        self.plateForces2 = {}
        self.compositePlateForces = {}

        self.shearForces = {}
        self.cshear_force = {}

        self.solidPressureForces = {}
        self.springForces = {}
        self.viscForces = {}

        self.force_VU = {}
        self.force_VU_2D = {}

        #OEF - Fluxes - tCode=4 thermal=1
        self.thermalLoad_CONV = {}
        self.thermalLoad_CHBDY = {}
        self.thermalLoad_1D = {}
        self.thermalLoad_2D_3D = {}
        self.thermalLoad_VU = {}
        self.thermalLoad_VU_3D = {}
        self.thermalLoad_VUBeam = {}
        #self.temperatureForces = {}

        # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
        #: OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
        self.celasStress = {}

        #: OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
        self.celasStrain = {}

        #: OES - CTRIAX6
        self.ctriaxStress = {}
        self.ctriaxStrain = {}

        #: OES - isotropic CROD/CONROD/CTUBE stress
        self.rodStress = {}

        #: OES - isotropic CROD/CONROD/CTUBE strain
        self.rodStrain = {}

        #: OES - nonlinear CROD/CONROD/CTUBE stress
        self.nonlinearRodStress = {}
        #: OES - nonlinear CROD/CONROD/CTUBE strain
        self.nonlinearRodStrain = {}
        #: OES - isotropic CBAR stress
        self.barStress = {}
        #: OES - isotropic CBAR strain
        self.barStrain = {}
        #: OES - isotropic CBEAM stress
        self.beamStress = {}
        #: OES - isotropic CBEAM strain
        self.beamStrain = {}
        #: OES - isotropic CBUSH stress
        self.bushStress = {}
        #: OES - isotropic CBUSH strain
        self.bushStrain = {}
         #: OES - isotropic CBUSH1D strain/strain
        self.bush1dStressStrain = {}

        #: OES - isotropic CTRIA3/CQUAD4 stress
        self.plateStress = {}
        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}
        #: OES - isotropic CTRIA3/CQUAD4 strain
        self.plateStrain = {}
        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

        #: OESNLXR - CTRIA3/CQUAD4 stress
        self.nonlinearPlateStress = {}
        #: OESNLXR - CTRIA3/CQUAD4 strain
        self.nonlinearPlateStrain = {}
        self.hyperelasticPlateStress = {}
        self.hyperelasticPlateStrain = {}

        #: OES - isotropic CTETRA/CHEXA/CPENTA stress
        self.solidStress = {}

        #: OES - isotropic CTETRA/CHEXA/CPENTA strain
        self.solidStrain = {}

        #: OES - composite CTRIA3/CQUAD4 stress
        self.compositePlateStress = {}
        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        #: OES - composite CTRIA3/CQUAD4 strain
        self.compositePlateStrain = {}
        self.cquad4_composite_strain = {}
        self.cquad8_composite_strain = {}
        self.ctria3_composite_strain = {}
        self.ctria6_composite_strain = {}


        #: OES - CSHEAR stress
        self.shearStress = {}
        self.cshear_stress = {}
        #: OES - CSHEAR strain
        self.shearStrain = {}
        self.cshear_strain = {}

        #: OES - CELAS1 224, CELAS3 225,
        self.nonlinearSpringStress = {}
        #: OES - GAPNL 86
        self.nonlinearGapStress = {}
        #: OES - CBUSH 226
        self.nolinearBushStress = {}

        # OQG - spc/mpc forces
        self.spcForces = {}  # tCode=3?
        self.mpcForces = {}  # tCode=39

        # OQG - thermal forces
        self.thermalGradientAndFlux = {}

        #: OGF - grid point forces
        self.gridPointForces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.gridPointStresses = {}       # tCode=26
        self.gridPointVolumeStresses = {}  # tCode=27

        #: OPG - summation of loads for each element
        self.loadVectors = {}       # tCode=2  thermal=0
        self.thermalLoadVectors = {}  # tCode=2  thermal=1
        self.appliedLoads = {}       # tCode=19 thermal=0
        self.forceVectors = {}       # tCode=12 thermal=0

        #: OEE - strain energy density
        self.strainEnergy = {}  # tCode=18

    def make_f06_header(self):
        """If this class is inherited, the F06 Header may be overwritten"""
        return make_f06_header()

    def make_stamp(self, Title, today):
        """If this class is inherited, the PAGE stamp may be overwritten"""
        return make_stamp(Title, today)

    def make_grid_point_singularity_table(self, failed):
        msg = ''
        if failed:
            msg += '0                                         G R I D   P O I N T   S I N G U L A R I T Y   T A B L E\n'
            msg += '0                             POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET\n'
            msg += '                               ID            DIRECTION      RATIO     EXCLUSIVE  UNION   EXCLUSIVE  UNION\n'
            for (nid, dof) in failed:
                msg += '                                %s        G      %s         0.00E+00          B        F         SB       SB   *\n' % (nid, dof)
        else:
            msg += 'No constraints have been applied...\n'

        page_stamp = self.make_stamp(self.Title, self.date)
        msg += page_stamp % self.page_num
        self.page_num += 1
        return msg

    def write_oload(self, model, Fg):
        msg = ''
        msg += '        *** USER INFORMATION MESSAGE 7310 (VECPRN)\n'
        msg += '            ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM WILL BE USED AS REFERENCE LOCATION.\n'
        msg += '            RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN SUPERELEMENT BASIC SYSTEM COORDINATES.\n'
        msg += '       0                                                  OLOAD    RESULTANT       \n'

        #nnodes = len(model.nodes)

        msg += '        SUBCASE/    LOAD\n'
        msg += '        DAREA ID    TYPE       T1            T2            T3            R1            R2            R3\n'
        msg += '      0        1     FX    3.000000E+03     ----          ----          ----       0.000000E+00 -6.000000E+04                             \n'
        msg += '                     FY       ----       5.000000E+03     ----       0.000000E+00     ----       2.000000E+05                             \n'
        msg += '                     FZ       ----          ----       0.000000E+00  0.000000E+00  0.000000E+00     ----                                  \n'
        msg += '                     MX       ----          ----          ----       1.300000E+04     ----          ----                                  \n'
        msg += '                     MY       ----          ----          ----          ----       0.000000E+00     ----                                  \n'
        msg += '                     MZ       ----          ----          ----          ----          ----       0.000000E+00                             \n'
        msg += '                   TOTALS  3.000000E+03  5.000000E+03  0.000000E+00  1.300000E+04  0.000000E+00  1.400000E+05\n'

        page_stamp = self.make_stamp(self.Title, self.date)
        msg += page_stamp % self.page_num
        self.page_num += 1
        return msg

    def write_summary(self, f06, card_count=None):
        summary_header = '                                        M O D E L   S U M M A R Y\n\n'
        summary = ''

        self.cardsToRead = set([

            # rigid elements
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3',

            # spc/mpc constraints
            'SPC', 'SPCADD', 'SPC1', 'SPCD', 'SPCAX',
            'MPC', 'MPCADD',
            'SUPORT', 'SUPORT1',

            # aero cards
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',

            # temperature cards
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'CONV',
        ])


        blocks = [
            ['POINTS', ['GRID', 'GRDSET', ]],
            ['ENTRIES', ['SPOINT']],

            ['ELEMENTS', [# these are sorted
                        # elements
                        'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

                        # springs
                        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CELAS5',

                        # bushings
                        'CBUSH', 'CBUSH1D', 'CBUSH2D',

                        # dampers
                        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',

                        # bar flags
                        'BAROR', 'CBARAO',
                        # bars
                        'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',

                        # shells
                        'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
                        'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',

                        # solids
                        'CTETRA', 'CPENTA', 'CHEXA',

                        # other
                        'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
                        'CGAP', 'CFAST',

                        # thermal
                        'CHBDYP', 'CHBDYG', 'CONV',
                          ]],
            ['ELEMENTS', ['RBE2', 'RBE3']],
        ]
        #print("self.card_count", self.card_count)
        if card_count is None:
            card_count = self.card_count

        for block in blocks:
            block_name, keys = block
            key_count = 0
            for key in sorted(keys):
                try:
                    value = card_count[key]
                    summary += '                                   NUMBER OF %-8s %-8s = %8s\n' % (key, block_name, value)
                    key_count += 1
                except KeyError:
                    pass
            if key_count:
                summary += ' \n'
        #sys.exit(summary)
        if summary:
            f06.write(summary_header)
            f06.write(summary)

            page_stamp = self.make_stamp(self.Title, self.date)
            f06.write(page_stamp % self.page_num)
            self.page_num += 1
            #print(summary)

    def write_f06(self, f06_outname, is_mag_phase=False,
                  delete_objects=True, end_flag=False):
        """
        Writes an F06 file based on the data we have stored in the object

        :param self:         the F06 object
        :param f06_outname:  the name of the F06 file to write
        :param is_mag_phase: should complex data be written using Magnitude/Phase
                         instead of Real/Imaginary (default=False; Real/Imag)
                         Real objects don't use this parameter.
        :param delete_objects: should objects be deleted after they're written
                         to reduce memory (default=True)
        :param end_flag: should a dummy Nastran "END" table be made
                         (default=False)
        """
        self._write_f06(f06_outname, is_mag_phase=is_mag_phase, delete_objects=delete_objects,
                       end_flag=end_flag, write_op2=False)

    def _write_f06(self, f06_outname, is_mag_phase=False,
                  delete_objects=True, end_flag=False, write_op2=False):
        if isinstance(f06_outname, str):
            f06 = open(f06_outname, 'wb')
            self.write_summary(f06)
        else:
            assert isinstance(f06_outname, file), 'type(f06_outname)= %s' % f06_outname
            f06 = f06_outname
            f06_outname = f06.name
            print('f06_outname =', f06_outname)
        if write_op2:
            op2 = open(f06_outname.replace('.f06', '.op2'), 'wb')

        page_stamp = self.make_stamp(self.Title, self.date)
        if self.grid_point_weight.reference_point is not None:
            print("grid_point_weight")
            self.page_num = self.grid_point_weight.write_f06(f06, page_stamp, self.page_num)
            assert isinstance(self.page_num, int), self.grid_point_weight.__class__.__name__

        #print "page_stamp = %r" % page_stamp
        #print "stamp      = %r" % stamp

        #is_mag_phase = False
        header = ['     DEFAULT                                                                                                                        \n',
                  '\n', '']

        # eigenvalues are written first
        for ikey, result in sorted(self.eigenvalues.iteritems()):
            header
            #print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
            self.page_num = result.write_f06(f06, header, page_stamp,
                                            page_num=self.page_num)
            assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
            if delete_objects:
                del result
            self.page_num += 1

        # then eigenvectors
        # has a special header
        for isubcase, result in sorted(self.eigenvectors.iteritems()):
            (subtitle, label) = self.iSubcaseNameMap[isubcase]
            subtitle = subtitle.strip()
            header[0] = '     %s\n' % subtitle
            header[1] = '0                                                                                                            SUBCASE %i\n' % (isubcase)
            #header[2] = complex/nonlinear
            print('%-18s SUBCASE=%i' % (result.__class__.__name__, isubcase))
            self.page_num = result.write_f06(header, page_stamp,
                                            self.page_num, f=f06, is_mag_phase=is_mag_phase)
            assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
            if delete_objects:
                del result
            self.page_num += 1

        # finally, we writte all the other tables
        # nastran puts the tables in order of the Case Control deck,
        # but we're lazy so we just hardcode the order

        # subcase name, subcase ID, transient word & value
        header_old = ['     DEFAULT                                                                                                                        \n',
                     '\n', ' \n']
        header = copy.deepcopy(header_old)
        res_types = [
            self.accelerations,
            self.displacements, self.displacementsPSD, self.displacementsATO, self.displacementsRMS,
            self.scaledDisplacements,  # ???

            self.forceVectors,
            self.loadVectors,
            self.temperatures,
            self.velocities, #self.eigenvectors,

            self.mpcForces,
            self.spcForces,
            self.thermalLoadVectors,


            #------------------------------------------
            # OEF - forces

            # alphabetical order...
            # bars
            self.barForces,

            # beam
            self.beamForces,
            self.bar100Forces,
            self.bendForces,

            # alphabetical
            self.conrod_force,
            self.cquad4_force,
            self.plateForces,   # centroidal elements
            self.plateForces2,  # bilinear elements

            self.crod_force,
            self.cshear_force,
            self.ctria3_force,
            self.ctube_force,

            # rods
            self.rodForces,

            # springs
            self.springForces,

            # dampers
            self.damperForces,

            # cshear,
            self.shearForces,
            # other
            self.bushForces, self.gapForces, self.solidPressureForces,

            #------------------------------------------
            # OES - strain
            # 1.  cbar
            # 2.  cbeam
            # 3.  crod/ctube/conrod

            # springs,
            self.celasStrain,

            # bars/beams
            self.barStrain, self.beamStrain,

            # plates
            self.plateStrain,
            self.shearStrain,
            self.compositePlateStrain,
            self.cquad4_composite_strain,
            self.cquad8_composite_strain,
            self.ctria3_composite_strain,
            self.ctria6_composite_strain,

            self.nonlinearPlateStrain,
            self.ctriaxStrain, self.hyperelasticPlateStress,


            # solids
            self.solidStrain,

            # rods
            self.rodStrain, self.nonlinearRodStrain,  # non-vectorized


            self.celas1_strain,
            self.celas2_strain,
            self.celas3_strain,
            self.celas4_strain,

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
            self.bushStrain,
            #------------------------------------------
            # cbars/cbeams
            self.barStress,
            self.beamStress,

            # bush
            self.bushStress, self.bush1dStressStrain,

            self.celasStress,
            self.shearStress,
            self.plateStress,
            self.solidStress,

            # rods
            self.rodStress, self.nonlinearRodStress,


            # shear
            # OES - stress
            self.celas1_stress,
            self.celas2_stress,
            self.celas3_stress,
            self.celas4_stress,

            self.chexa_stress,
            self.conrod_stress,
            self.cpenta_stress,
            self.cquad4_stress,
            self.cquad8_stress,
            self.cquadr_stress,
            self.crod_stress,
            self.cshear_stress,
            self.ctetra_stress,
            self.ctria3_stress,
            self.ctria6_stress,
            self.ctriar_stress,
            self.ctube_stress,

            self.compositePlateStress,
            self.cquad4_composite_stress,
            self.cquad8_composite_stress,
            self.ctria3_composite_stress,
            self.ctria6_composite_stress,

            self.nonlinearPlateStress,
            self.ctriaxStress, self.hyperelasticPlateStrain,

            #------------------------------------------

            self.gridPointStresses, self.gridPointVolumeStresses, self.gridPointForces,
        ]

        if 1:
            iSubcases = sorted(self.iSubcaseNameMap.keys())
            #print("self.iSubcaseNameMap = %s" %(self.iSubcaseNameMap))
            for isubcase in iSubcases:
                title = self.Title
                (subtitle, label) = self.iSubcaseNameMap[isubcase]
                subtitle = subtitle.strip()
                label = label.strip()
                #print "label = ",label

                (subtitle, label) = self.iSubcaseNameMap[isubcase]
                label = label.strip()
                subtitle = subtitle.strip()

                #header[0] = '     %-127s\n' % subtitle
                #header[1] = '0    %-72s                                SUBCASE %-15i\n' % (label, isubcase)
                #header[1] = '0    %-72s                                SUBCASE %-15i\n' % ('',isubcase)

                res_length = 0
                for res_type in res_types:
                    if isubcase in res_type:
                        result = res_type[isubcase]
                        if hasattr(result, 'write_op2'):
                            result.write_op2(op2)
                        res_length = max(len(result.__class__.__name__), res_length)
                        continue
                if res_length == 0:
                    return
                res_format = '%%-%is SUBCASE=%%i%%s' % res_length

                print("F06:")
                for res_type in res_types:
                    #print("res_type ", res_type)
                    header = ['', '']
                    #header[0] = '     %s\n' % subtitle
                    header[0] = '      %-126s\n' % subtitle
                    header[1] = '0     %-32s                                                                       SUBCASE %-15i\n \n' % (label, isubcase)
                    #print "res_type = ", res_type
                    if isubcase in res_type:
                        #header = copy.deepcopy(headerOld)  # fixes bug in case
                        #print("isubcase ", isubcase)
                        result = res_type[isubcase]
                        if result.nonlinear_factor is not None:
                            header.append('')
                        try:
                            element_name = ''
                            if hasattr(result, 'element_name'):
                                element_name = ' - ' + result.element_name
                            print(res_format % (result.__class__.__name__, isubcase, element_name))
                            self.page_num = result.write_f06(header, page_stamp, page_num=self.page_num, f=f06, is_mag_phase=False)
                            assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
                        except:
                            #print "result name = %r" % result.name()
                            raise
                        if delete_objects:
                            del result
                        self.page_num += 1
        f06.write(make_end(end_flag))
        f06.close()
        if write_op2:
            op2.close()
