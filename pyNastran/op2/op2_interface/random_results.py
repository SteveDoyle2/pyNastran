from typing import Dict, Any
import numpy as np

class RandomObjects:
    prefix = ''
    postfix = ''
    def __init__(self):
        self.displacements = {}
        self.velocities = {}
        self.accelerations = {}
        self.load_vectors = {}
        self.spc_forces = {}
        self.mpc_forces = {}

        self.crod_force = {}
        self.conrod_force = {}
        self.ctube_force = {}

        self.cbar_force = {}
        self.cbeam_force = {}

        self.cbush_stress = {}
        self.cbush_strain = {}

        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}
        self.cbar_stress = {}
        self.cbeam_stress = {}

        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}
        self.cbar_strain = {}
        self.cbeam_strain = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}

        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.celas1_stress = {}
        self.celas2_stress = {}
        self.celas3_stress = {}
        self.celas4_stress = {}

        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        self.celas1_force = {}
        self.celas2_force = {}
        self.celas3_force = {}
        self.celas4_force = {}

        self.ctria3_force = {}
        self.ctria6_force = {}
        self.ctriar_force = {}
        self.cquad4_force = {}
        self.cquad8_force = {}
        self.cquadr_force = {}

        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}

        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

        self.cbend_stress = {}
        self.cbend_strain = {}
        self.cbend_force = {}

        self.cshear_stress = {}
        self.cshear_strain = {}
        self.cshear_force = {}

        self.cbush_force = {}
        self.cdamp1_force = {}
        self.cdamp2_force = {}
        self.cdamp3_force = {}
        self.cdamp4_force = {}
        self.cvisc_force = {}

        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

        self.cquad4_composite_strain = {}
        self.cquad8_composite_strain = {}
        self.cquadr_composite_strain = {}
        self.ctria3_composite_strain = {}
        self.ctria6_composite_strain = {}
        self.ctriar_composite_strain = {}

    def get_table_types(self):
        tables = [
            'displacements', 'velocities', 'accelerations',
            'load_vectors', 'spc_forces', 'mpc_forces',

            'celas1_force', 'celas2_force', 'celas3_force', 'celas4_force',
            'crod_force', 'conrod_force', 'ctube_force',
            'cbar_force', 'cbeam_force',
            'cquad4_force', 'cquad8_force', 'cquadr_force',
            'ctria3_force', 'ctria6_force', 'ctriar_force',

            'celas1_stress', 'celas2_stress', 'celas3_stress', 'celas4_stress',
            'crod_stress', 'conrod_stress', 'ctube_stress',
            'cbar_stress', 'cbeam_stress',
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',
            'ctetra_stress', 'cpenta_stress', 'chexa_stress',

            'celas1_strain', 'celas2_strain', 'celas3_strain', 'celas4_strain',
            'crod_strain', 'conrod_strain', 'ctube_strain',
            'cbar_strain', 'cbeam_strain',
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',
            'ctetra_strain', 'cpenta_strain', 'chexa_strain',

            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',

            'cbend_stress', 'cbend_strain', 'cbend_force',
            'cbush_stress', 'cbush_strain',
            'cshear_stress', 'cshear_strain', 'cshear_force',

            'cbush_force',
            'cdamp1_force', 'cdamp2_force', 'cdamp3_force', 'cdamp4_force',
            'cvisc_force',

        ]
        return [self.prefix + table + self.postfix for table in tables]

class PSDObjects():
    """storage class for the ATO objects"""
    prefix = 'psds.'
    postfix = ''
    def __init__(self):
        self.displacements = {}
        self.velocities = {}
        self.accelerations = {}
        self.spc_forces = {}
        self.load_vectors = {}
        self.force = {}
        self.stress = {}
        self.strain = {}

    def get_table_types(self):
        tables = self._tables()
        return [self.prefix + table + self.postfix for table in tables]

    def _tables(self):
        tables = [
            'displacements', 'velocities', 'accelerations',
            'spc_forces', 'load_vectors',
            'force', 'stress', 'strain',
        ]
        return tables

    def get_results(self):
        tables = self._tables()
        results = {}
        for table in tables:
            result = getattr(self, table)
            if result:
                results[table] = result
        return results

    def get_stats(self, short=True):
        msg = ''
        psds_dict = self.get_results()

        for result_type, slot in psds_dict.items():
            npsds = len(slot)

            if short:
                msg += f'op2_results.psds.{result_type}; n={npsds}\n'
            else:
                ipsd = 0
                msg += f'op2_results.psds.{result_type}:\n'
                msg += f'  # (subtitle, analysis_code, stress_strain_flag, node, dof)\n'
                for key in slot:
                    msg += f'  {key}\n'
                    if ipsd == 10:
                        msg += f'  ... npsds={npsds}\n'
                        break
                    ipsd += 1
                msg += '\n'
        return msg

    def get_psds_by_subtitles(self) -> Dict[Any, Any]:
        psd_results = self.get_results()
        if not psd_results:
            return {}

        from collections import defaultdict
        psds_subtitle = defaultdict(dict)
        for res_type, psds in psd_results.items():
            for key, psd in psds.items():
                (subtitle, nid, dof) = key
                psds_subtitle[subtitle][(res_type, nid, dof)] = psd
        return psds_subtitle

    def plot(self):
        psds_subtitle = self.get_psds_by_subtitles()
        if not psds_subtitle:
            return

        import matplotlib.pyplot as plt
        for subtitle, psds in psds_subtitle.items():
            fig = plt.figure(1)
            for (res_type, nid, dof), psd in psds.items():
                freqs, psd = psd[:, 0], psd[:, 1]
                plt.plot(freqs, psd, name=f'(restype,nid,dof)=({res_type}, {nid}, {dof})')
            plt.legend()
        plt.show()

    def write_f06(self, f06):
        psds_subtitle = self.get_psds_by_subtitles()
        if not psds_subtitle:
            return

        psd_type_map = {
            'displacements' : 'DISP',
            'velocities' : 'VELO',
            'accelerations' : 'ACCE',
            'load_vectors' : 'OLOAD',
            'spc_forces' : 'SPCF',
            'force' : 'EL FOR',
            'stress' : 'EL STR',
            'strain' : 'STRAIN',
        }
        from scipy.integrate import trapz
        for subtitle, psds in psds_subtitle.items():
            f06.write(subtitle + '\n')
            f06.write('0                             X Y - O U T P U T  S U M M A R Y  ( A U T O  O R  P S D F )\n')
            f06.write('0 PLOT  CURVE FRAME    CURVE ID./       RMS       NO. POSITIVE   XMIN FOR   XMAX FOR   YMIN FOR    X FOR     YMAX FOR    X FOR*\n')
            f06.write('  TYPE   TYPE   NO.  PANEL  : GRID ID    VALUE        CROSSINGS   ALL DATA   ALL DATA   ALL DATA     YMIN     ALL DATA     YMAX\n')

            #fig = plt.figure(1)
            for (res_type, nid, dof), psd in psds.items():
                try:
                    psd_type = psd_type_map[res_type]
                except KeyError:
                    raise NotImplementedError(f'res_type = {res_type}')
                    #psd_type = analysis_code

                #rms_value = 2.879461E+00
                #no_crossings = 2.879461E+00
                #no_crossings = np.nan
                freqs, psd = psd[:, 0], psd[:, 1]
                #plt.plot(freqs, psd, name=f'(restype,nid,dof)=({res_type}, {nid}, {dof})')
                ymin = psd.min()
                ymax = psd.max()
                imin = np.where(psd == ymin)[0][0]
                imax = np.where(psd == ymax)[0][0]
                xmin = freqs[imin]
                xmax = freqs[imax]
                fmin = freqs.min()
                fmax = freqs.max()

                # If you want the RMS value, this is computed as RMS = SQRT(SUM(PSD*DF)) and,
                # where DF is the spectral resolution, where you integarate from Fmin to Fmax,
                # i.e. your lowest and highest analysis frequency of interest, respectively.
                psd_f = trapz(psd, freqs)
                rms = psd_f ** 0.5
                if psd_f == 0.0:
                    # really this is nan, but that's Nastran for you
                    no_crossings = 0.0
                else:
                    f2_psd_f = trapz(freqs**2 * psd, freqs)
                    no_crossings = (f2_psd_f / psd_f) ** 0.5  # Hz

                #print('ymin=%s ymax=%s xmin=%s xmax=%s fmin=%s fmax=%s' % (ymin, ymax, xmin, xmax, fmin, fmax))
                #'0                             X Y - O U T P U T  S U M M A R Y  ( A U T O  O R  P S D F )'
                #'0 PLOT  CURVE FRAME    CURVE ID./       RMS       NO. POSITIVE   XMIN FOR   XMAX FOR   YMIN FOR    X FOR     YMAX FOR    X FOR*'
                #'    TYPE   TYPE   NO.  PANEL  : GRID ID    VALUE        CROSSINGS   ALL DATA   ALL DATA   ALL DATA     YMIN     ALL DATA     YMAX'
                #'    PSDF ACCE       0  9400703(  5)    2.879461E+00  8.191217E+02  2.000E+01  2.000E+03  4.476E-06  7.900E+01  1.474E+00  3.980E+01'
                f06.write('0                                      \n')
                f06.write(f'  PSDF {psd_type:6s}     0 {nid:8d}( {dof:2d})    {rms:8.6E}  {no_crossings:9.6E}  {fmin:9.3E}  {fmax:9.3E}  {ymin:9.3E}  {xmin:9.3E}  {ymax:9.3E}  {xmax:9.3E}\n')
            #plt.legend()
        #plt.show()

class AutoCorrelationObjects(RandomObjects):
    """storage class for the ATO objects"""
    prefix = 'ato.'
    #postfix = ''

class PowerSpectralDensityObjects(RandomObjects):
    """storage class for the PSD objects"""
    prefix = 'psd.'
    #postfix = ''

class RootMeansSquareObjects(RandomObjects):
    """storage class for the RMS objects"""
    prefix = 'rms.'
    #postfix = ''

class CumulativeRootMeansSquareObjects(RandomObjects):
    """storage class for the CRMS objects"""
    prefix = 'crm.'
    #postfix = ''

class NumberOfCrossingsObjects(RandomObjects):
    """storage class for the NO objects"""
    prefix = 'no.'
    #postfix = ''

class RAECONS:
    """storage class for the RAECONS objects"""
    def __init__(self):
        self.ctria3_strain = {}
        self.cquad4_strain = {}
        self.chexa_strain = {}

    def get_table_types(self):
        tables = [
            'chexa_strain',
            'ctria3_strain', 'cquad4_strain',
        ]
        return ['RAECONS.' + table for table in tables]

class RASCONS:
    """storage class for the RASCONS objects"""
    def __init__(self):
        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}

        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}

        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

    def get_table_types(self):
        tables = [
            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',

            'ctetra_stress', 'chexa_stress', 'cpenta_stress',
            'ctetra_strain', 'chexa_strain', 'cpenta_strain',
        ]
        return ['RASCONS.' + table for table in tables]

class RAPCONS:
    """storage class for the RAPCONS objects"""
    def __init__(self):
        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

    def get_table_types(self):
        tables = [
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'cquadr_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',
            'ctriar_composite_stress',
            #'cquad4_composite_strain',
            #'cquad8_composite_strain',
            #'cquadr_composite_strain',
            #'ctria3_composite_strain',
            #'ctria6_composite_strain',
            #'ctriar_composite_strain',
        ]
        return ['RAPCONS.' + table for table in tables]

class RAPEATC:
    """storage class for the RAPEATC objects"""
    def __init__(self):
        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

    def get_table_types(self):
        tables = [
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'cquadr_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',
            'ctriar_composite_stress',

            #'cquad4_composite_strain',
            #'cquad8_composite_strain',
            #'cquadr_composite_strain',
            #'ctria3_composite_strain',
            #'ctria6_composite_strain',
            #'ctriar_composite_strain',
        ]
        return ['RAPEATC.' + table for table in tables]

class RAFCONS:
    """storage class for the RAFCONS objects"""
    def __init__(self):
        self.cbar_force = {}
        self.cquad4_force = {}
        self.cbush_force = {}

    def get_table_types(self):
        tables = [
            'cbar_force',
            'cquad4_force',
            'cbush_force',
        ]
        return ['RAFCONS.' + table for table in tables]

class RAGCONS:
    """storage class for the RAGCONS objects"""
    def __init__(self):
        self.grid_point_forces = {}
    def get_table_types(self):
        tables = [
            'grid_point_forces',
        ]
        return ['RAGCONS.' + table for table in tables]

class RAGEATC:
    """storage class for the RAGEATC objects"""
    def __init__(self):
        self.grid_point_forces = {}

    def get_table_types(self):
        tables = [
            'grid_point_forces',
        ]
        return ['RAGEATC.' + table for table in tables]


class RANCONS:
    """storage class for the RANCONS objects"""
    def __init__(self):
        self.cbar_strain_energy = {}
        self.cbush_strain_energy = {}
        self.chexa_strain_energy = {}
        self.ctria3_strain_energy = {}
        self.cquad4_strain_energy = {}

    def get_table_types(self):
        tables = [
            'cbar_strain_energy', 'cbush_strain_energy',
            'chexa_strain_energy',
            'ctria3_strain_energy', 'cquad4_strain_energy',
        ]
        return ['RANCONS.' + table for table in tables]

class RADEFFM:
    """storage class for the RADEFFM objects"""
    def __init__(self):
        self.eigenvectors = {}
    def get_table_types(self):
        tables = [
            'eigenvectors',
        ]
        return ['RADEFFM.' + table for table in tables]


class RADCONS:
    def __init__(self):
        self.eigenvectors = {}

    def get_table_types(self):
        tables = [
            'eigenvectors',
        ]
        return ['RADCONS.' + table for table in tables]


class RADEATC:
    """storage class for the RADEATC objects"""
    def __init__(self):
        self.eigenvectors = {}

    def get_table_types(self):
        tables = [
            'eigenvectors',
        ]
        return ['RADEATC.' + table for table in tables]


class RANEATC:
    """storage class for the RANEATC objects"""
    def __init__(self):
        self.cbar_strain_energy = {}
        self.cbush_strain_energy = {}
        self.chexa_strain_energy = {}
        self.ctria3_strain_energy = {}
        self.cquad4_strain_energy = {}

    def get_table_types(self):
        tables = [
            'cbar_strain_energy', 'cbush_strain_energy',
            'chexa_strain_energy',
            'ctria3_strain_energy', 'cquad4_strain_energy',
        ]
        return ['RANEATC.' + table for table in tables]


class ROUGV1:
    """storage class for the ROUGV1 objects"""
    def __init__(self):
        self.displacements = {}
        self.velocities = {}
        self.accelerations = {}
        self.eigenvectors = {}

    def get_table_types(self):
        tables = [
            'displacements', 'velocities', 'accelerations', 'eigenvectors',
        ]
        return ['ROUGV1.' + table for table in tables]

class RAFEATC:
    """storage class for the RAFEATC objects"""
    def __init__(self):
        self.cbar_force = {}
        self.cquad4_force = {}
        self.cbush_force = {}

    def get_table_types(self):
        tables = [
            'cbar_force',
            'cquad4_force',
            'cbush_force',
        ]
        return ['RAFEATC.' + table for table in tables]


class RASEATC:
    """storage class for the RASEATC objects"""
    def __init__(self):
        self.chexa_stress = {}
        self.cquad4_stress = {}

    def get_table_types(self):
        tables = [
            'chexa_stress',
            'cquad4_stress',
        ]
        return ['RASEATC.' + table for table in tables]

class RAEEATC:
    """storage class for the RAEEATC objects"""
    def __init__(self):
        self.chexa_strain = {}
        self.ctria3_strain = {}
        self.cquad4_strain = {}

    def get_table_types(self):
        tables = [
            'chexa_strain',
            'ctria3_strain', 'cquad4_strain',
        ]
        return ['RAEEATC.' + table for table in tables]
