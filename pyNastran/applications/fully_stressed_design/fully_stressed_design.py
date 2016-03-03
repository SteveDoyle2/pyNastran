import os
import copy
import shutil
from six import itervalues, iteritems

import numpy as np
from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import read_op2
from pyNastran.applications.aero_panel_buckling.run_patch_buckling_helper import run_nastran

def fully_stressed_design(bdf_filename, regions, keywords=None):
    """
    Optimizes shell thickness for minimum weight (ONLY shells)

    Parameters
    ----------
    bdf_filename : str
        the BDF filename
    regions : dict[pid]=list[float,float,float,float]
        the optimization dictionary

    Returns
    -------
    regions2 : dict[pid]=list[float,float,float,float]
        the "optimized" version of regions

    region = {
        pid : [tmin, tmax, ovm_min, ovm_max]
    }
    """
    force = True
    regions2 = copy.deepcopy(regions)

    iteration = 0
    niterations_max = 10
    model = read_bdf(bdf_filename)
    pid_to_eid = model.get_property_id_to_element_ids_map()

    bdf_filename2 = 'fem_baseline.bdf'
    op2_filename2 = 'fem_baseline.op2'
    shutil.copyfile(bdf_filename, bdf_filename2)

    while iteration < niterations_max:
        if not os.path.exists(op2_filename2) or force:
            run_nastran(bdf_filename2, keywords=keywords)
        results = read_op2(op2_filename2, combine=True, log=None,
                          debug=False,
                          debug_file=None,
                          build_dataframe=False,
                          skip_undefined_matrices=True,
                          mode='msc')

        isubcase = 1
        itime = 0

        stress_per_region = {}
        nopt = 0
        for pid, region in iteritems(regions):
            print('pid=%s region=%s' % (pid, region))
            (tmin, tmax, ovm_min, ovm_max) = region
            prop = model.properties[pid]
            told = prop.t
            eids_requested = pid_to_eid[pid]
            print('eids_requested[pid=%s] = %s' % (pid, eids_requested))

            stress = []
            eid_node = []
            for res in [results.cquad4_stress, results.ctria3_stress]:
                resi = res[isubcase]
                eid_nodei = resi.element_node
                #print('eid_nodei = %s' % (eid_nodei))
                eid = eid_nodei[:, 0]
                stress_data = resi.data
                eid_node.append(eid_nodei)


                # A
                #i = np.where(eid == eids_requested)

                # B
                #j = np.searchsorted(eid, eids_requested)
                #i = np.where(eid[j] == eids_requested)
                #j = np.in1d(eids_requested, eid) # A in B
                j = np.in1d(eid, eids_requested) # A in B
                i = np.where(j)
                #print('i = %s' % i) #[0]
                #print('j = %s' % j) #[0]
                #print('eid = %s' % eid) #[0]
                #print('eids_requested = %s' % eids_requested) #[0]
                if len(i) == 0:
                    continue
                #print('i=%s; ni=%s' % (i, len(i)))
                stress_datai = stress_data[itime, i, 7]
                #print('eids = %s' % eid[i])
                #print('stress_datai = %s' % stress_datai)
                stress.append(stress_datai)

            #print('stressA = %s' % stress)
            stress = np.hstack(stress)
            #print('stressB = %s' % stress)
            eid_node = np.vstack(eid_node)
            stress_max = stress.max()
            stress_min = stress.min()
            print('stress_min=%s' % stress_min)
            print('stress_max=%s' % stress_max)
            istress_max = np.where(stress == stress_max)[0]
            istress_min = np.where(stress == stress_min)[0]
            eid_max = eid_node[istress_max, 0]
            eid_min = eid_node[istress_min, 0]
            peak_stress = max(abs(stress_max), abs(stress_min))
            tnew = told * stress_max / ovm_max
            tnew = min(tmax, tnew)
            tnew = max(tmin, tnew)
            tratio = tnew / told
            if np.allclose(tratio, 1.):
                continue
            nopt += 1
            stress_per_region = [stress_min, stress_max,
                                 eid_min, eid_max]
            print('pid=%s' % pid)
            print('  stress_per_region (ovm_min, ovm_max, eid_min, eid_max) => %s' % stress_per_region)
            print('  told=%s tnew=%s tratio=%s\n' %(told, tnew, tratio))

            prop.t *= tratio
            prop.z1 *= tratio
            prop.z2 *= tratio

            regions2[pid] = [tnew, peak_stress, stress_max, eid_max, stress_min, eid_min] # t_new, ovm_new
        if nopt == 0:
            break

        #eid_node = np.hstack(eid_node)
        #s = np.hstack(eid_node)

        iteration += 1
        bdf_filename2 = 'fem_%i.bdf' % iteration
        op2_filename2 = 'fem_%i.op2' % iteration
        model.write_bdf(bdf_filename2)

    print('regions2 = %s' % regions2)
    return regions2

def main():
    regions = {
        # pid : [tmin, tmax, ovm_min, ovm_max]
        2 : [0.05, 1.],
        9 : [0.05, 1.],
        1 : [0.05, 1.],
    }
    bdf_filename = 'fem.bdf'
    fully_stressed_design(bdf_filename, regions, keywords=keywords)