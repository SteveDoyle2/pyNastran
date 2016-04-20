import os
import copy
import shutil
from six import itervalues, iteritems

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.op2.op2 import read_op2
from pyNastran.utils.nastran_utils import run_nastran

def fully_stressed_design(bdf_filename, keywords=None,
                          niterations_max=2, alpha=0.9):
    """
    Optimizes shell thickness for minimum weight (ONLY shells)

    Parameters
    ----------
    bdf_filename : str; BDF()
        the BDF filename or model

    Returns
    -------
    desvars : dict[id]=values
        the "optimization history of the design variables
    """
    force = True

    iteration = 0
    niterations_max = 10

    if isinstance(bdf_filename, str):
        model = read_bdf(bdf_filename)
    elif isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        raise TypeError(bdf_filename)

    doptparm = model.doptprm
    if doptparm is not None:
        if 'FSDALP' in doptparm.params:
            alpha = doptparm.params['FSDALP']
        else:
            alpha = doptparm.defaults['FSDALP']
        if not isinstance(alpha, float):
            msg = 'FSDALP on DOPTPARM must be an integer; FSDALP=%r' % (alpha)
            raise TypeError(msg)
        if not(0. < niterations_max <= 1.):
            msg = 'FSDALP on DOPTPARM must be between (0. < n <= 1.0); FSDALP=%s' % (alpha)
            raise ValueError(msg)

        if 'FSDMAX' in doptparm.params:
            niterations_max = doptparm.params['FSDMAX']
        else:
            niterations_max = doptparm.defaults['FSDMAX']

        if not isinstance(niterations_max, int):
            msg = 'FSDMAX on DOPTPARM must be an integer; FSDMAX=%r' % (niterations_max)
            raise TypeError(msg)
        if niterations_max <= 0:
            msg = 'FSDMAX on DOPTPARM must be > 0; FSDMAX=%s' % (niterations_max)
            raise ValueError(msg)
    else:
        niterations_max = 2
        alpha = 0.9


    dresps_to_consider, desvars_to_consider, dvprels_to_consider = get_inputs(model)

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

            #def compute_critical_stress(results, subcases)
        #def compute_critical_stress(results, subcases):
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

            # PROD area
            # PSHELL/PCOMP thickness
            # PSHEAR thickness
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
            #tnew = (oi/omax)**alpha * ti_old
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

def get_inputs(model):
    # PROD area
    # PSHELL/PCOMP thickness
    # PSHEAR thickness
    dresps_to_consider = []
    dresp_rtypes_to_skip = ['WEIGHT', 'DISP', 'VOLUME']
    dresp_rtypes_to_consider = ['VON-MIS']

    assert len(model.dresps) >= 1, len(model.dresps)
    for dresp in model.dresps:
        if dresp.type == 'DRESP1':
            if dresp.response_type in dresp_types_to_skip:
                continue
            elif dresp.response_type in ['STRESS', 'STRAIN']:
                pass
            else:
                raise RuntimeError(dresp)
        else:
            raise RuntimeError(dresp)

    assert len(model.dconadds) == 0, len(model.dconadds)
    #for dconstr in model.dconstrs:
        #if dconstr.

    for dvprel in model.dvprels:
        if dvprel.type == 'DVPREL1':
            if dvprel.Type == 'PROD':
                if dvprel.pNameFid == 'A':
                    dvprels_to_consider.append(dvprel)
                else:
                    raise RuntimeError(dvprel)
            elif dvprel.Type == 'PSHELL':
                if dvprel.pNameFid == 'T':
                    dvprels_to_consider.append(dvprel)
                else:
                    raise RuntimeError(dvprel)
            elif dvprel.Type == 'PSHEAR':
                if dvprel.pNameFid == 'T':
                    dvprels_to_consider.append(dvprel)
                else:
                    raise RuntimeError(dvprel)
            elif dvprel.Type == 'PCOMP':
                if dvprel.pNameFid.startswith('T'): # thickness; TODO: what about theta?
                    dvprels_to_consider.append(dvprel)
                elif dvprel.pNameFid in ['ZO']:
                    continue
                else:
                    raise RuntimeError(dvprel)
            elif dvprel.Type == 'PBARL':
                if 'DIM' in dvprel.pNameFid: # increased inertia -> lower stress.  TODO: is this always true?
                    pass
                else:
                    raise RuntimeError(dvprel)
            else:
                raise RuntimeError(dvprel)

            if len(dvprel.dvids) == 1:
                desvars_to_consider.append(dvprel.dvids[1])
                dvprels_to_consider.append(dvprel)
            #for dvid, coeff in zip(dvprel.dvids, dvprel.coeffs):
                #if coeff == 1.0

        elif dvprel.type == 'DVPREL2':
            continue
        else:
            raise RuntimeError(dvprel)
    return dresps_to_consider, desvars_to_consider, dvprels_to_consider

def main():
    regions = {
        # pid : [tmin, tmax, ovm_min, ovm_max]
        2 : [0.05, 1.],
        9 : [0.05, 1.],
        1 : [0.05, 1.],
    }
    bdf_filename = 'fem.bdf'
    fully_stressed_design(bdf_filename, regions, keywords=keywords)