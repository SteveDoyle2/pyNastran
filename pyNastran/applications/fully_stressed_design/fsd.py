from six import itervalues
import shutil

from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import read_op2

def fully_stressed_design(bdf_filename, regions, keywords=keywords):
    """
    Optimizes shell thickness for minimum weight (ONLY shells)

    region = {
        pid : [tmin, tmax, ovm_min, ovm_max]
    }
    """
    i = 0
    model = read_bdf(bdf_filename)
    pid_to_eid = model.get_property_id_to_element_ids_map()
    shutil.copyfile(bdf_filename, bdf_filename2)

    bdf_filename2 = 'fem_baseline.bdf'
    op2_filename2 = 'fem_baseline.op2'
    while i < niterations_max:
        bdf_filename2 = 'fem_%i.bdf' % i
        op2_filename2 = 'fem_%i.op2' % i
        results = read_op2(op2_filename, combine=True, log=None,
                          debug=False,
                          debug_file=None,
                          build_dataframe=False,
                          skip_undefined_matrices=True,
                          mode='msc')

        isubcase = 1
        itime = 0

        stress_per_region = {}
        for pid, (tmin, tmax, ovm_min, ovm_max) in itervalues(regions):
            prop = bdf.properties[pid]
            told = prop.t
            eids_requested = pid_to_eid[pid]

            stress = []
            eid_node = []
            for resi in [results.cquad4_stress, results.ctria3_stress]:
                eid_nodei = resi.element_node
                eid = eid_nodei[:, 0]
                stres_data = res.data
                eid_node.append(eid_nodei)
                i = where(eid == eids_requested)
                stress.append(stress_data[itime, i, 7])

            stress = np.vstack(stress)
            eid_node = np.vstack(eid_node)
            stress_max = stress.max()
            stress_min = stress.min()
            istress_max = where(stress == stress_max)[0]
            istress_min = where(stress == stress_min)[0]
            eid_max = eid_node[istress_max, 0]
            eid_min = eid_node[istress_min, 0]
            peak_stress = max(abs(stress_max), abs(stress_min))
            tnew = told * stress_max / ovm_max
            tnew = min(tmax, tnew)
            tnew = max(tmin, tnew)
            tratio = tnew / told

            stress_per_region = [stress_min, stress_max,
                                 eid_min, eid_max]
            print('pid=%s' % pid)
            print('  stress_per_region (ovm_min, ovm_max, eid_min, eid_max) => %s' % stress_per_region)
            print('  told=%s tnew=%s tratio=%s\n' %(told, tnew, tratio))

            prop.t *= tratio
            prop.z1 *= tratio
            prop.z2 *= tratio




        #eid_node = np.hstack(eid_node)
        #s = np.hstack(eid_node)

        i += 1
        model.write_bdf(bdf_filename2)

def main():
    regions = {
        # pid : [tmin, tmax, ovm_min, ovm_max]
        2 : [0.05, 1.],
        9 : [0.05, 1.],
        1 : [0.05, 1.],
    }
    bdf_filename = 'fem.bdf'
    fully_stressed_design(bdf_filename, regions, keywords=keywords)