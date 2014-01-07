import os
from collections import defaultdict

from pyNastran.converters.usm3d.usm3d_reader import Usm3dReader

def get_flo_files(dirname, model_name, include_dirname_in_path=True):
    """
    get the flo files in ascending order
    """
    if dirname == '':
        dirname = os.getcwd()
    flo_filenames = os.listdir(dirname)

    n_list = []
    # get the max N value
    for flo_filename in flo_filenames:
        base, ext = os.path.splitext(flo_filename)
        if ext == '.flo':
            n = int(base.split('_')[-1])
            n_list.append(n)

    n_list.sort()
    if include_dirname_in_path:
        flo_filenames = [os.path.join(dirname, '%s_%i.flo' % (model_name, n)) for n in n_list]
    else:
        flo_filenames = ['%s_%i.flo' % (model_name, n) for n in n_list]
    return flo_filenames



def main():
    model = Usm3dReader()
    node_ids = [1000]

    dirname = r'W:\usm3d\usm3d\spw4\nan_2'
    model_name = 'spw4'

    flo_filenames = get_flo_files(dirname, model_name)
    print "flo_filenames =", flo_filenames

    Cp = {}
    Mach = {}
    T = {}
    U = {}
    V = {}
    W = {}
    p = {}
    rhoU = {}
    for node_id in node_ids:
        Cp[node_id] = []
        Mach[node_id] = []
        T[node_id] = []
        U[node_id] = []
        V[node_id] = []
        W[node_id] = []
        p[node_id] = []
        rhoU[node_id] = []

    for flo_filename in flo_filenames:
        #print('flo_filename =', flo_filename)
        node_ids2, loads = model.read_flo(flo_filename, node_ids=node_ids)
        for i, node_id in enumerate(node_ids):
            Cp[node_id].append(loads['Cp'][i])
            Mach[node_id].append(loads['Mach'][i])
            T[node_id].append(loads['T'][i])
            U[node_id].append(loads['U'][i])
            V[node_id].append(loads['V'][i])
            W[node_id].append(loads['W'][i])
            p[node_id].append(loads['p'][i])
            rhoU[node_id].append(loads['rhoU'][i])

    print "loads.keys() = ", sorted(loads.keys())
    for node_id, Cpi in sorted(Cp.iteritems()):
        print "node_id=%i Cp=%s" % (node_id, Cpi)

if __name__ == '__main__':
    main()