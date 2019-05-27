import os
import multiprocessing as mp
from numpy import savetxt, arange, zeros

from pyNastran.converters.usm3d.usm3d_reader import Usm3d

def get_n_list(dirname, model_name):
    """
    gets files of the form:
     - modelname + '_xxx.flo'

    """
    flo_filenames = os.listdir(dirname)

    n_list = []
    # get the max N value
    for flo_filename in flo_filenames:
        base, ext = os.path.splitext(flo_filename)
        if ext == '.flo' and '.aux.' not in flo_filename:
            if '_' not in flo_filename:
                print(flo_filename)
                continue
            n = int(base.split('_')[-1])
            n_list.append(n)
    n_list.sort()
    return n_list

def get_flo_files(dirname, model_name, nstart=0, nlimit=None, include_dirname_in_path=True):
    """get the flo files in ascending order"""
    if dirname == '':
        dirname = os.getcwd()
    n_list = get_n_list(dirname, model_name)
    if nlimit is None:
        nlimit = n_list[-1]  # inclusive

    # handles case of user didn't define every data point
    n_list2 = n_list[nstart:nlimit+1]
    #for ni in n_list:
        #if nstart <= ni <= nlimit:
            #n_list2.append(ni)

    if include_dirname_in_path:
        flo_filenames = [os.path.join(dirname, '%s_%i.flo' % (model_name, n)) for n in n_list2]
    else:
        flo_filenames = ['%s_%i.flo' % (model_name, n) for n in n_list2]
    return n_list2, flo_filenames

def get_flo_files_from_n(dirname, model_name, n_list, include_dirname_in_path=True):
    """
    get the flo files in ascending order
    """
    if dirname == '':
        dirname = os.getcwd()

    if include_dirname_in_path:
        flo_filenames = [os.path.join(dirname, '%s_%i.flo' % (model_name, n)) for n in n_list]
    else:
        flo_filenames = ['%s_%i.flo' % (model_name, n) for n in n_list]
    return flo_filenames


def run_time_acc(dirname, model_name, node_ids, nstart=0, nlimit=None, num_cpus=8):
    """
    node ids start at index=0
    """
    model = Usm3d()

    n_list, flo_filenames = get_flo_files(dirname, model_name, nstart, nlimit)
    nlimit = len(flo_filenames)
    assert nlimit > 0, 'nfiles=%s' % (nlimit)
    #print("nmax =", nlimit)

    Cp = {}
    Mach = {}
    T = {}
    U = {}
    V = {}
    W = {}
    p = {}
    rhoU = {}

    # initialize the arrays
    for node_id in node_ids:
        Cp[node_id] = zeros(nlimit, 'float64')
        Mach[node_id] = zeros(nlimit, 'float64')
        T[node_id] = zeros(nlimit, 'float64')
        U[node_id] = zeros(nlimit, 'float64')
        V[node_id] = zeros(nlimit, 'float64')
        W[node_id] = zeros(nlimit, 'float64')
        p[node_id] = zeros(nlimit, 'float64')
        rhoU[node_id] = zeros(nlimit, 'float64')

    #flo_filenames = flo_filenames[:10]
    #num_cpus = 1
    if num_cpus == 1:
        model = Usm3d()
        for flo_filename in flo_filenames:
            node_ids2, loads = model.read_flo(flo_filename, node_ids=node_ids)
            for i, node_id in enumerate(node_ids):
                Cp[node_id][i] = loads['Cp'][i]
                Mach[node_id][i] = loads['Mach'][i]
                T[node_id][i] = loads['T'][i]
                U[node_id][i] = loads['U'][i]
                V[node_id][i] = loads['V'][i]
                W[node_id][i] = loads['W'][i]
                p[node_id][i] = loads['p'][i]
                rhoU[node_id][i] = loads['rhoU'][i]
    else:
        pool = mp.Pool(num_cpus)
        result = pool.imap(_loads_func,
                           [(flo_filename, node_ids) for flo_filename in flo_filenames])
        n = 0
        for j, loads in enumerate(result):
            if j % 500 == 0:
                print("n =", j, flo_filenames[j])
            for i, node_id in enumerate(node_ids):
                #print("Cp[node=%s] =%s" % (node_id, loads['Cp'][i]))
                Cp[node_id][j] = loads['Cp'][i]
                Mach[node_id][j] = loads['Mach'][i]
                T[node_id][j] = loads['T'][i]
                U[node_id][j] = loads['U'][i]
                V[node_id][j] = loads['V'][i]
                W[node_id][j] = loads['W'][i]
                p[node_id][j] = loads['p'][i]
                rhoU[node_id][j] = loads['rhoU'][i]
                n += 1
        pool.close()
        pool.join()
    return n_list, Cp, Mach, T, U, V, W, p, rhoU

def _loads_func(data):
    flo_filename, node_ids = data
    model = Usm3d()
    node_ids2, loads = model.read_flo(flo_filename, node_ids=node_ids)
    return loads

def write_loads(csv_filename, loads, node_id):
    (Cp, Mach, T, U, V, W, p, rhoU) = loads
    #print("loads.keys() = ", sorted(loads.keys()))
    f = open(csv_filename, 'wb')
    dt = 1.0
    t = arange(len(Cp[node_id])) * dt  # broken...
    f.write('time\t')
    savetxt(f, t, delimiter='', newline=',')
    f.write('\n')
    for node_id, Cpi in sorted(Cp.items()):
        f.write("\nnode_id=%i\n" % node_id)

        f.write('Cp[%s],' % node_id)
        savetxt(f, Cpi, delimiter='', newline=',')

        f.write('\np[%s],' % node_id)
        savetxt(f, p[node_id], delimiter='', newline=',')
        f.write('\n\n')
    f.close()


def main():
    dirname = r'W:\usm3d\usm3d\spw4'
    model_name = 'spw4'
    node_ids = [517]
    loads = run_time_acc(dirname, model_name, node_ids)
    write_loads('usm3d.csv', loads, node_ids)


if __name__ == '__main__':  # pragma: no cover
    main()
