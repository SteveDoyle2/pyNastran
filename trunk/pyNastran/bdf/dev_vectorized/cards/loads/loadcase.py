from six.moves import zip, StringIO
from numpy import unique, where

from collections import defaultdict

class LoadCase(object):
    def __init__(self):
        self.loads = defaultdict(list)

    def __getitem__(self, key):
        print("loadcase %s" % key)
        aaa
        return self.resolve(key)

    def add(self, load):
        if load.n:
            #print("load.n = %s" % load.n)
            load_ids = unique(load.load_id)
            for load_id in load_ids:
                if load.type in ['FORCE', 'FORCE1', 'FORCE2',
                                'MOMENT', 'MOMENT1', 'MOMENT2',
                                'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOADX1']:  # PLOAD4, RFORCE, GRAV
                    i = where(load_id == load.load_id)[0]
                    #print("i** =", i)
                    if len(i):
                        self.loads[load_id].append(load[i])
                else:
                    print("load.type %s" % load.type)
                    raise NotImplementedError(load.type)

    def add_reference(self, load):
        for load_id, loads in load.iteritems():
            #print("ref", loads)
            for load in loads:
                self.loads[load_id].append(load)

    def resolve(self, i):
        #print("key %s" % i)
        all_loads = self.loads[i]  # list of load objs

        #print("**********")
        f = StringIO.StringIO()
        for load in all_loads:
            load.write_bdf(f)
        #print(f.getvalue())
        #print("**********")

        all_loads_out = []
        all_loads_loop = all_loads

        i = 0
        while all_loads_loop:
            #print("--------------------------")
            #print("i = %s" % i)
            all_loads_loop = []
            #print "all_loads", all_loads
            for loads in all_loads:
                #print("a-loads =", loads)
                if isinstance(loads, tuple):
                    (sf, loads2) = loads
                    for load in loads2:
                        if load.type in ['FORCE', 'MOMENT']:
                            all_loads_out.append(sf * load)
                        else:
                            asdfff
                elif loads.type in ['FORCE', 'MOMENT', 'PLOAD1']:
                    all_loads_out.append(loads)

                else:
                    #print("*loads", loads)  # LOAD, DLOAD
                    load = loads
                    #for load in loads:
                    scale = load.scale
                    scale_factors = load.scale_factors
                    new_load_ids = load.load_ids
                    #print("new_load_ids =", new_load_ids)

                    for scale_factor, new_load_id in zip(scale_factors, new_load_ids):
                        sf = scale * scale_factor
                        new_all_loads = self.resolve(new_load_id)
                        all_loads_loop.append((sf, new_all_loads))
                        #print("added sf=%s new_load_id=%s" % (sf, new_load_id))
            #print("len(all_loads_loop) =", len(all_loads_loop))
            all_loads = all_loads_loop
            i += 1
            if i > 100:
                raise RuntimeError('i > 100')
        #print('------------------')
        #print("resolved LoadCase i=", i)
        f = StringIO.StringIO()
        for load in all_loads_out:
            load.write_bdf(f)
        #print(f.getvalue())
        #print('------------------')
        return all_loads_out