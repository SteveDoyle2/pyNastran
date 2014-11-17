from six import  iteritems
from six.moves import StringIO
from collections import defaultdict
from numpy import array, intersect1d, union1d, asarray

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter

from .load import LOAD
from .dload import DLOAD
from .dload import DLOAD as LSEQ
from .dload import DLOAD as SLOAD
from .grav import GRAV
from .force import FORCE
from .moment import MOMENT

from .force1 import FORCE1
from .force2 import FORCE2
#from .moment1 import MOMENT1
#from .moment2 import MOMENT2

# ACCEL1
# PLOAD3
# DAREA
# TLOAD1
# TLOAD2
# RLOAD1
# RLOAD2
# RANDPS

# loads
from .pload  import PLOAD
from .pload1 import PLOAD1
from .pload2 import PLOAD2
#from .pload3 import PLOAD3
from .pload4 import PLOAD4

from .ploadx1 import PLOADX1
#from .grav import GRAV

from .rforce import RFORCE
#from .sload import SLOAD

#from .loadcase import LoadCase
#from .loadset import LOADSET

#class LOAD(object):
    #def __init__(self, model):
        #self.load = defaultdict(list)

class VectorizedCardDict(object):
    def __len__(self):
        return self.n

    def allocate(self, ncards):
        pass

    def __init__(self, model):
        """
        Defines the object.

        :param self: the object
        :param model: the BDF object
        :param cards: the list of cards
        """
        self.type = type
        self.n = 0
        self.model = model
        self.keys = None
        self._cards = []
        self._comments = []

    #def add(self, card, comment):
        #prop = vPCOMP(card, comment=comment)
        #self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    def _get_objs(self, obj_dict, obj_keys_default, obj_keys):
        objs = []
        if obj_keys is None:
            ids = obj_keys_default
        for id in obj_dict:
            obj = obj_dict[id]
            objs.append(obj)
        return objs

    #def __getitem__(self, load_id):
        #load_id, int_flag = slice_to_iter(load_id)
        #return self.slice_by_load_id(load_id)

    def __contains__(self, key):
        """TODO: should check against unique values"""
        if key in self._objs:
            return True
        return False

    def __getitem__(self, load_id):
        return self._objs[load_id]

    def iteritems(self):
        for key, value in iteritems(self._objs):
            yield key, value

    def __repr__(self):
        f = StringIO()
        f.write('<%s object> n=%s\n' % (self.type, self.n))
        self.write_bdf(f)
        return f.getvalue()

class LOADs(VectorizedCardDict):
    type = 'LOAD'
    def __init__(self, model):
        VectorizedCardDict.__init__(self, model)
        self._objs = defaultdict(list)

    def get_load_ids(self):
        return array(self._objs.keys(), dtype='int32')

    def slice_by_load_id(self, load_id):
        return

    def add(self, card, comment):
        load = LOAD(self.model)
        load.add_from_bdf(card, comment)
        #print('load.load_id =%s' % load.load_id)
        #print('type(self.loads.load) =%s' % type(self.loads.load))
        self._objs[load.load_id].append(load)
        #self.loads.load.add(card, comment=comment)
        self.n += 1

    def write_bdf(self, f, size=8, is_double=False, load_id=None):
        if load_id is None:
            load_id = self.get_load_ids()
        #if isinstance(load_id, int):
            #load_id = array(load_id)
        if len(load_id.shape) == 0:
            load_id = array([load_id])

        #print('load_ids xx = %s' % load_id, type(load_id), load_id.shape)
        for lid in sorted(load_id):
            for load in self._objs[lid]:
                load.write_bdf(f, size=size, is_double=is_double, load_id=lid)

class Loads(object):
    def __init__(self, model):
        self.model = model

        #: stores LOAD, FORCE, MOMENT, etc.
        self.load = LOADs(model)
        self.dload = LOADs(model)
        #self.dload = defaultdict(list)
        #self.loadset = LOADSET(model)

        self.force = FORCE(model)
        self.force1 = FORCE1(model)
        self.force2 = FORCE2(model)
        self.moment = MOMENT(model)
        #self.moment1 = MOMENT1(model)
        #self.moment2 = MOMENT2(model)
        self.grav = GRAV(model)
        self.rforce = RFORCE(model)

        self.pload = PLOAD(model)
        self.pload1 = PLOAD1(model)
        self.pload2 = PLOAD2(model)
        #self.pload3 = PLOAD3(model)
        self.pload4 = PLOAD4(model)

        self.ploadx1 = PLOADX1(model)

    def _get_load_types(self, nlimit):
        types = [
            self.load,  #self.dload,
            self.force, self.force1, self.force2,
            self.grav, self.rforce,
            self.pload, self.pload1, self.pload2, self.pload4, self.ploadx1,
        ]
        if nlimit:
            types2 = []
            for ltype in types:
                if ltype.n > 0:
                    types2.append(ltype)
            types = types2
        return types

    def allocate(self, card_count):
        load_types = self._get_load_types(nlimit=False)
        for load_type in load_types:
            if load_type.type in card_count:
                load_type.allocate(card_count[load_type.type])
                del card_count[load_type.type]

    def __getitem__(self, value):
        pass

    def write_bdf(self, f, size=8, is_double=False, sort_data=False):
        sort_data = True
        if sort_data:
            types = self._get_load_types(nlimit=True)
            all_load_ids = array([], dtype='int32')
            for load_type in types:
                #print('*load-type = %s' % load_type.type)
                all_load_ids = union1d(all_load_ids, load_type.get_load_ids())
            #print('load_ids= %s' % all_load_ids)

            for load_id in sorted(all_load_ids):
                for load_type in types:
                    if load_id in load_type:
                        #print('******** load-type = %s' % load_type.type)
                        load_type.write_bdf(f, size=size, is_double=is_double, load_id=load_id)

        else:
            #self.loadcase.write_bdf(f, size)
            for load_id, loads in sorted(iteritems(self.load)):
                for load in loads:
                    load.write_bdf(f, size)

            for load_id, loads in sorted(iteritems(self.dload)):
                for load in loads:
                    load.write_bdf(f, size)

            #for load_id, loads in sorted(iteritems(self.sload)):
                #for load in loads:
                    #load.write_bdf(f, size)

            #for load_id, loads in sorted(iteritems(self.lseq)):
                #for load in loads:
                    #load.write_bdf(f, size)

            #self.loadset.write_bdf(f, size)
            self.force.write_bdf(f, size)
            #self.force1.write_bdf(f, size)
            #self.force2.write_bdf(f, size)
            self.moment.write_bdf(f, size)
            #self.moment1.write_bdf(f, size)
            #self.moment2.write_bdf(f, size)

            self.pload.write_bdf(f, size)
            self.pload1.write_bdf(f, size)
            self.pload2.write_bdf(f, size)
            #self.pload3.write_bdf(f, size)
            self.pload4.write_bdf(f, size)

            self.ploadx1.write_bdf(f, size)
            self.grav.write_bdf(f, size)
            self.rforce.write_bdf(f, size)

            #self.accel1.write_bdf(f, size)

            #self.tload1.write_bdf(f, size)
            #self.tload2.write_bdf(f, size)
            #self.rload1.write_bdf(f, size)
            #self.rload2.write_bdf(f, size)
            #self.randps.write_bdf(f, size)

            # DAREA
