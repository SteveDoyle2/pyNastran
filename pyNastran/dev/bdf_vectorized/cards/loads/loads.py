from collections import defaultdict
from io import StringIO
from numpy import array, union1d
from pyNastran.dev.bdf_vectorized.cards.loads.load import LOAD


#class LOAD:
    #def __init__(self, model):
        #self.load = defaultdict(list)

class VectorizedCardDict:
    def __len__(self):
        return self.n

    def allocate(self, ncards):
        pass

    def __init__(self, model):
        """
        Defines the object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.type = type
        self.n = 0
        self.model = model
        self.keys = None
        self._cards = []
        self._comments = []

    #def add_card(self, card, comment=''):
        #prop = vPCOMP(card, comment=comment)
        #self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    def _get_objs(self, obj_dict, obj_keys_default, obj_keys):
        objs = []
        if obj_keys is None:
            ids = obj_keys_default
        for idi in obj_dict:
            obj = obj_dict[idi]
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

    def items(self):
        for key, value in self._objs.items():
            yield key, value

    def __repr__(self):
        file_obj = StringIO()
        file_obj.write('<%s object> n=%s\n' % (self.type, self.n))
        self.write_card(file_obj)
        return file_obj.getvalue()

class LOADs(VectorizedCardDict):
    type = 'LOAD'
    def __init__(self, model):
        VectorizedCardDict.__init__(self, model)
        self._objs = defaultdict(list)

    def get_load_ids(self):
        return array(self._objs.keys(), dtype='int32')

    def slice_by_load_id(self, load_id):
        return

    def add_card(self, card, comment=''):
        load = LOAD(self.model)
        load.add_from_bdf(card, comment)
        #print('load.load_id =%s' % load.load_id)
        #print('type(self.loads.load) =%s' % type(self.loads.load))
        self._objs[load.load_id].append(load)
        #self.loads.load.add(card, comment=comment)
        self.n += 1

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if load_id is None:
            load_id = self.get_load_ids()
        #if isinstance(load_id, integer_types):
            #load_id = array(load_id)
        if len(load_id.shape) == 0:
            load_id = array([load_id])

        #print('load_ids xx = %s' % load_id, type(load_id), load_id.shape)
        for lid in sorted(load_id):
            self.model.log.debug('load_id = %s' % lid)
            for load in self._objs[lid]:
                load.write_card(bdf_file, size=size, is_double=is_double, load_id=lid)

class Loads:
    def __init__(self, model):
        self.model = model

        #: stores LOAD, FORCE, MOMENT, etc.
        self.load = model.load
        self.dload = model.dload
        #self.dload = defaultdict(list)
        #self.loadset = LOADSET(model)

        self.force = model.force
        self.force1 = model.force1
        self.force2 = model.force2
        self.moment = model.moment
        self.moment1 = model.moment1
        self.moment2 = model.moment2
        self.grav = model.grav
        self.rforce = model.rforce

        self.pload = model.pload
        self.pload1 = model.pload1
        self.pload2 = model.pload2
        #self.pload3 = model.pload3
        self.pload4 = model.pload4

        self.ploadx1 = model.ploadx1

        self.tload1 = model.tload1
        self.tload2 = model.tload2
        self.delay = model.delay

        self.rload1 = model.rload1
        #self.rload2 = model.rload2
        self.dphase = model.dphase
        self.darea = model.darea

    def _get_load_types(self, nlimit):
        types = [
            # static
            self.load,  #self.dload,
            self.force, self.force1, self.force2,
            self.moment, self.moment1, self.moment2,
            self.grav, self.rforce,
            self.pload, self.pload1, self.pload2, self.pload4, self.ploadx1,
            #self.accel, self.accel1,

            # dynamic
            self.tload1, self.tload2,
            self.rload1, #self.rload2,

            # thermal
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

    def write_card(self, bdf_file, size=8, is_double=False, sort_data=False):
        sort_data = True
        bdf_file.write('$ LOADS\n')
        if sort_data:
            types = self._get_load_types(nlimit=True)
            all_load_ids = array([], dtype='int32')
            for load_type in types:
                #print('*load-type = %s' % load_type.type)
                all_load_ids = union1d(all_load_ids, load_type.get_load_ids())
            self.model.log.debug('all_load_ids= %s' % all_load_ids)

            for load_id in sorted(all_load_ids):
                for load_type in types:
                    #self.model.log.debug('load_type = %s' % load_type)
                    #self.model.log.debug('load_type.load_id = %s' % load_type.load_id)
                    if load_id in load_type.load_id:
                        #self.model.log.debug('  load_id = %s' % load_id)
                        assert load_id > 0, load_id
                        #print('******** load-type = %s' % load_type.type)
                        load_type.write_card(bdf_file, size=size, is_double=is_double,
                                             load_id=load_id)

        else:
            #self.loadcase.write_card(bdf_file, size)
            for load_id, loads in sorted(self.load.items()):
                for load in loads:
                    load.write_card(bdf_file, size)

            for load_id, loads in sorted(self.dload.items()):
                for load in loads:
                    load.write_card(bdf_file, size)

            #for load_id, loads in sorted(self.sload.items()):
                #for load in loads:
                    #load.write_card(bdf_file, size)

            #for load_id, loads in sorted(self.lseq.items()):
                #for load in loads:
                    #load.write_card(bdf_file, size)

            #self.loadset.write_card(bdf_file, size)
            self.force.write_card(bdf_file, size)
            #self.force1.write_card(bdf_file, size)
            #self.force2.write_card(bdf_file, size)
            self.moment.write_card(bdf_file, size)
            #self.moment1.write_card(bdf_file, size)
            #self.moment2.write_card(bdf_file, size)

            self.pload.write_card(bdf_file, size)
            self.pload1.write_card(bdf_file, size)
            self.pload2.write_card(bdf_file, size)
            #self.pload3.write_card(bdf_file, size)
            self.pload4.write_card(bdf_file, size)

            self.ploadx1.write_card(bdf_file, size)
            self.grav.write_card(bdf_file, size)
            self.rforce.write_card(bdf_file, size)

            #self.accel1.write_card(bdf_file, size)

            self.tload1.write_card(bdf_file, size)
            self.tload2.write_card(bdf_file, size)

            #self.rload1.write_card(bdf_file, size)
            #self.rload2.write_card(bdf_file, size)
        self.delay.write_card(bdf_file, size)
        self.dphase.write_card(bdf_file, size)
        self.darea.write_card(bdf_file, size)

            #self.randps.write_card(bdf_file, size)

            # DAREA
        #self.model.log.debug('done with loads')
