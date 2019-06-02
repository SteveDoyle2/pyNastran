import numpy as np

from pyNastran.femutils.utils import unique2d
#from pyNastran.dev.bdf_vectorized.cards.elements.solid.ctetra4 import volume4
#from pyNastran.dev.bdf_vectorized.cards.elements.solid.chexa8 import quad_area_centroid
#from pyNastran.dev.bdf_vectorized.cards.elements.solid.cpenta6 import tri_area_centroid

#from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquad4 import _cquad4_normal_A
#from pyNastran.dev.bdf_vectorized.cards.elements.shell.ctria3 import _ctria3_normal_A
from pyNastran.dev.bdf_vectorized.cards.elements.utils import build_groups #, asarray

class Properties:
    def __init__(self, model):
        """
        Defines the Properties object.

        Parameters
        ----------
        model : BDF
           the BDF object

        """
        self.model = model
        self.nproperties = 0

        #: stores PSHELL, PCOMP, PCOMPG
        self.properties_shell = model.properties_shell

        # shear
        #: stores PSHEAR
        self.pshear = model.pshear

        # spring
        self.pelas = model.pelas

        # bush
        self.pbush = model.pbush

        # rods
        #self.conrod = model.conrod
        #self.crod = model.crod
        self.prod = model.prod

        # mass
        #: stores CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4, CMASS5, PMASS
        self.mass = model.mass

        # bars
        #: stores PBAR, PBARL
        self.properties_bar = model.properties_bar

        # beams
        #: stores PBEAM, PBEAML
        self.properties_beam = model.properties_beam

        # solids
        #: stores PSOLID, PLSOLID
        self.properties_solid = model.properties_solid

        # created by this class
        self.property_ids = None
        self.n = None
        self.property_groups = None

    def build(self):
        ptypes = self._get_property_types(nlimit=False)
        self.n = 0
        for props in ptypes:
            assert props is not None, props
            #props.build()
            self.nproperties += props.n

        pids = check_duplicate('property_id', ptypes, self.model.log)

        self.property_ids = np.array(list(pids), dtype='int32')
        self.property_ids.sort()
        self.property_groups = build_groups(ptypes, 'property_id')

    def get_properties(self, property_ids=None):
        return self.model.elements.get_properties(property_ids)

    def _get_property_types(self, nlimit=True):
        """
        Parameters
        ----------
        nlimit : bool; default=True
            limit the outputs to objects with data
        """
        types = [
            self.prod, self.pelas, self.pbush,
            self.properties_bar.pbar, self.properties_bar.pbarl,
            self.properties_beam.pbeam, self.properties_beam.pbeaml,
            self.pshear,

            #self.properties_shell,
            self.properties_shell.pshell,
            self.properties_shell.pcomp,
            self.properties_shell.pcompg,

            #self.properties_solid,
            self.properties_solid.psolid,
            #self.properties_solid.plsolid,
        ]
        if nlimit:
            types2 = []
            for etype in types:
                if etype.n > 0:
                    types2.append(etype)
            types = types2
        return types

    #def get_property_typemap(self):
        #TypeMap = {
            #'PELAS'  : self.pelas,
            #'PROD'  : self.prod,
            #'PSHEAR'  : self.pshear,

            #'PBAR'  : self.properties_bar.pbar,
            #'PBARL'  : self.properties_bar.pbarl,

            #'PBEAM'  : self.properties_beam.pbeam,
            #'PBEAML'  : self.properties_beam.pbeaml,
            ##'PBUSH'  : self.pbush,

            #'PSHELL'  : self.properties_shell.pshell,
            #'PCOMP'  : self.properties_shell.pcomp,
            #'PCOMPG'  : self.properties_shell.pcompg,

            #'PSOLID' : self.properties_solid.psolid,
        #}
        #return TypeMap

    def __len__(self):
        return self.model.elements.np

    def __iter__(self):
        pids = self.model.elements.property_ids
        for pid in pids:
            yield pid

    def values(self):
        pids = self.model.elements.property_ids
        for pid in pids:
            yield self.__getitem__(pid)

    def items(self):
        pids = self.model.elements.property_ids
        for pid in pids:
            yield pid, self.__getitem__(pid)

    def __getitem__(self, property_ids):
        return self.model.elements.get_properties(property_ids)

def check_duplicate(name, objs, log):
    unique_vals = set()
    for obj in objs:
        if hasattr(obj, name):
            vals = getattr(obj, name)
            if len(vals):
                #self.model.log.debug("%s vals = %s for class %s" % (
                    #name, vals, obj.__class__.__name__))
                unique_vals.update(list(vals))
            #print unique_vals
        else:
            #print("  %s has no %s"  % (obj.__class__.__name__, name))
            pass
    #print("unique %s = %s\n" %(name, unique_vals))
    if len(unique_vals) == 0:
        log.info("unique %s = %s" %(name, unique_vals)) # fails for CONRODs
        #raise RuntimeError
    #print('unique %s = %s' % (name, unique_vals))
    return unique_vals

def group_elements_by_property_type_and_element_type(elements, pid_data):
    """
    group elements of the same type by property type
       same element type & different property id (e.g. CTRIA3 PSHELL/PCOMP)  -> different group
       different element type & same property id (e.g. CTRIA3/CQUAD4 PSHELL) -> different group
       same element & same type -> same group

    we do this in order to think about one property at a time and not
    have to do a lot of special work to handle different methods for
    getting the mass
    """
    # find unique groups
    #print("pid_data = \n%s\n" % str(pid_data))
    pid_elementnum = unique2d(pid_data[:, 1:])

    data2 = {}
    etype_map = {
        1 : 'CROD', 5: 'CONROD',
        2 : 'CBEAM', 3 : 'CBAR',
        4 : 'CSHEAR',
        10 : 'CELAS1', 11 : 'CELAS2', 12 : 'CELAS3', 13 : 'CELAS4',
        73 : 'CTRIA3', 144 : 'CQUAD4',

        60 : 'CTETRA4', 61 : 'CTETRA10',
        62 : 'CPENTA6', 63 : 'CPENTA15',
        64 : 'CHEXA8', 65 : 'CHEXA20',
    }

    #self.model.log.debug("pid_elementnum = \n%s\n" % str(pid_elementnum))
    for (pid, element_num) in pid_elementnum:
        if pid not in elements.property_ids:
            print('Property pid=%s does not exist' % pid)
            #continue
        i = np.where(pid_data[:, 1] == pid)[0]
        #self.model.log.debug("pid=%i element_num=%s Step #1=> \n%s\n" % (
            #pid, element_num, pid_data[i, :]))
        j = np.where(pid_data[i, 2] == element_num)[0]
        eids = pid_data[i[j], 0]
        #self.model.log.debug("pid=%i element_num=%s eids=%s Step #2=> \n%s\n" % (
            #pid, element_num, eids, pid_data[i[j], :]))
        element_type = etype_map[element_num]
        data2[(pid, element_type)] = eids
    return data2
