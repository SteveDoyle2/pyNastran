from __future__ import print_function
from numpy import (array,
                   where, where)

from pyNastran.bdf.dev_vectorized.utils import unique2d
#from pyNastran.bdf.dev_vectorized.cards.elements.solid.ctetra4 import volume4
#from pyNastran.bdf.dev_vectorized.cards.elements.solid.chexa8 import quad_area_centroid
#from pyNastran.bdf.dev_vectorized.cards.elements.solid.cpenta6 import tri_area_centroid

#from pyNastran.bdf.dev_vectorized.cards.elements.shell.cquad4 import _cquad4_normal_A
#from pyNastran.bdf.dev_vectorized.cards.elements.shell.ctria3 import _ctria3_normal_A

class Properties(object):
    def __init__(self, model):
        """
        Defines the Properties object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.np = 0

        #: stores PSHELL, PCOMP, PCOMPG
        self.properties_shell = model.properties_shell

        # shear
        #: stores PSHEAR
        self.pshear = model.pshear

        # spring
        self.pelas = model.pelas

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

    def build(self):
        ptypes = self._get_property_types(nlimit=False)
        self.n = 0
        for props in ptypes:
            props.build()
            self.np += props.n

        pids = check_duplicate('property_id', ptypes)

        self.property_ids = array(list(pids), dtype='int32')
        self.property_ids.sort()
        self.property_groups = self.build_groups(ptypes, 'property_id')

    def get_properties(self, property_ids=None):
        return self.model.elements.get_properties(property_ids)

    #def _get_property_types(self, nlimit=True):
        #"""
        #:param nlimit: limit the outputs to objects with data
        #"""
        #types = [self.prod, self.pelas,
                 #self.properties_bar.pbar, self.properties_bar.pbarl,
                 #self.properties_beam.pbeam, self.properties_beam.pbeaml,
                 #self.pshear,

                 ##self.properties_shell,
                 #self.properties_shell.pshell,
                 #self.properties_shell.pcomp,
                 #self.properties_shell.pcompg,

                 ##self.properties_solid,
                 #self.properties_solid.psolid,
                 ##self.properties_solid.plsolid,
                 #]
        #if nlimit:
            #types2 = []
            #for etype in types:
                #if etype.n > 0:
                    #types2.append(etype)
            #types = types2
        #return types

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

def check_duplicate(name, objs):
    unique_vals = set([])
    for obj in objs:
        if hasattr(obj, name):
            vals = getattr(obj, name)
            if len(vals):
                #self.model.log.debug("%s vals = %s for class %s" % (name, vals, obj.__class__.__name__))
                unique_vals.update(list(vals))
            #print unique_vals
        else:
            #print "  %s has no %s"  % (obj.__class__.__name__, name)
            pass
    #print "unique %s = %s\n" %(name, unique_vals)
    if len(unique_vals) == 0:
        raise RuntimeError("unique %s = %s" %(name, unique_vals))
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
    pid_eType = unique2d(pid_data[:, 1:])

    data2 = {}
    eTypeMap = {
        1 : 'CROD', 5: 'CONROD',
        2 : 'CBEAM', 3 : 'CBAR',
        4 : 'CSHEAR',
        10 : 'CELAS1', 11 : 'CELAS2', 12 : 'CELAS3', 13 : 'CELAS4',
        73 : 'CTRIA3', 144 : 'CQUAD4',

        60 : 'CTETRA4', 61 : 'CTETRA10',
        62 : 'CPENTA6', 63 : 'CPENTA15',
        64 : 'CHEXA8', 65 : 'CHEXA20',
    }

    #self.model.log.debug("pid_eType = \n%s\n" % str(pid_eType))
    for (pid, eType) in pid_eType:
        if pid not in elements.property_ids:
            print('Property pid=%s does not exist' % pid)
            #continue
        i = where(pid_data[:, 1] == pid)[0]
        #self.model.log.debug("pid=%i eType=%s Step #1=> \n%s\n" % (pid, eType, pid_data[i, :]))
        j = where(pid_data[i, 2] == eType)[0]
        eids = pid_data[i[j], 0]
        #self.model.log.debug("pid=%i eType=%s eids=%s Step #2=> \n%s\n" % (pid, eType, eids, pid_data[i[j], :]))
        eType = eTypeMap[eType]
        data2[(pid, eType)] = eids
    return data2
