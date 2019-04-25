from __future__ import print_function
from six import string_types
#from six.moves import StringIO

from pyNastran.bdf.bdf import BDF

#def get_model(bdf_filename):
    #if isinstance(bdf_filename, str):
        #model = read_bdf(bdf_filename=bdf_filename, validate=True, xref=True,
                        #punch=False, skip_cards=None,
                        #read_cards=None,
                        #encoding=None, log=None,
                        #debug=True, mode='msc')
    #else:
        #model = bdf_filename
    #return model

def get_bdf_model(bdf_filename, xref=True, log=None, debug=False):
    # type: (Union[str, BDF], bool, Any, bool) -> BDF
    if isinstance(bdf_filename, string_types):
        model = BDF(log=log, debug=debug)
        model.read_bdf(bdf_filename, xref=xref)
    elif isinstance(bdf_filename, BDF):
        model = bdf_filename
        if xref:
            model.cross_reference(xref=xref)
    #elif isinstance(bdf_filename, StringIO):
        #model = BDF(log=log, debug=debug)
        #model.read_bdf(bdf_filename, xref=xref)
    else:
        raise NotImplementedError(bdf_filename)
    return model
