from io import StringIO
from pathlib import PurePath
from pyNastran.bdf.bdf import BDF, read_bdf

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

BDF_FILETYPE = BDF | str | StringIO | PurePath
def get_bdf_model(bdf_filename: BDF_FILETYPE,
                  xref: bool=True,
                  cards_to_skip=None,
                  cards_to_include=None,
                  validate: bool=True,
                  log=None, debug: bool=False) -> BDF:
    if isinstance(bdf_filename, (str, StringIO, PurePath)):
        #model = read_bdf(bdf_filename=bdf_filename, validate=True, xref=True,
                        #punch=False, skip_cards=None,
                        #read_cards=None,
                        #encoding=None, log=None,
                        #debug=True, mode='msc')
        model = read_bdf(
            bdf_filename,
            validate=validate,
            xref=xref,
            skip_cards=cards_to_skip,
            read_cards=cards_to_include,
            punch=False,
            #read_includes=True,
            save_file_structure=False,
            log=log, debug=debug)
    elif isinstance(bdf_filename, BDF):
        model = bdf_filename
        if xref:
            model.cross_reference(xref=xref)
    #elif isinstance(bdf_filename, StringIO):
        #model = BDF(log=log, debug=debug)
        #model.read_bdf(bdf_filename, xref=xref)
    else:  # pragma: no cover
        raise NotImplementedError(bdf_filename)
    return model
