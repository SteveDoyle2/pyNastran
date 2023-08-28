from typing import Optional
import numpy as np
from pyNastran.converters.aflr.ugrid.ugrid_reader import UGRID, read_ugrid
from pyNastran.converters.tecplot.tecplot import Tecplot, Zone

def get_ugrid_model(ugrid_filename: str, log=None, debug: bool=False) -> UGRID:
    """helper method for loading UGRID models

    Parameters
    ----------
    ugrid_filename : varies
        str : the input UGRID filename
        UGRID : the UGRID object

    Returns
    -------
    ugrid_model : UGRID()
        the UGRID object

    """
    if isinstance(ugrid_filename, str):
        #assert os.path.exists(ugrid_filename), '%r doesnt exist' % ugrid_filename
        model = read_ugrid(ugrid_filename=ugrid_filename,
                           encoding=None, log=log, debug=debug,
                           read_shells=True, read_solids=True, check=True)
        #UGRID(log=log, debug=debug)
        #model.read_ugrid(ugrid_filename)
    else:
        model = ugrid_filename
        assert isinstance(model, UGRID), 'expected UGRID; type(model)=%s' % type(model)
    return model

def ugrid3d_to_tecplot_filename(ugrid_filename: str, tecplot_filename: str,
                                log=None, debug: bool=False) -> Tecplot:
    """
    Converts a UGRID to a Tecplot ASCII file.

    Parameters
    ----------
    ugrid_filename : varies
        str : the input UGRID filename
        UGRID : the UGRID object
    tecplot_filename : str
        the output Tecplot filename
    log : logger; default=None
        a logger object
    debug : bool; default=False
        developer debug

    Returns
    -------
    tecplot_model : Tecplot()
        the Tecplot object
    """
    ugrid_model = get_ugrid_model(ugrid_filename, log=log, debug=debug)
    tecplot = write_tecplot(ugrid_model, tecplot_filename)
    return tecplot

def write_tecplot(ugrid_model: UGRID, tecplot_filename: str) -> Tecplot:
    ugrid_model.check_hanging_nodes()
    tecplot, zone = ugrid_to_tecplot(ugrid_model)
    tecplot.write_tecplot(tecplot_filename, adjust_nids=True)  # is adjust correct???
    return tecplot

def ugrid_to_tecplot(ugrid_filename: str,
                     tecplot_filename: Optional[str]=None,
                     log=None, debug: bool=False) -> tuple[Tecplot, Zone]:
    """
    Converts a UGRID to a Tecplot ASCII file.

    Parameters
    ----------
    ugrid_filename : varies
        str : the input UGRID filename
        UGRID : the UGRID object
    tecplot_filename : str
        the output Tecplot filename
    log : logger; default=None
        a logger object
    debug : bool; default=False
        developer debug

    Returns
    -------
    tecplot_model : Tecplot()
        the Tecplot object
    """
    ugrid_model = get_ugrid_model(ugrid_filename, log=log, debug=debug)
    #nnodes = len(ugrid_model.nodes)
    #nodes = zeros((nnodes, 3), dtype='float64')
    ugrid_model.check_hanging_nodes()

    ntets = len(ugrid_model.tets)
    nhexas = len(ugrid_model.hexas)
    non_tets = len(ugrid_model.penta5s) + len(ugrid_model.penta6s) + len(ugrid_model.hexas)
    non_hexas = len(ugrid_model.tets) + len(ugrid_model.penta5s) + len(ugrid_model.penta6s)
    assert ntets + non_tets > 0, 'nsolids=%s' % (ntets + non_tets)

    tecplot = Tecplot(log=ugrid_model.log, debug=debug)
    zone = Zone(ugrid_model.log)
    zone.headers_dict['VARIABLES'] = ['X', 'Y', 'Z']
    zone.zone_data = ugrid_model.nodes

    if ntets and non_tets == 0:
        elements = ugrid_model.tets
        zone.tet_elements = elements - 1
        zone_type = 'FETETRAHEDRON'
    elif nhexas and non_hexas == 0:
        elements = ugrid_model.hexas
        zone.hexa_elements = elements - 1
        zone_type = 'FEBRICK'
    elif non_tets:
        elements_list = []
        for element in ugrid_model.tets:
            n1, n2, n3, n4 = element
            elements_list.append([n1, n2, n3, n4,
                                  n4, n4, n4, n4])
        for element in ugrid_model.penta5s:
            n1, n2, n3, n4, n5 = element
            elements_list.append([n1, n2, n3, n4,
                                  n5, n5, n5, n5])
        for element in ugrid_model.penta6s:
            n1, n2, n3, n4, n5, n6 = element
            elements_list.append([n1, n2, n3, n4,
                                  n5, n6, n6, n6])
        for element in ugrid_model.hexas:
            n1, n2, n3, n4, n5, n6, n7, n8 = element
            elements_list.append([n1, n2, n3, n4,
                                  n5, n6, n7, n8])
        elements = np.array(elements_list, dtype='int32') - 1
        zone.hexa_elements = elements
        zone_type = 'FEBRICK'
    else:
        raise RuntimeError()
    zone.headers_dict['ZONETYPE'] = zone_type

    if tecplot_filename is not None:
        tecplot.write_tecplot(tecplot_filename)
    tecplot.zones = [zone]
    str(zone)
    return tecplot, zone


def main(): # pragma: no cover
    import sys
    if len(sys.argv) != 3:
        msg = ('number of arguments must be 2; ugrid_filename, tecplot_filename;'
               ' nargs=%s; args=%s' % (
                   len(sys.argv[1:]), sys.argv[1:]))
        raise RuntimeError(msg)
    ugrid_filename = sys.argv[1]
    tecplot_filename = sys.argv[2]
    ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
