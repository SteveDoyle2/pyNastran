from pyNastran.converters.cart3d.cart3d import Cart3D, read_cart3d
from pyNastran.converters.tecplot.tecplot import Tecplot
from pyNastran.converters.tecplot.zone import Zone, TecplotDict

def cart3d_to_tecplot(cart3d_filename, tecplot_filename, log=None, debug=False):
    """
    Converts Cart3d to Tecplot
    """
    if isinstance(cart3d_filename, Cart3D):
        model = cart3d_filename
    else:
        model = read_cart3d(cart3d_filename, log=log, debug=debug)

    tecplot = Tecplot()
    tecplot.log = model.log
    variables = ['X', 'Y', 'Z', 'Region']
    zonetype = 'FETRIANGLE'
    header_dict = TecplotDict({
        'VARIABLES': variables,
        'ZONETYPE': zonetype,
    })
    if 1:  # pragma: no cover
        print(model.points.shape)
        print(model.elements.shape)
        print(model.regions.shape)
        nelement = len(model.elements)
        element_results = model.regions.reshape(nelement, 1)
        name = '???'
        strand_id = -1
        data_packing = 'BLOCK'
        zone = Zone.set_zone_from_360(log,
                              header_dict,
                              variables,
                              name,
                              zonetype,
                              data_packing,
                              strand_id,
                              tris=model.elements+1,
                              quads=None,
                              tets=None, hexas=None,
                              zone_data=model.points,
                              element_data=element_results)
    else:
        zone = Zone(model.log)
        #zone.headers_dict['VARIABLES'] = variables
        #zone.headers_dict['ZONETYPE'] = zonetype
        zone.headers_dict = header_dict
        zone.zone_data = model.points
        zone.tri_elements = model.elements + 1
    tecplot.zones = [zone]

    tecplot.write_tecplot(tecplot_filename, adjust_nids=False)
    return tecplot
