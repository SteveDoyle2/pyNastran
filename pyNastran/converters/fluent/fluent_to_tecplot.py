#import os
import numpy as np
from pyNastran.utils import PathLike
from pyNastran.converters.fluent.fluent import read_fluent
from pyNastran.converters.tecplot.tecplot import Tecplot
from pyNastran.converters.tecplot.zone import Zone, TecplotDict

def fluent_to_tecplot(fluent_filename: PathLike,
                      tecplot_filename: PathLike='') -> Tecplot:
    """
    TODO: not tested
    TODO: probably need to consider non-sequential node ids
    TODO: filter unused nodes in
    """
    fluent_model = read_fluent(fluent_filename)
    log = fluent_model.log
    node_id = fluent_model.node_id
    xyz = fluent_model.xyz
    element_id = fluent_model.element_id
    titles = fluent_model.titles
    results = fluent_model.results
    quads = fluent_model.quads  # eid, pid, n1, n2, n3, n4
    tris = fluent_model.tris    # eid, pid, n1, n2, n3

    tecplot = Tecplot()
    tecplot.log = fluent_model.log

    # drop the element id
    #variables = ['Region'] + titles[1:].tolist()
    variables = ['X', 'Y', 'Z', 'Region'] + titles[1:].tolist()
    nvars = len(variables)
    #print('variables = ', variables)
    #print('nvars = ', nvars)

    data_packing = 'BLOCK'  # don't remember what this is
    name = '???'
    strand_id = -1

    ntri = len(tris)
    if len(tris):
        zonetype = 'FETRIANGLE'
        header_dict = TecplotDict({
            'VARIABLES': variables,
            'ZONETYPE': zonetype,
        })
        tri_eids = tris[:, 0]
        tri_regions = tris[:, 1].reshape(ntri, 1)
        tri_nodes = tris[:, 2:]
        ires_tri = np.searchsorted(element_id, tri_eids)
        res_tri = results[ires_tri, :]
        assert res_tri.ndim == 2, res_tri.shape
        #print('ires_tri.shape = ', ires_tri.shape)
        #print('res_tri.shape = ', res_tri.shape)
        #print('tri_regions.shape = ', tri_regions.shape)
        res_tri = np.column_stack((tri_regions, res_tri))
        assert res_tri.shape == (ntri, nvars-3) #, (tri_regions.shape, res_tri.shape, res_tri2.shape)
        del tris
        tri_zone = Zone.set_zone_from_360(
            log, header_dict, variables,
            name, zonetype, data_packing,
            strand_id,
            tris=tri_nodes + 1,
            quads=None,
            tets=None, hexas=None,
            zone_data=xyz,
            element_data=res_tri)
        tecplot.zones.append(tri_zone)

    nquad = len(quads)
    if nquad:
        zonetype = 'FEQUADRILATERAL'
        header_dict = TecplotDict({
            'VARIABLES': variables,
            'ZONETYPE': zonetype,
        })
        quad_eids = quads[:, 0]
        quad_regions = quads[:, 1].reshape(nquad, 1)
        quad_nodes = quads[:, 2:]
        ires_quad = np.searchsorted(element_id, quad_eids)
        res_quad = results[ires_quad, :]
        assert res_quad.ndim == 2, res_quad.shape
        res_quad = np.column_stack((quad_regions, res_quad))
        assert res_quad.shape == (nquad, nvars-3), res_quad.shape

        quad_zone = Zone.set_zone_from_360(
            log, header_dict, variables,
            name, zonetype, data_packing,
            strand_id,
            tris=None,
            quads=quad_nodes + 1,
            tets=None, hexas=None,
            zone_data=xyz,
            element_data=res_quad)
        tecplot.zones.append(quad_zone)

    if tecplot_filename:
        # TODO: fails b/c nodes aren't sequential?
        tecplot.write_tecplot(tecplot_filename)
    return tecplot
