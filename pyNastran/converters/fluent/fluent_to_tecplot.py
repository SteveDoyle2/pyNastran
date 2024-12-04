#import os
import numpy as np
from pyNastran.utils import PathLike
from pyNastran.converters.fluent.fluent import read_fluent
from pyNastran.converters.tecplot.tecplot import Tecplot
from pyNastran.converters.tecplot.zone import Zone, TecplotDict

def fluent_to_tecplot(fluent_filename: PathLike,
                      tecplot_filename: PathLike='') -> Tecplot:
    """
    not well tested

    considers:
      - non-sequential node ids
      - filters unused nodes
    """
    fluent_model = read_fluent(fluent_filename)
    log = fluent_model.log
    node_id = fluent_model.node_id
    xyz = fluent_model.xyz
    result_element_id = fluent_model.result_element_id
    titles = fluent_model.titles
    results = fluent_model.results
    quads = fluent_model.quads  # eid, pid, n1, n2, n3, n4
    tris = fluent_model.tris    # eid, pid, n1, n2, n3

    tecplot = Tecplot()
    tecplot.log = fluent_model.log

    # drop the element id
    variables = ['X', 'Y', 'Z', 'Region'] + titles[1:].tolist()
    nvars = len(variables)

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

        unodes = np.unique(tri_nodes.ravel())
        inid = np.searchsorted(node_id, unodes)
        jtri_nodes = np.searchsorted(unodes, tri_nodes)
        #-------------

        ires_tri = np.searchsorted(result_element_id, tri_eids)
        res_tri = results[ires_tri, :]
        assert res_tri.ndim == 2, res_tri.shape
        res_tri = np.column_stack((tri_regions, res_tri))
        assert res_tri.shape == (ntri, nvars-3)
        del tris
        tri_zone = Zone.set_zone_from_360(
            log, header_dict, variables,
            name, zonetype, data_packing,
            strand_id,
            tris=jtri_nodes,
            quads=None,
            tets=None, hexas=None,
            zone_data=xyz[inid, :],
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

        unodes = np.unique(quad_nodes.ravel())
        inid = np.searchsorted(node_id, unodes)
        jquad_nodes = np.searchsorted(unodes, quad_nodes)
        #-------------

        ires_quad = np.searchsorted(result_element_id, quad_eids)
        res_quad = results[ires_quad, :]
        assert res_quad.ndim == 2, res_quad.shape
        res_quad = np.column_stack((quad_regions, res_quad))
        assert res_quad.shape == (nquad, nvars-3), res_quad.shape

        quad_zone = Zone.set_zone_from_360(
            log, header_dict, variables,
            name, zonetype, data_packing,
            strand_id,
            tris=None,
            quads=jquad_nodes,
            tets=None, hexas=None,
            zone_data=xyz[inid, :],
            element_data=res_quad)
        tecplot.zones.append(quad_zone)

    if tecplot_filename:
        tecplot.write_tecplot(tecplot_filename)
    return tecplot
