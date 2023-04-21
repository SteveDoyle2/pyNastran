from typing import Optional
from numpy import vstack
from pyNastran.converters.tecplot.tecplot import Tecplot

def merge_tecplot_files(tecplot_filenames: list[str],
                        tecplot_filename_out: Optional[str]=None,
                        log=None) -> Tecplot:
    """merges one or more tecplot files"""
    assert isinstance(tecplot_filenames, (list, tuple)), type(tecplot_filenames)
    assert len(tecplot_filenames) > 0, tecplot_filenames

    xyz = []
    tri_elements = []
    quad_elements = []
    tet_elements = []
    hexa_elements = []
    results = []
    nnodes = 0

    model = Tecplot(log=log)
    if len(tecplot_filenames) == 1:
        model.read_tecplot(tecplot_filenames[0])
        if tecplot_filename_out is not None:
            model.write_tecplot(tecplot_filename_out)
        return model

    for tecplot_filename in tecplot_filenames:
        model.log.info('reading %s' % tecplot_filename)
        model.read_tecplot(tecplot_filename)
        for zone in model.zones:
            xyz.append(zone.xyz)
            if len(zone.tri_elements):
                tri_elements.append(zone.tri_elements + nnodes)
            if len(zone.quad_elements):
                quad_elements.append(zone.quad_elements + nnodes)
            if len(zone.tet_elements):
                tet_elements.append(zone.tet_elements + nnodes)
            if len(zone.hexa_elements):
                hexa_elements.append(zone.hexa_elements + nnodes)
            results.append(zone.nodal_results)
            nnodes += zone.nnodes

    model.xyz = vstack(xyz)
    if tri_elements:
        model.tri_elements = vstack(tri_elements)
    if quad_elements:
        model.quad_elements = vstack(quad_elements)
    if tet_elements:
        model.tet_elements = vstack(tet_elements)
    if hexa_elements:
        model.hexa_elements = vstack(hexa_elements)
    model.results = vstack(results)
    if tecplot_filename_out is not None:
        model.write_tecplot(tecplot_filename_out)
    return model
