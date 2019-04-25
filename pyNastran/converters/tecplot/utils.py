from numpy import vstack
from pyNastran.converters.tecplot.tecplot import Tecplot

def merge_tecplot_files(tecplot_filenames, tecplot_filename_out=None, log=None):
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
        xyz.append(model.xyz)
        if len(model.tri_elements):
            tri_elements.append(model.tri_elements + nnodes)
        if len(model.quad_elements):
            quad_elements.append(model.quad_elements + nnodes)
        if len(model.tet_elements):
            tet_elements.append(model.tet_elements + nnodes)
        if len(model.hexa_elements):
            hexa_elements.append(model.hexa_elements + nnodes)
        results.append(model.results)
        nnodes += model.nnodes

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


