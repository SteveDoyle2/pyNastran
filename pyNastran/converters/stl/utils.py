from numpy import vstack
from cpylog import get_logger2
from pyNastran.converters.stl.stl import STL


# merge_tecplot_files(tecplot_filenames, tecplot_filename_out=None, log=None):
def merge_stl_files(stl_filenames, stl_out_filename=None, remove_bad_elements=False,
                    is_binary=True, float_fmt='%6.12f', log=None):
    """
    Combines multiple STLs into a single file

    Parameters
    ----------
    stl_filenames : List[str, str, ...]
        list of stl filenames or a string filename
        (useful for removing bad elements)
    remove_bad_elements : bool; default=False
        should elements with invalid normals be removed?
    stl_out_filename : str; default=None -> no writing
        string of stl output filename
    is_binary : bool; default=True
        should the output file be binary
    float_fmt : str; default='%6.12f'
        the ascii float format

    Returns
    -------
    stl : STL()
        the stl object
    """
    if isinstance(stl_filenames, str):
        stl_filenames = [stl_filenames]
    assert isinstance(stl_filenames, (list, tuple)), type(stl_filenames)
    assert len(stl_filenames) > 0, stl_filenames

    log = get_logger2(log=log, debug=False, encoding='utf-8')
    if len(stl_filenames) == 1:
        model = STL(log=log)
        model.read_stl(stl_filenames[0])
        if remove_bad_elements:
            model.remove_elements_with_bad_normals()
        if stl_out_filename is not None:
            model.write_stl(stl_out_filename, is_binary=is_binary)
        return model

    nodes = []
    elements = []

    n0 = 0
    for i, fname in enumerate(stl_filenames):
        log.debug(f'reading file {i}: {fname}')
        model = STL(log=log)  # TODO: you shouldn't need to to reinstantiate the STL
        model.read_stl(fname)
        nnodes = model.nodes.shape[0]
        nodes.append(model.nodes)
        elements.append(model.elements + n0)
        n0 += nnodes

    model.nodes = vstack(nodes)
    model.elements = vstack(elements)

    if remove_bad_elements:
        model.remove_elements_with_bad_normals()

    if stl_out_filename is not None:
        model.write_stl(stl_out_filename, is_binary=is_binary, float_fmt=float_fmt)
    return model
