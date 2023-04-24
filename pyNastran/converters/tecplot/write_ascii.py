from typing import TextIO
import numpy as np
from cpylog import SimpleLogger
from .zone import Zone

def write_ascii_header(title: str,
                       is_x: bool, is_y: bool, is_z: bool,
                       res_types: list[str],
                       variables: list[str]) -> tuple[str, np.ndarray]:
    """"
    "tecplot geometry and solution file"
    """
    if '"' in title or "'" in title:
        msg = 'TITLE = %s\n' % title
    else:
        msg = 'TITLE = "%s"\n' % title

    variables_ = []
    result_indices_to_write = []
    #if is_x:
        #ivar = variables.index('X')
        #variables_.append('X')
        #result_indices_to_write.append(ivar)
    #if is_y:
        #ivar = variables.index('Y')
        #variables_.append('Y')
        #result_indices_to_write.append(ivar)
    #if is_z:
        #ivar = variables.index('Z')
        #variables_.append('Z')
        #result_indices_to_write.append(ivar)


    #msg += 'ZONE T="%s"\n' % r'\"processor 1\"'
    #print(f'res_types = {res_types}')
    #print(f'vars = {variables}')
    for ivar, var in enumerate(res_types):
        if var not in variables:
            raise RuntimeError(f'var={var!r} not in variables={variables}')
        #print(f'adding {var}')
        result_indices_to_write.append(variables.index(var))
    ivars = np.unique(result_indices_to_write)
    ivars.sort()
    assert len(result_indices_to_write) == len(ivars)

    for ivar in ivars:
        var = variables[ivar]
        variables_.append(var)

    #if len(variables_):
    msg += 'VARIABLES = ' + ''.join(f'"{var}"\n' for var in variables_)
    #print(f'ivars = {ivars}')
    assert len(ivars) > 0, ivars
    return msg, ivars

def write_ascii_tecplot_zone(tecplot_file: TextIO,
                             zone: Zone,
                             ivars: np.ndarray,
                             log: SimpleLogger,
                             adjust_nids: bool) -> None:
    nnodes = zone.nnodes
    nelements = zone.nelements
    (is_structured, is_unstructured, is_points, zone_type,
     is_tris, is_quads, is_tets, is_hexas) = zone.determine_element_type()
    #print(is_structured, is_unstructured, is_points, zone_type)
    #print(is_tris, is_quads, is_tets, is_hexas)

    if is_unstructured:
        zone.write_unstructured_zone(
            tecplot_file, ivars, is_points, nnodes, nelements, zone_type, log,
            is_tris, is_quads, is_tets, is_hexas, adjust_nids=adjust_nids)
    elif is_structured:
        zone.write_structured_zone(
            tecplot_file, ivars, log,
            zone.headers_dict, adjust_nids=adjust_nids)
    else:  # pragma: no cover
        raise RuntimeError('only structured/unstructured')
