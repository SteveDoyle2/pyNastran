from __future__ import annotations
from typing import Any, cast, TYPE_CHECKING
from numpy import int32, int64
from pyNastran.op2.result_objects.op2_objects import ScalarObject

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2
    from cpylog import SimpleLogger
COMPARE_KEYS = (int, int32, int64, str, bytes)


def get_op2_stats(model: OP2, short: bool=False) -> str:
    """see OP2.get_op2_stats(...)"""
    msg = []
    #msg += model.op2_results.responses.get_stats(short=short)
    op2_results = model.op2_results

    msg.extend(write_params_stats(model.params))
    for key, weight in model.grid_point_weight.items():
        msg += weight.get_stats(key, short=short)

    msg += op2_results.psds.get_stats(short=short)

    table_types = model._get_table_types_testing()

    if short:
        msg += _get_op2_stats_short(model, table_types, model.log)
    else:
        msg += _get_op2_stats_full(model, table_types, model.log)

    if model.matrices:
        msg.append('matrices:\n')
        for unused_name, matrix in sorted(model.matrices.items()):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append('  ' + str(matrix) + '\n')

    if model.matdicts:
        msg.append('matdicts:\n')
        for unused_name, matrix_dict in sorted(model.matdicts.items()):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append('  ' + str(matrix_dict) + '\n')
    try:
        return ''.join(msg)
    except TypeError:
        for msgi in msg:
            print('TypeError...%r' % msgi.rstrip())
            assert isinstance(msgi, str), msgi
        raise
    except UnicodeDecodeError:
        for msgi in msg:
            print('UnicodeDecodeError...%r' % msgi.rstrip())
            assert isinstance(msgi, str), msgi
        raise


def _get_op2_stats_short(model: OP2, table_types: list[str],
                         log: SimpleLogger) -> list[str]:
    """helper for get_op2_stats(...)"""
    msg = []
    handled_previously = ['params', 'grid_point_weight', 'psds', 'cstm']
    no_data_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
    for table_type in table_types:
        #table_type_print = ''
        if table_type in handled_previously:
            continue
        if table_type in ['gpdt', 'bgpdt', 'eqexin', 'monitor1', 'monitor3'] or table_type.startswith('responses.'):
            obj = model.get_result(table_type)
            if obj is None:
                continue
            elif isinstance(obj, dict):
                msg.extend(_get_op2_results_stats_dict(obj, table_type, short=True))
                continue

            obj = cast(ScalarObject, obj)
            stats: str = obj.get_stats(short=True)
            msg.extend(f'op2_results.{table_type}: ' + stats)  # TODO: a hack...not quite right...
            continue

        # # and not table_type.startswith('responses.')
        table_type_print = 'op2_results.' + table_type if '.' in table_type else table_type
        table = model.get_result(table_type)
        if table_type == 'superelement_tables':
            for key in table:
                msg.append(f'{table_type_print}[{key}]\n')
            continue

        try:
            sorted_tables = sorted(table.items(), key=_compare)
        except AttributeError:
            log.warning(f'table_type={table_type}; type(table)={type(table)}')
            raise

        for isubcase, subcase in sorted_tables:
            class_name = subcase.__class__.__name__
            if class_name in no_data_classes:
                msg.append('%s[%r]\n' % (table_type_print, isubcase))
            elif hasattr(subcase, 'data'):
                #data = subcase.data
                #shape = [int(i) for i in subcase.data.shape]
                #headers = subcase.get_headers()
                #headers_str = str(', '.join(headers))
                #msg.append('%s[%s]; %s; %s; [%s]\n' % (
                #table_type, isubcase, class_name, shape, headers_str))
                msg.append('%s[%s]\n' % (table_type_print, isubcase))
            elif table_type == 'params':  #  TODO: remove
                msgi = str(subcase)
            elif hasattr(subcase, 'get_stats'):
                msgi = '%s[%s] # unvectorized\n' % (table_type_print, isubcase)
                msg.append(msgi)
            else:
                msgi = 'skipping %r %s[%s]\n' % (class_name, table_type_print, isubcase)
                msg.append(msgi)
                #raise RuntimeError(msgi)
    return msg

def _get_op2_results_stats_dict(obj: dict[Any, Any],
                                table_type: str, short: bool) -> list[str]:
    msg = []
    for key, obji in obj.items():
        if isinstance(obji, list):
            for iobj, objii in enumerate(obji):
                stats = objii.get_stats(short=short)
                msg.extend(f'op2_results.{table_type}[{key}][{iobj}]: ' + stats)
        else:
            stats = obji.get_stats(short=short)
            msg.extend(f'op2_results.{table_type}[{key}]: ' + stats)
    return msg


def _get_op2_stats_full(model: OP2, table_types: list[str],
                        log: SimpleLogger) -> list[str]:
    """helper for get_op2_stats(...)"""
    msg = []
    handled_previously = ['params', 'grid_point_weight', 'psds']
    for table_type in table_types:
        table = model.get_result(table_type)
        if table_type in handled_previously:
            continue

        skip_results = ('gpdt', 'bgpdt', 'eqexin', 'monitor1', 'monitor3', 'cstm')
        if table_type in skip_results or table_type.startswith('responses.'):
            obj = model.get_result(table_type)
            if obj is None:
                continue
            elif isinstance(obj, dict):
                msg.extend(_get_op2_results_stats_dict(obj, table_type, short=False))
                continue

            obj = cast(ScalarObject, obj)
            stats = obj.get_stats(short=False)
            msg.extend(f'op2_results.{table_type}: ' + stats)  # TODO: a hack...not quite right...
            continue

        table_type_print = 'op2_results.' + table_type if '.' in table_type else table_type
        op2_result_tables = [
            #'vg_vf_response',
            'superelement_tables',
        ]
        if table_type in op2_result_tables:
            for key in table:
                msg.append(f'{table_type_print}[{key}]\n')
            continue

        try:
            for isubcase, subcase in sorted(table.items(), key=_compare):
                class_name = subcase.__class__.__name__
                if hasattr(subcase, 'get_stats'):
                    try:
                        stats = subcase.get_stats() # short=short
                    except Exception:
                        msgi = 'errored reading %s %s[%s]\n\n' % (
                            class_name, table_type_print, isubcase)
                        msg.append(msgi)
                        raise
                    else:
                        msg.append(f'{table_type_print}[{isubcase}]\n')
                        msg.extend(stats)
                        msg.append('\n')
                else:
                    msgi = 'skipping %s %s[%s]\n\n' % (class_name, table_type_print, isubcase)
                    msg.append(msgi)
                    raise RuntimeError(msgi)
        except Exception:
            log.warning(f'table_type={table_type}; type(table)={type(table)}')
            log.warning(str(table))
            raise
    return msg


def write_params_stats(params: dict[str, Any]) -> list[str]:
    """helper for get_op2_stats(...)"""
    if not params:
        return []
    msg = ['params:\n']
    iparam = 0
    for key, param in sorted(params.items()):
        if len(param.values) == 1:
            msg.append(f'  {key} = {param.values[0]!r}\n')
        else:
            msg.append(f'  {key} = {param.values}\n')
        iparam += 1
        if iparam > 10:
            msg.append(f'  ...\n')
            break
    return msg


def _compare(key_value):
    """helper for get_op2_stats(...)"""
    key = key_value[0]
    if isinstance(key, COMPARE_KEYS):
        return key
    #print('key=%s type=%s' % (key, type(key)))
    return key[0]
