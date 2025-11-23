import numpy as np
from typing import Optional, Callable


def get_setting(settings, setting_keys: list[str],
                setting_names: list[str],
                default: int | float | str | np.ndarray,
                auto_type: Optional[Callable]=None):
    """
    helper method for ``reapply_settings``

    does this, but for a variable number of input names, but one output name:
        screen_shape = settings.value("screen_shape", screen_shape_default)

    If the key is not defined, the default is used.
    """
    unused_set_name = setting_names[0]
    pull_name = None
    for key in setting_names:
        if key in setting_keys:
            pull_name = key
            break

    if pull_name is None:
        value = default
    else:
        try:
            value = settings.value(pull_name, default)
        except (TypeError, ValueError, RuntimeError):
            # print('couldnt load %s; using default' % pull_name)
            value = default

    if value is None and default is not None:
        #print('couldnt load %s; using default' % pull_name)
        value = default
    if default is not None:
        assert value is not None, pull_name

    value = autotype_value(value, auto_type)
    return value

def autotype_value(value, auto_type):
    """casts the values to the requested type"""
    if auto_type is None:
        return value

    if isinstance(value, np.ndarray):
        value = np.array([auto_type(valuei) for valuei in value])
    elif isinstance(value, list):
        value = [auto_type(valuei) for valuei in value]
    elif isinstance(value, tuple):
        value = [auto_type(valuei) for valuei in value]
        value = tuple(value)
    else:
        if auto_type == bool:
            if value == 'true':
                value = True
            elif value == 'false':
                value = False
            else:
                value = auto_type(value)
        else:
            value = auto_type(value)
    return value
