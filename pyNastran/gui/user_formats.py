"""
User-defined pyNastranGUI formats.

A user can add (or override) a format in ``CLASS_MAP`` without editing
pyNastran by listing it in ``~/pynastrangui_formats.json``::

    # pyNastranGUI user-defined formats
    # '#' starts a comment (outside of a quoted string)
    {
        # name used by "pyNastranGUI -f cart3d model.tri"
        "cart3d": {
            # a *.py file or a dotted module name (pyNastran.converters.cart3d.cart3d_io)
            "path": "C:/work/code/pyNastran/pyNastran/converters/cart3d/cart3d_io.py",
            "class": "Cart3dIO",
            # optional; used to guess the format from the file extension
            "extensions": [".tri", ".triq"]
        }
    }

The class must follow the standard pyNastranGUI format interface
(``__init__(self, gui)`` and
``get_<name>_wildcard_geometry_results_functions()``).

defines:
 - USER_FORMATS_FILENAME
 - strip_json_comments(text)
 - read_user_formats(filename=None, log=None)
 - write_user_formats(user_formats, filename=None)
 - load_user_format_classes(user_formats, log)
 - load_plugin_class(module_name, plugin_file, class_name, log)

"""
import os
import sys
import json
import importlib
import importlib.util
from typing import Optional, Any

from cpylog import SimpleLogger

USER_FORMATS_FILENAME = os.path.join(os.path.expanduser('~'), 'pynastrangui_formats.json')
UserFormats = dict[str, dict[str, Any]]


def strip_json_comments(text: str) -> str:
    """
    Removes ``#`` comments from a JSON string.

    Parameters
    ----------
    text : str
        the JSON text with comments

    Returns
    -------
    text2 : str
        the JSON text without comments

    Notes
    -----
    A ``#`` inside a double-quoted string (e.g., a path with a ``#``)
    is not a comment.  Line breaks are kept, so json error line numbers
    still match the original file.

    """
    out: list[str] = []
    in_string = False
    is_escaped = False
    in_comment = False
    for char in text:
        if in_comment:
            if char == '\n':
                in_comment = False
                out.append(char)
            continue
        if in_string:
            if is_escaped:
                is_escaped = False
            elif char == '\\':
                is_escaped = True
            elif char == '"':
                in_string = False
        elif char == '"':
            in_string = True
        elif char == '#':
            in_comment = True
            continue
        out.append(char)
    return ''.join(out)


def _format_error(filename: str, name: str, msg: str) -> str:
    return f'{filename}: user format {name!r} {msg}'


def read_user_formats(filename: Optional[str]=None,
                      log: Optional[SimpleLogger]=None) -> UserFormats:
    """
    Reads the user-defined formats file.

    Parameters
    ----------
    filename : str; default=None -> USER_FORMATS_FILENAME
        the json file (``#`` comments are allowed)
    log : SimpleLogger; default=None
        logs invalid entries; None -> invalid entries are silently skipped

    Returns
    -------
    user_formats : dict[str, dict[str, Any]]
        name -> {'path': str, 'class': str, 'extensions': list[str]}
        The name is lowercased.  A missing file returns {}.

    Notes
    -----
    A bad file/entry is logged and skipped rather than raised, so a typo
    in the json doesn't keep the GUI from opening.

    """
    if filename is None:
        filename = USER_FORMATS_FILENAME
    if not os.path.exists(filename):
        return {}

    with open(filename, 'r', encoding='utf-8') as json_file:
        text = json_file.read()
    text = strip_json_comments(text)
    if not text.strip():
        # everything is commented out
        return {}
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        # don't crash the gui on a typo; tell the user and move on
        if log is not None:
            log.error(f'{filename}: invalid json; {error}\n'
                      'Note: backslashes in paths must be doubled (or use /)')
        return {}

    if not isinstance(data, dict):
        if log is not None:
            log.error(f'{filename}: expected a dictionary of formats; got {type(data).__name__}')
        return {}

    user_formats: UserFormats = {}
    for name, value in data.items():
        name_lower = name.strip().lower()
        if not isinstance(value, dict):
            if log is not None:
                log.error(_format_error(filename, name, 'must be a dictionary with "path" and "class"'))
            continue
        path = value.get('path', '')
        class_name = value.get('class', '')
        extensions = value.get('extensions', [])
        if isinstance(extensions, str):
            extensions = [extensions]
        if not isinstance(path, str) or not path.strip():
            if log is not None:
                log.error(_format_error(filename, name, 'is missing "path"'))
            continue
        if not isinstance(class_name, str) or not class_name.strip():
            if log is not None:
                log.error(_format_error(filename, name, 'is missing "class"'))
            continue
        user_formats[name_lower] = {
            'path': path.strip(),
            'class': class_name.strip(),
            'extensions': _clean_extensions(extensions),
        }
    return user_formats


def _clean_extensions(extensions: list[str]) -> list[str]:
    """'tri' -> '.tri'; lowercased"""
    out = []
    for ext in extensions:
        ext = str(ext).strip().lower()
        if not ext:
            continue
        if not ext.startswith('.'):
            ext = '.' + ext
        out.append(ext)
    return out


def write_user_formats(user_formats: UserFormats,
                       filename: Optional[str]=None) -> None:
    """
    Writes the user-defined formats file.

    Parameters
    ----------
    user_formats : dict[str, dict[str, Any]]
        name -> {'path': str, 'class': str, 'extensions': list[str]}
    filename : str; default=None -> USER_FORMATS_FILENAME
        the json file

    Notes
    -----
    Comments in an existing file are not preserved (a header is written).

    """
    if filename is None:
        filename = USER_FORMATS_FILENAME
    data = {}
    for name, value in user_formats.items():
        data[name.strip().lower()] = {
            'path': value['path'],
            'class': value['class'],
            'extensions': _clean_extensions(value.get('extensions', [])),
        }
    header = (
        '# pyNastranGUI user-defined formats\n'
        '# "#" starts a comment; a user format overwrites a built-in format with the same name\n'
        '# name: {"path": <*.py file or dotted.module>, "class": <class name>, "extensions": [".ext"]}\n'
    )
    with open(filename, 'w', encoding='utf-8') as json_file:
        json_file.write(header)
        json.dump(data, json_file, indent=4)
        json_file.write('\n')


def _is_file_path(path: str) -> bool:
    """a *.py file (vs. a dotted module name)"""
    return path.lower().endswith('.py') or '/' in path or '\\' in path


def load_user_format_classes(user_formats: UserFormats,
                             log: SimpleLogger) -> dict[str, type]:
    """
    Imports the user-defined format classes.

    Parameters
    ----------
    user_formats : dict[str, dict[str, Any]]
        the output from ``read_user_formats``
    log : SimpleLogger
        logs failed imports

    Returns
    -------
    class_map : dict[str, type]
        name -> class; failed imports are skipped

    Notes
    -----
    A ``*.py`` file is loaded as module ``pynastrangui_user_format_<name>``,
    so it must use absolute imports (no ``from .foo import bar``).

    """
    class_map: dict[str, type] = {}
    for name, value in user_formats.items():
        path = value['path']
        class_name = value['class']
        if _is_file_path(path):
            if not os.path.exists(path):
                log.error(f'user format {name!r}: {path} does not exist')
                continue
            module_name = f'pynastrangui_user_format_{name}'
            cls = load_plugin_class(module_name, path, class_name, log)
        else:
            cls = _load_module_class(path, class_name, log)
        if cls is None:
            continue
        geom_results_funcs = f'get_{name}_wildcard_geometry_results_functions'
        if not hasattr(cls, geom_results_funcs):
            log.error(f'user format {name!r}: class {class_name} is missing '
                      f'{geom_results_funcs}()')
            continue
        class_map[name] = cls
    return class_map


def _load_module_class(module_name: str, class_name: str,
                       log: SimpleLogger) -> Optional[type]:
    """imports class_name from a dotted module (e.g., pyNastran.converters.cart3d.cart3d_io)"""
    try:
        module = importlib.import_module(module_name)
    except Exception as error:
        # any error in the user's module (ImportError, SyntaxError, ...)
        # shouldn't take down the gui
        log.warning(f'Could not import {module_name}')
        log.warning(str(error))
        return None
    try:
        return getattr(module, class_name)
    except AttributeError:
        log.warning(f"Class '{class_name}' not found in {module_name}")
        return None


def load_plugin_class(module_name: str, plugin_file: str, class_name: str,
                      log: SimpleLogger) -> Optional[type]:
    """Dynamically loads a Python file and returns a specific class from it.

    Compatible with Python 3.5-3.15+ (avoids deprecated load_module).
    """
    # 1. Generate the module metadata (spec)
    spec = importlib.util.spec_from_file_location(module_name, plugin_file)
    if spec is None or spec.loader is None:
        log.warning(f"Could not load spec for {plugin_file}")
        return None

    # 2. Create the empty module object
    module = importlib.util.module_from_spec(spec)

    # 3. Register it in sys.modules so internal imports work correctly
    sys.modules[module_name] = module

    try:
        # 4. Execute the file contents to populate the module
        spec.loader.exec_module(module)
    except Exception as error:
        # Clean up sys.modules if the plugin has a syntax/runtime error on load
        sys.modules.pop(module_name, None)
        log.warning(f"Could not load module for {plugin_file}")
        log.warning(str(error))
        return None

    # 5. Extract and return the requested class
    try:
        my_class = getattr(module, class_name)
    except AttributeError as error:
        log.warning(f"Class '{class_name}' not found in {plugin_file}")
        log.warning(str(error))
        return None
    return my_class
