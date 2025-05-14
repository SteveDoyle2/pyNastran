"""
Defines utilities for parsing include files:
 - get_include_filename(card_lines, include_dir='', is_windows=None)

"""
import os
import ntpath
import posixpath
from typing import Optional
from pathlib import Path, PurePosixPath, PureWindowsPath
from pyNastran.utils import print_bad_path, PathLike
from pyNastran.bdf.errors import EnvironmentVariableError

IS_WINDOWS = 'nt' in os.name
#is_linux = 'posix' in os.name
#is_mac = 'darwin' in os.name


def get_include_filename(log,
                         include_lines: list[str],
                         include_dirs: list[str],
                         replace_includes: dict[str, str],
                         source_filename: str='',
                         is_windows: Optional[bool]=None,
                         debug: bool=False,
                         write_env_on_error: bool=False) -> str:
    """
    Parses an INCLUDE file split into multiple lines (as a list).

    Parameters
    ----------
    include_lines : list[str]
        the list of lines in the include card (all the lines!)
    include_dirs : list[str]
        the include directory
    is_windows: Optional[bool]
        None: set it dynamically
        bool:
            Windows/Linux/Mac have different enviornment variable forms
            (%PATH% for Windows and $PATH for Linux/Mac). There are
            also different characters that are allowed in paths.

    Returns
    -------
    filename : str
        the INCLUDE filename

    """
    for line in include_lines:
        if len(line) > 72:
            msg = '\n - '.join(include_lines)
            log.warning(f'INCLUDE line={line!r} is too long (n={len(line)} for:\n - {msg}')
    #print(f'card_lines={card_lines};\nsource_filename={source_filename!r}')
    if not isinstance(include_dirs, list):
        assert isinstance(include_dirs, PathLike), include_dirs
        #if str(include_dirs) == '':
            #include_dirs = os.getcwd()
        include_dirs = [include_dirs]
    if is_windows is None:
        is_windows = IS_WINDOWS

    filename_raw = parse_include_lines(include_lines)
    if filename_raw in replace_includes:
        filename_raw = replace_includes[filename_raw]
        if len(filename_raw) == 0:
            return filename_raw

    ninclude_dirs = len(include_dirs)
    for include_dir in include_dirs:
        #print(f'filename_raw = {filename_raw}')
        #print(f'dir = {include_dir}')
        tokens = split_filename_into_tokens(
            include_dir, filename_raw, is_windows, debug=debug)
        filename = str(tokens)
        if os.path.exists(filename):
            if debug:
                print(f'found {filename}')
            break
        elif debug:
            print(f'could not find {filename}')
            print(print_bad_path(filename))
        if ninclude_dirs == 1:
            break
    else:
    #if 1:
        msg = f'Could not find INCLUDE line:\n{include_lines}\n'
        msg += f'  filename: {os.path.abspath(filename_raw)}\n'
        if source_filename:
            msg += f'  source file: {os.path.abspath(source_filename)}\n'
        msg += f'  include_dirs:\n - ' + '\n - '.join(repr(val) for val in include_dirs) + '\n'
        if write_env_on_error:
            msg += '  environment:'
            skip_keys = [
                'LESSOPEN', 'LOGNAME', 'LS_COLORS',
                'MAIL', 'NAME', 'XDG_SESSION_CLASS', 'XDG_DATA_DIRS',
                'XDG_RUNTIME_DIR', 'XDG_SESSION_ID', 'XDG_SESSION_TYPE', 'DISPLAY',
                'LANG', 'HOSTTYPE', 'MOTD_SHOWN', 'DEBUGINFOD_URLS', 'DISPLAY', 'USER',
                'TERM', 'WSL_DISTRO_NAME', 'WSL_INTEROP', 'WT_PROFILE_ID', 'WT_SESSION',
                'WSL2_GUI_APPS_ENABLED', 'WAYLAND_DISPLAY', 'PULSE_SERVER',
            ]
            for key, value in sorted(os.environ.items()):
                if key in skip_keys:
                    continue
                msg += f'  {key}: {value!r}\n'
            #else:
                #msg += f'  environment_keys: {list(key for key in os.environ)}'
        raise IOError(msg)
    return filename

def parse_include_lines(card_lines: list[str]) -> str:
    """handles splitting out the INCLUDE lines"""
    card_lines2 = []
    for line in card_lines:
        line = line.strip('\t\r\n ')
        card_lines2.append(line)

    # combine the lines
    # -----------------
    # initial:
    # INCLUDE '/dir123
    #          /dir456
    #          /dir789/filename.dat'
    #
    # final (not the extra single quotes):
    # "'/dir123/dir456/dir789/filename.dat'"
    card_lines2[0] = card_lines2[0][7:].strip()  # remove the INCLUDE
    filename = ''.join(card_lines2)

    # drop the single quotes:
    #     "'path1/path2/model.inc'" to "path1/path2/model.inc"
    filename = filename.strip('"').strip("'")

    if len(filename.strip()) == 0:
        raise SyntaxError(f'INCLUDE file is empty...card_lines={card_lines}\n'
                          'there is a $ sign in the INCLUDE card')

    # not handled...
    #    include '/mydir' /level1 /level2/ 'myfile.x'
    #
    # -> /proj/dept123/sect 456/joe/flange.bdf

    # probably not handled...
    # include c:\project,
    # $ A comment line
    # '\Data Files' \subdir\thisfile
    #
    # -> C:\PROJECT\Data Files\SUBDIR\THISFILE
    return filename

def split_filename_into_tokens(include_dir: str, filename: str,
                               is_windows: bool,
                               debug: bool=False) -> Path:
    r"""
    Tokens are the individual components of paths

    Invalid Linux Tokens
    '\0' (NUL)

    Invalid Windows Tokens
    < (less than)
    > (greater than)
    : (colon - sometimes works, but is actually NTFS Alternate Data Streams)
    " (double quote)
    / (forward slash)
    \ (backslash)
    | (vertical bar or pipe)
    ? (question mark)
    * (asterisk)
    All control codes (<= 31)
    """
    if is_windows:
        inc = PureWindowsPath(include_dir)
        pth = PureWindowsPath(filename).parts

        # fails if the comment has stripped out the file
        # (e.g., "INCLUDE '$ENV/model.bdf'")
        pth0 = pth[0]

        # Linux style paths
        # /work/model.bdf
        if len(pth0) == 1 and pth0[0] == '\\':
            # utterly breaks os.path.join
            raise SyntaxError(f'filename={filename!r} cannot start with / on Windows')
    else:
        inc = PurePosixPath(include_dir)
        pth = PurePosixPath(filename).parts

        # fails if the comment has stripped out the file (e.g., "INCLUDE '$ENV/model.bdf'")
        pth0 = pth[0]
        if len(pth0) >= 2 and pth0[:2] == r'\\':
            # Windows network paths
            # \\nas3\work\model.bdf
            raise SyntaxError(f'filename={filename!r} cannot start with \\\\ on Linux')

    pth2 = split_tokens(pth, is_windows, debug=debug)
    if is_windows:
        pth3 = ntpath.join(*pth2)
    else:
        pth3 = posixpath.join(*pth2)
    pth_out = inc / pth3
    return pth_out

def split_tokens(tokens: tuple[str], is_windows: bool,
                 debug: bool=False) -> list[str]:
    """converts a series of path tokens into a joinable path"""
    tokens2: list[str] = []
    is_mac_linux = not is_windows
    for itoken, token in enumerate(tokens):
        # this is technically legal...
        #   INCLUDE '/testdir/dir1/dir2/*/myfile.dat'
        assert '*' not in token, f'* in path not supported; tokens={tokens}'
        if is_windows:
            assert '$' not in token, f'$ in path not supported; tokens={tokens}'
        else:
            assert '%' not in token, f'%% in path not supported; tokens={tokens}'

        if itoken == 0 and is_mac_linux and ':' in token:
            tokensi, stokens = split_drive_token(
                itoken, tokens, is_windows, debug=debug)
            tokens2.extend(tokensi)
            tokens2.append(stokens[1])

        elif ':' in token:
            # Logical symbols provide you with a way of specifying file
            # locations with a convenient shorthand. This feature also allows
            # input files containing filename specifications to be moved
            # between computers without requiring modifications to the input
            # files. Only the logical symbol definitions that specify actual
            # file locations need to be modified.

            # Windows
            # this has an environment variable or a drive letter
            stokens = token.split(':')

            if len(stokens[0]) == 1:
                # this is a drive letter
                if len(stokens[1]) not in [0, 1]:
                    raise SyntaxError(
                        f'tokens={tokens!r} token={token!r} stokens={stokens} '
                        f'stokens[1]={stokens[1]!r}; len={len(stokens[1]):d}')
                # drive letter
                if itoken != 0:
                    raise SyntaxError('the drive letter is in the wrong place; '
                                      f'itoken={itoken:d}; token={token!r}; '
                                      f'stokend={stokens}; tokens={tokens}')
                tokens2.append(token)
            else:
                # variables in Windows are not case sensitive; not handled?
                environment_vars_to_expand = stokens[:-1]
                if len(environment_vars_to_expand) != 1:
                    raise SyntaxError(
                        'Only 1 environment variable can be expanded; '
                        f'environment_vars_to_expand = {environment_vars_to_expand!r}')
                for env_var in environment_vars_to_expand:
                    if env_var.strip('$ %') not in os.environ:
                        environment_variables = list(os.environ.keys())
                        environment_variables.sort()
                        raise EnvironmentVariableError(
                            f"Can't find environment variable={repr(env_var)}\n"
                            f'environ={environment_variables}\n'
                            f'which is required for {repr(tokens)}')

                    env_vari = os.path.expandvars('$' + env_var.strip('%'))
                    if debug:
                        print(f'expanded env {env_var!r} -> {env_vari!r}')
                    if '$' in env_vari:
                        raise SyntaxError(f'env_vari={env_vari!r} has a $ in it after expanding (token0={env_var!r})...\n'
                                          f'tokens={tokens} stokens={stokens}')
                    if is_windows:
                        tokensi = PureWindowsPath(env_vari).parts
                    else:
                        tokensi = PurePosixPath(env_vari).parts
                    tokens2.extend(tokensi)
                tokens2.append(stokens[-1])
        else:
            # standard
            tokens2.append(token)
    if debug:
        print(f'  split_tokens(is_windows={is_windows}):')
        print(f'    tokens:  {tokens}')
        print(f'    tokens2: {tokens2}')
    return tokens2


def split_drive_token(itoken: int, tokens: tuple[str, ...],
                      is_windows: bool,
                      debug: bool=False) -> tuple[tuple[str, ...], list[str]]:
    """
    splits the drive off:
       C:/work/file.bdf (Windows)
       ENVVAR:file.bdf  (Windows/Linux)

    """
    ## no C:/dir/model.bdf on linux/mac
    # raise SyntaxError('token cannot include colons (:); token=%r; tokens=%s' % (
    # token, str(tokens)))

    # this has an environment variable or a drive letter
    # print(token)
    assert isinstance(tokens, tuple), tokens
    token = tokens[itoken]
    stokens = token.split(':')
    nstokens = len(stokens)
    token0 = stokens[0]
    if nstokens != 2:
        msg = f'len(stokens)={nstokens:d} must be 2; stokens={stokens}'
        raise SyntaxError(msg)
    if len(token0) == 1:
        if len(stokens[1]) not in [0, 1]:
            msg = (f'tokens={tokens} tokens[{itoken:d}={token!r} stokens={stokens} '
                   f'stoken[1]={stokens[1]!r}; len={len(stokens[1]):d}')
            raise SyntaxError(msg)

    if len(token0) < 2:
        raise SyntaxError('token cannot include colons (:); '
                          f'token={token!r}; tokens={tokens}')

    # variables in Windows are not case sensitive; not handled?
    if is_windows:
        assert '$' not in token0, token0
        assert '%' not in token0, token0

        if '%' in token0:
            assert token0[0] == '%', token0
            assert token0[-1] == '%', token0
            token0 = '%' + token0 + '%'
        else:
            token0 = '$' + token0

        # tokeni = os.path.expandvars('$' + token0)
        tokeni = os.path.expandvars(token0)
        if debug:
            print(f'expanded env {token0!r} -> {tokeni!r}')
        if '$' in tokeni:
            raise SyntaxError(f'tokeni={tokeni!r} has a $ in it after '
                              f'expanding (token0={token0!r})...\n'
                              f'tokens={tokens} stokens={stokens}')
        tokensi = PureWindowsPath(tokeni).parts
    else:
        if '$' in token0:
            assert token0[0] == '$', token0
        else:
            token0 = '$' + token0
        assert '%' not in stokens[0], token0
        tokeni = os.path.expandvars(token0)
        if debug:
            print(f'expanded env {token0!r} -> {tokeni!r}')
        tokensi = PurePosixPath(tokeni).parts
    assert isinstance(tokensi, tuple), type(tokensi)
    assert isinstance(stokens, list), type(stokens)
    return tokensi, stokens
