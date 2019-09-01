"""
Defines utilities for parsing include files:
 - get_include_filename(card_lines, include_dir='', is_windows=None)

"""
import os
import ntpath
import posixpath
from typing import List, Tuple, Optional, Any
from pathlib import PurePosixPath, PureWindowsPath
from pyNastran.bdf.errors import EnvironmentVariableError

IS_WINDOWS = 'nt' in os.name
#is_linux = 'posix' in os.name
#is_mac = 'darwin' in os.name


def get_include_filename(card_lines: List[str], include_dir: str='',
                         is_windows: Optional[bool]=None) -> str:
    """
    Parses an INCLUDE file split into multiple lines (as a list).

    Parameters
    ----------
    card_lines : List[str]
        the list of lines in the include card (all the lines!)
    include_dir : str; default=''
        the include directory

    Returns
    -------
    filename : str
        the INCLUDE filename

    """
    if is_windows is None:
        is_windows = IS_WINDOWS

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
        raise SyntaxError('INCLUDE file is empty...card_lines=%s\n'
                          'there is a $ sign in the INCLUDE card' % card_lines)

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
    tokens = split_filename_into_tokens(include_dir, filename, is_windows)
    filename = str(tokens)
    return filename


def split_filename_into_tokens(include_dir: str, filename: str,
                               is_windows: bool) -> Any:
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

        # fails if the comment has stripped out the file (e.g., "INCLUDE '$ENV/model.bdf'")
        pth0 = pth[0]

        # Linux style paths
        # /work/model.bdf
        if len(pth0) == 1 and pth0[0] == '\\':
            # utterly breaks os.path.join
            raise SyntaxError('filename=%r cannot start with / on Windows' % filename)
    else:
        inc = PurePosixPath(include_dir)
        pth = PurePosixPath(filename).parts

        # fails if the comment has stripped out the file (e.g., "INCLUDE '$ENV/model.bdf'")
        pth0 = pth[0]
        if len(pth0) >= 2 and pth0[:2] == r'\\':
            # Windows network paths
            # \\nas3\work\model.bdf
            raise SyntaxError("filename=%r cannot start with \\\\ on Linux" % filename)

    pth2 = split_tokens(pth, is_windows)
    if is_windows:
        pth3 = ntpath.join(*pth2)
    else:
        pth3 = posixpath.join(*pth2)

    pth_out = inc / pth3
    return pth_out

def split_tokens(tokens, is_windows):
    # type: (Tuple[str], bool) -> List[str]
    """converts a series of path tokens into a joinable path"""
    tokens2 = [] # type: List[str]
    is_mac_linux = not is_windows
    for itoken, token in enumerate(tokens):
        # this is technically legal...
        #   INCLUDE '/testdir/dir1/dir2/*/myfile.dat'
        assert '*' not in token, '* in path not supported; tokens=%s' % tokens
        if is_windows:
            assert '$' not in token, '$ in path not supported; tokens=%s' % tokens
        else:
            assert '%' not in token, '%% in path not supported; tokens=%s' % tokens

        if itoken == 0 and is_mac_linux and ':' in token:
            ## no C:/dir/model.bdf on linux/mac
            #raise SyntaxError('token cannot include colons (:); token=%r; tokens=%s' % (
                #token, str(tokens)))

            # this has an environment variable or a drive letter
            #print(token)
            stokens = token.split(':')
            if len(stokens) != 2:
                msg = "len(stokens)=%s must be 2; stokens=%s" % (len(stokens), stokens)
                raise SyntaxError(msg)
            if len(stokens[0]) == 1:
                if len(stokens[1]) not in [0, 1]:
                    raise SyntaxError('tokens=%r token=%r stokens=%s stoken[1]=%r len=%r' % (
                        tokens, token, stokens, stokens[1], len(stokens[1])))

            if len(stokens[0]) < 2:
                raise SyntaxError('token cannot include colons (:); token=%r; tokens=%s' % (
                    token, str(tokens)))
            # variables in Windows are not case sensitive; not handled?
            token0 = stokens[0]
            if is_windows:
                assert '$' not in stokens[0], token0
                assert '%' not in stokens[0], token0

                if '%' in token0:
                    assert token0[0] == '%', token0
                    assert token0[-1] == '%', token0
                    token0 = '%' + token0 + '%'
                else:
                    token0 = '$' + token0

                #tokeni = os.path.expandvars('$' + stokens[0])
                tokeni = os.path.expandvars(token0)
                if '$' in tokeni:
                    raise SyntaxError('tokeni=%r has a $ in it after expanding (token0=%r)...\n'
                                      'tokens=%s stokens=%s' % (tokeni, token0, tokens, stokens))

                tokensi = PureWindowsPath(tokeni).parts
            else:
                if '$' in token0:
                    assert token0[0] == '$', token0
                else:
                    token0 = '$' + token0
                assert '%' not in stokens[0], token0
                tokeni = os.path.expandvars(token0)
                tokensi = PurePosixPath(tokeni).parts

            tokens2.extend(tokensi)
            tokens2.append(stokens[1])

        elif ':' in token:
            # Windows

            # this has an environment variable or a drive letter
            stokens = token.split(':')

            if len(stokens[0]) == 1:
                if len(stokens[1]) not in [0, 1]:
                    raise SyntaxError('tokens=%r token=%r stokens=%s stoken[1]=%r len=%r' % (
                        tokens, token, stokens, stokens[1], len(stokens[1])))
                # drive letter
                if itoken != 0:
                    raise SyntaxError('the drive letter is in the wrong place; '
                                      'itoken=%s; token=%r; stoken=%s; tokens=%s' % (
                                          itoken, token, stokens, tokens))
                tokens2.append(token)
            else:

                # variables in Windows are not case sensitive; not handled?
                environment_vars_to_expand = stokens[:-1]
                if len(environment_vars_to_expand) != 1:
                    raise SyntaxError(
                        'Only 1 environment variable can be expanded; '
                        'environment_vars_to_expand = %r' % environment_vars_to_expand)
                for env_var in environment_vars_to_expand:
                    if env_var.strip('$ %') not in os.environ:
                        environment_variables = list(os.environ.keys())
                        environment_variables.sort()
                        raise EnvironmentVariableError(
                            f"Can't find environment variable={repr(env_var)}"
                            f'\nenviron={environment_variables}\n'
                            f'which is required for {repr(tokens)}')

                    env_vari = os.path.expandvars('$' + env_var.strip('%'))
                    if '$' in env_vari:
                        raise SyntaxError('env_vari=%r has a $ in it after expanding (token0=%r)...\n'
                                          'tokens=%s stokens=%s' % (env_vari, env_var, tokens, stokens))
                    if is_windows:
                        tokensi = PureWindowsPath(env_vari).parts
                    else:
                        tokensi = PurePosixPath(env_vari).parts
                    tokens2.extend(tokensi)
                tokens2.append(stokens[-1])
        else:
            # standard
            tokens2.append(token)

    return tokens2
