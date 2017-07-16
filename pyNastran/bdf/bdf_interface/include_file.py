"""
Defines utilities for parsing include files:
 - get_include_filename(card_lines, include_dir='', is_windows=None)
"""
from __future__ import print_function, unicode_literals
import os
import ntpath
import posixpath
from six import PY2
if PY2:
    from pathlib2 import PurePosixPath, PureWindowsPath
else:
    from pathlib import PurePosixPath, PureWindowsPath

IS_WINDOWS = 'nt' in os.name
#is_linux = 'posix' in os.name
#is_mac = 'darwin' in os.name

#class PurePosixPath(object):
    #def __init__(self, path):
        #self.path = path

#class PureWindowsPath(object):
    #def __init__(self, path):
        ##self.path = path
        #if isinstance(path, list):
            #self._parts = path
        #else:
            #self._parts = self._split_path_by_slash(path)

    #@classmethod
    #def path_from_parts(parts):
        #pth = PureWindowsPath(parts)
        #print(pth)
        #return pth

    #def _split_path_by_slash(self, path):
        #print(path)
        #parts = []
        #partsi = os.path.splitunc(path)
        #if not partsi[0]:
            #partsi = os.path.splitdrive(path)
        #parts.append(partsi[0])
        #pth = os.path.normcase(partsi[1])
        ##pth2 = os.path.split(pth)
        #parts.append(pth)
        ##print(pth2)
        ##print(parts)
        ##print(partsi)
        #return parts

    #@property
    #def parts(self):
        #return self._parts

    ##def __truediv__ (self, pth):
        ##print(pth)
        ##print(self._parts)
    ##def __rdiv__(self, pth):
        ##print(pth)
        ##print(self._parts)
    #def __div__(self, right_path):
        #print('div...')
        #print('right_path =', right_path)
        #print(self._parts)
        #parts = os.path.join(self._parts,  right_path)
        #return PureWindowsPath.path_from_parts(parts)


    #def __repr__(self):
        #print("repr parts = %s" % self._parts)
        #return '%r' % str(os.path.join(*self._parts))

def get_include_filename(card_lines, include_dir='', is_windows=None):
    # type: (List[str], str) -> str
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
    print("filename = %r" % filename)

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
    print("tokens = ", tokens)
    filename = str(tokens)
    #print("filename = ", filename)
    return filename

#def _split_to_tokens(pth, is_windows):
    #if is_windows:
        #tokens = pth.split(ntpath.sep)
    #else:
        #tokens = pth.split(posixpath.sep)
    #return tokens


def split_filename_into_tokens(include_dir, filename, is_windows):
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
    #print('--------------')
    #print('  include_dir_in = %r' % include_dir)
    #print('  filename_in = %r' % filename)
    if is_windows:
        inc = PureWindowsPath(include_dir)#.parts
        pth = PureWindowsPath(filename).parts
        #print("pth[0] = %r" % pth[0])
        if len(pth[0]) == 1 and pth[0][0] == '\\':
            # utterly breaks os.path.join
            raise SyntaxError('filename=%r cannot start with / on Windows' % filename)
    else:
        inc = PurePosixPath(include_dir)#.parts
        pth = PurePosixPath(filename).parts
        pth0 = pth[0]
        #print("pth[:1] = %r" % pth0, '%r' % pth0[:2])
        #print('pth1 = %r' % pth0[:2], len(pth0[:2]))
        if len(pth0) >= 2 and pth0[:2] == r'\\':
            # network path
            raise SyntaxError("filename=%r cannot start with \\\\ on Linux" % filename)
    print('inc =', inc)
    print('pth.parts =', pth)

    #inc2 = split_tokens(inc, is_windows)
    pth2 = split_tokens(pth, is_windows)

    #print('inc2 =', inc2)
    print('pth2 (tokens) =', pth2)

    if is_windows:
        pth3 = ntpath.join(*pth2)
    else:
        pth3 = posixpath.join(*pth2)

    #print('pth3 =', pth3)
    pth_out = inc / pth3
    #print('pth_out =', pth_out)
    return pth_out

def split_tokens(tokens, is_windows):
    """converts a series of path tokens into a joinable path"""
    tokens2 = []
    is_mac_linux = not is_windows
    for itoken, token in enumerate(tokens):
        print('tokens[%i] = %r' % (itoken, token))

        # this is technically legal...
        #   INCLUDE '/testdir/dir1/dir2/*/myfile.dat'
        assert '*' not in token, '* in path not supported; tokens=%s' % tokens
        assert '$' not in token, '$ in path not supported; tokens=%s' % tokens
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
            else:
                # variables in Windows are not case sensitive; not handled?
                tokeni = os.path.expandvars('$' + stokens[0])
                if is_windows:
                    tokensi = PureWindowsPath(tokeni).parts
                else:
                    tokensi = PurePosixPath(tokeni).parts

                tokens2.extend(tokensi)
                tokens2.append(stokens[1])

        elif ':' in token:
            # this has an environment variable or a drive letter
            stokens = token.split(':')
            print('  stokens = %r' % stokens)
            if len(stokens) != 2:
                msg = 'stokens=%s h', stokens
            #stokens2 = []
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
                    env_vari = os.path.expandvars('$' + env_var)
                    if is_windows:
                        tokensi = PureWindowsPath(env_vari).parts
                    else:
                        tokensi = PurePosixPath(env_vari).parts
                    tokens2.extend(tokensi)
                    print("    expanding env_var=%r to %r" % (env_var, tokensi))
                tokens2.append(stokens[-1])
        else:
            # standard
            tokens2.append(token)

    print("tokens2 =", tokens2)
    #print(os.path.join(*tokens2))
    return tokens2
