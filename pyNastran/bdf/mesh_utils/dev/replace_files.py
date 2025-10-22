import os
import sys
import shutil
import filecmp
import argparse
from itertools import count

import pyNastran
from pyNastran.utils import PathLike
from pyNastran.utils.dev import get_files_of_types


def replace_files(dirname_source: PathLike, dirname_target: PathLike,
                  style: str='revision', test: bool=True) -> None:
    """WARNING: does not handle fixing up model traces"""
    extensions = ['.bdf', '.dat', '.ecd', '.blk', '.inc', '.pch']
    _replace_files(dirname_source, dirname_target, extensions, style=style, test=test)


def _replace_files(dirname_source: PathLike, dirname_target: PathLike,
                   extensions: list[str],
                   style: 'revision',
                   test: bool=True) -> None:
    """
    Files in dirname_source are assumed to follow:
     - wing/wing.bdf
     - fuselage/fuselage_core.bdf
    Let's update wing.bdf in dirname_target to get:
     - wing/wing_r1.bdf
     - fuselage/fuselage_core.bdf
    Let's do a second update with a new wing.bdf and fuselage_core.bdf to get:
     - wing/wing_r2.bdf
     - fuselage/fuselage_core_r2.bdf
    
    Parameters
    ----------
    dirname_source: Path/str
        source directory
    dirname_target: Path/str
        target directory
    style: str; default='revision'
        Does nothing
        revision: ???
     test: bool; default=True
        should the files actually be updated?

    """
    assert style in ['revision'], f'Harry {style!r}.title() is not a valid style'
    run = not test

    # full path to files
    dirname_fnames_source = get_files_of_types(dirname_source, extensions, max_size_mb=0.)
    dirname_fnames_target = get_files_of_types(dirname_target, extensions, max_size_mb=0.)

    # make the files relative to their parent
    fnames_source = [os.path.relpath(fname, dirname_source) for fname in dirname_fnames_source]
    fnames_target = [os.path.relpath(fname, dirname_target) for fname in dirname_fnames_target]

    # simplify the files for easy comparision
    fnames_target_no_revs, revisions = remove_revisions(fnames_target)

    assert run is False

    for isource, fname_source, filename_source in zip(count(), fnames_source, dirname_fnames_source):

        if fname_source not in fnames_target_no_revs:
            # if the new revision-less file doesn't exist, copy the new file
            print(f'copy file {isource}: {fname_source} b/c its a new file')
            filename_target = os.path.join(dirname_target, fname_source)
            if run:
                asdf
                shutil.copyfile(filename_source, filename_target)
            continue

        # fname_source exists in the target list
        itarget = fnames_target_no_revs.index(fname_source)

        # get the current target file with a revision
        filename_target0 = dirname_fnames_target[itarget]
        fname_target0 = fnames_target[itarget]
        revision = revisions[itarget]

        # let's check to see if the files are the same
        is_same = filecmp.cmp(filename_source, filename_target0)
        if is_same:
            #print(f' -> {fname_source} is same')
            pass
        else:
            rev_new = revision + 1
            base, ext = os.path.splitext(fname_source)
            fname_target = f'{base}_r{rev_new}{ext}'
            print(f' -> {fname_source} is different; {fname_target0}; r{rev_new} -> {fname_target}')
            filename_target_rev_new = os.path.join(dirname_target, fname_source)
            if run:
                asdf
                shutil.copyfile(filename_source, filename_target_rev_new)


def remove_revisions(fnames: list[str]) -> tuple[list[str], list[int]]:
    """

    Parameters
    ----------
    fnames: list[str]
        paths to input files of the form:
        - wing.bdf
        - fuselage_core_r2.bdf

    Returns
    -------
    fnames_no_revs: list[str]
        paths to input files of the form:
        - wing.bdf
        - fuselage_core.bdf
    revisions: list[int]
        revision numbers; [0, 2] in this case

    """
    fnames_no_revs = []
    revisions = []
    for fname in fnames:
        base, ext = os.path.splitext(fname)
        if '_' not in base:
            fnames_no_revs.append(fname)
            revisions.append(0)
            continue
        base_no_rev, rev = base.rsplit('_', 1)
        if not rev.startswith('r'):
            fnames_no_revs.append(fname)
            revisions.append(0)
            continue

        # validate that it's a revison
        ## TODO: try-except not a revision...
        num = int(rev[1:])
        fname2 = base_no_rev + ext
        #print(f' - to update: {fname} to {fname2}')
        fnames_no_revs.append(fname2)
        revisions.append(num)
    return fnames_no_revs, revisions


def cmd_line_replace_files(argv=None, quiet: bool=False) -> None:
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [FILE] + argv[2:]  # ['run_jobs'] + sys.argv[2:]
    if not quiet:
        print(f'argv = {argv}')
    #print(f'argv = {argv}')
    parser = argparse.ArgumentParser(prog='run_jobs')
    parser.add_argument('dirname_source', help='path to source directory')
    parser.add_argument('dirname_target', help='path to source directory')
    parser.add_argument('--test', default=False, help='trial run')
    #parser.add_argument('--debug', action='store_true', help='more debugging')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                        version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    show_folder = not args.nofolder
    dirname_source = args.dirname_source
    dirname_target = args.dirname_target
    replace_files(dirname_source, dirname_target)
