from pathlib import Path
import unittest

import pyNastran
from pyNastran.bdf.mesh_utils.run_jobs import get_bdf_filenames_to_run, cmd_line_run_jobs
from pyNastran.bdf.mesh_utils.host_jobs import cmd_line_host_jobs
from pyNastran.utils.nastran_utils import _get_keywords_list

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'


class TestRunHostJobs(unittest.TestCase):
    def test_keyworks(self):
        keywords_strs = [
            'scr=yes', 'old=no', 'news=no', 'mem=16gb',
            'parallel=8', 'auth=123@host', 'endian=little',
            'bat=no', 'fake=cat',
        ]
        for keywords_str in keywords_strs:
            keywords_list = _get_keywords_list(keywords_str)
        keywords_list = _get_keywords_list(keywords_strs)

        keyword_dict = {
            'scr': 'yes',
            'old': 'no',
            'news': "no",
        }
        keywords_list = _get_keywords_list(keyword_dict)

    def test_host_jobs(self):
        str_model_path = str(MODEL_PATH)
        args = ['bdf', 'host_jobs', str_model_path, '--test', '--nmax', '2']
        nfiles = cmd_line_host_jobs(args, quiet=True)

    def test_run_jobs_path(self):
        str_model_path = str(MODEL_PATH)
        args = ['bdf', 'run_jobs', str_model_path, '--cleanup', '-r', '--test']
        nfiles = cmd_line_run_jobs(args)
        nfiles = cmd_line_run_jobs(args, quiet=True)
        assert nfiles >= 1, nfiles  # 105

    def _test_run_jobs_path2(self):
        """doesn't work remotely"""
        str_model_path = str(MODEL_PATH)
        extensions = ['.dat', '.bdf']
        bdf_files = get_bdf_filenames_to_run(MODEL_PATH, extensions, recursive=True)
        assert len(bdf_files) >= 10, len(bdf_files)  # 105

        #------------------------------
        bdf_files = get_bdf_filenames_to_run(MODEL_PATH, extensions, recursive=False)
        assert len(bdf_files) == 1, len(bdf_files)

    def test_run_jobs_str(self):
        extensions = ['.dat', '.bdf']
        str_model_path = str(MODEL_PATH)

        nfiles = cmd_line_run_jobs([
            'bdf', 'run_jobs', str_model_path,
            '--cleanup', '-r', '--test'], quiet=True)

        bdf_files = get_bdf_filenames_to_run(str_model_path, extensions, recursive=True)
        assert len(bdf_files) >= 10, len(bdf_files)  # 105

        # lists
        bdf_files = get_bdf_filenames_to_run([str_model_path], extensions, recursive=True)
        assert len(bdf_files) >= 10, len(bdf_files)  # 105

    def test_run_jobs_out_in(self):
        #extensions = ['.dat', '.bdf']
        str_model_path = str(MODEL_PATH)
        # out_filename = str(MODEL_PATH / 'run_files.out')
        out_filename = 'run_files.out'

        args_out = ['bdf', 'run_jobs', str_model_path, '-r', '--outfile', out_filename, '--test']
        args_in = ['bdf', 'run_jobs', str_model_path, '-r', '--infile', out_filename, '--test']
        nfiles = cmd_line_run_jobs(args_out, quiet=True)
        nfiles = cmd_line_run_jobs(args_in, quiet=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
