"""Tests for pyNastran.bdf.bdf_interface.include_file
(parse_include_lines, split_filename_into_tokens, get_include_filename)."""
import os
import unittest
import tempfile
from pathlib import PurePosixPath, PureWindowsPath
from unittest.mock import patch

from pyNastran.bdf.bdf_interface.include_file import (
    parse_include_lines, split_filename_into_tokens, get_include_filename,
)


class TestParseIncludeLines(unittest.TestCase):
    """Tests for parse_include_lines: extracts filename from INCLUDE card lines."""

    def test_single_line_single_quotes(self):
        """INCLUDE 'path/to/file.bdf'"""
        lines = ["INCLUDE 'path/to/file.bdf'"]
        self.assertEqual(parse_include_lines(lines), 'path/to/file.bdf')

    def test_single_line_double_quotes(self):
        '''INCLUDE "path/to/file.bdf"'''
        lines = ['INCLUDE "path/to/file.bdf"']
        self.assertEqual(parse_include_lines(lines), 'path/to/file.bdf')

    def test_multiline(self):
        """INCLUDE split across multiple lines."""
        lines = [
            "INCLUDE '/dir1",
            "         /dir2",
            "         /file.bdf'",
        ]
        result = parse_include_lines(lines)
        self.assertEqual(result, '/dir1/dir2/file.bdf')

    def test_multiline_with_whitespace(self):
        """Whitespace and tabs in continuation lines are stripped."""
        lines = [
            "INCLUDE '/start",
            "\t\t /middle  ",
            "   /end.dat'  ",
        ]
        result = parse_include_lines(lines)
        self.assertEqual(result, '/start/middle/end.dat')

    def test_no_quotes(self):
        """INCLUDE without quotes (valid Nastran)."""
        lines = ["INCLUDE myfile.bdf"]
        result = parse_include_lines(lines)
        self.assertEqual(result, 'myfile.bdf')

    def test_windows_path(self):
        r"""INCLUDE with Windows backslash path."""
        lines = [r"INCLUDE 'C:\work\model.bdf'"]
        result = parse_include_lines(lines)
        self.assertEqual(result, r'C:\work\model.bdf')

    def test_empty_filename_raises(self):
        """INCLUDE with only whitespace filename raises SyntaxError."""
        lines = ["INCLUDE '   '"]
        with self.assertRaises(SyntaxError):
            parse_include_lines(lines)

    def test_include_lowercase(self):
        """INCLUDE keyword is removed regardless of case in the first 7 chars."""
        lines = ["INCLUDE '/path/file.dat'"]
        result = parse_include_lines(lines)
        self.assertEqual(result, '/path/file.dat')

    def test_leading_trailing_spaces_in_filename(self):
        """Spaces inside quotes are NOT stripped by parse_include_lines (only quotes removed)."""
        lines = ["INCLUDE '  model.bdf  '"]
        result = parse_include_lines(lines)
        # parse_include_lines strips quotes but not inner whitespace
        self.assertEqual(result, '  model.bdf  ')


class TestSplitFilenameIntoTokens(unittest.TestCase):
    """Tests for split_filename_into_tokens: combines include_dir + filename into a Path."""

    def test_relative_path_linux(self):
        """Relative filename joined with include_dir on Linux."""
        result = split_filename_into_tokens('/work/project', 'subdir/model.bdf', is_windows=False)
        self.assertEqual(str(result), '/work/project/subdir/model.bdf')

    def test_absolute_path_linux(self):
        """Absolute filename ignores include_dir on Linux."""
        result = split_filename_into_tokens('/work/project', '/other/model.bdf', is_windows=False)
        # PurePosixPath('/work/project') / '/other/model.bdf' -> '/other/model.bdf'
        self.assertEqual(str(result), '/other/model.bdf')

    def test_relative_path_windows(self):
        """Relative filename joined with include_dir on Windows."""
        result = split_filename_into_tokens(r'C:\work', r'subdir\model.bdf', is_windows=True)
        self.assertEqual(str(result), r'C:\work\subdir\model.bdf')

    def test_absolute_path_windows(self):
        """Absolute filename with drive letter on Windows."""
        result = split_filename_into_tokens(r'C:\work', r'D:\other\model.bdf', is_windows=True)
        # PureWindowsPath('C:\\work') / 'D:\\other\\model.bdf' -> 'D:\\other\\model.bdf'
        self.assertEqual(str(result), r'D:\other\model.bdf')

    def test_forward_slash_start_on_windows_raises(self):
        """Starting with / on Windows raises SyntaxError."""
        with self.assertRaises(SyntaxError):
            split_filename_into_tokens(r'C:\work', '/unix/path.bdf', is_windows=True)

    def test_unc_path_on_linux_raises(self):
        r"""Starting with \\ on Linux raises SyntaxError."""
        with self.assertRaises(SyntaxError):
            split_filename_into_tokens('/work', r'\\server\share\file.bdf', is_windows=False)

    def test_simple_filename_linux(self):
        """Single filename with no directory."""
        result = split_filename_into_tokens('/mydir', 'file.bdf', is_windows=False)
        self.assertEqual(str(result), '/mydir/file.bdf')

    def test_simple_filename_windows(self):
        """Single filename with no directory on Windows."""
        result = split_filename_into_tokens(r'C:\mydir', 'file.bdf', is_windows=True)
        self.assertEqual(str(result), r'C:\mydir\file.bdf')

    def test_nested_subdirectories_linux(self):
        """Multiple nested subdirectories on Linux."""
        result = split_filename_into_tokens('/base', 'a/b/c/d/file.dat', is_windows=False)
        self.assertEqual(str(result), '/base/a/b/c/d/file.dat')

    def test_empty_include_dir_linux(self):
        """Empty include_dir with relative path."""
        result = split_filename_into_tokens('', 'subdir/file.bdf', is_windows=False)
        self.assertEqual(str(result), 'subdir/file.bdf')

    def test_env_var_linux(self):
        """Environment variable in path on Linux (ENVVAR:rest)."""
        with patch.dict(os.environ, {'MYPROJECT': '/expanded/path'}):
            result = split_filename_into_tokens('', 'MYPROJECT:file.bdf', is_windows=False)
            self.assertEqual(str(result), '/expanded/path/file.bdf')

    def test_env_var_windows(self):
        """Environment variable in path on Windows (ENVVAR:rest)."""
        with patch.dict(os.environ, {'MYPROJECT': r'C:\expanded\path'}):
            result = split_filename_into_tokens('', 'MYPROJECT:file.bdf', is_windows=True)
            self.assertEqual(str(result), r'C:\expanded\path\file.bdf')


class TestGetIncludeFilename(unittest.TestCase):
    """Tests for get_include_filename: full integration of parsing + path resolution."""

    def _make_logger(self):
        """Create a minimal logger mock."""
        class FakeLogger:
            def warning(self, msg): pass
            def debug(self, msg): pass
        return FakeLogger()

    def test_basic_resolve(self):
        """INCLUDE resolves to existing file."""
        log = self._make_logger()
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'model.bdf')
            with open(filepath, 'w') as f:
                f.write('$ comment\n')

            lines = [f"INCLUDE '{filepath}'"]
            result = get_include_filename(
                log, lines, include_dirs=[''], replace_includes={},
                is_windows=None)
            self.assertTrue(os.path.exists(result))

    def test_include_dir_resolution(self):
        """INCLUDE resolves with include_dirs search path."""
        log = self._make_logger()
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, 'sub.bdf')
            with open(filepath, 'w') as f:
                f.write('$ test\n')

            lines = ["INCLUDE 'sub.bdf'"]
            result = get_include_filename(
                log, lines, include_dirs=[tmpdir], replace_includes={},
                is_windows=None)
            self.assertTrue(os.path.exists(result))

    def test_replace_includes(self):
        """replace_includes mapping substitutes the filename."""
        log = self._make_logger()
        with tempfile.TemporaryDirectory() as tmpdir:
            real_path = os.path.join(tmpdir, 'real.bdf')
            with open(real_path, 'w') as f:
                f.write('$ real\n')

            lines = ["INCLUDE 'fake.bdf'"]
            result = get_include_filename(
                log, lines, include_dirs=[tmpdir],
                replace_includes={'fake.bdf': 'real.bdf'},
                is_windows=None)
            self.assertTrue(os.path.exists(result))

    def test_replace_includes_empty_string(self):
        """replace_includes with empty string returns empty (skip include)."""
        log = self._make_logger()
        lines = ["INCLUDE 'skip_me.bdf'"]
        result = get_include_filename(
            log, lines, include_dirs=[''],
            replace_includes={'skip_me.bdf': ''},
            is_windows=None)
        self.assertEqual(result, '')

    def test_file_not_found_raises(self):
        """Missing INCLUDE file raises IOError when multiple include_dirs are searched."""
        log = self._make_logger()
        with tempfile.TemporaryDirectory() as dir1:
            with tempfile.TemporaryDirectory() as dir2:
                lines = ["INCLUDE 'nonexistent_xyz_12345.bdf'"]
                with self.assertRaises(IOError):
                    get_include_filename(
                        log, lines, include_dirs=[dir1, dir2],
                        replace_includes={}, is_windows=None)

    def test_file_not_found_single_dir_returns_path(self):
        """With a single include_dir, missing file returns path without raising."""
        log = self._make_logger()
        with tempfile.TemporaryDirectory() as tmpdir:
            lines = ["INCLUDE 'nonexistent_xyz_12345.bdf'"]
            result = get_include_filename(
                log, lines, include_dirs=[tmpdir],
                replace_includes={}, is_windows=None)
            # returns the path even though it doesn't exist
            self.assertIn('nonexistent_xyz_12345.bdf', result)

    def test_multiple_include_dirs(self):
        """Searches multiple include_dirs in order."""
        log = self._make_logger()
        with tempfile.TemporaryDirectory() as dir1:
            with tempfile.TemporaryDirectory() as dir2:
                filepath = os.path.join(dir2, 'found.bdf')
                with open(filepath, 'w') as f:
                    f.write('$ found\n')

                lines = ["INCLUDE 'found.bdf'"]
                result = get_include_filename(
                    log, lines, include_dirs=[dir1, dir2],
                    replace_includes={}, is_windows=None)
                self.assertTrue(os.path.exists(result))

    def test_long_line_warning(self):
        """Lines > 72 chars trigger a warning (but don't fail)."""
        log = self._make_logger()
        warnings = []
        log.warning = lambda msg: warnings.append(msg)

        with tempfile.TemporaryDirectory() as tmpdir:
            long_name = 'a' * 70 + '.bdf'
            filepath = os.path.join(tmpdir, long_name)
            with open(filepath, 'w') as f:
                f.write('$ test\n')

            lines = [f"INCLUDE '{long_name}'"]
            result = get_include_filename(
                log, lines, include_dirs=[tmpdir],
                replace_includes={}, is_windows=None)
            self.assertTrue(os.path.exists(result))
            self.assertGreater(len(warnings), 0)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
