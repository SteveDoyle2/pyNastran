import os
from typing import TextIO, BinaryIO


class DynamicFileWriter:
    def __init__(self, filename: str,
                 file_obj: TextIO | BinaryIO = None):
        ext = os.path.splitext(filename)[1][1:].lower()
        assert ext in ['h5', 'op4'], ext
        self.filename = filename
        file_obj = None
        self.file_obj = file_obj
        self.write = getattr(self, f'write_{ext}')
        self.ext = ext
        self.matrices_written = set([])

    def __repr__(self):
        opened = self.file_obj is not None
        return f'DynamicFileWriter({self.ext}, opened={opened})'

    #def __del__(self):
    #    self.close()
    #    del self.write

    def write_h5(self, name: str, matrix):
        if self._exists(name, matrix):
            return
        pass

    def write_f06(self, name: str, matrix):
        if self._exists(name, matrix):
            return
        if self.file_obj is None:
            self.file_obj = open(self.filename, 'w')
        #write_f06_matrix(matrix)

    def _exists(self, name, matrix) -> bool:
        if matrix is None or name in self.matrices_written:
            return True
        self.matrices_written.add(name)
        return False

    def write_op4(self, name: str, matrix):
        if self._exists(name, matrix):
            return
        if self.file_obj is None:
            self.file_obj = open(self.filename, 'w')
        #write_op4_matrix(matrix)

    def close(self):
        if self.file_obj is None:
            return
        self.file_obj.close()
        self.file_obj = None
