"""
Defines:
 - FortranFormat

"""
from __future__ import print_function
import sys
from copy import deepcopy
from struct import unpack
from six import iteritems
from six.moves import range

from pyNastran.op2.errors import FortranMarkerError, SortCodeError
from pyNastran.utils import object_attributes

# this is still a requirement, but disabling it so readthedocs works
if sys.version_info < (2, 7, 7):
    IMAJOR, MINOR1, MINOR2 = sys.version_info[:3]
    # makes sure we don't get the following bug:
    #   Issue #19099: The struct module now supports Unicode format strings.
    raise ImportError('Upgrade your Python to >= 2.7.7; version=(%s.%s.%s)' % (
        IMAJOR, MINOR1, MINOR2))

class FortranFormat(object):
    """defines basic methods for reading Fortran formatted data files"""
    def __init__(self):
        self.n = 0
        self.f = None
        #self.obj = None
        #self.data_code = None
        #self.table_name = None
        #self.isubcase = None
        self.binary_debug = None
        self.read_mode = 1
        self._endian = None
        self._table_mapper = {}

        #: stores if the user entered [] for isubcases
        self.is_all_subcases = True
        self.valid_subcases = []

    def show(self, n, types='ifs', endian=None):  # pragma: no cover
        """Shows binary data"""
        return self.op2_reader.show(n, types=types, endian=endian)

    def show_data(self, data, types='ifs', endian=None):  # pragma: no cover
        """Shows binary data"""
        return self.op2_reader.show_data(data, types=types, endian=endian)

    def show_ndata(self, n, types='ifs'):  # pragma: no cover
        self.op2_reader.show_ndata(n, types=types)

    #def passer(self, data):
        #"""
        #dummy function used for unsupported tables
        #"""
        #pass

    def _get_table_mapper(self):
        raise NotImplementedError('this should be overwritten')

    def _finish(self):
        raise NotImplementedError('overwrite this')
