from __future__ import print_function

try:
    import numexpr as ne
except ImportError:
    pass
from numpy import array, zeros


class ScalarResult(object):
    def __init__(self, data, uname='name', title='title', dim=1, location='node', fmt='%i'):
        self.uname = uname
        self.data = data
        self.title = title
        self.dim = dim
        self.location = location
        self.fmt = fmt
        assert location in ['node', 'centroid'], location
        assert dim in [1, 3], dim

    def get_methods(self, in_args):
        return [self.location]

    def get_data(self, in_args, method):
        return self.data

    def __repr__(self):
        msg = 'ScalarResult\n'
        msg += '    uname=%r\n' % self.uname
        msg += '    method=%r' % self.location
        return msg

#class ComboResult(object):

class FuncResult(object):
    def __init__(self, get_data, uname='name', title='title', dim=1, location='node', fmt='%i'):
        self.get_data = get_data
        self.uname = uname
        self.title = title
        self.dim = dim
        self.location = location
        self.fmt = fmt
        assert location in ['node', 'centroid'], location
        assert dim in [1, 3], dim

    def get_methods(self, in_args):
        return [self.location]

    def __repr__(self):
        msg = 'FuncResult\n'
        msg += '    uname=%r\n' % self.uname
        msg += '    location=%r' % self.location
        return msg


class Results(object):
    def __init__(self, *results):
        self.results = {}
        self.imap = {}
        for i, result in enumerate(results):
            print(i, result)
            uname = result.uname
            if hasattr(result, 'n'):
                nres = result.n
                for ii in range(nres):
                    self.imap[i + ii] = uname
            else:
                self.imap[i] = uname
            self.results[uname] = result

    def get_data(self, i, method):
        ii, name = i
        res = self.results[self.imap[ii]]
        return res.get_data(i, method)

    def get_methods(self, i):
        #res = self.results[i]
        ii, name = i
        res = self.results[self.imap[ii]]
        return res.get_methods(i)

    def get_uname(self, i):
        ii, name = i
        return self.imap[ii]

    def __repr__(self):
        msg = 'Results\n'
        msg += '    keys=%s' % self.results.keys()
        return msg


class Result(object):
    def __init__(self, uname='name', title='title', dim=1, location='node', fmt='%i'):
        self.uname = uname
        self.title = title
        self.dim = dim
        self.location = location
        self.fmt = fmt
        assert location in ['node', 'centroid'], location
        assert dim in [1, 3], dim

    def get_data(self, i, method):
        raise NotImplementedError()
        #res = self.results[self.imap[i]]
        #return res.get_data(self.results, method)

    def get_methods(self, i):
        raise NotImplementedError()
        #res = self.results[self.imap[i]]
        #return res.get_methods(i)

    def get_uname(self, i):
        return self.uname

    def __repr__(self):
        raise NotImplementedError()
        #msg = 'Results\n'
        #msg += '    keys=%s' % self.results.keys()
        #return msg



