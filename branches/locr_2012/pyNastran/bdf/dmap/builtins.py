class Number(object):
    def __init__(self, name):
        self.name = name

    def __add__(self, value):
        return '%s + %s' % (self.real, value)


class ZEROS(object):
    def __init__(self, shape, dtype='i'):
        self.name = None
        self.shape = shape
        self.dtype = dtype

    def shape(self):
        return self.shape


class ONES(object):
    def __init__(self, name, shape, dtype='i'):
        self.name = None
        self.shape = shape
        self.dtype = dtype

    def shape(self):
        return self.shape


class EYE(object):
    def __init__(self, name, shape, dtype='i'):
        self.name = None
        self.shape = (shape, shape)
        self.dtype = dtype

    def shape(self):
        return self.shape


class STRING(object):
    def __init__(self, name, word):
        self.name = name
        self.word = word  # .strip("'")

    def __add__(self, value):
        return '%s + %s' % (self.word, value)

    def length(self, globalsDict):
        #print "self.word = ",self.word
        val = globalsDict[self.name]
        #print "val = |%s|" %(val)
        return len(val.__repr__())

    def __repr__(self):
        return str(self.word)


class INT(Number):
    def __init__(self, name, real):
        Number.__init__(self, name)
        self.real = real


class FLOAT(Number):
    def __init__(self, name, real):
        Number.__init__(self, name)
        self.real = real


class DOUBLE(Number):
    def __init__(self, real):
        Number.__init__(self, None)
        self.real = real


class COMPLEX(Number):
    def __init__(self, name, real, imag):
        Number.__init__(self, name)
        self.real = real
        self.imag = imag

    def __add__(self, value):
        complexAdd


class ARRAY(object):
    def __init__(self, name, nrows, ncols):
        self.name = name
        self.nrows = nrows
        self.ncols = ncols

    def shape(self):
        return (self.nrows, self.ncols)
