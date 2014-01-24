class Zeros(object):
    def __init__(self, size, dtype='i'):
        if isinstance(size, int):
            size = [size]
        else:
            size = list(size)
        size = tuple([sizei+1 for sizei in size])
        self.data = zeros(size, dtype=dtype)

    def __getitem__(self, key):
        print "key_get =", str(key)
        key0 = key[0] - 1
        key1 = key[1]
        key = (key0, key1)
        print "key_get", key
        print "---"
        return self.data[key0, key1]

    def __setitem__(self, key, value):
        print "key_set1 =", str(key)
        key0 = key[0] - 1
        key1 = key[1]
        key = (key0, key1)
        print "key_set2", key
        print "value_set =", value - 1
        print "---"
        self.data[key] = value - 1

