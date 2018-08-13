from collections import OrderedDict


class MyObj(object):
    b = 1
    a = 2

    def __init__(self):
        object.__setattr__(self, '_attrs', OrderedDict())
        self.c = 1
        self.d = 2

    def __setattr__(self, key, value):
        assert key != '_attrs'
        self._attrs[key] = value

    def __getattr__(self, item):
        try:
            return self._attrs[item]
        except KeyError:
            return self.__class__.__dict__[item]

    @property
    def __dict__(self):
        return self._attrs


a = MyObj()

a.e = 3

print(a.__dict__)
print(MyObj.__dict__)

print(a.a)