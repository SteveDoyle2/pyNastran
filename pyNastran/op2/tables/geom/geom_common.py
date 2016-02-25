
class SuppressLogging(object):
    def __init__(self):
        pass
    def debug(self, msg):
        pass
    def info(self, msg):
        pass
    def flush(self):
        pass

class SuppressFileIO(object):
    def __init__(self):
        pass
    def open(self, fname):
        pass
    def close(self):
        pass
    def write(self, msg):
        pass
    def flush(self):
        pass

class GeomCommon(object):
    def __init__(self):
        self.card_count = {}
        self.is_debug_file = False
        self._endian = ''
        self.binary_debug = SuppressFileIO()
        self.log = SuppressLogging()

    def _read_fake(self, data, n):
        return len(data)
