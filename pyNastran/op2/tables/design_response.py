class WeightResponse(object):
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
    def append(self, internal_id, dresp_id, response_label, region,
               subcase, type_flag, seid,
               row_id, column_id):
        self.internal_id.append(internal_id)
        #self.response_label.append(response_label)
        #self.subcase.append(subcase)
        #self.mode.append(mode)
        #self.mach.append(mach)
        #self.velocity.append(velocity)
        #self.density.append(density)
        #self.flutter_id.append(flutter_id)
        #self.subcase.append(subcase)
        self._n += 1
        #if self.n == self._n:
            #print(self)
    def __repr__(self):
        msg = 'WeightResponse()\n'
        msg += '  n=%s\n' % self.n
        #msg += '  velocity=%s\n' % (velocity)
        return msg

class FlutterResponse(object):
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
        self.response_label = []
        self.subcase = []
        self.mode = []
        self.mach = []
        self.velocity = []
        self.density = []
        self.flutter_id = []
        self.subcase = []

    def append(self, internal_id, dresp_id, response_label, region, subcase, type_flag, seid,
               mode, mach, velocity, density, flutter_id):
        self.internal_id.append(internal_id)
        self.response_label.append(response_label)
        self.subcase.append(subcase)
        self.mode.append(mode)
        self.mach.append(mach)
        self.velocity.append(velocity)
        self.density.append(density)
        self.flutter_id.append(flutter_id)
        self.subcase.append(subcase)
        self._n += 1
        #if self.n == self._n:
            #print(self)
    def __repr__(self):
        msg = 'FlutterResponse()\n'
        msg += '  velocity=%s\n' % (self.velocity)
        return msg

