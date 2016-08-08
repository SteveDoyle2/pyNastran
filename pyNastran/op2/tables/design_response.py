import numpy as np


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


class Convergence(object):
    def __init__(self, ndesign_variables):
        self.n = 1
        self._n = 0
        self.is_built = False
        self.ndesign_variables = ndesign_variables
        self.design_iter = []
        self.iconvergence = []  #      1-soft, 2-hard
        self.conv_result = []   # 0-no, 1-soft, 2-hard
        self.obj_initial = []
        self.obj_final = []
        self.constraint_max = []
        self.row_constraint_max = []
        self.desvar_values = []

    def get_convergence(self):
        iconvergence = self.iconvergence[-1] # hard, soft (2)
        conv_result = self.conv_result[-1] # no, hard, soft (3)
        # 2*3 = 6
        if (conv_result, iconvergence) == ('hard', 'hard'):  # sure about this
            convergence = 'HARD'
        elif (conv_result, iconvergence) == ('no', 'hard'):  # pretty sure about this...
            convergence = 'MAX DESIGN CYCLES'
        elif (conv_result, iconvergence) == ('no', 'soft'):  # is this possible???
            convergence = 'MAX DESIGN CYCLES'

        # are these possible???
        # we'll just assume SOFT convergence
        if (conv_result, iconvergence) == ('hard', 'soft'):  # pretty sure about this...
            convergence = 'SOFT'
        elif (conv_result, iconvergence) == ('soft', 'hard'): # is this possible???
            convergence = 'SOFT'
        elif (conv_result, iconvergence) == ('soft', 'soft'): # is this possible???
            convergence = 'SOFT'
        else:
            raise NotImplementedError('conv_result=%r iconvergence=%r' % (conv_result, iconvergence))
        return convergence, (conv_result, iconvergence)


    def append(self, design_iter, iconvergence, conv_result, obj_initial, obj_final,
                 constraint_max, row_constraint_max, desvar_values):
        if not self.is_built:
            self.design_iter = np.zeros(self.n, dtype='int32')
            self.iconvergence = np.zeros(self.n, dtype=object)
            self.conv_result = np.zeros(self.n, dtype=object)
            self.obj_initial = np.zeros(self.n, dtype='float32')
            self.obj_final = np.zeros(self.n, dtype='float32')
            self.constraint_max = np.zeros(self.n, dtype='float32')
            self.row_constraint_max = np.zeros(self.n, dtype='int32')
            self.desvar_values = np.zeros((self.n, self.ndesign_variables), dtype='float32')
            self.is_built = True

        n = self._n
        self.design_iter[n] = design_iter
        self.iconvergence[n] = iconvergence
        self.conv_result[n] = conv_result
        self.obj_initial[n] = obj_initial
        self.obj_final[n] = obj_final
        self.constraint_max[n] = constraint_max
        self.row_constraint_max[n] = row_constraint_max
        self.desvar_values[n, :] = desvar_values
        self._n += 1
        #if self.n == self._n:
            #print(self)

    def __repr__(self):
        msg = 'Convergence()\n'
        msg += '  shape=(%s, %s)\n' % (self.n, self.ndesign_variables)
        msg += '  design_iter = %s\n' % self.design_iter
        msg += '  icovergence = %s\n' % self.iconvergence
        msg += '  conv_result = %s\n' % self.conv_result
        msg += '  obj_initial = %s\n' % self.obj_initial
        msg += '  constraint_max = %s\n' % self.constraint_max
        msg += '  row_constraint_max = %s\n' % self.row_constraint_max
        return msg
