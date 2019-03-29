from __future__ import print_function
import numpy as np

class Responses(object):
    """Defines SOL 200 responses"""
    def __init__(self):
        self.convergence_data = None
        self.weight_response = None
        self.stress_response = None
        self.strain_response = None
        self.force_response = None
        self.composite_stress_response = None
        self.composite_strain_response = None
        self.flutter_response = None

    def get_stats(self, short=False):
        objects = [
            self.convergence_data,
            self.weight_response,
            self.stress_response,
            self.strain_response,
            self.force_response,
            self.composite_stress_response,
            self.composite_strain_response,
            self.flutter_response,
        ]
        msg = []
        for obj in objects:
            if obj is not None:
                msg += obj.get_stats(short=short) + '\n'
        return msg

class WeightResponse(object):
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []

    def add_from_op2(self, out, log):
        """
        Weight Response

        1 IRID         I Internal response identification number
        2 RID          I External response identification number
        3 TYPE(C)      I Response type
        4 LABEL(2) CHAR4 Label
        6 REGION       I Region identifier
        7 SCID         I Subcase identification number
        8 UNDEF(2)     I Not used
        10 SEID        I Superelement identification number or ALL
        11 UNDEF(2)    I Not used

        13 UNDEF       I Not used
        14 TYFLG       I Flag to indicate how response is referenced
        15 SEID        I Superelement identificaiton number

        ---> 3i 8s 7i 3i

        #                             -----  WEIGHT RESPONSE  -----
        # ---------------------------------------------------------------------------------
        #  INTERNAL  DRESP1  RESPONSE  ROW  COLUMN  LOWER     INPUT      OUTPUT     UPPER
        #     ID       ID     LABEL     ID    ID    BOUND     VALUE       VALUE     BOUND
        # ---------------------------------------------------------------------------------
        #       1       1    WEIGHT     3     3       N/A   2.9861E+05  2.9852E+05   N/A

        # ?  ?     ?  LABEL?      ?  ?     ROW_ID? COL_ID? ?  ?  ?          ?          ?  ?
        #(1, 1,    1, 'WEIGHT  ', 0, 1011, 3,      3,      0, 0, 0,         0,         0, 0)
        #(1, 1000, 1, 'W       ', 0, 1,    3,      3,      0, 0, 0,         0,         0, 0)
        #
        # per dev forum; 538976288 is probably just '    '
        #(1, 1,    1, 'WEIGHT  ', 0, 1,    3,      3,      0, 0, 538976288, 538976288, 0, 0)
        """
        # F:\work\pyNastran\examples\Dropbox\move_tpl\mbcgen.op2
        # (1,  15, 1, 'W       ', -1,    1, 3, 3, 0, 0, 0, 0, 0, 10)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\i2002.op2
        # (20, 15, 1, 'W       ', -1,    1, 3, 3, 0, 0, 0, 0, 0, 20)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\edr2n.op2
        # (1, 201, 1, 'WEIGHT  ', -1,    0, 3, 3, 0, 0, 0, 0, 0, 100)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\ss200m2.op2
        # (1,   1, 1, 'W       ', 0, 30001, 3, 3, 0, 0, 0, 0, 0, 300)

        # F:\work\pyNastran\examples\Dropbox\move_tpl\mcso43.op2
        # (1,  10, 1, 'W       ', 0,     1, 3, 3, 0, 0, 0, 0, 1, 0)
        #--------------------------
        # common
        #
        # per the R1TAB DMAP page:
        #   all indicies are downshift by 1
        #   indices above out[3] are off by +2 because of the 2 field response_label
        internal_id = out[0]
        dresp_id = out[1]
        response_type = out[2]
        response_label = out[3].strip()
        # -1 for 2 field wide response_label
        region = out[4]
        subcase = out[5]
        type_flag = out[12]  # no meaning per MSC DMAP 2005
        seid = out[13]

        #--------------------------------------------------
        #row_id = out[4]

        # these should be blank?
        row_id = out[6]
        column_id = out[7]
        seid_weight = out[8]

        if np.abs(out[8:-1]).sum() != 0.0:
            msg = 'WEIGHT response sum error; out=%s 8=%s' % (out, out[8:-1])
            log.warning(msg)
        assert seid == out[-1]
        #assert out[-1] in [0, 1, 2, 3, 4, 5, 10, 20, 100, 200, 300], out
        #dunno_8 = out[8]
        #dunno_9 = out[9]
        #dunno_10 = out[10]
        #dunno_11 = out[11]
        #dunno_12 = out[12]
        #msg = ('WEIGHT - response_type=%r response_label=%r row_id=%r column_id=%r '
            #'6=%r 7=%r 8=%r 9=%r 10=%r 11=%r 12=%r seid=%r' % (
            #response_type, response_label, row_id, column_id,
            #dunno_6, dunno_7, dunno_8, dunno_9, dunno_10, dunno_11, dunno_12, seid))
        #out = unpack(self._endian + b'iii 8s iiff f fffff', data)
        #msg = 'WEIGHT - label=%r region=%s subcase=%s row_id=%r column_id=%r' % (
            #response_label, region, subcase, row_id, column_id)
        self.append(internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    row_id, column_id)
        #print(msg)
        #self.log.debug(msg)

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

    def __repr__(self):
        msg = 'responses.WeightResponse()\n'
        msg += '  n=%s\n' % self.n
        #msg += '  velocity=%s\n' % (velocity)
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.weight_response (%s)' % (self.n)
        else:
            return self.__repr__()


class GeneralResponse(object):
    """common class for StressResponse, StrainResponse, and ForceResponse"""
    def __init__(self):
        self.n = 1
        self._n = 0
        self._itable = 0
        self.is_built = False
        self.internal_id = []
        self.response_label = []
        self.subcase = []
        self.item_code = []
        self.pid = []

    def append(self, internal_id, dresp_id, response_label, region, subcase, type_flag, seid,
               item_code, pid):
        self.internal_id.append(internal_id)
        self.response_label.append(response_label)
        self.subcase.append(subcase)
        self.item_code.append(item_code)
        self.pid.append(pid)
        self._n += 1

    def __repr__(self):
        name = self.__class__.__name__
        msg = 'responses.%s()\n' % name
        msg += '  item_code=%s\n' % (self.item_code)
        msg += '  pid=%s\n' % (self.pid)
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.%s_response (%s)' % (self.name, self.n)
        else:
            return self.__repr__()

class ForceResponse(GeneralResponse):
    name = 'force'
class StressResponse(GeneralResponse):
    name = 'stress'
class StrainResponse(GeneralResponse):
    name = 'strain'

class FlutterResponse(object):
    name = 'flutter'
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
        self._n += 1

    def __repr__(self):
        msg = 'FlutterResponse()\n'
        msg += '  velocity=%s\n' % (self.velocity)
        return msg

    def get_stats(self, short=False):
        if short:
            return 'responses.%s_response (%s)' % (self.name, self.n)
        else:
            return self.__repr__()

class Convergence(object):
    def __init__(self, ndesign_variables):
        self.n = 1
        self._n = 0
        self.is_built = False
        self.ndesign_variables = ndesign_variables
        self.design_iter = []
        self.iconvergence = []  #       1-soft, 2-hard
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
            msg = 'conv_result=%r iconvergence=%r' % (conv_result, iconvergence)
            raise NotImplementedError(msg)
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

    def get_stats(self, short=False):
        if short:
            return 'responses.convergence_data (%s, %s)' % (self.n, self.ndesign_variables)
        else:
            return self.__repr__()
