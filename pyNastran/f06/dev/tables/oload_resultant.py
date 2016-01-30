from numpy import array

class OLOAD_Resultant(object):
    def __init__(self):
        self.Fx = 2.3e4

    def add(self, F):
        pass

    def write_f06(self, f, page_stamp, page_num):
        #msg = ''
        #msg += '        *** USER INFORMATION MESSAGE 7310 (VECPRN)\n'
        #msg += '            ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM WILL BE USED AS REFERENCE LOCATION.\n'
        #msg += '            RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN SUPERELEMENT BASIC SYSTEM COORDINATES.\n'
        #msg += '       0                                                  OLOAD    RESULTANT       \n'

        """
        0                                                  OLOAD    RESULTANT
          SUBCASE/    LOAD
          DAREA ID    TYPE       T1            T2            T3            R1            R2            R3
        0        1     FX    2.300000E+04     ----          ----          ----       3.320987E+04 -2.280395E+04
                       FY       ----       0.000000E+00     ----       0.000000E+00     ----       0.000000E+00
                       FZ       ----          ----       0.000000E+00  0.000000E+00  0.000000E+00     ----
                       MX       ----          ----          ----       0.000000E+00     ----          ----
                       MY       ----          ----          ----          ----       0.000000E+00     ----
                       MZ       ----          ----          ----          ----          ----       0.000000E+00
                     TOTALS  2.300000E+04  0.000000E+00  0.000000E+00  0.000000E+00  3.320987E+04 -2.280395E+04
        #1    MSC.NASTRAN JOB CREATED ON 28-JAN-12 AT 12:52:32                       OCTOBER  22, 2014  MSC.NASTRAN  6/17/05   PAGE     8
        """
        isubcase = 1
        Fx = 2.3e4
        msg = '0                                                  OLOAD    RESULTANT\n'
        msg += '  SUBCASE/    LOAD\n'
        msg += '  DAREA ID    TYPE       T1            T2            T3            R1            R2            R3\n'
        F = array([[2.3e4, 0., 0.],
                   [0., 0., 0.],
                   [0., 0., 0.]])
        M = array([[0., 3.320987e4, -2.280395e4],
                   [0., 0., 0.],
                   [0., 0., 0.]])
        msg += '0 %8i     FX    %12.6E     ----          ----       3.320987E+04 -2.280395E+04\n' % (isubcase, Fx)
        msg += '  %8s     FY    %12s    0.000000E+00     ----       0.000000E+00     ----       0.000000E+00\n' % ('', '----')
        msg += '  %8s     FZ    %12s       ----       0.000000E+00  0.000000E+00  0.000000E+00     ----\n' % ('', '----')
        msg += '  %8s     MX    %12s       ----          ----       0.000000E+00     ----          ----\n' % ('', '----')
        msg += '  %8s     MY    %12s       ----          ----          ----       0.000000E+00     ----\n' % ('', '----')
        msg += '  %8s     MZ    %12s       ----          ----          ----          ----       0.000000E+00\n' % ('', '----')
        msg += '  %8s     TOTALS  2.300000E+04  0.000000E+00  0.000000E+00  0.000000E+00  3.320987E+04 -2.280395E+04\n' % ('')

        msg += page_stamp % page_num
        f.write('\n'.join(msg))
        return page_num + 1

