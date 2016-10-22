# pylint: disable=C0301,C0326
"""
Defines:
  - Resultant
"""
from numpy import array

class Resultant(object):
    """interface for making the OLOAD Resultant Table"""
    def __init__(self, table_name, fxyz, isubcase):
        """
        common init for:
          - SPC_Force_Resultant
          - OLOAD_Resultant
         """
        self.table_name = table_name
        assert isubcase is not None
        self.fxyz = {}
        if isubcase is not None:
            self.fxyz[isubcase] = fxyz

    def add(self, fxyz, isubcase):
        """adds a result"""
        self.fxyz[isubcase] = fxyz

    def write_f06(self, f06_file, page_stamp, page_num):
        """writes an f06 file"""
        #msg = ''
        #msg += '        *** USER INFORMATION MESSAGE 7310 (VECPRN)\n'
        #msg += '            ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM WILL BE USED AS REFERENCE LOCATION.\n'
        #msg += '            RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN SUPERELEMENT BASIC SYSTEM COORDINATES.\n'
        #msg += '       0                                                  OLOAD    RESULTANT       \n'

        msg = ''
        assert len(self.fxyz) > 0, self.fxyz
        for isubcase, fxyz in sorted(self.fxyz.items()):
            msg += self._print_case(isubcase, fxyz)
            msg += page_stamp % page_num
            f06_file.write(''.join(msg))
        return page_num + 1

    def _print_case(self, isubcase, fxyz):
        """prints a single resultant case"""
        fx, fy, fz, mx, my, mz = fxyz
        msg = '0                                                  %-8s RESULTANT       \n' % self.table_name
        msg += '  SUBCASE/    LOAD\n'
        msg += '  DAREA ID    TYPE       T1            T2            T3            R1            R2            R3\n'
        force = array([[fx, 0., 0.],
                       [0., fy, 0.],
                       [0., 0., fz]])
        fxs, fys, fzs = force.sum(axis=0)
        moment = array([[mx, 0., 0.],
                        [0., my, 0.],
                        [0., 0., mz]])
        mxs, mys, mzs = moment.sum(axis=0)

        # I think the rest is just an rxF...
        dash = '   ----   '
        msg += '0 %8i     FX    %12.6E  %10s    %10s    %10s   %13.6E %13.6E\n' % (isubcase, fx, dash, dash,   dash, 0.0, 0.0)
        msg += '  %8s     FY    %10s    %12.6E  %10s   %13.6E  %10s   %13.6E\n'   % ('',     dash, fy, dash,   0.0, dash, 0.0)
        msg += '  %8s     FZ    %10s    %10s    %12.6E %13.6E %13.6E  %10s\n'     % ('',     dash, dash, fz,   0.0, 0.0, dash)
        msg += '  %8s     MX    %10s    %10s    %10s   %13.6E  %10s    %10s\n'    % ('',     dash, dash, dash, 0.0, dash, dash)
        msg += '  %8s     MY    %10s    %10s    %10s    %10s   %13.6E  %10s\n'    % ('',     dash, dash, dash, dash, 0.0, dash)
        msg += '  %8s     MZ    %10s    %10s    %10s    %10s    %10s   %13.6E\n'  % ('',     dash, dash, dash, dash, dash, 0.0)
        msg += '  %8s   TOTALS  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E\n' % ('',     fxs,  fys,  fzs,  mxs, mys, mzs)
        return msg

    def __repr__(self):
        """the str(obj)"""
        msg = ''
        assert len(self.fxyz) > 0, self.fxyz
        for isubcase, fxyz in sorted(self.fxyz.items()):
            msg += self._print_case(isubcase, fxyz) + '\n'
        return msg.rstrip()


#class OLOAD_Resultant(Resultant):
    #r"""
    #0                                                  OLOAD    RESULTANT
      #SUBCASE/    LOAD
      #DAREA ID    TYPE       T1            T2            T3            R1            R2            R3
    #0        1     FX    2.300000E+04     ----          ----          ----       3.320987E+04 -2.280395E+04
                   #FY       ----       0.000000E+00     ----       0.000000E+00     ----       0.000000E+00
                   #FZ       ----          ----       0.000000E+00  0.000000E+00  0.000000E+00     ----
                   #MX       ----          ----          ----       0.000000E+00     ----          ----
                   #MY       ----          ----          ----          ----       0.000000E+00     ----
                   #MZ       ----          ----          ----          ----          ----       0.000000E+00
                 #TOTALS  2.300000E+04  0.000000E+00  0.000000E+00  0.000000E+00  3.320987E+04 -2.280395E+04
    ##1    MSC.NASTRAN JOB CREATED ON 28-JAN-12 AT 12:52:32                       OCTOBER  22, 2014  MSC.NASTRAN  6/17/05   PAGE     8
    #"""
    #def __init__(self, fxyz=None, isubcase=None):
        #Resultant.__init__(self, 'OLOAD', fxyz, isubcase)

#class SPC_Force_Resultant(Resultant):
    #"""defines an 'SPCFORCE RESULTANT' table"""
    #def __init__(self, fxyz=None, isubcase=None):
        #Resultant.__init__(self, 'SPCFORCE', fxyz, isubcase)
