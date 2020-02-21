"""
Defines:
  - Resultant

"""
from numpy import array
from pyNastran.nptyping import NDArray6float

class Resultant:
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

    def _print_case(self, isubcase, fxyz: NDArray6float):
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
        msg += (
            '0 %8i     FX    %12.6E  %10s    %10s    %10s   %13.6E %13.6E\n'
            '  %8s     FY    %10s    %12.6E  %10s   %13.6E  %10s   %13.6E\n'
            '  %8s     FZ    %10s    %10s    %12.6E %13.6E %13.6E  %10s\n'
            '  %8s     MX    %10s    %10s    %10s   %13.6E  %10s    %10s\n'
            '  %8s     MY    %10s    %10s    %10s    %10s   %13.6E  %10s\n'
            '  %8s     MZ    %10s    %10s    %10s    %10s    %10s   %13.6E\n'
            '  %8s   TOTALS  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E\n' % (
                isubcase,
                fx, dash, dash, dash, 0.0, 0.0,
                '', dash, fy, dash, 0.0, dash, 0.0,
                '', dash, dash, fz, 0.0, 0.0, dash,
                '', dash, dash, dash, 0.0, dash, dash,
                '', dash, dash, dash, dash, 0.0, dash,
                '', dash, dash, dash, dash, dash, 0.0,
                '', fxs, fys, fzs, mxs, mys, mzs,)
        )
        return msg

    def __repr__(self):
        """the str(obj)"""
        msg = ''
        assert len(self.fxyz) > 0, self.fxyz
        for isubcase, fxyz in sorted(self.fxyz.items()):
            msg += self._print_case(isubcase, fxyz) + '\n'
        return msg.rstrip()
