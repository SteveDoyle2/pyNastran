from pyNastran.op2.tables.lama_eigenvalues.lama_objects import RealEigenvalues, ComplexEigenvalues

class LAMA(object):

    def __init__(self):
        pass

    def _real_eigenvalues(self):
        """
        ::

                                                     R E A L   E I G E N V A L U E S
           MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
            NO.       ORDER                                                                       MASS              STIFFNESS
                1         1        6.158494E+07        7.847607E+03        1.248985E+03        1.000000E+00        6.158494E+07
        """
        self.title = None
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        Title = None
        line1 = self.infile.readline().strip(); self.i += 1
        if line1 != 'MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED':
            Title = line1
            line1 = self.infile.readline().strip(); self.i += 1
        line2 = self.infile.readline().strip(); self.i += 1

        #MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
        # NO.       ORDER                                                                       MASS              STIFFNESS
        #     1         1        1.018377E-03        3.191203E-02        5.078956E-03        1.000000E+00        1.018377E-03
        #print(line1)
        #print(line2)
        #headers = self.skip(2)
        #print(headers)
        data = self._read_f06_table([int, int, float, float, float, float, float])

        self.eigenvalues[self.title] = RealEigenvalues(Title)
        self.eigenvalues[self.title].add_f06_data(data)

    def _complex_eigenvalue_summary(self):
        """
        ::

                                 C O M P L E X   E I G E N V A L U E   S U M M A R Y
          ROOT     EXTRACTION                  EIGENVALUE                     FREQUENCY              DAMPING
           NO.        ORDER             (REAL)           (IMAG)                (CYCLES)            COEFFICIENT
               1           6          0.0              6.324555E+01          1.006584E+01          0.0
               2           5          0.0              6.324555E+01          1.006584E+01          0.0
        """
        #(subcaseName,isubcase,transient,dt,analysis_code,is_sort1) = self.readSubcaseNameID()
        isubcase = 1  # .. todo:: fix this...

        headers = self.skip(2)
        data = self._read_f06_table([int, int, float, float, float, float])

        if self.title in self.eigenvalues:
            self.eigenvalues[self.title].add_f06_data(data)
        else:
            self.eigenvalues[self.title] = ComplexEigenvalues(self.title)
            self.eigenvalues[self.title].add_f06_data(data)
        self.iSubcases.append(isubcase)
