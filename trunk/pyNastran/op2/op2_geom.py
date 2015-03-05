from pyNastran.op2.tables.geom.geom1 import GEOM1
from pyNastran.op2.tables.geom.geom2 import GEOM2
from pyNastran.op2.tables.geom.geom3 import GEOM3
from pyNastran.op2.tables.geom.geom4 import GEOM4

from pyNastran.op2.tables.geom.ept import EPT
from pyNastran.op2.tables.geom.mpt import MPT

from pyNastran.op2.tables.geom.dit import DIT
from pyNastran.op2.tables.geom.dynamics import DYNAMICS

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_vectorized import OP2_Vectorized

class OP2Geom(BDF,
              GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT, DIT, DYNAMICS,
              OP2):
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param debug: enables the debug log and sets the debug in the logger (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        :param debug_file: sets the filename that will be written to (default=None -> no debug)
        """
        assert make_geom == True, make_geom

        BDF.__init__(self, debug=debug, log=log)
        GEOM1.__init__(self)
        GEOM2.__init__(self)
        GEOM3.__init__(self)
        GEOM4.__init__(self)

        EPT.__init__(self)
        MPT.__init__(self)
        DIT.__init__(self)
        DYNAMICS.__init__(self)

        OP2.__init__(self, debug, log=log, debug_file=debug_file)
        self.make_geom = make_geom

    def _get_table_mapper(self):
        table_mapper = OP2._get_table_mapper(self)

        # geometry
        table_mapper['GEOM1'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4'] = [self._read_geom4_4, self._read_geom4_4]

        # superelements
        table_mapper['GEOM1S'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2S'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3S'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4S'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper['GEOM1N'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2N'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3N'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4N'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper['GEOM1OLD'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper['GEOM2OLD'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper['GEOM3OLD'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper['GEOM4OLD'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper['EPT']    = [self._read_ept_4, self._read_ept_4]
        table_mapper['EPTS']   = [self._read_ept_4, self._read_ept_4]
        table_mapper['EPTOLD'] = [self._read_ept_4, self._read_ept_4]

        table_mapper['MPT']  = [self._read_mpt_4, self._read_mpt_4]
        table_mapper['MPTS'] = [self._read_mpt_4, self._read_mpt_4]

        table_mapper['DYNAMIC']  = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper['DYNAMICS'] = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper['DIT']      = [self._read_dit_4, self._read_dit_4]   # table objects (e.g. TABLED1)

        return table_mapper

        table_mapper = {
            #=======================
            # OEF
            # element forces
            'OEFIT' : [self._read_oef1_3, self._read_oef1_4],
            'OEF1X' : [self._read_oef1_3, self._read_oef1_4],
            'OEF1'  : [self._read_oef1_3, self._read_oef1_4],
            'DOEF1' : [self._read_oef1_3, self._read_oef1_4],
            #=======================
            # OQG
            # spc forces
            'OQG1'  : [self._read_oqg1_3, self._read_oqg1_4],  # spc forces in the nodal frame
            'OQGV1' : [self._read_oqg1_3, self._read_oqg1_4],  # spc forces in the nodal frame
            # mpc forces
            'OQMG1' : [self._read_oqg1_3, self._read_oqg1_4],

            # ???? - passer
            #'OQP1': [self._table_passer, self._table_passer],
            #=======================
            # OPG
            # applied loads
            'OPG1'  : [self._read_opg1_3, self._read_opg1_4],  # applied loads in the nodal frame
            'OPGV1' : [self._read_opg1_3, self._read_opg1_4],
            'OPNL1' : [self._read_opg1_3, self._read_opg1_4],

            # OGPFB1
            # grid point forces
            'OGPFB1' : [self._read_ogpf1_3, self._read_ogpf1_4],  # grid point forces

            # ONR/OEE
            # strain energy density
            'ONRGY1' : [self._read_onr1_3, self._read_onr1_4],  # strain energy density
            #=======================
            # OES
            # stress
            'OES1X1'  : [self._read_oes1_3, self._read_oes1_4],  # stress
            'OES1'    : [self._read_oes1_3, self._read_oes1_4],  # stress
            'OES1X'   : [self._read_oes1_3, self._read_oes1_4],  # stress
            'OES1C'   : [self._read_oes1_3, self._read_oes1_4],  # stress - composite
            'OESCP'   : [self._read_oes1_3, self._read_oes1_4],
            'OESNLXR' : [self._read_oes1_3, self._read_oes1_4],
            'OESNLXD' : [self._read_oes1_3, self._read_oes1_4],
            'OESNLBR' : [self._read_oes1_3, self._read_oes1_4],
            'OESTRCP' : [self._read_oes1_3, self._read_oes1_4],
            'OESNL1X' : [self._read_oes1_3, self._read_oes1_4],
            'OESRT'   : [self._read_oes1_3, self._read_oes1_4],

            # strain
            'OSTR1X'  : [self._read_oes1_3, self._read_oes1_4],  # strain - isotropic
            'OSTR1C'  : [self._read_oes1_3, self._read_oes1_4],  # strain - composite

            #=======================
            # OUG
            # displacement/velocity/acceleration/eigenvector/temperature
            'OUG1'    : [self._read_oug1_3, self._read_oug1_4],  # displacements in nodal frame
            'OUGV1'   : [self._read_oug1_3, self._read_oug1_4],  # displacements in nodal frame
            'BOUGV1'  : [self._read_oug1_3, self._read_oug1_4],  # OUG1 on the boundary???
            'OUGV1PAT': [self._read_oug1_3, self._read_oug1_4],  # OUG1 + coord ID
            'OUPV1'   : [self._read_oug1_3, self._read_oug1_4],

            #=======================
            # OGPWG
            # grid point weight
            'OGPWG'  : [self._read_ogpwg_3, self._read_ogpwg_4],  # grid point weight
            'OGPWGM' : [self._read_ogpwg_3, self._read_ogpwg_4],  # modal? grid point weight

            #=======================
            # OGS
            # grid point stresses
            'OGS1' : [self._read_ogs1_3, self._read_ogs1_4],  # grid point stresses
            #=======================
            # eigenvalues
            'BLAMA': [self._read_buckling_eigenvalue_3, self._read_buckling_eigenvalue_4],  # buckling eigenvalues
            'CLAMA': [self._read_complex_eigenvalue_3,  self._read_complex_eigenvalue_4],   # complex eigenvalues
            'LAMA' : [self._read_real_eigenvalue_3,     self._read_real_eigenvalue_4],      # eigenvalues

            # ===geom passers===
            # geometry
            b'GEOM1': [self._table_passer, self._table_passer],
            b'GEOM2': [self._table_passer, self._table_passer],
            b'GEOM3': [self._table_passer, self._table_passer],
            b'GEOM4': [self._table_passer, self._table_passer],

            # superelements
            b'GEOM1S': [self._table_passer, self._table_passer],  # GEOMx + superelement
            b'GEOM2S': [self._table_passer, self._table_passer],
            b'GEOM3S': [self._table_passer, self._table_passer],
            b'GEOM4S': [self._table_passer, self._table_passer],

            b'GEOM1N': [self._table_passer, self._table_passer],
            b'GEOM2N': [self._table_passer, self._table_passer],
            b'GEOM3N': [self._table_passer, self._table_passer],
            b'GEOM4N': [self._table_passer, self._table_passer],

            b'GEOM1OLD': [self._table_passer, self._table_passer],
            b'GEOM2OLD': [self._table_passer, self._table_passer],
            b'GEOM3OLD': [self._table_passer, self._table_passer],
            b'GEOM4OLD': [self._table_passer, self._table_passer],

            b'EPT' : [self._table_passer, self._table_passer],  # elements
            b'EPTS': [self._table_passer, self._table_passer],  # elements - superelements
            b'EPTOLD' : [self._table_passer, self._table_passer],

            'MPT' : [self._table_passer, self._table_passer],  # materials
            'MPTS': [self._table_passer, self._table_passer],  # materials - superelements

            'DYNAMIC': [self._table_passer, self._table_passer],
            'DYNAMICS': [self._table_passer, self._table_passer],
            'DIT': [self._table_passer, self._table_passer],

            # geometry
            #'GEOM1': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4': [self._read_geom4_4, self._read_geom4_4],

            # superelements
            #'GEOM1S': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2S': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3S': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4S': [self._read_geom4_4, self._read_geom4_4],

            #'GEOM1N': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2N': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3N': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4N': [self._read_geom4_4, self._read_geom4_4],

            #'GEOM1OLD': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2OLD': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3OLD': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4OLD': [self._read_geom4_4, self._read_geom4_4],

            #'EPT' : [self._read_ept_4, self._read_ept_4],
            #'EPTS': [self._read_ept_4, self._read_ept_4],
            #'EPTOLD' : [self._read_ept_4, self._read_ept_4],

            #'MPT' : [self._read_mpt_4, self._read_mpt_4],
            #'MPTS': [self._read_mpt_4, self._read_mpt_4],

            #'DYNAMIC': [self._read_dynamics_4, self._read_dynamics_4],
            #'DYNAMICS': [self._read_dynamics_4, self._read_dynamics_4],
            #'DIT': [self._read_dit_4, self._read_dit_4],   # table objects (e.g. TABLED1)

            # ===passers===
            'EQEXIN': [self._table_passer, self._table_passer],
            'EQEXINS': [self._table_passer, self._table_passer],

            'GPDT': [self._table_passer, self._table_passer],     # grid points?
            'BGPDT': [self._table_passer, self._table_passer],
            'BGPDTS': [self._table_passer, self._table_passer],
            'BGPDTOLD': [self._table_passer, self._table_passer],

            'PVT0': [self._table_passer, self._table_passer],
            'DESTAB': [self._table_passer, self._table_passer],
            'STDISP': [self._table_passer, self._table_passer],
            'R1TABRG': [self._table_passer, self._table_passer],
            'CASECC': [self._table_passer, self._table_passer],  # case control deck

            'HISADD': [self._table_passer, self._table_passer],  # optimization history (SOL200)
            'EDTS': [self._table_passer, self._table_passer],
            'FOL': [self._table_passer, self._table_passer],
            'MONITOR': [self._table_passer, self._table_passer],  # monitor points
            'PERF': [self._table_passer, self._table_passer],
            'VIEWTB': [self._table_passer, self._table_passer],   # view elements

            #==================================
            'OUGATO2': [self._table_passer, self._table_passer],
            'OUGCRM2': [self._table_passer, self._table_passer],
            'OUGNO2': [self._table_passer, self._table_passer],
            'OUGPSD2': [self._table_passer, self._table_passer],  # psd
            'OUGRMS2': [self._table_passer, self._table_passer],  # rms

            'OQGATO2': [self._table_passer, self._table_passer],
            'OQGCRM2': [self._table_passer, self._table_passer],

            'OQGNO2': [self._table_passer, self._table_passer],
            'OQGPSD2': [self._table_passer, self._table_passer],
            'OQGRMS2': [self._table_passer, self._table_passer],

            'OFMPF2M': [self._table_passer, self._table_passer],
            'OLMPF2M': [self._table_passer, self._table_passer],
            'OPMPF2M': [self._table_passer, self._table_passer],
            'OSMPF2M': [self._table_passer, self._table_passer],
            'OGPMPF2M': [self._table_passer, self._table_passer],

            'OEFATO2': [self._table_passer, self._table_passer],
            'OEFCRM2': [self._table_passer, self._table_passer],
            'OEFNO2': [self._table_passer, self._table_passer],
            'OEFPSD2': [self._table_passer, self._table_passer],
            'OEFRMS2': [self._table_passer, self._table_passer],

            'OESATO2': [self._table_passer, self._table_passer],
            'OESCRM2': [self._table_passer, self._table_passer],
            'OESNO2': [self._table_passer, self._table_passer],
            'OESPSD2': [self._table_passer, self._table_passer],
            'OESRMS2': [self._table_passer, self._table_passer],

            'OVGATO2': [self._table_passer, self._table_passer],
            'OVGCRM2': [self._table_passer, self._table_passer],
            'OVGNO2': [self._table_passer, self._table_passer],
            'OVGPSD2': [self._table_passer, self._table_passer],
            'OVGRMS2': [self._table_passer, self._table_passer],

            #==================================
            #'GPL': [self._table_passer, self._table_passer],
            'OMM2': [self._table_passer, self._table_passer],
            'ERRORN': [self._table_passer, self._table_passer],
            #==================================

            'OCRPG': [self._table_passer, self._table_passer],
            'OCRUG': [self._table_passer, self._table_passer],

            'EDOM': [self._table_passer, self._table_passer],

            'OAGPSD2': [self._table_passer, self._table_passer],
            'OAGATO2': [self._table_passer, self._table_passer],
            'OAGRMS2': [self._table_passer, self._table_passer],
            'OAGNO2': [self._table_passer, self._table_passer],
            'OAGCRM2': [self._table_passer, self._table_passer],

            'OPGPSD2': [self._table_passer, self._table_passer],
            'OPGATO2': [self._table_passer, self._table_passer],
            'OPGRMS2': [self._table_passer, self._table_passer],
            'OPGNO2': [self._table_passer, self._table_passer],
            'OPGCRM2': [self._table_passer, self._table_passer],

            'OSTRPSD2': [self._table_passer, self._table_passer],
            'OSTRATO2': [self._table_passer, self._table_passer],
            'OSTRRMS2': [self._table_passer, self._table_passer],
            'OSTRNO2': [self._table_passer, self._table_passer],
            'OSTRCRM2': [self._table_passer, self._table_passer],

            'OQMPSD2': [self._table_passer, self._table_passer],
            'OQMATO2': [self._table_passer, self._table_passer],
            'OQMRMS2': [self._table_passer, self._table_passer],
            'OQMNO2': [self._table_passer, self._table_passer],
            'OQMCRM2': [self._table_passer, self._table_passer],

            'AAA': [self._table_passer, self._table_passer],
            'AAA': [self._table_passer, self._table_passer],
            'AAA': [self._table_passer, self._table_passer],
            'AAA': [self._table_passer, self._table_passer],
        }
        return table_mapper

class OP2Geom_Vectorized(OP2Geom, OP2_Vectorized):
    def __init__(self, debug=False, log=None, debug_file=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param debug: enables the debug log and sets the debug in the logger (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        :param debug_file: sets the filename that will be written to (default=None -> no debug)
        """
        # make_geom=False, debug=True, log=None, debug_file=None
        OP2_Vectorized.__init__(self, debug=debug, log=log, debug_file=debug_file)
        OP2Geom.__init__(self, debug=debug, log=log, debug_file=debug_file)
        self.make_geom = True
