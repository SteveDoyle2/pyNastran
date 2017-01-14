from __future__ import print_function
from six import iteritems
from six.moves import zip
from math import sqrt, log10

from collections import defaultdict
import numpy as np
from numpy import array
from numpy.linalg import norm

from pyNastran.op2.op2 import OP2


class MarginChecker(object):
    def __init__(self, op2_filenames=None, subcases=None, IDs=None):
        """
        Performs load case combination for:

        Parameters
        ----------
        op2_filenames : List[str]; default=['fem.op2']
            list of op2_filenames that the subcase result will come from
        subcases : List[int]; default=[1]
            list of subcases to grab from from each filename

        Assumptions:
            * linear static analysis
            * nodes numbers are consistent across different OP2s/subcases
            * only does static loading & max_deflection
            * solid elements stress ONLY

        ::
          User Info:
            UnitCase  Tension      Compression   Bending
            Filename  tension.op2  compBend.op2  compBend.op2
            isubcase  1            1             2

            unit_load 100.         1.            1.
            req_load  100.         100.          100.
            ratio     1.           100.          100.

          Code Inputs:
            op2_filenames =  ['tension.op2','comp_bend.op2','comp_bend.op2']
            subcases  =  [1,            1,              2,             ]
            # only for VonMises
            vm_factors = [[1.,           100.,          100.          ]]
            IDs        =  ['tension',    'comp',       'bend'         ]

        .. note:: vm_factors and IDs may have multiple levels...

        .. code :: python

             vm_factors = [[1.,100.,100.],
                           [2., 50., 75.],]
             case_names = ['Case1','Case2']  # only for VM``
        """
        if op2_filenames is None:
            op2_filenames = ['fem.op2']
        if subcases is None:
            subcases = [1]
        if IDs is None:
            IDs = [None]

        self.op2_filenames = op2_filenames
        self.subcases = subcases
        self.IDs = IDs
        assert len(self.op2_filenames) == len(self.subcases) == len(self.IDs)

        print("op2_filenames = ", self.op2_filenames)
        print("subcases  = ", self.subcases)
        print("IDs       = ", self.IDs)

        #self.load_factors = [[1.,2.,],  # fatigue
        #                     [4.,5.,]]

        self.op2s = {}
        self.displacement_results = {}
        self.solid_stress_results = {}
        self.plate_stress_results = {}

        #self.composite_plate-stress_results = {}
        #self.composite_plate_strain_results = {}

        self.cases = defaultdict(list)
        for icase, op2_filename, isubcase in zip(count(), self.op2_filenames, self.subcases):
            self.cases[op2_filename] = (icase, isubcase)

        for key, values in sorted(iteritems(self.cases)):
            print("case[%s] = %s" % (key, values))

    def read_files(self):
        for op2_filename, (icase, subcase_ids) in sorted(iteritems(self.cases)):
            subcase_ids = list(set(subcase_ids))
            print("subcase_ids[%s] = %s" % (op2_filename, subcase_ids))

            op2 = OP2(debug=False)
            op2.set_subcases(subcase_ids)
            results = [
                'deflection',
                'ctria3_plate_stress', 'cquad4_plate_stress',
                'ctetra_solid_stress', 'cpenta_solid_stress', 'chexa_solid_stress',
            ]
            op2.set_results(results)
            #op2.get_all_results()
            op2.read_op2(op2_filename)
            for subcase_id in subcase_ids:
                self.op2s[icase] = op2

    def check_deflections(self, max_deflection):
        """
        Parameters
        ----------
        max_deflection : float?
            [x,y,z,combined] abs values are used
        """
        self.max_deflection = max_deflection
        res = []
        for op2 in self.op2_filenames:
            deflection = op2.deflection[isubcase]
            data = deflection[0, :, :3] * factor
            res.append(data)
            nids = deflection.node_gridtype[:, 0]

        if len(res) == 1:
            res = data
        else:
            res = np.hstack(data)

        txyz = norm(res, axis=1)
        deflection_margin = max_deflection / txyz - 1
        return nids, deflection_margin

    def check_von_mises(self, vm_factors=None, case_names=None, Fty=100.):
        r"""
        currently only handles von mises stress for solid elements...

        Parameters
        ----------
        vm_factors
            see ``self.__init__``
        case_names : List[str]
            see ``self.__init__``
        Fty : float
            max allowable yeild stress (same units as FEM stress!)

        \f[ \sigma_v^2 = \tfrac{1}{2}[(\sigma_{11} - \sigma_{22})^2 + (\sigma_{22} - \sigma_{33})^2 + (\sigma_{11} - \sigma_{33})^2 + 6(\sigma_{23}^2 + \sigma_{31}^2 + \sigma_{12}^2)]  \f]
        \f[ \sigma_v   = \sqrt{\frac{(\sigma_1 - \sigma_2)^2 + (\sigma_2 - \sigma_3)^2 + (\sigma_1 - \sigma_3)^2 } {2}} \f]
        """
        if vm_factors is None:
            vm_factors = [[1.]]
        if case_names is None:
            case_names = ['Case1']
        self.vm_factors = vm_factors
        self.case_names = case_names
        self.Fty = Fty  # ksi
        assert len(self.case_names) == len(self.vm_factors)

        print("vm_factors = ", self.vm_factors)
        print("case_names = ", self.case_names)
        (stressP, eid_list) = self.process_solid_stress()
        (stressP, eid_list) = self.process_plate_stress()
        #print(stressP)
        self.stressP = stressP

        Fty = self.Fty  # ksi
        min_margins = {}
        for eid in eid_list:
            margins = []
            for icase, vm_factor in enumerate(self.vm_factors):
                margin = stressP[icase][eid] / Fty - 1
                margins.append(margin)
            min_margin = min(margins)
            min_margins[eid] = min_margin
            case = margins.index(min_margin)
            print("case=%-6s min_margin[%s] = %g" % (
                self.case_names[case], eid, min_margin))

    def process_plate_stress(self):
        r"""
        \f[ \sigma_v = \sqrt{\sigma_1^2- \sigma_1\sigma_2+ \sigma_2^2 + 3\sigma_{12}^2} \f]
        ovm^2 = o1^2 - o1*o2 + o2^2 + 3*o12^2
        """
        # ovm^2 = o1^2 - o1*o2 + o2^2 + 3*o12^2
        stress = []

        #op2_filenames =  ['tension.op2','comp_bend.op2','comp_bend.op2']
        #subcases  =  [1,            1,              2,             ]
        #vm_factors = [[1.,           100.,          100.          ]]
        #IDs        =  ['tension',    'comp',       'bend'         ]

        ovm = []
        #   0          1    2    3
        # [fiber_dist, oxx, oyy, txy, angle, omax, omin, ovm]
        isolt = [1, 2, 3]
        res_types = [
            'ctria3_plate_stress',
            'cquad4_plate_stress',
        ]
        for icase, vm_factor in enumerate(self.vm_factors):
            op2 = self.op2s[icase]
            datas = []
            eids = []
            for res_type in res_types:
                res = getattr(op2.res_type)
                if icase in res:
                    plate = res[icase]
                    data = plate[icase].data[0, :, islot] * vm_factor
                    eid = plate.element_node
                    datas.append(data)
                    eids.append(eid)

            if len(datas) == 1:
                data = datai
                eid = eid
            else:
                data = np.vstack(data)
                eid = np.vstack(eids)
            oxx = data[:, 0]
            oyy = data[:, 1]
            txy = data[:, 2]
            ovmi = np.sqrt(oxx**2 - oxx*oyy + oyy**2 + 3*txy**2)
            ovm.append(ovmi)

        if len(ovm) == 1:
            ovm = data
        else:
            ovm = np.hstack(ovm)

        return ovm, eids

    def process_solid_stress(self):  # stressP[icase][eid] = min(eid_results)
        ovm = []
        #   0          1    2    3
        # [fiber_dist, oxx, oyy, txy, angle, omax, omin, ovm]

        #headers = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, von_mises]
        isolt = [0, 1, 2, 3, 4, 5]
        res_types = [
            'ctetra_solid_stress',
            'cpenta_solid_stress',
            'chexa_solid_stress',
        ]
        for icase, vm_factor in enumerate(self.vm_factors):
            op2 = self.op2s[icase]
            datas = []
            eids = []
            for res_type in res_types:
                res = getattr(op2.res_type)
                if icase in res:
                    solid = res[icase]
                    data = solid[icase].data[0, :, islot] * vm_factor
                    eid = solid.element_node
                    datas.append(data)
                    eids.append(eid)

            if len(datas) == 1:
                data = datai
                eid = eid
            else:
                data = np.vstack(data)
                eid = np.vstack(eids)
            oxx = data[:, 0]
            oyy = data[:, 1]
            ozz = data[:, 2]
            txy = data[:, 3]
            tyz = data[:, 4]
            txz = data[:, 5]

            ovmi = np.sqrt((oxx - oyy) ** 2 + (oyy - ozz) ** 2 + (oxx - ozz) ** 2 + 6 * (txy ** 2 + tyz ** 2 + txz ** 2))  # 3d stress
            ovm.append(ovmi)

        if len(ovm) == 1:
            ovm = data
        else:
            ovm = np.hstack(ovm)
        return ovm, eids

    def rainflow(self):
        """required for fatigue"""
        pass

    def damage_count(self):
        """required for fatigue"""
        # 1.  get principal stresses at different loading combinations
        # 2.  rainflow count
        # 3.  find damage caused by each load case
        # 4.  sum damage (Miner's Rule)
        # 5.  find equivalent stress
        # 6.  get damage margin per ratio of life used
        pass

    def SNcurve(self, stress):
        """
        .. code-block:: python

          stress amplitude = fatigue strength coefficient* (2 N(f))^b
          N(f) is the cycles to failure
          2N(f) is the number of load reversals to failure
          b is the fatigue strength exponent
          For an AISI Type 1015 steel
          b = -0.11
          fatigue strength coefficient where 2N(f) =1,120 Ksi

          Log Sf = Log a + b Log N = Log (1.62 Sut) + Log N-0.0851
          a = (0.9 Sut)2 / Se = 1.62 Sut
          b = -(Log (0.9 Sut / Se)) / 3 (= -0.0851 when Se = .5 Sut)
          Se = 0.5 Ftu  = Endurance Limit = Stress corresponding to 'infinite' life of 1,000,000 or more cycles.
          Sf = Stress corresponding to a fatigue life, N, of 1000 to 1,000,000 cycles inclusive.
        """
        Sut = 120.   # ksi
        Se = 0.5 * Sut  # ksi
        a = (0.9 * Sut) * 2 / Se
        b = -(log10(0.9 * Sut / Se)) / 3.
        #logSf = log10(a) + b * log10(N)
        #logN = (logSf - log10(a) )/b
        ncycles = 10. ** logN
        return ncycles


def main():
    op2_filenames = ['fem.op2', 'fem.op2']
    subcases = [1, 2]
    Fty = 150.  # ksi
    IDs = ['tension', 'compression']

    vm_factors = [
        # [1.,0.,],  # tension
        #[0.,1.,],  # bending
        #[0.,-1.,], # minus bending
        #[-1.,0.],  # compression
        #[-1.,1.],  # compression+bend
        #[1.,1.,],  # tension + bending
        #[-1.,-1.], # compression+mbend
        [1., -1., ], ]  # tension + mbending

    #caseNames = ['tens','bend','mBend','comp','cBend','tBend','cmBend','tmBend']
    case_names = ['tmBend']
    analysis = MarginChecker(op2_filenames, subcases, IDs)
    analysis.read_files()
    analysis.check_von_mises(vm_factors, case_names, Fty)
    #a.checkDeflections(0.018)

if __name__ == '__main__':  # pragma: no cover
    main()
