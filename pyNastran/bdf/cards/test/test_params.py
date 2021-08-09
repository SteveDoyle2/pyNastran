import unittest

#import pyNastran
from pyNastran.bdf.bdf import BDF, PARAM, BDFCard, MDLPRM
from pyNastran.bdf.cards.test.utils import save_load_deck

#ROOT_PATH = pyNastran.__path__[0]
#test_path = os.path.join(ROOT_PATH, 'bdf', 'cards', 'test')

class TestParams(unittest.TestCase):
    def test_mdlprm_1(self):
        model = BDF(debug=True, log=None, mode='msc')
        mdlprm_dict = {'HDF5': 1}
        model.add_mdlprm(mdlprm_dict, comment='mdlprm')
        save_load_deck(model, xref='standard', op2_log_level='warning')

    def test_param_01(self):
        card = PARAM('NOCOMP', [-1])
        #print('%r' % card)
        assert str(card) == 'PARAM     NOCOMP      -1\n', '%r' % str(card)

        cardi = BDFCard(['PARAM', 'NOCOMP', -1])
        card = PARAM.add_card(cardi)
        assert str(card) == 'PARAM     NOCOMP      -1\n', '%r' % str(card)
        #print('%r' % card)

    def test_param_02(self):
        model = BDF(debug=False)
        param = model.add_param('POST', -1, comment='param post')
        param.raw_fields()
        param.write_card(size=8, is_double=False)
        param.write_card(size=16, is_double=False)
        param.write_card(size=16, is_double=True)
        param.write_card_16(is_double=False)
        param.write_card_16(is_double=True)

    def test_param_mystran(self):
        model = BDF(debug=True, log=None, mode='mystran')
        param1 = model.add_param('POST', -1, comment='')
        str(param1)

        fields = ['PARAM', 'EQCHECK', '0', 3, 3, None, None, None, -1.E10]
        param2 = model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)
        str(param2)

        fields = ['PARAM', 'AUTOSPC', 'Y', 1.0E-9, None, 'Y', 'Y']
        param3 = model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)
        str(param3)

        fields = ['PARAM', 'ART_MASS', 'Y', '0.', 1.E-6]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'EIGNORM2', 'Y']
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'GRDPNT', 0]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'MEFMLOC', 'CG']
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'RCONDK', 'Y']
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'SOLLIB', 'LAPACK']
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'WTMASS', .002591]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        # basic flags
        fields = ['PARAM', 'PRTBASIC', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTCGLTM', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTDLR', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTGMN', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTGOA', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTSCP', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTTSET', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTHMN', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTIFLTM', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTKXX', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTMXX', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTPHIXA', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTPHIZL', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTQSYS', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTUO0', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTYS', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        # not basic flags
        fields = ['PARAM', 'PRTCORD', '2']
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTDOF', 1]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTRMG', 3]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        # fancy flags
        fields = ['PARAM', 'PRTSTIFF', 1, 3, 3, 3, 3]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTSTIFD', 1, 3, 3, 3, 3]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTFOR', 1, 3, 3, 3, 3]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)

        fields = ['PARAM', 'PRTMASS', 1, 3, 3, 3, 3]
        model.add_card(fields, 'PARAM', comment='param', ifile=None, is_list=True, has_none=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
