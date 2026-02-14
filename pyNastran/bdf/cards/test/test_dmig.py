import unittest
import os
import copy
import numpy as np

from cpylog import SimpleLogger
import pyNastran
from pyNastran.bdf.bdf import BDF, BDFCard, read_bdf, DMI, DMIG, fill_dmigs
from pyNastran.bdf.cards.test.utils import save_load_deck, get_matrices
from pyNastran.bdf.bdf_interface.utils import fill_dmigs


PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'cards', 'test')

class TestDTI(unittest.TestCase):

    def test_dti_units(self):
        """tests DTI,UNITS"""
        model = BDF(debug=False)
        fields = {
            'mass' : 'mass',
            'length' : 'length',
            'force' : 'force',
            'time' : 'time',
            'temp_stress' : 'temp_str',
        }
        dti = model.add_dti('UNITS', fields, comment='dti,units')
        dti2 = copy.deepcopy(dti)
        dti.raw_fields()
        #print(dti.write_card())
        save_load_deck(model)

    def test_dti_units2(self):
        """tests DTI,UNITS"""
        model = BDF(debug=False)
        fields = {
            'mass' : 'mass',
            'length' : 'length',
            'force' : 'force',
            'time' : 'time',
            'temp_stress' : None,
        }
        dti = model.add_dti('UNITS', fields, comment='dti,units')
        dti2 = copy.deepcopy(dti)
        dti.raw_fields()
        #print(dti.write_card())
        save_load_deck(model)

class TestMatrix(unittest.TestCase):
    def test_matrix_real(self):
        from pyNastran.op2.result_objects.matrix import Matrix
        eye = np.eye(3)
        mat = Matrix('EYER', 'square', eye)
        mat.write_dmi()

    def test_matrix_complex(self):
        from pyNastran.op2.result_objects.matrix import Matrix
        eye = np.eye(3) + np.eye(3) * 2j
        mat = Matrix('EYEC', 'square', eye)
        mat.write_dmi()

class TestDMI(unittest.TestCase):

    def test_dmi_square_mult(self):
        model = BDF(debug=True)
        name = 'WKK'
        form = 'square'
        tin = 'float32'
        tout = tin
        nrows = 4
        ncols = 4
        GCi = [1, 2, 2, 3, 3, 4, 4]
        GCj = [1, 1, 2, 2, 3, 3, 4]
        # Real = np.random.random(len(GCi))
        Real = np.ones(len(GCi))
        dmi1 = model.add_dmi(
            name, form, tin, nrows, ncols,
            GCi, GCj, Real, tout=tout)
        data1, rows1, cols1 = dmi1.get_matrix(
            is_sparse=False, apply_symmetry=False)
        # print(dmi1)
        GCid = [1, 2, 3, 4]
        GCjd = [1, 2, 3, 4]
        Reald = [3., 2., 7., 5.]
        dmi2 = model.add_dmi(
            'MYDIAG', 'square', tin, nrows, ncols,
            GCid, GCjd, Reald, tout=tout,
        )
        diag_data2, rows2, cols2 = dmi2.get_matrix(
            is_sparse=False, apply_symmetry=False)
        #print()
        #print(data1)
        #print(diag_data2)
        out = diag_data2 @ data1
        #print(out)

        updated_dmi1 = model.add_dense_dmi(
            'TEMP', out,
            'square', tin, tout=tout,
        )
        updated_dmi1.name = dmi1.name



    def test_dmi_alternate(self):
        model = BDF(debug=True)
        lines1 = ['DMI,ALTERNTE,0,3,1,1,,12,1']
        model.add_card(lines1, 'DMI', is_list=False)
        lines2 = [
            'DMI', 'ALTERNTE', '1', '2', '1.0', '0.0',
            '2.0', '0.0', '3.0', '0.0', '4.0', '0.0', '5.0',
            '0.0', '6.0',
        ]
        model.add_card(lines2, 'DMI')
        fill_dmigs(model)
        #print(model.dmi)
        dmi = model.dmi['ALTERNTE']
        alternate = dmi.get_matrix()[0]
        assert alternate.shape == (12, 1), alternate.shape
        alternate = alternate.flatten()
        alternate_expected = np.array([
            0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0, 0.0, 6.0])
        assert np.array_equal(alternate, alternate_expected), f'alternate={alternate};\nalternate_expected={alternate_expected}'
        print(str(dmi))

    def test_dmi_thru_real(self):
        model = BDF(debug=True)
        lines1 = ['DMI,RRR,0,3,1,1,,12,1']
        model.add_card(lines1, 'DMI', is_list=False)
        lines2 = ['DMI,RRR,1,2,1.0,THRU,10,12,2.0']
        model.add_card(lines2, 'DMI', is_list=False)
        fill_dmigs(model)
        #print(model.dmi)
        dmi = model.dmi['RRR']
        rrr = dmi.get_matrix()[0]
        assert rrr.shape == (12, 1), rrr.shape
        rrr = rrr.flatten()
        rrr_expected = np.array([0.] + [1.] * 9 + [0., 2.])
        assert np.array_equal(rrr, rrr_expected), f'rrr={rrr};\nrrr_expected={rrr_expected}'
        print(str(dmi))

    def test_dmi_thru_complex(self):
        model = BDF(debug=True)
        lines1 = ['DMI,RRR,0,3,3,3,,12,1']
        model.add_card(lines1, 'DMI', is_list=False)
        lines2 = ['DMI,RRR,1,2,1.0,3.0,THRU,10,12,',
                  ',2.0,4.0']
        model.add_card(lines2, 'DMI', is_list=False)
        fill_dmigs(model)
        #print(model.dmi)
        dmi = model.dmi['RRR']
        rrr = dmi.get_matrix()[0]
        assert rrr.shape == (12, 1), rrr.shape
        rrr = rrr.flatten()
        rrr_expected_real = np.array([0.] + [1.] * 9 + [0., 2.])
        rrr_expected_imag = np.array([0.] + [3.] * 9 + [0., 4.])
        #print(rrr.real)
        #print(rrr.imag)
        assert np.array_equal(rrr.real, rrr_expected_real), f'rrr={rrr.real};\nrrr_expected={rrr_expected_real}'
        assert np.array_equal(rrr.imag, rrr_expected_imag), f'rrr={rrr.imag};\nrrr_expected={rrr_expected_imag}'
        print(str(dmi))

    def test_dmi_01(self):
        """tests a DMI card"""
        lines = ['DMI,Q,0,6,1,0,,4,4']
        model = BDF(debug=False)
        card = model._process_card(lines)
        card_obj = BDFCard(card)

        size = 8
        dmi = DMI.add_card(card_obj)
        dmi.write_card(size, 'dummy')
        #dmi.raw_fields()

    def test_dmi_w2gj_1(self):
        data = """
DMI         W2GJ       0       2       1       0            1200       1
DMI         W2GJ       1       1 1.54685.1353939.1312423.0986108.0621382
        .0369645.0257369.0234453.0255875.05652071.561626.1205361.1125278
        .0846353.0557613.0355958.0250237.0222578.0246823.05448051.532335
        .1076103.0988449.0754307.0525759.0348965.0245642.0227323.0263744
        .05396551.515989 .095201.0867341.0672852.0484031  .03261 .022431
        .0212566.0260411.05388151.487622.0829165.0753736.0586076.0422512
        .0278508.0174376.0167388 .023505.05292421.448751 .070779.0644553
        .0497634.0353554.0217712.0110981 .010631.0189924.0513639 1.25276
        .0588486.0540163.0412702.0286076.0158065.0047222.0041635.0132007
         .050093.6120345.0477668.0442399.0331648 .022064 .009834-1.828-3
        -2.503-3.0071487.0496528-.041134.0378177.0352355.0256046.0157419
        .0039734-8.065-3-8.842-3.0020994.05098221.446484.0286787.0268212
        .0183474.0095515-1.486-3 -.01338-.013646-6.974-4 .0543671.548042
        .0209779.0194154.0117029.0038052-6.315-3-.017763 -.01672-1.465-3
        .0581247-.993292.0142634 .012689.0053179-2.231-3-.012116-.022461
        -.020295-3.444-3.0594645-1.02737.0076634.0058206-1.478-3-9.128-3
         -.01927 -.02811-.025353-7.069-3.0581478-1.06258.0013462-7.586-4
        -7.972-3-.015515-.025103-.033856-.030101-.010728  .05707-1.08126
        -4.495-3-6.589-3-.013483-.020721-.029771-.037558-.032107-.012514
        .0580607-1.09659-9.834-3-.011544 -.01795-.024825-.033374 -.04009
        -.033981-.013454.05842491.327594-.014595-.015787 -.02156-.027993
        -.036203 -.04211-.035644-.014236.0587609-1.12656-.018867  -.0194
        -.024477-.030468-.038234-.043564-.036789-.014931.05923291.406868
        -.022712-.022633-.026994-.032539-.039873-.044645-.037676-.015477
        .0600661 -1.1403-.026208-.025684-.029346-.034437-.041619-.046053
        -.038965-.016419.0602822-1.15545-.029382-.028471 -.03156 -.03624
         -.04336-.048018-.040796-.018028.0592796-1.16853-.032213-.030919
        -.033614-.038048-.045083-.049868-.042604-.019698.0580287-1.17282
        -.034734-.032994-.035557-.039919-.046795-.051645 -.04442-.021465
        .0565324-1.18209-.036935-.034901-.037274-.041597 -.04831-.053297
        -.046237-.023328.05519341.553182-.038891-.036845-.038764-.042911
        -.049524-.054687-.048015-.025159.0543401-1.17264-.040588  -.0387
        -.040155 -.04408-.050704 -.05603-.049735-.027008.0538784 1.55343
        -.041944-.040277-.041628-.045311-.052116-.057595-.051508-.028825
        .0533671-1.17067-.043074-.041669-.043193 -.04673-.053791-.059407
        -.053284-.030358.0528413 -1.1734-.044205-.043085-.044977 -.04865
        -.055744 -.06141 -.05495-.031233.0523855-1.16752-.045483-.044569
        -.047071-.051162 -.05808-.063396-.056434-.031236.0518605-1.16843
         -.04698-.046239-.049258-.053844-.060565-.065508-.057742-.031123
         .051474-1.15854-.048786-.048225-.051408-.056377-.063191-.067552
        -.059057-.031243.0501356-1.12756-.050892-.050449-.053679-.058932
        -.065748-.069427-.060094 -.03156.0482373-1.11765-.053268-.052849
        -.056289-.061745-.068092-.070729-.060566-.031679.0464637-1.08369
        -.055984-.055508-.059297-.064792-.070378-.071549-.060739-.031831
        .0438902-1.03624 -.05879-.058529-.062614-.067938 -.07253-.072523
        -.061477-.033638.03692041.523722-.062845-.062397-.066421-.071335
        -.074786-.073668-.062967-.036846.0276397 -.11955-.066722-.066324
        -.070425-.074913-.077237-.075171 -.06493-.041043.0165293.9973973
        -.062269-.064098-.070488  -.0769-.080234-.079264-.071769-.051909
         -.00519.5332272-.043435-.050199 -.06278-.075336-.083821-.088398
        -.088075-.075685-.044054     601 1.54685.1353939.1312423.0986108
        .0621382.0369645.0257369.0234453.0255875.05652071.561626.1205361
        .1125278.0846353.0557613.0355958.0250237.0222578.0246823.0544805
        1.532335.1076103.0988449.0754307.0525759.0348965.0245642.0227323
        .0263744.05396551.515989 .095201.0867341.0672852.0484031  .03261
         .022431.0212566.0260411.05388151.487622.0829165.0753736.0586076
        .0422512.0278508.0174376.0167388 .023505.05292421.448751 .070779
        .0644553.0497634.0353554.0217712.0110981 .010631.0189924.0513639
         1.25276.0588486.0540163.0412702.0286076.0158065.0047222.0041635
        .0132007 .050093.6120345.0477668.0442399.0331648 .022064 .009834
        -1.828-3-2.503-3.0071487.0496528-.041134.0378177.0352355.0256046
        .0157419.0039734-8.065-3-8.842-3.0020994.05098221.446484.0286787
        .0268212.0183474.0095515-1.486-3 -.01338-.013646-6.974-4 .054367
        1.548042.0209779.0194154.0117029.0038052-6.315-3-.017763 -.01672
        -1.465-3.0581247-.993292.0142634 .012689.0053179-2.231-3-.012116
        -.022461-.020295-3.444-3.0594645-1.02737.0076634.0058206-1.478-3
        -9.128-3 -.01927 -.02811-.025353-7.069-3.0581478-1.06258.0013462
        -7.586-4-7.972-3-.015515-.025103-.033856-.030101-.010728  .05707
        -1.08126-4.495-3-6.589-3-.013483-.020721-.029771-.037558-.032107
        -.012514.0580607-1.09659-9.834-3-.011544 -.01795-.024825-.033374
         -.04009-.033981-.013454.05842491.327594-.014595-.015787 -.02156
        -.027993-.036203 -.04211-.035644-.014236.0587609-1.12656-.018867
          -.0194-.024477-.030468-.038234-.043564-.036789-.014931.0592329
        1.406868-.022712-.022633-.026994-.032539-.039873-.044645-.037676
        -.015477.0600661 -1.1403-.026208-.025684-.029346-.034437-.041619
        -.046053-.038965-.016419.0602822-1.15545-.029382-.028471 -.03156
         -.03624 -.04336-.048018-.040796-.018028.0592796-1.16853-.032213
        -.030919-.033614-.038048-.045083-.049868-.042604-.019698.0580287
        -1.17282-.034734-.032994-.035557-.039919-.046795-.051645 -.04442
        -.021465.0565324-1.18209-.036935-.034901-.037274-.041597 -.04831
        -.053297-.046237-.023328.05519341.553182-.038891-.036845-.038764
        -.042911-.049524-.054687-.048015-.025159.0543401-1.17264-.040588
          -.0387-.040155 -.04408-.050704 -.05603-.049735-.027008.0538784
         1.55343-.041944-.040277-.041628-.045311-.052116-.057595-.051508
        -.028825.0533671-1.17067-.043074-.041669-.043193 -.04673-.053791
        -.059407-.053284-.030358.0528413 -1.1734-.044205-.043085-.044977
         -.04865-.055744 -.06141 -.05495-.031233.0523855-1.16752-.045483
        -.044569-.047071-.051162 -.05808-.063396-.056434-.031236.0518605
        -1.16843 -.04698-.046239-.049258-.053844-.060565-.065508-.057742
        -.031123 .051474-1.15854-.048786-.048225-.051408-.056377-.063191
        -.067552-.059057-.031243.0501356-1.12756-.050892-.050449-.053679
        -.058932-.065748-.069427-.060094 -.03156.0482373-1.11765-.053268
        -.052849-.056289-.061745-.068092-.070729-.060566-.031679.0464637
        -1.08369-.055984-.055508-.059297-.064792-.070378-.071549-.060739
        -.031831.0438902-1.03624 -.05879-.058529-.062614-.067938 -.07253
        -.072523-.061477-.033638.03692041.523722-.062845-.062397-.066421
        -.071335-.074786-.073668-.062967-.036846.0276397 -.11955-.066722
        -.066324-.070425-.074913-.077237-.075171 -.06493-.041043.0165293
        .9973973-.062269-.064098-.070488  -.0769-.080234-.079264-.071769
        -.051909 -.00519.5332272-.043435-.050199 -.06278-.075336-.083821
        -.088398-.088075-.075685-.044054
        """
        with open('dmi.bdf', 'w') as bdf_file:
            bdf_file.write(data)
        model = BDF(debug=False)
        model.read_bdf('dmi.bdf', punch=True)
        w2gj = model.dmi['W2GJ']
        assert w2gj.shape == (1200, 1), w2gj.shape
        w2gj.get_matrix()

        real2 = []
        for i, unused_real in enumerate(w2gj.Real):
            real2.append(0.1  * i)
        #w2gj.Real = real2
        #print(w2gj.GCi)  # varying (rows)
        #print(w2gj.GCj)  # constant (cols)

        model.write_bdf('dmi_out.bdf')

        model2 = BDF(debug=False)
        model2.read_bdf('dmi_out.bdf')
        w2gj_new = model.dmi['W2GJ']
        assert w2gj_new.shape == (1200, 1), w2gj_new.shape

        assert np.allclose(w2gj.GCi, w2gj_new.GCi)
        assert np.allclose(w2gj.GCj, w2gj_new.GCj)
        assert np.allclose(w2gj.Real, w2gj_new.Real)
        os.remove('dmi.bdf')
        os.remove('dmi_out.bdf')
        save_load_deck(model2)

    def test_dmi_wkk_diagonal(self):
        log = SimpleLogger(level='debug')
        model = BDF(log=log)
        nrows = 100
        GCi = np.arange(1, nrows+1, dtype='int32')
        GCj = GCi.copy()
        form = 'square'
        Real = 0.5 * np.ones(nrows, dtype='float32')
        dmi = model.add_dmi(
            name='WKK', form=form,
            tin='float32', tout='float32',
            nrows=nrows, ncols=nrows,
            GCj=GCj, GCi=GCi,
            Real=Real, Complex=None, comment='wkk')
        # print(dmi)
        model_out = save_load_deck(
            model, run_op2_writer=False, run_convert=False,
            run_mirror=False
        )
        dmi = model_out.dmi['WKK']
        str(dmi)

    def test_dmi_complex(self):
        """tests a complex DMI"""
        #DMI QQQ 0 2 3 3 4 2
        #DMI QQQ 1 1 1.0 2.0 3.0 0.0 3
        #5.0 6.0
        # DMI QQQ 2 2 6.0 7.0 4 8.0 9.0
        #[QQQ]
        #1.0 + 2.0i , 0.0 + 0.0i
        #3.0 + 0.0i , 6.0 + 7.0i
        #5.0 + 6.0i , 0.0 + 0.0i
        #0.0 + 0.0i , 8.0 + 9.0i
        j = [1, 1, 1, 1,
             2, 2, 2, 2]  # col
        i = [1, 2, 3, 4,
             1, 2, 3, 4] # rows
        reals = [1., 3., 5., 0.,
                 0., 6., 0., 8.]
        complexs = [2., 0., 6., 0.,
                    0., 7., 0., 9.]
        model = BDF(debug=True, log=None, mode='msc')
        form = 2 # square
        nrows = 4
        ncols = 2

        name = 'QQQR'
        tin = 1 # real
        tout = 1 # real
        dmi_real = model.add_dmi(
            name, form, tin, nrows, ncols,
            j, i, reals, Complex=None, tout=tout, comment='dmi_real')
        str(dmi_real)

        name = 'QQQI'
        tin = 3 # complex
        tout = 3 # complex
        dmi_imag = model.add_dmi(
            name, form, tin, nrows, ncols,
            j, i, reals, Complex=complexs, tout=tout,
            comment='dmi_complex')
        str(dmi_imag)

        matrix1r, unused_blank, unused_blank = dmi_real.get_matrix(is_sparse=False, apply_symmetry=False)
        matrix2r, unused_blank, unused_blank = dmi_imag.get_matrix(is_sparse=False, apply_symmetry=False)

        matrix1s, unused_blank, unused_blank = dmi_real.get_matrix(is_sparse=True, apply_symmetry=False)
        matrix2s, unused_blank, unused_blank = dmi_imag.get_matrix(is_sparse=True, apply_symmetry=False)
        str(matrix1r)
        str(matrix1s)
        str(matrix2r)
        str(matrix2s)
        #print(matrix1r)
        #print(matrix1s)


class TestDMIGReal(unittest.TestCase):
    def test_dmig_1(self):
        """Tests DMIG reading"""
        model = BDF(debug=False)
        bdf_name = os.path.join(TEST_PATH, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)
        out = model.dmig['REALS'].get_matrix(is_sparse=False)

        reals_actual, unused_rows_reversed, unused_cols_reversed = out
        #print "---reals_actual---\n", reals_actual
        #print "---out---\n", out

        reals_expected = [
            [1.0, 0.5, 0.25],
            [0.5, 2.0, 0.75],
            [0.25, 0.75, 3.0],
        ]
        a_matrix = model.dmig['REALS']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)
        self.assertTrue(np.allclose(reals_expected, reals_actual))
        a_matrix.get_matrix()

    def test_dmig_2(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(TEST_PATH, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)

        out = model.dmig['REAL'].get_matrix(is_sparse=False)
        real_actual, unused_rows_reversed, unused_cols_reversed = out
        #print "---REAL_actual---\n", REAL_actual
        real_expected = [
            [1.0, 0.5, 0.25],
            [0.0, 2.0, 0.75],
            [0.0, 0.0, 3.0],
        ]
        a_matrix = model.dmig['REAL']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        self.assertTrue(np.allclose(real_expected, real_actual))
        a_matrix.get_matrix()

        #model2 = BDF(debug=False)
        #bdf_name = os.path.join(TEST_PATH, 'include_dir', 'include.inc')
        #model2.read_bdf(bdf_name, xref=False, punch=True)
        save_load_deck(model)

    def test_dmig_3(self):
        model = BDF(debug=False)
        bdf_filename = os.path.join(TEST_PATH, 'dmig.bdf')
        model.read_bdf(bdf_filename, xref=False, punch=True)
        out = model.dmig['POLE'].get_matrix(is_sparse=False)

        pole_actual, unused_rows_reversed, unused_cols_reversed = out
        #print("---pole_actual---\n", pole_actual)
        mag_expected = np.array([
            [1.0, 4.0, 5.0],
            [0.0, 2.0, 6.0],
            [0.0, 0.0, 3.0],
        ])

        a_matrix = model.dmig['POLE']
        #print('GCi = ', a_matrix.GCi)
        #print('GCj = ', a_matrix.GCj)
        #print('Real = ', a_matrix.Real)
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        A_expected = mag_expected * np.cos(np.radians(45))
        B_expected = mag_expected * np.sin(np.radians(45))
        pole_expected = A_expected + B_expected * 1j

        np.set_printoptions(precision=20)
        msg = '\n%s_actual\n%s\n\n----' % ('POLE', pole_actual)
        msg += '\n%s_expected\n%s\n----' % ('POLE', pole_expected)
        msg += '\n%s_delta\n%s\n----' % ('POLE', pole_actual - pole_expected)
        #print(msg)
        self.assertTrue(np.allclose(pole_expected, pole_actual), msg)
        np.set_printoptions(precision=6)
        a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_4(self):
        lines = ['DMIG    ENFORCE 0       1       1       0']
        model = BDF(debug=False)
        card = model._process_card(lines)
        card_obj = BDFCard(card)

        size = 8
        dmi = DMIG.add_card(card_obj)
        dmi.write_card(size, 'dummy')
        #dmi.raw_fields()

    def test_dmig_5(self):
        cards = [
            ['DMIG, A, 0, 9, 1, 1,  ,    , 1'],
            ['DMIG, A, 1, 0,  , 2, 1, 1.0,'],
            ['DMIG, A, 1, 0,  , 2, 2, 1.0,'],
            ['DMIG, A, 1, 0,  , 2, 3, 1.0,'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        fill_dmigs(model)

        a_matrix = model.dmig['A']
        assert len(a_matrix.GCi) == 3, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 3, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)
        #a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_6(self):
        cards = [
            ['DMIG, A, 1, 0,  , 2, 1, 1.0,'],
            ['DMIG, A, 0, 9, 1, 1,  ,    , 1'],
            ['DMIG, A, 1, 0,  , 2, 2, 1.0,'],
            ['DMIG, A, 1, 0,  , 2, 3, 1.0,'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        fill_dmigs(model)
        get_matrices(model)

        a_matrix = model.dmig['A']
        assert len(a_matrix.GCi) == 3, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 3, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)
        #a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_7(self):
        cards = [
            ['DMIG, A, 0, 9, 1, 1,  ,    , 1'],
            ['DMIG, A, 1, ,  , 2, 1, 1.0,'],
            ['DMIG, A, 1, ,  , 2, 2, 1.0,'],
            ['DMIG, A, 1, ,  , 2, 3, 1.0,'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        fill_dmigs(model)
        get_matrices(model)

        a_matrix = model.dmig['A']
        assert len(a_matrix.GCi) == 3, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 3, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)
        #a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_8(self):
        """symmetric"""
        cards = [
            ['DMIG,AMTRXX,0,6,1,0'],
            ['DMIG,AMTRXX,2,1, ,2,1,201.0, ,+DM1',
             '+DM1,2,3,203.0'],
            ['DMIG,AMTRXX,3,1, ,3,1,301.0, ,+DM2',
             '+DM2,3,3,303.0'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        fill_dmigs(model)
        a_matrix = model.dmig['AMTRXX']
        assert len(a_matrix.GCi) == 4, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 4, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)
        assert a_matrix.shape == (4, 4), 'shape=%s' % str(a_matrix.shape)
        a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_9(self):
        pch_filename = os.path.join(TEST_PATH, 'dmig.pch')
        model = read_bdf(pch_filename, debug=False, punch=True)
        vax = model.dmig['VAX']
        vax_array, vax_dict_row, vax_dict_col = vax.get_matrix()
        assert vax_array.shape == (15, 1), vax_array
        vax_dict_row_expected = {
            0: (101, 1), 1: (102, 1), 2: (105, 1), 3: (106, 1), 4: (107, 1), 5: (108, 1),
            6: (109, 1), 7: (201, 1), 8: (301, 1), 9: (302, 1), 10: (305, 1), 11: (306, 1),
            12: (307, 1), 13: (308, 1), 14: (309, 1)
        }
        vax_dict_col_expected = {0: (1, 0)}
        assert list(sorted(vax_dict_col)) == list(sorted(vax_dict_col_expected)), 'vax_dict_col=%s vax_dict_col_expected=%s' % (vax_dict_col, vax_dict_col_expected)
        assert list(sorted(vax_dict_row)) == list(sorted(vax_dict_row_expected)), 'vax_dict_row=%s vax_dict_row_expected=%s' % (vax_dict_row, vax_dict_row_expected)
        #save_load_deck(model)

    def test_dmig_10(self):
        """tests the add card method with a real DMIG"""
        model = BDF(debug=False)
        name = 'DMIG_1'
        matrix_form = 6
        tin = 1
        ncols = None
        reals = [1.0, 2.0, 3.0]
        GCj = [[1, 1],  # grid, component
               [2, 1],
               [3, 1]]
        GCi = [[1, 1],  # grid, component
               [4, 1],
               [5, 1]]
        dmig = model.add_dmig(name, matrix_form, tin, ncols, GCj, GCi,
                              Real=reals, Complex=None,
                              comment='dmig')
        assert dmig.is_real is True, dmig.is_real
        assert dmig.is_complex is False, dmig.is_complex
        assert dmig.is_polar is False, dmig.is_polar
        dmig.get_matrix()

        name = 'DMIK_1'
        nrows = None
        dmik = model.add_dmik(name, matrix_form, tin, ncols, GCj, GCi,
                              Real=reals, Complex=None,
                              comment='dmik')
        dmik.get_matrix()

        dmiji = model.add_dmiji(name, matrix_form, tin, nrows, ncols, GCj, GCi,
                                reals, Complex=None, comment='dmiji')
        dmiji.get_matrix()

        #dmi = model.add_dmi(name, matrix_form, tin, tout, nrows, ncols, GCj, GCi,
                            #reals, Complex=None, comment='dmi')

        save_load_deck(model, run_test_bdf=False)

    def test_dmig_sparse(self):
        """tests dmig_sparse.pch"""
        bdf_filename = os.path.join(TEST_PATH, 'dmig_sparse.pch')
        model = read_bdf(bdf_filename, validate=True, xref=True, punch=True,
                         save_file_structure=False, skip_cards=None, read_cards=None,
                         encoding=None, log=None, debug=False, mode='msc')
        get_matrices(model)
        kaax = model.dmig['KAAX']
        #print(kaax.get_stats())
        get_matrices(model)

        # make the matrix complex to test that
        kaax.tin = 3   # complex
        kaax.tout = 3  # complex
        kaax.Complex = kaax.Real
        get_matrices(model)
        #kaax = model.dmig['KAAX'].get_matrix(is_sparse=True)

    def test_dmig_sparse2(self):
        """tests dmig_sparse.pch"""
        model = BDF(debug=True, log=None, mode='msc')
        name = 'TEST'
        ifo = 6  # symmetic
        tin = 1  # float32
        tout = 1 # float32
        polar = 0
        ncols = 2
        GCi = [ [1, 1],]  # rows
        GCj = [[1, 2], ]  # cols
        Real = [1.]
        model.add_dmig(name, ifo, tin, ncols, GCj, GCi, Real, Complex=None,
                       tout=tout, polar=polar, comment='')
        get_matrices(model)
        kaax = model.dmig[name]
        dense = kaax.get_matrix(is_sparse=False, apply_symmetry=True)[0]
        dense1 = kaax.get_matrix(is_sparse=True, apply_symmetry=True)[0].toarray()
        dense2 = kaax.get_matrix(is_sparse=True, apply_symmetry=True)[0].toarray()
        dense_expected = [
            [0., 1.],
            [1., 0.],
        ]
        assert np.array_equal(dense, dense_expected)
        assert np.array_equal(dense1, dense_expected)
        assert np.array_equal(dense2, dense_expected)
        #print(kaax.get_stats())
        get_matrices(model)

        # make the matrix complex to test that
        kaax.tin = 3   # complex
        kaax.tout = 3  # complex
        kaax.Complex = kaax.Real
        get_matrices(model)
        #kaax = model.dmig['KAAX'].get_matrix(is_sparse=True)

    def test_dmig_symmetry(self):
        """testing symmetric DMIGs"""
        model = BDF(debug=False, log=None, mode='msc')

        tin = 1 # real, float32
        tout = 1 # real, float32

        # [
        #   1., 2., 3.,
        #       4., 5.,
        #           6.,
        # ]
        GCi = [
            [1, 1], [1, 2], [1, 3],
            [1, 2], [1, 3],
            [1, 3],
            [1, 3],
        ]
        GCj = [
            [1, 1], [1, 1], [1, 1],
            [1, 2], [1, 2],
            [1, 3],
            [1, 3],
        ]
        polar = 0
        Real = np.array([1., 2., 3., 4., 5., 6., 7.])
        ncols = 0  # ????
        ifo = 6 # symmetric

        name = 'TEST'
        model.add_dmig(name, ifo, tin, ncols,
                       GCj, GCi, Real, Complex=None,
                       tout=tout, polar=polar, comment='')
        get_matrices(model)
        test = model.dmig[name]
        A1, _, _ = test.get_matrix(is_sparse=False, apply_symmetry=True)
        A2, _, _ = test.get_matrix(is_sparse=True, apply_symmetry=True)

        A_expected = np.array([
            [1., 2., 3.],
            [2., 4., 5.],
            [3., 5., 13.],
        ])
        assert np.array_equal(A1, A_expected)
        assert np.array_equal(A2.toarray(), A_expected)

        test.tin = 3   # complex
        test.tout = 3  # complex
        test.Complex = 10 * test.Real
        B1, _, _ = test.get_matrix(is_sparse=False, apply_symmetry=True)
        B2, _, _ = test.get_matrix(is_sparse=True, apply_symmetry=True)
        B_expected = A_expected + 10j * A_expected
        assert np.array_equal(B1, B_expected)
        assert np.array_equal(B2.toarray(), B_expected)

        get_matrices(model)
        #kaax = model.dmig['KAAX'].get_matrix(is_sparse=True)

    def test_dmig_rectangular(self):
        """testing symmetric DMIGs"""
        model = BDF(debug=False, log=None, mode='msc')

        tin = 1 # real, float32
        tout = 1 # real, float32

        # [
        #   1., 2., 3., 0.
        #       4., 5., 0.
        #           13., 8.
        # ]
        GCj = [ # cols
            [1, 1], [1, 2], [1, 3],
            [1, 2], [1, 3],
            [1, 3], [1, 3],
            [1, 4],
        ]
        GCi = [ # rows
            [1, 1], [1, 1], [1, 1],
            [1, 2], [1, 2],
            [1, 3], [1, 3],
            [1, 3],
        ]
        polar = 0
        Real = np.array([1., 2., 3.,
                         4., 5.,
                         6., 7.,
                         8.])
        ncols = 4  # ????
        ifo = 2 # rectangular (also 9...)

        name = 'TEST'
        model.add_dmig(name, ifo, tin, ncols,
                       GCj, GCi, Real, Complex=None,
                       tout=tout, polar=polar, comment='')
        get_matrices(model)
        test = model.dmig[name]
        A1, _, _ = test.get_matrix(is_sparse=False, apply_symmetry=True)
        A2, _, _ = test.get_matrix(is_sparse=True, apply_symmetry=True)

        A_expected = np.array([
            [1., 2., 3., 0.],
            [0., 4., 5., 0.],
            [0., 0., 13., 8.],
        ])
        assert np.array_equal(A1, A_expected)
        assert np.array_equal(A2.toarray(), A_expected)

        test.tin = 3   # complex
        test.tout = 3  # complex
        test.Complex = 10 * test.Real
        B1, _, _ = test.get_matrix(is_sparse=False, apply_symmetry=True)
        B2, _, _ = test.get_matrix(is_sparse=True, apply_symmetry=True)
        B_expected = A_expected + 10j * A_expected
        assert np.array_equal(B1, B_expected)
        assert np.array_equal(B2.toarray(), B_expected)

        get_matrices(model)
        #kaax = model.dmig['KAAX'].get_matrix(is_sparse=True)

    def test_dmig_column(self):
        """testing symmetric DMIGs"""
        model = BDF(debug=False, log=None, mode='msc')

        tin = 1 # real, float32
        tout = 1 # real, float32

        # [
        #   1., 2., 3., 0.
        #       4., 5., 0.
        #           13., 8.
        # ]
        GCj = [ # cols
            [1, 1], [1, 1], [1, 1], [1, 1],
        ]
        GCi = [ # rows
            [1, 1],
            [1, 2],
            [1, 3],
            [1, 3],
        ]
        polar = 0
        Real = np.array([1., 2., 3., 4.])
        ncols = 1  # ????
        ifo = 2 # rectangular (also 9...)

        name = 'TEST'
        model.add_dmig(name, ifo, tin, ncols,
                       GCj, GCi, Real, Complex=None,
                       tout=tout, polar=polar, comment='')
        get_matrices(model)
        test = model.dmig[name]
        A1, _, _ = test.get_matrix(is_sparse=False, apply_symmetry=True)
        A2, _, _ = test.get_matrix(is_sparse=True, apply_symmetry=True)

        A_expected = np.array([
            [1.],
            [2.],
            [7.],
        ])
        assert np.array_equal(A1, A_expected)
        assert np.array_equal(A2.toarray(), A_expected)

        test.tin = 3   # complex
        test.tout = 3  # complex
        test.Complex = 10 * test.Real
        B1, _, _ = test.get_matrix(is_sparse=False, apply_symmetry=True)
        B2, _, _ = test.get_matrix(is_sparse=True, apply_symmetry=True)
        B_expected = A_expected + 10j * A_expected
        assert np.array_equal(B1, B_expected)
        assert np.array_equal(B2.toarray(), B_expected)

        get_matrices(model)
        #kaax = model.dmig['KAAX'].get_matrix(is_sparse=True)

    def test_dmig_uaccel(self):
        """tests DMIG,UACCEL"""
        model = BDF(debug=None)

        lines = [
            ['DMIG', 'UACCEL', '0', 9, 1, None, None, None, 4],
            ['DMIG', 'UACCEL', 2, None, None, 2, 3, 386.4],
            ['DMIG', 'UACCEL', 3, None, None, 2, 4, 3.0],
            ['DMIG', 'UACCEL', 4, None, None, 2, 6, 1.0],
        ]
        for line in lines:
            model.add_card(line, 'DMIG')
        fill_dmigs(model)
        str(model.dmig)
        str(model.get_bdf_stats())

        tin = 1
        ncol = 3
        load_sequences = {
            # lseq : [ncx]
            10 : [[1001, 4, 1003.0]],
            20 : [[1, 2, 3.], [11, 12, 13.], ],
            30 : [[1, 2, 3.], [11, 12, 13.], [101, 4, 103.]],
        }
        uaccel = model.add_dmig_uaccel(tin, ncol, load_sequences, comment='uaccel')
        uaccel.raw_fields()
        str(uaccel)
        uaccel.tin = 2
        str(uaccel)
        save_load_deck(model)



class TestDMIGImag(unittest.TestCase):
    def test_dmig_1(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(TEST_PATH, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)
        out = model.dmig['IMAG'].get_matrix(is_sparse=False)

        imag_actual, unused_rows_reversed, unused_cols_reversed = out
        #print "---IMAG_actual---\n", IMAG_actual
        imag_expected_real = [
            [1.0, 0.5, 0.25],
            [0.0, 2.0, 0.75],
            [0.0, 0.0, 3.0],
        ]
        imag_expected_imag = [
            [1.1, 0.51, 0.251],
            [0.0, 2.1, 0.751],
            [0.0, 0.0, 3.1],
        ]
        a_matrix = model.dmig['IMAG']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        imag_expected = np.array(imag_expected_real) + np.array(imag_expected_imag)*1j
        self.assertTrue(np.allclose(imag_expected, imag_actual))
        a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_2(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(TEST_PATH, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)

        out = model.dmig['IMAGS'].get_matrix(is_sparse=False)
        imags_actual, unused_rows_reversed, unused_cols_reversed = out
        #print("---imag_actual---\n", imag_actual)
        imags_expected_real = [
            [1.0, 0.5, 0.25],
            [0.5, 2.0, 0.75],
            [0.25, 0.75, 3.0],
        ]
        imags_expected_imag = [
            [1.1, 0.51, 0.251],
            [0.51, 2.1, 0.751],
            [0.251, 0.751, 3.1],
        ]
        a_matrix = model.dmig['IMAGS']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        imags_expected = np.array(imags_expected_real) + np.array(imags_expected_imag)*1j
        msg = '\n%s_actual\n%s\n\n----' % ('IMAGS', imags_actual)
        msg += '\n%s_expected\n%s\n----' % ('IMAGS', imags_expected)
        msg += '\n%s_delta\n%s\n----' % ('IMAGS', imags_actual-imags_expected)
        self.assertTrue(np.allclose(imags_expected, imags_actual), msg)
        a_matrix.get_matrix()
        save_load_deck(model)

    def test_dmig_3(self):
        """tests the add card method with a complex DMIG"""
        model = BDF(debug=False)
        name = 'DMIG_1'
        ifo = 6
        tin = 3
        ncols = None
        reals = np.array([1.0, 2.0, 3.0])
        #complexs = reals
        GCj = [[1, 1],  # grid, component
               [2, 1],
               [3, 1]]
        GCi = [[1, 1],  # grid, component
               [4, 1],
               [5, 1]]
        dmig = model.add_dmig(name, ifo, tin, ncols, GCj, GCi,
                              Real=reals, Complex=10*reals,
                              comment='dmig')
        model.pop_parse_errors()
        assert dmig.is_real is False, dmig.is_real
        assert dmig.is_complex is True, dmig.is_complex
        assert dmig.is_polar is False, dmig.is_polar
        dmig.get_matrix()

        name = 'DMIK_1'
        unused_nrows = None
        unused_form = None
        dmik = model.add_dmik(name, ifo, tin, ncols, GCj, GCi,
                              Real=reals, Complex=10*reals,
                              # tout=tout, polar=polar,
                              comment='dmik')
        dmik.get_matrix()
        save_load_deck(model)

    def test_dmig_4(self):
        """tests the add card method with a complex polar DMIG"""
        model = BDF(debug=False)
        name = 'DMIG_1'
        ifo = 6
        tin = 3
        ncols = None
        reals = [1.0, 2.0, 3.0]
        #complexs = reals
        GCj = [[1, 1],  # grid, component
               [2, 1],
               [3, 1]]
        GCi = [[1, 1],  # grid, component
               [4, 1],
               [5, 1]]
        dmig = model.add_dmig(name, ifo, tin, ncols, GCj, GCi,
                              Real=reals, Complex=reals,
                              polar=1, comment='dmig')
        model.pop_parse_errors()
        assert dmig.is_real is False, dmig.is_real
        assert dmig.is_complex is True, dmig.is_complex
        assert dmig.is_polar is True, dmig.is_polar
        dmig.get_matrix()

        name = 'DMIK_1'
        unused_nrows = None
        unused_form = None
        dmik = model.add_dmik(name, ifo, tin, ncols, GCj, GCi,
                              Real=reals, Complex=reals,
                              # tout=tout, polar=polar,
                              comment='dmik')
        dmik.get_matrix()
        save_load_deck(model)


class TestDMIAX(unittest.TestCase):
    def test_dmiax(self):
        """tests DMIAX"""
        model = BDF(debug=None)

        #ifo = 1  # square
        tin = 1
        tout = None
        #ncol = None
        #gj = 2
        #cj = 8
        #g1 = 10
        #c1 = 20
        #a1 = 30
        #b1 = 40
        lines = [

            ['DMIAX', 'B2PP', 0, 1, 3],
            ['DMIAX', 'B2PP', 32, None, None, None, None, None, None,
             1027, 3, None, 4.25+6, '2.27+3'],
            #['DMIAX', 'AX', '0', ifo, tin, tout, polar, None, ncol],
            #['DMIAX', 'AX', gj, cj, None, g1, c1, a1, b1],
        ]
        for line in lines:
            model.add_card(line, 'DMIAX')
        fill_dmigs(model)
        str(model.dmiax)
        str(model.get_bdf_stats())
        str(model.dmiax['B2PP'])

        name = 'AXREAL1'
        tin = 1
        tout = 2
        ncols = 3
        GCNj = [(1, 2, 3)]
        GCNi = [[(10, 6, 30)]]
        Real = [0.1]
        matrix_form = 1
        dmiax_real = model.add_dmiax(name, matrix_form, tin, ncols,
                                     GCNj, GCNi, Real, Complex=None,
                                     tout=tout, comment='dmiax')
        assert dmiax_real.is_real is True, dmiax_real.is_real
        assert dmiax_real.is_complex is False, dmiax_real.is_complex
        assert dmiax_real.is_polar is False, dmiax_real.is_polar

        name = 'AXREAL1'
        tin = 2
        tout = 1
        dmiax_real = model.add_dmiax(name, matrix_form, tin, ncols,
                                     GCNj, GCNi, Real, Complex=None,
                                     tout=tout, comment='dmiax')
        assert dmiax_real.is_real is True, dmiax_real.is_real
        assert dmiax_real.is_complex is False, dmiax_real.is_complex
        assert dmiax_real.is_polar is False, dmiax_real.is_polar

        name = 'AXIMAG1'
        Complex = [0.2]
        tin = 3
        tout = 4
        dmiax_mag = model.add_dmiax(name, matrix_form, tin, ncols,
                                    GCNj, GCNi, Real, Complex=Complex,
                                    tout=tout, comment='dmiax')
        assert dmiax_mag.is_real is False, dmiax_mag.is_real
        assert dmiax_mag.is_complex is True, dmiax_mag.is_complex
        assert dmiax_mag.is_polar is False, dmiax_mag.is_polar

        name = 'AXIMAG2'
        Complex = [0.2]
        tout = 4
        tin = 3
        dmiax_imag = model.add_dmiax(name, matrix_form, tin, ncols,
                                     GCNj, GCNi, Real, Complex=Complex,
                                     tout=tout, comment='dmiax')
        assert dmiax_imag.is_real is False, dmiax_imag.is_real
        assert dmiax_imag.is_complex is True, dmiax_imag.is_complex
        assert dmiax_imag.is_polar is False, dmiax_imag.is_polar

        dmiax_real.raw_fields()
        dmiax_imag.raw_fields()
        str(dmiax_real)
        str(dmiax_imag)
        save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
