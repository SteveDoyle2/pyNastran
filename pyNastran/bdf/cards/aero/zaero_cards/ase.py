# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * TRIM

All cards are BaseCard objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string, string_or_blank,
    double_or_blank, integer_or_string, integer_or_double,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class ASE(BaseCard):
    type = 'ASE'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, ase_id: int, asecont_id: int, flutter_id: int,
                 mldstat_id: int, minstat_id: int, cmargin_id: int,
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.ase_id = ase_id
        self.asecont_id = asecont_id
        self.flutter_id = flutter_id
        self.mldstat_id = mldstat_id
        self.minstat_id = minstat_id
        self.cmargin_id = cmargin_id
        self.asecont_ref = None
        self.flutter_ref = None
        self.mldstat_ref = None
        self.minstat_ref = None
        self.asecont_ref = None
        self.cmargin_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # SE  SID CONID FLTID STATES RAAID MARID
        # ASE 5   20    30    1      40    50
        ase_id = integer(card, 1, 'ase_id')
        asecont_id = integer_or_blank(card, 2, 'conid, asecont_id', default=0)
        flutter_id = integer(card, 3, 'fltid, flutter_id')
        mldstat_id = integer_or_blank(card, 4, 'states, mldstat_id', default=0)
        minstat_id = integer_or_blank(card, 5, 'raaid, minstat_id', default=0)
        cmargin_id = integer_or_blank(card, 6, 'marid, cmargin_id', default=0)

        assert 4 <= len(card) <= 7, f'len(ASE card) = {len(card):d}\ncard={card}'
        return ASE(ase_id, asecont_id, flutter_id, mldstat_id, minstat_id, cmargin_id, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        zaero = model.zaero
        if self.asecont_id:
            self.asecont_ref = zaero.asecont[self.asecont_id]
        self.flutter_ref = model.flutters[self.flutter_id]
        if self.mldstat_id:
            self.mldstat_ref = zaero.mldstat[self.mldstat_id]
        if self.minstat_id:
            self.minstat_ref = zaero.minstat[self.minstat_id]
        if self.cmargin_id and 0:
            self.cmargin_ref = zaero.cmargin[self.cmargin_id]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.asecont_ref = None
        self.flutter_ref = None
        self.mldstat_ref = None
        self.minstat_ref = None
        self.cmargin_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'ASE', self.ase_id, self.asecont_id, self.flutter_id,
            self.mldstat_id, self.minstat_id, self.cmargin_id]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class ASECONT(BaseCard):
    type = 'ASECONT'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, asecont_id: int, surf_id: int, sens_id: int,
                 tf_id: int, gain_id: int, conct_id: int,
                 extinp_set_id: int, extout_set_id: int, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.asecont_id = asecont_id
        self.surf_id = surf_id
        self.sens_id = sens_id
        self.tf_id = tf_id
        self.gain_id = gain_id
        self.conct_id = conct_id
        self.extinp_set_id = extinp_set_id
        self.extout_set_id = extout_set_id

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a ASECONT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # ASECONT SID SURFID SENSID TFID GAINID CONCTID EINPID EOUTID
        # ASECONT 10  20     30     40   50     60

        asecont_id = integer(card, 1, 'asecont_id')
        surf_id = integer(card, 2, 'surf_id')
        sens_id = integer(card, 3, 'sens_id')
        tf_id = integer_or_blank(card, 4, 'tf_id', default=0)
        gain_id = integer_or_blank(card, 5, 'gain_id', default=0)
        conct_id = integer_or_blank(card, 6, 'conct_id', default=0)
        extinp_set_id = integer_or_blank(card, 7, 'extinp_set_id', default=0)
        extout_set_id = integer_or_blank(card, 8, 'extout_set_id', default=0)

        assert len(card) <= 9, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASECONT(asecont_id, surf_id, sens_id, tf_id, gain_id,
                          conct_id, extinp_set_id, extout_set_id, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by ASECONT={self.asecont_id}'
        # ASE, MLOADS, ELOADS, GLOADS, DFS, or NLFLTR
        # zaero = model.zaero
        # CNCTSET
        # self.conct_ref = model.conct[self.conct_id]
        if self.extinp_set_id:
            self.extinp_set_ref = model.Set(self.extinp_set_id, msg)
            for idi in self.extinp_set_ref.ids():
                asdf
        if self.extout_set_id:
            self.extout_set_ref = model.Set(self.extout_set_id, msg)
            for idi in self.extout_set_ref.ids():
                asdf


    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self) -> list:
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'ASECONT', self.asecont_id, self.surf_id, self.sens_id, self.tf_id, self.gain_id,
            self.conct_id, self.extinp_set_id, self.extout_set_id]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class ASEGAIN(BaseCard):
    type = 'ASEGAIN'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, asegain_id: int,
                 otf_id: int, c_out: int,
                 itf_id: int, c_in: int,
                 gain: float, gain_type: str, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.asegain_id = asegain_id
        self.otf_id = otf_id
        self.c_out = c_out
        self.itf_id = itf_id
        self.c_in = c_in
        self.gain = gain
        self.gain_type = gain_type
        self.input_ref = None
        self.output_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a ASECONT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # ASEGAIN ID OTFID CO ITFID CI GAIN TYPE
        # ASEGAIN 10 1     1  2     2  0.6  H
        # ['ASEGAIN', '301', '100', '1', '209', '1', '425', 'DEN']
        asegain_id = integer(card, 1, 'asegain_id')
        otf_id = integer(card, 2, 'otfid')
        c_out = integer(card, 3, 'CO')
        itf_id = integer(card, 4, 'itfid')
        c_in = integer(card, 5, 'CI')
        gain = integer_or_double(card, 6, 'gain')
        gain_type = string_or_blank(card, 7, 'gain_type', default='Q')

        assert len(card) < 9, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASEGAIN(asegain_id, otf_id, c_out, itf_id, c_in,
                          gain, gain_type, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by ASEGAIN={self.asegain_id}'
        # ASE, MLOADS, ELOADS, GLOADS, DFS, or NLFLTR
        zaero = model.zaero
        # CNCTSET
        # self.conct_ref = model.conct[self.conct_id]

        # CJUNCT, MIMOSS, SISOTF or ACTU
        self.input_ref = zona_cjunct_mimoss_sisotf_actu(
            'input', self.itf_id, zaero, f'ASEGAIN={self.asegain_id}')

        # CJUNCT, MIMOSS, SISOTF, ASESNSR, or ASESNS1
        self.output_ref = zona_cjunct_mimoss_sisotf_asesnsr_asesns1(
            'Output', self.otf_id, zaero, f'ASEGAIN={self.asegain_id}')

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        # ASEGAIN ID OTFID CO ITFID CI GAIN TYPE
        # ASEGAIN 10 1     1  2     2  0.6  H
        list_fields = [
            'ASEGAIN', self.asegain_id,
            self.otf_id, self.c_out,
            self.itf_id, self.c_in,
            self.gain, self.gain_type]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class ASESNSR(BaseCard):
    type = 'ASESNSR'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, asesnsr_id: int, sensor_type: int,
                 sgid: int, component: int,
                 factor: float, sum_method: str, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.asesnsr_id = asesnsr_id
        self.sensor_type = sensor_type
        self.sgid = sgid
        self.component = component
        self.factor = factor
        self.sum_method = sum_method

    @property
    def name(self):
        # {1-6}
        mapper = {
            # sensor_type
            0: {
                1: 'x', 2: 'y', 3: 'z',
                4: 'ϕ', 5: 'θ', 6: 'ψ',
            },
            1: {
                1: 'dx/dt', 2: 'dy/dt', 3: 'dz/dt',
                4: 'dϕ/dt', 5: 'dθ/dt', 6: 'dψ/dt',
            },
            2: {
                1: 'd²x/dt²', 2: 'd²y/dt²', 3: 'd²z/dt²',
                4: 'd²ϕ/dt²', 5: 'd²θ/dt²', 6: 'd²ψ/dt²',
            }
        }
        return mapper[self.sensor_type][self.component]

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a ASECONT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # ASESNSR ID TYPE SGID SC FACTOR SOF
        # ASESNSR 10 2    7    6  0.0176 YES

        asesnsr_id = integer(card, 1, 'asesnsr_id')
        sensor_type = integer(card, 2, 'sensor_type')
        sgid = integer(card, 3, 'sens_id')
        component = integer(card, 4, 'component')
        factor = double_or_blank(card, 5, 'factor', default=1.0)
        sum_method = string_or_blank(card, 6, 'conct_id', default='NO')

        assert len(card) <= 7, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASESNSR(asesnsr_id, sensor_type, sgid, component, factor,
                          sum_method, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        # SGID
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'ASESNSR', self.asesnsr_id, self.sensor_type,
            self.sgid, self.component, self.factor, self.sum_method]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class ASESNS1(BaseCard):
    type = 'ASESNS1'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, asesns1_id: int,
                 label: str, ikey: int | str,
                 factor: float, sum_method: str='NO', comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.asesns1_id = asesns1_id
        self.label = label
        self.ikey = ikey
        self.factor = factor
        self.sum_method = sum_method

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a ASECONT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # ASESNS1 ID LABEL   IKEY SOF FACTOR
        # ASESNS1 10 GRIDGT3 7   NO   0.0176
        asesns1_id = integer(card, 1, 'asesns1_id')
        label = string(card, 2, 'label')
        ikey = integer_or_string(card, 3, 'ikey')
        sum_method = string_or_blank(card, 6, 'sum_method', default='NO')
        factor = double_or_blank(card, 5, 'factor', default=1.0)

        assert len(card) < 7, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASESNS1(asesns1_id, label, ikey, factor,
                          sum_method=sum_method, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        # SGID
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'ASESNS1', self.asesns1_id, self.label,
            self.ikey, self.sum_method, self.factor]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CJUNCT(BaseCard):
    type = 'CJUNCT'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, cjunct_id: int, nu: int, ny: int,
                 values: list[float], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        values = np.asarray(values)
        self.cjunct_id = cjunct_id
        self.nu = nu
        self.ny = ny
        self.values = values.reshape(self.nu, self.ny) # (ninput, noutput)

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # CJUNCT ID  NU NY D11 D21 -etc- D(NY,1) D12 CONT
        #        D22 -etc- D(NY,NU)
        # CJUNCT 90 3 1 1.0 -1.0 0.5
        cjunct_id = integer(card, 1, 'cjunct_id')
        nu = integer(card, 2, 'nu')
        ny = integer(card, 3, 'ny')
        values = []
        for ifield in range(4, len(card)):
            di = double(card, ifield, 'di')
            values.append(di)

        nexpected = nu * ny
        assert len(values) == nexpected, f'len(values) ={len(values):d}; expected={nu}*{ny}={nexpected}'
        assert len(card) >= 5, f'len(CJUNCT card) = {len(card):d}\ncard={card}'
        return CJUNCT(cjunct_id, nu, ny, values, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'CJUNCT', self.cjunct_id, self.nu,
            self.ny] + self.values.ravel().tolist()
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CONCT(BaseCard):
    type = 'CONCT'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, conct_id: int,
                 output_tf_id: int, output_component: int,
                 input_tf_id: int, input_component: int,
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.conct_id = conct_id
        self.output_tf_id = output_tf_id
        self.output_component = output_component
        self.input_tf_id = input_tf_id
        self.input_component = input_component

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # CONCT ID OTFID CO ITFID CI
        # CONCT 10 1     1  2     2
        conct_id = integer(card, 1, 'conct_id')

        # CJUNCT, MIMOSS, SISOTF, ASESNSR, or ASESNS1
        output_tf_id = integer(card, 2, 'OTFID, output_tf_id')
        output_component = integer(card, 3, 'C0, output_component')
        # CJUNCT, MIMOSS, SISOTF or ACTU
        input_tf_id = integer(card, 4, 'ITFID, input_tf_id')
        input_component = integer(card, 5, 'CI, input_component')
        assert len(card) == 6, f'len(CJUNCT card) = {len(card):d}\ncard={card}'
        return CONCT(conct_id, output_tf_id, output_component, input_tf_id, input_component, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        zaero = model.zaero
        log = model.log

        # CJUNCT, MIMOSS, SISOTF or ACTU
        self.input_ref = zona_cjunct_mimoss_sisotf_actu(
            'input', self.input_tf_id,
            zaero, f'CONCT={self.conct_id}')

        #-------------
        idi = self.output_tf_id
        # CJUNCT, MIMOSS, SISOTF, ASESNSR, or ASESNS1
        self.output_ref = zona_cjunct_mimoss_sisotf_asesnsr_asesns1(
            'output', self.output_tf_id,
            zaero, f'CONCT={self.conct_id}')

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'CONCT', self.conct_id,
             self.output_tf_id, self.output_component,
             self.input_tf_id, self.input_component]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TFSET(BaseCard):
    type = 'TFSET'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, tfset_id: int, ids: list[int], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.tfset_id = tfset_id
        self.ids = ids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # TFSET SID TF1 TF2 TF3 TF4
        # TFSET 10  1   2   3
        tfset_id = integer(card, 1, 'cjunct_id')

        # CJUNCT, MIMOSS, SISOTF
        ids = []
        for ifield in range(2, len(card)):
            idi = integer(card, ifield, 'OTFID, id')
            ids.append(idi)
        assert len(card) >= 3, f'len(TFSET card) = {len(card):d}\ncard={card}'
        return TFSET(tfset_id, ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        ids_ref = []
        log = model.log
        zaero = model.zaero
        for idi in self.ids:
            if idi in zaero.cjunct:
                id_ref = zaero.cjunct[idi]
            elif idi in zaero.mimoss:
                id_ref = zaero.mimoss[idi]
            elif idi in zaero.sisotf:
                id_ref = zaero.sisotf[idi]
            else:
                cjunct = list(zaero.cjunct)
                mimoss = list(zaero.mimoss)
                sisotf = list(zaero.sisotf)
                msg = (
                    f'TFSET={self.tfset_id}: id={idi} is not [CJUNCT, MIMOSS, SISOTF]\n'
                    f' - cjunct = {cjunct}\n'
                    f' - mimoss = {mimoss}\n'
                    f' - sisotf = {sisotf}\n'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TFSET', self.tfset_id] + self.ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SENSET(BaseCard):
    type = 'SENSET'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, senset_id: int, ids: list[int], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.senset_id = senset_id
        self.ids = ids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SENSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # TFSET SID TF1 TF2 TF3 TF4
        # TFSET 10  1   2   3
        senset_id = integer(card, 1, 'senset_id')

        # CJUNCT, MIMOSS, SISOTF
        ids = []
        for ifield in range(2, len(card)):
            idi = integer(card, ifield, 'OTFID, id')
            ids.append(idi)
        assert len(card) >= 3, f'len(SENSET card) = {len(card):d}\ncard={card}'
        return SENSET(senset_id, ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        ids_ref = []
        log = model.log
        zaero = model.zaero
        for idi in self.ids:
            # ASESNSR or ASESNS1
            if idi in zaero.asesnsr:
                id_ref = zaero.asesnsr[idi]
            elif idi in zaero.asesns1:
                id_ref = zaero.asesns1[idi]
            else:
                asesnsr = list(zaero.asesnsr)
                asesns1 = list(zaero.asesns1)
                # sisotf = list(zaero.sisotf)
                msg = (
                    f'SENSET={self.senset_id}: id={idi} is not [ASESNSR, ASESNS1]\n'
                    f' - asesnsr = {asesnsr}\n'
                    f' - asesns1 = {asesns1}'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SENSET', self.senset_id] + self.ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class GAINSET(BaseCard):
    type = 'GAINSET'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, gainset_id: int, ids: list[int], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.gainset_id = gainset_id
        self.ids = ids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a GAINSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # GAINSET SID TF1 TF2 TF3 TF4
        # GAINSET 10  1   2   3
        gainset_id = integer(card, 1, 'gainset_id')

        # CJUNCT, MIMOSS, SISOTF
        ids = []
        for ifield in range(2, len(card)):
            idi = integer(card, ifield, 'OTFID, id')
            ids.append(idi)
        assert len(card) >= 3, f'len(GAINSET card) = {len(card):d}\ncard={card}'
        return GAINSET(gainset_id, ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        ids_ref = []
        log = model.log
        zaero = model.zaero
        for idi in self.ids:
            if idi in zaero.asegain:
                id_ref = zaero.asegain[idi]
            else:
                gainset = list(zaero.gainset)
                asegain = list(zaero.asegain)
                gainset.sort()
                asegain.sort()
                msg = (
                    f'GAINSET={self.gainset_id}: id={idi} is not [ASEGAIN]\n'
                    f' - gainset = {gainset}\n'
                    f' - asegain = {asegain}\n'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['GAINSET', self.gainset_id] + self.ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MIMOSS(BaseCard):
    type = 'MIMOSS'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, mimoss_id: int,
                 ntf: int, nu: int, ny: int,
                 dmi_label: str, print_flag: int,
                 values: list[float], labels: list[str],
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mimoss_id = mimoss_id
        self.ntf = ntf
        self.nu = nu
        self.ny = ny
        self.dmi_label = dmi_label
        self.print_flag = print_flag
        self.values = values
        self.labels = labels
        assert len(self.values) == len(self.labels)

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MIMOSS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MIMOSS ID  NTF    NU  NY LMIMO TYPE PRNT
        #        X1  L1     X2  L2 X3    L3   -etc-
        # MIMOSS 50  6      7   4        DEN  1
        #        0.0 MATRC1 1.0 MATRC2
        mimoss_id = integer(card, 1, 'mimoss_id')
        ntf = integer(card, 2, 'ntf')
        nu = integer(card, 3, 'nu')
        ny = integer(card, 4, 'ny')
        dmi_label = string_or_blank(card, 5, 'LMIMO, dmi_label', default='')
        mimoss_type = string_or_blank(card, 6, 'mimoss_type', default='Q')
        assert mimoss_type in ['VEL', 'DEN', 'Q', 'H'], mimoss_type
        # TODO: not sure on default
        print_flag = integer_or_blank(card, 7, 'print_flag', default=0)

        values = []
        labels = []
        if dmi_label: # blank
            # type, print, labels, values not required
            assert len(card) >= 5, f'len(MIMOSS card) = {len(card):d}\ncard={card}'
        else:
            # print(f'dmi_label = {dmi_label}')
            dfields = len(card) - 8
            assert dfields // 2 == 0, dfields
            assert dfields > 0, dfields
            for ifield in range(8, len(card), 2):
                di = double(card, ifield, 'di')
                labeli = string(card, ifield+1, 'labeli')
                values.append(di)
                labels.append(labeli)
            assert len(values) > 0, values
        # nexpected = nu * ny
        # assert len(values) == nexpected, f'len(values) ={len(values):d}; expected={nu}*{ny}={nexpected}'
        return MIMOSS(mimoss_id, ntf, nu, ny,
                      dmi_label, print_flag,
                      values, labels, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'CJUNCT', self.mimoss_id, self.ntf, self.nu, self.ny,
            self.dmi_label, self.print_flag]
        assert len(self.values) == len(self.labels)
        for value, label in zip(self.values, self.labels):
            list_fields.append(value)
            list_fields.append(label)
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SURFSET(BaseCard):
    type = 'SURFSET'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, surfset_id: int, ids: list[int], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.surfset_id = surfset_id
        self.ids = ids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SURFSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        surfset_id = integer(card, 1, 'surfset_id')

        ids = []
        j = 1
        for ifield in range(2, len(card)):
            idi = string(card, ifield, 'SURF, id{j}')
            ids.append(idi)
            j += 1
        assert len(card) >= 3, f'len(SURFSET card) = {len(card):d}\ncard={card}'
        return SURFSET(surfset_id, ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        ids_ref = []
        log = model.log
        zaero = model.zaero
        for idi in self.ids:
            # AESURFZ, AESLINK, PZTMODE or JETFRC
            if idi in model.aesurf:
                id_ref = model.aesurf[idi]
            # elif idi in zaero.aesurfz:
            #     id_ref = zaero.aesurfz[idi]
            elif idi in zaero.aeslink:
                id_ref = zaero.aeslink[idi]
            else:
                aesurf = list(model.aesurf)
                aeslink = list(zaero.aeslink)
                # aesurfz = list(zaero.aesurfz)
                asesnsr = list(zaero.asesnsr)
                # asesns1 = list(zaero.asesns1)
                msg = (
                    f'SURFSET={self.surfset_id}: id={idi} is not [AESURFZ, AESLINK, PZTMODE, JETFRC]\n'
                    f' - aesurf  = {aesurf}\n'
                    f' - aeslink = {aeslink}\n'
                    # f' - aesurfz = {aesurfz}\n'
                    f' - asesnsr = {asesnsr}\n'
                    # f' - asesns1 = {asesns1}'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SURFSET', self.surfset_id] + self.ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CNCTSET(BaseCard):
    type = 'CNCTSET'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, cnctset_id: int, ids: list[int], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.cnctset_id = cnctset_id
        self.ids = ids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a CNCTSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        cnctset_id = integer(card, 1, 'cnctset_id')
        ids = []
        for ifield in range(2, len(card)):
            idi = integer(card, ifield, 'CONCT, id')
            ids.append(idi)
        assert len(card) >= 3, f'len(CNCTSET card) = {len(card):d}\ncard={card}'
        return CNCTSET(cnctset_id, ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        ids_ref = []
        log = model.log
        zaero = model.zaero
        for idi in self.ids:
            # AESURFZ, AESLINK, PZTMODE or JETFRC
            if idi in zaero.cnctset:
                id_ref = zaero.cnctset[idi]
            elif idi in zaero.conct:
                id_ref = zaero.conct[idi]
            else:
                conct = list(zaero.conct)
                cnctset = list(zaero.cnctset)
                msg = (
                    f'CNCTSET={self.cnctset_id}: id={idi} is not [CONCT]\n'
                    f' - cnctset = {cnctset}\n'
                    f' - conct = {conct}\n'
                    # f' - asesns1 = {asesns1}'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['CNCTSET', self.cnctset_id] + self.ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SISOTF(BaseCard):
    type = 'SISOTF'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, sisotf_id: int,
                 nnumerator: int, ndenominator: int,
                 b: list[float], a: list[float],
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.sisotf_id = sisotf_id
        self.nnumerator = nnumerator
        self.ndenominator = ndenominator
        self.b = b # numerator
        self.a = a # denominator
        assert 0 <= self.nnumerator <= self.ndenominator, f'nnumerator={nnumerator}; ndenominator={ndenominator}'
        assert len(b) == nnumerator+1, f'nb={len(b)}; nnumerator={nnumerator}; b={b} a={a}'
        assert len(a) == ndenominator, f'na={len(a)}; ndenominator={ndenominator}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MIMOSS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # SISOTF ID NDEN  NNUM      A0 A1 A2 A3 A4
        #        A5 -etc- A(NDEN-1) B0 B1 B2 B3 B4
        #        B5 -etc- B(NNUM)
        sisotf_id = integer(card, 1, 'sisotf_id')
        ndenominator = integer(card, 2, 'NDEN')
        nnumerator = integer(card, 3, 'NNUM')
        j = 1
        a = []
        ifield0 = 4
        # nvalues = na + (nb + 1)
        # print(f'na={ndenominator} nb={nnumerator}')
        for ifield in range(ifield0, ifield0+ndenominator):
            ai = double(card, ifield, f'a{j}')
            j += 1
            a.append(ai)

        j = 1
        b = []
        ifield0 += ndenominator
        for ifield in range(ifield0, ifield0+nnumerator+1):
            bi = double(card, ifield, f'b{j}')
            j += 1
            b.append(bi)
        return SISOTF(sisotf_id, nnumerator, ndenominator,
                      b, a, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'SISOTF', self.sisotf_id, self.ndenominator, self.nnumerator] + \
            self.a + self.b
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEROLAG(BaseCard):
    type = 'AEROLAG'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, aerolag_id: int, nlag: int, lag_values: list[float],
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.aerolag_id = aerolag_id
        self.nlag = nlag
        self.lag_values = lag_values

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # AEROLAG SID NLAG R1    R2   R3    R4 R5 R6
        #         R7  R8   -etc-
        # AEROLAG 10  4  -0.2    -0.5 -1.0 -2.0
        aerolag_id = integer(card, 1, 'aerolag_id')
        nlag = integer(card, 1, 'nlag')

        lag_values = []
        for ifield in range(3, len(card)):
            idi = double(card, ifield, 'OTFID, id')
            lag_values.append(idi)
        # assert len(lag_values) == nlag, (lag_values, nlag)
        assert len(card) >= 3, f'len(AEROLAG card) = {len(card):d}\ncard={card}'
        return AEROLAG(aerolag_id, nlag, lag_values, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['AEROLAG', self.aerolag_id, self.nlag] + self.lag_values
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


def zona_cjunct_mimoss_sisotf_asesnsr_asesns1(input_output: str, idi: int, zaero, msg: str):
    if idi in zaero.sisotf:
        input_ref = zaero.sisotf[idi]
        assert idi not in zaero.cjunct, f'idi={idi} in sisotf (already in sisotf)'
        assert idi not in zaero.mimoss, f'idi={idi} in sisotf (already in mimoss)'
        assert idi not in zaero.asesnsr, f'idi={idi} in sisotf (already in asesnsr)'
        assert idi not in zaero.asesns1, f'idi={idi} in sisotf (already in asesns1)'
    elif idi in zaero.cjunct:
        input_ref = zaero.cjunct[idi]
        assert idi not in zaero.mimoss, f'idi={idi} in cjunct (already in mimoss)'
        assert idi not in zaero.asesnsr, f'idi={idi} in cjunct (already in asesnsr)'
        assert idi not in zaero.asesns1, f'idi={idi} in cjunct (already in asesns1)'
    elif idi in zaero.mimoss:
        input_ref = zaero.mimoss[idi]
        assert idi not in zaero.asesnsr, f'idi={idi} in mimoss (already in asesnsr)'
        assert idi not in zaero.asesns1, f'idi={idi} in mimoss (already in asesns1)'
    elif idi in zaero.asesnsr:
        input_ref = zaero.asesnsr[idi]
        assert idi not in zaero.asesns1, f'idi={idi} in asesns1 (already in asesnsr)'
    elif idi in zaero.asesns1:
        input_ref = zaero.asesns1[idi]
    else:
        cjunct = list(zaero.cjunct)
        mimoss = list(zaero.mimoss)
        sisotf = list(zaero.sisotf)
        asesnsr = list(zaero.asesnsr)
        asesns1 = list(zaero.asesns1)
        cjunct.sort()
        mimoss.sort()
        sisotf.sort()
        asesnsr.sort()
        asesns1.sort()
        msg = (
            f'{msg}: {input_output}={idi} is not [CJUNCT, MIMOSS, SISOTF, ASESNSR, ASESNS1]\n'
            f' - cjunct  = {cjunct}\n'
            f' - mimoss  = {mimoss}\n'
            f' - sisotf  = {sisotf}\n'
            f' - asesnsr = {asesnsr}\n'
            f' - asesns1 = {asesns1}\n')
        zaero.model.log.warning(msg)
        input_ref = None
    return input_ref


def zona_cjunct_mimoss_sisotf_actu(input_output: str, idi: int, zaero, msg: str,
                                   allow_actu: bool=True):
    if idi in zaero.sisotf:
        input_ref = zaero.sisotf[idi]
        assert idi not in zaero.cjunct, f'idi={idi} in sisotf (already in cjunct)'
        assert idi not in zaero.mimoss, f'idi={idi} in sisotf (already in mimoss)'
        assert idi not in zaero.actu, f'idi={idi} in sisotf (already in actu)'
    elif idi in zaero.cjunct:
        input_ref = zaero.cjunct[idi]
        assert idi not in zaero.mimoss, f'idi={idi} in cjunct (already in mimoss)'
        assert idi not in zaero.actu, f'idi={idi} in cjunct (already in actu)'
    elif idi in zaero.mimoss:
        input_ref = zaero.mimoss[idi]
        assert idi not in zaero.actu, f'idi={idi} in mimoss (already in actu)'
    elif idi in zaero.actu and allow_actu:
        input_ref = zaero.actu[idi]
    else:
        cjunct = list(zaero.cjunct)
        mimoss = list(zaero.mimoss)
        sisotf = list(zaero.sisotf)
        actu = list(zaero.actu)
        cjunct.sort()
        mimoss.sort()
        sisotf.sort()
        actu.sort()
        msg = (
            f'{msg}: {input_output}={idi} is not [CJUNCT, MIMOSS, SISOTF, ACTU]\n'
            f' - cjunct = {cjunct}\n'
            f' - mimoss = {mimoss}\n'
            f' - sisotf = {sisotf}\n'
            f' - actu   = {actu}\n')
        zaero.model.log.warning(msg)
        input_ref = None
    return input_ref
