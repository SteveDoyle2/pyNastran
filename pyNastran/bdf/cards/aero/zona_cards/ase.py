# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * TRIM

All cards are BaseCard objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.cards.aero import zona
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string, string_or_blank,
    double_or_blank, integer_or_string,
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
                 comment: str = ''):
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

        assert 4 <= len(card) < 8, f'len(ASE card) = {len(card):d}\ncard={card}'
        return ASE(ase_id, asecont_id, flutter_id, mldstat_id, minstat_id, cmargin_id, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        zona = model.zona
        if self.asecont_id:
            self.asecont_ref = zona.asecont[self.asecont_id]
        self.flutter_ref = model.flutters[self.flutter_id]
        if self.mldstat_id:
            self.mldstat_ref = zona.mldstat[self.mldstat_id]
        if self.minstat_id:
            self.minstat_ref = zona.minstat[self.minstat_id]
        if self.cmargin_id:
            self.asecont_ref = zona.cmargin[self.cmargin_id]

    def safe_cross_reference(self, model: BDF):
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

        assert len(card) < 9, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASECONT(asecont_id, surf_id, sens_id, tf_id, gain_id,
                          conct_id, extinp_set_id, extout_set_id, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by ASECONT={self.asecont_id}'
        # ASE, MLOADS, ELOADS, GLOADS, DFS, or NLFLTR
        zona = model.zona
        # CNCTSET
        # self.conct_ref = model.conct[self.conct_id]
        if self.extinp_set_id:
            self.extinp_set_ref = model.Set(self.extinp_set_id, msg)
            for id in self.extinp_set_ref.ids():
                asdf
        if self.extout_set_id:
            self.extout_set_ref = model.Set(self.extout_set_id, msg)
            for id in self.extout_set_ref.ids():
                asdf


    def safe_cross_reference(self, model: BDF):
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
        asegain_id = integer(card, 1, 'asegain_id')
        otf_id = integer(card, 2, 'otfid')
        c_out = integer(card, 3, 'CO')
        itf_id = integer(card, 4, 'itfid')
        c_in = integer(card, 5, 'CI')
        gain = double(card, 6, 'gain')
        gain_type = string_or_blank(card, 7, 'gain_type', default='Q')

        assert len(card) < 9, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASEGAIN(asegain_id, otf_id, c_out, itf_id, c_in,
                          gain, gain_type, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by ASECONT={self.asecont_id}'
        # ASE, MLOADS, ELOADS, GLOADS, DFS, or NLFLTR
        zona = model.zona
        # CNCTSET
        # self.conct_ref = model.conct[self.conct_id]
        if self.extinp_set_id:
            self.extinp_set_ref = model.Set(self.extinp_set_id, msg)
            for id in self.extinp_set_ref.ids():
                asdf
        if self.extout_set_id:
            self.extout_set_ref = model.Set(self.extout_set_id, msg)
            for id in self.extout_set_ref.ids():
                asdf


    def safe_cross_reference(self, model: BDF):
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
            'ASECONT', self.asegain_id,
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

        assert len(card) < 7, f'len(ASECONT card) = {len(card):d}\ncard={card}'
        asecont = ASESNSR(asesnsr_id, sensor_type, sgid, component, factor,
                          sum_method, comment=comment)
        return asecont

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        # SGID
        pass

    def safe_cross_reference(self, model: BDF):
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

    def safe_cross_reference(self, model: BDF):
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

        self.cjunct_id = cjunct_id
        self.nu = nu
        self.ny = ny
        self.values = values

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

    def safe_cross_reference(self, model: BDF):
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
            self.ny] + self.values
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

    def __init__(self, cjunct_id, output_tf_id, output_component,
                 input_tf_id, input_component, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.cjunct_id = cjunct_id
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
        cjunct_id = integer(card, 1, 'cjunct_id')

        # CJUNCT, MIMOSS, SISOTF, ASESNSR, or ASESNS1
        output_tf_id = integer(card, 2, 'OTFID, output_tf_id')
        output_component = integer(card, 3, 'C0, output_component')
        # CJUNCT, MIMOSS, SISOTF or ACTU
        input_tf_id = integer(card, 4, 'ITFID, input_tf_id')
        input_component = integer(card, 5, 'CI, input_component')
        assert len(card) == 6, f'len(CJUNCT card) = {len(card):d}\ncard={card}'
        return CONCT(cjunct_id, output_tf_id, output_component, input_tf_id, input_component, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF):
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
            'CONCT', self.cjunct_id,
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
        zona = model.zona
        for idi in self.ids:
            if idi in zona.cjunct:
                id_ref = zona.cjunct[idi]
            elif idi in zona.mimoss:
                id_ref = zona.mimoss[idi]
            elif idi in zona.sisotf:
                id_ref = zona.sisotf[idi]
            else:
                cjunct = list(zona.cjunct)
                mimoss = list(zona.mimoss)
                sisotf = list(zona.sisotf)
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

    def safe_cross_reference(self, model: BDF):
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
        zona = model.zona
        for idi in self.ids:
            # ASESNSR or ASESNS1
            if idi in zona.asesnsr:
                id_ref = zona.asesnsr[idi]
            elif idi in zona.asesns1:
                id_ref = zona.asesns1[idi]
            else:
                asesnsr = list(zona.asesnsr)
                asesns1 = list(zona.asesns1)
                # sisotf = list(zona.sisotf)
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

    def safe_cross_reference(self, model: BDF):
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
        zona = model.zona
        for idi in self.ids:
            if idi in zona.asegain:
                id_ref = zona.asegain[idi]
            # elif idi in zona.mimoss:
            #     id_ref = zona.mimoss[idi]
            # elif idi in zona.sisotf:
            #     id_ref = zona.sisotf[idi]
            # else:
                asegain = list(zona.asegain)
                # mimoss = list(zona.mimoss)
                # sisotf = list(zona.sisotf)
                msg = (
                    f'GAINSET={self.gainset_id}: id={idi} is not [ASEGAIN]\n'
                    f' - asegain = {asegain}\n'
                    # f' - mimoss = {mimoss}\n'
                    # f' - sisotf = {sisotf}\n'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF):
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
            print(f'dmi_label = {dmi_label}')
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

    def safe_cross_reference(self, model: BDF):
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
        zona = model.zona
        for idi in self.ids:
            # AESURFZ, AESLINK, PZTMODE or JETFRC
            if idi in model.aesurf:
                id_ref = model.aesurf[idi]
            else:
                aesurf = list(model.aesurf)
                asesnsr = list(zona.asesnsr)
                # asesns1 = list(zona.asesns1)
                msg = (
                    f'SURFSET={self.surfset_id}: id={idi} is not [AESURFZ, AESLINK, PZTMODE, JETFRC]\n'
                    f' - aesurf  = {aesurf}\n'
                    f' - asesnsr = {asesnsr}\n'
                    # f' - asesns1 = {asesns1}'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF):
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
        zona = model.zona
        for idi in self.ids:
            # AESURFZ, AESLINK, PZTMODE or JETFRC
            if idi in zona.cnctset:
                id_ref = zona.cnctset[idi]
            elif idi in zona.conct:
                id_ref = zona.conct[idi]
            else:
                conct = list(zona.conct)
                cnctset = list(zona.cnctset)
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

    def safe_cross_reference(self, model: BDF):
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
