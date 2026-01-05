from __future__ import annotations
from itertools import count
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, double_or_string,
    double_string_or_blank,
)
from pyNastran.bdf.cards.aero.aero import AELINK
from pyNastran.bdf.cards.aero.static_loads import TRIM
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.bdf.cards.aero.zaero import AEROZ


class TRIM_ZAERO(BaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    """
    type = 'TRIM_ZAERO'
    _field_map = {
        1: 'sid', 2: 'mkaeroz', 3: 'q',
        9: 'wtmass', 10: 'weight',
        17: 'true_g',
    }

    def __init__(self, sid: int, mkaeroz: int, q: float,
                 trimobj_id: int, trimcon_id: int,
                 weight: float, dcg: list[float], inertia: list[float], true_g: str,
                 nxyz: list[float], pqr_dot: list[float], loadset: Optional[int],
                 trimvar_ids: list[int], uxs: list[float],
                 wtmass: float=1.0, comment: str=''):
        """
        Creates a TRIM card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        q : float
            dynamic pressure
        true_g : str
            'TRUE': nxyz given in model units
            'G':    nxyz given in g's
        nxyz : list[float]
            g loading in xyz directions
        pqr_dot : list[float]
            [roll_accel, pitch_accel, yaw_accel]
        loadset : int
            Identification number of a SET1 or SETADD bulk data card that
            specifies a set of identification numbers of TRIMFNC or
            TRIMADD bulk data card.  All values of the trim functions
            defined by the TRIMFNC or TRIMADD bulk data card are computed
            and printed out.
        trimvar_ids : list[str]
            points to a TRIMVAR
            names of the fixed variables; TODO: why are these integers???
        uxs : list[float]
            values corresponding to labels
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Trim set identification number. (Integer > 0)
        self.sid = sid
        self.mkaeroz = mkaeroz
        #: Dynamic pressure. (Real > 0.0)
        self.q = q

        self.trimobj_id = trimobj_id
        self.trimcon_id = trimcon_id

        self.weight = weight
        self.wtmass = wtmass
        self.dcg = np.asarray(dcg)
        self.inertia = np.asarray(inertia)

        self.nxyz = nxyz
        self.true_g = true_g
        self.pqr_dot = pqr_dot
        self.loadset = loadset

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        # points to a TRIMVAR
        self.trimvar_ids = trimvar_ids
        assert len(trimvar_ids) > 0, trimvar_ids

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = uxs

        self.mkaeroz_ref = None
        self.trimvar_refs = None

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
        # TRIM TRIMID IDMK   QINF IDOBJ IDCONS RHOX RHOY RHOZ
        #      WTMASS WEIGHT IXX  IXY   IYY    IXZ  IYZ  IZZ
        #      TRNACC NX     NY   NZ    PDOT   QDOT RDOT LOADSET
        #      IDVAR1 VAL1 IDVAR2 VAL2
        sid = integer(card, 1, 'sid')
        mkaeroz = integer(card, 2, 'mkaeroz')
        qinf = double(card, 3, 'dynamic_pressure')

        # over-determined problem
        #       id    mkaero q     ???
        # TRIM, 100,  101,   42.,  ALPHA, 5.,   0.,   0., 0., 0.,
        #       1e4,  1e3,   1e3,  1e5,   1e3,  1e3,  1e4,
        #       TRUE, FREE, NONE, 32.,   FREE, NONE, 42., None, 17, 1.0
        trimobj_id = integer_or_blank(card, 4, 'trimobj_id', default=0)
        trimcon_id = integer_or_blank(card, 5, 'trimcon_id', default=0)

        dcg = [
            double(card, 6, 'dcg-x'),
            double(card, 7, 'dcg-y'),
            double(card, 8, 'dcg-z'),
        ]

        wtmass = double(card, 9, 'wtmass')
        weight = double(card, 10, 'weight')
        inertia = [
            double(card, 11, 'Ixx'),
            double(card, 12, 'Ixy'),
            double(card, 13, 'Iyy'),
            double(card, 14, 'Ixz'),
            double(card, 15, 'Iyz'),
            double(card, 16, 'Izz'),
        ]
        #  TRUE/G  NX      NY      NZ      PDOT    QDOT    RDOT    LOADSET
        true_g = string(card, 17, 'TRUE/G')
        nx = double_or_string(card, 18, 'NX')
        ny = double_or_string(card, 19, 'NY')
        nz = double_or_string(card, 20, 'NZ')
        nxyz = [nx, ny, nz]

        pdot = double_or_string(card, 21, 'P')
        qdot = double_or_string(card, 22, 'Q')
        rdot = double_or_string(card, 23, 'R')
        pqr_dot = [pdot, qdot, rdot]
        loadset = integer_or_blank(card, 24, 'loadset')

        uxs = []
        trimvar_ids = []

        i = 25
        n = 1
        # print(f'trim_id={sid} ncard={len(card)}')
        while i < len(card):
            trimvar_id = integer(card, i, f'label{n:d}')
            ux = double_or_string(card, i + 1, f'ux{n:d}')
            if isinstance(ux, str):
                assert ux == 'FREE', 'ux=%r' % ux
            # print('  label=%s ux=%s' % (trimvar_id, ux))
            trimvar_ids.append(trimvar_id)
            uxs.append(ux)
            i += 2
            n += 1
        # print(f'trimvar_ids = {trimvar_ids}')
        assert len(card) >= 25, f'len(TRIM card) = {len(card):d}\ncard={card}'
        assert len(trimvar_ids) > 0, trimvar_ids
        assert len(uxs) > 0, uxs
        return TRIM_ZAERO(sid, mkaeroz, qinf,
                          trimobj_id, trimcon_id,
                          weight, dcg, inertia,
                          true_g, nxyz, pqr_dot, loadset,
                          trimvar_ids, uxs, wtmass=wtmass, comment=comment)

    def validate(self):
        assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

        assert isinstance(self.nxyz[0], float) or self.nxyz[0] in ['FREE', 'NONE'], 'nx=%r' % self.nxyz[0]
        assert isinstance(self.nxyz[1], float) or self.nxyz[1] in ['FREE', 'NONE'], 'ny=%r' % self.nxyz[1]
        assert isinstance(self.nxyz[2], float) or self.nxyz[2] in ['FREE', 'NONE'], 'nz=%r' % self.nxyz[2]

        assert isinstance(self.pqr_dot[0], float) or self.pqr_dot[0] in ['FREE', 'NONE'], 'pdot=%r' % self.pqr_dot[0]
        assert isinstance(self.pqr_dot[1], float) or self.pqr_dot[1] in ['FREE', 'NONE'], 'qdot=%r' % self.pqr_dot[1]
        assert isinstance(self.pqr_dot[2], float) or self.pqr_dot[2] in ['FREE', 'NONE'], 'rdot=%r' % self.pqr_dot[2]

        assert self.q > 0.0, 'q=%s\n%s' % (self.q, str(self))
        if len(set(self.trimvar_ids)) != len(self.trimvar_ids):
            msg = 'not all labels are unique; labels=%s' % str(self.trimvar_ids)
            raise RuntimeError(msg)
        if len(self.trimvar_ids) != len(self.uxs):
            msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                len(self.trimvar_ids), len(self.uxs), str(self.trimvar_ids), str(self.uxs))
            raise RuntimeError(msg)

    def cross_reference(self, model: BDF) -> None:
        zaero = model.zaero
        self.mkaeroz_ref = zaero.mkaeroz[self.mkaeroz]

        assert self.trimobj_id == 0, self.trimobj_id
        assert self.trimcon_id == 0, self.trimcon_id
        weight_unit = '???'
        mass_unit = '???'
        length_unit = '???'
        inertia_unit = '???'
        aeroz = model.aeros
        if aeroz is not None:
            mass_unit = aeroz.fm_mass_unit
            weight_unit = aeroz.weight_unit
            length_unit = aeroz.fm_length_unit

        pressure_unit = f'{weight_unit}/{length_unit}^2'
        if self.wtmass == 1.0:
            inertia_unit = f'{mass_unit}*{length_unit}^2'
            weight_unit = mass_unit
        else:
            inertia_unit = f'{weight_unit}*{length_unit}^2'

        aeroz: AEROZ = model.aeros
        ref = aeroz.xyz_ref
        cg = self.dcg + ref
        mach = self.mkaeroz_ref.mach
        msg = (
            f'trim_id = {self.sid}\n'
            f'  mach={mach:g}; q={self.q} ({pressure_unit})\n'
            f'  weight={self.weight:g} ({weight_unit})\n'
            f'  mass={self.weight*self.wtmass:g} ({mass_unit})\n'
            f'  ref={ref} ({length_unit}); per AEROZ\n'
            f'  dcg={self.dcg} ({length_unit}); per TRIM\n'
            f'  cg={cg} ({length_unit})\n'
            f'  inertia={self.inertia*self.wtmass} ({inertia_unit})\n'
            f'  true/g={self.true_g}\n\n'
            f'  nxyz={self.nxyz}\n'
            f'  pqr_dot={self.pqr_dot}\n'
        )
        free_variables = []
        fixed_variables = []
        linked_variables = []
        rb_state_variables = []
        trim_dofs = []
        trim_variables = []
        for name, nxyzi in zip(('NX', 'NY', 'NZ'), self.nxyz):
            if isinstance(nxyzi, str):
                if nxyzi == 'FREE':
                    msg += f'Trim DOF (Free): {name} = {nxyzi}\n'
                    trim_dofs.append(name)
                    free_variables.append(name)
                elif nxyzi == 'NONE':
                    msg += f'Trim DOF (N/A): {name} = {nxyzi}\n'
                    pass
                else:
                    raise RuntimeError(f'{name}={nxyzi} is not [FREE, NONE]')
            elif isinstance(nxyzi, float):
                msg += f'Trim DOF (Fixed): {name} = {nxyzi}\n'
                trim_dofs.append(f'{name}={nxyzi}')
                fixed_variables.append(f'{name}={nxyzi}')
            else:
                raise RuntimeError(f'{name}={nxyzi!r} is not a valid type; type={type(nxyzi)}')

        for name, pqrdi in zip(('PDOT', 'QDOT', 'RDOT'), self.pqr_dot):
            # msg += f'Trim DOF: {name} = {pqrdi}\n'
            if isinstance(pqrdi, str):
                if pqrdi == 'FREE':
                    msg += f'Trim DOF (Free): {name} = {pqrdi}\n'
                    trim_dofs.append(name)
                    free_variables.append(name)
                elif pqrdi == 'NONE':
                    msg += f'Trim DOF (N/A): {name} = {pqrdi}\n'
                else:
                    raise RuntimeError(f'{name}={pqrdi!r} is not [FREE, NONE]')
            elif isinstance(pqrdi, float):
                msg += f'Trim DOF (Fixed): {name} = {pqrdi}\n'
                trim_dofs.append(f'{name}={pqrdi}')
                fixed_variables.append(f'{name}={pqrdi}')
            else:
                raise RuntimeError(f'{name}={pqrdi!r} is not FREE')

        assert len(self.trimvar_ids)
        assert len(self.trimvar_ids) == len(self.uxs)
        trimvar_refs = get_trimvars(model, self.trimvar_ids)
        assert len(trimvar_refs)

        # print(self.get_stats())
        aesurfz = [key for key in model.aesurf]
        aesurfz.sort()
        set_aesurfz = set(aesurfz)
        zona_state_vars = {'ALPHA', 'BETA', 'PRATE', 'QRATE', 'RRATE', 'THKCAM'}
        for trimvar_id, trimvar_ref, ux in zip(self.trimvar_ids, trimvar_refs, self.uxs):
            # print(trimvar_ref.get_stats())
            if trimvar_ref is None:
                free_variables.append('trimvar_id=???')
                trim_variables.append('trimvar_id=???')
                msg += f'Trim Variable: ???={ux} ({trimvar_id})\n'
                continue

            label = trimvar_ref.label
            trimlnk_id = trimvar_ref.trimlnk_id
            if trimlnk_id == 0:
                linkage = ''
            else:
                if trimvar_ref.trimlnk_ref is None:
                    linkage = f'; linked/missing ({trimlnk_id})'
                else:
                    linkage = f'; linked ({trimlnk_id})'
                linked_variables.append(f'{label}={trimlnk_id}')

            if isinstance(ux, str):
                if ux == 'FREE':
                    trim_variables.append(label)
                    free_variables.append(label)
                    if label in zona_state_vars:
                        rb_state_variables.append(f'{label}')
                        msg += f'Trim Variable (Free RB State): {label!r}={ux} ({trimvar_id}){linkage}\n'
                        continue
                    if label not in zona_state_vars and label in set_aesurfz:
                        set_aesurfz.remove(label)
                    msg += f'Trim Variable (Free): {label!r}={ux} ({trimvar_id}){linkage}\n'
                elif ux == 'NONE':
                    msg += f'Trim DOF (N/A): {label!r}={ux} ({trimvar_id}){linkage}\n'
                else:
                    raise RuntimeError(f'{label!r}={ux} ({trimvar_id}) is not [FREE, NONE]')
            elif isinstance(ux, float):
                msg += f'Trim Variable (Fixed RB State): {label!r}={ux} ({trimvar_id}){linkage}\n'
                rb_state_variables.append(f'{label}={ux}')
                # fixed_variables.append(f'{label}={ux}')  # TODO: not sure; yes
                trim_variables.append(f'{label}={ux}')
                if label not in zona_state_vars and label in set_aesurfz:
                    set_aesurfz.remove(label)
            else:
                raise RuntimeError(f'{label!r}={ux} ({trimvar_id}) is not [FREE, NONE]')
        self.trimvar_refs = trimvar_refs

        nfree = len(free_variables)
        nfixed = len(fixed_variables)
        nlinked = len(linked_variables)
        nstate = len(rb_state_variables)
        nunused_aesurfz = len(set_aesurfz)
        msg += f'\nSummary:\n'
        msg += f'  trim_dofs = {trim_dofs}; ntrim_dof={len(trim_dofs)}\n'
        msg += f'  trim_variables = {trim_variables}; ntrimvar={len(trim_variables)}\n\n'

        msg += f'  fixed_variables    = {fixed_variables}; nfixed={nfixed}\n'
        msg += f'  free_variables     = {free_variables}; nfree={nfree}\n\n'
        msg += f'  rb_state_variables = {rb_state_variables}; nstate={nstate}\n'
        msg += f'  aesurfz = {aesurfz}; n={len(aesurfz)}\n\n'
        msg += f'  linked_variables = {linked_variables}; nlinked={nlinked}\n'
        msg += f'  unused_aesurfz = {list(set_aesurfz)}; n={nunused_aesurfz} (should be linked, optimized, or unused)\n'
        # ndelta1 = nfree - (nfixed + nlinked)
        # msg += f'ndelta1 = nfree - (nfixed + nlinked) = {nfree} - ({nfixed} + {nlinked}) = {nfree} - {nfixed + nlinked} = {ndelta1}\n'

        ndelta = nfree - nfixed
        if ndelta == 0:
            determination = 'solvable'
        elif ndelta > 0:
            determination = 'over-determined; reduce ndelta'
        else:
            determination = 'under-determined; increase ndelta'
        msg += f'ndelta = nfree - nfixed = {nfree} - {nfixed} = {ndelta} ({determination})\n'
        # print(msg)
        assert nfixed == nfree, msg
        # assert nfixed == nfree, msg
        # assert nlinked == nunused_aesurfz, msg

        # ] + self.cg + [self.wtmass, self.weight] + self.inertia + [
        #     self.true_g] + self.nxyz + self.pqr_dot + [self.loadset]

        # nlabels = len(self.trimvar_ids)
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)
        del xref_errors

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model: BDF):
        mkaeroz_id = self.mkaeroz
        mkaeroz = model.zaero.mkaeroz[mkaeroz_id]
        #print(mkaeroz)
        mach = mkaeroz.mach
        labels = []
        uxs = []
        comment = str(self)
        assert len(self.trimvar_ids) == len(self.uxs)
        for trimvar_id, ux in zip(self.trimvar_ids, self.uxs):
            if ux != 'FREE':
                trimvar = model.zaero.trimvar[trimvar_id]
                label = trimvar.label
                assert isinstance(label, str), 'label=%r' % label
                comment += str(trimvar)
                labels.append(label)
                uxs.append(ux)

        assert self.q is not None
        if self.q == 'NONE':
            self.q = 1.
        assert isinstance(self.q, float), str(self)
        trim = TRIM(self.sid, mach, self.q, labels, uxs,
                    aeqr=1.0, comment=comment)
        trim.validate()
        return trim

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'TRIM', self.sid, self.mkaeroz, self.q, self.trimobj_id, self.trimcon_id,
        ] + list(self.dcg) + [self.wtmass, self.weight] + list(self.inertia) + [
            self.true_g] + self.nxyz + self.pqr_dot + [self.loadset]

        nlabels = len(self.trimvar_ids)
        assert nlabels > 0, self.trimvar_ids
        for (i, label, ux) in zip(count(), self.trimvar_ids, self.uxs):
            list_fields += [label, ux]
        return list_fields

    def repr_fields(self):
        # trimobj_id = set_blank_if_default(self.trimobj_id, 0)
        # trimcon_id = set_blank_if_default(self.trimcon_id, 0)
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIMLNK(BaseCard):
    """
    Defines a set of coefficient and trim variable identification
    number pairs for trim variable linking.

    +=========+========+========+========+========+========+========+========+========+
    |    1    |    2   |    3   |    4   |   5   |    6    |    7   |   8    |    9   |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+
    | TRIMLNK | IDLINK |   SYM  | COEFF1 | IDVAR1 | COEFF2 | IDVAR2 | COEFF3 | IDVAR3 |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+
    |         | COEFF4 | IDVAR4 |  etc.  |        |        |        |        |        |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+
    | TRIMLNK |   10   |   SYM  |   1.0  |   100  |   0.5  |   200  |         |        |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+
    """
    type = 'TRIMLNK'

    def __init__(self, link_id: int, sym: str,
                 coeffs: list[float], var_ids: list[int],
                 comment: str=''):
        """
        Creates a TRIMLNK card

        Parameters
        ----------
        link_id : int
            the TRIMLNK id
        sym : ???
            ???
        coeffs : list[float]
            linkage coefficients
        var_ids : list[int]
            pointer to the TRIMVAR
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.link_id = link_id
        self.sym = sym
        self.coeffs = coeffs
        self.var_ids = var_ids
        assert sym in ['SYM', 'ASYM', 'ANTI'], sym
        self.trimvar_refs = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIMLNK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        link_id = integer(card, 1, 'var_id')
        sym = string_or_blank(card, 2, 'sym')

        nfields = len(card) - 3
        assert nfields % 2 == 0, card
        icoeff = 1
        coeffs = []
        var_ids = []
        for ifield in range(3, len(card), 2):
            coeff = double(card, ifield, 'coeff_%i' % icoeff)
            var_id = integer(card, ifield + 1, 'var_%i' % icoeff)
            coeffs.append(coeff)
            var_ids.append(var_id)
            icoeff += 1
        assert len(card) >= 5, f'len(TRIMLNK card) = {len(card):d}\ncard={card}'
        return TRIMLNK(link_id, sym, coeffs, var_ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        trimvar_refs = []
        zaero = model.zaero
        for var_id in self.var_ids:
            if var_id in zaero.trimvar:
                ref_id = zaero.trimvar[var_id]
            else:
                raise RuntimeError(f'var_id {var_id} not in trimvar')
            trimvar_refs.append(ref_id)
        self.trimvar_refs = trimvar_refs
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)
        del xref_errors

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model):
        label = 'LNK_%s' % self.link_id
        trimvars = model.zaero.trimvar

        comment = str(self)
        independent_labels = []
        for var_id in self.var_ids:
            trimvar = trimvars[var_id]
            label = trimvar.label
            comment += str(trimvar)
            independent_labels.append(label)

        Cis = self.coeffs
        aelink = AELINK(self.link_id, label, independent_labels, Cis, comment=comment)
        aelink.validate()
        return aelink

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIMLNK', self.link_id, self.sym]
        for coeff, var in zip(self.coeffs, self.var_ids):
            list_fields.append(coeff)
            list_fields.append(var)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIMVAR(BaseCard):
    """
    Specifies a trim variable for static aeroelastic trim variables.

    """
    type = 'TRIMVAR'

    def __init__(self, var_id: int, label: str, lower: float, upper: float,
                 trimlnk_id: int, dmi: None, sym: int, initial: Optional[float]=None,
                 dcd: float|str='NONE', dcy: float|str='NONE', dcl: float|str='NONE',
                 dcr: float|str='NONE', dcm: float|str='NONE', dcn: float|str='NONE', comment: str=''):
        """
        Creates a TRIMVAR card for a static aero (144) analysis.

        Parameters
        ----------
        var_id : int
            the trim id; referenced by the Case Control TRIM field
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.var_id = var_id
        self.label = label
        self.lower = lower
        self.upper = upper
        self.trimlnk_id = trimlnk_id
        self.dmi = dmi
        self.sym = sym
        self.initial = initial
        self.dcd = dcd
        self.dcy = dcy
        self.dcl = dcl
        self.dcr = dcr
        self.dcm = dcm
        self.dcn = dcn
        self.trimlnk_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment:  str=''):
        """
        Adds a TRIMVAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        var_id = integer(card, 1, 'var_id')
        label = string(card, 2, 'label')
        lower = double_or_blank(card, 3, 'lower')
        upper = double_or_blank(card, 4, 'upper')
        trimlnk_id = integer_or_blank(card, 5, 'TRIMLNK', default=0)
        dmi = string_or_blank(card, 6, 'DMI', default='')
        sym = string_or_blank(card, 7, 'sym')
        initial = double_or_blank(card, 8, 'initial', default=None)
        dcd = double_string_or_blank(card, 9, 'DCD', default='NONE')
        dcy = double_string_or_blank(card, 10, 'DCY', default='NONE')
        dcl = double_string_or_blank(card, 11, 'DCL', default='NONE')
        dcr = double_string_or_blank(card, 12, 'DCR', default='NONE')
        dcm = double_string_or_blank(card, 13, 'DCM', default='NONE')
        dcn = double_string_or_blank(card, 14, 'DCN', default='NONE')
        return TRIMVAR(var_id, label, lower, upper, trimlnk_id, dmi, sym,
                       initial, dcd, dcy, dcl, dcr, dcm,
                       dcn, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        if self.trimlnk_id:
            self.trimlnk_ref = model.zaero.trimlnk[self.trimlnk_id]
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)
        del xref_errors

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model: BDF):
        raise NotImplementedError()
        #mkaeroz_id = self.mkaeroz
        #mkaeroz = model.zaero.mkaeroz[mkaeroz_id]
        #mach = mkaeroz.mach
        #labels = []
        #uxs = []
        #for label_id, ux in zip(self.labels, self.uxs):
            #if ux != 'FREE':
                #label = model.zaero.trimvar[label_id]
                #labels.append(label)
                #uxs.append(ux)
        #trim = TRIM(self.sid, mach, self.q, labels, uxs,
                    #aeqr=1.0, comment=str(self))
        #trim.validate()
        #return trim

    def raw_fields(self) -> list:
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIMVAR', self.var_id, self.label, self.lower, self.upper,
                       self.trimlnk_id, self.dmi, self.sym, self.initial,
                       self.dcd, self.dcy, self.dcl, self.dcr, self.dcm, self.dcn]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

def get_trimvars(model: BDF, trimvar_ids) -> list[TRIMVAR]:
    trimvar_refs = []
    for trimvar_id in trimvar_ids:
        try:
            trimvar_ref = model.zaero.trimvar[trimvar_id]
        except KeyError as error:
            trimvar_refs.append(None)
            continue
        trimvar_refs.append(trimvar_ref)
    return trimvar_refs
