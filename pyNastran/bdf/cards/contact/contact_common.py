from __future__ import annotations
from typing import Optional, TYPE_CHECKING

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.internal_get import coord_id
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
    string_or_blank, string_choice_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class BCPARA(BaseCard):
    """
    Defines contact parameters used in SOL 600.

    +--------+---------+--------+--------+--------+--------+---------+--------+
    |   1    |    2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +========+=========+========+========+========+========+=========+========+
    | BCPARA |  CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +--------+---------+--------+--------+--------+--------+---------+--------+
    |        | Param4  | Value4 | Param5 | Value5 |  etc.  |         |        |
    +--------+---------+--------+--------+--------+--------+---------+--------+
    | BCPARA | NBODIES |   4    |  BIAS  |   0.5  |        |         |        |
    +--------+---------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BCPARA'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'NBODIES' : 4}
        return BCPARA(csid, params)

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        keys, values = self.params

        self.params = {}
        for key, value in zip(keys, values):
            assert isinstance(key, str), (key, value)
            if isinstance(value, bytes):
                value = value.decode('ascii')
            assert isinstance(value, (int, float, str)), (key, value)
            self.params[key] = value

    def __init__(self, csid: int, params: dict[str, Optional[int | float]],
                 comment: str=''):
        """
        Creates a BCPARA card

        Parameters
        ----------
        csid : int
            ID is not used and should be set to zero. Only one BCPARA should be
            entered and it applies to all subcases.
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : int/float
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer_or_blank(card, 1, 'csid', default=0)
        i = 2
        j = 1
        params = {}
        ivalue_line = 0
        allowed_params = [
            'BIAS', 'DDULMT', 'ERRBAS', 'ERROR', 'FKIND',
            'FNTOL', 'FTYPE', 'IBSEP', 'ISPLIT', 'ICSEP',
            'LINQUAD', 'LINCNT', 'MAXNOD', 'METHOD', 'MAXENT',
            'NBODIES', 'NLGLUE', 'SLDLMT', 'SEGSYM', 'THKOFF',
        ]
        while i < card.nfields:
            param = string_or_blank(card, i, f'param{j}')
            i += 1

            if param is None:
                blank(card, i, f'blank_value{j}')
                i += 1
                j += 1
                continue
            #print('param', param)
            assert param in allowed_params, f'param={param!r} is not supported; allowed={allowed_params}'
            if param == 'BIAS':
                # Contact tolerance bias factor. (0.0 ≤ Real ≤ 1.0.;
                # Default = 0.9 for IGLUE=0, if field left blank or 0.0 (to obtain a near zero value,
                # enter 1.0E-16).
                # Default = 0.0 for IGLUE <>0. Note that when IGLUE<>0, BIAS can only be given
                # by the BCTABLE or BCONPRG.
                # Default = 0.0 for BEHAVE=SYMM on BCBODY, if field left blank or 0.0.
                value = double_or_blank(card, i, f'{param}_value{j}', default=0.0)
            elif param == 'DDULMT':
                # Maximum value of DDU in a Newton-Raphson cycle. (Real ≥ 0.0;
                # Default = 0.0, no limitation)
                value = double_or_blank(card, i, f'{param}_value{j}', default=0.0)
            elif param == 'ERRBAS':
                # Error computation option.
                # Integer
                # 0 = compute error globally or
                # 1 = calculate error based on each pair of slave-master;
                # Default = 0)
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1], f'ERRBAS must be [0, 1]; ERRBAS={value}'
            elif param == 'ERROR':
                # Distance below which a node is considered touching a body. (Real; Default = blank).
                # Automatically calculated if left blank. If left blank, the code calculates ERROR as the
                # smallest value resulting from: Either dividing the smallest nonzero element
                # dimension (plates or solids) in the contact body by 20. Or dividing the thinnest shell
                # thickness in the contact body by 4. This value is then used for all contact pairs.
                value = double_or_blank(card, i, f'{param}_value{j}')
            elif param == 'FKIND':
                # FKIND
                # (2,5)
                # Friction kind. (Integer 0 or 1)
                # 0 Friction based on nodal stress.
                # 1 Default if friction is present and FKIND is not entered.
                # Friction based on nodal force.
                value = integer_or_blank(card, i, f'value{j}', 1)
                assert value in [0, 1], f'FKIND must be [6]; FKIND={value}'
            elif param == 'FNTOL':
                # FNTOL Separation force (or stress if separation is controlled by stress as determined by
                # IBSEP) above which a node separates from a body. Automatically calculated if left
                # blank. (Real; Default = blank)
                value = double_or_blank(card, i, f'{param}_value{j}')

            elif param == 'FTYPE':
                # FTYPE Friction type. See Remark 5. (Integer)
                # 0 No friction. (Default)
                # 6 Bilinear Coulomb friction.
                # 7 Bilinear Shear friction.
                #
                # FTYPE
                # Friction type. (Integer)
                # 0 No friction. (Default)
                # 1 Shear friction.
                # 2 Coulomb Friction.
                # 3 Shear friction for rolling.
                # 4 Coulomb friction for rolling.
                # 5 Stick-slip Coulomb friction.
                # 6 Bilinear Coulomb friction. (Default if friction is present and FTYPE is not entered.)
                # 7 Bilinear Shear friction.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1, 2, 3, 4, 5, 6], f'FTYPE must be [0, 1, 2, 3, 4, 5, 6]; FTYPE={value}'
            elif param == 'IBSEP':
                # Flag for separation based on stresses or forces. (Integer > 0; Default = 0)
                # 0 Separation based on forces.
                # 1 Separation based on absolute stresses (force/area)
                # 2 Separation based on absolute stress (extrapolating integration point stresses)
                # 3 Relative nodal stress (force/area)
                # 4 Separation based on relative stress (extrapolating integration point stresses)
                # Only option 2 and 4 can be used with mid-side node elements where the mid-side
                # nodes contact (LINQUAD=-1). For segment to segment contact, the program will
                # set IBSEP to 2 internally. See Remarks 6. and 10.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1, 2, 3, 4], f'IBSEP must be [0, 1, 2, 3, 4]; IBSEP={value}'
            elif param == 'ISPLIT':
                # ISPLIT (2,7) Flag for increment splitting procedure. (Integer > 0; Default = 3 for
                # statics and 0 for dynamics)
                # 0 Uses increment splitting procedures for the fixed time step
                # procedures.
                # BCPARA 1115
                # Contact Parameters in SOL 600
                # Main Index
                # 1 Suppresses splitting for the fixed time step procedures. Any
                # penetration that occurred in the previous increment is
                # adjusted for equilibrium at the start of the next increment.
                # This method may require smaller time steps then the other
                # methods
                # 2 Suppresses splitting for adaptive time step procedures. Any
                # penetration that occurred in the previous increment is
                # adjusted for equilibrium at the start of the next increment.
                # This method may require smaller time steps then the other
                # methods.
                # 3 Uses contact procedure which does not require increment
                # splitting (3 is not available for dynamics). If a run does not
                # converge due to an excessive number of “iterative
                # penetration checking” messages, ISPLIT=2 may help,
                # however the time steps may need to be set to smaller values.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=3)
                assert value in [0, 1, 2, 3], f'ISPLIT must be [0, 1, 2, 3]; ISPLIT={value}'
            elif param == 'ICSEP':
                # ICSEP Flag to control separation. Not used for segment-to-segment contact. (Integer > 0;
                # Default = 0)
                # 0 The node separates and an iteration occurs if the force on the node is
                # greater than the separation force.
                # 1 If a node which was in contact at the end of the previous increment
                # has a force greater than the separation force, the node does NOT
                # separate in this increment, but separates at the beginning of the next
                # increment.
                # 2 If a new node comes into contact during this increment, it is not
                # allowed to separate during this increment, prevents chattering.
                # 3 Both 1 and 2 are in effect.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [0, 1, 2, 3], f'ICSEP must be [0, 1, 2, 3]; ICSEP={value}'
            elif param == 'LINQUAD':
                # LINQUAD (2,14)
                # Higher order element contact flag (Integer; Default = 1).
                # 1 The outer boundary of a contact body is described by the
                # corner nodes only and mid-side nodes can’t come into
                # contact.
                # -1 The outer boundary is described by a quadratic field and
                # both corner and mid-side nodes are considered in contact.
                # If this flag is set to -1 and IBSEP is blank, IBSEP will be reset
                # to 2. This option is only available with Marc 2003 and
                # subsequent releases.
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [-1, 1], f'LINQUAD must be [-1, 1]; LINQUAD={value}'
            elif param == 'LINCNT':
                value = integer_or_blank(card, i, f'{param}_value{j}', default=0)
                assert value in [-1, 0, 1], f'LINCNT must be [-1, 0, 1]; LINQUAD={value}'
            elif param == 'METHOD':
                #METHOD Flag to select Contact methods. (Character)
                #NODESURF Node to segment contact. (Default)
                #SEGTOSEG Segment to segment contact.
                value = string_choice_or_blank(card, i, f'{param}_value{j}',
                                       ('NODESURF', 'SEGTOSEG', 'SEGSMALL', 'SEGLARGE'),
                                        default='NODESURF')
            elif param == 'MAXENT':
                # MAXENT (2,2) Maximum number of entities created for any contact body. (Integer >
                # 0 or blank; default is max element number or 1.5 times the number of
                # nodes whichever is smaller)
                value = integer_or_blank(card, i, f'{param}_value{j}')
            elif param == 'MAXNOD':
                # MAXNOD (2,3) Maximum number of nodes that lie on the periphery of any
                # deformable contact body. (Integer > 0 or blank; default is the number
                # of nodes)
                value = integer_or_blank(card, i, f'{param}_value{j}')
            elif param == 'NBODIES':
                # NBODIES (2,1) Number of contact bodies defined in the analysis. (Integer > 0 or blank)
                value = integer(card, i, f'{param}_value{j}')
            elif param == 'NLGLUE':
                #(SOLs 101 and 400 only)
                #If all slave's for the BCTABLE or BCONPRG corresponding to the first loadcase
                #(first subcase and first step) contain IGLUE >0, permanent glued contact with small
                #rotation condition will be used for all SLAVE entries in all subcases and all steps
                #unless BCPARA,0,,1 is specified. If IGLUE < 0 exists, permanent glued
                value = integer_or_blank(card, i, f'{param}_value{j}', default=1)
                assert value in [0, 1], f'NLGLUE must be [0, 1]; NLGLUE={value}'
            elif param == 'SLDLMT':
                #SLDLMT Maximum allowed sliding distance, beyond it the contact segments are
                # redefined, for segment to segment contact analysis with large deformation.
                # (Real ≥ 0.0; Default = 0.0) See remark 9.
                value = double_or_blank(card, i, f'{param}_value{j}', default=0.0)
            elif param == 'SEGSYM':
                #Specify symmetric or non-symmetric friction matrix in segment to segment contact
                #analysis. (Integer 0 = symmetric matrix or 1 = non-symmetric matrix; Default = 0)
                value = integer_or_blank(card, i, f'{param}_value{j}', default=1)
                assert value in [0, 1], f'SEGSYM must be [0, 1]; SEGSYM={value}'
            elif param == 'THKOFF':
                #Ignore thickness from the tolerance used by ISEARCH=2 in node-to-surface contact
                #or from the characteristic length (for PENALT and AUGDIST) in segment-tosegment
                #contact. (Integer 0 = do not ignore thickness or 1 = remove thickness;
                #Default = 0)
                value = integer_or_blank(card, i, f'{param}_value{j}', default=1)
                assert value in [0, 1], f'THKOFF must be [0, 1]; THKOFF={value}'
            else:
                raise NotImplementedError(f'param={param} card={card}')

            params[param] = value
            i += 1
            j += 1

            ivalue_line += 1
            if ivalue_line == 3:
                if i == 8:
                    ivalue_line = 0
                    #print('*1', i, j, ivalue_line)
                    i += 1
                else:
                    #print('*2', i, j, ivalue_line)
                    i += 2
        return BCPARA(csid, params, comment=comment)

    def raw_fields(self):
        fields = ['BCPARA', self.csid]
        ivalue = 1
        ivalue_line = 1
        for key, value in sorted(self.params.items()):
            #print(ivalue, key, value)
            fields.append(key)
            fields.append(value)
            if ivalue_line == 3:
                if ivalue == 3:
                    #print('*')
                    fields.append(None)
                else:
                    #print('**')
                    fields.append(None)
                    fields.append(None)
                ivalue_line = 0
            ivalue += 1
            ivalue_line += 1
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCONP(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+----+-------+--------+-----+------+--------+-------+-----+
    |   1   |  2 |   3   |   4    |  5  |   6  |   7    |   8   |  9  |
    +=======+====+=======+========+=====+======+========+=======+=====+
    | BCONP | ID | SLAVE | MASTER |     | SFAC | FRICID | PTYPE | CID |
    +-------+----+-------+--------+-----+------+--------+-------+-----+
    | BCONP | 95 |   10  |   15   |     |  1.0 |   33   |   1   |     |
    +-------+----+-------+--------+-----+------+--------+-------+-----+

    """
    type = 'BCONP'
    @classmethod
    def _init_from_empty(cls):
        contact_id = 1
        slave = 2
        master = 3
        sfac = 4
        friction_id = 5
        ptype = 'cat'
        cid = 0
        return BCONP(contact_id, slave, master, sfac, friction_id, ptype, cid)

    def __init__(self, contact_id, slave, master, sfac, friction_id, ptype, cid, comment: str=''):
        if comment:
            self.comment = comment

        self.contact_id = contact_id
        self.slave = slave
        self.master = master
        self.sfac = sfac
        self.friction_id = friction_id
        self.ptype = ptype
        self.cid = cid
        self.cid_ref = None
        self.friction_id_ref = None
        self.slave_ref = None
        self.master_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        contact_id = integer(card, 1, 'contact_id')
        slave = integer(card, 2, 'slave')
        master = integer(card, 3, 'master')
        sfac = double_or_blank(card, 5, 'sfac', default=1.0)
        friction_id = integer_or_blank(card, 6, 'fric_id')
        ptype = integer_or_blank(card, 7, 'ptype', default=1)
        cid = integer_or_blank(card, 8, 'cid', default=0)
        return BCONP(contact_id, slave, master, sfac, friction_id, ptype, cid,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by BCONP line_id={self.contact_id}'
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)
        if self.friction_id is not None:
            self.friction_id_ref = model.bfric[self.friction_id]
        self.slave_ref = model.blseg[self.slave]
        self.master_ref = model.blseg[self.master]

    def uncross_reference(self) -> None:
        self.cid = self.Cid()
        self.friction_id = self.FrictionId()
        self.cid_ref = None
        self.friction_id_ref = None
        self.slave_ref = None
        self.master_ref = None

    def Cid(self) -> int:
        return coord_id(self.cid_ref, self.cid)

    def FrictionId(self) -> int:
        if self.friction_id_ref is not None:
            return self.friction_id_ref.friction_id
        return self.friction_id

    def Slave(self) -> int:
        if self.slave_ref is not None:
            return self.slave_ref.line_id
        return self.slave

    def Master(self) -> int:
        if self.master_ref is not None:
            return self.master_ref.line_id
        return self.master

    def raw_fields(self):
        list_fields = [
            'BCONP', self.contact_id, self.Slave(), self.Master(), None, self.sfac,
            self.FrictionId(), self.ptype, self.Cid()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
