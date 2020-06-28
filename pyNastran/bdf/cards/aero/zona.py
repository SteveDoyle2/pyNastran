# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * TRIM

All cards are BaseCard objects.

"""
from __future__ import annotations
from itertools import count
from typing import List, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils import object_attributes, object_methods
from pyNastran.utils.numpy_utils import integer_types

from pyNastran.bdf.cards.aero.dynamic_loads import Aero
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    loose_string_or_blank, string_or_blank, double_or_string, blank,
)
from pyNastran.bdf.cards.aero.aero import (Spline, CAERO1, CAERO2, PAERO2, # PAERO1,
                                           SPLINE1, AESURF, AELIST, # SPLINE2, SPLINE3,
                                           AELINK, AEFACT)
from pyNastran.bdf.cards.aero.static_loads import TRIM, AEROS
from pyNastran.bdf.cards.aero.dynamic_loads import AERO # MKAERO1,
from pyNastran.bdf.cards.aero.utils import (
    elements_from_quad, points_elements_from_quad_points, create_ellipse)
from pyNastran.bdf.cards.coordinate_systems import Coord
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class ZONA:
    def __init__(self, model):
        self.model = model
        self.caero_to_name_map = {}

        #: store PANLST1,PANLST2,PANLST3
        self.panlsts = {}
        self.mkaeroz = {}
        self.trimvar = {}
        self.trimlnk = {}
        #: store PAFOIL7/PAFOIL8
        self.pafoil = {}

    @classmethod
    def _init_from_self(cls, model):
        """helper method for dict_to_h5py"""
        return cls(model)

    def clear(self):
        """clears out the ZONA object"""
        self.panlsts = {}
        self.mkaeroz = {}
        self.trimvar = {}
        self.trimlnk = {}
        self.pafoil = {}
        #self.aeroz = {}

    def object_attributes(self, mode:str='public', keys_to_skip: Optional[List[str]]=None,
                          filter_properties: bool=False):
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : List[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong
        """
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            'log', 'model',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public',
                       keys_to_skip: Optional[List[str]]=None) -> List[str]:
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : List[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong
        """
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []  # type: List[str]

        my_keys_to_skip = ['log',]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def verify(self, xref):
        if self.model.nastran_format != 'zona':
            return
        for panlst in self.panlsts.values():
            panlst._verify(xref)
        for mkaeroz in self.mkaeroz.values():
            mkaeroz._verify(xref)
        for trimvar in self.trimvar.values():
            trimvar._verify(xref)
        for trimlnk in self.trimlnk.values():
            trimlnk._verify(xref)
        for pafoil in self.pafoil.values():
            pafoil._verify(xref)

    def validate(self):
        if self.model.nastran_format != 'zona':
            return
        for panlst in self.panlsts.values():
            panlst.validate()
        for mkaeroz in self.mkaeroz.values():
            mkaeroz.validate()
        for trimvar in self.trimvar.values():
            trimvar.validate()
        for trimlnk in self.trimlnk.values():
            trimlnk.validate()
        for pafoil in self.pafoil.values():
            pafoil.validate()

    def PAFOIL(self, pid, msg=''):
        """gets a pafoil profile (PAFOIL7/PAFOIL8)"""
        try:
            return self.pafoil[pid]
        except KeyError:
            raise KeyError('pid=%s not found%s.  Allowed pafoils=%s'
                           % (pid, msg, np.unique(list(self.pafoil.keys()))))

    def update_for_zona(self):
        card_parser = self.model._card_parser
        card_parser['TRIM'] = (TRIM_ZONA, self.model._add_trim_object)
        card_parser['CAERO7'] = (CAERO7, self.model._add_caero_object)
        card_parser['AEROZ'] = (AEROZ, self.model._add_aeros_object)
        card_parser['AESURFZ'] = (AESURFZ, self._add_aesurfz_object)
        card_parser['FLUTTER'] = (FLUTTER_ZONA, self.model._add_flutter_object)
        card_parser['SPLINE1'] = (SPLINE1_ZONA, self.model._add_spline_object)
        card_parser['SPLINE2'] = (SPLINE2_ZONA, self.model._add_spline_object)
        card_parser['SPLINE3'] = (SPLINE3_ZONA, self.model._add_spline_object)
        card_parser['PANLST1'] = (PANLST1, self._add_panlst_object)
        card_parser['PANLST3'] = (PANLST3, self._add_panlst_object)
        card_parser['PAFOIL7'] = (PAFOIL7, self._add_pafoil_object)
        card_parser['MKAEROZ'] = (MKAEROZ, self._add_mkaeroz_object)
        card_parser['SEGMESH'] = (SEGMESH, self.model._add_paero_object)
        card_parser['BODY7'] = (BODY7, self.model._add_caero_object)
        card_parser['ACOORD'] = (ACOORD, self.model._add_coord_object)
        card_parser['TRIMVAR'] = (TRIMVAR, self._add_trimvar_object)
        card_parser['TRIMLNK'] = (TRIMLNK, self._add_trimlnk_object)
        cards = [
            'CAERO7', 'AEROZ', 'AESURFZ', 'PANLST1', 'PANLST3', 'PAFOIL7',
            'SEGMESH', 'BODY7', 'ACOORD', 'MKAEROZ',
            'TRIMVAR', 'TRIMLNK', 'FLUTTER']
        self.model.cards_to_read.update(set(cards))

    def _add_panlst_object(self, panlst: Union[PANLST1, PANLST3]) -> None:
        """adds an PANLST1/PANLST2/PANLST3 object"""
        assert panlst.eid not in self.panlsts
        assert panlst.eid > 0
        key = panlst.eid
        self.panlsts[key] = panlst
        self.model._type_to_id_map[panlst.type].append(key)

    def _add_pafoil_object(self, pafoil: PAFOIL7) -> None:
        """adds an PAFOIL7/PAFOIL8 object"""
        assert pafoil.pid not in self.pafoil
        assert pafoil.pid > 0
        key = pafoil.pid
        self.pafoil[key] = pafoil
        self.model._type_to_id_map[pafoil.type].append(key)

    def _add_aesurfz_object(self, aesurf: AESURFZ) -> None:
        """adds an AESURFZ object"""
        key = aesurf.aesid
        model = self.model
        assert key not in model.aesurf, '\naesurf=\n%s old=\n%s' % (
            aesurf, model.aesurf[key])
        model.aesurf[key] = aesurf
        model._type_to_id_map[aesurf.type].append(key)

    def _add_mkaeroz_object(self, mkaeroz: MKAEROZ) -> None:
        """adds an MKAEROZ object"""
        assert mkaeroz.sid not in self.mkaeroz
        assert mkaeroz.sid > 0
        key = mkaeroz.sid
        self.mkaeroz[key] = mkaeroz
        self.model._type_to_id_map[mkaeroz.type].append(key)

    def _add_trimvar_object(self, trimvar: TRIMVAR) -> None:
        """adds an TRIMVAR object"""
        assert trimvar.var_id not in self.trimvar
        assert trimvar.var_id > 0
        key = trimvar.var_id
        self.trimvar[key] = trimvar
        self.model._type_to_id_map[trimvar.type].append(key)

    def _add_trimlnk_object(self, trimlnk: TRIMLNK) -> None:
        """adds an TRIMLNK object"""
        assert trimlnk.link_id not in self.trimlnk
        assert trimlnk.link_id > 0
        key = trimlnk.link_id
        self.trimlnk[key] = trimlnk
        self.model._type_to_id_map[trimlnk.type].append(key)

    def cross_reference(self):
        if self.model.nastran_format != 'zona':
            return
        for mkaeroz in self.mkaeroz.values():
            mkaeroz.cross_reference(self.model)
        for trimvar in self.trimvar.values():
            trimvar.cross_reference(self.model)
        for trimlnk in self.trimlnk.values():
            trimlnk.cross_reference(self.model)
        for unused_id, pafoil in self.pafoil.items():
            pafoil.cross_reference(self.model)
        #for aeroz in self.aeroz.values():
            #aeroz.cross_reference(self.model)

        for caero in self.model.caeros.values():
            #print('%s uses CAERO eid=%s' % (caero.label, caero.eid))
            self.caero_to_name_map[caero.label] = caero.eid

    def safe_cross_reference(self):
        self.cross_reference()

    def write_bdf(self, bdf_file, size=8, is_double=False):
        #if self.model.nastran_format != 'zona':
            #return
        for unused_id, panlst in self.panlsts.items():
            bdf_file.write(panlst.write_card(size=size, is_double=is_double))
        for unused_id, mkaeroz in self.mkaeroz.items():
            bdf_file.write(mkaeroz.write_card(size=size, is_double=is_double))
        for unused_id, trimvar in self.trimvar.items():
            bdf_file.write(trimvar.write_card(size=size, is_double=is_double))
        for unused_id, trimlnk in self.trimlnk.items():
            bdf_file.write(trimlnk.write_card(size=size, is_double=is_double))
        for unused_id, pafoil in self.pafoil.items():
            bdf_file.write(pafoil.write_card(size=size, is_double=is_double))

    def convert_to_nastran(self, save=True):
        """Converts a ZONA model to Nastran"""
        if self.model.nastran_format != 'zona':
            caeros = {}
            caero2s = []
            make_paero1 = False
            return caeros, caero2s, make_paero1

        caeros, caero2s, make_paero1 = self._convert_caeros()
        splines = self._convert_splines()
        aesurf, aelists = self._convert_aesurf_aelist()

        trims = self._convert_trim()
        aeros, aero = self.model.aeros.convert_to_zona(self.model)

        aelinks = self._convert_trimlnk()

        if save:
            self.clear()
            self.model.splines = splines
            self.model.aesurf = aesurf
            self.model.aelists = aelists
            self.model.aelinks = aelinks
            self.model.trims = trims
            self.model.aeros = aeros
            self.model.aero = aero
        return caeros, caero2s, make_paero1

    def _convert_caeros(self):
        """Converts ZONA CAERO7/BODY7 to CAERO1/CAERO2"""
        model = self.model
        caeros = {}
        caero2s = []
        make_paero1 = False
        for caero_id, caero in sorted(model.caeros.items()):
            if caero.type == 'CAERO7':
                caero_new = caero.convert_to_nastran()
                make_paero1 = True
            elif caero.type == 'BODY7':
                caero2s.append(caero)
                continue
            else:
                raise NotImplementedError(caero)
            caeros[caero_id] = caero_new

        self._add_caero2s(caero2s, add=False)
        return caeros, caero2s, make_paero1

    def _add_caero2s(self, caero2s, add=False):
        """Converts ZONA BODY7 to CAERO2/PAERO2/AEFACT"""
        model = self.model
        caero_body_ids = []
        for caero2 in caero2s:
            caero_id = caero2.eid
            out = caero2.convert_to_nastran(model)
            caero_new, paero2, aefact_xs, aefact_width, aefact_theta1, aefact_theta2 = out
            caero_body_ids.append(caero_id)
            if add:
                model._add_aefact_object(aefact_xs)
                model._add_aefact_object(aefact_width)
                model._add_aefact_object(aefact_theta1)
                model._add_aefact_object(aefact_theta2)
                model._add_paero_object(paero2)
                model._add_caero_object(caero_new)
        return

    def _convert_splines(self):
        """Converts ZONA splines to splines"""
        splines = {}
        for unused_spline_id, spline in self.model.splines.items():
            #print(spline)
            if spline.type == 'SPLINE1_ZONA':
                splines_new = spline.convert_to_nastran(self.model)
            elif spline.type == 'SPLINE3_ZONA':
                splines_new = spline.convert_to_nastran(self.model)
            else:
                raise NotImplementedError(spline)
            for spline_new in splines_new:
                splines[spline.eid] = spline_new
        return splines

    def _convert_aesurf_aelist(self):
        """
        Converts ZONA AESURFZ to AESURF/AELIST

        +---------+--------+-------+-------+-------+--------+--------+
        |    1    |   2    |   3   |   4   |   5   |   6    |    7   |
        +=========+========+=======+=======+=======+========+========+
        | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |
        +---------+--------+-------+-------+-------+--------+--------+
        | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |
        +---------+--------+-------+-------+-------+--------+--------+
        """
        model = self.model
        aelist_id = max(model.aelists) + 1 if model.aelists else 1
        aesurf_id = aelist_id
        aesurf = {}
        aelists = {}
        for unused_aesurf_name, aesurfi in sorted(model.aesurf.items()):
            aelist, aesurfi2 = aesurfi.convert_to_nastran(model, aesurf_id, aelist_id)
            aelists[aelist.sid] = aelist
            aesurf[aesurfi2.aesid] = aesurfi2
            aesurf_id += 1
            aelist_id += 1
        return aesurf, aelists

    def _convert_trim(self):
        """Converts ZONA TRIM to TRIM"""
        trims = {}
        model = self.model
        for trim_id, trim in sorted(model.trims.items()):
            trim_new = trim.convert_to_nastran(model)
            trims[trim_id] = trim_new
        return trims

    def _convert_trimlnk(self):
        """Converts ZONA TRIMLNK to AELINK"""
        model = self.model
        assert isinstance(model.aelinks, dict), model.aelinks
        aelinks = {}
        for trim_id, trimlnk in sorted(self.trimlnk.items()):
            aelink = trimlnk.convert_to_nastran(model)
            aelinks[trim_id] = aelink
        return aelinks

    def __repr__(self):
        msg = '<ZONA>; nPANLSTs=%s nmkaeroz=%s' % (
            len(self.panlsts), len(self.mkaeroz),
        )
        return msg


class ACOORD(Coord):  # not done
    """
    Defines a general coordinate system using three rotational angles as
    functions of coordinate values in the reference coordinate system.
    The CORD3G entry is used with the MAT9 entry to orient material principal
    axes for 3-D composite analysis.

    +--------+-----+--------+--------+--------+-------+-------+--------+
    |    1   |  2  |    3   |    4   |    5   |   6   |   7   |    8   |
    +========+=====+========+========+========+=======+=======+========+
    | ACOORD |  ID | XORIGN | YORIGN | ZORIGN | DELTA | THETA |        |
    +--------+-----+--------+--------+--------+-------+-------+--------+
    | ACOORD |  10 |  250.0 |  52.5  |  15.0  |  0.0  |  0.0  |        |
    +--------+-----+--------+--------+--------+-------+-------+--------+
    """
    type = 'ACOORD'
    Type = 'R'
    @property
    def rid(self):
        return None

    def __init__(self, cid, origin, delta, theta, comment=''):
        """
        Defines the CORD3G card

        Parameters
        ----------
        cid : int
            coordinate system id
        origin : List[float]
            the xyz origin
        delta : float
            pitch angle
        theta : float
            roll angle
        comment : str; default=''
            a comment for the card

        """
        Coord.__init__(self)
        if comment:
            self.comment = comment
        self.cid = cid
        self.origin = origin
        self.delta = delta
        self.theta = theta

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a ACOORD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        cid = integer(card, 1, 'cid')
        origin_x = double(card, 2, 'origin_x')
        origin_y = double(card, 3, 'origin_y')
        origin_z = double(card, 4, 'origin_z')
        origin = [origin_x, origin_y, origin_z]
        delta = double(card, 5, 'delta')
        theta = double(card, 6, 'theta')
        assert len(card) <= 7, 'len(ACOORD card) = %i\ncard=%s' % (len(card), card)
        return ACOORD(cid, origin, delta, theta, comment=comment)

    def setup(self):
        self.i = np.array([1., 0., 0.])
        self.j = np.array([0., 1., 0.])
        self.k = np.array([0., 0., 1.])

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def coord_to_xyz(self, p):
        return p
        #return self.acoord_transform_to_global(p)

    def acoord_transform_to_global(self, p):
        """
        Parameters
        ----------
        p : (3,) float ndarray
            the point to transform

        .. warning:: not done, just setting up how you'd do this
        .. note::    per http://en.wikipedia.org/wiki/Euler_angles

         "This means for example that a convention named (YXZ) is the result
         of performing first an intrinsic Z rotation, followed by X and
         Y rotations, in the moving axes (Note: the order of multiplication
         of matrices is the opposite of the order in which they're
         applied to a vector)."

        """
        ct = np.cos(np.radians(self.theta))
        st = np.sin(np.radians(self.theta))
        #if rotation == 1:
            #p = self.rotation_x(ct, st) @ p
        #elif rotation == 2:
        p = self.rotation_y(ct, st) @ p
        #elif rotation == 3:
            #p = self.rotation_z(ct, st) @ p
        #else:
            #raise RuntimeError('rotation=%s rotations=%s' % (rotation, rotations))
        return p

    def rotation_x(self, ct, st):
        matrix = np.array([[1., 0., 0.],
                           [ct, 0., -st],
                           [-st, 0., ct]])
        return matrix

    def rotation_y(self, ct, st):
        matrix = np.array([[ct, 0., st],
                           [0., 1., 0.],
                           [-st, 0., ct]])
        return matrix

    def rotation_z(self, ct, st):
        matrix = np.array([[ct, st, 0.],
                           [-st, ct, 0.],
                           [0., 0., 1.]])
        return matrix

    def raw_fields(self):
        list_fields = ['ACOORD', self.cid] + self.origin + [self.delta, self.theta]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESURFZ(BaseCard):
    """
    Specifies an aerodynamic control surface for aeroservoelastic, static
    aeroelastic/trim analysis, or the transient response analysis.

    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    |    1    |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
    +=========+========+=======+=======+=======+========+========+========+========+
    | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |        |        |
    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |        |        |
    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    """
    type = 'AESURFZ'

    @property
    def aesid(self):
        return self.label
    @property
    def alid1_ref(self):
        return None

    def __init__(self, label, surface_type, cid, panlst, setg, actuator_tf,
                 comment=''):
        """
        Creates an AESURF card, which defines a control surface

        Parameters
        ----------
        label : str
            controller name
        surface_type : str
            defines the control surface type {SYM, ASYM}
        cid : int
            coordinate system id to define the hinge axis
        panlst : int
            aero panels defined by PANLST
        setg : int
           ???
        actuator_tf : int
            ???
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Controller name.
        self.label = label
        self.surface_type = surface_type
        self.panlst = panlst
        self.setg = setg
        self.actuator_tf = actuator_tf

        #: Identification number of a rectangular coordinate system with a
        #: y-axis that defines the hinge line of the control surface
        #: component.
        self.cid = cid

        self.cid_ref = None
        self.panlst_ref = None
        self.aero_element_ids = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        label = string(card, 1, 'label')
        surface_type = string(card, 2, 'TYPE')
        cid = integer(card, 3, 'CID')
        panlst = integer(card, 4, 'PANLST/SETK') # PANLST1, PANLST2, PANLST3
        setg = integer(card, 5, 'SETG') # SET1, SETADD
        actuator_tf = integer_or_blank(card, 6, 'ACTID') # ACTU card
        assert len(card) <= 7, 'len(AESURFZ card) = %i\ncard=%s' % (len(card), card)
        assert surface_type in ['SYM', 'ANTISYM', 'ASYM']
        return AESURFZ(label, surface_type, cid, panlst, setg, actuator_tf, comment=comment)

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def SetK(self):
        if self.panlst_ref is not None:
            return self.panlst_ref.eid
        return self.panlst

    #def aelist_id1(self):
        #if self.alid1_ref is not None:
            #return self.alid1_ref.sid
        #return self.alid1

    #def aelist_id2(self):
        #if self.alid2_ref is not None:
            #return self.alid2_ref.sid
        #return self.alid2

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.cid_ref = model.Coord(self.cid)

        #self.alid1_ref = model.AELIST(self.alid1)
        #if self.alid2:
            #self.alid2_ref = model.AELIST(self.alid2)
        #if self.tqllim is not None:
            #self.tqllim_ref = model.TableD(self.tqllim)
        #if self.tqulim is not None:
            #self.tqulim_ref = model.TableD(self.tqulim)
        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.cross_reference(model)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by AESURF aesid=%s' % self.aesid
        self.cid_ref = model.safe_coord(self.cid, self.aesid, xref_errors, msg=msg)
        #if self.cid2 is not None:
            #self.cid2_ref = model.safe_coord(self.cid2, self.aesid, xref_errors, msg=msg)
        #try:
            #self.alid1_ref = model.AELIST(self.alid1)
        #except KeyError:
            #pass
        #if self.alid2:
            #try:
                #self.alid2_ref = model.AELIST(self.alid2)
            #except KeyError:
                #pass
        #if self.tqllim is not None:
            #try:
                #self.tqllim_ref = model.TableD(self.tqllim)
            #except KeyError:
                #pass
        #if self.tqulim is not None:
            #try:
                #self.tqulim_ref = model.TableD(self.tqulim)
            #except KeyError:
                #pass
        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.cross_reference(model)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid = self.Cid()
        self.cid_ref = None

        self.panlst = self.SetK()
        self.panlst_ref = None

    def convert_to_nastran(self, model, aesurf_id, aelist_id):
        """
        +--------+--------+-------+-------+-------+--------+--------+--------+--------+
        |    1   |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
        +========+========+=======+=======+=======+========+========+========+========+
        | AESURF |   ID   | LABEL | CID1  | ALID1 |  CID2  | ALID2  |  EFF   |  LDW   |
        +--------+--------+-------+-------+-------+--------+--------+--------+--------+
        |        |  CREFC | CREFS | PLLIM | PULIM | HMLLIM | HMULIM | TQLLIM | TQULIM |
        +--------+--------+-------+-------+-------+--------+--------+--------+--------+

        +---------+--------+-------+-------+-------+--------+--------+
        |    1    |   2    |   3   |   4   |   5   |   6    |    7   |
        +=========+========+=======+=======+=======+========+========+
        | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |
        +---------+--------+-------+-------+-------+--------+--------+
        | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |
        +---------+--------+-------+-------+-------+--------+--------+
        """
        assert self.surface_type == 'ASYM', str(self)
        aelist = AELIST(aelist_id, self.aero_element_ids)
        aesurf = AESURF(aesurf_id, self.label, self.cid, aelist_id, cid2=None, alid2=None,
                        eff=1.0, ldw='LDW', crefc=1.0,
                        crefs=1.0, pllim=-np.pi/2.,
                        pulim=np.pi/2., hmllim=None,
                        hmulim=None, tqllim=None,
                        tqulim=None, comment=self.comment)
        aesurf.validate()
        aelist.validate()
        return aelist, aesurf

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fieldsreset_camera[int/float/str]
            the fields that define the card

        """
        list_fields = ['AESURFZ', self.label, self.surface_type, self.cid,
                       self.panlst, self.setg, self.actuator_tf]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        Writes the card with the specified width and precision

        Parameters
        ----------
        size : int (default=8)
            size of the field; {8, 16}
        is_double : bool (default=False)
            is this card double precision

        Returns
        -------
        msg : str
            the string representation of the card

        """
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEROZ(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +-------+-------+-------+------+------+-------+-------+-------+
    |   1   |   2   |   3   |  4   |  5   |   6   |   7   |   8   |
    +=======+=======+=======+======+======+=======+=======+=======+
    | AEROS | ACSID | RCSID | REFC | REFB | REFS  | SYMXZ | SYMXY |
    +-------+-------+-------+------+------+-------+-------+-------+
    | AEROS |   10  |   20  | 10.  | 100. | 1000. |   1   |       |
    +-------+-------+-------+------+------+-------+-------+-------+
    """
    type = 'AEROZ'
    #_field_map = {
        #1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        #6:'symXZ', 7:'symXY',
    #}

    def __init__(self, fm_mass_unit, fm_length_unit,
                 cref, bref, sref,
                 flip='NO', acsid=0, rcsid=0, sym_xz=0, xyz_ref=None, comment=''):
        """
        Creates an AEROZ card

        Parameters
        ----------
        cref : float
            the aerodynamic chord
        bref : float
            the wing span
            for a half model, this should be the full span
            for a full model, this should be the full span
        sref : float
            the wing area
            for a half model, this should be the half area
            for a full model, this should be the full area
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        rcsid : int; default=0
            coordinate system for rigid body motions
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

        self.fm_mass_unit = fm_mass_unit
        self.fm_length_unit = fm_length_unit
        self.flip = flip

        #: Aerodynamic coordinate system identification.
        self.acsid = acsid

        #: Reference coordinate system identification for rigid body motions.
        self.rcsid = rcsid

        #: Reference chord length
        self.cref = cref

        #: Reference span
        self.bref = bref

        #: Reference wing area
        self.sref = sref

        #: Symmetry key for the aero coordinate x-z plane. See Remark 6.
        #: (Integer = +1 for symmetry, 0 for no symmetry, and -1 for antisymmetry;
        #: Default = 0)
        self.sym_xz = sym_xz

        self.xyz_ref = xyz_ref

        if self.acsid is None:
            self.acsid = 0
        if self.rcsid is None:
            self.rcsid = 0
        if self.sym_xz is None:
            self.sym_xz = 0
        if self.sym_xy is None:
            self.sym_xy = 0
        self.acsid_ref = None
        self.rcsid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEROZ card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card


        $       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        $+ABC   REFX  REFY  REFZ
        AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                59.53 0.0   0.0

        """
        acsid = integer_or_blank(card, 1, 'acsid', 0)
        sym_xz = string(card, 2, 'sym_xz')
        flip = string(card, 3, 'flip')
        fm_mass_unit = string(card, 4, 'fm_mass_unit')
        fm_length_unit = string(card, 5, 'fm_length_unit')

        # YES-aero=half,structure=half
        # NO-aero=full; structure=full
        # H2F-aero=full; structure=half
        assert sym_xz in ['YES', 'NO', 'H2F'], 'sym_xz=%r' % flip


        # YES-structure=left,aero=right
        assert flip in ['YES', 'NO'], 'flip=%r' % flip
        assert fm_mass_unit in ['SLIN', 'LBM'], 'fm_mass_unit=%r' % fm_mass_unit
        assert fm_length_unit in ['IN'], 'fm_length_unit=%r' % fm_length_unit

        #rcsid = integer_or_blank(card, 2, 'rcsid', 0)

        cref = double_or_blank(card, 6, 'cRef', 1.)
        bref = double_or_blank(card, 7, 'bRef', 1.)
        sref = double_or_blank(card, 8, 'Sref', 1.)

        xref = double_or_blank(card, 9, 'xRef', 0.)
        yref = double_or_blank(card, 10, 'yRef', 0.)
        zref = double_or_blank(card, 11, 'zref', 0.)
        xyz_ref = [xref, yref, zref]

        assert len(card) <= 12, 'len(AEROZ card) = %i\ncard=%s' % (len(card), card)

        # faking data to not change gui
        rcsid = 0
        #sym_xy = 0
        return AEROZ(fm_mass_unit, fm_length_unit,
                     cref, bref, sref, acsid=acsid, rcsid=rcsid,
                     sym_xz=sym_xz, flip=flip, xyz_ref=xyz_ref,
                     comment=comment)

    def Acsid(self):
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    def Rcsid(self):
        try:
            return self.rcsid_ref.cid
        except AttributeError:
            return self.rcsid

    #def validate(self):
        #msg = ''
        #if not isinstance(self.acsid, integer_types):
            #msg += 'acsid=%s must be an integer; type=%s\n' % (self.acsid, type(self.acsid))
        #if not isinstance(self.rcsid, integer_types):
            #msg += 'rcsid=%s must be an integer; type=%s\n' % (self.rcsid, type(self.rcsid))
        #if not isinstance(self.cref, float):
            #msg += 'cref=%s must be an float; type=%s\n' % (self.cref, type(self.cref))
        #if not isinstance(self.bref, float):
            #msg += 'bref=%s must be an float; type=%s\n' % (self.bref, type(self.bref))
        #if not isinstance(self.sref, float):
            #msg += 'sref=%s must be an float; type=%s\n' % (self.sref, type(self.sref))
        #if not isinstance(self.sym_xz, integer_types):
            #msg += 'sym_xz=%s must be an integer; type=%s\n' % (self.sym_xz, type(self.sym_xz))
        #if not isinstance(self.sym_xy, integer_types):
            #msg += 'sym_xy=%s must be an integer; type=%s\n' % (self.sym_xy, type(self.sym_xy))
        #if msg:
            #raise TypeError('There are errors on the AEROS card:\n%s%s' % (msg, self))

    def cross_reference(self, model: BDF) -> None:
        """
        Cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROZ'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)
        self.rcsid_ref = model.Coord(self.rcsid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Safe cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROZ'
        self.acsid_ref = model.safe_coord(self.acsid, None, xref_errors, msg=msg)
        self.rcsid_ref = model.safe_coord(self.rcsid, None, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.acsid_ref = None
        self.rcsid_ref = None

    def convert_to_zona(self, unused_model):
        #$       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        #$+ABC   REFX  REFY  REFZ
        #AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                #59.53 0.0   0.0
        cref = self.cref
        bref = self.bref
        sref = self.sref
        acsid = self.acsid
        rho_ref = 1.0
        if self.sym_xz == 'NO':
            sym_xz = 0
        elif self.sym_xz == 'YES':
            sym_xz = 1
        else:
            raise NotImplementedError(self.sym_xz)
        assert sym_xz in [0, 1], sym_xz
        aeros = AEROS(cref, bref, sref, acsid=acsid, rcsid=0, sym_xz=sym_xz, sym_xy=0,
                      comment=str(self))

        velocity = 1.
        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=sym_xz, sym_xy=0,
                    comment='')
        return aeros, aero

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        raise NotImplementedError()
        #list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       #self.bref, self.sref, self.sym_xz, self.sym_xy]
        #return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card

        """
        unused_sym_xz = set_blank_if_default(self.sym_xz, 0)
        unused_sym_xy = set_blank_if_default(self.sym_xy, 0)
        #$       ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS
        #$+ABC   REFX  REFY  REFZ
        #AEROZ   0     YES   NO   SLIN    IN       22.73 59.394 1175.8
                #59.53 0.0   0.0

        list_fields = ['AEROZ', self.Acsid(), self.sym_xz, self.flip,
                       self.fm_mass_unit, self.fm_length_unit,
                       self.cref, self.bref, self.sref] + list(self.xyz_ref)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MKAEROZ(BaseCard):
    type = 'MKAEROZ'
    def __init__(self, sid, mach, flt_id, filename, print_flag, freqs,
                 method=0, save=None, comment=''):
        """
        Parameters
        ==========
        sid : int
            the MKAEROZ id
        mach : float
            the mach number for the TRIM solution
        save : str
            save the AIC data to the filename
            SAVE    save the AICs
            ACQUIRE load an AIC database
            ADD     append the new acids to the existing AIC database
            RESTART continue an analysis
        filename : str
            the length of the file must be at most 56 characters
        print_flag : int
            ???
        freqs : List[float]
            ???
        method : int
            ???
        save : ???
            ???
        comment : str; default=''
             a comment for the card
        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        self.sid = sid
        self.mach = mach
        self.method = method
        self.flt_id = flt_id
        self.save = save
        self.freqs = freqs
        self.filename = filename
        self.print_flag = print_flag

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'IDMK')
        mach = double(card, 2, 'MACH')
        method = integer(card, 3, 'METHOD')
        flt_id = integer(card, 4, 'IDFLT')
        save = string_or_blank(card, 5, 'SAVE')
        filename_a = loose_string_or_blank(card, 6, 'FILENAMEA', '')
        filename_b = loose_string_or_blank(card, 7, 'FILENAMEB', '')
        #print(filename_a, filename_b)
        filename = (filename_a + filename_b).rstrip()
        print_flag = integer_or_blank(card, 8, 'PRINT_FLAG', 0)
        freqs = []
        ifreq = 1
        for ifield in range(9, len(card)):
            freq = double(card, ifield, 'FREQ%i'%  ifreq)
            freqs.append(freq)
            ifreq += 1
        return MKAEROZ(sid, mach, flt_id, filename, print_flag, freqs,
                       method=method, save=save, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card

        """
        filename_a = self.filename[:8]
        filename_b = self.filename[8:]
        list_fields = ['MKAEROZ', self.sid, self.mach, self.method, self.flt_id,
                       self.save, filename_a, filename_b, self.print_flag] + self.freqs
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class PANLST1(Spline):
    """
    Defines a set of aerodynamic boxes by the LABEL entry in CAERO7 or BODY7
    bulk data cards.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    +---------+-------+---------+------+------+------+----+-----+-------+
    |    1    |   2   |    3    |   4  |  5   |   6  |  7 |  8  |   9   |
    +=========+=======+=========+======+======+======+====+=====+=======+
    | PANLST1 | SETID | MACROID | BOX1 | BOX2 |      |    |     |       |
    +---------+-------+---------+------+------+------+----+-----+-------+
    | PANLST1 |  100  |   111   |  111 |  118 |      |    |     |       |
    +---------+-------+---------+------+------+------+----+-----+-------+

    PANLST1 is referred to by SPLINEi, ATTACH, LOADMOD, CPFACT, JETFRC, and/or
    AESURFZ bulk data card.
    """
    type = 'PANLST1'

    def __init__(self, eid, macro_id, box1, box2, comment=''):
        """
        Creates a PANLST1 card

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.macro_id = macro_id # points to CAERO7 / BODY7
        self.box1 = box1
        self.box2 = box2
        self.aero_element_ids = []
        self.caero_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PANLST3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        macro_id = integer(card, 2, 'macro_id')
        box1 = integer(card, 3, 'box1')
        box2 = integer(card, 4, 'box2')
        assert len(card) == 5, 'len(PANLST1 card) = %i\ncard=%s' % (len(card), card)
        return PANLST1(eid, macro_id, box1, box2, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PANLST1 eid=%s' % self.eid
        self.caero_ref = model.CAero(self.macro_id, msg=msg)
        self.aero_element_ids = np.arange(self.box1, self.box2)

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def raw_fields(self):
        list_fields = ['PANLST1', self.eid, self.macro_id, self.box1, self.box2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PANLST3(Spline):
    """
    Defines a set of aerodynamic boxes by the LABEL entry in CAERO7 or BODY7
    bulk data cards.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    +---------+-------+--------+--------+--------+-----+----+-----+-------+
    |    1    |   2   |   3    |    4   |    5   |  6  |  7 |  8  |   9   |
    +=========+=======+========+========+========+=====+====+=====+=======+
    | PANLST3 | SETID | LABEL1 | LABEL2 | LABEL3 | etc |    |     |       |
    +---------+-------+--------+--------+--------+-----+----+-----+-------+
    | PANLST3 | 100   | WING   | HTAIL  |        |     |    |     |       |
    +---------+-------+--------+--------+--------+-----+----+-----+-------+

    PANLST3 is referred to by SPLINEi, ATTACH, LOADMOD, CPFACT, JETFRC, and/or
    AESURFZ bulk data card.

    """
    type = 'PANLST3'

    def __init__(self, eid, panel_groups, comment=''):
        """
        Creates a PANLST3 card

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.panel_groups = panel_groups # points to CAERO7 / BODY7
        self.aero_element_ids = []
        self.caero_refs = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PANLST3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        group_id = 1
        panel_groups = []
        for ifield in range(2, len(card)):
            name = string(card, ifield, 'group_%i'%  (group_id))
            panel_groups.append(name)
        assert len(card) > 2, 'len(PANLST3 card) = %i; no panel_groups were defined\ncard=%s' % (len(card), card)
        return PANLST3(eid, panel_groups, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PANLST3 eid=%s' % self.eid
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        caero_refs = []
        aero_element_ids = []
        for caero_label in self.panel_groups:
            caero_eid = model.zona.caero_to_name_map[caero_label]
            caero_ref = model.CAero(caero_eid, msg=msg)
            caero_refs.append(caero_ref)
            eid = caero_ref.eid
            npanels = caero_ref.npanels
            if npanels == 0:
                model.log.warning('skipping PANLST3 because there are 0 panels in:\n%r' % caero_ref)
                continue
            aero_element_ids2 = range(eid, eid + npanels)
            assert len(aero_element_ids2) == npanels, npanels
            aero_element_ids += aero_element_ids2
        self.caero_refs = caero_refs
        self.aero_element_ids = aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def raw_fields(self):
        list_fields = ['PANLST3', self.eid] + self.panel_groups
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class PAFOIL7(BaseCard):
    """
    Defines an aerodynamic body macroelement of a body-like component.
    Similar to Nastran's CAERO2.

    +---------+----+------+------+-------+------+------+-------+------+
    |   1     |  2 |   3  |  4   |   5   |   6  |   7  |   8   |  9   |
    +=========+====+======+======+=======+======+======+=======+======+
    | PAFOIL7 | ID | ITAX | ITHR | ICAMR | RADR | ITHT | ICAMT | RADT |
    +---------+----+------+------+-------+------+------+-------+------+
    | PAFOIL7 |  1 | -201 |  202 |  203  |  0.1 |  211 |  212  |  0.1 |
    +---------+----+------+------+-------+------+------+-------+------+

    """
    type = 'PAFOIL7'

    def __init__(self, pid, i_axial,
                 i_thickness_root, i_camber_root, le_radius_root,
                 i_thickness_tip, i_camber_tip, le_radius_tip,
                 comment=''):
        """
        Defines a BODY7 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        pid : int
            PAFOIL7 identification number.
        i_axial : str
            Identification number of an AEFACT bulk data card used to
            specify the xcoordinate locations, in percentage of the
            chord length, where the thickness and camber are specified.
            ITAX can be a negative number (where ABS (ITAX) = AEFACT
            bulk data card identification number) to request linear
            interpolation.
        i_thickness_root / i_thickness_tip : int
            Identification number of an AEFACT bulk data card used to
            specify the half thickness of the airfoil at the wing
            root/tip.
        i_camber : int; default=0
            Identification number of an AEFACT bulk data card used to
            specify the camber of the airfoil at the wing root.
        le_radius_root / le_radius_root: float
            Leading edge radius at the root/tip normalized by the
            root/tip chord.
        i_thickness_tip : int
            Identification number of an AEFACT bulk data card used to
            specify the half thickness at the wing tip.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        self.pid = pid
        self.i_axial = i_axial

        self.i_thickness_root = i_thickness_root
        self.i_camber_root = i_camber_root
        self.le_radius_root = le_radius_root

        self.i_camber_tip = i_camber_tip
        self.le_radius_tip = le_radius_tip
        self.i_thickness_tip = i_thickness_tip

        self.i_thickness_root_ref = None
        self.i_camber_root_ref = None
        self.i_thickness_tip_ref = None
        self.i_camber_tip_ref = None
        self.i_axial_ref = None

    #@property
    #def cp(self):
        #return self.acoord
    #@property
    #def cp_ref(self):
        #return self.acoord_ref

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAFOIL7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """

        pid = integer(card, 1, 'pid')
        i_axial = integer(card, 2, 'i_axial')

        i_thickness_root = integer(card, 3, 'i_thickness_root')
        i_camber_root = integer(card, 4, 'i_camber_root')
        le_radius_root = double_or_blank(card, 5, 'le_radius_root')

        i_thickness_tip = integer(card, 6, 'i_thickness_tip')
        i_camber_tip = integer(card, 7, 'i_camber_tip')
        le_radius_tip = double_or_blank(card, 8, 'le_radius_tip')

        assert len(card) <= 9, 'len(PAFOIL7 card) = %i\ncard=%s' % (len(card), card)
        return PAFOIL7(pid, i_axial,
                       i_thickness_root, i_camber_root, le_radius_root,
                       i_thickness_tip, i_camber_tip, le_radius_tip,
                       comment=comment)

    #def ACoord(self):
        #if self.acoord_ref is not None:
            #return self.acoord_ref.cid
        #return self.acoord

    #def Pid(self):
        #if self.pid_ref is not None:
            #return self.pid_ref.pid
        #return self.pid

    #def Lsb(self):  # AEFACT
        #if self.lsb_ref is not None:
            #return self.lsb_ref.sid
        #return self.lsb

    #def Lint(self):  # AEFACT
        #if self.lint_ref is not None:
            #return self.lint_ref.sid
        #return self.lint

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PAFOIL7 pid=%s' % self.pid
        self.i_axial_ref = model.AEFact(abs(self.i_axial), msg=msg)

        self.i_thickness_root_ref = model.AEFact(self.i_thickness_root, msg=msg)
        self.i_camber_root_ref = model.AEFact(self.i_camber_root, msg=msg)

        self.i_thickness_tip_ref = model.AEFact(self.i_thickness_tip, msg=msg)
        self.i_camber_tip_ref = model.AEFact(self.i_camber_tip, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.i_thickness_root_ref = None
        self.i_camber_root_ref = None
        self.i_thickness_tip_ref = None
        self.i_camber_tip_ref = None

    def convert_to_nastran(self, model):
        """
        Should this be converted to a DMIG?

        +---------+----+------+------+-------+------+------+-------+------+
        |   1     |  2 |   3  |  4   |   5   |   6  |   7  |   8   |  9   |
        +=========+====+======+======+=======+======+======+=======+======+
        | PAFOIL7 | ID | ITAX | ITHR | ICAMR | RADR | ITHT | ICAMT | RADT |
        +---------+----+------+------+-------+------+------+-------+------+
        | PAFOIL7 |  1 | -201 |  202 |  203  |  0.1 |  211 |  212  |  0.1 |
        +---------+----+------+------+-------+------+------+-------+------+

        """
        raise NotImplementedError('PAFOIL7: convert_to_nastran')



    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        #pid = integer(card, 1, 'pid')
        #i_axial = integer(card, 2, 'i_axial')

        #i_thickness_root = integer(card, 3, 'i_thickness_root')
        #i_camber_root = integer(card, 4, 'i_camber_root')
        #le_radius_root = double_or_blank(card, 5, 'le_radius_root')

        #i_thickness_tip = integer(card, 6, 'i_thickness_tip')
        #le_radius_tip = integer(card, 7, 'le_radius_tip')
        #i_camber_tip = double_or_blank(card, 8, 'i_camber_tip')

        list_fields = [
            'PAFOIL7', self.pid, self.i_axial,
            self.i_thickness_root, self.i_camber_root, self.le_radius_root,
            self.i_thickness_tip, self.i_camber_tip, self.le_radius_tip,
        ]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class BODY7(BaseCard):
    """
    Defines an aerodynamic body macroelement of a body-like component.
    Similar to Nastran's CAERO2.

    +--------+-----+-----+----+-----+------+-----+------+------+
    |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
    +========+=====+=====+====+=====+======+=====+======+======+
    | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
    +--------+-----+-----+----+-----+------+-----+------+------+
    |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
    +--------+-----+-----+----+-----+------+-----+------+------+

    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    |   1   |    2    |   3   |     4   |    5   |   6  |     7   |    8    |    9    |
    +=======+=========+=======+=========+========+======+=========+=========+=========+
    | BODY7 |   BID   | LABEL | IPBODY7 | ACOORD | NSEG | IDMESH1 | IDMESH2 | IDMESH3 |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    |       | IDMESH4 |  etc  |         |        |      |         |         |         |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    | BODY7 |    4    |  BODY |    2    |    8   |   4  |    20   |    21   |    22   |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    |       |    23   |       |         |        |      |         |         |         |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+

    """
    type = 'BODY7'

    def __init__(self, eid, label, pid, nseg, idmeshes, acoord=0, comment=''):
        """
        Defines a BODY7 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            body id
        label : str
            An arbitrary character string used to define the body.
        pid  : int; default=0
            Identification number of PBODY7 bulk data card
            (specifying body wake and/or inlet aerodynamic boxes)
        acoord : int; default=0
            Identification number of ACOORD bulk data card
            (specifying body center line location and orientation)
        nseg : int
            Number of body segments
        idmeshes : List[int]
            Identification number of SEGMESH bulk data card (specifying body segments).
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid
        self.label = label

        #: Property identification number of a PAERO7 entry.
        self.pid = pid
        self.nseg = nseg
        self.idmeshes = idmeshes
        self.acoord = acoord
        self.pid_ref = None
        self.acoord_ref = None
        self.ascid_ref = None
        self.segmesh_refs = None

    #@property
    #def cp(self):
        #return self.acoord
    #@property
    #def cp_ref(self):
        #return self.acoord_ref

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a BODY7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        label = string(card, 2, 'label')
        assert len(card) >= 3, 'len(BODY7 card) = %i\ncard=%s' % (len(card), card)
        pid = integer_or_blank(card, 3, 'pid')
        acoord = integer_or_blank(card, 4, 'acoord', 0)
        nseg = integer_or_blank(card, 5, 'nseg')

        idmeshes = []
        for i, ifield in enumerate(range(6, len(card))):
            segmesh = integer(card, ifield, 'idmesh_%i' % (i+1))
            idmeshes.append(segmesh)
        assert len(card) <= 13, 'len(BODY7 card) = %i\ncard=%s' % (len(card), card)
        return BODY7(eid, label, pid, nseg, idmeshes, acoord=acoord, comment=comment)

    def ACoord(self):
        if self.acoord_ref is not None:
            return self.acoord_ref.cid
        return self.acoord

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    @property
    def nboxes(self):
        if self.nsb > 0:
            return self.nsb
        return len(self.lsb_ref.fractions) # AEFACT

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by BODY7 eid=%s' % self.eid
        self.segmesh_refs = []
        for segmesh_id in self.idmeshes:
            segmesh_ref = model.PAero(segmesh_id, msg=msg)  # links to SEGMESH/PAERO7
            self.segmesh_refs.append(segmesh_ref)

        #if self.pid is not None:
            #self.pid_ref = model.PAero(self.pid, msg=msg)  # links to PAERO7
        #self.acoord_ref = model.Coord(self.acoord, msg=msg)
        #if self.nsb == 0:
            #self.lsb_ref = model.AEFact(self.lsb, msg=msg)
        #if self.nint == 0:
            #self.lint_ref = model.AEFact(self.lint, msg=msg)
        if self.acoord is not None:
            self.acoord_ref = model.Coord(self.acoord, msg=msg)
        #self.ascid_ref = model.Acsid(msg=msg)
        self.ascid_ref = model.Coord(0, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.acoord = self.ACoord()
        self.pid_ref = None
        self.acoord_ref = None

    def convert_to_nastran(self, model):
        """
        +--------+-----+-----+----+-----+------+-----+------+------+
        |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
        +========+=====+=====+====+=====+======+=====+======+======+
        | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
        +--------+-----+-----+----+-----+------+-----+------+------+
        |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
        +--------+-----+-----+----+-----+------+-----+------+------+

        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        |   1   |    2    |   3   |     4   |    5   |   6  |     7   |    8    |    9    |
        +=======+=========+=======+=========+========+======+=========+=========+=========+
        | BODY7 |   BID   | LABEL | IPBODY7 | ACOORD | NSEG | IDMESH1 | IDMESH2 | IDMESH3 |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        |       | IDMESH4 |  etc  |         |        |      |         |         |         |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        | BODY7 |    4    |  BODY |    2    |    8   |   4  |    20   |    21   |    22   |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        |       |    23   |       |         |        |      |         |         |         |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+

        """
        pid = max(model.paeros) + 1000
        igroup = 1
        orient = 'ZY'
        cp = 0
        #width = 1

        #nsb  : AEFACT id for defining the location of the slender body elements
        #lsb  : AEFACT id for defining the location of interference elements
        #nint : Number of slender body elements
        #lint : Number of interference elements
        aefact_id = len(model.aefacts) + 1
        xs_id = aefact_id
        half_width_id = aefact_id + 1
        theta1_id = aefact_id + 2
        theta2_id = aefact_id + 3

        lsb = xs_id
        lint = xs_id

        #+---------+--------+-------+------+---------+-------+------+------+-----+
        #|    1    |    2   |   3   |   4  |    5    |   6   |   7  |   8  |  9  |
        #+=========+========+=======+======+=========+=======+======+======+=====+
        #| SEGMESH | IDMESH | NAXIS | NRAD | NOSERAD | IAXIS |      |      |     |
        #|         | ITYPE1 |   X1  | CAM1 |   YR1   |  ZR1  | IDY1 | IDZ1 |     |
        #|         | ITYPE2 |   X2  | CAM2 |   YR2   |  ZR2  | IDY2 | IDZ2 |     |
        #|         | ITYPE3 |   X3  | CAM3 |   YR3   |  ZR3  | IDY3 | IDZ3 |     |
        #+---------+--------+-------+------+---------+-------+------+------+-----+
        xpoints = []
        half_widths = []

        try:
            origin_x, origin_y, origin_z = self.acoord_ref.origin
        except AttributeError:  # pragma: no cover
            print(self.get_stats())
            raise
        #x_offset = origin_x + x
        #y_offset = origin_y + y
        #z_offset = origin_z + z
        nsegmesh = len(self.segmesh_refs)
        if nsegmesh == 0:
            raise RuntimeError('Number of SEGMESH  references on BODY7=0\n%s' % str(self))
        for isegmesh, segmesh in enumerate(self.segmesh_refs):
            itypes = segmesh.itypes

            #xs = segmesh.xs
            #idys_ref = segmesh.idys_ref
            #idzs_ref = segmesh.idzs_ref
            nitypes = len(itypes)
            idys_ref = [None] * nitypes if segmesh.idys_ref is None else segmesh.idys_ref
            idzs_ref = [None] * nitypes if segmesh.idzs_ref is None else segmesh.idzs_ref

            cambers = segmesh.cambers
            yrads = segmesh.ys
            zrads = segmesh.zs

            # what????
            if isegmesh in [0, nsegmesh - 1]:
                xs2 = segmesh.xs
                idys_ref2 = idys_ref
                idzs_ref2 = idzs_ref
            else:
                xs2 = segmesh.xs[1:]
                idys_ref2 = idys_ref[1:]
                idzs_ref2 = idzs_ref[1:]

            xpoints += xs2
            yz_mean = []


            thetas = self._get_thetas()
            for itype, camber, yrad, zrad, idy_ref, idz_ref in zip(
                    itypes, cambers, yrads, zrads, idys_ref2, idzs_ref2):
                if itype == 1:
                    # Body of Revolution
                    # Xi, CAMi, YRi
                    radius = yrad
                    aspect_ratio = 1.
                    yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
                    ypoints = yz[:, 0]
                    zpoints = camber + yz[:, 1]
                elif itype == 2:
                    # Elliptical body
                    height = zrad
                    width = yrad
                    aspect_ratio = height / width
                    radius = height
                    yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
                    ypoints = yz[:, 0]
                    zpoints = yz[:, 1]

                elif itype == 3:
                    # Arbitrary body using AEFACTss
                    try:
                        ypoints = idy_ref.fractions
                        zpoints = idz_ref.fractions
                    except AttributeError:  # pragma: no cover
                        print('idy_ref = %s' % idy_ref)
                        print('idz_ref = %s' % idz_ref)
                        print(self.get_stats())
                        raise
                else:  # pramga: no cover
                    msg = 'Unsupported itype=%s (must be 1/2/3)\n%s' % (itype, str(self))
                    raise NotImplementedError(msg)

                width = ypoints.max() - ypoints.min()
                height = zpoints.max() - zpoints.min()
                #elliptical_area = pi * width * height
                average_radius = (width + height) / 4.
                half_widths.append(average_radius)

                ymeani = ypoints.mean()
                zmeani = zpoints.mean()
                yz_mean.append([ymeani, zmeani])

            # I think you could area weight this and get a better mean...
            ymean, zmean = np.mean(yz_mean, axis=0)

        xpoints_local = [xi for xi in xpoints]
        assert len(half_widths) == len(xpoints_local)

        half_width = max(half_widths)
        AR = 1.0

        p1 = [origin_x, origin_y + ymean, origin_z + zmean]
        x12 = max(xpoints) - min(xpoints)
        dash = '-' * 80 + '\n'
        comment = dash
        comment += self.comment
        caero2 = CAERO2(self.eid, pid, igroup, p1, x12,
                        cp=cp,
                        nsb=0, lsb=lsb,
                        nint=0, lint=lint, comment=comment)

        #
        lrsb = half_width_id # slender body
        lrib = half_width_id # interference


        # theta arrays (AEFACTs)
        #lth1 = theta1_id
        #lth2 = theta2_id

        # 0-57 is excluded
        angles_body = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0,
                       200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0]
        angles_fin = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0,
                      200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0]


        aefact_xs = AEFACT(xs_id, xpoints_local, comment=dash+'Xs')
        aefact_width = AEFACT(half_width_id, half_widths, comment='half_widths')
        aefact_theta1 = AEFACT(theta1_id, angles_body, comment='angles_body')
        aefact_theta2 = AEFACT(theta2_id, angles_fin, comment='angles_fin')

        # which segments use theta1 array
        lth = [1, 10] #nsegments] # t
        thi = [1]
        thn = [1]
        paero2 = PAERO2(pid, orient, half_width, AR, thi, thn,
                        lrsb=lrsb, lrib=lrib,
                        lth=lth, comment='')
        caero2.validate()
        paero2.validate()
        return caero2, paero2, aefact_xs, aefact_width, aefact_theta1, aefact_theta2

    def _get_nthetas(self):
        return self.segmesh_refs[0].nradial  # npoints
        #nthetas = 17
        #for itype, idy_ref, unused_idz_ref in zip(itypes, idys_ref2, idzs_ref2):
            #if itype == 3:
                #fractions = idy_ref.fractions
                #nthetas = len(fractions)
                #break
        #return nthetas

    def _get_thetas(self):
        nthetas = self._get_nthetas()
        thetas = np.radians(np.linspace(0., 360., nthetas))
        return thetas

    def get_points(self):
        """creates a 1D representation of the BODY7"""
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))

        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return [p1, p2]

    @property
    def npanels(self):
        nz = len(self.segmesh_refs)
        unused_segmesh = self.segmesh_refs[0]
        nthetas = self._get_nthetas()
        npanels = nz * (nthetas - 1)
        return npanels

    def get_points_elements_3d(self):
        """
        Gets the points/elements in 3d space as CQUAD4s
        The idea is that this is used by the GUI to display CAERO panels.

        TODO: doesn't support the aero coordinate system

        """
        #paero2 = self.pid_ref
        xyz = []
        element = []
        npoints = 0
        for segmesh in self.segmesh_refs:
            #print(segmesh)
            xyzi, elementi = self._get_points_elements_3di(segmesh)
            xyz.append(xyzi)
            element.append(elementi + npoints)
            npoints += xyzi.shape[0]

        xyzs = np.vstack(xyz)
        elements = np.vstack(element)
        assert xyzs is not None, str(self)
        assert elements is not None, str(self)

        return xyzs, elements

    def _get_points_elements_3di(self, segmesh):
        """
        points (nchord, nspan) float ndarray; might be backwards???
            the points
        elements (nquads, 4) int ndarray
            series of quad elements
            nquads = (nchord-1) * (nspan-1)
        """
        #lengths_y = []
        #lengths_z = []
        nx = segmesh.naxial
        ny = segmesh.nradial

        xs = []
        ys = []
        zs = []

        origin_x, origin_y, origin_z = self.acoord_ref.origin

        nthetas = segmesh.nradial
        thetas = np.radians(np.linspace(0., 360., nthetas))
        for itype, x, yrad, zrad, camber, idy_ref, idz_ref in zip(
                segmesh.itypes,
                segmesh.xs, segmesh.ys, segmesh.zs,
                segmesh.cambers,
                segmesh.idys_ref, segmesh.idzs_ref):
            y = 0.
            z = 0.
            if itype == 1:
                # Body of Revolution
                # Xi, CAMi, YRi
                ## TODO: doesn't consider camber
                radius = yrad
                aspect_ratio = 1.
                yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
                ypoints = yz[:, 0]
                zpoints = camber + yz[:, 1]
            elif itype == 2:
                # Elliptical body
                #  Xi, YRi, ZRi
                height = zrad
                width = yrad
                aspect_ratio = height / width
                radius = height
                yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
                ypoints = yz[:, 0]
                zpoints = yz[:, 1]

            elif itype == 3:
                # Arbitrary body using AEFACTss
                #  Xi, IDYi, IDZi
                ypoints = idy_ref.fractions
                zpoints = idz_ref.fractions
                y = yrad
                z = zrad
            else:  # pramga: no cover
                msg = 'Unsupported itype=%s (must be 1/2/3)\n%s' % (itype, str(self))
                raise NotImplementedError(msg)

            assert len(ypoints) == len(zpoints), 'len(ypoints)=%s len(zpoints)=%s' % (len(ypoints), len(zpoints))
            nnodes = len(ypoints)

            x_offset = origin_x + x
            y_offset = origin_y + y
            z_offset = origin_z + z
            xs.append([x_offset] * nnodes)
            ys.append(y_offset + ypoints)
            zs.append(z_offset + zpoints)

        xyz = np.vstack([
            np.hstack(xs),
            np.hstack(ys),
            np.hstack(zs),
        ]).T
        elements = elements_from_quad(nx, ny, dtype='int32')  # nx,ny are points
        return xyz, elements

    #def set_points(self, points):
        #self.p1 = np.asarray(points[0])
        #p2 = np.asarray(points[1])
        #x12 = p2 - self.p1
        #self.x12 = x12[0]

    #def shift(self, dxyz):
        #"""shifts the aero panel"""
        #self.p1 += dxyz

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = ['BODY7', self.eid, self.label, self.Pid(), self.ACoord(),
                       self.nseg] + self.idmeshes
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SEGMESH(BaseCard):
    """
    Defines a grid system for a body segment; referenced by the BODY7 bulk data card.

    +---------+--------+-------+------+---------+-------+------+------+-----+
    |    1    |    2   |   3   |   4  |    5    |   6   |   7  |   8  |  9  |
    +=========+========+=======+======+=========+=======+======+======+=====+
    | SEGMESH | IDMESH | NAXIS | NRAD | NOSERAD | IAXIS |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE1 |   X1  | CAM1 |   YR1   |  ZR1  | IDY1 | IDZ1 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE2 |   X2  | CAM2 |   YR2   |  ZR2  | IDY2 | IDZ2 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE3 |   X3  | CAM3 |   YR3   |  ZR3  | IDY3 | IDZ3 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    | SEGMESH |    2   |  3    |  6   |         |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    1   |  0.0  |  0.0 | 0.0     |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    1   |  1.0  |  0.0 | 0.5     |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    3   |  2.0  |      |         |       |  103 | 104  |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+

    """
    type = 'SEGMESH'

    @property
    def pid(self):
        return self.segmesh_id
    @pid.setter
    def pid(self, segmesh_id):
        self.segmesh_id = segmesh_id

    def __init__(self, segmesh_id, naxial, nradial, nose_radius, iaxis,
                 itypes, xs, cambers, ys, zs, idys, idzs, comment=''):
        """
        Defines a SEGMESH card, which defines a cross-section for a PBODY7.

        Parameters
        ----------
        segmesh_id : int
            Body segment mesh identification number.
        naxial : int
            Number of axial stations (i.e., divisions) of the segment. (min=2).
        nradial : int
            Number of circumferential points of the segment (min=3).
        nose_radius : float
            Nose radius of blunt body.
            NOSERAD is active only if ZONA7U (Hypersonic Aerodynamic Method)
            is used (the METHOD entry of the MKAEROZ Bulk Data equals 2 or –2).
            Furthermore, NOSERAD is used only if the SEGMESH bulk data card is
            the first segment defined in the BODY7 bulk data card.
        iaxis : int
            The index of the axial station where the blunt nose ends.
            IAXIS is active only if ZONA7U (Hypersonic Aerodynamic
            Method) is used.
        ITYPEi : int
            Type of input used to define the circumferential box cuts
            - 1 body of revolution
            - 2 elliptical body
            - 3 arbitrary body
        Xi : List[float]
            X-location of the axial station; Xi must be in ascending
            order. (i.e., Xi+1 > Xi)
        cambers : List[float]
            Body camber at the Xi axial station. (Real)
        YRi : List[float]
            Body cross-sectional radius if ITYPEi = 1 or the semi-axis length
            of the elliptical body parallel to the Y-axis if ITYPEi=2.
        ZRi : List[float]
            The semi-axis length of the elliptical body parallel to the Z-axis.
            Used only if ITYPEi=2. (Real)
        IDYi : int
            Identification number of AEFACT bulk data card that specifies
            NRAD number of the Y-coordinate locations of the circumferential
            points at the Xi axial station. Use only if ITYPEi=3.
        IDZi : int
            Identification number of AEFACT bulk data card that specifies
            NRAD number of the Z-coordinate locations of the circumferential
            points at the Xi axial station. Use only if ITYPEi=3.
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        BaseCard.__init__(self)
        self.segmesh_id = segmesh_id
        self.naxial = naxial
        self.nradial = nradial
        self.nose_radius = nose_radius
        self.iaxis = iaxis
        self.itypes = itypes
        self.cambers = cambers
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.idys = idys
        self.idzs = idzs

        self.idys_ref = None
        self.idzs_ref = None
        self.pid_ref = None

    def validate(self):
        for i, itype in enumerate(self.itypes):
            assert itype in [1, 2, 3], 'itypes[%i]=%s is invalid; itypes=%s' % (i, itype, self.itypes)

        xi_old = self.xs[0]
        for i, xi in enumerate(self.xs[1:]):
            if xi <= xi_old:
                raise RuntimeError('xs=%s must be in ascending order\nx%i=%s x%i=%s (old)\n%s' % (
                    self.xs, i+2, xi, i+1, xi_old, str(self)))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SEGMESH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        segmesh_id = integer(card, 1, 'segmesh_id')
        naxial = integer(card, 2, 'naxial')
        nradial = integer(card, 3, 'nradial')
        nose_radius = double_or_blank(card, 4, 'nose_radius')
        iaxis = integer_or_blank(card, 5, 'iaxis')

        itypes = []
        xs = []
        ys = []
        zs = []
        cambers = []
        idys = []
        idzs = []
        assert len(card) >= 9, 'len(SEGMESH card) = %i\ncard=%s' % (len(card), card)

        for counter, ifield in enumerate(range(9, len(card), 8)):
            itype = integer(card, ifield, 'itype%i' % (counter+1))
            x = double_or_blank(card, ifield+1, 'itype%i' % (counter+1), 0.)
            camber = double_or_blank(card, ifield+2, 'camber%i' % (counter+1), 0.)
            y = double_or_blank(card, ifield+3, 'y%i' % (counter+1), 0.)
            z = double_or_blank(card, ifield+4, 'z%i' % (counter+1), 0.)
            idy = integer_or_blank(card, ifield+5, 'idy%i' % (counter+1))
            idz = integer_or_blank(card, ifield+6, 'idz%i' % (counter+1))
            itypes.append(itype)
            xs.append(x)
            ys.append(y)
            zs.append(z)
            cambers.append(camber)
            idys.append(idy)
            idzs.append(idz)
        assert len(itypes) == naxial, 'naxial=%s nradial=%s len(itypes)=%s' % (naxial, nradial, len(itypes))
        return SEGMESH(segmesh_id, naxial, nradial, nose_radius, iaxis,
                       itypes, xs, cambers, ys, zs, idys, idzs, comment=comment)

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SEGMESH eid=%s' % self.pid
        idys_ref = []
        idzs_ref = []
        for idy in self.idys:
            idy_ref = None
            if idy is not None and isinstance(idy, integer_types):
                idy_ref = model.AEFact(idy, msg=msg)
                assert len(idy_ref.fractions) > 2, 'idy_ref=%s' % idy_ref
            idys_ref.append(idy_ref)

        for idz in self.idzs:
            idz_ref = None
            if idz is not None and isinstance(idz, integer_types):
                idz_ref = model.AEFact(idz, msg=msg)
                assert len(idz_ref.fractions) > 2, 'idz_ref=%s' % idz_ref
            idzs_ref.append(idz_ref)
        self.idys_ref = idys_ref
        self.idzs_ref = idzs_ref
        #print(self.idys_ref)

    def safe_cross_reference(self, model, xref_errors):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        #self.cp = self.Cp()
        #self.idys = idys_ref
        #self.idzs = idzs_ref

        self.pid_ref = None
        #self.cp_ref = None
        self.idys_ref = None
        self.idzs_ref = None

    #def set_points(self, points):
        #self.p1 = np.asarray(points[0])
        #p2 = np.asarray(points[1])
        #x12 = p2 - self.p1
        #self.x12 = x12[0]

    def shift(self, dxyz):
        """shifts the aero panel"""
        self.p1 += dxyz

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = [
            'SEGMESH', self.segmesh_id, self.naxial, self.nradial, self.nose_radius,
            self.iaxis, None, None, None]
        for itype, x, camber, y, z, idy, idz in zip(self.itypes, self.xs, self.cambers,
                                                    self.ys, self.zs, self.idys, self.idzs):
            list_fields += [itype, x, camber, y, z, idy, idz, None]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class CAERO7(BaseCard):
    """
    Totally wrong...

    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |   1    |  2  |   3   |   4    |   5   |   6    |    7   |   8    |    9    |
    +========+=====+=======+========+=======+========+========+========+=========+
    | CAERO7 | WID | LABEL | ACOORD | NSPAN | NCHORD |  LSPAN |  ZTAIC | PAFOIL7 |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |        | XRL |  YRL  |   ZRL  |  RCH  |  LRCHD | ATTCHR | ACORDR |         |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |        | XTL |  YTL  |   ZTL  |  TCH  |  LTCHD | ATTCHT | ACORDT |         |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int, PAERO1
        int : PAERO1 ID
        PAERO1 : PAERO1 object (xref)
    igid : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int, CORDx
        int : coordinate system
        CORDx : Coordinate object (xref)
    nspan : int
        int > 0 : N spanwise boxes distributed evenly
        int = 0 : use lchord
    nchord : int
        int > 0 : N chordwise boxes distributed evenly
        int = 0 : use lchord
    lspan : int, AEFACT
        int > 0 : AEFACT reference for non-uniform nspan
        int = 0 : use nspan
        AEFACT : AEFACT object  (xref)
    lchord : int, AEFACT
        int > 0 : AEFACT reference for non-uniform nchord
        int = 0 : use nchord
        AEFACT : AEFACT object  (xref)
    comment : str; default=''
         a comment for the card

    """
    type = 'CAERO7'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nspan', 5:'nchord',
        6:'lspan', 7:'lchord', 8:'igid', 12:'x12', 16:'x43',
    }

    def __init__(self, eid, label, p1, x12, p4, x43,
                 cp=0, nspan=0, nchord=0, lspan=0, p_airfoil=None, ztaic=None, comment=''):
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO1
            int : PAERO1 ID
            PAERO1 : PAERO1 object (xref)
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card

        """
        BaseCard.__init__(self)
        if cp is None:
            cp = 0
        if nspan is None:
            nspan = 0
        if nchord is None:
            nchord = 0
        p1 = np.asarray(p1)
        p4 = np.asarray(p4)

        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.label = label

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.nchord = nchord
        self.lspan = lspan
        self.p_airfoil = p_airfoil
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43
        self.ztaic = ztaic

        self.pid_ref = None
        self.cp_ref = None
        self.lchord_ref = None
        self.lspan_ref = None
        self.ascid_ref = None
        self.box_ids = None
        self.pafoil_ref = None
        #self._init_ids() #TODO: make this work here?

    def validate(self):
        msg = ''
        is_failed = False
        if self.nspan == 0 and self.lspan == 0:
            msg += 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True
        if self.nspan <= 0:
            msg += 'NSPAN must be greater than 0; nspan=%r\n' % (
                self.nspan)
            is_failed = True

        if self.nchord <= 0:
            msg += 'NCHORD must be greater than 0; nchord=%r\n' % (
                self.nchord)
            is_failed = True
        if is_failed:
            msg += str(self)
            raise ValueError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4
        assert self.nspan < 100, 'nspan=%s\n%s' % (self.nspan, str(self))
        assert self.nchord < 100, 'nchord=%s\n%s' % (self.nchord, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        CAERO7  100101  RIBW1           2       25
                998.904  39.821 230.687 1298.159        310001
                1121.821 61.134 236.560 1175.243
        CAERO7  100201  RIBW2           2       25
                1121.821 61.134 236.560 1175.243
                1244.258 84.704 243.625 1052.805

        """
        eid = integer(card, 1, 'eid')
        name = string(card, 2, 'name')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        nchord = integer_or_blank(card, 5, 'nchord', 0)
        lspan = integer_or_blank(card, 6, 'aefact_lchord', 0)
        if lspan:
            lspan = 0
        ztaic = integer_or_blank(card, 7, 'ztaic')
        p_airfoil = integer_or_blank(card, 8, 'aefact')
        #assert cp == 0
        #igroup = integer(card, 8, 'igid')

        x1 = double_or_blank(card, 9, 'x1', 0.0)
        y1 = double_or_blank(card, 10, 'y1', 0.0)
        z1 = double_or_blank(card, 11, 'z1', 0.0)
        p1 = np.array([x1, y1, z1])
        x12 = double_or_blank(card, 12, 'x12', 0.)
        unused_lchord_root = integer_or_blank(card, 13, 'lchord_root')
        unused_attach_root = integer_or_blank(card, 14, 'attach_root')
        unused_achord_root = integer_or_blank(card, 15, 'achord_root')

        x4 = double_or_blank(card, 17, 'x4', 0.0)
        y4 = double_or_blank(card, 18, 'y4', 0.0)
        z4 = double_or_blank(card, 19, 'z4', 0.0)
        p4 = np.array([x4, y4, z4])
        x43 = double_or_blank(card, 20, 'x43', 0.)
        unused_lchord_tip = integer_or_blank(card, 21, 'lchord_tip')
        unused_attach_tip = integer_or_blank(card, 22, 'attach_tip')
        unused_achord_tip = integer_or_blank(card, 23, 'achord_tip')

        assert len(card) <= 23, 'len(CAERO7 card) = %i\ncard=%s' % (len(card), card)
        return CAERO7(eid, name, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, nchord=nchord, lspan=lspan,
                      p_airfoil=p_airfoil, ztaic=ztaic,
                      comment=comment)

    def flip_normal(self):
        self.p1, self.p4 = self.p4, self.p1
        self.x12, self.x43 = self.x43, self.x12

    def _init_ids(self, dtype='int32'):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)

        """
        nchord, nspan = self.shape
        assert nchord >= 1, 'nchord=%s' % nchord
        assert nspan >= 1, 'nspan=%s' % nspan
        self.box_ids = np.zeros((nchord, nspan), dtype=dtype)

        npanels = nchord * nspan
        try:
            self.box_ids = np.arange(self.eid, self.eid + npanels,
                                     dtype=dtype).reshape(nspan, nchord).T
        except OverflowError:
            if dtype == 'int64':
                # we already tried int64
                msg = 'eid=%s nspan=%s nchord=%s' % (
                    self.eid, nspan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    #def Pid(self):
        #if self.pid_ref is not None:
            #return self.pid_ref.pid
        #return self.pid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO1 eid=%s' % self.eid
        #self.pid_ref = model.PAero(self.pid, msg=msg)
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.ascid_ref = model.Acsid(msg=msg)

        #if self.nchord == 0:
            #assert isinstance(self.lchord, integer_types), self.lchord
            #self.lchord_ref = model.AEFact(self.lchord, msg)
        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.AEFact(self.lspan, msg)

        if self.p_airfoil:
            self.pafoil_ref = model.zona.PAFOIL(self.p_airfoil, msg)
        self._init_ids()

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO1 eid=%s' % self.eid
        #try:
            #self.pid_ref = model.PAero(self.pid, msg=msg)
        #except KeyError:
            #pass

        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)
        self.ascid_ref = model.safe_acsid(msg=msg)

        #if self.nchord == 0:
            #assert isinstance(self.lchord, integer_types), self.lchord
            #self.lchord_ref = model.safe_aefact(self.lchord, self.eid, xref_errors, msg)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.safe_aefact(self.lspan, self.eid, xref_errors, msg)

        self._init_ids()

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        #self.pid = self.Pid()
        self.cp = self.Cp()
        self.lspan = self.get_LSpan()
        #self.pid_ref = None
        self.cp_ref = None
        self.lspan_ref = None
        self.ascid_ref = None

    def convert_to_nastran(self):
        """
        +--------+-----+-----+----+-------+--------+--------+--------+------+
        |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
        +========+=====+=====+====+=======+========+========+========+======+
        | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
        +--------+-----+-----+----+-------+--------+--------+--------+------+
        |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
        +--------+-----+-----+----+-------+--------+--------+--------+------+

        """
        pid = 1
        igroup = 1
        caero = CAERO1(self.eid, pid, igroup, self.p1, self.x12,
                       self.p4, self.x43, cp=self.cp,
                       nspan=self.nspan, lspan=self.lspan,
                       nchord=self.nchord, lchord=0,
                       comment=self.comment)
        caero.validate()
        return caero

    @property
    def min_max_eid(self):
        """
        Gets the min and max element ids of the CAERO card

        Returns
        -------
        min_max_eid : (2, ) list
            The [min_eid, max_eid]

        """
        nchord, nspan = self.shape
        return [self.eid, self.eid + nchord * nspan]

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        if self.cp_ref is None and self.cp == 0:
            p1 = self.p1
            p4 = self.p4
        else:
            p1 = self.cp_ref.transform_node_to_global(self.p1)
            p4 = self.cp_ref.transform_node_to_global(self.p4)

        if self.ascid_ref is None:
            # yes, this really does list + array addition
            p2 = p1 + np.array([self.x12, 0., 0.])
            p3 = p4 + np.array([self.x43, 0., 0.])
        else:
            p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))
            p3 = p4 + self.ascid_ref.transform_vector_to_global(np.array([self.x43, 0., 0.]))
        return [p1, p2, p3, p4]

    def get_box_index(self, box_id):
        """
        Get the index of ``self.box_ids`` that coresponds to the given box id.

        Parameters
        -----------
        box_id : int
            Box id to ge tthe index of.

        Returns
        --------
        index : tuple
            Index of ``self.box_ids`` that coresponds to the given box id.

        """
        if box_id not in self.box_ids:
            self._box_id_error(box_id)
        index = np.where(self.box_ids == box_id)
        index = (index[0][0], index[1][0])
        return index

    def get_box_quarter_chord_center(self, box_id):
        """
        The the location of the quarter chord of the box along the centerline.

        Parameters
        -----------
        box_id : int
            Box id.

        Returns
        --------
        xyz_quarter_chord : ndarray
            Location of box quater chord in global.

        """
        return self._get_box_x_chord_center(box_id, 0.25)

    def get_box_mid_chord_center(self, box_id):
        """
        The the location of the mid chord of the box along the centerline.

        Parameters
        -----------
        box_id : int
            Box id.

        Returns
        --------
        xyz_mid_chord : ndarray
            Location of box mid chord in global.

        """
        return self._get_box_x_chord_center(box_id, 0.5)

    def _get_box_x_chord_center(self, box_id, x_chord):
        """The the location of the x_chord of the box along the centerline."""
        raise NotImplementedError('CAERO7: _get_box_x_chord_center')
        #if self.lchord != 0 or self.lspan != 0:
            #raise NotImplementedError()
        #ichord, ispan = self.get_box_index(box_id)

        #le_vector = self.p4 - self.p1
        #delta_xyz = le_vector * ((ispan + 0.5)/self.nspan)
        #yz = delta_xyz[1:3] + self.p1[1:3]
        #chord = ((ispan + 0.5)/self.nspan) * (self.x43 - self.x12) + self.x12
        #x = (ichord + x_chord)/self.nchord * chord + self.p1[0] + delta_xyz[0]
        #return np.array([x, yz[0], yz[1]])

    def _box_id_error(self, box_id):
        """
        Raise box_id IndexError.
        """
        msg = '%i not in range of aero box ids\nRange: %i to %i' % (box_id, self.box_ids[0, 0],
                                                                    self.box_ids[-1, -1])
        raise IndexError(msg)

    @property
    def npanels(self):
        nchord, nspan = self.shape
        return nchord * nspan

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        if self.nchord == 0:
            x = self.lchord_ref.fractions
            nchord = len(x) - 1
        else:
            nchord = self.nchord

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nchord < 1 or nspan < 1:
            msg = 'CAERO7 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

        Returns
        -------
        npoints : int
            The number of nodes for the CAERO
        nelmements : int
            The number of elements for the CAERO

        """
        nchord, nspan = self.shape
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        if self.nchord == 0:
            x = self.lchord_ref.fractions
            nchord = len(x) - 1
        else:
            nchord = self.nchord
            x = np.linspace(0., 1., nchord + 1)

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = 'CAERO7 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

        Returns
        -------
        points : (nnodes,3) ndarray of floats
            the array of points
        elements : (nelements,4) ndarray of integers
            the array of point ids

        """
        p1, p2, p3, p4 = self.get_points()
        x, y = self.xy
        # We're reordering the points so we get the node ids and element ids
        # to be consistent with Nastran.  This is only useful if you're plotting
        # aero panel forces
        #
        # this gives us chordwise panels and chordwise nodes
        return points_elements_from_quad_points(p1, p4, p3, p2, y, x, dtype='int32')

    def set_points(self, points):
        self.p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        self.p4 = points[3]
        self.x12 = p2[0] - self.p1[0]
        self.x43 = p3[0] - self.p4[0]
        assert self.x12 >= 0., 'p1=%s p2=%s' % (self.p1, p2)
        assert self.x43 >= 0., 'p4=%s p3=%s' % (self.p4, p3)
        assert self.x12 > 0. or self.x43 > 0., 'points=%s' % (points)

    def shift(self, dxyz):
        """shifts the aero panel"""
        self.p1 += dxyz
        self.p4 += dxyz

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
          the fields that define the card

        """
        cp = self.Cp()
        nspan = self.nspan
        nchord = self.nchord
        lspan = self.get_LSpan()
        list_fields = (
            ['CAERO7', self.eid, self.label, cp, nspan, nchord, lspan, self.ztaic,
             self.p_airfoil,] +
            list(self.p1) + [self.x12, None, None, None, None] +
            list(self.p4) + [self.x43, None, None, None, None])
        return list_fields

    def get_LSpan(self):
        if self.lspan_ref is not None:
            return self.lspan_ref.sid
        return self.lspan

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : LIST
            The fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        nspan = set_blank_if_default(self.nspan, 0)
        nchord = set_blank_if_default(self.nchord, 0)
        #lchord = set_blank_if_default(self.get_LChord(), 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)
        #lchord = 0
        lspan = 0
        list_fields = (
            ['CAERO7', self.eid, self.label, cp, nspan, nchord, lspan, self.ztaic,
             self.p_airfoil,] +
            list(self.p1) + [self.x12, None, None, None, None] +
            list(self.p4) + [self.x43, None, None, None, None])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIM_ZONA(BaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    """
    type = 'TRIM_ZONA'
    _field_map = {
        1: 'sid', 2:'mach', 3:'q', 8:'aeqr',
    }

    def __init__(self, sid, mkaeroz, q, cg, true_g, nxyz, pqr, loadset,
                 labels, uxs, comment=''):
        """
        Creates a TRIM card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        q : float
            dynamic pressure
        true_g : float
            ???
        nxyz : List[float]
            ???
        pqr : List[float]
            [roll_rate, pitch_rate, yaw_rate]
        loadset : int
            Identification number of a SET1 or SETADD bulk data card that
            specifies a set of identification numbers of TRIMFNC or
            TRIMADD bulk data card.  All values of the trim functions
            defined by the TRIMFNC or TRIMADD bulk data card are computed
            and printed out.
        labels : List[str]
            names of the fixed variables
        uxs : List[float]
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

        self.cg = cg
        self.nxyz = nxyz
        self.true_g = true_g
        self.pqr = pqr
        self.loadset = loadset

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        self.labels = labels

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = uxs

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mkaeroz = integer(card, 2, 'mkaeroz')
        qinf = double(card, 3, 'dynamic_pressure')
        # 5
        # 6
        cg = [
            double(card, 7, 'cg-x'),
            double(card, 8, 'cg-y'),
            double(card, 9, 'cg-z'),
        ]

        unused_wtmass = double(card, 9, 'wtmass')
        unused_weight = double(card, 10, 'weight')
        unused_inertia = [
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

        p = double_or_string(card, 21, 'P')
        q = double_or_string(card, 22, 'Q')
        r = double_or_string(card, 23, 'R')
        pqr = [p, q, r]
        loadset = integer_or_blank(card, 24, 'loadset')

        labels = []
        uxs = []

        i = 25
        n = 1
        while i < len(card):
            label = integer(card, i, 'label%i' % n)
            ux = double_or_string(card, i + 1, 'ux%i' % n)
            if isinstance(ux, str):
                assert ux == 'FREE', 'ux=%r' % ux
            #print('  label=%s ux=%s' % (label, ux))
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        assert len(card) >= 25, 'len(TRIM card) = %i\ncard=%s' % (len(card), card)
        return TRIM_ZONA(sid, mkaeroz, qinf, cg, true_g, nxyz, pqr, loadset,
                         labels, uxs, comment=comment)

    def validate(self):
        assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

        assert isinstance(self.nxyz[0], float) or self.nxyz[0] in ['FREE', 'NONE'], 'nx=%r' % self.nxyz[0]
        assert isinstance(self.nxyz[1], float) or self.nxyz[1] in ['FREE', 'NONE'], 'ny=%r' % self.nxyz[1]
        assert isinstance(self.nxyz[2], float) or self.nxyz[2] in ['FREE', 'NONE'], 'nz=%r' % self.nxyz[2]

        assert isinstance(self.pqr[0], float) or self.pqr[0] in ['FREE', 'NONE'], 'p=%r' % self.pqr[0]
        assert isinstance(self.pqr[1], float) or self.pqr[1] in ['FREE', 'NONE'], 'q=%r' % self.pqr[1]
        assert isinstance(self.pqr[2], float) or self.pqr[2] in ['FREE', 'NONE'], 'r=%r' % self.pqr[2]

        assert self.q > 0.0, 'q=%s\n%s' % (self.q, str(self))
        if len(set(self.labels)) != len(self.labels):
            msg = 'not all labels are unique; labels=%s' % str(self.labels)
            raise RuntimeError(msg)
        if len(self.labels) != len(self.uxs):
            msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            raise RuntimeError(msg)

    def cross_reference(self, model: BDF) -> None:
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model):
        mkaeroz_id = self.mkaeroz
        mkaeroz = model.zona.mkaeroz[mkaeroz_id]
        #print(mkaeroz)
        mach = mkaeroz.mach
        labels = []
        uxs = []
        comment = str(self)
        for label_id, ux in zip(self.labels, self.uxs):
            if ux != 'FREE':
                trimvar = model.zona.trimvar[label_id]
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
        mach = 1.0
        aeqr = 0.0
        list_fields = ['TRIM', self.sid, mach, self.q]
        nlabels = len(self.labels)
        assert nlabels > 0, self.labels
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            list_fields += [label, ux]
            if i == 1:
                list_fields += [aeqr]
        if nlabels == 1:
            list_fields += [None, None, aeqr]
        return list_fields

    def repr_fields(self):
        # fixes a Nastran bug
        #aeqr = set_blank_if_default(self.aeqr, 1.0)
        aeqr = 0.

        list_fields = self.raw_fields()
        list_fields[8] = aeqr
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        return ''
        #card = self.repr_fields()
        #return self.comment + print_card_8(card)

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

    """
    type = 'TRIMLNK'
    def __init__(self, link_id, sym, coeffs, var_ids, comment=''):
        """
        Creates a TRIMLNK card

        Parameters
        ----------
        link_id : int
            the TRIMLNK id
        sym : ???
            ???
        coeffs : ???
            ???
        var_ids : ???
            ???
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

    @classmethod
    def add_card(cls, card, comment=''):
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
        assert len(card) >= 5, 'len(TRIMLNK card) = %i\ncard=%s' % (len(card), card)
        return TRIMLNK(link_id, sym, coeffs, var_ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model):
        label = 'LNK_%s' % self.link_id
        trimvars = model.zona.trimvar

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
    def __init__(self, var_id, label, lower, upper,
                 trimlnk, dmi, sym, initial,
                 dcd, dcy, dcl, dcr, dcm, dcn, comment=''):
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
        self.trimlnk = trimlnk
        self.dmi = dmi
        self.sym = sym
        self.initial = initial
        self.dcd = dcd
        self.dcy = dcy
        self.dcl = dcl
        self.dcr = dcr
        self.dcm = dcm
        self.dcn = dcn

    @classmethod
    def add_card(cls, card, comment=''):
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
        trimlnk = integer_or_blank(card, 5, 'TRIMLNK')
        dmi = blank(card, 6, 'DMI')
        sym = string_or_blank(card, 7, 'sym')
        initial = blank(card, 8, 'initial')
        dcd = double_or_blank(card, 9, 'DCD')
        dcy = double_or_blank(card, 10, 'DCY')
        dcl = double_or_blank(card, 11, 'DCL')
        dcr = double_or_blank(card, 12, 'DCR')
        dcm = double_or_blank(card, 13, 'DCM')
        dcn = double_or_blank(card, 14, 'DCN')
        return TRIMVAR(var_id, label, lower, upper, trimlnk, dmi, sym,
                       initial, dcd, dcy, dcl, dcr, dcm,
                       dcn, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model):
        raise NotImplementedError()
        #mkaeroz_id = self.mkaeroz
        #mkaeroz = model.zona.mkaeroz[mkaeroz_id]
        #mach = mkaeroz.mach
        #labels = []
        #uxs = []
        #for label_id, ux in zip(self.labels, self.uxs):
            #if ux != 'FREE':
                #label = model.zona.trimvar[label_id]
                #labels.append(label)
                #uxs.append(ux)
        #trim = TRIM(self.sid, mach, self.q, labels, uxs,
                    #aeqr=1.0, comment=str(self))
        #trim.validate()
        #return trim

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIMVAR', self.var_id, self.label, self.lower, self.upper,
                       self.trimlnk, self.dmi, self.sym, self.initial,
                       self.dcd, self.dcy, self.dcl, self.dcr, self.dcm, self.dcn]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FLUTTER_ZONA(Spline):
    """
    Defines data needed to perform flutter, ASE, or a transient response analysis.

    +---------+-----+--------+------+------+-------+-------+-------------+------+
    |    1    |  2  |   3    |  4   |  5   |   6   |   7   |      8      |  9   |
    +=========+=====+========+======+======+=======+=======+=============+======+
    | FLUTTER | SID | METHOD | DENS | MACH | RFREQ | IMETH | NVALUE/OMAX | EPS  |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    | FLUTTER | 19  |   K    | 119  | 219  | 319   |   S   |      5      | 1.-4 |
    +---------+-----+--------+------+------+-------+-------+-------------+------+

    +---------+-------+-------+-------+-------+--------+-------+---------+--------+
    |    1    |   2   |   3   |   4   |   5   |    6   |   7   |    8    |    9   |
    +=========+=======+=======+=======+=======+========+=======+=========+========+
    | FLUTTER | SETID |  SYM  |  FIX  | NMODE | TABDMP | MLIST | CONMLST | NKSTEP |
    +---------+-------+-------+-------+-------+--------+-------+---------+--------+
    | FLUTTER |  100  | SYMM3 |   1   |   0   |  30    |  100  |    0    |   50   |
    +---------+-------+-------+-------+-------+--------+-------+---------+--------+

    """
    type = 'FLUTTER_ZONA'

    def __init__(self, sid, sym, fix, nmode, tabdmp, mlist, conmlst, nkstep=25, comment=''):
        """
        Creates a FLUTTER card, which is required for a flutter (SOL 145)
        analysis.

        Parameters
        ----------
        sid : int
            Unique set identification number. (Integer > 0)
        sym : str
           Character string up to 8 characters with no embedded blanks.
           The first 4 characters can be either 'SYMM' (or 'SYM'), 'ANTI',
           or 'ASYM' that defines the boundary condition of the structural
           finite element model as well as the unsteady aerodynamics, where:
            - SYMM Symmetric boundary condition.
            - ANTI Antisymmetric boundary condition.
            - ASYM Asymmetric boundary condition.
          The last 4 characters are used to specify the interpolation scheme
          for the generalized aerodynamic matrices. They can be either:
           - blank for a cubic spline
           - L for a linear interpolation.
             (such as SYM = 'SYMML', 'ANTIL', or 'ASYML')
           - P for a second-order-polynomial interpolation.
             (such as SYM = 'SYMMP', 'ANTIP', or 'ASYMP')
           - integer for a hybrid cubic spline and linear interpolation scheme.
             (such as SYM = 'SYMM1', 'SYMM2', 'ANTI3', etc.)
           - (Default = SYMML)
        fix : int
           Identification number of a FIXHATM, FIXMATM, FIXMACH, or FIXMDEN
           bulk data card. (Integer > 0)
        nmode : int
            Number of structural modes used in the flutter analysis. (Integer >= 0)
        tabdmp : int
            Identification number of a TABDMP1 bulk data card specifying modal damping as
           a function of natural frequency. (Integer ≥ 0)
        mlist : int
            Identification number of a SET1 or SETADD bulk data card
            specifying a list of normal modes to be omitted from the
            flutter analysis. (Integer >= 0)
        conmlst : int
            Identification number of a CONMLST bulk data card specifying
            a list of identification numbers of the CONM1 bulk data cards
            for mass perturbation. (Integer >= 0) (See Remark 8)
        nkstep : int; default=25
            Number of reduced frequency steps for the reduced-frequency-sweep
            technique of the g-Method flutter solution. (Integer >= 0)

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.sym = sym
        self.fix = fix
        self.nmode = nmode
        self.tabdmp = tabdmp
        self.mlist = mlist
        self.conmlst = conmlst
        self.nkstep = nkstep

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FLUTTER card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        sym = string(card, 2, 'sym')
        fix = integer(card, 3, 'fix')
        nmode = integer(card, 4, 'nmode')
        tabdmp = integer(card, 5, 'tabdmp')
        mlist = integer(card, 6, 'mlist')
        conmlst = integer(card, 7, 'conmlst')
        nkstep = integer_or_blank(card, 8, 'nkstep', 25)
        assert len(card) <= 9, 'len(FLUTTER card) = %i\ncard=%s' % (len(card), card)
        return FLUTTER_ZONA(sid, sym, fix, nmode, tabdmp, mlist, conmlst, nkstep,
                            comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return
        #msg = ', which is required by SPLINE1 eid=%s' % self.eid
        #self.setg_ref = model.Set(self.setg, msg=msg)
        #self.setg_ref.cross_reference_set(model, 'Node', msg=msg)

        #self.panlst_ref = model.zona.panlsts[self.panlst]
        #self.panlst_ref.cross_reference(model)
        #self.aero_element_ids = self.panlst_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors=None):
        return
        #msg = ', which is required by SPLINE1 eid=%s' % self.eid
        #try:
            #self.setg_ref = model.Set(self.setg, msg=msg)
            #self.setg_ref.safe_cross_reference(model, 'Node', msg=msg)
        #except KeyError:
            #model.log.warning('failed to find SETx set_id=%s%s; allowed_sets=%s' % (
                #self.setg, msg, np.unique(list(model.sets.keys()))))

        #try:
            #self.panlst_ref = model.zona.panlsts[self.panlst]
            #self.panlst_ref.safe_cross_reference(model, xref_errors)
            #self.aero_element_ids = self.panlst_ref.aero_element_ids
        #except KeyError:
            #pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        return
        #self.panlst_ref = None
        #self.setg_ref = None

    def convert_to_nastran(self, model):
        raise NotImplementedError()

    def raw_fields(self):
        raise NotImplementedError('FLUTTER - raw_fields')
        #list_fields = ['FLUTTER', self.sid, self.sym, self.fix, self.nmode,
                       #self.tabdmp, self.mlist, self.conmlst, self.nkstep]
        #return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE1_ZONA(Spline):
    """
    Defines an infinite plate spline method for displacements and loads
    transferal between CAERO7 macroelement and structural grid points.

    +---------+-------+-------+------+------+------+----+------+-------+
    |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+======+======+======+====+======+=======+
    | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
    +---------+-------+-------+------+------+------+----+------+-------+
    |         | NELEM | MELEM |      |      |      |    |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+
    | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    """
    type = 'SPLINE1_ZONA'

    def __init__(self, eid, panlst, setg, model=None, cp=None,
                 dz=None, eps=0.01, comment=''):
        """
        Creates a SPLINE1 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.model = model
        self.cp = cp
        self.panlst = panlst
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.panlst_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = string_or_blank(card, 2, 'model')
        cp = blank(card, 3, 'cp')

        panlst = integer(card, 4, 'panlst/setk')
        setg = integer(card, 5, 'setg')
        dz = blank(card, 6, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        return SPLINE1_ZONA(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                            comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)

        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.cross_reference(model)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.safe_cross_reference(model, 'Node', msg=msg)
        except KeyError:
            model.log.warning('failed to find SETx set_id=%s%s; allowed_sets=%s' % (
                self.setg, msg, np.unique(list(model.sets.keys()))))

        try:
            self.panlst_ref = model.zona.panlsts[self.panlst]
            self.panlst_ref.safe_cross_reference(model, xref_errors)
            self.aero_element_ids = self.panlst_ref.aero_element_ids
        except KeyError:
            pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.panlst_ref = None
        self.setg_ref = None

    def convert_to_nastran(self, model):
        """
        +---------+-------+-------+------+------+------+----+------+-------+
        |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
        +=========+=======+=======+======+======+======+====+======+=======+
        | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
        +---------+-------+-------+------+------+------+----+------+-------+
        |         | NELEM | MELEM |      |      |      |    |      |       |
        +---------+-------+-------+------+------+------+----+------+-------+
        | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
        +---------+-------+-------+------+------+------+----+------+-------+

        +---------+------+-------+-------+------+------+----+-----+-------+
        |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
        +=========+======+=======+=======+======+======+====+=====+=======+
        | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
        +---------+------+-------+-------+------+------+----+-----+-------+
        | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
        +---------+------+-------+-------+------+------+----+-----+-------+

        """
        #panlst = 100 # set_aero
        #return SPLINE1(self.eid, panlst, self.setg, model=None, cp=self.cp, dz=self.dz,
                       #eps=0.01, comment=self.comment)
        splines = []
        if not hasattr(self, '_comment'):
            self._comment = ''
        comment = '-' * 72 + '\n' #+ self._comment
        #self._comment = ''
        comment += str(self)
        for panel_groups in self.panlst_ref.panel_groups:
            eid = model.zona.caero_to_name_map[panel_groups]
            caero = model.caeros[eid]
            caero_id = eid
            box1 = caero.eid
            box2 = box1 + caero.npanels - 1
            assert caero.npanels > 0, caero
            #assert box1 > 0 and box2 > 0, 'box1=%s box2=%s' % (box1, box2)
            spline = SPLINE1(eid, caero_id, box1, box2, self.setg, dz=self.dz,
                             method='IPS', usage='BOTH',
                             nelements=10, melements=10, comment=comment)
            spline.validate()
            splines.append(spline)
            comment = ''
        return splines

    def raw_fields(self):
        list_fields = ['SPLINE1', self.eid, self.model, self.cp, self.panlst, self.setg,
                       self.dz, self.eps]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SPLINE2_ZONA(Spline):
    """
    Defines an infinite plate spline method for displacements and loads
    transferal between CAERO7 macroelement and structural grid points.

    +---------+------+-------+------+------+----+-----+-------+-------+
    |    1    |  2   |   3   |  5   |   6  |  6 |  7  |   8   |   9   |
    +=========+======+=======+======+======+====+=====+=======+=======+
    | SPLINE2 | EID  | MODEL | SETK | SETG | DZ | EPS |  CP   | CURV  |
    +---------+------+-------+------+------+----+-----+-------+-------+
    | SPLINE2 | 100  |       |  1   |  10  | 0. |     |       |       |
    +---------+------+-------+------+------+----+-----+-------+-------+

    """
    type = 'SPLINE2_ZONA'

    def __init__(self, eid, panlst, setg, model=None, dz=None, eps=0.01,
                 cp=None, curvature=None, comment=''):
        """
        Creates a SPLINE1 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.model = model
        self.cp = cp
        self.panlst = panlst
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.curvature = curvature
        self.panlst_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = string_or_blank(card, 2, 'model')
        panlst = integer(card, 3, 'panlst/setk')
        setg = integer(card, 4, 'setg')
        dz = blank(card, 5, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        cp = integer_or_blank(card, 7, 'cp', 0)
        curvature = double_or_blank(card, 8, 'curvature', 1.0)
        return SPLINE2_ZONA(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                            curvature=curvature, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.cross_reference(model)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        try:
            msg = ', which is required by SPLINE1 eid=%s' % self.eid
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        except:
            pass
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.safe_cross_reference(model, xref_errors)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.panlst_ref = None
        self.setg_ref = None

    def raw_fields(self):
        list_fields = ['SPLINE2', self.eid, self.model, self.panlst, self.setg,
                       self.dz, self.eps, self.cp, self.curvature]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SPLINE3_ZONA(Spline):
    """
    Defines a 3-D spline for the BODY7 and CAERO7 macroelement.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE3 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE3 | 100  |       |  N/A  |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    """
    type = 'SPLINE3_ZONA'

    def __init__(self, eid, panlst, setg, model=None, cp=None,
                 dz=None, eps=0.01, comment=''):
        """
        Creates a SPLINE3 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.model = model
        self.cp = cp
        self.panlst = panlst
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.panlst_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = blank(card, 2, 'model')
        cp = blank(card, 3, 'cp')

        panlst = integer(card, 4, 'panlst/setk')
        setg = integer(card, 5, 'setg')
        dz = blank(card, 6, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        return SPLINE3_ZONA(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                            comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.cross_reference(model)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        except:
            pass
        self.panlst_ref = model.zona.panlsts[self.panlst]
        self.panlst_ref.safe_cross_reference(model, xref_errors)
        self.aero_element_ids = self.panlst_ref.aero_element_ids

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.panlst_ref = None
        self.setg_ref = None

    def convert_to_nastran(self, model):
        return []

    def raw_fields(self):
        """
        +---------+------+-------+-------+------+----+----+-----+-------+
        |    1    |  2   |   3   |   4   |  5   |  6 |  7 |  8  |   9   |
        +=========+======+=======+=======+======+====+====+=====+=======+
        | SPLINE3 | EID  | CAERO | BOXID | COMP | G1 | C1 | A1  | USAGE |
        +---------+------+-------+-------+------+----+----+-----+-------+
        |         |  G2  |  C2   |  A2   | ---- | G3 | C3 | A2  |  ---  |
        +---------+------+-------+-------+------+----+----+-----+-------+
        |         |  G4  |  C4   |  A4   | etc. |    |    |     |       |
        +---------+------+-------+-------+------+----+----+-----+-------+

        """
        list_fields = ['SPLINE3', self.eid, self.model, self.cp, self.panlst, self.setg,
                       self.dz, self.eps]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
