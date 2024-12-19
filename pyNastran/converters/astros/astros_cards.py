from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.cards.coordinate_systems import Coord

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.cards.nodes import GRID, EPOINT, GRDSET, SEQGP, SPOINT
from pyNastran.bdf.cards.coordinate_systems import CORD1C, CORD1R, CORD1S, CORD2C, CORD2R, CORD2S
from pyNastran.bdf.cards.constraints import MPC, MPCADD
from pyNastran.bdf.cards.materials import MAT1, MAT2, MAT8, MAT9
from pyNastran.bdf.cards.elements.mass import CMASS1, CMASS2
from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2
from pyNastran.bdf.cards.elements.shell import CQUAD4, CSHEAR, CTRIA3, QuadShell
from pyNastran.bdf.cards.elements.bars import CBAR, BAROR
from pyNastran.bdf.cards.elements.rods import CROD, CONROD, RodElement
from pyNastran.bdf.cards.elements.mass import CONM1, CONM2, PointMassElement
from pyNastran.bdf.cards.elements.rigid import RBE1, RBE3, RBE2, RigidElementBase, RBAR, RROD
from pyNastran.bdf.cards.elements.elements import GENEL
from pyNastran.bdf.cards.properties.shell import PSHELL, Property, PCOMP, PSHEAR
from pyNastran.bdf.cards.properties.bars import PBAR
from pyNastran.bdf.cards.properties.springs import PELAS
from pyNastran.bdf.cards.properties.mass import PMASS
from pyNastran.bdf.cards.properties.rods import PROD
from pyNastran.bdf.cards.optimization import OptConstraint, DESVAR
from pyNastran.bdf.cards.constraints import SPC, SPC1, SPCADD, SUPORT1
from pyNastran.bdf.cards.aero.aero import AEFACT, AESURF, SPLINE1, SPLINE2, CAERO1, CAERO2, Spline, PAERO1, PAERO2
from pyNastran.bdf.cards.aero.static_loads import TRIM, AEROS
from pyNastran.bdf.cards.aero.dynamic_loads import AERO, Aero, FLFACT, GUST, MKAERO1, MKAERO2
from pyNastran.bdf.cards.loads.dloads import DLOAD
from pyNastran.bdf.cards.loads.static_loads import FORCE, FORCE1, GRAV, LOAD, MOMENT, MOMENT1, PLOAD, PLOAD2, PLOAD4
from pyNastran.bdf.cards.dynamic import FREQ, FREQ1, FREQ2, TF, TSTEP
from pyNastran.bdf.cards.methods import EIGC, EIGR 
from pyNastran.bdf.cards.bdf_sets import SET1, ASET, ASET1, OMIT, OMIT1
from pyNastran.bdf.cards.bdf_tables import TABDMP1, TABLED1
from pyNastran.bdf.cards.thermal.loads import TEMP, TEMPD

class GRID(GRID):
    """
    Copy of NASTRAN GRID class
    """
    def __init__(self, nid, xyz, cp = 0, cd = 0, ps = ''):
        if cp is None:
            cp = 0
        if cd is None:
            cd = 0
        if ps is None:
            ps = 0
        super().__init__(nid, xyz, cp, cd, ps)

class GRDSET(GRDSET):
    """
    Modified Nastran GRDSET
    """
    def __init__(self, cp, cd, ps):
        super().__init__(cp, cd, ps)

class SPOINT(SPOINT):
    """
    Copy of Nastran SPOINT
    """

class EPOINT(EPOINT):
    """
    Modified Nastran EPOINT
    """
    def __init__(self, setid, nids):
        super().__init__(nids)
        self.setid = setid

class MPC(MPC):
    """
    Copy of Nastran MPC
    """

class MPCADD(MPCADD):
    """
    Copy of Nastran MPCADD
    """

class GENEL(GENEL):
    """
    Look into more
    """

class GRAV(GRAV):
    """
    Modified Nastran GRAV
    """
    def __init__(self, sid, scale, N, cid = 0):
        super().__init__(sid, scale, N, cid)

class MKAERO1(MKAERO1):
    """
    Modified Nastran MKAERO1
    """
    def __init__(self, machs, reduced_freqs, symxz, symxy):
        super().__init__(machs, reduced_freqs)
        self.symxy = symxy
        self.symxz = symxz

class MKAERO2(MKAERO2):
    """
    Modified Nastran MKAERO2
    """
    def __init__(self, machs, reduced_freqs, symxz, symxy):
        super().__init__(machs,reduced_freqs)
        self.symxy = symxy
        self.symxz = symxz

class MAT1(MAT1):
    """
    Copy of NASTRAN MAT1 Class
    """

class MAT2(MAT2):
    """
    Copy of NASTRAN MAT2 Class
    """

class MAT8(MAT8):
    """
    Copy of NASTRAN MAT8 Class
    """

class MAT9(MAT9):
    """
    Copy of NASTRAN MAT9 Class
    """

class CMASS1(CMASS1):
    """
    Modified Nastran CMASS1
    """
    def __init__(self, eid, pid, nids, c1 = 0, c2 = 0, tmax:float=1E4):
        super().__init__(eid, pid, nids, c1, c2)
        self.tmax = tmax

class CMASS2(CMASS2):
    """
    Modified Nastran CMASS2
    """
    def __init__(self, eid, mass, nids, c1 = 0, c2 = 0, tmin:float=1E-4, tmax:float=1E4):
        super().__init__(eid, mass, nids, c1, c2)
        self.tmin = tmin
        self.tmax = tmax

class CQUAD4(CQUAD4):
    """
    Copy of NASTRAN CQUAD4 Class
    """
    def __init__(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                 tflag=0, T1=None, T2=None, T3=None, T4=None, comment=''):
        """
        Creates a CQUAD4 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        assert len(nids) == 4, nids
        self.nodes = self.prepare_node_ids(nids)
        self.zoffset = zoffset
        self.theta_mcid = theta_mcid
        self.tflag = tflag
        if self.theta_mcid == None:
            self.theta_mcid = 0.0
        if self.zoffset == None:
            self.zoffset = 0.
        if self.tflag == None:
            self.tflag = 0
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.theta_mcid_ref = None

class BAROR(BAROR):
    """
    Modified Nastran BAROR
    """
    def __init__(self,pid,is_g0,g0,x):
        super().__init__(pid,is_g0,g0,x)

class CBAR(CBAR):
    """
    Copy of NASTRAN CBAR Class
    """

class CROD(CROD):
    """
    Copy of NASTRAN CROD Class
    """
    def __init__(self, eid: int, pid: int, nids: list[int], tmax: Optional[float]=None, comment: str=''):
        """
        Creates a CROD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PROD)
        nids : list[int]
            node ids
        tmax : float; default=None
            Max allowable rod area
        comment : str; default=''
            a comment for the card
        """
        RodElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 2
        self.tmax = tmax
        self.nodes_ref = None
        self.pid_ref = None

class CSHEAR(CSHEAR):
    """
    Modified Nastran CSHEAR
    """
    def __init__(self, eid, pid, nids, tmax:Optional[float]=None):
        super().__init__(eid, pid, nids)
        self.tmax = tmax

class CTRIA3(CTRIA3):
    """
    Copy of Nastran CTRIA3
    """

class CONM1(CONM1):
    """
    Copy of Nastran CONM1
    """

class CONM2(CONM2):
    """
    Modified NASTRAN CONM2 Class
    """
    def __init__(self, eid: int, nid: int, mass: float,
                 cid: int=0,
                 X: Optional[list[int]]=None,
                 I: Optional[list[int]]=None, 
                 tmin: Optional[float]=None,
                 tmax: Optional[float]=None, comment: str=''):
        """
        Creates a CONM2 card

        Parameters
        ----------
        eid : int
           element id
        nid : int
           node id
        mass : float
           the mass of the CONM2
        cid : int; default=0
           coordinate frame of the offset (-1=absolute coordinates)
        X : (3, ) list[float]; default=None -> [0., 0., 0.]
            xyz offset vector relative to nid
        I : (6, ) list[float]; default=None -> [0., 0., 0., 0., 0., 0.]
            mass moment of inertia matrix about the CG
            I11, I21, I22, I31, I32, I33 = I
        comment : str; default=''
            a comment for the card

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        #: Element identification number. (0 < Integer < 100,000,000)
        self.eid = eid

        #: Grid point identification number. (Integer > 0)
        self.nid = nid

        #: Coordinate system identification number.
        #: For CID of -1; see X1, X2, X3 below.
        #: (Integer > -1; Default = 0)
        self.cid = cid

        #: Mass value. (Real)
        self.mass = mass

        if X is None:
            X = np.zeros(3)
        #: Offset distances from the grid point to the center of gravity of
        #: the mass in the coordinate system defined in field 4, unless
        #: CID = -1, in which case X1, X2, X3 are the coordinates, not
        #: offsets, of the center of gravity of the mass in the basic
        #: coordinate system. (Real)
        self.X = np.asarray(X)

        if I is None:
            I = np.zeros(6)
        #: Mass moments of inertia measured at the mass center of gravity in
        #: the coordinate system defined by field 4. If CID = -1, the basic
        #: coordinate system is implied. (Real)
        #: I11, I21, I22, I31, I32, I33 = I
        self.I = np.asarray(I)
        self.tmin = tmin
        self.tmax = tmax
        self.nid_ref = None
        self.cid_ref = None

class CONROD(CONROD):
    """
    Modified Nastran CONROD
    """
    def __init__(self, eid, mid, nids, A=0, j=0, c=0, nsm=0, tmin:float=1E-4, tmax:float=1E4):
        super().__init__(eid, mid, nids, A, j, c, nsm)
        self.tmin = tmin
        self.tmax = tmax

class CELAS1(CELAS1):
    """
    Modified Nastran CELAS1
    """
    def __init__(self,eid,pid,nids,c1:int=0,c2:int=0,tmax:float=1.0E4):
        super().__init__(eid,pid,nids,c1,c2)
        self.tmax = tmax

class CELAS2(CELAS2):
    """
    Modified Nastran CELAS2
    """
    def __init__(self, eid, k, nids, c1=0, c2=0, ge=0, s=0, tmin=None, tmax=None):
        super().__init__(eid, k, nids, c1, c2, ge, s)
        self.tmin = tmin
        self.tmax = tmax

class RBE3(RBE3):
    """
    Copy of NASTRAN RBE3 Class
    """
    def __init__(self, setid: int, eid: int, refgrid: int, refc: str,
                 weights: list[float], comps: list[str], Gijs: list[int],
                 Gmi=None, Cmi=None,
                 comment: str=''):
        """
        Creates an RBE3 element

        Parameters
        ----------
        eid : int
            element id
        refgrid : int
            dependent node
        refc : str
            dependent components for refgrid???
        weights : list[float]
            independent weights for the importance of the DOF
        comps : list[str]
            independent components
            len(comps) = len(weights)
        GiJs : varies
            independent nodes
            list[list[int]]:
                allows for different nodes for the different weights
                len(GiJs) = len(weights)
            list[int, ..., int]:
                intended for a single weight
                This will be expanded into list[list[int]]
        Gmi : list[int]; default=None -> []
            dependent nodes / UM Set
        Cmi : list[str]; default=None -> []
            dependent components / UM Set
        alpha : float; default=0.0
            thermal expansion coefficient
        tref : float; default=0.0
            reference temperature (new in MSC 2021)
        comment : str; default=''
            a comment for the card

        """
        RigidElementBase.__init__(self)
        if comment:
            self.comment = comment
        if Gmi is None:
            Gmi = []
        if Cmi is None:
            Cmi = []

        self.setid = setid
        self.eid = eid
        self.refgrid = refgrid
        self.refc = refc
        self.refgrid_ref = None
        self.Gmi_ref = None
        self.Gijs_ref = None

        if not len(weights) == len(comps) and len(weights) == len(Gijs):
            msg = 'len(weights)=%s len(comps)=%s len(Gijs)=%s' % (
                len(weights), len(comps), len(Gijs))
            raise RuntimeError(msg)

        self.weights = weights
        self.comps = comps
        # allow for Gijs as a list or list of lists
        if len(Gijs) == 0:
            raise RuntimeError(f'RBE3 eid={eid}; Gijs is be empty')

        if isinstance(Gijs[0], integer_types):
            Gijs2 = []
            for Gij in Gijs:
                assert isinstance(Gij, integer_types), 'Gij=%s type=%s' % (Gij, type(Gij))
                Gijs2.append([Gij])
            self.Gijs = Gijs2
        else:
            # default
            self.Gijs = Gijs

        if not len(Gmi) == len(Cmi):
            raise RuntimeError(f'len(Gmi)={len(Gmi):d} len(Cmi)={len(Cmi):d}')
        self.Gmi = Gmi
        self.Cmi = Cmi

        self.nodes_ref = None
        self.pid_ref = None

class RBE2(RBE2):
    """
    Copy of NASTRAN RBE3 Class
    """
    def __init__(self, setid: int, eid: int,
                 gn: int,  # independent
                 cm: str, Gmi: list[int], # dependent
                 comment: str=''):
        """
        Creates an RBE2 element

        Parameters
        ----------
        eid : int
            element id
        gn : int
           Identification number of grid point to which all six independent
           degrees-of-freedom for the element are assigned.
        cm : str
            Component numbers of the dependent degrees-of-freedom in the
            global coordinate system at grid points GMi.
        Gmi : list[int]
            dependent nodes
        alpha : float; default=0.0
            thermal expansion coefficient
        tref : float; default=0.0
            reference temperature
            TREF was added in MSC 2021
        """
        RigidElementBase.__init__(self)
        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid
        self.setid = setid
        #: Identification number of grid point to which all six independent
        #: degrees-of-freedom for the element are assigned. (Integer > 0)
        self.gn = gn

        #: Component numbers of the dependent degrees-of-freedom in the
        #: global coordinate system at grid points GMi. (Integers 1 through
        #: 6 with no embedded blanks.)
        self.cm = cm

        #: Grid point identification numbers at which dependent
        #: degrees-of-freedom are assigned. (Integer > 0)
        if isinstance(Gmi, integer_types):
            Gmi = [Gmi]
        elif isinstance(Gmi, np.ndarray):
            Gmi = Gmi.tolist()
        self.Gmi = Gmi
        #self.nodes_ref = None
        self.Gmi_ref = None
        self.gn_ref = None

class RBE1(RBE1):
    """
    Modified Nastran RBE1
    """
    def __init__(self, setid, eid, Gni, Cni, Gmi, Cmi):
        super().__init__(eid, Gni, Cni, Gmi, Cmi)
        self.setid = setid

class RBAR(RBAR):
    """
    Modified Nastran RBAR
    """
    def __init__(self, setid, eid, nids, cna, cnb, cma, cmb):
        super().__init__(eid, nids, cna, cnb, cma, cmb)
        self.setid = setid

class RROD(RROD):
    """
    Modified Nastran RROD
    """
    def __init__(self, setid, eid, nids, cma = None, cmb = None):
        super().__init__(eid, nids, cma, cmb)
        self.setid = setid

class PBAR(PBAR):
    """
    Look into further
    """

class PCOMP(PCOMP):
    """
    Copy of Nastran PCOMP
    """

class PELAS(PELAS):
    """
    Modified Nastran PELAS
    """
    def __init__(self, pid, k, ge=0, s=0, tmin:float=0.0001):
        super().__init__(pid, k, ge, s)
        self.tmin = tmin

class PSHELL(PSHELL):
    """
    Copy of NASTRAN PSHELL Class
    """
    def __init__(self, pid: int,
                 mid1: Optional[int]=None, t: Optional[float]=None,
                 mid2: Optional[int]=None, twelveIt3: float=1.0,
                 mid3: Optional[int]=None, tst: float=0.833333, nsm: float=0.0,
                 z1: Optional[float]=None, z2: Optional[float]=None,
                 mid4: Optional[int]=None, mcsid: Optional[int]=None, 
                 scsid: Optional[int]=None, zoff: float=0.0, tmin: float=0.0001,
                 comment: str=''):
        """
        Creates a PSHELL card

        Parameters
        ----------
        pid : int
            property id
        mid1 : int; default=None
            defines membrane material
            defines element density (unless blank)
        mid2 : int; default=None
            defines bending material
            defines element density if mid1=None
        mid3 : int; default=None
            defines transverse shear material
            (only defined if mid2 > 0)
        mid4 : int; default=None
            defines membrane-bending coupling material
            (only defined if mid1 > 0 and mid2 > 0; can't be mid1/mid2)
        twelveIt3 : float; default=1.0
            Bending moment of inertia ratio, 12I/T^3. Ratio of the actual
            bending moment inertia of the shell, I, to the bending
            moment of inertia of a homogeneous shell, T^3/12. The default
            value is for a homogeneous shell.
        nsm : float; default=0.0
            non-structural mass per unit area
        z1 / z2 : float; default=None
            fiber distance location 1/2 for stress/strain calculations
            z1 default : -t/2 if thickness is defined
            z2 default : t/2 if thickness is defined
        mcsid : int; default=None
            ID number of material coordiante system
        scsid : int; default=None
            ID number of stress coordinate system
        zoff : float; default=0
            Offset of element reference plane from plane of grid points
        tmin : float; default=0.0001
            Minimum thickness for design
        comment : str; default=''
            a comment for the card

        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        if mid2 == -1:
            mid2 = None

        #: Property ID
        self.pid = pid
        self.mid1 = mid1
        #: Material identification number for bending
        #: -1 for plane strin
        self.mid2 = mid2
        self.mid3 = mid3
        self.mid4 = mid4

        #: thickness
        self.t = t

        #: Scales the moment of interia of the element based on the
        #: moment of interia for a plate
        #:
        #: ..math:: I = \frac{12I}{t^3} I_{plate}
        self.twelveIt3 = twelveIt3

        #: Transverse shear thickness ratio, . Ratio of the shear thickness,
        #: ts/t, to the membrane thickness of the shell, t.
        #: The default value is for a homogeneous shell.
        self.tst = tst

        #: Non-structural Mass
        self.nsm = nsm

        if z1 is None and self.t is not None:
            z1 = -self.t / 2.
        if z2 is None and self.t is not None:
            z2 = self.t / 2.

        self.z1 = z1
        self.z2 = z2

        self.mcsid = mcsid
        self.scsid = scsid

        self.zoff = zoff

        self.tmin = tmin

        if self.t is not None:
            assert self.t >= 0.0, 'PSHELL pid=%s Thickness=%s must be >= 0' % (self.pid, self.t)

        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None
        self.mid4_ref = None

class PSHEAR(PSHEAR):
    """
    Modified Nastran PSHEAR
    """
    def __init__(self, pid, mid, t, nsm = 0):
        super().__init__(pid, mid, t, nsm)

class PROD(PROD):
    """
    Modified Nastran PROD
    """
    def __init__(self, pid, mid, A, j=0, c=0, nsm=0, tmin:float=0.0001):
        super().__init__(pid, mid, A, j, c, nsm)
        self.tmin = tmin

class PLIST(BaseCard):
    """
    ASTROS PLIST card
    """
    def __init__(self, linkid: int, ptype: str, pids: list[int]):
        self.linkid = linkid
        self.ptype = ptype
        self.pids = pids

    def cross_reference(self, model) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PFAST pid=%s' % self.pid
        if self.mcid != -1:
            self.mcid_ref = model.Coord(self.Mcid(), msg)

    def Mcid(self):
        if self.mcid_ref is None:
            return self.mcid
        return self.mcid_ref.cid
    
    def raw_fields(self):
        list_fields = ['PLIST', self.linkid]
        i = 0
        for key, value in sorted(self.params.items()):
            list_fields += [key, value]
            i += 1
            if i == 3:
                list_fields.append(None)
        return list_fields
    
class PLOAD(PLOAD):
    """
    Copy of Nastran PLOAD
    """

class PLOAD2(PLOAD2):
    """
    Copy of Nastran PLOAD2
    """

class PLOAD4(PLOAD4):
    """
    Modified Nastran PLOAD4
    """
    def __init__(self, sid, eid, pressures, cid=0, nvector=None):
        super().__init__(sid, eid, pressures, None, None, cid, nvector)

class PAERO1(PAERO1):
    """
    Copy of Nastran PAERO1
    """

class PAERO2(PAERO2):
    """
    Copy of Nastran PAERO2
    """

class PMASS(PMASS):
    """
    Modified Nastran PMASS
    """
    def __init__(self, pid, mass, tmin):
        super().__init__(pid, mass)
        self.tmin = tmin

class SPC(SPC):
    """
    Copy of Nastran SPC
    """

class SPC1(SPC1):
    """
    Copy of NASTRAN SPC1
    """

class SPCADD(SPCADD):
    """
    Copy of Nastran SPCADD
    """

class SUPORT(SUPORT1):
    """
    NASTRAN suport class
    """

class DCONVMP(OptConstraint):
    """
    Astros DCONVMP card
    """
    def __init__(self, sid: int, st: float, sc: float, ss: float, ptype: str, layrnum: int, pids: list[int], comment: str=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.st = st
        if sc == None:
            sc = st
        else:
            self.sc = sc
        self.ss = ss
        self.ptype = ptype
        self.layrnum = layrnum
        self.pids = pids

class DLOAD(DLOAD):
    """
    Copy of Nastran DLOAD
    """

class LOAD(LOAD):
    """
    Copy of Nastran LOAD
    """

class MOMENT(MOMENT):
    """
    Copy of Nastran MOMENT
    """

class MOMENT1(MOMENT1):
    """
    Copy of Nastran MOMENT1
    """

class DESVARP(DESVAR):
    """
    Modified NASTRAN desvar card
    """
    def __init__(self, desvar_id: int, linkid: int, label: str, layernum: int, layrlst: int,
                 xinit: float, xlb: float=0.001, xub: float=1000.0, comment: str=''):
        """
        Creates a DESVAR card

        Parameters
        ----------
        desvar_id : int
            design variable id
        label : str
            name of the design variable
        xinit : float
            the starting point value for the variable
        xlb : float; default=-1.e20
            the lower bound
        xub : float; default=1.e20
            the lower bound
        delx : float; default=1.e20
            fractional change allowed for design variables during
            approximate optimization
            NX  if blank : take from DOPTPRM; otherwise 1.0
            MSC if blank : take from DOPTPRM; otherwise 0.5
        ddval : int; default=0
            int : DDVAL id
                  allows you to set discrete values
            0/None : continuous
        comment : str; default=''
            a comment for the card

        """
        OptConstraint.__init__(self)
        if comment:
            self.comment = comment
        self.desvar_id = desvar_id
        self.linkid = linkid
        #: user-defined name for printing purposes
        self.label = label
        #xinit = np.clip(xinit, xlb, xub)
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        self.layernum = layernum
        self.layrlst = layrlst
        #assert ' ' not in label.rstrip(), self.get_stats()

class DESVARS(DESVAR):
    """
    Modified NASTRAN desvar card
    """
    def __init__(self, desvar_id: int, shapeid: int, label: str, layernum: int, layrlst: int,
                 xinit: float, xlb: float=0.001, xub: float=1000.0, comment: str=''):
        """
        Creates a DESVAR card

        Parameters
        ----------
        desvar_id : int
            design variable id
        label : str
            name of the design variable
        xinit : float
            the starting point value for the variable
        xlb : float; default=-1.e20
            the lower bound
        xub : float; default=1.e20
            the lower bound
        delx : float; default=1.e20
            fractional change allowed for design variables during
            approximate optimization
            NX  if blank : take from DOPTPRM; otherwise 1.0
            MSC if blank : take from DOPTPRM; otherwise 0.5
        ddval : int; default=0
            int : DDVAL id
                  allows you to set discrete values
            0/None : continuous
        comment : str; default=''
            a comment for the card

        """
        OptConstraint.__init__(self)
        if comment:
            self.comment = comment
        self.desvar_id = desvar_id
        self.shapeid = shapeid
        # self.linkid = shapeid

        #: user-defined name for printing purposes
        self.label = label
        #xinit = np.clip(xinit, xlb, xub)
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        self.layernum = layernum
        self.layrlst = layrlst
        #assert ' ' not in label.rstrip(), self.get_stats()

class FLFACT(FLFACT):
    """
    Copy of Nastran FLFACT
    """

class FORCE(FORCE):
    """
    Copy of Nastran Force
    """

class FORCE1(FORCE1):
    """
    Copy of Nastran Force1
    """

class GUST(GUST):
    """
    Look into more
    """

class AEFACT(AEFACT):
    """
    Copy of NASTRAN AEFACT
    """

class AESURF(AESURF):
    """
    Modified NASTRAN AESURF
    """
    def __init__(self,label: str, type: str, acid: int, fboxid: int, lboxid: int, cid: Optional[int] = None):
        AESURF.__init__(self,label,type,acid,fboxid,lboxid,cid)

class AIRFOIL(BaseCard):
    """
    ASTROS AIRFOIL Class
    """
    def __init__(self, acid: int, cmpnt: str, chord: int, cam: int, radius: float,
                 x1:float, y1:float, z1: float, x12: Optional[float]=None, cp: Optional[int]=None,
                 usothk: Optional[int]=None, lso: Optional[int]=None, ipanel: Optional[int]=None):
        self.acid = acid
        self.cmpnt = cmpnt
        self.chord = chord
        self.usothk = usothk
        self.ls0 = lso
        self.cam = cam
        self.radius = radius
        self.xyz = [x1,y1,z1]
        self.x12 = x12
        self.ipanel = ipanel

    def raw_fields(self):
        list_fields = ['AIRFOIL', self.acid]
        i = 0
        for key, value in sorted(self.params.items()):
            list_fields += [key, value]
            i += 1
            if i == 3:
                list_fields.append(None)
        return list_fields

class ASET(ASET):
    """
    Modified Nastran ASET
    """
    def __init__(self,setid,ids,comps):
        super().__init__(ids,comps)
        self.setid = setid

class ASET1(ASET1):
    """
    Modified Nastran ASET1
    """
    def __init__(self,setid,comp,ids):
        super().__init__(ids,comp)
        self.setid = setid

class OMIT(OMIT):
    """
    Modified Nastran OMIT
    """
    def __init__(self, setid, ids, components):
        super().__init__(ids, components)
        self.setid = setid

class OMIT1(OMIT1):
    """
    Modified Nastran OMIT1
    """
    def __init__(self, setid, ids, components):
        super().__init__(ids, components)
        self.setid = setid

class CAERO(CAERO1):
    """
    Used to convert CAERO6, AIRFOIL, and AEFACT to Nastran CAERO1
    """

class CAERO1(CAERO1):
    """
    Copy of Nastran CAERO1
    """

class CAERO2(CAERO2):
    """
    Copy of Nastran CAERO2
    """

class CAERO6(BaseCard):
    """
    ASTROS CAERO6 card
    """
    def __init__(self, acid: int, cmpnt: str, igrp: int, cp: Optional[int] = None,
                 lchord: Optional[int] = None, lspan: Optional[int] = None):
        super().__init__()
        self.acid = acid
        self.cmpnt = cmpnt
        self.igrp = igrp
        self.cp = cp
        self.lchord = lchord
        self.lspan = lspan

    def raw_fields(self):
        list_fields = ['CAERO6', self.acid]
        i = 0
        for key, value in sorted(self.params.items()):
            list_fields += [key, value]
            i += 1
            if i == 3:
                list_fields.append(None)
        return list_fields

class AERO(AERO):
    """
    Modified Nastran Aero card
    """
    def __init__(self, cref, rhoref, acsid=0):
        Aero.__init__(self)
        velocity = 0
        super().__init__(velocity,cref,rhoref,acsid,sym_xz=0,sym_xy=0)
        self.cref = cref
        self.rho_ref = rhoref
        self.acsid = acsid

class AEROS(AEROS):
    """
    Copy of NASTRAN AEROS
    """

class SPLINE1(SPLINE1):
    """
    Modified NASTRAN SPLINE1
    """

    def __init__(self, eid: int, caero: int, box1: int, box2: int, setg: int,
                 dz: float=0., method: str='IPS',
                 usage: str='BOTH', nelements: int=10,
                 melements: int=10, comment: str=''):
        """
        Creates a SPLINE1, which defines a surface spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements : int; default=10
            The number of FE elements along the local spline x-axis if
            using the FPS option
        melements : int; default=10
            The number of FE elements along the local spline y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        Spline.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.caero = caero
        self.box1 = box1
        self.box2 = box2
        self.setg = setg
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelements = nelements
        self.melements = melements
        if self.dz == None:
            self.dz = 0.0
        self.caero_ref = None
        self.setg_ref = None

class SPLINE2(SPLINE2):
    """
    Modified Nastran SPLINE2
    """
    def __init__(self, eid, caero, box1, box2, setg, dz = 0, dtor = 1, cid = 0, dthx = 0, dthy = 0):
        super().__init__(eid, caero, box1, box2, setg, dz, dtor, cid, dthx, dthy)

class FREQ(FREQ):
    """
    Copy of Nastran FREQ
    """

class FREQ1(FREQ1):
    """
    Copy of Nastran FREQ
    """

class FREQ2(FREQ2):
    """
    Copy of Nastran FREQ
    """

class TF(TF):
    """
    Copy of Nastran TF
    """

class TSTEP(TSTEP):
    """
    Copy of Nastran TSTEP
    """

class CORD1C(CORD1C):
    """
    Copy of Nastran CORD1C
    """

class CORD1R(CORD1R):
    """
    Copy of Nastran CORD1R
    """

class CORD1S(CORD1S):
    """
    Copy of Nastran CORD1S
    """

class CORD2C(CORD2C):
    """
    Copy of Nastran CORD2C
    """

class CORD2R(CORD2R):
    """
    Copy of Nastran CORD2R
    """

class CORD2S(CORD2S):
    """
    Copy of Nastran CORD2S
    """

class TRIM(TRIM):
    """
    Modified Nastran TRIM
    """
    def __init__(self, sid, mach, q, labels, uxs, trmtyp:Optional[str]=None, effid:Optional[int]=None, V0:Optional[float]=None):
        super().__init__(sid, mach, q, labels, uxs)
        self.trmtyp = trmtyp
        self.effid = effid
        self.V0 = V0

class EIGC(EIGC):
    """
    Copy of Nastran EIGC
    """

class EIGR(EIGR):
    """
    Modified Nastran EIGR
    """
    def __init__(self, sid, method='LAN', f1=None, f2=None, ne=None, nd=None, crit=None, norm='MASS', G=None, C=None, E:float=1E-10):
        super().__init__(sid, method, f1, f2, ne, nd, crit, norm, G, C)
        self.E = E

class TABDMP1(TABDMP1):
    """
    Copy of Nastran TABDMP1
    """

class TABLED1(TABLED1):
    """
    amodified Nastran TABLED1
    """
    def __init__(self, tid, x, y):
        super().__init__(tid, x, y)

class TEMP(TEMP):
    """
    Copy of Nastran TEMP
    """

class TEMPD(TEMPD):
    """
    Copy of Nastran TEMPD
    """

class SEQGP(SEQGP):
    """
    Copy of Nastran SEQGP
    """

class SET1(SET1):
    """
    Copy of NASTRAN SET1
    """