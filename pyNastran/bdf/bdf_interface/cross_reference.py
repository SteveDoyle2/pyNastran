"""
Links up the various cards in the BDF.

For example, with cross referencing...

.. code-block:: python

   >>> model = BDF()
   >>> model.read_bdf(bdf_filename, xref=True)

   >>> nid1 = 1
   >>> node1 = model.nodes[nid1]
   >>> node.nid
   1

   >>> node.xyz
   [1., 2., 3.]

   >>> node.Cid()
   3

   >>> node.cid
   3

   >>> node.cid_ref
   CORD2S, 3, 1, 0., 0., 0., 0., 0., 1.,
           1., 0., 0.
   # get the position in the global frame
   >>> node.get_position()
   [4., 5., 6.]

   # get the position with respect to another frame
   >>> node.get_position_wrt(model, cid=2)
   [4., 5., 6.]


Without cross referencing...

.. code-block:: python

   >>> model = BDF()
   >>> model.read_bdf(bdf_filename, xref=True)

   >>> nid1 = 1
   >>> node1 = model.nodes[nid1]
   >>> node.nid
   1

   >>> node.xyz
   [1., 2., 3.]

   >>> node.Cid()
   3

   >>> node.cid
   3

   >>> node.cid_ref
   None

   # get the position in the global frame
   >>> node.get_position()
   Error!

Cross-referencing allows you to easily jump across cards and also helps
with calculating things like position, area, and mass.  The BDF is designed
around the idea of cross-referencing, so it's recommended that you use it.

"""
# pylint: disable=R0902,R0904,R0914
import traceback
from typing import Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.cross_reference_obj import CrossReference


class XrefMesh(BDFAttributes):
    """Links up the various cards in the BDF."""
    def __init__(self) -> None:
        """The main BDF class defines all the parameters that are used."""
        BDFAttributes.__init__(self)

    # def geom_check(self):
        # """
        # Performs various geometry checks
          # 1.  nodal uniqueness on elements
        # """
        # for elem in model.elements:
            # elem.check_unique_nodes()

    def cross_reference(self,
                        xref: bool=True,
                        xref_nodes: bool=True,
                        xref_elements: bool=True,
                        xref_nodes_with_elements: bool=False,
                        xref_properties: bool=True,
                        xref_masses: bool=True,
                        xref_materials: bool=True,
                        xref_loads: bool=True,
                        xref_constraints: bool=True,
                        xref_aero: bool=True,
                        xref_sets: bool=True,
                        xref_optimization: bool=True,
                        word: str='') -> None:
        """
        Links up all the cards to the cards they reference

        Parameters
        ----------
        xref : bool; default=True
           cross references the model
        xref_nodes : bool; default=True
           set cross referencing of nodes/coords
        xref_element : bool; default=True
           set cross referencing of elements
        xref_properties : bool; default=True
           set cross referencing of properties
        xref_masses : bool; default=True
           set cross referencing of CMASS/PMASS
        xref_materials : bool; default=True
           set cross referencing of materials
        xref_loads : bool; default=True
            set cross referencing of loads
        xref_constraints : bool; default=True
            set cross referencing of constraints
        xref_aero : bool; default=True
            set cross referencing of CAERO/SPLINEs
        xref_sets : bool; default=True
            set cross referencing of SETx
        word : str; default=''
            model flag

        To only cross-reference nodes:

        .. code-block:: python

           model = BDF()
           model.read_bdf(bdf_filename, xref=False)
           model.cross_reference(xref=True, xref_loads=False, xref_constraints=False,
                                            xref_materials=False, xref_properties=False,
                                            xref_aero=False, xref_masses=False,
                                            xref_sets=False)

        .. warning:: be careful if you call this method with False values
        """
        if not xref:
            return
        self.log.debug("Cross Referencing%s..." % word)
        xref_obj: CrossReference = self.xref_obj
        #xref_obj.model
        #xref_obj.model.log
        if xref_nodes:
            xref_obj.cross_reference_nodes()
            xref_obj.cross_reference_coordinates()

        if xref_elements:
            xref_obj.cross_reference_bolts()
            xref_obj.cross_reference_elements()
            xref_obj.cross_reference_rigid_elements()
        if xref_properties:
            xref_obj.cross_reference_properties()
        if xref_masses:
            xref_obj.cross_reference_masses()
        if xref_materials:
            xref_obj.cross_reference_materials()

        if xref_aero:
            xref_obj.cross_reference_aero()
        if xref_constraints:
            xref_obj.cross_reference_constraints()
        if xref_loads:
            xref_obj.cross_reference_loads()
        if xref_sets:
            xref_obj.cross_reference_sets()
        if xref_optimization:
            xref_obj.cross_reference_optimization()
        if xref_nodes_with_elements:
            xref_obj.cross_reference_nodes_with_elements()
        xref_obj.cross_reference_contact()
        xref_obj.cross_reference_superelements()
        #self.case_control_deck.cross_reference(self)
        self.pop_xref_errors()

        for super_key, superelement in sorted(self.superelement_models.items()):
            word = ''
            if isinstance(super_key, int):
                word = f' (Superelement {super_key})'
            else:
                word = f'{super_key[0]}={super_key[1]}'

            superelement.cross_reference(
                xref=xref, xref_nodes=xref_nodes, xref_elements=xref_elements,
                xref_nodes_with_elements=xref_nodes_with_elements,
                xref_properties=xref_properties, xref_masses=xref_masses,
                xref_materials=xref_materials, xref_loads=xref_loads,
                xref_constraints=xref_constraints, xref_aero=xref_aero,
                xref_sets=xref_sets, xref_optimization=xref_optimization,
                word=word)

    def _create_superelement_from_sebulk(self, sebulk,
                                         seid: int, rseid: int) -> None:
        """helper for sebulk"""
        #C:\MSC.Software\MSC.Nastran\msc20051\nast\tpl\see103q4.dat
        model = self
        ref_model = model.superelement_models[('SUPER', rseid, '')]
        if sebulk.superelement_type == 'MIRROR':
            from pyNastran.bdf.mesh_utils.mirror_mesh import bdf_mirror_plane
            #print('creating superelement %s from %s' % (seid, rseid))
            sempln = model.sempln[seid]
            plane = np.array([node.get_position() for node in sempln.nodes_ref])

            # What about seloc on the primary and sempln+seloc on the secondary?
            #  - move the primary
            #  - then apply the mirror to make the secondary
            #  - then move the secondary
            #
            # Or what about sempln+seloc on the tertiary?
            #
            # this is fine for the secondary
            if rseid in model.seloc:
                # I think this is wrong...
                seloc = model.seloc[rseid]
                plane = seloc.transform(model, plane)

            ref_model, mirror_model, unused_nid_offset, unused_eid_offset = bdf_mirror_plane(
                ref_model, plane, mirror_model=None, log=None, debug=True, use_nid_offset=False)
            mirror_model.properties = ref_model.properties
            mirror_model.materials = ref_model.materials
            new_model = mirror_model
        elif sebulk.Type in ['MANUAL', 'PRIMARY', 'COLLCTR', 'EXTERNAL']:
            self.log.info('skipping:\n%s' % sebulk)
            new_model = None
        else:  # pragma: no cover
            raise NotImplementedError(sebulk)
        return new_model

    def _uncross_reference_superelements(self) -> None:
        """cross references the superelement objects"""
        model: BDF = self
        for unused_seid, csuper in model.csuper.items():
            csuper.uncross_reference()
        for unused_seid, csupext in model.csupext.items():
            csupext.uncross_reference()

        for unused_seid, sebulk in model.sebulk.items():
            sebulk.uncross_reference()
        for unused_seid, sebndry in model.sebndry.items():
            sebndry.uncross_reference()
        for unused_seid, seconct in model.seconct.items():
            seconct.uncross_reference()
        for unused_seid, seelt in model.seelt.items():
            seelt.uncross_reference()
        for unused_seid, seexcld in model.seexcld.items():
            seexcld.uncross_reference()
        for unused_seid, selabel in model.selabel.items():
            selabel.uncross_reference()
        for unused_seid, seloc in self.seloc.items():
            seloc.uncross_reference()
        for unused_seid, seload in model.seload.items():
            seload.uncross_reference()
        for unused_seid, sempln in model.sempln.items():
            sempln.uncross_reference()
        for unused_seid, setree in model.setree.items():
            setree.uncross_reference()

    def get_point_grids(self, nodes: list[Any], msg: str='') -> None:
        """gets GRID, POINT cards"""
        nodes_ref = []
        missing_nids = []
        model: BDF = self
        for nid in nodes:
            if nid in model.nodes:
                node = model.nodes[nid]
            elif nid in model.points:
                node = model.points[nid]
            else:
                missing_nids.append(nid)
                continue
            nodes_ref.append(node)
        if missing_nids:
            raise KeyError('missing GRID/POINT nids=%s%s' % (missing_nids, msg))
        return nodes_ref

    def superelement_nodes(self, seid: int, nodes: list[Any], msg: str='') -> None:
        if seid == 0:
            return self.Nodes(nodes, msg=msg)
        try:
            superelement = self.superelement_models[('SUPER', seid, '')]
        except KeyError:
            keys = list(self.superelement_models.keys())
            raise KeyError('cant find superelement=%i%s; seids=%s' % (seid, msg, keys))
        return superelement.Nodes(nodes, msg=msg)

    def geom_check(self, geom_check: bool, xref: bool) -> None:  # pragma: no cover
        """
        what about xref?
        """
        if geom_check:
            if xref:
                for unused_eid, element in self.elements.values():
                    #element.Mass()
                    element._verify(xref=True)
                #if 'GEOMCHECK' in self.params:  # should this be an executive control parameter?
                    #for eid, element in model.elements:
                        #element._verify()
            else:
                for unused_eid, element in self.elements.values():
                    element.verify_unique_node_ids()
                    element._verify(xref=False)

            # aspect ratio - ratio between element edges
            # warping - how planar is a face
            # taper - split a quad into 2 triangles and compare the area
            # skew - an angle, measures how skewed an element face is by drawing lines
            #        between midpoints of elements edges, finding the smallest angle
            #        between the intersecting lines and subtracting that from 90 degrees
            # Jacobian - how much does element deviate from the ideal shape by taking the
            #            determinant of the Jacobian matrix
            # quad skew <= 30.
            # quad warp >= 0.05
            # quad taper >= 0.5
            # quad iamin <= 30.
            # quad iamax >= 150.

            # tria skew <= 10.
            # tria iamax <= 160.

            # tetra ar >= 100.
            # tetra elpr <= 0.5
            # tetra detj <= 0.

            # hex ar >= 100.
            # hex elpr <= 0.5
            # hex detj <= 0.
            # hex warp <= 0.707

            # penta ar >= 100.
            # penta elpr <= 0.5
            # penta detj <= 0.
            # penta warp <= 0.707

            # pyram ar >= 100.
            # pyram elpr <= 0.5
            # pyram detj <= 0.
            # pyram warp <= 0.707
