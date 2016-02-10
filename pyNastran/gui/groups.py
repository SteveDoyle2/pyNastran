# Groups - not done
class GroupImport(object):

    def __init__(self):
        pass

    def clear_groups(self):
        all_groups = self.groups # set
        self.remove_groups(all_groups)

    def remove_groups(self, groups):
        assert isinstance(groups, list), type(groups)
        assert len(groups) >= 0, 'groups is empty'
        assert len(all_groups) >= 0, 'all_groups is empty'

        all_groups = self.groups # set
        for group in all_groups:
            if group in groups:
                self.remove_group(group)

    def remove_group(self, group):
        if group not in all_groups:
            raise RuntimeError('group=%r not found' % group)

    def show_group(self, name):
        self._group_shown[name] = True

    def hide_group(self, name):
        self._group_shown[name] = False

    def post_groups(self, groups):
        assert isinstance(groups, list), type(groups)
        assert len(groups) >= 0, 'groups is empty'

        all_groups = self.groups # set
        assert len(all_groups) >= 0, 'all_groups is empty'

        for group in all_groups:
            if group in groups:
                self.show_group(group)
            else:
                self.show_group(group)

    def _check_add(self, Format, name, element_id=None, property_id=None, coord_id=None):
        if element_id is None and property_id is None:
            raise RuntimeError('either element_id or property_id must be set')
        if isinstance(element_id, int):
            element_id = [element_id]
        if isinstance(property_id, int):
            property_id = [property_id]

        if Format == 'nastran':
            if property_id:
                element_id = self.model.get_element_id_by_property_id(property_id)

        elif Format == 'cart3d':
            if property_id:
                element_id = self.model.get_gelement_id_by_region_id(property_id)
        elif Format == 'panair':
            if element_id is None:
                raise RuntimeError('element_id must be set for panair')
        else:
            msg = "Format=%r is not supported; use 'nastran', 'cart3d', 'panair'" % Format
            raise NotImplementedError(msg)

        if coord_id is not None and Format != 'nastran':
            raise RuntimeError('coord_id must be None for format=%r' % Format)

        element_id = asarray(element_id)
        return element_id

    def _add_coord_id(self, name, coord_id):
        if coord_id is None:
            coord_id = set([0])
        elif isinstance(coord_id, int):
            coord_id = set([coord_id])
        else:
            for cid in coord_id:
                assert isinstance(cid, int), type(cid)
        if name in self._group_coords:
            self._group_coords[name].union(set(coord_id))
        else:
            self._group_coords[name] = set(coord_id)

    def _create_grid_mapper(self, name):
        self.grid = vtk.vtkUnstructuredGrid()
        self.grid_mapper = vtk.vtkDataSetMapper()
        self.grid_mapper.SetInput(self.grid)

        geometry_actor = vtk.vtkActor()
        geometry_actor.SetMapper(self.grid_mapper)
        geometry_actor.GetProperty().SetDiffuseColor(1, 0, 0)  # red
        self.rend.AddActor(geometry_actor)


class Groups(object):
    def __init__(self):
        self.nNodes = None
        self.model = None
        self.grid = None

    def add_to_group(self, Format, name, element_id=None, property_id=None, coord_id=None):
        assert name in self._group_elements
        element_id = self._check_add(Format, name,
                                     element_id=element_id, property_id=property_id,
                                     coord_id=coord_id)
        self._group_elements[name] = hstack([self._group_elements[name],
                                             element_id])
        self._add_coord_id(name, coord_id)

    def create_group(self, Format, name,
                     element_id=None, property_id=None, coord_id=None, show=True):
        element_id = self._check_add(Format, name,
                                     element_id=element_id, property_id=property_id,
                                     coord_id=coord_id)

        self.groups.add(name)
        self._group_elements[name] = element_id
        self._add_coord_id(name, coord_id)
        self._group_shown[name] = show
        if Format == 'nastran':
            self._create_nastran_group(name, self.model, element_id)
        elif Format == 'cart3d':
            self._create_cart3d_group(name, self.model, element_id)
        #elif Format == 'panair':
            #self._create_panair_group(name, self.model, element_id)

    def _create_nastran_group(self, name, model, element_id):
        pass

    def _create_cart3d_group(self, name, model, element_id):
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(self.nNodes)
        nelements = len(element_id)
        self.grid.Allocate(nelements, 1000)

        nodes = self.model.get_nodes_associated_with_elements(element_id)
        nodes.sort()

        nid = 0
        all_nodes = self.nodes
        for i in all_nodes:
            #if nid in nodes:
            points.InsertPoint(nid, all_nodes[i, :])
            nid += 1

        from vtk import vtkTriangle
        for eid in element_id:
            elem = vtkTriangle()
            node_ids = elements[eid, :]
            elem_nodes = searchsorted(nodes, node_ids)
            elem.GetPointIds().SetId(0, elem_nodes[0])
            elem.GetPointIds().SetId(1, elem_nodes[1])
            elem.GetPointIds().SetId(2, elem_nodes[2])
            self.grid.InsertNextCell(5, elem.GetPointIds())

        self.grid[name].SetPoints(points)
        self.grid[name].Modified()
        self.grid[name].Update()

