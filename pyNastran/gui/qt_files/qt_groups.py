
"""
group A
 - elements
 - properties/regions

group B


Nastran
  - type = 'elements', 'properties' (pick one)
  - elements = [1,2,3,4, 4:10]
  - properties = [2,4,6]
  - coords = [1,2,3]

Cart3d
  - type = 'elements', 'regions' (pick one)
  - elements = [1,2,3,4, 4:10]
  - regions = [2,4,6]

Panair
  - type = 'panels'
  - panels = [2,4,6]
"""

self.create_group('nastran', '0', element_id=element_id)  # panel ID
self.create_group('nastran', 'A-Nastran', element_id=element_id, coord_id=3)
self.create_group('nastran', 'B-Nastran', property_id=property_id, coord_id=[3,4])

self.create_group('cart3d', '0', element_id=element_id)  # panel ID
self.create_group('cart3d', 'A', element_id=element_id)
self.create_group('cart3d', 'B', property_id=property_id)  # regions

self.create_group('panair', '0', element_id=element_id)  # panel ID
self.create_group('panair', 'A', element_id=element_id)  # panel ID
self.create_group('panair', 'B', element_id=element_id)  # panel ID

self.remove_groups(['A', 'B'])

#self.set_active_group('0')
self.post_groups(['A', 'B'])

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
    assert len(all_groups) >= 0, 'all_groups is empty'

    all_groups = self.groups # set
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
    self._group_coords[name].union(set(coord_id))

def add_to_group(self, Format, name, element_id=None, property_id=None, coord_id=None):
    assert name in self._group_data
    element_id = self._check_add(Format, name,
                                 element_id=element_id, property_id=property_id,
                                 coord_id=coord_id)
    self._group_elements[name].extend(element_id)
    self._add_coord_id(name, coord_id)

def create_group(self, Format, name,
                 element_id=None, property_id=None, coord_id=None, show=True):
    element_id = self._check_add(Format, name,
                                 element_id=element_id, property_id=property_id,
                                 coord_id=coord_id)

    self.groups.add(name)
    self._group_data[name] = element_id
    self._add_coord_id(name, coord_id)
    self._group_shown[name] = show
