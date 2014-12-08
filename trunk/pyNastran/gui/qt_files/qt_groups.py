
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
  - type = 'panels' (pick one)
  - panels = [2,4,6]


"""

self.define_group('nastran', '0', element_id=element_id)  # panel ID
self.define_group('nastran', 'A-Nastran', element_id=element_id)
self.define_group('nastran', 'B-Nastran', property_id=property_id)

self.define_group('cart3d', '0', element_id=element_id)  # panel ID
self.define_group('cart3d', 'A', element_id=element_id)
self.define_group('cart3d', 'B', property_id=property_id)  # regions

self.define_group('panair', '0', element_id=element_id)  # panel ID
self.define_group('panair', 'A', element_id=element_id)  # panel ID
self.define_group('panair', 'B', element_id=element_id)  # panel ID

self.remove_groups(['A', 'B'])
self.set_active_group('0')

def define_group(self, Format, name, element_id=None, property_id=None):
    pass