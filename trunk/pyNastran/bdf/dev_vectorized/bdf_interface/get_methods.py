class GetMethods(object):
    def __init__(self):
        pass

    def Coord(self, cid, msg=''):
        try:
            return self.coords[cid]
        except KeyError:
            raise KeyError('cid=%s not found%s.  Allowed Cids=%s'
                           % (cid, msg, self.coordIDs()))

    def get_elements(self, *element_ids):
        elements = self.elements.get_elements(element_ids)
        if len(elements) != len(element_ids):
            eids = set([element.element_id if element is not None else None for element in elements])
            diff = set(element_ids).difference(eids)
            msg = 'The following element IDs could not be found %s' % list(diff)
            raise RuntimeError(msg)
        return elements

    def get_properties(self, *property_ids):
        properties = self.elements.get_properties(property_ids)
        if len(properties) != len(property_ids):
            pids = set([prop.property_id if prop is not None else None for prop in properties])
            diff = set(property_ids).difference(pids)
            msg = 'The following property IDs could not be found %s' % list(diff)
            raise RuntimeError(msg)
        return properties
