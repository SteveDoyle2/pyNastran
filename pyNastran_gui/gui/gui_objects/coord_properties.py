class CoordProperties(object):
    def __init__(self, label, coord_type, is_visible, scale):
        # type: (str, str, bool, float) -> None
        self.label = label
        #self.axes = axes
        self.coord_type = coord_type
        self.is_visible = is_visible
        self.representation = 'coord'
        self.scale = scale

    def __repr__(self):
        # type: () -> str
        msg = 'CoordProperties(label=%r, coord_type%r, is_visible=%s, scale=%s)' % (
            self.label, self.coord_type, self.is_visible, self.scale)
        return msg
