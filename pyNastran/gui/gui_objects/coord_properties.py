class CoordProperties:
    def __init__(self, label: str, coord_type: str, is_visible: bool,
                 scale: float) -> None:
        self.label = label
        #self.axes = axes
        self.coord_type = coord_type
        self.is_visible = is_visible
        self.representation = 'coord'
        self.scale = scale

    def __repr__(self) -> str:
        msg = 'CoordProperties(label=%r, coord_type%r, is_visible=%s, scale=%s)' % (
            self.label, self.coord_type, self.is_visible, self.scale)
        return msg
