class CoordProperties(object):
    def __init__(self, label, Type, is_visible, scale):
        self.label = label
        #self.axes = axes
        self.Type = Type
        self.is_visible = is_visible
        self.representation = 'coord'
        self.scale = scale
