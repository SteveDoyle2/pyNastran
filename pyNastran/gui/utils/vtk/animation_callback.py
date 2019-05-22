"""
defines:
 - AnimationCallback
"""
from itertools import cycle

class AnimationCallback:
    """
    http://www.vtk.org/Wiki/VTK/Examples/Python/Animation
    """
    def __init__(self, parent, scales, phases,
                 icases_fringe, icases_disp, icases_vector,
                 animate_fringe, animate_vector,
                 min_value, max_value):
        """
        creates AnimationCallback
        """
        self.parent = parent
        self.timer_count = 0
        self.cycler = cycle(range(len(icases_disp)))

        self.icase_fringe0 = -1
        self.icase_disp0 = -1
        self.icase_vector0 = -1
        self.ncases = len(icases_disp)

        self.scales = scales
        self.phases = phases
        self.icases_fringe = icases_fringe
        self.icases_disp = icases_disp
        self.icases_vector = icases_vector
        self.animate_fringe = animate_fringe
        self.animate_vector = animate_vector

        self.min_value = min_value
        self.max_value = max_value
        self.scale_max = max(abs(self.scales.max()), abs(self.scales.min()))
        #self.isteps = isteps

    def execute(self, obj, unused_event):
        """creates the ith frame"""
        unused_iren = obj
        i = self.timer_count % self.ncases
        #j = next(self.cycler)
        icase_fringe = self.icases_fringe[i]
        icase_disp = self.icases_disp[i]
        icase_vector = self.icases_vector[i]
        scale = self.scales[i]
        phase = self.phases[i]
        normalized_frings_scale = scale / self.scale_max
        is_valid = self.parent.animation_update(
            self.icase_fringe0, self.icase_disp0, self.icase_vector0,
            icase_fringe, icase_disp, icase_vector,
            scale, phase,
            self.animate_fringe, self.animate_vector,
            normalized_frings_scale,
            self.min_value, self.max_value)
        if not is_valid:
            self.parent.stop_animation()

        self.icase_disp0 = icase_disp
        self.icase_fringe0 = icase_fringe
        self.icase_vector0 = icase_vector

        self.parent.vtk_interactor.Render()
        self.timer_count += 1
