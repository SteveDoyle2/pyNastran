"""
Defines:
 - write_gif(gif_filename, png_filenames, time=2.0, fps=30,
              onesided=True, nrepeat=0,
              delete_images=False, make_gif=True)
"""
from __future__ import print_function
import os
import numpy as np
try:
    import imageio
    IS_IMAGEIO = True
except ImportError:
    IS_IMAGEIO = False
from pyNastran.utils import integer_types


def setup_animation(scale, istep=None,
                    animate_scale=True, animate_phase=False, animate_time=False,
                    icase=None,
                    icase_start=None, icase_end=None, icase_delta=None,
                    time=2.0, onesided=True,
                    fps=30):
    """
    helper method for ``make_gif``

    Returns
    -------
    phases : List[float]
        List[float] : the phase angles
    icases : List[int]
        List[int] : the icases to run
    isteps : List[int]
        nominal isteps = [0, 1, 2, 3, 4, ..., nframes]
        we can analyze pictures [1, 3, 4] by providing a subset
    scales : List[float]
        the displacement scale factor; true scale
    analysis_time : float
        the time that needs to be simulated for the analysis; not the runtime
    """
    if animate_time or animate_phase:
        onesided = True
    analysis_time = get_analysis_time(time, onesided)

    nframes = int(analysis_time * fps)
    scales = None
    phases = None
    if animate_scale:
        # TODO: we could start from 0 deflection, but that's more work
        # TODO: we could do a sine wave, but again, more work
        icases = icase
        scales = np.linspace(-scale, scale, num=nframes, endpoint=True)
        isteps = np.linspace(0, nframes, num=nframes, endpoint=True, dtype='int32')
        phases = [0.] * nframes
        assert len(scales) == len(isteps), 'nscales=%s nsteps=%s' % (len(scales), len(isteps))
        assert len(phases) == len(isteps), 'nphases=%s nsteps=%s' % (len(phases), len(isteps))
    elif animate_phase:
        # animate phase
        icases = icase
        phases = np.linspace(0., 360, num=nframes, endpoint=False)
        isteps = np.linspace(0, nframes, num=nframes, endpoint=False, dtype='int32')
        scales = [scale] * nframes
        assert len(phases) == len(isteps), 'nphases=%s nsteps=%s' % (len(phases), len(isteps))
        assert len(scales) == len(isteps), 'nscales=%s nsteps=%s' % (len(scales), len(isteps))
    elif animate_time:
        icases = np.arange(icase_start, icase_end+1, icase_delta)
        #min_value = 0.
        #max_value = 1.46862
        nfiles = len(icases)

        # specifying fps and dt makes the problem overdefined
        # assuming dt
        #
        # TDOO: change this to stepping similar to icase_delta
        #       icases = icases[::5]
        fps = nfiles / time

        # our scale will be constant
        # phases is just None
        scales = [scale] * nfiles
        assert len(icases) == nfiles, 'len(icases)=%s nfiles=%s' % (len(icases), nfiles)
        phases = None
        assert len(scales) == len(icases), 'nscales=%s len(icases)=%s' % (len(scales), len(icases))

        # TODO: this isn't maintained...
        #assert nframes == nfiles, 'nframes=%s nfiles=%s' % (nframes, nfiles)

        isteps = np.linspace(0, nfiles, num=nfiles, endpoint=True, dtype='int32')

    else:
        raise NotImplementedError('animate_scale=%s animate_phase=%s animate_time=%s' % (
            animate_scale, animate_phase, animate_time))
    if istep is not None:
        assert isinstance(istep, integer_types), 'istep=%r' % istep
        scales = (scales[istep],)
        phases = (phases[istep],)
        isteps = (istep,)

    phases2, icases2, isteps2, scales2 = update_animation_inputs(
        phases, icases, isteps, scales,
        analysis_time, fps)
    return phases2, icases2, isteps2, scales2, analysis_time

def get_analysis_time(time, onesided=True):
    """
    The analysis time is the time that needs to be simulated for the analysis.

    TODO: could we define time as 1/2-sided time so we can do less work?
    TODO: we could be more accurate regarding dt
          Nonesided = 5
          Ntwosided = 2 * Nonesided - 1 = 9
          Nonesided = (Ntwosided + 1) / 2

          Nframes = int(fps * t)
          Nonesided = Nframes
          Ntwosided = 2 * Nonesided - 1 = 9
          Nonesided = (Ntwosided + 1) / 2
    """
    if onesided:
        analysis_time = time / 2.
    else:
        analysis_time = time
    return analysis_time

def update_animation_inputs(phases, icases, isteps, scales, analysis_time, fps):
    """
    Simplifies the format of phases, icases, steps, scales to make them
    into lists of the correct length.

    Parameters
    ----------
    phases : List[float] or None
        List[float] : the phase angles
        None : real result (same as [0., 0., ...])
    icases : List[int] or int
        List[int] : the icases to run
        int : single icase (e.g., SOL 101, 103, 145; same as [icase, icase, ...]
    isteps : List[int]
        nominal isteps = [0, 1, 2, 3, 4, ..., nframes]
        we can analyze pictures [1, 3, 4] by providing a subset
    scales : List[float]
        the displacement scale factor; true scale
    """
    if phases is not None:
        pass
    elif phases is None:
        phases = [0.] * len(scales)
    else:
        raise RuntimeError('phases=%r' % phases)

    if isinstance(icases, integer_types):
        icases = [icases] * len(scales)


    if len(icases) != len(scales):
        msg = 'ncases=%s nscales=%s' % (len(icases), len(scales))
        print(msg)
        raise ValueError(msg)

    if len(icases) != len(phases):
        msg = 'ncases=%s nphases=%s' % (len(icases), len(phases))
        print(msg)
        raise ValueError(msg)

    if isteps is None:
        isteps = np.linspace(0, len(scales), endpoint=False, dtype='int32')
        print("setting isteps in make_gif")

    if len(scales) != len(isteps):
        msg = 'len(scales)=%s len(isteps)=%s analysis_time=%s fps=%s' % (
            len(scales), len(isteps), analysis_time, fps)
        print(msg)
        raise ValueError(msg)
    assert isinstance(isteps[0], integer_types), 'isteps=%s, must be integers' % isteps
    return phases, icases, isteps, scales

def write_gif(gif_filename, png_filenames, time=2.0,
              onesided=True, nrepeat=0,
              delete_images=False, make_gif=True):
    """
    Makes an animated gif

    Parameters
    ----------
    gif_filename : str
        path to the output gif & png folder
    png_filenames : List[str]
        the pictures to make the gif from
    time : float; default=2.0
        the runtime of the gif (seconds)

    Options
    -------
    onesided : bool; default=True
        should the animation go up and back down
    nrepeat : int; default=0
        0 : loop infinitely
        1 : loop 1 time
        2 : loop 2 times

    Final Control Options
    ---------------------
    delete_images : bool; default=False
        cleanup the png files at the end
    make_gif : bool; default=True
        actually make the gif at the end

    Other local variables
    ---------------------
    duration : float
        frame time (seconds)
    """
    if not IS_IMAGEIO:
        return False

    #assert fps >= 1, fps
    nframes = len(png_filenames)
    duration = time / nframes

    gif_dirname = os.path.dirname(os.path.abspath(gif_filename))
    if not os.path.exists(gif_dirname):
        os.makedirs(gif_dirname)

    if not onesided:
        # drop the duplicate middle frame
        # >>> a = [1, 2, 3, 4, 5]
        # >>> a + a[-2::-1]
        # [1, 2, 3, 4, 5, 4, 3, 2, 1]
        png_filenames = png_filenames + png_filenames[-2::-1]

    if make_gif and IS_IMAGEIO:
        images = []
        for png_filename in png_filenames:
            images.append(imageio.imread(png_filename))
        if nrepeat is True:
            nrepeat = 0
        imageio.mimsave(gif_filename, images, duration=duration,
                        loop=nrepeat)

    if delete_images:
        for png_filename in png_filenames:
            try:
                os.remove(png_filename)
            except OSError:
                pass
    return True
