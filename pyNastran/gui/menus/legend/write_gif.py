"""
Defines:
 - write_gif(gif_filename, png_filenames, time=2.0, fps=30,
              onesided=True, nrepeat=0,
              delete_images=False, make_gif=True)

"""
import os
from typing import Optional

import numpy as np
from pyNastran.utils import int_version
try:
    import imageio
    IMAGEIO_VERSION = int_version('imageio', imageio.__version__)
    if IMAGEIO_VERSION >= [2, 16, 2]:
        import imageio.v2 as imageio
    import PIL
    IS_IMAGEIO = True
except ModuleNotFoundError:
    IS_IMAGEIO = False

from pyNastran import is_pynastrangui_exe
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils import remove_files, PathLike

if is_pynastrangui_exe:  # pragma: no cover
    assert IS_IMAGEIO, 'imageio is not installed for the exe'

if IS_IMAGEIO:
    from pyNastran.utils import int_version
    iver = int_version('pillow', PIL.__version__)
    assert iver != [7, 1, 0], 'pillow=7.1.0 is not supported'


def setup_animation(scale, istep=None,
                    animate_scale: bool=True,
                    animate_phase: bool=False,
                    animate_time: bool=False,
                    icase_fringe=None, icase_disp=None, icase_vector=None,

                    icase_fringe_start=None, icase_fringe_end=None, icase_fringe_delta=None,
                    icase_disp_start=None, icase_disp_end=None, icase_disp_delta=None,
                    icase_vector_start=None, icase_vector_end=None, icase_vector_delta=None,

                    time: float=2.0, animation_profile: str='0 to scale',
                    fps: int=30, animate_in_gui: bool=False):
    """
    helper method for ``make_gif``

    Returns
    -------
    phases : list[float]
        list[float] : the phase angles
    icases_fringe/disp/vector : list[int]
        list[int] : the icases to run
    isteps : list[int]
        nominal isteps = [0, 1, 2, 3, 4, ..., nframes]
        we can analyze pictures [1, 3, 4] by providing a subset
    scales : list[float]
        the displacement scale factor; true scale
    analysis_time : float
        the time that needs to be simulated for the analysis; not the runtime

    """
    if animate_scale or animate_phase:
        # ignored for time
        assert isinstance(fps, integer_types), f'fps={fps:d} must be an integer'

    phases = None
    onesided = False
    is_symmetric = False
    endpoint = False
    if animate_scale:
        out = setup_animate_scale(
            scale,
            icase_fringe, icase_disp, icase_vector,
            time, animation_profile, fps)
        icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, is_symmetric, endpoint = out
    elif animate_phase:
        onesided = True
        endpoint = True
        is_symmetric = False
        out = setup_animate_phase(
            scale,
            icase_fringe, icase_disp, icase_vector,
            time, fps)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time = out
    elif animate_time:
        onesided = True
        is_symmetric = False
        endpoint = True
        out = setup_animate_time(
            scale, time, fps,
            icase_fringe_start, icase_fringe_end, icase_fringe_delta,
            icase_disp_start, icase_disp_end, icase_disp_delta,
            icase_vector_start, icase_vector_end, icase_vector_delta)
        icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, fps = out
    else:
        raise NotImplementedError('animate_scale=%s animate_phase=%s animate_time=%s' % (
            animate_scale, animate_phase, animate_time))

    out = update_animation_inputs(
        phases,
        icases_fringe, icases_disp, icases_vector,
        isteps, scales,
        analysis_time, fps)
    phases2, icases_fringe2, icases_disp2, icases_vector2, isteps2, scales2 = out

    if istep is None and animate_scale: # and not animate_in_gui:
        scales2, phases2, icases_fringe2, icases_disp2, icases_vector2, isteps2 = make_symmetric(
            scales2, phases2, icases_fringe2, icases_disp2, icases_vector2, isteps2, is_symmetric, endpoint)
        #if animate_in_gui:
            # double the number of frames
            # drop the duplicate end frame if necessary

    # grab the ith step
    if istep is not None:
        assert isinstance(istep, integer_types), f'istep={istep!r}'
        scales = (scales2[istep],)
        phases = (phases2[istep],)
        isteps = (istep,)
    #print('scales_final=%s' % scales)
    return phases2, icases_fringe2, icases_disp2, icases_vector2, isteps2, scales2, analysis_time, onesided, endpoint


def fix_nframes(nframes: int, profile: str) -> int:
    """
    make sure we break at the "true scale" max
    should this be in terms of fps or nframes? (I think nframes)
    """
    DIV2_PROFILES = {
        '0 to scale to 0',
        '-scale to scale to -scale',
    }
    DIV4_PROFILES = {
        '0 to scale to -scale to 0',
        'sinusoidal: 0 to scale to -scale to 0',
        'sinusoidal: scale to -scale to scale',
    }

    fix_nframes_even = nframes == 1 or nframes % 2 == 0
    fix_nframes_sin = nframes == 1 or nframes % 4 != 1
    is_div_two_profile = profile in DIV2_PROFILES
    is_div_four_profile = profile in DIV4_PROFILES
    if is_div_two_profile and fix_nframes_even:
        nframes_div_2 = nframes // 2
        nframes = 2 * (nframes_div_2 + 1) + 1
    if is_div_four_profile and fix_nframes_sin:
        nframes_div_4 = nframes // 4
        nframes = 4 * (nframes_div_4 + 1) + 1
    return nframes

def _get_scale_profile_info(profile: str) -> tuple[str, bool, bool, bool]:
    if isinstance(profile, str):
        profile = profile.lower()
        if profile == '0 to scale':
            # doesn't end at start point
            # not symmetric
            onesided = False
            endpoint = False
            is_symmetric = False
        elif profile == '0 to scale to 0':
            # ends at start point
            # symmetric
            onesided = True
            endpoint = True
            is_symmetric = True
        elif profile == '-scale to scale':
            # not symmetric
            onesided = False
            endpoint = False
            is_symmetric = False
        elif profile == '-scale to scale to -scale':
            # ends at start point
            # symmetric
            onesided = True
            endpoint = True
            is_symmetric = True
        elif profile in ['0 to scale to -scale to 0',
                         'sinusoidal: 0 to scale to -scale to 0']:
            # not symmetric
            onesided = False
            endpoint = True
            is_symmetric = False
        elif profile == 'sinusoidal: scale to -scale to scale':
            # symmetric
            #onesided = True
            #endpoint = False
            #is_symmetric = True
            onesided = False
            endpoint = True
            is_symmetric = False
        else:
            msg = (
                f'profile={profile!r} is not supported:\n'
                "  '0 to scale'\n"
                "   '0 to scale to 0'\n"
                "  '-scale to scale'\n"
                "  '-scale to scale to -scale'\n"
                "  '0 to scale to -scale to 0'\n"
                "  'sinusoidal: 0 to scale to -scale to 0'\n"
                "  'sinusoidal: scale to -scale to scale'\n")
            raise NotImplementedError(msg.rstrip())
    else:
        msg = f'profile={profile!r} is not supported'
        raise NotImplementedError(msg)
    return profile, onesided, endpoint, is_symmetric

def setup_animate_scale(scale: float,
                        icase_fringe: int,
                        icase_disp: int,
                        icase_vector: Optional[int],
                        time: float, profile: str, fps: int):
    """
    Gets the inputs for a displacement scale/real modal animation

    A onesided animation is an animation that:
     - does not need to loop back on itself because it ends back at the
       start point
     - is symmetric

    We want to set it to true because if an animation is onesided, we can skip
    making half the images.
    """
    profile, onesided, endpoint, is_symmetric = _get_scale_profile_info(profile)

    analysis_time = get_analysis_time(time, onesided)
    #fps = fix_nframes(fps, profile)
    nframes = int(analysis_time * fps)

    nframes = fix_nframes(nframes, profile)
    #print('analysis_time*fps=%s*%s=%s nframes=%s profile=%r' % (
        #analysis_time, fps, analysis_time * fps, nframes, profile))

    # endpoint=False if scale[0] == scale[-1]; otherwise True
    nframes_interp = nframes - 1

    # TODO: add more explicit flag for endpoint?
    is_linspace_endpoint = not endpoint
    isteps = np.linspace(0, nframes, num=nframes, endpoint=is_linspace_endpoint, dtype='int32')

    if isinstance(profile, str):
        profile = profile.lower()
        if profile == '0 to scale':
            yp = np.array([0, scale])
            xp = np.array([0., nframes_interp])
            scales = np.interp(isteps, xp, yp)

        elif profile == '-scale to scale':
            yp = np.array([-scale, scale])
            xp = np.array([0., nframes_interp])
            scales = np.interp(isteps, xp, yp)

        elif profile == '0 to scale to 0':
            yp = np.array([0, scale, 0.])
            xp = np.array([0., nframes_interp / 2., nframes_interp])
            scales = np.interp(isteps, xp, yp)

        elif profile == '-scale to scale to -scale':
            yp = np.array([-scale, scale, -scale])
            xp = np.array([0., nframes_interp / 2., nframes_interp])
            scales = np.interp(isteps, xp, yp)

        elif profile == '0 to scale to -scale to 0':
            yp = np.array([0., 1., 0., -1., 0.])
            xp = np.array([0., nframes_interp * 0.25, nframes_interp * 0.5,
                           nframes_interp * 0.75, nframes_interp])
            scales = scale * np.interp(isteps, xp, yp)
            #print('x=%s xp=%s yp=%s' % (x, xp, yp))
            #print('scales=%s' % (scales))
            #print('isteps* = ', isteps)

        elif profile == 'sinusoidal: 0 to scale to -scale to 0':
            theta = np.pi * np.array([0., 0.5, 1., 1.5, 2.])
            xp = np.array([0., nframes_interp * 0.25, nframes_interp * 0.5,
                           nframes_interp * 0.75, nframes_interp])
            scales = scale * np.sin(np.interp(isteps, xp, theta))
        elif profile == 'sinusoidal: scale to -scale to scale':
            theta = np.pi * np.array([0., 1., 2.])
            xp = np.array([0., nframes_interp / 2., nframes_interp])
            scales = scale * np.cos(np.interp(isteps, xp, theta))

        else:
            msg = (
                f'profile={profile!r} is not:\n'
                "  '0 to scale'\n"
                "  '0 to scale to 0'\n"
                "  '-scale to scale'\n"
                "  '-scale to scale to -scale'\n"
                "  '0 to scale to -scale to 0'\n"
                "  'sinusoidal: 0 to scale to -scale to 0'\n"
                "  'sinusoidal: scale to -scale to scale'\n")
            raise NotImplementedError(msg.rstrip())
    #elif isinstance(profile, list):
        #yp = np.array(profile)
        #xp = np.linspace(0., nframes_interp, num=len(yp), endpoint=True, dtype='float64')
    else:
        msg = f'profile={profile!r} is not supported'
        raise NotImplementedError(msg)

    icases_fringe = icase_fringe
    icases_disp = icase_disp
    icases_vector = icase_vector
    assert len(scales) == len(isteps), f'nscales={len(scales):d} nsteps={len(isteps):d}'
    #assert len(scales) == nframes, 'len(scales)=%s nframes=%s' % (len(scales), nframes)

    # TODO: this can hit
    #if profile == '0 to scale to 0' and len(scales) % 2 == 0:
        #raise RuntimeError('nscales=%s scales=%s' % (len(scales), scales))
    out = (
        icases_fringe, icases_disp, icases_vector,
        isteps, scales, analysis_time, onesided, is_symmetric, endpoint,
    )
    return out

def setup_animate_phase(scale: float,
                        icase_fringe: int, icase_disp: int, icase_vector: int,
                        time: float, fps: int):
    """Gets the inputs for a phase animation"""
    nframes = int(time * fps)
    icases_fringe = icase_fringe
    icases_disp = icase_disp
    icases_vector = icase_vector
    phases = np.linspace(0., 360., num=nframes, endpoint=False)
    isteps = np.linspace(0, nframes, num=nframes, endpoint=False, dtype='int32')
    scales = [scale] * len(isteps)
    assert len(phases) == len(isteps), f'nphases={len(phases):d} nsteps={len(isteps):d}'
    assert len(scales) == len(isteps), f'nscales={len(scales):d} nsteps={len(isteps):d}'
    #assert len(phases) == nframes, 'len(phases)=%s nframes=%s' % (len(phases), nframes)

    icases_fringe = icase_fringe
    icases_disp = icase_disp
    icases_vector = icase_vector
    return phases, icases_fringe, icases_disp, icases_vector, isteps, scales, time

def setup_animate_time(scale: float, time: float, fps: int,
                       icase_fringe_start=None, icase_fringe_end=None, icase_fringe_delta=None,
                       icase_disp_start=None, icase_disp_end=None, icase_disp_delta=None,
                       icase_vector_start=None, icase_vector_end=None, icase_vector_delta=None):
    """Gets the inputs for a transient animation"""
    analysis_time = time
    icases_tuple = [
        ('fringe', icase_fringe_start, icase_fringe_end, icase_fringe_delta),
        ('disp', icase_disp_start, icase_disp_end, icase_disp_delta),
        ('vector', icase_vector_start, icase_vector_end, icase_vector_delta),
    ]
    icases = {}
    for (name, icase_start, icase_end, icase_delta) in icases_tuple:
        if icase_start is None:
            continue
        assert isinstance(icase_start, integer_types), f'icase_start={icase_start:d}; type={type(icase_start)}'
        assert isinstance(icase_end, integer_types), f'icase_end={icase_end:d}; type={type(icase_end)}'
        assert isinstance(icase_delta, integer_types), f'icase_delta={icase_delta:d}; type={type(icase_delta)}'
        icases[name] = np.arange(icase_start, icase_end+1, icase_delta)

    icases_disp = icases.get('disp', [])
    icases_fringe = icases.get('fringe', [])
    icases_vector = icases.get('vector', [])

    # a fringe is always available
    if icase_fringe_start is None:
        if icase_disp_start is not None:
            icases_fringe = icases_disp
        elif icase_vector_start is not None:
            icases_fringe = icases_vector
    if len(icases_disp) == 0:
        icases_disp = None
    if len(icases_vector) == 0:
        icases_vector = None

    #min_value = 0.
    #max_value = 1.46862
    nfiles = len(icases_fringe)

    # specifying fps and dt makes the problem overdefined
    # assuming dt
    #
    # TDOO: change this to stepping similar to icase_delta
    #       icases = icases[::5]
    fps = nfiles / time

    # our scale will be constant
    # phases is just None b/c time is real
    scales = [scale] * nfiles
    if len(icases_fringe) != nfiles:
        msg = f'len(icases)={len(icases_fringe):d} nfiles={nfiles}'
        raise ValueError(msg)
    if len(scales) != len(icases_fringe):
        msg = f'nscales={len(scales):d} len(icases_fringe)={len(icases_fringe):d}'
        raise ValueError(msg)

    # TODO: this isn't maintained...
    #assert nframes == nfiles, 'nframes=%s nfiles=%s' % (nframes, nfiles)

    isteps = np.linspace(0, nfiles, num=nfiles, endpoint=True, dtype='int32')
    return icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, fps

def get_analysis_time(time: float, onesided: bool=True) -> float:
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

def update_animation_inputs(phases, icases_fringe, icases_disp, icases_vector,
                            isteps: list[int], scales: list[float],
                            analysis_time: float, fps: int):
    """
    Simplifies the format of phases, icases, steps, scales to make them
    into ndarrays of the correct length.

    Parameters
    ----------
    phases : list[float] or None
        list[float] : the phase angles
        None : real result (same as [0., 0., ...])
    icases_fringe/disp/vector : list[int] or int
        list[int] : the icases to run
        int : single icase (e.g., SOL 101, 103, 145; same as [icase, icase, ...]
    isteps : list[int]
        nominal isteps = [0, 1, 2, 3, 4, ..., nframes]
        we can analyze pictures [1, 3, 4] by providing a subset
    scales : list[float]
        the displacement scale factor; true scale
    analysis_time: float
        used for logging
    fps : int
        used for logging

    """
    if phases is not None:
        pass
    elif phases is None:
        phases = [0.] * len(scales)
    else:
        raise RuntimeError(f'phases={phases!r}')

    if icases_fringe is None or isinstance(icases_fringe, integer_types):
        icases_fringe = [icases_fringe] * len(scales)
    if icases_disp is None or isinstance(icases_disp, integer_types):
        icases_disp = [icases_disp] * len(scales)
    if icases_vector is None or isinstance(icases_vector, integer_types):
        icases_vector = [icases_vector] * len(scales)

    if not(len(icases_disp) or len(icases_fringe) or len(icases_vector)):
        raise RuntimeError('No displacement, fringe, or vector cases have been defined')

    assert icases_fringe is not None
    if len(icases_fringe) != len(scales):
        raise ValueError(f'ncases_fringe={len(icases_fringe):d} nscales={len(scales):d}')

    assert icases_disp is not None
    if len(icases_disp) != len(scales):
        raise ValueError(f'ncases_disp={len(icases_disp)} nscales={len(scales)}')

    assert icases_vector is not None
    if len(icases_vector) != len(scales):
        raise ValueError(f'ncases_vector={len(icases_vector):d} nscales={len(scales):d}')

    if len(icases_fringe) != len(phases):
        raise ValueError(f'ncases_fringe={len(icases_fringe):d} nphases={len(phases):d}')

    if isteps is None:
        isteps = np.linspace(0, len(scales), endpoint=False, dtype='int32')
        #print("setting isteps in make_gif")

    if len(scales) != len(isteps):
        msg = 'len(scales)=%s len(isteps)=%s analysis_time=%s fps=%s' % (
            len(scales), len(isteps), analysis_time, fps)
        print(msg)
        raise ValueError(msg)
    #print('scales=%s' % scales)
    assert isinstance(isteps[0], integer_types), f'isteps={isteps}, must be integers'

    phases2 = np.array(phases)
    icases_fringe2 = np.array(icases_fringe)
    icases_disp2 = np.array(icases_disp)
    icases_vector2 = np.array(icases_vector)
    isteps2 = np.array(isteps)
    scales2 = np.array(scales)
    return phases2, icases_fringe2, icases_disp2, icases_vector2, isteps2, scales2

def make_symmetric(scales, phases, icases_fringe, icases_disp, icases_vector,
                   isteps, is_symmetric: bool, endpoint: bool):
    """
    Chop the frames in half at the middle frame

    Examples
    --------
    >>> isteps = [1, 2, 3, 4, 5, 4, 3, 2, 1]
    >>> make_symmetric(scales, phases, icases_fringe, icases_disp, icases_vector,
                       isteps, is_symmetric)
    isteps = [1, 2, 3, 4, 5]

    Examples
    --------
    >>> scales = [0., 1., 0.]
    >>> isteps = [1, 2, 3]
    >>> make_symmetric(scales, phases, icases_fringe, icases_disp, icases_vector,
                       isteps, is_symmetric)
    scales = [0., 1.]
    isteps = [1, 2]
    """
    if not is_symmetric:
        if endpoint:
            out = (
                scales[:-1], phases[:-1],
                icases_fringe[:-1], icases_disp[:-1], icases_vector[:-1],
                isteps[:-1]
            )
        else:
            out = scales, phases, icases_fringe, icases_disp, icases_vector, isteps
        return out

    # if twosided
    i = None
    if is_symmetric:
        nframes = len(scales)
        ihalf_frame = nframes // 2 + 1
        assert nframes % 2 == 1, nframes
        i = ihalf_frame
    elif endpoint:
        i = -1

    if i is not None:
        if isinstance(phases, np.ndarray):
            phases = phases[:i]
        if isinstance(icases_fringe, np.ndarray):
            icases_fringe = icases_fringe[:i]
        if isinstance(icases_disp, np.ndarray):
            icases_disp = icases_disp[:i]
        if isinstance(icases_vector, np.ndarray):
            icases_vector = icases_vector[:i]
        isteps = isteps[:i]
        scales = scales[:i]
    return scales, phases, icases_fringe, icases_disp, icases_vector, isteps

def make_two_sided(scales, phases, icases_fringe, icases_disp, icases_vector, isteps,
                   onesided: bool):
    """
    Drop the duplicate middle frame if we're onesided.  Only for the GUI.

    Examples
    --------
    this is a onesided "mountain", we only go up
    >>> a = [1, 2, 3, 4, 5]
    >>> a + a[-2::-1]

    this is a two sided "mountain"; we have to go back down
    [1, 2, 3, 4, 5, 4, 3, 2, 1]
    """
    is_endpoint = np.allclose(scales[0], scales[-1])
    end = None

    if is_endpoint:
        end = -1

    # if twosided
    #if not onesided:
        #print('A-end=%s scales=%s, isteps=%s' % (end, scales, isteps))
        #print(scales[:end], scales[-1::-1])
        #phases = np.hstack([phases[:end], phases[-1::-1]])
        #icases_fringe = np.hstack([icases_fringe[:end], icases_fringe[-1::-1]])
        #icases_disp = np.hstack([icases_disp[:end], icases_disp[-1::-1]])
        #icases_vector = np.hstack([icases_vector[:end], icases_vector[-1::-1]])
        #isteps = np.hstack([isteps[:end], isteps[-1::-1]])
        #scales = np.hstack([scales[:end], scales[-1::-1]])
        #print('B-scales=%s n=%s, isteps=%s' % (scales, len(scales), isteps))
        ##analysis_time *= 2
        #end = None

    # make sure that we don't screw up onesided by putting this on repeat
    #print('scales[:end]=%s' % scales[:end])
    phases = np.hstack([phases[:end], phases[:end]])
    icases_fringe = np.hstack([icases_fringe[:end], icases_fringe[:end]])
    icases_disp = np.hstack([icases_disp[:end], icases_disp[:end]])
    icases_vector = np.hstack([icases_vector[:end], icases_vector[:end]])
    isteps = np.hstack([isteps[:end], isteps[:end]])
    scales = np.hstack([scales[:end], scales[:end]])
    #print('scales2     =%s n=%s, isteps2=%s' % (scales, len(scales), isteps))
    return scales, phases, icases_fringe, icases_disp, icases_vector, isteps

def write_gif(gif_filename: PathLike, png_filenames: list[PathLike], time: float=2.0,
              onesided: bool=True, nrepeat: int=0,
              delete_images: bool=False, make_gif: bool=True) -> bool:
    """
    Makes an animated gif

    Parameters
    ----------
    gif_filename : str
        path to the output gif & png folder
    png_filenames : list[str]
        the pictures to make the gif from
    time : float; default=2.0
        the runtime of the gif (seconds)
    onesided : bool; default=True
        should the animation go up and back down
        True : the video will use images [0...N]
        False : the video will use images [0...N...0]
    nrepeat : int; default=0
        0 : loop infinitely
        1 : loop 1 time
        2 : loop 2 times
    delete_images : bool; default=False
        cleanup the png files at the end
    make_gif : bool; default=True
        actually make the gif at the end

    Returns
    -------
    success : bool
        was the gif made

    """
    if not IS_IMAGEIO:
        return False

    #assert fps >= 1, fps
    nframes = len(png_filenames)
    assert nframes > 0, png_filenames

    # duration : float
    # frame time (seconds)
    duration = time / nframes

    gif_dirname = os.path.dirname(os.path.abspath(gif_filename))
    if not os.path.exists(gif_dirname):
        os.makedirs(gif_dirname)

    #if not onesided:
        # drop the duplicate middle frame
        # >>> a = [1, 2, 3, 4, 5]
        # >>> a + a[-2::-1]
        # [1, 2, 3, 4, 5, 4, 3, 2, 1]
        #png_filenames = png_filenames + png_filenames[-2::-1]

    if make_gif and IS_IMAGEIO:
        images = []
        for png_filename in png_filenames:
            if not isinstance(png_filename, str) and os.path.exists(png_filename):
                raise TypeError(f'png_filename={png_filename!r} is invalid')
            imagei = imageio.imread(png_filename)
            images.append(imagei)
        if nrepeat is True:
            nrepeat = 0
        try:
            imageio.mimsave(gif_filename, images, duration=duration,
                            loop=nrepeat)
        except IOError:  # file is open
            raise IOError(f'{gif_filename} is likely open')

    if delete_images:
        remove_files(png_filenames)
    return True
