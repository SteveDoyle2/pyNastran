"""
Defines:
 - write_gif(gif_filename, png_filenames, time=2.0, fps=30,
              onesided=True, nrepeat=0,
              delete_images=False, make_gif=True)
"""
import os
try:
    import imageio
    IS_IMAGEIO = True
except ImportError:
    IS_IMAGEIO = False


def write_gif(gif_filename, png_filenames, time=2.0, fps=30,
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
    fps : int; default=30
        the frames/second

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

    assert fps >= 1, fps
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
