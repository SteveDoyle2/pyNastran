import sys
from cmath import rect
from numpy import array, linspace, abs, angle, degrees, radians
from scipy import fft, ifft
from matplotlib.pylab import figure, show, grid, xlabel, ylabel

from scipy.interpolate import splrep, splev

if sys.version_info < (3, 0):
    # "fixes" bug where scipy screws up return code handling
    import scipy.weave

def resample_time_signal(dt, npoints, times, amplitudes):
    spline = splrep(times, amplitudes)
    tmax = npoints * dt
    times2 = linspace(0., tmax, npoints)
    amplitudes2 = splev(times2, spline)
    return times2, amplitudes2
    
def resample_frequency_signal(df, npoints, frequencies, magnitudes, phases):
    spline = splrep(frequencies, magnitudes)
    spline2 = splrep(frequencies, phases)
    fmax = npoints * df
    frequencies2 = linspace(0., fmax, npoints)
    magnitudes2 = splev(frequencies2, spline)
    phases2 = splev(frequencies2, spline2)
    return frequencies2, magnitudes2, phases2

def fft_1sided(times, amplitudes):
    """
    assumes times are evenly spaced
    """
    npoints = len(times)
    tmax = times[-1]
    tmin = times[0]
    dt = (tmax-tmin)/npoints
    fmax = 1./dt
    df = 1./tmax
    fmin = fmax - npoints*df
    fNyquist = fmax/2.
    
    n2 = npoints // 2
    (frequencies, magnitudes, phases, fNyquist, nNyquist) = fft_2sided(times, amplitudes)
    return (frequencies[:n2], magnitudes[:n2], phases[:n2], fNyquist, nNyquist)

def fft_2sided(times, amplitudes):
    """
    assumes times are evenly spaced
    """
    npoints = len(times)
    assert npoints == len(amplitudes)
    z = fft(amplitudes, npoints)
    (magnitudes) = abs(z)
    phases = angle(z)

    tmax = times[-1]
    tmin = times[0]
    dt = (tmax-tmin)/(npoints-1)
    fmax = 1./dt
    df = 1./tmax
    print "tmax = ",tmax
    print "dt = ",dt
    print "npoints = ",npoints
    print "df = ",df
    print "fmax = ",fmax
    fmin = fmax - npoints*df
    frequencies = linspace(fmin, fmax, npoints, endpoint=True)
    fNyquist = fmax/2.
    if npoints % 2 == 1: # odd
        nNyquist = 'iLast/2+1'
    else:
        nNyquist = 'iLast/2'
    return (frequencies, magnitudes, phases, fNyquist, nNyquist)

def power_spectral_density(frequencies, magnitudes):
    npoints = len(frequencies)
    assert npoints == len(magnitudes)
    fmax = frequencies[-1]
    fmin = frequencies[0]
    nfreq = len(frequencies)
    df = (fmax-fmin)/float(nfreq)
    psds = magnitudes*magnitudes/df
    return (psds)

def plot_fourier_transform(frequencies, magnitudes, phases):
    npoints = len(frequencies)
    assert npoints == len(magnitudes)
    assert npoints == len(phases)

    fig = plt.figure(1)
    ax = fig.add_subplot(3,1,1)
    ax.semilogy(frequencies, magnitudes)
    ylabel('Magnitude')
    grid(True)

    ax = fig.add_subplot(3,1,2)
    ax.plot(frequencies, degrees(phases))
    ylabel('Phase (degrees)')
    grid(True)

    ax = fig.add_subplot(3,1,3)
    ax.semilogy(frequencies,
            power_spectral_density(frequencies, magnitudes))
    ylabel('PSD (units$^2$/Hz)')
    xlabel('Frequency (Hz)')
    grid(True)
    show()

def plot_times(times1, amps1, times2, amps2):
    fig = figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.plot(times1, amps1, '-r', times2, amps2, '-b')
    ylabel('Amplitude')
    xlabel('Time (sec)')
    grid(True)
    show()

def ifft_2sided(frequencies, magnitudes, phases, nNyquist):
    """
    assumes times are evenly spaced
    """
    npoints = len(frequencies)
    assert npoints == len(magnitudes)
    assert npoints == len(phases)

    rphases = radians(phases)
    z = []
    for mag,phase in zip(magnitudes, rphases):
        z.append(rect(mag, phase))

    (amplitudes) = ifft(z, npoints)

    fmax = frequencies[-1]
    fmin = frequencies[0]
    df = (fmax-fmin)/(npoints-1)
    
    tmax = 1./df
    dt = 1./fmax
    tmin = tmax - npoints*dt
    print "tmin=%s tmax=%s" %(tmin, tmax)
    times = linspace(tmin, tmax, npoints, endpoint=True)
    assert len(times) == len(amplitudes)
    return (times, amplitudes)

def f(t,dt):
    return t/dt

def main():
    t = linspace(0.,10.,101)
    dt = t[1]-t[0]
    y = array([f(ti, dt) for ti in t])
    (frequencies, magnitudes, phases, fNyquist, nNyquist) = fft_2sided(t, y)
    (t2, y2) = ifft_2sided(frequencies, magnitudes, phases, nNyquist)
    #plot_fourier_transform(frequencies, magnitudes, phases)
    plot_times(t, y, t2, y2)

if __name__ == '__main__':
    main()
    