#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2007


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on pyaudiolab.
# """


# seem to need to have these imports here. . . to make sure names are defined
# from numpy import *
# from pyaudiolab import *
from numpy import *
from muse import *
from filters import *
from transforms import *
from encoders import *
from decoders import *

# import muse defined constants
import constants as C


def vrms(x, t, T, zi = None):
    """vrms(x, t, T, zi = None)
    
    Args:
      x -- Input signal
      t -- RMS averaging time, in sec.
      T -- Sampling period, in sec/sample.

      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.

    Outputs: (y, {zf})
    
      y -- The output, an rms tracked envelope.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            filter state.
    
    Track the envelope of a signal using the RMS time averaging method.
    """
    # Time varying RMS envelope tracking
    # see Zolzer, DAFX, p 98
    
    TAV = 1 - numpy.exp(-2.2 * T / t)
    b = array([TAV])
    a = array([1.,  (TAV - 1.)])

    res = ffilter(b, a, numpy.square(x), zi)

    return res


def vpeak(a, N, beta = 5, mode = 'z', kind = 'fft', zi = None):
    """vpeak(a, N, beta = 5, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - a    -- Input signal
        - N    -- filter (number of taps), should be odd
        - beta -- beta for Kaiser window FIR design.
                  5 = similiar to Hamming.
        - mode -- 'z' or 'full'. If mode is 'z', acts as a filter
                  with state 'z', and returns a vector of length
                  len(x). If mode is 'full', returns the full
                  convolution.
        - kind -- 'direct' or 'fft', for direct or fft convolution

        - zi   -- Initial state. An array of shape (len(kernel) - 1, 3).
                  If zi=None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y -- The output of the decoder.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            filter state.
    
    Track the envelope of a signal using the instantaneous method.
    """


    # convolve with hilbert kernel
    if zi is not None:
        hb, zf = convfilt(
            a,
            fir_hb(N, beta=5),      # generate hilbert kernel
            mode,
            kind,
            zi
            )
        over_dub(hb, zi, write_over = True)
    else:
        hb = convfilt(
            a,
            fir_hb(N, beta=5),      # generate hilbert kernel
            mode,
            kind
            )

    res = numpy.sqrt(numpy.square(hb.real) + numpy.square(hb.imag))

    if zi is not None:
        return res, zf
    else:
        return res


def venv(a, N, t, T, width=pi, mode = 'z', kind = 'fft', zi = None):
    """venv(a, N, t, T, width=pi, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - a    -- Input signal
        - N    -- filter (number of taps), should be odd
        - t    -- envelope attack/release (averaging) time, in sec.
        - T    -- Sampling period, in sec/sample.
        - width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

        - mode -- 'z' or 'full'. If mode is 'z', acts as a filter
                  with state 'z', and returns a vector of length
                  len(x). If mode is 'full', returns the full
                  convolution.
        - kind -- 'direct' or 'fft', for direct or fft convolution

        - zi   -- Initial state. An array of shape (len(kernel) - 1, 3).
                  If zi=None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y -- The output, tracked envelope.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            filter state.
    
    Track the envelope of a signal abs-->lp method.
    """

    Wn = freq_to_Wn(1./t, T)

    x = abs(a)

    # convolve with low-pass kernel
    if zi is not None:
        lp, zf = convfilt(
            x,
            fir_lp(N, Wn, width),      # generate lp kernel
            mode,
            kind,
            zi
            )
        over_dub(lp, zi, write_over = True)
    else:
        lp = convfilt(
            x,
            fir_lp(N, Wn, width),      # generate lp kernel
            mode,
            kind
            )

    res = lp

    if zi is not None:
        return res, zf
    else:
        return res


# NOTE: vaed() appears to return incorrect azimuth angles, 
#       at least for single sample scalar input.
#       Core testing is needed! Compare with aed()
# vaed (azimuth, elevation, directivity)
def vaed(a, zi = None):
    """vaed(a, zi = None)
    
    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity.
    
    Inputs:
        - b         : Input b-format signal
      zi -- Initial conditions for the filter delays.  It is a vector
            of shape (4, 2).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: ([a, e, d])
    
      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.

    """
    # normalise W
    b = copy(a)
    b[:, 0] *= sqrt(2)

    # pv & b**2
    if zi is None:
        pv = integ_filt(
            b * interleave(b[:, 0])
            )
        b_sqrd = integ_filt(
            b**2
            )
    else:
        zf = zeros_like(zi)
        
        pv, zf[:, :1] = integ_filt(
            b * interleave(b[:, 0]),
            zi[:, :1]
            )
        b_sqrd, zf[:, 1:] = integ_filt(
            b**2,
            zi[:, 1:]
            )

    # p**2 and v**2
    p_sqrd = b_sqrd[:, 0]
    v_sqrd = sum(b_sqrd[:, 1:], 1)

    # calculate azimuth, elevation
    a = unwrap(arctan2(pv[:, 2], pv[:, 1]))
    e = unwrap(arctan2(pv[:, 3], sqrt((pv[:, 1])**2 + (pv[:, 2])**2)))

    # calculate directivity
    # pi/2 - 2 * arctan(v/p)
    d = unwrap(pi/2 - 2 * arctan2(sqrt(v_sqrd), sqrt(p_sqrd)))

    res = interleave(array([a, e, d]))

    if zi is None:
        return res
    else:
        return res, zf


# aed (azimuth, elevation, directivity)
# return constants
def aed(a):
    """aed(a)
    
    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity.
    
    Inputs:
        - b         : Input b-format signal

    Outputs: ([a, e, d])
    
      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)

    """

    # normalise W
    b = copy(a)
    b[:, 0] *= sqrt(2)

    # pv & b**2 mean
    pv_mean = mean(b * interleave(b[:, 0]), 0)
    b_sqrd_mean = mean((b**2), 0)

    # p**2 and v**2
    p_sqrd = b_sqrd_mean[0]
    v_sqrd = sum(b_sqrd_mean[1:])

    # calculate azimuth, elevation
    a = arctan2(pv_mean[2], pv_mean[1])
    e = arctan2(pv_mean[3], sqrt((pv_mean[1])**2 + (pv_mean[2])**2))

    # calculate directivity
    # pi/2 - 2 * arctan(v/p)
    d = pi/2 - 2 * arctan2(sqrt(v_sqrd), sqrt(p_sqrd))

    # return [azimuth, elevation, directivity]
    res = array([a, e, d])

    return res
