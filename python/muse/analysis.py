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


# vaed (azimuth, elevation, directivity)
def vaed(b):
    """vaed(a)
    
    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity. Returns values as time varying.
    
    Inputs:
        - b         : Input b-format signal

    Outputs: ([a, e, d])
    
      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)

    """

    # calculate directivity
    b_sqrd = b**2

    d = 2 * arctan2(
        C.rec_sqrt2 * sqrt(
            (b_sqrd[:, 1] + b_sqrd[:, 2] + b_sqrd[:, 3])
            ),
        sqrt(b_sqrd[:, 0])
        )

    # adjust directivity so that it is idealized
    b_dir = direct(
        b,
        2. * arctan(1 / tan(d / 2.))
        )

    # calculate azimuth, elevation
    b_abs = a_to_b(
        abs(
            b_to_a(b_dir, weight = 'car')
            ),
        weight = 'car'
        )

    # translated to [r, theta, phi]
    spher = cart_to_spher(b_abs[:, 1:])

    # return [theta, phi, d]
    res = hstack(
        (spher[:, 1:],
        interleave(d))
        )

    return res


# aed (azimuth, elevation, directivity)
# return constants
# may want to make more like aed_rms
def aed(b):
    """aed(a)
    
    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity.
    
    Inputs:
        - b         : Input b-format signal

    Outputs: ([a, e, d])
    
      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)

    """

    # calculate directivity
    b_sqd_mean = (b**2).mean(0)

    d = 2 * arctan2(
        C.rec_sqrt2 * sqrt(
            (b_sqd_mean[1] + b_sqd_mean[2] + b_sqd_mean[3])
            ),
        sqrt(b_sqd_mean[0])
        )

    # adjust directivity so that it is idealized
    b_dir = direct(
        b,
        2. * arctan(1 / tan(d / 2.))
        )

    # calculate azimuth, elevation
    b_abs = a_to_b(
        abs(
            b_to_a(b_dir, weight = 'car')
            ),
        weight = 'car'
        )

    # translated to [r, theta, phi]
    spher = mean(
        cart_to_spher(b_abs[:, 1:]),
        0
        )

    # return [theta, phi, d]
    res = append(spher[1:], d)

    return res


# aed (azimuth, elevation)
# return constants
def aed_rms(b):
    """aed_rms(a)
    
    Analyze an ambisonic B-format sound field, returning azimuth,
    elevation and directivity. (Using RMS.)
    
    Inputs:
        - b         : Input b-format signal

    Outputs: ([a, e, d])
    
      [a, e, d] -- Azimuth, elevation and directivity in radians.
                   (See direct for details on directivity.)

    """

    # calculate directivity
    b_sqd_mean = (b**2).mean(0)

    d = 2 * arctan2(
        C.rec_sqrt2 * sqrt(
            (b_sqd_mean[1] + b_sqd_mean[2] + b_sqd_mean[3])
            ),
        sqrt(b_sqd_mean[0])
        )

    # adjust directivity so that it is idealized
    b_dir = direct(
        b,
        2. * arctan(1 / tan(d / 2.))
        )

    # calculate azimuth, elevation
    b_rms = a_to_b(
        rms(
            b_to_a(b_dir, weight = 'car'),
            0
            ),
        weight = 'car'
        )

    spher = cart_to_spher(
        b_rms[1:],
        )

    # return [theta, phi, d]
    res = append(spher[1:], d)

    return res
