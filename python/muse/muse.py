#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2007


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on numpy an pyaudiolab.
# """

from constants import *         # may not need this here
from numpy import *

#=========================
# Functions - value converters
# 
#=========================

def db_to_amp(a):
    """db_to_amp(a)
    
    Return amp from db.

    """

    return pow(10., .05 * asarray(a))


def amp_to_db(a):
    """amp_to_db(a)
    
    Return db from amp.

    """

    return 20. * log10(asarray(abs(a)))


def freq_to_period(a):
    """freq_to_period(a)
    
    Return period (in seconds) from freq (in hz).

    """
    return reciprocal(asarray(a).astype(float))


def period_to_freq(a):
    """period_to_freq(a)
    
    Return freq (in hz) from period (in seconds).

    """
    return reciprocal(asarray(a).astype(float))


def dur_to_nframes(a, sr, frac = False):
    """dur_to_nframes(a, sr, frac = False)
    
    Return nframes from duration (in sec).

    sr is the sampling rate.

    If frac is True, returns fractional value.

    """
    res = sr * asarray(a)

    if frac is False:
        res = around(res).astype(int)
    
    return res


def nframes_to_dur(a, T):
    """nframes_to_dur(a, T)
    
    Return duration (in sec) from nframes .

    T is the sampling period.

    """
    return T * asarray(a)


def freq_to_Wn(a, T):        # use with lfilter, normalized to nyquist
    """freq_to_Wn(a, T)
    
    Return normalized frequency from freq (in hz).

    T is the sampling period.

    """
    return 2. * T * asarray(a)


def Wn_to_freq(a, sr):       # use with lfilter, normalized to nyquist
    """Wn_to_freq(a, sr)
    
    Return freq (in hz) from normalized frequency.

    sr is the sampling rate.

    """
    return .5 * sr * asarray(a)


def freq_to_oct(a):
    """freq_to_oct(a)
    
    Return an octave-point-decimal value from freq (in hz).

    C4 = 8.00
    """
    A4 = 440.
    freq_ref = A4/pow(2, 8.75)
    
    return log2(asarray(a) / freq_ref)


def oct_to_freq(a):
    """oct_to_freq(a)
    
    Return freq (in hz) from an octave-point-decimal value.

    C4 = 8.00
    """
    A4 = 440.
    freq_ref = A4/pow(2, 8.75)
    
    return freq_ref * pow(2, asarray(a))


def pch_to_oct(a):
    """pch_to_oct(a)
    
    Return an octave-point-decimal value from octave-point-pitch.

    C4 = 8.00
    """
    octave = floor(a)
    decimal  = (a - octave) / .12
    
    return add(octave, decimal)


def oct_to_pch(a):
    """oct_to_pch(a)
    
    Return an octave-point-pitch value from octave-point-decimal.

    C4 = 8.00
    """
    octave = floor(a)
    pitch  = (a - octave) * .12
    
    return add(octave, pitch)


def freq_to_pch(a):
    """freq_to_pch(a)
    
    Return an octave-point-pitch value from freq (in hz).

    C4 = 8.00
    """
    return oct_to_pch(freq_to_oct(a))


def pch_to_freq(a):
    """pch_to_freq(a)
    
    Return freq (in hz) from an octave-point-pitch value.

    C4 = 8.00
    """
    return oct_to_freq(pch_to_oct(a))


def Wn_warp(Wn, c, Wn_w = False):
    """Wn_warp(Wn, c, Wn_w = False)
    
    Bilinear warping. Return warped Wn.

    Warping may be specified by c, (allpass) warping
    coefficient OR by Wn_w, the warped value of input Wn = .5.

    c = 0, no warping
    Wn_w = .5, no warping

    Note, c is -1 * value used by lagt.

    """

    if Wn_w:
        # this one here matches behaviour w/ calculating c and Wn
#         Wn1 = arctan2(
#             sin(pi * c) * sin(pi * Wn),
#             -cos(pi * c) + cos(pi * Wn)
#             ) / pi
        # this one acts as expecting w/ mapping
        # there may be a necessary modification of signs
        # w/ AP filtering?
        # ALSO NOTE: win warping is now dependent on this function
        Wn1 = arctan2(
            sin(pi * c) * sin(pi * Wn),
            cos(pi * c) + cos(pi * Wn)
            ) / pi
    else:
        Wn1 = arctan2(
            ((1 - c**2) * sin(pi * Wn)),
            (-2 * c + (1 + c**2) * cos(pi * Wn))
            ) / pi

    return Wn1


def deg_to_rad(a):
    """deg_to_rad(a)
    
    Return radians from degrees value.

    """
    return a * pi / 180.


def rad_to_deg(a):
    """rad_to_deg(a)
    
    Return degrees from radians value.

    """
    return a * 180. / pi


def pol_to_cart(a):
    """pol_to_cart(a)
    
    Return Cartesian coordinates from polar.
    Angles in radians.

    [ r, theta ] ---> [ x, y ]

    """
    if a.ndim > 2:
        raise "rank > 2 not supported"
    else:
        if a.ndim is 2:
            a = deinterleave(a)

    cosTheta, sinTheta = cos(a[1]), sin(a[1])

    res = array([
            a[0] * cosTheta,
            a[0] * sinTheta,
            ])

    if a.ndim is 2:
        return interleave(res)
    else:
        return res


def cart_to_pol(a):
    """cart_to_pol(a)
    
    Return polar coordinates from Cartesian.
    Angles in radians.

    [ x, y ] ---> [ r, theta ]

    """
    if a.ndim > 2:
        raise "rank > 2 not supported"
    else:
        if a.ndim is 2:
            a = deinterleave(a)

    r = sqrt(sum(asarray(a)**2, axis = 0))

    if r is 0.:
        theta = 0.
    elif r is zeros_like(r):
        theta = zeros_like(r)
    else:
#         theta = arctan2( a[1], a[0] )
        theta = nan_to_num(arctan2( a[1], a[0] ))

#     res = array([ r, theta ])

    if a.ndim is 2:
#         return interleave(res)
        return interleave(array([ r, unwrap(theta) ]))
    else:
        return array([ r, theta ])


def spher_to_cart(a):
    """spher_to_cart(a)
    
    Return Cartesian coordinates from spherical.
    Angles in radians.

    [ r, theta, phi ] ---> [ x, y, z ]

    """
    if a.ndim > 2:
        raise "rank > 2 not supported"
    else:
        if a.ndim is 2:
            a = deinterleave(a)

    cosTheta, sinTheta = cos(a[1]), sin(a[1])
    cosPhi, sinPhi = cos(a[2]), sin(a[2])

    res = array([
            a[0] * cosPhi * cosTheta,
            a[0] * cosPhi * sinTheta,
            a[0] * sinPhi
            ])

    if a.ndim is 2:
        return interleave(res)
    else:
        return res


def cart_to_spher(a):
    """cart_to_spher(a)
    
    Return spherical coordinates from Cartesian.
    Angles in radians.

    [ x, y, z ] ---> [ r, theta, phi ]

    """
    if a.ndim > 2:
        raise "rank > 2 not supported"
    else:
        if a.ndim is 2:
            a = deinterleave(a)

    r = sqrt(sum(asarray(a)**2, axis = 0))

    if r is 0.:
        theta = phi = 0.
    elif r is zeros_like(r):
        theta = phi = zeros_like(r)
    else:
#         theta = arctan2( a[1], a[0] )
#         phi = arcsin( a[2] / r )
        theta = nan_to_num(arctan2( a[1], a[0] ))
        phi = nan_to_num(arcsin( a[2] / r ))

#     res = array([ r, theta, phi ])

    if a.ndim is 2:
        return interleave(array([ r, unwrap(theta), unwrap(phi) ]))
    else:
        return array([ r, theta, phi ])


#=========================
# Unary Functions
# 
# To operate on soundfile arrays
# 
# Note: to broadcast correctly, these functions expect
# numpy arrays of shape (nframes, nchannels) OR (nframes)
#=========================

def nchannels(a):
    """nchannels(a)
    
    Returns number of channels of input.
    """
    if a.ndim > 2:
        raise "rank > 2 not supported"

    else:
        if a.ndim is 2:
            nchannels   = a.shape[1]
        else:
            nchannels   = 1

    return nchannels


def nframes(a):
    """nframes(a)
    
    Returns length of input (in nframes).
    """
    return a.shape[0]


def interleave(a, flat = True):
    """interleave(a)
    
    Interleave array of shape (nchannels, nframes)
    OR shape (nframes).

    Returning an array of shape (nframes, nchannels)

    If flat is True, mono will remain as flat array

    """
    if a.ndim > 2:
        raise "rank > 2 not supported"

    elif (nchannels(a) is 1) and flat:
        return reshape(a, (-1, 1))

    else:
        return transpose(a)


def deinterleave(a, flat = True):
    """deinterleave(a)
    
    De-interleave array of shape (nframes, nchannels)

    Returning an array of shape (nchannels, nframes)
    OR shape (nframes).

    If flat is True, deinterleaves mono to flat array

    """
    if a.ndim > 2:
        raise "rank > 2 not supported"

    elif (nchannels(a) is 1) and flat:
            return reshape(a, (-1))

    else:
        return transpose(a)


def peak(a, axis = None):
    """peak(a, axis = None)
    
    Return positive or negative peak value along dimension axis.

    """
    pp = a.max(axis)
    pn = a.min(axis)

    if axis is None or a.ndim is 1:
        if pp >= abs(pn):
            p = pp
        else:
            p = pn
    else:
        p = where(pp >= abs(pn), pp, pn)

    if (nchannels(a) is 1) and (axis is None or axis is 0):
        p = float(p)

    return p


def argpeak(a, axis = None):
    """argpeak(a, axis = None)
    
    Returns the indices of the peak value of the
    1-D arrays along the given axis.

    If axis = None, the index is into the flattened array.

    """
    p = peak(a, axis)

    app = a.argmax(axis)
    apn = a.argmin(axis)

    if axis is None or a.ndim is 1:
        if p >= 0:
            i = app
        else:
            i = apn
    else:
        i = where(p >= 0, app, apn)

    return i


def rms(a, axis = None):
    """rms(a, axis = None)
    
    Return rms value along dimension axis.

    If axis = None, returns the peak rms along
    the 0 axis.

    """
    if axis is None:
        axis = 0
        res = sqrt((a**2).mean(axis)).max()

    else:
        res = sqrt((a**2).mean(axis))

    if (nchannels(a) is 1) and (axis is None or axis is 0):
        res = float(res)

    return res


def scale_to(a, scale = 1., method = 'peak', axis = None):
    """scale_to(a, scale = 1., method = 'peak', axis = None)
    
    Normalize array to scale along axis using 'peak' or 'rms' method.

    If axis = None, normalize all channels to a single value.
    If axis = axis, normalize each channel separately along axis.
    """

    if method is 'rms':
        scale *=reciprocal(rms(a, axis))
    else:
        scale *=reciprocal(abs(peak(a, axis)))

    if axis is None or 0:
        res = scale * a
    else:
        res = rollaxis(
            scale * rollaxis(a, axis),
            0, axis + 1)

    return res


#=========================
# Binary Functions
# 
# To operate on TWO soundfile arrays
# 
# Note: to broadcast correctly, these functions expect
# numpy arrays of shape (nframes, nchannels) OR (nframes)
#=========================


def over_dub(a, b, index = 0, write_over = False):
    """over_dub(a, b, index = 0, write_over = False)
    
    Add b to a starting at the index.
    If b is too long only the first part is overdubbed.
    If write_over = True, write over input a.

    """
    if index >= nframes(a):         # error check
        raise ValueError, ("index must be < nframes(a)")
    elif nframes(b) > nframes(a):
        raise ValueError, ("nframes(b) must be < nframes(a)")
    else:                       # set slices
        if nframes(b) + index > nframes(a):
            s_b = slice(0, nframes(a) - index)
            s_a = slice(index, None)
        else:
            s_b = Ellipsis
            s_a = slice(index, index + nframes(b))

        if write_over is True:
            res = a
        else:
            res = a.copy()

        res[s_a] = a[s_a] + b[s_b]

    return res


def over_write(a, b, index = 0, write_over = False):
    """over_write(a, b, index = 0, write_over = False)
    
    Over write b into a starting at the index.
    If b is too long only the first part is written.
    If write_over = True, write over input a.

    """
    if index >= nframes(a):         # error check
        raise ValueError, ("index must be < nframes(a)")
    elif nframes(b) > nframes(a):
        raise ValueError, ("nframes(b) must be < nframes(a)")
    else:                       # set slices
        if nframes(b) + index > nframes(a):
            s_b = slice(0, nframes(a) - index)
            s_a = slice(index, None)
        else:
            s_b = Ellipsis
            s_a = slice(index, index + nframes(b))

        if write_over is True:
            res = a
        else:
            res = a.copy()

        res[s_a] = b[s_b]

    return res


#=========================
# Matrix Utility Functions
# 
#=========================

def mmix(a, b):
    """
    Args:
        - a         : Input multi-channel signal (interleaved)
        - b         : a matrix, may be a single matrix, or an array of
                      nframes(b) = nframes(a)

    Matrix mix a multi-channel signal. (Encode, Transform, Decode)
    
    Returs dot product, a[0] . b

    """
    
    # test shape
    dim_a = len(shape(a))
    dim_b = len(shape(b))
    
    if (dim_a is 2) and (dim_b is 2):
        res = einsum(a, [0,1], b, [2,1])

    elif (dim_a is 2) and (dim_b is 3):
        res = einsum(a, [0,1], b, [0,3,1], [0, 3])
        
    return res
    

def mmul(a, b):
    """
    Args:
        - a         : a matrix, may be a single matrix, or an array of
                      nframes(a)
        - b         : a matrix, may be a single matrix, or an array of
                      nframes(b)

    Returns dot product, b . a
    
    NOTE: if nframes(a) or nframes(b) !=, broadcasting will fail
    """
    
    # test shape
    dim_a = len(shape(a))
    dim_b = len(shape(b))
    
    if (dim_a is 2) and (dim_b is 2):
        res = einsum(b, [0, 1], a, [1, 3])

    elif (dim_a is 3) and (dim_b is 2):
        res = einsum(b, [0, 1], a, [2, 1, 4], [2, 0, 4])

    elif (dim_a is 2) and (dim_b is 3):
        res = einsum(b, [0, 1, 2], a, [2, 4], [0, 1, 4])

    elif (dim_a is 3) and (dim_b is 3):
        res = einsum(b, [0, 1, 2], a, [0, 2, 4], [0, 1, 4])
  
    return res


# consider creating a muse signal subclass of ndarray
# advantage would be to include the above as methods of
# the new signal class--could be very useful!
# this may be slightly involved-->review _wrapit
# see: "/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/numpy/core/fromnumeric.py"
# this is about wrapping methods to become functions,
# but could do the other way around

# class msig(ndarray):

#     pass
#     def __init__(self, vals, nframes, sr, gen = True):
#         MusObj.__init__(self, sr) # Run superclass init

#         self.vals = asarray(vals) # store vals

#         if isscalar(nframes):
#             self.nframes = nframes
#         else:
#             self.nframes = asarray(nframes)

#         self.nchannels = nchannels(asarray(vals)) # store nchannels

#         if gen is True:
#             self._gen_set()            # initialize / reset generator
