#! /usr/bin/env python

# vim:syntax=python

# Joseph Anderson 2013


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on pyaudiolab.
# """


# seem to need to have these imports here. . . to make sure names are defined
from muse import *
#from math import factorial
from scipy.misc import factorial

from muse.hoa import *

# import muse defined constants
import muse.constants as C



#=========================
# Definition of constants
#=========================

# scaling to convert from SN3D to maxN
# indexed by FuMa channel number
# NOTE: should not be directly touched by user
_fuma_maxN_weights = array([
    1,              # W : order 0
    1,              # X : order 1
    1,              # Y
    1,              # Z
    1,              # R : order 2
    2. / sqrt(3),   # S
    2. / sqrt(3),   # T
    2. / sqrt(3),   # U
    2. / sqrt(3),   # V
    1,              # K : order 3
    sqrt(45./32.),  # L
    sqrt(45./32.),  # M
    3. / sqrt(5),   # N
    3. / sqrt(5),   # O
    sqrt(8./5.),    # P
    sqrt(8./5.),    # Q
])
_fuma_MaxN_weights = concatenate((
    array([1./sqrt(2)]),
    _fuma_maxN_weights[1:]
))

#=========================
# Functions
#=========================


#=========================
# Spherical Harmonics - Coefficients

# Generate un-normalised spherical harmonic coefficients
#
# SEE: http://en.wikipedia.org/wiki/Ambisonic_data_exchange_formats
#
# AND:
#
# Nachbar, C. et al., 2011. AmbiX - A Suggested Ambisonics Format.
# In Proceedings of the 3rd International Symposium on Ambisonics and
# Spherical Acoustics. 3rd International Symposium on Ambisonics and
# Spherical Acoustics. Lexington, Kentucky.
#
#                                               / sin(abs(m) * theta)  -> m < 0
# sph_harm = norm(l, m) * alp(l, m, sin(phi)) * | 1                    -> m = 0
#                                               \ cos(abs(m) * theta)  -> m > 0
#
# associated legendre function (alp)
# scipy.special.lpmv(m,v,x)
#     returns the associated legendre function of
#     integer order m and real degree v
#
# NOTE: includes Condon-Shortley phase factor
#
# l = degree (classic ambisonic order)
# m = order (classic ambisonic index)
#
# azimuth --> 0 - 2pi
# elevation --> -pi/2 - pi/2
#
# returns unnormalised coefficients
def spher_harm(lm = (0, 0), theta = 0, phi = 0):
    """Args:
        - l         : Associated Legendre degree (ambisonic 'order')
        - m         : Associated Legendre order (ambisonic 'index')
        - theta     : azimuth, counter-clockwise, in radians
        - phi       : elevation, upwards, in radians

    The harmonic $Y^l_m$ sampled at theta and phi.

    Returns the un-normalised (real) spherical harmonic coefficient
    in the form required for Higher Order Ambisonics.
    """
    l, m = lm
    
    # for testing for, and constructing vectors if need be
    theta_scalar = isscalar(theta)
    phi_scalar = isscalar(phi)
    
    if not theta_scalar and phi_scalar:
        phi = repeat(array([phi]), len(theta))
    elif theta_scalar and not phi_scalar:
        theta = repeat(array([theta]), len(phi))
                
#    # remap.... via math library
#    theta = arctan2(
#        sin(theta) * cos(phi),
#        cos(theta) * cos(phi)
#    )
#    phi = arcsin(sin(phi))

    # remap... via 'hand'
    theta = theta - (pi * clip(sign(-cos(phi)), 0, 1))
    theta = mod(theta + pi, twoPi) - pi
    phi = (pi/2) - absolute(
        mod(phi + (pi/2), twoPi) - pi
    )

    # if theta / phi are arrays interleave here...
    if (not theta_scalar or not phi_scalar) and not isscalar(l):
        theta = interleave(theta)
        phi = interleave(phi)

    a = pow(-1, m) * scipy.special.lpmv(
        absolute(m),
        l,
        sin(phi)
    )

    boolean = less(m, 0)
    
    b = (sin(absolute(m) * theta) * boolean) + \
        (cos(absolute(m) * theta) * logical_not(boolean))
    
    res = a * b
    
    return res


#=========================
# Spherical Harmonics - Normalisation

# Schmidt semi-normalisation (3D)
def sn3d(lm = (0, 0)):
    """Args:
        - l              : Associated Legendre degree (ambisonic 'order')
        - m              : Associated Legendre order (ambisonic 'index')

    SN3D: returns Schmidt semi-normalisation (3D) weighting coefficients
    for spherical harmonics in the form required for Higher Order Ambisonics.
    
    This scheme has been adopted by the proposed AmbiX format.
    """
    l, m = lm
    
    dm = equal(m, 0).astype(int)
    
    res = sqrt(
        (2 - dm) * \
        (factorial(l - absolute(m)) / factorial(l + absolute(m)))
    )
        
    return res


# full 3D normalisation
def n3d(lm = (0, 0)):
    """Args:
        - l              : Associated Legendre degree (ambisonic 'order')
        - m              : Associated Legendre order (ambisonic 'index')

    N3D: returns full 3D normalisation weighting coefficients
    for spherical harmonics in the form required for Higher Order Ambisonics.
    
    Daniel describes it as follows: "Orthonormal basis for 3D decomposition.
    Simple relationship to SN3D [..]. Ensures equal power of the encoded
    components in the case of a perfectly diffuse 3D field. [..]
    Obvious significance for solving decoding problems [..]
    (3D reconstruction)."
    """
    l, m = lm
    
    res = sqrt(2 * l + 1) * sn3d(lm)
    
    return res


# full 2D normalisation
def n2d(lm = (0, 0)):
    """Args:
        - l              : Associated Legendre degree (ambisonic 'order')
        - m              : Associated Legendre order (ambisonic 'index')

    N2D: returns full 2D normalisation weighting coefficients
    for spherical harmonics in the form required for Higher Order Ambisonics.
    
    Daniel states: "2D-restricted formalism (cylindrical harmonics), but
    may also apply to 3D spherical harmonics provided that clear extension
    rules are given..."

    Useful for solving decoding problems [..]
    (2D reconstruction).
    """
    l, m = lm
    
    res = sqrt(
        pow(2, 2 * l) * pow(factorial(l), 2) / factorial(2 * l + 1)
    ) * n3d(lm)
    
    return res


# Semi-normalisation (2D)
def sn2d(lm = (0, 0)):
    """Args:
        - l              : Associated Legendre degree (ambisonic 'order')
        - m              : Associated Legendre order (ambisonic 'index')

    SN2D: returns Semi-normalisation (2D) weighting coefficients
    for spherical harmonics in the form required for Higher Order Ambisonics.
    
    Daniel states: "2D-restricted formalism (cylindrical harmonics), but
    may also apply to 3D spherical harmonics provided that clear extension
    rules are given..."

    Useful for solving decoding problems [..]
    (2D reconstruction).
    """
    l, m = lm
    
    lne0 = not_equal(0, l).astype(int)
    
    res = pow(2, -1./2 * lne0) * n2d(lm)
        
    return res


# Furse (29 June 2014) says:
# For MaxN, you need the maxima of the underlying spherical harmonics you're
# using. To find these, you can use the fact that the partial derivatives
# need to be zero. The azimuth derivative is trivial; the tricky bit comes
# from the ALPs, though this can be done by differentiating wrt sin(elev)
# and solving the resulting polynomial equations. Not particularly fun ;-)

# Heller (2 July 2014):
# Delivered file "MaxN Normalization.nb.pdf" listing coefficients up to
# 8th-order. Solved via Mathematica.

# maxN & MaxN normalisation
def maxN(lm = (0, 0)):
    """Args:
        - l              : Associated Legendre degree (ambisonic 'order')
        - m              : Associated Legendre order (ambisonic 'index')

    maxN: returns 'maxN' normalisation weighting coefficients
    for spherical harmonics in the form required for Higher Order Ambisonics.
    0th harmonic (W) scaled to 1.
    
    Normalises the maximum value for for each harmonic to 1. This is the
    normalisation scaling for the Furse-Malham convention.
    """
    res = sn3d(lm)  * _fuma_maxN_weights[lm_to_fuma(lm)]
    
    return res


def MaxN(lm = (0, 0)):
    """Args:
        - l              : Associated Legendre degree (ambisonic 'order')
        - m              : Associated Legendre order (ambisonic 'index')

    maxN: returns 'maxN' or 'MaxN' normalisation weighting coefficients
    for spherical harmonics in the form required for Higher Order Ambisonics.
    0th harmonic (W) scaled to 1./sqrt(2)
    
    Normalises the maximum value for for each harmonic to 1. This is the
    normalisation scaling for the Furse-Malham convention.
    """
    res = sn3d(lm)  * _fuma_MaxN_weights[lm_to_fuma(lm)]
    
    return res


#=========================
# Spherical Harmonics - Channel Order

# ACN (Ambisonic Channel Number) Utilities
def lm_to_acn(lm = (0, 0)):
    """Args:
        - l         : Associated Legendre degree (ambisonic 'order')
        - m         : Associated Legendre order (ambisonic 'index')

    Returns Ambisonic Channel Number (ACN) for Higher Order Ambisonics.
    """
    l, m = lm
    
    res = l**2 + l + m
    
    return res
    

def acn_to_lm(acn = 0):
    """Args:
        - acn       : Ambisonic Channel Number

    Returns spherical harmonic degree l (ambisonic 'order') and
    order m (ambisonic 'index') from Ambisonic Channel Number (ACN).
    """

    l = (floor(sqrt(acn))).astype(int)
    m = acn - l**2 - l
    
    res = (l, m)

    return res


# SID (Single Index Designation) Utilities
def lm_to_sid(lm = (0, 0)):
    """Args:
        - l         : Associated Legendre degree (ambisonic 'order')
        - m         : Associated Legendre order (ambisonic 'index')

    Returns Ambisonic Single Index Designation (SID) for Higher Order Ambisonics.
    """
    l, m = lm
    
    res = l**2 + (2 * (l - absolute(m)))

    res = res - clip(sign(m), -1, 0)


    return res


def sid_to_lm(sid = 0):
    """Args:
        - sid       : Ambisonic Single Index Designation

    Returns spherical harmonic degree l (ambisonic 'order') and
    order m (ambisonic 'index') from Ambisonic
    Single Index Designation (SID).
    """

    l = (floor(sqrt(sid))).astype(int)

    m0 = (floor(((l + 1)**2 - (sid + 1)) / 2.)).astype(int)
    m1 = -1 * (floor(((l + 1)**2 - (sid)) / 2.)).astype(int)
    
    boolean = equal(m0, absolute(m1))
    
    m = (m0 * boolean) + (m1 * logical_not(boolean))
    
    res = (l, m)

    return res


# FuMa (Furse-Malham) Utilities
def lm_to_fuma(lm = (0, 0)):
    """Args:
        - l         : Associated Legendre degree (ambisonic 'order')
        - m         : Associated Legendre order (ambisonic 'index')

    Returns Ambisonic Furse-Malham ordering for Higher Order Ambisonics.
    """
    l, m = lm
    
    # ordering for <=1st order
    res1 = lm_to_sid(lm)

    # ordering for >1st order
    res2 = l**2 + (2 * absolute(m))

    res2 = res2 - clip(sign(m), 0, 1)

    boolean = less(l, 2)
    
    res = (res1 * boolean) + (res2 * logical_not(boolean))

    return res


def fuma_to_lm(fuma = 0):
    """Args:
        - sid       : Ambisonic Furse-Malham channel number

    Returns spherical harmonic degree l (ambisonic 'order') and
    order m (ambisonic 'index') from Furse-Malham channel number (FuMa).
    """

    l = (floor(sqrt(fuma))).astype(int)

    m0 = -1 * (floor((fuma - l**2) / 2.)).astype(int)
    m1 = (floor((fuma + 1 - l**2) / 2.)).astype(int)

    boolean = equal(m1, absolute(m0))
    
    m = (m0 * boolean) + (m1 * logical_not(boolean))

    # correct m for <=1st order
    boolean2 = less(l, 2)

    m = (sid_to_lm(fuma)[1] * boolean2) + (m * logical_not(boolean2))
    
    res = (l, m)

    return res


#=========================
# Encoding - Matricies

# NOTE: suitable for both encoding, transforming, decoding
#       may want to move to 'utilities' 
def encoding_convert_matrix(format_in, format_out, order = 1):
    """Args:
        - format_in     : (ordering, normalisation)
        - format_out    : (ordering, normalisation)
        - order         : HOA order

        ordering      : 'acn', 'sid', 'fuma'
        normalisation : 'sn3d', 'n3d', 'sn2d', 'n2d', 'maxN', 'MaxN'


    Generate a matrix to encode / transcode an HOA signal.
    NOTE: AmbiX format = ('acn', 'sn3d')
    
    Use in conjunction with mmix().
    """
    N = order

    ordering_in, norm_in = format_in
    ordering_out, norm_out = format_out

    # ordering: acn, sid, fuma
    ordering_to_lm_dict = {
        'acn': acn_to_lm,
        'sid': sid_to_lm,
        'fuma': fuma_to_lm,
    }

    lm_to_ordering_dict = {
        'acn': lm_to_acn,
        'sid': lm_to_sid,
        'fuma': lm_to_fuma,
    }

    # normalisation: sn3d, n3d, maxN, MaxN
    lm_to_norma_dict = {
        'sn3d': sn3d,
        'n3d': n3d,
        'sn2d': sn2d,
        'n2d': n2d,
        'maxN': maxN,
        'MaxN': MaxN,
    }
    
    nchans = (N+1)**2
    lm = ordering_to_lm_dict[ordering_in](arange(nchans))
    norm = lm_to_norma_dict[norm_out](lm) / lm_to_norma_dict[norm_in](lm)

    res = identity(nchans)[:, lm_to_ordering_dict[ordering_out](lm)] * norm

    return res
    

def planewave_matrix(direction = array([0, 0]), order = 1):
    """
    Args:
        - direction : [azimuth, elevation]
        - order     : HOA order
    
    Direction argument may take several forms:
        
        shape(direction) == (2,): 
            encode one input as a TI, single planewave

        shape(direction) == (nframes, 2):
            encode one input as a TV, single planewave

        shape(direction) == (1, npws, 2):
            encode npws input as TI, multiple planewaves

        shape(direction) == (nframes, npws, 2):
            encode npws input as TV, multiple planewaves

    
    Generate a matrix to encode a signal (mono or multi-channel)
    into the B-format domain (AmbiX).
    
    Use in conjunction with mmix() to encode a signal.
    """
    # case 1: one input TI, single planewave
    #           - shape(dir) = (2,)
    #
    # case 2: one input TV, single planewave
    #           - shape(dir) = (nframes, 2)
    #
    # case 3: npws input TI, multiple planewaves
    #           - shape(dir) = (1, npws, 2)
    #
    # case 4: npws input TV, multiple planewaves
    #           - shape(dir) = (nframes, npws, 2)

    N = order
    
    theta, phi = transpose(direction)

    lm = acn_to_lm(arange((1+N)**2))
    norm = sn3d(lm)
    
    # case 4
    if len(shape(direction)) == 3 and shape(direction)[0] != 1:
        
        nframes, npws, dim = shape(direction)
        res = zeros((npws, nframes, (N+1)**2))
        
        # itterate by planewave
        for i in range(npws):
            res[i] = norm * spher_harm(lm, theta[i], phi[i])
    
        # reshape to return (nframes, out, in)
        res = transpose(res, (1, 2, 0))

    # cases 1, 2, 3
    else:
        res = norm * spher_harm(lm, theta, phi)
    
        if len(shape(direction)) == 3 and shape(direction)[0] == 1:
            # case 3
            res = transpose(res)
    
        else:        
            # case 1, 2
            res = reshape(
                res,
                shape(res) + (1,)
            )

    return res


# NOTE: may want to include pressurewave_matrix (omni) encoder
#       AND have a parameter to normalise according to:
#           'amp', 'rms', 'energy' - see decoders


#---------------------------------------------
# A to B encoder
#---------------------------------------------

# NOTE: For A to B encoding will need access to spherical sampling.
#       Incorporate platonic solids, t-design and other (polyhedral)
#       sampling methods.