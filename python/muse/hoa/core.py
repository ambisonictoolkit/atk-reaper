#! /usr/bin/env python
# -*- coding: utf-8 -*-

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

from numpy.polynomial import Chebyshev as T
from numpy.polynomial import Legendre as L
from numpy import math
#from math import factorial
from scipy.misc import factorial

from numpy.linalg.linalg import pinv


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
# Definition of dictionaries
#=========================

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
lm_to_normal_dict = {
    'sn3d': sn3d,
    'n3d': n3d,
    'sn2d': sn2d,
    'n2d': n2d,
    'maxN': maxN,
    'MaxN': MaxN,
}


#=========================
# Encoding - Matricies

def format_matrix(order, format_in, format_out):
    """Args:
        - order         : HOA order
        - format_in     : (ordering, normalisation)
        - format_out    : (ordering, normalisation)

        ordering      : 'acn', 'sid', 'fuma'
        normalisation : 'sn3d', 'n3d', 'sn2d', 'n2d', 'maxN', 'MaxN'


    Generate a matrix to encode / transcode an HOA signal from one
    HOA format to another.

    NOTE: AmbiX format = ('acn', 'sn3d')
    
    Use in conjunction with mmix().
    """
    N = order

    ordering_in, norm_in = format_in
    ordering_out, norm_out = format_out
    
    nchans = (N+1)**2
    lm = ordering_to_lm_dict[ordering_in](arange(nchans))
    norm = lm_to_normal_dict[norm_out](lm) / lm_to_normal_dict[norm_in](lm)

    res = identity(nchans)[:, lm_to_ordering_dict[ordering_out](lm)] * norm

    return res


def planewave_matrix(direction = array([0, 0]), order = 1, \
    format = ('acn', 'sn3d')):
    """
    Args:
        - direction : [azimuth, elevation]
        - order     : HOA order
        - format    : (ordering, normalisation)

    ordering      : 'acn', 'sid', 'fuma'
    normalisation : 'sn3d', 'n3d', 'sn2d', 'n2d', 'maxN', 'MaxN'
        
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
    into the B-format domain.
    
    format = ('acn', 'sn3d') returns AmbiX, ATK's convention.
    
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
    ordering, normalisation = format
    
    theta, phi = transpose(direction)

    nchans = (N+1)**2

    lm = ordering_to_lm_dict[ordering](arange(nchans))
    norm = lm_to_normal_dict[normalisation](lm)
    
    # case 4
    if len(shape(direction)) == 3 and shape(direction)[0] != 1:
        
        nframes, npws, dim = shape(direction)
        res = zeros((npws, nframes, nchans))
        
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



#=========================
# Decoding - Matricies

# --------------------------------------------------------------------------
#
# A set of utilities derived from BLaH and Daniel to calculate values
# required for HOA decoding. These are broken down into the following
# functions:
#
#     decoder_rV                : calculate rV for a 'regular' decoder
#     decoder_rE                : calculate rE for a 'regular' decoder
#     decoder_order_gains       : calculate 'order gains', g_m
#     decoder_E                 : decoder 'E', 'reduced energy'
#     decoder_matching_gain     : calculate gain to 'match' decoder types

## Heller, Aaron J., Eric M. Benjamin, and Richard Lee. 2012. “A Toolkit for
# the Design of Ambisonic Decoders.” In Proceedings of the Linux Audio
# Conference 2012. Stanford, CA.
#
# Daniel, Jérôme. 2001. “Représentation de champs acoustiques, application
# à la transmission et à la reproduction de scènes sonores complexes dans
# un contexte multimédia”. PhD Thesis, Paris: Université Paris 6.
#
# --------------------------------------------------------------------------


#=========================
# Functions - decoding utilities
#=========================

# decoder_rV            : calculate rV for a 'regular' decoder
def decoder_rV(order = 1, dec_type = 'basic', dim = 2):
    """
    Args:
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - dim        : Decoder dimensions: 2D or 3D


    Compute maximum average rV for an Ambisonic decoder.
    """
    
    M = order
    
    if dec_type == 'basic':
        res = 1

    elif dec_type == 'energy':
        res = decoder_rE(order, dec_type, dim)

    elif dec_type == 'controlled':
        if dim == 2:
            res = M / (M + 1.)
        else:
            res = M / (M + 2.)
    
    return res


# decoder_rE            : calculate rE for a 'regular' decoder
def decoder_rE(order = 1, dec_type = 'basic', dim = 2):
    """
    Args:
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - dim        : Decoder dimensions: 2D or 3D


    Compute maximum average rE for an Ambisonic decoder.
    """
    
    M = order
    
    if dec_type == 'basic' or dec_type == 'controlled':
        if dim == 2:
            res = 2*M / (2*M + 1.)
        else:
            res = M / (M + 1.)

    elif dec_type == 'energy':
        if dim == 2:
            res = max(T.basis(M + 1).roots()) # BLaH
            #res = cos(pi / (2 * M + 2)) # Daniel
        else:
            res = max(L.basis(M + 1).roots())
    
    return res


# decoder_order_gains       : calculate 'order gains', g_m
def decoder_order_gains(order = 1, dec_type = 'basic', dim = 2):
    """    
    Args:
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - dim        : Decoder dimensions: 2D or 3D


    Compute 'order gains' for a given Ambisonic decoder.
    
    Returns array
        G = [g_0, g_1, g_2, ... g_M]
        
    Where
        M = order
        g_m = scaling gain (scale) for each order coefficients
    """
    
    M = order

    res = ones(M + 1)
    
    if dec_type == 'basic':
        #res = ones(M + 1)
        pass

    elif dec_type == 'energy':
        if dim == 2:
            for m in range(M + 1):
                res[m] = T.basis(m)(decoder_rE(M, 'energy', 2)) # BLaH
                #res[m] = cos((m*pi) / (2*M +2)) # Daniel
        else:
            for m in range(M + 1):
                res[m] = L.basis(m)(decoder_rE(M, 'energy', 3))

    elif dec_type == 'controlled':
        if dim == 2:
            for m in range(M + 1):
                res[m] = 1./(factorial(M+m) * factorial(M-m))
            res *= factorial(M)**2
        else:
            for m in range(M + 1):
                res[m] = 1./(factorial(M+m+1) * factorial(M-m))
            res *= factorial(M)*factorial(M+1)
    
    return res


# decoder_E                 : decoder 'E'
def decoder_E(order = 1, dec_type = 'basic', dim = 2):
    """
    Args:
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - dim        : Decoder dimensions: 2D or 3D


    Compute 'l’énergie réduite E' for a given Ambisonic decoder.
    """
    
    M = order

    g_m = decoder_order_gains(M, dec_type, dim)
    
    if dim == 2:
        res = g_m[0]**2
        for m in range(1, M + 1):
            res += 2*(g_m[m]**2)

    elif dim == 3:
        res = 0.
        for m in range(M + 1):
            res += (2*m +1)*(g_m[m]**2)

    return res


# decoder_matching_gain     : calculate gain to 'match' decoder types
def decoder_matching_gain(order = 1, dec_type = 'basic', dim = 2, \
    match = 'amp', num_spkrs = nan):
    """
    Args:
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - dim        : Decoder dimensions: 2D or 3D
        
        - match      : Decoder matching criteria
                        'amp'    --> amplitude
                        'rms'    --> RMS (Gerzon / classic)
                        'energy' --> energy

        - num_spkrs : Number of decoder speakers


    Compute 'matching gain' (scale) for a given Ambisonic decoder.
    """
        
    M = order

    if match == 'amp':
        res = 1.

    elif match == 'rms':
        if dim == 2: 
            N = 2*M + 1

        elif dim == 3:
            N = (M + 1)**2

        res = sqrt(N/decoder_E(M, dec_type, dim))

    elif match == 'energy':
        N = num_spkrs
        res = sqrt(N/decoder_E(M, dec_type, dim))

    return res


#=========================
# Decoding - Matricies (utilities)

def order_gains_matrix(order = 1, dec_type = 'basic', dim = 2):
    """    
    Args:
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - dim        : Decoder dimensions: 2D or 3D


    Compute 'order gains' for a given Ambisonic decoder.
    
    Returns order gain matrix. See decoder_order_gains()
    """
    M = order

    res = diag(
        decoder_order_gains(
            M,
            dec_type,
            dim
        )[acn_to_lm(arange((M+1)**2))[0]])

    return res
    

def resize_order_matrix(order_in, order_out):
    """    
    Args:
        - order_in   : Ambisonic order, e.g., 1, 2, 3...

        - order_out  : Ambisonic order, e.g., 1, 2, 3...


    Resize Ambisonic order:
        order_in > order_out --> discard harmonics
        order_in < order_out --> insert zeros
    
    Returns matrix.
    """
    M_in = order_in
    M_out = order_out
    
    nchans_in = (M_in+1)**2
    nchans_out = (M_out+1)**2
    
    res = eye(nchans_out, nchans_in)
    
    return res


def zero_order_matrix(order_in, order_out):
    """
    Args:
        - order_in   : Ambisonic order, e.g., 1, 2, 3...

        - order_out  : Ambisonic order, e.g., 1, 2, 3...


    Zero (truncate) Ambisonic order:
        order_in > order_out --> insert zeros
        order_in < order_out --> invalid, returns identity matrix
    
    Returns (square) matrix.
    """
    M_in = order_in
    M_out = order_out
    
    nchans_in = (M_in+1)**2
    nchans_out = (M_out+1)**2
    
    if nchans_out > nchans_in:
        nchans_out = nchans_in
    
    res = diag(
        concatenate((
            ones(nchans_out),
            zeros(nchans_in - nchans_out)
        ))
    )
    
    return res


def peri_to_panto_matrix(order, ordering = 'acn'):
    """
    Args:
        - order     : Ambisonic order, e.g., 1, 2, 3...
        - ordering  : 'acn', 'sid', 'fuma'


    Discard periphonic (3D) harmonics. No (cylindrical) weighting
    is applied.
    
    Returns matrix.
    """
    M = order
    
    nchans_in = (M+1)**2
    
    l, m = ordering_to_lm_dict[ordering](arange(nchans_in))
    
    res = identity(nchans_in)[((l == abs(m)))]

    return res


def panto_to_peri_matrix(order, ordering = 'acn'):
    """
    Args:
        - order   : Ambisonic order, e.g., 1, 2, 3...
        - ordering  : 'acn', 'sid', 'fuma'


    Insert periphonic (3D) harmonics, as zeros. No (spherical) weighting
    is applied.
    
    Returns matrix.
    """
    M = order
    
    res = transpose(peri_to_panto_matrix(M, ordering))

    return res


def zero_peri_matrix(order, ordering = 'acn'):
    """
    Args:
        - order   : Ambisonic order, e.g., 1, 2, 3...
        - ordering  : 'acn', 'sid', 'fuma'


    Zero periphonic (3D) harmonics. No (cylindrical) weighting
    is applied.
    
    Returns (square) matrix.
    """
    M = order
    
    nchans_in = (M+1)**2
    
    l, m = ordering_to_lm_dict[ordering](arange(nchans_in))
    
    res = diag((l == abs(m)).astype(int))

    return res


#=========================
# Decoding - Matricies (decoders)

def panto_pinv_decoding_matrix(directions, order = 1, dec_type = 'basic', \
    match = 'amp', format = ('acn', 'sn3d')):
    """
    Args:
        - direction : [[azimuth, elevation], ... ]
    
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - match      : Decoder matching criteria
                        'amp'    --> amplitude
                        'rms'    --> RMS (Gerzon / classic)
                        'energy' --> energy

        - format    : (ordering, normalisation)

    ordering      : 'acn', 'sid', 'fuma'
    normalisation : 'sn3d', 'n3d', 'sn2d', 'n2d', 'maxN', 'MaxN'

    Returns pseudo-inverse matrix Ambisonic decoder.
    """

    dim = 2 # 2D
    num_spkrs = len(directions)
    
    N_in = order

    # calculate decoder order (N_out)
    if num_spkrs >= (2*N_in + 1):
        N_out = N_in
    else:
        N_out = (num_spkrs - 1) / 2

    # --------------------------------
    # prototype encoding matrix
    encoding_matrix = planewave_matrix(
        array([directions]),
        N_out,
        format
    )
        
    # discard 3D harmonics
    encoding_matrix = mmul(
        encoding_matrix,
        peri_to_panto_matrix(N_out, format[0])
    )
    
    # decoder: pseudo inverse - via pinv()
    decoding_matrix = pinv(encoding_matrix)
    
    # (re-)insert 3D harmonics
    decoding_matrix = mmul(
        peri_to_panto_matrix(N_out, format[0]),
        decoding_matrix
    )
    
    # apply order gains - for decoder type
    decoding_matrix = mmul(
        order_gains_matrix(N_out, dec_type, dim),
        decoding_matrix
    )
    
    # apply matching gain
    decoding_matrix = decoder_matching_gain(
        N_out,
        dec_type,
        dim,
        match,
        num_spkrs
        ) * decoding_matrix
        
    # expand to match input order (if necessary)
    if N_out != N_in:
        decoding_matrix = mmul(
            resize_order_matrix(N_in, N_out),
            decoding_matrix
        )
    
    return decoding_matrix


def panto_sad_decoding_matrix(directions, order = 1, dec_type = 'basic', \
    match = 'amp', format = ('acn', 'sn3d')):
    """
    Args:
        - direction : [[azimuth, elevation], ... ]
    
        - order      : Ambisonic order, e.g., 1, 2, 3...

        - dec_type   : Decoder type
                        'basic'      --> basic, velocity, rV
                        'energy'     --> energy, rE
                        'controlled' --> controlled opposites, in phase
                            
        - match      : Decoder matching criteria
                        'amp'    --> amplitude
                        'rms'    --> RMS (Gerzon / classic)
                        'energy' --> energy

        - format    : (ordering, normalisation)

    ordering      : 'acn', 'sid', 'fuma'
    normalisation : 'sn3d', 'n3d', 'sn2d', 'n2d', 'maxN', 'MaxN'

    Returns projection, aka 'simple Ambisonic decoder' (SAD),
    matrix Ambisonic decoder.
    """

    dim = 2 # 2D
    num_spkrs = len(directions)

    N_in = order

    # calculate decoder order (N_out)
    if num_spkrs >= (2*N_in + 1):
        N_out = N_in
    else:
        N_out = (num_spkrs - 1) / 2

    # --------------------------------
    # prototype encoding matrix ((ordering), ('n2d'))
    encoding_matrix = planewave_matrix(
        array([directions]),
        N_out,
        (format[0], 'n2d')
    )
    
    # zero 3D harmonics
    encoding_matrix = mmul(
        encoding_matrix,
        zero_peri_matrix(N_out, format[0])
    )
    
    # transpose and scale
    decoding_matrix = (transpose(encoding_matrix)) / num_spkrs
    
    # convert to n2d scaling - for input
    decoding_matrix = mmul(
        format_matrix(
            N_out,
            format,
            (format[0], 'n2d')
        ),
        decoding_matrix
    )    

    # apply order gains - for decoder type
    decoding_matrix = mmul(
        order_gains_matrix(N_out, dec_type, dim),
        decoding_matrix
    )
    
    # apply matching gain
    decoding_matrix = decoder_matching_gain(
        N_out,
        dec_type,
        dim,
        match,
        num_spkrs
        ) * decoding_matrix
    
    # expand to match input order (if necessary)
    if N_out != N_in:
        decoding_matrix = mmul(
            resize_order_matrix(N_in, N_out),
            decoding_matrix
        )
    
    return decoding_matrix