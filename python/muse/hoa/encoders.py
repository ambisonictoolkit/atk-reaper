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
from math import factorial

# import muse defined constants
import muse.constants as C



 #=========================
 # Definition of constants
 #=========================


#=========================
# Functions
#=========================

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
# m = index (order)
#
# returns unnormalised coefficients

# BUG: l, m = 1, 1 returns incorrect values! --> abs

def sph_harm(l = 0, m = 0, azimuth = 0., elevation = 0.):
    """Args:
        - l         : degree (classic ambisonic 'order')
        - m         : index (order)
        - azimuth   : counter-clockwise, in radians
        - elevation : upwards, in radians

    The harmonic $Y^l_m$ sampled at azimuth and elevation.

    Returns the un-normalised (real) spherical harmonic coefficient
    in the form required for Higher Order Ambisonics.
    """

    # for testing for and constructing vectors
    azim_scalar = isscalar(azimuth)
    elev_scalar = isscalar(elevation)
    
    if azim_scalar and not elev_scalar:
        azimuth = repeat(array([azimuth]), len(elevation))

    elif not azim_scalar and elev_scalar:
        elevation = repeat(array([elevation]), len(azimuth))


    # compute spherical harmonic
    a = pow(-1, m) * scipy.special.lpmv(
        abs(m),
        l,
        sin(elevation)
    )

    if m < 0:
        b = sin(abs(m) * azimuth)
    elif m == 0:
        b = 1
    elif m > 0:
        b = cos(abs(m) * azimuth)
    
    print 'azimuth coeff = ', b
    print 'elevation coeff = ', a
    
    res = a * b
    
    return res


#def mono_to_b(a, azimuth = 0., elevation = 0.):
#    """mono_to_b(a, azimuth = 0., elevation = 0.)
#    
#    Args:
#        - a         : Input mono signal
#        - azimuth   : counter-clockwise, in radians
#        - elevation : upwards, in radian
#
#    Encode a mono signal into the B-format domain.
#
#    """
#
#    # compute cosines and sines
#    cos_azim, sin_azim = cos(azimuth), sin(azimuth)
#    cos_elev, sin_elev = cos(elevation), sin(elevation)
#
#    # compute scalars
#    w_scale = C.rec_sqrt2
#    x_scale = cos_elev * cos_azim
#    y_scale = cos_elev * sin_azim
#    z_scale = sin_elev
#
#    # for testing for and constructing vectors
#    azim_scalar = isscalar(azimuth)
#    elev_scalar = isscalar(elevation)
#    n = len(a)
#
#    # construct appropriate encoder
#    if not elev_scalar:         # case 2, 4: x, y, z vectors
#        encoder = interleave(
#            array([
#                repeat(w_scale, n),
#                x_scale,
#                y_scale,
#                z_scale
#            ]))
#    elif not azim_scalar:       # case 3: x, y vectors
#        encoder = interleave(
#            array([
#                repeat(w_scale, n),
#                x_scale,
#                y_scale,
#                repeat(z_scale, n)
#            ]))
#    else:                       # case 1: all scalars
#        encoder = array([
#            w_scale,
#            x_scale,
#            y_scale,
#            z_scale
#            ])
#
#    # copy input mono array into four channel array
#    b = repeat(a, 4).reshape(-1,4)
#
#    # encode here!
#    return b * encoder


#---------------------------------------------
# A to B encoder
#---------------------------------------------






#-------------------------------------
# spherical harmonic test
#
# write up encoding using scipy's associated Legendre polynomials
#
# SEE: http://en.wikipedia.org/wiki/Ambisonic_data_exchange_formats
#
# AND:
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
# OK! This now works!! NOW, need to make broadcasting work correctly (yes?)
# to generate coefficients. (Following muse conventions.) AND add an
# HOA plane-wave encoder, where you specify the degree (ambisonic order).


#from math import factorial
#import pylab
#
#
## ---------------------------
## functions
#
## l = degree (classic ambisonic order)
## m = index (order)
##
## returns unnormalised coefficients
#def spher_harm(l, m, azimuth, elevation):
#    
#    a = pow(-1, m) * scipy.special.lpmv(
#        abs(m),
#        l,
#        sin(elevation)
#    )
#
#    if m < 0:
#        b = sin(abs(m) * azimuth)
#    elif m == 0:
#        b = 1
#    elif m > 0:
#        b = cos(abs(m) * azimuth)
#    
#    res = a * b
#    
#    return res
#
#
## Schmidt semi-normalisation (3D)
#def sn3d(l, m):
#    
#    if m == 0:
#        dm = 1
#    else:
#        dm = 0
#
#    res = sqrt(
#        ((2 - dm) / (4 * pi)) * \
#        factorial(l - abs(m)) / factorial(l + abs(m))
#    )
#    
#    return res
#
#
#def n3d(l, m):
#    
#    res = sqrt(2 * l + 1) * sn3d(l, m)
#    
#    return res
#
#
#
## ---------------------------
## parameters
#
#N = 2**8
#
#theta = twoPi * arange(N)/N
#phi = zeros(N)
#
##theta = 0
##phi = 0
#
#
## ---------------------------
## generate --> sn3d(l, m) * spher_harm(l, m, azimuth, elevation)
#
#
##b = interleave(array([
##    sn3d(0, 0) * spher_harm(0, 0, theta, phi),
##    sn3d(1, -1) * spher_harm(1, -1, theta, phi),
##    sn3d(1, 0) * spher_harm(1, 0, theta, phi),
##    sn3d(1, 1) * spher_harm(1, 1, theta, phi)
##]))
#
#b = interleave(array([
#    n3d(0, 0) / n3d(0, 0) * spher_harm(0, 0, theta, phi),
#    n3d(1, -1) / n3d(0, 0) * spher_harm(1, -1, theta, phi),
#    n3d(1, 0) / n3d(0, 0) * spher_harm(1, 0, theta, phi),
#    n3d(1, 1) / n3d(0, 0) * spher_harm(1, 1, theta, phi)
#]))
#
#
## ---------------------------
## plot
#pylab.plot(b)
#
#pylab.show()