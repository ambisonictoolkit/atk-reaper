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
from muse import *
# from filters import *

# import muse defined constants
import constants as C



# #=========================
# # Definition of constants
# #=========================


#=========================
# Functions
#=========================


#=========================
# frequency independent
#=========================


#=========================
# axis primitives first
#=========================
def _rotate(a, theta = 0.):
    """_rotate(a, theta = 0.)
    
    Soundfield rotate primitive.
    
    Args:
        - a         : Input signal (two axes)
        - theta     : roataion angle, in radians

    See rotate.
    """

    # compute cosines and sines
    cos_theta, sin_theta = cos(theta), sin(theta)

    # construct appropriate transform
    transform = array([
        [cos_theta, -sin_theta],
        [sin_theta,  cos_theta]
        ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)

    return (a[:, newaxis, :] * transform).sum(axis = -1)


def _dominate(a, gain = 0.):
    """_dominate(a, gain = 0.)
    
    Soundfield dominance primitive.
    
    Args:
        - a         : Input b-format signal
        - gain      : dominance gain, in dB

    See dominate.
    """

    # compute gamma, dominance scale
    g = db_to_amp(gain)
    rec_g = 1./g

    k0 =          .5 * (g + rec_g)
    k1 = C.rec_sqrt2 * (g - rec_g)

    # construct appropriate transform
    transform = array([
        [k0, .5 * k1],
        [k1,      k0]
        ])

    if not(isscalar(gain)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)

    return (a[:, newaxis, :] * transform).sum(axis = -1)


# old zoom: used for graz poster
# def _zoom(a, theta = 0.):
#     """_zoom(a, theta = 0.)
    
#     Soundfield zoom primitive.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians

#     See zoom.
#     """

#     # compute gamma, dominance scale
#     k0 =      1. / cos(theta)
#     k1 = C.sqrt2 * tan(theta)

#     # construct appropriate transform
#     transform = array([
#             [k0, .5 * k1],
#             [k1,      k0]
#             ])

#     if not(isscalar(theta)):                       # reshape for vectors
#         transform = transform.transpose(2, 0, 1)

#     return (a[:, newaxis, :] * transform).sum(axis = -1)

# new zoom
def _zoom(a, theta = 0.):
    """_zoom(a, theta = 0.)
    
    Soundfield zoom primitive.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    See zoom.
    """

    # compute stuff
    k0 = sin(theta)
    k1 = cos(theta)

    # construct appropriate transform
    transform = array([
            [ones_like(theta), C.rec_sqrt2 * k0],
            [C.sqrt2 * k0, ones_like(theta)]
            ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)
        k1 = interleave(k1)

    return (a[:, newaxis, :] * transform).sum(axis = -1), k1 # return tuple



# def _focus(a, theta = 0.):
#     """_focus(a, theta = 0.)
    
#     Soundfield focus primitive.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians

#     See focus.
#     """

#     # compute gamma, dominance scale
#     sin_theta = sin(theta)
#     abs_sin_theta = abs(sin_theta)

#     k0 = 1. / (1. + abs_sin_theta)
#     k1 = C.sqrt2 * k0 * sin_theta
#     k2 = sqrt(k0 * (1. - abs_sin_theta))

#     # construct appropriate transform
#     transform = array([
#             [k0, .5 * k1],
#             [k1,      k0]
#             ])

#     if not(isscalar(theta)):                       # reshape for vectors
#         transform = transform.transpose(2, 0, 1)
#         k2 = interleave(k2)

#     return (a[:, newaxis, :] * transform).sum(axis = -1), k2 # return tuple
def _focus(a, theta = 0.):
    """_focus(a, theta = 0.)
    
    Soundfield focus primitive.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    See focus.
    """

    # compute gamma, dominance scale
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    abs_sin_theta = abs(sin_theta)

    k0 = 1. / (1. + abs_sin_theta)
    k1 = C.sqrt2 * k0 * sin_theta
    k2 = k0 * cos_theta

    # construct appropriate transform
    transform = array([
            [k0, .5 * k1],
            [k1,      k0]
            ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)
        k2 = interleave(k2)

    return (a[:, newaxis, :] * transform).sum(axis = -1), k2 # return tuple


def _push(a, theta = 0.):
    """_push(a, theta = 0.)
    
    Soundfield push primitive.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    See push.
    """

    # compute stuff
    sin_theta = sin(theta)

    k0 = C.sqrt2 * sin_theta * abs(sin_theta)
    k1 = pow(cos(theta), 2.)

    # construct appropriate transform
    transform = array([
            [ones_like(theta), zeros_like(theta)],
            [k0, k1]
            ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)
        k1 = interleave(k1)

    return (a[:, newaxis, :] * transform).sum(axis = -1), k1 # return tuple



def _press(a, theta = 0.):
    """_press(a, theta = 0.)
    
    Soundfield press primitive.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    See press
    """

    # compute stuff
    sin_theta = sin(theta)

    k0 = C.sqrt2 * sin_theta * abs(sin_theta)
    k1 = cos(theta)
    k2 = k1**2

    # construct appropriate transform
    transform = array([
            [ones_like(theta), zeros_like(theta)],
            [k0, k2]
            ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)
        k1 = interleave(k1)

    return (a[:, newaxis, :] * transform).sum(axis = -1), k1 # return tuple


# old press, depricated
# def _press(a, theta = 0.):
#     """_press(a, theta = 0.)
    
#     Soundfield press primitive.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians

#     See press
#     """

#     # compute stuff
#     sin_theta = sin(theta)
#     abs_sin_theta = abs(sin_theta)

#     k0 = C.sqrt2 * sin_theta
#     k1 = (1. - abs_sin_theta) / (1. + abs_sin_theta)
#     k2 = .5 * k0 * k1
#     k3 = cos(theta)

#     # construct appropriate transform
#     transform = array([
#             [ones_like(theta), k2],
#             [k0, k1]
#             ])

#     if not(isscalar(theta)):                       # reshape for vectors
#         transform = transform.transpose(2, 0, 1)
#         k3 = interleave(k3)

#     return (a[:, newaxis, :] * transform).sum(axis = -1), k3 # return tuple

def _asymmetry(a, theta = 0.):
    """_asymmetry(a, theta = 0.)
    
    Soundfield asymmetry primitive.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    See push.
    """

    # compute stuff
    sin_theta = sin(theta)

    k0 = C.sqrt2 * sin_theta * abs(sin_theta)
    k1 = pow(cos(theta), 2.)


    sin_theta = sin(theta)
    cos_theta = cos(theta)

    k0 = -1 * sin_theta
    k1 = pow(sin_theta, 2.)
    k2 = pow(cos_theta, 2.)
    k3 = cos_theta * sin_theta
    k4 = cos_theta

    # construct appropriate transform
    transform = array([
            [ones_like(theta), zeros_like(theta), C.rec_sqrt2 * k0],
            [C.sqrt2 * k1, k2, k0],
            [-1 * C.sqrt2 * k3, k3, k4]
            ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)
        k4 = interleave(k4)

    return (a[:, newaxis, :] * transform).sum(axis = -1), k4 # return tuple



#=========================
# mirroring
#=========================


def mirror_x(a):
    """mirror_x(a)
    
    Mirror an ambisonic B-format across the x-axis.
    
    Args:
        - a         : Input b-format signal

    """

    chnls = array([False, True, False, False]) # X

    res = a.copy()              # operate on a copy
    
    res[:, chnls] = -res[:, chnls]  # mirror

    return res


def mirror_y(a):
    """mirror_y(a)
    
    Mirror an ambisonic B-format across the y-axis.
    
    Args:
        - a         : Input b-format signal

    """

    chnls = array([False, False, True, False]) # Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls] = -res[:, chnls]  # mirror

    return res


def mirror_z(a):
    """mirror_z(a)
    
    Mirror an ambisonic B-format across the z-axis.
    
    Args:
        - a         : Input b-format signal

    """

    chnls = array([False, False, False, True]) # Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls] = -res[:, chnls]  # mirror

    return res


def mirror_o(a):
    """mirror_z(a)
    
    Mirror an ambisonic B-format across the origin.
    
    Args:
        - a         : Input b-format signal

    """

    chnls = array([False, True, True, True]) # X, Y, Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls] = -res[:, chnls]  # mirror

    return res


#=========================
# rotations
#=========================


def rotate(a, theta = 0.):
    """rotate(a, theta = 0.)
    
    Rotate an ambisonic B-format sound field counter-clockwise around the z-axis.
    
    (Rotation around z-axis)

    Args:
        - a         : Input b-format signal
        - theta     : roataion angle, in radians

    Theta = pi/2 rotates what was front center to hard left and
    -pi/2 rotates to hard right. Both pi and -pi will rotate what
    was front center to back center. The default, 0, results in no change.

    """

    chnls = array([False, True, True, False]) # X and Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls] =  _rotate(a[:, chnls], theta) # rotate channels

    return res


def tilt(a, theta = 0.):
    """tilt(a, theta = 0.)
    
    Tilt an ambisonic B-format sound field counter-clockwise around the x-axis.
    
    (Rotation around x-axis)

    Args:
        - a         : Input b-format signal
        - theta     : roataion angle, in radians

    Theta = pi/2 rotates what was hard left to up and -pi/2
    rotates to down. Both pi and -pi will rotate what was
    hard left to hard right, also resulting in a vertically
    inverted image. The default, 0, results in no change.

    """

    chnls = array([False, False, True, True]) # Y and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls] =  _rotate(a[:, chnls], theta) # rotate channels

    return res


def tumble(a, theta = 0.):
    """tumble(a, theta = 0.)
    
    Tumble an ambisonic B-format sound field counter-clockwise around the x-axis.
    
    (Rotation around y-axis)

    Args:
        - a         : Input b-format signal
        - theta     : roataion angle, in radians

    Theta = pi/2 rotates what was hard left to up and -pi/2
    rotates to down. Both pi and -pi will rotate what was
    hard left to hard right, also resulting in a vertically
    inverted image. The default, 0, results in no change.

    """

    chnls = array([False, True, False, True]) # X and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls] =  _rotate(a[:, chnls], theta) # rotate channels

    return res


#=========================
# directivity (and squishing!)
#=========================


def direct(a, theta = 0):
    """direct(a, theta = 0)
    
    Adjust directivity of an ambisonic B-format sound field.

    Args:
        - a         : Input b-format signal
        - theta     : directivity, in radians (-pi/2 to pi/2)


    Theta = 0 retains the current directivity of the soundfield.
    Increasing theta towards pi/2 decreases the directivity,
    reducing the gains on the directional compenents (x, y, z)
    to zero--and is equivalent to a 'spatial lo-pass filter'.
    The resulting image becomes 'directionless'.

    Decreasing theta towards -pi/2 decreases the gain on the w,
    the zero order component, and can be regarded as a kind of
    'spatial sharpening' filter.

    Standard use of direct() is with theta >=0,
    
    """

    # compute gains
    g0 = sqrt(1 + sin(theta))
    g1 = sqrt(1 - sin(theta))

    # construct appropriate transform
    transform = array([
            g0, 
            g1,
            g1,
            g1
        ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = interleave(transform)

    return a * transform


def direct_x(a, theta = 0):
    """direct_x(a, theta = 0)
    
    Adjust the directivy of an ambisonic B-format sound field on the x-axis.

    Args:
        - a         : Input b-format signal
        - theta     : directivity, in radians (-pi/2 to pi/2)


    Theta = 0 retains the current directivity of the soundfield.
    Increasing theta towards pi/2 decreases the directivity on
    the x-axis, reducing the gains on this axis to zero--and is
    equivalent to a 'spatial lo-pass filter'. The resulting image
    becomes 'directionless' on the x-axis.

    Decreasing theta towards -pi/2 decreases the gain on w, y, z,
    and can be regarded as a kind of 'spatial sharpening' filter
    on the x-axis.

    Standard use of direct_x() is with theta >=0,
    
    """

    # compute gains
    g0 = sqrt(1 + sin(theta))
    g1 = sqrt(1 - sin(theta))

    # construct appropriate transform
    transform = array([
            g0, 
            g1,
            g0,
            g0
        ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = interleave(transform)

    return a * transform


def direct_y(a, theta = 0):
    """direct_y(a, theta = 0)
    
    Adjust the directivy of an ambisonic B-format sound field on the y-axis.

    Args:
        - a         : Input b-format signal
        - theta     : directivity, in radians (-pi/2 to pi/2)


    Theta = 0 retains the current directivity of the soundfield.
    Increasing theta towards pi/2 decreases the directivity on
    the y-axis, reducing the gains on this axis to zero--and is
    equivalent to a 'spatial lo-pass filter'. The resulting image
    becomes 'directionless' on the y-axis.

    Decreasing theta towards -pi/2 decreases the gain on w, x, z,
    and can be regarded as a kind of 'spatial sharpening' filter
    on the y-axis.

    Standard use of direct_y() is with theta >=0,

    """

    # compute gains
    g0 = sqrt(1 + sin(theta))
    g1 = sqrt(1 - sin(theta))

    # construct appropriate transform
    transform = array([
            g0, 
            g0,
            g1,
            g0
        ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = interleave(transform)

    return a * transform


def direct_z(a, theta = 0):
    """direct_z(a, theta = 0)
    
    Adjust the directivy of an ambisonic B-format sound field on the z-axis.

    Args:
        - a         : Input b-format signal
       - theta     : directivity, in radians (-pi/2 to pi/2)


    Theta = 0 retains the current directivity of the soundfield.
    Increasing theta towards pi/2 decreases the directivity on
    the z-axis, reducing the gains on this axis to zero--and is
    equivalent to a 'spatial lo-pass filter'. The resulting image
    becomes 'directionless' on the z-axis.

    Decreasing theta towards -pi/2 decreases the gain on w, x, y,
    and can be regarded as a kind of 'spatial sharpening' filter
    on the z-axis.

    Standard use of direct_z() is with theta >=0,

    """

    # compute gains
    g0 = sqrt(1 + sin(theta))
    g1 = sqrt(1 - sin(theta))

    # construct appropriate transform
    transform = array([
            g0, 
            g0,
            g0,
            g1
        ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = interleave(transform)

    return a * transform


# may want to update this so it inculdes a matrix calc
# rather than acting as a wrapper
def direct_a(a, theta = 0, azimuth = 0., elevation = 0.):
    """direct_a(a, theta = 0, azimuth = 0., elevation = 0.)
    
    Adjust the directivy of an ambisonic B-format sound field along an axis
    oriented along azimuth, elevation.

    Args:
        - a         : Input b-format signal
       - theta     : directivity, in radians (-pi/2 to pi/2)
        - azimuth   : azimuth to direct along
        - elevation : elevation to direct along

    Theta = 0 retains the current directivity of the soundfield.
    Increasing theta towards pi/2 decreases the directivity on
    the axis determined by (azimuth, elevation), reducing the gains
    on this axis to zero--and is equivalent to a 'spatial lo-pass filter'.
    The resulting image becomes 'directionless' on this axis.

    Decreasing theta towards -pi/2 decreases the gain on w and the
    plane perpendicular to the normal defined by (azimuth, elevation)
    and can be regarded as a kind of 'spatial sharpening' filter
    on the specified axis.

    Standard use of direct_a() is with theta >=0,
    """

    # transform here!
    return rotate(
        tumble(
            direct_x(
                tumble(
                    rotate(
                        a,
                        -azimuth),
                    -elevation),
                theta),
            elevation),
        azimuth)

###=========================
### squishing
###=========================
##
##def squish_x(a, theta = pi/2.):
##    """squish_x(a, theta = pi/2.0)
##    
##    Squish an ambisonic B-format sound field on the x-axis.
##
##    Args:
##        - a         : Input b-format signal
##        - theta     : angle of distortion, in radians (-pi to pi)
##
##    Theta = 0 squishes the x-axis to the center, bringing what
##    was at front center and back center to the middle of the
##    image, squishing the image to the y-z plane. Values greater
##    than pi/2 increase the gain of x while reducing w, y, z;
##    this range of values can be used to exaggerate the front/back
##    depth of the image. Values between pi/2 and -pi/2 squish the
##    image from unchanged, through the y-z plane, to front/back
##    reversed. The default, pi/2, results in no change.
##
##    """
##
##    # compute cosines and sines
##    cos_theta, sin_theta = cos(.5 * theta), sin(.5 * theta)
##
##    # construct appropriate transform
##    transform = C.sqrt2 * array([
##            cos_theta, 
##            sin_theta,
##            cos_theta,
##            cos_theta
##        ])
##
##    if not(isscalar(theta)):                       # reshape for vectors
##        transform = interleave(transform)
##
##    return a * transform
##
##
##def squish_y(a, theta = pi/2.):
##    """squish_y(a, theta = pi/2.)
##    
##    Squish an ambisonic B-format sound field on the y-axis.
##
##    Args:
##        - a         : Input b-format signal
##        - theta     : angle of distortion, in radians (-pi to pi)
##
##    Theta = 0 squishes the y-axis to the center, bringing what
##    was at hard left and hard right to the middle of the image,
##    squishing the image to the x-z plane. Values greater than
##    pi/2 increase the gain of y while reducing w, x, z; this
##    range of values can be used to exaggerate the left/right
##    width of the image. Values between pi/2 and -pi/2 squish
##    the image from unchanged, through the x-z plane, to
##    left/right reversed. The default, pi/2, results in no change.
##
##    """
##
##    # compute cosines and sines
##    cos_theta, sin_theta = cos(.5 * theta), sin(.5 * theta)
##
##    # construct appropriate transform
##    transform = C.sqrt2 * array([
##            cos_theta, 
##            cos_theta,
##            sin_theta,
##            cos_theta
##        ])
##
##    if not(isscalar(theta)):                       # reshape for vectors
##        transform = interleave(transform)
##
##    return a * transform
##
##
##def squish_z(a, theta = pi/2.):
##    """squish_z(a, theta = pi/2.)
##    
##    Squish an ambisonic B-format sound field on the z-axis.
##
##    Args:
##        - a         : Input b-format signal
##        - theta     : angle of distortion, in radians (-pi to pi)
##
##    Theta = 0 squishes the z-axis to the center, bringing what
##    was at up and down to the middle of the image, squishing the
##    image to the x-y plane. Values greater than pi/2 increase the
##    gain of z while reducing w, x, y; this range of values can be
##    used to exaggerate the height of the image. Values between pi/2
##    and -pi/2 squish the image from unchanged, through the x-y
##    plane, to up/down inverted. The default, pi/2, results in no change.
##
##    """
##
##    # compute cosines and sines
##    cos_theta, sin_theta = cos(.5 * theta), sin(.5 * theta)
##
##    # construct appropriate transform
##    transform = C.sqrt2 * array([
##            cos_theta, 
##            cos_theta,
##            cos_theta,
##            sin_theta
##        ])
##
##    if not(isscalar(theta)):                       # reshape for vectors
##        transform = interleave(transform)
##
##    return a * transform
##
##
### may want to update this so it inculdes a matrix calc
### rather than acting as a wrapper
##def squish(a, theta = pi/2., azimuth = 0., elevation = 0.):
##    """squish(a, theta = pi/2., azimuth = 0., elevation = 0.)
##    
##    Squish an ambisonic B-format sound field along an axis
##    oriented along azimuth, elevation.
##
##    Args:
##        - a         : Input b-format signal
##        - theta     : angle of distortion, in radians (-pi to pi)
##        - azimuth   : azimuth to squish along
##        - elevation : elevation to squish along
##
##    Theta = 0 squishes the axis defined by (azimuth, elevation)
##    to the center, bringing what was at (azimuth, elevation)
##    and (-azimuth, -elevation) to the middle of the image,
##    squishing the image to the plane normal to (azimuth, elevation).
##    Values greater than pi/2 increase the gain of axis defined by
##    (azimuth, elevation) while reducing w, and the gain along the
##    plane normal to (azimuth, elevation); this range of values can
##    be used to exaggerate the depth of the image along the axis
##    defined by (azimuth, elevation). Values between pi/2 and -pi/2
##    squish the image from unchanged, through the plane normal to
##    (azimuth, elevation), to axis defined by (azimuth, elevation)
##    reversed. The default, pi/2, results in no change.
##
##    """
##
##    # transform here!
##    return rotate(
##        tumble(
##            squish_x(
##                tumble(
##                    rotate(
##                        a,
##                        -azimuth),
##                    -elevation),
##                theta),
##            elevation),
##        azimuth)
##


#=========================
# dominance, zoom, focus, bloom (not coded yet), push
#=========================


def dominate_x(a, gain = 0.):
    """dominate_x(a, gain = 0.)
    
    Apply dominance to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - gain      : dominance gain, in dB

    Gain: the dominance gain, in dB, applied on the x-axis. Positive values
    increase the gain at front center to +gain, while decreasing the gain at
    back center to -gain, simultaneously distorting the image towards front
    center. Negative values of gain invert this distortion, distorting the
    image towards back center. The default, 0, results in no change.

    """

    chnls = array([True, True, False, False]) # W and X

    res = a.copy()              # operate on a copy
    
    res[:, chnls] =  _dominate(a[:, chnls], gain) # dominate channels

    return res


def dominate_y(a, gain = 0.):
    """dominate_y(a, gain = 0.)
    
    Apply dominance to an ambisonic B-format sound field along the y-axis.
    
    Args:
        - a         : Input b-format signal
        - gain      : dominance gain, in dB

    Gain: the dominance gain, in dB, applied on the y-axis. Positive values
    increase the gain at hard left to +gain, while decreasing the gain at
    hard right to -gain, simultaneously distorting the image towards hard
    left. Negative values of gain invert this distortion, distorting the
    image towards hard right. The default, 0, results in no change.

    """

    chnls = array([True, False, True, False]) # W and Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls] =  _dominate(a[:, chnls], gain) # dominate channels

    return res


def dominate_z(a, gain = 0.):
    """dominate_z(a, gain = 0.)
    
    Apply dominance to an ambisonic B-format sound field along the z-axis.
    
    Args:
        - a         : Input b-format signal
        - gain      : dominance gain, in dB

    Gain: the dominance gain, in dB, applied on the z-axis. Positive values
    increase the gain at up to +gain, while decreasing the gain at down to
    -gain, simultaneously distorting the image towards up. Negative values
    of gain invert this distortion, distorting the image towards down. The
    default, 0, results in no change.

    """

    chnls = array([True, False, False, True]) # W and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls] =  _dominate(a[:, chnls], gain) # dominate channels

    return res


def dominate(a, gain = 0., azimuth = 0., elevation = 0.):
    """dominate(a, gain = 0., azimuth = 0., elevation = 0.)
    
    Apply dominance to an ambisonic B-format sound field along an axis
    oriented along azimuth, elevation.
    
    Args:
        - a         : Input b-format signal
        - gain      : dominance gain, in dB
        - azimuth   : azimuth to apply dominance along
        - elevation : elevation to apply dominance along


    Gain: the dominance gain, in dB, applied the axis defined by
    (azimuth, elevation). Positive values increase the gain at
    (azimuth, elevation) to +gain, while decreasing the gain at
    (-azimuth, -elevation) to -gain, simultaneously distorting
    the image towards (azimuth, elevation). Negative values of
    gain invert this distortion, distorting the image towards
    (-azimuth, -elevation). The default, 0, results in no change.

    """

    # transform here!
    return rotate(
        tumble(
            dominate_x(
                tumble(
                    rotate(
                        a,
                        -azimuth),
                    -elevation),
                gain),
            elevation),
        azimuth)


# def zoom_x(a, theta = 0.):
#     """zoom_x(a, theta = 0.)
    
#     Apply dominance, in terms of distortion angle, to an ambisonic B-format sound field along the x-axis.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians

#     Theta: the angle of distortion in radians, from -pi/2 to pi/2.
#     Positive values increase gain in the front center of the image,
#     reducing the gain at back center. Negative values do the inverse.
#     The default, 0, results in no change. Note, +/-pi/2 generate infinite
#     gain, so the usable range is between these values.

#     """

#     chnls = array([True, True, False, False]) # W and X

#     res = a.copy()              # operate on a copy
    
#     res[:, chnls] =  _zoom(a[:, chnls], theta) # zoom channels

#     return res


# def zoom_y(a, theta = 0.):
#     """zoom_y(a, theta = 0.)
    
#     Apply dominance, in terms of distortion angle, to an ambisonic B-format sound field along the y-axis.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians
    
#     Gain: the angle of distortion in radians, from -pi/2 to pi/2.
#     Positive values increase gain at the hard left of the image,
#     reducing the gain at hard right. Negative values do the inverse.
#     The default, 0, results in no change. Note, +/-pi/2 generate
#     infinite gain, so the usable range is between these values.

#     """

#     chnls = array([True, False, True, False]) # W and Y

#     res = a.copy()              # operate on a copy
    
#     res[:, chnls] =  _zoom(a[:, chnls], theta) # zoom channels

#     return res


# def zoom_z(a, theta = 0.):
#     """zoom_z(a, theta = 0.)
    
#     Apply dominance, in terms of distortion angle, to an ambisonic B-format sound field along the z-axis.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians
    
#     Gain: the angle of distortion in radians, from -pi/2 to pi/2. Positive
#     values increase gain at up of the image, reducing the gain at down.
#     Negative values do the inverse. The default, 0, results in no change.
#     Note, +/-pi/2 generate infinite gain, so the usable range is between
#     these values.

#     """

#     chnls = array([True, False, False, True]) # W and Z

#     res = a.copy()              # operate on a copy
    
#     res[:, chnls] =  _zoom(a[:, chnls], theta) # zoom channels

#     return res


# def zoom(a, theta = 0., azimuth = 0., elevation = 0.):
#     """zoom(a, theta = 0., azimuth = 0., elevation = 0.)
    
#     Apply dominance, in terms of distortion angle, to an ambisonic
#     B-format sound field an axis oriented along azimuth, elevation.
    
#     Args:
#         - a         : Input b-format signal
#         - theta     : angular distortion, in radians
#         - azimuth   : azimuth to apply zoom along
#         - elevation : elevation to apply zoom along

#     Theta: the angle of distortion in radians, from -pi/2 to pi/2.
#     Positive values increase the gain at (azimuth, elevation) to +gain,
#     while decreasing the gain at (-azimuth, -elevation) to -gain,
#     simultaneously distorting the image towards (azimuth, elevation).
#     Negative values of gain invert this distortion, distorting the
#     image towards (-azimuth, -elevation). The default, 0, results
#     in no change.

#     """

#     # transform here!
#     return rotate(
#         tumble(
#             zoom_x(
#                 tumble(
#                     rotate(
#                         a,
#                         -azimuth),
#                     -elevation),
#                 theta),
#             elevation),
#         azimuth)

def zoom_x(a, theta = 0.):
    """zoom_x(a, theta = 0.)
    
    Apply dominance, in terms of distortion angle, to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values increase gain in the front center of the image,
    reducing the gain at back center. Negative values do the inverse.
    The default, 0, results in no change.

    """

    chnls = array([True, True, False, False]) # W and X

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k1 =  _zoom(a[:, chnls], theta) # zoom channels
    res[:, invert(chnls)] =  k1 * res[:, invert(chnls)] # scale remaining channels

    return res


def zoom_y(a, theta = 0.):
    """zoom_y(a, theta = 0.)
    
    Apply dominance, in terms of distortion angle, to an ambisonic B-format sound field along the y-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians
    
    Gain: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values increase gain at the hard left of the image,
    reducing the gain at hard right. Negative values do the inverse.
    The default, 0, results in no change.

    """

    chnls = array([True, False, True, False]) # W and Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k1 =  _zoom(a[:, chnls], theta) # zoom channels
    res[:, invert(chnls)] =  k1 * res[:, invert(chnls)] # scale remaining channels

    return res


def zoom_z(a, theta = 0.):
    """zoom_z(a, theta = 0.)
    
    Apply dominance, in terms of distortion angle, to an ambisonic B-format sound field along the z-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians
    
    Gain: the angle of distortion in radians, from -pi/2 to pi/2. Positive
    values increase gain at up of the image, reducing the gain at down.
    Negative values do the inverse.

    """

    chnls = array([True, False, False, True]) # W and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k1 =  _zoom(a[:, chnls], theta) # zoom channels
    res[:, invert(chnls)] =  k1 * res[:, invert(chnls)] # scale remaining channels

    return res


def zoom(a, theta = 0., azimuth = 0., elevation = 0.):
    """zoom(a, theta = 0., azimuth = 0., elevation = 0.)
    
    Apply dominance, in terms of distortion angle, to an ambisonic
    B-format sound field an axis oriented along azimuth, elevation.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians
        - azimuth   : azimuth to apply zoom along
        - elevation : elevation to apply zoom along

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values increase the gain at (azimuth, elevation) to +gain,
    while decreasing the gain at (-azimuth, -elevation) to -gain,
    simultaneously distorting the image towards (azimuth, elevation).
    Negative values of gain invert this distortion, distorting the
    image towards (-azimuth, -elevation). The default, 0, results
    in no change.

    """

    # transform here!
    return rotate(
        tumble(
            zoom_x(
                tumble(
                    rotate(
                        a,
                        -azimuth),
                    -elevation),
                theta),
            elevation),
        azimuth)


def focus_x(a, theta = 0.):
    """focus_x(a, theta = 0.)
    
    Apply focus to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians. Positive values focus
    on the front center of the image, and at pi/2 collapse the soundfield
    to mono, reducing the gain at back center to -inf dB. Negative values
    focus on back center. The default, 0, results in no change.

    """
    chnls = array([True, True, False, False]) # W and X

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k2 =  _focus(a[:, chnls], theta) # focus channels
    res[:, invert(chnls)] =  k2 * res[:, invert(chnls)] # scale remaining channels

    return res


def focus_y(a, theta = 0.):
    """focus_y(a, theta = 0.)
    
    Apply focus to an ambisonic B-format sound field along the y-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values focus on the hard left of the image, and at pi/2
    collapse the soundfield to mono, reducing the gain at hard right
    to -inf dB. Negative values focus on hard left. The default, 0,
    results in no change.

    """

    chnls = array([True, False, True, False]) # W and Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k2 =  _focus(a[:, chnls], theta) # focus channels
    res[:, invert(chnls)] =  k2 * res[:, invert(chnls)] # scale remaining channels

    return res


def focus_z(a, theta = 0.):
    """focus_z(a, theta = 0.)
    
    Apply focus to an ambisonic B-format sound field along the z-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values focus on up of the image, and at pi/2 collapse
    the soundfield to mono, reducing the gain at down to -inf dB.
    Negative values focus on down. The default, 0, results in no change.

    """

    chnls = array([True, False, False, True]) # W and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k2 =  _focus(a[:, chnls], theta) # focus channels
    res[:, invert(chnls)] =  k2 * res[:, invert(chnls)] # scale remaining channels

    return res


def focus(a, theta = 0., azimuth = 0., elevation = 0.):
    """focus(a, theta = 0., azimuth = 0., elevation = 0.)
    
    Apply focus to an ambisonic B-format sound field on an axis
    oriented along azimuth, elevation.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians
        - azimuth   : azimuth to apply focus along
        - elevation : elevation to apply focus along

    Theta: the angle of distortion in radians. Positive values focus
    on (azimuth, elevation) of the image, and at pi/2 collapse the soundfield
    to mono, reducing the gain at (-azimuth, -elevation) to -inf dB.
    Negative values focus on (-azimuth, -elevation). The default, 0,
    results in no change.

    """

    # transform here!
    return rotate(
        tumble(
            focus_x(
                tumble(
                    rotate(
                        a,
                        -azimuth),
                    -elevation),
                theta),
            elevation),
        azimuth)


# also. . . add bloom, with parameter to control +/-

#=========================
# push, press
#=========================

def push_x(a, theta = 0.):
    """push_x(a, theta = 0.)
    
    Apply push to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values push to the front center of the image, and at
    pi/2 collapse the soundfield to mono. Negative values push to
    back center. The default, 0, results in no change.

    """
    chnls = array([True, True, False, False]) # W and X

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k1 =  _push(a[:, chnls], theta) # push channels
    res[:, invert(chnls)] =  k1 * res[:, invert(chnls)] # scale remaining channels

    return res


def push_y(a, theta = 0.):
    """push_y(a, theta = 0.)
    
    Apply push to an ambisonic B-format sound field along the y-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values push to hard left of the image, and at pi/2
    collapse the soundfield to mono. Negative values push to hard
    left. The default, 0, results in no change.

    """

    chnls = array([True, False, True, False]) # W and Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k1 =  _push(a[:, chnls], theta) # push channels
    res[:, invert(chnls)] =  k1 * res[:, invert(chnls)] # scale remaining channels

    return res


def push_z(a, theta = 0.):
    """push_z(a, theta = 0.)
    
    Apply push to an ambisonic B-format sound field along the z-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values push to up of the image, and at pi/2 collapse
    the soundfield to mono. Negative values push to down. The default,
    0, results in no change.

    """

    chnls = array([True, False, False, True]) # W and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k1 =  _push(a[:, chnls], theta) # push channels
    res[:, invert(chnls)] =  k1 * res[:, invert(chnls)] # scale remaining channels

    return res


def push(a, theta = 0., azimuth = 0., elevation = 0.):
    """push(a, theta = 0., azimuth = 0., elevation = 0.)
    
    Apply push to an ambisonic B-format sound field on an axis
    oriented along azimuth, elevation.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians
        - azimuth   : azimuth to apply push along
        - elevation : elevation to apply push along

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values push to (azimuth, elevation) of the image, and at
    pi/2 collapse the soundfield to mono. Negative values push to
    (-azimuth, -elevation). The default, 0, results in no change.

    """

    # transform here!
    return rotate(
        tumble(
            push_x(
                tumble(
                    rotate(
                        a,
                        -azimuth),
                    -elevation),
                theta),
            elevation),
        azimuth)


def press_x(a, theta = 0.):
    """press_x(a, theta = 0.)
    
    Apply press to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values press to the front center of the image, and at
    pi/2 collapse the soundfield to mono. Negative values press to
    back center. The default, 0, results in no change.

    """

    chnls = array([True, True, False, False]) # W and X

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k3 =  _press(a[:, chnls], theta) # press channels
    res[:, invert(chnls)] =  k3 * res[:, invert(chnls)] # scale remaining channels

    return res


def press_y(a, theta = 0.):
    """press_y(a, theta = 0.)
    
    Apply press to an ambisonic B-format sound field along the y-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values press to hard left of the image, and at pi/2
    collapse the soundfield to mono. Negative values press to hard
    left. The default, 0, results in no change.

    """

    chnls = array([True, False, True, False]) # W and Y

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k3 =  _press(a[:, chnls], theta) # press channels
    res[:, invert(chnls)] =  k3 * res[:, invert(chnls)] # scale remaining channels

    return res


def press_z(a, theta = 0.):
    """press_z(a, theta = 0.)
    
    Apply press to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values press to up of the image, and at pi/2 collapse
    the soundfield to mono. Negative values press to down. The default,
    0, results in no change.

    """

    chnls = array([True, False, False, True]) # W and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k2 =  _press(a[:, chnls], theta) # press channels
    res[:, invert(chnls)] =  k2 * res[:, invert(chnls)] # scale remaining channels

    return res


def press(a, theta = 0., azimuth = 0., elevation = 0.):
    """press(a, theta = 0., azimuth = 0., elevation = 0.)
    
    Apply press to an ambisonic B-format sound field on an axis
    oriented along azimuth, elevation.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians
        - azimuth   : azimuth to apply press along
        - elevation : elevation to apply press along

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values press to (azimuth, elevation) of the image, and at
    pi/2 collapse the soundfield to mono. Negative values press to
    (-azimuth, -elevation). The default, 0, results in no change.

    """

    # transform here!
    return rotate(
        tumble(
            press_x(
                tumble(
                    rotate(
                        a,
                        -azimuth),
                    -elevation),
                theta),
            elevation),
        azimuth)

#=========================
# asymmetry
#=========================

def asymmetry_x(a, theta = 0.):
    """asymmetry_x(a, theta = 0.)
    
    Apply asymmetry to an ambisonic B-format sound field along the x-axis.
    
    Args:
        - a         : Input b-format signal
        - theta     : angular distortion, in radians

    Theta: the angle of distortion in radians, from -pi/2 to pi/2.
    Positive values push -Y to the front center of the image, and at
    pi/2 collapse the soundfield to mono. Negative values push +Y to
    front center. The default, 0, results in no change.

    """
    chnls = array([True, True, True, False]) # W and Z

    res = a.copy()              # operate on a copy
    
    res[:, chnls], k4 =  _asymmetry(a[:, chnls], theta) # push channels
    res[:, invert(chnls)] =  k4 * res[:, invert(chnls)] # scale remaining channels

    return res


#=========================
# a_to_a
#=========================


def a_to_a(a, in_orientation = 'flu', out_orientation = 'flu'):
    """a_to_a(a, in_orientation = 'flu', out_orientation = 'flu')
    
    Args:
        - a               : Input A-format signal
        - in_orientation  : Orientation of the A-format channel tetrahedron
            flu = front left up:          FLU, FRD, BLD, BRU
            fld = front left down:        FLD, FRU, BLU, BRD
            flr = front left-right:       FL, FR, BU, BD
            fud = front up-down:          FU, FD, BL, BR
            fbd = front and back down:    F, BD, BLU, BRU 
            fbu = front and back up:      F, BU, BLD, BRD
            flru = front left-right up:   FLU, FRU, FD, B
            flrd = front left-right down: FLD, FRD, FU, B

        - out_orientation : Orientation of the A-format channel tetrahedron
            as above


    Reorient an A-format soundfield.

    """
    decoder_dict = C.b_to_a_dict

    # construct encoder
    transform = dot(
        array(decoder_dict[out_orientation]),
        array(decoder_dict[in_orientation]).transpose()
        )

    # transform here!
    return inner(a, transform)


#=========================
# frequency dependent
#=========================


##def distance(x, r, T, zi = None):
##    """distance(x, r, T, zi = None)
##    
##    "Distance compensation" filter an ambisonic B-format sound field.
##    (Inverse of "proximity filter".)
##    
##    (1st order highpass on X, Y, Z)
##
##    Args:
##        - x         : Input b-format signal
##        - r         : Distance, in meters, to compensate for
##        - T         : Sampling period, 1./sr
##        - zi        : Initial conditions for the filter delays.  A vector
##            (or array of vectors for an N-dimensional input) of length
##            max(len(a),len(b)).  If zi=None or is not given then initial
##            rest is assumed.  SEE signal.lfiltic for more information.
##    
##    Outputs: (y, {zf})
##    
##      y -- The "distance filtered" output.
##      zf -- If zi is None, this is not returned, otherwise, zf holds the
##            final filter delay values.
##    
##    Algorithm:
##      See butter and lfilter
##
##    """
##    # Calculate Wn
##    Wn = freq_to_Wn(C.speed_of_sound / (C.twoPi * r), T)
##
##    # generate b, a for 1st order highpass
##    b, a = butter(1, Wn, 'highpass')
##
##    # X, Y, Z mask
##    chnls = array([False, True, True, True])
##
##    # operate on a copy
##    res = x.copy()
##
##    # filter!
##    # Note: could use fiir_hp, but choose
##    # to generate coefficients directly to
##    # make the inverse filter more explicit
##    if zi is None:
##        res[:, chnls] =  ffilter(b, a, x[:, chnls])
##        
##        return res
##
##    else:
##        res[:, chnls], zf =  ffilter(b, a, x[:, chnls], zi)
##        
##        return res, zf



