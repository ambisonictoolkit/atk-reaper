#! /usr/bin/env python
# -*- coding: utf-8 -*-

# vim:syntax=python

# Joseph Anderson 2014


# Muse
# ==========

# This is the start of the Muse project, a Music V type language.

# Built on pyaudiolab.
# """


# seem to need to have these imports here. . . to make sure names are defined
from muse import *
#from filters import *
#from transforms import *

from muse.hoa import *

from numpy.polynomial import Chebyshev as T
from numpy.polynomial import Legendre as L
from numpy import math
#from math import factorial
from scipy.misc import factorial

from numpy.linalg.linalg import pinv


# import muse defined constants
import muse.constants as C  # may not need this here...



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
# Functions
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
def decoder_matching_gain(order = 1, dec_type = 'basic', dim = 3, \
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

def order_gains_matrix(order = 1, dec_type = 'basic', dim = 3):
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


def peri_to_panto_matrix(order):
    """
    Args:
        - order   : Ambisonic order, e.g., 1, 2, 3...


    Discard periphonic (3D) harmonics. No (cylindrical) weighting
    is applied.
    
    Returns matrix.
    """
    M = order
    
    nchans_in = (M+1)**2
    
    l, m = acn_to_lm(arange(nchans_in))
    
    res = identity(nchans_in)[((l == abs(m)))]

    return res


def panto_to_peri_matrix(order):
    """
    Args:
        - order   : Ambisonic order, e.g., 1, 2, 3...


    Insert periphonic (3D) harmonics, as zeros. No (spherical) weighting
    is applied.
    
    Returns matrix.
    """
    M = order
    
    res = transpose(peri_to_panto_matrix(M))

    return res


def zero_peri_matrix(order):
    """
    Args:
        - order   : Ambisonic order, e.g., 1, 2, 3...


    Zero periphonic (3D) harmonics. No (cylindrical) weighting
    is applied.
    
    Returns (square) matrix.
    """
    M = order
    
    l, m = acn_to_lm(arange((M+1)**2))
    
    res = diag((l == abs(m)).astype(int))

    return res


#=========================
# Decoding - Matricies (decoders)

def panto_pinv_decoding_matrix(directions, order = 1, dec_type = 'basic', \
    match = 'amp'):
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
    # prototype encoding matrix (('acn'), ('sn3d'))
    encoding_matrix = planewave_matrix(
        array([directions]),
        N_out
    )
        
    # discard 3D harmonics
    encoding_matrix = mmul(
        encoding_matrix,
        peri_to_panto_matrix(N_out)
    )
    
    # decoder: pseudo inverse - via pinv()
    decoding_matrix = pinv(encoding_matrix)
    
    # (re-)insert 3D harmonics
    decoding_matrix = mmul(
        peri_to_panto_matrix(N_out),
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
    match = 'amp'):
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
    # prototype encoding matrix (('acn'), ('sn3d'))
    encoding_matrix = planewave_matrix(
        array([directions]),
        N_out
    )
        
    # convert to n2d scaling - for projection
    encoding_matrix = mmul(
        encoding_matrix,
        encoding_convert_matrix(
            (('acn'), ('sn3d')),
            (('acn'), ('n2d')),
            N_out
        )
    )
    
    # zero 3D harmonics
    encoding_matrix = mmul(
        encoding_matrix,
        zero_peri_matrix(N_out)
    )
    
    # transpose and scale
    decoding_matrix = (transpose(encoding_matrix)) / num_spkrs
    
    # convert to n2d scaling - for input
    decoding_matrix = mmul(
        encoding_convert_matrix(
            (('acn'), ('sn3d')),
            (('acn'), ('n2d')),
            N_out
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