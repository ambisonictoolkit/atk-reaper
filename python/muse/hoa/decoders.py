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

from muse import *

from numpy.polynomial import Chebyshev as T
from numpy.polynomial import Legendre as L
from numpy import math
#from math import factorial
from scipy.misc import factorial



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
