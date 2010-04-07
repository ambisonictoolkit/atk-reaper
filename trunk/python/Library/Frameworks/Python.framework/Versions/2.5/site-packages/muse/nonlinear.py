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
# from generators import *
# from scipy.signal import *
# from scipy.signal.signaltools import *
from numpy import round as nround

# #=========================
# # Definition of constants
# #=========================


#=========================
# non-linear distortion
#=========================



def quant(x, qnsr):
    """quant(x, qnsr)

    Introduce quantization noise.
    
    Description
    
      Process a data sequence, x, by introducing quantization noise.
    
    Inputs:
    
      x -- An N-dimensional input array.
      qnsr -- Quantization-error-noise-to-signal ratio, in dB.
    
    Outputs: y
    
      y -- The output of the process.

      """
    k = db_to_amp(qnsr)
    recip_k = 1. / k


#     y = k * numpy.round(recip_k * x)
    y = k * nround(recip_k * x)
    
    return y
