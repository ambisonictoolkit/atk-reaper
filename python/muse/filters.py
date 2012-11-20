#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
from generators import *
from sndfile import *          # uses sndfile rather than scikits.audiolab

import scipy.io as sio
import os
from scipy.signal import *
from scipy.signal.signaltools import *
from numpy.fft import fft, ifft, rfft, irfft
import scipy.fftpack

# import muse defined constants
import constants as C


# #=========================
# # Definition of constants
# #=========================


#=========================
# Functions
#=========================


#=========================
# FIR Filter Designers
#=========================

# should add IIR methods. . .

# ceptral methods
def rceps(b, min = -120.):
    """rceps(b, min = -120.)

    Return the real part of the cepstrum of kernel b.
    
    Inputs:
    
      b  -- signal (or filter) kernel
      min -- minimum dB value to clip amplitude response to.
               Reduces time aliasing.

    Outputs:
    
      b -- resulting kernel.
    """
    res = real(
        ifft(
            numpy.log(
                clip(
                    abs(fft(b)),
                    db_to_amp(min),
                    inf))))

    return res


def irceps(b):
    """irceps(b)

    Return the real part of the inverse cepstrum of kernel b.
    
    Inputs:
    
      b  -- signal (or filter) kernel

    Outputs:
    
      b -- resulting kernel.
    """
    res = real(
        ifft(
            numpy.exp(
                fft(b))
            ))

    return res


# Laguerre
def lagt(x, c, M):
    """lagt(x, c, M)

    Laguerre transform of input.
    
    Inputs:
    
      x -- input signal
      c -- Laguerre parameter (warping, -1 to +1, 0 = none)
      M -- number of terms of the Laguerre transform to compute.
    
    Outputs:
    
      y -- transformed output.


    Notes: 

    For for approximate length to match N:

    M = int((1 + abs(c)) / (1 - abs(c)) * N)


    To calculate c in terms of Wn:

    c = -(tan(pi / 2 * Wn) - 1.) / (tan(pi / 2 * Wn) + 1.)

    OR

    b, a = butter(1, Wn, 'lowpass')
    c = -a[1]

    """

    N = len(x)
    y = empty(M)                # empty output
    rx = flipud(x)              # time reverse x

    yy = ffilter(              # filter by normalizing filter lambda_0
        array([sqrt(1 - c**2)]),
        array([1., c]),
        x
        )
    y[0] = yy[-1]               # keep last sample for 1st output

    for k in range(1, M):        # allpass loop
        yy = ffilter(
            array([c, 1.]),
            array([1., c]),
            yy
            )
        y[k] = yy[-1]           # keep last sample for kth output

    return y


def lagtun(x, c, M):
    """lagtun(x, c, M)

    Laguerre transform of input, un-normalized.

    Inputs:
    
      x -- input signal
      c -- Laguerre parameter (warping, -1 to +1, 0 = none)
      M -- number of terms of the Laguerre transform to compute.
    
    Outputs:
    
      y -- transformed output.


    Notes: 

    For for approximate length to match N:

    M = int((1 + abs(c)) / (1 - abs(c)) * N)


    To calculate c in terms of Wn:

    c = -(tan(pi / 2 * Wn) - 1.) / (tan(pi / 2 * Wn) + 1.)

    OR

    b, a = butter(1, Wn, 'lowpass')
    c = -a[1]

    """

    N = len(x)
    y = empty(M)                # empty output
    rx = flipud(x)              # time reverse x

    yy = x                      # un-normalized input
    y[0] = yy[-1]               # keep last sample for 1st output

    for k in range(1, M):        # allpass loop
        yy = ffilter(
            array([c, 1.]),
            array([1., c]),
            yy
            )
        y[k] = yy[-1]           # keep last sample for kth output

    return y


# ************** New allpass normalised lagt
def lagtapn(x, c, M):
    """lagtapn(x, c, M)

    Laguerre transform of input, all-pass normalized.
    
    Inputs:
    
      x -- input signal
      c -- Laguerre parameter (warping, -1 to +1, 0 = none)
      M -- number of terms of the Laguerre transform to compute.
    
    Outputs:
    
      y -- transformed output.


    Notes: 

    For for approximate length to match N:

    M = int((1 + abs(c)) / (1 - abs(c)) * N)


    To calculate c in terms of Wn:

    c = -(tan(pi / 2 * Wn) - 1.) / (tan(pi / 2 * Wn) + 1.)

    OR

    b, a = butter(1, Wn, 'lowpass')
    c = -a[1]

    """

    N = len(x)
    y = empty(M)                # empty output
    rx = flipud(x)              # time reverse x

    yy = ffilter(              # filter by normalizing filter lambda_0
        array([sqrt(1 - c**2)]),
        array([1., c]),
        x
        )
    y[0] = yy[-1]               # keep last sample for 1st output

    for k in range(1, M):        # allpass loop
        yy = ffilter(
            array([c, 1.]),
            array([1., c]),
            yy
            )
        y[k] = yy[-1]           # keep last sample for kth output

    y_nor = lfilter(            # all-pass normalize the result
        array([1., -c]),
        array([sqrt(1 - c**2)]),
        y
        )

    return y_nor


# FIR methods. . . 

def reff(b, Wn):
    """reff(b, Wn)

    Return spectral reflection of kernel b about Wn.
    
    Inputs:
    
      b  -- filter kernel
      Wn -- center frequency of reflection (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)

    Outputs:
    
      b -- resulting kernel.
    """
    N = len(b)

    if N % 2 is 1:
        phi = (1 - N) * pi / 2 * Wn
        reff_b = 2 * fcososc(Wn, phi, N) * b

    else:
        phi = ((1 - N) * Wn + cos(pi / 2 * Wn)) * pi / 2
        reff_b = 2 * fsinosc(Wn, phi, N) * b

    return reff_b


def mirf(b):
    """mirf(b)

    Return spectral mirror of kernel b.
    (Reflection about Nyquist)
    """
    mir_b = .5 * reff(b, 1.)
    
    return mir_b


def invf(b):
    """invf(b)

    Return spectral inversion of a normalised linear phase kernel b.

    b + invf(b) returns allpass

    Note: invf is not optimal for even length b
    """
    N = len(b)

    ap_b = sinc(lin([-(N-1)/2., (N-1)/2.], N))
    ap_b *= hamming(N)

    inv_b = ap_b - b
    
    return inv_b

# # NOTE: FFT VERSION NEEDS REVISION, THERE IS A PROBLEM W/ PHASE
# #       IN THAT b + invf(b) != ap(b)
# def invf(b):
#     """invf(b)

#     Return spectral inversion of kernel b.

#     b + invf(b) returns allpass
#     """

#     # take real fft
#     fft_b = rfft(b)

#     # find magnitude
#     mag = abs(fft_b)

#     # find phase
#     phase = angle(fft_b)

#     # invert mag
#     inv_mag = 1. - mag

#     # take the ifft (real)
#     res = irfft(inv_mag * (cos(phase) + sin(phase) * 1j))
    
#     return res


def apf(b):
    """apf(b)

    Return allpass kernel, preserving phase of b.
    """
    
    # take real fft
    fft_b = rfft(b)

    # find phase
    phase = angle(fft_b)

    # take the ifft (real)
    res = irfft(cos(phase) + sin(phase) * 1j, len(b))

    return res


def linf(b):
    """linf(b)

    Return linear phase kernel, preserving magnitude of b.
    """

    # N, input/output kernel length
    N = len(b)
    
    # take real fft
    fft_b = rfft(b)

    # find magnitude
    mag = abs(fft_b)

    # set phase to linear
    M = len(mag)
    phase = lin([0, -pi * (M-1) * (N-1) / N], M)

    # take the ifft (real)
    res = irfft(mag * (cos(phase) + sin(phase) * 1j), N)

    return res


def norf(b):
    """norf(b)

    Return peak gain normalised kernel b.

    Normalised to 0dB.

    """
    N = len(b)

    res = reciprocal(
        peak(
            abs(
                rfft(b)
                )
            )
        ) * b

    return res


# minimum phase
# see: http://www.sfr-fresh.com/unix/privat/uade-2.09.tar.gz:a/uade-2.09/contrib/sinc-integral.py
# UADE (Unix Amiga Delitracker Emulator) plays old Amiga tunes with UAE emulation.
# Jens Schleusener 
# (T-Systems SfR)
# Bunsenstr. 10 
# D-37083 Gottingen 
# E-Mail: Jens.Schleusener@t-systems-sfr.com
# E-Mail: info@sfr-fresh.com

# also see. . .
# JOS:
# http://ccrma-www.stanford.edu/~jos/fp/Creating_Minimum_Phase_Filters.htmlhttp://ccrma-www.stanford.edu/~jos/fp/Creating_Minimum_Phase_Filters.html
# http://ccrma.stanford.edu/~jos/sasp/Minimum_Phase_Filter_Design.html
# http://ccrma.stanford.edu/~jos/fp/Matlab_Utilities.html
# http://ccrma.stanford.edu/~jos/sasp/Minimum_Phase_Causal_Cepstra.html#sec:laurent

# Cain:
# http://www.music.columbia.edu/pipermail/music-dsp/2004-February/059372.html

# see also:
# http://cnx.org/content/m12469/latest/

# for cepstrum see:
# http://www.dsprelated.com/showmessage/48073/1.php

def minf(b, min = -120., oversampling = 8):
    """minf(b, min = -120., oversampling = 8)

    Return a magnitude equivalent minimum phase kernel from b.
    
    Inputs:
    
      b  -- filter kernel
      min -- minimum dB value to clip amplitude response to.
               Reduces time aliasing.
      oversampling -- fft * size oversampling.
               Reduces time aliasing.

    Outputs:
    
      b -- resulting kernel.
    """
    n_b = len(b)                # length of input

    p = append(b,
               zeros(n_b * (oversampling - 1))
               )                # zero padded input

    n_p = len(p)                # length of padded

    # compute the real cepstrum
    x = rceps(p, min)

    # window the cepstrum so anticausal components are rejected
    w = zeros(n_p)
    w[0] = 1.
    w[1:n_p / 2] = 2.

    x *= w

    # take the inverse real cepstrum to return minimum phase
    res = irceps(x)[:n_b]

    return res


# kirkeby normalisation filter
# kernel is kernel to normalise
# Wns are LF and HF corner freqs defining band for normalisation
# roll_off is roll off in octaves from Wns (defines transition band)
# Es is the regularization parameter E(f)
# Note: kirkeby in the form described by Farina results in a 
#       normalised band pass filter between Wns
def kirkeby(b, Wns, roll_off = 1./3, Es = array([.01, 10.])):
    """kirkeby(b, Wns, roll_off = 1./3, Es = array([.01, 10.]))

    Return a Kirkeby normalisation filter kernel to normalise input kernel.
    
    Inputs:
    
      b         -- filter kernel (expected to be odd length)
      Wns       -- frequencies to normalise between, [Wn0, Wn1]
      roll_off  -- filter rolloff
      Es        -- Kirkeby regularisation parameter E(f)

    Outputs:
    
      resulting normalisation kernel.
    """
    
    N = nframes(b)         # kernel size
    M = N/2 + 1            # freqs (fft, has +-freqs)

    # Wn transition bands, for roll off
    c_Wns0 = empty(2)
    c_Wns0[0] = 2**(-roll_off/2) * Wns[0]
    c_Wns0[1] = 2**(roll_off/2) * Wns[0]

    c_Wns1 = empty(2)
    c_Wns1[0] = 2**(-roll_off/2) * Wns[1]
    c_Wns1[1] = 2**(roll_off/2) * Wns[1]

    # design regularization parameter E(f)
    kirkebyE = ones(M)

    fftWns = (scipy.fftpack.fftfreq(N) * 2)[:M] # only need zero + positive vals

    # test for mode...
    #          mode 1: [True, True], normalise all Wns
    #          mode 2: [False, False], normalise between Wns[0] and Wns[1]
    #          mode 3: [False, True], normalise above Wns[0] only
    #          mode 4: [True, False], normalise below Wns[1] only
    nanWns = isnan(Wns)

    if all(nanWns):             # mode 1
        kirkebyE = Es[0] * kirkebyE
        
    elif all(logical_not(nanWns)): # mode 2
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns0[0]:
                kirkebyE[n] = Es[1]

            elif c_Wns0[0] <= k_Wn and k_Wn < c_Wns0[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns0[0]) * pi) / (c_Wns0[1] - c_Wns0[0])) + 1) / 2 + Es[0]

            elif c_Wns0[1] <= k_Wn and k_Wn < c_Wns1[0]:
                kirkebyE[n] = Es[0]

            elif c_Wns1[0] <= k_Wn and k_Wn < c_Wns1[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns1[0]) * pi) / (c_Wns1[1] - c_Wns1[0]) + pi) \
                         + 1) / 2 + Es[0]

            else:               # c_Wns1[1] <= k_Wn:
                kirkebyE[n] = Es[1]

    elif logical_not(nanWns[0]):    # mode 3
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns0[0]:
                kirkebyE[n] = Es[1]

            elif c_Wns0[0] <= k_Wn and k_Wn < c_Wns0[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns0[0]) * pi) / (c_Wns0[1] - c_Wns0[0])) + 1) / 2 + Es[0]

            else:               # c_Wns0[1] <= k_Wn
                kirkebyE[n] = Es[0]

    else:                       # mode 4
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns1[0]:
                kirkebyE[n] = Es[0]

            elif c_Wns1[0] <= k_Wn and k_Wn < c_Wns1[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns1[0]) * pi) / (c_Wns1[1] - c_Wns1[0]) + pi) \
                         + 1) / 2 + Es[0]

            else:               # c_Wns1[1] <= k_Wn
                kirkebyE[n] = Es[1]

    kirkebyE = concatenate((    # complete building kirkebyE by mirroring
            kirkebyE,
            kirkebyE[::-1][:-1]
            ))

    # compute kirkeby "packing" (normalisation) filter

    # 1) take fft, transform to freq domain
    fftH = scipy.fftpack.fft(b)

    # 2) compute inverse filter in freq domain
    conjH = conjugate(fftH)
    fftK = conjH / (conjH * fftH + kirkebyE)

    # 3) take ifft, transform back to time domain (and window)
    res = real(scipy.fftpack.ifft(fftK))
    res *= hann(N)

    return res


# *************************
# filter kernel generation functions, below

# consider replacing FIR design methods with sinc methods
# see also: http://www.nicholson.com/rhn/dsp.html
def fir_lp(N, Wn, width=pi):
    """fir_lp(N, Wn, width=pi)

    Lowpass FIR Filter Design using windowed ideal filter method, with the
    Kaiser window.
    
    Inputs:
    
      N  -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    return firwin(N, Wn, window = ('kaiser', width))


def fir_hp(N, Wn, width=pi):
    """fir_hp(N, Wn, width=pi)

    Highpass FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N  -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    lp_b = fir_lp(N, Wn, width)
    hp_b = invf(lp_b)

    return hp_b


def fir_ls(N, Wn, k, width=pi):
    """fir_ls(N, Wn, k, width=pi)

    Low shelf FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N  -- order of filter (number of taps)

      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at low frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    lp_b = fir_lp(N, Wn, width)
    hp_b = invf(lp_b)

    ls_b = k * lp_b + hp_b

    return ls_b


def fir_hs(N, Wn, k, width=pi):
    """fir_hs(N, Wn, k, width=pi)

    High shelf FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N  -- order of filter (number of taps)

      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at high frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    lp_b = fir_lp(N, Wn, width)
    hp_b = invf(lp_b)

    hs_b = lp_b + k * hp_b

    return hs_b


def fir_bp(N, Wn, bw, width=pi):
    """fir_bp(N, Wn, bw, width=pi)

    Bandpass FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- center frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    if N % 2 is 1:
        lp_b = fir_lp(N, bw / 2, width)
        
        bp_b = reff(lp_b, Wn)
        
    else:
        bp_b = fir_lp(N, Wn + bw / 2, width) - fir_lp(N, Wn - bw / 2, width)

    return bp_b


def fir_bs(N, Wn, bw, width=pi):
    """fir_bs(N, Wn, bw, width=pi)

    Bandstop FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    if N % 2 is 1:
        bp_b = fir_bp(N, Wn, bw, width)
        
        bs_b = invf(bp_b)
        
    else:
        bs_b = fir_hp(N, Wn + bw / 2, width) + fir_lp(N, Wn - bw / 2, width)

    return bs_b


def fir_pk(N, Wn, bw, k, width=pi):
    """fir_pk(N, Wn, bw, k, width=pi)

    Peaking FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at peaking frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bp_b = fir_bp(N, Wn, bw, width)
    bs_b = fir_bs(N, Wn, bw, width)

    pk_b = k * bp_b + bs_b

    return pk_b


def fir_sk(N, Wn, bw, k, width=pi):
    """fir_sk(N, Wn, bw, k, width=pi)

    Skirting FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at peaking frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bp_b = fir_bp(N, Wn, bw, width)
    bs_b = fir_bs(N, Wn, bw, width)

    sk_b = bp_b + k * bs_b

    return sk_b


def fir_bps(N, Wn, bw, k, width=pi):
    """fir_bps(N, Wn, bw, k, width=pi)

    Bandpass bank FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- an array of center frequencies of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      bw -- an array of cutoff bandwidths of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      k -- an array of scales at center frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bps_b = zeros(N)

    for n in range(len(Wn)):
        bps_b += k[n] * fir_bp(N, Wn[n], bw[n], width)

    return bps_b


def fir_pks(N, Wn, bw, k, width=pi):
    """fir_pks(N, Wn, bw, k, width=pi)

    Peaking bank FIR Filter Design using windowed ideal filter method.

    Note: this algorithm doesn't perform well for even orders.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- an array of center frequencies of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      bw -- an array of cutoff bandwidths of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      k -- an array of scales at center frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bps_b = fir_bps(N, Wn, bw, 1. - k, width)

    pks_b = invf(bps_b)

    return pks_b


# fn to generate hilbert coeficients
# **************************************
# Gibson, D. (1996, April). "Desigining an SSB Outphaser, Part 1."
# Electronics World, 306-310.
# Gibson, D. (1996, May). "Desigining an SSB Outphaser, Part 2."
# Electronics World, 392-394.
# **************************************

def fir_hb(N, beta=5):
    """fir_hb(N, beta=5)

    Hilbert FIR Filter Design:

        N = odd, using method demonstrated by Gibson
        N = even, using Nyquist mirrored sinc()

    Windowed with the Kaiser window.

    Gibson, D. (1996, April). "Desigining an SSB Outphaser, Part 1."
    Electronics World, 306-310.
    Gibson, D. (1996, May). "Desigining an SSB Outphaser, Part 2."
    Electronics World, 392-394.
    
    Inputs:
    
      N     : order of filter (number of taps), should be odd
      beta  : beta for Kaiser window FIR design.
                  5 = similiar to Hamming.
    
    Outputs:
    
      b     : coefficients of length N FIR filter, returned as complex
              coefficients. Unit response is found in 'real' and Hilbert
              response is found in 'imag'.

    """

    # real response
    x_real = sinc(lin([-(N-1)/2., (N-1)/2.], N))

    # imag response
    if N % 2 is 1:                          # N odd
        x_imag = zeros(N)

        for i in range(N):
            if i == (N -1) / 2:
                pass
            else:
                x_imag[i] = (1 - (cos((i - (N - 1) / 2) * pi))) / \
                            ((i - (N - 1) / 2) * pi)

    else:                                   # N even
        x_imag = mirf(x_real)

    # sum real and imag (cast as complex)
    res = x_real + 1j * x_imag

    # window
    res *= kaiser(N, beta)

    return res


def fir_ap(N, width = pi):
    """fir_ap(N, width = pi)

    Allpass FIR Filter design using frequency sampling.
    
    Inputs:
    
      N         -- order of filter (number of taps)
      width     -- random phase width
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """

    M = N/2 + 1

    phase = white(M)
    phase *= width

    # convert to 'linear' phase
    phase += lin([0, -pi * (M-1) * (N-1) / N], M)

    # normalise phase at DC
    phase[0] = 0

    # normalise phase at Nyquist
    if mod(N, 2) == 0:
        phase[-1] = 0.
    else:
        phase[-1] = pi/2

    # take the ifft (real)
    ap_b = irfft(cos(phase) + sin(phase) * 1j, N)

    return ap_b


##def fir_ap(N, Wn=.5, width=pi):
##    """fir_ap(N, Wn=.5, width=pi)
##
##    Allpass FIR Filter Design using ideal filter method, un-windowed.
##    
##    Inputs:
##    
##      N  -- order of filter (number of taps)
##      Wn -- frequency warping
##      width -- phase range
##    
##    Outputs:
##    
##      b      -- coefficients of length N FIR filter.
##    """
##
##    Nr = N/2 + 1
##
##    mag = ones(Nr)
##    phase = white(Nr)
##    
##    if Wn != .5:                # warp
##
##        x0 = lin(nframes = Nr)      # linear index
##
##        phase = interp(
##            Wn_warp(x0, 1. - Wn, True),
##            x0,
##            phase
##            )
##
##    phase *= width
##    phase += lin([0., -pi * Nr], Nr)
##    phase[0] = 0
##
##    # take the ifft (real)
##    ap_b = irfft(mag * (cos(phase) + sin(phase) * 1j))
##
##    return ap_b
##

##def fir_ap(N, width=pi):
##    """fir_ap(N, width=pi)
##
##    Allpass FIR Filter Design using ideal filter method, un-windowed.
##    
##    Inputs:
##    
##      N  -- order of filter (number of taps)
##      width -- phase range
##    
##    Outputs:
##    
##      b      -- coefficients of length N FIR filter.
##    """
##
##    Nr = N/2 + 1
##
##    mag = ones(Nr)
##    phase = white(Nr)
##
##    phase *= width
##    phase += lin([0., -pi * Nr], Nr)
##    phase[0] = 0
##
##    # take the ifft (real)
##    ap_b = irfft(mag * (cos(phase) + sin(phase) * 1j))
##
##    return ap_b


### NOTE: we may want to revisit this in light of recent
###       comb-warping [S(PACE)] developed methods
##def fir_apw(N, Wn=.5, width=pi, minp=True):
##    """fir_apw(N, Wn=.5, width=pi, minp=True)
##
##    Allpass FIR Filter Design using ideal filter method, warped and windowed.
##    
##    Inputs:
##    
##      N  -- order of filter (number of taps)
##      Wn -- frequency warping
##      width -- kaiser window beta
##      minp -- minimum phase?
##    
##    Outputs:
##    
##      b      -- coefficients of length N FIR filter.
##    """
##
##    N /= 2
##
##    c = -(tan(pi / 2 * Wn) - 1.) / (tan(pi / 2 * Wn) + 1.)
##    Nr = int((1 - abs(c)) / (1 + abs(c)) * N)/2 + 1
##
##    phase = white(Nr)
##    phase[0] = 0
##
##    # take the ifft (real): generates ap filter
##    ap_b = irfft(cos(phase) + sin(phase) * 1j)
##
##    # warp: allpass normalised
##    ap_b = lagtapn(ap_b, c, N)
##
##    # window: reduce 'tail'
##    ap_b *= kaiser(2 * N, width)[N:]
##
##    # re-normalise gain for allpass: through convolution with inverse filter
##    rfft_b = rfft(ap_b)
##    mag = 1./absolute(rfft_b)
##    phase = lin([0., -pi * len(mag)], len(mag))
##    inv_b = irfft(mag * (cos(phase) + sin(phase) * 1j))
##    if minp:
##        inv_b = minf(inv_b)
##    ap_b = convfilt(ap_b, inv_b, 'full')
##
##    # append zeros as necessary
##    if len(ap_b) < 2 * N:
##        ap_b = concatenate((ap_b, zeros(2 * N - len(ap_b))))
##
##    return ap_b



# New all pass warping
def fir_lap(M, Wn, N = None, width = pi, phase = pi):
    """fir_lap(M, Wn, N, width=pi)

    Laguerre warped allpass FIR Filter Design using ideal filter method,
        un-windowed.
    
    Inputs:
    
        M -- order of filter (number of terms of the Laguerre
              transform to compute.)
        Wn -- Laguerre warping frequency. Wn = 0.5 gives no warping.
              See notes below
        N  -- order of FFT. If is None, computes as 'approximate length'
              FFT. See notes below.
        width     -- random phase width
        phase -- kernel phase. pi centers the response, pi/2 moves
                    response to 1/4 of the window
    
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.


    Notes: 

    For for approximate length to match N:

        N = int((1 - abs(c)) / (1 + abs(c)) * M)


    To calculate c in terms of Wn:

        c = -(tan(pi / 2 * Wn) - 1.) / (tan(pi / 2 * Wn) + 1.)

        OR

        b, a = butter(1, Wn, 'lowpass')
        c = -a[1]
    """

    # Laguerre warping coefficient
    b, a = butter(1, Wn, 'lowpass')
    c = -a[1]

    # compute N, if not given
    if N is None:
        N = int((1 - abs(c)) / (1 + abs(c)) * M)

    # warp
    ap_b = lagtapn(
            fir_ap(N, width),
            c,
            M
            )

    # restore response to allpass, and adjust 'kernel phase'
    ap_b = roll(
        apf(ap_b),
        int(M * phase / C.twoPi)
        )

    return ap_b


#=========================
# Convolution Function
#=========================

def fftpackrconvolve(in1, in2, fftN):
    # returns complete fft

    res = real(
        scipy.fftpack.ifft(
            scipy.fftpack.fft(in1, fftN) * scipy.fftpack.fft(in2, fftN)
            )
        )

    return res


# NOTE: fconvolve is principally to be used as a primitive
# doesn't handle all cases well--use convfilt! (this is a primitive for convfilt)
def fconvolve(in1, in2):
    """fconvolve(in1, in2)

    Convolve two N-dimensional arrays, using fftpack.
    
    Description:
    
       Convolve in1 and in2.
    
    Inputs:
    
      in1 -- an N-dimensional array.
      in2 -- an array with the same number of dimensions as in1.
    
    Outputs:  (out,)
    
      out -- an N-dimensional array containing a subset of the discrete linear
             convolution of in1 with in2. The output is the full discrete linear
             convolution of the inputs.
    """
    
    # size and channels for in1 and in2
    in1N = nframes(in1)
    in2N = nframes(in2)

    in1Ch = nchannels(in1)
    in2Ch = nchannels(in2)

    # resulting convolution length
    cN = in1N + in2N - 1

    # resulting fft size should be a power of 2 for speed
    fftN = 2**int(ceil(log2(cN)))

    # case 2: in1 is multichannel and in2 is single channel (interleaved)
    if (in1Ch > 1) and (in2Ch is 1):
        in2 = tile(
            in2,
            in1Ch
            )

    # case 3: in1 is single channel (interleaved) and in2 is multichannel
    elif (in1Ch is 1) and (in2Ch > 1):
        in1 = tile(
            in1,
            in2Ch
            )

    if in1Ch is 1 and in2Ch is 1:
        res = fftpackrconvolve(in1, in2, fftN)[:cN]

    else:
        maxCh = max(in1Ch, in2Ch)
        res = empty([maxCh, in1N + in2N - 1]) # deinterleaved, empty res

        in1_d = deinterleave(in1)
        in2_d = deinterleave(in2)

        for n in range(maxCh):
            res[n] = fftpackrconvolve(in1_d[n], in2_d[n], fftN)[:cN]
        res = interleave(res)

    return res


# may want to add a convolution function that
# takes kernel as a spectrum
def convfilt(x, kernel, mode = 'z', kind = 'fft', zi = None):
    """convfilt(x, kernel, mode = 'z', zi = None)

    Convolve input x by kernel. 
    
    Description
    
      Convolve input x by kernel along the 0 axis.
      Operates in two different modes, returning the complete
      convolution, or acting as a filter with state.

      Wraps convolve and fftconvolve.
    
    Inputs:
    
      x -- A N-dimensional input array.
      kernel -- A one or N-dimensional input array.
      mode -- 'z' or 'full'. If mode is 'z', acts as a filter
               with state 'z', and returns a vector of length len(x).
               If mode is 'full', returns the full convolution.
      kind -- 'direct' or 'fft', for direct or fft convolution
      zi -- Initial state.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            len(kernel) - 1.  If zi=None or is not given then initial
            rest is assumed.
    
    Outputs: (y, {zf})
    
      y -- The output of the delay.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            delay values.
    
      """
    if kind is 'direct':
        convfun = convolve
    else:
        convfun = fconvolve

    x_chans = nchannels(x)
    k_chans = nchannels(kernel)

    # case 1: x and kernel are are both single channel
    if (x_chans is 1) and (k_chans is 1):
        y = convfun(x, kernel)

    # case 2: x is multichannel and kernel is single channel
    elif (x_chans > 1) and (k_chans is 1):
        kernel = interleave(kernel)
        y = convfun(x, kernel)

    # case 3: x is single channel and kernel is multichannel
    elif (x_chans is 1) and (k_chans > 1):
        x = interleave(x)
        y = convfun(x, kernel)

    # case 4: x and kernel are multichannel, test if equal
    elif x_chans == k_chans:

        y = empty([x_chans, nframes(x) + nframes(kernel) - 1]) # deinterleaved, empty y

        x_d = deinterleave(x)
        k_d = deinterleave(kernel)

        for n in range(x_chans):
            y[n] = convfun(x_d[n], k_d[n])
        y = interleave(y)

    else:                       # raise error here
        raise ValueError, ("Doh!!, x and kernel don't broadcast!")

    # now return the result. . .
    if mode is 'z':

        y, zf = split(y, [nframes(x)])

        if zi is None:
            return y

        else:
            if nframes(zi) == (nframes(kernel) - 1):
                over_dub(y, zi, write_over = True)
                return y, zf

            else:
                print nframes(zi), (nframes(kernel) - 1)
                raise ValueError, ("Doh!!, zi is wrong size!!")

    else:
        return y

# # may want to add a convolution function that
# # takes kernel as a spectrum
# def convfilt(x, kernel, mode = 'z', kind = 'fft', zi = None):
#     """convfilt(x, kernel, mode = 'z', zi = None)

#     Convolve input x by kernel. 
    
#     Description
    
#       Convolve input x by kernel along the 0 axis.
#       Operates in two different modes, returning the complete
#       convolution, or acting as a filter with state.

#       Wraps convolve and fftconvolve.
    
#     Inputs:
    
#       x -- A N-dimensional input array.
#       kernel -- A one or N-dimensional input array.
#       mode -- 'z' or 'full'. If mode is 'z', acts as a filter
#                with state 'z', and returns a vector of length len(x).
#                If mode is 'full', returns the full convolution.
#       kind -- 'direct' or 'fft', for direct or fft convolution
#       zi -- Initial state.  It is a vector
#             (or array of vectors for an N-dimensional input) of length
#             len(kernel) - 1.  If zi=None or is not given then initial
#             rest is assumed.
    
#     Outputs: (y, {zf})
    
#       y -- The output of the delay.
#       zf -- If zi is None, this is not returned, otherwise, zf holds the
#             delay values.
    
#       """
#     if kind is 'direct':
#         convfun = convolve
#     else:
#         convfun = fftconvolve

#     x_chans = nchannels(x)
#     k_chans = nchannels(kernel)

#     # case 1: x and kernel are are both single channel
#     if (x_chans is 1) and (k_chans is 1):
#         y = convfun(x, kernel)

#     # case 2: x is multichannel and kernel is single channel
#     elif (x_chans > 1) and (k_chans is 1):
#         kernel = interleave(kernel)
#         y = convfun(x, kernel)

#     # case 3: x is single channel and kernel is multichannel
#     elif (x_chans is 1) and (k_chans > 1):
#         x = interleave(x)
#         y = convfun(x, kernel)

#     # case 4: x and kernel are multichannel, test if equal
#     elif x_chans == k_chans:

#         y = empty([x_chans, nframes(x) + nframes(kernel) - 1]) # deinterleaved, empty y

#         x_d = deinterleave(x)
#         k_d = deinterleave(kernel)

#         for n in range(x_chans):
#             y[n] = convfun(x_d[n], k_d[n])
#         y = interleave(y)

#     else:                       # raise error here
#         raise ValueError, ("Doh!!, x and kernel don't broadcast!")

#     # now return the result. . .
#     if mode is 'z':

#         y, zf = split(y, [nframes(x)])

#         if zi is None:
#             return y

#         else:
#             if nframes(zi) == (nframes(kernel) - 1):
#                 over_dub(y, zi, write_over = True)
#                 return y, zf

#             else:
#                 print nframes(zi), (nframes(kernel) - 1)
#                 raise ValueError, ("Doh!!, zi is wrong size!!")

#     else:
#         return y


#=========================
# IIR Filter Functions
#=========================


# consider adding:
#    ffilter_ar - accepts multi, channel b, a (following appropriate broadcast rules)
#    vfilter_ar - time varying b, a (built on vfilter_ar)
def ffilter(b, a, x, zi=None):
    """ffilter(b, a, x, zi=None)

    Filter data along one-dimension with an IIR or FIR filter.
    (ffilter is an lfilter wrapper to correct for 2-dimensional inputs.)
    N-dim >= 2. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using a digital filter.  This works for many
      fundamental data types (including Object type).  The filter is a direct
      form II transposed implementation of the standard difference equation
        (see "Algorithm").
    
    Inputs:
    
      b -- The numerator coefficient vector in a 1-D sequence.
      a -- The denominator coefficient vector in a 1-D sequence.  If a[0]
           is not 1, then both a and b are normalized by a[0].
      x -- An N-dimensional input array.
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See lfilter

      """
    axis = 0                    # filter along the 0 axis

    if (x.ndim is 2) and (zi != None): # test for 2-dim case where lfilter fails

        chans = nchannels(x)
        
        x = deinterleave(x, False) # deinterleave x

        y = empty_like(x) # empty outputs
        zf = empty_like(zi)

        for (n, x_n, zi_n) in zip(range(chans), x, zi): # iterate, by channel

            y[n], zf[n] = lfilter(b, a, x_n, axis, zi_n)

        y = interleave(y, False)   # deinterleave y

        return y, zf

    else:
        return lfilter(b, a, x, axis, zi)


# create a convenience / private function to do section filtering
# then ffos and fsos can use the section function
# **possibly use c-style indexing for coefs and states  

def _f_section_cascade(b, a, x, section_order = 1, zi=None):
    """_f_section_cascade(b, a, x, section_order = 1, zi=None)

    Filter data along one-dimension with an IIR or FIR filter.
    N-dim >= 2. Filter along the 0 axis

    Helper function to create a cascade of first or second order filters.
    See ffos or fsos.
        
    Algorithm:
      See ffilter & lfilter

      """
    # catch number of sections a == b
    if len(b) != len(a):
        raise ValueError, ("len(b) must equal len(a)")
    else:
        nos = len(b) / (section_order + 1) # number of sections

    # for convenience reshape b, a for iteration
    b.shape = (nos, section_order + 1)
    a.shape = (nos, section_order + 1)
    
    y = x.copy()             # set output to input, for cascading loop

    if zi is None:              # no zi
        for (b_n, a_n) in zip(b, a): # cascade
            y = ffilter(b_n, a_n, y) # filters

        return y
    else:
        if x.ndim is 1:
            zi.shape = (nos, section_order) # sections, section order
        else:
            zi = swapaxes(
                reshape(
                    zi,
                    (nos, shape(x)[1], section_order)), # nos, channels, section order
                0, 1)

        zf = zeros_like(zi)

        for (b_n, a_n, zi_n, n) in zip(b, a, zi, range(nos)): # cascade

            y, zf[n] = ffilter(b_n, a_n, y, zi_n)      # filters

        if x.ndim is 1:
            zf.shape = (section_order * nos,)
        else:
            zf = reshape(
                swapaxes(
                    zf, 0, 1),
                (nos, section_order * shape(x)[1])) # nos, order * channels

        return y, zf


def ffos(b, a, x, zi=None):
    """ffos(b, a, x, zi=None)

    Filter data along one-dimension with an IIR or FIR filter.
    N-dim >= 2. Filter along the 0 axis

    ffos is a cascade of first order sections.
    
    Description
    
      Filter a data sequence, x, using a digital filter.  This works for many
      fundamental data types (including Object type).  The filter is a direct
      form II transposed implementation of the standard difference equation
        (see "Algorithm").
    
    Inputs:
    
      b -- The collection of numerator coefficients in a 1-D sequence.
      a -- The collection of denominator coefficients in a 1-D sequence.
           len(a) must == len(b)
      x -- An N-dimensional input array.
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See ffilter & lfilter

      """
    return _f_section_cascade(b, a, x, 1, zi)


def fsos(b, a, x, zi=None):
    """fsos(b, a, x, zi=None)

    Filter data along one-dimension with an IIR or FIR filter.
    N-dim >= 2. Filter along the 0 axis

    fsos is a cascade of second order sections.
    
    Description
    
      Filter a data sequence, x, using a digital filter.  This works for many
      fundamental data types (including Object type).  The filter is a direct
      form II transposed implementation of the standard difference equation
        (see "Algorithm").
    
    Inputs:
    
      b -- The collection of numerator coefficients in a 1-D sequence.
      a -- The collection of denominator coefficients in a 1-D sequence.
           len(a) must == len(b)
      x -- An N-dimensional input array.
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See ffilter & lfilter

      """
    return _f_section_cascade(b, a, x, 2, zi)


# EXPECT NEED TO READDRESS THIS FUNCTION
def vfilter(b, a, x, zi=None):
    """vfilter(b, a, x, zi=None)

    Filter data along one-dimension with an IIR or FIR filter.
    vfilter is a time varying version of ffilter.
    N-dim >= 2. Filter along the 0 axis.
    
    Description
    
      Filter a data sequence, x, using a digital filter.  This works for many
      fundamental data types (including Object type).  The filter is a direct
      form II transposed implementation of the standard difference equation
      (see "Algorithm").

    Inputs:
    
      b -- The numerator coefficient vector in a 1-D sequence.
      a -- The denominator coefficient vector in a 1-D sequence.  If a[0]
           is not 1, then both a and b are normalized by a[0].
           Both b and a are time varying arrays of length order * sample length.
      x -- An N-dimensional input array.
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(shape(b)[1], shape(a)[1]).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See lfilter

      """
    # find order
    N = max(shape(b)[1], shape(a)[1]) - 1

    # find nframes
    nframes = shape(x)[0]

    # set initial/final state to zeros
    if zi is None:

        # set initial/final state to zeros
        if x.ndim is 1:
            zf = zeros(N)

        else:
            # find channels - only works for = 2-D input
            channels = shape(x)[1]
            zf = zeros((channels, N))

        zi_set = False

    # otherwise set to zi
    else:
        zf = zi
        zi_set = True

    y = empty_like(x)   # set empty y

    # iterate by frame
    for (n, b_n, a_n, x_n) in zip(range(nframes), b, a, x): # iterate, by sample frame
        y[n], zf = ffilter(b_n, a_n, array([x_n]), zf)

    y.shape = shape(x)          # this line is necessary for the 1-d case

    if zi_set:
        return y, zf
    else:
        return y


# **************************************
# IIR Filter Functions
# **************************************

def fiir_lp(x, N, Wn, zi = None):
    """fiir_lp(x, N, Wn, zi = None)

    Filter data along one-dimension with a fixed frequency butterworth
    low pass IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using a low pass digital butterworth filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      N -- order
      Wn -- cutoff
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    b, a = butter(N, Wn, 'lowpass')

    y = ffilter(b, a, x, zi)
    
    return y


def fiir_hp(x, N, Wn, zi = None):
    """fiir_hp(x, N, Wn, zi = None)

    Filter data along one-dimension with a fixed frequency butterworth
    high pass IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using a high pass digital butterworth filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      N -- order
      Wn -- cutoff
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    b, a = butter(N, Wn, 'highpass')

    y = ffilter(b, a, x, zi)
    
    return y


def fiir_bp(x, N, Wn, q = 1., zi = None):
    """fiir_bp(x, N, Wn, q = 1., zi = None)

    Filter data along one-dimension with a fixed frequency butterworth
    band pass IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using a band pass digital butterworth filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      N -- order
      Wn -- center Wn
      q -- filter Q
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    bw = Wn / q

    Wn = array([Wn - .5 * bw, Wn + .5 * bw])

    b, a = butter(N, Wn, 'bandpass')

    y = ffilter(b, a, x, zi)
    
    return y


def fiir_bs(x, N, Wn, q = 1., zi = None):
    """fiir_bs(x, N, Wn, q = 1., zi = None)

    Filter data along one-dimension with a fixed frequency butterworth
    band stop IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using a band stop digital butterworth filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      N -- order
      Wn -- center Wn
      q -- filter Q
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    bw = Wn / q

    Wn = array([Wn - .5 * bw, Wn + .5 * bw])

    b, a = butter(N, Wn, 'bandstop')

    y = ffilter(b, a, x, zi)
    
    return y


def fiir_ap(x, N, Wn, zi = None):
    """fiir_ap(x, N, Wn, zi = None)

    Filter data along one-dimension with a fixed frequency butterworth
    allpass IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using an all pass digital butterworth filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      N -- order
      Wn -- cutoff
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    b, a = butter(N, Wn, 'lowpass')
    b = flipud(a)

    y = ffilter(b, a, x, zi)
    
    return y


def diff_filt(x, zi = None):
    """diff_filt(x, zi = None)

    Filter data along one-dimension with a differentiator
    high pass IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using a differentiator high pass filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    b = array([1., -1.])
    a = array([1.])

    y = ffilter(b, a, x, zi)
    
    return y


def integ_filt(x, zi = None):
    """integ_filt(x, zi = None)

    Filter data along one-dimension with an integrator
    ow pass IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using an integrator low pass filter.
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    b = array([1.])             # integrator coeffs
    a = array([1., -1.])

    y = ffilter(b, a, x, zi)
    
    return y


# --------------------------------------
# Regalia-Mitra filters
# --------------------------------------

def fiir_rmhs2(x, Wn, k, zi = None):
    """fiir_rmhs2(x, Wn, k, zi = None)

    Filter data along one-dimension with a fixed frequency Regalia-Mitra
    2nd-order high-shelf IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using Regalia-Mitra 2nd-order high-shelf
      filter. This filter has a constant phase for all k, matching the
      phase of a 1st-order all-pass filter.
      
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      Wn -- cutoff
      k -- scale at high frequencies
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    # NOTE: for simplicity of handling of state (zi, zf) this
    #       filter is implemented directly as a 2nd order filter
    #       rather than through Regalia-Mitra architecture

    # generate allpass coefficients
    b, a = butter(1, Wn, 'lowpass')

    c = a[1]

    # generate 2nd-order shelf coefficients
    b = array([
        ((1.-k)/4 * (1+c**2)) + ((1.+k)/2 * c),
        ((1.-k) * c) + ((1.+k)/2 * (1.+c**2)),
        ((1.-k)/4 * (1+c**2)) + ((1.+k)/2 * c),
    ])
    a = array([
        1,
        2*c,
        c**2
    ])

    y = ffilter(b, a, x, zi)

    return y


def fiir_rmls2(x, Wn, k, zi = None):
    """fiir_rmls2(x, Wn, k, zi = None)

    Filter data along one-dimension with a fixed frequency Regalia-Mitra
    2nd-order low-shelf IIR filter. Filter along the 0 axis
    
    Description
    
      Filter a data sequence, x, using Regalia-Mitra 2nd-order low-shelf
      filter. This filter has a constant phase for all k, matching the
      phase of a 1st-order all-pass filter.
      
      The filter is a direct form II transposed implementation of the standard
      difference equation
      (see "Algorithm").
    
    Inputs:
    
      x -- An N-dimensional input array.
      Wn -- cutoff
      k -- scale at low frequencies
      zi -- Initial conditions for the filter delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

      """
    # NOTE: for simplicity of handling of state (zi, zf) this
    #       filter is implemented directly as a 2nd order filter
    #       rather than through Regalia-Mitra architecture

    # generate allpass coefficients
    b, a = butter(1, Wn, 'lowpass')

    c = a[1]

    # generate 2nd-order shelf coefficients
    b = array([
        ((k-1.)/4 * (1+c**2)) + ((k+1.)/2 * c),
        ((k-1.) * c) + ((k+1.)/2 * (1.+c**2)),
        ((k-1.)/4 * (1+c**2)) + ((k+1.)/2 * c),
    ])
    a = array([
        1,
        2*c,
        c**2
    ])

    y = ffilter(b, a, x, zi)

    return y


# --------------------------------------
# Ambisonic / Soundfield filters
# --------------------------------------

def nfc(x, r, T, zi = None):
    """nfc(x, r, T, zi = None)
    
    "Near-field distance compensation" filter an ambisonic B-format sound field.
    (Inverse of "proximity filter".)
    
    (1st order highpass on X, Y, Z)

    Args:
        - x         : Input b-format signal
        - r         : Distance, in meters, to compensate for
        - T         : Sampling period, 1./sr
        - zi        : Initial conditions for the filter delays.  A vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The near-field compensated output.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

    """
    # Calculate Wn
    Wn = freq_to_Wn(C.speed_of_sound / (C.twoPi * r), T)

    # generate b, a for 1st order highpass
    b, a = butter(1, Wn, 'highpass')

    # X, Y, Z mask
    chnls = array([False, True, True, True])

    # operate on a copy
    res = x.copy()

    # filter!
    if zi is None:
        res[:, chnls] =  ffilter(b, a, x[:, chnls])
        
        return res

    else:
        res[:, chnls], zf =  ffilter(b, a, x[:, chnls], zi)
        
        return res, zf


def proximity(x, r, T, zi = None):
    """proximity(x, r, T, zi = None)
    
    "Proximity filter" an ambisonic B-format sound field.
    (Inverse of "distance filter".)
    
    (Integrate, sum on X, Y, Z)

    Args:
        - x         : Input b-format signal
        - r         : Distance, in meters, to generate cues for
        - T         : Sampling period, 1./sr
        - zi        : Initial conditions for the filter delays.  A vector
            (or array of vectors for an N-dimensional input) of length
            max(len(a),len(b)).  If zi=None or is not given then initial
            rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The "proximity" filtered output.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See butter and lfilter

    NOTE: As X, Y, Z are integrated, it is best to initially prepare
          the signal with a highpass filter.

    """
    # Calculate Wn
    Wn = freq_to_Wn(C.speed_of_sound / (C.twoPi * r), T)

    # generate b, a for 1st order highpass
    b, a = butter(1, Wn, 'highpass')

    # X, Y, Z mask
    chnls = array([False, True, True, True])

    # operate on a copy
    res = x.copy()

    # filter!
    # b, a coefficients are inverted
    # see distance filter above
    if zi is None:
        res[:, chnls] =  ffilter(a, b, x[:, chnls])
        
        return res

    else:
        res[:, chnls], zf =  ffilter(a, b, x[:, chnls], zi)
        
        return res, zf


def psycho_shelf(x, Wn, k, zi = None):
    """psycho_shelf(x, Wn, k, zi = None)

    Psychoacoustic shelf filtering for Ambisonic Dual Band decoding.
        
    Args:
        x           : input b-format signal
        Wn          : shelf filter corner, set to ~400 Hz
        k           : an array [k0, k1] where k0 is the W scale at high
                        frequencies, and k1 is the X, Y, Z scale at high
                        frequencies
        zi          : Initial conditions for the filter delays, a
                        4 x 2 vector. (4 channels of 2nd order) If
                        zi = None or is not given then initial rest is
                        assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})

        y           : The output of the digital filter.
        zf          : If zi is None, this is not returned, otherwise, zf
                        holds the final filter delay values.
    
    Algorithm:
      See fiir_rmhs2, ffilter and lfilter

    """
    # channel masks
    deg = array([
        [True, False, False, False],    # W (degree 0) mask
        [False, True, True, True]       # X, Y, Z (degree 1) mask
        ])

    # operate on a copy
    y = x.copy()

    # filter!
    if zi is None:
        for i in range(size(k)):
            y[:, deg[i]] = fiir_rmhs2(x[:, deg[i]], Wn, k[i])
        return y

    else:
        zf = zeros_like(zi)             # init zf

        for i in range(size(k)):
            y[:, deg[i]], zf[deg[i]] = \
                 fiir_rmhs2(x[:, deg[i]], Wn, k[i], zi[deg[i]])
        return y, zf
        

# --------------------------------------
# Spherical and HRIR filters
# --------------------------------------

def sphere(N, theta, T, \
           r = 0.0875, k_0 = 0.1, theta_0 = 5./6 * pi, width = pi):
    """sphere(x, N, theta, T, \
           r = 0.0875, k_0 = 0.1, theta_0 = 5./6 * pi, width = pi)

    Spherical model FIR Filter Design using method demonstrated by Brown
    and Duda, windowed with the Kaiser window.

    Brown, C. Phillip and Richard O. Duda. "An Efficient HRTF Model for
    3-D Sound." In proceedings: IEEE Workshop on Applications of Signal
    Processing to Audio and Acoustics (WASPAA), New Paltz, NY 1997.
    
    Args:
    
        - N         : order of filter (number of taps)
        - theta     : incidence angle (-pi, pi)
        - T         : sampling period, 1./sr 
        - r         : sphere radius (default is Brown/Duda value)
        - k_0       : filter minimum scale
        - theta_0   : incidence angle for minimimum scale
        - width     : beta for Kaiser window FIR design.
                      pi = minimum ripple for steepest cutoff.    
    Outputs:
    
        - b         : coefficients of length N FIR filter.

    """

    # precondition angle, constrain to -pi, pi
    if theta > pi:
        theta -= C.twoPi
    elif theta > pi:
        theta += C.twoPi

    # compute delay from centre of sphere
    if abs(theta) < C.halfPi:
        delay = -r / C.speed_of_sound * cos(theta)
    else:
        delay = r / C.speed_of_sound * (abs(theta) - C.halfPi)

    delayN = delay * reciprocal(T)


    # compute filter coefficients
    w_0 = C.speed_of_sound / r * T
    k = (1 + k_0 / 2) + \
              (1 - k_0 / 2) * cos(theta / theta_0 * pi)

    b = array([w_0 + k, w_0 - k])
    a = array([w_0 + 1, w_0 - 1])


    # generate delay response
    x =  sinc(
        lin([-(N-1)/2., (N-1)/2.], N) - delayN
        )

    # apply frequency dependent gain, and window
    x = lfilter(b, a, x)
    x *= kaiser(N, width)

    return x


def open_sphere(N, theta, T, \
           r = 0.13125, k = 0., width = pi):
    """open_sphere(N, theta, T, \
           r = 0.13125, k = 0., width = pi)

    Open spherical model FIR Filter Design windowed with the Kaiser window.
    
    Args:
    
        - N         : order of filter (number of taps)
        - theta     : incidence angle
        - T         : sampling period, 1./sr 
        - r         : sphere radius
        - k         : Microphone response pattern, k. k = 1 returns
                        a bi-directional pattern (fig-8). k = 0.5
                        returns a cardioid pattern. I.e.:

                        (1-k) + k * cos
        - width     : beta for Kaiser window FIR design.
                      pi = minimum ripple for steepest cutoff.    
    Outputs:
    
      b      -- coefficients of length N FIR filter.

    """

    # compute delay from centre of sphere
    delay = -r / C.speed_of_sound * cos(theta)
    delayN = delay * reciprocal(T)

    # compute directional scale
    scale = (1.-k) + (k * cos(theta))

    # generate delay response
    x =  sinc(
        lin([-(N-1)/2., (N-1)/2.], N) - delayN
        )

    # apply gain, and window
    x *= scale
    x *= kaiser(N, width)

    return x


def sHRIR(N, azimuth, elevation, T, \
           r = 0.0875, theta_e = 5./9*pi, width = pi):
    """sHRIR(N, azimuth, elevation, T, \
           r = 0.0875, theta_e = 5./9*pi, width = pi)

    Spherical HRIR model FIR Filter Design using method demonstrated by Brown
    and Duda, windowed with the Kaiser window.

    Brown, C. Phillip and Richard O. Duda. "An Efficient HRTF Model for
    3-D Sound." In proceedings: IEEE Workshop on Applications of Signal
    Processing to Audio and Acoustics (WASPAA), New Paltz, NY 1997.
    
    Args:
    
        - N         : order of filter (number of taps)
        - azimuth   : source azimuth (-pi, pi)
        - elevation : sourze elevation (-pi/2, pi/2)
        - T         : sampling period, 1./sr 
        - r         : sphere radius (default is Brown/Duda value)
        - theta_e   : +/- ear angle (default is Brown/Duda value)
        - width     : beta for Kaiser window FIR design.
                      pi = minimum ripple for steepest cutoff.    
    Outputs:
    
        - b         : coefficients of length N FIR filter. Interleaved
                        [left_FIR, right_FIR]

    Note: As as simple sphereical model, the generated HRIR does not
            include torso or pinnae effects. Because of this, the
            simple spherical model does not generate elevation cues.

    """

    # calculate angles to ears
    theta_l = arccos(cos(azimuth - theta_e) * cos(elevation))
    theta_r = arccos(cos(azimuth + theta_e) * cos(elevation))

    # calculate sHRIR
    res = interleave(
        array([
            sphere(N, theta_l, T, r, width = width),
            sphere(N, theta_r, T, r, width = width)
            ])
        )

    return res

# NOTE: Consider adding general open and closed sphere arrays.
#       These would be useful for testing and simulation of
#       prototype real arrays before deployment.


def lHRIR(azimuth, elevation, subject_id, database_dir, status = 'C'):
    """lHRIR(azimuth, elevation, subject_id, database_dir, status = 'C')

    Return measured HRIRs from the IRCAM hosted Listen HRTF database.
    HRIR measurements were taken in blocked-meatus conditions.

    See: http://recherche.ircam.fr/equipes/salles/listen/
    

    Args:
    
        - azimuth       : source azimuth (-pi, pi)
        - elevation     : sourze elevation (-pi/2, pi/2)
        - subject_id    : 1002 - 1059 (as string)
        - database_dir  : path to local directory where subject directories
                          for the Listen database are located                          
        - status        : compensated ('C') or raw ('R') HRIRs
                            len('C' ) = 512, len('R') = 8192

    Outputs:
    
        - b         : coefficients of measured FIR filter. Interleaved
                        [left_FIR, right_FIR]


    The following table shows HRIR measurement points :

    Elevation (deg)         Azimuth increment (deg)     Points per elevation
    ------------------------------------------------------------------------
    -45                     15                          24
    -30                     15                          24
    -15                     15                          24
      0                     15                          24
     15                     15                          24
     30                     15                          24
     45                     15                          24
     60                     30                          12
     75                     60                           6
     90                     360                          1


    The setup consists of 10 elevation angles starting at -45deg ending at
    +90deg in 15deg steps vertical resolution. The steps per rotation varies
    from 24 to only 1 (90deg elevation). Measurement points are always
    located at the 15deg grid, but with increasing elevation only every
    second or fourth measurement point is taken into account. As a whole,
    there are 187 measurement points, hence 187 stereo audio files.


    Note: Returns the complete (asymmetric) HRIR

    """

    # lHRIR dirs, constants
    status_dir = {
        'C' : '/COMPENSATED/MAT/HRIR',
        'R' : '/RAW/MAT//HRIR'
        }
    status_key = {
        'C' : 'eq_hrir_S',
        'R' : 'hrir_S'
        }

    # generate lHRIR path
    HRIR_file = os.path.join(
        database_dir, 'IRC_' + subject_id + status_dir[status], \
        'IRC_' + subject_id + '_' + status + '_HRIR.mat'
        )

    # Read matlab file:
    HRIR_matfile = sio.loadmat(HRIR_file)

    # Read measured azimuths and elevations (from left)...
    azs = HRIR_matfile['l_' + status_key[status]]['azim_v'][0][0]
    els = HRIR_matfile['l_' + status_key[status]]['elev_v'][0][0]

    # convert azimuth, elevation to degrees
    # ... and determine index of desired HRIR
    az = int(round(rad_to_deg(azimuth))) % 360
    el = int(round(rad_to_deg(elevation))) % 360

    index = argmax(equal(az, azs) & equal(el, els))

    # retrieve HRIR
    res = interleave(array([
        HRIR_matfile['l_' + status_key[status]]['content_m'][0][0][index],
        HRIR_matfile['r_' + status_key[status]]['content_m'][0][0][index]
        ]))

    return res


def cHRIR(azimuth, elevation, subject_id, database_dir):
    """cHRIR(azimuth, elevation, subject_id, database_dir)

    Return measured HRIRs from the CIPIC HRTF database.
    HRIR measurements were taken in blocked-meatus conditions.

    See: http://interface.cipic.ucdavis.edu/sound/hrtf.html
    

    Args:
    
        - azimuth       : source azimuth (-pi, pi)
        - elevation     : sourze elevation (-pi/2, pi/2)
        - subject_id    : 003 - 165 (as string)
                          subjects '021' and '165' are the KEMAR head
        - database_dir  : path to local directory where subject directories
                          for the CIPIC database are located                          

    Outputs:
    
        - b         : coefficients of measured FIR filter. Interleaved
                        [left_FIR, right_FIR]


    Please see "Documentation for the UCD HRIR Files" avaliable at the
    above link for measurement details. As a whole, there are 1250
    measurement points.


    Note: Returns the complete (asymmetric) HRIR

    """

    # constants: cHRIR azimuth and elevation (in degrees)
    c_azimuths = concatenate((linspace(0, 45, 10), array([55, 65, 80])))
    c_azimuths = concatenate((-c_azimuths[1:][::-1], c_azimuths))
    c_elevations = -45 + (360. / 64 * arange(50))

    # map from ambisonic to CIPIC coordinates (in degrees)
    c_azimuth = round(
        rad_to_deg(
            -arcsin(sin(azimuth) * cos(elevation))
            ),
        3
        )
    c_elevation = round(
        rad_to_deg(
            arctan2(tan(elevation), cos(azimuth))
            ),
        3
        ) % 360
    if c_elevation > 230.625:
        c_elevation -= 360

    # generate az, el (indices) from azimuth, elevation
    # az = 12                 # 0 deg
    # el = 8                  # 0 deg
    az = where(c_azimuths == c_azimuth)[0][0]
    el = where(c_elevations == c_elevation)[0][0]


    # generate cHRIR path
    HRIR_file = os.path.join(
        database_dir, 'subject_' + subject_id, \
        'hrir_final.mat'
        )

    # Read matlab file:
    HRIR_matfile = sio.loadmat(HRIR_file)

    res = interleave(array([
        HRIR_matfile['hrir_l'][az, el],
        HRIR_matfile['hrir_r'][az, el]
        ]))


    # pad with zeros to nearest power of 2
    res = concatenate((
        res,
        zeros([2**int(ceil(log2(nframes(res)))) - nframes(res), nchannels(res)])
        ))

    return res


# --------------------------------------
# RIR filters
# --------------------------------------


def isophonicsRIR(subject_id, database_dir):
    """isophonicsRIR(subject_id, database_dir)

    Return measured RIRs from the Isophonics RIR database.

    See: http://isophonics.net/content/room-impulse-response-data-set
    

    Args:
    
        - subject_id    : 'xXXyYY', e.g., 'x00y00' (as string)
                          
        - database_dir  : path to local directory where subject directories
                          for the Isophonics database are located
                          (include desired space, e.g. 'greathall', 'octagon',
                          'classroom')

    Outputs:
    
        - b         : coefficients of measured FIR filter. Interleaved
                        [W_FIR, X_FIR, Y_FIR, Z_FIR]


    Please see documentation avaliable at the above link for measurement
    details. Measurement points vary for each space.


    Note: Returns the complete (asymmetric) HRIR

    """

    # harmonics
    harms = [ 'W', 'X', 'Y', 'Z' ]


    # generate RIR paths
    RIR_files = []
    for harm in harms:
        RIR_files += [os.path.join(
            database_dir, harm, \
            harm + subject_id + '.wav'
            )]

    # read in
    res = []        # start as list
    for RIR_file in RIR_files:
        RIR_sndfile = Sndfile(RIR_file)
        res += [RIR_sndfile.read_frames(RIR_sndfile.nframes)]
        RIR_sndfile.close()

    # convert to array and interleave

    res = interleave(array(res))

    return res



# **************************************
# osc. . .
# **************************************

def phasor(Wn, phase = 0., zi = None):
    """phasor(Wn, phase = 0., zi = None)

    Args:
        - Wn    : Normalized frequency array (may be multichannel)
        - phase : In radians (may be an array of len(Wn)).
        - zi    : Initial conditions for the filter state. It is
        a vector (or array of vectors for an N-dimensional input)
        of length 1.  If zi=None or is not given then initial
        rest is assumed.  SEE signal.lfiltic for more information.

    Returns (y, {zf}).

        - y     : The output of the phasor (digital filter).
        - zf:   : If zi is None, this is not returned, otherwise,
                  zf holds the final filter delay values.

    Return len(Wn) of a variable frequency un-wrapped phasor.
    Wn is the normalized frequency (may be an array of nchannels, for multichannel).
    Phase in radians (constant, or may be an array of nframes of nchannels).

    """
    # the phasor is actually a scaled integrator

    b = array([0., 1.]) # integrator coeffs
    a = array([1., -1.])

    if zi is None:
        return ffilter(b, a, (pi * Wn)) + phase
    
    else:
        y, zf = ffilter(b, a, (pi * Wn), zi)
        return y + phase, zf


def sinosc(Wn, phase = 0., zi = None):
    """phasor(Wn, phase = 0., zi = None)

    Args:
        - Wn    : Normalized frequency array (may be multichannel)
        - phase : In radians (may be an array of len(Wn)).
        - zi    : Initial conditions for the filter state. It is
        a vector (or array of vectors for an N-dimensional input)
        of length 1.  If zi=None or is not given then initial
        rest is assumed.  SEE signal.lfiltic for more information.

    Returns (y, {zf}).

        - y     : The output of the oscillator.
        - zf:   : If zi is None, this is not returned, otherwise,
                  zf holds the final phasor filter delay values.

    Return len(Wn) of a variable frequency sine oscillator.
    Wn is the normalized frequency (may be an array of nchannels, for multichannel).
    Phase in radians (constant, or may be an array of nframes of nchannels).


    """
    if zi is None:
        return sin(phasor(Wn, phase, zi))
    else:
        y, zf = phasor(Wn, phase, zi)
        return sin(y), zf


# **************************************
# low noise. . .
# **************************************

def low_noise(nframes, N, Wn, zi = None):
    """low_noise(nframes, N, Wn, zi = None)

    Generate low frequency noise by filtering white noise with
    a fixed frequency butterworth low pass IIR filter. Retains
    constant power with cutoff.
        
    Inputs:
    
      nframes -- nframes to generate, may be a shape tuple for multichannel.
      N       -- low pass filter order
      Wn      -- cutoff
      zi      -- Initial conditions for the filter delays.  It is a vector
                 (or array of vectors for an N-dimensional input) of length
                 max(len(a),len(b)).  If zi=None or is not given then initial
                 rest is assumed.  SEE signal.lfiltic for more information.
    
    Outputs: (y, {zf})
    
      y -- The output of the digital filter.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            final filter delay values.
    
    Algorithm:
      See white, butter and lfilter

      """
    if Wn == 1.:
        b = zeros(N+1)
        b[0] = 1.
        a = b
    else:
        b, a = butter(N, Wn, 'lowpass')

    y = ffilter(b, a, white(nframes) / sqrt(Wn), zi)
    
    return y


# **************************************
# delays. . .
# **************************************

# NOTE FOR DELAYS:
# review DAFX for other delay algorithms--particularly for feedback delays!!

# similar to SC3's BufDelayN
# f means fixed, n means no interpolation
# consider adding vdelay, delayi, etc.
def fdelayn(x, nframes, axis = 0, zi = None):
    """fdelayn(x, nframes, axis = 0, zi=None)

    Delay input x by nframes. 
    
    Description
    
      Delay input x by nframes along the given axis.
      If nframes is negative, results in abs(nframes) advance.
      
      Note: max delay = size(x, axis)
    
    Inputs:
    
      x -- A N-dimensional input array.
      nframes -- Number of frames to delay by.
      axis -- axis of choice
      zi -- Initial conditions for the delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            nframes.  If zi=None or is not given then initial rest is assumed.
    
    Outputs: (y, {zf})
    
      y -- The output of the delay.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            delay values.
    
      """
    # delay / advance
    y = roll(x, nframes, axis)

    # extract zf, which has wrapped around
    if nframes < 0:
        zf = take(y, arange(size(y, axis) + nframes, size(y, axis)), axis) # advance
    else:
        zf = take(y, arange(nframes), axis) # delay

    # generate zeros if none exist
    if zi is None:
        zi = zeros(shape(zf))
        zi_init = False
    else:
        zi_init = True

    # put appropriate zi in place--using slice tuple
    # (the put function isn't appropriate)
    s_r = (zi.ndim * [Ellipsis])

    if nframes < 0:
        s_r[axis] = slice(size(y, axis) + nframes, size(y, axis)) # advance
    else:
        s_r[axis] = slice(0, nframes) # delay

    y[s_r] = zi

    if zi_init:
        return y, zf
    else:
        return y


def fadvancen(x, nframes, axis = 0, zi = None):
    """fdelayn(x, nframes, axis = 0, zi=None)

    Advance input x by nframes. 
    
    Description
    
      Advance input x by nframes along the given axis.
      If nframes is negative, results in abs(nframes) delay.
      
      Note: max advance = size(x, axis)
    
    Inputs:
    
      x -- A N-dimensional input array.
      nframes -- Number of frames to advance by.
      axis -- axis of choice
      zi -- Initial conditions for the delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            nframes.  If zi=None or is not given then initial rest is assumed.
    
    Outputs: (y, {zf})
    
      y -- The output of the advance.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            delay values.
    
      """
    return fdelayn(x, -nframes, axis, zi)


# comb filters. . . using circular delays

# f means fixed, n means no interpolation
# plan to add damped combs!
# consider adding vcomb, combi, etc.
def faz_combn(x, nframes, b0, b1, axis = 0, zi = None):
    """faz_combn(x, nframes, b0, b1, axis = 0, zi = None)

    Filter input x by ALL ZERO comb filter of length nframes. 
    
    Description
    
      Filter input x along the given axis.
    
    Inputs:
    
      x -- A N-dimensional input array.
      nframes -- Length of comb filter.
      b0 -- The direct feed-forward coefficient.
      b1 -- The delay feed-forward coefficient.
      axis -- axis of choice
      zi -- Initial conditions for the delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            nframes.  If zi=None or is not given then initial rest is assumed.
    
    Outputs: (y, {zf})
    
      y -- The output of the delay.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            delay values.
    
      """

    if axis is not 0:
        raise 'currently only supports axis = 0'

    else:

        y = zeros_like(x)       # empty output

        if zi is None:
            if x.ndim is 1:
                z = zeros(nframes)  # nframes of a delay
            else:
                z = zeros((nframes, nchannels(x)))
        else:
            z = copy(zi)

        if isscalar(b0) and isscalar(b1): # b0, b1 as scalars

            for n in range(size(x, axis)):
                                  # n is input read index
                n_z = n % nframes # n_z is z read/write index

                y[n] = b0 * x[n] + b1 * z[n_z]
                z[n_z] = x[n]

        else:                   # upgrade for time varying b0, b1
            if isscalar(b0):
                b0 = repeat(b0, size(x, axis))
            if isscalar(b1):
                b1 = repeat(b1, size(x, axis))

            for n in range(size(x, axis)):

                n_z = n % nframes

                y[n] = b0[n] * x[n] + b1[n] * z[n_z]
                z[n_z] = x[n]

        if zi is None:
            return y
        else:
            zf = roll(z, -n_z)[::-1]
            return y, zf


def fap_combn(x, nframes, b0, a1, axis = 0, zi = None):
    """fap_combn(x, nframes, b0, a1, axis = 0, zi = None)

    Filter input x by ALL POLE comb filter of length nframes. 
    
    Description
    
      Filter input x along the given axis.
    
    Inputs:
    
      x -- A N-dimensional input array.
      nframes -- Length of comb filter.
      b0 -- The direct feed-forward coefficient.
      a1 -- The delay feed-back coefficient. (neg, to conform to transfer fn)
      axis -- axis of choice
      zi -- Initial conditions for the delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            nframes.  If zi=None or is not given then initial rest is assumed.
    
    Outputs: (y, {zf})
    
      y -- The output of the delay.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            delay values.
    
      """

    if axis is not 0:
        raise 'currently only supports axis = 0'

    else:

        y = zeros_like(x)       # empty output

        if zi is None:
            if x.ndim is 1:
                z = zeros(nframes)  # nframes of a delay
            else:
                z = zeros((nframes, nchannels(x)))
        else:
            z = copy(zi)

        if isscalar(b0) and isscalar(a1): # b0, a1 as scalars

            for n in range(size(x, axis)):
                                  # n is input read index
                n_z = n % nframes # n_z is z read/write index

                xh = x[n] - a1 * z[n_z]
                y[n] = b0 * xh
                z[n_z] = xh

        else:                   # upgrade for time varying b0, b1
            if isscalar(b0):
                b0 = repeat(b0, size(x, axis))
            if isscalar(a1):
                a1 = repeat(a1, size(x, axis))

            for n in range(size(x, axis)):

                n_z = n % nframes

                xh = x[n] - a1[n] * z[n_z]
                y[n] = b0[n] * xh
                z[n_z] = xh

        if zi is None:
            return y
        else:
            zf = roll(z, -n_z)[::-1]
            return y, zf


def fpz_combn(x, nframes, b0, b1, a1, axis = 0, zi = None):
    """fpz_combn(x, nframes, b0, b1, a1, axis = 0, zi = None)

    Filter input x by POLE ZERO comb filter of length nframes. 
    
    Description
    
      Filter input x along the given axis.
    
    Inputs:
    
      x -- A N-dimensional input array.
      nframes -- Length of comb filter.
      b0 -- The direct feed-forward coefficient.
      a1 -- The delay feed-back coefficient. (neg, to conform to transfer fn)
      axis -- axis of choice
      zi -- Initial conditions for the delays.  It is a vector
            (or array of vectors for an N-dimensional input) of length
            nframes.  If zi=None or is not given then initial rest is assumed.
    
    Outputs: (y, {zf})
    
      y -- The output of the delay.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            delay values.
    
      """

    if axis is not 0:
        raise 'currently only supports axis = 0'

    else:

        y = zeros_like(x)       # empty output

        if zi is None:
            if x.ndim is 1:
                z = zeros(nframes)  # nframes of a delay
            else:
                z = zeros((nframes, nchannels(x)))
        else:
            z = copy(zi)

        if isscalar(b0) and isscalar(b1) and isscalar(a1): # b0, a1 as scalars

            for n in range(size(x, axis)):
                                  # n is input read index
                n_z = n % nframes # n_z is z read/write index

                xh = x[n] - a1 * z[n_z]
                y[n] = b1 * z[n_z] + b0 * xh
                z[n_z] = xh

        else:                   # upgrade for time varying b0, b1
            if isscalar(b0):
                b0 = repeat(b0, size(x, axis))
            if isscalar(b1):
                b1 = repeat(b1, size(x, axis))
            if isscalar(a1):
                a1 = repeat(a1, size(x, axis))

            for n in range(size(x, axis)):

                n_z = n % nframes

                xh = x[n] - a1[n] * z[n_z]
                y[n] = b1[n] * z[n_z] + b0[n] * xh
                z[n_z] = xh

        if zi is None:
            return y
        else:
            zf = roll(z, -n_z)[::-1]
            return y, zf


# ****************************************************************************
# REVISIT EVERYTHING BELOW HERE
# ****************************************************************************


# ******************************************************************************
# 
# Classes. . .
# 
# **********************

# CLASSES NEED TO BE UPDATED TO REMOVE CHANNELS--CONFORM TO GENERATORS

# ******************************************************************************
# osc objects: filters returning signal from a freq array
#    .ar() returns complete array
#          depends on freq array assigned

# consider adding table and tablei osc
# AS OSC_OBJ IS A SPECIALIZED CLASS DON'T BRING INTO FILT_OBJ
class OscObj(MusObj):

    def __init__(self, phase_init = None, sr = None):
        MusObj.__init__(self, sr) # Run superclass init

        self.z = None    # state of the integrator
        self.phase_init = phase_init  # phase init'ed?


    def _phase_init(self, freq, phase): # initializes phase appropriately
        if phase is None:       # init phase, if set
            if self.phase_init is None:
                phase = 0.
            else:
                phase = self.phase_init
        else:
            if self.phase_init is None:
                pass
            else:
                phase += self.phase_init

        return asarray(phase)


    def _zi(self, freq, phase): # initializes state of the integrator appropriately

        if self.z is None:      # set up the state if it is None

            # determine channels of freq and phase to generate

            # set z_chans
            z_chans = max(
                nchannels(freq), # frequency channels
                ((isscalar(phase) and [None]) or [nchannels(phase)])[0] # phase channels
                )

            # create zeros. . .
            z = zeros(z_chans)
            z.shape = (z_chans, 1)

            # set z
            self.z = z


class SinOsc(OscObj):

    def __init__(self, phase_init = None, sr = None):
        OscObj.__init__(self, phase_init, sr = None) # Run superclass init

    def ar(self, freq, phase = None):
        
        # check phase is init'ed, and assign if need be
        phase = self._phase_init(freq, phase)

        # check if state is set, and assign if need be
        self._zi(freq, phase)

        # run sinosc()
        y, self.z = sinosc(
            freq_to_Wn(freq, self.T),
            phase,
            self.z)
        return y



# ******************************************************************************
# Delay objects: filters returning signal from a signal array
#    .ar() returns complete array
# consider making state a tuple, so that states for internal filters (e.g., lp)
# can be held too!
# also--consider making delay_time possible to be an array
# class DelObj(MusObj):


# ADD A FILT_OBJ CLASS HERE?. . . WHICH WILL HANDLE MULTICHANNEL EXPANSION
# AS DOES GEN_OBJ??

class FDelObj(MusObj):

    def __init__(self, delay, zi = None, sr = None, asnframes = False):
        MusObj.__init__(self, sr) # Run superclass init

        if asnframes:           # convert and set nframes
            nframes = delay
        else:
            nframes = dur_to_nframes(delay, self.sr)

        if isscalar(nframes):
            self.nframes = nframes
        else:
            self.nframes = asarray(nframes)

        self.nchannels = nchannels(asarray([delay])) # store nchannels

        self.z = zi        # state of the delay, may be a tuple, as needed


    def ar(self, x): # array return

        # check if state is set, and assign if need be
        self._zi(x)

        if self.nchannels is 1: # run fdelayn()
            y, self.z = self.fun(
                x,
                self.nframes,
                0,
                self.z)

        else:                         # multichannel, run fdelayn()
            y = range(self.nchannels) # set up two empty lists to fill
            z = range(self.nchannels) # with results of delay

            for n, x_n in zip(y, deinterleave(x)): # count through nchannels
                y[n], z[n] = self.fun(
                    x_n,
                    self.nframes[n],
                    0,
                    self.z[n])
            y = interleave(asarray(y))
            self.z = z

        return y


    def _zi(self, x):   # initializes state of the delay appropriately

        if self.z is None:      # set up the state if it is None

            if self.nchannels is 1: # create zeros - as a block. . .
                zi = zeros((abs(self.nframes), nchannels(x))) # generate zi zeros
                self.z = zi

            else:               # create zeros - as an array of nchannel arrays. . .
                if self.nchannels != nchannels(x):
                    raise ValueError, 'delay argument nchannels != nchannels(x)'
                else:
                    zi = range(self.nchannels)
                    for n in range(self.nchannels):
                        zi[n] = zeros(abs(self.nframes[n]))
                    self.z = zi


class FDelay(FDelObj):

    def __init__(self, delay = 1., zi = None, sr = None, asnframes = False):
        self.fun = fdelayn      # define function
        FDelObj.__init__(self, delay, zi, sr, asnframes) # Run superclass init


# combs. . .
# might be good to work out how to combine DEL_OBJ and CMB_OBJ into one class
class FCmbObj(MusObj):

    def __init__(self, freq, zi = None, sr = None):
        MusObj.__init__(self, sr) # Run superclass init

        # set flags for cos or sin response
        self.sign = sign(freq)  # negative freq --> sin response

        # set up nframes
        nframes = dur_to_nframes(
            freq_to_period(abs(freq)),
            self.sr) # convert and set nframes
        
        self.nframes = nframes

        self.nchannels = nchannels(asarray([nframes])) # store nchannels

        self.z = zi        # state of the delay, may be a tuple, as needed


    def _ar(self, x, *args): # array return

        # check if state is set, and assign if need be
        self._zi(x)

        if self.nchannels is 1: # run comb()
            if self.fun is fpz_combn:
                y, self.z = self.fun(
                    x,
                    self.nframes,
                    args[0],
                    args[1],
                    args[2],
                    0,
                    self.z)
            else:
                y, self.z = self.fun(
                    x,
                    self.nframes,
                    args[0],
                    args[1],
                    0,
                    self.z)
        else:                         # multichannel, run comb()
            y = range(self.nchannels) # set up two empty lists to fill
            z = range(self.nchannels) # with results of delay

            for n, x_n in zip(y, deinterleave(x)): # count through nchannels
                if self.fun is fpz_combn:
                    y[n], z[n] = self.fun(
                        x_n,
                        self.nframes[n],
                        args[0][n],
                        args[1][n],
                        args[2][n],
                        0,
                        self.z[n])
                else:
                    y[n], z[n] = self.fun(
                        x_n,
                        self.nframes[n],
                        args[0][n],
                        args[1][n],
                        0,
                        self.z[n])
            y = interleave(asarray(y))
            self.z = z

        return y


    def _zi(self, x):   # initializes state of the delay appropriately

        if self.z is None:      # set up the state if it is None

            if self.nchannels is 1: # create zeros - as a block. . .
                zi = zeros((abs(self.nframes), nchannels(x))) # generate zi zeros
                self.z = zi

            else:               # create zeros - as an array of nchannel arrays. . .
                if self.nchannels != nchannels(x):
                    raise ValueError, 'delay argument nchannels != nchannels(x)'
                else:
                    zi = range(self.nchannels)
                    for n in range(self.nchannels):
                        zi[n] = zeros(abs(self.nframes[n]))
                    self.z = zi


class FAzComb(FCmbObj):

    def __init__(self, freq = 440., zi = None, sr = None):
        self.fun = faz_combn      # define function
        FCmbObj.__init__(self, freq, zi, sr) # Run superclass init

    def ar(self, x, gain):      # gain in db


        # CONSIDER MOVING THIS TO _AR, ABOVE
        # massage gain to broadcast correctly
        # may want to rework this to take advantage of
        # numpy broadcasting
        if self.nchannels is 1:
            if not isscalar(gain):  # is gain an array?
                gain = asarray(gain)

                if nframes(gain) is nchannels(x): # upgrade gain
                    gain = tile(gain, (nframes(x), self.nchannels))

        else:               # multichannel broadcasting--e.g., need to
                            # use more than one comb function call
            if isscalar(gain):
                gain = repeat(gain, self.nchannels)

            elif nchannels(asarray(gain)) is 1 and nframes(asarray(gain)) is nframes(x):
                gain = tile(gain, self.nchannels)


        k = db_to_amp(gain)     # calculate coeffs here!
        b0 = .5 * (1 + k)
        b1 = .5 * (1 - k) * self.sign

        # need to reshape for multichannel broadcasing
        if shape(k) == shape(x) and self.nchannels > 1:
            b0 = deinterleave(b0)
            b1 = deinterleave(b1)

        return self._ar(x, b0, b1)


class FApComb(FCmbObj):

    def __init__(self, freq = 440., zi = None, sr = None, asT60 = False):
        self.fun = fap_combn      # define function
        self.asT60 = asT60
        FCmbObj.__init__(self, freq, zi, sr) # Run superclass init

    def ar(self, x, gain):      # gain in db


        # CONSIDER MOVING THIS TO _AR, ABOVE
        # massage gain to broadcast correctly
        # may want to rework this to take advantage of
        # numpy broadcasting
        if self.nchannels is 1:
            if not isscalar(gain):  # is gain an array?
                gain = asarray(gain)

                if nframes(gain) is nchannels(x): # upgrade gain
                    gain = tile(gain, (nframes(x), self.nchannels))

        else:               # multichannel broadcasting--e.g., need to
                            # use more than one comb function call
            if isscalar(gain):
                gain = repeat(gain, self.nchannels)

            elif nchannels(asarray(gain)) is 1 and nframes(asarray(gain)) is nframes(x):
                gain = tile(gain, self.nchannels)

        # calculate coeffs here!
        if not self.asT60:      # as gain
            k = db_to_amp(gain)

            a1 = (k - 1.)/(k + 1.) * self.sign
            b0 = (1 - abs(a1))

            # need to reshape for multichannel broadcasing
            if shape(k) == shape(x) and self.nchannels > 1:
                b0 = deinterleave(b0)
                b1 = deinterleave(b1)

        else:                   # as T60
            
            a1 = -pow(10, -3 * nframes_to_dur(self.nframes, self.sr) / gain) * self.sign
            b0 = (1 - abs(a1))

        return self._ar(x, b0, a1)

