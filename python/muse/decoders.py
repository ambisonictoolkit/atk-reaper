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
from filters import *
from transforms import *

# import muse defined constants
import constants as C

#=========================
# Functions
#=========================

#------------------------------------------------------------------------
# (Gerzon's) Diametric Decoder Theorem (DDT)
#------------------------------------------------------------------------
#
# Much of the code below is a transcoding of Aaron Heller's Octave
# code available at: http://www.ai.sri.com/ajh/ambisonics/
#
# Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"
# Preprint from AES-121, 10/2006, San Francisco
#
# Heller's original functions are noted through comments in each
# functions help field.
#
# Transcoding to Python/Numpy for use in muse/ATK by
# Joseph Anderson <josephlloydanderson@mac.com>
#
# aes_paper.m (expanded version of speaker_matrix.m) contains the
# following functions:
#
#   velocity_gain_matrix()**            : compute alpha, beta, and gamma
#   speaker_matrix()                    : compute alpha, beta, and gamma
#   decoder_gain_matrix()               : compute decoder matrix
#
#   rV()                                : compute the Makita direction and rV
#   rE()                                : compute rE (and direction?)
#
#   _virtual_mic()                      : virtual mic angle and directivity
#   decoder_matrix_to_virtual_mic()     : computes loudspeaker 'virtual mics' 
#
# ----------------------------------------
# the following functions are not included
# as they duplicate muse/ATK functionality
#
#   az2dir()                            : convert azimuth to directon cosines
#   degrees()                           : convert radians to degrees
#   radians()                           : convert degrees to radians
#   gain_to_db()
#
#   rectangular_speaker_arrays()        : example decodes
#   hexagonal_speaker_arrays()          : example decodes
#
#
# NOTE: speaker_matrix() and velocity_gain_matrix() are the same code.
#       It appears these two separate names are used (in error) in Heller's
#       code. The expanded version, aes_paper.m defines velocity_gain_matrix(),
#       but calls speaker_matrix().
#
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# DDT and related decoder matrix gains
#
#   NOTE:   These are the functions that compute gains to generate
#           loudspeaker feeds, and are not the functions which return
#           decoded B-format. See decoders, below.
#
#
#   speaker_matrix                  Heller's DDT (helper function)
#   decoder_gain_matrix             Heller's DDT (returns decoder gains)
#   panto_reg_decoder_gain_matrix   pantophonic
#   peri_reg_decoder_gain_matrix    periphonic
#   quad_decoder_gain_matrix        quad
#
#------------------------------------------------------------------------

def speaker_matrix(positions, k):
    """speaker_matrix(positions, k)
    
    Args:
        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
        
                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites


    Compute alpha, beta, and gamma. Return values are the weights for X, Y,
    and Z as columns of the matrix.

    retval:  alpha_0 ... alpha_i  ... alpha_n-1
             beta_0  ... beta_i   ... beta_n-1
             gamma_0 ... gamma_i  ... gamma_n-1 

    where the signal for the i'th speaker pair is:

    S_i = W +/- ( alpha_i*X + 
                  beta_i*Y +
                  gamma_i*Z )

    Note:   This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    

    Transcoding of Aaron Heller's Octave code available at:
    http://www.ai.sri.com/ajh/ambisonics/

    Transcoding to Python/Numpy for use in muse/ATK by
    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # allow entry of positions as
    # transpose for convenience
    # e.g., speaker positions are now in columns
    # rather than rows
    positions = transpose(positions)

    # n = number of speaker pairs
    # m = number of dimensions,
    #        2=horizontal, 3=periphonic 
    m, n = shape(positions)     

    # scatter matrix accumulator
    s = zeros((m, m))

    # speaker directions matrix
    directions = zeros((m, n))

    for i in range(n):

        # get the i'th speaker position
        # e.g., select the i'th column
        pos = positions[:,i]

        # normalize to get direction cosines
        dir = pos / sqrt(sum(pos**2))

        # form scatter matrix and accumulate
        s += outer(dir, dir)

        # form matrix of speaker directions
        directions[:,i] = dir

    res = sqrt(1./2) * n * k * dot(matrix(s).I, directions)

    return res


def decoder_gain_matrix(positions, k):
    """decoder_gain_matrix(positions, k)
    
    Args:
        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
        
                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
        

    Compute decoder gain matrix. Returns matrix G such that:
        [S_0, S_1, ... S_n-1]' = G [W, X, Y, Z]'


    Note:   This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    

    Transcoding of Aaron Heller's Octave code available at:
    http://www.ai.sri.com/ajh/ambisonics/

    Transcoding to Python/Numpy for use in muse/ATK by
    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # list all of the speakers
    # i.e., expand to actual pairs
    positions2 = vstack((positions, -positions))

    # get velocity gains
    # NOTE: this comment from Heller seems to be slightly
    #       misleading, in that the gains returned will be
    #       scaled by k, which may not request a velocity
    #       gain. I.e., k = 1 isn't necessarily true, as it
    #       is assigned as an argument to this function.
    sm = speaker_matrix(positions2, k)

    # n = number of speakers (NOTE: Heller says 'speaker pairs')
    # m = number of dimensions,
    #        2=horizontal, 3=periphonic 
    m, n = shape(sm)     

    # build decoder matrix 
    # rows are W, X, and Y gains
    # NOTE: this matrix construction can be simplified
    #       with a concatenation (hstack) of a column
    #       of ones and sm
    # ALSO: the below code calls for the complex conjugate
    #       of decoder_matrix. As we are expecting real vaules,
    #       we may regard this call as redundant.
    decoder_matrix = ones((m + 1, n))
    for i in range(n):
        for j in range(m):
            decoder_matrix[j + 1, i] = sm[j, i]

    res = sqrt(2)/n * decoder_matrix.conj().transpose()

    return res


def panto_reg_decoder_gain_matrix(num_speakers, orientation, k):
    """panto_reg_decoder_gain_matrix(num_speakers, orientation, k)
    
    Args:
        - num_speakers  : number of loudspeakers

        - orientation   : "flat" if the front bisects a side of
                            the polygon. The first speaker will
                            be the one left of center. "point" if
                            the front is a vertex of the polygon.
                            The first speaker will be directly in front.

                            Loudspeakers are returned in counter-clockwise
                            orientation.
                            
        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites


    Compute decoder gain matrix for a pantophonic (horizontal) regular
    polygon loudspeaker array. Returned gains are in counter-clockwise
    order from front centre.

    Returns matrix G such that:
        [S_0, S_1, ... S_n-1]' = G [W, X, Y]'


    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """

    # define function to return theta
    def theta(speaker):
        if orientation == 'point':
            theta = ((2. * speaker)/num_speakers) * pi # for 'point' case
        else :
            theta = ((1. + (2. * speaker))/num_speakers) * pi # for 'flat' case, default

        return theta


    # define constants
    g0 = 1.
    g1 = C.sqrt2

    # calculate decoding matrix
    decoder_matrix = []                # start with empty list
    speakers = range(num_speakers)

    for speaker in speakers:
        decoder_matrix.append([
                g0,
                k * g1 * cos(theta(speaker)),
                k * g1 * sin(theta(speaker))
                ])

    res = C.sqrt2/num_speakers * asarray(decoder_matrix)

    return res


def peri_reg_decoder_gain_matrix(num_speaker_pairs, elevation, orientation, k):
    """peri_reg_decoder_gain_matrix(num_speaker_pairs, orientation, k)
    
    Args:
        - num_speaker_pairs  : number of loudspeaker pairs

        - elevation     : elevation (radians) of the upper polygon

        - orientation   : "flat" if the front bisects a side of
                            the upper polygon. The first speaker will
                            be the one left of center in the upper
                            polygon. "point" if the front is a vertex of
                            the upper polygon. The first speaker will be
                            upper-front.

        - k         : W gain factor, corresponding to directivity

                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
                            Loudspeakers are returned in counter-clockwise
                            orientation.
                            

    Compute decoder gain matrix for a periphonic (full 3D) regular
    polygon loudspeaker array. Returned gains are in counter-clockwise
    order from front centre. The upper array is returned first,
    followed by the lower array. 

    Returns matrix G such that:
        [S_0, S_1, ... S_n-1]' = G [W, X, Y, Z]'


    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """

    # generate speaker pair positions
    # start with polar positions. . .
    theta = array([])
    for spkr_pr in range(num_speaker_pairs):
        theta = append(theta, 2 * pi * spkr_pr / num_speaker_pairs)

    if orientation is 'flat':
        theta += pi / num_speaker_pairs

    # convert to [ r, theta, phi ]
    spher = interleave(
        array([
            ones(num_speaker_pairs),
            theta,
            repeat(elevation, num_speaker_pairs)
            ])
        )

    # . . . then convert from spherical to cartesian
    positions = spher_to_cart(spher)

    # compute the decoder
    decoder_matrix = decoder_gain_matrix(positions, k)

    # reorder the bottom polygon
    top, bottom = split(decoder_matrix, 2)

    if orientation == 'flat' and mod(num_speaker_pairs, 2) == 1:
        bottom = roll(bottom, num_speaker_pairs/2 + 1, 0) #odd, flat
    else:
        bottom = roll(bottom, num_speaker_pairs/2, 0)   # odd, point
                                                        # even, flat
                                                        # even, point
    decoder_matrix = vstack((top, bottom))

    res = decoder_matrix

    return res


def quad_decoder_gain_matrix(angle, k):
    """quad_decoder_gain_matrix(angle, k)
    
    Args:
        - angle              : Front pair 1/2 angle
                            
        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites


    Compute decoder gain matrix for a quad (horizontal)
    loudspeaker array. The layout is adjustable via angle arg.
    Returned gains are in counter-clockwise order from front centre.

    Returns matrix G such that:
        [S_0, S_1, ... S_3]' = G [W, X, Y]'


    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """

    # calculate alpha, beta (scaled by k)
    alpha   = k / (C.sqrt2 * cos(angle))
    beta    = k / (C.sqrt2 * sin(angle))

    # fill decoding matrix
    decoder_matrix = array([
        [1.,  alpha,  beta],
        [1., -alpha,  beta],
        [1., -alpha, -beta],
        [1.,  alpha, -beta]
        ])

    res = C.sqrt2/4 * decoder_matrix

    return res


#---------------------------------------------
# UHJ and binaural decoder kernels
#---------------------------------------------

def uhj_decoder_kernel(N):
    """uhj_decoder_kernel(N)
    
    Generate a filter kernel suitable for b-format to UHJ decoding.

    See:

    Audio Engineering Society E-Library
    Ambisonic Decoders for HDTV
    Preprint Number:   3345    Convention:   92 (February 1992)
    Authors:   Gerzon, Michael A.; Barton, Geoffrey J.
    E-library Location: (CD aes12)   /pp9193/pp9203/3405.pdf

    Args:
        - N         : order of filter (number of taps)

    Outputs:
    
        - b         : coefficients of length N FIR filter: 
                        [[W_left_FIR, W_right_FIR],
                         [X_left_FIR, X_right_FIR],
                         [Y_left_FIR, Y_right_FIR]]

    Joseph Anderson <josephlloydanderson@mac.com>

    """

    #     S = 0.9396926*W + 0.1855740*X
    #     D = j(-0.3420201*W + 0.5098604*X) + 0.6554516*Y
    #     Left = (S + D)/2.0
    #     Right = (S - D)/2.0

    #---------------------------------
    # UHJ coefficients
    c_0 = 0.9396926
    c_1 = 0.1855740
    c_2 = -0.3420201
    c_3 = 0.5098604
    c_4 = 0.6554516

    gains = 0.5 * array([
        [complex(c_0, c_2), complex(c_0, -c_2)],
        [complex(c_1, c_3), complex(c_1, -c_3)],
        [complex(c_4, 0), complex(-c_4, 0)]
        ])

    m = 3               # harmonics (W, X, Y)

    #---------------------------------
    # calculate kernels

    hilbert = fir_hb(N)

    decoder_kernels = zeros((m, N, 2))  # harmonics, N, stereo

    # collect decoder kernel
    for i in range(m):          # i is harmonic number
        decoder_kernels[i] += gains[i].real * interleave(hilbert.real)
        decoder_kernels[i] += gains[i].imag * interleave(hilbert.imag)

    return decoder_kernels


def sHRIR_decoder_kernel(positions, k, N, T, \
                             r = 0.0875, theta_e = 5./9*pi, width = pi, \
                            phase = True):
    """sHRIR_decoder_kernel(positions, k, N, T, \
                           r = 0.0875, theta_e = 5./9*pi, width = pi)
    
    DDT, spherical HRIR model FIR Filter decoder kernel using spherical
    HRIR demonstrated by Brown and Duda, windowed with the Kaiser window.

    Decoder is Aaron Heller's DDT implementation. 

    See:

    Brown, C. Phillip and Richard O. Duda. "An Efficient HRTF Model for
    3-D Sound." In proceedings: IEEE Workshop on Applications of Signal
    Processing to Audio and Acoustics (WASPAA), New Paltz, NY 1997.

    http://www.ai.sri.com/ajh/ambisonics/

    Args:
        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
        
                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
        
        - N         : order of filter (number of taps)
        - T         : sampling period, 1./sr 
        - r         : sphere radius (default is Brown/Duda value)
        - theta_e   : +/- ear angle (default is Brown/Duda value)
        - width     : beta for Kaiser window FIR design.
                      pi = minimum ripple for steepest cutoff.
        - phase     : retain phase?

    Outputs:
    
        - b         : coefficients of length N FIR filter: 
                        [[W_left_FIR, W_right_FIR],
                        [[X_left_FIR, X_right_FIR],
                        [[Y_left_FIR, Y_right_FIR],
                        [[Z_left_FIR, Z_right_FIR]]

    Notes:  This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    
            As as simple sphereical model, the generated HRIR does not
            include torso or pinnae effects. Because of this, the
            simple spherical model does not generate elevation cues.


    Joseph Anderson <josephlloydanderson@mac.com>

    """
    gains = decoder_gain_matrix(positions, k)
    n, m = shape(positions)         # n = number of speaker pairs
                                    # m = number of dimensions,
                                    #        2=panto, 3=peri 
    positions2 = vstack((positions, -positions))

    speaker_kernels = empty((2 * n, N, 2))  # speakers, N, sHRIR channels
    decoder_kernels = zeros((m + 1, N, 2))  # harmonics, N, sHRIR channels

    # collect decoder kernel
    for i in range(2 * n):          # i is speaker number

        # collect speaker sHRIR kernels
        if m is 2:              # (2D)
            radius, azimuth = cart_to_pol(positions2)[i]
            elevation = 0.
        else:                   # (3D)
            radius, azimuth, elevation = cart_to_spher(positions2)[i]
                
        # find speaker kernel
        speaker_kernels[i] = sHRIR(N, azimuth, elevation, T, \
                                   r, theta_e, width)

        # discard phase for speaker..?
        if not phase:
            for chan in range(nchannels(speaker_kernels[i])):
                speaker_kernels[i][:, chan] = linf(speaker_kernels[i][:, chan])

        # sum to decoder kernels
        for j in range(m + 1):          # j is harmonic number
            decoder_kernels[j] += gains[i, j] * speaker_kernels[i]

    return decoder_kernels


def lHRIR_decoder_kernel(positions, k, subject_id, database_dir, status = 'C', \
                            phase = True):
    """lHRIR_decoder_kernel(positions, k, subject_id, database_dir, \
                                status = 'C', phase = True)
    
    DDT / HRIR FIR Filter decoder using measured HRIRs from the
    IRCAM hosted Listen HRTF database.

    Decoder is Aaron Heller's DDT implementation. 

    See: http://recherche.ircam.fr/equipes/salles/listen/
         http://www.ai.sri.com/ajh/ambisonics/

    Args:
        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
        
                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
        
        - subject_id    : 1002 - 1059 (as string)
        - database_dir  : path to local directory where subject directories
                          for the Listen database are located                          
        - status        : compensated ('C') or raw ('R') HRIRs
                            len('C' ) = 512, len('R') = 8192
        - phase         : retain phase?

    Outputs:
    
        - b         : coefficients of length N FIR filter: 
                        [[W_left_FIR, W_right_FIR],
                        [[X_left_FIR, X_right_FIR],
                        [[Y_left_FIR, Y_right_FIR],
                        [[Z_left_FIR, Z_right_FIR]]

    Notes:  This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    
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

    Source is at a distance of 1.95 meter. Use this value for NFC filtering.


    Note: Returns the complete (asymmetric) HRIR, measured at SR = 44.1kHz


    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # set kernel length (fixed for Listen HRIR database)
    if status is 'R':
        N = 8192
    else:
        N = 512
    
    gains = decoder_gain_matrix(positions, k)
    n, m = shape(positions)         # n = number of speaker pairs
                                    # m = number of dimensions,
                                    #        2=panto, 3=peri 
    positions2 = vstack((positions, -positions))

    speaker_kernels = empty((2 * n, N, 2))  # speakers, N, lHRIR channels
    decoder_kernels = zeros((m + 1, N, 2))  # harmonics, N, lHRIR channels

    # collect decoder kernel
    for i in range(2 * n):          # i is speaker number

        # collect speaker lHRIR kernels
        if m is 2:              # (2D)
            radius, azimuth = cart_to_pol(positions2)[i]
            elevation = 0.
        else:                   # (3D)
            radius, azimuth, elevation = cart_to_spher(positions2)[i]
                
        # find speaker kernel
        speaker_kernels[i] = lHRIR(azimuth, elevation, \
                                   subject_id, database_dir, status)

        # discard phase for speaker..?
        if not phase:
            for chan in range(nchannels(speaker_kernels[i])):
                speaker_kernels[i][:, chan] = linf(speaker_kernels[i][:, chan])

        # sum to decoder kernels
        for j in range(m + 1):          # j is harmonic number
            decoder_kernels[j] += gains[i, j] * speaker_kernels[i]

    return decoder_kernels


def cHRIR_decoder_kernel(positions, k, subject_id, database_dir, phase = True):
    """cHRIR_decoder_kernel(positions, k, subject_id, database_dir, \
                                phase = True)
    
    DDT / HRIR FIR Filter decoder using measured HRIRs from the
    CIPIC HRTF database.

    Decoder is Aaron Heller's DDT implementation. 

    See: http://interface.cipic.ucdavis.edu/sound/hrtf.html
         http://www.ai.sri.com/ajh/ambisonics/

    Args:
        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
        
                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
        
        - subject_id    : 003 - 165 (as string)
                          subjects '021' and '165' are the KEMAR head
        - database_dir  : path to local directory where subject directories
                          for the Listen database are located                          
        - phase         : retain phase?

    Outputs:
    
        - b         : coefficients of length N FIR filter: 
                        [[W_left_FIR, W_right_FIR],
                        [[X_left_FIR, X_right_FIR],
                        [[Y_left_FIR, Y_right_FIR],
                        [[Z_left_FIR, Z_right_FIR]]

    Notes:  This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    
    Please see "Documentation for the UCD HRIR Files" avaliable at the
    above link for measurement details. As a whole, there are 1250
    measurement points.

    Source is at a distance of 1.0 meter. Use this value for NFC filtering.

    Note: Returns the complete (asymmetric) HRIR, measured at SR = 44.1kHz


    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # set kernel length (fixed for CIPIC HRIR database)
    N = 256
    
    gains = decoder_gain_matrix(positions, k)
    n, m = shape(positions)         # n = number of speaker pairs
                                    # m = number of dimensions,
                                    #        2=panto, 3=peri 
    positions2 = vstack((positions, -positions))

    speaker_kernels = empty((2 * n, N, 2))  # speakers, N, cHRIR channels
    decoder_kernels = zeros((m + 1, N, 2))  # harmonics, N, cHRIR channels

    # collect decoder kernel
    for i in range(2 * n):          # i is speaker number

        # collect speaker cHRIR kernels
        if m is 2:              # (2D)
            radius, azimuth = cart_to_pol(positions2)[i]
            elevation = 0.
        else:                   # (3D)
            radius, azimuth, elevation = cart_to_spher(positions2)[i]
                
        # find speaker kernel
        speaker_kernels[i] = cHRIR(azimuth, elevation, \
                                   subject_id, database_dir)

        # discard phase for speaker..?
        if not phase:
            for chan in range(nchannels(speaker_kernels[i])):
                speaker_kernels[i][:, chan] = linf(speaker_kernels[i][:, chan])

        # sum to decoder kernels
        for j in range(m + 1):          # j is harmonic number
            decoder_kernels[j] += gains[i, j] * speaker_kernels[i]

    return decoder_kernels


#---------------------------------------------
# rV and rE analysis (in XY plane)
#---------------------------------------------

def rV2(decoder_matrix, positions2, source_angle):
    """rV2(decoder_matrix, positions2, source_angle)
    
    Args:
        - decoder_matrix : input decoder matrix, e.g.,
                            values returned by decoder_gain_matrix()

        - positions2    : complete set of loudspeaker positions
                            for DDT, this is the mirrored set

        - source_angle  : incidence angle (radians) of incoming source
        

    Compute Makita direction and rV (in the horizontal plane)


    Note:   This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    

    Transcoding of Aaron Heller's Octave code available at:
    http://www.ai.sri.com/ajh/ambisonics/

    Transcoding to Python/Numpy for use in muse/ATK by
    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # a signal from direction angle (in radians) in B-format
    # NOTE: only in the horizontal plane!
    #       Additionally, we could use Muse's encoder here
    source_WXY = array([
            sqrt(1./2),                      # W
	    cos(source_angle),              # X
	    sin(source_angle)               # Y
            ])

    # the gains from the signal to each speaker
    source_to_speaker_gains = dot(decoder_matrix, source_WXY)

    # the Velocity gains gains multiplied by 
    #  the directions to the speakers
    #  normalized by the Pressure gain
    source_to_speaker_gains_XY = \
        dot(source_to_speaker_gains.conj().T, positions2) / \
	sum(source_to_speaker_gains, 0)

    # convert to length and direction
    #  direction is the Makita direction
    #  mag is rV
    rVpolar = cart_to_pol(source_to_speaker_gains_XY)

    res = rVpolar

    return res


def rE2(decoder_matrix, positions2, source_angle):
    """rE2(decoder_matrix, positions2, source_angle)
    
    Args:
        - decoder_matrix : input decoder matrix, e.g.,
                            values returned by decoder_gain_matrix()

        - positions2    : complete set of loudspeaker positions
                            for DDT, this is the mirrored set

        - source_angle  : incidence angle (radians) of incoming source
        

    Compute energy direction and rE (in the horizontal plane)


    Note:   This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    

    Transcoding of Aaron Heller's Octave code available at:
    http://www.ai.sri.com/ajh/ambisonics/

    Transcoding to Python/Numpy for use in muse/ATK by
    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # a signal from direction angle (in radians) in B-format
    # NOTE: only in the horizontal plane!
    source_WXY = array([
            sqrt(1./2),                      # W
	    cos(source_angle),              # X
	    sin(source_angle)               # Y
            ])

    # the gains from the signal to each speaker
    source_to_speaker_gains = dot(decoder_matrix, source_WXY)

    # power gain
    source_to_speaker_gains2 = \
        abs( source_to_speaker_gains )**2

    # the Velocity gains gains squared (=energy) 
    #  multiplied by the directions to the speakers
    #  normalized by the Pressure gain
    source_to_speaker_gains_XY2 = \
        dot(source_to_speaker_gains2.conj().T, positions2) / \
	sum(source_to_speaker_gains2, 0)

    # convert to length and direction
    #  direction is the Makita direction
    #  mag is rV
    rEpolar = cart_to_pol(source_to_speaker_gains_XY2)

    res = rEpolar

    return res


#------------------------------------------------------------------------
# 'virtual mic' decoder analysis (in XY plane)
#------------------------------------------------------------------------

def _virtual_mic2(w, x, y):
    """_virtual_mic2( w, x, y )
    
    Args:
        - w, x, y : corresponding w, x, y coefficients
                    from a decoder matrix

    Helper function for decoder_matrix_to_virtual_mic


    Note:   This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    

    Transcoding of Aaron Heller's Octave code available at:
    http://www.ai.sri.com/ajh/ambisonics/

    Transcoding to Python/Numpy for use in muse/ATK by
    Joseph Anderson <josephlloydanderson@mac.com>

    """

    polar = cart_to_pol(array([x, y]))
    v = sqrt(2) * polar[0]
    res = array([ polar[1], v / (v + w) ])

    return res


def decoder_matrix_to_virtual_mic2(decoder_matrix):
    """decoder_matrix_to_virtual_mic2(decoder_matrix)
    
    Args:
        - decoder_matrix : decoder matrix returned by
                            decoder_matrix()


    Note:   This assumes standard B format
            definitions for W, X, Y, and Z, i.e., W
            is sqrt(2) lower than X, Y, and Z.
    

    Transcoding of Aaron Heller's Octave code available at:
    http://www.ai.sri.com/ajh/ambisonics/

    Transcoding to Python/Numpy for use in muse/ATK by
    Joseph Anderson <josephlloydanderson@mac.com>

    """

    # m speakers, n coefficients
    m, n = shape(decoder_matrix)     

    # virtual mic matrix
    vm_matrix = zeros((m, 2))

    # for loop
    # NOTE: there is likely a more Pythonic/Numpinic
    #       way of doing things here!
    for i in range(m):
        vm_matrix[i,:] = _virtual_mic2(
            decoder_matrix[i, 0],
            decoder_matrix[i, 1],
            decoder_matrix[i, 2]
            )

    res = vm_matrix

    return res


#------------------------------------------------------------------------
# Muse/ATK decoders built using the above decoder gain matrices
#   Developed and coded by J Anderson
#
#
#   ++++ Single Band Regular Decoders
#
#   panto_sbr
#   peri_sbr
#
#   ++++ Single Band Diametric Decoders
#
#   decode_sbd
#   quad_sbd
#
#
#   ++++ Dual Band Regular Decoders (Shelf-filtered)
#
#   panto_dbr
#   peri_dbr
#
#   ++++ Dual Band Diametric Decoders (Shelf-filtered)
#
#   decode_dbd
#   quad_dbd        (See special decoders)
#
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# Single Band Regular Decoders
#
#   panto_sbr           pantophonic
#   peri_sbr            periphonic
#------------------------------------------------------------------------

def panto_sbr(a, num_speakers = 4, orientation = 'flat', k = 0.7071):
    """panto_sbr(a, num_speakers = 4, orientation = 'flat', k = 0.7071)
    
    Pantophonic Single Band Regular Array Decoder
    
    Args:
        - num_speakers  : number of loudspeakers

        - orientation   : "flat" if the front bisects a side of
                            the polygon. The first speaker will
                            be the one left of center. "point" if
                            the front is a vertex of the polygon.
                            The first speaker will be directly in front.

                            Loudspeakers are returned in counter-clockwise
                            orientation.
                            
        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites

                         "velocity" returns "strict soundfield" decoding
                         "energy" is preferred for multiple listeners
                         "controlled opposites" is preferred for
                             acoustically live spaces (concert halls)


    Decode a three dimensional ambisonic B-format signal to a set of
    loudspeakers in a regular, horizontal polygon (a ring array). The
    outputs are in counter-clockwise order, and the position of the first
    speaker is either left of centre or centre.

    This decoder does not employ psychoacoustic shelf filtering.

    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """
    # return decoder gain matrix
    decoder = panto_reg_decoder_gain_matrix(num_speakers, orientation, k)

    # decode here!
    return inner(a[:, 0:-1], array(decoder))


def peri_sbr(a, num_speaker_pairs = 4, elevation = 0.6155, orientation = 'flat', \
             k = 0.5774):
    """peri_sbr(a, num_speaker_pairs = 4, elevation = 0.6155,
                                    orientation = 'flat', k = 0.5774)
    
    Periphonic Single Band Regular Array Decoder
    
    Args:
        - num_speaker_pairs  : number of loudspeaker pairss

        - orientation   : "flat" if the front bisects a side of
                            the upper polygon. The first speaker will
                            be the one left of center in the upper
                            polygon. "point" if the front is a vertex of
                            the upper polygon. The first speaker will be
                            upper-front.

                            Loudspeakers are returned in counter-clockwise
                            orientation.
                            
        - k         : W gain factor, corresponding to directivity

                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
                         
                         "velocity" returns "strict soundfield" decoding
                         "energy" is preferred for multiple listeners
                         "controlled opposites" is preferred for
                             acoustically live spaces (concert halls)


    Decode a three dimensional ambisonic B-format signal to a set of
    loudspeakers in a two regular, horizontal polygons: an upper and
    lower ring array. The outputs are in counter-clockwise order, and
    the position of the first speaker is either left of centre or centre.
    The upper array is returned first, followed by the lower array.

    This decoder does not employ psychoacoustic shelf filtering.

    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """
    # return decoder gain matrix
    decoder = peri_reg_decoder_gain_matrix(num_speaker_pairs, \
                                           elevation, orientation, k)

    # decode here!
    return inner(a, array(decoder))


#------------------------------------------------------------------------
# Single Band Diametric Decoder
#
#   decode_sbd
#   quad_sbd
#------------------------------------------------------------------------

def decode_sbd(a, positions, k = 0.7071):
    """decode_sbd(a, positions, k = 0.7071)
    
    Diametric Pairs Single Band Decoder (DDT)
    
    Args:
        - a                  : Input B-format signal

        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
        
                      For periphonic (3D)
                      k: 1         => velocity,
                         sqrt(1/3) => energy, 
                         1/3       => controlled opposites
                         
                         "velocity" returns "strict soundfield" decoding
                         "energy" is preferred for multiple listeners
                         "controlled opposites" is preferred for
                             acoustically live spaces (concert halls)


    Decode a three dimensional ambisonic B-format signal to a set of
    loudspeakers in an array as specified by positions arg.

    This decoder does not employ psychoacoustic shelf filtering.

    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """
    # return decoder gain matrix
    decoder = decoder_gain_matrix(positions, k)

    # pantophonic or periphonic?
    if shape(positions)[1] is 2:
        a = a[:, 0:-1]

    # decode here!
    return inner(a, array(decoder))


def quad_sbd(a, angle = 0.7854, k = 0.7071):
    """quad_sbd(a, angle = 0.7854, k = 0.7071)
    
    Quadraphonic Single Band Decoder
    
    Args:
        - a                  : Input B-format signal

        - angle              : Front pair 1/2 angle

        - k         : W gain factor, corresponding to directivity

                      For pantophonic (2D)
                      k: 1         => velocity,
                         sqrt(1/2) => energy, 
                         1/2       => controlled opposites
                         
    Decode a three dimensional ambisonic B-format signal to a
    quadraphonic (horizontal only) loudspeaker array. The layout
    is adjustable via the angle arg. Returned gains are in counter-
    clockwise order from front centre.

    The outputs will be:
    
         [ Front Left, Back Left, Back Right, Front Right ]


    This decoder does not employ psychoacoustic shelf filtering.

    """
    # return decoder gain matrix
    decoder = quad_decoder_gain_matrix(angle, k)

    # decode here!
    return inner(a[:, 0:-1], array(decoder))


#------------------------------------------------------------------------
# Dual Band Regular Decoders
#
#   panto_dbr           pantophonic
#   peri_dbr            periphonic
#------------------------------------------------------------------------

def panto_dbr(a, num_speakers = 4, orientation = 'flat', Wn = None, zi = None):
    """panto_dbr(a, num_speakers = 4, orientation = 'flat', Wn = None,
                                                                zi = None)
    
    Pantophonic Dual Band Regular Array Decoder
    
    Args:
        - num_speakers  : number of loudspeakers

        - orientation   : "flat" if the front bisects a side of
                            the polygon. The first speaker will
                            be the one left of center. "point" if
                            the front is a vertex of the polygon.
                            The first speaker will be directly in front.

                            Loudspeakers are returned in counter-clockwise
                            orientation.
                            
        - Wn            : shelf filter corner, set to ~400 Hz

        - zi            : Initial conditions for the filter delays, a
                            4 x 2 vector. (4 channels of 2nd order) If
                            zi = None or is not given then initial rest is
                            assumed.  SEE signal.lfiltic for more information.


    Decode a three dimensional ambisonic B-format signal to a set of
    loudspeakers in a regular, horizontal polygon (a ring array). The
    outputs are in counter-clockwise order, and the position of the first
    speaker is either left of centre or centre.

    This decoder employs psychoacoustic shelf filtering.

    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """
    # return decoder gain matrix (velocity decoder)
    decoder = panto_reg_decoder_gain_matrix(num_speakers, orientation, 1)

    # shelf filter
    if zi is None:
        b = psycho_shelf(a, Wn, C.k_2D)

        # decode here!
        res = inner(b, array(decoder))
        return res

    else:
        b, zf = psycho_shelf(a, Wn, C.k_2D, zi)

        # decode here!
        res = inner(b[:, 0:-1], array(decoder))
        return res, zf


def peri_dbr(a, num_speaker_pairs = 4, elevation = 0.6155, orientation = 'flat', \
             Wn = None, zi = None):
    """peri_dbr(a, num_speaker_pairs = 4, elevation = 0.6155,
                                orientation = 'flat', Wn = None, zi = None)
    
    Periphonic Dual Band Regular Array Decoder
    
    Args:
        - num_speaker_pairs  : number of loudspeaker pairss

        - orientation   : "flat" if the front bisects a side of
                            the upper polygon. The first speaker will
                            be the one left of center in the upper
                            polygon. "point" if the front is a vertex of
                            the upper polygon. The first speaker will be
                            upper-front.

                            Loudspeakers are returned in counter-clockwise
                            orientation.
                            
        - Wn            : shelf filter corner, set to ~400 Hz

        - zi            : Initial conditions for the filter delays, a
                            4 x 2 vector. (4 channels of 2nd order) If
                            zi = None or is not given then initial rest is
                            assumed.  SEE signal.lfiltic for more information.


    Decode a three dimensional ambisonic B-format signal to a set of
    loudspeakers in a two regular, horizontal polygons: an upper and
    lower ring array. The outputs are in counter-clockwise order, and
    the position of the first speaker is either left of centre or centre.
    The upper array is returned first, followed by the lower array.

    This decoder employs psychoacoustic shelf filtering.

    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """
    # return decoder gain matrix (velocity decoder)
    decoder = peri_reg_decoder_gain_matrix(num_speaker_pairs, \
                                           elevation, orientation, 1)

    # shelf filter
    if zi is None:
        b = psycho_shelf(a, Wn, C.k_3D)

        # decode here!
        res = inner(b, array(decoder))
        return res

    else:
        b, zf = psycho_shelf(a, Wn, C.k_3D, zi)

        # decode here!
        res = inner(b, array(decoder))
        return res, zf


#------------------------------------------------------------------------
# Dual Band Diametric Decoder
#
#   decode_dbd
#   quad_dbd
#------------------------------------------------------------------------


def decode_dbd(a, positions, Wn, zi = None):
    """decode_dbd(a, positions, Wn, zi = None)
    
    Diametric Pairs Dual Band Decoder (DDT)
    
    Args:
        - a         : Input B-format signal

        - positions : XYZ positions of the speaker pairs,
                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
                      If Z positions of speaker pairs are omitted,
                      it does a horizontal decode (otherwise Z
                      gain is infinite).

                      positions: [ [ x_0   y_0   z_0 ]
                                   ...
                                   [ x_i   y_i   z_i ]
                                   ...
                                   [ x_n-1 y_n-1 z_n-1 ] ]

        - Wn        : shelf filter corner, set to ~400 Hz

        - zi        : Initial conditions for the filter delays, a
                            4 x 2 vector. (4 channels of 2nd order) If
                            zi = None or is not given then initial rest is
                            assumed.  SEE signal.lfiltic for more information.


    Decode a three dimensional ambisonic B-format signal to a set of
    loudspeakers in an array as specified by positions arg.

    This decoder employs psychoacoustic shelf filtering.

    Note:   This assumes standard B format
            definitions for W, X, and Y, i.e., W
            is sqrt(2) lower than X, and Y.

    """

    # return decoder gain matrix (velocity decoder)
    decoder = decoder_gain_matrix(positions, 1)

    # determin k for pantophonic or periphonic shelving
    if shape(positions)[1] is 2:
        k = C.k_2D
    else:
        k = C.k_3D
    
    
    # shelf filter
    if zi is None:
        b = psycho_shelf(a, Wn, k)

        # decode here!
        if shape(positions)[1] is 2:
            b = b[:, 0:-1]
            
        res = inner(b, array(decoder))
        return res

    else:
        b, zf = psycho_shelf(a, Wn, k, zi)

        # decode here!
        if shape(positions)[1] is 2:
            b = b[:, 0:-1]
            
        res = inner(b, array(decoder))
        return res, zf


def quad_dbd(a, angle = 0.7854, Wn = None, zi = None):
    """quad_dbd(a, angle = 0.7854, Wn = None, zi = None)
    
    Quadraphonic Dual Band Decoder
    
    Args:
        - a         : Input B-format signal

        - angle     : Front pair 1/2 angle

        - Wn        : shelf filter corner, set to ~400 Hz

        - zi        : Initial conditions for the filter delays, a
                            4 x 2 vector. (4 channels of 2nd order) If
                            zi = None or is not given then initial rest is
                            assumed.  SEE signal.lfiltic for more information.
                         
    Decode a three dimensional ambisonic B-format signal to a
    quadraphonic (horizontal only) loudspeaker array. The layout
    is adjustable via the angle arg. Returned gains are in counter-
    clockwise order from front centre.

    The outputs will be:
    
         [ Front Left, Back Left, Back Right, Front Right ]


    This decoder employs psychoacoustic shelf filtering.

    """
    # return decoder gain matrix (velocity decoder)
    decoder = quad_decoder_gain_matrix(angle, 1)

    # shelf filter
    if zi is None:
        b = psycho_shelf(a, Wn, C.k_2D)

        # decode here!
        res = inner(b[:, 0:-1], array(decoder))
        return res

    else:
        b, zf = psycho_shelf(a, Wn, C.k_2D, zi)

        # decode here!
        res = inner(b[:, 0:-1], array(decoder))
        return res, zf


#------------------------------------------------------------------------
# Special decoders
#
#   b_to_ITU5           Wiggins coefficients (Single Band Decoders)
#
#------------------------------------------------------------------------

def b_to_ITU5(a, kind = 'foc'):
    """b_to_ITU5(a, kind = 'foc')
    
    Args:
        - a         : Input B-format signal
        - decoder : Three decoders are available. 

                'foc': LF and HF optimised with angular distortions
                'equ': Another potential full range decoder, less angular distortion
                'fou': Four speakers only (no centre), less angular distortion

    Decode a three dimensional ambisonic B-format signal to ITU-5.0.
    Embeded Wiggins coefficients used by agreement between:
    ******* University of Derby & University of Hull
    ******* ONLY to be used for Juan Pampin DVD-A project
    ******* further agreement necessary for incorporation into ATK
    ******* and other uses
    ******* NO further use or disclosure agreed

    The outputs will be in counter-clockwise order: FC, FL, BL, BR, FR. This "farfield" decode
    does not employ shelf filters.

    """

    # LF and HF optimised with angular distortions
    #         C       FL      BL      BR      FR

    # W       0.2000  0.4250  0.4700  0.4700  0.4250
    # C1 (X)  0.1600  0.3600  -0.3300 -0.3300 0.3600
    # S1 (Y)  0.0000  0.4050  0.4150  -0.4150 -0.4050


    # Another potential full range decoder, less angular distortion
    #         C       FL      BL      BR      FR

    # W       0.0000  0.3650  0.5550  0.5550  0.3650
    # C1 (X)  0.0850  0.4350  -0.2850 -0.2850 0.4350
    # S1 (Y)  0.0000  0.3400  0.4050  -0.4050 -0.3400

    # Hi Joseph, here are some more 5.1 coeffs
    #         C       FL      BL      BR      FR

    # W       0.0000  0.4250  0.6300  0.6300  0.4250
    # C1 (X)  0.0000  0.3850  -0.2750 -0.2750 0.3850
    # S1 (Y)  0.0000  0.3300  0.2850  -0.2850 -0.3300

    # define decoders, as a dictionary
    decoder_dict = {
        'foc' : [[.2000, .4250,  .4700,  .4700,  .4250],
                 [.1600, .3600, -.3300, -.3300,  .3600],
                 [.0000, .4050,  .4150, -.4150, -.4050]],

        'equ' : [[.0000, .3650,  .5550,  .5550,  .3650],
                 [.0850, .4350, -.2850, -.2850,  .4350],
                 [.0000, .3400,  .4050, -.4050, -.3400]],

        'fou' : [[.0000, .4250,  .6300,  .6300,  .4250],
                 [.0000, .3850, -.2750, -.2750,  .3850],
                 [.0000, .3300,  .2850, -.2850, -.3300]]
        }

    # construct decoder
    decoder = transpose(array(decoder_dict[kind]))

    # decode here!
    return inner(a[:, 0:-1], array(decoder))


#------------------------------------------------------------------------
# Stereo decoders
#
#   b_to_uhj            "Ambisonic Decoders for HDTV" (1992)
#   b_to_stereo         virtual stereo microphone decoding
#   b_to_binaural       HRTF decoding
#
#------------------------------------------------------------------------

def b_to_uhj(b, decoder_kernels, mode = 'z', kind = 'fft', zi = None):
    """b_to_uhj(b, decoder_kernels, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - b                 : Input B-format signal

        - decoder_kernels   : UHJ decoder kernels:
        
                              [[W_UHJ_l, W_UHJ_r],
                               [X_UHJ_l, X_UHJ_r],
                               [Y_UHJ_l, Y_UHJ_r]]

                               shape = (3, UHJ_kernel_size, 2)

        - mode  : 'z' or 'full'. If mode is 'z', acts as a filter with state
                  'z', and returns a vector of length nframes(b). If mode is
                  'full', returns the full convolution.

        - kind  : 'direct' or 'fft', for direct or fft convolution

        - zi    : Initial state. An array of shape...

                  (3, len(uhj_kernel_size) - 1, 2)
                  
                  If zi = None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y         : The output of the decoder.
      zf        : If zi is None, this is not returned, otherwise, zf holds
                  the filter state.
    
    Decode a three dimensional ambisonic B-format signal to two channel
    UHJ stereo using the supplied UHJ decoder kernels.

    """

    M = nframes(b)                      # length of input b-format
    N = shape(decoder_kernels)[1]       # length of UHJ kernels
    m = 3                               # number of harmonics

    # initialise result to correct size
    if mode is 'z':
        res = zeros((M, 2))
    else:
        res = zeros((M+N-1, 2))

    # convolve with UHJ kernel
    if zi is not None:
        zf = zeros_like(zi)
        
        for i in range(m):
            res_i, zf_i = convfilt(
                interleave(b[:, i]),
                decoder_kernels[i],
                mode,
                kind,
                zi[i]
                )

            res += res_i
            zf[i] = zf_i
    else:
        for i in range(m):
            res += convfilt(
                interleave(b[:, i]),
                decoder_kernels[i],
                mode,
                kind
                )

    if zi is not None:
        return res, zf
    else:
        return res


def b_to_stereo(a, angle = 0.7854, k = 1.0):
    """b_to_stereo(a, angle = 0.7854, k = 1.0)
    
    Virtual Stereo Microphone Decoder
    
    Args:
        - a         : Input B-format signal

        - angle     : Microphone pair 1/2 angle

        - pattern   : Microphone response pattern, k. k = 1 returns
                        a bi-directional pattern (fig-8). k = 0.5
                        returns a cardioid pattern. I.e.:

                        (1-k) + k * cos

                         
    Decode a three dimensional ambisonic B-format signal to two channel
    stereo.


    This decoder does not employ psychoacoustic shelf filtering.

    """
    # compute cosines and sines
    cos_angle, sin_angle = cos(angle), sin(angle)

    # construct decoder gain matrix
    decoder = array([
        [(1.-k) * C.sqrt2, k * cos_angle, k * sin_angle],
        [(1.-k) * C.sqrt2, k * cos_angle, -k * sin_angle]
        ])

    # decode here!
    res = inner(a[:, 0:-1], decoder)
    return res


def b_to_binaural(b, decoder_kernels, mode = 'z', kind = 'fft', zi = None):
    """b_to_binaural(b, decoder_kernels, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - b                 : Input B-format signal

        - decoder_kernels   : HRIR decoder kernels, for periphonic (3D):
        
                              [[W_HRIR_l, W_HRIR_r],
                               [X_HRIR_l, X_HRIR_r],
                               [Y_HRIR_l, Y_HRIR_r],
                               [Z_HRIR_l, Z_HRIR_r]]

                               shape = (4, HRIR_kernel_size, 2)

                               for pantophonic (2D):
        
                              [[W_HRIR_l, W_HRIR_r],
                               [X_HRIR_l, X_HRIR_r],
                               [Y_HRIR_l, Y_HRIR_r]]

                               shape = (3, HRIR_kernel_size, 2)

        - mode  : 'z' or 'full'. If mode is 'z', acts as a filter with state
                  'z', and returns a vector of length nframes(b). If mode is
                  'full', returns the full convolution.

        - kind  : 'direct' or 'fft', for direct or fft convolution

        - zi    : Initial state. An array of shape...

                  for periphonic (3D): (4, len(HRIR_kernel_size) - 1, 2)
                  for pantophonic (2D): (3, len(HRIR_kernel_size) - 1, 2)
                  
                  If zi = None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y         : The output of the decoder.
      zf        : If zi is None, this is not returned, otherwise, zf holds
                  the filter state.
    
    Decode a three dimensional ambisonic B-format signal to two channel
    binaural stereo using the supplied HRIR decoder kernels.

    """

    M = nframes(b)                      # length of input b-format
    N = shape(decoder_kernels)[1]       # length of HRIR kernels
    m = shape(decoder_kernels)[0]       # number of harmonics

    # initialise result to correct size
    # NOTE: At the moment, the assumption is asymmetrical HRIR,
    #       i.e., separate convolutions for L, R
    #
    #       In the symmetrical case, the algorithm can be optimised to
    #       use 3 convolutions for 2D and 4 for 3D. I.e., half the 
    #       convolution load. Such optimisation may be useful for
    #       real-time implementations. HOWEVER, there is evidence
    #       for better perceptial performance to use real-world
    #       assymetric HRIRs
    if mode is 'z':
        res = zeros((M, 2))
    else:
        res = zeros((M+N-1, 2))

    # convolve with HRIR kernel
    if zi is not None:
        zf = zeros_like(zi)
        
        for i in range(m):
            res_i, zf_i = convfilt(
                interleave(b[:, i]),
                decoder_kernels[i],
                mode,
                kind,
                zi[i]
                )

            res += res_i
            zf[i] = zf_i
    else:
        for i in range(m):
            res += convfilt(
                interleave(b[:, i]),
                decoder_kernels[i],
                mode,
                kind
                )

    if zi is not None:
        return res, zf
    else:
        return res


# ------------------------------------------------------------
# B-format to A-format
# ------------------------------------------------------------

def b_to_a(a, orientation = 'flu', weight = 'can'):
    """b_to_a(a, orientation = 'flu', weight = 'can')
    
    Args:
        - a              : Input B-format signal
        - orientation    : Orientation of the A-format channel tetrahedron
            flu = front left up:          FLU, FRD, BLD, BRU
            fld = front left down:        FLD, FRU, BLU, BRD
            flr = front left-right:       FL, FR, BU, BD
            fud = front up-down:          FU, FD, BL, BR
            fbd = front and back down:    F, BD, BLU, BRU 
            fbu = front and back up:      F, BU, BLD, BRD
            flru = front left-right up:   FLU, FRU, FD, B
            flrd = front left-right down: FLD, FRD, FU, B

        - weight : The scaling on W in the output a-format signal:
            can = canonical scaling of 1/sqrt(2)
            dec = scaling for decorrelated soundfields, weight of 1/sqrt(3) on W
            car = scaling for cardioid response, weight of sqrt(3) on W
            uns = unscaled, weight of 1 on W

            dec is the usual choice for use in reverberators

    Transform a B-format signal to the A-format domain.

    """

    decoder_dict = C.b_to_a_dict

    weight_dict = C.b_to_a_weight_dict

    # construct decoder
    decoder = array(decoder_dict[orientation]) * array(weight_dict[weight])

    # decode here!
    return inner(a, decoder)


# ------------------------------------------------------------
# Redundant
# ------------------------------------------------------------

### The function below, decoder_gain_matrix is a
### transcoding of Aaron Heller's Octave code available at:
##
### http://www.ai.sri.com/ajh/ambisonics/
### Aaron J. Heller <heller@ai.sri.com>
##
### Transcoding to Python/Numpy by Joseph Anderson
##
### Calculations for sec 2.4 Loudspeaker Arrays and sec 2.5 Decoding
### Equations of
##
### Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"
### Preprint from AES-121, 10/2006, San Francisco
##
##def decoder_gain_matrix(positions, k):
##    """decoder_gain_matrix(positions, k)
##    
##    Args:
##        - positions : XYZ positions of the speaker pairs,
##                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
##                      If Z positions of speaker pairs are omitted,
##                      it does a horizontal decode (otherwise Z
##                      gain is infinite).
##
##                      positions: [ [ x_0   y_0   z_0 ]
##                                   ...
##                                   [ x_i   y_i   z_i ]
##                                   ...
##                                   [ x_n-1 y_n-1 z_n-1 ] ]
##
##                      Note: radii of all positions must be equal
##
##        - k         : W gain factor, corresponding to directivity
##                      k: 1         => velocity,
##                         sqrt(1/2) => energy 2D, 
##                         sqrt(1/3) => energy 3D, 
##                         1/2       => controlled opposites
##        
##
##    Compute alpha, beta, and gamma. Return values are the weights for X, Y,
##    and Z as columns of the matrix.
##
##    retval:  alpha_0 ... alpha_i  ... alpha_n-1
##             beta_0  ... beta_i   ... beta_n-1
##             gamma_0 ... gamma_i  ... gamma_n-1 
##
##    where the signal for the i'th speaker pair is:
##
##    S_i = W +/- ( alpha_i*X + 
##                  beta_i*Y +
##                  gamma_i*Z )
##
##    Note: This assumes standard B format
##    definitions for W, X, Y, and Z, i.e., W
##    is sqrt(2) lower than X, Y, and Z.
##
##
##    """
##
##    # allow entry of positions as
##    # transpose for convenience
##    # e.g., speaker positions are now in columns
##    # rather than rows
##    positions = transpose(positions)
##
##    # n = number of speaker pairs
##    # m = number of dimensions,
##    #        2=horizontal, 3=periphonic 
##    m, n = shape(positions)     
##
##    # scatter matrix accumulator
##    s = zeros((m, m))
##
##    # speaker directions matrix
##    directions = zeros((m, n))
##
##    for i in range(n):
##
##        # get the i'th speaker position
##        # e.g., select the i'th column
##        pos = positions[:,i]
##
##        # normalize to get direction cosines
##        dir = pos / sqrt(sum(pos**2))
##
##        # form scatter matrix and accumulate
##        s += outer(dir, dir)
##
##        # form matrix of speaker directions
##        directions[:,i] = dir
##
##    res = C.rec_sqrt2 * n * k * dot(matrix(s).I, directions)
##
##    return res


##def pantof(a, num_speakers = 4, orientation = 'flat', directivity = 1):
##    """pantof(a, num_speakers = 4, orientation = 'flat', directivity = 1)
##    
##    Args:
##        - a         : Input B-format signal
##        - num_speakers      : Number of loudspeakers
##        - orientation        : Should be "flat" if the front bisects a side of
##            the polygon. The first speaker will be the one left of center. Should
##            be "point" if the front is a vertex of the polygon. The first speaker
##            will be directly in front.
##        - directivity        : The weighting of the ambisonic decoding equations.
##            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
##            for controlled opposites, in-phase decoding (cardioid). A value of 0
##            will give optimized energy decoding (max rE). Controlled opposites is
##            the option preferred for acoustically live spaces.
##
##    Decode a three dimensional ambisonic B-format signal to a farfield set
##    of loudspeakers in a regular, horizontal polygon (a ring array). The
##    "farfield" decode is preferred for large scale and concert hall decoding.
##    The outputs will be in counter-clockwise order. The position of the first
##    speaker is either left of center or center. The "farfield" decode does not
##    employ shelf filters.
##
##    """
##    # define function to return theta
##    def theta(speaker):
##        if orientation == 'point':
##            theta = ((2. * speaker)/num_speakers) * pi # for 'point' case
##        else :
##            theta = ((1. + (2. * speaker))/num_speakers) * pi # for 'flat' case, default
##
##        return theta
##
##
##    # define constants
##    g0 = 1.
##    g1 = C.rec_sqrt2
##
##    direct = pow(2., (1. - directivity)/2.) # directivity constant
##
##    # calculate decoding matrix
##    decoder = []                # start with empty list
##    speakers = range(num_speakers)
##
##    for speaker in speakers:
##        decoder.append([
##                g0,
##                direct * g1 * cos(theta(speaker)),
##                direct * g1 * sin(theta(speaker)),
##                0.
##                ])
##
##    # decode here!
##    return inner(a, array(decoder))


### *******************************
### DDT version of pantof
### (gives same result as pantof, above)
### *******************************
##
### def pantof(a, num_speaker_pairs = 2, orientation = 'flat', directivity = 1):
###     """pantof(a, num_speaker_pairs = 2, orientation = 'flat', directivity = 1)
##    
###     Args:
###         - a                  : Input B-format signal
###         - num_speaker_pairs  : Number of loudspeaker pairs
###         - orientation        : Should be "flat" if the front bisects a side of
###             the polygon. The first speaker will be the one left of center. Should
###             be "point" if the front is a vertex of the polygon. The first speaker
###             will be directly in front.
###         - directivity        : The weighting of the ambisonic decoding equations.
###             Varies between -1 for idealized, strict soundfield decoding (rV) to +1
###             for controlled opposites, in-phase decoding (cardioid). A value of 0
###             will give optimized energy decoding (max rE). Controlled opposites is
###             the option preferred for acoustically live spaces.
##
###     Decode a three dimensional ambisonic B-format signal to a farfield set
###     of speaker pairs in a regular, even sided, horizontal polygon. The "farfield"
###     decode is preferred for large scale and concert hall decoding. The outputs will
###     be in counter-clockwise order. The position of the first speaker is either left
###     of center or center. The "farfield" decode does not employ shelf filters.
##
###     This decoder is generated via the Diametric Decoder Theorem.
##
###     """
###     # map directivity to k
###     # k = 1          : velocity
###     # k = 1/sqrt(2)  : energy, 2d
###     # k = 1/2        : cardiod
###     k = 1. / pow(2., (1. + directivity) / 2.)
##    
###     # generate speaker pair positions
###     # start with polar positions. . .
###     theta = array([])
###     for spkr_pr in range(num_speaker_pairs):
###         theta = append(theta, pi * spkr_pr / num_speaker_pairs)
##
###     if orientation is 'flat':
###         theta += .5 * pi / num_speaker_pairs
##
###     polar = array([ ones(num_speaker_pairs), theta ]) # [ r, theta ]
##
###     positions = interleave(pol_to_cart(polar)) # . . . then convert from polar to cartesian
##
###     # generate decoder pairs
###     decoder = decoder_gain_matrix(positions, k)
##
###     # append opposing coefficients,
###     # so that decoder now contains all speakers
###     decoder = hstack((decoder, -decoder))
##
###     # add W and Z coefficients
###     decoder = vstack((ones(2 * num_speaker_pairs), decoder, zeros(2 * num_speaker_pairs)))
##
###     # interleave
###     decoder = interleave(decoder)
##
###     # decode here!
###     return inner(a, array(decoder))


##def quadf(a, angle = 0.78539816339744828, directivity = 1):
##    """quadf(a, angle = 0.7854, directivity = 1)
##    
##    Args:
##        - a                  : Input B-format signal
##        - angle              : Front pair 1/2 angle
##        - directivity        : The weighting of the ambisonic decoding equations.
##            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
##            for controlled opposites, in-phase decoding (cardioid). A value of 0
##            will give optimized energy decoding (max rE). Controlled opposites is
##            the option preferred for acoustically live spaces.
##
##    Decode a three dimensional ambisonic B-format signal to a farfield quadraphonic
##    set speaker pairs. The "farfield" decode is preferred for large scale and concert
##    hall decoding. The outputs will be:
##    
##         [ Left Front, Right Front, Left Back, Right Back ]
##
##    The "farfield" decode does not employ shelf filters.
##
##    This decoder is generated via the Diametric Decoder Theorem.
##
##    """
##    # map directivity to k
##    # k = 1          : velocity
##    # k = 1/sqrt(2)  : energy, 2d
##    # k = 1/2        : cardiod
##    k = 1. / pow(2., (1. + directivity) / 2.)
##    
##    # generate speaker pair positions
##    # start with polar positions. . .
##    num_speaker_pairs = 2
##
##    polar = array([[1, angle], [1, -angle]]) # [ r, theta ]
##    positions = pol_to_cart(polar) # . . . then convert from polar to cartesian
##
##    # generate decoder pairs
##    decoder = decoder_gain_matrix(positions, k)
##
##    # append opposing coefficients,
##    # so that decoder now contains all speakers
##    decoder = hstack(
##        (
##            decoder,
##            -roll(decoder, 1, 1)
##            )
##        )
##
##    # add W and Z coefficients
##    decoder = vstack((ones(2 * num_speaker_pairs), decoder, zeros(2 * num_speaker_pairs)))
##
##    # interleave
##    decoder = interleave(decoder)
##
##    # decode here!
##    return inner(a, array(decoder))


##def perif(a, num_speaker_pairs = 4, orientation = 'flat', elevation = 0.61547971, directivity = 1):
##    """perif(a, num_speaker_pairs = 4, orientation = 'flat', elevation = 0.61547971, directivity = 1)
##    
##    Args:
##        - a                  : Input B-format signal
##        - num_speaker_pairs  : Number of loudspeaker PAIRS
##        - orientation        : Should be "flat" if the front of the upper ring
##            bisects a side of the upper polygon. The first speaker will be the
##            one left of center in the upper ring. Should be "point" if the front
##            of the upper ring is a vertex of the polygon. The first speaker
##            will be directly in front of the upper ring.
##        - elevation          : The elevation half-angle in radians. The default, 
##            0.61547971, specifies the half-angle elevation for a cube
##        - directivity        : The weighting of the ambisonic decoding equations.
##            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
##            for controlled opposites, in-phase decoding (cardioid). A value of 0
##            will give optimized energy decoding (max rE). Controlled opposites is
##            the option preferred for acoustically live spaces.
##
##    Decode a three dimensional ambisonic B-format signal to a farfield set
##    of loudspeakers in two regular polygons, one above and one below. The "farfield"
##    decode is preferred for large scale and concert hall decoding. The outputs will
##    be in counter-clockwise order, beginning with the upper ring. The position of the
##    first speaker is either left of center upper or center upper. The "farfield"
##    decode does not employ shelf filters.
##
##    This decoder is generated via the Diametric Decoder Theorem.
##
##    """
##    # map directivity to k (simple quadratic mapping)
##    # k = 1          : velocity
##    # k = 1/sqrt(3)  : energy, 3d
##    # k = 1/2        : cardiod
##    c_c = 1. / sqrt(3)
##    c_b = -.25
##    c_a = .75 - c_c
##    k = c_a * directivity**2 + c_b * directivity + c_c
##    
##    # generate speaker pair positions
##    # start with spherical positions. . .
##    theta = array([])
##    for spkr_pr in range(num_speaker_pairs):
##        theta = append(theta, C.twoPi * spkr_pr / num_speaker_pairs)
##
##    if orientation is 'flat':
##        theta += pi / num_speaker_pairs
##
##    spher = array([
##            ones(num_speaker_pairs),
##            theta,
##            repeat(elevation, num_speaker_pairs)
##            ]) # [ r, theta, phi ]
##
##    positions = interleave(spher_to_cart(spher)) # . . . then convert from spherical to cartesian
##
##    # generate decoder pairs
##    # (positive output gives upper ring)
##    decoder = decoder_gain_matrix(positions, k)
##
##    # append opposing coefficients,
##    # so that decoder now contains all speakers
##    if orientation is 'point':
##        decoder = hstack(       # point. . .
##            (
##                decoder,                                 # upper ring
##                -roll(decoder, num_speaker_pairs / 2, 1) # lower ring (opposing pairs)
##                )                                        # rolled so first lower speaker
##            )                                            # is below first upper speaker
##    else:
##        decoder = hstack(       # . . . or flat
##            (
##                decoder,                                 # upper ring
##                -roll(decoder, (1 + num_speaker_pairs) / 2, 1) # lower ring (opposing pairs)
##                )                                        # rolled so first lower speaker
##            )                                            # is below first upper speaker
##
##    # add W coefficients
##    decoder = vstack((ones(2 * num_speaker_pairs), decoder))
##
##    # interleave
##    decoder = interleave(decoder)
##
##    # decode here!
##    return inner(a, array(decoder))


##def diametricf(a, positions, directivity = 1):
##    """diametricf(a, positions, directivity = 1)
##    
##    Args:
##        - a                  : Input B-format signal
##        - positions : XYZ positions of the speaker pairs,
##                      one speaker pair per row, i.e., [[1 1 1], [1 -1 -1]]
##                      If Z positions of speaker pairs are omitted,
##                      it does a horizontal decode (otherwise Z
##                      gain is infinite).
##
##                      positions: [ [ x_0   y_0   z_0 ]
##                                   ...
##                                   [ x_i   y_i   z_i ]
##                                   ...
##                                   [ x_n-1 y_n-1 z_n-1 ] ]
##
##                      Note: radii of all positions must be equal
##
##        - directivity        : The weighting of the ambisonic decoding equations.
##            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
##            for controlled opposites, in-phase decoding (cardioid). A value of 0
##            will give optimized energy decoding (max rE). Controlled opposites is
##            the option preferred for acoustically live spaces.
##
##    Decode a three dimensional ambisonic B-format signal to a farfield set
##    of loudspeakers specified in opposing speaker pairs. The "farfield"
##    decode is preferred for large scale and concert hall decoding. The outputs will
##    be in input position order, followed by opposing pairs, in the same order. It is
##    left to the user distribute output channels as appropriate.
##
##    The "farfield" decode does not employ shelf filters.
##
##    This decoder is generated via the Diametric Decoder Theorem.
##
##    """
##    # read number of speaker pairs
##    num_speaker_pairs = shape(positions)[0]
##
##    # test for periphonic (3D) or pantophonic (2D) decode
##    if shape(positions)[-1] is 3:
##        periphonic = True
##    else:
##        periphonic = False
##
##    # map directivity to k
##    if periphonic:
##        # (simple quadratic mapping)
##        # k = 1          : velocity
##        # k = 1/sqrt(3)  : energy, 3d
##        # k = 1/2        : cardiod
##        c_c = 1. / sqrt(3)
##        c_b = -.25
##        c_a = .75 - c_c
##        k = c_a * directivity**2 + c_b * directivity + c_c
##    
##    else:
##        # k = 1          : velocity
##        # k = 1/sqrt(2)  : energy, 2d
##        # k = 1/2        : cardiod
##        k = 1. / pow(2., (1. + directivity) / 2.)
##    
##    # generate decoder pairs
##    decoder = decoder_gain_matrix(positions, k)
##
##    # append opposing coefficients,
##    # so that decoder now contains all speakers
##    decoder = hstack((decoder, -decoder)) # ( specified positions, opposing pairs)
##
##    # add W coefficients, and Z, if necessary
##    if periphonic: 
##        decoder = vstack((ones(2 * num_speaker_pairs), decoder))
##    else:
##        decoder = vstack((ones(2 * num_speaker_pairs), decoder, zeros(2 * num_speaker_pairs)))
##
##    # interleave
##    decoder = interleave(decoder)
##
##    print decoder
##
##    # decode here!
##    return inner(a, array(decoder))
