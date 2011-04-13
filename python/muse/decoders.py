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


# #=========================
# # Definition of constants
# #=========================


#=========================
# Functions
#=========================

# The function below, decoder_gain_matrix is a
# transcoding of Aaron Heller's Octave code available at:

# http://www.ai.sri.com/ajh/ambisonics/
# Aaron J. Heller <heller@ai.sri.com>

# Transcoding to Python/Numpy by Joseph Anderson

# Calculations for sec 2.4 Loudspeaker Arrays and sec 2.5 Decoding
# Equations of

# Benjamin, et al., "Localization in Horizontal-Only Ambisonic Systems"
# Preprint from AES-121, 10/2006, San Francisco

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

                      Note: radii of all positions must be equal

        - k         : W gain factor, corresponding to directivity
                      k: 1         => velocity,
                         sqrt(1/2) => energy 2D, 
                         sqrt(1/3) => energy 3D, 
                         1/2       => controlled opposites
        

    Compute alpha, beta, and gamma. Return values are the weights for X, Y,
    and Z as columns of the matrix.

    retval:  alpha_0 ... alpha_i  ... alpha_n-1
             beta_0  ... beta_i   ... beta_n-1
             gamma_0 ... gamma_i  ... gamma_n-1 

    where the signal for the i'th speaker pair is:

    S_i = W +/- ( alpha_i*X + 
                  beta_i*Y +
                  gamma_i*Z )

    Note: This assumes standard B format
    definitions for W, X, Y, and Z, i.e., W
    is sqrt(2) lower than X, Y, and Z.


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

    res = C.rec_sqrt2 * n * k * dot(matrix(s).I, directions)

    return res


def pantof(a, num_speakers = 4, orientation = 'flat', directivity = 1):
    """pantof(a, num_speakers = 4, orientation = 'flat', directivity = 1)
    
    Args:
        - a         : Input B-format signal
        - num_speakers      : Number of loudspeakers
        - orientation        : Should be "flat" if the front bisects a side of
            the polygon. The first speaker will be the one left of center. Should
            be "point" if the front is a vertex of the polygon. The first speaker
            will be directly in front.
        - directivity        : The weighting of the ambisonic decoding equations.
            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
            for controlled opposites, in-phase decoding (cardioid). A value of 0
            will give optimized energy decoding (max rE). Controlled opposites is
            the option preferred for acoustically live spaces.

    Decode a three dimensional ambisonic B-format signal to a farfield set
    of loudspeakers in a regular, horizontal polygon (a ring array). The
    "farfield" decode is preferred for large scale and concert hall decoding.
    The outputs will be in counter-clockwise order. The position of the first
    speaker is either left of center or center. The "farfield" decode does not
    employ shelf filters.

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
    g1 = C.rec_sqrt2

    direct = pow(2., (1. - directivity)/2.) # directivity constant

    # calculate decoding matrix
    decoder = []                # start with empty list
    speakers = range(num_speakers)

    for speaker in speakers:
        decoder.append([
                g0,
                direct * g1 * cos(theta(speaker)),
                direct * g1 * sin(theta(speaker)),
                0.
                ])

    # decode here!
    return inner(a, array(decoder))


# *******************************
# DDT version of pantof
# (gives same result as pantof, above)
# *******************************

# def pantof(a, num_speaker_pairs = 2, orientation = 'flat', directivity = 1):
#     """pantof(a, num_speaker_pairs = 2, orientation = 'flat', directivity = 1)
    
#     Args:
#         - a                  : Input B-format signal
#         - num_speaker_pairs  : Number of loudspeaker pairs
#         - orientation        : Should be "flat" if the front bisects a side of
#             the polygon. The first speaker will be the one left of center. Should
#             be "point" if the front is a vertex of the polygon. The first speaker
#             will be directly in front.
#         - directivity        : The weighting of the ambisonic decoding equations.
#             Varies between -1 for idealized, strict soundfield decoding (rV) to +1
#             for controlled opposites, in-phase decoding (cardioid). A value of 0
#             will give optimized energy decoding (max rE). Controlled opposites is
#             the option preferred for acoustically live spaces.

#     Decode a three dimensional ambisonic B-format signal to a farfield set
#     of speaker pairs in a regular, even sided, horizontal polygon. The "farfield"
#     decode is preferred for large scale and concert hall decoding. The outputs will
#     be in counter-clockwise order. The position of the first speaker is either left
#     of center or center. The "farfield" decode does not employ shelf filters.

#     This decoder is generated via the Diametric Decoder Theorem.

#     """
#     # map directivity to k
#     # k = 1          : velocity
#     # k = 1/sqrt(2)  : energy, 2d
#     # k = 1/2        : cardiod
#     k = 1. / pow(2., (1. + directivity) / 2.)
    
#     # generate speaker pair positions
#     # start with polar positions. . .
#     theta = array([])
#     for spkr_pr in range(num_speaker_pairs):
#         theta = append(theta, pi * spkr_pr / num_speaker_pairs)

#     if orientation is 'flat':
#         theta += .5 * pi / num_speaker_pairs

#     polar = array([ ones(num_speaker_pairs), theta ]) # [ r, theta ]

#     positions = interleave(pol_to_cart(polar)) # . . . then convert from polar to cartesian

#     # generate decoder pairs
#     decoder = decoder_gain_matrix(positions, k)

#     # append opposing coefficients,
#     # so that decoder now contains all speakers
#     decoder = hstack((decoder, -decoder))

#     # add W and Z coefficients
#     decoder = vstack((ones(2 * num_speaker_pairs), decoder, zeros(2 * num_speaker_pairs)))

#     # interleave
#     decoder = interleave(decoder)

#     # decode here!
#     return inner(a, array(decoder))


def quadf(a, angle = 0.78539816339744828, directivity = 1):
    """quadf(a, angle = 0.7854, directivity = 1)
    
    Args:
        - a                  : Input B-format signal
        - angle              : Front pair 1/2 angle
        - directivity        : The weighting of the ambisonic decoding equations.
            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
            for controlled opposites, in-phase decoding (cardioid). A value of 0
            will give optimized energy decoding (max rE). Controlled opposites is
            the option preferred for acoustically live spaces.

    Decode a three dimensional ambisonic B-format signal to a farfield quadraphonic
    set speaker pairs. The "farfield" decode is preferred for large scale and concert
    hall decoding. The outputs will be:
    
         [ Left Front, Right Front, Left Back, Right Back ]

    The "farfield" decode does not employ shelf filters.

    This decoder is generated via the Diametric Decoder Theorem.

    """
    # map directivity to k
    # k = 1          : velocity
    # k = 1/sqrt(2)  : energy, 2d
    # k = 1/2        : cardiod
    k = 1. / pow(2., (1. + directivity) / 2.)
    
    # generate speaker pair positions
    # start with polar positions. . .
    num_speaker_pairs = 2

    polar = array([[1, angle], [1, -angle]]) # [ r, theta ]
    positions = pol_to_cart(polar) # . . . then convert from polar to cartesian

    # generate decoder pairs
    decoder = decoder_gain_matrix(positions, k)

    # append opposing coefficients,
    # so that decoder now contains all speakers
    decoder = hstack(
        (
            decoder,
            -roll(decoder, 1, 1)
            )
        )

    # add W and Z coefficients
    decoder = vstack((ones(2 * num_speaker_pairs), decoder, zeros(2 * num_speaker_pairs)))

    # interleave
    decoder = interleave(decoder)

    # decode here!
    return inner(a, array(decoder))


def perif(a, num_speaker_pairs = 4, orientation = 'flat', elevation = 0.61547971, directivity = 1):
    """perif(a, num_speaker_pairs = 4, orientation = 'flat', elevation = 0.61547971, directivity = 1)
    
    Args:
        - a                  : Input B-format signal
        - num_speaker_pairs  : Number of loudspeaker PAIRS
        - orientation        : Should be "flat" if the front of the upper ring
            bisects a side of the upper polygon. The first speaker will be the
            one left of center in the upper ring. Should be "point" if the front
            of the upper ring is a vertex of the polygon. The first speaker
            will be directly in front of the upper ring.
        - elevation          : The elevation half-angle in radians. The default, 
            0.61547971, specifies the half-angle elevation for a cube
        - directivity        : The weighting of the ambisonic decoding equations.
            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
            for controlled opposites, in-phase decoding (cardioid). A value of 0
            will give optimized energy decoding (max rE). Controlled opposites is
            the option preferred for acoustically live spaces.

    Decode a three dimensional ambisonic B-format signal to a farfield set
    of loudspeakers in two regular polygons, one above and one below. The "farfield"
    decode is preferred for large scale and concert hall decoding. The outputs will
    be in counter-clockwise order, beginning with the upper ring. The position of the
    first speaker is either left of center upper or center upper. The "farfield"
    decode does not employ shelf filters.

    This decoder is generated via the Diametric Decoder Theorem.

    """
    # map directivity to k (simple quadratic mapping)
    # k = 1          : velocity
    # k = 1/sqrt(3)  : energy, 3d
    # k = 1/2        : cardiod
    c_c = 1. / sqrt(3)
    c_b = -.25
    c_a = .75 - c_c
    k = c_a * directivity**2 + c_b * directivity + c_c
    
    # generate speaker pair positions
    # start with spherical positions. . .
    theta = array([])
    for spkr_pr in range(num_speaker_pairs):
        theta = append(theta, C.twoPi * spkr_pr / num_speaker_pairs)

    if orientation is 'flat':
        theta += pi / num_speaker_pairs

    spher = array([
            ones(num_speaker_pairs),
            theta,
            repeat(elevation, num_speaker_pairs)
            ]) # [ r, theta, phi ]

    positions = interleave(spher_to_cart(spher)) # . . . then convert from spherical to cartesian

    # generate decoder pairs
    # (positive output gives upper ring)
    decoder = decoder_gain_matrix(positions, k)

    # append opposing coefficients,
    # so that decoder now contains all speakers
    if orientation is 'point':
        decoder = hstack(       # point. . .
            (
                decoder,                                 # upper ring
                -roll(decoder, num_speaker_pairs / 2, 1) # lower ring (opposing pairs)
                )                                        # rolled so first lower speaker
            )                                            # is below first upper speaker
    else:
        decoder = hstack(       # . . . or flat
            (
                decoder,                                 # upper ring
                -roll(decoder, (1 + num_speaker_pairs) / 2, 1) # lower ring (opposing pairs)
                )                                        # rolled so first lower speaker
            )                                            # is below first upper speaker

    # add W coefficients
    decoder = vstack((ones(2 * num_speaker_pairs), decoder))

    # interleave
    decoder = interleave(decoder)

    # decode here!
    return inner(a, array(decoder))


def diametricf(a, positions, directivity = 1):
    """diametricf(a, positions, directivity = 1)
    
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

                      Note: radii of all positions must be equal

        - directivity        : The weighting of the ambisonic decoding equations.
            Varies between -1 for idealized, strict soundfield decoding (rV) to +1
            for controlled opposites, in-phase decoding (cardioid). A value of 0
            will give optimized energy decoding (max rE). Controlled opposites is
            the option preferred for acoustically live spaces.

    Decode a three dimensional ambisonic B-format signal to a farfield set
    of loudspeakers specified in opposing speaker pairs. The "farfield"
    decode is preferred for large scale and concert hall decoding. The outputs will
    be in input position order, followed by opposing pairs, in the same order. It is
    left to the user distribute output channels as appropriate.

    The "farfield" decode does not employ shelf filters.

    This decoder is generated via the Diametric Decoder Theorem.

    """
    # read number of speaker pairs
    num_speaker_pairs = shape(positions)[0]

    # test for periphonic (3D) or pantophonic (2D) decode
    if shape(positions)[-1] is 3:
        periphonic = True
    else:
        periphonic = False

    # map directivity to k
    if periphonic:
        # (simple quadratic mapping)
        # k = 1          : velocity
        # k = 1/sqrt(3)  : energy, 3d
        # k = 1/2        : cardiod
        c_c = 1. / sqrt(3)
        c_b = -.25
        c_a = .75 - c_c
        k = c_a * directivity**2 + c_b * directivity + c_c
    
    else:
        # k = 1          : velocity
        # k = 1/sqrt(2)  : energy, 2d
        # k = 1/2        : cardiod
        k = 1. / pow(2., (1. + directivity) / 2.)
    
    # generate decoder pairs
    decoder = decoder_gain_matrix(positions, k)

    # append opposing coefficients,
    # so that decoder now contains all speakers
    decoder = hstack((decoder, -decoder)) # ( specified positions, opposing pairs)

    # add W coefficients, and Z, if necessary
    if periphonic: 
        decoder = vstack((ones(2 * num_speaker_pairs), decoder))
    else:
        decoder = vstack((ones(2 * num_speaker_pairs), decoder, zeros(2 * num_speaker_pairs)))

    # interleave
    decoder = interleave(decoder)

    print decoder

    # decode here!
    return inner(a, array(decoder))


def b_to_ITU5(a, kind = 'foc'):
    """b_to_ITU5(a, kind = 'foc')
    
    Args:
        - a         : Input B-format signal
        - decoder : Two decoders are available. 

                'foc': LF and HF optimised with angular distortions
                'equ': Another potential full range decoder, less angular distortion

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


    # define decoders, as a dictionary
    decoder_dict = {
        'foc' : [[.2000, .4250,  .4700,  .4700,  .4250],
                 [.1600, .3600, -.3300, -.3300,  .3600],
                 [.0000, .4050,  .4150, -.4150, -.4050],
                 [.0000, .0000,  .0000,  .0000,  .0000]],

        'equ' : [[.0000, .3650,  .5550,  .5550,  .3650],
                 [.0850, .4350, -.2850, -.2850,  .4350],
                 [.0000, .3400,  .4050, -.4050, -.3400],
                 [.0000, .0000,  .0000,  .0000,  .0000]]
        }

    # construct decoder
    decoder = transpose(array(decoder_dict[kind]))

    # decode here!
    return inner(a, array(decoder))


def b_to_uhj(a, N, beta = 5, mode = 'z', kind = 'fft', zi = None):
    """b_to_uhj(a, N, beta = 5, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - a    -- Input B-format signal
        - N    -- filter (number of taps), should be odd
        - beta -- beta for Kaiser window FIR design.
                  5 = similiar to Hamming.
        - mode -- 'z' or 'full'. If mode is 'z', acts as a filter
                  with state 'z', and returns a vector of length
                  len(x). If mode is 'full', returns the full
                  convolution.
        - kind -- 'direct' or 'fft', for direct or fft convolution

        - zi   -- Initial state. An array of shape (len(kernel) - 1, 3).
                  If zi=None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y -- The output of the decoder.
      zf -- If zi is None, this is not returned, otherwise, zf holds the
            filter state.
    
    Decode a three dimensional ambisonic B-format signal to two channel
    stereo ambisonic UHJ using a linear phase hilbert transform filter.

    """
#     Audio Engineering Society E-Library
#     Ambisonic Decoders for HDTV
#     Preprint Number:   3345    Convention:   92 (February 1992)
#     Authors:   Gerzon, Michael A.; Barton, Geoffrey J.
#     E-library Location: (CD aes12)   /pp9193/pp9203/3405.pdf

#     S = 0.9396926*W + 0.1855740*X
#     D = j(-0.3420201*W + 0.5098604*X) + 0.6554516*Y
#     Left = (S + D)/2.0
#     Right = (S - D)/2.0

    # convolve with hilbert kernel
    if zi is not None:
        hb, zf = convfilt(
            a[:, :3],               # strip Z
            fir_hb(N, beta=5),      # generate hilbert kernel
            mode,
            kind,
            zi
            )
        over_dub(hb, zi)
    else:
        hb = convfilt(
            a[:, :3],               # strip Z
            fir_hb(N, beta=5),      # generate hilbert kernel
            mode,
            kind
            )

    s = (array([.9396926, .1855740]) * hb.real[:, :2]).sum(axis = -1)
    d = (array([-.3420201, .5098604]) * hb.imag[:, :2]).sum(axis = -1) + \
        .6554516 * hb.real[:, 2]
    
    res = .5 * interleave(array([s + d, s - d]))

    if zi is not None:
        return res, zf
    else:
        return res


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
