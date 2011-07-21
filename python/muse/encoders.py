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

# import muse defined constants
import constants as C



# #=========================
# # Definition of constants
# #=========================


#=========================
# Functions
#=========================

def mono_to_b(a, azimuth = 0., elevation = 0.):
    """mono_to_b(a, azimuth = 0., elevation = 0.)
    
    Args:
        - a         : Input mono signal
        - azimuth   : counter-clockwise, in radians
        - elevation : upwards, in radian

    Encode a mono signal into the B-format domain.

    """

    # compute cosines and sines
    cos_azim, sin_azim = cos(azimuth), sin(azimuth)
    cos_elev, sin_elev = cos(elevation), sin(elevation)

    # compute scalars
    w_scale = C.rec_sqrt2
    x_scale = cos_elev * cos_azim
    y_scale = cos_elev * sin_azim
    z_scale = sin_elev

    # for testing for and constructing vectors
    azim_scalar = isscalar(azimuth)
    elev_scalar = isscalar(elevation)
    n = len(a)

    # construct appropriate encoder
    if not elev_scalar:         # case 2, 4: x, y, z vectors
        encoder = interleave(
            array([
                repeat(w_scale, n),
                x_scale,
                y_scale,
                z_scale
            ]))
    elif not azim_scalar:       # case 3: x, y vectors
        encoder = interleave(
            array([
                repeat(w_scale, n),
                x_scale,
                y_scale,
                repeat(z_scale, n)
            ]))
    else:                       # case 1: all scalars
        encoder = array([
            w_scale,
            x_scale,
            y_scale,
            z_scale
            ])

    # copy input mono array into four channel array
    b = repeat(a, 4).reshape(-1,4)

    # encode here!
    return b * encoder


#---------------------------------------------
# A to B encoder
#---------------------------------------------
def a_to_b(a, orientation = 'flu', weight = 'can'):
    """a_to_b(a, orientation = 'flu', weight = 'can')
    
    Args:
        - a              : Input A-format signal
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

    Transform an A-format signal to the B-format domain.

    """
    decoder_dict = C.b_to_a_dict

    weight_dict = C.b_to_a_weight_dict

    # construct encoder
    encoder = (array(decoder_dict[orientation]) * \
               reciprocal(array(weight_dict[weight]))).transpose()

    # encode here!
    return inner(a, encoder)


#---------------------------------------------
# UHJ and SuperStereo encoder kernels
#---------------------------------------------

def uhj_encoder_kernel(N, Wn, k = 0.):
    """uhj_encoder_kernel(N, Wn, k = 0.)
    
    Generate a filter kernel suitable for UHJ to b-format encoding.
    The kernel is shelf filtered, with the assumption the resulting
    b-format signal will be decoded to a pantophonic (horizontal)
    array with the decoder using an appropriate psychoacoustic
    shelf filter.

    See:

    Audio Engineering Society E-Library
    Ambisonic Decoders for HDTV
    Preprint Number:   3345    Convention:   92 (February 1992)
    Authors:   Gerzon, Michael A.; Barton, Geoffrey J.
    E-library Location: (CD aes12)   /pp9193/pp9203/3405.pdf

    Args:
        - N     : order of filter (number of taps)
        - Wn    : shelf filter corner, set to ~400 Hz
        - k     : Forward preference, pushes phasiness to the rear
                  of the sound field where it is less offensive. With
                  no forward preference (= 0), the phasiness is distributed
                  equally around 360 degrees. 0. < k < .7

    Outputs:
    
        - b         : coefficients of length N FIR filter: 
                        [[ left_W_FIR,  left_X_FIR,  left_Y_FIR],
                         [right_W_FIR, right_X_FIR, right_Y_FIR]]

    Joseph Anderson <josephlloydanderson@mac.com>

    """

    #     S = (Left + Right)
    #     D = (Left - Right)

    #     W' = 0.982*S + j*0.164*D
    #     X' = 0.419*S - j*0.828*D
    #     Y' = 0.763*D + j*0.385*S
    #     B' = 0.116*D - j*0.694*S

    #     W" = k1*W'
    #     X" = k2*X'
    #     Y" = k2*Y' + k'*k3*B'         (where k' = forward preference)

    #   ------------------------
    #           LF      HF              (decoding shelfs)
    #   k1      0.646   1.000           NOTE: These shelfs gains should be 
    #   k2      1.263   1.000                   compensated for on
    #   k3      0.775   1.000                   encoding into b-format.
    #                                           Use a 2D psycho-acoustic
    #                                           shelf filter to do so.

    # Martin Lese on Forward Preference:
    # (http://members.tripod.com/martin_leese/Ambisonic/diy.html)
    #
    # Now, what value of "a" to use?  In theory, 0 < a < 1.  Gerzon's
    # patent 4081606 suggests 0.333 < a < 0.5.  Gerzon 1977 (refs at end)
    # uses a = 0.3.  Gerzon 1985 suggests 0 < k' < 0.7 which is equivalent
    # to 0 < a < 0.49.  I would either allow 0 < a < 1.0 and experiment or
    # go with the most recent reference, Gerzon 1985.


    # J Anderson on UHJ encode with Forward Preference:
    # 
    # The above can be combined to give the following: 
    #
    #     W = k1*(0.982*S + j*0.164*D)
    #     X = k2*(0.419*S - j*0.828*D)
    #     Y = k2*((0.763 + k'*k3/k2*0.116)*D + j*(0.385 - k'*k3/k2*0.694)*S)
    # 
    # k' = sqrt(2)/4 for 'optimum', e.g., no phasiness at front

    #---------------------------------
    # UHJ coefficients

    c_0 = 0.982
    c_1 = 0.164
    c_2 = 0.419
    c_3 = 0.828
    c_4 = 0.763
    c_5 = 0.116
    c_6 = 0.385
    c_7 = 0.694

    k_1 = 0.646                         # LF gains
    k_2 = 1.263                         # HF gains = 1
    k_3 = 0.775

    gains_lf = array([                  # LF complex gains
        [k_1*complex(c_0, c_1),
         k_2*complex(c_2, -c_3),
         k_2*complex((c_4 + (k*k_3/k_2)*c_5), (c_6 - (k*k_3/k_2)*c_7))],

        [k_1*complex(c_0, -c_1),
         k_2*complex(c_2, c_3),
         k_2*complex(-(c_4 + (k*k_3/k_2)*c_5), (c_6 - (k*k_3/k_2)*c_7))]
        ])

    gains_hf = array([                      # HF complex gains, normalised to
        [1./C.k_2D[0]*complex(c_0, c_1),    # 2D psycho_acoustic shelf gains
         1./C.k_2D[1]*complex(c_2, -c_3),
         1./C.k_2D[1]*complex((c_4 + k*c_5), (c_6 - k*c_7))],

        [1./C.k_2D[0]*complex(c_0, -c_1),
         1./C.k_2D[1]*complex(c_2, c_3),
         1./C.k_2D[1]*complex(-(c_4 + k*c_5), (c_6 - k*c_7))]
        ])

    m = 3               # harmonics (W, X, Y)

    #---------------------------------
    # calculate kernels

    hilbert = fir_hb(N)

    encoder_kernels = zeros((2, N, m))  # stereo, N, harmonics

    # collect decoder kernel... shelving
    for i in range(2):          # i is UHJ stereo channel number
        encoder_kernels[i] += gains_lf[i].real * fiir_rmhs2(
            interleave(hilbert.real),
            Wn,
            0.
            )
        encoder_kernels[i] += gains_lf[i].imag * fiir_rmhs2(
            interleave(hilbert.imag),
            Wn,
            0.
            )

        encoder_kernels[i] += gains_hf[i].real * fiir_rmls2(
            interleave(hilbert.real),
            Wn,
            0.
            )
        encoder_kernels[i] += gains_hf[i].imag * fiir_rmls2(
            interleave(hilbert.imag),
            Wn,
            0.
            )

    return encoder_kernels


def super_encoder_kernel(N, width = .593, k = 0.):
    """super_encoder_kernel(N, width = .593, k = 0.)
    
    Generate a filter kernel suitable for "super stereo" to b-format encoding.

    See:

    G. Barton, "RE: [Sursound] UHJ Coding/Decoding," 19-Dec-2003.

    Args:
        - N     : order of filter (number of taps)
        - width : Width setting, the optimal value of which is about .593.
        - k     : Forward preference, pushes phasiness to the rear
                  of the sound field where it is less offensive. With
                  no forward preference (= 0), the phasiness is distributed
                  equally around 360 degrees. .354 gives no phasiness at
                  front. 0. < k < .7

    Outputs:
    
        - b         : coefficients of length N FIR filter: 
                        [[ left_W_FIR,  left_X_FIR,  left_Y_FIR],
                         [right_W_FIR, right_X_FIR, right_Y_FIR]]

    Joseph Anderson <josephlloydanderson@mac.com>

    """

    #     S = (Left + Right)
    #     D = (Left - Right)

    # Geoff Barton on SuperStereo:
    # 
    # as you say; it was already out of date when it was published.
    # There were various later 'stereo decodes' we tried. By about 1982
    # we were using something like:-
    # 
    # W' = 0.6098637*S - 0.6896511*j*w*D
    # X' = 0.8624776*S + 0.7626955*j*w*D
    # Y' = 1.6822415*w*D - 0.2156194*j*S
    # 
    # where 'w' is a width setting, the optimal (in some senses) value
    # of which is about 0.593, S & D are as defined above. NB, in MAG's
    # notation W'' etc is the signal after the shelf filter, W' before.
    # 
    # There are various options, -j can be substituted for j and
    # 'forward preference' can be added.


    # Martin Lese on Forward Preference:
    # 
    # It pushes phasiness to the rear of the sound field where it is less 
    # offensive.  With no forward preference (a = 0), the phasiness is 
    # distributed equally around 360 degrees.  The transformation is:
    # 
    # new W = W
    # new X = X
    # new Y = Y - a.(jW)
    # 
    # where a = a constant
    # j = 90-degree phase shift
    # 
    # Now, what value of "a" to use?  In theory, 0 < a < 1.  Gerzon's
    # patent 4081606 suggests 0.333 < a < 0.5.  Gerzon 1977 (refs at end)
    # uses a = 0.3.  Gerzon 1985 suggests 0 < k' < 0.7 which is equivalent
    # to 0 < a < 0.49.  I would either allow 0 < a < 1.0 and experiment or
    # go with the most recent reference, Gerzon 1985.

    # J Anderson on Super Stereo encode with Forward Preference:
    # 
    # The above can be combined to give the following: 
    #
    #     W = 0.6098637*S - j*0.6896511*w*D
    #     X = 0.8624776*S + j*0.7626955*w*D
    #     Y = (1.6822415 + (k**2)*0.6896511)*w*D
    #           - j*(0.2156194 - (k**2)*0.6098637)*S
    # 
    # k' = sqrt(2)/4 for 'optimum', e.g., no phasiness at front (k' = sqrt(a))
    # Substitute -a to get expected behaviour.

    #---------------------------------
    # Super Stereo coefficients

    c_0 = 0.6098637
    c_1 = 0.6896511
    c_2 = 0.8624776
    c_3 = 0.7626955
    c_4 = 1.6822415
    c_5 = 0.2156194

    gains = array([
        [complex(c_0, -c_1*width),
         complex(c_2, c_3*width),
         complex(((c_4 + (k**2)*c_1)*width), -(c_5 - (k**2)*c_0))],

        [complex(c_0, c_1*width),
         complex(c_2, -c_3*width),
         complex(-((c_4 + (k**2)*c_1)*width), -(c_5 - (k**2)*c_0))]
        ])

    m = 3               # harmonics (W, X, Y)

    #---------------------------------
    # calculate kernels

    hilbert = fir_hb(N)

    encoder_kernels = zeros((2, N, m))  # stereo, N, harmonics

    # collect decoder kernel
    for i in range(2):          # i is stereo channel number
        encoder_kernels[i] += gains[i].real * interleave(hilbert.real)
        encoder_kernels[i] += gains[i].imag * interleave(hilbert.imag)

    return encoder_kernels


#---------------------------------------------
# UHJ, SuperStereo and SimpleStereo encoders
#---------------------------------------------

def uhj_to_b(a, encoder_kernels, mode = 'z', kind = 'fft', zi = None):
    """uhj_to_b(a, encoder_kernels, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - a                 : Input UHJ signal

        - decoder_kernels   : UHJ encoder kernels:
        
                        [[ left_W_FIR,  left_X_FIR,  left_Y_FIR],
                         [right_W_FIR, right_X_FIR, right_Y_FIR]]

                               shape = (2, UHJ_kernel_size, 3)

        - mode  : 'z' or 'full'. If mode is 'z', acts as a filter with state
                  'z', and returns a vector of length nframes(b). If mode is
                  'full', returns the full convolution.

        - kind  : 'direct' or 'fft', for direct or fft convolution

        - zi    : Initial state. An array of shape...

                  (2, len(uhj_kernel_size) - 1, 3)
                  
                  If zi = None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y         : The output of the encoder.
      zf        : If zi is None, this is not returned, otherwise, zf holds
                  the filter state.
    
    Encode a two dimensional UHJ stereo signal to horizontal
    ambisonic B-format using the supplied UHJ encoder kernels.

    """

    M = nframes(a)                      # length of input UHJ
    N = shape(encoder_kernels)[1]       # length of UHJ kernels
    m = 2                               # number of stereo channels

    # initialise result to correct size
    if mode is 'z':
        res = zeros((M, 3))
    else:
        res = zeros((M+N-1, 3))

    # convolve with UHJ kernel
    if zi is not None:
        zf = zeros_like(zi)
        
        for i in range(m):
            res_i, zf_i = convfilt(
                interleave(a[:, i]),
                encoder_kernels[i],
                mode,
                kind,
                zi[i]
                )

            res += res_i
            zf[i] = zf_i
    else:
        for i in range(m):
            res += convfilt(
                interleave(a[:, i]),
                encoder_kernels[i],
                mode,
                kind
                )

    res = hstack((res, zeros((nframes(res), 1))))

    if zi is not None:
        return res, zf
    else:
        return res


def superstereo(a, encoder_kernels, mode = 'z', kind = 'fft', zi = None):
    """superstereo(a, encoder_kernels, mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - a                 : Stereo signal

        - decoder_kernels   : Super Stereo encoder kernels:
        
                        [[ left_W_FIR,  left_X_FIR,  left_Y_FIR],
                         [right_W_FIR, right_X_FIR, right_Y_FIR]]

                               shape = (2, UHJ_kernel_size, 3)

        - mode  : 'z' or 'full'. If mode is 'z', acts as a filter with state
                  'z', and returns a vector of length nframes(b). If mode is
                  'full', returns the full convolution.

        - kind  : 'direct' or 'fft', for direct or fft convolution

        - zi    : Initial state. An array of shape...

                  (2, len(uhj_kernel_size) - 1, 3)
                  
                  If zi = None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y         : The output of the encoder.
      zf        : If zi is None, this is not returned, otherwise, zf holds
                  the filter state.
    
    Encode a two dimensional stereo signal to horizontal
    ambisonic B-format using the supplied Super Stereo encoder kernels.

    """
    
    # Wrapper for uhj_to_b()
    
    return uhj_to_b(a, encoder_kernels, mode, kind, zi)


def simplestereo(a, theta = 0.):
    """simplestereo(a, theta = 0.)
    
    Args:
        - a     : Stereo signal

        - theta : Stereo encoding angle, -pi/2 to pi/2. The default value,
                    0, place input stereo left/right at +/- pi/2 (hard
                    left/right). Theta = pi/2 collapses the image to front
                    centre, while -pi/2 collapses to back centre.
        
    Outputs: (y)
    
      y         : The output of the encoder.
    
    Encode a two dimensional stereo signal to horizontal
    ambisonic B-format using the Simple Stereo encoding method.

    """

    # coefficients
    k0 = 0.5 * ones_like(theta)
    k1 = C.rec_sqrt2 * sin(theta)
    k2 = C.rec_sqrt2 * cos(theta)

    # construct appropriate transform
    transform = array([
            [k0,  k0],
            [k1,  k1],
            [k2, -k2]
            ])

    if not(isscalar(theta)):                       # reshape for vectors
        transform = transform.transpose(2, 0, 1)

    res = (a[:, newaxis, :] * transform).sum(axis = -1)

    res = hstack((res, zeros((nframes(a), 1))))

    return res


#---------------------------------------------
# Special encoders
#
#   ZoomH2
#---------------------------------------------

def zoomH2_to_b(a):
    """zoomH2_to_b(a)
    
    Args:
        - a         : Input 4 channel signal from Zoom H2 recorder.

    Encode a Zoom H2 signal into the B-format domain.
    Channels should be in order: FL, FR, BL, BR

    """

    # construct appropriate encoder
    encoder = array(
        [[C.rec_sqrt2sum2, C.rec_sqrt2sum2, C.sqrt2div_sqrt2sum2, C.sqrt2div_sqrt2sum2],
         [C.twodiv_sqrt2sum1, C.twodiv_sqrt2sum1, -C.twodiv_sqrt2sum1, -C.twodiv_sqrt2sum1],
         [C.twodiv_sqrt2sum_sqrt3, -C.twodiv_sqrt2sum_sqrt3, C.twodiv_sqrt2sum_sqrt3, -C.twodiv_sqrt2sum_sqrt3],
         [0., 0., 0., 0.]]
        )

    # encode here!
    return inner(a, encoder)
