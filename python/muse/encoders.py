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


# it may be desirable to add "forward preference" to the uhj encoder
def uhj_to_b(a, hilbert_kernel, fpref = .354, \
             mode = 'z', kind = 'fft', zi = None):
    """uhj_to_b(a, hilbert_kernel, fpref = .354, \
                mode = 'z', kind = 'fft', zi = None)
    
    Args:
        - a                 : Input UHJ signal
        - hilbert_kernel    : Complex Hilbert transform kernel.
                              Real contains real resonse, Complex
                              contains complex response.
        - fpref : Forward preference, pushes phasiness to the rear
                  of the sound field where it is less offensive. With
                  no forward preference (= 0), the phasiness is distributed
                  equally around 360 degrees. .354 gives no phasiness at
                  front. 0. < fpref < .5
        - mode  : 'z' or 'full'. If mode is 'z', acts as a filter
                  with state 'z', and returns a vector of length
                  len(x). If mode is 'full', returns the full
                  convolution.
        - kind  : 'direct' or 'fft', for direct or fft convolution

        - zi    : Initial state. An array of shape (len(kernel) - 1, 2).
                  If zi=None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y         : The output of the encoder.
      zf        : If zi is None, this is not returned, otherwise, zf holds the
                  filter state.
    
    Decode a two dimensional ambisonic UHJ signal to four channel
    ambisonic b-format using a linear phase hilbert transform filter.
    (Channel Z is zeros.)

    """
    #     Audio Engineering Society E-Library
    #     Ambisonic Decoders for HDTV
    #     Preprint Number:   3345    Convention:   92 (February 1992)
    #     Authors:   Gerzon, Michael A.; Barton, Geoffrey J.
    #     E-library Location: (CD aes12)   /pp9193/pp9203/3405.pdf

    #     S = (Left + Right)/2.0
    #     D = (Left - Right)/2.0

    #     W = 0.982*S + j*0.164*D
    #     X = 0.419*S - j*0.828*D
    #     Y = 0.763*D + j*0.385*S

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


    # J Anderson on SuperStereo with Forward Preference:
    # 
    # The above two transforms can be combined
    # to give the following: 
    #
    #     W = 0.982*S + j*0.164*D
    #     X = 0.419*S - j*0.828*D
    #     Y = 0.763*D - j*a*0.763*D + a*0.385*S + j*0.385*S
    # 
    # a = .35355342513417343 for 'optimum', e.g., no phasiness at front

    s = a.sum(axis = -1)
    d = (array([1., -1.]) * a).sum(axis = -1)

    # convolve with hilbert kernel
    if zi is not None:
        hb_real, zf_real = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.real,
            mode,
            kind,
            zi.real
            )
        hb_imag, zf_imag = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.imag,
            mode,
            kind,
            zi.imag
            )
        zf = zf_real + 1j * zf_imag

    else:
        hb_real = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.real,
            mode,
            kind
            )
        hb_imag = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.imag,
            mode,
            kind
            )

    w = .982 * hb_real[:, 0] + .164 * hb_imag[:, 1]
    x = .419 * hb_real[:, 0] - .828 * hb_imag[:, 1]
    y = .763 * (hb_real[:, 1] - fpref * hb_imag[:, 1]) + \
        .385 * (fpref * hb_real[:, 0] + hb_imag[:, 0])
    z = zeros(nframes(hb_real))
    
    res = .5 * interleave(array([w, x, y, z]))

    if zi is not None:
        return res, zf
    else:
        return res


def superstereo(a, hilbert_kernel, width = .593, fpref = .354, mode = 'z', \
                    kind = 'fft', zi = None):
    """superstereo(a, hilbert_kernel, width = .593, fpref = .354, mode = 'z', \
                    kind = 'fft', zi = None)
    
    Args:
        - a                 : Input stereo signal (2 channel L/R)
        - hilbert_kernel    : Complex Hilbert transform kernel.
                              Real contains real resonse, Complex
                              contains complex response.
        - width     : Width setting, the optimal (in some senses)
                      value of which is about .593.
        - fpref     : Forward preference, pushes phasiness to the rear
                      of the sound field where it is less offensive. With
                      no forward preference (= 0), the phasiness is distributed
                      equally around 360 degrees. .354 gives no phasiness at
                      front. 0. < fpref < .5
        - mode      : 'z' or 'full'. If mode is 'z', acts as a filter
                      with state 'z', and returns a vector of length
                      len(x). If mode is 'full', returns the full
                      convolution.
        - kind      : 'direct' or 'fft', for direct or fft convolution

        - zi        : Initial state. An array of shape (len(kernel) - 1, 2).
                      If zi=None or is not given then initial rest is assumed.

    Outputs: (y, {zf})
    
      y             : The output of the encoder.
      zf            : If zi is None, this is not returned, otherwise, zf
                      holds the filter state.
    
    Encode a two channel L/R stereo signal via the ambisonic
    "super stereo" method to four channel ambisonic b-format
    using a linear phase hilbert transform filter.
    (Channel Z is zeros.)

    """
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


    # J Anderson on SuperStereo with Forward Preference:
    # 
    # Substituting -a allows the forward preference to act as expected (or
    # using -j for super stereo). The above two transforms can be combined
    # to give the following: 
    #
    # W = 0.6098637 * S - j * 0.6896511 * w * D
    # X = 0.8624776 * S + j * 0.7626955 * w * D
    # Y = (1.6822415 + a * 0.6896511) * w * D - j * (0.2156194 - a * 0.6098637) * S
    # 
    # a = .35355342513417343 for 'optimum', e.g., no phasiness at front


    s = a.sum(axis = -1)
    d = (array([1., -1.]) * a).sum(axis = -1)

    # convolve with hilbert kernel
    if zi is not None:
        hb_real, zf_real = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.real,
            mode,
            kind,
            zi.real
            )
        hb_imag, zf_imag = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.imag,
            mode,
            kind,
            zi.imag
            )
        zf = zf_real + 1j * zf_imag

    else:
        hb_real = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.real,
            mode,
            kind
            )
        hb_imag = convfilt(
            interleave(array([s, d])), # interleave s & d
            hilbert_kernel.imag,
            mode,
            kind
            )

    w = .6098637 * hb_real[:, 0] - .6896511 * width * hb_imag[:, 1]
    x = .8624776 * hb_real[:, 0] + .7626955 * width * hb_imag[:, 1]
    y = (1.6822415 + fpref * .6896511) * width * hb_real[:, 1] \
        - (.2156194 - fpref * .6098637) * hb_imag[:, 0]
    z = zeros(nframes(hb_real))
    
    res = .5 * interleave(array([w, x, y, z]))

    if zi is not None:
        return res, zf
    else:
        return res


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
    encoder = (array(decoder_dict[orientation]) * reciprocal(array(weight_dict[weight]))).transpose()

    # encode here!
    return inner(a, encoder)


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
