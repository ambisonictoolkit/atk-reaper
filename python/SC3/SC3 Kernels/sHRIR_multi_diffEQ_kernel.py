# -*- coding: utf-8 -*-
# ****************************************************************************
# sHRIR (Duda model) decoder kernels
# 
# Generate a set of kernels suitable for decoding b-format
# to binaural, using the Sphereical HRIR decoder.
#
# Kernels are classified by kernel size, N, subject, and stored in directory: 
#
#       'ATK_kernels/FOA/decoders/Spherical_HRIR/SR_[sr]/N_[kernel_size]/ \
#           [subject_id]'
#
# Within, four [W,X,Y,Z] two channel [L,R] kernels are found, named:
#
#       'HRIR_W', 'HRIR_X', 'HRIR_Y', 'HRIR_Z'
#
# A resulting binaural decode can be generated from these kernels as follows:
#
#       binaural = 'HRIR_W' * W + 'HRIR_X' * X \
#                   + 'HRIR_Y' * Y + 'HRIR_Z' * Z
#
#   where * is convolution, + is sum.
#
# ****************************************************************************

from muse import *
import os


# ****************************************************************************
# 
# Multi-band hrtf decoder using the synthetic spherical HRIR
#
# Built using three decoders. These are:
#
#   1) Classic NFC + shelf filtered decoder
#   2) Linear phase decoder
#   3) HF decoder - cardioid pairs
#
# ****************************************************************************


# -----functions

# Linkwitz-Riley Band-pass filter
def link_ril_bp(a, order, W0, W1):

    # HP
    if W0 != 0. and not(isnan(W0)):
        for i in range(2):      #Linkwitz-Riley
            a = fiir_hp(a[::-1], order, W0)
    # LP
    if W1 != 1. and not(isnan(W1)):
        for i in range(2):      #Linkwitz-Riley
            a = fiir_lp(a[::-1], order, W1)

    return a


# cosine probe (for single channel, linear phase, odd length kernels only!)
def cos_probe(a, Wn):

    M = len(a)

    phasor = fphasor(Wn, 0, M)
    phasor -= phasor[-1]/2
    prob_sig = cos(phasor)  # probe signal

    # --> extract magnitude at desired frequency
    res = sum(prob_sig * a)

    return res


# stereo kernels (2D)
def stereo_kernels(M, angle, k):

    proto_kernels = zeros([4, M, 2]) # NOTE: full 3D

    # --> generate prototype impulse
    impulse = zeros(M)
    impulse[M/2] = 1.0

    # --> generate decoder, cardioid pairs
    for harm in range(4):
        proto_kernels[harm] = b_to_stereo(
            roll(
                interleave(
                    array([
                        impulse,
                        zeros(M),
                        zeros(M),
                        zeros(M)
                        ])
                    ),
                harm,
                1
                ),
            angle,
            k
            )

    # remove zeros for Z, as the resulting kernels shouldn't include Z

    proto_kernels = transpose(proto_kernels)
    proto_kernels = array([
        proto_kernels[0, :, :-1],
        proto_kernels[1, :, :-1],
        ])
    res = transpose(proto_kernels)


    return res


# classic decoder using NFC and shelf filter -- W is linearised
def classic_kernels(M, pos, k, mod_r, psycho_freq, sr):

    proto_kernels = sHRIR_decoder_kernel(pos, k, (M+1)/2, 1./sr, \
                                             mod_r, \
                                             phase = True)
    num_harms = shape(proto_kernels)[0]
    num_chans = shape(proto_kernels)[2]

    res = zeros([num_harms, M, num_chans])

    
    # ---> shelf filter
    proto_kernels = transpose(proto_kernels)

    # add zeros for Z, for shelf filtering, NFC,
    # as the 2D generated kernels don't include Z
    proto_kernels = array([
        hstack((
            proto_kernels[0],
            zeros((shape(proto_kernels)[1],1))
            )),
        hstack((
            proto_kernels[1],
            zeros((shape(proto_kernels)[1],1))
            ))
        ])

    
    for chan in range(num_chans):
        proto_kernels[chan] = psycho_shelf(
            proto_kernels[chan],
            freq_to_Wn(psycho_freq, 1./sr),
            C.k_2D                  # 2D shelf
            )

    # remove zeros for Z, as the resulting kernels shouldn't include Z
    proto_kernels = array([
        proto_kernels[0, :, :-1],
        proto_kernels[1, :, :-1],
        ])

    proto_kernels = transpose(proto_kernels)


    # ---> linearise phase of W for decoder_1_kernel
    #      (keep X, Y w/ same phase relation)
    w_ap = zeros([(M+1)/2, num_chans])  # W linearising kernel,
                                        # for left and right ear W

    for w, chan in zip(deinterleave(proto_kernels[0]), \
                       range(num_chans)):
        w_ap[:, chan] = apf(w[::-1])

    # --> linearise phase here (size expanded to M)
    for harm in range(num_harms):
        res[harm] = convfilt(
            proto_kernels[harm],
            w_ap,
            'full'
            )

    return res


# linear decoder
def linear_kernels(M, pos, k, mod_r, sr):

    res = sHRIR_decoder_kernel(pos, k, M, 1./sr, \
                                             mod_r, \
                                             phase = False)

    return res



# Diffuse field normaliser
def normalise_diffuse_field(kern, norm_freqs, norm_coeffs):

    # NOTE: Kirkeby expects an odd length kernel as input

    # inits
    num_harms = shape(kern)[0]
    num_chans = shape(kern)[2]
    N = (shape(kern)[1] + 1) / 2
    

    # --> calculate normalisation kernel
    norm_kern = zeros([2*N + 1, num_chans])     # normalising kernel,
                                                # for left and right ear

    # find diffuse field magnitude (kernel size = 2N)
    diff_mag = zeros([N + 1, num_chans])
    diff_pha = lin([0, -pi/2 * (2*N-1)], N + 1)

    for harm in range(num_harms):
        for chan in range(num_chans):
            diff_mag[:, chan] += abs(
                rfft(
                    norm_coeffs[harm] * kern[harm][:, chan],
                    2*N
                    )
                )**2 / num_harms
    diff_mag = sqrt(diff_mag)

    # construct diffuse field kernel
    diff_kern = zeros([2*N, num_chans])
    for chan in range(num_chans):
        diff_kern[:, chan] = irfft(
            diff_mag[:, chan] * (cos(diff_pha) + sin(diff_pha) * 1j)
            )

    # construct normalisation kernel (kernel size = 2N + 1)
    for d, chan in zip(deinterleave(diff_kern), range(num_chans)):
        norm_kern[:, chan] = linf(
            kirkeby(
                2 * append(d, 0),          # to match UHJ normal
                freq_to_Wn(norm_freqs, 1./sr)
                )
            )

    # --> normalise here (and trim, kernel size = 2N -1)
    for harm in range(num_harms):
        kern[harm] = convfilt(
                        kern[harm],
                        norm_kern,
                        'full'
                        )[N:N + (2*N - 1)]

    return kern


# -----parameters
##srs         = array([44100, 48000, 88200, 96000, 192000]) # sample rates
##Ns          = array([256, 512])    # kernel lengths

# NOTE: for distribution Ns should be correlated to SR
#       44.1kHz: 512; 88.2kHz: 1024; 176.4kHz: 2048
#
# -->[sr, N]
#    [44100, 512]
#    [48000, 512]
#    [88200, 1024]
#    [96000, 1024]
#    [192000, 2048]

##srs         = array([44100]) # sample rates
##Ns          = array([512])    # kernel lengths

##srs         = array([48000]) # sample rates
##Ns          = array([512])    # kernel lengths

##srs         = array([88200]) # sample rates
##Ns          = array([1024])    # kernel lengths

##srs         = array([96000]) # sample rates
##Ns          = array([1024])    # kernel lengths

srs         = array([192000]) # sample rates
Ns          = array([2048])    # kernel lengths



# --> Duda data
a           = 2./3              # Duda anthro model values
var         = 0.123494002733    # used to generate radius map
med_r       = 0.0864243630344   # to generate various subjects

num_subjects = 10               # number of subjects to generate


# --> Decoders...
gain = 3.0 # final gain... similar to UHJ
k = 1.0     # velocity for all decoders
k2 = 0.5    # cardioid for HF

angle = 5./9*pi         # 1/2 angle for HF cardioids (HF decoder)
num_pos = 360           # number of positions to sample

# See: J. Daniel, "Evolving Views on HOA: from Technological
#      to Pragmatic Concerns," presented at the Ambisonics Symposium 2009,
#      Graz, 2009.
#
#      J. Daniel, "Représentation de champs acoustiques, application
#      à la transmission et à la reproduction de scènes sonores complexes
#      dans un contexte multimédia," PhD Thesis, Université Paris 6,
#      Paris, 2001., pg. 169, pg. 171
#
#                   1       2       3       4       5       Order
#                   1       2       3       4       5       Order
# free (cylinder)   714     1297    1899    2509    3125    Hz
# hrtf (cylinder)   543     1042    1575    2127    2691
psycho_freq = 700.
band_freqs = [NaN, 2000., 7000., NaN];   # filter bands

order = 2               # crossover filter order, x2 (Linkwitz-Riley)


norm_coeffs = sqrt([2, 2, 2]) # normalisation for 2D
norm_freqs = [NaN, NaN] # freqs to normalise between


num_harms = 3 # 2D
num_chans = 2 # stereo output (binaural)

width = pi  # window width for final windowing

M = 1023 # working kernel size - must be odd: 2**x -1


# --> file format and output
file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/decoders/spherical'

file_names  = ['HRIR_W', 'HRIR_X', 'HRIR_Y']

# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)


# -----initial calcs

# --> set up positions for decoder
pos = pol_to_cart(                  # horizontal loudspeaker layout
    interleave(                     #   (M+1) points
        array([ones((M+1)/2),
               arange(0, pi, pi/((M+1)/2))]))
    )

# --> generate subjects and subject data
mod_r = (var * (                        # model to match duda's measured rs
    (1-a) * lin([-1, 1], num_subjects) + \
    a * (lin([-1, 1], num_subjects))**3) + 1) \
    * med_r

subject_ids = []                        # subject ids
for i in range(num_subjects):
    subject_ids += [str(i).zfill(4)]



# ----- loop
for sr in srs:                          # SR
    for N in Ns:                        # kernel sizes

        for subject in arange(num_subjects):  # subject index


            # ----- 1) generate decoder kernels (individual for each subject)
            #          build low-mid freq decoder kernels
            #          velocity (strict soundfield) with shelf            
            decoder_1_kernels = classic_kernels(M, pos, k, mod_r[subject],
                                psycho_freq, sr)


            # ----- 2) generate decoder kernels (individual for each subject)
            #          build mid-high freq decoder kernels
            #          linear phase velocity (strict soundfield)            
            decoder_2_kernels = linear_kernels(M, pos, k, mod_r[subject],
                                sr)


            # ----- 3) generate decoder kernels (individual for each subject)
            #          build high-freq decoder kernels
            #          linear phase cardioid pairs            
            decoder_3_kernels = stereo_kernels(M, angle, k2)


            # ----- 5) crossover, sum and normalise

            # --> analyse W (l+r) at crossover frequencies
            scale_1 = 1.0
            scale_2 = scale_1 * cos_probe(
                sum(decoder_1_kernels[0], 1)/2,
                freq_to_Wn(band_freqs[1], 1./sr)
                ) / \
                cos_probe(
                sum(decoder_2_kernels[0], 1)/2,
                freq_to_Wn(band_freqs[1], 1./sr)
                )
            scale_3 = scale_2 * cos_probe(
                sum(decoder_2_kernels[0], 1)/2,
                freq_to_Wn(band_freqs[2], 1./sr)
                ) / \
                cos_probe(
                sum(decoder_3_kernels[0], 1)/2,
                freq_to_Wn(band_freqs[2], 1./sr)
                )

            # --> equalise gain at crossover frequencies
            decoder_1_kernels *= scale_1
            decoder_2_kernels *= scale_2
            decoder_3_kernels *= scale_3

            # --> bandpass filter
            for harm in range(num_harms):
                decoder_1_kernels[harm] = link_ril_bp(
                    decoder_1_kernels[harm],
                    order,
                    freq_to_Wn(band_freqs[0], 1./sr),
                    freq_to_Wn(band_freqs[1], 1./sr)
                    )
                decoder_2_kernels[harm] = link_ril_bp(
                    decoder_2_kernels[harm],
                    order,
                    freq_to_Wn(band_freqs[1], 1./sr),
                    freq_to_Wn(band_freqs[2], 1./sr)
                    )
                decoder_3_kernels[harm] = link_ril_bp(
                    decoder_3_kernels[harm],
                    order,
                    freq_to_Wn(band_freqs[2], 1./sr),
                    freq_to_Wn(band_freqs[3], 1./sr)
                    )

            # --> sum and normalise
            proto_kernels = decoder_1_kernels + decoder_2_kernels + \
                              decoder_3_kernels

            proto_kernels = normalise_diffuse_field(
                proto_kernels,
                norm_freqs,
                norm_coeffs
                )


            # ----- 6) centre, trim and window
            decoder_kernels = zeros([num_harms, N, num_chans])

            for harm in range(num_harms):
                decoder_kernels[harm] = convfilt(
                    proto_kernels[harm],
                    sinc(lin([-(M)/2., (M)/2.], M+1)) * hann(M+1),
                    'full'
                    )[(2*M-N)/2:(2*M-N)/2 + N] * \
                    interleave(kaiser(N, width)) * db_to_amp(gain)


            # ----- generate file names
            write_files = [] 
            for name in file_names:
                write_files += [
                    os.path.join(library_dir + database_dir, str(sr), str(N), \
                                 subject_ids[subject], \
                                 name + '.' + file_type[:3])
                    ]

            # ----- write out decoder kernels
            for i in range(len(write_files)):

                # ************************************************************
                # Set up sndfiles for writing:

                if not os.path.exists(os.path.dirname(write_files[i])):
                    os.makedirs(os.path.dirname(write_files[i]))
                
                write_sndfile = Sndfile(
                    write_files[i],
                    'w',
                    Format(file_type, encoding, endianness),
                    nchannels(decoder_kernels[i]),
                    int(sr or 44100)                    # sr defaults to 44100
                    )                                   # if none is supplied

                # ----- write out!
                write_sndfile.write_frames(decoder_kernels[i])

                # ----- close file
                write_sndfile.close()

