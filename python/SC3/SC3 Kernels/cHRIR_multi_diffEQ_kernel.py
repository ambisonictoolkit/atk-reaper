# -*- coding: utf-8 -*-
# ****************************************************************************
# cHRIR (CIPIC HRIR database) decoder kernels
# 
# Generate a set of kernels suitable for decoding b-format
# to binaural, using the Listen HRIR decoder.
#
# Kernels are classified by kernel size, N, subject, and stored in directory: 
#
#       'ATK_kernels/FOA/decoders/Listen_HRIR/SR_[sr]/N_[kernel_size]/ \
#           [subject_ID]'
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
# Multi-band hrtf decoder using the Listen database
#
# Built using four decoders. These are:
#
#   1) VLF - an omni (W only), to correct frequency response down to 0Hz
#   2) Classic NFC + shelf filtered decoder
#   3) Linear phase decoder
#   4) VHF decoder - extrapolated from decoder #3
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


# W only kernels (use for very low frequency)
def w_only_kernels(M, num_harms, num_chans):

    res = zeros([num_harms, M, num_chans])

    # --> generate decoder, just W
    for chan in range(num_chans):
        res[0][M/2, chan] = 1.0

    return res


# coefficient kernels (use for very high frequency)
def coeff_kernels(M, coeffs):

    num_harms = shape(coeffs)[0]
    num_chans = shape(coeffs)[1]

    res = zeros([num_harms, M, num_chans])

    # --> generate decoder, just W
    for harm in range(num_harms):
        for chan in range(num_chans):
            res[harm][M/2, chan] = coeffs[harm, chan]

    return res


# classic decoder using NFC and shelf filter -- W is linearised
def classic_kernels(M, pos, k, subject_id, hrir_database_dir,
                    nfc_r, psycho_freq, sr):

    proto_kernels = cHRIR_decoder_kernel(
        pos,
        k,
        subject_id,
        hrir_database_dir,
        phase = True
        )
    num_harms = shape(proto_kernels)[0]
    num_chans = shape(proto_kernels)[2]

    res = zeros([num_harms, M, num_chans])

    # --> pad with zeros, new size = (M+1)/2
    tmp = zeros([num_harms, (M+1)/2, num_chans])
    for harm in range(num_harms):
        tmp[harm] = over_dub(
            tmp[harm],
            proto_kernels[harm],
            ((M+1)/2 - shape(proto_kernels)[1]) / 2
            )
    proto_kernels = tmp
    
    # ---> shelf filter and near-field compensate
    proto_kernels = transpose(proto_kernels)
    for chan in range(num_chans):
        proto_kernels[chan] = psycho_shelf(
            nfc(proto_kernels[chan], nfc_r, 1./sr),
                freq_to_Wn(psycho_freq, 1./sr),
                C.k_3D                  # 3D shelf
                )
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
def linear_kernels(M, pos, k, subject_id, hrir_database_dir):

    proto_kernels = cHRIR_decoder_kernel(
        pos,
        k,
        subject_id,
        hrir_database_dir,
        phase = False
        )
    num_harms = shape(proto_kernels)[0]
    num_chans = shape(proto_kernels)[2]
    N = shape(proto_kernels)[1]

    res = zeros([num_harms, M, num_chans])

    # --> expanded size to M
    L = M - N + 1
    for harm in range(num_harms):
        res[harm] = convfilt(
            proto_kernels[harm],
            sinc(lin([-(L-1)/2., (L-1)/2.], L)) * hann(L),
            'full'
            )

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
srs         = array([44100])    # sample rates
##Ns          = array([256]) # target kernel lengths
Ns          = array([256, 512]) # target kernel lengths
##Ns          = array([512]) # target kernel lengths

# -------> edit target kernel lengths --> set = to 512 only
#          need to modify SC code to do so!!!!


# --> Decoders...
gain = 3.0 # final gain... similar to UHJ
k = 1.0 # velocity for all decoders

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
# free (sphere)     690     1250    1831    2423    3022    Hz
# hrtf (sphere)     528     1011    1532    2072    2626
psycho_freq = 700.
band_freqs = [NaN, 35., 2000., 18000., NaN];   # filter bands

order = 2               # crossover filter order, x2 (Linkwitz-Riley)

nfc_r = 1.95                # NFC distance (for compensation)


norm_coeffs = sqrt([2, 3, 3, 3]) # normalisation for 3D
norm_freqs = [NaN, NaN] # freqs to normalise between


num_harms = 4 # 3D
num_chans = 2 # stereo output (binaural)

width = pi  # window width for final windowing

# HRIR dirs
hrir_database_dir = '/Library/Application Support/HRTF/CIPIC/standard_hrir_database'
M = 4095 # working kernel size - must be odd: 2**x -1

# --> file format and output
file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/decoders/cipic'

file_names  = ['HRIR_W', 'HRIR_X', 'HRIR_Y', 'HRIR_Z']


# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)
    hrir_database_dir = os.path.expanduser('~' + hrir_database_dir)


# -----initial calcs

# --> set up positions for decoder

# constants: cHRIR azimuth and elevation (in degrees)
c_azimuths = deg_to_rad(concatenate((linspace(0, 45, 10), array([55, 65, 80]))))
c_azimuths = concatenate((-c_azimuths[1:][::-1], c_azimuths))
c_elevations = deg_to_rad(-45 + (360. / 64 * arange(17)))

pos_h = zeros([len(c_azimuths), 3])
pos = zeros([len(c_azimuths) * len(c_elevations), 3])

for i in range(len(c_azimuths)):
    pos_h[i] = spher_to_cart(array([1, c_azimuths[i], 0.]))

for i in range(len(c_elevations)):
    for j in range(len(c_azimuths)):
        pos[i + j * len(c_elevations)] = array([
            (cos(c_elevations[i]) * pos_h[j][0]) - (sin(c_elevations[i]) * pos_h[j][2]),
            pos_h[j][1],
            (sin(c_elevations[i]) * pos_h[j][0]) + (cos(c_elevations[i]) * pos_h[j][2])
            ])

# --> generate subject IDs
subject_ids = []
for dir_name in os.listdir(hrir_database_dir):
    if os.path.isdir(os.path.join(hrir_database_dir, dir_name)):
        name = dir_name[8:]
        if name.isdigit():
            subject_ids += [name]


# ----- loop
for sr in srs:                          # SR
    for N in Ns:                        # kernel sizes

        for subject_id in subject_ids:


            # ----- 1) build very low-freq decoder kernels
            #          linear phase omni pairs            
            decoder_1_kernels = w_only_kernels(M, num_harms, num_chans)


            # ----- 2) generate decoder kernels (individual for each subject)
            #          build low-mid freq decoder kernels
            #          velocity (strict soundfield) with shelf and NFC            
            decoder_2_kernels = classic_kernels(
                M, pos, k, subject_id, hrir_database_dir,
                    nfc_r, psycho_freq, sr
                )

            # ----- 3) generate decoder kernels (individual for each subject)
            #          build mid-high freq decoder kernels
            #          linear phase velocity (strict soundfield)            
            decoder_3_kernels = linear_kernels(
                M, pos, k, subject_id, hrir_database_dir
                )


            # ----- 4) generate decoder kernels (individual for each subject)
            #          build high-freq decoder kernels
            #          extract responses at crossover point from decoder_3            
            hfreq_coeffs = zeros([num_harms, num_chans])

            for harm in range(num_harms):
                for chan in range(num_chans):
                    hfreq_coeffs[harm, chan] = cos_probe(
                        decoder_3_kernels[harm][:, chan],
                        freq_to_Wn(band_freqs[3], 1./sr)
                        )

            decoder_4_kernels = coeff_kernels(M, hfreq_coeffs)


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
            scale_4 = scale_3 * cos_probe(
                sum(decoder_3_kernels[0], 1)/2,
                freq_to_Wn(band_freqs[3], 1./sr)
                ) / \
                cos_probe(
                sum(decoder_4_kernels[0], 1)/2,
                freq_to_Wn(band_freqs[3], 1./sr)
                )

            # --> equalise gain at crossover frequencies
            decoder_1_kernels *= scale_1
            decoder_2_kernels *= scale_2
            decoder_3_kernels *= scale_3
            decoder_4_kernels *= scale_4

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
                decoder_4_kernels[harm] = link_ril_bp(
                    decoder_4_kernels[harm],
                    order,
                    freq_to_Wn(band_freqs[3], 1./sr),
                    freq_to_Wn(band_freqs[4], 1./sr)
                    )

            # --> sum and normalise
            proto_kernels = decoder_1_kernels + decoder_2_kernels + \
                              decoder_3_kernels + decoder_4_kernels

            proto_kernels = normalise_diffuse_field(
                proto_kernels,
                norm_freqs,
                norm_coeffs
                )


            # ----- 6) centre, trim, window and gain
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
                                 '0' + subject_id, \
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
