# ****************************************************************************
# UHJ encoder kernels
# 
# Generate a set of kernels suitable for encoding a stereo UHJ signal
# to b-format.
#
# Kernels are classified by kernel size, N, and stored in directory: 
#
#       'ATK_kernels/FOA/encoders/UHJ/SR_00None/N_[kernel_size]/'
#
# Within, two [L, R] three channel [W, X, Y] kernels are found, named:
#
#       'UHJ_L', 'UHJ_R'
#
# A resulting b-format encode can be generated from these kernels as follows:
#
#       b-format = 'UHJ_L' * L + 'UHJ_R' * R
#
#   where * is convolution, + is sum.
#
# ****************************************************************************

from muse import *
import os

# params
srs         = array([44100, 48000, 88200, 96000, 192000]) # sample rates
Ns          = array([512, 1024, 2048, 4096, 8192])    # kernel lengths

psycho_freq = 400.          # freq (Hz) for psychoacoutic shelf filtering
                            # this is required for UHJ encoding!


file_type   = 'wav'         # write file...
#encoding    = 'pcm24'
encoding    = 'pcm32'
endianness  = 'file'


target_dir  = '/Volumes/Audio/test'      #temp write dir
file_dir    = '/ATK_kernels/FOA/encoders/UHJ'

file_names  = ['UHJ_L', 'UHJ_R']

subject_ids = ['0000']                  # only one subject
                                        # could add more, for 'preference'

# ----- loop
for sr in srs:                          # SR
    for N in Ns:                        # kernel sizes

        for subject_id in subject_ids:

            # ----- generate encoder kernels
            encoder_kernels = uhj_encoder_kernel(
                N, freq_to_Wn(psycho_freq, 1./sr)
                )

            # ----- generate file names
            write_files = []
            for name in file_names:
                write_files += [
                    target_dir + file_dir + \
                    '/SR_' + str(sr).zfill(6) + '/N_' + str(N).zfill(4) + '/' + \
                    subject_id + '/' + name + '.' + file_type[:3]
                    ]

            # ----- write out encoder kernels
            for i in range(len(write_files)):

                # ************************************************************
                # Set up sndfiles for writing:

                if not os.path.exists(os.path.dirname(write_files[i])):
                    os.makedirs(os.path.dirname(write_files[i]))
                
                write_sndfile = Sndfile(
                    write_files[i],
                    'w',
                    Format(file_type, encoding, endianness),
                    nchannels(encoder_kernels[i]),
                    int(sr or 44100)                    # sr defaults to 44100
                    )                                   # if none is supplied

                # ----- write out!
                write_sndfile.write_frames(encoder_kernels[i])

                # ----- close file
                write_sndfile.close()

