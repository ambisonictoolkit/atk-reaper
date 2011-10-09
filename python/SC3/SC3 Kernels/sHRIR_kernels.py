# ****************************************************************************
# sHRIR (Duda model) decoder kernels
# 
# Generate a set of kernels suitable for decoding b-format
# to binaural, using the Sphereical HRIR decoder.
#
# Kernels are classified by kernel size, N, and stored in directory: 
#
#       '(root or ~)/Library/Application Support/ATK/kernels/FOA/decoders/ \
#           spherical/[sr]/[kernel_size]/[subject_id]'
#
# Within, three [W,X,Y] two channel [L,R] kernels are found, named:
#
#       'HRIR_W', 'HRIR_X', 'HRIR_Y', 'HRIR_Z'
#
#   NOTE: The spherical model is symmetric about Z, so no Z kernel
#           is generated.
#
# A resulting binaural decode can be generated from these kernels as follows:
#
#       binaural = 'HRIR_W' * W + 'HRIR_X' * X + 'HRIR_Y' * Y
#
#   where * is convolution, + is sum.
#
# Subjects with varying head radii are generated using a model
# of the data found within anthro.mat. See 'chrir_anthro.py for
# details.
#
# ****************************************************************************

from muse import *
import os

# params
srs         = array([44100, 48000, 88200, 96000, 192000])   # sample rates
Ns          = array([256, 512, 1024, 2048, 4096, 8192])     # kernel lengths

a           = 2./3              # Duda anthro model values
var         = 0.123494002733    # used to generate radius map
med_r       = 0.0864243630344   # to generate various subjects

num_subjects = 10               # number of subjects to generate


k = 1.                      # velocity decode (shelf later)
psycho_freq = 700.          # freq (Hz) for psychoacoutic shelf filtering

pos = pol_to_cart(                    # horzontal loudspeaker layout
    array([[1, deg_to_rad(0)],        #     a cross
           [1, deg_to_rad(90)]])
    )

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


# generate subjects and subject data
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

            # ----- generate decoder kernels (individual for each subject)
            decoder_kernels = sHRIR_decoder_kernel(pos, k, N, 1./sr, \
                                                   mod_r[subject])

            # ----- shelf filter and NFC decoder kernels
            b_format_kernels = transpose(decoder_kernels)

            # add zeros for Z, for shelf filtering, NFC,
            # as the 2D generated kernels don't include Z
            b_format_kernels = array([
                hstack((
                    b_format_kernels[0],
                    zeros((shape(b_format_kernels)[1],1))
                    )),
                hstack((
                    b_format_kernels[1],
                    zeros((shape(b_format_kernels)[1],1))
                    ))
                ])

            # psychoacoustic shelf and NFC filtering
            for i in range(2):
                b_format_kernels[i] = psycho_shelf(
                    b_format_kernels[i],
                    freq_to_Wn(psycho_freq, 1./sr),
                    C.k_2D                  # 2D shelf
                    )

            # remove zeros for Z, as the resulting kernels shouldn't include Z
            b_format_kernels = array([
                b_format_kernels[0, :, :-1],
                b_format_kernels[1, :, :-1],
                ])

            decoder_kernels = transpose(b_format_kernels)


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

