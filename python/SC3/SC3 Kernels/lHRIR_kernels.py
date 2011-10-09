# ****************************************************************************
# lHRIR (Listen HRIR database) decoder kernels
# 
# Generate a set of kernels suitable for decoding b-format
# to binaural, using the Listen HRIR decoder.
#
# Kernels are classified by kernel size, N, subject, and stored in directory: 
#
#       '(root or ~)/Library/Application Support/ATK/kernels/FOA/decoders/ \
#           listen/[sr]/[kernel_size]/[subject_id]'
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

# params
srs         = array([44100])    # sample rates
Ns          = array([512])      # kernel lengths


k = 1.                      # velocity decode (shelf later)
psycho_freq = 700.          # freq (Hz) for psychoacoutic shelf filtering
nfc_r = 1.95                # NFC distance (for compensation)

pos = spher_to_cart(                    # B3
    array([[1, deg_to_rad(0), deg_to_rad(-30)],
           [1, deg_to_rad(60), deg_to_rad(30)],
           [1, deg_to_rad(-60), deg_to_rad(30)]])
    )

# HRIR dirs
hrir_database_dir = '/Library/Application Support/HRTF/Listen/hrir_database'


file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/decoders/listen'

file_names  = ['HRIR_W', 'HRIR_X', 'HRIR_Y', 'HRIR_Z']


# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)
    hrir_database_dir = os.path.expanduser('~' + hrir_database_dir)

# generate subject IDs (retrieved from directories in hrir_database_dir)
subject_ids = []
for dir_name in os.listdir(hrir_database_dir):
    if os.path.isdir(os.path.join(hrir_database_dir, dir_name)):
        subject_ids += [dir_name[4:]]


# ----- loop
for sr in srs:                          # SR
    for N in Ns:                        # kernel sizes

        # generate reference to normalise spherical against
        sphere_rms = rms(sHRIR_decoder_kernel(pos, k, N, 1./sr)[0])

        for subject_id in subject_ids:

            # ----- generate decoder kernels
            decoder_kernels = lHRIR_decoder_kernel(pos, k, \
                                                   subject_id, \
                                                   hrir_database_dir)
            # ----- normalise against spherical HRIR
            decoder_kernels *= sphere_rms / rms(decoder_kernels[0])

            # ----- shelf filter and NFC decoder kernels
            b_format_kernels = transpose(decoder_kernels)

            for i in range(2):
                b_format_kernels[i] = psycho_shelf(
                    nfc(b_format_kernels[i], nfc_r, 1./sr),
                    freq_to_Wn(psycho_freq, 1./sr),
                    C.k_3D                  # 3D shelf
                    )

            decoder_kernels = transpose(b_format_kernels)


            # ----- generate file names
            write_files = [] 
            for name in file_names:
                write_files += [
                    os.path.join(library_dir + database_dir, str(sr), str(N), \
                                 subject_id, \
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

