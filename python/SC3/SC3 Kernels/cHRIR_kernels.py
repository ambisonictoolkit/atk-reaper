# ****************************************************************************
# cHRIR (CIPIC HRIR database) decoder kernels
# 
# Generate a set of kernels suitable for decoding b-format
# to binaural, using the CIPIC HRIR decoder.
#
# Kernels are classified by kernel size, N, subject, and stored in directory: 
#

#       'ATK_kernels/FOA/decoders/CIPIC_HRIR/SR_[sr]/N_[kernel_size]/ \
#           [subject_ID]'

#       '(root or ~)/Library/Application Support/ATK/kernels/FOA/decoders/ \
#           cipic_raw/[sr]/[kernel_size]/[subject_id]'
#

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
srs         = array([44100]) # sample rates
Ns          = array([256])    # kernel lengths


k = 1.                      # velocity decode (shelf later)
psycho_freq = 700.          # freq (Hz) for psychoacoutic shelf filtering
nfc_r = 1.00                # NFC distance (for compensation)

#                                       B3c
#   [[  0,     -33.750],
#    [  57.609, 33.134],
#    [ -57.609, 33.134],
az = rad_to_deg(arctan2(
    -sin(deg_to_rad(-45)),
    cos(deg_to_rad(-45)) * cos(deg_to_rad(50.625))
    ))
el = rad_to_deg(arcsin(cos(deg_to_rad(-45)) * sin(deg_to_rad(50.625))))
pos = spher_to_cart(
    array([[1, deg_to_rad(0), deg_to_rad(-33.75)],
           [1, deg_to_rad(az), deg_to_rad(el)],
           [1, deg_to_rad(-az), deg_to_rad(el)]])
    )

# HRIR dirs
#database_dir = '/Users/josephla/Documents/Developer/CIPIC_hrtf_database/standard_hrir_database/'
hrir_database_dir = '/Library/Application Support/HRTF/CIPIC/standard_hrir_database'


file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


##target_dir  = '/Volumes/Audio/test'      # temp write dir
##file_dir    = '/ATK_kernels/FOA/decoders/CIPIC_HRIR'
user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/decoders/cipic_raw'

file_names  = ['HRIR_W', 'HRIR_X', 'HRIR_Y', 'HRIR_Z']


# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)
    hrir_database_dir = os.path.expanduser('~' + hrir_database_dir)

# generate subject IDs
subject_ids = []
for dir_name in os.listdir(hrir_database_dir):
    if os.path.isdir(os.path.join(hrir_database_dir, dir_name)):
        name = dir_name[8:]
        if name.isdigit():
            subject_ids += [name]


# ----- loop
for sr in srs:                          # SR
    for N in Ns:                        # kernel sizes

        # generate reference to normalise spherical against
        sphere_rms = rms(sHRIR_decoder_kernel(pos, k, N, 1./sr)[0])

        for subject_id in subject_ids:

            # ----- generate decoder kernels
            decoder_kernels = cHRIR_decoder_kernel(pos, k, \
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
##                write_files += [
##                    target_dir + file_dir + \
##                    '/SR_' + str(sr).zfill(6) + '/N_' + str(N).zfill(4) + '/' + \
##                    '0' + subject_id + '/' + name + '.' + file_type[:3]
##                    ]
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

