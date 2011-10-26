# ****************************************************************************
# Isophonics (Centre for Digital Music, Queen Mary University, London)
#       RIR database) encoder kernels
#       RIR ---> greathall
# 
# Generate a set of kernels suitable for encoding a mono signal to b-format,
# using Isophonics RIR encoding.
#
# Kernels are classified by kernel size, N, subject, and stored in directory: 
#
#       '(root or ~)/Library/Application Support/ATK/kernels/FOA/encoders/ \
#           greathall/[sr]/None/[subject_id]'
#
# Within, one [mono] four channel [W, X, Y, Z] kernels are found, named:
#
#       'greathall'
#
# A resulting b-format encode can be generated from these kernels as follows:
#
#       b-format = 'greathall' * M
#
#   where * is convolution
#
# ****************************************************************************

from muse import *
import os

# params
srs         = array([96000])    # sample rates
Ns          = array([None])      # kernel lengths


# RIR dirs
rir_database_dir = '/Library/Application Support/RIR/Isophonics/greathall'


file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/encoders/greathall'

file_names  = ['greathall']


# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)
    rir_database_dir = os.path.expanduser('~' + rir_database_dir)

# generate subject IDs (retrieved from directories in rir_database_dir)
subject_ids = []
for file_name in os.listdir(os.path.join(rir_database_dir, 'Omni')):
    if os.path.isfile(os.path.join(rir_database_dir, 'Omni', file_name)):
        if file_name[0] != '.':
            subject_ids += [os.path.splitext(file_name)[0]]

### ----- loop
for sr in srs:                          # SR
    for N in Ns:                        # kernel sizes

        for subject_id in subject_ids:

            # ----- generate encoder kernel
            encoder_kernel = isophonicsRIR(subject_id, rir_database_dir)
            encoder_kernel = array([encoder_kernel])

            # ----- generate file name
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
                    nchannels(encoder_kernel[i]),
                    int(sr or 44100)                    # sr defaults to 44100
                    )                                   # if none is supplied

                # ----- write out!
                write_sndfile.write_frames(encoder_kernel[i])

                # ----- close file
                write_sndfile.close()

