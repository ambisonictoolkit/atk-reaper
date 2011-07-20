# **************************************
# 
# UHJ decoder
# 
# decode a b-format signal to stereo (virtual mics)
# 
# the resulting decoded signal is generated 'all at once'
# using functional syntax
#
# **************************************


from muse import *


# params
sr = 44100          # sample rate
k = 1.              # velocity decode (shelf later)

gain        = -9            # in dB

angle       = 0.7854                # +-45 deg (Blumlein)
k           = 1.0                   # fig-8
##angle       = deg_to_rad(110/2)     # +-55 deg (Hyper-cardioid coincident)
##k           = 0.75                  # hyper-cardioid
##angle       = deg_to_rad(135./2)    # +-67.5 deg (Cardioid coincident)
##k           = 0.5                   # cardioid
##angle       = pi/3                  # +-60 deg (Uniform reverb)
##k           = 0.5                   # cardioid
##angle       = pi/2                  # +-90 deg (Sphere)
##k           = 0.5                   # cardioid

file_type   = 'wav'         # write file...
encoding    = 'pcm24'
endianness  = 'file'

file_dir    =   '/Volumes/Audio/test/'      #write file dir

read_file = '/Volumes/Audio/SinS/sweet_1.wav'

file_name = 'stereo_decoder_test'


# calculate
scale = db_to_amp(gain)

# ----- generate file name
write_file  =   file_dir + file_name + '.' + file_type[:3]


# ************************************************************
# Set up sndfiles for reading:
read_sndfile = Sndfile(read_file)

#--------------------
# read test signal
b = read_sndfile.read_frames(read_sndfile.nframes)

#--------------------
# decode

res = b_to_stereo(b, angle, k)

res *= scale

# ************************************************************
# Set up sndfiles for writing:
write_sndfile = Sndfile(
    write_file,
    'w',
    Format(file_type, encoding, endianness),
    nchannels(res),
    sr
    )

# ----- write out!
write_sndfile.write_frames(res)

# ----- close file
read_sndfile.close()
write_sndfile.close()

