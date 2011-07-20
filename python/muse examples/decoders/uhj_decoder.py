# **************************************
# 
# UHJ decoder
# 
# decode a b-format signal to UHJ
# 
# the resulting decoded signal is generated 'all at once'
# using functional syntax
#
# **************************************


from muse import *


# params
sr = 44100          # sample rate
k = 1.              # velocity decode (shelf later)

gain        = -6           # in dB

file_type   = 'wav'         # write file...
encoding    = 'pcm24'
endianness  = 'file'

file_dir    =   '/Volumes/Audio/test/'      #write file dir

read_file = '/Volumes/Audio/SinS/sweet_1.wav'

# UHJ kernel length
N = 2048

file_name = 'UHJ_decoder_test'


# calculate
scale = db_to_amp(gain)

# ----- generate file name
write_file  =   file_dir + file_name + '.' + file_type[:3]


# ----- generate UHJ kernel
decoder_kernels = uhj_decoder_kernel(N)


# ************************************************************
# Set up sndfiles for reading:
read_sndfile = Sndfile(read_file)

#--------------------
# read test signal
b = read_sndfile.read_frames(read_sndfile.nframes)

#--------------------
# convolve

res = b_to_uhj(b, decoder_kernels, 'full')

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

