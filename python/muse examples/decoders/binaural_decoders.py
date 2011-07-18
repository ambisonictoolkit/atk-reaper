# **************************************
# 
# binaural decoders
# 
# decode a b-format signal using
# either the Listen HRIR decoder
# or the synthetic Spherical HRIR decoder
# 
# the resulting decoded signal is generated 'all at once'
# using functional syntax
# 
# **************************************

from muse import *
import pylab as pl


# params
sr = 44100          # sample rate
k = 1.              # velocity decode (shelf later)


# speaker positions (pairs): un-comment your choice
# NOTE: prefer 'hex' arrangements (3D), especially sD & sH


# ------> 2D

##pos = pol_to_cart(                    # 1st choice
##    array([[1, deg_to_rad(0)],
##           [1, deg_to_rad(90)]])
##    )
##pos = pol_to_cart(                    # 2nd choice
##    array([[1, deg_to_rad(45)],
##           [1, deg_to_rad(-45)]])
##    )


# ------> 3D

# --- Octrahedron arrays
##pos = spher_to_cart(                    # sAO
##    array([[1, deg_to_rad(0), deg_to_rad(45)],
##           [1, deg_to_rad(0), deg_to_rad(-45)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # sBO
##    array([[1, deg_to_rad(0), deg_to_rad(0)],
##           [1, deg_to_rad(90), deg_to_rad(45)],
##           [1, deg_to_rad(90), deg_to_rad(-45)]])
##    )

# --- 'Regular' arrays
##pos = spher_to_cart(                    # sA
##    array([[1, deg_to_rad(0), deg_to_rad(30)],
##           [1, deg_to_rad(0), deg_to_rad(-30)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # sB
##    array([[1, deg_to_rad(0), deg_to_rad(0)],
##           [1, deg_to_rad(90), deg_to_rad(30)],
##           [1, deg_to_rad(90), deg_to_rad(-30)]])
##    )
##pos = spher_to_cart(                    # sC
##    array([[1, deg_to_rad(0), deg_to_rad(30)],
##           [1, deg_to_rad(60), deg_to_rad(-30)],
##           [1, deg_to_rad(-60), deg_to_rad(-30)]])
##    )
##pos = spher_to_cart(                    # sD
##    array([[1, deg_to_rad(0), deg_to_rad(-30)],
##           [1, deg_to_rad(60), deg_to_rad(30)],
##           [1, deg_to_rad(-60), deg_to_rad(30)]])
##    )

# --- vertically compressed arrays
##pos = spher_to_cart(                    # sE
##    array([[1, deg_to_rad(0), deg_to_rad(15)],
##           [1, deg_to_rad(0), deg_to_rad(-15)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # sF
##    array([[1, deg_to_rad(0), deg_to_rad(0)],
##           [1, deg_to_rad(90), deg_to_rad(15)],
##           [1, deg_to_rad(90), deg_to_rad(-15)]])
##    )
##pos = spher_to_cart(                    # sG
##    array([[1, deg_to_rad(0), deg_to_rad(15)],
##           [1, deg_to_rad(60), deg_to_rad(-15)],
##           [1, deg_to_rad(-60), deg_to_rad(-15)]])
##    )
pos = spher_to_cart(                    # sH
    array([[1, deg_to_rad(0), deg_to_rad(-15)],
	   [1, deg_to_rad(60), deg_to_rad(15)],
	   [1, deg_to_rad(-60), deg_to_rad(15)]])
    )

# --- Tilted regular arrays
##pos = spher_to_cart(                    # sI (good for horizontal plane)
##    array([[1, deg_to_rad(0), deg_to_rad(15)],
##           [1, deg_to_rad(0), deg_to_rad(-45)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # sJ
##    array([[1, deg_to_rad(0), deg_to_rad(-15)],
##           [1, deg_to_rad(0), deg_to_rad(45)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )


gain        = -8            # in dB

file_type   = 'wav'         #write file...
encoding    = 'pcm24'
endianness  = 'file'

#kind        = 'listen'          #HRIR kind: 'listen' or 'sphere'
kind        = 'sphere'          #HRIR kind: 'listen' or 'sphere'

file_dir    =   '/Volumes/Audio/test/'      #write file dir

read_file = '/Volumes/Audio/SinS/sweet_1.wav'

# lHRIR
database_dir = '/Users/josephla/Documents/Developer/Listen_hrtf_database/'

subject_id  = '1002'
nfc_r       = 1.95        # distance of lHRIR measurement

# sHRIR
N = 512

if kind is 'listen':
    file_name = 'lHRIR_decoder_test'
elif kind is 'sphere':
    file_name = 'sHRIR_decoder_test'
    

#---------------------------------
# calculate
T = 1./sr
scale = db_to_amp(gain)

# ----- generate file name
write_file  =   file_dir + file_name + '.' + file_type[:3]


# ----- generate lHRIR kernel
if kind is 'listen':
    decoder_kernels = decoder_lHRIR_kernel(pos, k, \
                                           subject_id, database_dir)
elif kind is 'sphere':
    decoder_kernels = sHRIR_decoder_kernel(pos, k, N, 1./sr)
    decoder_kernels *= db_to_amp(-3.5)  #scale to match lHRIR


m = shape(decoder_kernels)[0]   # number of harmonics to convolve
                                # determined from kernel

# ----- shelf filter and NFC decoder kernels
b_format_kernels = transpose(decoder_kernels)

# add zeros for Z, for shelf filtering, etc., if 2D
if m is 3:
    b_format_kernels = array([
        hstack((b_format_kernels[0], zeros((N,1)))),
        hstack((b_format_kernels[1], zeros((N,1))))
        ])

# psychoacoustic shelf and NFC filtering
#       psycho_shelf(x, Wn, k, zi = None)
#       nfc(x, r, T, zi = None)
for i in range(2):
    b_format_kernels[i] = psycho_shelf(
        nfc(b_format_kernels[i], nfc_r, T),
        freq_to_Wn(400., 1./sr),
        ((C.k_2D, C.k_3D)[m is 3])      # switch for 2D or 3D
        )

# remove zeros for Z, if 2D
if m is 3:
    b_format_kernels = array([
        b_format_kernels[0, :, :-1],
        b_format_kernels[1, :, :-1],
        ])

decoder_kernels = transpose(b_format_kernels)

# ************************************************************
# Set up sndfiles for reading:
read_sndfile = Sndfile(read_file)

#--------------------
# read test signal
b = read_sndfile.read_frames(read_sndfile.nframes)

#--------------------
# convolve

res = b_to_binaural(b, decoder_kernels)
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
