# **************************************
# 
# binaural decoders
# 
# decode a b-format signal to full 3D
# (periphonic) using one of 3 binaural decoders:
#       Listen HRIR decoder
#       CIPIC HRIR decoder
#       Spherical HRIR decoder
# 
# the resulting decoded signal is generated 'all at once'
# using functional syntax
#
# NOTE: While the measured HRIRs (Listen, CIPIC)
#       are described as 'free field' equalised
#       by their respective authors, there are
#       noticable tonal differences between all
#       three results.
#
# **************************************

from muse import *
import pylab as pl


# params
sr = 44100          # sample rate
k = 1.              # velocity decode (shelf later)


# speaker positions (pairs): un-comment your choice
# NOTE: prefer 'hex' arrangements (3D), especially B3 & C3
#       Also! While the CIPIC database is densely sampled,
#       due to limited sampling on the Y-Z plane (and +-Y)
#       limited arrays corresponding to the Listen
#       database may be easily found. Corresponding
#       arrays are listed below, in the secion labelled
#       '3D (for CIPIC)'.

# ------> 3D (for Sphere and Listen)

# --- Octrahedron arrays
##pos = spher_to_cart(                    # A0
##    array([[1, deg_to_rad(0), deg_to_rad(45)],
##           [1, deg_to_rad(0), deg_to_rad(-45)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # A1
##    array([[1, deg_to_rad(0), deg_to_rad(0)],
##           [1, deg_to_rad(90), deg_to_rad(45)],
##           [1, deg_to_rad(90), deg_to_rad(-45)]])
##    )

# --- 'Regular' arrays
##pos = spher_to_cart(                    # B0
##    array([[1, deg_to_rad(0), deg_to_rad(30)],
##           [1, deg_to_rad(0), deg_to_rad(-30)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # B1
##    array([[1, deg_to_rad(0), deg_to_rad(0)],
##           [1, deg_to_rad(90), deg_to_rad(30)],
##           [1, deg_to_rad(90), deg_to_rad(-30)]])
##    )
##pos = spher_to_cart(                    # B2
##    array([[1, deg_to_rad(0), deg_to_rad(30)],
##           [1, deg_to_rad(60), deg_to_rad(-30)],
##           [1, deg_to_rad(-60), deg_to_rad(-30)]])
##    )
pos = spher_to_cart(                    # B3
    array([[1, deg_to_rad(0), deg_to_rad(-30)],
           [1, deg_to_rad(60), deg_to_rad(30)],
           [1, deg_to_rad(-60), deg_to_rad(30)]])
    )

# --- vertically compressed arrays
##pos = spher_to_cart(                    # C0
##    array([[1, deg_to_rad(0), deg_to_rad(15)],
##           [1, deg_to_rad(0), deg_to_rad(-15)],
##           [1, deg_to_rad(90), deg_to_rad(0)]])
##    )
##pos = spher_to_cart(                    # C1
##    array([[1, deg_to_rad(0), deg_to_rad(0)],
##           [1, deg_to_rad(90), deg_to_rad(15)],
##           [1, deg_to_rad(90), deg_to_rad(-15)]])
##    )
##pos = spher_to_cart(                    # C2
##    array([[1, deg_to_rad(0), deg_to_rad(15)],
##           [1, deg_to_rad(60), deg_to_rad(-15)],
##           [1, deg_to_rad(-60), deg_to_rad(-15)]])
##    )
##pos = spher_to_cart(                    # C3
##    array([[1, deg_to_rad(0), deg_to_rad(-15)],
##	   [1, deg_to_rad(60), deg_to_rad(15)],
##	   [1, deg_to_rad(-60), deg_to_rad(15)]])
##    )


# ------> 3D (for CIPIC)

# --- 'Regular' arrays

###                                       B3c
###   [[  0,     -33.750],
###    [  57.609, 33.134],
###    [ -57.609, 33.134],
##az = rad_to_deg(arctan2(
##    -sin(deg_to_rad(-45)),
##    cos(deg_to_rad(-45)) * cos(deg_to_rad(50.625))
##    ))
##el = rad_to_deg(arcsin(cos(deg_to_rad(-45)) * sin(deg_to_rad(50.625))))
##pos = spher_to_cart(
##    array([[1, deg_to_rad(0), deg_to_rad(-33.75)],
##           [1, deg_to_rad(az), deg_to_rad(el)],
##           [1, deg_to_rad(-az), deg_to_rad(el)]])
##    )
##
# --- vertically compressed arrays
###                                       C2c
###   [[  0,       16.875],
###    [  58.304, -15.687],
###    [ -58.304, -15.687],
##az = rad_to_deg(arctan2(
##    -sin(deg_to_rad(-55)),
##    cos(deg_to_rad(-55)) * cos(deg_to_rad(28.125))
##    ))
##el = rad_to_deg(arcsin(cos(deg_to_rad(-55)) * sin(deg_to_rad(28.125))))
##pos = spher_to_cart(
##    array([[1, deg_to_rad(0), deg_to_rad(16.875)],
##           [1, deg_to_rad(az), deg_to_rad(-el)],
##           [1, deg_to_rad(-az), deg_to_rad(-el)]])
##    )
###                                       C3c
###   [[  0,      -16.875],
###    [  58.304,  15.687],
###    [ -58.304,  15.687],
##pos = spher_to_cart(
##    array([[1, deg_to_rad(0), deg_to_rad(-16.875)],
##           [1, deg_to_rad(az), deg_to_rad(el)],
##           [1, deg_to_rad(-az), deg_to_rad(el)]])
##    )


gain        = -12           # in dB

file_type   = 'wav'         # write file...
encoding    = 'pcm24'
endianness  = 'file'

# set kind for Listen, CIPIC or Sphere HRIR
#kind        = 'listen'
#kind        = 'cipic'
kind        = 'sphere'

file_dir    =   '/Volumes/Audio/test/'      #write file dir

read_file = '/Volumes/Audio/SinS/sweet_1.wav'


# HRIR dirs
l_database_dir = '/Users/josephla/Documents/Developer/Listen_hrtf_database/'
c_database_dir = '/Users/josephla/Documents/Developer/CIPIC_hrtf_database/standard_hrir_database/'

l_subject_id  = '1002'
c_subject_id = '021'        # KEMAR head

nfc_r = {                   # NFC distance (for compensation)
    'listen'    : 1.95,     # distance of lHRIR measurement
    'cipic'     : 1.00,     # distance of lHRIR measurement
    }

# sHRIR
N = 512

file_name = {
    'listen'    : 'lHRIR_3D_decoder_test',
    'cipic'     : 'cHRIR_3D_decoder_test',
    'sphere'    : 'sHRIR_3D_decoder_test',
    }


# calculate
T = 1./sr
scale = db_to_amp(gain)

# ----- generate file name
write_file  =   file_dir + file_name[kind] + '.' + file_type[:3]


# ----- generate lHRIR kernel
if kind is 'listen':
    decoder_kernels = lHRIR_decoder_kernel(pos, k, \
                                           l_subject_id, l_database_dir)
    decoder_kernels *= db_to_amp(5)     # scale to match other HRIRs

elif kind is 'cipic':
    decoder_kernels = cHRIR_decoder_kernel(pos, k, \
                                           c_subject_id, c_database_dir)
elif kind is 'sphere':
    decoder_kernels = sHRIR_decoder_kernel(pos, k, N, 1./sr)


m = shape(decoder_kernels)[0]   # number of harmonics to convolve
                                # determined from kernel

# ----- shelf filter and NFC decoder kernels
b_format_kernels = transpose(decoder_kernels)

# psychoacoustic shelf and NFC filtering
if kind is 'sphere':                # psycho_shelf only for sphere
    for i in range(2):
        b_format_kernels[i] = psycho_shelf(
            b_format_kernels[i],
            freq_to_Wn(700., 1./sr),
            C.k_3D                  # 3D shelf
            )
else:
    for i in range(2):              # psycho_shelf and nfc for measured HRIRs
        b_format_kernels[i] = psycho_shelf(
            nfc(b_format_kernels[i], nfc_r[kind], T),
            freq_to_Wn(700., 1./sr),
            C.k_3D                  # 3D shelf
            )

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

##zi = zeros((m, N - 1, 2))
##res, zf = b_to_binaural(b, decoder_kernels, zi = zi)

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
