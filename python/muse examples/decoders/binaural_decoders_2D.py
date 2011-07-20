# **************************************
# 
# binaural decoders
# 
# decode a b-format signal to horizontal 2D
# (pantophonic) using one of 3 binaural decoders:
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


# ------> 2D
# speaker positions (pairs): un-comment your choice

##pos = pol_to_cart(                    # 1st choice, but incompatable 
##    array([[1, deg_to_rad(0)],        # with cHRIR
##           [1, deg_to_rad(90)]])
##    )
pos = pol_to_cart(                    # 2nd choice
    array([[1, deg_to_rad(45)],
           [1, deg_to_rad(-45)]])
    )
##pos = pol_to_cart(                    # 3rd choice
##    array([[1, deg_to_rad(30)],
##           [1, deg_to_rad(-30)]])
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

# sHRIR kernel length
N = 512

file_name = {
    'listen'    : 'lHRIR_2D_decoder_test',
    'cipic'     : 'cHRIR_2D_decoder_test',
    'sphere'    : 'sHRIR_2D_decoder_test',
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

# add zeros for Z, for shelf filtering, NFC,
# as the 2D generated kernels don't include Z
b_format_kernels = array([
    hstack((b_format_kernels[0], zeros((shape(b_format_kernels)[1],1)))),
    hstack((b_format_kernels[1], zeros((shape(b_format_kernels)[1],1))))
    ])

# psychoacoustic shelf and NFC filtering
if kind is 'sphere':                # psycho_shelf only for sphere
    for i in range(2):
        b_format_kernels[i] = psycho_shelf(
            b_format_kernels[i],
            freq_to_Wn(700., 1./sr),
            C.k_2D                  # 2D shelf
            )
else:
    for i in range(2):              # psycho_shelf and nfc for measured HRIRs
        b_format_kernels[i] = psycho_shelf(
            nfc(b_format_kernels[i], nfc_r[kind], T),
            freq_to_Wn(700., 1./sr),
            C.k_2D                  # 2D shelf
            )

# remove zeros for Z, as the resulting kernels shouldn't include Z
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
