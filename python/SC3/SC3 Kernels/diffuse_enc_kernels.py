# ****************************************************************************
# Frequency diffusion encoder kernels
# 
# Generate a set of kernels suitable for encoding mono to
# b-format diffuse kernel encoder.
#
# Kernels are classified by kernel size, N, and stored in directory: 
#
#       '(root or ~)/Library/Application Support/ATK/kernels/FOA/encoders/ \
#           diffuse/[None]/[kernel_size]/[subject_id]'
#
# Within, one four channel [W,X,Y,Z] kernel is found, named:
#
#       'diffuse'
#
# A resulting b-format encode can be generated from these kernels as follows:
#
#       b-format = 'diffuse' * mono
#
#   where * is convolution, + is sum.
#
#
# ****************************************************************************

from muse import *
import os


# ----- functions

def resamp_freq(a):

    res = resample(a/2, len(a)*2)

    return res


def resamp_time(a):

    res = concatenate((
        zeros(len(a)/2),
        a,
        zeros(len(a)/2),
        ))

    return res


def rotate_phase(a, phase):

    tmp_cos = convfilt(a, real(fir_hb(len(a)+1)), 'full') 
    tmp_sin = convfilt(a, imag(fir_hb(len(a)+1)), 'full') 

    res = (cos(phase) * tmp_cos) + (sin(phase) * tmp_sin)
    res = res[len(a)/2:len(a)/2 + len(a)]

    return res
    

def allpass(N, bands, octaves):

    tmp = zeros(bands)

    for octave in range(octaves):
        tmp2 = mirf(
            resamp_freq(
                fir_ap(bands, 2*pi)
                )
            )
        tmp2 *= hann(bands*2)

        for i in range(octave):
            tmp2 = resamp_freq(tmp2)

        tmp = resamp_time(tmp) + tmp2

    # trim
    tmp = tmp[len(tmp)/2 - N/2: len(tmp)/2 + N/2]

    # normalise frequencies
    res = apf(tmp)

    return res


# ----- parameters
srs         = array([None])   # sample rates
Ns          = array([512, 1024, 2048, 4096, 8192, 16384])   # kernel lengths
bands       = array([4, 8, 16, 32, 64, 128])

octaves     = 16

# --
file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/encoders/diffuse'

file_name  = 'diffuse'

# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)


# generate subjects and subject data
subject_ids = []                        # subject ids
for i in range(len(bands)):
    subject_ids += [str(i).zfill(4)]


# ----- loop
for sr in srs:                          # SR

    for N in Ns:                        # kernel sizes
        for subject in arange(len(subject_ids)):  # subject index

            # ----- generate prototype kernels
            pulse = real(fir_hb(N))

            tmp = empty([N, 3])
            for i in range(3):
                tmp[:, i] = rotate_phase(
                    allpass(
                        N,
                        bands[subject],
                        octaves
                        ),
                    uniform(0, 2*pi)
                    )
                

            # ----- create b-format kernel - normalised for diffuse field
            encoder_kernel = interleave(
                array([
                    1./sqrt(2) * pulse,
                    sqrt(3)/2. * tmp[:, 0],
                    sqrt(3)/2. * tmp[:, 1],
                    sqrt(3)/2. * tmp[:, 2]
                    ])
                )


            # ----- generate file name
            write_file = os.path.join(library_dir + database_dir, \
                                      str(sr),  str(N), \
                                      subject_ids[subject], \
                                      file_name + '.' + file_type[:3])


            # ************************************************************
            # Set up sndfiles for writing:

            if not os.path.exists(os.path.dirname(write_file)):
                os.makedirs(os.path.dirname(write_file))
            
            write_sndfile = Sndfile(
                write_file,
                'w',
                Format(file_type, encoding, endianness),
                nchannels(encoder_kernel),
                int(sr or 44100)                    # sr defaults to 44100
                )                                   # if none is supplied

            # ----- write out!
            write_sndfile.write_frames(encoder_kernel)

            # ----- close file
            write_sndfile.close()

