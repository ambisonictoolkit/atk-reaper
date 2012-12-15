# ****************************************************************************
# Frequency spreading encoder kernels
# 
# Generate a set of kernels suitable for encoding mono to
# b-format Spread kernel encoder.
#
# Kernels are classified by kernel size, N, and stored in directory: 
#
#       '(root or ~)/Library/Application Support/ATK/kernels/FOA/encoders/ \
#           spread/[sr]/[kernel_size]/[subject_id]'
#
# Within, one four channel [W,X,Y,Z] kernel is found, named:
#
#       'spread'
#
# A resulting b-format encode can be generated from these kernels as follows:
#
#       b-format = 'spread' * mono
#
#   where * is convolution, + is sum.
#
#
# ****************************************************************************

from muse import *
import os


# ----- functions

# 2nd order allpass
# NOTE: q = 1./sqrt(2) gives same as fiir_ap(x, 2...)
def fiir_ap2(x, Wn, q, zi = None):

    Wb = Wn / q

    c = (tan(pi*Wb/2.) - 1) / (tan(pi*Wb/2.) + 1)
    d = -cos(pi*Wn)

    a = array([1, d*(1-c), -c])
    b = array([-c, d*(1-c), 1])
    
    res = ffilter(b, a, x, zi)

    return res


def cos_filts(x, Wn, q, N):

    ap = x
    for i in range(N):
        ap = fiir_ap2(ap, Wn, q)
    res = 0.5 * (real(ap[::-1] + ap) - 1j * imag(ap[::-1] - ap))

    return res

def sin_filts(x, Wn, q, N):

    ap = x
    for i in range(N):
        ap = fiir_ap2(ap, Wn, q)
    res = 0.5 * (imag(ap[::-1] + ap) + 1j * real(ap[::-1] - ap))

    return res


def normalise_xyz(x):

    nor_fac = array([sqrt(2), 1, 1, 1])

    fft_x = empty(shape(x), 'complex128')
    fft_x = fft_x[:nframes(x)/2 + 1]


    for i, k in zip(range(nchannels(x)), nor_fac):
        fft_x[:, i] = rfft(k * x[:, i])

    mag_dir = zeros(nframes(x)/2 + 1)
    for i in range(1, nchannels(x)):
        mag_dir += abs(fft_x[:, i])**2
    mag_dir = sqrt(mag_dir)

    if mag_dir[-1] == 0:
        mag_dir[-1] = 1

    mag_nor = abs(fft_x[:, 0]) / mag_dir
    mag_nor[-1] = 1

    res = x
    for i in range(1, nchannels(x)):
        res[:, i] = irfft(
            mag_nor * fft_x[:, i]
            )

    return res


def ae(x, Wn):

    phase = fphasor(Wn, 0, nframes(x))
    phase -= phase[-1]/2
    probe = cos(phase)

    res = aed((x * interleave(probe)).sum(0) * ones_like(x))[:-1]

    return res


def centre(x, Wn):

    x_aed = ae(x, Wn)

    res = tumble(
        rotate(
            x,
            -x_aed[0]
            ),
        -x_aed[1]
        )

    return res
    

# ----- parameters
srs         = array([44100, 48000, 88200, 96000, 192000])   # sample rates
Ns          = array([512, 1024, 2048, 4096, 8192, 16384])   # kernel lengths
bpos        = arange(1, 13 + 1)

freq0 = 20.     # LF band edge, in Hz

M = Ns[-1]       # working kernel size

# --
file_type   = 'wav'         # write file...
encoding    = 'pcm32'
endianness  = 'file'


user_dir        = True                              # write library to user dir?
library_dir     = '/Library/Application Support/ATK'      # library location
database_dir    = '/kernels/FOA/encoders/spread'

file_name  = 'spread'

# is user dir set?
if user_dir:
    library_dir = os.path.expanduser('~' + library_dir)


# generate subjects and subject data
subject_ids = []                        # subject ids
for i in range(len(bpos)):
    subject_ids += [str(i).zfill(4)]


# ----- loop
for sr in srs:                          # SR

    # ----- set up bandwidth (in Wn) & determine octaves
    Wn_band = array([freq_to_Wn(freq0, 1./sr), 1.])
    octaves = (freq_to_oct(Wn_to_freq(Wn_band, sr)) * array([-1, 1])).sum()

    # ----- find cf (in Wn) & determine q
    Wn = freq_to_Wn(
        oct_to_freq(
            octaves/2 + freq_to_oct(Wn_to_freq(Wn_band[0], sr))
            ),
        1./sr
        )
    q = 1./(octaves / 2) 

    for N in Ns:                        # kernel sizes
        for subject in arange(len(subject_ids)):  # subject index

            # -----initial calcs
            rotations = (bpos[subject] * octaves/4)
            theta_rots = int(ceil(rotations/2))
            phi_rots = theta_rots - 1


            # ----- generate prototype kernels (complex)
            pulse_hb = fir_hb(M)

            cos_phi_hb = cos_filts(pulse_hb, Wn, q, phi_rots)
            sin_phi_hb = sin_filts(pulse_hb, Wn, q, phi_rots)

            cos_phi_cos_theta_hb = cos_filts(cos_phi_hb, Wn, q, theta_rots)
            cos_phi_sin_theta_hb = sin_filts(cos_phi_hb, Wn, q, theta_rots)

            # ----- take real and create b-format kernel
            encoder_kernel = interleave(
                array([
                    1./sqrt(2) * real(pulse_hb),
                    real(cos_phi_cos_theta_hb),
                    real(cos_phi_sin_theta_hb),
                    real(sin_phi_hb)
                    ])
                )

            # ----- trim to N
            encoder_kernel = encoder_kernel[(M-N)/2:(M-N)/2 + N]

            # ----- normalise freq response (restore LFs)
            encoder_kernel = normalise_xyz(encoder_kernel)

            # ----- centre response about Wn
            encoder_kernel = tilt(encoder_kernel, pi/2)
            encoder_kernel = centre(encoder_kernel, Wn)
            encoder_kernel = mirror_y(encoder_kernel)


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

