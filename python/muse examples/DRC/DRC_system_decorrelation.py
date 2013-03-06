# -*- coding: utf-8 -*-
# ******************************************************
# 
# Deconvolve log or linear sine sweeps (ESS, LSS) for RIR measurement
# includes the use of DUT (mic/LS/measurement equipment) normalisation.
# 
# See: Angelo Farina, "Simultaneous Measurement of Impulse Response and Distortion
# with a Swept-Sine Technique," in  (presented at the Audio Engineering Society
# Convention no. 108, Audio Engineering Society, 2000).
# 
# AND
# 
# Angelo Farina, “Advancements in Impulse Response Measurements by Sine Sweeps,” in
# (presented at the Audio Engineering Society 122nd Convention, Vienna, Austria:
# Audio Engineering Society, 2007).

# decorrelation filters



from muse import *
import os.path
import shutil
import pylab



# A simple linear phase IIR filter
# order = 2*N
def lowpass(x, N, Wn):

    res = x
    for i in range(2):
        res = fiir_lp(res[::-1], N, Wn)

    return res


def hipass(x, N, Wn):

    res = x
    for i in range(2):
        res = fiir_hp(res[::-1], N, Wn)

    return res


# NOTE: windowing only handles even lengths
def octave_scaled_hann(N, octave):

    res = over_dub(
        zeros(N),
        hann(N * 2**-octave),
        N/2 - (N/2) * 2**-octave,
        True
        )

    return res


# Apply an octave scaled Hann window to x 
# NOTE: windowing only handles even lengths at the moment
#       can use conv x 2 w/ imag(fir_hb)
def octave_hann_windowed(x, N, octaves):

    # append zeros - for filtering
    M = nframes(x)
    
    x = concatenate((
        zeros(M),
        x,
        zeros(M)
        ))


    # cascade filter and window
    res = zeros(M)
    lp0 = x

    for octave in arange(octaves - 1, -1, -1):

        Wn = 2**-(octaves - octave)
        lp1 = lowpass(lp0, N, Wn)

        res += octave_scaled_hann(M, octave) * \
            (lp0 - lp1)[M:2*M]

        lp0 = lp1

    # add last octave...
    res += hann(M) * lp0[M:2*M]

    return res


# N should be odd! (2**x + 1)
def allpass(N, order, octaves, phase):

    res = apf(
        octave_hann_windowed(
            fir_ap(N-1),
            order,
            octaves
            )
        )

    # rotate phase
    tmp_cos = convfilt(res, real(fir_hb(N+1)), 'full') 
    tmp_sin = convfilt(res, imag(fir_hb(N+1)), 'full') 

    res = (cos(phase) * tmp_cos) + (sin(phase) * tmp_sin)
    res = res[N/2:N/2 + N]

    return res


# display spectrum
def plot_spec(x, min_db, max_db, sr):

    pylab.plot(
        (scipy.fftpack.fftfreq(nframes(x)) * sr)[:nframes(x)/2],
        clip(amp_to_db(abs(scipy.fftpack.fft(x))), min_db, max_db)[:nframes(x)/2]
        )


# ******************************************************
# parameters
##sr = 192000                     # sr! (highest of Fireface 800)
sr = 96000                     # sr! (highest of Fireface 800)
#sr = 48000                     # sr! (highest of Fireface 800)
file_format = 'wav'             # 
file_encoding = 'pcm32'         # note: applies to the normalisation filter


xfreq = 700.                 # crossover freq, in Hz
#xfreq = None                # single band

xN = 1                       # crossover order (* 2)

#apN = 2**15+1                 # allpass N (or read from file length?)
#apN = 2**7 + 1

winN = 2                    # octave windowing order
#winN = 1                    # octave windowing order

#octaves = 12                # number of windowed octaves
octaves = 11                # number of windowed octaves
#octaves = 10                # number of windowed octaves



file_name = 'drc'           # file name & dir
working_dir = '/Users/josephla/Sound/Seattle Feb-Mar 2013/Turnkey DRC'
filter_dir = 'xover'
decor_dir = 'decor'

#num_speakers = 2
#num_speakers = 4
#num_speakers = 12
#num_speakers = 24
num_speakers = 28



# ******************************************************
# generate file names
file_ext = file_format[:3]


# ******************************************************
# update output director if dual band...
# dual band?
if xfreq is not None:
    decor_dir = decor_dir + "_" + str(int(xfreq))


# ******************************************************
# copy data output text files

# check if output path exists, and create if no
if not os.path.exists(os.path.join(working_dir, str(sr), decor_dir)):
    os.makedirs(os.path.join(working_dir, str(sr), decor_dir))

shutil.copy(
    os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_distances" + os.extsep + "txt"
        ),
    os.path.join(
        working_dir, str(sr), decor_dir
        ),
    )

shutil.copy(
    os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_delays" + os.extsep + "txt"
        ),
    os.path.join(
        working_dir, str(sr), decor_dir
        ),
    )

shutil.copy(
    os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_gains" + os.extsep + "txt"
        ),
    os.path.join(
        working_dir, str(sr), decor_dir
        ),
    )


# ******************************************************
# loop over each speaker!


for spkr in range(num_speakers):

    print "\nProcessing speaker :", spkr


    # ******************************************************
    # generate file names
    drc_file = os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_" + str(spkr).zfill(2) + os.extsep + file_ext
        )

    decor_file = os.path.join(
        working_dir, str(sr), decor_dir, \
        file_name + "_" + str(spkr).zfill(2) + \
        os.extsep + file_ext
        )


    # ******************************************************
    # read in:
    #   drc
    # ... and set up output soundfile for decorrelated version

    # instantiate input and output soundfiles

    print 'Reading: ' + drc_file

    drc_sfile = Sndfile(
        drc_file,
        'r'
        )

    decor_sfile =  Sndfile(
        decor_file,
        'w',
        drc_sfile.format,        # these values should match
        drc_sfile.channels,      # those found in 'parameters'
        drc_sfile.samplerate     # above
        )


    # read in here:
    drc = drc_sfile.read_frames(
        drc_sfile.nframes            # read all frames
        )

    # ... and close
    drc_sfile.close()


    # ******************************************************
    #  decorrelate

    # generate full-band windowed allpass filters
    apN = len(drc) + 1

    ap = allpass(apN, winN, octaves, random.uniform(-pi, pi))

    # dual band?
    if xfreq is not None:

        Wn = freq_to_Wn(xfreq, 1./sr)

        impulse = zeros(apN)
        impulse[apN/2] = 1

        ap = apf(
            lowpass(impulse, xN, Wn) +
            hipass(ap, xN, Wn)
            )


    # convolve & trim
    decor = convfilt(drc, ap, 'full')
    decor = decor[len(drc)/2 : len(drc)/2 + len(drc)]

    
    # ******************************************************
    # write to output

    decor_sfile.write_frames(decor)

    # remember to close!!
    decor_sfile.close()



# ******************************************************
# display
min_db = -24.
max_db = 18.
min_hz = 10.

### display vals*******************
##max_Wn = 2 * scipy.fftpack.fftfreq(kirN)[argmax(abs(scipy.fftpack.fft(norm_filter)))]
##
##print "Max Kirkeby gain = ", amp_to_db(max_amp), "dB at", Wn_to_freq(max_Wn, sr), "Hz"
##
##
### Plot IRS***************
##
### normalise the measured RIR
##normalised_response = convfilt(response, norm_filter, 'full')
##normalised_response = normalised_response[kirN/2:kirN/2 + kirN]
##
##
### impulse response...
##ax1 = pylab.subplot(231)
##pylab.plot(response, color = 'b')
##ax1.set_ylim(-1., 1.)
##pylab.title('RIR')
##
### normalisation filter impulse response...
##ax2 = pylab.subplot(232)
##pylab.plot(norm_filter, color = 'g')
##ax2.set_ylim(-1., 1.)
##pylab.title('normalisation IR')
##
### DUT normalised impulse response...
##ax3 = pylab.subplot(233)
##pylab.plot(normalised_response, color = 'r')
##ax3.set_ylim(-1., 1.)
##pylab.title('normalised RIR')
##
##
##
### Plot spectra***************
##
##ax4 = pylab.subplot(212)
##ax4.set_xscale('log', basex = 2, nonposx='clip')
##
### spectrum of DUT...
##plot_spec(response, min_db, max_db, sr)
##
### spectrum of normalisation filter
##plot_spec(norm_filter, min_db, max_db, sr)
##
### spectrum of the normalised DUT
##plot_spec(normalised_response, min_db, max_db, sr)
##
##ax4.grid(True)
##
##pylab.title('frequency response')
##ax4.set_ylim(min_db, max_db)
##ax4.set_xlim(min_hz)
##
##pylab.show()



# ******************************************************
# finished!

print "\nDone!"
