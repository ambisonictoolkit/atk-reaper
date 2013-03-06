# -*- coding: utf-8 -*-
# ******************************************************
# 
# Deconvolve log sine sweep (LSS) for DRC
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


# NOTE: average gain is weighted by Q


from muse import *
import os.path
#import pylab

# kirkeby normalisation filter
# kernel is kernel to normalise
# Wns are LF and HF corner freqs defining band for normalisation
# roll_off is roll off in octaves from Wns (defines transition band)
# Es is the regularization parameter ε(f)
# Note: kirkeby in the form described by Farina results in a 
#       normalised band pass filter between Wns
#
#       Include this in muse.filters

# NOTE: Looks like the kirkeby filter is expecting an odd length
#       input. Need to document this before adding to Muse


# included in muse


# A simple linear phase IIR filter
# order = 2*N
def lowpass(x, N, Wn):

    res = x
    for i in range(2):
        res = fiir_lp(res[::-1], N, Wn)

    return res


# An octave scaled Hann window
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


# display spectrum
def plot_spec(x, min_db, max_db, sr):

    pylab.plot(
        (scipy.fftpack.fftfreq(nframes(x)) * sr)[:nframes(x)/2],
        clip(amp_to_db(abs(scipy.fftpack.fft(x))), min_db, max_db)[:nframes(x)/2]
        )


# ******************************************************
# parameters
##sr = 192000                     # sr!
sr = 96000                     # sr!
file_format = 'wav'             # 
file_encoding = 'pcm32'         # note: applies to the normalisation filter

                                # (make this an odd number)
#kirN = 2**14 + 1                   # size (samps) of DUT IR analysis window
kirN = 2**15 + 1                   # size (samps) of DUT IR analysis window
#kirN = 2**12 + 1                   # size (samps) of DUT IR analysis window
#drcN = kirN - 1
drcN = (kirN - 1) / 2           # temporary!!! bug in Convolution2.ar


#kfreqs = array([20., sr / 2])   # kirkeby cutoff freqs, in Hz (for mains)
#kfreqs = array([10., 150.])   # kirkeby cutoff freqs, in Hz (for subs)
#kfreqs = array([10., 200.])   # kirkeby cutoff freqs, in Hz (for subs)
kfreqs = concatenate((
    interleave(ones(24)) * array([20., sr / 2]),   # kirkeby cutoff freqs, in Hz (for mains)
    interleave(ones(4)) * array([10., 200.]),   # kirkeby cutoff freqs, in Hz (for subs)
    ))

roll_off = 1./3                 # kirkeby filter roll off, in octaves

file_name = 'drc'           # file name & dir
dut_norm_name = 'dut_norm'
dut_delay_name = 'dut_delay'    # txt file with delay in samples
drc_measure_gains_name = 'drc_measure_gains' # txt file with measurement gains
working_dir = '/Users/josephla/Sound/Seattle Feb-Mar 2013/Turnkey DRC'
filter_dir = 'filters'
signal_dir = 'signals'


#num_speakers = 2
#num_speakers = 4
#num_speakers = 12
#num_speakers = 24
num_speakers = 28


#normalise_dut = True            # use DUT normalisation filter?
normalise_dut = False            # use DUT normalisation filter?


norm_ref = False                # normalise to reference? (or to weighted band average)

#order = 1                       # order of octave windowing (*2)
order = 2                       # order of octave windowing (*2)

#octaves = 12                    # number of octaves to window
#octaves = 11                    # number of octaves to window
octaves = 10                    # number of octaves to window
#octaves = 9                    # number of octaves to window


# ******************************************************
# generate file names
file_ext = file_format[:3]

# measurement sweeps and dut filter
dut_norm_file = os.path.join(
    working_dir, str(sr), filter_dir, \
    dut_norm_name + os.extsep + file_ext
    )
signal_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_signal" + os.extsep + file_ext
    )
decon_filter_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_decon" + os.extsep + file_ext
    )

# data input text files
dut_delay_file = os.path.join(
    working_dir, str(sr), filter_dir, \
    dut_delay_name + os.extsep + "txt"
    ) 
drc_measure_gains_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    drc_measure_gains_name + os.extsep + "txt"
    ) 

# data output text files
distances_file = os.path.join(
    working_dir, str(sr), filter_dir, \
    file_name + "_distances" + os.extsep + "txt"
    ) 
delays_file = os.path.join(
    working_dir, str(sr), filter_dir, \
    file_name + "_delays" + os.extsep + "txt"
    ) 
gains_file = os.path.join(
    working_dir, str(sr), filter_dir, \
    file_name + "_gains" + os.extsep + "txt"
    ) 


# ******************************************************
# read in:
#       DUT delay and measured gains...

dut_delay = int(loadtxt(dut_delay_file))
measurement_gains = loadtxt(drc_measure_gains_file)


# ******************************************************
# read in:
#       DUT normalisation kernel
#       signal (sweep)
#       deconvolution filter (sweep)

# instantiate input and output soundfiles

print 'Reading: ' + dut_norm_file
print 'Reading: ' + signal_file
print 'Reading: ' + decon_filter_file

dut_norm_sfile = Sndfile(
    dut_norm_file,
    'r'
    )
signal_sfile = Sndfile(
    signal_file,
    'r'
    )
decon_filter_sfile = Sndfile(
    decon_filter_file,
    'r'
    )

# read in here:
dut_norm = dut_norm_sfile.read_frames(
    dut_norm_sfile.nframes            # read all frames
    )
signal = signal_sfile.read_frames(
    signal_sfile.nframes            # read all frames
    )
decon_filter = decon_filter_sfile.read_frames(
    decon_filter_sfile.nframes          # read all frames
    )

# ... and close
dut_norm_sfile.close()
signal_sfile.close()
decon_filter_sfile.close()


# ******************************************************
# deconvolve to find DRC reference

print "De-convolving reference..."
reference = convfilt(signal, decon_filter, 'full')
print "... now finished de-convolution."


# ******************************************************
# set up drc delay, gain and distance mesurements
distances = zeros(num_speakers)
normal_delays = zeros(num_speakers)
normal_dbs = zeros(num_speakers)


# ******************************************************
# loop over each speaker!


for spkr in range(num_speakers):

    print "\nProcessing speaker :", spkr


    # ******************************************************
    # generate file names
    measured_file = os.path.join(
        working_dir, str(sr), signal_dir, \
        file_name + "_measure_" + str(spkr).zfill(2) + os.extsep + file_ext
        )

    norm_file = os.path.join(
        working_dir, str(sr), filter_dir, \
        "drc_" + str(spkr).zfill(2) + os.extsep + file_ext
        )


    # ******************************************************
    # read in:
    #   measured sweep
    # ... and set up output soundfile for normalisation filter

    # instantiate input and output soundfiles

    print 'Reading: ' + measured_file

    measured_sfile = Sndfile(
        measured_file,
        'r'
        )

    norm_sfile =  Sndfile(
        norm_file,
        'w',
        Format(file_format, file_encoding),        # these values should match
        signal_sfile.channels,      # those found in 'parameters'
        signal_sfile.samplerate     # above
        )


    # read in here:
    measured = measured_sfile.read_frames(
        measured_sfile.nframes            # read all frames
        )

    # ... and close
    measured_sfile.close()


    # ******************************************************
    # remove DUT delay
    measured = fadvancen(           # remove DUT delay (hardware)
        measured,
        dut_delay
    )



    # ******************************************************
    # deconvolve to find RIR response

    print "De-convolving speaker", spkr, "..."

    if normalise_dut:           # correct response via DUT filter (if used!)
        measured = convfilt(
            fadvancen(           # remove DUT correction filter delay (if used!)
                measured,
                int(dut_norm_sfile.nframes)/2
            ),
            dut_norm,
            'full'
            )

    response = convfilt(
        measured,
        decon_filter,
        'full'
        )

    print "... now finished de-convolution."



    # ******************************************************
    #  calc vals for delay / distance
    delay = float(argpeak(response) - argpeak(reference)) / sr
    distance = speed_of_sound * delay

    normal_delays[spkr] = delay
    distances[spkr] = distance

    print "Delay :", delay, "seconds"
    print "Distance :", distance, "meters"



    # ******************************************************
    #  trim and window response

    response = response[
        argpeak(response) - (kirN / 2):argpeak(response) + (kirN / 2) + 1
        ]

    response = octave_hann_windowed(response[:-1], order, octaves)
    response = concatenate((response, zeros(1)))



    # ******************************************************
    # normalise gain

    # 1st, normalise test gain
    response *= reciprocal(db_to_amp(measurement_gains[spkr]))

    # Kirkeby normalisation freqs, in Wn
    k_Wns = freq_to_Wn(kfreqs[spkr], 1./sr)


    # to peak
    if norm_ref:
        normalisation_amp = reciprocal(peak(reference))

    # scale to normalise the average gain across kirkeby bandwidth
    else:
        response_rfft = scipy.fftpack.rfft(response)
        response_bin_mask = logical_and(
            greater_equal(scipy.fftpack.rfftfreq(kirN) * 2, k_Wns[0] * ones(kirN)),
            less_equal(scipy.fftpack.rfftfreq(kirN) * 2, k_Wns[1] * ones(kirN)))

        Wn0_bin = argpeak(greater_equal(scipy.fftpack.rfftfreq(kirN) * 2, k_Wns[0] * ones(kirN)))
        Wn1_bin = argmin(less_equal(scipy.fftpack.rfftfreq(kirN) * 2, k_Wns[1] * ones(kirN))) - 1
        if Wn1_bin == -1:
            Wn1_bin = kirN - 1

        average_amp = average(
            extract(
                response_bin_mask,
                2 * abs(response_rfft)
                ),
            0,
            1./(arange(Wn0_bin, Wn1_bin + 1))
            )

        normalisation_amp = reciprocal(average_amp)

    response *= normalisation_amp


    # ******************************************************
    # design kirkeby normalisation filter
    norm_filter = kirkeby(response, k_Wns, roll_off)


    # ******************************************************
    # write to output (scaled by gain)

    max_amp = max(abs(scipy.fftpack.fft(norm_filter)))

    # convert to length drcN
    norm_sfile.write_frames(
        convfilt(
            norm_filter * reciprocal(max_amp),
            real(fir_hb(drcN)),
            'full'
            )[drcN / 2: drcN / 2 + drcN]
        )

    # remember to close!!
    norm_sfile.close()


    # NOTE: final amplitude scaling... convert to dB and track
    #       this is what we'll need to write out...
    normal_db = amp_to_db(normalisation_amp * max_amp)
    normal_dbs[spkr] = normal_db

    print "Normalisation gain :", normal_db, "dB"


# ******************************************************
# scale normalisation gains and delays
normal_dbs -= max(normal_dbs) 
normal_delays = max(normal_delays) - normal_delays


# ******************************************************
# write out distance, gain, delay
savetxt(distances_file, distances, '%3.12f')
savetxt(delays_file, normal_delays, '%3.12f')
savetxt(gains_file, normal_dbs, '%3.12f')


### ******************************************************
### display
##min_db = -24.
##max_db = 18.
##min_hz = 10.
##
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
