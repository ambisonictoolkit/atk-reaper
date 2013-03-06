# -*- coding: utf-8 -*-
# ******************************************************
# 
# Generate DUT normalisation filter for use with
# log sine sweeps (LSS) for RIR measurement.
# The measured response and the convolution filter for the DUT
# should be presented. Normally this will be a test SS of a known
# mic and loudspeaker in an anecoic chamber (the ones used for the
# RIR measurement).
#
# In a practical case, this measurement may just be a 'feed through'
# measurement of the DA/AD itself
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

# reports hardware delay...
#
# NOTE: gain is normalised to peak of DUT response

from muse import *
import os.path
import pylab


# display spectrum
def plot_spec(x, min_db, max_db, sr):

    pylab.plot(
        (scipy.fftpack.fftfreq(nframes(x)) * sr)[:nframes(x)/2],
        clip(amp_to_db(abs(scipy.fftpack.fft(x))), min_db, max_db)[:nframes(x)/2]
        )


# ******************************************************
# parameters
sr = 96000                     # sr!
file_format = 'wav'             # 
file_encoding = 'pcm32'         # note: applies to the normalisation filter

dutN = 2**14 + 1                   # size (samps) of DUT IR analysis window
                                # NOTE: must be an odd number
kirN = dutN                                

#kfreqs = array([20., sr / 2])   # kirkeby cutoff freqs, in Hz
kfreqs = array([15., sr / 2])   # kirkeby cutoff freqs, in Hz
#kfreqs = array([10., sr / 2])   # kirkeby cutoff freqs, in Hz
roll_off = 1./3                 # kirkeby filter roll off, in octaves

file_name = 'dut'           # file name & dir
working_dir = '/Users/josephla/Sound/Seattle Feb-Mar 2013/Turnkey DRC'
filter_dir = 'filters'
signal_dir = 'signals'


gain = 0                       # scale DUT normalisationi filter by
                                # gain in dB on writing normalisation filter out

norm_ref = True                # normalise to reference? (or to weighted band average)


# ******************************************************
# generate file names
file_ext = file_format[:3]

signal_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_signal" + os.extsep + file_ext
    )
measured_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_measure" + os.extsep + file_ext
    )
decon_filter_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_decon" + os.extsep + file_ext
    )
norm_file = os.path.join(
    working_dir, str(sr), filter_dir, \
    file_name + "_norm" + os.extsep + file_ext
    )
delay_file = os.path.join(
    working_dir, str(sr), filter_dir,
    file_name + "_delay" + os.extsep + "txt"
    ) 

# check if output path exists, and create if no
if not os.path.exists(os.path.join(working_dir, str(sr), filter_dir)):
    os.makedirs(os.path.join(working_dir, str(sr), filter_dir))


# ******************************************************
# read in:
#   test signal sweep
#   deconvolution filter (sweep)
#   measured sweep
#
# ... and set up output soundfile for normalisation filter

# instantiate input and output soundfiles

print 'Reading: ' + signal_file
print 'Reading: ' + measured_file
print 'Reading: ' + decon_filter_file

signal_sfile = Sndfile(
    signal_file,
    'r'
    )
measured_sfile = Sndfile(
    measured_file,
    'r'
    )
decon_filter_sfile = Sndfile(
    decon_filter_file,
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
signal = signal_sfile.read_frames(
    signal_sfile.nframes            # read all frames
    )
measured = measured_sfile.read_frames(
    measured_sfile.nframes            # read all frames
    )
decon_filter = decon_filter_sfile.read_frames(
    decon_filter_sfile.nframes          # read all frames
    )

# ... and close
signal_sfile.close()
measured_sfile.close()
decon_filter_sfile.close()


# ******************************************************
# deconvolve to find (measurement reference and) DUT response

print "Actually starting to convolve..."

ref = convfilt(signal, decon_filter, 'full')
dut = convfilt(measured, decon_filter, 'full')

print "... now finished convolving."



# ******************************************************
#  calc vals
k_Wns = freq_to_Wn(kfreqs, 1./sr)

dutDel = argpeak(dut) - argpeak(ref)

print "Delay of DUT at SR =", sr, ":", dutDel

# write out delay
savetxt(delay_file, array([dutDel]), '%d')


# ******************************************************
#  trim and window DUT response
dut_response = dut[argpeak(dut) - (dutN / 2):argpeak(dut) + (dutN / 2) + 1]
dut_response *= win(dutN)


# ******************************************************
# normalise gain

# to peak
if norm_ref:
    dut_response *= reciprocal(peak(ref))

# scale to normalise the average gain across kirkeby bandwidth
else:
    dut_response_rfft = scipy.fftpack.rfft(dut_response)
    dut_response_bin_mask = logical_and(
        greater_equal(scipy.fftpack.rfftfreq(dutN) * 2, k_Wns[0] * ones(dutN)),
        less_equal(scipy.fftpack.rfftfreq(dutN) * 2, k_Wns[1] * ones(dutN)))

    Wn0_bin = argpeak(greater_equal(scipy.fftpack.rfftfreq(dutN) * 2, k_Wns[0] * ones(dutN)))
    Wn1_bin = argmin(less_equal(scipy.fftpack.rfftfreq(dutN) * 2, k_Wns[1] * ones(dutN))) - 1
    if Wn1_bin == -1:
        Wn1_bin = kirN - 1

    average_amp = average(
        extract(
            dut_response_bin_mask,
            2 * abs(dut_response_rfft)
            ),
        0,
        1./(arange(Wn0_bin, Wn1_bin + 1))
        )

    dut_response *= reciprocal(average_amp)


# ******************************************************
# design kirkeby normalisation filter


# NOTE: The kirkeby filter is expecting an odd length input.
norm_filter = kirkeby(dut_response, k_Wns, roll_off)


# ******************************************************
# write to output  scaled so peak gain = 0db + gain
#norm_sfile.write_frames(db_to_amp(gain) * norm_filter)

max_amp = max(abs(scipy.fftpack.fft(norm_filter)))

norm_sfile.write_frames(db_to_amp(gain)/max_amp * norm_filter)


# remember to close!!
norm_sfile.close()




# ******************************************************
# display
min_db = -24.
max_db = 18.
min_hz = 10.

# display vals*******************
max_Wn = 2 * scipy.fftpack.fftfreq(kirN)[argmax(abs(scipy.fftpack.fft(norm_filter)))]

print "Max normalisation gain = ", amp_to_db(max_amp), "dB at", Wn_to_freq(max_Wn, sr), "Hz"


# Plot IRS***************

# normalise the measured DUT IR
normalised_dut = convfilt(dut_response, norm_filter, 'full')


# DUT impulse response...
ax1 = pylab.subplot(231)
pylab.plot(dut_response, color = 'b')
ax1.set_ylim(-1., 1.)
pylab.title('DUT IR')

# normalisation filter impulse response...
ax2 = pylab.subplot(232)
pylab.plot(norm_filter, color = 'g')
ax2.set_ylim(-1., 1.)
pylab.title('normalisation IR')

# DUT normalised impulse response...
ax3 = pylab.subplot(233)
pylab.plot(normalised_dut, color = 'r')
ax3.set_ylim(-1., 1.)
pylab.title('normalised DUT IR')



# Plot spectra***************

ax4 = pylab.subplot(212)
ax4.set_xscale('log', basex = 2, nonposx='clip')

# spectrum of DUT...
plot_spec(dut_response, min_db, max_db, sr)

# spectrum of normalisation filter
plot_spec(norm_filter, min_db, max_db, sr)

# spectrum of the normalised DUT
plot_spec(normalised_dut, min_db, max_db, sr)

ax4.grid(True)

pylab.title('frequency response')
ax4.set_ylim(min_db, max_db)
ax4.set_xlim(min_hz)

pylab.show()



# ******************************************************
# finished!

print "\nDone!"
