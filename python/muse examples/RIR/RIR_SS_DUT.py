# ******************************************************
# 
# Generate DUT normalisation filter for use with
# log or linear sine sweeps (ESS, LSS) for RIR measurement.
# The measured response and the convolution filter for the DUT
# should be presented. Normally this will be a test SS of a known
# mic and loudspeaker (the ones used for the RIR measurement).
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


from muse import *
import pylab

# kirkeby normalisation filter
# kernel is kernel to normalise
# Wns are LF and HF corner freqs defining band for normalisation
# roll_off is roll off in octaves from Wns (defines transition band)
# Es is the regularization parameter ε(f)
# Note: kirkeby in the form described by Farina results in a 
#       normalised band pass filter between Wns
def kirkeby(kernel, Wns, roll_off = 1./3, Es = array([.01, 10.])):
    
    N = nframes(kernel)         # kernel size
    M = N/2 + 1                 # freqs (fft, has +-freqs)

    # Wn transition bands, for roll off
    c_Wns0 = empty(2)
    c_Wns0[0] = 2**(-roll_off/2) * Wns[0]
    c_Wns0[1] = 2**(roll_off/2) * Wns[0]

    c_Wns1 = empty(2)
    c_Wns1[0] = 2**(-roll_off/2) * Wns[1]
    c_Wns1[1] = 2**(roll_off/2) * Wns[1]

    # design regularization parameter ε(f)
    kirkebyE = ones(M)

    fftWns = (scipy.fftpack.fftfreq(N) * 2)[:M] # only need zero + positive vals

    # test for mode...
    #          mode 1: [True, True], normalise all Wns
    #          mode 2: [False, False], normalise between Wns[0] and Wns[1]
    #          mode 3: [False, True], normalise above Wns[0] only
    #          mode 4: [True, False], normalise below Wns[1] only
    nanWns = isnan(Wns)

    if all(nanWns):             # mode 1
        kirkebyE = Es[0] * kirkebyE
        
    elif all(logical_not(nanWns)): # mode 2
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns0[0]:
                kirkebyE[n] = Es[1]

            elif c_Wns0[0] <= k_Wn and k_Wn < c_Wns0[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns0[0]) * pi) / (c_Wns0[1] - c_Wns0[0])) + 1) / 2 + Es[0]

            elif c_Wns0[1] <= k_Wn and k_Wn < c_Wns1[0]:
                kirkebyE[n] = Es[0]

            elif c_Wns1[0] <= k_Wn and k_Wn < c_Wns1[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns1[0]) * pi) / (c_Wns1[1] - c_Wns1[0]) + pi) \
                         + 1) / 2 + Es[0]

            else:               # c_Wns1[1] <= k_Wn:
                kirkebyE[n] = Es[1]

    elif logical_not(nanWns[0]):    # mode 3
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns0[0]:
                kirkebyE[n] = Es[1]

            elif c_Wns0[0] <= k_Wn and k_Wn < c_Wns0[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns0[0]) * pi) / (c_Wns0[1] - c_Wns0[0])) + 1) / 2 + Es[0]

            else:               # c_Wns0[1] <= k_Wn
                kirkebyE[n] = Es[0]

    else:                       # mode 4
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns1[0]:
                kirkebyE[n] = Es[0]

            elif c_Wns1[0] <= k_Wn and k_Wn < c_Wns1[1]:
                kirkebyE[n] = (Es[1] - Es[0]) * \
                    (cos(((k_Wn - c_Wns1[0]) * pi) / (c_Wns1[1] - c_Wns1[0]) + pi) \
                         + 1) / 2 + Es[0]

            else:               # c_Wns1[1] <= k_Wn
                kirkebyE[n] = Es[1]

    kirkebyE = concatenate((    # complete building kirkebyE by mirroring
            kirkebyE,
            kirkebyE[::-1][:-1]
            ))

    # compute kirkeby "packing" (normalisation) filter

    # 1) take fft, transform to freq domain
    fftH = scipy.fftpack.fft(kernel)

    # 2) compute inverse filter in freq domain
    conjH = conjugate(fftH)
    fftK = conjH / (conjH * fftH + kirkebyE)

    # 3) take ifft, transform back to time domain (and window)
    res = real(scipy.fftpack.ifft(fftK))
    res *= hann(N)

    return res


# modified kirkeby normalisation filter
# kernel is kernel to normalise
# Wns are LF and HF corner freqs defining band for normalisation
# roll_off is roll off in octaves from Wns (defines transition band)
# Note: kirkeby is modified from the form described by Farina, and 
#       instead results in a filter where freq are normalised between
#       LF and HF, but left un-normalised outside the band of interest
#       Do expect that kernel should be length = 2**x+1
def kirkebyM(kernel, Wns, roll_off = 1./3):
    
    N = nframes(kernel)         # kernel size
    M = N/2 + 1                 # freqs (fft, has +-freqs)

    # Wn transition bands, for roll off
    c_Wns0 = empty(2)
    c_Wns0[0] = 2**(-roll_off/2) * Wns[0]
    c_Wns0[1] = 2**(roll_off/2) * Wns[0]

    c_Wns1 = empty(2)
    c_Wns1[0] = 2**(-roll_off/2) * Wns[1]
    c_Wns1[1] = 2**(roll_off/2) * Wns[1]

    # design regularization parameter ε(f)
    kirkebyE = ones(M)

    fftWns = (scipy.fftpack.fftfreq(N) * 2)[:M] # only need zero + positive vals

    # test for mode...
    #          mode 1: [True, True], normalise all Wns
    #          mode 2: [False, False], normalise between Wns[0] and Wns[1]
    #          mode 3: [False, True], normalise above Wns[0] only
    #          mode 4: [True, False], normalise below Wns[1] only
    nanWns = isnan(Wns)

    if all(nanWns):             # mode 1
        kirkebyE = kirkebyE
        
    elif all(logical_not(nanWns)): # mode 2
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns0[0]:
                kirkebyE[n] = 0.

            elif c_Wns0[0] <= k_Wn and k_Wn < c_Wns0[1]:
                kirkebyE[n] = (cos(((k_Wn - c_Wns0[0]) * pi) / (c_Wns0[1] - c_Wns0[0]) + pi) + 1) / 2

            elif c_Wns0[1] <= k_Wn and k_Wn < c_Wns1[0]:
                kirkebyE[n] = 1.

            elif c_Wns1[0] <= k_Wn and k_Wn < c_Wns1[1]:
                kirkebyE[n] = (cos(((k_Wn - c_Wns1[0]) * pi) / (c_Wns1[1] - c_Wns1[0])) + 1) / 2

            else:               # c_Wns1[1] <= k_Wn:
                kirkebyE[n] = 0.

    elif logical_not(nanWns[0]):    # mode 3
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns0[0]:
                kirkebyE[n] = 0.

            elif c_Wns0[0] <= k_Wn and k_Wn < c_Wns0[1]:
                kirkebyE[n] = (cos(((k_Wn - c_Wns0[0]) * pi) / (c_Wns0[1] - c_Wns0[0]) + pi) + 1) / 2

            else:               # c_Wns0[1] <= k_Wn
                kirkebyE[n] = 1.

    else:                       # mode 4
        for n, k_Wn in zip(range(M), fftWns):
            if 0. <= k_Wn and k_Wn < c_Wns1[0]:
                kirkebyE[n] = 1.

            elif c_Wns1[0] <= k_Wn and k_Wn < c_Wns1[1]:
                kirkebyE[n] = (cos(((k_Wn - c_Wns1[0]) * pi) / (c_Wns1[1] - c_Wns1[0])) + 1) / 2

            else:               # c_Wns1[1] <= k_Wn
                kirkebyE[n] = 0.

    kirkebyE = concatenate((    # complete building kirkebyE by mirroring
            kirkebyE,
            kirkebyE[::-1][:-1]
            ))

    # compute kirkeby "packing" (normalisation) filter

    # 1) take fft, transform to freq domain
    fftH = scipy.fftpack.fft(kernel)

    # 2) compute inverse filter in freq domain
    conjH = conjugate(fftH)

    # now interpolate amp to get a smooth transition
    # between the normalised and un-normalised bands
    fftC = conjH / (conjH * fftH + .001)

    magK = abs(fftC) * kirkebyE + ones(N) * (1 - kirkebyE)
    phaK = angle(fftC)

    fftK = magK * (cos(phaK) + sin(phaK) * 1j)

    # 3) take ifft, transform back to time domain (and window)
    res = real(scipy.fftpack.ifft(fftK))
    res *= hann(N)

    return res


# display spectrum
def plot_spec(x, min_db, max_db, sr):

    pylab.plot(
        (scipy.fftpack.fftfreq(nframes(x)) * sr)[:nframes(x)/2],
        clip(amp_to_db(abs(scipy.fftpack.fft(x))), min_db, max_db)[:nframes(x)/2]
        )


# parameters
sr = 44100                      # sr!
bit_depth = 24                  # bit depth
file_format = 'aiff'            # 


file_name = "neu"               # file name, dir & extension
# file_name = "neu-7"             # file name, dir & extension
file_dir = "/Users/josephla/Sound/test/midd/"
file_ext = ".aif"

measured_dur = None              # if none, read all
# measured_dur = 20.              # seconds to read in from the measurement
                                # (Could be tied to rtime and Wns, but
                                # this would also require the sweep method
                                # to be described. Best to determine by
                                # inspection.)
measured_sktm = 0.              # skiptime into the file (secs)

# dutN = 424 + 1                  # just excludes echo
dutN = 2**9 + 1                 # actually better than above (freq of analysis window?)
# dutN = 600 + 1                  # just includes echo: HF comb, +-2dB
# dutN = 674 + 1                  # just includes echo: HF comb, +-2dB
# dutN = 2**10 + 1                # just includes echo: HF comb, +-2dB (similar, but more LF)
# dutN = 2**11 + 1                # just includes echo: HF comb, +-5dB (v big HF and LF ringing)
                                # make this an odd number!

# kfreqs = array([50., 20000.])   # kirkeby cutoff freqs, in Hz
# kfreqs = array([500., 19000.])   # kirkeby cutoff freqs, in Hz
kfreqs = array([50., 19000.])   # kirkeby cutoff freqs, in Hz
roll_off = 1./3                 # kirkeby filter roll off, in octaves
# kmode = 'mod'                   # modified kirkeby filter?
kmode = 'std'

# ******************************************************
# generate file names
measured_file = file_dir + file_name + "_measured" + file_ext
filter_file = file_dir + file_name + "_filter" + file_ext
normal_file = file_dir + file_name + "_normal" + file_ext


# ******************************************************
#  calc vals
k_Wns = freq_to_Wn(kfreqs, 1./sr)
kirN = int(2**ceil(log2(dutN - 1))) + 1 # normalisation filter size


# ******************************************************
# read in measured and filter (signals/sweeps)

# instantiate input soundfiles
measured_sfile =  sndfile(
    measured_file,
    'read'
    )

filter_sfile =  sndfile(
    filter_file,
    'read'
    )

# read in
if measured_dur is None:
    measured = measured_sfile.read_frames(
        measured_sfile.get_nframes()                       # read all frames
        )                                                  # if no dur given

else:
    measured_sfile.seek(dur_to_nframes(measured_sktm, sr)) # skip forward
    measured = measured_sfile.read_frames(                 # then read dur frames
        dur_to_nframes(measured_dur, sr)
        )

filter = filter_sfile.read_frames(
    filter_sfile.get_nframes()                              # read all frames
    )

# ... and close
measured_sfile.close()
filter_sfile.close()


# ******************************************************
# deconvolve (trim to size and window...)

decon = convfilt(measured, filter, 'full')

peakframe = argpeak(decon)      # find peak, for centering

decon = decon[ peakframe - dutN/2 : peakframe + dutN/2 + 1] # trim to dutN around center

decon *= hann(dutN)             # window


# scale to normalise the average gain across kirkeby bandwidth
# decon = scale_to(decon)         # scale to +1
decon *= reciprocal(
    average(
        extract(
            logical_and(
                greater_equal(scipy.fftpack.rfftfreq(dutN) * 2, k_Wns[0] * ones(dutN)),
                less_equal(scipy.fftpack.rfftfreq(dutN) * 2, k_Wns[1] * ones(dutN))),
            2 * abs(scipy.fftpack.rfft(decon))
            )
        )
    )

# zero pad to bring dutN up to kirN
if dutN < kirN:
    decon = concatenate((
            zeros((kirN-dutN) / 2),
            decon,
            zeros((kirN-dutN) / 2)
            ))


# ******************************************************
# kirkeby normalise

if kmode is 'mod':
    normal = kirkebyM(decon, k_Wns, roll_off)

else:
    normal = kirkeby(decon, k_Wns, roll_off)


# ******************************************************
# display
min_db = -24.
max_db = 18.
min_hz = 10.

# display vals*******************
max_amp = max(abs(scipy.fftpack.fft(normal)))
max_Wn = 2 * scipy.fftpack.fftfreq(kirN)[argmax(abs(scipy.fftpack.fft(normal)))]

print "Max normalisation gain = ", amp_to_db(max_amp), "dB at", Wn_to_freq(max_Wn, sr), "Hz"


# Plot IRS***************

# normalise the measured DUT IR
decon_norm = convfilt(decon, normal, 'full')


# DUT impulse response...
ax1 = pylab.subplot(231)
pylab.plot(decon, color = 'b')
ax1.set_ylim(-1., 1.)
pylab.title('DUT IR')

# normalisation filter impulse response...
ax2 = pylab.subplot(232)
pylab.plot(normal, color = 'g')
ax2.set_ylim(-1., 1.)
pylab.title('normalisation IR')

# DUT normalised impulse response...
ax3 = pylab.subplot(233)
pylab.plot(decon_norm, color = 'r')
ax3.set_ylim(-1., 1.)
pylab.title('normalised DUT IR')



# Plot spectra***************

ax4 = pylab.subplot(212)
ax4.set_xscale('log', basex = 2, nonposx='clip')

# spectrum of DUT...
plot_spec(decon, min_db, max_db, sr)

# spectrum of normalisation filter
plot_spec(normal, min_db, max_db, sr)

# spectrum of the normalised DUT
plot_spec(decon_norm, min_db, max_db, sr)

ax4.grid(True)

pylab.title('frequency response')
ax4.set_ylim(min_db, max_db)
ax4.set_xlim(min_hz)

pylab.show()


# ******************************************************
# write normalisation filter out

normal_sfile =  sndfile(
    normal_file,
    'write',
    formatinfo(
        file_format,
        'float32'               # filter written as float
        ),
    1,
    sr
)

# write to output
normal_sfile.write_frames(
    normal,
    kirN
)


# now close files
normal_sfile.close()


# ******************************************************
# finished!

print "\nDone!"
