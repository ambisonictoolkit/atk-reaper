# ******************************************************
# 
# Generate log or linear sine sweeps (ESS, LSS) for RIR measurement
# along with inverse/deconvolution filter (sweep)
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

# parameters
sr = 44100                      # sr!
bit_depth = 24                  # bit depth
file_format = 'aiff'            # 
encoding = 'pcm24'              # note: only applies to the signal, not the deco filter



freq0 = sr/2.**12               # choice of freq0 determines log or lin (= 0.) sweep
# freq0 = 0.

# roll_off = None                   # upper and lower freq roll off in octaves
roll_off = 1./7                   # for log, applies only to lower freq
                                  # for lin, applies only to upper freq
                                  # used for time domain windowing (filter)
                                  # roll_off may be set to None, for no envelope
                                  # note: enhanced pre and post ringing may result
                                  #       if roll_off isn't used

rtime = 1.                      # expected rt, in secs

gain = -3.                      # gain, in dB

repeats = 1                     # number of repetitions of the signal
# repeats = 2                     # number of repetitions of the signal
# repeats = 3                     # number of repetitions of the signal
                                # note: averaging of repetition is considered
                                # to be deprecated for the sine sweep method

file_name = "test"              # file name, dir & extension
file_dir = "/Users/josephla/Sound/test/"
file_ext = ".aif"


# ******************************************************
# generate file names
signal_file = file_dir + file_name + "_signal" + file_ext
filter_file = file_dir + file_name + "_filter" + file_ext


# ******************************************************
# calculate signal parms

freqs = array([freq0, sr/2.])   # signal start and end freqs
Wns = freq_to_Wn(freqs, 1./sr)  # signal freqs as Wns

if freq0 is 0.:
    method = 'lin'
else:
    method = 'log'

# signal duration
if method is 'lin':
    signal_range = freqs[1] - freqs[0]

    signal_dur = rtime * log(signal_range) # this choice is ad hoc
    print "\nsignal duration =", signal_dur, "secs"

elif method is 'log':
    signal_octs = log2(freqs[1] / freqs[0])

    signal_dur = rtime * signal_octs
    print "\nsignal duration =", signal_dur, "secs"

signal_nframes = dur_to_nframes(signal_dur, sr)


# ******************************************************
# generate signal

if method is 'lin':
    signal = sinosc(lin(Wns, signal_nframes))

elif method is 'log':
    signal = sinosc(
        generators.exp(Wns[::-1], signal_nframes)
        )[::-1] # generate reversed, to guarantee 0 as last val
                # a la Farina, to reduce HF pre-ring

# scale it to gain
signal *= db_to_amp(gain)


# ******************************************************
# apply roll off filter to signal

if roll_off is not None:        # apply roll off filter?

    # roll off corner frequencies
    c_Wns = array([
            2**roll_off * Wns[0],
            2**-roll_off * Wns[1]
            ])

    # generate the filter, in the time domain, and apply (envelope)
    if method is 'lin':
        n1 = ceil((c_Wns[1] - Wns[1]) / (Wns[0] - Wns[1]) * (signal_nframes - 1))
        dec_nframes = ceil(pi * n1 / arccos(1. - C.sqrt2) + 1)

        dec = ((1 - cos(lin([0., pi], dec_nframes))) / 2)[::-1]

        signal *= concatenate(
                (
                    ones(signal_nframes - dec_nframes),
                    dec
                    )
                )

    elif method is 'log':
        n0 = ceil((log2(c_Wns[0]) - log2(Wns[0])) / (log2(Wns[1]) - log2(Wns[0])) * (signal_nframes - 1))
        ris_nframes = ceil(pi * n0 / arccos(1. - C.sqrt2) + 1)
        dec_nframes = ris_nframes

        ris = (1 - cos(lin([0., pi], ris_nframes))) / 2
        dec = ones(ris_nframes) # replace w/ no fade-out, a la Farina to reduce HF pre-ring

        signal *= concatenate(
                (
                    ris,
                    ones(signal_nframes - (ris_nframes + dec_nframes)),
                    dec
                    )
                )


# ******************************************************
# generate deconvolution filter
#       Kirkeby time-packing filter should be used in the DUT regularisation
#       process, but not here....

# time reverse
filter = copy(signal[::-1])

# -6dB/oct filter for log
if method is 'log':
    filter *= db_to_amp(
        lin(
            array([0., -6 * signal_octs]),
            signal_nframes
            )
        )


# ******************************************************
# zero pad

signal = concatenate((
        signal,
        zeros(dur_to_nframes(rtime, sr))
        ))


# ******************************************************
# repeats
# note: signal_nframes isn't updated here,
#       so lists the frames of a single
#       sweep w/out end zero padding or
#       repeats

signal = tile(signal, repeats)



# ******************************************************
# display

pylab.plot(signal)
pylab.plot(filter)
pylab.show()


# ******************************************************
# write signal and filter out

# instantiate soundfiles for writing
signal_sfile =  sndfile(
    signal_file,
    'write',
    formatinfo(
        file_format,
        encoding
        ),
    1,
    sr
)

filter_sfile =  sndfile(
    filter_file,
    'write',
    formatinfo(
        file_format,
        'float32'               # filter written as float
        ),
    1,
    sr
)

# write to output
signal_sfile.write_frames(
    signal,
    nframes(signal)
)

filter_sfile.write_frames(
    filter,
    nframes(filter)
)


# now close files
signal_sfile.close()
filter_sfile.close()


# ******************************************************
# finished!

print "\nDone!"
