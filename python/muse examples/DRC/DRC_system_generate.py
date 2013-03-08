# -*- coding: utf-8 -*-
# ******************************************************
# 
# Generate log sine sweep (LSS) for DRC measurement
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

# NOTES (check if this is true!!):
#   Also, for some reason this is only 'working' with EPD 6.3-2 (32-bit). Need
#   to have a look to see why it isn't behaving w/ newer EPD.


from muse import *
import os.path
import pylab

# parameters
sr = 96000                      # sr!
bit_depth = 24                  # bit depth
file_format = 'wav'             # 
file_encoding = 'pcm24'         # note: only applies to the signal
                                # not the deco filter. Deco filter
                                # is set to 'pcm32

freqs = array([sr/2.**13, sr/2.])   # signal start and end freqs
#freqs = array([sr/2.**14, sr/2.])

# roll_off = None                   # upper and lower freq roll off in octaves
roll_off = 1./7                   # if freqs[1] != sr/2, applies to freqs[0]
                                  # used for time domain windowing (filter)
                                  # roll_off may be set to None, for no envelope
                                  # note: enhanced pre and post ringing may result
                                  #       if roll_off isn't used

rtime = 1.5                     # expected rt, in secs

gain = -3.                      # gain, in dB

repeats = 1                     # number of repetitions of the signal
                                # note: averaging of repetition is considered
                                # to be deprecated for the sine sweep method
                                # See Farina for comments on averaging!

sltN = 2048                     # number of samples of pre/post noise burst
# sltN = 0                      # 'slate'. =0 --> NO pre/post slate


file_name = 'drc'           # file name & dir
working_dir = '/Users/josephla/Sound/Seattle Feb-Mar 2013/Turnkey DRC'
filter_dir = 'filters'
signal_dir = 'signals'


display = False                  # display waveforms?
#display = True                  # display waveforms?
                                # NOTE: EPD 7.2-2 doesn't successfully display
                                # OverflowError: Agg rendering complexity
                                # exceeded. Consider downsampling or decimating
                                # your data. 

# ******************************************************
# generate file names
file_ext = file_format[:3]

signal_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_signal" + os.extsep + file_ext
    )
decon_filter_file = os.path.join(
    working_dir, str(sr), signal_dir, \
    file_name + "_decon" + os.extsep + file_ext
    )


# ******************************************************
# calculate signal parms

Wns = freq_to_Wn(freqs, 1./sr)  # signal freqs as Wns

# signal duration
signal_octs = log2(freqs[1] / freqs[0])

signal_dur = rtime * signal_octs
print "\nsignal duration =", signal_dur, "secs"

signal_nframes = dur_to_nframes(signal_dur, sr)


# ******************************************************
# generate signal

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
    n0 = ceil((log2(c_Wns[0]) - log2(Wns[0])) / (log2(Wns[1]) - log2(Wns[0])) * (signal_nframes - 1))
    ris_nframes = ceil(pi * n0 / arccos(1. - C.sqrt2) + 1)
    dec_nframes = ris_nframes

    ris = (1 - cos(lin([0., pi], ris_nframes))) / 2

    if Wns[1] < 1.:
        dec = ris[::-1]         # apply a roll off for HFs

    else:
        dec = ones(ris_nframes) # replace w/ no fade-out,
                                # a la Farina to reduce HF pre-ring

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
decon_filter = copy(signal[::-1])

# -6dB/oct filter for log
decon_filter *= db_to_amp(
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
# slate with a noise burst
if sltN != 0:
    signal = concatenate((
        db_to_amp(gain) * white(sltN),
        zeros(dur_to_nframes(rtime, sr)),
        signal,
        zeros(dur_to_nframes(rtime, sr)),
        db_to_amp(gain) * white(sltN),
        ))


# ******************************************************
# display

if display:
    pylab.plot(signal)
    pylab.plot(decon_filter)
    pylab.show()


# ***********# ******************************************************
# write signal and filter out

# check if output path exists, and create if no
if not os.path.exists(os.path.join(working_dir, filter_dir)):
    os.mkdir(os.path.join(working_dir, filter_dir))

signal_sfile =  Sndfile(
    signal_file,
    'w',
    Format(
        file_format,
        file_encoding
        ),
    1,
    sr
)

decon_filter_sfile =  Sndfile(
    decon_filter_file,
    'w',
    Format(
        file_format,
        'pcm32'                 # decon_filter written as pcm32
        ),
    1,
    sr
)


# write to output
signal_sfile.write_frames(signal)
decon_filter_sfile.write_frames(decon_filter)


# now close files
signal_sfile.close()
decon_filter_sfile.close()


# ******************************************************
# finished!

print "\nDone!"
