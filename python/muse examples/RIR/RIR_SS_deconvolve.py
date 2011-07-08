# ******************************************************
# 
# Deconvolve log or linear sine sweeps (ESS, LSS) for RIR measurement
# includes the optional use of DUT (mic/LS/measurement equipment) normalisation.
# 
# The following files are required ( and exptected to be named):
# 
#       measured SS signal               :          NAME_measured
#       computed deconv filter           :          NAME_filter
#       computed DUT normalisation filter:          DUT_NAME_normal (optional)
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


rtime = 1.5                     # expected rt, in secs
ris = .01                       # envelope rise time
dec = .5                        # envelope decay time
                                # env looks like [ris, rtime, dec]

gain = -3.                      # gain, in dB

measured_dur = None              # if none, read all
# measured_dur = 20.              # seconds to read in from the measurement
                                # (Could be tied to rtime and Wns, but
                                # this would also require the sweep method
                                # to be described. Best to determine by
                                # inspection.)
measured_sktm = 0.              # skiptime into the file (secs)

dut_phase = 'nor'               # phase of DUT filter
                                # 'nor': corrects phase and gain
                                # 'lin': corrects gain only, phase response is linear
                                # 'min': corrects gain only, phase response is minimum

file_name = "A_BL"              # file name, dir & extension
dut_file_name = "neu-7"         # DUT, may be assigned to None

file_dir = "/Users/josephla/Sound/test/midd/"
file_ext = ".aif"


# ******************************************************
# generate file names
measured_file = file_dir + file_name + "_measured" + file_ext
filter_file = file_dir + file_name + "_filter" + file_ext
rir_file = file_dir + file_name + "_rir" + file_ext

if dut_file_name:
    normal_file = file_dir + dut_file_name + "_normal" + file_ext


# ******************************************************
# calculate signal parms


# ******************************************************
# read in measured and filter (signals/sweeps)

# instantiate input soundfiles and read in... and close
measured_sfile =  sndfile(
    measured_file,
    'read'
    )

filter_sfile =  sndfile(
    filter_file,
    'read'
    )

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


measured_sfile.close()
filter_sfile.close()


# do the same for DUT normal, if it exists

if dut_file_name:
    normal_sfile =  sndfile(
        normal_file,
        'read'
        )

    normal = normal_sfile.read_frames(
        normal_sfile.get_nframes() # read all frames
        )

    normal_sfile.close()


# ******************************************************
# adjust DUT normalisation filter phase response
if dut_file_name:
    if dut_phase is 'lin':
        normal = linf(normal)
        print "DUT normalise gain only, linear phase response."

    elif dut_phase is 'min':
        normal = minf(linf(normal))
        print "DUT normalise gain only, minimum phase response."

    else :
        print "DUT normalise phase and gain."


# ******************************************************
# deconvolve to find rir

if dut_file_name:
    filter = convfilt(filter, normal, 'full')

rir = scale_to(
    convfilt(measured, filter, 'full'),
    db_to_amp(gain)
)


# ******************************************************
# trim and envelope rir

peakframe = argpeak(rir) / nchannels(rir) # find peak

ris_nframes = dur_to_nframes(ris, sr)
dec_nframes = dur_to_nframes(dec, sr)
rtime_nframes = dur_to_nframes(rtime, sr)
dur_nframes = ris_nframes + dec_nframes + rtime_nframes

rir = fadvancen(rir, peakframe - ris_nframes)[:dur_nframes]

if nchannels(rir) is 1:
    rir *= linen([ris_nframes, dur_nframes - dec_nframes, dur_nframes], 'hann')
else:
    rir *= interleave(linen([ris_nframes, dur_nframes - dec_nframes, dur_nframes], 'hann'))


# ******************************************************
# write room impulse response (rir)

rir_sfile =  sndfile(
    rir_file,
    'write',
    formatinfo(
        file_format,
        encoding
        ),
    nchannels(rir),
    sr
)

# write to output
rir_sfile.write_frames(
    rir,
    nframes(rir)
)


# now close files
rir_sfile.close()


# ******************************************************
# finished!

print "\nDone!"
