# -*- coding: utf-8 -*-
# ******************************************************
# 
# Batch cross-over for DRC filters
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

# Cross-over system responses



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


xfreq = 80.                 # crossover freq, in Hz

lp_hp = concatenate((       # flags for lp or hp
    zeros(24),              # highpass
    ones(4)                 # lowpass
    ))

N = 2                       # crossover order (* 2)

file_name = 'drc'           # file name & dir
working_dir = '/Users/josephla/Sound/Seattle Feb-Mar 2013/Turnkey DRC'
filter_dir = 'filters'
xover_dir = 'xover'

#num_speakers = 2
#num_speakers = 4
#num_speakers = 12
#num_speakers = 24
num_speakers = 28



# ******************************************************
# generate file names
file_ext = file_format[:3]


# ******************************************************
# copy data output text files

# check if output path exists, and create if no
if not os.path.exists(os.path.join(working_dir, str(sr), xover_dir)):
    os.makedirs(os.path.join(working_dir, str(sr), xover_dir))

shutil.copy(
    os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_distances" + os.extsep + "txt"
        ),
    os.path.join(
        working_dir, str(sr), xover_dir
        ),
    )

shutil.copy(
    os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_delays" + os.extsep + "txt"
        ),
    os.path.join(
        working_dir, str(sr), xover_dir
        ),
    )

shutil.copy(
    os.path.join(
        working_dir, str(sr), filter_dir, \
        file_name + "_gains" + os.extsep + "txt"
        ),
    os.path.join(
        working_dir, str(sr), xover_dir
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

    xover_file = os.path.join(
        working_dir, str(sr), xover_dir, \
        file_name + "_" + str(spkr).zfill(2) + os.extsep + file_ext
        )


    # ******************************************************
    # read in:
    #   measured sweep
    # ... and set up output soundfile for normalisation filter

    # instantiate input and output soundfiles

    print 'Reading: ' + drc_file

    drc_sfile = Sndfile(
        drc_file,
        'r'
        )

    xover_sfile =  Sndfile(
        xover_file,
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
    #  crossover

    # pad with zeros
    xover = concatenate((
        zeros(int(drc_sfile.nframes)),
        drc,
        zeros(int(drc_sfile.nframes))
        ))
    
    if lp_hp[spkr] == 0:
        xover = hipass(xover, N, freq_to_Wn(xfreq, 1./sr))
        print "hp"
    else:
        xover = lowpass(xover, N, freq_to_Wn(xfreq, 1./sr))
        print "lp"

    # trim
    xover = xover[int(drc_sfile.nframes) : 2 * int(drc_sfile.nframes)]

    
    # ******************************************************
    # write to output

    xover_sfile.write_frames(xover)

    # remember to close!!
    xover_sfile.close()



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
