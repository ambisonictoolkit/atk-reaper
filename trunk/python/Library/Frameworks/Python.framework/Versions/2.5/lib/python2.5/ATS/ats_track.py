# python 2.5

# ATS
# by Juan Pampin
# pampin@u.washington.edu
#
# ported to python by
# Johnathan Lyon
# jgl@u.washington.edu
# 2007-06-06

#
# FILE: ats_track.py
#
# This file contains the tracker function
# definitions and constants


# standard libraries
from math import floor, sin, atan, cos

# external libraries
import pyaudiolab
from numpy import array, zeros
from numpy.numarray.fft import fft, inverse_fft

# libraries from this package
from ats_util import db_amp, ppp2, compute_frames, init_sound, find_peak_by_track, \
						optimize_sound, add_sound, simple_plot, TWO_PI, max_array
from ats_struct import ats_sound, ats_fft
from ats_win import make_blackman_window, make_fft_window, window_norm
from ats_peak_detect import peak_detection, sort_peaks_by_freq
from ats_critical_bands import evaluate_smr
from ats_peak_track import update_tracks, peak_tracking
from ats_residual import compute_residual



def tracker (file, snd,
			 	start=0.0,
			 	duration=None,
			 	lowest_frequency=20, # must be > 0
			 	highest_frequency=20000.0,
			 	frequency_deviation=0.1,
			 	window_cycles=4,
			 	window_type='blackman-harris-4-1',
			 	hop_size=0.25,
			 	fft_size=None,
			 	lowest_magnitude=db_amp(-60),
			 	track_length=3,
			 	min_segment_length=3,
			 	last_peak_contribution=0.0,
			 	SMR_continuity=0.0,
			 	SMR_threshold=None,
			 	amp_threshold=None,
			 	residual=None,
			 	par_energy=True,
			 	optimize=True,
			 	debug=None,
			 	verbose=None,
			 	force_M=None,
			 	force_window=None,
			 	window_mu=0.0,
			 	window_beta=1.0):
	#input file
	fil = pyaudiolab.sndfile(file,'read')
	#ATS sound
	sound = ats_sound(snd)
	# file sampling rate
	file_sampling_rate = fil.get_samplerate()
	file_format = fil.get_file_format()
	file_encoding = fil.get_encoding()
	# make sure file is mono
	if (fil.get_channels() != 1): 
		print 'ERROR in Tracker: not a mono file'
		return
	# index of first sample to read
	st = int(floor(start * file_sampling_rate))
	# index of last sample to read
	nd = 0
	if (duration==None):
		nd = fil.get_nframes()
	else:
		nd = int(st + (duration * file_sampling_rate))
	# convert input file to indexable array
	filtemp = fil.read_frames(fil.get_nframes())
	fil.close()
	fil = filtemp
	# number of samples to read
	total_samps = nd - st
	# file duration
	file_duration = (1.0 * total_samps) / file_sampling_rate
	# number of samples in a cycle
	cycle_samps = int(floor( (1.0 / lowest_frequency) * window_cycles * file_sampling_rate))
	# we want an odd lengthed window centered at time 0.0
	M=0
	if (force_M==None):
		if ((cycle_samps%2)==0):
			M = cycle_samps + 1
		else:
			M = cycle_samps
	else:
		M = force_M
	# fft size is next power of 2 or forced by user
	N=0
	if (fft_size==None):
		N = ppp2(2 * M)
	else:
		N = fft_size
	# fft structure
	fft_struct = ats_fft(N, file_sampling_rate, zeros(N, 'Complex64'))
	# window array
	window = []
	if (force_window==None):
		if (window_type.startswith('blackman')):
			window = make_blackman_window(window_type, M)
		else:
			window = make_fft_window(window_type, M, window_beta, window_mu)
	else:
		window = force_window
	# window normalization
	norm = window_norm(window)
	# hop in samples
	hop = int(floor(M * hop_size))
	# number of analysis frames
	frames = compute_frames(total_samps, hop, st, nd)
	# we keep sample numbers of central points of the windows
	win_samps = [0] * frames
	# magic number for fft frequencies (freq resolution)
	fft_mag = float(file_sampling_rate) / N
	# lowest frequency to analyze
	l_Frq = lowest_frequency
	if (l_Frq < 0.0):
		l_Frq = 0.0
	# highest frequency to analyze
	h_Frq = highest_frequency
	if ((h_Frq > (file_sampling_rate / 2)) or (h_Frq <= l_Frq)):
		h_Frq = int(floor(file_sampling_rate / 2))
	# lowest bin to read
	lowest_bin = int(floor(l_Frq / fft_mag))
	# highest bin to read
	highest_bin = int(floor(h_Frq / fft_mag))
	# Lists for data
	# list of lists for peaks
	ana_frames = []
	for i in range(0, frames):
		ana_frames.append([])
	# misc vars
	# timer
	tmp = 0.0
	smp = 0
	# central point of the window
	M_over_2 = int(floor( (M - 1) / 2))
	# first point in fft buffer where to write
	first_point = N - M_over_2
	# set file pointer half a window from the first sample
	filptr = st - M_over_2
	# minimum SMR
	min_smr = 0.0
	if (SMR_threshold != None):
		min_smr = SMR_threshold
	n_partials = 0
	tracks = []
	peaks = []
	unmatched_peaks = []
	
	# POST INIT STATUS
	# tell user we are starting partial tracking
	print "frames = " + str(frames)
	print "M = " + str(M)
	print "N = " + str(N)
	print "tracking..."
	
	# MAIN LOOP
	frame_n = 0
	while(frame_n < frames):
		# clear fft lists
		fft_struct.fd = zeros(N, 'Complex64')
		# multiply by window
		for k in range(0, M):
			#if ((filptr >= 0) and (filptr < len(fil))):			
			ind = int((k + first_point) % N)
			if ((filptr >= 0) and (filptr < len(fil))):
				fft_struct.fd[ind] = window[k] * fil[filptr]	
			else:
				fft_struct.fd[ind] = 0.0				
			filptr = filptr + 1
		# note that after the loop filptr=M and not M-1
		# the sample at the middle of the window is:
		smp = filptr - M_over_2 - 1
		if (debug):
			print "smp = " + str(smp)
		# we keep sample numbers of window midpoints in an array
		win_samps[frame_n] = smp
		# set timer
		tmp = float(smp - st) / file_sampling_rate
		# get the dft
		fft_struct.fd = fft(fft_struct.fd)
		
		# =====================================
		# Peak Detection:
		# get peaks (amplitudes normalized by window norm)
		# list of peaks is sorted by frequency
		peaks = peak_detection(fft_struct, 
								lowest_bin=lowest_bin, 
								highest_bin=highest_bin, 
								lowest_magnitude=lowest_magnitude, 
								norm=norm)
						
		# process peaks
		if (len(peaks) > 0):
			evaluate_smr(peaks)
			
			# ======================================
			# Peak Tracking:
			# try to match peaks
			# only if we have at least 2 frames
			# and if we have active tracks			
			if (frame_n > 0):
				tracks = update_tracks(tracks, track_length, frame_n, ana_frames, last_peak_contribution)
				if (len(tracks)>0):
					cpy_peak = None
					# tracks peaks and get leftover
					unmatched_peaks = peak_tracking(tracks, peaks, frequency_deviation, SMR_continuity, unmatched=[])
					# kill unmatched peaks from previous frame
					for k in unmatched_peaks[0]:
						cpy_peak = k
						cpy_peak.amp = 0.0
						cpy_peak.smr = 0.0
						peaks.insert(0, cpy_peak)
					# give birth to peaks from new frame
					for k in unmatched_peaks[1]:
						k.track = n_partials
						n_partials = n_partials + 1
						cpy_peak = k
						cpy_peak.amp = 0.0
						cpy_peak.smr = 0.0
						i = int(frame_n - 1)
						ana_frames[i].insert(0, cpy_peak)
						tracks.insert(0, k)	
			else:
				# give number to all peaks
				for p in peaks:
					p.track = n_partials
					n_partials = n_partials + 1
			ana_frames[frame_n] = peaks		
		filptr = filptr - M + hop
		if (verbose):
			print "<Frame: " + str(frame_n) + " Time: " + str(tmp) + " Tracks: " + str(n_partials)
		frame_n = frame_n + 1	
	# INIT ATS SOUND
	init_sound(sound, file_sampling_rate, hop, M, frames, file_duration, n_partials)
	# fill it up with data
	for k in range(0, n_partials):
		for frame in range(0, frames):
			pe = find_peak_by_track(k, ana_frames[frame])
			if (pe != None):
				sound.amp[k][frame] = pe.amp
				sound.frq[k][frame] = pe.frq
				sound.pha[k][frame] = pe.pha
			# set time anyways
			sound.time[k][frame] = float(win_samps[frame] - st) / file_sampling_rate
			
			
			######## ^
			######## |
			######## bug here
		##########
	
	# finally optimize and declare new sound in ATS
	if (optimize):
		optimize_sound(sound, min_frq=lowest_frequency, max_frq=highest_frequency, 
						min_length=min_segment_length, amp_threshold=amp_threshold, 
						verbose=verbose)
	if (verbose):
		print "Partials: " + str(sound.partials) + " Frames: " + str(sound.frames)
	# register sound in the system
	add_sound(sound)
	# now get the residual
	if (residual!=None):
		format = pyaudiolab.formatinfo(file_format, file_encoding)
		compute_residual(fil, residual, sound, win_samps, file_sampling_rate, verbose, fmt=format)
	
	#DEBUGS
	
"""
		residual_analysis(residual, sound, par_energy=par_energy, verbose=verbose, \
							debug=debug, equalize=True)
	close_input(fil)
	output to global variable
	"""
	
	
	
# TESTS

def test(dur=2.0):
	tracker('/clarinet.aif','test', duration=dur, residual=None, debug=False, verbose=True, optimize=True)
	
	# residual synthesis is fucked up
	# optimize is broken -> [ [list], int] instead of [ [list], [list]...]
	# output not interpolating frames
	# fix add sounds as dictionary
	
def interp(ind, arr):
	if (ind >= len(arr)):
		return arr[len(arr)]
	elif (ind < 0):
		return arr[0]
	else:
		lo = int(floor(ind))
		hi = lo + 1
		loval = arr[lo]
		hival = arr[hi]
		rang = hival - loval
		scale = ind - lo
		return (scale * rang) + loval
	
	
# interpolation of values in <env> of
# format [index0, value0, { index1, value1 }* ]
# based on a given index <val>
def envelope_interp(val, env):
	#check to see if env is in pairs
	if (((len(env)%2) == 1) or (len(env)==0)):
		print "ERROR: incorrect envelope format"
		return None
	# check if exceeds bounds of envelope
	if (val <= env[0]):
		return env[1]
	elif (val >= env[int(len(env) - 2)]):
		return env[int(len(env) - 1)]
	else:
		for i in range(2, int(len(env)), 2):
			if (env[i] == val):
				return env[i + 1]
			elif (env[i] > val):
				ind = val - env[int(i-2)]
				x = env[i] - env[int(i-2)]
				y = env[int(i + 1)] - env[int(i - 1)]
				return ((ind / x) * y) + env[int(i - 1)]	
	
