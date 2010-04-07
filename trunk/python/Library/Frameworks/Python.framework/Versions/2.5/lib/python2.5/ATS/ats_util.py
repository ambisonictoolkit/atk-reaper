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
# FILE: ats_util.py
#
# This file contains the utilities function
# definitions and constants
#
# includes:
#	constants:
#		TWO_PI
#		ATS_MAX_DB_SPL
#		ATS_CRITICAL_BANDS
#		ATS_MIN_SEGMENT_LENGTH
#		ATS_AMP_THRESHOLD
#
#	functions:
#		ppp2
#		clone_list
#		init_sound
#		find_peak_by_track
#		db_amp
#		amp_db
#		amp_db_spl
#		max_array
#		un_norm_array
#		prom_array
#		compute_frames
#		optimize_sound
#		get_ampmax
#		get_frqmax
#		fill_sound_gaps
#		get_gaps_arr
# 		segments_arr
#		find_next_val_arr
#		envelope_interp
#		trim_partials
#		remove_short_segments
#		set_amp_av
#		set_frq_av
#		get_valid_partials
#		simplify_sound
#		sort_partials_by_freq
#		simple_plot
#		plot_2
#		plot_3





from math import pi, floor, log
from pylab import plot, xlabel, ylabel, grid, show
from numpy import zeros



# ===================================
# CONSTANTS

TWO_PI = 2 * pi
ATS_MAX_DB_SPL = 100.0
ATS_CRITICAL_BANDS = 25
ATS_MIN_SEGMENT_LENGTH = 3
ATS_AMP_THRESHOLD = -60

# variable to keep names of loaded sounds
ATS_SOUNDS = []


# ====================================
# FUNCTIONS

# the closest power of 2
def ppp2(num,dep=0):
	tmp = pow(2,dep)
	if (tmp >= num):
		return tmp
	else:
		dep = dep + 1
		return ppp2(num, dep)

# clones a list		
def clone_list(l):
	n = zeros(1,'Float64')
	if (type(n)==type(l)):
		if (n.dtype == l.dtype):
			n = zeros(len(l), 'Float64')
		else:
			n = zeros(len(l), 'Complex64')
		for k in range(0, len(l)):
			n[k] = l[k]
		return n
	else:
		n = []
		for k in l:
			n.append(k)
		return n
						
# initializes an ATS sound
def init_sound (sound, sampling_rate, frame_size, window_size, frames, duration, partials,
				has_phase=True, has_noise=False, bands = ATS_CRITICAL_BANDS):
	sound.sampling_rate = sampling_rate
	sound.frame_size = frame_size
	sound.window_size = window_size
	sound.partials = partials
	sound.frames = frames
	sound.dur = duration
	sound.time = [0] * partials
	sound.frq_av = [0.0] * partials
	sound.amp_av = [0.0] * partials
	sound.frq = [0] * partials
	sound.amp = [0] * partials
	if (has_phase):
		sound.pha = [0] * partials
	if (has_noise):
		sound.band_energy = [0] * bands
	for tr in range(0, partials):
		sound.time[tr] = zeros(frames, 'Float64')
		sound.frq[tr] = zeros(frames, 'Float64')
		sound.amp[tr] = zeros(frames, 'Float64')
		if (has_phase):
			sound.pha[tr] = zeros(frames, 'Float64')
		if (has_noise and (tr < bands)):
			sound.band_energy[tr] = zeros(frames, 'Float64')
				
# left-to-right search for a track number in a peaklist		
def find_peak_by_track (track, peaklist):
	for k in peaklist:
		if (track==k.track):
			return k
	return None

def add_sound(sound):
	if ( ATS_SOUNDS.count(sound.name) < 1):
		ATS_SOUNDS.append(sound.name)

# ====================================
# functions for amplitude conversion

def db_amp(db):
	return pow(10.0, (db / 20.0))
	
def amp_db(amp):
	return 20 * log(amp, 10)

def amp_db_spl (db_spl):
	return db_amp(db_spl - ATS_MAX_DB_SPL)


# ====================================
# c-fun.cl
# all these functions were supported as C functions
# using Lisp's ffi. In this new version of ATS we
# use their Python version instead.

# array.c

# returns the maximum abs value in <a>, an array or list
def max_array(a):
    max = 0.0
    for i in a:
	i = abs(i)
	if (i>max): 
	    max = i
    return i
    
# multiplies all the values in the array or list <a> by <m>
def un_norm_array(a, m):
	count = 0
	for i in a:
		a[count] = i * m
		count = count + 1
	return a
    
# returns the average value of the array or list
# {!!Does not consider values less than 0 in the count of the average!!}
def prom_array(a):
    sum = 0.0
    count = 0
    for i in a:
		sum = sum + i
		if (i > 0.0):
			count = count + 1
    if (count > 0):
		return sum / count
    else:
		return 0.0


# ====================================
# ana-fun.cl
# analysis functions


# computes the number of frames in the analysis
# we want to have an extra frame at the end to prevent
# chopping the eneding
def compute_frames(total_samps, hop, st, nd):
	tmp = int(floor(total_samps / hop))
	tmp2 = (tmp * hop) - hop + st
	if (tmp2 > nd):
		return tmp
	else:
		return tmp + 1


# ====================================
# sound optimization

# optimize_sound <sound>
# initializes severl slots of <sound>
def optimize_sound (sound, verbose=False, get_max_values=True, fill_gaps=True,
					min_length=ATS_MIN_SEGMENT_LENGTH, trim=True,
					amp_threshold=ATS_AMP_THRESHOLD, simplify=True,
					min_frq=0.0, max_frq=20000.0):
	# double check parameters
	if (min_length==None):
		min_length=ATS_MIN_SEGMENT_LENGTH
	if (amp_threshold==None):
		amp_threshold=ATS_AMP_THRESHOLD
	# check previous optimization
	if (sound.optimized):
		print "Sound already optimized!"
	if (get_max_values):
		# get max amplitude of sound
		if (verbose):
			print "Getting Max. Amplitude..."
		sound.ampmax = get_ampmax(sound)
		if (verbose):
			print "Max. Amplitude: " + str(sound.ampmax)
		# get max frequency of sound
		if (verbose):
			print "Getting Max. Frequency..."
		sound.frqmax = get_frqmax(sound)
		if (verbose):
			print "Max. Frequency: " + str(sound.frqmax)
	if (fill_gaps):
		if (verbose):
			print "Filling out the sound gaps..."
		fill_sound_gaps(sound, min_length)
	if (trim):
		if (verbose):
			print "Trimming short segments off..."
		from ats_io import ats_save			
		ats_save(sound, "/pretrim.ats", True, False)
		trim_partials(sound, min_length)
	if (verbose):
		print "Getting amplitude and frequency averages..."
		ats_save(sound, "/preav.ats", True, False)
	set_amp_av(sound)
	set_frq_av(sound)
	if (simplify):
		if (verbose):
			print "Simplifying sound..."
		simplify_sound(sound, get_valid_partials(sound, min_frq, max_frq, amp_threshold))
	sound.optimized = True

		
# get_ampmax <sound>
# gets the maximum amplitude of <sound>
def get_ampmax(sound):
	ampmax = 0
	tmp = 0
	for h in range(0, sound.partials):
		tmp = max_array(sound.amp[h])
		if (tmp > ampmax):
			ampmax = tmp
	return ampmax

# get_frqmax <sound>
# gets the maximum frequency of <sound>
def get_frqmax(sound):
	frqmax = 0
	tmp = 0
	for h in range(0, sound.partials):
		tmp = max_array(sound.frq[h])
		if (tmp > frqmax):
			frqmax = tmp
	return frqmax

# fill up paramter for gaps shorter or equal to min_length
def fill_sound_gaps(sound, min_length):
	srate = sound.sampling_rate
	mag = TWO_PI / srate
	frame_size = sound.frame_size
	for par in range(0, sound.partials):
		amp_arr = sound.amp[par]
		frq_arr = sound.frq[par]
		next_val = 0
		gaps = get_gaps_arr(amp_arr)
		# first we fix the freq gap before attack
		next_val = find_next_val_arr(frq_arr, 0)
		if (next_val > 0):
			for f in range(0, next_val):
				frq_arr[f] = frq_arr[next_val]
		for k in gaps:
			if (k[1] <= min_length):
				left = 0
				if (k[0] != 0):
					left = k[0] - 1
				right = k[0] + k[1]
				if (right > sound.frames):
					right = sound.frames - 1
				# we know the boundaries of the gap, now let's fill it out...
				# frq
				for j in range(k[0], k[0] + k[1]):
					if (sound.frq[par][left] == 0.0):
						sound.frq[par][j] = sound.frq[par][right]
					elif (sound.amp[par][right] == 0.0):
						sound.frq[par][j] = sound.amp[par][left]
					else:
						sound.frq[par][j] = envelope_interp(j, [left, sound.frq[par][left],
																right, sound.frq[par][right]])
				# pha
				if (sound.amp[par][left] == 0.0):
					j = right - 1
					while (j > left):
						sound.pha[par][j] = sound.pha[par][j+1] + \
											(sound.frq[par][j] * mag * frame_size)
						j = j - 1
				else:
					for j in range(int(left + 1), right):
						sound.pha[par][j] = sound.pha[par][j-1] + \
											(sound.frq[par][j] * mag * frame_size)
				# and finally the amps
				for j in range(k[0], k[0]+k[1]):
					sound.amp[par][j] = envelope_interp(j, [left, sound.amp[par][left],
														right, sound.amp[par][right]])
					
# get gap positions and lengths
def get_gaps_arr(arr):
	segs = segments_arr(arr, l=[])
	st = 0
	l = []
	for k in segs:
		if (k[0] == 0):
			st = st + k[1]
		else:
			l.append([st, (k[0] - st)])
	return l
	
# returns a list of lists containing
# the index of the segment and its length
def segments_arr(arr, k=0, in_seg=False, st=0, l=[]):
	if (k == len(arr)):
		if (in_seg):
			l.append([st, (k - st)])
		else:
			l.reverse()
		return l
	else:
		if ((arr[k]==0.0) and in_seg):
			l.insert(0, [st, (k - st)])
			in_seg = False
		elif ((arr[k] != 0.0) and (not in_seg)):
			st = k
			in_seg = True
		return segments_arr(arr, k+1, in_seg, st, l)
	
# auxillary functions to fill functions	
def find_next_val_arr(arr, j):
	if ((arr[j] == 0) or (j == (len(arr) - 1))):
		return j
	else:
		return find_next_val_arr(arr, j+1)

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

				
# removes short segments from all partials
# this could take account of SMR, coming soon!
def trim_partials(sound, min_length=ATS_MIN_SEGMENT_LENGTH):
	for i in range(0, sound.partials):
		remove_short_segments(sound, i, min_length)

# removes short segments from a partial
def remove_short_segments(sound, par, min_length=ATS_MIN_SEGMENT_LENGTH):
	segments = segments_arr(sound.amp[par],l=[])
	for seg in segments:
		if (seg[1] < min_length):
			for i in range(seg[0], int(seg[0] + seg[1])):
				sound.amp[par][i] = 0.0

# sets the average amplitude for each partial of <sound>
def set_amp_av(sound):
	if (not sound.amp_av):
		sound.amp_av = zeros(sound.partials, 'Float64')
	frames = sound.frames
	for i in range(0, sound.partials):
		if (max_array(sound.amp[i]) > 0.0):
			sound.amp_av[i] = prom_array(sound.amp[i])

# sets the average frequency for each partial of <sound>
def set_frq_av(sound):
	if (not sound.frq_av):
		sound.frq_av = zeros(sound.partials, 'Float64')
	frames = sound.frames
	for i in range(0, sound.partials):
		if (max_array(sound.frq[i]) > 0.0):
			sound.frq_av[i] = prom_array(sound.frq[i])

# returns a list with the valid partial numbers
# valid partials are those with frq-av >= *ats-amp-threshold*
# and frq-av within min-frq and max-frq
def get_valid_partials(sound, min_frq=0.0, max_frq=20000.0, amp_threshold=ATS_AMP_THRESHOLD):
	l = []
	print db_amp(amp_threshold)
	for i in range(0, sound.partials):
		if ((sound.amp_av[i] >= db_amp(amp_threshold)) and \
				(min_frq <= sound.frq_av[i]) and \
				(sound.frq_av[i] <= max_frq)):
			l.append(i)
	print len(l)
	return l
		
# eliminates unvalid partials from <sound>
# valid partials in <valid> list		
def simplify_sound(sound, valid):
	n_partials = len(valid)
	print n_partials
	n_time = [0] * n_partials
	n_amp = [0] * n_partials
	n_frq = [0] * n_partials
	n_pha = []
	if (sound.pha):
		n_pha = [0] * n_partials
	n_noi = []
	if (sound.energy):
		n_noi = [0] * n_partials
	n_amp_av = zeros(n_partials, 'Float64')
	n_frq_av = zeros(n_partials, 'Float64')
	if (valid):
		tmp = []
		for i in valid:
			tmp.append([i,sound.frq_av[i]])
		sorted_valid = sort_partials_by_freq(tmp)
		i = 0
		for sv in sorted_valid:
			j = sv[0]
			n_time[i] = clone_list(sound.time[j])
			n_amp[i] = clone_list(sound.amp[j])
			n_frq[i] = clone_list(sound.frq[j])
			n_amp_av[i] = sound.amp_av[j]
			n_frq_av[i] = sound.frq_av[j]
			if (n_pha):
				n_pha[i] = clone_list(sound.pha[j])
			if (n_noi):
				n_noi[i] = clone_list(sound.energy[j])
		sound.time = n_time
		sound.amp = n_amp
		sound.frq = n_frq
		sound.amp_av = n_amp_av
		sound.frq_av = n_frq_av
		sound.partials = n_partials
		if (n_pha):
			sound.pha = n_pha
		if (n_noi):
			sound.energy = n_noi
						
# Quickstort implementation for sorting by partial frequency
def sort_partials_by_freq(arr):
	less = []
	pivotList = []
	greater = []
	if (len(arr) <= 1):
		return arr
	pivot = arr[int(len(arr) / 2)][1]
	for x in arr:
		if (x[1] < pivot):
			less.append(x)
		elif (x[1] == pivot):
			pivotList.append(x)
		else:
			greater.append(x)
	return sort_partials_by_freq(less) + pivotList + sort_partials_by_freq(greater)		
		
# ====================================
# plotting utilities

# simple plotting function: plots a graph of values (y-axis) along the index # of the array (x-axis)

def simple_plot(vals, x='indices', y='values', gridbool=False):
	a = [0.0] * len(vals)
	for count in range(0,len(vals)):
		a[count] = count
		count = count + 1
	plot(a, vals)
	xlabel(x)
	ylabel(y)
	grid(gridbool)
	show()
	
def plot_2(indices, values, x='indices', y='values', gridbool=False):
	plot(indices, values)
	xlabel(x)
	ylabel(y)
	grid(gridbool)
	show()

def plot_3(indices, values, x='indices', y='values', gridbool=False, marks='.'):
			plot(indices,values,marker=marks, ls=' ')
			bar(indices,values)
			xlabel(x)
			ylabel(y)
			show()
		

		
