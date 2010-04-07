# python 2.5

# ATS
# by Juan Pampin
# pampin@u.washington.edu
#
# ported to python by
# Johnathan Lyon
# jgl@u.washington.edu
# 2007-05-23

#
# FILE: ats_peak_track.py
#
# This file contains the implementation
# of the ATS's peak tracking algorithm
#
# includes:
#	constants:
#		
#
#	functions:
#		update_tracks
#		find_peak_by_track
#		peak_tracking
#		sort_tracks_by_smr

from ats_peak_detect import sort_peaks_by_freq
from ats_util import clone_list

def update_tracks (tracks, track_length, frame_n, ana_frames, beta=0.0):
	if (len(tracks) > 0):
		frames = min(frame_n, track_length)
		first_frame = frame_n - frames
		for g in tracks:
			track = g.track
			frq_acc = 0
			f = 0
			amp_acc = 0
			a = 0
			smr_acc = 0
			s = 0
			last_frq = 0
			last_amp = 0
			last_smr = 0
			for i in range(first_frame,frame_n):
				l_peaks = ana_frames[i]
				peak = find_peak_by_track(track, l_peaks)
				if (peak != None):
					if (peak.frq > 0.0):
						last_frq = peak.frq
						frq_acc = frq_acc + peak.frq
						f = f + 1
					if (peak.amp > 0.0):
						last_amp = peak.amp
						amp_acc = amp_acc + peak.amp
						a = a + 1
					if (peak.smr > 0.0):
						last_smr = peak.smr
						smr_acc = smr_acc + peak.smr
						s = s + 1
			if (f > 0):
				g.frq = ((1 - beta) * (frq_acc / f)) + (beta * last_frq)
			if (a > 0):
				g.amp = ((1 - beta) * (amp_acc / a)) + (beta * last_amp)
			if (s > 0):
				g.smr = ((1 - beta) * (smr_acc / s)) + (beta * last_smr)
		return tracks
			
	else:
		return clone_list(ana_frames[int(frame_n - 1)])
				

# left-to-right search for a track in a list of peaks

def find_peak_by_track(track, peaklist):
	for p in peaklist:
		if (track==p.track):
			return p
	return None

# left-to-right search for a frq in a list of peaks

def find_index_by_frq(frq, peaklist):
	i = 0
	for p in peaklist:
		if (frq==p.frq):
			return i
		else:
			i = i + 1
	
	
# peak-tracking <tracks> <peaks-b>
# given a list of tracks <tracks> and a list of peaks 
# <peaks-b> this function tracks stable 
# sinusoidal trajectories between the two frames.
# It returns  a list with two lists
# of untracked peaks: first the untracked ones from
# <tracks>, second the untracked ones from <peaks-b>.
# The <frq-deviation> parameter is used to decide
# which peaks in <peaks-b> are candidates for a 
# particular peak and track combination, this value is multiplied 
# by the frequency of the track that is used as a central 
# axis for the candidate search in the <peaks-b> pool.
# NOTE: this function assumes that the client passes
# peaks in <peaks-a> sorted by masking value

def peak_tracking (tracks, peaks_b, frq_deviation=0.45, alpha=0.0, unmatched=[]):
	tracks_t = clone_list(tracks)
	peaki = clone_list(peaks_b)
	tracks_t = sort_tracks_by_smr(tracks_t)
	if ((len(tracks_t) < 1) or (len(peaki) < 1)):
		return [unmatched + tracks_t, peaki]
	peak = tracks_t[0]
	peak_frq = peak.frq
	peak_smr = peak.smr
	# find the frq limits for candidates
	frq_limits = get_limits(peak_frq, frq_deviation)
	# get possible candidates
	peaki = sort_peaks_by_freq(peaki)
	peak_candidates = find_candidates(peaki, frq_limits, []) 
	# find best candidate
	matched_peak = find_best_candidate(peak_candidates, peak_frq, peak_smr, alpha)
	if (matched_peak != None):
		peaki.remove(matched_peak)
		matched_peak.track = peak.track	
	else:
		unmatched.insert(0,peak)
	return peak_tracking(tracks_t[1:len(tracks_t)], peaki, frq_deviation, alpha, unmatched)
	
# get-limits <peak-frq> <frq-dev>
# returns a list with the high and low frq
# limit for <peak-frq> given <ftq-dev>

def get_limits(peak_frq, frq_dev):
	half_band = 0.5 * peak_frq * frq_dev
	lo = peak_frq - half_band
	hi = peak_frq + half_band
	return [lo, hi]
	
	
# find-candidates <peaks> <frq-limits>
# returns a list with candidate peaks that 
# fall within the frq-limits in <peaks>
# NOTE: this function assumes that <peaks>
# is sorted by frequency

def find_candidates (peaks, frq_limits, l=[]):
	lo = frq_limits[0]
	hi = frq_limits[1]
	if ((len(peaks) < 1) or (peaks[0].frq > hi)):	
		return l
	if ((lo <= peaks[0].frq) and (peaks[0].frq <= hi)):
		l.append(peaks[0])
	return find_candidates(peaks[1:len(peaks)], frq_limits, l)

# find-best-candidate <peak-candidates> <peak-frq>
# returns the peak from the <peak-candidates> list
# that is closer in frequency to <peak-frq>

def find_best_candidate (peak_candidates, peak_frq, peak_smr, alpha=0.0):
	# delta will have the special value 'inf'
	delta = 1e30000
	best_peak = None
	local_delta = None
	for p in peak_candidates:
		local_delta = (abs(p.frq - peak_frq) + (alpha * abs(p.smr - peak_smr))) / (alpha + 1)
		if (local_delta < delta):
			best_peak = p
			delta = local_delta
	return best_peak
	

# quicksort implementation for sorting a list of peaks <tracks> by smr values

def sort_tracks_by_smr(tracks):
	less = []
	pivotList = []
	greater = []
	if (len(tracks) <= 1):
		return tracks
	pivot = tracks[int(len(tracks) / 2)].smr
	for x in tracks:
		if (x.smr < pivot):
			less.append(x)
		elif (x.smr == pivot):
			pivotList.append(x)
		else:
			greater.append(x)
	return sort_tracks_by_smr(less) + pivotList + sort_tracks_by_smr(greater)