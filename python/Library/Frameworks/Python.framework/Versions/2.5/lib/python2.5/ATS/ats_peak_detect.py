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
# FILE: ats_peak_detect.py
#
# This file contains the implementation
# of the ATS's peak detection algorithm
#
# includes:
#
#	functions:
#		peak_detection
#		to_polar
#		parabolic_interp
#		phase_interp
#		sort_peaks_by_freq

from math import sqrt, atan, pi
from numpy import zeros
from ats_struct import ats_peak
from ats_util import amp_db, db_amp, TWO_PI


def peak_detection(ats_fft, lowest_bin=None, highest_bin=None, 
					lowest_magnitude=0.0, norm=1.0):

	# empty list of peaks
	peaks = []
	# last bin to search
	if (highest_bin==None):
		N = ats_fft.size
	else:
		N = highest_bin
	# first bin to start search
	if ((lowest_bin==None) or (lowest_bin < 2)):
		first_bin = 2
	else:
		first_bin = lowest_bin
	# frequency scaling
	fft_mag = float(ats_fft.rate) / ats_fft.size
	# list of magnitudes
	fftmags = zeros(N, 'Float64')
	# list of phases
	fftphase = zeros(N, 'Float64')
	# peak bins
	right_bin = 0.0 #first point
	left_bin = 0.0 # third point
	central_bin = 0.0 # central point
	#temp vars for data
	frq = 0.0
	pha = 0.0
	# convert spectrum to polar coordinates
	to_polar(ats_fft, fftmags, fftphase, N, norm)
	central_bin = fftmags[int(first_bin - 2)]
	right_bin = fftmags[int(first_bin - 1)]
	# peak detection:
	# move by one bin and analyze by groups of 3
	for k in range(first_bin, N):
		atspeak = ats_peak()
		left_bin = central_bin
		central_bin = right_bin
		right_bin = fftmags[k]
		if ((central_bin > lowest_magnitude) and \
			(central_bin > right_bin) and \
			(central_bin > left_bin)):
			offset, amp = parabolic_interp(left_bin, central_bin, right_bin)
			frq = fft_mag * (k - 1 + offset) # actual frq of peak
			if (offset > 0):
				pha = phase_interp(fftphase[int(k - 2)], fftphase[k], offset)
			else:
				pha = phase_interp(fftphase[int(k - 1)], fftphase[k], offset)
			atspeak.amp = amp
			atspeak.frq = frq
			atspeak.pha = pha
			peaks.append(atspeak)
	return sort_peaks_by_freq(peaks)
			
			
		
	
	
	
# ==========================================================	
# HELPER FUNCTIONS

# to-polar <ats-fft> <mags> <phase> <N>
# returns the magnitude and phases of <ats-fft>
# into the passed arrays <mags> and <phase>.
# <N> is the highest bin where we look for spectral
# data in the fft. We normalize mag values by [norm=1.0]

def to_polar(ats_fft, mags, phase, N, norm=1.0):
	for k in range(0,N):
		x = ats_fft.fd[k].real
		y = ats_fft.fd[k].imag
		mags[k] = norm * sqrt((x * x) + (y * y))
		if ((x==0.0) and (y==0.0)):
			phase[k] = 0.0
		else:
			phase[k] = atan( -y / x)


# parabolic-interp <alpha> <beta> <gamma>
# does parabolic interpolation of 3 points
# returns the x offset and heigth
# of the interpolated peak		

def parabolic_interp(alpha, beta, gamma):
	dB_alpha = amp_db(alpha)
	dB_beta = amp_db(beta)
	dB_gamma = amp_db(gamma)
	offset  = 0.5 * ((dB_alpha - dB_gamma) / (dB_alpha + (-2 * dB_beta) + dB_gamma))
	height = db_amp(dB_beta - (0.25 * (dB_alpha - dB_gamma) * offset))
	return [offset, height]
	

def phase_interp(peak_phase, right_phase, offset):
	if ((peak_phase - right_phase) > (1.5 * pi)):
		right_phase = right_phase + TWO_PI
	elif ((right_phase - peak_phase) > (1.5 * pi)):
		peak_phase = peak_phase + TWO_PI
	return (peak_phase + (offset * (right_phase - peak_phase)))
	
	
# Quickstort implementation for sorting by peak frequency
def sort_peaks_by_freq(array_of_peaks):
	less = []
	pivotList = []
	greater = []
	if (len(array_of_peaks) <= 1):
		return array_of_peaks
	pivot = array_of_peaks[int(len(array_of_peaks) / 2)].frq
	for x in array_of_peaks:
		if (x.frq < pivot):
			less.append(x)
		elif (x.frq == pivot):
			pivotList.append(x)
		else:
			greater.append(x)
	return sort_peaks_by_freq(less) + pivotList + sort_peaks_by_freq(greater)
