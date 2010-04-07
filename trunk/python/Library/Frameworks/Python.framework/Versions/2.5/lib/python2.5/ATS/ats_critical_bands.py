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
# FILE: ats_critical_bands.py
#
# This file contains the implementation
# of a masking curve evaluation
# algorithm using a critical band based model
#
# includes:
#	constants:
#		ATS_CRITICAL_BAND_EDGES
#
#	functions:
#		evalutate_smr
#		clear_mask
#		frq_to_bark
#		find_band
#		compute_slope_r

from math import log

from ats_util import amp_db_spl


ATS_CRITICAL_BAND_EDGES = [0.0, 100.0, 200.0, 300.0, 400.0, 510.0, 630.0, 
							770.0, 920.0, 1080.0, 1270.0, 1480.0, 1720.0, 
							2000.0, 2320.0, 2700.0, 3150.0, 3700.0, 4400.0, 
							5300.0, 6400.0, 7700.0, 9500.0, 12000.0, 15500.0, 
							20000.0]


# evaluates masking values (SMR) for peaks in list <peaks>
# [slope-r] and [slope-l] are the slopes of the mask
# in dBs/bark, <delta-db> is the dB treshold for
# the masking curves (must be <= 0dB) 

def evaluate_smr(peaks, slope_l=-27.0, delta_dB=-50, debug=False):
	clear_mask(peaks)
	if (len(peaks)==1):
		peaks[0].smr = amp_db_spl(peaks[0].amp)
	else:
		for i in range(0,len(peaks)):
			maskee = peaks[i]
			frq_maskee = frq_to_bark(maskee.frq)
			amp_maskee = amp_db_spl(maskee.amp)
			if (debug):
				print "frq-maskee: " + str(frq_maskee) + " amp-maskee: " + str(amp_maskee)
			pp = peaks[0:i] + peaks[i+1:len(peaks)]
			for j in range(0, len(pp)):
				frq_masker = frq_to_bark(pp[j].frq)
				amp_masker = amp_db_spl(pp[j].amp)
				slope_r = compute_slope_r(amp_masker)
				mask_term = 0
				if (debug):
					print "frq-masker: " + str(frq_masker) + " amp_masker: " + str(amp_masker)
					print "slope-r: " + str(slope_r)
				if (frq_masker < frq_maskee):
					mask_term = (amp_masker + delta_dB + ((frq_maskee - frq_masker) * slope_r))
				else:
					mask_term = (amp_masker + delta_dB + ((frq_masker - frq_maskee) * slope_l))
				if (mask_term > maskee.smr):
					maskee.smr = mask_term
			if (debug):
				print "Maskee SMR: " + str(maskee.smr)
			maskee.smr = amp_maskee - maskee.smr
	
	
# sets mask values of all peaks in list <peaks> to 0.0

def clear_mask(peaks):
	for x in peaks:
		x.smr = 0.0


# NOTE WILL RETURN NONE???

# frq_to_bark <frq>
# converts <frq> into bark scale

def frq_to_bark(frq):
	if (frq <= 400.0):
		return 0.01 * frq
	elif (frq >= 20000.0):
		return None
	else:
		band = find_band(frq)
		lo_frq = ATS_CRITICAL_BAND_EDGES[band]
		hi_frq = ATS_CRITICAL_BAND_EDGES[int(band + 1)]
		return (1 + band + abs(log((frq / lo_frq), 10) / log((lo_frq / hi_frq), 10)))
		

# find_band <frq>
# finds the critical band for <frq>

def find_band (frq, l=ATS_CRITICAL_BAND_EDGES, b=0):
	if ((len(l)==0) or (l[0] > frq)):
		return (b - 1)
	else:
		return find_band(frq, l[1:len(l)], (b + 1))


def compute_slope_r(masker_amp_db):
	temp = masker_amp_db - 40.0
	if (temp < 0.0):
		temp = 0.0
	return (-27.0 + (temp * 0.37))

