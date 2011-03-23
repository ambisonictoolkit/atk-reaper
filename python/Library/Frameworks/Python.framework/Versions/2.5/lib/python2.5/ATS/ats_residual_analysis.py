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
# FILE: ats_residual_analysis.py
#
# This file contains the residual analysis function
# definitions and constants
#
# includes:
#	constants:
#
#	functions:
# 		energy_to_band
#		get_band_partials


from ats_critical_bands import ATS_CRITICAL_BAND_EDGES


# transfers energy from partials to a band
def energy_to_band (sound, band, frame):
	lo_frq = ATS_CRITICAL_BAND_EDGES[band]
	hi_frq = ATS_CRITICAL_BAND_EDGES[band + 1]
	par = get_band_partials(lo_frq, hi_frq, sound, frame)
	sum = 0
	for p in par:
		sum = sum + sound.energy[p][frame]
	return sum

# returns a list of partial numbers that fall
# in frequency between lo and hi
def get_band_partials(lo, hi, sound, frame):
	par = []
	for k in range(0, int(sound.partials)):
		if (lo <= sound.frq[k][frame] <= hi):
			par.append(k)
	return par