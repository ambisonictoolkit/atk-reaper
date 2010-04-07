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
# FILE: ats_io.py
#
# This file contains the i/o function
# definitions and constants
#
# includes:
#	constants:
#		ATS_HEADER_SIZE
#		ATS_MAGIC_NUMBER
#
#	functions:
#		ats_save
#		write_array_double_to_binary

from numpy import zeros
from struct import pack

from ats_util import ATS_CRITICAL_BANDS
from ats_residual_analysis import energy_to_band


# ===================================
# CONSTANTS

ATS_HEADER_SIZE = 10
ATS_MAGIC_NUMBER = 123.0



"""
-ATS header consists of (all double floats):

ATS_MAGIC_NUMBER
sampling-rate (samples/sec)
frame-size (samples)
window-size (samples)
partials (number)
frames (number)
ampmax (max. amplitude)
frqmax (max. frequecny)
dur (duration)
type (number, see below)

-ATS frames can be of four different types:

1) without phase or noise:
==========================
time (frame starting time)
amp (par#0 amplitude)
frq (par#0 frequency)
...
amp (par#n amplitude)
frq (par#n frequency)


2) with phase but not noise:
============================
time (frame starting time)
amp (par#0 amplitude)
frq (par#0 frequency)
pha (par#0 phase)
...
amp (par#n amplitude)
frq (par#n frequency)
pha (par#n phase)


3) with noise but not phase:
============================
time (frame starting time)
amp (par#0 amplitude)
frq (par#0 frequency)
...
amp (par#n amplitude)
frq (par#n frequency)

energy (band#0 energy)
...
energy (band#n energy)

4) with phase and noise:
========================
time (frame starting time)
amp (par#0 amplitude)
frq (par#0 frequency)
pha (par#0 phase)
...
amp (par#n amplitude)
frq (par#n frequency)
pha (par#n phase)

noise (band#0 energy)
...
noise (band#n energy)
"""


# saves <sound> into <file> in binary format
# <file> must be string with a file name
# in case the file already exists it gets overwritten

def ats_save (sound, file, save_phase=True, save_noise=True):
		fd = open(file, 'w')
		sr = float(sound.sampling_rate)
		frame_size = float(sound.frame_size)
		window_size = float(sound.window_size)
		partials = float(sound.partials)
		frames = float(sound.frames)
		max_frq = float(sound.frqmax)
		max_amp = float(sound.ampmax)
		dur = float(sound.dur)
		has_pha = False
		if (save_phase and sound.pha):
			has_pha = True
		has_noi = False
		if (save_noise and (sound.energy or sound.band_energy)):
			has_noi = True
		type = 1
		if (has_pha and has_noi):
			type = 4
		elif ((not has_pha) and has_noi):
			type = 3
		elif (has_pha and (not has_noi)):
			type = 2
		time_arr = zeros(1, 'Float64')
		header_arr = zeros(ATS_HEADER_SIZE, 'Float64')
		hop_arr = zeros(sound.frames, 'Float64')
		band_l = None
		if (has_noi):
			band_l = sound.bands
		data_arr = None
		noi_arr = None
    	# Header:
    	# [mag, sr, fs, ws, #par, #frm, MaxAmp, MaxFrq, dur, type]
   		# Simple mag word for now
    	# would read 123.0 if read with the correct byte order
		header_arr[0] = ATS_MAGIC_NUMBER
		header_arr[1] = sr
		header_arr[2] = frame_size
		header_arr[3] = window_size
		header_arr[4] = partials
		header_arr[5] = frames
		header_arr[6] = max_amp
		header_arr[7] = max_frq
		header_arr[8] = dur
		header_arr[9] = type
		print "Mag: " + str(header_arr[0]) + " SR: " + str(header_arr[1]) 	
		print "FS: " + str(header_arr[2]) + " WS: " + str(header_arr[3]) 
		print "Partials: " + str(header_arr[4]) + " Frames: " + str(header_arr[5]) 
		print "MaxAmp: " + str(header_arr[6]) + " MaxFrq: " + str(header_arr[7]) 
		print "Dur: " + str(header_arr[8]) + " Type: " + str(header_arr[9]) 
		# create array for data
		si = 2
		if (has_pha):
			si = 3
		data_arr = zeros(si, 'Float64')
		if (has_noi):
			noi_arr = zeros(ATS_CRITICAL_BANDS, 'Float64')
		# Store all times in array.
		# For now we consider all partials
		# have the same time structure
		# (need a different file format for multirate)
		for h in range(0, int(frames)):
			hop_arr[h] = sound.time[0][h]
		print "Saving sound..."
		for i in range(0, int(frames)):
			# write header
			if (i==0):
				write_array_double_to_binary(fd, header_arr)
			# write time
			time_arr[0] = hop_arr[i]
			write_array_double_to_binary(fd, time_arr)
			# loop for all partials and write blocks
			# of [par#, amp, frq, pha]
			for j in range(0, int(partials)):
				data_arr[0] = sound.amp[j][i]
				data_arr[1] = sound.frq[j][i]
				if (has_pha):
					data_arr[2] = sound.pha[j][i]
				# write buffer
				write_array_double_to_binary(fd, data_arr)
			# noise part
			if (has_noi):
				# NOTE:
				# for now critical band energy is stored as an array
				# of <frames> arrays of 25 elements each
				for k in range(0, ATS_CRITICAL_BANDS):
					if (band_l!=None and (k in band_l)):
						noi_arr[k] = float(sound.energy[index(k)][i])
					else:
						noi_arr[k] = float(energy_to_band(sound, k, i))
				write_array_double_to_binary(fd, noi_arr)
		# close the file
		fd.close()





def write_array_double_to_binary (file, arr):
	for i in range(0, len(arr)):
		file.write(pack('d',arr[i]))











