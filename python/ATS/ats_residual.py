# python 2.5

# ATS
# by Juan Pampin
# pampin@u.washington.edu
#
# ported to python by
# Johnathan Lyon
# jgl@u.washington.edu
# 2007-06-21

#
# FILE: ats_residual.py
#
# This file contains the implementation
# of the ATS's residual analysis and 
# computation
#
# includes:
#	constants:
#		
#
#	functions:
#		compute_residual
#		read_frame
#		synth_buffer
#		compute_M
#		compute_aux
#		compute_alpha
#		compute_beta
#		interp_phase
#		wrap_phase



from pyaudiolab import sndfile
from numpy import zeros
from math import floor, cos

from ats_util import TWO_PI, max_array, un_norm_array


# Computes the difference between the synthesis and the original sound. 
# The function is passed an I/O struct (fil) pointing to the analyzed sound
# the <win-samps> array contains the sample numbers in the input file 
# corresponding to each frame.
def compute_residual(fil, output_fil, sound, win_samps, file_sampling_rate, \
						verbose=False, equalize=False, first_frame=0, \
						last_frame=None, first_par=0, last_par=None, \
						fmt=['aiff', 'pcm16']):

	out_fil = sndfile(output_fil, 'write', fmt, 2, file_sampling_rate)
	frames = last_frame
	if (frames == None):
		frames = sound.frames
	partials = last_par
	if (partials == None):
		partials = sound.partials
	frm_samps = win_samps[1] - win_samps[0]
	# 2D array to send to the output file at close
	outsize = int((frames * frm_samps) - first_frame)
	output = zeros([outsize ,2], 'Float64')
	in_buff = zeros(frm_samps, 'Float64')
	synth_buff = zeros(frm_samps, 'Float64')
	mag = TWO_PI / file_sampling_rate
	out_smp = 0
	# now we go over the whole sound computing the synthesis frame by frame
    # we store the residual on channel A and the synthesis on channel B
	if (verbose):
		print "Computing residual..."
	for frm in range((first_frame + 1), frames):
		frm_1 = frm - 1
		frm_2 = frm
		samp_1 = win_samps[frm_1]
		samp_2 = win_samps[frm_2]
		max_in = 0.0
		max_synth = 0.0
		gain = 1.0
		if (verbose):
			print "<Frame: " + str(frm_1) + " s1: " + str(samp_1) + " s2: " + str(samp_2)
		# read samples from input
		read_frame(fil, samp_1, samp_2, in_buff)
		# now we have to compute one synthesis frame
		# we clear the array first
		synth_buff = zeros(frm_samps, 'Float64')
		for par in range(first_par, partials):
			a1 = sound.amp[par][frm_1]
			a2 = sound.amp[par][frm_2]
			# must convert to radians per sample
			f1 = sound.frq[par][frm_1] * mag
			f2 = sound.frq[par][frm_2] * mag
			p1 = sound.pha[par][frm_1]
			p2 = sound.pha[par][frm_2]
			if (not ((a1 <= 0.0) and (a2 <= 0.0))):
				# check form amp 0 in frame 1
				if (a1 <= 0.0):
					delta = p2 - (f2 * frm_samps)
					f1 = f2
					p1 = wrap_phase(delta)
				if (a2 <= 0.0):
					delta = p1 + (f1 * frm_samps)
					f2 = f1
					p2 = wrap_phase(delta)
				# synthesize partial and store in buffer
				synth_buffer(a1,a2,f1,f2,p1,p2,synth_buff,frm_samps)
		if (equalize):
			max_synth = max_array(synth_buff)
			max_in = max_array(in_buff)
			gain = 1.0
			if (max_synth > 0.0):
				gain = max_in / max_synth
			if (gain != 1.0):
				un_norm_array(synth_buff, gain)
		for i in range(0, frm_samps):
			# now we write the residual and synthesis to the output
			output[out_smp,0] = in_buff[i] - synth_buff[i]
			output[out_smp,1] = synth_buff[i]
			out_smp = out_smp + 1
	# write to file and close
	out_fil.write_frames(output, len(output))
	out_fil.close()
		
		
# read one frame into a buffer		
def read_frame (fil, samp_1, samp_2, in_buffer):
	in_buffer = zeros(len(in_buffer), 'Float64')
	if ((samp_2 - samp_1) != len(in_buffer)):
		print "ERROR: wrong number of samples"
		return
	i =0
	for ptr in range(samp_1, samp_2):
		if (ptr >= len(fil)):
			in_buffer[i] = 0.0
		else:
			in_buffer[i] = fil[ptr]
		i = i + 1

# synthesizes a buffer using phase interpolation
def synth_buffer (a1, a2, f1, f2, p1, p2, buffer, frame_samps):
	a_inc = (a2 - a1) / frame_samps
	M = compute_M(p1, f1, p2, f2, frame_samps)
	aux = compute_aux(p1, p2, f1, frame_samps, M)
	alpha = compute_alpha(aux, f1, f2, frame_samps)
	beta = compute_beta(aux, f1, f2, frame_samps)
	amp = a1
	for k in range(0, frame_samps):
		buffer[k] = buffer[k] + (amp * cos(interp_phase(p1, f1, alpha, beta, k)))
		amp = a1 + a_inc

# ============================
# Utility Macros

def compute_M (pha_1, frq_1, pha, frq, buffer_size):
	return floor((((pha_1 + (frq_1 * buffer_size) - pha) + (0.5 * buffer_size * (frq - frq_1))) / TWO_PI) + 0.5)

def compute_aux (pha_1, pha, frq_1, buffer_size, M):
	return ((pha + (TWO_PI * M)) - (pha_1 + (frq_1 * buffer_size)))

def compute_alpha (aux, frq_1, frq, buffer_size):
	return (((3.0 / (buffer_size * buffer_size)) * aux) - ((frq - frq_1) / buffer_size) )

def compute_beta (aux, frq_1, frq, buffer_size):
	return (((-2.0 / (buffer_size * buffer_size * buffer_size)) * aux) + ((frq - frq_1) / (buffer_size * buffer_size)))

def interp_phase (pha_1, frq_1, alpha, beta, i):
	return (pha_1 + (frq_1 * i) + (alpha * i * i) + (beta * i * i * i))

def wrap_phase (phase):
	return (phase - (floor(phase / TWO_PI) * TWO_PI))
		
	
	
	
	
	
	
	