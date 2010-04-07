# python 2.5

##### NOTES: check infinity on exp window


# ATS
# by Juan Pampin
# pampin@u.washington.edu
#
# ported to python by
# Johnathan Lyon
# jgl@u.washington.edu
# 2007-05-23

#
# FILE: ats_win.py
#
# This file contains the windowing function
# definitions and constants
#
# includes:
#	constants:
#		ATS_BLACKMAN_WINDOW_COEFF
#
#	functions:
#		make_blackman_window
#		make_fft_window
#		window_norm
#		bes_i0

from math import cos, exp, sqrt, sin, floor, log, pi, cosh
import cmath
from ats_util import TWO_PI
from numpy.fft import ifft

# All data coming form Harris' famous paper:
# "On the Use Of windows For Harmonic Analysis 
#  With The Discrete Fourier Transform"
# Proceedings of the IEEE, Vol. 66, No. 1 (pg. 51 to 84)
# January 1978
# Albert H. Nuttall, "Some Windows with Very Good Sidelobe Behaviour", 
# IEEE Transactions of Acoustics, Speech, and Signal Processing, Vol. ASSP-29,
# No. 1, February 1981, pp 84-91

# Window coefficients (a0, a1, a2, a3)
ATS_BLACKMAN_WINDOW_COEFF = [ [0.42659, -0.49656, 0.07685], # Exact Blackman (-51 dB)
			      [0.42, -0.5, 0.08], # Blackman (rounded coeffs) (-58 dB)
			      [0.42323, -0.49755, 0.07922], # 3-term Blackman-Harris 1 (-67 dB)
			      [0.44959, -0.49364, 0.05677], # 3-term Blackman-Harris 2 (-61 dB)
			      [0.35875, -0.48829, 0.14128, -0.01168], # 4-term Blackman-Harris 1 (-92 dB)
			      [0.40217, -0.49703, 0.09392, -0.00183]] # 4-term Blackman-Harris 2 (-71 dB)

# Window creation
# we generate a float list with the window values for each case
# @return: the window
def make_blackman_window(type, M):
    coeffs = {
	'blackman-exact': ATS_BLACKMAN_WINDOW_COEFF[0],
	'blackman': ATS_BLACKMAN_WINDOW_COEFF[1],
	'blackman-harris-3-1': ATS_BLACKMAN_WINDOW_COEFF[2],
	'blackman-harris-3-2': ATS_BLACKMAN_WINDOW_COEFF[3],
	'blackman-harris-4-1': ATS_BLACKMAN_WINDOW_COEFF[4],
	'blackman-harris-4-2': ATS_BLACKMAN_WINDOW_COEFF[5]} [type]
    two_pi_over_M = TWO_PI / M
    four_pi_over_M = 2 * two_pi_over_M
    six_pi_over_M = 3 * two_pi_over_M
    a0 = coeffs[0]
    a1 = coeffs[1]
    a2 = coeffs[2]
    if (len(coeffs)==4):
    	a3 = coeffs[3]
    else:
    	a3 = None
    win = [0.0] * M
    for count in range(0,M):
    	win[count] = a0 + (a1 * cos(two_pi_over_M * count)) + (a2 * cos(four_pi_over_M * count))
    	if (a3!=None):
    		win[count] = win[count] + (a3 * cos(six_pi_over_M * count))
    return win
    
    
    
def make_fft_window(type, M, beta=1.0, mu=0.0):
	window = [0.0] * M
	switch = {
		'rectangular' : [0, 'val=1.0'],
		'parzen' : [0, 'val=(1.0 - abs( (i - midn) / midp1))'],
		'welch' : [0, 'val=(1.0 - pow((float(i - midn) / midp1), 2))'],
		'kaiser' : [0, 'val=(bes_i0((beta * (sqrt(1.0 - pow(float(midn - i) / midn, 2))))) / I0beta)'],
		'gaussian' : [0, 'val=(exp( -0.5 * pow((beta * (float(midn - i)/ midn)), 2)))'],
		'poisson' : [0, 'val=(exp( -beta * (float(midn - i) / midn)))'],
		'cauchy' : [0, 'val=(1.0 / (1.0 + pow(((beta * float(midn - i)) / midn), 2)))'],
		'connes' : [0, 'val=pow((1.0 - pow( (float(i - midn) / midp1), 2)), 2)'],	
		### infinity check
		'exponential' : [1, 'val=(expsum - 1.0)', 'expsum=(expsum * expn)'],
		'bartlett' : [1, 'val=angle', 'angle=(angle + rate)'],
		'riemann' : [2, 'midn==i', 'val=1.0', 'val=(sin(sr * (midn - i)) / (sr * (midn - i)))'],
		'tukey' : [3, 'pos=(midn * (1.0 - beta))', 'i >= pos', 'val=1.0', 'val=(0.5 * (1.0 - cos( (pi * i) / pos)))'],
		'hamming' : [4, 'val=(0.54 - (0.46 * cx))'],
		'hann' : [4, 'val=(0.5 - (0.5 * cx))'],
		'hann-poisson' : [4,'val=((0.5 - (0.5 * cx)) * exp( -beta * (float(midn - i) / midn)))']
		} [type]		
	midn = int(floor(M / 2))
	midp1 = (M + 1) / 2
	freq = TWO_PI / M
	rate = 1.0 / midn
	sr = TWO_PI / M
	angle = 0.0
	expn = (1.0 + (log(2) / midn))
	expsum = 1.0
	I0beta = bes_i0(beta)
	val = 0.0
	j = M - 1
	if (switch[0]==0):
		for i in range(0,midn+1):
			exec(switch[1])
			window[i] = val
			window[j] = val
			j = j - 1
	elif (switch[0]==1):
		for i in range(0,midn+1):
			exec(switch[1])
			exec(switch[2])
			window[i] = val
			window[j] = val
			j = j - 1
	elif (switch[0]==2):
		for i in range(0,midn+1):
			if (eval(switch[1])):
				exec(switch[2])
			else:
				exec(switch[3])
			window[i] = val
			window[j] = val
			j = j - 1
	elif (switch[0]==3):
		exec(switch[1])
		for i in range(0,midn+1):
			if (eval(switch[2])):
				exec(switch[3])
			else:
				exec(switch[4])
			window[i] = val
			window[j] = val
			j = j - 1
	elif (switch[0]==4):
		for i in range(0,midn+1):
			cx = cos(angle)
			exec(switch[1])
			window[i] = val
			window[j] = val
			j = j - 1
			angle = angle + freq
	return window


# Returns the norm of the window
def window_norm (window):
	acc = 0
	M = len(window)
	count = 0
	while (count < M):
		acc = acc + abs(window[count])
		count = count + 1
	return (2.0 / acc)
	
	
# Modified Bessel Function of the First Kind
# from "Numerical Recipes in C"
def bes_i0 (x):
	if (abs(x) < 3.75):
		y = pow( (x / 3.75), 2)
		return (1.0 + (y * (3.5156229 + (y * (3.0899414 + \
			(y * (1.2067492 + (y * (0.2659732 + (y * \
			(0.360768e-1 + (y * 0.45813e-2))))))))))))
	else:
		ax = abs(x)
		y = 3.75 / ax
		return ( (exp(ax) / sqrt(ax)) * (0.39894228 + \
			(y * (0.1328592e-1 + (y * (0.225319e-2 + (y * \
			(-0.157565e-2 + (y * (0.916281e-2 + (y * \
			(-0.2057706e-1 + (y * (0.2635537e-1 + (y * \
			(-0.1647633e-1 + (y * 0.392377e-2)))))))))))))))))


